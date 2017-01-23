module SIS_tracer_advect
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of SIS2.                                        *
!*                                                                     *
!* SIS2 is free software; you can redistribute it and/or modify it and *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* SIS2 is distributed in the hope that it will be useful, but WITHOUT *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, 1996 - 2012, adapted for SIS2 in 2014-2016.    *
!*                                                                     *
!*    This program contains the subroutines that advect tracers        *
!*  horizontally (i.e. along layers).  This code was modified from the *
!*  corresponding MOM6 / GOLD code to work with the snow and ice       *
!*  tracers of SIS2.                                                   *
!*                                                                     *
!*    advect_SIS_tracers advects tracer concentrations using the       *
!*  modified flux advection scheme from Easter (Mon. Wea. Rev., 1993)  *
!*  with tracer distributions that are piecewise constant,             *
!*  piecewise linear (given by the monotonic scheme proposed by        *
!*  Lin et al. (Mon. Wea. Rev., 1994)), or the montonic piecewise      *
!*  parabolic method, as described in Carpenter et al. (MWR, 1990).    *
!*  This detects the mass of ice or snow in a grid cell and thickness  *
!*  category at the previous instance when the tracer concentration    *
!*  was changed is consistent with the mass fluxes and the new masses. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v, vh                                    *
!*    j    x ^ x ^ x   At >:  u, uh                                    *
!*    j    > o > o >   At o:  tr, h                                    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_MODULE, CLOCK_ROUTINE
use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_SIS_diag_field, safe_alloc_ptr, time_type
use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, max_across_PEs
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser, only : get_param, log_version, param_file_type
use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type
use SIS_tracer_registry, only : SIS_tracer_registry_type, SIS_tracer_type, SIS_tracer_chksum

implicit none ; private

#include <SIS2_memory.h>

public advect_SIS_tracers, advect_tracers_thicker, advect_scalar
public SIS_tracer_advect_init, SIS_tracer_advect_end

type, public :: SIS_tracer_advect_CS ; private
  real    :: dt             ! The baroclinic dynamics time step, in s.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                            ! timing of diagnostic output.
  logical :: debug          ! If true, write verbose checksums for debugging purposes.
  logical :: use_upwind2d   ! If true, use the non-split upwind scheme that was
                            ! was used in older versions of SIS.
  logical :: usePPM         ! If true, use PPM tracer advection instead of PLM.
  logical :: usePCM         ! If true, use PCM tracer advection instead of PLM.
end type SIS_tracer_advect_CS

logical :: first_call = .true.
integer :: id_clock_advect, id_clock_pass, id_clock_sync

contains

subroutine advect_SIS_tracers(h_prev, h_end, uhtr, vhtr, dt, G, IG, CS, TrReg, snow_tr ) ! (, OBC)
  type(SIS_hor_grid_type),                     intent(inout) :: G
  type(ice_grid_type),                         intent(inout) :: IG
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(in) :: h_prev, h_end
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), intent(in) :: uhtr
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), intent(in) :: vhtr
  real,                                        intent(in)    :: dt
  type(SIS_tracer_advect_CS),                  pointer       :: CS
  type(SIS_tracer_registry_type),              pointer       :: TrReg
  logical,                                     intent(in) :: snow_tr

! Arguments: h_prev - Category thickness times fractional coverage before advection, in m or kg m-2.
!  (in)      h_end - Layer thickness times fractional coverage after advection, in m or kg m-2.
!  (in)      uhtr - Accumulated volume or mass fluxes through zonal faces,
!                   in m3 s-1 or kg s-1.
!  (in)      vhtr - Accumulated volume or mass fluxes through meridional faces,
!                   in m3 s-1 or kg s-1.
!!  (in)      OBC - This open boundary condition type specifies whether, where,
!!                  and what open boundary conditions are used.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 SIS_tracer_advect_init.
!  (in)      TrReg - A pointer to the tracer registry.
!  (in)      snow_tr - If true, advect the snow tracers, otherwise advect the
!                      ice tracers.

  integer ntr

  if (.not. associated(CS)) call SIS_error(FATAL, "SIS_tracer_advect: "// &
       "SIS_tracer_advect_init must be called before advect_tracer.")
  if (.not. associated(TrReg)) call SIS_error(FATAL, "SIS_tracer_advect: "// &
       "register_tracer must be called before advect_tracer.")
  ntr = TrReg%ntr
  if (ntr==0) return

  call cpu_clock_begin(id_clock_advect)
  if (snow_tr) then
    if (CS%use_upwind2d) then
      call advect_upwind_2d(TrReg%Tr_snow, h_prev, h_end, uhtr, vhtr, ntr, dt, G, IG)
    else
      call advect_tracer(TrReg%Tr_snow, h_prev, h_end, uhtr, vhtr, ntr, dt, G, IG, CS)
    endif
  else
    if (CS%use_upwind2d) then
      call advect_upwind_2d(TrReg%Tr_ice, h_prev, h_end, uhtr, vhtr, ntr, dt, G, IG)
    else
      call advect_tracer(TrReg%Tr_ice, h_prev, h_end, uhtr, vhtr, ntr, dt, G, IG, CS)
    endif
  endif
  call cpu_clock_end(id_clock_advect)

end subroutine advect_SIS_tracers

subroutine advect_tracer(Tr, h_prev, h_end, uhtr, vhtr, ntr, dt, G, IG, CS) ! (, OBC)
  type(SIS_tracer_type), dimension(ntr),       intent(inout) :: Tr
  type(SIS_hor_grid_type),                     intent(inout) :: G
  type(ice_grid_type),                         intent(inout) :: IG
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(in)    :: h_prev, h_end
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), intent(in)    :: uhtr
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), intent(in)    :: vhtr
  real,                                        intent(in)    :: dt
  integer,                                     intent(in)    :: ntr
  type(SIS_tracer_advect_CS),                  pointer       :: CS
!!  type(ocean_OBC_type),                      pointer       :: OBC
!    This subroutine time steps the tracer concentration.
!  A monotonic, conservative, weakly diffusive scheme is used.

! Arguments: h_prev - Category thickness times fractional coverage before advection, in m or kg m-2.
!  (in)      h_end - Layer thickness times fractional coverage after advection, in m or kg m-2.
!  (in)      uhtr - Accumulated volume or mass fluxes through zonal faces,
!                   in m3 s-1 or kg s-1.
!  (in)      vhtr - Accumulated volume or mass fluxes through meridional faces,
!                   in m3 s-1 or kg s-1.
!!  (in)      OBC - This open boundary condition type specifies whether, where,
!!                  and what open boundary conditions are used.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 SIS_tracer_advect_init.
!  (in)      TrReg - A pointer to the tracer registry.

  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)) :: &
    hprev           ! The cell volume at the end of the previous tracer
                    ! change, in m3.
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)) :: &
    uhr             ! The remaining zonal thickness flux, in m3.
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)) :: &
    vhr             ! The remaining meridional thickness fluxes, in m3.
  real :: uh_neglect(SZIB_(G),SZJ_(G)) ! uh_neglect and vh_neglect are the
  real :: vh_neglect(SZI_(G),SZJB_(G)) ! magnitude of remaining transports that
                                ! can be simply discarded, in m3 or kg.

  real :: landvolfill         ! An arbitrary? nonzero cell volume, m3.
  real :: Idt                 ! 1/dt in s-1.
  real :: h_neglect
  logical :: domore_u(SZJ_(G),SZCAT_(IG))  ! domore__ indicate whether there is more
  logical :: domore_v(SZJB_(G),SZCAT_(IG)) ! advection to be done in the corresponding
                                ! row or column.
  logical :: x_first            ! If true, advect in the x-direction first.
  integer :: max_iter           ! The maximum number of iterations in
                                ! each layer.
  integer :: domore_k(SZCAT_(IG))
  integer :: stensil            ! The stensil of the advection scheme.
  integer :: nsten_halo         ! The number of stensils that fit in the halos.
  integer :: i, j, k, l, m, is, ie, js, je, isd, ied, jsd, jed
  integer :: ncat, nL_max, itt, do_any
  integer :: isv, iev, jsv, jev ! The valid range of the indices.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  landvolfill = 1.0e-20         ! This is arbitrary, but must be positive.
  stensil = 2                   ! The scheme's stensil; 2 for PLM.

  if (.not. associated(CS)) call SIS_error(FATAL, "SIS_tracer_advect: "// &
       "SIS_tracer_advect_init must be called before advect_tracer.")
  if (ntr==0) return

  x_first = (MOD(G%first_direction,2) == 0)

  Idt = 1.0/dt

  nL_max = 0
  do m=1,ntr ; nL_max = max(Tr(m)%nL,nL_max) ; enddo

  max_iter = 3
  if (CS%dt > 0.0) max_iter = 2*INT(CEILING(dt/CS%dt)) + 1

! This initializes the halos of uhr and vhr because pass_vector might do
! calculations on them, even though they are never used.
  uhr(:,:,:) = 0.0 ; vhr(:,:,:) = 0.0
  hprev(:,:,:) = landvolfill
  h_neglect = IG%H_subroundoff
! Initialize domore_u and domore_v to .false.; they will be reevaluated later.
  domore_u(:,:) = .false. ; domore_v(:,:) = .false.

!$OMP parallel default(none) shared(ncat,is,ie,js,je,domore_k,uhr,vhr,uhtr,vhtr,dt, &
!$OMP                               hprev,G,h_prev,h_end,isd,ied,jsd,jed,uh_neglect, &
!$OMP                               h_neglect,vh_neglect,ntr,Tr,domore_u,domore_v)
!$OMP do
  do k=1,ncat
    domore_k(k)=1
!  Put the remaining (total) thickness fluxes into uhr and vhr.
    do j=js,je ; do I=is-1,ie ; uhr(I,j,k) = dt*uhtr(I,j,k) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhr(i,J,k) = dt*vhtr(i,J,k) ; enddo ; enddo
    ! Find the previous total mass (or volume) of ice, but in the case that this
    ! category is now dramatically thinner than it was previously, add a tiny
    ! bit of extra mass to avoid nonsensical tracer concentrations.  This will
    ! lead rarely to a very slight non-conservation of tracers, but not mass.
    do j=js,je; do i=is,ie
      hprev(i,j,k) = G%areaT(i,j) * (h_prev(i,j,k) + &
                       max(0.0, 1.0e-13*h_prev(i,j,k) - h_end(i,j,k)))
      if (h_end(i,j,k) - h_prev(i,j,k) + ((uhr(I,j,k) - uhr(I-1,j,k)) + &
                            (vhr(i,J,k) - vhr(i,J-1,k))) * G%IareaT(i,j) > &
          1e-10*(h_end(i,j,k) + h_prev(i,j,k))) then
!$OMP critical
        call SIS_error(WARNING, "Apparently inconsistent h_prev, h_end, uhr and vhr in advect_tracer.")
!$OMP end critical
      endif
    enddo ; enddo
  enddo
!$OMP end do nowait
!$OMP do
  do j=jsd,jed ; do I=isd,ied-1
    uh_neglect(I,j) = h_neglect*MIN(G%areaT(i,j),G%areaT(i+1,j))
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsd,jed-1 ; do i=isd,ied
    vh_neglect(i,J) = h_neglect*MIN(G%areaT(i,j),G%areaT(i,j+1))
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do m=1,ntr
    if (associated(Tr(m)%ad2d_x)) then
      do j=jsd,jed ; do i=isd,ied ; Tr(m)%ad2d_x(I,j) = 0.0 ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad2d_y)) then
      do J=jsd,jed ; do i=isd,ied ; Tr(m)%ad2d_y(i,J) = 0.0 ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad3d_x)) then
      do k=1,ncat ; do j=jsd,jed ; do i=isd,ied
        Tr(m)%ad3d_x(I,j,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad3d_y)) then
      do k=1,ncat ; do J=jsd,jed ; do i=isd,ied
        Tr(m)%ad3d_y(i,J,k) = 0.0
      enddo ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad4d_x)) then
      do l=1,Tr(m)%nL ; do k=1,ncat ; do j=jsd,jed ; do i=isd,ied
        Tr(m)%ad4d_x(I,j,k,l) = 0.0
      enddo ; enddo ; enddo ; enddo
    endif
    if (associated(Tr(m)%ad4d_y)) then
      do l=1,Tr(m)%nL ; do k=1,ncat ; do J=jsd,jed ; do i=isd,ied
        Tr(m)%ad4d_y(i,J,k,l) = 0.0
      enddo ; enddo ; enddo ; enddo
    endif
  enddo
!$OMP end parallel

  isv = is ; iev = ie ; jsv = js ; jev = je

  do itt=1,max_iter

    if (isv > is-stensil) then
      call cpu_clock_begin(id_clock_pass)
      call pass_vector(uhr, vhr, G%Domain)
      do m=1,ntr ; do l=1,Tr(m)%nL
        call pass_var(Tr(m)%t(:,:,:,l), G%Domain, complete=.false.)
      enddo ; enddo
      call pass_var(hprev, G%Domain)
      call cpu_clock_end(id_clock_pass)

      nsten_halo = min(is-isd,ied-ie,js-jsd,jed-je)/stensil
      isv = is-nsten_halo*stensil ; jsv = js-nsten_halo*stensil
      iev = ie+nsten_halo*stensil ; jev = je+nsten_halo*stensil
      ! Reevaluate domore_u & domore_v unless the valid range is the same size as
      ! before.  Also, do this if there is Strang splitting.
      if ((nsten_halo > 1) .or. (itt==1)) then
!$OMP parallel do default(none) shared(ncat,domore_k,isv,iev,jsv,jev,stensil, &
!$OMP                                  domore_u,uhr,vhr,domore_v)
        do k=1,ncat ; if (domore_k(k) > 0) then
          do j=jsv,jev ; if (.not.domore_u(j,k)) then
            do i=isv+stensil-1,iev-stensil; if (uhr(I,j,k) /= 0.0) then
              domore_u(j,k) = .true. ; exit
            endif ; enddo ! i-loop
          endif ; enddo
          do J=jsv+stensil-1,jev-stensil ; if (.not.domore_v(J,k)) then
            do i=isv+stensil,iev-stensil; if (vhr(i,J,k) /= 0.0) then
              domore_v(J,k) = .true. ; exit
            endif ; enddo ! i-loop
          endif ; enddo

          !   At this point, domore_k is global.  Change it so that it indicates
          ! whether any work is needed on a layer on this processor.
          domore_k(k) = 0
          do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
          do J=jsv+stensil-1,jev-stensil ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo

        endif ; enddo ! k-loop
      endif
    endif

    ! Set the range of valid points after this iteration.
    isv = isv + stensil ; iev = iev - stensil
    jsv = jsv + stensil ; jev = jev - stensil
!$OMP parallel do default(none) shared(ncat,domore_k,x_first,Tr,hprev,uhr,uh_neglect, &
!$OMP                                  domore_u,ntr,nL_max,Idt,isv,iev,jsv,jev,       &
!$OMP                                  stensil,G,CS,vhr,vh_neglect,domore_v,IG)
    do k=1,ncat ; if (domore_k(k) > 0) then
!    To ensure positive definiteness of the thickness at each iteration, the
!  mass fluxes out of each layer are checked each step, and limited to keep
!  the thicknesses positive.  This means that several iteration may be required
!  for all the transport to happen.  The sum over domore_k keeps the processors
!  synchronized.  This may not be very efficient, but it should be reliable.

      if (x_first) then
  !    First, advect zonally.
        call advect_x(Tr, hprev, uhr, uh_neglect, domore_u, ntr, nL_max, Idt, &
                      isv, iev, jsv-stensil, jev+stensil, k, G, IG, &
                      CS%usePPM, CS%usePCM) !(, OBC)

  !    Next, advect meridionally.
        call advect_y(Tr, hprev, vhr, vh_neglect, domore_v, ntr, nL_max, Idt, &
                      isv, iev, jsv, jev, k, G, IG, CS%usePPM, CS%usePCM) !(, OBC)

        domore_k(k) = 0
        do j=jsv-stensil,jev+stensil ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
      else
  !    First, advect meridionally.
        call advect_y(Tr, hprev, vhr, vh_neglect, domore_v, ntr, nL_max, Idt, &
                      isv-stensil, iev+stensil, jsv, jev, k, G, IG, &
                      CS%usePPM, CS%usePCM) !(, OBC)

  !    Next, advect zonally.
        call advect_x(Tr, hprev, uhr, uh_neglect, domore_u, ntr, nL_max, Idt, &
                      isv, iev, jsv, jev, k, G, IG, CS%usePPM, CS%usePCM) !(, OBC)

        domore_k(k) = 0
        do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
      endif

    endif ; enddo ! End of k-loop

    ! If the advection just isn't finishing after max_iter, move on.
    if (itt >= max_iter) exit

    ! Exit if there are no layers that need more iterations.
    if (isv > is-stensil) then
      do_any = 0
      call cpu_clock_begin(id_clock_sync)
      call sum_across_PEs(domore_k(:), ncat)
      call cpu_clock_end(id_clock_sync)
      do k=1,ncat ; do_any = do_any + domore_k(k) ; enddo
      if (do_any == 0) exit
    endif

  enddo ! Iterations loop

end subroutine advect_tracer

subroutine advect_scalar(scalar, h_prev, h_end, uhtr, vhtr, dt, G, IG, CS) ! (, OBC)
  type(SIS_hor_grid_type),                      intent(inout) :: G
  type(ice_grid_type),                          intent(inout) :: IG
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(inout) :: scalar !< Scalar field to be advected
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(in)    :: h_prev, h_end
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), intent(in)    :: uhtr
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), intent(in)    :: vhtr
  real,                                         intent(in)    :: dt
  type(SIS_tracer_advect_CS),                   pointer       :: CS
! Arguments: scalar - The scalar tracer field to be advected, in arbitrary units.
!  (in)      h_prev - Category thickness times fractional coverage before advection, in m or kg m-2.
!  (in)      h_end - Layer thickness times fractional coverage after advection, in m or kg m-2.
!  (in)      uhtr - Accumulated volume or mass fluxes through zonal faces,
!                   in m3 s-1 or kg s-1.
!  (in)      vhtr - Accumulated volume or mass fluxes through meridional faces,
!                   in m3 s-1 or kg s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 SIS_tracer_advect_init.

  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)) :: &
    hprev           ! The cell volume at the end of the previous tracer
                    ! change, in m3.
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)) :: &
    uhr             ! The remaining zonal thickness flux, in m3.
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)) :: &
    vhr             ! The remaining meridional thickness fluxes, in m3.
  real :: uh_neglect(SZIB_(G),SZJ_(G)) ! uh_neglect and vh_neglect are the
  real :: vh_neglect(SZI_(G),SZJB_(G)) ! magnitude of remaining transports that
                                ! can be simply discarded, in m3 or kg.

  real :: landvolfill         ! An arbitrary? nonzero cell volume, m3.
  real :: Idt                 ! 1/dt in s-1.
  real :: h_neglect
  logical :: domore_u(SZJ_(G),SZCAT_(IG))  ! domore__ indicate whether there is more
  logical :: domore_v(SZJB_(G),SZCAT_(IG)) ! advection to be done in the corresponding
                                ! row or column.
  logical :: x_first            ! If true, advect in the x-direction first.
  integer :: max_iter           ! The maximum number of iterations in
                                ! each layer.

  real, dimension(SZIB_(G),SZJ_(G)) :: flux_U2d_x  ! x-direction tracer fluxes, in conc * kg
  real, dimension(SZI_(G),SZJB_(G)) :: flux_U2d_y  ! y-direction tracer fluxes, in conc * kg
  real    :: tr_up              ! Upwind tracer concentrations, in conc.
  real    :: vol_end, Ivol_end  ! Cell volume at the end of a step and its inverse.

  integer :: domore_k(SZCAT_(IG))
  integer :: stensil            ! The stensil of the advection scheme.
  integer :: nsten_halo         ! The number of stensils that fit in the halos.
  integer :: i, j, k, l, m, is, ie, js, je, isd, ied, jsd, jed
  integer :: ncat, nL_max, itt, do_any
  integer :: isv, iev, jsv, jev ! The valid range of the indices.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  landvolfill = 1.0e-20         ! This is arbitrary, but must be positive.
  stensil = 2                   ! The scheme's stensil; 2 for PLM.

  if (.not. associated(CS)) call SIS_error(FATAL, "advect_scalar: "// &
       "SIS_tracer_advect_init must be called before advect_scalar.")

  if (CS%use_upwind2d) then
!$OMP parallel do default(none) shared(is,ie,js,je,ncat,uhtr,vhtr,dt,G,h_end, &
!$OMP                                  h_prev,scalar) &
!$OMP                          private(tr_up,flux_U2d_x,flux_U2d_y,vol_end,Ivol_end)
    do k=1,ncat
      do j=js,je ; do I=is-1,ie
        if (uhtr(I,j,k) >= 0.0) then ; tr_up = scalar(i,j,k)
        else ; tr_up = scalar(i+1,j,k) ; endif
        flux_U2d_x(I,j) = (dt*uhtr(I,j,k)) * tr_up
      enddo ; enddo

      do J=js-1,je ; do i=is,ie
        if (vhtr(i,J,k) >= 0.0) then ; tr_up = scalar(i,j,k)
        else ; tr_up = scalar(i,j+1,k) ; endif
        flux_U2d_y(i,J) = (dt*vhtr(i,J,k)) * tr_up
      enddo ; enddo

      do j=js,je ; do i=is,ie
        vol_end = (G%areaT(i,j) * h_end(i,j,k))
        Ivol_end = 0.0 ; if (vol_end > 0.0) Ivol_end = 1.0 / vol_end
        scalar(i,j,k) = ( (G%areaT(i,j)*h_prev(i,j,k))*scalar(i,j,k) - &
                         ((flux_U2d_x(I,j) - flux_U2d_x(I-1,j)) + &
                          (flux_U2d_y(i,J) - flux_U2d_y(i,J-1))) ) * Ivol_end
      enddo ; enddo
    enddo
  else
    x_first = (MOD(G%first_direction,2) == 0)

    Idt = 1.0/dt

    max_iter = 3
    if (CS%dt > 0.0) max_iter = 2*INT(CEILING(dt/CS%dt)) + 1

    ! This initializes the halos of uhr and vhr because pass_vector might do
    ! calculations on them, even though they are never used.
    uhr(:,:,:) = 0.0 ; vhr(:,:,:) = 0.0
    hprev(:,:,:) = landvolfill
    h_neglect = IG%H_subroundoff

    ! Initialize domore_u and domore_v.  Curiously, the value used for
    ! initialization does not matter to the solutions, because if .false.
    ! they are reevaluated after the first halo update (and always on the first
    ! iteration, and if .true. a number of fluxes are exactly 0 anyway.  Of the
    ! two choices, .false. is more efficient in that it avoids extra
    ! calculations of 0 fluxes.
    domore_u(:,:) = .false. ; domore_v(:,:) = .false.

!$OMP parallel default(none) shared(is,ie,js,je,ncat,domore_k,uhr,vhr,uhtr,vhtr,dt,G, &
!$OMP                               hprev,h_prev,h_end,isd,ied,jsd,jed,uh_neglect,    &
!$OMP                               h_neglect,vh_neglect,domore_u,domore_v)
!$OMP do
    do k=1,ncat
      domore_k(k)=1

      ! Put the remaining (total) thickness fluxes into uhr and vhr.
      do j=js,je ; do I=is-1,ie ; uhr(I,j,k) = dt*uhtr(I,j,k) ; enddo ; enddo
      do J=js-1,je ; do i=is,ie ; vhr(i,J,k) = dt*vhtr(i,J,k) ; enddo ; enddo
      ! Find the previous total mass (or volume) of ice, but in the case that this
      ! category is now dramatically thinner than it was previously, add a tiny
      ! bit of extra mass to avoid nonsensical tracer concentrations.  This will
      ! lead rarely to a very slight non-conservation of tracers, but not mass.
      do i=is,ie ; do j=js,je
        hprev(i,j,k) = G%areaT(i,j) * (h_prev(i,j,k) + &
                         max(0.0, 1.0e-13*h_prev(i,j,k) - h_end(i,j,k)))
        if (h_end(i,j,k) - h_prev(i,j,k) + ((uhr(I,j,k) - uhr(I-1,j,k)) + &
                              (vhr(i,J,k) - vhr(i,J-1,k))) * G%IareaT(i,j) > &
            1e-10*(h_end(i,j,k) + h_prev(i,j,k))) then
!$OMP critical
          call SIS_error(WARNING, "Apparently inconsistent h_prev, h_end, uhr and vhr in advect_tracer.")
!$OMP end critical
        endif
      enddo ; enddo
    enddo
!$OMP end do nowait
!$OMP do
    do j=jsd,jed ; do I=isd,ied-1
      uh_neglect(I,j) = h_neglect*MIN(G%areaT(i,j),G%areaT(i+1,j))
    enddo ; enddo
!$OMP end do nowait
!$OMP do
    do J=jsd,jed-1 ; do i=isd,ied
      vh_neglect(i,J) = h_neglect*MIN(G%areaT(i,j),G%areaT(i,j+1))
    enddo ; enddo
!$OMP end parallel

    isv = is ; iev = ie ; jsv = js ; jev = je

    do itt=1,max_iter

      if (isv > is-stensil) then
        call cpu_clock_begin(id_clock_pass)
        call pass_vector(uhr, vhr, G%Domain)
        call pass_var(scalar, G%Domain, complete=.false.)
        call pass_var(hprev, G%Domain)
        call cpu_clock_end(id_clock_pass)

        nsten_halo = min(is-isd,ied-ie,js-jsd,jed-je)/stensil
        isv = is-nsten_halo*stensil ; jsv = js-nsten_halo*stensil
        iev = ie+nsten_halo*stensil ; jev = je+nsten_halo*stensil
        ! Reevaluate domore_u & domore_v unless the valid range is the same size as
        ! before.  Also, do this if there is Strang splitting.
        if ((nsten_halo > 1) .or. (itt==1)) then
!$OMP parallel do default(none) shared(isv,iev,jsv,jev,ncat,domore_k,domore_u,domore_v, &
!$OMP                                  stensil,uhr,vhr)
          do k=1,ncat ; if (domore_k(k) > 0) then
            do j=jsv,jev ; if (.not.domore_u(j,k)) then
              do i=isv+stensil-1,iev-stensil; if (uhr(I,j,k) /= 0.0) then
                domore_u(j,k) = .true. ; exit
              endif ; enddo ! i-loop
            endif ; enddo
            do J=jsv+stensil-1,jev-stensil ; if (.not.domore_v(J,k)) then
              do i=isv+stensil,iev-stensil; if (vhr(i,J,k) /= 0.0) then
                domore_v(J,k) = .true. ; exit
              endif ; enddo ! i-loop
            endif ; enddo

            !   At this point, domore_k is global.  Change it so that it indicates
            ! whether any work is needed on a layer on this processor.
            domore_k(k) = 0
            do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
            do J=jsv+stensil-1,jev-stensil ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo

          endif ; enddo ! k-loop
        endif
      endif

      ! Set the range of valid points after this iteration.
      isv = isv + stensil ; iev = iev - stensil
      jsv = jsv + stensil ; jev = jev - stensil
!$OMP parallel do default(none) shared(ncat,isv,iev,jsv,jev,x_first,domore_k,scalar, &
!$OMP                                  hprev,uhr,uh_neglect,domore_u,Idt,stensil,G,  &
!$OMP                                  CS,vhr,vh_neglect,domore_v,IG)
      do k=1,ncat ; if (domore_k(k) > 0) then
  !    To ensure positive definiteness of the thickness at each iteration, the
  !  mass fluxes out of each layer are checked each step, and limited to keep
  !  the thicknesses positive.  This means that several iteration may be required
  !  for all the transport to happen.  The sum over domore_k keeps the processors
  !  synchronized.  This may not be very efficient, but it should be reliable.

        if (x_first) then
    !    First, advect zonally.
          call advect_scalar_x(scalar, hprev, uhr, uh_neglect, domore_u, Idt, &
                        isv, iev, jsv-stensil, jev+stensil, k, G, IG, CS%usePPM, CS%usePCM) !(, OBC)

    !    Next, advect meridionally.
          call advect_scalar_y(scalar, hprev, vhr, vh_neglect, domore_v, Idt, &
                        isv, iev, jsv, jev, k, G, IG, CS%usePPM, CS%usePCM) !(, OBC)

          domore_k(k) = 0
          do j=jsv-stensil,jev+stensil ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
          do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
        else
    !    First, advect meridionally.
          call advect_scalar_y(scalar, hprev, vhr, vh_neglect, domore_v, Idt, &
                        isv-stensil, iev+stensil, jsv, jev, k, G, IG, CS%usePPM, CS%usePCM) !(, OBC)

    !    Next, advect zonally.
          call advect_scalar_x(scalar, hprev, uhr, uh_neglect, domore_u, Idt, &
                        isv, iev, jsv, jev, k, G, IG, CS%usePPM, CS%usePCM) !(, OBC)

          domore_k(k) = 0
          do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
          do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
        endif

      endif ; enddo ! End of k-loop

      ! If the advection just isn't finishing after max_iter, move on.
      if (itt >= max_iter) exit

      ! Exit if there are no layers that need more iterations.
      if (isv > is-stensil) then
        do_any = 0
        call cpu_clock_begin(id_clock_sync)
        call sum_across_PEs(domore_k(:), ncat)
        call cpu_clock_end(id_clock_sync)
        do k=1,ncat ; do_any = do_any + domore_k(k) ; enddo
        if (do_any == 0) exit
      endif

    enddo ! Iterations loop
  endif

end subroutine advect_scalar

subroutine advect_scalar_x(scalar, hprev, uhr, uh_neglect, domore_u, Idt, &
                    is, ie, js, je, k, G, IG, usePPM, usePCM) ! (, OBC)
  type(SIS_hor_grid_type),                      intent(inout) :: G
  type(ice_grid_type),                          intent(inout) :: IG
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(inout) :: scalar
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(inout) :: hprev
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), intent(inout) :: uhr
  real, dimension(SZIB_(G),SZJ_(G)),      intent(inout) :: uh_neglect
!!  type(ocean_OBC_type),                   pointer       :: OBC
  logical, dimension(SZJ_(G),SZCAT_(IG)), intent(inout) :: domore_u
  real,                                   intent(in)    :: Idt
  integer,                                intent(in)    :: is, ie, js, je, k
  logical,                                intent(in)    :: usePPM, usePCM
  !   This subroutine does 1-d flux-form advection in the zonal direction using
  ! a monotonic piecewise linear scheme.
  real, dimension(SZI_(G)) :: &
    slope_x         ! The concentration slope per grid point in units of
                    ! concentration (nondim.).
  real, dimension(SZIB_(G)) :: &
    Tr_x            ! The tracer concentration averaged over the water flux
                    ! across a zonal boundary in conc.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    mass_mask       ! A multiplicative mask at velocity points that is 1 if
                    ! both neighboring cells have any mass, and 0 otherwise.
  real :: maxslope  ! The maximum concentration slope per grid point consistent
                    ! with monotonicity, in conc. (nondim.).
  real :: hup, hlos ! hup is the upwind volume, hlos is the part of that volume
                    ! that might be lost due to advection out the other side of
                    ! the grid box, both in m3 or kg.
  real, dimension(SZIB_(G)) :: &
    uhh, &          ! The zonal flux that occurs during the current iteration, in m3 or kg.
    CFL             ! A nondimensional work variable.
  real, dimension(SZI_(G)) :: &
    hlst, Ihnew, &  ! Work variables with units of m3 or kg and m-3 or kg-1.
    haddE, haddW    ! Tiny amounts of thickness that should be added to the
                    ! tracer update with concentrations that match the average
                    ! over the fluxes through the faces to the nominal east
                    ! and west of the present cell, in m3 or kg.
  real :: hnew      ! The projected thickness, in m3 or kg.
  real :: h_add     ! A tiny thickness to add to keep the new tracer calculation
                    ! well defined in the limit of vanishing layers in m3 or kg.
  real :: I_htot    ! The inverse of the sum of thickness within or passing or
                    ! out of a cell, in m3 or kg.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  logical :: do_i(SZI_(G))  ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j
!  real :: aR, aL, dMx, dMn, Tp, Tc, Tm, dA, mA, a6
  logical :: usePLMslope

  usePLMslope = .not.(usePCM .or. usePPM)

  h_neglect = IG%H_subroundoff

  do I=is-1,ie ; CFL(I) = 0.0 ; enddo
  if (usePCM) then ; do i=is-1,ie+1 ; slope_x(i) = 0.0 ; enddo ; endif

  do j=js,je ; if (domore_u(j,k)) then
    domore_u(j,k) = .false.

    if (usePPM .or. usePLMslope) then ; do I=is-2,ie+1
      mass_mask(I,j) = 0.0
      if (G%mask2dCu(I,j) * hprev(i,j,k)*hprev(i+1,j,k) > 0.0) mass_mask(I,j) = 1.0
    enddo ; endif

    ! Calculate the i-direction profiles (slopes) of each tracer that is being advected.
    if (usePLMslope) then
      call kernel_PLM_slope_x(G, is-1, ie+1, j, scalar(:,:,k), mass_mask, slope_x(:))
    endif ! usePLMslope

    call kernel_uhh_CFL_x(G, is-1, ie, j, hprev(:,:,k), uhr(:,:,k), uhh, CFL, &
                          domore_u(j,k), h_neglect)

    if (usePPM) then
      call kernel_PPMH3_Tr_x(G, is-1, ie, j, &
             scalar(:,:,k), mass_mask, uhh, CFL, Tr_x(:))
    else ! PLM
      call kernel_PLM_Tr_x(G, is-1, ie, j, &
             scalar(:,:,k), uhh, CFL, slope_x(:), Tr_x(:))
    endif ! usePPM

    ! Calculate new tracer concentration in each cell after accounting for the i-direction fluxes.
!    call kernel_uhr_x(G, is, ie, j, uh_neglect, uhh, uhr(:,:,k), hprev(:,:,k), hlst, Ihnew, do_i)
    do I=is-1,ie
      uhr(I,j,k) = uhr(I,j,k) - uhh(I)
      if (abs(uhr(I,j,k)) < uh_neglect(I,j)) uhr(I,j,k) = 0.0
    enddo
    do i=is,ie
      if ((uhh(I) /= 0.0) .or. (uhh(I-1) /= 0.0)) then
        do_i(i) = .true.
        hlst(i) = max(hprev(i,j,k), 0.0) ! This max is here just for safety.
        hnew = hprev(i,j,k) - (uhh(I) - uhh(I-1))
        haddW(i) = 0.0 ; haddE(i) = 0.0
        if (hnew <= 0.0) then
          hnew = 0.0 ; do_i(i) = .false.
        elseif (hnew < h_neglect*G%areaT(i,j)) then
          ! Add a bit of thickness with tracer concentrations that are
          ! proportional to the mass associated with fluxes and the previous
          ! mass in the cell.
          h_add = h_neglect*G%areaT(i,j) - hnew
          I_htot = 1.0 / (hlst(i) + (abs(uhh(I)) + abs(uhh(I-1))))
          hlst(i) = hlst(i) + h_add*(hlst(i)*I_htot)
          haddW(i) = h_add * (abs(uhh(I-1))*I_htot)
          haddE(i) = h_add * (abs(uhh(I))*I_htot)

          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else
          Ihnew(i) = 1.0 / hnew
        endif
        ! Store hnew as hprev for the next iteration.
        hprev(i,j,k) = hnew
      else ! Nothing changes in this cell, so skip it.
        do_i(i) = .false.
      endif
    enddo
    do i=is,ie ; if ((do_i(i)) .and. (Ihnew(i) > 0.0)) then
      scalar(i,j,k) = (scalar(i,j,k) * hlst(i) - &
                       ((uhh(I)-haddE(i))*Tr_x(I) - &
                        (uhh(I-1)+haddW(i))*Tr_x(I-1))) * Ihnew(i)
    endif ; enddo
!    call kernel_tracer_div_x(G, is, ie, j, do_i, hlst, Ihnew, flux_x(:), scalar(:,:,k))

  endif ; enddo ! End of j-loop.

end subroutine advect_scalar_x

subroutine advect_x(Tr, hprev, uhr, uh_neglect, domore_u, ntr, nL_max, Idt, &
                    is, ie, js, je, k, G, IG, usePPM, usePCM) ! (, OBC)
  type(SIS_hor_grid_type),                      intent(inout) :: G
  type(ice_grid_type),                          intent(inout) :: IG
  type(SIS_tracer_type), dimension(ntr),        intent(inout) :: Tr
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(inout) :: hprev
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), intent(inout) :: uhr
  real, dimension(SZIB_(G),SZJ_(G)),      intent(inout) :: uh_neglect
!!  type(ocean_OBC_type),                   pointer       :: OBC
  logical, dimension(SZJ_(G),SZCAT_(IG)), intent(inout) :: domore_u
  real,                                   intent(in)    :: Idt
  integer,                                intent(in)    :: ntr, nL_max, is, ie, js, je,k
  logical,                                intent(in)    :: usePPM, usePCM
  !   This subroutine does 1-d flux-form advection in the zonal direction using
  ! a monotonic piecewise linear scheme.
  real, dimension(SZI_(G),nL_max,ntr) :: &
    slope_x         ! The concentration slope per grid point in units of
                    ! concentration (nondim.).
  real, dimension(SZIB_(G),nL_max,ntr) :: &
    Tr_x            ! The tracer concentration averaged over the water flux
                    ! across a zonal boundary in conc.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    mass_mask       ! A multiplicative mask at velocity points that is 1 if
                    ! both neighboring cells have any mass, and 0 otherwise.
  real :: maxslope  ! The maximum concentration slope per grid point consistent
                    ! with monotonicity, in conc. (nondim.).
  real :: hup, hlos ! hup is the upwind volume, hlos is the part of that volume
                    ! that might be lost due to advection out the other side of
                    ! the grid box, both in m3 or kg.
  real, dimension(SZIB_(G)) :: &
    uhh, &          ! The zonal flux that occurs during the current iteration, in m3 or kg.
    CFL             ! A nondimensional work variable.
  real, dimension(SZI_(G)) :: &
    hlst, Ihnew, &  ! Work variables with units of m3 or kg and m-3 or kg-1.
    haddE, haddW    ! Tiny amounts of thickness that should be added to the
                    ! tracer update with concentrations that match the average
                    ! over the fluxes through the faces to the nominal east
                    ! and west of the present cell, in m3 or kg.
  real :: hnew      ! The projected thickness, in m3 or kg.
  real :: h_add     ! A tiny thickness to add to keep the new tracer calculation
                    ! well defined in the limit of vanishing layers in m3 or kg.
  real :: I_htot    ! The inverse of the sum of thickness within or passing or
                    ! out of a cell, in m3 or kg.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  logical :: do_i(SZI_(G))  ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, l, m
  logical :: usePLMslope

  usePLMslope = .not.(usePCM .or. usePPM)

  h_neglect = IG%H_subroundoff

  do I=is-1,ie ; CFL(I) = 0.0 ; enddo

  if (usePCM) then ; do m=1,ntr ; do l=1,Tr(m)%nL ; do i=is-1,ie+1
    slope_x(i,l,m) = 0.0
  enddo ; enddo ; enddo ; endif

  do j=js,je ; if (domore_u(j,k)) then
    domore_u(j,k) = .false.

    if (usePPM .or. usePLMslope) then ; do I=is-2,ie+1
      mass_mask(I,j) = 0.0
      if (G%mask2dCu(I,j)*hprev(i,j,k)*hprev(i+1,j,k) > 0.0) mass_mask(I,j) = 1.0
    enddo ; endif

    ! Calculate the i-direction profiles (slopes) of each tracer that is being advected.
    if (usePLMslope) then
      do m=1,ntr ; do l=1,Tr(m)%nL
        call kernel_PLM_slope_x(G, is-1, ie+1, j, Tr(m)%t(:,:,k,l), mass_mask, slope_x(:,l,m))
      enddo ; enddo
    endif ! usePLMslope

    call kernel_uhh_CFL_x(G, is-1, ie, j, hprev(:,:,k), uhr(:,:,k), uhh, CFL, &
                          domore_u(j,k), h_neglect)

    if (usePPM) then
      do m=1,ntr ; do l=1,Tr(m)%nL
        call kernel_PPMH3_Tr_x(G, is-1, ie, j, &
               Tr(m)%t(:,:,k,l), mass_mask, uhh, CFL, Tr_x(:,l,m))
      enddo ; enddo
    else ! PLM
      do m=1,ntr ; do l=1,Tr(m)%nL
        call kernel_PLM_Tr_x(G, is-1, ie, j, &
               Tr(m)%t(:,:,k,l), uhh, CFL, slope_x(:,l,m), Tr_x(:,l,m))
      enddo ; enddo
    endif ! usePPM

    ! Calculate new tracer concentration in each cell after accounting for the i-direction fluxes.
!    call kernel_uhr_x(G, is, ie, j, uh_neglect, uhh, uhr(:,:,k), hprev(:,:,k), hlst, Ihnew, do_i)
    do I=is-1,ie
      uhr(I,j,k) = uhr(I,j,k) - uhh(I)
      if (abs(uhr(I,j,k)) < uh_neglect(I,j)) uhr(I,j,k) = 0.0
    enddo
    do i=is,ie
      if ((uhh(I) /= 0.0) .or. (uhh(I-1) /= 0.0)) then
        do_i(i) = .true.
        hlst(i) = max(hprev(i,j,k), 0.0) ! This max is here just for safety.
        hnew = hprev(i,j,k) - (uhh(I) - uhh(I-1))
        haddW(i) = 0.0 ; haddE(i) = 0.0
        if (hnew <= 0.0) then
          hnew = 0.0 ; do_i(i) = .false.
        elseif (hnew < h_neglect*G%areaT(i,j)) then
          ! Add a bit of thickness with tracer concentrations that are
          ! proportional to the mass associated with fluxes and the previous
          ! mass in the cell.
          h_add = h_neglect*G%areaT(i,j) - hnew
          I_htot = 1.0 / (hlst(i) + (abs(uhh(I)) + abs(uhh(I-1))))
          hlst(i) = hlst(i) + h_add*(hlst(i)*I_htot)
          haddW(i) = h_add * (abs(uhh(I-1))*I_htot)
          haddE(i) = h_add * (abs(uhh(I))*I_htot)

          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else
          Ihnew(i) = 1.0 / hnew
        endif
        ! Store hnew as hprev for the next iteration.
        hprev(i,j,k) = hnew
      else ! Nothing changes in this cell, so skip it.
        do_i(i) = .false.
      endif
    enddo
    do m=1,ntr ; do l=1,Tr(m)%nL
!       call kernel_tracer_div_x(G, is, ie, j, do_i, hlst, Ihnew, flux_x(:,l,m), Tr(m)%t(:,:,k,l))
      do i=is,ie ; if ((do_i(i)) .and. (Ihnew(i) > 0.0)) then
        Tr(m)%t(i,j,k,l) = (Tr(m)%t(i,j,k,l) * hlst(i) - &
                            ((uhh(I)-haddE(i))*Tr_x(I,l,m) - &
                             (uhh(I-1)+haddW(i))*Tr_x(I-1,l,m))) * Ihnew(i)
      endif ; enddo
      ! Diagnostics
      if (associated(Tr(m)%ad4d_x)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad4d_x(I,j,k,l) = Tr(m)%ad4d_x(I,j,k,l) + uhh(I)*Tr_x(I,l,m)*Idt
      endif ; enddo ; endif
      if (associated(Tr(m)%ad3d_x)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad3d_x(I,j,k) = Tr(m)%ad3d_x(I,j,k) + uhh(I)*Tr_x(I,l,m)*Idt
      endif ; enddo ; endif
      if (associated(Tr(m)%ad2d_x)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad2d_x(I,j) = Tr(m)%ad2d_x(I,j) + uhh(I)*Tr_x(I,l,m)*Idt
      endif ; enddo ; endif
    enddo ; enddo

  endif ; enddo ! End of j-loop.

end subroutine advect_x

!>  Calculate the mass flux and CFL such that the flux of tracer uses as much
!! the minimum of the remaining mass flux (uhr) and the half the mass
!! in the cell plus whatever part of its half of the mass flux that
!! the flux through the other side does not require.
subroutine kernel_uhh_CFL_x(G, is, ie, j, hprev, uhr, uhh, CFL, domore_u, h_neglect)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: hprev
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: uhr
  real, dimension(SZIB_(G)),         intent(inout) :: uhh, CFL
  logical,                           intent(inout) :: domore_u
  real,                              intent(in)    :: h_neglect
  ! Local
  integer :: i
  real :: hup, hlos

  do I=is,ie
    if (uhr(I,j) == 0.0) then
      uhh(I) = 0.0
      CFL(I) = 0.0
    elseif (uhr(I,j) < 0.0) then
      hup = hprev(i+1,j)
      hlos = MAX(0.0,uhr(I+1,j))
      if (((hup + uhr(I,j) - hlos) < 0.0) .and. &
          ((0.5*hup + uhr(I,j)) < 0.0)) then
        uhh(I) = MIN(-0.5*hup,-hup+hlos,0.0)
        domore_u = .true.
      else
        uhh(I) = uhr(I,j)
      endif
      CFL(I) = - uhh(I)/(hprev(i+1,j)+h_neglect) ! CFL is positive
    else
      hup = hprev(i,j)
      hlos = MAX(0.0,-uhr(I-1,j))
      if (((hup - uhr(I,j) - hlos) < 0.0) .and. &
          ((0.5*hup - uhr(I,j)) < 0.0)) then
        uhh(I) = MAX(0.5*hup,hup-hlos,0.0)
        domore_u = .true.
      else
        uhh(I) = uhr(I,j)
      endif
      CFL(I) = uhh(I)/(hprev(i,j)+h_neglect) ! CFL is positive
    endif
  enddo

end subroutine kernel_uhh_CFL_x

subroutine kernel_PLM_slope_x(G, is, ie, j, scalar, uMask, slope_x)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: scalar
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: uMask
  real, dimension(SZI_(G)),          intent(inout) :: slope_x
  ! Local
  integer :: i
  real :: Tp, Tc, Tm, dMx, dMn

  do i = is, ie
    Tp = scalar(i+1,j) ; Tc = scalar(i,j) ; Tm = scalar(i-1,j)
    dMx = max( Tp, Tc, Tm ) - Tc
    dMn= Tc - min( Tp, Tc, Tm )
    slope_x(i) = uMask(I,j)*uMask(I-1,j) * &
        sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
  enddo

end subroutine kernel_PLM_slope_x

subroutine kernel_PLM_Tr_x(G, is, ie, j, scalar, uhh, CFL, slope_x, Tr_x)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: scalar
  real, dimension(SZIB_(G)),         intent(in)    :: uhh, CFL
  real, dimension(SZI_(G)),          intent(in)    :: slope_x
  real, dimension(SZIB_(G)),         intent(inout) :: Tr_x
  ! Local
  integer :: i
  real :: Tc

  do I=is,ie
    if (uhh(I) >= 0.0) then
      Tc = scalar(i,j)
      Tr_x(I) = Tc + 0.5 * slope_x(i) * ( 1. - CFL(I) )
    else
      Tc = scalar(i+1,j)
      Tr_x(I) = Tc - 0.5 * slope_x(i+1) * ( 1. - CFL(I) )
    endif
  enddo

end subroutine kernel_PLM_Tr_x

subroutine kernel_PPMH3_Tr_x(G, is, ie, j, scalar, uMask, uhh, CFL, Tr_x)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: scalar
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: uMask
  real, dimension(SZIB_(G)),         intent(in)    :: uhh, CFL
  real, dimension(SZIB_(G)),         intent(inout) :: Tr_x
  ! Local
  integer :: i
  real :: Tp, Tc, Tm, aL, aR, dA, a6, mA

  do I=is,ie
    if (uhh(I) >= 0.0) then
      ! Implementation of PPM-H3
      Tp = scalar(i+1,j) ; Tc = scalar(i,j) ; Tm = scalar(i-1,j)
      aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
      aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
      aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
      aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
      dA = aR - aL ; mA = 0.5*( aR + aL )

      ! These expressions are uglier than they might be, but they are less
      ! sensitive to underflow than the alternatives would be.
      if ((uMask(I,j)*uMask(I-1,j) == 0.0) .or. (Tp == Tc) .or. (Tc == Tm) .or. &
          (sign(1.,Tp-Tc)*sign(1.,Tc-Tm) <= 0.)) then
        aL = Tc ; aR = Tc ! PCM for local extrema and boundary cells
        a6 = 0.0 ! Curvature
      elseif ( 6.*sign(1.,dA)*(Tc-mA) > abs(dA) ) then
        ! aL = 3.*Tc - 2.*aR
        aL = Tc + 2.*(Tc - aR)
        a6 = 3.*(aR - Tc) ! Curvature
      elseif ( -6.*sign(1.,dA)*(Tc-mA) > abs(dA) ) then
        ! aR = 3.*Tc - 2.*aL
        aR = Tc + 2.*(Tc - aL)
        a6 = 3.*(aL - Tc) ! Curvature
      else
        a6 = 3.*((Tc - aR) + (Tc - aL)) ! Curvature
      endif
      ! a6 = 6.*Tc - 3. * (aR + aL) ! Curvature

      Tr_x(I) = ( aR - 0.5 * CFL(I) * ( &
              ( aR - aL ) - a6 * ( 1. - 2./3. * CFL(I) ) ) )
    else
      ! Implementation of PPM-H3
      Tp = scalar(i+2,j) ; Tc = scalar(i+1,j) ; Tm = scalar(i,j)
      aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
      aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
      aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
      aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
      dA = aR - aL ; mA = 0.5*( aR + aL )

      if ((uMask(I,j)*uMask(I+1,j) == 0.0) .or. (Tp == Tc) .or. (Tc == Tm) .or. &
          (sign(1.,Tp-Tc)*sign(1.,Tc-Tm) <= 0.)) then
        aL = Tc ; aR = Tc ! PCM for local extrema and boundary cells
        a6 = 0.0 ! Curvature
      elseif ( 6.*sign(1.,dA)*(Tc-mA) > abs(dA) ) then
        ! aL = 3.*Tc - 2.*aR
        aL = Tc + 2.*(Tc - aR)
        a6 = 3.*(aR - Tc) ! Curvature
      elseif ( -6.*sign(1.,dA)*(Tc-mA) > abs(dA) ) then
        ! aR = 3.*Tc - 2.*aL
        aR = Tc + 2.*(Tc - aL)
        a6 = 3.*(aL - Tc) ! Curvature
      else
        a6 = 3.*((Tc - aR) + (Tc - aL)) ! Curvature
      endif
      ! a6 = 6.*Tc - 3. * (aR + aL) ! Curvature

      Tr_x(I) = ( aL + 0.5 * CFL(I) * ( &
              ( aR - aL ) + a6 * ( 1. - 2./3. * CFL(I) ) ) )
    endif
  enddo

end subroutine kernel_PPMH3_Tr_x

subroutine kernel_uhr_x(G, is, ie, j, uh_neglect, uhh, uhr, hprev, hlst, Ihnew, do_i, h_neglect)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: uh_neglect
  real, dimension(SZIB_(G)),         intent(in)    :: uhh
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: uhr
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: hprev
  real, dimension(SZI_(G)),          intent(inout) :: hlst, Ihnew
  logical, dimension(SZI_(G)),       intent(inout) :: do_i
  real,                              intent(in)    :: h_neglect
  ! Local
  integer :: i

  do I=is-1,ie
    uhr(I,j) = uhr(I,j) - uhh(I)
    if (abs(uhr(I,j)) < uh_neglect(I,j)) uhr(I,j) = 0.0
  enddo
  do i=is,ie
    if ((uhh(I) /= 0.0) .or. (uhh(I-1) /= 0.0)) then
      do_i(i) = .true.
      hlst(i) = hprev(i,j)
      hprev(i,j) = hprev(i,j) - (uhh(I) - uhh(I-1))
      if (hprev(i,j) <= 0.0) then
        do_i(i) = .false.
      elseif (hprev(i,j) < h_neglect*G%areaT(i,j)) then
        hlst(i) = hlst(i) + (h_neglect*G%areaT(i,j) - hprev(i,j))
        Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
      else
        Ihnew(i) = 1.0 / hprev(i,j)
      endif
    else
      do_i(i) = .false.
    endif
  enddo

end subroutine kernel_uhr_x

!> Updates a scalar with the divergence of x-flux
subroutine kernel_tracer_div_x(G, is, ie, j, do_i, hlst, Ihnew, flux_x, scalar)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  logical, dimension(SZI_(G)),       intent(in)    :: do_i
  real, dimension(SZI_(G)),          intent(in)    :: hlst, Ihnew
  real, dimension(SZIB_(G)),         intent(in)    :: flux_x
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: scalar
  ! Local
  integer :: i

  do i=is,ie ; if ((do_i(i)) .and. (Ihnew(i) > 0.0)) then
    scalar(i,j) = (scalar(i,j) * hlst(i) - &
                      (flux_x(I) - flux_x(I-1))) * Ihnew(i)
  endif ; enddo

end subroutine kernel_tracer_div_x

subroutine advect_scalar_y(scalar, hprev, vhr, vh_neglect, domore_v, Idt, &
                    is, ie, js, je, k, G, IG, usePPM, usePCM) ! (, OBC)
  type(SIS_hor_grid_type),                      intent(inout) :: G
  type(ice_grid_type),                          intent(inout) :: IG
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(inout) :: scalar
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(inout) :: hprev
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), intent(inout) :: vhr
  real, dimension(SZI_(G),SZJB_(G)),            intent(inout) :: vh_neglect
!  type(ocean_OBC_type),                        pointer       :: OBC
  logical, dimension(SZJB_(G),SZCAT_(IG)),      intent(inout) :: domore_v
  real,                                         intent(in)    :: Idt
  integer,                                      intent(in)    :: is, ie, js, je, k
  logical,                                      intent(in)    :: usePPM, usePCM
  !   This subroutine does 1-d flux-form advection using a monotonic piecewise
  ! linear scheme.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    slope_y         ! The concentration slope per grid point in units of
                    ! concentration (nondim.).
  real, dimension(SZI_(G),SZJB_(G)) :: &
    Tr_y            ! The tracer concentration averaged over the water flux
                    ! across a meridional boundary in conc.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    mass_mask, &    ! A multiplicative mask at velocity points that is 1 if
                    ! both neighboring cells have any mass, and 0 otherwise.
    vhh             ! The meridional flux that occurs during the current
                    ! iteration, in m3 or kg.
  real :: maxslope  ! The maximum concentration slope per grid point consistent
                    ! with monotonicity, in conc. (nondim.).
  real :: hup, hlos ! hup is the upwind volume, hlos is the part of that volume
                    ! that might be lost due to advection out the other side of
                    ! the grid box, both in m3 or kg.
  real, dimension(SZI_(G)) :: &
    hlst, Ihnew, &  ! Work variables with units of m3 or kg and m-3 or kg-1.
    haddN, haddS, & ! Tiny amounts of thickness that should be added to the
                    ! tracer update with concentrations that match the average
                    ! over the fluxes through the faces to the nominal north
                    ! and south of the present cell, in m3 or kg.
    CFL             ! A nondimensional work variable.
  real :: hnew      ! The projected thickness, in m3 or kg.
  real :: h_add     ! A tiny thickness to add to keep the new tracer calculation
                    ! well defined in the limit of vanishing layers in m3 or kg.
  real :: I_htot    ! The inverse of the sum of thickness within or passing or
                    ! out of a cell, in m3 or kg.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  logical :: do_j_tr(SZJ_(G))  ! If true, calculate the tracer profiles.
  logical :: do_i(SZI_(G))     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, l, m
!  real :: aR, aL, dMx, dMn, Tp, Tc, Tm, dA, mA, a6
  logical :: usePLMslope

  usePLMslope = .not.(usePCM .or. usePPM)

  h_neglect = IG%H_subroundoff

  do_j_tr(js-1) = domore_v(js-1,k) ; do_j_tr(je+1) = domore_v(je,k)
  do j=js,je ; do_j_tr(j) = (domore_v(J-1,k) .or. domore_v(J,k)) ; enddo

  if (usePPM .or. usePLMslope) then ; do J=js-2,je+1 ; do i=is,ie
    mass_mask(i,J) = 0.0
    if (G%mask2dCv(i,J)*hprev(i,j,k)*hprev(i,j+1,k) > 0.0) mass_mask(i,J) = 1.0
  enddo ; enddo ; endif

  ! Calculate the j-direction profiles (slopes) of each tracer that is being advected.
  if (usePLMslope) then
    do j=js-1,je+1 ; if (do_j_tr(j)) then
      call kernel_PLM_slope_y(G, is, ie, j, scalar(:,:,k), mass_mask, slope_y(:,j))
    endif ; enddo
  elseif (usePCM) then
    do j=js-1,je+1 ; do i=is,ie ; slope_y(i,j) = 0.0 ; enddo ; enddo
  endif ! usePLMslope

  do J=js-1,je ; if (domore_v(J,k)) then
    call kernel_vhh_CFL_y(G, is, ie, J, hprev(:,:,k), vhr(:,:,k), vhh, CFL, &
                          domore_v(:,k), h_neglect)
    if (usePPM) then
      call kernel_PPMH3_Tr_y(G, is, ie, J, &
             scalar(:,:,k), mass_mask, vhh, CFL, Tr_y(:,J))
    else ! PLM
      call kernel_PLM_Tr_y(G, is, ie, J, &
             scalar(:,:,k), vhh, CFL, slope_y(:,:), Tr_y(:,J))
    endif ! usePPM

  else ! not domore_v.
    do i=is,ie ; vhh(i,J) = 0.0 ; Tr_y(i,J) = 0.0 ; enddo
  endif ; enddo ! End of j-loop

  do J=js-1,je ; do i=is,ie
    vhr(i,J,k) = vhr(i,J,k) - vhh(i,J)
    if (abs(vhr(i,J,k)) < vh_neglect(i,J)) vhr(i,J,k) = 0.0
  enddo ; enddo

  ! Calculate new tracer concentration in each cell after accounting for the j-direction fluxes.
  do j=js,je ; if (do_j_tr(j)) then
!    call kernel_hlst_y(G, is, ie, j, vh_neglect, vhh, hprev(:,:,k), hlst, Ihnew, do_i, h_neglect)
!    call kernel_tracer_div_y(G, is, ie, j, do_i, hlst, Ihnew, flux_y(:,:), scalar(:,:,k))
    do i=is,ie
      if ((vhh(i,J) /= 0.0) .or. (vhh(i,J-1) /= 0.0)) then
        do_i(i) = .true.
        hlst(i) = max(hprev(i,j,k), 0.0) ! This max is here just for safety.
        hnew = hprev(i,j,k) - (vhh(i,J) - vhh(i,J-1))
        haddS(i) = 0.0 ; haddN(i) = 0.0
        if (hnew <= 0.0) then
          hnew = 0.0 ; do_i(i) = .false.
        elseif (hnew < h_neglect*G%areaT(i,j)) then
          ! Add a tiny bit of thickness with tracer concentrations that are
          ! proportional to the mass associated with fluxes and the previous
          ! mass in the cell.
          h_add = h_neglect*G%areaT(i,j) - hnew
          I_htot = 1.0 / (hlst(i) + (abs(vhh(i,J)) + abs(vhh(i,J-1))))
          hlst(i) = hlst(i) + h_add*(hlst(i)*I_htot)
          haddS(i) = h_add * (abs(vhh(i,J-1))*I_htot)
          haddN(i) = h_add * (abs(vhh(i,J))*I_htot)

          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else
          Ihnew(i) = 1.0 / hnew
        endif
        ! Store hnew as hprev for the next iteration.
        hprev(i,j,k) = hnew
      else ! Nothing changes in this cell, so skip it.
        do_i(i) = .false.
      endif
    enddo
    do i=is,ie ; if (do_i(i)) then
      scalar(i,j,k) = (scalar(i,j,k) * hlst(i) - &
                       ((vhh(i,J)-haddN(i))*Tr_y(i,J) - &
                        (vhh(i,J-1)+haddS(i))*Tr_y(i,J-1))) * Ihnew(i)
    endif ; enddo
  endif ; enddo ! End of j-loop.

end subroutine advect_scalar_y

subroutine advect_y(Tr, hprev, vhr, vh_neglect, domore_v, ntr, nL_max, Idt, &
                    is, ie, js, je, k, G, IG, usePPM, usePCM) ! (, OBC)
  type(SIS_hor_grid_type),                      intent(inout) :: G
  type(ice_grid_type),                          intent(inout) :: IG
  type(SIS_tracer_type), dimension(ntr),        intent(inout) :: Tr
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(inout) :: hprev
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), intent(inout) :: vhr
  real, dimension(SZI_(G),SZJB_(G)),            intent(inout) :: vh_neglect
!  type(ocean_OBC_type),                        pointer       :: OBC
  logical, dimension(SZJB_(G),SZCAT_(IG)),      intent(inout) :: domore_v
  real,                                         intent(in)    :: Idt
  integer,                                      intent(in)    :: ntr, nL_max, is, ie, js, je, k
  logical,                                      intent(in)    :: usePPM, usePCM
  !   This subroutine does 1-d flux-form advection using a monotonic piecewise
  ! linear scheme.
  real, dimension(SZI_(G),SZJ_(G),nL_max,ntr) :: &
    slope_y         ! The concentration slope per grid point in units of
                    ! concentration (nondim.).
  real, dimension(SZI_(G),SZJB_(G),nL_max,ntr) :: &
    Tr_y            ! The tracer concentration averaged over the water flux
                    ! across a meridional boundary in conc.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    mass_mask, &    ! A multiplicative mask at velocity points that is 1 if
                    ! both neighboring cells have any mass, and 0 otherwise.
    vhh             ! The meridional flux that occurs during the current
                    ! iteration, in m3 or kg.
  real :: maxslope  ! The maximum concentration slope per grid point consistent
                    ! with monotonicity, in conc. (nondim.).
  real :: hup, hlos ! hup is the upwind volume, hlos is the part of that volume
                    ! that might be lost due to advection out the other side of
                    ! the grid box, both in m3 or kg.
  real, dimension(SZI_(G)) :: &
    hlst, Ihnew, &  ! Work variables with units of m3 or kg and m-3 or kg-1.
    haddN, haddS, & ! Tiny amounts of thickness that should be added to the
                    ! tracer update with concentrations that match the average
                    ! over the fluxes through the faces to the nominal north
                    ! and south of the present cell, in m3 or kg.
    CFL             ! A nondimensional work variable.
  real :: hnew      ! The projected thickness, in m3 or kg.
  real :: h_add     ! A tiny thickness to add to keep the new tracer calculation
                    ! well defined in the limit of vanishing layers in m3 or kg.
  real :: I_htot    ! The inverse of the sum of thickness within or passing or
                    ! out of a cell, in m3 or kg.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in m.
  logical :: do_j_tr(SZJ_(G))  ! If true, calculate the tracer profiles.
  logical :: do_i(SZI_(G))     ! If true, work on given points.
  logical :: do_any_i
  logical :: usePLMslope
  integer :: i, j, l, m

  usePLMslope = .not.(usePCM .or. usePPM)

  h_neglect = IG%H_subroundoff

  do_j_tr(js-1) = domore_v(js-1,k) ; do_j_tr(je+1) = domore_v(je,k)
  do j=js,je ; do_j_tr(j) = (domore_v(J-1,k) .or. domore_v(J,k)) ; enddo

  if (usePPM .or. usePLMslope) then ; do J=js-2,je+1 ; do i=is,ie
    mass_mask(i,J) = 0.0
    if (G%mask2dCv(i,J)*hprev(i,j,k)*hprev(i,j+1,k) > 0.0) mass_mask(i,J) = 1.0
  enddo ; enddo ; endif

  ! Calculate the j-direction profiles (slopes) of each tracer that is being advected.
  if (usePLMslope) then
    do j=js-1,je+1 ; if (do_j_tr(j)) then ; do m=1,ntr ; do l=1,Tr(m)%nL
      call kernel_PLM_slope_y(G, is, ie, j, Tr(m)%t(:,:,k,l), mass_mask, slope_y(:,j,l,m))
    enddo ; enddo ; endif ; enddo ! End of l-, m-, & j- loops.
  elseif (usePCM) then
    do m=1,ntr ; do l=1,Tr(m)%nL ; do j=js-1,je+1 ; do i=is,ie
      slope_y(i,j,l,m) = 0.0
    enddo ; enddo ; enddo ; enddo
  endif ! usePLMslope

  do J=js-1,je ; if (domore_v(J,k)) then
    call kernel_vhh_CFL_y(G, is, ie, J, hprev(:,:,k), vhr(:,:,k), vhh, CFL, &
                          domore_v(:,k), h_neglect)
    if (usePPM) then
      do m=1,ntr ; do l=1,Tr(m)%nL
        call kernel_PPMH3_Tr_y(G, is, ie, J, &
               Tr(m)%t(:,:,k,l), mass_mask, vhh, CFL, Tr_y(:,J,l,m))
      enddo ; enddo
    else ! PLM
      do m=1,ntr ; do l=1,Tr(m)%nL
        call kernel_PLM_Tr_y(G, is, ie, J, &
               Tr(m)%t(:,:,k,l), vhh, CFL, slope_y(:,:,l,m), Tr_y(:,J,l,m))
      enddo ; enddo
    endif ! usePPM

  else ! not domore_v.
    do i=is,ie ; vhh(i,J) = 0.0 ; enddo
    do m=1,ntr ; do l=1,Tr(m)%nL ; do i=is,ie ; Tr_y(i,J,l,m) = 0.0 ; enddo ; enddo ; enddo
  endif ; enddo ! End of j-loop

  do J=js-1,je ; do i=is,ie
    vhr(i,J,k) = vhr(i,J,k) - vhh(i,J)
    if (abs(vhr(i,J,k)) < vh_neglect(i,J)) vhr(i,J,k) = 0.0
  enddo ; enddo

  ! Calculate new tracer concentration in each cell after accounting for the j-direction fluxes.
  do j=js,je ; if (do_j_tr(j)) then
!    call kernel_hlst_y(G, is, ie, j, vh_neglect, vhh, hprev(:,:,k), hlst, Ihnew, do_i, h_neglect)
    do i=is,ie
      if ((vhh(i,J) /= 0.0) .or. (vhh(i,J-1) /= 0.0)) then
        do_i(i) = .true.
        hlst(i) = max(hprev(i,j,k), 0.0) ! This max is here just for safety.
        hnew = hprev(i,j,k) - (vhh(i,J) - vhh(i,J-1))
        haddS(i) = 0.0 ; haddN(i) = 0.0
        if (hnew <= 0.0) then
          hnew = 0.0 ; do_i(i) = .false.
        elseif (hnew < h_neglect*G%areaT(i,j)) then
          ! Add a tiny bit of thickness with tracer concentrations that are
          ! proportional to the mass associated with fluxes and the previous
          ! mass in the cell.
          h_add = h_neglect*G%areaT(i,j) - hnew
          I_htot = 1.0 / (hlst(i) + (abs(vhh(i,J)) + abs(vhh(i,J-1))))
          hlst(i) = hlst(i) + h_add*(hlst(i)*I_htot)
          haddS(i) = h_add * (abs(vhh(i,J-1))*I_htot)
          haddN(i) = h_add * (abs(vhh(i,J))*I_htot)

          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else
          Ihnew(i) = 1.0 / hnew
        endif
        ! Store hnew as hprev for the next iteration.
        hprev(i,j,k) = hnew
      else ! Nothing changes in this cell, so skip it.
        do_i(i) = .false.
      endif
    enddo
    do m=1,ntr ; do l=1,Tr(m)%nL
!      call kernel_tracer_div_y(G, is, ie, j, do_i, hlst, Ihnew, flux_y(:,:,l,m), Tr(m)%t(:,:,k,l))
      do i=is,ie ; if (do_i(i)) then
        Tr(m)%t(i,j,k,l) = (Tr(m)%t(i,j,k,l) * hlst(i) - &
               ((vhh(i,J)-haddN(i))*Tr_y(i,J,l,m) - &
                (vhh(i,J-1)+haddS(i))*Tr_y(i,J-1,l,m))) * Ihnew(i)
      endif ; enddo
      ! Diagnostics
      if (associated(Tr(m)%ad4d_y)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad4d_y(i,J,k,l) = Tr(m)%ad4d_y(i,J,k,l) + vhh(i,J)*Tr_y(i,J,l,m)*Idt
      endif ; enddo ; endif
      if (associated(Tr(m)%ad3d_y)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad3d_y(i,J,k) = Tr(m)%ad3d_y(i,J,k) + vhh(i,J)*Tr_y(i,J,l,m)*Idt
      endif ; enddo ; endif
      if (associated(Tr(m)%ad2d_y)) then ; do i=is,ie ; if (do_i(i)) then
        Tr(m)%ad2d_y(i,J) = Tr(m)%ad2d_y(i,J) + vhh(i,J)*Tr_y(i,J,l,m)*Idt
      endif ; enddo ; endif
    enddo ; enddo
  endif ; enddo ! End of j-loop.
  ! Diagnostics (on southern edge?)
  J = js-1
  do m=1,ntr
    if (associated(Tr(m)%ad4d_y) .or. associated(Tr(m)%ad3d_y) .or. associated(Tr(m)%ad3d_y)) then
      if (associated(Tr(m)%ad4d_y)) then ; do l=1,Tr(m)%nL ; do i=is,ie
        Tr(m)%ad4d_y(i,J,k,l) = Tr(m)%ad4d_y(i,J,k,l) + vhh(i,J)*Tr_y(i,J,l,m)*Idt
      enddo ; enddo ; endif
      if (associated(Tr(m)%ad3d_y)) then ; do l=1,Tr(m)%nL ; do i=is,ie
        Tr(m)%ad3d_y(i,J,k) = Tr(m)%ad3d_y(i,J,k) + vhh(i,J)*Tr_y(i,J,l,m)*Idt
      enddo ; enddo ; endif
      if (associated(Tr(m)%ad2d_y)) then ; do l=1,Tr(m)%nL ; do i=is,ie
        Tr(m)%ad2d_y(i,J) = Tr(m)%ad2d_y(i,J) + vhh(i,J)*Tr_y(i,J,l,m)*Idt
      enddo ; enddo ; endif
    endif
  enddo ! m

end subroutine advect_y

!>  Calculate the mass flux and CFL such that the flux of tracer uses as much
!! the minimum of the remaining mass flux (vhr) and the half the mass
!! in the cell plus whatever part of its half of the mass flux that
!! the flux through the other side does not require.
subroutine kernel_vhh_CFL_y(G, is, ie, J, hprev, vhr, vhh, CFL, domore_v, h_neglect)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, J
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: hprev
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: vhr
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: vhh
  real, dimension(SZI_(G)),          intent(inout) :: CFL
  logical, dimension(SZJB_(G)),      intent(inout) :: domore_v
  real,                              intent(in)    :: h_neglect
  ! Local
  integer :: i
  real :: hup, hlos

  domore_v(J) = .false.
  do i=is,ie
    if (vhr(i,J) == 0.0) then
      vhh(i,J) = 0.0
      CFL(i) = 0.0
    elseif (vhr(i,J) < 0.0) then
      hup = hprev(i,j+1)
      hlos = MAX(0.0,vhr(i,J+1))
      if ((((hup - hlos) + vhr(i,J)) < 0.0) .and. &
          ((0.5*hup + vhr(i,J)) < 0.0)) then
        vhh(i,J) = MIN(-0.5*hup,-hup+hlos,0.0)
        domore_v(J) = .true.
      else
        vhh(i,J) = vhr(i,J)
      endif
      CFL(i) = - vhh(i,J) / (hprev(i,j+1)+h_neglect) ! CFL is positive
    else
      hup = hprev(i,j)
      hlos = MAX(0.0,-vhr(i,J-1))
      if ((((hup - hlos) - vhr(i,J)) < 0.0) .and. &
          ((0.5*hup - vhr(i,J)) < 0.0)) then
        vhh(i,J) = MAX(0.5*hup,hup-hlos,0.0)
        domore_v(J) = .true.
      else
        vhh(i,J) = vhr(i,J)
      endif
      CFL(i) = vhh(i,J) / (hprev(i,j)+h_neglect) ! CFL is positive
    endif
  enddo

end subroutine kernel_vhh_CFL_y

subroutine kernel_PLM_slope_y(G, is, ie, j, scalar, vMask, slope_y)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: scalar
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: vMask
  real, dimension(SZI_(G)),          intent(inout) :: slope_y
  ! Local
  integer :: i
  real :: Tp, Tc, Tm, dMx, dMn

  do i = is, ie
    Tp = scalar(i,j+1) ; Tc = scalar(i,j) ; Tm = scalar(i,j-1)
    dMx = max( Tp, Tc, Tm ) - Tc
    dMn= Tc - min( Tp, Tc, Tm )
    slope_y(i) = vMask(i,J)*vMask(i,J-1) * &
        sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
  enddo

end subroutine kernel_PLM_slope_y

subroutine kernel_PLM_Tr_y(G, is, ie, J, scalar, vhh, CFL, slope_y, Tr_y)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, J
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: scalar
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: vhh
  real, dimension(SZI_(G)),          intent(in)    :: CFL
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: slope_y
  real, dimension(SZI_(G)),          intent(inout) :: Tr_y
  ! Local
  integer :: i
  real :: Tc

  do i=is,ie
    if (vhh(i,J) >= 0.0) then
      Tc = scalar(i,j)
      Tr_y(i) = Tc + 0.5 * slope_y(i,j) * ( 1. - CFL(i) )
    else
      Tc = scalar(i,j+1)
      Tr_y(i) = Tc - 0.5 * slope_y(i,j+1) * ( 1. - CFL(i) )
    endif
  enddo

end subroutine kernel_PLM_Tr_y

subroutine kernel_PPMH3_Tr_y(G, is, ie, J, scalar, vMask, vhh, CFL, Tr_y)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, J
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: scalar
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: vMask
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: vhh
  real, dimension(SZI_(G)),          intent(in)    :: CFL
  real, dimension(SZI_(G)),          intent(inout) :: Tr_y
  ! Local variables, all with the same units as scalar.
  real :: Tp, Tc, Tm, aL, aR, dA, a6, mA
  integer :: i

  do i=is,ie
    if (vhh(i,J) >= 0.0) then
      ! Implementation of PPM-H3
      Tp = scalar(i,j+1) ; Tc = scalar(i,j) ; Tm = scalar(i,j-1)
      aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
      aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
      aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
      aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
      dA = aR - aL ; mA = 0.5*( aR + aL )

      ! These expressions are uglier than they might be, but they are less
      ! sensitive to underflow than the alternatives would be.
      if ((vMask(i,J)*vMask(i,J-1) == 0.0) .or. (Tp == Tc) .or. (Tc == Tm) .or. &
          (sign(1.,Tp-Tc)*sign(1.,Tc-Tm) <= 0.)) then
        aL = Tc ; aR = Tc ! PCM for local extrema and boundary cells
        a6 = 0.0 ! Curvature
      elseif ( 6.*sign(1.,dA)*(Tc-mA) > abs(dA) ) then
        ! aL = 3.*Tc - 2.*aR
        aL = Tc + 2.*(Tc - aR)
        a6 = 3.*(aR - Tc) ! Curvature
      elseif ( -6.*sign(1.,dA)*(Tc-mA) > abs(dA) ) then
        ! aR = 3.*Tc - 2.*aL
        aR = Tc + 2.*(Tc - aL)
        a6 = 3.*(aL - Tc) ! Curvature
      else
        a6 = 3.*((Tc - aR) + (Tc - aL)) ! Curvature
      endif
      ! a6 = 6.*Tc - 3. * (aR + aL) ! Curvature

      Tr_y(i) = ( aR - 0.5 * CFL(i) * ( &
            ( aR - aL ) - a6 * ( 1. - 2./3. * CFL(i) ) ) )
    else
      ! Implementation of PPM-H3
      Tp = scalar(i,j+2) ; Tc = scalar(i,j+1) ; Tm = scalar(i,j)
      aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
      aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
      aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
      aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
      dA = aR - aL ; mA = 0.5*( aR + aL )
      ! These expressions are uglier than they might be, but they are less
      ! sensitive to underflow than the alternatives would be.
      if ((vMask(i,J)*vMask(i,J+1) == 0.0) .or. (Tp == Tc) .or. (Tc == Tm) .or. &
          (sign(1.,Tp-Tc)*sign(1.,Tc-Tm) <= 0.)) then
        aL = Tc ; aR = Tc ! PCM for local extrema and boundary cells
        a6 = 0.0 ! Curvature
      elseif ( 6.*sign(1.,dA)*(Tc-mA) > abs(dA) ) then
        ! aL = 3.*Tc - 2.*aR
        aL = Tc + 2.*(Tc - aR)
        a6 = 3.*(aR - Tc) ! Curvature
      elseif ( -6.*sign(1.,dA)*(Tc-mA) > abs(dA) ) then
        ! aR = 3.*Tc - 2.*aL
        aR = Tc + 2.*(Tc - aL)
        a6 = 3.*(aL - Tc) ! Curvature
      else
        a6 = 3.*((Tc - aR) + (Tc - aL)) ! Curvature
      endif
      ! a6 = 6.*Tc - 3. * (aR + aL) ! Curvature
      Tr_y(i) = ( aL + 0.5 * CFL(i) * ( &
            ( aR - aL ) + a6 * ( 1. - 2./3. * CFL(i) ) ) )
    endif
  enddo

end subroutine kernel_PPMH3_Tr_y

subroutine kernel_hlst_y(G, is, ie, j, vh_neglect, vhh, hprev, hlst, Ihnew, do_i, h_neglect)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: vh_neglect, vhh
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: hprev
  real, dimension(SZI_(G)),          intent(inout) :: hlst, Ihnew
  logical, dimension(SZI_(G)),       intent(inout) :: do_i
  real,                              intent(in)    :: h_neglect
  ! Local
  integer :: i

  do i=is,ie
    if ((vhh(i,J) /= 0.0) .or. (vhh(i,J-1) /= 0.0)) then
      do_i(i) = .true.
      hlst(i) = hprev(i,j)
      hprev(i,j) = max(hprev(i,j) - (vhh(i,J) - vhh(i,J-1)), 0.0)
      if (hprev(i,j) <= 0.0) then
        do_i(i) = .false.
      elseif (hprev(i,j) < h_neglect*G%areaT(i,j)) then
        hlst(i) = hlst(i) + (h_neglect*G%areaT(i,j) - hprev(i,j))
        Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
      else
        Ihnew(i) = 1.0 / hprev(i,j)
      endif
    else
      do_i(i) = .false.
    endif
  enddo

end subroutine kernel_hlst_y

!> Updates a scalar with the divergence of y-flux
subroutine kernel_tracer_div_y(G, is, ie, j, do_i, hlst, Ihnew, flux_y,  scalar)
  type(SIS_hor_grid_type),           intent(in)    :: G
  integer,                           intent(in)    :: is, ie, j
  logical, dimension(SZI_(G)),       intent(in)    :: do_i
  real, dimension(SZI_(G)),          intent(in)    :: hlst, Ihnew
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: flux_y
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: scalar
  ! Local
  integer :: i

  do i=is,ie ; if (do_i(i)) then
    scalar(i,j) = (scalar(i,j) * hlst(i) - &
                      (flux_y(i,J) - flux_y(i,J-1))) * Ihnew(i)
  endif ; enddo

end subroutine kernel_tracer_div_y

subroutine advect_upwind_2d(Tr, h_prev, h_end, uhtr, vhtr, ntr, dt, G, IG)
  type(SIS_hor_grid_type),                     intent(inout) :: G
  type(ice_grid_type),                         intent(inout) :: IG
  type(SIS_tracer_type), dimension(ntr),       intent(inout) :: Tr
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),  intent(in) :: h_prev, h_end
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), intent(in) :: uhtr
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), intent(in) :: vhtr
  real,                                   intent(in)    :: dt
  integer,                                intent(in)    :: ntr
! Arguments: tr - The arrays of tracer concentration being worked on.
!  (in)      h_prev - Category thickness times fractional coverage before advection, in m or kg m-2.
!  (in)      h_end - Layer thickness times fractional coverage after advection, in m or kg m-2.
!  (in)      uhtr - Accumulated volume or mass fluxes through zonal faces,
!                   in m3 s-1 or kg s-1.
!  (in)      vhtr - Accumulated volume or mass fluxes through meridional faces,
!                   in m3 s-1 or kg s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 tracer_advect_init.
!  (in)      TrReg - A pointer to the tracer registry.

  real, dimension(SZIB_(G),SZJ_(G)) :: flux_x  ! x-direction tracer fluxes, in conc * kg
  real, dimension(SZI_(G),SZJB_(G)) :: flux_y  ! y-direction tracer fluxes, in conc * kg
  real    :: tr_up              ! Upwind tracer concentrations, in conc.
  real    :: Idt
  real    :: vol_end, Ivol_end  ! Cell volume at the end of a step and its inverse.
  integer :: i, j, k, l, m, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  Idt = 1.0/dt

  ! Reconstruct the old value of h ???
  ! if (h_prev(i,j,k) > 0.0) then
  ! h_last(i,j,k) = h_end(i,j,k) + dt * G%IareaT(i,j) * &
  !        ((uh(I,j,k) - uh(I-1,j,k)) + (vh(i,J,k) - vh(i,J-1,k)))

  ! For now this is just non-directionally split upwind advection.
  do m=1,ntr ; do l=1,Tr(m)%nL ; do k=1,IG%CatIce
    do j=js,je ; do I=is-1,ie
      if (uhtr(I,j,k) >= 0.0) then ; tr_up = Tr(m)%t(i,j,k,l)
      else ; tr_up = Tr(m)%t(i+1,j,k,l) ; endif
      flux_x(I,j) = (dt*uhtr(I,j,k)) * tr_up
    enddo ; enddo

    do J=js-1,je ; do i=is,ie
      if (vhtr(i,J,k) >= 0.0) then ; tr_up = Tr(m)%t(i,j,k,l)
      else ; tr_up = Tr(m)%t(i,j+1,k,l) ; endif
      flux_y(i,J) = (dt*vhtr(i,J,k)) * tr_up
    enddo ; enddo

    do j=js,je ; do i=is,ie
      vol_end = (G%areaT(i,j) * h_end(i,j,k))
      Ivol_end = 0.0 ; if (vol_end > 0.0) Ivol_end = 1.0 / vol_end
      Tr(m)%t(i,j,k,l) = ( (G%areaT(i,j)*h_prev(i,j,k))*Tr(m)%t(i,j,k,l) - &
                       ((flux_x(I,j) - flux_x(I-1,j)) + &
                        (flux_y(i,J) - flux_y(i,J-1))) ) * Ivol_end
    enddo ; enddo

    if (associated(Tr(m)%ad4d_x)) then ; do j=js,je ; do I=is-1,ie
      Tr(m)%ad4d_x(I,j,k,l) = Tr(m)%ad4d_x(I,j,k,l) + flux_x(I,j)*Idt
    enddo ; enddo ; endif
    if (associated(Tr(m)%ad3d_x)) then ; do j=js,je ; do I=is-1,ie
      Tr(m)%ad3d_x(I,j,k) = Tr(m)%ad3d_x(I,j,k) + flux_x(I,j)*Idt
    enddo ; enddo ; endif
    if (associated(Tr(m)%ad2d_x)) then ; do j=js,je ; do I=is-1,ie
      Tr(m)%ad2d_x(I,j) = Tr(m)%ad2d_x(I,j) + flux_x(I,j)*Idt
    enddo ; enddo ; endif

    if (associated(Tr(m)%ad4d_y)) then ; do J=js-1,je ; do i=is,ie
      Tr(m)%ad4d_y(i,J,k,l) = Tr(m)%ad4d_y(i,J,k,l) + flux_y(i,J)*Idt
    enddo ; enddo ; endif
    if (associated(Tr(m)%ad3d_y)) then ; do J=js-1,je ; do i=is,ie
      Tr(m)%ad3d_y(i,J,k) = Tr(m)%ad3d_y(i,J,k) + flux_y(i,J)*Idt
    enddo ; enddo ; endif
    if (associated(Tr(m)%ad2d_y)) then ; do J=js-1,je ; do i=is,ie
      Tr(m)%ad2d_y(i,J) = Tr(m)%ad2d_y(i,J) + flux_y(i,J)*Idt
    enddo ; enddo ; endif
  enddo ; enddo ; enddo

end subroutine advect_upwind_2d

subroutine advect_tracers_thicker(vol_start, vol_trans, G, IG, CS, &
                                  TrReg, snow_tr, j, is, ie)
  type(SIS_hor_grid_type),             intent(in) :: G
  type(ice_grid_type),                 intent(inout) :: IG
  real, dimension(SZI_(G),SZCAT_(IG)), intent(in) :: vol_start, vol_trans
  type(SIS_tracer_advect_CS),          pointer    :: CS
  type(SIS_tracer_registry_type),      pointer    :: TrReg
  logical,                             intent(in) :: snow_tr
  integer,                             intent(in) :: j, is, ie

  real, dimension(SZI_(G),SZCAT_(IG)) :: vol
  type(SIS_tracer_type), dimension(:), pointer :: Tr=>NULL()
  real :: Ivol_new
  integer :: i, k, m, n, ncat

  if (.not. associated(CS)) call SIS_error(FATAL, "SIS_tracer_advect: "// &
       "SIS_tracer_advect_init must be called before advect_tracers_thicker.")
  if (.not. associated(TrReg)) call SIS_error(FATAL, "SIS_tracer_advect: "// &
       "register_tracer must be called before advect_tracers_thicker.")
  if (TrReg%ntr==0) return

  ncat = IG%CatIce

  if (snow_tr) then
    Tr => TrReg%Tr_snow
  else
    Tr => TrReg%Tr_ice
  endif

  do k=1,ncat ; do i=is,ie ; vol(i,k) = vol_start(i,k) ; enddo ; enddo
  do K=1,ncat-1 ; do i=is,ie ; if (vol_trans(i,K) > 0.0) then
    Ivol_new = 1.0 / (vol(i,k+1) + vol_trans(i,K))
    ! This is upwind advection across categories.  Improve it later.
    do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
      Tr(n)%t(i,j,k+1,m) = (vol_trans(i,K)*Tr(n)%t(i,j,k,m) + &
                       vol(i,k+1)*Tr(n)%t(i,j,k+1,m)) * Ivol_new
    enddo ; enddo
    vol(i,k+1) = vol(i,k+1) + vol_trans(i,K)
    vol(i,k) = vol(i,k) - vol_trans(i,K)
  endif ; enddo ; enddo

  do K=ncat-1,1,-1 ; do i=is,ie ; if (vol_trans(i,K) < 0.0) then
    Ivol_new = 1.0 / (vol(i,k) - vol_trans(i,K))
    ! This is upwind advection across categories.  Improve it later.
    do n=1,TrReg%ntr ; do m=1,Tr(n)%nL
      Tr(n)%t(i,j,k,m) = (vol(i,k)*Tr(n)%t(i,j,k,m) - &
                         vol_trans(i,K)*Tr(n)%t(i,j,k+1,m)) * Ivol_new
    enddo ; enddo
    vol(i,k+1) = vol(i,k+1) + vol_trans(i,K)
    vol(i,k) = vol(i,k) - vol_trans(i,K)
  endif ; enddo ; enddo

end subroutine advect_tracers_thicker

subroutine SIS_tracer_advect_init(Time, G, param_file, diag, CS, scheme)
  type(time_type), target, intent(in)    :: Time
  type(SIS_hor_grid_type), intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(SIS_tracer_advect_CS),  pointer       :: CS
  character(len=*), optional, intent(in) :: scheme
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer to the control structure for this module
  integer, save :: init_calls = 0
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_tracer_advect" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_tracer_advect_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  if ((first_call) .or. .not.present(scheme)) &
    call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "DT_ICE_DYNAMICS", CS%dt, &
                 "The time step used for the slow ice dynamics, including "//&
                 "stepping the continuity equation and interactions between "//&
                 "the ice mass field and velocities.", units="s", &
                 default=-1.0, do_not_log=.true.)
  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false.)
  if (present(scheme)) then ; mesg = scheme ; else
    call get_param(param_file, mod, "SIS_TRACER_ADVECTION_SCHEME", mesg, &
          desc="The horizontal transport scheme for tracers:\n"//&
          "  UPWIND_2D - Non-directionally split upwind\n"//&
          "  PCM    - Directionally split piecewise constant\n"//&
          "  PLM    - Piecewise Linear Method\n"//&
          "  PPM:H3 - Piecewise Parabolic Method (Huyhn 3rd order)", &
          default='UPWIND_2D')
  endif
  CS%use_upwind2d = .false. ; CS%usePPM = .false. ; CS%usePCM = .false.
  select case (trim(mesg))
    case ("UPWIND_2D")
      CS%use_upwind2d = .true.
    case ("PCM")
      CS%usePCM = .true.
    case ("PLM")
      CS%usePPM = .false.
    case ("PPM:H3")
      CS%usePPM = .true.
    case default
      if (present(scheme)) then
        call SIS_error(FATAL, "SIS_tracer_advect, SIS_tracer_advect_init: "//&
           "Unknown input scheme "//trim(mesg))
      else
        call SIS_error(FATAL, "SIS_tracer_advect, SIS_tracer_advect_init: "//&
           "Unknown SIS_TRACER_ADVECTION_SCHEME = "//trim(mesg))
      endif
  end select

  if (first_call) then
    id_clock_advect = cpu_clock_id('(Ocean advect tracer)', grain=CLOCK_MODULE)
    id_clock_pass = cpu_clock_id('(Ocean tracer halo updates)', grain=CLOCK_ROUTINE)
    id_clock_sync = cpu_clock_id('(Ocean tracer global synch)', grain=CLOCK_ROUTINE)
    first_call = .false.
  endif

end subroutine SIS_tracer_advect_init

subroutine SIS_tracer_advect_end(CS)
  type(SIS_tracer_advect_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine SIS_tracer_advect_end

end module SIS_tracer_advect
