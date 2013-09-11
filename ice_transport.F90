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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   This module does the transport and redistribution between thickness        !
! categories for the SIS2 sea ice model.                                       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_transport_mod

use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser, only : get_param, log_param, read_param, log_version, param_file_type
use MOM_domains,     only : pass_var, pass_vector, BGRID_NE, CGRID_NE

use ice_grid_mod, only : sea_ice_grid_type
use ice_thm_mod,  only : get_thermo_coefs

implicit none ; private

#include <SIS2_memory.h>

public :: ice_transport_init, ice_transport, ice_transport_end !, ice_transport_register_restarts

type, public :: ice_transport_CS ; private

  ! parameters for doing advective and parameterized advection.
  logical :: SLAB_ICE = .false. ! should we do old style GFDL slab ice?
  real    :: chan_visc  = 0.     ! viscosity used in one-cell wide channels to parameterize transport (m^2/s)
  real    :: smag_ocn           = 0.15   ! Smagorinksy coefficient for viscosity (dimensionless)
  real    :: chan_cfl_limit     = 0.25   ! CFL limit for channel viscosity parameterization (dimensionless)
  real :: Rho_ocean = 1030.0  ! The nominal density of sea water, in kg m-3.
  real :: Rho_ice = 905.0     ! The nominal density of sea ice, in kg m-3.
  real :: Rho_snow = 330.0    ! The nominal density of snow on sea ice, in
                              ! kg m-3.
  logical :: specified_ice    ! If true, the sea ice is specified and there is
                              ! no need for ice dynamics.
  integer :: adv_sub_steps    ! The number of advective iterations for each slow
                              ! time step.
  type(time_type), pointer :: Time ! A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_ustar = -1, id_uocean = -1, id_uchan = -1
  integer :: id_vstar = -1, id_vocean = -1, id_vchan = -1
  integer :: id_ix_trans = -1, id_iy_trans = -1
end type ice_transport_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! transport - do ice transport and thickness class redistribution              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_transport (part_sz, h_ice, h_snow, uc, vc, t_ice, t_snow, &
                          sea_lev, hlim, dt_slow, G, CS)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT0_(G)), intent(inout) :: part_sz
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: h_ice, h_snow
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G),SZK_ICE_(G)), intent(inout) :: t_ice
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: t_snow
  real, dimension(SZIB_(G),SZJ_(G)),           intent(inout) :: uc
  real, dimension(SZI_(G),SZJB_(G)),           intent(inout) :: vc
  real, dimension(SZI_(G),SZJ_(G)),            intent(in)    :: sea_lev
  real, dimension(:),      intent(in) :: hlim  ! Move to grid type?
  real,                    intent(in) :: dt_slow
  type(ice_transport_CS), pointer :: CS
! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1, in/out.
!  (inout)   h_ice - The thickness of the ice in each category in m.
!  (inout)   h_snow - The thickness of the snow atop the ice in each category
!                     in m.
!  (in)      uc - The zonal ice velocity, in m s-1.
!  (in)      vc - The meridional ice velocity, in m s-1.
!  (inout)   t_ice - The temperature of the ice in each category and layer
!                    within the ice in degC.
!  (inout)   t_snow - The temperature of the snow atop the ice in each category
!                     in degC.
!  (in)      hlim - The lower thickness limit of each category, in m.
!  (in)      sea_lev - The height of the sea level, including contributions
!                      from non-levitating ice from an earlier time step, in m.
!  (in)      dt_slow - The amount of time over which the ice dynamics are to be
!                      advanced, in s.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.
  
  real, dimension(G%NkIce) :: pocket_enth ! Coefficients relating to the enthalpy
                                          ! contribution from melting in brine
                                          ! pockets, in degC2.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    uf0, uf, & ! Zonal fluxes in m3 s-1 and kg s-1.
    ustar, ustaro, ustarv ! Local variables, transporting velocities
  real, dimension(SZI_(G),SZJB_(G)) :: &
    vf0, vf, & ! Zonal fluxes in m3 s-1 and kg s-1.
    vstar, vstaro, vstarv ! Local variables, transporting velocities
  real, dimension(SZI_(G),SZJ_(G)) :: tmp1 ! Local variables, 2D ice concentration
  real :: u_visc, u_ocn, cnn, grad_eta ! Variables for channel parameterization

  real :: dt_adv
  integer :: i, j, k, l, bad, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (CS%slab_ice) then
    call pass_vector(uc, vc, G%Domain, stagger=CGRID_NE)
    call slab_ice_advect(uc, vc, h_ice(:,:,2), 4.0, dt_slow, G, CS)
    call pass_var(h_ice(:,:,2), G%Domain)
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      if (h_ice(i,j,2) > 0.0) then
        part_sz(i,j,2) = 1.0
      else 
        part_sz(i,j,2) = 0.0
      endif
    enddo ; enddo
    return
  endif

  call get_thermo_coefs(layer_coefs=pocket_enth)
  !   Convert from ice temperature (which is not conserved) to enthalpy, which
  ! includes the heat requirements for melting of brine pockets associated with 
  ! temperature changes.  These expressions stem from the assumptions that
  ! brine pockets will shrink or grow until their salinity gives a freezing
  ! point that matches the local temperature.
  ! This was prevously the subroutine thm_pack.
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    if (h_ice(i,j,k)>0.0) then
      h_ice(i,j,k) = part_sz(i,j,k)*h_ice(i,j,k)
      do l=1,G%NkIce
        t_ice(i,j,k,l) = (t_ice(i,j,k,l) - pocket_enth(l)/t_ice(i,j,k,l)) * &
                          h_ice(i,j,k)
      enddo
      h_snow(i,j,k) = part_sz(i,j,k)*h_snow(i,j,k)
      t_snow(i,j,k) = t_snow(i,j,k)*h_snow(i,j,k)
    else 
      part_sz(i,j,k) = 0.0 ; h_ice(i,j,k) = 0.0
      do l=1,G%NkIce ; t_ice(i,j,k,l) = 0.0 ; enddo
      h_snow(i,j,k) = 0.0
      t_snow(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo
  
  call pass_var(part_sz, G%Domain) ! cannot be combined with updates below
  do l=1,G%NkIce ! The do loop allows these to be combined with the following.
    call pass_var(t_ice(:,:,:,l), G%Domain, complete=.false.)
  enddo
  call pass_var(t_snow, G%Domain, complete=.false.)
  call pass_var(h_ice,  G%Domain, complete=.false.)
  call pass_var(h_snow, G%Domain, complete=.true.)

  if (CS%chan_visc>0. .and. CS%adv_sub_steps>0) then
  ! This block of code is a parameterization of either (or both)
  ! i) the pressure driven oceanic flow in a narrow channel, or
  ! ii) pressure driven flow of ice itself.
  ! The latter is a speculative but both are missing due to
  ! masking of velocities to zero in a single-cell wide channel.
    dt_adv = dt_slow/CS%adv_sub_steps

    tmp1=1.-max(1.-sum(part_sz(:,:,1:G%CatIce),dim=3),0.0)
    ustar(:,:)=0.; vstar(:,:)=0.
    ustaro(:,:)=0.; vstaro(:,:)=0.
    ustarv(:,:)=0.; vstarv(:,:)=0.
    do j=jsc,jec ; do I=isc-1,iec
      if ((uc(I,j)==0.) .and. & ! this is a redundant test due to following line
          (G%mask2dBu(I,J)+G%mask2dBu(I,J-1)==0.) .and. &  ! =0 => no vels
          (G%mask2dT(i,j)*G%mask2dT(i+1,j)>0.)) then ! >0 => open for transport
        grad_eta=(sea_lev(i+1,j)-sea_lev(i,j)) * G%IdxCu(I,j) ! delta_i eta/dx
        u_visc=-G%g_Earth*((G%dyCu(I,j)*G%dyCu(I,j))/(12.*CS%chan_visc)) & ! -g*dy^2/(12*visc)
                 *grad_eta                                                  ! d/dx eta
        u_ocn=sqrt( G%g_Earth*G%dyCu(I,j)*abs(grad_eta)/(36.*CS%smag_ocn) ) ! Magnitude of ocean current
        u_ocn=sign(u_ocn, -grad_eta) ! Direct down the ssh gradient
        cnn=max(tmp1(i,j),tmp1(i+1,j))**2. ! Use the larger concentration
        uc(I,j)=cnn*u_visc+(1.-cnn)*u_ocn
        ! Limit flow to be stable for fully divergent flow
        if (uc(I,j)>0.) then
          uc(I,j)=min( uc(I,j), CS%chan_cfl_limit*G%dxT(i,j)/dt_adv)
        else
          uc(I,j)=max( uc(I,j),(-1*CS%chan_cfl_limit)*G%dxT(i+1,j)/dt_adv)
        endif
        if (CS%id_ustar>0) ustar(I,j)=uc(I,j)
        if (CS%id_uocean>0) ustaro(I,j)=u_ocn
        if (CS%id_uchan>0) ustarv(I,j)=u_visc
      endif
    enddo ; enddo
    do J=jsc-1,jec ; do i=isc,iec
      if ((vc(i,J)==0.) .and. & ! this is a redundant test due to following line
          (G%mask2dBu(I,J)+G%mask2dBu(I-1,J)==0.) .and. &  ! =0 => no vels
          (G%mask2dT(i,j)*G%mask2dT(i,j+1)>0.)) then ! >0 => open for transport
        grad_eta=(sea_lev(i,j+1)-sea_lev(i,j)) * G%IdyCv(i,J) ! delta_i eta/dy
        u_visc=-G%g_Earth*((G%dxCv(i,J)*G%dxCv(i,J))/(12.*CS%chan_visc)) & ! -g*dy^2/(12*visc)
                *grad_eta                                                  ! d/dx eta
        u_ocn=sqrt( G%g_Earth*G%dxCv(i,J)*abs(grad_eta)/(36.*CS%smag_ocn) ) ! Magnitude of ocean current
        u_ocn=sign(u_ocn, -grad_eta) ! Direct down the ssh gradient
        cnn=max(tmp1(i,j),tmp1(i,j+1))**2. ! Use the larger concentration
        vc(i,J)=cnn*u_visc+(1.-cnn)*u_ocn
        ! Limit flow to be stable for fully divergent flow
        if (vc(i,J)>0.) then
         vc(i,J)=min( vc(i,J), CS%chan_cfl_limit*G%dyT(i,j)/dt_adv)
        else
         vc(i,J)=max( vc(i,J),(-1*CS%chan_cfl_limit)*G%dyT(i,j+1)/dt_adv)
        endif
        if (CS%id_vstar>0) vstar(i,J)=vc(i,J)
        if (CS%id_vocean>0) vstaro(i,J)=u_ocn
        if (CS%id_vchan>0) vstarv(i,J)=u_visc
      endif
    enddo ; enddo
  endif

  call pass_vector(uc, vc, G%Domain, stagger=CGRID_NE)

  uf(:,:) = 0.0; vf(:,:) = 0.0
  do k=1,G%CatIce
    call ice_advect(uc, vc, part_sz(:,:,k), dt_slow, G, CS)
    call ice_advect(uc, vc, h_snow(:,:,k), dt_slow, G, CS, uf0, vf0)
    uf = uf + CS%Rho_snow*uf0; vf = vf + CS%Rho_snow*vf0
    call ice_advect(uc, vc, h_ice(:,:,k), dt_slow, G, CS, uf0, vf0)
    uf = uf + CS%Rho_ice*uf0; vf = vf + CS%Rho_ice*vf0
    call ice_advect(uc, vc, t_snow(:,:,k), dt_slow, G, CS)
    do l=1,G%NkIce
      call ice_advect(uc, vc, t_ice(:,:,k,l), dt_slow, G, CS)
    enddo
  enddo
  call post_SIS_data(CS%id_ix_trans, uf, CS%diag)
  call post_SIS_data(CS%id_iy_trans, vf, CS%diag)

  do j=jsc, jec
     do i=isc, iec
        if (sum(h_ice(i,j,:))>0) &
          call ice_redistribute(part_sz(i,j,1:G%CatIce), G, &
             h_snow(i,j,:), t_snow(i,j,:), h_ice (i,j,:), &
             t_ice(i,j,:,1), t_ice(i,j,:,2), &
             t_ice(i,j,:,3), t_ice(i,j,:,4), hlim)
     end do
  end do

  !   Convert from enthalpy back to ice temperature. These expressions stem from
  ! the assumptions that brine pockets will shrink or grow until their salinity
  ! gives a freezing point that matches the local temperature.
  ! This was prevously the subroutine thm_unpack.
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    if (part_sz(i,j,k)<1e-10) h_ice(i,j,k) = 0.0
    if (h_ice(i,j,k)>0.0) then
      do l=1,G%NkIce
        t_ice(i,j,k,l) = t_ice(i,j,k,l) / h_ice(i,j,k)
        t_ice(i,j,k,l) = 0.5*(t_ice(i,j,k,l) - &
                              sqrt(t_ice(i,j,k,l)*t_ice(i,j,k,l) + 4.0*pocket_enth(l)))
      enddo
      h_ice(i,j,k) = h_ice(i,j,k)/part_sz(i,j,k)
      if (h_snow(i,j,k)>0.0) then
        t_snow(i,j,k) = t_snow(i,j,k)/h_snow(i,j,k)
      else
        t_snow(i,j,k) = 0.0
      endif
      h_snow(i,j,k) = h_snow(i,j,k)/part_sz(i,j,k)
      
      ! Check for bad values.
      bad = 0
      if (h_snow(i,j,k) <0.0 .or. h_snow(i,j,k) > 1e3 ) bad = bad+1
      if (h_ice(i,j,k) <0.0 .or. h_ice(i,j,k) > 1e3 ) bad = bad+1
      if (t_snow(i,j,k)>0.0.or.t_snow(i,j,k)<-100.0) bad = bad+1
      do l=1,G%NkIce
        if (t_ice(i,j,k,l) > 0.0 .or. t_ice(i,j,k,l) < -100.0) bad = bad+1
      enddo      
!### USE BETTER ERROR HANDLING LATER.
      if (bad>0) then
        print *, 'BAD ICE AFTER UNPACK ', 'hs/hi=',h_snow(i,j,k),h_ice(i,j,k),'tsn/tice=', &
                          t_snow(i,j,k),t_ice(i,j,k,1),t_ice(i,j,k,2),t_ice(i,j,k,3), &
                          t_ice(i,j,k,4),'cn=',part_sz(i,j,k)
      endif
    else
      part_sz(i,j,k) = 0.0 ; h_ice(i,j,k) = 0.0
      do l=1,G%NkIce ; t_ice(i,j,k,l) = 0.0 ; enddo
      h_snow(i,j,k) = 0.0
      t_snow(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo

  call pass_var(part_sz, G%Domain) ! cannot be combined with the two updates below
  call pass_var(h_snow, G%Domain, complete=.false.)
  call pass_var(h_ice, G%Domain, complete=.true.)

  if (CS%id_ustar >0) call post_SIS_data(CS%id_ustar, ustar, CS%diag)
  if (CS%id_vstar >0) call post_SIS_data(CS%id_vstar , vstar, CS%diag)
  if (CS%id_vocean>0) call post_SIS_data(CS%id_vocean, vstaro, CS%diag)
  if (CS%id_uocean>0) call post_SIS_data(CS%id_uocean, ustaro, CS%diag)
  if (CS%id_vchan>0)  call post_SIS_data(CS%id_vchan,  vstarv, CS%diag)
  if (CS%id_uchan>0)  call post_SIS_data(CS%id_uchan,  ustarv, CS%diag)

end subroutine ice_transport


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_advect - take adv_sub_steps upstream advection timesteps                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_advect(uc, vc, trc, dt_slow, G, CS, uf, vf)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: uc  ! x-face advecting velocity
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: vc  ! y-face advecting velocity
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: trc ! tracer to advect
  real,                              intent(in   ) :: dt_slow
  type(ice_transport_CS),            pointer       :: CS
  real, dimension(SZIB_(G),SZJ_(G)), optional, intent(  out) :: uf
  real, dimension(SZI_(G),SZJB_(G)), optional, intent(  out) :: vf
! Arguments: uc - The zonal ice velocity, in m s-1.
!  (in)      vc - The meridional ice velocity, in m s-1.
!  (inout)   trc - A tracer concentration times thickness, in m kg kg-1 or
!                  other units.
!  (in)      dt_slow - The amount of time over which the ice dynamics are to be
!                      advanced, in s.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.
!  (out)     vf - The averaged zonal tracer flux, in m3 kg kg-1 s-1.
!  (out)     vf - The averaged meridional tracer flux, in m3 kg kg-1 s-1.

  real, dimension(SZIB_(G),SZJ_(G)) :: uflx
  real, dimension(SZI_(G),SZJB_(G)) :: vflx
  real :: dt_adv, I_adv_steps
  integer :: l, i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (CS%adv_sub_steps==0) return;
  dt_adv = dt_slow/CS%adv_sub_steps

  if (present(uf)) uf(:,:) = 0.0
  if (present(vf)) vf(:,:) = 0.0

  uflx(:,:) = 0.0
  vflx(:,:) = 0.0

  do l=1,CS%adv_sub_steps
    do j=jsc,jec ; do I=isc-1,iec
      if( uc(I,j) > 0.0 ) then
         uflx(I,j) = uc(I,j) * trc(i,j) * G%dyCu(I,j)
      else
         uflx(I,j) = uc(I,j) * trc(i+1,j) * G%dyCu(I,j)
      endif
    enddo ; enddo

    do J=jsc-1,jec ; do i=isc,iec
      if( vc(i,J) > 0.0 ) then
         vflx(i,J) = vc(i,J) * trc(i,j) * G%dxCv(i,J)
      else
         vflx(i,J) = vc(i,J) * trc(i,j+1) * G%dxCv(i,J)
      endif
    enddo ; enddo

    do j=jsc,jec ; do i=isc,iec
      trc(i,j) = trc(i,j) + dt_adv * ( (uflx(I-1,j) - uflx(I,j)) + &
                 (vflx(i,J-1) - vflx(i,J)) ) * G%IareaT(i,j)
    enddo ; enddo

    call pass_var(trc, G%Domain)

    if (present(uf)) then ; do j=jsc,jec ; do I=isc-1,iec       
      uf(I,j) = uf(I,j) + uflx(I,j)
    enddo ; enddo ; endif

    if (present(vf)) then ;  do J=jsc-1,jec ; do i=isc,iec       
      vf(i,J) = vf(i,J) + vflx(i,J)
    enddo ; enddo ; endif
  enddo

  if (CS%adv_sub_steps>1) then
    I_adv_steps = 1.0/CS%adv_sub_steps
    if (present(uf)) then ; do j=jsc,jec ; do I=isc-1,iec       
      uf(I,j) = uf(I,j) * I_adv_steps
    enddo ; enddo ; endif
    if (present(vf)) then ;  do J=jsc-1,jec ; do i=isc,iec       
      vf(i,J) = vf(i,J) * I_adv_steps
    enddo ; enddo ; endif
  endif

end subroutine ice_advect

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine slab_ice_advect(uc, vc, trc, stop_lim, dt_slow, G, CS)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: uc  ! x-face advecting velocity
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: vc  ! y-face advecting velocity
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: trc ! tracer to advect
  real,                              intent(in   ) :: stop_lim
  real,                              intent(in   ) :: dt_slow
  type(ice_transport_CS),            pointer       :: CS
! Arguments: uc - The zonal ice velocity, in m s-1.
!  (in)      vc - The meridional ice velocity, in m s-1.
!  (inout)   trc - A tracer concentration times thickness, in m kg kg-1 or
!                  other units.
!  (in)      stop_lim - ?
!  (in)      dt_slow - The amount of time over which the ice dynamics are to be
!                      advanced, in s.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.

  real, dimension(SZIB_(G),SZJ_(G)) :: uflx
  real, dimension(SZI_(G),SZJB_(G)) :: vflx
  real :: avg, dif
  real :: dt_adv
  integer :: l, i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec


  if (CS%adv_sub_steps==0) return;
  dt_adv = dt_slow/CS%adv_sub_steps


  do l=1,CS%adv_sub_steps
    do j=jsc,jec ; do I=isc-1,iec
      avg = ( trc(i,j) + trc(i+1,j) )/2
      dif = trc(i+1,j) - trc(i,j)
      if( avg > stop_lim .and. uc(I,j) * dif > 0.0) then
        uflx(I,j) = 0.0
      else if( uc(i,j) > 0.0 ) then
        uflx(I,j) = uc(I,j) * trc(i,j) * G%dyCu(I,j)
      else
        uflx(I,j) = uc(I,j) * trc(i+1,j) * G%dyCu(I,j)
      endif
    enddo ; enddo

    do J=jsc-1,jec ; do i=isc,iec
      avg = ( trc(i,j) + trc(i,j+1) )/2
      dif = trc(i,j+1) - trc(i,j)
      if( avg > stop_lim .and. vc(i,J) * dif > 0.0) then
        vflx(i,J) = 0.0
      else if( vc(i,J) > 0.0 ) then
        vflx(i,J) = vc(i,J) * trc(i,j) * G%dxCv(i,J)
      else
        vflx(i,J) = vc(i,J) * trc(i,j+1) * G%dxCv(i,J)
      endif
    enddo ; enddo

    do j=jsc,jec ; do i=isc,iec
      trc(i,j) = trc(i,j) + dt_adv * ((uflx(I-1,j) - uflx(I,j)) + &
                                      (vflx(i,J-1) - vflx(i,J)) ) * G%IareaT(i,j)
    enddo ; enddo

    call pass_var(trc, G%Domain)
  enddo

end subroutine slab_ice_advect

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_redistribute - a simple ice redistribution scheme from Igor Polyakov     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_redistribute(cn, G, hs, tsn, hi, t1, t2, t3, t4, hlim)
!!$    real, intent(inout), dimension(1:G%CatIce)           :: cn, hs, hi, t1, t2
  type(sea_ice_grid_type), intent(inout) :: G
  real, intent(inout), dimension(:) :: cn, hs, tsn, hi, t1, t2, t3, t4
  real, dimension(:), intent(in) :: hlim

  real    :: cw                                 ! open water concentration
  integer :: k

  cw = 1-sum(cn)
  if (cw<0) cn(1) = cw + cn(1) ! open water has been eliminated by convergence

  do k=1,G%CatIce-1 ; if (cn(k)<0) then
    cn(k+1) = cn(k+1)+cn(k) ; cn(k) = 0 ! pass concentration deficit up to
    hs(k+1) = hs(k+1)+hs(k) ; hs(k) = 0 ! next thicker category
    tsn(k+1) = tsn(k+1)+tsn(k) ; tsn(k) = 0
    hi(k+1) = hi(k+1)+hi(k) ; hi(k) = 0
    t1(k+1) = t1(k+1)+t1(k) ; t1(k) = 0 ! NOTE:  here between the thm_pack and
    t2(k+1) = t2(k+1)+t2(k) ; t2(k) = 0 ! thm_unpack calls, all quantities are
    t3(k+1) = t3(k+1)+t3(k) ; t3(k) = 0 ! extensive, so we add instead of
    t4(k+1) = t4(k+1)+t4(k) ; t4(k) = 0 ! averaging
  endif ; enddo

  do k=1,G%CatIce-1 ; if (hi(k)>hlim(k+1)*cn(k)) then
    cn(k+1) = cn(k+1)+cn(k) ; cn(k) = 0 ! upper thickness limit exceeded
    hs(k+1) = hs(k+1)+hs(k) ; hs(k) = 0 ! move ice up to next thicker category
    tsn(k+1) = tsn(k+1)+tsn(k) ; tsn(k) = 0
    hi(k+1) = hi(k+1)+hi(k) ; hi(k) = 0
    t1(k+1) = t1(k+1)+t1(k) ; t1(k) = 0
    t2(k+1) = t2(k+1)+t2(k) ; t2(k) = 0
    t3(k+1) = t3(k+1)+t3(k) ; t3(k) = 0
    t4(k+1) = t4(k+1)+t4(k) ; t4(k) = 0
  endif ; enddo

  do k=G%CatIce,2,-1 ; if (hi(k)<hlim(k)*cn(k)) then
    cn(k-1) = cn(k-1)+cn(k) ; cn(k) = 0  ! lower thickness limit exceeded;
    hs(k-1) = hs(k-1)+hs(k) ; hs(k) = 0  ! move ice down to thinner category
    tsn(k-1) = tsn(k-1)+tsn(k) ; tsn(k) = 0
    hi(k-1) = hi(k-1)+hi(k) ; hi(k) = 0
    t1(k-1) = t1(k-1)+t1(k) ; t1(k) = 0
    t2(k-1) = t2(k-1)+t2(k) ; t2(k) = 0
    t3(k-1) = t3(k-1)+t3(k) ; t3(k) = 0
    t4(k-1) = t4(k-1)+t4(k) ; t4(k) = 0
  endif ; enddo

end subroutine ice_redistribute

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_transport_init - initialize the ice transport and set parameters.        !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_transport_init(Time, G, param_file, diag, CS)
  type(time_type),     target, intent(in)    :: Time
  type(sea_ice_grid_type),     intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(ice_transport_CS),      pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.

!   This subroutine sets the parameters and registers the diagnostics associated
! with the ice dynamics.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "ice_transport" ! This module's name.
  real, parameter :: missing = -1e34

  if (associated(CS)) then
    call SIS_error(WARNING, "ice_transport_init called with an associated control structure.")
!    call SIS_error(FATAL, "ice_transport_init called with an unassociated control structure. \n"//&
!                    "ice_transport_register_restarts must be called before ice_transport_init.")
    return
  endif
  allocate(CS)

  CS%diag => diag
  CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "SPECIFIED_ICE", CS%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  if ( CS%specified_ice ) then
    CS%adv_sub_steps = 0
    call log_param(param_file, mod, "NSTEPS_ADV", CS%adv_sub_steps, &
                 "The number of advective iterations for each slow time \n"//&
                 "step.  With SPECIFIED_ICE this is always 0.")
  else
    call get_param(param_file, mod, "NSTEPS_ADV", CS%adv_sub_steps, &
                 "The number of advective iterations for each slow time \n"//&
                 "step.", default=1)
  endif

  call get_param(param_file, mod, "ICE_CHANNEL_VISCOSITY", CS%chan_visc, &
                 "A viscosity used in one-cell wide channels to \n"//&
                 "parameterize transport, especially with B-grid sea ice \n"//&
                 "coupled to a C-grid ocean model.", units="m2 s-1", default=0.0)
  call get_param(param_file, mod, "ICE_CHANNEL_SMAG_COEF", CS%smag_ocn, &
                 "A Smagorinsky coefficient for viscosity in channels.", &
                 units="Nondim", default=0.15)
  call get_param(param_file, mod, "ICE_CHANNEL_CFL_LIMIT", CS%chan_cfl_limit, &
                 "The CFL limit that is applied to the parameterized \n"//&
                 "viscous transport in single-point channels.", &
                 units="Nondim", default=0.25)

  call get_param(param_file, mod, "RHO_OCEAN", CS%Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0)
  call get_param(param_file, mod, "RHO_ICE", CS%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0)
  call get_param(param_file, mod, "RHO_SNOW", CS%Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0)

  call get_param(param_file, mod, "USE_SLAB_ICE", CS%SLAB_ICE, &
                 "If true, use the very old slab-style ice.", default=.false.)

  CS%id_ustar = register_diag_field('ice_model', 'U_STAR', diag%axesCu1, Time, &
              'channel transport velocity - x component', 'm/s', missing_value=missing)
  CS%id_vstar = register_diag_field('ice_model', 'V_STAR', diag%axesCv1, Time, &
              'channel transport velocity - y component', 'm/s', missing_value=missing)
  CS%id_uocean = register_diag_field('ice_model', 'U_CHAN_OCN', diag%axesCu1, Time, &
              'ocean component of channel transport - x', 'm/s', missing_value=missing)
  CS%id_vocean = register_diag_field('ice_model', 'V_CHAN_OCN', diag%axesCv1, Time, &
              'ocean component of channel transport - y', 'm/s', missing_value=missing)
  CS%id_uchan = register_diag_field('ice_model', 'U_CHAN_VISC', diag%axesCu1, Time, &
              'viscous component of channel transport - x', 'm/s', missing_value=missing)
  CS%id_vchan = register_diag_field('ice_model', 'V_CHAN_VISC', diag%axesCv1, Time, &
              'viscous component of channel transport - y', 'm/s', missing_value=missing)
  CS%id_ix_trans = register_diag_field('ice_model', 'IX_TRANS', diag%axesCu1, Time, &
               'x-direction ice transport', 'kg/s', missing_value=missing)
  CS%id_iy_trans = register_diag_field('ice_model', 'IY_TRANS', diag%axesCv1, Time, &
               'y-direction ice transport', 'kg/s', missing_value=missing)

end subroutine ice_transport_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_transport_end - deallocate the memory associated with this module.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_transport_end(CS)
  type(ice_transport_CS), pointer :: CS

  deallocate(CS)
end subroutine ice_transport_end

end module ice_transport_mod
