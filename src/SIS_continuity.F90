!> The sea ice continuity solver
module SIS_continuity

! This file is a part of SIS2.  See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg and Alistair Adcroft, September 2006 - .        *
!*                                                                     *
!*    This program contains the subroutine that advects layer          *
!*  thickness.  The scheme here uses a Piecewise-Parabolic method with *
!*  a positive definite limiter.                                       *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use ice_grid,          only : ice_grid_type
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_obsolete_params, only : obsolete_logical
use MOM_unit_scaling,  only : unit_scale_type
use SIS_diag_mediator, only : time_type, SIS_diag_ctrl
use SIS_hor_grid,      only : SIS_hor_grid_type

implicit none ; private

#include <SIS2_memory.h>

public ice_continuity, SIS_continuity_init, SIS_continuity_end
public summed_continuity, proportionate_continuity, ice_cover_transport

integer :: id_clock_update  !< A CPU time clock ID
integer :: id_clock_correct !< A CPU time clock ID

!> The control structure with parameters regulating the continuity solver
type, public :: SIS_continuity_CS ; private
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.
  logical :: use_upwind2d    !< If true, use the non-split upwind scheme that was
                             !! used in older versions of SIS.
  logical :: upwind_1st      !< If true, use a directionally-split first-order upwind scheme.
  logical :: monotonic       !< If true, use the Colella & Woodward monotonic limiter;
                             !! otherwise use a simple positive definite limiter.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic mean) interpolation
                             !! of the edge values instead of the higher order interpolation.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths to the tracer
                             !! cell areas when estimating CFL numbers.
  real :: h_neglect_cont     !< The category ice mass per ocean cell area below which the
                             !! transport within this thickness category of out of a cell is
                             !! set to zero [R Z ~> kg m-2]
  real :: frac_neglect       !< When the total fluxes are distributed between categories, any
                             !! category whose ice is less than this fraction of the total mass
                             !! contributes no flux [nondim]
end type SIS_continuity_CS

!> This type is used to specify the active loop bounds
type :: loop_bounds_type ; private
  !>@{ The active index range
  integer :: ish, ieh, jsh, jeh
  !!@}
end type loop_bounds_type

contains

!> ice_continuity time steps the category thickness changes due to advection,
!! using a monotonically limited, directionally split PPM scheme.
subroutine ice_continuity(u, v, hin, h, uh, vh, dt, G, US, IG, CS, use_h_neg, masking_uh, masking_vh)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: u   !< Zonal ice velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: v   !< Meridional ice velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(in)    :: hin !< Initial ice or snow thickness by category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(inout) :: h   !< Final ice or snow thickness by category [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(out)   :: uh  !< Volume flux through zonal faces = u*h*dy
                                                !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                           intent(out)   :: vh  !< Volume flux through meridional faces = v*h*dx
                                                !! [R Z L2 T-1 ~> kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_continuity_CS), pointer       :: CS  !< The control structure returned by a
                                                !! previous call to SIS_continuity_init.
  logical,       optional, intent(in)    :: use_h_neg !< If true, only move mass out of a cell if it
                                                !! is thicker than CS%h_neglect_cont.
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: masking_uh  !< If this is 0, uh = 0.  Often this is another
                                                !! zonal mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: masking_vh  !< If this is 0, vh = 0.  Often this is another
                                                !! meridional mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].

!    This subroutine time steps the category thicknesses, using a monotonically
!  limit, directionally split PPM scheme, based on Lin (1994).  In the following
!  documentation, H is used for the units of thickness (usually m or kg m-2.)

  ! Local variables
  type(loop_bounds_type) :: LB  ! A structure with the active loop bounds.
  real    :: h_up ! The upwind thickness [R Z ~> kg m-2]
  logical :: apply_h_neg
  integer :: is, ie, js, je, nCat, stensil
  integer :: i, j, k

  logical :: x_first
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nCat = IG%CatIce

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_continuity: Module must be initialized before it is used.")
  x_first = (MOD(G%first_direction,2) == 0)

  apply_h_neg = .false.
  if (present(use_h_neg)) then
    apply_h_neg = ((use_h_neg) .and. (CS%h_neglect_cont > 0.0))
  endif

  stensil = 3 ; if (CS%simple_2nd) stensil = 2 ; if (CS%upwind_1st) stensil = 1

  do k=1,nCat ; do j=js,je ; do i=is,ie ; if (h(i,j,k) < 0.0) then
    call SIS_error(FATAL, 'Negative thickness input to ice_continuity().')
  endif ; enddo ; enddo ; enddo

  if (CS%use_upwind2d) then
    ! This reproduces the scheme that was originally used in SIS1.
    !$OMP parallel default(shared) private(h_up)
    !$OMP do
    do j=js,je ; do k=1,nCat ; do I=is-1,ie
      if (u(I,j) >= 0.0) then ; h_up = hin(i,j,k)
      else ; h_up = hin(i+1,j,k) ; endif
      if (apply_h_neg .and. (h_up < CS%h_neglect_cont)) h_up = 0.0

      uh(I,j,k) = G%dy_Cu(I,j) * u(I,j) * h_up
      if (present(masking_uh)) then
        if (masking_uh(I,j,k) == 0.0) uh(I,j,k) = 0.0
        if (abs(uh(I,j,k)) < CS%frac_neglect*abs(masking_uh(I,j,k))) uh(I,j,k) = 0.0
      endif
    enddo ; enddo ; enddo
    !$OMP do
    do J=js-1,je ; do k=1,nCat ; do i=is,ie
      if (v(i,J) >= 0.0) then ; h_up = hin(i,j,k)
      else ; h_up = hin(i,j+1,k) ; endif
      if (apply_h_neg .and. (h_up < CS%h_neglect_cont)) h_up = 0.0

      vh(i,J,k) = G%dx_Cv(i,J) * v(i,J) * h_up
      if (present(masking_vh)) then
        if (masking_vh(i,J,k) == 0.0) vh(i,J,k) = 0.0
        if (abs(vh(i,J,k)) < CS%frac_neglect*abs(masking_vh(i,J,k))) vh(i,J,k) = 0.0
      endif
    enddo ; enddo ; enddo
    !$OMP do
    do j=js,je ; do k=1,nCat ; do i=is,ie
      h(i,j,k) = hin(i,j,k) - dt* G%IareaT(i,j) * &
           ((uh(I,j,k) - uh(I-1,j,k)) + (vh(i,J,k) - vh(i,J-1,k)))

      if (h(i,j,k) < 0.0) then
        call SIS_error(FATAL, 'Negative thickness encountered in ice_continuity().')
      endif
    enddo ; enddo ; enddo
    !$OMP end parallel
  elseif (x_first) then
  !    First, advect zonally.
    LB%ish = G%isc ; LB%ieh = G%iec
    LB%jsh = G%jsc-stensil ; LB%jeh = G%jec+stensil
    if (apply_h_neg) then
      call zonal_mass_flux(u, dt, G, US, IG, CS, LB, hin, uh, &
                           h_mobilize=CS%h_neglect_cont, masking_uh=masking_uh)
    else
      call zonal_mass_flux(u, dt, G, US, IG, CS, LB, hin, uh, masking_uh=masking_uh)
    endif

    call cpu_clock_begin(id_clock_update)
    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
      h(i,j,k) = hin(i,j,k) - dt* G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k))
      if (h(i,j,k) < 0.0) then
        call SIS_error(FATAL, &
        'Negative thickness encountered in u-pass of ice_continuity().')
      endif
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_update)

    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec

  !    Now advect meridionally, using the updated thicknesses to determine
  !  the fluxes.
    if (apply_h_neg) then
      call meridional_mass_flux(v, dt, G, US, IG, CS, LB, h, vh, &
                                h_mobilize=CS%h_neglect_cont, masking_vh=masking_vh)
    else
      call meridional_mass_flux(v, dt, G, US, IG, CS, LB, h, vh, masking_vh=masking_vh)
    endif

    call cpu_clock_begin(id_clock_update)
    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
      h(i,j,k) = h(i,j,k) - dt*G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k))
      if (h(i,j,k) < 0.0) then
        call SIS_error(FATAL, &
        'Negative thickness encountered in v-pass of ice_continuity().')
      endif
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_update)

  else  ! .not. x_first
  !    First, advect meridionally, so set the loop bounds accordingly.
    LB%ish = G%isc-stensil ; LB%ieh = G%iec+stensil
    LB%jsh = G%jsc ; LB%jeh = G%jec

    if (apply_h_neg) then
      call meridional_mass_flux(v, dt, G, US, IG, CS, LB, hin, vh, &
                                h_mobilize=CS%h_neglect_cont, masking_vh=masking_vh)
    else
      call meridional_mass_flux(v, dt, G, US, IG, CS, LB, hin, vh, masking_vh=masking_vh)
    endif

    call cpu_clock_begin(id_clock_update)
    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
      h(i,j,k) = hin(i,j,k) - dt*G%IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k))
      if (h(i,j,k) < 0.0) then
        call SIS_error(FATAL, &
        'Negative thickness encountered in v-pass of ice_continuity().')
      endif
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_update)

  !    Now advect zonally, using the updated thicknesses to determine
  !  the fluxes.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    if (apply_h_neg) then
      call zonal_mass_flux(u, dt, G, US, IG, CS, LB, h, uh, &
                           h_mobilize=CS%h_neglect_cont, masking_uh=masking_uh)
    else
      call zonal_mass_flux(u, dt, G, US, IG, CS, LB, h, uh, masking_uh=masking_uh)
    endif

    call cpu_clock_begin(id_clock_update)
    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
      h(i,j,k) = h(i,j,k) - dt* G%IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k))
      if (h(i,j,k) < 0.0) then
        call SIS_error(FATAL, &
        'Negative thickness encountered in u-pass of ice_continuity().')
      endif
    enddo ; enddo ; enddo
    call cpu_clock_end(id_clock_update)

  endif  ! End of x_first block.

end subroutine ice_continuity


!> ice_cover_transport advects the total fractional ice cover and limits them not to exceed 1.
subroutine ice_cover_transport(u, v, cvr, dt, G, US, IG, CS, masking_uhtot, masking_vhtot)
  type(SIS_hor_grid_type),           intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),               intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: u   !< Zonal ice velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: v   !< Meridional ice velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: cvr !< Fractional ice cover [nondim].
  real,                              intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_continuity_CS),           pointer       :: CS  !< The control structure returned by a
                                                          !! previous call to SIS_continuity_init.
  real, dimension(SZIB_(G),SZJ_(G)), &
                           optional, intent(in)    :: masking_uhtot !< If this zonal mass flux is 0, the
                                                          !! ice-cover u-transport is 0 [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           optional, intent(in)    :: masking_vhtot !< If this meridional mass flux is 0, the
                                                          !! ice-cover v-transport is 0 [R Z L2 T-1 ~> kg s-1].

  ! Local variables
  type(loop_bounds_type) :: LB  ! A structure with the active loop bounds.
  real, dimension(SZIB_(G),SZJ_(G)) :: ucvr ! Ice cover flux through zonal faces = u*cvr*dy [L2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJB_(G)) :: vcvr ! Ice cover flux through meridional faces = v*cvr*dx [L2 T-1 ~> m2 s-1].
  real    :: cvr_up
  integer :: is, ie, js, je, stensil
  integer :: i, j

  logical :: x_first
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_continuity: Module must be initialized before it is used.")
  x_first = (MOD(G%first_direction,2) == 0)

  stensil = 3 ; if (CS%simple_2nd) stensil = 2 ; if (CS%upwind_1st) stensil = 1

  do j=js,je ; do i=is,ie ; if (cvr(i,j) < 0.0) then
    call SIS_error(FATAL, 'Negative mass input to ice_cover_transport().')
  endif ; enddo ; enddo

  if (CS%use_upwind2d) then
    ! This reproduces the scheme that was originally used in SIS1.
    !$OMP parallel default(shared) private(cvr_up)
    !$OMP do
    do j=js,je ; do I=is-1,ie
      if (u(I,j) >= 0.0) then ; cvr_up = cvr(i,j)
      else ; cvr_up = cvr(i+1,j) ; endif
      ucvr(I,j) = G%dy_Cu(I,j) * u(I,j) * cvr_up
      if (present(masking_uhtot)) then ; if (masking_uhtot(I,j) == 0.0) ucvr(I,j) = 0.0 ; endif
    enddo ; enddo
    !$OMP do
    do J=js-1,je ; do i=is,ie
      if (v(i,J) >= 0.0) then ; cvr_up = cvr(i,j)
      else ; cvr_up = cvr(i,j+1) ; endif
      vcvr(i,J) = G%dx_Cv(i,J) * v(i,J) * cvr_up
      if (present(masking_vhtot)) then ; if (masking_vhtot(i,J) == 0.0) vcvr(i,J) = 0.0 ; endif
    enddo ; enddo
    !$OMP do
    do j=js,je ; do i=is,ie
      cvr(i,j) = cvr(i,j) - dt * G%IareaT(i,j) * &
           ((ucvr(I,j) - ucvr(I-1,j)) + (vcvr(i,J) - vcvr(i,J-1)))
      if (cvr(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative ice cover encountered in ice_cover_transport().')
    enddo ; enddo
    !$OMP end parallel
  elseif (x_first) then
    ! First, advect zonally.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc-stensil ; LB%jeh = G%jec+stensil
    call zonal_mass_flux(u, dt, G, US, IG, CS, LB, htot_in=cvr, uh_tot=ucvr, masking_uhtot=masking_uhtot)

    call cpu_clock_begin(id_clock_update)
    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      cvr(i,j) = cvr(i,j) - G%IareaT(i,j) * (dt*(ucvr(I,j) - ucvr(I-1,j)))
      if (cvr(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative ice cover encountered in u-pass of ice_cover_transport().')
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    ! Now advect meridionally, using the updated ice covers to determine the fluxes.
    call meridional_mass_flux(v, dt, G, US, IG, CS, LB, htot_in=cvr, vh_tot=vcvr, masking_vhtot=masking_vhtot)

    call cpu_clock_begin(id_clock_update)
    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      cvr(i,j) = max(1.0, cvr(i,j) - dt*G%IareaT(i,j) * (vcvr(i,J) - vcvr(i,J-1)))
      if (cvr(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative ice cover encountered in v-pass of ice_cover_transport().')
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

  else  ! .not. x_first
    !  First, advect meridionally, so set the loop bounds accordingly.
    LB%ish = G%isc-stensil ; LB%ieh = G%iec+stensil ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call meridional_mass_flux(v, dt, G, US, IG, CS, LB, htot_in=cvr, vh_tot=vcvr, masking_vhtot=masking_vhtot)

    call cpu_clock_begin(id_clock_update)
    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      cvr(i,j) = cvr(i,j) - dt*G%IareaT(i,j) * (vcvr(i,J) - vcvr(i,J-1))
      if (cvr(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative ice cover encountered in v-pass of ice_cover_transport().')
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

    ! Now advect zonally, using the updated ice covers to determine the fluxes.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call zonal_mass_flux(u, dt, G, US, IG, CS, LB, htot_in=cvr, uh_tot=ucvr, masking_uhtot=masking_uhtot)

    call cpu_clock_begin(id_clock_update)
    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      cvr(i,j) = max(1.0, cvr(i,j) - dt* G%IareaT(i,j) * (ucvr(I,j) - ucvr(I-1,j)))
      if (cvr(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative ice cover encountered in u-pass of ice_cover_transport().')
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

  endif  ! End of x_first block.

end subroutine ice_cover_transport


!> summed_continuity time steps the total ice, water, and snow mass changes summed across all the
!! thickness categories due to advection, using a monotonically limited, directionally split PPM
!! scheme or simple upwind 2-d scheme.  It may also update the ice thickness, using fluxes that are
!! proportional to the total fluxes times the ice mass divided by the total mass in the upwind cell.
subroutine summed_continuity(u, v, h_in, h, uh, vh, dt, G, US, IG, CS, h_ice)
  type(SIS_hor_grid_type),           intent(inout) :: G  !< The horizontal grid type
  type(ice_grid_type),               intent(inout) :: IG !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in)    :: u  !< Zonal ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(in)    :: v  !< Meridional ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)    :: h_in !< Initial total ice and snow mass per
                                                         !! unit cell area [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: h  !< Total ice and snow mass per unit cell
                                                         !! area [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), intent(out)   :: uh !< Total mass flux through zonal faces
                                                         !! = u*h*dy [R Z L2 T-1 ~> kg s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(out)   :: vh !< Total mass flux through meridional faces
                                                         !! = v*h*dx [R Z L2 T-1 ~> kg s-1]
  real,                              intent(in)    :: dt !< Time increment [T ~> s]
  type(unit_scale_type),             intent(in)    :: US !< A structure with unit conversion factors
  type(SIS_continuity_CS),           pointer       :: CS !< The control structure returned by a
                                                         !! previous call to SIS_continuity_init.
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(inout) :: h_ice  !< Total ice mass per unit cell
                                                         !! area [R Z ~> kg m-2].  h_ice must not exceed h.

  ! Local variables
  type(loop_bounds_type) :: LB  ! A structure with the active loop bounds.
  real, dimension(SZIB_(G),SZJ_(G)) :: uh_ice ! Ice mass flux through zonal faces = u*h*dy
                                              ! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G)) :: vh_ice ! Ice mass flux through meridional faces = v*h*dx
                                              ! [R Z L2 T-1 ~> kg s-1].
  real    :: h_up
  integer :: is, ie, js, je, stensil
  integer :: i, j

  logical :: x_first
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_continuity: Module must be initialized before it is used.")
  x_first = (MOD(G%first_direction,2) == 0)

  stensil = 3 ; if (CS%simple_2nd) stensil = 2 ; if (CS%upwind_1st) stensil = 1

  do j=js,je ; do i=is,ie ; if (h_in(i,j) < 0.0) then
    call SIS_error(FATAL, 'Negative mass input to summed_continuity().')
  endif ; enddo ; enddo

  if (present(h_ice)) then ; do j=js,je ; do i=is,ie ; if (h_ice(i,j) > h_in(i,j)) then
    call SIS_error(FATAL, 'ice mass exceeds total mass in summed_continuity().')
  endif ; enddo ; enddo ; endif

  if (CS%use_upwind2d) then
    ! This reproduces the scheme that was originally used in SIS1.
    !$OMP parallel default(shared) private(h_up)
    !$OMP do
    do j=js,je ; do I=is-1,ie
      if (u(I,j) >= 0.0) then ; h_up = h_in(i,j)
      else ; h_up = h_in(i+1,j) ; endif
      if (h_up < IG%CatIce*CS%h_neglect_cont) h_up = 0.0
      uh(I,j) = G%dy_Cu(I,j) * u(I,j) * h_up
    enddo ; enddo
    !$OMP do
    do J=js-1,je ; do i=is,ie
      if (v(i,J) >= 0.0) then ; h_up = h_in(i,j)
      else ; h_up = h_in(i,j+1) ; endif
      if (h_up < IG%CatIce*CS%h_neglect_cont) h_up = 0.0
      vh(i,J) = G%dx_Cv(i,J) * v(i,J) * h_up
    enddo ; enddo
    if (present(h_ice)) then
      !$OMP do
      do j=js,je ; do I=is-1,ie
        if (uh(I,j) < 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i+1,j) / h_in(i+1,j))
        elseif (uh(I,j) > 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i,j) / h_in(i,j))
        else ; uh_ice(I,j) = 0.0 ; endif
      enddo ; enddo
      !$OMP do
      do J=js-1,je ; do i=is,ie
        if (vh(i,J) < 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j+1) / h_in(i,j+1))
        elseif (vh(i,J) > 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j) / h_in(i,j))
        else ; vh_ice(i,J) = 0.0 ; endif
      enddo ; enddo
      !$OMP do
      do j=js,je ; do i=is,ie
        h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * &
             ((uh_ice(I,j) - uh_ice(I-1,j)) + (vh_ice(i,J) - vh_ice(i,J-1)))
      enddo ; enddo
    endif
    !$OMP do
    do j=js,je ; do i=is,ie
      h(i,j) = h_in(i,j) - (dt * G%IareaT(i,j)) * &
           ((uh(I,j) - uh(I-1,j)) + (vh(i,J) - vh(i,J-1)))
      ! if (h(i,j) < 0.0) call SIS_error(FATAL, &
      !   'Negative thickness encountered in ice_total_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() 2d.')
      ! endif ; endif
    enddo ; enddo
    !$OMP end parallel
  elseif (x_first) then
    ! First, advect zonally.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc-stensil ; LB%jeh = G%jec+stensil
    call zonal_mass_flux(u, dt, G, US, IG, CS, LB, htot_in=h_in, uh_tot=uh, h_mobilize=CS%h_neglect_cont)

    call cpu_clock_begin(id_clock_update)

    if (present(h_ice)) then
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh
        do I=LB%ish-1,LB%ieh
          if (uh(I,j) < 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i+1,j) / h_in(i+1,j))
          elseif (uh(I,j) > 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i,j) / h_in(i,j))
          else ; uh_ice(I,j) = 0.0 ; endif
        enddo
        do i=LB%ish,LB%ieh
          h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * (uh_ice(I,j) - uh_ice(I-1,j))
        enddo
      enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j) = h_in(i,j) - (dt * G%IareaT(i,j)) * (uh(I,j) - uh(I-1,j))
      ! if (h(i,j) < 0.0) call SIS_error(FATAL, &
      !   'Negative thickness encountered in u-pass of ice_total_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() x-1.')
      ! endif ; endif
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    ! Now advect meridionally, using the updated thicknesses to determine the fluxes.
    call meridional_mass_flux(v, dt, G, US, IG, CS, LB, htot_in=h, vh_tot=vh, h_mobilize=CS%h_neglect_cont)

    call cpu_clock_begin(id_clock_update)
    if (present(h_ice)) then
      !$OMP parallel do default(shared)
      do J=LB%jsh-1,LB%jeh ; do i=LB%ish,LB%ieh
        if (vh(i,J) < 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j+1) / h(i,j+1))
        elseif (vh(i,J) > 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j) / h(i,j))
        else ; vh_ice(i,J) = 0.0 ; endif
      enddo ; enddo
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
        h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * (vh_ice(i,J) - vh_ice(i,J-1))
      enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j) = h(i,j) - (dt * G%IareaT(i,j)) * (vh(i,J) - vh(i,J-1))
      if (h(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative thickness encountered in v-pass of summed_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() x-2.')
      ! endif ; endif
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

  else  ! .not. x_first
    !  First, advect meridionally, so set the loop bounds accordingly.
    LB%ish = G%isc-stensil ; LB%ieh = G%iec+stensil ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call meridional_mass_flux(v, dt, G, US, IG, CS, LB, htot_in=h_in, vh_tot=vh, h_mobilize=CS%h_neglect_cont)

    call cpu_clock_begin(id_clock_update)
    if (present(h_ice)) then
      !$OMP parallel do default(shared)
      do J=LB%jsh-1,LB%jeh ; do i=LB%ish,LB%ieh
        if (vh(i,J) < 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j+1) / h_in(i,j+1))
        elseif (vh(i,J) > 0.0) then ; vh_ice(i,J) = vh(i,J) * (h_ice(i,j) / h_in(i,j))
        else ; vh_ice(i,J) = 0.0 ; endif
      enddo ; enddo
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
        h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * (vh_ice(i,J) - vh_ice(i,J-1))
      enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j) = h_in(i,j) - (dt * G%IareaT(i,j)) * (vh(i,J) - vh(i,J-1))
      ! if (h(i,j) < 0.0) call SIS_error(FATAL, &
      !   'Negative thickness encountered in v-pass of ice_total_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() y-1.')
      ! endif ; endif
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

    !  Now advect zonally, using the updated thicknesses to determine the fluxes.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    call zonal_mass_flux(u, dt, G, US, IG, CS, LB, htot_in=h, uh_tot=uh, h_mobilize=CS%h_neglect_cont)

    call cpu_clock_begin(id_clock_update)

    if (present(h_ice)) then
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh
        do I=LB%ish-1,LB%ieh
          if (uh(I,j) < 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i+1,j) / h(i+1,j))
          elseif (uh(I,j) > 0.0) then ; uh_ice(I,j) = uh(I,j) * (h_ice(i,j) / h(i,j))
          else ; uh_ice(I,j) = 0.0 ; endif
        enddo
        do i=LB%ish,LB%ieh
          h_ice(i,j) = h_ice(i,j) - (dt * G%IareaT(i,j)) * (uh_ice(I,j) - uh_ice(I-1,j))
        enddo
      enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h(i,j) = h(i,j) - (dt * G%IareaT(i,j)) * (uh(I,j) - uh(I-1,j))
      if (h(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative thickness encountered in u-pass of summed_continuity().')
      ! if (present(h_ice)) then ; if (h_ice(i,j) > h(i,j)) then
      !   call SIS_error(FATAL, 'ice mass exceeds total mass in ice_total_continuity() y-2.')
      ! endif ; endif
    enddo ; enddo
    call cpu_clock_end(id_clock_update)

  endif  ! End of x_first block.

end subroutine summed_continuity

!> proportionate_continuity time steps the category thickness changes due to advection,
!! using input total mass fluxes with the fluxes proportionate to the relative upwind
!! thicknesses.
subroutine proportionate_continuity(h_tot_in, uh_tot, vh_tot, dt, G, US, IG, CS, &
                                    h1, uh1, vh1, h2, uh2, vh2, h3, uh3, vh3)
  type(SIS_hor_grid_type),        intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),            intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: h_tot_in !< Initial total ice and snow mass per unit
                                                       !! cell area [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: uh_tot !< Total mass flux through zonal faces
                                                       !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: vh_tot !< Total mass flux through meridional faces
                                                       !! [R Z L2 T-1 ~> kg s-1].
  real,                           intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type),          intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_continuity_CS),        pointer       :: CS  !< The control structure returned by a
                                                       !! previous call to SIS_continuity_init.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                        optional, intent(inout) :: h1  !< Updated mass of medium 1 (often ice) by
                                                       !! category [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                        optional, intent(out)   :: uh1 !< Zonal mass flux of medium 1 by category
                                                       !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                        optional, intent(out)   :: vh1 !< Meridional mass flux of medium 1 by category
                                                       !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                        optional, intent(inout) :: h2  !< Updated mass of medium 2 (often snow) by
                                                       !! category [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                        optional, intent(out)   :: uh2 !< Zonal mass flux of medium 2 by category
                                                       !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                        optional, intent(out)   :: vh2 !< Meridional mass flux of medium 2 by category
                                                       !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                        optional, intent(inout) :: h3  !< Updated mass of medium 3 (pond water?) by
                                                       !! category [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                        optional, intent(out)   :: uh3 !< Zonal mass flux of medium 3 by category
                                                       !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                        optional, intent(out)   :: vh3 !< Meridional mass flux of medium 3 by category
                                                       !! [R Z L2 T-1 ~> kg s-1].

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: h_tot  ! Total thicknesses [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G)) :: I_htot ! The Adcroft reciprocal of the total thicknesses [R-1 Z-1 ~> m2 kg-1].
  type(loop_bounds_type) :: LB  ! A structure with the active loop bounds.
  logical :: do_mask  ! If true, mask very small fluxes.
  integer :: is, ie, js, je, nCat, stensil
  integer :: i, j, k

  logical :: x_first
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nCat = IG%CatIce

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_continuity: Module must be initialized before it is used.")
  x_first = (MOD(G%first_direction,2) == 0)

  do_mask = ((CS%h_neglect_cont > 0.0) .or. (CS%frac_neglect > 0.0))

  do j=js,je ; do i=is,ie ; if (h_tot_in(i,j) < 0.0) then
    call SIS_error(FATAL, 'Negative thickness input to proportionate_continuity().')
  endif ; enddo ; enddo

  !$OMP parallel do default(shared)
  do j=js-1,je+1 ; do i=is-1,ie+1
    I_htot(i,j) = 0.0 ; if (h_tot_in(i,j) > 0.0) I_htot(i,j) = 1.0 / h_tot_in(i,j)
  enddo ; enddo

  if (CS%use_upwind2d) then
    ! Both directions are updated based on the original thicknesses.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec

    if (present(h1)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h1, uh1, G, IG, LB, frac_neglect=CS%frac_neglect)
        call merid_proportionate_fluxes(vh_tot, I_htot, h1, vh1, G, IG, LB, frac_neglect=CS%frac_neglect)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h1, uh1, G, IG, LB)
        call merid_proportionate_fluxes(vh_tot, I_htot, h1, vh1, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h1(i,j,k) = h1(i,j,k) - G%IareaT(i,j) * (dt * &
             ((uh1(I,j,k) - uh1(I-1,j,k)) + (vh1(i,J,k) - vh1(i,J-1,k))))
      enddo ; enddo ; enddo
    endif
    if (present(h2)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h2, uh2, G, IG, LB, masking_uh=uh1)
        call merid_proportionate_fluxes(vh_tot, I_htot, h2, vh2, G, IG, LB, masking_vh=vh1)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h2, uh2, G, IG, LB)
        call merid_proportionate_fluxes(vh_tot, I_htot, h2, vh2, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h2(i,j,k) = h2(i,j,k) - G%IareaT(i,j) * (dt * &
             ((uh2(I,j,k) - uh2(I-1,j,k)) + (vh2(i,J,k) - vh2(i,J-1,k))))
      enddo ; enddo ; enddo
    endif
    if (present(h3)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h3, uh3, G, IG, LB, masking_uh=uh1)
        call merid_proportionate_fluxes(vh_tot, I_htot, h3, vh3, G, IG, LB, masking_vh=vh1)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h3, uh3, G, IG, LB)
        call merid_proportionate_fluxes(vh_tot, I_htot, h3, vh3, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h3(i,j,k) = h3(i,j,k) - G%IareaT(i,j) * (dt * &
             ((uh3(I,j,k) - uh3(I-1,j,k)) + (vh3(i,J,k) - vh3(i,J-1,k))))
      enddo ; enddo ; enddo
    endif

  elseif (x_first) then
    ! First, advect zonally.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc-1 ; LB%jeh = G%jec+1
    if (present(h1)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h1, uh1, G, IG, LB, frac_neglect=CS%frac_neglect)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h1, uh1, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
        h1(i,j,k) = h1(i,j,k) - G%IareaT(i,j) * (dt * (uh1(I,j,k) - uh1(I-1,j,k)))
      enddo ; enddo ; enddo
    endif
    if (present(h2)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h2, uh2, G, IG, LB, masking_uh=uh1)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h2, uh2, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
        h2(i,j,k) = h2(i,j,k) - G%IareaT(i,j) * (dt * (uh2(I,j,k) - uh2(I-1,j,k)))
      enddo ; enddo ; enddo
    endif
    if (present(h3)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h3, uh3, G, IG, LB, masking_uh=uh1)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h3, uh3, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
        h3(i,j,k) = h3(i,j,k) - G%IareaT(i,j) * (dt * (uh3(I,j,k) - uh3(I-1,j,k)))
      enddo ; enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h_tot(i,j) = h_tot_in(i,j) - dt* G%IareaT(i,j) * (uh_tot(I,j) - uh_tot(I-1,j))
      if (h_tot(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative thickness encountered in u-pass of proportionate_continuity().')
      I_htot(i,j) = 0.0 ; if (h_tot(i,j) > 0.0) I_htot(i,j) = 1.0 / h_tot(i,j)
    enddo ; enddo

    !  Now advect meridionally, using the updated thicknesses to determine the fluxes.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    if (present(h1)) then
      if (do_mask) then
        call merid_proportionate_fluxes(vh_tot, I_htot, h1, vh1, G, IG, LB, frac_neglect=CS%frac_neglect)
      else
        call merid_proportionate_fluxes(vh_tot, I_htot, h1, vh1, G, IG, LB)
      endif
     !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h1(i,j,k) = h1(i,j,k) - G%IareaT(i,j) * (dt * (vh1(i,J,k) - vh1(i,J-1,k)) )
      enddo ; enddo ; enddo
    endif
    if (present(h2)) then
      if (do_mask) then
        call merid_proportionate_fluxes(vh_tot, I_htot, h2, vh2, G, IG, LB, masking_vh=vh1)
      else
        call merid_proportionate_fluxes(vh_tot, I_htot, h2, vh2, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h2(i,j,k) = h2(i,j,k) - G%IareaT(i,j) * (dt * (vh2(i,J,k) - vh2(i,J-1,k)) )
      enddo ; enddo ; enddo
    endif
    if (present(h3)) then
      if (do_mask) then
        call merid_proportionate_fluxes(vh_tot, I_htot, h3, vh3, G, IG, LB, masking_vh=vh1)
      else
        call merid_proportionate_fluxes(vh_tot, I_htot, h3, vh3, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h3(i,j,k) = h3(i,j,k) - G%IareaT(i,j) * (dt * (vh3(i,J,k) - vh3(i,J-1,k)) )
      enddo ; enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h_tot(i,j) = h_tot(i,j) - dt* G%IareaT(i,j) * (vh_tot(i,J) - vh_tot(i,J-1))
      if (h_tot(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative thickness encountered in v-pass of proportionate_continuity().')
      ! I_htot(i,j) = 0.0 ; if (h_tot(i,j) > 0.0) I_htot(i,j) = 1.0 / h_tot(i,j)
    enddo ; enddo

  else  ! .not. x_first
    ! First, advect meridionally, so set the loop bounds accordingly.
    LB%ish = G%isc-1 ; LB%ieh = G%iec+1 ; LB%jsh = G%jsc ; LB%jeh = G%jec

    if (present(h1)) then
      if (do_mask) then
        call merid_proportionate_fluxes(vh_tot, I_htot, h1, vh1, G, IG, LB, frac_neglect=CS%frac_neglect)
      else
        call merid_proportionate_fluxes(vh_tot, I_htot, h1, vh1, G, IG, LB)
      endif
     !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h1(i,j,k) = h1(i,j,k) - G%IareaT(i,j) * (dt * (vh1(i,J,k) - vh1(i,J-1,k)) )
      enddo ; enddo ; enddo
    endif
    if (present(h2)) then
      if (do_mask) then
        call merid_proportionate_fluxes(vh_tot, I_htot, h2, vh2, G, IG, LB, masking_vh=vh1)
      else
        call merid_proportionate_fluxes(vh_tot, I_htot, h2, vh2, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h2(i,j,k) = h2(i,j,k) - G%IareaT(i,j) * (dt * (vh2(i,J,k) - vh2(i,J-1,k)) )
      enddo ; enddo ; enddo
    endif
    if (present(h3)) then
      if (do_mask) then
        call merid_proportionate_fluxes(vh_tot, I_htot, h3, vh3, G, IG, LB, frac_neglect=CS%frac_neglect)
      else
        call merid_proportionate_fluxes(vh_tot, I_htot, h3, vh3, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=js,je ; do k=1,nCat ; do i=is,ie
        h3(i,j,k) = h3(i,j,k) - G%IareaT(i,j) * (dt * (vh3(i,J,k) - vh3(i,J-1,k)) )
      enddo ; enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h_tot(i,j) = h_tot(i,j) - dt* G%IareaT(i,j) * (vh_tot(i,J) - vh_tot(i,J-1))
      if (h_tot(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative thickness encountered in v-pass of proportionate_continuity().')
      I_htot(i,j) = 0.0 ; if (h_tot(i,j) > 0.0) I_htot(i,j) = 1.0 / h_tot(i,j)
    enddo ; enddo

    ! Now advect zonally, using the updated thicknesses to determine the fluxes.
    LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
    if (present(h1)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h1, uh1, G, IG, LB, frac_neglect=CS%frac_neglect)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h1, uh1, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
        h1(i,j,k) = h1(i,j,k) - G%IareaT(i,j) * (dt * (uh1(I,j,k) - uh1(I-1,j,k)))
      enddo ; enddo ; enddo
    endif
    if (present(h2)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h2, uh2, G, IG, LB, masking_uh=uh1)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h2, uh2, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
        h2(i,j,k) = h2(i,j,k) - G%IareaT(i,j) * (dt * (uh2(I,j,k) - uh2(I-1,j,k)))
      enddo ; enddo ; enddo
    endif
    if (present(h3)) then
      if (do_mask) then
        call zonal_proportionate_fluxes(uh_tot, I_htot, h3, uh3, G, IG, LB, masking_uh=uh1)
      else
        call zonal_proportionate_fluxes(uh_tot, I_htot, h3, uh3, G, IG, LB)
      endif
      !$OMP parallel do default(shared)
      do j=LB%jsh,LB%jeh ; do k=1,nCat ; do i=LB%ish,LB%ieh
        h3(i,j,k) = h3(i,j,k) - G%IareaT(i,j) * (dt * (uh3(I,j,k) - uh3(I-1,j,k)))
      enddo ; enddo ; enddo
    endif

    !$OMP parallel do default(shared)
    do j=LB%jsh,LB%jeh ; do i=LB%ish,LB%ieh
      h_tot(i,j) = h_tot_in(i,j) - dt* G%IareaT(i,j) * (uh_tot(I,j) - uh_tot(I-1,j))
      if (h_tot(i,j) < 0.0) call SIS_error(FATAL, &
        'Negative thickness encountered in u-pass of proportionate_continuity().')
      ! I_htot(i,j) = 0.0 ; if (h_tot(i,j) > 0.0) I_htot(i,j) = 1.0 / h_tot(i,j)
    enddo ; enddo

  endif  ! End of x_first block.

end subroutine proportionate_continuity

!> Calculate zonal fluxes by category that are proportionate to the relative masses in the upwind cell.
subroutine zonal_proportionate_fluxes(uh_tot, I_htot, h, uh, G, IG, LB, masking_uh, frac_neglect)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: uh_tot !< Total mass flux through zonal faces
                                                !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: I_htot !< Adcroft reciprocal of the total mass per unit
                                                !! cell area [R-1 Z-1 ~> m2 kg-1].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(inout) :: h   !< Mass per unit cell area by category [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(out)   :: uh  !< Category mass flux through zonal faces = u*h*dy.
                                                !! [R Z L2 T-1 ~> kg s-1].
  type(loop_bounds_type),  intent(in)    :: LB  !< A structure with the active loop bounds.
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: masking_uh  !< If this is 0, uh = 0.  Often this is another
                                                !! zonal mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
  real,          optional, intent(in)    :: frac_neglect  !< Any category with less than this fraction
                                                !! of the total transport has no transport [nondim].

  ! Local variables
  real :: frac ! The fraction of the total transport in a category [nondim].
  integer :: i, j, k, ish, ieh, jsh, jeh, nCat

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nCat = IG%CatIce
  !$OMP parallel do default(shared)
  do j=jsh,jeh ; do k=1,nCat ; do I=ish-1,ieh
    frac = 0.0
    if (uh_tot(I,j) < 0.0) then ; frac = (h(i+1,j,k) * I_htot(i+1,j))
    elseif (uh_tot(I,j) > 0.0) then ; frac = (h(i,j,k) * I_htot(i,j)) ; endif
    if (present(masking_uh)) then ; if (masking_uh(I,j,k) == 0.0) frac = 0.0 ; endif
    if (present(frac_neglect)) then ; if (frac < frac_neglect) frac = 0.0 ; endif

    uh(I,j,k) = frac * uh_tot(I,j)
  enddo ; enddo ; enddo

end subroutine zonal_proportionate_fluxes

!> Calculate meridional mass fluxes by category that are proportionate to the relative masses in the upwind cell.
subroutine merid_proportionate_fluxes(vh_tot, I_htot, h, vh, G, IG, LB, masking_vh, frac_neglect)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: vh_tot !< Total mass flux through meridional faces
                                                !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: I_htot !< Adcroft reciprocal of the total mass per unit
                                                !! cell area [R-1 Z-1 ~> m2 kg-1].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(inout) :: h   !< Mass per unit cell area by category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                           intent(out)   :: vh  !< Category mass flux through meridional faces = v*h*dx
                                                !! [R Z L2 T-1 ~> kg s-1].
  type(loop_bounds_type),  intent(in)    :: LB  !< A structure with the active loop bounds.
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: masking_vh  !< If this is 0, vh = 0.  Often this is another
                                                !! meridional mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
  real,          optional, intent(in)    :: frac_neglect  !< Any category with less than this fraction
                                                !! of the total transport has no transport [nondim].

  ! Local variables
  real :: frac ! The fraction of the total transport in a category [nondim].
  integer :: i, j, k, ish, ieh, jsh, jeh, nCat

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nCat = IG%CatIce
  !$OMP parallel do default(shared)
  do J=jsh-1,jeh ; do k=1,nCat ; do i=ish,ieh
    frac = 0.0
    if (vh_tot(i,J) < 0.0) then ; frac = (h(i,j+1,k) * I_htot(i,j+1))
    elseif (vh_tot(i,J) > 0.0) then ; frac = (h(i,j,k) * I_htot(i,j)) ; endif
    if (present(masking_vh)) then ; if (masking_vh(i,J,k) == 0.0) frac = 0.0 ; endif
    if (present(frac_neglect)) then ; if (frac < frac_neglect) frac = 0.0 ; endif

    vh(i,J,k) = frac * vh_tot(i,J)
  enddo ; enddo ; enddo

end subroutine merid_proportionate_fluxes

!> Calculates the mass or volume fluxes through the zonal
!! faces, and other related quantities.
subroutine zonal_mass_flux(u, dt, G, US, IG, CS, LB, h_in, uh, htot_in, uh_tot, h_mobilize, &
                           masking_uh, masking_uhtot)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: u   !< Zonal ice velocity [L T-1 ~> m s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_continuity_CS), pointer       :: CS  !< The control structure returned by a
                                                !! previous call to SIS_continuity_init.
  type(loop_bounds_type),  intent(in)    :: LB  !< A structure with the active loop bounds.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: h_in !< Category thickness used to calculate the fluxes [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(out)   :: uh  !< Category volume flux through zonal faces = u*h*dy
                                                !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: htot_in !< Total thicknesses used to calculate the fluxes [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(out)   :: uh_tot !< Total mass flux through zonal faces = u*htot*dy
                                                !! [R Z L2 T-1 ~> kg s-1].
  real,          optional, intent(in)    :: h_mobilize !< The minimum ice thickness per category that
                                                !! is able to move out of a cell [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: masking_uh  !< If this is 0, uh = 0.  Often this is another
                                                !! zonal mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: masking_uhtot !< If this is 0, uh_tot = 0.  Often this is another
                                                !! total zonal mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
!   This subroutine calculates the mass or volume fluxes through the zonal
! faces, and other related quantities.

  ! Local variables
!  real, dimension(SZIB_(G)) :: &
!    duhdu      ! Partial derivative of uh with u [R Z L ~> kg m-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &    ! The total thickness summed across categories [R Z ~> kg m-2].
    I_htot, &  ! The inverse of htot or 0 [R-1 Z-1 ~> m2 kg-1].
    hl, hr      ! Left and right face thicknesses [R Z ~> kg m-2].
  real, dimension(SZIB_(G)) :: &
    uhtot      ! The total transports [R Z L2 T-1 ~> kg s-1].
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim].
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
!  real :: h_marg ! The marginal thickness of a flux [R Z ~> kg m-2].
!  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west [L ~> m].
  integer :: i, j, k, ish, ieh, jsh, jeh, nCat

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nCat = IG%CatIce

  call cpu_clock_begin(id_clock_update)

  htot(:,:) = 0.0
  if (present(htot_in)) then
    if (present(h_mobilize)) then
      !$OMP parallel do default(shared)
      do j=jsh,jeh ; do i=G%isd,G%ied
        if (htot_in(i,j) >= nCat*h_mobilize) htot(i,j) = htot(i,j) + htot_in(i,j)
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=jsh,jeh ; do i=G%isd,G%ied
        htot(i,j) = htot(i,j) + htot_in(i,j)
      enddo ; enddo
    endif
  elseif (present(h_in)) then
    if (present(h_mobilize)) then
      !$OMP parallel do default(shared)
      do j=jsh,jeh ; do k=1,nCat ; do i=G%isd,G%ied
        if (h_in(i,j,k) >= h_mobilize) htot(i,j) = htot(i,j) + h_in(i,j,k)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=jsh,jeh ; do k=1,nCat ; do i=G%isd,G%ied
        htot(i,j) = htot(i,j) + h_in(i,j,k)
      enddo ; enddo ; enddo
    endif
  else
    call SIS_error(FATAL, "Either h_in or htot_in must be present in call to zonal_mass_flux.")
  endif

  ! This sets hl and hr.
  if (CS%upwind_1st) then
    do j=jsh,jeh ; do i=ish-1,ieh+1
      hl(i,j) = htot(i,j) ; hr(i,j) = htot(i,j)
    enddo ; enddo
  else
    call PPM_reconstruction_x(htot, hl, hr, G, LB, 0.0, CS%monotonic, &
                              simple_2nd=CS%simple_2nd)
  endif
  call cpu_clock_end(id_clock_update)

  call cpu_clock_begin(id_clock_correct)

  !$OMP parallel do default(shared) private(uhtot)
  do j=jsh,jeh
    ! Set uhtot and duhdu.
    do I=ish-1,ieh
      ! Set new values of uh and duhdu.
      if (u(I,j) > 0.0) then
        if (CS%vol_CFL) then ; CFL = (u(I,j) * dt) * (G%dy_Cu(I,j) * G%IareaT(i,j))
        else ; CFL = u(I,j) * dt * G%IdxT(i,j) ; endif
        curv_3 = hL(i,j) + hR(i,j) - 2.0*htot(i,j)
        uhtot(I) = G%dy_Cu(I,j) * u(I,j) * &
            (hR(i,j) + CFL * (0.5*(hL(i,j) - hR(i,j)) + curv_3*(CFL - 1.5)))
!        h_marg = hR(i,j) + CFL * ((hL(i,j) - hR(i,j)) + 3.0*curv_3*(CFL - 1.0))
      elseif (u(I,j) < 0.0) then
        if (CS%vol_CFL) then ; CFL = (-u(I,j) * dt) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
        else ; CFL = -u(I,j) * dt * G%IdxT(i+1,j) ; endif
        curv_3 = hL(i+1,j) + hR(i+1,j) - 2.0*htot(i+1,j)
        uhtot(I) = G%dy_Cu(I,j) * u(I,j) * &
            (hL(i+1,j) + CFL * (0.5*(hR(i+1,j)-hL(i+1,j)) + curv_3*(CFL - 1.5)))
!        h_marg = hL(i+1) + CFL * ((hR(i+1,j)-hL(i+1,j)) + 3.0*curv_3*(CFL - 1.0))
      else
        uhtot(I) = 0.0
!        h_marg = 0.5 * (hl(i+1,j) + hr(i,j))
      endif
!      duhdu(I,j) = G%dy_Cu(I,j) * h_marg ! * visc_rem(I)
    enddo

    ! Partition the transports by category in proportion to their relative masses.
    if (present(uh)) then
      do i=ish-1,ieh+1
        I_htot(i,j) = 0.0 ; if (htot(i,j) > 0.0) I_htot(i,j) = 1.0 / htot(i,j)
      enddo
      if (present(h_mobilize)) then
        do k=1,nCat ; do I=ish-1,ieh
          uh(I,j,k) = 0.0
          if (u(I,j) >= 0.0) then
            if (h_in(i,j,k) >= h_mobilize) uh(I,j,k) = uhtot(I) * (h_in(i,j,k) * I_htot(i,j))
          else
            if (h_in(i+1,j,k) >= h_mobilize) uh(I,j,k) = uhtot(I) * (h_in(i+1,j,k) * I_htot(i+1,j))
          endif
        enddo ; enddo
      else
        do k=1,nCat ; do I=ish-1,ieh
          if (u(I,j) >= 0.0) then
            uh(I,j,k) = uhtot(I) * (h_in(i,j,k) * I_htot(i,j))
          else
            uh(I,j,k) = uhtot(I) * (h_in(i+1,j,k) * I_htot(i+1,j))
          endif
        enddo ; enddo
      endif
    endif

    ! Block mass fluxes in categories where a related flux (e.g. of ice) is zero.
    if (present(masking_uh) .and. present(uh)) then
      do k=1,nCat ; do I=ish-1,ieh
        if (masking_uh(I,j,k) == 0.0) uh(I,j,k) = 0.0
        if (abs(uh(I,j,k)) < CS%frac_neglect*abs(masking_uh(I,j,k))) uh(I,j,k) = 0.0
      enddo ; enddo
    endif

    if (present(uh_tot) .and. present(masking_uhtot)) then
      do I=ish-1,ieh
        uh_tot(I,j) = uhtot(I)
        if (masking_uhtot(I,j) == 0.0) uh_tot(I,j) = 0.0
      enddo
    elseif (present(uh_tot)) then
      do I=ish-1,ieh
        uh_tot(I,j) = uhtot(I)
      enddo
    endif

  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

end subroutine zonal_mass_flux

!> Calculates the mass or volume fluxes through the meridional
!! faces, and other related quantities.
subroutine meridional_mass_flux(v, dt, G, US, IG, CS, LB, h_in, vh, htot_in, vh_tot, &
                                h_mobilize, masking_vh, masking_vhtot)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: v   !< Meridional ice velocity [L T-1 ~> m s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s]
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_continuity_CS), pointer       :: CS  !< The control structure returned by a
                                                !! previous call to SIS_continuity_init.
  type(loop_bounds_type),  intent(in)    :: LB  !< A structure with the active loop bounds.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: h_in !< Category thickness used to calculate the fluxes [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                 optional, intent(out)   :: vh  !< Category volume flux through meridional faces = v*h*dx
                                                !! [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G)), &
                 optional, intent(in)    :: htot_in !< Total thicknesses used to calculate the fluxes [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(out)   :: vh_tot !< Total mass flux through meridional faces = v*htot*dx
                                                !! [R Z L2 T-1 ~> kg s-1].
  real,          optional, intent(in)    :: h_mobilize !< The minimum ice thickness per category that
                                                !! is able to move out of a cell [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                 optional, intent(in)    :: masking_vh  !< If this is 0, vh = 0.  Often this is another
                                                !! meridional mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(in)    :: masking_vhtot !< If this is 0, vh_tot = 0.  Often this is another
                                                !! total meridional mass flux, e.g. of ice [R Z L2 T-1 ~> kg s-1].

!   This subroutine calculates the mass or volume fluxes through the meridional
! faces, and other related quantities.

  ! Local variables
  real, dimension(SZI_(G)) :: &
    dvhdv      ! Partial derivative of vh with v [R Z L ~> kg m-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &    ! The total thickness summed across categories [R Z ~> kg m-2].
    I_htot, &  ! The inverse of htot or 0 [R-1 Z-1 ~> m2 kg-1].
    hl, hr     ! Left and right face thicknesses [R Z ~> kg m-2].
  real, dimension(SZI_(G)) :: &
    vhtot      ! The total transports [R Z L2 s-1 ~> kg s-1].
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim].
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  real :: h_marg ! The marginal thickness of a flux [R Z ~> kg m-2].
!  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south [L ~> m].
  integer :: i, j, k, ish, ieh, jsh, jeh, nCat

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nCat = IG%CatIce

  call cpu_clock_begin(id_clock_update)

  htot(:,:) = 0.0

  if (present(htot_in)) then
    if (present(h_mobilize)) then
      !$OMP parallel do default(shared)
      do j=G%jsd,G%jed ; do i=ish,ieh
        if (htot_in(i,j) >= nCat*h_mobilize) htot(i,j) = htot(i,j) + htot_in(i,j)
      enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=G%jsd,G%jed ; do i=ish,ieh
        htot(i,j) = htot(i,j) + htot_in(i,j)
      enddo ; enddo
    endif
  elseif (present(h_in)) then
    if (present(h_mobilize)) then
      !$OMP parallel do default(shared)
      do j=G%jsd,G%jed ; do k=1,nCat ; do i=ish,ieh
        if (h_in(i,j,k) >= h_mobilize) htot(i,j) = htot(i,j) + h_in(i,j,k)
      enddo ; enddo ; enddo
    else
      !$OMP parallel do default(shared)
      do j=G%jsd,G%jed ; do k=1,nCat ; do i=ish,ieh
        htot(i,j) = htot(i,j) + h_in(i,j,k)
      enddo ; enddo ; enddo
    endif
  else
    call SIS_error(FATAL, "Either h_in or htot_in must be present in call to meridional_mass_flux.")
  endif
  if (present(vh)) then ; do j=jsh-1,jeh+1 ; do i=ish,ieh
    I_htot(i,j) = 0.0 ; if (htot(i,j) > 0.0) I_htot(i,j) = 1.0 / htot(i,j)
  enddo ; enddo ; endif

  ! This sets hl and hr.
  if (CS%upwind_1st) then
    do j=jsh-1,jeh+1 ; do i=ish,ieh
      hl(i,j) = htot(i,j) ; hr(i,j) = htot(i,j)
    enddo ; enddo
  else
    call PPM_reconstruction_y(htot, hl, hr, G, LB, 0.0, CS%monotonic, &
                              simple_2nd=CS%simple_2nd)
  endif
  call cpu_clock_end(id_clock_update)

  call cpu_clock_begin(id_clock_correct)
  !$OMP parallel do default(shared) private(vhtot)
  do J=jsh-1,jeh
    ! This sets vh and dvhdv.
    do i=ish,ieh
      if (v(i,J) > 0.0) then
        if (CS%vol_CFL) then ; CFL = (v(i,J) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j))
        else ; CFL = v(i,J) * dt * G%IdyT(i,j) ; endif
        curv_3 = hL(i,j) + hR(i,j) - 2.0*htot(i,j)
        vhtot(i) = G%dx_Cv(i,J) * v(i,J) * ( hR(i,j) + CFL * &
            (0.5*(hL(i,j) - hR(i,j)) + curv_3*(CFL - 1.5)) )
       ! h_marg = hR(i,j) + CFL * ((hL(i,j) - hR(i,j)) + 3.0*curv_3*(CFL - 1.0))
      elseif (v(i,J) < 0.0) then
        if (CS%vol_CFL) then ; CFL = (-v(i,J) * dt) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
        else ; CFL = -v(i,J) * dt *  G%IdyT(i,j+1) ; endif
        curv_3 = hL(i,j+1) + hR(i,j+1) - 2.0*htot(i,j+1)
        vhtot(i) = G%dx_Cv(i,J) * v(i,J) * ( hL(i,j+1) + CFL * &
            (0.5*(hR(i,j+1)-hL(i,j+1)) + curv_3*(CFL - 1.5)) )
       ! h_marg = hL(i,j+1) + CFL * ((hR(i,j+1)-hL(i,j+1)) + 3.0*curv_3*(CFL - 1.0))
      else
        vhtot(i) = 0.0
        ! h_marg = 0.5 * (hl(i,j+1) + hr(i,j))
      endif
      ! dvhdv(i) = G%dx_Cv(i,J) * h_marg ! * visc_rem(i)
    enddo

    ! Partition the transports by category in proportion to their relative masses.
    if (present(vh)) then
      if (present(h_mobilize)) then
        do k=1,nCat ; do i=ish,ieh
          vh(i,J,k) = 0.0
          if (v(i,J) >= 0.0) then
            if (h_in(i,j,k) >= h_mobilize) vh(i,J,k) = vhtot(i) * (h_in(i,j,k) * I_htot(i,j))
          else
            if (h_in(i,j+1,k) >= h_mobilize) vh(i,J,k) = vhtot(i) * (h_in(i,j+1,k) * I_htot(i,j+1))
          endif
        enddo ; enddo
      else
        do k=1,nCat ; do i=ish,ieh
          if (v(i,J) >= 0.0) then
            vh(i,J,k) = vhtot(i) * (h_in(i,j,k) * I_htot(i,j))
          else
            vh(i,J,k) = vhtot(i) * (h_in(i,j+1,k) * I_htot(i,j+1))
          endif
        enddo ; enddo
      endif
    endif

    ! Block mass fluxes in categories where a related flux (e.g. of ice) is zero.
    if (present(masking_vh) .and. present(vh)) then
      do k=1,nCat ; do i=ish,ieh
        if (masking_vh(i,J,k) == 0.0) vh(i,J,k) = 0.0
        if (abs(vh(i,J,k)) < CS%frac_neglect*abs(masking_vh(i,J,k))) vh(i,J,k) = 0.0
      enddo ; enddo
    endif

    if (present(vh_tot) .and. present(masking_vhtot)) then
      do i=ish,ieh
        vh_tot(i,J) = vhtot(I)
        if (masking_vhtot(i,J) == 0.0) vh_tot(i,J) = 0.0
      enddo
    elseif (present(vh_tot)) then
      do i=ish,ieh
        vh_tot(i,J) = vhtot(I)
      enddo
    endif

  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

end subroutine meridional_mass_flux

!> Calculate a piecewise parabolic thickness reconstruction in the x-direction.
subroutine PPM_reconstruction_x(h_in, h_l, h_r, G, LB, h_min, monotonic, simple_2nd)
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Initial thickness of a category [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l !< Left edge value of thickness reconstruction [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r !< Right edge value of thickness reconstruction [R Z ~> kg m-2]
  type(loop_bounds_type),           intent(in)  :: LB  !< A structure with the active loop bounds.
  real,                             intent(in)  :: h_min !< The minimum thickness that can be
                                                       !! obtained by a concave parabolic fit [R Z ~> kg m-2].
  logical, optional,                intent(in)  :: monotonic !< If true, use the Colella & Woodward monotonic limiter.
                                                       !! Otherwise use a simple positive-definite limiter.
  logical, optional,                intent(in)  :: simple_2nd !< If true, use the arithmetic mean thicknesses as the
                                                       !! default edge values for a simple 2nd order scheme.
! This subroutine calculates left/right edge values for PPM reconstruction.

  ! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_ip1, h_im1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, stensil

  use_CW84 = .false. ; if (present(monotonic)) use_CW84 = monotonic
  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish-1 ; iel = LB%ieh+1 ; jsl = LB%jsh ; jel = LB%jeh

  ! This is the stensil of the reconstruction, not the scheme overall.
  stensil = 2 ; if (use_2nd) stensil = 1

  if ((isl-stensil < G%isd) .or. (iel+stensil > G%ied)) then
    write(mesg,'("In SIS_continuity, PPM_reconstruction_x called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               stensil + max(G%isd-isl,iel-G%ied)
    call SIS_error(FATAL,mesg)
  endif
  if ((jsl < G%jsd) .or. (jel > G%jed)) then
    write(mesg,'("In SIS_continuity, PPM_reconstruction_x called with a ", &
               & "y-halo that needs to be increased by ",i2,".")') &
               max(G%jsd-jsl,jel-G%jed)
    call SIS_error(FATAL,mesg)
  endif

  if (use_2nd) then
    !$OMP parallel do default(shared) private(h_im1,h_ip1)
    do j=jsl,jel ; do i=isl,iel
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) )
    enddo ; enddo
  else
    !$OMP parallel do default(shared) private(dMx,dMn,h_im1,h_ip1)
    do j=jsl,jel
      do i=isl-1,iel+1
        if ((G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j)) == 0.0) then
          slp(i,j) = 0.0
        else
          ! This uses a simple 2nd order slope.
          slp(i,j) = 0.5 * (h_in(i+1,j) - h_in(i-1,j))
          ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
          dMx = max(h_in(i+1,j), h_in(i-1,j), h_in(i,j)) - h_in(i,j)
          dMn = h_in(i,j) - min(h_in(i+1,j), h_in(i-1,j), h_in(i,j))
          slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                  ! * (G%mask2dT(i-1,j) * G%mask2dT(i,j) * G%mask2dT(i+1,j))
        endif
      enddo

      do i=isl,iel
        ! Neighboring values should take into account any boundaries.  The 3
        ! following sets of expressions are equivalent.
      ! h_im1 = h_in(i-1,j,k) ; if (G%mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
      ! h_ip1 = h_in(i+1,j,k) ; if (G%mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
        h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
        h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
        ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
        h_l(i,j) = 0.5*( h_im1 + h_in(i,j) ) + oneSixth*( slp(i-1,j) - slp(i,j) )
        h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i+1,j) )
      enddo
    enddo
  endif

  if (use_CW84) then
    call PPM_limit_CW84(h_in, h_l, h_r, G, isl, iel, jsl, jel)
  else
    call PPM_limit_pos(h_in, h_l, h_r, h_min, G, isl, iel, jsl, jel)
  endif

  return
end subroutine PPM_reconstruction_x

!> Calculate a piecewise parabolic thickness reconstruction in the y-direction.
subroutine PPM_reconstruction_y(h_in, h_l, h_r, G, LB, h_min, monotonic, simple_2nd)
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Initial thickness of a category [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l !< Left edge value of thickness reconstruction [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r !< Right edge value of thickness reconstruction [R Z ~> kg m-2]
  type(loop_bounds_type),           intent(in)  :: LB  !< A structure with the active loop bounds.
  real,                             intent(in)  :: h_min !< The minimum thickness that can be
                                                       !! obtained by a concave parabolic fit [R Z ~> kg m-2].
  logical, optional,                intent(in)  :: monotonic !< If true, use the Colella & Woodward monotonic limiter.
                                                       !! Otherwise use a simple positive-definite limiter.
  logical, optional,                intent(in)  :: simple_2nd !< If true, use the arithmetic mean thicknesses as the
                                                       !! default edge values for a simple 2nd order scheme.
! This subroutine calculates left/right edge values for PPM reconstruction.

! Local variables with useful mnemonic names.
  real, dimension(SZI_(G),SZJ_(G))  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_jp1, h_jm1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  character(len=256) :: mesg
  integer :: i, j, isl, iel, jsl, jel, stensil

  use_CW84 = .false. ; if (present(monotonic)) use_CW84 = monotonic
  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  isl = LB%ish ; iel = LB%ieh ; jsl = LB%jsh-1 ; jel = LB%jeh+1

  ! This is the stensil of the reconstruction, not the scheme overall.
  stensil = 2 ; if (use_2nd) stensil = 1

  if ((isl < G%isd) .or. (iel > G%ied)) then
    write(mesg,'("In SIS_continuity, PPM_reconstruction_y called with a ", &
               & "x-halo that needs to be increased by ",i2,".")') &
               max(G%isd-isl,iel-G%ied)
    call SIS_error(FATAL,mesg)
  endif
  if ((jsl-stensil < G%jsd) .or. (jel+stensil > G%jed)) then
    write(mesg,'("In SIS_continuity, PPM_reconstruction_y called with a ", &
                 & "y-halo that needs to be increased by ",i2,".")') &
                 stensil + max(G%jsd-jsl,jel-G%jed)
    call SIS_error(FATAL,mesg)
  endif

  if (use_2nd) then
    !$OMP parallel do default(shared) private(h_jm1,h_jp1)
    do j=jsl,jel ; do i=isl,iel
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) )
    enddo ; enddo
  else
    !$OMP parallel do default(shared)  private(dMx,dMn)
    do j=jsl-1,jel+1 ; do i=isl,iel
      if ((G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1)) == 0.0) then
        slp(i,j) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j) = 0.5 * (h_in(i,j+1) - h_in(i,j-1))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i,j+1), h_in(i,j-1), h_in(i,j)) - h_in(i,j)
        dMn = h_in(i,j) - min(h_in(i,j+1), h_in(i,j-1), h_in(i,j))
        slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dMx, dMn))
                ! * (G%mask2dT(i,j-1) * G%mask2dT(i,j) * G%mask2dT(i,j+1))
      endif
    enddo ; enddo
    !$OMP parallel do default(shared) private(h_jm1,h_jp1)
    do j=jsl,jel ; do i=isl,iel
      ! Neighboring values should take into account any boundaries.
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) ) + oneSixth*( slp(i,j-1) - slp(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) ) + oneSixth*( slp(i,j) - slp(i,j+1) )
    enddo ; enddo
  endif

  if (use_CW84) then
    call PPM_limit_CW84(h_in, h_l, h_r, G, isl, iel, jsl, jel)
  else
    call PPM_limit_pos(h_in, h_l, h_r, h_min, G, isl, iel, jsl, jel)
  endif

  return
end subroutine PPM_reconstruction_y

!> Limit the left or right edge values of the PPM reconstruction to give a
!! reconstruction that is positive-definite.
subroutine PPM_limit_pos(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
  type(SIS_hor_grid_type),          intent(in)    :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: h_in !< Initial thickness of a category [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: h_L !< Left edge value of thickness reconstruction [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: h_R !< Right edge value of thickness reconstruction [R Z ~> kg m-2]
  real,                             intent(in)    :: h_min !< The minimum thickness that can be
                                                           !! obtained by a concave parabolic fit [R Z ~> kg m-2].
  integer,                          intent(in)    :: iis !< The starting i-index to work on
  integer,                          intent(in)    :: iie !< The ending i-index to work on
  integer,                          intent(in)    :: jis !< The starting j-index to work on
  integer,                          intent(in)    :: jie !< The ending j-index to work on
! This subroutine limits the left/right edge values of the PPM reconstruction
! to give a reconstruction that is positive-definite.  Here this is
! reinterpreted as giving a constant thickness if the mean thickness is less
! than h_min, with a minimum of h_min otherwise.

! Local variables
  real    :: curv, dh, scale
  character(len=256) :: mesg
  integer :: i,j

  do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*(h_L(i,j) + h_R(i,j) - 2.0*h_in(i,j))
    if (curv > 0.0) then ! Only minima are limited.
      dh = h_R(i,j) - h_L(i,j)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
        if (h_in(i,j) <= h_min) then
          h_L(i,j) = h_in(i,j) ; h_R(i,j) = h_in(i,j)
        elseif (12.0*curv*(h_in(i,j) - h_min) < (curv**2 + 3.0*dh**2)) then
          ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
          ! be limited in this case.  0 < scale < 1.
          scale = 12.0*curv*(h_in(i,j) - h_min) / (curv**2 + 3.0*dh**2)
          h_L(i,j) = h_in(i,j) + scale*(h_L(i,j) - h_in(i,j))
          h_R(i,j) = h_in(i,j) + scale*(h_R(i,j) - h_in(i,j))
        endif
      endif
    endif
  enddo ; enddo

end subroutine PPM_limit_pos

!> Limit the left or right edge values of the PPM reconstruction to be monotonic
!! using prescription of Colella and Woodward, 1984.
subroutine PPM_limit_CW84(h_in, h_l, h_r, G, iis, iie, jis, jie)
  type(SIS_hor_grid_type),          intent(in)    :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: h_in !< Initial thickness of a category [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: h_L !< Left edge value of thickness reconstruction [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: h_R !< Right edge value of thickness reconstruction [R Z ~> kg m-2]
  integer,                          intent(in)    :: iis !< The starting i-index to work on
  integer,                          intent(in)    :: iie !< The ending i-index to work on
  integer,                          intent(in)    :: jis !< The starting j-index to work on
  integer,                          intent(in)    :: jie !< The ending j-index to work on
! This subroutine limits the left/right edge values of the PPM reconstruction
! according to the monotonic prescription of Colella and Woodward, 1984.

! Local variables
  real    :: h_i, RLdiff, RLdiff2, RLmean, FunFac
  character(len=256) :: mesg
  integer :: i,j

  do j=jis,jie ; do i=iis,iie
    ! This limiter monotonizes the parabola following
    ! Colella and Woodward, 1984, Eq. 1.10
    h_i = h_in(i,j)
    if ( ( h_r(i,j) - h_i ) * ( h_i - h_l(i,j) ) <= 0. ) then
      h_l(i,j) = h_i ; h_r(i,j) = h_i
    else
      RLdiff = h_r(i,j) - h_l(i,j)            ! Difference of edge values
      RLmean = 0.5 * ( h_r(i,j) + h_l(i,j) )  ! Mean of edge values
      FunFac = 6. * RLdiff * ( h_i - RLmean ) ! Some funny factor
      RLdiff2 = RLdiff * RLdiff               ! Square of difference
      if ( FunFac >  RLdiff2 ) h_l(i,j) = 3. * h_i - 2. * h_r(i,j)
      if ( FunFac < -RLdiff2 ) h_r(i,j) = 3. * h_i - 2. * h_l(i,j)
    endif
  enddo ; enddo

  return
end subroutine PPM_limit_CW84

!> Initializes the sea ice continuity module
subroutine SIS_continuity_init(Time, G, US, param_file, diag, CS, CS_cvr)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model time.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid type
  type(unit_scale_type),       intent(in)    :: US   !< A structure with unit conversion factors
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(SIS_continuity_CS),     pointer       :: CS   !< The control structure for mass transport that
                                                     !! is carried out by this module; it is allocated
                                                     !! and populated here.
  type(SIS_continuity_CS),     optional, pointer :: CS_cvr !< A secondary control structure for the
                                                     !! transport of ice cover.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40) :: mdl = "SIS_continuity" ! This module's name.
  character(len=40) :: mesg    ! Message for error messages.

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_continuity_init called with associated control structure.")
    return
  endif
  allocate(CS)

! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "SIS_CONTINUITY_SCHEME", mesg, &
          desc="The horizontal transport scheme used in continuity:\n"//&
          "  UPWIND_2D - Non-directionally split upwind\n"//&
          "  PCM       - Directionally split piecewise constant\n"//&
          "  PPM:C2PD  - Positive definite PPM with 2nd order edge values\n"//&
          "  PPM:C2MO  - Monotonic PPM with 2nd order edge values\n", &
          default='UPWIND_2D')
  CS%use_upwind2d = .false. ; CS%upwind_1st = .false. ; CS%simple_2nd = .false.
  CS%monotonic = .false.
  select case (trim(mesg))
    case ("UPWIND_2D") ; CS%use_upwind2d = .true.
    case ("PCM")      ; CS%upwind_1st = .true.
    case ("PPM:C2PD") ; CS%simple_2nd = .true.
    case ("PPM:C2MO") ; CS%simple_2nd = .true. ; CS%monotonic = .true.
    case default ; call SIS_error(FATAL, "SIS_continuity, SIS_continuity_init: "//&
           "Unknown SIS_CONTINUITY_SCHEME = "//trim(mesg))
  end select
  call obsolete_logical(param_file, "MONOTONIC_CONTINUITY", &
       hint="Use SIS_CONTINUITY_SCHEME instead.")
  call obsolete_logical(param_file, "UPWIND_2D_CONTINUITY", &
       hint="Use SIS_CONTINUITY_SCHEME instead.")
  call obsolete_logical(param_file, "SIMPLE_2ND_PPM_CONTINUITY", &
       hint="Use SIS_CONTINUITY_SCHEME instead.")
  call obsolete_logical(param_file, "UPWIND_1ST_CONTINUITY", &
       hint="Use SIS_CONTINUITY_SCHEME instead.")

  call get_param(param_file, mdl, "CONT_PPM_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the "//&
                 "tracer cell areas when estimating CFL numbers.", &
                 default=.false.)
  call get_param(param_file, mdl, "CONTINUITY_H_NEGLECT", CS%h_neglect_cont, &
                 "The category ice mass per ocean cell area below which the transport "//&
                 "within this thickness category of out of a cell is set to zero.  A suggested "//&
                 "non-default value might be of order 3e-32 kg m-2, which is one molecule of "//&
                 "ice per square kilometer.", &
                 default=0.0, units="kg m-2", scale=US%kg_m3_to_R*US%m_to_Z)
  call get_param(param_file, mdl, "CONTINUITY_FRAC_NEGLECT", CS%frac_neglect, &
                 "When the total fluxes are distributed between categories with "//&
                 "MERGED_CONTINUITY, any category whose ice is less than this fraction of the "//&
                 "total mass contributes no flux.  Without MERGED_CONTINUITY, any snow or "//&
                 "melt pond transport that is less than this fraction of the ice transport "//&
                 "is zeroed out.  A suggested non-default value might be of order 1e-20.", &
                 default=0.0, units="nondim")

  CS%diag => diag

  if (present(CS_cvr)) then
    allocate(CS_cvr)

    call get_param(param_file, mdl, "SIS_COVER_TRANSPORT_SCHEME", mesg, &
            desc="The horizontal transport scheme used for projections of ice cover:\n"//&
            "  UPWIND_2D - Non-directionally split upwind\n"//&
            "  PCM       - Directionally split piecewise constant\n"//&
            "  PPM:C2PD  - Positive definite PPM with 2nd order edge values\n"//&
            "  PPM:C2MO  - Monotonic PPM with 2nd order edge values\n", &
            default='UPWIND_2D')
    CS_cvr%use_upwind2d = .false. ; CS_cvr%upwind_1st = .false. ; CS_cvr%simple_2nd = .false.
    CS_cvr%monotonic = .false.
    select case (trim(mesg))
      case ("UPWIND_2D") ; CS_cvr%use_upwind2d = .true.
      case ("PCM")      ; CS_cvr%upwind_1st = .true.
      case ("PPM:C2PD") ; CS_cvr%simple_2nd = .true.
      case ("PPM:C2MO") ; CS_cvr%simple_2nd = .true. ; CS_cvr%monotonic = .true.
      case default ; call SIS_error(FATAL, "SIS_continuity, SIS_continuity_init: "//&
             "Unknown SIS_COVER_TRANSPORT_SCHEME = "//trim(mesg))
    end select

    call get_param(param_file, mdl, "COVER_PPM_VOLUME_BASED_CFL", CS_cvr%vol_CFL, &
                   "If true, use the ratio of the open face lengths to the cell "//&
                   "areas when estimating CFL numbers in ice cover transport.", &
                   default=.false.)
  endif

  id_clock_update = cpu_clock_id('(Ocean continuity update)', grain=CLOCK_ROUTINE)
  id_clock_correct = cpu_clock_id('(Ocean continuity correction)', grain=CLOCK_ROUTINE)

end subroutine SIS_continuity_init

!> Deallocate memory associated with the control structure for the SIS_continuity module
subroutine SIS_continuity_end(CS)
  type(SIS_continuity_CS), pointer :: CS  !< The control structure returned by a previous call
                                          !! to SIS_continuity_init that is deallocated here
  deallocate(CS)
end subroutine SIS_continuity_end

end module SIS_continuity
