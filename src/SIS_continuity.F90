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

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_obsolete_params, only : obsolete_logical
use SIS_diag_mediator, only : time_type, SIS_diag_ctrl
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type
! use MOM_variables, only : ocean_OBC_type, OBC_SIMPLE
! use MOM_variables, only : OBC_FLATHER_E, OBC_FLATHER_W, OBC_FLATHER_N, OBC_FLATHER_S

implicit none ; private

#include <SIS2_memory.h>

public ice_continuity, SIS_continuity_init, SIS_continuity_end

integer :: id_clock_update  !< A CPU time clock ID
integer :: id_clock_correct !< A CPU time clock ID

!> The control structure with parameters regulating the continuity solver
type, public :: SIS_continuity_CS ; private
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.
  logical :: use_upwind2d    !< If true, use the non-split upwind scheme that was
                             !! used in older versions of SIS.
  logical :: upwind_1st      !< If true, use a directionally-split first-order
                             !! upwind scheme.
  logical :: monotonic       !< If true, use the Colella & Woodward monotonic
                             !! limiter; otherwise use a simple positive
                             !! definite limiter.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic
                             !! mean) interpolation of the edge values instead
                             !! of the higher order interpolation.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths
                             !! to the tracer cell areas when estimating CFL
                             !! numbers.
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
subroutine ice_continuity(u, v, hin, h, uh, vh, dt, G, IG, CS)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: u   !< Zonal ice velocity, in m s-1.
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: v   !< Meridional ice velocity, in m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(in)    :: hin !< Initial ice or snow thickness by category, in H.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(inout) :: h   !< Final ice or snow thickness by category, in H.
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(out)   :: uh  !< Volume flux through zonal faces = u*h*dy, H m2 s-1.
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                           intent(out)   :: vh  !< Volume flux through meridional faces = v*h*dx, in H m2 s-1.
  real,                    intent(in)    :: dt  !< Time increment in s.
  type(SIS_continuity_CS), pointer       :: CS  !< The control structure returned by a
                                                !! previous call to SIS_continuity_init.
!    This subroutine time steps the category thicknesses, using a monotonically
!  limit, directionally split PPM scheme, based on Lin (1994).  In the following
!  documentation, H is used for the units of thickness (usually m or kg m-2.)

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)) :: &
    h_input      ! Left and right face thicknesses, in H.
  type(loop_bounds_type) :: LB  ! A structure with the active loop bounds.
  real    :: h_up
  integer :: is, ie, js, je, nCat, stensil
  integer :: i, j, k

  logical :: x_first
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nCat = IG%CatIce

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_continuity: Module must be initialized before it is used.")
  x_first = (MOD(G%first_direction,2) == 0)

  stensil = 3 ; if (CS%simple_2nd) stensil = 2 ; if (CS%upwind_1st) stensil = 1

  do k=1,nCat ; do j=js,je ; do i=is,ie ; if (h(i,j,k) < 0.0) then
    call SIS_error(FATAL, 'Negative thickness input to ice_continuity().')
  endif ; enddo ; enddo ; enddo

  if (CS%use_upwind2d) then
    ! This reproduces the scheme that was originally used in SIS1.
!$OMP parallel default(none) shared(G,is,ie,js,je,u,v,hin,uh,vh,h,dt,nCat) &
!$OMP                       private(h_up)
!$OMP do
    do j=js,je ; do k=1,nCat ; do I=is-1,ie
      if (u(I,j) >= 0.0) then ; h_up = hin(i,j,k)
      else ; h_up = hin(i+1,j,k) ; endif
      uh(I,j,k) = G%dy_Cu(I,j) * u(I,j) * h_up
    enddo ; enddo ; enddo
!$OMP do
    do J=js-1,je ; do k=1,nCat ; do i=is,ie
      if (v(i,J) >= 0.0) then ; h_up = hin(i,j,k)
      else ; h_up = hin(i,j+1,k) ; endif
      vh(i,J,k) = G%dx_Cv(i,J) * v(i,J) * h_up
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
    call zonal_mass_flux(u, hin, uh, dt, G, IG, CS, LB)

    call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) shared(LB,nCat,G,uh,hin,dt,h)
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
    call meridional_mass_flux(v, h, vh, dt, G, IG, CS, LB)

    call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) shared(nCat,LB,h,dt,G,vh)
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

    call meridional_mass_flux(v, hin, vh, dt, G, IG, CS, LB)

    call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) shared(nCat,LB,h,hin,dt,G,vh)
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
    call zonal_mass_flux(u, h, uh, dt, G, IG, CS, LB)

    call cpu_clock_begin(id_clock_update)
!$OMP parallel do default(none) shared(nCat,LB,h,dt,G,uh)
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

!> Calculates the mass or volume fluxes through the zonal
!! faces, and other related quantities.
subroutine zonal_mass_flux(u, h_in, uh, dt, G, IG, CS, LB)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: u   !< Zonal ice velocity, in m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(in)    :: h_in !< Category thickness used to calculate the fluxes, in H.
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(out)   :: uh  !< Volume flux through zonal faces = u*h*dy, H m2 s-1.
  real,                    intent(in)    :: dt  !< Time increment in s
  type(SIS_continuity_CS), pointer       :: CS  !< The control structure returned by a
                                                !! previous call to SIS_continuity_init.
  type(loop_bounds_type),  intent(in)    :: LB  !< A structure with the active loop bounds.
!   This subroutine calculates the mass or volume fluxes through the zonal
! faces, and other related quantities.

  ! Local variables
!  real, dimension(SZIB_(G)) :: &
!    duhdu      ! Partial derivative of uh with u, in H m.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &    ! The total thickness summed across categories, in H.
    I_htot, &  ! The inverse of htot or 0, in H-1.
    hl, hr      ! Left and right face thicknesses, in H.
  real, dimension(SZIB_(G)) :: &
    uhtot      ! The total transports in H m2 s-1.
  real :: CFL  ! The CFL number based on the local velocity and grid spacing, ND.
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
!  real :: h_marg ! The marginal thickness of a flux, in H.
  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west, in m.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = IG%CatIce

  call cpu_clock_begin(id_clock_update)

  htot(:,:) = 0.0
!$OMP parallel do default(none) shared(jsh,jeh,nz,G,htot,h_in,I_htot)
  do j=jsh,jeh
    do k=1,nz ; do i=G%isd,G%ied
      htot(i,j) = htot(i,j) + h_in(i,j,k)
    enddo ; enddo
    do i=G%isd,G%ied
      I_htot(i,j) = 0.0 ; if (htot(i,j) > 0.0) I_htot(i,j) = 1.0 / htot(i,j)
    enddo
  enddo

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
    do k=1,nz ; do I=ish-1,ieh
      if (u(I,j) >= 0.0) then
        uh(I,j,k) = uhtot(I) * (h_in(i,j,k) * I_htot(i,j))
      else
        uh(I,j,k) = uhtot(I) * (h_in(i+1,j,k) * I_htot(i+1,j))
      endif
    enddo ; enddo

  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

end subroutine zonal_mass_flux

!> Calculates the mass or volume fluxes through the meridional
!! faces, and other related quantities.
subroutine meridional_mass_flux(v, h_in, vh, dt, G, IG, CS, LB)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: v   !< Meridional ice velocity, in m s-1.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(in)    :: h_in !< Category thickness used to calculate the fluxes, in H.
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)), &
                           intent(out)   :: vh  !< Volume flux through meridional faces = u*h*dy, H m2 s-1.
  real,                    intent(in)    :: dt  !< Time increment in s
  type(SIS_continuity_CS), pointer       :: CS  !< The control structure returned by a
                                                !! previous call to SIS_continuity_init.
  type(loop_bounds_type),  intent(in)    :: LB  !< A structure with the active loop bounds.

!   This subroutine calculates the mass or volume fluxes through the meridional
! faces, and other related quantities.

  ! Local variables
  real, dimension(SZI_(G)) :: &
    dvhdv      ! Partial derivative of vh with v, in m2.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    htot, &    ! The total thickness summed across categories, in H.
    I_htot, &  ! The inverse of htot or 0, in H-1.
    hl, hr     ! Left and right face thicknesses, in m.
  real, dimension(SZI_(G)) :: &
    vhtot      ! The total transports in H m2 s-1.
  real :: CFL ! The CFL number based on the local velocity and grid spacing, ND.
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h_in.
  real :: h_marg ! The marginal thickness of a flux, in m.
  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south, in m.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh ; nz = IG%CatIce

  call cpu_clock_begin(id_clock_update)

  htot(:,:) = 0.0

!$OMP parallel do default(none) shared(ish,ieh,G,nz,htot,h_in,I_htot)
  do j=G%jsd,G%jed
    do k=1,nz ; do i=ish,ieh
      htot(i,j) = htot(i,j) + h_in(i,j,k)
    enddo ; enddo
    do i=ish,ieh
      I_htot(i,j) = 0.0 ; if (htot(i,j) > 0.0) I_htot(i,j) = 1.0 / htot(i,j)
    enddo
  enddo

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
        else ; CFL = -v(i,J) * dt * G%IdyT(i,j+1) ; endif
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
    do k=1,nz ; do i=ish,ieh
      if (v(i,J) >= 0.0) then
        vh(i,J,k) = vhtot(i) * (h_in(i,j,k) * I_htot(i,j))
      else
        vh(i,J,k) = vhtot(i) * (h_in(i,j+1,k) * I_htot(i,j+1))
      endif
    enddo ; enddo

  enddo ! j-loop
  call cpu_clock_end(id_clock_correct)

end subroutine meridional_mass_flux

!> Calculate a piecewise parabolic thickness reconstruction in the x-direction.
subroutine PPM_reconstruction_x(h_in, h_l, h_r, G, LB, h_min, monotonic, simple_2nd)
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Initial thickness of a category, in H
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l !< Left edge value of thickness reconstruction, in H
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r !< Right edge value of thickness reconstruction, in H
  type(loop_bounds_type),           intent(in)  :: LB  !< A structure with the active loop bounds.
  real,                             intent(in)  :: h_min !< The minimum thickness that can be
                                                       !! obtained by a concave parabolic fit, in H.
  logical, optional,                intent(in)  :: monotonic !< If true, use the Colella & Woodward monotonic limiter.
                                                       !! Otherwise use a simple positive-definite limiter.
  logical, optional,                intent(in)  :: simple_2nd !< If true, use the arithmetic mean thicknesses as the
                                                       !! default edge values for a simple 2nd order scheme.
! This subroutine calculates left/right edge valus for PPM reconstruction.

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
!$OMP parallel do default(none) shared(isl,iel,jsl,jel,G,h_in,h_l,h_r) &
!$OMP                          private(h_im1,h_ip1)
    do j=jsl,jel ; do i=isl,iel
      h_im1 = G%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-G%mask2dT(i-1,j)) * h_in(i,j)
      h_ip1 = G%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-G%mask2dT(i+1,j)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_im1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) )
    enddo ; enddo
  else
!$OMP parallel do default(none) shared(isl,iel,jsl,jel,G,h_in,h_l,h_r,slp) &
!$OMP                          private(dMx,dMn,h_im1,h_ip1)
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
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: h_in !< Initial thickness of a category, in H
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_l !< Left edge value of thickness reconstruction, in H
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: h_r !< Right edge value of thickness reconstruction, in H
  type(loop_bounds_type),           intent(in)  :: LB  !< A structure with the active loop bounds.
  real,                             intent(in)  :: h_min !< The minimum thickness that can be
                                                       !! obtained by a concave parabolic fit, in H.
  logical, optional,                intent(in)  :: monotonic !< If true, use the Colella & Woodward monotonic limiter.
                                                       !! Otherwise use a simple positive-definite limiter.
  logical, optional,                intent(in)  :: simple_2nd !< If true, use the arithmetic mean thicknesses as the
                                                       !! default edge values for a simple 2nd order scheme.
! This subroutine calculates left/right edge valus for PPM reconstruction.

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
!$OMP parallel do default(none) shared(isl,iel,jsl,jel,G,h_in,h_l,h_r) &
!$OMP                          private(h_jm1,h_jp1)
    do j=jsl,jel ; do i=isl,iel
      h_jm1 = G%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-G%mask2dT(i,j-1)) * h_in(i,j)
      h_jp1 = G%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-G%mask2dT(i,j+1)) * h_in(i,j)
      h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) )
      h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) )
    enddo ; enddo
  else
!$OMP parallel do default(none) shared(isl,iel,jsl,jel,G,h_in,slp) &
!$OMP                          private(dMx,dMn)
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
!$OMP parallel do default(none) shared(isl,iel,jsl,jel,G,h_in,h_l,h_r,slp) &
!$OMP                          private(h_jm1,h_jp1)
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
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: h_in !< Initial thickness of a category, in H
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: h_L !< Left edge value of thickness reconstruction, in H
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: h_R !< Right edge value of thickness reconstruction, in H
  real,                             intent(in)    :: h_min !< The minimum thickness that can be
                                                           !! obtained by a concave parabolic fit, in H.
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
  real, dimension(SZI_(G),SZJ_(G)), intent(in)    :: h_in !< Initial thickness of a category, in H
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: h_L !< Left edge value of thickness reconstruction, in H
  real, dimension(SZI_(G),SZJ_(G)), intent(inout) :: h_R !< Right edge value of thickness reconstruction, in H
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
subroutine SIS_continuity_init(Time, G, param_file, diag, CS)
  type(time_type),     target, intent(in)    :: Time   !< The sea-ice model's clock,
                                                       !! set with the current model time.
  type(SIS_hor_grid_type),     intent(in)    :: G   !< The horizontal grid type
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(SIS_continuity_CS),     pointer       :: CS  !< The control structure for this module that
                                                    !! is allocated and populated here.

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
    case ("UPWIND_2D")
      CS%use_upwind2d = .true.
    case ("PCM")
      CS%upwind_1st = .true.
    case ("PPM:C2PD")
      CS%simple_2nd = .true.
    case ("PPM:C2MO")
      CS%simple_2nd = .true.
      CS%monotonic = .true.
    case default
      call SIS_error(FATAL, "SIS_continuity, SIS_continuity_init: "//&
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
                 "If true, use the ratio of the open face lengths to the \n"//&
                 "tracer cell areas when estimating CFL numbers.", &
                 default=.false.)

  CS%diag => diag

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
