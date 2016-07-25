!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!                       FIVE-LAYER VERTICAL THERMODYNAMICS                     !
!                                                                              !
! Reference:                                                                   !
!   Winton, M., 2011:  A conservative non-iterative n-layer sea ice            !
!   temperature solver, in prep.                                               !
!                                                                              !
!                                                                              !
!         ->+---------+ <- ts - diagnostic surface temperature ( <= 0C )       !
!        /  |         |                                                        !
!      hs   |  snow   | <- tsn   One snow layer with heat capacity
!        \  |         |                                                        !
!         =>+---------+                                                        !
!        /  |         |                                                        !
!       /   |         | <- t1    Four salty ice layers with heat capacity      !
!      /    |         |                                                        !
!     /     |         | <- t2                                                  !
!   hi      |...ice...|                                                        !
!     \     |         | <- t3                                                  !
!      \    |         |                                                        !
!       \   |         | <- t4                                                  !
!        \  |         |                                                        !
!         ->+---------+ <- base of ice fixed at seawater freezing temp.        !
!                                                                              !
!                                         Mike Winton (Michael.Winton@noaa.gov)!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

! TK mod:  SLAB_ICE treatment modified to follow supersource
!          (after Bryan 1969).  The conductive heat
!           flux from ice to atmosphere is computed based on
!           an effective ice thickness which ensures a minimum
!           thickness of 1.7cm for the calculation.

module ice_thm_mod

implicit none ; private

public :: slab_ice_optics, get_thermo_coefs

! salinities from S(z) = 0.5*3.2*(1.0-cos(3.1416*z**(0.407/(z+0.573))))
! z=[1 3 5 7]/8 ; ref: Hunke et al: CICE V. 4.0, 2008, p. 26
real, parameter :: SI1   = 0.65      ! salinity of sea ice top layer
real, parameter :: SI2   = 2.35      ! salinity of sea ice second layer
real, parameter :: SI3   = 3.03      ! salinity of sea ice third layer
real, parameter :: SI4   = 3.19      ! salinity of sea ice bottom layer

!
! slab ice optics specific parameters
!
real, parameter :: CRIT_THICKNESS       = 1.00
real, parameter :: T_RANGE              = 10.0
real, parameter :: MIN_ICE_ALB          = 0.55   ! coupled model uses 0.55
real, parameter :: MAX_ICE_ALB          = 0.80
real, parameter :: ALB_OCEAN            = 0.10

contains

subroutine slab_ice_optics(hs, hi, ts, tfw, albedo)
  real, intent(in   ) :: hs  ! snow thickness (m-snow)
  real, intent(in   ) :: hi  ! ice thickness (m-ice)
  real, intent(in   ) :: ts  ! surface temperature
  real, intent(in   ) :: tfw ! seawater freezing temperature
  real, intent(  out) :: albedo ! ice surface albedo (0-1)
  real :: alb, as, ai, cs
  real :: thick_ice_alb, tcrit, fh

  tcrit = tfw - T_RANGE
  if (ts <= tcrit) then
    thick_ice_alb = MAX_ICE_ALB
  else if (ts >= tfw) then
    thick_ice_alb = MIN_ICE_ALB
  else
    thick_ice_alb = MAX_ICE_ALB + (MIN_ICE_ALB-MAX_ICE_ALB)*(ts-tcrit)/T_RANGE
  endif

  if (hi >= crit_thickness) then
    albedo = thick_ice_alb
  else
    albedo = ALB_OCEAN + (thick_ice_alb-ALB_OCEAN)*sqrt(hi/CRIT_THICKNESS)
  endif

end subroutine slab_ice_optics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_thermo_coefs - return various thermodynamic coefficients.                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine get_thermo_coefs(ice_salinity)
  real, dimension(:), optional, intent(out) :: ice_salinity
! Arguments: ice_salinity - The specified salinity of each layer when the
!                           thermodynamic salinities are pre-specified.
  integer k, nk

  if (present(ice_salinity)) then
    nk = size(ice_salinity)
    if (nk >= 1) ice_salinity(1) = SI1
    if (nk >= 2) ice_salinity(2) = SI2
    if (nk >= 3) ice_salinity(3) = SI3
    if (nk >= 4) ice_salinity(4) = SI4
    do k=5,nk ;  ice_salinity(k) = SI4 ; enddo
  endif

end subroutine get_thermo_coefs

end module ice_thm_mod
