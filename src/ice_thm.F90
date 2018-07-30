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

!> A module for specifying certain sea-ice properties
module ice_thm_mod

implicit none ; private

public :: get_thermo_coefs

!>@{ salinities from S(z) = 0.5*3.2*(1.0-cos(3.1416*z**(0.407/(z+0.573))))
! z=[1 3 5 7]/8 ; ref: Hunke et al: CICE V. 4.0, 2008, p. 26
real, parameter :: SI1   = 0.65      ! salinity of sea ice top layer
real, parameter :: SI2   = 2.35      ! salinity of sea ice second layer
real, parameter :: SI3   = 3.03      ! salinity of sea ice third layer
real, parameter :: SI4   = 3.19      ! salinity of sea ice bottom layer
!!@}

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> get_thermo_coefs returns various thermodynamic coefficients.
subroutine get_thermo_coefs(ice_salinity)
  real, dimension(:), optional, intent(out) :: ice_salinity !< The specified salinity of each layer
                                            !! when the thermodynamic salinities are pre-specified.
  integer k, nk

  if (present(ice_salinity)) then
    nk = size(ice_salinity)
    if (nk >= 1) ice_salinity(1) = SI1
    if (nk >= 2) ice_salinity(2) = SI2
    if (nk >= 3) ice_salinity(3) = SI3
    if (nk >= 4) ice_salinity(4) = SI4
    do k=5, nk ;  ice_salinity(k) = SI4 ; enddo
  endif

end subroutine get_thermo_coefs

end module ice_thm_mod
