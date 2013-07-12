module ice_dyn_mod
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
!                                                                              !
! SEA ICE DYNAMICS using ELASTIC-VISCOUS-PLASTIC RHEOLOGY adapted from         !
! Hunke and Dukowicz (JPO 1997, H&D hereafter) -Mike Winton (Michael.Winton@noaa.gov) !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser, only : get_param, read_param, log_version, param_file_type
  use diag_manager_mod, only:  send_data

  use mpp_domains_mod, only: mpp_update_domains, BGRID_NE
  use constants_mod,   only: pi
  use ice_grid_mod,    only: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed
  use ice_grid_mod,    only: dt_evp, evp_sub_steps
use ice_grid_mod,    only: sea_ice_grid_type
!  use ice_thm_mod,     only: DI, DS, DW

implicit none ; private

public :: ice_dyn_init, ice_dynamics

type, public :: ice_dyn_CS ; private
  ! parameters for calculating water drag and internal ice stresses
  logical :: SLAB_ICE = .false. ! should we do old style GFDL slab ice?
  real :: p0 = 2.75e4         ! pressure constant (Pa)
  real :: c0 = 20.0           ! another pressure constant
  real :: cdw = 3.24e-3       ! ice/water drag coef.
  real :: blturn = 25.0       ! air/water surf. turning angle (NH) 25
  real :: EC = 2.0            ! yield curve axis ratio
  real :: MIV_MIN =  1.0      ! min ice mass to do dynamics (kg/m^2)
  real :: Rho_ocean = 1030.0  ! The nominal density of sea water, in kg m-3.
  real :: Rho_ice = 905.0     ! The nominal density of sea ice, in kg m-3.
  real :: Rho_snow = 330.0    ! The nominal density of snow on sea ice, in
                              ! kg m-3.
  type(time_type), pointer :: Time ! A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_fix = -1, id_fiy = -1, id_fcx = -1, id_fcy = -1                  
  integer :: id_fwx = -1, id_fwy = -1, id_sigi = -1, id_sigii = -1                  
  integer :: id_stren = -1, id_ui = -1, id_vi = -1
end type ice_dyn_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_dyn_param - set ice dynamic parameters                                   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_dyn_init(Time, G, param_file, diag, CS, p0_in, c0_in, cdw_in, wd_turn_in, slab_ice_in)
  type(time_type), target, intent(in)    :: Time
  type(sea_ice_grid_type), intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(ice_dyn_CS),        pointer    :: CS
    real,    intent(in)   :: p0_in, c0_in, cdw_in, wd_turn_in
    logical, intent(in)   :: slab_ice_in

! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (inout)   AD - A structure pointing to the various accelerations in
!                 the momentum equations.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module

!   This subroutine sets the parameters and registers the diagnostics associated
! with the ice dynamics.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "ice_dyn" ! This module's name.

    real, parameter       :: missing = -1e34

  if (associated(CS)) then
    call SIS_error(WARNING, "ice_dyn_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag
  CS%Time => Time

  CS%p0 = p0_in
  CS%c0 = c0_in
  CS%cdw = cdw_in
  CS%blturn = wd_turn_in
  CS%SLAB_ICE = slab_ice_in
  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "ICE_STRENGTH_PSTAR", CS%p0, &
                 "A constant in the expression for the ice strength, \n"//&
                 "P* in Hunke & Dukowicz 1997.", units="Pa", default=2.75e4)
  call get_param(param_file, mod, "ICE_STRENGTH_CSTAR", CS%c0, &
                 "A constant in the exponent of the expression for the \n"//&
                 "ice strength, c* in Hunke & Dukowicz 1997.", &
                 units="nondim", default=20.)
  call get_param(param_file, mod, "ICE_CDRAG_WATER", CS%cdw, &
                 "The drag coefficient between the sea ice and water.", &
                 units="nondim", default=3.24e-3)

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
  call get_param(param_file, mod, "AIR_WATER_STRESS_TURN_ANGLE", CS%blturn, &
                 "An angle by which to rotate the velocities at the air- \n"//&
                 "water boundary in calculating stresses.", units="degrees", &
                 default=0.0)


  CS%id_sigi  = register_diag_field('ice_model','SIGI' ,G%axesT1, Time,         &
            'first stress invariant', 'none', missing_value=missing)
  CS%id_sigii = register_diag_field('ice_model','SIGII' ,G%axesT1, Time,        &
            'second stress invariant', 'none', missing_value=missing)
  CS%id_stren = register_diag_field('ice_model','STRENGTH' ,G%axesT1, Time,     &
            'ice strength', 'Pa*m', missing_value=missing)
  CS%id_fix   = register_diag_field('ice_model', 'FI_X', G%axesB1, Time,        &
            'ice internal stress - x component', 'Pa', missing_value=missing)
  CS%id_fiy   = register_diag_field('ice_model', 'FI_Y', G%axesB1, Time,        &
            'ice internal stress - y component', 'Pa', missing_value=missing)
  CS%id_fcx   = register_diag_field('ice_model', 'FC_X', G%axesB1, Time,        &
            'coriolis force - x component', 'Pa', missing_value=missing)
  CS%id_fcy   = register_diag_field('ice_model', 'FC_Y', G%axesB1, Time,        &
            'coriolis force - y component', 'Pa', missing_value=missing)
  CS%id_fwx   = register_diag_field('ice_model', 'FW_X', G%axesB1, Time,        &
            'water stress on ice - x component', 'Pa', missing_value=missing)
  CS%id_fwy   = register_diag_field('ice_model', 'FW_Y', G%axesB1, Time,        &
            'water stress on ice - y component', 'Pa', missing_value=missing)
  CS%id_ui    = register_diag_field('ice_model', 'UI', G%axesB1, Time,          &
            'ice velocity - x component', 'm/s', missing_value=missing)
  CS%id_vi    = register_diag_field('ice_model', 'VI', G%axesB1, Time,          &
            'ice velocity - y component', 'm/s', missing_value=missing)

end subroutine ice_dyn_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! find_ice_strength - magnitude of force on ice in plastic deformation         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine find_ice_strength(hi, ci, ice_strength, G, CS) ! ??? may change to do loop
  real, dimension(isd:ied,jsd:jed), intent(in)  :: hi, ci
  real, dimension(isd:ied,jsd:jed), intent(out) :: ice_strength
  type(sea_ice_grid_type),          intent(in)  :: G
  type(ice_dyn_CS),                 pointer     :: CS

  integer :: i, j

  do j=jsc,jec ; do i=isc,iec
    ice_strength(i,j) = CS%p0*hi(i,j)*ci(i,j)*exp(-CS%c0*(1-ci(i,j)))
  enddo ; enddo

end subroutine find_ice_strength

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_dynamics - take a single dynamics timestep with EVP subcycles            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_dynamics(ci, hs, hi, ui, vi, sig11, sig22, sig12, uo, vo,       &
     fxat, fyat, sea_lev, fxoc, fyoc, G, CS)

!!$    real, intent(in   ), dimension(isd:ied,jsd:jed) :: ci, hs, hi  ! ice properties
    real, intent(in   ), dimension(isd:ied,jsd:jed) :: ci, hs, hi  ! ice properties
    real, intent(inout), dimension(isd:ied,jsd:jed) :: ui, vi      ! ice velocity
    real, intent(inout), dimension(isd:ied,jsd:jed) :: sig11, sig22, sig12       ! stress tensor
    real, intent(in   ), dimension(isd:ied,jsd:jed) :: uo, vo      ! ocean velocity
!!$    real, intent(in   ), dimension(isc:iec,jsc:jec) :: fxat, fyat  ! air stress on ice
    real, intent(in   ), dimension(isc:,jsc:) :: fxat, fyat  ! air stress on ice
    real, intent(in   ), dimension(isd:ied,jsd:jed) :: sea_lev     ! sea level
    real, intent(  out), dimension(isc:iec,jsc:jec) :: fxoc, fyoc  ! ice stress on ocean
  type(sea_ice_grid_type), intent(in) :: G
  type(ice_dyn_CS),        pointer    :: CS

  real, dimension(isc:iec,jsc:jec) :: fxic, fyic  ! ice int. stress
  real, dimension(isc:iec,jsc:jec) :: fxco, fyco  ! coriolis force

  real, dimension(isd:ied,jsd:jed)    :: prs                    ! ice pressure
  real                                :: zeta, eta              ! bulk/shear viscosities
  real, dimension(isc:iec,jsc:jec)    :: strn11, strn12, strn22 ! strain tensor

    real,    dimension(isd:ied,jsd:jed) :: mit                 ! mass on t-points
    real,    dimension(isd:ied,jsd:jed) :: miv                 ! mass on v-points
    real,    dimension(isd:ied,jsd:jed) :: civ                 ! conc. on v-points
    real,    dimension(isd:ied,jsd:jed) :: diag_val            ! A temporary diagnostic array
    complex                             :: rr                  ! linear drag coefficient
    real                                :: fxic_now, fyic_now  ! ice internal stress

    ! temporaries for strain calculation
    real, dimension(isd:ied,jsd:jed) :: &
      grid_fac1, grid_fac2, grid_fac3, grid_fac4
    
    ! temporaries for ice stress calculation
    real                             :: del2, a, b, tmp
    real, dimension(isc:iec,jsc:jec) :: edt, mp4z, t0, t1, t2
    real                             :: f11, f22
    real, dimension(isd:ied,jsd:jed) :: sldx, sldy
    real, dimension(isd:ied,jsd:jed) :: dydx, dxdy
    real   :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7

    ! for velocity calculation
    real,    dimension(isc:iec,jsc:jec) :: dtmiv, rpart, fpart, uvfac
  real :: EC2I  ! 1/EC^2, where EC is the yield curve axis ratio.
    complex                             :: newuv

  logical :: sent
  integer :: i,j,l

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "ice_dynamics: Module must be initialized before it is used.")

  EC2I = 1.0/(CS%EC*CS%EC)

  fxoc(:,:) = 0.0 ; fyoc(:,:) = 0.0 ! zero these for summing later
  fxic(:,:) = 0.0 ; fyic(:,:) = 0.0
  fxco(:,:) = 0.0 ; fyco(:,:) = 0.0

  if (CS%SLAB_ICE) then
    ui(:,:) = uo(:,:) ; vi(:,:) = vo(:,:)
    fxoc(:,:) = fxat(:,:) ; fyoc(:,:) = fyat(:,:)
    return
  end if

  if (evp_sub_steps==0) return;

  !### ADD PARENTHESIS FOR REPRODUCIBILITY.
  do J=jsc-1,jec ; do I=isc-1,iec
    dydx(I,J) = 0.5*(G%dyT(i+1,j+1) - G%dyT(i,j+1) + G%dyT(i+1,j) - G%dyT(i,j) )
    dxdy(I,J) = 0.5*(G%dxT(i+1,j+1) - G%dxT(i+1,j) + G%dxT(i,j+1) - G%dxT(i,j) )
  enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    grid_fac1(i,j) = (G%dxCv(i,J)-G%dxCv(i,J-1))/G%dyT(i,j) !### Use G%IdyT?
    grid_fac2(i,j) = (G%dyCu(I,j)-G%dyCu(I-1,j))/G%dxT(i,j)
    grid_fac3(i,j) = 0.5*G%dyT(i,j)/G%dxT(i,j) !### Use G%IdxT?
    grid_fac4(i,j) = 0.5*G%dxT(i,j)/G%dyT(i,j)
  enddo ; enddo

  ! sea level slope force
  ! ### Add parentheses for rotational consistency.
  ! ### Multiply by G%IdxBu, etc.
  do j=jsc,jec ; do i=isc,iec ! ###RESIZE  do J=jsc-1,jec ; do I=isc-1,iec
    sldx(I,J) = -dt_evp*G%g_Earth*(0.5*(sea_lev(i+1,j+1)-sea_lev(i,j+1) &
              + sea_lev(i+1,j)-sea_lev(i,j))) / G%dxBu(i,J)
    sldy(I,J) = -dt_evp*G%g_Earth*(0.5*(sea_lev(i+1,j+1)-sea_lev(i+1,j) &
              + sea_lev(i,j+1)-sea_lev(i,j))) / G%dyBu(I,J)
  enddo ; enddo

  ! put ice/snow mass and concentration on v-grid, first finding mass on t-grid.
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    mit(i,j) = ci(i,j)*(hi(i,j)*CS%Rho_ice + hs(i,j)*CS%Rho_snow)
  enddo ; enddo
  !### ADD PARENTHESIS FOR REPRODUCIBILITY.
  do J=jsc-1,jec ; do I=isc-1,iec ; if (G%mask2dBu(i,j) > 0.5 ) then
    miv(I,J) = 0.25*( mit(i+1,j+1) + mit(i+1,j) + mit(i,j+1) + mit(i,j) )
    civ(I,J) = 0.25*( ci(i+1,j+1) + ci(i+1,j) + ci(i,j+1) + ci(i,j) )
  else
    miv(I,J) = 0.0 ; civ(I,J) = 0.0
  endif ; enddo ; enddo
    
  ! precompute prs, elastic timestep parameter, and linear drag coefficient
  !
  call find_ice_strength(hi, ci, prs, G, CS)

  do j=jsc,jec ; do i=isc,iec
    if(G%dxT(i,j) < G%dyT(i,j) ) then
      edt(i,j) = (CS%Rho_ice*G%dxT(i,j)*G%dxT(i,j)*ci(i,j)*hi(i,j))/(2*dt_evp)
    else
      edt(i,j) = (CS%Rho_ice*G%dyT(i,j)*G%dyT(i,j)*ci(i,j)*hi(i,j))/(2*dt_evp)
    endif
  enddo ; enddo

  do j=jsc,jec ; do i=isc,iec ! ###RESIZE do J=jsc-1,jec ; do I=isc-1,iec
    if ((G%mask2dBu(I,J)>0.5) .and. (miv(I,J) > CS%MIV_MIN) ) then ! values for velocity calculation (on v-grid)
      dtmiv(i,j) = dt_evp/miv(i,j)
    else
      ui(I,J) = 0.0 ; vi(I,J) = 0.0
    endif
  enddo ; enddo

  do l=1,evp_sub_steps

    ! calculate strain tensor for viscosities and forcing elastic eqn.
    call mpp_update_domains(ui, vi, Domain, gridtype=BGRID_NE)

  !### ADD PARENTHESES FOR CLARITY
  !### MULTIPLY BY GRID METRIC INVERSES FOR EFFICIENCY.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_strain - calculate generalized orthogonal coordinate strain tensor       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    do j=jsc,jec ; do i=isc,iec
      strn11(i,j) = (0.5*(ui(i,j)-ui(i-1,j)+ui(i,j-1)-ui(i-1,j-1))        &
                  + 0.25*(vi(i,j)+vi(i,j-1)+vi(i-1,j)+vi(i-1,j-1))*grid_fac1(i,j)) / &
                  G%dxT(i,j)   !### Use G%IdxT?
      strn22(i,j) = (0.5*(vi(i,j)-vi(i,j-1)+vi(i-1,j)-vi(i-1,j-1))        &
                  + 0.25*(ui(i,j)+ui(i,j-1)+ui(i-1,j)+ui(i-1,j-1))*grid_fac2(i,j)) / &
                  G%dyT(i,j)
      strn12(i,j) = grid_fac3(i,j)*(0.5*(vi(i,j)/G%dyBu(I,J) - vi(i-1,j)/G%dyBu(I-1,J) + &
                                    vi(i,j-1)/G%dyBu(i,j-1) - vi(i-1,j-1)/G%dyBu(i-1,j-1))) &
                  + grid_fac4(i,j)*(0.5*(ui(i,j)/G%dxBu(i,j)-ui(i,j-1)/G%dxBu(i,j-1) + &
                                    ui(i-1,j)/G%dxBu(i-1,j)-ui(i-1,j-1)/G%dxBu(i-1,j-1)))
    enddo ; enddo

    ! calculate viscosities - how often should we do this ?
    if (l>=1) then
      do j=jsc,jec ; do i=isc,iec
        del2 = (strn11(i,j)*strn11(i,j) + strn22(i,j)*strn22(i,j)) * (1+EC2I)     &
              + 4*EC2I*strn12(i,j)*strn12(i,j) + 2*strn11(i,j)*strn22(i,j)*(1-EC2I)  ! H&D eqn 9

        if ( del2 > 4e-18 ) then
          zeta = 0.5*prs(i,j)/sqrt(del2)
        else
          zeta = 2.5e8*prs(i,j)
        endif

        if (zeta<4e8) zeta = 4e8 ! Hibler uses to prevent nonlinear instability

        eta = zeta*EC2I
        !
        ! some helpful temporaries
        !
        if ( hi(i,j) > 0.0 ) then
          mp4z(i,j) = -prs(i,j)/(4*zeta)
          t0(i,j)   = 2*eta/(2*eta+edt(i,j))
          tmp       = 1/(4*eta*zeta)
          a         = 1/edt(i,j) + (zeta+eta)*tmp ! = 1/edt(i,j) + (1+EC2I)/(4*eta)
          b         = (zeta-eta)*tmp              ! = (1-EC2I)/(4*eta)
          t1(i,j)   = b/a
          t2(i,j)   = a - b*b/a
        endif
      enddo ; enddo
    endif

   ! timestep stress tensor (H&D eqn 21)
    do j=jsc,jec ; do i=isc,iec
      if( (G%mask2dT(i,j)>0.5) .and. &
          (ci(i,j)*(CS%Rho_ice*hi(i,j)+CS%Rho_snow*hs(i,j))>CS%MIV_MIN) ) then
        f11   = mp4z(i,j) + sig11(i,j)/edt(i,j) + strn11(i,j)
        f22   = mp4z(i,j) + sig22(i,j)/edt(i,j) + strn22(i,j)
        sig11(i,j) = (t1(i,j)*f22 + f11) / t2(i,j)
        sig22(i,j) = (t1(i,j)*f11 + f22) / t2(i,j)
        sig12(i,j) = t0(i,j) * (sig12(i,j) + edt(i,j)*strn12(i,j))
      else
        sig11(i,j) = 0.0
        sig22(i,j) = 0.0
        sig12(i,j) = 0.0 ! eliminate internal ice forces 
      endif
    enddo ; enddo

    call mpp_update_domains(sig11, Domain, complete=.false.)
    call mpp_update_domains(sig22, Domain, complete=.false.)
    call mpp_update_domains(sig12, Domain, complete=.true.)

    do j=jsc,jec ; do i=isc,iec ! ###RESIZE  do J=jsc-1,jec ; do I=isc-1,iec
      if( (G%mask2dBu(i,j)>0.5).and.(miv(i,j)>CS%MIV_MIN)) then ! timestep ice velocity (H&D eqn 22)
        rr = CS%cdw*CS%Rho_ocean*abs(cmplx(ui(i,j)-uo(i,j),vi(i,j)-vo(i,j))) * &
             exp(sign(CS%blturn*pi/180,G%CoriolisBu(i,j))*(0.0,1.0))
        !
        ! first, timestep explicit parts (ice, wind & ocean part of water stress)
        !
  !### ADD PARENTHESIS FOR REPRODUCIBILITY.
        tmp1 = 0.5*(sig12(i+1,j+1)*G%dxT(i+1,j+1) - sig12(i+1,j)*G%dxT(i+1,j) + &
                    sig12(i,j+1)*G%dxT(i,j+1) - sig12(i,j)*G%dxT(i,j) )
        tmp2 = 0.5*(sig11(i+1,j+1)*G%dyT(i+1,j+1) - sig11(i,j+1)*G%dyT(i,j+1) + &
                    sig11(i+1,j)*G%dyT(i+1,j) - sig11(i,j)*G%dyT(i,j) )
        tmp6 = 0.5*(sig12(i+1,j+1)*G%dyT(i+1,j+1) - sig12(i,j+1)*G%dyT(i,j+1) + &
                    sig12(i+1,j)*G%dyT(i+1,j) - sig12(i,j)*G%dyT(i,j) )
        tmp7 = 0.5*(sig22(i+1,j+1)*G%dxT(i+1,j+1) - sig22(i+1,j)*G%dxT(i+1,j) + &
                    sig22(i,j+1)*G%dxT(i,j+1) - sig22(i,j)*G%dxT(i,j) )
        tmp3 = 0.25*(sig12(i+1,j+1)+sig12(i+1,j)+sig12(i,j+1)+sig12(i,j) )
        tmp4 = 0.25*(sig22(i+1,j+1)+sig22(i+1,j)+sig22(i,j+1)+sig22(i,j) )
        tmp5 = 0.25*(sig11(i+1,j+1)+sig11(i+1,j)+sig11(i,j+1)+sig11(i,j) )

  !### ADD PARENTHESIS FOR REPRODUCIBILITY.
        fxic_now = ( tmp1 + tmp2 + tmp3*dxdy(i,j) - tmp4*dydx(i,j) ) / (G%dxBu(i,j)*G%dyBu(I,J)) 
        fyic_now = ( tmp6 + tmp7 + tmp3*dydx(i,j) - tmp5*dxdy(i,j) ) / (G%dxBu(i,j)*G%dyBu(I,J)) 

        ui(i,j) = ui(i,j) + (fxic_now + civ(i,j)*fxat(i,j) + &
                             real(civ(i,j)*rr*cmplx(uo(i,j),vo(i,j)))) * dtmiv(i,j) + sldx(i,j)
        vi(i,j) = vi(i,j) + (fyic_now + civ(i,j)*fyat(i,j) + &
                             aimag(civ(i,j)*rr*cmplx(uo(i,j),vo(i,j)))) * dtmiv(i,j) + sldy(i,j)
        !
        ! second, timestep implicit parts (Coriolis and ice part of water stress)
        !
        newuv = cmplx(ui(i,j),vi(i,j)) / &
            (1 + dt_evp*(0.0,1.0)*G%CoriolisBu(i,j) + civ(i,j)*rr*dtmiv(i,j))
        ui(i,j) = real(newuv); vi(i,j) = aimag(newuv)
        !
        ! sum for averages
        !
        fxic(i,j) = fxic(i,j) + fxic_now
        fyic(i,j) = fyic(i,j) + fyic_now
        fxoc(i,j) = fxoc(i,j) +  real(civ(i,j)*rr*cmplx(ui(i,j)-uo(i,j), vi(i,j)-vo(i,j)))
        fyoc(i,j) = fyoc(i,j) + aimag(civ(i,j)*rr*cmplx(ui(i,j)-uo(i,j), vi(i,j)-vo(i,j)))
        fxco(i,j) = fxco(i,j) - miv(i,j)*real ((0.0,1.0)*G%CoriolisBu(i,j) * cmplx(ui(i,j),vi(i,j)))
        fyco(i,j) = fyco(i,j) - miv(i,j)*aimag((0.0,1.0)*G%CoriolisBu(i,j) * cmplx(ui(i,j),vi(i,j)))              
      endif
    enddo ; enddo
  enddo ! l=1,evp_sub_steps

  ! make averages
  ! ### Multiply by reciprocal of evp_sub_steps?
  do j=jsc,jec ; do i=isc,iec ! ###RESIZE  do J=jsc-1,jec ; do I=isc-1,iec
    if( (G%mask2dBu(i,j)>0.5) .and. miv(i,j)>CS%MIV_MIN ) then
       fxoc(i,j) = fxoc(i,j)/evp_sub_steps;  fyoc(i,j) = fyoc(i,j)/evp_sub_steps
       fxic(i,j) = fxic(i,j)/evp_sub_steps;  fyic(i,j) = fyic(i,j)/evp_sub_steps
       fxco(i,j) = fxco(i,j)/evp_sub_steps;  fyco(i,j) = fyco(i,j)/evp_sub_steps             
    endif
  enddo ; enddo

  ! Write out diagnostics associated with the ice dynamics.
!  The diagnistics of fxat and fyat are supposed to be taken over all partitions
!  (ocean & ice), whereas fxat and fyat here are only averaged over ice.
!## if (CS%id_fax>0) &
!##      sent = send_data(CS%id_fax, all_avg(CS%flux_u_top_bgrid(isc:iec,jsc:jec,:),CS%part_size_uv), CS%Time)
!## if (CS%id_fay>0) &
!##      sent = send_data(CS%id_fay, all_avg(CS%flux_v_top_bgrid(isc:iec,jsc:jec,:),CS%part_size_uv), CS%Time)

  if (CS%id_fix>0) sent = send_data(CS%id_fix, fxic, CS%Time)
  if (CS%id_fiy>0) sent = send_data(CS%id_fiy, fyic, CS%Time)
  if (CS%id_fcx>0) sent = send_data(CS%id_fcx, fxco, CS%Time)
  if (CS%id_fcy>0) sent = send_data(CS%id_fcy, fyco, CS%Time)
  if (CS%id_fwx>0) sent = send_data(CS%id_fwx, -fxoc, CS%Time) ! water force on ice
  if (CS%id_fwy>0) sent = send_data(CS%id_fwy, -fyoc, CS%Time) ! ...= -ice on water

  if (CS%id_sigi>0) then
    diag_val(:,:) =  sigI(hi, ci, sig11, sig22, sig12, G, CS)
    sent = send_data(CS%id_sigi, diag_val(isc:iec,jsc:jec), CS%Time) !### , mask=CS%mask)
  endif
  if (CS%id_sigii>0) then
    diag_val(:,:) = sigII(hi, ci, sig11, sig22, sig12, G, CS)
    sent = send_data(CS%id_sigii, diag_val(isc:iec,jsc:jec), CS%Time) !### , mask=Ice%mask)
  endif
  if (CS%id_stren>0) then
    call find_ice_strength(hi, ci, diag_val, G, CS)
    sent = send_data(CS%id_stren, diag_val(isc:iec,jsc:jec), CS%Time) !, mask=Ice%mask)
  endif

  if (CS%id_ui>0) sent = send_data(CS%id_ui, ui(isc:iec,jsc:jec), CS%Time)
  if (CS%id_vi>0) sent = send_data(CS%id_vi, vi(isc:iec,jsc:jec), CS%Time)

end subroutine ice_dynamics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! sigI - first stress invariant                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function sigI(hi, ci, sig11, sig22, sig12, G, CS)
  real, dimension(isd:ied,jsd:jed), intent(in) :: hi, ci, sig11, sig22, sig12
  real, dimension(isd:ied,jsd:jed)             :: sigI
  type(sea_ice_grid_type), intent(in)    :: G
  type(ice_dyn_CS),        pointer       :: CS

  integer :: i, j

  call find_ice_strength(hi, ci, sigI, G, CS)

  do j=jsc,jec ; do i=isc,iec
    if (sigI(i,j) > 0.0) sigI(i,j) = (sig11(i,j) + sig22(i,j)) / sigI(i,j)
  enddo ; enddo

end function sigI

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! sigII - second stress invariant                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function sigII(hi, ci, sig11, sig22, sig12, G, CS)
  real, dimension(isd:ied,jsd:jed), intent(in) :: hi, ci, sig11, sig22, sig12
  real, dimension(isd:ied,jsd:jed)             :: sigII
  type(sea_ice_grid_type), intent(in)    :: G
  type(ice_dyn_CS),        pointer       :: CS

  integer :: i, j

  call find_ice_strength(hi, ci, sigII, G, CS)

  do j=jsc,jec ; do i=isc,iec
    if (sigII(i,j) > 0.0) sigII(i,j) = (((sig11(i,j)-sig22(i,j))**2+4*sig12(i,j)*sig12(i,j))/(sigII(i,j)**2))**0.5
  enddo ; enddo

end function sigII

end module ice_dyn_mod
