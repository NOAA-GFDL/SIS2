!> Update sea-ice dynamics using elastic-viscous-plastic rheology with a B-grid discretization
module SIS_dyn_bgrid

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
! B-grid SEA ICE DYNAMICS using ELASTIC-VISCOUS-PLASTIC RHEOLOGY adapted from  !
! Hunke and Dukowicz (JPO 1997, H&D hereafter) -Mike Winton (Michael.Winton@noaa.gov) !
! Recoded for consistency with MOM6 and symmetry by R. Hallberg, 7/2013.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_domains,       only : pass_var, pass_vector, BGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, log_param, read_param, log_version, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_unit_scaling,  only : unit_scale_type
use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use SIS_debugging,     only : chksum, Bchksum, hchksum, check_redundant_B
use SIS_debugging,     only : Bchksum_pair
use SIS_restart,       only : register_restart_field, SIS_restart_CS
use SIS_framework,     only : safe_alloc_ptr
use SIS_hor_grid,      only : SIS_hor_grid_type
use ice_ridging_mod,   only : ridge_rate

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_B_dyn_init, SIS_B_dynamics, SIS_B_dyn_end, SIS_B_dyn_register_restarts

!> The control structure with parameters regulating B-grid ice dynamics
type, public :: SIS_B_dyn_CS ; private
  real, dimension(:,:), pointer :: &
    sig11 => NULL(), &  !< The xx component of the stress tensor [R Z L2 T-2 ~> Pa m] (or N m-1).
    sig12 => NULL(), &  !< The xy and yx component of the stress tensor [R Z L2 T-2 ~> Pa m] (or N m-1).
    sig22 => NULL()     !< The yy component of the stress tensor [RZ Z L2 T-2 ~> Pa m] (or N m-1).

  ! parameters for calculating water drag and internal ice stresses
  real :: p0                  !< Hibbler rheology pressure constant [R L2 T-2 ~> Pa]
  real :: p0_rho              !< The pressure constant divided by ice density [L2 T-2 ~> N m kg-1].
  real :: c0                  !< another pressure constant [nondim]
  real :: cdw                 !< ice/water drag coef. [nondim]
  real :: blturn              !< air/water surf. turning angle (degrees)
  real :: EC = 2.0            !< yield curve axis ratio [nondim]
  real :: MIV_MIN             !< min ice mass to do dynamics [R Z ~> kg m-2]
  real :: Rho_ocean           !< The nominal density of sea water [R ~> kg m-3].
  real :: Rho_ice             !< The nominal density of sea ice [R ~> kg m-3].
  logical :: debug            !< If true, write verbose checksums for debugging purposes.
  logical :: debug_redundant  !< If true, debug redundant points
  integer :: evp_sub_steps    !< The number of iterations in the EVP dynamics
                              !! for each slow time step.
  real    :: dt_Rheo          !< The maximum sub-cycling time step for the rheology
                              !! and momentum equations [T ~> s].
  type(time_type), pointer :: Time => NULL() !< A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                             !! timing of diagnostic output.
  !>@{ Diagnostic IDs
  integer :: id_fix = -1, id_fiy = -1, id_fcx = -1, id_fcy = -1
  integer :: id_fwx = -1, id_fwy = -1, id_sigi = -1, id_sigii = -1
  integer :: id_stren = -1, id_ui = -1, id_vi = -1
  !!@}
end type SIS_B_dyn_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_B_dyn_init initializes the ice dynamics and sets parameters.
subroutine SIS_B_dyn_init(Time, G, US, param_file, diag, CS)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model time.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid type
  type(unit_scale_type),       intent(in)    :: US   !< A structure with unit conversion factors
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(SIS_B_dyn_CS),          pointer       :: CS   !< The control structure for this module

!   This subroutine sets the parameters and registers the diagnostics associated
! with the ice dynamics.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40) :: mdl = "SIS_dyn_bgrid" ! This module's name.
  logical           :: debug
  real, parameter   :: missing = -1e34

  if (.not.associated(CS)) then
    call SIS_error(FATAL, "SIS_B_dyn_init called with an unassociated control structure. \n"//&
                    "SIS_B_dyn_register_restarts must be called before SIS_B_dyn_init.")
    return
  endif

  CS%diag => diag
  CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "DT_RHEOLOGY", CS%dt_Rheo, &
                 "The sub-cycling time step for iterating the rheology "//&
                 "and ice momentum equations. If DT_RHEOLOGY is negative, "//&
                 "the time step is set via NSTEPS_DYN.", units="seconds", scale=US%s_to_T, &
                 default=-1.0)
  CS%evp_sub_steps = -1
  if (CS%dt_Rheo <= 0.0) &
    call get_param(param_file, mdl, "NSTEPS_DYN", CS%evp_sub_steps, &
                 "The number of iterations of the rheology and ice "//&
                 "momentum equations for each slow ice time step.", default=432)

  call get_param(param_file, mdl, "ICE_STRENGTH_PSTAR", CS%p0, &
                 "A constant in the expression for the ice strength, "//&
                 "P* in Hunke & Dukowicz 1997.", &
                 units="Pa", scale=US%kg_m3_to_R*US%m_s_to_L_T**2, default=2.75e4)
  call get_param(param_file, mdl, "ICE_STRENGTH_CSTAR", CS%c0, &
                 "A constant in the exponent of the expression for the "//&
                 "ice strength, c* in Hunke & Dukowicz 1997.", &
                 units="nondim", default=20.)
  call get_param(param_file, mdl, "ICE_CDRAG_WATER", CS%cdw, &
                 "The drag coefficient between the sea ice and water.", &
                 units="nondim", default=3.24e-3)

  call get_param(param_file, mdl, "RHO_OCEAN", CS%Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RHO_ICE", CS%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0, scale=US%kg_m3_to_R)
  CS%p0_rho = CS%p0 / CS%Rho_ice
  !### There should be a new get_param call for this variable
  CS%MIV_MIN = 1.0*US%kg_m3_to_R*US%m_to_Z


  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", default=.false., &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_SLOW_ICE", CS%debug, &
                 "If true, write out verbose debugging data on the slow ice PEs.", &
                 default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_REDUNDANT", CS%debug_redundant, &
                 "If true, debug redundant data points.", default=CS%debug, &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "AIR_WATER_STRESS_TURN_ANGLE", CS%blturn, &
                 "An angle by which to rotate the velocities at the air-water "//&
                 "boundary in calculating stresses.", units="degrees", &
                 default=0.0)


  CS%id_sigi  = register_diag_field('ice_model','SIGI' ,diag%axesT1, Time,     &
            'first stress invariant', 'none', missing_value=missing)
  CS%id_sigii = register_diag_field('ice_model','SIGII' ,diag%axesT1, Time,    &
            'second stress invariant', 'none', missing_value=missing)
  CS%id_stren = register_diag_field('ice_model','STRENGTH' ,diag%axesT1, Time, &
            'ice strength', 'Pa*m', missing_value=missing, conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
  CS%id_fix   = register_diag_field('ice_model', 'FI_X', diag%axesB1, Time,    &
            'ice internal stress - x component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fiy   = register_diag_field('ice_model', 'FI_Y', diag%axesB1, Time,    &
            'ice internal stress - y component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fcx   = register_diag_field('ice_model', 'FC_X', diag%axesB1, Time,    &
            'coriolis force - x component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fcy   = register_diag_field('ice_model', 'FC_Y', diag%axesB1, Time,    &
            'coriolis force - y component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fwx   = register_diag_field('ice_model', 'FW_X', diag%axesB1, Time,    &
            'water stress on ice - x component', 'Pa', conversion=-US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fwy   = register_diag_field('ice_model', 'FW_Y', diag%axesB1, Time,    &
            'water stress on ice - y component', 'Pa', conversion=-US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_ui    = register_diag_field('ice_model', 'UI', diag%axesB1, Time,      &
            'ice velocity - x component', 'm/s', conversion=US%L_T_to_m_s,     &
            missing_value=missing, interp_method='none')
  CS%id_vi    = register_diag_field('ice_model', 'VI', diag%axesB1, Time,      &
            'ice velocity - y component', 'm/s', conversion=US%L_T_to_m_s,     &
            missing_value=missing, interp_method='none')

end subroutine SIS_B_dyn_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> find_ice_strength determines the magnitude of force on ice in plastic deformation
subroutine find_ice_strength(mi, ci, ice_strength, G, US, CS) !, nCat)
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  type(unit_scale_type),            intent(in)  :: US  !< A structure with unit conversion factors
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: mi  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: ice_strength  !< The ice strength [R Z L2 T-2 ~> N m-1]
  type(SIS_B_dyn_CS),               pointer     :: CS  !< The control structure for this module
  ! integer, intent(in) :: nCat !< The number of sea ice categories.

  ! Local variables
  logical :: prs_rothrock = .false.
  !Niki: TOM has a new option for calculating ice strength. If you want to use it set
  !      prs_rothrock=.true.  In that case we need to review and fix the new code
  !      inside if (prs_rothrock), particularly work on getting hi3,ci3,hi3v,Cp & Cf.
  logical :: rdg_lipscomb = .true.
  !Niki: What are ci3,hi3,Cp, Cf?
!  real, dimension(SZI_(G),SZJ_(G),nCat) :: hi3,ci3
!  real, dimension(nCat) :: hi3v
  real :: Cp, Cf

  integer :: i, j, isc, iec, jsc, jec, k
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ! ice strength derived after Rothrock et al. (JGR, 1975)
  if (prs_rothrock) then
    call SIS_error(FATAL, "B-grid find_ice_strength: "//&
                   "The Rothrock et al. strength scheme has not been coded yet.")
!   do j=jsc,jec ; do i=isc,iec

!     ! preparations
!     ice_strength(i,j) = 0.0
!     do k=1,nCat
!       hi3v(k) = hi3(i,j,k)*ci3(i,j,k)
!     enddo
!Niki: The arguments for the following call are missing
!           call ice_ridging_init(nCat, ci3(i,j,:), hi3v, part_undef, part_undef_sum, &
!                hmin, hmax, efold, rdg_ratio)
!      !
!     if (rdg_lipscomb) then   ! based on exponential distribution after Lipscomb
!       do k=1,nCat
!Niki                tmp(k) = hmin(k)*hmin(k) + 2.*hmin(k)*efold(k) + 2.*efold(k)*efold(k)
!       enddo
!     else   ! based on uniform distribution after Thorndike and Hibler
!       do k=1,nCat
!Niki      tmp(k) = ( hmax(k)*hmax(k)*hmax(k) - hmin(k)*hmin(k)*hmin(k) ) / &
!                      ( 3 * (hmax(k)-hmin(k)) )
!       enddo
!     endif
!     !
!     ! ice strength calculation
!     do k=1,nCat
!Niki   ice_strength(i,j) = ice_strength(i,j) + &
!           part_undef(k)/rdg_ratio(k) * (tmp(k)) - &  ! potential energy of new ridges
!           part_undef(k)*hi3(i,j,k)*hi3(i,j,k)        ! potential energy of ridging
!     enddo
!     !if (part_undef_sum>0.0 .and. part_undef_sum<=1.0) then
!     !  ice_strength(i,j) = ice_strength(i,j) * Cf * Cp / part_undef_sum
!     !
!   enddo ; enddo

  else
    ! ice strength derived after Hibler, JGR, 1979 (default!)
    do j=jsc,jec ; do i=isc,iec
      ice_strength(i,j) = CS%p0_rho*mi(i,j)*exp(-CS%c0*(1-ci(i,j)))
    enddo ; enddo
  endif

end subroutine find_ice_strength

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_B_dynamics takes a single dynamics timestep with EVP subcycles
subroutine SIS_B_dynamics(ci, misp, mice, ui, vi, uo, vo, &
     fxat, fyat, sea_lev, fxoc, fyoc, do_ridging, rdg_rate, dt_slow, G, US, CS)

  type(SIS_hor_grid_type),            intent(inout) :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),   intent(in   ) :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)),   intent(in   ) :: misp  !< Mass per unit ocean area of sea ice,
                                                             !! snow and melt pond water [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),   intent(in   ) :: mice  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZIB_(G),SZJB_(G)), intent(inout) :: ui    !< Zonal ice velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJB_(G)), intent(inout) :: vi    !< Meridional ice velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in   ) :: uo    !< Zonal ocean velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in   ) :: vo    !< Meridional ocean velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in   ) :: fxat  !< Zonal air stress on ice [R Z L T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJB_(G)), intent(in   ) :: fyat  !< Meridional air stress on ice [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)),   intent(in   ) :: sea_lev !< The height of the sea level, including
                                                             !! contributions from non-levitating ice from
                                                             !! an earlier time step [Z ~> m].
  real, dimension(SZIB_(G),SZJB_(G)), intent(  out) :: fxoc  !< Zonal ice stress on ocean [R Z L T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJB_(G)), intent(  out) :: fyoc  !< Meridional ice stress on ocean [R Z L T-2 ~> Pa]
  logical,                            intent(in   ) :: do_ridging !< If true, the ice can ridge
  real, dimension(SZI_(G),SZJ_(G)),   intent(  out) :: rdg_rate !< ridging rate from drift state [T-1 ~> s-1]
  real,                               intent(in   ) :: dt_slow !< The amount of time over which the ice
                                                             !! dynamics are to be advanced [T ~> s].
  type(unit_scale_type),              intent(in)    :: US    !< A structure with unit conversion factors
  type(SIS_B_dyn_CS),                 pointer       :: CS    !< The control structure for this module

  ! Local variables
  real, dimension(SZIB_(G),SZJB_(G)) :: fxic, fyic  ! ice internal stresses [R Z L2 T-2 ~> Pa m]
  real, dimension(SZIB_(G),SZJB_(G)) :: fxco, fyco  ! Coriolis force [R Z L T-2 ~> kg m-1 s-2 = Pa]

  real, dimension(SZI_(G),SZJ_(G)) :: prs                    ! Ice internal pressure [R Z L2 T-2 ~> Pa m]
  real                             :: zeta, eta              ! bulk/shear viscosities [R Z L2 T-1 ~> Pa m s]
  real, dimension(SZI_(G),SZJ_(G)) :: strn11, strn12, strn22 ! strain tensor [T-1 ~> s-1]

  real, dimension(SZIB_(G),SZJB_(G)) :: miv                 ! mass on v-points [R Z ~> kg m-2]
  real, dimension(SZIB_(G),SZJB_(G)) :: civ                 ! conc. on v-points [nondim]
  real, dimension(SZI_(G),SZJ_(G))   :: diag_val            ! A temporary diagnostic array
  complex                            :: rr                  ! linear drag rate [T-1 ~> s-1]
  real                               :: fxic_now, fyic_now  ! ice internal stresses [R Z L2 T-2 ~> Pa m]

  logical, dimension(SZI_(G),SZJ_(G))   :: ice_present
  logical :: evp_new = .false.
  real    :: edt_new,e0 !Niki: What is e0
  ! temporaries for strain calculation
  real, dimension(SZI_(G),SZJ_(G)) :: &
       grid_fac1, grid_fac2, grid_fac3, grid_fac4

  ! temporaries for ice stress calculation
  real                             :: del2  ! The squared magnitude of the strain rate ellipse [T-2 ~> s-2]
  real, dimension(SZI_(G),SZJ_(G)) :: edt   ! The elasticity (E) times a time-step [R Z L2 T-1 ~> Pa m s]
  real, dimension(SZI_(G),SZJ_(G)) :: mp4z  ! [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: t0, t1 ! Ratios of viscosities [nondim]
  real, dimension(SZI_(G),SZJ_(G)) :: It2   ! A scaled viscosity [R Z L2 T-1 ~> Pa m s]
  real                             :: a, b  ! Temporary inverse elasticity variables [T R-1 Z-1 L-2 ~> Pa-1 m-1 s-1]
  real                             :: tmp   ! An inverse squared viscosity [T2 R-2 Z-2 L-4 ~> Pa-2 m-2 s-2]
  real                             :: f11, f22 ! Temporary variables [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJB_(G)) :: sldx, sldy  ! Pressure accelerations due to the sloping
                                                    ! sea surface integrated over a timestep [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJB_(G)) :: dydx, dxdy  ! Differences in grid spacings between adjacent points [L ~> m].
  real   :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9

  ! for velocity calculation
  real,    dimension(SZIB_(G),SZJB_(G)) :: dtmiv  ! The timestep divided by the ice mass [T R-1 Z-1 ~> s m2 kg-1]
  real :: dt_Rheo  ! The short timestep associated with the rheology [T ~> s].
  real :: I_2dt_Rheo ! 1.0 / (2*dt_Rheo) [T-1 ~> s-1]
  integer :: EVP_steps ! The number of EVP sub-steps that will actually be taken.
  real :: I_sub_steps
  real :: EC2I    ! 1/EC^2, where EC is the yield curve axis ratio.
  complex                             :: newuv
  real :: pi

  logical :: sent
  integer :: l
  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (.not.associated(CS)) call SIS_error(FATAL, &
       "SIS_B_dynamics: Module must be initialized before it is used.")

  pi = 4.0*atan(1.0)
  EC2I = 1.0/(CS%EC*CS%EC)

  !TOM> derive 2D fields of ice concentration, and snow & ice thicknesses
!Niki: The following are undefined
!  ci = 1-ci3(:,:,1)
!  hs = ice_avg(hs3,ci3)
!  hi = ice_avg(hi3,ci3)

  fxoc(:,:) = 0.0 ; fyoc(:,:) = 0.0 ! zero these for summing later
  fxic(:,:) = 0.0 ; fyic(:,:) = 0.0
  fxco(:,:) = 0.0 ; fyco(:,:) = 0.0

  if ((CS%evp_sub_steps<=0) .and. (CS%dt_Rheo<=0.0)) return

  if (CS%dt_Rheo > 0.0) then
    EVP_steps = max(CEILING(dt_slow/CS%dt_Rheo - 0.0001), 1)
  else
    EVP_steps = CS%evp_sub_steps
  endif
  dt_Rheo = dt_slow/EVP_steps

  do J=jsc-1,jec ; do I=isc-1,iec
    dydx(I,J) = 0.5*((G%dyT(i+1,j+1) - G%dyT(i,j+1)) + (G%dyT(i+1,j) - G%dyT(i,j)))
    dxdy(I,J) = 0.5*((G%dxT(i+1,j+1) - G%dxT(i+1,j)) + (G%dxT(i,j+1) - G%dxT(i,j)))
  enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    grid_fac1(i,j) = (G%dxCv(i,J)-G%dxCv(i,J-1))*G%IdyT(i,j)
    grid_fac2(i,j) = (G%dyCu(I,j)-G%dyCu(I-1,j))*G%IdxT(i,j)
    grid_fac3(i,j) = 0.5*G%dyT(i,j) * G%IdxT(i,j)
    grid_fac4(i,j) = 0.5*G%dxT(i,j) * G%IdyT(i,j)
  enddo ; enddo

  !TOM> check where ice is present
  do j=jsc,jec ; do i=isc,iec
    ice_present(i,j) = ( (G%mask2dT(i,j)>0.5) .and. (misp(i,j) > CS%MIV_MIN) )
  enddo ; enddo

  ! sea level slope force
  do J=jsc-1,jec ; do I=isc-1,iec
    sldx(I,J) = -dt_Rheo * G%g_Earth*(0.5*((sea_lev(i+1,j+1)-sea_lev(i,j+1)) + &
                                           (sea_lev(i+1,j)-sea_lev(i,j)))) * G%IdxBu(i,J)
    sldy(I,J) = -dt_Rheo * G%g_Earth*(0.5*((sea_lev(i+1,j+1)-sea_lev(i+1,j)) + &
                                           (sea_lev(i,j+1)-sea_lev(i,j)))) * G%IdyBu(I,J)
  enddo ; enddo

  ! put ice/snow mass and concentration on v-grid, first finding mass on t-grid.
  do J=jsc-1,jec ; do I=isc-1,iec ; if (G%mask2dBu(i,j) > 0.5 ) then
    miv(I,J) = 0.25*( (misp(i+1,j+1) + misp(i,j)) + (misp(i+1,j) + misp(i,j+1)) )
    civ(I,J) = 0.25*( (ci(i+1,j+1) + ci(i,j)) + (ci(i+1,j) + ci(i,j+1)) )
  else
    miv(I,J) = 0.0 ; civ(I,J) = 0.0
  endif ; enddo ; enddo

  ! precompute prs, elastic timestep parameter, and linear drag coefficient
  !
  call find_ice_strength(mice, ci, prs, G, US, CS)

  !TOM> towards a leaner calculation of the ice stress
  if (evp_new) then
    ! calculate elastic parameter:
    ! E=zeta/(E_0*dt) => E*dt_Rheo=zeta/(E_0*N_evp), where dt_Rheo=dt/N_evp
    ! here, edt_new is 2*zeta/(E*dt_Rheo) = 2*E_0*N_evp for computational reasons
    edt_new = 2. * e0 * REAL(EVP_steps)
  else
    ! This is H&D97, Eq. 44, with their E_0 = 0.25.
    I_2dt_Rheo = 1.0 / (2.0*dt_Rheo)
    do j=jsc,jec ; do i=isc,iec
      if (G%dxT(i,j) < G%dyT(i,j) ) then
        edt(i,j) = I_2dt_Rheo * (G%dxT(i,j)**2 * mice(i,j))
      else
        edt(i,j) = I_2dt_Rheo * (G%dyT(i,j)**2 * mice(i,j))
      endif
    enddo ; enddo
  endif

  do J=jsc-1,jec ; do I=isc-1,iec
    if ((G%mask2dBu(I,J)>0.5) .and. (miv(I,J) > CS%MIV_MIN) ) then ! values for velocity calculation (on v-grid)
      dtmiv(I,J) = dt_Rheo/miv(I,J)
    else
      ui(I,J) = 0.0 ; vi(I,J) = 0.0
    endif
  enddo ; enddo

  if (CS%debug .or. CS%debug_redundant) then
    call Bchksum_pair("sld[xy] in SIS_B_dynamics", sldx, sldy, G, symmetric=.true., scale=US%L_T_to_m_s)
    call Bchksum_pair("f[xy]at in SIS_B_dynamics", fxat, fyat, G, symmetric=.true., &
                      scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    call Bchksum_pair("[uv]i pre-steps SIS_B_dynamics", ui, vi, G, symmetric=.true., scale=US%L_T_to_m_s)
    call Bchksum_pair("[uv]o in SIS_B_dynamics", uo, vo, G, symmetric=.true., scale=US%L_T_to_m_s)
    call Bchksum_pair("d[yx]d[xy] in SIS_B_dynamics", dydx, dxdy, G, scalars=.true., scale=US%L_to_m)
  endif
  if (CS%debug_redundant) then
    call check_redundant_B("civ in SIS_B_dynamics", civ, G)
    call check_redundant_B("miv in SIS_B_dynamics", miv, G)
  endif

  do l=1,EVP_steps

    ! calculate strain tensor for viscosities and forcing elastic eqn.
    call pass_vector(ui, vi, G%Domain, stagger=BGRID_NE)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ! set_strain - calculate generalized orthogonal coordinate strain tensor       !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    do j=jsc,jec ; do i=isc,iec
      strn11(i,j) = (0.5*((ui(I,J)-ui(I-1,J)) + (ui(I,J-1)-ui(I-1,J-1))) + &
           0.25*((vi(I,J)+vi(I-1,J-1)) + (vi(I,J-1)+vi(I-1,J))) * &
           grid_fac1(i,j)) * G%IdxT(i,j)
      strn22(i,j) = (0.5*((vi(I,J)-vi(I,J-1)) + (vi(I-1,J)-vi(I-1,J-1))) + &
           0.25*((ui(I,J)+ui(I,J-1)) + (ui(I-1,J)+ui(I-1,J-1))) * &
           grid_fac2(i,j)) * G%IdyT(i,j)
      strn12(i,j) = 0.5*grid_fac3(i,j) * &
           ( (vi(I,J)*G%IdyBu(I,J) - vi(I-1,J)*G%IdyBu(I-1,J)) + &
           (vi(I,J-1)*G%IdyBu(I,J-1) - vi(I-1,J-1)*G%IdyBu(I-1,J-1)) ) + &
           0.5*grid_fac4(i,j) * &
           ( (ui(I,J)*G%IdxBu(I,J) - ui(I,J-1)*G%IdxBu(I,J-1)) + &
           (ui(I-1,J)*G%IdxBu(I-1,J) - ui(I-1,J-1)*G%IdxBu(I-1,J-1)) )
    enddo ; enddo

    ! calculate viscosities - how often should we do this ?
    if (l>=1) then
      do j=jsc,jec ; do i=isc,iec
        del2 = (strn11(i,j)*strn11(i,j) + strn22(i,j)*strn22(i,j)) * (1+EC2I)     &
             + 4*EC2I*strn12(i,j)*strn12(i,j) + 2*strn11(i,j)*strn22(i,j)*(1-EC2I)  ! H&D eqn 9

        if ( del2 > 4.0e-18*US%T_to_s**2 ) then
          zeta = 0.5*prs(i,j) / sqrt(del2)
        else
          zeta = 2.50e8*US%s_to_T * prs(i,j)
        endif

        ! Hibler uses the following dimensonal constant ceiling to prevent nonlinear instability.
        zeta = max(zeta, 4.0e8*US%kg_m2s_to_RZ_T*US%m_to_L**2)

        eta = zeta*EC2I
        !
        ! some helpful temporaries
        !
        if ( mice(i,j) > 0.0 ) then
          mp4z(i,j) = -prs(i,j)/(4*zeta)
          t0(i,j)   = 2*eta / (2*eta + edt(i,j))
          tmp       = 1/(4*eta*zeta)
          a         = 1/edt(i,j) + (zeta+eta)*tmp ! = 1/edt(i,j) + (1+EC2I)/(4*eta)
          b         = (zeta-eta)*tmp              ! = (1-EC2I)/(4*eta)
          t1(i,j)   = b/a                         ! = (1-EC2I)*edt(i,j) / (4*eta + (1+EC2I)*edt(i,j))
          It2(i,j)  = a / (a**2 - b**2)           ! 1/t2 = a / (a*a - b*b)
        endif
      enddo ; enddo
    endif

    ! timestep stress tensor (H&D eqn 21)
    do j=jsc,jec ; do i=isc,iec
      if( (G%mask2dT(i,j)>0.5) .and. (misp(i,j) > CS%MIV_MIN) ) then
        f11   = mp4z(i,j) + CS%sig11(i,j)/edt(i,j) + strn11(i,j)
        f22   = mp4z(i,j) + CS%sig22(i,j)/edt(i,j) + strn22(i,j)
        CS%sig11(i,j) = (t1(i,j)*f22 + f11) * It2(i,j)
        CS%sig22(i,j) = (t1(i,j)*f11 + f22) * It2(i,j)
        CS%sig12(i,j) = t0(i,j) * (CS%sig12(i,j) + edt(i,j)*strn12(i,j))
      else
        CS%sig11(i,j) = 0.0
        CS%sig22(i,j) = 0.0
        CS%sig12(i,j) = 0.0 ! eliminate internal ice forces
      endif
    enddo ; enddo

    !Niki:
    !Instead of the above block for calculating stresses TOM uses subroutine calls below
    !I have ported these subroutines to this file for completeness.
    !TOM> viscosity and ice stress calculation in subroutines
    !     a) new: EVP solver as documented for CICE 4.0
    !     b) old: solver as in H&D and omsk_2008_03
    !  if (evp_new) then
    !    call ice_stress_new(prs,strn11,strn22,strn12,edt_new, CS%EC,        &
    !                        sig11(isc:iec,jsc:jec), sig22(isc:iec,jsc:jec), &
    !                        sig12(isc:iec,jsc:jec), del2, ice_present       )
    !  else
    !    call ice_stress_old(prs,strn11,strn22,strn12,edt, CS%EC,            &
    !                        sig11(isc:iec,jsc:jec), sig22(isc:iec,jsc:jec), &
    !                        sig12(isc:iec,jsc:jec), del2, ice_present       )
    !  endif
    !


    ! ### SIG11 and SIG22  SHOULD BE PAIRED ON A CUBED SPHERE.
    call pass_var(CS%sig11, G%Domain, complete=.false.)
    call pass_var(CS%sig22, G%Domain, complete=.false.)
    call pass_var(CS%sig12, G%Domain, complete=.true.)

    do J=jsc-1,jec ; do I=isc-1,iec
      if( (G%mask2dBu(i,j)>0.5).and.(miv(i,j)>CS%MIV_MIN)) then ! timestep ice velocity (H&D eqn 22)
        rr = CS%cdw*US%L_to_Z*CS%Rho_ocean * abs(cmplx(ui(i,j)-uo(i,j),vi(i,j)-vo(i,j))) * &
             exp(sign(CS%blturn*pi/180, G%CoriolisBu(i,j)) * (0.0,1.0))
        !
        ! first, timestep explicit parts (ice, wind & ocean part of water stress)
        !
        tmp1 = 0.5*((CS%sig12(i+1,j+1)*G%dxT(i+1,j+1) - CS%sig12(i+1,j)*G%dxT(i+1,j)) + &
                    (CS%sig12(i,j+1)*G%dxT(i,j+1) - CS%sig12(i,j)*G%dxT(i,j)) )
        tmp2 = 0.5*((CS%sig11(i+1,j+1)*G%dyT(i+1,j+1) - CS%sig11(i,j+1)*G%dyT(i,j+1)) + &
                    (CS%sig11(i+1,j)*G%dyT(i+1,j) - CS%sig11(i,j)*G%dyT(i,j)) )
        tmp6 = 0.5*((CS%sig12(i+1,j+1)*G%dyT(i+1,j+1) - CS%sig12(i,j+1)*G%dyT(i,j+1)) + &
                    (CS%sig12(i+1,j)*G%dyT(i+1,j) - CS%sig12(i,j)*G%dyT(i,j)) )
        tmp7 = 0.5*((CS%sig22(i+1,j+1)*G%dxT(i+1,j+1) - CS%sig22(i+1,j)*G%dxT(i+1,j)) + &
                    (CS%sig22(i,j+1)*G%dxT(i,j+1) - CS%sig22(i,j)*G%dxT(i,j)))
        tmp3 = 0.25*((CS%sig12(i+1,j+1)+CS%sig12(i,j)) + (CS%sig12(i+1,j)+CS%sig12(i,j+1)) )
        tmp4 = 0.25*((CS%sig22(i+1,j+1)+CS%sig22(i,j)) + (CS%sig22(i+1,j)+CS%sig22(i,j+1)) )
        tmp5 = 0.25*((CS%sig11(i+1,j+1)+CS%sig11(i,j)) + (CS%sig11(i+1,j)+CS%sig11(i,j+1)) )

        fxic_now = ( (tmp1 + tmp2) + (tmp3*dxdy(I,J) - tmp4*dydx(I,J)) ) * G%IareaBu(I,J)
        fyic_now = ( (tmp6 + tmp7) + (tmp3*dydx(I,J) - tmp5*dxdy(I,J)) ) * G%IareaBu(I,J)

        !### REWRITE TO AVOID COMPLEX EXPRESSIONS.
        ui(I,J) = ui(I,J) + (fxic_now + civ(I,J)*fxat(I,J) + &
             real(civ(I,J)*rr*cmplx(uo(I,J),vo(I,J)))) * dtmiv(I,J) + sldx(I,J)
        vi(I,J) = vi(I,J) + (fyic_now + civ(I,J)*fyat(I,J) + &
             aimag(civ(I,J)*rr*cmplx(uo(I,J),vo(I,J)))) * dtmiv(I,J) + sldy(I,J)

        !
        ! second, timestep implicit parts (Coriolis and ice part of water stress)
        !
        newuv = cmplx(ui(I,J),vi(I,J)) / &
             (1 + dt_Rheo*(0.0,1.0)*G%CoriolisBu(I,J) + civ(I,J)*rr*dtmiv(I,J))
        ui(I,J) = real(newuv); vi(I,J) = aimag(newuv)
        !
        ! sum for averages
        !
        fxic(I,J) = fxic(I,J) + fxic_now
        fyic(I,J) = fyic(I,J) + fyic_now
        fxoc(I,J) = fxoc(I,J) + real( civ(I,J)*rr*cmplx(ui(I,J)-uo(I,J), vi(I,J)-vo(I,J)))
        fyoc(I,J) = fyoc(I,J) + aimag(civ(I,J)*rr*cmplx(ui(I,J)-uo(I,J), vi(I,J)-vo(I,J)))
        fxco(I,J) = fxco(I,J) - miv(I,J)*real ((0.0,1.0)*G%CoriolisBu(I,J) * cmplx(ui(I,J),vi(I,J)))
        fyco(I,J) = fyco(I,J) - miv(I,J)*aimag((0.0,1.0)*G%CoriolisBu(I,J) * cmplx(ui(I,J),vi(I,J)))

      endif
    enddo ; enddo

    if (CS%debug) then
      call hchksum(CS%sig11, "sig11 in SIS_B_dynamics", G%HI, haloshift=1, scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
      call hchksum(CS%sig22, "sig22 in SIS_B_dynamics", G%HI, haloshift=1, scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
      call hchksum(CS%sig12, "sig12 in SIS_B_dynamics", G%HI, haloshift=1, scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)

      call Bchksum(fxic, "fxic in SIS_B_dynamics", G%HI, symmetric=.true., scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call Bchksum(fyic, "fyic in SIS_B_dynamics", G%HI, symmetric=.true., scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call Bchksum(fxoc, "fxoc in SIS_B_dynamics", G%HI, symmetric=.true., scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call Bchksum(fyoc, "fyoc in SIS_B_dynamics", G%HI, symmetric=.true., scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call Bchksum(fxco, "fxco in SIS_B_dynamics", G%HI, symmetric=.true., scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call Bchksum(fyco, "fyco in SIS_B_dynamics", G%HI, symmetric=.true., scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call Bchksum(ui, "ui in SIS_B_dynamics", G%HI, symmetric=.true., scale=US%L_T_to_m_s)
      call Bchksum(vi, "vi in SIS_B_dynamics", G%HI, symmetric=.true., scale=US%L_T_to_m_s)
    endif
    if (CS%debug_redundant) then
      call check_redundant_B("fxic/fyic in SIS_B_dynamics steps",fxic, fyic, G)
      call check_redundant_B("fxco in SIS_B_dynamics steps", fxco, fyco, G)
      call check_redundant_B("fxoc in SIS_B_dynamics steps", fxoc, fyoc, G)
      call check_redundant_B("ui/vi in SIS_B_dynamics steps", ui, vi, G)
    endif

  enddo ! l=1,EVP_steps

  if (CS%debug) then
    call Bchksum(ui, "ui end SIS_B_dynamics", G%HI, symmetric=.true., scale=US%L_T_to_m_s)
    call Bchksum(vi, "vi end SIS_B_dynamics", G%HI, symmetric=.true., scale=US%L_T_to_m_s)
  endif
  if (CS%debug_redundant) &
    call check_redundant_B("ui/vi end SIS_B_dynamics", ui, vi, G)

  ! make averages
  I_sub_steps = 1.0/EVP_steps
  do J=jsc-1,jec ; do I=isc-1,iec
    if ( (G%mask2dBu(i,j)>0.5) .and. miv(i,j)>CS%MIV_MIN ) then
      fxoc(i,j) = fxoc(i,j)*I_sub_steps ; fyoc(i,j) = fyoc(i,j)*I_sub_steps
      fxic(i,j) = fxic(i,j)*I_sub_steps ; fyic(i,j) = fyic(i,j)*I_sub_steps
      fxco(i,j) = fxco(i,j)*I_sub_steps ; fyco(i,j) = fyco(i,j)*I_sub_steps
    endif
  enddo ; enddo

  ! Write out diagnostics associated with the ice dynamics.
  if (query_SIS_averaging_enabled(CS%diag)) then
    if (CS%id_fix>0) call post_SIS_data(CS%id_fix, fxic, CS%diag)
    if (CS%id_fiy>0) call post_SIS_data(CS%id_fiy, fyic, CS%diag)
    if (CS%id_fcx>0) call post_SIS_data(CS%id_fcx, fxco, CS%diag)
    if (CS%id_fcy>0) call post_SIS_data(CS%id_fcy, fyco, CS%diag)
    if (CS%id_fwx>0) call post_SIS_data(CS%id_fwx, fxoc, CS%diag) ! water force on ice
    if (CS%id_fwy>0) call post_SIS_data(CS%id_fwy, fyoc, CS%diag) ! ...= -ice on water
    !  The diagnistics of fxat and fyat are supposed to be taken over all partitions
    !  (ocean & ice), whereas fxat and fyat here are only averaged over the ice.

    if (CS%id_sigi>0) then
      diag_val(:,:) = sigI(mice, ci, CS%sig11, CS%sig22, CS%sig12, G, US, CS)
      call post_SIS_data(CS%id_sigi, diag_val, CS%diag)
    endif
    if (CS%id_sigii>0) then
      diag_val(:,:) = sigII(mice, ci, CS%sig11, CS%sig22, CS%sig12, G, US, CS)
      call post_SIS_data(CS%id_sigii, diag_val, CS%diag)
    endif
    if (CS%id_stren>0) then
      call find_ice_strength(mice, ci, diag_val, G, US, CS)
      call post_SIS_data(CS%id_stren, diag_val, CS%diag)
    endif

    if (CS%id_ui>0) call post_SIS_data(CS%id_ui, ui, CS%diag)
    if (CS%id_vi>0) call post_SIS_data(CS%id_vi, vi, CS%diag)
  endif

  ! Compute the deformation rate for ridging
  if (do_ridging) then
    do j=jsc,jec ; do i=isc,iec ; if (ice_present(i,j)) then
      del2 = (strn11(i,j)*strn11(i,j) + strn22(i,j)*strn22(i,j)) * (1+EC2I) +  &
             4.0*EC2I*strn12(i,j)*strn12(i,j) + 2.0*strn11(i,j)*strn22(i,j)*(1-EC2I) ! H&D eqn 9
      rdg_rate(i,j) = ridge_rate(del2, (strn11(i,j)+strn22(i,j)))
    else
      rdg_rate(i,j) = 0.0
    endif ; enddo ; enddo
  endif

end subroutine SIS_B_dynamics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> sigI evaluates the first stress invariant
function sigI(mi, ci, sig11, sig22, sig12, G, US, CS)
  type(SIS_hor_grid_type),          intent(in) :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: mi  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: sig11 !< The xx component of the stress tensor [RZ Z L2 T-2 ~> Pa m]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: sig22 !< The yy component of the stress tensor [RZ Z L2 T-2 ~> Pa m]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: sig12 !< The xy & yx component of the stress
                                                        !! tensor [RZ Z L2 T-2 ~> Pa m] (or N m-1)
  real, dimension(SZI_(G),SZJ_(G))             :: sigI !< The first stress invariant [nondim]
  type(unit_scale_type),            intent(in) :: US   !< A structure with unit conversion factors
  type(SIS_B_dyn_CS),               pointer    :: CS  !< The control structure for this module

  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call find_ice_strength(mi, ci, sigI, G, US, CS)

  do j=jsc,jec ; do i=isc,iec
    if (sigI(i,j) > 0.0) sigI(i,j) = (sig11(i,j) + sig22(i,j)) / sigI(i,j)
  enddo ; enddo

end function sigI

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> sigII evaluates the second stress invariant
function sigII(mi, ci, sig11, sig22, sig12, G, US, CS)
  type(SIS_hor_grid_type),          intent(in) :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: mi  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: sig11 !< The xx component of the stress tensor [RZ Z L2 T-2 ~> Pa m]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: sig22 !< The yy component of the stress tensor [RZ Z L2 T-2 ~> Pa m]
  real, dimension(SZI_(G),SZJ_(G)), intent(in) :: sig12 !< The xy & yx component of the stress
                                                        !! tensor [RZ Z L2 T-2 ~> Pa m] (or N m-1)
  real, dimension(SZI_(G),SZJ_(G))             :: sigII !< The second stress invariant [nondim]
  type(unit_scale_type),            intent(in) :: US   !< A structure with unit conversion factors
  type(SIS_B_dyn_CS),               pointer    :: CS  !< The control structure for this module

  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call find_ice_strength(mi, ci, sigII, G, US, CS)

  do j=jsc,jec ; do i=isc,iec
    if (sigII(i,j) > 0.0) sigII(i,j) = (((sig11(i,j)-sig22(i,j))**2+4*sig12(i,j)*sig12(i,j))/(sigII(i,j)**2))**0.5
  enddo ; enddo

end function sigII

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_B_dyn_register_restarts allocates and registers any variables for this
!!      module that need to be included in the restart files.
subroutine SIS_B_dyn_register_restarts(HI, param_file, CS, Ice_restart)
  type(hor_index_type),    intent(in) :: HI    !< The horizontal index type describing the domain
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS_B_dyn_CS),      pointer    :: CS    !< The control structure for this module that
                                               !! will be allocated here
  type(SIS_restart_CS),    pointer    :: Ice_restart !< The control structure for the ice restarts

!   This subroutine registers the restart variables associated with the
! the ice dynamics.

  integer :: isd, ied, jsd, jed, id
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_B_dyn_register_restarts called with an "//&
                            "associated control structure.")
    return
  endif
  allocate(CS)

  call safe_alloc_ptr(CS%sig11, isd, ied, jsd, jed)
  call safe_alloc_ptr(CS%sig12, isd, ied, jsd, jed)
  call safe_alloc_ptr(CS%sig22, isd, ied, jsd, jed)

  if (associated(Ice_restart)) then
    call register_restart_field(Ice_restart, 'sig11', CS%sig11, mandatory=.false.)
    call register_restart_field(Ice_restart, 'sig22', CS%sig22, mandatory=.false.)
    call register_restart_field(Ice_restart, 'sig12', CS%sig12, mandatory=.false.)
  endif
end subroutine SIS_B_dyn_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_B_dyn_end - deallocates the memory associated with this module.
subroutine SIS_B_dyn_end(CS)
  type(SIS_B_dyn_CS), pointer :: CS !< The control structure for this module
                                    !! that will be deallocated here.

  deallocate(CS%sig11) ; deallocate(CS%sig12) ; deallocate(CS%sig22)

  deallocate(CS)
end subroutine SIS_B_dyn_end

!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_stress_old derives the ice stress as in SIS of CM2.1 (after Hunke and Dukowicz, 1997)
subroutine ice_stress_old(isc,iec,jsc,jec,prs,strn11,strn22,strn12,edt,EC, &
                          sig11,sig22,sig12,del2,ice_present)
  integer,                          intent(in   ) :: isc !< The starting i-index to work on
  integer,                          intent(in   ) :: iec !< The ending i-index to work on
  integer,                          intent(in   ) :: jsc !< The starting i-index to work on
  integer,                          intent(in   ) :: jec !< The ending j-index to work on
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: prs !< The internal ice pressure [Pa m].
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: strn11 !< The xx component of the strain rate [s-1]
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: strn22 !< The yy component of the strain rate [s-1]
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: strn12 !< The xy & yx component of the strain rate [s-1]
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: edt   !< The ice elasticity times a time-step [Pa m s].
  real,                             intent(in   ) :: EC    !< The yeild curve axis ratio
  real, dimension(isc:iec,jsc:jec), intent(inout) :: sig11 !< The xx component of the stress tensor [N m-1]
  real, dimension(isc:iec,jsc:jec), intent(inout) :: sig22 !< The yy component of the stress tensor [N m-1]
  real, dimension(isc:iec,jsc:jec), intent(inout) :: sig12 !< The xy & yx component of the stress tensor [N m-1]
  real, dimension(isc:iec,jsc:jec), intent(  out) :: del2  !< An elipticity modulated estimate of
                                                           !! the squared strain rate [s-2].
  logical, dimension(isc:iec,jsc:jec), intent(in) :: ice_present !< True where there is any ice present in a cell
  !
  integer                          :: i, j
  real, dimension(isc:iec,jsc:jec) :: mp4z, t0, t1, t2
  real, dimension(isc:iec,jsc:jec) :: f11, f22
  real                             :: a, b, tmp
  real                             :: zeta, eta               ! bulk/shear viscosities
  !
  real :: EC2I

  EC2I = 1.0 / (EC*EC)

  do j=jsc,jec ; do i=isc,iec
    del2(i,j) = (strn11(i,j)*strn11(i,j)+strn22(i,j)*strn22(i,j))*(1+EC2I)     &
                    + 4*EC2I*strn12(i,j)*strn12(i,j) +2*strn11(i,j)*strn22(i,j)*(1-EC2I)  ! H&D eqn 9
    !
    ! calculate viscosities
    if (del2(i,j) > 4e-18 ) then
      zeta = 0.5*prs(i,j)/sqrt(del2(i,j))
    else
      zeta = 2.5e8*prs(i,j)
    endif
    if (zeta<4e8) zeta = 4e8 ! Hibler uses to prevent nonlinear instability

    eta = zeta*EC2I
    !
    ! timestep stress tensor (H&D eqn 21)
    !
    if(ice_present(i,j)) then
 ! some helpful temporaries
      mp4z(i,j) = -prs(i,j)/(4*zeta)
      t0(i,j)   = 2*eta/(2*eta+edt(i,j))
      tmp       = 1/(4*eta*zeta)
      a         = 1/edt(i,j)+(zeta+eta)*tmp
      b         = (zeta-eta)*tmp
      t1(i,j)   = b/a
      t2(i,j)   = a-b*b/a

      f11(i,j)   = mp4z(i,j)+sig11(i,j)/edt(i,j)+strn11(i,j)
      f22(i,j)   = mp4z(i,j)+sig22(i,j)/edt(i,j)+strn22(i,j)
      ! stresses
      sig11(i,j) = (t1(i,j)*f22(i,j)+f11(i,j))/t2(i,j)
      sig22(i,j) = (t1(i,j)*f11(i,j)+f22(i,j))/t2(i,j)
      sig12(i,j) = t0(i,j)*(sig12(i,j)+edt(i,j)*strn12(i,j))
    else
      sig11(i,j) = 0.0
      sig22(i,j) = 0.0
      sig12(i,j) = 0.0 ! eliminate internal ice forces
    endif
  !
  enddo ; enddo
  !
end subroutine ice_stress_old

!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_stress_new derives the ice stress as in CICE 4.0
subroutine ice_stress_new(isc,iec,jsc,jec,prs,strn11,strn22,strn12,edt, EC, &
                          sig11,sig22,sig12,del2,ice_present)
  integer,                          intent(in   ) :: isc !< The starting i-index to work on
  integer,                          intent(in   ) :: iec !< The ending i-index to work on
  integer,                          intent(in   ) :: jsc !< The starting i-index to work on
  integer,                          intent(in   ) :: jec !< The ending j-index to work on
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: prs !< The internal ice pressure [Pa m].
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: strn11 !< The xx component of the strain rate
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: strn22 !< The yy component of the strain rate
  real, dimension(isc:iec,jsc:jec), intent(in   ) :: strn12 !< The xy & yx component of the strain rate
  real,                             intent(in   ) :: edt   !< The ice elasticity times a time-step [Pa m s].
  real,                             intent(in   ) :: EC    !< The yeild curve axis ratio
  real, dimension(isc:iec,jsc:jec), intent(inout) :: sig11 !< The xx component of the stress tensor
  real, dimension(isc:iec,jsc:jec), intent(inout) :: sig22 !< The yy component of the stress tensor
  real, dimension(isc:iec,jsc:jec), intent(inout) :: sig12 !< The xy & yx component of the stress tensor
  real, dimension(isc:iec,jsc:jec), intent(  out) :: del2  !< An elipticity modulated estimate of
                                                           !! the squared strain rate [s-2].
  logical, dimension(isc:iec,jsc:jec), intent(in) :: ice_present !< True where there is any ice present in a cell
  !
  integer :: i, j
  real    :: zeta, eta               ! bulk/shear viscosities
  real    :: strn1, strn2, sig1, sig2
  real    :: tmp
  real :: EC2I

  EC2I = 1.0 / (EC*EC)

  do j=jsc,jec ; do i=isc,iec
    del2(i,j) = (strn11(i,j)*strn11(i,j) + strn22(i,j)*strn22(i,j))*(1+EC2I) + &
                 4*EC2I*strn12(i,j)*strn12(i,j) + 2*strn11(i,j)*strn22(i,j)*(1-EC2I)  ! H&D eqn 9
    !
    ! calculate viscosities
    if (del2(i,j) > 4e-18 ) then
      zeta = 0.5*prs(i,j)/sqrt(del2(i,j))
    else
      zeta = 2.5e8*prs(i,j)
    endif
    if(zeta<4e8) zeta = 4e8 ! Hibler uses to prevent nonlinear instability
    zeta = 2.*zeta
    eta = zeta*EC2I

    ! calculate stresses
    if (ice_present(i,j)) then
      strn1      = strn11(i,j) + strn22(i,j)
      strn2      = strn11(i,j) - strn22(i,j)
      sig1       = sig11(i,j)  + sig22(i,j)
      sig2       = sig11(i,j)  - sig22(i,j)

      sig1       =      edt*sig1       + zeta*strn1 - prs(i,j)
      sig2       = EC2I*edt*sig2       +  eta*strn2
      sig12(i,j) = EC2I*edt*sig12(i,j) +  eta*strn12(i,j)

      sig1       = sig1       * 1./(1.+     edt)
      tmp        =              1./(1.+EC2I*edt)
      sig2       = sig2       * tmp
      sig12(i,j) = sig12(i,j) * tmp

      sig11(i,j) = 0.5*(sig1+sig2)
      sig22(i,j) = 0.5*(sig1-sig2)
    else
      sig11(i,j) = 0.0
      sig22(i,j) = 0.0
      sig12(i,j) = 0.0 ! eliminate internal ice forces
    endif

  enddo ; enddo

end subroutine ice_stress_new

end module SIS_dyn_bgrid
