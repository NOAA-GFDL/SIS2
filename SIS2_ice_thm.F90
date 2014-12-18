!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!                       N-LAYER VERTICAL THERMODYNAMICS                        !
!                                                                              !
! References:                                                                  !
!   Hallberg, R., and M. Winton, 2014:  The SIS2.0 sea-ice model, in prep.     !
!                                                                              !
!   Winton, M., 2011:  A conservative non-iterative n-layer sea ice            !
!     temperature solver, in prep.                                             !
!                                                                              !
!                                                                              !
!         ->+---------+ <- ts - diagnostic surface temperature ( <= 0C )       !
!        /  |         |                                                        !
!      hs   |  snow   | <- tsn   One snow layer with heat capacity             !
!        \  |         |                                                        !
!         =>+---------+                                                        !
!        /  |         |                                                        !
!       /   |         | <- t1    N salty ice layers with heat capacity         !
!      /    |         |                                                        !
!     /     |         | <- t2                                                  !
!   hi      |...ice...|                                                        !
!     \     |         | <- tN-1                                                !
!      \    |         |                                                        !
!       \   |         | <- tN                                                  !
!        \  |         |                                                        !
!         ->+---------+ <- base of ice fixed at seawater freezing temp.        !
!                                                                              !
!                                       Bob Hallberg (Robert.Hallberg@noaa.gov)!
!                                       Mike Winton  (Michael.Winton@noaa.gov) !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module SIS2_ice_thm

! for calling delta-Eddington shortwave from ice_optics
use ice_shortwave_dEdd, only : shortwave_dEdd0_set_snow, shortwave_dEdd0_set_pond
use ice_shortwave_dEdd, only : shortwave_dEdd0, shortwave_dEdd0_set_params
use ice_shortwave_dEdd, only : dbl_kind, int_kind, nilyr, nslyr
use ice_thm_mod, only : get_thermo_coefs
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,  only : get_param, log_param, read_param, log_version, param_file_type

  use constants_mod, only : hlv, hlf ! latent heats of vaporization and fusion.

implicit none ; private

public :: get_thermo_coefs, get_SIS2_thermo_coefs, SIS2_ice_thm_end
public :: SIS2_ice_thm_init, ice_optics_SIS2, ice_temp_SIS2, ice_resize_SIS2
public :: Temp_from_Enth_S, Temp_from_En_S, enth_from_TS, enthalpy_from_TS
public :: enthalpy_liquid_freeze, T_Freeze, calculate_T_Freeze, enthalpy_liquid
public :: e_to_melt_TS, energy_melt_enthS

type, public :: ice_thermo_type ; private
  real :: Cp_ice            ! The heat capacity of ice, in J kg-1 K-1.
  real :: Cp_SeaWater       ! The heat capacity of liquid seawater, in J/(kg K).
  real :: Cp_water          ! The heat capacity of liquid water in the ice model,
                            ! but not in the brine pockets, in J/(kg K).
  real :: Cp_brine          ! The heat capacity of liquid water in the brine
                            ! pockets within the ice, in J/(kg K).  Cp_brine
                            ! should be set equal to Cp_SeaWater, but for
                            ! algorithmic convenience has often been
                            ! set equal to Cp_ice.
  real :: rho_ice, rho_snow, rho_water  ! The nominal densities of ice and water in kg m-3.
  real :: LI                ! The latent heat of fusion, in J kg-1.
  real :: dTf_dS            ! The derivative of the freezing point with salinity, 
                            ! in degC per PSU.  (dTf_dS is negative.)
  real :: mu_TS             ! Negative the derivative of the freezing point with
                            ! salinity, in degC per PSU.

  real :: enth_liq_0 = 0.0     ! The value of enthalpy for liquid fresh
                               ! water at 0 C, in J kg-1.
  real :: enth_unit = 1.0      ! A conversion factor for enthalpy from Joules kg-1.
end type ice_thermo_type

type, public :: SIS2_ice_thm_CS ; private
  ! properties of ice, snow, and seawater (NCAR CSM values)
  real :: KS   ! conductivity of snow, often 0.31 W/(mK)
  real :: KI   ! conductivity of ice, often 2.03 W/(mK)

  ! albedos are from CSIM4 assumming 0.53 visible and 0.47 near-ir insolation
  real :: alb_snow        ! albedo of snow (not melting)
  real :: alb_ice         ! albedo of ice (not melting)
  real :: pen_ice         ! ice surface penetrating solar fraction
  real :: opt_dep_ice     ! ice optical depth (m)
  real :: t_range_melt    ! melt albedos scaled in below melting T
  real :: temp_ice_freeze ! The freezing temperature of the top ice layer, in C.
  real :: temp_range_est  ! An estimate of the range of snow and ice temperatures
                          ! that is used to evaluate whether an explicit
                          ! diffusive form of the heat fluxes or an inversion
                          ! based on the layer heat budget is more likely to
                          ! be the most accurate.

  real :: h_lo_lim        ! hi/hs lower limit for temp. calc.
  real :: frazil_temp_offset = 0.5 ! An offset between the temperature with
                          ! which frazil ice forms and the freezing point of
                          ! each sublayer of the ice.  This functionality could
                          ! later be accounted for using liq_lim instead.

  ! In the ice temperature calculation we place a limit to below (salinity
  ! dependent) freezing point on the prognosed temperatures.  For ice_resize
  ! it is better to make a slightly more restrictive limit that requires the
  ! temperature to be such that the brine content is less than "liq_lim" of
  ! the total mass.  That is T_f/T < liq_lim implying T<T_f/liq_lim
  real :: liq_lim = .99

  logical :: do_deltaEdd = .true.  ! If true, use a delta-Eddington radiative
                          ! transfer calculation for the shortwave radiation
                          ! within the sea-ice and snow.
end type SIS2_ice_thm_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_thm_param - set ice thermodynamic parameters                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS2_ice_thm_init(param_file, CS, ITV )

  use constants_mod, only : hlf ! latent heat of fusion - 334e3 J/(kg-ice)

  type(param_file_type),       intent(in)    :: param_file
  type(SIS2_ice_thm_CS), pointer :: CS
  type(ice_thermo_type), pointer :: ITV ! A pointer to the ice thermodynamic parameter structure.

  real :: sal_ice_top(1)  ! A specified surface salinity of ice.

  real :: deltaEdd_R_ice  ! Mysterious delta-Eddington tuning parameters, unknown.
  real :: deltaEdd_R_snow ! Mysterious delta-Eddington tuning parameters, unknown.
  real :: deltaEdd_R_pond ! Mysterious delta-Eddington tuning parameters, unknown.
  character(len=40)  :: mod = "SIS2_ice_thm" ! This module's name.

  if (.not.associated(CS)) allocate(CS)
  if (.not.associated(ITV)) allocate(ITV)


  ! LI must be taken from the constants mod for internal consistency in the
  ! coupled climate model.
  ITV%LI = hlf

  call get_param(param_file, mod, "RHO_OCEAN", ITV%Rho_water, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0)
  call get_param(param_file, mod, "RHO_ICE", ITV%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0)
  call get_param(param_file, mod, "RHO_SNOW", ITV%Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0)
  call get_param(param_file, mod, "CP_ICE", ITV%Cp_ice, &
                 "The heat capacity of fresh ice, approximated as a \n"//&
                 "constant.", units="J kg-1 K-1", default=2100.0)
  call get_param(param_file, mod, "CP_SEAWATER", ITV%Cp_SeaWater, &
                 "The heat capacity of sea water, approximated as a \n"//&
                 "constant.", units="J kg-1 K-1", default=4200.0)
  call get_param(param_file, mod, "CP_WATER", ITV%Cp_water, &
                 "The heat capacity of water in sea-ice, approximated as \n"//&
                 "a constant.  CP_WATER and CP_SEAWATER should be equal, \n"//&
                 "but for computational convenience CP_WATER has often \n"//&
                 "been set equal to CP_ICE instead.", units="J kg-1 K-1", &
                 default=ITV%Cp_SeaWater)
  call get_param(param_file, mod, "CP_BRINE", ITV%Cp_brine, &
                 "The heat capacity of water in brine pockets within the \n"//&
                 "sea-ice, approximated as a constant.  CP_BRINE and \n"//&
                 "CP_WATER should be equal, but for computational \n"//&
                 "convenience CP_BRINE has often been set equal to CP_ICE.", &
                 units="J kg-1 K-1", default=ITV%Cp_ice)  !### CHANGE LATER TO default=CP_WATER)
  call get_param(param_file, mod, "DTFREEZE_DS", ITV%dTf_dS, &
                 "The derivative of the freezing temperature with salinity.", &
                 units="deg C PSU-1", default=-0.054)
  ITV%mu_TS = -ITV%dTf_dS

  call get_thermo_coefs(ice_salinity=sal_ice_top)
  CS%temp_ice_freeze = T_freeze(sal_ice_top(1), ITV)

  call get_param(param_file, mod, "ENTHALPY_LIQUID_0", ITV%enth_liq_0, &
                 "The enthalpy of liquid fresh water at 0 C.  The solutions \n"//&
                 "should be physically consistent when this is adjusted, \n"//&
                 "because only the relative value is of physical meaning, \n"//&
                 "but roundoff errors can change the solution.", units="J kg-1", &
                 default=0.0)
  call get_param(param_file, mod, "ENTHALPY_UNITS", ITV%enth_unit, &
                 "A constant that rescales enthalpy from J/kg to a \n"//&
                 "different scale in its internal representation.  Changing \n"//&
                 "this by a power of 2 is useful for debugging, as answers \n"//&
                 "should not change.  A negative values is taken as an inverse.", &
                 units="J kg-1", default=1.0)
  if (ITV%enth_unit < 0.) ITV%enth_unit = -1.0 / ITV%enth_unit

  call get_param(param_file, mod, "SNOW_CONDUCTIVITY", CS%Ks, &
                 "The conductivity of heat in snow.", units="W m-1 K-1", &
                 default=0.31)
  call get_param(param_file, mod, "ICE_CONDUCTIVITY", CS%Ki, &
                 "The conductivity of heat in ice.", units="W m-1 K-1", &
                 default=2.03)
  call get_param(param_file, mod, "MIN_H_FOR_TEMP_CALC", CS%h_lo_lim, &
                 "The minimum ice thickness at which to do temperature \n"//&
                 "calculations.", units="m", default=0.0)
 
  call get_param(param_file, mod, "DO_DELTA_EDDINGTON_SW", CS%do_deltaEdd, &
                 "If true, a delta-Eddington radiative transfer calculation \n"//&
                 "for the shortwave radiation within the sea-ice.", default=.true.)
  call get_param(param_file, mod, "ICE_TEMP_RANGE_ESTIMATE", CS%temp_range_est,&
                 "An estimate of the range of snow and ice temperatures \n"//&
                 "that is used to evaluate whether an explicit diffusive \n"//&
                 "form of the heat fluxes or an inversion based on the \n"//&
                 "layer heat budget is more likely to be more accurate. \n"//&
                 "Setting this to 0 causes the explicit diffusive form. \n"//&
                 "to always be used.", units="degC", default=40.0)
  CS%temp_range_est = abs(CS%temp_range_est)

  if (CS%do_deltaEdd) then
    call get_param(param_file, mod, "ICE_DELTA_EDD_R_ICE", deltaEdd_R_ice, &
                   "A dreadfully documented tuning parameter for the radiative \n"//&
                   "propeties of sea ice with the delta-Eddington radiative \n"//&
                   "transfer calculation.", units="perhaps nondimensional?", default=0.0)
    call get_param(param_file, mod, "ICE_DELTA_EDD_R_SNOW", deltaEdd_R_snow, &
                   "A dreadfully documented tuning parameter for the radiative \n"//&
                   "propeties of snow on sea ice with the delta-Eddington \n"//&
                   "radiative transfer calculation.", &
                   units="perhaps nondimensional?", default=0.0)
    call get_param(param_file, mod, "ICE_DELTA_EDD_R_POND", deltaEdd_R_pond, &
                   "A dreadfully documented tuning parameter for the radiative \n"//&
                   "propeties of meltwater ponds on sea ice with the delta-Eddington \n"//&
                   "radiative transfer calculation.", units="perhaps nondimensional?", &
                   default=0.0)
    call shortwave_dEdd0_set_params(deltaEdd_R_ice, deltaEdd_R_snow, deltaEdd_R_pond)

  else
    call get_param(param_file, mod, "SNOW_ALBEDO", CS%alb_snow, &
                   "The albedo of dry snow atop sea ice.", units="nondim", &
                   default=0.85)
    call get_param(param_file, mod, "ICE_ALBEDO", CS%alb_ice, &
                   "The albedo of dry bare sea ice.", units="nondim", &
                   default=0.5826)
    call get_param(param_file, mod, "ICE_SW_PEN_FRAC", CS%pen_ice, &
                   "The fraction of the unreflected shortwave radiation that \n"//&
                   "penetrates into the ice.", units="Nondimensional", default=0.3)
    call get_param(param_file, mod, "ICE_OPTICAL_DEPTH", CS%opt_dep_ice, &
                   "The optical depth of shortwave radiation in sea ice.", &
                   units="m", default=0.67)
    call get_param(param_file, mod, "ALBEDO_T_MELT_RANGE", CS%t_range_melt, &
                   "The temperature range below freezing over which the \n"//&
                   "albedos are changed by partial melting.", units="degC", &
                   default=1.0)
  endif

end subroutine SIS2_ice_thm_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_optics - set albedo, penetrating solar, and ice/snow transmissivity      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_optics_SIS2(hs, hi, ts, tfw, NkIce, alb_vis_dir, alb_vis_dif, &
                    alb_nir_dir, alb_nir_dif, abs_sfc, abs_snow, abs_ice_lay, &
                    abs_ocn, abs_int, CS, coszen_in)
  real, intent(in   ) :: hs  ! snow thickness (m-snow)
  real, intent(in   ) :: hi  ! ice thickness (m-ice)
  real, intent(in   ) :: ts  ! surface temperature
  real, intent(in   ) :: tfw ! seawater freezing temperature
  integer, intent(in) :: NkIce
  real, intent(  out) :: alb_vis_dir ! ice surface albedo (0-1)
  real, intent(  out) :: alb_vis_dif ! ice surface albedo (0-1)
  real, intent(  out) :: alb_nir_dir ! ice surface albedo (0-1)
  real, intent(  out) :: alb_nir_dif ! ice surface albedo (0-1)
  real, intent(  out) :: abs_sfc  ! frac abs sw abs at surface
  real, intent(  out) :: abs_snow ! frac abs sw abs in snow
  real, intent(  out) :: abs_ice_lay(NkIce) ! frac abs sw abs by each ice layer
  real, intent(  out) :: abs_ocn  ! frac abs sw abs in ocean
  real, intent(  out) :: abs_int  ! frac abs sw abs in ice interior
  type(SIS2_ice_thm_CS), intent(in) :: CS
  real, intent(in),optional :: coszen_in

  real :: alb, as, ai, snow_cover, fh
  real :: coalb, I_coalb  ! The coalbedo and its reciprocal.
  real :: SW_frac_top     ! The fraction of the SW at the top of the snow that
                          ! is still present at the top of each ice layer (ND).
  real :: opt_decay_lay   ! The optical extinction in each ice layer (ND).
  real :: pen      ! frac sw passed below the surface (frac 1-pen absorbed at the surface)
  integer :: m
  character(len=200) :: mesg

  integer (kind=int_kind) :: &
    nx_block, ny_block, & ! block dimensions
    icells                ! number of ice-covered grid cells

  integer (kind=int_kind), dimension (1) :: &
    indxi   , & ! compressed indices for ice-covered cells
    indxj

  ! inputs
  real (kind=dbl_kind), dimension (1,1) :: &
    aice   , & ! concentration of ice
    vice   , & ! volume of ice
    vsno   , & ! volume of snow
    Tsfc   , & ! surface temperature
    coszen , & ! cosine of solar zenith angle
    tarea  , & ! cell area - not used
    swvdr  , & ! sw down, visible, direct  (W/m^2)
    swvdf  , & ! sw down, visible, diffuse (W/m^2)
    swidr  , & ! sw down, near IR, direct  (W/m^2)
    swidf      ! sw down, near IR, diffuse (W/m^2)

  ! outputs
  real (kind=dbl_kind), dimension (1,1) :: &
    fs     , & ! horizontal coverage of snow
    fp     , & ! pond fractional coverage (0 to 1)
    hp         ! pond depth (m)

  real (kind=dbl_kind), dimension (1,1,1) :: &
    rhosnw , & ! density in snow layer (kg/m3)
    rsnw       ! grain radius in snow layer (micro-meters)

  real (kind=dbl_kind), dimension (1,1,18) :: &
    trcr        ! aerosol tracers

  real (kind=dbl_kind), dimension (1,1) :: &
    alvdr   , & ! visible, direct, albedo (fraction) 
    alvdf   , & ! visible, diffuse, albedo (fraction) 
    alidr   , & ! near-ir, direct, albedo (fraction) 
    alidf   , & ! near-ir, diffuse, albedo (fraction) 
    fswsfc  , & ! SW absorbed at snow/bare ice/pondedi ice surface (W m-2)
    fswint  , & ! SW interior absorption (below surface, above ocean,W m-2)
    fswthru     ! SW through snow/bare ice/ponded ice into ocean (W m-2)

  real (kind=dbl_kind), dimension (1,1,1) :: &
    Sswabs      ! SW absorbed in snow layer (W m-2)

  real (kind=dbl_kind), dimension (1,1,NkIce) :: &
    Iswabs      ! SW absorbed in ice layer (W m-2)

  real (kind=dbl_kind), dimension (1,1) :: &
    albice  , & ! bare ice albedo, for history  
    albsno  , & ! snow albedo, for history  
    albpnd      ! pond albedo, for history  

  if (CS%do_deltaEdd) then

    if (nilyr /= NkIce) then
      write(mesg, '("The Delta-Eddington sea-ice radiation is hard-coded to use ",(I4),&
                   &" ice layers, not ",(I4),".")') nilyr, NkIce
      call SIS_error(FATAL, mesg)
    endif

    ! temporary for delta-Eddington shortwave call
    nx_block = 1 ; ny_block = 1
    icells = 1 ; indxi(1) = 1 ; indxj(1) = 1
    aice(1,1) = 1.0
    tarea(1,1) = 1.0 ! not used

    ! stuff that matters
    coszen(1,1) = cos(3.14*67.0/180.0) ! NP summer solstice
    if(present(coszen_in))  coszen(1,1) = max(0.01,coszen_in)
    Tsfc(1,1) = ts
    vsno(1,1) = hs
    vice(1,1) = hi
    swvdr(1,1) = 0.25
    swvdf(1,1) = 0.25
    swidr(1,1) = 0.25
    swidf(1,1) = 0.25

    call shortwave_dEdd0_set_snow(nx_block, ny_block, icells, indxi, indxj, &
             aice, vsno, Tsfc, fs, rhosnw, rsnw) ! out: fs, rhosnw, rsnw

    call shortwave_dEdd0_set_pond(nx_block, ny_block, icells, indxi, indxj, &
             aice, Tsfc, fs, fp, hp) ! out: fp, hp
    call shortwave_dEdd0  (nx_block, ny_block, icells, indxi, indxj, coszen, &
             aice, vice, vsno, fs, rhosnw, rsnw, fp, hp, swvdr, swvdf, &
             swidr, swidf, alvdf, alvdr, alidr, alidf, fswsfc, fswint, &
             fswthru, Sswabs, Iswabs, albice, albsno, albpnd)
    ! out: alvdf, alvdr, and subsequent.

    ! Note: fswint = Sswabs + sum(Iswabs)
    alb = 1.0 - (fswsfc(1,1) + (fswint(1,1) + fswthru(1,1)))
    coalb = fswsfc(1,1) + (fswint(1,1) + fswthru(1,1))
    I_coalb = 0.0 ; if (coalb > 0.0) I_coalb = 1.0 / coalb
    abs_sfc  = fswsfc(1,1)   * I_coalb
    abs_snow = Sswabs(1,1,1) * I_coalb
    do m=1,NkIce ; abs_ice_lay(m) = Iswabs(1,1,m) * I_coalb ; enddo
    abs_ocn  = fswthru(1,1)  * I_coalb

    alb_vis_dir = alvdr(1,1)
    alb_vis_dif = alvdf(1,1)
    alb_nir_dir = alidr(1,1)
    alb_nir_dif = alidf(1,1)

    pen = (fswint(1,1) + fswthru(1,1)) * I_coalb
    abs_int = fswint(1,1) * I_coalb

  else
    as = CS%alb_snow ; ai = CS%alb_ice
    snow_cover = hs/(hs+0.02)                ! thin snow partially covers ice

    fh = min(atan(5.0*hi)/atan(5.0*0.5),1.0) ! use this form from CSIM4 to
    ! reduce albedo for thin ice
    if (ts+CS%T_RANGE_MELT > CS%temp_ice_freeze) then        ! reduce albedo for melting as in
       ! CSIM4 assuming 0.53/0.47 vis/ir
       as = as-0.1235*min((ts+CS%T_RANGE_MELT-CS%temp_ice_freeze)/CS%T_RANGE_MELT,1.0)
       ai = ai-0.075 *min((ts+CS%T_RANGE_MELT-CS%temp_ice_freeze)/CS%T_RANGE_MELT,1.0)
    endif
    ai = fh*ai+(1-fh)*0.06                 ! reduce albedo for thin ice

    alb = snow_cover*as + (1-snow_cover)*ai
    alb_vis_dir = alb ; alb_vis_dif = alb
    alb_nir_dir = alb ; alb_nir_dif = alb

    pen = (1-snow_cover)*CS%pen_ice
    opt_decay_lay = exp(-hi/(NkIce*CS%opt_dep_ice))
    abs_ocn = pen * exp(-hi/CS%opt_dep_ice)
    abs_sfc  = 1.0 - pen
    abs_snow = 0.0
    SW_frac_top = pen
    do m=1,NkIce
      abs_ice_lay(m) = SW_frac_top * (1.0 - opt_decay_lay)
      SW_frac_top = SW_frac_top * opt_decay_lay
    enddo
    abs_int = pen - SW_frac_top

     !! check for ice albedos out of range (0 to 1)
     ! if (alb.lt.0.0 .or. alb.gt.1.0) then
     !    print *,'ice_optics: albedo out of range, alb=',alb
     !    print *,'cs=',cs,  'as=',as, 'ai=',ai
     !    print *,'ts=',ts,  'fh=',fh, 'hs=',hs, 'hi=',hi, 'tfw=',tfw
     !    print *,'ALB_SNO=',ALB_SNO,  'ALB_ICE=',ALB_ICE, 'T_RANGE_MELT,=',T_RANGE_MELT, 'TFI=',TFI
     !    stop
     ! end if
  endif

end subroutine ice_optics_SIS2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_temp_SIS2 - A subroutine that calculates the snow and ice enthalpy       !
!    changes due to surface forcing.                                           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_temp_SIS2(m_snow, m_ice, enthalpy, sice, sh_T0, B, sol, tfw, fb, &
                         tsurf, dtt, NkIce, tmelt, bmelt, CS, ITV, check_conserve)

  real, intent(in   ) :: m_snow  ! snow mass per unit area (H, usually kg m-2)
  real, intent(in   ) :: m_ice   ! ice mass per unit area (H, usually kg m-2)
  real, dimension(0:NkIce) , &
        intent(inout) :: enthalpy ! The enthalpy of each layer in a column of
                                  ! snow and ice, in enth_unit (J kg-1).
  real, dimension(NkIce), &
        intent(in)    :: Sice  ! ice salinity by layer (g/kg)
  real, intent(in   ) :: sh_T0 ! net surface heat flux (+ up) at ts=0 (W/m^2)
  real, intent(in   ) :: B     ! d(sfc heat flux)/d(ts) [W/(m^2 deg-C)]
  real, dimension(0:NkIce), &
        intent(in)    :: sol   ! Solar heating of the snow and ice layers (W m-2)
  real, intent(in   ) :: tfw   ! seawater freezing temperature (deg-C)
  real, intent(in   ) :: fb    ! heat flux upward from ocean to ice bottom (W/m^2)
  real, intent(  out) :: tsurf ! surface temperature (deg-C)
  real, intent(in   ) :: dtt   ! timestep (sec)
  integer, intent(in   ) :: NkIce ! The number of ice layers.
  real, intent(inout) :: tmelt ! accumulated top melting energy  (J/m^2)
  real, intent(inout) :: bmelt ! accumulated bottom melting energy (J/m^2)
  type(SIS2_ice_thm_CS), intent(in) :: CS
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  logical, optional, intent(in) :: check_conserve ! If true, check for local heat conservation.
!
! variables for temperature calculation [see Winton (1999) section II.A.]
! note:  here equations are multiplied by hi to improve thin ice accuracy
!
  real :: A ! Net downward surface heat flux from the atmosphere at 0C (W/m^2)
!  real, dimension(0:NkIce) :: &
!    temp_est, &    ! An estimated snow and ice temperature, in degC.
!    temp_IC, &     ! The temperatures of the snow and ice based on the initial
!                   ! enthalpy, in degC.
!    temp_new       ! The updated temperatures, in degC.
  real, dimension(0:NkIce) :: temp_est   ! An estimated snow and ice temperature, in degC.
  real, dimension(0:NkIce) :: temp_IC    ! The temperatures of the snow and ice based on the initial
                                         ! enthalpy, in degC.
  real, dimension(0:NkIce) :: temp_new   ! The updated temperatures, in degC.
  real, dimension(NkIce) :: tfi  ! The ice freezing temperatures, in degC.
  real :: mL_ice   ! The mass-per-unit-area of each ice layer in kg m-2 (not H).
  real :: mL_snow  ! The mass-per-unit-area of each snow layer in kg m-2 (not H).
  real :: e_extra
  real, dimension(0:NkIce) :: m_lay     ! Masses of all layers in kg m-2.
  real :: enth_fp  ! The enthalpy at the freezing point (solid for fresh ice).
  real :: kk, k10, k0a, k0skin
  real :: I_bb, b_denom_1
  real :: comp_rat ! The complement of rat, going from 0 to 1.
  real :: tsf      ! The surface freezing temperature in degC.
  real :: k0a_x_ta, tsno_est, salt_part, rat
  real :: tsurf_est  ! An estimate of the surface temperature in degC.
  real, dimension(0:NkIce+1) :: cc ! Interfacial coupling coefficients.
  real, dimension(0:NkIce) :: bb   ! Effective layer heat capacities.
  real, dimension(0:NkIce) :: cc_bb ! Remaining coupling ratios.
  real, dimension(-1:NkIce) :: heat_flux_int ! The downward heat fluxes at the
                                ! interfaces between layers, in W m-2.
                                ! heat_flux_int uses the index convention from
                                ! MOM6 that interface K is below layer k.
  real :: I_liq_lim     ! The inverse of CS%liq_lim.
  real :: heat_flux_err_rat
  real :: col_enth1, col_enth2, col_enth2b, col_enth3
  real :: d_e_extra, e_extra_sum
  real :: tflux_bot, tflux_bot_diff, tflux_bot_resid
  real :: tfb_diff_err, tfb_resid_err
  real :: tflux_sfc, sum_sol, d_tflux_bot
  real :: hsnow_eff
  real :: snow_temp_new, snow_temp_max
  real :: hL_ice_eff
  real :: enth_liq_lim
  real :: enth_prev
  real :: I_enth_unit
  logical :: col_check
  integer :: k

  col_check = .false. ; if (present(check_conserve)) col_check = check_conserve

  temp_IC(0) = temp_from_En_S(enthalpy(0), 0.0, ITV)
  call temp_from_Enth_S(enthalpy(1:), sice(1:), temp_IC(1:), ITV)

  A = -sh_T0

  I_enth_unit = 1.0 / ITV%enth_unit
  mL_ice = m_ice / NkIce   ! ice mass per unit area of each layer
  mL_snow = m_snow         ! snow mass per unit area (in kg m-2).
  call calculate_T_Freeze(sice, tfi, ITV)    ! freezing temperature of ice layers

  ! Set the effective thickness of each ice and snow layer, limited to avoid
  ! instabilities for thin layers.
  hL_ice_eff = max(mL_ice / ITV%Rho_ice, CS%H_LO_LIM) 
  hsnow_eff = mL_snow / ITV%Rho_snow + max(1e-35, 1e-20*CS%H_LO_LIM)

  kk = CS%KI/hL_ice_eff       ! full ice layer conductivity

  tsf = tfi(1)                ! surface freezing temperature
  if (mL_snow>0.0) tsf = 0.0
  
  k10 = 2.0*(CS%KS*CS%KI) / (hL_ice_eff*CS%KS + hsnow_eff*CS%KI) ! coupling ice layer 1 to snow
  k0a = (CS%KS*B) / (0.5*B*hsnow_eff + CS%KS)      ! coupling snow to "air"
  k0skin = 2.0*CS%KS / hsnow_eff
  k0a_x_ta = (CS%KS*A) / (0.5*B*hsnow_eff + CS%KS) ! coupling times "air" temperture

  enth_liq_lim = Enth_from_TS(0.0, 0.0, ITV)
  
  ! Determine the enthalpy for conservation checks.
  m_lay(0) = mL_snow ; do k=1,NkIce ; m_lay(k) = mL_ice ; enddo
 
  if (col_check) then
    col_enth1 = 0.0
    do k=0,NkIce ; col_enth1 = col_enth1 + m_lay(k)*enthalpy(k) ; enddo
  endif
  e_extra_sum = 0.0

  !   First get a non-conservative estimate with a linearized implicit treatment
  ! of layer coupling.
  
  ! Determine the effective layer heat capacities.
  !   bb = dheat/dTemp.  This should be a proper linearization of the enthalpy equation.
  bb(0) = mL_snow*ITV%Cp_ice
  do k=1,NkIce   ! load bb with heat capacity term.
    salt_part = 0.0
    if (sice(k)>0.0) salt_part = tfi(k)*ITV%LI/(temp_IC(k)*temp_IC(k))
    !### ADD TERM TO ACCOUNT FOR DIFFERENCE BETWEEN CP_ICE AND CP_Brine.
    bb(k) = mL_ice*(ITV%Cp_ice-salt_part) ! add coupling to this later
  enddo

  cc(0) = k0a*dtt  ! Atmosphere-snow coupling
  cc(1) = k10*dtt  ! Snow-ice coupling
  do k=2,NkIce ; cc(k) = kk*dtt ; enddo
  cc(NkIce+1) = 2*kk*dtt ! bottom of ice coupling with ocean at freezing point.

  !   This is a version of a tridiagonal solver where the relationship between
  ! the diagonal elements and the coupling between layers has been used to
  ! ensure that there are no differences ever taken, so there will be minimal
  ! truncation errors.  This closely follows schemes that have previously been
  ! found to work very well in GOLD and MOM6.

  ! Go UP the ice column.
  b_denom_1 = (bb(NkIce) + cc(NkIce+1))
  I_bb = 1.0 / (b_denom_1 + cc(NkIce))
  temp_est(NkIce) = ( (sol(NkIce)*dtt + bb(NkIce)*temp_IC(NkIce)) + &
                      cc(NkIce+1)*tfw ) * I_bb
  comp_rat = b_denom_1 * I_bb
  cc_bb(NkIce) = cc(NkIce) * I_bb

  do k=NkIce-1,1,-1
    b_denom_1 = bb(k) + comp_rat*cc(k+1)
    I_bb =  1.0 / (b_denom_1 + cc(k))
    temp_est(k) = ((sol(k)*dtt + bb(k)*temp_IC(k)) + cc(k+1)*temp_est(k+1)) * I_bb

    comp_rat = b_denom_1 * I_bb ! 1.0 >= comp_rat >= 0.0
    cc_bb(k) = cc(k) * I_bb
  end do

  b_denom_1 = bb(0) + comp_rat*cc(1)
  I_bb =  1.0 / (b_denom_1 + cc(0))
  !   This is a complete calculation of temp_est(0), assuming that the surface
  ! flux is given by A - B*tsurf, with tsurf estimated as part of this calculation.
  temp_est(0) = (((sol(0)*dtt + bb(0)*temp_IC(0)) + k0a_x_ta*dtt) + cc(1)*temp_est(1)) * I_bb
  
  ! Diagnose the surface skin temperature by matching the diffusive fluxes in
  ! the snow with the atmospheric fluxes.  I.e. solve the following for tsurf_est:
  !  (A - B*tsurf_est) = k0skin * (tsurf_est - temp_est(0))
  tsurf_est = (A + k0skin*temp_est(0)) / (B + k0skin)

  if (tsurf_est > tsf) then
    ! The surface is melting, set tsurf to melt temp. and recalculate I_bb.
    tsurf_est = tsf
    ! cc(0) = k0skin*dtt

    I_bb =  1.0 / (b_denom_1 + k0skin*dtt)
    temp_est(0) = min(tsf, &
      (((sol(0)*dtt + bb(0)*temp_IC(0)) + k0skin*dtt*tsf) + cc(1)*temp_est(1)) * I_bb)
  endif
  ! Go back DOWN the ice column to get the estimated temperatures, subject to the
  ! limitation that all temperatures are assumed to be at or below freezing.
  do k=1,NkIce
    temp_est(k) = min(temp_est(k) + cc_bb(k) * temp_est(k-1), tfi(k))
  end do

  !
  ! Quasi-conservative iterative pass going UP the ice column
  !
  temp_est(NkIce) = laytemp_SIS2(mL_ice, tfi(NkIce), sol(NkIce) + kk*(2*tfw+temp_est(NkIce-1)), &
                                 3*kk, temp_IC(NkIce), dtt, ITV)
  do k=NkIce-1,2,-1
    temp_est(k) = laytemp_SIS2(mL_ice, tfi(k), sol(k) + kk*(temp_est(k-1)+temp_est(k+1)), &
                               2*kk, temp_IC(k), dtt, ITV)
  enddo
  temp_est(1) = laytemp_SIS2(mL_ice, tfi(1), sol(1) + (kk*temp_est(2) + k10*temp_est(0)), &
                             kk + k10, temp_IC(1), dtt, ITV)
  
  ! Calculate the bulk snow temperature and surface skin temperature together.
  temp_est(0) = laytemp_SIS2(mL_snow, 0.0, sol(0) + (k10*temp_est(1)+k0a_x_ta), &
                          k10+k0a, temp_IC(0), dtt, ITV)
  tsurf = (A + k0skin*temp_est(0)) / (B + k0skin)  ! diagnose surface skin temp.

  !
  !   The following conservative update pass going DOWN the ice column is where
  ! the layer enthalpies are actually updated and extra heat that should drive
  ! melting is accumulated and stored.
  !

  !   Pre-calculate a conversion factor that can be used to figure out whether it
  ! is more accurate to use the explicit form for the fluxes or to invert enthalpy
  ! conservation for the fluxes, depending on the relative layer masses and the 
  ! strength of the coupling between layers.  The expression below is a simplified
  ! form that assume that temperatures are measured in Celsius and are less than
  ! about CS%temp_err_est, and that the heat budget has 4 terms of comparable
  ! magnitude to the enthalpy content of the layer.  If there were a larger
  ! tolerance for the iterative estimate of a new temperature, that would go into
  ! the numerator but not the denominator.  The present form is based on the
  ! assumption that floating-point roundoff dominates the errors.
  heat_flux_err_rat = 0.7071 * dtt * (CS%temp_range_est) / (CS%temp_range_est * ITV%Cp_Ice + ITV%LI)

  e_extra = 0.0
  if (tsurf > tsf) then ! The surface is melting: update enthalpy with the surface at melt temp.
    tsurf = tsf
    ! Accumulate surface melt energy.
    if (mL_snow>0.0) then
      heat_flux_int(-1) = k0skin * tsf
      heat_flux_int(0) = -k10*temp_est(1)
      call update_lay_enth(mL_snow, 0.0, enthalpy(0), heat_flux_int(-1), &
                           sol(0), heat_flux_int(0), -k0skin, k10, dtt, &
                           heat_flux_err_rat, ITV, e_extra)

      tmelt = tmelt + e_extra + dtt*((A-B*tsf) - heat_flux_int(-1))
      tflux_sfc = dtt*heat_flux_int(-1)
      e_extra_sum = e_extra_sum + e_extra
    else
      ! There is no snow mass, so just convert tsnow = tsf to enthalpy.
      enthalpy(0) = enth_from_TS(tsf, 0.0, ITV)

      heat_flux_int(0) = k10*(tsf - temp_est(1))
      heat_flux_int(-1) = heat_flux_int(0)

      ! Replace use tsurf in the calculation of tmelt
      tmelt = tmelt + dtt*((sol(0)+(A-B*tsf)) - heat_flux_int(0))
      tflux_sfc = dtt*heat_flux_int(0)

    endif
  else
    heat_flux_int(-1) = k0a_x_ta
    heat_flux_int(0) = -k10*temp_est(1)
    snow_temp_max = (tsf*(B + k0skin) - A) / k0skin
    call update_lay_enth(mL_snow, 0.0, enthalpy(0), heat_flux_int(-1), &
                         sol(0), heat_flux_int(0), -k0a, k10, dtt, &
                         heat_flux_err_rat, ITV, e_extra, &
                         temp_new=snow_temp_new, temp_max=snow_temp_max)
    tsurf = (A + k0skin*snow_temp_new) / (B + k0skin)  ! diagnose surface skin temp.
    ! This is equivalent to, but safer than, tsurf = (A - heat_flux_int(-1)) / B

    tflux_sfc = dtt*heat_flux_int(-1)
    e_extra_sum = e_extra_sum + e_extra
    tmelt = tmelt + e_extra
  endif

  do k=1,NkIce-1 ! flux from above is fixed, only have downward feedback
    heat_flux_int(K) = -kk*temp_est(k+1)
    call update_lay_enth(mL_ice, Sice(k), enthalpy(k), heat_flux_int(K-1), &
                         sol(k), heat_flux_int(K), 0.0, kk, dtt, &
                         heat_flux_err_rat, ITV, e_extra)

    e_extra_sum = e_extra_sum + e_extra
    if (k <= NkIce/2) then ; tmelt = tmelt + e_extra
    else ; bmelt = bmelt + e_extra ; endif
  enddo
  
  heat_flux_int(NkIce) = -2.0*kk*tfw
  call update_lay_enth(mL_ice, Sice(NkIce), enthalpy(NkIce), heat_flux_int(NkIce-1), &
                       sol(NkIce), heat_flux_int(NkIce), 0.0, 2.0*kk, dtt, &
                       heat_flux_err_rat, ITV, e_extra)
  e_extra_sum = e_extra_sum + e_extra 
  bmelt = bmelt + e_extra
  !
  ! END of the conservative update of enthalpy.
  !

  if (col_check) then
    col_enth2 = e_extra_sum*ITV%enth_unit ; sum_sol = 0.0 ; col_enth2b = 0.0
    do k=0,NkIce
      col_enth2 = col_enth2 + m_lay(k)*enthalpy(k)
      col_enth2b = col_enth2b + m_lay(k)*enthalpy(k)
      sum_sol = sum_sol + dtt*sol(k)
    enddo 
    !   tflux_bot_resid and tflux_bot_diff are two mathematically equivalent
    ! estimates of the heat flux at the base of the ice.
    tflux_bot_resid = (col_enth2 - col_enth1) - (sum_sol + tflux_sfc)
    tflux_bot_diff = -heat_flux_int(NkIce)*dtt

    ! Estimate the errors with these two expressions from 64-bit roundoff.
    tfb_diff_err = 1e-15*2.0*kk*dtt * sqrt(tfw**2 + 10.0**2)  ! The -10 deg is arbitrary but good enough?
    tfb_resid_err = 1e-15*sqrt(col_enth2**2 + col_enth1**2 + sum_sol**2 + tflux_sfc**2)

    d_tflux_bot = tflux_bot_diff - tflux_bot_resid
    if (abs(d_tflux_bot) > 1.0e-9*(abs(tflux_bot_resid) + abs(tflux_bot_diff) + &
                                   abs(col_enth2 - col_enth1))) then
      d_tflux_bot = tflux_bot_diff - tflux_bot_resid
    endif
  endif

  tflux_bot = -heat_flux_int(NkIce)*dtt

  !   Accumulate the difference between the ocean's heat flux to the ice-ocean
  ! interface and the sea-ice heat flux to evaluate bottom melting/freezing.
  bmelt = bmelt + (dtt*fb - tflux_bot)

  e_extra_sum = 0.0

  if (enthalpy(0) > enth_liq_lim) then ! put excess snow energy into top melt.
    e_extra = (enthalpy(0) - enth_liq_lim) * mL_snow * I_enth_unit
    tmelt = tmelt + e_extra
    e_extra_sum = e_extra_sum + e_extra
    enthalpy(0) = enth_liq_lim
  endif

  I_liq_lim = 1.0 / CS%liq_lim
  do k=1,NkIce
    enth_liq_lim = enth_from_TS(tfi(k)*I_liq_lim, sice(k), ITV)
    if (enthalpy(k) > enth_liq_lim) then ! push excess energy to closer of top or bottom melt
      e_extra = (enthalpy(k) - enth_liq_lim) * mL_ice * I_enth_unit
      e_extra_sum = e_extra_sum + e_extra
      enthalpy(k) = enth_liq_lim
      if (k<=NkIce/2) then
        tmelt = tmelt+e_extra
      else
        bmelt = bmelt+e_extra
      endif
    endif
  enddo

  if (col_check) then
    col_enth3 = 0.0
    do k=0,NkIce ; col_enth3 = col_enth3 + m_lay(k)*enthalpy(k) ; enddo

    d_e_extra = col_enth3 - (col_enth2b - e_extra_sum*ITV%enth_unit)
    if (abs(d_e_extra) > 1.0e-12*(abs(col_enth3) + abs(col_enth2b) + &
                                  abs(e_extra_sum*ITV%enth_unit))) then
      d_e_extra = (col_enth3 - col_enth2b) - e_extra_sum*ITV%enth_unit
    endif
  endif

  call ice_check(mL_snow, NkIce*mL_ice, enthalpy, sice, NkIce, &
           "at end of ice_temp_SIS2", ITV, bmelt=bmelt, tmelt=tmelt, t_sfc=tsurf)

end subroutine ice_temp_SIS2

!
! laytemp_SIS2 - implicit calculation of new layer temperature
!
function laytemp_SIS2(m, tfi, f, b, tp, dtt, ITV) result (new_temp)
  real ::  new_temp
  real, intent(in) :: m    ! mass of ice - kg/m2
  real, intent(in) :: tfi  ! ice freezing temp. (determined by salinity)
  real, intent(in) :: f    ! Inward forcing - W/m2
  real, intent(in) :: b    ! response of outward heat flux to local temperature - W/m2/K
  real, intent(in) :: tp   ! prior step temperature
  real, intent(in) :: dtt  ! timestep in s.
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.

  real :: E0   ! Starting heat relative to salinity dependent freezing.
  real :: AA, BB, CC
  real :: Cp_Ice, LI
  Cp_Ice = ITV%Cp_Ice ; LI = ITV%LI

  if ( tfi == 0.0 ) then
    ! For fresh water, avoid the degeneracy of the enthalpy-temperature
    ! relationship by extending the linear expression for frozen water, then
    ! limit it later to be at or below freezing.
    !
    !   m * {Cp_Ice} * (tn-tp) = dtt * (f - b*tn)
    !
    new_temp = (m*Cp_Ice*tp + f*dtt) / (m*Cp_Ice + b*dtt) ! = -BB/AA

  elseif (ITV%Cp_ice == ITV%Cp_brine) then
    if (tp >= tfi) then
      E0 = ITV%Cp_Water*(tp - tfi)  ! >= 0
    else
      E0 = Cp_Ice*(tp - tfi) - LI*(1 - tfi/tp)  ! < 0
    endif
    ! Determine whether the new solution will be above or below freezing.
    
    if (m*E0 + dtt * (f - b*tfi) >= 0) then
      ! This layer will be completely melted.
      new_temp = tfi + (m*E0 + dtt* (f - b*tfi)) / (Cp_Ice*m + dtt*b)
    else
      ! This layer will be partly melted.
      ! Solve a quadratic equation for the new layer temperature, tn:
      !
      !   m * {Cp_Ice-LI*tfi/(tn*tp)} * (tn-tp) = dtt * (f - b*tn)
      !
      AA = m*Cp_Ice + b*dtt
      BB = -(m*((E0 + LI) + Cp_Ice*tfi) + f*dtt)
      CC = m*LI*tfi
      ! This form avoids round-off errors.
      if (BB >= 0) then
        new_temp = -(BB + sqrt(BB*BB - 4*AA*CC)) / (2*AA)
      else
        new_temp = (2*CC) / (-BB + sqrt(BB*BB - 4*AA*CC))
      endif
    endif
  else
    call SIS_error(FATAL, "Write laytemp_SIS2 for Cp_ice /= Cp_brine.")
  endif

  ! Only return temperatures that are at or below the freezing point.
  new_temp = min(new_temp, tfi)

end function laytemp_SIS2

!
! update_lay_enth - implicit calculation of new layer enthalpy
!
subroutine update_lay_enth(m_lay, sice, enth, ftop, ht_body, fbot, dftop_dT, &
                           dfbot_dT, dtt, hf_err_rat, ITV, extra_heat, temp_new, temp_max)
  real, intent(in) :: m_lay    ! This layers mass of ice in kg/m2
  real, intent(in) :: sice     ! ice salinity in g/kg
  real, intent(inout) :: enth  ! ice enthalpy in enth_units (proportional to J kg-1).
  real, intent(inout) :: ftop  ! Downward heat flux atop the layer in W/m2.
  real, intent(in) :: ht_body  ! Body forcing  to layer in W/m2
  real, intent(inout) :: fbot  ! Downward heat below the layer in W/m2.
  real, intent(in) :: dftop_dT ! The linearization of ftop with layer temperature in W m-2 K-1.
  real, intent(in) :: dfbot_dT ! The linearization of fbot with layer temperature in W m-2 K-1.
  real, intent(in) :: dtt      ! The timestep in s.
  real, intent(in) :: hf_err_rat  ! A conversion factor for comparing the errors
                               ! in explicit and implicit estimates of the updated
                               ! heat fluxes, in (kg m-2) / (W m-2 K-1).
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real, intent(out) :: extra_heat ! The heat above the melt point, in J.
  real, optional, intent(out) :: temp_new ! The new temperature, in degC.
  real, optional, intent(in)  :: temp_max ! The maximum new temperature, in degC.

  real :: htg      ! The rate of heating of the layer in W m-2.
  real :: new_temp ! The new layer temperature, in degC.
  real :: max_temp ! The maximum new layer temperature, in degC.
  real :: max_enth ! The maximum new layer enthalpy, in degC.
  real :: fb       ! The negative of the dependence of layer heating on
                   ! temperature, in W m-2 K-1. fb > 0.
  real :: extra_enth ! Excess enthalpy above the melt point, in kg enth_units.
  real :: enth_in  ! The initial enthalpy, in enth_units.
  real :: enth_fp  ! The enthalpy at the freezing point, in enth_units.
  real :: AA, BB, CC ! Temporary variables used to solve a quadratic equation.
  real :: dtEU     ! The timestep times the unit conversion from J to Enth_units, in s?
  real :: dT_dEnth ! The partial derivative of temperature with enthalpy,
                   ! in units of K / Enth_unit.
  real :: En_J     ! The enthalpy in Joules with 0 offset for liquid at 0 C.
  real :: tfi      ! ice freezing temp. (determined by salinity)
  real :: fbot_in, ftop_in ! Input values of fbot and ftop in W m-2.
  real :: dflux_dtot_dT  ! A temporary work array in units of degC.

  ! Solve m_lay*(enth - enth_in) + extra_heat = dt * (ht_body + ftop - fbot)
  !  ftop = ftop_in + temp*dftop_dT
  !  fbot = fbot_in + temp*dfbot_dT
  !    enth <= enth_fp and extra_heat >= 0

  ftop_in = ftop ; fbot_in = fbot
  htg = (ht_body + ftop_in) - fbot_in
  fb = -(dftop_dT - dfbot_dT)   !  = -dhgt_dt > 0

  extra_heat = 0.0 ; extra_enth = 0.0
  if (sice > 0.0) then
    tfi = T_freeze(sice, ITV)
    enth_fp = enthalpy_liquid_freeze(sice, ITV)
  else
    tfi = 0.0
    enth_fp = enth_from_TS(0.0, 0.0, ITV)
  endif
  max_temp = tfi ; max_enth = enth_fp
  if (present(temp_max)) then ; if (temp_max < tfi) then
    max_temp = temp_max ; max_enth = enth_from_TS(temp_max, sice, ITV)
  endif ; endif
  enth_in = enth
  dtEU = ITV%enth_unit * dtt
  
  ! Solve m_lay * (enth_new - enth) = dtEU * (htg - fb*t_new)
  !       t_new = Temp_from_En_S(enth_new, Sice, ITV)
  if (m_lay == 0.0) then
    new_temp = min(htg / fb, max_temp)
    enth = enth_from_TS(new_temp, sice, ITV)
  elseif (dtEU * (htg - fb*max_temp) >= m_lay*(max_enth - enth_in)) then
    ! There is enough heat being applied here that the ice would be above the
    ! maximum temperature (often the freezing point).  The ice should be set to
    ! the maximum temperature and enthalpy, and the extra heat stored for later
    ! use in melting.
    extra_enth = m_lay*(enth_in - max_enth) + dtEU * (htg - fb*max_temp)
    extra_heat = extra_enth / ITV%enth_unit
    new_temp = max_temp
    enth = max_enth
  elseif ( sice == 0.0 ) then  ! Note that tfi = 0.
    ! dT_dEnth is 0 for enth > enth_fp.
    !   dT_dEnth = dTemp_dEnth(enth_in, Sice, ITV)
    dT_dEnth = 1.0 / (ITV%Cp_Ice * ITV%enth_unit)

    ! Solve for enth:  m_lay  * (enth - enth_in) = 
    !       dtEU * (htg - fb*tfi - fb*dT_dEnth*(enth - enth_fp))
    !  enth = enth_in + dtEU * (htg - fb*(0.0 - dT_dEnth*(enth_fp-enth_in))) / &
    !                          (m_lay + dtEU*b*dT_dEnth)
    ! Or equivalently...  (noting that tfi = 0.0)
    enth = enth_fp + (dtEU * (htg - fb*0.0) + m_lay * (enth_in-enth_fp)) / &
                     (m_lay  + dtEU*(fb*dT_dEnth))
    ! The following is equivalent to new_temp = Temp_from_En_S(enth, 0.0, ITV)
    !     or  new_temp = dT_dEnth * (enth - enth_fp) ! + tfi==0.
    ! but it avoids serious roundoff issues later on when b is large.
    new_temp = dT_dEnth * ((dtEU * (htg - fb*0.0) + m_lay * (enth_in-enth_fp)) / &
                           (m_lay  + dtEU*(fb*dT_dEnth)))
  elseif (ITV%Cp_ice == ITV%Cp_brine) then
    En_J = enth_in  / ITV%enth_unit - ITV%enth_liq_0
    ! Solve a quadratic equation for the new layer temperature, tn:
    !
    !   m * (En_J - (ITV%Cp_Water-Cp_Ice)*T_fr + L + htg*dt/m) = 
    !        (m*Cp_Ice + b*dt) *tn + m*LI*tfi/tn
    !
    AA = m_lay *ITV%Cp_Ice + fb*dtt
    BB = -(m_lay*((En_J - (ITV%Cp_Water-ITV%Cp_Ice)*tfi) + ITV%LI) + htg*dtt)
    CC = m_lay *ITV%LI*tfi
    ! This form avoids round-off errors.
    if (BB >= 0) then
      new_temp = -(BB + sqrt(BB*BB - 4*AA*CC)) / (2*AA)
    else
      new_temp = (2*CC) / (-BB + sqrt(BB*BB - 4*AA*CC))
    endif
!  These should be equivalent.
!   enth = enth_in + (dtEU/m) * (f - b*new_temp)
    enth = enth_from_TS(new_temp, sice, ITV)
  else
    call SIS_error(FATAL, "Write update_lay_enth for Cp_ice /= Cp_brine.")
  endif


!    Figure out whether it is more accurate to use the explicit form for the
!  fluxes or to invert enthalpy conservation for the fluxes.
  if (abs(hf_err_rat*dftop_dT) <= m_lay) then
    ! This branch currently applies in all interior ice layers.
    ftop = ftop_in + dftop_dT*new_temp
    ! An explicit estimate (or no update) is used for the top flux.
    if (hf_err_rat*dfbot_dT <= m_lay) then  ! Use the explicit expression for fbot.
      fbot = fbot_in + dfbot_dT*new_temp
    else ! Use conservation to invert for fbot.
      fbot = (ht_body + ftop) - (m_lay*(enth - enth_in) + extra_enth)/dtEU
    endif
  elseif (hf_err_rat*dfbot_dT <= m_lay) then
    ! Use the explicit expression for fbot and invert for ftop.
    fbot = fbot_in + dfbot_dT*new_temp
    ftop = (fbot - ht_body) + (m_lay*(enth - enth_in) + extra_enth)/dtEU
  else
    ! Conservation is used to invert for both fbot and ftop, partitioning
    ! the changes in proportion to their sensitivities.
    !   dflux = (htg - (m*(enth - enth_in) + extra_enth)/dtEU)
    !   dflux = dfbot - dftop ; dftop / dfbot = dftop_dT / dfbot_dT

    if (dfbot_dT - dftop_dT > 0.0) then
      dflux_dtot_dT = (htg - (m_lay*(enth - enth_in) + extra_enth)/dtEU) / &
                      (dfbot_dT - dftop_dT)
    else
      dflux_dtot_dT = 0.0 ! This should never occur.
    endif

    ftop = ftop_in + dftop_dT * dflux_dtot_dT
    fbot = fbot_in + dfbot_dT * dflux_dtot_dT
  endif

  if (present(temp_new)) temp_new = new_temp

end subroutine update_lay_enth

subroutine ice_check(ms, mi, enthalpy, s_ice, NkIce, msg_part, ITV, &
                      bmelt, tmelt, t_sfc)
  real, intent(in) :: ms, mi
  real, dimension(0:NkIce), intent(in) :: enthalpy
  real, dimension(NkIce), intent(in) :: s_ice
  integer, intent(in) :: NkIce
  character(len=*), intent(in) :: msg_part
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real, optional, intent(in) :: bmelt, tmelt, t_sfc

  character(len=300) :: mesg
  character(len=80) :: msg2
  integer :: k, bad

  real, dimension(0:NkIce) :: t_col

  bad = 0
  if (present(t_sfc)) then
    if ((t_sfc > 1.e-14) .or. (t_sfc < -100.0)) bad = bad+1
  endif
  if (ms <0.0 .or. ms > 3.30e5  ) bad = bad+1
  if (mi <0.0 .or. mi > 1.0e6  ) bad = bad+1

  t_col(0) = temp_from_En_S(enthalpy(0), 0.0, ITV)
  call temp_from_Enth_S(enthalpy(1:), s_ice(:), t_col(1:), ITV)
  do k=0,NkIce ; if ((t_col(k) > 0.0) .or. (t_col(k) < -100.0)) bad = bad+1 ; enddo

  if (bad>0) then
    mesg = "BAD ICE "//trim(msg_part)
    write (msg2,'(" ms,mi=",2(ES11.3))') ms*1e-3,mi*1e-3 ; mesg = trim(mesg)//trim(msg2)
    if (present(t_sfc)) then
      write (msg2,'(" t_sfc,tsn,tice=",ES11.3)') t_sfc ; mesg = trim(mesg)//trim(msg2)
    else
      mesg = trim(mesg)//"tsn,tice="
    endif
    do k=0,NkIce ; write (msg2,'(ES11.3)') t_col(k) ; mesg = trim(mesg)//trim(msg2) ; enddo
    if (present(bmelt)) then
      write (msg2,'(" bmelt=",ES11.3)') bmelt ; mesg = trim(mesg)//trim(msg2)
    endif
    if (present(tmelt)) then
      write (msg2,'(" tmelt=",ES11.3)') tmelt ; mesg = trim(mesg)//trim(msg2)
    endif
    call SIS_error(WARNING, mesg, all_print=.true.)
  endif

end subroutine ice_check

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! T_Freeze - Return the freezing temperature as a function of salinity (and    !
!            possibly later pressure).                                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function T_Freeze(S, ITV)
  real, intent(in) :: S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real :: T_Freeze

  T_Freeze = 0.0 + ITV%dTf_dS * S

end function T_Freeze

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! calculate_T_Freeze - Calculate an array of freezing temperatures for an      !
!            an array of salinities (and maybe later pressures).               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine calculate_T_Freeze(S, T_Freeze, ITV)
  real, dimension(:), intent(in) :: S
  real, dimension(:), intent(out) :: T_Freeze
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.

  integer :: k, nk_ice
  nk_ice = size(S)

  do k=1,nk_ice ; T_Freeze(k) = 0.0 + ITV%dTf_dS * S(k) ; enddo

end subroutine calculate_T_Freeze

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! enthalpy_from_TS - Set a column of enthalpies from temperature and salinity. !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine enthalpy_from_TS(T, S, enthalpy, ITV)
  real, dimension(:), intent(in) :: T, S
  real, dimension(:), intent(out) :: enthalpy
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.

  integer :: k, nk_ice
  nk_ice = size(T)
 
  do k=1,nk_ice ; enthalpy(k) = enth_from_TS(T(k), S(k), ITV) ; enddo

end subroutine enthalpy_from_TS

function enth_from_TS(T, S, ITV) result(enthalpy)
  real, intent(in)  :: T, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real :: enthalpy

  real :: T_fr  ! The freezing temperature in deg C.  
  real :: Cp_Ice, Enth_liq_0, LI, enth_unit
  Cp_Ice = ITV%Cp_Ice ; LI = ITV%LI
  Enth_liq_0 = ITV%Enth_liq_0 ; enth_unit = ITV%enth_unit

  T_fr = -ITV%mu_TS*S

  if ((S == 0.0) .and. (T <= 0.0)) then
    ! Note that at the freezing point, fresh water is assumed to be all ice,
    ! due to the degeneracy in inverting temperature for enthalpy.
    enthalpy = enth_unit * ((ENTH_LIQ_0 - LI) + Cp_Ice*T)
  elseif (T >= T_fr) then ! This layer is already melted, so just warm or cool it to 0 C.
    enthalpy = enth_unit * (ENTH_LIQ_0 + ITV%Cp_Water*T)
  elseif (ITV%Cp_Ice == ITV%Cp_Brine) then
    enthalpy = enth_unit * ((ENTH_LIQ_0 - LI * (1.0 - T_fr/T)) + &
                  (Cp_Ice*T + (ITV%Cp_Water-Cp_Ice)*T_fr))
  else
    call SIS_error(FATAL, "Write enth_from_TS for Cp_ice /= Cp_brine.")
  endif

end function enth_from_TS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! enthalpy_liquid_freeze - Return the enthalpy of liquid water at the freezing !
!    point for a given salinity.                                               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function enthalpy_liquid_freeze(S, ITV)
  real, intent(in)  :: S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real :: enthalpy_liquid_freeze

  enthalpy_liquid_freeze = ITV%enth_unit * &
    (ITV%Cp_water*(-ITV%mu_TS*S) + ITV%ENTH_LIQ_0)

end function enthalpy_liquid_freeze

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! enthalpy_liquid - Returns the enthalpy of liquid water at the given          !
!     temperature and salinity, in enth_unit.                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function enthalpy_liquid(T, S, ITV)
  real, intent(in)  :: T, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real :: enthalpy_liquid

  enthalpy_liquid = ITV%enth_unit * (ITV%ENTH_LIQ_0 + ITV%CP_SeaWater*T)

end function enthalpy_liquid

! This returns the enthalpy change associated with melting water of
! a given temperature (T, in C) and salinity (S), in enth_unit.
function enth_melt(T, S, ITV) result (emelt)
  real, intent(in) :: T, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real             :: emelt

  emelt = enthalpy_liquid_freeze(S, ITV) - enth_from_TS(T, S, ITV)
end function enth_melt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Temp_from_Enth_S - Set a column of temperatures from enthalpy and salinity.  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine Temp_from_Enth_S(En, S, Temp, ITV)
  real, dimension(:), intent(in) :: En, S
  real, dimension(:), intent(out) :: Temp
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.

  integer :: k, nk_ice
  
  nk_ice = size(Temp)
 
  do k=1,nk_ice
    Temp(k) = Temp_from_En_S(En(k), S(k), ITV)
  enddo

end subroutine Temp_from_Enth_S

function dTemp_dEnth(En, S, ITV) result(dT_dE)
  real, intent(in)  :: En, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real :: dT_dE  ! Partial derivative of temperature with enthalpy in degC/Enth_unit.

  real :: I_Cp_Ice, BB, I_CpI_Eu, I_CpW_Eu
  real :: I_enth_unit
  real :: Cp_Ice, LI, Mu_TS
  real :: T_fr  ! The freezing temperature in deg C.  
  real :: En_J  ! Enthalpy in Joules with 0 offset.
  Cp_Ice = ITV%Cp_Ice ; LI = ITV%LI ; Mu_TS = ITV%mu_TS

  I_Cp_Ice = 1.0 / Cp_Ice ; I_enth_unit = 1.0 / ITV%enth_unit
  I_CpI_Eu = 1.0 / (Cp_Ice * ITV%enth_unit)
  I_CpW_Eu = 1.0 / (ITV%Cp_Water * ITV%enth_unit)
  ! I_Cp_Water = 1.0 / ITV%CP_Water

  T_fr = -ITV%mu_TS*S

  En_J = En * I_enth_unit - ITV%enth_liq_0
  if (S <= 0.0) then ! There is a step function for fresh water.
    if (En_J >= 0.0) then ; dT_dE = I_CpW_Eu
    elseif (En_J > -LI) then ; dT_dE = 0.0
    else ; dT_dE = I_CpI_Eu ; endif
  elseif (ITV%Cp_Ice == ITV%Cp_Brine) then
    ! This makes the assumption that all water in the ice and snow categories,
    ! both fluid and in pockets, has the same heat capacity.
    if (En_J < T_fr * ITV%Cp_water) then
      BB = 0.5*((En_J - T_fr*(ITV%Cp_water-ITV%Cp_ice)) + LI)
      dT_dE = I_Cp_Ice * 0.5 * (1.0 - BB / sqrt(BB**2 - T_fr*Cp_Ice*LI))
    else  ! This layer is already melted, so just warm it to 0 C.
      dT_dE = I_CpW_Eu
    endif
  else
    call SIS_error(FATAL, "Write dTemp_dEnth for Cp_ice /= Cp_brine.")
  endif

end function dTemp_dEnth

function Temp_from_En_S(En, S, ITV) result(Temp)
  real, intent(in)  :: En, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real :: Temp  ! Temperature in deg C.

  real :: I_Cp_Ice, I_Cp_Water  ! Inverse heat capacities, in kg K J-1.
  real :: BB
  real :: I_enth_unit
  real :: T_fr  ! The freezing temperature in deg C.  
  real :: Cp_Ice, Cp_Water, LI, Mu_TS
  real :: En_J  ! Enthalpy in Joules with 0 offset.
  Cp_Ice = ITV%Cp_Ice ; Cp_water = ITV%Cp_water ; LI = ITV%LI ; Mu_TS = ITV%mu_TS
  
  I_Cp_Ice = 1.0 / Cp_Ice ; I_enth_unit = 1.0 / ITV%enth_unit
  I_Cp_Water = 1.0 / CP_Water

  T_fr = -ITV%mu_TS*S
  En_J = En * I_enth_unit - ITV%enth_liq_0

  if (S <= 0.0) then ! There is a step function for fresh water.
    if (En_J >= 0.0) then ; Temp = En_J * I_Cp_Water
    elseif (En_J >= -LI) then ; Temp = 0.0
    else ; Temp = I_Cp_Ice * (En_J + LI) ; endif
  elseif (Cp_Ice == ITV%Cp_brine) then
    ! This makes the assumption that all water in the ice and snow categories,
    ! both fluid and in pockets, has the same heat capacity.

    !   LI * (T_fr/T) + Cp_Ice*T = En_J - (ITV%Cp_Water-Cp_Ice)*T_fr + LI 
    !   LI * (T_fr/T) + Cp_Ice*T = 2.0*BB
    !   Cp_Ice*T**2 - 2.0*BB*T + LI * T_fr = 0.0

    if (En_J < T_fr*Cp_Water) then
      BB = 0.5*((En_J - T_fr*(ITV%Cp_water-ITV%Cp_ice)) + LI)
      Temp = I_Cp_Ice * (BB - sqrt(BB**2 - T_fr*Cp_Ice*LI))
    else  ! This layer is completely melted.
      Temp = En_J * I_Cp_Water
    endif
  else
    call SIS_error(FATAL, "Write Temp_from_En_S for Cp_ice /= Cp_brine.")
  endif

end function Temp_from_En_S

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! e_to_melt_TS - return the energy needed to melt a given snow/ice             !
!      configuration, in J kg-1.                                               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function e_to_melt_TS(T, S, ITV) result(e_to_melt)
  real, intent(in) :: T, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.

  real :: e_to_melt  ! The energy required to melt this mixture of ice and brine
                     ! and warm it to its bulk freezing temperature, in J kg-1.

  real :: T_fr  ! The freezing temperature in deg C.  
  T_fr = -ITV%mu_TS*S

  if (T >= T_Fr) then ! This layer is already melted and has excess heat.
    e_to_melt = ITV%Cp_Water * (-T + T_Fr)
  elseif (ITV%Cp_Ice == ITV%Cp_brine) then
    e_to_melt = (ITV%LI - ITV%Cp_ice*T) * (1.0 - T_Fr/T)
  else
    call SIS_error(FATAL, "Write e_to_melt_TS for Cp_ice /= Cp_brine.")
  endif

end function e_to_melt_TS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! energy_melt_enthS - return the energy needed to melt a given snow/ice        !
!      configuration, in J kg-1.                                               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function energy_melt_enthS(En, S, ITV) result(e_to_melt)
  real, intent(in) :: En, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.

  real :: e_to_melt  ! The energy required to melt this mixture of ice and brine
                     ! and warm it to its bulk freezing temperature, in J kg-1.

  e_to_melt = ITV%enth_unit * (enthalpy_liquid_freeze(S, ITV) - En)

end function energy_melt_enthS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_thermo_coefs - return various thermodynamic coefficients.                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine get_SIS2_thermo_coefs(ITV, ice_salinity, Cp_Ice, enthalpy_units, &
                                 specified_thermo_salinity)
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real, dimension(:), optional, intent(out) :: ice_salinity
  real, optional, intent(out) :: Cp_Ice, enthalpy_units
  logical, optional, intent(out) :: specified_thermo_salinity
! Arguments: ITV - The ice_thermo_type that contains all sea-ice thermodynamic
!                  parameters.
!            ice_salinity - The specified salinity of each layer when the
!                           thermodynamic salinities are pre-specified.
!            enthalpy_units - A unit conversion factor for ethalpy from Joules.
!            Cp_Ice - The heat capacity of ice in J kg-1 K-1.
!            specified_thermo_salinity - If true, all thermodynamic calculations
!                  are done with a specified salinity profile that may be
!                  independent of the ice bulk salinity.

  call get_thermo_coefs(ice_salinity=ice_salinity)
  
  if (present(Cp_Ice)) Cp_Ice = ITV%Cp_Ice
  if (present(enthalpy_units)) enthalpy_units = ITV%enth_unit
  if (present(specified_thermo_salinity)) specified_thermo_salinity = .true.

end subroutine get_SIS2_thermo_coefs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_resize_SIS2 - An n-layer code for applying snow and ice thickness and    !
!    temperature changes due to thermodynamic forcing.                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_resize_SIS2(mH_snow, mH_ice, H_to_kg_m2, Enthalpy, Sice_therm, &
                           Salin, snow, frazil, evap, tmlt, bmlt, tfw, NkIce, &
                           heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                           snow_to_ice, salt_to_ice, ITV, CS, bablt, &
                           enthalpy_evap, enthalpy_melt, enthalpy_freeze)
  real, intent(inout) :: mH_snow     ! snow mass per unit area in units of H (often kg m-2).
  real, intent(inout) :: mH_ice      ! ice mass per unit area in units of H (often kg m-2).
  real, intent(in)    :: H_to_kg_m2  ! The conversion factor from the thickness
                                     ! units to kg m-2.
  real, dimension(0:NkIce+1), &
        intent(inout) :: Enthalpy    ! Snow, ice, and ocean enthalpy by layer in enth_units
                                     ! (which might be J/kg).
  real, dimension(NkIce), &
        intent(in)    :: Sice_therm  ! ice salinity by layer, as used for thermodynamics (g/kg)
  real, dimension(NkIce+1), &
        intent(inout) :: Salin       ! Conserved ice bulk salinity by layer (g/kg)
  real, intent(in   ) :: snow        ! new snow (kg/m^2-snow)
  real, intent(in   ) :: frazil      ! frazil in energy units
  real, intent(in   ) :: evap        ! ice evaporation (kg/m^2)
  real, intent(in   ) :: tmlt        ! top melting energy (J/m^2)
  real, intent(in   ) :: bmlt        ! bottom melting energy (J/m^2)
  real, intent(in   ) :: tfw         ! seawater freezing temperature (deg-C)
  integer, intent(in) :: NkIce       ! The number of ice layers.
  real, intent(  out) :: heat_to_ocn ! energy left after ice all melted (J/m^2)
  real, intent(  out) :: h2o_ice_to_ocn ! liquid water flux to ocean (kg/m^2)
  real, intent(  out) :: h2o_ocn_to_ice ! liquid water flux from ocean (kg/m^2)
  real, intent(  out) :: evap_from_ocn! evaporation flux from ocean (kg/m^2)
  real, intent(  out) :: snow_to_ice ! snow below waterline becomes ice
  real, intent(  out) :: salt_to_ice ! Net flux of salt to the ice, in g m-2.
  type(SIS2_ice_thm_CS), intent(in) :: CS
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.

  real, intent(  out), optional :: bablt ! bottom ablation (kg/m^2)
  real, intent(  out), optional :: enthalpy_evap ! The enthalpy loss due to the
                                     ! mass loss by evaporation / sublimation.
  real, intent(  out), optional :: enthalpy_melt ! The enthalpy loss due to the
                                     ! mass loss by melting, in J m-2.
  real, intent(  out), optional :: enthalpy_freeze ! The enthalpy gain due to the
                                     ! mass gain by freezing, in J m-2.

  real :: top_melt, bot_melt, melt_left ! Heating amounts, all in melt_unit.
  real, dimension(0:NkIce) :: m_lay ! temporary ice masses in kg m-2.
  real :: mtot_ice    ! The summed ice mass in kg m-2.
  real :: enth_frazil ! The enthalpy of newly formed frazil ice, in enth_unit.
  real :: enth_freeze ! The enthalpy of newly formed congelation ice, in enth_unit.
  real, dimension(0:NkIce) :: enth_fr ! The snow and ice layers' freezing point
                                      ! enthalpy, in units of enth_unit.
  real, dimension(NkIce) :: mlay_new, enth_ice_new, sal_ice_new
  real :: frazil_per_layer    ! The frazil heat sink from each of the sublayers of
                              ! of the ice, in units of enth_unit.
  real :: t_frazil  ! The temperature which with the frazil-ice is created, in C.
  real :: m_frazil  ! The newly-formed mass per unit area of frazil ice, in kg m-2.
  real :: hw                  ! waterline height above ice base.
  real :: ablation  ! The mass loss from bottom melt, in kg m-2.
  real :: min_dEnth_freeze    ! The minimum enthalpy change that must occur when
                              ! freezing water, usually enough to account for
                              ! the latent heat of fusion in a small fraction of
                              ! the water, in Enth_unit kg-1 (perhaps J kg-1).
  real :: m_k1_to_k2, m_ice_avg
  real :: m_freeze            ! The newly formed ice from freezing, in kg m-2.
  real :: M_melt              ! The ice mass lost to melting, in kg m-2.
  real :: evap_left           ! The remaining evaporation, in kg m-2.
  real :: evap_here           ! The evaporation from the current layer, in kg m-2.
  real :: m_submerged         ! The submerged mass of ice, in kg m-2.
  real :: salin_freeze        ! The salinity of newly frozen ice, in g kg-1.
  real :: enthM_evap, enthM_melt, enthM_freezing, enthM_snowfall
  real :: etot
  real :: enth_unit, LI
  real :: h2o_to_ocn, h2o_orig, h2o_imb
  integer :: k, k1, k2, kold, knew

  enth_unit = ITV%enth_unit ; LI = ITV%LI
  min_dEnth_freeze = (LI*enth_unit) * (1.0-CS%liq_lim)

  top_melt = tmlt*enth_unit ; bot_melt = bmlt*enth_unit

  ! set mass mark; will subtract mass at end for melt flux to ocean
  h2o_orig = H_to_kg_m2 * (mH_snow + mH_ice)

  enth_fr(0) = enthalpy_liquid_freeze(0.0, ITV)
  m_lay(0) = mH_snow * H_to_kg_m2
  do k=1,NkIce ! break out individual layers
    enth_fr(k) = enthalpy_liquid_freeze(sice_therm(k), ITV)
    m_lay(k) = mH_ice * H_to_kg_m2 / NkIce
  enddo

  heat_to_ocn = 0.0   ! for excess melt energy

  evap_from_ocn = 0.0 ! for excess evap-melt
  h2o_ocn_to_ice = 0.0 ; h2o_ice_to_ocn = 0.0 ; snow_to_ice  = 0.0
  salt_to_ice = 0.0
  enthM_freezing = 0.0 ; enthM_melt = 0.0 ; enthM_evap = 0.0 ; enthM_snowfall = 0.0

  if (mH_ice == 0.0) m_lay(0) = 0.0  ! This should already be true! Trap an error?  Convert the snow to ice?

  m_lay(0) = m_lay(0) + snow ! add snow
  enthM_snowfall = snow*enthalpy(0)
  if (evap < 0.0) then
    m_lay(0) = m_lay(0) - evap ! Treat frost formation like snow.
    enthM_snowfall = enthM_snowfall - evap*enthalpy(0)
  endif
  ! Assume that Salin(NkIce+1) already is the freezing salinity.
  salin_freeze = Salin(NkIce+1)

  ! add frazil
  ! Frazil mostly forms in leads, so add its heat uniformly over all of the
  ! layers rather than just adding it to the ice bottom.
  if (frazil > 0.0) then
    frazil_per_layer = (enth_unit*frazil)/NkIce
    do k=1,NkIce
      t_frazil = min(tfw, -ITV%mu_TS*sice_therm(k) - CS%Frazil_temp_offset)
      enth_frazil = min(enth_from_TS(t_frazil, sice_therm(k), ITV), &
                        enthalpy(NkIce+1) - min_dEnth_freeze)
      m_frazil = frazil_per_layer / (enthalpy(NkIce+1) - enth_frazil)

      Enthalpy(k) = (m_lay(k)*Enthalpy(k) + m_frazil*enth_frazil) / &
                    (m_lay(k) + m_frazil)
      Salin(k) = (m_lay(k)*Salin(k) + m_frazil*salin_freeze) / &
                 (m_lay(k) + m_frazil)
      Salt_to_ice = Salt_to_ice + m_frazil*salin_freeze

      m_lay(k) = m_lay(k) + m_frazil
      h2o_ocn_to_ice = h2o_ocn_to_ice + m_frazil

      enthM_freezing = enthM_freezing + m_frazil*enthalpy(NkIce+1)
    enddo
  endif

  if (top_melt < 0.0) then  ! this usually shouldn't happen
    bot_melt = bot_melt + top_melt
    top_melt = 0.0
  endif

  if (bot_melt < 0.0 ) then ! add freezing to bottom layer at tice and salin_freeze.
    ! Enth_freeze is based on the colder of the properties of the existing ice
    ! or the heat content of ocean water after a small fraction has frozen.
    enth_freeze = min(Enthalpy(NkIce), enthalpy(NkIce+1) - min_dEnth_freeze)
    m_freeze = -bot_melt / (Enthalpy(NkIce+1) - enth_freeze)

    Enthalpy(NkIce) = (m_lay(NkIce)*Enthalpy(NkIce) + m_freeze*enth_freeze) / &
                      (m_lay(NkIce) + m_freeze)
    Salin(NkIce) = (m_lay(NkIce)*Salin(NkIce) + m_freeze*salin_freeze) / &
                   (m_lay(NkIce) + m_freeze)
    Salt_to_ice = Salt_to_ice + m_freeze*salin_freeze

    m_lay(NkIce) = m_lay(NkIce) + m_freeze
    h2o_ocn_to_ice = h2o_ocn_to_ice + m_freeze
    enthM_freezing = enthM_freezing + m_freeze*enthalpy(NkIce+1)

    bot_melt = 0.0
  endif

  ! Apply mass losses from evaporation.

  if (evap > 0.0) then ! apply evaporation mass flux
    evap_left = evap
    do k=0,NkIce
      evap_here = min(evap_left, m_lay(k))
      evap_left = evap_left - evap_here
      m_lay(k) = m_lay(k) - evap_here
      ! Assume that evaporation does not make ice salty?
      if (k>0) Salt_to_ice = Salt_to_ice - Salin(k) * evap_here
      enthM_evap = enthM_evap + evap_here * enthalpy(k)
      
      if (evap_left <= 0.0) exit
    enddo

    evap_from_ocn = evap_left
    !   The energy required to evaporate was already taken into account in
    ! ice_thm, but there is excess energy that has not been used here that needs
    ! to be passed on to the ocean.
    heat_to_ocn = heat_to_ocn + evap_left*(hlv+hlf)
  endif

  if (top_melt > 0.0 ) then ! apply top melt heat flux
    melt_left = top_melt
    do k=0,NkIce
      if (melt_left < m_lay(k) * (enth_fr(k) - Enthalpy(k))) then
        M_melt = melt_left / (enth_fr(k) - Enthalpy(k)) ! melt part of this layer
        melt_left = 0.0
      else
        M_melt = m_lay(k) ! melt this whole layer
        melt_left = melt_left - M_melt*(enth_fr(k) - Enthalpy(k))
      endif
      m_lay(k) = m_lay(k) - M_melt
      if (k>0) Salt_to_ice = Salt_to_ice - Salin(k) * M_melt
      h2o_ice_to_ocn = h2o_ice_to_ocn + M_melt
      enthM_melt = enthM_melt + M_melt*enth_fr(k)

      if (melt_left <= 0.0) exit ! All melt energy has been used.
    enddo

    heat_to_ocn = heat_to_ocn + melt_left/enth_unit ! melt heat left after snow & ice gone
  endif

  ! apply bottom melt heat flux

  melt_left = bot_melt ; ablation = 0.0
  if (melt_left > 0.0) then ! melt ice and snow from below
    do k=NkIce,0,-1
      if (melt_left < m_lay(k) * (enth_fr(k) - Enthalpy(k))) then
        M_melt = melt_left / (enth_fr(k) - Enthalpy(k)) ! melt part of this layer
        melt_left = 0.0
      else
        M_melt = m_lay(k) ! melt this whole layer
        melt_left = melt_left - M_melt*(enth_fr(k) - Enthalpy(k))
      endif
      m_lay(k) = m_lay(k) - M_melt
      if (k>0) Salt_to_ice = Salt_to_ice - Salin(k) * M_melt
      h2o_ice_to_ocn = h2o_ice_to_ocn + M_melt
      enthM_melt = enthM_melt + M_melt*enth_fr(k)
      ablation = ablation + M_melt

      if (melt_left <= 0.0) exit ! All melt energy has been used.
    enddo

    heat_to_ocn = heat_to_ocn + melt_left/enth_unit
  endif

  if (present(bablt)) bablt = ablation

  ! There are no further heat or mass losses or gains by the ice+snow.
  if (present(Enthalpy_evap)) Enthalpy_evap = enthM_evap
  if (present(Enthalpy_melt)) Enthalpy_melt = enthM_melt
  if (present(Enthalpy_freeze)) Enthalpy_freeze = enthM_freezing

  ! Make the snow below waterline adjustment.
  mtot_ice = 0.0 ; do k=1,NkIce ; mtot_ice = mtot_ice + m_lay(k) ; enddo

  m_submerged = (mtot_ice+m_lay(0))* (ITV%Rho_ice/ITV%Rho_water) ! The mass of ice that will
                ! be submerged when floating according to Archimede's principle.
  if (m_submerged > mtot_ice) then ! convert snow to ice to maintain ice top at waterline
    snow_to_ice = m_submerged - mtot_ice ! need this much ice mass from snow

    m_lay(0) = m_lay(0) - snow_to_ice

    ! Add ice to the topmost layer and dilute its salinity.
    Enthalpy(1) = (m_lay(1)*Enthalpy(1) + snow_to_ice*Enthalpy(0)) / &
                  (m_lay(1) + snow_to_ice)
    Salin(1) = Salin(1) * m_lay(1) / (m_lay(1) + snow_to_ice)
    ! ### THIS NEEDS TO WORK WITH THE TRACER REGISTRY TO ADD SNOW TRACERS TO
    ! ### THE CORRESPONDING ICE TRACER.

    m_lay(1) = m_lay(1) + snow_to_ice
  else
    snow_to_ice = 0.0
  endif

  ! Even up ice layer thicknesses.
  ! ### THIS NEEDS TO WORK WITH THE TRACER REGISTRY SO THAT IT WORKS ON ALL
  ! ### TRACERS WITH MORE THAN 1 LAYER.
  mtot_ice = 0.0 ; do k=1,NkIce ; mtot_ice = mtot_ice + m_lay(k) ; enddo
  if (mtot_ice == 0.0) then ! There is no ice, so quit.
    do k=1,NkIce ; Enthalpy(k) = enth_fr(k) ; enddo
  else
    do k=1,NkIce ; mlay_new(k) = 0.0 ; enth_ice_new(k) = 0.0 ; enddo
    do k=1,NkIce ; sal_ice_new(k) = 0.0 ; enddo
    m_ice_avg = mtot_ice / NkIce
    k1 = 1 ; k2 = 1
    do  ! Add ice from k1 to k2 to even up layer thicknesses.
      ! m_k1_to_k2 = min(m_ice_avg - mlay_new(k2), m_lay(k1))
      if ((mlay_new(k2) >= m_ice_avg) .and. (k2 < NkIce)) then ; k2 = k2+1
      elseif (m_lay(k1) <= 0.0) then ; k1 = k1+1
      elseif ((m_ice_avg - mlay_new(k2) > m_lay(k1)) .or. (k2 == NkIce)) then
        ! Move all remaining ice from k1 to k2.
        m_k1_to_k2 = m_lay(k1)
        enth_ice_new(k2) = enth_ice_new(k2) + m_k1_to_k2 * Enthalpy(k1)
        sal_ice_new(k2) = sal_ice_new(k2) + m_k1_to_k2 * Salin(k1)
        mlay_new(k2) = mlay_new(k2) + m_k1_to_k2
        m_lay(k1) = 0.0 ! = m_lay(k1) - m_k1_to_k2
        k1 = k1+1
      else
        ! Move some of the ice from k1 to k2.
        m_k1_to_k2 = m_ice_avg - mlay_new(k2)
        enth_ice_new(k2) = enth_ice_new(k2) + m_k1_to_k2 * Enthalpy(k1)
        sal_ice_new(k2) = sal_ice_new(k2) + m_k1_to_k2 * Salin(k1)
        mlay_new(k2) = m_ice_avg ! = mlay_new(k2) + m_k1_to_k2
        m_lay(k1) = m_lay(k1) - m_k1_to_k2
        k2 = k2+1
      endif
      if (k1 > NkIce) exit
    enddo
    do k=1,NkIce ; if (mlay_new(k) > 0.0) then
      Enthalpy(k) = enth_ice_new(k)/mlay_new(k)
      Salin(k) = sal_ice_new(k)/mlay_new(k)
    endif ; enddo
  endif

  mH_ice = mtot_ice / H_to_kg_m2
  mH_snow = m_lay(0) / H_to_kg_m2

  h2o_to_ocn = h2o_orig + snow - (evap-evap_from_ocn) - H_to_kg_m2 * (mH_snow + mH_ice)

  h2o_imb = h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)
  if (abs(h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)) > &
      max(1e-10, 1e-12*h2o_orig, 1e-12*(abs(h2o_ice_to_ocn)+abs(h2o_ocn_to_ice)))) then
    h2o_imb = h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)
  endif

  call ice_check(mH_snow*H_to_kg_m2, mH_ice*H_to_kg_m2, enthalpy, Sice_therm, &
                    NkIce, "at end of ice_resize_SIS2", ITV, &
                    bmelt=bot_melt/enth_unit, tmelt=top_melt/enth_unit)

end subroutine ice_resize_SIS2

subroutine SIS2_ice_thm_end(CS, ITV)
  type(SIS2_ice_thm_CS), pointer :: CS
  type(ice_thermo_type), pointer :: ITV ! A pointer to the ice thermodynamic parameter structure.

  deallocate(ITV)
  deallocate(CS)
end subroutine SIS2_ice_thm_end

end module SIS2_ice_thm
