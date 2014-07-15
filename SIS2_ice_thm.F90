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

type, public :: ice_thermo_type ; private
  real :: Cp_ice            ! The heat capacity of ice, in J kg-1 K-1.
  real :: Cp_Water = 4.2e3  ! The heat capacity of liquid seawater 4200 J/(kg K)
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

  real :: h_lo_lim        ! hi/hs lower limit for temp. calc.
  real :: frazil_temp_offset = 0.5 ! A temperature offset between the 

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

  use ice_thm_mod, only : TFI
  use constants_mod, only : hlf ! latent heat of fusion - 334e3 J/(kg-ice)

  type(param_file_type),       intent(in)    :: param_file
  type(SIS2_ice_thm_CS), pointer :: CS
  type(ice_thermo_type), pointer :: ITV ! A pointer to the ice thermodynamic parameter structure.

  real :: deltaEdd_R_ice  ! Mysterious delta-Eddington tuning parameters, unknown.
  real :: deltaEdd_R_snow ! Mysterious delta-Eddington tuning parameters, unknown.
  real :: deltaEdd_R_pond ! Mysterious delta-Eddington tuning parameters, unknown.
  character(len=40)  :: mod = "SIS2_ice_thm" ! This module's name.

  if (.not.associated(CS)) allocate(CS)
  if (.not.associated(ITV)) allocate(ITV)

  CS%temp_ice_freeze = TFI

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
  call get_param(param_file, mod, "CP_WATER", ITV%Cp_water, &
                 "The heat capacity of sea water, approximated as a \n"//&
                 "constant.", units="J kg-1 K-1", default=4200.0)
  call get_param(param_file, mod, "CP_ICE", ITV%Cp_ice, &
                 "The heat capacity of fresh ice, approximated as a \n"//&
                 "constant.", units="J kg-1 K-1", default=2100.0)
  call get_param(param_file, mod, "DTFREEZE_DS", ITV%dTf_dS, &
                 "The derivative of the freezing temperature with salinity.", &
                 units="deg C PSU-1", default=-0.054)
  ITV%mu_TS = -ITV%dTf_dS

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
! ice5lay_temp - A subroutine that calculates the snow and ice enthalpy        !
!    changes due to surface forcing.                                           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_temp_SIS2(hsno, tsn, hice, tice, sice, sh_T0, B, sol, tfw, fb, &
                         tsurf, dtt, NkIce, tmelt, bmelt, CS, ITV, check_conserve)

  real, intent(in   ) :: hsno  ! snow thickness (m)
  real, intent(inout) :: tsn   ! snow temperature (deg-C)
  real, intent(in   ) :: hice  ! ice thickness (m)
  real, dimension(NkIce), &
        intent(inout) :: tice ! ice temperature by layer (deg-C)
  real, dimension(NkIce), &
        intent(in)    :: Sice ! ice salinity by layer (g/kg)
  real, intent(in   ) :: sh_T0 ! net surface heat flux (+ up) at ts=0 (W/m^2)
  real, intent(in   ) :: B     ! d(sfc heat flux)/d(ts) [W/(m^2 deg-C)]
  real, dimension(0:NkIce), &
        intent(in)    :: sol   ! Solar heating of the snow and ice layers (W m-2)
  real, intent(in   ) :: tfw   ! seawater freezing temperature (deg-C)
  real, intent(in   ) :: fb    ! heat flux from ocean to ice bottom (W/m^2)
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
  real, dimension(0:NkIce) :: temp_est
  real, dimension(NkIce) :: tfi, tice_est ! estimated new ice temperatures
  real :: m_ice, m_snow, e_extra
  real, dimension(NkIce) :: m_lay
  real :: kk, k10, k0a, k0skin
  real :: I_bb, b_denom_1
  real :: comp_rat ! The complement of rat, going from 0 to 1.
  real :: tsf, k0a_x_ta, tsno_est, hie, salt_part, rat, tsurf_est
  real, dimension(0:NkIce+1) :: cc ! Interfacial coupling coefficients.
  real, dimension(0:NkIce) :: bb   ! Effective layer heat capacities.
  real, dimension(0:NkIce) :: cc_bb ! Remaining coupling ratios.
  real :: I_liq_lim     ! The inverse of CS%liq_lim.
  real :: col_enth1, col_enth2, col_enth3
  real :: d_e_extra, e_extra_sum
  real :: tflux_bot, tflux_bot_diff, tflux_sfc, sum_sol, d_tflux_bot
  real :: hsnow_eff ! , Ks_h
  logical :: col_check
  integer :: k

  col_check = .false. ; if (present(check_conserve)) col_check = check_conserve

  A = -sh_T0

  m_ice = ITV%Rho_ice*hice/NkIce ! ice mass of each layer
  m_snow = ITV%Rho_snow*hsno     ! full snow layer mass
  call calculate_T_Freeze(sice, tfi, ITV)    ! freezing temperature of ice layers
  hie = max(hice, CS%H_LO_LIM); ! prevent thin ice inaccuracy (mw)
  kk = NkIce*CS%KI/hie                     ! full ice layer conductivity

  tsf = tfi(1)                          ! surface freezing temperature
  if (hsno>0.0) tsf = 0.0
  hsnow_eff = hsno + max(1e-35, 1e-20*CS%H_LO_LIM)
  
  k10 = 2.0*(CS%KS*(NkIce*CS%KI)) / (hice*CS%KS + hsnow_eff*(NkIce*CS%KI)) ! coupling ice layer 1 to snow
  k0a = (CS%KS*B) / (0.5*B*hsnow_eff + CS%KS)      ! coupling snow to "air"
  k0skin = 2.0*CS%KS / hsnow_eff
  k0a_x_ta = (CS%KS*A) / (0.5*B*hsnow_eff + CS%KS) ! coupling times "air" temperture

  ! Determine the enthalpy for conservation checks.
  if (col_check) then
    do k=1,NkIce ; m_lay(k) = m_ice ; enddo
    col_enth1 = column_enthalpy(m_snow, m_lay, tsn, tice, sice, ITV)
  endif

  ! First get non-conservative estimate with implicit treatment of layer coupling.
  
  ! This is the start of what was ice_temp_est.

  ! Determine the effective layer heat capacities.
  !   bb = dheat/dTemp.  It should be a proper linearization of the enthalpy equation.
  bb(0) = m_snow*ITV%Cp_ice
  do k=1,NkIce   ! load bb with heat capacity term.
    salt_part = 0.0
    if (sice(k)>0.0) salt_part = tfi(k)*ITV%LI/(tice(k)*tice(k))
    bb(k) = m_ice*(ITV%Cp_ice-salt_part) ! add coupling to this later
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
  temp_est(NkIce) = ( (sol(NkIce)*dtt + bb(NkIce)*tice(NkIce)) + &
                      cc(NkIce+1)*tfw ) * I_bb
  comp_rat = b_denom_1 * I_bb
  cc_bb(NkIce) = cc(NkIce) * I_bb

  do k=NkIce-1,1,-1
    b_denom_1 = bb(k) + comp_rat*cc(k+1)
    I_bb =  1.0 / (b_denom_1 + cc(k))
    temp_est(k) = ((sol(k)*dtt + bb(k)*tice(k)) + cc(k+1)*temp_est(k+1)) * I_bb

    comp_rat = b_denom_1 * I_bb ! 1.0 >= comp_rat >= 0.0
    cc_bb(k) = cc(k) * I_bb
  end do

  b_denom_1 = bb(0) + comp_rat*cc(1)
  I_bb =  1.0 / (b_denom_1 + cc(0))
  !   This is a complete calculation of temp_est(0), assuming that the surface
  ! flux is given by A - B*tsurf, with tsurf estimated as part of this calculation.
  temp_est(0) = ((bb(0)*tsn + k0a_x_ta*dtt) + cc(1)*temp_est(1)) * I_bb
  ! Note that the SIS1 code omits sol(0) at its equivalent of this point.
  temp_est(0) = temp_est(0) + (sol(0)*dtt) * I_bb
  
  ! Diagnose the surface skin temperature by matching the diffusive fluxes in
  ! the snow with the atmospheric fluxes.  I.e. solve the following for tsurf_est:
  !  (A - B*tsurf_est) = k0skin * (tsurf_est - temp_est(0))
  tsurf_est = (A + k0skin*temp_est(0)) / (B + k0skin)

  if (tsurf_est > tsf) then
    ! The surface is melting, set tsurf to melt temp. and recalculate I_bb.
    tsurf_est = tsf
    ! cc(0) = k0skin*dtt

    I_bb =  1.0 / (b_denom_1 + (k0skin*dtt))
    temp_est(0) = ((bb(0)*tsn + k0skin*dtt*tsf) + cc(1)*temp_est(1)) * I_bb
    temp_est(0) = temp_est(0) + (sol(0)*dtt) * I_bb
  endif
  ! Go back DOWN the ice column to get the final temperatures.
  do k=1,NkIce
    ! Probably this should be limited so that temp_est(k) <= t_freeze(k)
    temp_est(k) = temp_est(k) + cc_bb(k) * temp_est(k-1)
  end do


  tice_est(:) = temp_est(1:NkIce)
  tsno_est = temp_est(0)

! The following line disables the implicit coupling estimate; in a 10 year CORE
! test the skin temperature without implicit estimate was mostly within .05K,
! but was cooler around topography in the NH (up to 1K near Baffin Island).
! The thickness looked very similar.
!   tice_est(:) = tice(:) ; tsno_est=tsn

  !
  ! conservative pass going UP the ice column
  !
  tice_est(NkIce) = laytemp_SIS2(m_ice, tfi(NkIce), sol(NkIce) + kk*(2*tfw+tice_est(NkIce-1)), &
                                 3*kk, tice(NkIce), dtt, ITV)
  do k=NkIce-1,2,-1
    tice_est(k) = laytemp_SIS2(m_ice, tfi(k), sol(k) + kk*(tice_est(k-1)+tice_est(k+1)), &
                               2*kk, tice(k), dtt, ITV)
  enddo
  tice_est(1) = laytemp_SIS2(m_ice, tfi(1), sol(1) + (kk*tice_est(2) + k10*tsno_est), &
                             kk + k10, tice(1), dtt, ITV)
  tsno_est = laytemp_SIS2(m_snow, 0.0, sol(0) + (k10*tice_est(1)+k0a_x_ta), &
                          k10+k0a, tsn, dtt, ITV)
  tsurf = (A + k0skin*tsno_est) / (B + k0skin)  ! diagnose surface skin temp.

  if (tsurf > tsf) then ! surface is melting, redo with surf. at melt temp.
    tsurf = tsf
    tsno_est = laytemp_SIS2(m_snow, 0.0, sol(0) + k10*tice_est(1) + k0skin*tsf,&
                            k10+k0skin, tsn, dtt, ITV)
    ! add in surf. melt
    if (hsno>0.0) then
      tmelt = tmelt+dtt*((A-B*tsurf)-2*CS%KS*(tsurf-tsno_est)/hsno)
      tflux_sfc = dtt*2*CS%KS*(tsurf-tsno_est)/hsno
    else
      tmelt = tmelt + dtt*((sol(0)+A-B*tsurf) - k10*(tsurf-tice_est(1))) ! tsno = tsurf
      tflux_sfc = dtt*(k10*(tsurf-tice_est(1)) - sol(0))
    endif
  else
    tflux_sfc = dtt*(A - B*tsurf)
  endif
  tsn = tsno_est ! finalize snow temperature

  !
  ! conservative pass going DOWN the ice column
  !
  tice(1)  = laytemp_SIS2(m_ice, tfi(1), sol(1)+kk*tice_est(1+1) &
                          + k10*(tsno_est-tice_est(1)), kk, tice(1), dtt, ITV)
  do k=2,NkIce-1 ! flux from above is fixed, only have downward feedback
    tice(k) = laytemp_SIS2(m_ice, tfi(k), sol(k)+kk*tice_est(k+1) &
                                       +kk*(tice(k-1)-tice_est(k)), kk, tice(k), dtt, ITV)
  enddo
  tice(NkIce) = laytemp_SIS2(m_ice, tfi(NkIce), sol(NkIce)+2*kk*tfw &
                             +kk*(tice(NkIce-1)-tice_est(NkIce)), 2*kk, tice(NkIce), dtt, ITV)
  !
  ! END conservative update
  !
  if (col_check) then
    col_enth2 = column_enthalpy(m_snow, m_lay, tsn, tice, sice, ITV)
    sum_sol = 0.0 ; do k=0,NkIce ; sum_sol = sum_sol + dtt*sol(k) ; enddo 
    tflux_bot = (col_enth2 - col_enth1) - (sum_sol + tflux_sfc)
    tflux_bot_diff = 2*kk*(tfw-tice(NkIce))*dtt

    d_tflux_bot = tflux_bot_diff - tflux_bot
    if (abs(d_tflux_bot) > 1.0e-9*(abs(tflux_bot) + abs(tflux_bot_diff) + &
                                   abs(col_enth2 - col_enth1))) then
      d_tflux_bot = tflux_bot_diff - tflux_bot
    endif
  endif

  ! The following is a dangerous calculation if kk is too large, in part because
  ! all of the other fluxes are calculated implicitly for the layer above, so
  ! the temperatures are well bounded.  Can this be rearranged or calculated
  ! as a residual of the heat changes in the ice and snow?
  bmelt = bmelt + dtt*(fb - 2*kk*(tfw-tice(NkIce))) ! add in bottom melting/freezing

  e_extra_sum = 0.0
  if (tsn > 0.0) then ! put excess snow energy into top melt
    e_extra = ITV%Cp_ice*tsn * m_snow
    tmelt = tmelt + e_extra
    e_extra_sum = e_extra_sum + e_extra
    tsn = 0.0
  endif

  I_liq_lim = 1.0 / CS%liq_lim
  do k=1,NkIce ; if (tice(k)>tfi(k)*I_liq_lim) then ! push excess energy to closer of top or bottom melt
    e_extra = (enth_from_TS(tice(k),sice(k), ITV) - &
               enth_from_TS(tfi(k)*I_liq_lim,sice(k), ITV)) * (m_ice / ITV%enth_unit)
    e_extra_sum = e_extra_sum + e_extra
    tice(k) = tfi(k)*I_liq_lim
    if (k<=NkIce/2) then
      tmelt = tmelt+e_extra
    else
      bmelt = bmelt+e_extra
    endif
  endif ; enddo
  
  if (col_check) then
    col_enth3 = column_enthalpy(m_snow, m_lay, tsn, tice, sice, ITV)

    d_e_extra = col_enth3 - (col_enth2 - e_extra_sum)
    if (abs(d_e_extra) > 1.0e-12*(abs(col_enth3) + abs(col_enth2) + &
                                   abs(e_extra_sum))) then
      d_e_extra = (col_enth3 - col_enth2) - e_extra_sum
    endif
  endif

  call temp_check(tsurf, hsno, tsn, hice, tice, NkIce, bmelt, tmelt)

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
    ! relationship by extending the linear expression for frozen water.
    !
    !   m * {Cp_Ice} * (tn-tp) = dtt * (f - b*tn)
    !
    new_temp = (m*Cp_Ice*tp + f*dtt) / (m*Cp_Ice + b*dtt) ! = -BB/AA
    ! 
    ! if (tp > 0.0) then
    !   E0 = Cp_Ice*tp + LI  ! >= LI
    ! else
    !   E0 = Cp_Ice*tp  ! < 0
    ! endif
    ! ! Determine whether the new solution will be above, at, or below freezing.
    ! if (m*E0 + dtt * (f - b*tfi) >= m*LI) then
    !   new_temp = (m*(E0 - LI) + f*dtt) / (m*Cp_Ice + b*dtt)
    !   extra_heat = LI
    ! elseif (m*E0 + dtt * (f - b*tfi) >= 0) then
    !   new_temp = 0.0 ; extra_heat = m*E0 + ddt * (f - b*tfi)
    ! else
    !   new_temp = (m*E0 + f*dtt) / (m*Cp_Ice + b*dtt)
    !   extra_heat = 0.0
    ! endif
  else
    if (tp >= tfi) then
      E0 = Cp_Ice*(tp - tfi)  ! >= 0
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
  endif

end function laytemp_SIS2

subroutine temp_check(ts, hs, tsn, hi, t_ice, NkIce, bmelt, tmelt)
  real, intent(in) :: ts, hs, tsn, hi, bmelt, tmelt
  real, dimension(NkIce), intent(in) :: t_ice
  integer, intent(in) :: NkIce
  integer :: k, bad

  bad = 0
  if (ts >0.0.or.ts <-100.0) bad = bad+1
  if (tsn>0.0.or.tsn<-100.0) bad = bad+1
  do k=1,NkIce ; if (t_ice(k) >0.0 .or. t_ice(k) < -100.0) bad = bad+1 ; enddo

  if (bad>0) then
    print *, 'BAD ICE AFTER TEMP ', 'hs/hi=',hs,hi,'ts/tsn/tice=',ts, &
                      tsn,t_ice(:),'tmelt/bmelt=',tmelt,bmelt
  endif
end subroutine temp_check

subroutine resize_check(ms, mi, enthalpy, s_ice, NkIce, bmelt, tmelt, ITV)
  real, intent(in) :: ms, mi, bmelt, tmelt
  real, dimension(0:NkIce), intent(in) :: enthalpy
  real, dimension(NkIce), intent(in) :: s_ice
  integer, intent(in) :: NkIce
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  integer :: k, bad

  real :: t_snow
  real, dimension(NkIce) :: t_ice

  bad = 0
  if (ms <0.0 .or. ms > 3.30e5  ) bad = bad+1
  if (mi <0.0 .or. mi > 1.0e6  ) bad = bad+1
  t_snow = temp_from_En_S(enthalpy(0), 0.0, ITV)

  if (t_snow>0.0 .or. t_snow<-100.0) bad = bad+1
  do k=1,NkIce
    t_ice(k) = temp_from_En_S(enthalpy(k), s_ice(k), ITV)
    if (t_ice(k) >0.0 .or. t_ice(k) < -100.0) bad = bad+1
  enddo

  if (bad>0) then
    print *, 'BAD ICE AFTER RESIZE ', 'hs/hi=',ms/330.,mi/905.,'tsn/tice=',&
                      t_snow, t_ice(:),'tmelt/bmelt=',tmelt,bmelt
  endif
end subroutine resize_check

subroutine unpack_check(hs, tsn, hi, t_ice, NkIce, cn)
  real, intent(in) :: hs, tsn, hi, cn
  real, dimension(NkIce), intent(in) :: t_ice
  integer, intent(in) :: NkIce
  integer :: k, bad

  bad = 0
  if (hs <0.0.or.hs > 1e3  ) bad = bad+1
  if (hi <0.0.or.hi > 1e3  ) bad = bad+1
  if (tsn>0.0.or.tsn<-100.0) bad = bad+1
  do k=1,NkIce ; if (t_ice(k) >0.0 .or. t_ice(k) < -100.0) bad = bad+1 ; enddo

  if (bad>0) print *, 'BAD ICE AFTER UNPACK ', 'hs/hi=',hs,hi,'tsn/tice=', &
                      tsn,t_ice(:),'cn=',cn
end subroutine unpack_check

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
  
  real :: Mu_TS, Cp_Ice, Enth_liq_0, LI, enth_unit
  Mu_TS = ITV%mu_TS ; Cp_Ice = ITV%Cp_Ice ; LI = ITV%LI
  Enth_liq_0 = ITV%Enth_liq_0 ; enth_unit = ITV%enth_unit

  ! This makes the assumption that all water in the ice and snow categories,
  ! both fluid and in pockets, has the same heat capacity.
  if ((S == 0.0) .and. (T <= 0.0)) then
    enthalpy = enth_unit * ((ENTH_LIQ_0 - LI) + Cp_Ice*T)
  elseif (T <= -MU_TS*S) then
    enthalpy = enth_unit * ((ENTH_LIQ_0 - LI * (1.0 + MU_TS*S/T)) + Cp_Ice*T)
       !### + enth_unit * (-MU_TS*S) * (CP_Water - Cp_Ice)
  else  ! This layer is already melted, so just warm it to 0 C.
    enthalpy = enth_unit * (ENTH_LIQ_0 + Cp_Ice*T)  !### Change Cp_Ice to CP_Water.
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
    ((-ITV%Cp_Ice*ITV%mu_TS)*S + ITV%ENTH_LIQ_0)  !### Change Cp_Ice to CP_Water.

end function enthalpy_liquid_freeze

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! enthalpy_liquid - Returns the enthalpy of liquid water at the given          !
!     temperature and salinity, in enth_unit.                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function enthalpy_liquid(T, S, ITV)
  real, intent(in)  :: T, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real :: enthalpy_liquid

  enthalpy_liquid = ITV%enth_unit * (ITV%ENTH_LIQ_0 + ITV%CP_Water*T)

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
  real :: I_Cp_Ice, BB
  real :: I_enth_unit
  real :: En_J  ! Enthalpy in Joules with 0 offset.
  
  nk_ice = size(Temp)
!  I_Cp_Ice = 1.0 / ITV%Cp_Ice
!  I_enth_unit = 1.0 / ITV%enth_unit
 
  do k=1,nk_ice
    Temp(k) = Temp_from_En_S(En(k), S(k), ITV)
!   En_J = En(k) * I_enth_unit - ENTH_LIQ_0
!   ! This makes the assumption that all water in the ice and snow categories,
!   ! both fluid and in pockets, has the same heat capacity.
!   if (S(k) <= 0.0) then ! There is a step function for fresh water.
!     if (En_J >= 0.0) then ; Temp(k) = En_J * I_Cp_Ice
!     elseif (En_J >= -ITV%LI) then ; Temp(k) = 0.0
!     else ; Temp(k) = I_Cp_Ice * (En_J + ITV%LI) ; endif
!   else
!     if (En_J < -ITV%MU_TS*S(k)*ITV%Cp_Ice) then
!       BB = 0.5*(En_J + ITV%LI)
!       Temp(k) = I_Cp_Ice * (BB - sqrt(BB**2 + ITV%MU_TS*S(k)*ITV%Cp_Ice*ITV%LI))
!     else  ! This layer is already melted, so just warm it to 0 C.
!       Temp(k) = En_J * I_Cp_Ice
!     endif
!   endif
  enddo

end subroutine Temp_from_Enth_S

function Temp_from_En_S(En, S, ITV) result(Temp)
  real, intent(in)  :: En, S
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real :: Temp  ! Temperature in deg C.

  real :: I_Cp_Ice, BB
  real :: I_enth_unit
  real :: Cp_Ice, LI, Mu_TS
  real :: En_J  ! Enthalpy in Joules with 0 offset.
  Cp_Ice = ITV%Cp_Ice ; LI = ITV%LI ; Mu_TS = ITV%mu_TS
  
  I_Cp_Ice = 1.0 / Cp_Ice ; I_enth_unit = 1.0 / ITV%enth_unit
  ! I_Cp_Water = 1.0 / CP_Water

  En_J = En * I_enth_unit - ITV%enth_liq_0
  ! This makes the assumption that all water in the ice and snow categories,
  ! both fluid and in pockets, has the same heat capacity.
  if (S <= 0.0) then ! There is a step function for fresh water.
    if (En_J >= 0.0) then ; Temp = En_J * I_Cp_Ice  ! ### Change to I_Cp_Water
    elseif (En_J >= -LI) then ; Temp = 0.0
    else ; Temp = I_Cp_Ice * (En_J + LI) ; endif
  else
    if (En_J < -MU_TS*S*Cp_Ice) then  ! ### Change Cp_Ice to Cp_Water
      BB = 0.5*(En_J + LI)        ! ### Change En_J to En_J - (-MU_TS*S) * (CP_Water - Cp_Ice)
      Temp = I_Cp_Ice * (BB - sqrt(BB**2 + MU_TS*S*Cp_Ice*LI))
    else  ! This layer is already melted, so just warm it to 0 C.
      Temp = En_J * I_Cp_Ice          ! ### Change to I_Cp_Water
    endif
  endif

end function Temp_from_En_S

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
subroutine ice_resize_SIS2(m_snow, m_ice, Enthalpy, Sice_therm, Salin, snow, &
                           frazil, evap, tmlt, bmlt, tfw, NkIce, heat_to_ocn, &
                           h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                           snow_to_ice, salt_to_ice, ITV, CS, bablt, &
                           enthalpy_evap, enthalpy_melt, enthalpy_freeze)
  real, intent(inout) :: m_snow      ! snow mass per unit area (kg m-2)
  real, intent(inout) :: m_ice       ! ice mass per unit area (kg m-2)
  real, dimension(0:NkIce+1), &
        intent(inout) :: Enthalpy    ! snow, ice, and ocean enthalpy by layer (J/kg)
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
  real, dimension(0:NkIce) :: m_lay ! temporary ice mass
  real :: enth_frazil ! The enthalpy of newly formed frazil ice, in enth_unit..
  real, dimension(0:NkIce) :: enth_fr ! The snow and ice layers' freezing point
                                      ! enthalpy, in units of enth_unit.
  real, dimension(NkIce) :: mlay_new, enth_ice_new, sal_ice_new
  real :: frazil_per_layer    ! The frazil heat sink from each of the sublayers of
                              ! of the ice, in units of enth_unit.
  real :: t_frazil  ! The temperature which with the frazil-ice is created, in C.
  real :: m_frazil  ! The newly-formed mass per unit area of frazil ice, in kg m-2.
  real :: hw                  ! waterline height above ice base.
  real :: ablation  ! The mass loss from bottom melt, in kg m-2.
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

  top_melt = tmlt*enth_unit ; bot_melt = bmlt*enth_unit

  ! set mass mark; will subtract mass at end for melt flux to ocean
  h2o_orig = m_snow + m_ice

  enth_fr(0) = enthalpy_liquid_freeze(0.0, ITV)
  m_lay(0) = m_snow
  do k=1,NkIce ! break out individual layers
    enth_fr(k) = enthalpy_liquid_freeze(sice_therm(k), ITV)
    m_lay(k) = m_ice / NkIce
  enddo

  heat_to_ocn = 0.0   ! for excess melt energy

  evap_from_ocn = 0.0 ! for excess evap-melt
  h2o_ocn_to_ice = 0.0 ; h2o_ice_to_ocn = 0.0 ; snow_to_ice  = 0.0
  salt_to_ice = 0.0
  enthM_freezing = 0.0 ; enthM_melt = 0.0 ; enthM_evap = 0.0 ; enthM_snowfall = 0.0

  if (m_ice == 0.0) m_lay(0) = 0.0  ! This should already be true! Trap an error?  Convert the snow to ice?

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
      enth_frazil = enth_from_TS(t_frazil, sice_therm(k), ITV)
      ! Enth_fr here should be based on the temperature of ocean water and the
      ! salinity of the newly formed ice.
      m_frazil = frazil_per_layer / (enth_fr(k) - enth_frazil)

      Enthalpy(k) = (m_lay(k)*Enthalpy(k) + m_frazil*enth_frazil) / &
                    (m_lay(k) + m_frazil)
      Salin(k) = (m_lay(k)*Salin(k) + m_frazil*salin_freeze) / &
                 (m_lay(k) + m_frazil)
      Salt_to_ice = Salt_to_ice + m_frazil*salin_freeze

      m_lay(k) = m_lay(k) + m_frazil
      h2o_ocn_to_ice = h2o_ocn_to_ice + m_frazil

      ! This should be based on the enthalpy of ocean water.
      enthM_freezing = enthM_freezing + m_frazil*enth_fr(k)
    enddo
  endif

  if (top_melt < 0.0) then  ! this usually shouldn't happen
    bot_melt = bot_melt + top_melt
    top_melt = 0.0
  endif

  if (bot_melt < 0.0 ) then ! add freezing to bottom layer at tice and salin_freeze.
    ! Enth_fr here should be based on the temperature of ocean water and the
    ! salinity of the newly formed ice.
    m_freeze = -bot_melt / (enth_fr(NkIce) - Enthalpy(NkIce))

    Salin(NkIce) = (m_lay(NkIce)*Salin(NkIce) + m_freeze*salin_freeze) / &
                   (m_lay(NkIce) + m_freeze)
    Salt_to_ice = Salt_to_ice + m_freeze*salin_freeze

    m_lay(NkIce) = m_lay(NkIce) + m_freeze
    h2o_ocn_to_ice = h2o_ocn_to_ice + m_freeze
    enthM_freezing = enthM_freezing + m_freeze*enth_fr(NkIce)

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
  m_ice = 0.0 ; do k=1,NkIce ; m_ice = m_ice + m_lay(k) ; enddo

  m_submerged = (m_ice+m_lay(0))* (ITV%Rho_ice/ITV%Rho_water) ! The mass of ice that will
                ! be submerged when floating according to Archimede's principle.
  if (m_submerged > m_ice) then ! convert snow to ice to maintain ice top at waterline
    snow_to_ice = m_submerged - m_ice ! need this much ice mass from snow

    m_lay(0) = m_lay(0) - snow_to_ice

    ! Add ice to the topmost layer.  Currently, the salinity of the target layer
    ! and the bulk ice salinity do not change.
    Enthalpy(1) = (m_lay(1)*Enthalpy(1) + snow_to_ice*Enthalpy(0)) / &
                  (m_lay(1) + snow_to_ice)

    ! Pick one of the following...
    ! Salt_to_ice = Salt_to_ice + Salin(k) * snow_to_ice
    Salin(1) = Salin(1) * m_lay(1) / (m_lay(1) + snow_to_ice)

    m_lay(1) = m_lay(1) + snow_to_ice
  else
    snow_to_ice = 0.0;
  endif

  ! Even up ice layer thicknesses.
  m_ice = 0.0 ; do k=1,NkIce ; m_ice = m_ice + m_lay(k) ; enddo
  if (m_ice == 0.0) then ! There is no ice, so quit.
    do k=1,NkIce ; Enthalpy(k) = enth_fr(k) ; enddo
  else
    do k=1,NkIce ; mlay_new(k) = 0.0 ; enth_ice_new(k) = 0.0 ; enddo
    do k=1,NkIce ; sal_ice_new(k) = 0.0 ; enddo
    m_ice_avg = m_ice / NkIce
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

   m_snow = m_lay(0)

  h2o_to_ocn = h2o_orig + snow - (evap-evap_from_ocn) - (m_snow + m_ice)

  h2o_imb = h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)
  if (abs(h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)) > &
      max(1e-10, 1e-12*h2o_orig, 1e-12*(abs(h2o_ice_to_ocn)+abs(h2o_ocn_to_ice)))) then
    h2o_imb = h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)
  endif

  call resize_check(m_snow, m_ice, enthalpy, Sice_therm, NkIce, bot_melt/enth_unit, &
                    top_melt/enth_unit, ITV)

end subroutine ice_resize_SIS2

function column_enthalpy(m_snow, m_ice, t_snow, t_ice, s_ice, ITV)
  real, dimension(:), intent(in) :: m_ice, T_ice, S_ice
  real, intent(in) :: m_snow, t_snow
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.

  real :: column_enthalpy

  integer :: k, nk_ice
  real :: I_nk_ice
  
  nk_ice = size(T_ice)

  column_enthalpy = m_snow*enth_from_TS(t_snow, 0.0, ITV)
  do k=1,nk_ice
    column_enthalpy = column_enthalpy + m_ice(k)*enth_from_TS(t_ice(k), s_ice(k), ITV)
  enddo

end function column_enthalpy

subroutine SIS2_ice_thm_end(CS, ITV)
  type(SIS2_ice_thm_CS), pointer :: CS
  type(ice_thermo_type), pointer :: ITV ! A pointer to the ice thermodynamic parameter structure.

  deallocate(ITV)
  deallocate(CS)
end subroutine SIS2_ice_thm_end

end module SIS2_ice_thm
