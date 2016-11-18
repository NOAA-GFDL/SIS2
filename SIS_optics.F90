module SIS_optics

! for calling delta-Eddington shortwave from ice_optics
use ice_shortwave_dEdd, only : shortwave_dEdd0_set_snow, shortwave_dEdd0_set_pond
use ice_shortwave_dEdd, only : shortwave_dEdd0, shortwave_dEdd0_set_params
use ice_shortwave_dEdd, only : dbl_kind, int_kind, nilyr, nslyr
use SIS2_ice_thm, only : ice_thermo_type, get_SIS2_thermo_coefs, T_freeze
! use MOM_EOS, only : EOS_type, EOS_init, EOS_end
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,  only : get_param, log_param, read_param, log_version, param_file_type

implicit none ; private

public :: ice_optics_SIS2, SIS_optics_init, SIS_optics_end

type, public :: SIS_optics_CS ; private

  ! albedos are from CSIM4 assumming 0.53 visible and 0.47 near-ir insolation
  real :: alb_snow        ! albedo of snow (not melting)
  real :: alb_ice         ! albedo of ice (not melting)
  real :: pen_ice         ! ice surface penetrating solar fraction
  real :: opt_dep_ice     ! ice optical depth (m)
  real :: t_range_melt    ! melt albedos scaled in below melting T

  logical :: do_deltaEdd = .true.  ! If true, use a delta-Eddington radiative
                          ! transfer calculation for the shortwave radiation
                          ! within the sea-ice and snow.

  logical :: do_pond = .false. ! activate melt pond scheme - mw/new
  real :: max_pond_frac = 0.5  ! pond water beyond this is dumped
  real :: min_pond_frac = 0.2  ! ponds below sea level don't drain

end type SIS_optics_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_optics_init initalizes the SIS2 optics module
subroutine SIS_optics_init(param_file, CS)

  type(param_file_type), intent(in)    :: param_file
  type(SIS_optics_CS), pointer :: CS

  !
  ! Albedo tuning parameters are documented in:
  !
  ! Briegleb, B.P., and B. Light, 2007:  A delta-Eddington multiple scattering
  !  parameterization for solar radiation in the sea ice component of the
  !  Community Climate System Model, NCAR/TN+472+STR.
  !
  real :: deltaEdd_R_ice  ! delta-Eddington ice albedo tuning, non-dim.
  real :: deltaEdd_R_snow ! delta-Eddington snow albedo tuning, non-dim.
  real :: deltaEdd_R_pond ! delta-Eddington pond albedo tuning, non-dim.
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_optics" ! This module's name.

  if (.not.associated(CS)) allocate(CS)

  call log_version(param_file, mod, version, &
     "This module calculates the albedo and absorption profiles for shortwave radiation.")

  call get_param(param_file, mod, "DO_DELTA_EDDINGTON_SW", CS%do_deltaEdd, &
                 "If true, a delta-Eddington radiative transfer calculation \n"//&
                 "for the shortwave radiation within the sea-ice.", default=.true.)
  call get_param(param_file, mod, "DO_POND", CS%do_pond, &
                 "If true, calculate melt ponds and use them for\n"//&
                 "shortwave radiation calculation.", default=.false.)
  call get_param(param_file, mod, "MIN_POND_FRAC", CS%min_pond_frac, &
                 "Minimum melt pond cover (by ponds at sea level)\n"//&
                 "pond drains to this when ice is porous.", default=0.2)
  call get_param(param_file, mod, "MAX_POND_FRAC", CS%max_pond_frac, &
                 "Maximum melt pond cover - associated with pond volume\n"//&
                 "that suppresses ice top to waterline", default=0.5)

  if (CS%do_deltaEdd) then
    call get_param(param_file, mod, "ICE_DELTA_EDD_R_ICE", deltaEdd_R_ice, &
                   "A dreadfully documented tuning parameter for the radiative \n"//&
                   "propeties of sea ice with the delta-Eddington radiative \n"//&
                   "transfer calculation.", units="nondimensional", default=0.0)
    call get_param(param_file, mod, "ICE_DELTA_EDD_R_SNOW", deltaEdd_R_snow, &
                   "A dreadfully documented tuning parameter for the radiative \n"//&
                   "propeties of snow on sea ice with the delta-Eddington \n"//&
                   "radiative transfer calculation.", &
                   units="nondimensional", default=0.0)
    call get_param(param_file, mod, "ICE_DELTA_EDD_R_POND", deltaEdd_R_pond, &
                   "A dreadfully documented tuning parameter for the radiative \n"//&
                   "propeties of meltwater ponds on sea ice with the delta-Eddington \n"//&
                   "radiative transfer calculation.", units="nondimensional", &
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

end subroutine SIS_optics_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_optics - set albedo, penetrating solar, and ice/snow transmissivity      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_optics_SIS2(mp, hs, hi, ts, tfw, NkIce, alb_vis_dir, alb_vis_dif, &
                    alb_nir_dir, alb_nir_dif, abs_sfc, abs_snow, abs_ice_lay, &
                    abs_ocn, abs_int, CS, ITV, coszen_in)
  real, intent(in   ) :: mp  ! pond mass (kg/m2) mw/new
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
  type(SIS_optics_CS), intent(in) :: CS
  type(ice_thermo_type), intent(in) :: ITV ! The ice thermodynamic parameter structure.
  real, intent(in),optional :: coszen_in

  real :: alb, as, ai, snow_cover, fh
  real :: coalb, I_coalb  ! The coalbedo and its reciprocal.
  real :: SW_frac_top     ! The fraction of the SW at the top of the snow that
                          ! is still present at the top of each ice layer (ND).
  real :: opt_decay_lay   ! The optical extinction in each ice layer (ND).
  real :: rho_ice  ! The nominal density of sea ice in kg m-3.
  real :: rho_snow ! The nominal density of snow in kg m-3.
  real :: rho_water ! The nominal density of sea water in kg m-3.
  real :: pen      ! frac sw passed below the surface (frac 1-pen absorbed at the surface)
  real :: sal_ice_top(1)  ! A specified surface salinity of ice.
  real :: temp_ice_freeze ! The freezing temperature of the top ice layer, in C.
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
    hprad      ! pond depth (m) for radiation code - may be diagnosed

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

  real (kind=dbl_kind) :: max_mp, hs_mask_pond, pond_decr

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

    if ( CS%do_pond ) then ! mw/new
      call get_SIS2_thermo_coefs(ITV, rho_ice=rho_ice, rho_snow=rho_snow, &
                                 rho_water=rho_water)

      max_mp = (Rho_water-Rho_ice)*hi  ! max pond allowed by waterline
      fp(1,1) = CS%max_pond_frac*sqrt(min(1.0,mp/max_mp))
      ! set average pond depth (max. = 2*average)
      hprad(1,1) = mp/(fp(1,1)*1000)  ! freshwater density = 1000 kg/m2
      fs(1,1) = fs(1,1)*(1-fp(1,1))   ! reduce fs to frac of pond-free ice
      ! decrement fp (increment fs) for snow masking of pond: pond is completely
      ! masked when snow depth contains 2*average_pond_depth in its pore space
      if (hs>0.0 .and. hprad(1,1)>0.0) then
        hs_mask_pond = 2*hprad(1,1)*Rho_ice/(Rho_ice-Rho_snow)
        pond_decr = fp(1,1)*min(1.0,hs/hs_mask_pond)
        fp(1,1) = fp(1,1) - pond_decr
        fs(1,1) = fs(1,1) + pond_decr
      endif
    else
      call shortwave_dEdd0_set_pond(nx_block, ny_block, icells, indxi, indxj, &
               aice, Tsfc, fs, fp, hprad) ! out: fp, hprad
    endif

    call shortwave_dEdd0  (nx_block, ny_block, icells, indxi, indxj, coszen, &
             aice, vice, vsno, fs, rhosnw, rsnw, fp, hprad, swvdr, swvdf, &
             swidr, swidf, alvdr, alvdf, alidr, alidf, fswsfc, fswint, &
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
    call get_SIS2_thermo_coefs(ITV, ice_salinity=sal_ice_top)
    temp_ice_freeze = T_freeze(sal_ice_top(1), ITV)

    fh = min(atan(5.0*hi)/atan(5.0*0.5),1.0) ! use this form from CSIM4 to
    ! reduce albedo for thin ice
    if (ts+CS%T_RANGE_MELT > temp_ice_freeze) then        ! reduce albedo for melting as in
       ! CSIM4 assuming 0.53/0.47 vis/ir
       as = as-0.1235*min((ts+CS%T_RANGE_MELT- temp_ice_freeze)/CS%T_RANGE_MELT,1.0)
       ai = ai-0.075 *min((ts+CS%T_RANGE_MELT- temp_ice_freeze)/CS%T_RANGE_MELT,1.0)
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


subroutine SIS_optics_end(CS)
  type(SIS_optics_CS), pointer :: CS

  deallocate(CS)

end subroutine SIS_optics_end

end module SIS_optics
