!> Specifies the sea-ice optical properties
module SIS_optics

! This file is a part of SIS2. See LICENSE.md for the license.

! for calling delta-Eddington shortwave from ice_optics
use ice_shortwave_dEdd, only : shortwave_dEdd0_set_snow, shortwave_dEdd0_set_pond
use ice_shortwave_dEdd, only : shortwave_dEdd0, shortwave_dEdd0_set_params
use ice_shortwave_dEdd, only : dbl_kind, int_kind, nilyr, nslyr
use SIS2_ice_thm, only : ice_thermo_type, get_SIS2_thermo_coefs, T_freeze
! use MOM_EOS, only : EOS_type, EOS_init, EOS_end
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,  only : get_param, log_param, read_param, log_version, param_file_type
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

public :: ice_optics_SIS2, SIS_optics_init, SIS_optics_end, bright_ice_temp
public :: VIS_DIR, VIS_DIF, NIR_DIR, NIR_DIF

! parameters facilitate the use of 4-D arrays for shortwave radiation and
! albedos within SIS2.
integer, parameter :: VIS_DIR=1 !< Indicates the visible direct band
integer, parameter :: VIS_DIF=2 !< Indicates the visible diffuse band
integer, parameter :: NIR_DIR=3 !< Indicates the near-infrared direct band
integer, parameter :: NIR_DIF=4 !< Indicates the near-infrared diffuse band

!> This type contains the parameters regulating sea-ice optics.
type, public :: SIS_optics_CS ; private

  ! albedos are from CSIM4 assumming 0.53 visible and 0.47 near-ir insolation
  real :: alb_snow        !< albedo of snow (not melting) [nondim]
  real :: alb_ice         !< albedo of ice (not melting) [nondim]
  real :: pen_ice         !< ice surface penetrating solar fraction [nondim]
  real :: opt_dep_ice     !< ice optical depth [Z ~> m]
  real :: t_range_melt    !< melt albedos scaled in below melting T [degC]

  logical :: do_deltaEdd = .true.  !< If true, use a delta-Eddington radiative
                          !! transfer calculation for the shortwave radiation
                          !! within the sea-ice and snow.

  logical :: do_pond = .false. !< activate melt pond scheme - mw/new
  real :: max_pond_frac = 0.5  !< pond water beyond this is dumped [nondim]
  real :: min_pond_frac = 0.2  !< ponds below sea level don't drain [nondim]

  logical :: slab_optics = .false. !< If true use the very old slab ice optics
                                   !! from the supersource model.
  real :: slab_crit_thick !< The thickness beyond which the slab ice optics no
                          !! longer exhibits a thickness dependencs on albedo [Z ~> m].
  real :: slab_alb_ocean  !< The ocean albedo as used in the slab ice optics [nondim].
  real :: slab_min_ice_alb !< The minimum thick ice albedo with the slab ice optics [nondim].

end type SIS_optics_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_optics_init initalizes the SIS2 optics module
subroutine SIS_optics_init(param_file, US, CS, slab_optics)

  type(param_file_type), intent(in) :: param_file  !< Parameter file handle
  type(unit_scale_type), intent(in) :: US          !< A structure with unit conversion factors
  type(SIS_optics_CS),   pointer    :: CS          !< A pointer to the SIS_optics control structure.
  logical, optional,     intent(in) :: slab_optics !< If true use the very old slab ice optics
                                                   !! from the supersource model.

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

  real :: T_range_dflt, alb_ice_dflt
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "SIS_optics" ! This module's name.

  if (.not.associated(CS)) allocate(CS)

  CS%slab_optics = .false. ; if (present(slab_optics)) CS%slab_optics = slab_optics

  call log_version(param_file, mdl, version, &
     "This module calculates the albedo and absorption profiles for shortwave radiation.")

  if (CS%slab_optics) then
    call log_param(param_file, mdl, "! USE_SLAB_ICE_OPTICS", CS%slab_optics, &
                 "Using the very old slab-style ice optics.")
  endif

  call get_param(param_file, mdl, "DO_DELTA_EDDINGTON_SW", CS%do_deltaEdd, &
                 "If true, a delta-Eddington radiative transfer calculation "//&
                 "for the shortwave radiation within the sea-ice.", &
                 default=.not.CS%slab_optics)
  call get_param(param_file, mdl, "DO_POND", CS%do_pond, &
                 "If true, calculate melt ponds and use them for "//&
                 "shortwave radiation calculation.", default=.false.)
  call get_param(param_file, mdl, "MIN_POND_FRAC", CS%min_pond_frac, &
                 "Minimum melt pond cover (by ponds at sea level) "//&
                 "pond drains to this when ice is porous.", default=0.2, units="nondim")
  call get_param(param_file, mdl, "MAX_POND_FRAC", CS%max_pond_frac, &
                 "Maximum melt pond cover - associated with pond volume "//&
                 "that suppresses ice top to waterline", default=0.5, units="nondim")

  call get_param(param_file, mdl, "ICE_DELTA_EDD_R_ICE", deltaEdd_R_ice, &
                 "A dreadfully documented tuning parameter for the radiative "//&
                 "propeties of sea ice with the delta-Eddington radiative "//&
                 "transfer calculation.", units="nondimensional", default=0.0, &
                 do_not_log=.not.CS%do_deltaEdd)
  call get_param(param_file, mdl, "ICE_DELTA_EDD_R_SNOW", deltaEdd_R_snow, &
                 "A dreadfully documented tuning parameter for the radiative "//&
                 "propeties of snow on sea ice with the delta-Eddington "//&
                 "radiative transfer calculation.", units="nondimensional", &
                 default=0.0, do_not_log=.not.CS%do_deltaEdd)
  call get_param(param_file, mdl, "ICE_DELTA_EDD_R_POND", deltaEdd_R_pond, &
                 "A dreadfully documented tuning parameter for the radiative "//&
                 "propeties of meltwater ponds on sea ice with the delta-Eddington "//&
                 "radiative transfer calculation.", units="nondimensional", &
                 default=0.0, do_not_log=.not.CS%do_deltaEdd)
  call shortwave_dEdd0_set_params(deltaEdd_R_ice, deltaEdd_R_snow, deltaEdd_R_pond)

  call get_param(param_file, mdl, "SNOW_ALBEDO", CS%alb_snow, &
                 "The albedo of dry snow atop sea ice.", units="nondim", &
                 default=0.85, do_not_log=CS%do_deltaEdd)
  alb_ice_dflt = 0.5826 ; if (CS%slab_optics) alb_ice_dflt = 0.8
  call get_param(param_file, mdl, "ICE_ALBEDO", CS%alb_ice, &
                 "The albedo of dry bare sea ice.", units="nondim", &
                 default=alb_ice_dflt, do_not_log=CS%do_deltaEdd)
  call get_param(param_file, mdl, "ICE_SW_PEN_FRAC", CS%pen_ice, &
                 "The fraction of the unreflected shortwave radiation that "//&
                 "penetrates into the ice.", units="Nondimensional", &
                 default=0.3, do_not_log=CS%do_deltaEdd)
  call get_param(param_file, mdl, "ICE_OPTICAL_DEPTH", CS%opt_dep_ice, &
                 "The optical depth of shortwave radiation in sea ice.", &
                 units="m", default=0.67, scale=US%m_to_Z, do_not_log=CS%do_deltaEdd)
  T_range_dflt = 1.0 ; if (CS%slab_optics) T_range_dflt = 10.0
  call get_param(param_file, mdl, "ALBEDO_T_MELT_RANGE", CS%t_range_melt, &
                 "The temperature range below freezing over which the "//&
                 "albedos are changed by partial melting.", units="degC", &
                 default=1.0, do_not_log=CS%do_deltaEdd)

  ! These parameters pertain only to the ancient slab ice optics parameterization.
  call get_param(param_file, mdl, "SLAB_OPTICS_CRITICAL_THICK", CS%slab_crit_thick, &
                 "The thickness beyond which the slab ice optics no longer "//&
                 "exhibits a thickness dependencs on albedo.", units="m", scale=US%m_to_Z, &
                 default=1.0, do_not_log=.not.CS%slab_optics)
  call get_param(param_file, mdl, "SLAB_OPTICS_OCEAN_ALBEDO", CS%slab_alb_ocean, &
                 "The ocean albedo as used in the slab ice optics.", units="nondim", &
                 default=0.1, do_not_log=.not.CS%slab_optics)
  call get_param(param_file, mdl, "SLAB_OPTICS_MIN_ICE_ALBEDO", CS%slab_min_ice_alb, &
                 "The minimum thick ice albedo with the slab ice optics.", &
                 units="nondim", default=0.55, do_not_log=.not.CS%slab_optics)

end subroutine SIS_optics_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_optics_SIS2 sets albedo, penetrating solar, and ice/snow transmissivity
subroutine ice_optics_SIS2(m_pond, m_snow, m_ice, ts, tfw, NkIce, albedos, abs_sfc, &
                    abs_snow, abs_ice_lay, abs_ocn, abs_int, US, CS, ITV, coszen_in)
  real, intent(in   ) :: m_pond   !< pond mass [R Z ~> kg m-2]
  real, intent(in   ) :: m_snow   !< snow mass per unit area [R Z ~> kg m-2]
  real, intent(in   ) :: m_ice    !< ice thickness [R Z ~> kg m-2]
  real, intent(in   ) :: ts       !< surface temperature [degC]
  real, intent(in   ) :: tfw      !< seawater freezing temperature [degC]
  integer, intent(in) :: NkIce    !< The number of sublayers in the ice
  real, dimension(:), intent(  out) :: albedos !< ice surface albedos (0-1) [nondim]
  real, intent(  out) :: abs_sfc  !< fraction of absorbed SW that is absorbed at surface [nondim]
  real, intent(  out) :: abs_snow !< fraction of absorbed SW that is absorbed in snow [nondim]
  real, intent(  out) :: abs_ice_lay(NkIce) !< fraction of absorbed SW that is absorbed by each ice layer [nondim]
  real, intent(  out) :: abs_ocn  !< fraction of absorbed SW that is absorbed in ocean [nondim]
  real, intent(  out) :: abs_int  !< fraction of absorbed SW that is absorbed in ice interior [nondim]
  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors
  type(SIS_optics_CS),   intent(in) :: CS  !< The ice optics control structure.
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real, intent(in), optional :: coszen_in !< The cosine of the solar zenith angle [nondim].

  ! Local variables
  real :: hs              ! snow thickness [Z ~> m]
  real :: hi              ! ice thickness [Z ~> m]
  real :: alb             ! The albedo for all bands, 0-1 [nondim].
  real :: as              ! A snow albedo, 0-1 [nondim].
  real :: ai              ! The ice albedo, 0-1 [nondim].
  real :: snow_cover      ! The fraction of the area covered by snow, 0-1 [nondim].
  real :: fh              ! A weighting fraction of the ice albedo (as compared
                          ! with the albedo of water) when ice is thin, 0-1 [nondim].
  real :: coalb, I_coalb  ! The coalbedo (0-1) and its reciprocal.
  real :: SW_frac_top     ! The fraction of the SW at the top of the snow that
                          ! is still present at the top of each ice layer [nondim].
  real :: opt_decay_lay   ! The optical extinction in each ice layer [nondim].
  real :: rho_ice         ! The nominal density of sea ice [R ~> kg m-3].
  real :: rho_snow        ! The nominal density of snow [R ~> kg m-3].
  real :: rho_water       ! The nominal density of sea water [R ~> kg m-3].
  real :: pen             ! The fraction of the shortwave flux that will pass below
                          ! the surface (frac 1-pen absorbed at the surface) [nondim]
  real :: sal_ice_top(1)  ! A specified surface salinity of ice [gSalt kg-1].
  real :: temp_ice_freeze ! The freezing temperature of the top ice layer [degC].
  real :: max_mp          ! The maximum melt pond mass at the waterline [R Z ~> kg m-2]
  integer :: m, b, nb
  character(len=200) :: mesg

  real :: tcrit, thick_ice_alb  ! Slab optics variables

  integer (kind=int_kind) :: &
    nx_block, ny_block, & ! block dimensions
    icells                ! number of ice-covered grid cells

  integer (kind=int_kind), dimension (1) :: &
    indxi   , & ! compressed indices for ice-covered cells
    indxj

  ! inputs
  real (kind=dbl_kind), dimension (1,1) :: &
    aice   , & ! concentration of ice
    vice   , & ! volume of ice [m]
    vsno   , & ! volume of snow [m]
    Tsfc   , & ! surface temperature
    coszen , & ! cosine of solar zenith angle
    tarea  , & ! cell area - not used
    swvdr  , & ! sw down, visible, direct  [W m-2]
    swvdf  , & ! sw down, visible, diffuse [W m-2]
    swidr  , & ! sw down, near IR, direct  [W m-2]
    swidf      ! sw down, near IR, diffuse [W m-2]

  ! outputs
  real (kind=dbl_kind), dimension (1,1) :: &
    fs     , & ! horizontal coverage of snow
    fp     , & ! pond fractional coverage (0 to 1) [nondim]
    hprad      ! pond depth [m] for radiation code - may be diagnosed

  real (kind=dbl_kind), dimension (1,1,1) :: &
    rhosnw , & ! density in snow layer [kg m-3]
    rsnw       ! grain radius in snow layer [micro-meters]

  real (kind=dbl_kind), dimension (1,1,18) :: &
    trcr        ! aerosol tracers

  real (kind=dbl_kind), dimension (1,1) :: &
    alvdr   , & ! visible, direct, albedo [nondim]
    alvdf   , & ! visible, diffuse, albedo [nondim]
    alidr   , & ! near-ir, direct, albedo [nondim]
    alidf   , & ! near-ir, diffuse, albedo [nondim]
    fswsfc  , & ! SW absorbed at snow/bare ice/pondedi ice surface [W m-2]
    fswint  , & ! SW interior absorption (below surface, above ocean) [W m-2]
    fswthru     ! SW through snow/bare ice/ponded ice into ocean [W m-2]

  real (kind=dbl_kind), dimension (1,1,1) :: &
    Sswabs      ! SW absorbed in snow layer [W m-2]

  real (kind=dbl_kind), dimension (1,1,NkIce) :: &
    Iswabs      ! SW absorbed in ice layer [W m-2]

  real (kind=dbl_kind), dimension (1,1) :: &
    albice  , & ! bare ice albedo, for history [nondim]
    albsno  , & ! snow albedo, for history [nondim]
    albpnd      ! pond albedo, for history [nondim]

  real (kind=dbl_kind) :: hs_mask_pond, pond_decr

  nb = size(albedos)

  call get_SIS2_thermo_coefs(ITV, rho_ice=rho_ice, rho_snow=rho_snow, rho_water=rho_water)
  hi = (1.0 / Rho_ice) * m_ice
  hs = (1.0 / Rho_snow) * m_snow

  if (CS%slab_optics) then
    ! This option uses a very old slab ice albedo parameterization, which was
    ! used in Supersource and other GFDL models from the 1990s and before.
    tcrit = tfw - CS%T_RANGE_MELT
    if (ts <= tcrit) then
      thick_ice_alb = CS%alb_ice
    else if (ts >= tfw) then
      thick_ice_alb = CS%slab_min_ice_alb
    else
      thick_ice_alb = CS%alb_ice + (CS%slab_min_ice_alb-CS%alb_ice) * &
                      (ts-tcrit) / CS%T_RANGE_MELT
    endif

    if (hi >= CS%slab_CRIT_THICK) then
      alb = thick_ice_alb
    else
      alb = CS%slab_alb_ocean + (thick_ice_alb - CS%slab_alb_ocean) * &
                    sqrt(hi/CS%slab_CRIT_THICK)
    endif

    do b=1,nb ; albedos(b) = alb ; enddo
    abs_sfc = 1.0
    abs_snow = 0.0 ; abs_ice_lay(:) = 0.0 ; abs_ocn = 0.0 ; abs_int = 0.0
  elseif (CS%do_deltaEdd) then

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
    if (present(coszen_in)) coszen(1,1) = max(0.01,coszen_in)
    Tsfc(1,1) = ts
    vsno(1,1) = US%Z_to_m*hs
    vice(1,1) = US%Z_to_m*hi
    swvdr(1,1) = 0.25
    swvdf(1,1) = 0.25
    swidr(1,1) = 0.25
    swidf(1,1) = 0.25

    call shortwave_dEdd0_set_snow(nx_block, ny_block, icells, indxi, indxj, &
             aice, vsno, Tsfc, fs, rhosnw, rsnw) ! out: fs, rhosnw, rsnw

    if ( CS%do_pond ) then ! mw/new

      max_mp = (rho_water - rho_ice) * hi  ! max pond allowed by waterline
      fp(1,1) = CS%max_pond_frac*sqrt(min(1.0,m_pond/max_mp))
      ! set average pond depth (max. = 2*average)
      hprad(1,1) = US%RZ_to_kg_m2*m_pond / (fp(1,1)*1000.0) ! freshwater density = 1000 kg/m2
      fs(1,1) = fs(1,1)*(1.0-fp(1,1))   ! reduce fs to frac of pond-free ice
      ! decrement fp (increment fs) for snow masking of pond: pond is completely
      ! masked when snow depth contains 2*average_pond_depth in its pore space
      if (hs>0.0 .and. hprad(1,1)>0.0) then
        hs_mask_pond = 2*hprad(1,1)*rho_ice / (rho_ice - rho_snow)
        pond_decr = fp(1,1)*min(1.0, US%Z_to_m*hs/hs_mask_pond)
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

    albedos(:) = 0.0
    if (nb < max(VIS_DIR, VIS_DIF, NIR_DIR, NIR_DIF)) call SIS_error(FATAL, &
      "The argument array albedos is not large enough in the call to ice_optics_SIS2.")
    albedos(VIS_DIR) = alvdr(1,1) ; albedos(VIS_DIF) = alvdf(1,1)
    albedos(NIR_DIR) = alidr(1,1) ; albedos(NIR_DIF) = alidf(1,1)

    pen = (fswint(1,1) + fswthru(1,1)) * I_coalb
    abs_int = fswint(1,1) * I_coalb

  else
    as = CS%alb_snow ; ai = CS%alb_ice
    snow_cover = hs / (hs + 0.02*US%m_to_Z)  ! thin snow partially covers ice
    call get_SIS2_thermo_coefs(ITV, ice_salinity=sal_ice_top)
    temp_ice_freeze = T_freeze(sal_ice_top(1), ITV)

    fh = min(atan(5.0*US%Z_to_m*hi)/atan(5.0*0.5),1.0) ! use this form from CSIM4 to
    ! reduce albedo for thin ice
    if (ts+CS%T_RANGE_MELT > temp_ice_freeze) then        ! reduce albedo for melting as in
       ! CSIM4 assuming 0.53/0.47 vis/ir
       as = as-0.1235*min((ts+CS%T_RANGE_MELT- temp_ice_freeze)/CS%T_RANGE_MELT,1.0)
       ai = ai-0.075 *min((ts+CS%T_RANGE_MELT- temp_ice_freeze)/CS%T_RANGE_MELT,1.0)
    endif
    ai = fh*ai+(1-fh)*0.06                 ! reduce albedo for thin ice

    alb = snow_cover*as + (1-snow_cover)*ai
    do b=1,nb ; albedos(b) = alb ; enddo

    pen = (1-snow_cover)*CS%pen_ice
    opt_decay_lay = exp(-hi/(NkIce*CS%opt_dep_ice))
    abs_ocn = pen * exp(-hi/(CS%opt_dep_ice))
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

!> bright_ice_temp returns the skin temperature (in degC) below which the snow
!! and ice attain their greatest brightness and albedo no longer varies, for
!! the highest attainable salinity.
function bright_ice_temp(CS, ITV) result(bright_temp)
  type(SIS_optics_CS), intent(in)   :: CS  !< The ice optics control structure
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: bright_temp

  real :: salin_max       ! The maximum attainable salinity [gSalt kg-1].
  real :: temp_freeze_min ! The freezing temperature of water at salin_max [degC].

  salin_max = 40.0

  temp_freeze_min = T_freeze(salin_max, ITV)

  if (CS%do_deltaEdd) then ! This is hard-coded for the delta-Eddington scheme.
    bright_temp = temp_freeze_min - 1.0
  else
    bright_temp = temp_freeze_min - CS%T_RANGE_MELT
  endif

end function bright_ice_temp

!> Deallocate memory associated with the SIS_optics module
subroutine SIS_optics_end(CS)
  type(SIS_optics_CS), pointer :: CS !< The ice optics control structure that is deallocated here

  deallocate(CS)

end subroutine SIS_optics_end

end module SIS_optics
