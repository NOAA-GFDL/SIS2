!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Contains a number of common SIS types, along with subroutines to perform various tasks on
!! these types, including allocation, deallocation, registration for restarts, and checksums.
module SIS_types


use ice_grid,          only : ice_grid_type
use MOM_coms,          only : PE_here
use MOM_domains,       only : MOM_domain_type, pass_vector, BGRID_NE, CGRID_NE, clone_MOM_domain
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_time_manager,  only : time_type, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type
use SIS_diag_mediator, only : SIS_diag_ctrl, post_data=>post_SIS_data
use SIS_diag_mediator, only : register_SIS_diag_field, register_static_field
use SIS_debugging,     only : chksum, Bchksum, Bchksum_pair, hchksum, uvchksum
use SIS_debugging,     only : check_redundant_B, check_redundant_C
use SIS_framework,     only : domain2D, CORNER, EAST_FACE, NORTH_FACE, redistribute_data
use SIS_restart,       only : register_restart_field, SIS_restart_CS, restore_SIS_state
use SIS_restart,       only : query_initialized=>query_inited, only_read_from_restarts
use SIS_framework,     only : safe_alloc, safe_alloc_ptr
use SIS_framework,     only : coupler_1d_bc_type, coupler_2d_bc_type, coupler_3d_bc_type
use SIS_framework,     only : coupler_type_spawn, coupler_type_initialized
use SIS_framework,     only : coupler_type_redistribute_data, coupler_type_copy_data
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_tracer_registry, only : SIS_tracer_registry_type
use SIS2_ice_thm,      only : ice_thermo_type, SIS2_ice_thm_CS, get_SIS2_thermo_coefs
use SIS2_ice_thm,      only : enth_from_TS, energy_melt_EnthS, temp_from_En_S

implicit none ; private

#include <SIS2_memory.h>

public :: ice_state_type, alloc_IST_arrays, ice_state_register_restarts
public :: IST_chksum, IST_bounds_check, copy_IST_to_IST, dealloc_IST_arrays
public :: ice_state_read_alt_restarts, register_fast_to_slow_restarts
public :: rescale_fast_to_slow_restart_fields, rescale_ice_state_restart_fields
public :: ice_ocean_flux_type, alloc_ice_ocean_flux, dealloc_ice_ocean_flux
public :: ocean_sfc_state_type, alloc_ocean_sfc_state, dealloc_ocean_sfc_state
public :: fast_ice_avg_type, alloc_fast_ice_avg, dealloc_fast_ice_avg, copy_FIA_to_FIA
public :: OSS_chksum, IOF_chksum, FIA_chksum, register_unit_conversion_restarts
public :: ice_rad_type, ice_rad_register_restarts, dealloc_ice_rad
public :: simple_OSS_type, alloc_simple_OSS, dealloc_simple_OSS, copy_sOSS_to_sOSS
public :: redistribute_IST_to_IST, redistribute_FIA_to_FIA, redistribute_sOSS_to_sOSS
public :: total_sfc_flux_type, alloc_total_sfc_flux, dealloc_total_sfc_flux
public :: copy_TSF_to_TSF, redistribute_TSF_to_TSF, TSF_chksum
public :: copy_Rad_to_Rad, redistribute_Rad_to_Rad, alloc_ice_Rad
public :: translate_OSS_to_sOSS

integer, parameter :: NBANDS=4 !< the number of 4-D arrays for shortwave radiation and
!! albedos within SIS2.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> This structure contains the ice model state, and is intended to be private
!! to SIS2.  It is not to be shared with other components and modules, and may
!! use different indexing conventions than other components.
type ice_state_type
  ! The 8 of the following 10 variables constitute the sea-ice state.
  real, allocatable, dimension(:,:,:) :: part_size !< The fractional coverage of a grid cell by
                !! each ice thickness category [nondim], 0 to 1.  Category 0 is open ocean.
                !!  The sum of part_size is 1.
  ! These velocities are only used on the slow ice processors
  real, allocatable, dimension(:,:) :: u_ice_B  !< The pseudo-zonal ice velocity along the
                !! along the grid directions on a B-grid [L T-1 ~> m s-1].
                !! All thickness categories are assumed to have the same velocities.
  real, allocatable, dimension(:,:) :: v_ice_B  !< The pseudo-meridional ice velocity along the
                !! along the grid directions on a B-grid [L T-1 ~> m s-1].
  real, allocatable, dimension(:,:) :: u_ice_C  !< The pseudo-zonal ice velocity along the
                !! along the grid directions on a C-grid [L T-1 ~> m s-1].
                !! All thickness categories are assumed to have the same velocities.
  real, allocatable, dimension(:,:) :: v_ice_C  !< The pseudo-meridional ice velocity along the
                !! along the grid directions on a C-grid [L T-1 ~> m s-1].

  real, allocatable, dimension(:,:,:) :: &
    mH_pond, &  !< The mass per unit area of the pond in each category [R Z ~> kg m-2].
    mH_snow, &  !< The mass per unit area of the snow in each category [R Z ~> kg m-2].
    mH_ice, &   !< The mass per unit area of the ice in each category [R Z ~> kg m-2].
    t_surf      !< The surface temperature [Kelvin].

  real, allocatable, dimension(:,:) :: &
    snow_to_ocn, & !< The mass per unit ocean area of snow that will be dumped into the
                   !! ocean due to recent mechanical activities like ridging or drifting [R Z ~> kg m-2].
    water_to_ocn, & !< The mass per unit ocean area of pond water that will be dumped into the
                   !! ocean due to recent mechanical activities like ridging or drifting [R Z ~> kg m-2].
    enth_snow_to_ocn !< The average enthalpy of the snow that will be dumped into the
                   !! ocean due to recent mechanical activities like ridging or drifting [Q ~> J kg-1].

  real, allocatable, dimension(:,:,:,:) :: sal_ice  !< The salinity of the sea ice
                !! in each category and fractional thickness layer [gSalt kg-1].
  real, allocatable, dimension(:,:,:,:) :: enth_ice !< The enthalpy of the sea ice
                !! in each category and fractional thickness layer [Q ~> J kg-1].
  real, allocatable, dimension(:,:,:,:) :: enth_snow !< The enthalpy of the snow
                !! in each category and snow thickness layer [Q ~> J kg-1].
  real, allocatable, dimension(:,:) :: rdg_rate !< The rate of fractional area loss by ridging [T-1 ~> s-1]
  real, allocatable, dimension(:,:,:) :: &
    rdg_mice    !< A diagnostic of the ice load that was formed by ridging [R Z ~> kg m-2].

  logical :: Cgrid_dyn !< If true use a C-grid discretization of the sea-ice dynamics.
  logical :: valid_IST !< If true, this is currently the valid state of the ice.  Otherwise the ice
                       !! is in the midst of a dynamics cycle where the evolving state has changes
                       !! that are not yet reflected here.

  type(SIS_tracer_registry_type), pointer :: TrReg => NULL() !< A pointer to the SIS tracer registry

  type(ice_thermo_type), pointer  :: ITV => NULL() !< A pointer to the ice thermodynamics type
end type ice_state_type

!> ocean_sfc_state_type contains variables that describe the ocean's surface
!! state as seen by the slowly evolving sea-ice, on the ice grid.
type ocean_sfc_state_type
  ! 7 of the following 9 variables describe the ocean state as seen by the sea ice.
  real, allocatable, dimension(:,:) :: &
    s_surf , &  !< The ocean's surface salinity [gSalt kg-1].
    SST_C  , &  !< The ocean's bulk surface temperature [degC].
    T_fr_ocn, & !< The freezing point temperature at the ocean's surface salinity [degC].
    u_ocn_B, &  !< The ocean's zonal velocity on B-grid points [L T-1 ~> m s-1].
    v_ocn_B, &  !< The ocean's meridional velocity on B-grid points [L T-1 ~> m s-1].
    u_ocn_C, &  !< The ocean's zonal velocity on C-grid points [L T-1 ~> m s-1].
    v_ocn_C     !< The ocean's meridional velocity on C-grid points [L T-1 ~> m s-1].
  real, allocatable, dimension(:,:) :: bheat !< The upward diffusive heat flux from the ocean
                !! to the ice at the base of the ice [Q R Z T-1 ~> W m-2].
  real, allocatable, dimension(:,:) :: frazil !< A downward heat flux from the ice into the ocean
                !! associated with the formation of frazil ice in the ocean integrated over a
                !! timestep [Q R Z ~> J m-2]. This is the input value and is not changed by the ice.
  real, allocatable, dimension(:,:) :: sea_lev !< The equivalent sea-level, after any non-levitating
                !! ice has been converted to sea-water, as determined by the ocean [Z ~> m].
                !! Sea-ice only contributes by applying pressure to the ocean that is then
                !! (partially) converted back to its equivalent by the ocean.

  type (coupler_2d_bc_type) :: &
    tr_fields   !< A structure of fields related to properties for additional tracers.

!   type(coupler_3d_bc_type)   :: ocean_fields       ! array of fields used for additional tracers

  real :: kmelt !< A constant that is used in the calculation of the ocean/ice basal heat flux,
                !! [Q R Z T-1 degC-1 ~> W m-2 degC-1].  This could be replaced with an array
                !! reflecting the turbulence in the under-ice ocean boundary layer and the effective
                !! depth of the reported value of t_ocn.

  logical :: Cgrid_dyn !< If true use a C-grid discretization of the sea-ice dynamics.

  !>@{ diagnostic IDs for ocean surface properties
  integer :: id_sst=-1, id_sss=-1, id_ssh=-1, id_uo=-1, id_vo=-1, id_frazil=-1
  !!@}
end type ocean_sfc_state_type

!> simple_OSS_type contains variables that describe the ocean's surface
!! state as seen by the fast sea-ice or atmosphere, on the ice grid.
type simple_OSS_type
  ! The following 5 variables describe the ocean state as seen by the
  ! atmosphere and use for the rapid thermodynamic sea ice changes.
  real, allocatable, dimension(:,:) :: &
    s_surf , &  !< The ocean's surface salinity [gSalt kg-1].
    SST_C  , &  !< The ocean's bulk surface temperature [degC].
    T_fr_ocn, & !< The freezing point temperature at the ocean's surface salinity [degC].
    u_ocn_A, &  !< The ocean's zonal surface velocity on A-grid points [L T-1 ~> m s-1].
    v_ocn_A, &  !< The ocean's meridional surface velocity on A-grid points [L T-1 ~> m s-1].
    u_ice_A, &  !< The sea ice's zonal velocity on A-grid points [L T-1 ~> m s-1].
    v_ice_A     !< The sea ice's meridional velocity on A-grid points [L T-1 ~> m s-1].
  real, allocatable, dimension(:,:) :: bheat !< The upward diffusive heat flux
                !! from the ocean to the ice at the base of the ice [Q R Z T-1 ~> W m-2].

  type (coupler_2d_bc_type) :: &
    tr_fields   !< A structure of fields related to properties for additional tracers.
end type simple_OSS_type


!> fast_ice_avg_type contains variables that describe the fluxes between the
!! atmosphere and the ice or that have been accumulated over fast thermodynamic
!! steps but will be applied to the slow (mass-changing) thermodynamics.  Some
!! of these are diagnostics, while others are averages of fluxes taken during
!! the fast ice thermodynamics and used during the slow ice thermodynamics or dynamics.
type fast_ice_avg_type
!FAST ONLY
  integer :: avg_count   !< The number of times that surface fluxes to the ice have been incremented.
  logical :: atmos_winds !< If true, the wind stresses come directly from the atmosphere model
                         !! and have the wrong sign.
  ! These are the arrays that are averaged over the fast thermodynamics.  They
  ! are either used to communicate to the slow thermodynamics or diagnostics or
  ! both.
  real, allocatable, dimension(:,:,:) :: &
    ! The 3rd dimension in each of the following is ice thickness category.
    flux_u_top  , & !< The downward flux of zonal momentum on an A-grid [R Z L T-2 ~> Pa].
    flux_v_top  , & !< The downward flux of meridional momentum on an A-grid [R Z L T-2 ~> Pa].
    flux_sh_top , & !< The upward sensible heat flux at the ice top [Q R Z T-1 ~> W m-2].
    evap_top    , & !< The upward evaporative moisture flux at top of the ice [R Z T-1 ~> kg m-2 s-1].
    flux_lw_top , & !< The net downward flux of longwave radiation at the top of the ice [Q R Z T-1 ~> W m-2].
    flux_lh_top , & !< The upward flux of latent heat at the top of the ice [Q R Z T-1 ~> W m-2].
    lprec_top   , & !< The downward flux of liquid precipitation at the top of the ice [R Z T-1 ~> kg m-2 s-1].
    fprec_top   , & !< The downward flux of frozen precipitation at the top of the ice [R Z T-1 ~> kg m-2 s-1].
    tmelt       , & !< Ice-top melt energy into the ice/snow [Q R Z ~> J m-2].
    bmelt       , & !< Ice-bottom melting energy into the ice [Q R Z ~> J m-2].
    Tskin_cat       !< The ice skin temperature by category [degC].
  real, allocatable, dimension(:,:,:) ::  sw_abs_ocn !< The fraction of the absorbed
                    !! shortwave radiation that is absorbed in the ocean, <=1, [nondim].
                    !! Equivalent sw_abs_ocn fields are in both the fast_ice_avg_type and the
                    !! ice_rad_type because it is used as a part of the slow thermodynamic updates.
  ! The last dimension in each of the following is angular and frequency radiation band.
  real, allocatable, dimension(:,:,:,:) :: flux_sw_top
                    !< The downward flux of shortwave radiation at the top of the sea-ice [Q R Z T-1 ~> W m-2].
                    !! The fourth dimension combines angular orientation (direct or diffuse) and
                    !! frequency (visible or near-IR) bands, with the integer parameters
                    !! from this module helping to distinguish them.
  real, allocatable, dimension(:,:,:) :: flux_sw_dn !< The total downward shortwave flux
                    !! by wavelength band, averaged across all thickness categories [Q R Z T-1 ~> W m-2].
  real, allocatable, dimension(:,:) :: &
    WindStr_x  , &  !< The zonal wind stress averaged over the ice categories on an A-grid [R L Z T-2 ~> Pa].
    WindStr_y  , &  !< The meridional wind stress averaged over the ice categories on an A-grid [R L Z T-2 ~> Pa].
    WindStr_ocn_x, & !< The zonal wind stress on open water on an A-grid [R L Z T-2 ~> Pa].
    WindStr_ocn_y, & !< The meridional wind stress on open water on an A-grid [R L Z T-2 ~> Pa].
    p_atm_surf , &  !< The atmospheric pressure at the top of the ice [R L Z T-2 ~> Pa].
    runoff, &       !< Liquid runoff into the ocean [R Z T-1 ~> kg m-2].
    calving         !< Calving of ice or runoff of frozen fresh  water into the ocean [R Z T-1 ~> kg m-2].
  real, allocatable, dimension(:,:) :: runoff_hflx !< The heat flux associated with runoff, based
                    !! on the temperature difference relative to a reference temperature [Q R Z T-1 ~> W m-2]
  real, allocatable, dimension(:,:) :: calving_hflx !< The heat flux associated with calving, based
                    !! on the temperature difference relative to a reference temperature [Q R Z T-1 ~> W m-2]
  real, allocatable, dimension(:,:) :: calving_preberg !< Calving of ice or runoff of frozen fresh
                    !! water into the ocean, exclusive of any iceberg contributions [R Z T-1 ~> kg m-2].
  real, allocatable, dimension(:,:) :: calving_hflx_preberg !< The heat flux associated with calving
                    !! exclusive of any iceberg contributions, based on the temperature difference
                    !! relative to a reference temperature [Q R Z T-1 ~> W m-2]
  real, allocatable, dimension(:,:) :: Tskin_avg !< The area-weighted average skin temperature
                    !! across all ice thickness categories [degC], or 0 if there is no ice.
  real, allocatable, dimension(:,:) :: ice_free  !< The fractional open water used in calculating
                    !! WindStr_[xy]_A, between 0 & 1 [nondim].
  real, allocatable, dimension(:,:) :: ice_cover !< The fractional ice coverage, summed across all
                    !! thickness categories, used in calculating WindStr_[xy]_A, between 0 & 1 [nondim].q

  integer :: copy_calls = 0 !< The number of times this structure has been
                    !! copied from the fast ice to the slow ice.
  type (coupler_3d_bc_type) :: &
    tr_flux         !< A structure of additional tracer fluxes at the top of the sea-ice

  ! These are the arrays that are averaged over the fast thermodynamics and
  ! then interpolated into unoccupied categories for the purpose of redoing
  ! the application of the fast thermodynamics
  real, allocatable, dimension(:,:,:) ::  flux_sh0 !< The upward sensible heat flux at the ice top
                !! extrapolated to a skin temperature of 0 degC [Q R Z T-1 ~> W m-2].
  real, allocatable, dimension(:,:,:) ::  evap0 !< The upward evaporative moisture flux
                !! at the top of the ice extrapolated to a skin temperature of 0 degC [R Z T-1 ~> kg m-2 s-1].
  real, allocatable, dimension(:,:,:) ::  flux_lw0 !< The net downward flux of longwave radiation
                !! at the top of the  ice extrapolated to a skin temperature of 0 degC [Q R Z T-1 ~> W m-2].
  real, allocatable, dimension(:,:,:) :: &
    dshdt, &    !< The partial derivative of flux_sh0 with ice skin temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].
    devapdt, &  !< The partial derivative of evap0 with ice skin temperature [R Z T-1 degC-1 ~> kg m-2 s-1 degC-1].
    dlwdt       !< The partial derivative of flux_lw0 with ice skin temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].

!SLOW ONLY
  real, allocatable, dimension(:,:) :: frazil_left !< The frazil heat flux that has not yet been
                    !! consumed in making ice [Q R Z ~> J m-2]. This array is decremented by the ice
                    !! model as the heat flux is used up.
!SLOW ONLY
  !!@{ Diagnostic IDs
  integer :: id_sh=-1, id_lh=-1, id_sw=-1, id_slp=-1
  integer :: id_lw=-1, id_snofl=-1, id_rain=-1,  id_evap=-1
  integer :: id_sw_vis_dir=-1, id_sw_vis_dif=-1, id_sw_nir_dir=-1, id_sw_nir_dif=-1
  integer :: id_sw_vis=-1, id_sw_dir=-1, id_sw_dif=-1, id_sw_dn=-1, id_albedo=-1
  integer :: id_runoff=-1, id_calving=-1, id_runoff_hflx=-1, id_calving_hflx=-1
  integer :: id_tmelt=-1, id_bmelt=-1, id_bheat=-1
  integer :: id_tsfc=-1, id_sitemptop=-1, id_sitemptop_CMOR=-1

  integer :: id_evap_cat=-1, id_lw_cat=-1, id_sh_cat=-1, id_tsfc_cat=-1
  integer :: id_evap0=-1, id_lw0=-1, id_sh0=-1
  integer :: id_devdt=-1, id_dlwdt=-1, id_dshdt=-1
  !!@}
end type fast_ice_avg_type

!> total_sfc_flux_type contains variables that describe the fluxes between the
!! atmosphere and the ice or ocean that have been accumulated over fast thermodynamic
!! steps and integrated across the part-size categories.
type total_sfc_flux_type

  ! These are the arrays that are averaged over the categories and in time over
  ! the fast thermodynamics.
  real, allocatable, dimension(:,:) :: &
    flux_u  , & !< The downward flux of zonal momentum on an A-grid [R L Z T-2 ~> Pa].
    flux_v  , & !< The downward flux of meridional momentum on an A-grid [R L Z T-2 ~> Pa].
    flux_sh , & !< The upward sensible heat flux at the ice top [Q R Z T-1 ~> W m-2].
    evap    , & !< The upward evaporative moisture flux at top of the ice [R Z T-1 ~> kg m-2 s-1].
    flux_lw , & !< The downward flux of longwave radiation at  the top of the ice [Q R Z T-1 ~> W m-2].
    flux_lh , & !< The upward flux of latent heat at the top of the ice [Q R Z T-1 ~> W m-2].
    lprec   , & !< The downward flux of liquid precipitation  at the top of the ice [R Z T-1 ~> kg m-2 s-1].
    fprec       !< The downward flux of frozen precipitation at the top of the ice [R Z T-1 ~> kg m-2 s-1].
  real, allocatable, dimension(:,:,:) :: flux_sw
                !< The downward flux of shortwave radiation at the top of the sea-ice [Q R Z T-1 ~> W m-2].
                !! The third dimension combines angular orientation (direct or diffuse) and
                !! frequency (visible or near-IR) bands, with the integer parameters
                !! from this module helping to distinguish them.
  integer :: copy_calls = 0  !< The number of times this structure has been
                !! copied from the fast ice to the slow ice.
  type (coupler_2d_bc_type) :: &
    tr_flux     !< A structure of additional tracer fluxes at the top of the sea-ice
end type total_sfc_flux_type


!> ice_rad_type contains variables that describe the absorption and reflection
!! of shortwave radiation in and around the sea ice.
type ice_rad_type

  ! The ice skin temperature that can next be used for radiation
  real, allocatable, dimension(:,:,:) :: &
    t_skin, &   !< The surface skin temperature as calculated by the most
                !! recent fast atmospheric timestep, or a value filled in
                !! from other ice categories or the local freezing point of
                !! seawater when there is no ice at all [degC].
    Tskin_Rad   !< The surface skin temperature that was most recently used in
                !! ice optics calculations [degC].
  ! Shortwave absorption parameters that are set in ice_optics.
  real, allocatable, dimension(:,:,:) :: &
    sw_abs_sfc , &  !< The fraction of the absorbed shortwave radiation that is
                    !! absorbed in a surface skin layer, <=1, [nondim].
    sw_abs_snow, &  !< The fraction of the absorbed shortwave radiation that is
                    !! absorbed in the snow, <=1, [nondim].
    sw_abs_ocn , &  !< The fraction of the absorbed shortwave radiation that is
                    !! absorbed in the ocean, <=1, [nondim].
                    !  Only sw_abs_ocn is used in the slow step.
    sw_abs_int      !< The fraction of the absorbed shortwave radiation that is
                    !! absorbed by all ice layers in aggregate, <=1, [nondim].
                    !  sw_abs_int is only used for diagnostics.
  real, allocatable, dimension(:,:,:,:) :: &
    sw_abs_ice      !< The fraction of the absorbed shortwave that is
                    !! absorbed in each of the ice layers, <=1, [nondim].

  real, allocatable, dimension(:,:)   :: &
    coszen_lastrad, & !< Cosine of the solar zenith angle averaged
                    !! over the last radiation timestep [nondim].
    coszen_nextrad  !< Cosine of the solar zenith angle averaged
                    !! over the next radiation timestep [nondim].

  logical :: add_diurnal_sw       !< If true, apply a synthetic diurnal cycle to
                                  !! the shortwave radiation.
  logical :: do_sun_angle_for_alb !< If true, find the sun angle for calculating
                                  !! the ocean albedo in the frame of the ice model.

  !!@{ Diagnostic IDs
  integer, allocatable, dimension(:)   :: id_sw_abs_ice
  integer :: id_sw_abs_sfc=-1, id_sw_abs_snow=-1, id_sw_pen=-1, id_sw_abs_ocn=-1
  integer :: id_alb=-1, id_coszen=-1, id_swdn=-1, id_lwdn=-1
  integer :: id_alb_vis_dir=-1, id_alb_vis_dif=-1, id_alb_nir_dir=-1, id_alb_nir_dif=-1
  integer :: id_tskin=-1, id_cn=-1, id_mi=-1
  !!@}

end type ice_rad_type

!> ice_ocean_flux_type contains variables that describe the fluxes between the
!! ice and the ocean, on the ice grid.
type ice_ocean_flux_type
  ! These variables describe the fluxes between ice or atmosphere and the ocean.
  real, allocatable, dimension(:,:)   :: &
    flux_sh_ocn_top, & !< The upward sensible heat flux from the ocean to the ice or atmosphere [Q R Z T-1 ~> W m-2].
    evap_ocn_top, &    !< The upward evaporative moisture flux at the ocean surface [R Z T-1 ~> kg m-2 s-1].
    flux_lw_ocn_top, & !< The downward flux of longwave radiation at the ocean surface [Q R Z T-1 ~> W m-2].
    flux_lh_ocn_top, & !< The upward flux of latent heat at the ocean surface [Q R Z T-1 ~> W m-2].
    lprec_ocn_top, &   !< The downward flux of liquid precipitation at the ocean surface [R Z T-1 ~> kg m-2 s-1].
    fprec_ocn_top, &   !< The downward flux of frozen precipitation at the ocean surface [R Z T-1 ~> kg m-2 s-1].
    flux_u_ocn, &      !< The flux of x-momentum into the ocean at locations given by
                       !! flux_uv_stagger [R Z L T-2 ~> Pa].
                       !! Note that regardless of the staggering, flux_u_ocn is allocated as though on an A-grid.
    flux_v_ocn, &      !< The flux of y-momentum into the ocean at locations given by
                       !! flux_uv_stagger [R Z L T-2 ~> Pa].
                       !! Note that regardless of the staggering, flux_v_ocn is allocated as though on an A-grid.
    stress_mag, &      !< The area-weighted time-mean of the magnitude of the stress on the ocean [R Z L T-2 ~> Pa].
    melt_nudge, &      !< A downward fresh water flux into the ocean that acts to nudge the ocean
                       !! surface salinity to facilitate the retention of sea ice [R Z T-1 ~> kg m-2 s-1].
    flux_salt, &       !< The flux of salt out of the ocean [kgSalt kg-1 R Z T-1 ~> kgSalt m-2 s-1].
    transmutation_salt_flux, & !< The difference between the salt flux extracted from the ice and the
                       !! salt flux added to the ocean when the ice is transmuted directly into seawater
                       !! as a form of open boundary condition [kgSalt kg-1 R Z T-1 ~> kgSalt m-2 s-1].
    mass_ice_sn_p, &   !< The combined mass per unit ocean area of ice, snow and pond water [R Z ~> kg m-2].
    pres_ocn_top       !< The hydrostatic pressure at the ocean surface due to the weight of ice,
                       !! snow and ponds, exclusive of atmospheric pressure [R Z L T-2 ~> Pa].
                       !### What about pressure from bergs?
  real, allocatable, dimension(:,:,:) :: flux_sw_ocn !< The downward flux of shortwave radiation
                       !! at the ocean surface [Q R Z T-1 ~> W m-2].  The third dimension combines
                       !! angular orientation (direct or diffuse) and frequency
                       !! (visible or near-IR) bands, with the integer parameters
                       !! from this module helping to distinguish them.

  ! Iceberg fields - these are passed unchanged from the icebergs module, so are not rescaled.
  real, pointer, dimension(:,:)   :: &
    ustar_berg => NULL(), & !< ustar contribution below icebergs [m s-1]
    area_berg => NULL(),  & !< fraction of grid cell covered by icebergs [m2 m-2]
    mass_berg => NULL()     !< mass of icebergs [kg m-2]

  ! These arrays are used for enthalpy change diagnostics in the slow thermodynamics.
  real, allocatable, dimension(:,:)   :: &
    ! These terms diagnose the enthalpy change associated with the addition or
    ! removal of water mass (liquid or frozen) from the ice model are required
    ! to close the enthalpy budget. Ice enthalpy is generally negative, so terms
    ! that add mass to the ice are generally negative.
    Enth_Mass_in_atm , & !< The enthalpy introduced to the ice by water fluxes from the atmosphere [Q R Z ~> J m-2].
    Enth_Mass_out_atm, & !< Negative of the enthalpy extracted from the ice by water fluxes to
                         !! the atmosphere [Q R Z ~> J m-2].
    Enth_Mass_in_ocn , & !< The enthalpy introduced to the ice by water fluxes from the ocean [Q R Z ~> J m-2].
    Enth_Mass_out_ocn, & !< Negative of the enthalpy extracted from the ice by water fluxes to
                         !! the ocean [Q R Z ~> J m-2].
    transmutation_enth   !< The difference between the enthalpy extracted from the ice and the
                         !! enthalpy added to the ocean when the ice is transmuted directly into
                         !! seawater as a form of open boundary condition [Q R Z ~> J m-2].

  integer :: stress_count !< The number of times that the stresses from the ice to the ocean have been incremented.
  integer :: flux_uv_stagger = -999 !< The staggering relative to the tracer points of the two wind
                        !! stress components. Valid entries include AGRID, BGRID_NE, CGRID_NE,
                        !! BGRID_SW, and CGRID_SW, following the Arakawa grid-staggering notation.
                        !! (These are named integers taken from mpp_parameter_mod.)
                        !! Following SIS, this is BGRID_NE by default when the sea ice is initialized,
                        !! but flux_uv_stagger is set to -999 here so that a global max across ice and
                        !! non-ice processors can be used to determine its value.
  logical :: slp2ocean  !< If true, apply sea level pressure to ocean surface.

  type (coupler_2d_bc_type) :: &
    tr_flux_ocn_top !< A structure of additional tracer fluxes at the top of the ocean

  !>@{ diagnostic IDs for ice-to-ocean fluxes.
  integer :: id_saltf=-1
  ! The following are diagnostic IDs for iceberg-related fields.  These are only
  ! used if the iceberg code is activated.
  integer ::  id_ustar_berg=-1, id_area_berg=-1, id_mass_berg=-1
  !!@}
end type ice_ocean_flux_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_IST_arrays allocates the arrays in an ice_state_type.
subroutine alloc_IST_arrays(HI, IG, IST, omit_velocities, omit_Tsurf, do_ridging)
  type(hor_index_type), intent(in)    :: HI  !< The horizontal index type describing the domain
  type(ice_grid_type),  intent(in)    :: IG  !< The sea-ice specific grid type
  type(ice_state_type), intent(inout) :: IST !< A type describing the state of the sea ice
  logical,    optional, intent(in)    :: omit_velocities !< If true, do not allocate velocity arrays
  logical,    optional, intent(in)    :: omit_Tsurf !< If true, do not allocate the surface temperature array
  logical,    optional, intent(in)    :: do_ridging !< If true, allocate arrays related to ridging

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  integer :: isd, ied, jsd, jed, CatIce, NkIce
  logical :: do_vel, do_Tsurf

  do_vel = .true. ; if (present(omit_velocities)) do_vel = .not.omit_velocities
  do_Tsurf = .true. ; if (present(omit_Tsurf)) do_Tsurf = .not.omit_Tsurf

  CatIce = IG%CatIce ; NkIce = IG%NkIce
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  IST%valid_IST = .true.
  allocate(IST%part_size(isd:ied, jsd:jed, 0:CatIce)) ; IST%part_size(:,:,:) = 0.0
  allocate(IST%mH_pond(  isd:ied, jsd:jed, CatIce)) ; IST%mH_pond(:,:,:) = 0.0
  allocate(IST%mH_snow(  isd:ied, jsd:jed, CatIce)) ; IST%mH_snow(:,:,:) = 0.0
  allocate(IST%enth_snow(isd:ied, jsd:jed, CatIce, 1)) ; IST%enth_snow(:,:,:,:) = 0.0
  allocate(IST%mH_ice(   isd:ied, jsd:jed, CatIce)) ; IST%mH_ice(:,:,:) = 0.0
  allocate(IST%enth_ice( isd:ied, jsd:jed, CatIce, NkIce)) ; IST%enth_ice(:,:,:,:) = 0.0
  allocate(IST%sal_ice(  isd:ied, jsd:jed, CatIce, NkIce)) ; IST%sal_ice(:,:,:,:) = 0.0

  if (present(do_ridging)) then ; if (do_ridging) then
    allocate(IST%snow_to_ocn(isd:ied, jsd:jed)) ; IST%snow_to_ocn(:,:) = 0.0
    allocate(IST%water_to_ocn(isd:ied, jsd:jed)) ; IST%water_to_ocn(:,:) = 0.0
    allocate(IST%enth_snow_to_ocn(isd:ied, jsd:jed)) ; IST%enth_snow_to_ocn(:,:) = 0.0
    allocate(IST%rdg_rate(isd:ied, jsd:jed)) ; IST%rdg_rate(:,:) = 0.0
    allocate(IST%rdg_mice(isd:ied, jsd:jed, CatIce)) ; IST%rdg_mice(:,:,:) = 0.0
  endif ; endif

  if (do_vel) then
    ! These velocities are only required for the slow ice processes, and hence
    ! can use the memory macros.
    allocate(IST%u_ice_C(SZIB_(HI), SZJ_(HI))) ; IST%u_ice_C(:,:) = 0.0
    allocate(IST%v_ice_C(SZI_(HI), SZJB_(HI))) ; IST%v_ice_C(:,:) = 0.0
    if (.not.IST%Cgrid_dyn) then
      allocate(IST%u_ice_B(SZIB_(HI), SZJB_(HI))) ; IST%u_ice_B(:,:) = 0.0
      allocate(IST%v_ice_B(SZIB_(HI), SZJB_(HI))) ; IST%v_ice_B(:,:) = 0.0
    endif
  endif

  if (do_Tsurf) then
    ! IST%tsurf is only used with some older options.
    allocate(IST%t_surf(isd:ied, jsd:jed, CatIce)) ; IST%t_surf(:,:,:) = T_0degC
  endif

end subroutine alloc_IST_arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_state_register_restarts registers any variables in the ice state type
!!     that need to be includedin the restart files.
subroutine ice_state_register_restarts(IST, G, IG, Ice_restart)
  type(ice_state_type),    intent(inout) :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(in)    :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  type(SIS_restart_CS),    pointer       :: Ice_restart !< The control structure for the ice restarts

  ! Now register some of these arrays to be read from the restart files.
  if (associated(Ice_restart)) then
    call register_restart_field(Ice_restart, 'part_size', IST%part_size, dim_3='cat0')
    if (allocated(IST%t_surf)) then
      call register_restart_field(Ice_restart, 't_surf_ice', IST%t_surf, &
                                  mandatory=.false., units="deg K")
    endif
    call register_restart_field(Ice_restart, 'h_pond', IST%mH_pond, &
                                mandatory=.false., units="H_to_kg_m2 kg m-2")
    call register_restart_field(Ice_restart, 'h_snow', IST%mH_snow, &
                                mandatory=.true., units="H_to_kg_m2 kg m-2")
    call register_restart_field(Ice_restart, 'enth_snow', IST%enth_snow, dim_4='z_snow', &
                                mandatory=.false.)
    call register_restart_field(Ice_restart, 'h_ice', IST%mH_ice, &
                                mandatory=.true., units="H_to_kg_m2 kg m-2")
    call register_restart_field(Ice_restart, 'H_to_kg_m2', IG%H_to_kg_m2, &
                                longname="The conversion factor from SIS2 mass-thickness units to kg m-2.", &
                                mandatory=.false.)

    call register_restart_field(Ice_restart, 'enth_ice', IST%enth_ice, &
                                mandatory=.false., units="J kg-1")
    call register_restart_field(Ice_restart, 'sal_ice', IST%sal_ice, &
                                mandatory=.false., units="kg/kg")

    if (allocated(IST%snow_to_ocn)) then
      call register_restart_field(Ice_restart, 'snow_to_ocn', IST%snow_to_ocn, &
                                  mandatory=.false., units="kg m-2")
      call register_restart_field(Ice_restart, 'enth_snow_to_ocn', IST%enth_snow_to_ocn, &
                                  mandatory=.false., units="J kg-1")
    endif

    if (IST%Cgrid_dyn) then
      if (G%symmetric) then
        call register_restart_field(Ice_restart, 'sym_u_ice_C', IST%u_ice_C, &
                                    position=EAST_FACE, mandatory=.false.)
        call register_restart_field(Ice_restart, 'sym_v_ice_C', IST%v_ice_C, &
                                    position=NORTH_FACE, mandatory=.false.)
      else
        call register_restart_field(Ice_restart, 'u_ice_C', IST%u_ice_C, &
                                    position=EAST_FACE, mandatory=.false.)
        call register_restart_field(Ice_restart, 'v_ice_C', IST%v_ice_C, &
                                    position=NORTH_FACE, mandatory=.false.)
      endif
    else
      if (G%symmetric) then
        call register_restart_field(Ice_restart, 'sym_u_ice_B', IST%u_ice_B, &
                                    position=CORNER, mandatory=.false.)
        call register_restart_field(Ice_restart, 'sym_v_ice_B', IST%v_ice_B, &
                                    position=CORNER, mandatory=.false.)
      else
        call register_restart_field(Ice_restart, 'u_ice', IST%u_ice_B, &
                                    position=CORNER, mandatory=.false.)
        call register_restart_field(Ice_restart, 'v_ice', IST%v_ice_B, &
                                    position=CORNER, mandatory=.false.)
      endif
    endif
  endif

end subroutine ice_state_register_restarts

subroutine register_unit_conversion_restarts(US, Ice_restart)
  type(unit_scale_type),   intent(inout) :: US    !< A structure with unit conversion factors
  type(SIS_restart_CS),    pointer       :: Ice_restart !< The control structure for the ice restarts

  ! Register scalar unit conversion factors.
  call register_restart_field(Ice_restart, "m_to_Z", US%m_to_Z_restart, &
                                 longname="The conversion factor from m to SIS2 height units.", &
                                 units="Z meter-1", mandatory=.false.)
  call register_restart_field(Ice_restart, "m_to_L", US%m_to_L_restart, &
                                 longname="The conversion factor from m to SIS2 length units.", &
                                 units="L meter-1", mandatory=.false.)
  call register_restart_field(Ice_restart, "s_to_T", US%s_to_T_restart, &
                                 longname="The conversion factor from s to SIS2 time units.", &
                                 units="T second-1", mandatory=.false.)
  call register_restart_field(Ice_restart, "kg_m3_to_R", US%kg_m3_to_R_restart, &
                                 longname="The conversion factor from kg m-3 to SIS2 density units.", &
                                 units="R m3 kg-1", mandatory=.false.)
  call register_restart_field(Ice_restart, "J_kg_to_Q", US%J_kg_to_Q_restart, &
                                 longname="The conversion factor from J kg-1 to SIS2 enthalpy units.", &
                                 units="Q kg J-1", mandatory=.false.)

end subroutine register_unit_conversion_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_state_read_alt_restarts reads in alternative variables that might have been in the restart
!! file, specifically dealing with changing between symmetric and non-symmetric memory restart files.
subroutine ice_state_read_alt_restarts(IST, G, IG, Ice_restart, restart_dir)
  type(ice_state_type),    intent(inout) :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(in)    :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  type(SIS_restart_CS),    pointer       :: Ice_restart !< The control structure for the ice restarts
  character(len=*),        intent(in)    :: restart_dir !< A directory in which to find the restart file

  ! These are temporary variables that will be used only here for reading and then discarded.
  real, allocatable, target, dimension(:,:) :: u_tmp, v_tmp
  type(MOM_domain_type),   pointer :: domain_tmp => NULL()
  logical :: u_set, v_set, read_u, read_v
  integer :: i, j

  if (.not.associated(Ice_restart)) return

  if (G%symmetric) then

    if (IST%Cgrid_dyn) then
      u_set = query_initialized(Ice_restart, 'sym_u_ice_C')
      v_set = query_initialized(Ice_restart, 'sym_v_ice_C')
    else
      u_set = query_initialized(Ice_restart, 'sym_u_ice_B')
      v_set = query_initialized(Ice_restart, 'sym_v_ice_B')
    endif
    if (u_set .and. v_set) return

    if (u_set .neqv. v_set) call SIS_error(FATAL, "ice_state_read_alt_restarts: "//&
      "Only one of the u and v input variables were successfully read from the restart file.")

    call clone_MOM_domain(G%domain, domain_tmp, symmetric=.false., &
                          domain_name="ice temporary domain")

    if (IST%Cgrid_dyn .and. (.not.u_set)) then
      call safe_alloc(u_tmp, G%isd, G%ied, G%jsd, G%jed)
      call safe_alloc(v_tmp, G%isd, G%ied, G%jsd, G%jed)
      call only_read_from_restarts(Ice_restart, 'u_ice_C', u_tmp, domain_tmp, position=EAST_FACE, &
                   directory=restart_dir, success=read_u)
      call only_read_from_restarts(Ice_restart, 'v_ice_C', v_tmp, domain_tmp, position=NORTH_FACE, &
                   directory=restart_dir, success=read_v)
      if (read_u .and. read_v) then
        ! The non-symmetric variant of this vector has been successfully read.
        call pass_vector(u_tmp, v_tmp, domain_tmp, stagger=CGRID_NE)
        do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
          IST%u_ice_C(I,j) = u_tmp(I,j)
        enddo ; enddo
        do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
          IST%v_ice_C(i,J) = v_tmp(i,J)
        enddo ; enddo
      endif
    endif
    if ((.not.IST%Cgrid_dyn) .and. (.not.u_set)) then
      call safe_alloc(u_tmp, G%isd, G%ied, G%jsd, G%jed)
      call safe_alloc(v_tmp, G%isd, G%ied, G%jsd, G%jed)
      call only_read_from_restarts(Ice_restart, 'u_ice', u_tmp, domain_tmp, position=CORNER, &
                                   directory=restart_dir, success=read_u)
      call only_read_from_restarts(Ice_restart, 'v_ice', v_tmp, domain_tmp, position=CORNER, &
                                   directory=restart_dir, success=read_v)
      if (read_u .and. read_v) then
        ! The non-symmetric variant of this variable has been successfully read.
        call pass_vector(u_tmp, v_tmp, domain_tmp, stagger=BGRID_NE)
        do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
          IST%u_ice_B(I,J) = u_tmp(I,J)
          IST%v_ice_B(I,J) = v_tmp(I,J)
        enddo ; enddo
      endif
    endif

  else  ! .not. symmetric
    if (IST%Cgrid_dyn) then
      u_set = query_initialized(Ice_restart, 'u_ice_C')
      v_set = query_initialized(Ice_restart, 'v_ice_C')
    else
      u_set = query_initialized(Ice_restart, 'u_ice')
      v_set = query_initialized(Ice_restart, 'v_ice')
    endif
    if (u_set .and. v_set) return

    if (u_set .neqv. v_set) call SIS_error(FATAL, "ice_state_read_alt_restarts: "//&
      "Only one of the u and v input variables were successfully read from the restart file.")

    call clone_MOM_domain(G%domain, domain_tmp, symmetric=.true., &
                          domain_name="ice temporary sym")

    if (IST%Cgrid_dyn .and. (.not.u_set)) then
      call safe_alloc(u_tmp, G%isd-1, G%ied, G%jsd, G%jed)
      call safe_alloc(v_tmp, G%isd, G%ied, G%jsd-1, G%jed)
      call only_read_from_restarts(Ice_restart, 'sym_u_ice_C', u_tmp, domain_tmp, &
                                   position=EAST_FACE, directory=restart_dir, success=read_u)
      call only_read_from_restarts(Ice_restart, 'sym_v_ice_C', v_tmp, domain_tmp, &
                                   position=NORTH_FACE, directory=restart_dir, success=read_v)
      if (read_u .and. read_v) then
        ! The symmetric variant of this vector has been successfully read.
        do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
          IST%u_ice_C(I,j) = u_tmp(I,j)
        enddo ; enddo
        do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
          IST%v_ice_C(i,J) = v_tmp(i,J)
        enddo ; enddo
      endif
    endif
    if ((.not.IST%Cgrid_dyn) .and. (.not.u_set)) then
      call safe_alloc(u_tmp, G%isd-1, G%ied, G%jsd-1, G%jed)
      call safe_alloc(v_tmp, G%isd-1, G%ied, G%jsd-1, G%jed)
      call only_read_from_restarts(Ice_restart, 'sym_u_ice_B', u_tmp, domain_tmp, position=CORNER, &
                                   directory=restart_dir, success=read_u)
      call only_read_from_restarts(Ice_restart, 'sym_v_ice_B', v_tmp, domain_tmp, position=CORNER, &
                                   directory=restart_dir, success=read_v)
      if (read_u .and. read_v) then
        ! The symmetric variant of this variable has been successfully read.
        do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
          IST%u_ice_B(I,J) = u_tmp(I,J)
          IST%v_ice_B(I,J) = v_tmp(I,J)
        enddo ; enddo
      endif
    endif
  endif

  deallocate(u_tmp, v_tmp)
  deallocate(domain_tmp%mpp_domain) ; deallocate(domain_tmp)

end subroutine ice_state_read_alt_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> rescale_ice_state_restart_fields handles any changes in dimensional rescaling of ice state
!! variables between what is stored in the restart file and what is done for the current run segment.
subroutine rescale_ice_state_restart_fields(IST, G, US, IG, H_to_kg_m2, Rho_ice, Rho_snow)
  type(ice_state_type),    intent(inout) :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(in)    :: G   !< The horizontal grid type
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  real,                    intent(in)    :: H_to_kg_m2 !< The mass conversion_factor that will be
                                                     !! used for the run [kg m-2 H-1 ~> 1].
  real,                    intent(in)    :: Rho_ice  !< The nominal density of ice [R ~> kg m-2]
  real,                    intent(in)    :: Rho_snow !< The nominal density of snow [R ~> kg m-2]

  ! Local variables
  real :: vel_rescale, Q_rescale, RZ_rescale, H_rescale_ice, H_rescale_snow
  integer :: i, j, k, m

  ! Redo the dimensional rescaling of the ice state type variables as necessary.
  ! The rescaling of the ice and snow thickness are dealt with in ice_model-init so that
  ! older (SIS1) sea ice restart files can be used.
  vel_rescale = 1.0
  if ((US%s_to_T_restart*US%m_to_L_restart /= 0.0) .and. &
      (US%m_to_L*US%s_to_T_restart) /= (US%m_to_L_restart*US%s_to_T)) &
    vel_rescale = (US%m_to_L*US%s_to_T_restart) / (US%m_to_L_restart*US%s_to_T)
  Q_rescale = 1.0
  if ((US%J_kg_to_Q_restart /= 0.0) .and. &
      (US%J_kg_to_Q /= US%J_kg_to_Q_restart)) &
    Q_rescale = US%J_kg_to_Q / US%J_kg_to_Q_restart
  RZ_rescale = 1.0
  if ((US%kg_m3_to_R_restart*US%m_to_Z_restart /= 0.0) .and. &
      (US%kg_m3_to_R*US%m_to_Z) /= (US%kg_m3_to_R_restart*US%m_to_Z_restart)) &
    RZ_rescale = (US%kg_m3_to_R*US%m_to_Z) / (US%kg_m3_to_R_restart*US%m_to_Z_restart)

  ! Determine the thickness rescaling factors that are needed.
  H_rescale_ice = 1.0 ; H_rescale_snow = 1.0
  if (IG%H_to_kg_m2 == -1.0) then
    ! This is an older restart file, and the snow and ice thicknesses are in m.
    H_rescale_ice = US%R_to_kg_m3*Rho_ice / H_to_kg_m2
    H_rescale_snow = US%R_to_kg_m3*Rho_snow / H_to_kg_m2
  elseif (IG%H_to_kg_m2 /= H_to_kg_m2) then
    H_rescale_ice = IG%H_to_kg_m2 / H_to_kg_m2
    H_rescale_snow = H_rescale_ice
  endif

  if (H_rescale_ice /= 1.0) then
    do k=1,IG%CatIce ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      IST%mH_ice(i,j,k) = H_rescale_ice * IST%mH_ice(i,j,k)
    enddo ; enddo ; enddo
  endif
  if (H_rescale_snow /= 1.0) then
    do k=1,IG%CatIce ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      IST%mH_snow(i,j,k) = H_rescale_snow * IST%mH_snow(i,j,k)
    enddo ; enddo ; enddo
  endif

  if (IST%Cgrid_dyn .and. (vel_rescale /= 1.0)) then
    do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
      IST%u_ice_C(I,j) = vel_rescale * IST%u_ice_C(I,j)
    enddo ; enddo
    do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
      IST%v_ice_C(i,J) = vel_rescale * IST%v_ice_C(i,J)
    enddo ; enddo
  endif
  if (.not.IST%Cgrid_dyn .and. (vel_rescale /= 1.0)) then
    do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
      IST%u_ice_B(I,J) = vel_rescale * IST%u_ice_B(I,J)
      IST%v_ice_B(I,J) = vel_rescale * IST%v_ice_B(I,J)
    enddo ; enddo
  endif

  if (Q_rescale /= 1.0) then
    do m=1,IG%NkIce ; do k=1,IG%CatIce ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      IST%enth_ice(i,j,k,m) = Q_rescale * IST%enth_ice(i,j,k,m)
      IST%enth_snow(i,j,k,m) = Q_rescale * IST%enth_snow(i,j,k,m)
    enddo ; enddo ; enddo ; enddo
  endif
  if (allocated(IST%snow_to_ocn) .and. (Q_rescale /= 1.0) .or. (RZ_rescale /= 1.0)) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      IST%snow_to_ocn(i,j) = RZ_rescale * IST%snow_to_ocn(i,j)
      IST%enth_snow_to_ocn(i,j) = Q_rescale * IST%enth_snow_to_ocn(i,j)
    enddo ; enddo
  endif

end subroutine rescale_ice_state_restart_fields

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> rescale_fast_to_slow_restart_fields redoes the dimensional rescaling of the restart fields
!!   that are required to be sent from the fast state to the slow state when
!!   the model is restart.  These are the fields that would be copied via the
!!   subroutines copy_FIA_to_FIA, copy_TSF_to_TSF and copy_Rad_to_Rad, and it
!!   should be called from the fast ice processors when redo_fast_update is true.
subroutine rescale_fast_to_slow_restart_fields(FIA, Rad, TSF, G, US, IG)
  type(fast_ice_avg_type),   pointer     :: FIA     !< The fast ice model's fast_ice_avg_type
  type(ice_rad_type),        pointer     :: Rad     !< The fast ice model's ice_rad_type
  type(total_sfc_flux_type), pointer     :: TSF     !< The fast ice model's total_sfc_flux_type
  type(SIS_hor_grid_type),   intent(in)  :: G       !< The horizontal grid type
  type(unit_scale_type),     intent(in)  :: US      !< A structure with unit conversion factors
  type(ice_grid_type),       intent(in)  :: IG  !< The sea-ice specific grid type

  real :: QRZ_T_rescale, RZ_T_rescale, RZL_T2_rescale ! Rescaling correction factors [all ~> 1.0]
  integer :: i, j, k, b

  QRZ_T_rescale = 1.0 ; RZ_T_rescale = 1.0 ; RZL_T2_rescale = 1.0
  if ((US%J_kg_to_Q_restart*US%kg_m3_to_R_restart*US%s_to_T_restart*US%m_to_Z_restart /= 0.0) .and. &
      ((US%J_kg_to_Q*US%kg_m3_to_R*US%m_to_Z*US%s_to_T_restart) /= &
       (US%J_kg_to_Q_restart*US%kg_m3_to_R_restart*US%m_to_Z_restart*US%s_to_T)) ) &
    QRZ_T_rescale = (US%J_kg_to_Q*US%kg_m3_to_R*US%m_to_Z*US%s_to_T_restart) / &
                    (US%J_kg_to_Q_restart*US%kg_m3_to_R_restart*US%m_to_Z_restart*US%s_to_T)

  if ((US%kg_m3_to_R_restart*US%s_to_T_restart*US%m_to_Z_restart /= 0.0) .and. &
      ((US%kg_m3_to_R*US%m_to_Z*US%s_to_T_restart) /= &
       (US%kg_m3_to_R_restart*US%m_to_Z_restart*US%s_to_T)) ) &
    RZ_T_rescale = (US%kg_m3_to_R*US%m_to_Z*US%s_to_T_restart) / &
                   (US%kg_m3_to_R_restart*US%m_to_Z_restart*US%s_to_T)

  if ((US%kg_m3_to_R_restart*US%s_to_T_restart*US%m_to_L_restart*US%m_to_Z_restart /= 0.0) .and. &
      ((US%kg_m3_to_R*US%m_to_Z*US%m_to_L*US%s_to_T_restart**2) /= &
       (US%kg_m3_to_R_restart*US%m_to_Z_restart*US%m_to_L_restart*US%s_to_T**2)) ) &
    RZL_T2_rescale = (US%kg_m3_to_R*US%m_to_Z*US%m_to_L*US%s_to_T_restart**2) / &
                     (US%kg_m3_to_R_restart*US%m_to_Z_restart*US%m_to_L_restart*US%s_to_T**2)

  if ((QRZ_T_rescale == 1.0) .and. (RZ_T_rescale == 1.0) .and. (RZL_T2_rescale == 1.0)) return

  do k=0,IG%CatIce ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    FIA%flux_sh_top(i,j,k) = QRZ_T_rescale * FIA%flux_sh_top(i,j,k) ! [Q R Z T-1 ~> W m-2]
    FIA%evap_top(i,j,k)    = RZ_T_rescale * FIA%evap_top(i,j,k) ! [R Z T-1 ~> kg m-2 s-1]
    FIA%flux_lw_top(i,j,k) = QRZ_T_rescale * FIA%flux_lw_top(i,j,k) ! [Q R Z T-1 ~> W m-2]
    FIA%flux_lh_top(i,j,k) = QRZ_T_rescale * FIA%flux_lh_top(i,j,k) ! [Q R Z T-1 ~> W m-2]
    FIA%lprec_top(i,j,k)   = RZ_T_rescale * FIA%lprec_top(i,j,k) ! [R Z T-1 ~> kg m-2 s-1]
    FIA%fprec_top(i,j,k)   = RZ_T_rescale * FIA%fprec_top(i,j,k) ! [R Z T-1 ~> kg m-2 s-1]
  enddo ; enddo ; enddo
  do b=1,size(FIA%flux_sw_top,4) ; do k=0,IG%CatIce ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    FIA%flux_sw_top(i,j,k,b) = QRZ_T_rescale * FIA%flux_sw_top(i,j,k,b) ! [Q R Z T-1 ~> W m-2]
  enddo ; enddo ; enddo ; enddo

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    FIA%WindStr_x(i,j) = RZL_T2_rescale * FIA%WindStr_x(i,j) ! [R Z L T-2 ~> Pa]
    FIA%WindStr_y(i,j) = RZL_T2_rescale * FIA%WindStr_y(i,j) ! [R Z L T-2 ~> Pa]
    FIA%WindStr_ocn_x(i,j) = RZL_T2_rescale * FIA%WindStr_ocn_x(i,j) ! [R Z L T-2 ~> Pa]
    FIA%WindStr_ocn_y(i,j) = RZL_T2_rescale * FIA%WindStr_ocn_y(i,j) ! [R Z L T-2 ~> Pa]
    FIA%p_atm_surf(i,j) = RZL_T2_rescale * FIA%p_atm_surf(i,j) ! [R Z L T-2 ~> Pa]
    FIA%runoff(i,j) = RZ_T_rescale * FIA%runoff(i,j) ! [R Z T-1 ~> kg m-2 s-1]
    FIA%calving(i,j) = RZ_T_rescale * FIA%calving(i,j) ! [R Z T-1 ~> kg m-2 s-1]
    FIA%runoff_hflx(i,j) = QRZ_T_rescale * FIA%runoff_hflx(i,j) ! [Q R Z T-1 ~> W m-2]
    FIA%calving_hflx(i,j) = QRZ_T_rescale * FIA%calving_hflx(i,j) ! [Q R Z T-1 ~> W m-2]
  enddo ; enddo
  do b=1,size(FIA%flux_sw_dn,3) ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    FIA%flux_sw_dn(i,j,b) = QRZ_T_rescale * FIA%flux_sw_dn(i,j,b) ! [Q R Z T-1 ~> W m-2]
  enddo ; enddo ; enddo


  if (allocated(FIA%flux_sh0)) then ; do k=0,IG%CatIce ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    FIA%flux_sh0(i,j,k) = QRZ_T_rescale * FIA%flux_sh0(i,j,k) ! [Q R Z T-1 ~> W m-2]
    FIA%flux_lw0(i,j,k) = QRZ_T_rescale * FIA%flux_lw0(i,j,k) ! [Q R Z T-1 ~> W m-2]
    FIA%evap0(i,j,k) = RZ_T_rescale * FIA%evap0(i,j,k) ! [R Z T-1 ~> kg m-2 s-1]
    FIA%dshdt(i,j,k) = QRZ_T_rescale * FIA%dshdt(i,j,k) ! [Q R Z T-1 degC-1 ~> W m-2 degC-1]
    FIA%dlwdt(i,j,k) = QRZ_T_rescale * FIA%dlwdt(i,j,k) ! [Q R Z T-1 degC-1 ~> W m-2 degC-1]
    FIA%devapdt(i,j,k) = RZ_T_rescale * FIA%devapdt(i,j,k) ! [Q R Z T-1 degC-1 ~> kg m-2 s-1 degC-1]
!    ! Do not rescale FIA%Tskin_cat(i,j,k) =  FIA%Tskin_cat(i,j,k)  ! [degC]
  enddo ; enddo ; enddo ; endif

 ! Do not rescale Rad%tskin_rad(i,j) = Rad%tskin_rad(i,j) ! [degC]
 ! Do not rescale Rad%coszen_lastrad(i,j) = Rad%coszen_lastrad(i,j) ! [nondim]

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    TSF%flux_sh(i,j) = QRZ_T_rescale * TSF%flux_sh(i,j) ! [Q R Z T-1 ~> W m-2]
    TSF%flux_lw(i,j) = QRZ_T_rescale * TSF%flux_lw(i,j) ! [Q R Z T-1 ~> W m-2]
    TSF%flux_lh(i,j) = QRZ_T_rescale * TSF%flux_lh(i,j) ! [Q R Z T-1 ~> W m-2]
    TSF%evap(i,j) = RZ_T_rescale * TSF%evap(i,j) ! [R Z T-1 ~> kg m-2 s-1]
    TSF%lprec(i,j) = RZ_T_rescale * TSF%lprec(i,j) ! [R Z T-1 ~> kg m-2 s-1]
    TSF%fprec(i,j) = RZ_T_rescale * TSF%fprec(i,j) ! [R Z T-1 ~> kg m-2 s-1]
  enddo ; enddo
  do b=1,size(TSF%flux_sw,3) ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    TSF%flux_sw(i,j,b) = QRZ_T_rescale * TSF%flux_sw(i,j,b) ! [Q R Z T-1 ~> W m-2]
  enddo ; enddo ; enddo

end subroutine rescale_fast_to_slow_restart_fields

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_fast_ice_avg allocates and zeros out the arrays in a fast_ice_avg_type.
subroutine alloc_fast_ice_avg(FIA, HI, IG, interp_fluxes, gas_fluxes)
  type(fast_ice_avg_type), pointer    :: FIA !< A type containing averages of fields
                                             !! (mostly fluxes) over the fast updates
  type(hor_index_type),    intent(in) :: HI  !< The horizontal index type describing the domain
  type(ice_grid_type),     intent(in) :: IG  !< The sea-ice specific grid type
  logical,                 intent(in) :: interp_fluxes !< If true, allocate fields to permit the
                                             !! interpolation of a linearized version of the
                                             !! fast fluxes into arealess categories.
  type(coupler_1d_bc_type), &
                 optional, intent(in) :: gas_fluxes !< If present, this type describes the
                                             !! additional gas or other tracer fluxes between the
                                             !! ocean, ice, and atmosphere.

  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, CatIce

  if (.not.associated(FIA)) allocate(FIA)
  CatIce = IG%CatIce
  isc = HI%isc ; iec = HI%iec ; jsc = HI%jsc ; jec = HI%jec
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  FIA%avg_count = 0
  allocate(FIA%flux_u_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_u_top(:,:,:) = 0.0
  allocate(FIA%flux_v_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_v_top(:,:,:) = 0.0
  allocate(FIA%flux_sh_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_sh_top(:,:,:) = 0.0
  allocate(FIA%evap_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%evap_top(:,:,:) = 0.0
  allocate(FIA%flux_sw_top(isd:ied, jsd:jed, 0:CatIce, NBANDS)) ; FIA%flux_sw_top(:,:,:,:) = 0.0
  allocate(FIA%flux_lw_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_lw_top(:,:,:) = 0.0
  allocate(FIA%flux_lh_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_lh_top(:,:,:) = 0.0
  allocate(FIA%lprec_top(isd:ied, jsd:jed, 0:CatIce)) ;  FIA%lprec_top(:,:,:) = 0.0
  allocate(FIA%fprec_top(isd:ied, jsd:jed, 0:CatIce)) ;  FIA%fprec_top(:,:,:) = 0.0
  allocate(FIA%runoff(isd:ied, jsd:jed)) ; FIA%runoff(:,:) = 0.0
  allocate(FIA%calving(isd:ied, jsd:jed)) ; FIA%calving(:,:) = 0.0
  allocate(FIA%calving_preberg(isd:ied, jsd:jed)) ; FIA%calving_preberg(:,:) = 0.0 ! diag
  allocate(FIA%runoff_hflx(isd:ied, jsd:jed)) ; FIA%runoff_hflx(:,:) = 0.0
  allocate(FIA%calving_hflx(isd:ied, jsd:jed)) ; FIA%calving_hflx(:,:) = 0.0
  allocate(FIA%calving_hflx_preberg(isd:ied, jsd:jed)) ; FIA%calving_hflx_preberg(:,:) = 0.0 ! diag

  allocate(FIA%frazil_left(isd:ied, jsd:jed)) ; FIA%frazil_left(:,:) = 0.0
  allocate(FIA%tmelt(isd:ied, jsd:jed, CatIce)) ; FIA%tmelt(:,:,:) = 0.0
  allocate(FIA%bmelt(isd:ied, jsd:jed, CatIce)) ; FIA%bmelt(:,:,:) = 0.0
  allocate(FIA%WindStr_x(isd:ied, jsd:jed)) ; FIA%WindStr_x(:,:) = 0.0
  allocate(FIA%WindStr_y(isd:ied, jsd:jed)) ; FIA%WindStr_y(:,:) = 0.0
  allocate(FIA%WindStr_ocn_x(isd:ied, jsd:jed)) ; FIA%WindStr_ocn_x(:,:) = 0.0
  allocate(FIA%WindStr_ocn_y(isd:ied, jsd:jed)) ; FIA%WindStr_ocn_y(:,:) = 0.0
  allocate(FIA%p_atm_surf(isd:ied, jsd:jed)) ; FIA%p_atm_surf(:,:) = 0.0
  allocate(FIA%Tskin_avg(isd:ied, jsd:jed)) ; FIA%Tskin_avg(:,:) = 0.0 ! diag
  allocate(FIA%ice_free(isd:ied, jsd:jed))  ; FIA%ice_free(:,:) = 0.0
  allocate(FIA%ice_cover(isd:ied, jsd:jed)) ; FIA%ice_cover(:,:) = 0.0

  if (interp_fluxes) then
    allocate(FIA%flux_sh0(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_sh0(:,:,:) = 0.0
    allocate(FIA%evap0(isd:ied, jsd:jed, 0:CatIce)) ; FIA%evap0(:,:,:) = 0.0
    allocate(FIA%flux_lw0(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_lw0(:,:,:) = 0.0
    allocate(FIA%dshdt(isd:ied, jsd:jed, 0:CatIce)) ; FIA%dshdt(:,:,:) = 0.0
    allocate(FIA%devapdt(isd:ied, jsd:jed, 0:CatIce)) ; FIA%devapdt(:,:,:) = 0.0
    allocate(FIA%dlwdt(isd:ied, jsd:jed, 0:CatIce)) ; FIA%dlwdt(:,:,:) = 0.0
    allocate(FIA%Tskin_cat(isd:ied, jsd:jed, 0:CatIce)) ; FIA%Tskin_cat(:,:,:) = 0.0
  endif

  allocate(FIA%flux_sw_dn(isd:ied, jsd:jed, NBANDS)) ; FIA%flux_sw_dn(:,:,:) = 0.0
  allocate(FIA%sw_abs_ocn(isd:ied, jsd:jed, CatIce)) ; FIA%sw_abs_ocn(:,:,:) = 0.0

  if (present(gas_fluxes)) &
    call coupler_type_spawn(gas_fluxes, FIA%tr_flux, (/isd, isc, iec, ied/), &
                            (/jsd, jsc, jec, jed/), (/0, CatIce/))

end subroutine alloc_fast_ice_avg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_total_sfc_flux allocates and zeros out the arrays in a total_sfc_flux_type.
subroutine alloc_total_sfc_flux(TSF, HI, gas_fluxes)
  type(total_sfc_flux_type), pointer    :: TSF  !< The total surface flux type being allocated
  type(hor_index_type),      intent(in) :: HI   !< The hor_index_type with information about the
                                                !! array extents to be allocated.
  type(coupler_1d_bc_type), &
                optional, intent(in)    :: gas_fluxes !< If present, this type describes the
                                              !! additional gas or other tracer fluxes between the
                                              !! ocean, ice, and atmosphere.

  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed

  if (.not.associated(TSF)) allocate(TSF)
  isc = HI%isc ; iec = HI%iec ; jsc = HI%jsc ; jec = HI%jec
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  allocate(TSF%flux_u(isd:ied, jsd:jed)) ; TSF%flux_u(:,:) = 0.0
  allocate(TSF%flux_v(isd:ied, jsd:jed)) ; TSF%flux_v(:,:) = 0.0
  allocate(TSF%flux_sh(isd:ied, jsd:jed)) ; TSF%flux_sh(:,:) = 0.0
  allocate(TSF%flux_sw(isd:ied, jsd:jed, NBANDS)) ; TSF%flux_sw(:,:,:) = 0.0
  allocate(TSF%flux_lw(isd:ied, jsd:jed)) ; TSF%flux_lw(:,:) = 0.0
  allocate(TSF%flux_lh(isd:ied, jsd:jed)) ; TSF%flux_lh(:,:) = 0.0
  allocate(TSF%evap(isd:ied, jsd:jed)) ; TSF%evap(:,:) = 0.0
  allocate(TSF%lprec(isd:ied, jsd:jed)) ;  TSF%lprec(:,:) = 0.0
  allocate(TSF%fprec(isd:ied, jsd:jed)) ;  TSF%fprec(:,:) = 0.0
  if (present(gas_fluxes)) &
    call coupler_type_spawn(gas_fluxes, TSF%tr_flux, (/isd, isc, iec, ied/), &
                            (/jsd, jsc, jec, jed/))

end subroutine alloc_total_sfc_flux


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_rad_register_restarts allocates the arrays in the ice_rad_type
!!     and registers any variables in the ice rad type that need to be included
!!     in the restart files.
subroutine ice_rad_register_restarts(HI, IG, param_file, Rad, Ice_restart)
  type(hor_index_type),    intent(in)    :: HI  !< The horizontal index type describing the domain
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ice_rad_type),      pointer       :: Rad !< A structure with fields related to the absorption,
                                                !! reflection and transmission of shortwave radiation.
  type(SIS_restart_CS),    pointer       :: Ice_restart !< The control structure for the ice restarts

  integer :: isd, ied, jsd, jed, CatIce, NkIce

  if (.not.associated(Rad)) allocate(Rad)
  CatIce = IG%CatIce ; NkIce = IG%NkIce
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  call safe_alloc(Rad%t_skin, isd, ied, jsd, jed, CatIce)
  call safe_alloc(Rad%Tskin_rad, isd, ied, jsd, jed, CatIce)

  call safe_alloc(Rad%sw_abs_sfc, isd, ied, jsd, jed, CatIce)
  call safe_alloc(Rad%sw_abs_snow, isd, ied, jsd, jed, CatIce)
  if (.not. allocated(Rad%sw_abs_ice)) then
    allocate(Rad%sw_abs_ice(isd:ied, jsd:jed, CatIce, NkIce)) ; Rad%sw_abs_ice(:,:,:,:) = 0.0
  endif
  call safe_alloc(Rad%sw_abs_ocn, isd, ied, jsd, jed, CatIce)
  call safe_alloc(Rad%sw_abs_int, isd, ied, jsd, jed, CatIce)

  call safe_alloc(Rad%coszen_nextrad, isd, ied, jsd, jed)
  call safe_alloc(Rad%coszen_lastrad, isd, ied, jsd, jed)

  call register_restart_field(Ice_restart, 'coszen', Rad%coszen_nextrad, mandatory=.false.)
  call register_restart_field(Ice_restart, 'T_skin', Rad%t_skin, mandatory=.false.)

end subroutine ice_rad_register_restarts

subroutine alloc_ice_rad(Rad, HI, IG)
  type(ice_rad_type),      pointer       :: Rad !< A structure with fields related to the absorption,
                                                !! reflection and transmission of shortwave radiation.
  type(hor_index_type),    intent(in)    :: HI  !< The horizontal index type describing the domain
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type

  integer :: isd, ied, jsd, jed, CatIce, NkIce

  if (.not.associated(Rad)) allocate(Rad)
  CatIce = IG%CatIce ; NkIce = IG%NkIce
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  allocate(Rad%t_skin(isd:ied, jsd:jed, CatIce)) ; Rad%t_skin(:,:,:) = 0.0
  allocate(Rad%Tskin_rad(isd:ied, jsd:jed, CatIce)) ; Rad%Tskin_rad(:,:,:) = 0.0

  allocate(Rad%sw_abs_sfc(isd:ied, jsd:jed, CatIce)) ; Rad%sw_abs_sfc(:,:,:) = 0.0
  allocate(Rad%sw_abs_snow(isd:ied, jsd:jed, CatIce)) ; Rad%sw_abs_snow(:,:,:) = 0.0
  allocate(Rad%sw_abs_ice(isd:ied, jsd:jed, CatIce, NkIce)) ; Rad%sw_abs_ice(:,:,:,:) = 0.0
  allocate(Rad%sw_abs_ocn(isd:ied, jsd:jed, CatIce)) ; Rad%sw_abs_ocn(:,:,:) = 0.0
  allocate(Rad%sw_abs_int(isd:ied, jsd:jed, CatIce)) ; Rad%sw_abs_int(:,:,:) = 0.0

  allocate(Rad%coszen_nextrad(isd:ied, jsd:jed)) ; Rad%coszen_nextrad(:,:) = 0.0
  allocate(Rad%coszen_lastrad(isd:ied, jsd:jed)) ; Rad%coszen_lastrad(:,:) = 0.0

end subroutine alloc_ice_rad

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_ice_ocean_flux allocates and zeros out the arrays in an ice_ocean_flux_type.
subroutine alloc_ice_ocean_flux(IOF, HI, do_stress_mag, do_iceberg_fields, do_transmute)
  type(ice_ocean_flux_type), pointer    :: IOF !< A structure containing fluxes from the ice to
                                               !! the ocean that are calculated by the ice model.
  type(hor_index_type),      intent(in) :: HI  !< The horizontal index type describing the domain
  logical,         optional, intent(in) :: do_stress_mag !< If true, allocate memory to use for
                                               !! the magnitude of the ice-ocean stress.
  logical,         optional, intent(in) :: do_iceberg_fields !< If true, allocate fields related
                                               !! to exchanges with icebergs
  logical,         optional, intent(in) :: do_transmute !< If true, allocate fields related to
                                               !! transmuting ice directly into seawater as a form
                                               !! of open boundary condition
  integer :: CatIce
  logical :: alloc_bergs, alloc_stress_mag

  alloc_bergs = .false. ; if (present(do_iceberg_fields)) alloc_bergs = do_iceberg_fields
  alloc_stress_mag = .false. ; if (present(do_stress_mag)) alloc_stress_mag = do_stress_mag

  if (.not.associated(IOF)) allocate(IOF)

  allocate(IOF%flux_salt(SZI_(HI), SZJ_(HI))) ; IOF%flux_salt(:,:) = 0.0
  if (do_transmute) then
    allocate(IOF%transmutation_salt_flux(SZI_(HI), SZJ_(HI))) ; IOF%transmutation_salt_flux(:,:) = 0.0
  endif

  allocate(IOF%flux_sh_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sh_ocn_top(:,:) = 0.0
  allocate(IOF%evap_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%evap_ocn_top(:,:) = 0.0
  allocate(IOF%flux_lw_ocn_top(SZI_(HI), SZJ_(HI))) ; IOF%flux_lw_ocn_top(:,:) = 0.0
  allocate(IOF%flux_lh_ocn_top(SZI_(HI), SZJ_(HI))) ; IOF%flux_lh_ocn_top(:,:) = 0.0
  allocate(IOF%flux_sw_ocn(SZI_(HI), SZJ_(HI), NBANDS)) ;  IOF%flux_sw_ocn(:,:,:) = 0.0
  allocate(IOF%lprec_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%lprec_ocn_top(:,:) = 0.0
  allocate(IOF%fprec_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%fprec_ocn_top(:,:) = 0.0
  allocate(IOF%flux_u_ocn(SZI_(HI), SZJ_(HI)))    ;  IOF%flux_u_ocn(:,:) = 0.0
  allocate(IOF%flux_v_ocn(SZI_(HI), SZJ_(HI)))    ;  IOF%flux_v_ocn(:,:) = 0.0
  if (alloc_stress_mag) then
    allocate(IOF%stress_mag(SZI_(HI), SZJ_(HI)))  ;  IOF%stress_mag(:,:) = 0.0
  endif
  allocate(IOF%pres_ocn_top(SZI_(HI), SZJ_(HI)))  ; IOF%pres_ocn_top(:,:) = 0.0
  allocate(IOF%mass_ice_sn_p(SZI_(HI), SZJ_(HI))) ; IOF%mass_ice_sn_p(:,:) = 0.0

  allocate(IOF%Enth_Mass_in_atm(SZI_(HI), SZJ_(HI)))  ; IOF%Enth_Mass_in_atm(:,:) = 0.0
  allocate(IOF%Enth_Mass_out_atm(SZI_(HI), SZJ_(HI))) ; IOF%Enth_Mass_out_atm(:,:) = 0.0
  allocate(IOF%Enth_Mass_in_ocn(SZI_(HI), SZJ_(HI)))  ; IOF%Enth_Mass_in_ocn(:,:) = 0.0
  allocate(IOF%Enth_Mass_out_ocn(SZI_(HI), SZJ_(HI))) ; IOF%Enth_Mass_out_ocn(:,:) = 0.0
  if (do_transmute) then
    allocate(IOF%transmutation_enth(SZI_(HI), SZJ_(HI))) ; IOF%transmutation_enth(:,:) = 0.0
  endif
  ! Allocating iceberg fields (only used if pass_iceberg_area_to_ocean=.True.)
  ! Please note that these are only allocated on the computational domain so that they
  ! can be passed conveniently to the iceberg code.
  if (alloc_bergs) then
    allocate(IOF%mass_berg(HI%isc:HI%iec, HI%jsc:HI%jec)) ; IOF%mass_berg(:,:) = 0.0
    allocate(IOF%ustar_berg(HI%isc:HI%iec, HI%jsc:HI%jec)) ; IOF%ustar_berg(:,:) = 0.0
    allocate(IOF%area_berg(HI%isc:HI%iec, HI%jsc:HI%jec)) ; IOF%area_berg(:,:) = 0.0
  endif

end subroutine alloc_ice_ocean_flux

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_ocean_sfc_state allocates and zeros out the arrays in an ocean_sfc_state_type.
subroutine alloc_ocean_sfc_state(OSS, HI, Cgrid_dyn, gas_fields_ocn)
  type(ocean_sfc_state_type), pointer    :: OSS  !< The ocean_sfc_state_type being allocated
  type(hor_index_type),       intent(in) :: HI   !< The hor_index_type with information about the
                                                 !! array extents to be allocated.
  logical,                    intent(in) :: Cgrid_dyn  !< A variable indicating whether the ice
                                                 !! ice dynamics are calculated on a C-grid (true)
                                                 !! or on a B-grid (false).
  type(coupler_1d_bc_type), &
           optional, intent(in)     :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes.
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  isc = HI%isc ; iec = HI%iec ; jsc = HI%jsc ; jec = HI%jec
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  if (.not.associated(OSS)) allocate(OSS)

  ! The ocean_sfc_state_type only occurs on slow ice PEs, so it can use the memory macros.
  allocate(OSS%s_surf(SZI_(HI), SZJ_(HI))) ; OSS%s_surf(:,:) = 0.0
  allocate(OSS%SST_C(SZI_(HI), SZJ_(HI)))  ; OSS%SST_C(:,:) = 0.0
  allocate(OSS%T_fr_ocn(SZI_(HI), SZJ_(HI))) ; OSS%T_fr_ocn(:,:) = 0.0
  allocate(OSS%sea_lev(SZI_(HI), SZJ_(HI))) ; OSS%sea_lev(:,:) = 0.0
  allocate(OSS%bheat(SZI_(HI), SZJ_(HI)))  ; OSS%bheat(:,:) = 0.0
  allocate(OSS%frazil(SZI_(HI), SZJ_(HI))) ; OSS%frazil(:,:) = 0.0

  if (Cgrid_dyn) then
    allocate(OSS%u_ocn_C(SZIB_(HI), SZJ_(HI))) ; OSS%u_ocn_C(:,:) = 0.0
    allocate(OSS%v_ocn_C(SZI_(HI), SZJB_(HI))) ; OSS%v_ocn_C(:,:) = 0.0
  else
    allocate(OSS%u_ocn_B(SZIB_(HI), SZJB_(HI))) ; OSS%u_ocn_B(:,:) = 0.0
    allocate(OSS%v_ocn_B(SZIB_(HI), SZJB_(HI))) ; OSS%v_ocn_B(:,:) = 0.0
  endif

  OSS%Cgrid_dyn = Cgrid_dyn

  if (present(gas_fields_ocn)) &
    call coupler_type_spawn(gas_fields_ocn, OSS%tr_fields, (/isd, isc, iec, ied/), &
                            (/jsd, jsc, jec, jed/))

end subroutine alloc_ocean_sfc_state


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_simple_ocean_sfc_state allocates and zeros out the arrays in a
!! simple_OSS_type.
subroutine alloc_simple_OSS(OSS, HI, gas_fields_ocn)
  type(simple_OSS_type), pointer    :: OSS    !< The simple_OSS_type being allocated
  type(hor_index_type),  intent(in) :: HI     !< The hor_index_type with information about the
                                              !! array extents to be allocated.
  type(coupler_1d_bc_type), &
           optional, intent(in)     :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes.

  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed

  if (.not.associated(OSS)) allocate(OSS)
  isc = HI%isc ; iec = HI%iec ; jsc = HI%jsc ; jec = HI%jec
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  allocate(OSS%s_surf(isd:ied, jsd:jed)) ; OSS%s_surf(:,:) = 0.0
  allocate(OSS%SST_C(isd:ied, jsd:jed))  ; OSS%SST_C(:,:) = 0.0
  allocate(OSS%T_fr_ocn(isd:ied, jsd:jed)) ; OSS%T_fr_ocn(:,:) = 0.0
  allocate(OSS%bheat(isd:ied, jsd:jed))   ; OSS%bheat(:,:) = 0.0
  allocate(OSS%u_ocn_A(isd:ied, jsd:jed)) ; OSS%u_ocn_A(:,:) = 0.0
  allocate(OSS%v_ocn_A(isd:ied, jsd:jed)) ; OSS%v_ocn_A(:,:) = 0.0
  allocate(OSS%u_ice_A(isd:ied, jsd:jed)) ; OSS%u_ice_A(:,:) = 0.0
  allocate(OSS%v_ice_A(isd:ied, jsd:jed)) ; OSS%v_ice_A(:,:) = 0.0
  if (present(gas_fields_ocn)) &
    call coupler_type_spawn(gas_fields_ocn, OSS%tr_fields, (/isd, isc, iec, ied/), &
                            (/jsd, jsc, jec, jed/))

end subroutine alloc_simple_OSS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> copy_IST_to_IST copies the computational domain of one ice state type into
!! the computational domain of another ice_state_type.  Both must use the same
!! domain decomposition and indexing convention (for now), but they may have
!! different halo sizes.
subroutine copy_IST_to_IST(IST_in, IST_out, HI_in, HI_out, IG)
  type(ice_state_type), intent(in)    :: IST_in !< The ice_state_type that is being copied from
  type(ice_state_type), intent(inout) :: IST_out !< The ice_state_type that is being copied into
  type(hor_index_type), intent(in)    :: HI_in  !< The horizontal index type for the input type
  type(hor_index_type), intent(in)    :: HI_out !< The horizontal index type for the output type
  type(ice_grid_type),  intent(in)    :: IG  !< The sea-ice specific grid type

  integer :: i, j, k, m, isc, iec, jsc, jec, ncat, NkIce
  integer :: i2, j2, i_off, j_off

  isc = HI_in%isc ; iec = HI_in%iec ; jsc = HI_in%jsc ; jec = HI_in%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  if ((HI_in%iec-HI_in%isc /= HI_out%iec-HI_out%isc) .or. &
      (HI_in%jec-HI_in%jsc /= HI_out%jec-HI_out%jsc)) then
    call SIS_error(FATAL, "copy_IST_to_IST called with inconsistent domain "//&
                          "decompositions of the two ice types.")
  endif
  i_off = HI_out%iec-HI_in%iec ;  j_off = HI_out%jec-HI_in%jec

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    IST_out%part_size(i2,j2,k) = IST_in%part_size(i,j,k)
  enddo ; enddo ; enddo

  if (allocated(IST_out%t_surf) .and. allocated(IST_in%t_surf)) then
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      IST_out%t_surf(i2,j2,k) = IST_in%t_surf(i,j,k)
    enddo ; enddo ; enddo
  endif

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    IST_out%mH_pond(i2,j2,k) = IST_in%mH_pond(i,j,k)
    IST_out%mH_snow(i2,j2,k) = IST_in%mH_snow(i,j,k)
    IST_out%mH_ice(i2,j2,k) = IST_in%mH_ice(i,j,k)

    IST_out%enth_snow(i2,j2,k,1) = IST_in%enth_snow(i,j,k,1)
  enddo ; enddo ; enddo

  do m=1,NkIce ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    IST_out%enth_ice(i2,j2,k,m) = IST_in%enth_ice(i,j,k,m)
    IST_out%sal_ice(i2,j2,k,m) = IST_in%sal_ice(i,j,k,m)
  enddo ; enddo ; enddo ; enddo

  ! The velocity components, rdg_mice, TrReg, and ITV are deliberately not being copied.

end subroutine copy_IST_to_IST

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> redistribute_IST_to_IST redistributes the computational domain of one ice state type into
!! the computational domain of another ice_state_type.
subroutine redistribute_IST_to_IST(IST_in, IST_out, domain_in, domain_out)
  type(ice_state_type), pointer    :: IST_in     !< The ice_state_type that is being copied from (intent in).
  type(ice_state_type), pointer    :: IST_out    !< The ice_state_type that is being copied into (intent inout).
  type(domain2d),       intent(in) :: domain_in  !< The source data domain.
  type(domain2d),       intent(in) :: domain_out !< The target data domain.

  real, pointer, dimension(:,:,:) :: null_ptr3D => NULL()
  real, pointer, dimension(:,:,:,:) :: null_ptr4D => NULL()

  ! The velocity components, rdg_mice, TrReg, and ITV are deliberately not being copied.
  if (associated(IST_out) .and. associated(IST_in)) then
    call redistribute_data(domain_in, IST_in%part_size, domain_out, &
                           IST_out%part_size, complete=.true.)

    if (allocated(IST_out%t_surf) .or. allocated(IST_in%t_surf)) then
      call redistribute_data(domain_in, IST_in%t_surf, domain_out, &
                             IST_out%t_surf, complete=.false.)
    endif
    call redistribute_data(domain_in, IST_in%mH_pond, domain_out, &
                           IST_out%mH_pond, complete=.false.)
    call redistribute_data(domain_in, IST_in%mH_snow, domain_out, &
                           IST_out%mH_snow, complete=.false.)
    call redistribute_data(domain_in, IST_in%mH_ice, domain_out, &
                           IST_out%mH_ice, complete=.false.)
    call redistribute_data(domain_in, IST_in%enth_snow, domain_out, &
                           IST_out%enth_snow, complete=.true.)

    call redistribute_data(domain_in, IST_in%enth_ice, domain_out, &
                           IST_out%enth_ice, complete=.false.)
    call redistribute_data(domain_in, IST_in%sal_ice, domain_out, &
                           IST_out%sal_ice, complete=.true.)
  elseif (associated(IST_out)) then
    ! Use the null pointers in place of the unneeded input arrays.
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           IST_out%part_size, complete=.true.)

    if (allocated(IST_out%t_surf)) then
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             IST_out%t_surf, complete=.false.)
    endif
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           IST_out%mH_pond, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           IST_out%mH_snow, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           IST_out%mH_ice, complete=.false.)
    call redistribute_data(domain_in, null_ptr4D, domain_out, &
                           IST_out%enth_snow, complete=.true.)

    call redistribute_data(domain_in, null_ptr4D, domain_out, &
                           IST_out%enth_ice, complete=.false.)
    call redistribute_data(domain_in, null_ptr4D, domain_out, &
                           IST_out%sal_ice, complete=.true.)
  elseif (associated(IST_in)) then
    ! Use the null pointers in place of the unneeded output arrays.
    call redistribute_data(domain_in, IST_in%part_size, domain_out, &
                           null_ptr3D, complete=.true.)

    if (allocated(IST_in%t_surf)) then
      call redistribute_data(domain_in, IST_in%t_surf, domain_out, &
                             null_ptr3D, complete=.false.)
    endif
    call redistribute_data(domain_in, IST_in%mH_pond, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, IST_in%mH_snow, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, IST_in%mH_ice, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, IST_in%enth_snow, domain_out, &
                           null_ptr4D, complete=.true.)

    call redistribute_data(domain_in, IST_in%enth_ice, domain_out, &
                           null_ptr4D, complete=.false.)
    call redistribute_data(domain_in, IST_in%sal_ice, domain_out, &
                           null_ptr4D, complete=.true.)

  else
    call SIS_error(FATAL, "redistribute_IST_to_IST called with "//&
                          "neither IST_in nor IST_out associated.")
  endif

end subroutine redistribute_IST_to_IST

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> translate_OSS_to_sOSS translates the full ocean surface state, as seen by the slow
!! ice processors into a simplified version with the fields that are shared with
!! the atmosphere and the fast ice thermodynamics.
subroutine translate_OSS_to_sOSS(OSS, IST, sOSS, G, US)
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(ice_state_type),       intent(in)    :: IST !< A type describing the state of the sea ice
  type(simple_OSS_type),      intent(inout) :: sOSS !< The simple ocean surface state type that is being copied into
  type(SIS_hor_grid_type),    intent(in)    :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,sOSS,OSS,IST)
  do j=jsc,jec ; do i=isc,iec
    sOSS%s_surf(i,j) = OSS%s_surf(i,j)
    sOSS%SST_C(i,j) = OSS%SST_C(i,j)
    sOSS%T_fr_ocn(i,j) = OSS%T_fr_ocn(i,j)

    if (G%mask2dT(i,j) > 0.5) then
      sOSS%bheat(i,j) = OSS%bheat(i,j)
      ! Interpolate the ocean and ice velocities onto tracer cells.
      if (OSS%Cgrid_dyn) then
        sOSS%u_ocn_A(i,j) = 0.5*(OSS%u_ocn_C(I,j) + OSS%u_ocn_C(I-1,j))
        sOSS%v_ocn_A(i,j) = 0.5*(OSS%v_ocn_C(i,J) + OSS%v_ocn_C(i,J-1))
      else
        sOSS%u_ocn_A(i,j) = 0.25*((OSS%u_ocn_B(I,J) + OSS%u_ocn_B(I-1,J-1)) + &
                                  (OSS%u_ocn_B(I,J-1) + OSS%u_ocn_B(I-1,J)) )
        sOSS%v_ocn_A(i,j) = 0.25*((OSS%v_ocn_B(I,J) + OSS%v_ocn_B(I-1,J-1)) + &
                                  (OSS%v_ocn_B(I,J-1) + OSS%v_ocn_B(I-1,J)) )
      endif
      if (IST%Cgrid_dyn) then
        sOSS%u_ice_A(i,j) = 0.5*(IST%u_ice_C(I,j) + IST%u_ice_C(I-1,j))
        sOSS%v_ice_A(i,j) = 0.5*(IST%v_ice_C(i,J) + IST%v_ice_C(i,J-1))
      else
        sOSS%u_ice_A(i,j) = 0.25*((IST%u_ice_B(I,J) + IST%u_ice_B(I-1,J-1)) + &
                                  (IST%u_ice_B(I,J-1) + IST%u_ice_B(I-1,J)) )
        sOSS%v_ice_A(i,j) = 0.25*((IST%v_ice_B(I,J) + IST%v_ice_B(I-1,J-1)) + &
                                  (IST%v_ice_B(I,J-1) + IST%v_ice_B(I-1,J)) )
      endif
    else ! This is a land point.
      sOSS%bheat(i,j) = 0.0
      sOSS%u_ocn_A(i,j) = 0.0 ; sOSS%v_ocn_A(i,j) = 0.0
      sOSS%u_ice_A(i,j) = 0.0 ; sOSS%v_ice_A(i,j) = 0.0
    endif
  enddo ; enddo

  call coupler_type_spawn(OSS%tr_fields, sOSS%tr_fields, (/isd, isc, iec, ied/), &
                          (/jsd, jsc, jec, jed/), as_needed=.true. )

  call coupler_type_copy_data(OSS%tr_fields, sOSS%tr_fields)

end subroutine translate_OSS_to_sOSS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> copy_sOSS_to_sOSS copies the computational domain of one simple_OSS_type into
!! the computational domain of another simple_OSS_type.  Both must use the same
!! domain decomposition and indexing convention, but they may have different
!! halo sizes.
subroutine copy_sOSS_to_sOSS(OSS_in, OSS_out, HI_in, HI_out)
  type(simple_OSS_type), intent(inout) :: OSS_in  !< The simple ocean surface state type that is being copied from
  type(simple_OSS_type), intent(inout) :: OSS_out !< The simple ocean surface state type that is being copied into
  type(hor_index_type),  intent(in)    :: HI_in  !< The horizontal index type for the input type
  type(hor_index_type),  intent(in)    :: HI_out !< The horizontal index type for the output type

  integer :: i, j, isc, iec, jsc, jec
  integer :: i2, j2, i_off, j_off

  isc = HI_in%isc ; iec = HI_in%iec ; jsc = HI_in%jsc ; jec = HI_in%jec

  if ((HI_in%iec-HI_in%isc /= HI_out%iec-HI_out%isc) .or. &
      (HI_in%jec-HI_in%jsc /= HI_out%jec-HI_out%jsc)) then
    call SIS_error(FATAL, "copy_sOSS_to_sOSS called with inconsistent domain "//&
                          "decompositions of the two ice types.")
  endif
  i_off = HI_out%iec-HI_in%iec ;  j_off = HI_out%jec-HI_in%jec

  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    OSS_out%SST_C(i2,j2) = OSS_in%SST_C(i,j)
    OSS_out%s_surf(i2,j2) = OSS_in%s_surf(i,j)
    OSS_out%T_fr_ocn(i2,j2) = OSS_in%T_fr_ocn(i,j)
    OSS_out%bheat(i2,j2) = OSS_in%bheat(i,j)
    OSS_out%u_ocn_A(i2,j2) = OSS_in%u_ocn_A(i,j)
    OSS_out%v_ocn_A(i2,j2) = OSS_in%v_ocn_A(i,j)
    OSS_out%u_ice_A(i2,j2) = OSS_in%u_ice_A(i,j)
    OSS_out%v_ice_A(i2,j2) = OSS_in%v_ice_A(i,j)
  enddo ; enddo

  call coupler_type_copy_data(OSS_in%tr_fields, OSS_out%tr_fields)

end subroutine copy_sOSS_to_sOSS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> redistribute_sOSS_to_sOSS copies the computational domain of one simple_OSS_type into
!! the computational domain of another simple_OSS_type.  When the source and target
!! simple_OSS_types are on different PE lists, one or the other may be unassociated.
subroutine redistribute_sOSS_to_sOSS(OSS_in, OSS_out, domain_in, domain_out, HI_out)
  type(simple_OSS_type),          pointer    :: OSS_in     !< The simple OSS type that is being copied from.
  type(simple_OSS_type),          pointer    :: OSS_out    !< The simple OSS type that is being copied into.
  type(domain2d),                 intent(in) :: domain_in  !< The source data domain.
  type(domain2d),                 intent(in) :: domain_out !< The target data domain.
  type(hor_index_type), optional, intent(in) :: HI_out     !< The hor_index_type on the target domain; HI_out
                                                           !! may be omitted if this is not a target PE.
  real, pointer, dimension(:,:) :: null_ptr => NULL()
  type(coupler_2d_bc_type)  :: null_bc

  if (.not. (associated(OSS_out) .or. associated(OSS_in))) &
    call SIS_error(FATAL, "redistribute_sOSS_to_sOSS called with "//&
                          "neither OSS_in nor OSS_out associated.")

  if (associated(OSS_out) .and. associated(OSS_in)) then
    ! This could have complete set to .false. if the halo sizes matched.
    call coupler_type_redistribute_data(OSS_in%tr_fields, domain_in, &
                          OSS_out%tr_fields, domain_out, complete=.false.)
    call redistribute_data(domain_in, OSS_in%SST_C, domain_out, &
                           OSS_out%SST_C, complete=.false.)
    call redistribute_data(domain_in, OSS_in%s_surf, domain_out, &
                           OSS_out%s_surf, complete=.false.)
    call redistribute_data(domain_in, OSS_in%T_fr_ocn, domain_out, &
                           OSS_out%T_fr_ocn, complete=.false.)
    call redistribute_data(domain_in, OSS_in%bheat, domain_out, &
                           OSS_out%bheat, complete=.false.)
    call redistribute_data(domain_in, OSS_in%u_ocn_A, domain_out, &
                           OSS_out%u_ocn_A, complete=.false.)
    call redistribute_data(domain_in, OSS_in%v_ocn_A, domain_out, &
                           OSS_out%v_ocn_A, complete=.false.)
    call redistribute_data(domain_in, OSS_in%u_ice_A, domain_out, &
                           OSS_out%u_ice_A, complete=.false.)
    call redistribute_data(domain_in, OSS_in%v_ice_A, domain_out, &
                           OSS_out%v_ice_A, complete=.true.)
  elseif (associated(OSS_out)) then
    ! Use the null pointer in place of the unneeded input arrays.

    call coupler_type_redistribute_data(null_bc, domain_in, &
                          OSS_out%tr_fields, domain_out, complete=.false.)
    call redistribute_data(domain_in, null_ptr, domain_out, &
                           OSS_out%SST_C, complete=.false.)
    call redistribute_data(domain_in, null_ptr, domain_out, &
                           OSS_out%s_surf, complete=.false.)
    call redistribute_data(domain_in, null_ptr, domain_out, &
                           OSS_out%T_fr_ocn, complete=.false.)
    call redistribute_data(domain_in, null_ptr, domain_out, &
                           OSS_out%bheat, complete=.false.)
    call redistribute_data(domain_in, null_ptr, domain_out, &
                           OSS_out%u_ocn_A, complete=.false.)
    call redistribute_data(domain_in, null_ptr, domain_out, &
                           OSS_out%v_ocn_A, complete=.false.)
    call redistribute_data(domain_in, null_ptr, domain_out, &
                           OSS_out%u_ice_A, complete=.false.)
    call redistribute_data(domain_in, null_ptr, domain_out, &
                           OSS_out%v_ice_A, complete=.true.)
  elseif (associated(OSS_in)) then
    ! Use the null pointer in place of the unneeded output arrays.
    call coupler_type_redistribute_data(OSS_in%tr_fields, domain_in, &
                          null_bc, domain_out, complete=.false.)
    call redistribute_data(domain_in, OSS_in%SST_C, domain_out, &
                           null_ptr, complete=.false.)
    call redistribute_data(domain_in, OSS_in%s_surf, domain_out, &
                           null_ptr, complete=.false.)
    call redistribute_data(domain_in, OSS_in%T_fr_ocn, domain_out, &
                           null_ptr, complete=.false.)
    call redistribute_data(domain_in, OSS_in%bheat, domain_out, &
                           null_ptr, complete=.false.)
    call redistribute_data(domain_in, OSS_in%u_ocn_A, domain_out, &
                           null_ptr, complete=.false.)
    call redistribute_data(domain_in, OSS_in%v_ocn_A, domain_out, &
                           null_ptr, complete=.false.)
    call redistribute_data(domain_in, OSS_in%u_ice_A, domain_out, &
                           null_ptr, complete=.false.)
    call redistribute_data(domain_in, OSS_in%v_ice_A, domain_out, &
                           null_ptr, complete=.true.)
  else
    call SIS_error(FATAL, "redistribute_sOSS_to_sOSS called with "//&
                          "neither OSS_in nor OSS_out associated.")
  endif

end subroutine redistribute_sOSS_to_sOSS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> copy_FIA_to_FIA copies the computational domain of one fast_ice_avg_type into
!! the computational domain of another fast_ice_avg_type.  Both must use the same
!! domain decomposition and indexing convention, but they may have different
!! halo sizes.
subroutine copy_FIA_to_FIA(FIA_in, FIA_out, HI_in, HI_out, IG)
  type(fast_ice_avg_type), intent(inout) :: FIA_in   !< The fast_ice_avg_type that is being copied from.
  type(fast_ice_avg_type), intent(inout) :: FIA_out  !< The fast_ice_avg_type that is being copied into.
  type(hor_index_type),    intent(in)    :: HI_in  !< The horizontal index type for the input type
  type(hor_index_type),    intent(in)    :: HI_out !< The horizontal index type for the output type
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type

  integer :: b, i, j, k, m, n, nb, isc, iec, jsc, jec, ncat, NkIce
  integer :: i2, j2, i_off, j_off
  integer :: isd, ied, jsd, jed

  isc = HI_in%isc ; iec = HI_in%iec ; jsc = HI_in%jsc ; jec = HI_in%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce ; nb = size(FIA_in%flux_sw_top,4)

  if ((HI_in%iec-HI_in%isc /= HI_out%iec-HI_out%isc) .or. &
      (HI_in%jec-HI_in%jsc /= HI_out%jec-HI_out%jsc)) then
    call SIS_error(FATAL, "copy_FIA_to_FIA called with inconsistent domain "//&
                          "decompositions of the two ice types.")
  endif
  i_off = HI_out%iec-HI_in%iec ;  j_off = HI_out%jec-HI_in%jec

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    FIA_out%flux_sh_top(i2,j2,k) = FIA_in%flux_sh_top(i,j,k)
    FIA_out%evap_top(i2,j2,k) = FIA_in%evap_top(i,j,k)
    FIA_out%flux_lw_top(i2,j2,k) = FIA_in%flux_lw_top(i,j,k)
    FIA_out%flux_lh_top(i2,j2,k) = FIA_in%flux_lh_top(i,j,k)
    FIA_out%lprec_top(i2,j2,k) = FIA_in%lprec_top(i,j,k)
    FIA_out%fprec_top(i2,j2,k) = FIA_in%fprec_top(i,j,k)
    do b=1,nb ; FIA_out%flux_sw_top(i2,j2,k,b) = FIA_in%flux_sw_top(i,j,k,b) ; enddo
  enddo ; enddo ; enddo

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    FIA_out%tmelt(i2,j2,k) = FIA_in%tmelt(i,j,k)
    FIA_out%bmelt(i2,j2,k) = FIA_in%bmelt(i,j,k)
    FIA_out%sw_abs_ocn(i2,j2,k) = FIA_in%sw_abs_ocn(i,j,k)
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    FIA_out%WindStr_x(i2,j2) = FIA_in%WindStr_x(i,j)
    FIA_out%WindStr_y(i2,j2) = FIA_in%WindStr_y(i,j)
    FIA_out%WindStr_ocn_x(i2,j2) = FIA_in%WindStr_ocn_x(i,j)
    FIA_out%WindStr_ocn_y(i2,j2) = FIA_in%WindStr_ocn_y(i,j)
    FIA_out%p_atm_surf(i2,j2) = FIA_in%p_atm_surf(i,j)
    FIA_out%runoff(i2,j2) = FIA_in%runoff(i,j)
    FIA_out%calving(i2,j2) =  FIA_in%calving(i,j)
    FIA_out%runoff_hflx(i2,j2) = FIA_in%runoff_hflx(i,j)
    FIA_out%calving_hflx(i2,j2) =  FIA_in%calving_hflx(i,j)
    FIA_out%Tskin_avg(i2,j2) = FIA_in%Tskin_avg(i,j)
    FIA_out%ice_free(i2,j2) = FIA_in%ice_free(i,j)
    FIA_out%ice_cover(i2,j2) = FIA_in%ice_cover(i,j)
    do b=1,nb ; FIA_out%flux_sw_dn(i2,j2,b) = FIA_in%flux_sw_dn(i,j,b) ; enddo
  enddo ; enddo
  !   FIA%flux_u_top and flux_v_top are deliberately not being copied, as they
  ! are only needed on the fast_ice_PEs
  !   FIA%frazil_left is deliberately not being copied, as it is only valid on
  ! the slow_ice_PEs.
  !   FIA%calving_preberg and FIA%calving_hflx_preberg are deliberately not
  ! being copied over.
  if (allocated(FIA_out%flux_sh0)) then
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off
      FIA_out%flux_sh0(i2,j2,k) = FIA_in%flux_sh0(i,j,k)
      FIA_out%evap0(i2,j2,k) = FIA_in%evap0(i,j,k)
      FIA_out%flux_lw0(i2,j2,k) = FIA_in%flux_lw0(i,j,k)
      FIA_out%dshdt(i2,j2,k) = FIA_in%dshdt(i,j,k)
      FIA_out%devapdt(i2,j2,k) = FIA_in%devapdt(i,j,k)
      FIA_out%dlwdt(i2,j2,k) = FIA_in%dlwdt(i,j,k)
      FIA_out%Tskin_cat(i2,j2,k) = FIA_in%Tskin_cat(i,j,k)
    enddo ; enddo ; enddo
  endif

  if (FIA_in%copy_calls /= FIA_out%copy_calls) call SIS_error(WARNING, &
    "copy_FIA_to_FIA called an inconsistent number of time for the input and output types.")

  FIA_in%copy_calls = FIA_in%copy_calls + 1 ; FIA_out%copy_calls = FIA_out%copy_calls + 1

  call coupler_type_copy_data(FIA_in%tr_flux, FIA_out%tr_flux)

  ! avg_count, atmos_winds, and the IDs are deliberately not being copied.
end subroutine copy_FIA_to_FIA

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> redistribute_FIA_to_FIA copies the computational domain of one fast_ice_avg_type into
!! the computational domain of another fast_ice_avg_type.
subroutine redistribute_FIA_to_FIA(FIA_in, FIA_out, domain_in, domain_out, G_out, IG)
  type(fast_ice_avg_type), pointer    :: FIA_in     !< The fast_ice_avg_type that is being copied from.
  type(fast_ice_avg_type), pointer    :: FIA_out    !< The fast_ice_avg_type that is being copied into.
  type(domain2d),          intent(in) :: domain_in  !< The source data domain.
  type(domain2d),          intent(in) :: domain_out !< The target data domain.
  type(SIS_hor_grid_type), optional, intent(in) :: G_out !< The horizontal grid on the target domain.
  type(ice_grid_type),     optional, intent(in) :: IG    !< The ice grid on the target domain.

  real, pointer, dimension(:,:) :: null_ptr2D => NULL()
  real, pointer, dimension(:,:,:) :: null_ptr3D => NULL()
  real, pointer, dimension(:,:,:,:) :: null_ptr4D => NULL()
  type(coupler_3d_bc_type)  :: null_bc
  integer :: copy_calls  ! The number of times these FIA_types have been copied.
  integer :: i, j, b, isd, ied, jsd, jed, ncat

  copy_calls = 0
  if (associated(FIA_out)) then
    FIA_out%copy_calls = FIA_out%copy_calls + 1
    copy_calls = FIA_out%copy_calls
  endif
  if (associated(FIA_in)) then
    FIA_in%copy_calls = FIA_in%copy_calls + 1
    copy_calls = max(copy_calls,FIA_in%copy_calls)
  endif

  !   FIA%flux_u_top and flux_v_top are deliberately not being copied, as they
  ! are only needed on the fast_ice_PEs
  !   FIA%frazil_left is deliberately not being copied, as it is only valid on
  ! the slow_ice_PEs.
  !   FIA%calving_preberg and FIA%calving_hflx_preberg are deliberately not
  ! being copied over.
  ! avg_count, atmos_winds, and the IDs are deliberately not being copied.

  if (associated(FIA_out) .and. associated(FIA_in)) then
    call coupler_type_redistribute_data(FIA_in%tr_flux, domain_in, &
                          FIA_out%tr_flux, domain_out, complete=.false.)
    do b=1,size(FIA_in%flux_sw_top,4)
      call redistribute_data(domain_in, FIA_in%flux_sw_top(:,:,:,b), domain_out, &
                             FIA_out%flux_sw_top(:,:,:,b), complete=.false.)
    enddo
    call  redistribute_data(domain_in, FIA_in%flux_sh_top, domain_out, &
                          FIA_out%flux_sh_top, complete=.false.)
    call redistribute_data(domain_in, FIA_in%evap_top, domain_out, &
                           FIA_out%evap_top, complete=.false.)
    call redistribute_data(domain_in, FIA_in%flux_lw_top, domain_out, &
                           FIA_out%flux_lw_top, complete=.false.)
    call redistribute_data(domain_in, FIA_in%flux_lh_top, domain_out, &
                           FIA_out%flux_lh_top, complete=.false.)
    call redistribute_data(domain_in, FIA_in%lprec_top, domain_out, &
                           FIA_out%lprec_top, complete=.false.)
    call redistribute_data(domain_in, FIA_in%fprec_top, domain_out, &
                           FIA_out%fprec_top, complete=.true.)

    call redistribute_data(domain_in, FIA_in%tmelt, domain_out, &
                           FIA_out%tmelt, complete=.false.)
    call redistribute_data(domain_in, FIA_in%bmelt, domain_out, &
                           FIA_out%bmelt, complete=.false.)
    call redistribute_data(domain_in, FIA_in%sw_abs_ocn, domain_out, &
                           FIA_out%sw_abs_ocn, complete=.true.)

    call redistribute_data(domain_in, FIA_in%WindStr_x, domain_out, &
                           FIA_out%WindStr_x, complete=.false.)
    call redistribute_data(domain_in, FIA_in%WindStr_y, domain_out, &
                           FIA_out%WindStr_y, complete=.false.)
    call redistribute_data(domain_in, FIA_in%WindStr_ocn_x, domain_out, &
                           FIA_out%WindStr_ocn_x, complete=.false.)
    call redistribute_data(domain_in, FIA_in%WindStr_ocn_y, domain_out, &
                           FIA_out%WindStr_ocn_y, complete=.false.)
    call redistribute_data(domain_in, FIA_in%p_atm_surf, domain_out, &
                           FIA_out%p_atm_surf, complete=.false.)
    call redistribute_data(domain_in, FIA_in%runoff, domain_out, &
                           FIA_out%runoff, complete=.false.)
    call redistribute_data(domain_in, FIA_in%calving, domain_out, &
                           FIA_out%calving, complete=.false.)
    call redistribute_data(domain_in, FIA_in%runoff_hflx, domain_out, &
                           FIA_out%runoff_hflx, complete=.false.)
    call redistribute_data(domain_in, FIA_in%calving_hflx, domain_out, &
                           FIA_out%calving_hflx, complete=.false.)
    call redistribute_data(domain_in, FIA_in%Tskin_avg, domain_out, &
                           FIA_out%Tskin_avg, complete=.false.)
    do b=1,size(FIA_in%flux_sw_dn,3)
      call redistribute_data(domain_in, FIA_in%flux_sw_dn(:,:,b), domain_out, &
                             FIA_out%flux_sw_dn(:,:,b), complete=.false.)
    enddo
    call redistribute_data(domain_in, FIA_in%ice_free, domain_out, &
                           FIA_out%ice_free, complete=.false.)
    call redistribute_data(domain_in, FIA_in%ice_cover, domain_out, &
                           FIA_out%ice_cover, complete=.true.)

    if (allocated(FIA_in%flux_sh0)) then
      call redistribute_data(domain_in, FIA_in%flux_sh0, domain_out, &
                             FIA_out%flux_sh0, complete=.false.)
      call redistribute_data(domain_in, FIA_in%evap0, domain_out, &
                             FIA_out%evap0, complete=.false.)
      call redistribute_data(domain_in, FIA_in%flux_lw0, domain_out, &
                             FIA_out%flux_lw0, complete=.false.)
      call redistribute_data(domain_in, FIA_in%dshdt, domain_out, &
                             FIA_out%dshdt, complete=.false.)
      call redistribute_data(domain_in, FIA_in%devapdt, domain_out, &
                             FIA_out%devapdt, complete=.false.)
      call redistribute_data(domain_in, FIA_in%dlwdt, domain_out, &
                             FIA_out%dlwdt, complete=.true.)
      call redistribute_data(domain_in, FIA_in%Tskin_cat, domain_out, &
                             FIA_out%Tskin_cat, complete=.true.)
    endif

  elseif (associated(FIA_out)) then
    ! Use the null pointers in place of the unneeded input arrays.
    call coupler_type_redistribute_data(null_bc, domain_in, &
                          FIA_out%tr_flux, domain_out, complete=.false.)
    do b=1,size(FIA_out%flux_sw_top,4)
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             FIA_out%flux_sw_top(:,:,:,b), complete=.false.)
    enddo
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%flux_sh_top, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%evap_top, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%flux_lw_top, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%flux_lh_top, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%lprec_top, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%fprec_top, complete=.true.)

    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%tmelt, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%bmelt, complete=.false.)
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           FIA_out%sw_abs_ocn, complete=.true.)

    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%WindStr_x, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%WindStr_y, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%WindStr_ocn_x, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%WindStr_ocn_y, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%p_atm_surf, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%runoff, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%calving, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%runoff_hflx, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%calving_hflx, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%Tskin_avg, complete=.false.)
    do b=1,size(FIA_out%flux_sw_dn,3)
      call redistribute_data(domain_in, null_ptr2D, domain_out, &
                             FIA_out%flux_sw_dn(:,:,b), complete=.false.)
    enddo
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%ice_free, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           FIA_out%ice_cover, complete=.true.)

    if (allocated(FIA_out%flux_sh0)) then
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             FIA_out%flux_sh0, complete=.false.)
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             FIA_out%evap0, complete=.false.)
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             FIA_out%flux_lw0, complete=.false.)
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             FIA_out%dshdt, complete=.false.)
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             FIA_out%devapdt, complete=.false.)
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             FIA_out%dlwdt, complete=.true.)
      call redistribute_data(domain_in, null_ptr3D, domain_out, &
                             FIA_out%Tskin_cat, complete=.true.)
    endif


  elseif (associated(FIA_in)) then
    ! Use the null pointers in place of the unneeded output arrays.
    call coupler_type_redistribute_data(FIA_in%tr_flux, domain_in, &
                          null_bc, domain_out, complete=.false.)
    do b=1,size(FIA_in%flux_sw_top,4)
      call redistribute_data(domain_in, FIA_in%flux_sw_top(:,:,:,b), domain_out, &
                             null_ptr3D, complete=.false.)
    enddo
    call redistribute_data(domain_in, FIA_in%flux_sh_top, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%evap_top, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%flux_lw_top, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%flux_lh_top, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%lprec_top, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%fprec_top, domain_out, &
                           null_ptr3D, complete=.true.)

    call redistribute_data(domain_in, FIA_in%tmelt, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%bmelt, domain_out, &
                           null_ptr3D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%sw_abs_ocn, domain_out, &
                           null_ptr3D, complete=.true.)

    call redistribute_data(domain_in, FIA_in%WindStr_x, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%WindStr_y, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%WindStr_ocn_x, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%WindStr_ocn_y, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%p_atm_surf, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%runoff, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%calving, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%runoff_hflx, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%calving_hflx, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%Tskin_avg, domain_out, &
                           null_ptr2D, complete=.false.)
    do b=1,size(FIA_in%flux_sw_dn,3)
      call redistribute_data(domain_in, FIA_in%flux_sw_dn(:,:,b), domain_out, &
                             null_ptr2D, complete=.false.)
    enddo
    call redistribute_data(domain_in, FIA_in%ice_free, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, FIA_in%ice_cover, domain_out, &
                           null_ptr2D, complete=.true.)

    if (allocated(FIA_in%flux_sh0)) then
      call redistribute_data(domain_in, FIA_in%flux_sh0, domain_out, &
                             null_ptr3D, complete=.false.)
      call redistribute_data(domain_in, FIA_in%evap0, domain_out, &
                             null_ptr3D, complete=.false.)
      call redistribute_data(domain_in, FIA_in%flux_lw0, domain_out, &
                             null_ptr3D, complete=.false.)
      call redistribute_data(domain_in, FIA_in%dshdt, domain_out, &
                             null_ptr3D, complete=.false.)
      call redistribute_data(domain_in, FIA_in%devapdt, domain_out, &
                             null_ptr3D, complete=.false.)
      call redistribute_data(domain_in, FIA_in%dlwdt, domain_out, &
                             null_ptr3D, complete=.true.)
      call redistribute_data(domain_in, FIA_in%Tskin_cat, domain_out, &
                             null_ptr3D, complete=.true.)
    endif

  else
    call SIS_error(FATAL, "redistribute_FIA_to_FIA called with "//&
                          "neither FIA_in nor FIA_out associated.")
  endif

end subroutine redistribute_FIA_to_FIA

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> copy_TSF_to_TSF copies the computational domain of one fast_ice_avg_type into
!! the computational domain of another fast_ice_avg_type.  Both must use the same
!! domain decomposition and indexing convention, but they may have different
!! halo sizes.
subroutine copy_TSF_to_TSF(TSF_in, TSF_out, HI_in, HI_out)
  type(total_sfc_flux_type), intent(inout) :: TSF_in  !< The TSF type that is being copied from
  type(total_sfc_flux_type), intent(inout) :: TSF_out !< The TSF type that is being copied into
  type(hor_index_type),      intent(in)    :: HI_in  !< The horizontal index type for the input type
  type(hor_index_type),      intent(in)    :: HI_out !< The horizontal index type for the output type

  integer :: b, i, j, k, m, n, nb, isc, iec, jsc, jec, ncat
  integer :: i2, j2, i_off, j_off
  integer :: isd, ied, jsd, jed

  isc = HI_in%isc ; iec = HI_in%iec ; jsc = HI_in%jsc ; jec = HI_in%jec
  nb = size(TSF_in%flux_sw,3)

  if ((HI_in%iec-HI_in%isc /= HI_out%iec-HI_out%isc) .or. &
      (HI_in%jec-HI_in%jsc /= HI_out%jec-HI_out%jsc)) then
    call SIS_error(FATAL, "copy_TSF_to_TSF called with inconsistent domain "//&
                          "decompositions of the two ice types.")
  endif
  i_off = HI_out%iec-HI_in%iec ;  j_off = HI_out%jec-HI_in%jec

  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    TSF_out%flux_sh(i2,j2) = TSF_in%flux_sh(i,j)
    TSF_out%flux_lw(i2,j2) = TSF_in%flux_lw(i,j)
    TSF_out%flux_lh(i2,j2) = TSF_in%flux_lh(i,j)
    TSF_out%evap(i2,j2) = TSF_in%evap(i,j)
    TSF_out%lprec(i2,j2) = TSF_in%lprec(i,j)
    TSF_out%fprec(i2,j2) = TSF_in%fprec(i,j)
    do b=1,nb ; TSF_out%flux_sw(i2,j2,b) = TSF_in%flux_sw(i,j,b) ; enddo
  enddo ; enddo

  if (TSF_in%copy_calls /= TSF_out%copy_calls) call SIS_error(WARNING, &
    "copy_TSF_to_TSF called an inconsistent number of time for the input and output types.")
  TSF_in%copy_calls = TSF_in%copy_calls + 1 ; TSF_out%copy_calls = TSF_out%copy_calls + 1

  call coupler_type_copy_data(TSF_in%tr_flux, TSF_out%tr_flux)

end subroutine copy_TSF_to_TSF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> redistribute_TSF_to_TSF redistributes the computational domain of one
!! total_sfc_flux_type into the computational domain of another total_sfc_flux_type.
subroutine redistribute_TSF_to_TSF(TSF_in, TSF_out, domain_in, domain_out, HI_out)
  type(total_sfc_flux_type), pointer    :: TSF_in     !< The total_sfc_flux_type that is being copied from
  type(total_sfc_flux_type), pointer    :: TSF_out    !< The total_sfc_flux_type that is being copied into
  type(domain2d),            intent(in) :: domain_in  !< The source data domain.
  type(domain2d),            intent(in) :: domain_out !< The target data domain.
  type(hor_index_type), optional, intent(in) :: HI_out !< The hor_index_type on the target domain; HI_out
                                                       !! may be omitted if this is not a target PE.

  real, pointer, dimension(:,:) :: null_ptr2D => NULL()
  type(coupler_2d_bc_type)  :: null_bc
  integer :: copy_calls  ! The number of times these TSF_types have been copied.
  integer :: b, m, num_tr

  if (.not. (associated(TSF_out) .or. associated(TSF_in))) &
    call SIS_error(FATAL, "redistribute_TSF_to_TSF called with "//&
                          "neither TSF_in nor TSF_out associated.")
  copy_calls = 0
  if (associated(TSF_out)) then
    TSF_out%copy_calls = TSF_out%copy_calls + 1
    copy_calls = TSF_out%copy_calls
  endif
  if (associated(TSF_in)) then
    TSF_in%copy_calls = TSF_in%copy_calls + 1
    copy_calls = max(copy_calls,TSF_in%copy_calls)
  endif

  if (associated(TSF_out) .and. associated(TSF_in)) then
    ! The extra tracer arrays are copied first so that they can all have
    ! complete=.false.
    call coupler_type_redistribute_data(TSF_in%tr_flux, domain_in, &
                          TSF_out%tr_flux, domain_out, complete=.false.)
    do b=1,size(TSF_in%flux_sw,3)
      call redistribute_data(domain_in, TSF_in%flux_sw(:,:,b), domain_out, &
                             TSF_out%flux_sw(:,:,b), complete=.false.)
    enddo
    call redistribute_data(domain_in, TSF_in%flux_sh, domain_out, &
                           TSF_out%flux_sh, complete=.false.)
    call redistribute_data(domain_in, TSF_in%flux_lw, domain_out, &
                           TSF_out%flux_lw, complete=.false.)
    call redistribute_data(domain_in, TSF_in%flux_lh, domain_out, &
                           TSF_out%flux_lh, complete=.false.)
    call redistribute_data(domain_in, TSF_in%evap, domain_out, &
                           TSF_out%evap, complete=.false.)
    call redistribute_data(domain_in, TSF_in%lprec, domain_out, &
                           TSF_out%lprec, complete=.false.)
    call redistribute_data(domain_in, TSF_in%fprec, domain_out, &
                           TSF_out%fprec, complete=.true.)
  elseif (associated(TSF_out)) then
    ! Use the null pointer in place of the unneeded input arrays.
    call coupler_type_redistribute_data(null_bc, domain_in, &
                          TSF_out%tr_flux, domain_out, complete=.false.)
    do b=1,size(TSF_out%flux_sw,3)
      call redistribute_data(domain_in, null_ptr2D, domain_out, &
                             TSF_out%flux_sw(:,:,b), complete=.false.)
    enddo
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           TSF_out%flux_sh, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           TSF_out%flux_lw, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           TSF_out%flux_lh, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           TSF_out%evap, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           TSF_out%lprec, complete=.false.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           TSF_out%fprec, complete=.true.)
  elseif (associated(TSF_in)) then
    ! Use the null pointer in place of the unneeded output arrays.
    call coupler_type_redistribute_data(TSF_in%tr_flux, domain_in, &
                          null_bc, domain_out, complete=.false.)
    do b=1,size(TSF_in%flux_sw,3)
      call redistribute_data(domain_in, TSF_in%flux_sw(:,:,b), domain_out, &
                             null_ptr2D, complete=.false.)
    enddo
    call redistribute_data(domain_in, TSF_in%flux_sh, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, TSF_in%flux_lw, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, TSF_in%flux_lh, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, TSF_in%evap, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, TSF_in%lprec, domain_out, &
                           null_ptr2D, complete=.false.)
    call redistribute_data(domain_in, TSF_in%fprec, domain_out, &
                           null_ptr2D, complete=.true.)
  else
    call SIS_error(FATAL, "redistribute_TSF_to_TSF called with "//&
                          "neither TSF_in nor TSF_out associated.")
  endif


end subroutine redistribute_TSF_to_TSF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> copy_Rad_to_Rad copies the computational domains of several fields from one
!! ice_rad_type into the computational domain of another ice_rad_type.  Both must
!! use the same domain decomposition and indexing convention, but they may have
!! different halo sizes.
subroutine copy_Rad_to_Rad(Rad_in, Rad_out, HI_in, HI_out, IG)
  type(ice_rad_type),   intent(inout) :: Rad_in !< A structure with fields related to the absorption,
                                                !! reflection and transmission of shortwave radiation
                                                !! that is being copied from
  type(ice_rad_type),   intent(inout) :: Rad_out !< A structure with fields related to the absorption,
                                                !! reflection and transmission of shortwave radiation
                                                !! that is being copied into
  type(hor_index_type), intent(in)    :: HI_in  !< The horizontal index type for the input type
  type(hor_index_type), intent(in)    :: HI_out !< The horizontal index type for the output type
  type(ice_grid_type),  intent(in)    :: IG     !< The sea-ice specific grid type

  integer :: b, i, j, k, m, n, nb, isc, iec, jsc, jec, ncat, NkIce
  integer :: i2, j2, i_off, j_off

  isc = HI_in%isc ; iec = HI_in%iec ; jsc = HI_in%jsc ; jec = HI_in%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  if ((HI_in%iec-HI_in%isc /= HI_out%iec-HI_out%isc) .or. &
      (HI_in%jec-HI_in%jsc /= HI_out%jec-HI_out%jsc)) then
    call SIS_error(FATAL, "copy_Rad_to_Rad called with inconsistent domain "//&
                          "decompositions of the two ice types.")
  endif
  i_off = HI_out%iec-HI_in%iec ;  j_off = HI_out%jec-HI_in%jec

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    Rad_out%tskin_rad(i2,j2,k) = Rad_in%tskin_rad(i,j,k)
  enddo ; enddo ; enddo
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    Rad_out%coszen_lastrad(i2,j2) = Rad_in%coszen_lastrad(i,j)
  enddo ; enddo

end subroutine copy_Rad_to_Rad

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> redistribute_Rad_to_Rad redistributes the computational domain of several fields
!! from one ice_rad_type into the computational domain of another ice_rad_type.
subroutine redistribute_Rad_to_Rad(Rad_in, Rad_out, domain_in, domain_out)
  type(ice_rad_type), pointer    :: Rad_in     !< The ice_rad_type that is being copied from (intent in).
  type(ice_rad_type), pointer    :: Rad_out    !< The ice_rad_type that is being copied into (intent inout).
  type(domain2d),     intent(in) :: domain_in  !< The source data domain.
  type(domain2d),     intent(in) :: domain_out !< The target data domain.

  real, pointer, dimension(:,:,:) :: null_ptr3D => NULL()
  real, pointer, dimension(:,:) :: null_ptr2D => NULL()
  integer :: m

  if (associated(Rad_out) .and. associated(Rad_in)) then
    call redistribute_data(domain_in, Rad_in%tskin_rad, domain_out, &
                           Rad_out%tskin_rad, complete=.true.)
    call redistribute_data(domain_in, Rad_in%coszen_lastrad, domain_out, &
                           Rad_out%coszen_lastrad, complete=.true.)
  elseif (associated(Rad_out)) then
    ! Use the null pointers in place of the unneeded input arrays.
    call redistribute_data(domain_in, null_ptr3D, domain_out, &
                           Rad_out%tskin_rad, complete=.true.)
    call redistribute_data(domain_in, null_ptr2D, domain_out, &
                           Rad_out%coszen_lastrad, complete=.true.)
  elseif (associated(Rad_in)) then
    ! Use the null pointers in place of the unneeded output arrays.
    call redistribute_data(domain_in, Rad_in%tskin_rad, domain_out, &
                           null_ptr3D, complete=.true.)
    call redistribute_data(domain_in, Rad_in%coszen_lastrad, domain_out, &
                           null_ptr2D, complete=.true.)
  else
    call SIS_error(FATAL, "redistribute_Rad_to_Rad called with "//&
                          "neither Rad_in nor Rad_out associated.")
  endif

end subroutine redistribute_Rad_to_Rad

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> register_fast_to_slow_restarts registers all of the restart fields that are
!!   required in order to be sent from the fast state to the slow state when
!!   the model is restart.  These are the fields that would be copied via the
!!   subroutines copy_FIA_to_FIA, copy_TSF_to_TSF and copy_Rad_to_Rad, and it
!!   should be called from the fast ice processors when redo_fast_update is true.
subroutine register_fast_to_slow_restarts(FIA, Rad, TSF, mpp_domain, Ice_restart, restart_file)
  type(fast_ice_avg_type),   pointer     :: FIA     !< The fast ice model's fast_ice_avg_type
  type(ice_rad_type),        pointer     :: Rad     !< The fast ice model's ice_rad_type
  type(total_sfc_flux_type), pointer     :: TSF     !< The fast ice model's total_sfc_flux_type
  type(domain2d),            intent(in)  :: mpp_domain !< The mpp domain descriptor
  type(SIS_restart_CS),      pointer     :: Ice_restart !< The control structure for the ice restarts
  character(len=*),          intent(in)  :: restart_file !< The name and path to the restart file

! These fields are needed because the open-water fluxes are not recalculated.  It might be
! possible to make the fast-to-slow restart file smaller by breaking out the open-ocean
! category.
  call register_restart_field(Ice_restart, 'flux_sh_top', FIA%flux_sh_top, dim_3="cat0", &
                              mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'evap_top', FIA%evap_top, dim_3="cat0", &
                              mandatory=.false., units="kg m-2 s-1")
  call register_restart_field(Ice_restart, 'flux_lw_top', FIA%flux_lw_top, dim_3="cat0", &
                              mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'flux_lh_top', FIA%flux_lh_top, dim_3="cat0", &
                              mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'lprec_top', FIA%lprec_top, dim_3="cat0", &
                              mandatory=.false., units="kg m-2 s-1")
  call register_restart_field(Ice_restart, 'fprec_top', FIA%fprec_top, dim_3="cat0", &
                              mandatory=.false., units="kg m-2 s-1")
  call register_restart_field(Ice_restart, 'flux_sw_top', FIA%flux_sw_top, dim_3="cat0", &
                              dim_4="band", mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'WindStr_x', FIA%WindStr_x, &
                              mandatory=.false., units="Pa")
  call register_restart_field(Ice_restart, 'WindStr_y', FIA%WindStr_y, &
                              mandatory=.false., units="Pa")
  call register_restart_field(Ice_restart, 'WindStr_ocn_x', FIA%WindStr_ocn_x, &
                              mandatory=.false., units="Pa")
  call register_restart_field(Ice_restart, 'WindStr_ocn_y', FIA%WindStr_ocn_y, &
                              mandatory=.false., units="Pa")
  call register_restart_field(Ice_restart, 'p_atm_surf', FIA%p_atm_surf, &
                              mandatory=.false., units="Pa")
  call register_restart_field(Ice_restart, 'runoff', FIA%runoff, &
                              mandatory=.false., units="kg m-2 s-1")
  call register_restart_field(Ice_restart, 'calving', FIA%calving, &
                              mandatory=.false., units="kg m-2 s-1")
  call register_restart_field(Ice_restart, 'runoff_hflx', FIA%runoff_hflx, &
                              mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'calving_hflx', FIA%calving_hflx, &
                              mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'Tskin_avg', FIA%Tskin_avg, &
                              mandatory=.false., units="degC")
  call register_restart_field(Ice_restart, 'ice_free', FIA%ice_free, &
                              mandatory=.false., units="nondim")
  call register_restart_field(Ice_restart, 'ice_cover', FIA%ice_cover, &
                              mandatory=.false., units="nondim")
  call register_restart_field(Ice_restart, 'flux_sw_dn', FIA%flux_sw_dn, dim_3="band", &
                              mandatory=.false., units="W m-2")

  if (allocated(FIA%flux_sh0)) then
    call register_restart_field(Ice_restart, 'flux_sh_T0', FIA%flux_sh0, dim_3="cat0", &
                                mandatory=.false., units="W m-2")
    call register_restart_field(Ice_restart, 'flux_lw_T0', FIA%flux_lw0, dim_3="cat0", &
                                mandatory=.false., units="W m-2")
    call register_restart_field(Ice_restart, 'evap_T0', FIA%evap0, dim_3="cat0", &
                                mandatory=.false., units="kg m-2 s-1")
    call register_restart_field(Ice_restart, 'dsh_dT', FIA%dshdt, dim_3="cat0", &
                                mandatory=.false., units="W m-2 degC-1")
    call register_restart_field(Ice_restart, 'dlw_dT', FIA%dlwdt, dim_3="cat0", &
                                mandatory=.false., units="W m-2 degC-1")
    call register_restart_field(Ice_restart, 'devap_dT', FIA%devapdt, dim_3="cat0", &
                                mandatory=.false., units="kg m-2 s-1 degC-1")
    call register_restart_field(Ice_restart, 'Tskin_can', FIA%Tskin_cat, dim_3="cat0", &
                                mandatory=.false., units="degC")
  endif

  call register_restart_field(Ice_restart, 'tskin_rad', Rad%tskin_rad, &
                              mandatory=.false., units="degC")
  call register_restart_field(Ice_restart, 'coszen_rad', Rad%coszen_lastrad, &
                              mandatory=.false., units="nondim")

  call register_restart_field(Ice_restart, 'total_flux_sh', TSF%flux_sh, &
                              mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'total_flux_lw', TSF%flux_lw, &
                              mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'total_flux_lh', TSF%flux_lh, &
                              mandatory=.false., units="W m-2")
  call register_restart_field(Ice_restart, 'total_evap', TSF%evap, &
                              mandatory=.false., units="kg m-2 s-1")
  call register_restart_field(Ice_restart, 'total_lprec', TSF%lprec, &
                              mandatory=.false., units="kg m-2 s-1")
  call register_restart_field(Ice_restart, 'total_fprec', TSF%fprec, &
                              mandatory=.false., units="kg m-2 s-1")
  call register_restart_field(Ice_restart, 'total_flux_sw', TSF%flux_sw, dim_3="band", &
                              longname="Total shortwave flux by frequency and angular band", &
                              mandatory=.false., units="W m-2")

  if (coupler_type_initialized(TSF%tr_flux) .and. &
      coupler_type_initialized(FIA%tr_flux)) then
    call register_restart_field(Ice_restart, TSF%tr_flux, "TSF_")
    call register_restart_field(Ice_restart, FIA%tr_flux, "FIA_")
  endif

end subroutine register_fast_to_slow_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_IST_arrays deallocates the arrays in an ice_state_type.
subroutine dealloc_IST_arrays(IST)
  type(ice_state_type), intent(inout) :: IST !< A type describing the state of the sea ice

  deallocate(IST%part_size, IST%mH_snow, IST%mH_ice)
  deallocate(IST%mH_pond) ! mw/new
  deallocate(IST%enth_snow, IST%enth_ice, IST%sal_ice)
  if (allocated(IST%snow_to_ocn)) deallocate(IST%snow_to_ocn)
  if (allocated(IST%enth_snow_to_ocn)) deallocate(IST%enth_snow_to_ocn)
  if (allocated(IST%t_surf)) deallocate(IST%t_surf)

  if (allocated(IST%u_ice_C)) deallocate(IST%u_ice_C)
  if (allocated(IST%v_ice_C)) deallocate(IST%v_ice_C)
  if (allocated(IST%u_ice_B)) deallocate(IST%u_ice_B)
  if (allocated(IST%v_ice_B)) deallocate(IST%v_ice_B)

end subroutine dealloc_IST_arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_ocean_sfc_state deallocates the arrays in an ocean_sfc_state_type.
subroutine dealloc_ocean_sfc_state(OSS)
  type(ocean_sfc_state_type), pointer :: OSS !< A structure containing the arrays that describe
                                        !! the ocean's surface state that is deallocated here.


  if (.not.associated(OSS)) then
    call SIS_error(WARNING, "dealloc_ocean_sfc_state called with an unassociated pointer.")
    return
  endif

  deallocate(OSS%s_surf, OSS%SST_C, OSS%sea_lev, OSS%T_fr_ocn, OSS%frazil, OSS%bheat)
  if (allocated(OSS%u_ocn_B)) deallocate(OSS%u_ocn_B)
  if (allocated(OSS%v_ocn_B)) deallocate(OSS%v_ocn_B)
  if (allocated(OSS%u_ocn_C)) deallocate(OSS%u_ocn_C)
  if (allocated(OSS%v_ocn_C)) deallocate(OSS%v_ocn_C)

  deallocate(OSS)
end subroutine dealloc_ocean_sfc_state

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_simple_OSS deallocates the arrays in a simple_OSS_type.
subroutine dealloc_simple_OSS(OSS)
  type(simple_OSS_type), pointer :: OSS !< A structure containing the arrays that describe
                                        !! the ocean's surface state that is deallocated here.

  if (.not.associated(OSS)) then
    call SIS_error(WARNING, "dealloc_ocean_sfc_state called with an unassociated pointer.")
    return
  endif

  deallocate(OSS%s_surf, OSS%SST_C, OSS%bheat, OSS%T_fr_ocn)
  deallocate(OSS%u_ocn_A, OSS%v_ocn_A, OSS%u_ice_A, OSS%v_ice_A)

  deallocate(OSS)
end subroutine dealloc_simple_OSS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_fast_ice_avg deallocates the arrays in a fast_ice_avg_type.
subroutine dealloc_fast_ice_avg(FIA)
  type(fast_ice_avg_type), pointer    :: FIA !< A type containing averages of fields
                                             !! that is being deallocated here

  if (.not.associated(FIA)) then
    call SIS_error(WARNING, "dealloc_fast_ice_avg called with an unassociated pointer.")
    return
  endif

  deallocate(FIA%flux_u_top, FIA%flux_v_top )
  deallocate(FIA%flux_sh_top, FIA%evap_top, FIA%flux_lw_top)
  deallocate(FIA%flux_lh_top, FIA%lprec_top, FIA%fprec_top)
  deallocate(FIA%flux_sw_top)
  deallocate(FIA%runoff, FIA%calving, FIA%runoff_hflx, FIA%calving_hflx)
  deallocate(FIA%calving_preberg, FIA%calving_hflx_preberg)

  deallocate(FIA%tmelt, FIA%bmelt, FIA%frazil_left)
  deallocate(FIA%WindStr_x, FIA%WindStr_y, FIA%p_atm_surf)
  deallocate(FIA%WindStr_ocn_x, FIA%WindStr_ocn_y)
  deallocate(FIA%ice_free, FIA%ice_cover, FIA%sw_abs_ocn, FIA%Tskin_avg)

  if (allocated(FIA%flux_sh0)) deallocate(FIA%flux_sh0)
  if (allocated(FIA%evap0)) deallocate(FIA%evap0)
  if (allocated(FIA%flux_lw0)) deallocate(FIA%flux_lw0)
  if (allocated(FIA%dshdt))  deallocate(FIA%dshdt)
  if (allocated(FIA%devapdt))  deallocate(FIA%devapdt)
  if (allocated(FIA%dlwdt)) deallocate(FIA%dlwdt)
  if (allocated(FIA%Tskin_cat)) deallocate(FIA%Tskin_cat)

  deallocate(FIA)
end subroutine dealloc_fast_ice_avg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_total_sfc_flux deallocates the arrays in a total_sfc_flux_type.
subroutine dealloc_total_sfc_flux(TSF)
  type(total_sfc_flux_type), pointer :: TSF !< A type with averaged surface fluxes
                                            !! that is to be deallocated here.

  if (.not.associated(TSF)) then
    call SIS_error(WARNING, "dealloc_total_sfc_flux called with an unassociated pointer.")
    return
  endif

  deallocate(TSF%flux_u, TSF%flux_v, TSF%flux_sh, TSF%evap)
  deallocate(TSF%flux_sw)
  deallocate(TSF%flux_lw, TSF%flux_lh, TSF%lprec, TSF%fprec)

end subroutine dealloc_total_sfc_flux

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_ice_rad deallocates the arrays in a ice_rad_type.
subroutine dealloc_ice_rad(Rad)
  type(ice_rad_type), pointer :: Rad !< A structure with fields related to
                            !! the absorption, reflection and transmission of
                            !! shortwave radiation that is deallocated here.

  if (.not.associated(Rad)) then
    call SIS_error(WARNING, "dealloc_ice_rad called with an unassociated pointer.")
    return
  endif

  deallocate(Rad%sw_abs_sfc, Rad%sw_abs_snow, Rad%sw_abs_ice)
  deallocate(Rad%sw_abs_ocn, Rad%sw_abs_int)
  deallocate(Rad%coszen_nextrad, Rad%coszen_lastrad)
  deallocate(Rad%T_skin, Rad%Tskin_rad)

  deallocate(Rad)
end subroutine dealloc_ice_rad

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_ice_ocean_flux deallocates the arrays in a ice_ocean_flux_type.
subroutine dealloc_ice_ocean_flux(IOF)
  type(ice_ocean_flux_type), pointer :: IOF !< A structure containing fluxes from the ice to
                                            !! the ocean that is deallocated here.

  if (.not.associated(IOF)) then
    call SIS_error(WARNING, "dealloc_ice_ocean_flux called with an unassociated pointer.")
    return
  endif

  deallocate(IOF%flux_sh_ocn_top, IOF%evap_ocn_top)
  deallocate(IOF%flux_lw_ocn_top, IOF%flux_lh_ocn_top)
  deallocate(IOF%flux_sw_ocn)
  deallocate(IOF%lprec_ocn_top, IOF%fprec_ocn_top, IOF%flux_salt)
  deallocate(IOF%flux_u_ocn, IOF%flux_v_ocn, IOF%pres_ocn_top, IOF%mass_ice_sn_p)
  if (allocated(IOF%stress_mag)) deallocate(IOF%stress_mag)
  if (allocated(IOF%transmutation_salt_flux)) deallocate(IOF%transmutation_salt_flux)

  deallocate(IOF%Enth_Mass_in_atm, IOF%Enth_Mass_out_atm)
  deallocate(IOF%Enth_Mass_in_ocn, IOF%Enth_Mass_out_ocn)
  if (allocated(IOF%transmutation_enth)) deallocate(IOF%transmutation_enth)

  !Deallocating iceberg fields
  if (associated(IOF%mass_berg)) deallocate(IOF%mass_berg)
  if (associated(IOF%ustar_berg)) deallocate(IOF%ustar_berg)
  if (associated(IOF%area_berg)) deallocate(IOF%area_berg)

  deallocate(IOF)
end subroutine dealloc_ice_ocean_flux

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays in an ice_ocean_flux_type.
subroutine IOF_chksum(mesg, IOF, G, US, mech_fluxes, thermo_fluxes)
  character(len=*),          intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(ice_ocean_flux_type), intent(in) :: IOF   !< The structure whose arrays are being checksummed.
  type(SIS_hor_grid_type),   intent(inout) :: G  !< The ice-model's horizontal grid type.
  type(unit_scale_type),     intent(in)    :: US !< A structure with unit conversion factors
  logical,         optional, intent(in)    :: mech_fluxes !< If true, do checksums of mechanical fluxes
  logical,         optional, intent(in)    :: thermo_fluxes !< If true, do checksums of thermodynamic fluxes

  logical :: do_mech, do_thermo

  ! Do all fluxes unless a subset of fluxes are specified, in which case only do those that are
  ! indicated by the arguments.
  do_mech = .not.present(thermo_fluxes) ; do_thermo = .not.present(mech_fluxes)
  if (present(mech_fluxes)) then ; do_mech = mech_fluxes ; endif
  if (present(thermo_fluxes)) then ; do_thermo = thermo_fluxes ; endif

  if (do_thermo) then
    call hchksum(IOF%flux_salt, trim(mesg)//" IOF%flux_salt", G%HI, scale=US%RZ_T_to_kg_m2s)
    if (allocated(IOF%transmutation_salt_flux)) call hchksum(IOF%transmutation_salt_flux, &
          trim(mesg)//" IOF%transmutation_salt_flux", G%HI, scale=US%RZ_T_to_kg_m2s)

    call hchksum(IOF%flux_sh_ocn_top, trim(mesg)//" IOF%flux_sh_ocn_top", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(IOF%flux_lw_ocn_top, trim(mesg)//" IOF%flux_lw_ocn_top", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(IOF%flux_lh_ocn_top, trim(mesg)//" IOF%flux_lh_ocn_top", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(IOF%flux_sw_ocn,     trim(mesg)//" IOF%flux_sw_ocn",     G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(IOF%evap_ocn_top,    trim(mesg)//" IOF%evap_ocn_top",  G%HI, scale=US%RZ_T_to_kg_m2s)
    call hchksum(IOF%lprec_ocn_top,   trim(mesg)//" IOF%lprec_ocn_top", G%HI, scale=US%RZ_T_to_kg_m2s)
    call hchksum(IOF%fprec_ocn_top,   trim(mesg)//" IOF%fprec_ocn_top", G%HI, scale=US%RZ_T_to_kg_m2s)
  endif
  if (do_mech) then
    call hchksum(IOF%flux_u_ocn,      trim(mesg)//" IOF%flux_u_ocn",   G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    call hchksum(IOF%flux_v_ocn,      trim(mesg)//" IOF%flux_v_ocn",   G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    call hchksum(IOF%pres_ocn_top,    trim(mesg)//" IOF%pres_ocn_top", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    call hchksum(IOF%mass_ice_sn_p,   trim(mesg)//" IOF%mass_ice_sn_p", G%HI, scale=US%RZ_to_kg_m2)
    if (allocated(IOF%stress_mag)) &
      call hchksum(IOF%stress_mag,    trim(mesg)//" IOF%stress_mag", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
  endif

  if (do_thermo) then
    call hchksum(IOF%Enth_Mass_in_atm,  trim(mesg)//" IOF%Enth_Mass_in_atm",  G%HI, scale=US%QRZ_T_to_W_m2*US%T_to_s)
    call hchksum(IOF%Enth_Mass_out_atm, trim(mesg)//" IOF%Enth_Mass_out_atm", G%HI, scale=US%QRZ_T_to_W_m2*US%T_to_s)
    call hchksum(IOF%Enth_Mass_in_ocn,  trim(mesg)//" IOF%Enth_Mass_in_ocn",  G%HI, scale=US%QRZ_T_to_W_m2*US%T_to_s)
    call hchksum(IOF%Enth_Mass_out_ocn, trim(mesg)//" IOF%Enth_Mass_out_ocn", G%HI, scale=US%QRZ_T_to_W_m2*US%T_to_s)

    if (allocated(IOF%transmutation_enth)) call hchksum(IOF%transmutation_enth, &
          trim(mesg)//" IOF%transmutation_enth", G%HI, scale=US%QRZ_T_to_W_m2*US%T_to_s)
  endif
end subroutine IOF_chksum

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays in a fast_ice_avg_type.
subroutine FIA_chksum(mesg, FIA, G, US, check_ocean)
  character(len=*),        intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(fast_ice_avg_type), intent(in) :: FIA   !< The structure whose arrays are being checksummed.
  type(SIS_hor_grid_type), intent(inout) :: G  !< The ice-model's horizontal grid type.
  type(unit_scale_type),   intent(in)    :: US !< A structure with unit conversion factors
  logical, optional,       intent(in) :: check_ocean !< If present and true, check the fluxes to the ocean.

  character(len=8) :: nstr
  integer :: b

  call hchksum(FIA%flux_sh_top(:,:,1:), trim(mesg)//" FIA%flux_sh_top", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%evap_top(:,:,1:), trim(mesg)//" FIA%evap_top", G%HI, scale=US%RZ_T_to_kg_m2s)
  do b=1,size(FIA%flux_sw_top,4)
    write(nstr, '(I4)') b ; nstr = adjustl(nstr)
    call hchksum(FIA%flux_sw_top(:,:,1:,b), &
                 trim(mesg)//" FIA%flux_sw_top("//trim(nstr)//")", G%HI, scale=US%QRZ_T_to_W_m2)
  enddo
  call hchksum(FIA%flux_lw_top(:,:,1:), trim(mesg)//" FIA%flux_lw_top", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%flux_lh_top(:,:,1:), trim(mesg)//" FIA%flux_lh_top", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%lprec_top(:,:,1:), trim(mesg)//" FIA%lprec_top", G%HI, scale=US%RZ_T_to_kg_m2s)
  call hchksum(FIA%fprec_top(:,:,1:), trim(mesg)//" FIA%fprec_top", G%HI, scale=US%RZ_T_to_kg_m2s)

  if (present(check_ocean)) then ; if (check_ocean) then
    call hchksum(FIA%flux_sh_top(:,:,0), trim(mesg)//" FIA%flux_sh_top0", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(FIA%evap_top(:,:,0), trim(mesg)//" FIA%evap_top0", G%HI, scale=US%RZ_T_to_kg_m2s)
    do b=1,size(FIA%flux_sw_top,4)
      write(nstr, '(I4)') b ; nstr = adjustl(nstr)
      call hchksum(FIA%flux_sw_top(:,:,0,b), &
                   trim(mesg)//" FIA%flux_sw_top0("//trim(nstr)//")", G%HI, scale=US%QRZ_T_to_W_m2)
    enddo
    call hchksum(FIA%flux_lw_top(:,:,0), trim(mesg)//" FIA%flux_lw_top0", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(FIA%flux_lh_top(:,:,0), trim(mesg)//" FIA%flux_lh_top0", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(FIA%lprec_top(:,:,0), trim(mesg)//" FIA%lprec_top0", G%HI, scale=US%RZ_T_to_kg_m2s)
    call hchksum(FIA%fprec_top(:,:,0), trim(mesg)//" FIA%fprec_top0", G%HI, scale=US%RZ_T_to_kg_m2s)
  endif ; endif

  call hchksum(FIA%tmelt, trim(mesg)//" FIA%tmelt", G%HI, scale=US%QRZ_T_to_W_m2*US%T_to_s)
  call hchksum(FIA%bmelt, trim(mesg)//" FIA%bmelt", G%HI, scale=US%QRZ_T_to_W_m2*US%T_to_s)
  call hchksum(FIA%sw_abs_ocn, trim(mesg)//" FIA%sw_abs_ocn", G%HI)

  call hchksum(FIA%WindStr_x, trim(mesg)//" FIA%WindStr_x", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
  call hchksum(FIA%WindStr_y, trim(mesg)//" FIA%WindStr_y", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
  call hchksum(FIA%WindStr_ocn_x, trim(mesg)//" FIA%WindStr_ocn_x", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
  call hchksum(FIA%WindStr_ocn_y, trim(mesg)//" FIA%WindStr_ocn_y", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
  call hchksum(FIA%p_atm_surf, trim(mesg)//" FIA%p_atm_surf", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
  call hchksum(FIA%runoff, trim(mesg)//" FIA%runoff", G%HI, scale=US%RZ_T_to_kg_m2s)
  call hchksum(FIA%calving, trim(mesg)//" FIA%calving", G%HI, scale=US%RZ_T_to_kg_m2s)
  call hchksum(FIA%runoff_hflx, trim(mesg)//" FIA%runoff_hflx", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%calving_hflx, trim(mesg)//" FIA%calving_hflx", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%ice_free, trim(mesg)//" FIA%ice_free", G%HI)
  call hchksum(FIA%ice_cover, trim(mesg)//" FIA%ice_cover", G%HI)
  call hchksum(FIA%flux_sw_dn, trim(mesg)//" FIA%flux_sw_dn", G%HI, scale=US%QRZ_T_to_W_m2)

end subroutine FIA_chksum


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays in an ocean_surface_state_type.
subroutine OSS_chksum(mesg, OSS, G, US, haloshift)
  character(len=*),           intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(ocean_sfc_state_type), intent(in) :: OSS   !< A structure containing the arrays that describe
                                                  !! the ocean's surface state for the ice model.
  type(SIS_hor_grid_type), intent(inout) :: G     !< The ice-model's horizontal grid type.
  type(unit_scale_type),      intent(in) :: US    !< A structure with unit conversion factors
  integer,          optional, intent(in) :: haloshift !< The width of halos to check, or 0 if missing.

  ! Local variables
  integer :: hs

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  hs=0 ; if (present(haloshift)) hs=haloshift

  call hchksum(OSS%s_surf, trim(mesg)//" OSS%s_surf", G%HI, haloshift=hs)
  call hchksum(OSS%SST_C, trim(mesg)//" OSS%SST_C", G%HI, haloshift=hs)
  call hchksum(OSS%T_fr_ocn, trim(mesg)//" OSS%T_fr_ocn", G%HI, haloshift=hs)
  call hchksum(OSS%sea_lev, trim(mesg)//" OSS%sea_lev", G%HI, haloshift=hs, scale=US%Z_to_m)
  call hchksum(OSS%bheat, trim(mesg)//" OSS%bheat", G%HI, haloshift=hs, scale=US%QRZ_T_to_W_m2)
  call hchksum(OSS%frazil, trim(mesg)//" OSS%frazil", G%HI, haloshift=hs, scale=US%QRZ_T_to_W_m2*US%T_to_s)

  if (OSS%Cgrid_dyn) then
    call uvchksum(mesg//" OSS%[uv]_ocn_C", OSS%u_ocn_C, OSS%v_ocn_C, G, halos=hs, scale=US%L_T_to_m_s)
    call check_redundant_C(mesg//" OSS%u/v_ocn_C", OSS%u_ocn_C, OSS%v_ocn_C, G)
  else
    call Bchksum_pair(mesg//" OSS%[uv]_ocn_B", OSS%u_ocn_B, OSS%v_ocn_B, G, halos=hs, scale=US%L_T_to_m_s)
    call check_redundant_B(mesg//" OSS%u/v_ocn", OSS%u_ocn_B, OSS%v_ocn_B, G)
  endif

end subroutine OSS_chksum



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays in an total_sfc_flux_type.
subroutine TSF_chksum(mesg, TSF, G, US, haloshift)
  character(len=*),           intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(total_sfc_flux_type),  pointer    :: TSF   !< The total_sfc_flux_type being checksummed.
  type(SIS_hor_grid_type), intent(inout) :: G     !< The ice-model's horizontal grid type.
  type(unit_scale_type),      intent(in) :: US    !< A structure with unit conversion factors
  integer,          optional, intent(in) :: haloshift !< The width of halos to check, or 0 if missing.

  ! Local variables
  integer :: hs, nb

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  hs=0 ; if (present(haloshift)) hs=haloshift

  call hchksum(TSF%flux_sh, trim(mesg)//" TSF%flux_sh", G%HI, haloshift=hs, scale=US%QRZ_T_to_W_m2)
  ! call hchksum(TSF%flux_sw, trim(mesg)//" TSF%flux_sw", G%HI, haloshift=hs, scale=US%QRZ_T_to_W_m2)
  call hchksum(TSF%flux_lw, trim(mesg)//" TSF%flux_lw", G%HI, haloshift=hs, scale=US%QRZ_T_to_W_m2)
  call hchksum(TSF%flux_lh, trim(mesg)//" TSF%flux_lh", G%HI, haloshift=hs, scale=US%QRZ_T_to_W_m2)
  call hchksum(TSF%evap, trim(mesg)//" TSF%evap", G%HI, haloshift=hs, scale=US%RZ_T_to_kg_m2s)
  call hchksum(TSF%lprec, trim(mesg)//" TSF%lprec", G%HI, haloshift=hs, scale=US%RZ_T_to_kg_m2s)
  call hchksum(TSF%fprec, trim(mesg)//" TSF%fprec", G%HI, haloshift=hs, scale=US%RZ_T_to_kg_m2s)
  call uvchksum(mesg//" TSF%flux_[uv]", TSF%flux_u, TSF%flux_v, G, halos=hs, scale=US%L_T_to_m_s)

end subroutine TSF_chksum


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays in an ice_state_type.
subroutine IST_chksum(mesg, IST, G, US, IG, haloshift)
  character(len=*),        intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(ice_state_type),    intent(in) :: IST   !< The structure whose arrays are being checksummed.
  type(SIS_hor_grid_type), intent(inout) :: G  !< The ice-model's horizontal grid type.
  type(unit_scale_type),   intent(in)    :: US !< A structure with unit conversion factors
  type(ice_grid_type),     intent(in) :: IG    !< The sea-ice grid type.
  integer, optional,       intent(in) :: haloshift !< The width of halos to check, or 0 if missing.

  ! Local variables
  character(len=20) :: k_str1, k_str
  integer :: hs, k

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  hs=0; if (present(haloshift)) hs=haloshift

  call hchksum(IST%part_size(:,:,0), trim(mesg)//" IST%part_size(0)", G%HI, haloshift=hs)
  call hchksum(IST%part_size(:,:,1:), trim(mesg)//" IST%part_size", G%HI, haloshift=hs)
  call hchksum(IST%mH_ice, trim(mesg)//" IST%mH_ice", G%HI, haloshift=hs, scale=US%RZ_to_kg_m2)
  do k=1,IG%NkIce
    write(k_str1,'(I8)') k ;  k_str = "("//trim(adjustl(k_str1))//")"
    call hchksum(IST%enth_ice(:,:,:,k), trim(mesg)//" IST%enth_ice("//trim(k_str), G%HI, &
                 haloshift=hs, scale=US%Q_to_J_kg)
    call hchksum(IST%sal_ice(:,:,:,k), trim(mesg)//" IST%sal_ice("//trim(k_str), G%HI, haloshift=hs)
  enddo
  call hchksum(IST%mH_snow, trim(mesg)//" IST%mH_snow", G%HI, haloshift=hs, scale=US%RZ_to_kg_m2)
  call hchksum(IST%enth_snow(:,:,:,1), trim(mesg)//" IST%enth_snow", G%HI, haloshift=hs, scale=US%Q_to_J_kg)
  call hchksum(IST%mH_pond, trim(mesg)//" IST%mH_pond", G%HI, haloshift=hs, scale=US%RZ_to_kg_m2)

  if (allocated(IST%u_ice_B) .and. allocated(IST%v_ice_B)) then
    call Bchksum_pair(mesg//" IST%[uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G, halos=hs, scale=US%L_T_to_m_s)
    call check_redundant_B(mesg//" IST%u/v_ice", IST%u_ice_B, IST%v_ice_B, G)
  endif
  if (allocated(IST%u_ice_C) .and. allocated(IST%v_ice_C)) then
    call uvchksum(mesg//" IST%[uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G, halos=hs, scale=US%L_T_to_m_s)
    call check_redundant_C(mesg//" IST%u/v_ice_C", IST%u_ice_C, IST%v_ice_C, G)
  endif

end subroutine IST_chksum

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Check the ice state for bad values of temperature, thickness, areal coverage,
!! enthalpy or ice mass, and write diagnostics about any offending columns
subroutine IST_bounds_check(IST, G, US, IG, msg, OSS, Rad)
  type(ice_state_type),    intent(in)    :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  character(len=*),        intent(in)    :: msg !< An identifying message
  type(ocean_sfc_state_type), optional, intent(in) :: OSS !< A structure containing the arrays that describe
                                                !! the ocean's surface state for the ice model.
  type(ice_rad_type),         optional, intent(in) :: Rad !< A structure with fields related to the
                                                !! absorption, reflection and transmission of shortwave radiation.

  character(len=512) :: mesg1, mesg2
  character(len=24) :: err
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: sum_part_sz
  real, dimension(IG%NkIce) :: S_col
  real    :: tsurf_min, tsurf_max, tice_min, tice_max, tOcn_min, tOcn_max
  real    :: enth_min, enth_max
  real    :: m_max ! Maximum mass per unit area [R Z ~> kg m-2]
  logical :: spec_thermo_sal
  integer :: i, j, k, m, isc, iec, jsc, jec, ncat, NkIce, i_off, j_off
  integer :: n_bad, i_bad, j_bad, k_bad

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  n_bad = 0 ; i_bad = 0 ; j_bad = 0 ; k_bad = 0 ; err = ":"

  m_max = 1.0e6*US%kg_m3_to_R*US%m_to_Z

  sum_part_sz(:,:) = 0.0
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    sum_part_sz(i,j) = sum_part_sz(i,j) + IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  tOcn_min = -100. ; tOcn_max = 60.
  if (present(OSS)) then
    do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
      if ((OSS%s_surf(i,j) < 0.0) .or. (OSS%s_surf(i,j) > 100.0) .or. &
          (OSS%SST_C(i,j) < tOcn_min) .or. (OSS%SST_C(i,j) > tOcn_max)) then
        n_bad = n_bad + 1
        if (n_bad == 1) then ; i_bad = i ; j_bad = j ; err = "t_ocn" ; endif
      endif
    endif ; enddo ; enddo
  endif
  do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
    if (abs(sum_part_sz(i,j) - 1.0) > 2.0*(ncat+1)*epsilon(sum_part_sz(i,j))) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; err = "sum_part_sz" ; endif
    endif
  endif ; enddo ; enddo

  tsurf_min = tOcn_min ; tsurf_max = tOcn_max
  tice_min = -100. ; tice_max = 1.0
  enth_min = enth_from_TS(tice_min, 0., IST%ITV)
  enth_max = enth_from_TS(tice_max, 0., IST%ITV)
  if (present(Rad)) then
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
      if ((Rad%t_skin(i,j,k) < tsurf_min) .or. (Rad%t_skin(i,j,k) > tsurf_max)) then
        n_bad = n_bad + 1
        if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "tsurf" ; endif
      endif
    endif ; enddo ; enddo ; enddo
  endif

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
    if ((IST%mH_ice(i,j,k) > m_max) .or. (IST%mH_snow(i,j,k) > m_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "large mass" ; endif
    endif
    if ((IST%enth_snow(i,j,k,1) < enth_min) .or. (IST%enth_snow(i,j,k,1) > enth_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "enth_snow" ; endif
    endif
  endif ; enddo ; enddo ; enddo

  do m=1,NkIce ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
    if ((IST%enth_ice(i,j,k,m) < enth_min) .or. (IST%enth_ice(i,j,k,m) > enth_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "enth_ice" ; endif
    endif
    if ((IST%sal_ice(i,j,k,m) < 0.0) .or. (IST%sal_ice(i,j,k,m) > 1000.0)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "sal_ice" ; endif
    endif
  endif ; enddo ; enddo ; enddo ; enddo

  if (n_bad > 0) then
    i = i_bad ; j=j_bad ; k = k_bad
    write(mesg1,'(" at ", 2(F6.1)," or i,j,k = ",3i4,"; nbad = ",i6," on pe ",i4)') &
           G%geolonT(i,j), G%geolatT(i,j), i_bad, j_bad, k_bad, n_bad, pe_here()
    if (k_bad > 0) then
      if (present(Rad)) then
        write(mesg2,'("T_skin = ",1pe12.4,", ps = ",1pe12.4)') Rad%t_skin(i,j,k), IST%part_size(i,j,k)
      else
        write(mesg2,'("part_size = ",1pe12.4)') IST%part_size(i,j,k)
      endif
    elseif (present(OSS)) then
      if (sum_part_sz(i,j) < 0.9999) then
        write(mesg2,'("T_ocn = ",1pe12.4,", S_sfc = ",1pe12.4,", sum_ps = ",1pe12.4)') &
              OSS%SST_C(i,j), OSS%s_surf(i,j), sum_part_sz(i,j)
      else
        write(mesg2,'("T_ocn = ",1pe12.4,", S_sfc = ",1pe12.4,", sum_ps = 1 - ",1pe12.4)') &
              OSS%SST_C(i,j), OSS%s_surf(i,j), 1.0-sum_part_sz(i,j)
      endif
    else
      if (sum_part_sz(i,j) < 0.9999) then
        write(mesg2,'("sum_part_sz = ",1pe12.4)') sum_part_sz(i,j)
      else
        write(mesg2,'("sum_part_sz = 1 - ",1pe12.4)') 1.0-sum_part_sz(i,j)
      endif
    endif
    call SIS_error(WARNING, "Bad ice state "//trim(err)//" "//trim(msg)//" ; "//trim(mesg1)//&
                            " ; "//trim(mesg2), all_print=.true.)
    if (k_bad > 0) then
      call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, spec_thermo_salin=spec_thermo_sal)
      if (.not.spec_thermo_sal) then
        do m=1,NkIce ; S_col(m) = IST%sal_ice(i,j,k,m) ; enddo
      endif
      write(mesg1,'("mi/ms = ", 2(1pe12.4)," ts = ",1pe12.4," ti = ",1pe12.4)') &
             IST%mH_ice(i,j,k)*US%RZ_to_kg_m2, IST%mH_snow(i,j,k)*US%RZ_to_kg_m2, &
             temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV), &
             temp_from_En_S(IST%enth_ice(i,j,k,1), S_col(1), IST%ITV)
      do m=2,NkIce
        write(mesg2,'(", ", 1pe12.4)') temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV)
        mesg1 = trim(mesg1)//trim(mesg2)
      enddo
      call SIS_error(WARNING, mesg1, all_print=.true.)
      write(mesg1,'("enth_snow = ",1pe12.4," enth_ice = ",1pe12.4)') &
             IST%enth_snow(i,j,k,1), IST%enth_ice(i,j,k,1)
      do m=2,NkIce
        write(mesg2,'(", ", 1pe12.4)') IST%enth_ice(i,j,k,m)
        mesg1 = trim(mesg1)//trim(mesg2)
      enddo
      call SIS_error(WARNING, mesg1, all_print=.true.)
      write(mesg1,'("salin_ice = ",1pe12.4)') IST%sal_ice(i,j,k,1)
      do m=2,NkIce
        write(mesg2,'(", ", 1pe12.4)') IST%sal_ice(i,j,k,m)
        mesg1 = trim(mesg1)//trim(mesg2)
      enddo
      call SIS_error(WARNING, mesg1, all_print=.true.)
    endif
  endif

end subroutine IST_bounds_check

end module SIS_types
