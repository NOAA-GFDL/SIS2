!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_types contains a number of common SIS types, along with subroutines to   !
!   perform various tasks on these types, including allocation, deallocation,  !
!   registration for restarts, and checksums.                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_types

! use mpp_mod,          only: mpp_sum, stdout, input_nml_file, PE_here => mpp_pe
! use mpp_domains_mod,  only: domain2D, mpp_get_compute_domain, CORNER, EAST, NORTH
use mpp_domains_mod,  only: domain2D, CORNER, EAST, NORTH, mpp_redistribute
! use mpp_parameter_mod, only: CGRID_NE, BGRID_NE, AGRID
! use fms_mod,          only: open_namelist_file, check_nml_error, close_file
! use fms_io_mod,       only: save_restart, restore_state, query_initialized
use fms_io_mod,       only: register_restart_field, restart_file_type
use time_manager_mod, only: time_type, time_type_to_real
use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type

use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type

use SIS2_ice_thm, only : ice_thermo_type, SIS2_ice_thm_CS, enth_from_TS, energy_melt_EnthS
use SIS2_ice_thm, only : get_SIS2_thermo_coefs, temp_from_En_S

use MOM_coms, only : PE_here, max_across_PEs
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg, is_root_pe
use MOM_file_parser, only : param_file_type
use MOM_hor_index,   only : hor_index_type
use SIS_diag_mediator, only : SIS_diag_ctrl, post_data=>post_SIS_data
use SIS_diag_mediator, only : register_SIS_diag_field, register_static_field
use SIS_debugging,   only : chksum, Bchksum, hchksum, uchksum, vchksum
use SIS_debugging,   only : check_redundant_B, check_redundant_C
use SIS_sum_output_type, only : SIS_sum_out_CS
use SIS_tracer_registry, only : SIS_tracer_registry_type

implicit none ; private

#include <SIS2_memory.h>

public :: ice_state_type, alloc_IST_arrays, ice_state_register_restarts, dealloc_IST_arrays
public :: IST_chksum, IST_bounds_check, copy_IST_to_IST
public :: ice_ocean_flux_type, alloc_ice_ocean_flux, dealloc_ice_ocean_flux
public :: ocean_sfc_state_type, alloc_ocean_sfc_state, dealloc_ocean_sfc_state
public :: fast_ice_avg_type, alloc_fast_ice_avg, dealloc_fast_ice_avg, copy_FIA_to_FIA
public :: IOF_chksum, FIA_chksum
public :: ice_rad_type, ice_rad_register_restarts, dealloc_ice_rad
public :: simple_OSS_type, alloc_simple_OSS, dealloc_simple_OSS, copy_sOSS_to_sOSS
public :: redistribute_IST_to_IST, redistribute_FIA_to_FIA, redistribute_sOSS_to_sOSS
public :: total_sfc_flux_type, alloc_total_sfc_flux, dealloc_total_sfc_flux
public :: translate_OSS_to_sOSS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! This structure contains the ice model state, and is intended to be private   !
! to SIS2.  It is not to be shared with other components and modules, and may  !
! use different indexing conventions than other components.                    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type ice_state_type
  ! The 8 of the following 10 variables constitute the sea-ice state.
  real, allocatable, dimension(:,:,:) :: &
    part_size   ! The fractional coverage of a grid cell by each ice
                ! thickness category, nondim, 0 to 1.  Category 0 is
                ! open ocean.  The sum of part_size is 1.
  ! These velocities are only used on the slow ice processors
  real, allocatable, dimension(:,:) :: &
    u_ice_B, &  ! The pseudo-zonal and pseudo-meridional ice velocities
    v_ice_B, &  ! along the model's grid directions on a B-grid, in m s-1.
                ! All thickness categories are assumed to have the same
                ! velocity.
    u_ice_C, &  ! The pseudo-zonal and pseudo-meridional ice velocities
    v_ice_C     ! along the model's grid directions on a C-grid, in m s-1.
                ! All thickness categories are assumed to have the same
                ! velocity.
  real, allocatable, dimension(:,:,:) :: &
    mH_pond, &  ! The mass per unit area of the pond in each category,
                ! in units of H (usually kg m-2). mw/new
    mH_snow, &  ! The mass per unit area of the snow in each category,
                ! in units of H (usually kg m-2).
    mH_ice, &   ! The mass per unit area of the ice in each category,
                ! in units of H (usually kg m-2).
    t_surf      ! The surface temperature, in Kelvin.

  real, allocatable, dimension(:,:,:,:) :: &
    sal_ice, &  ! The salinity of the sea ice in each category and
                ! fractional thickness layer, in g/kg.
    enth_ice, & ! The enthalpy of the sea ice in each category and
                ! fractional thickness layer, in enth_unit (J/kg or rescaled).
    enth_snow   ! The enthalpy of the snow in each category, in enth_unit.

  real, allocatable, dimension(:,:,:) :: &
    rdg_mice    ! A diagnostic of the ice load that was formed by
                ! ridging, in H (usually kg m-2).

  logical :: Cgrid_dyn ! If true use a C-grid discretization of the
                       ! sea-ice dynamics.

  type(SIS_tracer_registry_type), pointer :: TrReg => NULL()

  type(ice_thermo_type), pointer  :: ITV => NULL()
end type ice_state_type

!> ocean_sfc_state_type contains variables that describe the ocean's surface
!! state as seen by the slowly evolving sea-ice, on the ice grid.
type ocean_sfc_state_type
  ! 6 of the following 8 variables describe the ocean state as seen by the sea ice.
  real, allocatable, dimension(:,:) :: &
    s_surf , &  ! The ocean's surface salinity in g/kg.
    SST_C  , &  ! The ocean's bulk surface temperature in degC.
    T_fr_ocn, & ! The freezing point temperature in degC at the ocean's surface salinity.
    u_ocn_B, &  ! The ocean's zonal velocity on B-grid points in m s-1.
    v_ocn_B, &  ! The ocean's meridional velocity on B-grid points in m s-1.
    u_ocn_C, &  ! The ocean's zonal and meridional velocity on C-grid
    v_ocn_C, &  ! points, both in m s-1.
    frazil , &  ! A downward heat flux from the ice into the ocean
                ! associated with the formation of frazil ice in
                ! the ocean integrated over a timestep, in J m-2.
                ! This is the input value and is not changed by the ice.
    sea_lev     ! The equivalent sea-level, after any non-levitating
                ! ice has been converted to sea-water, as determined
                ! by the ocean, in m.  Sea-ice only contributes by
                ! applying pressure to the ocean that is then
                ! (partially) converted back to its equivalent by the
                ! ocean.

  real, allocatable, dimension(:,:,:) :: &
    tr_array    ! An array of fields related to properties for additional tracers.

  integer :: num_tr = -1  ! The number of additional tracer-related arrays.
!   type(coupler_3d_bc_type)   :: ocean_fields       ! array of fields used for additional tracers

  real :: kmelt ! A constant that is used in the calculation of the ocean/ice
                ! basal heat flux, in W m-2 K-1.  This could be replaced with
                ! an array reflecting the turbulence in the under-ice ocean
                ! boundary layer and the effective depth of the reported value
                ! of t_ocn.

  logical :: Cgrid_dyn ! If true use a C-grid discretization of the
                       ! sea-ice dynamics.

  ! diagnostic IDs for ocean surface  properties.
  integer :: id_sst=-1, id_sss=-1, id_ssh=-1, id_uo=-1, id_vo=-1, id_frazil=-1
end type ocean_sfc_state_type

!> simple_OSS_type contains variables that describe the ocean's surface
!! state as seen by the fast sea-ice or atmosphere, on the ice grid.
type simple_OSS_type
  ! The following 5 variables describe the ocean state as seen by the
  ! atmosphere and use for the rapid thermodynamic sea ice changes.
  real, allocatable, dimension(:,:) :: &
    s_surf , &  ! The ocean's surface salinity in g/kg.
    SST_C  , &  ! The ocean's bulk surface temperature in degC.
    T_fr_ocn, & ! The freezing point temperature in degC at the ocean's surface salinity.
    u_ocn_A, &  ! The ocean's zonal surface velocity on A-grid points in m s-1.
    v_ocn_A, &  ! The ocean's meridional surface velocity on A-grid points in m s-1.
    u_ice_A, &  ! The sea ice's zonal velocity on A-grid points in m s-1.
    v_ice_A, &  ! The sea ice's meridional velocity on A-grid points in m s-1.
    bheat       ! The upward diffusive heat flux from the ocean
                ! to the ice at the base of the ice, in W m-2.

  real, allocatable, dimension(:,:,:) :: &
    tr_array    ! An array of fields related to properties for additional tracers.
  integer :: num_tr = -1  ! The number of additional tracer-related arrays.
  logical :: first_copy = .true.

end type simple_OSS_type


!> fast_ice_avg_type contains variables that describe the fluxes between the
!! atmosphere and the ice or that have been accumulated over fast thermodynamic
!! steps but will be applied to the slow (mass-changing) thermodynamics.  Some
!! of these are diagnostics, while others are averages of fluxes taken during
!! the fast ice thermodynamics and used during the slow ice thermodynamics or dynamics.
type fast_ice_avg_type
!FAST ONLY
  integer :: avg_count  ! The number of times that surface fluxes to the ice
                        ! have been incremented.
  logical :: atmos_winds ! The wind stresses come directly from the atmosphere
                         ! model and have the wrong sign.
  ! These are the arrays that are averaged over the fast thermodynamics.  They
  ! are either used to communicate to the slow thermodynamics or diagnostics or
  ! both.
  real, allocatable, dimension(:,:,:) :: &
    ! The 3rd dimension in each of the following is ice thickness category.
    flux_u_top         , & ! The downward flux of zonal and meridional
    flux_v_top         , & ! momentum on an A-grid in Pa.
    flux_t_top         , & ! The upward sensible heat flux at the ice top
                           ! in W m-2.
    flux_q_top         , & ! The upward evaporative moisture flux at
                           ! top of the ice, in kg m-2 s-1.
    flux_lw_top        , & ! The net downward flux of longwave radiation at the
                           ! top of the ice, in W m-2.
    flux_sw_vis_dir_top, & ! The downward diffuse flux of direct (dir)
    flux_sw_vis_dif_top, & ! and diffuse (dif) shortwave radiation in
    flux_sw_nir_dir_top, & ! the visible (vis) and near-infrared (nir)
    flux_sw_nir_dif_top, & ! bands at the top of the ice, in W m-2.
    flux_lh_top        , & ! The upward flux of latent heat at the top
                           ! of the ice, in W m-2.
    lprec_top          , & ! The downward flux of liquid precipitation
                           ! at the top of the ice, in kg m-2 s-1.
    fprec_top          , & ! The downward flux of frozen precipitation
                           ! at the top of the ice, in kg m-2 s-1.
    tmelt              , & ! Ice-top melt energy into the ice/snow in J m-2.
    bmelt              , & ! Ice-bottom melting energy into the ice in J m-2.
    Tskin_cat          , & ! The ice skin temperature by category, in degC.
    sw_abs_ocn      !  The fraction of the absorbed shortwave radiation that is
                    !  absorbed in the ocean, nondim and <=1.
                    !  Equivalent sw_abs_ocn fields are in both the fast_ice_avg_type
                    !  and the ice_rad_type because it is used as a part of the slow
                    !  thermodynamic updates.
  real, allocatable, dimension(:,:) :: &
    bheat      , &   ! The upward diffusive heat flux from the ocean
                     ! to the ice at the base of the ice, in W m-2.
    WindStr_x  , &   ! The zonal wind stress averaged over the ice
                     ! categories on an A-grid, in Pa.
    WindStr_y  , &   ! The meridional wind stress averaged over the
                     ! ice categories on an A-grid, in Pa.
    WindStr_ocn_x, & ! The zonal wind stress on open water on an A-grid, in Pa.
    WindStr_ocn_y, & ! The meridional wind stress on open water on an A-grid, in Pa.
    p_atm_surf , &   ! The atmospheric pressure at the top of the ice, in Pa.
    flux_sw_dn, &    ! The total downward shortwave flux, summed across all
                     ! wavelengths and averaged across all thickness categories
                     ! in W m-2.
    runoff, &        ! Liquid runoff into the ocean, in kg m-2.
    calving, &       ! Calving of ice or runoff of frozen fresh
                     ! water into the ocean, in kg m-2.
    runoff_hflx, &   ! The heat flux associated with runoff, based
                     ! on the temperature difference relative to a
                     ! reference temperature, in ???.
    calving_hflx, &  ! The heat flux associated with calving, based
                     ! on the temperature difference relative to a
                     ! reference temperature, in ???.
    calving_preberg, &  ! Calving of ice or runoff of frozen fresh
                        ! water into the ocean, exclusive of any
                        ! iceberg contributions, in kg m-2.
    calving_hflx_preberg, & ! The heat flux associated with calving,
                        ! exclusive of any iceberg contributions, based on
                        ! the temperature difference relative to a
                        ! reference temperature, in ???.
    Tskin_avg, &     ! The area-weighted average skin temperature across all
                     ! ice thickness categories, in deg C, or 0 if there is
                     ! no ice.
    ice_free   , &   ! The fractional open water used in calculating
                     ! WindStr_[xy]_A; nondimensional, between 0 & 1.
    ice_cover        ! The fractional ice coverage, summed across all
                     ! thickness categories, used in calculating
                     ! WindStr_[xy]_A; nondimensional, between 0 & 1.

  logical :: first_copy = .true.
  integer :: num_tr_fluxes = -1   ! The number of tracer flux fields
  real, allocatable, dimension(:,:,:,:) :: &
    tr_flux_top    ! An array of tracer fluxes at the top of the
                   ! sea ice.

  ! These are the arrays that are averaged over the fast thermodynamics and
  ! then interpolated into unoccupied categories for the purpose of redoing
  ! the application of the fast thermodynamics
  real, allocatable, dimension(:,:,:) :: &
    flux_t0, &  ! The upward sensible heat flux at the ice top
                ! extrapolated to a skin temperature of 0 deg C, in W m-2.
    flux_q0, &  ! The upward evaporative moisture flux at the top of the ice
                ! extrapolated to a skin temperature of 0 deg C, in kg m-2 s-1.
    flux_lw0, & ! The net downward flux of longwave radiation at the top of the
                ! ice extrapolated to a skin temperature of 0 deg C, in W m-2.
    dhdt, &     ! The partial derivative of flux_t0 with ice skin temperature
                ! in W m-2 K-1.
    dedt, &     ! The partial derivative of flux_q0 with ice skin temperature
                ! in kg m-2 s-1 K-1.
    dlwdt       ! The partial derivative of flux_lw0 with ice skin temperature
                ! in W m-2 K-1.

!SLOW ONLY
  real, allocatable, dimension(:,:) :: &
    frazil_left    ! The frazil heat flux that has not yet been
                   ! consumed in making ice, in J m-2. This array
                   ! is decremented by the ice model as the heat
                   ! flux is used up.
!SLOW ONLY
  integer :: id_sh=-1, id_lh=-1, id_sw=-1, id_slp=-1
  integer :: id_lw=-1, id_snofl=-1, id_rain=-1,  id_evap=-1
  integer :: id_sw_vis_dir=-1, id_sw_vis_dif=-1, id_sw_nir_dir=-1, id_sw_nir_dif=-1
  integer :: id_sw_vis=-1, id_sw_dir=-1, id_sw_dif=-1, id_sw_dn=-1, id_albedo=-1
  integer :: id_runoff=-1, id_calving=-1, id_runoff_hflx=-1, id_calving_hflx=-1
  integer :: id_tmelt=-1, id_bmelt=-1, id_bheat=-1
  integer :: id_tsfc=-1, id_sitemptop=-1

  integer :: id_evap_cat=-1, id_lw_cat=-1, id_sh_cat=-1, id_tsfc_cat=-1
  integer :: id_evap0=-1, id_lw0=-1, id_sh0=-1
  integer :: id_dedt=-1, id_dlwdt=-1, id_dshdt=-1
end type fast_ice_avg_type

!> total_sfc_flux_type contains variables that describe the fluxes between the
!! atmosphere and the ice or ocean that have been accumulated over fast thermodynamic
!! steps and integrated across the part-size categories.
type total_sfc_flux_type

  ! These are the arrays that are averaged over the categories and in time over
  ! the fast thermodynamics.
  real, allocatable, dimension(:,:) :: &
    flux_u         , & ! The downward flux of zonal and meridional
    flux_v         , & ! momentum on an A-grid in Pa.
    flux_t         , & ! The upward sensible heat flux at the ice top
                       ! in W m-2.
    flux_q         , & ! The upward evaporative moisture flux at
                       ! top of the ice, in kg m-2 s-1.
    flux_lw        , & ! The downward flux of longwave radiation at
                       ! the top of the ice, in W m-2.
    flux_sw_vis_dir, & ! The downward diffuse flux of direct (dir)
    flux_sw_vis_dif, & ! and diffuse (dif) shortwave radiation in
    flux_sw_nir_dir, & ! the visible (vis) and near-infrared (nir)
    flux_sw_nir_dif, & ! bands at the top of the ice, in W m-2.
    flux_lh        , & ! The upward flux of latent heat at the top
                       ! of the ice, in W m-2.
    lprec          , & ! The downward flux of liquid precipitation
                       ! at the top of the ice, in kg m-2 s-1.
    fprec              ! The downward flux of frozen precipitation
                       ! at the top of the ice, in kg m-2 s-1.
!  logical :: first_copy = .true.
  integer :: num_tr_fluxes = -1   ! The number of tracer flux fields
  real, allocatable, dimension(:,:,:) :: &
    tr_flux        ! An array of tracer fluxes at the top of the
                   ! sea ice.
end type total_sfc_flux_type


!> ice_rad_type contains variables that describe the absorption and reflection
!! of shortwave radiation in and around the sea ice.
type ice_rad_type

  ! The ice skin temperature that can next be used for radiation
  real, allocatable, dimension(:,:,:) :: &
    t_skin      ! The surface skin temperature as calculated by the most
                ! recent fast atmospheric timestep, or a value filled in
                ! from other ice categories or the local freezing point of
                ! seawater when there is no ice at all, in degrees Celsius.
  ! Shortwave absorption parameters that are set in ice_optics.
  real, allocatable, dimension(:,:,:) :: &
    sw_abs_sfc , &  !< The fraction of the absorbed shortwave radiation that is
                    !! absorbed in a surface skin layer, nondim and <=1.
    sw_abs_snow, &  !< The fraction of the absorbed shortwave radiation that is
                    !! absorbed in the snow, nondim and <=1.
    sw_abs_ocn , &  !< The fraction of the absorbed shortwave radiation that is
                    !! absorbed in the ocean, nondim and <=1.
                    !  Only sw_abs_ocn is used in the slow step.
    sw_abs_int      !< The fraction of the absorbed shortwave radiation that is
                    !! absorbed by all ice layers in aggregate, nondim and <=1.
                    !  sw_abs_int is only used for diagnostics.
  real, allocatable, dimension(:,:,:,:) :: &
    sw_abs_ice      !< The fraction of the absorbed shortwave that is
                    !! absorbed in each of the ice layers, nondim, <=1.

  real, allocatable, dimension(:,:)   :: &
    coszen_nextrad  !< Cosine of the solar zenith angle averaged
                    !! over the next radiation timestep, nondim.

  logical :: add_diurnal_sw       !< If true, apply a synthetic diurnal cycle to
                                  !! the shortwave radiation.
  logical :: do_sun_angle_for_alb !< If true, find the sun angle for calculating
                                  !! the ocean albedo in the frame of the ice model.
  logical :: frequent_albedo_update !< If true, update the ice and ocean albedos
                                  !! within the fast ice model update.  Otherwise,
                                  !! the albedos are only updated within
                                  !! set_ice_surface_state.

  integer, allocatable, dimension(:)   :: id_sw_abs_ice
  integer :: id_sw_abs_sfc=-1, id_sw_abs_snow=-1, id_sw_pen=-1, id_sw_abs_ocn=-1
  integer :: id_alb=-1, id_coszen=-1, id_swdn=-1, id_lwdn=-1
  integer :: id_alb_vis_dir=-1, id_alb_vis_dif=-1, id_alb_nir_dir=-1, id_alb_nir_dif=-1
  integer :: id_tskin=-1, id_cn=-1, id_mi=-1

end type ice_rad_type

!> ice_ocean_flux_type contains variables that describe the fluxes between the
!! ice and the ocean, on the ice grid.
type ice_ocean_flux_type
  ! These variables describe the fluxes between ice or atmosphere and the ocean.
  real, allocatable, dimension(:,:)   :: &
    flux_t_ocn_top , &  ! The upward sensible heat flux from the ocean
                        ! to the ice or atmosphere, in W m-2.
    flux_q_ocn_top , &  ! The upward evaporative moisture flux at
                        ! the ocean surface, in kg m-2 s-1.
    flux_lw_ocn_top, &  ! The downward flux of longwave radiation at
                        ! the ocean surface, in W m-2.
    flux_sw_vis_dir_ocn, & ! The downward diffuse flux of direct (dir)
    flux_sw_vis_dif_ocn, & ! and diffuse (dif) shortwave radiation in
    flux_sw_nir_dir_ocn, & ! the visible (vis) and near-infrared (nir)
    flux_sw_nir_dif_ocn, & ! bands at the ocean surface, in W m-2.
    flux_lh_ocn_top, &  ! The upward flux of latent heat at the
                        ! ocean surface, in W m-2.
    lprec_ocn_top, &    ! The downward flux of liquid precipitation at
                        ! the ocean surface, in kg m-2 s-1.
    fprec_ocn_top, &    ! The downward flux of frozen precipitation at
                        ! the ocean surface, in kg m-2 s-1.
    flux_u_ocn, &       ! The flux of x-momentum into the ocean, in Pa,
                        ! at locations determined by flux_uv_stagger,
                        ! but allocated as though on an A-grid.
    flux_v_ocn, &       ! The flux of y-momentum into the ocean, in Pa,
                        ! at locations determined by flux_uv_stagger,
                        ! but allocated as though on an A-grid.
    melt_nudge, &       ! A downward fresh water flux into the ocean that
                        ! acts to nudge the ocean surface salinity to
                        ! facilitate the retention of sea ice, in kg m-2 s-1.
    flux_salt           ! The flux of salt out of the ocean in kg m-2.

  !Iceberg fields
  real, pointer, dimension(:,:)   :: &
    ustar_berg =>NULL(), &  ! ustar contribution below icebergs in m/s
    area_berg =>NULL(),  &  ! fraction of grid cell covered by icebergs in m2/m2
    mass_berg =>NULL()      ! mass of icebergs in km/m^2

  ! These arrays are used for enthalpy change diagnostics in the slow thermodynamics.
  real, allocatable, dimension(:,:)   :: &
    ! These terms diagnose the enthalpy change associated with the addition or
    ! removal of water mass (liquid or frozen) from the ice model are required
    ! to close the enthalpy budget. Ice enthalpy is generally negative, so terms
    ! that add mass to the ice are generally negative.
    Enth_Mass_in_atm , & ! The enthalpy introduced to the ice by water
                         ! fluxes from the atmosphere, in J m-2.
    Enth_Mass_out_atm, & ! Negative of the enthalpy extracted from the
                         ! ice by water fluxes to the atmosphere, in J m-2.
    Enth_Mass_in_ocn , & ! The enthalpy introduced to the ice by water
                         ! fluxes from the ocean, in J m-2.
    Enth_Mass_out_ocn    ! Negative of the enthalpy extracted from the
                         ! ice by water fluxes to the ocean, in J m-2.

  integer :: stress_count ! The number of times that the stresses from the ice
                        ! to the ocean have been incremented.
  integer :: flux_uv_stagger = -999 ! The staggering relative to the tracer points
                    ! points of the two wind stress components. Valid entries
                    ! include AGRID, BGRID_NE, CGRID_NE, BGRID_SW, and CGRID_SW,
                    ! corresponding to the community-standard Arakawa notation.
                    ! (These are named integers taken from mpp_parameter_mod.)
                    ! Following SIS, this is BGRID_NE by default when the sea
                    ! ice is initialized, but here it is set to -999 so that a
                    ! global max across ice and non-ice processors can be used
                    ! to determine its value.
  logical :: slp2ocean  ! If true, apply sea level pressure to ocean surface.

!   type(coupler_2d_bc_type)   :: ocean_fluxes       ! array of fluxes used for additional tracers

  integer :: num_tr_fluxes = -1 ! The number of tracer flux fields
  real, allocatable, dimension(:,:,:) :: &
    tr_flux_ocn_top     ! An array of tracer fluxes at the ocean's surface.

  ! diagnostic IDs for ice-to-ocean fluxes.
  integer :: id_saltf=-1
  ! The following are diagnostic IDs for iceberg-related fields.  These are only
  ! used if the iceberg code is activated.
  integer ::  id_ustar_berg=-1, id_area_berg=-1, id_mass_berg=-1
end type ice_ocean_flux_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_IST_arrays allocates the arrays in an ice_state_type.
subroutine alloc_IST_arrays(HI, IG, IST, omit_velocities, omit_Tsurf)
  type(hor_index_type),    intent(in)    :: HI
  type(ice_grid_type),     intent(in)    :: IG
  type(ice_state_type),    intent(inout) :: IST
  logical, optional,       intent(in)    :: omit_velocities
  logical, optional,       intent(in)    :: omit_Tsurf

  integer :: isd, ied, jsd, jed, CatIce, NkIce, idr
  logical :: do_vel, do_Tsurf

  do_vel = .true. ; if (present(omit_velocities)) do_vel = .not.omit_velocities
  do_Tsurf = .true. ; if (present(omit_Tsurf)) do_Tsurf = .not.omit_Tsurf

  CatIce = IG%CatIce ; NkIce = IG%NkIce
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  allocate(IST%part_size(isd:ied, jsd:jed, 0:CatIce)) ; IST%part_size(:,:,:) = 0.0
  allocate(IST%mH_pond(  isd:ied, jsd:jed, CatIce)) ; IST%mH_pond(:,:,:) = 0.0
  allocate(IST%mH_snow(  isd:ied, jsd:jed, CatIce)) ; IST%mH_snow(:,:,:) = 0.0
  allocate(IST%enth_snow(isd:ied, jsd:jed, CatIce, 1)) ; IST%enth_snow(:,:,:,:) = 0.0
  allocate(IST%mH_ice(   isd:ied, jsd:jed, CatIce)) ; IST%mH_ice(:,:,:) = 0.0
  allocate(IST%enth_ice( isd:ied, jsd:jed, CatIce, NkIce)) ; IST%enth_ice(:,:,:,:) = 0.0
  allocate(IST%sal_ice(  isd:ied, jsd:jed, CatIce, NkIce)) ; IST%sal_ice(:,:,:,:) = 0.0

  if (do_vel) then
    ! These velocities are only required for the slow ice processes, and hence
    ! can use the memory macros.
    if (IST%Cgrid_dyn) then
      allocate(IST%u_ice_C(SZIB_(HI), SZJ_(HI))) ; IST%u_ice_C(:,:) = 0.0
      allocate(IST%v_ice_C(SZI_(HI), SZJB_(HI))) ; IST%v_ice_C(:,:) = 0.0
    else
      allocate(IST%u_ice_B(SZIB_(HI), SZJB_(HI))) ; IST%u_ice_B(:,:) = 0.0
      allocate(IST%v_ice_B(SZIB_(HI), SZJB_(HI))) ; IST%v_ice_B(:,:) = 0.0
    endif

    ! ### THESE ARE DIAGNOSTICS.  PERHAPS THEY SHOULD ONLY BE ALLOCATED IF USED.
    allocate(IST%rdg_mice(isd:ied, jsd:jed, CatIce)) ; IST%rdg_mice(:,:,:) = 0.0
  endif

  if (do_Tsurf) then
    ! IST%tsurf is only used with some older options.
    allocate(IST%t_surf(isd:ied, jsd:jed, CatIce)) ; IST%t_surf(:,:,:) = 0.0
  endif

end subroutine alloc_IST_arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_state_register_restarts registers any variables in the ice state type
!!     that need to be includedin the restart files.
subroutine ice_state_register_restarts(mpp_domain, IST, IG, Ice_restart, restart_file)
  type(domain2d),          intent(in)    :: mpp_domain
  type(ice_state_type),    intent(inout) :: IST
  type(ice_grid_type),     intent(in)    :: IG
  type(restart_file_type), pointer       :: Ice_restart
  character(len=*),        intent(in)    :: restart_file

  integer :: idr

  ! Now register some of these arrays to be read from the restart files.
  if (associated(Ice_restart)) then
    idr = register_restart_field(Ice_restart, restart_file, 'part_size', IST%part_size, domain=mpp_domain)
    if (allocated(IST%t_surf)) then
      idr = register_restart_field(Ice_restart, restart_file, 't_surf_ice', IST%t_surf, &
                                 domain=mpp_domain, mandatory=.false., units="deg K")
    endif
    idr = register_restart_field(Ice_restart, restart_file, 'h_pond', IST%mH_pond, & ! mw/new
                                 domain=mpp_domain, mandatory=.false., units="H_to_kg_m2 kg m-2")
    idr = register_restart_field(Ice_restart, restart_file, 'h_snow', IST%mH_snow, &
                                 domain=mpp_domain, mandatory=.true., units="H_to_kg_m2 kg m-2")
    idr = register_restart_field(Ice_restart, restart_file, 'enth_snow', IST%enth_snow, &
                                 domain=mpp_domain, mandatory=.false.)
    idr = register_restart_field(Ice_restart, restart_file, 'h_ice',  IST%mH_ice, &
                                 domain=mpp_domain, mandatory=.true., units="H_to_kg_m2 kg m-2")
    idr = register_restart_field(Ice_restart, restart_file, 'H_to_kg_m2', IG%H_to_kg_m2, &
                                 longname="The conversion factor from SIS2 mass-thickness units to kg m-2.", &
                                 no_domain=.true., mandatory=.false.)

    idr = register_restart_field(Ice_restart, restart_file, 'enth_ice', IST%enth_ice, &
                                 domain=mpp_domain, mandatory=.false., units="J kg-1")
    idr = register_restart_field(Ice_restart, restart_file, 'sal_ice', IST%sal_ice, &
                                 domain=mpp_domain, mandatory=.false., units="kg/kg")

    if (IST%Cgrid_dyn) then
      idr = register_restart_field(Ice_restart, restart_file, 'u_ice_C', IST%u_ice_C, &
                                   domain=mpp_domain, position=EAST, mandatory=.false.)
      idr = register_restart_field(Ice_restart, restart_file, 'v_ice_C', IST%v_ice_C, &
                                   domain=mpp_domain, position=NORTH, mandatory=.false.)
    else
      idr = register_restart_field(Ice_restart, restart_file, 'u_ice',   IST%u_ice_B, &
                                   domain=mpp_domain, position=CORNER, mandatory=.false.)
      idr = register_restart_field(Ice_restart, restart_file, 'v_ice',   IST%v_ice_B, &
                                   domain=mpp_domain, position=CORNER, mandatory=.false.)
    endif
  endif

end subroutine ice_state_register_restarts


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_fast_ice_avg allocates and zeros out the arrays in a fast_ice_avg_type.
subroutine alloc_fast_ice_avg(FIA, HI, IG, interp_fluxes)
  type(fast_ice_avg_type), pointer    :: FIA
  type(hor_index_type),    intent(in) :: HI
  type(ice_grid_type),     intent(in) :: IG
  logical,                 intent(in) :: interp_fluxes

  integer :: isd, ied, jsd, jed, CatIce

  if (.not.associated(FIA)) allocate(FIA)
  CatIce = IG%CatIce
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  FIA%avg_count = 0
  allocate(FIA%flux_u_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_u_top(:,:,:) = 0.0
  allocate(FIA%flux_v_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_v_top(:,:,:) = 0.0
  allocate(FIA%flux_t_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_t_top(:,:,:) = 0.0
  allocate(FIA%flux_q_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_q_top(:,:,:) = 0.0
  allocate(FIA%flux_sw_vis_dir_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_sw_vis_dir_top(:,:,:) = 0.0
  allocate(FIA%flux_sw_vis_dif_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_sw_vis_dif_top(:,:,:) = 0.0
  allocate(FIA%flux_sw_nir_dir_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_sw_nir_dir_top(:,:,:) = 0.0
  allocate(FIA%flux_sw_nir_dif_top(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_sw_nir_dif_top(:,:,:) = 0.0
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
  allocate(FIA%bheat(isd:ied, jsd:jed)) ; FIA%bheat(:,:) = 0.0
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
    allocate(FIA%flux_t0(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_t0(:,:,:) = 0.0
    allocate(FIA%flux_q0(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_q0(:,:,:) = 0.0
    allocate(FIA%flux_lw0(isd:ied, jsd:jed, 0:CatIce)) ; FIA%flux_lw0(:,:,:) = 0.0
    allocate(FIA%dhdt(isd:ied, jsd:jed, 0:CatIce)) ; FIA%dhdt(:,:,:) = 0.0
    allocate(FIA%dedt(isd:ied, jsd:jed, 0:CatIce)) ; FIA%dedt(:,:,:) = 0.0
    allocate(FIA%dlwdt(isd:ied, jsd:jed, 0:CatIce)) ; FIA%dlwdt(:,:,:) = 0.0
    allocate(FIA%Tskin_cat(isd:ied, jsd:jed, 0:CatIce)) ; FIA%Tskin_cat(:,:,:) = 0.0
  endif

  allocate(FIA%flux_sw_dn(isd:ied, jsd:jed))  ; FIA%flux_sw_dn(:,:) = 0.0
  allocate(FIA%sw_abs_ocn(isd:ied, jsd:jed, CatIce)) ; FIA%sw_abs_ocn(:,:,:) = 0.0

end subroutine alloc_fast_ice_avg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_total_sfc_flux allocates and zeros out the arrays in a total_sfc_flux_type.
subroutine alloc_total_sfc_flux(TSF, HI)
  type(total_sfc_flux_type), pointer    :: TSF
  type(hor_index_type),      intent(in) :: HI

  integer :: isd, ied, jsd, jed

  if (.not.associated(TSF)) allocate(TSF)
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  allocate(TSF%flux_u(isd:ied, jsd:jed)) ; TSF%flux_u(:,:) = 0.0
  allocate(TSF%flux_v(isd:ied, jsd:jed)) ; TSF%flux_v(:,:) = 0.0
  allocate(TSF%flux_t(isd:ied, jsd:jed)) ; TSF%flux_t(:,:) = 0.0
  allocate(TSF%flux_q(isd:ied, jsd:jed)) ; TSF%flux_q(:,:) = 0.0
  allocate(TSF%flux_sw_vis_dir(isd:ied, jsd:jed)) ; TSF%flux_sw_vis_dir(:,:) = 0.0
  allocate(TSF%flux_sw_vis_dif(isd:ied, jsd:jed)) ; TSF%flux_sw_vis_dif(:,:) = 0.0
  allocate(TSF%flux_sw_nir_dir(isd:ied, jsd:jed)) ; TSF%flux_sw_nir_dir(:,:) = 0.0
  allocate(TSF%flux_sw_nir_dif(isd:ied, jsd:jed)) ; TSF%flux_sw_nir_dif(:,:) = 0.0
  allocate(TSF%flux_lw(isd:ied, jsd:jed)) ; TSF%flux_lw(:,:) = 0.0
  allocate(TSF%flux_lh(isd:ied, jsd:jed)) ; TSF%flux_lh(:,:) = 0.0
  allocate(TSF%lprec(isd:ied, jsd:jed)) ;  TSF%lprec(:,:) = 0.0
  allocate(TSF%fprec(isd:ied, jsd:jed)) ;  TSF%fprec(:,:) = 0.0

end subroutine alloc_total_sfc_flux


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_rad_register_restarts allocates the arrays in the ice_rad_type
!!     and registers any variables in the ice rad type that need to be included
!!     in the restart files.
subroutine ice_rad_register_restarts(mpp_domain, HI, IG, param_file, Rad, &
                                       Ice_restart, restart_file)
  type(domain2d),          intent(in)    :: mpp_domain
  type(hor_index_type),    intent(in)    :: HI
  type(ice_grid_type),     intent(in)    :: IG
  type(param_file_type),   intent(in)    :: param_file
  type(ice_rad_type),      pointer       :: Rad
  type(restart_file_type), intent(inout) :: Ice_restart
  character(len=*),        intent(in)    :: restart_file

  integer :: isd, ied, jsd, jed, CatIce, NkIce, idr

  if (.not.associated(Rad)) allocate(Rad)
  CatIce = IG%CatIce ; NkIce = IG%NkIce
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  allocate(Rad%t_skin(isd:ied, jsd:jed, CatIce)) ; Rad%t_skin(:,:,:) = 0.0

  allocate(Rad%sw_abs_sfc(isd:ied, jsd:jed, CatIce)) ; Rad%sw_abs_sfc(:,:,:) = 0.0
  allocate(Rad%sw_abs_snow(isd:ied, jsd:jed, CatIce)) ; Rad%sw_abs_snow(:,:,:) = 0.0
  allocate(Rad%sw_abs_ice(isd:ied, jsd:jed, CatIce, NkIce)) ; Rad%sw_abs_ice(:,:,:,:) = 0.0
  allocate(Rad%sw_abs_ocn(isd:ied, jsd:jed, CatIce)) ; Rad%sw_abs_ocn(:,:,:) = 0.0
  allocate(Rad%sw_abs_int(isd:ied, jsd:jed, CatIce)) ; Rad%sw_abs_int(:,:,:) = 0.0

  allocate(Rad%coszen_nextrad(isd:ied, jsd:jed)) ; Rad%coszen_nextrad(:,:) = 0.0

  idr = register_restart_field(Ice_restart, restart_file, 'coszen', Rad%coszen_nextrad, &
                               domain=mpp_domain, mandatory=.false.)
  idr = register_restart_field(Ice_restart, restart_file, 'T_skin', Rad%t_skin, &
                               domain=mpp_domain, mandatory=.false.)

end subroutine ice_rad_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_ice_ocean_flux allocates and zeros out the arrays in an ice_ocean_flux_type.
subroutine alloc_ice_ocean_flux(IOF, HI, do_iceberg_fields)
  type(ice_ocean_flux_type), pointer    :: IOF
  type(hor_index_type),      intent(in) :: HI
  logical,         optional, intent(in) :: do_iceberg_fields

  integer :: CatIce
  logical :: alloc_bergs

  alloc_bergs = .false. ; if (present(do_iceberg_fields)) alloc_bergs = do_iceberg_fields

  if (.not.associated(IOF)) allocate(IOF)

  allocate(IOF%flux_salt(SZI_(HI), SZJ_(HI))) ; IOF%flux_salt(:,:) = 0.0

  allocate(IOF%flux_t_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%flux_t_ocn_top(:,:) = 0.0
  allocate(IOF%flux_q_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%flux_q_ocn_top(:,:) = 0.0
  allocate(IOF%flux_lw_ocn_top(SZI_(HI), SZJ_(HI))) ; IOF%flux_lw_ocn_top(:,:) = 0.0
  allocate(IOF%flux_lh_ocn_top(SZI_(HI), SZJ_(HI))) ; IOF%flux_lh_ocn_top(:,:) = 0.0
  allocate(IOF%flux_sw_vis_dir_ocn(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sw_vis_dir_ocn(:,:) = 0.0
  allocate(IOF%flux_sw_vis_dif_ocn(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sw_vis_dif_ocn(:,:) = 0.0
  allocate(IOF%flux_sw_nir_dir_ocn(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sw_nir_dir_ocn(:,:) = 0.0
  allocate(IOF%flux_sw_nir_dif_ocn(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sw_nir_dif_ocn(:,:) = 0.0
  allocate(IOF%lprec_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%lprec_ocn_top(:,:) = 0.0
  allocate(IOF%fprec_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%fprec_ocn_top(:,:) = 0.0
  allocate(IOF%flux_u_ocn(SZI_(HI), SZJ_(HI)))    ;  IOF%flux_u_ocn(:,:) = 0.0
  allocate(IOF%flux_v_ocn(SZI_(HI), SZJ_(HI)))    ;  IOF%flux_v_ocn(:,:) = 0.0

  allocate(IOF%Enth_Mass_in_atm(SZI_(HI), SZJ_(HI)))  ; IOF%Enth_Mass_in_atm(:,:) = 0.0
  allocate(IOF%Enth_Mass_out_atm(SZI_(HI), SZJ_(HI))) ; IOF%Enth_Mass_out_atm(:,:) = 0.0
  allocate(IOF%Enth_Mass_in_ocn(SZI_(HI), SZJ_(HI)))  ; IOF%Enth_Mass_in_ocn(:,:) = 0.0
  allocate(IOF%Enth_Mass_out_ocn(SZI_(HI), SZJ_(HI))) ; IOF%Enth_Mass_out_ocn(:,:) = 0.0

  !Allocating iceberg fields (only used if pass_iceberg_area_to_ocean=.True.)
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
subroutine alloc_ocean_sfc_state(OSS, HI, Cgrid_dyn)
  type(ocean_sfc_state_type), pointer    :: OSS
  type(hor_index_type),       intent(in) :: HI
  logical,                    intent(in) :: Cgrid_dyn

  if (.not.associated(OSS)) allocate(OSS)

  ! The ocean_sfc_state_type only occurs on slow ice PEs, so it can use the memory macros.
  allocate(OSS%s_surf(SZI_(HI), SZJ_(HI))) ; OSS%s_surf(:,:) = 0.0
  allocate(OSS%SST_C(SZI_(HI), SZJ_(HI)))  ; OSS%SST_C(:,:) = 0.0
  allocate(OSS%T_fr_ocn(SZI_(HI), SZJ_(HI))) ; OSS%T_fr_ocn(:,:) = 0.0
  allocate(OSS%sea_lev(SZI_(HI), SZJ_(HI))) ; OSS%sea_lev(:,:) = 0.0
  allocate(OSS%frazil(SZI_(HI), SZJ_(HI))) ; OSS%frazil(:,:) = 0.0

  if (Cgrid_dyn) then
    allocate(OSS%u_ocn_C(SZIB_(HI), SZJ_(HI))) ; OSS%u_ocn_C(:,:) = 0.0
    allocate(OSS%v_ocn_C(SZI_(HI), SZJB_(HI))) ; OSS%v_ocn_C(:,:) = 0.0
  else
    allocate(OSS%u_ocn_B(SZIB_(HI), SZJB_(HI))) ; OSS%u_ocn_B(:,:) = 0.0
    allocate(OSS%v_ocn_B(SZIB_(HI), SZJB_(HI))) ; OSS%v_ocn_B(:,:) = 0.0
  endif

  OSS%Cgrid_dyn = Cgrid_dyn

end subroutine alloc_ocean_sfc_state


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_simple_ocean_sfc_state allocates and zeros out the arrays in a
!! simple_OSS_type.
subroutine alloc_simple_OSS(OSS, HI)
  type(simple_OSS_type), pointer    :: OSS
  type(hor_index_type),  intent(in) :: HI

  integer :: isd, ied, jsd, jed

  if (.not.associated(OSS)) allocate(OSS)
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  allocate(OSS%s_surf(isd:ied, jsd:jed)) ; OSS%s_surf(:,:) = 0.0
  allocate(OSS%SST_C(isd:ied, jsd:jed))  ; OSS%SST_C(:,:) = 0.0
  allocate(OSS%T_fr_ocn(isd:ied, jsd:jed)) ; OSS%T_fr_ocn(:,:) = 0.0
  allocate(OSS%bheat(isd:ied, jsd:jed))   ; OSS%bheat(:,:) = 0.0
  allocate(OSS%u_ocn_A(isd:ied, jsd:jed)) ; OSS%u_ocn_A(:,:) = 0.0
  allocate(OSS%v_ocn_A(isd:ied, jsd:jed)) ; OSS%v_ocn_A(:,:) = 0.0
  allocate(OSS%u_ice_A(isd:ied, jsd:jed)) ; OSS%u_ice_A(:,:) = 0.0
  allocate(OSS%v_ice_A(isd:ied, jsd:jed)) ; OSS%v_ice_A(:,:) = 0.0

end subroutine alloc_simple_OSS


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> copy_IST_to_IST copies the computational domain of one ice state type into
!! the computational domain of another ice_state_type.  Both must use the same
!! domain decomposition and indexing convention (for now), but they may have
!! different halo sizes.
subroutine copy_IST_to_IST(IST_in, IST_out, HI_in, HI_out, IG)
  type(ice_state_type), intent(in)    :: IST_in
  type(ice_state_type), intent(inout) :: IST_out
  type(hor_index_type), intent(in)    :: HI_in, HI_out
  type(ice_grid_type),  intent(in)    :: IG

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

  ! The velocity components, rdg_mice, TrReg, and ITV are deliberately not being
  ! copied.

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

  ! The velocity components, rdg_mice, TrReg, and ITV are deliberately not being
  ! copied.
  if (associated(IST_out) .and. associated(IST_in)) then
    call mpp_redistribute(domain_in, IST_in%part_size, domain_out, &
                          IST_out%part_size, complete=.true.)

    if (allocated(IST_out%t_surf) .or. allocated(IST_in%t_surf)) then
      call mpp_redistribute(domain_in, IST_in%t_surf, domain_out, &
                          IST_out%t_surf, complete=.false.)
    endif
    call mpp_redistribute(domain_in, IST_in%mH_pond, domain_out, &
                          IST_out%mH_pond, complete=.false.)
    call mpp_redistribute(domain_in, IST_in%mH_snow, domain_out, &
                          IST_out%mH_snow, complete=.false.)
    call mpp_redistribute(domain_in, IST_in%mH_ice, domain_out, &
                          IST_out%mH_ice, complete=.false.)
    call mpp_redistribute(domain_in, IST_in%enth_snow, domain_out, &
                          IST_out%enth_snow, complete=.true.)

    call mpp_redistribute(domain_in, IST_in%enth_ice, domain_out, &
                          IST_out%enth_ice, complete=.false.)
    call mpp_redistribute(domain_in, IST_in%sal_ice, domain_out, &
                          IST_out%sal_ice, complete=.true.)
  elseif (associated(IST_out)) then
    ! Use the null pointers in place of the unneeded input arrays.
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          IST_out%part_size, complete=.true.)

    if (allocated(IST_out%t_surf)) then
      call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          IST_out%t_surf, complete=.false.)
    endif
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          IST_out%mH_pond, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          IST_out%mH_snow, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          IST_out%mH_ice, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr4D, domain_out, &
                          IST_out%enth_snow, complete=.true.)

    call mpp_redistribute(domain_in, null_ptr4D, domain_out, &
                          IST_out%enth_ice, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr4D, domain_out, &
                          IST_out%sal_ice, complete=.true.)
  elseif (associated(IST_in)) then
    ! Use the null pointers in place of the unneeded output arrays.
    call mpp_redistribute(domain_in, IST_in%part_size, domain_out, &
                          null_ptr3D, complete=.true.)

    if (allocated(IST_in%t_surf)) then
      call mpp_redistribute(domain_in, IST_in%t_surf, domain_out, &
                          null_ptr3D, complete=.false.)
    endif
    call mpp_redistribute(domain_in, IST_in%mH_pond, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, IST_in%mH_snow, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, IST_in%mH_ice, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, IST_in%enth_snow, domain_out, &
                          null_ptr4D, complete=.true.)

    call mpp_redistribute(domain_in, IST_in%enth_ice, domain_out, &
                          null_ptr4D, complete=.false.)
    call mpp_redistribute(domain_in, IST_in%sal_ice, domain_out, &
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
subroutine translate_OSS_to_sOSS(OSS, IST, sOSS, G)
  type(ocean_sfc_state_type), intent(in)    :: OSS
  type(ice_state_type),       intent(in)    :: IST
  type(simple_OSS_type),      intent(inout) :: sOSS
  type(SIS_hor_grid_type),    intent(in)    :: G

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  !$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,sOSS,OSS,IST)
  do j=jsc,jec ; do i=isc,iec
    sOSS%s_surf(i,j) = OSS%s_surf(i,j)
    sOSS%SST_C(i,j) = OSS%SST_C(i,j)
    sOSS%T_fr_ocn(i,j) = OSS%T_fr_ocn(i,j)

    if (G%mask2dT(i,j) > 0.5) then
      sOSS%bheat(i,j) = OSS%kmelt*(OSS%SST_C(i,j) - sOSS%T_fr_ocn(i,j))
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

  if (sOSS%num_tr<0) then
    sOSS%num_tr = OSS%num_tr
    if (sOSS%num_tr > 0) then
      allocate(sOSS%tr_array(G%isd:G%ied,G%jsd:G%jed,sOSS%num_tr)) ; sOSS%tr_array(:,:,:) = 0.0
    endif
  endif
  do m=1,OSS%num_tr ; do j=jsc,jec ; do i=isc,iec
    sOSS%tr_array(i,j,m) = OSS%tr_array(i,j,m)
  enddo ; enddo ; enddo

end subroutine translate_OSS_to_sOSS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> copy_sOSS_to_sOSS copies the computational domain of one simple_OSS_type into
!! the computational domain of another simple_OSS_type.  Both must use the same
!! domain decomposition and indexing convention (for now), but they may have
!! different halo sizes.
subroutine copy_sOSS_to_sOSS(OSS_in, OSS_out, HI_in, HI_out)
  type(simple_OSS_type), intent(inout) :: OSS_in
  type(simple_OSS_type), intent(inout) :: OSS_out
  type(hor_index_type),  intent(in)    :: HI_in, HI_out

  integer :: i, j, m, isc, iec, jsc, jec
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

  if (OSS_out%first_copy) then
    OSS_in%first_copy = .false. ; OSS_out%first_copy = .false.
    OSS_out%num_tr = OSS_in%num_tr
    if (OSS_out%num_tr > 0) then
      allocate(OSS_out%tr_array(HI_out%isd:HI_out%ied,HI_out%jsd:HI_out%jed,OSS_out%num_tr))
      OSS_out%tr_array(:,:,:) = 0.0
    endif
  endif

  do m=1,OSS_in%num_tr ; do j=jsc,jec ; do i=isc,iec ; i2=i+i_off ; j2=j+j_off
    OSS_out%tr_array(i2,j2,m) = OSS_in%tr_array(i,j,m)
  enddo ; enddo ; enddo

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
  logical :: first_copy
  integer :: m, num_tr

  if (.not. (associated(OSS_out) .or. associated(OSS_in))) &
    call SIS_error(FATAL, "redistribute_sOSS_to_sOSS called with "//&
                          "neither OSS_in nor OSS_out associated.")
  first_copy = .false.
  if (associated(OSS_out)) first_copy = OSS_out%first_copy
  if (associated(OSS_in)) first_copy = first_copy .or. OSS_in%first_copy

  if (first_copy) then
    ! Determine the number of fluxes.
    num_tr = 0 ; if (associated(OSS_in)) num_tr = OSS_in%num_tr
    call max_across_PEs(num_tr)

    if (associated(OSS_out)) then
      if (.not. present(HI_out)) &
        call SIS_error(FATAL, "redistribute_sOSS_to_sOSS called with an "//&
                              "associated OSS_out but without HI_out.")
      OSS_out%num_tr = num_tr
      if ((num_tr > 0) .and. .not.allocated(OSS_out%tr_array)) then
        allocate(OSS_out%tr_array(HI_out%isd:HI_out%ied,HI_out%jsd:HI_out%jed,num_tr))
        OSS_out%tr_array(:,:,:) = 0.0
      endif
      OSS_out%first_copy = .false.
    endif

    if (associated(OSS_in)) OSS_in%first_copy = .false.
  endif

  if (associated(OSS_out) .and. associated(OSS_in)) then
    ! The extra tracer arrays are copied first so that they can all have
    ! complete=.false.
    do m=1,OSS_in%num_tr
      call mpp_redistribute(domain_in, OSS_in%tr_array(:,:,m), domain_out, &
                            OSS_out%tr_array(:,:,m), complete=.false.)
    enddo

    call mpp_redistribute(domain_in, OSS_in%SST_C, domain_out, &
                          OSS_out%SST_C, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%s_surf, domain_out, &
                          OSS_out%s_surf, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%T_fr_ocn, domain_out, &
                          OSS_out%T_fr_ocn, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%bheat, domain_out, &
                          OSS_out%bheat, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%u_ocn_A, domain_out, &
                          OSS_out%u_ocn_A, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%v_ocn_A, domain_out, &
                          OSS_out%v_ocn_A, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%u_ice_A, domain_out, &
                          OSS_out%u_ice_A, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%v_ice_A, domain_out, &
                          OSS_out%v_ice_A, complete=.true.)
  elseif (associated(OSS_out)) then
    ! Use the null pointer in place of the unneeded input arrays.
    do m=1,OSS_out%num_tr
      call mpp_redistribute(domain_in, null_ptr, domain_out, &
                            OSS_out%tr_array(:,:,m), complete=.false.)
    enddo

    call mpp_redistribute(domain_in, null_ptr, domain_out, &
                          OSS_out%SST_C, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr, domain_out, &
                          OSS_out%s_surf, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr, domain_out, &
                          OSS_out%T_fr_ocn, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr, domain_out, &
                          OSS_out%bheat, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr, domain_out, &
                          OSS_out%u_ocn_A, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr, domain_out, &
                          OSS_out%v_ocn_A, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr, domain_out, &
                          OSS_out%u_ice_A, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr, domain_out, &
                          OSS_out%v_ice_A, complete=.true.)
  elseif (associated(OSS_in)) then
    ! Use the null pointer in place of the unneeded output arrays.
    do m=1,OSS_in%num_tr
      call mpp_redistribute(domain_in, OSS_in%tr_array(:,:,m), domain_out, &
                            null_ptr, complete=.false.)
    enddo

    call mpp_redistribute(domain_in, OSS_in%SST_C, domain_out, &
                          null_ptr, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%s_surf, domain_out, &
                          null_ptr, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%T_fr_ocn, domain_out, &
                          null_ptr, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%bheat, domain_out, &
                          null_ptr, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%u_ocn_A, domain_out, &
                          null_ptr, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%v_ocn_A, domain_out, &
                          null_ptr, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%u_ice_A, domain_out, &
                          null_ptr, complete=.false.)
    call mpp_redistribute(domain_in, OSS_in%v_ice_A, domain_out, &
                          null_ptr, complete=.true.)
  else
    call SIS_error(FATAL, "redistribute_sOSS_to_sOSS called with "//&
                          "neither OSS_in nor OSS_out associated.")
  endif

end subroutine redistribute_sOSS_to_sOSS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> copy_FIA_to_FIA copies the computational domain of one fast_ice_avg_type into
!! the computational domain of another fast_ice_avg_type.  Both must use the same
!! domain decomposition and indexing convention (for now), but they may have
!! different halo sizes.
subroutine copy_FIA_to_FIA(FIA_in, FIA_out, HI_in, HI_out, IG)
  type(fast_ice_avg_type), intent(inout) :: FIA_in
  type(fast_ice_avg_type), intent(inout) :: FIA_out
  type(hor_index_type),    intent(in)    :: HI_in, HI_out
  type(ice_grid_type),     intent(in)    :: IG

  integer :: i, j, k, m, n, isc, iec, jsc, jec, ncat, NkIce ! , i_off, j_off
  integer :: i2, j2, i_off, j_off
  integer :: isd, ied, jsd, jed

  isc = HI_in%isc ; iec = HI_in%iec ; jsc = HI_in%jsc ; jec = HI_in%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  if ((HI_in%iec-HI_in%isc /= HI_out%iec-HI_out%isc) .or. &
      (HI_in%jec-HI_in%jsc /= HI_out%jec-HI_out%jsc)) then
    call SIS_error(FATAL, "copy_FIA_to_FIA called with inconsistent domain "//&
                          "decompositions of the two ice types.")
  endif
  i_off = HI_out%iec-HI_in%iec ;  j_off = HI_out%jec-HI_in%jec

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    FIA_out%flux_t_top(i2,j2,k) = FIA_in%flux_t_top(i,j,k)
    FIA_out%flux_q_top(i2,j2,k) = FIA_in%flux_q_top(i,j,k)
    FIA_out%flux_sw_vis_dir_top(i2,j2,k) = FIA_in%flux_sw_vis_dir_top(i,j,k)
    FIA_out%flux_sw_vis_dif_top(i2,j2,k) = FIA_in%flux_sw_vis_dif_top(i,j,k)
    FIA_out%flux_sw_nir_dir_top(i2,j2,k) = FIA_in%flux_sw_nir_dir_top(i,j,k)
    FIA_out%flux_sw_nir_dif_top(i2,j2,k) = FIA_in%flux_sw_nir_dif_top(i,j,k)
    FIA_out%flux_lw_top(i2,j2,k) = FIA_in%flux_lw_top(i,j,k)
    FIA_out%flux_lh_top(i2,j2,k) = FIA_in%flux_lh_top(i,j,k)
    FIA_out%lprec_top(i2,j2,k) = FIA_in%lprec_top(i,j,k)
    FIA_out%fprec_top(i2,j2,k) = FIA_in%fprec_top(i,j,k)
  enddo ; enddo ; enddo

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    FIA_out%tmelt(i2,j2,k) = FIA_in%tmelt(i,j,k)
    FIA_out%bmelt(i2,j2,k) = FIA_in%bmelt(i,j,k)
    FIA_out%sw_abs_ocn(i2,j2,k) = FIA_in%sw_abs_ocn(i,j,k)
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    FIA_out%bheat(i2,j2) = FIA_in%bheat(i,j)
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
    FIA_out%flux_sw_dn(i2,j2) = FIA_in%flux_sw_dn(i,j)
  enddo ; enddo
  !   FIA%flux_u_top and flux_v_top are deliberately not being copied, as they
  ! are only needed on the fast_ice_PEs
  !   FIA%frazil_left is deliberately not being copied, as it is only valid on
  ! the slow_ice_PEs.
  !   FIA%calving_preberg and FIA%calving_hflx_preberg are deliberately not
  ! being copied over.
  if (allocated(FIA_out%flux_t0)) then
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off
      FIA_out%flux_t0(i2,j2,k) = FIA_in%flux_t0(i,j,k)
      FIA_out%flux_q0(i2,j2,k) = FIA_in%flux_q0(i,j,k)
      FIA_out%flux_lw0(i2,j2,k) = FIA_in%flux_lw0(i,j,k)
      FIA_out%dhdt(i2,j2,k) = FIA_in%dhdt(i,j,k)
      FIA_out%dedt(i2,j2,k) = FIA_in%dedt(i,j,k)
      FIA_out%dlwdt(i2,j2,k) = FIA_in%dlwdt(i,j,k)
      FIA_out%Tskin_cat(i2,j2,k) = FIA_in%Tskin_cat(i,j,k)
    enddo ; enddo ; enddo
  endif

  if (FIA_in%first_copy .or. FIA_out%first_copy) then ; if (FIA_in%num_tr_fluxes >= 0) then
    if (FIA_out%num_tr_fluxes < 0) then
      ! Allocate the tr_flux_top arrays to accommodate the size of the input
      ! fluxes.  This only occurs the first time FIA_out is copied from a fully
      ! initialized FIA_in.
      FIA_out%num_tr_fluxes = FIA_in%num_tr_fluxes
      if (FIA_out%num_tr_fluxes > 0) then
        isd = HI_out%isd ; ied = HI_out%ied ; jsd = HI_out%jsd ; jed = HI_out%jed
        allocate(FIA_out%tr_flux_top(isd:ied, jsd:jed, 0:ncat, FIA_out%num_tr_fluxes))
        FIA_out%tr_flux_top(:,:,:,:) = 0.0
      endif
    endif
    FIA_in%first_copy = .false. ; FIA_out%first_copy = .false.
  endif ; endif

  if (FIA_in%num_tr_fluxes >= 0) then
    if (FIA_in%num_tr_fluxes /= FIA_out%num_tr_fluxes) &
      call SIS_error(FATAL, "copy_FIA_to_FIA called with different num_tr_fluxes.")

    do n=1,FIA_in%num_tr_fluxes ; do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off
      FIA_out%tr_flux_top(i2,j2,k,n) = FIA_in%tr_flux_top(i,j,k,n)
    enddo ; enddo ; enddo ; enddo
  endif

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
  logical :: first_copy
  integer :: i, j, isd, ied, jsd, jed, ncat
  integer :: num_tr


  first_copy = .false.
  if (associated(FIA_out)) first_copy = FIA_out%first_copy
  if (associated(FIA_in)) first_copy = first_copy .or. FIA_in%first_copy

  if (first_copy) then
    ! Determine the number of fluxes.
    num_tr = 0 ; if (associated(FIA_in)) num_tr = FIA_in%num_tr_fluxes
    call max_across_PEs(num_tr)

    if (associated(FIA_out)) then
      if (.not. present(G_out)) &
        call SIS_error(FATAL, "redistribute_sFIA_to_sFIA called with an "//&
                              "associated FIA_out but without G_out.")
      if (.not. present(IG)) &
        call SIS_error(FATAL, "redistribute_sFIA_to_sFIA called with an "//&
                              "associated FIA_out but without IG.")
      FIA_out%num_tr_fluxes = num_tr
      if ((num_tr > 0) .and. .not.allocated(FIA_out%tr_flux_top)) then
        isd = G_out%isd ; ied = G_out%ied ; jsd = G_out%jsd ; jed = G_out%jed
        ncat = IG%CatIce
        allocate(FIA_out%tr_flux_top(isd:ied, jsd:jed, 0:ncat, num_tr))
        FIA_out%tr_flux_top(:,:,:,:) = 0.0
      endif
      FIA_out%first_copy = .false.
    endif

    if (associated(FIA_in)) FIA_in%first_copy = .false.
  endif

  !   FIA%flux_u_top and flux_v_top are deliberately not being copied, as they
  ! are only needed on the fast_ice_PEs
  !   FIA%frazil_left is deliberately not being copied, as it is only valid on
  ! the slow_ice_PEs.
  !   FIA%calving_preberg and FIA%calving_hflx_preberg are deliberately not
  ! being copied over.
  ! avg_count, atmos_winds, and the IDs are deliberately not being copied.

  if (associated(FIA_out) .and. associated(FIA_in)) then
    call mpp_redistribute(domain_in, FIA_in%flux_t_top, domain_out, &
                          FIA_out%flux_t_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_q_top, domain_out, &
                          FIA_out%flux_q_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_vis_dir_top, domain_out, &
                          FIA_out%flux_sw_vis_dir_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_vis_dif_top, domain_out, &
                          FIA_out%flux_sw_vis_dif_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_nir_dir_top, domain_out, &
                          FIA_out%flux_sw_nir_dir_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_nir_dif_top, domain_out, &
                          FIA_out%flux_sw_nir_dif_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_lw_top, domain_out, &
                          FIA_out%flux_lw_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_lh_top, domain_out, &
                          FIA_out%flux_lh_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%lprec_top, domain_out, &
                          FIA_out%lprec_top, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%fprec_top, domain_out, &
                          FIA_out%fprec_top, complete=.true.)

    call mpp_redistribute(domain_in, FIA_in%tmelt, domain_out, &
                          FIA_out%tmelt, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%bmelt, domain_out, &
                          FIA_out%bmelt, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%sw_abs_ocn, domain_out, &
                          FIA_out%sw_abs_ocn, complete=.true.)

    call mpp_redistribute(domain_in, FIA_in%bheat, domain_out, &
                          FIA_out%bheat, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%WindStr_x, domain_out, &
                          FIA_out%WindStr_x, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%WindStr_y, domain_out, &
                          FIA_out%WindStr_y, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%WindStr_ocn_x, domain_out, &
                          FIA_out%WindStr_ocn_x, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%WindStr_ocn_y, domain_out, &
                          FIA_out%WindStr_ocn_y, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%p_atm_surf, domain_out, &
                          FIA_out%p_atm_surf, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%runoff, domain_out, &
                          FIA_out%runoff, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%calving, domain_out, &
                          FIA_out%calving, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%runoff_hflx, domain_out, &
                          FIA_out%runoff_hflx, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%calving_hflx, domain_out, &
                          FIA_out%calving_hflx, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%Tskin_avg, domain_out, &
                          FIA_out%Tskin_avg, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%ice_free, domain_out, &
                          FIA_out%ice_free, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%ice_cover, domain_out, &
                          FIA_out%ice_cover, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_dn, domain_out, &
                          FIA_out%flux_sw_dn, complete=.true.)

    if (allocated(FIA_in%flux_t0)) then
      call mpp_redistribute(domain_in, FIA_in%flux_t0, domain_out, &
                            FIA_out%flux_t0, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%flux_q0, domain_out, &
                            FIA_out%flux_q0, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%flux_lw0, domain_out, &
                            FIA_out%flux_lw0, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%dhdt, domain_out, &
                            FIA_out%dhdt, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%dedt, domain_out, &
                            FIA_out%dedt, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%dlwdt, domain_out, &
                            FIA_out%dlwdt, complete=.true.)
      call mpp_redistribute(domain_in, FIA_in%Tskin_cat, domain_out, &
                            FIA_out%Tskin_cat, complete=.true.)
    endif

    if (FIA_in%num_tr_fluxes > 0) then
      call mpp_redistribute(domain_in, FIA_in%tr_flux_top, domain_out, &
                            FIA_out%tr_flux_top)
    endif
  elseif (associated(FIA_out)) then
    ! Use the null pointers in place of the unneeded input arrays.
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%flux_t_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%flux_q_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%flux_sw_vis_dir_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%flux_sw_vis_dif_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%flux_sw_nir_dir_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%flux_sw_nir_dif_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%flux_lw_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%flux_lh_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%lprec_top, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%fprec_top, complete=.true.)

    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%tmelt, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%bmelt, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                          FIA_out%sw_abs_ocn, complete=.true.)

    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%bheat, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%WindStr_x, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%WindStr_y, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%WindStr_ocn_x, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%WindStr_ocn_y, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%p_atm_surf, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%runoff, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%calving, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%runoff_hflx, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%calving_hflx, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%Tskin_avg, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%ice_free, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%ice_cover, complete=.false.)
    call mpp_redistribute(domain_in, null_ptr2D, domain_out, &
                          FIA_out%flux_sw_dn, complete=.true.)

    if (allocated(FIA_out%flux_t0)) then
      call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                            FIA_out%flux_t0, complete=.false.)
      call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                            FIA_out%flux_q0, complete=.false.)
      call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                            FIA_out%flux_lw0, complete=.false.)
      call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                            FIA_out%dhdt, complete=.false.)
      call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                            FIA_out%dedt, complete=.false.)
      call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                            FIA_out%dlwdt, complete=.true.)
      call mpp_redistribute(domain_in, null_ptr3D, domain_out, &
                            FIA_out%Tskin_cat, complete=.true.)
    endif


    if (FIA_out%num_tr_fluxes > 0) then
      call mpp_redistribute(domain_in, null_ptr4D, domain_out, &
                            FIA_out%tr_flux_top)
    endif
  elseif (associated(FIA_in)) then
    ! Use the null pointers in place of the unneeded output arrays.
    call mpp_redistribute(domain_in, FIA_in%flux_t_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_q_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_vis_dir_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_vis_dif_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_nir_dir_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_nir_dif_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_lw_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_lh_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%lprec_top, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%fprec_top, domain_out, &
                          null_ptr3D, complete=.true.)

    call mpp_redistribute(domain_in, FIA_in%tmelt, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%bmelt, domain_out, &
                          null_ptr3D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%sw_abs_ocn, domain_out, &
                          null_ptr3D, complete=.true.)

    call mpp_redistribute(domain_in, FIA_in%bheat, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%WindStr_x, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%WindStr_y, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%WindStr_ocn_x, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%WindStr_ocn_y, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%p_atm_surf, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%runoff, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%calving, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%runoff_hflx, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%calving_hflx, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%Tskin_avg, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%ice_free, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%ice_cover, domain_out, &
                          null_ptr2D, complete=.false.)
    call mpp_redistribute(domain_in, FIA_in%flux_sw_dn, domain_out, &
                          null_ptr2D, complete=.true.)

    if (allocated(FIA_in%flux_t0)) then
      call mpp_redistribute(domain_in, FIA_in%flux_t0, domain_out, &
                            null_ptr3D, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%flux_q0, domain_out, &
                            null_ptr3D, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%flux_lw0, domain_out, &
                            null_ptr3D, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%dhdt, domain_out, &
                            null_ptr3D, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%dedt, domain_out, &
                            null_ptr3D, complete=.false.)
      call mpp_redistribute(domain_in, FIA_in%dlwdt, domain_out, &
                            null_ptr3D, complete=.true.)
      call mpp_redistribute(domain_in, FIA_in%Tskin_cat, domain_out, &
                            null_ptr3D, complete=.true.)
    endif

    if (FIA_in%num_tr_fluxes > 0) then
      call mpp_redistribute(domain_in, FIA_in%tr_flux_top, domain_out, &
                            null_ptr4D)
    endif
  else
    call SIS_error(FATAL, "redistribute_FIA_to_FIA called with "//&
                          "neither FIA_in nor FIA_out associated.")
  endif

end subroutine redistribute_FIA_to_FIA


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_IST_arrays deallocates the arrays in an ice_state_type.
subroutine dealloc_IST_arrays(IST)
  type(ice_state_type), intent(inout) :: IST

  deallocate(IST%part_size, IST%mH_snow, IST%mH_ice)
  deallocate(IST%mH_pond) ! mw/new
  deallocate(IST%enth_snow, IST%enth_ice, IST%sal_ice)
  if (allocated(IST%t_surf)) deallocate(IST%t_surf)

  if (allocated(IST%u_ice_C)) deallocate(IST%u_ice_C)
  if (allocated(IST%v_ice_C)) deallocate(IST%v_ice_C)
  if (allocated(IST%u_ice_B)) deallocate(IST%u_ice_B)
  if (allocated(IST%v_ice_B)) deallocate(IST%v_ice_B)

end subroutine dealloc_IST_arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_ocean_sfc_state deallocates the arrays in an ocean_sfc_state_type.
subroutine dealloc_ocean_sfc_state(OSS)
  type(ocean_sfc_state_type), pointer :: OSS

  if (.not.associated(OSS)) then
    call SIS_error(WARNING, "dealloc_ocean_sfc_state called with an unassociated pointer.")
    return
  endif

  deallocate(OSS%s_surf, OSS%SST_C, OSS%sea_lev, OSS%T_fr_ocn, OSS%frazil)
  if (allocated(OSS%u_ocn_B)) deallocate(OSS%u_ocn_B)
  if (allocated(OSS%v_ocn_B)) deallocate(OSS%v_ocn_B)
  if (allocated(OSS%u_ocn_C)) deallocate(OSS%u_ocn_C)
  if (allocated(OSS%v_ocn_C)) deallocate(OSS%v_ocn_C)

  deallocate(OSS)
end subroutine dealloc_ocean_sfc_state

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_simple_OSS deallocates the arrays in a simple_OSS_type.
subroutine dealloc_simple_OSS(OSS)
  type(simple_OSS_type), pointer :: OSS

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
  type(fast_ice_avg_type), pointer    :: FIA

  if (.not.associated(FIA)) then
    call SIS_error(WARNING, "dealloc_fast_ice_avg called with an unassociated pointer.")
    return
  endif

  deallocate(FIA%flux_u_top, FIA%flux_v_top )
  deallocate(FIA%flux_t_top, FIA%flux_q_top, FIA%flux_lw_top)
  deallocate(FIA%flux_lh_top, FIA%lprec_top, FIA%fprec_top)
  deallocate(FIA%flux_sw_vis_dir_top, FIA%flux_sw_vis_dif_top)
  deallocate(FIA%flux_sw_nir_dir_top, FIA%flux_sw_nir_dif_top)
  deallocate(FIA%runoff, FIA%calving, FIA%runoff_hflx, FIA%calving_hflx)
  deallocate(FIA%calving_preberg, FIA%calving_hflx_preberg)

  deallocate(FIA%bheat, FIA%tmelt, FIA%bmelt, FIA%frazil_left)
  deallocate(FIA%WindStr_x, FIA%WindStr_y, FIA%p_atm_surf)
  deallocate(FIA%WindStr_ocn_x, FIA%WindStr_ocn_y)
  deallocate(FIA%ice_free, FIA%ice_cover, FIA%sw_abs_ocn, FIA%Tskin_avg)

  if (allocated(FIA%flux_t0)) deallocate(FIA%flux_t0)
  if (allocated(FIA%flux_q0)) deallocate(FIA%flux_q0)
  if (allocated(FIA%flux_lw0)) deallocate(FIA%flux_lw0)
  if (allocated(FIA%dhdt))  deallocate(FIA%dhdt)
  if (allocated(FIA%dedt))  deallocate(FIA%dedt)
  if (allocated(FIA%dlwdt)) deallocate(FIA%dlwdt)
  if (allocated(FIA%Tskin_cat)) deallocate(FIA%Tskin_cat)

  deallocate(FIA)
end subroutine dealloc_fast_ice_avg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_total_sfc_flux deallocates the arrays in a total_sfc_flux_type.
subroutine dealloc_total_sfc_flux(TSF)
  type(total_sfc_flux_type), pointer    :: TSF

  if (.not.associated(TSF)) then
    call SIS_error(WARNING, "dealloc_total_sfc_flux called with an unassociated pointer.")
    return
  endif

  deallocate(TSF%flux_u, TSF%flux_v, TSF%flux_t, TSF%flux_q)
  deallocate(TSF%flux_sw_vis_dir, TSF%flux_sw_vis_dif)
  deallocate(TSF%flux_sw_nir_dir, TSF%flux_sw_nir_dif)
  deallocate(TSF%flux_lw, TSF%flux_lh, TSF%lprec, TSF%fprec)

end subroutine dealloc_total_sfc_flux

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_ice_rad deallocates the arrays in a ice_rad_type.
subroutine dealloc_ice_rad(Rad)
  type(ice_rad_type), pointer    :: Rad

  if (.not.associated(Rad)) then
    call SIS_error(WARNING, "dealloc_ice_rad called with an unassociated pointer.")
    return
  endif

  deallocate(Rad%sw_abs_sfc, Rad%sw_abs_snow, Rad%sw_abs_ice)
  deallocate(Rad%sw_abs_ocn, Rad%sw_abs_int)
  deallocate(Rad%coszen_nextrad)
  deallocate(Rad%T_skin)

  deallocate(Rad)
end subroutine dealloc_ice_rad

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_ice_ocean_flux deallocates the arrays in a ice_ocean_flux_type.
subroutine dealloc_ice_ocean_flux(IOF)
  type(ice_ocean_flux_type), pointer    :: IOF

  if (.not.associated(IOF)) then
    call SIS_error(WARNING, "dealloc_ice_ocean_flux called with an unassociated pointer.")
    return
  endif

  deallocate(IOF%flux_t_ocn_top, IOF%flux_q_ocn_top)
  deallocate(IOF%flux_lw_ocn_top, IOF%flux_lh_ocn_top)
  deallocate(IOF%flux_sw_vis_dir_ocn, IOF%flux_sw_vis_dif_ocn)
  deallocate(IOF%flux_sw_nir_dir_ocn, IOF%flux_sw_nir_dif_ocn)
  deallocate(IOF%lprec_ocn_top, IOF%fprec_ocn_top)
  deallocate(IOF%flux_u_ocn, IOF%flux_v_ocn, IOF%flux_salt)

  deallocate(IOF%Enth_Mass_in_atm, IOF%Enth_Mass_out_atm)
  deallocate(IOF%Enth_Mass_in_ocn, IOF%Enth_Mass_out_ocn)

  !Deallocating iceberg fields
  if (associated(IOF%mass_berg)) deallocate(IOF%mass_berg)
  if (associated(IOF%ustar_berg)) deallocate(IOF%ustar_berg)
  if (associated(IOF%area_berg)) deallocate(IOF%area_berg)

  deallocate(IOF)
end subroutine dealloc_ice_ocean_flux

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays in an ice_ocean_flux_type.
subroutine IOF_chksum(mesg, IOF, G)
  character(len=*),          intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(ice_ocean_flux_type), intent(in) :: IOF   !< The structure whose arrays are being checksummed.
  type(SIS_hor_grid_type),   intent(inout) :: G  !< The ice-model's horizonal grid type.

  call hchksum(IOF%flux_salt, trim(mesg)//" IOF%flux_salt", G%HI)

  call hchksum(IOF%flux_t_ocn_top, trim(mesg)//"  IOF%flux_t_ocn_top", G%HI)
  call hchksum(IOF%flux_q_ocn_top, trim(mesg)//"  IOF%flux_q_ocn_top", G%HI)
  call hchksum(IOF%flux_lw_ocn_top, trim(mesg)//" IOF%flux_lw_ocn_top", G%HI)
  call hchksum(IOF%flux_lh_ocn_top, trim(mesg)//" IOF%flux_lh_ocn_top", G%HI)
  call hchksum(IOF%flux_sw_vis_dir_ocn, trim(mesg)//"  IOF%flux_sw_vis_dir_ocn", G%HI)
  call hchksum(IOF%flux_sw_vis_dif_ocn, trim(mesg)//"  IOF%flux_sw_vis_dif_ocn", G%HI)
  call hchksum(IOF%flux_sw_nir_dir_ocn, trim(mesg)//"  IOF%flux_sw_nir_dir_ocn", G%HI)
  call hchksum(IOF%flux_sw_nir_dif_ocn, trim(mesg)//"  IOF%flux_sw_nir_dif_ocn", G%HI)
  call hchksum(IOF%lprec_ocn_top, trim(mesg)//"  IOF%lprec_ocn_top", G%HI)
  call hchksum(IOF%fprec_ocn_top, trim(mesg)//"  IOF%fprec_ocn_top", G%HI)
  call hchksum(IOF%flux_u_ocn, trim(mesg)//"  IOF%flux_u_ocn", G%HI)
  call hchksum(IOF%flux_v_ocn, trim(mesg)//"  IOF%flux_v_ocn", G%HI)

  call hchksum(IOF%Enth_Mass_in_atm, trim(mesg)//" IOF%Enth_Mass_in_atm", G%HI)
  call hchksum(IOF%Enth_Mass_out_atm, trim(mesg)//" IOF%Enth_Mass_out_atm", G%HI)
  call hchksum(IOF%Enth_Mass_in_ocn, trim(mesg)//" IOF%Enth_Mass_in_ocn", G%HI)
  call hchksum(IOF%Enth_Mass_out_ocn, trim(mesg)//" IOF%Enth_Mass_out_ocn", G%HI)
end subroutine IOF_chksum

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays in a fast_ice_avg_type.
subroutine FIA_chksum(mesg, FIA, G, check_ocean)
  character(len=*),        intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(fast_ice_avg_type), intent(in) :: FIA   !< The structure whose arrays are being checksummed.
  type(SIS_hor_grid_type), intent(inout) :: G  !< The ice-model's horizonal grid type.
  logical, optional,       intent(in) :: check_ocean !< If present and true, check the fluxes to the ocean.

  call hchksum(FIA%flux_t_top(:,:,1:), trim(mesg)//" FIA%flux_t_top", G%HI)
  call hchksum(FIA%flux_q_top(:,:,1:), trim(mesg)//" FIA%flux_q_top", G%HI)
  call hchksum(FIA%flux_sw_vis_dir_top(:,:,1:), trim(mesg)//" FIA%flux_sw_vis_dir_top", G%HI)
  call hchksum(FIA%flux_sw_vis_dif_top(:,:,1:), trim(mesg)//" FIA%flux_sw_vis_dif_top", G%HI)
  call hchksum(FIA%flux_sw_nir_dir_top(:,:,1:), trim(mesg)//" FIA%flux_sw_nir_dir_top", G%HI)
  call hchksum(FIA%flux_sw_nir_dif_top(:,:,1:), trim(mesg)//" FIA%flux_sw_nir_dif_top", G%HI)
  call hchksum(FIA%flux_lw_top(:,:,1:), trim(mesg)//" FIA%flux_lw_top", G%HI)
  call hchksum(FIA%flux_lh_top(:,:,1:), trim(mesg)//" FIA%flux_lh_top", G%HI)
  call hchksum(FIA%lprec_top(:,:,1:), trim(mesg)//" FIA%lprec_top", G%HI)
  call hchksum(FIA%fprec_top(:,:,1:), trim(mesg)//" FIA%fprec_top", G%HI)

  if (present(check_ocean)) then ; if (check_ocean) then
    call hchksum(FIA%flux_t_top(:,:,0), trim(mesg)//" FIA%flux_t_top0", G%HI)
    call hchksum(FIA%flux_q_top(:,:,0), trim(mesg)//" FIA%flux_q_top0", G%HI)
    call hchksum(FIA%flux_sw_vis_dir_top(:,:,0), trim(mesg)//" FIA%flux_sw_vis_dir_top0", G%HI)
    call hchksum(FIA%flux_sw_vis_dif_top(:,:,0), trim(mesg)//" FIA%flux_sw_vis_dif_top0", G%HI)
    call hchksum(FIA%flux_sw_nir_dir_top(:,:,0), trim(mesg)//" FIA%flux_sw_nir_dir_top0", G%HI)
    call hchksum(FIA%flux_sw_nir_dif_top(:,:,0), trim(mesg)//" FIA%flux_sw_nir_dif_top0", G%HI)
    call hchksum(FIA%flux_lw_top(:,:,0), trim(mesg)//" FIA%flux_lw_top0", G%HI)
    call hchksum(FIA%flux_lh_top(:,:,0), trim(mesg)//" FIA%flux_lh_top0", G%HI)
    call hchksum(FIA%lprec_top(:,:,0), trim(mesg)//" FIA%lprec_top0", G%HI)
    call hchksum(FIA%fprec_top(:,:,0), trim(mesg)//" FIA%fprec_top0", G%HI)
  endif ; endif

  call hchksum(FIA%tmelt, trim(mesg)//" FIA%tmelt", G%HI)
  call hchksum(FIA%bmelt, trim(mesg)//" FIA%bmelt", G%HI)
  call hchksum(FIA%sw_abs_ocn, trim(mesg)//" FIA%sw_abs_ocn", G%HI)

  call hchksum(FIA%bheat, trim(mesg)//" FIA%bheat", G%HI)
  call hchksum(FIA%WindStr_x, trim(mesg)//" FIA%WindStr_x", G%HI)
  call hchksum(FIA%WindStr_y, trim(mesg)//" FIA%WindStr_y", G%HI)
  call hchksum(FIA%WindStr_ocn_x, trim(mesg)//" FIA%WindStr_ocn_x", G%HI)
  call hchksum(FIA%WindStr_ocn_y, trim(mesg)//" FIA%WindStr_ocn_y", G%HI)
  call hchksum(FIA%p_atm_surf, trim(mesg)//" FIA%p_atm_surf", G%HI)
  call hchksum(FIA%runoff, trim(mesg)//" FIA%runoff", G%HI)
  call hchksum(FIA%calving, trim(mesg)//" FIA%calving", G%HI)
  call hchksum(FIA%runoff_hflx, trim(mesg)//" FIA%runoff_hflx", G%HI)
  call hchksum(FIA%calving_hflx, trim(mesg)//" FIA%calving_hflx", G%HI)
  call hchksum(FIA%ice_free, trim(mesg)//" FIA%ice_free", G%HI)
  call hchksum(FIA%ice_cover, trim(mesg)//" FIA%ice_cover", G%HI)
  call hchksum(FIA%flux_sw_dn, trim(mesg)//" FIA%flux_sw_dn", G%HI)

end subroutine FIA_chksum

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays in an ice_state_type.
subroutine IST_chksum(mesg, IST, G, IG, haloshift)
  character(len=*),        intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(ice_state_type),    intent(in) :: IST   !< The structure whose arrays are being checksummed.
  type(SIS_hor_grid_type), intent(inout) :: G  !< The ice-model's horizonal grid type.
  type(ice_grid_type),     intent(in) :: IG    !< The sea-ice grid type.
  integer, optional,       intent(in) :: haloshift !< The width of halos to check, or 0 if missing.
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      IST - The ice state type variable to be checked.
!  (in)      G - The ocean's grid structure.  (Inout due to halo updates.)
!  (in,opt)  haloshift - If present, check halo points out this far.
  character(len=20) :: k_str1, k_str
  integer :: hs, k

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  hs=0; if (present(haloshift)) hs=haloshift

  call hchksum(IST%part_size(:,:,0), trim(mesg)//" IST%part_size(0)", G%HI, haloshift=hs)
  call hchksum(IST%part_size(:,:,1:), trim(mesg)//" IST%part_size", G%HI, haloshift=hs)
  call hchksum(IST%mH_ice*IG%H_to_kg_m2, trim(mesg)//" IST%mH_ice", G%HI, haloshift=hs)
  do k=1,IG%NkIce
    write(k_str1,'(I8)') k
    k_str = "("//trim(adjustl(k_str1))//")"
    call hchksum(IST%enth_ice(:,:,:,k), trim(mesg)//" IST%enth_ice("//trim(k_str), G%HI, haloshift=hs)
    call hchksum(IST%sal_ice(:,:,:,k), trim(mesg)//" IST%sal_ice("//trim(k_str), G%HI, haloshift=hs)
  enddo
  call hchksum(IST%mH_snow*IG%H_to_kg_m2, trim(mesg)//" IST%mH_snow", G%HI, haloshift=hs)
  call hchksum(IST%enth_snow(:,:,:,1), trim(mesg)//" IST%enth_snow", G%HI, haloshift=hs)

  if (allocated(IST%u_ice_B) .and. allocated(IST%v_ice_B)) then
    if (allocated(IST%u_ice_B)) call Bchksum(IST%u_ice_B, mesg//" IST%u_ice_B", G%HI, haloshift=hs)
    if (allocated(IST%v_ice_B)) call Bchksum(IST%v_ice_B, mesg//" IST%v_ice_B", G%HI, haloshift=hs)
    call check_redundant_B(mesg//" IST%u/v_ice", IST%u_ice_B, IST%v_ice_B, G)
  endif
  if (allocated(IST%u_ice_C) .and. allocated(IST%v_ice_C)) then
    call uchksum(IST%u_ice_C, mesg//" IST%u_ice_C", G%HI, haloshift=hs)
    call vchksum(IST%v_ice_C, mesg//" IST%v_ice_C", G%HI, haloshift=hs)
    call check_redundant_C(mesg//" IST%u/v_ice_C", IST%u_ice_C, IST%v_ice_C, G)
  endif

end subroutine IST_chksum

subroutine IST_bounds_check(IST, G, IG, msg, OSS, Rad)
  type(ice_state_type),    intent(in)    :: IST
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(in)    :: IG
  character(len=*),        intent(in)    :: msg
  type(ocean_sfc_state_type), optional, intent(in) :: OSS
  type(ice_rad_type),         optional, intent(in) :: Rad

  character(len=512) :: mesg1, mesg2
  character(len=24) :: err
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: sum_part_sz
  real, dimension(IG%NkIce) :: S_col
  real    :: tsurf_min, tsurf_max, tice_min, tice_max, tOcn_min, tOcn_max
  real    :: enth_min, enth_max, m_max
  logical :: spec_thermo_sal
  integer :: i, j, k, m, isc, iec, jsc, jec, ncat, NkIce, i_off, j_off
  integer :: n_bad, i_bad, j_bad, k_bad

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  n_bad = 0 ; i_bad = 0 ; j_bad = 0 ; k_bad = 0 ; err = ":"

  m_max = 1.0e6*IG%kg_m2_to_H

  sum_part_sz(:,:) = 0.0
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    sum_part_sz(i,j) = sum_part_sz(i,j) + IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  tOcn_min = -100. ; tOcn_max = 60.
  if (present(OSS)) then
    do j=jsc,jec ; do i=isc,iec
      if ((OSS%s_surf(i,j) < 0.0) .or. (OSS%s_surf(i,j) > 100.0) .or. &
          (OSS%SST_C(i,j) < tOcn_min) .or. (OSS%SST_C(i,j) > tOcn_max)) then
        n_bad = n_bad + 1
        if (n_bad == 1) then ; i_bad = i ; j_bad = j ; err = "t_ocn" ; endif
      endif
    enddo ; enddo
  endif
  do j=jsc,jec ; do i=isc,iec
    if (abs(sum_part_sz(i,j) - 1.0) > 1.0e-5) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; err = "sum_part_sz" ; endif
    endif
  enddo ; enddo

  tsurf_min = tOcn_min ; tsurf_max = tOcn_max
  tice_min = -100. ; tice_max = 1.0
  enth_min = enth_from_TS(tice_min, 0., IST%ITV)
  enth_max = enth_from_TS(tice_max, 0., IST%ITV)
  if (present(Rad)) then
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      if ((Rad%t_skin(i,j,k) < tsurf_min) .or. (Rad%t_skin(i,j,k) > tsurf_max)) then
        n_bad = n_bad + 1
        if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "tsurf" ; endif
      endif
    enddo ; enddo ; enddo
  endif

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if ((IST%mH_ice(i,j,k) > m_max) .or. (IST%mH_snow(i,j,k) > m_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "large mass" ; endif
    endif
    if ((IST%enth_snow(i,j,k,1) < enth_min) .or. (IST%enth_snow(i,j,k,1) > enth_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "enth_snow" ; endif
    endif
  enddo ; enddo ; enddo

  do m=1,NkIce ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if ((IST%enth_ice(i,j,k,m) < enth_min) .or. (IST%enth_ice(i,j,k,m) > enth_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "enth_ice" ; endif
    endif
    if ((IST%sal_ice(i,j,k,m) < 0.0) .or. (IST%sal_ice(i,j,k,m) > 1000.0)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "sal_ice" ; endif
    endif
  enddo ; enddo ; enddo ; enddo

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
      write(mesg2,'("T_ocn = ",1pe12.4,", S_sfc = ",1pe12.4,", sum_ps = ",1pe12.4)') &
            OSS%SST_C(i,j), OSS%s_surf(i,j), sum_part_sz(i,j)
    else
      write(mesg2,'("sum_part_sz = ",1pe12.4)') sum_part_sz(i,j)
    endif
    call SIS_error(WARNING, "Bad ice state "//trim(err)//" "//trim(msg)//" ; "//trim(mesg1)//&
                            " ; "//trim(mesg2), all_print=.true.)
    if (k_bad > 0) then
      call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, &
                                 specified_thermo_salinity=spec_thermo_sal)
      if (.not.spec_thermo_sal) then
        do m=1,NkIce ; S_col(m) = IST%sal_ice(i,j,k,m) ; enddo
      endif
      write(mesg1,'("mi/ms = ", 2(1pe12.4)," ts = ",1pe12.4," ti = ",1pe12.4)') &
             IST%mH_ice(i,j,k)*IG%H_to_kg_m2, IST%mH_snow(i,j,k)*IG%H_to_kg_m2, &
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
