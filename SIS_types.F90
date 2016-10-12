!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_types contains a number of common SIS types, along with subroutines to   !
!   perform various tasks on these types, including allocation, deallocation,  !
!   registration for restarts, and checksums.                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_types

! use mpp_mod,          only: mpp_sum, stdout, input_nml_file, PE_here => mpp_pe
! use mpp_domains_mod,  only: domain2D, mpp_get_compute_domain, CORNER, EAST, NORTH
use mpp_domains_mod,  only: domain2D, CORNER, EAST, NORTH
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

use MOM_coms, only : PE_here
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg, is_root_pe
use MOM_file_parser, only : param_file_type
use MOM_hor_index,   only : hor_index_type
use SIS_diag_mediator, only : SIS_diag_ctrl, post_data=>post_SIS_data
use SIS_diag_mediator, only : register_SIS_diag_field, register_static_field
use MOM_checksums,      only : chksum, Bchksum, hchksum, uchksum, vchksum
use SIS_error_checking, only : check_redundant_B, check_redundant_C
use SIS_sum_output_type, only : SIS_sum_out_CS
use SIS_tracer_registry, only : SIS_tracer_registry_type
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS

implicit none ; private

#include <SIS2_memory.h>

public :: ice_state_type, ice_state_register_restarts, dealloc_IST_arrays
public :: IST_chksum, IST_bounds_check
public :: ice_ocean_flux_type, alloc_ice_ocean_flux, dealloc_ice_ocean_flux
public :: ocean_sfc_state_type, alloc_ocean_sfc_state, dealloc_ocean_sfc_state
public :: fast_ice_avg_type, alloc_fast_ice_avg, dealloc_fast_ice_avg
public :: ice_rad_type, ice_rad_register_restarts, dealloc_ice_rad
! public :: ice_diagnostics_init

public :: fast_thermo_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! This structure contains the ice model state, and is intended to be private   !
! to SIS2.  It is not to be shared with other components and modules, and may  !
! use different indexing conventions than other components.                    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type ice_state_type
  type(time_type) :: Time_Init, Time
  type(time_type) :: Time_step_fast, Time_step_slow

  ! The 8 of the following 10 variables constitute the sea-ice state.
  real, allocatable, dimension(:,:,:) :: &
    part_size   ! The fractional coverage of a grid cell by each ice
                ! thickness category, nondim, 0 to 1.  Category 0 is
                ! open ocean.  The sum of part_size is 1.
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

  ! State type
  logical :: slab_ice  ! If true, do the old style GFDL slab ice.
  ! State type
  logical :: Cgrid_dyn ! If true use a C-grid discretization of the
                       ! sea-ice dynamics.

  real :: Rho_ice      ! The nominal density of sea ice, in kg m-3.
  real :: Rho_snow     ! The nominal density of snow on sea ice, in kg m-3.

  logical :: specified_ice  ! If true, the sea ice is specified and there is
                            ! no need for ice dynamics.
!  logical :: bounds_check    ! If true, check for sensible values of thicknesses
                             ! temperatures, fluxes, etc.
!  logical :: debug           ! If true, write verbose checksums for debugging purposes.

  integer, dimension(:), allocatable :: id_t, id_sal
  integer :: id_cn=-1, id_hi=-1, id_hp = -1, id_hs=-1, id_tsn=-1, id_tsfc=-1, id_ext=-1 ! id_hp mw/new
  integer :: id_t_iceav=-1, id_s_iceav=-1, id_e2m=-1
  
  integer :: id_rdgr=-1 ! These do not exist yet: id_rdgf=-1, id_rdgo=-1, id_rdgv=-1

  integer :: id_slp=-1

  type(SIS_tracer_registry_type), pointer :: TrReg => NULL()

  type(ice_thermo_type), pointer  :: ITV => NULL()
  type(SIS2_ice_thm_CS), pointer  :: ice_thm_CSp => NULL()
  type(SIS_diag_ctrl)             :: diag ! A structure that regulates diagnostics.
end type ice_state_type

type fast_thermo_CS ! To be made ; private
  ! These two arrarys are used with column_check when evaluating the enthalpy
  ! conservation with the fast thermodynamics code. 
  real, pointer, dimension(:,:,:) :: &
    enth_prev, heat_in

  logical :: debug        ! If true, write verbose checksums for debugging purposes.
  logical :: column_check ! If true, enable the heat check column by column.
  real    :: imb_tol      ! The tolerance for imbalances to be flagged by
                          ! column_check, nondim.
  logical :: bounds_check ! If true, check for sensible values of thicknesses
                          ! temperatures, fluxes, etc.

  integer :: n_fast = 0   ! The number of times update_ice_model_fast
                          ! has been called.
end type fast_thermo_CS

!> ocean_sfc_state_type contains variables that describe the ocean's surface
!! state as seen by the sea-ice or atmosphere, on the ice grid.
type ocean_sfc_state_type
  ! 5 of the following 7 variables describe the ocean state as seen by the sea ice.
  real, allocatable, dimension(:,:) :: &
    s_surf , &  ! The ocean's surface salinity in g/kg.
    t_ocn  , &  ! The ocean's bulk surface temperature in degC.
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

!   type(coupler_3d_bc_type)   :: ocean_fields       ! array of fields used for additional tracers

  real :: kmelt ! A constant that is used in the calculation of the ocean/ice
                ! basal heat flux, in W m-2 K-1.  This could be replaced with
                ! an array reflecting the turbulence in the under-ice ocean
                ! boundary layer and the effective depth of the reported value
                ! of t_ocn.
 
  ! diagnostic IDs for ocean surface  properties.
  integer :: id_sst=-1, id_sss=-1, id_ssh=-1, id_uo=-1, id_vo=-1, id_frazil=-1
end type ocean_sfc_state_type

!> fast_ice_avg_type contains variables that describe the fluxes between the
!! atmosphere and the ice or that have been accumlated over fast thermodynamic
!! steps but will be applied to the slow (mass-changing) thermodynamics.  Some
!! of these are diagnostics, while others are averages of fluxes taken during
!! the fast ice thermodynamics and used during the slow ice thermodynamics or dynamics.
type fast_ice_avg_type
  ! These are the arrays that are averaged over the fast thermodynamics.  They
  ! are either used to communicate to the slow thermodynamics or diagnostics or
  ! both.
  integer :: avg_count  ! The number of times that surface fluxes to the ice
                        ! have been incremented.
  logical :: atmos_winds ! The wind stresses come directly from the atmosphere
                         ! model and have the wrong sign.
  real, allocatable, dimension(:,:,:) :: &
    ! The 3rd dimension in each of the following is ice thickness category.
    flux_u_top         , & ! The downward flux of zonal and meridional
    flux_v_top         , & ! momentum on an A-grid in ???.
    flux_t_top         , & ! The upward sensible heat flux at the ice top
                           ! in W m-2.
    flux_q_top         , & ! The upward evaporative moisture flux at
                           ! top of the ice, in kg m-2 s-1.
    flux_lw_top        , & ! The downward flux of longwave radiation at
                           ! the top of the ice, in W m-2.
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
    sw_abs_ocn      !  The fraction of the absorbed shortwave radiation that is
                    !  absorbed in the ocean, nondim and <=1.
                    !  Equivalent sw_abs_ocn fields are in both the fast_ice_avg_type
                    !  and the ice_rad_type because it is used as a part of the slow
                    !  thermodynamic updates.
  real, allocatable, dimension(:,:) :: &
    bheat      , & ! The upward diffusive heat flux from the ocean
                   ! to the ice at the base of the ice, in W m-2.
    frazil_left, & ! The frazil heat flux that has not yet been
                   ! consumed in making ice, in J m-2. This array
                   ! is decremented by the ice model as the heat
                   ! flux is used up.
    WindStr_x  , & ! The zonal wind stress averaged over the ice
                   ! categories on an A-grid, in Pa.
    WindStr_y  , & ! The meridional wind stress averaged over the
                   ! ice categories on an A-grid, in Pa.
    ice_free   , & ! The fractional open water used in calculating
                   ! WindStr_[xy]_A; nondimensional, between 0 & 1.
    ice_cover      ! The fractional ice coverage, summed across all
                   ! thickness categories, used in calculating
                   ! WindStr_[xy]_A; nondimensional, between 0 & 1.

  integer :: num_tr_fluxes = -1   ! The number of tracer flux fields
  real, allocatable, dimension(:,:,:,:) :: &
    tr_flux_top    ! An array of tracer fluxes at the top of the
                   ! sea ice.
  integer, allocatable, dimension(:,:) :: tr_flux_index

  integer :: id_sh=-1, id_lh=-1, id_sw=-1
  integer :: id_lw=-1, id_snofl=-1, id_rain=-1,  id_evap=-1
  integer :: id_sw_vis_dir=-1, id_sw_vis_dif=-1, id_sw_nir_dir=-1, id_sw_nir_dif=-1
  integer :: id_sw_vis=-1, id_sw_dir=-1, id_sw_dif=-1
  integer :: id_tmelt=-1, id_bmelt=-1, id_bheat=-1

end type fast_ice_avg_type

!> ice_rad_type contains variables that describe the absorption and reflection
!! of shortwave radiation in and around the sea ice.
type ice_rad_type

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
    runoff, &           ! Liquid runoff into the ocean, in kg m-2.
    calving, &          ! Calving of ice or runoff of frozen fresh
                        ! water into the ocean, in kg m-2.
    calving_preberg, &  ! Calving of ice or runoff of frozen fresh
                        ! water into the ocean, exclusive of any
                        ! iceberg contributions, in kg m-2.
    runoff_hflx, &      ! The heat flux associated with runoff, based
                        ! on the temperature difference relative to a
                        ! reference temperature, in ???.
    calving_hflx, &     ! The heat flux associated with calving, based
                        ! on the temperature difference relative to a
                        ! reference temperature, in ???.
    calving_hflx_preberg, & ! The heat flux associated with calving,
                        ! exclusive of any iceberg contributions, based on
                        ! the temperature difference relative to a
                        ! reference temperature, in ???.
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
  integer, allocatable, dimension(:,:) :: tr_flux_index

  ! diagnostic IDs for ice-to-ocean fluxes.
  integer :: id_runoff=-1, id_calving=-1, id_runoff_hflx=-1, id_calving_hflx=-1
  integer :: id_saltf=-1
  ! The following are diagnostic IDs for iceberg-related fields.  These are only
  ! used if the iceberg code is activated.
  integer ::  id_ustar_berg=-1, id_area_berg=-1, id_mass_berg=-1
end type ice_ocean_flux_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_state_register_restarts allocates the arrays in the ice_state_type
!!     and registers any variables in the ice state type that need to be included
!!     in the restart files.
subroutine ice_state_register_restarts(mpp_domain, HI, IG, param_file, IST, &
                                       Ice_restart, restart_file)
  type(domain2d),          intent(in)    :: mpp_domain
  type(hor_index_type),    intent(in)    :: HI
  type(ice_grid_type),     intent(in)    :: IG
  type(param_file_type),   intent(in)    :: param_file
  type(ice_state_type),    intent(inout) :: IST
  type(restart_file_type), intent(inout) :: Ice_restart
  character(len=*),        intent(in)    :: restart_file

  integer :: CatIce, NkIce, idr, n
  character(len=8) :: nstr
 
  CatIce = IG%CatIce ; NkIce = IG%NkIce
  allocate(IST%part_size(SZI_(HI), SZJ_(HI), 0:CatIce)) ; IST%part_size(:,:,:) = 0.0
  allocate(IST%mH_pond(SZI_(HI), SZJ_(HI), CatIce)) ; IST%mH_pond(:,:,:) = 0.0 !  mw/new
  allocate(IST%mH_snow(SZI_(HI), SZJ_(HI), CatIce)) ; IST%mH_snow(:,:,:) = 0.0
  allocate(IST%enth_snow(SZI_(HI), SZJ_(HI), CatIce, 1)) ; IST%enth_snow(:,:,:,:) = 0.0
  allocate(IST%mH_ice(SZI_(HI), SZJ_(HI), CatIce)) ; IST%mH_ice(:,:,:) = 0.0
  allocate(IST%enth_ice(SZI_(HI), SZJ_(HI), CatIce, NkIce)) ; IST%enth_ice(:,:,:,:) = 0.0
  allocate(IST%sal_ice(SZI_(HI), SZJ_(HI), CatIce, NkIce)) ; IST%sal_ice(:,:,:,:) = 0.0
  if (IST%Cgrid_dyn) then
    allocate(IST%u_ice_C(SZIB_(HI), SZJ_(HI))) ; IST%u_ice_C(:,:) = 0.0
    allocate(IST%v_ice_C(SZI_(HI), SZJB_(HI))) ; IST%v_ice_C(:,:) = 0.0
  else
    allocate(IST%u_ice_B(SZIB_(HI), SZJB_(HI))) ; IST%u_ice_B(:,:) = 0.0
    allocate(IST%v_ice_B(SZIB_(HI), SZJB_(HI))) ; IST%v_ice_B(:,:) = 0.0
  endif

  allocate(IST%t_surf(SZI_(HI), SZJ_(HI), 0:CatIce)) ; IST%t_surf(:,:,:) = 0.0 !X

  ! ### THESE ARE DIAGNOSTICS.  PERHAPS THEY SHOULD ONLY BE ALLOCATED IF USED.
  allocate(IST%rdg_mice(SZI_(HI), SZJ_(HI), CatIce)) ; IST%rdg_mice(:,:,:) = 0.0


  ! Now register some of these arrays to be read from the restart files.
  idr = register_restart_field(Ice_restart, restart_file, 'part_size', IST%part_size, domain=mpp_domain)
  idr = register_restart_field(Ice_restart, restart_file, 't_surf', IST%t_surf, &
                               domain=mpp_domain)
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

end subroutine ice_state_register_restarts


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_fast_ice_avg allocates and zeros out the arrays in a fast_ice_avg_type.
subroutine alloc_fast_ice_avg(FIA, HI, IG)
  type(fast_ice_avg_type), pointer    :: FIA
  type(hor_index_type),    intent(in) :: HI
  type(ice_grid_type),     intent(in) :: IG

  integer :: CatIce, NkIce

  if (.not.associated(FIA)) allocate(FIA)
  CatIce = IG%CatIce ; NkIce = IG%NkIce

  FIA%avg_count = 0
  allocate(FIA%flux_u_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ; FIA%flux_u_top(:,:,:) = 0.0 !NR
  allocate(FIA%flux_v_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ; FIA%flux_v_top(:,:,:) = 0.0 !NR
  allocate(FIA%flux_t_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ;  FIA%flux_t_top(:,:,:) = 0.0 !NI
  allocate(FIA%flux_q_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ;  FIA%flux_q_top(:,:,:) = 0.0 !NI
  allocate(FIA%flux_sw_vis_dir_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ; FIA%flux_sw_vis_dir_top(:,:,:) = 0.0 !NI
  allocate(FIA%flux_sw_vis_dif_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ; FIA%flux_sw_vis_dif_top(:,:,:) = 0.0 !NI
  allocate(FIA%flux_sw_nir_dir_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ; FIA%flux_sw_nir_dir_top(:,:,:) = 0.0 !NI
  allocate(FIA%flux_sw_nir_dif_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ; FIA%flux_sw_nir_dif_top(:,:,:) = 0.0 !NI
  allocate(FIA%flux_lw_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ; FIA%flux_lw_top(:,:,:) = 0.0 !NI
  allocate(FIA%flux_lh_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ; FIA%flux_lh_top(:,:,:) = 0.0 !NI
  allocate(FIA%lprec_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ;  FIA%lprec_top(:,:,:) = 0.0 !NI
  allocate(FIA%fprec_top(SZI_(HI), SZJ_(HI), 0:CatIce)) ;  FIA%fprec_top(:,:,:) = 0.0 !NI

  allocate(FIA%frazil_left(SZI_(HI), SZJ_(HI))) ; FIA%frazil_left(:,:) = 0.0 !NR
  allocate(FIA%bheat(SZI_(HI), SZJ_(HI))) ; FIA%bheat(:,:) = 0.0 !NI
  allocate(FIA%tmelt(SZI_(HI), SZJ_(HI), CatIce)) ; FIA%tmelt(:,:,:) = 0.0 !NR
  allocate(FIA%bmelt(SZI_(HI), SZJ_(HI), CatIce)) ; FIA%bmelt(:,:,:) = 0.0 !NR
  allocate(FIA%WindStr_x(SZI_(HI), SZJ_(HI))) ; FIA%WindStr_x(:,:) = 0.0 !NI
  allocate(FIA%WindStr_y(SZI_(HI), SZJ_(HI))) ; FIA%WindStr_y(:,:) = 0.0 !NI
  allocate(FIA%ice_free(SZI_(HI), SZJ_(HI)))  ; FIA%ice_free(:,:) = 0.0 !NI
  allocate(FIA%ice_cover(SZI_(HI), SZJ_(HI))) ; FIA%ice_cover(:,:) = 0.0 !NI

  allocate(FIA%sw_abs_ocn(SZI_(HI), SZJ_(HI), CatIce)) ; FIA%sw_abs_ocn(:,:,:) = 0.0

end subroutine alloc_fast_ice_avg


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

  integer :: CatIce, NkIce, idr

  if (.not.associated(Rad)) allocate(Rad)
  CatIce = IG%CatIce ; NkIce = IG%NkIce

  allocate(Rad%sw_abs_sfc(SZI_(HI), SZJ_(HI), CatIce)) ; Rad%sw_abs_sfc(:,:,:) = 0.0
  allocate(Rad%sw_abs_snow(SZI_(HI), SZJ_(HI), CatIce)) ; Rad%sw_abs_snow(:,:,:) = 0.0
  allocate(Rad%sw_abs_ice(SZI_(HI), SZJ_(HI), CatIce, NkIce)) ; Rad%sw_abs_ice(:,:,:,:) = 0.0
  allocate(Rad%sw_abs_ocn(SZI_(HI), SZJ_(HI), CatIce)) ; Rad%sw_abs_ocn(:,:,:) = 0.0
  allocate(Rad%sw_abs_int(SZI_(HI), SZJ_(HI), CatIce)) ; Rad%sw_abs_int(:,:,:) = 0.0

  allocate(Rad%coszen_nextrad(SZI_(HI), SZJ_(HI))) ; Rad%coszen_nextrad(:,:) = 0.0

  idr = register_restart_field(Ice_restart, restart_file, 'coszen', Rad%coszen_nextrad, &
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

  allocate(IOF%runoff(SZI_(HI), SZJ_(HI))) ; IOF%runoff(:,:) = 0.0 !NI
  allocate(IOF%calving(SZI_(HI), SZJ_(HI))) ; IOF%calving(:,:) = 0.0 !NI
  allocate(IOF%calving_preberg(SZI_(HI), SZJ_(HI))) ; IOF%calving_preberg(:,:) = 0.0 !NI, diag
  allocate(IOF%runoff_hflx(SZI_(HI), SZJ_(HI))) ; IOF%runoff_hflx(:,:) = 0.0 !NI
  allocate(IOF%calving_hflx(SZI_(HI), SZJ_(HI))) ; IOF%calving_hflx(:,:) = 0.0 !NI
  allocate(IOF%calving_hflx_preberg(SZI_(HI), SZJ_(HI))) ; IOF%calving_hflx_preberg(:,:) = 0.0 !NI, diag
  allocate(IOF%flux_salt(SZI_(HI), SZJ_(HI))) ; IOF%flux_salt(:,:) = 0.0 !NI

  allocate(IOF%flux_t_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%flux_t_ocn_top(:,:) = 0.0 !NI
  allocate(IOF%flux_q_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%flux_q_ocn_top(:,:) = 0.0 !NI
  allocate(IOF%flux_lw_ocn_top(SZI_(HI), SZJ_(HI))) ; IOF%flux_lw_ocn_top(:,:) = 0.0 !NI
  allocate(IOF%flux_lh_ocn_top(SZI_(HI), SZJ_(HI))) ; IOF%flux_lh_ocn_top(:,:) = 0.0 !NI
  allocate(IOF%flux_sw_vis_dir_ocn(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sw_vis_dir_ocn(:,:) = 0.0 !NI
  allocate(IOF%flux_sw_vis_dif_ocn(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sw_vis_dif_ocn(:,:) = 0.0 !NI
  allocate(IOF%flux_sw_nir_dir_ocn(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sw_nir_dir_ocn(:,:) = 0.0 !NI
  allocate(IOF%flux_sw_nir_dif_ocn(SZI_(HI), SZJ_(HI))) ;  IOF%flux_sw_nir_dif_ocn(:,:) = 0.0 !NI
  allocate(IOF%lprec_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%lprec_ocn_top(:,:) = 0.0 !NI
  allocate(IOF%fprec_ocn_top(SZI_(HI), SZJ_(HI))) ;  IOF%fprec_ocn_top(:,:) = 0.0 !NI
  allocate(IOF%flux_u_ocn(SZI_(HI), SZJ_(HI)))    ;  IOF%flux_u_ocn(:,:) = 0.0 !NI
  allocate(IOF%flux_v_ocn(SZI_(HI), SZJ_(HI)))    ;  IOF%flux_v_ocn(:,:) = 0.0 !NI

  allocate(IOF%Enth_Mass_in_atm(SZI_(HI), SZJ_(HI)))  ; IOF%Enth_Mass_in_atm(:,:) = 0.0 !NR
  allocate(IOF%Enth_Mass_out_atm(SZI_(HI), SZJ_(HI))) ; IOF%Enth_Mass_out_atm(:,:) = 0.0 !NR
  allocate(IOF%Enth_Mass_in_ocn(SZI_(HI), SZJ_(HI)))  ; IOF%Enth_Mass_in_ocn(:,:) = 0.0 !NR
  allocate(IOF%Enth_Mass_out_ocn(SZI_(HI), SZJ_(HI))) ; IOF%Enth_Mass_out_ocn(:,:) = 0.0 !NR

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

  allocate(OSS%s_surf(SZI_(HI), SZJ_(HI))) ; OSS%s_surf(:,:) = 0.0
  allocate(OSS%t_ocn(SZI_(HI), SZJ_(HI)))  ; OSS%t_ocn(:,:) = 0.0 
  allocate(OSS%sea_lev(SZI_(HI), SZJ_(HI))) ; OSS%sea_lev(:,:) = 0.0
  allocate(OSS%frazil(SZI_(HI), SZJ_(HI))) ; OSS%frazil(:,:) = 0.0


  if (Cgrid_dyn) then
    allocate(OSS%u_ocn_C(SZIB_(HI), SZJ_(HI))) ; OSS%u_ocn_C(:,:) = 0.0
    allocate(OSS%v_ocn_C(SZI_(HI), SZJB_(HI))) ; OSS%v_ocn_C(:,:) = 0.0
  else
    allocate(OSS%u_ocn_B(SZIB_(HI), SZJB_(HI))) ; OSS%u_ocn_B(:,:) = 0.0
    allocate(OSS%v_ocn_B(SZIB_(HI), SZJB_(HI))) ; OSS%v_ocn_B(:,:) = 0.0
  endif

end subroutine alloc_ocean_sfc_state

subroutine dealloc_IST_arrays(IST)
  type(ice_state_type), intent(inout) :: IST

  deallocate(IST%part_size, IST%mH_snow, IST%mH_ice)
  deallocate(IST%mH_pond) ! mw/new
  deallocate(IST%enth_snow, IST%enth_ice, IST%sal_ice, IST%t_surf)
  if (IST%Cgrid_dyn) then
    deallocate(IST%u_ice_C, IST%v_ice_C)
  else
    deallocate(IST%u_ice_B, IST%v_ice_B)
  endif

end subroutine dealloc_IST_arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_ocean_sfc_state deallocates the arrays in an ocean_sfc_state_type.
subroutine dealloc_ocean_sfc_state(OSS)
  type(ocean_sfc_state_type), pointer :: OSS

  if (.not.associated(OSS)) then
    call SIS_error(WARNING, "dealloc_ocean_sfc_state called with an unassociated pointer.")
    return
  endif

  deallocate(OSS%s_surf, OSS%t_ocn, OSS%sea_lev, OSS%frazil)
  if (allocated(OSS%u_ocn_B)) deallocate(OSS%u_ocn_B)
  if (allocated(OSS%v_ocn_B)) deallocate(OSS%v_ocn_B)
  if (allocated(OSS%u_ocn_C)) deallocate(OSS%u_ocn_C)
  if (allocated(OSS%v_ocn_C)) deallocate(OSS%v_ocn_C)

  deallocate(OSS)
end subroutine dealloc_ocean_sfc_state

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

  deallocate(FIA%bheat, FIA%tmelt, FIA%bmelt, FIA%frazil_left)
  deallocate(FIA%WindStr_x, FIA%WindStr_y, FIA%ice_free, FIA%ice_cover)
  deallocate(FIA%sw_abs_ocn)

  deallocate(FIA)
end subroutine dealloc_fast_ice_avg

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
  deallocate(IOF%runoff, IOF%calving, IOF%runoff_hflx, IOF%calving_hflx)
  deallocate(IOF%calving_preberg, IOF%calving_hflx_preberg)
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
subroutine IST_chksum(mesg, IST, G, IG, haloshift)
  character(len=*),        intent(in) :: mesg
  type(ice_state_type),    intent(in) :: IST
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(in) :: IG
  integer, optional,       intent(in) :: haloshift
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

  call hchksum(IST%part_size, trim(mesg)//" IST%part_size", G%HI, haloshift=hs)
  call hchksum(IST%mH_ice*IG%H_to_kg_m2, trim(mesg)//" IST%mH_ice", G%HI, haloshift=hs)
  do k=1,IG%NkIce
    write(k_str1,'(I8)') k
    k_str = "("//trim(adjustl(k_str1))//")"
    call hchksum(IST%enth_ice(:,:,:,k), trim(mesg)//" IST%enth_ice("//trim(k_str), G%HI, haloshift=hs)
    call hchksum(IST%sal_ice(:,:,:,k), trim(mesg)//" IST%sal_ice("//trim(k_str), G%HI, haloshift=hs)
  enddo
  call hchksum(IST%mH_snow*IG%H_to_kg_m2, trim(mesg)//" IST%mH_snow", G%HI, haloshift=hs)
  call hchksum(IST%enth_snow(:,:,:,1), trim(mesg)//" IST%enth_snow", G%HI, haloshift=hs)
  if (allocated(IST%u_ice_B)) call Bchksum(IST%u_ice_B, mesg//" IST%u_ice_B", G%HI, haloshift=hs)
  if (allocated(IST%v_ice_B)) call Bchksum(IST%v_ice_B, mesg//" IST%v_ice_B", G%HI, haloshift=hs)
  call check_redundant_B(mesg//" IST%u/v_ice", IST%u_ice_B, IST%v_ice_B, G)
  if (IST%Cgrid_dyn) then
    call uchksum(IST%u_ice_C, mesg//" IST%u_ice_C", G%HI, haloshift=hs)
    call vchksum(IST%v_ice_C, mesg//" IST%v_ice_C", G%HI, haloshift=hs)
    call check_redundant_C(mesg//" IST%u/v_ice_C", IST%u_ice_C, IST%v_ice_C, G)
  endif

end subroutine IST_chksum

subroutine IST_bounds_check(IST, G, IG, msg, OSS)
  type(ice_state_type),    intent(in)    :: IST
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(in)    :: IG
  character(len=*),        intent(in)    :: msg
  type(ocean_sfc_state_type), optional, intent(in) :: OSS

  character(len=512) :: mesg1, mesg2
  character(len=24) :: err
  real, dimension(SZI_(G),SZJ_(G)) :: sum_part_sz
  real, dimension(IG%NkIce) :: S_col
  real    :: tsurf_min, tsurf_max, tice_min, tice_max, tOcn_min, tOcn_max
  real    :: enth_min, enth_max, m_max
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
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
          (OSS%t_ocn(i,j) < tOcn_min) .or. (OSS%t_ocn(i,j) > tOcn_max)) then
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

  tsurf_min = tOcn_min + T_0degC ; tsurf_max = tOcn_max + T_0degC
  tice_min = -100. ; tice_max = 1.0
  enth_min = enth_from_TS(tice_min, 0., IST%ITV)
  enth_max = enth_from_TS(tice_max, 0., IST%ITV)
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    if ((IST%t_surf(i,j,k) < tsurf_min) .or. (IST%t_surf(i,j,k) > tsurf_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; err = "tsurf" ; endif
    endif
  enddo ; enddo ; enddo

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
      write(mesg2,'("T_sfc = ",1pe12.4,", ps = ",1pe12.4)') IST%t_surf(i,j,k), IST%part_size(i,j,k)
    elseif (present(OSS)) then
      write(mesg2,'("T_ocn = ",1pe12.4,", S_sfc = ",1pe12.4,", sum_ps = ",1pe12.4)') &
            OSS%t_ocn(i,j), OSS%s_surf(i,j), sum_part_sz(i,j)
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
