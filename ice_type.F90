!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_type_mod - maintains the sea ice data, reads/writes restarts, reads the  !
!                namelist and initializes diagnostics. - Mike Winton           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_type_mod

use mpp_mod,          only: mpp_sum, stdout, input_nml_file, PE_here => mpp_pe
use mpp_domains_mod,  only: domain2D, mpp_get_compute_domain, CORNER, EAST, NORTH
use mpp_parameter_mod, only: CGRID_NE, BGRID_NE, AGRID
use fms_mod,          only: open_namelist_file, check_nml_error, close_file
use fms_io_mod,       only: save_restart, restore_state, query_initialized
use fms_io_mod,       only: register_restart_field, restart_file_type
use time_manager_mod, only: time_type, time_type_to_real
use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type

use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type

use SIS_dyn_bgrid,    only: SIS_B_dyn_CS
use SIS_dyn_cgrid,    only: SIS_C_dyn_CS
use ice_transport_mod, only: ice_transport_CS
use SIS2_ice_thm, only : ice_thermo_type, SIS2_ice_thm_CS, enth_from_TS, energy_melt_EnthS
use SIS2_ice_thm, only : get_SIS2_thermo_coefs, temp_from_En_S
use ice_bergs, only: icebergs, icebergs_stock_pe, icebergs_save_restart

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

public :: ice_data_type, ice_state_type
public :: ice_ocean_flux_type, alloc_ice_ocean_flux, dealloc_ice_ocean_flux
public :: ocean_sfc_state_type, alloc_ocean_sfc_state, dealloc_ocean_sfc_state
public :: fast_ice_avg_type, alloc_fast_ice_avg, dealloc_fast_ice_avg
public :: ice_model_restart, dealloc_ice_arrays, dealloc_IST_arrays
public :: ice_data_type_register_restarts, ice_state_register_restarts
public :: ice_diagnostics_init, ice_stock_pe
public :: ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum, ice_data_type_chksum
public :: IST_chksum, Ice_public_type_chksum, Ice_public_type_bounds_check, IST_bounds_check

public :: fast_thermo_CS, slow_thermo_CS, dyn_trans_CS

real, parameter :: missing = -1e34
integer, parameter :: miss_int = -9999

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

  ! These two arrarys are used with column_check when evaluating the enthalpy
  ! conservation with the fast thermodynamics code. 
  real, pointer, dimension(:,:,:) :: &
    enth_prev, heat_in


  ! Shortwave absorption parameters that are set in ice_optics.
  real, allocatable, dimension(:,:,:) :: &
    sw_abs_sfc , &  ! The fractions of the absorbed shortwave radiation
    sw_abs_snow, &  ! that are absorbed in a surface skin layer (_sfc),
    sw_abs_ocn , &  ! the snow (_snow), by the ocean (_ocn), or integrated
    sw_abs_int      ! across all of the ice layers (_int), all nondim
                    ! and <=1.  sw_abs_int is only used for diagnostics.
                    ! Only sw_abs_ocn is used in the slow step.
  real, allocatable, dimension(:,:,:,:) :: &
    sw_abs_ice      ! The fraction of the absorbed shortwave that is
                    ! absorbed in each of the ice layers, nondim, <=1.

  real, allocatable, dimension(:,:)   :: &
    coszen_nextrad  ! Cosine of the solar zenith angle averaged
                    ! over the next radiation timestep, nondim.


  ! State type
  logical :: slab_ice  ! If true, do the old style GFDL slab ice.
  ! State type
  logical :: Cgrid_dyn ! If true use a C-grid discretization of the
                       ! sea-ice dynamics.
  ! Delete rho_ocean?
  real :: Rho_ocean    ! The nominal density of sea water, in kg m-3.
  real :: Rho_ice      ! The nominal density of sea ice, in kg m-3.
  real :: Rho_snow     ! The nominal density of snow on sea ice, in kg m-3.
  logical :: do_icebergs    ! If true, use the Lagrangian iceberg code, which
                            ! modifies the calving field among other things.
  ! SLOW THERMO (mostly)
  logical :: do_ridging     ! If true, use the ridging code

  logical :: specified_ice  ! If true, the sea ice is specified and there is
                            ! no need for ice dynamics.
  logical :: column_check   ! If true, enable the heat check column by column.
  real    :: imb_tol        ! The tolerance for imbalances to be flagged by
                            ! column_check, nondim.
  logical :: bounds_check    ! If true, check for sensible values of thicknesses
                             ! temperatures, fluxes, etc.
  logical :: debug           ! If true, write verbose checksums for debugging purposes.

  ! SLOW DYNAMICS?
  type(time_type) :: ice_stats_interval ! The interval between writes of the
                             ! globally summed ice statistics and conservation checks.
  type(time_type) :: write_ice_stats_time ! The next time to write out the ice statistics.

  ! FAST THERMO
  real :: kmelt          ! A constant that is used in the calculation of the
                         ! ocean/ice basal heat flux, in W m-2 K-1.

  ! Set_ocean_top
  logical :: slp2ocean  ! If true, apply sea level pressure to ocean surface.

  ! top level fast
  logical :: add_diurnal_sw ! If true, apply a synthetic diurnal cycle to the shortwave radiation.
  logical :: do_sun_angle_for_alb ! If true, find the sun angle for calculating
                                  ! the ocean albedo in the frame of the ice model.
  logical :: frequent_albedo_update ! If true, update the ice and ocean albedos
                                  ! within the fast ice model update.  Otherwise,
                                  ! the albedos are only updated within
                                  ! set_ice_surface_state.

!   type(coupler_3d_bc_type)   :: ocean_fields       ! array of fields used for additional tracers
!   type(coupler_2d_bc_type)   :: ocean_fluxes       ! array of fluxes used for additional tracers

  integer, dimension(:), allocatable :: id_t, id_sw_abs_ice, id_sal
  integer :: id_cn=-1, id_hi=-1, id_hs=-1, id_tsn=-1, id_tsfc=-1, id_ext=-1
  integer :: id_t_iceav=-1, id_s_iceav=-1, id_e2m=-1, id_swdn=-1, id_lwdn=-1
  
  integer :: id_rdgr=-1 ! These do not exist yet: id_rdgf=-1, id_rdgo=-1, id_rdgv=-1

  integer :: id_slp=-1
  !### THESE DIAGNOSTICS ARE NEVER SENT!
  ! integer :: id_ta=-1, id_obi=-1
  integer :: id_coszen=-1, id_alb=-1
  integer :: id_alb_vis_dir=-1, id_alb_vis_dif=-1, id_alb_nir_dir=-1, id_alb_nir_dif=-1
  integer :: id_sw_abs_sfc=-1, id_sw_abs_snow=-1, id_sw_pen=-1, id_sw_abs_ocn=-1

  !### THESE DIAGNOSTICS ARE NEVER SENT!
  ! integer :: id_strna=-1

  type(SIS_tracer_registry_type), pointer :: TrReg => NULL()
  type(SIS_tracer_flow_control_CS), pointer :: SIS_tracer_flow_CSp => NULL()

  type(fast_thermo_CS), pointer :: fast_thermo_CSp => NULL()
  type(slow_thermo_CS), pointer :: slow_thermo_CSp => NULL()
  type(dyn_trans_CS), pointer :: dyn_trans_CSp => NULL()

  type(ice_thermo_type), pointer  :: ITV => NULL()
  type(SIS2_ice_thm_CS), pointer  :: ice_thm_CSp => NULL()
  type(SIS_diag_ctrl)             :: diag ! A structure that regulates diagnostics.
!   type(icebergs), pointer     :: icebergs => NULL()
end type ice_state_type

type fast_thermo_CS ! To be made ; private
  logical :: debug        ! If true, write verbose checksums for debugging purposes.
  logical :: column_check ! If true, enable the heat check column by column.
  real    :: imb_tol      ! The tolerance for imbalances to be flagged by
                          ! column_check, nondim.
  logical :: bounds_check ! If true, check for sensible values of thicknesses
                          ! temperatures, fluxes, etc.

  integer :: n_fast = 0   ! The number of times update_ice_model_fast
                          ! has been called.
end type fast_thermo_CS

type slow_thermo_CS ! To be made ; private
  real :: ice_bulk_salin ! The globally constant sea ice bulk salinity, in g/kg
                         ! that is used to calculate the ocean salt flux.
  real :: ice_rel_salin  ! The initial bulk salinity of sea-ice relative to the
                         ! salinity of the water from which it formed, nondim.

  logical :: filling_frazil  ! If true, apply frazil to fill as many categories
                             ! as possible to fill in a uniform (minimum) amount
                             ! of frazil in all the thinnest categories.
                             ! Otherwise the frazil is always assigned to a
                             ! single category with part size > 0.01.
  real    :: fraz_fill_time  ! A timescale with which the filling frazil causes
                             ! the thinest cells to attain similar thicknesses,
                             ! or a negative number to apply the frazil flux
                             ! uniformly, in s.

  logical :: do_ice_restore ! If true, restore the sea-ice toward climatology
                            ! by applying a restorative heat flux.
  real    :: ice_restore_timescale ! The time scale for restoring ice when
                            ! do_ice_restore is true, in days.

  logical :: do_ice_limit   ! Limit the sea ice thickness to max_ice_limit.
  real    :: max_ice_limit  ! The maximum sea ice thickness, in m, when
                            ! do_ice_limit is true.

  logical :: nudge_sea_ice = .false. ! If true, nudge sea ice concentrations towards observations.
  real    :: nudge_sea_ice_rate = 0.0 ! The rate of cooling of ice-free water that
                              ! should be ice  covered in order to constrained the
                              ! ice concentration to track observations.  A suggested
                              ! value is of order 10000 W m-2.
  real    :: nudge_stab_fac   ! A factor that determines whether the buoyancy
                              ! flux associated with the sea ice nudging of
                              ! warm water includes a freshwater flux so as to
                              ! be destabilizing on net (<1), stabilizing (>1),
                              ! or neutral (=1).  The default is 1.
  real    :: nudge_conc_tol   ! The tolerance for mismatch in the sea ice concentations
                              ! before nudging begins to be applied.

  logical :: debug        ! If true, write verbose checksums for debugging purposes.
  logical :: column_check ! If true, enable the heat check column by column.
  real    :: imb_tol      ! The tolerance for imbalances to be flagged by
                          ! column_check, nondim.
  logical :: bounds_check ! If true, check for sensible values of thicknesses
                          ! temperatures, fluxes, etc.

  integer :: n_calls = 0  ! The number of times update_ice_model_slow_down
                          ! has been called.

  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                                   ! timing of diagnostic output.
  type(ice_transport_CS), pointer :: ice_transport_CSp => NULL()
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()

  integer :: id_qflim=-1, id_qfres=-1, id_fwnudge=-1
  integer :: id_lsrc=-1, id_lsnk=-1, id_bsnk=-1, id_sn2ic=-1
end type slow_thermo_CS

type dyn_trans_CS ! To be made ; private
  logical :: Cgrid_dyn ! If true use a C-grid discretization of the
                       ! sea-ice dynamics.
  real    :: dt_ice_dyn   ! The time step used for the slow ice dynamics, including
                          ! stepping the continuity equation and interactions
                          ! between the ice mass field and velocities, in s. If
                          ! 0 or negative, the coupling time step will be used.
  logical :: debug        ! If true, write verbose checksums for debugging purposes.
  logical :: column_check ! If true, enable the heat check column by column.
  real    :: imb_tol      ! The tolerance for imbalances to be flagged by
                          ! column_check, nondim.
  logical :: bounds_check ! If true, check for sensible values of thicknesses
                          ! temperatures, fluxes, etc.
  logical :: verbose      ! A flag to control the printing of an ice-diagnostic
                          ! message.  When true, this will slow the model down.

  integer :: ntrunc = 0   ! The number of times the velocity has been truncated
                          ! since the last call to write_ice_statistics.

  integer :: n_calls = 0  ! The number of times SIS_dynamics_trans has been called.

  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                                   ! timing of diagnostic output.

  integer :: id_fax=-1, id_fay=-1, id_xprt=-1, id_mib=-1, id_mi=-1
  type(SIS_B_dyn_CS), pointer     :: SIS_B_dyn_CSp => NULL()
  type(SIS_C_dyn_CS), pointer     :: SIS_C_dyn_CSp => NULL()
  type(ice_transport_CS), pointer :: ice_transport_CSp => NULL()
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()
  logical :: module_is_initialized = .false.
end type dyn_trans_CS

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
    bmelt                  ! Ice-bottom melting energy into the ice in J m-2.
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

  integer :: num_tr_fluxes = -1 ! The number of tracer flux fields
  real, allocatable, dimension(:,:,:) :: &
    tr_flux_ocn_top     ! An array of tracer fluxes at the ocean's surface.
  integer, allocatable, dimension(:,:) :: tr_flux_index

  ! diagnostic IDs for ice-to-ocean fluxes.
  integer :: id_runoff=-1, id_calving=-1, id_runoff_hflx=-1, id_calving_hflx=-1
  integer :: id_saltf=-1
end type ice_ocean_flux_type

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! This structure contains the ice model data (some used by calling routines);  !
! the third index is partition (1 is open water; 2 is ice cover)               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type ice_data_type !  ice_public_type
  type(domain2D)                     :: Domain
  type(time_type)                    :: Time
  logical                            :: pe
  integer, pointer, dimension(:)     :: pelist   =>NULL() ! Used for flux-exchange.
  logical, pointer, dimension(:,:)   :: ocean_pt =>NULL() ! An array that indicates all ocean points as true.

  ! These fields are used to provide information about the ice surface to the
  ! atmosphere, and contain separate values for each ice thickness category.
  real, pointer, dimension(:,:,:) :: &
    part_size => NULL(), &    ! The fractional coverage of a grid cell by each ice
                              ! thickness category, nondim, 0 to 1.  Category 1 is
                              ! open ocean.  The sum of part_size is 1.
    albedo    => NULL(), &    ! The surface albedo averaged across all wavelength
                              ! and orientation bands within each ice-thickness
                              ! category.  Nondimensional, between 0 and 1.
    albedo_vis_dir => NULL(), &  ! The surface albedos for visible (_vis) or
    albedo_nir_dir => NULL(), &  ! near-infrared (_nir) wavelengths of direct (_dir)
    albedo_vis_dif => NULL(), &  ! diffuse (_dif) shortwave radiation in each
    albedo_nir_dif => NULL(), &  ! ice-thickness category. Nondim, between 0 and 1.
    rough_mom   => NULL(), &  ! The roughnesses for momentum, heat, and moisture
    rough_heat  => NULL(), &  ! at the ocean surface, as provided by ocean_rough_mod,
    rough_moist => NULL(), &  ! apparently in m.
    t_surf      => NULL(), &  ! The surface temperature for the ocean or for
                              ! each ice-thickness category, in Kelvin.
    u_surf      => NULL(), &  ! The eastward (u_) and northward (v_) surface
    v_surf      => NULL()     ! velocities of the ocean (:,:,1) or sea-ice, in m s-1.
  real, pointer, dimension(:,:)   :: &
    s_surf         =>NULL()   ! The ocean's surface salinity, in g/kg.

  ! These arrays will be used to set the forcing for the ocean.
  real, pointer, dimension(:,:) :: &
    flux_u => NULL(), &   ! The flux of x-momentum into the ocean, in Pa.
    flux_v => NULL(), &   ! The flux of y-momentum into the ocean, in Pa.
    flux_t => NULL(), &   ! The flux of sensible heat out of the ocean, in W m-2.
    flux_q => NULL(), &   ! The evaporative moisture flux out of the ocean, in kg m-2 s-1.
    flux_lw => NULL(), &  ! The sensible heat flux out of the ocean, in W m-2.
    flux_sw_vis_dir => NULL(), &  ! The direct (dir) or diffuse (dif) shortwave
    flux_sw_vis_dif => NULL(), &  ! heat fluxes into the ocean in the visible
    flux_sw_nir_dir => NULL(), &  ! (vis) or near-infrared (nir) band, all
    flux_sw_nir_dif => NULL(), &  ! in W m-2.
    flux_lh => NULL(), &  ! The latent heat flux out of the ocean, in W m-2.
    lprec => NULL(), &    ! The liquid precipitation flux into the ocean, in kg m-2.
    fprec => NULL(), &    ! The frozen precipitation flux into the ocean, in kg m-2.
    p_surf => NULL(), &   ! The pressure at the ocean surface, in Pa.  This may
                          ! or may not include atmospheric pressure.
    runoff => NULL(), &   ! Liquid runoff into the ocean, in kg m-2.
    calving => NULL(), &  ! Calving of ice or runoff of frozen fresh water into
                          ! the ocean, in kg m-2.
    runoff_hflx => NULL(), &  ! The heat flux associated with runoff, based on
                              ! the temperature difference relative to a
                              ! reference temperature, in ???.
    calving_hflx => NULL(), & ! The heat flux associated with calving, based on
                              ! the temperature difference relative to a
                              ! reference temperature, in ???.
    flux_salt  => NULL()  ! The flux of salt out of the ocean in kg m-2.

  real, pointer, dimension(:,:) :: &
    area => NULL() , &    ! The area of ocean cells, in m2.  Land cells have
                          ! a value of 0, so this could also be used as a mask.
    mi   => NULL()        ! The total ice+snow mass, in kg m-2.
             ! mi is needed for the wave model. It is introduced here,
             ! because flux_ice_to_ocean cannot handle 3D fields. This may be
             ! removed, if the information on ice thickness can be derived from
             ! eventually from h_ice outside the ice module.
  integer, dimension(3)    :: axes
  type(coupler_3d_bc_type) :: ocean_fields       ! array of fields used for additional tracers
  type(coupler_2d_bc_type) :: ocean_fluxes       ! array of fluxes used for additional tracers
  type(coupler_3d_bc_type) :: ocean_fluxes_top   ! ###THIS IS ARCHAIC AND COULD BE DELETED!
  integer :: flux_uv_stagger = -999 ! The staggering relative to the tracer points
                    ! points of the two wind stress components. Valid entries
                    ! include AGRID, BGRID_NE, CGRID_NE, BGRID_SW, and CGRID_SW,
                    ! corresponding to the community-standard Arakawa notation.
                    ! (These are named integers taken from mpp_parameter_mod.)
                    ! Following SIS, this is BGRID_NE by default when the sea
                    ! ice is initialized, but here it is set to -999 so that a
                    ! global max across ice and non-ice processors can be used
                    ! to determine its value.

      type(icebergs), pointer     :: icebergs => NULL()
  type(SIS_hor_grid_type), pointer :: G => NULL() ! A structure containing metrics and grid info.
  type(ice_grid_type),  pointer :: IG => NULL() ! A structure containing sea-ice specific grid info.
  type(ice_state_type), pointer :: Ice_state => NULL() ! A structure containing the internal
                               ! representation of the ice state.
  type(ice_ocean_flux_type), pointer :: IOF => NULL()  ! A structure containing fluxes from
                               ! the ice to the ocean that are calculated by the ice model.
  type(ocean_sfc_state_type), pointer :: OSS => NULL() ! A structure containing the arrays
                               ! that describe the ocean's surface state, as it is revealed
                               ! to the ice model.
  type(fast_ice_avg_type), pointer :: FIA => NULL()    ! A structure of the fluxes and other
                               ! fields that are calculated during the fast ice step but
                               ! stored for later use by the slow ice step or the ocean.
  type(restart_file_type), pointer :: Ice_restart => NULL()
end type ice_data_type !  ice_public_type

!   The following three types are for data exchange with the FMS coupler
! they are defined here but declared in coupler_main and allocated in flux_init.
type :: ocean_ice_boundary_type
  real, dimension(:,:),   pointer :: &
    u      => NULL(), &  ! The x-direction ocean velocity at a position
                         ! determined by stagger, in m s-1.
    v      => NULL(), &  ! The y-direction ocean velocity at a position
                         ! determined by stagger, in m s-1.
    t      => NULL(), &  ! The ocean's surface temperature in Kelvin.
    s      => NULL(), &  ! The ocean's surface temperature in g/kg.
    frazil => NULL(), &  ! The frazil heat rejected by the ocean, in J m-2.
    sea_level => NULL()  ! The sea level after adjustment for any surface
                         ! pressure that the ocean allows to be expressed, in m.
  real, dimension(:,:,:), pointer :: data      =>NULL() ! collective field for "named" fields above
  integer                         :: stagger = BGRID_NE
  integer                         :: xtype              ! REGRID, REDIST or DIRECT used by coupler
  type(coupler_2d_bc_type)        :: fields     ! array of fields used for additional tracers
end type

type :: atmos_ice_boundary_type
  real, dimension(:,:,:), pointer :: &
    u_flux  => NULL(), & ! The true-eastward stresses (momentum fluxes) from the atmosphere
                         ! to the ocean or ice in each category, discretized on an A-grid,
                         ! and _not_ rotated to align with the model grid, in Pa.
    v_flux  => NULL(), & ! The true-northward stresses (momentum fluxes) from the atmosphere
                         ! to the ocean or ice in each category, discretized on an A-grid,
                         ! and _not_ rotated to align with the model grid, in Pa.
    u_star  => NULL(), & ! The atmospheric friction velocity on an A-grid, in Pa.
    t_flux  => NULL(), & ! The sensible heat flux flux from the ocean or ice into the
                         ! atmosphere at the surface, in W m-2.
    q_flux  => NULL(), & ! The flux of moisture from the ice or ocean to the
                         ! atmosphere due to evaporation or sublimation, in kg m-2 s-1.
    lw_flux => NULL(), & ! The flux longwave radiation from the atmosphere into the
                         ! ice or ocean, in W m-2.
    sw_flux_vis_dir => NULL(), &  ! The visible (_vis) or near-infrared (_nir),
    sw_flux_vis_dif => NULL(), &  ! direct (_dir) or diffuse (_dif) shortwave
    sw_flux_nir_dir => NULL(), &  ! radiation fluxes from the atmosphere into
    sw_flux_nir_dif => NULL(), &  ! the ice or ocean, in W m-2.
    lprec   => NULL(), & ! The liquid precipitation from the atmosphere onto the
                         ! atmosphere or ice in each thickness category, in kg m-2 s-1.
                         ! Rain falling on snow is currently assumed to pass or drain
                         ! directly through the ice into the ocean; this should be
                         ! revisitied!
    fprec   => NULL(), & ! The frozen precipitation (snowfall) from the atmosphere
                         ! to the ice or ocean, in kg m-2 s-1.  Currently in SIS2
                         ! all frozen precipitation, including snow, sleet, hail
                         ! or graupel, are all treated as snow.
    dhdt    => NULL(), & ! The derivative of the upward sensible heat flux with the
                         ! surface temperature in W m-2 K-1.
    dedt    => NULL(), & ! The derivative of the sublimation and evaporation rate
                         ! with the surface temperature, in kg m-2 s-1 K-1.
    drdt    => NULL(), & ! The derivative of the downward longwave radiative heat
                         ! flux with surface temperature, in W m-2 K-1.
    coszen  => NULL(), & ! The cosine of the solar zenith angle averged over the
                         ! next radiation timestep (not the one that was used to
                         ! calculate the sw_flux fields), nondim and <=1.
    p       => NULL(), & ! The atmospheric surface pressure, in Pa, often ~1e5 Pa.
    data    => NULL()
  integer                   :: xtype  ! DIRECT or REDIST - used by coupler.
  type(coupler_3d_bc_type)  :: fluxes ! array of fluxes used for additional tracers
end type

type :: land_ice_boundary_type
  real, dimension(:,:),   pointer :: &
    runoff  =>NULL(), &  ! The liquid runoff into the ocean, in kg m-2.
    calving =>NULL(), &  ! The frozen runoff into each cell, that is offered
                         ! first to the icebergs (if any), where it might be
                         ! used or modified before being passed to the ocean,
                         ! in kg m-2.
    runoff_hflx  =>NULL(), & ! The heat flux associated with the temperature of
                             ! of the liquid runoff, relative to liquid water
                             ! at 0 deg C, in W m-2.
    calving_hflx =>NULL()    ! The heat flux associated with the temperature of
                             ! of the frozen runoff, relative to liquid? (or frozen?) water
                             ! at 0 deg C, in W m-2.
  real, dimension(:,:,:), pointer :: data    =>NULL() ! collective field for "named" fields above
  integer                         :: xtype  ! REGRID, REDIST or DIRECT - used by coupler.
end type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_data_type_register_restarts - allocate the arrays in the ice_data_type   !
!     and register any variables in the ice data type that need to be included !
!     in the restart files.                                                    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_data_type_register_restarts(domain, CatIce, param_file, Ice, &
                                           Ice_restart, restart_file)
  type(domain2d),          intent(in)    :: domain
  integer,                 intent(in)    :: CatIce
  type(param_file_type),   intent(in)    :: param_file
  type(ice_data_type),     intent(inout) :: Ice
  type(restart_file_type), intent(inout) :: Ice_restart
  character(len=*),        intent(in)    :: restart_file

  ! This subroutine allocates the externally visible ice_data_type's arrays and
  ! registers the appopriate ones for inclusion in the restart file.
  integer :: isc, iec, jsc, jec, km, idr

  call mpp_get_compute_domain(domain, isc, iec, jsc, jec )
  km = CatIce + 1

  allocate(Ice%ocean_pt(isc:iec, jsc:jec)) ; Ice%ocean_pt(:,:) = .false. !derived
  allocate(Ice%t_surf(isc:iec, jsc:jec, km)) ; Ice%t_surf(:,:,:) = 0.0
  allocate(Ice%s_surf(isc:iec, jsc:jec)) ; Ice%s_surf(:,:) = 0.0 !NI
  allocate(Ice%u_surf(isc:iec, jsc:jec, km)) ; Ice%u_surf(:,:,:) = 0.0 !NI
  allocate(Ice%v_surf(isc:iec, jsc:jec, km)) ; Ice%v_surf(:,:,:) = 0.0 !NI
  allocate(Ice%part_size(isc:iec, jsc:jec, km)) ; Ice%part_size(:,:,:) = 0.0
  allocate(Ice%rough_mom(isc:iec, jsc:jec, km)) ; Ice%rough_mom(:,:,:) = 0.0
  allocate(Ice%rough_heat(isc:iec, jsc:jec, km)) ; Ice%rough_heat(:,:,:) = 0.0
  allocate(Ice%rough_moist(isc:iec, jsc:jec, km)) ; Ice%rough_moist(:,:,:) = 0.0

  allocate(Ice%albedo(isc:iec, jsc:jec, km)) ; Ice%albedo(:,:,:) = 0.0  ! Derived?
  allocate(Ice%albedo_vis_dir(isc:iec, jsc:jec, km)) ; Ice%albedo_vis_dir(:,:,:) = 0.0
  allocate(Ice%albedo_nir_dir(isc:iec, jsc:jec, km)) ; Ice%albedo_nir_dir(:,:,:) = 0.0
  allocate(Ice%albedo_vis_dif(isc:iec, jsc:jec, km)) ; Ice%albedo_vis_dif(:,:,:) = 0.0
  allocate(Ice%albedo_nir_dif(isc:iec, jsc:jec, km)) ; Ice%albedo_nir_dif(:,:,:) = 0.0

  allocate(Ice%flux_u(isc:iec, jsc:jec)) ; Ice%flux_u(:,:) = 0.0
  allocate(Ice%flux_v(isc:iec, jsc:jec)) ; Ice%flux_v(:,:) = 0.0
  allocate(Ice%flux_t(isc:iec, jsc:jec)) ; Ice%flux_t(:,:) = 0.0
  allocate(Ice%flux_q(isc:iec, jsc:jec)) ; Ice%flux_q(:,:) = 0.0
  allocate(Ice%flux_sw_vis_dir(isc:iec, jsc:jec)) ; Ice%flux_sw_vis_dir(:,:) = 0.0
  allocate(Ice%flux_sw_vis_dif(isc:iec, jsc:jec)) ; Ice%flux_sw_vis_dif(:,:) = 0.0
  allocate(Ice%flux_sw_nir_dir(isc:iec, jsc:jec)) ; Ice%flux_sw_nir_dir(:,:) = 0.0
  allocate(Ice%flux_sw_nir_dif(isc:iec, jsc:jec)) ; Ice%flux_sw_nir_dif(:,:) = 0.0
  allocate(Ice%flux_lw(isc:iec, jsc:jec)) ; Ice%flux_lw(:,:) = 0.0
  allocate(Ice%flux_lh(isc:iec, jsc:jec)) ; Ice%flux_lh(:,:) = 0.0 !NI
  allocate(Ice%lprec(isc:iec, jsc:jec)) ; Ice%lprec(:,:) = 0.0
  allocate(Ice%fprec(isc:iec, jsc:jec)) ; Ice%fprec(:,:) = 0.0
  allocate(Ice%p_surf(isc:iec, jsc:jec)) ; Ice%p_surf(:,:) = 0.0
  allocate(Ice%runoff(isc:iec, jsc:jec)) ; Ice%runoff(:,:) = 0.0
  allocate(Ice%calving(isc:iec, jsc:jec)) ; Ice%calving(:,:) = 0.0
  allocate(Ice%runoff_hflx(isc:iec, jsc:jec)) ; Ice%runoff_hflx(:,:) = 0.0
  allocate(Ice%calving_hflx(isc:iec, jsc:jec)) ; Ice%calving_hflx(:,:) = 0.0
  allocate(Ice%flux_salt(isc:iec, jsc:jec)) ; Ice%flux_salt(:,:) = 0.0

  allocate(Ice%area(isc:iec, jsc:jec)) ; Ice%area(:,:) = 0.0 !derived
  allocate(Ice%mi(isc:iec, jsc:jec)) ; Ice%mi(:,:) = 0.0 !NR


  ! Now register some of these arrays to be read from the restart files.
  idr = register_restart_field(Ice_restart, restart_file, 'rough_mom',   Ice%rough_mom,   domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'rough_heat',  Ice%rough_heat,  domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'rough_moist', Ice%rough_moist, domain=domain)

  idr = register_restart_field(Ice_restart, restart_file, 'flux_u',      Ice%flux_u,       domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_v',      Ice%flux_v,       domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_t',      Ice%flux_t,       domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_q',      Ice%flux_q,       domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_salt',   Ice%flux_salt,    domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_lw',     Ice%flux_lw,      domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'lprec',       Ice%lprec,        domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'fprec',       Ice%fprec,        domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'runoff',      Ice%runoff,       domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'calving',     Ice%calving,      domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'runoff_hflx', Ice%runoff_hflx,  domain=domain, mandatory=.false.)
  idr = register_restart_field(Ice_restart, restart_file, 'calving_hflx',Ice%calving_hflx, domain=domain, mandatory=.false.)
  idr = register_restart_field(Ice_restart, restart_file, 'p_surf',      Ice%p_surf,       domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_sw_vis_dir', Ice%flux_sw_vis_dir, &
                                      domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_sw_vis_dif', Ice%flux_sw_vis_dif, &
                                      domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_sw_nir_dir', Ice%flux_sw_nir_dir, &
                                      domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'flux_sw_nir_dif', Ice%flux_sw_nir_dif, &
                                      domain=domain)
end subroutine ice_data_type_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_state_register_restarts allocates the arrays in the ice_state_type
!!     and registers any variables in the ice data type that need to be included
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

  allocate(IST%coszen_nextrad(SZI_(HI), SZJ_(HI))) ; IST%coszen_nextrad(:,:) = 0.0 !NR X

  allocate(IST%sw_abs_sfc(SZI_(HI), SZJ_(HI), CatIce)) ; IST%sw_abs_sfc(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_snow(SZI_(HI), SZJ_(HI), CatIce)) ; IST%sw_abs_snow(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_ice(SZI_(HI), SZJ_(HI), CatIce, NkIce)) ; IST%sw_abs_ice(:,:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_ocn(SZI_(HI), SZJ_(HI), CatIce)) ; IST%sw_abs_ocn(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_int(SZI_(HI), SZJ_(HI), CatIce)) ; IST%sw_abs_int(:,:,:) = 0.0 !NR

  allocate(IST%enth_prev(SZI_(HI), SZJ_(HI), CatIce)) ; IST%enth_prev(:,:,:) = 0.0
  allocate(IST%heat_in(SZI_(HI), SZJ_(HI), CatIce)) ; IST%heat_in(:,:,:) = 0.0

  ! ### THESE ARE DIAGNOSTICS.  PERHAPS THEY SHOULD ONLY BE ALLOCATED IF USED.
  allocate(IST%rdg_mice(SZI_(HI), SZJ_(HI), CatIce)) ; IST%rdg_mice(:,:,:) = 0.0


  ! Now register some of these arrays to be read from the restart files.
  idr = register_restart_field(Ice_restart, restart_file, 'part_size', IST%part_size, domain=mpp_domain)
  idr = register_restart_field(Ice_restart, restart_file, 't_surf', IST%t_surf, &
                               domain=mpp_domain)
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
  idr = register_restart_field(Ice_restart, restart_file, 'coszen', IST%coszen_nextrad, &
                               domain=mpp_domain, mandatory=.false.)

end subroutine ice_state_register_restarts


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_fast_ice_avg allocates and zeros out the arrays in a fast_ice_avg_type.
subroutine alloc_fast_ice_avg(FIA, HI, IG)
  type(fast_ice_avg_type), pointer    :: FIA
  type(hor_index_type),    intent(in) :: HI
  type(ice_grid_type),     intent(in) :: IG

  integer :: CatIce

  if (.not.associated(FIA)) allocate(FIA)
  CatIce = IG%CatIce

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

end subroutine alloc_fast_ice_avg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> alloc_ice_ocean_flux allocates and zeros out the arrays in an ice_ocean_flux_type.
subroutine alloc_ice_ocean_flux(IOF, HI)
  type(ice_ocean_flux_type), pointer    :: IOF
  type(hor_index_type),      intent(in) :: HI

  integer :: CatIce

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

subroutine dealloc_Ice_arrays(Ice)
  type(ice_data_type), intent(inout) :: Ice

  deallocate(Ice%ocean_pt, Ice%t_surf, Ice%s_surf)
  deallocate(Ice%u_surf, Ice%v_surf, Ice%part_size)
  deallocate(Ice%rough_mom, Ice%rough_heat, Ice%rough_moist)
  deallocate(Ice%albedo, Ice%albedo_vis_dir, Ice%albedo_nir_dir)
  deallocate(Ice%albedo_vis_dif, Ice%albedo_nir_dif)

  deallocate(Ice%flux_u, Ice%flux_v, Ice%flux_t, Ice%flux_q, Ice%flux_lw)
  deallocate(Ice%flux_lh, Ice%lprec, Ice%fprec, Ice%p_surf, Ice%runoff)
  deallocate(Ice%calving, Ice%runoff_hflx, Ice%calving_hflx)
  deallocate(Ice%flux_salt)
  deallocate(Ice%flux_sw_vis_dir, Ice%flux_sw_vis_dif)
  deallocate(Ice%flux_sw_nir_dir, Ice%flux_sw_nir_dif)
  deallocate(Ice%area, Ice%mi)
end subroutine dealloc_Ice_arrays

subroutine dealloc_IST_arrays(IST)
  type(ice_state_type), intent(inout) :: IST

  deallocate(IST%part_size, IST%mH_snow, IST%mH_ice)
  deallocate(IST%enth_snow, IST%enth_ice, IST%sal_ice, IST%t_surf)
  if (IST%Cgrid_dyn) then
    deallocate(IST%u_ice_C, IST%v_ice_C)
  else
    deallocate(IST%u_ice_B, IST%v_ice_B)
  endif

  deallocate(IST%sw_abs_sfc, IST%sw_abs_snow, IST%sw_abs_ice)
  deallocate(IST%sw_abs_ocn, IST%sw_abs_int, IST%coszen_nextrad)

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

  deallocate(FIA)
end subroutine dealloc_fast_ice_avg

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

subroutine Ice_public_type_chksum(mesg, Ice)
  character(len=*),    intent(in) :: mesg
  type(ice_data_type), intent(in) :: Ice
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      Ice - An ice_data_type structure whose elements are to be
!                  checksummed.

  ! Note that the publicly visible ice_data_type has no halos, so it is not
  ! possible do check their values.

  call chksum(Ice%part_size, trim(mesg)//" Ice%part_size")
  call chksum(Ice%albedo, trim(mesg)//" Ice%albedo")
  call chksum(Ice%albedo_vis_dir, trim(mesg)//" Ice%albedo_vis_dir")
  call chksum(Ice%albedo_nir_dir, trim(mesg)//" Ice%albedo_nir_dir")
  call chksum(Ice%albedo_vis_dif, trim(mesg)//" Ice%albedo_vis_dif")
  call chksum(Ice%albedo_nir_dif, trim(mesg)//" Ice%albedo_nir_dif")
  call chksum(Ice%rough_mom, trim(mesg)//" Ice%rough_mom")
  call chksum(Ice%rough_mom, trim(mesg)//" Ice%rough_mom")
  call chksum(Ice%rough_moist, trim(mesg)//" Ice%rough_moist")

  call chksum(Ice%t_surf, trim(mesg)//" Ice%t_surf")
  call chksum(Ice%u_surf, trim(mesg)//" Ice%u_surf")
  call chksum(Ice%v_surf, trim(mesg)//" Ice%v_surf")
  call chksum(Ice%s_surf, trim(mesg)//" Ice%s_surf")

  call chksum(Ice%flux_u, trim(mesg)//" Ice%flux_u")
  call chksum(Ice%flux_v, trim(mesg)//" Ice%flux_v")
  call chksum(Ice%flux_t, trim(mesg)//" Ice%flux_t")
  call chksum(Ice%flux_q, trim(mesg)//" Ice%flux_q")
  call chksum(Ice%flux_lw, trim(mesg)//" Ice%flux_lw")
  call chksum(Ice%flux_sw_vis_dir, trim(mesg)//" Ice%flux_sw_vis_dir")
  call chksum(Ice%flux_sw_nir_dir, trim(mesg)//" Ice%flux_sw_nir_dir")
  call chksum(Ice%flux_sw_vis_dif, trim(mesg)//" Ice%flux_sw_vis_dif")
  call chksum(Ice%flux_sw_nir_dif, trim(mesg)//" Ice%flux_sw_nir_dif")
  call chksum(Ice%flux_lh, trim(mesg)//" Ice%flux_lh")
  call chksum(Ice%lprec, trim(mesg)//" Ice%lprec")
  call chksum(Ice%fprec, trim(mesg)//" Ice%fprec")
  call chksum(Ice%p_surf, trim(mesg)//" Ice%p_surf")
  call chksum(Ice%calving, trim(mesg)//" Ice%calving")
  call chksum(Ice%runoff, trim(mesg)//" Ice%runoff")

end subroutine Ice_public_type_chksum

subroutine Ice_public_type_bounds_check(Ice, G, msg)
  type(ice_data_type),     intent(in)    :: Ice
  type(SIS_hor_grid_type), intent(inout) :: G
  character(len=*),        intent(in)    :: msg

  character(len=512) :: mesg1, mesg2
  integer :: i, j, k, l, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  integer :: n_bad, i_bad, j_bad, k_bad
  real    :: t_min, t_max
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = Ice%IG%CatIce
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc

  n_bad = 0 ; i_bad = 0 ; j_bad = 0 ; k_bad = 0

  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    if ((Ice%s_surf(i2,j2) < 0.0) .or. (Ice%s_surf(i2,j2) > 100.0)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; endif
    endif
    if ((abs(Ice%flux_t(i2,j2)) > 1e4) .or. (abs(Ice%flux_lw(i2,j2)) > 1e4)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; endif
    endif
  enddo ; enddo
  t_min = T_0degC-100. ; t_max = T_0degC+60.
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    if ((Ice%t_surf(i2,j2,k2) < t_min) .or. (Ice%t_surf(i2,j2,k2) > t_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; endif
    endif
  enddo ; enddo ; enddo

  if (n_bad > 0) then
    i2 = i_bad+i_off ; j2 = j_bad+j_off ; k2 = k_bad+1
    write(mesg1,'(" at ", 2(F6.1)," or i,j,k = ",3i4,"; nbad = ",i6," on pe ",i4)') &
           G%geolonT(i_bad,j_bad), G%geolatT(i_bad,j_bad), i_bad, j_bad, k_bad, n_bad, pe_here()
    write(mesg2,'("T_sfc = ",1pe12.4,", ps = ",1pe12.4,", flux_t,lw,q = ",3(1pe12.4))') &
       Ice%t_surf(i2,j2,k2), Ice%part_size(i2,j2,k2), Ice%flux_t(i2,j2), Ice%flux_lw(i2,j2), Ice%flux_q(i2,j2)
    call SIS_error(WARNING, "Bad ice data "//trim(msg)//" ; "//trim(mesg1)//" ; "//trim(mesg2), all_print=.true.)
  endif

end subroutine Ice_public_type_bounds_check

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

!=======================================================================
! <SUBROUTINE NAME="ice_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ice_model_restart(Ice, time_stamp)
  type(ice_data_type), intent(inout) :: Ice
  character(len=*),    intent(in), optional :: time_stamp

  call save_restart(Ice%Ice_restart, time_stamp)
  call icebergs_save_restart(Ice%icebergs)

end subroutine ice_model_restart
! </SUBROUTINE>
!=======================================================================

subroutine ice_diagnostics_init(Ice, IST, IOF, OSS, FIA, G, diag, Time)
  type(ice_data_type),        intent(inout) :: Ice
  type(ice_state_type),       intent(inout) :: IST
  type(ice_ocean_flux_type),  intent(inout) :: IOF
  type(ocean_sfc_state_type), intent(inout) :: OSS
  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(SIS_diag_ctrl),        intent(in)    :: diag
  type(time_type),            intent(inout) :: Time

  real, parameter       :: missing = -1e34
  integer               :: id_geo_lon, id_geo_lat, id_sin_rot, id_cos_rot, id_cell_area
  logical               :: sent
  integer :: i, j, k, isc, iec, jsc, jec, n, nLay
  character(len=8) :: nstr

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  nLay = Ice%IG%NkIce

  Ice%axes(1:2) = diag%axesTc%handles(1:2)

  id_sin_rot   = register_static_field('ice_model', 'SINROT', diag%axesT1, &
                 '-SINROT,COSROT points north', 'none')
  id_cos_rot   = register_static_field('ice_model', 'COSROT', diag%axesT1, &
                 'COSROT,SINROT points east','none')
  id_geo_lon   = register_static_field('ice_model', 'GEOLON', diag%axesT1, 'longitude', &
                 'degrees')
  id_geo_lat   = register_static_field('ice_model', 'GEOLAT', diag%axesT1, 'latitude', &
                 'degrees')
  id_cell_area = register_static_field('ice_model', 'CELL_AREA', diag%axesT1, &
                 'cell area', 'sphere')

  IST%id_ext = register_SIS_diag_field('ice_model', 'EXT', diag%axesT1, Time, &
               'ice modeled', '0 or 1', missing_value=missing)
  IST%id_cn       = register_SIS_diag_field('ice_model', 'CN', diag%axesTc, Time, &
               'ice concentration', '0-1', missing_value=missing)
  IST%id_hs       = register_SIS_diag_field('ice_model', 'HS', diag%axesT1, Time, &
               'snow thickness', 'm-snow', missing_value=missing)
  IST%id_tsn      = register_SIS_diag_field('ice_model', 'TSN', diag%axesT1, Time, &
               'snow layer temperature', 'C',  missing_value=missing)
  IST%id_hi       = register_SIS_diag_field('ice_model', 'HI', diag%axesT1, Time, &
               'ice thickness', 'm-ice', missing_value=missing)

  IST%id_t_iceav = register_SIS_diag_field('ice_model', 'T_bulkice', diag%axesT1, Time, &
               'Volume-averaged ice temperature', 'C', missing_value=missing)
  IST%id_s_iceav = register_SIS_diag_field('ice_model', 'S_bulkice', diag%axesT1, Time, &
               'Volume-averaged ice salinity', 'g/kg', missing_value=missing)
  call safe_alloc_ids_1d(IST%id_t, nLay)
  call safe_alloc_ids_1d(IST%id_sal, nLay)
  do n=1,nLay
    write(nstr, '(I4)') n ; nstr = adjustl(nstr)
    IST%id_t(n)   = register_SIS_diag_field('ice_model', 'T'//trim(nstr), &
                 diag%axesT1, Time, 'ice layer '//trim(nstr)//' temperature', &
                 'C',  missing_value=missing)
    IST%id_sal(n)   = register_SIS_diag_field('ice_model', 'Sal'//trim(nstr), &
               diag%axesT1, Time, 'ice layer '//trim(nstr)//' salinity', &
               'g/kg',  missing_value=missing)
  enddo
  IST%id_tsfc     = register_SIS_diag_field('ice_model', 'TS', diag%axesT1, Time, &
               'surface temperature', 'C', missing_value=missing)
  FIA%id_sh       = register_SIS_diag_field('ice_model','SH' ,diag%axesT1, Time, &
               'sensible heat flux', 'W/m^2',  missing_value=missing)
  FIA%id_lh       = register_SIS_diag_field('ice_model','LH' ,diag%axesT1, Time, &
               'latent heat flux', 'W/m^2', missing_value=missing)
  FIA%id_sw       = register_SIS_diag_field('ice_model','SW' ,diag%axesT1, Time, &
               'short wave heat flux', 'W/m^2', missing_value=missing)
  FIA%id_lw       = register_SIS_diag_field('ice_model','LW' ,diag%axesT1, Time, &
               'long wave heat flux over ice', 'W/m^2', missing_value=missing)
  FIA%id_snofl    = register_SIS_diag_field('ice_model','SNOWFL' ,diag%axesT1, Time, &
               'rate of snow fall', 'kg/(m^2*s)', missing_value=missing)
  FIA%id_rain     = register_SIS_diag_field('ice_model','RAIN' ,diag%axesT1, Time, &
               'rate of rain fall', 'kg/(m^2*s)', missing_value=missing)
  IOF%id_runoff   = register_SIS_diag_field('ice_model','RUNOFF' ,diag%axesT1, Time, &
               'liquid runoff', 'kg/(m^2*s)', missing_value=missing)

  IOF%id_calving  = register_SIS_diag_field('ice_model','CALVING',diag%axesT1, Time, &
               'frozen runoff', 'kg/(m^2*s)', missing_value=missing)
  IOF%id_runoff_hflx   = register_SIS_diag_field('ice_model','RUNOFF_HFLX' ,diag%axesT1, Time, &
               'liquid runoff sensible heat flux', 'W/m^2', missing_value=missing)
  IOF%id_calving_hflx  = register_SIS_diag_field('ice_model','CALVING_HFLX',diag%axesT1, Time, &
               'frozen runoff sensible heat flux', 'W/m^2', missing_value=missing)
  FIA%id_evap     = register_SIS_diag_field('ice_model','EVAP',diag%axesT1, Time, &
               'evaporation', 'kg/(m^2*s)', missing_value=missing)
  IOF%id_saltf    = register_SIS_diag_field('ice_model','SALTF' ,diag%axesT1, Time, &
               'ice to ocean salt flux', 'kg/(m^2*s)', missing_value=missing)
  FIA%id_tmelt    = register_SIS_diag_field('ice_model','TMELT'  ,diag%axesT1, Time, &
               'upper surface melting energy flux', 'W/m^2', missing_value=missing)
  FIA%id_bmelt    = register_SIS_diag_field('ice_model','BMELT'  ,diag%axesT1, Time, &
               'bottom surface melting energy flux', 'W/m^2', missing_value=missing)
  FIA%id_bheat    = register_SIS_diag_field('ice_model','BHEAT'  ,diag%axesT1, Time, &
               'ocean to ice heat flux', 'W/m^2', missing_value=missing)
  IST%id_e2m      = register_SIS_diag_field('ice_model','E2MELT' ,diag%axesT1, Time, &
               'heat needed to melt ice', 'J/m^2', missing_value=missing)
  OSS%id_frazil   = register_SIS_diag_field('ice_model','FRAZIL' ,diag%axesT1, Time, &
               'energy flux of frazil formation', 'W/m^2', missing_value=missing)
  IST%id_alb      = register_SIS_diag_field('ice_model','ALB',diag%axesT1, Time, &
               'surface albedo','0-1', missing_value=missing )
  IST%id_coszen   = register_SIS_diag_field('ice_model','coszen',diag%axesT1, Time, &
               'cosine of the solar zenith angle for the next radiation step','-1:1', missing_value=missing )
  IST%id_sw_abs_sfc= register_SIS_diag_field('ice_model','sw_abs_sfc',diag%axesT1, Time, &
               'SW frac. abs. at the ice surface','0:1', missing_value=missing )
  IST%id_sw_abs_snow= register_SIS_diag_field('ice_model','sw_abs_snow',diag%axesT1, Time, &
               'SW frac. abs. in snow','0:1', missing_value=missing )

  call safe_alloc_ids_1d(IST%id_sw_abs_ice, nLay)
  do n=1,nLay
    write(nstr, '(I4)') n ; nstr = adjustl(nstr)
    IST%id_sw_abs_ice(n) = register_SIS_diag_field('ice_model','sw_abs_ice'//trim(nstr), &
                 diag%axesT1, Time, 'SW frac. abs. in ice layer '//trim(nstr), &
                 '0:1', missing_value=missing )
  enddo
  IST%id_sw_pen= register_SIS_diag_field('ice_model','sw_pen',diag%axesT1, Time, &
               'SW frac. pen. surf.','0:1', missing_value=missing )
  IST%id_sw_abs_ocn= register_SIS_diag_field('ice_model','sw_abs_ocn',diag%axesT1, Time, &
               'SW frac. sent to the ocean','0:1', missing_value=missing )


  IST%id_alb_vis_dir = register_SIS_diag_field('ice_model','alb_vis_dir',diag%axesT1, Time, &
               'ice surface albedo vis_dir','0-1', missing_value=missing )
  IST%id_alb_vis_dif = register_SIS_diag_field('ice_model','alb_vis_dif',diag%axesT1, Time, &
               'ice surface albedo vis_dif','0-1', missing_value=missing )
  IST%id_alb_nir_dir = register_SIS_diag_field('ice_model','alb_nir_dir',diag%axesT1, Time, &
               'ice surface albedo nir_dir','0-1', missing_value=missing )
  IST%id_alb_nir_dif = register_SIS_diag_field('ice_model','alb_nir_dif',diag%axesT1, Time, &
               'ice surface albedo nir_dif','0-1', missing_value=missing )

!### THIS DIAGNOSTIC IS MISSING.
!  IST%id_strna    = register_SIS_diag_field('ice_model','STRAIN_ANGLE', diag%axesT1,Time, &
!               'strain angle', 'none', missing_value=missing)
  if (IST%Cgrid_dyn) then
    OSS%id_uo     = register_SIS_diag_field('ice_model', 'UO', diag%axesCu1, Time, &
               'surface current - x component', 'm/s', missing_value=missing)
    OSS%id_vo     = register_SIS_diag_field('ice_model', 'VO', diag%axesCv1, Time, &
               'surface current - y component', 'm/s', missing_value=missing)
  else
    OSS%id_uo     = register_SIS_diag_field('ice_model', 'UO', diag%axesB1, Time, &
               'surface current - x component', 'm/s', missing_value=missing)
    OSS%id_vo     = register_SIS_diag_field('ice_model', 'VO', diag%axesB1, Time, &
               'surface current - y component', 'm/s', missing_value=missing)
  endif
  FIA%id_sw_vis   = register_SIS_diag_field('ice_model','SW_VIS' ,diag%axesT1, Time, &
               'visible short wave heat flux', 'W/m^2', missing_value=missing)
  FIA%id_sw_dir   = register_SIS_diag_field('ice_model','SW_DIR' ,diag%axesT1, Time, &
               'direct short wave heat flux', 'W/m^2', missing_value=missing)
  FIA%id_sw_dif   = register_SIS_diag_field('ice_model','SW_DIF' ,diag%axesT1, Time, &
               'diffuse short wave heat flux', 'W/m^2', missing_value=missing)
  FIA%id_sw_vis_dir = register_SIS_diag_field('ice_model','SW_VIS_DIR' ,diag%axesT1, Time, &
               'visible direct short wave heat flux', 'W/m^2', missing_value=missing)
  FIA%id_sw_vis_dif = register_SIS_diag_field('ice_model','SW_VIS_DIF' ,diag%axesT1, Time, &
               'visible diffuse short wave heat flux', 'W/m^2', missing_value=missing)
  FIA%id_sw_nir_dir = register_SIS_diag_field('ice_model','SW_NIR_DIR' ,diag%axesT1, Time, &
               'near IR direct short wave heat flux', 'W/m^2', missing_value=missing)
  FIA%id_sw_nir_dif = register_SIS_diag_field('ice_model','SW_NIR_DIF' ,diag%axesT1, Time, &
               'near IR diffuse short wave heat flux', 'W/m^2', missing_value=missing)

  !
  ! diagnostics for quantities produced outside the ice model
  !
  IST%id_swdn  = register_SIS_diag_field('ice_model','SWDN' ,diag%axesT1, Time, &
             'downward shortwave flux', 'W/m^2', missing_value=missing)
  IST%id_lwdn  = register_SIS_diag_field('ice_model','LWDN' ,diag%axesT1, Time, &
             'downward longwave flux', 'W/m^2', missing_value=missing)
!### THIS DIAGNOSTIC IS MISSING.
! IST%id_ta    = register_SIS_diag_field('ice_model', 'TA', diag%axesT1, Time, &
!            'surface air temperature', 'C', missing_value=missing)
  IST%id_slp   = register_SIS_diag_field('ice_model', 'SLP', diag%axesT1, Time, &
             'sea level pressure', 'Pa', missing_value=missing)
  OSS%id_sst   = register_SIS_diag_field('ice_model', 'SST', diag%axesT1, Time, &
             'sea surface temperature', 'deg-C', missing_value=missing)
  OSS%id_sss   = register_SIS_diag_field('ice_model', 'SSS', diag%axesT1, Time, &
             'sea surface salinity', 'psu', missing_value=missing)
  OSS%id_ssh   = register_SIS_diag_field('ice_model', 'SSH', diag%axesT1, Time, &
             'sea surface height', 'm', missing_value=missing)
!### THIS DIAGNOSTIC IS MISSING.
!  IST%id_obi   = register_SIS_diag_field('ice_model', 'OBI', diag%axesT1, Time, &
!       'ice observed', '0 or 1', missing_value=missing)


  IST%id_rdgr    = register_SIS_diag_field('ice_model','RDG_RATE' ,diag%axesT1, Time, &
               'ice ridging rate', '1/sec', missing_value=missing)
!### THESE DIAGNOSTICS DO NOT EXIST YET.
!  IST%id_rdgf    = register_SIS_diag_field('ice_model','RDG_FRAC' ,diag%axesT1, Time, &
!               'ridged ice fraction', '0-1', missing_value=missing)
!  IST%id_rdgo    = register_SIS_diag_field('ice_model','RDG_OPEN' ,diag%axesT1, Time, &
!               'opening due to ridging', '1/s', missing_value=missing)
!  IST%id_rdgv    = register_SIS_diag_field('ice_model','RDG_VOSH' ,diag%axesT1, Time, &
!               'volume shifted from level to ridged ice', 'm^3/s', missing_value=missing)

  if (id_sin_rot>0) call post_data(id_sin_rot, G%sin_rot, diag, is_static=.true.)
  if (id_cos_rot>0) call post_data(id_cos_rot, G%cos_rot, diag, is_static=.true.)
  if (id_geo_lon>0) call post_data(id_geo_lon, G%geoLonT, diag, is_static=.true.)
  if (id_geo_lat>0) call post_data(id_geo_lat, G%geoLatT, diag, is_static=.true.)
  if (id_cell_area>0) call post_data(id_cell_area, &
            Ice%area / (16.0*atan(1.0)*G%Rad_Earth**2), diag, is_static=.true.)


end subroutine ice_diagnostics_init

subroutine safe_alloc_ids_1d(ids, nids)
  integer, allocatable :: ids(:)
  integer, intent(in)  :: nids

  if (.not.ALLOCATED(ids)) then
    allocate(ids(nids)) ; ids(:) = -1
  endif
end subroutine safe_alloc_ids_1d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_stock_pe - returns stocks of heat, water, etc. for conservation checks   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_stock_pe(Ice, index, value)

  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT

  type(ice_data_type) :: Ice
  integer, intent(in) :: index
  real, intent(out)   :: value
  type(ice_state_type), pointer :: IST => NULL()
  type(ice_grid_type),  pointer :: IG => NULL()

  integer :: i, j, k, m, isc, iec, jsc, jec, ncat
  real :: icebergs_value
  real :: LI
  real :: part_wt, I_NkIce, kg_H, kg_H_Nk

  value = 0.0
  if(.not.Ice%pe) return

  IST => Ice%Ice_state
  IG => Ice%IG

  isc = Ice%G%isc ; iec = Ice%G%iec ; jsc = Ice%G%jsc ; jec = Ice%G%jec
  ncat = IG%CatIce ; I_NkIce = 1.0 / IG%NkIce
  kg_H = IG%H_to_kg_m2 ; kg_H_Nk = IG%H_to_kg_m2 / IG%NkIce
  call get_SIS2_thermo_coefs(Ice%Ice_State%ITV, Latent_fusion=LI)

  select case (index)

    case (ISTOCK_WATER)
      value = 0.0
      do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + kg_H * (IST%mH_ice(i,j,k) + IST%mH_snow(i,j,k)) * &
               IST%part_size(i,j,k) * (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j))
      enddo ; enddo ; enddo

    case (ISTOCK_HEAT)
      value = 0.0
      if (IST%slab_ice) then
        do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
          if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
              value = value - (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j)) * IST%part_size(i,j,k) * &
                              (kg_H * IST%mH_ice(i,j,1)) * LI
          endif
        enddo ; enddo ; enddo
      else
        do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
          part_wt = (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j)) * IST%part_size(i,j,k)
          if (part_wt*IST%mH_ice(i,j,k) > 0.0) then
            value = value - (part_wt * (kg_H * IST%mH_snow(i,j,k))) * &
                Energy_melt_enthS(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
            do m=1,IG%NkIce
              value = value - (part_wt * (kg_H_Nk * IST%mH_ice(i,j,k))) * &
                  Energy_melt_enthS(IST%enth_ice(i,j,k,m), IST%sal_ice(i,j,k,m), IST%ITV)
            enddo
          endif
        enddo ; enddo ; enddo
      endif

    case (ISTOCK_SALT)
      !There is no salt in the snow.
      value = 0.0
      do m=1,IG%NkIce ; do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + (IST%part_size(i,j,k) * (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j))) * &
            (0.001*(kg_H_Nk*IST%mH_ice(i,j,k))) * IST%sal_ice(i,j,k,m)
      enddo ; enddo ; enddo ; enddo

    case default

      value = 0.0

  end select

  if (IST%do_icebergs) then
    call icebergs_stock_pe(Ice%icebergs, index, icebergs_value)
    value = value + icebergs_value
  endif

end subroutine ice_stock_pe

subroutine ice_data_type_chksum(id, timestep, Ice)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

  character(len=*), intent(in) :: id
  integer         , intent(in) :: timestep
  type(ice_data_type), intent(in) :: Ice
  integer ::   n, m, outunit

  outunit = stdout()
  write(outunit,*) "BEGIN CHECKSUM(ice_data_type):: ", id, timestep
  write(outunit,100) 'ice_data_type%part_size          ',mpp_chksum(Ice%part_size          )
  write(outunit,100) 'ice_data_type%albedo             ',mpp_chksum(Ice%albedo             )
  write(outunit,100) 'ice_data_type%albedo_vis_dir     ',mpp_chksum(Ice%albedo_vis_dir     )
  write(outunit,100) 'ice_data_type%albedo_nir_dir     ',mpp_chksum(Ice%albedo_nir_dir     )
  write(outunit,100) 'ice_data_type%albedo_vis_dif     ',mpp_chksum(Ice%albedo_vis_dif     )
  write(outunit,100) 'ice_data_type%albedo_nir_dif     ',mpp_chksum(Ice%albedo_nir_dif     )
  write(outunit,100) 'ice_data_type%rough_mom          ',mpp_chksum(Ice%rough_mom          )
  write(outunit,100) 'ice_data_type%rough_heat         ',mpp_chksum(Ice%rough_heat         )
  write(outunit,100) 'ice_data_type%rough_moist        ',mpp_chksum(Ice%rough_moist        )

  write(outunit,100) 'ice_data_type%t_surf             ',mpp_chksum(Ice%t_surf             )
  write(outunit,100) 'ice_data_type%u_surf             ',mpp_chksum(Ice%u_surf             )
  write(outunit,100) 'ice_data_type%v_surf             ',mpp_chksum(Ice%v_surf             )
  write(outunit,100) 'ice_data_type%s_surf             ',mpp_chksum(Ice%s_surf             )
  write(outunit,100) 'ice_data_type%flux_u             ',mpp_chksum(Ice%flux_u             )
  write(outunit,100) 'ice_data_type%flux_v             ',mpp_chksum(Ice%flux_v             )
  write(outunit,100) 'ice_data_type%flux_t             ',mpp_chksum(Ice%flux_t             )
  write(outunit,100) 'ice_data_type%flux_q             ',mpp_chksum(Ice%flux_q             )
  write(outunit,100) 'ice_data_type%flux_lw            ',mpp_chksum(Ice%flux_lw            )
  write(outunit,100) 'ice_data_type%flux_sw_vis_dir    ',mpp_chksum(Ice%flux_sw_vis_dir    )
  write(outunit,100) 'ice_data_type%flux_sw_vis_dif    ',mpp_chksum(Ice%flux_sw_vis_dif    )
  write(outunit,100) 'ice_data_type%flux_sw_nir_dir    ',mpp_chksum(Ice%flux_sw_nir_dir    )
  write(outunit,100) 'ice_data_type%flux_sw_nir_dif    ',mpp_chksum(Ice%flux_sw_nir_dif    )
  write(outunit,100) 'ice_data_type%flux_lh            ',mpp_chksum(Ice%flux_lh            )
  write(outunit,100) 'ice_data_type%lprec              ',mpp_chksum(Ice%lprec              )
  write(outunit,100) 'ice_data_type%fprec              ',mpp_chksum(Ice%fprec              )
  write(outunit,100) 'ice_data_type%p_surf             ',mpp_chksum(Ice%p_surf             )
  write(outunit,100) 'ice_data_type%runoff             ',mpp_chksum(Ice%runoff             )
  write(outunit,100) 'ice_data_type%calving            ',mpp_chksum(Ice%calving            )
  write(outunit,100) 'ice_data_type%flux_salt          ',mpp_chksum(Ice%flux_salt          )

  do n=1,Ice%ocean_fields%num_bcs ; do m=1,Ice%ocean_fields%bc(n)%num_fields
    write(outunit,101) 'ice%', trim(Ice%ocean_fields%bc(n)%name), &
                       trim(Ice%ocean_fields%bc(n)%field(m)%name), &
                       mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
  enddo ; enddo

100 FORMAT("   CHECKSUM::",A32," = ",Z20)
101 FORMAT("   CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ice_data_type_chksum


subroutine ocn_ice_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

  character(len=*), intent(in) :: id
  integer         , intent(in) :: timestep
  type(ocean_ice_boundary_type), intent(in) :: bnd_type
  integer ::   n, m, outunit

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(ocean_ice_boundary_type):: ', id, timestep
  write(outunit,100) 'ocn_ice_bnd_type%u        ',mpp_chksum(bnd_type%u        )
  write(outunit,100) 'ocn_ice_bnd_type%v        ',mpp_chksum(bnd_type%v        )
  write(outunit,100) 'ocn_ice_bnd_type%t        ',mpp_chksum(bnd_type%t        )
  write(outunit,100) 'ocn_ice_bnd_type%s        ',mpp_chksum(bnd_type%s        )
  write(outunit,100) 'ocn_ice_bnd_type%frazil   ',mpp_chksum(bnd_type%frazil   )
  write(outunit,100) 'ocn_ice_bnd_type%sea_level',mpp_chksum(bnd_type%sea_level)
  !    write(outunit,100) 'ocn_ice_bnd_type%data     ',mpp_chksum(bnd_type%data     )
  100 FORMAT("CHECKSUM::",A32," = ",Z20)

  do n = 1, bnd_type%fields%num_bcs  !{
    do m = 1, bnd_type%fields%bc(n)%num_fields  !{
        write(outunit,101) 'oibt%',trim(bnd_type%fields%bc(n)%name), &
             trim(bnd_type%fields%bc(n)%field(m)%name), &
             mpp_chksum(bnd_type%fields%bc(n)%field(m)%values)
    enddo  !} m
  enddo  !} n
  101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ocn_ice_bnd_type_chksum

subroutine atm_ice_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

  character(len=*), intent(in) :: id
  integer         , intent(in) :: timestep
  type(atmos_ice_boundary_type), intent(in) :: bnd_type
  integer ::   n, outunit

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(atmos_ice_boundary_type):: ', id, timestep
  write(outunit,100) 'atm_ice_bnd_type%u_flux          ',mpp_chksum(bnd_type%u_flux)
  write(outunit,100) 'atm_ice_bnd_type%v_flux          ',mpp_chksum(bnd_type%v_flux)
  write(outunit,100) 'atm_ice_bnd_type%u_star          ',mpp_chksum(bnd_type%u_star)
  write(outunit,100) 'atm_ice_bnd_type%t_flux          ',mpp_chksum(bnd_type%t_flux)
  write(outunit,100) 'atm_ice_bnd_type%q_flux          ',mpp_chksum(bnd_type%q_flux)
  write(outunit,100) 'atm_ice_bnd_type%lw_flux         ',mpp_chksum(bnd_type%lw_flux)
  write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dir ',mpp_chksum(bnd_type%sw_flux_vis_dir)
  write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dif ',mpp_chksum(bnd_type%sw_flux_vis_dif)
  write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dir ',mpp_chksum(bnd_type%sw_flux_nir_dir)
  write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dif ',mpp_chksum(bnd_type%sw_flux_nir_dif)
  write(outunit,100) 'atm_ice_bnd_type%lprec           ',mpp_chksum(bnd_type%lprec)
  write(outunit,100) 'atm_ice_bnd_type%fprec           ',mpp_chksum(bnd_type%fprec)
  write(outunit,100) 'atm_ice_bnd_type%dhdt            ',mpp_chksum(bnd_type%dhdt)
  write(outunit,100) 'atm_ice_bnd_type%dedt            ',mpp_chksum(bnd_type%dedt)
  write(outunit,100) 'atm_ice_bnd_type%drdt            ',mpp_chksum(bnd_type%drdt)
  write(outunit,100) 'atm_ice_bnd_type%coszen          ',mpp_chksum(bnd_type%coszen)
  write(outunit,100) 'atm_ice_bnd_type%p               ',mpp_chksum(bnd_type%p)
!    write(outunit,100) 'atm_ice_bnd_type%data            ',mpp_chksum(bnd_type%data)
100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine atm_ice_bnd_type_chksum

subroutine lnd_ice_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

  character(len=*), intent(in) :: id
  integer         , intent(in) :: timestep
  type(land_ice_boundary_type), intent(in) :: bnd_type
  integer ::   n, outunit

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(land_ice_boundary_type):: ', id, timestep
  write(outunit,100) 'lnd_ice_bnd_type%runoff  ',mpp_chksum(bnd_type%runoff)
  write(outunit,100) 'lnd_ice_bnd_type%calving ',mpp_chksum(bnd_type%calving)
  !    write(outunit,100) 'lnd_ice_bnd_type%data    ',mpp_chksum(bnd_type%data)
  100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine lnd_ice_bnd_type_chksum

end module ice_type_mod
