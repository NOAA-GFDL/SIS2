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
use constants_mod,    only: T_0degC=>Tfreeze

use ice_grid_mod,     only: sea_ice_grid_type, cell_area

use ice_thm_mod,      only: e_to_melt
use ice_dyn_bgrid,    only: ice_B_dyn_CS
use ice_dyn_cgrid,    only: ice_C_dyn_CS
use ice_transport_mod, only: ice_transport_CS
use SIS2_ice_thm, only : ice_thermo_type, SIS2_ice_thm_CS, enth_from_TS
use constants_mod,    only: radius, pi, LI => hlf ! latent heat of fusion - 334e3 J/(kg-ice)
use ice_bergs, only: icebergs, icebergs_stock_pe, icebergs_save_restart

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg, is_root_pe
use MOM_file_parser, only : param_file_type
use SIS_diag_mediator, only : SIS_diag_ctrl, post_data=>post_SIS_data
use SIS_diag_mediator, only : register_SIS_diag_field, register_static_field
use SIS_error_checking, only : chksum, Bchksum, hchksum, uchksum, vchksum
use SIS_error_checking, only : check_redundant_B, check_redundant_C
use SIS_get_input, only : archaic_nml_check
use SIS_sum_output_type, only : SIS_sum_out_CS
use SIS_tracer_registry, only : SIS_tracer_registry_type

implicit none ; private

#include <SIS2_memory.h>

public :: ice_data_type, ice_state_type
public :: ice_model_restart, dealloc_ice_arrays, dealloc_IST_arrays
public :: ice_data_type_register_restarts, ice_state_register_restarts
public :: ice_diagnostics_init, ice_stock_pe, Ice_restart, check_ice_model_nml
public :: ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum, ice_data_type_chksum
public :: IST_chksum, Ice_public_type_chksum, Ice_public_type_bounds_check, IST_bounds_check

public  :: earth_area

  real, parameter :: earth_area = 4*PI*RADIUS*RADIUS !5.10064471909788E+14 m^2
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
  integer :: avg_count
!   logical                            :: pe

!   logical, pointer, dimension(:,:)   :: mask                =>NULL() ! where ice can be
  real, pointer, dimension(:,:,:) :: &
    part_size =>NULL()     ! The fractional coverage of a grid cell by each ice
                           ! thickness category, nondim, 0 to 1.  Category 0 is
                           ! open ocean.  The sum of part_size is 1.

  ! The following are the 6 variables that constitute the sea-ice state.
  real, pointer, dimension(:,:) :: &
    u_ice =>NULL(), & ! The pseudo-zonal and pseudo-meridional ice velocities
    v_ice =>NULL(), & ! along the model's grid directions on a B-grid, in m s-1.
                      ! All thickness categories are assumed to have the same
                      ! velocity.
    u_ice_C =>NULL(), & ! The pseudo-zonal and pseudo-meridional ice velocities
    v_ice_C =>NULL()  ! along the model's grid directions on a C-grid, in m s-1.
                      ! All thickness categories are assumed to have the same
                      ! velocity.
  real, pointer, dimension(:,:,:) :: &
    h_snow =>NULL(), &  ! The thickness of the snow in each category, in m.
    h_ice =>NULL(), &   ! The thickness of the ice in each category, in m.
    t_snow =>NULL()     ! The temperture of the snow in each category, in degC.
  real, pointer, dimension(:,:,:,:) :: &
    t_ice =>NULL(), &   ! The temperature of the sea ice in each category and
                        ! fractional thickness layer, in degC.
    sal_ice =>NULL(), & ! The salinity of the sea ice in each category and
                        ! fractional thickness layer, in g/kg.
    enth_ice =>NULL(), & ! The enthalpy of the sea ice in each category and
                        ! fractional thickness layer, in enth_unit (J or rescaled).
    enth_snow =>NULL()  ! The enthalpy of the snow in each category, in enth_unit.

  real,    pointer, dimension(:,:) :: &
    s_surf  =>NULL(), &    ! The ocean's surface salinity in g/kg.
    t_ocn   =>NULL(), &    ! The ocean's bulk surface temperature in degC.
    u_ocn   =>NULL(), &    ! The ocean's zonal velocity on B-grid points in m s-1.
    v_ocn   =>NULL(), &    ! The ocean's meridional velocity on B-grid points in m s-1.
    u_ocn_C =>NULL(), &    ! The ocean's zonal and meridional velocity on C-grid
    v_ocn_C =>NULL(), &    ! points, both in m s-1.
    sea_lev =>NULL()       ! The equivalent sea-level, after any non-levitating
                           ! ice has been converted to sea-water, as determined
                           ! by the ocean, in m.  Sea-ice only contributes by
                           ! applying pressure to the ocean that is then
                           ! (partially) converted back to its equivalent by the
                           ! ocean. 
  real, pointer, dimension(:,:,:) :: &
    enth_prev, heat_in
  real,    pointer, dimension(:,:,:) :: &
    ! The 3rd dimension in each of the following is ice thickness category.
    t_surf              =>NULL(), & ! The surface temperature, in Kelvin.
    flux_u_top          =>NULL(), & ! The downward? flux of zonal and meridional
    flux_v_top          =>NULL(), & ! momentum on an A-grid in ???.
    flux_u_top_Bgrid    =>NULL(), & ! The downward? flux of zonal and meridional
    flux_v_top_Bgrid    =>NULL(), & ! momentum on a B-grid velocity point in ???.
    flux_u_top_Cu       =>NULL(), & ! The downward? flux of zonal and meridional
    flux_v_top_Cv       =>NULL(), & ! momentum on a B-grid velocity point in ???.
    flux_t_top          =>NULL(), & ! The upward sensible heat flux at the ice top
                                    ! in W m-2.
    flux_q_top          =>NULL(), & ! The upward evaporative moisture flux at
                                    ! top of the ice, in kg m-2 s-1.
    flux_lw_top         =>NULL(), & ! The downward flux of longwave radiation at
                                    ! the top of the ice, in W m-2.
    flux_sw_vis_dir_top =>NULL(), & ! The downward diffuse flux of direct (dir)
    flux_sw_vis_dif_top =>NULL(), & ! and diffuse (dif) shortwave radiation in
    flux_sw_nir_dir_top =>NULL(), & ! the visible (vis) and near-infrared (nir)
    flux_sw_nir_dif_top =>NULL(), & ! bands at the top of the ice, in W m-2.
    flux_lh_top         =>NULL(), & ! The upward flux of latent heat at the top
                                    ! of the ice, in W m-2.
    lprec_top           =>NULL(), & ! The downward flux of liquid precipitation
                                    ! at the top of the ice, in kg m-2 s-1.
    fprec_top           =>NULL()    ! The downward flux of frozen precipitation
                                    ! at the top of the ice, in kg m-2 s-1.

  real, pointer, dimension(:,:)   :: &
    ! These terms diagnose the enthalpy change associated with the addition or
    ! removal of water mass (liquid or frozen) from the ice model are required
    ! to close the enthalpy budget. Ice enthalpy is generally negative, so terms
    ! that add mass to the ice are generally negative.
    Enth_Mass_in_atm  =>NULL(), & ! The enthalpy introduced to the ice by water
                                  ! fluxes from the atmosphere, in J m-2.
    Enth_Mass_out_atm =>NULL(), & ! Negative of the enthalpy extracted from the
                                  ! ice by water fluxes to the atmosphere, in J m-2.
    Enth_Mass_in_ocn  =>NULL(), & ! The enthalpy introduced to the ice by water
                                  ! fluxes from the ocean, in J m-2.
    Enth_Mass_out_ocn =>NULL(), & ! Negative of the enthalpy extracted from the
                                  ! ice by water fluxes to the ocean, in J m-2.

    flux_t_ocn_top => NULL(), &   ! The upward sensible heat flux from the ocean
                                  ! to the ice or atmosphere, in W m-2.
    flux_q_ocn_top => NULL(), &   ! The upward evaporative moisture flux at
                                  ! the ocean surface, in kg m-2 s-1.
    flux_lw_ocn_top =>NULL(), &   ! The downward flux of longwave radiation at
                                  ! the ocean surface, in W m-2.
    flux_sw_vis_dir_ocn =>NULL(), & ! The downward diffuse flux of direct (dir)
    flux_sw_vis_dif_ocn =>NULL(), & ! and diffuse (dif) shortwave radiation in
    flux_sw_nir_dir_ocn =>NULL(), & ! the visible (vis) and near-infrared (nir)
    flux_sw_nir_dif_ocn =>NULL(), & ! bands at the ocean surface, in W m-2.
    flux_lh_ocn_top =>NULL(), &   ! The upward flux of latent heat at the
                                  ! ocean surface, in W m-2.
    lprec_ocn_top => NULL(), &    ! The downward flux of liquid precipitation at
                                  ! the ocena surface, in kg m-2 s-1.
    fprec_ocn_top => NULL(), &    ! The downward flux of frozen precipitation at
                                  ! the ocena surface, in kg m-2 s-1.
  !  ### ADD BETTER COMMENTS, WITH UNITS.
    lwdn         =>NULL(), &      ! Accumulated diagnostics of downward long-
    swdn         =>NULL()         ! and short-wave radiation <WHERE?> in <UNITS?>.

  real, pointer, dimension(:,:,:) :: sw_abs_sfc   =>NULL() ! frac abs sw abs @ surf.
  real, pointer, dimension(:,:,:) :: sw_abs_snow  =>NULL() ! frac abs sw abs in snow
  real, pointer, dimension(:,:,:,:) :: sw_abs_ice =>NULL() ! frac abs sw abs in ice layers
  real, pointer, dimension(:,:,:) :: sw_abs_ocn   =>NULL() ! frac abs sw abs in ocean
  real, pointer, dimension(:,:,:) :: sw_abs_int   =>NULL() ! frac abs sw abs in ice interior
  real, pointer, dimension(:,:)   :: coszen       =>NULL()
  real, pointer, dimension(:,:,:) :: tmelt        =>NULL()
  real, pointer, dimension(:,:,:) :: bmelt        =>NULL()

  real, pointer, dimension(:,:)   :: frazil       =>NULL()
  real, pointer, dimension(:,:)   :: frazil_input =>NULL()
  real, pointer, dimension(:,:)   :: bheat        =>NULL()
  real, pointer, dimension(:,:)   :: mi           =>NULL() ! The total ice+snow mass, in kg m-2.
  logical :: slab_ice  ! If true, do the old style GFDL slab ice.
  logical :: Cgrid_dyn ! If true use a C-grid discretization of the
                       ! sea-ice dynamics.
  logical :: SIS1_5L_thermo ! If true, the thermodynamic calculations inhereted
                       ! from the 5-layer version of SIS1. Otherwise, use the
                       ! newer SIS2 version.                  
  real :: Rho_ocean    ! The nominal density of sea water, in kg m-3.
  real :: Rho_ice      ! The nominal density of sea ice, in kg m-3.
  real :: Rho_snow     ! The nominal density of snow on sea ice, in kg m-3.
  logical :: do_icebergs    ! If true, use the Lagrangian iceberg code, which
                            ! modifies the calving field among other things.
  logical :: specified_ice  ! If true, the sea ice is specified and there is
                            ! no need for ice dynamics.
  logical :: column_check   ! If true, enable the heat check column by column.
  real    :: imb_tol        ! The tolerance for imbalances to be flagged by
                            ! column_check, nondim.
  logical :: bounds_check    ! If true, check for sensible values of thicknesses
                             ! temperatures, fluxes, etc.
  logical :: debug           ! If true, write verbose checksums for debugging purposes.
  type(time_type) :: ice_stats_interval ! The interval between writes of the
                             ! globally summed ice statistics and conservation checks.
  type(time_type) :: write_ice_stats_time ! The next time to write out the ice statistics.
  

  logical :: atmos_winds ! The wind stresses come directly from the atmosphere
                         ! model and have the wrong sign.
  real :: kmelt          ! A constant that is used in the calculation of the 
                         ! ocean/ice basal heat flux, in W m-2 K-1.
  real :: ice_bulk_salin ! The globally constant sea ice bulk salinity, in g/kg
                         ! that is used to calculate the ocean salt flux.
  real :: ice_rel_salin  ! The initial bulk salinity of sea-ice relative to the
                         ! salinity of the water from which it formed, nondim.
  logical :: do_ice_restore ! If true, restore the sea-ice toward climatology
                            ! by applying a restorative heat flux.
  real    :: ice_restore_timescale ! The time scale for restoring ice when
                            ! do_ice_restore is true, in days.
  logical :: do_ice_limit   ! Limit the sea ice thickness to max_ice_limit.
  real    :: max_ice_limit  ! The maximum sea ice thickness, in m, when
                            ! do_ice_limit is true.
  logical :: slp2ocean  ! If true, apply sea level pressure to ocean surface.
  logical :: verbose    ! A flag to control the printing of an ice-diagnostic 
                        ! message.  When true, this will slow the model down.
  logical :: add_diurnal_sw ! If true, apply a synthetic diurnal cycle to the shortwave radiation.
  logical :: do_sun_angle_for_alb ! If true, find the sun angle for calculating
                                  ! the ocean albedo in the frame of the ice model.

  integer :: n_calls = 0     ! The number of times update_ice_model_slow_down
                             ! has been called.
  integer :: n_fast = 0      ! The number of times update_ice_model_fast
                             ! has been called.
  logical :: do_init = .false. ! If true, there is still some initialization
                               ! that needs to be done.
  logical :: first_time = .true. ! If true, this is the first call to 
                               ! update_ice_model_slow_up 

!   type(coupler_3d_bc_type)   :: ocean_fields       ! array of fields used for additional tracers
!   type(coupler_2d_bc_type)   :: ocean_fluxes       ! array of fluxes used for additional tracers
!   type(coupler_3d_bc_type)   :: ocean_fluxes_top   ! array of fluxes for averaging

  integer, dimension(:), allocatable :: id_t, id_sw_abs_ice, id_sal
  integer :: id_cn=-1, id_hi=-1, id_hs=-1, id_tsn=-1
  integer :: id_ts=-1, id_t_iceav=-1, id_s_iceav=-1, id_hio=-1, id_mi=-1, id_sh=-1
  integer :: id_lh=-1, id_sw=-1, id_lw=-1, id_snofl=-1, id_rain=-1, id_runoff=-1
  integer :: id_calving=-1, id_runoff_hflx=-1, id_calving_hflx=-1, id_evap=-1
  integer :: id_saltf=-1, id_tmelt=-1, id_bmelt=-1, id_bheat=-1, id_e2m=-1
  integer :: id_frazil=-1, id_alb=-1, id_xprt=-1, id_lsrc=-1, id_lsnk=-1, id_bsnk=-1
  integer :: id_strna=-1, id_fax=-1, id_fay=-1, id_swdn=-1, id_lwdn=-1, id_sn2ic=-1
  integer :: id_slp=-1, id_ext=-1, id_sst=-1, id_sss=-1, id_ssh=-1, id_uo=-1, id_vo=-1
  integer :: id_ta=-1, id_obi=-1, id_qfres=-1, id_qflim=-1, id_ix_trans=-1
  integer :: id_iy_trans=-1, id_sw_vis=-1, id_sw_dir=-1, id_sw_dif=-1
  integer :: id_sw_vis_dir=-1, id_sw_vis_dif=-1, id_sw_nir_dir=-1, id_sw_nir_dif=-1
  integer :: id_mib=-1, id_coszen=-1
  integer :: id_alb_vis_dir=-1, id_alb_vis_dif=-1, id_alb_nir_dir=-1, id_alb_nir_dif=-1
  integer :: id_abs_int=-1, id_sw_abs_sfc=-1, id_sw_abs_snow=-1
  integer :: id_sw_pen=-1, id_sw_abs_ocn=-1

  type(SIS_tracer_registry_type), pointer :: TrReg => NULL()

  type(ice_B_dyn_CS), pointer     :: ice_B_dyn_CSp => NULL()
  type(ice_C_dyn_CS), pointer     :: ice_C_dyn_CSp => NULL()
  type(ice_transport_CS), pointer :: ice_transport_CSp => NULL()
  type(ice_thermo_type), pointer  :: ITV => NULL()
  type(SIS2_ice_thm_CS), pointer  :: ice_thm_CSp => NULL()
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()
  type(SIS_diag_ctrl)             :: diag ! A structure that regulates diagnostis.
!   type(icebergs), pointer     :: icebergs => NULL()
end type ice_state_type

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! This structure contains the ice model data (some used by calling routines);  !
! the third index is partition (1 is open water; 2 is ice cover)               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type ice_data_type !  ice_public_type
  type(domain2D)                     :: Domain
  type(time_type)                    :: Time
  logical                            :: pe
  integer, pointer, dimension(:)     :: pelist              =>NULL() ! Used for flux-exchange.
     logical, pointer, dimension(:,:)   :: mask                =>NULL() ! where ice can be
  logical, pointer, dimension(:,:,:) :: ice_mask            =>NULL() ! where ice actually is (Used for k-size only?)

  ! These fields are used to provide information about the ice surface to the
  ! atmosphere, and contain separate values for each ice thickness category.
  real, pointer, dimension(:,:,:) :: &
    !### ADD COMMENTS DESCRIBING EACH FIELD.
    part_size => NULL(), &
    albedo    => NULL(), &
    albedo_vis_dir => NULL(), &
    albedo_nir_dir => NULL(), &
    albedo_vis_dif => NULL(), &
    albedo_nir_dif => NULL(), &
    rough_mom   => NULL(), &
    rough_heat  => NULL(), &
    rough_moist => NULL(), &
    t_surf      => NULL(), &
    u_surf      => NULL(), &
    v_surf      => NULL()
  real, pointer, dimension(:,:)   :: s_surf         =>NULL()

  ! These arrays will be used to set the forcing for the ocean.
  real, pointer, dimension(:,:) :: &
    flux_u => NULL(), &   ! The flux of x-momentum into the ocean, in Pa.
    flux_v => NULL(), &   ! The flux of y-momentum into the ocean, in Pa.
    flux_t => NULL(), &   ! The flux of sensible heat out of the ocean, in W m-2.
    flux_q => NULL(), &   ! The evaporative moisture flux out of the ocean, in kg m-2 s-1.
    flux_lw => NULL(), &  ! The sensible heat flux out of the ocena, in W m-2.
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

  real, pointer, dimension(:,:) :: area => NULL()
  real, pointer, dimension(:,:) :: mi   => NULL() ! The total ice+snow mass, in kg m-2.
             ! mi is needed for the wave model. It is introduced here,
             ! because flux_ice_to_ocean cannot handle 3D fields. This may be
			       ! removed, if the information on ice thickness can be derived from 
			       ! eventually from h_ice outside the ice module.
  integer, dimension(3)    :: axes
  type(coupler_3d_bc_type) :: ocean_fields       ! array of fields used for additional tracers
  type(coupler_2d_bc_type) :: ocean_fluxes       ! array of fluxes used for additional tracers
  type(coupler_3d_bc_type) :: ocean_fluxes_top   ! array of fluxes for averaging
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
  type(sea_ice_grid_type), pointer :: G ! A structure containing metrics and grid info.
  type(ice_state_type), pointer :: Ice_state => NULL() ! A structure containing the internal
                               ! representation of the ice state.
  type(restart_file_type), pointer :: Ice_restart
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
    frazil => NULL(), &  ! The frazil heat rejected by the ocean, in J.
    sea_level => NULL()  ! The sea level after adjustment for any surface
                         ! pressure that the ocean allows to be expressed, in m.
  real, dimension(:,:,:), pointer :: data      =>NULL() ! collective field for "named" fields above
  integer                         :: stagger = BGRID_NE
  integer                         :: xtype              ! REGRID, REDIST or DIRECT used by coupler
  type(coupler_2d_bc_type)        :: fields     ! array of fields used for additional tracers
end type 

type :: atmos_ice_boundary_type 
!### These arrays need to be described in comments with units and directionality.
  real, dimension(:,:,:), pointer :: u_flux  =>NULL()
  real, dimension(:,:,:), pointer :: v_flux  =>NULL()
  real, dimension(:,:,:), pointer :: u_star  =>NULL()
  real, dimension(:,:,:), pointer :: t_flux  =>NULL()
  real, dimension(:,:,:), pointer :: q_flux  =>NULL()
  real, dimension(:,:,:), pointer :: lw_flux =>NULL()
  real, dimension(:,:,:), pointer :: sw_flux_vis_dir =>NULL()
  real, dimension(:,:,:), pointer :: sw_flux_vis_dif =>NULL()
  real, dimension(:,:,:), pointer :: sw_flux_nir_dir =>NULL()
  real, dimension(:,:,:), pointer :: sw_flux_nir_dif =>NULL()
  real, dimension(:,:,:), pointer :: lprec   =>NULL()
  real, dimension(:,:,:), pointer :: fprec   =>NULL()
  real, dimension(:,:,:), pointer :: dhdt    =>NULL()
  real, dimension(:,:,:), pointer :: dedt    =>NULL()
  real, dimension(:,:,:), pointer :: drdt    =>NULL()
  real, dimension(:,:,:), pointer :: coszen  =>NULL()
  real, dimension(:,:,:), pointer :: p       =>NULL()
  real, dimension(:,:,:), pointer :: data    =>NULL()
  integer                         :: xtype
  type(coupler_3d_bc_type)        :: fluxes     ! array of fluxes used for additional tracers
end type

type :: land_ice_boundary_type
  real, dimension(:,:),   pointer :: runoff  =>NULL()
  real, dimension(:,:),   pointer :: calving =>NULL()
  real, dimension(:,:),   pointer :: runoff_hflx  =>NULL()
  real, dimension(:,:),   pointer :: calving_hflx =>NULL()
  real, dimension(:,:,:), pointer :: data    =>NULL() ! collective field for "named" fields above
  integer                         :: xtype            ! REGRID, REDIST or DIRECT used by coupler
end type

type(restart_file_type), pointer, save :: Ice_restart

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

  allocate(Ice%mask(isc:iec, jsc:jec)) ; Ice%mask(:,:) = .false. !derived
  allocate(Ice%ice_mask(isc:iec, jsc:jec, km)) ; Ice%ice_mask(:,:,:) = .false. !NI
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
  idr = register_restart_field(Ice_restart, restart_file, 'albedo',    Ice%albedo,    domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'albedo_vis_dir', Ice%albedo_vis_dir, &
                                      domain=domain, mandatory=.false.)
  idr = register_restart_field(Ice_restart, restart_file, 'albedo_nir_dir', Ice%albedo_nir_dir, &
                                      domain=domain, mandatory=.false.)
  idr = register_restart_field(Ice_restart, restart_file, 'albedo_vis_dif', Ice%albedo_vis_dif, &
                                      domain=domain, mandatory=.false.)
  idr = register_restart_field(Ice_restart, restart_file, 'albedo_nir_dif', Ice%albedo_nir_dif, &
                                      domain=domain, mandatory=.false.)
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
! ice_state_register_restarts - allocate the arrays in the ice_state_type      !
!     and register any variables in the ice data type that need to be included !
!     in the restart files.                                                    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_state_register_restarts(G, param_file, IST, Ice_restart, restart_file)
  type(sea_ice_grid_type), intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(ice_state_type),    intent(inout) :: IST
  type(restart_file_type), intent(inout) :: Ice_restart
  character(len=*),        intent(in)    :: restart_file
  
  type(domain2d), pointer :: domain
  integer :: CatIce, idr, n
  character(len=8) :: nstr

  CatIce = G%CatIce
  allocate(IST%t_surf(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%t_surf(:,:,:) = 0.0 !X
  allocate(IST%s_surf(SZI_(G), SZJ_(G))) ; IST%s_surf(:,:) = 0.0 !NI X
  allocate(IST%t_ocn(SZI_(G), SZJ_(G))) ; IST%t_ocn(:,:) = 0.0   !NI X
  allocate(IST%sea_lev(SZI_(G), SZJ_(G))) ; IST%sea_lev(:,:) = 0.0 !NR 
  allocate(IST%part_size(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%part_size(:,:,:) = 0.0
  allocate(IST%coszen(SZI_(G), SZJ_(G))) ; IST%coszen(:,:) = 0.0 !NR X

  allocate(IST%flux_u_top(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%flux_u_top(:,:,:) = 0.0 !NR
  allocate(IST%flux_v_top(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%flux_v_top(:,:,:) = 0.0 !NR 
  allocate(IST%flux_t_top(SZI_(G), SZJ_(G), 0:CatIce)) ;  IST%flux_t_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_q_top(SZI_(G), SZJ_(G), 0:CatIce)) ;  IST%flux_q_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_sw_vis_dir_top(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%flux_sw_vis_dir_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_sw_vis_dif_top(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%flux_sw_vis_dif_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_sw_nir_dir_top(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%flux_sw_nir_dir_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_sw_nir_dif_top(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%flux_sw_nir_dif_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_lw_top(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%flux_lw_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_lh_top(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%flux_lh_top(:,:,:) = 0.0 !NI
  allocate(IST%lprec_top(SZI_(G), SZJ_(G), 0:CatIce)) ;  IST%lprec_top(:,:,:) = 0.0 !NI
  allocate(IST%fprec_top(SZI_(G), SZJ_(G), 0:CatIce)) ;  IST%fprec_top(:,:,:) = 0.0 !NI

  allocate(IST%Enth_Mass_in_atm(SZI_(G), SZJ_(G)))  ; IST%Enth_Mass_in_atm(:,:) = 0.0 !NR
  allocate(IST%Enth_Mass_out_atm(SZI_(G), SZJ_(G))) ; IST%Enth_Mass_out_atm(:,:) = 0.0 !NR
  allocate(IST%Enth_Mass_in_ocn(SZI_(G), SZJ_(G)))  ; IST%Enth_Mass_in_ocn(:,:) = 0.0 !NR
  allocate(IST%Enth_Mass_out_ocn(SZI_(G), SZJ_(G))) ; IST%Enth_Mass_out_ocn(:,:) = 0.0 !NR
  allocate(IST%flux_t_ocn_top(SZI_(G), SZJ_(G))) ;  IST%flux_t_ocn_top(:,:) = 0.0 !NI
  allocate(IST%flux_q_ocn_top(SZI_(G), SZJ_(G))) ;  IST%flux_q_ocn_top(:,:) = 0.0 !NI
  allocate(IST%flux_lw_ocn_top(SZI_(G), SZJ_(G))) ; IST%flux_lw_ocn_top(:,:) = 0.0 !NI
  allocate(IST%flux_lh_ocn_top(SZI_(G), SZJ_(G))) ; IST%flux_lh_ocn_top(:,:) = 0.0 !NI
  allocate(IST%flux_sw_vis_dir_ocn(SZI_(G), SZJ_(G))) ;  IST%flux_sw_vis_dir_ocn(:,:) = 0.0 !NI
  allocate(IST%flux_sw_vis_dif_ocn(SZI_(G), SZJ_(G))) ;  IST%flux_sw_vis_dif_ocn(:,:) = 0.0 !NI
  allocate(IST%flux_sw_nir_dir_ocn(SZI_(G), SZJ_(G))) ;  IST%flux_sw_nir_dir_ocn(:,:) = 0.0 !NI
  allocate(IST%flux_sw_nir_dif_ocn(SZI_(G), SZJ_(G))) ;  IST%flux_sw_nir_dif_ocn(:,:) = 0.0 !NI
  allocate(IST%lprec_ocn_top(SZI_(G), SZJ_(G))) ;  IST%lprec_ocn_top(:,:) = 0.0 !NI
  allocate(IST%fprec_ocn_top(SZI_(G), SZJ_(G))) ;  IST%fprec_ocn_top(:,:) = 0.0 !NI

  allocate(IST%lwdn(SZI_(G), SZJ_(G))) ; IST%lwdn(:,:) = 0.0 !NR
  allocate(IST%swdn(SZI_(G), SZJ_(G))) ; IST%swdn(:,:) = 0.0 !NR
  allocate(IST%frazil(SZI_(G), SZJ_(G))) ; IST%frazil(:,:) = 0.0 !NR
  allocate(IST%frazil_input(SZI_(G), SZJ_(G))) ; IST%frazil_input(:,:) = 0.0 !NR
  allocate(IST%bheat(SZI_(G), SZJ_(G))) ; IST%bheat(:,:) = 0.0 !NI
  allocate(IST%tmelt(SZI_(G), SZJ_(G), CatIce)) ; IST%tmelt(:,:,:) = 0.0 !NR
  allocate(IST%bmelt(SZI_(G), SZJ_(G), CatIce)) ; IST%bmelt(:,:,:) = 0.0 !NR

  allocate(IST%sw_abs_sfc(SZI_(G), SZJ_(G), CatIce)) ; IST%sw_abs_sfc(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_snow(SZI_(G), SZJ_(G), CatIce)) ; IST%sw_abs_snow(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_ice(SZI_(G), SZJ_(G), CatIce, G%NkIce)) ; IST%sw_abs_ice(:,:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_ocn(SZI_(G), SZJ_(G), CatIce)) ; IST%sw_abs_ocn(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_int(SZI_(G), SZJ_(G), CatIce)) ; IST%sw_abs_int(:,:,:) = 0.0 !NR

  allocate(IST%h_snow(SZI_(G), SZJ_(G), CatIce)) ; IST%h_snow(:,:,:) = 0.0
  allocate(IST%t_snow(SZI_(G), SZJ_(G), CatIce)) ; IST%t_snow(:,:,:) = 0.0
  allocate(IST%enth_snow(SZI_(G), SZJ_(G), CatIce, 1)) ; IST%enth_snow(:,:,:,:) = 0.0
  allocate(IST%h_ice(SZI_(G), SZJ_(G), CatIce)) ; IST%h_ice(:,:,:) = 0.0
  allocate(IST%t_ice(SZI_(G), SZJ_(G), CatIce, G%NkIce)) ; IST%t_ice(:,:,:,:) = 0.0
  allocate(IST%enth_ice(SZI_(G), SZJ_(G), CatIce, G%NkIce)) ; IST%enth_ice(:,:,:,:) = 0.0
  allocate(IST%sal_ice(SZI_(G), SZJ_(G), CatIce, G%NkIce)) ; IST%sal_ice(:,:,:,:) = 0.0

  allocate(IST%enth_prev(SZI_(G), SZJ_(G), CatIce)) ; IST%enth_prev(:,:,:) = 0.0
  allocate(IST%heat_in(SZI_(G), SZJ_(G), CatIce)) ; IST%heat_in(:,:,:) = 0.0


  if (IST%Cgrid_dyn) then
    allocate(IST%u_ice_C(SZIB_(G), SZJ_(G))) ; IST%u_ice_C(:,:) = 0.0
    allocate(IST%v_ice_C(SZI_(G), SZJB_(G))) ; IST%v_ice_C(:,:) = 0.0
    allocate(IST%u_ocn_C(SZIB_(G), SZJ_(G))) ; IST%u_ocn_C(:,:) = 0.0 !NR
    allocate(IST%v_ocn_C(SZI_(G), SZJB_(G))) ; IST%v_ocn_C(:,:) = 0.0 !NR
    allocate(IST%flux_u_top_Cu(SZIB_(G), SZJ_(G), 0:CatIce)) ; IST%flux_u_top_Cu(:,:,:) = 0.0 !NR
    allocate(IST%flux_v_top_Cv(SZI_(G), SZJB_(G), 0:CatIce)) ; IST%flux_v_top_Cv(:,:,:) = 0.0 !NR
  else
    allocate(IST%u_ice(SZIB_(G), SZJB_(G))) ; IST%u_ice(:,:) = 0.0
    allocate(IST%v_ice(SZIB_(G), SZJB_(G))) ; IST%v_ice(:,:) = 0.0
    allocate(IST%u_ocn(SZIB_(G), SZJB_(G))) ; IST%u_ocn(:,:) = 0.0 !NR
    allocate(IST%v_ocn(SZIB_(G), SZJB_(G))) ; IST%v_ocn(:,:) = 0.0 !NR
    allocate(IST%flux_u_top_bgrid(SZIB_(G), SZJB_(G), 0:CatIce)) ; IST%flux_u_top_bgrid(:,:,:) = 0.0 !NR
    allocate(IST%flux_v_top_bgrid(SZIB_(G), SZJB_(G), 0:CatIce)) ; IST%flux_v_top_bgrid(:,:,:) = 0.0 !NR
  endif

  ! Now register some of these arrays to be read from the restart files.
  domain => G%domain%mpp_domain
  idr = register_restart_field(Ice_restart, restart_file, 'part_size', IST%part_size, domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 't_surf', IST%t_surf, &
                               domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 'h_snow', IST%h_snow, &
                               domain=domain)
  idr = register_restart_field(Ice_restart, restart_file, 't_snow', IST%t_snow, &
                               domain=domain, mandatory=.false.)
  idr = register_restart_field(Ice_restart, restart_file, 'h_ice',  IST%h_ice, &
                               domain=domain)
  do n=1,G%NkIce
    write(nstr, '(I4)') n ; nstr = adjustl(nstr)
    idr = register_restart_field(Ice_restart, restart_file, 't_ice'//trim(nstr), &
                                 IST%t_ice(:,:,:,n), domain=domain, mandatory=(n==1))
    idr = register_restart_field(Ice_restart, restart_file, 'sal_ice'//trim(nstr), &
                                 IST%sal_ice(:,:,:,n), domain=domain, mandatory=.false.)
  enddo

  if (IST%Cgrid_dyn) then
    idr = register_restart_field(Ice_restart, restart_file, 'u_ice_C', IST%u_ice_C, &
                                 domain=domain, position=EAST, mandatory=.false.)
    idr = register_restart_field(Ice_restart, restart_file, 'v_ice_C', IST%v_ice_C, &
                                 domain=domain, position=NORTH, mandatory=.false.)
  else
    idr = register_restart_field(Ice_restart, restart_file, 'u_ice',   IST%u_ice, &
                                 domain=domain, position=CORNER, mandatory=.false.)
    idr = register_restart_field(Ice_restart, restart_file, 'v_ice',   IST%v_ice, &
                                 domain=domain, position=CORNER, mandatory=.false.)
  endif
  idr = register_restart_field(Ice_restart, restart_file, 'coszen', IST%coszen, &
                               domain=domain, mandatory=.false.)

end subroutine ice_state_register_restarts

subroutine dealloc_Ice_arrays(Ice)
  type(ice_data_type), intent(inout) :: Ice

  deallocate(Ice%mask, Ice%ice_mask, Ice%t_surf, Ice%s_surf)
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

  deallocate(IST%t_surf, IST%s_surf, IST%t_ocn, IST%sea_lev)
  deallocate(IST%part_size)
  if (IST%Cgrid_dyn) then
    deallocate(IST%u_ice_C, IST%v_ice_C, IST%u_ocn_C, IST%v_ocn_C)
    deallocate(IST%flux_u_top_Cu, IST%flux_v_top_Cv)
  else
    deallocate(IST%u_ocn, IST%v_ocn, IST%u_ice, IST%v_ice)
    deallocate(IST%flux_u_top_bgrid, IST%flux_v_top_bgrid)
  endif

  deallocate(IST%flux_u_top, IST%flux_v_top )
  deallocate(IST%flux_t_top, IST%flux_q_top, IST%flux_lw_top)
  deallocate(IST%flux_lh_top, IST%lprec_top, IST%fprec_top)
  deallocate(IST%flux_sw_vis_dir_top, IST%flux_sw_vis_dif_top)
  deallocate(IST%flux_sw_nir_dir_top, IST%flux_sw_nir_dif_top)

  deallocate(IST%Enth_Mass_in_atm, IST%Enth_Mass_out_atm)
  deallocate(IST%Enth_Mass_in_ocn, IST%Enth_Mass_out_ocn)
  deallocate(IST%flux_t_ocn_top, IST%flux_q_ocn_top)
  deallocate(IST%flux_lw_ocn_top, IST%flux_lh_ocn_top)
  deallocate(IST%flux_sw_vis_dir_ocn, IST%flux_sw_vis_dif_ocn)
  deallocate(IST%flux_sw_nir_dir_ocn, IST%flux_sw_nir_dif_ocn)
  deallocate(IST%lprec_ocn_top, IST%fprec_ocn_top)

  deallocate(IST%lwdn, IST%swdn, IST%coszen, IST%frazil, IST%frazil_input)
  deallocate(IST%bheat, IST%tmelt, IST%bmelt)
  deallocate(IST%sw_abs_sfc, IST%sw_abs_snow, IST%sw_abs_ice)
  deallocate(IST%sw_abs_ocn, IST%sw_abs_int)

  deallocate(IST%h_snow, IST%t_snow, IST%h_ice, IST%t_ice)
  deallocate(IST%enth_snow, IST%enth_ice, IST%sal_ice)

end subroutine dealloc_IST_arrays

subroutine IST_chksum(mesg, IST, G, haloshift)
  character(len=*),        intent(in)    :: mesg
  type(ice_state_type),    intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G
  integer, optional,       intent(in)    :: haloshift
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      IST - The ice state type variable to be checked.
!  (in)      G - The ocean's grid structure.
!  (in,opt)  haloshift - If present, check halo points out this far.
  character(len=20) :: k_str1, k_str
  integer :: hs, k

  ! Note that for the chksum calls to be useful for reproducing across PE
  ! counts, there must be no redundant points, so all variables use is..ie
  ! and js...je as their extent.
  hs=0; if (present(haloshift)) hs=haloshift

  call hchksum(IST%part_size, trim(mesg)//" IST%part_size",G,haloshift=hs)
  call hchksum(IST%h_ice, trim(mesg)//" IST%h_ice",G,haloshift=hs)
  do k=1,G%NkIce
    write(k_str1,'(I8)') k
    k_str = "("//trim(adjustl(k_str1))//")"
    call hchksum(IST%t_ice(:,:,:,k), trim(mesg)//" IST%h_ice("//trim(k_str),G,haloshift=hs)
    call hchksum(IST%enth_ice(:,:,:,k), trim(mesg)//" IST%enth_ice("//trim(k_str),G,haloshift=hs)
    call hchksum(IST%sal_ice(:,:,:,k), trim(mesg)//" IST%sal_ice("//trim(k_str),G,haloshift=hs)
  enddo
  call hchksum(IST%h_snow, trim(mesg)//" IST%h_snow",G,haloshift=hs)
  call hchksum(IST%t_snow, trim(mesg)//" IST%t_snow",G,haloshift=hs)
  call hchksum(IST%enth_snow(:,:,:,1), trim(mesg)//" IST%enth_snow",G,haloshift=hs)
  if (associated(IST%u_ice)) call Bchksum(IST%u_ice, mesg//" IST%u_ice",G,haloshift=hs)
  if (associated(IST%v_ice)) call Bchksum(IST%v_ice, mesg//" IST%v_ice",G,haloshift=hs)
  call check_redundant_B(mesg//" IST%u/v_ice", IST%u_ice, IST%v_ice, G)
  if (IST%Cgrid_dyn) then
    call uchksum(IST%u_ice_C, mesg//" IST%u_ice_C",G,haloshift=hs)
    call vchksum(IST%v_ice_C, mesg//" IST%v_ice_C",G,haloshift=hs)
    call check_redundant_C(mesg//" IST%u/v_ice_C", IST%u_ice_C, IST%v_ice_C, G)
  endif

end subroutine IST_chksum

subroutine Ice_public_type_chksum(mesg, Ice)
  character(len=*),        intent(in)    :: mesg
  type(ice_data_type),     intent(inout) :: Ice
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
  type(ice_data_type),     intent(in) :: Ice
  type(sea_ice_grid_type), intent(inout) :: G
  character(len=*),        intent(in) :: msg

  character(len=512) :: mesg1, mesg2
  integer :: i, j, k, l, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  integer :: n_bad, i_bad, j_bad, k_bad
  real    :: t_min, t_max

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc

  n_bad = 0 ; i_bad = 0 ; j_bad = 0 ; k_bad = 0

  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    if ((Ice%s_surf(i2,j2) < 0.0) .or. (Ice%s_surf(i2,j2) > 100.0)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; endif
    endif
    if ((abs(Ice%flux_t(i2,j2)) > 1e4) .or. (abs(Ice%flux_lw(i2,j2)) > 1e4)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; endif
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

subroutine IST_bounds_check(IST, G, msg)
  type(ice_state_type),    intent(in) :: IST
  type(sea_ice_grid_type), intent(inout) :: G
  character(len=*),        intent(in) :: msg

  character(len=512) :: mesg1, mesg2
  real, dimension(SZI_(G),SZJ_(G)) :: sum_part_sz
  real    :: tsurf_min, tsurf_max, tice_min, tice_max, tOcn_min, tOcn_max
  real    :: enth_min, enth_max
  integer :: i, j, k, l, isc, iec, jsc, jec, ncat, i_off, j_off
  integer :: n_bad, i_bad, j_bad, k_bad


  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce

  n_bad = 0 ; i_bad = 0 ; j_bad = 0 ; k_bad = 0

  sum_part_sz(:,:) = 0.0
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    sum_part_sz(i,j) = sum_part_sz(i,j) + IST%part_size(i,j,k) 
  enddo ; enddo ; enddo

  tOcn_min = -100. ; tOcn_max = 60.
  do j=jsc,jec ; do i=isc,iec
    if ((abs(sum_part_sz(i,j) - 1.0) > 1.0e-5) .or. &
        (IST%s_surf(i,j) < 0.0) .or. (IST%s_surf(i,j) > 100.0) .or. &
        (IST%t_ocn(i,j) < tOcn_min) .or. (IST%t_ocn(i,j) > tOcn_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; endif
    endif
  enddo ; enddo
  tsurf_min = tOcn_min + T_0degC ; tsurf_max = tOcn_max + T_0degC
  tice_min = -100. ; tice_max = 1.0
  enth_min = enth_from_TS(tice_min, 0., IST%ITV)
  enth_max = enth_from_TS(tice_max, 0., IST%ITV)
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    if ((IST%t_surf(i,j,k) < tsurf_min) .or. (IST%t_surf(i,j,k) > tsurf_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; endif
    endif
  enddo ; enddo ; enddo

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if ((IST%h_ice(i,j,k) > 1000.0) .or. (IST%h_snow(i,j,k) > 1000.0) .or. &
        (IST%t_snow(i,j,k) < tice_min) .or. (IST%t_snow(i,j,k) > tice_max) .or. &
        (IST%enth_snow(i,j,k,1) < enth_min) .or. (IST%enth_snow(i,j,k,1) > enth_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; endif
    endif
  enddo ; enddo ; enddo

  do l=1,G%NkIce ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if ((IST%t_ice(i,j,k,l) < tice_min) .or. (IST%t_ice(i,j,k,l) > tice_max) .or. &
        (IST%enth_ice(i,j,k,l) < enth_min) .or. (IST%enth_ice(i,j,k,l) > enth_max) .or. &
        (IST%sal_ice(i,j,k,l) < 0.0) .or. (IST%sal_ice(i,j,k,l) > 1000.0)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; endif
    endif
  enddo ; enddo ; enddo ; enddo

  if (n_bad > 0) then
    i = i_bad ; j=j_bad ; k = k_bad
    write(mesg1,'(" at ", 2(F6.1)," or i,j,k = ",3i4,"; nbad = ",i6," on pe ",i4)') &
           G%geolonT(i,j), G%geolatT(i,j), i_bad, j_bad, k_bad, n_bad, pe_here()
    if (k_bad > 0) then
      write(mesg2,'("T_sfc = ",1pe12.4,", ps = ",1pe12.4)') IST%t_surf(i,j,k), IST%part_size(i,j,k)
    else
      write(mesg2,'("T_ocn = ",1pe12.4,", S_sfc = ",1pe12.4,", sum_ps = ",1pe12.4)') &
            IST%t_ocn(i,j), IST%s_surf(i,j), sum_part_sz(i,j)
    endif
    call SIS_error(WARNING, "Bad ice state "//trim(msg)//" ; "//trim(mesg1)//" ; "//trim(mesg2), all_print=.true.)
    if (k_bad > 0) then
      write(mesg1,'("hi/hs = ", 2(1pe12.4)," ts = ",1pe12.4," ti = ",1pe12.4)') &
             IST%h_ice(i,j,k), IST%h_snow(i,j,k), IST%T_snow(i,j,k), IST%T_ice(i,j,k,1)
      do l=2,G%NkIce
        write(mesg2,'(", ", 1pe12.4)') IST%T_ice(i,j,k,l)
        mesg1 = trim(mesg1)//trim(mesg2)
      enddo
      call SIS_error(WARNING, mesg1, all_print=.true.)
    endif
  endif

end subroutine IST_bounds_check

!#######################################################################
! <SUBROUTINE NAME="ice_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ice_model_restart(Ice, time_stamp)
  type(ice_data_type), intent(inout), optional :: Ice
  character(len=*),    intent(in), optional :: time_stamp

  if (present(Ice)) then
    call save_restart(Ice%Ice_restart, time_stamp)
    call icebergs_save_restart(Ice%icebergs)
  else
    ! This option is here only to accomodate an old and inappropriate interface.
    ! Redo the order of the arguments when done right.
    call save_restart(Ice_restart, time_stamp)
  endif

end subroutine ice_model_restart
! </SUBROUTINE>
!#######################################################################

subroutine ice_diagnostics_init(Ice, IST, G, diag, Time)
  type(ice_data_type),     intent(inout)    :: Ice
  type(ice_state_type),    intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G
  type(SIS_diag_ctrl),     intent(in)    :: diag
  type(time_type),         intent(inout) :: Time

  real, parameter       :: missing = -1e34
  integer               :: id_geo_lon, id_geo_lat, id_sin_rot, id_cos_rot, id_cell_area
  logical               :: sent
  integer :: i, j, k, isc, iec, jsc, jec, n
  character(len=8) :: nstr

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

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
  IST%id_mi       = register_SIS_diag_field('ice_model', 'MI', diag%axesT1, Time, &
               'ice mass', 'kg/m^2', missing_value=missing)
  IST%id_mib      = register_SIS_diag_field('ice_model', 'MIB', diag%axesT1, Time, &
               'ice + bergs mass', 'kg/m^2', missing_value=missing)
  IST%id_cn       = register_SIS_diag_field('ice_model', 'CN', diag%axesTc, Time, &
               'ice concentration', '0-1', missing_value=missing)
  IST%id_hs       = register_SIS_diag_field('ice_model', 'HS', diag%axesT1, Time, &
               'snow thickness', 'm-snow', missing_value=missing)
  IST%id_tsn      = register_SIS_diag_field('ice_model', 'TSN', diag%axesT1, Time, &
               'snow layer temperature', 'C',  missing_value=missing)
  IST%id_hi       = register_SIS_diag_field('ice_model', 'HI', diag%axesT1, Time, &
               'ice thickness', 'm-ice', missing_value=missing)
  IST%id_hio      = register_SIS_diag_field('ice_model', 'HIO', diag%axesT1, Time, &
               'ice thickness', 'm-ice', missing_value=missing)
  
  IST%id_t_iceav = register_SIS_diag_field('ice_model', 'T_bulkice', diag%axesT1, Time, &
               'Volume-averaged ice temperature', 'C', missing_value=missing)
  IST%id_s_iceav = register_SIS_diag_field('ice_model', 'S_bulkice', diag%axesT1, Time, &
               'Volume-averaged ice salinity', 'g/kg', missing_value=missing)
  call safe_alloc_ids_1d(IST%id_t, G%NkIce)
  call safe_alloc_ids_1d(IST%id_sal, G%NkIce)
  do n=1,G%NkIce
    write(nstr, '(I4)') n ; nstr = adjustl(nstr)
    IST%id_t(n)   = register_SIS_diag_field('ice_model', 'T'//trim(nstr), &
                 diag%axesT1, Time, 'ice layer '//trim(nstr)//' temperature', &
                 'C',  missing_value=missing)
    IST%id_sal(n)   = register_SIS_diag_field('ice_model', 'Sal'//trim(nstr), &
               diag%axesT1, Time, 'ice layer '//trim(nstr)//' salinity', &
               'g/kg',  missing_value=missing)
  enddo
  IST%id_ts       = register_SIS_diag_field('ice_model', 'TS', diag%axesT1, Time, &
               'surface temperature', 'C', missing_value=missing)
  IST%id_sh       = register_SIS_diag_field('ice_model','SH' ,diag%axesT1, Time, &
               'sensible heat flux', 'W/m^2',  missing_value=missing)
  IST%id_lh       = register_SIS_diag_field('ice_model','LH' ,diag%axesT1, Time, &
               'latent heat flux', 'W/m^2', missing_value=missing)
  IST%id_sw       = register_SIS_diag_field('ice_model','SW' ,diag%axesT1, Time, &
               'short wave heat flux', 'W/m^2', missing_value=missing)
  IST%id_lw       = register_SIS_diag_field('ice_model','LW' ,diag%axesT1, Time, &
               'long wave heat flux over ice', 'W/m^2', missing_value=missing)
  IST%id_snofl    = register_SIS_diag_field('ice_model','SNOWFL' ,diag%axesT1, Time, &
               'rate of snow fall', 'kg/(m^2*s)', missing_value=missing)
  IST%id_rain     = register_SIS_diag_field('ice_model','RAIN' ,diag%axesT1, Time, &
               'rate of rain fall', 'kg/(m^2*s)', missing_value=missing)
  IST%id_runoff   = register_SIS_diag_field('ice_model','RUNOFF' ,diag%axesT1, Time, &
               'liquid runoff', 'kg/(m^2*s)', missing_value=missing)
  IST%id_calving  = register_SIS_diag_field('ice_model','CALVING',diag%axesT1, Time, &
               'frozen runoff', 'kg/(m^2*s)', missing_value=missing)
  IST%id_runoff_hflx   = register_SIS_diag_field('ice_model','RUNOFF_HFLX' ,diag%axesT1, Time, &
               'liquid runoff sensible heat flux', 'W/m^2', missing_value=missing)
  IST%id_calving_hflx  = register_SIS_diag_field('ice_model','CALVING_HFLX',diag%axesT1, Time, &
               'frozen runoff sensible heat flux', 'W/m^2', missing_value=missing)
  IST%id_evap     = register_SIS_diag_field('ice_model','EVAP',diag%axesT1, Time, &
               'evaporation', 'kg/(m^2*s)', missing_value=missing)
  IST%id_saltf    = register_SIS_diag_field('ice_model','SALTF' ,diag%axesT1, Time, &
               'ice to ocean salt flux', 'kg/(m^2*s)', missing_value=missing)
  IST%id_sn2ic    = register_SIS_diag_field('ice_model','SN2IC'  ,diag%axesT1,Time, &
               'rate of snow to ice conversion', 'kg/(m^2*s)', missing_value=missing)
  IST%id_tmelt    = register_SIS_diag_field('ice_model','TMELT'  ,diag%axesT1, Time, &
               'upper surface melting energy flux', 'W/m^2', missing_value=missing)
  IST%id_bmelt    = register_SIS_diag_field('ice_model','BMELT'  ,diag%axesT1, Time, &
               'bottom surface melting energy flux', 'W/m^2', missing_value=missing)
  IST%id_bheat    = register_SIS_diag_field('ice_model','BHEAT'  ,diag%axesT1, Time, &
               'ocean to ice heat flux', 'W/m^2', missing_value=missing)
  IST%id_e2m      = register_SIS_diag_field('ice_model','E2MELT' ,diag%axesT1, Time, &
               'heat needed to melt ice', 'J/m^2', missing_value=missing)
  IST%id_frazil   = register_SIS_diag_field('ice_model','FRAZIL' ,diag%axesT1, Time, &
               'energy flux of frazil formation', 'W/m^2', missing_value=missing)
  IST%id_alb      = register_SIS_diag_field('ice_model','ALB',diag%axesT1, Time, &
               'surface albedo','0-1', missing_value=missing )
  IST%id_coszen   = register_SIS_diag_field('ice_model','coszen',diag%axesT1, Time, &
               'cosine of zenith','-1:1', missing_value=missing )
  IST%id_sw_abs_sfc= register_SIS_diag_field('ice_model','sw_abs_sfc',diag%axesT1, Time, &
               'SW frac. abs. at the ice surface','0:1', missing_value=missing )
  IST%id_sw_abs_snow= register_SIS_diag_field('ice_model','sw_abs_snow',diag%axesT1, Time, &
               'SW frac. abs. in snow','0:1', missing_value=missing )

  call safe_alloc_ids_1d(IST%id_sw_abs_ice, G%NkIce)
  do n=1,G%NkIce
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
  IST%id_xprt     = register_SIS_diag_field('ice_model','XPRT',diag%axesT1, Time, &
               'frozen water transport convergence', 'kg/(m^2*yr)', missing_value=missing)
  IST%id_lsrc     = register_SIS_diag_field('ice_model','LSRC', diag%axesT1, Time, &
               'frozen water local source', 'kg/(m^2*yr)', missing_value=missing)
  IST%id_lsnk     = register_SIS_diag_field('ice_model','LSNK',diag%axesT1, Time, &
               'frozen water local sink', 'kg/(m^2*yr)', missing_value=missing)
  IST%id_bsnk     = register_SIS_diag_field('ice_model','BSNK',diag%axesT1, Time, &
               'frozen water local bottom sink', 'kg/(m^2*yr)', missing_value=missing)
  IST%id_qfres    = register_SIS_diag_field('ice_model', 'QFLX_RESTORE_ICE', diag%axesT1, Time, &
               'Ice Restoring heat flux', 'W/m^2', missing_value=missing)
  IST%id_qflim    = register_SIS_diag_field('ice_model', 'QFLX_LIMIT_ICE', diag%axesT1, Time, &
               'Ice Limit heat flux', 'W/m^2', missing_value=missing)
  IST%id_strna    = register_SIS_diag_field('ice_model','STRAIN_ANGLE', diag%axesT1,Time, &
               'strain angle', 'none', missing_value=missing)
  IST%id_uo       = register_SIS_diag_field('ice_model', 'UO', diag%axesB1, Time, &
               'surface current - x component', 'm/s', missing_value=missing)
  IST%id_vo       = register_SIS_diag_field('ice_model', 'VO', diag%axesB1, Time, &
               'surface current - y component', 'm/s', missing_value=missing)
  IST%id_sw_vis   = register_SIS_diag_field('ice_model','SW_VIS' ,diag%axesT1, Time, &
               'visible short wave heat flux', 'W/m^2', missing_value=missing)
  IST%id_sw_dir   = register_SIS_diag_field('ice_model','SW_DIR' ,diag%axesT1, Time, &
               'direct short wave heat flux', 'W/m^2', missing_value=missing)
  IST%id_sw_dif   = register_SIS_diag_field('ice_model','SW_DIF' ,diag%axesT1, Time, &
               'diffuse short wave heat flux', 'W/m^2', missing_value=missing)
  IST%id_sw_vis_dir = register_SIS_diag_field('ice_model','SW_VIS_DIR' ,diag%axesT1, Time, &
               'visible direct short wave heat flux', 'W/m^2', missing_value=missing)
  IST%id_sw_vis_dif = register_SIS_diag_field('ice_model','SW_VIS_DIF' ,diag%axesT1, Time, &
               'visible diffuse short wave heat flux', 'W/m^2', missing_value=missing)
  IST%id_sw_nir_dir = register_SIS_diag_field('ice_model','SW_NIR_DIR' ,diag%axesT1, Time, &
               'near IR direct short wave heat flux', 'W/m^2', missing_value=missing)
  IST%id_sw_nir_dif = register_SIS_diag_field('ice_model','SW_NIR_DIF' ,diag%axesT1, Time, &
               'near IR diffuse short wave heat flux', 'W/m^2', missing_value=missing)

  !
  ! diagnostics for quantities produced outside the ice model
  !
  IST%id_swdn  = register_SIS_diag_field('ice_model','SWDN' ,diag%axesT1, Time, &
             'downward shortwave flux', 'W/m^2', missing_value=missing)
  IST%id_lwdn  = register_SIS_diag_field('ice_model','LWDN' ,diag%axesT1, Time, &
             'downward longwave flux', 'W/m^2', missing_value=missing)
  IST%id_ta    = register_SIS_diag_field('ice_model', 'TA', diag%axesT1, Time, &
             'surface air temperature', 'C', missing_value=missing)
  IST%id_slp   = register_SIS_diag_field('ice_model', 'SLP', diag%axesT1, Time, &
             'sea level pressure', 'Pa', missing_value=missing)
  IST%id_sst   = register_SIS_diag_field('ice_model', 'SST', diag%axesT1, Time, &
             'sea surface temperature', 'deg-C', missing_value=missing)
  IST%id_sss   = register_SIS_diag_field('ice_model', 'SSS', diag%axesT1, Time, &
             'sea surface salinity', 'psu', missing_value=missing)
  IST%id_ssh   = register_SIS_diag_field('ice_model', 'SSH', diag%axesT1, Time, &
             'sea surface height', 'm', missing_value=missing)
  IST%id_obi   = register_SIS_diag_field('ice_model', 'OBI', diag%axesT1, Time, &
       'ice observed', '0 or 1', missing_value=missing)

  !
  ! diagnostics that are specific to C-grid dynamics of the ice model
  !
  if (IST%Cgrid_dyn) then
    IST%id_fax = register_SIS_diag_field('ice_model', 'FA_X', diag%axesCu1, Time, &
               'Air stress on ice on C-grid - x component', 'Pa', missing_value=missing)
    IST%id_fay = register_SIS_diag_field('ice_model', 'FA_Y', diag%axesCv1, Time, &
               'Air stress on ice on C-grid - y component', 'Pa', missing_value=missing)
  else
    IST%id_fax = register_SIS_diag_field('ice_model', 'FA_X', diag%axesB1, Time, &
               'air stress on ice - x component', 'Pa', missing_value=missing)
    IST%id_fay = register_SIS_diag_field('ice_model', 'FA_Y', diag%axesB1, Time, &
               'air stress on ice - y component', 'Pa', missing_value=missing)
  endif

  if (id_sin_rot>0) call post_data(id_sin_rot, G%sin_rot, diag, is_static=.true.)
  if (id_cos_rot>0) call post_data(id_cos_rot, G%cos_rot, diag, is_static=.true.)
  if (id_geo_lon>0) call post_data(id_geo_lon, G%geoLonT, diag, is_static=.true.)
  if (id_geo_lat>0) call post_data(id_geo_lat, G%geoLatT, diag, is_static=.true.)
  if (id_cell_area>0) call post_data(id_cell_area, cell_area, diag, is_static=.true.)

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

  integer :: i, j, k, m, isc, iec, jsc, jec, ncat
  real :: icebergs_value

  value = 0.0
  if(.not.Ice%pe) return

  IST => Ice%Ice_state

  isc = Ice%G%isc ; iec = Ice%G%iec ; jsc = Ice%G%jsc ; jec = Ice%G%jec
  ncat = Ice%G%CatIce

  select case (index)

    case (ISTOCK_WATER)

      value = 0.0
      do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + (IST%Rho_ice*IST%h_ice(i,j,k) + IST%Rho_snow*IST%h_snow(i,j,k)) * &
               IST%part_size(i,j,k) * (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j))
      enddo ; enddo ; enddo

    case (ISTOCK_HEAT)
      !### Fix the heat stock calculation.
      value = 0.0
      do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
        if ((IST%part_size(i,j,k)>0.0.and.IST%h_ice(i,j,k)>0.0)) then
          if (IST%slab_ice) then
            value = value - (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j)) * IST%part_size(i,j,k) * &
                           IST%h_ice(i,j,1)*IST%Rho_ice*LI
          else
            value = value - (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j)) * IST%part_size(i,j,k) * &
                            e_to_melt(IST%h_snow(i,j,k), IST%t_snow(i,j,k), &
                                      IST%h_ice(i,j,k), IST%t_ice(i,j,k,1),  &
                                      IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3), &
                                      IST%t_ice(i,j,k,4) )
          endif
        endif
      enddo ; enddo ; enddo

    case (ISTOCK_SALT)
      !There is no salt in the snow.
      value = 0.0
      do m=1,Ice%G%NkIce ; do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + (IST%part_size(i,j,k) * (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j))) * &
            (0.001*IST%Rho_ice*IST%h_ice(i,j,k)) * IST%sal_ice(i,j,k,m)
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



subroutine check_ice_model_nml(param_file)

  type(param_file_type), intent(in) :: param_file
  !--- namelist interface --------------
  real    :: h_lo_lim       = 0.0        ! min ice thickness for temp. calc.
  integer :: num_part       = 6          ! number of ice grid partitions
                                         ! partition 1 is open water
                                         ! partitions 2 to num_part-1 are
                                         !   thickness limited ice categories
                                         ! partition num_part is unlimited ice
  logical :: atmos_winds = .true.        ! wind stress from atmosphere model over t points and has wrong sign
  logical :: slab_ice    = .false.       ! do old-style GFDL slab ice?
  logical :: spec_ice    = .false.       ! old-style GFDL slab ice with SST, ice thickness and conc. from data
  logical :: do_ice_restore  = .false.   ! restore sea-ice toward climatology
  logical :: do_ice_limit    = .false.   ! limit sea ice to max_ice_limit
  real    :: max_ice_limit   = 4.0       ! maximum sea ice height(m),
                                         ! if do_ice_limit is true
                                         ! TK: default chosen based on observed
                                         !     ice thickness data used by climate
                                         !     group, which range up to 7.31 m
  real    :: ice_restore_timescale = 5.0 ! time scale for restoring ice (days)
  logical :: conservation_check = .true. ! check for heat and h2o conservation
  logical :: slp2ocean          = .false.! apply sea level pressure to ocean surface
  logical :: verbose            = .false.! control printing message, will slow model down when turn true
  logical :: do_icebergs        = .false.! call iceberg code to modify calving field
  logical :: add_diurnal_sw     = .false.! apply an additional diurnal cycle to shortwave radiation
  logical :: do_sun_angle_for_alb = .false.! find the sun angle for ocean albed in the frame of the ice model
  integer :: layout(2)          = (/0, 0/)
  integer :: io_layout(2)       = (/0, 0/)

  real    :: R_ice=0., R_snw=0., R_pnd=0.
  logical :: do_deltaEdd = .true.
  character(len=128) :: mask_table = "INPUT/ice_mask_table"

  ! The following are archaic namelist variables. They are here only for error
  ! checking of attempts to use out-of-date namelist values.
  real    :: p0             = missing    ! ice strength parameter
  real    :: c0             = missing    ! another ice strength parameter
  real    :: cdw            = missing    ! water/ice drag coefficient
  real    :: wd_turn        = missing    ! water/ice drag turning angle
  real    :: channel_viscosity = missing ! viscosity used in one-cell wide channels to parameterize transport (m^2/s)
  real    :: smag_ocn          = missing ! Smagorinksy coefficient for viscosity (dimensionless)
  real    :: ssh_gravity       = missing ! Gravity parameter used in channel viscosity parameterization (m/s^2)
  real    :: chan_cfl_limit    = missing ! CFL limit for channel viscosity parameterization (dimensionless)
  integer :: nsteps_dyn     = miss_int   ! dynamics steps per slow timestep
  integer :: nsteps_adv     = miss_int   ! advection steps per slow timestep
  real    :: mom_rough_ice  = missing    ! momentum same, cd10=(von_k/ln(10/z0))^2
  real    :: heat_rough_ice = missing    ! heat roughness length
  real    :: kmelt          = missing    ! ocean/ice heat flux constant
  real    :: ks             = missing    ! snow conductivity (W/mK)
  real    :: alb_sno        = missing    ! snow albedo (less if melting)
  real    :: alb_ice        = missing    ! ice albedo (less if melting)
  real    :: pen_ice        = missing    ! part unreflected solar penetrates ice
  real    :: opt_dep_ice    = missing    ! ice optical depth
  real    :: t_range_melt   = missing    ! melt albedos scaled in over T range
  real    :: ice_bulk_salin = missing    ! ice bulk salinity (for ocean salt flux)

  namelist /ice_model_nml/ mom_rough_ice, heat_rough_ice, p0, c0, cdw, wd_turn,  &
                           kmelt, alb_sno, alb_ice, pen_ice, opt_dep_ice,        &
                           nsteps_dyn, nsteps_adv, num_part, atmos_winds,        &
                           slab_ice, spec_ice, ice_bulk_salin, layout,           &
                           do_ice_restore, do_ice_limit, max_ice_limit,          &
                           ice_restore_timescale, slp2ocean, conservation_check, &
                           t_range_melt, ks, h_lo_lim, verbose,        &
                           do_icebergs, add_diurnal_sw, io_layout, channel_viscosity,&
                           smag_ocn, ssh_gravity, chan_cfl_limit, do_sun_angle_for_alb, &
                           mask_table, do_deltaEdd,R_ice,R_snw,R_pnd
  integer :: io, unit, ierr

#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ice_model_nml, iostat=io)
#else
  unit = open_namelist_file()
  read  (unit, ice_model_nml,iostat=io)
  call close_file (unit)
#endif
  ierr = check_nml_error(io,'ice_model_nml')

  ! Check for parameters that are still being set via the namelist but which
  ! should be set via the param_file instead.

  call archaic_nml_check(param_file, "USE_SLAB_ICE", "SLAB_ICE", slab_ice, .false.)
  call archaic_nml_check(param_file, "DO_ICEBERGS", "do_icebergs", do_icebergs, .false.)
  call archaic_nml_check(param_file, "SPECIFIED_ICE", "spec_ice", spec_ice, .false.)
  call archaic_nml_check(param_file, "NSTEPS_DYN", "nsteps_dyn", nsteps_dyn, miss_int, 432)
  call archaic_nml_check(param_file, "ICE_STRENGTH_PSTAR", "p0", p0, missing)
  call archaic_nml_check(param_file, "ICE_STRENGTH_CSTAR", "c0", c0, missing)
  call archaic_nml_check(param_file, "ICE_CDRAG_WATER", "cdw", cdw, missing)
  call archaic_nml_check(param_file, "AIR_WATER_STRESS_TURN_ANGLE", "wd_turn", wd_turn, missing)
  call archaic_nml_check(param_file, "NSTEPS_ADV", "nsteps_adv", nsteps_adv, miss_int, 1)
  call archaic_nml_check(param_file, "ICE_CHANNEL_VISCOSITY", &
                         "channel_viscosity", channel_viscosity, missing, 0.0)
  call archaic_nml_check(param_file, "ICE_CHANNEL_SMAG_COEF", "smag_ocn", smag_ocn, missing)
  call archaic_nml_check(param_file, "ICE_CHANNEL_CFL_LIMIT", "chan_cfl_limit", chan_cfl_limit, missing)

  call archaic_nml_check(param_file, "ICE_BULK_SALINITY", "ice_bulk_salin", ice_bulk_salin, missing)
  call archaic_nml_check(param_file, "SNOW_ALBEDO", "alb_snow", alb_sno, missing)
  call archaic_nml_check(param_file, "ICE_ALBEDO", "alb_ice", alb_ice, missing)
  call archaic_nml_check(param_file, "MOMENTUM_ROUGH_ICE", "mom_rough_ice", mom_rough_ice, missing)
  call archaic_nml_check(param_file, "HEAT_ROUGH_ICE", "heat_rough_ice", heat_rough_ice, missing)
  call archaic_nml_check(param_file, "ICE_KMELT", "kmelt", kmelt, missing)
  call archaic_nml_check(param_file, "SNOW_CONDUCT", "ks", ks, missing)
  call archaic_nml_check(param_file, "ICE_SW_PEN_FRAC", "pen_ice", pen_ice, missing)
  call archaic_nml_check(param_file, "ICE_OPTICAL_DEPTH", "opt_dep_ice", opt_dep_ice, missing)
  call archaic_nml_check(param_file, "ALBEDO_T_MELT_RANGE", "t_range_melt", t_range_melt, missing)
  call archaic_nml_check(param_file, "ICE_SEES_ATMOS_WINDS", "atmos_winds", atmos_winds, .true.)
  call archaic_nml_check(param_file, "DO_ICE_RESTORE", "do_ice_restore", do_ice_restore, .false.)
  call archaic_nml_check(param_file, "APPLY_ICE_LIMIT", "do_ice_limit", do_ice_limit, .false.)
  call archaic_nml_check(param_file, "MAX_ICE_THICK_LIMIT", "max_ice_limit", max_ice_limit, 4.0)
  call archaic_nml_check(param_file, "ICE_RESTORE_TIMESCALE", "ice_restore_timescale", ice_restore_timescale, 5.0)
  call archaic_nml_check(param_file, "APPLY_SLP_TO_OCEAN", "slp2ocean", slp2ocean, .false.)
  call archaic_nml_check(param_file, "MIN_H_FOR_TEMP_CALC", "h_lo_lim", h_lo_lim, 0.0)
  call archaic_nml_check(param_file, "ADD_DIURNAL_SW", "add_diurnal_sw", add_diurnal_sw, .false.)
  call archaic_nml_check(param_file, "DO_SUN_ANGLE_FOR_ALB", "do_sun_angle_for_alb", do_sun_angle_for_alb, .false.)
  call archaic_nml_check(param_file, "DO_DELTA_EDDINGTON_SW", "do_deltaEdd", do_deltaEdd, .true.)
  call archaic_nml_check(param_file, "ICE_DELTA_EDD_R_ICE", "R_ice", R_ice, 0.0)
  call archaic_nml_check(param_file, "ICE_DELTA_EDD_R_SNOW", "R_snw", R_snw, 0.0)
  call archaic_nml_check(param_file, "ICE_DELTA_EDD_R_POND", "R_pnd", R_pnd, 0.0)

end subroutine check_ice_model_nml

end module ice_type_mod
