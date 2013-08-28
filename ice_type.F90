!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_type_mod - maintains the sea ice data, reads/writes restarts, reads the  !
!                namelist and initializes diagnostics. - Mike Winton           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_type_mod

  use mpp_mod,          only: mpp_sum, mpp_clock_id, CLOCK_COMPONENT, &
                              CLOCK_LOOP, CLOCK_ROUTINE, stdout,input_nml_file
  use mpp_domains_mod,  only: domain2D, CORNER
  use mpp_domains_mod,  only: CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
  use mpp_domains_mod,  only: mpp_get_compute_domain
  use fms_mod,          only: file_exist, open_namelist_file, check_nml_error, write_version_number,&
                              read_data, close_file, field_exist, &
                              stderr, stdlog, error_mesg, FATAL, WARNING, NOTE, clock_flag_default
  use fms_io_mod,       only: save_restart, restore_state, query_initialized, &
                              register_restart_field, restart_file_type, set_domain, nullify_domain, &
                              parse_mask_table
  use time_manager_mod, only: time_type, time_type_to_real
  use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type
  use constants_mod,    only: Tfreeze, radius, pi

use ice_grid_mod,     only: set_ice_grid, ice_grid_end, sea_ice_grid_type
  use ice_grid_mod,     only: Domain
  use ice_grid_mod,     only: cell_area

  use ice_thm_mod,      only: ice_thm_param, DI, DS, e_to_melt
use ice_dyn_mod,       only: ice_dyn_init, ice_dyn_CS, ice_dyn_register_restarts, ice_dyn_end
use ice_transport_mod, only: ice_transport_init, ice_transport_CS, ice_transport_end ! , ice_transport_register_restarts
  use constants_mod,    only: LI => hlf ! latent heat of fusion - 334e3 J/(kg-ice)
  use ice_bergs,        only: icebergs_init, icebergs_end, icebergs, icebergs_stock_pe
  use ice_bergs,        only: icebergs_save_restart
  use astronomy_mod,    only: astronomy_init, astronomy_end
  use ice_shortwave_dEdd, only: shortwave_dEdd0_set_params

use MOM_domains,     only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : open_param_file, close_param_file
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg, is_root_pe
use SIS_diag_mediator, only : SIS_diag_ctrl, set_SIS_axes_info, SIS_diag_mediator_init
use SIS_diag_mediator, only : post_data=>post_SIS_data, diag_axis_init, register_static_field
use SIS_diag_mediator, only : register_SIS_diag_field
use SIS_get_input, only : Get_SIS_input, directories, archaic_nml_check

implicit none ; private

#include <SIS2_memory.h>

public :: ice_data_type, ice_state_type, ice_model_init, ice_model_end, ice_stock_pe, &
          ice_model_restart, ice_data_type_chksum
public :: do_sun_angle_for_alb

public  :: iceClock,iceClock1,iceClock2,iceClock3,iceClock4,iceClock5,iceClock6,iceClock7,iceClock8,iceClock9
public  :: iceClocka,iceClockb,iceClockc
public  :: earth_area

  real, parameter :: earth_area = 4*PI*RADIUS*RADIUS !5.10064471909788E+14 m^2
  real, parameter :: missing = -1e34
  integer, parameter :: miss_int = -9999

  character(len=128) :: version = '$Id: ice_type.F90,v 1.1.2.1.6.1.2.2.2.1 2013/06/18 22:24:14 nnz Exp $'
  character(len=128) :: tagname = '$Name: siena_201305_ice_sis2_5layer_dEdd_nnz $'

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
  real    :: ice_bulk_salin = missing    ! ice bulk salinity (for ocean salt flux)!CICE value

  ! mask_table contains information for masking domain ( n_mask, layout and mask_list).
  !   A text file to specify n_mask, layout and mask_list to reduce number of processor
  !   usage by masking out some domain regions which contain all land points. 
  !   The default file name of mask_table is "INPUT/ice_mask_table". Please note that 
  !   the file name must begin with "INPUT/". The first 
  !   line of mask_table will be number of region to be masked out. The second line 
  !   of the mask_table will be the layout of the model. User need to set ice_model_nml
  !   variable layout to be the same as the second line of the mask table.
  !   The following n_mask line will be the position of the processor to be masked out.
  !   The mask_table could be created by tools check_mask. 
  !   For example the mask_table will be as following if n_mask=2, layout=4,6 and 
  !   the processor (1,2) and (3,6) will be masked out. 
  !     2
  !     4,6
  !     1,2
  !     3,6
  character(len=128) :: mask_table = "INPUT/ice_mask_table"


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

  real    :: hlim_dflt(8) = (/ 0.0, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! thickness limits 1...num_part-1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! This structure contains the ice model state, and is intended to be private   !
! to SIS2.  It is not to be shared with other components and modules, and may  !
! use different indexing conventions to other modules..                        !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type ice_state_type
  type (time_type)                   :: Time_Init, Time
  type (time_type)                   :: Time_step_fast, Time_step_slow
  integer                            :: avg_count
   logical                            :: pe

   logical, pointer, dimension(:,:)   :: mask                =>NULL() ! where ice can be
  real, pointer, dimension(:,:,:) :: &
    part_size =>NULL(), &  ! The fractional coverage of a grid cell by each ice
                           ! thickness category, nondim, 0 to 1.  Category 0 is
                           ! open ocean.  The sum of part_size is 1.
    part_size_uv =>NULL()  ! The equivalent of part_size for B-grid velocity
                           ! cells.  Nondim., and 0 to 1, sums to 1 on a cell.

  ! The following are the 6 variables that constitute the sea-ice state.
  real, pointer, dimension(:,:) :: &
    u_ice =>NULL(), & ! The pseudo-zonal and pseudo-meridional ice velocities
    v_ice =>NULL()    ! along the model's grid directions, in m s-1.  All
                      ! thickness categories are assumed to have the same
                      ! velocity.
  real, pointer, dimension(:,:,:) :: &
    h_snow =>NULL(), &  ! The thickness of the snow in each category, in m.
    h_ice =>NULL(), &   ! The thickness of the ice in each category, in m.
    t_snow =>NULL()     ! The temperture of the snow in each category, in degC.
  real, pointer, dimension(:,:,:,:) :: &
    t_ice =>NULL()      ! The temperature of the sea ice in each category and
                        ! fractional thickness layer, in degC.
  
  real,    pointer, dimension(:,:) :: &
    s_surf  =>NULL(), &    ! The ocean surface salinity in g/kg.
    u_ocn   =>NULL(), &    ! The ocean's zonal velocity in m s-1.
    v_ocn   =>NULL(), &    ! The ocean's meridional velocity in m s-1.
    sea_lev =>NULL()       ! The equivalent sea-level, after any non-levitating
                           ! ice has been converted to sea-water, as determined
                           ! by the ocean, in m.  Sea-ice only contributes by
                           ! applying pressure to the ocean that is then
                           ! (partially) converted back to its equivalent by the
                           ! ocean. 
  
  real,    pointer, dimension(:,:,:) :: &
    ! The 3rd dimension in each of the following is ice thickness category.
    t_surf              =>NULL(), & ! The surface temperature, in Kelvin.
    flux_u_top          =>NULL(), & ! The downward? flux of zonal and meridional
    flux_v_top          =>NULL(), & ! momentum on an A-grid in ???.
    flux_u_top_bgrid    =>NULL(), & ! The downward? flux of zonal and meridional
    flux_v_top_bgrid    =>NULL(), & ! momentum on a B-grid velocity point in ???.
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

  real,    pointer, dimension(:,:)   :: lwdn                =>NULL() ! Accumulated diagnostics of
  real,    pointer, dimension(:,:  ) :: swdn                =>NULL() ! downward long/shortwave
  real,    pointer, dimension(:,:,:) :: pen                 =>NULL()
  real,    pointer, dimension(:,:,:) :: trn                 =>NULL() ! ice optical parameters
  real,    pointer, dimension(:,:,:) :: sw_abs_sfc          =>NULL() ! frac abs sw abs @ surf.
  real,    pointer, dimension(:,:,:) :: sw_abs_snow         =>NULL() ! frac abs sw abs in snow
  real,    pointer, dimension(:,:,:,:) :: sw_abs_ice        =>NULL() ! frac abs sw abs in ice layers
  real,    pointer, dimension(:,:,:) :: sw_abs_ocn          =>NULL() ! frac abs sw abs in ocean
  real,    pointer, dimension(:,:,:) :: sw_abs_int          =>NULL() ! frac abs sw abs in ice interior
  real,    pointer, dimension(:,:)   :: coszen              =>NULL()
  real,    pointer, dimension(:,:,:) :: tmelt               =>NULL()
  real,    pointer, dimension(:,:,:) :: bmelt               =>NULL()

  real,    pointer, dimension(:,:)   :: frazil              =>NULL()
  real,    pointer, dimension(:,:)   :: bheat               =>NULL()
   real,    pointer, dimension(:,:)   :: mi                  =>NULL() ! This is needed for the wave model. It is introduced here,
                                                                      ! because flux_ice_to_ocean cannot handle 3D fields. This may be
								! removed, if the information on ice thickness can be derived from 
								! eventually from h_ice outside the ice module.
   logical, pointer, dimension(:,:)   :: maskmap             =>NULL() ! A pointer to an array indicating which
                                                                      ! logical processors are actually used for
                                                                      ! the ocean code. The other logical
                                                                      ! processors would be all land points and
                                                                      ! are not assigned to actual processors.
                                                                      ! This need not be assigned if all logical
                                                                      ! processors are used

  logical :: slab_ice  ! If true, do the old style GFDL slab ice.
  real :: Rho_ocean    ! The nominal density of sea water, in kg m-3.
  real :: Rho_ice      ! The nominal density of sea ice, in kg m-3.
  real :: Rho_snow     ! The nominal density of snow on sea ice, in kg m-3.
  logical :: do_icebergs    ! If true, use the Lagrangian iceberg code, which
                            ! modifies the calving field among other things.
  logical :: specified_ice  ! If true, the sea ice is specified and there is
                            ! no need for ice dynamics.
  logical :: conservation_check ! If true, check for heat, salt and h2o conservation.

                  !### FIX THESE COMMENTS.
  logical :: atmos_winds ! wind stress from atmosphere model over t points and has wrong sign
  real :: mom_rough_ice  ! momentum same, cd10=(von_k/ln(10/z0))^2, in m.
  real :: heat_rough_ice ! heat roughness length, in m.
  real :: kmelt          ! ocean/ice heat flux constant, W m-2 K-1.
  real :: k_snow         ! snow conductivity (W/mK)
  real :: ice_bulk_salin ! ice bulk salinity (for ocean salt flux), in kg/kg.
  real :: alb_snow       ! snow albedo (less if melting), nondim.
  real :: alb_ice        ! ice albedo (less if melting), nondim.
  real :: pen_ice      ! part unreflected solar penetrates ice, nondim.
  real :: opt_dep_ice  ! ice optical depth, in m-1.
  real :: t_range_melt ! melt albedos scaled in over T range, in deg C.
  logical :: do_ice_restore   ! restore sea-ice toward climatology
  real    :: ice_restore_timescale ! time scale for restoring ice (days)
  logical :: do_ice_limit   ! limit sea ice to max_ice_limit
  real    :: max_ice_limit  ! The maximum sea ice height, in m,
                            ! if do_ice_limit is true
  logical :: slp2ocean  ! If true, apply sea level pressure to ocean surface.
  real    :: h_lo_lim   ! The min ice thickness for temp. calc, in m.
  logical :: verbose    ! control printing message, will slow model down when turn true
  logical :: add_diurnal_sw ! If true, apply a synthetic diurnal cycle to the shortwave radiation.
  logical :: do_sun_angle_for_alb ! If true, find the sun angle for calculating
                                  ! the ocean albedo in the frame of the ice model.
  logical :: do_deltaEdd  ! If true, a delta-Eddington radiative transfer calculation
                          ! for the shortwave radiation within the sea-ice.
  real :: deltaEdd_R_ice  ! Mysterious delta-Eddington tuning parameters, unknown.
  real :: deltaEdd_R_snow ! Mysterious delta-Eddington tuning parameters, unknown.
  real :: deltaEdd_R_pond ! Mysterious delta-Eddington tuning parameters, unknown.

  real    :: h2o(4), heat(4), salt(4) ! for conservation analysis
                             ! 1 - initial ice h2o/heat content
                             ! 2 - h2o/heat flux down at top of ice
                             ! 3 - h2o/heat flux down at bottom of ice
                             ! 4 - final ice h2o/heat content

  logical :: do_init = .false. ! If true, there is still some initialization
                               ! that needs to be done.

!   type(coupler_3d_bc_type)   :: ocean_fields       ! array of fields used for additional tracers
!   type(coupler_2d_bc_type)   :: ocean_fluxes       ! array of fluxes used for additional tracers
!   type(coupler_3d_bc_type)   :: ocean_fluxes_top   ! array of fluxes for averaging

  integer :: id_cn=-1, id_hi=-1, id_hs=-1, id_tsn=-1, id_t1=-1
  integer :: id_t2=-1, id_t3=-1, id_t4=-1, id_ts=-1, id_hio=-1, id_mi=-1, id_sh=-1
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
  integer :: id_abs_int=-1, id_sw_abs_snow=-1, id_sw_abs_ice1=-1, id_sw_abs_ice2=-1
  integer :: id_sw_abs_ice3=-1, id_sw_abs_ice4=-1, id_sw_pen=-1, id_sw_trn=-1

  type(ice_dyn_CS), pointer       :: ice_dyn_CSp => NULL()
  type(ice_transport_CS), pointer :: ice_transport_CSp => NULL()
  type(SIS_diag_ctrl)         :: diag
  type(icebergs), pointer     :: icebergs => NULL()
!  type(sea_ice_grid_type), pointer :: G ! A structure containing metrics and grid info.
end type ice_state_type

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! This structure contains the ice model data (some used by calling routines);  !
! the third index is partition (1 is open water; 2 is ice cover)               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type ice_data_type !  ice_public_type
  type(domain2D)                     :: Domain
  type (time_type)                   :: Time
  logical                            :: pe
  integer, pointer, dimension(:)     :: pelist              =>NULL() ! Used for flux-exchange.
     logical, pointer, dimension(:,:)   :: mask                =>NULL() ! where ice can be
  logical, pointer, dimension(:,:,:) :: ice_mask            =>NULL() ! where ice actually is (Size only?)
  real,    pointer, dimension(:,:,:) :: part_size           =>NULL()
  real,    pointer, dimension(:,:,:) :: albedo              =>NULL()
  real,    pointer, dimension(:,:,:) :: albedo_vis_dir      =>NULL()
  real,    pointer, dimension(:,:,:) :: albedo_nir_dir      =>NULL()
  real,    pointer, dimension(:,:,:) :: albedo_vis_dif      =>NULL()
  real,    pointer, dimension(:,:,:) :: albedo_nir_dif      =>NULL()
  real,    pointer, dimension(:,:,:) :: rough_mom           =>NULL()
  real,    pointer, dimension(:,:,:) :: rough_heat          =>NULL()
  real,    pointer, dimension(:,:,:) :: rough_moist         =>NULL()
  real,    pointer, dimension(:,:,:) :: t_surf              =>NULL()
  real,    pointer, dimension(:,:,:) :: u_surf              =>NULL()
  real,    pointer, dimension(:,:,:) :: v_surf              =>NULL()
  real,    pointer, dimension(:,:)   :: s_surf              =>NULL()

  ! These arrays will be used to set the forcing for the ocean.
  real,    pointer, dimension(:,:  ) :: flux_u              =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_v              =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_t              =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_q              =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_lw             =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_sw_vis_dir     =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_sw_vis_dif     =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_sw_nir_dir     =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_sw_nir_dif     =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_lh             =>NULL()
  real,    pointer, dimension(:,:  ) :: lprec               =>NULL()
  real,    pointer, dimension(:,:  ) :: fprec               =>NULL()
  real,    pointer, dimension(:,:  ) :: p_surf              =>NULL()
  real,    pointer, dimension(:,:  ) :: runoff              =>NULL()
  real,    pointer, dimension(:,:  ) :: calving             =>NULL()
  real,    pointer, dimension(:,:  ) :: runoff_hflx         =>NULL()
  real,    pointer, dimension(:,:  ) :: calving_hflx        =>NULL()
  real,    pointer, dimension(:,:  ) :: flux_salt           =>NULL()

  real,    pointer, dimension(:,:)   :: area                =>NULL()
  real,    pointer, dimension(:,:)   :: mi                  =>NULL() ! This is needed for the wave model. It is introduced here,
                                                                     ! because flux_ice_to_ocean cannot handle 3D fields. This may be
								! removed, if the information on ice thickness can be derived from 
								! eventually from h_ice outside the ice module.
     logical, pointer, dimension(:,:)   :: maskmap           =>NULL() ! A pointer to an array indicating which
                                                                      ! logical processors are actually used for
                                                                      ! the ocean code. The other logical
                                                                      ! processors would be all land points and
                                                                      ! are not assigned to actual processors.
                                                                      ! This need not be assigned if all logical
                                                                      ! processors are used
  integer, dimension(3)              :: axes
  type(coupler_3d_bc_type)           :: ocean_fields       ! array of fields used for additional tracers
  type(coupler_2d_bc_type)           :: ocean_fluxes       ! array of fluxes used for additional tracers
  type(coupler_3d_bc_type)           :: ocean_fluxes_top   ! array of fluxes for averaging

!      type(ice_dyn_CS), pointer       :: ice_dyn_CSp => NULL()
!      type(ice_transport_CS), pointer :: ice_transport_CSp => NULL()
!      type(SIS_diag_ctrl)         :: diag
      type(icebergs), pointer     :: icebergs => NULL()
  type(sea_ice_grid_type), pointer :: G ! A structure containing metrics and grid info.
  type(ice_state_type), pointer :: Ice_state => NULL() ! A structure containing the internal
                               ! representation of the ice state.
end type ice_data_type !  ice_public_type

  integer :: iceClock, iceClock1, iceCLock2, iceCLock3, iceClock4, iceClock5, &
             iceClock6, iceClock7, iceClock8, iceClock9, iceClocka, iceClockb, iceClockc
  type(restart_file_type), save :: Ice_restart

  contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_stock_pe - returns stocks of heat, water, etc. for conservation checks   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_stock_pe(Ice, index, value)

  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT

  type(ice_data_type) :: Ice
  integer, intent(in) :: index
  real, intent(out)   :: value
  type(ice_state_type), pointer :: IST => NULL()

  integer :: i, j, k, isc, iec, jsc, jec, ncat
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
        value = value + (DI*IST%h_ice(i,j,k) + DS*IST%h_snow(i,j,k)) * &
               IST%part_size(i,j,k) * (ICE%G%areaT(i,j)*Ice%G%mask2dT(i,j))
      enddo ; enddo ; enddo

    case (ISTOCK_HEAT)

      value = 0.0
      do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
        if ((IST%part_size(i,j,k)>0.0.and.IST%h_ice(i,j,k)>0.0)) then
          if (slab_ice) then
            value = value - (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j)) * IST%part_size(i,j,k) * &
                           IST%h_ice(i,j,2)*DI*LI
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
      !No salt in the h_snow component.
      value = 0.0
      do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + (DI*IST%h_ice(i,j,k)) * IST%ice_bulk_salin * &
               IST%part_size(i,j,k) * (Ice%G%areaT(i,j)*Ice%G%mask2dT(i,j))
      enddo ; enddo ; enddo

    case default

      value = 0.0

  end select

  if (IST%do_icebergs) then
    call icebergs_stock_pe(Ice%icebergs, index, icebergs_value)
    value = value + icebergs_value
  endif

end subroutine ice_stock_pe

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_init - initializes ice model data, parameters and diagnostics      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_init (Ice, Time_Init, Time, Time_step_fast, Time_step_slow )

  type (ice_data_type), intent(inout) :: Ice
  type (time_type)    , intent(in)    :: Time_Init      ! starting time of model integration
  type (time_type)    , intent(in)    :: Time           ! current time
    type (time_type)    , intent(in)    :: Time_step_fast ! time step for the ice_model_fast
    type (time_type)    , intent(in)    :: Time_step_slow ! time step for the ice_model_slow

  logical :: x_cyclic, tripolar_grid
    integer           :: io, ierr, nlon, nlat, npart, unit, log_unit, k
  integer :: sc, dy, i, j, l, i2, j2, k2, i_off, j_off
  integer :: isc, iec, jsc, jec
  integer :: CatIce
    integer           :: id_restart, id_restart_flux_sw
    real              :: dt_slow
    character(len=64) :: restart_file
  character(len=40)  :: mod = "ice_model" ! This module's name.
    integer           :: stdlogunit, stdoutunit
  type(param_file_type) :: param_file
  type(ice_state_type),    pointer :: IST => NULL()
  type(sea_ice_grid_type), pointer :: G => NULL()

    stdlogunit=stdlog()
    stdoutunit = stdout()

  if (associated(Ice%Ice_state)) then
    call SIS_error(WARNING, "ice_model_init called with an associated "// &
                    "Ice%Ice_state structure. Model is already initialized.")
    return
  endif
  allocate(Ice%Ice_state)
  IST => Ice%Ice_state
  allocate(Ice%G)
  G => Ice%G

  ! read namelist and write to logfile
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, nml=ice_model_nml, iostat=io)
#else
  unit = open_namelist_file()
  read  (unit, ice_model_nml,iostat=io)
  call close_file (unit)
#endif
  ierr = check_nml_error(io,'ice_model_nml')
  write (stdoutunit,'(/)')
  write (stdoutunit, ice_model_nml)
  write (stdlogunit, ice_model_nml)

  call Get_SIS_Input(param_file)

  call write_version_number( version, tagname )

  ! Check for parameters that are still being set via the namelist but which
  ! should be set via the param_file instead.  This has to be somewhere that
  ! has all of the namelist variables in scope.
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
  call archaic_nml_check(param_file, "ICE_CONSERVATION_CHECK", "conservation_check", conservation_check, .true.)
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

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "SPECIFIED_ICE", IST%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  call get_param(param_file, mod, "RHO_OCEAN", IST%Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0)
  call get_param(param_file, mod, "RHO_ICE", IST%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0)
  call get_param(param_file, mod, "RHO_SNOW", IST%Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0)
  call get_param(param_file, mod, "USE_SLAB_ICE", IST%slab_ice, &
                 "If true, use the very old slab-style ice.", default=.false.)
  call get_param(param_file, mod, "MOMENTUM_ROUGH_ICE", IST%mom_rough_ice, &
                 "The default momentum roughness length scale for the ocean.", &
                 units="m", default=1.0e-4)
  call get_param(param_file, mod, "HEAT_ROUGH_ICE", IST%heat_rough_ice, &
                 "The default roughness length scale for the turbulent \n"//&
                 "transfer of heat into the ocean.", units="m", default=1.0e-4)
  call get_param(param_file, mod, "ICE_KMELT", IST%kmelt, &
                 "A constant giving the proportionality of the ocean/ice \n"//&
                 "base heat flux to the tempature difference, given by \n"//&
                 "the product of the heat capacity per unit volume of sea \n"//&
                 "water times a molecular diffusive piston velocity.", &
                 units="W m-2 K-1", default=6e-5*4e6)
  call get_param(param_file, mod, "SNOW_CONDUCT", IST%k_snow, &
                 "The conductivity of heat in snow.", units="W m-1 K-1", &
                 default=0.31)
  call get_param(param_file, mod, "SNOW_ALBEDO", IST%alb_snow, &
                 "The albedo of dry snow atop sea ice.", units="nondim", &
                 default=0.85)
  call get_param(param_file, mod, "ICE_ALBEDO", IST%alb_ice, &
                 "The albedo of dry bare sea ice.", units="nondim", &
                 default=0.5826)
  call get_param(param_file, mod, "ICE_SW_PEN_FRAC", IST%pen_ice, &
                 "The fraction of the unreflected shortwave radiation that \n"//&
                 "penetrates into the ice.", units="Nondimensional", default=0.3)
  call get_param(param_file, mod, "ICE_OPTICAL_DEPTH", IST%opt_dep_ice, &
                 "The optical depth of shortwave radiation in sea ice.", &
                 units="m", default=0.67)
  call get_param(param_file, mod, "ALBEDO_T_MELT_RANGE", IST%t_range_melt, &
                 "The temperature range below freezing over which the \n"//&
                 "albedos are changed by partial melting.", units="degC", &
                 default=1.0)
  call get_param(param_file, mod, "ICE_CONSERVATION_CHECK", IST%conservation_check, &
                 "If true, do additional calculations to check for \n"//&
                 "internal conservation of heat, salt, and water mass in \n"//&
                 "the sea ice model.  This does not change answers, but \n"//&
                 "can increase model run time.", default=.true.)
  call get_param(param_file, mod, "ICE_SEES_ATMOS_WINDS", IST%atmos_winds, &
                 "If true, the sea ice is being given wind stresses with \n"//&
                 "the atmospheric sign convention, and need to have their \n"//&
                 "sign changed.", default=.true.)
  call get_param(param_file, mod, "ICE_BULK_SALINITY", IST%ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "kg/kg", default=0.004)
!  call get_param(param_file, mod, "ICE_LAYOUT?",  ###HANDLE THIS LIKE MOM6?
!  call get_param(param_file, mod, "ICE_IO_LAYOUT?",  ###HANDLE THIS LIKE MOM6?
!  call get_param(param_file, mod, "ICE_MASK_TABLE?",  ###HANDLE THIS LIKE MOM6?
  call get_param(param_file, mod, "DO_ICE_RESTORE", IST%do_ice_restore, &
                 "If true, restore the sea ice state toward climatology.", &
                 default=.false.)
  if (IST%do_ice_restore) &
    call get_param(param_file, mod, "ICE_RESTORE_TIMESCALE", IST%ice_restore_timescale, &
                 "The restoring timescale when DO_ICE_RESTORE is true.", &
                 units="days", default=5.0)
  call get_param(param_file, mod, "APPLY_ICE_LIMIT", IST%do_ice_limit, &
                 "If true, restore the sea ice state toward climatology.", &
                 default=.false.)
  if (IST%do_ice_limit) &
    call get_param(param_file, mod, "MAX_ICE_THICK_LIMIT", IST%max_ice_limit, &
                 "The maximum permitted sea ice thickness when \n"//&
                 "APPLY_ICE_LIMIT is true.", units="m", default=4.0)
  call get_param(param_file, mod, "APPLY_SLP_TO_OCEAN", IST%slp2ocean, &
                 "If true, apply the atmospheric sea level pressure to \n"//&
                 "the ocean.", default=.false.)
  call get_param(param_file, mod, "MIN_H_FOR_TEMP_CALC", IST%h_lo_lim, &
                 "The minimum ice thickness at which to do temperature \n"//&
                 "calculations.", units="m", default=0.0)
  call get_param(param_file, mod, "VERBOSE", IST%verbose, &
                 "If true, write out verbose diagnostics.", default=.false.)
  call get_param(param_file, mod, "DO_ICEBERGS", IST%do_icebergs, &
                 "If true, call the iceberg module.", default=.false.)
  call get_param(param_file, mod, "ADD_DIURNAL_SW", IST%add_diurnal_sw, &
                 "If true, add a synthetic diurnal cycle to the shortwave \n"//&
                 "radiation.", default=.false.)
  call get_param(param_file, mod, "DO_SUN_ANGLE_FOR_ALB", IST%do_sun_angle_for_alb, &
                 "If true, find the sun angle for calculating the ocean \n"//&
                 "albedo within the sea ice model.", default=.false.)
  call get_param(param_file, mod, "DO_DELTA_EDDINGTON_SW", IST%do_deltaEdd, &
                 "If true, a delta-Eddington radiative transfer calculation \n"//&
                 "for the shortwave radiation within the sea-ice.", default=.true.)
  call get_param(param_file, mod, "ICE_DELTA_EDD_R_ICE", IST%deltaEdd_R_ice, &
                 "A dreadfully documented tuning parameter for the radiative \n"//&
                 "propeties of sea ice with the delta-Eddington radiative \n"//&
                 "transfer calculation.", units="perhaps nondimensional?", default=0.0)
  call get_param(param_file, mod, "ICE_DELTA_EDD_R_SNOW", IST%deltaEdd_R_snow, &
                 "A dreadfully documented tuning parameter for the radiative \n"//&
                 "propeties of snow on sea ice with the delta-Eddington \n"//&
                 "radiative transfer calculation.", &
                 units="perhaps nondimensional?", default=0.0)
  call get_param(param_file, mod, "ICE_DELTA_EDD_R_POND", IST%deltaEdd_R_pond, &
                 "A dreadfully documented tuning parameter for the radiative \n"//&
                 "propeties of meltwater ponds on sea ice with the delta-Eddington \n"//&
                 "radiative transfer calculation.", units="perhaps nondimensional?", &
                 default=0.0)

  if (IST%specified_ice) then
    IST%slab_ice = .true.
       nsteps_dyn = 0
       nsteps_adv = 0
  end if
  if (IST%slab_ice) num_part = 2 ! open water and ice ... but never in same place

  dt_slow = time_type_to_real(Time_step_slow)
!    call get_time(Time_step_slow, sc, dy); dt_slow=864e2*dy+sc

  if (file_exist(mask_table)) then
     call SIS_mesg(' ice_model_init:  reading maskmap information from '//trim(mask_table))
     if(layout(1) == 0 .OR. layout(2) == 0 ) call error_mesg('ice_model_init', &
        'ice_model_nml layout should be set when file '//trim(mask_table)//' exists', FATAL)

     allocate(Ice%maskmap(layout(1), layout(2)))
     call parse_mask_table(mask_table, Ice%maskmap, "Ice model")
  endif

  if( ASSOCIATED(Ice%maskmap) ) then
    call set_ice_grid(Ice%G, param_file, Ice%domain, num_part, layout, io_layout, Ice%maskmap  )
  else
    call set_ice_grid(Ice%G, param_file, Ice%domain, num_part, layout, io_layout )
  endif
  if (IST%slab_ice) then
    G%CatIce = 1 ! open water and ice ... but never in same place
  endif

  ! Initialize G%H_cat_lim here.  ###This needs to be extended to add more options.
  do k=1,min(G%CatIce+1,size(hlim_dflt(:)))
    G%H_cat_lim(k) = hlim_dflt(k)
  enddo
  if ((G%CatIce+1 > size(hlim_dflt(:))) .and. (size(hlim_dflt(:)) > 1)) then
    do k=min(G%CatIce+1,size(hlim_dflt(:))) + 1, G%CatIce+1
      G%H_cat_lim(k) =  2.0*G%H_cat_lim(k-1) - G%H_cat_lim(k-2)
    enddo
  endif

  call set_domain(G%Domain%mpp_domain)
  CatIce = G%CatIce

!  call allocate_ice_data_type_arrays(Ice, G%Domain%mpp_domain, G%CatIce)

  ! Allocate the internally visible ice_state_type's arrays.
!  call allocate_state_data_type_arrays(IST, G)

  ! Allocate and register fields for restarts.
  restart_file = 'ice_model.res.nc'
  call ice_data_type_register_restarts(G%Domain%mpp_domain, G%CatIce, param_file, Ice, Ice_restart, restart_file)

  call ice_state_register_restarts(G, param_file, IST, Ice_restart, restart_file)

  call ice_dyn_register_restarts(Ice%G, param_file, IST%ice_dyn_CSp, Ice_restart, restart_file)
!    call ice_transport_register_restarts(Ice%G, param_file, IST%ice_transport_CSp, Ice_restart, restart_file)



  ! Redefine the computational domain sizes to use the ice model's indexing convention.
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc

  Ice%area(:,:)   = cell_area(:,:) * 4*PI*RADIUS*RADIUS  ! ### Eliminate later
  IST%coszen(:,:) = cos(3.14*67.0/180.0) ! NP summer solstice.

  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%mask(i2,j2) = ( Ice%G%mask2dT(i,j) > 0.5 )
!###   Ice%area(i2,j2) = G%areaT(i,j) * G%mask2dT(i,j)
  enddo ; enddo
 
  Ice%Time           = Time
  IST%Time           = Time
  IST%Time_Init      = Time_Init
  IST%Time_step_fast = Time_step_fast
  IST%Time_step_slow = Time_step_slow

  IST%avg_count      = 0

  !
  ! read restart
  !
  restart_file = 'INPUT/ice_model.res.nc'
  if (file_exist(restart_file)) then
    call restore_state(Ice_restart)

!    if ( .NOT.query_initialized(Ice_restart, id_restart_flux_sw) ) then
!      call error_mesg ('ice_model_init', &
!         'Restart file does not contain flux_sw_* subcomponents!', WARNING)
!      ! for compatibility with preN restarts, we should check for total SW flux 
!      if (field_exist(restart_file,'flux_sw')) then
!        call read_data( restart_file, 'flux_sw', Ice%flux_sw_vis_dir, domain) 
!        ! simplest way to break the total flux to 4 components
!        Ice%flux_sw_vis_dir(:,:) = Ice%flux_sw_vis_dir(:,:) / 4
!        Ice%flux_sw_vis_dif(:,:) = Ice%flux_sw_vis_dir(:,:)
!        Ice%flux_sw_nir_dir(:,:) = Ice%flux_sw_vis_dir(:,:)
!        Ice%flux_sw_nir_dif(:,:) = Ice%flux_sw_vis_dir(:,:)
!      else
!        call error_mesg ('ice_model_init', &
!             'Restart file does not contain flux_sw total or its components!', FATAL)
!      endif
!    endif

    !--- update the halo values.
    call pass_var(IST%part_size, Ice%G%Domain, complete=.false.)
    call pass_var(IST%h_ice, Ice%G%Domain, complete=.false.)
    call pass_var(IST%h_snow, Ice%G%Domain, complete=.false.)
    do l=1,G%NkIce
      call pass_var(IST%t_ice(:,:,:,l), Ice%G%Domain, complete=.false.)
    enddo
    call pass_var(IST%t_snow, Ice%G%Domain, complete=.true.)

    call pass_vector(IST%u_ice, IST%v_ice, Ice%G%Domain, stagger=BGRID_NE)
  else ! no restart implies initialization with no ice
    IST%part_size(:,:,:) = 0.0
    IST%part_size(:,:,0) = 1.0

    Ice%rough_mom(:,:,:)   = IST%mom_rough_ice
    Ice%rough_heat(:,:,:)  = IST%heat_rough_ice
    Ice%rough_moist(:,:,:) = IST%heat_rough_ice
    IST%t_surf(:,:,:) = Tfreeze-5.0
    IST%t_snow(:,:,:) = -5.0
    IST%t_ice(:,:,:,:) = -5.0

    IST%do_init = .true. ! Some more initilization needs to be done in ice_model.
  endif ! file_exist(restart_file)

  do k=0,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
    Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  IST%part_size_uv(:,:,:) = 0.0
  IST%part_size_uv(:,:,0) = 1.0
  !### ADD PARENTHESIS FOR REPRODUCIBILITY.
  do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
    if(Ice%G%mask2dBu(i,j) > 0.5 ) then
       IST%part_size_uv(i,j,k) = 0.25*(IST%part_size(i+1,j+1,k) + IST%part_size(i+1,j,k) + &
                                       IST%part_size(i,j+1,k) + IST%part_size(i,j,k))
    else
       IST%part_size_uv(i,j,k) = 0.0
    endif
    IST%part_size_uv(i,j,1) = IST%part_size_uv(i,j,1) - IST%part_size_uv(i,j,k)
  enddo ; enddo ; enddo

    
  call SIS_diag_mediator_init(Ice%G, param_file, IST%diag, component="SIS")
  call set_SIS_axes_info(Ice%G, param_file, IST%diag)

  call ice_diagnostics_init(Ice, IST, Ice%G, IST%diag, IST%Time)

  call ice_dyn_init(IST%Time, Ice%G, param_file, IST%diag, IST%ice_dyn_CSp)
  call ice_transport_init(IST%Time, Ice%G, param_file, IST%diag, IST%ice_transport_CSp)
  call ice_thm_param(IST%alb_snow, IST%alb_ice, IST%pen_ice, IST%opt_dep_ice, IST%slab_ice, &
                     IST%t_range_melt, IST%k_snow, IST%h_lo_lim, IST%do_deltaEdd)

  call close_param_file(param_file)

  iceClock = mpp_clock_id( 'Ice', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  iceClock1 = mpp_clock_id( 'Ice: bot to top', flags=clock_flag_default, grain=CLOCK_ROUTINE )
  iceClock2 = mpp_clock_id( 'Ice: update slow (dn)', flags=clock_flag_default, grain=CLOCK_ROUTINE )
  iceClock7 = mpp_clock_id( '  Ice: slow: conservation check', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock4 = mpp_clock_id( '  Ice: slow: dynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClocka = mpp_clock_id( '       slow: ice_dynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClockb = mpp_clock_id( '       slow: comm/cut check ', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClockc = mpp_clock_id( '       slow: diags', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock5 = mpp_clock_id( '  Ice: slow: thermodynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock6 = mpp_clock_id( '  Ice: slow: restore/limit', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock8 = mpp_clock_id( '  Ice: slow: salt to ocean', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock9 = mpp_clock_id( '  Ice: slow: thermodyn diags', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock3 = mpp_clock_id( 'Ice: update fast', flags=clock_flag_default, grain=CLOCK_ROUTINE )

  ! Initialize icebergs
  x_cyclic = (Ice%G%Domain%X_FLAGS == CYCLIC_GLOBAL_DOMAIN)
  tripolar_grid = (Ice%G%Domain%Y_FLAGS == FOLD_NORTH_EDGE)
  if (IST%do_icebergs) call icebergs_init(Ice%icebergs, &
           Ice%G%Domain%niglobal, Ice%G%Domain%njglobal, Ice%G%Domain%layout, Ice%G%Domain%io_layout, Ice%axes(1:2), Ice%maskmap, x_cyclic, tripolar_grid, &
           dt_slow, Time, Ice%G%geoLonBu(isc:iec,jsc:jec), Ice%G%geoLatBu(isc:iec,jsc:jec), &
           Ice%G%mask2dT, Ice%G%dxCv, Ice%G%dyCu, cell_area, Ice%G%cos_rot, Ice%G%sin_rot )

  if (IST%add_diurnal_sw .or. IST%do_sun_angle_for_alb) call astronomy_init

  call shortwave_dEdd0_set_params(IST%deltaEdd_R_ice,IST%deltaEdd_R_snow,IST%deltaEdd_R_pond)

  call nullify_domain()

end subroutine ice_model_init


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
  integer :: isc, iec, jsc, jec, km, id_restart

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
  id_restart = register_restart_field(Ice_restart, restart_file, 'albedo',    Ice%albedo,    domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'albedo_vis_dir', Ice%albedo_vis_dir, &
                                      domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'albedo_nir_dir', Ice%albedo_nir_dir, &
                                      domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'albedo_vis_dif', Ice%albedo_vis_dif, &
                                      domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'albedo_nir_dif', Ice%albedo_nir_dif, &
                                      domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'rough_mom',   Ice%rough_mom,   domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'rough_heat',  Ice%rough_heat,  domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'rough_moist', Ice%rough_moist, domain=domain)

  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_u',      Ice%flux_u,       domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_v',      Ice%flux_v,       domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_t',      Ice%flux_t,       domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_q',      Ice%flux_q,       domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_salt',   Ice%flux_salt,    domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_lw',     Ice%flux_lw,      domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'lprec',       Ice%lprec,        domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'fprec',       Ice%fprec,        domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'runoff',      Ice%runoff,       domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'calving',     Ice%calving,      domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'runoff_hflx', Ice%runoff_hflx,  domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'calving_hflx',Ice%calving_hflx, domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'p_surf',      Ice%p_surf,       domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_sw_vis_dir', Ice%flux_sw_vis_dir, &
                                      domain=domain)    
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_sw_vis_dif', Ice%flux_sw_vis_dif, &
                                      domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_sw_nir_dir', Ice%flux_sw_nir_dir, &
                                      domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_sw_nir_dif', Ice%flux_sw_nir_dif, &
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
  integer :: CatIce, id_restart

  CatIce = G%CatIce
  allocate(IST%t_surf(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%t_surf(:,:,:) = 0.0 !X
  allocate(IST%s_surf(SZI_(G), SZJ_(G))) ; IST%s_surf(:,:) = 0.0 !NI X
  allocate(IST%sea_lev(SZI_(G), SZJ_(G))) ; IST%sea_lev(:,:) = 0.0 !NR 
  allocate(IST%part_size(SZI_(G), SZJ_(G), 0:CatIce)) ; IST%part_size(:,:,:) = 0.0
  allocate(IST%part_size_uv(SZIB_(G), SZJB_(G), 0:CatIce)) ; IST%part_size_uv(:,:,:) = 0.0 !NR X
  allocate(IST%u_ocn(SZI_(G), SZJ_(G))) ; IST%u_ocn(:,:) = 0.0 !NR
  allocate(IST%v_ocn(SZI_(G), SZJ_(G))) ; IST%v_ocn(:,:) = 0.0 !NR
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
  allocate(IST%flux_u_top_bgrid(SZIB_(G), SZJB_(G), 0:CatIce)) ; IST%flux_u_top_bgrid(:,:,:) = 0.0 !NR
  allocate(IST%flux_v_top_bgrid(SZIB_(G), SZJB_(G), 0:CatIce)) ; IST%flux_v_top_bgrid(:,:,:) = 0.0 !NR

  allocate(IST%lwdn(SZI_(G), SZJ_(G))) ; IST%lwdn(:,:) = 0.0 !NR
  allocate(IST%swdn(SZI_(G), SZJ_(G))) ; IST%swdn(:,:) = 0.0 !NR
  allocate(IST%frazil(SZI_(G), SZJ_(G))) ; IST%frazil(:,:) = 0.0 !NR
  allocate(IST%bheat(SZI_(G), SZJ_(G))) ; IST%bheat(:,:) = 0.0 !NI
  allocate(IST%tmelt(SZI_(G), SZJ_(G), CatIce)) ; IST%tmelt(:,:,:) = 0.0 !NR
  allocate(IST%bmelt(SZI_(G), SZJ_(G), CatIce)) ; IST%bmelt(:,:,:) = 0.0 !NR

  allocate(IST%pen(SZI_(G), SZJ_(G), CatIce)) ; IST%pen(:,:,:) = 0.0 !NI
  allocate(IST%trn(SZI_(G), SZJ_(G), CatIce)) ; IST%trn(:,:,:) = 0.0 !NI
  allocate(IST%sw_abs_sfc(SZI_(G), SZJ_(G), CatIce)) ; IST%sw_abs_sfc(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_snow(SZI_(G), SZJ_(G), CatIce)) ; IST%sw_abs_snow(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_ice(SZI_(G), SZJ_(G), CatIce, G%NkIce)) ; IST%sw_abs_ice(:,:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_ocn(SZI_(G), SZJ_(G), CatIce)) ; IST%sw_abs_ocn(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_int(SZI_(G), SZJ_(G), CatIce)) ; IST%sw_abs_int(:,:,:) = 0.0 !NR

  allocate(IST%u_ice(SZIB_(G), SZJB_(G))) ; IST%u_ice(:,:) = 0.0
  allocate(IST%v_ice(SZIB_(G), SZJB_(G))) ; IST%v_ice(:,:) = 0.0
  allocate(IST%h_snow(SZI_(G), SZJ_(G), CatIce)) ; IST%h_snow(:,:,:) = 0.0
  allocate(IST%t_snow(SZI_(G), SZJ_(G), CatIce)) ; IST%t_snow(:,:,:) = 0.0
  allocate(IST%h_ice(SZI_(G), SZJ_(G), CatIce)) ; IST%h_ice(:,:,:) = 0.0
  allocate(IST%t_ice(SZI_(G), SZJ_(G), CatIce, G%NkIce)) ; IST%t_ice(:,:,:,:) = 0.0
  

  ! Now register some of these arrays to be read from the restart files.
  domain => G%domain%mpp_domain
  id_restart = register_restart_field(Ice_restart, restart_file, 'part_size', IST%part_size, domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_surf',    IST%t_surf,    domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'h_snow',    IST%h_snow,    domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_snow',    IST%t_snow,    domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'h_ice',     IST%h_ice,     domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_ice1',    IST%t_ice(:,:,:,1), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_ice2',    IST%t_ice(:,:,:,2), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_ice3',    IST%t_ice(:,:,:,3), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_ice4',    IST%t_ice(:,:,:,4), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'u_ice',     IST%u_ice,     domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'v_ice',     IST%v_ice,     domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'coszen',    IST%coszen,    domain=domain, mandatory=.false.)

end subroutine ice_state_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_end - writes the restart file and deallocates memory               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_end (Ice)
  type (ice_data_type), intent(inout) :: Ice
  integer           :: k

  integer           :: unit
  character(len=22) :: restart='RESTART/ice_model.res'
  type(ice_state_type), pointer :: IST => NULL()

  IST => Ice%Ice_state
  if (IST%conservation_check) call ice_print_budget(IST)

  call ice_model_restart()

  !--- release memory ------------------------------------------------
  call ice_grid_end(Ice%G)

  deallocate(Ice%mask, Ice%ice_mask, Ice%t_surf, Ice%s_surf, IST%t_surf, IST%s_surf, IST%sea_lev )
  deallocate(Ice%part_size, IST%part_size, IST%part_size_uv, Ice%u_surf, Ice%v_surf )
  deallocate(IST%u_ocn, IST%v_ocn ,  Ice%rough_mom, Ice%rough_heat )
  deallocate(Ice%rough_moist, Ice%albedo, IST%flux_u_top, IST%flux_v_top )
  deallocate(IST%flux_u_top_bgrid, IST%flux_v_top_bgrid )
  deallocate(IST%flux_t_top, IST%flux_q_top, IST%flux_lw_top )
  deallocate(IST%flux_lh_top, IST%lprec_top, IST%fprec_top, Ice%flux_u )
  deallocate(Ice%flux_v, Ice%flux_t, Ice%flux_q, Ice%flux_lw )
  deallocate(Ice%flux_lh, Ice%lprec, Ice%fprec, Ice%p_surf, Ice%runoff ) 
  deallocate(Ice%calving, Ice%runoff_hflx, Ice%calving_hflx )
  deallocate(Ice%flux_salt)
  deallocate(IST%lwdn, IST%swdn, IST%coszen)
  deallocate(IST%frazil )
  deallocate(IST%bheat, IST%u_ice, IST%v_ice )
  deallocate(IST%tmelt, IST%bmelt, IST%pen, IST%trn )
  deallocate(IST%h_snow, IST%t_snow, IST%h_ice )
  deallocate(IST%t_ice)
  deallocate(Ice%flux_sw_vis_dir, Ice%flux_sw_vis_dif )
  deallocate(Ice%flux_sw_nir_dir, Ice%flux_sw_nir_dif )

  call ice_dyn_end(IST%ice_dyn_CSp)
  call ice_transport_end(IST%ice_transport_CSp)

  ! End icebergs
  if (IST%do_icebergs) call icebergs_end(Ice%icebergs)
  
  deallocate(Ice%Ice_state)

  if (add_diurnal_sw .or. do_sun_angle_for_alb) call astronomy_end

end subroutine ice_model_end

subroutine ice_print_budget(IST)
  type(ice_state_type), intent(inout) :: IST
  integer :: k
 
  do k=1,4
    call mpp_sum(IST%h2o(k))
    call mpp_sum(IST%heat(k))
    call mpp_sum(IST%salt(k))
!          call mpp_sum(tracer(k))
  end do
  if (is_root_pe()) then
    print *, 'ICE MODEL BUDGET' ! PER EARTH AREA'
    print '(a10,5a22)',   'ICE MODEL ','   AT START  ', &
         ' TOP FLUX DN.', ' BOT FLUX DN.', '   AT END    ', '   ERROR     '
    print '(a10,5es22.14)','WATER(Kg) ', IST%h2o(:), &
                 -(IST%h2o(4) -IST%h2o(1) -IST%h2o(2) +IST%h2o(3))/(IST%h2o(4) +1.0) 
    print '(a10,5es22.14)','HEAT(J)   ', IST%heat(:), &
                                       -(IST%heat(4)-IST%heat(1)-IST%heat(2)+IST%heat(3))/(IST%heat(4)+1.0)
    print '(a10,5es22.14)','SALT(sal) ', IST%salt(:), &
                                       -(IST%salt(4)-IST%salt(1)-IST%salt(2)+IST%salt(3))/(IST%salt(4)+1.0)
!          print '(a10,5es22.14)','TRACER      ', tracer, &
!                            -(tracer(4)-tracer(1)-tracer(2)+tracer(3))/(tracer(4)+1.0)
    print *
  endif
end subroutine ice_print_budget


!#######################################################################
! <SUBROUTINE NAME="ice_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine ice_model_restart(time_stamp)
  character(len=*),         intent(in), optional :: time_stamp

  call save_restart(Ice_restart, time_stamp)
 !call icebergs_save_restart(Ice%icebergs)
 ! This should go here but since "Ice" is not available we have to
 ! rely on the restart written via ice_model_end() -AJA

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
!  integer, dimension(2) :: axt, axv, axtv, axvt
!  integer, dimension(3) :: axt2
  integer               :: id_geo_lon, id_geo_lat, id_sin_rot, id_cos_rot, id_cell_area
  logical               :: sent
!  integer               :: id_xb, id_xt, id_yb, id_yt, id_ct, id_xv, id_yv
  integer :: i, j, k, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  !
  ! diagnostics MUST use a domain without halos otherwise same as the
  ! regular domain:  Domain (see ice_grid.f90)
  !
! id_xv = diag_axis_init('xv', G%gridLonB(G%isg:G%ieg), 'degrees_E', 'X','longitude', set_name='ice', Domain2=Domain )
! id_yv = diag_axis_init('yv', G%gridLatB(G%jsg:G%jeg), 'degrees_N', 'Y','latitude',  set_name='ice', Domain2=Domain )
! id_xb = diag_axis_init('xb', G%gridLonB, 'degrees_E', 'X', 'longitude', set_name='ice', Domain2=Domain )
! id_yb = diag_axis_init('yb', G%gridLatB, 'degrees_N', 'Y', 'latitude', set_name='ice', Domain2=Domain )
! id_xt = diag_axis_init('xt', G%gridLonT, 'degrees_E', 'X', &
!         'longitude',set_name='ice',edges=id_xb,Domain2=Domain)
! id_yt = diag_axis_init('yt', G%gridLatT, 'degrees_N', 'Y', &
!         'latitude',set_name='ice', edges=id_yb,Domain2=Domain)
! id_ct = diag_axis_init('ct', G%H_cat_lim(1:G%CatIce), 'meters','Z', 'thickness')

!  axv  = (/ id_xv, id_yv  /)
!  axt  = (/ id_xt, id_yt  /)
!  axt2 = (/ id_xt, id_yt, id_ct/)
  Ice%axes(:) = diag%axesTc(:)
!  axtv = (/ id_xt, id_yv /); ! for north faces of t-cells
!  axvt = (/ id_xv, id_yt /); ! for east  faces of t-cells

!  G%axesT1(:) = (/ id_xt, id_yt  /)
!  G%axesB1(:) = (/ id_xv, id_yv  /)
!  G%axesCv1(:) = (/ id_xt, id_yv /)
!  G%axesCu1(:) = (/ id_xv, id_yt /)

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
  IST%id_t1       = register_SIS_diag_field('ice_model', 'T1', diag%axesT1, Time, &
               'top ice layer temperature', 'C',  missing_value=missing)
  IST%id_t2       = register_SIS_diag_field('ice_model', 'T2', diag%axesT1, Time, &
               'second ice layer temperature', 'C',  missing_value=missing)
  IST%id_t3       = register_SIS_diag_field('ice_model', 'T3', diag%axesT1, Time, &
               'third ice layer temperature', 'C',  missing_value=missing)
  IST%id_t4       = register_SIS_diag_field('ice_model', 'T4', diag%axesT1, Time, &
               'bottom ice layer temperature', 'C',  missing_value=missing)
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
  IST%id_sw_abs_snow= register_SIS_diag_field('ice_model','sw_abs_snow',diag%axesT1, Time, &
               'SW frac. abs. in snow','0:1', missing_value=missing )
  IST%id_sw_abs_ice1= register_SIS_diag_field('ice_model','sw_abs_ice1',diag%axesT1, Time, &
               'SW frac. abs. in ice1','0:1', missing_value=missing )
  IST%id_sw_abs_ice2= register_SIS_diag_field('ice_model','sw_abs_ice2',diag%axesT1, Time, &
               'SW frac. abs. in ice2','0:1', missing_value=missing )
  IST%id_sw_abs_ice3= register_SIS_diag_field('ice_model','sw_abs_ice3',diag%axesT1, Time, &
               'SW frac. abs. in ice3','0:1', missing_value=missing )
  IST%id_sw_abs_ice4= register_SIS_diag_field('ice_model','sw_abs_ice4',diag%axesT1, Time, &
               'SW frac. abs. in ice4','0:1', missing_value=missing )
  IST%id_sw_pen= register_SIS_diag_field('ice_model','sw_pen',diag%axesT1, Time, &
               'SW frac. pen. surf.','0:1', missing_value=missing )
  IST%id_sw_trn= register_SIS_diag_field('ice_model','sw_trn',diag%axesT1, Time, &
               'SW frac. trans. to ice bot.','0:1', missing_value=missing )


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
  IST%id_fax      = register_SIS_diag_field('ice_model', 'FA_X', diag%axesB1, Time, &
               'air stress on ice - x component', 'Pa', missing_value=missing)
  IST%id_fay      = register_SIS_diag_field('ice_model', 'FA_Y', diag%axesB1, Time, &
               'air stress on ice - y component', 'Pa', missing_value=missing)
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

  if (id_sin_rot>0) call post_data(id_sin_rot, G%sin_rot, diag, is_static=.true.)
  if (id_cos_rot>0) call post_data(id_cos_rot, G%cos_rot, diag, is_static=.true.)
  if (id_geo_lon>0) call post_data(id_geo_lon, G%geoLonT, diag, is_static=.true.)
  if (id_geo_lat>0) call post_data(id_geo_lat, G%geoLatT, diag, is_static=.true.)
  if (id_cell_area>0) call post_data(id_cell_area, cell_area, diag, is_static=.true.)

end subroutine ice_diagnostics_init

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
! write(outunit,100) 'ice_data_type%part_size_uv       ',mpp_chksum(IST%part_size_uv       )
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
! write(outunit,100) 'ice_data_type%sea_lev            ',mpp_chksum(IST%sea_lev            )
  write(outunit,100) 'ice_data_type%s_surf             ',mpp_chksum(Ice%s_surf             )
! write(outunit,100) 'ice_data_type%u_ocn              ',mpp_chksum(IST%u_ocn              )
! write(outunit,100) 'ice_data_type%v_ocn              ',mpp_chksum(IST%v_ocn              )
! write(outunit,100) 'ice_data_type%flux_u_top         ',mpp_chksum(IST%flux_u_top         )
! write(outunit,100) 'ice_data_type%flux_v_top         ',mpp_chksum(IST%flux_v_top         )
! write(outunit,100) 'ice_data_type%flux_t_top         ',mpp_chksum(IST%flux_t_top         )
! write(outunit,100) 'ice_data_type%flux_q_top         ',mpp_chksum(IST%flux_q_top         )
! write(outunit,100) 'ice_data_type%flux_lw_top        ',mpp_chksum(IST%flux_lw_top        )
! write(outunit,100) 'ice_data_type%flux_sw_vis_dir_top',mpp_chksum(IST%flux_sw_vis_dir_top)
! write(outunit,100) 'ice_data_type%flux_sw_vis_dif_top',mpp_chksum(IST%flux_sw_vis_dif_top)
! write(outunit,100) 'ice_data_type%flux_sw_nir_dir_top',mpp_chksum(IST%flux_sw_nir_dir_top)
! write(outunit,100) 'ice_data_type%flux_sw_nir_dif_top',mpp_chksum(IST%flux_sw_nir_dif_top)
! write(outunit,100) 'ice_data_type%flux_lh_top        ',mpp_chksum(IST%flux_lh_top        )
! write(outunit,100) 'ice_data_type%lprec_top          ',mpp_chksum(IST%lprec_top          )
! write(outunit,100) 'ice_data_type%fprec_top          ',mpp_chksum(IST%fprec_top          )
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
!  write(outunit,100) 'ice_data_type%lwdn               ',mpp_chksum(IST%lwdn               )
!  write(outunit,100) 'ice_data_type%swdn               ',mpp_chksum(IST%swdn               )
!  write(outunit,100) 'ice_data_type%pen                ',mpp_chksum(IST%pen                )
!  write(outunit,100) 'ice_data_type%trn                ',mpp_chksum(IST%trn                )
!  write(outunit,100) 'ice_data_type%tmelt              ',mpp_chksum(IST%tmelt              )
!  write(outunit,100) 'ice_data_type%bmelt              ',mpp_chksum(IST%bmelt              )
!  write(outunit,100) 'ice_data_type%h_snow             ',mpp_chksum(IST%h_snow             )
!  write(outunit,100) 'ice_data_type%t_snow             ',mpp_chksum(IST%t_snow             )
!  write(outunit,100) 'ice_data_type%h_ice              ',mpp_chksum(IST%h_ice              )
!  write(outunit,100) 'ice_data_type%t_ice(1)           ',mpp_chksum(IST%t_ice(:,:,:,1)     )
!  write(outunit,100) 'ice_data_type%t_ice(2)           ',mpp_chksum(IST%t_ice(:,:,:,2)     )
!  write(outunit,100) 'ice_data_type%t_ice(3)           ',mpp_chksum(IST%t_ice(:,:,:,3)     )
!  write(outunit,100) 'ice_data_type%t_ice(4)           ',mpp_chksum(IST%t_ice(:,:,:,4)     )
!  write(outunit,100) 'ice_data_type%u_ice              ',mpp_chksum(IST%u_ice              )
!  write(outunit,100) 'ice_data_type%v_ice              ',mpp_chksum(IST%v_ice              )
!  write(outunit,100) 'ice_data_type%frazil             ',mpp_chksum(IST%frazil)
!  write(outunit,100) 'ice_data_type%bheat              ',mpp_chksum(IST%bheat)

  do n=1,Ice%ocean_fields%num_bcs ; do m=1,Ice%ocean_fields%bc(n)%num_fields
    write(outunit,101) 'ice%', trim(Ice%ocean_fields%bc(n)%name), &
                       trim(Ice%ocean_fields%bc(n)%field(m)%name), &
                       mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
  enddo ; enddo

100 FORMAT("   CHECKSUM::",A32," = ",Z20)
101 FORMAT("   CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ice_data_type_chksum

end module ice_type_mod
