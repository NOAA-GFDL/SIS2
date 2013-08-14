!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_type_mod - maintains the sea ice data, reads/writes restarts, reads the  !
!                namelist and initializes diagnostics. - Mike Winton           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_type_mod

  use mpp_mod,          only: mpp_pe, mpp_root_pe, mpp_sum, mpp_clock_id, CLOCK_COMPONENT, &
                              CLOCK_LOOP, CLOCK_ROUTINE, stdout,input_nml_file
  use mpp_domains_mod,  only: domain2D, mpp_update_domains, CORNER, BGRID_NE
  use fms_mod,          only: file_exist, open_namelist_file, check_nml_error, write_version_number,&
                              read_data, close_file, field_exist, &
                              stderr, stdlog, error_mesg, FATAL, WARNING, NOTE, clock_flag_default
  use fms_io_mod,       only: save_restart, restore_state, query_initialized, &
                              register_restart_field, restart_file_type, set_domain, nullify_domain, &
                              parse_mask_table
  use diag_manager_mod, only: diag_axis_init, register_diag_field, &
                              register_static_field, send_data
  use time_manager_mod, only: time_type, get_time
  use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type
  use constants_mod,    only: Tfreeze, radius, pi
  use ice_grid_mod,     only: set_ice_grid, t_to_uv, ice_grid_end
  use ice_grid_mod,     only: sea_ice_grid_type
  use ice_grid_mod,     only: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, im, jm, km
  use ice_grid_mod,     only: cell_area, sin_rot, cos_rot, xb1d, yb1d
  use ice_grid_mod,     only: grid_x_t,grid_y_t
  use ice_grid_mod,     only: x_cyclic, tripolar_grid
  use ice_thm_mod,      only: ice_thm_param, DI, DS, e_to_melt
use ice_dyn_mod,       only: ice_dyn_init, ice_dyn_CS, ice_dyn_register_restarts, ice_dyn_end
use ice_transport_mod, only: ice_transport_init, ice_transport_CS, ice_transport_end ! , ice_transport_register_restarts
  use constants_mod,    only: LI => hlf ! latent heat of fusion - 334e3 J/(kg-ice)
  use ice_bergs,        only: icebergs_init, icebergs_end, icebergs, icebergs_stock_pe
  use ice_bergs,        only: icebergs_save_restart
  use astronomy_mod,    only: astronomy_init, astronomy_end
  use ice_shortwave_dEdd,only: shortwave_dEdd0_set_params

use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : open_param_file, close_param_file
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use SIS_diag_mediator, only: SIS_diag_ctrl, set_SIS_axes_info, SIS_diag_mediator_init
use SIS_get_input, only : Get_SIS_input, directories, archaic_nml_check

implicit none ; private

public :: ice_data_type, ice_state_type, ice_model_init, ice_model_end, ice_stock_pe,  &
          mom_rough_ice, heat_rough_ice, atmos_winds, hlim, slab_ice, kmelt,  &
          spec_ice, verbose, ice_bulk_salin, do_ice_restore, do_ice_limit,    &
          max_ice_limit, ice_restore_timescale, do_init, h2o, heat, salt, slp2ocean,&
          conservation_check, do_icebergs, ice_model_restart,       &
          add_diurnal_sw, ice_data_type_chksum
public :: do_sun_angle_for_alb

public  :: id_cn, id_hi, id_hs, id_tsn, id_t1, id_t2, id_t3, id_t4, id_ts,id_hio
public  :: id_mi, id_sh, id_lh, id_sw, id_lw, id_snofl, id_rain, id_runoff,    &
           id_calving, id_runoff_hflx, id_calving_hflx,                        &
           id_evap, id_saltf, id_tmelt, id_bmelt, id_bheat, id_e2m,            &
           id_frazil, id_alb, id_xprt, id_lsrc, id_lsnk, id_bsnk, id_strna,    &
           id_fax, id_fay, id_swdn, id_lwdn, id_sn2ic,                         &
           id_slp, id_ext, id_sst, id_sss, id_ssh, id_uo, id_vo, id_ta, id_obi,&
           id_qfres, id_qflim, id_ix_trans, id_iy_trans,                       &
           id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir, id_sw_vis_dif,      &
           id_sw_nir_dir, id_sw_nir_dif, id_mib

public  :: id_alb_vis_dir, id_alb_vis_dif,id_alb_nir_dir, id_alb_nir_dif, id_coszen
public  :: id_abs_int,id_sw_abs_snow,id_sw_abs_ice1,id_sw_abs_ice2,id_sw_abs_ice3,id_sw_abs_ice4,id_sw_pen,id_sw_trn
public  :: iceClock,iceClock1,iceClock2,iceClock3,iceClock4,iceClock5,iceClock6,iceClock7,iceClock8,iceClock9
public  :: iceClocka,iceClockb,iceClockc
public  :: earth_area

  real, parameter :: earth_area = 4*PI*RADIUS*RADIUS !5.10064471909788E+14 m^2
  real, parameter :: missing = -1e34
  integer, parameter :: miss_int = -9999
  !---- id for diagnositics -------------------
  integer :: id_xb, id_xt, id_yb, id_yt, id_ct, id_xv, id_yv
  integer :: id_cn, id_hi, id_hs, id_tsn, id_t1, id_t2, id_t3, id_t4, id_ts,id_hio
  integer :: id_mi, id_sh, id_lh, id_sw, id_lw, id_snofl, id_rain
  integer :: id_runoff, id_calving, id_runoff_hflx, id_calving_hflx
  integer :: id_evap, id_saltf, id_tmelt, id_bmelt, id_bheat, id_e2m
  integer :: id_frazil, id_alb, id_xprt, id_lsrc, id_lsnk, id_bsnk, id_strna
  integer :: id_sigi, id_sigii, id_stren, id_ui, id_vi, id_fax, id_fay, id_fix
  integer :: id_fiy, id_fcx, id_fcy, id_fwx, id_fwy, id_swdn, id_lwdn, id_sn2ic
  integer :: id_slp, id_ext, id_sst, id_sss, id_ssh, id_uo, id_vo, id_ta, id_obi
  integer :: id_qfres, id_qflim, id_ix_trans, id_iy_trans
  integer :: id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir, id_sw_vis_dif
  integer :: id_sw_nir_dir, id_sw_nir_dif, id_mib
  integer :: id_alb_vis_dir, id_alb_vis_dif, id_alb_nir_dir, id_alb_nir_dif, id_coszen 
  integer :: id_abs_int,id_sw_abs_snow,id_sw_abs_ice1,id_sw_abs_ice2,id_sw_abs_ice3,id_sw_abs_ice4,id_sw_pen,id_sw_trn

  character(len=128) :: version = '$Id: ice_type.F90,v 1.1.2.1.6.1.2.2.2.1 2013/06/18 22:24:14 nnz Exp $'
  character(len=128) :: tagname = '$Name: siena_201305_ice_sis2_5layer_dEdd_nnz $'

  !--- namelist interface --------------
  real    :: mom_rough_ice  = 1.0e-4     ! momentum same, cd10=(von_k/ln(10/z0))^2
  real    :: heat_rough_ice = 1.0e-4     ! heat roughness length
  real    :: kmelt          = 6e-5*4e6   ! ocean/ice heat flux constant
  real    :: ks             = 0.31       ! snow conductivity (W/mK)
  real    :: alb_sno        = 0.85       ! snow albedo (less if melting)
  real    :: alb_ice        = 0.5826     ! ice albedo (less if melting)
  real    :: pen_ice        = 0.3        ! part unreflected solar penetrates ice
  real    :: opt_dep_ice    = 0.67       ! ice optical depth
  real    :: t_range_melt   = 1.0        ! melt albedos scaled in over T range
  real    :: ice_bulk_salin = 0.004      ! ice bulk salinity (for ocean salt flux)!CICE value
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

  logical :: do_init = .false.
  real    :: hlim(8) = (/ 0.0, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! thickness limits 1...num_part-1
  real    :: h2o(4), heat(4), salt(4) ! for conservation analysis
                             ! 1 - initial ice h2o/heat content
                             ! 2 - h2o/heat flux down at top of ice
                             ! 3 - h2o/heat flux down at bottom of ice
                             ! 4 - final ice h2o/heat content

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
  real,    pointer, dimension(:,:,:) :: part_size           =>NULL()
  real,    pointer, dimension(:,:,:) :: part_size_uv        =>NULL()
  real,    pointer, dimension(:,:)   :: sea_lev             =>NULL()

  real,    pointer, dimension(:,:,:)   :: t_surf            =>NULL()
  real,    pointer, dimension(:,:)   :: s_surf              =>NULL()
  real,    pointer, dimension(:,:)   :: u_ocn               =>NULL()
  real,    pointer, dimension(:,:)   :: v_ocn               =>NULL()
  
  real,    pointer, dimension(:,:,:) :: flux_u_top          =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_v_top          =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_u_top_bgrid    =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_v_top_bgrid    =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_t_top          =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_q_top          =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_lw_top         =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_sw_vis_dir_top =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_sw_vis_dif_top =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_sw_nir_dir_top =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_sw_nir_dif_top =>NULL()
  real,    pointer, dimension(:,:,:) :: flux_lh_top         =>NULL()
  real,    pointer, dimension(:,:,:) :: lprec_top           =>NULL()
  real,    pointer, dimension(:,:,:) :: fprec_top           =>NULL()

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
  real,    pointer, dimension(:,:,:) :: h_snow              =>NULL()
  real,    pointer, dimension(:,:,:) :: t_snow              =>NULL()
  real,    pointer, dimension(:,:,:) :: h_ice               =>NULL()
  real,    pointer, dimension(:,:,:,:) :: t_ice             =>NULL()
  real,    pointer, dimension(:,:)   :: u_ice               =>NULL()
  real,    pointer, dimension(:,:)   :: v_ice               =>NULL()
  real,    pointer, dimension(:,:)   :: frazil              =>NULL()
  real,    pointer, dimension(:,:)   :: bheat               =>NULL()
  real,    pointer, dimension(:,:)   :: qflx_lim_ice        =>NULL()
  real,    pointer, dimension(:,:)   :: qflx_res_ice        =>NULL()
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
!   integer, dimension(3)              :: axes
   type(coupler_3d_bc_type)           :: ocean_fields       ! array of fields used for additional tracers
   type(coupler_2d_bc_type)           :: ocean_fluxes       ! array of fluxes used for additional tracers
   type(coupler_3d_bc_type)           :: ocean_fluxes_top   ! array of fluxes for averaging

  type(ice_dyn_CS), pointer       :: ice_dyn_CSp => NULL()
  type(ice_transport_CS), pointer :: ice_transport_CSp => NULL()
  type(SIS_diag_ctrl)         :: diag
  type(icebergs), pointer     :: icebergs => NULL()
  type(sea_ice_grid_type) :: G ! A structure containing metrics and grid info.
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
    type(sea_ice_grid_type) :: G ! A structure containing metrics and grid info.
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

  use ice_grid_mod, only : all_avg

  type(ice_data_type)    :: Ice
  integer, intent(in) :: index
  real, intent(out)   :: value
  type(ice_state_type), pointer :: IST => NULL()

  integer :: i, j, k
  real :: icebergs_value

  value = 0.0
  if(.not.Ice%pe) return

  IST => Ice%Ice_state

  select case (index)

    case (ISTOCK_WATER)

      value = sum(cell_area(:,:)*all_avg(DI*IST%h_ice(isc:iec,jsc:jec,:)+   &
          DS*IST%h_snow(isc:iec,jsc:jec,:),  &
          IST%part_size(isc:iec,jsc:jec,:))) &
          *4*pi*radius*radius 

    case (ISTOCK_HEAT)

      value = 0.0
      do k=2,km ; do j=jsc,jec ; do i=isc,iec
        if ((IST%part_size(i,j,k)>0.0.and.IST%h_ice(i,j,k)>0.0)) then
          if (slab_ice) then
            value = value - cell_area(i,j) * IST%part_size(i,j,k)*IST%h_ice(i,j,2)*DI*LI
          else
            value = value - cell_area(i,j) * IST%part_size(i,j,k)           &
                           *e_to_melt(IST%h_snow(i,j,k), IST%t_snow(i,j,k), &
                                      IST%h_ice(i,j,k), IST%t_ice(i,j,k,1),  &
                                      IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3), &
                                      IST%t_ice(i,j,k,4) )
          endif
        endif
      enddo ; enddo ; enddo
      value = value*4*pi*radius*radius

    case (ISTOCK_SALT)
       !No salt in the h_snow component.
      value =  sum(cell_area(:,:)*all_avg(DI*IST%h_ice(isc:iec,jsc:jec,:),IST%part_size(isc:iec,jsc:jec,:))) &
              *ice_bulk_salin*4*pi*radius*radius
    case default

      value = 0.0

  end select

  if (do_icebergs) then
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

    integer           :: io, ierr, nlon, nlat, npart, unit, log_unit, k
    integer           :: sc, dy, i, j, l
    integer           :: id_restart, id_restart_albedo, id_restart_flux_sw
    real              :: dt_slow
    character(len=64) :: restart_file
    integer           :: stdlogunit, stdoutunit
    type(param_file_type) :: param_file
  type(ice_state_type), pointer :: IST => NULL()

    stdlogunit=stdlog()
    stdoutunit = stdout()

  if (associated(Ice%Ice_state)) then
    call SIS_error(WARNING, "ice_model_init called with an associated "// &
                    "Ice%Ice_state structure. Model is already initialized.")
    return
  endif
  allocate(Ice%Ice_state)
  IST => Ice%Ice_state
  

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
  call archaic_nml_check(param_file, "NSTEPS_DYN", "nsteps_dyn", nsteps_dyn, miss_int, 432)
  call archaic_nml_check(param_file, "ICE_STRENGTH_PSTAR", "p0", p0, missing)
  call archaic_nml_check(param_file, "ICE_STRENGTH_CSTAR", "c0", c0, missing)
  call archaic_nml_check(param_file, "ICE_CDRAG_WATER", "cdw", cdw, missing)
  ! call archaic_nml_check(param_file, "USE_SLAB_ICE", "SLAB_ICE", slab_ice, .false.)
  call archaic_nml_check(param_file, "AIR_WATER_STRESS_TURN_ANGLE", "wd_turn", wd_turn, missing)
  ! call archaic_nml_check(param_file, "SPECIFIED_ICE", "spec_ice", spec_ice, .false.)
  call archaic_nml_check(param_file, "NSTEPS_ADV", "nsteps_adv", nsteps_adv, miss_int, 1)
  call archaic_nml_check(param_file, "ICE_CHANNEL_VISCOSITY", &
                         "channel_viscosity", channel_viscosity, missing, 0.0)
  call archaic_nml_check(param_file, "ICE_CHANNEL_SMAG_COEF", "smag_ocn", smag_ocn, missing)
  call archaic_nml_check(param_file, "ICE_CHANNEL_CFL_LIMIT", "chan_cfl_limit", chan_cfl_limit, missing)

    if (spec_ice) then
       slab_ice = .true.
       nsteps_dyn = 0
       nsteps_adv = 0
    end if
    if (slab_ice) num_part = 2 ! open water and ice ... but never in same place
    if (num_part>size(hlim(:))+1) &
         call error_mesg ('ice_model_init', 'not enough thickness limits', FATAL)

    call get_time(Time_step_slow, sc, dy); dt_slow=864e2*dy+sc

  if (file_exist(mask_table)) then
     write(stdoutunit, *) '==> NOTE from ice_model_init:  reading maskmap information from '//trim(mask_table)
     if(layout(1) == 0 .OR. layout(2) == 0 ) call error_mesg('ice_model_init', &
        'ice_model_nml layout should be set when file '//trim(mask_table)//' exists', FATAL)

     allocate(Ice%maskmap(layout(1), layout(2)))
     call parse_mask_table(mask_table, Ice%maskmap, "Ice model")
  endif

  if( ASSOCIATED(Ice%maskmap) ) then
     call set_ice_grid(Ice%G, param_file, Ice%domain, num_part, layout, io_layout, Ice%maskmap  )
  else
     call set_ice_grid(Ice%G, param_file, Ice%domain, num_part, layout, io_layout )
  end if
  call set_domain(domain)

  allocate(Ice%mask(isc:iec, jsc:jec)) ; Ice%mask(:,:) = .false. !derived
  allocate(Ice%ice_mask(isc:iec, jsc:jec, km)) ; Ice%ice_mask(:,:,:) = .false. !NI
  allocate(Ice%t_surf(isc:iec, jsc:jec, km)) ; Ice%t_surf(:,:,:) = 0.0
  allocate(Ice%s_surf(isc:iec, jsc:jec)) ; Ice%s_surf(:,:) = 0.0 !NI
  allocate(IST%t_surf(isc:iec, jsc:jec, km)) ; IST%t_surf(:,:,:) = 0.0
  allocate(IST%s_surf(isc:iec, jsc:jec)) ; IST%s_surf(:,:) = 0.0 !NI
  allocate(IST%sea_lev(isd:ied, jsd:jed)) ; IST%sea_lev(:,:) = 0.0 !NR
  allocate(Ice%part_size(isd:ied, jsd:jed, km)) ; Ice%part_size(:,:,:) = 0.0
  allocate(IST%part_size(isd:ied, jsd:jed, km)) ; IST%part_size(:,:,:) = 0.0
  allocate(IST%part_size_uv(isc:iec, jsc:jec, km)) ; IST%part_size_uv(:,:,:) = 0.0 !NR
  allocate(Ice%u_surf(isc:iec, jsc:jec, km)) ; Ice%u_surf(:,:,:) = 0.0 !NI
  allocate(Ice%v_surf(isc:iec, jsc:jec, km)) ; Ice%v_surf(:,:,:) = 0.0 !NI
  allocate(IST%u_ocn(isd:ied, jsd:jed)) ; IST%u_ocn(:,:) = 0.0 !NR
  allocate(IST%v_ocn(isd:ied, jsd:jed)) ; IST%v_ocn(:,:) = 0.0 !NR
  allocate(Ice%rough_mom(isc:iec, jsc:jec, km)) ; Ice%rough_mom(:,:,:) = 0.0
  allocate(Ice%rough_heat(isc:iec, jsc:jec, km)) ; Ice%rough_heat(:,:,:) = 0.0
  allocate(Ice%rough_moist(isc:iec, jsc:jec, km)) ; Ice%rough_moist(:,:,:) = 0.0
  allocate(IST%coszen(isc:iec, jsc:jec)) ; IST%coszen(:,:) = 0.0 !NR
  allocate(Ice%albedo(isc:iec, jsc:jec, km)) ; Ice%albedo(:,:,:) = 0.0
  allocate(Ice%albedo_vis_dir(isc:iec, jsc:jec, km)) ; Ice%albedo_vis_dir(:,:,:) = 0.0
  allocate(Ice%albedo_nir_dir(isc:iec, jsc:jec, km)) ; Ice%albedo_nir_dir(:,:,:) = 0.0
  allocate(Ice%albedo_vis_dif(isc:iec, jsc:jec, km)) ; Ice%albedo_vis_dif(:,:,:) = 0.0
  allocate(Ice%albedo_nir_dif(isc:iec, jsc:jec, km)) ; Ice%albedo_nir_dif(:,:,:) = 0.0

  allocate(IST%flux_u_top(isd:ied, jsd:jed, km)) ; IST%flux_u_top(:,:,:) = 0.0 !NR
  allocate(IST%flux_v_top(isd:ied, jsd:jed, km)) ; IST%flux_v_top(:,:,:) = 0.0 !NR 
  allocate(IST%flux_t_top(isc:iec, jsc:jec, km)) ;  IST%flux_t_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_q_top(isc:iec, jsc:jec, km)) ;  IST%flux_q_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_sw_vis_dir_top(isc:iec, jsc:jec, km)) ; IST%flux_sw_vis_dir_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_sw_vis_dif_top(isc:iec, jsc:jec, km)) ; IST%flux_sw_vis_dif_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_sw_nir_dir_top(isc:iec, jsc:jec, km)) ; IST%flux_sw_nir_dir_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_sw_nir_dif_top(isc:iec, jsc:jec, km)) ; IST%flux_sw_nir_dif_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_lw_top(isc:iec, jsc:jec, km)) ; IST%flux_lw_top(:,:,:) = 0.0 !NI
  allocate(IST%flux_lh_top(isc:iec, jsc:jec, km)) ; IST%flux_lh_top(:,:,:) = 0.0 !NI
  allocate(IST%lprec_top(isc:iec, jsc:jec, km)) ;  IST%lprec_top(:,:,:) = 0.0 !NI
  allocate(IST%fprec_top(isc:iec, jsc:jec, km)) ;  IST%fprec_top(:,:,:) = 0.0 !NI

  allocate(IST%flux_u_top_bgrid(isd:ied, jsd:jed, km)) ; IST%flux_u_top_bgrid(:,:,:) = 0.0 !NR
  allocate(IST%flux_v_top_bgrid(isd:ied, jsd:jed, km)) ; IST%flux_v_top_bgrid(:,:,:) = 0.0 !NR

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
  allocate(IST%lwdn(isc:iec, jsc:jec)) ; IST%lwdn(:,:) = 0.0 !NR
  allocate(IST%swdn(isc:iec, jsc:jec)) ; IST%swdn(:,:) = 0.0 !NR
  allocate(IST%frazil(isc:iec, jsc:jec)) ; IST%frazil(:,:) = 0.0 !NR
  allocate(IST%bheat(isc:iec, jsc:jec)) ; IST%bheat(:,:) = 0.0 !NI
  allocate(IST%u_ice(isd:ied, jsd:jed)) ; IST%u_ice(:,:) = 0.0
  allocate(IST%v_ice(isd:ied, jsd:jed)) ; IST%v_ice(:,:) = 0.0
  allocate(IST%tmelt(isc:iec, jsc:jec, 2:km)) ; IST%tmelt(:,:,:) = 0.0 !NR
  allocate(IST%bmelt(isc:iec, jsc:jec, 2:km)) ; IST%bmelt(:,:,:) = 0.0 !NR
  allocate(IST%pen(isc:iec, jsc:jec, 2:km)) ; IST%pen(:,:,:) = 0.0 !NI
  allocate(IST%trn(isc:iec, jsc:jec, 2:km)) ; IST%trn(:,:,:) = 0.0 !NI
  allocate(IST%sw_abs_sfc(isc:iec, jsc:jec, 2:km)) ; IST%sw_abs_sfc(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_snow(isc:iec, jsc:jec, 2:km)) ; IST%sw_abs_snow(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_ice(isc:iec, jsc:jec, 2:km, Ice%G%NkIce)) ; IST%sw_abs_ice(:,:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_ocn(isc:iec, jsc:jec, 2:km)) ; IST%sw_abs_ocn(:,:,:) = 0.0 !NR
  allocate(IST%sw_abs_int(isc:iec, jsc:jec, 2:km)) ; IST%sw_abs_int(:,:,:) = 0.0 !NR
  allocate(IST%h_snow(isd:ied, jsd:jed, 2:km)) ; IST%h_snow(:,:,:) = 0.0
  allocate(IST%t_snow(isd:ied, jsd:jed, 2:km)) ; IST%t_snow(:,:,:) = 0.0
  allocate(IST%h_ice(isd:ied, jsd:jed, 2:km)) ; IST%h_ice(:,:,:) = 0.0
  allocate(IST%t_ice(isd:ied, jsd:jed, 2:km, Ice%G%NkIce)) ; IST%t_ice(:,:,:,:) = 0.0
  allocate(IST%qflx_lim_ice(isc:iec, jsc:jec)) ; IST%qflx_lim_ice(:,:) = 0.0 !NR
  allocate(IST%qflx_res_ice(isc:iec, jsc:jec)) ; IST%qflx_res_ice(:,:) = 0.0 !NR

  allocate(Ice%area(isc:iec, jsc:jec)) ; Ice%area(:,:) = 0.0 !derived
  allocate(Ice%mi(isc:iec, jsc:jec)) ; Ice%mi(:,:) = 0.0 !NR

  Ice%area(:,:)       = cell_area(:,:) * 4*PI*RADIUS*RADIUS
  IST%coszen(:,:) = cos(3.14*67.0/180.0) ! NP summer solstice.

  do j=jsc,jec ; do i=isc,iec
    Ice%mask(i,j) = ( Ice%G%mask2dT(i,j) > 0.5 )
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
  restart_file = 'ice_model.res.nc'
  id_restart = register_restart_field(Ice_restart, restart_file, 'part_size', IST%part_size, domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'albedo',    Ice%albedo,    domain=domain)
  id_restart_albedo = register_restart_field(Ice_restart, restart_file, 'albedo_vis_dir', Ice%albedo_vis_dir, &
                                             domain=domain, mandatory=.false.)
  id_restart        = register_restart_field(Ice_restart, restart_file, 'albedo_nir_dir', Ice%albedo_nir_dir, &
                                             domain=domain, mandatory=.false.)
  id_restart        = register_restart_field(Ice_restart, restart_file, 'albedo_vis_dif', Ice%albedo_vis_dif, &
                                             domain=domain, mandatory=.false.)
  id_restart        = register_restart_field(Ice_restart, restart_file, 'albedo_nir_dif', Ice%albedo_nir_dif, &
                                             domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'rough_mom',   Ice%rough_mom,        domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'rough_heat',  Ice%rough_heat,       domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'rough_moist', Ice%rough_moist,      domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_surf',      IST%t_surf,           domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'h_snow',      IST%h_snow(:,:,2:km), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_snow',      IST%t_snow(:,:,2:km), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'h_ice',       IST%h_ice(:,:,2:km),  domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_ice1',      IST%t_ice(:,:,2:km,1), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_ice2',      IST%t_ice(:,:,2:km,2), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_ice3',      IST%t_ice(:,:,2:km,3), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 't_ice4',      IST%t_ice(:,:,2:km,4), domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'u_ice',       IST%u_ice,            domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'v_ice',       IST%v_ice,            domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_u',      Ice%flux_u,           domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_v',      Ice%flux_v,           domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_t',      Ice%flux_t,           domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_q',      Ice%flux_q,           domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_salt',   Ice%flux_salt,        domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'flux_lw',     Ice%flux_lw,          domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'lprec',       Ice%lprec,            domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'fprec',       Ice%fprec,            domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'runoff',      Ice%runoff,           domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'calving',     Ice%calving,          domain=domain)
  id_restart = register_restart_field(Ice_restart, restart_file, 'runoff_hflx', Ice%runoff_hflx,      domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'calving_hflx',Ice%calving_hflx,     domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'p_surf',      Ice%p_surf,           domain=domain)
  id_restart         = register_restart_field(Ice_restart, restart_file, 'flux_sw_vis_dir', Ice%flux_sw_vis_dir, &
                                              domain=domain, mandatory=.false.)    
  id_restart         = register_restart_field(Ice_restart, restart_file, 'flux_sw_vis_dif', Ice%flux_sw_vis_dif, &
                                              domain=domain, mandatory=.false.)
  id_restart_flux_sw = register_restart_field(Ice_restart, restart_file, 'flux_sw_nir_dir', Ice%flux_sw_nir_dir, &
                                              domain=domain, mandatory=.false.)
  id_restart         = register_restart_field(Ice_restart, restart_file, 'flux_sw_nir_dif', Ice%flux_sw_nir_dif, &
                                              domain=domain, mandatory=.false.)
  id_restart = register_restart_field(Ice_restart, restart_file, 'coszen',    IST%coszen,    domain=domain, mandatory=.false.)

  call ice_dyn_register_restarts(Ice%G, param_file, IST%ice_dyn_CSp, Ice_restart, restart_file)
!    call ice_transport_register_restarts(Ice%G, param_file, IST%ice_transport_CSp, Ice_restart, restart_file)

  restart_file = 'INPUT/ice_model.res.nc'
  if (file_exist(restart_file)) then
    call restore_state(Ice_restart)

    if ( .NOT.query_initialized(Ice_restart, id_restart_flux_sw) ) then
       call error_mesg ('ice_model_init', &
          'Restart file does not contain flux_sw_* subcomponents!', WARNING)
       ! for compatibility with preN restarts, we should check for total SW flux 
       if (field_exist(restart_file,'flux_sw')) then
         call read_data( restart_file, 'flux_sw', Ice%flux_sw_vis_dir, domain) 
         ! simplest way to break the total flux to 4 components
         Ice%flux_sw_vis_dir(:,:) = Ice%flux_sw_vis_dir(:,:) / 4
         Ice%flux_sw_vis_dif(:,:) = Ice%flux_sw_vis_dir(:,:)
         Ice%flux_sw_nir_dir(:,:) = Ice%flux_sw_vis_dir(:,:)
         Ice%flux_sw_nir_dif(:,:) = Ice%flux_sw_vis_dir(:,:)
       else
         call error_mesg ('ice_model_init', &
              'Restart file does not contain flux_sw total or its components!', FATAL)
       endif
     endif

     !--- update to data domain
     call mpp_update_domains(IST%part_size, Domain)
     call mpp_update_domains(IST%h_snow(:,:,2:km), Domain )
     call mpp_update_domains(IST%t_snow(:,:,2:km), Domain )
     call mpp_update_domains(IST%h_ice (:,:,2:km), Domain )

     do l=1,Ice%G%NkIce
       call mpp_update_domains(IST%t_ice(:,:,2:km,l), Domain )
     enddo

    call mpp_update_domains(IST%u_ice, IST%v_ice, Domain, gridtype=BGRID_NE )
  else ! no restart => no ice
    IST%part_size(:,:,:) = 0.0
    IST%part_size(:,:,1) = 1.0

    Ice%rough_mom(:,:,:)   = mom_rough_ice
    Ice%rough_heat(:,:,:)  = heat_rough_ice
    Ice%rough_moist(:,:,:) = heat_rough_ice
    IST%t_surf(:,:,:) = Tfreeze-5.0
    IST%t_snow(:,:,:) = -5.0
    IST%t_ice(:,:,:,:) = -5.0

    do_init = .true. ! done in ice_model
  endif ! file_exist(restart_file)

  Ice%part_size(:,:,:) = IST%part_size(:,:,:)
  Ice%t_surf(:,:,:) = IST%t_surf(:,:,:)

  IST%part_size_uv(:,:,1) = 1.0
  do k=2,km
    IST%part_size_uv(:,:,k) = 0.0
    call t_to_uv(IST%part_size(:,:,k), IST%part_size_uv(:,:,k), Ice%G)
    IST%part_size_uv (:,:,1) = IST%part_size_uv(:,:,1)-IST%part_size_uv (:,:,k)
  enddo

    
  call SIS_diag_mediator_init(Ice%G, param_file, IST%diag, component="SIS")
  call set_SIS_axes_info(Ice%G, param_file, IST%diag)

  call ice_diagnostics_init(Ice, Ice%G, IST%Time)

  call ice_dyn_init(IST%Time, Ice%G, param_file, IST%diag, IST%ice_dyn_CSp)
  call ice_transport_init(IST%Time, Ice%G, param_file, IST%diag, IST%ice_transport_CSp)
  call ice_thm_param(alb_sno, alb_ice, pen_ice, opt_dep_ice, slab_ice, &
                     t_range_melt, ks, h_lo_lim,do_deltaEdd)

  call close_param_file(param_file)

  !Balaji
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
  if (do_icebergs) call icebergs_init(Ice%icebergs, &
           im, jm, layout, io_layout, Ice%axes(1:2), Ice%maskmap, x_cyclic, tripolar_grid, &
           dt_slow, Time, Ice%G%geoLonBu(isc:iec,jsc:jec), Ice%G%geoLatBu(isc:iec,jsc:jec), &
           Ice%G%mask2dT, Ice%G%dxCv, Ice%G%dyCu, cell_area, cos_rot, sin_rot )

  if (add_diurnal_sw .or. do_sun_angle_for_alb) call astronomy_init

  call shortwave_dEdd0_set_params(R_ice,R_snw,R_pnd)

  call nullify_domain()

end subroutine ice_model_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_end - writes the restart file and deallocates memory               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_end (Ice)
  type (ice_data_type), intent(inout) :: Ice
  integer           :: k

  integer           :: unit
  character(len=22) :: restart='RESTART/ice_model.res'
  type(ice_state_type), pointer :: IST => NULL()

  if (conservation_check) call ice_print_budget()

  call ice_model_restart()
  IST => Ice%Ice_state

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
  deallocate(IST%qflx_lim_ice, IST%qflx_res_ice )
  deallocate(Ice%flux_sw_vis_dir, Ice%flux_sw_vis_dif )
  deallocate(Ice%flux_sw_nir_dir, Ice%flux_sw_nir_dif )

  call ice_dyn_end(IST%ice_dyn_CSp)
  call ice_transport_end(IST%ice_transport_CSp)

  ! End icebergs
  if (do_icebergs) call icebergs_end(Ice%icebergs)
  
  deallocate(Ice%Ice_state)

  if (add_diurnal_sw .or. do_sun_angle_for_alb) call astronomy_end

end subroutine ice_model_end

    subroutine ice_print_budget
      integer :: k
       do k=1,4
          call mpp_sum(h2o(k))
          call mpp_sum(heat(k))
          call mpp_sum(salt(k))
!          call mpp_sum(tracer(k))
       end do
       if (mpp_pe()==mpp_root_pe()) then
          print *, 'ICE MODEL BUDGET' ! PER EARTH AREA'
          print '(a10,5a22)',   'ICE MODEL ','   AT START  ', &
               ' TOP FLUX DN.', &
               ' BOT FLUX DN.', &
               '   AT END    ', &
               '   ERROR     '
          print '(a10,5es22.14)','WATER(Kg) ', h2o * earth_area , &
                                             -(h2o(4) -h2o(1) -h2o(2) +h2o(3))/(h2o(4) +1.0/earth_area) 
          print '(a10,5es22.14)','HEAT(J)   ', heat* earth_area , &
                                             -(heat(4)-heat(1)-heat(2)+heat(3))/(heat(4)+1.0/earth_area)
          print '(a10,5es22.14)','SALT(sal) ', salt * earth_area, &
                                             -(salt(4)-salt(1)-salt(2)+salt(3))/(salt(4)+1.0/earth_area)
!          print '(a10,5es22.14)','TRACER      ', tracer * earth_area, &
!                                             -(tracer(4)-tracer(1)-tracer(2)+tracer(3))/(tracer(4)+1.0/earth_area)
          print *
       end if
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

subroutine ice_diagnostics_init(Ice, G, Time)
  type(ice_data_type),     intent(inout)    :: Ice
  type(sea_ice_grid_type), intent(inout) :: G
  type(time_type),         intent(inout) :: Time

  real, parameter       :: missing = -1e34
  integer, dimension(2) :: axt, axv, axtv, axvt, axto
  integer, dimension(3) :: axt2
  integer               :: id_geo_lon, id_geo_lat, id_sin_rot, id_cos_rot, id_cell_area
  logical               :: sent
  integer               :: id_xto,id_yto

  !
  ! diagnostics MUST use a domain without halos otherwise same as the
  ! regular domain:  Domain (see ice_grid.f90)
  !
  id_xv = diag_axis_init('xv', xb1d(2:im+1), 'degrees_E', 'X','longitude', set_name='ice', Domain2=Domain )
  id_yv = diag_axis_init('yv', yb1d(2:jm+1), 'degrees_N', 'Y','latitude',  set_name='ice', Domain2=Domain )
  id_xb = diag_axis_init('xb', xb1d, 'degrees_E', 'X', 'longitude', set_name='ice', Domain2=Domain )
  id_yb = diag_axis_init('yb', yb1d, 'degrees_N', 'Y', 'latitude', set_name='ice', Domain2=Domain )
  id_xt = diag_axis_init('xt', (xb1d(1:im)+xb1d(2:im+1))/2, 'degrees_E', 'X', &
          'longitude',set_name='ice',edges=id_xb,Domain2=Domain)
  id_yt = diag_axis_init('yt', (yb1d(1:jm)+yb1d(2:jm+1))/2, 'degrees_N', 'Y', &
          'latitude',set_name='ice', edges=id_yb,Domain2=Domain)
  id_ct = diag_axis_init('ct', hlim(1:num_part-1), 'meters','Z', 'thickness')

  id_xto = diag_axis_init ('xt_ocean',grid_x_t,'degrees_E','x','tcell longitude',&
           set_name='ice', Domain2=Domain, aux='geolon_t')
  id_yto = diag_axis_init ('yt_ocean',grid_y_t,'degrees_N','y','tcell latitude',&
           set_name='ice', Domain2=Domain, aux='geolat_t')
  axto = (/ id_xto, id_yto /)
  axv  = (/ id_xv, id_yv  /)
  axt  = (/ id_xt, id_yt  /)
  axt2 = (/ id_xt, id_yt, id_ct/)
  Ice%axes(:) = axt2(:)
  axtv = (/ id_xt, id_yv /); ! for north faces of t-cells
  axvt = (/ id_xv, id_yt /); ! for east  faces of t-cells

  G%axesT1(:) = (/ id_xt, id_yt  /)
  G%axesB1(:) = (/ id_xv, id_yv  /)
  G%axesCv1(:) = (/ id_xt, id_yv /)
  G%axesCu1(:) = (/ id_xv, id_yt /)

  id_sin_rot   = register_static_field('ice_model', 'SINROT', axt,              &
                 '-SINROT,COSROT points north', 'none')
  id_cos_rot   = register_static_field('ice_model', 'COSROT', axt,              &
                 'COSROT,SINROT points east','none')
  id_geo_lon   = register_static_field('ice_model', 'GEOLON', axt, 'longitude', &
                 'degrees')
  id_geo_lat   = register_static_field('ice_model', 'GEOLAT', axt, 'latitude',  &
                 'degrees')
  id_cell_area = register_static_field('ice_model', 'CELL_AREA', axt,           &
                 'cell area', 'sphere')
  id_ext       = register_diag_field('ice_model', 'MOI', axt, Time,   &
                 'ice modeled', '0 or 1', missing_value=missing)
  if (id_ext > 0 ) then
     call error_mesg ('ice_model_init', &
          'Diagnostic MOI has been renamed EXT.  Change your diag_table.', WARNING)
  else
     id_ext = register_diag_field('ice_model', 'EXT', axt, Time, &
              'ice modeled', '0 or 1', missing_value=missing)
  end if
  id_mi       = register_diag_field('ice_model', 'MI', axt, Time,                  &
               'ice mass', 'kg/m^2', missing_value=missing)
  id_mib      = register_diag_field('ice_model', 'MIB', axt, Time,                 &
               'ice + bergs mass', 'kg/m^2', missing_value=missing)
  id_cn       = register_diag_field('ice_model', 'CN', axt2, Time,                 &
               'ice concentration', '0-1', missing_value=missing)
  id_hs       = register_diag_field('ice_model', 'HS', axt, Time,                  &
               'snow thickness', 'm-snow', missing_value=missing)
  id_tsn      = register_diag_field('ice_model', 'TSN', axt, Time,                 &
               'snow layer temperature', 'C',  missing_value=missing)
  id_hi       = register_diag_field('ice_model', 'HI', axt, Time,                  &
               'ice thickness', 'm-ice', missing_value=missing)
  id_hio      = register_diag_field('ice_model', 'HIO', axto, Time,                &
               'ice thickness', 'm-ice', missing_value=missing)
  id_t1       = register_diag_field('ice_model', 'T1', axt, Time,                  &
               'top ice layer temperature', 'C',  missing_value=missing)
  id_t2       = register_diag_field('ice_model', 'T2', axt, Time,                  &
               'second ice layer temperature', 'C',  missing_value=missing)
  id_t3       = register_diag_field('ice_model', 'T3', axt, Time,                  &
               'third ice layer temperature', 'C',  missing_value=missing)
  id_t4       = register_diag_field('ice_model', 'T4', axt, Time,                  &
               'bottom ice layer temperature', 'C',  missing_value=missing)
  id_ts       = register_diag_field('ice_model', 'TS', axt, Time,                  &
               'surface temperature', 'C', missing_value=missing)
  id_sh       = register_diag_field('ice_model','SH' ,axt, Time,                   &
               'sensible heat flux', 'W/m^2',  missing_value=missing)
  id_lh       = register_diag_field('ice_model','LH' ,axt, Time,                   &
               'latent heat flux', 'W/m^2', missing_value=missing)
  id_sw       = register_diag_field('ice_model','SW' ,axt, Time,                   &
               'short wave heat flux', 'W/m^2', missing_value=missing)
  id_lw       = register_diag_field('ice_model','LW' ,axt, Time,                   &
               'long wave heat flux over ice', 'W/m^2', missing_value=missing)
  id_snofl    = register_diag_field('ice_model','SNOWFL' ,axt, Time,               &
               'rate of snow fall', 'kg/(m^2*s)', missing_value=missing)
  id_rain     = register_diag_field('ice_model','RAIN' ,axt, Time,                 &
               'rate of rain fall', 'kg/(m^2*s)', missing_value=missing)
  id_runoff   = register_diag_field('ice_model','RUNOFF' ,axt, Time,               &
               'liquid runoff', 'kg/(m^2*s)', missing_value=missing)
  id_calving  = register_diag_field('ice_model','CALVING',axt, Time,               &
               'frozen runoff', 'kg/(m^2*s)', missing_value=missing)
  id_runoff_hflx   = register_diag_field('ice_model','RUNOFF_HFLX' ,axt, Time,               &
               'liquid runoff sensible heat flux', 'W/m^2', missing_value=missing)
  id_calving_hflx  = register_diag_field('ice_model','CALVING_HFLX',axt, Time,               &
               'frozen runoff sensible heat flux', 'W/m^2', missing_value=missing)
  id_evap     = register_diag_field('ice_model','EVAP',axt, Time,                  &
               'evaporation', 'kg/(m^2*s)', missing_value=missing)
  id_saltf    = register_diag_field('ice_model','SALTF' ,axt, Time,                &
               'ice to ocean salt flux', 'kg/(m^2*s)', missing_value=missing)
  id_sn2ic    = register_diag_field('ice_model','SN2IC'  ,axt,Time,                &
               'rate of snow to ice conversion', 'kg/(m^2*s)', missing_value=missing)
  id_tmelt    = register_diag_field('ice_model','TMELT'  ,axt, Time,               &
               'upper surface melting energy flux', 'W/m^2', missing_value=missing)
  id_bmelt    = register_diag_field('ice_model','BMELT'  ,axt, Time,               &
               'bottom surface melting energy flux', 'W/m^2', missing_value=missing)
  id_bheat    = register_diag_field('ice_model','BHEAT'  ,axt, Time,               &
               'ocean to ice heat flux', 'W/m^2', missing_value=missing)
  id_e2m      = register_diag_field('ice_model','E2MELT' ,axt, Time,               &
               'heat needed to melt ice', 'J/m^2', missing_value=missing)
  id_frazil   = register_diag_field('ice_model','FRAZIL' ,axt, Time,               &
               'energy flux of frazil formation', 'W/m^2', missing_value=missing)
  id_alb      = register_diag_field('ice_model','ALB',axt, Time,                   &
               'surface albedo','0-1', missing_value=missing )
  id_coszen   = register_diag_field('ice_model','coszen',axt, Time,                   &
               'cosine of zenith','-1:1', missing_value=missing )
  id_sw_abs_snow= register_diag_field('ice_model','sw_abs_snow',axt, Time,&
               'SW frac. abs. in snow','0:1', missing_value=missing )
  id_sw_abs_ice1= register_diag_field('ice_model','sw_abs_ice1',axt, Time,&
               'SW frac. abs. in ice1','0:1', missing_value=missing )
  id_sw_abs_ice2= register_diag_field('ice_model','sw_abs_ice2',axt, Time,&
               'SW frac. abs. in ice2','0:1', missing_value=missing )
  id_sw_abs_ice3= register_diag_field('ice_model','sw_abs_ice3',axt, Time,&
               'SW frac. abs. in ice3','0:1', missing_value=missing )
  id_sw_abs_ice4= register_diag_field('ice_model','sw_abs_ice4',axt, Time,&
               'SW frac. abs. in ice4','0:1', missing_value=missing )
  id_sw_pen= register_diag_field('ice_model','sw_pen',axt, Time,&
               'SW frac. pen. surf.','0:1', missing_value=missing )
  id_sw_trn= register_diag_field('ice_model','sw_trn',axt, Time,&
               'SW frac. trans. to ice bot.','0:1', missing_value=missing )



  id_alb_vis_dir = register_diag_field('ice_model','alb_vis_dir',axt, Time,                &
               'ice surface albedo vis_dir','0-1', missing_value=missing )
  id_alb_vis_dif = register_diag_field('ice_model','alb_vis_dif',axt, Time,                &
               'ice surface albedo vis_dif','0-1', missing_value=missing )
  id_alb_nir_dir = register_diag_field('ice_model','alb_nir_dir',axt, Time,                &
               'ice surface albedo nir_dir','0-1', missing_value=missing )
  id_alb_nir_dif = register_diag_field('ice_model','alb_nir_dif',axt, Time,                &
               'ice surface albedo nir_dif','0-1', missing_value=missing )
  id_xprt     = register_diag_field('ice_model','XPRT',axt, Time,                  &
               'frozen water transport convergence', 'kg/(m^2*yr)', missing_value=missing)
  id_lsrc     = register_diag_field('ice_model','LSRC', axt, Time,                 &
               'frozen water local source', 'kg/(m^2*yr)', missing_value=missing)
  id_lsnk     = register_diag_field('ice_model','LSNK',axt, Time,                  &
               'frozen water local sink', 'kg/(m^2*yr)', missing_value=missing)
  id_bsnk     = register_diag_field('ice_model','BSNK',axt, Time,                  &
               'frozen water local bottom sink', 'kg/(m^2*yr)', missing_value=missing)
  id_qfres    = register_diag_field('ice_model', 'QFLX_RESTORE_ICE', axt, Time,    &
               'Ice Restoring heat flux', 'W/m^2', missing_value=missing)
  id_qflim    = register_diag_field('ice_model', 'QFLX_LIMIT_ICE', axt, Time,      &
               'Ice Limit heat flux', 'W/m^2', missing_value=missing)
  id_strna    = register_diag_field('ice_model','STRAIN_ANGLE', axt,Time,          &
               'strain angle', 'none', missing_value=missing)
  id_fax      = register_diag_field('ice_model', 'FA_X', axv, Time,                &
               'air stress on ice - x component', 'Pa', missing_value=missing)
  id_fay      = register_diag_field('ice_model', 'FA_Y', axv, Time,                &
               'air stress on ice - y component', 'Pa', missing_value=missing)
  id_uo       = register_diag_field('ice_model', 'UO', axv, Time,                  &
               'surface current - x component', 'm/s', missing_value=missing)
  id_vo       = register_diag_field('ice_model', 'VO', axv, Time,                  &
               'surface current - y component', 'm/s', missing_value=missing)
  id_sw_vis   = register_diag_field('ice_model','SW_VIS' ,axt, Time,               &
               'visible short wave heat flux', 'W/m^2', missing_value=missing)
  id_sw_dir   = register_diag_field('ice_model','SW_DIR' ,axt, Time,               &
               'direct short wave heat flux', 'W/m^2', missing_value=missing)
  id_sw_dif   = register_diag_field('ice_model','SW_DIF' ,axt, Time,               &
               'diffuse short wave heat flux', 'W/m^2', missing_value=missing)
  id_sw_vis_dir = register_diag_field('ice_model','SW_VIS_DIR' ,axt, Time,         &
               'visible direct short wave heat flux', 'W/m^2', missing_value=missing)
  id_sw_vis_dif = register_diag_field('ice_model','SW_VIS_DIF' ,axt, Time,         &
               'visible diffuse short wave heat flux', 'W/m^2', missing_value=missing)
  id_sw_nir_dir = register_diag_field('ice_model','SW_NIR_DIR' ,axt, Time,         &
               'near IR direct short wave heat flux', 'W/m^2', missing_value=missing)
  id_sw_nir_dif = register_diag_field('ice_model','SW_NIR_DIF' ,axt, Time,         &
               'near IR diffuse short wave heat flux', 'W/m^2', missing_value=missing)

  !
  ! diagnostics for quantities produced outside the ice model
  !
  id_swdn  = register_diag_field('ice_model','SWDN' ,axt, Time,       &
             'downward shortwave flux', 'W/m^2', missing_value=missing)
  id_lwdn  = register_diag_field('ice_model','LWDN' ,axt, Time,       &
             'downward longwave flux', 'W/m^2', missing_value=missing)
  id_ta    = register_diag_field('ice_model', 'TA', axt, Time,        &
             'surface air temperature', 'C', missing_value=missing)
  id_slp   = register_diag_field('ice_model', 'SLP', axt, Time,       &
             'sea level pressure', 'Pa', missing_value=missing)
  id_sst   = register_diag_field('ice_model', 'SST', axt, Time,       &
             'sea surface temperature', 'deg-C', missing_value=missing)
  id_sss   = register_diag_field('ice_model', 'SSS', axt, Time,       &
             'sea surface salinity', 'psu', missing_value=missing)
  id_ssh   = register_diag_field('ice_model', 'SSH', axt, Time,       &
             'sea surface height', 'm', missing_value=missing)
  id_obi   = register_diag_field('ice_model', 'OBI', axt, Time,       &
       'ice observed', '0 or 1', missing_value=missing)

  if (id_sin_rot>0)   sent=send_data(id_sin_rot, sin_rot(isc:iec,jsc:jec), Time);
  if (id_cos_rot>0)   sent=send_data(id_cos_rot, cos_rot(isc:iec,jsc:jec), Time);
  if (id_geo_lon>0)   sent=send_data(id_geo_lon, Ice%G%geoLonT(isc:iec,jsc:jec), Time);
  if (id_geo_lat>0)   sent=send_data(id_geo_lat, Ice%G%geoLatT(isc:iec,jsc:jec), Time);
  if (id_cell_area>0) sent=send_data(id_cell_area, cell_area, Time);

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
!  write(outunit,100) 'ice_data_type%qflx_lim_ice       ',mpp_chksum(IST%qflx_lim_ice)
!  write(outunit,100) 'ice_data_type%qflx_res_ice       ',mpp_chksum(IST%qflx_res_ice)

  do n=1,Ice%ocean_fields%num_bcs ; do m=1,Ice%ocean_fields%bc(n)%num_fields
    write(outunit,101) 'ice%', trim(Ice%ocean_fields%bc(n)%name), &
                       trim(Ice%ocean_fields%bc(n)%field(m)%name), &
                       mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
  enddo ; enddo

100 FORMAT("   CHECKSUM::",A32," = ",Z20)
101 FORMAT("   CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ice_data_type_chksum

end module ice_type_mod
