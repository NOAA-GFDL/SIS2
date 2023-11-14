!> maintains the sea ice data, reads/writes restarts, reads the namelist and initializes diagnostics.
module ice_type_mod

use ice_bergs,         only : icebergs, icebergs_stock_pe, icebergs_save_restart
use ice_grid,          only : ice_grid_type
use MOM_coms,          only : PE_here
use MOM_domains,       only : CGRID_NE, BGRID_NE, AGRID
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : is_root_pe, stdout
use MOM_file_parser,   only : param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_time_manager,  only : time_type, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type
use SIS_ctrl_types,    only : SIS_fast_CS, SIS_slow_CS
use SIS_debugging,     only : chksum
use SIS_diag_mediator, only : SIS_diag_ctrl, post_data=>post_SIS_data, register_SIS_diag_field
use SIS_framework,     only : domain2D, SIS_chksum, get_domain_extent, safe_alloc, safe_alloc_ptr
use SIS_restart,       only : register_restart_field, save_restart, SIS_restart_CS, query_initialized
use SIS_framework,     only : coupler_1d_bc_type, coupler_2d_bc_type, coupler_3d_bc_type
use SIS_framework,     only : coupler_type_spawn, coupler_type_write_chksums
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_open_boundary, only : ice_OBC_type
use SIS_types,         only : ice_state_type, fast_ice_avg_type
use SIS2_ice_thm,      only : ice_thermo_type, energy_0degC, get_SIS2_thermo_coefs
use iso_fortran_env,   only : int64

implicit none ; private

public :: ice_data_type, dealloc_ice_arrays
public :: ice_type_slow_reg_restarts, ice_type_fast_reg_restarts
public :: ice_model_restart, ice_stock_pe, ice_data_type_chksum
public :: Ice_public_type_chksum, Ice_public_type_bounds_check

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> This structure contains the ice model data (some used by calling routines);   !
!! the third index is partition (1 is open water; 2... are ice cover by category)!
type ice_data_type !  ice_public_type
  type(domain2D)          :: Domain         !< A copy of the fast ice domain without halos.
  type(domain2D)          :: slow_Domain_NH !< A copy of the slow ice domain without halos.
  type(domain2D), pointer :: &
    fast_domain => NULL(), & !< A pointer to the fast ice mpp domain or a copy
                             !! on slow ice PEs.
    slow_domain => NULL()    !< A pointer to the fast ice mpp domain or a copy
                             !! on slow ice PEs.
  type(time_type) :: Time    !< The sea-ice model's clock, that
                             !! set with the current model time.
  logical  :: pe             !< If true, there is ice on this PE.
  logical  :: slow_ice_pe = .false. !< If true, this is a slow ice PE
  logical  :: fast_ice_pe = .false. !< If true, this is a fast ice PE
  logical  :: shared_slow_fast_PEs = .true. !< If true, the fast and slow ice use the same processors
                                    !! and domain decomposition
  integer  :: xtype          !< An integer specifying the type for the exchange
  integer, pointer, dimension(:)   :: slow_pelist =>NULL() !< Used for flux-exchange with slow processes.
  integer, pointer, dimension(:)   :: fast_pelist =>NULL() !< Used for flux-exchange with fast processes.
  integer, pointer, dimension(:)   :: pelist   =>NULL() !< Used for flux-exchange.
  logical, pointer, dimension(:,:) :: ocean_pt =>NULL() !< An array that indicates ocean points as true.

  ! These fields are used to provide information about the ice surface to the
  ! atmosphere, and contain separate values for each ice thickness category.
  real, pointer, dimension(:,:,:) :: &
    part_size => NULL(), &    !< The fractional coverage of a grid cell by each ice
                              !! thickness category [nondim], 0 to 1.  Category 1 is
                              !! open ocean.  The sum of part_size is 1.
    albedo    => NULL(), &    !< The surface albedo averaged across all wavelength
                              !! and orientation bands within each ice-thickness
                              !! category [nondim], between 0 and 1.
    albedo_vis_dir => NULL(), & !< The surface albedo for direct visible shortwave radiation
                                !! in each ice-thickness category [nondim], between 0 and 1.
    albedo_nir_dir => NULL(), & !< The surface albedo for diffuse visible shortwave radiation
                                !! in each ice-thickness category [nondim], between 0 and 1.
    albedo_vis_dif => NULL(), & !< The surface albedo for direct near-infrared shortwave radiation
                                !! in each ice-thickness category [nondim], between 0 and 1.
    albedo_nir_dif => NULL(), & !< The surface albedo for diffuse near-infrared shortwave radiation
                                !! in each ice-thickness category [nondim], between 0 and 1.
    rough_mom   => NULL(), &  !< The roughness for momentum at the ocean surface, as provided by
                              !! ocean_rough_mod, apparently [m].
    rough_heat  => NULL(), &  !< The roughness for heat at the ocean surface, as provided by
                              !! ocean_rough_mod, apparently [m].
    rough_moist => NULL(), &  !< The roughness for moisture at the ocean surface, as provided by
                              !! ocean_rough_mod, apparently [m].
    t_surf      => NULL(), &  !< The surface temperature for the ocean or for
                              !! each ice-thickness category [Kelvin].
    u_surf      => NULL(), &  !< The eastward surface velocities of the ocean (:,:,1) or sea-ice [m s-1].
    v_surf      => NULL()     !< The northward surface velocities of the ocean (:,:,1) or sea-ice [m s-1].
  real, pointer, dimension(:,:)   :: &
    s_surf         =>NULL()   !< The ocean's surface salinity [gSalt kg-1].

  ! These arrays will be used to set the forcing for the ocean.
  real, pointer, dimension(:,:) :: &
    SST_C => NULL(), &    !< The ocean surface temperature [degC].
    flux_u => NULL(), &   !< The flux of x-momentum into the ocean [Pa].
    flux_v => NULL(), &   !< The flux of y-momentum into the ocean [Pa].
    flux_t => NULL(), &   !< The flux of sensible heat out of the ocean [W m-2].
    flux_q => NULL(), &   !< The evaporative moisture flux out of the ocean [kg m-2 s-1].
    flux_lw => NULL(), &  !< The longwave flux out of the ocean [W m-2].
    flux_sw_vis_dir => NULL(), & !< The direct visible shortwave heat flux into the ocean [W m-2].
    flux_sw_vis_dif => NULL(), & !< The diffuse visible shortwave heat flux into the ocean [W m-2].
    flux_sw_nir_dir => NULL(), & !< The direct near-infrared heat flux into the ocean [W m-2].
    flux_sw_nir_dif => NULL(), & !< The diffuse near-infrared heat flux into the ocean [W m-2].
    flux_lh => NULL(), &  !< The latent heat flux out of the ocean [W m-2].
    lprec => NULL(), &    !< The liquid precipitation flux into the ocean [kg m-2].
    fprec => NULL(), &    !< The frozen precipitation flux into the ocean [kg m-2].
    p_surf => NULL(), &   !< The pressure at the ocean surface [Pa].  This may
                          !! or may not include atmospheric pressure.
    runoff => NULL(), &   !< Liquid runoff into the ocean [kg m-2].
    calving => NULL(), &  !< Calving of ice or runoff of frozen fresh water into
                          !! the ocean [kg m-2].
    stress_mag => NULL(), & !< The time-mean magnitude of the stress on the ocean [Pa].
    ustar_berg => NULL(), &  !< ustar contribution below icebergs [m s-1]
    area_berg => NULL(),  &  !< fraction of grid cell covered by icebergs in [m2 m-2]
    mass_berg => NULL(),  &  !< mass of icebergs in [kg m-2]
    runoff_hflx => NULL(), &  !< The heat flux associated with runoff, based on
                              !! the temperature difference relative to a
                              !! reference temperature, in ???.
    calving_hflx => NULL(), & !< The heat flux associated with calving, based on
                              !! the temperature difference relative to a
                              !! reference temperature, in ???.
    flux_salt  => NULL()  !< The flux of salt out of the ocean [kg m-2 s-1].

  real, pointer, dimension(:,:) :: &
    area => NULL() , &    !< The area of ocean cells [m2].  Land cells have
                          !! a value of 0, so this could also be used as a mask.
    mi   => NULL()        !< The total ice+snow mass [kg m-2].
             ! mi is needed for the wave model. It is introduced here,
             ! because flux_ice_to_ocean cannot handle 3D fields. This may be
             ! removed, if the information on ice thickness can be derived from
             ! h_ice outside the ice module.
  integer, dimension(3)    :: axes  !< The sea ice surface field axes.
  type(coupler_3d_bc_type) :: ocean_fields !< array of fields used for additional tracers
                                           !! whose surface state is shared with the atmosphere.
  type(coupler_2d_bc_type) :: ocean_fluxes !< array of fluxes from the ice to the ocean used
                                           !! for additional tracers
  type(coupler_3d_bc_type) :: ocean_fluxes_top  !< An ARCHAIC element that should eventually be deleted.
                               ! ###THIS IS ARCHAIC AND COULD BE DELETED!
  integer :: flux_uv_stagger = -999 !< The staggering relative to the tracer points of the two
                          !! wind stress components. Valid entries include AGRID, BGRID_NE,
                          !! CGRID_NE, BGRID_SW, and CGRID_SW, corresponding to the community-
                          !! standard Arakawa notation. (These are named integers taken from
                          !! mpp_parameter_mod.) Following SIS, this is BGRID_NE by default when the
                          !! sea ice is initialized, but here it is set to -999 so that a global max
                          !! across ice and non-ice processors can be used to determine its value.

  ! The following are actually private to SIS2, and are not used elsewhere by other FMS modules.
  type(icebergs),    pointer :: icebergs => NULL()
          !< A pointer to the icebergs control structure
  type(SIS_fast_CS), pointer :: fCS => NULL()
          !< A pointer to the SIS fast ice update control structure
  type(SIS_slow_CS), pointer :: sCS => NULL()
          !< A pointer to the SIS slow ice update control structure
  type(unit_scale_type), pointer :: US => NULL()
          !< structure containing various unit conversion factors
  type(SIS_restart_CS), pointer :: Ice_restart => NULL()
          !< A pointer to the slow ice restart control structure
  type(SIS_restart_CS), pointer :: Ice_fast_restart => NULL()
          !< A pointer to the fast ice restart control structure
  type(ice_OBC_type), pointer :: OBC => NULL()
          !< A pointer to the ice OBC control structure
  character(len=240) :: restart_output_dir = './RESTART/'
          !< The directory into which to write restart files.
end type ice_data_type !  ice_public_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_type_slow_reg_restarts allocates the arrays in the ice_data_type that are
!! predominantly associated with the slow processors, and register any variables
!! in the ice data type that need to be included in the slow ice restart files.
subroutine ice_type_slow_reg_restarts(domain, CatIce, param_file, Ice, &
                                      Ice_restart, gas_fluxes)
  type(domain2d),          intent(in)    :: domain   !< The ice models' FMS domain type
  integer,                 intent(in)    :: CatIce   !< The number of ice thickness categories
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ice_data_type),     intent(inout) :: Ice      !< The publicly visible ice data type.
  type(SIS_restart_CS),    pointer       :: Ice_restart !< The control structure for the ice restarts
  type(coupler_1d_bc_type), &
                 optional, intent(in)    :: gas_fluxes !< If present, this type describes the
                                              !! additional gas or other tracer fluxes between the
                                              !! ocean, ice, and atmosphere.

  ! This subroutine allocates the externally visible ice_data_type's arrays and
  ! registers the appropriate ones for inclusion in the restart file.
  integer :: isc, iec, jsc, jec, km, idr

  call get_domain_extent(domain, isc, iec, jsc, jec )
  km = CatIce + 1

  ! The fields t_surf, s_surf, and part_size are only available on fast PEs.

  call safe_alloc_ptr(Ice%flux_u, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_v, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_t, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_q, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_sw_vis_dir, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_sw_vis_dif, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_sw_nir_dir, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_sw_nir_dif, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_lw, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_lh, isc, iec, jsc, jec)  !NI
  call safe_alloc_ptr(Ice%lprec, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%fprec, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%p_surf, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%runoff, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%calving, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%runoff_hflx, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%calving_hflx, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%flux_salt, isc, iec, jsc, jec)

  call safe_alloc_ptr(Ice%SST_C, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%area, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%mi, isc, iec, jsc, jec)  !NR

  if (Ice%sCS%pass_stress_mag) then
    call safe_alloc_ptr(Ice%stress_mag, isc, iec, jsc, jec)
  endif

  if (Ice%sCS%pass_iceberg_area_to_ocean) then
    call safe_alloc_ptr(Ice%ustar_berg, isc, iec, jsc, jec)
    call safe_alloc_ptr(Ice%area_berg, isc, iec, jsc, jec)
    call safe_alloc_ptr(Ice%mass_berg, isc, iec, jsc, jec)
  endif

  if (present(gas_fluxes)) &
    call coupler_type_spawn(gas_fluxes, Ice%ocean_fluxes, (/isc,isc,iec,iec/), &
                            (/jsc,jsc,jec,jec/),  suffix = '_ice')

  ! These are used by the ocean model, and need to be in the slow PE restarts.
  if (associated(Ice_restart)) then
    call register_restart_field(Ice_restart, 'flux_u',      Ice%flux_u)
    call register_restart_field(Ice_restart, 'flux_v',      Ice%flux_v)
    call register_restart_field(Ice_restart, 'flux_t',      Ice%flux_t)
    call register_restart_field(Ice_restart, 'flux_q',      Ice%flux_q)
    call register_restart_field(Ice_restart, 'flux_salt',   Ice%flux_salt)
    call register_restart_field(Ice_restart, 'flux_lw',     Ice%flux_lw)
    call register_restart_field(Ice_restart, 'lprec',       Ice%lprec)
    call register_restart_field(Ice_restart, 'fprec',       Ice%fprec)
    call register_restart_field(Ice_restart, 'runoff',      Ice%runoff)
    call register_restart_field(Ice_restart, 'calving',     Ice%calving)
    call register_restart_field(Ice_restart, 'runoff_hflx', Ice%runoff_hflx, mandatory=.false.)
    call register_restart_field(Ice_restart, 'calving_hflx',Ice%calving_hflx, mandatory=.false.)
    call register_restart_field(Ice_restart, 'p_surf',      Ice%p_surf)
    call register_restart_field(Ice_restart, 'flux_sw_vis_dir', Ice%flux_sw_vis_dir)
    call register_restart_field(Ice_restart, 'flux_sw_vis_dif', Ice%flux_sw_vis_dif)
    call register_restart_field(Ice_restart, 'flux_sw_nir_dir', Ice%flux_sw_nir_dir)
    call register_restart_field(Ice_restart, 'flux_sw_nir_dif', Ice%flux_sw_nir_dif)
    if (Ice%sCS%pass_stress_mag) &
      call register_restart_field(Ice_restart, 'stress_mag', Ice%stress_mag, mandatory=.false.)
    if (Ice%sCS%pass_iceberg_area_to_ocean) then
      call register_restart_field(Ice_restart, 'ustar_berg', Ice%ustar_berg, mandatory=.false.)
      call register_restart_field(Ice_restart, 'area_berg', Ice%area_berg, mandatory=.false.)
      call register_restart_field(Ice_restart, 'mass_berg', Ice%mass_berg, mandatory=.false.)
    endif
  endif
end subroutine ice_type_slow_reg_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_type_fast_reg_restarts allocates the arrays in the ice_data_type that are
!! predominantly associated with the fast processors, and registers any variables
!! in the ice data type that need to be included in the fast ice restart files.
subroutine ice_type_fast_reg_restarts(domain, CatIce, param_file, Ice, &
                                      Ice_restart, gas_fields_ocn)
  type(domain2d),          intent(in)    :: domain   !< The ice models' FMS domain type
  integer,                 intent(in)    :: CatIce   !< The number of ice thickness categories
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ice_data_type),     intent(inout) :: Ice !< The publicly visible ice data type.
  type(SIS_restart_CS),    pointer       :: Ice_restart !< The control structure for the ice restarts
  type(coupler_1d_bc_type), &
                 optional, intent(in)    :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes.

  ! This subroutine allocates the externally visible ice_data_type's arrays and
  ! registers the appropriate ones for inclusion in the restart file.
  integer :: isc, iec, jsc, jec, km, idr

  call get_domain_extent(domain, isc, iec, jsc, jec )
  km = CatIce + 1

  call safe_alloc_ptr(Ice%t_surf, isc, iec, jsc, jec, km)
  call safe_alloc_ptr(Ice%s_surf, isc, iec, jsc, jec)
  call safe_alloc_ptr(Ice%part_size, isc, iec, jsc, jec, km)

  call safe_alloc_ptr(Ice%u_surf, isc, iec, jsc, jec, km)
  call safe_alloc_ptr(Ice%v_surf, isc, iec, jsc, jec, km)
  if (.not.associated(Ice%ocean_pt)) then
    allocate(Ice%ocean_pt(isc:iec, jsc:jec), source=.false.) !derived
  endif

  call safe_alloc_ptr(Ice%rough_mom, isc, iec, jsc, jec, km)
  call safe_alloc_ptr(Ice%rough_heat, isc, iec, jsc, jec, km)
  call safe_alloc_ptr(Ice%rough_moist, isc, iec, jsc, jec, km)

  call safe_alloc_ptr(Ice%albedo, isc, iec, jsc, jec, km)  ! Derived?
  call safe_alloc_ptr(Ice%albedo_vis_dir, isc, iec, jsc, jec, km)
  call safe_alloc_ptr(Ice%albedo_nir_dir, isc, iec, jsc, jec, km)
  call safe_alloc_ptr(Ice%albedo_vis_dif, isc, iec, jsc, jec, km)
  call safe_alloc_ptr(Ice%albedo_nir_dif, isc, iec, jsc, jec, km)

  if (present(gas_fields_ocn)) &
    call coupler_type_spawn(gas_fields_ocn, Ice%ocean_fields, (/isc,isc,iec,iec/), &
                            (/jsc,jsc,jec,jec/), (/1, km/), suffix = '_ice')

  ! Now register some of these arrays to be read from the restart files.
  ! These are used by the atmospheric model, and need to be in the fast PE restarts.
  if (associated(Ice_restart)) then
    call register_restart_field(Ice_restart, 'rough_mom',   Ice%rough_mom, dim_3="cat0")
    call register_restart_field(Ice_restart, 'rough_heat',  Ice%rough_heat, dim_3="cat0")
    call register_restart_field(Ice_restart, 'rough_moist', Ice%rough_moist, dim_3="cat0")
  endif

end subroutine ice_type_fast_reg_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dealloc_Ice_arrays deallocates the memory in the publicly visible ice data type.
subroutine dealloc_Ice_arrays(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type.

  if (associated(Ice%ocean_pt)) deallocate(Ice%ocean_pt)
  if (associated(Ice%t_surf)) deallocate(Ice%t_surf)
  if (associated(Ice%s_surf)) deallocate(Ice%s_surf)
  if (associated(Ice%u_surf)) deallocate(Ice%u_surf)
  if (associated(Ice%v_surf)) deallocate(Ice%v_surf)
  if (associated(Ice%part_size)) deallocate(Ice%part_size)
  if (associated(Ice%rough_mom)) deallocate(Ice%rough_mom)
  if (associated(Ice%rough_heat)) deallocate(Ice%rough_heat)
  if (associated(Ice%rough_moist)) deallocate(Ice%rough_moist)
  if (associated(Ice%albedo)) deallocate(Ice%albedo)
  if (associated(Ice%albedo_vis_dir)) deallocate(Ice%albedo_vis_dir)
  if (associated(Ice%albedo_nir_dir)) deallocate(Ice%albedo_nir_dir)
  if (associated(Ice%albedo_vis_dif)) deallocate(Ice%albedo_vis_dif)
  if (associated(Ice%albedo_nir_dif)) deallocate(Ice%albedo_nir_dif)

  if (associated(Ice%flux_u)) deallocate(Ice%flux_u)
  if (associated(Ice%flux_v)) deallocate(Ice%flux_v)
  if (associated(Ice%flux_t)) deallocate(Ice%flux_t)
  if (associated(Ice%flux_q)) deallocate(Ice%flux_q)
  if (associated(Ice%flux_lw)) deallocate(Ice%flux_lw)
  if (associated(Ice%flux_lh)) deallocate(Ice%flux_lh)
  if (associated(Ice%lprec)) deallocate(Ice%lprec)
  if (associated(Ice%fprec)) deallocate(Ice%fprec)
  if (associated(Ice%p_surf)) deallocate(Ice%p_surf)
  if (associated(Ice%runoff)) deallocate(Ice%runoff)
  if (associated(Ice%calving)) deallocate(Ice%calving)
  if (associated(Ice%runoff_hflx)) deallocate(Ice%runoff_hflx)
  if (associated(Ice%calving_hflx)) deallocate(Ice%calving_hflx)
  if (associated(Ice%stress_mag)) deallocate(Ice%stress_mag)

  if (associated(Ice%flux_salt)) deallocate(Ice%flux_salt)
  if (associated(Ice%flux_sw_vis_dir)) deallocate(Ice%flux_sw_vis_dir)
  if (associated(Ice%flux_sw_vis_dif)) deallocate(Ice%flux_sw_vis_dif)
  if (associated(Ice%flux_sw_nir_dir)) deallocate(Ice%flux_sw_nir_dir)
  if (associated(Ice%flux_sw_nir_dif)) deallocate(Ice%flux_sw_nir_dif)
  if (associated(Ice%area)) deallocate(Ice%area)
  if (associated(Ice%mi)) deallocate(Ice%mi)

  if (associated(Ice%ustar_berg)) deallocate(Ice%ustar_berg)
  if (associated(Ice%area_berg)) deallocate(Ice%area_berg)
  if (associated(Ice%mass_berg)) deallocate(Ice%mass_berg)

end subroutine dealloc_Ice_arrays

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Ice_public_type_chksum writes out checksums of the variables in a publicly
!! visible ice data type
subroutine Ice_public_type_chksum(mesg, Ice, check_fast, check_slow, check_rough)
  character(len=*),    intent(in) :: mesg !< An identifying message
  type(ice_data_type), intent(in) :: Ice !< The publicly visible ice data type.
  logical, optional,   intent(in) :: check_fast !< If true, check the fast ice fields
  logical, optional,   intent(in) :: check_slow !< If true, check the slow ice fields
  logical, optional,   intent(in) :: check_rough !< If true, check the roughness fields
!   This subroutine writes out chksums for the model's basic state variables.

  ! Note that the publicly visible ice_data_type has no halos, so it is not
  ! possible do check their values.

  logical :: fast_fields, slow_fields, roughnesses

  fast_fields = Ice%fast_ice_PE ; slow_fields = Ice%slow_ice_PE

  if (present(check_fast)) then
    fast_fields = Ice%fast_ice_PE .and. check_fast
    if (.not.present(check_slow)) slow_fields = .false.
  elseif (present(check_slow)) then
    slow_fields = Ice%slow_ice_PE .and. check_slow
  endif

  roughnesses = fast_fields ; if (present(check_rough)) roughnesses = Ice%fast_ice_PE .and. check_rough

  if (fast_fields) then ! This is a fast-ice PE.
    call chksum(Ice%part_size, trim(mesg)//" Ice%part_size")
    call chksum(Ice%albedo, trim(mesg)//" Ice%albedo")
    call chksum(Ice%albedo_vis_dir, trim(mesg)//" Ice%albedo_vis_dir")
    call chksum(Ice%albedo_nir_dir, trim(mesg)//" Ice%albedo_nir_dir")
    call chksum(Ice%albedo_vis_dif, trim(mesg)//" Ice%albedo_vis_dif")
    call chksum(Ice%albedo_nir_dif, trim(mesg)//" Ice%albedo_nir_dif")

    call chksum(Ice%t_surf, trim(mesg)//" Ice%t_surf")
    call chksum(Ice%s_surf, trim(mesg)//" Ice%s_surf")
    call chksum(Ice%u_surf, trim(mesg)//" Ice%u_surf")
    call chksum(Ice%v_surf, trim(mesg)//" Ice%v_surf")
  endif
  if (roughnesses) then
    call chksum(Ice%rough_mom, trim(mesg)//" Ice%rough_mom")
    call chksum(Ice%rough_heat, trim(mesg)//" Ice%rough_heat")
    call chksum(Ice%rough_moist, trim(mesg)//" Ice%rough_moist")
  endif

  if (slow_fields) then ! This is a slow-ice PE.
    call chksum(Ice%SST_C, trim(mesg)//" Ice%SST_C")
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
    if (associated(Ice%sCS)) then ; if (Ice%sCS%pass_stress_mag) then
      call chksum(Ice%stress_mag, trim(mesg)//" Ice%stress_mag")
    endif ; endif
  endif

  if (slow_fields .and. associated(Ice%sCS)) then ; if (Ice%sCS%pass_iceberg_area_to_ocean) then
    call chksum(Ice%ustar_berg, trim(mesg)//" Ice%ustar_berg")
    call chksum(Ice%area_berg, trim(mesg)//" Ice%area_berg")
    call chksum(Ice%mass_berg, trim(mesg)//" Ice%mass_berg")
  endif ; endif
end subroutine Ice_public_type_chksum

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Ice_public_type_bounds_check checks for unphysical values in a publicly
!! visible ice data type, and writes out diagnostics for any offending columns
subroutine Ice_public_type_bounds_check(Ice, G, msg)
  type(ice_data_type),     intent(in)    :: Ice !< The publicly visible ice data type.
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  character(len=*),        intent(in)    :: msg !< An identifying message

  character(len=512) :: mesg1, mesg2
  integer :: i, j, k, l, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  integer :: n_bad, i_bad, j_bad, k_bad
  logical :: fluxes_avail
  real    :: t_min, t_max
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc
  ncat = SIZE(Ice%t_surf,3) - 1

  fluxes_avail = .false. ! (associated(Ice%flux_t) .and. associated(Ice%flux_lw))

  n_bad = 0 ; i_bad = 0 ; j_bad = 0 ; k_bad = 0

  t_min = T_0degC-100. ; t_max = T_0degC+60.
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    if ((Ice%t_surf(i2,j2,k2) < t_min) .or. (Ice%t_surf(i2,j2,k2) > t_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; endif
    endif
    endif ; enddo ; enddo ; enddo
  do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then ; i2 = i+i_off ; j2 = j+j_off
    if ((Ice%s_surf(i2,j2) < 0.0) .or. (Ice%s_surf(i2,j2) > 100.0)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; endif
    endif
    endif ; enddo ; enddo
  if (fluxes_avail) then ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then ; i2 = i+i_off ; j2 = j+j_off
    if ((abs(Ice%flux_t(i2,j2)) > 1e4) .or. (abs(Ice%flux_lw(i2,j2)) > 1e4)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; endif
    endif
    endif ; enddo ; enddo ; endif

  if (n_bad > 0) then
    i2 = i_bad+i_off ; j2 = j_bad+j_off ; k2 = k_bad+1
    write(mesg1,'(" at ", 2(F6.1)," or i,j,k = ",3i4,"; nbad = ",i6," on pe ",i4)') &
           G%geolonT(i_bad,j_bad), G%geolatT(i_bad,j_bad), i_bad, j_bad, k_bad, n_bad, pe_here()
    if (fluxes_avail) then
      write(mesg2,'("T_sfc = ",1pe12.4,", ps = ",1pe12.4,", flux_t,lw,q = ",3(1pe12.4))') &
         Ice%t_surf(i2,j2,k2), Ice%part_size(i2,j2,k2), Ice%flux_t(i2,j2), Ice%flux_lw(i2,j2), Ice%flux_q(i2,j2)
    else
      write(mesg2,'("T_sfc = ",1pe12.4,", ps = ",1pe12.4,", S_sfc = ",1pe12.4)') &
       Ice%t_surf(i2,j2,k2), Ice%part_size(i2,j2,k2), Ice%s_surf(i2,j2)
    endif
    call SIS_error(WARNING, "Bad ice data "//trim(msg)//" ; "//trim(mesg1)//" ; "//trim(mesg2), all_print=.true.)
  endif

end subroutine Ice_public_type_bounds_check


!=======================================================================
! <SUBROUTINE NAME="ice_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!>  Write out restart files registered through register_restart_file
subroutine ice_model_restart(Ice, time_stamp)
  type(ice_data_type),        intent(inout) :: Ice !< The publicly visible ice data type.
  character(len=*), optional, intent(in)    :: time_stamp !< A date stamp to include in the restart file name

  if (associated(Ice%Ice_restart) .and. associated(Ice%sCS)) then
    call save_restart(Ice%restart_output_dir, Ice%Time, Ice%sCS%G, Ice%Ice_restart, IG=Ice%sCS%IG, &
                      time_stamp=time_stamp)
    if (associated(Ice%Ice_fast_restart)) then
      if (.not.associated(Ice%Ice_fast_restart, Ice%Ice_restart)) &
        call save_restart(Ice%restart_output_dir, Ice%Time, Ice%fCS%G, Ice%Ice_fast_restart, &
                          IG=Ice%fCS%IG, time_stamp=time_stamp)
    endif
  elseif (associated(Ice%Ice_fast_restart)) then
    call save_restart(Ice%restart_output_dir, Ice%Time, Ice%fCS%G, Ice%Ice_fast_restart, &
                      IG=Ice%fCS%IG, time_stamp=time_stamp)
  endif
  call icebergs_save_restart(Ice%icebergs, time_stamp)

end subroutine ice_model_restart
! </SUBROUTINE>
!=======================================================================

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_stock_pe returns stocks of heat, water, etc. for conservation checks
subroutine ice_stock_pe(Ice, index, value)

  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT

  type(ice_data_type), intent(in) :: Ice !< The publicly visible ice data type.
  integer, intent(in) :: index !< A coded integer indicating which stock to find
  real, intent(out)   :: value !< The integrated stock quantity

  type(ice_state_type), pointer :: IST => NULL()
  real :: icebergs_value
  real :: LI  ! Latent heat of fusion [Q ~> J kg-1]
  real :: part_area ! The area of an ice thickness partition in a cell [m2]
  real :: kg_H    ! A conversion factor from the ice thickness units to kg m-2 [kg m-2 H-1 ~> 1]
  real :: kg_H_Nk ! The ice thickness unit conversion factor divided by the number of ice layers [kg m-2 H-1 ~> 1]
  integer :: i, j, k, m, isc, iec, jsc, jec, ncat, NkIce
  logical :: slab_ice    ! If true, use the very old slab ice thermodynamics,
                         ! with effectively zero heat capacity of ice and snow.
  type(SIS_hor_grid_type), pointer :: G => NULL()

  value = 0.0
  if(.not.Ice%pe) return

  if (associated(Ice%sCS)) then
    IST => Ice%sCS%IST
    G => Ice%sCS%G
    ncat = Ice%sCS%IG%CatIce ; NkIce = Ice%sCS%IG%NkIce
  elseif (associated(Ice%fCS)) then
    IST => Ice%fCS%IST
    G => Ice%fCS%G
    ncat = Ice%fCS%IG%CatIce ; NkIce = Ice%fCS%IG%NkIce
  else
    call SIS_error(WARNING, "ice_stock_pe called with an ice_data_type "//&
                   "without either sCS or fCS associated.")
    return
  endif

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  kg_H = G%US%RZ_to_kg_m2 ; kg_H_Nk = G%US%RZ_to_kg_m2 / NkIce
  call get_SIS2_thermo_coefs(IST%ITV, Latent_fusion=LI, slab_ice=slab_ice)

  value = 0.0

  select case (index)

    case (ISTOCK_WATER)
      do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + kg_H * (IST%mH_ice(i,j,k) + (IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k))) * &
               IST%part_size(i,j,k) * (G%US%L_to_m**2*G%areaT(i,j)*G%mask2dT(i,j))
      enddo ; enddo ; enddo

    case (ISTOCK_HEAT)
      if (slab_ice) then
        do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
          if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
              value = value - (G%US%L_to_m**2*G%areaT(i,j)*G%mask2dT(i,j)) * IST%part_size(i,j,k) * &
                              (kg_H * IST%mH_ice(i,j,k)) * LI*G%US%Q_to_J_kg
          endif
        enddo ; enddo ; enddo
      else
        do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
          part_area = (G%US%L_to_m**2*G%areaT(i,j)*G%mask2dT(i,j)) * IST%part_size(i,j,k)
          if (part_area*IST%mH_ice(i,j,k) > 0.0) then
            value = value - (part_area * kg_H * IST%mH_snow(i,j,k)) * &
                  Energy_0degC(IST%enth_snow(i,j,k,1), IST%ITV)
            ! The pond contribution here is 0 because ponds are assumed be at 0 degC already.
            ! Otherwise add something like:
            ! value = value - (part_area * kg_H * IST%mH_pond(i,j,k)) * &
            !     Energy_0degC(enthalpy_liquid(IST%Temperature_pond(i,j,k), 0.0, IST%ITV), IST%ITV)
            do m=1,NkIce
              value = value - (part_area * (kg_H_Nk * IST%mH_ice(i,j,k))) * &
                  Energy_0degC(IST%enth_ice(i,j,k,m), IST%ITV)
            enddo
          endif
        enddo ; enddo ; enddo
      endif

    case (ISTOCK_SALT)
      !There is no salt in the snow or in the ponds.
      do m=1,NkIce ; do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + (IST%part_size(i,j,k) * (G%US%L_to_m**2*G%areaT(i,j)*G%mask2dT(i,j))) * &
            (0.001*(kg_H_Nk*IST%mH_ice(i,j,k))) * G%US%S_to_ppt*IST%sal_ice(i,j,k,m)
      enddo ; enddo ; enddo ; enddo

  end select

  if (associated(Ice%icebergs)) then
    call icebergs_stock_pe(Ice%icebergs, index, icebergs_value)
    value = value + icebergs_value
  endif

end subroutine ice_stock_pe

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Write chksums of fields in the ice data type, using interfaces that shared
!! with older sea ice models
subroutine ice_data_type_chksum(mesg, timestep, Ice, init_call)
  character(len=*),    intent(in) :: mesg      !< An identifying message
  integer         ,    intent(in) :: timestep  !< The timestep number
  type(ice_data_type), intent(in) :: Ice       !< The publicly visible ice data type.
  logical, optional,   intent(in) :: init_call !< If true omit checksums that do not make sense
                                               !! to output during initialization.

  ! Local variables
  integer(kind=int64) :: chks ! A checksum for the field
  logical :: root    ! True only on the root PE.
  logical :: init    ! If true, omit checksums that do not make sense to output
                     ! during initialization.
  integer :: outunit ! The output unit to write to.

  outunit = stdout()
  root = is_root_pe()
  init = .false. ; if (present(init_call)) init = init_call

  if (root) write(outunit,*) "BEGIN CHECKSUM(ice_data_type):: ", mesg, timestep

  if (Ice%fast_ice_PE) then
    ! These fields are only valid on fast ice PEs.
    if (.not.init) then
      chks = SIS_chksum(Ice%part_size      ) ; if (root) write(outunit,100) 'ice_data_type%part_size      ', chks
      chks = SIS_chksum(Ice%t_surf         ) ; if (root) write(outunit,100) 'ice_data_type%t_surf         ', chks
      chks = SIS_chksum(Ice%s_surf         ) ; if (root) write(outunit,100) 'ice_data_type%s_surf         ', chks
      chks = SIS_chksum(Ice%albedo         ) ; if (root) write(outunit,100) 'ice_data_type%albedo         ', chks
      chks = SIS_chksum(Ice%albedo_vis_dir ) ; if (root) write(outunit,100) 'ice_data_type%albedo_vis_dir ', chks
      chks = SIS_chksum(Ice%albedo_nir_dir ) ; if (root) write(outunit,100) 'ice_data_type%albedo_nir_dir ', chks
      chks = SIS_chksum(Ice%albedo_vis_dif ) ; if (root) write(outunit,100) 'ice_data_type%albedo_vis_dif ', chks
      chks = SIS_chksum(Ice%albedo_nir_dif ) ; if (root) write(outunit,100) 'ice_data_type%albedo_nir_dif ', chks
    endif
    chks = SIS_chksum(Ice%rough_mom   ) ; if (root) write(outunit,100)   'ice_data_type%rough_mom  ', chks
    chks = SIS_chksum(Ice%rough_heat  ) ; if (root) write(outunit,100)   'ice_data_type%rough_heat ', chks
    chks = SIS_chksum(Ice%rough_moist ) ; if (root) write(outunit,100)   'ice_data_type%rough_moist', chks

    if (.not.init) then
      chks = SIS_chksum(Ice%u_surf) ; if (root) write(outunit,100) 'ice_data_type%u_surf ', chks
      chks = SIS_chksum(Ice%v_surf) ; if (root) write(outunit,100) 'ice_data_type%v_surf ', chks
    endif

    call coupler_type_write_chksums(Ice%ocean_fields, outunit, 'ice%')
  endif

  if (Ice%slow_ice_PE) then
    ! These fields are only valid on slow ice PEs.
    chks = SIS_chksum(Ice%flux_u          ) ; if (root) write(outunit,100) 'ice_data_type%flux_u          ', chks
    chks = SIS_chksum(Ice%flux_v          ) ; if (root) write(outunit,100) 'ice_data_type%flux_v          ', chks
    chks = SIS_chksum(Ice%flux_t          ) ; if (root) write(outunit,100) 'ice_data_type%flux_t          ', chks
    chks = SIS_chksum(Ice%flux_q          ) ; if (root) write(outunit,100) 'ice_data_type%flux_q          ', chks
    chks = SIS_chksum(Ice%flux_lw         ) ; if (root) write(outunit,100) 'ice_data_type%flux_lw         ', chks
    chks = SIS_chksum(Ice%flux_sw_vis_dir ) ; if (root) write(outunit,100) 'ice_data_type%flux_sw_vis_dir ', chks
    chks = SIS_chksum(Ice%flux_sw_vis_dif ) ; if (root) write(outunit,100) 'ice_data_type%flux_sw_vis_dif ', chks
    chks = SIS_chksum(Ice%flux_sw_nir_dir ) ; if (root) write(outunit,100) 'ice_data_type%flux_sw_nir_dir ', chks
    chks = SIS_chksum(Ice%flux_sw_nir_dif ) ; if (root) write(outunit,100) 'ice_data_type%flux_sw_nir_dif ', chks
    chks = SIS_chksum(Ice%flux_lh         ) ; if (root) write(outunit,100) 'ice_data_type%flux_lh         ', chks
    chks = SIS_chksum(Ice%lprec           ) ; if (root) write(outunit,100) 'ice_data_type%lprec           ', chks
    chks = SIS_chksum(Ice%fprec           ) ; if (root) write(outunit,100) 'ice_data_type%fprec           ', chks
    chks = SIS_chksum(Ice%p_surf          ) ; if (root) write(outunit,100) 'ice_data_type%p_surf          ', chks
    chks = SIS_chksum(Ice%runoff          ) ; if (root) write(outunit,100) 'ice_data_type%runoff          ', chks
    chks = SIS_chksum(Ice%calving         ) ; if (root) write(outunit,100) 'ice_data_type%calving         ', chks
    chks = SIS_chksum(Ice%flux_salt       ) ; if (root) write(outunit,100) 'ice_data_type%flux_salt       ', chks

    if (associated(Ice%sCS)) then ; if (Ice%sCS%pass_iceberg_area_to_ocean) then
      chks = SIS_chksum(Ice%ustar_berg    ) ; if (root) write(outunit,100) 'ice_data_type%ustar_berg      ', chks
      chks = SIS_chksum(Ice%area_berg     ) ; if (root) write(outunit,100) 'ice_data_type%area_berg       ', chks
      chks = SIS_chksum(Ice%mass_berg     ) ; if (root) write(outunit,100) 'ice_data_type%mass_berg       ', chks
    endif ; endif
  endif

100 FORMAT("   CHECKSUM::",A32," = ",Z20)

end subroutine ice_data_type_chksum

end module ice_type_mod
