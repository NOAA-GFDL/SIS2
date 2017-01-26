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

use SIS2_ice_thm, only : ice_thermo_type, enth_from_TS, energy_melt_EnthS
use SIS2_ice_thm, only : get_SIS2_thermo_coefs, temp_from_En_S
use ice_bergs, only: icebergs, icebergs_stock_pe, icebergs_save_restart

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg, is_root_pe
use MOM_file_parser, only : param_file_type
use MOM_hor_index,   only : hor_index_type
use SIS_debugging,     only : chksum
use SIS_diag_mediator, only : SIS_diag_ctrl, post_data=>post_SIS_data
use SIS_diag_mediator, only : register_SIS_diag_field

use SIS_types, only : ice_state_type, fast_ice_avg_type
use SIS_ctrl_types, only : SIS_fast_CS, SIS_slow_CS

implicit none ; private

public :: ice_data_type, dealloc_ice_arrays
public :: ice_type_slow_reg_restarts, ice_type_fast_reg_restarts
public :: ice_model_restart, ice_stock_pe, ice_data_type_chksum
public :: Ice_public_type_chksum, Ice_public_type_bounds_check

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! This structure contains the ice model data (some used by calling routines);   !
! the third index is partition (1 is open water; 2... are ice cover by category)!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
type ice_data_type !  ice_public_type
  type(domain2D)          :: Domain         ! A copy of the fast ice domain without halos.
  type(domain2D)          :: slow_Domain_NH ! A copy of the slow ice domain without halos.
  type(domain2D), pointer :: &
    fast_domain => NULL(), & ! A pointer to the fast ice mpp domain or a copy
                             ! on slow ice PEs.
    slow_domain => NULL()    ! A pointer to the fast ice mpp domain or a copy
                             ! on slow ice PEs.
  type(time_type)                  :: Time
  logical                          :: pe
  logical                          :: slow_ice_pe = .false.
  logical                          :: fast_ice_pe = .false.
  logical                          :: shared_slow_fast_PEs = .true.
  integer                          :: xtype
  integer, pointer, dimension(:)   :: slow_pelist =>NULL() ! Used for flux-exchange with slow processes.
  integer, pointer, dimension(:)   :: fast_pelist =>NULL() ! Used for flux-exchange with fast processes.
  integer, pointer, dimension(:)   :: pelist   =>NULL() ! Used for flux-exchange.
  logical, pointer, dimension(:,:) :: ocean_pt =>NULL() ! An array that indicates all ocean points as true.

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
    SST_C => NULL(), &    ! The ocean surface temperature, in deg C.
    flux_u => NULL(), &   ! The flux of x-momentum into the ocean, in Pa.
    flux_v => NULL(), &   ! The flux of y-momentum into the ocean, in Pa.
    flux_t => NULL(), &   ! The flux of sensible heat out of the ocean, in W m-2.
    flux_q => NULL(), &   ! The evaporative moisture flux out of the ocean, in kg m-2 s-1.
    flux_lw => NULL(), &  ! The longwave flux out of the ocean, in W m-2.
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
    ustar_berg => NULL(), &  !ustar contribution below icebergs in m/s
    area_berg => NULL(),  &  !fraction of grid cell covered by icebergs in m2/m2
    mass_berg => NULL(),  &  !mass of icebergs in km/m^2
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
             ! h_ice outside the ice module.
  integer, dimension(3)    :: axes
  type(coupler_3d_bc_type) :: ocean_fields ! array of fields used for additional tracers
                                           ! whose surface state is shared with the atmosphere.
  type(coupler_2d_bc_type) :: ocean_fluxes ! array of fluxes from the ice to the ocean used
                                           ! for additional tracers
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

  ! The following are actually private to SIS2, and are not used elsewhere by
  ! other FMS modules.
  type(icebergs),    pointer :: icebergs => NULL()
  type(SIS_fast_CS), pointer :: fCS => NULL()
  type(SIS_slow_CS), pointer :: sCS => NULL()
  type(restart_file_type), pointer :: Ice_restart => NULL()
  type(restart_file_type), pointer :: Ice_fast_restart => NULL()
end type ice_data_type !  ice_public_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_type_slow_reg_restarts - allocate the arrays in the ice_data_type        !
!     that are predominantly associated with the slow processors, and register !
!     any variables in the ice data type that need to be included in the slow  !
!     ice restart files.                                                       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_type_slow_reg_restarts(domain, CatIce, param_file, Ice, &
                                           Ice_restart, restart_file)
  type(domain2d),          intent(in)    :: domain
  integer,                 intent(in)    :: CatIce
  type(param_file_type),   intent(in)    :: param_file
  type(ice_data_type),     intent(inout) :: Ice
  type(restart_file_type), pointer       :: Ice_restart
  character(len=*),        intent(in)    :: restart_file

  ! This subroutine allocates the externally visible ice_data_type's arrays and
  ! registers the appopriate ones for inclusion in the restart file.
  integer :: isc, iec, jsc, jec, km, idr

  call mpp_get_compute_domain(domain, isc, iec, jsc, jec )
  km = CatIce + 1

  ! The fields t_surf, s_surf, and part_size are only available on fast PEs.

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

  allocate(Ice%SST_C(isc:iec, jsc:jec)) ; Ice%SST_C(:,:) = 0.0
  allocate(Ice%area(isc:iec, jsc:jec)) ; Ice%area(:,:) = 0.0
  allocate(Ice%mi(isc:iec, jsc:jec)) ; Ice%mi(:,:) = 0.0 !NR

  if (associated(Ice%sCS)) then ; if (Ice%sCS%pass_iceberg_area_to_ocean) then
    allocate(Ice%ustar_berg(isc:iec, jsc:jec)) ; Ice%ustar_berg(:,:) = 0.0
    allocate(Ice%area_berg(isc:iec, jsc:jec)) ; Ice%area_berg(:,:) = 0.0
    allocate(Ice%mass_berg(isc:iec, jsc:jec)) ; Ice%mass_berg(:,:) = 0.0
  endif ; endif

  ! These are used by the ocean model, and need to be in the slow PE restarts.
  if (associated(Ice_restart)) then
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
  endif
end subroutine ice_type_slow_reg_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_type_slow_reg_restarts - allocate the arrays in the ice_data_type        !
!     that are predominantly associated with the fast processors, and register !
!     any variables in the ice data type that need to be included in the fast  !
!     ice restart files.                                                       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_type_fast_reg_restarts(domain, CatIce, param_file, Ice, &
                                           Ice_restart, restart_file)
  type(domain2d),          intent(in)    :: domain
  integer,                 intent(in)    :: CatIce
  type(param_file_type),   intent(in)    :: param_file
  type(ice_data_type),     intent(inout) :: Ice
  type(restart_file_type), pointer       :: Ice_restart
  character(len=*),        intent(in)    :: restart_file

  ! This subroutine allocates the externally visible ice_data_type's arrays and
  ! registers the appopriate ones for inclusion in the restart file.
  integer :: isc, iec, jsc, jec, km, idr

  call mpp_get_compute_domain(domain, isc, iec, jsc, jec )
  km = CatIce + 1

  allocate(Ice%t_surf(isc:iec, jsc:jec, km)) ; Ice%t_surf(:,:,:) = 0.0
  allocate(Ice%s_surf(isc:iec, jsc:jec)) ; Ice%s_surf(:,:) = 0.0
  allocate(Ice%part_size(isc:iec, jsc:jec, km)) ; Ice%part_size(:,:,:) = 0.0

  allocate(Ice%u_surf(isc:iec, jsc:jec, km)) ; Ice%u_surf(:,:,:) = 0.0
  allocate(Ice%v_surf(isc:iec, jsc:jec, km)) ; Ice%v_surf(:,:,:) = 0.0
  allocate(Ice%ocean_pt(isc:iec, jsc:jec)) ; Ice%ocean_pt(:,:) = .false. !derived

  allocate(Ice%rough_mom(isc:iec, jsc:jec, km)) ; Ice%rough_mom(:,:,:) = 0.0
  allocate(Ice%rough_heat(isc:iec, jsc:jec, km)) ; Ice%rough_heat(:,:,:) = 0.0
  allocate(Ice%rough_moist(isc:iec, jsc:jec, km)) ; Ice%rough_moist(:,:,:) = 0.0

  allocate(Ice%albedo(isc:iec, jsc:jec, km)) ; Ice%albedo(:,:,:) = 0.0  ! Derived?
  allocate(Ice%albedo_vis_dir(isc:iec, jsc:jec, km)) ; Ice%albedo_vis_dir(:,:,:) = 0.0
  allocate(Ice%albedo_nir_dir(isc:iec, jsc:jec, km)) ; Ice%albedo_nir_dir(:,:,:) = 0.0
  allocate(Ice%albedo_vis_dif(isc:iec, jsc:jec, km)) ; Ice%albedo_vis_dif(:,:,:) = 0.0
  allocate(Ice%albedo_nir_dif(isc:iec, jsc:jec, km)) ; Ice%albedo_nir_dif(:,:,:) = 0.0

  ! Now register some of these arrays to be read from the restart files.
  ! These are used by the atmospheric model, and need to be in the fast PE restarts.
  if (associated(Ice_restart)) then
    idr = register_restart_field(Ice_restart, restart_file, 'rough_mom',   Ice%rough_mom,   domain=domain)
    idr = register_restart_field(Ice_restart, restart_file, 'rough_heat',  Ice%rough_heat,  domain=domain)
    idr = register_restart_field(Ice_restart, restart_file, 'rough_moist', Ice%rough_moist, domain=domain)
  endif

end subroutine ice_type_fast_reg_restarts


subroutine dealloc_Ice_arrays(Ice)
  type(ice_data_type), intent(inout) :: Ice

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

subroutine Ice_public_type_chksum(mesg, Ice, check_fast, check_slow)
  character(len=*),    intent(in) :: mesg
  type(ice_data_type), intent(in) :: Ice
  logical, optional,   intent(in) :: check_fast, check_slow
!   This subroutine writes out chksums for the model's basic state variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      Ice - An ice_data_type structure whose elements are to be
!                  checksummed.

  ! Note that the publicly visible ice_data_type has no halos, so it is not
  ! possible do check their values.

  logical :: fast_fields, slow_fields

  fast_fields = Ice%fast_ice_PE ; slow_fields = Ice%slow_ice_PE

  if (present(check_fast)) then
    fast_fields = Ice%fast_ice_PE .and. check_fast
    if (.not.present(check_slow)) slow_fields = .false.
  elseif (present(check_slow)) then
    slow_fields = Ice%slow_ice_PE .and. check_slow
  endif

  ! These fields are on all PEs.
  if (fast_fields .or. slow_fields) &
  call chksum(Ice%part_size, trim(mesg)//" Ice%part_size")

  if (fast_fields) then ! This is a fast-ice PE.
    call chksum(Ice%albedo, trim(mesg)//" Ice%albedo")
    call chksum(Ice%albedo_vis_dir, trim(mesg)//" Ice%albedo_vis_dir")
    call chksum(Ice%albedo_nir_dir, trim(mesg)//" Ice%albedo_nir_dir")
    call chksum(Ice%albedo_vis_dif, trim(mesg)//" Ice%albedo_vis_dif")
    call chksum(Ice%albedo_nir_dif, trim(mesg)//" Ice%albedo_nir_dif")
    call chksum(Ice%rough_mom, trim(mesg)//" Ice%rough_mom")
    call chksum(Ice%rough_heat, trim(mesg)//" Ice%rough_heat")
    call chksum(Ice%rough_moist, trim(mesg)//" Ice%rough_moist")

    call chksum(Ice%t_surf, trim(mesg)//" Ice%t_surf")
    call chksum(Ice%s_surf, trim(mesg)//" Ice%s_surf")
    call chksum(Ice%u_surf, trim(mesg)//" Ice%u_surf")
    call chksum(Ice%v_surf, trim(mesg)//" Ice%v_surf")
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
  endif

  if (slow_fields .and. associated(Ice%sCS)) then ; if (Ice%sCS%pass_iceberg_area_to_ocean) then
    call chksum(Ice%ustar_berg, trim(mesg)//" Ice%ustar_berg")
    call chksum(Ice%area_berg, trim(mesg)//" Ice%area_berg")
    call chksum(Ice%mass_berg, trim(mesg)//" Ice%mass_berg")
  endif ; endif
end subroutine Ice_public_type_chksum

subroutine Ice_public_type_bounds_check(Ice, G, msg)
  type(ice_data_type),     intent(in)    :: Ice
  type(SIS_hor_grid_type), intent(inout) :: G
  character(len=*),        intent(in)    :: msg

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
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    if ((Ice%t_surf(i2,j2,k2) < t_min) .or. (Ice%t_surf(i2,j2,k2) > t_max)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; k_bad = k ; endif
    endif
  enddo ; enddo ; enddo
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    if ((Ice%s_surf(i2,j2) < 0.0) .or. (Ice%s_surf(i2,j2) > 100.0)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; endif
    endif
  enddo ; enddo
  if (fluxes_avail) then ; do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    if ((abs(Ice%flux_t(i2,j2)) > 1e4) .or. (abs(Ice%flux_lw(i2,j2)) > 1e4)) then
      n_bad = n_bad + 1
      if (n_bad == 1) then ; i_bad = i ; j_bad = j ; endif
    endif
  enddo ; enddo ; endif

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
subroutine ice_model_restart(Ice, time_stamp)
  type(ice_data_type), intent(inout) :: Ice
  character(len=*),    intent(in), optional :: time_stamp

  if (associated(Ice%Ice_restart)) then
    call save_restart(Ice%Ice_restart, time_stamp)
    if (associated(Ice%Ice_fast_restart)) then
      if (.not.associated(Ice%Ice_fast_restart,Ice%Ice_restart)) &
        call save_restart(Ice%Ice_fast_restart, time_stamp)
    endif
  elseif (associated(Ice%Ice_fast_restart)) then
    call save_restart(Ice%Ice_fast_restart, time_stamp)
  endif
  call icebergs_save_restart(Ice%icebergs)

end subroutine ice_model_restart
! </SUBROUTINE>
!=======================================================================

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_stock_pe - returns stocks of heat, water, etc. for conservation checks   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_stock_pe(Ice, index, value)

  use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT

  type(ice_data_type) :: Ice
  integer, intent(in) :: index
  real, intent(out)   :: value
  type(ice_state_type), pointer :: IST => NULL()

  real :: icebergs_value
  real :: LI
  real :: part_wt, I_NkIce, kg_H, kg_H_Nk
  integer :: i, j, k, m, isc, iec, jsc, jec, ncat, NkIce
  logical :: slab_ice    ! If true, use the very old slab ice thermodynamics,
                         ! with effectively zero heat capacity of ice and snow.
  type(SIS_hor_grid_type), pointer :: G => NULL()

  value = 0.0
  if(.not.Ice%pe) return

  if (associated(Ice%sCS)) then
    IST => Ice%sCS%IST
    G => Ice%sCS%G
    ncat = Ice%sCS%IG%CatIce ; NkIce = Ice%sCS%IG%NkIce ; kg_H = Ice%sCS%IG%H_to_kg_m2
  elseif (associated(Ice%fCS)) then
    IST => Ice%fCS%IST
    G => Ice%fCS%G
    ncat = Ice%fCS%IG%CatIce ; NkIce = Ice%fCS%IG%NkIce ; kg_H = Ice%fCS%IG%H_to_kg_m2
  else
    call SIS_error(WARNING, "ice_stock_pe called with an ice_data_type "//&
                   "without either sCS or fCS associated.")
    return
  endif

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  I_NkIce = 1.0 / NkIce  ; kg_H_Nk = kg_H / NkIce
  call get_SIS2_thermo_coefs(IST%ITV, Latent_fusion=LI, slab_ice=slab_ice)

  select case (index)

    case (ISTOCK_WATER)
      value = 0.0
      do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + kg_H * (IST%mH_ice(i,j,k) + IST%mH_snow(i,j,k)) * &
               IST%part_size(i,j,k) * (G%areaT(i,j)*G%mask2dT(i,j))
      enddo ; enddo ; enddo

    case (ISTOCK_HEAT)
      value = 0.0
      if (slab_ice) then
        do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
          if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
              value = value - (G%areaT(i,j)*G%mask2dT(i,j)) * IST%part_size(i,j,k) * &
                              (kg_H * IST%mH_ice(i,j,k)) * LI
          endif
        enddo ; enddo ; enddo
      else
        do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
          part_wt = (G%areaT(i,j)*G%mask2dT(i,j)) * IST%part_size(i,j,k)
          if (part_wt*IST%mH_ice(i,j,k) > 0.0) then
            value = value - (part_wt * (kg_H * IST%mH_snow(i,j,k))) * &
                Energy_melt_enthS(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
            do m=1,NkIce
              value = value - (part_wt * (kg_H_Nk * IST%mH_ice(i,j,k))) * &
                  Energy_melt_enthS(IST%enth_ice(i,j,k,m), IST%sal_ice(i,j,k,m), IST%ITV)
            enddo
          endif
        enddo ; enddo ; enddo
      endif

    case (ISTOCK_SALT)
      !There is no salt in the snow.
      value = 0.0
      do m=1,NkIce ; do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
        value = value + (IST%part_size(i,j,k) * (G%areaT(i,j)*G%mask2dT(i,j))) * &
            (0.001*(kg_H_Nk*IST%mH_ice(i,j,k))) * IST%sal_ice(i,j,k,m)
      enddo ; enddo ; enddo ; enddo

    case default

      value = 0.0

  end select

  if (associated(Ice%icebergs)) then
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
  ! These fields are on all PEs.
  write(outunit,100) 'ice_data_type%part_size          ',mpp_chksum(Ice%part_size         )
  write(outunit,100) 'ice_data_type%t_surf             ',mpp_chksum(Ice%t_surf            )
  write(outunit,100) 'ice_data_type%s_surf             ',mpp_chksum(Ice%s_surf            )

  if (Ice%fast_ice_PE) then
    ! These fields are only valid on fast ice PEs.
    write(outunit,100) 'ice_data_type%albedo             ',mpp_chksum(Ice%albedo          )
    write(outunit,100) 'ice_data_type%albedo_vis_dir     ',mpp_chksum(Ice%albedo_vis_dir  )
    write(outunit,100) 'ice_data_type%albedo_nir_dir     ',mpp_chksum(Ice%albedo_nir_dir  )
    write(outunit,100) 'ice_data_type%albedo_vis_dif     ',mpp_chksum(Ice%albedo_vis_dif  )
    write(outunit,100) 'ice_data_type%albedo_nir_dif     ',mpp_chksum(Ice%albedo_nir_dif  )
    write(outunit,100) 'ice_data_type%rough_mom          ',mpp_chksum(Ice%rough_mom       )
    write(outunit,100) 'ice_data_type%rough_heat         ',mpp_chksum(Ice%rough_heat      )
    write(outunit,100) 'ice_data_type%rough_moist        ',mpp_chksum(Ice%rough_moist     )

    write(outunit,100) 'ice_data_type%u_surf             ',mpp_chksum(Ice%u_surf          )
    write(outunit,100) 'ice_data_type%v_surf             ',mpp_chksum(Ice%v_surf          )

    do n=1,Ice%ocean_fields%num_bcs ; do m=1,Ice%ocean_fields%bc(n)%num_fields
      write(outunit,101) 'ice%', trim(Ice%ocean_fields%bc(n)%name), &
                         trim(Ice%ocean_fields%bc(n)%field(m)%name), &
                         mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
    enddo ; enddo
  endif

  if (Ice%slow_ice_PE) then
    ! These fields are only valid on slow ice PEs.
    write(outunit,100) 'ice_data_type%flux_u             ',mpp_chksum(Ice%flux_u          )
    write(outunit,100) 'ice_data_type%flux_v             ',mpp_chksum(Ice%flux_v          )
    write(outunit,100) 'ice_data_type%flux_t             ',mpp_chksum(Ice%flux_t          )
    write(outunit,100) 'ice_data_type%flux_q             ',mpp_chksum(Ice%flux_q          )
    write(outunit,100) 'ice_data_type%flux_lw            ',mpp_chksum(Ice%flux_lw         )
    write(outunit,100) 'ice_data_type%flux_sw_vis_dir    ',mpp_chksum(Ice%flux_sw_vis_dir )
    write(outunit,100) 'ice_data_type%flux_sw_vis_dif    ',mpp_chksum(Ice%flux_sw_vis_dif )
    write(outunit,100) 'ice_data_type%flux_sw_nir_dir    ',mpp_chksum(Ice%flux_sw_nir_dir )
    write(outunit,100) 'ice_data_type%flux_sw_nir_dif    ',mpp_chksum(Ice%flux_sw_nir_dif )
    write(outunit,100) 'ice_data_type%flux_lh            ',mpp_chksum(Ice%flux_lh         )
    write(outunit,100) 'ice_data_type%lprec              ',mpp_chksum(Ice%lprec           )
    write(outunit,100) 'ice_data_type%fprec              ',mpp_chksum(Ice%fprec           )
    write(outunit,100) 'ice_data_type%p_surf             ',mpp_chksum(Ice%p_surf          )
    write(outunit,100) 'ice_data_type%runoff             ',mpp_chksum(Ice%runoff          )
    write(outunit,100) 'ice_data_type%calving            ',mpp_chksum(Ice%calving         )
    write(outunit,100) 'ice_data_type%flux_salt          ',mpp_chksum(Ice%flux_salt       )

    if (associated(Ice%sCS)) then ; if (Ice%sCS%pass_iceberg_area_to_ocean) then
      write(outunit,100) 'ice_data_type%ustar_berg         ',mpp_chksum(Ice%ustar_berg    )
      write(outunit,100) 'ice_data_type%area_berg          ',mpp_chksum(Ice%area_berg     )
      write(outunit,100) 'ice_data_type%mass_berg          ',mpp_chksum(Ice%mass_berg     )
    endif ; endif

  endif

100 FORMAT("   CHECKSUM::",A32," = ",Z20)
101 FORMAT("   CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ice_data_type_chksum

end module ice_type_mod
