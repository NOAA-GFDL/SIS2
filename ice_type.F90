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
use SIS_dyn_trans, only : dyn_trans_CS
use SIS_slow_thermo, only : slow_thermo_CS
use SIS_sum_output_type, only : SIS_sum_out_CS
use SIS_tracer_registry, only : SIS_tracer_registry_type
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS

use SIS_types, only : ice_state_type, fast_ice_avg_type, ice_rad_type
use SIS_types, only : ice_ocean_flux_type, ocean_sfc_state_type
use SIS_types, only : fast_thermo_CS

implicit none ; private

#include <SIS2_memory.h>

public :: ice_data_type, ice_data_type_register_restarts, dealloc_ice_arrays 
public :: ice_model_restart, ice_stock_pe
public :: ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum, ice_data_type_chksum
public :: Ice_public_type_chksum, Ice_public_type_bounds_check

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

  ! The following are actually private to SIS2, and are not used elsewhere by
  ! other FMS modules.
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
  type(ice_rad_type), pointer :: Rad => NULL()    ! A structure with fields related to
                               ! the absorption, reflection and transmission of
                               ! shortwave radiation.
  type(fast_thermo_CS), pointer :: fast_thermo_CSp => NULL()
  type(slow_thermo_CS), pointer :: slow_thermo_CSp => NULL()
  type(dyn_trans_CS), pointer :: dyn_trans_CSp => NULL()
  type(SIS_tracer_flow_control_CS), pointer :: SIS_tracer_flow_CSp => NULL()

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
