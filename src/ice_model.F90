!> This is the central module for the SIS2 sea ice model.
module ice_model_mod

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   SIS2 is a SEA ICE MODEL for coupling through the GFDL exchange grid. SIS2  !
! is a revision of the original SIS with have extended capabilities, including !
! the option of using a B-grid or C-grid spatial discretization.  The SIS2     !
! software has been extensively reformulated from SIS for greater consistency  !
! with the Modular Ocean Model, version 6 (MOM6), and to permit might tighter  !
! dynamical coupling between the ocean and sea-ice.                            !
!   This module manages fluxes between sub-modules, many diagnostics, and the  !
! overall time stepping of the sea ice. Sea ice dynamics are handled in        !
! ice_dyn_bgrid.F90 or ice_dyn_cgrid.F90, while the transport of mass, heat,   !
! and tracers occurs in ice_transport.F90.  Sea ice thermodynamics is treated  !
! in ice_thm.F90 and other modules that are subsequently called from there.    !
! The Lagrangian icebergs code of Adcroft and Martin is called from SIS.       !
!   The original SIS was developed by Mike Winton (Michael.Winton@noaa.gov).   !
! SIS2 has been developed by Robert Hallberg and Mike Winton, with             !
! contributions from many people at NOAA/GFDL, including Alistair Adcroft and  !
! Niki Zadeh.                                                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use MOM_domains,       only : MOM_domain_type
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : fill_symmetric_edges, MOM_domains_init, clone_MOM_domain
use MOM_dyn_horgrid,   only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_param, log_version, read_param, param_file_type
use MOM_file_parser,   only : open_param_file, close_param_file
use MOM_hor_index,     only : hor_index_type, hor_index_init
use MOM_io,            only : file_exists
use MOM_obsolete_params, only : obsolete_logical, obsolete_real
use MOM_string_functions, only : uppercase
use MOM_time_manager,  only : time_type, time_type_to_real, real_to_time
use MOM_time_manager,  only : operator(+), operator(-)
use MOM_time_manager,  only : operator(>), operator(*), operator(/), operator(/=)
use MOM_unit_scaling,  only : unit_scale_type, unit_scaling_init
use MOM_unit_scaling,  only : unit_scaling_end, fix_restart_unit_scaling

use astronomy_mod, only : astronomy_init, astronomy_end
use astronomy_mod, only : universal_time, orbital_time, diurnal_solar, daily_mean_solar
use ocean_albedo_mod, only : compute_ocean_albedo            ! ice sets ocean surface
use ocean_rough_mod,  only : compute_ocean_roughness         ! properties over water

use ice_bergs,          only : icebergs, icebergs_run, icebergs_init, icebergs_end
use ice_boundary_types, only : ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
use ice_boundary_types, only : ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
use ice_boundary_types, only : lnd_ice_bnd_type_chksum
use ice_grid,           only : set_ice_grid, ice_grid_end, ice_grid_type
use ice_spec_mod,       only : get_sea_surface
use ice_type_mod,       only : ice_data_type, dealloc_ice_arrays
use ice_type_mod,       only : ice_type_slow_reg_restarts, ice_type_fast_reg_restarts
use ice_type_mod,       only : Ice_public_type_chksum, Ice_public_type_bounds_check
use ice_type_mod,       only : ice_model_restart, ice_stock_pe, ice_data_type_chksum

use SIS_ctrl_types,    only : SIS_slow_CS, SIS_fast_CS
use SIS_ctrl_types,    only : ice_diagnostics_init, ice_diags_fast_init
use SIS_debugging,     only : chksum, uvchksum, Bchksum, SIS_debugging_init
use SIS_diag_mediator, only : set_SIS_axes_info, SIS_diag_mediator_init, SIS_diag_mediator_end
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_dyn_trans,     only : SIS_dynamics_trans, SIS_multi_dyn_trans, update_icebergs
use SIS_dyn_trans,     only : slab_ice_dyn_trans
use SIS_dyn_trans,     only : SIS_dyn_trans_register_restarts, SIS_dyn_trans_init, SIS_dyn_trans_end
use SIS_dyn_trans,     only : SIS_dyn_trans_read_alt_restarts, stresses_to_stress_mag
use SIS_dyn_trans,     only : SIS_dyn_trans_transport_CS, SIS_dyn_trans_sum_output_CS
use SIS_fast_thermo,   only : accumulate_deposition_fluxes, convert_frost_to_snow
use SIS_fast_thermo,   only : do_update_ice_model_fast, avg_top_quantities, total_top_quantities
use SIS_fast_thermo,   only : redo_update_ice_model_fast, find_excess_fluxes
use SIS_fast_thermo,   only : infill_array, SIS_fast_thermo_init, SIS_fast_thermo_end
use SIS_framework,     only : set_domain, nullify_domain, broadcast_domain
use SIS_restart,       only : restore_SIS_state, query_initialized=>query_inited, SIS_restart_init
use SIS_restart,       only : determine_is_new_run, is_new_run
use SIS_framework,     only : coupler_1d_bc_type, coupler_2d_bc_type, coupler_3d_bc_type
use SIS_framework,     only : coupler_type_spawn, coupler_type_initialized
use SIS_framework,     only : coupler_type_rescale_data, coupler_type_copy_data
use SIS_fixed_initialization, only : SIS_initialize_fixed
use SIS_get_input,     only : Get_SIS_input, directories
use SIS_hor_grid,      only : SIS_hor_grid_type, set_hor_grid, SIS_hor_grid_end, set_first_direction
use SIS_open_boundary, only : ice_OBC_type
use SIS_optics,        only : ice_optics_SIS2, SIS_optics_init, SIS_optics_end, SIS_optics_CS
use SIS_optics,        only : VIS_DIR, VIS_DIF, NIR_DIR, NIR_DIF
use SIS_slow_thermo,   only : slow_thermodynamics, SIS_slow_thermo_init, SIS_slow_thermo_end
use SIS_slow_thermo,   only : SIS_slow_thermo_set_ptrs
use SIS_state_initialization, only : read_archaic_thermo_restarts, initialize_ice_categories
use SIS_state_initialization, only : ice_state_mass_init, ice_state_thermo_init
use SIS_sum_output,    only : SIS_sum_output_init,  write_ice_statistics
use SIS_tracer_flow_control, only : SIS_call_tracer_register, SIS_tracer_flow_control_init
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_end
use SIS_tracer_registry, only : register_SIS_tracer, register_SIS_tracer_pair
use SIS_transcribe_grid, only : copy_dyngrid_to_SIS_horgrid, copy_SIS_horgrid_to_dyngrid
use SIS_transport,     only : adjust_ice_categories
use SIS_types,         only : ice_ocean_flux_type, alloc_ice_ocean_flux, dealloc_ice_ocean_flux
use SIS_types,         only : ocean_sfc_state_type, alloc_ocean_sfc_state, dealloc_ocean_sfc_state, OSS_chksum
use SIS_types,         only : fast_ice_avg_type, alloc_fast_ice_avg, dealloc_fast_ice_avg
use SIS_types,         only : total_sfc_flux_type, alloc_total_sfc_flux, dealloc_total_sfc_flux
use SIS_types,         only : ice_rad_type, ice_rad_register_restarts, dealloc_ice_rad, alloc_ice_rad
use SIS_types,         only : simple_OSS_type, alloc_simple_OSS, dealloc_simple_OSS
use SIS_types,         only : ice_state_type, alloc_IST_arrays, dealloc_IST_arrays
use SIS_types,         only : IST_chksum, IST_bounds_check, ice_state_register_restarts
use SIS_types,         only : ice_state_read_alt_restarts, register_fast_to_slow_restarts
use SIS_types,         only : register_unit_conversion_restarts
use SIS_types,         only : rescale_fast_to_slow_restart_fields, rescale_ice_state_restart_fields
use SIS_types,         only : copy_IST_to_IST, copy_FIA_to_FIA, copy_sOSS_to_sOSS
use SIS_types,         only : copy_TSF_to_TSF, redistribute_TSF_to_TSF, TSF_chksum
use SIS_types,         only : copy_Rad_to_Rad, redistribute_Rad_to_Rad
use SIS_types,         only : redistribute_IST_to_IST, redistribute_FIA_to_FIA
use SIS_types,         only : redistribute_sOSS_to_sOSS, FIA_chksum, IOF_chksum, translate_OSS_to_sOSS
use SIS_utils,         only : post_avg, ice_grid_chksum
use SIS2_ice_thm,      only : ice_temp_SIS2, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,      only : ice_thermo_init, ice_thermo_end, T_freeze, ice_thermo_type
use specified_ice,     only : specified_ice_dynamics, specified_ice_init, specified_ice_CS
use specified_ice,     only : specified_ice_end, specified_ice_sum_output_CS

implicit none ; private

#include <SIS2_memory.h>

public :: ice_data_type, ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
public :: ice_model_init, share_ice_domains, ice_model_end, ice_stock_pe
public :: update_ice_model_fast
public :: ice_model_restart  ! for intermediate restarts
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum, ice_data_type_chksum
public :: update_ice_atm_deposition_flux
public :: unpack_ocean_ice_boundary, exchange_slow_to_fast_ice, set_ice_surface_fields
public :: ice_model_fast_cleanup, unpack_land_ice_boundary
public :: exchange_fast_to_slow_ice, update_ice_model_slow
public :: update_ice_slow_thermo, update_ice_dynamics_trans

!>@{ CPU time clock IDs
integer :: iceClock
integer :: ice_clock_slow, ice_clock_fast, ice_clock_exchange
!!@}

integer, parameter :: REDIST=2 !< Redistribute for exchange
integer, parameter :: DIRECT=3 !< Use direct exchange

contains


!-----------------------------------------------------------------------
!> Update the sea-ice state due to slow processes, including dynamics,
!! freezing and melting, precipitation, and transport.
subroutine update_ice_model_slow(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type.

  call update_ice_slow_thermo(Ice)

  call update_ice_dynamics_trans(Ice)

end subroutine update_ice_model_slow


!-----------------------------------------------------------------------
!> Update the sea-ice state due to slow thermodynamic processes, including
!! freezing and melting, precipitation, and brine drainage, and possibly also
!! also the accumulated effects of faster thermodynamic processes
subroutine update_ice_slow_thermo(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type.

  ! These pointers are used to simplify the code below.
  type(ice_grid_type),     pointer :: sIG => NULL()
  type(SIS_hor_grid_type), pointer :: sG => NULL()
  type(ice_state_type),    pointer :: sIST => NULL()
  type(fast_ice_avg_type), pointer :: FIA => NULL()
  type(ice_rad_type),      pointer :: Rad => NULL()
  type(unit_scale_type),   pointer :: US => NULL()
  real :: dt_slow  ! The time step over which to advance the model [T ~> s].
  integer :: i, j, i2, j2, i_off, j_off

  if (.not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "The pointer to Ice%sCS must be associated in update_ice_slow_thermo.")

  sIST => Ice%sCS%IST ; sIG => Ice%sCS%IG ; sG => Ice%sCS%G ; FIA => Ice%sCS%FIA
  Rad => Ice%sCS%Rad ; US => Ice%sCS%US
  call cpu_clock_begin(iceClock) ; call cpu_clock_begin(ice_clock_slow)

  ! Advance the slow PE clock to give the end time of the slow timestep.  There
  ! is a separate clock inside the fCS that is advanced elsewhere.
  Ice%sCS%Time = Ice%sCS%Time + Ice%sCS%Time_step_slow
  if (.not.associated(Ice%fCS)) then
    Ice%Time = Ice%sCS%Time
  endif
  dt_slow = US%s_to_T*time_type_to_real(Ice%sCS%Time_step_slow)

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Start update_ice_slow_thermo", Ice, check_slow=.true.)
    call FIA_chksum("Start update_ice_slow_thermo", FIA, sG, US)
    ! call IOF_chksum("Start update_ice_slow_thermo", Ice%sCS%IOF, sG, US)
  endif

  ! Store some diagnostic fluxes...
  !$OMP parallel do default(none) shared(sG, FIA)
  do j=sG%jsc,sG%jec ; do i=sG%isc,sG%iec
    FIA%calving_preberg(i,j) = FIA%calving(i,j)
    FIA%calving_hflx_preberg(i,j) = FIA%calving_hflx(i,j)
  enddo ; enddo

  if (Ice%sCS%redo_fast_update) then
    call redo_update_ice_model_fast(sIST, Ice%sCS%sOSS, Ice%sCS%Rad, FIA, Ice%sCS%TSF, &
              Ice%sCS%optics_CSp, Ice%sCS%Time_step_slow, Ice%sCS%fast_thermo_CSp, sG, US, sIG)

    call find_excess_fluxes(FIA, Ice%sCS%TSF, Ice%sCS%XSF, sIST%part_size, sG, US, sIG)
  endif

  call convert_frost_to_snow(FIA, sG, US, sIG)

  if (Ice%sCS%do_icebergs) then
    if (Ice%sCS%berg_windstress_bug) then
      ! This code is only required to reproduce an old bug.
      i_off = LBOUND(Ice%flux_t,1) - sG%isc
      j_off = LBOUND(Ice%flux_t,2) - sG%jsc
      !$OMP parallel do default(none) shared(Ice,sG,US,i_off,j_off) private(i2,j2)
      do j=sG%jsc,sG%jec ; do i=sG%isc,sG%iec
        i2 = i+i_off ; j2 = j+j_off
        Ice%sCS%IOF%flux_u_ocn(i,j) = US%kg_m2s_to_RZ_T*US%m_s_to_L_T*Ice%flux_u(i2,j2)
        Ice%sCS%IOF%flux_v_ocn(i,j) = US%kg_m2s_to_RZ_T*US%m_s_to_L_T*Ice%flux_v(i2,j2)
      enddo ; enddo
    endif

    call cpu_clock_end(ice_clock_slow) ; call cpu_clock_end(iceClock)
    call update_icebergs(sIST, Ice%sCS%OSS, Ice%sCS%IOF, FIA, Ice%icebergs, US%T_to_s*dt_slow, &
                         sG, US, sIG, Ice%sCS%dyn_trans_CSp)
    call cpu_clock_begin(iceClock) ; call cpu_clock_begin(ice_clock_slow)

    if (Ice%sCS%debug) then
      call FIA_chksum("After update_icebergs", FIA, sG, US)
    endif
  endif

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Before slow_thermodynamics", Ice, check_slow=.true.)
    call FIA_chksum("Before slow_thermodynamics", FIA, sG, US)
    call IST_chksum("Before slow_thermodynamics", sIST, sG, US, sIG)
    call OSS_chksum("Before slow_thermodynamics", Ice%sCS%OSS, sG, US)
    if (associated(Ice%sCS%XSF)) &
      call TSF_chksum("Before slow_thermodynamics XSF", Ice%sCS%XSF, sG, US)
    ! call IOF_chksum("Before slow_thermodynamics", Ice%sCS%IOF, sG, US)
  endif

  call slow_thermodynamics(sIST, dt_slow, Ice%sCS%slow_thermo_CSp, Ice%sCS%OSS, FIA, &
                           Ice%sCS%XSF, Ice%sCS%IOF, sG, US, sIG)
  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Before set_ocean_top_fluxes", Ice, check_slow=.true.)
    call IOF_chksum("Before set_ocean_top_fluxes", Ice%sCS%IOF, sG, US, thermo_fluxes=.true.)
    call IST_chksum("Before set_ocean_top_fluxes", sIST, sG, US, sIG)
  endif
  ! Set up the thermodynamic fluxes in the externally visible structure Ice.
  call set_ocean_top_fluxes(Ice, sIST, Ice%sCS%IOF, FIA, Ice%sCS%OSS, sG, US, sIG, Ice%sCS)

  call cpu_clock_end(ice_clock_slow) ; call cpu_clock_end(iceClock)

end subroutine update_ice_slow_thermo

!-----------------------------------------------------------------------
!> Update the sea-ice state due to dynamics and ice transport.
subroutine update_ice_dynamics_trans(Ice, time_step, start_cycle, end_cycle, cycle_length)
  type(ice_data_type),       intent(inout) :: Ice !< The publicly visible ice data type.
  type(time_type), optional, intent(in)    :: time_step !< The amount of time to cover in this update.
  logical,         optional, intent(in)    :: start_cycle !< This indicates whether this call is to be
                                                  !! treated as the first call to update_ice_dynamics_trans
                                                  !! in a time-stepping cycle; missing is like true.
  logical,         optional, intent(in)    :: end_cycle   !< This indicates whether this call is to be
                                                  !! treated as the last call to update_ice_dynamics_trans
                                                  !! in a time-stepping cycle; missing is like true.
  real,            optional, intent(in)    :: cycle_length !< The duration of a coupled time stepping cycle [s].

  ! These pointers are used to simplify the code below.
  type(ice_grid_type),     pointer :: sIG => NULL()
  type(SIS_hor_grid_type), pointer :: sG => NULL()
  type(ice_state_type),    pointer :: sIST => NULL()
  type(fast_ice_avg_type), pointer :: FIA => NULL()
  type(unit_scale_type),   pointer :: US => NULL()
  real :: dt_slow  ! The time step over which to advance the model [T ~> s].
  logical :: do_multi_trans, cycle_start

  if (.not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "The pointer to Ice%sCS must be associated in update_ice_dynamics_trans.")

  sIST => Ice%sCS%IST ; sIG => Ice%sCS%IG ; sG => Ice%sCS%G ; FIA => Ice%sCS%FIA ; US => Ice%sCS%US
  dt_slow = US%s_to_T*time_type_to_real(Ice%sCS%Time_step_slow)
  if (present(time_step)) dt_slow = US%s_to_T*time_type_to_real(time_step)
  cycle_start = .true. ; if (present(start_cycle)) cycle_start = start_cycle

  call cpu_clock_begin(iceClock) ; call cpu_clock_begin(ice_clock_slow)

  ! Do halo updates on the forcing fields, as necessary.  This must occur before
  ! the call to SIS_dynamics_trans, because update_icebergs does its own halo
  ! updates, and slow_thermodynamics only works on the computational domain.
  if (cycle_start) then
    call pass_vector(FIA%WindStr_x, FIA%WindStr_y, sG%Domain, stagger=AGRID, complete=.false.)
    call pass_vector(FIA%WindStr_ocn_x, FIA%WindStr_ocn_y, sG%Domain, stagger=AGRID)
    call pass_var(FIA%ice_cover, sG%Domain, complete=.false.)
    call pass_var(FIA%ice_free,  sG%Domain, complete=.true.)
  endif
  if (sIST%valid_IST) then
    call pass_var(sIST%part_size, sG%Domain)
    call pass_var(sIST%mH_ice, sG%Domain, complete=.false.)
    call pass_var(sIST%mH_pond, sG%Domain, complete=.false.)
    call pass_var(sIST%mH_snow, sG%Domain, complete=.true.)
  endif

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Before SIS_dynamics_trans", Ice, check_slow=.true.)
  endif

  do_multi_trans = (present(start_cycle) .or. present(end_cycle) .or. present(cycle_length))

  if (Ice%sCS%specified_ice) then ! There is no ice dynamics or transport.
    call specified_ice_dynamics(sIST, Ice%sCS%OSS, FIA, Ice%sCS%IOF, dt_slow, &
                                Ice%sCS%specified_ice_CSp, sG, US, sIG)
  elseif (do_multi_trans) then
    call SIS_multi_dyn_trans(sIST, Ice%sCS%OSS, FIA, Ice%sCS%IOF, dt_slow, Ice%sCS%dyn_trans_CSp, &
                             Ice%icebergs, sG, US, sIG, Ice%sCS%SIS_tracer_flow_CSp, &
                             Ice%OBC, start_cycle, end_cycle, cycle_length)
  elseif (Ice%sCS%slab_ice) then ! Use a very old slab ice model.
    call slab_ice_dyn_trans(sIST, Ice%sCS%OSS, FIA, Ice%sCS%IOF, dt_slow, Ice%sCS%dyn_trans_CSp, &
                            sG, US, sIG, Ice%sCS%SIS_tracer_flow_CSp, Ice%OBC)
  else ! This is the typical branch used by SIS2.
    call SIS_dynamics_trans(sIST, Ice%sCS%OSS, FIA, Ice%sCS%IOF, dt_slow, Ice%sCS%dyn_trans_CSp, &
                            Ice%icebergs, sG, US, sIG, Ice%sCS%SIS_tracer_flow_CSp, Ice%OBC)
  endif

 ! Set up the stresses and surface pressure in the externally visible structure Ice.
  if (sIST%valid_IST) call ice_mass_from_IST(sIST, Ice%sCS%IOF, sG, sIG)

  if (Ice%sCS%debug) then
    call IOF_chksum("Before set_ocean_top_dyn_fluxes", Ice%sCS%IOF, sG, US, mech_fluxes=.true.)
  endif
  call set_ocean_top_dyn_fluxes(Ice, Ice%sCS%IOF, FIA, sG, US, Ice%sCS)

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("End update_ice_dynamics_trans", Ice, check_slow=.true.)
  endif

  !### THIS NO LONGER WORKS ON SLOW ICE PES.
!  if (Ice%sCS%bounds_check) then
!    call Ice_public_type_bounds_check(Ice, sG, "End update_ice_slow")
!  endif

  call cpu_clock_end(ice_clock_slow) ; call cpu_clock_end(iceClock)

end subroutine update_ice_dynamics_trans

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_model_fast_cleanup performs the final steps in the fast ice update cycle
!! and prepares data to drive the slow ice updates.  This includes finding the
!! averaged fluxes and unpacking the land to ice forcing.
subroutine ice_model_fast_cleanup(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type.

  if (.not.associated(Ice%fCS)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS must be associated in ice_model_fast_cleanup.")

  ! average fluxes from update_ice_model_fast
  call avg_top_quantities(Ice%fCS%FIA, Ice%fCS%Rad, Ice%fCS%IST, Ice%fCS%G, Ice%fCS%US, Ice%fCS%IG)

  call total_top_quantities(Ice%fCS%FIA, Ice%fCS%TSF, Ice%fCS%IST%part_size, Ice%fCS%G, Ice%fCS%US, Ice%fCS%IG)

  if (allocated(Ice%fCS%IST%t_surf)) &
    Ice%fCS%IST%t_surf(:,:,1:) = Ice%fCS%Rad%T_skin(:,:,:) + Ice%fCS%IST%T_0degC
  call infill_array(Ice%fCS%IST, Ice%fCS%sOSS%T_fr_ocn, Ice%fCS%Rad%T_skin, Ice%fCS%G, Ice%fCS%IG)

end subroutine ice_model_fast_cleanup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> unpack_land_ice_bdry converts the information in a publicly visible
!! land_ice_boundary_type into an internally visible fast_ice_avg_type variable.
subroutine unpack_land_ice_boundary(Ice, LIB)
  type(ice_data_type),          intent(inout) :: Ice !< The publicly visible ice data type.
  type(land_ice_boundary_type), intent(in)    :: LIB !< The land ice boundary type that is being unpacked.

  type(fast_ice_avg_type), pointer :: FIA => NULL()
  type(SIS_hor_grid_type), pointer :: G => NULL()
  type(unit_scale_type),   pointer :: US => NULL()

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off

  if (.not.associated(Ice%fCS)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS must be associated in unpack_land_ice_boundary.")
  if (.not.associated(Ice%fCS%FIA)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS%FIA must be associated in unpack_land_ice_boundary.")
  if (.not.associated(Ice%fCS%G)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS%G must be associated in unpack_land_ice_boundary.")

  FIA => Ice%fCS%FIA ; G => Ice%fCS%G
  US => Ice%fCS%US

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  ! Store liquid runoff and other fluxes from the land to the ice or ocean.
  i_off = LBOUND(LIB%runoff,1) - G%isc ; j_off = LBOUND(LIB%runoff,2) - G%jsc
  !$OMP parallel do default(none) shared(isc,iec,jsc,jec,FIA,LIB,i_off,j_off,G,US) &
  !$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j) > 0.0) then
    i2 = i+i_off ; j2 = j+j_off
    FIA%runoff(i,j)  = US%kg_m2s_to_RZ_T*LIB%runoff(i2,j2)
    FIA%calving(i,j) = US%kg_m2s_to_RZ_T*LIB%calving(i2,j2)
    FIA%runoff_hflx(i,j)  = US%W_m2_to_QRZ_T*LIB%runoff_hflx(i2,j2)
    FIA%calving_hflx(i,j) = US%W_m2_to_QRZ_T*LIB%calving_hflx(i2,j2)
  else
    ! This is a land point from the perspective of the sea-ice.
    ! At some point it might make sense to check for non-zero fluxes, which
    ! might indicate regridding errors.  However, bad-data values are also
    ! non-zero and should not be flagged.
    FIA%runoff(i,j)  = 0.0 ; FIA%calving(i,j) = 0.0
    FIA%runoff_hflx(i,j)  = 0.0 ; FIA%calving_hflx(i,j) = 0.0
  endif ; enddo ; enddo

  if (Ice%fCS%debug) then
    call FIA_chksum("End of unpack_land_ice_boundary", FIA, G, Ice%fCS%US)
  endif

end subroutine unpack_land_ice_boundary

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> This subroutine copies information (mostly fluxes and the updated temperatures)
!! from the fast part of the sea-ice to the  slow part of the sea ice.
subroutine exchange_fast_to_slow_ice(Ice)
  type(ice_data_type), &
    intent(inout) :: Ice            !< The publicly visible ice data type whose fast
                                    !! part is to be exchanged with the slow part.
  type(fast_ice_avg_type),   pointer :: FIA_null => NULL()
  type(ice_state_type),      pointer :: IST_null => NULL()
  type(ice_rad_type),        pointer :: Rad_null => NULL()
  type(total_sfc_flux_type), pointer :: TSF_null => NULL()

  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  logical :: redo_fast_update

  redo_fast_update = .false.
  if (associated(Ice%fCS)) redo_fast_update = Ice%fCS%redo_fast_update
  if (associated(Ice%sCS)) redo_fast_update = Ice%sCS%redo_fast_update

  if (associated(Ice%fCS)) then
    isc = Ice%fCS%G%isc ; iec = Ice%fCS%G%iec ; jsc = Ice%fCS%G%jsc ; jec = Ice%fCS%G%jec
    isd = Ice%fCS%G%isd ; ied = Ice%fCS%G%ied ; jsd = Ice%fCS%G%jsd ; jed = Ice%fCS%G%jed

    ! Propagate the coupler_type info to Ice%fCS%FIA%tr_flux and allocate its arrays.
    call coupler_type_spawn(Ice%ocean_fluxes, Ice%fCS%FIA%tr_flux, &
                            (/isd, isc, iec, ied/),  (/jsd, jsc, jec, jed/), &
                            (/0, Ice%fCS%IG%CatIce/), as_needed=.true.)

    if (redo_fast_update) &
      ! Propagate the coupler_type info to Ice%fCS%TSF%tr_flux and allocate its arrays.
      call coupler_type_spawn(Ice%ocean_fluxes, Ice%fCS%TSF%tr_flux, &
                              (/isd, isc, iec, ied/),  (/jsd, jsc, jec, jed/), as_needed=.true.)
  endif

  if (associated(Ice%sCS)) then
    isc = Ice%sCS%G%isc ; iec = Ice%sCS%G%iec ; jsc = Ice%sCS%G%jsc ; jec = Ice%sCS%G%jec
    isd = Ice%sCS%G%isd ; ied = Ice%sCS%G%ied ; jsd = Ice%sCS%G%jsd ; jed = Ice%sCS%G%jed

    ! Propagate the coupler_type info to Ice%sCS%FIA%tr_flux and allocate its arrays.
    call coupler_type_spawn(Ice%ocean_fluxes, Ice%sCS%FIA%tr_flux, &
                            (/isd, isc, iec, ied/),  (/jsd, jsc, jec, jed/), &
                            (/0, Ice%sCS%IG%CatIce/), as_needed=.true.)

    if (redo_fast_update) &
      ! Propagate the coupler_type info to Ice%sCS%TSF%tr_flux and allocate its arrays.
      call coupler_type_spawn(Ice%ocean_fluxes, Ice%sCS%TSF%tr_flux, &
                              (/isd, isc, iec, ied/),  (/jsd, jsc, jec, jed/), as_needed=.true.)
  endif

  if (Ice%xtype == DIRECT) then
    if (.not.associated(Ice%fCS) .or. .not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "With xtype=DIRECT, both the pointer to Ice%sCS and the pointer to Ice%fCS must be "//&
      "associated (although perhaps not with each other) in exchange_fast_to_slow_ice.")

    if (.not.associated(Ice%fCS%FIA, Ice%sCS%FIA)) then
      call copy_FIA_to_FIA(Ice%fCS%FIA, Ice%sCS%FIA, Ice%fCS%G%HI, Ice%sCS%G%HI, Ice%sCS%IG)
    endif

    if (redo_fast_update) then
      if (.not.associated(Ice%fCS%TSF, Ice%sCS%TSF)) &
        call copy_TSF_to_TSF(Ice%fCS%TSF, Ice%sCS%TSF, Ice%fCS%G%HI, Ice%sCS%G%HI)
      if (.not.associated(Ice%fCS%Rad, Ice%sCS%Rad)) &
        call copy_Rad_to_Rad(Ice%fCS%Rad, Ice%sCS%Rad, Ice%fCS%G%HI, Ice%sCS%G%HI, Ice%fCS%IG)
    else
      if (.not.associated(Ice%fCS%IST, Ice%sCS%IST)) &
        call copy_IST_to_IST(Ice%fCS%IST, Ice%sCS%IST, Ice%fCS%G%HI, Ice%sCS%G%HI, Ice%fCS%IG)
    endif
  elseif (Ice%xtype == REDIST) then
    if (.not.associated(Ice%fCS) .and. .not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "Either the pointer to Ice%sCS or the pointer to Ice%fCS must be "//&
      "associated in exchange_fast_to_slow_ice.")

    if (associated(Ice%fCS) .and. associated(Ice%sCS)) then
      if (.not.associated(Ice%fCS%FIA, Ice%sCS%FIA)) &
        call redistribute_FIA_to_FIA(Ice%fCS%FIA, Ice%sCS%FIA, Ice%fast_domain, &
                                     Ice%slow_domain, Ice%sCS%G, Ice%sCS%IG)

      if (redo_fast_update) then
        call redistribute_TSF_to_TSF(Ice%fCS%TSF, Ice%sCS%TSF, Ice%fast_domain, &
                                     Ice%slow_domain, Ice%sCS%G%HI)
        call redistribute_Rad_to_Rad(Ice%fCS%Rad, Ice%sCS%Rad, Ice%fast_domain, Ice%slow_domain)
      else
        if (.not.associated(Ice%fCS%IST, Ice%sCS%IST)) &
          call redistribute_IST_to_IST(Ice%fCS%IST, Ice%sCS%IST, Ice%fast_domain, Ice%slow_domain)
      endif
    elseif (associated(Ice%fCS)) then
      call redistribute_FIA_to_FIA(Ice%fCS%FIA, FIA_null, Ice%fast_domain, Ice%slow_domain)
      if (redo_fast_update) then
        call redistribute_TSF_to_TSF(Ice%fCS%TSF, TSF_null, Ice%fast_domain, Ice%slow_domain)
        call redistribute_Rad_to_Rad(Ice%fCS%Rad, Rad_null, Ice%fast_domain, Ice%slow_domain)
      else
        call redistribute_IST_to_IST(Ice%fCS%IST, IST_null, Ice%fast_domain, Ice%slow_domain)
      endif
    elseif (associated(Ice%sCS)) then
      call redistribute_FIA_to_FIA(FIA_null, Ice%sCS%FIA, Ice%fast_domain, &
                                   Ice%slow_domain, Ice%sCS%G, Ice%sCS%IG)
      if (redo_fast_update) then
        call redistribute_TSF_to_TSF(TSF_null, Ice%sCS%TSF, Ice%fast_domain, &
                                     Ice%slow_domain, Ice%sCS%G%HI)
        call redistribute_Rad_to_Rad(Rad_null, Ice%sCS%Rad, Ice%fast_domain, Ice%slow_domain)
      else
        call redistribute_IST_to_IST(IST_null, Ice%sCS%IST, Ice%fast_domain, Ice%slow_domain)
      endif
    else
      call SIS_error(FATAL, "Either the pointer to Ice%sCS or the pointer to "//&
                     "Ice%fCS must be associated in exchange_fast_to_slow_ice.")
    endif
  else
    call SIS_error(FATAL, "exchange_fast_to_slow_ice called with an unrecognized Ice%xtype value.")
  endif

end subroutine exchange_fast_to_slow_ice


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ocean_top_fluxes translates ice-bottom fluxes of heat, mass, salt, and
!!  tracers from the ice model's internal state to the public ice data type
!!  for use by the ocean model.
subroutine set_ocean_top_fluxes(Ice, IST, IOF, FIA, OSS, G, US, IG, sCS)
  type(ice_data_type),        intent(inout) :: Ice !< The publicly visible ice data type.
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ice_ocean_flux_type),  intent(in)    :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  type(fast_ice_avg_type),    intent(in)    :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(in)    :: IG  !< The sea-ice specific grid type
  type(SIS_slow_CS),          intent(in)    :: sCS !< The slow ice control structure

  real :: I_count
  integer :: i, j, k, isc, iec, jsc, jec, m, n
  integer :: i2, j2, i_off, j_off, ncat, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  if (sCS%debug) then
    call Ice_public_type_chksum("Start set_ocean_top_fluxes", Ice, check_slow=.true.)
    call IOF_chksum("Start set_ocean_top_fluxes", IOF, G, sCS%US, thermo_fluxes=.true.)
    call FIA_chksum("Start set_ocean_top_fluxes", FIA, G, US)
  endif

!   It is possible that the ice mass and surface pressure will be needed after
! the thermodynamic step, in which case this should be uncommented.
!  ! Sum the concentration weighted mass.
!  Ice%mi(:,:) = 0.0
!  i_off = LBOUND(Ice%mi,1) - G%isc ; j_off = LBOUND(Ice%mi,2) - G%jsc
!  !$OMP parallel do default(shared) private(i2,j2)
!  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
!    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
!    Ice%mi(i2,j2) = Ice%mi(i2,j2) + IST%part_size(i,j,k) * &
!        (G%US%RZ_to_kg_m2 * ((IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k)) + IST%mH_ice(i,j,k)))
!  enddo ; enddo ; enddo

  ! This block of code is probably unnecessary.
  Ice%flux_t(:,:) = 0.0 ; Ice%flux_q(:,:) = 0.0
  Ice%flux_sw_nir_dir(:,:) = 0.0 ; Ice%flux_sw_nir_dif(:,:) = 0.0
  Ice%flux_sw_vis_dir(:,:) = 0.0 ; Ice%flux_sw_vis_dif(:,:) = 0.0
  Ice%flux_lw(:,:) = 0.0 ; Ice%flux_lh(:,:) = 0.0
  Ice%fprec(:,:) = 0.0 ; Ice%lprec(:,:) = 0.0
  call coupler_type_rescale_data(Ice%ocean_fluxes, 0.0)

  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc
  !$OMP parallel do default(none) shared(isc,iec,jsc,jec,Ice,IST,IOF,FIA,i_off,j_off,G,US,OSS) &
  !$OMP                           private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
    Ice%flux_t(i2,j2) = US%QRZ_T_to_W_m2*IOF%flux_sh_ocn_top(i,j)
    Ice%flux_q(i2,j2) = US%RZ_T_to_kg_m2s*IOF%evap_ocn_top(i,j)
    Ice%flux_sw_vis_dir(i2,j2) = US%QRZ_T_to_W_m2*IOF%flux_sw_ocn(i,j,VIS_DIR)
    Ice%flux_sw_vis_dif(i2,j2) = US%QRZ_T_to_W_m2*IOF%flux_sw_ocn(i,j,VIS_DIF)
    Ice%flux_sw_nir_dir(i2,j2) = US%QRZ_T_to_W_m2*IOF%flux_sw_ocn(i,j,NIR_DIR)
    Ice%flux_sw_nir_dif(i2,j2) = US%QRZ_T_to_W_m2*IOF%flux_sw_ocn(i,j,NIR_DIF)
    Ice%flux_lw(i2,j2) = US%QRZ_T_to_W_m2*IOF%flux_lw_ocn_top(i,j)
    Ice%flux_lh(i2,j2) = US%QRZ_T_to_W_m2*IOF%flux_lh_ocn_top(i,j)
    Ice%fprec(i2,j2) = US%RZ_T_to_kg_m2s*IOF%fprec_ocn_top(i,j)
    Ice%lprec(i2,j2) = US%RZ_T_to_kg_m2s*IOF%lprec_ocn_top(i,j)
    Ice%runoff(i2,j2)  = US%RZ_T_to_kg_m2s*FIA%runoff(i,j)
    Ice%calving(i2,j2) = US%RZ_T_to_kg_m2s*FIA%calving(i,j)
    Ice%runoff_hflx(i2,j2)  = US%QRZ_T_to_W_m2*FIA%runoff_hflx(i,j)
    Ice%calving_hflx(i2,j2) = US%QRZ_T_to_W_m2*FIA%calving_hflx(i,j)
    Ice%flux_salt(i2,j2) = US%S_to_ppt*US%RZ_T_to_kg_m2s*IOF%flux_salt(i,j)
    Ice%SST_C(i2,j2) = US%C_to_degC*OSS%SST_C(i,j)

!   It is possible that the ice mass and surface pressure will be needed after
! the thermodynamic step, in which case this should be uncommented.
!  if (IOF%slp2ocean) then
!     Ice%p_surf(i2,j2) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*FIA%p_atm_surf(i,j) - 1e5 ! SLP - 1 std. atmosphere [Pa].
!   else
!     Ice%p_surf(i2,j2) = 0.0
!   endif
!   Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) + US%L_T_to_m_s**2*US%m_to_Z*G%g_Earth*Ice%mi(i2,j2)
  enddo ; enddo
  if (allocated(IOF%melt_nudge)) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
      Ice%lprec(i2,j2) = Ice%lprec(i2,j2) + US%RZ_T_to_kg_m2s*IOF%melt_nudge(i,j)
    enddo ; enddo
  endif

  ! This copy may need to be skipped in the first step of a cold-start run with lagged ice
  ! coupling, but otherwise if it is skipped may indicate a problem that should be trapped.
  if (coupler_type_initialized(IOF%tr_flux_ocn_top)) &
    call coupler_type_copy_data(IOF%tr_flux_ocn_top, Ice%ocean_fluxes)

  if (sCS%debug) then
    call Ice_public_type_chksum("End set_ocean_top_fluxes", Ice, check_slow=.true.)
  endif

end subroutine set_ocean_top_fluxes

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_mass_from_IST stores the total ice mass determined from IST in the IOF type.
subroutine ice_mass_from_IST(IST, IOF, G, IG)
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),        intent(in)    :: IG  !< The sea-ice specific grid type

  integer :: i, j, k, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  ! Sum the concentration weighted mass.
  IOF%mass_ice_sn_p(:,:) = 0.0
  !$OMP parallel do default(shared)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    IOF%mass_ice_sn_p(i,j) = IOF%mass_ice_sn_p(i,j) + IST%part_size(i,j,k) * &
          ((IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k)) + IST%mH_ice(i,j,k))
  enddo ; enddo ; enddo

end subroutine ice_mass_from_IST


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ocean_top_dyn_fluxes translates ice-bottom stresses and mass from the ice
!!  model's ice-ocean flux type and the fast-ice average type to the public
!!  ice data type for use by the ocean model.
subroutine set_ocean_top_dyn_fluxes(Ice, IOF, FIA, G, US, sCS)
  type(ice_data_type),        intent(inout) :: Ice !< The publicly visible ice data type.
  type(ice_ocean_flux_type),  intent(in)    :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  type(fast_ice_avg_type),    intent(in)    :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_slow_CS),          intent(in)    :: sCS !< The slow ice control structure

  real :: I_count
  integer :: i, j, k, isc, iec, jsc, jec
  integer :: i2, j2, i_off, j_off, ind
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (sCS%debug) then
    call Ice_public_type_chksum("Start set_ocean_top_dyn_fluxes", Ice, check_slow=.true.)
    call IOF_chksum("Start set_ocean_top_dyn_fluxes", IOF, G, US, mech_fluxes=.true.)
  endif

  ! Sum the concentration weighted mass.
  Ice%mi(:,:) = 0.0
  i_off = LBOUND(Ice%mi,1) - G%isc ; j_off = LBOUND(Ice%mi,2) - G%jsc
  !$OMP parallel do default(shared) private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
    Ice%mi(i2,j2) = Ice%mi(i2,j2) + US%RZ_to_kg_m2*IOF%mass_ice_sn_p(i,j)
  enddo ; enddo

  if (sCS%do_icebergs .and. associated(IOF%mass_berg)) then
    ! Note that the IOF berg fields and Ice fields are only allocated on the
    ! computational domains, although they may use different indexing conventions.
    Ice%mi(:,:) = Ice%mi(:,:) + IOF%mass_berg(:,:)
    if (sCS%pass_iceberg_area_to_ocean) then
      Ice%mass_berg(:,:) = IOF%mass_berg(:,:)
      if (associated(IOF%ustar_berg)) Ice%ustar_berg(:,:) = IOF%ustar_berg(:,:)
      if (associated(IOF%area_berg))  Ice%area_berg(:,:) = IOF%area_berg(:,:)
    endif
  endif

  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc
  !$OMP parallel do default(shared) private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
    Ice%flux_u(i2,j2) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*IOF%flux_u_ocn(i,j)
    Ice%flux_v(i2,j2) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*IOF%flux_v_ocn(i,j)

    if (IOF%slp2ocean) then
      Ice%p_surf(i2,j2) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*FIA%p_atm_surf(i,j) - 1e5 ! SLP - 1 std. atmosphere [Pa].
    else
      Ice%p_surf(i2,j2) = 0.0
    endif
    Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) + US%L_T_to_m_s**2*US%m_to_Z*G%g_Earth*Ice%mi(i2,j2)
  enddo ; enddo
  if (associated(Ice%stress_mag) .and. allocated(IOF%stress_mag)) then
    i_off = LBOUND(Ice%stress_mag,1) - G%isc ; j_off = LBOUND(Ice%stress_mag,2) - G%jsc
    !$OMP parallel do default(shared) private(i2,j2)
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      Ice%stress_mag(i2,j2) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*IOF%stress_mag(i,j)
    enddo ; enddo
  endif

  if (sCS%debug) then
    call Ice_public_type_chksum("End set_ocean_top_dyn_fluxes", Ice, check_slow=.true.)
  endif

end subroutine set_ocean_top_dyn_fluxes


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> This subroutine copies information from the slow part of the sea-ice to the
!! fast part of the sea ice.
subroutine exchange_slow_to_fast_ice(Ice)
  type(ice_data_type), &
    intent(inout) :: Ice            !< The publicly visible ice data type whose slow
                                    !! part is to be exchanged with the fast part.
  type(simple_OSS_type), pointer :: sOSS_null => NULL()
  type(ice_state_type),  pointer :: IST_null => NULL()

  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed

  call cpu_clock_begin(iceClock) ; call cpu_clock_begin(ice_clock_exchange)

  if (associated(Ice%fCS)) then
    isc = Ice%fCS%G%isc ; iec = Ice%fCS%G%iec ; jsc = Ice%fCS%G%jsc ; jec = Ice%fCS%G%jec
    isd = Ice%fCS%G%isd ; ied = Ice%fCS%G%ied ; jsd = Ice%fCS%G%jsd ; jed = Ice%fCS%G%jed

    ! Propagate the coupler_type info to Ice%fCS%sOSS%tr_fields and allocate its arrays.
    call coupler_type_spawn(Ice%ocean_fields, Ice%fCS%sOSS%tr_fields, &
                            (/isd, isc, iec, ied/),  (/jsd, jsc, jec, jed/), as_needed=.true. )
  endif

  if (Ice%xtype == DIRECT) then
    if (.not.associated(Ice%fCS) .or. .not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "With xtype=DIRECT, both the pointer to Ice%sCS and the pointer to Ice%fCS must be "//&
      "associated (although perhaps not with each other) in exchange_slow_to_fast_ice.")

    if (.not.associated(Ice%fCS%sOSS, Ice%sCS%sOSS)) then
      call copy_sOSS_to_sOSS(Ice%sCS%sOSS, Ice%fCS%sOSS, Ice%sCS%G%HI, Ice%fCS%G%HI)
    endif

    if (.not.associated(Ice%fCS%IST, Ice%sCS%IST)) then
      call copy_IST_to_IST(Ice%sCS%IST, Ice%fCS%IST, Ice%sCS%G%HI, Ice%fCS%G%HI, &
           Ice%sCS%IG)
    endif

  elseif (Ice%xtype == REDIST) then
    if (.not.associated(Ice%fCS) .and. .not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "Either the pointer to Ice%sCS or the pointer to Ice%fCS must be "//&
      "associated in exchange_slow_to_fast_ice.")

    if (associated(Ice%fCS) .and. associated(Ice%sCS)) then
      if (.not.associated(Ice%fCS%sOSS, Ice%sCS%sOSS)) &
        call redistribute_sOSS_to_sOSS(Ice%sCS%sOSS, Ice%fCS%sOSS, Ice%slow_domain, &
                                       Ice%fast_domain, Ice%fCS%G%HI)

      if (.not.associated(Ice%fCS%IST, Ice%sCS%IST)) &
        call redistribute_IST_to_IST(Ice%sCS%IST, Ice%fCS%IST, Ice%slow_domain, &
                                     Ice%fast_domain)
    elseif (associated(Ice%fCS)) then
      call redistribute_sOSS_to_sOSS(sOSS_null, Ice%fCS%sOSS, Ice%slow_domain, &
                                     Ice%fast_domain, HI_out=Ice%fCS%G%HI)
      call redistribute_IST_to_IST(IST_null, Ice%fCS%IST, Ice%slow_domain, &
                                   Ice%fast_domain)
    elseif (associated(Ice%sCS)) then
      call redistribute_sOSS_to_sOSS(Ice%sCS%sOSS, sOSS_null, Ice%slow_domain, &
                                     Ice%fast_domain)
      call redistribute_IST_to_IST(Ice%sCS%IST, IST_null, Ice%slow_domain, &
                                   Ice%fast_domain)
    else
      call SIS_error(FATAL, "Either the pointer to Ice%sCS or the pointer to "//&
                     "Ice%fCS must be associated in exchange_slow_to_fast_ice.")
    endif
  else
    call SIS_error(FATAL, "exchange_slow_to_fast_ice called with an unrecognized Ice%xtype value.")
  endif

  call cpu_clock_end(ice_clock_exchange) ; call cpu_clock_end(iceClock)

end subroutine exchange_slow_to_fast_ice

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> This subroutine copies information from an ocean_ice_boundary_type into the
!! slow part of an ice_data type, using a coupler-friendly interface.
subroutine unpack_ocean_ice_boundary(Ocean_boundary, Ice)
  type(ocean_ice_boundary_type), &
    intent(inout) :: Ocean_boundary !< A structure containing information about
                                    !! the ocean that is being shared with the sea-ice.
  type(ice_data_type), &
    intent(inout) :: Ice            !< The publicly visible ice data type in the slow part
                                    !! of which the ocean surface information is to be stored.

  if (.not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "The pointer to Ice%sCS must be associated in unpack_ocean_ice_boundary.")

  call unpack_ocn_ice_bdry(Ocean_boundary, Ice%sCS%OSS, Ice%sCS%IST%ITV, Ice%sCS%G, Ice%sCS%US, &
                           Ice%sCS%specified_ice, Ice%ocean_fields)

  call translate_OSS_to_sOSS(Ice%sCS%OSS, Ice%sCS%IST, Ice%sCS%sOSS, Ice%sCS%G, Ice%sCS%US)

end subroutine unpack_ocean_ice_boundary

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> This subroutine converts the information in a publicly visible
!! ocean_ice_boundary_type into an internally visible ocean_sfc_state_type
!! variable.
subroutine unpack_ocn_ice_bdry(OIB, OSS, ITV, G, US, specified_ice, ocean_fields)
  type(ocean_ice_boundary_type), intent(in)    :: OIB !< A type containing ocean surface fields that
                                                      !! are used to drive the sea ice
  type(ocean_sfc_state_type),    intent(inout) :: OSS !< A structure containing the arrays that describe
                                                      !! the ocean's surface state for the ice model.
  type(ice_thermo_type),         intent(in)    :: ITV !< The ice thermodynamics parameter structure.
  type(SIS_hor_grid_type),       intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),         intent(in)    :: US  !< A structure with unit conversion factors
  logical,                       intent(in)    :: specified_ice !< If true, use specified ice properties.
  type(coupler_3d_bc_type),      intent(inout) :: ocean_fields  !< A structure of ocean fields, often
                                                                !! related to passive tracers.

  real, dimension(G%isd:G%ied, G%jsd:G%jed) :: u_nonsym, v_nonsym ! Nonsymmetric velocities [L T-1 ~> m s-1]
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  logical :: Cgrid_ocn
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, index
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = LBOUND(OIB%t,1) - G%isc ; j_off = LBOUND(OIB%t,2) - G%jsc

  call cpu_clock_begin(iceClock) ; call cpu_clock_begin(ice_clock_slow)

  ! Pass the ocean state through ice on partition 0, unless using specified ice.
  if (.not. specified_ice) then
    !$OMP parallel do default(shared) private(i2,j2)
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      OSS%SST_C(i,j) = US%degC_to_C * (OIB%t(i2,j2) - T_0degC)
    enddo ; enddo
  endif

  !$OMP parallel do default(shared) private(i2,j2)
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    OSS%s_surf(i,j) = US%ppt_to_S*OIB%s(i2,j2)
    OSS%T_fr_ocn(i,j) = T_Freeze(OSS%s_surf(i,j), ITV)
    OSS%bheat(i,j) = OSS%kmelt*(OSS%SST_C(i,j) - OSS%T_fr_ocn(i,j))
    OSS%frazil(i,j) = US%W_m2_to_QRZ_T*US%s_to_T*OIB%frazil(i2,j2)
    OSS%sea_lev(i,j) = US%m_to_Z*OIB%sea_level(i2,j2)
  enddo ; enddo

  Cgrid_ocn = (allocated(OSS%u_ocn_C) .and. allocated(OSS%v_ocn_C))

  ! Unpack the ocean surface velocities.
  if (OIB%stagger == AGRID) then
    u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      u_nonsym(i,j) = US%m_s_to_L_T*OIB%u(i2,j2) ; v_nonsym(i,j) = US%m_s_to_L_T*OIB%v(i2,j2)
    enddo ; enddo
    call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, stagger=AGRID)

    if (Cgrid_ocn) then
      do j=jsc,jec ; do I=isc-1,iec
        OSS%u_ocn_C(I,j) = 0.5*(u_nonsym(i,j) + u_nonsym(i+1,j))
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        OSS%v_ocn_C(i,J) = 0.5*(v_nonsym(i,j) + v_nonsym(i,j+1))
      enddo ; enddo
    else
      do J=jsc-1,jec ; do I=isc-1,iec
        OSS%u_ocn_B(I,J) = 0.25*((u_nonsym(i,j) + u_nonsym(i+1,j+1)) + &
                               (u_nonsym(i+1,j) + u_nonsym(i,j+1)))
        OSS%v_ocn_B(I,J) = 0.25*((v_nonsym(i,j) + v_nonsym(i+1,j+1)) + &
                               (v_nonsym(i+1,j) + v_nonsym(i,j+1)))
      enddo ; enddo
    endif

  elseif (OIB%stagger == BGRID_NE) then
    if (Cgrid_ocn) then
      u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
      do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        u_nonsym(i,j) = US%m_s_to_L_T*OIB%u(i2,j2) ; v_nonsym(i,j) = US%m_s_to_L_T*OIB%v(i2,j2)
      enddo ; enddo
      call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, stagger=BGRID_NE)

      do j=jsc,jec ; do I=isc-1,iec
        OSS%u_ocn_C(I,j) = 0.5*(u_nonsym(I,J) + u_nonsym(I,J-1))
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        OSS%v_ocn_C(i,J) = 0.5*(v_nonsym(I,J) + v_nonsym(I-1,J))
      enddo ; enddo
    else
      do J=jsc,jec ; do I=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        OSS%u_ocn_B(I,J) = US%m_s_to_L_T*OIB%u(i2,j2)
        OSS%v_ocn_B(I,J) = US%m_s_to_L_T*OIB%v(i2,j2)
      enddo ; enddo
      if (G%symmetric) &
        call fill_symmetric_edges(OSS%u_ocn_B, OSS%v_ocn_B, G%Domain, stagger=BGRID_NE)
    endif

  elseif (OIB%stagger == CGRID_NE) then
    if (Cgrid_ocn) then
      do j=jsc,jec ; do I=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        OSS%u_ocn_C(I,j) = US%m_s_to_L_T*OIB%u(i2,j2)
      enddo ; enddo
      do J=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        OSS%v_ocn_C(i,J) = US%m_s_to_L_T*OIB%v(i2,j2)
      enddo ; enddo
      if (G%symmetric) &
        call fill_symmetric_edges(OSS%u_ocn_C, OSS%v_ocn_C, G%Domain, stagger=CGRID_NE)
    else
      u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
      do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        u_nonsym(I,j) = US%m_s_to_L_T*OIB%u(i2,j2) ; v_nonsym(i,J) = US%m_s_to_L_T*OIB%v(i2,j2)
      enddo ; enddo
      call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, stagger=CGRID_NE)
      do J=jsc-1,jec ; do I=isc-1,iec
        OSS%u_ocn_B(I,J) = 0.5*(u_nonsym(I,j) + u_nonsym(I,j+1))
        OSS%v_ocn_B(I,J) = 0.5*(v_nonsym(i,J) + v_nonsym(i+1,J))
      enddo ; enddo
    endif
  else
    call SIS_error(FATAL, "unpack_ocn_ice_bdry: Unrecognized OIB%stagger.")
  endif

  ! Fill in the halo values.
  if (Cgrid_ocn) then
    call pass_vector(OSS%u_ocn_C, OSS%v_ocn_C, G%Domain, stagger=CGRID_NE)
  else
    call pass_vector(OSS%u_ocn_B, OSS%v_ocn_B, G%Domain, stagger=BGRID_NE)
  endif
  call pass_var(OSS%sea_lev, G%Domain)

! Transfer the ocean state for extra tracer fluxes.
  call coupler_type_spawn(OIB%fields, OSS%tr_fields, (/isd, isc, iec, ied/), &
                          (/jsd, jsc, jec, jed/), as_needed=.true. )
  call coupler_type_copy_data(OIB%fields, OSS%tr_fields)

  call cpu_clock_end(ice_clock_slow) ; call cpu_clock_end(iceClock)

end subroutine unpack_ocn_ice_bdry

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ice_surface_fields prepares the ice surface state for atmosphere fast
!! physics and does precalculation of ice radiative properties.
subroutine set_ice_surface_fields(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type whose
                                            !! surface properties are being set.

  call cpu_clock_begin(iceClock) ; call cpu_clock_begin(ice_clock_fast)
  if (.not.associated(Ice%fCS)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS must be associated in set_ice_surface_fields.")

  call set_ice_surface_state(Ice, Ice%fCS%IST, Ice%fCS%sOSS, Ice%fCS%Rad, &
                             Ice%fCS%FIA, Ice%fCS%G, Ice%fCS%US, Ice%fCS%IG, Ice%fCS )

  call cpu_clock_end(ice_clock_fast) ; call cpu_clock_end(iceClock)
end subroutine set_ice_surface_fields

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ice_surface_state prepares the surface state for atmosphere fast physics
subroutine set_ice_surface_state(Ice, IST, OSS, Rad, FIA, G, US, IG, fCS)
  type(ice_data_type),        intent(inout) :: Ice !< The publicly visible ice data type.
  type(ice_state_type),       intent(in)    :: IST !< A type describing the state of the sea ice
  type(simple_OSS_type),      intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(ice_rad_type),         intent(inout) :: Rad !< A structure with fields related to the absorption,
                                                   !! reflection and transmission of shortwave radiation.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(in)    :: IG  !< The sea-ice specific grid type
  type(SIS_fast_CS),          intent(inout) :: fCS !< The fast ice thermodynamics control structure

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: m_ice_tot !< The total mass of ice in a cell [R Z ~> kg m-2]
  real, dimension(IG%NkIce) :: sw_abs_lay ! The fraction of the absorbed shortwave that is
                                          ! absorbed in each of the ice layers, <=1, [nondim].
  real, dimension(size(FIA%flux_sw_top,4)) :: &
    albedos        ! The albedos for the various wavelength and direction bands
                   ! for the current partition, non-dimensional and 0 to 1.
  real :: u, v     ! Ice velocity components [m s-1]

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  integer :: index
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc

  if (fCS%bounds_check) &
    call IST_bounds_check(IST, G, US, IG, "Start of set_ice_surface_state", Rad=Rad) !, OSS=OSS)

  if (fCS%debug) then
    call IST_chksum("Start set_ice_surface_state", IST, G, fCS%US, IG)
    call Ice_public_type_chksum("Start set_ice_surface_state", Ice, &
                                check_fast=.false., check_rough=.true.)
  endif

  m_ice_tot(:,:) = 0.0
  !$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,IST,OSS,FIA,ncat,m_ice_tot)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    FIA%tmelt(i,j,k) = 0.0 ; FIA%bmelt(i,j,k) = 0.0
    m_ice_tot(i,j) = m_ice_tot(i,j) + IST%mH_ice(i,j,k) * IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  if (.not.fCS%Eulerian_tsurf) then
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      Rad%t_skin(i,j,k) = IST%t_surf(i,j,k) - IST%T_0degC
    enddo ; enddo ; enddo
  endif

  ! Determine the sea-ice optical properties.

  !   These initialization calls for ice-free categories are not really
  ! needed because these arrays are only used where there is ice.
  !   The following lines can be uncommented without changing answers.
  ! Rad%sw_abs_sfc(:,:,:) = 0.0 ; Rad%sw_abs_snow(:,:,:) = 0.0
  ! Rad%sw_abs_ice(:,:,:,:) = 0.0 ; Rad%sw_abs_ocn(:,:,:) = 0.0
  ! Rad%sw_abs_int(:,:,:) = 0.0
  ! Ice%albedo(:,:,:) = 0.0
  ! Ice%albedo_vis_dir(:,:,:) = 0.0 ; Ice%albedo_vis_dif(:,:,:) = 0.0
  ! Ice%albedo_nir_dir(:,:,:) = 0.0 ; Ice%albedo_nir_dif(:,:,:) = 0.0

  ! Set the initial ocean albedos, either using coszen_nextrad or a synthetic sun angle.
  if (Rad%do_sun_angle_for_alb) then
    call set_ocean_albedo_from_astronomy(Ice, G, fCS%Time, fCS%Time + fCS%Time_step_fast)
  else
    call set_ocean_albedo_from_coszen(Ice, G, Rad%coszen_nextrad)
  endif

  !$OMP parallel do default(shared) private(i2,j2,k2,sw_abs_lay,albedos)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    if (IST%part_size(i,j,k)*IST%MH_ice(i,j,k) > 0.0) then
      call ice_optics_SIS2(IST%mH_pond(i,j,k), IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), &
               Rad%t_skin(i,j,k), OSS%T_fr_ocn(i,j), IG%NkIce, albedos, &
               Rad%sw_abs_sfc(i,j,k), Rad%sw_abs_snow(i,j,k), &
               sw_abs_lay, Rad%sw_abs_ocn(i,j,k), Rad%sw_abs_int(i,j,k), &
               US, fCS%optics_CSp, IST%ITV, coszen_in=Rad%coszen_nextrad(i,j))
      Ice%albedo_vis_dir(i2,j2,k2) = albedos(VIS_DIR)
      Ice%albedo_vis_dif(i2,j2,k2) = albedos(VIS_DIF)
      Ice%albedo_nir_dir(i2,j2,k2) = albedos(NIR_DIR)
      Ice%albedo_nir_dif(i2,j2,k2) = albedos(NIR_DIF)

      do m=1,IG%NkIce ; Rad%sw_abs_ice(i,j,k,m) = sw_abs_lay(m) ; enddo

      !Niki: Is the following correct for diagnostics?
      !###  This calculation of the "average" albedo should be replaced
      ! with a calculation that weights the averaging by the fraction of the
      ! shortwave radiation in each wavelength and orientation band.  However,
      ! since this is only used for diagnostic purposes, making this change
      ! might not be too urgent. -RWH
      Ice%albedo(i2,j2,k2) = (Ice%albedo_vis_dir(i2,j2,k2)+Ice%albedo_nir_dir(i2,j2,k2)&
                        +Ice%albedo_vis_dif(i2,j2,k2)+Ice%albedo_nir_dif(i2,j2,k2))/4
    else ! Zero the albedos and absorbed shortwave radiation out for debugging purposes.
      Rad%sw_abs_sfc(i,j,k) = 0.0
      Rad%sw_abs_snow(i,j,k) = 0.0
      Rad%sw_abs_int(i,j,k) = 0.0
      Rad%sw_abs_ocn(i,j,k) = 0.0
      Ice%albedo_vis_dir(i2,j2,k2) = 0.0
      Ice%albedo_vis_dif(i2,j2,k2) = 0.0
      Ice%albedo_nir_dir(i2,j2,k2) = 0.0
      Ice%albedo_nir_dif(i2,j2,k2) = 0.0
      Ice%albedo(i2,j2,k2) = 0.0
    endif
  enddo ; enddo ; enddo

  !$OMP parallel do default(shared)
  do j=jsc,jec
    do k=1,ncat ; do i=isc,iec
      Rad%Tskin_rad(i,j,k) = Rad%t_skin(i,j,k)
    enddo ; enddo
    do i=isc,iec
      Rad%coszen_lastrad(i,j) = Rad%coszen_nextrad(i,j)
    enddo
  enddo

  if (fCS%bounds_check) &
    call IST_bounds_check(IST, G, US, IG, "Midpoint set_ice_surface_state", Rad=Rad) !, OSS=OSS)

  ! Copy the surface temperatures into the externally visible data type.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,ncat,i_off,j_off,OSS,US) &
!$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    Ice%t_surf(i2,j2,1) = US%C_to_degC*OSS%SST_C(i,j) + T_0degC
    Ice%part_size(i2,j2,1) = IST%part_size(i,j,0)
  enddo ; enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Rad,Ice,ncat,i_off,j_off,OSS,US) &
!$OMP                          private(i2,j2,k2)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      Ice%t_surf(i2,j2,k2) = US%C_to_degC*Rad%t_skin(i,j,k) + T_0degC
      Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
    enddo ; enddo
  enddo

  ! Put ocean salinity and ocean and ice velocities into Ice%u_surf/v_surf on t-cells.
  !$OMP parallel do default(shared) private(i2,j2)
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%u_surf(i2,j2,1) = US%L_T_to_m_s*OSS%u_ocn_A(i,j)
    Ice%v_surf(i2,j2,1) = US%L_T_to_m_s*OSS%v_ocn_A(i,j)
    Ice%u_surf(i2,j2,2) = US%L_T_to_m_s*OSS%u_ice_A(i,j)
    Ice%v_surf(i2,j2,2) = US%L_T_to_m_s*OSS%v_ice_A(i,j)
    Ice%s_surf(i2,j2) = US%S_to_ppt*OSS%s_surf(i,j)
  enddo ; enddo

  if (fCS%debug) then
    call chksum(Ice%u_surf(:,:,1), "Intermed Ice%u_surf(1)")
    call chksum(Ice%v_surf(:,:,1), "Intermed Ice%v_surf(1)")
    call chksum(Ice%u_surf(:,:,2), "Intermed Ice%u_surf(2)")
    call chksum(Ice%v_surf(:,:,2), "Intermed Ice%v_surf(2)")
    call chksum(G%mask2dT(isc:iec,jsc:jec), "Intermed G%mask2dT")
!   if (allocated(OSS%u_ocn_C) .and. allocated(OSS%v_ocn_C)) &
!     call uvchksum(OSS%u_ocn_C, "OSS%u_ocn_C", &
!                   OSS%v_ocn_C, "OSS%v_ocn_C", G%HI, haloshift=1, scale=US%L_T_to_m_s)
!   if (allocated(OSS%u_ocn_B)) &
!     call Bchksum(OSS%u_ocn_B, "OSS%u_ocn_B", G%HI, haloshift=1, scale=US%L_T_to_m_s)
!   if (allocated(OSS%v_ocn_B)) &
!     call Bchksum(OSS%v_ocn_B, "OSS%v_ocn_B", G%HI, haloshift=1)
    call chksum(G%sin_rot(isc:iec,jsc:jec), "G%sin_rot")
    call chksum(G%cos_rot(isc:iec,jsc:jec), "G%cos_rot")
  endif

  ! Rotate the velocities from the ocean coordinates to lat/lon coordinates.
  do k2=1,2 ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    u = Ice%u_surf(i2,j2,k2) ; v = Ice%v_surf(i2,j2,k2)
    Ice%u_surf(i2,j2,k2) =  u*G%cos_rot(i,j)+v*G%sin_rot(i,j)
    Ice%v_surf(i2,j2,k2) =  v*G%cos_rot(i,j)-u*G%sin_rot(i,j)
  enddo ; enddo ; enddo
  do k2=3,ncat+1
    Ice%u_surf(:,:,k2) = Ice%u_surf(:,:,2)  ! same ice flow on all ice partitions
    Ice%v_surf(:,:,k2) = Ice%v_surf(:,:,2)  !
  enddo
  if (fCS%debug) then
    do k2=1,ncat+1
      call chksum(Ice%u_surf(:,:,k2), "End Ice%u_surf(k2)")
      call chksum(Ice%v_surf(:,:,k2), "End Ice%v_surf(k2)")
    enddo
  endif

  ! Copy over the additional tracer fields into the open-ocean category of
  ! the Ice%ocean_fields structure.
  call coupler_type_copy_data(OSS%tr_fields, Ice%ocean_fields, ind3_start=1, ind3_end=1)

  if (fCS%debug) then
    call IST_chksum("End set_ice_surface_state", IST, G, fCS%US, IG)
    call Ice_public_type_chksum("End set_ice_surface_state", Ice, check_fast=.true.)
  endif

  if (fCS%bounds_check) &
    call Ice_public_type_bounds_check(Ice, G, "End set_ice_surface_state")

end subroutine set_ice_surface_state

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> update_ice_model_fast records fluxes (in Ice) and calculates ice temperature
!!    on the (fast) atmospheric timestep
subroutine update_ice_model_fast( Atmos_boundary, Ice )
  type(ice_data_type),           intent(inout) :: Ice !< The publicly visible ice data type.
  type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary !< A type containing atmospheric boundary
                                                !! forcing fields that are used to drive the ice

  type(time_type) :: Time_start, Time_end, dT_fast

  call cpu_clock_begin(iceClock) ; call cpu_clock_begin(ice_clock_fast)

  if (Ice%fCS%debug) &
    call Ice_public_type_chksum("Pre do_update_ice_model_fast", Ice, check_fast=.true.)

  dT_fast = Ice%fCS%Time_step_fast
  Time_start = Ice%fCS%Time
  Time_end = Time_start + dT_fast

  if (Ice%fCS%Rad%add_diurnal_sw) &
    call add_diurnal_sw(Atmos_boundary, Ice%fCS%G, Time_start, Time_end)

  call do_update_ice_model_fast(Atmos_boundary, Ice%fCS%IST, Ice%fCS%sOSS, Ice%fCS%Rad, &
                                Ice%fCS%FIA, dT_fast, Ice%fCS%fast_thermo_CSp, &
                                Ice%fCS%G, Ice%fCS%US, Ice%fCS%IG )

  ! Advance the master sea-ice time.
  Ice%fCS%Time = Ice%fCS%Time + dT_fast

  Ice%Time = Ice%fCS%Time

  call fast_radiation_diagnostics(Atmos_boundary, Ice, Ice%fCS%IST, Ice%fCS%Rad, Ice%fCS%FIA, &
                                  Ice%fCS%G, Ice%fCS%US, Ice%fCS%IG, Ice%fCS, Time_start, Time_end)

  ! Set some of the evolving ocean properties that will be seen by the
  ! atmosphere in the next time-step.
  call set_fast_ocean_sfc_properties(Atmos_boundary, Ice, Ice%fCS%IST, Ice%fCS%Rad, Ice%fCS%FIA, &
                                     Ice%fCS%G, Ice%fCS%US, Ice%fCS%IG, Time_end, Time_end + dT_fast)

  if (Ice%fCS%debug) &
    call Ice_public_type_chksum("End do_update_ice_model_fast", Ice, check_fast=.true.)
  if (Ice%fCS%bounds_check) &
    call Ice_public_type_bounds_check(Ice, Ice%fCS%G, "End update_ice_fast")

  call cpu_clock_end(ice_clock_fast) ; call cpu_clock_end(iceClock)

end subroutine update_ice_model_fast

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_fast_ocean_sfc_properties updates the ocean surface properties like
!! roughness and albedo for the rapidly evolving atmospheric updates
subroutine set_fast_ocean_sfc_properties( Atmos_boundary, Ice, IST, Rad, FIA, &
                                          G, US, IG, Time_start, Time_end)
  type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary !< A type containing atmospheric boundary
                                                      !! forcing fields that are used to drive the ice
  type(ice_data_type),           intent(inout) :: Ice !< The publicly visible ice data type.
  type(ice_state_type),          intent(inout) :: IST !< A type describing the state of the sea ice
  type(ice_rad_type),            intent(inout) :: Rad !< A structure with fields related to the absorption,
                                                      !! reflection and transmission of shortwave radiation.
  type(fast_ice_avg_type),       intent(inout) :: FIA !< A type containing averages of fields
                                                      !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),       intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),         intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),           intent(inout) :: IG  !< The sea-ice specific grid type
  type(time_type),               intent(in)    :: Time_start !< The start of the time covered by this call
  type(time_type),               intent(in)    :: Time_end   !< The end of the time covered by this call

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  logical :: coszen_changed
  integer :: i, j, k, i2, j2, k2, i3, j3, isc, iec, jsc, jec, ncat
  integer :: io_A, jo_A, io_I, jo_I  ! Offsets for indexing conventions.

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  io_A = LBOUND(Atmos_boundary%t_flux,1) - G%isc
  jo_A = LBOUND(Atmos_boundary%t_flux,2) - G%jsc
  io_I = LBOUND(Ice%t_surf,1) - G%isc
  jo_I = LBOUND(Ice%t_surf,2) - G%jsc

  call compute_ocean_roughness (Ice%ocean_pt, Atmos_boundary%u_star(:,:,1), Ice%rough_mom(:,:,1), &
                                Ice%rough_heat(:,:,1), Ice%rough_moist(:,:,1)  )

  ! Update publicly visible ice_data_type variables..
  coszen_changed = .false.
  !$OMP parallel do default(shared) private(i3,j3)
  do j=jsc,jec ; do i=isc,iec
    i3 = i+io_A ; j3 = j+jo_A
    Rad%coszen_nextrad(i,j) = Atmos_boundary%coszen(i3,j3,1)
    FIA%p_atm_surf(i,j) = US%kg_m2s_to_RZ_T*US%m_s_to_L_T*Atmos_boundary%p(i3,j3,1)
    if (Rad%coszen_nextrad(i,j) /= Rad%coszen_lastrad(i,j)) coszen_changed = .true.
  enddo ; enddo

  !$OMP parallel do default(shared) private(i2,j2,k2)
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+io_I ; j2 = j+jo_I ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = US%C_to_degC*Rad%t_skin(i,j,k) + T_0degC
  enddo ; enddo ; enddo

  ! set_ocean_albedo only needs to be called if do_sun_angle_for_alb is true or
  ! if the coupled model's radiation timestep is shorter than the slow coupling
  ! timestep.  However, it is safe (if wasteful) to call it more frequently.
  if (Rad%do_sun_angle_for_alb) then
    call set_ocean_albedo_from_astronomy(Ice, G, Time_start, Time_end)
  elseif (coszen_changed) then
    call set_ocean_albedo_from_coszen(Ice, G, Rad%coszen_nextrad)
  endif

end subroutine set_fast_ocean_sfc_properties

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ocean_albedo uses the time and astronomical calculates to set the solar
!! zenith angle and calculate the ocean albedo.
subroutine set_ocean_albedo_from_astronomy(Ice, G, Time_start, Time_end)
  type(ice_data_type),     intent(inout) :: Ice !< The publicly visible ice data type.
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(time_type),         intent(in)    :: Time_start !< The start of the time covered by this call
  type(time_type),         intent(in)    :: Time_end   !< The end of the time covered by this call

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: &
    dummy, &  ! A dummy array that is not used again.
    cosz_alb  ! The cosine of the solar zenith angle for calculating albedo [nondim].
  real :: rad
  real :: rrsun_dt_ice
  type(time_type) :: dT_ice   ! The time interval for this update.
  integer :: i, j, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  rad = acos(-1.)/180.
  dT_ice = Time_end - Time_start

  call diurnal_solar(G%geoLatT(isc:iec,jsc:jec)*rad, G%geoLonT(isc:iec,jsc:jec)*rad, &
               Time_start, cosz=cosz_alb, fracday=dummy, rrsun=rrsun_dt_ice, &
               dt_time=dT_ice)

  call compute_ocean_albedo(Ice%ocean_pt, cosz_alb(:,:), Ice%albedo_vis_dir(:,:,1), &
                            Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1), &
                            Ice%albedo_nir_dif(:,:,1), rad*G%geoLatT(isc:iec,jsc:jec), &
                            Ice%flux_u(:,:), Ice%flux_v(:,:) )

end subroutine set_ocean_albedo_from_astronomy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ocean_albedo uses either the input cosine of the solar zenith
!! angle to calculate the ocean albedo.
subroutine set_ocean_albedo_from_coszen(Ice, G, coszen)
  type(ice_data_type),     intent(inout) :: Ice !< The publicly visible ice data type.
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  real, dimension(G%isd:G%ied, G%jsd:G%jed), &
                           intent(in)    :: coszen !< Cosine of the solar zenith angle for this step

  real :: rad
  integer :: i, j, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  rad = acos(-1.)/180.

  call compute_ocean_albedo(Ice%ocean_pt, coszen(isc:iec,jsc:jec), Ice%albedo_vis_dir(:,:,1),&
                            Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                            Ice%albedo_nir_dif(:,:,1), rad*G%geoLatT(isc:iec,jsc:jec), Ice%flux_u(:,:), Ice%flux_v(:,:))

end subroutine set_ocean_albedo_from_coszen

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> fast_radiation_diagnostics offers diagnostics of the rapidly changing shortwave
!! radiative and other properties of the ice, and it accumulates the shortwave radiation.
subroutine fast_radiation_diagnostics(ABT, Ice, IST, Rad, FIA, G, US, IG, CS, &
                                      Time_start, Time_end)
  type(atmos_ice_boundary_type), &
                           intent(in)    :: ABT !< A type containing atmospheric boundary
                                                !! forcing fields that are used to drive the ice
  type(ice_data_type),     intent(in)    :: Ice !< The publicly visible ice data type.
  type(ice_state_type),    intent(in)    :: IST !< A type describing the state of the sea ice
  type(ice_rad_type),      intent(in)    :: Rad !< A structure with fields related to the absorption,
                                                !! reflection and transmission of shortwave radiation.
  type(fast_ice_avg_type), intent(inout) :: FIA !< A type containing averages of fields
                                                !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type), intent(in)    :: G   !< The horizontal grid type
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  type(SIS_fast_CS),       intent(inout) :: CS  !< The fast ice thermodynamics control structure
  type(time_type),         intent(in)    :: Time_start !< The start time of the diagnostics in this call
  type(time_type),         intent(in)    :: Time_end   !< The end time of the diagnostics in this call

  real, dimension(G%isd:G%ied, G%jsd:G%jed) :: tmp_diag, sw_dn, net_sw, avg_alb
  real, dimension(G%isd:G%ied, G%jsd:G%jed,size(FIA%flux_sw_dn,3)) :: &
    sw_dn_bnd  ! The downward shortwave radiation by frequency and angular band
               ! averaged over all of the ice thickness categories [W m-2].
  real, dimension(G%isd:G%ied) :: Tskin_avg ! Average skin temperature [C ~> degC]
  real, dimension(G%isd:G%ied) :: ice_conc
  real :: dt_diag
  real    :: Stefan ! The Stefan-Boltzmann constant [W m-2 degK-4] as used for
                    ! strictly diagnostic purposes.
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  integer :: i, j, k, m, i2, j2, k2, i3, j3, isc, iec, jsc, jec, ncat, NkIce
  integer :: b
  integer :: io_A, jo_A, io_I, jo_I  ! Offsets for indexing conventions.

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  NkIce = IG%NkIce
  io_A = LBOUND(ABT%lw_flux,1) - G%isc ; jo_A = LBOUND(ABT%lw_flux,2) - G%jsc
  io_I = LBOUND(Ice%albedo_vis_dir,1) - G%isc
  jo_I = LBOUND(Ice%albedo_vis_dir,2) - G%jsc

  dt_diag = time_type_to_real(Time_end - Time_start)

  call enable_SIS_averaging(dt_diag, Time_end, CS%diag)

  if (Rad%id_alb_vis_dir>0) call post_avg(Rad%id_alb_vis_dir, Ice%albedo_vis_dir, &
                              IST%part_size(isc:iec,jsc:jec,:), CS%diag)
  if (Rad%id_alb_vis_dif>0) call post_avg(Rad%id_alb_vis_dif, Ice%albedo_vis_dif, &
                              IST%part_size(isc:iec,jsc:jec,:), CS%diag)
  if (Rad%id_alb_nir_dir>0) call post_avg(Rad%id_alb_nir_dir, Ice%albedo_nir_dir, &
                              IST%part_size(isc:iec,jsc:jec,:), CS%diag)
  if (Rad%id_alb_nir_dif>0) call post_avg(Rad%id_alb_nir_dif, Ice%albedo_nir_dif, &
                              IST%part_size(isc:iec,jsc:jec,:), CS%diag)
  if (Rad%id_alb>0)         call post_avg(Rad%id_alb, Ice%albedo, &
                              IST%part_size(isc:iec,jsc:jec,:), CS%diag)

  if (Rad%id_tskin>0) call post_data(Rad%id_tskin, Rad%t_skin, CS%diag)
  if (Rad%id_cn>0) call post_data(Rad%id_cn, IST%part_size(:,:,1:), CS%diag)
  if (Rad%id_mi>0) call post_data(Rad%id_mi, IST%mH_ice(:,:,:), CS%diag)

  if (Rad%id_sw_abs_sfc>0) call post_avg(Rad%id_sw_abs_sfc, Rad%sw_abs_sfc, &
                                   IST%part_size(:,:,1:), CS%diag, G=G)
  if (Rad%id_sw_abs_snow>0) call post_avg(Rad%id_sw_abs_snow, Rad%sw_abs_snow, &
                                   IST%part_size(:,:,1:), CS%diag, G=G)
  if (allocated(Rad%id_sw_abs_ice)) then ; do m=1,NkIce
    if (Rad%id_sw_abs_ice(m)>0) call post_avg(Rad%id_sw_abs_ice(m), Rad%sw_abs_ice(:,:,:,m), &
                                     IST%part_size(:,:,1:), CS%diag, G=G)
  enddo ; endif
  if (Rad%id_sw_abs_ocn>0) call post_avg(Rad%id_sw_abs_ocn, Rad%sw_abs_ocn, &
                                   IST%part_size(:,:,1:), CS%diag, G=G)

  if (Rad%id_sw_pen>0) then
    tmp_diag(:,:) = 0.0
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * &
                     (Rad%sw_abs_ocn(i,j,k) + Rad%sw_abs_int(i,j,k))
    enddo ; enddo ; enddo
    call post_data(Rad%id_sw_pen, tmp_diag, CS%diag)
  endif

  if (Rad%id_lwdn > 0) then
    tmp_diag(:,:) = 0.0
    Stefan = 5.6734e-8  ! Set the Stefan-Bolzmann constant [W m-2 degK-4].
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then
      i3 = i+io_A ; j3 = j+jo_A ; k2 = k+1
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * &
                           (ABT%lw_flux(i3,j3,k2) + Stefan*(US%C_to_degC*Rad%t_skin(i,j,k)+T_0degC)**4)
    endif ; enddo ; enddo ; enddo
    call post_data(Rad%id_lwdn, tmp_diag, CS%diag)
  endif

  sw_dn(:,:) = 0.0 ; net_sw(:,:) = 0.0 ; avg_alb(:,:) = 0.0 ; sw_dn_bnd(:,:,:) = 0.0
  !$OMP parallel do default(shared) private(i2,j2,k2,i3,j3)
  do j=jsc,jec ; do k=0,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then
    i2 = i+io_I ; j2 = j+jo_I ; i3 = i+io_A ; j3 = j+jo_A ; k2 = k+1
    if (associated(ABT%sw_down_vis_dir)) then
      sw_dn(i,j) = sw_dn(i,j) + IST%part_size(i,j,k) * ( &
            (ABT%sw_down_vis_dir(i3,j3,k2) + ABT%sw_down_vis_dif(i3,j3,k2)) + &
            (ABT%sw_down_nir_dir(i3,j3,k2) + ABT%sw_down_nir_dif(i3,j3,k2)) )
      sw_dn_bnd(i,j,VIS_DIR) = sw_dn_bnd(i,j,VIS_DIR) + &
                     IST%part_size(i,j,k) * ABT%sw_down_vis_dir(i3,j3,k2)
      sw_dn_bnd(i,j,VIS_DIF) = sw_dn_bnd(i,j,VIS_DIF) + &
                     IST%part_size(i,j,k) * ABT%sw_down_vis_dif(i3,j3,k2)
      sw_dn_bnd(i,j,NIR_DIR) = sw_dn_bnd(i,j,NIR_DIR) + &
                     IST%part_size(i,j,k) * ABT%sw_down_nir_dir(i3,j3,k2)
      sw_dn_bnd(i,j,NIR_DIF) = sw_dn_bnd(i,j,NIR_DIF) + &
                     IST%part_size(i,j,k) * ABT%sw_down_nir_dif(i3,j3,k2)
    else
      sw_dn(i,j) = sw_dn(i,j) + IST%part_size(i,j,k) * ( &
            (ABT%sw_flux_vis_dir(i3,j3,k2)/(1-Ice%albedo_vis_dir(i2,j2,k2)) + &
             ABT%sw_flux_vis_dif(i3,j3,k2)/(1-Ice%albedo_vis_dif(i2,j2,k2))) + &
            (ABT%sw_flux_nir_dir(i3,j3,k2)/(1-Ice%albedo_nir_dir(i2,j2,k2)) + &
             ABT%sw_flux_nir_dif(i3,j3,k2)/(1-Ice%albedo_nir_dif(i2,j2,k2))) )
      sw_dn_bnd(i,j,VIS_DIR) = sw_dn_bnd(i,j,VIS_DIR) + IST%part_size(i,j,k) * &
            (ABT%sw_flux_vis_dir(i3,j3,k2)/(1.0-Ice%albedo_vis_dir(i2,j2,k2)))
      sw_dn_bnd(i,j,VIS_DIF) = sw_dn_bnd(i,j,VIS_DIF) + IST%part_size(i,j,k) * &
            (ABT%sw_flux_vis_dif(i3,j3,k2)/(1.0-Ice%albedo_vis_dif(i2,j2,k2)))
      sw_dn_bnd(i,j,NIR_DIR) = sw_dn_bnd(i,j,NIR_DIR) + IST%part_size(i,j,k) * &
            (ABT%sw_flux_nir_dir(i3,j3,k2)/(1.0-Ice%albedo_nir_dir(i2,j2,k2)))
      sw_dn_bnd(i,j,NIR_DIF) = sw_dn_bnd(i,j,NIR_DIF) + IST%part_size(i,j,k) * &
            (ABT%sw_flux_nir_dif(i3,j3,k2)/(1.0-Ice%albedo_nir_dif(i2,j2,k2)))
    endif

    net_sw(i,j) = net_sw(i,j) + IST%part_size(i,j,k) * ( &
          (ABT%sw_flux_vis_dir(i3,j3,k2) + ABT%sw_flux_vis_dif(i3,j3,k2)) + &
          (ABT%sw_flux_nir_dir(i3,j3,k2) + ABT%sw_flux_nir_dif(i3,j3,k2)) )
    avg_alb(i,j) = avg_alb(i,j) + IST%part_size(i,j,k) * 0.25 * ( &
            (Ice%albedo_vis_dir(i2,j2,k2) + Ice%albedo_vis_dif(i2,j2,k2)) + &
            (Ice%albedo_nir_dir(i2,j2,k2) + Ice%albedo_nir_dif(i2,j2,k2)) )
    ! Consider recalculating this as avg_alb(i,j) = 1.0 - net_sw(i,j) / sw_dn(i,j) ? -RWH
  endif ; enddo ; enddo ; enddo

  !$OMP parallel do default(shared) private(Tskin_avg,ice_conc)
  do j=jsc,jec
    Tskin_avg(:) = 0.0 ; ice_conc(:) = 0.0
    do k=1,ncat ; do i=isc,iec
      Tskin_avg(i) = Tskin_avg(i) + Rad%t_skin(i,j,k) * IST%part_size(i,j,k)
      ice_conc(i) = ice_conc(i) + IST%part_size(i,j,k)
    enddo ; enddo
    do i=isc,iec
      if (ice_conc(i)>0.0) then
        FIA%Tskin_avg(i,j) = FIA%Tskin_avg(i,j) + (Tskin_avg(i) / ice_conc(i))
      ! else there is nothing to add, because Tskin_avg = 0.
      endif
    enddo
  enddo

  if (Rad%id_swdn > 0) call post_data(Rad%id_swdn, sw_dn, CS%diag)

  if (Rad%id_alb > 0) then
    do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then
      if (sw_dn(i,j) > 0.0) &
        avg_alb(i,j) = (sw_dn(i,j) - net_sw(i,j)) / sw_dn(i,j)
      ! Otherwise keep the simple average that was set above.
    endif ; enddo ; enddo
    call post_data(Rad%id_alb, avg_alb, CS%diag)
  endif

  do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then
    do b=1,size(FIA%flux_sw_dn,3)
      FIA%flux_sw_dn(i,j,b) = FIA%flux_sw_dn(i,j,b) + US%W_m2_to_QRZ_T*sw_dn_bnd(i,j,b)
    enddo
  endif ; enddo ; enddo

  if (Rad%id_coszen>0) call post_data(Rad%id_coszen, Rad%coszen_nextrad, CS%diag)

  call disable_SIS_averaging(CS%diag)

end subroutine fast_radiation_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> add_diurnal_SW adjusts the shortwave fluxes in an atmos_boundary_type variable
!! to add a synthetic diurnal cycle.
subroutine add_diurnal_SW(ABT, G, Time_start, Time_end)
  type(atmos_ice_boundary_type), intent(inout) :: ABT !< The type with atmospheric fluxes to be adjusted.
  type(SIS_hor_grid_type),       intent(in)    :: G   !< The sea-ice lateral grid type.
  type(time_type),               intent(in)    :: Time_start !< The start time for this step.
  type(time_type),               intent(in)    :: Time_end   !< The ending time for this step.

  real :: diurnal_factor, time_since_ae, rad
  real :: fracday_dt, fracday_day
  real :: cosz_day, cosz_dt, rrsun_day, rrsun_dt
  type(time_type) :: dt_here

  integer :: i, j, k, i2, j2, isc, iec, jsc, jec, ncat, i_off, j_off

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = size(ABT%sw_flux_nir_dir,3)
  i_off = LBOUND(ABT%t_flux,1) - G%isc ; j_off = LBOUND(ABT%t_flux,2) - G%jsc

  !   Orbital_time extracts the time of year relative to the northern
  ! hemisphere autumnal equinox from a time_type variable.
  time_since_ae = orbital_time(Time_start)
  dt_here = Time_end - Time_start
  rad = acos(-1.)/180.

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,rad,Time_start,dt_here,time_since_ae, &
!$OMP                                  ncat,ABT,i_off,j_off) &
!$OMP                          private(i,j,i2,j2,k,cosz_dt,fracday_dt,rrsun_dt, &
!$OMP                                  fracday_day,cosz_day,rrsun_day,diurnal_factor)
  do j=jsc,jec ; do i=isc,iec
!    Per Rick Hemler:
!      Call diurnal_solar with dtime=dt_here to get cosz averaged over dt_here.
!      Call daily_mean_solar to get cosz averaged over a day.  Then
!      diurnal_factor = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice /
!                       cosz_day*fracday_day*rrsun_day

    call diurnal_solar(G%geoLatT(i,j)*rad, G%geoLonT(i,j)*rad, Time_start, cosz=cosz_dt, &
                       fracday=fracday_dt, rrsun=rrsun_dt, dt_time=dt_here)
    call daily_mean_solar (G%geoLatT(i,j)*rad, time_since_ae, cosz_day, fracday_day, rrsun_day)
    diurnal_factor = cosz_dt*fracday_dt*rrsun_dt / &
                     max(1e-30, cosz_day*fracday_day*rrsun_day)

    i2 = i+i_off ; j2 = j+j_off
    do k=1,ncat
      ABT%sw_flux_nir_dir(i2,j2,k) = ABT%sw_flux_nir_dir(i2,j2,k) * diurnal_factor
      ABT%sw_flux_nir_dif(i2,j2,k) = ABT%sw_flux_nir_dif(i2,j2,k) * diurnal_factor
      ABT%sw_flux_vis_dir(i2,j2,k) = ABT%sw_flux_vis_dir(i2,j2,k) * diurnal_factor
      ABT%sw_flux_vis_dif(i2,j2,k) = ABT%sw_flux_vis_dif(i2,j2,k) * diurnal_factor
    enddo
    if (associated(ABT%sw_down_nir_dir)) then ; do k=1,ncat
      ABT%sw_down_nir_dir(i2,j2,k) = ABT%sw_down_nir_dir(i2,j2,k) * diurnal_factor
    enddo ; endif
    if (associated(ABT%sw_down_nir_dif)) then ; do k=1,ncat
      ABT%sw_down_nir_dif(i2,j2,k) = ABT%sw_down_nir_dif(i2,j2,k) * diurnal_factor
    enddo ; endif
    if (associated(ABT%sw_down_vis_dir)) then ; do k=1,ncat
      ABT%sw_down_vis_dir(i2,j2,k) = ABT%sw_down_vis_dir(i2,j2,k) * diurnal_factor
    enddo ; endif
    if (associated(ABT%sw_down_vis_dif)) then ; do k=1,ncat
      ABT%sw_down_vis_dif(i2,j2,k) = ABT%sw_down_vis_dif(i2,j2,k) * diurnal_factor
    enddo ; endif
  enddo ; enddo

end subroutine add_diurnal_sw

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_model_init - initializes ice model data, parameters and diagnostics. It
!! might operate on the fast ice processors, the slow ice processors or both.
subroutine ice_model_init(Ice, Time_Init, Time, Time_step_fast, Time_step_slow, &
                          Verona_coupler, Concurrent_atm, Concurrent_ice, gas_fluxes, gas_fields_ocn )

  type(ice_data_type), intent(inout) :: Ice            !< The ice data type that is being initialized.
  type(time_type)    , intent(in)    :: Time_Init      !< The starting time of the model integration
  type(time_type)    , intent(in)    :: Time           !< The current time
  type(time_type)    , intent(in)    :: Time_step_fast !< The time step for the ice_model_fast
  type(time_type)    , intent(in)    :: Time_step_slow !< The time step for the ice_model_slow
  logical,   optional, intent(in)    :: Verona_coupler !< If false or not present, use the input values
                                              !! in Ice to determine whether this is a fast or slow
                                              !! ice processor or both.  SIS2 will now throw a fatal
                                              !! error if this is present and true.
  logical,   optional, intent(in)    :: Concurrent_atm
  logical,   optional, intent(in)    :: Concurrent_ice !< If present and true, use sea ice model
                                              !! settings appropriate for running the atmosphere and
                                              !! slow ice simultaneously, including embedding the
                                              !! slow sea-ice time stepping in the ocean model.
  type(coupler_1d_bc_type), &
             optional, intent(in)     :: gas_fluxes !< If present, this type describes the
                                              !! additional gas or other tracer fluxes between the
                                              !! ocean, ice, and atmosphere, and can be used to
                                              !! spawn related internal variables in the ice model.
  type(coupler_1d_bc_type), &
             optional, intent(in)     :: gas_fields_ocn !< If present, this type describes the
                                              !! ocean and surface-ice fields that will participate
                                              !! in the calculation of additional gas or other
                                              !! tracer fluxes, and can be used to spawn related
                                              !! internal variables in the ice model.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  real :: enth_spec_snow, enth_spec_ice
  real :: pi ! pi = 3.1415926... calculated as 4*atan(1)
  integer :: i, j, k, l, i2, j2, k2, i_off, j_off, n
  integer :: isc, iec, jsc, jec, nCat_dflt
  character(len=120) :: restart_file, fast_rest_file
  character(len=40)  :: mdl = "ice_model" ! This module's name.
  character(len=8)   :: nstr
  type(directories)  :: dirs   ! A structure containing several relevant directory paths.

  type(param_file_type) :: param_file !< A structure to parse for run-time parameters
  type(hor_index_type)  :: fHI  !  A hor_index_type for array extents on fast_ice_PEs.
  type(hor_index_type)  :: sHI  !  A hor_index_type for array extents on slow_ice_PEs.

  type(dyn_horgrid_type),  pointer :: dG => NULL()
  type(unit_scale_type),   pointer :: US => NULL()
  ! These pointers are used only for coding convenience on slow PEs.
  type(SIS_hor_grid_type), pointer :: sG => NULL()
  type(MOM_domain_type),   pointer :: sGD => NULL()
  type(ice_state_type),    pointer :: sIST => NULL()
  type(ice_grid_type),     pointer :: sIG => NULL()

  ! These pointers are used only for coding convenience on fast PEs.
  type(SIS_hor_grid_type), pointer :: fG => NULL()
  type(MOM_domain_type), pointer :: fGD => NULL()

  ! Parameters that are read in and used to initialize other modules.  If those
  ! other modules had control states, these would be moved to those modules.
  real :: mom_rough_ice  ! momentum same, cd10=(von_k/ln(10/z0))^2 [m].
  real :: heat_rough_ice ! heat roughness length [m].
  real :: dt_Rad_real    ! The radiation timestep [s].
  type(time_type) :: dt_Rad ! The radiation timestep, used initializing albedos.
  real :: rad            ! The conversion factor from degrees to radians.
  real :: rrsun          ! An unused temporary factor related to the Earth-sun distance.

  ! Parameters that properly belong exclusively to ice_thm.
  real :: massless_ice_enth, massless_snow_enth ! Enthalpy fill values [Q ~> J kg-1]
  real :: massless_ice_salin  ! A salinity fill value [S ~> ppt]

  real, allocatable, dimension(:,:) :: &
    h_ice_input, dummy   ! Temporary arrays.
  real, allocatable, dimension(:,:) :: &
    str_x, str_y, stress_mag ! Temporary stress arrays

  real :: g_Earth        !   The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: ice_bulk_salin ! The globally constant sea ice bulk salinity [S ~> gSalt kg-1] = [S ~> ppt]
                         ! that is used to calculate the ocean salt flux.
  real :: ice_rel_salin  ! The initial bulk salinity of sea-ice relative to the
                         ! salinity of the water from which it formed [nondim].
  real :: coszen_IC      ! A constant value that is used to initialize
                         ! coszen if it is not read from a restart file, or a
                         ! negative number to use the time and geometry.
  real :: rho_ice        ! The nominal density of sea ice [R ~> kg m-3].
  real :: rho_snow       ! The nominal density of snow [R ~> kg m-3].
  real :: rho_Ocean      ! The nominal density of seawater [R ~> kg m-3].
  real :: kmelt          ! A constant that is used in the calculation of the ocean/ice basal heat
                         ! flux [Q R Z T-1 C-1 ~> W m-2 degC-1]. This could be changed to reflect
                         ! the turbulence in the under-ice ocean boundary layer and the effective
                         ! depth of the reported value of t_ocn.
  real :: opm_dflt       ! The default value for ocean_part_min, which is the
                         ! minimum value for the fractional open-ocean
                         ! area.  With redo_fast_update, this is set to 1e-40 so
                         ! that the open ocean fluxes can be used in
                         ! interpolation across categories; otherwise it is 0
                         ! to reduce the number of categories being evaluated.
                         ! In ice/ocean models, sufficiently small values do not
                         ! change answers, but in coupled models they do change
                         ! answers at roundoff.

  integer :: CatIce, NkIce, isd, ied, jsd, jed
  integer :: write_geom     ! A flag indicating whether to write the grid geometry files only for
                            ! new runs (1), both new runs and restarts (2) or neither (0).
  logical :: nudge_sea_ice  ! If true, nudge sea ice concentrations towards observations.
  logical :: transmute_ice  ! If true, allow ice to be transmuted directly into seawater with a
                            ! spatially varying rate as a form of outflow open boundary condition.
  logical :: atmos_winds, slp2ocean
  logical :: do_icebergs, pass_iceberg_area_to_ocean
  logical :: pass_stress_mag
  logical :: do_ridging
  logical :: specified_ice    ! If true, the ice is specified and there is no dynamics.
  logical :: Cgrid_dyn
  logical :: new_sim     ! If true, this is a new simulation, based on the contents of dirs and
                         ! the presence or absence of a named restart file.
  logical :: slab_ice    ! If true, use the very old slab ice thermodynamics,
                         ! with effectively zero heat capacity of ice and snow.
  logical :: debug, debug_slow, debug_fast, bounds_check
  logical :: do_sun_angle_for_alb, add_diurnal_sw
  logical :: init_coszen, init_Tskin, init_rough
  logical :: Eulerian_tsurf   ! If true, use previous calculations of the ice-top
                              ! surface skin temperature for tsurf at the start of
                              ! atmospheric time stepping, including interpolating between
                              ! tsurf values from other categories in the same location.
  logical :: write_geom_files ! If true, write out the grid geometry files.
  logical :: symmetric        ! If true, use symmetric memory allocation.
  logical :: global_indexing  ! If true use global horizontal index values instead
                              ! of having the data domain on each processor start at 1.
  integer :: first_direction  ! An integer that indicates which direction is to be
                              ! updated first in directionally split parts of the
                              ! calculation.  This can be altered during the course
                              ! of the run via calls to set_first_direction.
  logical :: fast_ice_PE      ! If true, fast ice processes are handled on this PE.
  logical :: slow_ice_PE      ! If true, slow ice processes are handled on this PE.
  logical :: single_IST       ! If true, fCS%IST and sCS%IST point to the same structure.
  logical :: interp_fluxes    ! If true, interpolate a linearized version of the
                              ! fast fluxes into arealess categories.
  logical :: redo_fast_update ! If true, recalculate the thermal updates from the fast
                              ! dynamics on the slowly evolving ice state, rather than
                              ! copying over the slow ice state to the fast ice state.
  logical :: do_mask_restart  ! If true, apply the scaling and masks to mH_snow, mH_ice, part_size
                              ! mH_pond, t_surf, t_skin, sal_ice, enth_ice and enth_snow
                              ! after a restart. However this may cause answers to diverge
                              ! after a restart.Provide a switch to turn this option off.
  logical :: recategorize_ice ! If true, adjust the distribution of the ice among thickness
                              ! categories after initialization.
  logical :: Verona
  logical :: Concurrent
  logical :: read_aux_restart
  logical :: split_restart_files
  logical :: is_restart = .false.
  character(len=16) :: stagger, dflt_stagger
  type(ice_OBC_type), pointer :: OBC_in => NULL()

  if (associated(Ice%sCS)) then ; if (associated(Ice%sCS%IST)) then
    call SIS_error(WARNING, "ice_model_init called with an associated "// &
                    "Ice%sCS%Ice_state structure. Model is already initialized.")
    return
  endif ; endif

  ! For now, both fast and slow processes occur on all sea-ice PEs.
  fast_ice_PE = .true. ; slow_ice_PE = .true.
  Verona = .false. ; if (present(Verona_coupler)) Verona = Verona_coupler
  if (Verona) call SIS_error(FATAL, "SIS2 no longer works with pre-Warsaw couplers.")
  fast_ice_PE = Ice%fast_ice_pe ; slow_ice_PE = Ice%slow_ice_pe
  Concurrent = .false. ; if (present(Concurrent_atm)) Concurrent = Concurrent_atm

  ! Open the parameter file.
  if (fast_ice_PE.eqv.slow_ice_PE) then
    call Get_SIS_Input(param_file, dirs, check_params=.true., component='SIS')  
    call Get_SIS_Input(param_file, dirs, check_params=.false., component='SIS_fast')
  elseif (slow_ice_PE) then
    call Get_SIS_Input(param_file, dirs, check_params=.true., component='SIS')
  elseif (fast_ice_PE) then
    call Get_SIS_Input(param_file, dirs, check_params=.false., component='SIS_fast')
  endif

  call callTree_enter("ice_model_init(), ice_model.F90")

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  ! Determining the internal unit scaling factors for this run.
  call unit_scaling_init(param_file, US)

  call get_param(param_file, mdl, "SPECIFIED_ICE", specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  call get_param(param_file, mdl, "CGRID_ICE_DYNAMICS", Cgrid_dyn, &
                 "If true, use a C-grid discretization of the sea-ice "//&
                 "dynamics; if false use a B-grid discretization.", &
                 default=.true.)
  if (specified_ice) then
    slab_ice = .true.
    call log_param(param_file, mdl, "USE_SLAB_ICE", slab_ice, &
                 "Use the very old slab-style ice.  With SPECIFIED_ICE, "//&
                 "USE_SLAB_ICE is always true.")
  else
    call get_param(param_file, mdl, "USE_SLAB_ICE", slab_ice, &
                 "If true, use the very old slab-style ice.", default=.false.)
  endif
  call get_param(param_file, mdl, "SINGLE_ICE_STATE_TYPE", single_IST, &
                 "If true, the fast and slow portions of the ice use a "//&
                 "single common ice_state_type.  Otherwise they point to "//&
                 "different ice_state_types that need to be explicitly "//&
                 "copied back and forth.", default=.true.)
  call get_param(param_file, mdl, "EULERIAN_TSURF", Eulerian_tsurf, &
                 "If true, use previous calculations of the ice-top surface "//&
                 "skin temperature for tsurf at the start of atmospheric "//&
                 "time stepping, including interpolating between tsurf "//&
                 "values from other categories in the same location.", default=.true.)

  call obsolete_logical(param_file, "SIS1_5L_THERMODYNAMICS", warning_val=.false.)
  call obsolete_logical(param_file, "INTERSPERSED_ICE_THERMO", warning_val=.false.)
  call obsolete_logical(param_file, "AREA_WEIGHTED_STRESSES", warning_val=.true.)

  dflt_stagger = "B" ; if (Cgrid_dyn) dflt_stagger = "C"
  call get_param(param_file, mdl, "ICE_OCEAN_STRESS_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the "//&
                 "staggering of the stress field on the ocean that is "//&
                 "returned to the coupler.  Valid values include "//&
                 "'A', 'B', or 'C', with a default that follows the "//&
                 "value of CGRID_ICE_DYNAMICS.", default=dflt_stagger)

  ! Rho_ocean is not actually used here, but it used from later get_param
  ! calls in other modules.  This call is here to avoid changing the order of
  ! the entries in the SIS_parameter_doc files.
  call get_param(param_file, mdl, "RHO_OCEAN", Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RHO_ICE", Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RHO_SNOW", Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0, scale=US%kg_m3_to_R)

  call get_param(param_file, mdl, "G_EARTH", g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80, scale=US%m_s_to_L_T**2*US%Z_to_m)

  call get_param(param_file, mdl, "MOMENTUM_ROUGH_ICE", mom_rough_ice, &
                 "The default momentum roughness length scale for the ocean.", &
                 units="m", default=1.0e-4)
  call get_param(param_file, mdl, "HEAT_ROUGH_ICE", heat_rough_ice, &
                 "The default roughness length scale for the turbulent "//&
                 "transfer of heat into the ocean.", units="m", default=1.0e-4)

  call get_param(param_file, mdl, "CONSTANT_COSZEN_IC", coszen_IC, &
                 "A constant value to use to initialize the cosine of "//&
                 "the solar zenith angle for the first radiation step, "//&
                 "or a negative number to use the current time and astronomy.", &
                 units="nondim", default=-1.0)
  call get_param(param_file, mdl, "DT_RADIATION", dt_Rad_real, &
                 "The time step with which the shortwave radiation and "//&
                 "fields like albedos are updated.  Currently this is only "//&
                 "used to initialize albedos when there is no restart file.", &
                 units="s", default=time_type_to_real(Time_step_slow))
  dt_Rad = real_to_time(dt_Rad_real)
  call get_param(param_file, mdl, "ICE_KMELT", kmelt, &
                 "A constant giving the proportionality of the ocean/ice "//&
                 "base heat flux to the tempature difference, given by "//&
                 "the product of the heat capacity per unit volume of sea "//&
                 "water times a molecular diffusive piston velocity.", &
                 units="W m-2 K-1", scale=US%W_m2_to_QRZ_T*US%C_to_degC, default=6e-5*4e6)
  call obsolete_real(param_file, "SNOW_CONDUCT", warning_val=0.31)
  call get_param(param_file, mdl, "ICE_BOUNDS_CHECK", bounds_check, &
                 "If true, periodically check the values of ice and snow "//&
                 "temperatures and thicknesses to ensure that they are "//&
                 "sensible, and issue warnings if they are not.  This "//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_SLOW_ICE", debug_slow, &
                 "If true, write out verbose debugging data on the slow ice PEs.", &
                 default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_FAST_ICE", debug_fast, &
                 "If true, write out verbose debugging data on the fast ice PEs.", &
                 default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "GLOBAL_INDEXING", global_indexing, &
                 "If true, use a global lateral indexing convention, so "//&
                 "that corresponding points on different processors have "//&
                 "the same index. This does not work with static memory.", &
                 default=.false., layoutParam=.true.)
#ifdef STATIC_MEMORY_
  if (global_indexing) call SIS_error(FATAL, "ice_model_init: "//&
       "GLOBAL_INDEXING can not be true with STATIC_MEMORY.")
#endif
  call get_param(param_file, mdl, "FIRST_DIRECTION", first_direction, &
                 "An integer that indicates which direction goes first "//&
                 "in parts of the code that use directionally split "//&
                 "updates, with even numbers (or 0) used for x- first "//&
                 "and odd numbers used for y-first.", default=0)

  call get_param(param_file, mdl, "ICE_SEES_ATMOS_WINDS", atmos_winds, &
                 "If true, the sea ice is being given wind stresses with "//&
                 "the atmospheric sign convention, and need to have their sign changed.", &
                 default=.true.)
  call get_param(param_file, mdl, "ICE_BULK_SALINITY", ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "g/kg", &
                 default=4.0, do_not_log=.true.)
  call get_param(param_file, mdl, "ICE_RELATIVE_SALINITY", ice_rel_salin, &
                 "The initial salinity of sea ice as a fraction of the "//&
                 "salinity of the seawater from which it formed.", &
                 units = "nondim", default=0.0, do_not_log=.true.)
  if ((ice_bulk_salin < 0.0) .or. (ice_rel_salin > 0.0)) ice_bulk_salin = 0.0

  call get_param(param_file, mdl, "APPLY_SLP_TO_OCEAN", slp2ocean, &
                 "If true, apply the atmospheric sea level pressure to the ocean.", &
                 default=.false.)
  call get_param(param_file, mdl, "PASS_STRESS_MAG_TO_OCEAN", pass_stress_mag, &
                 "If true, provide the time and area weighted mean magnitude "//&
                 "of the stresses on the ocean to the ocean.", default=.false.)
  call get_param(param_file, mdl, "DO_ICEBERGS", do_icebergs, &
                 "If true, call the iceberg module.", default=.false.)
  if (do_icebergs) then
    call get_param(param_file, mdl, "PASS_ICEBERG_AREA_TO_OCEAN", pass_iceberg_area_to_ocean, &
                 "If true, iceberg area is passed through coupler", default=.false.)
  else ; pass_iceberg_area_to_ocean = .false. ; endif

  call get_param(param_file, mdl, "ADD_DIURNAL_SW", add_diurnal_sw, &
                 "If true, add a synthetic diurnal cycle to the shortwave radiation.", &
                 default=.false.)
  call get_param(param_file, mdl, "DO_SUN_ANGLE_FOR_ALB", do_sun_angle_for_alb, &
                 "If true, find the sun angle for calculating the ocean "//&
                 "albedo within the sea ice model.", default=.false.)
  call get_param(param_file, mdl, "DO_RIDGING", do_ridging, &
                 "If true, call the ridging routines.", default=.false.)

  call get_param(param_file, mdl, "RESTARTFILE", restart_file, &
                 "The name of the restart file.", default="ice_model.res.nc")
  if (fast_ice_PE.eqv.slow_ice_PE) then
    call get_param(param_file, mdl, "FAST_ICE_RESTARTFILE", fast_rest_file, &
                   "The name of the restart file for those elements of the "//&
                   "the sea ice that are handled by the fast ice PEs.", default=restart_file)
  else
    call get_param(param_file, mdl, "FAST_ICE_RESTARTFILE", fast_rest_file, &
                   "The name of the restart file for those elements of the "//&
                   "the sea ice that are handled by the fast ice PEs.", &
                   default="ice_model_fast.res.nc")
  endif
  call get_param(param_file, mdl, "APPLY_MASKS_AFTER_RESTART", do_mask_restart, &
                 "If true, applies masks to mH_ice,mH_snow and part_size after a restart.", &
                  default=.true.)

  call get_param(param_file, mdl, "MASSLESS_ICE_ENTH", massless_ice_enth, &
                 "The ice enthalpy fill value for massless categories.", &
                 units="J kg-1", default=0.0, scale=US%J_kg_to_Q, do_not_log=.true.)
  call get_param(param_file, mdl, "MASSLESS_SNOW_ENTH", massless_snow_enth, &
                 "The snow enthalpy fill value for massless categories.", &
                 units="J kg-1", default=0.0, scale=US%J_kg_to_Q, do_not_log=.true.)
  call get_param(param_file, mdl, "MASSLESS_ICE_SALIN", massless_ice_salin, &
                 "The ice salinity fill value for massless categories.", &
                 units="g kg-1", default=0.0, scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, "MOM", "WRITE_GEOM", write_geom, &
                 "If =0, never write the geometry and vertical grid files. "//&
                 "If =1, write the geometry and vertical grid files only for "//&
                 "a new simulation. If =2, always write the geometry and "//&
                 "vertical grid files. Other values are invalid.", default=1)
  if (write_geom<0 .or. write_geom>2) call SIS_error(FATAL,"SIS2: "//&
         "WRITE_GEOM must be equal to 0, 1 or 2.")
  call get_param(param_file, "MOM", "INTERPOLATE_FLUXES", interp_fluxes, &
                 "If true, interpolate a linearized version of the fast "//&
                 "fluxes into arealess categories.", default=.true.)
  call get_param(param_file, "MOM", "REDO_FAST_ICE_UPDATE", redo_fast_update, &
                 "If true, recalculate the thermal updates from the fast "//&
                 "dynamics on the slowly evolving ice state, rather than "//&
                 "copying over the slow ice state to the fast ice state.", default=Concurrent_ice)

  call get_param(param_file, mdl, "NUDGE_SEA_ICE", nudge_sea_ice, &
                 "If true, constrain the sea ice concentrations using observations.", &
                 default=.false., do_not_log=.true.) ! Defer logging to SIS_slow_thermo.
  call get_param(param_file, mdl, "TRANSMUTE_SEA_ICE", transmute_ice, &
                 "If true, allow ice to be transmuted directly into seawater with a spatially "//&
                 "varying rate as a form of outflow open boundary condition.", &
                 default=.false., do_not_log=.true.) ! Defer logging to SIS_slow_thermo.


  nCat_dflt = 5 ; if (slab_ice) nCat_dflt = 1
  opm_dflt = 0.0 ; if (redo_fast_update) opm_dflt = 1.0e-40
#ifdef SYMMETRIC_MEMORY_
  symmetric = .true.
#else
  symmetric = .false.
#endif

  call SIS_debugging_init(param_file)

  ! Interpret and do error checking on some of the parameters.
  split_restart_files = (trim(restart_file) /= trim(fast_rest_file))
  if ((fast_ice_PE.neqv.slow_ice_PE) .and. .not.split_restart_files) then
    call SIS_error(FATAL, "The fast ice restart file must be separate from the "//&
           "standard ice restart file when there are separate fast and slow ice PEs. "//&
           "Choose different values of RESTARTFILE and FAST_ICE_RESTARTFILE.")
  endif

  if (fast_ice_PE.neqv.slow_ice_PE) single_IST = .false.

  if (uppercase(stagger(1:1)) == 'A') then ; Ice%flux_uv_stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; Ice%flux_uv_stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; Ice%flux_uv_stagger = CGRID_NE
  else ; call SIS_error(FATAL,"ice_model_init: ICE_OCEAN_STRESS_STAGGER = "//&
                        trim(stagger)//" is invalid.") ; endif

  Ice%Time = Time

  !   Now that all top-level sea-ice parameters have been read, allocate the
  ! various structures and register fields for restarts.
  if (slow_ice_PE) then
    if (.not.associated(Ice%sCS)) allocate(Ice%sCS)
    if (.not.associated(Ice%sCS%IG)) allocate(Ice%sCS%IG)
    if (.not.associated(Ice%sCS%IST)) allocate(Ice%sCS%IST)
    Ice%sCS%US => US
    Ice%sCS%Time = Time

    ! Set some pointers for convenience.
    sIST => Ice%sCS%IST ; sIG => Ice%sCS%IG
    sIST%Cgrid_dyn = Cgrid_dyn

    Ice%sCS%do_icebergs = do_icebergs
    Ice%sCS%pass_iceberg_area_to_ocean = pass_iceberg_area_to_ocean
    Ice%sCS%pass_stress_mag = pass_stress_mag
    Ice%sCS%slab_ice = slab_ice
    Ice%sCS%specified_ice = specified_ice
    Ice%sCS%Cgrid_dyn = Cgrid_dyn
    Ice%sCS%redo_fast_update = redo_fast_update
    Ice%sCS%bounds_check = bounds_check
    Ice%sCS%debug = debug_slow

    ! Set up the ice-specific grid describing categories and ice layers.
    call set_ice_grid(sIG, US, param_file, nCat_dflt, ocean_part_min_dflt=opm_dflt)
    if (slab_ice) sIG%CatIce = 1 ! open water and ice ... but never in same place
    CatIce = sIG%CatIce ; NkIce = sIG%NkIce
    call initialize_ice_categories(sIG, Rho_ice, US, param_file)


    ! Set up the domains and lateral grids.
    if (.not.associated(Ice%sCS%G)) allocate(Ice%sCS%G)
    sG => Ice%sCS%G

    ! Set up the MOM_domain_type structures.
#ifdef STATIC_MEMORY_
    call MOM_domains_init(Ice%sCS%G%domain, param_file, symmetric=symmetric, &
              static_memory=.true., NIHALO=NIHALO_, NJHALO=NJHALO_, &
              NIGLOBAL=NIGLOBAL_, NJGLOBAL=NJGLOBAL_, NIPROC=NIPROC_, &
              NJPROC=NJPROC_, domain_name="ice model", include_name="SIS2_memory.h")
#else
    call MOM_domains_init(Ice%sCS%G%domain, param_file, symmetric=symmetric, &
             domain_name="ice model", include_name="SIS2_memory.h")
#endif
    sGD => Ice%sCS%G%Domain

    call callTree_waypoint("domains initialized (ice_model_init)")
    call hor_index_init(sGD, sHI, param_file, &
                        local_indexing=.not.global_indexing)

    call create_dyn_horgrid(dG, sHI) !, bathymetry_at_vel=bathy_at_vel)
    call clone_MOM_domain(sGD, dG%Domain)

    ! Set up the restart file and determine whether this is a new simulation.
    call set_domain(sGD%mpp_domain)
    if (.not.associated(Ice%Ice_restart)) &
      call SIS_restart_init(Ice%Ice_restart, restart_file, sGD, param_file)
    new_sim = determine_is_new_run(dirs%input_filename, dirs%restart_input_dir, sG, Ice%Ice_restart)
    write_geom_files = ((write_geom==2) .or. ((write_geom==1) .and. new_sim))

    ! Set the bathymetry, Coriolis parameter, open channel widths and masks.
    call SIS_initialize_fixed(dG, US, param_file, write_geom_files, dirs%output_directory, OBC_in)

    call set_hor_grid(sG, param_file, global_indexing=global_indexing)
    call copy_dyngrid_to_SIS_horgrid(dG, sG)
    call destroy_dyn_horgrid(dG)
    Ice%sCS%G%US => US
    Ice%OBC => OBC_in

  ! Allocate and register fields for restarts.

    call ice_type_slow_reg_restarts(sGD%mpp_domain, CatIce, &
                      param_file, Ice, Ice%Ice_restart)

    call alloc_IST_arrays(sHI, sIG, US, sIST, omit_tsurf=Eulerian_tsurf, do_ridging=do_ridging)
    call ice_state_register_restarts(sIST, sG, sIG, US, Ice%Ice_restart)
    call register_unit_conversion_restarts(Ice%sCS%US, Ice%Ice_restart)

    call alloc_ocean_sfc_state(Ice%sCS%OSS, sHI, sIST%Cgrid_dyn, gas_fields_ocn)
    Ice%sCS%OSS%kmelt = kmelt

    call alloc_simple_OSS(Ice%sCS%sOSS, sHI, gas_fields_ocn)

    call alloc_ice_ocean_flux(Ice%sCS%IOF, sHI, do_stress_mag=Ice%sCS%pass_stress_mag, &
                              do_iceberg_fields=Ice%sCS%do_icebergs, do_transmute=transmute_ice)
    Ice%sCS%IOF%slp2ocean = slp2ocean
    Ice%sCS%IOF%flux_uv_stagger = Ice%flux_uv_stagger
    call alloc_fast_ice_avg(Ice%sCS%FIA, sHI, sIG, interp_fluxes, gas_fluxes)

    if (Ice%sCS%redo_fast_update) then
      call alloc_total_sfc_flux(Ice%sCS%TSF, sHI, gas_fluxes)
      call alloc_total_sfc_flux(Ice%sCS%XSF, sHI, gas_fluxes)
      call alloc_ice_rad(Ice%sCS%Rad, sHI, sIG)
    endif

    if (.not.specified_ice) &
      call SIS_dyn_trans_register_restarts(sHI, sIG, param_file, Ice%sCS%dyn_trans_CSp, US, &
                                           Ice%Ice_restart)

    call SIS_diag_mediator_init(sG, sIG, param_file, Ice%sCS%diag, component="SIS", &
                                doc_file_dir = dirs%output_directory)
    call set_SIS_axes_info(sG, sIG, param_file, Ice%sCS%diag)

    call ice_thermo_init(param_file, sIST%ITV, US, init_EOS=nudge_sea_ice)

    ! Register tracers that will be advected around.
    call register_SIS_tracer_pair(sIST%enth_ice, NkIce, "enth_ice", &
                                  sIST%enth_snow, 1, "enth_snow", &
                                  sG, sIG, param_file, sIST%TrReg, &
                                  massless_iceval=massless_ice_enth, &
                                  massless_snowval=massless_snow_enth)

    if (ice_rel_salin > 0.0) then
      call register_SIS_tracer(sIST%sal_ice, sG, sIG, NkIce, "salin_ice", param_file, &
                               sIST%TrReg, snow_tracer=.false., massless_val=massless_ice_salin, &
                               nonnegative=.true., conc_scale=US%S_to_ppt)
    endif

  !   Register any tracers that will be handled via tracer flow control for
  ! restarts and advection.
    call SIS_call_tracer_register(sG, sIG, param_file, Ice%sCS%SIS_tracer_flow_CSp, &
                                  Ice%sCS%diag, sIST%TrReg, Ice%Ice_restart)

    ! Set a few final things to complete the setup of the grid.
    sG%g_Earth = g_Earth
    call set_first_direction(sG, first_direction)
    call clone_MOM_domain(sGD, sG%domain_aux, symmetric=.false., &
                          domain_name="ice model aux")

    ! Copy the ice model's domain into one with no halos that can be shared
    ! publicly for use by the exchange grid.
    call clone_MOM_domain(sGD, Ice%slow_domain_NH, halo_size=0, symmetric=.false., &
                          domain_name="ice_nohalo")

    ! Set the computational domain sizes using the ice model's indexing convention.
    isc = sHI%isc ; iec = sHI%iec ; jsc = sHI%jsc ; jec = sHI%jec
    i_off = LBOUND(Ice%area,1) - sHI%isc ; j_off = LBOUND(Ice%area,2) - sHI%jsc
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      Ice%area(i2,j2) = US%L_to_m**2 * sG%areaT(i,j) * sG%mask2dT(i,j)
    enddo ; enddo

  endif ! slow_ice_PE

  !   Allocate the various structures and register fields for restarts connected
  ! to the fast ice processes.  This is interspersed between the slow ice
  ! registration calls and the actual reading of the restart files because there
  ! might be a single common restart file being used, and because the fast ice
  ! state might use some structures that point to their counterparts in the slow
  ! ice state.
  if (fast_ice_PE) then
    if (.not.associated(Ice%fCS)) allocate(Ice%fCS)
    if (.not.associated(Ice%fCS%IG)) allocate(Ice%fCS%IG)
    Ice%fCS%Time = Time

    if (single_IST .and. .not.redo_fast_update) then
      Ice%fCS%IST => Ice%sCS%IST
    else
      if (.not.associated(Ice%fCS%IST)) allocate(Ice%fCS%IST)
    endif

    Ice%fCS%US => US

    if (single_IST) then
      Ice%fCS%G => Ice%sCS%G
      fG => Ice%fCS%G
      fGD => Ice%fCS%G%Domain
      fHI = sHI
      Ice%fCS%FIA => Ice%sCS%FIA
      Ice%fCS%sOSS => Ice%sCS%sOSS
    else
      ! Set up the domains and lateral grids.
      Ice%fCS%IST%Cgrid_dyn = Cgrid_dyn
      if (.not.associated(Ice%fCS%G)) allocate(Ice%fCS%G)
      fG => Ice%fCS%G

      ! Set up the MOM_domain_type structures.
      call MOM_domains_init(Ice%fCS%G%domain, param_file, symmetric=.true., &
               domain_name="ice model fast", include_name="SIS2_memory.h", &
               static_memory=.false., NIHALO=0, NJHALO=0, param_suffix="_FAST")
      fGD => Ice%fCS%G%Domain

      call callTree_waypoint("domains initialized (ice_model_init)")
      call hor_index_init(fGD, fHI, param_file, &
                          local_indexing=.not.global_indexing)

      call create_dyn_horgrid(dG, fHI) !, bathymetry_at_vel=bathy_at_vel)
      call clone_MOM_domain(fGD, dG%Domain)

      ! Set the bathymetry, Coriolis parameter, open channel widths and masks.
      call SIS_initialize_fixed(dG, US, param_file, .false., dirs%output_directory, OBC_in)

      call set_hor_grid(Ice%fCS%G, param_file, global_indexing=global_indexing)
      call copy_dyngrid_to_SIS_horgrid(dG, Ice%fCS%G)
      call destroy_dyn_horgrid(dG)
      Ice%fCS%G%US => US
      Ice%OBC => OBC_in
    endif

    Ice%fCS%bounds_check = bounds_check
    Ice%fCS%debug = debug_fast
    Ice%fCS%Eulerian_tsurf = Eulerian_tsurf
    Ice%fCS%redo_fast_update = redo_fast_update

    ! Set up the ice-specific grid describing categories and ice layers.
    call set_ice_grid(Ice%fCS%IG, US, param_file, nCat_dflt, ocean_part_min_dflt=opm_dflt)
    if (slab_ice) Ice%fCS%IG%CatIce = 1 ! open water and ice ... but never in same place
    CatIce = Ice%fCS%IG%CatIce ; NkIce = Ice%fCS%IG%NkIce

    call initialize_ice_categories(Ice%fCS%IG, Rho_ice, US, param_file)

  ! Allocate and register fields for restarts.

    if (.not.slow_ice_PE) call set_domain(fGD%mpp_domain)
    ! if this is just a fast ice PE then read restart filename
    if (.not.associated(Ice%Ice_restart)) &
        call SIS_restart_init(Ice%Ice_restart, restart_file, sGD, param_file)

    if (split_restart_files) then
      if (.not.associated(Ice%Ice_fast_restart)) &
        call SIS_restart_init(Ice%Ice_fast_restart, fast_rest_file, fGD, param_file)
    else
      Ice%Ice_fast_restart => Ice%Ice_restart
    endif

  ! These allocation routines are called on all PEs; whether or not the variables
  ! they allocate are registered for inclusion in restart files is determined by
  ! whether the Ice%Ice...restart types are associated.
    call ice_type_fast_reg_restarts(fGD%mpp_domain, CatIce, &
                      param_file, Ice, Ice%Ice_fast_restart)
    if (split_restart_files) &
      call register_unit_conversion_restarts(Ice%fCS%US, Ice%Ice_fast_restart)

    if (redo_fast_update .or. .not.single_IST) then
      call alloc_IST_arrays(fHI, Ice%fCS%IG, US, Ice%fCS%IST, &
                            omit_velocities=.true., omit_tsurf=Eulerian_tsurf)
    endif
    if (.not.single_IST) then
      call alloc_fast_ice_avg(Ice%fCS%FIA, fHI, Ice%fCS%IG, interp_fluxes, gas_fluxes)

      call alloc_simple_OSS(Ice%fCS%sOSS, fHI, gas_fields_ocn)
    endif
    call alloc_total_sfc_flux(Ice%fCS%TSF, fHI, gas_fluxes)
    Ice%fCS%FIA%atmos_winds = atmos_winds

    call ice_rad_register_restarts(fHI, Ice%fCS%IG, US, param_file, Ice%fCS%Rad, Ice%Ice_fast_restart)
    Ice%fCS%Rad%do_sun_angle_for_alb = do_sun_angle_for_alb
    Ice%fCS%Rad%add_diurnal_sw = add_diurnal_sw

    if (Concurrent_ice) then
      call register_fast_to_slow_restarts(Ice%fCS%FIA, Ice%fCS%Rad, Ice%fCS%TSF, &
                       fGD%mpp_domain, US, Ice%Ice_fast_restart, fast_rest_file)
    endif

    allocate(Ice%fCS%diag)
    call SIS_diag_mediator_init(fG, Ice%fCS%IG, param_file, Ice%fCS%diag, component="SIS_fast", &
                                doc_file_dir = dirs%output_directory)
    call set_SIS_axes_info(fG, Ice%fCS%IG, param_file, Ice%fCS%diag, axes_set_name="ice_fast")

    if (redo_fast_update .or. .not.single_IST) then
      call ice_thermo_init(param_file, Ice%fCS%IST%ITV, US, init_EOS=nudge_sea_ice)
    endif

    if (.not.single_IST) then
      ! Set a few final things to complete the setup of the grid.
      fG%g_Earth = g_Earth
      call set_first_direction(fG, first_direction)
      call clone_MOM_domain(fGD, fG%domain_aux, symmetric=.false., &
                            domain_name="ice model aux")

    endif

    ! Copy the ice model's domain into one with no halos that can be shared
    ! publicly for use by the exchange grid.
    call clone_MOM_domain(Ice%fCS%G%Domain, Ice%domain, halo_size=0, &
                          symmetric=.false., domain_name="ice_nohalo")
  endif


  ! Read the restart file, if it exists and reading it is indicated by the value of dirs%input_filename,
  ! and initialize the ice arrays to default values using other methods if it does not.
  if (slow_ice_PE) then
    ! Set some pointers for convenience.
    sIST => Ice%sCS%IST ; sIG => Ice%sCS%IG ; sG => Ice%sCS%G

    new_sim = is_new_run(Ice%Ice_restart)
    if (.not.new_sim) then
      call callTree_enter("ice_model_init():restore_from_restart_files "//trim(restart_file))
      ! Set a value of IG%H_to_kg_m2 that will permit its absence from the restart file to be
      ! detected, as a way to detect an archaic restart file format.
      sIG%H_to_kg_m2 = -1.0
      is_restart = .true.
      recategorize_ice = .false. ! Assume that the ice is already in the right thickness categories.

      call restore_SIS_state(Ice%Ice_restart, dirs%restart_input_dir, dirs%input_filename, sG)

      ! If the velocity and other fields have not been initialized, check for
      ! the fields that would have been read if symmetric were toggled.
      call ice_state_read_alt_restarts(sIST, sG, sIG, US, Ice%Ice_restart, dirs%restart_input_dir)
      if (.not.specified_ice) &
        call SIS_dyn_trans_read_alt_restarts(Ice%sCS%dyn_trans_CSp, sG, US, Ice%Ice_restart, &
                                             dirs%restart_input_dir)

      call rescale_ice_state_restart_fields(sIST, sG, US, sIG, Rho_ice, Rho_snow)
      sIG%H_to_kg_m2 = 1.0

      if ((.not.query_initialized(Ice%Ice_restart, 'enth_ice')) .or. &
          (.not.query_initialized(Ice%Ice_restart, 'enth_snow')) .or. &
          (.not.query_initialized(Ice%Ice_restart, 'sal_ice'))) then
        ! Approximately initialize state fields that are not present
        ! in SIS1 restart files.  This is obsolete and can probably be eliminated.
        call read_archaic_thermo_restarts(Ice, sIST, sG, sIG, US, param_file, dirs, restart_file)
      endif

      if (Ice%sCS%pass_stress_mag .and. .not.query_initialized(Ice%Ice_restart, 'stress_mag')) then
        ! Determine the magnitude of the stresses from the (non-symmetric-memory) stresses
        ! in the Ice type, which will have been read from the restart files.
        allocate(str_x(sG%isd:sG%ied, sG%jsd:sG%jed), source=0.0)
        allocate(str_y(sG%isd:sG%ied, sG%jsd:sG%jed), source=0.0)
        allocate(stress_mag(sG%isd:sG%ied, sG%jsd:sG%jed), source=0.0)

        i_off = LBOUND(Ice%stress_mag,1) - sG%isc ; j_off = LBOUND(Ice%stress_mag,2) - sG%jsc
        do j=sG%jsc,sG%jec ; do i=sG%isc,sG%iec ; i2 = i+i_off ; j2 = j+j_off ! Correct for indexing differences.
          str_x(i,j) = Ice%flux_u(i2,j2) ; str_y(i,j) = Ice%flux_v(i2,j2)
        enddo ; enddo
        ! This serves to fill in the symmetric-edge stress points
        if ((Ice%flux_uv_stagger == BGRID_NE) .or. (Ice%flux_uv_stagger == CGRID_NE)) &
          call pass_vector(str_x, str_y, sG%Domain_aux, stagger=Ice%flux_uv_stagger) ! Breaks gnu:, halo=1)

        call stresses_to_stress_mag(sG, str_x, str_y, Ice%flux_uv_stagger, stress_mag)

        do j=sG%jsc,sG%jec ; do i=sG%isc,sG%iec ; i2 = i+i_off ; j2 = j+j_off ! Correct for indexing differences.
          Ice%stress_mag(i2,j2) = stress_mag(i,j)
        enddo ; enddo

        deallocate(str_x, str_y, stress_mag)
      endif

      if (fast_ice_PE .and. .not.split_restart_files) then
        init_coszen = .not.query_initialized(Ice%Ice_fast_restart, 'coszen')
        init_Tskin  = .not.query_initialized(Ice%Ice_fast_restart, 'T_skin')
        init_rough  = .not.(query_initialized(Ice%Ice_fast_restart, 'rough_mom') .and. &
                            query_initialized(Ice%Ice_fast_restart, 'rough_heat') .and. &
                            query_initialized(Ice%Ice_fast_restart, 'rough_moist'))
      endif

      call callTree_leave("ice_model_init():restore_from_restart_files")
    endif ! End of (.not.new_sim)

    ! If there is not a restart file, initialize the ice another way, perhaps with no ice.
    ! If there is a restart file, the following two calls are just here to read and log parameters.
    call ice_state_mass_init(sIST, Ice, sG, sIG, US, param_file, Ice%sCS%Time, just_read_params=is_restart)

    call ice_state_thermo_init(sIST, Ice, sG, sIG, US, param_file, Ice%sCS%Time, &
                               just_read_params=is_restart)

    if (.not.is_restart) then
      ! Record the need to transfer ice to the correct thickness category.
      recategorize_ice = .true.
      init_coszen = .true. ; init_Tskin = .true. ; init_rough = .true.
    endif

    ! The restart files have now been read or the variables that would have been in the restart
    ! files have been initialized, although some corrections and halo updates still need to be done.
    ! Now call the initialization routines for any dependent sub-modules.

    call ice_diagnostics_init(Ice%sCS%IOF, Ice%sCS%OSS, Ice%sCS%FIA, sG, US, sIG, &
                              Ice%sCS%diag, Ice%sCS%Time, Cgrid=sIST%Cgrid_dyn)
    Ice%axes(1:3) = Ice%sCS%diag%axesTc0%handles(1:3)

    Ice%sCS%Time_step_slow = Time_step_slow

    call SIS_slow_thermo_init(Ice%sCS%Time, sG, US, sIG, param_file, Ice%sCS%diag, &
                              Ice%sCS%slow_thermo_CSp, Ice%sCS%SIS_tracer_flow_CSp)

    if (specified_ice) then
      recategorize_ice = .false.
      call specified_ice_init(Ice%sCS%Time, sG, sIG, param_file, Ice%sCS%diag, &
                              Ice%sCS%specified_ice_CSp, dirs%output_directory, Time_Init)
      call SIS_slow_thermo_set_ptrs(Ice%sCS%slow_thermo_CSp, &
                   sum_out_CSp=specified_ice_sum_output_CS(Ice%sCS%specified_ice_CSp))

      ! When SPECIFIED_ICE=True, the variable Ice%sCS%OSS%SST_C is used for the skin temperature
      ! and needs to be updated for each run segment, regardless of whether a restart file is used.
      call get_sea_surface(Ice%sCS%Time, sG%HI, SST=Ice%sCS%OSS%SST_C, ice_domain=Ice%slow_domain_NH, US=US)
      if (.not.is_restart) then
        ! Perhaps ice_conc and h_ice_input should also be read with the get_sea_surface limits.
        allocate(h_ice_input(sG%isd:sG%ied, sG%jsd:sG%jed), source=0.0)
        call get_sea_surface(Ice%sCS%Time, sG%HI, SST=Ice%sCS%OSS%SST_C, ice_conc=sIST%part_size(:,:,1), &
                             ice_thick=h_ice_input, ice_domain=Ice%slow_domain_NH, US=US)
        do j=jsc,jec ; do i=isc,iec
          sIST%part_size(i,j,0) = 1.0 - sIST%part_size(i,j,1)
          sIST%mH_ice(i,j,1) = h_ice_input(i,j)*US%m_to_Z * Rho_ice
        enddo ; enddo
        deallocate(h_ice_input)
      endif
    else
      call SIS_dyn_trans_init(Ice%sCS%Time, sG, US, sIG, param_file, Ice%sCS%diag, &
                              Ice%sCS%dyn_trans_CSp, dirs%output_directory, Time_Init, &
                              slab_ice=slab_ice)
      call SIS_slow_thermo_set_ptrs(Ice%sCS%slow_thermo_CSp, &
                   transport_CSp=SIS_dyn_trans_transport_CS(Ice%sCS%dyn_trans_CSp), &
                   sum_out_CSp=SIS_dyn_trans_sum_output_CS(Ice%sCS%dyn_trans_CSp))
    endif


    ! Apply corrections to the ice state, like readjusting ice categories and fixing values on land.
    ! These corrections occur here so that they can use adjust_ice_categories.

    if (do_mask_restart) then
      ! Deal with any ice masses, thicknesses and other properties over land.
      if (allocated(sIST%t_surf)) then ; do j=jsc,jec ; do i=isc,iec
        if (sG%mask2dT(i,j) < 0.5) sIST%t_surf(i,j,:) = sIST%T_0degC
      enddo ; enddo ; endif
      do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
        sIST%sal_ice(i,j,k,n) = sIST%sal_ice(i,j,k,n) * sG%mask2dT(i,j)
        sIST%enth_ice(i,j,k,n) = sIST%enth_ice(i,j,k,n) * sG%mask2dT(i,j)
      enddo ; enddo ; enddo ; enddo
      do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
        sIST%mH_snow(i,j,k) = sIST%mH_snow(i,j,k) * sG%mask2dT(i,j)
        sIST%enth_snow(i,j,k,1) = sIST%enth_snow(i,j,k,1) * sG%mask2dT(i,j)
        sIST%mH_ice(i,j,k) = sIST%mH_ice(i,j,k) * sG%mask2dT(i,j)
        sIST%mH_pond(i,j,k) = sIST%mH_pond(i,j,k) * sG%mask2dT(i,j)
        sIST%part_size(i,j,k) = sIST%part_size(i,j,k) * sG%mask2dT(i,j)
      enddo ; enddo ; enddo
      ! Since we masked out the part_size on land we should set
      ! part_size(i,j,0) = 1. on land to satisfy the summation check
      do j=jsc,jec ; do i=isc,iec
        if (sG%mask2dT(i,j) < 0.5) sIST%part_size(i,j,0) = 1.
      enddo ; enddo
    endif

    if (recategorize_ice) then
      ! Transfer ice to the correct thickness categories.  This does not change the thicknesses or
      ! other properties where the ice is already in the correct categories.
      call adjust_ice_categories(sIST%mH_ice, sIST%mH_snow, sIST%mH_pond, sIST%part_size, &
                            sIST%TrReg, sG, sIG, SIS_dyn_trans_transport_CS(Ice%sCS%dyn_trans_CSp))
    endif

    if (sIG%ocean_part_min > 0.0) then ; do j=jsc,jec ; do i=isc,iec
      sIST%part_size(i,j,0) = max(sIST%part_size(i,j,0), sIG%ocean_part_min)
    enddo ; enddo ; endif

    !--- update the halo values for the physical ice state.
    call pass_var(sIST%part_size, sGD, complete=.true.)
    call pass_var(sIST%mH_ice, sGD, complete=.false.)
    call pass_var(sIST%mH_snow, sGD, complete=.false.)
    call pass_var(sIST%mH_pond, sGD, complete=.false.)
    do l=1,NkIce
      call pass_var(sIST%enth_ice(:,:,:,l), sGD, complete=.false.)
    enddo
    call pass_var(sIST%enth_snow(:,:,:,1), sGD, complete=.true.)
    if (Cgrid_dyn) then
      call pass_vector(sIST%u_ice_C, sIST%v_ice_C, sGD, stagger=CGRID_NE)
    else
      call pass_vector(sIST%u_ice_B, sIST%v_ice_B, sGD, stagger=BGRID_NE)
    endif

    ! The slow physical ice properties do not change after this point.

    if (Ice%sCS%redo_fast_update) then
      call SIS_fast_thermo_init(Ice%sCS%Time, sG, sIG, param_file, Ice%sCS%diag, &
                                Ice%sCS%fast_thermo_CSp)
      call SIS_optics_init(param_file, US, Ice%sCS%optics_CSp, slab_optics=slab_ice)
    endif

  !   Initialize any tracers that will be handled via tracer flow control.
    call SIS_tracer_flow_control_init(Ice%sCS%Time, sG, sIG, param_file, &
                                      Ice%sCS%SIS_tracer_flow_CSp, is_restart)

  ! Initialize icebergs
    if (Ice%sCS%do_icebergs) then
      call get_param(param_file, mdl, "ICEBERG_WINDSTRESS_BUG", Ice%sCS%berg_windstress_bug, &
                 "If true, use older code that applied an old ice-ocean "//&
                 "stress to the icebergs in place of the current air-ocean "//&
                 "stress.  This option is here for backward compatibility, "//&
                 "but should be avoided.", default=.false.)

      isc = sG%isc ; iec = sG%iec ; jsc = sG%jsc ; jec = sG%jec

      if (ASSOCIATED(sGD%maskmap)) then
        call icebergs_init(Ice%icebergs, sGD%niglobal, sGD%njglobal, &
                sGD%layout, sGD%io_layout, Ice%axes(1:2), &
                sGD%X_flags, sGD%Y_flags, time_type_to_real(Time_step_slow), &
                Time, sG%geoLonBu(isc:iec,jsc:jec), sG%geoLatBu(isc:iec,jsc:jec), &
                sG%mask2dT(isc-1:iec+1,jsc-1:jec+1), &
                US%L_to_m*sG%dxCv(isc-1:iec+1,jsc-1:jec+1), US%L_to_m*sG%dyCu(isc-1:iec+1,jsc-1:jec+1), &
                Ice%area,  sG%cos_rot(isc-1:iec+1,jsc-1:jec+1), &
                sG%sin_rot(isc-1:iec+1,jsc-1:jec+1), maskmap=sGD%maskmap )
      else
        call icebergs_init(Ice%icebergs, sGD%niglobal, sGD%njglobal, &
                 sGD%layout, sGD%io_layout, Ice%axes(1:2), &
                 sGD%X_flags, sGD%Y_flags, time_type_to_real(Time_step_slow), &
                 Time, sG%geoLonBu(isc:iec,jsc:jec), sG%geoLatBu(isc:iec,jsc:jec), &
                 sG%mask2dT(isc-1:iec+1,jsc-1:jec+1), &
                 US%L_to_m*sG%dxCv(isc-1:iec+1,jsc-1:jec+1), US%L_to_m*sG%dyCu(isc-1:iec+1,jsc-1:jec+1), &
                 Ice%area, sG%cos_rot(isc-1:iec+1,jsc-1:jec+1), &
                 sG%sin_rot(isc-1:iec+1,jsc-1:jec+1) )
      endif
    endif

    ! Do any error checking here.
    if (Ice%sCS%debug) call ice_grid_chksum(sG, US, haloshift=1)

    if (specified_ice) then
      call write_ice_statistics(sIST, Ice%sCS%Time, 0, sG, US, sIG, &
               specified_ice_sum_output_CS(Ice%sCS%specified_ice_CSp))
    else
      call write_ice_statistics(sIST, Ice%sCS%Time, 0, sG, US, sIG, &
               SIS_dyn_trans_sum_output_CS(Ice%sCS%dyn_trans_CSp))
    endif
  endif  ! slow_ice_PE


  if (fast_ice_PE) then
    ! Read the fast_restart file and initialize the subsidiary modules of the
    ! fast ice processes.

    ! Set some pointers for convenience.
    fG => Ice%fCS%G ; fGD => Ice%fCS%G%Domain

    if ((.not.slow_ice_PE) .or. split_restart_files) then
      ! Read the fast restart file, if it exists and this is indicated by the value of dirs%input_filename.
      new_sim = determine_is_new_run(dirs%input_filename, dirs%restart_input_dir, fG, Ice%Ice_fast_restart)
      if (.not.new_sim) then
        call restore_SIS_state(Ice%Ice_fast_restart, dirs%restart_input_dir, dirs%input_filename, fG)
        init_coszen = .not.query_initialized(Ice%Ice_fast_restart, 'coszen')
        init_Tskin  = .not.query_initialized(Ice%Ice_fast_restart, 'T_skin')
        init_rough  = .not.(query_initialized(Ice%Ice_fast_restart, 'rough_mom') .and. &
                            query_initialized(Ice%Ice_fast_restart, 'rough_heat') .and. &
                            query_initialized(Ice%Ice_fast_restart, 'rough_moist'))
      else
        init_coszen = .true. ; init_Tskin = .true. ; init_rough = .true.
      endif
    endif

    if (Concurrent_ice) then
      call rescale_fast_to_slow_restart_fields(Ice%fCS%FIA, Ice%fCS%Rad, Ice%fCS%TSF, &
                                               Ice%fCS%G, US, Ice%fCS%IG)
    endif


!  if (Ice%fCS%Rad%add_diurnal_sw .or. Ice%fCS%Rad%do_sun_angle_for_alb) then
!    call set_domain(fGD%mpp_domain)
    call astronomy_init()
!    call nullify_domain()
!  endif

    if (init_coszen) then
      if (coszen_IC >= 0.0) then
        Ice%fCS%Rad%coszen_nextrad(:,:) = coszen_IC
      else
        rad = acos(-1.)/180.
        allocate(dummy(fG%isd:fG%ied,fG%jsd:fG%jed))
        call diurnal_solar(fG%geoLatT(:,:)*rad, fG%geoLonT(:,:)*rad, &
                           Ice%fCS%Time, cosz=Ice%fCS%Rad%coszen_nextrad, fracday=dummy, &
                           rrsun=rrsun, dt_time=dT_rad)
        deallocate(dummy)
      endif
    endif
    if (init_Tskin) then
      Ice%fCS%Rad%t_skin(:,:,:) = 0.0
    elseif (do_mask_restart) then
      do k=1,CatIce
        Ice%fCS%Rad%t_skin(:,:,k) = Ice%fCS%Rad%t_skin(:,:,k) * fG%mask2dT(:,:)
      enddo
    endif
    if (init_rough) then
      Ice%rough_mom(:,:,:)   = mom_rough_ice
      Ice%rough_heat(:,:,:)  = heat_rough_ice
      Ice%rough_moist(:,:,:) = heat_rough_ice
    endif

    call ice_diags_fast_init(Ice%fCS%Rad, fG, Ice%fCS%IG, Ice%fCS%diag, &
                             Ice%fCS%Time, component="ice_model_fast")

    call SIS_fast_thermo_init(Ice%fCS%Time, fG, Ice%fCS%IG, param_file, Ice%fCS%diag, &
                              Ice%fCS%fast_thermo_CSp)
    call SIS_optics_init(param_file, US, Ice%fCS%optics_CSp, slab_optics=slab_ice)

    Ice%fCS%Time_step_fast = Time_step_fast
    Ice%fCS%Time_step_slow = Time_step_slow

    isc = fHI%isc ; iec = fHI%iec ; jsc = fHI%jsc ; jec = fHI%jec
    i_off = LBOUND(Ice%ocean_pt,1) - fHI%isc ; j_off = LBOUND(Ice%ocean_pt,2) - fHI%jsc
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      Ice%ocean_pt(i2,j2) = ( fG%mask2dT(i,j) > 0.0 )
    enddo ; enddo
    if (.not.slow_ice_PE) then
      Ice%axes(1:3) = Ice%fCS%diag%axesTc0%handles(1:3)
    endif
  endif ! fast_ice_PE
   
  call fix_restart_unit_scaling(US, unscaled=.true.)

  !nullify_domain perhaps could be called somewhere closer to set_domain
  !but it should be called after restore_SIS_state() otherwise it causes a restart mismatch
  call nullify_domain()

  ! Close the parameter file, supplying information for logging the default values.
  if (slow_ice_PE) then
    call close_param_file(param_file, component='SIS')
  else
    call close_param_file(param_file, component='SIS_fast')
  endif

  ! Ice%xtype can be REDIST or DIRECT, depending on the relationship between
  ! the fast and slow ice PEs.  REDIST should always work but may be slower.
  if (fast_ice_PE .neqv. slow_ice_PE) then
    Ice%xtype = REDIST
  elseif (single_IST .or. ((fGD%layout(1) == sGD%layout(1)) .and. &
                           (fGD%layout(2) == sGD%layout(2))) ) then
    Ice%xtype = DIRECT
  else
    Ice%xtype = REDIST
  endif
  
  if (fast_ice_PE .eqv. slow_ice_PE) then
    call exchange_fast_to_slow_ice(Ice)
  endif 

  if (Ice%shared_slow_fast_PEs) then
    iceClock = cpu_clock_id( 'Ice', grain=CLOCK_COMPONENT )
    ice_clock_fast = cpu_clock_id('Ice Fast', grain=CLOCK_SUBCOMPONENT )
    ice_clock_slow = cpu_clock_id('Ice Slow', grain=CLOCK_SUBCOMPONENT )
  else
    iceClock = 0 ! The comprehensive ice clock can not be used with separate fast and slow ice PEs.
    if (fast_ice_PE) then
      ice_clock_fast = cpu_clock_id('Ice Fast', grain=CLOCK_COMPONENT )
    elseif (slow_ice_PE) then
      ice_clock_slow = cpu_clock_id('Ice Slow', grain=CLOCK_COMPONENT )
    endif
  endif

  call callTree_leave("ice_model_init()")

end subroutine ice_model_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> share_ice_domains exchanges domain information between the fast and slow ice PEs
subroutine share_ice_domains(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type.

  ! This code has to be called for all of the ice processors using the
  ! union of the fast and slow ice PE_lists.

  if (associated(Ice%sCS)) then
    Ice%slow_domain => Ice%sCS%G%domain%mpp_domain
  else
    allocate(Ice%slow_domain)
  endif
  if (associated(Ice%fCS)) then
    Ice%fast_domain => Ice%fCS%G%domain%mpp_domain
  else
    allocate(Ice%fast_domain)
  endif
  call broadcast_domain(Ice%Domain)
  call broadcast_domain(Ice%slow_domain_NH)
  call broadcast_domain(Ice%slow_domain)
  call broadcast_domain(Ice%fast_domain)

  if (Ice%shared_slow_fast_PEs) then
    ice_clock_exchange = cpu_clock_id('Ice Fast/Slow Exchange', grain=CLOCK_SUBCOMPONENT )
  else
    ice_clock_exchange = cpu_clock_id('Ice Fast/Slow Exchange', grain=CLOCK_COMPONENT )
  endif

end subroutine share_ice_domains


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> update_ice_atm_deposition_flux updates the values of any fluxes that are
!! labeled as having the type 'air_sea_deposition'.  With the FMS coupler,
!! these fluxes were not available when update_ice_model_fast was called.
subroutine update_ice_atm_deposition_flux( Atmos_boundary, Ice )
  type(ice_data_type),           intent(inout) :: Ice !< The publicly visible ice data type.
  type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary !< A type containing atmospheric boundary
                                                      !! forcing fields that are used to drive the ice

  call accumulate_deposition_fluxes(Atmos_boundary, Ice%fCS%FIA, Ice%fCS%G, Ice%fCS%IG)

end subroutine update_ice_atm_deposition_flux


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_model_end writes the restart file and deallocates memory
subroutine ice_model_end(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type.

  logical :: fast_ice_PE       ! If true, fast ice processes are handled on this PE.
  logical :: slow_ice_PE       ! If true, slow ice processes are handled on this PE.

  call ice_model_restart(Ice=Ice)

  !--- release memory ------------------------------------------------

  fast_ice_PE = associated(Ice%fCS)
  slow_ice_PE = associated(Ice%sCS)

  if (fast_ice_PE) then
    call SIS_fast_thermo_end(Ice%fCS%fast_thermo_CSp)

    call SIS_optics_end(Ice%fCS%optics_CSp)

    if (Ice%fCS%Rad%add_diurnal_sw .or. Ice%fCS%Rad%do_sun_angle_for_alb) &
      call astronomy_end

    call dealloc_ice_rad(Ice%fCS%Rad)
    call dealloc_total_sfc_flux(Ice%fCS%TSF)
    call ice_grid_end(Ice%fCS%IG)

    if (.not.associated(Ice%sCS)) then
      call dealloc_IST_arrays(Ice%fCS%IST)
      deallocate(Ice%fCS%IST)

      call dealloc_fast_ice_avg(Ice%fCS%FIA)
      call dealloc_simple_OSS(Ice%fCS%sOSS)
      call SIS_hor_grid_end(Ice%fCS%G)
    else
      if (.not.associated(Ice%fCS%IST,Ice%sCS%IST)) then
        call dealloc_IST_arrays(Ice%fCS%IST)
        deallocate(Ice%fCS%IST)
      endif
      if (.not.associated(Ice%fCS%FIA,Ice%sCS%FIA)) &
        call dealloc_fast_ice_avg(Ice%fCS%FIA)
      if (.not.associated(Ice%fCS%sOSS,Ice%sCS%sOSS)) &
        call dealloc_simple_OSS(Ice%fCS%sOSS)
      if (.not.associated(Ice%fCS%G,Ice%sCS%G)) &
        call SIS_hor_grid_end(Ice%fCS%G)
    endif

    if (associated(Ice%Ice_fast_restart) .and. &
        (.not.associated(Ice%Ice_fast_restart, Ice%Ice_restart))) &
      deallocate(Ice%Ice_fast_restart)

  endif

  if (slow_ice_PE) then

    if (associated(Ice%sCS%dyn_trans_CSp)) &
      call SIS_dyn_trans_end(Ice%sCS%dyn_trans_CSp)

    if (associated(Ice%sCS%specified_ice_CSp)) &
      call specified_ice_end(Ice%sCS%specified_ice_CSp)

    call SIS_slow_thermo_end(Ice%sCS%slow_thermo_CSp)

    call ice_thermo_end(Ice%sCS%IST%ITV)

    ! End icebergs
    if (Ice%sCS%do_icebergs) call icebergs_end(Ice%icebergs)

    call SIS_tracer_flow_control_end(Ice%sCS%SIS_tracer_flow_CSp)

    call dealloc_ice_ocean_flux(Ice%sCS%IOF)

    if (Ice%sCS%redo_fast_update) then
      call SIS_fast_thermo_end(Ice%sCS%fast_thermo_CSp)
      call SIS_optics_end(Ice%sCS%optics_CSp)

      call dealloc_total_sfc_flux(Ice%sCS%TSF)
      call dealloc_total_sfc_flux(Ice%sCS%XSF)
      call dealloc_ice_rad(Ice%sCS%Rad)
    endif

    call dealloc_ocean_sfc_state(Ice%sCS%OSS)

    call dealloc_simple_OSS(Ice%sCS%sOSS)

    call ice_grid_end(Ice%sCS%IG)

    call dealloc_IST_arrays(Ice%sCS%IST)
    deallocate(Ice%sCS%IST)

    call dealloc_fast_ice_avg(Ice%sCS%FIA)

    call SIS_hor_grid_end(Ice%sCS%G)
  endif

  call dealloc_Ice_arrays(Ice)

  if (associated(Ice%Ice_restart)) deallocate(Ice%Ice_restart)

  if (slow_ice_PE) then
    call SIS_diag_mediator_end(Ice%sCS%Time, Ice%sCS%diag)
  else
    call SIS_diag_mediator_end(Ice%fCS%Time, Ice%fCS%diag)
  endif

  if (associated(Ice%fCS)) deallocate(Ice%fCS)
  if (associated(Ice%sCS)) deallocate(Ice%sCS)

end subroutine ice_model_end

end module ice_model_mod
