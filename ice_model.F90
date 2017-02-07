!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of SIS2.                                        *
!*                                                                     *
!* SIS2 is free software; you can redistribute it and/or modify it and *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* SIS2 is distributed in the hope that it will be useful, but WITHOUT *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

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
module ice_model_mod

use SIS_debugging,     only : chksum, uchksum, vchksum, Bchksum, SIS_debugging_init
use SIS_diag_mediator, only : set_SIS_axes_info, SIS_diag_mediator_init, SIS_diag_mediator_end
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
! use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
! use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_get_input, only : Get_SIS_input, directories
use SIS_sum_output, only : SIS_sum_output_init,  write_ice_statistics
use SIS_transcribe_grid, only : copy_dyngrid_to_SIS_horgrid, copy_SIS_horgrid_to_dyngrid

use MOM_domains,       only : MOM_domain_type
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : fill_symmetric_edges, MOM_domains_init, clone_MOM_domain
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, log_param, log_version, read_param, param_file_type
use MOM_file_parser, only : open_param_file, close_param_file
use MOM_hor_index, only : hor_index_type, hor_index_init
use MOM_obsolete_params, only : obsolete_logical
use MOM_string_functions, only : uppercase
use MOM_time_manager, only : time_type, time_type_to_real, real_to_time_type
use MOM_time_manager, only : set_date, set_time, operator(+), operator(-)
use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)

use fms_mod, only : file_exist, clock_flag_default
use fms_io_mod, only : set_domain, nullify_domain, restore_state, query_initialized
use fms_io_mod, only : restore_state, query_initialized
use fms_io_mod, only : register_restart_field, restart_file_type
use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
use mpp_domains_mod, only : mpp_broadcast_domain

use astronomy_mod, only : astronomy_init, astronomy_end
use astronomy_mod, only : universal_time, orbital_time, diurnal_solar, daily_mean_solar
use coupler_types_mod, only : coupler_3d_bc_type
use ocean_albedo_mod, only : compute_ocean_albedo            ! ice sets ocean surface
use ocean_rough_mod,  only : compute_ocean_roughness         ! properties over water

use ice_type_mod, only : ice_data_type, dealloc_ice_arrays
use ice_type_mod, only : ice_type_slow_reg_restarts, ice_type_fast_reg_restarts
use ice_type_mod, only : Ice_public_type_chksum, Ice_public_type_bounds_check
use ice_type_mod, only : ice_model_restart, ice_stock_pe, ice_data_type_chksum
use ice_boundary_types, only : ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
use ice_boundary_types, only : ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
use ice_boundary_types, only : lnd_ice_bnd_type_chksum
use SIS_ctrl_types, only : SIS_slow_CS, SIS_fast_CS
use SIS_ctrl_types, only : ice_diagnostics_init, ice_diags_fast_init
use SIS_types, only : ice_ocean_flux_type, alloc_ice_ocean_flux, dealloc_ice_ocean_flux
use SIS_types, only : ocean_sfc_state_type, alloc_ocean_sfc_state, dealloc_ocean_sfc_state
use SIS_types, only : fast_ice_avg_type, alloc_fast_ice_avg, dealloc_fast_ice_avg
use SIS_types, only : total_sfc_flux_type, alloc_total_sfc_flux, dealloc_total_sfc_flux
use SIS_types, only : ice_rad_type, ice_rad_register_restarts, dealloc_ice_rad
use SIS_types, only : simple_OSS_type, alloc_simple_OSS, dealloc_simple_OSS
use SIS_types, only : ice_state_type, alloc_IST_arrays, dealloc_IST_arrays
use SIS_types, only : IST_chksum, IST_bounds_check, ice_state_register_restarts
use SIS_types, only : copy_IST_to_IST, copy_FIA_to_FIA, copy_sOSS_to_sOSS
use SIS_types, only : redistribute_IST_to_IST, redistribute_FIA_to_FIA
use SIS_types, only : redistribute_sOSS_to_sOSS, FIA_chksum, IOF_chksum, translate_OSS_to_sOSS
use SIS_utils, only : post_avg, ice_grid_chksum
use SIS_hor_grid, only : SIS_hor_grid_type, set_hor_grid, SIS_hor_grid_end, set_first_direction
use SIS_fixed_initialization, only : SIS_initialize_fixed

use ice_grid, only : set_ice_grid, ice_grid_end, ice_grid_type
use ice_spec_mod, only : get_sea_surface

use SIS_tracer_registry, only : register_SIS_tracer, register_SIS_tracer_pair
use SIS_tracer_flow_control, only : SIS_call_tracer_register, SIS_tracer_flow_control_init
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_end

use ice_thm_mod,     only : slab_ice_optics
use SIS_dyn_trans,   only : SIS_dynamics_trans, update_icebergs
use SIS_dyn_trans,   only : SIS_dyn_trans_register_restarts, SIS_dyn_trans_init, SIS_dyn_trans_end
use SIS_dyn_trans,   only : SIS_dyn_trans_transport_CS, SIS_dyn_trans_sum_output_CS
use SIS_slow_thermo, only : slow_thermodynamics, SIS_slow_thermo_init, SIS_slow_thermo_end
use SIS_slow_thermo, only : SIS_slow_thermo_set_ptrs
use SIS_fast_thermo, only : do_update_ice_model_fast, avg_top_quantities, total_top_quantities
use SIS_fast_thermo, only : infill_array, SIS_fast_thermo_init, SIS_fast_thermo_end
use SIS_optics,      only : ice_optics_SIS2, SIS_optics_init, SIS_optics_end

use SIS2_ice_thm,  only : ice_temp_SIS2, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,  only : ice_thermo_init, ice_thermo_end, get_SIS2_thermo_coefs
use SIS2_ice_thm,  only : enth_from_TS, Temp_from_En_S, T_freeze, ice_thermo_type
use ice_bergs,     only : icebergs, icebergs_run, icebergs_init, icebergs_end

implicit none ; private

public :: ice_data_type, ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
public :: ice_model_init, share_ice_domains, ice_model_end, ice_stock_pe
public :: update_ice_model_fast
public :: update_ice_model_slow_up, update_ice_model_slow_dn ! The old Verona interfaces.
public :: ice_model_restart  ! for intermediate restarts
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum, ice_data_type_chksum
public :: unpack_ocean_ice_boundary, exchange_slow_to_fast_ice, set_ice_surface_fields
public :: ice_model_fast_cleanup, unpack_land_ice_boundary
public :: exchange_fast_to_slow_ice, update_ice_model_slow

integer :: iceClock
integer :: ice_clock_slow, ice_clock_fast, ice_clock_exchange

integer, parameter :: REDIST=2, DIRECT=3

contains

!-----------------------------------------------------------------------
!> Update the sea-ice state due to slow processes, including dynamics,
!! freezing and melting, precipitation, and transport.
subroutine update_ice_model_slow_dn ( Atmos_boundary, Land_boundary, Ice )
  type(atmos_ice_boundary_type), &
    intent(in)    :: Atmos_boundary !< Atmos_boundary is not actually used, and
                                   !! is still here only for backward compatibilty with the
                                   !! interface to Verona and earlier couplers.
  type(land_ice_boundary_type), &
    intent(in)    :: Land_boundary !< A structure containing information about
                                   !! the fluxes from the land that is being shared with the
                                   !! sea-ice.  If this argument is not present, it is assumed
                                   !! that this information has already been exchanged.
  type(ice_data_type), &
    intent(inout) :: Ice           !< The publicly visible ice data type; this must always be
                                   !! present, but is optional because of an unfortunate
                                   !! order of arguments.

  if (.not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "The pointer to Ice%sCS must be associated in update_ice_model_slow_dn.")

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(ice_clock_slow)

  call ice_model_fast_cleanup(Ice)

  call unpack_land_ice_boundary(Ice, Land_boundary)

  call exchange_fast_to_slow_ice(Ice)

  call mpp_clock_end(ice_clock_slow) ; call mpp_clock_end(iceClock)

  call update_ice_model_slow(Ice)

end subroutine update_ice_model_slow_dn


!-----------------------------------------------------------------------
!> Update the sea-ice state due to slow processes, including dynamics,
!! freezing and melting, precipitation, and transport.
subroutine update_ice_model_slow(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type.

  real :: dt_slow  ! The time step over which to advance the model.
  integer :: i, j, i2, j2, i_off, j_off

  if (.not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "The pointer to Ice%sCS must be associated in update_ice_model_slow.")

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(ice_clock_slow)

  ! Advance the slow PE clock to give the end time of the slow timestep.  There
  ! is a separate clock inside the fCS that is advanced elsewhere.
  Ice%sCS%Time = Ice%sCS%Time + Ice%sCS%Time_step_slow
  if (.not.associated(Ice%fCS)) then
    Ice%Time = Ice%sCS%Time
  endif
  dt_slow = time_type_to_real(Ice%sCS%Time_step_slow)

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Start update_ice_model_slow", Ice, check_slow=.true.)
    call FIA_chksum("Start update_ice_model_slow", Ice%sCS%FIA, Ice%sCS%G)
    call IOF_chksum("Start update_ice_model_slow", Ice%sCS%IOF, Ice%sCS%G)
  endif

  ! Store some diagnostic fluxes...
!$OMP parallel do default(none) shared(Ice)
  do j=Ice%sCS%G%jsc,Ice%sCS%G%jec ; do i=Ice%sCS%G%isc,Ice%sCS%G%iec
    Ice%sCS%FIA%calving_preberg(i,j) = Ice%sCS%FIA%calving(i,j)
    Ice%sCS%FIA%calving_hflx_preberg(i,j) = Ice%sCS%FIA%calving_hflx(i,j)
  enddo ; enddo

  if (Ice%sCS%do_icebergs) then
    if (Ice%sCS%berg_windstress_bug) then
      ! This code is only required to reproduce an old bug.
      i_off = LBOUND(Ice%flux_t,1) - Ice%sCS%G%isc
      j_off = LBOUND(Ice%flux_t,2) - Ice%sCS%G%jsc
!$OMP parallel do default(none) shared(Ice,i_off,j_off) private(i2,j2)
      do j=Ice%sCS%G%jsc,Ice%sCS%G%jec ; do i=Ice%sCS%G%isc,Ice%sCS%G%iec
        i2 = i+i_off ; j2 = j+j_off
        Ice%sCS%IOF%flux_u_ocn(i,j) = Ice%flux_u(i2,j2)
        Ice%sCS%IOF%flux_v_ocn(i,j) = Ice%flux_v(i2,j2)
      enddo ; enddo
    endif

    call mpp_clock_end(ice_clock_slow) ; call mpp_clock_end(iceClock)
    call update_icebergs(Ice%sCS%IST, Ice%sCS%OSS, Ice%sCS%IOF, Ice%sCS%FIA, Ice%icebergs, &
                         dt_slow, Ice%sCS%G, Ice%sCS%IG, Ice%sCS%dyn_trans_CSp)
    call mpp_clock_begin(iceClock) ; call mpp_clock_begin(ice_clock_slow)

    if (Ice%sCS%debug) then
      call FIA_chksum("After update_icebergs", Ice%sCS%FIA, Ice%sCS%G)
    endif
  endif

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Before slow_thermodynamics", Ice, check_slow=.true.)
    call IOF_chksum("Before slow_thermodynamics", Ice%sCS%IOF, Ice%sCS%G)
  endif

  call slow_thermodynamics(Ice%sCS%IST, dt_slow, Ice%sCS%slow_thermo_CSp, &
                           Ice%sCS%OSS, Ice%sCS%FIA, Ice%sCS%IOF, Ice%sCS%G, Ice%sCS%IG)

  ! Do halo updates on the forcing fields, as necessary.  This must occur before
  ! the call to SIS_dynamics_trans, because update_icebergs does its own halo
  ! updates, and slow_thermodynamics only works on the computational domain.
  call pass_vector(Ice%sCS%FIA%WindStr_x, Ice%sCS%FIA%WindStr_y, &
                   Ice%sCS%G%Domain, stagger=AGRID, complete=.false.)
  call pass_vector(Ice%sCS%FIA%WindStr_ocn_x, Ice%sCS%FIA%WindStr_ocn_y, &
                   Ice%sCS%G%Domain, stagger=AGRID)
  call pass_var(Ice%sCS%FIA%ice_cover, Ice%sCS%G%Domain, complete=.false.)
  call pass_var(Ice%sCS%FIA%ice_free,  Ice%sCS%G%Domain, complete=.true.)
  call pass_var(Ice%sCS%IST%part_size, Ice%sCS%G%Domain)
  call pass_var(Ice%sCS%IST%mH_ice, Ice%sCS%G%Domain, complete=.false.)
  call pass_var(Ice%sCS%IST%mH_pond, Ice%sCS%G%Domain, complete=.false.)
  call pass_var(Ice%sCS%IST%mH_snow, Ice%sCS%G%Domain, complete=.true.)

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Before SIS_dynamics_trans", Ice, check_slow=.true.)
    call IOF_chksum("Before SIS_dynamics_trans", Ice%sCS%IOF, Ice%sCS%G)
  endif

  call SIS_dynamics_trans(Ice%sCS%IST, Ice%sCS%OSS, Ice%sCS%FIA, Ice%sCS%IOF, &
                          dt_slow, Ice%sCS%dyn_trans_CSp, Ice%icebergs, Ice%sCS%G, &
                          Ice%sCS%IG, Ice%sCS%SIS_tracer_flow_CSp)

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Before set_ocean_top_fluxes", Ice, check_slow=.true.)
    call IOF_chksum("Before set_ocean_top_fluxes", Ice%sCS%IOF, Ice%sCS%G)
    call IST_chksum("Before set_ocean_top_fluxes", Ice%sCS%IST, Ice%sCS%G, Ice%sCS%IG)
  endif
  ! Set up the thermodynamic fluxes in the externally visible structure Ice.
  call set_ocean_top_fluxes(Ice, Ice%sCS%IST, Ice%sCS%IOF, Ice%sCS%FIA, Ice%sCS%OSS, &
                            Ice%sCS%G, Ice%sCS%IG, Ice%sCS)

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("End update_ice_model_slow", Ice, check_slow=.true.)
  endif

  !### THIS NO LONGER WORKS ON SLOW ICE PES.
!  if (Ice%sCS%bounds_check) then
!    call Ice_public_type_bounds_check(Ice, Ice%sCS%G, "End update_ice_slow")
!  endif

  call mpp_clock_end(ice_clock_slow) ; call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_model_fast_cleanup performs the final steps in the fast ice update cycle
!! and prepares data to drive the slow ice updates.  This includes finding the
!! averaged fluxes and unpacking the land to ice forcing.
subroutine ice_model_fast_cleanup(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type.

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  if (.not.associated(Ice%fCS)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS must be associated in ice_model_fast_cleanup.")

  ! average fluxes from update_ice_model_fast
  call avg_top_quantities(Ice%fCS%FIA, Ice%fCS%Rad, Ice%fCS%IST, &
                          Ice%fCS%G, Ice%fCS%IG)

  call total_top_quantities(Ice%fCS%FIA, Ice%fCS%TSF, Ice%fCS%IST%part_size, &
                            Ice%fCS%G, Ice%fCS%IG)

  if (allocated(Ice%fCS%IST%t_surf)) &
    Ice%fCS%IST%t_surf(:,:,1:) = Ice%fCS%Rad%T_skin(:,:,:) + T_0degC
  call infill_array(Ice%fCS%IST, Ice%fCS%sOSS%T_fr_ocn, Ice%fCS%Rad%T_skin, &
                    Ice%fCS%G, Ice%fCS%IG)

end subroutine ice_model_fast_cleanup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> unpack_land_ice_bdry converts the information in a publicly visible
!! land_ice_boundary_type into an internally visible fast_ice_avg_type variable.
subroutine unpack_land_ice_boundary(Ice, LIB)
  type(ice_data_type),          intent(inout) :: Ice !< The publicly visible ice data type.
  type(land_ice_boundary_type), intent(in)    :: LIB !< The land ice boundary type that is being unpacked.

  type(fast_ice_avg_type), pointer :: FIA => NULL()
  type(SIS_hor_grid_type), pointer :: G => NULL()

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off

  if (.not.associated(Ice%fCS)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS must be associated in unpack_land_ice_boundary.")
  if (.not.associated(Ice%fCS%FIA)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS%FIA must be associated in unpack_land_ice_boundary.")
  if (.not.associated(Ice%fCS%G)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS%G must be associated in unpack_land_ice_boundary.")

  FIA => Ice%fCS%FIA ; G => Ice%fCS%G

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  ! Store liquid runoff and other fluxes from the land to the ice or ocean.
  i_off = LBOUND(LIB%runoff,1) - G%isc ; j_off = LBOUND(LIB%runoff,2) - G%jsc
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,FIA,LIB,i_off,j_off,G) &
!$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j) > 0.0) then
    i2 = i+i_off ; j2 = j+j_off
    FIA%runoff(i,j)  = LIB%runoff(i2,j2)
    FIA%calving(i,j) = LIB%calving(i2,j2)
    FIA%runoff_hflx(i,j)  = LIB%runoff_hflx(i2,j2)
    FIA%calving_hflx(i,j) = LIB%calving_hflx(i2,j2)
  else
    ! This is a land point from the perspective of the sea-ice.
    ! At some point it might make sense to check for non-zero fluxes, which
    ! might indicate regridding errors.  However, bad-data values are also
    ! non-zero and should not be flagged.
    FIA%runoff(i,j)  = 0.0 ; FIA%calving(i,j) = 0.0
    FIA%runoff_hflx(i,j)  = 0.0 ; FIA%calving_hflx(i,j) = 0.0
  endif ; enddo ; enddo

  if (Ice%fCS%debug) then
    call FIA_chksum("End of unpack_land_ice_boundary", FIA, G)
  endif

end subroutine unpack_land_ice_boundary

!> This subroutine copies information (mostly fluxes and the updated temperatures)
!! from the fast part of the sea-ice to the  slow part of the sea ice.
subroutine exchange_fast_to_slow_ice(Ice)
  type(ice_data_type), &
    intent(inout) :: Ice            !< The publicly visible ice data type whose fast
                                    !! part is to be exchanged with the slow part.
  type(fast_ice_avg_type), pointer :: FIA_null => NULL()
  type(ice_state_type),    pointer :: IST_null => NULL()

  if(Ice%xtype == DIRECT) then
    if (.not.associated(Ice%fCS) .or. .not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "With xtype=DIRECT, both the pointer to Ice%sCS and the pointer to Ice%fCS must be "//&
      "associated (although perhaps not with each other) in exchange_fast_to_slow_ice.")

    if (.not.associated(Ice%fCS%FIA, Ice%sCS%FIA)) then
      call copy_FIA_to_FIA(Ice%fCS%FIA, Ice%sCS%FIA, Ice%fCS%G%HI, Ice%sCS%G%HI, &
                           Ice%sCS%IG)
    endif

    if (.not.associated(Ice%fCS%IST, Ice%sCS%IST)) then
      call copy_IST_to_IST(Ice%fCS%IST, Ice%sCS%IST, Ice%fCS%G%HI, Ice%sCS%G%HI, &
           Ice%fCS%IG)
    endif
  elseif (Ice%xtype == REDIST) then
    if (.not.associated(Ice%fCS) .and. .not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "Either the pointer to Ice%sCS or the pointer to Ice%fCS must be "//&
      "associated in exchange_fast_to_slow_ice.")

    if (associated(Ice%fCS) .and. associated(Ice%sCS)) then
      if (.not.associated(Ice%fCS%FIA, Ice%sCS%FIA)) &
        call redistribute_FIA_to_FIA(Ice%fCS%FIA, Ice%sCS%FIA, Ice%fast_domain, &
                                     Ice%slow_domain, Ice%sCS%G, Ice%sCS%IG)

      if (.not.associated(Ice%fCS%IST, Ice%sCS%IST)) &
        call redistribute_IST_to_IST(Ice%fCS%IST, Ice%sCS%IST, Ice%fast_domain, &
                                     Ice%slow_domain)
    elseif (associated(Ice%fCS)) then
      call redistribute_FIA_to_FIA(Ice%fCS%FIA, FIA_null, Ice%fast_domain, &
                                   Ice%slow_domain)
      call redistribute_IST_to_IST(Ice%fCS%IST, IST_null, Ice%fast_domain, &
                                   Ice%slow_domain)
    elseif (associated(Ice%sCS)) then
      call redistribute_FIA_to_FIA(FIA_null, Ice%sCS%FIA, Ice%fast_domain, &
                                   Ice%slow_domain, Ice%sCS%G, Ice%sCS%IG)
      call redistribute_IST_to_IST(IST_null, Ice%sCS%IST, Ice%fast_domain, &
                                   Ice%slow_domain)
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
subroutine set_ocean_top_fluxes(Ice, IST, IOF, FIA, OSS, G, IG, sCS)
  type(ice_data_type),        intent(inout) :: Ice
  type(ice_state_type),       intent(inout) :: IST
  type(ice_ocean_flux_type),  intent(in)    :: IOF
  type(fast_ice_avg_type),    intent(in)    :: FIA
  type(ocean_sfc_state_type), intent(in)    :: OSS
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(in)    :: IG
  type(SIS_slow_CS),          intent(in)    :: sCS

  real :: I_count
  integer :: i, j, k, isc, iec, jsc, jec, m, n
  integer :: i2, j2, i_off, j_off, ind, ncat, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  if (sCS%debug) then
    call Ice_public_type_chksum("Start set_ocean_top_fluxes", Ice, check_slow=.true.)
    call IOF_chksum("Start set_ocean_top_fluxes", IOF, G)
    call FIA_chksum("Start set_ocean_top_fluxes", FIA, G)
  endif

  ! This block of code is probably unneccessary.
  Ice%flux_t(:,:) = 0.0 ; Ice%flux_q(:,:) = 0.0
  Ice%flux_sw_nir_dir(:,:) = 0.0 ; Ice%flux_sw_nir_dif(:,:) = 0.0
  Ice%flux_sw_vis_dir(:,:) = 0.0 ; Ice%flux_sw_vis_dif(:,:) = 0.0
  Ice%flux_lw(:,:) = 0.0 ; Ice%flux_lh(:,:) = 0.0
  Ice%fprec(:,:) = 0.0 ; Ice%lprec(:,:) = 0.0
  do n=1,Ice%ocean_fluxes%num_bcs ; do m=1,Ice%ocean_fluxes%bc(n)%num_fields
    Ice%ocean_fluxes%bc(n)%field(m)%values(:,:) = 0.0
  enddo ; enddo

  ! Sum the concentration weighted mass.
  Ice%mi(:,:) = 0.0
  i_off = LBOUND(Ice%mi,1) - G%isc ; j_off = LBOUND(Ice%mi,2) - G%jsc
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,Ice,IST,IG,i_off,j_off) &
!$OMP                           private(i2,j2)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
    Ice%mi(i2,j2) = Ice%mi(i2,j2) + IST%part_size(i,j,k) * &
        (IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k)))
  enddo ; enddo ; enddo

  if (sCS%do_icebergs .and. associated(IOF%mass_berg)) then
    ! Note that the IOF berg fields and Ice fields are only allocated on the
    ! computational domains, although they may use different indexing conventions.
    Ice%mi(:,:) = Ice%mi(:,:) + IOF%mass_berg(:,:)
    if  (sCS%pass_iceberg_area_to_ocean) then
      Ice%mass_berg(:,:) = IOF%mass_berg(:,:)
      if (associated(IOF%ustar_berg)) Ice%ustar_berg(:,:) = IOF%ustar_berg(:,:)
      if (associated(IOF%area_berg))  Ice%area_berg(:,:) = IOF%area_berg(:,:)
    endif
  endif

  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,Ice,IST,IOF,FIA,i_off,j_off,G,OSS) &
!$OMP                           private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
    Ice%flux_u(i2,j2) = IOF%flux_u_ocn(i,j)
    Ice%flux_v(i2,j2) = IOF%flux_v_ocn(i,j)
    Ice%flux_t(i2,j2) = IOF%flux_t_ocn_top(i,j)
    Ice%flux_q(i2,j2) = IOF%flux_q_ocn_top(i,j)
    Ice%flux_sw_vis_dir(i2,j2) = IOF%flux_sw_vis_dir_ocn(i,j)
    Ice%flux_sw_vis_dif(i2,j2) = IOF%flux_sw_vis_dif_ocn(i,j)
    Ice%flux_sw_nir_dir(i2,j2) = IOF%flux_sw_nir_dir_ocn(i,j)
    Ice%flux_sw_nir_dif(i2,j2) = IOF%flux_sw_nir_dif_ocn(i,j)
    Ice%flux_lw(i2,j2) = IOF%flux_lw_ocn_top(i,j)
    Ice%flux_lh(i2,j2) = IOF%flux_lh_ocn_top(i,j)
    Ice%fprec(i2,j2) = IOF%fprec_ocn_top(i,j)
    Ice%lprec(i2,j2) = IOF%lprec_ocn_top(i,j)
    Ice%runoff(i2,j2)  = FIA%runoff(i,j)
    Ice%calving(i2,j2) = FIA%calving(i,j)
    Ice%runoff_hflx(i2,j2)  = FIA%runoff_hflx(i,j)
    Ice%calving_hflx(i2,j2) = FIA%calving_hflx(i,j)
    Ice%flux_salt(i2,j2) = IOF%flux_salt(i,j)
    Ice%SST_C(i2,j2) = OSS%SST_C(i,j)

    if (IOF%slp2ocean) then
      Ice%p_surf(i2,j2) = FIA%p_atm_surf(i,j) - 1e5 ! SLP - 1 std. atmosphere, in Pa.
    else
      Ice%p_surf(i2,j2) = 0.0
    endif
    Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) + G%G_Earth*Ice%mi(i2,j2)
  enddo ; enddo
  if (allocated(IOF%melt_nudge)) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
      Ice%lprec(i2,j2) = Ice%lprec(i2,j2) + IOF%melt_nudge(i,j)
    enddo ; enddo
  endif

  ind = 0
  do n=1,Ice%ocean_fluxes%num_bcs ; do m=1,Ice%ocean_fluxes%bc(n)%num_fields
    ind = ind + 1
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off  ! Use these to correct for indexing differences.
        Ice%ocean_fluxes%bc(n)%field(m)%values(i2,j2) = IOF%tr_flux_ocn_top(i,j,ind)
    enddo ; enddo
  enddo ; enddo


! This extra block is required with the Verona and earlier versions of the coupler.
  i_off = LBOUND(Ice%part_size,1) - G%isc ; j_off = LBOUND(Ice%part_size,2) - G%jsc
  if (Ice%shared_slow_fast_PEs) then
    if ((Ice%fCS%G%iec-Ice%fCS%G%isc==iec-isc) .and. &
        (Ice%fCS%G%jec-Ice%fCS%G%jsc==jec-jsc)) then
      ! The fast and slow ice PEs are using the same PEs and layout, so the
      ! part_size arrays can be copied directly from the fast ice PEs.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,Ice,IST,i_off,j_off) &
!$OMP                           private(i2,j2)
      do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
        i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
        Ice%part_size(i2,j2,k+1) = IST%part_size(i,j,k)
      enddo ; enddo ; enddo
    endif
  endif

  if (sCS%debug) then
    call Ice_public_type_chksum("End set_ocean_top_fluxes", Ice, check_slow=.true.)
    call IST_chksum("Start set_ocean_top_fluxes", IST, G, IG)
  endif

end subroutine set_ocean_top_fluxes

! Coupler interface to provide ocean surface data to atmosphere.
!
!> update_ice_model_slow_up prepares the ice surface data for forcing the atmosphere
!! and also unpacks the data from the ocean and shares it between the fast and
!! slow ice structures.
subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
  type(ocean_ice_boundary_type), &
    intent(inout) :: Ocean_boundary  !< A structure containing information about
                                   !! the ocean that is being shared with the sea-ice.
  type(ice_data_type), &
    intent(inout) :: Ice           !< The publicly visible ice data type.

  if (.not.associated(Ice%fCS)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS must be associated in update_ice_model_slow_up.")
  if (.not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "The pointer to Ice%sCS must be associated in update_ice_model_slow_up.")

  call unpack_ocn_ice_bdry(Ocean_boundary, Ice%sCS%OSS, Ice%sCS%IST%ITV, Ice%sCS%G, &
                           Ice%sCS%specified_ice, Ice%ocean_fields)

  call translate_OSS_to_sOSS(Ice%sCS%OSS, Ice%sCS%IST, Ice%sCS%sOSS, Ice%sCS%G)

  call exchange_slow_to_fast_ice(Ice)

  call set_ice_surface_fields(Ice)

end subroutine update_ice_model_slow_up

!> This subroutine copies information from the slow part of the sea-ice to the
!! fast part of the sea ice.
subroutine exchange_slow_to_fast_ice(Ice)
  type(ice_data_type), &
    intent(inout) :: Ice            !< The publicly visible ice data type whose slow
                                    !! part is to be exchanged with the fast part.
  type(simple_OSS_type), pointer :: sOSS_null => NULL()
  type(ice_state_type),  pointer :: IST_null => NULL()

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(ice_clock_exchange)


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

  call mpp_clock_end(ice_clock_exchange) ; call mpp_clock_end(iceClock)

end subroutine exchange_slow_to_fast_ice

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

  call unpack_ocn_ice_bdry(Ocean_boundary, Ice%sCS%OSS, Ice%sCS%IST%ITV, Ice%sCS%G, &
                           Ice%sCS%specified_ice, Ice%ocean_fields)

  call translate_OSS_to_sOSS(Ice%sCS%OSS, Ice%sCS%IST, Ice%sCS%sOSS, Ice%sCS%G)

end subroutine unpack_ocean_ice_boundary

!> This subroutine converts the information in a publicly visible
!! ocean_ice_boundary_type into an internally visible ocean_sfc_state_type
!! variable.
subroutine unpack_ocn_ice_bdry(OIB, OSS, ITV, G, specified_ice, ocean_fields)
  type(ocean_ice_boundary_type), intent(in)    :: OIB
  type(ocean_sfc_state_type),    intent(inout) :: OSS
  type(ice_thermo_type),         intent(in)    :: ITV
  type(SIS_hor_grid_type),       intent(inout) :: G
  logical,                       intent(in)    :: specified_ice ! If true, use specified ice properties.
  type(coupler_3d_bc_type),      intent(inout) :: ocean_fields  ! A structure of ocean fields, often
                                                                ! related to passive tracers.

  real, dimension(G%isd:G%ied, G%jsd:G%jed) :: u_nonsym, v_nonsym
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  logical :: Cgrid_ocn
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, index
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = LBOUND(OIB%t,1) - G%isc ; j_off = LBOUND(OIB%t,2) - G%jsc

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(ice_clock_slow)

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,OSS,OIB,ITV,i_off,j_off) &
!$OMP                           private(i2,j2)
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    OSS%s_surf(i,j) = OIB%s(i2,j2)
    OSS%T_fr_ocn(i,j) = T_Freeze(OSS%s_surf(i,j), ITV)
    OSS%frazil(i,j) = OIB%frazil(i2,j2)
    OSS%sea_lev(i,j) = OIB%sea_level(i2,j2)
  enddo ; enddo

  ! Pass the ocean state through ice on partition 0, unless using specified ice.
  if (.not. specified_ice) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,OSS,OIB,i_off,j_off) &
!$OMP                           private(i2,j2)
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      OSS%SST_C(i,j) = OIB%t(i2,j2) - T_0degC
    enddo ; enddo
  endif

  Cgrid_ocn = (allocated(OSS%u_ocn_C) .and. allocated(OSS%v_ocn_C))

  ! Unpack the ocean surface velocities.
  if (OIB%stagger == AGRID) then
    u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      u_nonsym(i,j) = OIB%u(i2,j2) ; v_nonsym(i,j) = OIB%v(i2,j2)
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
        u_nonsym(i,j) = OIB%u(i2,j2) ; v_nonsym(i,j) = OIB%v(i2,j2)
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
        OSS%u_ocn_B(I,J) = OIB%u(i2,j2)
        OSS%v_ocn_B(I,J) = OIB%v(i2,j2)
      enddo ; enddo
      if (G%symmetric) &
        call fill_symmetric_edges(OSS%u_ocn_B, OSS%v_ocn_B, G%Domain, stagger=BGRID_NE)
    endif

  elseif (OIB%stagger == CGRID_NE) then
    if (Cgrid_ocn) then
      do j=jsc,jec ; do I=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        OSS%u_ocn_C(I,j) = OIB%u(i2,j2)
      enddo ; enddo
      do J=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        OSS%v_ocn_C(i,J) = OIB%v(i2,j2)
      enddo ; enddo
      if (G%symmetric) &
        call fill_symmetric_edges(OSS%u_ocn_C, OSS%v_ocn_C, G%Domain, stagger=CGRID_NE)
    else
      u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
      do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        u_nonsym(I,j) = OIB%u(i2,j2) ; v_nonsym(i,J) = OIB%v(i2,j2)
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
  if (OSS%num_tr<0) then
    ! Determine the number of tracer arrays and allocate them.
    OSS%num_tr = 0
    do n=1,OIB%fields%num_bcs
      OSS%num_tr = OSS%num_tr + OIB%fields%bc(n)%num_fields
    enddo
    if (OSS%num_tr > 0) then
      allocate(OSS%tr_array(G%isd:G%ied,G%jsd:G%jed,OSS%num_tr)) ; OSS%tr_array(:,:,:) = 0.0
    endif
  endif
  index = 0
  do n=1,OIB%fields%num_bcs  ; do m=1,OIB%fields%bc(n)%num_fields
    index = index+1
    if (index == 1) then
      i_off = LBOUND(OIB%fields%bc(n)%field(m)%values,1) - G%isc
      j_off = LBOUND(OIB%fields%bc(n)%field(m)%values,2) - G%jsc
    endif
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      OSS%tr_array(i,j,index) = OIB%fields%bc(n)%field(m)%values(i2,j2)
    enddo ; enddo
  enddo ; enddo

  call mpp_clock_end(ice_clock_slow) ; call mpp_clock_end(iceClock)

end subroutine unpack_ocn_ice_bdry

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ice_surface_fields prepares the ice surface state for atmosphere fast
!! physics and does precalculation of ice radiative properties.
subroutine set_ice_surface_fields(Ice)
  type(ice_data_type), intent(inout) :: Ice !< The publicly visible ice data type whose
                                            !! surface properties are being set.

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(ice_clock_fast)
  if (.not.associated(Ice%fCS)) call SIS_error(FATAL, &
      "The pointer to Ice%fCS must be associated in set_ice_surface_fields.")

  call set_ice_surface_state(Ice, Ice%fCS%IST, Ice%fCS%sOSS, Ice%fCS%Rad, &
                             Ice%fCS%FIA, Ice%fCS%G, Ice%fCS%IG, Ice%fCS )

  call mpp_clock_end(ice_clock_fast) ; call mpp_clock_end(iceClock)
end subroutine set_ice_surface_fields

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ice_surface_state - prepare surface state for atmosphere fast physics    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ice_surface_state(Ice, IST, OSS, Rad, FIA, G, IG, fCS)
  type(ice_data_type),        intent(inout) :: Ice
  type(ice_state_type),       intent(inout) :: IST
  type(simple_OSS_type),      intent(in)    :: OSS
  type(ice_rad_type),         intent(inout) :: Rad
  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(in)    :: IG
  type(SIS_fast_CS),          intent(inout) :: fCS

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: m_ice_tot
  real, dimension(IG%NkIce) :: sw_abs_lay
  real :: u, v
  real :: area_pt
  real :: rho_ice  ! The nominal density of sea ice in kg m-3.
  real :: rho_snow ! The nominal density of snow in kg m-3.
  type(time_type) :: dt_r   ! A temporary radiation timestep.

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  integer :: index
  real :: H_to_m_ice     ! The specific volumes of ice and snow times the
  real :: H_to_m_snow    ! conversion factor from thickness units, in m H-1.
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  logical :: slab_ice    ! If true, use the very old slab ice thermodynamics,
                         ! with effectively zero heat capacity of ice and snow.

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc

  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice, rho_snow=rho_snow, slab_ice=slab_ice)
  H_to_m_snow = IG%H_to_kg_m2 / Rho_snow ; H_to_m_ice = IG%H_to_kg_m2 / Rho_ice


  if (fCS%bounds_check) &
    call IST_bounds_check(IST, G, IG, "Start of set_ice_surface_state", Rad=Rad) !, OSS=OSS)

  if (fCS%debug) then
    call IST_chksum("Start set_ice_surface_state", IST, G, IG)
    call Ice_public_type_chksum("Start set_ice_surface_state", Ice, check_fast=.true.)
  endif

  m_ice_tot(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,IST,OSS,FIA,ncat,m_ice_tot)
  do j=jsc,jec

    do k=1,ncat ; do i=isc,iec
      FIA%tmelt(i,j,k) = 0.0 ; FIA%bmelt(i,j,k) = 0.0
      m_ice_tot(i,j) = m_ice_tot(i,j) + IST%mH_ice(i,j,k) * IST%part_size(i,j,k)
    enddo ; enddo

    do i=isc,iec
      if (m_ice_tot(i,j) > 0.0) then
        FIA%bheat(i,j) = OSS%bheat(i,J)
      else
        FIA%bheat(i,j) = 0.0
      endif
    enddo
  enddo

  if (.not.fCS%Eulerian_tsurf) then
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      Rad%t_skin(i,j,k) = IST%t_surf(i,j,k) - T_0degC
    enddo ; enddo ; enddo
  endif

  ! Determine the sea-ice optical properties.

  !   These initialization calls for ice-free categories are not really
  ! needed because these arrays are only used where there is ice.
  ! In the case of slab_ice, the various Rad%sw_abs arrays are initialized
  ! to 0 when they are allocated, and this never changes.
  ! The following lines can be uncommented without changing answers.
  ! Rad%sw_abs_sfc(:,:,:) = 0.0 ; Rad%sw_abs_snow(:,:,:) = 0.0
  ! Rad%sw_abs_ice(:,:,:,:) = 0.0 ; Rad%sw_abs_ocn(:,:,:) = 0.0
  ! Rad%sw_abs_int(:,:,:) = 0.0
  ! Ice%albedo(:,:,:) = 0.0
  ! Ice%albedo_vis_dir(:,:,:) = 0.0 ; Ice%albedo_vis_dif(:,:,:) = 0.0
  ! Ice%albedo_nir_dir(:,:,:) = 0.0 ; Ice%albedo_nir_dif(:,:,:) = 0.0

  ! Set the initial ocean albedos, either using coszen_nextrad or a
  ! synthetic sun angle.
  dT_r = fCS%Time_step_slow
  if (Rad%frequent_albedo_update) dT_r = fCS%Time_step_fast
  call set_ocean_albedo(Ice, Rad%do_sun_angle_for_alb, G, fCS%Time, &
                        fCS%Time + dT_r, Rad%coszen_nextrad)

  if (slab_ice) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Rad,Ice,i_off,j_off, &
!$OMP                                  H_to_m_snow,H_to_m_ice,OSS) &
!$OMP                          private(i2,j2,k2)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k) > 0.0) then
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      call slab_ice_optics(IST%mH_snow(i,j,k)*H_to_m_snow, IST%mH_ice(i,j,k)*H_to_m_ice, &
               Rad%t_skin(i,j,k), OSS%T_fr_ocn(i,j), &
               Ice%albedo(i2,j2,k2))

      Ice%albedo_vis_dir(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_vis_dif(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_nir_dir(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_nir_dif(i2,j2,k2) = Ice%albedo(i2,j2,k2)
    endif ; enddo ; enddo ; enddo
  else
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Ice,G,IG,i_off,j_off, &
!$OMP                                  H_to_m_snow,H_to_m_ice,OSS,Rad,fCS) &
!$OMP                          private(i2,j2,k2,sw_abs_lay)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k) > 0.0) then
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      call ice_optics_SIS2(IST%mH_pond(i,j,k), IST%mH_snow(i,j,k)*H_to_m_snow, &
               IST%mH_ice(i,j,k)*H_to_m_ice, &
               Rad%t_skin(i,j,k), OSS%T_fr_ocn(i,j), IG%NkIce, &
               Ice%albedo_vis_dir(i2,j2,k2), Ice%albedo_vis_dif(i2,j2,k2), &
               Ice%albedo_nir_dir(i2,j2,k2), Ice%albedo_nir_dif(i2,j2,k2), &
               Rad%sw_abs_sfc(i,j,k),  Rad%sw_abs_snow(i,j,k), &
               sw_abs_lay, Rad%sw_abs_ocn(i,j,k), Rad%sw_abs_int(i,j,k), &
               fCS%optics_CSp, IST%ITV, coszen_in=Rad%coszen_nextrad(i,j))

      do m=1,IG%NkIce ; Rad%sw_abs_ice(i,j,k,m) = sw_abs_lay(m) ; enddo

      !Niki: Is the following correct for diagnostics?
      !   Probably this calculation of the "average" albedo should be replaced
      ! with a calculation that weights the averaging by the fraction of the
      ! shortwave radiation in each wavelength and orientation band.  However,
      ! since this is only used for diagnostic purposes, making this change
      ! might not be too urgent. -RWH
      Ice%albedo(i2,j2,k2) = (Ice%albedo_vis_dir(i2,j2,k2)+Ice%albedo_nir_dir(i2,j2,k2)&
                        +Ice%albedo_vis_dif(i2,j2,k2)+Ice%albedo_nir_dif(i2,j2,k2))/4

    endif ; enddo ; enddo ; enddo
  endif

  if (fCS%bounds_check) &
    call IST_bounds_check(IST, G, IG, "Midpoint set_ice_surface_state", Rad=Rad) !, OSS=OSS)

  ! Copy the surface temperatures into the externally visible data type.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,ncat,i_off,j_off,OSS) &
!$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    Ice%t_surf(i2,j2,1) = OSS%SST_C(i,j) + T_0degC
    Ice%part_size(i2,j2,1) = IST%part_size(i,j,0)
  enddo ; enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Rad,Ice,ncat,i_off,j_off,OSS) &
!$OMP                          private(i2,j2,k2)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      Ice%t_surf(i2,j2,k2) = Rad%t_skin(i,j,k) + T_0degC
      Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
    enddo ; enddo
  enddo

  ! Put ocean salinity and ocean and ice velocities into Ice%u_surf/v_surf on t-cells.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,G,i_off,j_off,OSS) &
!$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%u_surf(i2,j2,1) = OSS%u_ocn_A(i,j) ; Ice%v_surf(i2,j2,1) = OSS%v_ocn_A(i,j)
    Ice%u_surf(i2,j2,2) = OSS%u_ice_A(i,j) ; Ice%v_surf(i2,j2,2) = OSS%v_ice_A(i,j)
    Ice%s_surf(i2,j2) = OSS%s_surf(i,j)
  enddo ; enddo

  if (fCS%debug) then
    call chksum(Ice%u_surf(:,:,1), "Intermed Ice%u_surf(1)")
    call chksum(Ice%v_surf(:,:,1), "Intermed Ice%v_surf(1)")
    call chksum(Ice%u_surf(:,:,2), "Intermed Ice%u_surf(2)")
    call chksum(Ice%v_surf(:,:,2), "Intermed Ice%v_surf(2)")
    call chksum(G%mask2dT(isc:iec,jsc:jec), "Intermed G%mask2dT")
!   if (allocated(OSS%u_ocn_C)) &
!     call uchksum(OSS%u_ocn_C, "OSS%u_ocn_C", G%HI, haloshift=1)
!   if (allocated(OSS%v_ocn_C)) &
!     call vchksum(OSS%v_ocn_C, "OSS%v_ocn_C", G%HI, haloshift=1)
!   if (allocated(OSS%u_ocn_B)) &
!     call Bchksum(OSS%u_ocn_B, "OSS%u_ocn_B", G%HI, haloshift=1)
!   if (allocated(OSS%v_ocn_B)) &
!     call Bchksum(OSS%v_ocn_B, "OSS%v_ocn_B", G%HI, haloshift=1)
    call chksum(G%sin_rot(isc:iec,jsc:jec), "G%sin_rot")
    call chksum(G%cos_rot(isc:iec,jsc:jec), "G%cos_rot")
  endif

  ! Rotate the velocities from the ocean coordinates to lat/lon coordiantes.
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

  ! Copy over the additional tracer fields into the ocean_fields structure.
  index = 0
  do n=1,Ice%ocean_fields%num_bcs  ; do m=1,Ice%ocean_fields%bc(n)%num_fields
    index = index+1
    if (index == 1) then
      i_off = LBOUND(Ice%ocean_fields%bc(n)%field(m)%values,1) - G%isc
      j_off = LBOUND(Ice%ocean_fields%bc(n)%field(m)%values,2) - G%jsc
    endif
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      Ice%ocean_fields%bc(n)%field(m)%values(i2,j2,1) = OSS%tr_array(i,j,index)
    enddo ; enddo
  enddo ; enddo

  if (fCS%debug) then
    call IST_chksum("End set_ice_surface_state", IST, G, IG)
    call Ice_public_type_chksum("End set_ice_surface_state", Ice, check_fast=.true.)
  endif

  if (fCS%bounds_check) &
    call Ice_public_type_bounds_check(Ice, G, "End set_ice_surface_state")

end subroutine set_ice_surface_state

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! update_ice_model_fast - records fluxes (in Ice) and calculates ice temp. on  !
!                         (fast) atmospheric timestep (see coupler_main.f90)   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine update_ice_model_fast( Atmos_boundary, Ice )
  type(ice_data_type),              intent(inout) :: Ice
  type(atmos_ice_boundary_type),    intent(inout) :: Atmos_boundary

  type(time_type) :: Time_start, Time_end, dT_fast

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(ice_clock_fast)

  if (Ice%fCS%debug) &
    call Ice_public_type_chksum("Pre do_update_ice_model_fast", Ice, check_fast=.true.)

  dT_fast = Ice%fCS%Time_step_fast
  Time_start = Ice%fCS%Time
  Time_end = Time_start + dT_fast

  if (Ice%fCS%Rad%add_diurnal_sw) &
    call add_diurnal_sw(Atmos_boundary, Ice%fCS%G, Time_start, Time_end)

  call do_update_ice_model_fast(Atmos_boundary, Ice%fCS%IST, Ice%fCS%sOSS, Ice%fCS%Rad, &
                                Ice%fCS%FIA, dT_fast, Ice%fCS%fast_thermo_CSp, &
                                Ice%fCS%G, Ice%fCS%IG )

  ! Advance the master sea-ice time.
  Ice%fCS%Time = Ice%fCS%Time + dT_fast

  Ice%Time = Ice%fCS%Time

  call fast_radiation_diagnostics(Atmos_boundary, Ice, Ice%fCS%IST, Ice%fCS%Rad, &
                                  Ice%fCS%FIA, Ice%fCS%G, Ice%fCS%IG, Ice%fCS, &
                                  Time_start, Time_end)

  ! Set some of the evolving ocean properties that will be seen by the
  ! atmosphere in the next time-step.
  call set_fast_ocean_sfc_properties(Atmos_boundary, Ice, Ice%fCS%IST, Ice%fCS%Rad, &
                                     Ice%fCS%FIA, Ice%fCS%G, Ice%fCS%IG, Time_end, Time_end + dT_fast)

  if (Ice%fCS%debug) &
    call Ice_public_type_chksum("End do_update_ice_model_fast", Ice, check_fast=.true.)
  if (Ice%fCS%bounds_check) &
    call Ice_public_type_bounds_check(Ice, Ice%fCS%G, "End update_ice_fast")

  call mpp_clock_end(ice_clock_fast) ; call mpp_clock_end(iceClock)

end subroutine update_ice_model_fast

subroutine set_fast_ocean_sfc_properties( Atmos_boundary, Ice, IST, Rad, FIA, &
                                          G, IG, Time_start, Time_end)
  type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary
  type(ice_data_type),           intent(inout) :: Ice
  type(ice_state_type),          intent(inout) :: IST
  type(ice_rad_type),            intent(inout) :: Rad
  type(fast_ice_avg_type),       intent(inout) :: FIA
  type(SIS_hor_grid_type),       intent(inout) :: G
  type(ice_grid_type),           intent(inout) :: IG
  type(time_type),               intent(in)    :: Time_start, Time_end

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
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
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,Ice,IST,Atmos_boundary,Rad,&
!$OMP                                  FIA,io_A,jo_A,io_I,jo_I ) &
!$OMP                          private(i3,j3)
  do j=jsc,jec ; do i=isc,iec
    i3 = i+io_A ; j3 = j+jo_A
    Rad%coszen_nextrad(i,j) = Atmos_boundary%coszen(i3,j3,1)
    FIA%p_atm_surf(i,j) = Atmos_boundary%p(i3,j3,1)
  enddo ; enddo

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,Rad,Ice,io_I,jo_I ) &
!$OMP                           private(i2,j2,k2)
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+io_I ; j2 = j+jo_I ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = Rad%t_skin(i,j,k) + T_0degC
  enddo ; enddo ; enddo

  ! set_ocean_albedo only needs to be called if do_sun_angle_for_alb is true or
  ! if the coupled model's radiation timestep is shorter than the slow coupling
  ! timestep.  However, it is safe (if wasteful) to call it more frequently.
  if (Rad%frequent_albedo_update) then
    call set_ocean_albedo(Ice, Rad%do_sun_angle_for_alb, G, Time_start, &
                          Time_end, Rad%coszen_nextrad)
  endif

end subroutine set_fast_ocean_sfc_properties

!> set_ocean_albedo uses either the time or the input cosine of solar zenith
!!   angle to calculate the ocean albedo.
subroutine set_ocean_albedo(Ice, recalc_sun_angle, G, Time_start, Time_end, coszen)
  type(ice_data_type),     intent(inout) :: Ice
  logical,                 intent(in)    :: recalc_sun_angle
  type(SIS_hor_grid_type), intent(inout) :: G
  type(time_type),         intent(in)    :: Time_start, Time_end
  real, dimension(G%isd:G%ied, G%jsd:G%jed), &
                           intent(in)    :: coszen

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: &
    dummy, &  ! A dummy array that is not used again.
    cosz_alb  ! The cosine of the solar zenith angle for calculating albedo, ND.
  real :: rad
  real :: rrsun_dt_ice
  type(time_type) :: dT_ice   ! The time interval for this update.
  integer :: i, j, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  rad = acos(-1.)/180.
  dT_ice = Time_end - Time_start

  if (recalc_sun_angle) then
    call diurnal_solar(G%geoLatT(isc:iec,jsc:jec)*rad, G%geoLonT(isc:iec,jsc:jec)*rad, &
                 Time_start, cosz=cosz_alb, fracday=dummy, rrsun=rrsun_dt_ice, &
                 dt_time=dT_ice)
  else
    do j=jsc,jec ; do i=isc,iec ; cosz_alb(i,j) = coszen(i,j) ; enddo ; enddo
  endif
  call compute_ocean_albedo(Ice%ocean_pt, cosz_alb(:,:), Ice%albedo_vis_dir(:,:,1),&
                            Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                            Ice%albedo_nir_dif(:,:,1), rad*G%geoLatT(isc:iec,jsc:jec) )

end subroutine set_ocean_albedo


subroutine fast_radiation_diagnostics(ABT, Ice, IST, Rad, FIA, G, IG, CS, &
                                      Time_start, Time_end)
  type(atmos_ice_boundary_type), intent(in)    :: ABT
  type(ice_data_type),           intent(in)    :: Ice
  type(ice_state_type),          intent(in)    :: IST
  type(ice_rad_type),            intent(in)    :: Rad
  type(fast_ice_avg_type),       intent(inout) :: FIA
  type(SIS_hor_grid_type),       intent(in)    :: G
  type(ice_grid_type),           intent(in)    :: IG
  type(SIS_fast_CS),             intent(inout) :: CS
  type(time_type),               intent(in)    :: Time_start, Time_end

  real, dimension(G%isd:G%ied, G%jsd:G%jed) :: tmp_diag, sw_dn, net_sw, avg_alb
  real, dimension(G%isd:G%ied) :: Tskin_avg, ice_conc
  real :: dt_diag
  real    :: Stefan ! The Stefan-Boltzmann constant in W m-2 K-4 as used for
                    ! strictly diagnostic purposes.
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  integer :: i, j, k, m, i2, j2, k2, i3, j3, isc, iec, jsc, jec, ncat, NkIce
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
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Rad,tmp_diag)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * &
                     (Rad%sw_abs_ocn(i,j,k) + Rad%sw_abs_int(i,j,k))
    enddo ; enddo ; enddo
    call post_data(Rad%id_sw_pen, tmp_diag, CS%diag)
  endif

  if (Rad%id_lwdn > 0) then
    tmp_diag(:,:) = 0.0
    Stefan = 5.6734e-8  ! Set the Stefan-Bolzmann constant, in W m-2 K-4.
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
      i3 = i+io_A ; j3 = j+jo_A ; k2 = k+1
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * &
                           (ABT%lw_flux(i3,j3,k2) + Stefan*(Rad%t_skin(i,j,k)+T_0degC)**4)
    endif ; enddo ; enddo ; enddo
    call post_data(Rad%id_lwdn, tmp_diag, CS%diag)
  endif

  sw_dn(:,:) = 0.0 ; net_sw(:,:) = 0.0 ; avg_alb(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,IST,Ice,ABT, &
!$OMP                                  io_I,jo_I,io_A,jo_A,sw_dn,net_sw,avg_alb) &
!$OMP                          private(i2,j2,k2,i3,j3)
  do j=jsc,jec ; do k=0,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
    i2 = i+io_I ; j2 = j+jo_I ; i3 = i+io_A ; j3 = j+jo_A ; k2 = k+1
    if (associated(ABT%sw_down_vis_dir)) then
      sw_dn(i,j) = sw_dn(i,j) + IST%part_size(i,j,k) * ( &
            (ABT%sw_down_vis_dir(i3,j3,k2) + ABT%sw_down_vis_dif(i3,j3,k2)) + &
            (ABT%sw_down_nir_dir(i3,j3,k2) + ABT%sw_down_nir_dif(i3,j3,k2)) )
    else
      sw_dn(i,j) = sw_dn(i,j) + IST%part_size(i,j,k) * ( &
            (ABT%sw_flux_vis_dir(i3,j3,k2)/(1-Ice%albedo_vis_dir(i2,j2,k2)) + &
             ABT%sw_flux_vis_dif(i3,j3,k2)/(1-Ice%albedo_vis_dif(i2,j2,k2))) + &
            (ABT%sw_flux_nir_dir(i3,j3,k2)/(1-Ice%albedo_nir_dir(i2,j2,k2)) + &
             ABT%sw_flux_nir_dif(i3,j3,k2)/(1-Ice%albedo_nir_dif(i2,j2,k2))) )
    endif

    net_sw(i,j) = net_sw(i,j) + IST%part_size(i,j,k) * ( &
          (ABT%sw_flux_vis_dir(i3,j3,k2) + ABT%sw_flux_vis_dif(i3,j3,k2)) + &
          (ABT%sw_flux_nir_dir(i3,j3,k2) + ABT%sw_flux_nir_dif(i3,j3,k2)) )
    avg_alb(i,j) = avg_alb(i,j) + IST%part_size(i,j,k) * 0.25 * ( &
            (Ice%albedo_vis_dir(i2,j2,k2) + Ice%albedo_vis_dif(i2,j2,k2)) + &
            (Ice%albedo_nir_dir(i2,j2,k2) + Ice%albedo_nir_dif(i2,j2,k2)) )
    ! Consider recalculating this as avg_alb(i,j) = 1.0 - net_sw(i,j) / sw_dn(i,j) ? -RWH
  endif ; enddo ; enddo ; enddo

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,Rad,IST,FIA) &
!$OMP                          private(Tskin_avg, ice_conc)
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
    do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
      if (sw_dn(i,j) > 0.0) &
        avg_alb(i,j) = (sw_dn(i,j) - net_sw(i,j)) / sw_dn(i,j)
      ! Otherwise keep the simple average that was set above.
    endif ; enddo ; enddo
    call post_data(Rad%id_alb, avg_alb, CS%diag)
  endif

  do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
    FIA%flux_sw_dn(i,j) = FIA%flux_sw_dn(i,j) + sw_dn(i,j)
  endif ; enddo ; enddo

  if (Rad%id_coszen>0) call post_data(Rad%id_coszen, Rad%coszen_nextrad, CS%diag)

  call disable_SIS_averaging(CS%diag)

end subroutine fast_radiation_diagnostics

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
subroutine ice_model_init(Ice, Time_Init, Time, Time_step_fast, Time_step_slow, Verona_coupler )

  type(ice_data_type), intent(inout) :: Ice            !< The ice data type that is being initialized.
  type(time_type)    , intent(in)    :: Time_Init      !< The starting time of the model integration
  type(time_type)    , intent(in)    :: Time           !< The current time
  type(time_type)    , intent(in)    :: Time_step_fast !< The time step for the ice_model_fast
  type(time_type)    , intent(in)    :: Time_step_slow !< The time step for the ice_model_slow
  logical,   optional, intent(in)    :: Verona_coupler !< If present and false, use the input values
                                              !! in Ice to determine whether this is a fast or slow
                                              !! ice processor or both.  Otherwise, carry out all of
                                              !! the sea ice iniatialization calls so that SIS2 will
                                              !! work with the Verona and earlier releases of the FMS
                                              !! coupler code in configurations that use the exchange
                                              !! grid to communicate with the atmosphere or land.

! This include declares and sets the variable "version".
#include "version_variable.h"
  real :: enth_spec_snow, enth_spec_ice
  real, allocatable :: S_col(:)
  real :: pi ! pi = 3.1415926... calculated as 4*atan(1)
  integer :: i, j, k, l, i2, j2, k2, i_off, j_off, n
  integer :: isc, iec, jsc, jec, nCat_dflt
  logical :: spec_thermo_sal
  character(len=120) :: restart_file, fast_rest_file
  character(len=240) :: restart_path, fast_rest_path
  character(len=40)  :: mod = "ice_model" ! This module's name.
  character(len=8)   :: nstr
  type(directories)  :: dirs   ! A structure containing several relevant directory paths.

  type(param_file_type) :: param_file
  type(hor_index_type)  :: fHI  !  A hor_index_type for array extents on fast_ice_PEs.
  type(hor_index_type)  :: sHI  !  A hor_index_type for array extents on slow_ice_PEs.

  type(dyn_horgrid_type),  pointer :: dG => NULL()
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
  real :: mom_rough_ice  ! momentum same, cd10=(von_k/ln(10/z0))^2, in m.
  real :: heat_rough_ice ! heat roughness length, in m.
  real :: dt_Rad_real    ! The radiation timestep, in s.
  type(time_type) :: dt_Rad ! The radiation timestep, used initializing albedos.
  real :: rad            ! The conversion factor from degrees to radians.
  real :: rrsun          ! An unused temporary factor related to the Earth-sun distance.

  ! Parameters that properly belong exclusively to ice_thm.
  real :: k_snow         ! snow conductivity (W/mK)
  real :: h_lo_lim       ! The min ice thickness for temp. calc, in m.
  real :: H_to_kg_m2_tmp ! A temporary variable for holding the intended value
                         ! of the thickness to mass-per-unit-area conversion
                         ! factor.
  real :: enth_unit      ! A conversion factor for enthalpy from Joules kg-1.
  real :: massless_ice_enth, massless_snow_enth, massless_ice_salin
  real :: H_rescale_ice, H_rescale_snow ! Rescaling factors to account for
                         ! differences in thickness units between the current
                         ! model run and the input files.

  real, allocatable, dimension(:,:) :: &
    h_ice_input, dummy  ! Temporary arrays.

  real, allocatable, target, dimension(:,:,:,:) :: t_ice_tmp, sal_ice_tmp
  real, allocatable, target, dimension(:,:,:) :: t_snow_tmp
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  real :: g_Earth        !   The gravitational acceleration in m s-2.
  real :: ice_bulk_salin ! The globally constant sea ice bulk salinity, in g/kg
                         ! that is used to calculate the ocean salt flux.
  real :: ice_rel_salin  ! The initial bulk salinity of sea-ice relative to the
                         ! salinity of the water from which it formed, nondim.
  real :: coszen_IC      ! A constant value that is used to initialize
                         ! coszen if it is not read from a restart file, or a
                         ! negative number to use the time and geometry.
  real :: rho_ice        ! The nominal density of sea ice in kg m-3.
  real :: rho_snow       ! The nominal density of snow in kg m-3.
  real :: rho_Ocean      ! The nominal density of seawater, in kg m-3.
  real :: kmelt          ! A constant that is used in the calculation of the
                         ! ocean/ice basal heat flux, in W m-2 K-1.  This could
                         ! be changed to reflect the turbulence in the under-ice
                         ! ocean boundary layer and the effective depth of the
                         ! reported value of t_ocn.
  real :: ocean_part_min ! The minimum value for the fractional open-ocean
                         ! area.  This can be 0, but for some purposes it
                         ! may be useful to set this to a miniscule value
                         ! (like 1e-40) that will be lost to roundoff
                         ! during any sums so that the open ocean fluxes
                         ! can be used in interpolation across categories.

  integer :: CatIce, NkIce, isd, ied, jsd, jed
  integer :: idr, id_sal
  integer :: write_geom
  logical :: nudge_sea_ice
  logical :: atmos_winds, slp2ocean
  logical :: do_icebergs, pass_iceberg_area_to_ocean
  logical :: do_ridging
  logical :: specified_ice
  logical :: Cgrid_dyn, slab_ice
  logical :: debug, bounds_check
  logical :: do_sun_angle_for_alb, add_diurnal_sw
  logical :: init_coszen, init_Tskin, init_rough
  logical :: write_error_mesg
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
  logical :: Verona
  logical :: read_aux_restart
  logical :: split_restart_files
  logical :: is_restart = .false.
  character(len=16)  :: stagger, dflt_stagger

  ! ### These are just here to keep the order of SIS_parameter_doc.
  logical :: column_check
  real :: imb_tol

  if (associated(Ice%sCS)) then ; if (associated(Ice%sCS%IST)) then
    call SIS_error(WARNING, "ice_model_init called with an associated "// &
                    "Ice%sCS%Ice_state structure. Model is already initialized.")
    return
  endif ; endif

  ! For now, both fast and slow processes occur on all sea-ice PEs.
  fast_ice_PE = .true. ; slow_ice_PE = .true.
  if (present(Verona_coupler)) then ; if (.not.Verona_coupler) then
    fast_ice_PE = Ice%fast_ice_pe ; slow_ice_PE = Ice%slow_ice_pe
  endif ; endif
  Verona = .true. ; if (present(Verona_coupler)) Verona = Verona_coupler

  ! Open the parameter file.
  call Get_SIS_Input(param_file, dirs, check_params=slow_ice_PE)

  call callTree_enter("ice_model_init(), ice_model.F90")

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "SPECIFIED_ICE", specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  call get_param(param_file, mod, "CGRID_ICE_DYNAMICS", Cgrid_dyn, &
                 "If true, use a C-grid discretization of the sea-ice \n"//&
                 "dynamics; if false use a B-grid discretization.", &
                 default=.false.)
  if (specified_ice) then
    slab_ice = .true.
    call log_param(param_file, mod, "USE_SLAB_ICE", slab_ice, &
                 "Use the very old slab-style ice.  With SPECIFIED_ICE, \n"//&
                 "USE_SLAB_ICE is always true.")
  else
    call get_param(param_file, mod, "USE_SLAB_ICE", slab_ice, &
                 "If true, use the very old slab-style ice.", default=.false.)
  endif
  call get_param(param_file, mod, "SINGLE_ICE_STATE_TYPE", single_IST, &
                 "If true, the fast and slow portions of the ice use a \n"//&
                 "single common ice_state_type.  Otherwise they point to \n"//&
                 "different ice_state_types that need to be explicitly \n"//&
                 "copied back and forth.", default=.true.)
  call get_param(param_file, mod, "EULERIAN_TSURF", Eulerian_tsurf, &
                 "If true, use previous calculations of the ice-top surface \n"//&
                 "skin temperature for tsurf at the start of atmospheric \n"//&
                 "time stepping, including interpolating between tsurf \n"//&
                 "values from other categories in the same location.", default=.true.)

  call obsolete_logical(param_file, "SIS1_5L_THERMODYNAMICS", warning_val=.false.)
  call obsolete_logical(param_file, "INTERSPERSED_ICE_THERMO", warning_val=.false.)
  call obsolete_logical(param_file, "AREA_WEIGHTED_STRESSES", warning_val=.true.)

  dflt_stagger = "B" ; if (Cgrid_dyn) dflt_stagger = "C"
  call get_param(param_file, mod, "ICE_OCEAN_STRESS_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the \n"//&
                 "staggering of the stress field on the ocean that is \n"//&
                 "returned to the coupler.  Valid values include \n"//&
                 "'A', 'B', or 'C', with a default that follows the \n"//&
                 "value of CGRID_ICE_DYNAMICS.", default=dflt_stagger)

  ! Rho_ocean is not actually used here, but it used from later get_param
  ! calls in other modules.  This call is here to avoid changing the order of
  ! the entries in the SIS_parameter_doc files.
  call get_param(param_file, mod, "RHO_OCEAN", Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0)
  call get_param(param_file, mod, "RHO_ICE", Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0)
  call get_param(param_file, mod, "RHO_SNOW", Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0)

  call get_param(param_file, mod, "G_EARTH", g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default = 9.80)

  call get_param(param_file, mod, "MOMENTUM_ROUGH_ICE", mom_rough_ice, &
                 "The default momentum roughness length scale for the ocean.", &
                 units="m", default=1.0e-4)
  call get_param(param_file, mod, "HEAT_ROUGH_ICE", heat_rough_ice, &
                 "The default roughness length scale for the turbulent \n"//&
                 "transfer of heat into the ocean.", units="m", default=1.0e-4)

  call get_param(param_file, mod, "CONSTANT_COSZEN_IC", coszen_IC, &
                 "A constant value to use to initialize the cosine of \n"//&
                 "the solar zenith angle for the first radiation step, \n"//&
                 "or a negative number to use the current time and astronomy.", &
                 units="nondim", default=-1.0)
  call get_param(param_file, mod, "DT_RADIATION", dt_Rad_real, &
                 "The time step with which the shortwave radiation and \n"//&
                 "fields like albedos are updated.  Currently this is only \n"//&
                 "used to initialize albedos when there is no restart file.", &
                 units="s", default=time_type_to_real(Time_step_slow))
  dt_Rad = real_to_time_type(dt_Rad_real)
  call get_param(param_file, mod, "ICE_KMELT", kmelt, &
                 "A constant giving the proportionality of the ocean/ice \n"//&
                 "base heat flux to the tempature difference, given by \n"//&
                 "the product of the heat capacity per unit volume of sea \n"//&
                 "water times a molecular diffusive piston velocity.", &
                 units="W m-2 K-1", default=6e-5*4e6)
  call get_param(param_file, mod, "SNOW_CONDUCT", k_snow, &
                 "The conductivity of heat in snow.", units="W m-1 K-1", &
                 default=0.31)
  call get_param(param_file, mod, "COLUMN_CHECK", column_check, &
                 "If true, add code to allow debugging of conservation \n"//&
                 "column-by-column.  This does not change answers, but \n"//&
                 "can increase model run time.", default=.false.)
  call get_param(param_file, mod, "IMBALANCE_TOLERANCE", imb_tol, &
                 "The tolerance for imbalances to be flagged by COLUMN_CHECK.", &
                 units="nondim", default=1.0e-9)
  call get_param(param_file, mod, "ICE_BOUNDS_CHECK", bounds_check, &
                 "If true, periodically check the values of ice and snow \n"//&
                 "temperatures and thicknesses to ensure that they are \n"//&
                 "sensible, and issue warnings if they are not.  This \n"//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)
  call get_param(param_file, mod, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call get_param(param_file, mod, "GLOBAL_INDEXING", global_indexing, &
                 "If true, use a global lateral indexing convention, so \n"//&
                 "that corresponding points on different processors have \n"//&
                 "the same index. This does not work with static memory.", &
                 default=.false., layoutParam=.true.)
#ifdef STATIC_MEMORY_
  if (global_indexing) call SIS_error(FATAL, "ice_model_init: "//&
       "GLOBAL_INDEXING can not be true with STATIC_MEMORY.")
#endif
  call get_param(param_file, mod, "FIRST_DIRECTION", first_direction, &
                 "An integer that indicates which direction goes first \n"//&
                 "in parts of the code that use directionally split \n"//&
                 "updates, with even numbers (or 0) used for x- first \n"//&
                 "and odd numbers used for y-first.", default=0)
  call log_param(param_file, mod, "! VERONA_COUPLER", Verona, &
                 "If true, carry out all of the sea ice calls so that SIS2 \n"//&
                 "will work with the Verona and earlier releases of the \n"//&
                 "FMS coupler code in configurations that use the exchange \n"//&
                 "grid to communicate with the atmosphere or land.", &
                 layoutParam=.true.)

  call get_param(param_file, mod, "ICE_SEES_ATMOS_WINDS", atmos_winds, &
                 "If true, the sea ice is being given wind stresses with \n"//&
                 "the atmospheric sign convention, and need to have their \n"//&
                 "sign changed.", default=.true.)
  call get_param(param_file, mod, "ICE_BULK_SALINITY", ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "g/kg", &
                 default=4.0, do_not_log=.true.)
  call get_param(param_file, mod, "ICE_RELATIVE_SALINITY", ice_rel_salin, &
                 "The initial salinity of sea ice as a fraction of the \n"//&
                 "salinity of the seawater from which it formed.", &
                 units = "nondim", default=0.0, do_not_log=.true.)
  if ((ice_bulk_salin < 0.0) .or. (ice_rel_salin > 0.0)) ice_bulk_salin = 0.0

  call get_param(param_file, mod, "APPLY_SLP_TO_OCEAN", slp2ocean, &
                 "If true, apply the atmospheric sea level pressure to \n"//&
                 "the ocean.", default=.false.)
  call get_param(param_file, mod, "MIN_H_FOR_TEMP_CALC", h_lo_lim, &
                 "The minimum ice thickness at which to do temperature \n"//&
                 "calculations.", units="m", default=0.0)
  call get_param(param_file, mod, "MIN_OCEAN_PARTSIZE", ocean_part_min, &
                 "The minimum value for the fractional open-ocean area. \n"//&
                 "This can be 0, but for some purposes it may be useful \n"//&
                 "to set this to a miniscule value (like 1e-40) that will \n"//&
                 "be lost to roundoff during any sums so that the open \n"//&
                 "ocean fluxes can be used in with new categories.", &
                 units="nondim", default=0.0)
  if ((ocean_part_min < 0.0) .or. (ocean_part_min > 1.0e-20)) &
    call SIS_error(FATAL, "MIN_OCEAN_PARTSIZE has been set outside of the valid"//&
                   "range of 0 to 1e-20.")
  call get_param(param_file, mod, "DO_ICEBERGS", do_icebergs, &
                 "If true, call the iceberg module.", default=.false.)
  if (do_icebergs) then
    call get_param(param_file, mod, "PASS_ICEBERG_AREA_TO_OCEAN", pass_iceberg_area_to_ocean, &
                 "If true, iceberg area is passed through coupler", default=.false.)
  else ; pass_iceberg_area_to_ocean = .false. ; endif

  call get_param(param_file, mod, "ADD_DIURNAL_SW", add_diurnal_sw, &
                 "If true, add a synthetic diurnal cycle to the shortwave \n"//&
                 "radiation.", default=.false.)
  call get_param(param_file, mod, "DO_SUN_ANGLE_FOR_ALB", do_sun_angle_for_alb, &
                 "If true, find the sun angle for calculating the ocean \n"//&
                 "albedo within the sea ice model.", default=.false.)
  call get_param(param_file, mod, "DO_RIDGING", do_ridging, &
                 "If true, call the ridging routines.", default=.false.)

  call get_param(param_file, mod, "RESTARTFILE", restart_file, &
                 "The name of the restart file.", default="ice_model.res.nc")
  if (fast_ice_PE.eqv.slow_ice_PE) then
    call get_param(param_file, mod, "FAST_ICE_RESTARTFILE", fast_rest_file, &
                   "The name of the restart file for those elements of the \n"//&
                   "the sea ice that are handled by the fast ice PEs.", &
                   default=restart_file)
  else
    call get_param(param_file, mod, "FAST_ICE_RESTARTFILE", fast_rest_file, &
                   "The name of the restart file for those elements of the \n"//&
                   "the sea ice that are handled by the fast ice PEs.", &
                   default="ice_model_fast.res.nc")
  endif

  call get_param(param_file, mod, "MASSLESS_ICE_ENTH", massless_ice_enth, &
                 "The ice enthalpy fill value for massless categories.", &
                 units="J kg-1", default=0.0, do_not_log=.true.)
  call get_param(param_file, mod, "MASSLESS_SNOW_ENTH", massless_snow_enth, &
                 "The snow enthalpy fill value for massless categories.", &
                 units="J kg-1", default=0.0, do_not_log=.true.)
  call get_param(param_file, mod, "MASSLESS_ICE_SALIN", massless_ice_salin, &
                 "The ice salinity fill value for massless categories.", &
                 units="g kg-1", default=0.0, do_not_log=.true.)
  call get_param(param_file, "MOM", "WRITE_GEOM", write_geom, &
                 "If =0, never write the geometry and vertical grid files.\n"//&
                 "If =1, write the geometry and vertical grid files only for\n"//&
                 "a new simulation. If =2, always write the geometry and\n"//&
                 "vertical grid files. Other values are invalid.", default=1)
  call get_param(param_file, "MOM", "INTERPOLATE_FLUXES", interp_fluxes, &
                 "If true, interpolate a linearized version of the fast \n"//&
                 "fluxes into arealess categories.", default=.true.)
  if (write_geom<0 .or. write_geom>2) call SIS_error(FATAL,"SIS2: "//&
         "WRITE_GEOM must be equal to 0, 1 or 2.")
  write_geom_files = ((write_geom==2) .or. ((write_geom==1) .and. &
     ((dirs%input_filename(1:1)=='n') .and. (LEN_TRIM(dirs%input_filename)==1))))

  nudge_sea_ice = .false. ; call read_param(param_file, "NUDGE_SEA_ICE", nudge_sea_ice)
  nCat_dflt = 5 ; if (slab_ice) nCat_dflt = 1
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
    Ice%sCS%Time = Time

    ! Set some pointers for convenience.
    sIST => Ice%sCS%IST ; sIG => Ice%sCS%IG
    sIST%Cgrid_dyn = Cgrid_dyn

    Ice%sCS%do_icebergs = do_icebergs
    Ice%sCS%pass_iceberg_area_to_ocean = pass_iceberg_area_to_ocean
    Ice%sCS%specified_ice = specified_ice
    Ice%sCS%Cgrid_dyn = Cgrid_dyn
    Ice%sCS%bounds_check = bounds_check
    Ice%sCS%debug = debug

    ! Set up the ice-specific grid describing categories and ice layers.
    call set_ice_grid(sIG, param_file, nCat_dflt)
    if (slab_ice) sIG%CatIce = 1 ! open water and ice ... but never in same place
    CatIce = sIG%CatIce ; NkIce = sIG%NkIce
    call initialize_ice_categories(sIG, Rho_ice, param_file)


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

    ! Set the bathymetry, Coriolis parameter, open channel widths and masks.
    call SIS_initialize_fixed(dG, param_file, write_geom_files, dirs%output_directory)

    call set_hor_grid(sG, param_file, global_indexing=global_indexing)
    call copy_dyngrid_to_SIS_horgrid(dG, sG)
    call destroy_dyn_horgrid(dG)

  ! Allocate and register fields for restarts.

    call set_domain(sGD%mpp_domain)
    if (.not.associated(Ice%Ice_restart)) allocate(Ice%Ice_restart)

    call ice_type_slow_reg_restarts(sGD%mpp_domain, CatIce, &
                      param_file, Ice, Ice%Ice_restart, restart_file)

    call alloc_IST_arrays(sHI, sIG, sIST, omit_tsurf=Eulerian_tsurf)
    call ice_state_register_restarts(sGD%mpp_domain, sIST, sIG, Ice%Ice_restart, restart_file)

    call alloc_ocean_sfc_state(Ice%sCS%OSS, sHI, sIST%Cgrid_dyn)
    Ice%sCS%OSS%kmelt = kmelt

    call alloc_simple_OSS(Ice%sCS%sOSS, sHI)

    call alloc_ice_ocean_flux(Ice%sCS%IOF, sHI, do_iceberg_fields=Ice%sCS%do_icebergs)
    Ice%sCS%IOF%slp2ocean = slp2ocean
    Ice%sCS%IOF%flux_uv_stagger = Ice%flux_uv_stagger
    call alloc_fast_ice_avg(Ice%sCS%FIA, sHI, sIG, interp_fluxes)

    call SIS_dyn_trans_register_restarts(sGD%mpp_domain, sHI, sIG, param_file,&
                                Ice%sCS%dyn_trans_CSp, Ice%Ice_restart, restart_file)

    call SIS_diag_mediator_init(sG, sIG, param_file, Ice%sCS%diag, component="SIS", &
                                doc_file_dir = dirs%output_directory)
    call set_SIS_axes_info(sG, sIG, param_file, Ice%sCS%diag)

    call ice_thermo_init(param_file, sIST%ITV, init_EOS=nudge_sea_ice)
    call get_SIS2_thermo_coefs(sIST%ITV, enthalpy_units=enth_unit)

    ! Register tracers that will be advected around.
    call register_SIS_tracer_pair(sIST%enth_ice, NkIce, "enth_ice", &
                                  sIST%enth_snow, 1, "enth_snow", &
                                  sG, sIG, param_file, sIST%TrReg, &
                                  massless_iceval=massless_ice_enth*enth_unit, &
                                  massless_snowval=massless_snow_enth*enth_unit)

    if (ice_rel_salin > 0.0) then
      call register_SIS_tracer(sIST%sal_ice, sG, sIG, NkIce, "salin_ice", param_file, &
                               sIST%TrReg, snow_tracer=.false., &
                               massless_val=massless_ice_salin, nonnegative=.true.)
    endif

  !   Register any tracers that will be handled via tracer flow control for
  ! restarts and advection.
    call SIS_call_tracer_register(sG, sIG, param_file, Ice%sCS%SIS_tracer_flow_CSp, &
                                  Ice%sCS%diag, sIST%TrReg, Ice%Ice_restart, restart_file)

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
      Ice%area(i2,j2) = sG%areaT(i,j) * sG%mask2dT(i,j)
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

    if (single_IST) then
      Ice%fCS%IST => Ice%sCS%IST
      Ice%fCS%G => Ice%sCS%G
      fG => Ice%fCS%G
      fGD => Ice%fCS%G%Domain
      fHI = sHI
      Ice%fCS%FIA => Ice%sCS%FIA
      Ice%fCS%sOSS => Ice%sCS%sOSS
    else
      ! Set up the domains and lateral grids.
      if (.not.associated(Ice%fCS%IST)) allocate(Ice%fCS%IST)
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
      call SIS_initialize_fixed(dG, param_file, .false., dirs%output_directory)

      call set_hor_grid(Ice%fCS%G, param_file, global_indexing=global_indexing)
      call copy_dyngrid_to_SIS_horgrid(dG, Ice%fCS%G)
      call destroy_dyn_horgrid(dG)
    endif

    Ice%fCS%bounds_check = bounds_check
    Ice%fCS%debug = debug
    Ice%fCS%Eulerian_tsurf = Eulerian_tsurf

    ! Set up the ice-specific grid describing categories and ice layers.
    call set_ice_grid(Ice%fCS%IG, param_file, nCat_dflt)
    if (slab_ice) Ice%fCS%IG%CatIce = 1 ! open water and ice ... but never in same place
    CatIce = Ice%fCS%IG%CatIce ; NkIce = Ice%fCS%IG%NkIce

    call initialize_ice_categories(Ice%fCS%IG, Rho_ice, param_file)

  ! Allocate and register fields for restarts.

    if (.not.slow_ice_PE) call set_domain(fGD%mpp_domain)
    if (split_restart_files) then
      if (.not.associated(Ice%Ice_fast_restart)) allocate(Ice%Ice_fast_restart)
    else
      Ice%Ice_fast_restart => Ice%Ice_restart
    endif

  ! These allocation routines are called on all PEs; whether or not the variables
  ! they allocate are registered for inclusion in restart files is determined by
  ! whether the Ice%Ice...restart types are associated.
    call ice_type_fast_reg_restarts(fGD%mpp_domain, CatIce, &
                      param_file, Ice, Ice%Ice_fast_restart, fast_rest_file)

    if (.not.single_IST) then
      call alloc_IST_arrays(fHI, Ice%fCS%IG, Ice%fCS%IST, &
                            omit_velocities=.true., omit_tsurf=Eulerian_tsurf)

      call alloc_fast_ice_avg(Ice%fCS%FIA, fHI, Ice%fCS%IG, interp_fluxes)

      call alloc_simple_OSS(Ice%fCS%sOSS, fHI)
    endif
    call alloc_total_sfc_flux(Ice%fCS%TSF, fHI)
    Ice%fCS%FIA%atmos_winds = atmos_winds

    call ice_rad_register_restarts(fGD%mpp_domain, fHI, Ice%fCS%IG, param_file, &
                                   Ice%fCS%Rad, Ice%Ice_fast_restart, fast_rest_file)
    Ice%fCS%Rad%do_sun_angle_for_alb = do_sun_angle_for_alb
    Ice%fCS%Rad%add_diurnal_sw = add_diurnal_sw
    Ice%fCS%Rad%frequent_albedo_update = .true.
    !### Instead perhaps this could be
    !###   Ice%fCS%Rad%frequent_albedo_update = Ice%fCS%Rad%do_sun_angle_for_alb .or. (Time_step_slow > dT_Rad)
    !### However this changes answers in coupled models.  I don't understand why. -RWH

    allocate(Ice%fCS%diag)
    call SIS_diag_mediator_init(fG, Ice%fCS%IG, param_file, Ice%fCS%diag, component="SIS_fast", &
                                doc_file_dir = dirs%output_directory)
    call set_SIS_axes_info(fG, Ice%fCS%IG, param_file, Ice%fCS%diag, axes_set_name="ice_fast")

    if (.not.single_IST) then
      call ice_thermo_init(param_file, Ice%fCS%IST%ITV, init_EOS=nudge_sea_ice)

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


  ! Read the restart file, if it exists, and initialize the ice arrays to
  ! to default values if it does not.
  if (slow_ice_PE) then
    ! Set some pointers for convenience.
    sIST => Ice%sCS%IST ; sIG => Ice%sCS%IG ; sG => Ice%sCS%G

    allocate(S_col(NkIce)) ; S_col(:) = 0.0
    call get_SIS2_thermo_coefs(sIST%ITV, ice_salinity=S_col, enthalpy_units=enth_unit, &
                               specified_thermo_salinity=spec_thermo_sal)

    restart_path = trim(dirs%restart_input_dir)//trim(restart_file)

    if (file_exist(restart_path)) then
      ! Set values of IG%H_to_kg_m2 that will permit its absence from the restart
      ! file to be detected, and its difference from the value in this run to
      ! be corrected for.
      H_to_kg_m2_tmp = sIG%H_to_kg_m2
      sIG%H_to_kg_m2 = -1.0
      is_restart = .true.

      call restore_state(Ice%Ice_restart, directory=dirs%restart_input_dir)

      ! Approximately initialize state fields that are not present
      ! in SIS1 restart files.  This is obsolete and can probably be eliminated.

      ! Initialize the ice salinity.
      if (.not.query_initialized(Ice%Ice_restart, 'sal_ice')) then
        allocate(sal_ice_tmp(sG%isd:sG%ied, sG%jsd:sG%jed, CatIce, NkIce)) ; sal_ice_tmp(:,:,:,:) = 0.0
        do n=1,NkIce
          write(nstr, '(I4)') n ; nstr = adjustl(nstr)
          id_sal = register_restart_field(Ice%Ice_restart, restart_file, 'sal_ice'//trim(nstr), &
                                       sal_ice_tmp(:,:,:,n), domain=sGD%mpp_domain, &
                                       mandatory=.false., read_only=.true.)
          call restore_state(Ice%Ice_restart, id_sal, directory=dirs%restart_input_dir)
        enddo

        if (query_initialized(Ice%Ice_restart, 'sal_ice1')) then
          do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
            sIST%sal_ice(i,j,k,1) = sal_ice_tmp(i,j,k,1)
          enddo ; enddo ; enddo
        else
          sIST%sal_ice(:,:,:,1) = ice_bulk_salin
        endif
        do n=2,NkIce
          write(nstr, '(I4)') n ; nstr = adjustl(nstr)
          if (query_initialized(Ice%Ice_restart, 'sal_ice'//trim(nstr))) then
            do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
              sIST%sal_ice(i,j,k,n) = sal_ice_tmp(i,j,k,n)
            enddo ; enddo ; enddo
          else
            sIST%sal_ice(:,:,:,n) = sIST%sal_ice(:,:,:,n-1)
          endif
        enddo

        deallocate(sal_ice_tmp)
      endif

      read_aux_restart = (.not.query_initialized(Ice%Ice_restart, 'enth_ice')) .or. &
                         (.not.query_initialized(Ice%Ice_restart, 'enth_snow'))
      if (read_aux_restart) then
        allocate(t_snow_tmp(sG%isd:sG%ied, sG%jsd:sG%jed, CatIce)) ; t_snow_tmp(:,:,:) = 0.0
        allocate(t_ice_tmp(sG%isd:sG%ied, sG%jsd:sG%jed, CatIce, NkIce)) ; t_ice_tmp(:,:,:,:) = 0.0

        idr = register_restart_field(Ice%Ice_restart, restart_file, 't_snow', t_snow_tmp, &
                                     domain=sGD%mpp_domain, mandatory=.false., read_only=.true.)
        call restore_state(Ice%Ice_restart, idr, directory=dirs%restart_input_dir)
        do n=1,NkIce
          write(nstr, '(I4)') n ; nstr = adjustl(nstr)
          idr = register_restart_field(Ice%Ice_restart, restart_file, 't_ice'//trim(nstr), &
                                       t_ice_tmp(:,:,:,n), domain=sGD%mpp_domain, &
                                       mandatory=.false., read_only=.true.)
          call restore_state(Ice%Ice_restart, idr, directory=dirs%restart_input_dir)
        enddo
      endif

      ! Initialize the ice enthalpy.
      if (.not.query_initialized(Ice%Ice_restart, 'enth_ice')) then
        if (.not.query_initialized(Ice%Ice_restart, 't_ice1')) then
          call SIS_error(FATAL, "Either t_ice1 or enth_ice must be present in the SIS2 restart file "//restart_path)
        endif

        if (spec_thermo_sal) then
          do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
            sIST%enth_ice(i,j,k,1) = Enth_from_TS(t_ice_tmp(i,j,k,1), S_col(1), sIST%ITV)
          enddo ; enddo ; enddo
        else
          do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
            sIST%enth_ice(i,j,k,1) = Enth_from_TS(t_ice_tmp(i,j,k,1), &
                     sIST%sal_ice(i,j,k,1), sIST%ITV)
          enddo ; enddo ; enddo
        endif

        do n=2,NkIce
          write(nstr, '(I4)') n ; nstr = adjustl(nstr)
          if (.not.query_initialized(Ice%Ice_restart, 't_ice'//trim(nstr))) &
            t_ice_tmp(:,:,:,n) = t_ice_tmp(:,:,:,n-1)

          if (spec_thermo_sal) then
            do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
              sIST%enth_ice(i,j,k,n) = Enth_from_TS(t_ice_tmp(i,j,k,n), &
                               S_col(n), sIST%ITV)
            enddo ; enddo ; enddo
          else
            do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
              sIST%enth_ice(i,j,k,n) = Enth_from_TS(t_ice_tmp(i,j,k,n), &
                               sIST%sal_ice(i,j,k,n), sIST%ITV)
            enddo ; enddo ; enddo
          endif
        enddo
      endif

      ! Initialize the snow enthalpy.
      if (.not.query_initialized(Ice%Ice_restart, 'enth_snow')) then
        if (.not.query_initialized(Ice%Ice_restart, 't_snow')) then
          if (query_initialized(Ice%Ice_restart, 't_ice1')) then
            t_snow_tmp(:,:,:) = t_ice_tmp(:,:,:,1)
          else
            do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
              t_snow_tmp(i,j,k) = Temp_from_En_S(sIST%enth_ice(i,j,k,1), &
                                    sIST%sal_ice(i,j,k,1), sIST%ITV)
            enddo ; enddo ; enddo
          endif
        endif
        do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          sIST%enth_snow(i,j,k,1) = Enth_from_TS(t_snow_tmp(i,j,k), 0.0, sIST%ITV)
        enddo ; enddo ; enddo
      endif

      if (read_aux_restart) deallocate(t_snow_tmp, t_ice_tmp)

      if (allocated(sIST%t_surf) .and. &
          .not.query_initialized(Ice%Ice_restart, 't_surf_ice')) then
        sIST%t_surf(:,:,:) = T_0degC
      endif

      H_rescale_ice = 1.0 ; H_rescale_snow = 1.0
      if (sIG%H_to_kg_m2 == -1.0) then
        ! This is an older restart file, and the snow and ice thicknesses are in
        ! m, and not a mass coordinate.
        H_rescale_ice = Rho_ice / H_to_kg_m2_tmp
        H_rescale_snow = Rho_snow / H_to_kg_m2_tmp
      elseif (sIG%H_to_kg_m2 /= H_to_kg_m2_tmp) then
        H_rescale_ice = sIG%H_to_kg_m2 / H_to_kg_m2_tmp
        H_rescale_snow = H_rescale_ice
      endif
      sIG%H_to_kg_m2 = H_to_kg_m2_tmp

      ! Deal with any ice masses or thicknesses over land, and rescale to
      ! account for differences between the current thickness units and whatever
      ! thickness units were in the input restart file.
      do k=1,CatIce
        sIST%mH_snow(:,:,k) = sIST%mH_snow(:,:,k) * H_rescale_snow * sG%mask2dT(:,:)
        sIST%mH_ice(:,:,k) = sIST%mH_ice(:,:,k) * H_rescale_ice * sG%mask2dT(:,:)
      enddo

      if (ocean_part_min > 0.0) then ; do j=jsc,jec ; do i=isc,iec
        sIST%part_size(i,j,0) = max(sIST%part_size(i,j,0), ocean_part_min)
      enddo ; enddo ; endif

      !--- update the halo values.
      call pass_var(sIST%part_size, sGD)
      call pass_var(sIST%mH_ice, sGD, complete=.false.)
      call pass_var(sIST%mH_snow, sGD, complete=.false.)
      do l=1,NkIce
        call pass_var(sIST%enth_ice(:,:,:,l), sGD, complete=.false.)
      enddo
      call pass_var(sIST%enth_snow(:,:,:,1), sGD, complete=.true.)

      if (Cgrid_dyn) then
        call pass_vector(sIST%u_ice_C, sIST%v_ice_C, sGD, stagger=CGRID_NE)
      else
        call pass_vector(sIST%u_ice_B, sIST%v_ice_B, sGD, stagger=BGRID_NE)
      endif

      if (fast_ice_PE .and. .not.split_restart_files) then
        init_coszen = .not.query_initialized(Ice%Ice_fast_restart, 'coszen')
        init_Tskin  = .not.query_initialized(Ice%Ice_fast_restart, 'T_skin')
        init_rough  = .not.(query_initialized(Ice%Ice_fast_restart, 'rough_mom') .and. &
                            query_initialized(Ice%Ice_fast_restart, 'rough_heat') .and. &
                            query_initialized(Ice%Ice_fast_restart, 'rough_moist'))
      endif

    else ! no restart file implies initialization with no ice
      sIST%part_size(:,:,:) = 0.0
      sIST%part_size(:,:,0) = 1.0

      if (allocated(sIST%t_surf)) sIST%t_surf(:,:,:) = T_0degC
      sIST%sal_ice(:,:,:,:) = ice_bulk_salin

      enth_spec_snow = Enth_from_TS(0.0, 0.0, sIST%ITV)
      sIST%enth_snow(:,:,:,1) = enth_spec_snow
      do n=1,NkIce
        enth_spec_ice = Enth_from_TS(0.0, S_col(n), sIST%ITV)
        sIST%enth_ice(:,:,:,n) = enth_spec_ice
      enddo

      allocate(h_ice_input(sG%isc:sG%iec,sG%jsc:sG%jec))
      call get_sea_surface(Ice%sCS%Time, Ice%sCS%OSS%SST_C(isc:iec,jsc:jec), &
                           sIST%part_size(isc:iec,jsc:jec,0:1), &
                           h_ice_input, ice_domain=Ice%slow_domain_NH, ts_in_K=.false. )
      do j=jsc,jec ; do i=isc,iec
        sIST%mH_ice(i,j,1) = h_ice_input(i,j)*(Rho_ice*sIG%kg_m2_to_H)
      enddo ; enddo

      !   Transfer ice to the correct thickness category.  If do_ridging=.false.,
      ! the first call to ice_redistribute has the same result.  At present, all
      ! tracers are initialized to their default values, and snow is set to 0,
      ! and so do not need to be updated here.
      if (do_ridging) then
        do j=jsc,jec ; do i=isc,iec ; if (sIST%mH_ice(i,j,1) > sIG%mH_cat_bound(1)) then
          do k=CatIce,2,-1 ; if (sIST%mH_ice(i,j,1) > sIG%mH_cat_bound(k-1)) then
            sIST%part_size(i,j,k) = sIST%part_size(i,j,1)
            sIST%part_size(i,j,1) = 0.0
            sIST%mH_ice(i,j,k) = sIST%mH_ice(i,j,1) ; sIST%mH_ice(i,j,1) = 0.0
            !  sIST%mH_snow(i,j,k) = sIST%mH_snow(i,j,1) ; sIST%mH_snow(i,j,1) = 0.0
            exit ! from k-loop
          endif ; enddo
        endif ; enddo ; enddo
      endif

      if (ocean_part_min > 0.0) then ; do j=jsc,jec ; do i=isc,iec
        sIST%part_size(i,j,0) = max(sIST%part_size(i,j,0), ocean_part_min)
      enddo ; enddo ; endif

      deallocate(h_ice_input)

      call pass_var(sIST%part_size, sGD, complete=.true. )
      call pass_var(sIST%mH_ice, sGD, complete=.true. )

      init_coszen = .true. ; init_Tskin = .true. ; init_rough = .true.

    endif ! file_exist(restart_path)

    deallocate(S_col)

  ! The restart files have now been read or the variables that would have been
  ! in the restart files have been initialized.  Now call the initialization
  ! routines for any dependent sub-modules.

    call ice_diagnostics_init(Ice%sCS%IOF, Ice%sCS%OSS, Ice%sCS%FIA, sG, sIG, &
                              Ice%sCS%diag, Ice%sCS%Time, Cgrid=sIST%Cgrid_dyn)
    Ice%axes(1:2) = Ice%sCS%diag%axesTc%handles(1:2)

    Ice%sCS%Time_step_slow = Time_step_slow

    call SIS_slow_thermo_init(Ice%sCS%Time, sG, sIG, param_file, Ice%sCS%diag, &
                              Ice%sCS%slow_thermo_CSp, Ice%sCS%SIS_tracer_flow_CSp)

    call SIS_dyn_trans_init(Ice%sCS%Time, sG, sIG, param_file, Ice%sCS%diag, &
                            Ice%sCS%dyn_trans_CSp, dirs%output_directory, Time_Init)

    call SIS_slow_thermo_set_ptrs(Ice%sCS%slow_thermo_CSp, &
             transport_CSp=SIS_dyn_trans_transport_CS(Ice%sCS%dyn_trans_CSp), &
             sum_out_CSp=SIS_dyn_trans_sum_output_CS(Ice%sCS%dyn_trans_CSp))

  !   Initialize any tracers that will be handled via tracer flow control.
    call SIS_tracer_flow_control_init(Ice%sCS%Time, sG, sIG, param_file, &
                                      Ice%sCS%SIS_tracer_flow_CSp, is_restart)

  ! Initialize icebergs
    if (Ice%sCS%do_icebergs) then
      call get_param(param_file, mod, "ICEBERG_WINDSTRESS_BUG", Ice%sCS%berg_windstress_bug, &
                 "If true, use older code that applied an old ice-ocean \n"//&
                 "stress to the icebergs in place of the current air-ocean \n"//&
                 "stress.  This option is here for backward compatibility, \n"//&
                 "but should be avoided.", default=.false.)

      isc = sG%isc ; iec = sG%iec ; jsc = sG%jsc ; jec = sG%jec

      if (ASSOCIATED(sGD%maskmap)) then
        call icebergs_init(Ice%icebergs, sGD%niglobal, sGD%njglobal, &
                sGD%layout, sGD%io_layout, Ice%axes(1:2), &
                sGD%X_flags, sGD%Y_flags, time_type_to_real(Time_step_slow), &
                Time, sG%geoLonBu(isc:iec,jsc:jec), sG%geoLatBu(isc:iec,jsc:jec), &
                sG%mask2dT(isc-1:iec+1,jsc-1:jec+1), &
                sG%dxCv(isc-1:iec+1,jsc-1:jec+1), sG%dyCu(isc-1:iec+1,jsc-1:jec+1), &
                Ice%area,  sG%cos_rot(isc-1:iec+1,jsc-1:jec+1), &
                sG%sin_rot(isc-1:iec+1,jsc-1:jec+1), maskmap=sGD%maskmap )
      else
        call icebergs_init(Ice%icebergs, sGD%niglobal, sGD%njglobal, &
                 sGD%layout, sGD%io_layout, Ice%axes(1:2), &
                 sGD%X_flags, sGD%Y_flags, time_type_to_real(Time_step_slow), &
                 Time, sG%geoLonBu(isc:iec,jsc:jec), sG%geoLatBu(isc:iec,jsc:jec), &
                 sG%mask2dT(isc-1:iec+1,jsc-1:jec+1), &
                 sG%dxCv(isc-1:iec+1,jsc-1:jec+1), sG%dyCu(isc-1:iec+1,jsc-1:jec+1), &
                 Ice%area, sG%cos_rot(isc-1:iec+1,jsc-1:jec+1), &
                 sG%sin_rot(isc-1:iec+1,jsc-1:jec+1) )
      endif
    endif

    if (Verona) then
      !   The Verona and earlier versions of the coupler code make calls to set
      ! up the exchange grid right at the start of the coupled timestep, before
      ! information about the part_size distribution can be copied from the slow
      ! processors to the fast processors.  This will cause coupled models with

      if (fast_ice_PE) then
        write_error_mesg = .not.((sHI%iec-sHI%isc==fHI%iec-fHI%isc) .and. &
                                 (sHI%jec-sHI%jsc==fHI%jec-fHI%jsc))
      else ; write_error_mesg = .true.
      endif

      if (write_error_mesg) call SIS_error(FATAL, &
          "The Verona coupler will not work unless the fast and slow portions "//&
          "of SIS2 use the same PEs and layout.")

      ! Set the computational domain sizes using the ice model's indexing convention.
      isc = sHI%isc ; iec = sHI%iec ; jsc = sHI%jsc ; jec = sHI%jec
      i_off = LBOUND(Ice%part_size,1) - sHI%isc ; j_off = LBOUND(Ice%part_size,2) - sHI%jsc
      do k=0,CatIce ; do j=jsc,jec ; do i=isc,iec
        i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
        Ice%part_size(i2,j2,k2) = sIST%part_size(i,j,k)
      enddo ; enddo ; enddo

    endif

    ! Do any error checking here.
    if (debug) call ice_grid_chksum(sG, haloshift=1)

    call write_ice_statistics(sIST, Ice%sCS%Time, 0, sG, sIG, &
                   SIS_dyn_trans_sum_output_CS(Ice%sCS%dyn_trans_CSp))
  endif  ! slow_ice_PE


  if (fast_ice_PE) then
    ! Read the fast_restart file and initialize the subsidiary modules of the
    ! fast ice processes.

    ! Set some pointers for convenience.
    fG => Ice%fCS%G ; fGD => Ice%fCS%G%Domain

    if ((.not.slow_ice_PE) .or. split_restart_files) then
      ! Read the fast restart file, if it exists.
      fast_rest_path = trim(dirs%restart_input_dir)//trim(fast_rest_file)
      if (file_exist(fast_rest_path)) then
        call restore_state(Ice%Ice_fast_restart, directory=dirs%restart_input_dir)
        init_coszen = .not.query_initialized(Ice%Ice_fast_restart, 'coszen')
        init_Tskin = .not.query_initialized(Ice%Ice_fast_restart, 'T_skin')
        init_rough  = .not.(query_initialized(Ice%Ice_fast_restart, 'rough_mom') .and. &
                            query_initialized(Ice%Ice_fast_restart, 'rough_heat') .and. &
                            query_initialized(Ice%Ice_fast_restart, 'rough_moist'))
      else
        init_coszen = .true. ; init_Tskin = .true. ; init_rough = .true.
      endif
    endif

!  if (Ice%fCS%Rad%add_diurnal_sw .or. Ice%fCS%Rad%do_sun_angle_for_alb) then
!    call set_domain(fGD%mpp_domain)
    call astronomy_init
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
    call SIS_optics_init(param_file, Ice%fCS%optics_CSp)

    Ice%fCS%Time_step_fast = Time_step_fast
    Ice%fCS%Time_step_slow = Time_step_slow

    isc = fHI%isc ; iec = fHI%iec ; jsc = fHI%jsc ; jec = fHI%jec
    i_off = LBOUND(Ice%ocean_pt,1) - fHI%isc ; j_off = LBOUND(Ice%ocean_pt,2) - fHI%jsc
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      Ice%ocean_pt(i2,j2) = ( fG%mask2dT(i,j) > 0.5 )
    enddo ; enddo
    if (.not.slow_ice_PE) then
      Ice%axes(1:2) = Ice%fCS%diag%axesTc%handles(1:2)
    endif
  endif ! fast_ice_PE

  !nullify_domain perhaps could be called somewhere closer to set_domain
  !but it should be called after restore_state() otherwise it causes a restart mismatch
  call nullify_domain()

  call close_param_file(param_file)

  ! In the post-Verona coupler, share_ice_domains is called by the coupler
  ! after it switches the current PE_list to the one with all ice PEs.
  if (Verona) call share_ice_domains(Ice)

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

  if (Ice%shared_slow_fast_PEs) then
    iceClock = mpp_clock_id( 'Ice', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    ice_clock_fast = mpp_clock_id('Ice Fast', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    ice_clock_slow = mpp_clock_id('Ice Slow', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  else
    iceClock = 0 ! The comprehensive ice clock can not be used with separate fast and slow ice PEs.
    if (fast_ice_PE) then
      ice_clock_fast = mpp_clock_id('Ice Fast', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    elseif (slow_ice_PE) then
      ice_clock_slow = mpp_clock_id('Ice Slow', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    endif
  endif

  call callTree_leave("ice_model_init()")

end subroutine ice_model_init

subroutine share_ice_domains(Ice)
  type(ice_data_type), intent(inout) :: Ice

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
  call mpp_broadcast_domain(Ice%Domain)
  call mpp_broadcast_domain(Ice%slow_domain_NH)
  call mpp_broadcast_domain(Ice%slow_domain)
  call mpp_broadcast_domain(Ice%fast_domain)

  if (Ice%shared_slow_fast_PEs) then
    ice_clock_exchange = mpp_clock_id('Ice Fast/Slow Exchange', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  else
    ice_clock_exchange = mpp_clock_id('Ice Fast/Slow Exchange', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  endif

end subroutine share_ice_domains

!> initialize_ice_categories sets the bounds of the ice thickness categories.
subroutine initialize_ice_categories(IG, Rho_ice, param_file, hLim_vals)
  type(ice_grid_type),          intent(inout) :: IG
  real,                         intent(in)    :: Rho_ice
  type(param_file_type),        intent(in)    :: param_file
  real, dimension(:), optional, intent(in)    :: hLim_vals

  ! Initialize IG%cat_thick_lim and IG%mH_cat_bound here.
  !  ###This subroutine should be extended to add more options.

  real :: hlim_dflt(8) = (/ 1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! lower thickness limits 1...CatIce
  integer :: k, CatIce, list_size

  CatIce = IG%CatIce
  list_size = -1
  if (present(hLim_vals)) then ; if (size(hLim_vals(:)) > 1) then
    list_size = size(hlim_vals(:))
    do k=1,min(CatIce+1,list_size) ; IG%cat_thick_lim(k) = hlim_vals(k) ; enddo
  endif ; endif
  if (list_size < 2) then  ! Use the default categories.
    list_size = size(hlim_dflt(:))
    do k=1,min(CatIce+1,list_size) ; IG%cat_thick_lim(k) = hlim_dflt(k) ; enddo
  endif

  if ((CatIce+1 > list_size) .and. (list_size > 1)) then
    do k=list_size+1, CatIce+1
      IG%cat_thick_lim(k) =  2.0*IG%cat_thick_lim(k-1) - IG%cat_thick_lim(k-2)
    enddo
  endif

  do k=1,IG%CatIce+1
    IG%mH_cat_bound(k) = IG%cat_thick_lim(k) * (Rho_ice*IG%kg_m2_to_H)
  enddo
end subroutine initialize_ice_categories

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_end - writes the restart file and deallocates memory               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_end (Ice)
  type(ice_data_type), intent(inout) :: Ice

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
    call SIS_dyn_trans_end(Ice%sCS%dyn_trans_CSp)

    call SIS_slow_thermo_end(Ice%sCS%slow_thermo_CSp)

    call ice_thermo_end(Ice%sCS%IST%ITV)

    ! End icebergs
    if (Ice%sCS%do_icebergs) call icebergs_end(Ice%icebergs)

    call SIS_tracer_flow_control_end(Ice%sCS%SIS_tracer_flow_CSp)

    call dealloc_ice_ocean_flux(Ice%sCS%IOF)

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
