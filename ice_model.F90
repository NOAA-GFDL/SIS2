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

use SIS_diag_mediator, only : set_SIS_axes_info, SIS_diag_mediator_init, SIS_diag_mediator_end
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
! use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
! use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_get_input, only : Get_SIS_input, directories
use SIS_sum_output, only : SIS_sum_output_init,  write_ice_statistics
use SIS_transcribe_grid, only : copy_dyngrid_to_SIS_horgrid, copy_SIS_horgrid_to_dyngrid

use MOM_checksums,     only : chksum, uchksum, vchksum, Bchksum
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
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_ROUTINE

use astronomy_mod, only : astronomy_init, astronomy_end
use astronomy_mod, only : universal_time, orbital_time, diurnal_solar, daily_mean_solar
use coupler_types_mod, only : coupler_3d_bc_type
use ocean_albedo_mod, only : compute_ocean_albedo            ! ice sets ocean surface
use ocean_rough_mod,  only : compute_ocean_roughness         ! properties over water

use ice_type_mod, only : ice_data_type, dealloc_ice_arrays, ice_data_type_register_restarts
use ice_type_mod, only : Ice_public_type_chksum, Ice_public_type_bounds_check
use ice_type_mod, only : ice_model_restart, ice_stock_pe
use ice_type_mod, only : ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
use ice_type_mod, only : ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
use ice_type_mod, only : lnd_ice_bnd_type_chksum, ice_data_type_chksum
use SIS_types, only : ice_ocean_flux_type, alloc_ice_ocean_flux, dealloc_ice_ocean_flux
use SIS_types, only : ocean_sfc_state_type, alloc_ocean_sfc_state, dealloc_ocean_sfc_state
use SIS_types, only : fast_ice_avg_type, alloc_fast_ice_avg, dealloc_fast_ice_avg
use SIS_types, only : ice_rad_type, ice_rad_register_restarts, dealloc_ice_rad
use SIS_types, only : ice_state_type, ice_state_register_restarts, dealloc_IST_arrays
use SIS_ctrl_types, only : ice_diagnostics_init, SIS_slow_CS, SIS_fast_CS
use SIS_types, only : IST_chksum, IST_bounds_check
use ice_utils_mod, only : post_avg, ice_grid_chksum
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
use SIS_fast_thermo, only : do_update_ice_model_fast, avg_top_quantities
use SIS_fast_thermo, only : SIS_fast_thermo_init, SIS_fast_thermo_end

use SIS2_ice_thm,  only : ice_temp_SIS2, ice_optics_SIS2, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,  only : get_SIS2_thermo_coefs, enth_from_TS, Temp_from_En_S, T_freeze
use ice_bergs,     only : icebergs, icebergs_run, icebergs_init, icebergs_end

implicit none ; private

#include <SIS2_memory.h>

public :: ice_data_type, ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
public :: ice_model_init, ice_model_end, update_ice_model_fast, ice_stock_pe
public :: update_ice_model_slow_up, update_ice_model_slow_dn
public :: ice_model_restart  ! for intermediate restarts
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum, ice_data_type_chksum

integer :: iceClock, iceClock1, iceCLock2, iceCLock3

contains

!-----------------------------------------------------------------------
!
! Coupler interface to do slow ice processes:  dynamics, transport, mass
!
subroutine update_ice_model_slow_dn ( Atmos_boundary, Land_boundary, Ice )
  type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
  type(land_ice_boundary_type),  intent(inout) :: Land_boundary
  type(ice_data_type),           intent(inout) :: Ice

  real :: dt_slow  ! The time step over which to advance the model.

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(iceClock2)
  dt_slow = time_type_to_real(Ice%Ice_state%Time_step_slow)

  ! average fluxes from update_ice_model_fast
  call avg_top_quantities(Ice%FIA, Ice%Rad, Ice%Ice_state%part_size, Ice%G, Ice%IG)

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Start update_ice_model_slow_dn", Ice)
  endif

  call set_ice_ocean_fluxes(Ice%IOF, Ice, Land_boundary, Ice%G, Ice%IG)

  if (Ice%sCS%do_icebergs) then
    call mpp_clock_end(iceClock2) ; call mpp_clock_end(iceClock)
    call update_icebergs(Ice%Ice_state, Ice%OSS, Ice%IOF, Ice%FIA, Ice%icebergs, &
                         dt_slow, Ice%G, Ice%IG, Ice%dyn_trans_CSp)
    call mpp_clock_begin(iceClock) ; call mpp_clock_begin(iceClock2)
  endif

  call slow_thermodynamics(Ice%Ice_state, dt_slow, Ice%slow_thermo_CSp, &
                           Ice%OSS, Ice%FIA, Ice%IOF, Ice%G, Ice%IG)

  call SIS_dynamics_trans(Ice%Ice_state, Ice%OSS, Ice%FIA, Ice%IOF, &
                          dt_slow, Ice%dyn_trans_CSp, Ice%icebergs, Ice%G, Ice%IG)

  if (Ice%sCS%debug) &
    call IST_chksum("Before set_ocean_top_fluxes", Ice%Ice_state, Ice%G, Ice%IG)
  ! Set up the thermodynamic fluxes in the externally visible structure Ice.
  call set_ocean_top_fluxes(Ice, Ice%Ice_state, Ice%IOF, Ice%G, Ice%IG, Ice%sCS)

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("End update_ice_model_slow_dn", Ice)
  endif
  if (Ice%sCS%bounds_check) then
    call Ice_public_type_bounds_check(Ice, Ice%G, "End update_ice_slow")
  endif

  call mpp_clock_end(iceClock2) ; call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow_dn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ice_ocean_fluxes copies the ice surface fluxes and any other fields into
!! the ice_ocean_flux_type.
subroutine set_ice_ocean_fluxes(IOF, Ice, LIB, G, IG)
  type(ice_ocean_flux_type),    intent(inout) :: IOF
  type(ice_data_type),          intent(in)    :: Ice
  type(land_ice_boundary_type), intent(in)    :: LIB
  type(SIS_hor_grid_type),      intent(in)    :: G
  type(ice_grid_type),          intent(in)    :: IG

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  ! Store liquid runoff and other fluxes from the land to the ice or ocean.
  i_off = LBOUND(LIB%runoff,1) - G%isc ; j_off = LBOUND(LIB%runoff,2) - G%jsc
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IOF,LIB,i_off,j_off) &
!$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    IOF%runoff(i,j)  = LIB%runoff(i2,j2)
    IOF%calving(i,j) = LIB%calving(i2,j2)
    IOF%runoff_hflx(i,j)  = LIB%runoff_hflx(i2,j2)
    IOF%calving_hflx(i,j) = LIB%calving_hflx(i2,j2)
    ! diagnostic fluxes...
    IOF%calving_preberg(i,j) = IOF%calving(i,j)
    IOF%calving_hflx_preberg(i,j) = IOF%calving_hflx(i,j)
  enddo ; enddo

  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IOF,Ice,i_off,j_off) &
!$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    IOF%flux_u_ocn(i,j) = Ice%flux_u(i2,j2)
    IOF%flux_v_ocn(i,j) = Ice%flux_v(i2,j2)
  enddo ; enddo

end subroutine set_ice_ocean_fluxes

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ocean_top_fluxes translates ice-bottom fluxes of heat, mass, salt, and
!!  tracers from the ice model's internal state to the public ice data type
!!  for use by the ocean model.
subroutine set_ocean_top_fluxes(Ice, IST, IOF, G, IG, sCS)
  type(ice_data_type),       intent(inout) :: Ice
  type(ice_state_type),      intent(inout) :: IST
  type(ice_ocean_flux_type), intent(in)    :: IOF
  type(SIS_hor_grid_type),   intent(inout) :: G
  type(ice_grid_type),       intent(in)    :: IG
  type(SIS_slow_CS),         intent(in)    :: sCS

  real :: I_count
  integer :: i, j, k, isc, iec, jsc, jec, m, n
  integer :: i2, j2, i_off, j_off, ind, ncat, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce
  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc

  if (sCS%debug) then
    call Ice_public_type_chksum("Start set_ocean_top_fluxes", Ice)
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

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,Ice,IST,i_off,j_off) &
!$OMP                           private(i2,j2)
  do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
    Ice%part_size(i2,j2,k+1) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo

! ### This call seems out of place here.
  if (IST%id_slp>0) then
    call enable_SIS_averaging(real(time_type_to_real(IST%Time_step_slow)), IST%Time, IST%diag)
    call post_data(IST%id_slp, Ice%p_surf, IST%diag)
    call disable_SIS_averaging(IST%diag)
  endif

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,Ice,IST,IOF,i_off,j_off,G) &
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
    Ice%runoff(i2,j2)  = IOF%runoff(i,j)
    Ice%calving(i2,j2) = IOF%calving(i,j)
    Ice%runoff_hflx(i2,j2)  = IOF%runoff_hflx(i,j)
    Ice%calving_hflx(i2,j2) = IOF%calving_hflx(i,j)
    Ice%flux_salt(i2,j2) = IOF%flux_salt(i,j)

    if (IOF%slp2ocean) then
      Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) - 1e5 ! SLP - 1 std. atmosphere, in Pa.
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

  do n=1,Ice%ocean_fluxes%num_bcs ; do m=1,Ice%ocean_fluxes%bc(n)%num_fields
    ind = IOF%tr_flux_index(m,n)
    if (ind < 1) call SIS_error(FATAL, "Bad boundary flux index in set_ocean_top_fluxes.")
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off  ! Use these to correct for indexing differences.
        Ice%ocean_fluxes%bc(n)%field(m)%values(i2,j2) = IOF%tr_flux_ocn_top(i,j,ind)
    enddo ; enddo
  enddo ; enddo

  if (sCS%debug) then
    call Ice_public_type_chksum("End set_ocean_top_fluxes", Ice)
  endif

end subroutine set_ocean_top_fluxes

!
! Coupler interface to provide ocean surface data to atmosphere.
!
subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
  type(ocean_ice_boundary_type), intent(inout) :: Ocean_boundary
  type(ice_data_type),           intent(inout) :: Ice

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(iceClock1)

  call unpack_ocn_ice_bdry(Ocean_boundary, Ice%OSS, Ice%G, &
                           Ice%ocean_fields)

  !### Exchange information from the slow ice processors to the fast ice processors.

  call set_ice_surface_state(Ice, Ice%Ice_state, Ocean_boundary%t, &
                             Ice%OSS, Ice%Rad, Ice%FIA, Ice%G, Ice%IG, Ice%fCS )

  call mpp_clock_end(iceClock1) ; call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow_up

!> This subroutine converts the information in a publicly visible
!! ocean_ice_boundary_type into an internally visible ocean_sfc_state_type
!! variable.
subroutine unpack_ocn_ice_bdry(OIB, OSS, G, ocean_fields)
  type(ocean_ice_boundary_type), intent(in)    :: OIB
  type(ocean_sfc_state_type),    intent(inout) :: OSS
  type(SIS_hor_grid_type),       intent(inout) :: G
  type(coupler_3d_bc_type),      intent(inout) :: ocean_fields

  real, dimension(SZI_(G),SZJ_(G)) :: u_nonsym, v_nonsym
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  logical :: Cgrid_ocn
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = LBOUND(OIB%t,1) - G%isc ; j_off = LBOUND(OIB%t,2) - G%jsc

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,OSS,OIB,i_off,j_off) &
!$OMP                           private(i2,j2)
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off

    OSS%t_ocn(i,j) = OIB%t(i2,j2) - T_0degC
    OSS%s_surf(i,j) = OIB%s(i2,j2)
    OSS%frazil(i,j) = OIB%frazil(i2,j2)
    OSS%sea_lev(i,j) = OIB%sea_level(i2,j2)
  enddo ; enddo

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
      call pass_vector(OSS%u_ocn_C, OSS%v_ocn_C, G%Domain, stagger=CGRID_NE)
    else
      do J=jsc-1,jec ; do I=isc-1,iec
        OSS%u_ocn_B(I,J) = 0.25*((u_nonsym(i,j) + u_nonsym(i+1,j+1)) + &
                               (u_nonsym(i+1,j) + u_nonsym(i,j+1)))
        OSS%v_ocn_B(I,J) = 0.25*((v_nonsym(i,j) + v_nonsym(i+1,j+1)) + &
                               (v_nonsym(i+1,j) + v_nonsym(i,j+1)))
      enddo ; enddo
      call pass_vector(OSS%u_ocn_B, OSS%v_ocn_B, G%Domain, stagger=BGRID_NE)
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
      call pass_vector(OSS%u_ocn_C, OSS%v_ocn_C, G%Domain, stagger=CGRID_NE)
    else
      do J=jsc,jec ; do I=isc,iec ; i2 = i+i_off ; j2 = j+j_off
        OSS%u_ocn_B(I,J) = OIB%u(i2,j2)
        OSS%v_ocn_B(I,J) = OIB%v(i2,j2)
      enddo ; enddo
      if (G%symmetric) &
        call fill_symmetric_edges(OSS%u_ocn_B, OSS%v_ocn_B, G%Domain, stagger=BGRID_NE)

      call pass_vector(OSS%u_ocn_B, OSS%v_ocn_B, G%Domain, stagger=BGRID_NE)
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

      call pass_vector(OSS%u_ocn_C, OSS%v_ocn_C, G%Domain, stagger=CGRID_NE)
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
      call pass_vector(OSS%u_ocn_B, OSS%v_ocn_B, G%Domain, stagger=BGRID_NE)
    endif
  else
    call SIS_error(FATAL, "set_ice_surface_state: Unrecognized OIB%stagger.")
  endif

  call pass_var(OSS%sea_lev, G%Domain)

! Transfer the ocean state for extra tracer fluxes.
  do n=1,OIB%fields%num_bcs  ; do m=1,OIB%fields%bc(n)%num_fields
    ocean_fields%bc(n)%field(m)%values(:,:,1) = OIB%fields%bc(n)%field(m)%values
  enddo ; enddo

end subroutine unpack_ocn_ice_bdry

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ice_surface_state - prepare surface state for atmosphere fast physics    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ice_surface_state(Ice, IST, t_surf_ice_bot, OSS, Rad, FIA, G, IG, fCS)
  type(ice_data_type),        intent(inout) :: Ice
  type(ice_state_type),       intent(inout) :: IST
  type(ocean_sfc_state_type), intent(in)    :: OSS
  type(ice_rad_type),         intent(inout) :: Rad
  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(in)    :: IG
  type(SIS_fast_CS),          intent(in)    :: fCS
  real, dimension(G%isc:G%iec,G%jsc:G%jec), &
                              intent(in)    :: t_surf_ice_bot

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: m_ice_tot
  real, dimension(IG%NkIce) :: sw_abs_lay
  real :: u, v
  real :: area_pt
  real :: I_Nk
  real :: kg_H_Nk  ! The conversion factor from units of H to kg/m2 over Nk.
  real :: dt_slow  ! The thermodynamic step, in s.
  real :: Idt_slow ! The inverse of the thermodynamic step, in s-1.
  type(time_type) :: dt_r   ! A temporary radiation timestep.

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  logical :: sent
  real :: H_to_m_ice     ! The specific volumes of ice and snow times the
  real :: H_to_m_snow    ! conversion factor from thickness units, in m H-1.
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc
  I_Nk = 1.0 / IG%NkIce ; kg_H_Nk = IG%H_to_kg_m2 * I_Nk

  H_to_m_snow = IG%H_to_kg_m2 / IST%Rho_snow ; H_to_m_ice = IG%H_to_kg_m2 / IST%Rho_ice

  ! Pass the ocean state through ice on partition 0, unless using specified ice.
  if (.not. IST%specified_ice) then
    IST%t_surf(isc:iec,jsc:jec,0) = t_surf_ice_bot(isc:iec,jsc:jec)
  endif

  if (fCS%bounds_check) &
    call IST_bounds_check(IST, G, IG, "Start of set_ice_surface_state", OSS=OSS)

  if (fCS%debug) then
    call IST_chksum("Start set_ice_surface_state", IST, G, IG)
    call Ice_public_type_chksum("Start set_ice_surface_state", Ice)
  endif

  m_ice_tot(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,IST,OSS,FIA,ncat,m_ice_tot,i_off,j_off)
  do j=jsc,jec

    do k=1,ncat ; do i=isc,iec
      FIA%tmelt(i,j,k) = 0.0 ; FIA%bmelt(i,j,k) = 0.0
      m_ice_tot(i,j) = m_ice_tot(i,j) + IST%mH_ice(i,j,k) * IST%part_size(i,j,k)
    enddo ; enddo

    do i=isc,iec
      if (m_ice_tot(i,j) > 0.0) then
        FIA%bheat(i,j) = OSS%kmelt*(OSS%t_ocn(i,j) - T_Freeze(OSS%s_surf(i,j), IST%ITV))
      else
        FIA%bheat(i,j) = 0.0
      endif
      FIA%frazil_left(i,j) = OSS%frazil(i,j)
    enddo
  enddo

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
  dT_r = IST%Time_step_slow
  if (Rad%frequent_albedo_update) dT_r = IST%Time_step_fast
  call set_ocean_albedo(Ice, Rad%do_sun_angle_for_alb, G, IST%Time, &
                        IST%Time + dT_r, Rad%coszen_nextrad)

  if (IST%slab_ice) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Ice,i_off,j_off, &
!$OMP                                  H_to_m_snow,H_to_m_ice,OSS) &
!$OMP                          private(i2,j2,k2)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      call slab_ice_optics(IST%mH_snow(i,j,k)*H_to_m_snow, IST%mH_ice(i,j,k)*H_to_m_ice, &
               IST%t_surf(i,j,k)-T_0degC, T_Freeze(OSS%s_surf(i,j),IST%ITV), &
               Ice%albedo(i2,j2,k2))

      Ice%albedo_vis_dir(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_vis_dif(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_nir_dir(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_nir_dif(i2,j2,k2) = Ice%albedo(i2,j2,k2)
    endif ; enddo ; enddo ; enddo
  else
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Ice,G,IG,i_off,j_off, &
!$OMP                                  H_to_m_snow,H_to_m_ice,OSS) &
!$OMP                          private(i2,j2,k2,sw_abs_lay)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      call ice_optics_SIS2(IST%mH_pond(i,j,k), IST%mH_snow(i,j,k)*H_to_m_snow, &
               IST%mH_ice(i,j,k)*H_to_m_ice, &
               IST%t_surf(i,j,k)-T_0degC, T_Freeze(OSS%s_surf(i,j),IST%ITV), IG%NkIce, &
               Ice%albedo_vis_dir(i2,j2,k2), Ice%albedo_vis_dif(i2,j2,k2), &
               Ice%albedo_nir_dir(i2,j2,k2), Ice%albedo_nir_dif(i2,j2,k2), &
               Rad%sw_abs_sfc(i,j,k),  Rad%sw_abs_snow(i,j,k), &
               sw_abs_lay, Rad%sw_abs_ocn(i,j,k), Rad%sw_abs_int(i,j,k), &
               IST%ice_thm_CSp, IST%ITV, coszen_in=Rad%coszen_nextrad(i,j))

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
    call IST_bounds_check(IST, G, IG, "Midpoint set_ice_surface_state", OSS=OSS)

  ! Copy the surface temperatures into the externally visible data type.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,ncat,i_off,j_off,OSS) &
!$OMP                          private(i2,j2,k2)
  do j=jsc,jec
    do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      Ice%s_surf(i2,j2) = OSS%s_surf(i,j)
    enddo
    do k=0,ncat ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
      Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
    enddo ; enddo
  enddo

  ! put ocean and ice velocities into Ice%u_surf/v_surf on t-cells
  if (IST%Cgrid_dyn) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,G,i_off,j_off,OSS) &
!$OMP                          private(i2,j2)
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      if (G%mask2dT(i,j) > 0.5 ) then
        Ice%u_surf(i2,j2,1) = 0.5*(OSS%u_ocn_C(I,j) + OSS%u_ocn_C(I-1,j))
        Ice%v_surf(i2,j2,1) = 0.5*(OSS%v_ocn_C(i,J) + OSS%v_ocn_C(i,J-1))
        Ice%u_surf(i2,j2,2) = 0.5*(IST%u_ice_C(I,j) + IST%u_ice_C(I-1,j))
        Ice%v_surf(i2,j2,2) = 0.5*(IST%v_ice_C(i,J) + IST%v_ice_C(i,J-1))
      else
        Ice%u_surf(i2,j2,1) = 0.0 ; Ice%v_surf(i2,j2,1) = 0.0
        Ice%u_surf(i2,j2,2) = 0.0 ; Ice%v_surf(i2,j2,2) = 0.0
      endif
    enddo ; enddo
  else ! B-grid discretization.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,G,i_off,j_off,OSS) &
!$OMP                          private(i2,j2)
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      if (G%mask2dT(i,j) > 0.5 ) then
        Ice%u_surf(i2,j2,1) = 0.25*((OSS%u_ocn_B(I,J) + OSS%u_ocn_B(I-1,J-1)) + &
                                    (OSS%u_ocn_B(I,J-1) + OSS%u_ocn_B(I-1,J)) )
        Ice%v_surf(i2,j2,1) = 0.25*((OSS%v_ocn_B(I,J) + OSS%v_ocn_B(I-1,J-1)) + &
                                    (OSS%v_ocn_B(I,J-1) + OSS%v_ocn_B(I-1,J)) )
        Ice%u_surf(i2,j2,2) = 0.25*((IST%u_ice_B(I,J) + IST%u_ice_B(I-1,J-1)) + &
                                    (IST%u_ice_B(I,J-1) + IST%u_ice_B(I-1,J)) )
        Ice%v_surf(i2,j2,2) = 0.25*((IST%v_ice_B(I,J) + IST%v_ice_B(I-1,J-1)) + &
                                    (IST%v_ice_B(I,J-1) + IST%v_ice_B(I-1,J)) )
      else
        Ice%u_surf(i2,j2,1) = 0.0 ; Ice%v_surf(i2,j2,1) = 0.0
        Ice%u_surf(i2,j2,2) = 0.0 ; Ice%v_surf(i2,j2,2) = 0.0
      endif
    enddo ; enddo
  endif

  if (fCS%debug) then
    call chksum(Ice%u_surf(:,:,1), "Intermed Ice%u_surf(1)")
    call chksum(Ice%v_surf(:,:,1), "Intermed Ice%v_surf(1)")
    call chksum(Ice%u_surf(:,:,2), "Intermed Ice%u_surf(2)")
    call chksum(Ice%v_surf(:,:,2), "Intermed Ice%v_surf(2)")
    call chksum(G%mask2dT(isc:iec,jsc:jec), "Intermed G%mask2dT")
    if (allocated(OSS%u_ocn_C)) &
      call uchksum(OSS%u_ocn_C, "OSS%u_ocn_C", G%HI, haloshift=1)
    if (allocated(OSS%v_ocn_C)) &
      call vchksum(OSS%v_ocn_C, "OSS%v_ocn_C", G%HI, haloshift=1)
    if (allocated(OSS%u_ocn_B)) &
      call Bchksum(OSS%u_ocn_B, "OSS%u_ocn_B", G%HI, haloshift=1)
    if (allocated(OSS%v_ocn_B)) &
      call Bchksum(OSS%v_ocn_B, "OSS%v_ocn_B", G%HI, haloshift=1)
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

  !
  ! Pre-timestep diagnostics
  !
  dt_slow = time_type_to_real(IST%Time_step_slow)
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)
  if (Rad%id_alb>0) call post_avg(Rad%id_alb, Ice%albedo, &
                     IST%part_size(isc:iec,jsc:jec,:), IST%diag)
  if (OSS%id_sst>0) call post_data(OSS%id_sst, OSS%t_ocn, IST%diag)
  if (OSS%id_sss>0) call post_data(OSS%id_sss, OSS%s_surf, IST%diag)
  if (OSS%id_ssh>0) call post_data(OSS%id_ssh, OSS%sea_lev, IST%diag)
  if (IST%Cgrid_dyn) then
    if (OSS%id_uo>0) call post_data(OSS%id_uo, OSS%u_ocn_C, IST%diag)
    if (OSS%id_vo>0) call post_data(OSS%id_vo, OSS%v_ocn_C, IST%diag)
  else
    if (OSS%id_uo>0) call post_data(OSS%id_uo, OSS%u_ocn_B, IST%diag)
    if (OSS%id_vo>0) call post_data(OSS%id_vo, OSS%v_ocn_B, IST%diag)
  endif
  if (OSS%id_frazil>0) &
    call post_data(OSS%id_frazil, OSS%frazil*Idt_slow, IST%diag)

  if (FIA%id_bheat>0) call post_data(FIA%id_bheat, FIA%bheat, IST%diag)
  call disable_SIS_averaging(IST%diag)

  if (fCS%debug) then
    call IST_chksum("End set_ice_surface_state", IST, G, IG)
    call Ice_public_type_chksum("End set_ice_surface_state", Ice)
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

  call mpp_clock_begin(iceClock) ; call mpp_clock_begin(iceClock3)

  if (Ice%fCS%debug) &
    call Ice_public_type_chksum("Pre do_update_ice_model_fast", Ice)

  dT_fast = Ice%Ice_state%Time_step_fast
  Time_start = Ice%Ice_state%Time
  Time_end = Time_start + dT_fast

  if (Ice%Rad%add_diurnal_sw) &
    call add_diurnal_sw(Atmos_boundary, Ice%G, Time_start, Time_end)

  call do_update_ice_model_fast(Atmos_boundary, Ice%Ice_state, Ice%OSS, Ice%Rad, &
                                Ice%FIA, Ice%fast_thermo_CSp, &
                                Ice%G, Ice%IG )

  Ice%Time = Ice%Ice_state%Time
  Time_end = Ice%Ice_state%Time ! Probably there is no change to Time_end.

  call fast_radiation_diagnostics(Atmos_boundary, Ice, Ice%Ice_state, Ice%Rad, &
                                  Ice%G, Ice%IG, Time_start, Time_end)

  ! Set some of the evolving ocean properties that will be seen by the
  ! atmosphere in the next time-step.
  call set_fast_ocean_sfc_properties(Atmos_boundary, Ice, Ice%Ice_state, Ice%Rad, &
                                     Ice%G, Ice%IG, Time_end, Time_end + dT_fast)

  if (Ice%fCS%debug) &
    call Ice_public_type_chksum("End do_update_ice_model_fast", Ice)
  if (Ice%fCS%bounds_check) &
    call Ice_public_type_bounds_check(Ice, Ice%G, "End update_ice_fast")

  call mpp_clock_end(iceClock3) ; call mpp_clock_end(iceClock)

end subroutine update_ice_model_fast

subroutine set_fast_ocean_sfc_properties( Atmos_boundary, Ice, IST, Rad, G, IG, Time_start, Time_end)
  type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary
  type(ice_data_type),           intent(inout) :: Ice
  type(ice_state_type),          intent(inout) :: IST
  type(ice_rad_type),            intent(inout) :: Rad
  type(SIS_hor_grid_type),       intent(inout) :: G
  type(ice_grid_type),           intent(inout) :: IG
  type(time_type),               intent(in)    :: Time_start, Time_end
   
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
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,Ice,IST,Atmos_boundary,io_A,jo_A,io_I,jo_I ) &
!$OMP                           private(i2,j2,i3,j3)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+io_I ; j2 = j+jo_I ; i3 = i+io_A ; j3 = j+jo_A
    Rad%coszen_nextrad(i,j) = Atmos_boundary%coszen(i3,j3,1)
    Ice%p_surf(i2,j2) = Atmos_boundary%p(i3,j3,1)
  enddo ; enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Ice,io_I,jo_I ) &
!$OMP                           private(i2,j2,k2)
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+io_I ; j2 = j+jo_I ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
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
  real, dimension(SZI_(G), SZJ_(G)), &
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


subroutine fast_radiation_diagnostics(ABT, Ice, IST, Rad, G, IG, Time_start, Time_end)
  type(atmos_ice_boundary_type), intent(in)    :: ABT
  type(ice_data_type),           intent(in)    :: Ice
  type(ice_state_type),          intent(inout) :: IST
  type(ice_rad_type),            intent(inout) :: Rad
  type(SIS_hor_grid_type),       intent(inout) :: G
  type(ice_grid_type),           intent(in)    :: IG
  type(time_type),               intent(in)    :: Time_start, Time_end

  real, dimension(SZI_(G), SZJ_(G)) :: tmp_diag
  real :: dt_diag
  real    :: Stefan ! The Stefan-Boltzmann constant in W m-2 K-4 as used for
                    ! strictly diagnostic purposes.
  integer :: i, j, k, m, i2, j2, k2, i3, j3, isc, iec, jsc, jec, ncat, NkIce
  integer :: io_A, jo_A, io_I, jo_I  ! Offsets for indexing conventions.

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  NkIce = IG%NkIce
  io_A = LBOUND(ABT%t_flux,1) - G%isc ; jo_A = LBOUND(ABT%t_flux,2) - G%jsc
  io_I = LBOUND(Ice%t_surf,1) - G%isc ; jo_I = LBOUND(Ice%t_surf,2) - G%jsc

  dt_diag = time_type_to_real(Time_end - Time_start)

  call enable_SIS_averaging(dt_diag, Time_end, IST%diag)

  if (Rad%id_alb_vis_dir>0) call post_avg(Rad%id_alb_vis_dir, Ice%albedo_vis_dir, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag)
  if (Rad%id_alb_vis_dif>0) call post_avg(Rad%id_alb_vis_dif, Ice%albedo_vis_dif, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag)
  if (Rad%id_alb_nir_dir>0) call post_avg(Rad%id_alb_nir_dir, Ice%albedo_nir_dir, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag)
  if (Rad%id_alb_nir_dif>0) call post_avg(Rad%id_alb_nir_dif, Ice%albedo_nir_dif, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag)

  if (Rad%id_sw_abs_sfc>0) call post_avg(Rad%id_sw_abs_sfc, Rad%sw_abs_sfc, &
                                   IST%part_size(:,:,1:), IST%diag, G=G)
  if (Rad%id_sw_abs_snow>0) call post_avg(Rad%id_sw_abs_snow, Rad%sw_abs_snow, &
                                   IST%part_size(:,:,1:), IST%diag, G=G)
  do m=1,NkIce
    if (Rad%id_sw_abs_ice(m)>0) call post_avg(Rad%id_sw_abs_ice(m), Rad%sw_abs_ice(:,:,:,m), &
                                     IST%part_size(:,:,1:), IST%diag, G=G)
  enddo
  if (Rad%id_sw_abs_ocn>0) call post_avg(Rad%id_sw_abs_ocn, Rad%sw_abs_ocn, &
                                   IST%part_size(:,:,1:), IST%diag, G=G)

  if (Rad%id_sw_pen>0) then
    tmp_diag(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,tmp_diag)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * &
                     (Rad%sw_abs_ocn(i,j,k) + Rad%sw_abs_int(i,j,k))
    enddo ; enddo ; enddo
    call post_data(Rad%id_sw_pen, tmp_diag, IST%diag)
  endif

  if (Rad%id_lwdn > 0) then
    tmp_diag(:,:) = 0.0
    Stefan = 5.6734e-8  ! Set the Stefan-Bolzmann constant, in W m-2 K-4.
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
      i3 = i+io_A ; j3 = j+jo_A ; k2 = k+1
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * &
                           (ABT%lw_flux(i3,j3,k2) + Stefan*IST%t_surf(i,j,k)**4)
    endif ; enddo ; enddo ; enddo
    call post_data(Rad%id_lwdn, tmp_diag, IST%diag)
  endif

  if (Rad%id_swdn > 0) then
    tmp_diag(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,IST,Ice,ABT, &
!$OMP                                  io_I,jo_I,io_A,jo_A,tmp_diag) &
!$OMP                          private(i2,j2,k2,i3,j3)
    do j=jsc,jec ; do k=0,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
      i2 = i+io_I ; j2 = j+jo_I ; i3 = i+io_A ; j3 = j+jo_A ; k2 = k+1
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * ( &
            (ABT%sw_flux_vis_dir(i3,j3,k2)/(1-Ice%albedo_vis_dir(i2,j2,k2)) + &
             ABT%sw_flux_vis_dif(i3,j3,k2)/(1-Ice%albedo_vis_dif(i2,j2,k2))) + &
            (ABT%sw_flux_nir_dir(i3,j3,k2)/(1-Ice%albedo_nir_dir(i2,j2,k2)) + &
             ABT%sw_flux_nir_dif(i3,j3,k2)/(1-Ice%albedo_nir_dif(i2,j2,k2))) )
    endif ; enddo ; enddo ; enddo
    call post_data(Rad%id_swdn, tmp_diag, IST%diag)
  endif

  if (Rad%id_coszen>0) call post_data(Rad%id_coszen, Rad%coszen_nextrad, IST%diag)

  call disable_SIS_averaging(IST%diag)

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
  enddo ; enddo

end subroutine add_diurnal_sw

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_init - initializes ice model data, parameters and diagnostics      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_init(Ice, Time_Init, Time, Time_step_fast, Time_step_slow )

  type(ice_data_type), intent(inout) :: Ice
  type(time_type)    , intent(in)    :: Time_Init      ! starting time of model integration
  type(time_type)    , intent(in)    :: Time           ! current time
  type(time_type)    , intent(in)    :: Time_step_fast ! time step for the ice_model_fast
  type(time_type)    , intent(in)    :: Time_step_slow ! time step for the ice_model_slow

! This include declares and sets the variable "version".
#include "version_variable.h"
  real :: hlim_dflt(8) = (/ 1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! lower thickness limits 1...CatIce
  real :: enth_spec_snow, enth_spec_ice
  real, allocatable :: S_col(:)
  real :: pi ! pi = 3.1415926... calculated as 4*atan(1)
  integer :: i, j, k, l, i2, j2, k2, i_off, j_off, n
  integer :: isc, iec, jsc, jec, nCat_dflt
  logical :: spec_thermo_sal
  character(len=120) :: restart_file
  character(len=240) :: restart_path
  character(len=40)  :: mod = "ice_model" ! This module's name.
  character(len=8)   :: nstr
  type(directories)  :: dirs   ! A structure containing several relevant directory paths.

  type(param_file_type) :: param_file
  type(hor_index_type)  :: HI  !  A hor_index_type for array extents
  type(ice_state_type),    pointer :: IST => NULL()
  type(SIS_hor_grid_type), pointer :: G => NULL()
  type(ice_grid_type),     pointer :: IG => NULL()
  type(dyn_horgrid_type),  pointer :: dG => NULL()

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
  real :: rho_Ocean      ! The nominal density of seawater, in kg m-3.
  real :: kmelt          ! A constant that is used in the calculation of the
                         ! ocean/ice basal heat flux, in W m-2 K-1.  This could
                         ! be changed to reflect the turbulence in the under-ice
                         ! ocean boundary layer and the effective depth of the
                         ! reported value of t_ocn.

  integer :: idr, id_sal
  integer :: write_geom
  logical :: test_grid_copy = .false.
  logical :: nudge_sea_ice
  logical :: atmos_winds, slp2ocean
  logical :: do_icebergs, pass_iceberg_area_to_ocean
  logical :: do_ridging
  logical :: debug, bounds_check
  logical :: do_sun_angle_for_alb, add_diurnal_sw
  logical :: write_geom_files  ! If true, write out the grid geometry files.
  logical :: symmetric         ! If true, use symmetric memory allocation.
  logical :: global_indexing   ! If true use global horizontal index values instead
                               ! of having the data domain on each processor start at 1.
  integer :: first_direction   ! An integer that indicates which direction is to be
                               ! updated first in directionally split parts of the
                               ! calculation.  This can be altered during the course
                               ! of the run via calls to set_first_direction.
  logical :: fast_ice_PE       ! If true, fast ice processes are handled on this PE.
  logical :: slow_ice_PE       ! If true, slow ice processes are handled on this PE.
  logical :: read_aux_restart
  logical :: is_restart = .false.
  character(len=16)  :: stagger, dflt_stagger

  ! ### These are just here to keep the order of SIS_parameter_doc.
  logical :: column_check
  real :: imb_tol

  if (associated(Ice%Ice_state)) then
    call SIS_error(WARNING, "ice_model_init called with an associated "// &
                    "Ice%Ice_state structure. Model is already initialized.")
    return
  endif
  ! For now, both fast and slow processes occur on all sea-ice PEs.
  fast_ice_PE = .true. ; slow_ice_PE = .true.

  if (fast_ice_PE) then
    if (.not.associated(Ice%fCS)) allocate(Ice%fCS)
  endif
  if (slow_ice_PE) then
    if (.not.associated(Ice%sCS)) allocate(Ice%sCS)
  endif

  if (.not.associated(Ice%Ice_state)) allocate(Ice%Ice_state) ; IST => Ice%Ice_state
  if (.not.associated(Ice%G)) allocate(Ice%G)
  if (test_grid_copy) then ; allocate(G)
  else ; G => Ice%G ; endif
  if (.not.associated(Ice%IG)) allocate(Ice%IG) ; IG => Ice%IG
  if (.not.associated(Ice%Ice_restart)) allocate(Ice%Ice_restart)

  ! Open the parameter file.
  call Get_SIS_Input(param_file, dirs)

  call callTree_enter("ice_model_init(), ice_model.F90")

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "SPECIFIED_ICE", IST%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  call get_param(param_file, mod, "CGRID_ICE_DYNAMICS", IST%Cgrid_dyn, &
                 "If true, use a C-grid discretization of the sea-ice \n"//&
                 "dynamics; if false use a B-grid discretization.", &
                 default=.false.)
  call get_param(param_file, mod, "USE_SLAB_ICE", IST%slab_ice, &
                 "If true, use the very old slab-style ice.", default=.false.)

  call obsolete_logical(param_file, "SIS1_5L_THERMODYNAMICS", warning_val=.false.)
  call obsolete_logical(param_file, "INTERSPERSED_ICE_THERMO", warning_val=.false.)
  call obsolete_logical(param_file, "AREA_WEIGHTED_STRESSES", warning_val=.true.)

  dflt_stagger = "B" ; if (IST%Cgrid_dyn) dflt_stagger = "C"
  call get_param(param_file, mod, "ICE_OCEAN_STRESS_STAGGER", stagger, &
                 "A case-insensitive character string to indicate the \n"//&
                 "staggering of the stress field on the ocean that is \n"//&
                 "returned to the coupler.  Valid values include \n"//&
                 "'A', 'B', or 'C', with a default that follows the \n"//&
                 "value of CGRID_ICE_DYNAMICS.", default=dflt_stagger)
  if (uppercase(stagger(1:1)) == 'A') then ; Ice%flux_uv_stagger = AGRID
  elseif (uppercase(stagger(1:1)) == 'B') then ; Ice%flux_uv_stagger = BGRID_NE
  elseif (uppercase(stagger(1:1)) == 'C') then ; Ice%flux_uv_stagger = CGRID_NE
  else ; call SIS_error(FATAL,"ice_model_init: ICE_OCEAN_STRESS_STAGGER = "//&
                        trim(stagger)//" is invalid.") ; endif

  ! Rho_ocean is not actually used here, but it used from later get_param
  ! calls in other modules.  This call is here to avoid changing the order of
  ! the entries in the SIS_parameter_doc files.
  call get_param(param_file, mod, "RHO_OCEAN", Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0)
  call get_param(param_file, mod, "RHO_ICE", IST%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0)
  call get_param(param_file, mod, "RHO_SNOW", IST%Rho_snow, &
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

!  call get_param(param_file, mod, "ICE_SEES_ATMOS_WINDS", Ice%FIA%atmos_winds, &
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
  call get_param(param_file, mod, "DO_ICEBERGS", do_icebergs, &
                 "If true, call the iceberg module.", default=.false.)
  if (do_icebergs) then
    call get_param(param_file, mod, "PASS_ICEBERG_AREA_TO_OCEAN", pass_iceberg_area_to_ocean, &
                 "If true, iceberg area is passed through coupler", default=.false.)
  else ; pass_iceberg_area_to_ocean = .false. ; endif
  if (slow_ice_PE) then
    Ice%sCS%do_icebergs = do_icebergs
    Ice%sCS%pass_iceberg_area_to_ocean = pass_iceberg_area_to_ocean
  endif
  
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
  if (write_geom<0 .or. write_geom>2) call SIS_error(FATAL,"SIS2: "//&
         "WRITE_GEOM must be equal to 0, 1 or 2.")
  write_geom_files = ((write_geom==2) .or. ((write_geom==1) .and. &
     ((dirs%input_filename(1:1)=='n') .and. (LEN_TRIM(dirs%input_filename)==1))))

  if (IST%specified_ice) IST%slab_ice = .true.

  ! Set up the ice-specific grid describing categories and ice layers.
  nCat_dflt = 5 ; if (IST%slab_ice)  nCat_dflt = 1 ! open water and ice ... but never in same place
  call set_ice_grid(Ice%IG, param_file, nCat_dflt)
  if (IST%slab_ice) IG%CatIce = 1 ! open water and ice ... but never in same place
  ! Initialize IG%cat_thick_lim here.  ###This needs to be extended to add more options.
  do k=1,min(IG%CatIce+1,size(hlim_dflt(:)))
    IG%cat_thick_lim(k) = hlim_dflt(k)
  enddo
  if ((IG%CatIce+1 > size(hlim_dflt(:))) .and. (size(hlim_dflt(:)) > 1)) then
    do k=min(IG%CatIce+1,size(hlim_dflt(:))) + 1, IG%CatIce+1
      IG%cat_thick_lim(k) =  2.0*IG%cat_thick_lim(k-1) - IG%cat_thick_lim(k-2)
    enddo
  endif
  do k=1,IG%CatIce+1
    IG%mH_cat_bound(k) = IG%cat_thick_lim(k) * (IST%Rho_ice*IG%kg_m2_to_H)
  enddo

  ! Set up the domains and lateral grids.

  ! Set up the MOM_domain_type structures.
#ifdef SYMMETRIC_MEMORY_
  symmetric = .true.
#else
  symmetric = .false.
#endif
#ifdef STATIC_MEMORY_
  call MOM_domains_init(G%domain, param_file, symmetric=symmetric, &
            static_memory=.true., NIHALO=NIHALO_, NJHALO=NJHALO_, &
            NIGLOBAL=NIGLOBAL_, NJGLOBAL=NJGLOBAL_, NIPROC=NIPROC_, &
            NJPROC=NJPROC_, domain_name="ice model", include_name="SIS2_memory.h")
#else
  call MOM_domains_init(G%domain, param_file, symmetric=symmetric, &
           domain_name="ice model", include_name="SIS2_memory.h")
#endif

  call callTree_waypoint("domains initialized (ice_model_init)")
  call hor_index_init(G%Domain, HI, param_file, &
                      local_indexing=.not.global_indexing)

  call create_dyn_horgrid(dG, HI) !, bathymetry_at_vel=bathy_at_vel)
  call clone_MOM_domain(G%Domain, dG%Domain)

  ! Set the bathymetry, Coriolis parameter, open channel widths and masks.
  call SIS_initialize_fixed(dG, param_file, write_geom_files, dirs%output_directory)

  call set_hor_grid(G, param_file, global_indexing=global_indexing)
  call copy_dyngrid_to_SIS_horgrid(dG, G)
  call destroy_dyn_horgrid(dG)

  call set_domain(G%Domain%mpp_domain)
  ! Allocate and register fields for restarts.
  call ice_data_type_register_restarts(G%domain%mpp_domain, IG%CatIce, &
                         param_file, Ice, Ice%Ice_restart, restart_file)

  call ice_state_register_restarts(G%domain%mpp_domain, HI, IG, param_file, &
                                   IST, Ice%Ice_restart, restart_file)

  call alloc_ocean_sfc_state(Ice%OSS, HI, IST%Cgrid_dyn)
  Ice%OSS%kmelt = kmelt

  if (slow_ice_PE) then
    call alloc_ice_ocean_flux(Ice%IOF, HI, do_iceberg_fields=Ice%sCS%do_icebergs)
    Ice%IOF%slp2ocean = slp2ocean
  endif

  call alloc_fast_ice_avg(Ice%FIA, HI, IG)
  Ice%FIA%atmos_winds = atmos_winds

  call ice_rad_register_restarts(G%domain%mpp_domain, HI, IG, param_file, &
                                 Ice%Rad, Ice%Ice_restart, restart_file)
  Ice%Rad%do_sun_angle_for_alb = do_sun_angle_for_alb
  Ice%Rad%add_diurnal_sw = add_diurnal_sw
  Ice%Rad%frequent_albedo_update = .true.
  !### Instead perhaps this could be
  !###   Ice%Rad%frequent_albedo_update = Ice%Rad%do_sun_angle_for_alb .or. (Time_step_slow > dT_Rad)
  !### However this changes answers in coupled models.  I don't understand why. -RWH

  ! Ice%IOF has now been set and can be used.
  Ice%IOF%flux_uv_stagger = Ice%flux_uv_stagger

  call SIS_dyn_trans_register_restarts(G%domain%mpp_domain, HI, IG, param_file,&
                              Ice%dyn_trans_CSp, Ice%Ice_restart, restart_file)

  call SIS_diag_mediator_init(G, IG, param_file, IST%diag, component="SIS", &
                              doc_file_dir = dirs%output_directory)
  call set_SIS_axes_info(G, IG, param_file, IST%diag)

  nudge_sea_ice = .false. ; call read_param(param_file, "NUDGE_SEA_ICE", nudge_sea_ice)
  call SIS2_ice_thm_init(param_file, IST%ice_thm_CSp, IST%ITV, &
                         init_EOS=nudge_sea_ice)

  call get_SIS2_thermo_coefs(IST%ITV, enthalpy_units=enth_unit)

  ! Register tracers that will be advected around.
  call register_SIS_tracer_pair(IST%enth_ice, IG%NkIce, "enth_ice", &
                                IST%enth_snow, 1, "enth_snow", &
                                G, IG, param_file, IST%TrReg, &
                                massless_iceval=massless_ice_enth*enth_unit, &
                                massless_snowval=massless_snow_enth*enth_unit)

  if (ice_rel_salin > 0.0) then
    call register_SIS_tracer(IST%sal_ice, G, IG, IG%NkIce, "salin_ice", param_file, &
                             IST%TrReg, snow_tracer=.false., &
                             massless_val=massless_ice_salin)
  endif

  !   Register any tracers that will be handled via tracer flow control for 
  ! restarts and advection.
  call SIS_call_tracer_register(G, IG, param_file, Ice%SIS_tracer_flow_CSp, &
                                IST%diag, IST%TrReg, Ice%Ice_restart, restart_file)

  ! Redefine the computational domain sizes to use the ice model's indexing convention.
  isc = HI%isc ; iec = HI%iec ; jsc = HI%jsc ; jec = HI%jec
  i_off = LBOUND(Ice%t_surf,1) - HI%isc ; j_off = LBOUND(Ice%t_surf,2) - HI%jsc

  if (test_grid_copy) then
    !  Copy the data from the temporary grid to the dyn_hor_grid to CS%G.
    call create_dyn_horgrid(dG, G%HI)
    call clone_MOM_domain(G%Domain, dG%Domain)

    call clone_MOM_domain(G%Domain, Ice%G%Domain)
    call set_hor_grid(Ice%G, param_file)

    call copy_SIS_horgrid_to_dyngrid(G, dG)
    call copy_dyngrid_to_SIS_horgrid(dG, Ice%G)

    call destroy_dyn_horgrid(dG)
    call SIS_hor_grid_end(G) ; deallocate(G)

    G => Ice%G
  endif

  ! Set a few final things to complete the setup of the grid. 
  G%g_Earth = g_Earth
  call set_first_direction(G, first_direction)
  call clone_MOM_domain(G%domain, G%domain_aux, symmetric=.false., &
                        domain_name="ice model aux")

  ! Copy the ice model's domain into one with no halos that can be shared
  ! publicly for use by the exchange grid.
  call clone_MOM_domain(G%domain, Ice%domain, halo_size=0, symmetric=.false., &
                        domain_name="ice_nohalo")

  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%ocean_pt(i2,j2) = ( G%mask2dT(i,j) > 0.5 )
    Ice%area(i2,j2) = G%areaT(i,j) * G%mask2dT(i,j)
  enddo ; enddo

  Ice%Time           = Time
  IST%Time           = Time
  IST%Time_Init      = Time_Init
  IST%Time_step_fast = Time_step_fast
  IST%Time_step_slow = Time_step_slow

!  if (Ice%Rad%add_diurnal_sw .or. Ice%Rad%do_sun_angle_for_alb) then
!    call set_domain(G%Domain%mpp_domain)
    call astronomy_init
!    call nullify_domain()
!  endif

  !
  ! Read the restart file, if it exists.
  !
  allocate(S_col(IG%NkIce)) ; S_col(:) = 0.0
  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_unit, &
                             specified_thermo_salinity=spec_thermo_sal)

!  restart_path = 'INPUT/'//trim(restart_file)
  restart_path = trim(dirs%restart_input_dir)//trim(restart_file)
  if (file_exist(restart_path)) then
    ! Set values of IG%H_to_kg_m2 that will permit its absence from the restart
    ! file to be detected, and its difference from the value in this run to
    ! be corrected for.
    H_to_kg_m2_tmp = IG%H_to_kg_m2
    IG%H_to_kg_m2 = -1.0
    is_restart = .true.

    call restore_state(Ice%Ice_restart, directory=dirs%restart_input_dir)

    ! Approximately initialize state fields that are not present
    ! in SIS1 restart files.  This is obsolete and can probably be eliminated.

    ! Initialize the ice salinity.
    if (.not.query_initialized(Ice%Ice_restart, 'sal_ice')) then
      allocate(sal_ice_tmp(SZI_(G), SZJ_(G), IG%CatIce, IG%NkIce)) ; sal_ice_tmp(:,:,:,:) = 0.0
      do n=1,IG%NkIce
        write(nstr, '(I4)') n ; nstr = adjustl(nstr)
        id_sal = register_restart_field(Ice%Ice_restart, restart_file, 'sal_ice'//trim(nstr), &
                                     sal_ice_tmp(:,:,:,n), domain=G%domain%mpp_domain, &
                                     mandatory=.false., read_only=.true.)
        call restore_state(Ice%Ice_restart, id_sal, directory=dirs%restart_input_dir)
      enddo

      if (query_initialized(Ice%Ice_restart, 'sal_ice1')) then
        do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%sal_ice(i,j,k,1) = sal_ice_tmp(i,j,k,1)
        enddo ; enddo ; enddo
      else
        IST%sal_ice(:,:,:,1) = ice_bulk_salin
      endif
      do n=2,IG%NkIce
        write(nstr, '(I4)') n ; nstr = adjustl(nstr)
        if (query_initialized(Ice%Ice_restart, 'sal_ice'//trim(nstr))) then
          do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
            IST%sal_ice(i,j,k,n) = sal_ice_tmp(i,j,k,n)
          enddo ; enddo ; enddo
        else
          IST%sal_ice(:,:,:,n) = IST%sal_ice(:,:,:,n-1)
        endif
      enddo

      deallocate(sal_ice_tmp)
    endif

    read_aux_restart = (.not.query_initialized(Ice%Ice_restart, 'enth_ice')) .or. &
                       (.not.query_initialized(Ice%Ice_restart, 'enth_snow'))
    if (read_aux_restart) then
      allocate(t_snow_tmp(SZI_(G), SZJ_(G), IG%CatIce)) ; t_snow_tmp(:,:,:) = 0.0
      allocate(t_ice_tmp(SZI_(G), SZJ_(G), IG%CatIce, IG%NkIce)) ; t_ice_tmp(:,:,:,:) = 0.0

      idr = register_restart_field(Ice%Ice_restart, restart_file, 't_snow', t_snow_tmp, &
                                   domain=G%domain%mpp_domain, mandatory=.false., read_only=.true.)
      call restore_state(Ice%Ice_restart, idr, directory=dirs%restart_input_dir)
      do n=1,IG%NkIce
        write(nstr, '(I4)') n ; nstr = adjustl(nstr)
        idr = register_restart_field(Ice%Ice_restart, restart_file, 't_ice'//trim(nstr), &
                                     t_ice_tmp(:,:,:,n), domain=G%domain%mpp_domain, &
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
        do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,1) = Enth_from_TS(t_ice_tmp(i,j,k,1), S_col(1), IST%ITV)
        enddo ; enddo ; enddo
      else
        do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,1) = Enth_from_TS(t_ice_tmp(i,j,k,1), IST%sal_ice(i,j,k,1), IST%ITV)
        enddo ; enddo ; enddo
      endif

      do n=2,IG%NkIce
        write(nstr, '(I4)') n ; nstr = adjustl(nstr)
        if (.not.query_initialized(Ice%Ice_restart, 't_ice'//trim(nstr))) &
          t_ice_tmp(:,:,:,n) = t_ice_tmp(:,:,:,n-1)

        if (spec_thermo_sal) then
          do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
            IST%enth_ice(i,j,k,n) = Enth_from_TS(t_ice_tmp(i,j,k,n), S_col(n), IST%ITV)
          enddo ; enddo ; enddo
        else
          do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
            IST%enth_ice(i,j,k,n) = Enth_from_TS(t_ice_tmp(i,j,k,n), IST%sal_ice(i,j,k,n), IST%ITV)
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
          do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
            t_snow_tmp(i,j,k) = Temp_from_En_S(IST%enth_ice(i,j,k,1), IST%sal_ice(i,j,k,1), IST%ITV)
          enddo ; enddo ; enddo
        endif
      endif
      do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
        IST%enth_snow(i,j,k,1) = Enth_from_TS(t_snow_tmp(i,j,k), 0.0, IST%ITV)
      enddo ; enddo ; enddo
    endif

    if (read_aux_restart) deallocate(t_snow_tmp, t_ice_tmp)

    if (.not.query_initialized(Ice%Ice_restart, 'coszen')) then
      if (coszen_IC >= 0.0) then
        Ice%Rad%coszen_nextrad(:,:) = coszen_IC
      else
        rad = acos(-1.)/180.
        allocate(dummy(G%isd:G%ied,G%jsd:G%jed))
        call diurnal_solar(G%geoLatT(:,:)*rad, G%geoLonT(:,:)*rad, &
                   IST%Time, cosz=Ice%Rad%coszen_nextrad, fracday=dummy, &
                   rrsun=rrsun, dt_time=dT_rad)
        deallocate(dummy)
      endif
    endif

    H_rescale_ice = 1.0 ; H_rescale_snow = 1.0
    if (IG%H_to_kg_m2 == -1.0) then
      ! This is an older restart file, and the snow and ice thicknesses are in
      ! m, and not a mass coordinate.
      H_rescale_ice = IST%Rho_ice / H_to_kg_m2_tmp
      H_rescale_snow = IST%Rho_snow / H_to_kg_m2_tmp
    elseif (IG%H_to_kg_m2 /= H_to_kg_m2_tmp) then
      H_rescale_ice = IG%H_to_kg_m2 / H_to_kg_m2_tmp
      H_rescale_snow = H_rescale_ice
    endif
    IG%H_to_kg_m2 = H_to_kg_m2_tmp

    ! Deal with any ice masses or thicknesses over land, and rescale to
    ! account for differences between the current thickness units and whatever
    ! thickness units were in the input restart file.
    do k=1,IG%CatIce
      IST%mH_snow(:,:,k) = IST%mH_snow(:,:,k) * H_rescale_snow * G%mask2dT(:,:)
      IST%mH_ice(:,:,k) = IST%mH_ice(:,:,k) * H_rescale_ice * G%mask2dT(:,:)
    enddo

    !--- update the halo values.
    call pass_var(IST%part_size, G%Domain)
    call pass_var(IST%mH_ice, G%Domain, complete=.false.)
    call pass_var(IST%mH_snow, G%Domain, complete=.false.)
    do l=1,IG%NkIce
      call pass_var(IST%enth_ice(:,:,:,l), G%Domain, complete=.false.)
    enddo
    call pass_var(IST%enth_snow(:,:,:,1), G%Domain, complete=.true.)

    if (IST%Cgrid_dyn) then
      call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
    else
      call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
    endif

  else ! no restart implies initialization with no ice
    IST%part_size(:,:,:) = 0.0
    IST%part_size(:,:,0) = 1.0

    Ice%rough_mom(:,:,:)   = mom_rough_ice
    Ice%rough_heat(:,:,:)  = heat_rough_ice
    Ice%rough_moist(:,:,:) = heat_rough_ice
    IST%t_surf(:,:,:) = T_0degC
    IST%sal_ice(:,:,:,:) = ice_bulk_salin

    enth_spec_snow = Enth_from_TS(0.0, 0.0, IST%ITV)
    IST%enth_snow(:,:,:,1) = enth_spec_snow
    do n=1,IG%NkIce
      enth_spec_ice = Enth_from_TS(0.0, S_col(n), IST%ITV)
      IST%enth_ice(:,:,:,n) = enth_spec_ice
    enddo

    allocate(h_ice_input(G%isc:G%iec,G%jsc:G%jec))
    call get_sea_surface(IST%Time, IST%t_surf(isc:iec,jsc:jec,0), IST%part_size(isc:iec,jsc:jec,0:1), &
                         h_ice_input, ice_domain=Ice%domain )
    do j=jsc,jec ; do i=isc,iec
      IST%mH_ice(i,j,1) = h_ice_input(i,j)*(IST%Rho_ice*IG%kg_m2_to_H)
    enddo ; enddo

    !   Transfer ice to the correct thickness category.  If do_ridging=.false.,
    ! the first call to ice_redistribute has the same result.  At present, all
    ! tracers are initialized to their default values, and snow is set to 0,
    ! and so do not need to be updated here.
    if (do_ridging) then
      do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,1) > IG%mH_cat_bound(1)) then
        do k=IG%CatIce,2,-1 ; if (IST%mH_ice(i,j,1) > IG%mH_cat_bound(k-1)) then
          IST%part_size(i,j,k) = IST%part_size(i,j,1)
          IST%part_size(i,j,1) = 0.0
          IST%mH_ice(i,j,k) = IST%mH_ice(i,j,1) ; IST%mH_ice(i,j,1) = 0.0
          !  IST%mH_snow(i,j,k) = IST%mH_snow(i,j,1) ; IST%mH_snow(i,j,1) = 0.0
          exit ! from k-loop
        endif ; enddo
      endif ; enddo ; enddo
    endif

    deallocate(h_ice_input)

    call pass_var(IST%part_size, G%Domain, complete=.true. )
    call pass_var(IST%mH_ice, G%Domain, complete=.true. )

    if (coszen_IC >= 0.0) then
      Ice%Rad%coszen_nextrad(:,:) = coszen_IC
    else
      rad = acos(-1.)/180.
      allocate(dummy(G%isd:G%ied,G%jsd:G%jed))
      call diurnal_solar(G%geoLatT(:,:)*rad, G%geoLonT(:,:)*rad, &
                         IST%Time, cosz=Ice%Rad%coszen_nextrad, fracday=dummy, &
                         rrsun=rrsun, dt_time=dT_rad)
      deallocate(dummy)
    endif

  endif ! file_exist(restart_path)
  deallocate(S_col)

  do k=0,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
    Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo


  call ice_diagnostics_init(IST, Ice%IOF, Ice%OSS, Ice%FIA, Ice%Rad, G, IG, IST%diag, IST%Time)
  Ice%axes(1:2) = IST%diag%axesTc%handles(1:2)

  if (fast_ice_PE) then
    call SIS_fast_thermo_init(Ice%Time, G, IG, param_file, IST%diag, Ice%fast_thermo_CSp)
  endif

  if (slow_ice_PE) then
    call SIS_slow_thermo_init(Ice%Time, G, IG, param_file, IST%diag, Ice%slow_thermo_CSp, &
                              Ice%SIS_tracer_flow_CSp)

    call SIS_dyn_trans_init(Ice%Time, G, IG, param_file, IST%diag, &
                            Ice%dyn_trans_CSp, dirs%output_directory, Time_Init)
  !  IST%ice_transport_CSp => SIS_dyn_trans_transport_CS(Ice%dyn_trans_CSp)

    call SIS_slow_thermo_set_ptrs(Ice%slow_thermo_CSp, &
             transport_CSp=SIS_dyn_trans_transport_CS(Ice%dyn_trans_CSp), &
             sum_out_CSp=SIS_dyn_trans_sum_output_CS(Ice%dyn_trans_CSp))

  !   Initialize any tracers that will be handled via tracer flow control.
    call SIS_tracer_flow_control_init(Ice%Time, G, IG, param_file, Ice%SIS_tracer_flow_CSp, is_restart)
  endif

  call close_param_file(param_file)

  iceClock = mpp_clock_id( 'Ice', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  iceClock1 = mpp_clock_id( 'Ice: bot to top', flags=clock_flag_default, grain=CLOCK_ROUTINE )
  iceClock2 = mpp_clock_id( 'Ice: update slow (dn)', flags=clock_flag_default, grain=CLOCK_ROUTINE )
  iceClock3 = mpp_clock_id( 'Ice: update fast', flags=clock_flag_default, grain=CLOCK_ROUTINE )

  ! Initialize icebergs
  if (slow_ice_PE) then ; if (Ice%sCS%do_icebergs) then
     if( ASSOCIATED(G%Domain%maskmap)) then
       call icebergs_init(Ice%icebergs, G%Domain%niglobal, G%Domain%njglobal, &
               G%Domain%layout, G%Domain%io_layout, Ice%axes(1:2), &
               G%Domain%X_flags, G%Domain%Y_flags, time_type_to_real(Time_step_slow), &
               Time, G%geoLonBu(isc:iec,jsc:jec), G%geoLatBu(isc:iec,jsc:jec), &
               G%mask2dT(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
               G%dxCv(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), G%dyCu(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
               Ice%area,  G%cos_rot(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
               G%sin_rot(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), maskmap=G%Domain%maskmap )
     else
       call icebergs_init(Ice%icebergs, G%Domain%niglobal, G%Domain%njglobal, &
                G%Domain%layout, G%Domain%io_layout, Ice%axes(1:2), &
                G%Domain%X_flags, G%Domain%Y_flags, time_type_to_real(Time_step_slow), &
                Time, G%geoLonBu(isc:iec,jsc:jec), G%geoLatBu(isc:iec,jsc:jec), &
                G%mask2dT(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
                G%dxCv(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), G%dyCu(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
                Ice%area, G%cos_rot(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
                G%sin_rot(G%isc-1:G%iec+1,G%jsc-1:G%jec+1) )
     endif
  endif ; endif
  !nullify_domain perhaps could be called somewhere closer to set_domain 
  !but it should be callled after restore_state() otherwise it causes a restart mismatch
  call nullify_domain()

  ! Duplicate what is currently in IST in Ice%fCS and/or Ice%sCS.  This
  ! will be moved up later.
  if (fast_ice_PE) then
    Ice%fCS%IST => IST
    Ice%fCS%diag => IST%diag

    Ice%fCS%slab_ice = IST%slab_ice
    Ice%fCS%Cgrid_dyn = IST%Cgrid_dyn
    Ice%fCS%Rho_ice = IST%Rho_ice
    Ice%fCS%Rho_snow = IST%Rho_snow
    Ice%fCS%specified_ice = IST%specified_ice
    Ice%fCS%bounds_check = bounds_check
    Ice%fCS%debug = debug
  endif

  ! Duplicate what is currently in IST in Ice%fCS and/or Ice%sCS.  This
  ! will be moved up later.
  if (slow_ice_PE) then
    Ice%sCS%IST => IST
    Ice%sCS%diag => IST%diag

    Ice%sCS%slab_ice = IST%slab_ice
    Ice%sCS%Cgrid_dyn = IST%Cgrid_dyn
    Ice%sCS%Rho_ice = IST%Rho_ice
    Ice%sCS%Rho_snow = IST%Rho_snow
    Ice%sCS%specified_ice = IST%specified_ice
    Ice%sCS%bounds_check = bounds_check
    Ice%sCS%debug = debug
  endif

  ! Do any error checking here.
  if (debug) then
    call ice_grid_chksum(G, haloshift=2)
  endif

  call write_ice_statistics(IST, IST%Time, 0, G, IG, &
                            SIS_dyn_trans_sum_output_CS(Ice%dyn_trans_CSp))

  call callTree_leave("ice_model_init()")

end subroutine ice_model_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_end - writes the restart file and deallocates memory               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_end (Ice)
  type(ice_data_type), intent(inout) :: Ice

  type(ice_state_type), pointer :: IST => NULL()
  logical :: fast_ice_PE       ! If true, fast ice processes are handled on this PE.
  logical :: slow_ice_PE       ! If true, slow ice processes are handled on this PE.

  ! For now, both fast and slow processes occur on all sea-ice PEs.
  fast_ice_PE = .true. ; slow_ice_PE = .true.
  IST => Ice%Ice_state

  call ice_model_restart(Ice=Ice)

  !--- release memory ------------------------------------------------

  if (fast_ice_PE) then
    call SIS_fast_thermo_end(Ice%fast_thermo_CSp)
  endif

  if (slow_ice_PE) then
    call SIS_dyn_trans_end(Ice%dyn_trans_CSp)

    call SIS_slow_thermo_end(Ice%slow_thermo_CSp)

    call SIS2_ice_thm_end(IST%ice_thm_CSp, IST%ITV)

    ! End icebergs
    if (Ice%sCS%do_icebergs) call icebergs_end(Ice%icebergs)
  endif

  call SIS_hor_grid_end(Ice%G)
  call ice_grid_end(Ice%IG)
  call dealloc_Ice_arrays(Ice)

  call SIS_tracer_flow_control_end(Ice%SIS_tracer_flow_CSp)

  call dealloc_ocean_sfc_state(Ice%OSS)

  call dealloc_fast_ice_avg(Ice%FIA)

  call dealloc_ice_ocean_flux(Ice%IOF)

  call dealloc_IST_arrays(IST)
  deallocate(Ice%Ice_restart)

  if (Ice%Rad%add_diurnal_sw .or. Ice%Rad%do_sun_angle_for_alb) call astronomy_end

  call dealloc_ice_rad(Ice%Rad)

  call SIS_diag_mediator_end(IST%Time, IST%diag)

  deallocate(Ice%Ice_state)
  if (associated(Ice%fCS)) deallocate(Ice%fCS)
  if (associated(Ice%sCS)) deallocate(Ice%sCS)

end subroutine ice_model_end

end module ice_model_mod
