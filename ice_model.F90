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

use MOM_checksums,     only : chksum
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : fill_symmetric_edges, MOM_domains_init, clone_MOM_domain
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : open_param_file, close_param_file
use MOM_hor_index, only : hor_index_type, hor_index_init
use MOM_obsolete_params, only : obsolete_logical
use MOM_string_functions, only : uppercase
use MOM_time_manager, only : time_type, time_type_to_real
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

use ice_type_mod, only : ice_data_type, ice_state_type
use ice_type_mod, only : ice_model_restart, dealloc_ice_arrays, dealloc_IST_arrays
use ice_type_mod, only : ice_data_type_register_restarts, ice_state_register_restarts
use ice_type_mod, only : ice_diagnostics_init, ice_stock_pe, check_ice_model_nml
use ice_type_mod, only : ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
use ice_type_mod, only : ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
use ice_type_mod, only : lnd_ice_bnd_type_chksum, ice_data_type_chksum
use ice_type_mod, only : IST_chksum, Ice_public_type_chksum
use ice_type_mod, only : IST_bounds_check, Ice_public_type_bounds_check
use ice_utils_mod, only : post_avg, ice_grid_chksum
use SIS_hor_grid, only : SIS_hor_grid_type, set_hor_grid, SIS_hor_grid_end, set_first_direction
use SIS_fixed_initialization, only : SIS_initialize_fixed

use ice_grid, only : set_ice_grid, ice_grid_end, ice_grid_type
use ice_spec_mod, only : get_sea_surface

use SIS_tracer_registry, only : register_SIS_tracer, register_SIS_tracer_pair
use SIS_tracer_flow_control, only : SIS_call_tracer_register, SIS_tracer_flow_control_init
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_end

use ice_thm_mod,   only : slab_ice_optics, ice_thm_param, ice5lay_temp, TFI, CI
use SIS_slow_mod,  only : update_ice_model_slow, SIS_slow_register_restarts
use SIS_slow_mod,  only : SIS_slow_init, SIS_slow_end
use SIS2_ice_thm,  only : ice_temp_SIS2, ice_optics_SIS2, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,  only : get_SIS2_thermo_coefs, enth_from_TS, Temp_from_En_S, T_freeze
use ice_bergs,     only : icebergs, icebergs_run, icebergs_init, icebergs_end, icebergs_incr_mass

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

  call mpp_clock_begin(iceClock)
  call mpp_clock_begin(iceClock2)
  ! average fluxes from update_ice_model_fast
  call avg_top_quantities(Ice%Ice_state, Ice%G, Ice%IG)

  if (Ice%Ice_state%debug) then
    call Ice_public_type_chksum("Start update_ice_model_slow_dn", Ice)
  endif

  call set_ice_state_fluxes(Ice%Ice_state, Ice, Land_boundary, Ice%G, Ice%IG)

  call update_ice_model_slow(Ice%Ice_state, Ice%icebergs, Ice%G, Ice%IG)

  ! Set up the thermodynamic fluxes in the externally visible structure Ice.
  call set_ocean_top_fluxes(Ice, Ice%Ice_state, Ice%G, Ice%IG)

  if (Ice%Ice_state%debug) then
    call Ice_public_type_chksum("End update_ice_model_slow_dn", Ice)
  endif
  if (Ice%Ice_state%bounds_check) then
    call Ice_public_type_bounds_check(Ice, Ice%G, "End update_ice_slow")
  endif

  call mpp_clock_end(iceClock2)
  call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow_dn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! sum_top_quantities - sum fluxes for later use by ice/ocean slow physics.     !
!   Nothing here will be exposed to other modules.                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine sum_top_quantities ( Ice, IST, Atmos_boundary_fluxes, flux_u, flux_v, flux_t, flux_q, &
       flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif,&
       flux_lw, lprec, fprec, flux_lh, G, IG)
  type(ice_data_type),              intent(inout) :: Ice
  type(ice_state_type),             intent(inout) :: IST
  type(coupler_3d_bc_type),         intent(inout) :: Atmos_boundary_fluxes
  type(SIS_hor_grid_type),          intent(inout) :: G
  type(ice_grid_type),              intent(inout) :: IG
  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:IG%CatIce), intent(in) :: &
    flux_u, flux_v, flux_t, flux_q, flux_lw, lprec, fprec, flux_lh
  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:IG%CatIce), intent(in) :: &
    flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif
  real    :: Stefan ! The Stefan-Boltzmann constant in W m-2 K-4 as used for
                    ! strictly diagnostic purposes.
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, ncat
  integer :: ind, max_num_fields, next_index

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  i_off = LBOUND(Ice%albedo_vis_dir,1) - G%isc
  j_off = LBOUND(Ice%albedo_vis_dir,2) - G%jsc

  if (IST%num_tr_fluxes < 0) then
    ! Determine how many atmospheric boundary fluxes have been passed in, and
    ! set up both an indexing array for these and a space to take their average.
    ! This code is only exercised the first time that sum_top_quantities is called.
    IST%num_tr_fluxes = 0
    if (Atmos_boundary_fluxes%num_bcs > 0) then
      max_num_fields = 0
      do n=1,Atmos_boundary_fluxes%num_bcs
        IST%num_tr_fluxes = IST%num_tr_fluxes + Atmos_boundary_fluxes%bc(n)%num_fields
        max_num_fields = max(max_num_fields, Atmos_boundary_fluxes%bc(n)%num_fields)
      enddo
      if (IST%num_tr_fluxes > 0) then
        allocate(IST%tr_flux_top(SZI_(G), SZJ_(G), 0:IG%CatIce, IST%num_tr_fluxes))
        allocate(IST%tr_flux_ocn_top(SZI_(G), SZJ_(G), IST%num_tr_fluxes))
        IST%tr_flux_top(:,:,:,:) = 0.0 ; IST%tr_flux_ocn_top(:,:,:) = 0.0
        allocate(IST%tr_flux_index(max_num_fields, Atmos_boundary_fluxes%num_bcs))
        IST%tr_flux_index(:,:) = -1 ; next_index = 1
        do n=1,Atmos_boundary_fluxes%num_bcs ; do m=1,Atmos_boundary_fluxes%bc(n)%num_fields
          IST%tr_flux_index(m, n) = next_index ; next_index = next_index + 1
        enddo ; enddo
      endif
    endif
  endif

  if (IST%avg_count == 0) then
    ! zero_top_quantities - zero fluxes to begin summing in ice fast physics.
    IST%flux_u_top(:,:,:) = 0.0 ; IST%flux_v_top(:,:,:) = 0.0
    IST%lwdn(:,:) = 0.0 ; IST%swdn(:,:) = 0.0
    IST%flux_t_top(:,:,:) = 0.0 ; IST%flux_q_top(:,:,:) = 0.0
    IST%flux_lw_top(:,:,:) = 0.0 ; IST%flux_lh_top(:,:,:) = 0.0
    IST%flux_sw_nir_dir_top(:,:,:) = 0.0 ; IST%flux_sw_nir_dif_top(:,:,:) = 0.0
    IST%flux_sw_vis_dir_top(:,:,:) = 0.0 ; IST%flux_sw_vis_dif_top(:,:,:) = 0.0
    IST%lprec_top(:,:,:) = 0.0 ; IST%fprec_top(:,:,:) = 0.0
    if (IST%num_tr_fluxes > 0) IST%tr_flux_top(:,:,:,:) = 0.0
  endif

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,flux_u,flux_v,flux_t, &
!$OMP                                  flux_q,flux_sw_nir_dir,flux_sw_nir_dif,        &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,flux_lw,       &
!$OMP                                  lprec,fprec,flux_lh)
  do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
    IST%flux_u_top(i,j,k)  = IST%flux_u_top(i,j,k)  + flux_u(i,j,k)
    IST%flux_v_top(i,j,k)  = IST%flux_v_top(i,j,k)  + flux_v(i,j,k)
    IST%flux_t_top(i,j,k)  = IST%flux_t_top(i,j,k)  + flux_t(i,j,k)
    IST%flux_q_top(i,j,k)  = IST%flux_q_top(i,j,k)  + flux_q(i,j,k)
    IST%flux_sw_nir_dir_top(i,j,k) = IST%flux_sw_nir_dir_top(i,j,k) + flux_sw_nir_dir(i,j,k)
    IST%flux_sw_nir_dif_top(i,j,k) = IST%flux_sw_nir_dif_top(i,j,k) + flux_sw_nir_dif(i,j,k)
    IST%flux_sw_vis_dir_top(i,j,k) = IST%flux_sw_vis_dir_top(i,j,k) + flux_sw_vis_dir(i,j,k)
    IST%flux_sw_vis_dif_top(i,j,k) = IST%flux_sw_vis_dif_top(i,j,k) + flux_sw_vis_dif(i,j,k)
    IST%flux_lw_top(i,j,k) = IST%flux_lw_top(i,j,k) + flux_lw(i,j,k)
    IST%lprec_top(i,j,k)   = IST%lprec_top(i,j,k)   + lprec(i,j,k)
    IST%fprec_top(i,j,k)   = IST%fprec_top(i,j,k)   + fprec(i,j,k)
    IST%flux_lh_top(i,j,k) = IST%flux_lh_top(i,j,k) + flux_lh(i,j,k)
  enddo ; enddo ; enddo

  do n=1,Atmos_boundary_fluxes%num_bcs ; do m=1,Atmos_boundary_fluxes%bc(n)%num_fields
    ind = IST%tr_flux_index(m,n)
    if (ind < 1) call SIS_error(FATAL, "Bad boundary flux index in sum_top_quantities.")
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      IST%tr_flux_top(i,j,k,ind) = IST%tr_flux_top(i,j,k,ind) + &
            Atmos_boundary_fluxes%bc(n)%field(m)%values(i2,j2,k2)
    enddo ; enddo ; enddo
  enddo ; enddo

  if (IST%id_lwdn > 0) then
    Stefan = 5.6734e-8  ! Set the Stefan-Bolzmann constant, in W m-2 K-4.
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
      IST%lwdn(i,j) = IST%lwdn(i,j) + IST%part_size(i,j,k) * &
                           (flux_lw(i,j,k) + Stefan*IST%t_surf(i,j,k)**4)
    endif ; enddo ; enddo ; enddo
  endif

  if (IST%id_swdn > 0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,IST,Ice,i_off,j_off, &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,            &
!$OMP                                  flux_sw_nir_dir,flux_sw_nir_dif)            &
!$OMP                          private(i2,j2,k2)
    do j=jsc,jec ; do k=0,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      IST%swdn(i,j) = IST%swdn(i,j) + IST%part_size(i,j,k) * ( &
            (flux_sw_vis_dir(i,j,k)/(1-Ice%albedo_vis_dir(i2,j2,k2)) + &
             flux_sw_vis_dif(i,j,k)/(1-Ice%albedo_vis_dif(i2,j2,k2))) + &
            (flux_sw_nir_dir(i,j,k)/(1-Ice%albedo_nir_dir(i2,j2,k2)) + &
             flux_sw_nir_dif(i,j,k)/(1-Ice%albedo_nir_dif(i2,j2,k2))) )
    endif ; enddo ; enddo ; enddo
  endif

  IST%avg_count = IST%avg_count + 1

end subroutine sum_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! avg_top_quantities - time average fluxes for ice and ocean slow physics      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine avg_top_quantities(IST, G, IG)
  type(ice_state_type),    intent(inout) :: IST
  type(SIS_hor_grid_type), intent(in)    :: G
  type(ice_grid_type),     intent(in)    :: IG

  real    :: u, v, divid, sign
  integer :: i, j, k, m, n, isc, iec, jsc, jec, ncat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  !
  ! compute average fluxes
  !
  if (IST%avg_count == 0) call SIS_error(FATAL,'avg_top_quantities: '//&
       'no ocean model fluxes have been averaged')

  ! Rotate the stress from lat/lon to ocean coordinates and possibly change the
  ! sign to positive for downward fluxes of positive momentum.
  sign = 1.0 ; if (IST%atmos_winds) sign = -1.0
  divid = 1.0/real(IST%avg_count)

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,sign,divid,G) private(u,v)
  do j=jsc,jec
    do k=0,ncat ;  do i=isc,iec
      u = IST%flux_u_top(i,j,k) * (sign*divid)
      v = IST%flux_v_top(i,j,k) * (sign*divid)
      IST%flux_u_top(i,j,k) = u*G%cos_rot(i,j)-v*G%sin_rot(i,j) ! rotate stress from lat/lon
      IST%flux_v_top(i,j,k) = v*G%cos_rot(i,j)+u*G%sin_rot(i,j) ! to ocean coordinates
      IST%flux_t_top(i,j,k)  = IST%flux_t_top(i,j,k)  * divid
      IST%flux_q_top(i,j,k)  = IST%flux_q_top(i,j,k)  * divid
      IST%flux_sw_nir_dir_top(i,j,k) = IST%flux_sw_nir_dir_top(i,j,k) * divid
      IST%flux_sw_nir_dif_top(i,j,k) = IST%flux_sw_nir_dif_top(i,j,k) * divid
      IST%flux_sw_vis_dir_top(i,j,k) = IST%flux_sw_vis_dir_top(i,j,k) * divid
      IST%flux_sw_vis_dif_top(i,j,k) = IST%flux_sw_vis_dif_top(i,j,k) * divid
      IST%flux_lw_top(i,j,k) = IST%flux_lw_top(i,j,k) * divid
      IST%fprec_top(i,j,k)   = IST%fprec_top(i,j,k)   * divid
      IST%lprec_top(i,j,k)   = IST%lprec_top(i,j,k)   * divid
      IST%flux_lh_top(i,j,k) = IST%flux_lh_top(i,j,k) * divid
      ! Convert frost forming atop sea ice into frozen precip.
      if ((k>0) .and. (IST%flux_q_top(i,j,k) < 0.0)) then
        IST%fprec_top(i,j,k) = IST%fprec_top(i,j,k) - IST%flux_q_top(i,j,k)
        IST%flux_q_top(i,j,k) = 0.0
      endif
      do n=1,IST%num_tr_fluxes
        IST%tr_flux_top(i,j,k,n) = IST%tr_flux_top(i,j,k,n) * divid
      enddo
    enddo ; enddo
    do i=isc,iec
      IST%lwdn(i,j) = IST%lwdn(i,j)* divid
      IST%swdn(i,j) = IST%swdn(i,j)* divid
    enddo
  enddo
  call pass_vector(IST%flux_u_top, IST%flux_v_top, G%Domain, stagger=AGRID)

  !
  ! set count to zero and fluxes will be zeroed before the next sum
  !
  IST%avg_count = 0
end subroutine avg_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ice_state_fluxes copies the ice surface fluxes and any other fields into
!! the ice_state_type.
subroutine set_ice_state_fluxes(IST, Ice, LIB, G, IG)
  type(ice_state_type),         intent(inout) :: IST
  type(ice_data_type),          intent(in)    :: Ice
  type(land_ice_boundary_type), intent(in)    :: LIB
  type(SIS_hor_grid_type),      intent(in)    :: G
  type(ice_grid_type),          intent(in)    :: IG

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  ! Store liquid runoff and other fluxes from the land to the ice or ocean.
  i_off = LBOUND(LIB%runoff,1) - G%isc ; j_off = LBOUND(LIB%runoff,2) - G%jsc
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,LIB,i_off,j_off) &
!$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    IST%runoff(i,j)  = LIB%runoff(i2,j2)
    IST%calving(i,j) = LIB%calving(i2,j2)
    IST%runoff_hflx(i,j)  = LIB%runoff_hflx(i2,j2)
    IST%calving_hflx(i,j) = LIB%calving_hflx(i2,j2)
  enddo ; enddo

  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,i_off,j_off) &
!$OMP                          private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    IST%flux_u_ocn(i,j) = Ice%flux_u(i2,j2)
    IST%flux_v_ocn(i,j) = Ice%flux_v(i2,j2)
  enddo ; enddo

end subroutine set_ice_state_fluxes

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! post_flux_diagnostics - write out any diagnostics of surface fluxes.         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine post_flux_diagnostics(IST, G, IG, Idt_slow)
  type(ice_state_type),    intent(in) :: IST
  type(SIS_hor_grid_type), intent(in) :: G
  type(ice_grid_type),     intent(in) :: IG
  real,                    intent(in) :: Idt_slow

  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: tmp2d
  integer :: i, j, k, m, n, isc, iec, jsc, jec, ncat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  ! Flux diagnostics
  !
  if (IST%id_runoff>0) &
    call post_data(IST%id_runoff, IST%runoff, IST%diag)
  if (IST%id_calving>0) &
    call post_data(IST%id_calving, IST%calving, IST%diag)
  if (IST%id_runoff_hflx>0) &
    call post_data(IST%id_runoff_hflx, IST%runoff_hflx, IST%diag)
  if (IST%id_calving_hflx>0) &
    call post_data(IST%id_calving_hflx, IST%calving_hflx, IST%diag)
  if (IST%id_frazil>0) &
    call post_data(IST%id_frazil, IST%frazil*Idt_slow, IST%diag)
  if (IST%id_sh>0) call post_avg(IST%id_sh, IST%flux_t_top, IST%part_size, &
                                 IST%diag, G=G)
  if (IST%id_lh>0) call post_avg(IST%id_lh, IST%flux_lh_top, IST%part_size, &
                                 IST%diag, G=G)
  if (IST%id_evap>0) call post_avg(IST%id_evap, IST%flux_q_top, IST%part_size, &
                                 IST%diag, G=G)
  if (IST%id_sw>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,tmp2d,IST)
    do j=jsc,jec
      do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * ( &
              IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k) + &
              IST%flux_sw_nir_dir_top(i,j,k) + IST%flux_sw_nir_dif_top(i,j,k) )
      enddo ; enddo
    enddo
    call post_data(IST%id_sw, tmp2d, IST%diag)
  endif
  if (IST%id_lw>0) call post_avg(IST%id_lw, IST%flux_lw_top, &
                                 IST%part_size, IST%diag, G=G)
  if (IST%id_snofl>0) call post_avg(IST%id_snofl, IST%fprec_top, &
                                    IST%part_size, IST%diag, G=G)
  if (IST%id_rain>0) call post_avg(IST%id_rain, IST%lprec_top, &
                                   IST%part_size, IST%diag, G=G)
  if (IST%id_lwdn>0) call post_data(IST%id_lwdn, IST%lwdn, IST%diag)
  if (IST%id_swdn>0) call post_data(IST%id_swdn, IST%swdn, IST%diag)
  if (IST%id_sw_vis>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,tmp2d,IST)
    do j=jsc,jec
      do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * ( &
              IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k) )
      enddo ; enddo
    enddo
    call post_data(IST%id_sw_vis, tmp2d, IST%diag)
  endif
  if (IST%id_sw_nir_dir>0) call post_avg(IST%id_sw_nir_dir, IST%flux_sw_nir_dir_top, &
                             IST%part_size, IST%diag, G=G)
  if (IST%id_sw_nir_dif>0) call post_avg(IST%id_sw_nir_dif, IST%flux_sw_nir_dif_top, &
                             IST%part_size, IST%diag, G=G)
  if (IST%id_sw_vis_dir>0) call post_avg(IST%id_sw_vis_dir, IST%flux_sw_vis_dir_top, &
                             IST%part_size, IST%diag, G=G)
  if (IST%id_sw_vis_dif>0) call post_avg(IST%id_sw_vis_dif, IST%flux_sw_vis_dif_top, &
                             IST%part_size, IST%diag, G=G)

  if (IST%nudge_sea_ice .and. IST%id_fwnudge>0) then
    call post_data(IST%id_fwnudge, IST%melt_nudge, IST%diag)
  endif

end subroutine post_flux_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ocean_top_fluxes - Translate ice-bottom fluxes of heat, mass, salt, and  !
!   tracers from the ice model's internal state to the public ice data type    !
!   for use by the ocean model.                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ocean_top_fluxes(Ice, IST, G, IG)
  type(ice_data_type),     intent(inout) :: Ice
  type(ice_state_type),    intent(in)    :: IST
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(in)    :: IG

  real :: I_count
  integer :: i, j, k, isc, iec, jsc, jec, m, n
  integer :: i2, j2, i_off, j_off, ind, ncat, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce
  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc

  if (IST%debug) then
    call IST_chksum("Start set_ocean_top_fluxes", IST, G, IG)
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
  if (IST%do_icebergs) call icebergs_incr_mass(Ice%icebergs, Ice%mi(:,:)) ! Add icebergs mass in kg/m^2

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,Ice,IST,i_off,j_off) &
!$OMP                           private(i2,j2)
  do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
    Ice%part_size(i2,j2,k+1) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  if (IST%id_slp>0) call post_data(IST%id_slp, Ice%p_surf, IST%diag)

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,Ice,IST,i_off,j_off) &
!$OMP                           private(i2,j2)
  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
    Ice%flux_u(i2,j2) = IST%flux_u_ocn(i,j)
    Ice%flux_v(i2,j2) = IST%flux_v_ocn(i,j)
    Ice%flux_t(i2,j2) = IST%flux_t_ocn_top(i,j)
    Ice%flux_q(i2,j2) = IST%flux_q_ocn_top(i,j)
    Ice%flux_sw_vis_dir(i2,j2) = IST%flux_sw_vis_dir_ocn(i,j)
    Ice%flux_sw_vis_dif(i2,j2) = IST%flux_sw_vis_dif_ocn(i,j)
    Ice%flux_sw_nir_dir(i2,j2) = IST%flux_sw_nir_dir_ocn(i,j)
    Ice%flux_sw_nir_dif(i2,j2) = IST%flux_sw_nir_dif_ocn(i,j)
    Ice%flux_lw(i2,j2) = IST%flux_lw_ocn_top(i,j)
    Ice%flux_lh(i2,j2) = IST%flux_lh_ocn_top(i,j)
    Ice%fprec(i2,j2) = IST%fprec_ocn_top(i,j)
    Ice%lprec(i2,j2) = IST%lprec_ocn_top(i,j)
    Ice%runoff(i2,j2)  = IST%runoff(i,j)
    Ice%calving(i2,j2) = IST%calving(i,j)
    Ice%runoff_hflx(i2,j2)  = IST%runoff_hflx(i,j)
    Ice%calving_hflx(i2,j2) = IST%calving_hflx(i,j)
    Ice%flux_salt(i2,j2) = IST%flux_salt(i,j)

    if (IST%slp2ocean) then
      Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) - 1e5 ! SLP - 1 std. atmosphere, in Pa.
    else
      Ice%p_surf(i2,j2) = 0.0
    endif
    Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) + G%G_Earth*Ice%mi(i2,j2)
  enddo ; enddo
  if (IST%nudge_sea_ice) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
      Ice%lprec(i2,j2) = Ice%lprec(i2,j2) + IST%melt_nudge(i,j)
    enddo ; enddo
  endif

  do n=1,Ice%ocean_fluxes%num_bcs ; do m=1,Ice%ocean_fluxes%bc(n)%num_fields
    ind = IST%tr_flux_index(m,n)
    if (ind < 1) call SIS_error(FATAL, "Bad boundary flux index in set_ocean_top_fluxes.")
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off  ! Use these to correct for indexing differences.
        Ice%ocean_fluxes%bc(n)%field(m)%values(i2,j2) = IST%tr_flux_ocn_top(i,j,ind)
    enddo ; enddo
  enddo ; enddo

  if (IST%debug) then
    call IST_chksum("End set_ocean_top_fluxes", IST, G, IG)
    call Ice_public_type_chksum("End set_ocean_top_fluxes", Ice)
  endif

end subroutine set_ocean_top_fluxes

!
! Coupler interface to provide ocean surface data to atmosphere.
!
subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
  type(ocean_ice_boundary_type), intent(inout) :: Ocean_boundary
  type(ice_data_type),           intent(inout) :: Ice

  call mpp_clock_begin(iceClock)
  call mpp_clock_begin(iceClock1)
  call set_ice_surface_state(Ice, Ice%Ice_state, Ocean_boundary%t, Ocean_boundary%u, Ocean_boundary%v, &
                             Ocean_boundary%frazil, Ocean_boundary, Ice%G, Ice%IG, &
                             Ocean_boundary%s, Ocean_boundary%sea_level )
  call mpp_clock_end(iceClock1)
  call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow_up

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ice_surface_state - prepare surface state for atmosphere fast physics    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ice_surface_state(Ice, IST, t_surf_ice_bot, u_surf_ice_bot, v_surf_ice_bot, &
                                 frazil_ice_bot, OIB, G, IG, s_surf_ice_bot, sea_lev_ice_bot )
  type(ice_data_type),                     intent(inout) :: Ice
  type(ice_state_type),                    intent(inout) :: IST
  type(SIS_hor_grid_type),                 intent(inout) :: G
  type(ice_grid_type),                     intent(inout) :: IG
  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: t_surf_ice_bot, u_surf_ice_bot
  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: v_surf_ice_bot, frazil_ice_bot
  type(ocean_ice_boundary_type),           intent(inout) :: OIB
  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: s_surf_ice_bot, sea_lev_ice_bot

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: m_ice_tot
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: h_ice_input
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: icec, icec_obs
  real, dimension(SZI_(G),SZJ_(G)) :: u_nonsym, v_nonsym
  real, dimension(IG%NkIce) :: sw_abs_lay
  real :: u, v
  real :: area_pt
  real :: I_Nk
  real :: kg_H_Nk  ! The conversion factor from units of H to kg/m2 over Nk.
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  logical :: sent
  real :: H_to_m_ice     ! The specific volumes of ice and snow times the
  real :: H_to_m_snow    ! conversion factor from thickness units, in m H-1.
  real :: dt_slow        ! The ice thermodynamics time step
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc
  I_Nk = 1.0 / IG%NkIce ; kg_H_Nk = IG%H_to_kg_m2 * I_Nk
  dt_slow = time_type_to_real(IST%Time_step_slow)

  H_to_m_snow = IG%H_to_kg_m2 / IST%Rho_snow ; H_to_m_ice = IG%H_to_kg_m2 / IST%Rho_ice

  ! pass ocean state through ice on first partition
  if (.not. IST%specified_ice) then ! otherwise, already set by update_ice_model_slow
    IST%t_surf(isc:iec,jsc:jec,0) = t_surf_ice_bot(isc:iec,jsc:jec)
  endif

  if (IST%do_init) then
    call get_sea_surface(IST%Time, IST%t_surf(isc:iec,jsc:jec,0), IST%part_size(isc:iec,jsc:jec,0:1), &
                         h_ice_input )
    do j=jsc,jec ; do i=isc,iec
      IST%mH_ice(i,j,1) = h_ice_input(i,j)*(IST%Rho_ice*IG%kg_m2_to_H)
    enddo ; enddo

    !   Transfer ice to the correct thickness category.  If do_ridging=.false.,
    ! the first call to ice_redistribute has the same result.  At present, all
    ! tracers are initialized to their default values, and snow is set to 0,
    ! and so do not need to be updated here.
    if (IST%do_ridging) then
      do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,1) > IG%mH_cat_bound(1)) then
        do k=ncat,2,-1 ; if (IST%mH_ice(i,j,1) > IG%mH_cat_bound(k-1)) then
          IST%part_size(i,j,k) = IST%part_size(i,j,1)
          IST%part_size(i,j,1) = 0.0
          IST%mH_ice(i,j,k) = IST%mH_ice(i,j,1) ; IST%mH_ice(i,j,1) = 0.0
          !  IST%mH_snow(i,j,k) = IST%mH_snow(i,j,1) ; IST%mH_snow(i,j,1) = 0.0
          exit ! from k-loop
        endif ; enddo
      endif ; enddo ; enddo
    endif

    call pass_var(IST%part_size, G%Domain, complete=.true. )
    call pass_var(IST%mH_ice, G%Domain, complete=.true. )
    IST%do_init = .false.
  endif

  ! Any special first-time initialization must be completed before this point.
  IST%first_time = .false.

  if (IST%bounds_check) &
    call IST_bounds_check(IST, G, IG, "Start of set_ice_surface_state")

  if (IST%debug) then
    call IST_chksum("Start set_ice_surface_state", IST, G, IG)
    call Ice_public_type_chksum("Start set_ice_surface_state", Ice)
    call chksum(u_surf_ice_bot(isc:iec,jsc:jec), "Start IB2IT u_surf_ice_bot")
    call chksum(v_surf_ice_bot(isc:iec,jsc:jec), "Start IB2IT v_surf_ice_bot")
  endif

! Transfer the ocean state for extra tracer fluxes.
  do n=1,OIB%fields%num_bcs  ; do m=1,OIB%fields%bc(n)%num_fields
    Ice%ocean_fields%bc(n)%field(m)%values(:,:,1) = OIB%fields%bc(n)%field(m)%values
  enddo ; enddo
  m_ice_tot(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,IST,t_surf_ice_bot,          &
!$OMP                                  s_surf_ice_bot,frazil_ice_bot,sea_lev_ice_bot, &
!$OMP                                  ncat,m_ice_tot,Ice,i_off,j_off)                &
!$OMP                          private(i2,j2,k2)
  do j=jsc,jec
    do i=isc,iec
      IST%t_ocn(i,j) = t_surf_ice_bot(i,j) - T_0degC
      IST%s_surf(i,j) = s_surf_ice_bot(i,j)
      IST%frazil(i,j) = frazil_ice_bot(i,j)
      IST%sea_lev(i,j) = sea_lev_ice_bot(i,j)
    enddo

    do k=1,ncat ; do i=isc,iec
      IST%tmelt(i,j,k) = 0.0 ; IST%bmelt(i,j,k) = 0.0
      m_ice_tot(i,j) = m_ice_tot(i,j) + IST%mH_ice(i,j,k) * IST%part_size(i,j,k)
    enddo ; enddo

    do i=isc,iec
      if (m_ice_tot(i,j) > 0.0) then
        IST%bheat(i,j) = IST%kmelt*(IST%t_ocn(i,j) - T_Freeze(IST%s_surf(i,j), IST%ITV))
      else
        IST%bheat(i,j) = 0.0
      endif
    enddo
  enddo

  if (IST%slab_ice) then
    IST%sw_abs_sfc(:,:,:) = 0.0 ; IST%sw_abs_snow(:,:,:) = 0.0
    IST%sw_abs_ice(:,:,:,:) = 0.0 ; IST%sw_abs_ocn(:,:,:) = 0.0
    IST%sw_abs_int(:,:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Ice,i_off,j_off, &
!$OMP                                  H_to_m_snow,H_to_m_ice) &
!$OMP                          private(i2,j2,k2)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      call slab_ice_optics(IST%mH_snow(i,j,k)*H_to_m_snow, IST%mH_ice(i,j,k)*H_to_m_ice, &
               IST%t_surf(i,j,k)-T_0degC, T_Freeze(IST%s_surf(i,j),IST%ITV), &
               Ice%albedo(i2,j2,k2))

      Ice%albedo_vis_dir(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_vis_dif(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_nir_dir(i2,j2,k2) = Ice%albedo(i2,j2,k2)
      Ice%albedo_nir_dif(i2,j2,k2) = Ice%albedo(i2,j2,k2)
    endif ; enddo ; enddo ; enddo
  else
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Ice,G,IG,i_off,j_off, &
!$OMP                                  H_to_m_snow,H_to_m_ice) &
!$OMP                          private(i2,j2,k2,sw_abs_lay)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      call ice_optics_SIS2(IST%mH_snow(i,j,k)*H_to_m_snow, IST%mH_ice(i,j,k)*H_to_m_ice, &
               IST%t_surf(i,j,k)-T_0degC, T_Freeze(IST%s_surf(i,j),IST%ITV), IG%NkIce, &
               Ice%albedo_vis_dir(i2,j2,k2), Ice%albedo_vis_dif(i2,j2,k2), &
               Ice%albedo_nir_dir(i2,j2,k2), Ice%albedo_nir_dif(i2,j2,k2), &
               IST%sw_abs_sfc(i,j,k),  IST%sw_abs_snow(i,j,k), &
               sw_abs_lay, IST%sw_abs_ocn(i,j,k), IST%sw_abs_int(i,j,k), &
               IST%ice_thm_CSp, coszen_in=IST%coszen(i,j))

      do m=1,IG%NkIce ; IST%sw_abs_ice(i,j,k,m) = sw_abs_lay(m) ; enddo

      !Niki: Is the following correct for diagnostics?
      !   Probably this calculation of the "average" albedo should be replaced
      ! with a calculation that weights the averaging by the fraction of the
      ! shortwave radiation in each wavelength and orientation band.  However,
      ! since this is only used for diagnostic purposes, making this change
      ! might not be too urgent. -RWH
      Ice%albedo(i2,j2,k2)=(Ice%albedo_vis_dir(i2,j2,k2)+Ice%albedo_nir_dir(i2,j2,k2)&
                        +Ice%albedo_vis_dif(i2,j2,k2)+Ice%albedo_nir_dif(i2,j2,k2))/4

    endif ; enddo ; enddo ; enddo
  endif

  if (Ice%flux_uv_stagger == AGRID) then
    u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      u_nonsym(i,j) = u_surf_ice_bot(i,j) ; v_nonsym(i,j) = v_surf_ice_bot(i,j)
    enddo ; enddo
    call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, stagger=AGRID)

    if (associated(IST%u_ocn) .and. associated(IST%v_ocn)) then
      do J=jsc-1,jec ; do I=isc-1,iec
        IST%u_ocn(I,J) = 0.25*((u_nonsym(i,j) + u_nonsym(i+1,j+1)) + &
                               (u_nonsym(i+1,j) + u_nonsym(i,j+1)))
        IST%v_ocn(I,J) = 0.25*((v_nonsym(i,j) + v_nonsym(i+1,j+1)) + &
                               (v_nonsym(i+1,j) + v_nonsym(i,j+1)))
      enddo ; enddo
      call pass_vector(IST%u_ocn, IST%v_ocn, G%Domain, stagger=BGRID_NE)
    endif

    if (associated(IST%u_ocn_C) .and. associated(IST%v_ocn_C)) then
      do j=jsc,jec ; do I=isc-1,iec
        IST%u_ocn_C(I,j) = 0.5*(u_nonsym(i,j) + u_nonsym(i+1,j))
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        IST%v_ocn_C(i,J) = 0.5*(v_nonsym(i,j) + v_nonsym(i,j+1))
      enddo ; enddo
      call pass_vector(IST%u_ocn_C, IST%v_ocn_C, G%Domain, stagger=CGRID_NE)
    endif

  elseif (OIB%stagger == BGRID_NE) then
    if (IST%Cgrid_dyn) then
        u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
        do j=jsc,jec ; do i=isc,iec
          u_nonsym(i,j) = u_surf_ice_bot(i,j) ; v_nonsym(i,j) = v_surf_ice_bot(i,j)
        enddo ; enddo
        call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, stagger=BGRID_NE)

      do j=jsc,jec ; do I=isc-1,iec
        IST%u_ocn_C(I,j) = 0.5*(u_nonsym(I,J) + u_nonsym(I,J-1))
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        IST%v_ocn_C(i,J) = 0.5*(v_nonsym(I,J) + v_nonsym(I-1,J))
      enddo ; enddo
      call pass_vector(IST%u_ocn_C, IST%v_ocn_C, G%Domain, stagger=CGRID_NE)
    else
      do J=jsc,jec ; do I=isc,iec
        IST%u_ocn(I,J) = u_surf_ice_bot(I,J) ! need under-ice current
        IST%v_ocn(I,J) = v_surf_ice_bot(I,J) ! for water drag term
      enddo ; enddo
      if (G%symmetric) &
        call fill_symmetric_edges(IST%u_ocn, IST%v_ocn, G%Domain, stagger=BGRID_NE)

      call pass_vector(IST%u_ocn, IST%v_ocn, G%Domain, stagger=BGRID_NE)
    endif

  elseif (OIB%stagger == CGRID_NE) then
    if (IST%Cgrid_dyn) then
      do j=jsc,jec ; do I=isc,iec
        IST%u_ocn_C(I,j) = u_surf_ice_bot(I,j)
      enddo ; enddo
      do J=jsc,jec ; do i=isc,iec
        IST%v_ocn_C(i,J) = v_surf_ice_bot(I,j)
      enddo ; enddo
      if (G%symmetric) &
        call fill_symmetric_edges(IST%u_ocn_C, IST%v_ocn_C, G%Domain, stagger=CGRID_NE)

      call pass_vector(IST%u_ocn_C, IST%v_ocn_C, G%Domain, stagger=CGRID_NE)
    else
      u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
      do j=jsc,jec ; do i=isc,iec
        u_nonsym(I,j) = u_surf_ice_bot(I,j) ; v_nonsym(i,J) = v_surf_ice_bot(i,J)
      enddo ; enddo
      call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, stagger=CGRID_NE)
      do J=jsc-1,jec ; do I=isc-1,iec
        IST%u_ocn(I,J) = 0.5*(u_nonsym(I,j) + u_nonsym(I,j+1))
        IST%v_ocn(I,J) = 0.5*(v_nonsym(i,J) + v_nonsym(i+1,J))
      enddo ; enddo
      call pass_vector(IST%u_ocn, IST%v_ocn, G%Domain, stagger=BGRID_NE)
    endif
  else
    call SIS_error(FATAL, "set_ice_surface_state: Unrecognized OIB%stagger.")
  endif

  call pass_var(IST%sea_lev, G%Domain)

  if (IST%debug) then
    if (associated(IST%u_ocn) .and. associated(IST%v_ocn)) then
      call chksum(IST%u_ocn(isc:iec,jsc:jec), "Post-pass IST%u_ocn(0,0)")
      call chksum(IST%v_ocn(isc:iec,jsc:jec), "Post-pass IST%v_ocn(0,0)")
    endif
    if (associated(IST%u_ocn_C) .and. associated(IST%v_ocn_C)) then
      call chksum(IST%u_ocn_C(isc:iec,jsc:jec), "Post-pass IST%u_ocn_C(0,0)")
      call chksum(IST%v_ocn_C(isc:iec,jsc:jec), "Post-pass IST%v_ocn_C(0,0)")
    endif
  endif

  if (IST%bounds_check) &
    call IST_bounds_check(IST, G, IG, "Midpoint set_ice_surface_state")

  ! Copy the surface temperatures into the externally visible data type.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,ncat,i_off,j_off) &
!$OMP                          private(i2,j2,k2)
  do j=jsc,jec
    do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      Ice%s_surf(i2,j2) = IST%s_surf(i,j)
    enddo
    do k=0,ncat ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
      Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
    enddo ; enddo
  enddo

  ! put ocean and ice velocities into Ice%u_surf/v_surf on t-cells
  if (IST%Cgrid_dyn) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,G,i_off,j_off) &
!$OMP                          private(i2,j2)
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      if (G%mask2dT(i,j) > 0.5 ) then
        Ice%u_surf(i2,j2,1) = 0.5*(IST%u_ocn_C(I,j) + IST%u_ocn_C(I-1,j))
        Ice%v_surf(i2,j2,1) = 0.5*(IST%v_ocn_C(i,J) + IST%v_ocn_C(i,J-1))
        Ice%u_surf(i2,j2,2) = 0.5*(IST%u_ice_C(I,j) + IST%u_ice_C(I-1,j))
        Ice%v_surf(i2,j2,2) = 0.5*(IST%v_ice_C(i,J) + IST%v_ice_C(i,J-1))
      else
        Ice%u_surf(i2,j2,1) = 0.0 ; Ice%v_surf(i2,j2,1) = 0.0
        Ice%u_surf(i2,j2,2) = 0.0 ; Ice%v_surf(i2,j2,2) = 0.0
      endif
    enddo ; enddo
  else ! B-grid discretization.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST,Ice,G,i_off,j_off) &
!$OMP                          private(i2,j2)
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      if (G%mask2dT(i,j) > 0.5 ) then
        Ice%u_surf(i2,j2,1) = 0.25*((IST%u_ocn(I,J) + IST%u_ocn(I-1,J-1)) + &
                                    (IST%u_ocn(I,J-1) + IST%u_ocn(I-1,J)) )
        Ice%v_surf(i2,j2,1) = 0.25*((IST%v_ocn(I,J) + IST%v_ocn(I-1,J-1)) + &
                                    (IST%v_ocn(I,J-1) + IST%v_ocn(I-1,J)) )
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

  if (IST%debug) then
    call chksum(Ice%u_surf(:,:,1), "Intermed Ice%u_surf(1)")
    call chksum(Ice%v_surf(:,:,1), "Intermed Ice%v_surf(1)")
    call chksum(Ice%u_surf(:,:,2), "Intermed Ice%u_surf(2)")
    call chksum(Ice%v_surf(:,:,2), "Intermed Ice%v_surf(2)")
    call chksum(G%mask2dT(isc:iec,jsc:jec), "Intermed G%mask2dT")
    if (associated(IST%u_ocn_C) .and. associated(IST%v_ocn_C)) then
      call chksum(IST%u_ocn_C(isc:iec,jsc:jec), "Intermed IST%u_ocn_C(0,0)")
      call chksum(IST%v_ocn_C(isc:iec,jsc:jec), "Intermed IST%v_ocn_C(0,0)")
    endif
    if (associated(IST%u_ocn) .and. associated(IST%v_ocn)) then
      call chksum(IST%u_ocn(isc:iec,jsc:jec), "Intermed IST%u_ocn(0,0)")
      call chksum(IST%u_ocn(isc-1:iec-1,jsc:jec), "Intermed IST%u_ocn(-,0)")
      call chksum(IST%u_ocn(isc:iec,jsc-1:jec-1), "Intermed IST%u_ocn(0,-)")
      call chksum(IST%u_ocn(isc-1:iec-1,jsc-1:jec-1), "Intermed IST%u_ocn(-,-)")
      call chksum(IST%v_ocn(isc:iec,jsc:jec), "Intermed IST%v_ocn(0,0)")
      call chksum(IST%v_ocn(isc-1:iec-1,jsc:jec), "Intermed IST%v_ocn(-,0)")
      call chksum(IST%v_ocn(isc:iec,jsc-1:jec-1), "Intermed IST%v_ocn(0,-)")
      call chksum(IST%v_ocn(isc-1:iec-1,jsc-1:jec-1), "Intermed IST%v_ocn(-,-)")
    endif
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
  if (IST%debug) then
    do k2=1,ncat+1
      call chksum(Ice%u_surf(:,:,k2), "End Ice%u_surf(k2)")
      call chksum(Ice%v_surf(:,:,k2), "End Ice%v_surf(k2)")
    enddo
  endif

  if (IST%column_check) then
    IST%enth_prev(:,:,:) = 0.0
    IST%heat_in(:,:,:) = 0.0
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,k)>0.0) then
      IST%enth_prev(i,j,k) = (IST%mH_snow(i,j,k)*IG%H_to_kg_m2) * IST%enth_snow(i,j,k,1)
      do m=1,IG%NkIce
        IST%enth_prev(i,j,k) = IST%enth_prev(i,j,k) + &
                               (IST%mH_ice(i,j,k)*kg_H_Nk) * IST%enth_ice(i,j,k,m)
      enddo
    endif ; enddo ; enddo ; enddo
  endif
  !
  ! Pre-timestep diagnostics
  !
  call enable_SIS_averaging(real(time_type_to_real(IST%Time_step_slow)), IST%Time, IST%diag)
  if (IST%id_alb>0) call post_avg(IST%id_alb, Ice%albedo, &
                     IST%part_size(isc:iec,jsc:jec,:), IST%diag)
  if (IST%id_sst>0) call post_data(IST%id_sst, IST%t_ocn, IST%diag)
  if (IST%id_sss>0) call post_data(IST%id_sss, IST%s_surf, IST%diag)
  if (IST%id_ssh>0) call post_data(IST%id_ssh, IST%sea_lev, IST%diag)
  if (IST%Cgrid_dyn) then
    if (IST%id_uo>0) call post_data(IST%id_uo, IST%u_ocn_C, IST%diag)
    if (IST%id_vo>0) call post_data(IST%id_vo, IST%v_ocn_C, IST%diag)
  else
    if (IST%id_uo>0) call post_data(IST%id_uo, IST%u_ocn, IST%diag)
    if (IST%id_vo>0) call post_data(IST%id_vo, IST%v_ocn, IST%diag)
  endif
  if (IST%id_bheat>0) call post_data(IST%id_bheat, IST%bheat, IST%diag)
  call disable_SIS_averaging(IST%diag)

  if (IST%debug) then
    call IST_chksum("End set_ice_surface_state", IST, G, IG)
    call Ice_public_type_chksum("End set_ice_surface_state", Ice)
  endif

  if (IST%bounds_check) &
    call Ice_public_type_bounds_check(Ice, G, "End set_ice_surface_state")

end subroutine set_ice_surface_state

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! update_ice_model_fast - records fluxes (in Ice) and calculates ice temp. on  !
!                         (fast) atmospheric timestep (see coupler_main.f90)   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine update_ice_model_fast( Atmos_boundary, Ice )

  type(ice_data_type),              intent(inout) :: Ice
  type(atmos_ice_boundary_type),    intent(inout) :: Atmos_boundary

  call mpp_clock_begin(iceClock)
  call mpp_clock_begin(iceClock3)

  call do_update_ice_model_fast( Atmos_boundary, Ice, Ice%Ice_state, Ice%G, Ice%IG )

  call mpp_clock_end(iceClock3)
  call mpp_clock_end(iceClock)

end subroutine update_ice_model_fast

subroutine do_update_ice_model_fast( Atmos_boundary, Ice, IST, G, IG )

  type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
  type(ice_data_type),           intent(inout) :: Ice
  type(ice_state_type),          intent(inout) :: IST
  type(SIS_hor_grid_type),       intent(inout) :: G
  type(ice_grid_type),           intent(inout) :: IG

  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:IG%CatIce) :: &
    flux_t, flux_q, flux_lh, flux_lw, &
    flux_sw_nir_dir, flux_sw_nir_dif, &
    flux_sw_vis_dir, flux_sw_vis_dif, &
    flux_u, flux_v, lprec, fprec, &
    dhdt, &   ! The derivative of the upward sensible heat flux with the surface
              ! temperature in W m-2 K-1.
    dedt, &   ! The derivative of the sublimation rate with the surface
              ! temperature, in kg m-2 s-1 K-1.
    drdt      ! The derivative of the upward radiative heat flux with surface
              ! temperature (i.e. d(flux)/d(surf_temp)) in W m-2 K-1.
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: &
    diurnal_factor, cosz_alb
  real, dimension(SZI_(G), SZJ_(G)) :: tmp_diag
  real, dimension(0:IG%NkIce) :: T_col ! The temperature of a column of ice and snow in degC.
  real, dimension(IG%NkIce)   :: S_col ! The thermodynamic salinity of a column of ice, in g/kg.
  real, dimension(0:IG%NkIce) :: enth_col   ! The enthalpy of a column of snow and ice, in enth_unit (J/kg?).
  real, dimension(0:IG%NkIce) :: SW_abs_col
  real :: dt_fast, ts_new, dts, hf, hfd, latent
  real :: hf_0    ! The positive upward surface heat flux when T_sfc = 0 C, in W m-2.
  real :: dhf_dt  ! The deriviative of the upward surface heat flux with Ts, in W m-2 C-1.
  real :: gmt, time_since_ae, cosz, rrsun, fracday, fracday_dt_ice, fracday_day
  real :: rad, cosz_day, cosz_dt_ice, rrsun_day, rrsun_dt_ice
  real :: flux_sw ! sum over dir/dif vis/nir components
  real :: T_freeze_surf ! The freezing temperature at the surface salinity of
                        ! the ocean, in deg C.
  real :: T_freeze_ice_top ! The freezing temperature at the salinity of the
                        ! upper layer of the ice, in deg C.
  real :: LatHtFus       ! The latent heat of fusion of ice in J/kg.
  real :: LatHtVap       ! The latent heat of vaporization of water at 0C in J/kg.
  real :: H_to_m_ice     ! The specific volumes of ice and snow times the
  real :: H_to_m_snow    ! conversion factor from thickness units, in m H-1.
  type(time_type) :: Dt_ice
  logical :: sent
  integer :: i, j, k, m, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off, NkIce

  real :: tot_heat_in, enth_here, enth_imb, norm_enth_imb, SW_absorbed
  real :: enth_liq_0 ! The value of enthalpy for liquid fresh water at 0 C, in
                     ! enthalpy units (sometimes J kg-1).
  real :: enth_units ! A conversion factor from Joules kg-1 to enthalpy units.
  real :: I_enth_unit  ! The inverse of enth_units, in J kg-1 enth_unit-1.
  real :: I_Nk
  real :: kg_H_Nk  ! The conversion factor from units of H to kg/m2 over Nk.
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  i_off = LBOUND(Atmos_boundary%t_flux,1) - G%isc
  j_off = LBOUND(Atmos_boundary%t_flux,2) - G%jsc
  NkIce = IG%NkIce ; I_Nk = 1.0 / NkIce ; kg_H_Nk = IG%H_to_kg_m2 * I_Nk

  rad = acos(-1.)/180.

  IST%n_fast = IST%n_fast + 1

  if (IST%debug) then
    call IST_chksum("Start do_update_ice_model_fast", IST, G, IG)
    call Ice_public_type_chksum("Start do_update_ice_model_fast", Ice)
  endif
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Atmos_boundary,i_off, &
!$OMP                                  j_off,Ice,flux_u,flux_v,flux_t,flux_q,flux_lw, &
!$OMP                                  flux_sw_nir_dir,flux_sw_nir_dif,               &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,               &
!$OMP                                  lprec,fprec,dhdt,dedt,drdt        )            &
!$OMP                          private(i2,j2,k2)
  do j=jsc,jec
    do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off
      IST%coszen(i,j) = Atmos_boundary%coszen(i2,j2,1)
      Ice%p_surf(i2,j2) = Atmos_boundary%p(i2,j2,1)
    enddo
    !   Set up local copies of fluxes.  The Atmos_boundary arrays may have
    ! different index conventions than are used internally in this component.
    do k=0,ncat ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      flux_u(i,j,k)  = Atmos_boundary%u_flux(i2,j2,k2)
      flux_v(i,j,k)  = Atmos_boundary%v_flux(i2,j2,k2)
      flux_t(i,j,k)  = Atmos_boundary%t_flux(i2,j2,k2)
      flux_q(i,j,k)  = Atmos_boundary%q_flux(i2,j2,k2)
      flux_lw(i,j,k) = Atmos_boundary%lw_flux(i2,j2,k2)
      flux_sw_nir_dir(i,j,k) = Atmos_boundary%sw_flux_nir_dir(i2,j2,k2)
      flux_sw_nir_dif(i,j,k) = Atmos_boundary%sw_flux_nir_dif(i2,j2,k2)
      flux_sw_vis_dir(i,j,k) = Atmos_boundary%sw_flux_vis_dir(i2,j2,k2)
      flux_sw_vis_dif(i,j,k) = Atmos_boundary%sw_flux_vis_dif(i2,j2,k2)
      lprec(i,j,k)   = Atmos_boundary%lprec(i2,j2,k2)
      fprec(i,j,k)   = Atmos_boundary%fprec(i2,j2,k2)
      dhdt(i,j,k) = Atmos_boundary%dhdt(i2,j2,k2)
      dedt(i,j,k) = Atmos_boundary%dedt(i2,j2,k2)
      drdt(i,j,k) = Atmos_boundary%drdt(i2,j2,k2)
    enddo ; enddo
  enddo

  if (IST%add_diurnal_sw .or. IST%do_sun_angle_for_alb) then
!---------------------------------------------------------------------
!    extract time of day (gmt) from time_type variable time with
!    function universal_time.
!---------------------------------------------------------------------
    gmt = universal_time(IST%Time)
!---------------------------------------------------------------------
!    extract the time of year relative to the northern hemisphere
!    autumnal equinox (time_since_ae) from time_type variable
!    time using the function orbital_time.
!---------------------------------------------------------------------
    time_since_ae = orbital_time(IST%Time)
!--------------------------------------------------------------------
!    call diurnal_solar_2d to calculate astronomy fields
!    convert G%geoLonT and G%geoLatT to radians
!    Per Rick Hemler:
!      call daily_mean_solar to get cosz (over a day)
!      call diurnal_solar with dtime=Dt_ice to get cosz over Dt_ice
!      diurnal_factor = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice/cosz_day*fracday_day*rrsun_day
!--------------------------------------------------------------------
    Dt_ice = IST%Time_step_fast
  endif
  if (IST%add_diurnal_sw) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,rad,IST,Dt_ice,time_since_ae, &
!$OMP                                  diurnal_factor) &
!$OMP                          private(cosz_dt_ice,fracday_dt_ice,rrsun_dt_ice, &
!$OMP                                  fracday_day,cosz_day,rrsun_day)
    do j=jsc,jec ; do i=isc,iec
      call diurnal_solar(G%geoLatT(i,j)*rad, G%geoLonT(i,j)*rad, IST%Time, cosz=cosz_dt_ice, &
                         fracday=fracday_dt_ice, rrsun=rrsun_dt_ice, dt_time=Dt_ice)
      call daily_mean_solar (G%geoLatT(i,j)*rad, time_since_ae, cosz_day, fracday_day, rrsun_day)
      diurnal_factor(i,j) = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice /   &
                                   max(1e-30, cosz_day*fracday_day*rrsun_day)
    enddo ; enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,flux_sw_nir_dir, &
!$OMP                                  flux_sw_nir_dif,flux_sw_vis_dir,      &
!$OMP                                  flux_sw_vis_dif,diurnal_factor)
    do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
      flux_sw_nir_dir(i,j,k) = flux_sw_nir_dir(i,j,k) * diurnal_factor(i,j)
      flux_sw_nir_dif(i,j,k) = flux_sw_nir_dif(i,j,k) * diurnal_factor(i,j)
      flux_sw_vis_dir(i,j,k) = flux_sw_vis_dir(i,j,k) * diurnal_factor(i,j)
      flux_sw_vis_dif(i,j,k) = flux_sw_vis_dif(i,j,k) * diurnal_factor(i,j)
    enddo ; enddo ; enddo
  endif

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             Latent_fusion=LatHtFus, Latent_vapor=LatHtVap)

  do j=jsc,jec ; do i=isc,iec
    flux_lh(i,j,0) = LatHtVap * flux_q(i,j,0)
  enddo ; enddo

  !
  ! implicit update of ice surface temperature
  !
  dt_fast = time_type_to_real(IST%Time_step_fast)

  if (IST%SIS1_5L_thermo) then
    if (NkIce /= 4) call SIS_error(FATAL, "SIS1_5L_thermodynamics requires that NK_ICE=4.")

    H_to_m_snow = IG%H_to_kg_m2 / IST%Rho_snow
    H_to_m_ice  = IG%H_to_kg_m2 / IST%Rho_ice
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,NkIce,IST,S_col,dhdt,dedt,  &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,flux_sw_nir_dir, &
!$OMP                                  flux_sw_nir_dif,drdt,flux_t,flux_q,flux_lw,      &
!$OMP                                  H_to_m_snow,H_to_m_ice,dt_fast,flux_lh,LatHtFus,LatHtVap) &
!$OMP                          private(T_Freeze_surf,latent,T_col,flux_sw,hfd,hf, &
!$OMP                                  ts_new,dts)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j), IST%ITV)
      if (IST%mH_ice(i,j,k) > 0.0) then
        T_col(0) = Temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
        do m=1,NkIce ; T_col(m) = Temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV) ; enddo

        if (IST%slab_ice) then
          latent         = LatHtVap+LatHtFus
        elseif (IST%mH_snow(i,j,k)>0.0) then
          latent         = LatHtVap + (LatHtFus-CI*T_col(0))
        else
          latent         = LatHtVap + LatHtFus*(1-TFI/T_col(1))
        endif
        flux_sw = (flux_sw_vis_dir(i,j,k) + flux_sw_vis_dif(i,j,k)) + &
                  (flux_sw_nir_dir(i,j,k) + flux_sw_nir_dif(i,j,k))
        hfd = dhdt(i,j,k) + dedt(i,j,k)*latent + drdt(i,j,k)
        hf  = flux_t(i,j,k) + flux_q(i,j,k)*latent - flux_lw(i,j,k)   &
              - IST%sw_abs_sfc(i,j,k)*flux_sw - hfd*(IST%t_surf(i,j,k)-T_0degC)
        !   This call updates the snow and ice temperatures and accumulates the
        ! surface and bottom melting/freezing energy.  The ice and snow do not
        ! actually lose or gain any mass from freezing or melting.
        call ice5lay_temp(IST%mH_snow(i,j,k)*H_to_m_snow, T_col(0), &
                          IST%mH_ice(i,j,k)*H_to_m_ice,    &
                          T_col(1), T_col(2), T_col(3),    &
                          T_col(4), ts_new, hf, hfd,       &
                          IST%sw_abs_snow(i,j,k)*flux_sw,  &
                          IST%sw_abs_ice(i,j,k,1)*flux_sw, &
                          IST%sw_abs_ice(i,j,k,2)*flux_sw, &
                          IST%sw_abs_ice(i,j,k,3)*flux_sw, &
                          IST%sw_abs_ice(i,j,k,4)*flux_sw, &
                          T_Freeze_surf, IST%bheat(i,j), dt_fast, &
                          IST%tmelt(i,j,k), IST%bmelt(i,j,k))
        dts                = ts_new - (IST%t_surf(i,j,k)-T_0degC)
        IST%t_surf(i,j,k)  = IST%t_surf(i,j,k) + dts
        flux_t(i,j,k)  = flux_t(i,j,k)  + dts * dhdt(i,j,k)
        flux_q(i,j,k)  = flux_q(i,j,k)  + dts * dedt(i,j,k)
        flux_lw(i,j,k) = flux_lw(i,j,k) - dts * drdt(i,j,k)
        flux_lh(i,j,k) = latent * flux_q(i,j,k)

        IST%enth_snow(i,j,k,1) = enth_from_TS(T_col(0), 0.0, IST%ITV)
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_from_TS(T_col(m), S_col(m), IST%ITV) ; enddo
      else ! IST%mH_ice <= 0
        flux_lh(i,j,k) = LatHtVap * flux_q(i,j,k)
      endif
    enddo ; enddo ; enddo
  else
    enth_liq_0 = Enth_from_TS(0.0, 0.0, IST%ITV) ; I_enth_unit = 1.0 / enth_units

    T_freeze_ice_top = T_Freeze(S_col(1), IST%ITV)
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,NkIce,IST,dhdt,dedt,drdt,   &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,flux_sw_nir_dir, &
!$OMP                                  flux_sw_nir_dif,flux_t,flux_q,flux_lw,enth_liq_0,&
!$OMP                                  dt_fast,flux_lh,I_enth_unit,G,S_col,kg_H_Nk,     &
!$OMP                                  enth_units,LatHtFus,LatHtVap,IG)                 &
!$OMP                          private(T_Freeze_surf,latent,enth_col,flux_sw,dhf_dt,    &
!$OMP                                  hf_0,ts_new,dts,SW_abs_col,SW_absorbed,enth_here,&
!$OMP                                  tot_heat_in,enth_imb,norm_enth_imb     )
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j), IST%ITV)
      if (IST%mH_ice(i,j,k) > 0.0) then
        enth_col(0) = IST%enth_snow(i,j,k,1)
        do m=1,NkIce ; enth_col(m) = IST%enth_ice(i,j,k,m) ; enddo

        ! In the case of sublimation of either snow or ice, the vapor is at 0 C.
        ! If the vapor should be at a different temperature, a correction would be
        ! made here.
        if (IST%slab_ice) then
          latent = LatHtVap + LatHtFus
        elseif (IST%mH_snow(i,j,k)>0.0) then
          latent = LatHtVap + (enth_liq_0 - IST%enth_snow(i,j,k,1)) * I_enth_unit
        else
          latent = LatHtVap + (enth_liq_0 - IST%enth_ice(i,j,k,1)) * I_enth_unit
        endif
        flux_sw = (flux_sw_vis_dir(i,j,k) + flux_sw_vis_dif(i,j,k)) + &
                  (flux_sw_nir_dir(i,j,k) + flux_sw_nir_dif(i,j,k))

        dhf_dt = (dhdt(i,j,k) + dedt(i,j,k)*latent) + drdt(i,j,k)
        hf_0 = ((flux_t(i,j,k) + flux_q(i,j,k)*latent) - &
                (flux_lw(i,j,k) + IST%sw_abs_sfc(i,j,k)*flux_sw)) - &
               dhf_dt * (IST%t_surf(i,j,k)-T_0degC)

        SW_abs_col(0) = IST%sw_abs_snow(i,j,k)*flux_sw
        do m=1,NkIce ; SW_abs_col(m) = IST%sw_abs_ice(i,j,k,m)*flux_sw ; enddo

        !   This call updates the snow and ice temperatures and accumulates the
        ! surface and bottom melting/freezing energy.  The ice and snow do not
        ! actually lose or gain any mass from freezing or melting.
        call ice_temp_SIS2(IST%mH_snow(i,j,k)*IG%H_to_kg_m2, IST%mH_ice(i,j,k)*IG%H_to_kg_m2, &
                          enth_col, S_col, hf_0, dhf_dt, SW_abs_col, &
                          T_Freeze_surf, IST%bheat(i,j), ts_new, &
                          dt_fast, NkIce, IST%tmelt(i,j,k), IST%bmelt(i,j,k), &
                          IST%ice_thm_CSp, IST%ITV, IST%column_check)
        IST%enth_snow(i,j,k,1) = enth_col(0)
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_col(m) ; enddo

        dts               = ts_new - (IST%t_surf(i,j,k)-T_0degC)
        IST%t_surf(i,j,k) = IST%t_surf(i,j,k) + dts
        flux_t(i,j,k)  = flux_t(i,j,k)  + dts * dhdt(i,j,k)
        flux_q(i,j,k)  = flux_q(i,j,k)  + dts * dedt(i,j,k)
        flux_lw(i,j,k) = flux_lw(i,j,k) - dts * drdt(i,j,k)
        flux_lh(i,j,k) = latent * flux_q(i,j,k)

        if (IST%column_check) then
          SW_absorbed = SW_abs_col(0)
          do m=1,NkIce ; SW_absorbed = SW_absorbed + SW_abs_col(m) ; enddo
          IST%heat_in(i,j,k) = IST%heat_in(i,j,k) + dt_fast * &
            ((flux_lw(i,j,k) + IST%sw_abs_sfc(i,j,k)*flux_sw) + SW_absorbed + &
             IST%bheat(i,j) - (flux_t(i,j,k) + flux_lh(i,j,k)))

          enth_here = (IG%H_to_kg_m2*IST%mH_snow(i,j,k)) * enth_col(0)
          do m=1,NkIce
            enth_here = enth_here + (IST%mH_ice(i,j,k)*kg_H_Nk) * enth_col(m)
          enddo
          tot_heat_in = enth_units * (IST%heat_in(i,j,k) - &
                                      (IST%bmelt(i,j,k) + IST%tmelt(i,j,k)))
          enth_imb = enth_here - (IST%enth_prev(i,j,k) + tot_heat_in)
          if (abs(enth_imb) > IST%imb_tol * (abs(enth_here) + &
                    abs(IST%enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
            norm_enth_imb = enth_imb / (abs(enth_here) + &
                    abs(IST%enth_prev(i,j,k)) + abs(tot_heat_in))
            enth_imb = enth_here - (IST%enth_prev(i,j,k) + tot_heat_in)
          endif
        endif

      else ! IST%mH_ice <= 0
        flux_lh(i,j,k) = LatHtVap * flux_q(i,j,k)
      endif
    enddo ; enddo ; enddo

  endif

  ! This routine works on the boundary exchange state.
  call compute_ocean_roughness (Ice%ocean_pt, Atmos_boundary%u_star(:,:,1), Ice%rough_mom(:,:,1), &
                                Ice%rough_heat(:,:,1), Ice%rough_moist(:,:,1)  )

  ! This routine works on the boundary exchange state.
  if (IST%do_sun_angle_for_alb) then
    call diurnal_solar(G%geoLatT(isc:iec,jsc:jec)*rad, G%geoLonT(isc:iec,jsc:jec)*rad, &
                 IST%time, cosz=cosz_alb, fracday=diurnal_factor, rrsun=rrsun_dt_ice, dt_time=Dt_ice)  !diurnal_factor as dummy
    call compute_ocean_albedo(Ice%ocean_pt, cosz_alb(:,:), Ice%albedo_vis_dir(:,:,1),&
                              Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                              Ice%albedo_nir_dif(:,:,1), rad*G%geoLatT(isc:iec,jsc:jec) )
  else
    call compute_ocean_albedo(Ice%ocean_pt, IST%coszen(isc:iec,jsc:jec), Ice%albedo_vis_dir(:,:,1),&
                              Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                              Ice%albedo_nir_dif(:,:,1), rad*G%geoLatT(isc:iec,jsc:jec) )
  endif

  call sum_top_quantities(Ice, IST, Atmos_boundary%fluxes, flux_u, flux_v, flux_t, &
    flux_q, flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif, &
    flux_lw, lprec, fprec, flux_lh, G, IG )

  IST%Time = IST%Time + IST%Time_step_fast ! advance time
  Ice%Time = IST%Time

  ! Copy the surface temperatures into the externally visible data type.
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%s_surf(i2,j2) = IST%s_surf(i,j)
  enddo ; enddo
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
  enddo ; enddo ; enddo

  call enable_SIS_averaging(dt_fast, IST%Time, IST%diag)
  if (IST%id_alb_vis_dir>0) call post_avg(IST%id_alb_vis_dir, Ice%albedo_vis_dir, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag)
  if (IST%id_alb_vis_dif>0) call post_avg(IST%id_alb_vis_dif, Ice%albedo_vis_dif, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag)
  if (IST%id_alb_nir_dir>0) call post_avg(IST%id_alb_nir_dir, Ice%albedo_nir_dir, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag)
  if (IST%id_alb_nir_dif>0) call post_avg(IST%id_alb_nir_dif, Ice%albedo_nir_dif, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag)

  if (IST%id_sw_abs_sfc>0) call post_avg(IST%id_sw_abs_sfc, IST%sw_abs_sfc, &
                                   IST%part_size(:,:,1:), IST%diag, G=G)
  if (IST%id_sw_abs_snow>0) call post_avg(IST%id_sw_abs_snow, IST%sw_abs_snow, &
                                   IST%part_size(:,:,1:), IST%diag, G=G)
  do m=1,NkIce
    if (IST%id_sw_abs_ice(m)>0) call post_avg(IST%id_sw_abs_ice(m), IST%sw_abs_ice(:,:,:,m), &
                                     IST%part_size(:,:,1:), IST%diag, G=G)
  enddo
  if (IST%id_sw_abs_ocn>0) call post_avg(IST%id_sw_abs_ocn, IST%sw_abs_ocn, &
                                   IST%part_size(:,:,1:), IST%diag, G=G)

  if (IST%id_sw_pen>0) then
    tmp_diag(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,tmp_diag)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * &
                     (IST%sw_abs_ocn(i,j,k) + IST%sw_abs_int(i,j,k))
    enddo ; enddo ; enddo
    call post_data(IST%id_sw_pen, tmp_diag, IST%diag)
  endif

  if (IST%id_coszen>0) call post_data(IST%id_coszen, IST%coszen, IST%diag)
  call disable_SIS_averaging(IST%diag)

  if (IST%debug) then
    call IST_chksum("End do_update_ice_model_fast", IST, G, IG)
    call Ice_public_type_chksum("End do_update_ice_model_fast", Ice)
  endif

  if (IST%bounds_check) &
    call Ice_public_type_bounds_check(Ice, G, "End update_ice_fast")
  if (IST%bounds_check) &
    call IST_bounds_check(IST, G, IG, "End of update_ice_fast")

end subroutine do_update_ice_model_fast

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

  ! Parameters that properly belong exclusively to ice_thm.
  real :: k_snow         ! snow conductivity (W/mK)
  real :: h_lo_lim       ! The min ice thickness for temp. calc, in m.
  real :: Time_unit      ! The time unit in seconds for ICE_STATS_INTERVAL.
  real :: H_to_kg_m2_tmp ! A temporary variable for holding the intended value
                         ! of the thickness to mass-per-unit-area conversion
                         ! factor.
  real :: enth_unit      ! A conversion factor for enthalpy from Joules kg-1.
  real :: massless_ice_enth, massless_snow_enth, massless_ice_salin
  real :: H_rescale_ice, H_rescale_snow ! Rescaling factors to account for
                         ! differences in thickness units between the current
                         ! model run and the input files.
  real, allocatable, target, dimension(:,:,:,:) :: t_ice_tmp, sal_ice_tmp
  real, allocatable, target, dimension(:,:,:) :: t_snow_tmp
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  real :: g_Earth !   The gravitational acceleration in m s-2.
  integer :: idr, id_sal
  integer :: write_geom
  logical :: test_grid_copy = .false.
  logical :: write_geom_files  ! If true, write out the grid geometry files.
  logical :: symmetric         ! If true, use symmetric memory allocation.
  logical :: global_indexing   ! If true use global horizontal index values instead
                               ! of having the data domain on each processor start at 1.
  integer :: first_direction   ! An integer that indicates which direction is to be
                               ! updated first in directionally split parts of the
                               ! calculation.  This can be altered during the course
                               ! of the run via calls to set_first_direction.
  logical :: read_aux_restart
  logical :: is_restart = .false.
  character(len=16)  :: stagger, dflt_stagger

  if (associated(Ice%Ice_state)) then
    call SIS_error(WARNING, "ice_model_init called with an associated "// &
                    "Ice%Ice_state structure. Model is already initialized.")
    return
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
  call get_param(param_file, mod, "SIS1_5L_THERMODYNAMICS", IST%SIS1_5L_thermo, &
                 "If true, use the thermodynamic calculations inhereted \n"//&
                 "from the SIS1 5 layer. Otherwise, use the newer SIS2 version.", &
                 default=.false.)
  
  call obsolete_logical(param_file, "INTERSPERSED_ICE_THERMO", warning_val=.false.)
  call get_param(param_file, mod, "AREA_WEIGHTED_STRESSES", IST%area_wtd_stress, &
                 "If true, use wind stresses that are weighted by the ice \n"//&
                 "areas in the neighboring cells.  The default (true) is \n"//&
                 "probably the right behavior, and this option will be \n"//&
                 "obsoleted as soon as it is verified to work properly.", &
                 default=.true.)

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
  IST%flux_uv_stagger = Ice%flux_uv_stagger

  call get_param(param_file, mod, "DT_ICE_DYNAMICS", IST%dt_ice_dyn, &
                 "The time step used for the slow ice dynamics, including \n"//&
                 "stepping the continuity equation and interactions \n"//&
                 "between the ice mass field and velocities.  If 0 or \n"//&
                 "negative the coupling time step will be used.", &
                 units="seconds", default=-1.0)

  call get_param(param_file, mod, "RHO_OCEAN", IST%Rho_ocean, &
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
  call get_param(param_file, mod, "ICE_KMELT", IST%kmelt, &
                 "A constant giving the proportionality of the ocean/ice \n"//&
                 "base heat flux to the tempature difference, given by \n"//&
                 "the product of the heat capacity per unit volume of sea \n"//&
                 "water times a molecular diffusive piston velocity.", &
                 units="W m-2 K-1", default=6e-5*4e6)
  call get_param(param_file, mod, "SNOW_CONDUCT", k_snow, &
                 "The conductivity of heat in snow.", units="W m-1 K-1", &
                 default=0.31)
  call get_param(param_file, mod, "COLUMN_CHECK", IST%column_check, &
                 "If true, add code to allow debugging of conservation \n"//&
                 "column-by-column.  This does not change answers, but \n"//&
                 "can increase model run time.", default=.false.)
  call get_param(param_file, mod, "IMBALANCE_TOLERANCE", IST%imb_tol, &
                 "The tolerance for imbalances to be flagged by COLUMN_CHECK.", &
                 units="nondim", default=1.0e-9)
  call get_param(param_file, mod, "ICE_BOUNDS_CHECK", IST%bounds_check, &
                 "If true, periodically check the values of ice and snow \n"//&
                 "temperatures and thicknesses to ensure that they are \n"//&
                 "sensible, and issue warnings if they are not.  This \n"//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)
  call get_param(param_file, mod, "DEBUG", IST%debug, &
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

  call get_param(param_file, mod, "ICE_SEES_ATMOS_WINDS", IST%atmos_winds, &
                 "If true, the sea ice is being given wind stresses with \n"//&
                 "the atmospheric sign convention, and need to have their \n"//&
                 "sign changed.", default=.true.)
  call get_param(param_file, mod, "ICE_BULK_SALINITY", IST%ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "g/kg", default=4.0)
  call get_param(param_file, mod, "ICE_RELATIVE_SALINITY", IST%ice_rel_salin, &
                 "The initial salinity of sea ice as a fraction of the \n"//&
                 "salinity of the seawater from which it formed.", &
                 units = "nondim", default=0.0)
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


  call get_param(param_file, mod, "NUDGE_SEA_ICE", IST%nudge_sea_ice, &
                 "If true, constrain the sea ice concentrations using observations.", &
                 default=.false.)
  if (IST%nudge_sea_ice) then
    call get_param(param_file, mod, "NUDGE_SEA_ICE_RATE", IST%nudge_sea_ice_rate, &
                 "The rate of cooling of ice-free water that should be ice \n"//&
                 "covered in order to constrained the ice concentration to \n"//&
                 "track observations.  A suggested value is ~10000 W m-2.", &
                 units = "W m-2", default=0.0)
    call get_param(param_file, mod, "NUDGE_SEA_ICE_TOLERANCE", IST%nudge_conc_tol, &
                 "The tolerance for mismatch in the sea ice concentations \n"//&
                 "before nudging begins to be applied.  Values of order 0.1\n"//&
                 "should work well.", units = "nondim", default=0.0)
    call get_param(param_file, mod, "NUDGE_SEA_ICE_STABILITY", IST%nudge_stab_fac, &
                 "A factor that determines whether the buoyancy flux \n"//&
                 "associated with the sea ice nudging of warm water includes \n"//&
                 "a freshwater flux so as to be destabilizing on net (<1), \n"//&
                 "stabilizing (>1), or neutral (=1).", units="nondim", default=1.0)
  endif

  call get_param(param_file, mod, "APPLY_SLP_TO_OCEAN", IST%slp2ocean, &
                 "If true, apply the atmospheric sea level pressure to \n"//&
                 "the ocean.", default=.false.)
  call get_param(param_file, mod, "MIN_H_FOR_TEMP_CALC", h_lo_lim, &
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
  call get_param(param_file, mod, "DO_RIDGING", IST%do_ridging, &
                 "If true, call the ridging routines.", default=.false.)
  if (.not.IST%SIS1_5L_thermo) then
    call get_param(param_file, mod, "SIS2_FILLING_FRAZIL", IST%filling_frazil, &
                 "If true, apply frazil to fill as many categories as \n"//&
                 "possible to fill in a uniform (minimum) amount of ice \n"//&
                 "in all the thinnest categories. Otherwise the frazil is \n"//&
                 "always assigned to a single category.", default=.false.) !###CHANGE DEFAULTS.
    if (IST%filling_frazil) then
      call get_param(param_file, mod, "FILLING_FRAZIL_TIMESCALE", IST%fraz_fill_time, &
                 "A timescale with which the filling frazil causes the \n"//&
                 "thinest cells to attain similar thicknesses, or a negative \n"//&
                 "number to apply the frazil flux uniformly.", default=0.0, units="s")
    endif
  endif
  call get_param(param_file, mod, "RESTARTFILE", restart_file, &
                 "The name of the restart file.", default="ice_model.res.nc")

  call get_param(param_file, mod, "TIMEUNIT", Time_unit, &
                 "The time unit for ICE_STATS_INTERVAL.", &
                 units="s", default=86400.0)
  call get_param(param_file, mod, "ICE_STATS_INTERVAL",IST%ice_stats_interval, &
                 "The interval in units of TIMEUNIT between writes of the \n"//&
                 "globally summed ice statistics and conservation checks.", &
                 default=set_time(0,1), timeunit=Time_unit)

  if ((IST%ice_bulk_salin > 0.0) .and. (IST%ice_rel_salin > 0.0)) &
    call SIS_error(FATAL, "It is inconsistent to have both ICE_BULK_SALINITY "//&
                   "and ICE_RELATIVE_SALINITY set to positive values.")
  if (IST%ice_bulk_salin < 0.0) IST%ice_bulk_salin = 0.0

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

  call check_ice_model_nml(param_file)

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

  ! Allocate and register fields for restarts.
  call ice_data_type_register_restarts(G%domain%mpp_domain, IG%CatIce, &
                         param_file, Ice, Ice%Ice_restart, restart_file)

  call ice_state_register_restarts(G%domain%mpp_domain, HI, IG, param_file, &
                                   IST, Ice%Ice_restart, restart_file)

  call SIS_slow_register_restarts(G%domain%mpp_domain, HI, IG, param_file, &
                                  IST, Ice%Ice_restart, restart_file)


  ! Register tracers that will be advected around.
  call register_SIS_tracer_pair(IST%enth_ice, IG%NkIce, "enth_ice", &
                                IST%enth_snow, 1, "enth_snow", &
                                G, IG, param_file, IST%TrReg, &
                                massless_iceval=massless_ice_enth*enth_unit, &
                                massless_snowval=massless_snow_enth*enth_unit)

  if (IST%ice_rel_salin > 0.0) then
    call register_SIS_tracer(IST%sal_ice, G, IG, IG%NkIce, "salin_ice", param_file, &
                             IST%TrReg, snow_tracer=.false., &
                             massless_val=massless_ice_salin)
  endif

  !   Register any tracers that will be handled via tracer flow control for 
  ! restarts and advection.
  call SIS_call_tracer_register(G, IG, param_file, IST%SIS_tracer_flow_CSp, &
                                IST%diag, IST%TrReg, Ice%Ice_restart, restart_file)

  ! Redefine the computational domain sizes to use the ice model's indexing convention.
  isc = HI%isc ; iec = HI%iec ; jsc = HI%jsc ; jec = HI%jec
  i_off = LBOUND(Ice%t_surf,1) - HI%isc ; j_off = LBOUND(Ice%t_surf,2) - HI%jsc

  ! This will likely be replaced later with information provided along
  ! with the shortwave fluxes.
  IST%coszen(:,:) = cos(3.14*67.0/180.0) ! NP summer solstice.

  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%ocean_pt(i2,j2) = ( G%mask2dT(i,j) > 0.5 )
    Ice%area(i2,j2) = G%areaT(i,j) * G%mask2dT(i,j)
  enddo ; enddo

  Ice%Time           = Time
  IST%Time           = Time
  IST%Time_Init      = Time_Init
  IST%Time_step_fast = Time_step_fast
  IST%Time_step_slow = Time_step_slow

  IST%avg_count = 0
  IST%n_calls = 0

  call SIS_diag_mediator_init(G, IG, param_file, IST%diag, component="SIS", &
                              doc_file_dir = dirs%output_directory)
  call set_SIS_axes_info(G, IG, param_file, IST%diag)

  call SIS2_ice_thm_init(param_file, IST%ice_thm_CSp, IST%ITV, &
                         init_EOS=IST%nudge_sea_ice)

  if (IST%nudge_sea_ice) then
    allocate(IST%cool_nudge(isc:iec,jsc:jec)) ; IST%cool_nudge(:,:) = 0.0
    allocate(IST%melt_nudge(isc:iec,jsc:jec)) ; IST%melt_nudge(:,:) = 0.0
  endif

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
        IST%sal_ice(:,:,:,1) = IST%ice_bulk_salin
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
    IST%sal_ice(:,:,:,:) = IST%ice_bulk_salin

    enth_spec_snow = Enth_from_TS(0.0, 0.0, IST%ITV)
    IST%enth_snow(:,:,:,1) = enth_spec_snow
    do n=1,IG%NkIce
      enth_spec_ice = Enth_from_TS(0.0, S_col(n), IST%ITV)
      IST%enth_ice(:,:,:,n) = enth_spec_ice
    enddo

    IST%do_init = .true. ! Some more initialization needs to be done in ice_model.
  endif ! file_exist(restart_path)
  deallocate(S_col)

  do k=0,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
    Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo

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

  ! Set a few final things to complete the  setup of the grid. 
  Ice%G%g_Earth = g_Earth
  call set_first_direction(G, first_direction)
  call clone_MOM_domain(G%domain, G%domain_aux, symmetric=.false., &
                        domain_name="ice model aux")

  ! Copy the ice model's domain into one with no halos that can be shared
  ! publicly for use by the exchange grid.
  call clone_MOM_domain(G%domain, Ice%domain, halo_size=0, symmetric=.false., &
                        domain_name="ice_nohalo")

  call ice_diagnostics_init(Ice, IST, G, IST%diag, IST%Time)

  call ice_thm_param(IST%slab_ice, k_snow, h_lo_lim)

  call SIS_slow_init(Ice%Time, G, IG, param_file, IST%diag, IST)

  call SIS_sum_output_init(G, param_file, dirs%output_directory, Time_Init, &
                           IST%sum_output_CSp, IST%ntrunc)

  !   Initialize any tracers that will be handled via tracer flow control.
  call SIS_tracer_flow_control_init(Ice%Time, G, IG, param_file, IST%SIS_tracer_flow_CSp, is_restart)

  call close_param_file(param_file)

  iceClock = mpp_clock_id( 'Ice', flags=clock_flag_default, grain=CLOCK_COMPONENT )
  iceClock1 = mpp_clock_id( 'Ice: bot to top', flags=clock_flag_default, grain=CLOCK_ROUTINE )
  iceClock2 = mpp_clock_id( 'Ice: update slow (dn)', flags=clock_flag_default, grain=CLOCK_ROUTINE )
  iceClock3 = mpp_clock_id( 'Ice: update fast', flags=clock_flag_default, grain=CLOCK_ROUTINE )

  ! Initialize icebergs
  if (IST%do_icebergs) then
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
  endif

  if (IST%add_diurnal_sw .or. IST%do_sun_angle_for_alb) then
    call set_domain(G%Domain%mpp_domain)
    call astronomy_init
    call nullify_domain()
  endif

  ! Do any error checking here.
  if (IST%debug) then
    call ice_grid_chksum(G, haloshift=2)
  endif

  call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IG, IST%sum_output_CSp)
  IST%write_ice_stats_time = Time_Init + IST%ice_stats_interval * &
      (1 + (IST%Time - Time_init) / IST%ice_stats_interval)

  call callTree_leave("ice_model_init()")

end subroutine ice_model_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_end - writes the restart file and deallocates memory               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_end (Ice)
  type(ice_data_type), intent(inout) :: Ice

  type(ice_state_type), pointer :: IST => NULL()

  IST => Ice%Ice_state

  call ice_model_restart(Ice=Ice)

  !--- release memory ------------------------------------------------

  call SIS_slow_end(IST)

  call SIS2_ice_thm_end(IST%ice_thm_CSp, IST%ITV)

  call SIS_hor_grid_end(Ice%G)
  call ice_grid_end(Ice%IG)
  call dealloc_Ice_arrays(Ice)

  call SIS_tracer_flow_control_end(IST%SIS_tracer_flow_CSp)

  call dealloc_IST_arrays(IST)
  deallocate(Ice%Ice_restart)

  ! End icebergs
  if (IST%do_icebergs) call icebergs_end(Ice%icebergs)
  if (IST%add_diurnal_sw .or. IST%do_sun_angle_for_alb) call astronomy_end

  call SIS_diag_mediator_end(IST%Time, IST%diag)

  deallocate(Ice%Ice_state)

end subroutine ice_model_end

end module ice_model_mod
