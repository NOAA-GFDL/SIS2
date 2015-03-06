!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

use SIS_diag_mediator, only : set_SIS_axes_info, SIS_diag_mediator_init
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_error_checking, only : chksum, Bchksum, hchksum, uchksum, vchksum
use SIS_error_checking, only : check_redundant_B, check_redundant_C
use SIS_get_input, only : Get_SIS_input
use SIS_sum_output, only : write_ice_statistics, SIS_sum_output_init
use SIS_sum_output, only : accumulate_bottom_input, accumulate_input_1, accumulate_input_2

use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : fill_symmetric_edges
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : open_param_file, close_param_file
use MOM_string_functions, only : uppercase

use fms_mod, only : file_exist, clock_flag_default
use fms_io_mod, only : set_domain, nullify_domain, restore_state, query_initialized
use fms_io_mod, only : register_restart_field, restart_file_type
use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_LOOP, CLOCK_ROUTINE

use time_manager_mod, only : time_type, time_type_to_real, get_date, get_time
use time_manager_mod, only : set_date, set_time, operator(+), operator(-)
use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)
use astronomy_mod, only: astronomy_init, astronomy_end
use astronomy_mod, only: universal_time, orbital_time, diurnal_solar, daily_mean_solar
  use coupler_types_mod,only: coupler_3d_bc_type
  use constants_mod,    only: hlv, hlf, T_0degC=>Tfreeze, grav, STEFAN
use ocean_albedo_mod, only: compute_ocean_albedo            ! ice sets ocean surface
use ocean_rough_mod,  only: compute_ocean_roughness         ! properties over water

use ice_type_mod, only : ice_data_type, ice_state_type
use ice_type_mod, only : ice_model_restart, dealloc_ice_arrays, dealloc_IST_arrays
use ice_type_mod, only : ice_data_type_register_restarts, ice_state_register_restarts
use ice_type_mod, only : ice_diagnostics_init, ice_stock_pe, check_ice_model_nml
use ice_type_mod, only : ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
use ice_type_mod, only : ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
use ice_type_mod, only : lnd_ice_bnd_type_chksum, ice_data_type_chksum
use ice_type_mod, only : IST_chksum, Ice_public_type_chksum
use ice_type_mod, only : IST_bounds_check, Ice_public_type_bounds_check
use ice_utils_mod, only : get_avg, post_avg, ice_line, ice_grid_chksum
use ice_grid_mod, only : sea_ice_grid_type, set_ice_grid, ice_grid_end, cell_area
use ice_spec_mod, only : get_sea_surface

use SIS_tracer_registry, only : register_SIS_tracer

use ice_thm_mod,   only: slab_ice_optics, ice_thm_param, ice5lay_temp, ice5lay_resize
  use ice_thm_mod,      only: MU_TS, TFI, CI, e_to_melt, get_thermo_coefs
use SIS2_ice_thm,  only: ice_temp_SIS2, ice_resize_SIS2, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,  only: ice_optics_SIS2, get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,  only: enthalpy_from_TS, enth_from_TS, Temp_from_En_S, Temp_from_Enth_S
use SIS2_ice_thm,  only: T_freeze, calculate_T_freeze, enthalpy_liquid, e_to_melt_TS
use ice_dyn_bgrid, only: ice_B_dynamics, ice_B_dyn_init, ice_B_dyn_register_restarts, ice_B_dyn_end
use ice_dyn_cgrid, only: ice_C_dynamics, ice_C_dyn_init, ice_C_dyn_register_restarts, ice_C_dyn_end
use ice_transport_mod, only : ice_transport, ice_transport_init, ice_transport_end
use ice_transport_mod, only : adjust_ice_categories
use ice_bergs,        only: icebergs_run, icebergs_init, icebergs_end, icebergs_incr_mass

implicit none ; private

#include <SIS2_memory.h>

public :: ice_data_type, ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
public :: ice_model_init, ice_model_end, update_ice_model_fast, ice_stock_pe, cell_area
public :: update_ice_model_slow_up, update_ice_model_slow_dn
public :: ice_model_restart  ! for intermediate restarts
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum, ice_data_type_chksum


integer :: iceClock, iceClock1, iceCLock2, iceCLock3, iceClock4, iceClock5, &
           iceClock6, iceClock7, iceClock8, iceClock9, iceClocka, iceClockb, iceClockc

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
  call update_ice_model_slow(Ice, Ice%Ice_state, Ice%G, Land_boundary%runoff, Land_boundary%calving, &
                             Land_boundary%runoff_hflx, Land_boundary%calving_hflx )
  call mpp_clock_end(iceClock2)
  call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow_dn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! sum_top_quantities - sum fluxes for later use by ice/ocean slow physics.     !
!   Nothing here will be exposed to other modules.                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine sum_top_quantities ( Ice, IST, Atmos_boundary_fluxes, flux_u,  flux_v, flux_t,  flux_q, &
       flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif,&
       flux_lw, lprec, fprec, flux_lh, G)
  type(ice_data_type),              intent(inout) :: Ice
  type(ice_state_type),             intent(inout) :: IST
  type(coupler_3d_bc_type),         intent(inout) :: Atmos_boundary_fluxes
  type(sea_ice_grid_type),          intent(inout) :: G
  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:G%CatIce), intent(in) :: &
    flux_u, flux_v, flux_t, flux_q, flux_lw, lprec, fprec, flux_lh
  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:G%CatIce), intent(in) :: &
    flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif

  real,dimension(G%isc:G%iec,G%jsc:G%jec)                 :: tmp
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, ncat
  integer :: ind, max_num_fields, next_index

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce

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
        allocate(IST%tr_flux_top(SZI_(G), SZJ_(G), 0:G%CatIce, IST%num_tr_fluxes))
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

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
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
    call get_avg(flux_lw(:,:,:) + STEFAN*IST%t_surf(isc:iec,jsc:jec,:)**4, &
                       IST%part_size(isc:iec,jsc:jec,:), tmp(:,:))
    do j=jsc,jec ; do i=isc,iec
      if (G%Lmask2dT(i,j)) IST%lwdn(i,j) = IST%lwdn(i,j) + tmp(i,j)
    enddo ; enddo
  endif

  if (IST%id_swdn > 0) then
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%Lmask2dT(i,j)) then
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
subroutine avg_top_quantities(Ice, IST, G)
  type(ice_data_type),     intent(inout) :: Ice
  type(ice_state_type),    intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G

  real    :: u, v, divid, sign
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: tmp2d
  integer :: i, j, k, m, n, isc, iec, jsc, jec, ncat
  logical :: sent

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce

  !
  ! compute average fluxes
  !
  if (IST%avg_count == 0) call SIS_error(FATAL,'avg_top_quantities: '//&
       'no ocean model fluxes have been averaged')

  ! Rotate the stress from lat/lon to ocean coordinates and possibly change the
  ! sign to positive for downward fluxes of positive momentum.
  sign = 1.0 ; if (IST%atmos_winds) sign = -1.0
  divid = 1.0/real(IST%avg_count)
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    u = IST%flux_u_top(i,j,k) * (sign*divid)
    v = IST%flux_v_top(i,j,k) * (sign*divid)
    IST%flux_u_top(i,j,k) = u*G%cos_rot(i,j)-v*G%sin_rot(i,j) ! rotate stress from lat/lon
    IST%flux_v_top(i,j,k) = v*G%cos_rot(i,j)+u*G%sin_rot(i,j) ! to ocean coordinates
  enddo ; enddo ; enddo
  call pass_vector(IST%flux_u_top, IST%flux_v_top, G%Domain, stagger=AGRID)

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
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
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    IST%lwdn(i,j) = IST%lwdn(i,j)* divid
    IST%swdn(i,j) = IST%swdn(i,j)* divid
  enddo ; enddo
  !
  ! Flux diagnostics
  !
  if (IST%id_sh>0) call post_avg(IST%id_sh, IST%flux_t_top, IST%part_size, &
                                 IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_lh>0) call post_avg(IST%id_lh, IST%flux_lh_top, IST%part_size, &
                                 IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_evap>0) call post_avg(IST%id_evap, IST%flux_q_top, IST%part_size, &
                                 IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw>0) then
    do j=jsc,jec ; do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo ; enddo
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * ( &
            IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k) + &
            IST%flux_sw_nir_dir_top(i,j,k) + IST%flux_sw_nir_dif_top(i,j,k) )
    enddo ; enddo ; enddo
    call post_data(IST%id_sw, tmp2d, IST%diag, mask=G%Lmask2dT)
  endif
  if (IST%id_lw>0) call post_avg(IST%id_lw, IST%flux_lw_top, &
                             IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_snofl>0) call post_avg(IST%id_snofl, IST%fprec_top, &
                             IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_rain>0) call post_avg(IST%id_rain, IST%lprec_top, &
                             IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_lwdn>0) call post_data(IST%id_lwdn, IST%lwdn, IST%diag, mask=G%Lmask2dT)
  if (IST%id_swdn>0) call post_data(IST%id_swdn, IST%swdn, IST%diag, mask=G%Lmask2dT)
  if (IST%id_sw_vis>0) then
    do j=jsc,jec ; do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo ; enddo
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * ( &
            IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k) )
    enddo ; enddo ; enddo
    call post_data(IST%id_sw_vis, tmp2d, IST%diag, mask=G%Lmask2dT)
  endif
  if (IST%id_sw_nir_dir>0) call post_avg(IST%id_sw_nir_dir, IST%flux_sw_nir_dir_top, &
                             IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_nir_dif>0) call post_avg(IST%id_sw_nir_dif, IST%flux_sw_nir_dif_top, &
                             IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_vis_dir>0) call post_avg(IST%id_sw_vis_dir, IST%flux_sw_vis_dir_top, &
                             IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_vis_dif>0) call post_avg(IST%id_sw_vis_dif, IST%flux_sw_vis_dif_top, &
                             IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  !
  ! set count to zero and fluxes will be zeroed before the next sum
  !
  IST%avg_count = 0

end subroutine avg_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ocean_top_fluxes - Translate ice-bottom fluxes of heat, mass, salt, and  !
!   tracers from the ice model's internal state to the public ice data type    !
!   for use by the ocean model.                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ocean_top_fluxes(Ice, IST, G)
  type(ice_data_type),  intent(inout) :: Ice
  type(ice_state_type), intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G

  real :: I_count
  integer :: i, j, isc, iec, jsc, jec, m, n, i2, j2, i_off, j_off, ind
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc

  if (IST%debug) then
    call IST_chksum("Start set_ocean_top_fluxes", IST, G)
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

  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off! Use these to correct for indexing differences.
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
  enddo ; enddo
  do n=1,Ice%ocean_fluxes%num_bcs ; do m=1,Ice%ocean_fluxes%bc(n)%num_fields
    ind = IST%tr_flux_index(m,n)
    if (ind < 1) call SIS_error(FATAL, "Bad boundary flux index in set_ocean_top_fluxes.")
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off  ! Use these to correct for indexing differences.
        Ice%ocean_fluxes%bc(n)%field(m)%values(i2,j2) = IST%tr_flux_ocn_top(i,j,ind)
    enddo ; enddo
  enddo ; enddo

  if (IST%debug) then
    call IST_chksum("End set_ocean_top_fluxes", IST, G)
    call Ice_public_type_chksum("End set_ocean_top_fluxes", Ice)
  endif

end subroutine set_ocean_top_fluxes

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! finish_ocean_top_stresses - Finish setting the ice-ocean stresses by dividing!
!   them through the stresses by the number of times they have been augmented. !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine finish_ocean_top_stresses(Ice, IST, G)
  type(ice_data_type),  intent(inout) :: Ice
  type(ice_state_type), intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G

  real :: I_count

  if (IST%stress_count > 1) then
    I_count = 1.0 / IST%stress_count
    Ice%flux_u(:,:) = Ice%flux_u(:,:) * I_count
    Ice%flux_v(:,:) = Ice%flux_v(:,:) * I_count
  endif

  if (IST%debug) then
    call Ice_public_type_chksum("finish_ocean_top_stresses", Ice)
  endif

end subroutine finish_ocean_top_stresses

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ocean_top_stress_Bgrid - Calculate the stresses on the ocean integrated  !
!   across all the thickness categories with the appropriate staggering, and   !
!   store them in the public ice data type for use by the ocean model.  This   !
!   version of the routine uses wind and ice-ocean stresses on a B-grid.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ocean_top_stress_Bgrid(Ice, IST, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, part_size, G)
  type(ice_data_type),  intent(inout) :: Ice
  type(ice_state_type), intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: windstr_x_water, str_ice_oce_x
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: windstr_y_water, str_ice_oce_y
  real, dimension (SZI_(G),SZJ_(G),0:G%CatIce), intent(in) :: part_size

  real    :: ps_vel ! part_size interpolated to a velocity point, nondim.
  integer :: i, j, k, isc, iec, jsc, jec, ncat, i2, j2, k2, i_off, j_off
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc

  if (IST%debug) then
    call IST_chksum("Start set_ocean_top_stress_Bgrid", IST, G)
    call Ice_public_type_chksum("Start set_ocean_top_stress_Bgrid", Ice)
  endif

  if (IST%stress_count == 0) then
    Ice%flux_u(:,:) = 0.0 ; Ice%flux_v(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Bgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.
  if (Ice%flux_uv_stagger == AGRID) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      ps_vel = G%mask2dT(i,j) * part_size(i,j,0)
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + ps_vel * 0.25 * &
          ((windstr_x_water(I,J) + windstr_x_water(I-1,J-1)) + &
           (windstr_x_water(I-1,J) + windstr_x_water(I,J-1)))
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + ps_vel * 0.25 * &
          ((windstr_y_water(I,J) + windstr_y_water(I-1,J-1)) + &
           (windstr_y_water(I-1,J) + windstr_y_water(I,J-1)))
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%Lmask2dT(i,j)) then
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + part_size(i,j,k) * 0.25 * &
          ((str_ice_oce_x(I,J) + str_ice_oce_x(I-1,J-1)) + &
           (str_ice_oce_x(I-1,J) + str_ice_oce_x(I,J-1)))
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + part_size(i,j,k) * 0.25 * &
          ((str_ice_oce_y(I,J) + str_ice_oce_y(I-1,J-1)) + &
           (str_ice_oce_y(I-1,J) + str_ice_oce_y(I,J-1)))
    endif ; enddo ; enddo ; enddo
  elseif (Ice%flux_uv_stagger == BGRID_NE) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      ps_vel = 1.0 ; if (G%Lmask2dBu(I,J)) ps_vel = &
                         0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                               (part_size(i+1,j,0) + part_size(i,j+1,0)) )
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + windstr_x_water(I,J) * ps_vel
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + windstr_y_water(I,J) * ps_vel
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%Lmask2dBu(I,J)) then
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                       (part_size(i+1,j,k) + part_size(i,j+1,k)) )
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + str_ice_oce_x(I,J) * ps_vel
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + str_ice_oce_y(I,J) * ps_vel
    endif ; enddo ; enddo ; enddo
  elseif (Ice%flux_uv_stagger == CGRID_NE) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      ps_vel = 1.0 ; if (G%Lmask2dCu(I,j)) ps_vel = &
                         0.5*(part_size(i+1,j,0) + part_size(i,j,0))
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + ps_vel * &
              0.5 * (windstr_x_water(I,J) + windstr_x_water(I,J-1))
      ps_vel = 1.0 ; if (G%Lmask2dCv(i,J)) ps_vel = &
                         0.5*(part_size(i,j+1,0) + part_size(i,j,0))
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + ps_vel * &
              0.5 * (windstr_y_water(I,J) + windstr_y_water(I-1,J))
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      if (G%Lmask2dCu(I,j)) then
        ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
        Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + ps_vel * &
            0.5 * (str_ice_oce_x(I,J) + str_ice_oce_x(I,J-1))
      endif
      if (G%Lmask2dCv(i,J)) then
        ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
        Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + ps_vel * &
                0.5 * (str_ice_oce_y(I,J) + str_ice_oce_y(I-1,J))
      endif
    enddo ; enddo ; enddo
  else
    call SIS_error(FATAL, "set_ocean_top_stress_Bgrid: Unrecognized flux_uv_stagger.")
  endif

  IST%stress_count = IST%stress_count + 1

  if (IST%debug) then
    call IST_chksum("End set_ocean_top_stress_Bgrid", IST, G)
    call Ice_public_type_chksum("End set_ocean_top_stress_Bgrid", Ice)
  endif

end subroutine set_ocean_top_stress_Bgrid

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ocean_top_stress_Cgrid - Calculate the stresses on the ocean integrated  !
!   across all the thickness categories with the appropriate staggering, and   !
!   store them in the public ice data type for use by the ocean model.  This   !
!   version of the routine uses wind and ice-ocean stresses on a C-grid.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ocean_top_stress_Cgrid(Ice, IST, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, part_size, G)
  type(ice_data_type),  intent(inout) :: Ice
  type(ice_state_type), intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: windstr_x_water, str_ice_oce_x
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: windstr_y_water, str_ice_oce_y
  real, dimension (SZI_(G),SZJ_(G),0:G%CatIce), intent(in) :: part_size

  real    :: ps_vel ! part_size interpolated to a velocity point, nondim.
  integer :: i, j, k, isc, iec, jsc, jec, ncat, i2, j2, k2, i_off, j_off
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc

  if (IST%debug) then
    call IST_chksum("Start set_ocean_top_stress_Cgrid", IST, G)
    call Ice_public_type_chksum("Start set_ocean_top_stress_Cgrid", Ice)
  endif

  if (IST%stress_count == 0) then
    Ice%flux_u(:,:) = 0.0 ; Ice%flux_v(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Cgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.
  if (Ice%flux_uv_stagger == AGRID) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      ps_vel = G%mask2dT(i,j) * part_size(i,j,0)
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + ps_vel * 0.5 * &
                          (windstr_x_water(I,j) + windstr_x_water(I-1,j))
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + ps_vel * 0.5 * &
                          (windstr_y_water(I,j) + windstr_y_water(i,J-1))
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%Lmask2dT(i,j)) then
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) +  part_size(i,j,k) * 0.5 * &
                          (str_ice_oce_x(I,j) + str_ice_oce_x(I-1,j))
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + part_size(i,j,k) * 0.5 * &
                          (str_ice_oce_y(I,j) + str_ice_oce_y(i,J-1))
    endif ; enddo ; enddo ; enddo
  elseif (Ice%flux_uv_stagger == BGRID_NE) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      ps_vel = 1.0 ; if (G%Lmask2dBu(I,J)) ps_vel = &
                         0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                               (part_size(i+1,j,0) + part_size(i,j+1,0)) )
      ! Consider deleting the masks here?
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
              (windstr_x_water(I,j) + windstr_x_water(I,j+1))
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
              (windstr_y_water(I,j) + windstr_y_water(i+1,J))
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (G%Lmask2dBu(I,J)) then
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                       (part_size(i+1,j,k) + part_size(i,j+1,k)) )
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + ps_vel * 0.5 * &
                          (str_ice_oce_x(I,j) + str_ice_oce_x(I,j+1))
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + ps_vel * 0.5 * &
                          (str_ice_oce_y(I,j) + str_ice_oce_y(i+1,J))
    endif ; enddo ; enddo ; enddo
  elseif (Ice%flux_uv_stagger == CGRID_NE) then
    do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      ps_vel = 1.0 ; if (G%Lmask2dCu(I,j)) ps_vel = &
                         0.5*(part_size(i+1,j,0) + part_size(i,j,0))
      Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + ps_vel * windstr_x_water(I,j)
      ps_vel = 1.0 ; if (G%Lmask2dCv(i,J)) ps_vel = &
                         0.5*(part_size(i,j+1,0) + part_size(i,j,0))
      Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + ps_vel * windstr_y_water(i,J)
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ! Use these to correct for indexing differences.
      if (G%Lmask2dCu(I,j)) then
        ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
        Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + ps_vel * str_ice_oce_x(I,j)
      endif
      if (G%Lmask2dCv(i,J)) then
        ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
        Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + ps_vel * str_ice_oce_y(I,j)
      endif
    enddo ; enddo ; enddo
  else
    call SIS_error(FATAL, "set_ocean_top_stress_Cgrid: Unrecognized flux_uv_stagger.")
  endif

  IST%stress_count = IST%stress_count + 1

  if (IST%debug) then
    call IST_chksum("End set_ocean_top_stress_Cgrid", IST, G)
    call Ice_public_type_chksum("End set_ocean_top_stress_Cgrid", Ice)
  endif

end subroutine set_ocean_top_stress_Cgrid

!
! Coupler interface to provide ocean surface data to atmosphere.
!
subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
  type(ocean_ice_boundary_type), intent(inout) :: Ocean_boundary
  type(ice_data_type),           intent(inout) :: Ice

  call mpp_clock_begin(iceClock)
  call mpp_clock_begin(iceClock1)
  call set_ice_surface_state(Ice, Ice%Ice_state, Ocean_boundary%t, Ocean_boundary%u, Ocean_boundary%v, &
                             Ocean_boundary%frazil, Ocean_boundary, Ice%G, &
                             Ocean_boundary%s, Ocean_boundary%sea_level )
  call mpp_clock_end(iceClock1)
  call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow_up

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ice_surface_state - prepare surface state for atmosphere fast physics    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ice_surface_state(Ice, IST, t_surf_ice_bot, u_surf_ice_bot, v_surf_ice_bot, &
                                 frazil_ice_bot, OIB, G, s_surf_ice_bot, sea_lev_ice_bot )
  type(ice_data_type),                     intent(inout) :: Ice
  type(ice_state_type),                    intent(inout) :: IST
  type(sea_ice_grid_type),                 intent(inout) :: G
  real, dimension(G%isc:G%iec,G%jsc:G%jec),   intent(in) :: t_surf_ice_bot, u_surf_ice_bot
  real, dimension(G%isc:G%iec,G%jsc:G%jec),   intent(in) :: v_surf_ice_bot, frazil_ice_bot
  type(ocean_ice_boundary_type),           intent(inout) :: OIB
  real, dimension(G%isc:G%iec,G%jsc:G%jec),   intent(in) :: s_surf_ice_bot, sea_lev_ice_bot

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: m_ice_tot
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: h_ice_input
  real, dimension(SZI_(G),SZJ_(G)) :: u_nonsym, v_nonsym
  real, dimension(G%NkIce) :: sw_abs_lay
  real :: u, v
  real :: area_pt
  real :: I_Nk
  real :: kg_H_Nk  ! The conversion factor from units of H to kg/m2 over Nk.
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  logical :: sent
  real :: H_to_m_ice     ! The specific volumes of ice and snow times the
  real :: H_to_m_snow    ! conversion factor from thickness units, in m H-1.
  real, parameter                  :: LI = hlf
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc
  I_Nk = 1.0 / G%NkIce ; kg_H_Nk = G%H_to_kg_m2 * I_Nk

  H_to_m_snow = G%H_to_kg_m2 / IST%Rho_snow ; H_to_m_ice = G%H_to_kg_m2 / IST%Rho_ice

  ! pass ocean state through ice on first partition
  if (.not. IST%specified_ice) then ! otherwise, already set by update_ice_model_slow
    IST%t_surf(isc:iec,jsc:jec,0) = t_surf_ice_bot(isc:iec,jsc:jec)
  endif

  if (IST%do_init) then
    call get_sea_surface(IST%Time, IST%t_surf(isc:iec,jsc:jec,0), IST%part_size(isc:iec,jsc:jec,0:1), &
                         h_ice_input )
    do j=jsc,jec ; do i=isc,iec
      IST%mH_ice(i,j,1) = h_ice_input(i,j)*(IST%Rho_ice*G%kg_m2_to_H)
    enddo ; enddo

    !   Transfer ice to the correct thickness category.  If do_ridging=.false.,
    ! the first call to ice_redistribute has the same result.  At present, all
    ! tracers are initialized to their default values, and snow is set to 0,
    ! and so do not need to be updated here.
    if (IST%do_ridging) then
      do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,1) > G%mH_cat_bound(1)) then
        do k=ncat,2,-1 ; if (IST%mH_ice(i,j,1) > G%mH_cat_bound(k-1)) then
          IST%part_size(i,j,k) = IST%part_size(i,j,1)
          IST%part_size(i,j,1) = 0.0
          IST%mH_ice(i,j,k) = IST%mH_ice(i,j,1) ; IST%mH_ice(i,j,1) = 0.0
          !	IST%mH_snow(i,j,k) = IST%mH_snow(i,j,1) ; IST%mH_snow(i,j,1) = 0.0
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
    call IST_bounds_check(IST, G, "Start of set_ice_surface_state")

  if (IST%debug) then
    call IST_chksum("Start set_ice_surface_state", IST, G)
    call Ice_public_type_chksum("Start set_ice_surface_state", Ice)
    call chksum(u_surf_ice_bot(isc:iec,jsc:jec), "Start IB2IT u_surf_ice_bot")
    call chksum(v_surf_ice_bot(isc:iec,jsc:jec), "Start IB2IT v_surf_ice_bot")
  endif

  do j=jsc,jec ; do i=isc,iec
    IST%t_ocn(i,j) = t_surf_ice_bot(i,j) - T_0degC

    IST%s_surf(i,j) = s_surf_ice_bot(i,j)
    IST%frazil(i,j) = frazil_ice_bot(i,j)
    IST%sea_lev(i,j) = sea_lev_ice_bot(i,j)
  enddo ; enddo


! Transfer the ocean state for extra tracer fluxes.
  do n=1,OIB%fields%num_bcs  ; do m=1,OIB%fields%bc(n)%num_fields
    Ice%ocean_fields%bc(n)%field(m)%values(:,:,1) = OIB%fields%bc(n)%field(m)%values
  enddo ; enddo

  m_ice_tot(:,:) = 0.0
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    IST%tmelt(i,j,k) = 0.0 ; IST%bmelt(i,j,k) = 0.0
    m_ice_tot(i,j) = m_ice_tot(i,j) + IST%mH_ice(i,j,k) * IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    if (m_ice_tot(i,j) > 0.0) then
      IST%bheat(i,j) = IST%kmelt*(IST%t_ocn(i,j) - T_Freeze(IST%s_surf(i,j), IST%ITV))
    else
      IST%bheat(i,j) = 0.0
    endif
  enddo ; enddo

  Ice%ice_mask(:,:,1) = .false.
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%ice_mask(i2,j2,k2) = (IST%mH_ice(i,j,k) > 0.0)
  enddo ; enddo ; enddo

  if (IST%slab_ice) then
    IST%sw_abs_sfc(:,:,:) = 0.0 ; IST%sw_abs_snow(:,:,:) = 0.0
    IST%sw_abs_ice(:,:,:,:) = 0.0 ; IST%sw_abs_ocn(:,:,:) = 0.0
    IST%sw_abs_int(:,:,:) = 0.0
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
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
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      call ice_optics_SIS2(IST%mH_snow(i,j,k)*H_to_m_snow, IST%mH_ice(i,j,k)*H_to_m_ice, &
               IST%t_surf(i,j,k)-T_0degC, T_Freeze(IST%s_surf(i,j),IST%ITV), G%NkIce, &
               Ice%albedo_vis_dir(i2,j2,k2), Ice%albedo_vis_dif(i2,j2,k2), &
               Ice%albedo_nir_dir(i2,j2,k2), Ice%albedo_nir_dif(i2,j2,k2), &
               IST%sw_abs_sfc(i,j,k),  IST%sw_abs_snow(i,j,k), &
               sw_abs_lay, IST%sw_abs_ocn(i,j,k), IST%sw_abs_int(i,j,k), &
               IST%ice_thm_CSp, coszen_in=IST%coszen(i,j))

      do m=1,G%NkIce ; IST%sw_abs_ice(i,j,k,m) = sw_abs_lay(m) ; enddo

      !Niki: Is the following correct for diagnostics?
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
      if (G%symmetric) then  ! This is a place-holder until the Tikal release.
        u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
        do j=jsc,jec ; do i=isc,iec
          u_nonsym(i,j) = u_surf_ice_bot(i,j) ; v_nonsym(i,j) = v_surf_ice_bot(i,j)
        enddo ; enddo
        call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, stagger=BGRID_NE)

        ! The under-ice current is needed for the water drag term.
        do J=jsc-1,jec ; do I=isc-1,iec
          IST%u_ocn(I,J) = u_nonsym(I,J) ; IST%v_ocn(I,J) = v_nonsym(I,J)
        enddo ; enddo
      else
        do J=jsc,jec ; do I=isc,iec
          IST%u_ocn(I,J) = u_surf_ice_bot(I,J) ! need under-ice current
          IST%v_ocn(I,J) = v_surf_ice_bot(I,J) ! for water drag term
        enddo ; enddo
      endif
   !   This will be used with Tikal and later shared code.  However, it does
   ! not appear to work properly yet.
   !   do J=jsc,jec ; do I=isc,iec
   !     IST%u_ocn(I,J) = u_surf_ice_bot(I,J) ! need under-ice current
   !     IST%v_ocn(I,J) = v_surf_ice_bot(I,J) ! for water drag term
   !   enddo ; enddo
   !   if (G%symmetric) &
   !     call fill_symmetric_edges(IST%u_ocn, IST%v_ocn, G%Domain, stagger=BGRID_NE)

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
      ! This can only be used with Tikal and later shared code.
      if (G%symmetric) &
        call fill_symmetric_edges(IST%u_ocn_C, IST%v_ocn_C, G%Domain, stagger=CGRID_NE)

      call pass_vector(IST%u_ocn_C, IST%v_ocn_C, G%Domain, stagger=CGRID_NE)
    else
      u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
      do j=jsc,jec ; do i=isc,iec
        u_nonsym(i,j) = u_surf_ice_bot(i,j) ; v_nonsym(i,j) = v_surf_ice_bot(i,j)
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
    call IST_bounds_check(IST, G, "Midpoint set_ice_surface_state")

  ! Copy the surface temperatures into the externally visible data type.
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%s_surf(i2,j2) = IST%s_surf(i,j)
  enddo ; enddo
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
    Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  ! put ocean and ice velocities into Ice%u_surf/v_surf on t-cells
  if (IST%Cgrid_dyn) then
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
      IST%enth_prev(i,j,k) = (IST%mH_snow(i,j,k)*G%H_to_kg_m2) * IST%enth_snow(i,j,k,1)
      do m=1,G%NkIce
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
                     IST%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (IST%id_sst>0) call post_data(IST%id_sst, IST%t_ocn, IST%diag, mask=G%Lmask2dT)
  if (IST%id_sss>0) call post_data(IST%id_sss, IST%s_surf, IST%diag, mask=G%Lmask2dT)
  if (IST%id_ssh>0) call post_data(IST%id_ssh, IST%sea_lev, IST%diag, mask=G%Lmask2dT)
  if (IST%Cgrid_dyn) then
    if (IST%id_uo>0) call post_data(IST%id_uo, IST%u_ocn_C, IST%diag, mask=G%Lmask2dCu)
    if (IST%id_vo>0) call post_data(IST%id_vo, IST%v_ocn_C, IST%diag, mask=G%Lmask2dCv)
    if (IST%id_uo_filt>0) call post_data(IST%id_uo_filt, IST%u_ocn_filt, IST%diag, mask=G%Lmask2dCu)
    if (IST%id_vo_filt>0) call post_data(IST%id_vo_filt, IST%v_ocn_filt, IST%diag, mask=G%Lmask2dCv)
  else
    if (IST%id_uo>0) call post_data(IST%id_uo, IST%u_ocn, IST%diag, mask=G%Lmask2dBu)
    if (IST%id_vo>0) call post_data(IST%id_vo, IST%v_ocn, IST%diag, mask=G%Lmask2dBu)
    if (IST%id_uo_filt>0) call post_data(IST%id_uo_filt, IST%u_ocn_filt, IST%diag, mask=G%Lmask2dBu)
    if (IST%id_vo_filt>0) call post_data(IST%id_vo_filt, IST%v_ocn_filt, IST%diag, mask=G%Lmask2dBu)
  endif
  if (IST%id_bheat>0) call post_data(IST%id_bheat, IST%bheat, IST%diag, mask=G%Lmask2dT)
  call disable_SIS_averaging(IST%diag)

  if (IST%debug) then
    call IST_chksum("End set_ice_surface_state", IST, G)
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

  call do_update_ice_model_fast( Atmos_boundary, Ice, Ice%Ice_state, Ice%G )

  call mpp_clock_end(iceClock3)
  call mpp_clock_end(iceClock)

end subroutine update_ice_model_fast

subroutine do_update_ice_model_fast( Atmos_boundary, Ice, IST, G )

  type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
  type(ice_data_type),           intent(inout) :: Ice
  type(ice_state_type),          intent(inout) :: IST
  type(sea_ice_grid_type),       intent(inout) :: G

  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:G%CatIce) :: &
    flux_t, flux_q, flux_lh, flux_lw, &
    flux_sw_nir_dir, flux_sw_nir_dif, &
    flux_sw_vis_dir, flux_sw_vis_dif, &
    flux_u, flux_v, lprec, fprec, &
    dhdt, &   ! The derivative of the upward sensible heat flux with the surface
              ! temperature in W m-2 K-1.
    dedt, &   ! The derivative of the sublimation rate with the surface
              ! temperature, in kg m-2 s-1 K-1 (I think).
    drdt      ! The derivative of the upward radiative heat flux with surface
              ! temperature (i.e. d(flux)/d(surf_temp) in W m-2 K-1.
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: &
    diurnal_factor, cosz_alb, tmp_diag
  real, dimension(0:G%NkIce) :: T_col ! The temperature of a column of ice and snow in degC.
  real, dimension(G%NkIce)   :: S_col ! The thermodynamic salinity of a column of ice, in g/kg.
  real, dimension(0:G%NkIce) :: enth_col   ! The enthalpy of a column of snow and ice, in enth_unit (J/kg?).
  real, dimension(0:G%NkIce) :: SW_abs_col
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

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Atmos_boundary%t_flux,1) - G%isc
  j_off = LBOUND(Atmos_boundary%t_flux,2) - G%jsc
  NkIce = G%NkIce ; I_Nk = 1.0 / G%NkIce ; kg_H_Nk = G%H_to_kg_m2 * I_Nk

  rad = acos(-1.)/180.

  IST%n_fast = IST%n_fast + 1

  if (IST%debug) then
    call IST_chksum("Start do_update_ice_model_fast", IST, G)
    call Ice_public_type_chksum("Start do_update_ice_model_fast", Ice)
  endif

  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    IST%coszen(i,j) = Atmos_boundary%coszen(i2,j2,1)
    Ice%p_surf(i2,j2) = Atmos_boundary%p(i2,j2,1)
  enddo ; enddo

  !   Set up local copies of fluxes.  The Atmos_boundary arrays may have
  ! different index conventions than are used internally in this component.
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
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
  enddo ; enddo ; enddo

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
    do j=jsc,jec ; do i=isc,iec
      call diurnal_solar(G%geoLatT(i,j)*rad, G%geoLonT(i,j)*rad, IST%Time, cosz=cosz_dt_ice, &
                         fracday=fracday_dt_ice, rrsun=rrsun_dt_ice, dt_time=Dt_ice)
      call daily_mean_solar (G%geoLatT(i,j)*rad, time_since_ae, cosz_day, fracday_day, rrsun_day)
      diurnal_factor(i,j) = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice /   &
                                   max(1e-30, cosz_day*fracday_day*rrsun_day)
    enddo ; enddo

    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      flux_sw_nir_dir(i,j,k) = flux_sw_nir_dir(i,j,k) * diurnal_factor(i,j)
      flux_sw_nir_dif(i,j,k) = flux_sw_nir_dif(i,j,k) * diurnal_factor(i,j)
      flux_sw_vis_dir(i,j,k) = flux_sw_vis_dir(i,j,k) * diurnal_factor(i,j)
      flux_sw_vis_dif(i,j,k) = flux_sw_vis_dif(i,j,k) * diurnal_factor(i,j)
    enddo ; enddo ; enddo
  endif

  do j=jsc,jec ; do i=isc,iec
    flux_lh(i,j,0) = hlv * flux_q(i,j,0)
  enddo ; enddo

  !
  ! implicit update of ice surface temperature
  !
  dt_fast = time_type_to_real(IST%Time_step_fast)

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units)

  if (IST%SIS1_5L_thermo) then
    if (G%NkIce /= 4) call SIS_error(FATAL, "SIS1_5L_thermodynamics requires that NK_ICE=4.")

    H_to_m_snow = G%H_to_kg_m2 / IST%Rho_snow
    H_to_m_ice  = G%H_to_kg_m2 / IST%Rho_ice
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j), IST%ITV)
      if (IST%mH_ice(i,j,k) > 0.0) then
        T_col(0) = Temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
        do m=1,NkIce ; T_col(m) = Temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV) ; enddo

        if (IST%slab_ice) then
          latent         = hlv+hlf
        elseif (IST%mH_snow(i,j,k)>0.0) then
          latent         = hlv + (hlf-CI*T_col(0))
        else
          latent         = hlv + hlf*(1-TFI/T_col(1))
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
        flux_lh(i,j,k) = hlv * flux_q(i,j,k)
      endif
    enddo ; enddo ; enddo
  else
    enth_liq_0 = Enth_from_TS(0.0, 0.0, IST%ITV) ; I_enth_unit = 1.0 / enth_units

    T_freeze_ice_top = T_Freeze(S_col(1), IST%ITV)
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j), IST%ITV)
      if (IST%mH_ice(i,j,k) > 0.0) then
        enth_col(0) = IST%enth_snow(i,j,k,1)
        do m=1,NkIce ; enth_col(m) = IST%enth_ice(i,j,k,m) ; enddo
        
        ! In the case of sublimation of either snow or ice, the vapor is at 0 C. 
        ! If the vapor should be at a different temperature, a correction would be
        ! made here.
        if (IST%slab_ice) then
          latent = hlv + hlf
        elseif (IST%mH_snow(i,j,k)>0.0) then
          latent = hlv + (enth_liq_0 - IST%enth_snow(i,j,k,1)) * I_enth_unit
        else
          latent = hlv + (enth_liq_0 - IST%enth_ice(i,j,k,1)) * I_enth_unit
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
        call ice_temp_SIS2(IST%mH_snow(i,j,k)*G%H_to_kg_m2, IST%mH_ice(i,j,k)*G%H_to_kg_m2, &
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

          enth_here = (G%H_to_kg_m2*IST%mH_snow(i,j,k)) * enth_col(0)
          do m=1,G%NkIce
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
        flux_lh(i,j,k) = hlv * flux_q(i,j,k)
      endif
    enddo ; enddo ; enddo

  endif

  ! This routine works on the boundary exchange state.
  call compute_ocean_roughness (Ice%mask, Atmos_boundary%u_star(:,:,1), Ice%rough_mom(:,:,1), &
                                Ice%rough_heat(:,:,1), Ice%rough_moist(:,:,1)  )

  ! This routine works on the boundary exchange state.
  if (IST%do_sun_angle_for_alb) then
    call diurnal_solar(G%geoLatT(isc:iec,jsc:jec)*rad, G%geoLonT(isc:iec,jsc:jec)*rad, &
                 IST%time, cosz=cosz_alb, fracday=diurnal_factor, rrsun=rrsun_dt_ice, dt_time=Dt_ice)  !diurnal_factor as dummy
    call compute_ocean_albedo(Ice%mask, cosz_alb(:,:), Ice%albedo_vis_dir(:,:,1),&
                              Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                              Ice%albedo_nir_dif(:,:,1), rad*G%geoLatT(isc:iec,jsc:jec) )
  else
    call compute_ocean_albedo(Ice%mask, IST%coszen(isc:iec,jsc:jec), Ice%albedo_vis_dir(:,:,1),&
                              Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                              Ice%albedo_nir_dif(:,:,1), rad*G%geoLatT(isc:iec,jsc:jec) )
  endif

  call sum_top_quantities(Ice, IST, Atmos_boundary%fluxes, flux_u, flux_v, flux_t, &
    flux_q, flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif, &
    flux_lw, lprec, fprec, flux_lh, G )

  IST%Time = IST%Time + IST%Time_step_fast ! advance time
  Ice%Time = IST%Time

  ! Copy the surface temperatures into the externally visible data type.
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%s_surf(i2,j2) = IST%s_surf(i,j)
  enddo ; enddo
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
    !   I do not know whether this is needed.  It should have been set in
    ! set_ice_surface_state and not changed since.
    Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  call enable_SIS_averaging(dt_fast, IST%Time, IST%diag)
  if (IST%id_alb_vis_dir>0) call post_avg(IST%id_alb_vis_dir, Ice%albedo_vis_dir, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (IST%id_alb_vis_dif>0) call post_avg(IST%id_alb_vis_dif, Ice%albedo_vis_dif, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (IST%id_alb_nir_dir>0) call post_avg(IST%id_alb_nir_dir, Ice%albedo_nir_dir, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (IST%id_alb_nir_dif>0) call post_avg(IST%id_alb_nir_dif, Ice%albedo_nir_dif, &
                             IST%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)

  if (IST%id_sw_abs_sfc>0) call post_avg(IST%id_sw_abs_sfc, IST%sw_abs_sfc, &
                                   IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_abs_snow>0) call post_avg(IST%id_sw_abs_snow, IST%sw_abs_snow, &
                                   IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  do m=1,G%NkIce
    if (IST%id_sw_abs_ice(m)>0) call post_avg(IST%id_sw_abs_ice(m), IST%sw_abs_ice(:,:,:,m), &
                                     IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  enddo
  if (IST%id_sw_abs_ocn>0) call post_avg(IST%id_sw_abs_ocn, IST%sw_abs_ocn, &
                                   IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)

  if (IST%id_sw_pen>0) then
    tmp_diag(:,:) = 0.0
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      tmp_diag(i,j) = tmp_diag(i,j) + IST%part_size(i,j,k) * &
                     (IST%sw_abs_ocn(i,j,k) + IST%sw_abs_int(i,j,k))
    enddo ; enddo ; enddo
    call post_data(IST%id_sw_pen, tmp_diag, IST%diag, mask=G%Lmask2dT)
  endif

  if (IST%id_coszen>0) call post_data(IST%id_coszen, IST%coszen, IST%diag, mask=G%Lmask2dT)
  call disable_SIS_averaging(IST%diag)

  if (IST%debug) then
    call IST_chksum("End do_update_ice_model_fast", IST, G)
    call Ice_public_type_chksum("End do_update_ice_model_fast", Ice)
  endif

  if (IST%bounds_check) &
    call Ice_public_type_bounds_check(Ice, G, "End update_ice_fast")
  if (IST%bounds_check) &
    call IST_bounds_check(IST, G, "End of update_ice_fast")

end subroutine do_update_ice_model_fast

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! update_ice_model_slow - do ice dynamics, transport, and mass changes         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine update_ice_model_slow(Ice, IST, G, runoff, calving, &
                                 runoff_hflx, calving_hflx)

  type(ice_data_type),                   intent(inout) :: Ice
  type(ice_state_type),                  intent(inout) :: IST
  type(sea_ice_grid_type),               intent(inout) :: G
  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: runoff, calving
  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: runoff_hflx, calving_hflx

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: h2o_chg_xprt, mass, tmp2d
  real, dimension(SZI_(G),SZJ_(G),G%CatIce,G%NkIce) :: &
    temp_ice    ! A diagnostic array with the ice temperature in degC.
  real, dimension(SZI_(G),SZJ_(G),G%CatIce) :: &
    temp_snow   ! A diagnostic array with the snow temperature in degC.
  real, dimension(G%isc:G%iec,G%jsc:G%jec,2) :: Obs_cn_ice      ! partition 2 = ice concentration
  real, dimension(SZI_(G),SZJ_(G))   :: &
    hi_avg, &           ! The area-weighted average ice thickness, in m.
    h_ice_input, &      ! The specified ice thickness, with specified_ice, in m.
    ms_sum, mi_sum, &   ! Masses of snow and ice per unit total area, in kg m-2.
    ice_free, &         ! The fractional open water; nondimensional, between 0 & 1.
    ice_cover, &        ! The fractional ice coverage, summed across all
                        ! thickness categories; nondimensional, between 0 & 1.
    WindStr_x_A, &      ! Zonal (_x_) and meridional (_y_) wind stresses 
    WindStr_y_A, &      ! averaged over the ice categories on an A-grid, in Pa.
    WindStr_x_ocn_A, &  ! Zonal (_x_) and meridional (_y_) wind stresses on the
    WindStr_y_ocn_A, &  ! ice-free ocean on an A-grid, in Pa.
    ice_free_in, &      ! The initial fractional open water; nondimensional, between 0 & 1.
    ice_cover_in, &     ! The initial fractional ice coverage, summed across all
                        ! thickness categories; nondimensional, between 0 & 1.
    WindStr_x_A_in, &   ! Initial onal (_x_) and meridional (_y_) wind stresses 
    WindStr_y_A_in      ! averaged over the ice categories on an A-grid, in Pa.
 real, dimension(SZIB_(G),SZJB_(G)) :: &
    WindStr_x_B, &      ! Zonal (_x_) and meridional (_y_) wind stresses 
    WindStr_y_B, &      ! averaged over the ice categories on a B-grid, in Pa.
    WindStr_x_ocn_B, WindStr_y_ocn_B, & ! Wind stresses on the ice-free ocean on a B-grid, in Pa.
    str_x_ice_ocn_B, str_y_ice_ocn_B  ! Ice-ocean stresses on a B-grid, in Pa.
  real, dimension(SZIB_(G),SZJ_(G))  :: &
    WindStr_x_Cu, &   ! Zonal wind stress averaged over the ice categores on C-grid u-points, in Pa.
    WindStr_x_ocn_Cu, & ! Zonal wind stress on the ice-free ocean on C-grid u-points, in Pa.
    str_x_ice_ocn_Cu   ! Zonal ice-ocean stress on C-grid u-points, in Pa.
  real, dimension(SZI_(G),SZJB_(G))  :: &
    WindStr_y_Cv, &   ! Meridional wind stress averaged over the ice categores on C-grid v-points, in Pa.
    WindStr_y_ocn_Cv, & ! Meridional wind stress on the ice-free ocean on C-grid v-points, in Pa.
    str_y_ice_ocn_Cv  ! Meridional ice-ocean stress on C-grid v-points, in Pa.
  real, dimension(SZIB_(G),SZJ_(G))  :: uc ! Ice velocities interpolated onto
  real, dimension(SZI_(G),SZJB_(G))  :: vc ! a C-grid, in m s-1.

  real, dimension(SZI_(G),SZJ_(G))   :: diagVar ! An temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBx ! An temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBy ! An temporary array for diagnostics.

  real, dimension(SZIB_(G),SZJB_(G)) :: wts  ! A sum of the weights by category.
  real :: weights  ! A sum of the weights around a point.
  real :: I_wts    ! 1.0 / wts or 0 if wts is 0, nondim.
  real :: ps_vel   ! The fractional thickness catetory coverage at a velocity point.

  real :: H_to_m_ice  ! The specific volume of ice times the conversion factor
                      ! from thickness units, in m H-1.
  real :: tot_frazil
  real :: area_h, area_pt
  real :: dt_slow
  real :: dt_slow_dyn
  integer :: ndyn_steps
  real :: Idt_slow
  real :: I_Nk
  integer :: i, j, k, l, m, isc, iec, jsc, jec, ncat, NkIce, nds
  integer :: isd, ied, jsd, jed
  integer :: i2, j2, k2, i_off, j_off
  integer ::iyr, imon, iday, ihr, imin, isec

  real, dimension(G%NkIce) :: S_col ! Specified thermodynamic salinity of each
                                    ! ice layer if spec_thermo_sal is true.
  real :: heat_fill_val   ! A value of enthalpy to use for massless categories.
  logical :: spec_thermo_sal
  logical :: do_temp_diags
  real :: enth_units, I_enth_units, alpha_filt, beta_filt

  real, dimension(SZI_(G),SZJ_(G),G%CatIce) :: &
    rdg_frac, & ! fraction of ridged ice per category
    mi_old      ! Ice mass per unit area before thermodynamics.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rdg_open, & ! formation rate of open water due to ridging
    rdg_vosh, & ! rate of ice volume shifted from level to ridged ice
!   rdg_s2o, &  ! snow volume [m] dumped into ocean during ridging
    rdg_rate, & ! Niki: Where should this come from?
    snow2ocn
  real    :: tmp3

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; NkIce = G%NkIce
  I_Nk = 1.0 / G%NkIce
  i_off = LBOUND(Ice%runoff,1) - G%isc ; j_off = LBOUND(Ice%runoff,2) - G%jsc
  dt_slow = time_type_to_real(IST%Time_step_slow) ; Idt_slow = 1.0/dt_slow

  if (IST%ocean_filter_dt>0.) then
    ! TBD: To initialize the running mean we should set alpha=1,beta=0 for the
    ! first time-step of a new run. -AJA
    ! if (query_initialized(Ice%Ice_restart, 'sea_lev_filt')) then ...
    alpha_filt = min(1., (2. * dt_slow) / ( dt_slow + IST%ocean_filter_dt ) )
    beta_filt = max(0., (IST%ocean_filter_dt - dt_slow) / ( dt_slow + IST%ocean_filter_dt ) )
  endif

  ndyn_steps = 1
  if ((IST%dt_ice_dyn > 0.0) .and. (IST%dt_ice_dyn < dt_slow)) &
    ndyn_steps = max(CEILING(dt_slow/IST%dt_ice_dyn - 0.000001), 1)

  dt_slow_dyn = dt_slow / ndyn_steps

  IST%n_calls = IST%n_calls + 1
  IST%stress_count = 0

  if (IST%debug) then
    call IST_chksum("Start update_ice_model_slow", IST, G)
    call Ice_public_type_chksum("Start update_ice_model_slow", Ice)
  endif

  if (IST%bounds_check) &
    call IST_bounds_check(IST, G, "Start of update_ice_model_slow")

  !
  ! Set up fluxes
  !

  ! save liquid runoff for ocean
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%runoff(i2,j2)  = runoff(i,j)
    Ice%calving(i2,j2) = calving(i,j)
    Ice%runoff_hflx(i2,j2)  = runoff_hflx(i,j)
    Ice%calving_hflx(i2,j2) = calving_hflx(i,j)
  enddo ; enddo

  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

  if (IST%id_runoff>0) &
    call post_data(IST%id_runoff, runoff, IST%diag, mask=Ice%mask )
  if (IST%id_calving>0) &
    call post_data(IST%id_calving, calving, IST%diag, mask=Ice%mask )
  if (IST%id_runoff_hflx>0) &
    call post_data(IST%id_runoff_hflx, runoff_hflx, IST%diag, mask=Ice%mask )
  if (IST%id_calving_hflx>0) &
    call post_data(IST%id_calving_hflx, calving_hflx, IST%diag, mask=Ice%mask )
  if (IST%id_frazil>0) &
    call post_data(IST%id_frazil, IST%frazil*Idt_slow, IST%diag, mask=G%Lmask2dT)

  call avg_top_quantities(Ice, IST, G) ! average fluxes from update_ice_model_fast

  do j=jsc,jec ; do i=isc,iec
    IST%frazil_input(i,j) = IST%frazil(i,j)

    IST%Enth_Mass_in_atm(i,j) = 0.0 ; IST%Enth_Mass_out_atm(i,j) = 0.0
    IST%Enth_Mass_in_ocn(i,j) = 0.0 ; IST%Enth_Mass_out_ocn(i,j) = 0.0
  enddo ; enddo

  !
  ! conservation checks: top fluxes
  !
  call mpp_clock_begin(iceClock7)
  call accumulate_input_1(IST, Ice, dt_slow, G, IST%sum_output_CSp)
  if (IST%column_check) &
    call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IST%sum_output_CSp, &
                              message="    Start of update", check_column=.true.)
  call mpp_clock_end(iceClock7)

  ! Determine the fractional ice coverage and the wind stresses averaged
  ! across all the ice thickness categories on an A-grid.  This is done
  ! over the entire data domain for safety.
  WindStr_x_A(:,:) = 0.0 ; WindStr_y_A(:,:) = 0.0 ; ice_cover(:,:) = 0.0
  do k=1,ncat ; do j=jsd,jed ; do i=isd,ied
    WindStr_x_A(i,j) = WindStr_x_A(i,j) + IST%part_size(i,j,k) * IST%flux_u_top(i,j,k)
    WindStr_y_A(i,j) = WindStr_y_A(i,j) + IST%part_size(i,j,k) * IST%flux_v_top(i,j,k)
    ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
  enddo ; enddo ; enddo
  do j=jsd,jed ; do i=isd,ied
    if (ice_cover(i,j) > 0.0) then
      I_wts = 1.0 / ice_cover(i,j)
      WindStr_x_A(i,j) = WindStr_x_A(i,j) * I_wts
      WindStr_y_A(i,j) = WindStr_y_A(i,j) * I_wts
      if (ice_cover(i,j) > 1.0) ice_cover(i,j) = 1.0

      ! The max with 0 in the following line is here for safety; the only known
      ! instance where it has been required is when reading a SIS-1-derived
      ! restart file with tiny negative concentrations. SIS2 should not need it.
      ice_free(i,j) = max(IST%part_size(i,j,0), 0.0)
  !    Rescale to add up to 1?
  !    I_wts = 1.0 / (ice_free(i,j) + ice_cover(i,j))
  !    ice_free(i,j) = ice_free(i,j) * I_wts ; ice_cover(i,j) = ice_cover(i,j) * I_wts
    else
      ice_free(i,j) = 1.0 ; ice_cover(i,j) = 0.0
    endif
    WindStr_x_ocn_A(i,j) = IST%flux_u_top(i,j,0)
    WindStr_y_ocn_A(i,j) = IST%flux_v_top(i,j,0)

    ice_cover_in(i,j) = ice_cover(i,j) ; ice_free_in(i,j) = ice_free(i,j)
    WindStr_x_A_in(i,j) = WindStr_x_A(i,j) ; WindStr_y_A_in(i,j) = WindStr_y_A(i,j)
  enddo ; enddo

  if (IST%ocean_filter_dt>0.) then
    if (IST%Cgrid_dyn) then
      IST%u_ocn_filt(:,:) = alpha_filt * IST%u_ocn_C(:,:) + beta_filt * IST%u_ocn_filt(:,:)
      IST%v_ocn_filt(:,:) = alpha_filt * IST%v_ocn_C(:,:) + beta_filt * IST%v_ocn_filt(:,:)
    else
      IST%u_ocn_filt(:,:) = alpha_filt * IST%u_ocn(:,:) + beta_filt * IST%u_ocn_filt(:,:)
      IST%v_ocn_filt(:,:) = alpha_filt * IST%v_ocn(:,:) + beta_filt * IST%v_ocn_filt(:,:)
    endif
    IST%sea_lev_filt(:,:) = alpha_filt * IST%sea_lev(:,:) + beta_filt * IST%sea_lev_filt(:,:)
  else
    if (IST%Cgrid_dyn) then
      IST%u_ocn_filt => IST%u_ocn_C
      IST%v_ocn_filt => IST%v_ocn_C
    else
      IST%u_ocn_filt => IST%u_ocn
      IST%v_ocn_filt => IST%v_ocn
    endif
    IST%sea_lev_filt => IST%sea_lev
  endif

  ! Calve off icebergs and integrate forward iceberg trajectories
  if (IST%do_icebergs) then
    call mpp_clock_end(iceClock2) ; call mpp_clock_end(iceClock) ! Stop the sea-ice clocks.
    H_to_m_ice = G%H_to_kg_m2 / IST%Rho_ice
    call get_avg(IST%mH_ice, IST%part_size(:,:,1:), hi_avg, wtd=.true.)
    hi_avg(:,:) = hi_avg(:,:) * H_to_m_Ice
    if (IST%Cgrid_dyn) then
      call icebergs_run( Ice%icebergs, IST%Time, &
              Ice%calving(:,:), IST%u_ocn_filt(isc-2:iec+1,jsc-1:jec+1), &
              IST%v_ocn_filt(isc-1:iec+1,jsc-2:jec+1), IST%u_ice_C(isc-2:iec+1,jsc-1:jec+1), &
              IST%v_ice_C(isc-1:iec+1,jsc-2:jec+1), &
              Ice%flux_u(:,:), Ice%flux_v(:,:), &
              IST%sea_lev_filt(isc-1:iec+1,jsc-1:jec+1), IST%t_surf(isc:iec,jsc:jec,0),  &
              Ice%calving_hflx(:,:), ice_cover, hi_avg, stagger=CGRID_NE, &
              stress_stagger=Ice%flux_uv_stagger)
    else
      call icebergs_run( Ice%icebergs, IST%Time, &
              Ice%calving(:,:), IST%u_ocn_filt(isc-1:iec+1,jsc-1:jec+1), &
              IST%v_ocn_filt(isc-1:iec+1,jsc-1:jec+1), IST%u_ice_B(isc-1:iec+1,jsc-1:jec+1), &
              IST%v_ice_B(isc-1:iec+1,jsc-1:jec+1), &
              Ice%flux_u(:,:), Ice%flux_v(:,:), &
              IST%sea_lev_filt(isc-1:iec+1,jsc-1:jec+1), IST%t_surf(isc:iec,jsc:jec,0),  &
              Ice%calving_hflx(:,:), ice_cover, hi_avg, stagger=BGRID_NE, &
              stress_stagger=Ice%flux_uv_stagger)
    endif
    call mpp_clock_begin(iceClock) ; call mpp_clock_begin(iceClock2) ! Restart the sea-ice clocks.
  endif

  !
  ! Thermodynamics
  !
  if (.not.IST%interspersed_thermo) then
    !TOM> Store old ice mass per unit area for calculating partial ice growth.
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      mi_old(i,j,k) = IST%mH_ice(i,j,k)
    enddo ; enddo ; enddo
    !TOM> derive ridged ice fraction prior to thermodynamic changes of ice thickness
    !     in order to subtract ice melt proportionally from ridged ice volume (see below)
    if (IST%do_ridging) then
      do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
        tmp3 = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
        rdg_frac(i,j,k) = 0.0 ; if (tmp3 > 0.0) &
            rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp3
      enddo ; enddo ; enddo
    endif

    call disable_SIS_averaging(IST%diag)

    ! The thermodynamics routines return updated values of the ice and snow
    ! masses-per-unit area and enthalpies.
    call accumulate_input_2(IST, Ice, IST%part_size, dt_slow, G, IST%sum_output_CSp)
    if (IST%SIS1_5L_thermo) then
      call SIS1_5L_thermodynamics(Ice, IST, G)
    else
      call SIS2_thermodynamics(Ice, IST, G)
    endif

    call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

    !TOM> calculate partial ice growth for ridging and aging.
    if (IST%do_ridging) then
      !     ice growth (Ice%mH_ice > mi_old) does not affect ridged ice volume
      !     ice melt   (ice%mH_ice < mi_old) reduces ridged ice volume proportionally
      do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
        if (IST%mH_ice(i,j,k) < mi_old(i,j,k)) &
          IST%rdg_mice(i,j,k) = IST%rdg_mice(i,j,k) + rdg_frac(i,j,k) * &
             (IST%mH_ice(i,j,k) - mi_old(i,j,k)) * IST%part_size(i,j,k)
        IST%rdg_mice(i,j,k) = max(IST%rdg_mice(i,j,k), 0.0)
      enddo ; enddo ; enddo
    endif

    !  Sea-ice age ... changes due to growth and melt of ice volume and aging (time stepping)
    if (IST%id_age>0) call ice_aging(G, IST%mH_ice, IST%age_ice, mi_old, dt_slow)
    !  Other routines that do thermodynamic vertical processes should be added here

    ! Set up the thermodynamic fluxes in the externally visible structure Ice.
    call set_ocean_top_fluxes(Ice, IST, G)
    call accumulate_bottom_input(IST, Ice, dt_slow, G, IST%sum_output_CSp)

    if (IST%column_check) &
      call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IST%sum_output_CSp, &
                                message="      Post_thermo A", check_column=.true.)
    call adjust_ice_categories(IST%mH_ice, IST%mH_snow, IST%part_size, &
                               IST%TrReg, G, IST%ice_transport_CSp) !Niki: add ridging?
    call pass_var(IST%part_size, G%Domain)
    call pass_var(IST%mH_ice, G%Domain, complete=.false.)
    call pass_var(IST%mH_snow, G%Domain, complete=.true.)

    if (IST%column_check) &
      call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IST%sum_output_CSp, &
                                message="      Post_thermo B ", check_column=.true.)
  endif

  if (IST%id_xprt>0) then
    ! Store values to determine the ice and snow mass change due to transport.
    h2o_chg_xprt(:,:) = 0.0
  endif

  do nds=1,ndyn_steps

    if (.not.IST%interspersed_thermo .or. nds>1) then
      ! Correct the wind stresses for changes in the fractional ice-coverage.
      ice_cover(:,:) = 0.0
      do k=1,ncat ; do j=jsd,jed ; do i=isd,ied
        ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
      enddo ; enddo ; enddo
      do j=jsd,jed ; do i=isd,ied
        ice_free(i,j) = IST%part_size(i,j,0)

        if (ice_cover(i,j) > ice_cover_in(i,j)) then
          WindStr_x_A(i,j) = ((ice_cover(i,j)-ice_cover_in(i,j))*IST%flux_u_top(i,j,0) + &
                              ice_cover_in(i,j)*WindStr_x_A_in(i,j)) / ice_cover(i,j)
          WindStr_y_A(i,j) = ((ice_cover(i,j)-ice_cover_in(i,j))*IST%flux_v_top(i,j,0) + &
                              ice_cover_in(i,j)*WindStr_y_A_in(i,j)) / ice_cover(i,j)
        elseif (ice_free(i,j) > ice_free_in(i,j)) then
          WindStr_x_ocn_A(i,j) = ((ice_free(i,j)-ice_free_in(i,j))*WindStr_x_A_in(i,j) + &
                              ice_free_in(i,j)*IST%flux_u_top(i,j,0)) / ice_free(i,j)
          WindStr_y_ocn_A(i,j) = ((ice_free(i,j)-ice_free_in(i,j))*WindStr_y_A_in(i,j) + &
                              ice_free_in(i,j)*IST%flux_v_top(i,j,0)) / ice_free(i,j)
        endif
      enddo ; enddo
    endif

    !
    ! Dynamics - update ice velocities.
    !
    call mpp_clock_begin(iceClock4)

    ms_sum(:,:) = 0.0 ; mi_sum(:,:) = 0.0
    do k=1,ncat ; do j=jsd,jed ; do i=isd,ied
      ms_sum(i,j) = ms_sum(i,j) + (G%H_to_kg_m2 * IST%mH_snow(i,j,k)) * IST%part_size(i,j,k)
      mi_sum(i,j) = mi_sum(i,j) + (G%H_to_kg_m2 * IST%mH_ice(i,j,k))  * IST%part_size(i,j,k)
    enddo ; enddo ; enddo

    ! In the dynamics code, only the ice velocities are changed, and the ice-ocean
    ! stresses are calculated.  The gravity wave dynamics (i.e. the continuity
    ! equation) are not included in the dynamics.  All of the thickness categories
    ! are merged together.
    if (IST%Cgrid_dyn) then
      if (IST%area_wtd_stress) then
        !   The j-loop extents here are larger than they would normally be in case
        ! the stresses are being passed to the ocean on a B-grid.
        do j=jsc-1,jec+1 ; do I=isc-1,iec
          weights = (G%areaT(i,j)*ice_cover(i,j) + G%areaT(i+1,j)*ice_cover(i+1,j))
          if (G%mask2dCu(I,j) * weights > 0.0) then ; I_wts = 1.0 / weights
            WindStr_x_Cu(I,j) = G%mask2dCu(I,j) * &
                (G%areaT(i,j) * ice_cover(i,j) * WindStr_x_A(i,j) + &
                 G%areaT(i+1,j)*ice_cover(i+1,j)*WindStr_x_A(i+1,j)) * I_wts
          else
            WindStr_x_Cu(I,j) = 0.0
          endif

          weights = (G%areaT(i,j)*ice_free(i,j) + G%areaT(i+1,j)*ice_free(i+1,j))
          if (G%mask2dCu(I,j) * weights > 0.0) then ; I_wts = 1.0 / weights
            WindStr_x_ocn_Cu(I,j) = G%mask2dCu(I,j) * &
                (G%areaT(i,j) * ice_free(i,j) * WindStr_x_ocn_A(i,j) + &
                 G%areaT(i+1,j)*ice_free(i+1,j)*WindStr_x_ocn_A(i+1,j)) * I_wts
          else
            WindStr_x_ocn_Cu(I,j) = 0.0
          endif
        enddo ; enddo

        do J=jsc-1,jec ; do i=isc-1,iec+1
          weights = (G%areaT(i,j)*ice_cover(i,j) + G%areaT(i,j+1)*ice_cover(i,j+1))
          if (G%mask2dCv(i,J) * weights > 0.0) then ; I_wts = 1.0 / weights
            WindStr_y_Cv(i,J) = G%mask2dCv(i,J) * &
                (G%areaT(i,j) * ice_cover(i,j) * WindStr_y_A(i,j) + &
                 G%areaT(i,j+1)*ice_cover(i,j+1)*WindStr_y_A(i,j+1)) * I_wts
          else
            WindStr_y_Cv(i,J) = 0.0
          endif

          weights = (G%areaT(i,j)*ice_free(i,j) + G%areaT(i,j+1)*ice_free(i,j+1))
          if (weights > 0.0) then ; I_wts = 1.0 / weights
            WindStr_y_ocn_Cv(i,J) = G%mask2dCv(i,J) * &
                (G%areaT(i,j) * ice_free(i,j) * WindStr_y_ocn_A(i,j) + &
                 G%areaT(i,j+1)*ice_free(i,j+1)*WindStr_y_ocn_A(i,j+1)) * I_wts
          else
            WindStr_y_ocn_Cv(i,J) = 0.0
          endif
        enddo ; enddo
      else
        WindStr_x_Cu(:,:) = 0.0 ; WindStr_x_ocn_Cu(:,:) = 0.0 ; wts(:,:) = 0.0
        do k=1,ncat ; do j=jsc-1,jec+1 ; do I=isc-1,iec
          ps_vel = 0.5*G%mask2dCu(I,j) * (IST%part_size(i+1,j,k) + IST%part_size(i,j,k))
          WindStr_x_Cu(I,j) = WindStr_x_Cu(I,j) + ps_vel * (G%mask2dCu(I,j) * &
                           0.5* (IST%flux_u_top(i,j,k) + IST%flux_u_top(i+1,j,k)) )
          wts(I,J) = wts(I,J) + ps_vel
        enddo ; enddo ; enddo
        do j=jsc-1,jec+1 ; do I=isc-1,iec
          if (wts(I,j) > 0.) WindStr_x_Cu(I,j) = WindStr_x_Cu(I,j) / wts(I,j)

          WindStr_x_ocn_Cu(I,j) = G%mask2dCu(I,j) * &
                     0.5 * (IST%flux_u_top(i,j,0) + IST%flux_u_top(i+1,j,0))
        enddo ; enddo

        WindStr_y_Cv(:,:) = 0.0 ; WindStr_y_ocn_Cv(:,:) = 0.0 ; wts(:,:) = 0.0
        do k=1,ncat ; do J=jsc-1,jec ; do i=isc-1,iec+1
          ps_vel = 0.5*G%mask2dCv(i,J) * (IST%part_size(i,j+1,k) + IST%part_size(i,j,k))
          WindStr_y_Cv(i,j) = WindStr_y_Cv(i,J) + ps_vel * ( G%mask2dCv(i,J) * &
                         0.5*(IST%flux_v_top(i,j,k) + IST%flux_v_top(i,j+1,k)) )
          wts(i,J) = wts(i,J) + ps_vel
        enddo ; enddo ; enddo
        do J=jsc-1,jec ; do i=isc-1,iec+1
          if (wts(i,J) > 0.) WindStr_y_Cv(i,J) = WindStr_y_Cv(i,J) / wts(i,J)

          WindStr_y_ocn_Cv(i,J) = G%mask2dCv(i,J) * &
                     0.5*(IST%flux_v_top(i,j,0) + IST%flux_v_top(i,j+1,0))
        enddo ; enddo
      endif

      if (IST%debug) then
        call IST_chksum("Before ice_C_dynamics", IST, G)
        call hchksum(IST%part_size(:,:,0), "ps(0) before ice_C_dynamics", G)
        call hchksum(ms_sum, "ms_sum before ice_C_dynamics", G)
        call hchksum(mi_sum, "mi_sum before ice_C_dynamics", G)
        call hchksum(IST%sea_lev, "sea_lev before ice_C_dynamics", G, haloshift=1)
        call uchksum(IST%u_ocn_C, "u_ocn_C before ice_C_dynamics", G)
        call vchksum(IST%v_ocn_C, "v_ocn_C before ice_C_dynamics", G)
        call uchksum(WindStr_x_Cu, "WindStr_x_Cu before ice_C_dynamics", G)
        call vchksum(WindStr_y_Cv, "WindStr_y_Cv before ice_C_dynamics", G)
        call check_redundant_C("WindStr before ice_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G)
      endif

      call mpp_clock_begin(iceClocka)
      !### Ridging needs to be added with C-grid dynamics.
      call ice_C_dynamics(1.0-IST%part_size(:,:,0), ms_sum, mi_sum, IST%u_ice_C, IST%v_ice_C, &
                        IST%u_ocn_filt, IST%v_ocn_filt, &
                        WindStr_x_Cu, WindStr_y_Cv, IST%sea_lev_filt, str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, &
                        dt_slow_dyn, G, IST%ice_C_dyn_CSp)
      call mpp_clock_end(iceClocka)

      if (IST%debug) then
        call IST_chksum("After ice_dynamics", IST, G)
      endif

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
      call pass_vector(str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, G%Domain, stagger=CGRID_NE)
      call mpp_clock_end(iceClockb)
      !
      ! Dynamics diagnostics
      !
      if (IST%id_fax>0) call post_data(IST%id_fax, WindStr_x_Cu, IST%diag)
      if (IST%id_fay>0) call post_data(IST%id_fay, WindStr_y_Cv, IST%diag)

      call set_ocean_top_stress_Cgrid(Ice, IST, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, &
                                      str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, IST%part_size, G)
      call mpp_clock_end(iceClockc)

    else ! B-grid dynamics.

      if (IST%area_wtd_stress) then
        do J=jsc-1,jec ; do I=isc-1,iec ; if (G%mask2dBu(I,J) > 0.0) then
          weights = ((G%areaT(i+1,j+1)*ice_cover(i+1,j+1) + G%areaT(i,j)*ice_cover(i,j)) + &
                     (G%areaT(i+1,j)*ice_cover(i+1,j) + G%areaT(i,j+1)*ice_cover(i,j+1)) )
          I_wts = 0.0 ; if (weights > 0.0) I_wts = 1.0 / weights
          WindStr_x_B(I,J) = G%mask2dBu(I,J) * &
                  ((G%areaT(i+1,j+1)*ice_cover(i+1,j+1)*WindStr_x_A(i+1,j+1) + &
                    G%areaT(i,j)   * ice_cover(i,j)   * WindStr_x_A(i,j)) + &
                   (G%areaT(i+1,j) * ice_cover(i+1,j) * WindStr_x_A(i+1,j) + &
                    G%areaT(i,j+1) * ice_cover(i,j+1) * WindStr_x_A(i,j+1)) ) * I_wts
          WindStr_y_B(I,J) = G%mask2dBu(I,J) * &
                  ((G%areaT(i+1,j+1)*ice_cover(i+1,j+1)*WindStr_y_A(i+1,j+1) + &
                    G%areaT(i,j)   * ice_cover(i,j)   * WindStr_y_A(i,j)) + &
                   (G%areaT(i+1,j) * ice_cover(i+1,j) * WindStr_y_A(i+1,j) + &
                    G%areaT(i,j+1) * ice_cover(i,j+1) * WindStr_y_A(i,j+1)) ) * I_wts


          weights = ((G%areaT(i+1,j+1)*ice_free(i+1,j+1) + G%areaT(i,j)*ice_free(i,j)) + &
                     (G%areaT(i+1,j)*ice_free(i+1,j) + G%areaT(i,j+1)*ice_free(i,j+1)) )
          I_wts = 0.0 ; if (weights > 0.0) I_wts = 1.0 / weights
          WindStr_x_ocn_B(I,J) = G%mask2dBu(I,J) * &
                  ((G%areaT(i+1,j+1)*ice_free(i+1,j+1)*WindStr_x_ocn_A(i+1,j+1) + &
                    G%areaT(i,j)   * ice_free(i,j)   * WindStr_x_ocn_A(i,j)) + &
                   (G%areaT(i+1,j) * ice_free(i+1,j) * WindStr_x_ocn_A(i+1,j) + &
                    G%areaT(i,j+1) * ice_free(i,j+1) * WindStr_x_ocn_A(i,j+1)) ) * I_wts
          WindStr_y_ocn_B(I,J) = G%mask2dBu(I,J) * &
                  ((G%areaT(i+1,j+1)*ice_free(i+1,j+1)*WindStr_y_ocn_A(i+1,j+1) + &
                    G%areaT(i,j)   * ice_free(i,j)   * WindStr_y_ocn_A(i,j)) + &
                   (G%areaT(i+1,j) * ice_free(i+1,j) * WindStr_y_ocn_A(i+1,j) + &
                    G%areaT(i,j+1) * ice_free(i,j+1) * WindStr_y_ocn_A(i,j+1)) ) * I_wts
        else
          WindStr_x_B(I,J) = 0.0 ; WindStr_y_B(I,J) = 0.0
          WindStr_x_ocn_B(I,J) = 0.0 ; WindStr_y_ocn_B(I,J) = 0.0
        endif ; enddo ; enddo
      else
        WindStr_x_B(:,:) = 0.0 ; WindStr_y_B(:,:) = 0.0 ! ; wts(:,:) = 0.0
        do k=1,ncat ; do J=jsc-1,jec ; do I=isc-1,iec
          ps_vel = 0.25*G%mask2dBu(I,J) * &
              ((IST%part_size(i+1,j+1,k) + IST%part_size(i,j,k)) + &
               (IST%part_size(i+1,j,k) + IST%part_size(i,j+1,k)) )
          WindStr_x_B(I,J) = WindStr_x_B(I,J) + ps_vel * 0.25*( &
                  (IST%flux_u_top(i+1,j+1,k) + IST%flux_u_top(i,j,k)) + &
                  (IST%flux_u_top(i+1,j,k) + IST%flux_u_top(i,j+1,k)) )
          WindStr_y_B(I,J) = WindStr_y_B(I,J) + ps_vel * 0.25*( &
                  (IST%flux_v_top(i+1,j+1,k) + IST%flux_v_top(i,j,k)) + &
                  (IST%flux_v_top(i+1,j,k) + IST%flux_v_top(i,j+1,k)) )
          wts(I,J) = wts(I,J) + ps_vel
        enddo ; enddo ; enddo
        do J=jsc-1,jec ; do I=isc-1,iec
          if (wts(i,j) > 0.) then
            WindStr_x_B(I,J) = WindStr_x_B(I,J) / wts(I,J)
            WindStr_y_B(I,J) = WindStr_y_B(I,J) / wts(I,J)
          endif
          WindStr_x_ocn_B(I,J) = G%mask2dBu(I,J) * 0.25*( &
                  (IST%flux_u_top(i+1,j+1,0) + IST%flux_u_top(i,j,0)) + &
                  (IST%flux_u_top(i+1,j,0) + IST%flux_u_top(i,j+1,0)) )
          WindStr_y_ocn_B(I,J) = G%mask2dBu(I,J) * 0.25*( &
                  (IST%flux_v_top(i+1,j+1,0) + IST%flux_v_top(i,j,0)) + &
                  (IST%flux_v_top(i+1,j,0) + IST%flux_v_top(i,j+1,0)) )
        enddo ; enddo
      endif

      if (IST%debug) then
        call IST_chksum("Before ice_dynamics", IST, G)
        call hchksum(IST%part_size(:,:,0), "ps(0) before ice_dynamics", G)
        call hchksum(ms_sum, "ms_sum before ice_dynamics", G)
        call hchksum(mi_sum, "mi_sum before ice_dynamics", G)
        call hchksum(IST%sea_lev, "sea_lev before ice_dynamics", G, haloshift=1)
        call Bchksum(IST%u_ocn, "u_ocn before ice_dynamics", G, symmetric=.true.)
        call Bchksum(IST%v_ocn, "v_ocn before ice_dynamics", G, symmetric=.true.)
        call Bchksum(WindStr_x_B, "WindStr_x_B before ice_dynamics", G, symmetric=.true.)
        call Bchksum(WindStr_y_B, "WindStr_y_B before ice_dynamics", G, symmetric=.true.)
        call check_redundant_B("WindStr before ice_dynamics",WindStr_x_B, WindStr_y_B, G)
      endif

      rdg_rate(:,:) = 0.0
      call mpp_clock_begin(iceClocka)
      call ice_B_dynamics(1.0-IST%part_size(:,:,0), ms_sum, mi_sum, IST%u_ice_B, IST%v_ice_B, &
                        IST%u_ocn_filt, IST%v_ocn_filt, WindStr_x_B, WindStr_y_B, IST%sea_lev_filt, &
                        str_x_ice_ocn_B, str_y_ice_ocn_B, IST%do_ridging, &
                        rdg_rate(isc:iec,jsc:jec), dt_slow_dyn, G, IST%ice_B_dyn_CSp)
      call mpp_clock_end(iceClocka)

      if (IST%debug) then
        call IST_chksum("After ice_dynamics", IST, G)
      endif

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
      call mpp_clock_end(iceClockb)

      call mpp_clock_begin(iceClockc)
      !
      ! Dynamics diagnostics
      !
      if ((IST%id_fax>0) .or. (IST%id_fay>0)) then
        do J=jsc-1,jec ; do I=isc-1,iec
          ps_vel = (1.0 - G%mask2dBu(I,J)) + 0.25*G%mask2dBu(I,J) * &
                ((IST%part_size(i+1,j+1,0) + IST%part_size(i,j,0)) + &
                 (IST%part_size(i+1,j,0) + IST%part_size(i,j+1,0)) )
          diagVarBx(I,J) = ps_vel *  WindStr_x_ocn_B(I,J) + &
                           (1.0-ps_vel) * WindStr_x_B(I,J)
          diagVarBy(I,J) = ps_vel * WindStr_y_ocn_B(I,J) + &
                           (1.0-ps_vel) * WindStr_y_B(I,J)
        enddo ; enddo

        if (IST%id_fax>0) call post_data(IST%id_fax, diagVarBx, IST%diag)
        if (IST%id_fay>0) call post_data(IST%id_fay, diagVarBy, IST%diag)
      endif

      call set_ocean_top_stress_Bgrid(Ice, IST, WindStr_x_ocn_B, WindStr_y_ocn_B, &
                                      str_x_ice_ocn_B, str_y_ice_ocn_B, IST%part_size, G)
      call mpp_clock_end(iceClockc)
    endif ! End of B-grid dynamics


    call mpp_clock_end(iceClock4)

    !
    ! Thermodynamics (The thermodynamic changes might have been applied above.)
    !
    if (IST%interspersed_thermo .and. nds==1) then

      !TOM> Store old ice mass per unit area for calculating partial ice growth.
      do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
        mi_old(i,j,k) = IST%mH_ice(i,j,k)
      enddo ; enddo ; enddo
      !TOM> derive ridged ice fraction prior to thermodynamic changes of ice thickness
      !     in order to subtract ice melt proportionally from ridged ice volume (see below)
      if (IST%do_ridging) then
        do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
          tmp3 = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
          rdg_frac(i,j,k) = 0.0 ; if (tmp3 > 0.0) &
              rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp3
        enddo ; enddo ; enddo
      endif

      call disable_SIS_averaging(IST%diag)

      ! The thermodynamics routines return updated values of the ice and snow
      ! masses-per-unit area and enthalpies.
      call accumulate_input_2(IST, Ice, IST%part_size, dt_slow, G, IST%sum_output_CSp)
      if (IST%SIS1_5L_thermo) then
        call SIS1_5L_thermodynamics(Ice, IST, G) !, runoff, calving, runoff_hflx, calving_hflx)
      else
        call SIS2_thermodynamics(Ice, IST, G) !, runoff, calving, runoff_hflx, calving_hflx)
      endif

      !TOM> calculate partial ice growth for ridging and aging.
      if (IST%do_ridging) then
        !     ice growth (Ice%mH_ice > mi_old) does not affect ridged ice volume
        !     ice melt   (ice%mH_ice < mi_old) reduces ridged ice volume proportionally
        do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
          if (IST%mH_ice(i,j,k) < mi_old(i,j,k)) &
            IST%rdg_mice(i,j,k) = IST%rdg_mice(i,j,k) + rdg_frac(i,j,k) * &
               (IST%mH_ice(i,j,k) - mi_old(i,j,k)) * IST%part_size(i,j,k)
          IST%rdg_mice(i,j,k) = max(IST%rdg_mice(i,j,k), 0.0)
        enddo ; enddo ; enddo
      endif

      !  Sea-ice age ... changes due to growth and melt of ice volume and aging (time stepping)
      if (IST%id_age>0) call ice_aging(G, IST%mH_ice, IST%age_ice, mi_old, dt_slow)
      !  Other routines that do thermodynamic vertical processes should be added here.

      call adjust_ice_categories(IST%mH_ice, IST%mH_snow, IST%part_size, &
                                 IST%TrReg, G, IST%ice_transport_CSp) !Niki: add ridging?

      ! Set up the thermodynamic fluxes in the externally visible structure Ice.
      call set_ocean_top_fluxes(Ice, IST, G)
      call accumulate_bottom_input(IST, Ice, dt_slow, G, IST%sum_output_CSp)

      if (IST%column_check) &
        call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IST%sum_output_CSp, &
                                  message="        Post_thermo", check_column=.true.)
    endif  ! Interspersed thermo


    call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

    !
    ! Do ice transport ... all ocean fluxes have been calculated by now.
    !
    call mpp_clock_begin(iceClock8)

    if (IST%id_xprt>0) then ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      h2o_chg_xprt(i,j) = h2o_chg_xprt(i,j) - IST%part_size(i,j,k) * &
                        G%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
    enddo ; enddo ; enddo ; endif

    if (IST%debug) then
      call IST_chksum("Before ice_transport", IST, G)
    endif

    if (IST%Cgrid_dyn) then
      call ice_transport(IST%part_size, IST%mH_ice, IST%mH_snow, IST%u_ice_C, IST%v_ice_C, &
                         IST%TrReg, IST%sea_lev, dt_slow_dyn, G, IST%ice_transport_CSp, &
                         IST%rdg_mice, IST%age_ice(:,:,:,1), snow2ocn, rdg_rate, &
                         rdg_open, rdg_vosh)
    else
      ! B-grid transport
      ! Convert the velocities to C-grid points for transport.
      uc(:,:) = 0.0; vc(:,:) = 0.0
      do j=jsc,jec ; do I=isc-1,iec
        uc(I,j) = 0.5 * ( IST%u_ice_B(I,J-1) + IST%u_ice_B(I,J) )
      enddo ; enddo
      do J=jsc-1,jec ; do i = isc,iec
        vc(i,J) = 0.5 * ( IST%v_ice_B(I-1,J) + IST%v_ice_B(I,J) )
      enddo ; enddo

      call ice_transport(IST%part_size, IST%mH_ice, IST%mH_snow, uc, vc, &
                         IST%TrReg, IST%sea_lev, dt_slow_dyn, G, IST%ice_transport_CSp, &
                         IST%rdg_mice, IST%age_ice(:,:,:,1), snow2ocn, rdg_rate, &
                         rdg_open, rdg_vosh)
    endif
    if (IST%column_check) &
      call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IST%sum_output_CSp, &
                                message="      Post_transport")! , check_column=.true.)

    if (IST%id_xprt>0) then ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      h2o_chg_xprt(i,j) = h2o_chg_xprt(i,j) + IST%part_size(i,j,k) * &
                        G%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
    enddo ; enddo ; enddo ; endif

    call mpp_clock_end(iceClock8)

  enddo ! nds=1,ndyn_steps

  ! Add snow volume dumped into ocean to flux of frozen precipitation:
  !### WARNING - rdg_s2o is never calculated!!!
!  if (IST%do_ridging) then ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
!    IST%fprec_top(i,j,k) = IST%fprec_top(i,j,k) + rdg_s2o(i,j)*(IST%Rho_snow/dt_slow)
!  enddo ; enddo ; enddo ; endif

  call mpp_clock_begin(iceClock8)

  call finish_ocean_top_stresses(Ice, IST, G)

  ! Set appropriate surface quantities in categories with no ice.
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%part_size(i,j,k)<1e-10) &
    IST%t_surf(i,j,k) = T_0degC + T_Freeze(IST%s_surf(i,j),IST%ITV)
  enddo ; enddo ; enddo

  if (IST%bounds_check) call IST_bounds_check(IST, G, "After ice_transport")
  if (IST%debug) call IST_chksum("After ice_transport", IST, G)

  ! Sum the concentration weighted mass.
  mass(:,:) = 0.0
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    mass(i,j) = mass(i,j) + (G%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))) * &
                IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  if (IST%id_mi>0) call post_data(IST%id_mi, mass(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))

  if (IST%do_icebergs) call icebergs_incr_mass(Ice%icebergs, mass(isc:iec,jsc:jec)) ! Add icebergs mass in kg/m^2

  if (IST%id_mib>0) call post_data(IST%id_mib, mass(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec)) ! Diagnose total mass
  if (IST%id_slp>0) call post_data(IST%id_slp, Ice%p_surf, IST%diag, mask=Ice%mask)


  if (IST%specified_ice) then   ! over-write changes with specifications.
    h_ice_input(:,:) = 0.0
    call get_sea_surface(IST%Time, IST%t_surf(isc:iec,jsc:jec,0), IST%part_size(isc:iec,jsc:jec,:), &
                         h_ice_input(isc:iec,jsc:jec))
    do j=jsc,jec ; do i=isc,iec
      IST%mH_ice(i,j,1) = h_ice_input(i,j) * (G%kg_m2_to_H * IST%Rho_ice)
    enddo ; enddo
    call pass_var(IST%part_size, G%Domain)
  endif

  call mpp_clock_end(iceClock8)

  !
  ! Thermodynamic state diagnostics
  !
  call mpp_clock_begin(iceClock9)
  if (IST%id_cn>0) call post_data(IST%id_cn, IST%part_size(:,:,1:ncat), IST%diag, &
                                  mask=spread(Ice%G%Lmask2dT,3,ncat) )
  ! TK Mod: 10/18/02
  !  if (IST%id_obs_cn>0) call post_data(IST%id_obs_cn, Obs_cn_ice(:,:,2), IST%diag, &
  !                                       mask=G%Lmask2dT(isc:iec,jsc:jec)       )
  ! TK Mod: 10/18/02: (commented out...does not compile yet... add later
  !  if (IST%id_obs_hi>0) &
  !    call post_avg(IST%id_obs_hi, Obs_h_ice(isc:iec,jsc:jec), IST%part_size(isc:iec,jsc:jec,1:), &
  !                  IST%diag, G=G, mask=G%Lmask2dT(isc:iec,jsc:jec), wtd=.true.)

  !   Convert from ice and snow enthalpy back to temperature for diagnostic purposes.
  do_temp_diags = (IST%id_tsn > 0)
  do m=1,G%NkIce ; if (IST%id_t(m)>0) do_temp_diags = .true. ; enddo
  do_temp_diags = .true.  !### DELETE THIS WHEN T_ICE BECOMES A READ-ONLY RESTART VARIABLE.
  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             specified_thermo_salinity=spec_thermo_sal)
  I_enth_units = 1.0 / enth_units

  if (do_temp_diags) then
    do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
      if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
        if (spec_thermo_sal) then ; do m=1,G%NkIce
          temp_ice(i,j,k,m) = temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV)
        enddo ; else ; do m=1,G%NkIce
          temp_ice(i,j,k,m) = temp_from_En_S(IST%enth_ice(i,j,k,m), &
                                              IST%sal_ice(i,j,k,m), IST%ITV)
        enddo ; endif
      else
        do m=1,G%NkIce ; temp_ice(i,j,k,m) = 0.0 ; enddo
      endif
      if (IST%part_size(i,j,k)*IST%mH_snow(i,j,k) > 0.0) then
        temp_snow(i,j,k) = temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
      else
        temp_snow(i,j,k) = 0.0 ! ### Should this be = temp_ice(i,j,k,1)?
      endif
    enddo ; enddo ; enddo
  endif

  if (IST%id_ext>0) then
    do j=jsc,jec ; do i=isc,iec
      diagVar(i,j) = 0.0 ; if (IST%part_size(i,j,0) < 0.85) diagVar(i,j) = 1.0
    enddo ; enddo
    call post_data(IST%id_ext, diagVar(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  if (IST%id_hs>0) call post_avg(IST%id_hs, IST%mH_snow, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, &
                                 scale=G%H_to_kg_m2/IST%Rho_snow, wtd=.true.)
  if (IST%id_hi>0) call post_avg(IST%id_hi, IST%mH_ice, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, &
                                 scale=G%H_to_kg_m2/IST%Rho_ice, wtd=.true.)
  if (IST%id_ts>0) call post_avg(IST%id_ts, IST%t_surf(:,:,1:), IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, offset=-T_0degC, wtd=.true.)
  if (IST%id_tsn>0) call post_avg(IST%id_tsn, temp_snow, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  do m=1,G%NkIce
    if (IST%id_t(m)>0) call post_avg(IST%id_t(m), temp_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
    if (IST%id_sal(m)>0) call post_avg(IST%id_sal(m), IST%sal_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  enddo
  if (IST%id_t_iceav>0) call post_avg(IST%id_t_iceav, temp_ice, IST%part_size(:,:,1:), &
                                    IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_S_iceav>0) call post_avg(IST%id_S_iceav, IST%sal_ice, IST%part_size(:,:,1:), &
                                    IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_age>0) call post_avg(IST%id_age, IST%age_ice, IST%part_size(:,:,1:), &
                                  IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)


  if (IST%id_xprt>0) then
    call post_data(IST%id_xprt, h2o_chg_xprt(isc:iec,jsc:jec)*864e2*365/dt_slow, &
                   IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  if (IST%id_e2m>0) then
    tmp2d(:,:) = 0.0
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)>0.0) then
      tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k)*IST%mH_snow(i,j,k)*G%H_to_kg_m2 * &
                       ((enthalpy_liquid_freeze(0.0, IST%ITV) - &
                         IST%enth_snow(i,j,k,1)) * I_enth_units)
      if (spec_thermo_sal) then ; do m=1,NkIce
        tmp2d(i,j) = tmp2d(i,j) + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*G%H_to_kg_m2*I_Nk) * &
                       ((enthalpy_liquid_freeze(S_col(m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
      enddo ; else ; do m=1,NkIce
        tmp2d(i,j) = tmp2d(i,j) + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*G%H_to_kg_m2*I_Nk) * &
                       ((enthalpy_liquid_freeze(IST%sal_ice(i,j,k,m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
      enddo ; endif
    endif ; enddo ; enddo ; enddo
    call post_data(IST%id_e2m,  tmp2d(:,:), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  call disable_SIS_averaging(IST%diag)

  !
  ! Ridging diagnostics
  !
  !TOM> preparing output field fraction of ridged ice rdg_frac = (ridged ice volume) / (total ice volume)
  !     in each category; Ice%rdg_mice is ridged ice mass per unit total
  !     area throughout the code.
  if (IST%do_ridging) then
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      tmp3 = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
      if (tmp3*G%H_to_kg_m2 > IST%Rho_Ice*1.e-5) then   ! 1 mm ice thickness x 1% ice concentration
        rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp3
      else
        rdg_frac(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo

    call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

    if (IST%id_rdgr>0) call post_data(IST%id_rdgr, rdg_rate(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
!    if (id_rdgf>0) sent = send_data(id_rdgf,     rdg_frac(isc:iec,jsc:jec,2:km), Ice%Time, mask=spread(Ice%mask,3,km-1))
!    if (id_rdgo>0) sent = send_data(id_rdgo,     rdg_open(isc:iec,jsc:jec),      Ice%Time, mask=Ice%mask)
!    if (id_rdgv>0) sent = send_data(id_rdgv,     rdg_vosh(isc:iec,jsc:jec)*cell_area(isc:iec,jsc:jec), &
  endif

  !   Copy the surface properties, fractional areas and other variables to the
  ! externally visible structure Ice.
  !   Ice and IST may use different indexing conventions.
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    if (IST%slp2ocean) then
      Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) - 1e5 ! SLP - 1 std. atmosphere, in Pa.
    else
      Ice%p_surf(i2,j2) = 0.0
    endif
    Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) + grav*mass(i,j)

    Ice%mi(i2,j2) = mass(i,j)
  enddo ; enddo
  do k=0,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
  ! This may not be needed here?
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
  enddo ; enddo ; enddo

  if (IST%verbose) then
    call get_date(IST%Time, iyr, imon, iday, ihr, imin, isec)
    call get_time(IST%Time-set_date(iyr,1,1,0,0,0),isec,iday)
    call ice_line(iyr, iday+1, isec, IST%part_size(isc:iec,jsc:jec,0), &
                              IST%t_surf(:,:,0)-T_0degC, G)
  endif

  call mpp_clock_end(iceClock9)

  if (IST%debug) then
    call IST_chksum("End UIMS", IST, G)
    call Ice_public_type_chksum("End UIMS", Ice)
  endif

  if (IST%bounds_check) then
    call IST_bounds_check(IST, G, "End of update_ice_slow")
    call Ice_public_type_bounds_check(Ice, G, "End update_ice_slow")
  endif

  if (IST%Time + (IST%Time_step_slow/2) > IST%write_ice_stats_time) then
    call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IST%sum_output_CSp)
    IST%write_ice_stats_time = IST%write_ice_stats_time + IST%ice_stats_interval
  elseif (IST%column_check) then
    call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IST%sum_output_CSp)
  endif

end subroutine update_ice_model_slow

subroutine SIS1_5L_thermodynamics(Ice, IST, G) !, runoff, calving, &
                               ! runoff_hflx, calving_hflx)
  type(ice_data_type),                intent(inout) :: Ice
  type(ice_state_type),               intent(inout) :: IST
  type(sea_ice_grid_type),            intent(inout) :: G
!  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: runoff, calving
!  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: runoff_hflx, calving_hflx

  ! This subroutine does the thermodynamic calculations following SIS1.

  real, dimension(G%isc:G%iec,G%jsc:G%jec)  :: mi_change, h2o_change, bsnk, tmp2d
  real, dimension(SZI_(G),SZJ_(G),1:G%CatIce) :: &
    snow_to_ice, &   ! The conversion from snow to ice in m.
    h_ice, &         ! The ice thickness in m.
    h_snow           ! The snow thickness in m.
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: dum1, Obs_h_ice ! for qflux calculation
  real, dimension(G%isc:G%iec,G%jsc:G%jec,2) :: Obs_cn_ice      ! partition 2 = ice concentration
  real, dimension(SZI_(G),SZJ_(G))   :: qflx_lim_ice, qflx_res_ice
  real, dimension(1:G%CatIce)        :: e2m

  real, dimension(0:G%NkIce) :: T_col ! The temperature of a column of ice and snow in degC.
  real, dimension(G%NkIce) :: S_col ! The thermodynamic salinity of a column of ice, in g/kg.
  real :: heat_fill_val   ! A value of enthalpy to use for massless categories.
  real :: dt_slow, Idt_slow, yr_dtslow
  integer :: i, j, k, l, m, n, isc, iec, jsc, jec, ncat, NkIce
  integer :: i2, j2, k2, i_off, j_off
  real :: heat_to_ocn, h2o_to_ocn, evap_from_ocn, sn2ic, bablt
  real            :: heat_limit_ice, heat_res_ice
  real            :: tot_heat, heating, tot_frazil
  real :: T_Freeze_surf
  real :: H_to_m_ice, H_to_m_snow  ! The specific volumes of ice and snow times the
                               ! conversion factor from thickness units, in m H-1.
  real, parameter :: LI = hlf

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%runoff,1) - G%isc ; j_off = LBOUND(Ice%runoff,2) - G%jsc
  NkIce = G%NkIce
  dt_slow = time_type_to_real(IST%Time_step_slow) ; Idt_slow = 1.0/dt_slow
  H_to_m_ice = G%H_to_kg_m2 / IST%Rho_Ice ; H_to_m_snow = G%H_to_kg_m2 / IST%Rho_Snow

  if (G%NkIce /= 4) call SIS_error(FATAL, "SIS1_5L_thermodynamics requires that NK_ICE=4.")

  call mpp_clock_begin(iceClock5)

  snow_to_ice(:,:,:) = 0.0
  bsnk(:,:) = 0.0

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col)

  mi_change(:,:) = 0.0
  h2o_change(:,:) = 0.0
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    mi_change(i,j) = mi_change(i,J) - G%H_to_kg_m2*IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
    h2o_change(i,j) = h2o_change(i,j) - IST%part_size(i,j,k) * &
                      G%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
    h_ice(i,j,k) = IST%mH_ice(i,j,k) * H_to_m_Ice
    h_snow(i,j,k) = IST%mH_snow(i,j,k) * H_to_m_Snow
  enddo ; enddo ; enddo

  ! Start accumulating certain fluxes at the ocean's surface.
  do j=jsc,jec ; do i=isc,iec
    IST%flux_t_ocn_top(i,j) = IST%part_size(i,j,0) * IST%flux_t_top(i,j,0)
    IST%flux_q_ocn_top(i,j) = IST%part_size(i,j,0) * IST%flux_q_top(i,j,0)
    IST%flux_lw_ocn_top(i,j) = IST%part_size(i,j,0) * IST%flux_lw_top(i,j,0)
    IST%flux_lh_ocn_top(i,j) = IST%part_size(i,j,0) * IST%flux_lh_top(i,j,0)
    IST%flux_sw_vis_dir_ocn(i,j) = IST%part_size(i,j,0) * IST%flux_sw_vis_dir_top(i,j,0)
    IST%flux_sw_vis_dif_ocn(i,j) = IST%part_size(i,j,0) * IST%flux_sw_vis_dif_top(i,j,0)
    IST%flux_sw_nir_dir_ocn(i,j) = IST%part_size(i,j,0) * IST%flux_sw_nir_dir_top(i,j,0)
    IST%flux_sw_nir_dif_ocn(i,j) = IST%part_size(i,j,0) * IST%flux_sw_nir_dif_top(i,j,0)
    IST%lprec_ocn_top(i,j) = IST%part_size(i,j,0) * IST%lprec_top(i,j,0)
    IST%fprec_ocn_top(i,j) = IST%part_size(i,j,0) * IST%fprec_top(i,j,0)
  enddo ; enddo
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + IST%part_size(i,j,k) * IST%lprec_top(i,j,k)
  enddo ; enddo ; enddo
  if (IST%num_tr_fluxes>0) then ; do n=1,IST%num_tr_fluxes
    do j=jsc,jec ; do i=isc,iec
      IST%tr_flux_ocn_top(i,j,n) = IST%part_size(i,j,0) * IST%tr_flux_top(i,j,0,n)
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      IST%tr_flux_ocn_top(i,j,n) = IST%tr_flux_ocn_top(i,j,n) + &
                   IST%part_size(i,j,k) * IST%tr_flux_top(i,j,k,n)
    enddo ; enddo ; enddo
  enddo ; endif

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if (G%mask2dT(i,j) > 0 .and. IST%part_size(i,j,k) > 0) then
      T_col(0) = Temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
      do m=1,NkIce ; T_col(m) = Temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV) ; enddo

      ! reshape the ice based on fluxes
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j),IST%ITV)

      evap_from_ocn = 0.0; h2o_to_ocn = 0.0; heat_to_ocn = 0.0
      call ice5lay_resize(h_snow(i,j,k), T_col(0), h_ice(i,j,k),&
                          T_col(1), T_col(2), T_col(3), T_col(4),            &
                          IST%fprec_top(i,j,k) *dt_slow, 0.0,                &
                          IST%flux_q_top(i,j,k)*dt_slow,                     &
                          IST%tmelt (i,j,k), IST%bmelt(i,j,k),               &
                          T_Freeze_surf,                                     &
                          heat_to_ocn, h2o_to_ocn, evap_from_ocn,            &
                          snow_to_ice(i,j,k), bablt                          )

      IST%flux_q_ocn_top(i,j) = IST%flux_q_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                  (evap_from_ocn*Idt_slow) ! no ice, evaporation left
      IST%flux_lh_ocn_top(i,j) = IST%flux_lh_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                 ((hlv*evap_from_ocn)*Idt_slow)
      IST%flux_t_ocn_top(i,j) = IST%flux_t_ocn_top(i,j) + IST%part_size(i,j,k) * &
            (IST%bheat(i,j) - (heat_to_ocn - hlf*evap_from_ocn)*Idt_slow)
      IST%flux_sw_vis_dif_ocn(i,j) = IST%flux_sw_vis_dif_ocn(i,j) + IST%part_size(i,j,k) * &
             (((IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k)) + &
               (IST%flux_sw_nir_dir_top(i,j,k) + IST%flux_sw_nir_dif_top(i,j,k))) * &
              IST%sw_abs_ocn(i,j,k))
      IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + IST%part_size(i,j,k) * &
              (h2o_to_ocn*Idt_slow)

      bsnk(i,j) = bsnk(i,j) - IST%part_size(i,j,k)*bablt ! bot. melt. ablation

      IST%enth_snow(i,j,k,1) = enth_from_TS(T_col(0), 0.0, IST%ITV)
      do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_from_TS(T_col(m), S_col(m), IST%ITV) ; enddo

    endif

    !
    ! absorb frazil in thinest ice partition available
    !
    if (IST%frazil(i,j)>0 .and. IST%part_size(i,j,0)+IST%part_size(i,j,k)>0.01) then
      !                                                           was ...>0.0
      ! raised above threshold from 0 to 0.01 to avert ocean-ice model blow-ups
      !
      T_col(0) = Temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
      do m=1,NkIce ; T_col(m) = Temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV) ; enddo

      T_Freeze_surf = T_Freeze(IST%s_surf(i,j),IST%ITV)
      h_snow(i,j,k) = h_snow(i,j,k) * IST%part_size(i,j,k)
      h_ice(i,j,k)  = h_ice(i,j,k)  * IST%part_size(i,j,k)
      IST%t_surf(i,j,k) = (IST%t_surf(i,j,k) * IST%part_size(i,j,k) + &
                           (T_0degC + T_Freeze_surf)*IST%part_size(i,j,0))
      IST%part_size(i,j,k) = IST%part_size(i,j,k) + IST%part_size(i,j,0)
      IST%part_size(i,j,0) = 0.0
      h_snow(i,j,k) = h_snow(i,j,k) / IST%part_size(i,j,k)
      h_ice(i,j,k)  = h_ice(i,j,k) / IST%part_size(i,j,k)
      IST%t_surf(i,j,k) = IST%t_surf(i,j,k) / IST%part_size(i,j,k)

      call ice5lay_resize(h_snow(i,j,k), T_col(0), h_ice(i,j,k), &
                          T_col(1), T_col(2), T_col(3), T_col(4), 0.0,           &
                          IST%frazil(i,j) / IST%part_size(i,j,k), 0.0, 0.0, 0.0, &
                          T_Freeze_surf, &
                          heat_to_ocn, h2o_to_ocn, evap_from_ocn, sn2ic)
      IST%frazil(i,j) = 0.0;
      !
      ! spread frazil salinification over all partitions
      !
      IST%lprec_top(i,j,:) = IST%lprec_top(i,j,:) + h2o_to_ocn*IST%part_size(i,j,k)/dt_slow
      IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + &
                               (h2o_to_ocn*IST%part_size(i,j,k)) / dt_slow

      IST%enth_snow(i,j,k,1) = enth_from_TS(T_col(0), 0.0, IST%ITV)
      do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_from_TS(T_col(m), S_col(m), IST%ITV) ; enddo

    endif

  enddo ; enddo ; enddo   ! i-, j-, and k-loops
  call mpp_clock_end(iceClock5)

  !
  ! Calculate QFLUX's from (1) restoring to obs and (2) limiting total ice.
  !
  call mpp_clock_begin(iceClock6)
  ! get observed ice thickness for ice restoring, if calculating qflux
  if (IST%do_ice_restore) &
    call get_sea_surface(IST%Time, dum1, Obs_cn_ice, Obs_h_ice)

  if (IST%do_ice_restore .or. IST%do_ice_limit) then
    qflx_lim_ice(:,:) = 0.0
    qflx_res_ice(:,:) = 0.0

    do j=jsc,jec ; do i=isc,iec
      T_col(0) = Temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
      do m=1,NkIce ; T_col(m) = Temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV) ; enddo

      heat_res_ice   = 0.0
      heat_limit_ice = 0.0
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j),IST%ITV)
      !
      ! calculate enthalpy
      !
      if (IST%slab_ice) then
        e2m(1) = h_ice(i,j,1)*IST%Rho_ice*LI
      else
        do k=1,ncat
          if ((IST%part_size(i,j,k)>0.0 .and. h_ice(i,j,k)>0.0)) then
             e2m(k) = e_to_melt(h_snow(i,j,k), T_col(0), h_ice(i,j,k), &
                      T_col(1), T_col(2), T_col(3), T_col(4)) * IST%part_size(i,j,k)
          else
             e2m(k) = 0.0
          endif
        enddo
      endif
      !
      ! calculate heat needed to constrain ice enthalpy
      !
      if (IST%do_ice_restore) then
        ! Restore to observed enthalpy, implying restoring toward
        ! thickness * concentration.
        if (IST%slab_ice) then
          heat_res_ice = -(LI*IST%Rho_ice*Obs_h_ice(i,j)-sum(e2m)) &
                         *dt_slow/(86400*IST%ice_restore_timescale)
        else
          heat_res_ice = -(LI*IST%Rho_ice*Obs_h_ice(i,j)*Obs_cn_ice(i,j,2)-sum(e2m)) &
                         *dt_slow/(86400*IST%ice_restore_timescale)
        endif
      endif

      if (IST%do_ice_limit .and. (sum(e2m) > IST%max_ice_limit*IST%Rho_ice*LI)) then
        heat_limit_ice = sum(e2m)-LI*IST%Rho_ice*IST%max_ice_limit
        ! should we "heat_ice_res = 0.0" ?
      endif

      !
      ! apply constraining heat to ice
      !
      tot_heat = heat_res_ice+heat_limit_ice
      if (IST%slab_ice) h_ice(i,j,1) = h_ice(i,j,1) - tot_heat/(IST%Rho_ice*LI)

      if (.not. IST%slab_ice .and. (tot_heat>0.0)) then  ! add like ocean-ice heat
        do k=0,ncat-1
          if (e2m(k) > 0.0) then
            heating = tot_heat/sum(IST%part_size(i,j,k:ncat))
            if (heating*IST%part_size(i,j,k) > e2m(k)) then ! cat. melts away
              h_ice (i,j,k) = 0.0
              h_snow(i,j,k) = 0.0
              tot_heat = tot_heat - e2m(k)
            else
              evap_from_ocn = 0.0; h2o_to_ocn = 0.0; heat_to_ocn = 0.0
              call ice5lay_resize(h_snow(i,j,k), T_col(0), h_ice(i,j,k), &
                                  T_col(1), T_col(2), T_col(3), T_col(4), &
                                  0.0, 0.0, 0.0, 0.0, heating, T_Freeze_surf, &
                                  heat_to_ocn, h2o_to_ocn, evap_from_ocn, &
                                  snow_to_ice(i,j,k), bablt              )
              tot_heat = tot_heat - heating*IST%part_size(i,j,k)
            endif
          endif
        enddo
      endif

      tot_heat = heat_res_ice+heat_limit_ice
      if (.not. IST%slab_ice .and. (tot_heat<0.0)) then ! add like frazil
        do k=1,ncat
          if (IST%part_size(i,j,0)+IST%part_size(i,j,k)>0) exit
        enddo
        ! k is thinnest ice partition that can recieve frazil
        h_snow(i,j,k) = h_snow(i,j,k) * IST%part_size(i,j,k)
        h_ice(i,j,k)  = h_ice(i,j,k)  * IST%part_size(i,j,k)
        IST%t_surf(i,j,k) = IST%t_surf(i,j,k) * IST%part_size(i,j,k) &
                       + (T_0degC + T_Freeze_surf)*IST%part_size(i,j,0)
        IST%part_size(i,j,k) = IST%part_size(i,j,k) + IST%part_size(i,j,0)
        IST%part_size(i,j,0) = 0.0
        h_snow(i,j,k) = h_snow(i,j,k) / IST%part_size(i,j,k)
        h_ice(i,j,k) =  h_ice(i,j,k)  / IST%part_size(i,j,k)
        IST%t_surf(i,j,k) = IST%t_surf(i,j,k) / IST%part_size(i,j,k)

        call ice5lay_resize(h_snow(i,j,k), T_col(0), h_ice(i,j,k), &
                            T_col(1), T_col(2), T_col(3), T_col(4), 0.0,   &
                            -tot_heat/IST%part_size(i,j,k), 0.0, 0.0, 0.0, &
                            T_Freeze_surf,                        &
                            heat_to_ocn, h2o_to_ocn, evap_from_ocn, sn2ic)
      endif

      ! Convert constraining heat from energy (J/m^2) to flux (W/m^2)
      qflx_lim_ice(i,j) = heat_limit_ice / dt_slow
      qflx_res_ice(i,j) = heat_res_ice / dt_slow
      !
      ! Check for energy conservation
      !
      if (IST%slab_ice) then
        e2m(1) = e2m(1) - h_ice(i,j,1)*IST%Rho_ice*LI
      else
        do k=1,ncat
          if (IST%part_size(i,j,k)>0.0 .and. h_ice(i,j,k)>0.0) &
            e2m(k) = e2m(k)-e_to_melt(h_snow(i,j,k), T_col(0), h_ice(i,j,k), &
                     T_col(1), T_col(2), T_col(3), T_col(4)) * IST%part_size(i,j,k)
        enddo
      endif
      ! if (abs(sum(e2m) - heat_res_ice - heat_limit_ice)>IST%Rho_ice*LI*1e-3) &
      !       print *, 'QFLUX conservation error at', i, j, 'heat2ice=',  &
      !             tot_heat, 'melted=', sum(e2m), 'h*part_size=', &
      !             h_ice(i,j,:)*IST%part_size(i,j,:)

      IST%enth_snow(i,j,k,1) = enth_from_TS(T_col(0), 0.0, IST%ITV)
      do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_from_TS(T_col(m), S_col(m), IST%ITV) ; enddo

    enddo ; enddo
  endif ! End of (IST%do_ice_restore .or. IST%do_ice_limit) block
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    IST%mH_ice(i,j,k) = h_ice(i,j,k) * (G%kg_m2_to_H * IST%Rho_ice)
    IST%mH_snow(i,j,k) = h_snow(i,j,k) * (G%kg_m2_to_H * IST%Rho_snow)
  enddo ; enddo ; enddo

  !   Convert from ice temperature (which is not conserved) to enthalpy, which
  ! includes the heat requirements for melting of brine pockets associated with
  ! temperature changes.
  heat_fill_val = Enth_from_TS(0.0, 0.0, IST%ITV)

  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) <= 0.0) then
      do m=1,G%NkIce ; IST%enth_ice(i,j,k,m) = heat_fill_val ; enddo
      IST%enth_snow(i,j,k,1) = heat_fill_val
    endif
  enddo ; enddo ; enddo

  call mpp_clock_end(iceClock6)

  ! Determine the salt fluxes to ocean
  ! Note that at this point mi_change and h2o_change are the negative of the masses.
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    mi_change(i,j) = mi_change(i,J) + (G%H_to_kg_m2*IST%mH_ice(i,j,k))*IST%part_size(i,j,k)
    h2o_change(i,j) = h2o_change(i,j) + IST%part_size(i,j,k) * &
                      (G%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k)))
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    ! Note the conversion here from g m-2 to kg m-2 s-1.
    Ice%flux_salt(i2,j2) = (0.001*IST%ice_bulk_salin) * &
                           (mi_change(i,j) * Idt_slow)
  enddo ; enddo

  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

  yr_dtslow = (864e2*365/dt_slow)
  if (IST%id_lsnk>0) then
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = min(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(IST%id_lsnk, tmp2d(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  if (IST%id_lsrc>0) then
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = max(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(IST%id_lsrc, tmp2d(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  if (IST%id_saltf>0) call post_data(IST%id_saltf, Ice%flux_salt, IST%diag, &
                                     mask=G%Lmask2dT(isc:iec,jsc:jec))
  if (IST%id_bsnk>0)  call post_data(IST%id_bsnk, bsnk(isc:iec,jsc:jec)*yr_dtslow, &
                                     IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  if (IST%id_tmelt>0) call post_avg(IST%id_tmelt, IST%tmelt, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_bmelt>0) call post_avg(IST%id_bmelt, IST%bmelt, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_sn2ic>0) call post_avg(IST%id_sn2ic, snow_to_ice, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, mask=G%Lmask2dT)
  if (IST%id_qflim>0) call post_data(IST%id_qflim, qflx_lim_ice, IST%diag, mask=G%Lmask2dT)
  if (IST%id_qfres>0) call post_data(IST%id_qfres, qflx_res_ice, IST%diag, mask=G%Lmask2dT)

  call disable_SIS_averaging(IST%diag)

end subroutine SIS1_5L_thermodynamics

subroutine SIS2_thermodynamics(Ice, IST, G) !, runoff, calving, &
                               ! runoff_hflx, calving_hflx)
  type(ice_data_type),                intent(inout) :: Ice
  type(ice_state_type),               intent(inout) :: IST
  type(sea_ice_grid_type),            intent(inout) :: G
!  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: runoff, calving
!  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: runoff_hflx, calving_hflx

  ! This subroutine does the thermodynamic calculations in the same order as SIS1,
  ! but with a greater emphasis on enthalpy as the dominant state variable.

  real, dimension(G%isc:G%iec,G%jsc:G%jec)  :: salt_change, h2o_change, bsnk, tmp2d
  real, dimension(SZI_(G),SZJ_(G),1:G%CatIce) :: snow_to_ice
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: dum1, Obs_h_ice ! for qflux calculation
  real, dimension(G%isc:G%iec,G%jsc:G%jec,2) :: Obs_cn_ice      ! partition 2 = ice concentration
  real, dimension(SZI_(G),SZJ_(G))   :: qflx_lim_ice, qflx_res_ice
  real, dimension(SZI_(G),SZJ_(G),1:G%CatIce)   :: heat_in, enth_prev, enth
  real, dimension(SZI_(G),SZJ_(G))   :: heat_in_col, enth_prev_col, enth_col, enth_mass_in_col
  real, dimension(1:G%CatIce)        :: e2m
  real, dimension(G%NkIce) :: S_col        ! The salinity of a column of ice, in g/kg.
  real, dimension(G%NkIce+1) :: Salin      ! The conserved bulk salinity of each
                                           ! layer in g/kg, with the salinity of
                                           ! newly formed ice in layer NkIce+1.
  real, dimension(0:G%NkIce)   :: Tcol0    ! The temperature of a column of ice and snow, in degC.
  real, dimension(0:G%NkIce)   :: S_col0   ! The salinity of a column of ice and snow, in g/kg.
  real, dimension(0:G%NkIce)   :: Tfr_col0 ! The freezing temperature of a column of ice and snow, in degC.
  real, dimension(0:G%NkIce+1) :: &
    enthalpy              ! The initial enthalpy of a column of ice and snow
                          ! and the surface ocean, in enth_units (often J/kg).
  real :: enthalpy_ocean  ! The enthalpy of the ocean surface waters, in Enth_units.
  real :: heat_fill_val   ! An enthalpy to use for massless categories, in enth_units.

  real :: T_freeze_surf ! The freezing temperature at the surface salinity of
                        ! the ocean, in deg C.
  real :: I_part   ! The inverse of a part_size, nondim.
  logical :: spec_thermo_sal  ! If true, use the specified salinities of the
                              ! various sub-layers of the ice for all thermodynamic
                              ! calculations; otherwise use the prognostic
                              ! salinity fields for these calculations.

  real :: dt_slow, Idt_slow, yr_dtslow
  integer :: i, j, k, l, m, n, isc, iec, jsc, jec, ncat, NkIce
  integer :: i2, j2, k2, i_off, j_off
  real :: heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, sn2ic, bablt
  real :: salt_to_ice
  real :: heat_limit_ice, heat_res_ice
  real :: enth_evap, enth_ice_to_ocn, enth_ocn_to_ice, enth_snowfall
  real :: tot_heat, heating, tot_frazil, heat_mass_in, heat_input
  real :: mass_in, mass_here, mass_prev, mass_imb
  real :: enth_units, I_enth_units
  real :: I_Nk
  real :: kg_H_Nk  ! The conversion factor from units of H to kg/m2 over Nk.
  real, parameter :: LI = hlf

  real :: tot_heat_in, enth_here, enth_imb, norm_enth_imb, emic2, tot_heat_in2, enth_imb2

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%runoff,1) - G%isc ; j_off = LBOUND(Ice%runoff,2) - G%jsc
  NkIce = G%NkIce ; I_Nk = 1.0 / G%NkIce ; kg_H_Nk = G%H_to_kg_m2 * I_Nk
  dt_slow = time_type_to_real(IST%Time_step_slow) ; Idt_slow = 1.0/dt_slow

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             specified_thermo_salinity=spec_thermo_sal)
  S_col0(0) = 0.0 ; do m=1,G%NkIce ; S_col0(m) = S_col(m) ; enddo
  call calculate_T_Freeze(S_col0(0:G%NkIce), Tfr_col0(0:G%NkIce), IST%ITV)
  I_enth_units = 1.0 / enth_units

  heat_fill_val = Enth_from_TS(0.0, 0.0, IST%ITV)

  if (.not.spec_thermo_sal) call SIS_error(FATAL, "SIS2_thermodynamics is not "//&
    "prepared for SPECIFIED_THERMO_SALINITY to be false.")

  if (IST%column_check) then
    enth_prev(:,:,:) = 0.0 ; heat_in(:,:,:) = 0.0

    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
      enth_prev(i,j,k) = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
      do m=1,G%NkIce
        enth_prev(i,j,k) = enth_prev(i,j,k) + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
      enddo
    endif ; enddo ; enddo ; enddo

    enth_prev_col(:,:) = 0.0 ; heat_in_col(:,:) = 0.0 ; enth_mass_in_col(:,:) = 0.0

    do j=jsc,jec ; do i=isc,iec
      heat_in_col(i,j) = heat_in_col(i,j) - IST%frazil(i,j)
      heat_in_col(i,j) = heat_in_col(i,j) - IST%part_size(i,j,0) * dt_slow*IST%flux_t_top(i,j,0)
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%part_size(i,j,k) > 0.0) then
        heat_in_col(i,j) = heat_in_col(i,j) + IST%part_size(i,j,k) * &
          (IST%tmelt(i,j,k) + IST%bmelt(i,j,k) - dt_slow*IST%bheat(i,j))
    endif ; enddo ; enddo ; enddo

    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
      enth_prev_col(i,j) = enth_prev_col(i,j) + &
        (IST%mH_snow(i,j,k)*IST%part_size(i,j,k)) * IST%enth_snow(i,j,k,1)
      do m=1,G%NkIce
        enth_prev_col(i,j) = enth_prev_col(i,j) + &
          (IST%mH_ice(i,j,k)*IST%part_size(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
      enddo
    endif ; enddo ; enddo ; enddo
  endif

  call mpp_clock_begin(iceClock5)

  snow_to_ice(:,:,:) = 0.0
  bsnk(:,:) = 0.0

  salt_change(:,:) = 0.0
  h2o_change(:,:) = 0.0
  if (IST%ice_rel_salin <= 0.0) then
    do m=1,G%NkIce ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      salt_change(i,j) = salt_change(i,j) - &
         (IST%sal_ice(i,j,k,m)*(IST%mH_ice(i,j,k)*kg_H_Nk)) * IST%part_size(i,j,k)
    enddo ; enddo ; enddo ; enddo
  endif
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    h2o_change(i,j) = h2o_change(i,j) - IST%part_size(i,j,k) * &
                      G%H_to_kg_m2*(IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
  enddo ; enddo ; enddo

  ! Start accumulating the fluxes at the ocean's surface.
  do j=jsc,jec ; do i=isc,iec
    IST%flux_t_ocn_top(i,j) = IST%part_size(i,j,0) * IST%flux_t_top(i,j,0)
    IST%flux_q_ocn_top(i,j) = IST%part_size(i,j,0) * IST%flux_q_top(i,j,0)
    IST%flux_lw_ocn_top(i,j) = IST%part_size(i,j,0) * IST%flux_lw_top(i,j,0)
    IST%flux_lh_ocn_top(i,j) = IST%part_size(i,j,0) * IST%flux_lh_top(i,j,0)
    IST%flux_sw_vis_dir_ocn(i,j) = IST%part_size(i,j,0) * IST%flux_sw_vis_dir_top(i,j,0)
    IST%flux_sw_vis_dif_ocn(i,j) = IST%part_size(i,j,0) * IST%flux_sw_vis_dif_top(i,j,0)
    IST%flux_sw_nir_dir_ocn(i,j) = IST%part_size(i,j,0) * IST%flux_sw_nir_dir_top(i,j,0)
    IST%flux_sw_nir_dif_ocn(i,j) = IST%part_size(i,j,0) * IST%flux_sw_nir_dif_top(i,j,0)
    IST%lprec_ocn_top(i,j) = IST%part_size(i,j,0) * IST%lprec_top(i,j,0)
    IST%fprec_ocn_top(i,j) = IST%part_size(i,j,0) * IST%fprec_top(i,j,0)
  enddo ; enddo
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + IST%part_size(i,j,k) * IST%lprec_top(i,j,k)
  enddo ; enddo ; enddo

  if (IST%num_tr_fluxes>0) then ; do n=1,IST%num_tr_fluxes
    do j=jsc,jec ; do i=isc,iec
      IST%tr_flux_ocn_top(i,j,n) = IST%part_size(i,j,0) * IST%tr_flux_top(i,j,0,n)
    enddo ; enddo
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      IST%tr_flux_ocn_top(i,j,n) = IST%tr_flux_ocn_top(i,j,n) + &
                   IST%part_size(i,j,k) * IST%tr_flux_top(i,j,k,n)
    enddo ; enddo ; enddo
  enddo ; endif

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if (G%mask2dT(i,j) > 0 .and. IST%part_size(i,j,k) > 0) then
      ! reshape the ice based on fluxes
      if (IST%column_check) then
        mass_prev = IST%mH_snow(i,j,k)
        mass_prev = mass_prev + IST%mH_ice(i,j,k)
      endif

 !     evap_from_ocn = 0.0; h2o_ice_to_ocn = 0.0; heat_to_ocn = 0.0

      if (IST%mH_snow(i,j,k) == 0.0) IST%enth_snow(i,j,k,1) = &
          enth_from_TS(Temp_from_En_S(IST%enth_ice(i,j,k,1), S_col0(1), IST%ITV), &
                       0.0, IST%ITV)
      enthalpy(0) = IST%enth_snow(i,j,k,1)
      do m=1,NkIce ; enthalpy(m) = IST%enth_ice(i,j,k,m) ; enddo
      enthalpy_ocean = enthalpy_liquid(IST%t_ocn(i,j), IST%s_surf(i,j), IST%ITV)
      enthalpy(NkIce+1) = enthalpy_ocean

      if (IST%ice_rel_salin > 0.0) then
        do m=1,NkIce ; Salin(m) = IST%sal_ice(i,j,k,m) ; enddo
        salin(NkIce+1) = IST%ice_rel_salin * IST%s_surf(i,j)
      else
        do m=1,NkIce+1 ; Salin(m) = IST%ice_bulk_salin ; enddo
      endif

      call ice_resize_SIS2(IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), G%H_to_kg_m2, enthalpy, &
                   S_col, Salin, IST%fprec_top(i,j,k)*dt_slow, 0.0, &
                   IST%flux_q_top(i,j,k)*dt_slow, IST%tmelt(i,j,k), &
                   IST%bmelt(i,j,k), T_Freeze(IST%s_surf(i,j),IST%ITV), NkIce, &
                   heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                   snow_to_ice(i,j,k), salt_to_ice, IST%ITV, IST%ice_thm_CSp, bablt, &
                   Enthalpy_evap=enth_evap, Enthalpy_melt=enth_ice_to_ocn, &
                   Enthalpy_freeze=enth_ocn_to_ice)
      if (IST%mH_ice(i,j,k) == 0.0) then
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy_liquid_freeze(S_col(m), IST%ITV) ; enddo
      else
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy(m) ; enddo
      endif
      if (IST%ice_rel_salin > 0.0) then
        if (IST%mH_ice(i,j,k) == 0.0) then
          do m=1,NkIce ; IST%sal_ice(i,j,k,m) = 0.0 ; enddo
        else
          do m=1,NkIce ; IST%sal_ice(i,j,k,m) = Salin(m) ; enddo
        endif
        salt_change(i,j) = salt_change(i,j) + IST%part_size(i,j,k) * salt_to_ice
      endif

      ! The snow enthalpy should not have changed.  This should do nothing.
      ! IST%enth_snow(i,j,k,1) = Enthalpy(0)

      enth_snowfall = ((dt_slow*IST%fprec_top(i,j,k)) * enthalpy(0))
      IST%Enth_Mass_in_atm(i,j) = IST%Enth_Mass_in_atm(i,j) + &
           IST%part_size(i,j,k) * enth_snowfall

!      IST%Enth_Mass_in_ocn(i,j) = IST%Enth_Mass_in_ocn(i,j) + &
!          IST%part_size(i,j,k) * enth_ocn_to_ice
      IST%Enth_Mass_in_ocn(i,j) = IST%Enth_Mass_in_ocn(i,j) + &
          IST%part_size(i,j,k) * (h2o_ocn_to_ice * enthalpy_ocean)
      IST%Enth_Mass_out_ocn(i,j) = IST%Enth_Mass_out_ocn(i,j) - &
          IST%part_size(i,j,k) * enth_ice_to_ocn
      IST%Enth_Mass_out_atm(i,j) = IST%Enth_Mass_out_atm(i,j) - &
          IST%part_size(i,j,k) * enth_evap


      if (IST%column_check) then
        heat_in(i,j,k) = heat_in(i,j,k) + IST%tmelt(i,j,k) + IST%bmelt(i,j,k) - &
                     (heat_to_ocn - (hlv+hlf)*evap_from_ocn)

        heat_input = IST%tmelt(i,j,k) + IST%bmelt(i,j,k) - (heat_to_ocn - (hlv+hlf)*evap_from_ocn)
        heat_mass_in = enth_snowfall + enth_ocn_to_ice - enth_ice_to_ocn - enth_evap
        mass_in = dt_slow*IST%fprec_top(i,j,k) + h2o_ocn_to_ice - h2o_ice_to_ocn - &
                 (dt_slow*IST%flux_q_top(i,j,k)-evap_from_ocn)

        mass_here = IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k)
        enth_here = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,G%NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
        tot_heat_in = G%kg_m2_to_H*(enth_units*heat_input + heat_mass_in)
        mass_in = mass_in*G%kg_m2_to_H

        enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        mass_imb = mass_here - (mass_prev + mass_in)
        if (abs(enth_imb) > IST%imb_tol * &
            (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
          norm_enth_imb = enth_imb / (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in))
          enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        endif
      endif

      if (IST%column_check) then
       enth_mass_in_col(i,j) = enth_mass_in_col(i,j) + IST%part_size(i,j,k) * &
          (enth_snowfall + enth_ocn_to_ice - enth_ice_to_ocn - enth_evap)
      endif

      IST%flux_q_ocn_top(i,j) = IST%flux_q_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                  (evap_from_ocn*Idt_slow) ! no ice, evaporation left
      IST%flux_lh_ocn_top(i,j) = IST%flux_lh_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                 ((hlv*evap_from_ocn)*Idt_slow)
      IST%flux_t_ocn_top(i,j) = IST%flux_t_ocn_top(i,j) + IST%part_size(i,j,k) * &
             (IST%bheat(i,j) - (heat_to_ocn - hlf*evap_from_ocn)*Idt_slow)
      IST%flux_sw_vis_dif_ocn(i,j) = IST%flux_sw_vis_dif_ocn(i,j) + IST%part_size(i,j,k) * &
             (((IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k)) + &
               (IST%flux_sw_nir_dir_top(i,j,k) + IST%flux_sw_nir_dif_top(i,j,k))) * &
               IST%sw_abs_ocn(i,j,k))
      IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + IST%part_size(i,j,k) * &
              ((h2o_ice_to_ocn-h2o_ocn_to_ice)*Idt_slow)

      bsnk(i,j) = bsnk(i,j) - IST%part_size(i,j,k)*bablt ! bot. melt. ablation

    endif

    !
    ! absorb frazil in thinest ice partition available
    !
    if (IST%frazil(i,j)>0 .and. IST%part_size(i,j,0)+IST%part_size(i,j,k)>0.01) then
      !                                                           was ...>0.0
      ! raised above threshold from 0 to 0.01 to avert ocean-ice model blow-ups

      T_Freeze_surf = T_Freeze(IST%s_surf(i,j),IST%ITV)

      if (IST%column_check) then
        enth_prev(i,j,k) = 0.0 ; heat_in(i,j,k) = 0.0

        if (IST%mH_ice(i,j,k) > 0.0) then
          enth_prev(i,j,k) = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
          do m=1,G%NkIce
            enth_prev(i,j,k) = enth_prev(i,j,k) + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
          enddo
          enth_prev(i,j,k) = enth_prev(i,j,k) * IST%part_size(i,j,k)
        endif

      endif

      IST%mH_snow(i,j,k) = IST%mH_snow(i,j,k) * IST%part_size(i,j,k)
      IST%mH_ice(i,j,k)  = IST%mH_ice(i,j,k)  * IST%part_size(i,j,k)
      IST%t_surf(i,j,k) = (IST%t_surf(i,j,k) * IST%part_size(i,j,k) &
                           +(T_0degC + T_Freeze_surf)*IST%part_size(i,j,0))
      IST%part_size(i,j,k) = IST%part_size(i,j,k) + IST%part_size(i,j,0)
      IST%part_size(i,j,0) = 0.0
      I_part = 1.0 / IST%part_size(i,j,k)
      IST%mH_snow(i,j,k) = IST%mH_snow(i,j,k) * I_part
      IST%mH_ice(i,j,k)  = IST%mH_ice(i,j,k)  * I_part
      IST%t_surf(i,j,k) = IST%t_surf(i,j,k) * I_part

      enthalpy(0) = IST%enth_snow(i,j,k,1)
      do m=1,NkIce ; enthalpy(m) = IST%enth_ice(i,j,k,m) ; enddo
      enthalpy_ocean = enthalpy_liquid(IST%t_ocn(i,j), IST%s_surf(i,j), IST%ITV)
      enthalpy(NkIce+1) = enthalpy_ocean

      if (IST%ice_rel_salin > 0.0) then
        do m=1,NkIce ; Salin(m) = IST%sal_ice(i,j,k,m) ; enddo
        salin(NkIce+1) = IST%ice_rel_salin * IST%s_surf(i,j)
      else
        do m=1,NkIce+1 ; Salin(m) = IST%ice_bulk_salin ; enddo
      endif

      call ice_resize_SIS2(IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), G%H_to_kg_m2, &
                   enthalpy, S_col, Salin, 0.0, IST%frazil(i,j) * I_part, &
                   0.0, 0.0, 0.0, T_Freeze_surf, NkIce, &
                   heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                   sn2ic, salt_to_ice, IST%ITV, IST%ice_thm_CSp, Enthalpy_freeze=enth_ocn_to_ice)
      if (IST%mH_ice(i,j,k) == 0.0) then
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy_liquid_freeze(S_col(m), IST%ITV) ; enddo
      else
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy(m) ; enddo
      endif
      if (IST%ice_rel_salin > 0.0) then
        if (IST%mH_ice(i,j,k) == 0.0) then
          do m=1,NkIce ; IST%sal_ice(i,j,k,m) = 0.0 ; enddo
        else
          do m=1,NkIce ; IST%sal_ice(i,j,k,m) = Salin(m) ; enddo
        endif
        salt_change(i,j) = salt_change(i,j) + IST%part_size(i,j,k) * salt_to_ice
      endif

      ! The snow temperature should not have changed.  This should do nothing.
      ! IST%enth_snow(i,j,k,1) = Enthalpy(0)

!      IST%Enth_Mass_in_ocn(i,j) = IST%Enth_Mass_in_ocn(i,j) + &
!          IST%part_size(i,j,k) * enth_ocn_to_ice
      IST%Enth_Mass_in_ocn(i,j) = IST%Enth_Mass_in_ocn(i,j) + &
          IST%part_size(i,j,k) * (h2o_ocn_to_ice * enthalpy_ocean)

      if (IST%column_check) then
        heat_in(i,j,k) = heat_in(i,j,k) - IST%frazil(i,j) * I_part

        enth_here = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,G%NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
        enth_here = enth_here * IST%part_size(i,j,k)
        tot_heat_in = (enth_units * heat_in(i,j,k) + enth_ocn_to_ice) * IST%part_size(i,j,k)
        enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        if (abs(enth_imb) > IST%imb_tol * &
            (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
          enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        endif

        enth_mass_in_col(i,j) = enth_mass_in_col(i,j) + IST%part_size(i,j,k) * enth_ocn_to_ice
      endif

      IST%frazil(i,j) = 0.0
      !
      ! spread frazil salinification of the ocean over all partitions
      !
      IST%lprec_top(i,j,:) = IST%lprec_top(i,j,:) + (h2o_ice_to_ocn-h2o_ocn_to_ice)* &
             IST%part_size(i,j,k) * Idt_slow
      IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + &
             ((h2o_ice_to_ocn-h2o_ocn_to_ice) * IST%part_size(i,j,k)) * Idt_slow
    endif

  enddo ; enddo ; enddo   ! i-, j-, and k-loops
  call mpp_clock_end(iceClock5)

  !
  ! Calculate QFLUX's from (1) restoring to obs and (2) limiting total ice.
  !
  call mpp_clock_begin(iceClock6)
  ! get observed ice thickness for ice restoring, if calculating qflux
  if (IST%do_ice_restore) &
    call get_sea_surface(IST%Time, dum1, Obs_cn_ice, Obs_h_ice)

  if (IST%do_ice_restore .or. IST%do_ice_limit) then
    qflx_lim_ice(:,:) = 0.0
    qflx_res_ice(:,:) = 0.0

    do j=jsc,jec ; do i=isc,iec
      heat_res_ice   = 0.0
      heat_limit_ice = 0.0
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j),IST%ITV)
      !
      ! calculate enthalpy
      !
      do k=1,ncat
        if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)>0.0) then
          e2m(k) = (IST%part_size(i,j,k)*IST%mH_snow(i,j,k)) * G%H_to_kg_m2 * &
                       ((enthalpy_liquid_freeze(0.0, IST%ITV) - &
                         IST%enth_snow(i,j,k,1)) * I_enth_units)
          if (spec_thermo_sal) then ; do m=1,NkIce
            e2m(k) = e2m(k) + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*kg_H_Nk) * &
                       ((enthalpy_liquid_freeze(S_col(m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
          enddo ; else ; do m=1,NkIce
            e2m(k) = e2m(k) + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*kg_H_Nk) * &
                       ((enthalpy_liquid_freeze(IST%sal_ice(i,j,k,m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
          enddo ; endif
        else
          e2m(k) = 0.0
        endif
      enddo
      !
      ! calculate heat needed to constrain ice enthalpy
      !
      if (IST%do_ice_restore) then
        ! Restore to observed enthalpy, implying restoring toward
        ! thickness * concentration.
        heat_res_ice = -(LI*IST%Rho_ice*Obs_h_ice(i,j)*Obs_cn_ice(i,j,2)-sum(e2m(:))) &
                       *dt_slow/(86400*IST%ice_restore_timescale)
      endif

      if (IST%do_ice_limit .and. (sum(e2m) > IST%max_ice_limit*IST%Rho_ice*LI)) then
        heat_limit_ice = sum(e2m(:)) - LI*IST%Rho_ice*IST%max_ice_limit
        ! should we "heat_ice_res = 0.0" ?
      endif

      !
      ! apply constraining heat to ice
      !
      tot_heat = heat_res_ice+heat_limit_ice

      if (tot_heat>0.0) then  ! add like ocean-ice heat
        do k=0,ncat-1
          if (e2m(k) > 0.0) then
            heating = tot_heat/sum(IST%part_size(i,j,k:ncat))
            if (heating*IST%part_size(i,j,k) > e2m(k)) then ! cat. melts away
              IST%mH_ice(i,j,k) = 0.0 ; IST%mH_snow(i,j,k) = 0.0
              tot_heat = tot_heat - e2m(k)
            else
              evap_from_ocn = 0.0; h2o_ice_to_ocn = 0.0; heat_to_ocn = 0.0

              enthalpy(0) = IST%enth_snow(i,j,k,1)
              do m=1,NkIce ; enthalpy(m) = IST%enth_ice(i,j,k,m) ; enddo
              enthalpy(NkIce+1) = enthalpy_liquid(IST%t_ocn(i,j), IST%s_surf(i,j), IST%ITV)

              if (IST%ice_rel_salin > 0.0) then
                do m=1,NkIce ; Salin(m) = IST%sal_ice(i,j,k,m) ; enddo
                salin(NkIce+1) = IST%ice_rel_salin * IST%s_surf(i,j)
              else
                do m=1,NkIce+1 ; Salin(m) = IST%ice_bulk_salin ; enddo
              endif

              call ice_resize_SIS2(IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), &
                          G%H_to_kg_m2, enthalpy, &
                          S_col, Salin, 0.0, 0.0, 0.0, 0.0, &
                          heating, T_Freeze_surf, NkIce, &
                          heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                          snow_to_ice(i,j,k), salt_to_ice, IST%ITV, IST%ice_thm_CSp, bablt )
              if (IST%mH_ice(i,j,k) == 0.0) then
                do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy_liquid_freeze(S_col(m), IST%ITV) ; enddo
              else
                do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy(m) ; enddo
              endif
              if (IST%ice_rel_salin > 0.0) then
                if (IST%mH_ice(i,j,k) == 0.0) then
                  do m=1,NkIce ; IST%sal_ice(i,j,k,m) = 0.0 ; enddo
                else
                  do m=1,NkIce ; IST%sal_ice(i,j,k,m) = Salin(m) ; enddo
                endif
                salt_change(i,j) = salt_change(i,j) + IST%part_size(i,j,k) * salt_to_ice
              endif

              ! The snow temperature should not have changed.  This should do nothing.
              ! IST%enth_snow(i,j,k,1) = Enthalpy(0)
              tot_heat = tot_heat - heating*IST%part_size(i,j,k)
            endif
          endif
        enddo
      endif

      tot_heat = heat_res_ice+heat_limit_ice
      if (tot_heat<0.0) then ! add like frazil
        do k=1,ncat
          if (IST%part_size(i,j,0)+IST%part_size(i,j,k)>0) exit
        enddo
        ! k is thinnest ice partition that can recieve frazil
        IST%mH_snow(i,j,k) = IST%mH_snow(i,j,k) * IST%part_size(i,j,k)
        IST%mH_ice(i,j,k)  = IST%mH_ice(i,j,k)  * IST%part_size(i,j,k)
        IST%t_surf(i,j,k) = IST%t_surf(i,j,k) * IST%part_size(i,j,k) &
                       + (T_0degC + T_Freeze_surf)*IST%part_size(i,j,0)
        IST%part_size(i,j,k) = IST%part_size(i,j,k) + IST%part_size(i,j,0)
        IST%part_size(i,j,0) = 0.0
        IST%mH_snow(i,j,k) = IST%mH_snow(i,j,k) / IST%part_size(i,j,k)
        IST%mH_ice(i,j,k) =  IST%mH_ice(i,j,k)  / IST%part_size(i,j,k)
        IST%t_surf(i,j,k) = IST%t_surf(i,j,k) / IST%part_size(i,j,k)

        enthalpy(0) = IST%enth_snow(i,j,k,1)
        do m=1,NkIce ; enthalpy(m) = IST%enth_ice(i,j,k,m) ; enddo
        enthalpy(NkIce+1) = enthalpy_liquid(IST%t_ocn(i,j), IST%s_surf(i,j), IST%ITV)
        if (IST%ice_rel_salin > 0.0) then
          do m=1,NkIce ; Salin(m) = IST%sal_ice(i,j,k,m) ; enddo
          salin(NkIce+1) = IST%ice_rel_salin * IST%s_surf(i,j)
        else
          do m=1,NkIce+1 ; Salin(m) = IST%ice_bulk_salin ; enddo
        endif

        call ice_resize_SIS2(IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), G%H_to_kg_m2, enthalpy, &
                            S_col, Salin, 0.0, -tot_heat/IST%part_size(i,j,k), &
                            0.0, 0.0, 0.0, T_Freeze_surf, NkIce, &
                            heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, &
                            evap_from_ocn, sn2ic, salt_to_ice, IST%ITV, IST%ice_thm_CSp)
        if (IST%mH_ice(i,j,k) == 0.0) then
          do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy_liquid_freeze(S_col(m), IST%ITV) ; enddo
        else
          do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy(m) ; enddo
        endif
        if (IST%ice_rel_salin > 0.0) then
          if (IST%mH_ice(i,j,k) == 0.0) then
            do m=1,NkIce ; IST%sal_ice(i,j,k,m) = 0.0 ; enddo
          else
            do m=1,NkIce ; IST%sal_ice(i,j,k,m) = Salin(m) ; enddo
          endif
          salt_change(i,j) = salt_change(i,j) + IST%part_size(i,j,k) * salt_to_ice
        endif
        ! The snow temperature should not have changed.  This should do nothing.
        ! IST%enth_snow(i,j,k,1) = Enthalpy(0)
      endif

      ! Convert constraining heat from energy (J/m^2) to flux (W/m^2)
      qflx_lim_ice(i,j) = heat_limit_ice * Idt_slow
      qflx_res_ice(i,j) = heat_res_ice * Idt_slow

    enddo ; enddo
  endif ! End of (IST%do_ice_restore .or. IST%do_ice_limit) block
  call mpp_clock_end(iceClock6)

  if (IST%column_check) then
    enth_col(:,:) = 0.0
    ! Add back any frazil that has not been used yet.
    do j=jsc,jec ; do i=isc,iec
      heat_in_col(i,j) = heat_in_col(i,j) + IST%frazil(i,j) + IST%flux_t_ocn_top(i,j)*dt_slow
    enddo ; enddo

    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
      enth_col(i,j) = enth_col(i,j) + &
        (IST%mH_snow(i,j,k)*IST%part_size(i,j,k)) * IST%enth_snow(i,j,k,1)
      do m=1,G%NkIce
        enth_col(i,j) = enth_col(i,j) + &
          (IST%mH_ice(i,j,k)*IST%part_size(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
      enddo
    endif ; enddo ; enddo ; enddo
    do j=jsc,jec ; do i=isc,iec
      enth_here = enth_col(i,j)
      tot_heat_in = enth_units*heat_in_col(i,j) + enth_mass_in_col(i,j)
      emic2 = (IST%Enth_Mass_in_ocn(i,j) + IST%Enth_Mass_in_atm(i,j) + &
              IST%Enth_Mass_out_ocn(i,j) + IST%Enth_Mass_out_atm(i,j))
      tot_heat_in2 = enth_units*heat_in_col(i,j) + emic2

      enth_imb = enth_here - (enth_prev_col(i,j) + tot_heat_in)
      if (abs(enth_imb) > IST%imb_tol * &
          (abs(enth_here) + abs(enth_prev_col(i,j)) + abs(tot_heat_in)) ) then
        norm_enth_imb = enth_imb / (abs(enth_here) + abs(enth_prev_col(i,j)) + abs(tot_heat_in))
        enth_imb = enth_here - (enth_prev_col(i,j) + tot_heat_in)
      endif
      enth_imb2 = enth_here - (enth_prev_col(i,j) + tot_heat_in2)
      if (abs(enth_imb2) > IST%imb_tol * &
          (abs(enth_here) + abs(enth_prev_col(i,j)) + abs(tot_heat_in2)) ) then
        norm_enth_imb = enth_imb2 / (abs(enth_here) + abs(enth_prev_col(i,j)) + abs(tot_heat_in2))
        enth_imb2 = enth_here - (enth_prev_col(i,j) + tot_heat_in2)
      endif
    enddo ; enddo
  endif

  ! Determine the salt fluxes to ocean
  ! Note that at this point salt_change and h2o_change are the negative of the masses.
  if (IST%ice_rel_salin <= 0.0) then
    do m=1,G%NkIce ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      salt_change(i,j) = salt_change(i,j) + &
         (IST%sal_ice(i,j,k,m)*(IST%mH_ice(i,j,k)*kg_H_Nk)) * IST%part_size(i,j,k)
    enddo ; enddo ; enddo ; enddo
  endif
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    h2o_change(i,j) = h2o_change(i,j) + IST%part_size(i,j,k) * &
                      G%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    ! Note the conversion here from g m-2 to kg m-2 s-1.
    Ice%flux_salt(i2,j2) = salt_change(i,j) * (0.001*Idt_slow)
  enddo ; enddo


  !   The remainder of this routine deals with any thermodynamics diagnostic
  ! output that has been requested.
  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

  yr_dtslow = (864e2*365*Idt_slow)
  if (IST%id_lsnk>0) then
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = min(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(IST%id_lsnk, tmp2d(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  if (IST%id_lsrc>0) then
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = max(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(IST%id_lsrc, tmp2d(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  if (IST%id_saltf>0) call post_data(IST%id_saltf, Ice%flux_salt, IST%diag, &
                                     mask=G%Lmask2dT(isc:iec,jsc:jec))
  if (IST%id_bsnk>0)  call post_data(IST%id_bsnk, bsnk(isc:iec,jsc:jec)*yr_dtslow, &
                                     IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  if (IST%id_tmelt>0) call post_avg(IST%id_tmelt, IST%tmelt, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_bmelt>0) call post_avg(IST%id_bmelt, IST%bmelt, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_sn2ic>0) call post_avg(IST%id_sn2ic, snow_to_ice, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, mask=G%Lmask2dT)
  if (IST%id_qflim>0) call post_data(IST%id_qflim, qflx_lim_ice, IST%diag, mask=G%Lmask2dT)
  if (IST%id_qfres>0) call post_data(IST%id_qfres, qflx_res_ice, IST%diag, mask=G%Lmask2dT)

  call disable_SIS_averaging(IST%diag)

end subroutine SIS2_thermodynamics

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
  integer :: i, j, k, l, i2, j2, k2, i_off, j_off, n
  integer :: isc, iec, jsc, jec, CatIce, nCat_dflt
  logical :: spec_thermo_sal
  character(len=128) :: restart_file, restart_path
  character(len=40)  :: mod = "ice_model" ! This module's name.
  character(len=8)   :: nstr

  type(param_file_type) :: param_file
  type(ice_state_type),    pointer :: IST => NULL()
  type(sea_ice_grid_type), pointer :: G => NULL()

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
  integer :: idr, id_sal
  logical :: read_aux_restart
  character(len=16)  :: stagger, dflt_stagger

  if (associated(Ice%Ice_state)) then
    call SIS_error(WARNING, "ice_model_init called with an associated "// &
                    "Ice%Ice_state structure. Model is already initialized.")
    return
  endif
  allocate(Ice%Ice_state) ; IST => Ice%Ice_state
  allocate(Ice%G) ; G => Ice%G
  allocate(Ice%Ice_restart)

  ! Open the parameter file.
  call Get_SIS_Input(param_file)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
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
  call get_param(param_file, mod, "INTERSPERSED_ICE_THERMO", IST%interspersed_thermo, &
                 "If true, the sea ice thermodynamic updates are applied \n"//&
                 "after the new velocities are determined, but before the \n"//&
                 "transport occurs.  Otherwise, the ice thermodynamic \n"//&
                 "updates occur at the start of the slow ice update and \n"//&
                 "dynamics and continuity can occur together.\n"//&
                 "The default should be changed to false.", default=.true.)
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
  call get_param(param_file, mod, "OCEAN_FILTER_DT", IST%ocean_filter_dt, &
                 "If >0, is the time-scale for a running mean of the\n"//&
                 "ocean surface seen by sea-ice dynamics and icebergs.", &
                 default=0.)

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

  call check_ice_model_nml(param_file)

  if (IST%specified_ice) IST%slab_ice = .true.

  nCat_dflt = 5
  if (IST%slab_ice)  nCat_dflt = 1 ! open water and ice ... but never in same place

  call set_ice_grid(Ice%G, param_file, Ice%domain, nCat_dflt )

  if (IST%slab_ice) G%CatIce = 1 ! open water and ice ... but never in same place
  ! Initialize G%cat_thick_lim here.  ###This needs to be extended to add more options.
  do k=1,min(G%CatIce+1,size(hlim_dflt(:)))
    G%cat_thick_lim(k) = hlim_dflt(k)
  enddo
  if ((G%CatIce+1 > size(hlim_dflt(:))) .and. (size(hlim_dflt(:)) > 1)) then
    do k=min(G%CatIce+1,size(hlim_dflt(:))) + 1, G%CatIce+1
      G%cat_thick_lim(k) =  2.0*G%cat_thick_lim(k-1) - G%cat_thick_lim(k-2)
    enddo
  endif
  do k=1,G%CatIce+1
    G%mH_cat_bound(k) = G%cat_thick_lim(k) * (IST%Rho_ice*G%kg_m2_to_H)
  enddo

  call set_domain(G%Domain%mpp_domain)
  CatIce = G%CatIce

  ! Allocate and register fields for restarts.
  call ice_data_type_register_restarts(G%Domain%mpp_domain, G%CatIce, &
                         param_file, Ice, Ice%Ice_restart, restart_file)

  call ice_state_register_restarts(G, param_file, IST, Ice%Ice_restart, &
                                   restart_file)

  if (IST%Cgrid_dyn) then
    call ice_C_dyn_register_restarts(G, param_file, IST%ice_C_dyn_CSp, &
                                     Ice%Ice_restart, restart_file)
  else
    call ice_B_dyn_register_restarts(G, param_file, IST%ice_B_dyn_CSp, &
                                     Ice%Ice_restart, restart_file)
  endif
!  call ice_transport_register_restarts(G, param_file, IST%ice_transport_CSp, &
!                                       Ice%Ice_restart, restart_file)

  ! Redefine the computational domain sizes to use the ice model's indexing convention.
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc

  IST%coszen(:,:) = cos(3.14*67.0/180.0) ! NP summer solstice.

  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%mask(i2,j2) = ( G%mask2dT(i,j) > 0.5 )
    Ice%area(i2,j2) = G%areaT(i,j) * G%mask2dT(i,j)
  enddo ; enddo

  Ice%Time           = Time
  IST%Time           = Time
  IST%Time_Init      = Time_Init
  IST%Time_step_fast = Time_step_fast
  IST%Time_step_slow = Time_step_slow

  IST%avg_count = 0
  IST%n_calls = 0

  call SIS_diag_mediator_init(G, param_file, IST%diag, component="SIS")
  call set_SIS_axes_info(G, param_file, IST%diag)

  call SIS2_ice_thm_init(param_file, IST%ice_thm_CSp, IST%ITV)

  !
  ! Read the restart file, if it exists.
  !
  allocate(S_col(G%NkIce)) ; S_col(:) = 0.0
  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_unit, &
                             specified_thermo_salinity=spec_thermo_sal)

  restart_path = 'INPUT/'//trim(restart_file)
  if (file_exist(restart_path)) then
    ! Set values of G%H_to_kg_m2 that will permit its absence from the restart
    ! file to be detected, and its difference from the value in this run to
    ! be corrected for.
    H_to_kg_m2_tmp = G%H_to_kg_m2
    G%H_to_kg_m2 = -1.0

    call restore_state(Ice%Ice_restart)

    ! Approximately initialize state fields that are not present
    ! in SIS1 restart files.

    ! Initialize the ice salinity.
    if (.not.query_initialized(Ice%Ice_restart, 'sal_ice')) then
      allocate(sal_ice_tmp(SZI_(G), SZJ_(G), G%CatIce, G%NkIce)) ; sal_ice_tmp(:,:,:,:) = 0.0
      do n=1,G%NkIce
        write(nstr, '(I4)') n ; nstr = adjustl(nstr)
        id_sal = register_restart_field(Ice%Ice_restart, restart_file, 'sal_ice'//trim(nstr), &
                                     sal_ice_tmp(:,:,:,n), domain=G%domain%mpp_domain, &
                                     mandatory=.false., read_only=.true.)
        call restore_state(Ice%Ice_restart, id_sal)
      enddo

      if (query_initialized(Ice%Ice_restart, 'sal_ice1')) then
        do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%sal_ice(i,j,k,1) = sal_ice_tmp(i,j,k,1)
        enddo ; enddo ; enddo
      else
        IST%sal_ice(:,:,:,1) = IST%ice_bulk_salin
      endif
      do n=2,G%NkIce
        write(nstr, '(I4)') n ; nstr = adjustl(nstr)
        if (query_initialized(Ice%Ice_restart, 'sal_ice'//trim(nstr))) then
          do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
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
      allocate(t_snow_tmp(SZI_(G), SZJ_(G), G%CatIce)) ; t_snow_tmp(:,:,:) = 0.0
      allocate(t_ice_tmp(SZI_(G), SZJ_(G), G%CatIce, G%NkIce)) ; t_ice_tmp(:,:,:,:) = 0.0

      idr = register_restart_field(Ice%Ice_restart, restart_file, 't_snow', t_snow_tmp, &
                                   domain=G%domain%mpp_domain, mandatory=.false., read_only=.true.)
      call restore_state(Ice%Ice_restart, idr)
      do n=1,G%NkIce
        write(nstr, '(I4)') n ; nstr = adjustl(nstr)
        idr = register_restart_field(Ice%Ice_restart, restart_file, 't_ice'//trim(nstr), &
                                     t_ice_tmp(:,:,:,n), domain=G%domain%mpp_domain, &
                                     mandatory=.false., read_only=.true.)
        call restore_state(Ice%Ice_restart, idr)
      enddo
    endif

    ! Initialize the ice enthalpy.
    if (.not.query_initialized(Ice%Ice_restart, 'enth_ice')) then
      if (.not.query_initialized(Ice%Ice_restart, 't_ice1')) then
        call SIS_error(FATAL, "Either t_ice1 or enth_ice must be present in the SIS2 restart file "//restart_path)
      endif

      if (spec_thermo_sal) then
        do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,1) = Enth_from_TS(t_ice_tmp(i,j,k,1), S_col(1), IST%ITV)
        enddo ; enddo ; enddo
      else
        do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,1) = Enth_from_TS(t_ice_tmp(i,j,k,1), IST%sal_ice(i,j,k,1), IST%ITV)
        enddo ; enddo ; enddo
      endif

      do n=2,G%NkIce
        write(nstr, '(I4)') n ; nstr = adjustl(nstr)
        if (.not.query_initialized(Ice%Ice_restart, 't_ice'//trim(nstr))) &
          t_ice_tmp(:,:,:,n) = t_ice_tmp(:,:,:,n-1)

        if (spec_thermo_sal) then
          do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
            IST%enth_ice(i,j,k,n) = Enth_from_TS(t_ice_tmp(i,j,k,n), S_col(n), IST%ITV)
          enddo ; enddo ; enddo
        else
          do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
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
          do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
            t_snow_tmp(i,j,k) = Temp_from_En_S(IST%enth_ice(i,j,k,1), IST%sal_ice(i,j,k,1), IST%ITV)
          enddo ; enddo ; enddo
        endif
      endif
      do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
        IST%enth_snow(i,j,k,1) = Enth_from_TS(t_snow_tmp(i,j,k), 0.0, IST%ITV)
      enddo ; enddo ; enddo
    endif

    if (read_aux_restart) deallocate(t_snow_tmp, t_ice_tmp)

    H_rescale_ice = 1.0 ; H_rescale_snow = 1.0
    if (G%H_to_kg_m2 == -1.0) then
      ! This is an older restart file, and the snow and ice thicknesses are in
      ! m, and not a mass coordinate.
      H_rescale_ice = IST%Rho_ice / H_to_kg_m2_tmp
      H_rescale_snow = IST%Rho_snow / H_to_kg_m2_tmp
    elseif (G%H_to_kg_m2 /= H_to_kg_m2_tmp) then
      H_rescale_ice = G%H_to_kg_m2 / H_to_kg_m2_tmp
      H_rescale_snow = H_rescale_ice
    endif
    G%H_to_kg_m2 = H_to_kg_m2_tmp

    ! Deal with any ice masses or thicknesses over land, and rescale to
    ! account for differences between the current thickness units and whatever
    ! thickness units were in the input restart file.
    do k=1,G%CatIce
      IST%mH_snow(:,:,k) = IST%mH_snow(:,:,k) * H_rescale_snow * G%mask2dT(:,:)
      IST%mH_ice(:,:,k) = IST%mH_ice(:,:,k) * H_rescale_ice * G%mask2dT(:,:)
    enddo

    !--- update the halo values.
    call pass_var(IST%part_size, G%Domain)
    call pass_var(IST%mH_ice, G%Domain, complete=.false.)
    call pass_var(IST%mH_snow, G%Domain, complete=.false.)
    do l=1,G%NkIce
      call pass_var(IST%enth_ice(:,:,:,l), G%Domain, complete=.false.)
    enddo
    call pass_var(IST%enth_snow(:,:,:,1), G%Domain, complete=.true.)

    if (IST%Cgrid_dyn) then
      call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
    else
      call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
    endif

    if (IST%ocean_filter_dt>0.) then
       if (IST%Cgrid_dyn) then
         call pass_vector(IST%u_ocn_filt, IST%v_ocn_filt, G%Domain, stagger=CGRID_NE)
       else
         call pass_vector(IST%u_ocn_filt, IST%v_ocn_filt, G%Domain, stagger=BGRID_NE)
       endif
       call pass_var(IST%sea_lev_filt(:,:), G%Domain, complete=.true.)
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
    do n=1,G%NkIce
      enth_spec_ice = Enth_from_TS(0.0, S_col(n), IST%ITV)
      IST%enth_ice(:,:,:,n) = enth_spec_ice
    enddo

    IST%do_init = .true. ! Some more initialization needs to be done in ice_model.
  endif ! file_exist(restart_path)
  deallocate(S_col)

  do k=0,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
    Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo


  call ice_diagnostics_init(Ice, IST, G, IST%diag, IST%Time)

  call ice_thm_param(IST%slab_ice, k_snow, h_lo_lim)

  if (IST%Cgrid_dyn) then
    call ice_C_dyn_init(IST%Time, G, param_file, IST%diag, IST%ice_C_dyn_CSp, IST%ntrunc)
  else
    call ice_B_dyn_init(IST%Time, G, param_file, IST%diag, IST%ice_B_dyn_CSp)
  endif
  call ice_transport_init(IST%Time, G, param_file, IST%diag, IST%ice_transport_CSp)

  ! Register tracers that will be advected around.
  call register_SIS_tracer(IST%enth_ice, G, G%NkIce, "enth_ice", param_file, &
                           IST%TrReg, snow_tracer=.false., &
                           massless_val=massless_ice_enth*enth_unit)
  call register_SIS_tracer(IST%enth_snow, G, 1, "enth_snow", param_file, &
                           IST%TrReg, snow_tracer=.true., &
                           massless_val=massless_snow_enth*enth_unit)
  if (IST%ice_rel_salin > 0.0) then
    call register_SIS_tracer(IST%sal_ice, G, G%NkIce, "salin_ice", param_file, &
                             IST%TrReg, snow_tracer=.false., &
                             massless_val=massless_ice_salin)
  endif
  if (IST%id_age>0) &
    call register_SIS_tracer(IST%age_ice, G, 1, "age_ice", param_file, &
                             IST%TrReg, snow_tracer=.false.)

  call SIS_sum_output_init(G, param_file, "./", Time_Init, IST%sum_output_CSp, IST%ntrunc)

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
  iceClock8 = mpp_clock_id( '  Ice: slow: transport', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock9 = mpp_clock_id( '  Ice: slow: thermodyn diags', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock3 = mpp_clock_id( 'Ice: update fast', flags=clock_flag_default, grain=CLOCK_ROUTINE )

  ! Initialize icebergs
  if (IST%do_icebergs) then
     if( ASSOCIATED(G%Domain%maskmap)) then
       call icebergs_init(Ice%icebergs, &
       G%Domain%niglobal, G%Domain%njglobal, G%Domain%layout, G%Domain%io_layout, &
       Ice%axes(1:2), G%Domain%X_flags, G%Domain%Y_flags, &
       time_type_to_real(Time_step_slow), Time, G%geoLonBu(isc:iec,jsc:jec), G%geoLatBu(isc:iec,jsc:jec), &
       G%mask2dT(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
       G%dxCv(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), G%dyCu(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
       Ice%area, &
       G%cos_rot(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), G%sin_rot(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
       maskmap=G%Domain%maskmap )
     else
       call icebergs_init(Ice%icebergs, &
       G%Domain%niglobal, G%Domain%njglobal, G%Domain%layout, G%Domain%io_layout, &
       Ice%axes(1:2), G%Domain%X_flags, G%Domain%Y_flags, &
       time_type_to_real(Time_step_slow), Time, G%geoLonBu(isc:iec,jsc:jec), G%geoLatBu(isc:iec,jsc:jec), &
       G%mask2dT(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
       G%dxCv(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), G%dyCu(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), &
       Ice%area, &
       G%cos_rot(G%isc-1:G%iec+1,G%jsc-1:G%jec+1), G%sin_rot(G%isc-1:G%iec+1,G%jsc-1:G%jec+1) )
     endif
  endif

  if (IST%add_diurnal_sw .or. IST%do_sun_angle_for_alb) call astronomy_init

  call nullify_domain()

  ! Do any error checking here.
  if (IST%debug) then
    call ice_grid_chksum(G)
  endif

  call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IST%sum_output_CSp)
  IST%write_ice_stats_time = Time_Init + IST%ice_stats_interval * &
      (1 + (IST%Time - Time_init) / IST%ice_stats_interval)

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

  if (IST%Cgrid_dyn) then
    call ice_C_dyn_end(IST%ice_C_dyn_CSp)
  else
    call ice_B_dyn_end(IST%ice_B_dyn_CSp)
  endif
  call ice_transport_end(IST%ice_transport_CSp)
  call SIS2_ice_thm_end(IST%ice_thm_CSp, IST%ITV)

  call ice_grid_end(Ice%G)
  call dealloc_Ice_arrays(Ice)
  call dealloc_IST_arrays(IST)
  deallocate(Ice%Ice_restart)

  ! End icebergs
  if (IST%do_icebergs) call icebergs_end(Ice%icebergs)
  if (IST%add_diurnal_sw .or. IST%do_sun_angle_for_alb) call astronomy_end

  deallocate(Ice%Ice_state)

end subroutine ice_model_end


!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_aging - step the age of sea ice with time step                           !
!             based on Harder, M. (1997) Roughness, age and drift trajectories !
!             of sea ice in large-scale simulations and their use in model     !
!             verifications, Annals of Glaciology 25, p. 237-240.              !
!           - adding a tracer for the oldest ice per category                  !
!             T. Martin, April/June 2008                                       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_aging(G, mi, age, mi_old, dt)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)), intent(in) :: mi, mi_old
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G),1), intent(inout) :: age
  real, intent(in) :: dt   ! has unit seconds

  integer :: i, j, k
  integer :: isc, iec, jsc, jec
  real    :: dt_day
  real    :: mi_min  ! A zero-age mass per unit area, in units of H (often kg m-2).
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  dt_day = dt / 86400.0
  mi_min = 1.0e-7*G%kg_m2_to_H

  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    ! age the ice by one time step, age has units of days.
    age(i,j,k,1) = age(i,j,k,1) + dt_day

    ! Ice age decreases when new ice is formed (mi>mi_old), but is not changed
    ! by ice melt (growth<0) because melt is assumed to occur uniformly to all
    ! ages of ice.
    if (mi(i,j,k) > mi_old(i,j,k)) &
      age(i,j,k,1) = age(i,j,k,1) * (mi_old(i,j,k) / mi(i,j,k))

    ! Ice age is at least 0.01 time step, and excessively thin is set to age 0.
    if ((mi(i,j,k) <= mi_min) .or. (age(i,j,k,1) < 0.01*dt_day)) age(i,j,k,1) = 0.0

  enddo ; enddo ; enddo

end subroutine ice_aging

end module ice_model_mod
