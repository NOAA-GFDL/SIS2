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
! A SEA ICE MODEL for coupling through the GFDL exchange grid;  this module    !
! manages fluxes, diagnostics, and ice timesteps; sea ice dynamics and         !
! thermodynamics are performed in ice_[dyn|thm].f90 - Mike Winton (Michael.Winton@noaa.gov)!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_model_mod

use SIS_diag_mediator, only : set_SIS_axes_info, SIS_diag_mediator_init
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_error_checking, only : chksum, Bchksum, hchksum, check_redundant_B
use SIS_get_input, only : Get_SIS_input

use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
! use MOM_domains,       only : fill_symmetric_edges
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_file_parser, only : open_param_file, close_param_file

use fms_mod, only: file_exist, clock_flag_default
use fms_io_mod, only : set_domain, nullify_domain, restore_state
use mpp_mod, only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only: CLOCK_COMPONENT, CLOCK_LOOP, CLOCK_ROUTINE

use time_manager_mod, only: time_type, time_type_to_real, get_date, get_time
use time_manager_mod, only: operator(+), operator(-), set_date
use astronomy_mod, only: astronomy_init, astronomy_end
use astronomy_mod, only: universal_time, orbital_time, diurnal_solar, daily_mean_solar
  use coupler_types_mod,only: coupler_3d_bc_type
  use constants_mod,    only: hlv, hlf, Tfreeze, grav, STEFAN
use ocean_albedo_mod, only: compute_ocean_albedo            ! ice sets ocean surface
use ocean_rough_mod,  only: compute_ocean_roughness         ! properties over water

use ice_type_mod, only : ice_data_type, ice_state_type, Ice_restart
use ice_type_mod, only : ice_model_restart, dealloc_ice_arrays, dealloc_IST_arrays
use ice_type_mod, only : ice_data_type_register_restarts, ice_state_register_restarts
use ice_type_mod, only : ice_diagnostics_init, ice_stock_pe, check_ice_model_nml
use ice_type_mod, only : ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
use ice_type_mod, only : ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
use ice_type_mod, only : lnd_ice_bnd_type_chksum, ice_data_type_chksum, ice_print_budget
use ice_type_mod, only : IST_chksum, Ice_public_type_chksum
use ice_utils_mod, only : get_avg, post_avg, ice_line, ice_grid_chksum
use ice_grid_mod, only: sea_ice_grid_type, set_ice_grid, ice_grid_end, cell_area
use ice_shortwave_dEdd, only: shortwave_dEdd0_set_params
use ice_spec_mod, only: get_sea_surface

use ice_thm_mod,      only: ice_optics, ice_thm_param, ice5lay_temp, ice5lay_resize
  use ice_thm_mod,      only: MU_TS, TFI, CI, e_to_melt
use ice_dyn_bgrid, only: ice_B_dynamics, ice_B_dyn_init, ice_B_dyn_register_restarts, ice_B_dyn_end
! use ice_dyn_cgrid, only: ice_C_dynamics, ice_C_dyn_init, ice_C_dyn_register_restarts, ice_C_dyn_end
use ice_transport_mod, only : ice_transport, ice_transport_init, ice_transport_end
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
! zero_top_quantities - zero fluxes to begin summing in ice fast physics.      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine zero_top_quantities ( Ice, IST )
  type(ice_data_type),  intent(inout)  :: Ice
  type(ice_state_type), intent(inout) :: IST

  integer :: n, m

  IST%avg_count = 0

  IST%flux_u_top(:,:,:) = 0.0
  IST%flux_v_top(:,:,:) = 0.0

  IST%flux_t_top(:,:,:)          = 0.0
  IST%flux_q_top(:,:,:)          = 0.0
  IST%flux_lw_top(:,:,:)         = 0.0
  IST%flux_lh_top(:,:,:)         = 0.0
  IST%flux_sw_nir_dir_top(:,:,:) = 0.0
  IST%flux_sw_nir_dif_top(:,:,:) = 0.0
  IST%flux_sw_vis_dir_top(:,:,:) = 0.0
  IST%flux_sw_vis_dif_top(:,:,:) = 0.0
  IST%lprec_top(:,:,:)           = 0.0
  IST%fprec_top(:,:,:)           = 0.0
  do n = 1, Ice%ocean_fluxes_top%num_bcs  !{
    do m = 1, Ice%ocean_fluxes_top%bc(n)%num_fields  !{
      Ice%ocean_fluxes_top%bc(n)%field(m)%values(:,:,:) = 0.0
    enddo  !} m
  enddo  !} n

  IST%lwdn(:,:) = 0.0
  IST%swdn(:,:) = 0.0

end subroutine zero_top_quantities

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

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce

  i_off = LBOUND(Ice%albedo_vis_dir,1) - G%isc
  j_off = LBOUND(Ice%albedo_vis_dir,2) - G%jsc

  if (IST%avg_count == 0) call zero_top_quantities (Ice, IST)

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

  do n = 1, Ice%ocean_fluxes_top%num_bcs  !{
    do m = 1, Ice%ocean_fluxes_top%bc(n)%num_fields  !{
      do k2=1,ncat+1 ; do j2=jsc+j_off,jec+j_off ; do i2=isc+i_off,iec+i_off
        Ice%ocean_fluxes_top%bc(n)%field(m)%values(i2,j2,k2) = &
              Ice%ocean_fluxes_top%bc(n)%field(m)%values(i2,j2,k2) + &
             Atmos_boundary_fluxes%bc(n)%field(m)%values(i2,j2,k2)
      enddo ; enddo ; enddo
    enddo  !} m
  enddo  !} n

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
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  logical :: sent

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = 0 ; j_off = 0
  if (Ice%ocean_fluxes_top%num_bcs>=1) then ; if (Ice%ocean_fluxes_top%bc(1)%num_fields>=1) then
    i_off = LBOUND(Ice%ocean_fluxes_top%bc(1)%field(1)%values,1) - G%isc
    j_off = LBOUND(Ice%ocean_fluxes_top%bc(1)%field(1)%values,2) - G%jsc
  endif ; endif
  !
  ! compute average fluxes
  !
  if (IST%avg_count == 0) call SIS_error(FATAL,'avg_top_quantities: '//&
       'no ocean model fluxes have been averaged')

  divid = 1.0/real(IST%avg_count)

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    u = IST%flux_u_top(i,j,k) * divid
    v = IST%flux_v_top(i,j,k) * divid
    IST%flux_u_top(i,j,k) = u*G%cos_rot(i,j)-v*G%sin_rot(i,j) ! rotate stress from lat/lon
    IST%flux_v_top(i,j,k) = v*G%cos_rot(i,j)+u*G%sin_rot(i,j) ! to ocean coordinates
  enddo ; enddo ; enddo

  ! Put wind stress on u,v points and change sign to +down
  call pass_vector(IST%flux_u_top, IST%flux_v_top, G%Domain, stagger=AGRID)
  sign = 1.0 ; if (IST%atmos_winds) sign = -1.0
  do k=0,ncat ; do J=jsc-1,jec ; do I=isc-1,iec
    if ( G%mask2dBu(i,j) > 0.5 ) then
      IST%flux_u_top_bgrid(I,J,k) = sign*0.25*( &
            (IST%flux_u_top(i+1,j+1,k) + IST%flux_u_top(i,j,k)) + &
            (IST%flux_u_top(i+1,j,k) + IST%flux_u_top(i,j+1,k)) )
      IST%flux_v_top_bgrid(I,J,k) = sign*0.25*( &
            (IST%flux_v_top(i+1,j+1,k) + IST%flux_v_top(i,j,k)) + &
            (IST%flux_v_top(i+1,j,k) + IST%flux_v_top(i,j+1,k)) )
    else
      IST%flux_u_top_bgrid(I,J,k) = 0.0
      IST%flux_v_top_bgrid(I,J,k) = 0.0
    endif
  enddo ; enddo ; enddo

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
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    do n=1,Ice%ocean_fluxes_top%num_bcs ; do m=1,Ice%ocean_fluxes_top%bc(n)%num_fields
      Ice%ocean_fluxes_top%bc(n)%field(m)%values(i2,j2,k2) = &
          Ice%ocean_fluxes_top%bc(n)%field(m)%values(i2,j2,k2) * divid
    enddo ; enddo
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
! ice_top_to_ice_bottom - Translate quantities from the ice model's internal   !
!   state to the public ice data type for use by the ocean model.              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_top_to_ice_bottom (Ice, IST, part_size, part_size_uv, G)
  type(ice_data_type),  intent(inout) :: Ice
  type(ice_state_type), intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension (SZI_(G),SZJ_(G),0:G%CatIce), intent(in) :: part_size
  real, dimension (SZIB_(G),SZJB_(G),0:G%CatIce), intent(in) :: part_size_uv

  integer :: i, j, k, isc, iec, jsc, jec, ncat, m, n, i2, j2, k2, i_off, j_off
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%flux_t,1) - G%isc ; j_off = LBOUND(Ice%flux_t,2) - G%jsc

  if (IST%debug) then
    call IST_chksum("Start ice_top_to_ice_bottom", IST, G)
    call Ice_public_type_chksum("Start ice_top_to_ice_bottom", Ice)
  endif

  Ice%flux_u(:,:) = 0.0 ; Ice%flux_v(:,:) = 0.0
  Ice%flux_t(:,:) = 0.0 ; Ice%flux_q(:,:) = 0.0
  Ice%flux_sw_nir_dir(:,:) = 0.0 ; Ice%flux_sw_nir_dif(:,:) = 0.0
  Ice%flux_sw_vis_dir(:,:) = 0.0 ; Ice%flux_sw_vis_dif(:,:) = 0.0
  Ice%flux_lw(:,:) = 0.0 ; Ice%flux_lh(:,:) = 0.0
  Ice%fprec(:,:) = 0.0 ; Ice%lprec(:,:) = 0.0
  do n=1,Ice%ocean_fluxes%num_bcs ; do m=1,Ice%ocean_fluxes%bc(n)%num_fields
    Ice%ocean_fluxes%bc(n)%field(m)%values(:,:) = 0.0
  enddo ; enddo

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1  ! Use these to correct for indexing differences.
    Ice%flux_u(i2,j2) = Ice%flux_u(i2,j2) + IST%flux_u_top_bgrid(I,J,k) * part_size_uv(I,J,k)
    Ice%flux_v(i2,j2) = Ice%flux_v(i2,j2) + IST%flux_v_top_bgrid(I,J,k) * part_size_uv(I,J,k)
    Ice%flux_t(i2,j2) = Ice%flux_t(i2,j2) + IST%flux_t_top(i,j,k) * part_size(i,j,k)
    Ice%flux_q(i2,j2) = Ice%flux_q(i2,j2) + IST%flux_q_top(i,j,k) * part_size(i,j,k)
    Ice%flux_sw_nir_dir(i2,j2) = Ice%flux_sw_nir_dir(i2,j2) + &
            IST%flux_sw_nir_dir_top(i,j,k) * part_size(i,j,k)
    Ice%flux_sw_nir_dif(i2,j2) = Ice%flux_sw_nir_dif(i2,j2) + &
            IST%flux_sw_nir_dif_top(i,j,k) * part_size(i,j,k)
    Ice%flux_sw_vis_dir(i2,j2) = Ice%flux_sw_vis_dir(i2,j2) + &
            IST%flux_sw_vis_dir_top(i,j,k) * part_size(i,j,k)
    Ice%flux_sw_vis_dif(i2,j2) = Ice%flux_sw_vis_dif(i2,j2) + &
            IST%flux_sw_vis_dif_top(i,j,k) * part_size(i,j,k)
    Ice%flux_lw(i2,j2) = Ice%flux_lw(i2,j2) + IST%flux_lw_top(i,j,k) * part_size(i,j,k)
    Ice%flux_lh(i2,j2) = Ice%flux_lh(i2,j2) + IST%flux_lh_top(i,j,k) * part_size(i,j,k)
    Ice%fprec(i2,j2) = Ice%fprec(i2,j2) + IST%fprec_top(i,j,k) * part_size(i,j,k)
    Ice%lprec(i2,j2) = Ice%lprec(i2,j2) + IST%lprec_top(i,j,k) * part_size(i,j,k)

    do n=1,Ice%ocean_fluxes%num_bcs ; do m=1,Ice%ocean_fluxes%bc(n)%num_fields
      Ice%ocean_fluxes%bc(n)%field(m)%values(i2,j2) = Ice%ocean_fluxes%bc(n)%field(m)%values(i2,j2) + &
            Ice%ocean_fluxes_top%bc(n)%field(m)%values(i2,j2,k2) * part_size(i,j,k)
    enddo ; enddo
  enddo ; enddo ; enddo

  if (IST%debug) then
    call IST_chksum("End ice_top_to_ice_bottom", IST, G)
    call Ice_public_type_chksum("End ice_top_to_ice_bottom", Ice)
  endif

end subroutine ice_top_to_ice_bottom

!
! Coupler interface to provide ocean surface data to atmosphere.
!
subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
  type(ocean_ice_boundary_type), intent(inout) :: Ocean_boundary
  type(ice_data_type),           intent(inout) :: Ice

  call mpp_clock_begin(iceClock)
  call mpp_clock_begin(iceClock1)
  call ice_bottom_to_ice_top(Ice, Ice%Ice_state, Ocean_boundary%t, Ocean_boundary%u, Ocean_boundary%v, &
                             Ocean_boundary%frazil, Ocean_boundary, Ice%G, &
                             Ocean_boundary%s, Ocean_boundary%sea_level )
  call mpp_clock_end(iceClock1)
  call mpp_clock_end(iceClock)

end subroutine update_ice_model_slow_up

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_bottom_to_ice_top - prepare surface state for atmosphere fast physics    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_bottom_to_ice_top(Ice, IST, t_surf_ice_bot, u_surf_ice_bot, v_surf_ice_bot, &
                                 frazil_ice_bot, OIB, G, s_surf_ice_bot, sea_lev_ice_bot )
  type(ice_data_type),                     intent(inout) :: Ice
  type(ice_state_type),                    intent(inout) :: IST
  type(sea_ice_grid_type),                 intent(inout) :: G
  real, dimension(G%isc:G%iec,G%jsc:G%jec),   intent(in) :: t_surf_ice_bot, u_surf_ice_bot
  real, dimension(G%isc:G%iec,G%jsc:G%jec),   intent(in) :: v_surf_ice_bot, frazil_ice_bot
  type(ocean_ice_boundary_type),           intent(inout) :: OIB
  real, dimension(G%isc:G%iec,G%jsc:G%jec),   intent(in) :: s_surf_ice_bot, sea_lev_ice_bot

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: sst, tmp
  real, dimension(SZI_(G),SZJ_(G)) :: u_nonsym, v_nonsym
  real :: u, v
  real :: area_pt
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  logical :: sent
  real, parameter                  :: LI = hlf
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%t_surf,1) - G%isc ; j_off = LBOUND(Ice%t_surf,2) - G%jsc
  !
  ! pass ocean state through ice on first partition
  !
  if (.not. IST%specified_ice) then ! otherwise, already set by update_ice_model_slow
    IST%t_surf(isc:iec,jsc:jec,0) = t_surf_ice_bot(isc:iec,jsc:jec)
  endif

  if (IST%do_init) then
    call get_sea_surface(IST%Time, IST%t_surf(isc:iec,jsc:jec,0), IST%part_size(isc:iec,jsc:jec,0:1), &
                         IST%h_ice(isc:iec,jsc:jec,1) )
    call pass_var(IST%part_size(:,:,0), G%Domain, complete=.false. )
    call pass_var(IST%part_size(:,:,1), G%Domain, complete=.false. )
    call pass_var(IST%h_ice(:,:,1), G%Domain, complete=.true. )
    IST%do_init = .false.
  endif

  if (IST%first_time .and. IST%conservation_check) then
    IST%h2o(:) = 0.0
    IST%salt(:) = 0.0
    IST%heat(:) = 0.0

    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      area_pt = G%areaT(i,j) * G%mask2dT(i,j) * IST%part_size(i,j,k)
      IST%h2o(1) = IST%h2o(1) + (IST%Rho_ice*IST%h_ice(i,j,k) + &
                                 IST%Rho_snow*IST%h_snow(i,j,k)) * area_pt
      IST%salt(1) = IST%salt(1) + (IST%Rho_ice*IST%ice_bulk_salin*IST%h_ice(i,j,k)) * area_pt

      if ((IST%part_size(i,j,k)>0.0) .and. (IST%h_ice(i,j,k)>0.0)) then
        if (IST%slab_ice) then
          IST%heat(1) = IST%heat(1) - area_pt * IST%h_ice(i,j,1)*IST%Rho_ice*LI
        else
          IST%heat(1) = IST%heat(1) - area_pt * &
                    e_to_melt(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k),         &
                    IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4))
        endif
      endif
    enddo ; enddo ; enddo
  endif

  ! Any special first-time initialization must be completed before this point.
  IST%first_time = .false.

  if (IST%debug) then
    call IST_chksum("Start ice_bottom_to_ice_top", IST, G)
    call Ice_public_type_chksum("Start ice_bottom_to_ice_top", Ice)
    call chksum(u_surf_ice_bot(isc:iec,jsc:jec), "Start IB2IT u_surf_ice_bot")
    call chksum(v_surf_ice_bot(isc:iec,jsc:jec), "Start IB2IT v_surf_ice_bot")
  endif

  do j=jsc,jec ; do i=isc,iec
    sst(i,j) = IST%t_surf(i,j,0) - Tfreeze
  enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    IST%s_surf(i,j) = s_surf_ice_bot(i,j)
    IST%frazil(i,j) = frazil_ice_bot(i,j)
    IST%sea_lev(i,j) = sea_lev_ice_bot(i,j)
  enddo ; enddo

  call pass_var(IST%sea_lev, G%Domain)

! Transfer the ocean state for extra tracer fluxes.
  do n=1,OIB%fields%num_bcs  ; do m=1,OIB%fields%bc(n)%num_fields
    Ice%ocean_fields%bc(n)%field(m)%values(:,:,1) = OIB%fields%bc(n)%field(m)%values
  enddo ; enddo

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    IST%tmelt(i,j,k) = 0.0
    IST%bmelt(i,j,k) = 0.0
  enddo ; enddo ; enddo

  call get_avg(IST%h_ice(isc:iec,jsc:jec,:), IST%part_size(isc:iec,jsc:jec,1:), &
               tmp(isc:iec,jsc:jec), wtd=.true.)
  do j=jsc,jec ; do i=isc,iec
    if ( tmp(i,j) > 0.0) then
       IST%bheat(i,j) = IST%kmelt*(sst(i,j) + MU_TS*IST%s_surf(i,j))
    else
       IST%bheat(i,j) = 0.0
    endif
  enddo ; enddo

  Ice%ice_mask(:,:,1) = .false.
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%ice_mask(i2,j2,k2) = (IST%h_ice(i,j,k) > 0.0)
  enddo ; enddo ; enddo

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%h_ice(i,j,k) > 0.0) then
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    call ice_optics(IST%h_snow(i,j,k), IST%h_ice(i,j,k), &
         IST%t_surf(i,j,k)-Tfreeze, -MU_TS*IST%s_surf(i,j), &
         Ice%albedo_vis_dir(i2,j2,k2), Ice%albedo_vis_dif(i2,j2,k2), &
         Ice%albedo_nir_dir(i2,j2,k2), Ice%albedo_nir_dif(i2,j2,k2), &
         IST%sw_abs_sfc(i,j,k), &
         IST%sw_abs_snow(i,j,k), IST%sw_abs_ice(i,j,k,1), &
         IST%sw_abs_ice(i,j,k,2), IST%sw_abs_ice(i,j,k,3), &
         IST%sw_abs_ice(i,j,k,4), IST%sw_abs_ocn(i,j,k) , &
         IST%sw_abs_int(i,j,k),  IST%pen(i,j,k), IST%trn(i,j,k), &
         coszen_in=IST%coszen(i,j))
    !
    !Niki: Is the following correct for diagnostics?
    Ice%albedo(i2,j2,k2)=(Ice%albedo_vis_dir(i2,j2,k2)+Ice%albedo_nir_dir(i2,j2,k2)&
                      +Ice%albedo_vis_dif(i2,j2,k2)+Ice%albedo_nir_dif(i2,j2,k2))/4

  endif ; enddo ; enddo ; enddo

  if (G%symmetric) then  ! This is a place-holder until the Tikal release.
    u_nonsym(:,:) = 0.0 ; v_nonsym(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      u_nonsym(i,j) = u_surf_ice_bot(i,j) ! need under-ice current
      v_nonsym(i,j) = v_surf_ice_bot(i,j) ! for water drag term
    enddo ; enddo
    call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, stagger=BGRID_NE)
    do j=jsc-1,jec ; do i=isc-1,iec
      IST%u_ocn(i,j) = u_nonsym(i,j) ! need under-ice current
      IST%v_ocn(i,j) = v_nonsym(i,j) ! for water drag term
    enddo ; enddo
  else
    do j=jsc,jec ; do i=isc,iec
      IST%u_ocn(i,j) = u_surf_ice_bot(i,j) ! need under-ice current
      IST%v_ocn(i,j) = v_surf_ice_bot(i,j) ! for water drag term
    enddo ; enddo
  endif
  !   This will be used with Tikal and later shared code.
  ! if (G%symmetric) &
  !   call fill_symmetric_edges(IST%u_ocn, IST%v_ocn, G%Domain, stagger=BGRID_NE)

  if (IST%debug) then
    call chksum(u_surf_ice_bot(isc:iec,jsc:jec), "Pre-pass u_surf_ice_bot")
    call chksum(v_surf_ice_bot(isc:iec,jsc:jec), "Pre-pass v_surf_ice_bot")
    call chksum(IST%u_ocn(isc:iec,jsc:jec), "Pre-pass IST%u_ocn(0,0)")
    call chksum(IST%v_ocn(isc:iec,jsc:jec), "Pre-pass IST%v_ocn(0,0)")
  endif

  call pass_vector(IST%u_ocn, IST%v_ocn, G%Domain, stagger=BGRID_NE)
  if (IST%debug) then
    call chksum(IST%u_ocn(isc:iec,jsc:jec), "Post-pass IST%u_ocn(0,0)")
    call chksum(IST%v_ocn(isc:iec,jsc:jec), "Post-pass IST%v_ocn(0,0)")
  endif

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
    call SIS_error(FATAL, "ice_bottom_to_ice_top: Cgrid dynamics are not yet available.")
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
        Ice%u_surf(i2,j2,2) = 0.25*((IST%u_ice(I,J) + IST%u_ice(I-1,J-1)) + &
                                    (IST%u_ice(I,J-1) + IST%u_ice(I-1,J)) )
        Ice%v_surf(i2,j2,2) = 0.25*((IST%v_ice(I,J) + IST%v_ice(I-1,J-1)) + &
                                    (IST%v_ice(I,J-1) + IST%v_ice(I-1,J)) )
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
    call chksum(IST%u_ocn(isc:iec,jsc:jec), "Intermed IST%u_ocn(0,0)")
    call chksum(IST%u_ocn(isc-1:iec-1,jsc:jec), "Intermed IST%u_ocn(-,0)")
    call chksum(IST%u_ocn(isc:iec,jsc-1:jec-1), "Intermed IST%u_ocn(0,-)")
    call chksum(IST%u_ocn(isc-1:iec-1,jsc-1:jec-1), "Intermed IST%u_ocn(-,-)")
    call chksum(IST%v_ocn(isc:iec,jsc:jec), "Intermed IST%v_ocn(0,0)")
    call chksum(IST%v_ocn(isc-1:iec-1,jsc:jec), "Intermed IST%v_ocn(-,0)")
    call chksum(IST%v_ocn(isc:iec,jsc-1:jec-1), "Intermed IST%v_ocn(0,-)")
    call chksum(IST%v_ocn(isc-1:iec-1,jsc-1:jec-1), "Intermed IST%v_ocn(-,-)")
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
  !
  ! Pre-timestep diagnostics
  !
  call enable_SIS_averaging(real(time_type_to_real(IST%Time_step_slow)), IST%Time, IST%diag)
  if (IST%id_alb>0) call post_avg(IST%id_alb, Ice%albedo, &
                     IST%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (IST%id_sst>0) call post_data(IST%id_sst, sst(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  if (IST%id_sss>0) call post_data(IST%id_sss, IST%s_surf, IST%diag, mask=G%Lmask2dT)
  if (IST%id_ssh>0) call post_data(IST%id_ssh, IST%sea_lev, IST%diag, mask=G%Lmask2dT)
  if (IST%id_uo >0) call post_data(IST%id_uo , IST%u_ocn, IST%diag, mask=G%Lmask2dBu)
  if (IST%id_vo >0) call post_data(IST%id_vo , IST%v_ocn, IST%diag, mask=G%Lmask2dBu)
  if (IST%id_bheat>0) call post_data(IST%id_bheat, IST%bheat, IST%diag, mask=G%Lmask2dT)
  call disable_SIS_averaging(IST%diag)

  if (IST%debug) then
    call IST_chksum("End ice_bottom_to_ice_top", IST, G)
    call Ice_public_type_chksum("End ice_bottom_to_ice_top", Ice)
  endif

end subroutine ice_bottom_to_ice_top

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
    dhdt, dedt, drdt, & ! d(flux)/d(surf_temp) (+ up)
    albedo_vis_dir, albedo_vis_dif, &
    albedo_nir_dir, albedo_nir_dif, &
    sw_abs_sfc,sw_abs_snow, &
    sw_abs_ocn, sw_abs_int, pen, trn
  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:G%CatIce) :: &
    sw_abs_ice1, sw_abs_ice2,sw_abs_ice3,sw_abs_ice4
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: &
    diurnal_factor, cosz_alb
  real :: dt_fast, ts_new, dts, hf, hfd, latent
  real :: gmt, time_since_ae, cosz, rrsun, fracday, fracday_dt_ice, fracday_day
  real :: rad, cosz_day, cosz_dt_ice, rrsun_day, rrsun_dt_ice
  real :: flux_sw ! sum over dir/dif vis/nir components
  type(time_type) :: Dt_ice
  logical :: sent
  integer :: i, j, k, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Atmos_boundary%t_flux,1) - G%isc
  j_off = LBOUND(Atmos_boundary%t_flux,2) - G%jsc

  rad = acos(-1.)/180.

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
    flux_lh(i,j,k) = hlv*Atmos_boundary%q_flux(i2,j2,k2)
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

  !
  ! implicit update of ice surface temperature
  !
  dt_fast = time_type_to_real(IST%Time_step_fast)

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if (IST%h_ice(i,j,k) > 0.0) then
      if (IST%slab_ice) then
        flux_lh(i,j,k) = flux_lh(i,j,k) + hlf*flux_q(i,j,k)
        latent             = hlv+hlf
      elseif (IST%h_snow(i,j,k)>0.0) then
        flux_lh(i,j,k) = flux_lh(i,j,k) + (hlf-CI*IST%t_snow(i,j,k))*flux_q(i,j,k)
        latent             = hlv + (hlf-CI*IST%t_snow(i,j,k))
      else
        flux_lh(i,j,k) = flux_lh(i,j,k) + hlf*flux_q(i,j,k)*(1-TFI/IST%t_ice(i,j,k,1))
        latent             = hlv + hlf*(1-TFI/IST%t_ice(i,j,k,1))
      endif
      flux_sw = (flux_sw_vis_dir(i,j,k) + flux_sw_vis_dif(i,j,k)) + &
                (flux_sw_nir_dir(i,j,k) + flux_sw_nir_dif(i,j,k))
      !### ADD PARENTHESES.
      hfd = dhdt(i,j,k) + dedt(i,j,k)*latent + drdt(i,j,k)
      hf  = flux_t(i,j,k) + flux_q(i,j,k)*latent - flux_lw(i,j,k)   &
            - (1-IST%pen(i,j,k))*flux_sw - hfd*(IST%t_surf(i,j,k)-Tfreeze)
      call ice5lay_temp(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k),    &
                        IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3),   &
                        IST%t_ice(i,j,k,4), ts_new, hf, hfd,                        &
                        IST%sw_abs_snow(i,j,k)*flux_sw, &
                        IST%sw_abs_ice(i,j,k,1)*flux_sw, &
                        IST%sw_abs_ice(i,j,k,2)*flux_sw, &
                        IST%sw_abs_ice(i,j,k,3)*flux_sw, &
                        IST%sw_abs_ice(i,j,k,4)*flux_sw, &
                        -MU_TS*IST%s_surf(i,j), IST%bheat(i,j), dt_fast, &
                        IST%tmelt(i,j,k), IST%bmelt(i,j,k))
      dts                = ts_new - (IST%t_surf(i,j,k)-Tfreeze)
      IST%t_surf(i,j,k)  = IST%t_surf(i,j,k) + dts
      flux_t(i,j,k)  = flux_t(i,j,k)  + dts * dhdt(i,j,k)
      flux_q(i,j,k)  = flux_q(i,j,k)  + dts * dedt(i,j,k)
      flux_lh(i,j,k) = flux_lh(i,j,k) + dts * dedt(i,j,k) * latent
      flux_lw(i,j,k) = flux_lw(i,j,k) - dts * drdt(i,j,k)
    endif
  enddo ; enddo ; enddo

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

  if (IST%id_sw_abs_snow>0) call post_avg(IST%id_sw_abs_snow, IST%sw_abs_snow, &
                                   IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_abs_ice1>0) call post_avg(IST%id_sw_abs_ice1, IST%sw_abs_ice(:,:,:,1), &
                                   IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_abs_ice2>0) call post_avg(IST%id_sw_abs_ice4, IST%sw_abs_ice(:,:,:,2), &
                                   IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_abs_ice3>0) call post_avg(IST%id_sw_abs_ice4, IST%sw_abs_ice(:,:,:,3), &
                                   IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_abs_ice4>0) call post_avg(IST%id_sw_abs_ice4, IST%sw_abs_ice(:,:,:,4), &
                                   IST%part_size, IST%diag, G=G, mask=G%Lmask2dT)

  if (IST%id_sw_pen>0) call post_avg(IST%id_sw_pen, IST%pen, IST%part_size, &
                                     IST%diag, G=G, mask=G%Lmask2dT)
  if (IST%id_sw_trn>0) call post_avg(IST%id_sw_trn, IST%trn, IST%part_size, &
                                     IST%diag, G=G, mask=G%Lmask2dT)

  if (IST%id_coszen>0) call post_data(IST%id_coszen, IST%coszen, IST%diag, mask=G%Lmask2dT)
  call disable_SIS_averaging(IST%diag)

  if (IST%debug) then
    call IST_chksum("End do_update_ice_model_fast", IST, G)
    call Ice_public_type_chksum("End do_update_ice_model_fast", Ice)
  endif

end subroutine do_update_ice_model_fast

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! update_ice_model_slow - do ice dynamics, transport, and mass changes         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine update_ice_model_slow(Ice, IST, G, runoff, calving, &
                                 runoff_hflx, calving_hflx)

  type(ice_data_type),                intent(inout) :: Ice
  type(ice_state_type),               intent(inout) :: IST
  type(sea_ice_grid_type),            intent(inout) :: G
    real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in), optional :: runoff, calving
    real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in), optional :: runoff_hflx, calving_hflx

  real, dimension(SZIB_(G),SZJB_(G))      :: fx_wat, fy_wat
    real, dimension(G%isc:G%iec,G%jsc:G%jec)      :: hi_change, h2o_change, bsnk, tmp2d, mass
  real, dimension(SZI_(G),SZJ_(G),1:G%CatIce) :: snow_to_ice
  real, dimension(SZI_(G),SZJ_(G),0:G%CatIce) :: &
    part_save
  real, dimension(SZIB_(G),SZJB_(G),0:G%CatIce) :: &
    part_size_uv
!  real, dimension(SZIB_(G),SZJ_(G),0:G%CatIce) :: &
!    part_size_u
!  real, dimension(SZI_(G),SZJB_(G),0:G%CatIce) :: &
!    part_size_v
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: dum1, Obs_h_ice ! for qflux calculation
  real, dimension(G%isc:G%iec,G%jsc:G%jec,2) :: Obs_cn_ice      ! partition 2 = ice concentration
  real, dimension(SZI_(G),SZJ_(G))   :: hs_avg, hi_avg
  real, dimension(SZI_(G),SZJ_(G))   :: tmp1
  real, dimension(SZI_(G),SZJ_(G))   :: qflx_lim_ice, qflx_res_ice
  real, dimension(SZIB_(G),SZJB_(G)) :: wind_stress_x, wind_stress_y
  real, dimension(1:G%CatIce)       :: e2m
  real, dimension(SZIB_(G),SZJ_(G)) :: uc ! Ice velocities interpolated onto
  real, dimension(SZI_(G),SZJB_(G)) :: vc ! a C-grid, in m s-1.
  real, dimension(SZI_(G),SZJ_(G))  :: diagVar ! An temporary array for diagnostics.

  real :: area_h, area_pt
  real :: Idt_slow, yr_dtslow
  integer :: i, j, k, l, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off
  integer ::iyr, imon, iday, ihr, imin, isec
    real                                  :: dt_slow, heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic, bablt
    real                                  :: heat_limit_ice, heat_res_ice
    real                                  :: tot_heat, heating, tot_frazil
    logical                               :: sent
    real, parameter                       :: LI = hlf

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce
  i_off = LBOUND(Ice%runoff,1) - G%isc ; j_off = LBOUND(Ice%runoff,2) - G%jsc
  dt_slow = time_type_to_real(IST%Time_step_slow) ; Idt_slow = 1.0/dt_slow

  if (IST%debug) then
    call IST_chksum("Start update_ice_model_slow", IST, G)
    call Ice_public_type_chksum("Start update_ice_model_slow", Ice)
  endif

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

  !TOM> assume that open water area is not up to date:
  call mpp_clock_end(iceClock)
  call mpp_clock_end(iceClock2)
  tmp1(:,:) = 1.-max(1.-sum(IST%part_size(:,:,1:ncat),dim=3),0.0)
  call get_avg(IST%h_ice, IST%part_size(:,:,1:), hi_avg, wtd=.true.)
  ! Calve off icebergs and integrate forward iceberg trajectories
  if (IST%do_icebergs) &
    call icebergs_run( Ice%icebergs, IST%Time, &
            Ice%calving(:,:), IST%u_ocn(isc-1:iec+1,jsc-1:jec+1), &
            IST%v_ocn(isc-1:iec+1,jsc-1:jec+1), IST%u_ice(isc-1:iec+1,jsc-1:jec+1), IST%v_ice(isc-1:iec+1,jsc-1:jec+1), &
            Ice%flux_u(:,:), Ice%flux_v(:,:), &
            IST%sea_lev(isc-1:iec+1,jsc-1:jec+1), IST%t_surf(isc:iec,jsc:jec,0),  &
            Ice%calving_hflx(:,:), tmp1, hi_avg)
  call mpp_clock_begin(iceClock2)
  call mpp_clock_begin(iceClock)

  call avg_top_quantities(Ice, IST, G) ! average fluxes from update_ice_model_fast


  do k=0,ncat ; do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    part_save(i,j,k) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  !
  ! conservation checks: top fluxes
  !
  call mpp_clock_begin(iceClock7)
  if (IST%conservation_check) then
    tot_frazil = 0.0
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      area_h = G%areaT(i,j) * G%mask2dT(i,j)
      IST%h2o(2) = IST%h2o(2) + (dt_slow * area_h) * &
           (Ice%runoff(i2,j2) + Ice%calving(i2,j2))
      IST%heat(2) = IST%heat(2) + (dt_slow * area_h) * &
           (-LI * Ice%calving(i2,j2))
      tot_frazil = tot_frazil + area_h*IST%frazil(i,j)
    enddo ; enddo
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      area_pt = G%areaT(i,j) * G%mask2dT(i,j) * IST%part_size(i,j,k)
      IST%h2o(2) = IST%h2o(2) + (dt_slow * area_pt) * &
          ( (IST%lprec_top(i,j,k) + IST%fprec_top(i,j,k)) - IST%flux_q_top(i,j,k) )
      IST%heat(2) = IST%heat(2) + (dt_slow * area_pt) * &
          ( ((IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k)) + &
             (IST%flux_sw_nir_dir_top(i,j,k) + IST%flux_sw_nir_dif_top(i,j,k))) + &
            (IST%flux_lw_top(i,j,k) - IST%flux_t_top(i,j,k)) + &
            (-IST%flux_lh_top(i,j,k) - LI*IST%fprec_top(i,j,k)) )
    enddo ; enddo ; enddo
    ! Runoff and calving do not bring in salt, so salt(2) = 0.0
  endif
  call mpp_clock_end(iceClock7)

  ! Dynamics
  !
  call mpp_clock_begin(iceClock4)

  call get_avg(IST%h_snow, IST%part_size(:,:,1:), hs_avg, wtd=.true.)
  call get_avg(IST%h_ice,  IST%part_size(:,:,1:), hi_avg, wtd=.true.)

  if (IST%Cgrid_dyn) then
    call SIS_error(FATAL, "update_ice_model_slow: Cgrid dynamics are not yet available.")

!    part_size_u(:,:,0) = 1.0 ; part_size_u(:,:,1:) = 0.0
!    part_size_v(:,:,0) = 1.0 ; part_size_v(:,:,1:) = 0.0
!    do k=1,ncat
!      do j=jsc,jec ; do I=isc-1,iec
!        part_size_u(I,j,k) = 0.5*G%mask2dCu(I,j) * &
!                             (IST%part_size(i+1,j,k) + IST%part_size(i,j,k))
!        part_size_u(I,j,0) = part_size_u(I,j,0) - part_size_u(I,j,k)
!      enddo ; enddo
!      do J=jsc-1,jec ; do i=isc-,iec
!        part_size_v(i,J,k) = 0.5*G%mask2dCv(i,J) * &
!                             (IST%part_size(i,j+1,k) + IST%part_size(i,j,k))
!        part_size_v(i,J,0) = part_size_v(i,J,0) - part_size_v(i,J,k)
!      enddo ; enddo
!    enddo
!
!    wind_stress_Cu(:,:) = 0.0 ; wind_stress_Cv(:,:) = 0.0
!    call get_avg(IST%flux_u_top_cgrid(isc-1:iec,jsc:jec,1:), part_size_u(isc-1:iec,jsc:jec,1:), &
!                 wind_stress_x(isc-1:iec,jsc:jec), wtd=.true.)
!    call get_avg(IST%flux_v_top_cgrid(isc:iec,jsc-1:jec,1:), part_size_v(isc:iec,jsc-1:jec,1:), &
!                 wind_stress_Cv(isc:iec,jsc-1:jec), wtd=.true.)
!
!    if (IST%debug) then
!      call IST_chksum("Before ice_C_dynamics", IST, G)
!      call hchksum(IST%part_size(:,:,0), "ps(0) before ice_C_dynamics", G)
!      call hchksum(hs_avg, "hs_avg before ice_C_dynamics", G)
!      call hchksum(hi_avg, "hi_avg before ice_C_dynamics", G)
!      call hchksum(IST%sea_lev, "sea_lev before ice_C_dynamics", G, haloshift=1)
!      call uchksum(IST%u_ocn_C, "u_ocn_C before ice_C_dynamics", G, symmetric=.true.)
!      call vchksum(IST%v_ocn_C, "v_ocn_C before ice_C_dynamics", G, symmetric=.true.)
!      call uchksum(wind_stress_Cu, "wind_stress_Cu before ice_C_dynamics", G, symmetric=.true.)
!      call vchksum(wind_stress_Cv, "wind_stress_Cv before ice_C_dynamics", G, symmetric=.true.)
!      call check_redundant_C("wind_stress before ice_C_dynamics", wind_stress_x, wind_stress_y, G)
!      call check_redundant_B("flux_u/v_top before ice_C_dynamics", IST%flux_u_top_bgrid, IST%flux_v_top_bgrid, G)
!      call check_redundant_C("part_size_uv before ice_C_dynamics", part_size_u, part_size_v, G)
!    endif

!    call mpp_clock_begin(iceClocka)
!    call ice_C_dynamics(1.0-IST%part_size(:,:,0), hs_avg, hi_avg, IST%u_ice_C, IST%v_ice_C, &
!                      IST%u_ocn_C, IST%v_ocn_C, &
!                      wind_stress_Cu, wind_stress_Cv, IST%sea_lev, fx_wat_C, fy_wat_C, &
!                      dt_slow, G, IST%ice_C_dyn_CSp)
!    call mpp_clock_end(iceClocka)

!    if (IST%debug) then
!      call IST_chksum("After ice_dynamics", IST, G)
!    endif

!    call mpp_clock_begin(iceClockb)
!    call pass_vector(IST%u_ice, IST%v_ice, G%Domain, stagger=BGRID_NE)
!    call mpp_clock_end(iceClockb)

!    do k=1,ncat ; do J=jsc-1,jec ; do I=isc-1,iec
!      IST%flux_u_top_bgrid(I,J,k) = 0.5*(fx_wat_C(I,j) + fx_wat_C(I,j+1)) ! stress of ice on ocean
!      IST%flux_v_top_bgrid(I,J,k) = 0.5*(fy_wat_C(i,J) + fy_wat_C(i+1,J)) !
!    enddo ; enddo ; enddo
!    call mpp_clock_end(iceClockc)

  else ! B-grid dynamics.

    part_size_uv(:,:,0) = 1.0 ; part_size_uv(:,:,1:) = 0.0
    do k=1,ncat
      do J=jsc-1,jec ; do I=isc-1,iec
        if (G%mask2dBu(I,J) > 0.5 ) then
          part_size_uv(I,J,k) = 0.25*((IST%part_size(i+1,j+1,k) + IST%part_size(i,j,k)) + &
                                      (IST%part_size(i+1,j,k) + IST%part_size(i,j+1,k)) )
        else
          part_size_uv(I,J,k) = 0.0
        endif
        part_size_uv(I,J,0) = part_size_uv(I,J,0) - part_size_uv(I,J,k)
      enddo ; enddo
    enddo

    wind_stress_x(:,:) = 0.0 ; wind_stress_y(:,:) = 0.0
    call get_avg(IST%flux_u_top_bgrid(isc-1:iec,jsc-1:jec,1:), part_size_uv(isc-1:iec,jsc-1:jec,1:), &
                 wind_stress_x(isc-1:iec,jsc-1:jec), wtd=.true.)
    call get_avg(IST%flux_v_top_bgrid(isc-1:iec,jsc-1:jec,1:), part_size_uv(isc-1:iec,jsc-1:jec,1:), &
                 wind_stress_y(isc-1:iec,jsc-1:jec), wtd=.true.)

    if (IST%debug) then
      call IST_chksum("Before ice_dynamics", IST, G)
      call hchksum(IST%part_size(:,:,0), "ps(0) before ice_dynamics", G)
      call hchksum(hs_avg, "hs_avg before ice_dynamics", G)
      call hchksum(hi_avg, "hi_avg before ice_dynamics", G)
      call hchksum(IST%sea_lev, "sea_lev before ice_dynamics", G, haloshift=1)
      call Bchksum(IST%u_ocn, "u_ocn before ice_dynamics", G, symmetric=.true.)
      call Bchksum(IST%v_ocn, "v_ocn before ice_dynamics", G, symmetric=.true.)
      call Bchksum(wind_stress_x, "wind_stress_x before ice_dynamics", G, symmetric=.true.)
      call Bchksum(wind_stress_y, "wind_stress_y before ice_dynamics", G, symmetric=.true.)
      call check_redundant_B("wind_stress before ice_dynamics",wind_stress_x, wind_stress_y, G)
      call check_redundant_B("flux_u/v_top before ice_dynamics",IST%flux_u_top_bgrid, IST%flux_v_top_bgrid, G)
      call check_redundant_B("part_size_uv before ice_dynamics",part_size_uv, G)
    endif

    call mpp_clock_begin(iceClocka)
    call ice_B_dynamics(1.0-IST%part_size(:,:,0), hs_avg, hi_avg, IST%u_ice, IST%v_ice, &
                      IST%u_ocn, IST%v_ocn, &
                      wind_stress_x, wind_stress_y, IST%sea_lev, fx_wat, fy_wat, &
                      dt_slow, G, IST%ice_B_dyn_CSp)
    call mpp_clock_end(iceClocka)

    if (IST%debug) then
      call IST_chksum("After ice_dynamics", IST, G)
    endif

    call mpp_clock_begin(iceClockb)
    call pass_vector(IST%u_ice, IST%v_ice, G%Domain, stagger=BGRID_NE)
    call mpp_clock_end(iceClockb)

    call mpp_clock_begin(iceClockc)
    !
    ! Dynamics diagnostics
    !
    if (IST%id_fax>0) call post_avg(IST%id_fax, IST%flux_u_top_bgrid, &
                                    part_size_uv, IST%diag, G=G)
    if (IST%id_fay>0) call post_avg(IST%id_fay, IST%flux_v_top_bgrid, &
                                    part_size_uv, IST%diag, G=G)

    do k=1,ncat ; do J=jsc-1,jec ; do I=isc-1,iec
      IST%flux_u_top_bgrid(I,J,k) = fx_wat(I,J)  ! stress of ice on ocean
      IST%flux_v_top_bgrid(I,J,k) = fy_wat(I,J)  !
    enddo ; enddo ; enddo
    call mpp_clock_end(iceClockc)
  endif ! End of B-grid dynamics

  call mpp_clock_end(iceClock4)

  !
  ! Thermodynamics
  !
  call mpp_clock_begin(iceClock5)
  if (IST%id_frazil>0) call post_data(IST%id_frazil, IST%frazil(isc:iec,jsc:jec)*Idt_slow, &
                                    IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  call disable_SIS_averaging(IST%diag)

  snow_to_ice(:,:,:) = 0.0
  bsnk(:,:) = 0.0

  hi_change(:,:) = 0.0
  h2o_change(:,:) = 0.0
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    hi_change(i,j) = hi_change(i,J) - IST%h_ice(i,j,k)*IST%part_size(i,j,k)
    h2o_change(i,j) = h2o_change(i,j) - IST%part_size(i,j,k) * &
                  (IST%Rho_snow*IST%h_snow(i,j,k) + IST%Rho_ice*IST%h_ice(i,j,k))
  enddo ; enddo ; enddo

  ! get observed ice thickness for ice restoring, if calculating qflux
  if (IST%do_ice_restore) &
    call get_sea_surface(IST%Time, dum1, Obs_cn_ice, Obs_h_ice)
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if (G%mask2dT(i,j) > 0 .and.IST%h_ice(i,j,k) > 0) then
      ! reshape the ice based on fluxes

      h2o_from_ocn = 0.0; h2o_to_ocn = 0.0; heat_to_ocn = 0.0
      call ice5lay_resize(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k),&
                          IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2),              &
                          IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4),              &
                          IST%fprec_top(i,j,k) *dt_slow, 0.0,                &
                          IST%flux_q_top(i,j,k)*dt_slow,                     &
                          IST%tmelt (i,j,k), IST%bmelt(i,j,k),               &
                          -MU_TS*IST%s_surf(i,j),                            &
                          heat_to_ocn, h2o_to_ocn, h2o_from_ocn,             &
                          snow_to_ice(i,j,k), bablt                          )

      ! modify above-ice to under-ice fluxes for passing to ocean
      IST%flux_q_top(i,j,k) = h2o_from_ocn*Idt_slow ! no ice, evaporation left
      IST%flux_lh_top(i,j,k) = hlv*IST%flux_q_top(i,j,k)
      IST%flux_lw_top(i,j,k) = 0.0
      IST%flux_t_top(i,j,k) = IST%bheat(i,j) - heat_to_ocn*Idt_slow
      IST%flux_sw_vis_dif_top(i,j,k) = (IST%flux_sw_vis_dir_top(i,j,k)+      &
            IST%flux_sw_vis_dif_top(i,j,k)+IST%flux_sw_nir_dir_top(i,j,k)+   &
            IST%flux_sw_nir_dif_top(i,j,k))*IST%pen(i,j,k)*IST%trn(i,j,k)
      IST%flux_sw_nir_dir_top(i,j,k) = 0.0
      IST%flux_sw_nir_dif_top(i,j,k) = 0.0
      IST%flux_sw_vis_dir_top(i,j,k) = 0.0
      IST%fprec_top(i,j,k) = 0.0
      IST%lprec_top(i,j,k) = IST%lprec_top(i,j,k) + h2o_to_ocn*Idt_slow

      bsnk(i,j) = bsnk(i,j) - IST%part_size(i,j,k)*bablt ! bot. melt. ablation

    endif

    !
    ! absorb frazil in thinest ice partition available
    !
    if (IST%frazil(i,j)>0 .and. IST%part_size(i,j,0)+IST%part_size(i,j,k)>0.01) then
      !                                                           was ...>0.0
      ! raised above threshold from 0 to 0.01 to avert ocean-ice model blow-ups
      !
      IST%h_snow(i,j,k) = IST%h_snow(i,j,k) * IST%part_size(i,j,k)
      IST%h_ice(i,j,k)  = IST%h_ice(i,j,k)  * IST%part_size(i,j,k)
      IST%t_surf(i,j,k) = (IST%t_surf(i,j,k) * IST%part_size(i,j,k) &
                           +(Tfreeze - MU_TS*IST%s_surf(i,j))*IST%part_size(i,j,0))
      IST%part_size(i,j,k) = IST%part_size(i,j,k) + IST%part_size(i,j,0)
      IST%part_size(i,j,0) = 0.0
      IST%h_snow(i,j,k) = IST%h_snow(i,j,k) / IST%part_size(i,j,k)
      IST%h_ice(i,j,k)  = IST%h_ice(i,j,k) / IST%part_size(i,j,k)
      IST%t_surf(i,j,k) = IST%t_surf(i,j,k) / IST%part_size(i,j,k)

      call ice5lay_resize(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k), &
                          IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2),                &
                          IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4), 0.0,           &
                          IST%frazil(i,j) / IST%part_size(i,j,k), 0.0, 0.0, 0.0, &
                          -MU_TS*IST%s_surf(i,j),                              &
                          heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic)
      IST%frazil(i,j) = 0.0;
      !
      ! spread frazil salinification over all partitions
      !
      IST%lprec_top  (i,j,:) = IST%lprec_top(i,j,:) + h2o_to_ocn*IST%part_size(i,j,k)/dt_slow
    endif

  enddo ; enddo ; enddo   ! i-, j-, and k-loops
  call mpp_clock_end(iceClock5)

  !
  ! Calculate QFLUX's from (1) restoring to obs and (2) limiting total ice.
  !
  call mpp_clock_begin(iceClock6)
  if (IST%do_ice_restore .or. IST%do_ice_limit) then
    qflx_lim_ice(:,:) = 0.0
    qflx_res_ice(:,:) = 0.0

    do j=jsc,jec ; do i=isc,iec
      heat_res_ice   = 0.0
      heat_limit_ice = 0.0
      !
      ! calculate enthalpy
      !
      if (IST%slab_ice) then
        e2m(1) = IST%h_ice(i,j,1)*IST%Rho_ice*LI
      else
        do k=1,ncat
          if ((IST%part_size(i,j,k)>0.0.and.IST%h_ice(i,j,k)>0.0)) then
             e2m(k) = e_to_melt(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k), &
                      IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3), &
                      IST%t_ice(i,j,k,4)) * IST%part_size(i,j,k)
          else
             e2m(k) = 0.0
          endif
        enddo
      endif
      !
      ! calculate heat needed to constrain ice enthalpy
      !
      if (IST%do_ice_restore) then
        ! TK Mod: restore to observed enthalpy (no change for slab, but for
        !         sis ice, results in restoring toward thickness * concentration

        ! TK Mod for test 10/18/02
        !   Concentration is 1.0 where there is ice for slab case.
        !   The input field Obs_cn_ice may have values less than 1,
        !     so put in if test...

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
      if (IST%slab_ice) IST%h_ice(i,j,1) = IST%h_ice(i,j,1) - tot_heat/(IST%Rho_ice*LI)

      if (.not. IST%slab_ice .and. (tot_heat>0.0)) then  ! add like ocean-ice heat
        do k=0,ncat-1
          if (e2m(k) > 0.0) then
            heating = tot_heat/sum(IST%part_size(i,j,k:ncat))
            if (heating*IST%part_size(i,j,k) > e2m(k)) then ! cat. melts away
              IST%h_ice (i,j,k) = 0.0
              IST%h_snow(i,j,k) = 0.0
              tot_heat = tot_heat - e2m(k)
            else
              h2o_from_ocn = 0.0; h2o_to_ocn = 0.0; heat_to_ocn = 0.0
              call ice5lay_resize(IST%h_snow(i,j,k), IST%t_snow(i,j,k),  &
                                  IST%h_ice(i,j,k), IST%t_ice(i,j,k,1),   &
                                  IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3),  &
                                  IST%t_ice(i,j,k,4),  0.0, 0.0, 0.0, 0.0,&
                                  heating, -MU_TS*IST%s_surf(i,j),       &
                                  heat_to_ocn, h2o_to_ocn, h2o_from_ocn, &
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
        IST%h_snow(i,j,k) = IST%h_snow(i,j,k) * IST%part_size(i,j,k)
        IST%h_ice(i,j,k)  = IST%h_ice(i,j,k)  * IST%part_size(i,j,k)
        IST%t_surf(i,j,k) = IST%t_surf(i,j,k) * IST%part_size(i,j,k) &
                       + (Tfreeze - MU_TS*IST%s_surf(i,j))*IST%part_size(i,j,0)
        IST%part_size(i,j,k) = IST%part_size(i,j,k) + IST%part_size(i,j,0)
        IST%part_size(i,j,0) = 0.0
        IST%h_snow(i,j,k) = IST%h_snow(i,j,k) / IST%part_size(i,j,k)
        IST%h_ice(i,j,k) =  IST%h_ice(i,j,k)  / IST%part_size(i,j,k)
        IST%t_surf(i,j,k) = IST%t_surf(i,j,k) / IST%part_size(i,j,k)

        call ice5lay_resize(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k), &
                            IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2),        &
                            IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4), 0.0,   &
                            -tot_heat/IST%part_size(i,j,k), 0.0, 0.0, 0.0, &
                            -MU_TS*IST%s_surf(i,j),                        &
                            heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic)
      endif

      ! Convert constraining heat from energy (J/m^2) to flux (W/m^2)
      qflx_lim_ice(i,j) = heat_limit_ice / dt_slow
      qflx_res_ice(i,j) = heat_res_ice / dt_slow
      !
      ! Check for energy conservation
      !
      if (IST%slab_ice) then
        e2m(1) = e2m(1) - IST%h_ice(i,j,1)*IST%Rho_ice*LI
      else
        do k=1,ncat
          if (IST%part_size(i,j,k)>0.0.and.IST%h_ice(i,j,k)>0.0) &
            e2m(k) = e2m(k)-e_to_melt(IST%h_snow(i,j,k), IST%t_snow(i,j,k),  &
                     IST%h_ice(i,j,k), IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2), &
                     IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4)) * IST%part_size(i,j,k)
        enddo
      endif
      ! if (abs(sum(e2m) - heat_res_ice - heat_limit_ice)>IST%Rho_ice*LI*1e-3) &
      !       print *, 'QFLUX conservation error at', i, j, 'heat2ice=',  &
      !             tot_heat, 'melted=', sum(e2m), 'h*part_size=', &
      !             IST%h_ice(i,j,:)*IST%part_size(i,j,:)

    enddo ; enddo
  endif ! End of (IST%do_ice_restore .or. IST%do_ice_limit) block
  call mpp_clock_end(iceClock6)
  !
  ! Salt fluxes to ocean
  !
  call mpp_clock_begin(iceClock8)
  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

  ! Note that at this point hi_change and h2o_change are the negative of the masses.
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    hi_change(i,j) = hi_change(i,J) + IST%h_ice(i,j,k)*IST%part_size(i,j,k)
    h2o_change(i,j) = h2o_change(i,j) + IST%part_size(i,j,k) * &
                (IST%Rho_snow*IST%h_snow(i,j,k) + IST%Rho_ice*IST%h_ice(i,j,k))
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off
    Ice%flux_salt(i2,j2) = IST%ice_bulk_salin*IST%Rho_ice*hi_change(i,j) * Idt_slow
  enddo ; enddo
  if (IST%id_saltf>0) call post_data(IST%id_saltf, Ice%flux_salt, IST%diag, &
                                     mask=G%Lmask2dT(isc:iec,jsc:jec))

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
  !
  ! Ice transport ... all ocean fluxes have been calculated by now
  !
  h2o_change(:,:) = 0.0
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    h2o_change(i,j) = h2o_change(i,j) + IST%part_size(i,j,k) * &
                (IST%Rho_snow*IST%h_snow(i,j,k)+IST%Rho_ice*IST%h_ice(i,j,k))
  enddo ; enddo ; enddo

  if (IST%debug) then
    call IST_chksum("Before ice_transport", IST, G)
  endif

  if (IST%Cgrid_dyn) then
    call SIS_error(FATAL, "update_ice_model_slow: Cgrid dynamics are not yet available.")
!    call ice_transport(IST%part_size, IST%h_ice, IST%h_snow, IST%u_ice_C, IST%v_ice_C, &
!                       IST%t_ice, IST%t_snow, IST%sea_lev, G%H_cat_lim, dt_slow, &
!                       G, IST%ice_transport_CSp)
  else
    ! B-grid transport 
    ! Convert the velocities to C-grid points for transport.
    uc(:,:) = 0.0; vc(:,:) = 0.0
    do j=jsc,jec ; do I=isc-1,iec
      uc(I,j) = 0.5 * ( IST%u_ice(I,J-1) + IST%u_ice(I,J) )
    enddo ; enddo
    do J=jsc-1,jec ; do i = isc,iec
      vc(i,J) = 0.5 * ( IST%v_ice(I-1,J) + IST%v_ice(I,J) )
    enddo ; enddo

    call ice_transport(IST%part_size, IST%h_ice, IST%h_snow, uc, vc, &
                       IST%t_ice, IST%t_snow, IST%sea_lev, G%H_cat_lim, dt_slow, &
                       G, IST%ice_transport_CSp)
  endif

  ! Set appropriate surface quantities in categories with no ice.
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%part_size(i,j,k)<1e-10) &
    IST%t_surf(i,j,k) = Tfreeze - MU_TS*IST%s_surf(i,j)
  enddo ; enddo ; enddo

  if (IST%debug) then
    call IST_chksum("After ice_transport", IST, G)
  endif

  ! Convert thickness and concentration to mass.
  mass(:,:) = 0.0
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    mass(i,j) = mass(i,j) + (IST%Rho_snow*IST%h_snow(i,j,k) + &
                             IST%Rho_ice*IST%h_ice(i,j,k)) * IST%part_size(i,j,k)
  enddo ; enddo ; enddo
  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    Ice%mi(i2,j2) = mass(i,j)
  enddo ; enddo
  if (IST%id_mi>0) call post_data(IST%id_mi, mass(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))

  if (IST%do_icebergs) call icebergs_incr_mass(Ice%icebergs, mass(isc:iec,jsc:jec)) ! Add icebergs mass in kg/m^2

  if (IST%id_mib>0) call post_data(IST%id_mib, mass(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec)) ! Diagnose total mass
  if (IST%id_slp>0) call post_data(IST%id_slp, Ice%p_surf, IST%diag, mask=Ice%mask)

  do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
    if (IST%slp2ocean) then
      Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) - 1e5 ! SLP - 1 std. atmosphere, in Pa.
    else
      Ice%p_surf(i2,j2) = 0.0
    endif
    Ice%p_surf(i2,j2) = Ice%p_surf(i2,j2) + grav*mass(i,j)

    h2o_change(i,j) = mass(i,j) - h2o_change(i,j)
  enddo ; enddo

  if (IST%specified_ice) then                  ! over-write changes with specifications
    call get_sea_surface(IST%Time, IST%t_surf(isc:iec,jsc:jec,0), IST%part_size(isc:iec,jsc:jec,:), &
                         IST%h_ice(isc:iec,jsc:jec,1))
    call pass_var(IST%part_size, G%Domain)
  endif

  IST%part_size(:,:,0) = 1.0

  do k=1,ncat
    IST%part_size(:,:,0) = IST%part_size(:,:,0) - IST%part_size(:,:,k)
  enddo
  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%part_size(i2,j2,k2) = IST%part_size(i,j,k)
  enddo ; enddo ; enddo

  call mpp_clock_end(iceClock8)
  !
  ! Thermodynamics diagnostics
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

  if (IST%id_ext>0) then
    do j=jsc,jec ; do i=isc,iec
      diagVar(i,j) = 0.0 ; if (IST%part_size(i,j,0) < 0.85) diagVar(i,j) = 1.0
    enddo ; enddo
    call post_data(IST%id_ext, diagVar(isc:iec,jsc:jec), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  if (IST%id_hs>0) call post_avg(IST%id_hs, IST%h_snow, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_hi>0) call post_avg(IST%id_hi, IST%h_ice, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_hio>0) call post_avg(IST%id_hio, IST%h_ice, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_ts>0) call post_avg(IST%id_ts, IST%t_surf(:,:,1:), IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, offset=-Tfreeze, wtd=.true.)
  if (IST%id_tsn>0) call post_avg(IST%id_tsn, IST%t_snow, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_t1>0) call post_avg(IST%id_t1, IST%t_ice(:,:,:,1), IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_t2>0) call post_avg(IST%id_t2, IST%t_ice(:,:,:,2), IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_t3>0) call post_avg(IST%id_t3, IST%t_ice(:,:,:,3), IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)
  if (IST%id_t4>0) call post_avg(IST%id_t4, IST%t_ice(:,:,:,4), IST%part_size(:,:,1:), &
                                 IST%diag, G=G, mask=G%Lmask2dT, wtd=.true.)

  if (IST%id_xprt>0) call post_data(IST%id_xprt,  h2o_change(isc:iec,jsc:jec)*864e2*365/dt_slow, &
               IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  if (IST%id_e2m>0) then
    tmp2d(:,:) = 0.0
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      if (IST%h_ice(i,j,k)>0.0) &
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k)*e_to_melt(IST%h_snow(i,j,k), &
                     IST%t_snow(i,j,k), IST%h_ice(i,j,k), IST%t_ice(i,j,k,1),    &
                     IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4)    )
    enddo ; enddo ; enddo
    call post_data(IST%id_e2m,  tmp2d(:,:), IST%diag, mask=G%Lmask2dT(isc:iec,jsc:jec))
  endif
  call disable_SIS_averaging(IST%diag)

  call ice_top_to_ice_bottom(Ice, IST, part_save, part_size_uv, G)
  !
  ! conservation checks:  bottom fluxes and final
  !
  if (IST%conservation_check) then
    IST%heat(3) = IST%heat(3) + tot_frazil
    do j=jsc,jec ; do i=isc,iec ; i2 = i+i_off ; j2 = j+j_off
      area_h = G%areaT(i,j) * G%mask2dT(i,j)
      IST%h2o(3) = IST%h2o(3) + dt_slow*area_h * &
           ( (Ice%runoff(i2,j2) + Ice%calving(i2,j2)) + &
             (Ice%lprec(i2,j2) + Ice%fprec(i2,j2)) - Ice%flux_q(i2,j2) )
      IST%heat(3) = IST%heat(3) + dt_slow*area_h * &
           ( (Ice%flux_sw_vis_dir(i2,j2) + Ice%flux_sw_vis_dif(i2,j2) + &
              Ice%flux_sw_nir_dir(i2,j2) + Ice%flux_sw_nir_dif(i2,j2)) + &
             ((Ice%flux_lw(i2,j2) - Ice%flux_t(i2,j2)) - Ice%flux_lh(i2,j2)) &
           - LI*(Ice%fprec(i2,j2) + Ice%calving(i2,j2)) )
      IST%salt(3)  = IST%salt(3) + dt_slow*area_h * Ice%flux_salt(i2,j2) !Kg salt added since start
    enddo ; enddo

    IST%h2o(4)  = 0.0
    IST%salt(4) = 0.0
    IST%heat(4) = 0.0

    do k=1,ncat ; do j=jsc,jec ;  do i=isc,iec
      area_pt = G%areaT(i,j) * G%mask2dT(i,j) * IST%part_size(i,j,k)
      IST%h2o(4) = IST%h2o(4) + (IST%Rho_ice*IST%h_ice(i,j,k) + &
                                 IST%Rho_snow*IST%h_snow(i,j,k)) * area_pt
      IST%salt(4) = IST%salt(4) + (IST%Rho_ice*IST%ice_bulk_salin*IST%h_ice(i,j,k)) * area_pt

      if ((IST%part_size(i,j,k)>0.0) .and. (IST%h_ice(i,j,k)>0.0)) then
        if (IST%slab_ice) then
          IST%heat(4) = IST%heat(4) - area_pt * IST%h_ice(i,j,1)*IST%Rho_ice*LI
        else
          IST%heat(4) = IST%heat(4) - area_pt * e_to_melt(IST%h_snow(i,j,k), &
                    IST%t_snow(i,j,k), IST%h_ice(i,j,k), IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2),   &
                                                         IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4)    )
        endif
      endif
    enddo ;  enddo ; enddo
  endif

  if (IST%verbose) then
    call get_date(IST%Time, iyr, imon, iday, ihr, imin, isec)
    call get_time(IST%Time-set_date(iyr,1,1,0,0,0),isec,iday)
    call ice_line(iyr, iday+1, isec, IST%part_size(isc:iec,jsc:jec,0), &
                              IST%t_surf(:,:,0)-Tfreeze, G)
  endif

  ! Copy the surface properties to the externally visible structure Ice.
  ! These may use different indexing conventions.  This may not be needed here.
  do k=0,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
    Ice%t_surf(i2,j2,k2) = IST%t_surf(i,j,k)
  enddo ; enddo ; enddo

  ! ### ADD BETTER ERROR HANDLING.
  do j=jsc,jec ; do i=isc,iec
    if (G%Lmask2dT(i,j) .and. (abs(sum(IST%part_size(i,j,:))-1.0)>1e-2)) &
         print *,'IST%part_size=',IST%part_size(i,j,:), 'DOES NOT SUM TO 1 AT', &
         G%geoLonT(i,j), G%geoLatT(i,j)
  enddo ; enddo
  call mpp_clock_end(iceClock9)

  if (IST%debug) then
    call IST_chksum("End UIMS", IST, G)
    call Ice_public_type_chksum("End UIMS", Ice)
  endif

end subroutine update_ice_model_slow


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_init - initializes ice model data, parameters and diagnostics      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_init (Ice, Time_Init, Time, Time_step_fast, Time_step_slow )

  type(ice_data_type), intent(inout) :: Ice
  type(time_type)    , intent(in)    :: Time_Init      ! starting time of model integration
  type(time_type)    , intent(in)    :: Time           ! current time
  type(time_type)    , intent(in)    :: Time_step_fast ! time step for the ice_model_fast
  type(time_type)    , intent(in)    :: Time_step_slow ! time step for the ice_model_slow

! This include declares and sets the variable "version".
#include "version_variable.h"
  real :: hlim_dflt(8) = (/ 0.0, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! lower thickness limits 1...NumCat
  integer :: i, j, k, l, i2, j2, k2, i_off, j_off
  integer :: isc, iec, jsc, jec, CatIce, nCat_dflt
  character(len=128) :: restart_file
  character(len=40)  :: mod = "ice_model" ! This module's name.
  type(param_file_type) :: param_file
  type(ice_state_type),    pointer :: IST => NULL()
  type(sea_ice_grid_type), pointer :: G => NULL()

  ! Parameters that are read in and used to initialize other modules.  If those
  ! other modules had control states, these would be moved to those modules.
  real :: mom_rough_ice  ! momentum same, cd10=(von_k/ln(10/z0))^2, in m.
  real :: heat_rough_ice ! heat roughness length, in m.
    ! Parameters that properly belong exclusively to ice_thm.
  real :: alb_snow       ! snow albedo (less if melting), nondim.
  real :: alb_ice        ! ice albedo (less if melting), nondim.
  real :: k_snow         ! snow conductivity (W/mK)
  real :: pen_ice      ! part unreflected solar penetrates ice, nondim.
  real :: opt_dep_ice  ! ice optical depth, in m-1.
  real :: t_range_melt ! melt albedos scaled in over T range, in deg C.
  real    :: h_lo_lim   ! The min ice thickness for temp. calc, in m.
  logical :: do_deltaEdd  ! If true, a delta-Eddington radiative transfer calculation
                          ! for the shortwave radiation within the sea-ice.
  real :: deltaEdd_R_ice  ! Mysterious delta-Eddington tuning parameters, unknown.
  real :: deltaEdd_R_snow ! Mysterious delta-Eddington tuning parameters, unknown.
  real :: deltaEdd_R_pond ! Mysterious delta-Eddington tuning parameters, unknown.


  if (associated(Ice%Ice_state)) then
    call SIS_error(WARNING, "ice_model_init called with an associated "// &
                    "Ice%Ice_state structure. Model is already initialized.")
    return
  endif
  allocate(Ice%Ice_state)
  IST => Ice%Ice_state
  allocate(Ice%G)
  allocate(Ice%Ice_restart)
  G => Ice%G

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
  call get_param(param_file, mod, "SNOW_ALBEDO", alb_snow, &
                 "The albedo of dry snow atop sea ice.", units="nondim", &
                 default=0.85)
  call get_param(param_file, mod, "ICE_ALBEDO", alb_ice, &
                 "The albedo of dry bare sea ice.", units="nondim", &
                 default=0.5826)
  call get_param(param_file, mod, "ICE_SW_PEN_FRAC", pen_ice, &
                 "The fraction of the unreflected shortwave radiation that \n"//&
                 "penetrates into the ice.", units="Nondimensional", default=0.3)
  call get_param(param_file, mod, "ICE_OPTICAL_DEPTH", opt_dep_ice, &
                 "The optical depth of shortwave radiation in sea ice.", &
                 units="m", default=0.67)
  call get_param(param_file, mod, "ALBEDO_T_MELT_RANGE", t_range_melt, &
                 "The temperature range below freezing over which the \n"//&
                 "albedos are changed by partial melting.", units="degC", &
                 default=1.0)
  call get_param(param_file, mod, "ICE_CONSERVATION_CHECK", IST%conservation_check, &
                 "If true, do additional calculations to check for \n"//&
                 "internal conservation of heat, salt, and water mass in \n"//&
                 "the sea ice model.  This does not change answers, but \n"//&
                 "can increase model run time.", default=.true.)
  call get_param(param_file, mod, "DEBUG", IST%debug, &
                 "If true, write out verbose debugging data.", default=.false.)

  call get_param(param_file, mod, "ICE_SEES_ATMOS_WINDS", IST%atmos_winds, &
                 "If true, the sea ice is being given wind stresses with \n"//&
                 "the atmospheric sign convention, and need to have their \n"//&
                 "sign changed.", default=.true.)
  call get_param(param_file, mod, "ICE_BULK_SALINITY", IST%ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "kg/kg", default=0.004)
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
  call get_param(param_file, mod, "DO_DELTA_EDDINGTON_SW", do_deltaEdd, &
                 "If true, a delta-Eddington radiative transfer calculation \n"//&
                 "for the shortwave radiation within the sea-ice.", default=.true.)
  call get_param(param_file, mod, "ICE_DELTA_EDD_R_ICE", deltaEdd_R_ice, &
                 "A dreadfully documented tuning parameter for the radiative \n"//&
                 "propeties of sea ice with the delta-Eddington radiative \n"//&
                 "transfer calculation.", units="perhaps nondimensional?", default=0.0)
  call get_param(param_file, mod, "ICE_DELTA_EDD_R_SNOW", deltaEdd_R_snow, &
                 "A dreadfully documented tuning parameter for the radiative \n"//&
                 "propeties of snow on sea ice with the delta-Eddington \n"//&
                 "radiative transfer calculation.", &
                 units="perhaps nondimensional?", default=0.0)
  call get_param(param_file, mod, "ICE_DELTA_EDD_R_POND", deltaEdd_R_pond, &
                 "A dreadfully documented tuning parameter for the radiative \n"//&
                 "propeties of meltwater ponds on sea ice with the delta-Eddington \n"//&
                 "radiative transfer calculation.", units="perhaps nondimensional?", &
                 default=0.0)

  call check_ice_model_nml(param_file)

  if (IST%specified_ice) IST%slab_ice = .true.

  nCat_dflt = 5
  if (IST%slab_ice)  nCat_dflt = 1 ! open water and ice ... but never in same place


  call set_ice_grid(Ice%G, param_file, Ice%domain, nCat_dflt )

  if (IST%slab_ice) G%CatIce = 1 ! open water and ice ... but never in same place
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

  ! Eliminate this after the Tikal release interface changes.
  Ice_restart => Ice%Ice_restart

  ! Allocate and register fields for restarts.
  restart_file = 'ice_model.res.nc'
  call ice_data_type_register_restarts(G%Domain%mpp_domain, G%CatIce, &
                         param_file, Ice, Ice%Ice_restart, restart_file)

  call ice_state_register_restarts(G, param_file, IST, Ice%Ice_restart, &
                                   restart_file)

  if (IST%Cgrid_dyn) then
    call SIS_error(FATAL, "ice_model_init: Cgrid dynamics are not yet available.")
    ! call ice_C_dyn_register_restarts(G, param_file, IST%ice_C_dyn_CSp, &
    !                                  Ice%Ice_restart, restart_file)
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

  IST%avg_count      = 0

  !
  ! read restart
  !
  restart_file = 'INPUT/ice_model.res.nc'
  if (file_exist(restart_file)) then
    call restore_state(Ice%Ice_restart)

    !--- update the halo values.
    call pass_var(IST%part_size, G%Domain)
    call pass_var(IST%h_ice, G%Domain, complete=.false.)
    call pass_var(IST%h_snow, G%Domain, complete=.false.)
    do l=1,G%NkIce
      call pass_var(IST%t_ice(:,:,:,l), G%Domain, complete=.false.)
    enddo
    call pass_var(IST%t_snow, G%Domain, complete=.true.)

    if (IST%Cgrid_dyn) then
      call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
    endif
    call pass_vector(IST%u_ice, IST%v_ice, G%Domain, stagger=BGRID_NE)
  else ! no restart implies initialization with no ice
    IST%part_size(:,:,:) = 0.0
    IST%part_size(:,:,0) = 1.0

    Ice%rough_mom(:,:,:)   = mom_rough_ice
    Ice%rough_heat(:,:,:)  = heat_rough_ice
    Ice%rough_moist(:,:,:) = heat_rough_ice
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

  call SIS_diag_mediator_init(G, param_file, IST%diag, component="SIS")
  call set_SIS_axes_info(G, param_file, IST%diag)

  call ice_diagnostics_init(Ice, IST, G, IST%diag, IST%Time)

  if (IST%Cgrid_dyn) then
    call SIS_error(FATAL, "ice_model_init: Cgrid dynamics are not yet available.")
    ! call ice_C_dyn_init(IST%Time, G, param_file, IST%diag, IST%ice_C_dyn_CSp)
  else
    call ice_B_dyn_init(IST%Time, G, param_file, IST%diag, IST%ice_B_dyn_CSp)
  endif
  call ice_transport_init(IST%Time, G, param_file, IST%diag, IST%ice_transport_CSp)
  call ice_thm_param(alb_snow, alb_ice, pen_ice, opt_dep_ice, IST%slab_ice, &
                     t_range_melt, k_snow, h_lo_lim, do_deltaEdd)

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
  if (IST%do_icebergs) call icebergs_init(Ice%icebergs, &
       G%Domain%niglobal, G%Domain%njglobal, G%Domain%layout, G%Domain%io_layout, &
       Ice%axes(1:2), G%Domain%maskmap, G%Domain%X_flags, G%Domain%Y_flags, &
       time_type_to_real(Time_step_slow), Time, G%geoLonBu(isc:iec,jsc:jec), G%geoLatBu(isc:iec,jsc:jec), &
       G%mask2dT, G%dxCv, G%dyCu, Ice%area, G%cos_rot, G%sin_rot )

  if (IST%add_diurnal_sw .or. IST%do_sun_angle_for_alb) call astronomy_init

  if (do_deltaEdd) &
    call shortwave_dEdd0_set_params(deltaEdd_R_ice, deltaEdd_R_snow, deltaEdd_R_pond)

  call nullify_domain()

  ! Do any error checking here.
  if (IST%debug) then
    call ice_grid_chksum(G)
  endif

end subroutine ice_model_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_model_end - writes the restart file and deallocates memory               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_model_end (Ice)
  type(ice_data_type), intent(inout) :: Ice

  type(ice_state_type), pointer :: IST => NULL()

  IST => Ice%Ice_state
  if (IST%conservation_check) call ice_print_budget(IST)

  call ice_model_restart(Ice=Ice)

  !--- release memory ------------------------------------------------

  if (IST%Cgrid_dyn) then
    call SIS_error(FATAL, "ice_model_end: Cgrid dynamics are not yet available.")
    ! call ice_C_dyn_end(IST%ice_C_dyn_CSp)
  else
    call ice_B_dyn_end(IST%ice_B_dyn_CSp)
  endif
  call ice_transport_end(IST%ice_transport_CSp)

  call ice_grid_end(Ice%G)
  call dealloc_Ice_arrays(Ice)
  call dealloc_IST_arrays(IST)
  deallocate(Ice%Ice_restart)

  ! End icebergs
  if (IST%do_icebergs) call icebergs_end(Ice%icebergs)
  if (IST%add_diurnal_sw .or. IST%do_sun_angle_for_alb) call astronomy_end

  deallocate(Ice%Ice_state)

end subroutine ice_model_end

end module ice_model_mod
