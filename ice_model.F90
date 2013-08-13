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

use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg

  use mpp_mod,          only: mpp_clock_begin, mpp_clock_end
  use mpp_domains_mod,  only: mpp_update_domains, BGRID_NE, CGRID_NE
  use fms_mod,          only: error_mesg
  use diag_manager_mod, only: send_data
  use time_manager_mod, only: time_type, operator(+), get_date, get_time
  use time_manager_mod, only: operator(-), set_date
  use astronomy_mod,    only: universal_time, orbital_time, diurnal_solar, daily_mean_solar
  use coupler_types_mod,only: coupler_2d_bc_type, coupler_3d_bc_type
  use constants_mod,    only: hlv, hlf, Tfreeze, grav, STEFAN
  use ocean_albedo_mod, only: compute_ocean_albedo            ! ice sets ocean surface
  use ocean_rough_mod,  only: compute_ocean_roughness         ! properties over water
use ice_type_mod,     only: ice_data_type, ice_state_type
  use ice_type_mod,     only: ice_model_init, ice_model_end, hlim,&
                              mom_rough_ice, heat_rough_ice, atmos_winds, kmelt, &
                              slab_ice, spec_ice, ice_bulk_salin, id_cn, id_hi, id_hio,  &
                              id_hs, id_tsn, id_t1, id_t2, id_t3, id_t4, id_ts,  &
                              id_sh, id_lh, id_sw,                               &
                              id_swdn, id_lw, id_lwdn, id_snofl, id_rain,        &
                              id_runoff_hflx, id_calving_hflx,                   &
                              id_runoff, id_calving, id_evap, id_saltf, id_tmelt,&
                              id_mi, id_bmelt, id_bheat, id_frazil, id_alb,      &
                              id_xprt, id_lsrc, id_lsnk, id_bsnk, id_strna,      &
                              id_fax, id_fay,                                    &
                              id_sn2ic,  id_ext, id_slp, id_sst, id_sss,         &
                              id_ssh, id_uo, id_vo, id_e2m, id_qflim, id_qfres,  &
                              id_ix_trans, id_iy_trans,                          &
                              do_ice_restore, do_ice_limit, max_ice_limit,       &
                              ice_restore_timescale, do_init, h2o, heat, salt,   &
                              conservation_check, slp2ocean, iceClock, verbose,  &
                              iceClock1, iceClock2, iceClock3, &
                              iceClock4, iceClock5, iceClock6, &
                              iceClock7, iceClock8, iceClock9, &
                              iceClocka, iceClockb, iceClockc, &
                              id_sw_vis, id_sw_dir, id_sw_dif, id_sw_vis_dir,    &
                              id_sw_vis_dif, id_sw_nir_dir, id_sw_nir_dif, id_coszen, &
                              ice_stock_pe, do_icebergs, ice_model_restart, &
                              add_diurnal_sw, id_mib, ice_data_type_chksum
  use ice_type_mod,     only: do_sun_angle_for_alb,              &
			      id_alb_vis_dir, id_alb_vis_dif,    &
	                      id_alb_nir_dir, id_alb_nir_dif
  use ice_type_mod,     only: id_sw_abs_snow,id_sw_abs_ice1,id_sw_abs_ice2,id_sw_abs_ice3,id_sw_abs_ice4,&
                              id_sw_pen,id_sw_trn

  use ice_grid_mod,     only: uv_to_t, t_to_uv, cut_check, tripolar_grid
  use ice_grid_mod,     only: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, im, jm, km
  use ice_grid_mod,     only: ice_avg, all_avg, get_avg, ice_line
  use ice_grid_mod,     only: cell_area, sin_rot, cos_rot, latitude
  use ice_grid_mod,     only: sea_ice_grid_type
  use ice_spec_mod,     only: get_sea_surface

  !
  ! the following four modules are the work horses of the sea ice model
  !
  use ice_thm_mod,      only: ice_optics, ice_thm_param, ice5lay_temp, ice5lay_resize
  use ice_thm_mod,      only: DI, DS, MU_TS, TFI, CI, e_to_melt
  use ice_dyn_mod,      only: ice_dynamics
  use ice_transport_mod, only : ice_transport
  use ice_bergs,        only: icebergs_run, icebergs_incr_mass

  implicit none ; private

  public :: ice_data_type, ocean_ice_boundary_type, atmos_ice_boundary_type, land_ice_boundary_type
  public :: ice_model_init, ice_model_end, update_ice_model_fast, ice_stock_pe, cell_area
  public :: update_ice_model_slow_up, update_ice_model_slow_dn
  public :: ice_model_restart  ! for intermediate restart
  public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum, &
            lnd_ice_bnd_type_chksum, ice_data_type_chksum

  !
  ! the following three types are for data exchange with the new coupler
  ! they are defined here but declared in coupler_main and allocated in flux_init
  !
  type :: ocean_ice_boundary_type
     real, dimension(:,:),   pointer :: u         =>NULL()
     real, dimension(:,:),   pointer :: v         =>NULL()
     real, dimension(:,:),   pointer :: t         =>NULL()
     real, dimension(:,:),   pointer :: s         =>NULL()
     real, dimension(:,:),   pointer :: frazil    =>NULL()
     real, dimension(:,:),   pointer :: sea_level =>NULL()
     real, dimension(:,:,:), pointer :: data      =>NULL() ! collective field for "named" fields above
     integer                         :: xtype              ! REGRID, REDIST or DIRECT used by coupler
     type(coupler_2d_bc_type)        :: fields     ! array of fields used for additional tracers
  end type 

  type :: atmos_ice_boundary_type 
     real, dimension(:,:,:), pointer :: u_flux  =>NULL()
     real, dimension(:,:,:), pointer :: v_flux  =>NULL()
     real, dimension(:,:,:), pointer :: u_star  =>NULL()
     real, dimension(:,:,:), pointer :: t_flux  =>NULL()
     real, dimension(:,:,:), pointer :: q_flux  =>NULL()
     real, dimension(:,:,:), pointer :: lw_flux =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_vis_dir =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_vis_dif =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_nir_dir =>NULL()
     real, dimension(:,:,:), pointer :: sw_flux_nir_dif =>NULL()
     real, dimension(:,:,:), pointer :: lprec   =>NULL()
     real, dimension(:,:,:), pointer :: fprec   =>NULL()
     real, dimension(:,:,:), pointer :: dhdt    =>NULL()
     real, dimension(:,:,:), pointer :: dedt    =>NULL()
     real, dimension(:,:,:), pointer :: drdt    =>NULL()
     real, dimension(:,:,:), pointer :: coszen  =>NULL()
     real, dimension(:,:,:), pointer :: p       =>NULL()
     real, dimension(:,:,:), pointer :: data    =>NULL()
     integer                         :: xtype
     type(coupler_3d_bc_type)        :: fluxes     ! array of fluxes used for additional tracers
  end type

  type :: land_ice_boundary_type
     real, dimension(:,:),   pointer :: runoff  =>NULL()
     real, dimension(:,:),   pointer :: calving =>NULL()
     real, dimension(:,:),   pointer :: runoff_hflx  =>NULL()
     real, dimension(:,:),   pointer :: calving_hflx =>NULL()
     real, dimension(:,:,:), pointer :: data    =>NULL() ! collective field for "named" fields above
     integer                         :: xtype            ! REGRID, REDIST or DIRECT used by coupler
  end type

  logical :: first_time = .true.               ! first time ice_bottom_to_ice_top

interface post_avg
  module procedure post_avg_all, post_avg_G
end interface post_avg

contains

!#######################################################################
!
! Coupler interface to do slow ice processes:  dynamics, transport, mass
!
subroutine update_ice_model_slow_dn ( Atmos_boundary, Land_boundary, Ice )
  type(atmos_ice_boundary_type), intent(inout) :: Atmos_boundary
  type(land_ice_boundary_type),  intent(inout) :: Land_boundary
  type (ice_data_type),          intent(inout) :: Ice

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
  type (ice_data_type),             intent(inout) :: Ice
  type(ice_state_type),             intent(inout) :: IST
  type(coupler_3d_bc_type),         intent(inout) :: Atmos_boundary_fluxes
  type(sea_ice_grid_type),          intent(inout) :: G
  real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_u,  flux_v, flux_t, flux_q
  real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_lw, lprec, fprec, flux_lh
  real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_sw_nir_dir
  real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_sw_nir_dif
  real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_sw_vis_dir
  real, dimension(isc:iec,jsc:jec,km), intent(in) :: flux_sw_vis_dif

  real,dimension(isc:iec,jsc:jec)                 :: tmp
  integer                                         :: i, j, k, m, n

  if (IST%avg_count == 0) call zero_top_quantities (Ice, Ice%Ice_state)

  do k=1,km ; do j=jsc,jec ; do i=isc,iec
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
      do k=1,km ; do j=jsc,jec ; do i=isc,iec
        Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) = &
              Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) +     &
             Atmos_boundary_fluxes%bc(n)%field(m)%values(i,j,k)
      enddo ; enddo ; enddo
    enddo  !} m
  enddo  !} n

  if (id_lwdn > 0) then
    call get_avg(flux_lw(:,:,:) + STEFAN*Ice%t_surf(:,:,:)**4, &
                       Ice%part_size(isc:iec,jsc:jec,:), tmp(:,:))
    do j=jsc,jec ; do i=isc,iec
      if (Ice%mask(i,j)) IST%lwdn(i,j) = IST%lwdn(i,j) + tmp(i,j)
    enddo ; enddo
  endif

  if (id_swdn > 0) then
    !### REWRITE TO MERGE LOOPS
    call get_avg(flux_sw_vis_dir(:,:,:)/(1-Ice%albedo_vis_dir(:,:,:)) + &
                       flux_sw_vis_dif(:,:,:)/(1-Ice%albedo_vis_dif(:,:,:)) + &
                       flux_sw_nir_dir(:,:,:)/(1-Ice%albedo_nir_dir(:,:,:)) + &
                       flux_sw_nir_dif(:,:,:)/(1-Ice%albedo_nir_dif(:,:,:)), &
         Ice%part_size(isc:iec,jsc:jec,:), tmp(:,:) )
    do j=jsc,jec ; do i=isc,iec
      if (Ice%mask(i,j)) IST%swdn(i,j) = IST%swdn(i,j) + tmp(i,j)
    enddo ; enddo
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
  integer :: i, j, k, m, n
  logical :: sent
  !
  ! compute average fluxes
  !
  if (IST%avg_count == 0) call error_mesg ('avg_top_quantities', &
       'no ocean model fluxes have been averaged', 3)

  divid = 1.0/real(IST%avg_count)

  do k=1,km ; do j=jsc,jec ; do i=isc,iec
     u = IST%flux_u_top(i,j,k) * divid
     v = IST%flux_v_top(i,j,k) * divid
     IST%flux_u_top(i,j,k) = u*cos_rot(i,j)-v*sin_rot(i,j) ! rotate stress from lat/lon
     IST%flux_v_top(i,j,k) = v*cos_rot(i,j)+u*sin_rot(i,j) ! to ocean coordinates
  enddo ; enddo ; enddo

  ! Put wind stress on u,v points and change sign to +down
  call mpp_update_domains(IST%flux_u_top, IST%flux_v_top, Domain  )
  !### ADD PARENTHESIS FOR REPRODUCIBILITY.
  sign = 1.0 ; if (atmos_winds) sign = -1.0
  do k=1,km ; do j=jsc,jec ; do i=isc,iec
    if ( Ice%G%mask2dBu(i,j) > 0.5 ) then
      IST%flux_u_top_bgrid(i,j,k) = sign*0.25*( &
            IST%flux_u_top(i+1,j+1,k) + IST%flux_u_top(i+1,j,k) + &
            IST%flux_u_top(i,j+1,k) + IST%flux_u_top(i,j,k) )
      IST%flux_v_top_bgrid(i,j,k) = sign*0.25*( &
            IST%flux_v_top(i+1,j+1,k) + IST%flux_v_top(i+1,j,k) + &
            IST%flux_v_top(i,j+1,k) + IST%flux_v_top(i,j,k) )
    else
      IST%flux_u_top_bgrid(i,j,k) = 0.0
      IST%flux_v_top_bgrid(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo

  do k = 1,km ; do j=jsc,jec ; do i=isc,iec
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
    do n = 1, Ice%ocean_fluxes_top%num_bcs  !{
      do m = 1, Ice%ocean_fluxes_top%bc(n)%num_fields  !{
        Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) = Ice%ocean_fluxes_top%bc(n)%field(m)%values(i,j,k) * divid
      enddo  !} m
    enddo  !} n
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    IST%lwdn(i,j) = IST%lwdn(i,j)* divid
    IST%swdn(i,j) = IST%swdn(i,j)* divid
  enddo ; enddo
  !
  ! Flux diagnostics
  !
  if (id_sh>0) call post_avg(id_sh, IST%flux_t_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_lh>0) call post_avg(id_lh, IST%flux_lh_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_evap>0) call post_avg(id_evap, IST%flux_q_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw>0) then
    do j=jsc,jec ; do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo ; enddo
    do k=1,km ; do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = tmp2d(i,j) + Ice%part_size(i,j,k) * ( &
            IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k) + &
            IST%flux_sw_nir_dir_top(i,j,k) + IST%flux_sw_nir_dif_top(i,j,k) )
    enddo ; enddo ; enddo
    sent = send_data(id_sw, tmp2d(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
  endif
  if (id_lw>0) call post_avg(id_lw, IST%flux_lw_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_snofl>0) call post_avg(id_snofl, IST%fprec_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_rain>0) call post_avg(id_rain, IST%lprec_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_lwdn>0) sent = send_data(id_lwdn, IST%lwdn(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
  if (id_swdn>0) sent = send_data(id_swdn, IST%swdn(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
  if (id_sw_vis>0) then
    do j=jsc,jec ; do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo ; enddo
    do k=1,km ; do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = tmp2d(i,j) + Ice%part_size(i,j,k) * ( &
            IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k) )
    enddo ; enddo ; enddo
    sent = send_data(id_sw_vis, tmp2d(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
!    sent = send_data(id_sw_vis, all_avg(IST%flux_sw_vis_dif_top(isc:iec,jsc:jec,:)+&
!                          IST%flux_sw_vis_dir_top(isc:iec,jsc:jec,:), &
!                          Ice%part_size(isc:iec,jsc:jec,:)), Ice%Time, mask=Ice%mask)
  endif
  if (id_sw_nir_dir>0) call post_avg(id_sw_nir_dir, IST%flux_sw_nir_dir_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw_nir_dif>0) call post_avg(id_sw_nir_dif, IST%flux_sw_nir_dif_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw_vis_dir>0) call post_avg(id_sw_vis_dir, IST%flux_sw_vis_dir_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw_vis_dif>0) call post_avg(id_sw_vis_dif, IST%flux_sw_vis_dif_top(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
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
  type (ice_data_type), intent(inout) :: Ice
  type(ice_state_type), intent(inout) :: IST
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension (:,:,:), intent(in) :: part_size, part_size_uv
  integer                             :: m, n

!    Ice%flux_u  = all_avg( IST%flux_u_top_bgrid(isc:iec,jsc:jec,:) , part_size_uv(isc:iec,jsc:jec,:) )
!    Ice%flux_v  = all_avg( IST%flux_v_top_bgrid(isc:iec,jsc:jec,:) , part_size_uv(isc:iec,jsc:jec,:) )

  Ice%flux_u  = all_avg( IST%flux_u_top_bgrid(isc:iec,jsc:jec,:) , part_size_uv )
  Ice%flux_v  = all_avg( IST%flux_v_top_bgrid(isc:iec,jsc:jec,:) , part_size_uv )
  Ice%flux_t  = all_avg( IST%flux_t_top , part_size )
  Ice%flux_q  = all_avg( IST%flux_q_top , part_size )
  Ice%flux_sw_nir_dir = all_avg( IST%flux_sw_nir_dir_top, part_size )
  Ice%flux_sw_nir_dif = all_avg( IST%flux_sw_nir_dif_top, part_size )
  Ice%flux_sw_vis_dir = all_avg( IST%flux_sw_vis_dir_top, part_size )
  Ice%flux_sw_vis_dif = all_avg( IST%flux_sw_vis_dif_top, part_size )
  Ice%flux_lw = all_avg( IST%flux_lw_top, part_size )
  Ice%fprec   = all_avg( IST%fprec_top  , part_size )
  Ice%lprec   = all_avg( IST%lprec_top  , part_size )
  Ice%flux_lh = all_avg( IST%flux_lh_top, part_size )
  do n = 1, Ice%ocean_fluxes%num_bcs  !{
    do m = 1, Ice%ocean_fluxes%bc(n)%num_fields  !{
      Ice%ocean_fluxes%bc(n)%field(m)%values =                &
           all_avg(Ice%ocean_fluxes_top%bc(n)%field(m)%values, part_size)
    enddo  !} m
  enddo  !} n

end subroutine ice_top_to_ice_bottom

!
! Coupler interface to provide ocean surface data to atmosphere.
!
subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
  type(ocean_ice_boundary_type), intent(inout) :: Ocean_boundary
  type (ice_data_type),          intent(inout) :: Ice

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
subroutine ice_bottom_to_ice_top (Ice, IST, t_surf_ice_bot, u_surf_ice_bot, v_surf_ice_bot, &
                                  frazil_ice_bot, Ocean_ice_boundary, G, &
                                  s_surf_ice_bot, sea_lev_ice_bot )
  type(ice_data_type),                     intent(inout) :: Ice
  type(ice_state_type),                    intent(inout) :: IST
  type(sea_ice_grid_type),                 intent(inout) :: G
  real, dimension(isc:iec,jsc:jec),           intent(in) :: t_surf_ice_bot, u_surf_ice_bot
  real, dimension(isc:iec,jsc:jec),           intent(in) :: v_surf_ice_bot, frazil_ice_bot
  type(ocean_ice_boundary_type),           intent(inout) :: Ocean_ice_boundary
  real, dimension(isc:iec,jsc:jec), intent(in), optional :: s_surf_ice_bot, sea_lev_ice_bot

  real, dimension(isc:iec,jsc:jec) :: sst, tmp
  real                             :: u, v
  integer                          :: i, j, k, m, n
  logical                          :: sent
  real, parameter                  :: LI = hlf

  !
  ! pass ocean state through ice on first partition
  !
  if (.not. spec_ice) & ! otherwise, already set by update_ice_model_slow
    Ice%t_surf(:,:,1) = t_surf_ice_bot

  if (do_init) then
    call get_sea_surface(Ice%Time, Ice%t_surf(isc:iec,jsc:jec,1), Ice%part_size(isc:iec,jsc:jec,1:2), &
                         IST%h_ice(isc:iec,jsc:jec,2) )
    call mpp_update_domains(Ice%part_size(:,:,1:2), Domain ) ! these two updates cannot be combined
    call mpp_update_domains(IST%h_ice(:,:,2), Domain )       ! these two updates cannot be combined
    IST%part_size_uv(:,:,:) = 0.0
    IST%part_size_uv(:,:,1) = 1.0
    do k=2,km
      call t_to_uv(Ice%part_size(:,:,k),IST%part_size_uv(:,:,k), Ice%G)
    enddo
    do k=2,km ; do j=jsc,jec ; do i=isc,iec
      IST%part_size_uv (i,j,1) = IST%part_size_uv(i,j,1)-IST%part_size_uv (i,j,k)
    enddo ; enddo ; enddo
    do_init = .false.
  endif

  if (first_time .and. conservation_check) then
    h2o    = 0.0
    h2o(1) = sum(cell_area*all_avg(DI*IST%h_ice(isc:iec,jsc:jec,:) +   &
                                   DS*IST%h_snow(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)))
    heat   = 0.0

    do k=2,km ; do j=jsc,jec ; do i=isc,iec
      if ((Ice%part_size(i,j,k)>0.0) .and. (IST%h_ice(i,j,k)>0.0)) then
        if (slab_ice) then
          heat(1) = heat(1) - cell_area(i,j) * Ice%part_size(i,j,k) * IST%h_ice(i,j,2)*DI*LI
        else
          heat(1) = heat(1) - cell_area(i,j) * Ice%part_size(i,j,k) *        &
                    e_to_melt(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k),         &
                    IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4))
        endif
      endif
    enddo ; enddo ; enddo

    salt = 0.0
    salt(1) = sum(cell_area*all_avg(DI*IST%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))) !*ice_bulk_salin 

    first_time = .false.
  endif

  sst(:,:) = Ice%t_surf(:,:,1)
!
! temporary check to find escomp problem
!
do j=jsc,jec
  do i=isc,iec     
    if (Ice%mask(i,j)) then
      if ((Ice%part_size(i,j,1) > 0.0).and.(is_NaN(Ice%t_surf(i,j,1))&
                                        .or.Ice%t_surf(i,j,1)<Tfreeze-5.0&
                                        .or.Ice%t_surf(i,j,1)>Tfreeze+40.0)) &
         print '(a,f10.2,2i5,1PG12.6)','BAD SST',Ice%t_surf(i,j,1),i,j,&
                                             Ice%part_size(i,j,1)
      do k = 2, km
       if ((Ice%part_size(i,j,k) > 0.0).and.(is_NaN(Ice%t_surf(i,j,k))&
                                        .or.Ice%t_surf(i,j,k)<Tfreeze-70.0&
                                        .or.Ice%t_surf(i,j,k)>Tfreeze+0.1)) &
         print '(a,f10.2,3i5,1PG12.6,2f10.2)','BAD ICE',Ice%t_surf(i,j,k),i,j,k,&
                                             Ice%part_size(i,j,k),&
                                             IST%h_ice(i,j,k),IST%h_snow(i,j,k)
      enddo
    endif
  enddo
enddo
! end temporary check

  if (present(s_surf_ice_bot)) then
    do j=jsc,jec ; do i=isc,iec     
      Ice%s_surf(i,j) = s_surf_ice_bot(i,j)
    enddo ; enddo
  else
    do j=jsc,jec ; do i=isc,iec  
      Ice%s_surf(i,j) = 0.0
    enddo ; enddo
  endif

  do j=jsc,jec ; do i=isc,iec
    IST%frazil(i,j) = frazil_ice_bot(i,j)
  enddo ; enddo

!       transfer the ocean state for extra tracer fluxes
!
  do n = 1, Ocean_ice_boundary%fields%num_bcs  !{
    do m = 1, Ocean_ice_boundary%fields%bc(n)%num_fields  !{
      Ice%ocean_fields%bc(n)%field(m)%values(:,:,1) = Ocean_ice_boundary%fields%bc(n)%field(m)%values
    enddo  !} m
  enddo  !} n

  if (present(sea_lev_ice_bot)) then
     do j=jsc,jec ; do i=isc,iec  
       IST%sea_lev(i,j) = sea_lev_ice_bot(i,j)
     enddo ; enddo
  else
     do j=jsc,jec ; do i=isc,iec
       IST%sea_lev(i,j) = 0.0
     enddo ; enddo
  endif

  call mpp_update_domains(IST%sea_lev, Domain)

  do k=2,km ; do j=jsc,jec ; do i=isc,iec  
    IST%tmelt(i,j,k) = 0.0
    IST%bmelt(i,j,k) = 0.0
  enddo ; enddo ; enddo

  tmp = ice_avg(IST%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:) )
  do j=jsc,jec ; do i=isc,iec
    if ( tmp(i,j) > 0.0) then
       IST%bheat(i,j) = kmelt*(sst(i,j)-Tfreeze+MU_TS*Ice%s_surf(i,j))
    else
       IST%bheat(i,j) = 0.0
    endif
  enddo ; enddo

  Ice%ice_mask(:,:,1) = .false.
  do k=2,km ; do j=jsc,jec ; do i=isc,iec
    Ice%ice_mask(i,j,k) = (IST%h_ice(i,j,k) > 0.0)
  enddo ; enddo ; enddo

  do k=2,km ; do j=jsc,jec ; do i=isc,iec ; if (IST%h_ice(i,j,k) > 0.0) then
    call ice_optics(IST%h_snow(i,j,k), IST%h_ice(i,j,k), &
         Ice%t_surf(i,j,k)-Tfreeze, -MU_TS*Ice%s_surf(i,j), &
         Ice%albedo_vis_dir(i,j,k), Ice%albedo_vis_dif(i,j,k), &
         Ice%albedo_nir_dir(i,j,k), Ice%albedo_nir_dif(i,j,k), &
         IST%sw_abs_sfc(i,j,k), &
         IST%sw_abs_snow(i,j,k), IST%sw_abs_ice(i,j,k,1), &
         IST%sw_abs_ice(i,j,k,2), IST%sw_abs_ice(i,j,k,3), &
         IST%sw_abs_ice(i,j,k,4), IST%sw_abs_ocn(i,j,k) , &
         IST%sw_abs_int(i,j,k),  IST%pen(i,j,k), IST%trn(i,j,k), &
         coszen_in=IST%coszen(i,j))
    !
    !Niki: Is the following correct for diagnostics?
    Ice%albedo(i,j,k)=(Ice%albedo_vis_dir(i,j,k)+Ice%albedo_nir_dir(i,j,k)&
                      +Ice%albedo_vis_dif(i,j,k)+Ice%albedo_nir_dif(i,j,k))/4

  endif ; enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    IST%u_ocn(i,j) = u_surf_ice_bot(i,j) ! need under-ice current
    IST%v_ocn(i,j) = v_surf_ice_bot(i,j) ! for water drag term
  enddo ; enddo

  call mpp_update_domains(IST%u_ocn, IST%v_ocn, Domain, gridtype=BGRID_NE)

  ! put ocean and ice velocities into Ice%u_surf/v_surf on t-cells
  call uv_to_t(IST%u_ocn, Ice%u_surf(:,:,1), Ice%G)
  call uv_to_t(IST%v_ocn, Ice%v_surf(:,:,1), Ice%G)

  call uv_to_t(IST%u_ice, Ice%u_surf(:,:,2), Ice%G)
  call uv_to_t(IST%v_ice, Ice%v_surf(:,:,2), Ice%G)

  do k=1,2 ; do j=jsc,jec ; do i=isc,iec
    u = Ice%u_surf(i,j,k)
    v = Ice%v_surf(i,j,k)
    Ice%u_surf(i,j,k) =  u*cos_rot(i,j)+v*sin_rot(i,j) ! rotate velocity from ocean
    Ice%v_surf(i,j,k) =  v*cos_rot(i,j)-u*sin_rot(i,j) ! coord. to lat/lon coord.
  enddo ; enddo ; enddo

  do k=3,km ; do j=jsc,jec ; do i=isc,iec
    Ice%u_surf(i,j,k) = Ice%u_surf(i,j,2)  ! same ice flow on all ice partitions
    Ice%v_surf(i,j,k) = Ice%v_surf(i,j,2)  !
  enddo ; enddo ; enddo
  !
  ! Pre-timestep diagnostics
  !
  if (id_sst>0) sent = send_data(id_sst, sst(isc:iec,jsc:jec)-Tfreeze, Ice%Time, mask=Ice%mask)
  if (id_sss>0) sent = send_data(id_sss, Ice%s_surf(isc:iec,jsc:jec) , Ice%Time, mask=Ice%mask)
  if (id_ssh>0) sent = send_data(id_ssh, IST%sea_lev(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
  if (id_uo >0) sent = send_data(id_uo , IST%u_ocn(isc:iec,jsc:jec)  , Ice%Time, mask=Ice%mask)
  if (id_vo >0) sent = send_data(id_vo , IST%v_ocn(isc:iec,jsc:jec)  , Ice%Time, mask=Ice%mask)
  if (id_bheat>0) sent = send_data(id_bheat, IST%bheat(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)

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
  
  real, dimension(isc:iec,jsc:jec,km) :: flux_t, flux_q, flux_lh, flux_lw
  real, dimension(isc:iec,jsc:jec,km) :: flux_sw_nir_dir
  real, dimension(isc:iec,jsc:jec,km) :: flux_sw_nir_dif
  real, dimension(isc:iec,jsc:jec,km) :: flux_sw_vis_dir
  real, dimension(isc:iec,jsc:jec,km) :: flux_sw_vis_dif
  real, dimension(isc:iec,jsc:jec,km) :: flux_u, flux_v, lprec, fprec
  real, dimension(isc:iec,jsc:jec,km) :: dhdt, dedt, drdt ! d(flux)/d(surf_temp) (+ up)
  integer                             :: dy, sc, i, j, k
  real                                :: dt_fast, ts_new, dts, hf, hfd, latent
  logical                             :: sent
  real                                :: gmt, time_since_ae, cosz, rrsun, fracday, fracday_dt_ice, fracday_day
  real, dimension(isc:iec,jsc:jec)    :: diurnal_factor
  real                                :: rad, cosz_day, cosz_dt_ice, rrsun_day, rrsun_dt_ice
  type (time_type)                    :: Dt_ice
  real, dimension(isc:iec,jsc:jec)    :: cosz_alb
  real                                :: flux_sw ! sum over dir/dif vis/nir components

  real, dimension(isc:iec,jsc:jec,km) :: albedo_vis_dir,albedo_vis_dif,albedo_nir_dir,albedo_nir_dif
  real, dimension(isc:iec,jsc:jec,km) :: sw_abs_sfc,sw_abs_snow,sw_abs_ice1,sw_abs_ice2,sw_abs_ice3,sw_abs_ice4,sw_abs_ocn,sw_abs_int,pen,trn

  if (id_alb>0) call post_avg(id_alb, Ice%albedo(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
!  if (id_alb>0) sent = send_data(id_alb, all_avg(Ice%albedo(isc:iec,jsc:jec,:), &
!       Ice%part_size(isc:iec,jsc:jec,:)), Ice%Time, mask=Ice%mask)

  do j=jsc,jec ; do i=isc,iec
    IST%coszen(i,j) = Atmos_boundary%coszen(i,j,1)
    Ice%p_surf(i,j) = Atmos_boundary%p(i,j,1)
  enddo ; enddo

  !   Set up local copies of fluxes.  The Atmos_boundary arrays may have
  ! different index conventions than are used internally in this component.
  do k=1,km ; do j=jsc,jec ; do i=isc,iec
    flux_u(i,j,k)  = Atmos_boundary%u_flux(i,j,k)
    flux_v(i,j,k)  = Atmos_boundary%v_flux(i,j,k)
    flux_t(i,j,k)  = Atmos_boundary%t_flux(i,j,k)
    flux_q(i,j,k)  = Atmos_boundary%q_flux(i,j,k)
    flux_lh(i,j,k) = hlv*Atmos_boundary%q_flux(i,j,k)
    flux_lw(i,j,k) = Atmos_boundary%lw_flux(i,j,k)
    flux_sw_nir_dir(i,j,k) = Atmos_boundary%sw_flux_nir_dir(i,j,k)
    flux_sw_nir_dif(i,j,k) = Atmos_boundary%sw_flux_nir_dif(i,j,k)
    flux_sw_vis_dir(i,j,k) = Atmos_boundary%sw_flux_vis_dir(i,j,k)
    flux_sw_vis_dif(i,j,k) = Atmos_boundary%sw_flux_vis_dif(i,j,k)
    lprec(i,j,k)   = Atmos_boundary%lprec(i,j,k)
    fprec(i,j,k)   = Atmos_boundary%fprec(i,j,k)
    dhdt(i,j,k) = Atmos_boundary%dhdt(i,j,k)
    dedt(i,j,k) = Atmos_boundary%dedt(i,j,k)
    drdt(i,j,k) = Atmos_boundary%drdt(i,j,k)
  enddo ; enddo ; enddo

  if (add_diurnal_sw .or. do_sun_angle_for_alb) then
!---------------------------------------------------------------------
!    extract time of day (gmt) from time_type variable time with
!    function universal_time.
!---------------------------------------------------------------------
    gmt = universal_time(Ice%Time)
!---------------------------------------------------------------------
!    extract the time of year relative to the northern hemisphere
!    autumnal equinox (time_since_ae) from time_type variable 
!    time using the function orbital_time.
!---------------------------------------------------------------------
    time_since_ae = orbital_time(Ice%Time)
!--------------------------------------------------------------------
!    call diurnal_solar_2d to calculate astronomy fields
!    convert G%geoLonT and G%geoLatT to radians
!    Per Rick Hemler:
!      call daily_mean_solar to get cosz (over a day)
!      call diurnal_solar with dtime=Dt_ice to get cosz over Dt_ice
!      diurnal_factor = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice/cosz_day*fracday_day*rrsun_day
!--------------------------------------------------------------------
    rad = acos(-1.)/180.
    Dt_ice = Ice%Time_step_fast
  endif
  if (add_diurnal_sw) then
    do j=jsc,jec ; do i=isc,iec
      call diurnal_solar(Ice%G%geoLatT(i,j)*rad, Ice%G%geoLonT(i,j)*rad, Ice%time, cosz=cosz_dt_ice, &
                         fracday=fracday_dt_ice, rrsun=rrsun_dt_ice, dt_time=Dt_ice)
      call daily_mean_solar (Ice%G%geoLatT(i,j)*rad, time_since_ae, cosz_day, fracday_day, rrsun_day)
      diurnal_factor(i,j) = cosz_dt_ice*fracday_dt_ice*rrsun_dt_ice /   &
                                   max(1e-30, cosz_day*fracday_day*rrsun_day)
    enddo ; enddo

    do k=1,km ; do j=jsc,jec ; do i=isc,iec
      flux_sw_nir_dir(i,j,k) = flux_sw_nir_dir(i,j,k) * diurnal_factor(i,j)
      flux_sw_nir_dif(i,j,k) = flux_sw_nir_dif(i,j,k) * diurnal_factor(i,j)
      flux_sw_vis_dir(i,j,k) = flux_sw_vis_dir(i,j,k) * diurnal_factor(i,j)
      flux_sw_vis_dif(i,j,k) = flux_sw_vis_dif(i,j,k) * diurnal_factor(i,j)
    enddo ; enddo ; enddo
  endif

  !
  ! implicit update of ice surface temperature
  !
  call get_time(Ice%Time_step_fast,sc,dy); dt_fast = 864e2*dy+sc

  do k=2,km ; do j=jsc,jec ; do i=isc,iec
    if (IST%h_ice(i,j,k) > 0.0) then
      if (slab_ice) then
        flux_lh(i,j,k) = flux_lh(i,j,k) + hlf*flux_q(i,j,k)
        latent             = hlv+hlf
      elseif (IST%h_snow(i,j,k)>0.0) then
        flux_lh(i,j,k) = flux_lh(i,j,k) + (hlf-CI*IST%t_snow(i,j,k))*flux_q(i,j,k)
        latent             = hlv + (hlf-CI*IST%t_snow(i,j,k))
      else
        flux_lh(i,j,k) = flux_lh(i,j,k)+hlf*flux_q(i,j,k)*(1-TFI/IST%t_ice(i,j,k,1))
        latent             = hlv + hlf*(1-TFI/IST%t_ice(i,j,k,1))
      endif
      !### ADD PARENTHESES.
      flux_sw = flux_sw_vis_dir(i,j,k)+flux_sw_vis_dif(i,j,k) &
               +flux_sw_nir_dir(i,j,k)+flux_sw_nir_dif(i,j,k)
      hfd = dhdt(i,j,k) + dedt(i,j,k)*latent + drdt(i,j,k)
      hf  = flux_t(i,j,k) + flux_q(i,j,k)*latent - flux_lw(i,j,k)   &
            - (1-IST%pen(i,j,k))*flux_sw - hfd*(Ice%t_surf(i,j,k)-Tfreeze)
      call ice5lay_temp(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k),    &
                        IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3),   &
                        IST%t_ice(i,j,k,4), ts_new, hf, hfd,                        &
                        IST%sw_abs_snow(i,j,k)*flux_sw, &
                        IST%sw_abs_ice(i,j,k,1)*flux_sw, &
                        IST%sw_abs_ice(i,j,k,2)*flux_sw, &
                        IST%sw_abs_ice(i,j,k,3)*flux_sw, &
                        IST%sw_abs_ice(i,j,k,4)*flux_sw, &
                        -MU_TS*Ice%s_surf(i,j), IST%bheat(i,j), dt_fast, &
                        IST%tmelt(i,j,k), IST%bmelt(i,j,k))
      dts                = ts_new-(Ice%t_surf(i,j,k)-Tfreeze)
      Ice%t_surf(i,j,k)  = Ice%t_surf(i,j,k)  + dts
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
  if (do_sun_angle_for_alb) then
    call diurnal_solar(Ice%G%geoLatT(isc:iec,jsc:jec)*rad, Ice%G%geoLonT(isc:iec,jsc:jec)*rad, &
                 Ice%time, cosz=cosz_alb, fracday=diurnal_factor, rrsun=rrsun_dt_ice, dt_time=Dt_ice)  !diurnal_factor as dummy
    call compute_ocean_albedo (Ice%mask, cosz_alb(:,:), Ice%albedo_vis_dir(:,:,1),&
                                 Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                                 Ice%albedo_nir_dif(:,:,1), latitude )
  else
    call compute_ocean_albedo (Ice%mask, IST%coszen(:,:), Ice%albedo_vis_dir(:,:,1),&
                               Ice%albedo_vis_dif(:,:,1), Ice%albedo_nir_dir(:,:,1),&
                               Ice%albedo_nir_dif(:,:,1), latitude )
  endif

  call sum_top_quantities(Ice, Ice%Ice_state, Atmos_boundary%fluxes, flux_u, flux_v, flux_t, &
    flux_q, flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif, &
    flux_lw, lprec, fprec, flux_lh, G )

  Ice%Time = Ice%Time + Ice%Time_step_fast ! advance time

  if (id_alb_vis_dir>0) call post_avg(id_alb_vis_dir, Ice%albedo_vis_dir(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_alb_vis_dif>0) call post_avg(id_alb_vis_dif, Ice%albedo_vis_dif(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_alb_nir_dir>0) call post_avg(id_alb_nir_dir, Ice%albedo_nir_dir(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_alb_nir_dif>0) call post_avg(id_alb_nir_dif, Ice%albedo_nir_dif(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)

  if (id_sw_abs_snow>0) call post_avg(id_sw_abs_snow, IST%sw_abs_snow(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw_abs_ice1>0) call post_avg(id_sw_abs_ice1, IST%sw_abs_ice(isc:iec,jsc:jec,:,1), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw_abs_ice2>0) call post_avg(id_sw_abs_ice4, IST%sw_abs_ice(isc:iec,jsc:jec,:,2), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw_abs_ice3>0) call post_avg(id_sw_abs_ice4, IST%sw_abs_ice(isc:iec,jsc:jec,:,3), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw_abs_ice4>0) call post_avg(id_sw_abs_ice4, IST%sw_abs_ice(isc:iec,jsc:jec,:,4), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)

  if (id_sw_pen>0) call post_avg(id_sw_pen, IST%pen(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)
  if (id_sw_trn>0) call post_avg(id_sw_trn, IST%trn(isc:iec,jsc:jec,:), &
                             Ice%part_size(isc:iec,jsc:jec,:), IST%diag, mask=Ice%mask)

  if (id_coszen>0) sent = send_data(id_coszen, IST%coszen(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)

end subroutine do_update_ice_model_fast

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! update_ice_model_slow - do ice dynamics, transport, and mass changes         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine update_ice_model_slow(Ice, IST, G, runoff, calving, &
                                 runoff_hflx, calving_hflx)

  type(ice_data_type),                intent(inout) :: Ice
  type(ice_state_type),               intent(inout) :: IST
  type(sea_ice_grid_type),            intent(inout) :: G
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: runoff, calving
    real, dimension(isc:iec,jsc:jec), intent(in), optional :: runoff_hflx, calving_hflx

    real, dimension(isd:ied,jsd:jed)      :: fx_wat, fy_wat
    real, dimension(isc:iec,jsc:jec)      :: hi_change, h2o_change, bsnk, tmp2d, mass
    real, dimension(isc:iec,jsc:jec,2:km) :: snow_to_ice
    real, dimension(isc:iec,jsc:jec,km)   :: part_save, part_save_uv
    real, dimension(isc:iec,jsc:jec)      :: dum1, Obs_h_ice ! for qflux calculation
    real, dimension(isc:iec,jsc:jec,2)    :: Obs_cn_ice      ! partition 2 = ice concentration
    real, dimension(isd:ied,jsd:jed)      :: tmp1, tmp2
    real, dimension(isd:ied,jsd:jed)      :: wind_stress_x, wind_stress_y
    real, dimension(2:km)                 :: e2m
  real, dimension(isd:ied,jsd:jed) :: uc, vc ! Ice velocities interpolated onto
                                             ! a C-grid, in m s-1.
  real, dimension(isd:ied,jsd:jed) :: diagVar ! An temporary array for diagnostics.

    integer                               :: i, j, k, l, sc, dy, iyr, imon, iday, ihr, imin, isec
    real                                  :: dt_slow, heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic, bablt
    real                                  :: heat_limit_ice, heat_res_ice
    real                                  :: tot_heat, heating, tot_frazil
    logical                               :: sent
    real, parameter                       :: LI = hlf

  call get_time(Ice%Time_step_slow,sc,dy); dt_slow = 864e2*dy+sc
  !
  ! Set up fluxes
  !
  if (present(runoff)) then ! save liquid runoff for ocean
    do j=jsc,jec ; do i=isc,iec
      Ice%runoff(i,j)  = runoff(i,j)
    enddo ; enddo
    if (id_runoff>0) &
      sent = send_data(id_runoff, Ice%runoff(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask )
  else
    do j=jsc,jec ; do i=isc,iec ; Ice%runoff(i,j) = 0.0 ; enddo ; enddo
  endif

  if (present(calving)) then ! save frozen runoff for ocean
    do j=jsc,jec ; do i=isc,iec
      Ice%calving(i,j) = calving(i,j)
    enddo ; enddo
    if (id_calving>0) &
      sent = send_data(id_calving, Ice%calving(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask )
  else
    do j=jsc,jec ; do i=isc,iec ; Ice%calving(i,j) = 0.0 ; enddo ; enddo
  endif

  if (present(runoff_hflx)) then ! save liquid runoff hflx for ocean
    do j=jsc,jec ; do i=isc,iec
      Ice%runoff_hflx(i,j)  = runoff_hflx(i,j)
    enddo ; enddo
    if (id_runoff_hflx>0) &
      sent = send_data( id_runoff_hflx, Ice%runoff_hflx(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask )
  else
    do j=jsc,jec ; do i=isc,iec ; Ice%runoff_hflx(i,j) = 0.0 ; enddo ; enddo
  endif

  if (present(calving_hflx)) then ! save frozen runoff hflx for ocean
    do j=jsc,jec ; do i=isc,iec
      Ice%calving_hflx(i,j) = calving_hflx(i,j)
    enddo ; enddo
    if (id_calving_hflx>0) &
      sent = send_data(id_calving_hflx, Ice%calving_hflx(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask )
  else
    do j=jsc,jec ; do i=isc,iec ; Ice%calving_hflx(i,j) = 0.0 ; enddo ; enddo
  endif

  !TOM> assume that open water area is not up to date:
  call mpp_clock_end(iceClock)
  call mpp_clock_end(iceClock2)
  tmp1(:,:) = 1.-max(1.-sum(Ice%part_size(:,:,2:km),dim=3),0.0)
  tmp2(:,:) = ice_avg(IST%h_ice,Ice%part_size)
  ! Calve off icebergs and integrate forward iceberg trajectories
  if (do_icebergs) call icebergs_run( Ice%icebergs, Ice%Time,                &
                    Ice%calving, IST%u_ocn, IST%v_ocn, IST%u_ice, IST%v_ice, &
                    Ice%flux_u, Ice%flux_v, IST%sea_lev, Ice%t_surf(:,:,1),  &
                    Ice%calving_hflx, tmp1, tmp2)
  call mpp_clock_begin(iceClock2)
  call mpp_clock_begin(iceClock)

  call enable_SIS_averaging(dt_slow, Ice%Time, IST%diag)
  call avg_top_quantities(Ice, Ice%Ice_state, G) ! average fluxes from update_ice_model_fast
  call disable_SIS_averaging(IST%diag)

  do k=1,km ; do j=jsc,jec ; do i=isc,iec
    part_save(i,j,k)    = Ice%part_size(i,j,k)
    part_save_uv(i,j,k) = IST%part_size_uv(i,j,k)
  enddo ; enddo ; enddo
  !
  ! conservation checks: top fluxes
  !
  call mpp_clock_begin(iceClock7)
  if (conservation_check) then
    ! Regroup for efficiency?
    h2o(2) = h2o(2)+dt_slow*sum(cell_area*(Ice%runoff(isc:iec,jsc:jec) &
                                          +Ice%calving(isc:iec,jsc:jec)&
            +all_avg(IST%lprec_top(isc:iec,jsc:jec,:)+IST%fprec_top(isc:iec,jsc:jec,:)&
                    -IST%flux_q_top(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))),&
             cell_area > 0 )
    heat(2) = heat(2)+dt_slow*sum(cell_area*(                              &
              all_avg(IST%flux_sw_vis_dir_top(isc:iec,jsc:jec,:)+IST%flux_sw_vis_dif_top(isc:iec,jsc:jec,:) &
             +IST%flux_sw_nir_dir_top(isc:iec,jsc:jec,:)+IST%flux_sw_nir_dif_top(isc:iec,jsc:jec,:), &
              Ice%part_size(isc:iec,jsc:jec,:))                            &
             +all_avg(IST%flux_lw_top(isc:iec,jsc:jec,:), Ice%part_size(isc:iec,jsc:jec,:))   &
             -all_avg(IST%flux_t_top(isc:iec,jsc:jec,:) , Ice%part_size(isc:iec,jsc:jec,:))   &
             -all_avg(IST%flux_lh_top(isc:iec,jsc:jec,:), Ice%part_size(isc:iec,jsc:jec,:))   &
             -LI*(all_avg(IST%fprec_top,Ice%part_size(isc:iec,jsc:jec,:))  &
             +Ice%calving(isc:iec,jsc:jec))), cell_area > 0 )
    tot_frazil = sum(cell_area*IST%frazil)
    !Niki: Does runoff or calving bring in salt? Or is salt(2) = 0.0
  endif
  call mpp_clock_end(iceClock7)

  ! Dynamics
  !

  call mpp_clock_begin(iceClock4)
  tmp1(:,:) = ice_avg(IST%h_snow,Ice%part_size)
  tmp2(:,:) = ice_avg(IST%h_ice,Ice%part_size)
  wind_stress_x(:,:) = 0.0 ; wind_stress_y(:,:) = 0.0
  wind_stress_x(isc:iec,jsc:jec) = ice_avg(IST%flux_u_top_bgrid(isc:iec,jsc:jec,:), IST%part_size_uv(isc:iec,jsc:jec,:))
  wind_stress_y(isc:iec,jsc:jec) = ice_avg(IST%flux_v_top_bgrid(isc:iec,jsc:jec,:), IST%part_size_uv(isc:iec,jsc:jec,:))

  call enable_SIS_averaging(dt_slow, Ice%Time, IST%diag)
  call mpp_clock_begin(iceClocka)
  call ice_dynamics(1-Ice%part_size(:,:,1), tmp1, tmp2, IST%u_ice, IST%v_ice, &
                    IST%u_ocn, IST%v_ocn, &
                    wind_stress_x, wind_stress_y, IST%sea_lev, fx_wat, fy_wat, &
                    dt_slow, Ice%G, IST%ice_dyn_CSp)
  call mpp_clock_end(iceClocka)

  call mpp_clock_begin(iceClockb)
  call mpp_update_domains(IST%u_ice, IST%v_ice, Domain, gridtype=BGRID_NE)
  call mpp_clock_end(iceClockb)

  call mpp_clock_begin(iceClockc)
  !
  ! Dynamics diagnostics
  !
  if (id_fax>0) call post_avg(id_fax, IST%flux_u_top_bgrid(isc:iec,jsc:jec,:), &
                              IST%part_size_uv(isc:iec,jsc:jec,:), IST%diag)
  if (id_fay>0) call post_avg(id_fay, IST%flux_v_top_bgrid(isc:iec,jsc:jec,:), &
                              IST%part_size_uv(isc:iec,jsc:jec,:), IST%diag)
  call disable_SIS_averaging(IST%diag)


  do k=2,km ; do j=jsc,jec ; do i=isc,iec
    IST%flux_u_top_bgrid(I,J,k) = fx_wat(I,J)  ! stress of ice on ocean
    IST%flux_v_top_bgrid(I,J,k) = fy_wat(I,J)  !
  enddo ; enddo ; enddo
  call mpp_clock_end(iceClockc)
  call mpp_clock_end(iceClock4)

  !
  ! Thermodynamics
  !
  call mpp_clock_begin(iceClock5)
  if (id_frazil>0) sent = send_data(id_frazil, IST%frazil(isc:iec,jsc:jec)/dt_slow,  Ice%Time, mask=Ice%mask)
  snow_to_ice(:,:,:) = 0.
  bsnk(:,:) = 0.

  call get_avg(IST%h_ice(isc:iec,jsc:jec,:), Ice%part_size(isc:iec,jsc:jec,2:), hi_change(:,:))
  h2o_change(:,:) = all_avg(DS*IST%h_snow(isc:iec,jsc:jec,:)+DI*IST%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))

  ! get observed ice thickness for ice restoring, if calculating qflux
  if (do_ice_restore) &
       call get_sea_surface(Ice%Time, dum1, Obs_cn_ice, Obs_h_ice)
  do k=2,km ; do j=jsc,jec ; do i=isc,iec 
    if (cell_area(i,j) > 0 .and.IST%h_ice(i,j,k) > 0) then
      ! reshape the ice based on fluxes

      h2o_from_ocn = 0.0; h2o_to_ocn = 0.0; heat_to_ocn = 0.0
      call ice5lay_resize(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k),&
                          IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2),              &
                          IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4),              &
                          IST%fprec_top(i,j,k) *dt_slow, 0.0,                &
                          IST%flux_q_top(i,j,k)*dt_slow,                     &
                          IST%tmelt (i,j,k), IST%bmelt(i,j,k),               &
                          -MU_TS*Ice%s_surf(i,j),                            &
                          heat_to_ocn, h2o_to_ocn, h2o_from_ocn,             &
                          snow_to_ice(i,j,k), bablt                          )

      ! modify above-ice to under-ice fluxes for passing to ocean
      IST%flux_q_top (i,j,k) = h2o_from_ocn/dt_slow ! no ice, evaporation left
      IST%flux_lh_top(i,j,k) = hlv*IST%flux_q_top(i,j,k)
      IST%flux_lw_top(i,j,k) = 0.0
      IST%flux_t_top (i,j,k) = IST%bheat(i,j)-heat_to_ocn/dt_slow
      IST%flux_sw_vis_dif_top(i,j,k) = (IST%flux_sw_vis_dir_top(i,j,k)+      &
            IST%flux_sw_vis_dif_top(i,j,k)+IST%flux_sw_nir_dir_top(i,j,k)+   &
            IST%flux_sw_nir_dif_top(i,j,k))*IST%pen(i,j,k)*IST%trn(i,j,k)
      IST%flux_sw_nir_dir_top(i,j,k) = 0.0
      IST%flux_sw_nir_dif_top(i,j,k) = 0.0
      IST%flux_sw_vis_dir_top(i,j,k) = 0.0
      IST%fprec_top  (i,j,k) = 0.0
      IST%lprec_top  (i,j,k) = IST%lprec_top(i,j,k) + h2o_to_ocn/dt_slow

      bsnk(i,j) = bsnk(i,j) - Ice%part_size(i,j,k)*bablt ! bot. melt. ablation

    endif

    !
    ! absorb frazil in thinest ice partition available
    !
    if (IST%frazil(i,j)>0 .and. Ice%part_size(i,j,1)+Ice%part_size(i,j,k)>0.01) then
      !                                                           was ...>0.0
      ! raised above threshold from 0 to 0.01 to avert ocean-ice model blow-ups
      !
      IST%h_snow(i,j,k) = IST%h_snow(i,j,k) * Ice%part_size(i,j,k)
      IST%h_ice(i,j,k)  = IST%h_ice(i,j,k)  * Ice%part_size(i,j,k)
      Ice%t_surf(i,j,k) = (Ice%t_surf(i,j,k) * Ice%part_size(i,j,k) &
                           +(Tfreeze - MU_TS*Ice%s_surf(i,j))*Ice%part_size(i,j,1))
      Ice%part_size(i,j,k) = Ice%part_size(i,j,k) + Ice%part_size(i,j,1)
      Ice%part_size(i,j,1) = 0.0
      IST%h_snow(i,j,k) = IST%h_snow(i,j,k) / Ice%part_size(i,j,k)
      IST%h_ice(i,j,k)  = IST%h_ice(i,j,k) / Ice%part_size(i,j,k)
      Ice%t_surf(i,j,k) = Ice%t_surf(i,j,k) / Ice%part_size(i,j,k)

      call ice5lay_resize(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k), &
                          IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2),                &
                          IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4), 0.0,           &
                          IST%frazil(i,j) / Ice%part_size(i,j,k), 0.0, 0.0, 0.0, &
                          -MU_TS*Ice%s_surf(i,j),                              &
                          heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic)
      IST%frazil(i,j) = 0.0;
      !
      ! spread frazil salinification over all partitions
      !
      IST%lprec_top  (i,j,:) = IST%lprec_top(i,j,:) + h2o_to_ocn*Ice%part_size(i,j,k)/dt_slow
    endif

  enddo ; enddo ; enddo   ! i-, j-, and k-loops
  call mpp_clock_end(iceClock5)

  !
  ! Calculate QFLUX's from (1) restoring to obs and (2) limiting total ice.
  !
  call mpp_clock_begin(iceClock6)
  if (do_ice_restore .or. do_ice_limit) then
    do j=jsc,jec ; do i=isc,iec
      IST%qflx_lim_ice(i,j) = 0.0
      IST%qflx_res_ice(i,j) = 0.0
    enddo ; enddo

    do j=jsc,jec ; do i=isc,iec
      heat_res_ice   = 0.0
      heat_limit_ice = 0.0
      !
      ! calculate enthalpy
      !
      if (slab_ice) then
        e2m(2) = IST%h_ice(i,j,2)*DI*LI
      else
        do k=2,km
          if ((Ice%part_size(i,j,k)>0.0.and.IST%h_ice(i,j,k)>0.0)) then
             e2m(k) = e_to_melt(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k), &
                      IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3), &
                      IST%t_ice(i,j,k,4)) * Ice%part_size(i,j,k)
          else
             e2m(k) = 0.0
          endif
        enddo
      endif
      !
      ! calculate heat needed to constrain ice enthalpy
      !
      if (do_ice_restore) then
        ! TK Mod: restore to observed enthalpy (no change for slab, but for
        !         sis ice, results in restoring toward thickness * concentration      

        ! TK Mod for test 10/18/02
        !   Concentration is 1.0 where there is ice for slab case.
        !   The input field Obs_cn_ice may have values less than 1,
        !     so put in if test...

        if (slab_ice) then
          heat_res_ice = -(LI*DI*Obs_h_ice(i,j)-sum(e2m)) &
                         *dt_slow/(86400*ice_restore_timescale)
        else                   
          heat_res_ice = -(LI*DI*Obs_h_ice(i,j)*Obs_cn_ice(i,j,2)-sum(e2m)) &
                         *dt_slow/(86400*ice_restore_timescale)
        endif

        !        heat_res_ice = -(LI*DI*Obs_h_ice(i,j)-sum(e2m)) &
        !                        *dt_slow/(86400*ice_restore_timescale)

      endif

      if (do_ice_limit .and. (sum(e2m) > max_ice_limit*DI*LI)) then
        heat_limit_ice = sum(e2m)-LI*DI*max_ice_limit
        ! should we "heat_ice_res = 0.0" ?
      endif

      !
      ! apply constraining heat to ice
      !
      tot_heat = heat_res_ice+heat_limit_ice
      if (slab_ice) IST%h_ice(i,j,2) = IST%h_ice(i,j,2) - tot_heat/(DI*LI)

      if (.not. slab_ice .and. (tot_heat>0.0)) then  ! add like ocean-ice heat
        do k=2,km
          if (e2m(k) > 0.0) then
            heating = tot_heat/sum(Ice%part_size(i,j,k:km))
            if (heating*Ice%part_size(i,j,k) > e2m(k)) then ! cat. melts away
              IST%h_ice (i,j,k) = 0.0
              IST%h_snow(i,j,k) = 0.0
              tot_heat = tot_heat - e2m(k)
            else
              h2o_from_ocn = 0.0; h2o_to_ocn = 0.0; heat_to_ocn = 0.0
              call ice5lay_resize(IST%h_snow(i,j,k), IST%t_snow(i,j,k),  &
                                  IST%h_ice(i,j,k), IST%t_ice(i,j,k,1),   &
                                  IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3),  &
                                  IST%t_ice(i,j,k,4),  0.0, 0.0, 0.0, 0.0,&
                                  heating, -MU_TS*Ice%s_surf(i,j),       &
                                  heat_to_ocn, h2o_to_ocn, h2o_from_ocn, &
                                  snow_to_ice(i,j,k), bablt              )
              tot_heat = tot_heat - heating*Ice%part_size(i,j,k)
            endif
          endif
        enddo
      endif

      tot_heat = heat_res_ice+heat_limit_ice
      if (.not. slab_ice .and. (tot_heat<0.0)) then ! add like frazil
        do k=2,km
          if (Ice%part_size(i,j,1)+Ice%part_size(i,j,k)>0) exit
        enddo
        ! k is thinnest ice partition that can recieve frazil
        IST%h_snow(i,j,k) = IST%h_snow(i,j,k) * Ice%part_size(i,j,k)
        IST%h_ice(i,j,k)  = IST%h_ice(i,j,k)  * Ice%part_size(i,j,k)
        Ice%t_surf(i,j,k) = Ice%t_surf(i,j,k) * Ice%part_size(i,j,k) &
                       + (Tfreeze - MU_TS*Ice%s_surf(i,j))*Ice%part_size(i,j,1)
        Ice%part_size(i,j,k) = Ice%part_size(i,j,k) + Ice%part_size(i,j,1)
        Ice%part_size(i,j,1) = 0.0
        IST%h_snow(i,j,k) = IST%h_snow(i,j,k) / Ice%part_size(i,j,k)
        IST%h_ice(i,j,k) =  IST%h_ice(i,j,k)  / Ice%part_size(i,j,k)
        Ice%t_surf(i,j,k) = Ice%t_surf(i,j,k) / Ice%part_size(i,j,k)

        call ice5lay_resize(IST%h_snow(i,j,k), IST%t_snow(i,j,k), IST%h_ice(i,j,k), &
                            IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2),        &
                            IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4), 0.0,   &
                            -tot_heat/Ice%part_size(i,j,k), 0.0, 0.0, 0.0, &
                            -MU_TS*Ice%s_surf(i,j),                        &
                            heat_to_ocn, h2o_to_ocn, h2o_from_ocn, sn2ic)
      endif

      ! Convert constraining heat from energy (J/m^2) to flux (W/m^2)
      IST%qflx_lim_ice(i,j) = heat_limit_ice / dt_slow
      IST%qflx_res_ice(i,j) = heat_res_ice / dt_slow
      !
      ! Check for energy conservation
      !
      if (slab_ice) then
        e2m(2) = e2m(2) - IST%h_ice(i,j,2)*DI*LI
      else
        do k=2,km
          if (Ice%part_size(i,j,k)>0.0.and.IST%h_ice(i,j,k)>0.0) &
            e2m(k) = e2m(k)-e_to_melt(IST%h_snow(i,j,k), IST%t_snow(i,j,k),  &
                     IST%h_ice(i,j,k), IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2), &
                     IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4)) * Ice%part_size(i,j,k)
        enddo
      endif
      ! if (abs(sum(e2m) - heat_res_ice - heat_limit_ice)>DI*LI*1e-3) &
      !       print *, 'QFLUX conservation error at', i, j, 'heat2ice=',  &
      !             tot_heat, 'melted=', sum(e2m), 'h*part_size=', &
      !             IST%h_ice(i,j,:)*Ice%part_size(i,j,:)

    enddo ; enddo
  endif
  call mpp_clock_end(iceClock6)
  !
  ! Salt fluxes to ocean
  !
  call mpp_clock_begin(iceClock8)
  !### Inlining sensibly would change answers.
  hi_change(:,:)  = all_avg(IST%h_ice(isc:iec,jsc:jec,:), Ice%part_size(isc:iec,jsc:jec,:))-hi_change(:,:)
  Ice%flux_salt(:,:) = ice_bulk_salin*DI*hi_change(:,:)/dt_slow
  if (id_saltf>0)  sent = send_data(id_saltf, Ice%flux_salt(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)

  h2o_change(:,:) = all_avg(DS*IST%h_snow(isc:iec,jsc:jec,:)+DI*IST%h_ice(isc:iec,jsc:jec,:), &
               Ice%part_size(isc:iec,jsc:jec,:))-h2o_change(:,:)

  if (id_lsnk>0) then
    tmp2d(:,:) = h2o_change(:,:)*864e2*365/dt_slow
    do j=jsc,jec ; do i=isc,iec
      if (tmp2d(i,j)>0.0) tmp2d(i,j) = 0.0
    enddo ; enddo
    sent = send_data(id_lsnk, tmp2d(:,:), Ice%Time, mask=Ice%mask)
  endif
  if (id_lsrc>0) then
    tmp2d(:,:) = h2o_change(:,:)*864e2*365/dt_slow
    do j=jsc,jec ; do i=isc,iec
      if (tmp2d(i,j)<0.0) tmp2d(i,j) = 0.0
    enddo ; enddo
    sent = send_data(id_lsrc,  tmp2d(:,:), Ice%Time, mask=Ice%mask)
  endif

  call enable_SIS_averaging(dt_slow, Ice%Time, IST%diag)
  if (id_bsnk>0)  sent = send_data(id_bsnk, bsnk(isc:iec,jsc:jec)*864e2*365/dt_slow, Ice%Time, mask=Ice%mask)
  if (id_tmelt>0) call post_avg(id_tmelt, IST%tmelt, Ice%part_size(isc:iec,jsc:jec,2:), IST%diag, &
                                scale=1.0/dt_slow, mask=Ice%mask, wtd=.true.)
  if (id_bmelt>0) call post_avg(id_bmelt, IST%bmelt, Ice%part_size(isc:iec,jsc:jec,2:), IST%diag, &
                                scale=1.0/dt_slow, mask=Ice%mask, wtd=.true.)
  if (id_sn2ic>0) call post_avg(id_sn2ic, snow_to_ice, Ice%part_size(isc:iec,jsc:jec,2:), IST%diag, scale=1.0/dt_slow, mask=Ice%mask)
  if (id_qflim>0) sent = send_data(id_qflim, IST%qflx_lim_ice(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
  if (id_qfres>0) sent = send_data(id_qfres, IST%qflx_res_ice(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
  !
  ! Ice transport ... all ocean fluxes have been calculated by now
  !
  h2o_change(isc:iec,jsc:jec) = all_avg(DS*IST%h_snow(isc:iec,jsc:jec,:) &
                                       +DI*IST%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))

  ! Convert the velocities to C-grid points for transport.
  uc(:,:) = 0.0; vc(:,:) = 0.0
  do j=jsc,jec ; do I=isc-1,iec
    uc(I,j) = 0.5 * ( IST%u_ice(i,j-1) + IST%u_ice(i,j) )
  enddo ; enddo
  do J=jsc-1,jec ; do i = isc,iec
    vc(i,J) = 0.5 * ( IST%v_ice(i-1,j) + IST%v_ice(i,j) )
  enddo ; enddo

  call ice_transport(Ice%part_size, IST%h_ice, IST%h_snow, uc, vc, &
                     IST%t_ice, IST%t_snow, IST%sea_lev, hlim, dt_slow, &
                     Ice%G, IST%ice_transport_CSp)
  ! Set appropriate surface quantities in categories with no ice.
  do k=2,km ; do j=jsc,jec ; do i=isc,iec ; if (Ice%part_size(i,j,k)<1e-10) &
    Ice%t_surf(i,j,k) = Tfreeze-MU_TS*Ice%s_surf(i,j)
  enddo ; enddo ; enddo


!  tmp2d(:,:) = all_avg(DS*IST%h_snow(isc:iec,jsc:jec,:)+DI*IST%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))
!  do j=jsc,jec ; do i=isc,iec
!    Ice%mi(i,j) = tmp2d(i,j) 
!  enddo ; enddo
  ! Convert thickness and concentration to mass.
  Ice%mi(:,:) = 0.0 
  do k=2,km ; do j=jsc,jec ; do i=isc,iec
    Ice%mi(i,j) = Ice%mi(i,j) + (DS*IST%h_snow(i,j,k) + DI*IST%h_ice(i,j,k)) * Ice%part_size(i,j,k)
  enddo ; enddo ; enddo
  do j=jsc,jec ; do i=isc,iec ; mass(i,j) = Ice%mi(i,j) ; enddo ; enddo
  if (id_mi>0) sent = send_data(id_mi, Ice%mi(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)

  call disable_SIS_averaging(IST%diag)

  if (do_icebergs) call icebergs_incr_mass(Ice%icebergs, mass(isc:iec,jsc:jec)) ! Add icebergs mass in kg/m^2
  if (id_mib>0) sent = send_data(id_mib, mass(:,:), Ice%Time, mask = Ice%mask) ! Diagnose total mass
  if (id_slp>0) sent = send_data(id_slp, Ice%p_surf(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
  do j=jsc,jec ; do i=isc,iec
    if (slp2ocean) then
      Ice%p_surf(i,j) = Ice%p_surf(i,j)-1e5
    else
      Ice%p_surf(i,j) = 0.0
    endif
    Ice%p_surf(i,j) = Ice%p_surf(i,j) + grav*mass(i,j)
  enddo ; enddo
  h2o_change(:,:) = mass(:,:)-h2o_change(:,:)

  if (spec_ice) then                  ! over-write changes with specifications
    call get_sea_surface(Ice%Time, Ice%t_surf(:,:,1), Ice%part_size(isc:iec,jsc:jec,:), IST%h_ice (isc:iec,jsc:jec,2))
    call mpp_update_domains(Ice%part_size, Domain)
  endif

  Ice%part_size(:,:,1) = 1.0

  IST%part_size_uv(:,:,1) = 1.0

  do k=2,km
    Ice%part_size(:,:,1) = Ice%part_size(:,:,1) - Ice%part_size(:,:,k)

    call t_to_uv(Ice%part_size(:,:,k),IST%part_size_uv(:,:,k), Ice%G)
    IST%part_size_uv(:,:,1) = IST%part_size_uv(:,:,1) - IST%part_size_uv(:,:,k)
  enddo
  call mpp_clock_end(iceClock8)
  !
  ! Thermodynamics diagnostics
  !
  call mpp_clock_begin(iceClock9)
  if (id_cn>0) sent = send_data(id_cn, Ice%part_size(isc:iec,jsc:jec,2:km), Ice%Time, mask=spread(Ice%mask,3,km-1)       )

  ! TK Mod: 10/18/02
  !  if (id_obs_cn>0) sent = send_data(id_obs_cn, Obs_cn_ice(:,:,2), Ice%Time, &
       !                                       mask=Ice%mask       )

  if (id_ext>0) then
    do j=jsc,jec ; do i=isc,iec
      diagVar(i,j) = 0.0 ; if (Ice%part_size(i,j,1) < 0.85) diagVar(i,j) = 1.0
    enddo ; enddo
    sent = send_data(id_ext, diagVar(isc:iec,jsc:jec), Ice%Time, mask=Ice%mask)
  endif
  if (id_hs>0) sent  = send_data(id_hs, ice_avg(IST%h_snow(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)), &
                       Ice%Time, mask=Ice%mask)
  if (id_hi>0) sent  = send_data(id_hi, ice_avg(IST%h_ice(isc:iec,jsc:jec,:) ,Ice%part_size(isc:iec,jsc:jec,:)), &
                       Ice%Time, mask=Ice%mask)

  if (id_hio>0) sent  = send_data(id_hio, ice_avg(IST%h_ice(isc:iec,jsc:jec,:) ,Ice%part_size(isc:iec,jsc:jec,:)), &
                       Ice%Time, mask=Ice%mask,is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)

  ! TK Mod: 10/18/02: (commented out...does not compile yet... add later
  !  if (id_obs_hi>0) &
       !    sent = send_data(id_obs_hi, ice_avg(Obs_h_ice,Ice%part_size), Ice%Time, &
  !                     mask=Ice%mask)

  if (id_ts>0) sent = send_data(id_ts, ice_avg(Ice%t_surf-Tfreeze,Ice%part_size(isc:iec,jsc:jec,:)), Ice%Time, mask=Ice%mask)
  if (id_tsn>0) sent = send_data(id_tsn, ice_avg(IST%t_snow(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)), &
               Ice%Time, mask=Ice%mask)
  if (id_t1>0) sent = send_data(id_t1, ice_avg(IST%t_ice(isc:iec,jsc:jec,:,1),Ice%part_size(isc:iec,jsc:jec,:)), &
               Ice%Time, mask=Ice%mask)
  if (id_t2>0) sent = send_data(id_t2, ice_avg(IST%t_ice(isc:iec,jsc:jec,:,2),Ice%part_size(isc:iec,jsc:jec,:)), &
               Ice%Time, mask=Ice%mask)
  if (id_t3>0) sent = send_data(id_t3, ice_avg(IST%t_ice(isc:iec,jsc:jec,:,3),Ice%part_size(isc:iec,jsc:jec,:)), &
               Ice%Time, mask=Ice%mask)
  if (id_t4>0) sent = send_data(id_t4, ice_avg(IST%t_ice(isc:iec,jsc:jec,:,4),Ice%part_size(isc:iec,jsc:jec,:)), &
               Ice%Time, mask=Ice%mask)
  if (id_xprt>0) sent = send_data(id_xprt,  h2o_change(isc:iec,jsc:jec)*864e2*365/dt_slow, Ice%Time, mask=Ice%mask)
  if (id_e2m>0) then
    tmp2d(:,:) = 0.0
    do k=2,km ; do j=jsc,jec ; do i=isc,iec
      if (IST%h_ice(i,j,k)>0.0) &
        tmp2d(i,j) = tmp2d(i,j) + Ice%part_size(i,j,k)*e_to_melt(IST%h_snow(i,j,k), &
                     IST%t_snow(i,j,k), IST%h_ice(i,j,k), IST%t_ice(i,j,k,1),    &
                     IST%t_ice(i,j,k,2), IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4)    )
    enddo ; enddo ; enddo
    sent = send_data(id_e2m,  tmp2d(:,:), Ice%Time, mask=Ice%mask)
  endif

  call ice_top_to_ice_bottom(Ice, Ice%Ice_state, part_save, part_save_uv, G)
  !
  ! conservation checks:  bottom fluxes and final
  !
  if (conservation_check) then
    h2o(3)  = h2o(3) + dt_slow*sum(cell_area *(Ice%lprec+Ice%fprec-Ice%flux_q+Ice%runoff+Ice%calving),&
              cell_area > 0 )
    heat(3) = heat(3) + tot_frazil + dt_slow*sum(cell_area*(Ice%flux_sw_vis_dir+Ice%flux_sw_vis_dif+  &
              Ice%flux_sw_nir_dir+Ice%flux_sw_nir_dif+Ice%flux_lw- &
              Ice%flux_t-Ice%flux_lh -LI*(Ice%fprec+Ice%calving)), cell_area > 0 )
    salt(3)  = salt(3) + dt_slow*sum(cell_area*(Ice%flux_salt))/ice_bulk_salin !Kg salt added since start / earth_area

    h2o(4)  = sum(cell_area*all_avg(DI*IST%h_ice(isc:iec,jsc:jec,:)+ &
              DS*IST%h_snow(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:)))

    salt(4) = sum(cell_area*all_avg(DI*IST%h_ice(isc:iec,jsc:jec,:),Ice%part_size(isc:iec,jsc:jec,:))) !*ice_bulk_salin

    heat(4) = 0.0

    do k=2,km ; do j=jsc,jec ;  do i=isc,iec
      if ((Ice%part_size(i,j,k)>0.0) .and. (IST%h_ice(i,j,k)>0.0)) then
        if (slab_ice) then
          heat(4) = heat(4) - cell_area(i,j) * Ice%part_size(i,j,k)*IST%h_ice(i,j,2)*DI*LI
        else
          heat(4) = heat(4) - cell_area(i,j) * Ice%part_size(i,j,k)*e_to_melt(IST%h_snow(i,j,k), &
                    IST%t_snow(i,j,k), IST%h_ice(i,j,k), IST%t_ice(i,j,k,1), IST%t_ice(i,j,k,2),   &
                                                         IST%t_ice(i,j,k,3), IST%t_ice(i,j,k,4)    )
        endif
      endif
    enddo ;  enddo ; enddo
  endif

  call get_date(Ice%Time, iyr, imon, iday, ihr, imin, isec)
  call get_time(Ice%Time-set_date(iyr,1,1,0,0,0),isec,iday)
  if(verbose)  call ice_line(iyr, iday+1, isec, Ice%part_size(isc:iec,jsc:jec,:), &
             Ice%t_surf(:,:,1)-Tfreeze, Ice%G) 

  ! ### ADD BETTER ERROR HANDLING.
  do j=jsc,jec ; do i=isc,iec 
    if (Ice%mask(i,j).and.(abs(sum(Ice%part_size(i,j,:))-1.0)>1e-2)) &
         print *,'ICE%PART_SIZE=',Ice%part_size(i,j,:), 'DOES NOT SUM TO 1 AT', &
         Ice%G%geoLonT(i,j), Ice%G%geoLatT(i,j)
  enddo ; enddo
  call mpp_clock_end(iceClock9)


end subroutine update_ice_model_slow

! =====================================================================

subroutine post_avg_all(id, val, part, diag, mask, scale, offset, wtd)
  integer, intent(in) :: id
  real, dimension(:,:,:), intent(in) :: val, part
  type(SIS_diag_ctrl),  intent(in) :: diag
  logical, dimension(:,:), optional, intent(in) :: mask
  real,                    optional, intent(in) :: scale, offset
  logical,                 optional, intent(in) :: wtd
  ! This subroutine determines the average of a quantity across thickness
  ! categories and does a send data on it.

  real :: avg(size(val,1),size(val,2)), wts(size(val,1),size(val,2))
  real :: scl, off
  logical :: do_wt
  integer :: i, j, k, ni, nj, nk

  ni = size(val,1) ; nj = size(val,2) ; nk = size(val,3)
  if (size(part,1) /= ni) call SIS_error(FATAL, &
    "Mismatched i-sizes in post_avg.")
  if (size(part,2) /= nj) call SIS_error(FATAL, &
    "Mismatched j-sizes in post_avg.")
  if (size(part,3) /= nk) call SIS_error(FATAL, &
    "Mismatched k-sizes in post_avg.")

  scl = 1.0 ; if (present(scale)) scl = scale
  off = 0.0 ; if (present(offset)) off = offset
  do_wt = .false. ; if (present(wtd)) do_wt = wtd

  if (do_wt) then
    avg(:,:) = 0.0 ; wts(:,:) = 0.0
    do k=1,nk ; do j=1,nj ; do i=1,ni
      avg(i,j) = avg(i,j) + part(i,j,k)*(scl*val(i,j,k) + off)
      wts(i,j) = wts(i,j) + part(i,j,k)
    enddo ; enddo ; enddo
    do j=1,nj ; do i=1,ni
      if (wts(i,j) > 0.) then
        avg(i,j) = avg(i,j) / wts(i,j)
      else
        avg(i,j) = 0.0
      endif
    enddo ; enddo
  else
    avg(:,:) = 0.0
    do k=1,nk ; do j=1,nj ; do i=1,ni
      avg(i,j) = avg(i,j) + part(i,j,k)*(scl*val(i,j,k) + off)
    enddo ; enddo ; enddo
  endif

  call post_SIS_data(id, avg, diag, mask=mask)

end subroutine post_avg_all

subroutine post_avg_G(id, G, val, part, diag, mask, scale, offset, wtd)
  integer, intent(in) :: id
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:G%CatIce), intent(in) :: val, part
  type(SIS_diag_ctrl),  intent(in) :: diag
  logical, dimension(:,:), optional, intent(in) :: mask
  real,                    optional, intent(in) :: scale, offset
  logical,                 optional, intent(in) :: wtd
  ! This subroutine determines the average of a quantity across thickness
  ! categories and does a send data on it.

  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: avg, wts
  real :: scl, off
  logical :: do_wt
  integer :: i, j, k, isc, iec, jsc, jec, ncat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = G%CatIce

  scl = 1.0 ; if (present(scale)) scl = scale
  off = 0.0 ; if (present(offset)) off = offset
  do_wt = .false. ; if (present(wtd)) do_wt = wtd

  if (do_wt) then
    avg(:,:) = 0.0 ; wts(:,:) = 0.0
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      avg(i,j) = avg(i,j) + part(i,j,k)*(scl*val(i,j,k) + off)
      wts(i,j) = wts(i,j) + part(i,j,k)
    enddo ; enddo ; enddo
    do j=jsc,jec ; do i=isc,iec
      if (wts(i,j) > 0.) then
        avg(i,j) = avg(i,j) / wts(i,j)
      else
        avg(i,j) = 0.0
      endif
    enddo ; enddo
  else
    avg(:,:) = 0.0
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      avg(i,j) = avg(i,j) + part(i,j,k)*(scl*val(i,j,k) + off)
    enddo ; enddo ; enddo
  endif

  call post_SIS_data(id, avg, diag, mask=mask)

end subroutine post_avg_G

function is_NaN(x)
  real, intent(in) :: x
  logical :: is_nan
! This subroutine returns .true. if x is a NaN, and .false. otherwise.

  is_nan = (((x < 0.0) .and. (x >= 0.0)) .or. &
            (.not.(x < 0.0) .and. .not.(x >= 0.0)))

end function is_nan

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

end module ice_model_mod
