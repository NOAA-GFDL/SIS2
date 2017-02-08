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
!   This module handles the rapid thermodynamic interactions between the ice   !
! and the atmosphere, including heating and the accumulation of fluxes, but    !
! not changes to the ice or snow mass.                                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_fast_thermo

use SIS_diag_mediator, only : SIS_diag_ctrl
! ! use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
! ! use SIS_diag_mediator, only : query_SIS_averaging_enabled, post_SIS_data
! ! use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field

use SIS_debugging,     only : hchksum
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
! use MOM_hor_index, only : hor_index_type, hor_index_init
! use MOM_obsolete_params, only : obsolete_logical
! use MOM_string_functions, only : uppercase
use MOM_time_manager, only : time_type, time_type_to_real
use MOM_time_manager, only : set_date, set_time, operator(+), operator(-)
use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)

use coupler_types_mod, only : coupler_3d_bc_type

use SIS_types, only : ice_state_type, IST_chksum, IST_bounds_check
use SIS_types, only : fast_ice_avg_type, ice_rad_type, simple_OSS_type, total_sfc_flux_type
use ice_boundary_types, only : atmos_ice_boundary_type ! , land_ice_boundary_type
use SIS_hor_grid, only : SIS_hor_grid_type

use ice_grid, only : ice_grid_type

use SIS2_ice_thm,  only : SIS2_ice_thm_CS, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,  only : ice_temp_SIS2
use SIS2_ice_thm,  only : get_SIS2_thermo_coefs, enth_from_TS, Temp_from_En_S

implicit none ; private

public :: do_update_ice_model_fast, SIS_fast_thermo_init, SIS_fast_thermo_end
public :: fast_thermo_CS, avg_top_quantities, total_top_quantities, infill_array

type fast_thermo_CS ; private
  ! These two arrarys are used with column_check when evaluating the enthalpy
  ! conservation with the fast thermodynamics code.
  real, pointer, dimension(:,:,:) :: &
    enth_prev, heat_in

  logical :: debug        ! If true, write verbose checksums for debugging purposes.
  logical :: column_check ! If true, enable the heat check column by column.
  real    :: imb_tol      ! The tolerance for imbalances to be flagged by
                          ! column_check, nondim.
  logical :: bounds_check ! If true, check for sensible values of thicknesses
                          ! temperatures, fluxes, etc.

  integer :: n_fast = 0   ! The number of times update_ice_model_fast
                          ! has been called.

  ! These are pointers to the control structures for subsidiary modules.
  type(SIS2_ice_thm_CS), pointer  :: ice_thm_CSp => NULL()
end type fast_thermo_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> sum_top_quantities does a running sum of fluxes for later use by the slow ice
!!   physics and the ocean.  Nothing here will be exposed to other modules until
!!   after it has passed through avg_top_quantities.
subroutine sum_top_quantities (FIA, ABT, flux_u, flux_v, flux_t, flux_q, &
       flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif, &
       flux_lw, lprec, fprec, flux_lh, t_skin, dhdt, dedt, &
       dlwdt, SST, G, IG)
  type(fast_ice_avg_type),       intent(inout) :: FIA
  type(atmos_ice_boundary_type), intent(in)    :: ABT
  type(SIS_hor_grid_type),       intent(in)    :: G
  type(ice_grid_type),           intent(in)    :: IG
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), intent(in) :: &
    flux_u, flux_v, flux_t, flux_q, &
    flux_lw, &  ! The net longwave heat flux from the atmosphere into the
                ! ice or ocean, in W m-2.
    lprec, fprec, flux_lh, &
    flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif, &
    dhdt, dedt, &
    dlwdt       ! The partial derivative of the longwave heat flux from the
                ! atmosphere into the ice or ocean with ice skin temperature,
                ! in W m-2 K-1.
  real, dimension(G%isd:G%ied,G%jsd:G%jed,IG%CatIce), intent(in) :: &
    t_skin
  real, dimension(G%isd:G%ied,G%jsd:G%jed), intent(in) :: &
    SST

  real :: t_sfc
  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, ncat
  integer :: ind

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  i_off = LBOUND(ABT%t_flux,1) - G%isc
  j_off = LBOUND(ABT%t_flux,2) - G%jsc

  if (FIA%num_tr_fluxes < 0) then
    ! Determine how many atmospheric boundary fluxes have been passed in, and
    ! set up both an indexing array for these and a space to take their average.
    ! This code is only exercised the first time that sum_top_quantities is called.
    FIA%num_tr_fluxes = 0
    if (ABT%fluxes%num_bcs > 0) then
      do n=1,ABT%fluxes%num_bcs
        FIA%num_tr_fluxes = FIA%num_tr_fluxes + ABT%fluxes%bc(n)%num_fields
      enddo

      if (FIA%num_tr_fluxes > 0) then
        allocate(FIA%tr_flux_top(G%isd:G%ied, G%jsd:G%jed, 0:IG%CatIce, FIA%num_tr_fluxes))
        FIA%tr_flux_top(:,:,:,:) = 0.0
      endif
    endif
  endif

  if (FIA%avg_count == 0) then
    ! zero_top_quantities - zero fluxes to begin summing in ice fast physics.
    FIA%flux_u_top(:,:,:) = 0.0 ; FIA%flux_v_top(:,:,:) = 0.0
    FIA%flux_t_top(:,:,:) = 0.0 ; FIA%flux_q_top(:,:,:) = 0.0
    FIA%flux_lw_top(:,:,:) = 0.0 ; FIA%flux_lh_top(:,:,:) = 0.0
    FIA%flux_sw_nir_dir_top(:,:,:) = 0.0 ; FIA%flux_sw_nir_dif_top(:,:,:) = 0.0
    FIA%flux_sw_vis_dir_top(:,:,:) = 0.0 ; FIA%flux_sw_vis_dif_top(:,:,:) = 0.0
    FIA%lprec_top(:,:,:) = 0.0 ; FIA%fprec_top(:,:,:) = 0.0
    FIA%flux_sw_dn(:,:) = 0.0 ; FIA%Tskin_avg(:,:) = 0.0

    if (allocated(FIA%flux_t0)) then
      FIA%dhdt(:,:,:) = 0.0 ; FIA%dedt(:,:,:) = 0.0 ; FIA%dlwdt(:,:,:) = 0.0
      FIA%flux_t0(:,:,:) = 0.0 ; FIA%flux_q0(:,:,:) = 0.0
      FIA%flux_lw0(:,:,:) = 0.0 ; FIA%Tskin_cat(:,:,:) = 0.0
    endif

    if (FIA%num_tr_fluxes > 0) FIA%tr_flux_top(:,:,:,:) = 0.0
  endif

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,flux_u,flux_v,flux_t, &
!$OMP                                  flux_q,flux_sw_nir_dir,flux_sw_nir_dif,        &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,flux_lw,       &
!$OMP                                  lprec,fprec,flux_lh,FIA)
  do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
    FIA%flux_u_top(i,j,k)  = FIA%flux_u_top(i,j,k)  + flux_u(i,j,k)
    FIA%flux_v_top(i,j,k)  = FIA%flux_v_top(i,j,k)  + flux_v(i,j,k)
    FIA%flux_t_top(i,j,k)  = FIA%flux_t_top(i,j,k)  + flux_t(i,j,k)
    FIA%flux_q_top(i,j,k)  = FIA%flux_q_top(i,j,k)  + flux_q(i,j,k)
    FIA%flux_sw_nir_dir_top(i,j,k) = FIA%flux_sw_nir_dir_top(i,j,k) + flux_sw_nir_dir(i,j,k)
    FIA%flux_sw_nir_dif_top(i,j,k) = FIA%flux_sw_nir_dif_top(i,j,k) + flux_sw_nir_dif(i,j,k)
    FIA%flux_sw_vis_dir_top(i,j,k) = FIA%flux_sw_vis_dir_top(i,j,k) + flux_sw_vis_dir(i,j,k)
    FIA%flux_sw_vis_dif_top(i,j,k) = FIA%flux_sw_vis_dif_top(i,j,k) + flux_sw_vis_dif(i,j,k)
    FIA%flux_lw_top(i,j,k) = FIA%flux_lw_top(i,j,k) + flux_lw(i,j,k)
    FIA%lprec_top(i,j,k)   = FIA%lprec_top(i,j,k)   + lprec(i,j,k)
    FIA%fprec_top(i,j,k)   = FIA%fprec_top(i,j,k)   + fprec(i,j,k)
    FIA%flux_lh_top(i,j,k) = FIA%flux_lh_top(i,j,k) + flux_lh(i,j,k)
  enddo ; enddo ; enddo
  ! FIA%flux_sw_dn is accumulated where the fast radiation diagnostics are output
  ! because it depends on arrays that are stored in the public ice_data_type.

  ind = 0
  do n=1,ABT%fluxes%num_bcs ; do m=1,ABT%fluxes%bc(n)%num_fields
    ind = ind + 1
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      FIA%tr_flux_top(i,j,k,ind) = FIA%tr_flux_top(i,j,k,ind) + &
            ABT%fluxes%bc(n)%field(m)%values(i2,j2,k2)
    enddo ; enddo ; enddo
  enddo ; enddo

  if (allocated(FIA%flux_t0)) then
    !$OMP parallel do default(shared) private(t_sfc)
    do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
      t_sfc = SST(i,j) ; if (k>0) t_sfc = t_skin(i,j,k)
      FIA%dhdt(i,j,k) = FIA%dhdt(i,j,k) + dhdt(i,j,k)
      FIA%dedt(i,j,k) = FIA%dedt(i,j,k) + dedt(i,j,k)
      FIA%dlwdt(i,j,k) = FIA%dlwdt(i,j,k) + dlwdt(i,j,k)
      FIA%flux_t0(i,j,k) = FIA%flux_t0(i,j,k) + (flux_t(i,j,k) - t_sfc*dhdt(i,j,k))
      FIA%flux_q0(i,j,k) = FIA%flux_q0(i,j,k) + (flux_q(i,j,k) - t_sfc*dedt(i,j,k))
      FIA%flux_lw0(i,j,k) = FIA%flux_lw0(i,j,k) + (flux_lw(i,j,k) - t_sfc*dlwdt(i,j,k))
      FIA%Tskin_cat(i,j,k) = FIA%Tskin_cat(i,j,k) + t_sfc
    enddo ; enddo ; enddo
  endif

  FIA%avg_count = FIA%avg_count + 1

end subroutine sum_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> avg_top_quantities determines time average fluxes for later use by the
!!    slow ice physics and by the ocean.
subroutine avg_top_quantities(FIA, Rad, IST, G, IG)
  type(fast_ice_avg_type), intent(inout) :: FIA
  type(ice_rad_type),      intent(in)    :: Rad
  type(ice_state_type),    intent(in)    :: IST
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(in)    :: IG

  real    :: u, v, divid, sign
  real    :: I_avc    ! The inverse of the number of contributions.
  real    :: I_wts    ! 1.0 / ice_cover or 0 if ice_cover is 0, nondim.
  integer :: i, j, k, m, n, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  !
  ! compute average fluxes
  !
  if (FIA%avg_count == 0) call SIS_error(FATAL,'avg_top_quantities: '//&
       'no ocean model fluxes have been averaged')

  ! Rotate the stress from lat/lon to ocean coordinates and possibly change the
  ! sign to positive for downward fluxes of positive momentum.
  sign = 1.0 ; if (FIA%atmos_winds) sign = -1.0
  I_avc = 1.0/real(FIA%avg_count)

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,sign,I_avc,G,FIA,Rad) &
!$OMP                          private(u,v)
  do j=jsc,jec
    do k=0,ncat ; do i=isc,iec
      u = FIA%flux_u_top(i,j,k) * (sign*I_avc)
      v = FIA%flux_v_top(i,j,k) * (sign*I_avc)
      FIA%flux_u_top(i,j,k) = u*G%cos_rot(i,j)-v*G%sin_rot(i,j) ! rotate stress from lat/lon
      FIA%flux_v_top(i,j,k) = v*G%cos_rot(i,j)+u*G%sin_rot(i,j) ! to ocean coordinates
      FIA%flux_t_top(i,j,k)  = FIA%flux_t_top(i,j,k)  * I_avc
      FIA%flux_q_top(i,j,k)  = FIA%flux_q_top(i,j,k)  * I_avc
      FIA%flux_sw_nir_dir_top(i,j,k) = FIA%flux_sw_nir_dir_top(i,j,k) * I_avc
      FIA%flux_sw_nir_dif_top(i,j,k) = FIA%flux_sw_nir_dif_top(i,j,k) * I_avc
      FIA%flux_sw_vis_dir_top(i,j,k) = FIA%flux_sw_vis_dir_top(i,j,k) * I_avc
      FIA%flux_sw_vis_dif_top(i,j,k) = FIA%flux_sw_vis_dif_top(i,j,k) * I_avc
      FIA%flux_lw_top(i,j,k) = FIA%flux_lw_top(i,j,k) * I_avc
      FIA%fprec_top(i,j,k)   = FIA%fprec_top(i,j,k)   * I_avc
      FIA%lprec_top(i,j,k)   = FIA%lprec_top(i,j,k)   * I_avc
      FIA%flux_lh_top(i,j,k) = FIA%flux_lh_top(i,j,k) * I_avc

      ! Copy radiation fields from the fast to the slow states.
      if (k>0) FIA%sw_abs_ocn(i,j,k) = Rad%sw_abs_ocn(i,j,k)
      ! Convert frost forming atop sea ice into frozen precip.
      if ((k>0) .and. (FIA%flux_q_top(i,j,k) < 0.0)) then
        FIA%fprec_top(i,j,k) = FIA%fprec_top(i,j,k) - FIA%flux_q_top(i,j,k)
        FIA%flux_q_top(i,j,k) = 0.0
      endif
      do n=1,FIA%num_tr_fluxes
        FIA%tr_flux_top(i,j,k,n) = FIA%tr_flux_top(i,j,k,n) * I_avc
      enddo
    enddo ; enddo
    do i=isc,iec
      FIA%flux_sw_dn(i,j) = FIA%flux_sw_dn(i,j)*I_avc
      FIA%Tskin_avg(i,j) = FIA%Tskin_avg(i,j) * I_avc
    enddo
  enddo

  ! Determine the fractional ice coverage and the wind stresses averaged
  ! across all the ice thickness categories on an A-grid.
  FIA%WindStr_x(:,:) = 0.0 ; FIA%WindStr_y(:,:) = 0.0 ; FIA%ice_cover(:,:) = 0.0

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,FIA,Rad,IST) &
!$OMP                           private(I_wts)
  do j=jsc,jec
    do k=1,ncat ; do i=isc,iec
      FIA%WindStr_x(i,j) = FIA%WindStr_x(i,j) + IST%part_size(i,j,k) * FIA%flux_u_top(i,j,k)
      FIA%WindStr_y(i,j) = FIA%WindStr_y(i,j) + IST%part_size(i,j,k) * FIA%flux_v_top(i,j,k)
      FIA%ice_cover(i,j) = FIA%ice_cover(i,j) + IST%part_size(i,j,k)
    enddo ; enddo
    do i=isc,iec
      FIA%WindStr_ocn_x(i,j) = FIA%flux_u_top(i,j,0)
      FIA%WindStr_ocn_y(i,j) = FIA%flux_v_top(i,j,0)
      if (FIA%ice_cover(i,j) > 0.0) then
        I_wts = 1.0 / FIA%ice_cover(i,j)
        FIA%WindStr_x(i,j) = FIA%WindStr_x(i,j) * I_wts
        FIA%WindStr_y(i,j) = FIA%WindStr_y(i,j) * I_wts
        if (FIA%ice_cover(i,j) > 1.0) FIA%ice_cover(i,j) = 1.0

        ! The max with 0 in the following line is here for safety; the only known
        ! instance where it has been required is when reading a SIS-1-derived
        ! restart file with tiny negative concentrations. SIS2 should not need it.
        FIA%ice_free(i,j) = max(IST%part_size(i,j,0), 0.0)

    !    Rescale to add up to 1?
    !    I_wts = 1.0 / (FIA%ice_free(i,j) + FIA%ice_cover(i,j))
    !    FIA%ice_free(i,j) = FIA%ice_free(i,j) * I_wts
    !    FIA%ice_cover(i,j) = FIA%ice_cover(i,j) * I_wts
      else
        FIA%ice_free(i,j) = 1.0 ; FIA%ice_cover(i,j) = 0.0
        FIA%WindStr_x(i,j) = FIA%flux_u_top(i,j,0)
        FIA%WindStr_y(i,j) = FIA%flux_v_top(i,j,0)
      endif
    enddo
  enddo

  if (allocated(FIA%flux_t0)) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
      FIA%dhdt(i,j,k) = FIA%dhdt(i,j,k) * I_avc
      FIA%dedt(i,j,k) = FIA%dedt(i,j,k) * I_avc
      FIA%dlwdt(i,j,k) = FIA%dlwdt(i,j,k) * I_avc
      FIA%flux_t0(i,j,k) = FIA%flux_t0(i,j,k) * I_avc
      FIA%flux_q0(i,j,k) = FIA%flux_q0(i,j,k) * I_avc
      FIA%flux_lw0(i,j,k) = FIA%flux_lw0(i,j,k) * I_avc
      FIA%Tskin_cat(i,j,k) = FIA%Tskin_cat(i,j,k) * I_avc
    enddo ; enddo ; enddo

    ! Fill in the information to reconstruct the fluxes for any area-less categories.
    ! The open-ocean category must always be calculated for this to work properly.
    call infill_array(IST, FIA%dhdt(:,:,0), FIA%dhdt(:,:,1:), G, IG)
    call infill_array(IST, FIA%dedt(:,:,0), FIA%dedt(:,:,1:), G, IG)
    call infill_array(IST, FIA%dlwdt(:,:,0), FIA%dlwdt(:,:,1:), G, IG)
    call infill_array(IST, FIA%flux_t0(:,:,0), FIA%flux_t0(:,:,1:), G, IG)
    call infill_array(IST, FIA%flux_q0(:,:,0), FIA%flux_q0(:,:,1:), G, IG)
    call infill_array(IST, FIA%flux_lw0(:,:,0), FIA%flux_lw0(:,:,1:), G, IG)
    call infill_array(IST, FIA%Tskin_cat(:,:,0), FIA%Tskin_cat(:,:,1:), G, IG)
  endif

  ! set count to zero and fluxes will be zeroed before the next sum
  FIA%avg_count = 0

end subroutine avg_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> total_top_quantities determines the sum across partitions of various fluxes
!! for later use on a potentially different ice state on the slow side.
subroutine total_top_quantities(FIA, TSF, part_size, G, IG)
  type(fast_ice_avg_type),   intent(in)    :: FIA
  type(total_sfc_flux_type), intent(inout) :: TSF
  type(SIS_hor_grid_type),   intent(inout) :: G
  type(ice_grid_type),       intent(in)    :: IG
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
                             intent(in)    :: part_size

  integer :: i, j, k, m, n, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (TSF%num_tr_fluxes < 0) then
    ! Allocate the arrays to hold the tracer fluxes. This code is only exercised
    ! the first time that total_top_quantities is called.
    TSF%num_tr_fluxes = FIA%num_tr_fluxes
    if (TSF%num_tr_fluxes > 0) then
      allocate(TSF%tr_flux(G%isd:G%ied, G%jsd:G%jed, TSF%num_tr_fluxes))
    endif
  endif

  TSF%flux_u(:,:) = 0.0 ; TSF%flux_v(:,:) = 0.0
  TSF%flux_t(:,:) = 0.0 ; TSF%flux_q(:,:) = 0.0
  TSF%flux_sw_nir_dir(:,:) = 0.0 ; TSF%flux_sw_nir_dif(:,:) = 0.0
  TSF%flux_sw_vis_dir(:,:) = 0.0 ; TSF%flux_sw_vis_dif(:,:) = 0.0

  TSF%flux_lw(:,:) = 0.0 ; TSF%flux_lh(:,:) = 0.0
  TSF%fprec(:,:) = 0.0 ; TSF%lprec(:,:) = 0.0
  if (TSF%num_tr_fluxes > 0) TSF%tr_flux(:,:,:) = 0.0

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    TSF%flux_u(i,j) = TSF%flux_u(i,j) + part_size(i,j,k) * FIA%flux_u_top(i,j,k)
    TSF%flux_v(i,j) = TSF%flux_v(i,j) + part_size(i,j,k) * FIA%flux_v_top(i,j,k)
    TSF%flux_t(i,j) = TSF%flux_t(i,j) + part_size(i,j,k) * FIA%flux_t_top(i,j,k)
    TSF%flux_q(i,j) = TSF%flux_q(i,j) + part_size(i,j,k) * FIA%flux_q_top(i,j,k)
    TSF%flux_sw_nir_dir(i,j) = TSF%flux_sw_nir_dir(i,j) + &
                                part_size(i,j,k) * FIA%flux_sw_nir_dir_top(i,j,k)
    TSF%flux_sw_nir_dif(i,j) = TSF%flux_sw_nir_dif(i,j) + &
                                part_size(i,j,k) * FIA%flux_sw_nir_dif_top(i,j,k)
    TSF%flux_sw_vis_dir(i,j) = TSF%flux_sw_vis_dir(i,j) + &
                                part_size(i,j,k) * FIA%flux_sw_vis_dir_top(i,j,k)
    TSF%flux_sw_vis_dif(i,j) = TSF%flux_sw_vis_dif(i,j) + &
                                part_size(i,j,k) * FIA%flux_sw_vis_dif_top(i,j,k)

    TSF%flux_lw(i,j) = TSF%flux_lw(i,j) + part_size(i,j,k) * FIA%flux_lw_top(i,j,k)
    TSF%flux_lh(i,j) = TSF%flux_lh(i,j) + part_size(i,j,k) * FIA%flux_lh_top(i,j,k)
    TSF%fprec(i,j) = TSF%fprec(i,j) + part_size(i,j,k) * FIA%fprec_top(i,j,k)
    TSF%lprec(i,j) = TSF%lprec(i,j) + part_size(i,j,k) * FIA%lprec_top(i,j,k)

    do n=1,TSF%num_tr_fluxes
      TSF%tr_flux(i,j,n) = TSF%tr_flux(i,j,n) + part_size(i,j,k) * FIA%tr_flux_top(i,j,k,n)
    enddo
  enddo ; enddo ; enddo

  !   If the sum of part_size across all the ice and ocean categories is not
  ! exactly 1, rescaling might be advisable, but for now it is assumed that
  ! part_size is properly scaled.

end subroutine total_top_quantities

!> infill_array fills in an array with actual, interpolated or plausible values
!! from other ice categories.
subroutine infill_array(IST, ta_ocn, ta, G, IG)
  type(ice_state_type),    intent(in)    :: IST  !< The ice state type whose values of part_size and
                                                 !! mH_ice are being used to control the infilling.
  type(SIS_hor_grid_type), intent(in)    :: G    !< The sea-ice lateral grid type.
  type(ice_grid_type),     intent(in)    :: IG   !< The sea-ice grid type.
  real, dimension(G%isd:G%ied, G%jsd:G%jed), &
                           intent(in)    :: ta_ocn !< The value of the array for ocean partitions.
  real, dimension(G%isd:G%ied, G%jsd:G%jed, IG%CatIce), &
                           intent(inout) :: ta   !< The array that is being infilled.

  integer, dimension(G%isd:G%ied) :: k_thick, k_thin
  real, dimension(G%isd:G%ied) :: &
    ta_thick, ta_thin, mH_thin
  real :: wt2, mH_cat_tgt
  logical, dimension(G%isd:G%ied) :: infill
  logical :: any_ice
  integer :: i, j, k, k2, k3, isc, iec, jsc, jec, ncat, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  do j=jsc,jec
    ! Determine which categories exist.
    do i=isc,iec ; k_thick(i) = -1 ; k_thin(i) = ncat+1 ; mH_thin(i) = 0.0 ; enddo
    any_ice = .false.
    do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k) > 0.0) then
      k_thick(i) = k
      ta_thick(i) = ta(i,j,k)
      any_ice = .true.
    endif ; enddo ; enddo

    do i=isc,iec
      if (k_thick(i) < 1) ta_thick(i) = G%mask2dT(i,j) * ta_ocn(i,j)
    enddo

    if (.not.any_ice) then
      ! This entire row is ice free.
      do k=1,ncat ; do i=isc,iec
        ta(i,j,k) = G%mask2dT(i,j)*ta_ocn(i,j)
      enddo ; enddo
    else
      do k=ncat,1,-1 ; do i=isc,iec ; if (IST%part_size(i,j,k) > 0.0) then
        k_thin(i) = k
        ta_thin(i) = ta(i,j,k)
        mH_thin(i) = IST%mH_ice(i,j,k)
      endif ; enddo ; enddo

      ! The logic for the k=1 (thinest) category is particularly simple.
      do i=isc,iec ; if (IST%part_size(i,j,1) <= 0.0) then
        ta(i,j,1) = G%mask2dT(i,j)*ta_ocn(i,j)
      endif ; enddo
      do k=2,ncat ; do i=isc,iec
        if (G%mask2dT(i,j) > 0.0) then
          if (k > k_thick(i)) then
            ! This is the most common case, since it applies to ice-free water.
            ta(i,j,k) = ta_thick(i)
          elseif (IST%part_size(i,j,k) > 0.0) then
            ! This is a valid value - no interpolation is needed.
            cycle
          elseif (k < k_thin(i)) then
            ! Linearly interpolate to the ocean's freezing temperature.
            mH_cat_tgt = 0.5*(IG%mH_cat_bound(K) + IG%mH_cat_bound(K+1))
            wt2 = 0.0  ! If in doubt, use ta_thin.  This should not happen.
            if (mH_thin(i) > mH_cat_tgt) wt2 = (mH_thin(i) - mH_cat_tgt) / mH_thin(i)

            ta(i,j,k) = wt2*ta_ocn(i,j) + (1.0-wt2)*ta_thin(i)
          else
            ! Linearly interpolate between bracketing neighboring
            ! categories with valid values.
            do k2=k+1,ncat ; if (IST%part_size(i,j,k2) > 0.0) exit ; enddo
            do k3=k-1,1,-1 ; if (IST%part_size(i,j,k3) > 0.0) exit ; enddo
            mH_cat_tgt = 0.5*(IG%mH_cat_bound(K) + IG%mH_cat_bound(K+1))
            wt2 = 0.5
            if ((IST%mH_ice(i,j,k3) > mH_cat_tgt) .and. &
                (IST%mH_ice(i,j,k2) < mH_cat_tgt)) &
              wt2 = (IST%mH_ice(i,j,k3) - mH_cat_tgt) / &
                    (IST%mH_ice(i,j,k3) - IST%mH_ice(i,j,k2))

            ta(i,j,k) = wt2*ta(i,j,k2) + (1.0-wt2) * ta(i,j,k3)
          endif
        else ! This is a land point.
          ta(i,j,k) = 0.0
        endif
      enddo ; enddo

    endif ! .not.any_ice
  enddo

end subroutine infill_array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> do_update_ice_model_fast applies the surface heat fluxes, shortwave radiation
!!   diffusion of heat to the sea-ice to implicitly determine a new temperature
!!   profile, subject to the constraint that ice and snow temperatures are never
!!   above freezing.  Melting and freezing occur elsewhere.
subroutine do_update_ice_model_fast(Atmos_boundary, IST, sOSS, Rad, FIA, &
                                    Time_step, CS, G, IG )

  type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary
  type(ice_state_type),          intent(inout) :: IST
  type(simple_OSS_type),         intent(in)    :: sOSS
  type(ice_rad_type),            intent(inout) :: Rad
  type(fast_ice_avg_type),       intent(inout) :: FIA
  type(time_type),               intent(in)    :: Time_step  ! The amount of time over which to advance the ice.
  type(fast_thermo_CS),          pointer       :: CS
  type(SIS_hor_grid_type),       intent(inout) :: G
  type(ice_grid_type),           intent(in)    :: IG

  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce) :: &
    flux_t, flux_q, flux_lh, &
    flux_lw, &  ! The net longwave heat flux into the ice, in W m-2.
    flux_sw_nir_dir, flux_sw_nir_dif, &
    flux_sw_vis_dir, flux_sw_vis_dif, &
    flux_u, flux_v, lprec, fprec, &
    dhdt, &   ! The derivative of the upward sensible heat flux with the surface
              ! temperature in W m-2 K-1.
    dedt, &   ! The derivative of the sublimation rate with the surface
              ! temperature, in kg m-2 s-1 K-1.
    dlwdt     ! The derivative of the downward radiative heat flux with surface
              ! temperature (i.e. d(flux_lw)/d(surf_temp)) in W m-2 K-1.
  real, dimension(0:IG%NkIce) :: T_col ! The temperature of a column of ice and snow in degC.
  real, dimension(IG%NkIce)   :: S_col ! The thermodynamic salinity of a column of ice, in g/kg.
  real, dimension(0:IG%NkIce) :: enth_col   ! The enthalpy of a column of snow and ice, in enth_unit (J/kg?).
  real, dimension(0:IG%NkIce) :: SW_abs_col
  real :: dt_fast, ts_new, dts, hf, hfd, latent
  real :: hf_0    ! The positive upward surface heat flux when T_sfc = 0 C, in W m-2.
  real :: dhf_dt  ! The deriviative of the upward surface heat flux with Ts, in W m-2 C-1.
  real :: flux_sw ! sum over dir/dif vis/nir components
  real :: LatHtFus       ! The latent heat of fusion of ice in J/kg.
  real :: LatHtVap       ! The latent heat of vaporization of water at 0C in J/kg.
  real :: H_to_m_ice     ! The specific volumes of ice and snow times the
  real :: H_to_m_snow    ! conversion factor from thickness units, in m H-1.
  logical :: slab_ice    ! If true, use the very old slab ice thermodynamics,
                         ! with effectively zero heat capacity of ice and snow.
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

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_fast_thermo: Module must be initialized before it is used.")

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  i_off = LBOUND(Atmos_boundary%t_flux,1) - G%isc
  j_off = LBOUND(Atmos_boundary%t_flux,2) - G%jsc
  NkIce = IG%NkIce ; I_Nk = 1.0 / NkIce ; kg_H_Nk = IG%H_to_kg_m2 * I_Nk

  CS%n_fast = CS%n_fast + 1

  if (CS%column_check .and. (FIA%avg_count==0)) then
    CS%heat_in(:,:,:) = 0.0
    CS%enth_prev(:,:,:) = 0.0
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,k)>0.0) then
      CS%enth_prev(i,j,k) = (IST%mH_snow(i,j,k)*IG%H_to_kg_m2) * IST%enth_snow(i,j,k,1)
      do m=1,IG%NkIce
        CS%enth_prev(i,j,k) = CS%enth_prev(i,j,k) + &
                               (IST%mH_ice(i,j,k)*kg_H_Nk) * IST%enth_ice(i,j,k,m)
      enddo
    endif ; enddo ; enddo ; enddo
  endif

  if (CS%debug) &
    call IST_chksum("Start do_update_ice_model_fast", IST, G, IG)

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Atmos_boundary,i_off, &
!$OMP                                  j_off,flux_u,flux_v,flux_t,flux_q,flux_lw, &
!$OMP                                  flux_sw_nir_dir,flux_sw_nir_dif,               &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,               &
!$OMP                                  lprec,fprec,dhdt,dedt,dlwdt        )            &
!$OMP                           private(i2,j2,k2)
  do j=jsc,jec
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
      dlwdt(i,j,k) = -1.*Atmos_boundary%drdt(i2,j2,k2)
    enddo ; enddo
  enddo

  if (CS%debug) then
    call hchksum(flux_u(:,:,1:), "Mid do_fast flux_u", G%HI)
    call hchksum(flux_v(:,:,1:), "Mid do_fast flux_v", G%HI)
    call hchksum(flux_t(:,:,1:), "Mid do_fast flux_t", G%HI)
    call hchksum(flux_q(:,:,1:), "Mid do_fast flux_q", G%HI)
    call hchksum(flux_lw(:,:,1:), "Mid do_fast flux_lw", G%HI)
    call hchksum(flux_sw_nir_dir(:,:,1:), "Mid do_fast flux_sw_nir_dir", G%HI)
    call hchksum(flux_sw_nir_dif(:,:,1:), "Mid do_fast flux_sw_nir_dif", G%HI)
    call hchksum(flux_sw_vis_dir(:,:,1:), "Mid do_fast flux_sw_vis_dir", G%HI)
    call hchksum(flux_sw_vis_dif(:,:,1:), "Mid do_fast flux_sw_vis_dif", G%HI)
    call hchksum(lprec(:,:,1:), "Mid do_fast lprec", G%HI)
    call hchksum(fprec(:,:,1:), "Mid do_fast fprec", G%HI)
    call hchksum(dhdt(:,:,1:), "Mid do_fast dhdt", G%HI)
    call hchksum(dedt(:,:,1:), "Mid do_fast dedt", G%HI)
    call hchksum(dlwdt(:,:,1:), "Mid do_fast dlwdt", G%HI)
  endif

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             Latent_fusion=LatHtFus, Latent_vapor=LatHtVap, slab_ice=slab_ice)

  do j=jsc,jec ; do i=isc,iec
    flux_lh(i,j,0) = LatHtVap * flux_q(i,j,0)
  enddo ; enddo

  !
  ! implicit update of ice surface temperature
  !
  dt_fast = time_type_to_real(Time_step)

  enth_liq_0 = Enth_from_TS(0.0, 0.0, IST%ITV) ; I_enth_unit = 1.0 / enth_units

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,NkIce,IST,dhdt,dedt,dlwdt,   &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,flux_sw_nir_dir, &
!$OMP                                  flux_sw_nir_dif,flux_t,flux_q,flux_lw,enth_liq_0,&
!$OMP                                  dt_fast,flux_lh,I_enth_unit,G,S_col,kg_H_Nk,slab_ice,&
!$OMP                                  enth_units,LatHtFus,LatHtVap,IG,sOSS,FIA,Rad,CS) &
!$OMP                          private(latent,enth_col,flux_sw,dhf_dt,                  &
!$OMP                                  hf_0,ts_new,dts,SW_abs_col,SW_absorbed,enth_here,&
!$OMP                                  tot_heat_in,enth_imb,norm_enth_imb     )
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    if (IST%part_size(i,j,k) > 0.0) then
      enth_col(0) = IST%enth_snow(i,j,k,1)
      do m=1,NkIce ; enth_col(m) = IST%enth_ice(i,j,k,m) ; enddo

      ! In the case of sublimation of either snow or ice, the vapor is at 0 C.
      ! If the vapor should be at a different temperature, a correction would be
      ! made here.
      if (slab_ice) then
        latent = LatHtVap + LatHtFus
      elseif (IST%mH_snow(i,j,k)>0.0) then
        latent = LatHtVap + (enth_liq_0 - IST%enth_snow(i,j,k,1)) * I_enth_unit
      else
        latent = LatHtVap + (enth_liq_0 - IST%enth_ice(i,j,k,1)) * I_enth_unit
      endif
      flux_sw = (flux_sw_vis_dir(i,j,k) + flux_sw_vis_dif(i,j,k)) + &
                (flux_sw_nir_dir(i,j,k) + flux_sw_nir_dif(i,j,k))

      dhf_dt = (dhdt(i,j,k) + dedt(i,j,k)*latent) - dlwdt(i,j,k)
      hf_0 = ((flux_t(i,j,k) + flux_q(i,j,k)*latent) - &
              (flux_lw(i,j,k) + Rad%sw_abs_sfc(i,j,k)*flux_sw)) - &
             dhf_dt * Rad%t_skin(i,j,k)

      SW_abs_col(0) = Rad%sw_abs_snow(i,j,k)*flux_sw
      do m=1,NkIce ; SW_abs_col(m) = Rad%sw_abs_ice(i,j,k,m)*flux_sw ; enddo

      !   This call updates the snow and ice temperatures and accumulates the
      ! surface and bottom melting/freezing energy.  The ice and snow do not
      ! actually lose or gain any mass from freezing or melting.
      ! mw/new - pass melt pond (surface temp fixed at freezing when present)
      call ice_temp_SIS2(IST%mH_pond(i,j,k)*IG%H_to_kg_m2, &
                         IST%mH_snow(i,j,k)*IG%H_to_kg_m2, &
                         IST%mH_ice(i,j,k)*IG%H_to_kg_m2, &
                         enth_col, S_col, hf_0, dhf_dt, SW_abs_col, &
                         sOSS%T_fr_ocn(i,j), FIA%bheat(i,j), ts_new, &
                         dt_fast, NkIce, FIA%tmelt(i,j,k), FIA%bmelt(i,j,k), &
                         CS%ice_thm_CSp, IST%ITV, CS%column_check)
      IST%enth_snow(i,j,k,1) = enth_col(0)
      do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_col(m) ; enddo

      dts               = ts_new - Rad%t_skin(i,j,k)
      Rad%t_skin(i,j,k) = ts_new
      flux_t(i,j,k)  = flux_t(i,j,k)  + dts * dhdt(i,j,k)
      flux_q(i,j,k)  = flux_q(i,j,k)  + dts * dedt(i,j,k)
      flux_lw(i,j,k) = flux_lw(i,j,k) + dts * dlwdt(i,j,k)
      flux_lh(i,j,k) = latent * flux_q(i,j,k)

      if (CS%column_check) then
        SW_absorbed = SW_abs_col(0)
        do m=1,NkIce ; SW_absorbed = SW_absorbed + SW_abs_col(m) ; enddo
        CS%heat_in(i,j,k) = CS%heat_in(i,j,k) + dt_fast * &
          ((flux_lw(i,j,k) + Rad%sw_abs_sfc(i,j,k)*flux_sw) + SW_absorbed + &
           FIA%bheat(i,j) - (flux_t(i,j,k) + flux_lh(i,j,k)))

        enth_here = (IG%H_to_kg_m2*IST%mH_snow(i,j,k)) * enth_col(0)
        do m=1,NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k)*kg_H_Nk) * enth_col(m)
        enddo
        tot_heat_in = enth_units * (CS%heat_in(i,j,k) - &
                                    (FIA%bmelt(i,j,k) + FIA%tmelt(i,j,k)))
        enth_imb = enth_here - (CS%enth_prev(i,j,k) + tot_heat_in)
        if (abs(enth_imb) > CS%imb_tol * (abs(enth_here) + &
                  abs(CS%enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
          norm_enth_imb = enth_imb / (abs(enth_here) + &
                  abs(CS%enth_prev(i,j,k)) + abs(tot_heat_in))
          enth_imb = enth_here - (CS%enth_prev(i,j,k) + tot_heat_in)
        endif
      endif

    else ! IST%mH_ice <= 0
      flux_lh(i,j,k) = LatHtVap * flux_q(i,j,k)
    endif
  enddo ; enddo ; enddo

  call sum_top_quantities(FIA, Atmos_boundary, flux_u, flux_v, flux_t, &
    flux_q, flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif, &
    flux_lw, lprec, fprec, flux_lh, Rad%t_skin, dhdt, dedt, dlwdt, sOSS%SST_C, &
    G, IG )

  if (CS%debug) &
    call IST_chksum("End do_update_ice_model_fast", IST, G, IG)

  if (CS%bounds_check) &
    call IST_bounds_check(IST, G, IG, "End of update_ice_fast", Rad=Rad) !, OSS=sOSS)

end subroutine do_update_ice_model_fast

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_fast_thermo_init - initializes the parameters and diagnostics associated
!!    with the SIS_fast_thermo module.
subroutine SIS_fast_thermo_init(Time, G, IG, param_file, diag, CS)
  type(time_type),     target, intent(in)    :: Time   ! current time
  type(SIS_hor_grid_type),     intent(in)    :: G      ! The horizontal grid structure
  type(ice_grid_type),         intent(in)    :: IG     ! The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(fast_thermo_CS),        pointer       :: CS

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_fast_thermo" ! This module's name.

  call callTree_enter("SIS_fast_thermo_init(), SIS_fast_thermo.F90")

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_fast_thermo_init called with associated control structure.")
!    return
  else
    allocate(CS)
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, &
     "This module applies rapidly varying heat fluxes to the ice and does an "//&
     "implicit surface temperature calculation.")

  call get_param(param_file, mod, "COLUMN_CHECK", CS%column_check, &
                 "If true, add code to allow debugging of conservation \n"//&
                 "column-by-column.  This does not change answers, but \n"//&
                 "can increase model run time.", default=.false.)
  call get_param(param_file, mod, "IMBALANCE_TOLERANCE", CS%imb_tol, &
                 "The tolerance for imbalances to be flagged by COLUMN_CHECK.", &
                 units="nondim", default=1.0e-9)
  call get_param(param_file, mod, "ICE_BOUNDS_CHECK", CS%bounds_check, &
                 "If true, periodically check the values of ice and snow \n"//&
                 "temperatures and thicknesses to ensure that they are \n"//&
                 "sensible, and issue warnings if they are not.  This \n"//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)
  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)

  call SIS2_ice_thm_init(param_file, CS%ice_thm_CSp)

  if (CS%column_check) then
    allocate(CS%enth_prev(G%HI%isd:G%HI%ied, G%HI%jsd:G%HI%jed, IG%CatIce)) ; CS%enth_prev(:,:,:) = 0.0
    allocate(CS%heat_in(G%HI%isd:G%HI%ied, G%HI%jsd:G%HI%jed, IG%CatIce)) ; CS%heat_in(:,:,:) = 0.0
  endif

  call callTree_leave("SIS_fast_thermo_init()")

end subroutine SIS_fast_thermo_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_fast_thermo_end deallocates any memory associated with this module.
subroutine SIS_fast_thermo_end(CS)
  type(fast_thermo_CS), pointer :: CS

  call SIS2_ice_thm_end(CS%ice_thm_CSp)

  if (associated(CS)) deallocate(CS)

end subroutine SIS_fast_thermo_end

end module SIS_fast_thermo
