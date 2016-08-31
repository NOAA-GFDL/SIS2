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

use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
! ! use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field

use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
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

use ice_type_mod, only : ice_state_type, IST_chksum, IST_bounds_check
use ice_type_mod, only : fast_ice_avg_type, ocean_sfc_state_type
use ice_type_mod, only : atmos_ice_boundary_type, land_ice_boundary_type
use ice_type_mod, only : fast_thermo_CS
use ice_utils_mod, only : post_avg
use SIS_hor_grid, only : SIS_hor_grid_type

use ice_grid, only : ice_grid_type

use SIS2_ice_thm,  only : ice_temp_SIS2
use SIS2_ice_thm,  only : get_SIS2_thermo_coefs, enth_from_TS, Temp_from_En_S, T_freeze

implicit none ; private

#include <SIS2_memory.h>

public :: do_update_ice_model_fast, SIS_fast_thermo_init, SIS_fast_thermo_end
public :: avg_top_quantities

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> sum_top_quantities does a running sum of fluxes for later use by the slow ice
!!   physics and the ocean.  Nothing here will be exposed to other modules until
!!   after it has passed through avg_top_quantities.
subroutine sum_top_quantities (FIA, ABT, flux_u, flux_v, flux_t, flux_q, &
       flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif,&
       flux_lw, lprec, fprec, flux_lh, G, IG)
  type(fast_ice_avg_type),       intent(inout) :: FIA
  type(atmos_ice_boundary_type), intent(in)    :: ABT
  type(SIS_hor_grid_type),       intent(in)    :: G
  type(ice_grid_type),           intent(in)    :: IG
  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:IG%CatIce), intent(in) :: &
    flux_u, flux_v, flux_t, flux_q, flux_lw, lprec, fprec, flux_lh
  real, dimension(G%isc:G%iec,G%jsc:G%jec,0:IG%CatIce), intent(in) :: &
    flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif

  integer :: i, j, k, m, n, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, ncat
  integer :: ind, max_num_fields, next_index

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  i_off = LBOUND(ABT%t_flux,1) - G%isc
  j_off = LBOUND(ABT%t_flux,2) - G%jsc

  if (FIA%num_tr_fluxes < 0) then
    ! Determine how many atmospheric boundary fluxes have been passed in, and
    ! set up both an indexing array for these and a space to take their average.
    ! This code is only exercised the first time that sum_top_quantities is called.
    FIA%num_tr_fluxes = 0
    if (ABT%fluxes%num_bcs > 0) then
      max_num_fields = 0
      do n=1,ABT%fluxes%num_bcs
        FIA%num_tr_fluxes = FIA%num_tr_fluxes + ABT%fluxes%bc(n)%num_fields
        max_num_fields = max(max_num_fields, ABT%fluxes%bc(n)%num_fields)
      enddo

      if (FIA%num_tr_fluxes > 0) then
        allocate(FIA%tr_flux_top(SZI_(G), SZJ_(G), 0:IG%CatIce, FIA%num_tr_fluxes))
        FIA%tr_flux_top(:,:,:,:) = 0.0

        allocate(FIA%tr_flux_index(max_num_fields, ABT%fluxes%num_bcs))
        FIA%tr_flux_index(:,:) = -1 ; next_index = 1
        do n=1,ABT%fluxes%num_bcs ; do m=1,ABT%fluxes%bc(n)%num_fields
          FIA%tr_flux_index(m, n) = next_index ; next_index = next_index + 1
        enddo ; enddo
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

  do n=1,ABT%fluxes%num_bcs ; do m=1,ABT%fluxes%bc(n)%num_fields
    ind = FIA%tr_flux_index(m,n)
    if (ind < 1) call SIS_error(FATAL, "Bad boundary flux index in sum_top_quantities.")
    do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      FIA%tr_flux_top(i,j,k,ind) = FIA%tr_flux_top(i,j,k,ind) + &
            ABT%fluxes%bc(n)%field(m)%values(i2,j2,k2)
    enddo ; enddo ; enddo
  enddo ; enddo

  FIA%avg_count = FIA%avg_count + 1

end subroutine sum_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> avg_top_quantities determines time average fluxes for later use by the
!!    slow ice physics and by the ocean.
subroutine avg_top_quantities(FIA, part_size, G, IG)
  type(fast_ice_avg_type), intent(inout) :: FIA
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(in)    :: IG
  real, dimension(SZI_(G),SZJ_(G),0:IG%CatIce), &
                           intent(in)    :: part_size

  real    :: u, v, divid, sign
  real :: I_wts    ! 1.0 / wts or 0 if wts is 0, nondim.
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
  divid = 1.0/real(FIA%avg_count)

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,sign,divid,G,FIA) private(u,v)
  do j=jsc,jec
    do k=0,ncat ;  do i=isc,iec
      u = FIA%flux_u_top(i,j,k) * (sign*divid)
      v = FIA%flux_v_top(i,j,k) * (sign*divid)
      FIA%flux_u_top(i,j,k) = u*G%cos_rot(i,j)-v*G%sin_rot(i,j) ! rotate stress from lat/lon
      FIA%flux_v_top(i,j,k) = v*G%cos_rot(i,j)+u*G%sin_rot(i,j) ! to ocean coordinates
      FIA%flux_t_top(i,j,k)  = FIA%flux_t_top(i,j,k)  * divid
      FIA%flux_q_top(i,j,k)  = FIA%flux_q_top(i,j,k)  * divid
      FIA%flux_sw_nir_dir_top(i,j,k) = FIA%flux_sw_nir_dir_top(i,j,k) * divid
      FIA%flux_sw_nir_dif_top(i,j,k) = FIA%flux_sw_nir_dif_top(i,j,k) * divid
      FIA%flux_sw_vis_dir_top(i,j,k) = FIA%flux_sw_vis_dir_top(i,j,k) * divid
      FIA%flux_sw_vis_dif_top(i,j,k) = FIA%flux_sw_vis_dif_top(i,j,k) * divid
      FIA%flux_lw_top(i,j,k) = FIA%flux_lw_top(i,j,k) * divid
      FIA%fprec_top(i,j,k)   = FIA%fprec_top(i,j,k)   * divid
      FIA%lprec_top(i,j,k)   = FIA%lprec_top(i,j,k)   * divid
      FIA%flux_lh_top(i,j,k) = FIA%flux_lh_top(i,j,k) * divid
      ! Convert frost forming atop sea ice into frozen precip.
      if ((k>0) .and. (FIA%flux_q_top(i,j,k) < 0.0)) then
        FIA%fprec_top(i,j,k) = FIA%fprec_top(i,j,k) - FIA%flux_q_top(i,j,k)
        FIA%flux_q_top(i,j,k) = 0.0
      endif
      do n=1,FIA%num_tr_fluxes
        FIA%tr_flux_top(i,j,k,n) = FIA%tr_flux_top(i,j,k,n) * divid
      enddo
    enddo ; enddo
  enddo
  call pass_vector(FIA%flux_u_top, FIA%flux_v_top, G%Domain, stagger=AGRID)

  ! Determine the fractional ice coverage and the wind stresses averaged
  ! across all the ice thickness categories on an A-grid.  This is done
  ! over the entire data domain for safety.
  FIA%WindStr_x(:,:) = 0.0 ; FIA%WindStr_y(:,:) = 0.0 ; FIA%ice_cover(:,:) = 0.0
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,ncat,FIA,part_size) &
!$OMP                           private(I_wts)
  do j=jsd,jed
    do k=1,ncat ; do i=isd,ied
      FIA%WindStr_x(i,j) = FIA%WindStr_x(i,j) + part_size(i,j,k) * FIA%flux_u_top(i,j,k)
      FIA%WindStr_y(i,j) = FIA%WindStr_y(i,j) + part_size(i,j,k) * FIA%flux_v_top(i,j,k)
      FIA%ice_cover(i,j) = FIA%ice_cover(i,j) + part_size(i,j,k)
    enddo ; enddo
    do i=isd,ied
      if (FIA%ice_cover(i,j) > 0.0) then
        I_wts = 1.0 / FIA%ice_cover(i,j)
        FIA%WindStr_x(i,j) = FIA%WindStr_x(i,j) * I_wts
        FIA%WindStr_y(i,j) = FIA%WindStr_y(i,j) * I_wts
        if (FIA%ice_cover(i,j) > 1.0) FIA%ice_cover(i,j) = 1.0

        ! The max with 0 in the following line is here for safety; the only known
        ! instance where it has been required is when reading a SIS-1-derived
        ! restart file with tiny negative concentrations. SIS2 should not need it.
        FIA%ice_free(i,j) = max(part_size(i,j,0), 0.0)

    !    Rescale to add up to 1?
    !    I_wts = 1.0 / (FIA%ice_free(i,j) + FIA%ice_cover(i,j))
    !    FIA%ice_free(i,j) = FIA%ice_free(i,j) * I_wts
    !    FIA%ice_cover(i,j) = FIA%ice_cover(i,j) * I_wts
      else
        FIA%ice_free(i,j) = 1.0 ; FIA%ice_cover(i,j) = 0.0
        FIA%WindStr_x(i,j) = FIA%flux_u_top(i,j,0)
        FIA%WindStr_y(i,j) = FIA%flux_u_top(i,j,0)
      endif
    enddo
  enddo

  !
  ! set count to zero and fluxes will be zeroed before the next sum
  !
  FIA%avg_count = 0
end subroutine avg_top_quantities


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> do_update_ice_model_fast applies the surface heat fluxes, shortwave radiation
!!   diffusion of heat to the sea-ice to implicitly determine a new temperature
!!   profile, subject to the constraint that ice and snow temperatures are never
!!   above freezing.  Melting and freezing occur elsewhere.
subroutine do_update_ice_model_fast( Atmos_boundary, IST, OSS, FIA, CS, G, IG )

  type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary
  type(ice_state_type),          intent(inout) :: IST
  type(ocean_sfc_state_type),    intent(in)    :: OSS
  type(fast_ice_avg_type),       intent(inout) :: FIA
  type(fast_thermo_CS),          pointer       :: CS
  type(SIS_hor_grid_type),       intent(inout) :: G
  type(ice_grid_type),           intent(in)    :: IG

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
  real, dimension(SZI_(G), SZJ_(G)) :: tmp_diag
  real, dimension(0:IG%NkIce) :: T_col ! The temperature of a column of ice and snow in degC.
  real, dimension(IG%NkIce)   :: S_col ! The thermodynamic salinity of a column of ice, in g/kg.
  real, dimension(0:IG%NkIce) :: enth_col   ! The enthalpy of a column of snow and ice, in enth_unit (J/kg?).
  real, dimension(0:IG%NkIce) :: SW_abs_col
  real :: dt_fast, ts_new, dts, hf, hfd, latent
  real :: hf_0    ! The positive upward surface heat flux when T_sfc = 0 C, in W m-2.
  real :: dhf_dt  ! The deriviative of the upward surface heat flux with Ts, in W m-2 C-1.
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

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_fast_thermo: Module must be initialized before it is used.")

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  i_off = LBOUND(Atmos_boundary%t_flux,1) - G%isc
  j_off = LBOUND(Atmos_boundary%t_flux,2) - G%jsc
  NkIce = IG%NkIce ; I_Nk = 1.0 / NkIce ; kg_H_Nk = IG%H_to_kg_m2 * I_Nk

  CS%n_fast = CS%n_fast + 1

  if (CS%debug) &
    call IST_chksum("Start do_update_ice_model_fast", IST, G, IG)

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,Atmos_boundary,i_off, &
!$OMP                                  j_off,flux_u,flux_v,flux_t,flux_q,flux_lw, &
!$OMP                                  flux_sw_nir_dir,flux_sw_nir_dif,               &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,               &
!$OMP                                  lprec,fprec,dhdt,dedt,drdt        )            &
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
      drdt(i,j,k) = Atmos_boundary%drdt(i2,j2,k2)
    enddo ; enddo
  enddo

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             Latent_fusion=LatHtFus, Latent_vapor=LatHtVap)

  do j=jsc,jec ; do i=isc,iec
    flux_lh(i,j,0) = LatHtVap * flux_q(i,j,0)
  enddo ; enddo

  !
  ! implicit update of ice surface temperature
  !
  dt_fast = time_type_to_real(IST%Time_step_fast)

  enth_liq_0 = Enth_from_TS(0.0, 0.0, IST%ITV) ; I_enth_unit = 1.0 / enth_units

  T_freeze_ice_top = T_Freeze(S_col(1), IST%ITV)
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,NkIce,IST,dhdt,dedt,drdt,   &
!$OMP                                  flux_sw_vis_dir,flux_sw_vis_dif,flux_sw_nir_dir, &
!$OMP                                  flux_sw_nir_dif,flux_t,flux_q,flux_lw,enth_liq_0,&
!$OMP                                  dt_fast,flux_lh,I_enth_unit,G,S_col,kg_H_Nk,     &
!$OMP                                  enth_units,LatHtFus,LatHtVap,IG,OSS,FIA,CS)      &
!$OMP                          private(T_Freeze_surf,latent,enth_col,flux_sw,dhf_dt,    &
!$OMP                                  hf_0,ts_new,dts,SW_abs_col,SW_absorbed,enth_here,&
!$OMP                                  tot_heat_in,enth_imb,norm_enth_imb     )
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    T_Freeze_surf = T_Freeze(OSS%s_surf(i,j), IST%ITV)
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
                        T_Freeze_surf, FIA%bheat(i,j), ts_new, &
                        dt_fast, NkIce, FIA%tmelt(i,j,k), FIA%bmelt(i,j,k), &
                        IST%ice_thm_CSp, IST%ITV, CS%column_check)
      IST%enth_snow(i,j,k,1) = enth_col(0)
      do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_col(m) ; enddo

      dts               = ts_new - (IST%t_surf(i,j,k)-T_0degC)
      IST%t_surf(i,j,k) = IST%t_surf(i,j,k) + dts
      flux_t(i,j,k)  = flux_t(i,j,k)  + dts * dhdt(i,j,k)
      flux_q(i,j,k)  = flux_q(i,j,k)  + dts * dedt(i,j,k)
      flux_lw(i,j,k) = flux_lw(i,j,k) - dts * drdt(i,j,k)
      flux_lh(i,j,k) = latent * flux_q(i,j,k)

      if (CS%column_check) then
        SW_absorbed = SW_abs_col(0)
        do m=1,NkIce ; SW_absorbed = SW_absorbed + SW_abs_col(m) ; enddo
        IST%heat_in(i,j,k) = IST%heat_in(i,j,k) + dt_fast * &
          ((flux_lw(i,j,k) + IST%sw_abs_sfc(i,j,k)*flux_sw) + SW_absorbed + &
           FIA%bheat(i,j) - (flux_t(i,j,k) + flux_lh(i,j,k)))

        enth_here = (IG%H_to_kg_m2*IST%mH_snow(i,j,k)) * enth_col(0)
        do m=1,NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k)*kg_H_Nk) * enth_col(m)
        enddo
        tot_heat_in = enth_units * (IST%heat_in(i,j,k) - &
                                    (FIA%bmelt(i,j,k) + FIA%tmelt(i,j,k)))
        enth_imb = enth_here - (IST%enth_prev(i,j,k) + tot_heat_in)
        if (abs(enth_imb) > CS%imb_tol * (abs(enth_here) + &
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

  call sum_top_quantities(FIA, Atmos_boundary, flux_u, flux_v, flux_t, &
    flux_q, flux_sw_nir_dir, flux_sw_nir_dif, flux_sw_vis_dir, flux_sw_vis_dif, &
    flux_lw, lprec, fprec, flux_lh, G, IG )

  IST%Time = IST%Time + IST%Time_step_fast ! advance time

  if (CS%debug) &
    call IST_chksum("End do_update_ice_model_fast", IST, G, IG)

  if (CS%bounds_check) &
    call IST_bounds_check(IST, G, IG, "End of update_ice_fast", OSS=OSS)

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

  call callTree_leave("SIS_fast_thermo_init()")

end subroutine SIS_fast_thermo_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_fast_thermo_end deallocates any memory associated with this module.
subroutine SIS_fast_thermo_end(CS)
  type(fast_thermo_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine SIS_fast_thermo_end

end module SIS_fast_thermo
