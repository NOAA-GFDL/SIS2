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
!   This module handles the main updates of the ice states at the slower time- !
! scales of the couplng or the interactions with the ocean, including the ice  !
! mass balance and related thermodynamics, the ice momentum balance, ice and   !
! tracer transport, salinity changes, and coupling with the ocean.  The        !
! radiative heating and diffusive temperature changes due to coupling with the !
! atmosphere are handled elsewhere.                                            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_slow_mod

use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use MOM_checksums,     only :  chksum, Bchksum, hchksum, uchksum, vchksum
use SIS_error_checking, only : check_redundant_B, check_redundant_C
! use SIS_get_input, only : Get_SIS_input, directories
use SIS_sum_output, only : write_ice_statistics! , SIS_sum_output_init
use SIS_sum_output, only : accumulate_bottom_input, accumulate_input_1, accumulate_input_2

use mpp_domains_mod,  only  : domain2D !, mpp_get_compute_domain, CORNER, EAST, NORTH
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : fill_symmetric_edges !, MOM_domains_init, clone_MOM_domain
! use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
! use MOM_file_parser, only : open_param_file, close_param_file
use MOM_hor_index, only : hor_index_type ! , hor_index_init
! use MOM_string_functions, only : uppercase
use MOM_EOS, only : EOS_type, calculate_density_derivs

use fms_mod, only : clock_flag_default !, file_exist
! use fms_io_mod, only : restore_state, query_initialized
use fms_io_mod, only : register_restart_field, restart_file_type
use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_LOOP, CLOCK_ROUTINE

use time_manager_mod, only : time_type, time_type_to_real, get_date, get_time
use time_manager_mod, only : set_date, set_time, operator(+), operator(-)
use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)
use data_override_mod, only : data_override

use ice_type_mod, only : ice_state_type
! use ice_type_mod, only : dealloc_IST_arrays, ice_state_register_restarts
! use ice_type_mod, only : ice_diagnostics_init
use ice_type_mod, only : IST_chksum,  IST_bounds_check
use ice_utils_mod, only : get_avg, post_avg, ice_line !, ice_grid_chksum
use SIS_hor_grid, only : SIS_hor_grid_type

use ice_grid, only : ice_grid_type
use ice_spec_mod, only : get_sea_surface

! use SIS_tracer_registry, only : register_SIS_tracer, register_SIS_tracer_pair
! use SIS_tracer_flow_control, only : SIS_call_tracer_register, SIS_tracer_flow_control_init
use SIS_tracer_flow_control, only : SIS_call_tracer_column_fns
! use SIS_tracer_flow_control, only : SIS_tracer_flow_control_end


use ice_thm_mod,   only: e_to_melt, ice5lay_resize
use SIS2_ice_thm,  only: get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,  only: ice_resize_SIS2, add_frazil_SIS2, rebalance_ice_layers
use SIS2_ice_thm,  only: enth_from_TS, Temp_from_En_S
use SIS2_ice_thm,  only: T_freeze, enthalpy_liquid, calculate_T_freeze
use SIS_dyn_bgrid, only: SIS_B_dynamics, SIS_B_dyn_init, SIS_B_dyn_register_restarts, SIS_B_dyn_end
use SIS_dyn_cgrid, only: SIS_C_dynamics, SIS_C_dyn_init, SIS_C_dyn_register_restarts, SIS_C_dyn_end
use ice_transport_mod, only : ice_transport, ice_transport_init, ice_transport_end
use ice_transport_mod, only : adjust_ice_categories
use ice_bergs,        only: icebergs, icebergs_run, icebergs_init, icebergs_end, icebergs_incr_mass

implicit none ; private

#include <SIS2_memory.h>

public :: update_ice_model_slow, SIS_slow_register_restarts, SIS_slow_init, SIS_slow_end

integer :: iceClock, iceCLock2, iceClock4, iceClock5, &
           iceClock6, iceClock7, iceClock8, iceClock9, iceClocka, iceClockb, iceClockc

contains

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
! update_ice_model_slow - do ice dynamics, transport, and mass changes         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine update_ice_model_slow(IST, icebergs_CS, G, IG)

  type(ice_state_type),    intent(inout) :: IST
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(inout) :: IG
  type(icebergs),          pointer       :: icebergs_CS

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: h2o_chg_xprt, mass, tmp2d
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,IG%NkIce) :: &
    temp_ice    ! A diagnostic array with the ice temperature in degC.
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
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
    WindStr_x_A_in, &   ! Initial zonal (_x_) and meridional (_y_) wind stresses
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
  integer ::iyr, imon, iday, ihr, imin, isec

  real, dimension(IG%NkIce) :: S_col ! Specified thermodynamic salinity of each
                                    ! ice layer if spec_thermo_sal is true.
  real :: heat_fill_val   ! A value of enthalpy to use for massless categories.
  logical :: spec_thermo_sal
  logical :: do_temp_diags
  real :: enth_units, I_enth_units
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    rdg_frac, & ! fraction of ridged ice per category
    mi_old      ! Ice mass per unit area before thermodynamics.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rdg_open, & ! formation rate of open water due to ridging
    rdg_vosh, & ! rate of ice volume shifted from level to ridged ice
!   rdg_s2o, &  ! snow volume [m] dumped into ocean during ridging
    rdg_rate, & ! Niki: Where should this come from?
    snow2ocn
  real    :: tmp3
  mi_old(:,:,:) = 0.0
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; NkIce = IG%NkIce
  I_Nk = 1.0 / NkIce
  dt_slow = time_type_to_real(IST%Time_step_slow) ; Idt_slow = 1.0/dt_slow

  if (IST%specified_ice) then
    ndyn_steps = 0.0 ; dt_slow_dyn = 0.0
  else
    ndyn_steps = 1
    if ((IST%dt_ice_dyn > 0.0) .and. (IST%dt_ice_dyn < dt_slow)) &
      ndyn_steps = max(CEILING(dt_slow/IST%dt_ice_dyn - 0.000001), 1)
    dt_slow_dyn = dt_slow / ndyn_steps
  endif

  IST%n_calls = IST%n_calls + 1
  IST%stress_count = 0

  if (IST%debug) then
    call IST_chksum("Start update_ice_model_slow", IST, G, IG)
  endif

  if (IST%bounds_check) &
    call IST_bounds_check(IST, G, IG, "Start of update_ice_model_slow")

  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

  call post_flux_diagnostics(IST, G, IG, Idt_slow) ! save out diagnostics of fluxes.

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IST)
  do j=jsc,jec ; do i=isc,iec
    IST%frazil_input(i,j) = IST%frazil(i,j)

    IST%Enth_Mass_in_atm(i,j) = 0.0 ; IST%Enth_Mass_out_atm(i,j) = 0.0
    IST%Enth_Mass_in_ocn(i,j) = 0.0 ; IST%Enth_Mass_out_ocn(i,j) = 0.0
  enddo ; enddo

  !
  ! conservation checks: top fluxes
  !
  call mpp_clock_begin(iceClock7)
  call accumulate_input_1(IST, dt_slow, G, IG, IST%sum_output_CSp)
  if (IST%column_check) &
    call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IG, IST%sum_output_CSp, &
                              message="    Start of update", check_column=.true.)
  call mpp_clock_end(iceClock7)

  ! Determine the fractional ice coverage and the wind stresses averaged
  ! across all the ice thickness categories on an A-grid.  This is done
  ! over the entire data domain for safety.
  WindStr_x_A(:,:) = 0.0 ; WindStr_y_A(:,:) = 0.0 ; ice_cover(:,:) = 0.0
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,ncat,WindStr_x_A,WindStr_y_A, &
!$OMP                                  IST,ice_cover,ice_free,WindStr_x_ocn_A,       &
!$OMP                                  WindStr_y_ocn_A,ice_cover_in,ice_free_in,     &
!$OMP                                  WindStr_x_A_in,WindStr_y_A_in)                &
!$OMP                          private(I_wts)
  do j=jsd,jed
    do k=1,ncat ; do i=isd,ied
      WindStr_x_A(i,j) = WindStr_x_A(i,j) + IST%part_size(i,j,k) * IST%flux_u_top(i,j,k)
      WindStr_y_A(i,j) = WindStr_y_A(i,j) + IST%part_size(i,j,k) * IST%flux_v_top(i,j,k)
      ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
    enddo ; enddo
    do i=isd,ied
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
    enddo
   enddo


  ! Calve off icebergs and integrate forward iceberg trajectories
  if (IST%do_icebergs) then
    call mpp_clock_end(iceClock2) ; call mpp_clock_end(iceClock) ! Stop the sea-ice clocks.
    H_to_m_ice = IG%H_to_kg_m2 / IST%Rho_ice
    call get_avg(IST%mH_ice, IST%part_size(:,:,1:), hi_avg, wtd=.true.)
    hi_avg(:,:) = hi_avg(:,:) * H_to_m_Ice
    
    !### I think that there is long-standing bug here, in that the old ice-ocean
    !###  stresses are being passed in place of the wind stresses on the icebergs. -RWH
    if (IST%Cgrid_dyn) then
      call icebergs_run( icebergs_CS, IST%Time, &
              IST%calving(isc:iec,jsc:jec), IST%u_ocn_C(isc-2:iec+1,jsc-1:jec+1), &
              IST%v_ocn_C(isc-1:iec+1,jsc-2:jec+1), IST%u_ice_C(isc-2:iec+1,jsc-1:jec+1), &
              IST%v_ice_C(isc-1:iec+1,jsc-2:jec+1), &
              IST%flux_u_ocn(isc:iec,jsc:jec), IST%flux_v_ocn(isc:iec,jsc:jec), &
              IST%sea_lev(isc-1:iec+1,jsc-1:jec+1), IST%t_surf(isc:iec,jsc:jec,0),  &
              IST%calving_hflx(isc:iec,jsc:jec), ice_cover(isc-1:iec+1,jsc-1:jec+1), &
              hi_avg(isc-1:iec+1,jsc-1:jec+1), stagger=CGRID_NE, &
              stress_stagger=IST%flux_uv_stagger)
    else
      call icebergs_run( icebergs_CS, IST%Time, &
              IST%calving(isc:iec,jsc:jec), IST%u_ocn(isc-1:iec+1,jsc-1:jec+1), &
              IST%v_ocn(isc-1:iec+1,jsc-1:jec+1), IST%u_ice_B(isc-1:iec+1,jsc-1:jec+1), &
              IST%v_ice_B(isc-1:iec+1,jsc-1:jec+1), &
              IST%flux_u_ocn(isc:iec,jsc:jec), IST%flux_v_ocn(isc:iec,jsc:jec), &
              IST%sea_lev(isc-1:iec+1,jsc-1:jec+1), IST%t_surf(isc:iec,jsc:jec,0),  &
              IST%calving_hflx(isc:iec,jsc:jec), ice_cover(isc-1:iec+1,jsc-1:jec+1), &
              hi_avg(isc-1:iec+1,jsc-1:jec+1), stagger=BGRID_NE, &
              stress_stagger=IST%flux_uv_stagger)
    endif
    call mpp_clock_begin(iceClock) ; call mpp_clock_begin(iceClock2) ! Restart the sea-ice clocks.
  endif

  !
  ! Thermodynamics
  !
  if (.not.IST%specified_ice) then
    !TOM> Store old ice mass per unit area for calculating partial ice growth.  
    mi_old = IST%mH_ice
    
    !TOM> derive ridged ice fraction prior to thermodynamic changes of ice thickness
    !     in order to subtract ice melt proportionally from ridged ice volume (see below)
    if (IST%do_ridging) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,rdg_frac) &
!$OMP                          private(tmp3)
      do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
        tmp3 = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
        rdg_frac(i,j,k) = 0.0 ; if (tmp3 > 0.0) &
            rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp3
      enddo ; enddo ; enddo
    endif

    call disable_SIS_averaging(IST%diag)

    ! The thermodynamics routines return updated values of the ice and snow
    ! masses-per-unit area and enthalpies.
    call accumulate_input_2(IST, IST%part_size, dt_slow, G, IG, IST%sum_output_CSp)
    if (IST%SIS1_5L_thermo) then
      call SIS1_5L_thermodynamics(IST, G, IG)
    else
      call SIS2_thermodynamics(IST, G, IG)
    endif

    call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

    !TOM> calculate partial ice growth for ridging and aging.
    if (IST%do_ridging) then
      !     ice growth (IST%mH_ice > mi_old) does not affect ridged ice volume
      !     ice melt   (IST%mH_ice < mi_old) reduces ridged ice volume proportionally
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,mi_old,rdg_frac)
      do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
        if (IST%mH_ice(i,j,k) < mi_old(i,j,k)) &
          IST%rdg_mice(i,j,k) = IST%rdg_mice(i,j,k) + rdg_frac(i,j,k) * &
             (IST%mH_ice(i,j,k) - mi_old(i,j,k)) * IST%part_size(i,j,k)
        IST%rdg_mice(i,j,k) = max(IST%rdg_mice(i,j,k), 0.0)
      enddo ; enddo ; enddo
    endif

    !  Other routines that do thermodynamic vertical processes should be added here


    ! Do tracer column physics
    call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)
    call SIS_call_tracer_column_fns(dt_slow, G, IG, IST%SIS_tracer_flow_CSp, IST%mH_ice, mi_old)
    call disable_SIS_averaging(IST%diag)

    call accumulate_bottom_input(IST, dt_slow, G, IG, IST%sum_output_CSp)

    if (IST%column_check) &
      call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IG, IST%sum_output_CSp, &
                                message="      Post_thermo A", check_column=.true.)
    call adjust_ice_categories(IST%mH_ice, IST%mH_snow, IST%part_size, &
                               IST%TrReg, G, IG, IST%ice_transport_CSp) !Niki: add ridging?
    call pass_var(IST%part_size, G%Domain)
    call pass_var(IST%mH_ice, G%Domain, complete=.false.)
    call pass_var(IST%mH_snow, G%Domain, complete=.true.)

    if (IST%column_check) &
      call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IG, IST%sum_output_CSp, &
                                message="      Post_thermo B ", check_column=.true.)
  endif

  if (IST%id_xprt>0) then
    ! Store values to determine the ice and snow mass change due to transport.
    h2o_chg_xprt(:,:) = 0.0
  endif

  do nds=1,ndyn_steps

    call enable_SIS_averaging(dt_slow_dyn, IST%Time - set_time(int((ndyn_steps-nds)*dt_slow_dyn)), IST%diag)

    ! Correct the wind stresses for changes in the fractional ice-coverage.
    ice_cover(:,:) = 0.0
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,ncat,ice_cover,IST,ice_free, &
!$OMP                                  ice_cover_in,WindStr_x_A,WindStr_x_A_in,     &
!$OMP                                  WindStr_y_A,WindStr_y_A_in,ice_free_in,      &
!$OMP                                  WindStr_x_ocn_A,WindStr_y_ocn_A)
    do j=jsd,jed
      do k=1,ncat ; do i=isd,ied
        ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
      enddo ; enddo
      do i=isd,ied
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
      enddo
    enddo

    !
    ! Dynamics - update ice velocities.
    !
    call mpp_clock_begin(iceClock4)

    ms_sum(:,:) = 0.0 ; mi_sum(:,:) = 0.0
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,ncat,ms_sum,mi_sum,G,IST,IG)
    do j=jsd,jed ; do k=1,ncat ; do i=isd,ied
      ms_sum(i,j) = ms_sum(i,j) + (IG%H_to_kg_m2 * IST%mH_snow(i,j,k)) * IST%part_size(i,j,k)
      mi_sum(i,j) = mi_sum(i,j) + (IG%H_to_kg_m2 * IST%mH_ice(i,j,k))  * IST%part_size(i,j,k)
    enddo ; enddo ; enddo

    ! In the dynamics code, only the ice velocities are changed, and the ice-ocean
    ! stresses are calculated.  The gravity wave dynamics (i.e. the continuity
    ! equation) are not included in the dynamics.  All of the thickness categories
    ! are merged together.
    if (IST%Cgrid_dyn) then
      if (IST%area_wtd_stress) then
        !   The j-loop extents here are larger than they would normally be in case
        ! the stresses are being passed to the ocean on a B-grid.
!$OMP parallel default(none) shared(isc,iec,jsc,jec,G,ice_cover,WindStr_x_Cu,ice_free, &
!$OMP                               WindStr_x_A,WindStr_x_ocn_Cu,WindStr_x_ocn_A,      &
!$OMP                               WindStr_y_Cv,WindStr_y_A,WindStr_y_ocn_Cv,         &
!$OMP                               WindStr_y_ocn_A) &
!$OMP                       private(weights,I_wts)
!$OMP do
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
!$OMP end do nowait
!$OMP do
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
!$OMP end parallel
      else
        WindStr_x_Cu(:,:) = 0.0 ; WindStr_x_ocn_Cu(:,:) = 0.0 ; wts(:,:) = 0.0
        WindStr_y_Cv(:,:) = 0.0 ; WindStr_y_ocn_Cv(:,:) = 0.0 ; wts(:,:) = 0.0
!$OMP parallel default(none) shared(isc,iec,jsc,jec,G,ncat,IST,wts,WindStr_x_Cu,    &
!$OMP                               WindStr_x_ocn_Cu,WindStr_y_Cv,WindStr_y_ocn_Cv) &
!$OMP                       private(ps_vel)
!$OMP do
        do j=jsc-1,jec+1
          do k=1,ncat ; do I=isc-1,iec
            ps_vel = 0.5*G%mask2dCu(I,j) * (IST%part_size(i+1,j,k) + IST%part_size(i,j,k))
            WindStr_x_Cu(I,j) = WindStr_x_Cu(I,j) + ps_vel * (G%mask2dCu(I,j) * &
                             0.5* (IST%flux_u_top(i,j,k) + IST%flux_u_top(i+1,j,k)) )
            wts(I,J) = wts(I,J) + ps_vel
          enddo ; enddo
          do I=isc-1,iec
            if (wts(I,j) > 0.) WindStr_x_Cu(I,j) = WindStr_x_Cu(I,j) / wts(I,j)

            WindStr_x_ocn_Cu(I,j) = G%mask2dCu(I,j) * &
                       0.5 * (IST%flux_u_top(i,j,0) + IST%flux_u_top(i+1,j,0))
          enddo
        enddo
!$OMP end do nowait
!$OMP do
        do J=jsc-1,jec
          do k=1,ncat ; do i=isc-1,iec+1
            ps_vel = 0.5*G%mask2dCv(i,J) * (IST%part_size(i,j+1,k) + IST%part_size(i,j,k))
            WindStr_y_Cv(i,j) = WindStr_y_Cv(i,J) + ps_vel * ( G%mask2dCv(i,J) * &
                           0.5*(IST%flux_v_top(i,j,k) + IST%flux_v_top(i,j+1,k)) )
            wts(i,J) = wts(i,J) + ps_vel
          enddo ; enddo
          do i=isc-1,iec+1
            if (wts(i,J) > 0.) WindStr_y_Cv(i,J) = WindStr_y_Cv(i,J) / wts(i,J)

            WindStr_y_ocn_Cv(i,J) = G%mask2dCv(i,J) * &
                       0.5*(IST%flux_v_top(i,j,0) + IST%flux_v_top(i,j+1,0))
          enddo
        enddo
!$OMP end parallel
      endif

      if (IST%debug) then
        call IST_chksum("Before SIS_C_dynamics", IST, G, IG)
        call hchksum(IST%part_size(:,:,0), "ps(0) before SIS_C_dynamics", G%HI)
        call hchksum(ms_sum, "ms_sum before SIS_C_dynamics", G%HI)
        call hchksum(mi_sum, "mi_sum before SIS_C_dynamics", G%HI)
        call hchksum(IST%sea_lev, "sea_lev before SIS_C_dynamics", G%HI, haloshift=1)
        call uchksum(IST%u_ocn_C, "u_ocn_C before SIS_C_dynamics", G%HI)
        call vchksum(IST%v_ocn_C, "v_ocn_C before SIS_C_dynamics", G%HI)
        call uchksum(WindStr_x_Cu, "WindStr_x_Cu before SIS_C_dynamics", G%HI)
        call vchksum(WindStr_y_Cv, "WindStr_y_Cv before SIS_C_dynamics", G%HI)
        call check_redundant_C("WindStr before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G)
      endif

      call mpp_clock_begin(iceClocka)
      !### Ridging needs to be added with C-grid dynamics.
      call SIS_C_dynamics(1.0-IST%part_size(:,:,0), ms_sum, mi_sum, IST%u_ice_C, IST%v_ice_C, &
                        IST%u_ocn_C, IST%v_ocn_C, &
                        WindStr_x_Cu, WindStr_y_Cv, IST%sea_lev, str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, &
                        dt_slow_dyn, G, IST%SIS_C_dyn_CSp)
      call mpp_clock_end(iceClocka)

      if (IST%debug) then
        call IST_chksum("After ice_dynamics", IST, G, IG)
      endif

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
      call pass_vector(str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, G%Domain, stagger=CGRID_NE)
      call mpp_clock_end(iceClockb)
      !
      ! Dynamics diagnostics
      !
      call mpp_clock_begin(iceClockc)
      if (IST%id_fax>0) call post_data(IST%id_fax, WindStr_x_Cu, IST%diag)
      if (IST%id_fay>0) call post_data(IST%id_fay, WindStr_y_Cv, IST%diag)

      call set_ocean_top_stress_Cgrid(IST, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, &
                                      str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, IST%part_size, G, IG)
      call mpp_clock_end(iceClockc)

    else ! B-grid dynamics.

      if (IST%area_wtd_stress) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,ice_cover,WindStr_x_B,ice_free, &
!$OMP                                  WindStr_x_A,WindStr_x_ocn_B,WindStr_x_ocn_A,      &
!$OMP                                  WindStr_y_ocn_B,WindStr_y_ocn_A,WindStr_y_B,      &
!$OMP                                  WindStr_y_A)                                      &
!$OMP                          private(weights,I_wts)
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
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,ncat,IST,WindStr_x_B,wts,  &
!$OMP                                  WindStr_y_B,WindStr_x_ocn_B,WindStr_y_ocn_B) &
!$OMP                          private(ps_vel)
        do J=jsc-1,jec
          do k=1,ncat ; do I=isc-1,iec
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
          enddo ; enddo
          do I=isc-1,iec
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
          enddo
        enddo
      endif

      if (IST%debug) then
        call IST_chksum("Before ice_dynamics", IST, G, IG)
        call hchksum(IST%part_size(:,:,0), "ps(0) before ice_dynamics", G%HI)
        call hchksum(ms_sum, "ms_sum before ice_dynamics", G%HI)
        call hchksum(mi_sum, "mi_sum before ice_dynamics", G%HI)
        call hchksum(IST%sea_lev, "sea_lev before ice_dynamics", G%HI, haloshift=1)
        call Bchksum(IST%u_ocn, "u_ocn before ice_dynamics", G%HI, symmetric=.true.)
        call Bchksum(IST%v_ocn, "v_ocn before ice_dynamics", G%HI, symmetric=.true.)
        call Bchksum(WindStr_x_B, "WindStr_x_B before ice_dynamics", G%HI, symmetric=.true.)
        call Bchksum(WindStr_y_B, "WindStr_y_B before ice_dynamics", G%HI, symmetric=.true.)
        call check_redundant_B("WindStr before ice_dynamics",WindStr_x_B, WindStr_y_B, G)
      endif

      rdg_rate(:,:) = 0.0
      call mpp_clock_begin(iceClocka)
      call SIS_B_dynamics(1.0-IST%part_size(:,:,0), ms_sum, mi_sum, IST%u_ice_B, IST%v_ice_B, &
                        IST%u_ocn, IST%v_ocn, WindStr_x_B, WindStr_y_B, IST%sea_lev, &
                        str_x_ice_ocn_B, str_y_ice_ocn_B, IST%do_ridging, &
                        rdg_rate(isc:iec,jsc:jec), dt_slow_dyn, G, IST%SIS_B_dyn_CSp)
      call mpp_clock_end(iceClocka)

      if (IST%debug) then
        call IST_chksum("After ice_dynamics", IST, G, IG)
      endif

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
      call mpp_clock_end(iceClockb)

      call mpp_clock_begin(iceClockc)
      !
      ! Dynamics diagnostics
      !
      if ((IST%id_fax>0) .or. (IST%id_fay>0)) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,IST,diagVarBx,diagVarBy, &
!$OMP                                  WindStr_x_ocn_B,WindStr_y_ocn_B,           &
!$OMP                                  WindStr_x_B,WindStr_y_B) &
!$OMP                          private(ps_vel)
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

      call set_ocean_top_stress_Bgrid(IST, WindStr_x_ocn_B, WindStr_y_ocn_B, &
                                      str_x_ice_ocn_B, str_y_ice_ocn_B, IST%part_size, G, IG)
      call mpp_clock_end(iceClockc)
    endif ! End of B-grid dynamics


    call mpp_clock_end(iceClock4)

    call enable_SIS_averaging(dt_slow_dyn, IST%Time - set_time(int((ndyn_steps-nds)*dt_slow_dyn)), IST%diag)
    !
    ! Do ice transport ... all ocean fluxes have been calculated by now.
    !
    call mpp_clock_begin(iceClock8)

    if (IST%id_xprt>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,h2o_chg_xprt,IST,G,IG)
      do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
        h2o_chg_xprt(i,j) = h2o_chg_xprt(i,j) - IST%part_size(i,j,k) * &
                          IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
      enddo ; enddo ; enddo
    endif

    if (IST%debug) then
      call IST_chksum("Before ice_transport", IST, G, IG)
    endif

    if (IST%Cgrid_dyn) then
      call ice_transport(IST%part_size, IST%mH_ice, IST%mH_snow, IST%u_ice_C, IST%v_ice_C, &
                         IST%TrReg, IST%sea_lev, dt_slow_dyn, G, IG, IST%ice_transport_CSp, &
                         IST%rdg_mice, snow2ocn, rdg_rate, &
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
                         IST%TrReg, IST%sea_lev, dt_slow_dyn, G, IG, IST%ice_transport_CSp, &
                         IST%rdg_mice, snow2ocn, rdg_rate, &
                         rdg_open, rdg_vosh)
    endif
    if (IST%column_check) &
      call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IG, IST%sum_output_CSp, &
                                message="      Post_transport")! , check_column=.true.)

    if (IST%id_xprt>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,h2o_chg_xprt,IST,G,IG)
      do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      h2o_chg_xprt(i,j) = h2o_chg_xprt(i,j) + IST%part_size(i,j,k) * &
                        IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
    enddo ; enddo ; enddo ; endif

    call mpp_clock_end(iceClock8)

  enddo ! nds=1,ndyn_steps
  call finish_ocean_top_stresses(IST, G%HI)

  ! Add snow volume dumped into ocean to flux of frozen precipitation:
  !### WARNING - rdg_s2o is never calculated!!!
!  if (IST%do_ridging) then ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
!    IST%fprec_top(i,j,k) = IST%fprec_top(i,j,k) + rdg_s2o(i,j)*(IST%Rho_snow/dt_slow)
!  enddo ; enddo ; enddo ; endif

  call mpp_clock_begin(iceClock8)

  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

  ! Set appropriate surface quantities in categories with no ice.  Change <1e-10 to == 0?
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)<1e-10) &
    IST%t_surf(i,j,k) = T_0degC + T_Freeze(IST%s_surf(i,j),IST%ITV)
  enddo ; enddo ; enddo

  if (IST%bounds_check) call IST_bounds_check(IST, G, IG, "After ice_transport")
  if (IST%debug) call IST_chksum("After ice_transport", IST, G, IG)

  ! Sum the concentration weighted mass for diagnostics.
  if (IST%id_mi>0 .or. IST%id_mib>0) then
    mass(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,mass,G,IST,IG)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      mass(i,j) = mass(i,j) + (IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))) * &
                  IST%part_size(i,j,k)
    enddo ; enddo ; enddo
    if (IST%id_mi>0) call post_data(IST%id_mi, mass(isc:iec,jsc:jec), IST%diag)

    if (IST%id_mib>0) then
      if (IST%do_icebergs) call icebergs_incr_mass(icebergs_CS, mass(isc:iec,jsc:jec)) ! Add icebergs mass in kg/m^2
      call post_data(IST%id_mib, mass(isc:iec,jsc:jec), IST%diag)
    endif
  endif

  if (IST%specified_ice) then   ! over-write changes with specifications.
    h_ice_input(:,:) = 0.0
    call get_sea_surface(IST%Time, IST%t_surf(isc:iec,jsc:jec,0), IST%part_size(isc:iec,jsc:jec,:), &
                         h_ice_input(isc:iec,jsc:jec))
    do j=jsc,jec ; do i=isc,iec
      IST%mH_ice(i,j,1) = h_ice_input(i,j) * (IG%kg_m2_to_H * IST%Rho_ice)
    enddo ; enddo
    call pass_var(IST%part_size, G%Domain)
  endif

  call mpp_clock_end(iceClock8)

  !
  ! Thermodynamic state diagnostics
  !
  call mpp_clock_begin(iceClock9)
  if (IST%id_cn>0) call post_data(IST%id_cn, IST%part_size(:,:,1:ncat), IST%diag)
  ! TK Mod: 10/18/02
  !  if (IST%id_obs_cn>0) call post_data(IST%id_obs_cn, Obs_cn_ice(:,:,2), IST%diag)
  ! TK Mod: 10/18/02: (commented out...does not compile yet... add later)
  !  if (IST%id_obs_hi>0) &
  !    call post_avg(IST%id_obs_hi, Obs_h_ice(isc:iec,jsc:jec), IST%part_size(isc:iec,jsc:jec,1:), &
  !                  IST%diag, G=G, wtd=.true.)

  !   Convert from ice and snow enthalpy back to temperature for diagnostic purposes.
  do_temp_diags = (IST%id_tsn > 0)
  do m=1,NkIce ; if (IST%id_t(m)>0) do_temp_diags = .true. ; enddo
  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             specified_thermo_salinity=spec_thermo_sal)
  I_enth_units = 1.0 / enth_units

  if (do_temp_diags) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,IST,spec_thermo_sal,temp_ice, &
!$OMP                                  S_col,temp_snow,ncat,NkIce)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
        if (spec_thermo_sal) then ; do m=1,NkIce
          temp_ice(i,j,k,m) = temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV)
        enddo ; else ; do m=1,NkIce
          temp_ice(i,j,k,m) = temp_from_En_S(IST%enth_ice(i,j,k,m), &
                                              IST%sal_ice(i,j,k,m), IST%ITV)
        enddo ; endif
      else
        do m=1,NkIce ; temp_ice(i,j,k,m) = 0.0 ; enddo
      endif
      if (IST%part_size(i,j,k)*IST%mH_snow(i,j,k) > 0.0) then
        temp_snow(i,j,k) = temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
      else
        temp_snow(i,j,k) = 0.0 ! ### Should this be = temp_ice(i,j,k,1)?
      endif
    enddo ; enddo ; enddo
  endif

  if (IST%id_ext>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (IST%part_size(i,j,0) < 0.85) diagVar(i,j) = 1.0
    enddo ; enddo
    call post_data(IST%id_ext, diagVar, IST%diag)
  endif
  if (IST%id_hs>0) call post_avg(IST%id_hs, IST%mH_snow, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, &
                                 scale=IG%H_to_kg_m2/IST%Rho_snow, wtd=.true.)
  if (IST%id_hi>0) call post_avg(IST%id_hi, IST%mH_ice, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, &
                                 scale=IG%H_to_kg_m2/IST%Rho_ice, wtd=.true.)
  if (IST%id_ts>0) call post_avg(IST%id_ts, IST%t_surf(:,:,1:), IST%part_size(:,:,1:), &
                                 IST%diag, G=G, offset=-T_0degC, wtd=.true.)
  if (IST%id_tsn>0) call post_avg(IST%id_tsn, temp_snow, IST%part_size(:,:,1:), &
                                 IST%diag, G=G, wtd=.true.)
  do m=1,NkIce
    if (IST%id_t(m)>0) call post_avg(IST%id_t(m), temp_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   IST%diag, G=G, wtd=.true.)
    if (IST%id_sal(m)>0) call post_avg(IST%id_sal(m), IST%sal_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   IST%diag, G=G, wtd=.true.)
  enddo
  if (IST%id_t_iceav>0) call post_avg(IST%id_t_iceav, temp_ice, IST%part_size(:,:,1:), &
                                    IST%diag, G=G, wtd=.true.)
  if (IST%id_S_iceav>0) call post_avg(IST%id_S_iceav, IST%sal_ice, IST%part_size(:,:,1:), &
                                    IST%diag, G=G, wtd=.true.)


  if (IST%id_xprt>0) then
    call post_data(IST%id_xprt, h2o_chg_xprt(isc:iec,jsc:jec)*864e2*365/dt_slow, &
                   IST%diag)
  endif
  if (IST%id_e2m>0) then
    tmp2d(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,G,tmp2d,I_enth_units, &
!$OMP                                  spec_thermo_sal,NkIce,I_Nk,S_col,IG)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)>0.0) then
      tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k)*IST%mH_snow(i,j,k)*IG%H_to_kg_m2 * &
                       ((enthalpy_liquid_freeze(0.0, IST%ITV) - &
                         IST%enth_snow(i,j,k,1)) * I_enth_units)
      if (spec_thermo_sal) then ; do m=1,NkIce
        tmp2d(i,j) = tmp2d(i,j) + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*IG%H_to_kg_m2*I_Nk) * &
                       ((enthalpy_liquid_freeze(S_col(m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
      enddo ; else ; do m=1,NkIce
        tmp2d(i,j) = tmp2d(i,j) + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*IG%H_to_kg_m2*I_Nk) * &
                       ((enthalpy_liquid_freeze(IST%sal_ice(i,j,k,m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
      enddo ; endif
    endif ; enddo ; enddo ; enddo
    call post_data(IST%id_e2m,  tmp2d(:,:), IST%diag)
  endif
  call disable_SIS_averaging(IST%diag)

  !
  ! Ridging diagnostics
  !
  !TOM> preparing output field fraction of ridged ice rdg_frac = (ridged ice volume) / (total ice volume)
  !     in each category; IST%rdg_mice is ridged ice mass per unit total
  !     area throughout the code.
  if (IST%do_ridging) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,G,rdg_frac,IG) &
!$OMP                          private(tmp3)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      tmp3 = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
      if (tmp3*IG%H_to_kg_m2 > IST%Rho_Ice*1.e-5) then   ! 1 mm ice thickness x 1% ice concentration
        rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp3
      else
        rdg_frac(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo

    call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

    if (IST%id_rdgr>0) call post_data(IST%id_rdgr, rdg_rate(isc:iec,jsc:jec), IST%diag)
!    if (IST%id_rdgf>0) call post_data(IST%id_rdgf, rdg_frac(isc:iec,jsc:jec), IST%diag)
!    if (IST%id_rdgo>0) call post_data(IST%id_rdgo, rdg_open(isc:iec,jsc:jec), IST%diag)
!    if (IST%id_rdgv>0) then
!      do j=jsc,jec ; do i=isc,iec
!        tmp2d(i,j) = rdg_vosh(i,j) * G%areaT(i,j) * G%mask2dT(i,j)
!      enddo ; enddo
!      call post_data(IST%id_rdgv, tmp2d, IST%diag)
!    endif
  endif

  if (IST%verbose) then
    call get_date(IST%Time, iyr, imon, iday, ihr, imin, isec)
    call get_time(IST%Time-set_date(iyr,1,1,0,0,0),isec,iday)
    call ice_line(iyr, iday+1, isec, IST%part_size(isc:iec,jsc:jec,0), &
                              IST%t_surf(:,:,0)-T_0degC, G)
  endif

  call mpp_clock_end(iceClock9)

  if (IST%debug) then
    call IST_chksum("End UIMS", IST, G, IG)
  endif

  if (IST%bounds_check) then
    call IST_bounds_check(IST, G, IG, "End of update_ice_slow")
  endif

  if (IST%Time + (IST%Time_step_slow/2) > IST%write_ice_stats_time) then
    call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IG, IST%sum_output_CSp)
    IST%write_ice_stats_time = IST%write_ice_stats_time + IST%ice_stats_interval
  elseif (IST%column_check) then
    call write_ice_statistics(IST, IST%Time, IST%n_calls, G, IG, IST%sum_output_CSp)
  endif

end subroutine update_ice_model_slow

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! finish_ocean_top_stresses - Finish setting the ice-ocean stresses by dividing!
!   them through the stresses by the number of times they have been augmented. !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine finish_ocean_top_stresses(IST, HI)
  type(hor_index_type), intent(in)    :: HI
  type(ice_state_type), intent(inout) :: IST

  real :: I_count
  integer :: i, j, isc, iec, jsc, jec
  isc = HI%isc ; iec = HI%iec ; jsc = HI%jsc ; jec = HI%jec

  if (IST%stress_count > 1) then
    I_count = 1.0 / IST%stress_count
    do j=jsc,jec ; do i=isc,iec
      IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) * I_count
      IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) * I_count
    enddo ; enddo
  endif

end subroutine finish_ocean_top_stresses

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ocean_top_stress_Bgrid - Calculate the stresses on the ocean integrated  !
!   across all the thickness categories with the appropriate staggering, and   !
!   store them in the public ice data type for use by the ocean model.  This   !
!   version of the routine uses wind and ice-ocean stresses on a B-grid.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ocean_top_stress_Bgrid(IST, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, part_size, G, IG)
  type(ice_state_type),    intent(inout) :: IST
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(inout) :: IG
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: windstr_x_water, str_ice_oce_x
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: windstr_y_water, str_ice_oce_y
  real, dimension (SZI_(G),SZJ_(G),0:IG%CatIce), intent(in) :: part_size

  real    :: ps_vel ! part_size interpolated to a velocity point, nondim.
  integer :: i, j, k, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  if (IST%debug) then
    call IST_chksum("Start set_ocean_top_stress_Bgrid", IST, G, IG)
  endif

  if (IST%stress_count == 0) then
    IST%flux_u_ocn(:,:) = 0.0 ; IST%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Bgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.
!$OMP parallel default(none) shared(isc,iec,jsc,jec,ncat,G,                    &
!$OMP                               part_size,windstr_x_water,windstr_y_water, &
!$OMP                               str_ice_oce_x,str_ice_oce_y)               &
!$OMP                       private(ps_vel)
  if (IST%flux_uv_stagger == AGRID) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = G%mask2dT(i,j) * part_size(i,j,0)
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + ps_vel * 0.25 * &
            ((windstr_x_water(I,J) + windstr_x_water(I-1,J-1)) + &
             (windstr_x_water(I-1,J) + windstr_x_water(I,J-1)))
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + ps_vel * 0.25 * &
            ((windstr_y_water(I,J) + windstr_y_water(I-1,J-1)) + &
             (windstr_y_water(I-1,J) + windstr_y_water(I,J-1)))
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + part_size(i,j,k) * 0.25 * &
            ((str_ice_oce_x(I,J) + str_ice_oce_x(I-1,J-1)) + &
             (str_ice_oce_x(I-1,J) + str_ice_oce_x(I,J-1)))
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + part_size(i,j,k) * 0.25 * &
            ((str_ice_oce_y(I,J) + str_ice_oce_y(I-1,J-1)) + &
             (str_ice_oce_y(I-1,J) + str_ice_oce_y(I,J-1)))
      endif ; enddo ; enddo
    enddo
  elseif (IST%flux_uv_stagger == BGRID_NE) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dBu(I,J)>0.5) ps_vel = &
                           0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                                 (part_size(i+1,j,0) + part_size(i,j+1,0)) )
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + windstr_x_water(I,J) * ps_vel
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + windstr_y_water(I,J) * ps_vel
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dBu(I,J)>0.5) then
        ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                         (part_size(i+1,j,k) + part_size(i,j+1,k)) )
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + str_ice_oce_x(I,J) * ps_vel
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + str_ice_oce_y(I,J) * ps_vel
      endif ; enddo ; enddo
    enddo
  elseif (IST%flux_uv_stagger == CGRID_NE) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.5) ps_vel = &
                           0.5*(part_size(i+1,j,0) + part_size(i,j,0))
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + ps_vel * &
                0.5 * (windstr_x_water(I,J) + windstr_x_water(I,J-1))
        ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.5) ps_vel = &
                           0.5*(part_size(i,j+1,0) + part_size(i,j,0))
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + ps_vel * &
                0.5 * (windstr_y_water(I,J) + windstr_y_water(I-1,J))
      enddo
      do k=1,ncat ; do i=isc,iec
        if (G%mask2dCu(I,j)>0.5) then
          ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
          IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + ps_vel * &
              0.5 * (str_ice_oce_x(I,J) + str_ice_oce_x(I,J-1))
        endif
        if (G%mask2dCv(i,J)>0.5) then
          ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
          IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + ps_vel * &
                  0.5 * (str_ice_oce_y(I,J) + str_ice_oce_y(I-1,J))
        endif
      enddo ; enddo
    enddo
  else
!$OMP single
    call SIS_error(FATAL, "set_ocean_top_stress_Bgrid: Unrecognized flux_uv_stagger.")
!$OMP end single
  endif
!$OMP end parallel
  IST%stress_count = IST%stress_count + 1

  if (IST%debug) then
    call IST_chksum("End set_ocean_top_stress_Bgrid", IST, G, IG)
  endif

end subroutine set_ocean_top_stress_Bgrid

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ocean_top_stress_Cgrid - Calculate the stresses on the ocean integrated  !
!   across all the thickness categories with the appropriate staggering, and   !
!   store them in the public ice data type for use by the ocean model.  This   !
!   version of the routine uses wind and ice-ocean stresses on a C-grid.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ocean_top_stress_Cgrid(IST, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, part_size, G, IG)
  type(ice_state_type),    intent(inout) :: IST
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(inout) :: IG
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: windstr_x_water, str_ice_oce_x
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: windstr_y_water, str_ice_oce_y
  real, dimension (SZI_(G),SZJ_(G),0:IG%CatIce), intent(in) :: part_size

  real    :: ps_vel ! part_size interpolated to a velocity point, nondim.
  integer :: i, j, k, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  if (IST%debug) then
    call IST_chksum("Start set_ocean_top_stress_Cgrid", IST, G, IG)
  endif

  if (IST%stress_count == 0) then
    IST%flux_u_ocn(:,:) = 0.0 ; IST%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Cgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.
!$OMP parallel default(none) shared(isc,iec,jsc,jec,ncat,G,    &
!$OMP                               part_size,windstr_x_water,windstr_y_water, &
!$OMP                               str_ice_oce_x,str_ice_oce_y)               &
!$OMP                       private(ps_vel)
  if (IST%flux_uv_stagger == AGRID) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = G%mask2dT(i,j) * part_size(i,j,0)
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + ps_vel * 0.5 * &
                            (windstr_x_water(I,j) + windstr_x_water(I-1,j))
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + ps_vel * 0.5 * &
                            (windstr_y_water(I,j) + windstr_y_water(i,J-1))
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) +  part_size(i,j,k) * 0.5 * &
                            (str_ice_oce_x(I,j) + str_ice_oce_x(I-1,j))
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + part_size(i,j,k) * 0.5 * &
                            (str_ice_oce_y(I,j) + str_ice_oce_y(i,J-1))
      endif ; enddo ; enddo
    enddo
  elseif (IST%flux_uv_stagger == BGRID_NE) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dBu(I,J)>0.5) ps_vel = &
                           0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                                 (part_size(i+1,j,0) + part_size(i,j+1,0)) )
        ! Consider deleting the masks here?
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
                (windstr_x_water(I,j) + windstr_x_water(I,j+1))
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
                (windstr_y_water(I,j) + windstr_y_water(i+1,J))
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dBu(I,J)>0.5) then
        ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                         (part_size(i+1,j,k) + part_size(i,j+1,k)) )
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + ps_vel * 0.5 * &
                            (str_ice_oce_x(I,j) + str_ice_oce_x(I,j+1))
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + ps_vel * 0.5 * &
                            (str_ice_oce_y(I,j) + str_ice_oce_y(i+1,J))
      endif ; enddo ; enddo
    enddo
  elseif (IST%flux_uv_stagger == CGRID_NE) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.5) ps_vel = &
                           0.5*(part_size(i+1,j,0) + part_size(i,j,0))
        IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + ps_vel * windstr_x_water(I,j)
        ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.5) ps_vel = &
                           0.5*(part_size(i,j+1,0) + part_size(i,j,0))
        IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + ps_vel * windstr_y_water(i,J)
      enddo
      do k=1,ncat ; do i=isc,iec
        if (G%mask2dCu(I,j)>0.5) then
          ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
          IST%flux_u_ocn(i,j) = IST%flux_u_ocn(i,j) + ps_vel * str_ice_oce_x(I,j)
        endif
        if (G%mask2dCv(i,J)>0.5) then
          ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
          IST%flux_v_ocn(i,j) = IST%flux_v_ocn(i,j) + ps_vel * str_ice_oce_y(I,j)
        endif
      enddo ; enddo
    enddo
  else
!$OMP single
    call SIS_error(FATAL, "set_ocean_top_stress_Cgrid: Unrecognized flux_uv_stagger.")
!$OMP end single
  endif
!$OMP end parallel

  IST%stress_count = IST%stress_count + 1

  if (IST%debug) then
    call IST_chksum("End set_ocean_top_stress_Cgrid", IST, G, IG)
  endif

end subroutine set_ocean_top_stress_Cgrid

subroutine SIS1_5L_thermodynamics(IST, G, IG)
  type(ice_state_type),               intent(inout) :: IST
  type(SIS_hor_grid_type),            intent(inout) :: G
  type(ice_grid_type),                intent(inout) :: IG

  ! This subroutine does the thermodynamic calculations following SIS1.

  real, dimension(G%isc:G%iec,G%jsc:G%jec)  :: mi_change, h2o_change, bsnk, tmp2d
  real, dimension(SZI_(G),SZJ_(G),1:IG%CatIce) :: &
    snow_to_ice, &   ! The conversion from snow to ice in m.
    h_ice, &         ! The ice thickness in m.
    h_snow           ! The snow thickness in m.
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: dum1, Obs_h_ice ! for qflux calculation
  real, dimension(G%isc:G%iec,G%jsc:G%jec,2) :: Obs_cn_ice      ! partition 2 = ice concentration
  real, dimension(SZI_(G),SZJ_(G))   :: qflx_lim_ice, qflx_res_ice
  real, dimension(1:IG%CatIce)        :: e2m

  real, dimension(0:IG%NkIce) :: T_col ! The temperature of a column of ice and snow in degC.
  real, dimension(0:IG%NkIce,0:IG%CatIce) :: T_col_k
  real, dimension(IG%NkIce) :: S_col ! The thermodynamic salinity of a column of ice, in g/kg.
  real :: heat_fill_val   ! A value of enthalpy to use for massless categories.
  real :: dt_slow, Idt_slow, yr_dtslow
  integer :: i, j, k, l, m, n, isc, iec, jsc, jec, ncat, NkIce

  real :: heat_to_ocn, h2o_to_ocn, evap_from_ocn, sn2ic, bablt
  real            :: heat_limit_ice, heat_res_ice
  real            :: tot_heat, heating, tot_frazil
  real :: T_Freeze_surf
  real :: H_to_m_ice, H_to_m_snow  ! The specific volumes of ice and snow times the
                               ! conversion factor from thickness units, in m H-1.
  real :: LatHtFus     ! The latent heat of fusion of ice in J/kg.
  real :: LatHtVap     ! The latent heat of vaporization of water at 0C in J/kg.
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  NkIce = IG%NkIce
  dt_slow = time_type_to_real(IST%Time_step_slow) ; Idt_slow = 1.0/dt_slow
  H_to_m_ice = IG%H_to_kg_m2 / IST%Rho_Ice ; H_to_m_snow = IG%H_to_kg_m2 / IST%Rho_Snow

  if (NkIce /= 4) call SIS_error(FATAL, "SIS1_5L_thermodynamics requires that NK_ICE=4.")

  call mpp_clock_begin(iceClock5)

  snow_to_ice(:,:,:) = 0.0
  bsnk(:,:) = 0.0

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, Latent_fusion=LatHtFus, &
                             Latent_vapor=LatHtVap)
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,IST,mi_change,h2o_change, &
!$OMP                                  h_ice,h_snow,H_to_m_Ice,H_to_m_Snow,IG)
  do j=jsc,jec
    do i=isc,iec
      mi_change(i,j) = 0.0
      h2o_change(i,j) = 0.0
    enddo
    do k=1,ncat ; do i=isc,iec
      mi_change(i,j) = mi_change(i,J) - IG%H_to_kg_m2*IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
      h2o_change(i,j) = h2o_change(i,j) - IST%part_size(i,j,k) * &
                        IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
      h_ice(i,j,k) = IST%mH_ice(i,j,k) * H_to_m_Ice
      h_snow(i,j,k) = IST%mH_snow(i,j,k) * H_to_m_Snow
    enddo ; enddo
    ! Start accumulating certain fluxes at the ocean's surface.
    do i=isc,iec
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
    enddo
    do k=1,ncat ; do i=isc,iec
      IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + IST%part_size(i,j,k) * IST%lprec_top(i,j,k)
    enddo ; enddo
  enddo

  if (IST%num_tr_fluxes>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST)
    do n=1,IST%num_tr_fluxes
      do j=jsc,jec ; do i=isc,iec
        IST%tr_flux_ocn_top(i,j,n) = IST%part_size(i,j,0) * IST%tr_flux_top(i,j,0,n)
      enddo ; enddo
      do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
        IST%tr_flux_ocn_top(i,j,n) = IST%tr_flux_ocn_top(i,j,n) + &
                     IST%part_size(i,j,k) * IST%tr_flux_top(i,j,k,n)
      enddo ; enddo ; enddo
    enddo
  endif

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,NkIce,G,IST,h_snow,h_ice, &
!$OMP                                  dt_slow,snow_to_ice,Idt_slow,bsnk,S_col,Ice, &
!$OMP                                  LatHtFus,LatHtVap) &
!$OMP                          private(T_col,T_Freeze_surf,evap_from_ocn,h2o_to_ocn, &
!$OMP                                  heat_to_ocn,bablt,sn2ic)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
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
                                 ((LatHtVap*evap_from_ocn)*Idt_slow)
      IST%flux_t_ocn_top(i,j) = IST%flux_t_ocn_top(i,j) + IST%part_size(i,j,k) * &
            (IST%bheat(i,j) - (heat_to_ocn - LatHtFus*evap_from_ocn)*Idt_slow)
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
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,NkIce,qflx_lim_ice,qflx_res_ice, &
!$OMP                                  IST,S_col,h_ice,ncat,h_snow,Obs_h_ice,dt_slow,   &
!$OMP                                  Obs_cn_ice,snow_to_ice,LatHtFus)  &
!$OMP                          private(T_col_k,heat_res_ice,heat_limit_ice,T_Freeze_surf, &
!$OMP                                  e2m,tot_heat,heating,evap_from_ocn,h2o_to_ocn,   &
!$OMP                                  heat_to_ocn,bablt,sn2ic)
    do j=jsc,jec ; do i=isc,iec
      do k=0,ncat
         T_col_k(0,k) = Temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
         do m=1,NkIce ; T_col_k(m,k) = Temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV) ; enddo
      enddo
      heat_res_ice   = 0.0
      heat_limit_ice = 0.0
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j),IST%ITV)
      !
      ! calculate enthalpy
      !
      if (IST%slab_ice) then
        e2m(1) = h_ice(i,j,1)*IST%Rho_ice*LatHtFus
      else
        do k=1,ncat
          if ((IST%part_size(i,j,k)>0.0 .and. h_ice(i,j,k)>0.0)) then
             e2m(k) = e_to_melt(h_snow(i,j,k), T_col_k(0,k), h_ice(i,j,k), &
                      T_col_k(1,k), T_col_k(2,k), T_col_k(3,k), T_col_k(4,k)) * IST%part_size(i,j,k)
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
          heat_res_ice = -(LatHtFus*IST%Rho_ice*Obs_h_ice(i,j)-sum(e2m)) &
                         *dt_slow/(86400*IST%ice_restore_timescale)
        else
          heat_res_ice = -(LatHtFus*IST%Rho_ice*Obs_h_ice(i,j)*Obs_cn_ice(i,j,2)-sum(e2m)) &
                         *dt_slow/(86400*IST%ice_restore_timescale)
        endif
      endif

      if (IST%do_ice_limit .and. (sum(e2m) > IST%max_ice_limit*IST%Rho_ice*LatHtFus)) then
        heat_limit_ice = sum(e2m)-LatHtFus*IST%Rho_ice*IST%max_ice_limit
        ! should we "heat_ice_res = 0.0" ?
      endif

      !
      ! apply constraining heat to ice
      !
      tot_heat = heat_res_ice+heat_limit_ice
      if (IST%slab_ice) h_ice(i,j,1) = h_ice(i,j,1) - tot_heat/(IST%Rho_ice*LatHtFus)

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
              call ice5lay_resize(h_snow(i,j,k), T_col_k(0,k), h_ice(i,j,k), &
                                  T_col_k(1,k), T_col_k(2,k), T_col_k(3,k), T_col_k(4,k), &
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

        call ice5lay_resize(h_snow(i,j,k), T_col_k(0,k), h_ice(i,j,k), &
                            T_col_k(1,k), T_col_k(2,k), T_col_k(3,k), T_col_k(4,k), 0.0,   &
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
        e2m(1) = e2m(1) - h_ice(i,j,1)*IST%Rho_ice*LatHtFus
      else
        do k=1,ncat
          if (IST%part_size(i,j,k)>0.0 .and. h_ice(i,j,k)>0.0) &
            e2m(k) = e2m(k)-e_to_melt(h_snow(i,j,k), T_col_k(0,k), h_ice(i,j,k), &
                     T_col_k(1,k), T_col_k(2,k), T_col_k(3,k), T_col_k(4,k)) * IST%part_size(i,j,k)
        enddo
      endif
      ! if (abs(sum(e2m) - heat_res_ice - heat_limit_ice)>IST%Rho_ice*LI*1e-3) &
      !       print *, 'QFLUX conservation error at', i, j, 'heat2ice=',  &
      !             tot_heat, 'melted=', sum(e2m), 'h*part_size=', &
      !             h_ice(i,j,:)*IST%part_size(i,j,:)

      do k=1,ncat
         IST%enth_snow(i,j,k,1) = enth_from_TS(T_col_k(0,k), 0.0, IST%ITV)
         do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_from_TS(T_col_k(m,k), S_col(m), IST%ITV) ; enddo
      enddo
    enddo ; enddo
  endif ! End of (IST%do_ice_restore .or. IST%do_ice_limit) block
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    IST%mH_ice(i,j,k) = h_ice(i,j,k) * (IG%kg_m2_to_H * IST%Rho_ice)
    IST%mH_snow(i,j,k) = h_snow(i,j,k) * (IG%kg_m2_to_H * IST%Rho_snow)
  enddo ; enddo ; enddo

  !   Convert from ice temperature (which is not conserved) to enthalpy, which
  ! includes the heat requirements for melting of brine pockets associated with
  ! temperature changes.
  heat_fill_val = Enth_from_TS(0.0, 0.0, IST%ITV)

  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) <= 0.0) then
      do m=1,NkIce ; IST%enth_ice(i,j,k,m) = heat_fill_val ; enddo
      IST%enth_snow(i,j,k,1) = heat_fill_val
    endif
  enddo ; enddo ; enddo

  call mpp_clock_end(iceClock6)

  ! Determine the salt fluxes to ocean
  ! Note that at this point mi_change and h2o_change are the negative of the masses.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,IST,mi_change,h2o_change, &
!$OMP                                  Idt_slow,IG)
  do j=jsc,jec
    do k=1,ncat ; do i=isc,iec
      mi_change(i,j) = mi_change(i,J) + (IG%H_to_kg_m2*IST%mH_ice(i,j,k))*IST%part_size(i,j,k)
      h2o_change(i,j) = h2o_change(i,j) + IST%part_size(i,j,k) * &
                        (IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k)))
    enddo ; enddo
    do i=isc,iec
      ! Note the conversion here from g m-2 to kg m-2 s-1.
      IST%flux_salt(i,j) = (0.001*IST%ice_bulk_salin) * &
                             (mi_change(i,j) * Idt_slow)
    enddo
  enddo

  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

  yr_dtslow = (864e2*365/dt_slow)
  if (IST%id_lsnk>0) then
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = min(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(IST%id_lsnk, tmp2d(isc:iec,jsc:jec), IST%diag)
  endif
  if (IST%id_lsrc>0) then
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = max(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(IST%id_lsrc, tmp2d(isc:iec,jsc:jec), IST%diag)
  endif
  if (IST%id_saltf>0) call post_data(IST%id_saltf, IST%flux_salt, IST%diag)
  if (IST%id_bsnk>0)  call post_data(IST%id_bsnk, bsnk(isc:iec,jsc:jec)*yr_dtslow, &
                                     IST%diag )
  if (IST%id_tmelt>0) call post_avg(IST%id_tmelt, IST%tmelt, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, wtd=.true.)
  if (IST%id_bmelt>0) call post_avg(IST%id_bmelt, IST%bmelt, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, wtd=.true.)
  if (IST%id_sn2ic>0) call post_avg(IST%id_sn2ic, snow_to_ice, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow)
  if (IST%id_qflim>0) call post_data(IST%id_qflim, qflx_lim_ice, IST%diag)
  if (IST%id_qfres>0) call post_data(IST%id_qfres, qflx_res_ice, IST%diag)

  call disable_SIS_averaging(IST%diag)

end subroutine SIS1_5L_thermodynamics

subroutine SIS2_thermodynamics(IST, G, IG)
  type(ice_state_type),               intent(inout) :: IST
  type(SIS_hor_grid_type),            intent(inout) :: G
  type(ice_grid_type),                intent(inout) :: IG

  ! This subroutine does the thermodynamic calculations in the same order as SIS1,
  ! but with a greater emphasis on enthalpy as the dominant state variable.

  real, dimension(G%isc:G%iec,G%jsc:G%jec)  :: salt_change, h2o_change, bsnk, tmp2d
  real, dimension(SZI_(G),SZJ_(G),1:IG%CatIce) :: snow_to_ice
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: Obs_sst, Obs_h_ice ! for qflux calculation
  real, dimension(G%isc:G%iec,G%jsc:G%jec,2) :: Obs_cn_ice      ! partition 2 = ice concentration
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: icec, icec_obs
  real, dimension(SZI_(G),SZJ_(G))   :: &
    qflx_lim_ice, qflx_res_ice, &
    net_melt              ! The net mass flux from the ice and snow into the
                          ! ocean due to melting and freezing integrated
                          ! across all categories, in kg m-2 s-1.
  real, dimension(SZI_(G),SZJ_(G),1:IG%CatIce)   :: heat_in, enth_prev, enth
  real, dimension(SZI_(G),SZJ_(G))   :: heat_in_col, enth_prev_col, enth_col, enth_mass_in_col

  real, dimension(IG%NkIce) :: S_col        ! The salinity of a column of ice, in g/kg.
  real, dimension(IG%NkIce+1) :: Salin      ! The conserved bulk salinity of each
                                           ! layer in g/kg, with the salinity of
                                           ! newly formed ice in layer NkIce+1.
  real, dimension(0:IG%NkIce) :: m_lay    ! The masses of a column of ice and snow, in kg m-2.
  real, dimension(0:IG%NkIce) :: Tcol0    ! The temperature of a column of ice and snow, in degC.
  real, dimension(0:IG%NkIce) :: S_col0   ! The salinity of a column of ice and snow, in g/kg.
  real, dimension(0:IG%NkIce) :: Tfr_col0 ! The freezing temperature of a column of ice and snow, in degC.
  real, dimension(0:IG%NkIce+1) :: &
    enthalpy              ! The initial enthalpy of a column of ice and snow
                          ! and the surface ocean, in enth_units (often J/kg).
  real, dimension(IG%CatIce) :: frazil_cat  ! The frazil heating applied to each thickness
                       ! category, averaged over the area of that category in J m-2.
  real :: enthalpy_ocean  ! The enthalpy of the ocean surface waters, in Enth_units.
  real :: heat_fill_val   ! An enthalpy to use for massless categories, in enth_units.

  real :: T_freeze_surf ! The freezing temperature at the surface salinity of
                        ! the ocean, in deg C.
  real :: I_part        ! The inverse of a part_size, nondim.
  logical :: spec_thermo_sal  ! If true, use the specified salinities of the
                              ! various sub-layers of the ice for all thermodynamic
                              ! calculations; otherwise use the prognostic
                              ! salinity fields for these calculations.

  type(EOS_type), pointer :: EOS
  real :: Cp_water
  real :: drho_dT(1), drho_dS(1), pres_0(1)

  real :: dt_slow     ! The thermodynamic step, in s.
  real :: Idt_slow    ! The inverse of the thermodynamic step, in s-1.
  real :: yr_dtslow   ! The ratio of 1 year to the thermodyamic time step, used
                      ! to change the units of several diagnostics to rate yr-1
  real :: heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, sn2ic, bablt
  real :: salt_to_ice ! The flux of salt from the ocean to the ice, in kg m-2 s-1.
                      ! This may be of either sign; in some places it is an
                      ! average over the whole cell, while in others just a partition.
  real :: mtot_ice    ! The total mass of ice ans snow in a cell, in kg m-2.
  real :: e2m_tot     ! The total enthalpy requred to melt all ice and snow, in J m-2.
  real :: enth_evap, enth_ice_to_ocn, enth_ocn_to_ice, enth_snowfall
  real :: tot_heat, heating, tot_frazil, heat_mass_in, heat_input
  real :: mass_in, mass_here, mass_prev, mass_imb
  real :: enth_units, I_enth_units ! The units of enthaply and their inverse.
  real :: frac_keep, frac_melt  ! The fraction of ice and snow to keep or remove, nd.
  real :: ice_melt_lay ! The amount of excess ice removed from each layer in kg/m2.
  real :: snow_melt    ! The amount of excess snow that is melted, in kg/m2.
  real :: enth_freeze  ! The freezing point enthalpy of a layer, in enth_units.
  real :: enth_to_melt ! The enthalpy addition required to melt the excess ice
                       ! and snow in enth_unit kg/m2.
  real :: I_Nk         ! The inverse of the number of layers in the ice, nondim.
  real :: kg_H_Nk      ! The conversion factor from units of H to kg/m2 over Nk.
  real :: part_sum     ! A running sum of partition sizes.
  real :: d_enth       ! The change in enthalpy between categories.
  real :: fill_frac    ! The fraction of the difference between the thicknesses
                       ! in thin categories that will be removed within a single
                       ! timestep with filling_frazil.
  integer :: i, j, k, l, m, n, isc, iec, jsc, jec, ncat, NkIce
  integer :: k_merge
  real :: LatHtFus     ! The latent heat of fusion of ice in J/kg.
  real :: LatHtVap     ! The latent heat of vaporization of water at 0C in J/kg.
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  real :: tot_heat_in, enth_here, enth_imb, norm_enth_imb, emic2, tot_heat_in2, enth_imb2

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  NkIce = IG%NkIce ; I_Nk = 1.0 / NkIce ; kg_H_Nk = IG%H_to_kg_m2 * I_Nk
  dt_slow = time_type_to_real(IST%Time_step_slow) ; Idt_slow = 1.0/dt_slow

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             specified_thermo_salinity=spec_thermo_sal, &
                             Latent_fusion=LatHtFus, Latent_vapor=LatHtVap)
  S_col0(0) = 0.0 ; do m=1,NkIce ; S_col0(m) = S_col(m) ; enddo
  call calculate_T_Freeze(S_col0(0:NkIce), Tfr_col0(0:NkIce), IST%ITV)
  I_enth_units = 1.0 / enth_units

  heat_fill_val = Enth_from_TS(0.0, 0.0, IST%ITV)

  if (.not.spec_thermo_sal) call SIS_error(FATAL, "SIS2_thermodynamics is not "//&
    "prepared for SPECIFIED_THERMO_SALINITY to be false.")

  call mpp_clock_begin(iceClock6)
  ! Add any heat fluxes to restore the sea-ice properties toward a prescribed
  ! state, potentially including freshwater fluxes to avoid driving oceanic
  ! convection.
  if (IST%nudge_sea_ice) then
    IST%cool_nudge(:,:) = 0.0 ; IST%melt_nudge(:,:) = 0.0
    icec(:,:) = 0.0
    call data_override('ICE','icec',icec_obs,IST%Time)

    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      icec(i,j) = icec(i,j) + IST%part_size(i,j,k)
    enddo ; enddo ; enddo
    pres_0(:) = 0.0
    call get_SIS2_thermo_coefs(IST%ITV, Cp_SeaWater=Cp_water, EOS=EOS)
    do j=jsc,jec ; do i=isc,iec
      if (icec(i,j) < icec_obs(i,j) - IST%nudge_conc_tol) then
        IST%cool_nudge(i,j) = IST%nudge_sea_ice_rate * &
             ((icec_obs(i,j)-IST%nudge_conc_tol) - icec(i,j))**2.0 ! W/m2
        if (IST%nudge_stab_fac /= 0.0) then
          if (IST%t_ocn(i,j) > T_Freeze(IST%s_surf(i,j), IST%ITV) ) then
            call calculate_density_derivs(IST%t_ocn(i:i,j),IST%s_surf(i:i,j),pres_0,&
                           drho_dT,drho_dS,1,1,EOS)
            IST%melt_nudge(i,j) = IST%nudge_stab_fac * (-IST%cool_nudge(i,j)*drho_dT(1)) / &
                                  ((Cp_Water*drho_dS(1)) * max(IST%s_surf(i,j), 1.0) )
          endif
        endif
      elseif (icec(i,j) > icec_obs(i,j) + IST%nudge_conc_tol) then
        ! Heat the ice but do not apply a fresh water flux.
        IST%cool_nudge(i,j) = -IST%nudge_sea_ice_rate * &
             (icec(i,j) - (icec_obs(i,j)+IST%nudge_conc_tol))**2.0 ! W/m2
      endif

      if (IST%cool_nudge(i,J) > 0.0) then
        IST%frazil(i,j) = IST%frazil(i,j) + IST%cool_nudge(i,j)*dt_slow
      elseif (IST%cool_nudge(i,J) < 0.0) then
        IST%bheat(i,j) = IST%bheat(i,j) - IST%cool_nudge(i,j)
      endif
    enddo ; enddo
  endif
  if (IST%do_ice_restore) then
    ! get observed ice thickness for ice restoring, if calculating qflux
    call get_sea_surface(IST%Time, Obs_SST, Obs_cn_ice, Obs_h_ice)
    call get_SIS2_thermo_coefs(IST%ITV, Latent_fusion=LatHtFus)
    qflx_res_ice(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      e2m_tot = 0.0
      ! Calculate the average ice and snow enthalpy relative to freezing per unit area.
      do k=1,ncat
        if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)>0.0) then
          e2m_tot = (IST%part_size(i,j,k)*IST%mH_snow(i,j,k)) * IG%H_to_kg_m2 * &
                       ((enthalpy_liquid_freeze(0.0, IST%ITV) - &
                         IST%enth_snow(i,j,k,1)) * I_enth_units)
          if (spec_thermo_sal) then ; do m=1,NkIce
            e2m_tot = e2m_tot + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*kg_H_Nk) * &
                       ((enthalpy_liquid_freeze(S_col(m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
          enddo ; else ; do m=1,NkIce
            e2m_tot = e2m_tot + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*kg_H_Nk) * &
                       ((enthalpy_liquid_freeze(IST%sal_ice(i,j,k,m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
          enddo ; endif
        endif
      enddo
      qflx_res_ice(i,j) = -(LatHtFus*IST%Rho_ice*Obs_h_ice(i,j)*Obs_cn_ice(i,j,2)-e2m_tot) / &
                           (86400.0*IST%ice_restore_timescale)
      if (qflx_res_ice(i,j) < 0.0) then
        IST%frazil(i,j) = IST%frazil(i,j) - qflx_res_ice(i,j)*dt_slow
      elseif (qflx_res_ice(i,j) >  0.0) then
        IST%bheat(i,j) = IST%bheat(i,j) + qflx_res_ice(i,j)
      endif
    enddo ; enddo
  endif
  call mpp_clock_end(iceClock6)


  if (IST%column_check) then
    enth_prev(:,:,:) = 0.0 ; heat_in(:,:,:) = 0.0
    enth_prev_col(:,:) = 0.0 ; heat_in_col(:,:) = 0.0 ; enth_mass_in_col(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,enth_prev,G,I_Nk, &
!$OMP                                  heat_in_col,dt_slow,enth_prev_col,NkIce)
    do j=jsc,jec
      do k=1,ncat ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
        enth_prev(i,j,k) = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,NkIce
          enth_prev(i,j,k) = enth_prev(i,j,k) + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
      endif ; enddo ; enddo

      do i=isc,iec
        heat_in_col(i,j) = heat_in_col(i,j) - IST%frazil(i,j)
        heat_in_col(i,j) = heat_in_col(i,j) - IST%part_size(i,j,0) * dt_slow*IST%flux_t_top(i,j,0)
      enddo

      do k=1,ncat ; do i=isc,iec
        if (IST%part_size(i,j,k) > 0.0) then
          heat_in_col(i,j) = heat_in_col(i,j) + IST%part_size(i,j,k) * &
              (IST%tmelt(i,j,k) + IST%bmelt(i,j,k) - dt_slow*IST%bheat(i,j))
        endif
        if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
          enth_prev_col(i,j) = enth_prev_col(i,j) + &
            (IST%mH_snow(i,j,k)*IST%part_size(i,j,k)) * IST%enth_snow(i,j,k,1)
          do m=1,NkIce
            enth_prev_col(i,j) = enth_prev_col(i,j) + &
              (IST%mH_ice(i,j,k)*IST%part_size(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
          enddo
        endif
      enddo ; enddo
    enddo
  endif


  call mpp_clock_begin(iceClock5)

  snow_to_ice(:,:,:) = 0.0 ; net_melt(:,:) = 0.0
  bsnk(:,:) = 0.0
  salt_change(:,:) = 0.0
  h2o_change(:,:) = 0.0
!$OMP parallel default(none) shared(isc,iec,jsc,jec,ncat,G,IST,salt_change,kg_H_Nk, &
!$OMP                               h2o_change,NkIce,IG)
  if (IST%ice_rel_salin <= 0.0) then
!$OMP do
    do j=jsc,jec ; do m=1,NkIce ; do k=1,ncat ; do i=isc,iec
      salt_change(i,j) = salt_change(i,j) - &
         (IST%sal_ice(i,j,k,m)*(IST%mH_ice(i,j,k)*kg_H_Nk)) * IST%part_size(i,j,k)
    enddo ; enddo ; enddo ; enddo
  endif
!$OMP do
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    h2o_change(i,j) = h2o_change(i,j) - IST%part_size(i,j,k) * &
                      IG%H_to_kg_m2*(IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
  enddo ; enddo ; enddo

  ! Start accumulating the fluxes at the ocean's surface.
!$OMP do
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
!$OMP do
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + IST%part_size(i,j,k) * IST%lprec_top(i,j,k)
  enddo ; enddo ; enddo

  if (IST%num_tr_fluxes>0) then
!$OMP do
    do n=1,IST%num_tr_fluxes
      do j=jsc,jec ; do i=isc,iec
        IST%tr_flux_ocn_top(i,j,n) = IST%part_size(i,j,0) * IST%tr_flux_top(i,j,0,n)
      enddo ; enddo
      do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
        IST%tr_flux_ocn_top(i,j,n) = IST%tr_flux_ocn_top(i,j,n) + &
                   IST%part_size(i,j,k) * IST%tr_flux_top(i,j,k,n)
      enddo ; enddo ; enddo
    enddo
  endif
!$OMP end parallel

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,IST,S_col0,NkIce,S_col, &
!$OMP                                  dt_slow,snow_to_ice,heat_in,I_NK,enth_units,   &
!$OMP                                  enth_prev,enth_mass_in_col,Idt_slow,bsnk,      &
!$OMP                                  salt_change,net_melt,kg_H_nk,LatHtFus,LatHtVap,IG) &
!$OMP                          private(mass_prev,enthalpy,enthalpy_ocean,Salin,     &
!$OMP                                  heat_to_ocn,h2o_ice_to_ocn,h2o_ocn_to_ice,   &
!$OMP                                  evap_from_ocn,salt_to_ice,bablt,enth_evap,   &
!$OMP                                  enth_ice_to_ocn,enth_ocn_to_ice,heat_input,  &
!$OMP                                  heat_mass_in,mass_in,mass_here,enth_here,    &
!$OMP                                  tot_heat_in,enth_imb,mass_imb,norm_enth_imb, &
!$OMP                                  m_lay, mtot_ice,                             &
!$OMP                                  T_Freeze_surf,I_part,sn2ic,enth_snowfall)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
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

      m_lay(0) = IST%mH_snow(i,j,k) * IG%H_to_kg_m2
      do m=1,NkIce ; m_lay(m) = IST%mH_ice(i,j,k) * kg_H_Nk ; enddo

      call ice_resize_SIS2(m_lay, enthalpy, S_col, Salin, &
                   IST%fprec_top(i,j,k)*dt_slow, IST%flux_q_top(i,j,k)*dt_slow, &
                   IST%tmelt(i,j,k), IST%bmelt(i,j,k), NkIce, &
                   heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                   snow_to_ice(i,j,k), salt_to_ice, IST%ITV, IST%ice_thm_CSp, bablt, &
                   enth_evap, enth_ice_to_ocn, enth_ocn_to_ice)

      IST%mH_snow(i,j,k) = m_lay(0) * IG%kg_m2_to_H
      call rebalance_ice_layers(m_lay, mtot_ice, Enthalpy, Salin, NkIce)
      IST%mH_ice(i,j,k) = mtot_ice * IG%kg_m2_to_H

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
                     (heat_to_ocn - (LatHtVap+LatHtFus)*evap_from_ocn)

        heat_input = IST%tmelt(i,j,k) + IST%bmelt(i,j,k) - (heat_to_ocn - (LatHtVap+LatHtFus)*evap_from_ocn)
        heat_mass_in = enth_snowfall + enth_ocn_to_ice - enth_ice_to_ocn - enth_evap
        mass_in = dt_slow*IST%fprec_top(i,j,k) + h2o_ocn_to_ice - h2o_ice_to_ocn - &
                 (dt_slow*IST%flux_q_top(i,j,k)-evap_from_ocn)

        mass_here = IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k)
        enth_here = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
        tot_heat_in = IG%kg_m2_to_H*(enth_units*heat_input + heat_mass_in)
        mass_in = mass_in*IG%kg_m2_to_H

        enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        mass_imb = mass_here - (mass_prev + mass_in)
        if (abs(enth_imb) > IST%imb_tol * &
            (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
          norm_enth_imb = enth_imb / (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in))
          enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        endif

        enth_mass_in_col(i,j) = enth_mass_in_col(i,j) + IST%part_size(i,j,k) * &
          (enth_snowfall + enth_ocn_to_ice - enth_ice_to_ocn - enth_evap)
      endif

      IST%flux_q_ocn_top(i,j) = IST%flux_q_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                  (evap_from_ocn*Idt_slow) ! no ice, evaporation left
      IST%flux_lh_ocn_top(i,j) = IST%flux_lh_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                 ((LatHtVap*evap_from_ocn)*Idt_slow)
      IST%flux_t_ocn_top(i,j) = IST%flux_t_ocn_top(i,j) + IST%part_size(i,j,k) * &
             (IST%bheat(i,j) - (heat_to_ocn - LatHtFus*evap_from_ocn)*Idt_slow)
      IST%flux_sw_vis_dif_ocn(i,j) = IST%flux_sw_vis_dif_ocn(i,j) + IST%part_size(i,j,k) * &
             (((IST%flux_sw_vis_dir_top(i,j,k) + IST%flux_sw_vis_dif_top(i,j,k)) + &
               (IST%flux_sw_nir_dir_top(i,j,k) + IST%flux_sw_nir_dif_top(i,j,k))) * &
               IST%sw_abs_ocn(i,j,k))
      net_melt(i,j) = net_melt(i,j) + IST%part_size(i,j,k) * &
              ((h2o_ice_to_ocn-h2o_ocn_to_ice)*Idt_slow)
      bsnk(i,j) = bsnk(i,j) - IST%part_size(i,j,k)*bablt ! bot. melt. ablation

    endif ! Applying surface fluxes to each category.
  enddo ; enddo ; enddo

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,IG,IST,S_col0,NkIce,S_col, &
!$OMP                                  dt_slow,snow_to_ice,heat_in,I_NK,enth_units,   &
!$OMP                                  enth_prev,enth_mass_in_col,Idt_slow,bsnk,      &
!$OMP                                  salt_change,net_melt,kg_h_Nk,LatHtFus)         &
!$OMP                          private(mass_prev,enthalpy,enthalpy_ocean,Salin,     &
!$OMP                                  heat_to_ocn,h2o_ice_to_ocn,h2o_ocn_to_ice,   &
!$OMP                                  evap_from_ocn,salt_to_ice,bablt,enth_evap,   &
!$OMP                                  enth_ice_to_ocn,enth_ocn_to_ice,heat_input,  &
!$OMP                                  heat_mass_in,mass_in,mass_here,enth_here,    &
!$OMP                                  tot_heat_in,enth_imb,mass_imb,norm_enth_imb, &
!$OMP                                  m_lay,mtot_ice,frazil_cat,k_merge,part_sum,  &
!$OMP                                  fill_frac,d_enth,                            &
!$OMP                                  T_Freeze_surf,I_part,sn2ic,enth_snowfall)
  do j=jsc,jec ; do i=isc,iec ; if (IST%frazil(i,j)>0.0) then

    frazil_cat(1:ncat) = 0.0
    k_merge = 1  ! Find the category that will be combined with the ice free category.
    if (.not.IST%filling_frazil) then
      do k=1,ncat ; if (IST%part_size(i,j,0)+IST%part_size(i,j,k)>0.01) then
        ! absorb frazil in thinest ice partition available    (was ...>0.0)
        ! raised above threshold from 0 to 0.01 to avert ocean-ice model blow-ups
        k_merge = k ; exit
      endif ; enddo
    endif

!   if (IST%part_size(i,j,0) > 0.0) then
      k = k_merge
      T_Freeze_surf = T_Freeze(IST%s_surf(i,j),IST%ITV)

      ! Combine the ice-free part size with one of the categories.
      I_part = 1.0 / (IST%part_size(i,j,k) + IST%part_size(i,j,0))
      IST%mH_snow(i,j,k) = (IST%mH_snow(i,j,k) * IST%part_size(i,j,k)) * I_part
      IST%mH_ice(i,j,k)  = (IST%mH_ice(i,j,k)  * IST%part_size(i,j,k)) * I_part
      IST%t_surf(i,j,k) = (IST%t_surf(i,j,k) * IST%part_size(i,j,k) + &
                       (T_0degC + T_Freeze_surf)*IST%part_size(i,j,0)) * I_part
      IST%part_size(i,j,k) = IST%part_size(i,j,k) + IST%part_size(i,j,0)
      IST%part_size(i,j,0) = 0.0
!   endif

    if (IST%filling_frazil) then
      if (IST%fraz_fill_time < 0.0) then
        frazil_cat(k) = IST%frazil(i,J)
        IST%frazil(i,j) = 0.0
      else
        part_sum = 0.0
        fill_frac = 1.0 ; if (IST%fraz_fill_time > 0.0) &
          fill_frac = dt_slow / (dt_slow + IST%fraz_fill_time)
        do k=1,ncat-1
          part_sum = part_sum + IST%part_size(i,j,k)
          d_enth = fill_frac * max(0.0, LatHtFus * IG%H_to_kg_m2 * &
                         (IG%mH_cat_bound(k+1) - IST%mH_ice(i,j,k)))
          if (d_enth*part_sum > IST%frazil(i,j)) then
            frazil_cat(k) = IST%frazil(i,j) / part_sum
            IST%frazil(i,j) = 0.0
            exit
          else
            frazil_cat(k) = d_enth
            IST%frazil(i,j) = IST%frazil(i,j) - frazil_cat(k)*part_sum
          endif
        enddo
        if (IST%frazil(i,j) > 0.0) &
          frazil_cat(ncat) = IST%frazil(i,j)
          ! Note that at this point we should have that part_sum = 1.0.
        do k=ncat-1,1 ; frazil_cat(k) = frazil_cat(k) + frazil_cat(k+1) ; enddo
      endif
    else
      ! Set the frazil that is absorbed in this category and remove it from
      ! the overall frazil energy.
      I_part = 1.0 / (IST%part_size(i,j,k))
      frazil_cat(k_merge) = IST%frazil(i,j) * I_part
      IST%frazil(i,j) = 0.0
    endif

    do k=1,ncat ; if (frazil_cat(k) > 0.0) then
      if (IST%column_check) then
        enth_prev(i,j,k) = 0.0 ; heat_in(i,j,k) = 0.0
        if (IST%mH_ice(i,j,k) > 0.0) then
          enth_prev(i,j,k) = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
          do m=1,NkIce
            enth_prev(i,j,k) = enth_prev(i,j,k) + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
          enddo
          enth_prev(i,j,k) = enth_prev(i,j,k) * IST%part_size(i,j,k)
        endif
      endif

      ! Set up the columns of enthalpy, salinity, and mass.
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

      m_lay(0) = IST%mH_snow(i,j,k) * IG%H_to_kg_m2
      do m=1,NkIce ; m_lay(m) = IST%mH_ice(i,j,k) * kg_H_Nk ; enddo

      call add_frazil_SIS2(m_lay, enthalpy, S_col, Salin, frazil_cat(k), &
                   T_Freeze_surf, NkIce, h2o_ocn_to_ice, salt_to_ice, IST%ITV, &
                   IST%ice_thm_CSp, Enthalpy_freeze=enth_ocn_to_ice)

      call rebalance_ice_layers(m_lay, mtot_ice, Enthalpy, Salin, NkIce)

      ! Unpack the columns of mass, enthalpy and salinity.
      IST%mH_snow(i,j,k) = m_lay(0) * IG%kg_m2_to_H
      IST%mH_ice(i,j,k) = mtot_ice * IG%kg_m2_to_H

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

!      IST%Enth_Mass_in_ocn(i,j) = IST%Enth_Mass_in_ocn(i,j) + &
!          IST%part_size(i,j,k) * enth_ocn_to_ice
      IST%Enth_Mass_in_ocn(i,j) = IST%Enth_Mass_in_ocn(i,j) + &
          IST%part_size(i,j,k) * (h2o_ocn_to_ice * enthalpy_ocean)
      net_melt(i,j) = net_melt(i,j) - &
             (h2o_ocn_to_ice * IST%part_size(i,j,k)) * Idt_slow

      if (IST%column_check) then
        heat_in(i,j,k) = heat_in(i,j,k) - frazil_cat(k)

        enth_here = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,NkIce
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
    endif ; enddo ! frazil_cat>0, k-loop

  endif ; enddo ; enddo   ! frazil>0, i-, and j-loops

  call mpp_clock_end(iceClock5)

  call mpp_clock_begin(iceClock6)
  if (IST%do_ice_limit) then
    qflx_lim_ice(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,NkIce,IST,G,I_enth_units,   &
!$OMP                                  spec_thermo_sal,kg_H_Nk,S_col,Obs_h_ice,dt_slow, &
!$OMP                                  Obs_cn_ice,snow_to_ice,salt_change,qflx_lim_ice, &
!$OMP                                  Idt_slow,net_melt,IG)                            &
!$OMP                          private(mtot_ice,frac_keep,frac_melt,salt_to_ice,  &
!$OMP                                  h2o_ice_to_ocn,enth_to_melt,enth_ice_to_ocn,   &
!$OMP                                  ice_melt_lay,snow_melt,enth_freeze)
    do j=jsc,jec ; do i=isc,iec
      mtot_ice = 0.0
      do k=1,ncat
         mtot_ice = mtot_ice + IST%part_size(i,j,k) * IG%H_to_kg_m2 * &
                     (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
      enddo
      if (mtot_ice > IST%max_ice_limit*IST%Rho_ice) then
        frac_keep = IST%max_ice_limit*IST%Rho_ice / mtot_ice
        frac_melt = 1.0 - frac_keep
        salt_to_ice = 0.0 ; h2o_ice_to_ocn = 0.0
        enth_to_melt = 0.0 ; enth_ice_to_ocn = 0.0
        do k=1,ncat
          ice_melt_lay = frac_melt * IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*kg_H_Nk
          snow_melt = frac_melt * IST%part_size(i,j,k)*IST%mH_snow(i,j,k)*IG%H_to_kg_m2
          enth_freeze = enthalpy_liquid_freeze(0.0, IST%ITV)
          enth_to_melt = enth_to_melt + snow_melt * &
                         (enth_freeze - IST%enth_snow(i,j,k,1))
          enth_ice_to_ocn = enth_ice_to_ocn + snow_melt * enth_freeze
          do m=1,NkIce
            if (spec_thermo_sal) then
              enth_freeze = enthalpy_liquid_freeze(S_col(m), IST%ITV)
            else
              enth_freeze = enthalpy_liquid_freeze(IST%sal_ice(i,j,k,m), IST%ITV)
            endif
            enth_to_melt = enth_to_melt + ice_melt_lay * &
                           (enth_freeze - IST%enth_ice(i,j,k,m))
            enth_ice_to_ocn = enth_ice_to_ocn + ice_melt_lay * enth_freeze
            salt_to_ice = salt_to_ice - ice_melt_lay * IST%sal_ice(i,j,k,m)
          enddo

          h2o_ice_to_ocn = h2o_ice_to_ocn + (snow_melt + NkIce*ice_melt_lay)
          IST%mH_ice(i,j,k) = frac_keep*IST%mH_ice(i,j,k)
          IST%mH_snow(i,j,k) = frac_keep*IST%mH_snow(i,j,k)
        enddo
        net_melt(i,j) = net_melt(i,j) + h2o_ice_to_ocn * Idt_slow
        qflx_lim_ice(i,j) = enth_to_melt * I_enth_units * Idt_slow
        IST%Enth_Mass_out_ocn(i,j) = IST%Enth_Mass_out_ocn(i,j) - enth_ice_to_ocn
        if (IST%ice_rel_salin > 0.0) then
          salt_change(i,j) = salt_change(i,j) + IST%part_size(i,j,k) * salt_to_ice
        endif
      endif

    enddo ; enddo
  endif ! End of (IST%do_ice_limit) block
  call mpp_clock_end(iceClock6)

  if (IST%column_check) then
    enth_col(:,:) = 0.0
    ! Add back any frazil that has not been used yet.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,heat_in_col,IST,dt_slow)
    do j=jsc,jec ; do i=isc,iec
      heat_in_col(i,j) = heat_in_col(i,j) + IST%frazil(i,j) + IST%flux_t_ocn_top(i,j)*dt_slow
    enddo ; enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,enth_col,IST,I_Nk,NkIce)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
      enth_col(i,j) = enth_col(i,j) + &
        (IST%mH_snow(i,j,k)*IST%part_size(i,j,k)) * IST%enth_snow(i,j,k,1)
      do m=1,NkIce
        enth_col(i,j) = enth_col(i,j) + &
          (IST%mH_ice(i,j,k)*IST%part_size(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
      enddo
    endif ; enddo ; enddo ; enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,enth_col,IST,enth_units,    &
!$OMP                                  heat_in_col,enth_mass_in_col,enth_prev_col) &
!$OMP                          private(enth_here,tot_heat_in,emic2,tot_heat_in2,   &
!$OMP                                  enth_imb,norm_enth_imb,enth_imb2)
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
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,ncat,salt_change,IST,kg_H_Nk,NkIce)
    do j=jsc,jec ; do m=1,NkIce ; do k=1,ncat ; do i=isc,iec
      salt_change(i,j) = salt_change(i,j) + &
         (IST%sal_ice(i,j,k,m)*(IST%mH_ice(i,j,k)*kg_H_Nk)) * IST%part_size(i,j,k)
    enddo ; enddo ; enddo ; enddo
  endif
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,G,h2o_change, &
!$OMP                                  salt_change,Idt_slow,IG) &
  do j=jsc,jec
    do k=1,ncat ; do i=isc,iec
      h2o_change(i,j) = h2o_change(i,j) + IST%part_size(i,j,k) * &
                        IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
    enddo ; enddo
    do i=isc,iec
      ! Note the conversion here from g m-2 to kg m-2 s-1.
      IST%flux_salt(i,j) = salt_change(i,j) * (0.001*Idt_slow)
    enddo
  enddo

  !   The remainder of this routine deals with any thermodynamics diagnostic
  ! output that has been requested.
  call enable_SIS_averaging(dt_slow, IST%Time, IST%diag)

  yr_dtslow = (864e2*365*Idt_slow)
  if (IST%id_lsnk>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,tmp2d,h2o_change,yr_dtslow)
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = min(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(IST%id_lsnk, tmp2d(isc:iec,jsc:jec), IST%diag)
  endif
  if (IST%id_lsrc>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,tmp2d,h2o_change,yr_dtslow)
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = max(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(IST%id_lsrc, tmp2d(isc:iec,jsc:jec), IST%diag)
  endif
  if (IST%id_saltf>0) call post_data(IST%id_saltf, IST%flux_salt, IST%diag)
  if (IST%id_bsnk>0)  call post_data(IST%id_bsnk, bsnk(isc:iec,jsc:jec)*yr_dtslow, &
                                     IST%diag)
  if (IST%id_tmelt>0) call post_avg(IST%id_tmelt, IST%tmelt, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, wtd=.true.)
  if (IST%id_bmelt>0) call post_avg(IST%id_bmelt, IST%bmelt, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow, wtd=.true.)
  if (IST%id_sn2ic>0) call post_avg(IST%id_sn2ic, snow_to_ice, IST%part_size(:,:,1:), IST%diag, G=G, &
                                    scale=Idt_slow)
  if (IST%id_qflim>0) call post_data(IST%id_qflim, qflx_lim_ice, IST%diag)
  if (IST%id_qfres>0) call post_data(IST%id_qfres, qflx_res_ice, IST%diag)

  call disable_SIS_averaging(IST%diag)

  ! Combine the liquid precipitation with the net melt of ice and snow for
  ! passing to the ocean. These may later be kept separate.
  do j=jsc,jec ; do i=isc,iec
    IST%lprec_ocn_top(i,j) = IST%lprec_ocn_top(i,j) + net_melt(i,j)
  enddo ; enddo

end subroutine SIS2_thermodynamics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_slow_register_restarts - allocate and register any variables for this    !
!      module that need to be included in the restart files.                   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_slow_register_restarts(mpp_domain, HI, IG, param_file, IST, &
                                      Ice_restart, restart_file)
  type(domain2d),          intent(in)    :: mpp_domain
  type(hor_index_type),    intent(in)    :: HI
  type(ice_grid_type),     intent(in)    :: IG     ! The sea-ice grid type
  type(param_file_type),   intent(in)    :: param_file
  type(ice_state_type),    pointer       :: IST
  type(restart_file_type), intent(inout) :: Ice_restart
  character(len=*),        intent(in)    :: restart_file

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
!

!   This subroutine registers the restart variables associated with the
! the slow ice dynamics and thermodynamics.

!  integer :: isd, ied, jsd, jed, id
!  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

!  if (associated(IST)) then
!    call SIS_error(WARNING, "SIS_slow_register_restarts called with an "//&
!                            "associated control structure.")
!    return
!  endif
!  allocate(IST)
  if (IST%Cgrid_dyn) then
    call SIS_C_dyn_register_restarts(mpp_domain, HI, param_file, &
                 IST%SIS_C_dyn_CSp, Ice_restart, restart_file)
  else
    call SIS_B_dyn_register_restarts(mpp_domain, HI, param_file, &
                 IST%SIS_B_dyn_CSp, Ice_restart, restart_file)
  endif
!  call ice_transport_register_restarts(G, param_file, IST%ice_transport_CSp, &
!                                       Ice_restart, restart_file)

end subroutine SIS_slow_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_slow_init - initializes ice model data, parameters and diagnostics       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_slow_init(Time, G, IG, param_file, diag, IST)
  type(time_type),     target, intent(in)    :: Time   ! current time
  type(SIS_hor_grid_type),     intent(in)    :: G      ! The horizontal grid structure
  type(ice_grid_type),         intent(in)    :: IG     ! The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(ice_state_type), intent(inout) :: IST

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_slow" ! This module's name.

  call callTree_enter("SIS_slow_init(), SIS_slow.F90")

  ! Read all relevant parameters and write them to the model log.
!   call log_version(param_file, mod, version, "")

  if (IST%Cgrid_dyn) then
    call SIS_C_dyn_init(IST%Time, G, param_file, IST%diag, IST%SIS_C_dyn_CSp, IST%ntrunc)
  else
    call SIS_B_dyn_init(IST%Time, G, param_file, IST%diag, IST%SIS_B_dyn_CSp)
  endif
  call ice_transport_init(IST%Time, G, param_file, IST%diag, IST%ice_transport_CSp)


  iceClock = mpp_clock_id( 'Ice', flags=clock_flag_default, grain=CLOCK_COMPONENT )
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

  call callTree_leave("SIS_slow_init()")

end subroutine SIS_slow_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_slow_end - deallocates memory                                            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_slow_end (IST)
  type(ice_state_type), pointer :: IST

  if (IST%Cgrid_dyn) then
    call SIS_C_dyn_end(IST%SIS_C_dyn_CSp)
  else
    call SIS_B_dyn_end(IST%SIS_B_dyn_CSp)
  endif
  call ice_transport_end(IST%ice_transport_CSp)

end subroutine SIS_slow_end

end module SIS_slow_mod
