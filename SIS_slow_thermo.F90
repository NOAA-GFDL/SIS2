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
! mass balance and related thermodynamics and salinity changes, and            !
! thermodynamic coupling with the ocean.  The radiative heating and diffusive  !
! temperature changes due to coupling with the atmosphere are handled elsewhere.                                            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_slow_thermo

use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_sum_output, only : SIS_sum_out_CS, write_ice_statistics! , SIS_sum_output_init
use SIS_sum_output, only : accumulate_bottom_input, accumulate_input_1, accumulate_input_2

! use MOM_domains,       only : pass_var
! ! use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_hor_index, only : hor_index_type
use MOM_EOS, only : EOS_type, calculate_density_derivs

use fms_mod, only : clock_flag_default
use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_LOOP, CLOCK_ROUTINE

use time_manager_mod, only : time_type, time_type_to_real! , get_date, get_time
! use time_manager_mod, only : set_date, set_time, operator(+), operator(-)
! use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)
use data_override_mod, only : data_override

use SIS_types, only : ice_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_types, only : ocean_sfc_state_type
use SIS_types, only : IST_chksum, IST_bounds_check
use SIS_utils, only : post_avg
use SIS_hor_grid, only : SIS_hor_grid_type

use ice_grid, only : ice_grid_type
use ice_spec_mod, only : get_sea_surface

use SIS_tracer_flow_control, only : SIS_call_tracer_column_fns

use SIS2_ice_thm, only : SIS2_ice_thm_CS, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm, only : ice_resize_SIS2, add_frazil_SIS2, rebalance_ice_layers
use SIS2_ice_thm, only : get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm, only : enth_from_TS, Temp_from_En_S
use SIS2_ice_thm, only : enthalpy_liquid, calculate_T_freeze
use SIS_transport, only : adjust_ice_categories, SIS_transport_CS
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS
use SIS_tracer_registry, only : SIS_unpack_passive_ice_tr, SIS_repack_passive_ice_tr
use SIS_tracer_registry, only : SIS_count_passive_tracers

implicit none ; private

#include <SIS2_memory.h>

public :: slow_thermodynamics, SIS_slow_thermo_init, SIS_slow_thermo_end
public :: slow_thermo_CS, SIS_slow_thermo_set_ptrs

type slow_thermo_CS ; private
  logical :: specified_ice  ! If true, the sea ice is specified and there is
                            ! no need for ice dynamics.
  real :: ice_bulk_salin ! The globally constant sea ice bulk salinity, in g/kg
                         ! that is used to calculate the ocean salt flux.
  real :: ice_rel_salin  ! The initial bulk salinity of sea-ice relative to the
                         ! salinity of the water from which it formed, nondim.

  logical :: filling_frazil ! If true, apply frazil to fill as many categories
                            ! as possible to fill in a uniform (minimum) amount
                            ! of frazil in all the thinnest categories.
                            ! Otherwise the frazil is always assigned to a
                            ! single category with part size > 0.01.
  real    :: fraz_fill_time ! A timescale with which the filling frazil causes
                            ! the thinest cells to attain similar thicknesses,
                            ! or a negative number to apply the frazil flux
                            ! uniformly, in s.
  real :: ocean_part_min    ! The minimum value for the fractional open-ocean
                            ! area.  This can be 0, but for some purposes it
                            ! may be useful to set this to a miniscule value
                            ! (like 1e-40) that will be lost to roundoff
                            ! during any sums so that the open ocean fluxes
                            ! can be used in interpolation across categories.

  logical :: do_ridging     !   If true, apply a ridging scheme to the convergent
                            ! ice.  The original SIS2 implementation is based on
                            ! work by Torge Martin.  Otherwise, ice is compressed
                            ! proportionately if the concentration exceeds 1.

  logical :: do_ice_restore ! If true, restore the sea-ice toward climatology
                            ! by applying a restorative heat flux.
  real    :: ice_restore_timescale ! The time scale for restoring ice when
                            ! do_ice_restore is true, in days.

  logical :: do_ice_limit   ! Limit the sea ice thickness to max_ice_limit.
  real    :: max_ice_limit  ! The maximum sea ice thickness, in m, when
                            ! do_ice_limit is true.

  logical :: nudge_sea_ice = .false. ! If true, nudge sea ice concentrations towards observations.
  real    :: nudge_sea_ice_rate = 0.0 ! The rate of cooling of ice-free water that
                              ! should be ice  covered in order to constrained the
                              ! ice concentration to track observations.  A suggested
                              ! value is of order 10000 W m-2.
  real    :: nudge_stab_fac   ! A factor that determines whether the buoyancy
                              ! flux associated with the sea ice nudging of
                              ! warm water includes a freshwater flux so as to
                              ! be destabilizing on net (<1), stabilizing (>1),
                              ! or neutral (=1).  The default is 1.
  real    :: nudge_conc_tol   ! The tolerance for mismatch in the sea ice concentations
                              ! before nudging begins to be applied.

  logical :: debug        ! If true, write verbose checksums for debugging purposes.
  logical :: column_check ! If true, enable the heat check column by column.
  real    :: imb_tol      ! The tolerance for imbalances to be flagged by
                          ! column_check, nondim.
  logical :: bounds_check ! If true, check for sensible values of thicknesses
                          ! temperatures, fluxes, etc.

  integer :: n_calls = 0  ! The number of times update_ice_model_slow_down
                          ! has been called.

  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                                   ! timing of diagnostic output.

  ! These are pointers to the control structures for subsidiary modules.
  type(SIS2_ice_thm_CS), pointer  :: ice_thm_CSp => NULL()

  type(SIS_transport_CS), pointer :: SIS_transport_CSp => NULL()
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()
  type(SIS_tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL()

  integer :: id_qflim=-1, id_qfres=-1, id_fwnudge=-1
  integer :: id_lsrc=-1, id_lsnk=-1, id_bsnk=-1, id_sn2ic=-1
end type slow_thermo_CS

integer :: iceClock5, iceClock6, iceClock7

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! post_flux_diagnostics - write out any diagnostics of surface fluxes.         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine post_flux_diagnostics(IST, FIA, IOF, CS, G, IG, Idt_slow)
  type(ice_state_type),      intent(in) :: IST
  type(fast_ice_avg_type),   intent(in) :: FIA
  type(ice_ocean_flux_type), intent(in) :: IOF
  type(slow_thermo_CS),      pointer    :: CS
  type(SIS_hor_grid_type),   intent(in) :: G
  type(ice_grid_type),       intent(in) :: IG
  real,                      intent(in) :: Idt_slow

  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: tmp2d, net_sw
  integer :: i, j, k, m, n, isc, iec, jsc, jec, ncat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  ! Flux diagnostics
  !
  if (FIA%id_runoff>0) &
    call post_data(FIA%id_runoff, FIA%runoff, CS%diag)
  if (FIA%id_calving>0) &
    call post_data(FIA%id_calving, FIA%calving_preberg, CS%diag)
  if (FIA%id_runoff_hflx>0) &
    call post_data(FIA%id_runoff_hflx, FIA%runoff_hflx, CS%diag)
  if (FIA%id_calving_hflx>0) &
    call post_data(FIA%id_calving_hflx, FIA%calving_hflx_preberg, CS%diag)
  ! The frazil diagnostic is with the other ocean surface diagnostics.
  ! if (IST%id_frazil>0) &
  !   call post_data(IST%id_frazil, FIA%frazil_left*Idt_slow, CS%diag)
  if (FIA%id_sh>0) call post_avg(FIA%id_sh, FIA%flux_t_top, IST%part_size, &
                                 CS%diag, G=G)
  if (FIA%id_lh>0) call post_avg(FIA%id_lh, FIA%flux_lh_top, IST%part_size, &
                                 CS%diag, G=G)
  if (FIA%id_evap>0) call post_avg(FIA%id_evap, FIA%flux_q_top, IST%part_size, &
                                 CS%diag, G=G)
  if (FIA%id_slp>0) &
    call post_data(FIA%id_slp, FIA%p_atm_surf, CS%diag)

  if ((FIA%id_sw>0) .or. (FIA%id_albedo>0)) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,net_sw,IST,FIA)
    do j=jsc,jec
      do i=isc,iec ; net_sw(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        net_sw(i,j) = net_sw(i,j) + IST%part_size(i,j,k) * ( &
              FIA%flux_sw_vis_dir_top(i,j,k) + FIA%flux_sw_vis_dif_top(i,j,k) + &
              FIA%flux_sw_nir_dir_top(i,j,k) + FIA%flux_sw_nir_dif_top(i,j,k) )
      enddo ; enddo
    enddo
    if (FIA%id_sw>0) call post_data(FIA%id_sw, net_sw, CS%diag)
    if (FIA%id_albedo>0) then
      do j=jsc,jec ; do i=isc,iec
        if (G%mask2dT(i,j)<=0.5) then
          tmp2d(i,j) = -1.0 ! This is land.
        elseif ((FIA%flux_sw_dn(i,j) > 0.0)) then
          ! The 10.0 below is deliberate.  An albedo of down to -9 can be reported
          ! for detecting inconsistent net_sw and sw_dn.
          tmp2d(i,j) = (FIA%flux_sw_dn(i,j) - min(net_sw(i,j), 10.0*FIA%flux_sw_dn(i,j))) / &
                       FIA%flux_sw_dn(i,j)
        else
          tmp2d(i,j) = 0.0 ! What does the albedo mean at night?
        endif
      enddo ; enddo
      call post_data(FIA%id_albedo, tmp2d, CS%diag)
    endif
  endif
  if (FIA%id_lw>0) call post_avg(FIA%id_lw, FIA%flux_lw_top, &
                                 IST%part_size, CS%diag, G=G)
  if (FIA%id_snofl>0) call post_avg(FIA%id_snofl, FIA%fprec_top, &
                                    IST%part_size, CS%diag, G=G)
  if (FIA%id_rain>0) call post_avg(FIA%id_rain, FIA%lprec_top, &
                                   IST%part_size, CS%diag, G=G)
  if (FIA%id_sw_dn>0) call post_data(FIA%id_sw_dn, FIA%flux_sw_dn, CS%diag)
  if (FIA%id_tsfc>0) call post_data(FIA%id_tsfc, FIA%Tskin_avg, CS%diag)
  if (FIA%id_sitemptop>0) call post_data(FIA%id_sitemptop, FIA%Tskin_avg, CS%diag)

  if (FIA%id_sh0>0) call post_data(FIA%id_sh0, FIA%flux_t0, CS%diag)
  if (FIA%id_evap0>0) call post_data(FIA%id_evap0, FIA%flux_q0, CS%diag)
  if (FIA%id_lw0>0) call post_data(FIA%id_lw0, FIA%flux_lw0, CS%diag)
  if (FIA%id_dshdt>0) call post_data(FIA%id_dshdt, FIA%dhdt, CS%diag)
  if (FIA%id_dedt>0) call post_data(FIA%id_dedt, FIA%dedt, CS%diag)
  if (FIA%id_dlwdt>0) call post_data(FIA%id_dlwdt, FIA%dlwdt, CS%diag)

  if (FIA%id_sh_cat>0) call post_data(FIA%id_sh_cat, FIA%flux_t_top, CS%diag)
  if (FIA%id_evap_cat>0) call post_data(FIA%id_evap_cat, FIA%flux_q_top, CS%diag)
  if (FIA%id_lw_cat>0) call post_data(FIA%id_lw_cat, FIA%flux_lw_top, CS%diag)
  if (FIA%id_Tsfc_cat>0) call post_data(FIA%id_Tsfc_cat, FIA%Tskin_cat, CS%diag)

  if (FIA%id_sw_vis>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,tmp2d,IST,FIA)
    do j=jsc,jec
      do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * ( &
              FIA%flux_sw_vis_dir_top(i,j,k) + FIA%flux_sw_vis_dif_top(i,j,k) )
      enddo ; enddo
    enddo
    call post_data(FIA%id_sw_vis, tmp2d, CS%diag)
  endif
  if (FIA%id_sw_dir>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,tmp2d,IST,FIA)
    do j=jsc,jec
      do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * ( &
              FIA%flux_sw_vis_dir_top(i,j,k) + FIA%flux_sw_nir_dir_top(i,j,k) )
      enddo ; enddo
    enddo
    call post_data(FIA%id_sw_dir, tmp2d, CS%diag)
  endif
  if (FIA%id_sw_dif>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,tmp2d,IST,FIA)
    do j=jsc,jec
      do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * ( &
              FIA%flux_sw_vis_dif_top(i,j,k) + FIA%flux_sw_nir_dif_top(i,j,k) )
      enddo ; enddo
    enddo
    call post_data(FIA%id_sw_dif, tmp2d, CS%diag)
  endif
  if (FIA%id_sw_nir_dir>0) call post_avg(FIA%id_sw_nir_dir, FIA%flux_sw_nir_dir_top, &
                             IST%part_size, CS%diag, G=G)
  if (FIA%id_sw_nir_dif>0) call post_avg(FIA%id_sw_nir_dif, FIA%flux_sw_nir_dif_top, &
                             IST%part_size, CS%diag, G=G)
  if (FIA%id_sw_vis_dir>0) call post_avg(FIA%id_sw_vis_dir, FIA%flux_sw_vis_dir_top, &
                             IST%part_size, CS%diag, G=G)
  if (FIA%id_sw_vis_dif>0) call post_avg(FIA%id_sw_vis_dif, FIA%flux_sw_vis_dif_top, &
                             IST%part_size, CS%diag, G=G)

  if (CS%nudge_sea_ice .and. CS%id_fwnudge>0) then
    call post_data(CS%id_fwnudge, IOF%melt_nudge, CS%diag)
  endif

end subroutine post_flux_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> slow_thermodynamics takes care of slow ice thermodynamics and mass changes
subroutine slow_thermodynamics(IST, dt_slow, CS, OSS, FIA, IOF, G, IG)

  type(ice_state_type),       intent(inout) :: IST
  real,                       intent(in)    :: dt_slow ! The thermodynamic step, in s.
  type(slow_thermo_CS),       pointer       :: CS
  type(ocean_sfc_state_type), intent(inout) :: OSS
  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(ice_ocean_flux_type),  intent(inout) :: IOF
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(inout) :: IG

  real, dimension(SZI_(G),SZJ_(G))   :: &
    h_ice_input         ! The specified ice thickness, with specified_ice, in m.

  real :: rho_ice  ! The nominal density of sea ice in kg m-3.
  real :: Idt_slow
  integer :: i, j, k, l, m, isc, iec, jsc, jec, ncat, NkIce
  integer :: isd, ied, jsd, jed

  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    rdg_frac, & ! fraction of ridged ice per category
    mi_old      ! Ice mass per unit area before thermodynamics.
  real    :: tmp3  ! This is a bad name - make it more descriptive!

  mi_old(:,:,:) = 0.0
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; NkIce = IG%NkIce
!  I_Nk = 1.0 / NkIce
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  CS%n_calls = CS%n_calls + 1

  if (CS%debug) then
    call IST_chksum("Start update_ice_model_slow", IST, G, IG)
  endif

  if (CS%bounds_check) &
    call IST_bounds_check(IST, G, IG, "Start of SIS_slow_thermo", OSS=OSS)

  ! Set the frazil heat flux that remains to be applied. This might need
  ! to be moved earlier in the algorithm, if there ever to be multiple calls to
  ! slow_thermodynamics per coupling timestep.
  do j=jsc,jec ; do i=isc,iec
    FIA%frazil_left(i,j) = OSS%frazil(i,j)
  enddo ; enddo

  !
  ! conservation checks: top fluxes
  !
  call mpp_clock_begin(iceClock7)
  call accumulate_input_1(IST, FIA, dt_slow, G, IG, CS%sum_output_CSp)
  if (CS%column_check) &
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
                              message="    SIS_slow_thermo", check_column=.true.)
  call mpp_clock_end(iceClock7)

  !
  ! Thermodynamics
  !
  if (CS%specified_ice) then   ! over-write changes with specifications.
    h_ice_input(:,:) = 0.0
    call get_sea_surface(CS%Time, OSS%SST_C(isc:iec,jsc:jec), IST%part_size(isc:iec,jsc:jec,:), &
                         h_ice_input(isc:iec,jsc:jec), ts_in_K=.false.)
    call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
    do j=jsc,jec ; do i=isc,iec
      IST%mH_ice(i,j,1) = h_ice_input(i,j) * (IG%kg_m2_to_H * rho_ice)
    enddo ; enddo

    do j=jsc,jec ; do i=isc,iec
      IOF%flux_t_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%flux_t_top(i,j,0)
      IOF%flux_q_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%flux_q_top(i,j,0)
      IOF%flux_lw_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%flux_lw_top(i,j,0)
      IOF%flux_lh_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%flux_lh_top(i,j,0)
      IOF%flux_sw_vis_dir_ocn(i,j) = IST%part_size(i,j,0) * FIA%flux_sw_vis_dir_top(i,j,0)
      IOF%flux_sw_vis_dif_ocn(i,j) = IST%part_size(i,j,0) * FIA%flux_sw_vis_dif_top(i,j,0)
      IOF%flux_sw_nir_dir_ocn(i,j) = IST%part_size(i,j,0) * FIA%flux_sw_nir_dir_top(i,j,0)
      IOF%flux_sw_nir_dif_ocn(i,j) = IST%part_size(i,j,0) * FIA%flux_sw_nir_dif_top(i,j,0)
      IOF%lprec_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%lprec_top(i,j,0)
      IOF%fprec_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%fprec_top(i,j,0)
    enddo ; enddo

  endif

  ! IOF must be updated regardless of whether the ice is specified or the prognostic model
  ! is being used
  if (FIA%num_tr_fluxes>0) then
!It is necessary and sufficient that only one OMP thread goes through the following block
!since IOF is shared between the threads (hence the block is not thread-safe).
!$OMP SINGLE
    if (IOF%num_tr_fluxes < 0) then
      ! This is the first call, and the IOF arrays need to be allocated.
      IOF%num_tr_fluxes = FIA%num_tr_fluxes

      allocate(IOF%tr_flux_ocn_top(SZI_(G), SZJ_(G), IOF%num_tr_fluxes))
      IOF%tr_flux_ocn_top(:,:,:) = 0.0
    endif
!$OMP END SINGLE
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,FIA,IOF)
    do m=1,FIA%num_tr_fluxes
      do j=jsc,jec ; do i=isc,iec
        IOF%tr_flux_ocn_top(i,j,m) = IST%part_size(i,j,0) * FIA%tr_flux_top(i,j,0,m)
      enddo ; enddo
      do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
        IOF%tr_flux_ocn_top(i,j,m) = IOF%tr_flux_ocn_top(i,j,m) + &
                   IST%part_size(i,j,k) * FIA%tr_flux_top(i,j,k,m)
      enddo ; enddo ; enddo
    enddo
  endif

  ! No other thermodynamics need to be done for ice that is specified,
  if(CS%specified_ice) return ;
  ! Otherwise, Continue with the remainder of the prognostic slow thermodynamics

  !TOM> Store old ice mass per unit area for calculating partial ice growth.
  mi_old = IST%mH_ice

  !TOM> derive ridged ice fraction prior to thermodynamic changes of ice thickness
  !     in order to subtract ice melt proportionally from ridged ice volume (see below)
  if (CS%do_ridging) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,rdg_frac) &
!$OMP                          private(tmp3)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      tmp3 = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
      rdg_frac(i,j,k) = 0.0 ; if (tmp3 > 0.0) &
          rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp3
    enddo ; enddo ; enddo
  endif

  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)

  ! Save out diagnostics of fluxes.  This must go before SIS2_thermodynamics.
  call post_flux_diagnostics(IST, FIA, IOF, CS, G, IG, Idt_slow)

  call disable_SIS_averaging(CS%diag)

  call accumulate_input_2(IST, FIA, IOF, IST%part_size, dt_slow, G, IG, CS%sum_output_CSp)
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IOF)
  do j=jsc,jec ; do i=isc,iec
    IOF%Enth_Mass_in_atm(i,j) = 0.0 ; IOF%Enth_Mass_out_atm(i,j) = 0.0
    IOF%Enth_Mass_in_ocn(i,j) = 0.0 ; IOF%Enth_Mass_out_ocn(i,j) = 0.0
  enddo ; enddo

  ! The thermodynamics routines return updated values of the ice and snow
  ! masses-per-unit area and enthalpies.
  call SIS2_thermodynamics(IST, dt_slow, CS, OSS, FIA, IOF, G, IG)

  !TOM> calculate partial ice growth for ridging and aging.
  if (CS%do_ridging) then
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
  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)
  call SIS_call_tracer_column_fns(dt_slow, G, IG, CS%tracer_flow_CSp, IST%mH_ice, mi_old)
  call disable_SIS_averaging(CS%diag)

  call accumulate_bottom_input(IST, OSS, FIA, IOF, dt_slow, G, IG, CS%sum_output_CSp)

  if (CS%column_check) &
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
                              message="      Post_thermo A", check_column=.true.)
  call adjust_ice_categories(IST%mH_ice, IST%mH_snow, IST%mH_pond, IST%part_size, &
                             IST%TrReg, G, IG, CS%SIS_transport_CSp) !Niki: add ridging?

  if (CS%column_check) &
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
                              message="      Post_thermo B ", check_column=.true.)


end subroutine slow_thermodynamics

subroutine SIS2_thermodynamics(IST, dt_slow, CS, OSS, FIA, IOF, G, IG)
  type(ice_state_type),       intent(inout) :: IST
  real,                       intent(in)    :: dt_slow ! The thermodynamic step, in s.
  type(slow_thermo_CS),       pointer       :: CS
  type(ocean_sfc_state_type), intent(in)    :: OSS
  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(ice_ocean_flux_type),  intent(inout) :: IOF
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(inout) :: IG

  ! This subroutine does the thermodynamic calculations in the same order as SIS1,
  ! but with a greater emphasis on enthalpy as the dominant state variable.

  real, dimension(G%isc:G%iec,G%jsc:G%jec)  :: salt_change, h2o_change, bsnk, tmp2d
  real, dimension(SZI_(G),SZJ_(G),1:IG%CatIce) :: snow_to_ice
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: Obs_sst, Obs_h_ice ! for qflux calculation
  real, dimension(G%isc:G%iec,G%jsc:G%jec,2) :: Obs_cn_ice      ! partition 2 = ice concentration
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: icec, icec_obs
  real, dimension(SZI_(G),SZJ_(G))   :: &
    qflx_lim_ice, qflx_res_ice, &
    cool_nudge, &         ! A heat flux out of the sea ice that
                          ! acts to create sea-ice, in W m-2.
    net_melt              ! The net mass flux from the ice and snow into the
                          ! ocean due to melting and freezing integrated
                          ! across all categories, in kg m-2 s-1.
  real, dimension(SZI_(G),SZJ_(G),1:IG%CatIce)   :: heat_in, enth_prev, enth
  real, dimension(SZI_(G),SZJ_(G))   :: heat_in_col, enth_prev_col, enth_col, enth_mass_in_col

  real, dimension(IG%NkIce) :: S_col      ! The salinity of a column of ice, in g/kg.
  real, dimension(IG%NkIce+1) :: Salin    ! The conserved bulk salinity of each
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

  real :: I_part        ! The inverse of a part_size, nondim.
  logical :: spec_thermo_sal  ! If true, use the specified salinities of the
                              ! various sub-layers of the ice for all thermodynamic
                              ! calculations; otherwise use the prognostic
                              ! salinity fields for these calculations.
  real, dimension(:,:), allocatable :: TrLay
                          ! Passive tracer slice through categories
                          ! By default, both the 0 (snow layer) boundary and
                          ! the NkIce+1 (surface ocean layer) are both set to 0
                          ! for all tracers

  type(EOS_type), pointer :: EOS
  real :: Cp_water
  real :: drho_dT(1), drho_dS(1), pres_0(1)
  real :: rho_ice     ! The nominal density of sea ice in kg m-3.

  real :: Idt_slow    ! The inverse of the thermodynamic step, in s-1.
  real :: yr_dtslow   ! The ratio of 1 year to the thermodyamic time step, used
                      ! to change the units of several diagnostics to rate yr-1
  real :: heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, sn2ic, bablt
  real :: salt_to_ice ! The flux of salt from the ocean to the ice, in kg m-2 s-1.
                      ! This may be of either sign; in some places it is an
                      ! average over the whole cell, while in others just a partition.
  real :: mtot_ice    ! The total mass of ice and snow in a cell, in kg m-2.
  real :: e2m_tot     ! The total enthalpy required to melt all ice and snow, in J m-2.
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
  real :: part_ocn     ! A slightly modified ocean part size.
  real :: d_enth       ! The change in enthalpy between categories.
  real :: fill_frac    ! The fraction of the difference between the thicknesses
                       ! in thin categories that will be removed within a single
                       ! timestep with filling_frazil.
  integer :: i, j, k, l, m, n, isc, iec, jsc, jec, ncat, NkIce, tr, npassive
  integer :: k_merge
  real :: LatHtFus     ! The latent heat of fusion of ice in J/kg.
  real :: LatHtVap     ! The latent heat of vaporization of water at 0C in J/kg.
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  real :: tot_heat_in, enth_here, enth_imb, norm_enth_imb, emic2, tot_heat_in2, enth_imb2

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  NkIce = IG%NkIce ; I_Nk = 1.0 / NkIce ; kg_H_Nk = IG%H_to_kg_m2 * I_Nk
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                   rho_ice=rho_ice, specified_thermo_salinity=spec_thermo_sal, &
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
  if (CS%nudge_sea_ice) then
    if (.not.allocated(IOF%melt_nudge)) allocate(IOF%melt_nudge(isc:iec,jsc:jec))

    cool_nudge(:,:) = 0.0 ; IOF%melt_nudge(:,:) = 0.0
    icec(:,:) = 0.0
    call data_override('ICE','icec',icec_obs,CS%Time)

    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      icec(i,j) = icec(i,j) + IST%part_size(i,j,k)
    enddo ; enddo ; enddo
    pres_0(:) = 0.0
    call get_SIS2_thermo_coefs(IST%ITV, Cp_Water=Cp_water, EOS=EOS)
    do j=jsc,jec ; do i=isc,iec
      if (icec(i,j) < icec_obs(i,j) - CS%nudge_conc_tol) then
        cool_nudge(i,j) = CS%nudge_sea_ice_rate * &
             ((icec_obs(i,j)-CS%nudge_conc_tol) - icec(i,j))**2.0 ! W/m2
        if (CS%nudge_stab_fac /= 0.0) then
          if (OSS%SST_C(i,j) > OSS%T_fr_ocn(i,j)) then
            call calculate_density_derivs(OSS%SST_C(i:i,j),OSS%s_surf(i:i,j),pres_0,&
                           drho_dT,drho_dS,1,1,EOS)
            IOF%melt_nudge(i,j) = CS%nudge_stab_fac * (-cool_nudge(i,j)*drho_dT(1)) / &
                                  ((Cp_Water*drho_dS(1)) * max(OSS%s_surf(i,j), 1.0) )
          endif
        endif
      elseif (icec(i,j) > icec_obs(i,j) + CS%nudge_conc_tol) then
        ! Heat the ice but do not apply a fresh water flux.
        cool_nudge(i,j) = -CS%nudge_sea_ice_rate * &
             (icec(i,j) - (icec_obs(i,j)+CS%nudge_conc_tol))**2.0 ! W/m2
      endif

      if (cool_nudge(i,J) > 0.0) then
        FIA%frazil_left(i,j) = FIA%frazil_left(i,j) + cool_nudge(i,j)*dt_slow
      elseif (cool_nudge(i,J) < 0.0) then
        FIA%bheat(i,j) = FIA%bheat(i,j) - cool_nudge(i,j)
      endif
    enddo ; enddo
  endif
  if (CS%do_ice_restore) then
    ! get observed ice thickness for ice restoring, if calculating qflux
    call get_sea_surface(CS%Time, Obs_SST, Obs_cn_ice, Obs_h_ice)
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
      qflx_res_ice(i,j) = -(LatHtFus*Rho_ice*Obs_h_ice(i,j)*Obs_cn_ice(i,j,2)-e2m_tot) / &
                           (86400.0*CS%ice_restore_timescale)
      if (qflx_res_ice(i,j) < 0.0) then
        FIA%frazil_left(i,j) = FIA%frazil_left(i,j) - qflx_res_ice(i,j)*dt_slow
      elseif (qflx_res_ice(i,j) >  0.0) then
        FIA%bheat(i,j) = FIA%bheat(i,j) + qflx_res_ice(i,j)
      endif
    enddo ; enddo
  endif
  call mpp_clock_end(iceClock6)


  if (CS%column_check) then
    enth_prev(:,:,:) = 0.0 ; heat_in(:,:,:) = 0.0
    enth_prev_col(:,:) = 0.0 ; heat_in_col(:,:) = 0.0 ; enth_mass_in_col(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,enth_prev,G,I_Nk, &
!$OMP                                  heat_in_col,dt_slow,enth_prev_col,NkIce,FIA)
    do j=jsc,jec
      do k=1,ncat ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
        enth_prev(i,j,k) = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,NkIce
          enth_prev(i,j,k) = enth_prev(i,j,k) + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
      endif ; enddo ; enddo

      do i=isc,iec
        heat_in_col(i,j) = heat_in_col(i,j) - FIA%frazil_left(i,j)
        heat_in_col(i,j) = heat_in_col(i,j) - IST%part_size(i,j,0) * dt_slow*FIA%flux_t_top(i,j,0)
      enddo

      do k=1,ncat ; do i=isc,iec
        if (IST%part_size(i,j,k) > 0.0) then
          heat_in_col(i,j) = heat_in_col(i,j) + IST%part_size(i,j,k) * &
              (FIA%tmelt(i,j,k) + FIA%bmelt(i,j,k) - dt_slow*FIA%bheat(i,j))
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
!$OMP                               h2o_change,NkIce,IG,CS,IOF,FIA) &
!$OMP                        private(part_ocn)
  if (CS%ice_rel_salin <= 0.0) then
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
    part_ocn = 0.0
    if (IST%part_size(i,j,0) > CS%ocean_part_min) part_ocn = IST%part_size(i,j,0)

    IOF%flux_t_ocn_top(i,j) = part_ocn * FIA%flux_t_top(i,j,0)
    IOF%flux_q_ocn_top(i,j) = part_ocn * FIA%flux_q_top(i,j,0)
    IOF%flux_lw_ocn_top(i,j) = part_ocn * FIA%flux_lw_top(i,j,0)
    IOF%flux_lh_ocn_top(i,j) = part_ocn * FIA%flux_lh_top(i,j,0)
    IOF%flux_sw_vis_dir_ocn(i,j) = part_ocn * FIA%flux_sw_vis_dir_top(i,j,0)
    IOF%flux_sw_vis_dif_ocn(i,j) = part_ocn * FIA%flux_sw_vis_dif_top(i,j,0)
    IOF%flux_sw_nir_dir_ocn(i,j) = part_ocn * FIA%flux_sw_nir_dir_top(i,j,0)
    IOF%flux_sw_nir_dif_ocn(i,j) = part_ocn * FIA%flux_sw_nir_dif_top(i,j,0)
    IOF%lprec_ocn_top(i,j) = part_ocn * FIA%lprec_top(i,j,0)
    IOF%fprec_ocn_top(i,j) = part_ocn * FIA%fprec_top(i,j,0)
  enddo ; enddo
! mw/new precip will eventually be intercepted by pond eliminating need for next 3 lines
!$OMP do
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    IOF%lprec_ocn_top(i,j) = IOF%lprec_ocn_top(i,j) + &
                             IST%part_size(i,j,k) * FIA%lprec_top(i,j,k)
  enddo ; enddo ; enddo
!$OMP end parallel

  ! Set up temporary tracer array
  npassive = SIS_count_passive_tracers(IST%TrReg)
  if(npassive>0) allocate(TrLay(0:NkIce+1,npassive))

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,IST,S_col0,NkIce,S_col, &
!$OMP                                  dt_slow,snow_to_ice,heat_in,I_NK,enth_units,   &
!$OMP                                  enth_prev,enth_mass_in_col,Idt_slow,bsnk,      &
!$OMP                                  salt_change,net_melt,kg_H_nk,LatHtFus,LatHtVap,&
!$OMP                                  IG,CS,OSS,FIA,IOF,npassive) &
!$OMP                          private(mass_prev,enthalpy,enthalpy_ocean,Salin,     &
!$OMP                                  heat_to_ocn,h2o_ice_to_ocn,h2o_ocn_to_ice,   &
!$OMP                                  evap_from_ocn,salt_to_ice,bablt,enth_evap,   &
!$OMP                                  enth_ice_to_ocn,enth_ocn_to_ice,heat_input,  &
!$OMP                                  heat_mass_in,mass_in,mass_here,enth_here,    &
!$OMP                                  tot_heat_in,enth_imb,mass_imb,norm_enth_imb, &
!$OMP                                  m_lay, mtot_ice, TrLay,                      &
!$OMP                                  I_part,sn2ic,enth_snowfall)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    if (G%mask2dT(i,j) > 0 .and. IST%part_size(i,j,k) > 0) then
      ! reshape the ice based on fluxes
      if (CS%column_check) then
        mass_prev = IST%mH_snow(i,j,k)
        mass_prev = mass_prev + IST%mH_ice(i,j,k)
      endif

 !     evap_from_ocn = 0.0; h2o_ice_to_ocn = 0.0; heat_to_ocn = 0.0

      if (IST%mH_snow(i,j,k) == 0.0) IST%enth_snow(i,j,k,1) = &
          enth_from_TS(Temp_from_En_S(IST%enth_ice(i,j,k,1), S_col0(1), IST%ITV), &
                       0.0, IST%ITV)
      enthalpy(0) = IST%enth_snow(i,j,k,1)
      do m=1,NkIce ; enthalpy(m) = IST%enth_ice(i,j,k,m) ; enddo
      enthalpy_ocean = enthalpy_liquid(OSS%SST_C(i,j), OSS%s_surf(i,j), IST%ITV)
      enthalpy(NkIce+1) = enthalpy_ocean

      ! Handle unpacking and BCs for passive tracers
      call SIS_unpack_passive_ice_tr(i, j, k, nkice, IST%TrReg, TrLay)

      if (CS%ice_rel_salin > 0.0) then
        do m=1,NkIce ; Salin(m) = IST%sal_ice(i,j,k,m) ; enddo
        salin(NkIce+1) = CS%ice_rel_salin * OSS%s_surf(i,j)
      else
        do m=1,NkIce+1 ; Salin(m) = CS%ice_bulk_salin ; enddo
      endif

      m_lay(0) = IST%mH_snow(i,j,k) * IG%H_to_kg_m2
      do m=1,NkIce ; m_lay(m) = IST%mH_ice(i,j,k) * kg_H_Nk ; enddo

      ! mw/new - melt pond size is now adjusted here (rain ignored in resize, for now)
      call ice_resize_SIS2(1-IST%part_size(i,j,0), IST%mH_pond(i,j,k), m_lay, &
                   enthalpy, S_col, Salin, FIA%fprec_top(i,j,k)*dt_slow, &
                   FIA%lprec_top(i,j,k)*dt_slow, FIA%flux_q_top(i,j,k)*dt_slow, &
                   FIA%tmelt(i,j,k), FIA%bmelt(i,j,k), NkIce, npassive, TrLay, &
                   heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                   snow_to_ice(i,j,k), salt_to_ice, IST%ITV, CS%ice_thm_CSp, bablt, &
                   enth_evap, enth_ice_to_ocn, enth_ocn_to_ice)

      IST%mH_snow(i,j,k) = m_lay(0) * IG%kg_m2_to_H
      call rebalance_ice_layers(m_lay, mtot_ice, Enthalpy, Salin, NkIce, npassive, TrLay)
      IST%mH_ice(i,j,k) = mtot_ice * IG%kg_m2_to_H

      if (IST%mH_ice(i,j,k) == 0.0) then
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy_liquid_freeze(S_col(m), IST%ITV) ; enddo
      else
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy(m) ; enddo
      endif
      if (CS%ice_rel_salin > 0.0) then
        if (IST%mH_ice(i,j,k) == 0.0) then
          do m=1,NkIce ; IST%sal_ice(i,j,k,m) = 0.0 ; enddo
        else
          do m=1,NkIce ; IST%sal_ice(i,j,k,m) = Salin(m) ; enddo
        endif
        salt_change(i,j) = salt_change(i,j) + IST%part_size(i,j,k) * salt_to_ice
      endif

      ! Copy back into the tracer array
      call SIS_repack_passive_ice_tr(i, j, k, nkice, IST%TrReg, TrLay)

      ! The snow enthalpy should not have changed.  This should do nothing.
      ! IST%enth_snow(i,j,k,1) = Enthalpy(0)

      enth_snowfall = ((dt_slow*FIA%fprec_top(i,j,k)) * enthalpy(0))
      IOF%Enth_Mass_in_atm(i,j) = IOF%Enth_Mass_in_atm(i,j) + &
           IST%part_size(i,j,k) * enth_snowfall

!      IOF%Enth_Mass_in_ocn(i,j) = IOF%Enth_Mass_in_ocn(i,j) + &
!          IST%part_size(i,j,k) * enth_ocn_to_ice

      IOF%Enth_Mass_in_ocn(i,j) = IOF%Enth_Mass_in_ocn(i,j) + &
          IST%part_size(i,j,k) * (h2o_ocn_to_ice * enthalpy_ocean)

      IOF%Enth_Mass_out_ocn(i,j) = IOF%Enth_Mass_out_ocn(i,j) - &
          IST%part_size(i,j,k) * enth_ice_to_ocn
      IOF%Enth_Mass_out_atm(i,j) = IOF%Enth_Mass_out_atm(i,j) - &
          IST%part_size(i,j,k) * enth_evap


      if (CS%column_check) then
        heat_in(i,j,k) = heat_in(i,j,k) + FIA%tmelt(i,j,k) + FIA%bmelt(i,j,k) - &
                     (heat_to_ocn - (LatHtVap+LatHtFus)*evap_from_ocn)

        heat_input = FIA%tmelt(i,j,k) + FIA%bmelt(i,j,k) - (heat_to_ocn - (LatHtVap+LatHtFus)*evap_from_ocn)
        heat_mass_in = enth_snowfall + enth_ocn_to_ice - enth_ice_to_ocn - enth_evap
        mass_in = dt_slow*FIA%fprec_top(i,j,k) & ! +FIA%lprec_top(i,j,k) <- eventually
                + h2o_ocn_to_ice - h2o_ice_to_ocn &
                - (dt_slow*FIA%flux_q_top(i,j,k)-evap_from_ocn)

        mass_here = IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k) + IST%mH_ice(i,j,k)
        enth_here = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
        tot_heat_in = IG%kg_m2_to_H*(enth_units*heat_input + heat_mass_in)
        mass_in = mass_in*IG%kg_m2_to_H

        enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        mass_imb = mass_here - (mass_prev + mass_in)
        if (abs(enth_imb) > CS%imb_tol * &
            (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
          norm_enth_imb = enth_imb / (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in))
          enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        endif

        enth_mass_in_col(i,j) = enth_mass_in_col(i,j) + IST%part_size(i,j,k) * &
          (enth_snowfall + enth_ocn_to_ice - enth_ice_to_ocn - enth_evap)
      endif

      IOF%flux_q_ocn_top(i,j) = IOF%flux_q_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                  (evap_from_ocn*Idt_slow) ! no ice, evaporation left
      IOF%flux_lh_ocn_top(i,j) = IOF%flux_lh_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                 ((LatHtVap*evap_from_ocn)*Idt_slow)
      IOF%flux_t_ocn_top(i,j) = IOF%flux_t_ocn_top(i,j) + IST%part_size(i,j,k) * &
             (FIA%bheat(i,j) - (heat_to_ocn - LatHtFus*evap_from_ocn)*Idt_slow)
      IOF%flux_sw_vis_dif_ocn(i,j) = IOF%flux_sw_vis_dif_ocn(i,j) + IST%part_size(i,j,k) * &
             (((FIA%flux_sw_vis_dir_top(i,j,k) + FIA%flux_sw_vis_dif_top(i,j,k)) + &
               (FIA%flux_sw_nir_dir_top(i,j,k) + FIA%flux_sw_nir_dif_top(i,j,k))) * &
               FIA%sw_abs_ocn(i,j,k))
      net_melt(i,j) = net_melt(i,j) + IST%part_size(i,j,k) * &
              ((h2o_ice_to_ocn-h2o_ocn_to_ice)*Idt_slow)
      bsnk(i,j) = bsnk(i,j) - IST%part_size(i,j,k)*bablt ! bot. melt. ablation

    endif ! Applying surface fluxes to each category.
  enddo ; enddo ; enddo

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,npassive,G,IG,IST,S_col0,NkIce,S_col, &
!$OMP                                  dt_slow,snow_to_ice,heat_in,I_NK,enth_units,   &
!$OMP                                  enth_prev,enth_mass_in_col,Idt_slow,bsnk,      &
!$OMP                                  salt_change,net_melt,kg_h_Nk,LatHtFus,FIA,CS,OSS,&
!$OMP                                  IOF) &
!$OMP                          private(mass_prev,enthalpy,enthalpy_ocean,Salin,     &
!$OMP                                  heat_to_ocn,h2o_ice_to_ocn,h2o_ocn_to_ice,   &
!$OMP                                  evap_from_ocn,salt_to_ice,bablt,enth_evap,   &
!$OMP                                  enth_ice_to_ocn,enth_ocn_to_ice,heat_input,  &
!$OMP                                  heat_mass_in,mass_in,mass_here,enth_here,    &
!$OMP                                  tot_heat_in,enth_imb,mass_imb,norm_enth_imb, &
!$OMP                                  m_lay,mtot_ice,frazil_cat,k_merge,part_sum,  &
!$OMP                                  fill_frac,d_enth,TrLay,I_part,sn2ic,enth_snowfall)
  do j=jsc,jec ; do i=isc,iec ; if (FIA%frazil_left(i,j)>0.0) then

    frazil_cat(1:ncat) = 0.0
    k_merge = 1  ! Find the category that will be combined with the ice free category.
    if (.not.CS%filling_frazil) then
      do k=1,ncat ; if (IST%part_size(i,j,0)+IST%part_size(i,j,k)>0.01) then
        ! absorb frazil in thinest ice partition available    (was ...>0.0)
        ! raised above threshold from 0 to 0.01 to avert ocean-ice model blow-ups
        k_merge = k ; exit
      endif ; enddo
    endif

    if (IST%part_size(i,j,0) > CS%ocean_part_min) then
      ! Combine the ice-free part size with one of the categories.
      !   Whether or not this is also applied when part_size(i,j,0)==0 changes
      ! answers at roundoff because (t*h)*(1/h) /= t.
      k = k_merge
      I_part = 1.0 / (IST%part_size(i,j,k) + IST%part_size(i,j,0))
      IST%mH_snow(i,j,k) = (IST%mH_snow(i,j,k) * IST%part_size(i,j,k)) * I_part
      IST%mH_ice(i,j,k)  = (IST%mH_ice(i,j,k)  * IST%part_size(i,j,k)) * I_part
      if (allocated(IST%t_surf)) then
        IST%t_surf(i,j,k) = (IST%t_surf(i,j,k) * IST%part_size(i,j,k) + &
                       (T_0degC + OSS%T_fr_ocn(i,j))*IST%part_size(i,j,0)) * I_part
        if (IST%part_size(i,j,k) + IST%part_size(i,j,0) == 0.0) &
          IST%t_surf(i,j,k) = OSS%T_fr_ocn(i,j) + T_0degC
      endif
      IST%part_size(i,j,k) = IST%part_size(i,j,k) + IST%part_size(i,j,0)
      IST%part_size(i,j,0) = CS%ocean_part_min
    endif

    if (CS%filling_frazil) then
      if (CS%fraz_fill_time < 0.0) then
        ! This will apply the frazil uniformly to all categories.
        frazil_cat(ncat) = FIA%frazil_left(i,J)
        FIA%frazil_left(i,j) = 0.0
      else
        part_sum = 0.0
        fill_frac = 1.0 ; if (CS%fraz_fill_time > 0.0) &
          fill_frac = dt_slow / (dt_slow + CS%fraz_fill_time)
        do k=1,ncat-1
          part_sum = part_sum + IST%part_size(i,j,k)
          d_enth = fill_frac * max(0.0, LatHtFus * IG%H_to_kg_m2 * &
                         (IG%mH_cat_bound(k+1) - &
                          max(IST%mH_ice(i,j,k),IG%mH_cat_bound(k))) )
          if (d_enth*part_sum > FIA%frazil_left(i,j)) then
            frazil_cat(k) = FIA%frazil_left(i,j) / part_sum
            FIA%frazil_left(i,j) = 0.0
            exit
          else
            frazil_cat(k) = d_enth
            FIA%frazil_left(i,j) = FIA%frazil_left(i,j) - frazil_cat(k)*part_sum
          endif
        enddo
        if (FIA%frazil_left(i,j) > 0.0) then
          ! Note that at this point we should have that part_sum = 1.0.
          frazil_cat(ncat) = FIA%frazil_left(i,j)
          FIA%frazil_left(i,j) = 0.0
        endif
      endif
      do k=ncat-1,1,-1 ; frazil_cat(k) = frazil_cat(k) + frazil_cat(k+1) ; enddo
    else  ! Not filling frazil.
      ! Set the frazil that is absorbed in this category and remove it from
      ! the overall frazil energy.
      I_part = 1.0 / (IST%part_size(i,j,k_merge))
      frazil_cat(k_merge) = FIA%frazil_left(i,j) * I_part
      FIA%frazil_left(i,j) = 0.0
    endif

    do k=1,ncat ; if (frazil_cat(k) > 0.0) then
      if (CS%column_check) then
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
      enthalpy_ocean = enthalpy_liquid(OSS%SST_C(i,j), OSS%s_surf(i,j), IST%ITV)
      enthalpy(NkIce+1) = enthalpy_ocean

      if (CS%ice_rel_salin > 0.0) then
        do m=1,NkIce ; Salin(m) = IST%sal_ice(i,j,k,m) ; enddo
        salin(NkIce+1) = CS%ice_rel_salin * OSS%s_surf(i,j)
      else
        do m=1,NkIce+1 ; Salin(m) = CS%ice_bulk_salin ; enddo
      endif

      ! Handle unpacking and BCs for passive tracers
      call SIS_unpack_passive_ice_tr(i, j, k, nkice, IST%TrReg, TrLay)

      m_lay(0) = IST%mH_snow(i,j,k) * IG%H_to_kg_m2
      do m=1,NkIce ; m_lay(m) = IST%mH_ice(i,j,k) * kg_H_Nk ; enddo

      call add_frazil_SIS2(m_lay, enthalpy, S_col, Salin, npassive, TrLay, frazil_cat(k), &
                   OSS%T_fr_ocn(i,j), NkIce, h2o_ocn_to_ice, salt_to_ice, IST%ITV, &
                   CS%ice_thm_CSp, Enthalpy_freeze=enth_ocn_to_ice)

      call rebalance_ice_layers(m_lay, mtot_ice, Enthalpy, Salin, NkIce, npassive, TrLay)

      ! Unpack the columns of mass, enthalpy and salinity.
      IST%mH_snow(i,j,k) = m_lay(0) * IG%kg_m2_to_H
      IST%mH_ice(i,j,k) = mtot_ice * IG%kg_m2_to_H

      if (IST%mH_ice(i,j,k) == 0.0) then
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy_liquid_freeze(S_col(m), IST%ITV) ; enddo
      else
        do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enthalpy(m) ; enddo
      endif
      if (CS%ice_rel_salin > 0.0) then
        if (IST%mH_ice(i,j,k) == 0.0) then
          do m=1,NkIce ; IST%sal_ice(i,j,k,m) = 0.0 ; enddo
        else
          do m=1,NkIce ; IST%sal_ice(i,j,k,m) = Salin(m) ; enddo
        endif
        salt_change(i,j) = salt_change(i,j) + IST%part_size(i,j,k) * salt_to_ice
      endif

      ! Copy back into the tracer array
      call SIS_repack_passive_ice_tr(i, j, k, nkice, IST%TrReg, TrLay)

!      IOF%Enth_Mass_in_ocn(i,j) = IOF%Enth_Mass_in_ocn(i,j) + &
!          IST%part_size(i,j,k) * enth_ocn_to_ice
      IOF%Enth_Mass_in_ocn(i,j) = IOF%Enth_Mass_in_ocn(i,j) + &
          IST%part_size(i,j,k) * (h2o_ocn_to_ice * enthalpy_ocean)
      net_melt(i,j) = net_melt(i,j) - &
             (h2o_ocn_to_ice * IST%part_size(i,j,k)) * Idt_slow

      if (CS%column_check) then
        heat_in(i,j,k) = heat_in(i,j,k) - frazil_cat(k)

        enth_here = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
        enth_here = enth_here * IST%part_size(i,j,k)
        tot_heat_in = (enth_units * heat_in(i,j,k) + enth_ocn_to_ice) * IST%part_size(i,j,k)
        enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        if (abs(enth_imb) > CS%imb_tol * &
            (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
          enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        endif

        enth_mass_in_col(i,j) = enth_mass_in_col(i,j) + IST%part_size(i,j,k) * enth_ocn_to_ice
      endif
    endif ; enddo ! frazil_cat>0, k-loop

  endif ; enddo ; enddo   ! frazil>0, i-, and j-loops

  call mpp_clock_end(iceClock5)

  call mpp_clock_begin(iceClock6)
  if (CS%do_ice_limit) then
    qflx_lim_ice(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,NkIce,IST,G,I_enth_units,   &
!$OMP                                  spec_thermo_sal,kg_H_Nk,S_col,Obs_h_ice,dt_slow, &
!$OMP                                  Obs_cn_ice,snow_to_ice,salt_change,qflx_lim_ice, &
!$OMP                                  Idt_slow,net_melt,IG,CS,IOF,FIA,Rho_ice)         &
!$OMP                          private(mtot_ice,frac_keep,frac_melt,salt_to_ice,  &
!$OMP                                  h2o_ice_to_ocn,enth_to_melt,enth_ice_to_ocn,   &
!$OMP                                  ice_melt_lay,snow_melt,enth_freeze)
    do j=jsc,jec ; do i=isc,iec
      mtot_ice = 0.0
      do k=1,ncat
         mtot_ice = mtot_ice + IST%part_size(i,j,k) * IG%H_to_kg_m2 * &
                     (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
      enddo
      if (mtot_ice > CS%max_ice_limit*Rho_ice) then
        frac_keep = CS%max_ice_limit*Rho_ice / mtot_ice
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
        IOF%Enth_Mass_out_ocn(i,j) = IOF%Enth_Mass_out_ocn(i,j) - enth_ice_to_ocn
        if (CS%ice_rel_salin > 0.0) then
          salt_change(i,j) = salt_change(i,j) + IST%part_size(i,j,k) * salt_to_ice
        endif
      endif

    enddo ; enddo
  endif ! End of (CS%do_ice_limit) block
  call mpp_clock_end(iceClock6)

  if (CS%column_check) then
    enth_col(:,:) = 0.0
    ! Add back any frazil that has not been used yet.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,heat_in_col,IST,dt_slow,FIA,IOF)
    do j=jsc,jec ; do i=isc,iec
      heat_in_col(i,j) = heat_in_col(i,j) + FIA%frazil_left(i,j) + IOF%flux_t_ocn_top(i,j)*dt_slow
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
!$OMP                                  heat_in_col,enth_mass_in_col,enth_prev_col,IOF,CS) &
!$OMP                          private(enth_here,tot_heat_in,emic2,tot_heat_in2,   &
!$OMP                                  enth_imb,norm_enth_imb,enth_imb2)
    do j=jsc,jec ; do i=isc,iec
      enth_here = enth_col(i,j)
      tot_heat_in = enth_units*heat_in_col(i,j) + enth_mass_in_col(i,j)
      emic2 = (IOF%Enth_Mass_in_ocn(i,j) + IOF%Enth_Mass_in_atm(i,j) + &
               IOF%Enth_Mass_out_ocn(i,j) + IOF%Enth_Mass_out_atm(i,j))
      tot_heat_in2 = enth_units*heat_in_col(i,j) + emic2

      enth_imb = enth_here - (enth_prev_col(i,j) + tot_heat_in)
      if (abs(enth_imb) > CS%imb_tol * &
          (abs(enth_here) + abs(enth_prev_col(i,j)) + abs(tot_heat_in)) ) then
        norm_enth_imb = enth_imb / (abs(enth_here) + abs(enth_prev_col(i,j)) + abs(tot_heat_in))
        enth_imb = enth_here - (enth_prev_col(i,j) + tot_heat_in)
      endif
      enth_imb2 = enth_here - (enth_prev_col(i,j) + tot_heat_in2)
      if (abs(enth_imb2) > CS%imb_tol * &
          (abs(enth_here) + abs(enth_prev_col(i,j)) + abs(tot_heat_in2)) ) then
        norm_enth_imb = enth_imb2 / (abs(enth_here) + abs(enth_prev_col(i,j)) + abs(tot_heat_in2))
        enth_imb2 = enth_here - (enth_prev_col(i,j) + tot_heat_in2)
      endif
    enddo ; enddo
  endif

  ! Determine the salt fluxes to ocean
  ! Note that at this point salt_change and h2o_change are the negative of the masses.
  if (CS%ice_rel_salin <= 0.0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,G,ncat,salt_change,IST,kg_H_Nk,NkIce)
    do j=jsc,jec ; do m=1,NkIce ; do k=1,ncat ; do i=isc,iec
      salt_change(i,j) = salt_change(i,j) + &
         (IST%sal_ice(i,j,k,m)*(IST%mH_ice(i,j,k)*kg_H_Nk)) * IST%part_size(i,j,k)
    enddo ; enddo ; enddo ; enddo
  endif
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,G,h2o_change, &
!$OMP                                  salt_change,Idt_slow,IG,IOF)
  do j=jsc,jec
    do k=1,ncat ; do i=isc,iec
      h2o_change(i,j) = h2o_change(i,j) + IST%part_size(i,j,k) * &
                        IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
    enddo ; enddo
    do i=isc,iec
      ! Note the conversion here from g m-2 to kg m-2 s-1.
      IOF%flux_salt(i,j) = salt_change(i,j) * (0.001*Idt_slow)
    enddo
  enddo

  !   The remainder of this routine deals with any thermodynamics diagnostic
  ! output that has been requested.
  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)

  yr_dtslow = (864e2*365*Idt_slow)
  if (CS%id_lsnk>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,tmp2d,h2o_change,yr_dtslow)
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = min(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(CS%id_lsnk, tmp2d(isc:iec,jsc:jec), CS%diag)
  endif
  if (CS%id_lsrc>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,tmp2d,h2o_change,yr_dtslow)
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = max(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(CS%id_lsrc, tmp2d(isc:iec,jsc:jec), CS%diag)
  endif
  if (IOF%id_saltf>0) call post_data(IOF%id_saltf, IOF%flux_salt, CS%diag)
  if (CS%id_bsnk>0)  call post_data(CS%id_bsnk, bsnk(isc:iec,jsc:jec)*yr_dtslow, &
                                    CS%diag)
  if (FIA%id_tmelt>0) call post_avg(FIA%id_tmelt, FIA%tmelt, IST%part_size(:,:,1:), CS%diag, G=G, &
                                    scale=Idt_slow, wtd=.true.)
  if (FIA%id_bmelt>0) call post_avg(FIA%id_bmelt, FIA%bmelt, IST%part_size(:,:,1:), CS%diag, G=G, &
                                    scale=Idt_slow, wtd=.true.)
  if (FIA%id_bheat>0) call post_data(FIA%id_bheat, FIA%bheat, CS%diag)
  if (CS%id_sn2ic>0) call post_avg(CS%id_sn2ic, snow_to_ice, IST%part_size(:,:,1:), CS%diag, G=G, &
                                    scale=Idt_slow)
  if (CS%id_qflim>0) call post_data(CS%id_qflim, qflx_lim_ice, CS%diag)
  if (CS%id_qfres>0) call post_data(CS%id_qfres, qflx_res_ice, CS%diag)

  call disable_SIS_averaging(CS%diag)

  ! Combine the liquid precipitation with the net melt of ice and snow for
  ! passing to the ocean. These may later be kept separate.
  do j=jsc,jec ; do i=isc,iec
    IOF%lprec_ocn_top(i,j) = IOF%lprec_ocn_top(i,j) + net_melt(i,j)
  enddo ; enddo

  ! Make sure TrLay is no longer allocated
  if(allocated(TrLay)) deallocate(TrLay)
end subroutine SIS2_thermodynamics


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_slow_thermo_init - initializes the parameters and diagnostics associated
!!    with the SIS_slow_thermo module.
subroutine SIS_slow_thermo_init(Time, G, IG, param_file, diag, CS, tracer_flow_CSp)
  type(time_type),     target, intent(in)    :: Time   ! current time
  type(SIS_hor_grid_type),     intent(in)    :: G      ! The horizontal grid structure
  type(ice_grid_type),         intent(in)    :: IG     ! The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(slow_thermo_CS),        pointer       :: CS
  type(SIS_tracer_flow_control_CS), pointer  :: tracer_flow_CSp

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_slow_thermo" ! This module's name.
  real, parameter    :: missing = -1e34

  call callTree_enter("SIS_slow_thermo_init(), SIS_slow_thermo.F90")

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_slow_thermo_init called with associated control structure.")
!    return
  else
    allocate(CS)
  endif

  CS%diag => diag ; CS%Time => Time
  if (associated(tracer_flow_CSp)) CS%tracer_flow_CSp => tracer_flow_CSp

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, &
     "This module calculates the slow evolution of the ice mass, heat, and salt budgets.")

  call get_param(param_file, mod, "SPECIFIED_ICE", CS%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  call get_param(param_file, mod, "DO_RIDGING", CS%do_ridging, &
                 "If true, apply a ridging scheme to the convergent ice. \n"//&
                 "Otherwise, ice is compressed proportionately if the \n"//&
                 "concentration exceeds 1.  The original SIS2 implementation \n"//&
                 "is based on work by Torge Martin.", default=.false.)
  call get_param(param_file, mod, "ICE_BULK_SALINITY", CS%ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "g/kg", default=4.0)
  call get_param(param_file, mod, "ICE_RELATIVE_SALINITY", CS%ice_rel_salin, &
                 "The initial salinity of sea ice as a fraction of the \n"//&
                 "salinity of the seawater from which it formed.", &
                 units = "nondim", default=0.0)
  if ((CS%ice_bulk_salin > 0.0) .and. (CS%ice_rel_salin > 0.0)) &
    call SIS_error(FATAL, "It is inconsistent to have both ICE_BULK_SALINITY "//&
                   "and ICE_RELATIVE_SALINITY set to positive values.")
  if (CS%ice_bulk_salin < 0.0) CS%ice_bulk_salin = 0.0

  call get_param(param_file, mod, "SIS2_FILLING_FRAZIL", CS%filling_frazil, &
               "If true, apply frazil to fill as many categories as \n"//&
               "possible to fill in a uniform (minimum) amount of ice \n"//&
               "in all the thinnest categories. Otherwise the frazil is \n"//&
               "always assigned to a single category.", default=.true.)
  call get_param(param_file, mod, "FILLING_FRAZIL_TIMESCALE", CS%fraz_fill_time, &
               "A timescale with which the filling frazil causes the \n"//&
               "thinest cells to attain similar thicknesses, or a negative \n"//&
               "number to apply the frazil flux uniformly.", default=0.0, &
               units="s", do_not_log=.not.CS%filling_frazil)
  call get_param(param_file, mod, "MIN_OCEAN_PARTSIZE", CS%ocean_part_min, &
                 "The minimum value for the fractional open-ocean area. \n"//&
                 "This can be 0, but for some purposes it may be useful \n"//&
                 "to set this to a miniscule value (like 1e-40) that will \n"//&
                 "be lost to roundoff during any sums so that the open \n"//&
                 "ocean fluxes can be used in with new categories.", &
                 units="nondim", default=0.0)

  call get_param(param_file, mod, "APPLY_ICE_LIMIT", CS%do_ice_limit, &
                 "If true, restore the sea ice state toward climatology.", &
                 default=.false.)
  call get_param(param_file, mod, "MAX_ICE_THICK_LIMIT", CS%max_ice_limit, &
                 "The maximum permitted sea ice thickness when \n"//&
                 "APPLY_ICE_LIMIT is true.", units="m", default=4.0, &
                 do_not_log=.not.CS%do_ice_limit)

  call get_param(param_file, mod, "DO_ICE_RESTORE", CS%do_ice_restore, &
                 "If true, restore the sea ice state toward climatology.", &
                 default=.false.)
  call get_param(param_file, mod, "ICE_RESTORE_TIMESCALE", CS%ice_restore_timescale, &
                 "The restoring timescale when DO_ICE_RESTORE is true.", &
                 units="days", default=5.0, do_not_log=.not.CS%do_ice_restore)

  call get_param(param_file, mod, "NUDGE_SEA_ICE", CS%nudge_sea_ice, &
                 "If true, constrain the sea ice concentrations using observations.", &
                 default=.false.)
  call get_param(param_file, mod, "NUDGE_SEA_ICE_RATE", CS%nudge_sea_ice_rate, &
                 "The rate of cooling of ice-free water that should be ice \n"//&
                 "covered in order to constrained the ice concentration to \n"//&
                 "track observations.  A suggested value is ~10000 W m-2.", &
                 units = "W m-2", default=0.0, do_not_log=.not.CS%nudge_sea_ice)
  call get_param(param_file, mod, "NUDGE_SEA_ICE_TOLERANCE", CS%nudge_conc_tol, &
                 "The tolerance for mismatch in the sea ice concentations \n"//&
                 "before nudging begins to be applied.  Values of order 0.1\n"//&
                 "should work well.", units = "nondim", default=0.0, &
                 do_not_log=.not.CS%nudge_sea_ice)
  call get_param(param_file, mod, "NUDGE_SEA_ICE_STABILITY", CS%nudge_stab_fac, &
                 "A factor that determines whether the buoyancy flux \n"//&
                 "associated with the sea ice nudging of warm water includes \n"//&
                 "a freshwater flux so as to be destabilizing on net (<1), \n"//&
                 "stabilizing (>1), or neutral (=1).", units="nondim", &
                 default=1.0, do_not_log=.not.CS%nudge_sea_ice)

  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
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

  CS%id_lsrc = register_diag_field('ice_model','LSRC', diag%axesT1, Time, &
               'frozen water local source', 'kg/(m^2*yr)', missing_value=missing)
  CS%id_lsnk = register_diag_field('ice_model','LSNK',diag%axesT1, Time, &
               'frozen water local sink', 'kg/(m^2*yr)', missing_value=missing)
  CS%id_bsnk = register_diag_field('ice_model','BSNK',diag%axesT1, Time, &
               'frozen water local bottom sink', 'kg/(m^2*yr)', missing_value=missing)
  CS%id_sn2ic = register_diag_field('ice_model','SN2IC'  ,diag%axesT1,Time, &
               'rate of snow to ice conversion', 'kg/(m^2*s)', missing_value=missing)

  if (CS%do_ice_restore) then
    CS%id_qfres = register_diag_field('ice_model', 'QFLX_RESTORE_ICE', diag%axesT1, Time, &
                 'Ice Restoring heat flux', 'W/m^2', missing_value=missing)
  endif
  if (CS%do_ice_limit) then
    CS%id_qflim = register_diag_field('ice_model', 'QFLX_LIMIT_ICE', diag%axesT1, Time, &
                 'Ice Limit heat flux', 'W/m^2', missing_value=missing)
  endif
  if (CS%nudge_sea_ice) then
    CS%id_fwnudge  = register_diag_field('ice_model','FW_NUDGE' ,diag%axesT1, Time, &
               'nudging freshwater flux', 'kg/(m^2*s)', missing_value=missing)
  endif

  call SIS2_ice_thm_init(param_file, CS%ice_thm_CSp)

  iceClock7 = mpp_clock_id( '  Ice: slow: conservation check', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock5 = mpp_clock_id( '  Ice: slow: thermodynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock6 = mpp_clock_id( '  Ice: slow: restore/limit', flags=clock_flag_default, grain=CLOCK_LOOP )

  call callTree_leave("SIS_slow_thermo_init()")

end subroutine SIS_slow_thermo_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_slow_thermo_set_ptrs can be used to set one of several pointers that
!! are in the slow_therm_CS.
subroutine SIS_slow_thermo_set_ptrs(CS, transport_CSp, sum_out_CSp)
  type(slow_thermo_CS), pointer :: CS
  type(SIS_transport_CS), optional, pointer :: transport_CSp
  type(SIS_sum_out_CS),   optional, pointer :: sum_out_CSp

  if (present(transport_CSp)) CS%SIS_transport_CSp => transport_CSp
  if (present(sum_out_CSp)) CS%sum_output_CSp => sum_out_CSp

end subroutine SIS_slow_thermo_set_ptrs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_slow_thermo_end deallocates any memory associated with this module.
subroutine SIS_slow_thermo_end (CS)
  type(slow_thermo_CS), pointer :: CS

  call SIS2_ice_thm_end(CS%ice_thm_CSp)

  if (associated(CS)) deallocate(CS)

end subroutine SIS_slow_thermo_end

end module SIS_slow_thermo
