!> Handles the main updates of the ice states at the slower time-scales of the coupling or the
!! interactions with the ocean, including the ice mass balance and related thermodynamics and
!! salinity changes, and thermodynamic coupling with the ocean.  The radiative heating and
!! diffusive temperature changes due to coupling with the atmosphere are handled elsewhere.
module SIS_slow_thermo

! This file is a part of SIS2.  See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   SIS2 is a SEA ICE MODEL for coupling through the GFDL exchange grid. SIS2  !
! is a revision of the original SIS with have extended capabilities, including !
! the option of using a B-grid or C-grid spatial discretization.  The SIS2     !
! software has been extensively reformulated from SIS for greater consistency  !
! with the Modular Ocean Model, version 6 (MOM6), and to permit might tighter  !
! dynamical coupling between the ocean and sea-ice.                            !
!   This module handles the main updates of the ice states at the slower time- !
! scales of the coupling or the interactions with the ocean, including the ice !
! mass balance and related thermodynamics and salinity changes, and            !
! thermodynamic coupling with the ocean.  The radiative heating and diffusive  !
! temperature changes due to coupling with the atmosphere are handled elsewhere.!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!



use ice_grid,          only : ice_grid_type
use ice_spec_mod,      only : get_sea_surface

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_COMPONENT, CLOCK_LOOP, CLOCK_ROUTINE
use MOM_data_override, only : data_override
use MOM_EOS,           only : EOS_type, calculate_density_derivs
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : file_exists, MOM_read_data, slasher
use MOM_time_manager,  only : time_type, time_type_to_real
use MOM_unit_scaling,  only : unit_scale_type

use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_framework,     only : coupler_type_spawn, coupler_type_initialized
use SIS_framework,     only : coupler_type_increment_data, coupler_type_rescale_data
use SIS_framework,     only : coupler_type_send_data
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_optics,        only : VIS_DIR, VIS_DIF, NIR_DIR, NIR_DIF
use SIS_sum_output,    only : SIS_sum_out_CS, write_ice_statistics! , SIS_sum_output_init
use SIS_sum_output,    only : accumulate_bottom_input, accumulate_input_1, accumulate_input_2
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS, SIS_call_tracer_column_fns
use SIS_tracer_registry, only : SIS_unpack_passive_ice_tr, SIS_repack_passive_ice_tr
use SIS_tracer_registry, only : SIS_count_passive_tracers
use SIS_transport,     only : adjust_ice_categories, SIS_transport_CS
use SIS_types,         only : ice_state_type, IST_chksum, IST_bounds_check, total_sfc_flux_type
use SIS_types,         only : ocean_sfc_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_utils,         only : post_avg
use SIS2_ice_thm,      only : SIS2_ice_thm_CS, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,      only : ice_resize_SIS2, add_frazil_SIS2, rebalance_ice_layers
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,      only : enth_from_TS, Temp_from_En_S, enthalpy_liquid, calculate_T_freeze

implicit none ; private

#include <SIS2_memory.h>

public :: slow_thermodynamics, SIS_slow_thermo_init, SIS_slow_thermo_end
public :: slow_thermo_CS, SIS_slow_thermo_set_ptrs

!> The control structure for the SIS slow thermodynamics module
type slow_thermo_CS ; private
  logical :: specified_ice  !< If true, the sea ice is specified and there is
                            !! no need for ice dynamics.
  real :: ice_bulk_salin    !< The globally constant sea ice bulk salinity [gSalt kg-1]
                            !! that is used to calculate the ocean salt flux.
  real :: ice_rel_salin     !< The initial bulk salinity of sea-ice relative to the
                            !! salinity of the water from which it formed [nondim].

  logical :: filling_frazil !< If true, apply frazil to fill as many categories
                            !! as possible to fill in a uniform (minimum) amount
                            !! of frazil in all the thinnest categories.
                            !! Otherwise the frazil is always assigned to a
                            !! single category with part size > 0.01.
  real    :: fraz_fill_time !< A timescale with which the filling frazil causes
                            !! the thinnest cells to attain similar thicknesses,
                            !! or a negative number to apply the frazil flux
                            !! uniformly [T ~> s].

  logical :: do_ridging     !<   If true, apply a ridging scheme to the convergent
                            !! ice.  The original SIS2 implementation is based on
                            !! work by Torge Martin.  Otherwise, ice is compressed
                            !! proportionately if the concentration exceeds 1.

  logical :: do_ice_restore !< If true, restore the sea-ice toward climatology
                            !! by applying a restorative heat flux.
  real    :: ice_restore_timescale !< The time scale for restoring ice when
                            !! do_ice_restore is true [T ~> s].

  logical :: do_ice_limit   !< Limit the sea ice thickness to max_ice_limit.
  real    :: max_ice_limit  !< The maximum sea ice thickness [Z ~> m], when do_ice_limit is true.

  logical :: nudge_sea_ice  !< If true, nudge sea ice concentrations towards observations.
  real    :: nudge_sea_ice_rate !< The rate of cooling of ice-free water that should be ice
                            !! covered in order to constrained the ice concentration to track
                            !! observations [Q R Z T-1 ~> W m-2].  A suggested value is of order 10000 W m-2.
  real    :: nudge_stab_fac !< A factor that determines whether the buoyancy flux associated with
                            !! the sea ice nudging of warm water includes a freshwater flux so as to
                            !! be destabilizing on net (<1), stabilizing (>1), or neutral (=1).
                            !!  The default is 1.
  real    :: nudge_conc_tol !< The tolerance for mismatch in the sea ice concentrations
                            !! before nudging begins to be applied.
  logical :: transmute_ice  !< If true, allow ice to be transmuted directly into seawater with a
                            !! spatially varying rate as a form of outflow open boundary condition.
  real, allocatable, dimension(:,:) :: transmutation_rate  !< A spatially varying rate with which
                            !! sea ice and snow are converted into sea-water [T-1 ~> s-1]

  logical :: debug          !< If true, write verbose checksums for debugging purposes.
  logical :: column_check   !< If true, enable the heat check column by column.
  real    :: imb_tol        !< The tolerance for imbalances to be flagged by column_check [nondim].
  logical :: bounds_check   !< If true, check for sensible values of thicknesses temperatures, fluxes, etc.

  integer :: n_calls = 0    !< The number of times slow_thermodynamics has been called.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.

  type(SIS2_ice_thm_CS), pointer  :: ice_thm_CSp => NULL()
                            !< A pointers to the control structures for a subsidiary module
  type(SIS_transport_CS), pointer :: SIS_transport_CSp => NULL()
                            !< A pointers to the control structures for a subsidiary module
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()
                            !< A pointers to the control structures for a subsidiary module
  type(SIS_tracer_flow_control_CS), pointer :: tracer_flow_CSp => NULL()
                            !< A pointers to the control structures for a subsidiary module

  !>@{ Diagnostic IDs
  integer :: id_qflim=-1, id_qfres=-1, id_fwnudge=-1, id_net_melt=-1, id_CMOR_melt=-1
  integer :: id_lsrc=-1, id_lsnk=-1, id_bsnk=-1, id_sn2ic=-1
  !!@}
end type slow_thermo_CS

!>@{ CPU time clock IDs
integer :: iceClock5, iceClock6, iceClock7
!!@}

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> write out any diagnostics of surface fluxes
subroutine post_flux_diagnostics(IST, FIA, IOF, CS, G, US, IG, Idt_slow)
  type(ice_state_type),      intent(in) :: IST !< A type describing the state of the sea ice
  type(fast_ice_avg_type),   intent(in) :: FIA !< A type containing averages of fields
                                               !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type), intent(in) :: IOF !< A structure containing fluxes from the ice to
                                               !! the ocean that are calculated by the ice model.
  type(slow_thermo_CS),      pointer    :: CS  !< The control structure for the SIS_slow_thermo module
  type(SIS_hor_grid_type),   intent(in) :: G   !< The horizontal grid type
  type(unit_scale_type),     intent(in) :: US  !< A structure with unit conversion factors
  type(ice_grid_type),       intent(in) :: IG  !< The sea-ice specific grid type
  real,                      intent(in) :: Idt_slow !< The inverse of the slow thermodynamic
                                               !! time step [T-1 ~> s-1]

  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: tmp2d, net_sw, sw_dn ! Shortwave fluxes[Q R Z T-1 ~> W m-2]
  real :: sw_cat ! [Q R Z T-1 ~> W m-2]
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin
  real, parameter :: missing = -1e34  ! A missing data fill value
  integer :: i, j, k, m, n, b, nb, isc, iec, jsc, jec, ncat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  nb = size(FIA%flux_sw_top,4)
  ! Flux diagnostics
  !
  if (FIA%id_runoff>0) call post_data(FIA%id_runoff, FIA%runoff, CS%diag)
  if (FIA%id_calving>0) call post_data(FIA%id_calving, FIA%calving_preberg, CS%diag)
  if (FIA%id_runoff_hflx>0) call post_data(FIA%id_runoff_hflx, FIA%runoff_hflx, CS%diag)
  if (FIA%id_calving_hflx>0) call post_data(FIA%id_calving_hflx, FIA%calving_hflx_preberg, CS%diag)
  ! The frazil diagnostic is with the other ocean surface diagnostics.
  ! if (IST%id_frazil>0) call post_data(IST%id_frazil, FIA%frazil_left*Idt_slow, CS%diag)
  if (FIA%id_sh>0) call post_avg(FIA%id_sh, FIA%flux_sh_top, IST%part_size, CS%diag, G=G)
  if (FIA%id_lh>0) call post_avg(FIA%id_lh, FIA%flux_lh_top, IST%part_size, CS%diag, G=G)
  if (FIA%id_evap>0) call post_avg(FIA%id_evap, FIA%evap_top, IST%part_size, CS%diag, G=G)
  if (FIA%id_slp>0) call post_data(FIA%id_slp, FIA%p_atm_surf, CS%diag)

  if ((FIA%id_sw>0) .or. (FIA%id_albedo>0)) then
    !$OMP parallel do default(shared) private(sw_cat)
    do j=jsc,jec
      do i=isc,iec ; net_sw(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        sw_cat = 0 ; do b=1,nb ; sw_cat = sw_cat + FIA%flux_sw_top(i,j,k,b) ; enddo
        net_sw(i,j) = net_sw(i,j) + IST%part_size(i,j,k) * sw_cat
      enddo ; enddo
    enddo
    if (FIA%id_sw>0) call post_data(FIA%id_sw, net_sw, CS%diag)
    if (FIA%id_albedo>0) then
      do j=jsc,jec ; do i=isc,iec
        sw_dn(i,j) = 0.0
        do b=1,size(FIA%flux_sw_dn,3)
          sw_dn(i,j) = sw_dn(i,j) + FIA%flux_sw_dn(i,j,b)
        enddo

        if (G%mask2dT(i,j)<=0.5) then
          tmp2d(i,j) = -1.0 ! This is land.
        elseif ((sw_dn(i,j) > 0.0)) then
          ! The 10.0 below is deliberate.  An albedo of down to -9 can be reported
          ! for detecting inconsistent net_sw and sw_dn.
          tmp2d(i,j) = (sw_dn(i,j) - min(net_sw(i,j), 10.0*sw_dn(i,j))) / &
                       sw_dn(i,j)
        else
          tmp2d(i,j) = 0.0 ! What does the albedo mean at night?
        endif
      enddo ; enddo
      call post_data(FIA%id_albedo, tmp2d, CS%diag)
    endif
  endif
  if (FIA%id_lw>0) call post_avg(FIA%id_lw, FIA%flux_lw_top, IST%part_size, CS%diag, G=G)
  if (FIA%id_snofl>0) call post_avg(FIA%id_snofl, FIA%fprec_top, IST%part_size, CS%diag, G=G)
  if (FIA%id_rain>0) call post_avg(FIA%id_rain, FIA%lprec_top, IST%part_size, CS%diag, G=G)
  if (FIA%id_sw_dn>0) then
    sw_dn(:,:) = 0.0
    do b=1,size(FIA%flux_sw_dn,3) ; do j=jsc,jec ; do i=isc,iec
      sw_dn(i,j) = sw_dn(i,j) + FIA%flux_sw_dn(i,j,b)
    enddo ; enddo ; enddo
    call post_data(FIA%id_sw_dn, sw_dn, CS%diag)
  endif
  if (FIA%id_tsfc>0) call post_data(FIA%id_tsfc, FIA%Tskin_avg, CS%diag)
  if (FIA%id_sitemptop>0) call post_data(FIA%id_sitemptop, FIA%Tskin_avg, CS%diag)
  if (FIA%id_sitemptop_CMOR>0) then
    tmp2d(:,:) = missing
    do j=jsc,jec ; do i=isc,iec
      if (FIA%Tskin_avg(i,j) /= missing) tmp2d(i,j) = FIA%Tskin_avg(i,j) + T_0degC
    enddo ; enddo
    call post_data(FIA%id_sitemptop_CMOR, tmp2d, CS%diag)
  endif

  if (FIA%id_sh0>0) call post_data(FIA%id_sh0, FIA%flux_sh0, CS%diag)
  if (FIA%id_evap0>0) call post_data(FIA%id_evap0, FIA%evap0, CS%diag)
  if (FIA%id_lw0>0) call post_data(FIA%id_lw0, FIA%flux_lw0, CS%diag)
  if (FIA%id_dshdt>0) call post_data(FIA%id_dshdt, FIA%dshdt, CS%diag)
  if (FIA%id_devdt>0) call post_data(FIA%id_devdt, FIA%devapdt, CS%diag)
  if (FIA%id_dlwdt>0) call post_data(FIA%id_dlwdt, FIA%dlwdt, CS%diag)

  if (FIA%id_sh_cat>0) call post_data(FIA%id_sh_cat, FIA%flux_sh_top, CS%diag)
  if (FIA%id_evap_cat>0) call post_data(FIA%id_evap_cat, FIA%evap_top, CS%diag)
  if (FIA%id_lw_cat>0) call post_data(FIA%id_lw_cat, FIA%flux_lw_top, CS%diag)
  if (FIA%id_Tsfc_cat>0) call post_data(FIA%id_Tsfc_cat, FIA%Tskin_cat, CS%diag)

  if (FIA%id_sw_vis>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,tmp2d,IST,FIA)
    do j=jsc,jec
      do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * &
             ( FIA%flux_sw_top(i,j,k,VIS_DIR) + FIA%flux_sw_top(i,j,k,VIS_DIF) )
      enddo ; enddo
    enddo
    call post_data(FIA%id_sw_vis, tmp2d, CS%diag)
  endif
  if (FIA%id_sw_dir>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,tmp2d,IST,FIA)
    do j=jsc,jec
      do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * &
             ( FIA%flux_sw_top(i,j,k,VIS_DIR) + FIA%flux_sw_top(i,j,k,NIR_DIR) )
      enddo ; enddo
    enddo
    call post_data(FIA%id_sw_dir, tmp2d, CS%diag)
  endif
  if (FIA%id_sw_dif>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,tmp2d,IST,FIA)
    do j=jsc,jec
      do i=isc,iec ; tmp2d(i,j) = 0.0 ; enddo
      do k=0,ncat ; do i=isc,iec
        tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k) * &
              ( FIA%flux_sw_top(i,j,k,VIS_DIF) + FIA%flux_sw_top(i,j,k,NIR_DIF) )
      enddo ; enddo
    enddo
    call post_data(FIA%id_sw_dif, tmp2d, CS%diag)
  endif
  if (FIA%id_sw_nir_dir>0) call post_avg(FIA%id_sw_nir_dir, FIA%flux_sw_top(:,:,:,NIR_DIR), &
                             IST%part_size, CS%diag, G=G)
  if (FIA%id_sw_nir_dif>0) call post_avg(FIA%id_sw_nir_dif, FIA%flux_sw_top(:,:,:,NIR_DIF), &
                             IST%part_size, CS%diag, G=G)
  if (FIA%id_sw_vis_dir>0) call post_avg(FIA%id_sw_vis_dir, FIA%flux_sw_top(:,:,:,VIS_DIR), &
                             IST%part_size, CS%diag, G=G)
  if (FIA%id_sw_vis_dif>0) call post_avg(FIA%id_sw_vis_dif, FIA%flux_sw_top(:,:,:,VIS_DIF), &
                             IST%part_size, CS%diag, G=G)

  if (CS%nudge_sea_ice .and. CS%id_fwnudge>0) then
    call post_data(CS%id_fwnudge, IOF%melt_nudge, CS%diag)
  endif

end subroutine post_flux_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> slow_thermodynamics takes care of slow ice thermodynamics and mass changes
subroutine slow_thermodynamics(IST, dt_slow, CS, OSS, FIA, XSF, IOF, G, US, IG)

  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  real,                       intent(in)    :: dt_slow !< The thermodynamic step [T ~> s].
  type(slow_thermo_CS),       pointer       :: CS  !< The control structure for the SIS_slow_thermo module
  type(ocean_sfc_state_type), intent(inout) :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(total_sfc_flux_type),  pointer       :: XSF !< A structure of the excess fluxes between the
                                                   !! atmosphere and the ice or ocean relative to
                                                   !! those stored in TSF
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    h_ice_input    ! The specified ice thickness, with specified_ice [m].

  real :: rho_ice  ! The nominal density of sea ice [R ~> kg m-3].
  real :: Idt_slow ! The inverse of the slow thermodynamic time step [T-1 ~> s-1]
  integer :: i, j, k, l, m, b, nb, isc, iec, jsc, jec, ncat, NkIce
  integer :: isd, ied, jsd, jed

  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    rdg_frac, & ! fraction of ridged ice per category [nondim]
    mi_old      ! Ice mass per unit area before thermodynamics [R Z ~> kg m-2].
  real :: mass_part  ! The mass per unit cell area in a thickness category [R Z ~> kg m-2]

  mi_old(:,:,:) = 0.0
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; NkIce = IG%NkIce
  nb = size(FIA%flux_sw_top,4)
!  I_Nk = 1.0 / NkIce
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  CS%n_calls = CS%n_calls + 1

  if (CS%debug) then
    call IST_chksum("Start update_ice_model_slow", IST, G, US, IG)
  endif

  if (CS%bounds_check) &
    call IST_bounds_check(IST, G, US, IG, "Start of SIS_slow_thermo", OSS=OSS)

  ! Set the frazil heat flux that remains to be applied. This might need
  ! to be moved earlier in the algorithm, if there ever to be multiple calls to
  ! slow_thermodynamics per coupling timestep.
  do j=jsc,jec ; do i=isc,iec
    FIA%frazil_left(i,j) = OSS%frazil(i,j)
  enddo ; enddo

  !
  ! conservation checks: top fluxes
  !
  call cpu_clock_begin(iceClock7)
  call accumulate_input_1(IST, FIA, OSS, dt_slow, G, US, IG, CS%sum_output_CSp)
  if (CS%column_check) &
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
                              message="    SIS_slow_thermo", check_column=.true.)
  call cpu_clock_end(iceClock7)

  !
  ! Thermodynamics
  !
  if (CS%specified_ice) then   ! over-write changes with specifications.
    h_ice_input(:,:) = 0.0
    call get_sea_surface(CS%Time, G%HI, SST=OSS%SST_C, ice_conc=IST%part_size, ice_thick=h_ice_input)
    call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)

    do j=jsc,jec ; do i=isc,iec
      IST%part_size(i,j,0) = 1.0 - IST%part_size(i,j,1)
      IST%mH_ice(i,j,1) = US%m_to_Z*h_ice_input(i,j) * rho_ice

      IOF%flux_sh_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%flux_sh_top(i,j,0)
      IOF%evap_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%evap_top(i,j,0)
      IOF%flux_lw_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%flux_lw_top(i,j,0)
      IOF%flux_lh_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%flux_lh_top(i,j,0)
      do b=1,nb
        IOF%flux_sw_ocn(i,j,b) = IST%part_size(i,j,0) * FIA%flux_sw_top(i,j,0,b)
      enddo
      IOF%lprec_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%lprec_top(i,j,0)
      IOF%fprec_ocn_top(i,j) = IST%part_size(i,j,0) * FIA%fprec_top(i,j,0)
    enddo ; enddo

  endif

  ! IOF must be updated regardless of whether the ice is specified or the
  ! prognostic model is being used.
  call coupler_type_spawn(FIA%tr_flux, IOF%tr_flux_ocn_top, (/isd, isc, iec, ied/), &
                          (/jsd, jsc, jec, jed/), as_needed=.true. )
  call coupler_type_rescale_data(IOF%tr_flux_ocn_top, 0.0)
  call coupler_type_increment_data(FIA%tr_flux, IST%part_size, IOF%tr_flux_ocn_top)

  ! No other thermodynamics need to be done for ice that is specified.
  if (CS%specified_ice) then
    if (associated(XSF)) call add_excess_fluxes(IOF, XSF, G, US)
    return
  endif

  !TOM> Store old ice mass per unit area for calculating partial ice growth.
  mi_old(:,:,:) = IST%mH_ice(:,:,:)

  !TOM> derive ridged ice fraction prior to thermodynamic changes of ice thickness
  !     in order to subtract ice melt proportionally from ridged ice volume (see below)
  if (CS%do_ridging) then
    !$OMP parallel do default(shared) private(mass_part)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      mass_part = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
      rdg_frac(i,j,k) = 0.0 ; if (mass_part > 0.0) &
          rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / mass_part
    enddo ; enddo ; enddo
  endif

  call enable_SIS_averaging(US%T_to_s*dt_slow, CS%Time, CS%diag)

  ! Save out diagnostics of fluxes.  This must go before SIS2_thermodynamics.
  call post_flux_diagnostics(IST, FIA, IOF, CS, G, US, IG, Idt_slow)

  call disable_SIS_averaging(CS%diag)

  call accumulate_input_2(IST, FIA, IOF, OSS, IST%part_size, dt_slow, G, US, IG, CS%sum_output_CSp)
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,IOF)
  do j=jsc,jec ; do i=isc,iec
    IOF%Enth_Mass_in_atm(i,j) = 0.0 ; IOF%Enth_Mass_out_atm(i,j) = 0.0
    IOF%Enth_Mass_in_ocn(i,j) = 0.0 ; IOF%Enth_Mass_out_ocn(i,j) = 0.0
  enddo ; enddo

  ! The thermodynamics routines return updated values of the ice and snow
  ! masses-per-unit area and enthalpies.
  call SIS2_thermodynamics(IST, dt_slow, CS, OSS, FIA, IOF, G, US, IG)

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
  call enable_SIS_averaging(US%T_to_s*dt_slow, CS%Time, CS%diag)
  call SIS_call_tracer_column_fns(dt_slow, G, IG, CS%tracer_flow_CSp, IST%mH_ice, mi_old)
  call disable_SIS_averaging(CS%diag)

  call accumulate_bottom_input(IST, OSS, FIA, IOF, dt_slow, G, US, IG, CS%sum_output_CSp)

  ! This needs to go after accumulate_bottom_input.
  if (associated(XSF)) call add_excess_fluxes(IOF, XSF, G, US)

  if (CS%column_check) &
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
                              message="      Post_thermo A", check_column=.true.)
  call adjust_ice_categories(IST%mH_ice, IST%mH_snow, IST%mH_pond, IST%part_size, &
                             IST%TrReg, G, IG, CS%SIS_transport_CSp) !Niki: add ridging?

  if (CS%column_check) &
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
                              message="      Post_thermo B ", check_column=.true.)

end subroutine slow_thermodynamics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Add in any excess fluxes due to ice state type into the ice_ocean_flux_type
subroutine add_excess_fluxes(IOF, XSF, G, US)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(total_sfc_flux_type), intent(in)    :: XSF !< A structure of the excess fluxes between the
                                                  !! atmosphere and the ice or ocean relative to
                                                  !! those stored in TSF
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors

  real :: sw_comb  ! A combination of two downward shortwave fluxes [Q R Z T-1 ~> W m-2].
  integer :: i, j, k, m, n, b, nb, isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  nb = size(XSF%flux_sw,3)

  do j=jsc,jec ; do i=isc,iec
    IOF%flux_sh_ocn_top(i,j) = IOF%flux_sh_ocn_top(i,j) - XSF%flux_sh(i,j)
    IOF%evap_ocn_top(i,j) = IOF%evap_ocn_top(i,j) - XSF%evap(i,j)
    IOF%flux_lw_ocn_top(i,j) = IOF%flux_lw_ocn_top(i,j) - XSF%flux_lw(i,j)
    IOF%flux_lh_ocn_top(i,j) = IOF%flux_lh_ocn_top(i,j) - XSF%flux_lh(i,j)

    IOF%lprec_ocn_top(i,j) = IOF%lprec_ocn_top(i,j) - XSF%lprec(i,j)
    IOF%fprec_ocn_top(i,j) = IOF%fprec_ocn_top(i,j) - XSF%fprec(i,j)

    call coupler_type_increment_data(XSF%tr_flux, IOF%tr_flux_ocn_top, scale_factor=-1.0)

    ! The shortwave fluxes are more complicated because there are multiple bands
    ! and none of these should have negative fluxes if it can be avoided.
    do b=2,nb,2
      ! Combine the direct and diffuse excess fluxes and convert them into
      ! diffuse fluxes, since the ice scatters any light passing through it.
      IOF%flux_sw_ocn(i,j,b) = IOF%flux_sw_ocn(i,j,b) - (XSF%flux_sw(i,j,b-1) + XSF%flux_sw(i,j,b))

      ! Rearrange shortwave fluxes to prevent any band from having a negative
      ! flux.  (The ocean does not glow.)
      if (IOF%flux_sw_ocn(i,j,b) < 0.0) then
        sw_comb = IOF%flux_sw_ocn(i,j,b) + IOF%flux_sw_ocn(i,j,b-1)
        if (sw_comb >= 0.0) then
          ! Borrow from the direct flux to bring the diffuse flux up to 0.
          IOF%flux_sw_ocn(i,j,b-1) = sw_comb
          IOF%flux_sw_ocn(i,j,b) = 0.0
        elseif ((b==VIS_DIF) .or. (b-1==VIS_DIF)) then
          ! The visible diffuse flux is the total in this band
          IOF%flux_sw_ocn(i,j,b) = sw_comb
          IOF%flux_sw_ocn(i,j,b-1) = 0.0
        else
          ! Borrow from the direct visible flux to bring both the fluxes at this
          ! wavelength up to zero.  (This might be a bad idea that should be revisited.)
          IOF%flux_sw_ocn(i,j,VIS_DIF) = IOF%flux_sw_ocn(i,j,VIS_DIF) + sw_comb
          IOF%flux_sw_ocn(i,j,b) = 0.0
          IOF%flux_sw_ocn(i,j,b-1) = 0.0
        endif
      endif
    enddo
    if (IOF%flux_sw_ocn(i,j,VIS_DIF) < 0.0) then
      ! Borrow from other positive bands if the visible diffuse flux is negative.
      do b=1,nb
        if ((b /= VIS_DIF) .and. (IOF%flux_sw_ocn(i,j,b) > 0.0)) then
          sw_comb = IOF%flux_sw_ocn(i,j,VIS_DIF) + IOF%flux_sw_ocn(i,j,b)
          if (sw_comb < 0.0) then
            IOF%flux_sw_ocn(i,j,VIS_DIF) = sw_comb
            IOF%flux_sw_ocn(i,j,b) = 0.0
          else
            IOF%flux_sw_ocn(i,j,b) = sw_comb
            IOF%flux_sw_ocn(i,j,VIS_DIF) = 0.0
            exit ! Break out of the do b loop.
          endif
        endif
      enddo
    endif
  enddo ; enddo

end subroutine add_excess_fluxes

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS2_thermodynamics does the slow thermodynamic update of the ice state,
!! including freezing or melting, and the accumulation of snow and frazil ice.
subroutine SIS2_thermodynamics(IST, dt_slow, CS, OSS, FIA, IOF, G, US, IG)
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  real,                       intent(in)    :: dt_slow !< The thermodynamic step [T ~> s].
  type(slow_thermo_CS),       pointer       :: CS  !< The control structure for the SIS_slow_thermo module
  type(ocean_sfc_state_type), intent(inout) :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type

  ! This subroutine does the thermodynamic calculations in the same order as SIS1,
  ! but with a greater emphasis on enthalpy as the dominant state variable.

  real, dimension(SZI_(G),SZJ_(G),1:IG%CatIce) :: snow_to_ice ! The flux of snow to the ice [R Z ~> kg m-2]
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: Obs_h_ice ! Observed ice thickness for qflux calculation
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: Obs_cn_ice ! Observed total ice concentration [nondim]
  real, dimension(G%isc:G%iec,G%jsc:G%jec)   :: icec  ! Total ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G))   :: &
    salt_change, &        ! The change in integrated salinity [R Z gSalt kg-1 ~> gSalt m-2]
    h2o_change, &         ! The change in water in the ice [R Z ~> kg m-2]
    bsnk, &               ! The bottom melting mass sink [R Z T-1 ~> kg m-2 s-1]
    tmp2d, &
    qflx_lim_ice, &       ! Ice limiting heat flux [Q R Z T-1 ~> W m-2]
    qflx_res_ice, &       ! Ice restoring heat flux [Q R Z T-1 ~> W m-2]
    cool_nudge, &         ! A heat flux out of the sea ice that
                          ! acts to create sea-ice [Q R Z T-1 ~> W m-2].
    net_melt              ! The net mass flux from the ice and snow into the
                          ! ocean due to melting and freezing integrated
                          ! across all categories [R Z T-1 ~> kg m-2 s-1].
  real, dimension(SZI_(G),SZJ_(G),1:IG%CatIce) :: &
    heat_in, &            ! The input heat [Q R Z ~> J m-2]
    enth_prev, &
    enth
  real, dimension(SZI_(G),SZJ_(G)) :: &
    heat_in_col, &        ! The total heat in each column [Q R Z ~> J m-2]
    enth_prev_col, &
    enth_col, &
    enth_mass_in_col

  real, dimension(IG%NkIce) :: S_col      ! The salinity of a column of ice [gSalt kg-1].
  real, dimension(IG%NkIce+1) :: Salin    ! The conserved bulk salinity of each
                                          ! layer [gSalt kg-1], with the salinity of
                                          ! newly formed ice in layer NkIce+1.
  real, dimension(0:IG%NkIce) :: m_lay    ! The masses of a column of ice and snow [R Z ~> kg m-2].
  real, dimension(0:IG%NkIce) :: Tcol0    ! The temperature of a column of ice and snow [degC].
  real, dimension(0:IG%NkIce) :: S_col0   ! The salinity of a column of ice and snow [gSalt kg-1].
  real, dimension(0:IG%NkIce) :: Tfr_col0 ! The freezing temperature of a column of ice and snow [degC].
  real, dimension(0:IG%NkIce+1) :: &
    enthalpy              ! The initial enthalpy of a column of ice and snow
                          ! and the surface ocean [Q ~> J kg-1].
  real, dimension(IG%CatIce) :: frazil_cat  ! The frazil heating applied to each thickness
                          ! category, averaged over the area of that category [Q R Z ~> J m-2].
  real :: enthalpy_ocean  ! The enthalpy of the ocean surface waters [Q ~> J kg-1].
  real :: heat_fill_val   ! An enthalpy to use for massless categories [Q ~> J kg-1].

  real :: I_part          ! The inverse of a part_size [nondim].
  logical :: spec_thermo_sal  ! If true, use the specified salinities of the
                              ! various sub-layers of the ice for all thermodynamic
                              ! calculations; otherwise use the prognostic
                              ! salinity fields for these calculations.
  real, dimension(:,:), allocatable :: TrLay
                          ! Passive tracer slice through categories
                          ! By default, both the 0 (snow layer) boundary and
                          ! the NkIce+1 (surface ocean layer) are both set to 0
                          ! for all tracers

  type(EOS_type), pointer :: EOS => NULL()
  real :: Cp_water    ! The heat capacity of sea water [Q degC-1 ~> J kg-1 degC-1]
  real :: drho_dT(1), drho_dS(1)
  real :: pres_0(1)
  real :: rho_ice     ! The nominal density of sea ice [R ~> kg m-3].

  real :: Idt_slow    ! The inverse of the thermodynamic step [T-1 ~> s-1].
  real :: yr_dtslow   ! The ratio of 1 year to the thermodynamic time step times some scaling
                      ! factors, used to change the units of several diagnostics to rate yr-1.
  real :: heat_to_ocn    ! The heat passed from the ice to the ocean [Q R Z ~> J m-2]
  real :: water_to_ocn   ! The water passed to the ocean [R Z ~> kg m-2]
  real :: salt_to_ocn    ! The salt passed to the ocean [R Z gSalt kg-1 ~> gSalt m-2]
  real :: heat_from_ice  ! The heat extracted from the ice [Q R Z ~> J m-2]
  real :: water_from_ice ! The water extracted from the ice [R Z ~> kg m-2]
  real :: salt_from_ice  ! The salt extracted from the ice [R Z gSalt kg-1 ~> gSalt m-2]
  real :: ice_loss       ! The loss of ice mass from transmutation [R Z ~> kg m-2]
  real :: snow_loss      ! The loss of snow mass from transmutation [R Z ~> kg m-2]
  real :: h2o_ice_to_ocn ! The downward water flux from the ice to the ocean [R Z ~> kg m-2]
  real :: h2o_ocn_to_ice ! The upward water flux from the ocean to the ice [R Z ~> kg m-2]
  real :: evap_from_ocn  ! The evaporation from the ocean [R Z ~> kg m-2]
  real :: bablt       ! The bottom ablation ice loss [R Z ~> kg m-2]
  real :: salt_to_ice ! The flux of salt from the ocean to the ice [R Z gSalt kg-1 ~> gSalt m-2].
                      ! This may be of either sign; in some places it is an
                      ! average over the whole cell, while in others just a partition.
  real :: mtot_ice    ! The total mass of ice and snow in a cell [R Z ~> kg m-2].
  real :: e2m_tot     ! The total enthalpy required to melt all ice and snow [Q R Z ~> J m-2].
  real :: enth_evap   ! The heat flux associated with sublimation [Q R Z ~> J m-2]
  real :: enth_ice_to_ocn ! The heat flux associated with melting at the ice-ocean interface [Q R Z ~> J m-2]
  real :: enth_ocn_to_ice ! The heat flux associated with freezing at the ice-ocean interface [Q R Z ~> J m-2]
  real :: enth_snowfall ! The heat flux associated with snowfall [Q R Z ~> J m-2]
  real :: tot_heat, heating, tot_frazil, heat_mass_in
  real :: heat_input   ! The input heat [Q R Z ~> J m-2]
  real :: mass_in, mass_here, mass_prev, mass_imb
  real :: frac_keep, frac_melt  ! The fraction of ice and snow to keep or remove [nondim].
  real :: ice_melt_lay ! The amount of excess ice removed from each layer [R Z ~> kg m-2].
  real :: snow_melt    ! The amount of excess snow that is melted [R Z ~> kg m-2].
  real :: enth_freeze  ! The freezing point enthalpy of a layer [Q ~> J kg-1].
  real :: enth_to_melt ! The enthalpy addition required to melt the excess ice
                       ! and snow [Q R Z ~> J m-2].
  real :: I_Nk         ! The inverse of the number of layers in the ice [nondim].
  real :: part_sum     ! A running sum of partition sizes [nondim].
  real :: part_ocn     ! A slightly modified ocean part size [nondim].
  real :: d_enth       ! The change in enthalpy between categories [Q R Z ~> J m-2].
  real :: fill_frac    ! The fraction of the difference between the thicknesses
                       ! in thin categories that will be removed within a single
                       ! timestep with filling_frazil [nondim].
  real :: sw_tot       ! The total shortwave radiation incident on a category [Q R Z T-1 ~> W m-2].
  integer :: i, j, k, l, m, n, b, nb, isc, iec, jsc, jec, ncat, NkIce, tr, npassive
  integer :: k_merge
  real :: LatHtFus     ! The latent heat of fusion of ice [Q ~> J kg-1].
  real :: LatHtVap     ! The latent heat of vaporization of water at 0C [Q ~> J kg-1].
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  real :: tot_heat_in, enth_here, enth_imb, norm_enth_imb, emic2, tot_heat_in2, enth_imb2

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  NkIce = IG%NkIce ; I_Nk = 1.0 / NkIce
  nb = size(FIA%flux_sw_top,4)
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0 / dt_slow

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, &
                   rho_ice=rho_ice, spec_thermo_salin=spec_thermo_sal, &
                   Latent_fusion=LatHtFus, Latent_vapor=LatHtVap)
  S_col0(0) = 0.0 ; do m=1,NkIce ; S_col0(m) = S_col(m) ; enddo
  call calculate_T_Freeze(S_col0(0:NkIce), Tfr_col0(0:NkIce), IST%ITV)

  heat_fill_val = Enth_from_TS(0.0, 0.0, IST%ITV)

  if (.not.spec_thermo_sal) call SIS_error(FATAL, "SIS2_thermodynamics is not "//&
    "prepared for SPECIFIED_THERMO_SALINITY to be false.")

  call cpu_clock_begin(iceClock6)
  ! Add any heat fluxes to restore the sea-ice properties toward a prescribed
  ! state, potentially including freshwater fluxes to avoid driving oceanic
  ! convection.
  if (CS%nudge_sea_ice) then
    if (.not.allocated(IOF%melt_nudge)) allocate(IOF%melt_nudge(isc:iec,jsc:jec))

    cool_nudge(:,:) = 0.0 ; IOF%melt_nudge(:,:) = 0.0

    Obs_cn_ice(:,:) = 0.0
    call data_override('ICE', 'icec', Obs_cn_ice, CS%Time)

    icec(:,:) = 0.0
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
      icec(i,j) = icec(i,j) + IST%part_size(i,j,k)
    enddo ; enddo ; enddo
    pres_0(:) = 0.0
    call get_SIS2_thermo_coefs(IST%ITV, Cp_Water=Cp_water, EOS=EOS)

    do j=jsc,jec ; do i=isc,iec
      if (icec(i,j) < Obs_cn_ice(i,j) - CS%nudge_conc_tol) then
        cool_nudge(i,j) = CS%nudge_sea_ice_rate * &
             ((Obs_cn_ice(i,j)-CS%nudge_conc_tol) - icec(i,j))**2.0 ! W/m2
        if (CS%nudge_stab_fac /= 0.0) then
          if (OSS%SST_C(i,j) > OSS%T_fr_ocn(i,j)) then
            call calculate_density_derivs(OSS%SST_C(i:i,j),OSS%s_surf(i:i,j),pres_0,&
                           drho_dT,drho_dS,1,1,EOS)
            IOF%melt_nudge(i,j) = CS%nudge_stab_fac * (-cool_nudge(i,j)*drho_dT(1)) / &
                                  ((Cp_water*drho_dS(1)) * max(OSS%s_surf(i,j), 1.0) )
          endif
        endif
      elseif (icec(i,j) > Obs_cn_ice(i,j) + CS%nudge_conc_tol) then
        ! Heat the ice but do not apply a fresh water flux.
        cool_nudge(i,j) = -CS%nudge_sea_ice_rate * &
             (icec(i,j) - (Obs_cn_ice(i,j)+CS%nudge_conc_tol))**2.0 ! W/m2
      endif

      if (cool_nudge(i,J) > 0.0) then
        FIA%frazil_left(i,j) = FIA%frazil_left(i,j) + cool_nudge(i,j)*dt_slow
      elseif (cool_nudge(i,J) < 0.0) then
        OSS%bheat(i,j) = OSS%bheat(i,j) - cool_nudge(i,j)
      endif
    enddo ; enddo
  endif
  if (CS%do_ice_restore) then
    ! Get observed ice thickness and concentration for ice restoring, if calculating qflux.
    ! There is no need to apply limits, and only the product is used in these calculations.
    Obs_cn_ice(:,:) = 0.0 ; Obs_h_ice(:,:) = 0.0
    call data_override('ICE', 'sic_obs', Obs_cn_ice, CS%Time)
    call data_override('ICE', 'sit_obs', Obs_h_ice, CS%Time)

    call get_SIS2_thermo_coefs(IST%ITV, Latent_fusion=LatHtFus)
    qflx_res_ice(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      e2m_tot = 0.0
      ! Calculate the average ice and snow enthalpy relative to freezing per unit area.
      do k=1,ncat
        if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)>0.0) then
          e2m_tot = (IST%part_size(i,j,k)*IST%mH_snow(i,j,k)) * &
                       ((enthalpy_liquid_freeze(0.0, IST%ITV) - &
                         IST%enth_snow(i,j,k,1)) )
          if (spec_thermo_sal) then ; do m=1,NkIce
            e2m_tot = e2m_tot + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) * I_Nk) * &
                       ((enthalpy_liquid_freeze(S_col(m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) )
          enddo ; else ; do m=1,NkIce
            e2m_tot = e2m_tot + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) * I_Nk) * &
                       ((enthalpy_liquid_freeze(IST%sal_ice(i,j,k,m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) )
          enddo ; endif
        endif
      enddo
      qflx_res_ice(i,j) = -(LatHtFus*rho_ice*Obs_h_ice(i,j)*Obs_cn_ice(i,j) - e2m_tot) / &
                           CS%ice_restore_timescale
      if (qflx_res_ice(i,j) < 0.0) then
        !There is less ice in model than Obs,
        !so make some ice by increasing frazil heat
        FIA%frazil_left(i,j) = FIA%frazil_left(i,j) - qflx_res_ice(i,j)*dt_slow
        !Note that ice should grow when frazil heat is positive
      elseif (qflx_res_ice(i,j) >  0.0) then
        !There is more ice in model than Obs,
        !so melt ice by increasing heat input to ice from ocean (bheat),
        !        OSS%bheat(i,j) = OSS%bheat(i,j) + qflx_res_ice(i,j)
        !Note that ice should melt when bheat increases.
        !BUT, here it's too late for the bheat to have a negative feedback on the ice thickness
        !since thickness is determined by the melting energies calculated in the fast ice
        !module call ice_temp_SIS2() before this point.
        !So, we should rather change the bottom melt energy directly here
        !(as prescribed in ice_temp_SIS2) to have a restoring effect on the ice thickness
        !later in the call ice_resize_SIS2() in this module.
        do k=1,ncat
          FIA%bmelt(i,j,k) = FIA%bmelt(i,j,k) + dt_slow*qflx_res_ice(i,j)
        enddo
      endif
    enddo ; enddo
  endif
  call cpu_clock_end(iceClock6)


  if (CS%column_check) then
    enth_prev(:,:,:) = 0.0 ; heat_in(:,:,:) = 0.0
    enth_prev_col(:,:) = 0.0 ; heat_in_col(:,:) = 0.0 ; enth_mass_in_col(:,:) = 0.0
    !$OMP parallel do default(shared)
    do j=jsc,jec
      do k=1,ncat ; do i=isc,iec ; if (IST%mH_ice(i,j,k) > 0.0) then
        enth_prev(i,j,k) = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,NkIce
          enth_prev(i,j,k) = enth_prev(i,j,k) + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
      endif ; enddo ; enddo

      do i=isc,iec
        heat_in_col(i,j) = heat_in_col(i,j) - FIA%frazil_left(i,j)
        heat_in_col(i,j) = heat_in_col(i,j) - IST%part_size(i,j,0) * dt_slow*FIA%flux_sh_top(i,j,0)
      enddo

      do k=1,ncat ; do i=isc,iec
        if (IST%part_size(i,j,k) > 0.0) then
          heat_in_col(i,j) = heat_in_col(i,j) + IST%part_size(i,j,k) * &
              (FIA%tmelt(i,j,k) + FIA%bmelt(i,j,k) - dt_slow*OSS%bheat(i,j))
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

  call cpu_clock_begin(iceClock5)

  snow_to_ice(:,:,:) = 0.0 ; net_melt(:,:) = 0.0
  bsnk(:,:) = 0.0
  salt_change(:,:) = 0.0
  h2o_change(:,:) = 0.0
  !$OMP parallel default(shared) private(part_ocn)
  if (CS%ice_rel_salin <= 0.0) then
    !$OMP do
    do j=jsc,jec ; do m=1,NkIce ; do k=1,ncat ; do i=isc,iec
      salt_change(i,j) = salt_change(i,j) - &
         (IST%sal_ice(i,j,k,m)*(IST%mH_ice(i,j,k) * I_Nk)) * IST%part_size(i,j,k)
    enddo ; enddo ; enddo ; enddo
  endif
  !$OMP do
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    h2o_change(i,j) = h2o_change(i,j) - IST%part_size(i,j,k) * &
                      (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
  enddo ; enddo ; enddo

  ! Start accumulating the fluxes at the ocean's surface.
  !$OMP do
  do j=jsc,jec ; do i=isc,iec
    part_ocn = 0.0
    if (IST%part_size(i,j,0) > IG%ocean_part_min) part_ocn = IST%part_size(i,j,0)

    IOF%flux_sh_ocn_top(i,j) = part_ocn * FIA%flux_sh_top(i,j,0)
    IOF%evap_ocn_top(i,j) = part_ocn * FIA%evap_top(i,j,0)
    IOF%flux_lw_ocn_top(i,j) = part_ocn * FIA%flux_lw_top(i,j,0)
    IOF%flux_lh_ocn_top(i,j) = part_ocn * FIA%flux_lh_top(i,j,0)
    do b=1,nb ; IOF%flux_sw_ocn(i,j,b) = part_ocn * FIA%flux_sw_top(i,j,0,b) ; enddo
    IOF%lprec_ocn_top(i,j) = part_ocn * FIA%lprec_top(i,j,0)
    IOF%fprec_ocn_top(i,j) = part_ocn * FIA%fprec_top(i,j,0)
  enddo ; enddo

! mw/new precip will eventually be intercepted by pond eliminating need for next 3 lines
  !$OMP do
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    IOF%lprec_ocn_top(i,j) = IOF%lprec_ocn_top(i,j) + IST%part_size(i,j,k) * FIA%lprec_top(i,j,k)
  enddo ; enddo ; enddo

  ! Add fluxes of snow and other properties to the ocean due to recent ridging or drifting events.
  if (allocated(IST%snow_to_ocn)) then
    !$OMP do
    do j=jsc,jec ; do i=isc,iec ; if (IST%snow_to_ocn(i,j) > 0.0) then
      IOF%fprec_ocn_top(i,j) = IOF%fprec_ocn_top(i,j) + IST%snow_to_ocn(i,j) * Idt_slow
      IOF%Enth_Mass_out_ocn(i,j) = IOF%Enth_Mass_out_ocn(i,j) - &
              IST%snow_to_ocn(i,j) * IST%enth_snow_to_ocn(i,j)
      ! h2o_change(i,j) = h2o_change(i,j) - IST%snow_to_ocn(i,j)
      IST%snow_to_ocn(i,j) = 0.0 ;  IST%enth_snow_to_ocn(i,j) = 0.0
    endif ; enddo ; enddo
  endif
!$OMP end parallel

  ! Set up temporary tracer array
  npassive = SIS_count_passive_tracers(IST%TrReg)
  if(npassive>0) allocate(TrLay(0:NkIce+1,npassive))

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,G,US,IST,S_col0,NkIce,S_col,dt_slow, &
!$OMP                                  snow_to_ice,heat_in,I_NK,enth_prev,enth_mass_in_col,bsnk, &
!$OMP                                  Idt_slow,salt_change,net_melt,LatHtFus,LatHtVap,IG,CS,OSS, &
!$OMP                                  FIA,IOF,npassive,nb) &
!$OMP                          private(mass_prev,enthalpy,enthalpy_ocean,Salin,     &
!$OMP                                  heat_to_ocn,h2o_ice_to_ocn,h2o_ocn_to_ice,   &
!$OMP                                  evap_from_ocn,salt_to_ice,bablt,enth_evap,   &
!$OMP                                  enth_ice_to_ocn,enth_ocn_to_ice,heat_input,  &
!$OMP                                  heat_mass_in,mass_in,mass_here,enth_here,    &
!$OMP                                  tot_heat_in,enth_imb,mass_imb,norm_enth_imb, &
!$OMP                                  m_lay, mtot_ice, TrLay,sw_tot,               &
!$OMP                                  I_part,enth_snowfall)
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

      m_lay(0) = IST%mH_snow(i,j,k)
      do m=1,NkIce ; m_lay(m) = IST%mH_ice(i,j,k) * I_Nk ; enddo

      ! mw/new - melt pond size is now adjusted here (rain ignored in resize, for now)
      call ice_resize_SIS2(1-IST%part_size(i,j,0), IST%mH_pond(i,j,k), m_lay, &
                   enthalpy, S_col, Salin, FIA%fprec_top(i,j,k)*dt_slow, &
                   FIA%lprec_top(i,j,k)*dt_slow, FIA%evap_top(i,j,k)*dt_slow, &
                   FIA%tmelt(i,j,k), FIA%bmelt(i,j,k), NkIce, npassive, TrLay, &
                   heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                   snow_to_ice(i,j,k), salt_to_ice, IST%ITV, US, CS%ice_thm_CSp, bablt, &
                   enth_evap, enth_ice_to_ocn, enth_ocn_to_ice)

      IST%mH_snow(i,j,k) = m_lay(0)
      call rebalance_ice_layers(m_lay, mtot_ice, Enthalpy, Salin, NkIce, npassive, TrLay)
      IST%mH_ice(i,j,k) = mtot_ice

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

      enth_snowfall = (dt_slow*FIA%fprec_top(i,j,k)) * enthalpy(0)
      if (FIA%evap_top(i,j,k) < 0.0) &
        enth_snowfall = enth_snowfall - (dt_slow*FIA%evap_top(i,j,k)) * enthalpy(0)
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

        heat_input = (FIA%tmelt(i,j,k) + FIA%bmelt(i,j,k)) - &
                     (heat_to_ocn - (LatHtVap+LatHtFus)*evap_from_ocn)
        heat_mass_in = (enth_snowfall + enth_ocn_to_ice - enth_ice_to_ocn - enth_evap)
        mass_in = dt_slow*FIA%fprec_top(i,j,k) & ! +FIA%lprec_top(i,j,k) <- eventually
                + h2o_ocn_to_ice - h2o_ice_to_ocn &
                - (dt_slow*FIA%evap_top(i,j,k) - evap_from_ocn)

        mass_here = IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k) + IST%mH_ice(i,j,k)
        enth_here = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
        do m=1,NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
        tot_heat_in = heat_input + heat_mass_in

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

      IOF%evap_ocn_top(i,j) = IOF%evap_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                  (evap_from_ocn*Idt_slow) ! no ice, evaporation left
      IOF%flux_lh_ocn_top(i,j) = IOF%flux_lh_ocn_top(i,j) + IST%part_size(i,j,k) * &
                                 ((LatHtVap*evap_from_ocn)*Idt_slow)
      IOF%flux_sh_ocn_top(i,j) = IOF%flux_sh_ocn_top(i,j) + IST%part_size(i,j,k) * &
             (OSS%bheat(i,j) - (heat_to_ocn - LatHtFus*evap_from_ocn)*Idt_slow)
      sw_tot = 0.0
      do b=2,nb,2 ! This sum combines direct and diffuse fluxes to preserve answers.
        sw_tot = sw_tot + (FIA%flux_sw_top(i,j,k,b-1) + FIA%flux_sw_top(i,j,k,b))
      enddo
      IOF%flux_sw_ocn(i,j,VIS_DIF) = IOF%flux_sw_ocn(i,j,VIS_DIF) + IST%part_size(i,j,k) * &
                                     (sw_tot * FIA%sw_abs_ocn(i,j,k))
      net_melt(i,j) = net_melt(i,j) + IST%part_size(i,j,k) * &
              ((h2o_ice_to_ocn-h2o_ocn_to_ice)*Idt_slow)
      bsnk(i,j) = bsnk(i,j) - IST%part_size(i,j,k)*bablt*Idt_slow ! bot. melt. ablation

    endif ! Applying surface fluxes to each category.
  enddo ; enddo ; enddo

  !$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,npassive,G,US,IG,IST,S_col0,NkIce,  &
  !$OMP                                  S_col,dt_slow,snow_to_ice,heat_in,I_NK,enth_mass_in_col, &
  !$OMP                                  enth_prev,Idt_slow,bsnk,salt_change,net_melt,LatHtFus,   &
  !$OMP                                  FIA,CS,OSS,IOF) &
  !$OMP                          private(mass_prev,enthalpy,enthalpy_ocean,Salin,heat_to_ocn,     &
  !$OMP                                  h2o_ice_to_ocn,h2o_ocn_to_ice,evap_from_ocn,salt_to_ice, &
  !$OMP                                  bablt,enth_evap,enth_ice_to_ocn,enth_ocn_to_ice,         &
  !$OMP                                  heat_input,heat_mass_in,mass_in,mass_here,enth_here,     &
  !$OMP                                  tot_heat_in,enth_imb,mass_imb,norm_enth_imb,m_lay,       &
  !$OMP                                  mtot_ice,frazil_cat,k_merge,part_sum,fill_frac,d_enth,   &
  !$OMP                                  TrLay,I_part,enth_snowfall)
  do j=jsc,jec ; do i=isc,iec ; if (FIA%frazil_left(i,j)>0.0) then

    frazil_cat(1:ncat) = 0.0
    k_merge = 1  ! Find the category that will be combined with the ice free category.
    if (.not.CS%filling_frazil) then
      do k=1,ncat ; if (IST%part_size(i,j,0)+IST%part_size(i,j,k)>0.01) then
        ! absorb frazil in thinnest ice partition available    (was ...>0.0)
        ! raised above threshold from 0 to 0.01 to avert ocean-ice model blow-ups
        k_merge = k ; exit
      endif ; enddo
    endif

    if (IST%part_size(i,j,0) > IG%ocean_part_min) then
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
      IST%part_size(i,j,0) = IG%ocean_part_min
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
          d_enth = fill_frac * max(0.0, LatHtFus * &
                         (IG%mH_cat_bound(k+1) - max(IST%mH_ice(i,j,k),IG%mH_cat_bound(k))) )
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

      m_lay(0) = IST%mH_snow(i,j,k)
      do m=1,NkIce ; m_lay(m) = IST%mH_ice(i,j,k) * I_Nk ; enddo

      call add_frazil_SIS2(m_lay, enthalpy, S_col, Salin, npassive, TrLay, frazil_cat(k), &
                   OSS%T_fr_ocn(i,j), NkIce, h2o_ocn_to_ice, salt_to_ice, IST%ITV, US, &
                   CS%ice_thm_CSp, Enthalpy_freeze=enth_ocn_to_ice)
      call rebalance_ice_layers(m_lay, mtot_ice, Enthalpy, Salin, NkIce, npassive, TrLay)

      ! Unpack the columns of mass, enthalpy and salinity.
      IST%mH_snow(i,j,k) = m_lay(0)
      IST%mH_ice(i,j,k) = mtot_ice

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
        tot_heat_in = (heat_in(i,j,k) + enth_ocn_to_ice) * IST%part_size(i,j,k)
        enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        if (abs(enth_imb) > CS%imb_tol * &
            (abs(enth_here) + abs(enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
          enth_imb = enth_here - (enth_prev(i,j,k) + tot_heat_in)
        endif

        enth_mass_in_col(i,j) = enth_mass_in_col(i,j) + IST%part_size(i,j,k) * enth_ocn_to_ice
      endif
    endif ; enddo ! frazil_cat>0, k-loop

  endif ; enddo ; enddo   ! frazil>0, i-, and j-loops

  call cpu_clock_end(iceClock5)

  call cpu_clock_begin(iceClock6)
  if (CS%do_ice_limit) then
    qflx_lim_ice(:,:) = 0.0
    !$OMP parallel do default(shared) private(mtot_ice,frac_keep,frac_melt,salt_to_ice,  &
    !$OMP                                  h2o_ice_to_ocn,enth_to_melt,enth_ice_to_ocn,   &
    !$OMP                                  ice_melt_lay,snow_melt,enth_freeze)
    do j=jsc,jec ; do i=isc,iec
      mtot_ice = 0.0
      do k=1,ncat
         mtot_ice = mtot_ice + IST%part_size(i,j,k) * &
                     (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
      enddo
      if (mtot_ice > CS%max_ice_limit*rho_ice) then
        frac_keep = CS%max_ice_limit*rho_ice / mtot_ice
        frac_melt = 1.0 - frac_keep
        salt_to_ice = 0.0 ; h2o_ice_to_ocn = 0.0
        enth_to_melt = 0.0 ; enth_ice_to_ocn = 0.0
        do k=1,ncat
          ice_melt_lay = frac_melt * IST%part_size(i,j,k)*IST%mH_ice(i,j,k) * I_Nk
          snow_melt = frac_melt * IST%part_size(i,j,k)*IST%mH_snow(i,j,k)
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
        qflx_lim_ice(i,j) = enth_to_melt * Idt_slow
        IOF%Enth_Mass_out_ocn(i,j) = IOF%Enth_Mass_out_ocn(i,j) - enth_ice_to_ocn
        if (CS%ice_rel_salin > 0.0) then
          salt_change(i,j) = salt_change(i,j) + salt_to_ice
        endif
      endif

    enddo ; enddo
  endif ! End of (CS%do_ice_limit) block
  call cpu_clock_end(iceClock6)

  if (CS%column_check) then
    enth_col(:,:) = 0.0
    ! Add back any frazil that has not been used yet.
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do i=isc,iec
      heat_in_col(i,j) = heat_in_col(i,j) + FIA%frazil_left(i,j) + &
                         IOF%flux_sh_ocn_top(i,j)*dt_slow
    enddo ; enddo
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
      enth_col(i,j) = enth_col(i,j) + &
        (IST%mH_snow(i,j,k)*IST%part_size(i,j,k)) * IST%enth_snow(i,j,k,1)
      do m=1,NkIce
        enth_col(i,j) = enth_col(i,j) + &
          (IST%mH_ice(i,j,k)*IST%part_size(i,j,k)*I_Nk) * IST%enth_ice(i,j,k,m)
      enddo
    endif ; enddo ; enddo ; enddo
    !$OMP parallel do default(shared) private(enth_here,tot_heat_in,emic2,tot_heat_in2, &
    !$OMP                                     enth_imb,norm_enth_imb,enth_imb2)
    do j=jsc,jec ; do i=isc,iec
      enth_here = enth_col(i,j)
      tot_heat_in = heat_in_col(i,j) + enth_mass_in_col(i,j)
      emic2 = (IOF%Enth_Mass_in_ocn(i,j) + IOF%Enth_Mass_in_atm(i,j) + &
               IOF%Enth_Mass_out_ocn(i,j) + IOF%Enth_Mass_out_atm(i,j))
      tot_heat_in2 = heat_in_col(i,j) + emic2

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
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do m=1,NkIce ; do k=1,ncat ; do i=isc,iec
      salt_change(i,j) = salt_change(i,j) + &
         (IST%sal_ice(i,j,k,m)*(IST%mH_ice(i,j,k)*I_Nk)) * IST%part_size(i,j,k)
    enddo ; enddo ; enddo ; enddo
  endif
  !$OMP parallel do default(shared)
  do j=jsc,jec
    do k=1,ncat ; do i=isc,iec
      h2o_change(i,j) = h2o_change(i,j) + IST%part_size(i,j,k) * &
                        (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
    enddo ; enddo
    do i=isc,iec
      ! Note the conversion here from g m-2 to kg m-2 s-1.
      IOF%flux_salt(i,j) = salt_change(i,j) * (0.001*Idt_slow)
    enddo
  enddo

  ! Optionally cause the ice to be converted into sea-water with the ocean properties at a
  ! spatially varying rate by reduction of the part size.  The thicknesses do not change.
  if (CS%transmute_ice .and. allocated(CS%transmutation_rate)) then
    do j=jsc,jec ; do i=isc,iec ; if (CS%transmutation_rate(i,j) > 0.0) then
      salt_from_ice = 0.0 ; heat_from_ice = 0.0 ; water_from_ice = 0.0

      frac_keep = exp(-dt_slow*CS%transmutation_rate(i,j))
      do k=1,ncat ; if (IST%part_size(i,j,k) > 0.0) then
        ice_loss = (1.0 - frac_keep) * IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
        do m=1,NkIce
          salt_from_ice = salt_from_ice + (ice_loss * I_Nk) * IST%sal_ice(i,j,k,m)
          heat_from_ice = heat_from_ice + (ice_loss * I_Nk) * IST%enth_ice(i,j,k,m)
        enddo
        snow_loss = (1.0 - frac_keep) * IST%mH_snow(i,j,k)*IST%part_size(i,j,k)
        heat_from_ice = heat_from_ice + snow_loss * IST%enth_snow(i,j,k,1)

        water_from_ice = ice_loss + snow_loss
        IST%part_size(i,j,k) = frac_keep * IST%part_size(i,j,k)
      endif ; enddo

      water_to_ocn = water_from_ice
      heat_to_ocn = water_from_ice * enthalpy_liquid(OSS%SST_C(i,j), OSS%s_surf(i,j), IST%ITV)
      salt_to_ocn = water_from_ice * OSS%s_surf(i,j)

      ! With transmuted ice, the ice is non-conservatively changed to match the ocean properties.
      IOF%flux_salt(i,j) = IOF%flux_salt(i,j) + salt_to_ocn * (0.001*Idt_slow)
      net_melt(i,j) = net_melt(i,j) + water_to_ocn * Idt_slow  ! This goes to IOF%lprec_ocn_top.
      IOF%Enth_mass_out_ocn(i,j) = IOF%Enth_mass_out_ocn(i,j) + heat_to_ocn

      ! With transmuted ice, the imbalances are stored to close the heat and salt budgets.
      IOF%transmutation_salt_flux(i,j) = (salt_from_ice - salt_to_ocn) * (0.001*Idt_slow)
      IOF%transmutation_enth(i,j) = (heat_from_ice - heat_to_ocn)
    endif ; enddo ; enddo
  endif

  !   The remainder of this routine deals with any thermodynamics diagnostic
  ! output that has been requested.
  call enable_SIS_averaging(US%T_to_s*dt_slow, CS%Time, CS%diag)

  yr_dtslow = US%RZ_T_to_kg_m2s*(864e2*365*Idt_slow)
  if (CS%id_lsnk>0) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = min(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(CS%id_lsnk, tmp2d(isc:iec,jsc:jec), CS%diag)
  endif
  if (CS%id_lsrc>0) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do i=isc,iec
      tmp2d(i,j) = max(h2o_change(i,j),0.0) * yr_dtslow
    enddo ; enddo
    call post_data(CS%id_lsrc, tmp2d(isc:iec,jsc:jec), CS%diag)
  endif
  if (IOF%id_saltf>0) call post_data(IOF%id_saltf, IOF%flux_salt, CS%diag)
  if (CS%id_bsnk>0)  call post_data(CS%id_bsnk, bsnk, CS%diag)
  if (FIA%id_tmelt>0) call post_avg(FIA%id_tmelt, FIA%tmelt, IST%part_size(:,:,1:), CS%diag, G=G, &
                                    scale=Idt_slow, wtd=.true.)
  if (FIA%id_bmelt>0) call post_avg(FIA%id_bmelt, FIA%bmelt, IST%part_size(:,:,1:), CS%diag, G=G, &
                                    scale=Idt_slow, wtd=.true.)
  if (FIA%id_bheat>0) call post_data(FIA%id_bheat, OSS%bheat, CS%diag)
  if (CS%id_sn2ic>0) call post_avg(CS%id_sn2ic, snow_to_ice, IST%part_size(:,:,1:), CS%diag, G=G, &
                                    scale=US%RZ_T_to_kg_m2s*Idt_slow)
  if (CS%id_qflim>0) call post_data(CS%id_qflim, qflx_lim_ice, CS%diag)
  if (CS%id_qfres>0) call post_data(CS%id_qfres, qflx_res_ice, CS%diag)
  if (CS%id_net_melt>0) call post_data(CS%id_net_melt, net_melt, CS%diag)
  if (CS%id_CMOR_melt>0) call post_data(CS%id_CMOR_melt, net_melt, CS%diag)

  if (coupler_type_initialized(IOF%tr_flux_ocn_top)) &
    call coupler_type_send_data(IOF%tr_flux_ocn_top, CS%Time)

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
subroutine SIS_slow_thermo_init(Time, G, US, IG, param_file, diag, CS, tracer_flow_CSp)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid structure
  type(unit_scale_type),       intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),         intent(in)    :: IG   !< The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(slow_thermo_CS),        pointer       :: CS   !< The control structure for the SIS_slow_thermo
                                                     !! module that is initialized here
  type(SIS_tracer_flow_control_CS), &
                               pointer       :: tracer_flow_CSp !< A structure that is used to
                                                     !! orchestrate the calling ice tracer packages

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40) :: mdl = "SIS_slow_thermo" ! This module's name.
  logical           :: debug
  real               :: transmute_scale ! A scaling factor to use when reading the transmutation rate.
  character(len=64)  :: transmute_var   ! Transmutation rate variable name in file
  character(len=200) :: filename, transmute_file, inputdir ! Strings for file/path
  integer :: i, j
  real, parameter   :: missing = -1e34

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
  call log_version(param_file, mdl, version, &
     "This module calculates the slow evolution of the ice mass, heat, and salt budgets.")

  call get_param(param_file, mdl, "SPECIFIED_ICE", CS%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  call get_param(param_file, mdl, "DO_RIDGING", CS%do_ridging, &
                 "If true, apply a ridging scheme to the convergent ice. "//&
                 "Otherwise, ice is compressed proportionately if the "//&
                 "concentration exceeds 1.  The original SIS2 implementation "//&
                 "is based on work by Torge Martin.", default=.false.)
  call get_param(param_file, mdl, "ICE_BULK_SALINITY", CS%ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "g/kg", default=4.0)
  call get_param(param_file, mdl, "ICE_RELATIVE_SALINITY", CS%ice_rel_salin, &
                 "The initial salinity of sea ice as a fraction of the "//&
                 "salinity of the seawater from which it formed.", &
                 units = "nondim", default=0.0)
  if ((CS%ice_bulk_salin > 0.0) .and. (CS%ice_rel_salin > 0.0)) &
    call SIS_error(FATAL, "It is inconsistent to have both ICE_BULK_SALINITY "//&
                   "and ICE_RELATIVE_SALINITY set to positive values.")
  if (CS%ice_bulk_salin < 0.0) CS%ice_bulk_salin = 0.0

  call get_param(param_file, mdl, "SIS2_FILLING_FRAZIL", CS%filling_frazil, &
               "If true, apply frazil to fill as many categories as "//&
               "possible to fill in a uniform (minimum) amount of ice "//&
               "in all the thinnest categories. Otherwise the frazil is "//&
               "always assigned to a single category.", default=.true.)
  call get_param(param_file, mdl, "FILLING_FRAZIL_TIMESCALE", CS%fraz_fill_time, &
               "A timescale with which the filling frazil causes the "//&
               "thinest cells to attain similar thicknesses, or a negative "//&
               "number to apply the frazil flux uniformly.", default=0.0, &
               units="s", scale=US%s_to_T, do_not_log=.not.CS%filling_frazil)

  call get_param(param_file, mdl, "APPLY_ICE_LIMIT", CS%do_ice_limit, &
                 "If true, restore the sea ice state toward climatology.", &
                 default=.false.)
  call get_param(param_file, mdl, "MAX_ICE_THICK_LIMIT", CS%max_ice_limit, &
                 "The maximum permitted sea ice thickness when "//&
                 "APPLY_ICE_LIMIT is true.", units="m", default=4.0, scale=US%m_to_Z, &
                 do_not_log=.not.CS%do_ice_limit)

  call get_param(param_file, mdl, "DO_ICE_RESTORE", CS%do_ice_restore, &
                 "If true, restore the sea ice state toward climatology.", &
                 default=.false.)
  call get_param(param_file, mdl, "ICE_RESTORE_TIMESCALE", CS%ice_restore_timescale, &
                 "The restoring timescale when DO_ICE_RESTORE is true.", &
                 units="days", default=5.0, scale=86400.0*US%s_to_T, do_not_log=.not.CS%do_ice_restore)

  call get_param(param_file, mdl, "NUDGE_SEA_ICE", CS%nudge_sea_ice, &
                 "If true, constrain the sea ice concentrations using observations.", &
                 default=.false.)
  call get_param(param_file, mdl, "NUDGE_SEA_ICE_RATE", CS%nudge_sea_ice_rate, &
                 "The rate of cooling of ice-free water that should be ice "//&
                 "covered in order to constrained the ice concentration to "//&
                 "track observations.  A suggested value is ~10000 W m-2.", &
                 units = "W m-2", default=0.0, scale=US%W_m2_to_QRZ_T, &
                 do_not_log=.not.CS%nudge_sea_ice)
  call get_param(param_file, mdl, "NUDGE_SEA_ICE_TOLERANCE", CS%nudge_conc_tol, &
                 "The tolerance for mismatch in the sea ice concentations "//&
                 "before nudging begins to be applied.  Values of order 0.1 "//&
                 "should work well.", units = "nondim", default=0.0, &
                 do_not_log=.not.CS%nudge_sea_ice)
  call get_param(param_file, mdl, "NUDGE_SEA_ICE_STABILITY", CS%nudge_stab_fac, &
                 "A factor that determines whether the buoyancy flux "//&
                 "associated with the sea ice nudging of warm water includes "//&
                 "a freshwater flux so as to be destabilizing on net (<1), "//&
                 "stabilizing (>1), or neutral (=1).", units="nondim", &
                 default=1.0, do_not_log=.not.CS%nudge_sea_ice)

  call get_param(param_file, mdl, "TRANSMUTE_SEA_ICE", CS%transmute_ice, &
                 "If true, allow ice to be transmuted directly into seawater with a spatially "//&
                 "varying rate as a form of outflow open boundary condition.", default=.false.)
  if (CS%transmute_ice) then

    allocate(CS%transmutation_rate(SZI_(G), SZJ_(G))) ; CS%transmutation_rate(:,:) = 0.0
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    call get_param(param_file, mdl, "TRANSMUTATION_RATE_FILE", transmute_file, &
                 "The file from which the transmutation rate should be read.", &
                 fail_if_missing=.true.)
    filename = trim(slasher(inputdir))//trim(transmute_file)
    call log_param(param_file, mdl, "INPUTDIR/TRANSMUTATION_RATE_FILE", filename)
    call get_param(param_file, mdl, "TRANSMUTATION_RATE_VAR", transmute_var, &
                 "The variable with the map of sea-ice transmutation rate.  No transmutation "//&
                 "occurs where this field is 0.", default="transmute_rate")
    call get_param(param_file, mdl, "TRANSMUTATION_RATE_RESCALE", transmute_scale, &
                 "A rescaling factor to use when reading the transmutation rate to scale it "//&
                 "to units of [s-1].  If this factor is negative, the field that is being read "//&
                 "is a timescale, so the rate is the Adcroft reciprocal of the field that is "//&
                 "read times this scaling factor", default=1.0, units="nondim")

    if (.not.file_exists(filename, G%Domain)) call SIS_error(FATAL, &
       " SIS_slow_thermo_init: Unable to open transmutation rate file "//trim(filename))
    if (transmute_scale >= 0.0) then
      call MOM_read_data(filename, transmute_var, CS%transmutation_rate, G%Domain, scale=transmute_scale*US%T_to_s)
    else ! When (transmute_scale < 0.0) read the scaled timescales, then take their Adcroft reciprocal.
      call MOM_read_data(filename, transmute_var, CS%transmutation_rate, G%Domain, scale=(-transmute_scale)*US%s_to_T)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        if (abs(CS%transmutation_rate(i,j)) > 0.0) &
          CS%transmutation_rate(i,j) = 1.0 / (abs(CS%transmutation_rate(i,j)))
      enddo ; enddo
    endif

  endif

  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_SLOW_ICE", CS%debug, &
                 "If true, write out verbose debugging data on the slow ice PEs.", &
                 default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "COLUMN_CHECK", CS%column_check, &
                 "If true, add code to allow debugging of conservation "//&
                 "column-by-column.  This does not change answers, but "//&
                 "can increase model run time.", default=.false., &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "IMBALANCE_TOLERANCE", CS%imb_tol, &
                 "The tolerance for imbalances to be flagged by COLUMN_CHECK.", &
                 units="nondim", default=1.0e-9, debuggingParam=.true.)
  call get_param(param_file, mdl, "ICE_BOUNDS_CHECK", CS%bounds_check, &
                 "If true, periodically check the values of ice and snow "//&
                 "temperatures and thicknesses to ensure that they are "//&
                 "sensible, and issue warnings if they are not.  This "//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)

  CS%id_lsrc = register_diag_field('ice_model','LSRC', diag%axesT1, Time, &
               'frozen water local source', 'kg/(m^2*yr)', missing_value=missing)
  CS%id_lsnk = register_diag_field('ice_model','LSNK',diag%axesT1, Time, &
               'frozen water local sink', 'kg/(m^2*yr)', missing_value=missing)
  CS%id_bsnk = register_diag_field('ice_model','BSNK',diag%axesT1, Time, &
               'frozen water local bottom sink', &
               'kg/(m^2*yr)', conversion= 864e2*365.*US%RZ_T_to_kg_m2s, &
               missing_value=missing)
  CS%id_sn2ic = register_diag_field('ice_model','SN2IC'  ,diag%axesT1,Time, &
               'rate of snow to ice conversion', 'kg/(m^2*s)', missing_value=missing)
  CS%id_net_melt = register_diag_field('ice_model','net_melt' ,diag%axesT1, Time, &
               'net mass flux from ice & snow to ocean due to melting & freezing', &
               'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, missing_value=missing)
  CS%id_CMOR_melt = register_diag_field('ice_model','fsitherm' ,diag%axesT1, Time, &
               'water_flux_into_sea_water_due_to_sea_ice_thermodynamics', &
               'kg m-2 s-1', conversion=US%RZ_T_to_kg_m2s, missing_value=missing)

  if (CS%do_ice_restore) then
    CS%id_qfres = register_diag_field('ice_model', 'QFLX_RESTORE_ICE', diag%axesT1, Time, &
                 'Ice Restoring heat flux', 'W/m^2', conversion=US%QRZ_T_to_W_m2, missing_value=missing)
  endif
  if (CS%do_ice_limit) then
    CS%id_qflim = register_diag_field('ice_model', 'QFLX_LIMIT_ICE', diag%axesT1, Time, &
                 'Ice Limit heat flux', 'W/m^2', conversion=US%QRZ_T_to_W_m2, missing_value=missing)
  endif
  if (CS%nudge_sea_ice) then
    CS%id_fwnudge  = register_diag_field('ice_model','FW_NUDGE' ,diag%axesT1, Time, &
               'nudging freshwater flux', 'kg/(m^2*s)', conversion=US%RZ_T_to_kg_m2s, missing_value=missing)
  endif

  call SIS2_ice_thm_init(US, param_file, CS%ice_thm_CSp)

  iceClock7 = cpu_clock_id( '  Ice: slow: conservation check', grain=CLOCK_LOOP )
  iceClock5 = cpu_clock_id( '  Ice: slow: thermodynamics', grain=CLOCK_LOOP )
  iceClock6 = cpu_clock_id( '  Ice: slow: restore/limit', grain=CLOCK_LOOP )

  call callTree_leave("SIS_slow_thermo_init()")

end subroutine SIS_slow_thermo_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_slow_thermo_set_ptrs can be used to set one of several pointers that
!! are in the slow_therm_CS.
subroutine SIS_slow_thermo_set_ptrs(CS, transport_CSp, sum_out_CSp)
  type(slow_thermo_CS), pointer :: CS   !< The control structure for the SIS_slow_thermo module
  type(SIS_transport_CS), &
              optional, pointer :: transport_CSp !< This pointer will be set to the control structure
                                        !! for ice transport
  type(SIS_sum_out_CS), &
              optional, pointer :: sum_out_CSp !< This pointer will be set to the control structure
                                        !! for globally summed diagnostics

  if (present(transport_CSp)) CS%SIS_transport_CSp => transport_CSp
  if (present(sum_out_CSp)) CS%sum_output_CSp => sum_out_CSp

end subroutine SIS_slow_thermo_set_ptrs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_slow_thermo_end deallocates any memory associated with this module.
subroutine SIS_slow_thermo_end (CS)
  type(slow_thermo_CS), pointer :: CS   !< The control structure for the SIS_slow_thermo module
                                        !! that is deallocated here

  call SIS2_ice_thm_end(CS%ice_thm_CSp)

  if (associated(CS)) deallocate(CS)

end subroutine SIS_slow_thermo_end

end module SIS_slow_thermo
