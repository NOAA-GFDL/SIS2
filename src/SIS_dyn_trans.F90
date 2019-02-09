!> Handles the main updates of the ice states at the slower time-scales of the couplng or
!! the interactions with the ocean due to ice dynamics and lateral transport.
module SIS_dyn_trans

! This file is part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   SIS2 is a SEA ICE MODEL for coupling through the GFDL exchange grid. SIS2  !
! is a revision of the original SIS with have extended capabilities, including !
! the option of using a B-grid or C-grid spatial discretization.  The SIS2     !
! software has been extensively reformulated from SIS for greater consistency  !
! with the Modular Ocean Model, version 6 (MOM6), and to permit might tighter  !
! dynamical coupling between the ocean and sea-ice.                            !
!   This module handles the main updates of the ice states at the slower time- !
! scales of the couplng or the interactions with the ocean due to ice dynamics !
! and lateral transport.                                                       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : fill_symmetric_edges !, MOM_domains_init, clone_MOM_domain
! use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, log_version, param_file_type
use MOM_hor_index, only : hor_index_type ! , hor_index_init
! use MOM_string_functions, only : uppercase
use MOM_time_manager, only : time_type, time_type_to_real, real_to_time
use MOM_time_manager, only : operator(+), operator(-)
use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)
use MOM_EOS, only : EOS_type, calculate_density_derivs

use coupler_types_mod, only: coupler_type_initialized, coupler_type_send_data
use fms_mod, only : clock_flag_default
! use fms_io_mod, only : restore_state, query_initialized
use fms_io_mod, only : register_restart_field, restart_file_type
use mpp_domains_mod,  only  : domain2D
use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_LOOP, CLOCK_ROUTINE

use SIS_continuity,    only : SIS_continuity_CS, summed_continuity
use SIS_debugging,     only : chksum, Bchksum, hchksum
use SIS_debugging,     only : hchksum_pair, Bchksum_pair, uvchksum
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl, safe_alloc_alloc
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_dyn_bgrid, only : SIS_B_dyn_CS, SIS_B_dynamics, SIS_B_dyn_init
use SIS_dyn_bgrid, only : SIS_B_dyn_register_restarts, SIS_B_dyn_end
use SIS_dyn_cgrid, only : SIS_C_dyn_CS, SIS_C_dynamics, SIS_C_dyn_init
use SIS_dyn_cgrid, only : SIS_C_dyn_register_restarts, SIS_C_dyn_end
use SIS_dyn_cgrid, only : SIS_C_dyn_read_alt_restarts
use SIS_hor_grid,  only : SIS_hor_grid_type
use SIS_sum_output, only : write_ice_statistics, SIS_sum_output_init, SIS_sum_out_CS
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS
use SIS_transport, only : SIS_transport_init, SIS_transport_end
use SIS_transport, only : SIS_transport_CS, adjust_ice_categories, cell_average_state_type
use SIS_transport, only : alloc_cell_average_state_type, dealloc_cell_average_state_type
use SIS_transport, only : cell_ave_state_to_ice_state, ice_state_to_cell_ave_state
use SIS_transport, only : ice_cat_transport, finish_ice_transport
use SIS_types,     only : ocean_sfc_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_types,     only : ice_state_type, IST_chksum, IST_bounds_check
use SIS_utils,     only : get_avg, post_avg, ice_line !, ice_grid_chksum
use SIS2_ice_thm,  only : get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,  only : enth_from_TS, Temp_from_En_S
use slab_ice,      only : slab_ice_advect, slab_ice_dynamics
use ice_bergs,     only : icebergs, icebergs_run, icebergs_init, icebergs_end
use ice_grid,      only : ice_grid_type

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_dynamics_trans, update_icebergs, dyn_trans_CS
public :: specified_ice_dynamics, slab_ice_dyn_trans
public :: SIS_dyn_trans_register_restarts, SIS_dyn_trans_init, SIS_dyn_trans_end
public :: SIS_dyn_trans_read_alt_restarts, stresses_to_stress_mag
public :: SIS_dyn_trans_transport_CS, SIS_dyn_trans_sum_output_CS
public :: post_ocean_sfc_diagnostics, post_ice_state_diagnostics

!> The control structure for the SIS_dyn_trans module
type dyn_trans_CS ; private
  logical :: Cgrid_dyn    !< If true use a C-grid discretization of the sea-ice dynamics.
  real    :: dt_ice_dyn   !< The time step used for the slow ice dynamics, including
                          !! stepping the continuity equation and interactions
                          !! between the ice mass field and velocities [s]. If
                          !! 0 or negative, the coupling time step will be used.
  logical :: merged_cont  !< If true, update the continuity equations for the ice, snow,
                          !! and melt pond water together with proportionate fluxes.
                          !! Otherwise the three media are updated separately.
  logical :: do_ridging   !<   If true, apply a ridging scheme to the convergent
                          !! ice.  The original SIS2 implementation is based on
                          !! work by Torge Martin.  Otherwise, ice is compressed
                          !! proportionately if the concentration exceeds 1.
  integer :: adv_substeps !< The number of advective iterations for each slow time step.
  logical :: berg_windstress_bug  !< If true, use older code that applied an old
                          !! ice-ocean stress to the icebergs in place of the
                          !! current air-ice stress.  This option is here for
                          !! backward compatibility, but should be avoided.
  logical :: Warsaw_sum_order !< If true, use the order of sums in the Warsaw version
                          !! of SIS2.  This option exists for backward compatibilty
                          !! but may eventually be obsoleted.

  logical :: debug        !< If true, write verbose checksums for debugging purposes.
  logical :: column_check !< If true, enable the heat check column by column.
  real    :: imb_tol      !< The tolerance for imbalances to be flagged by
                          !! column_check [nondim].
  logical :: bounds_check !< If true, check for sensible values of thicknesses
                          !! temperatures, fluxes, etc.
  logical :: verbose      !< A flag to control the printing of an ice-diagnostic
                          !! message.  When true, this will slow the model down.

  integer :: max_nts      !< The maximum number of transport steps that can be stored
                          !! before they are carried out.
  integer :: nts = 0      !< The number of accumulated transport steps since the last update.
  integer :: ntrunc = 0   !< The number of times the velocity has been truncated
                          !! since the last call to write_ice_statistics.

  integer :: n_calls = 0  !< The number of times SIS_dynamics_trans has been called.
  type(time_type) :: ice_stats_interval !< The interval between writes of the
                          !! globally summed ice statistics and conservation checks.
  type(time_type) :: write_ice_stats_time !< The next time to write out the ice statistics.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  real, allocatable, dimension(:,:,:) :: mca_step !< The total mass per unit total area of snow
                          !! and ice summed across thickness categories in a cell, before each
                          !! transportation substep [H ~> kg m-2].
  real, allocatable, dimension(:,:,:) :: uh_step !< The total zonal mass fluxes during each
                          !! transportation substep [H m2 s-1 ~> kg s-1].
  real, allocatable, dimension(:,:,:) :: vh_step !< The total meridional mass fluxes during each
                          !! transportation substep [H m2 s-1 ~> kg s-1].

  !>@{ Diagnostic IDs
  integer :: id_fax=-1, id_fay=-1, id_mib=-1, id_mi=-1

  ! These are the diagnostic ids for describing the ice state.
  integer, dimension(:), allocatable :: id_t, id_sal
  integer :: id_cn=-1, id_hi=-1, id_hp=-1, id_hs=-1, id_tsn=-1, id_ext=-1 ! id_hp mw/new
  integer :: id_t_iceav=-1, id_s_iceav=-1, id_e2m=-1

  integer :: id_simass=-1, id_sisnmass=-1, id_sivol=-1
  integer :: id_siconc=-1, id_sithick=-1, id_sisnconc=-1, id_sisnthick=-1
  integer :: id_siu=-1, id_siv=-1, id_sispeed=-1, id_sitimefrac=-1
  !!@}

  type(cell_average_state_type), pointer :: CAS => NULL()
          !< A structure with ocean-cell averaged masses.
  type(SIS_B_dyn_CS), pointer     :: SIS_B_dyn_CSp => NULL()
      !< Pointer to the control structure for the B-grid dynamics module
  type(SIS_C_dyn_CS), pointer     :: SIS_C_dyn_CSp => NULL()
      !< Pointer to the control structure for the C-grid dynamics module
  type(SIS_transport_CS), pointer :: SIS_transport_CSp => NULL()
      !< Pointer to the control structure for the ice transport module
  type(SIS_continuity_CS),    pointer :: continuity_CSp => NULL()
      !< The control structure for the SIS continuity module
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()
     !< Pointer to the control structure for the summed diagnostics module
  logical :: module_is_initialized = .false. !< If true, this module has been initialized.
end type dyn_trans_CS

!>@{ CPU time clock IDs
integer :: iceClock4, iceClock8, iceClock9, iceClocka, iceClockb, iceClockc
!!@}

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> update_icebergs calls icebergs_run and offers diagnostics of some of the
!! iceberg fields that might drive the sea ice or ocean
subroutine update_icebergs(IST, OSS, IOF, FIA, icebergs_CS, dt_slow, G, IG, CS)
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow !< The slow ice dynamics timestep [s].
  type(icebergs),             pointer       :: icebergs_CS !< A control structure for the iceberg model.
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module

  real, dimension(SZI_(G),SZJ_(G))   :: &
    hi_avg            ! The area-weighted average ice thickness [m].
  real, dimension(G%isc:G%iec, G%jsc:G%jec)   :: &
    windstr_x, &      ! The area-weighted average ice thickness [Pa].
    windstr_y         ! The area-weighted average ice thickness [Pa].
  real :: rho_ice     ! The nominal density of sea ice [kg m-3].
  real :: H_to_m_ice  ! The specific volume of ice times the conversion factor
                      ! from thickness units [m H-1 ~> m3].
  integer :: stress_stagger
  integer :: i, j, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  H_to_m_ice = IG%H_to_kg_m2 / rho_ice
  call get_avg(IST%mH_ice, IST%part_size(:,:,1:), hi_avg, wtd=.true.)
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    hi_avg(i,j) = hi_avg(i,j) * H_to_m_Ice
  enddo ; enddo

  if (CS%berg_windstress_bug) then
    ! This code reproduces a long-standing bug, in that the old ice-ocean
    ! stresses are being passed in place of the wind stresses on the icebergs.
    do j=jsc,jec ; do i=isc,iec
      windstr_x(i,j) = IOF%flux_u_ocn(i,j)
      windstr_y(i,j) = IOF%flux_v_ocn(i,j)
    enddo ; enddo
    stress_stagger = IOF%flux_uv_stagger
  else
    do j=jsc,jec ; do i=isc,iec
      windstr_x(i,j) = FIA%WindStr_ocn_x(i,j)
      windstr_y(i,j) = FIA%WindStr_ocn_y(i,j)
    enddo ; enddo
    stress_stagger = AGRID
  endif

  if (IST%Cgrid_dyn) then
    call icebergs_run( icebergs_CS, CS%Time, &
            FIA%calving(isc:iec,jsc:jec), OSS%u_ocn_C(isc-2:iec+1,jsc-1:jec+1), &
            OSS%v_ocn_C(isc-1:iec+1,jsc-2:jec+1), IST%u_ice_C(isc-2:iec+1,jsc-1:jec+1), &
            IST%v_ice_C(isc-1:iec+1,jsc-2:jec+1), windstr_x, windstr_y, &
            OSS%sea_lev(isc-1:iec+1,jsc-1:jec+1), OSS%SST_C(isc:iec,jsc:jec),  &
            FIA%calving_hflx(isc:iec,jsc:jec), FIA%ice_cover(isc-1:iec+1,jsc-1:jec+1), &
            hi_avg(isc-1:iec+1,jsc-1:jec+1), stagger=CGRID_NE, &
            stress_stagger=stress_stagger,sss=OSS%s_surf(isc:iec,jsc:jec), &
            mass_berg=IOF%mass_berg, ustar_berg=IOF%ustar_berg, &
            area_berg=IOF%area_berg )
  else
    call icebergs_run( icebergs_CS, CS%Time, &
            FIA%calving(isc:iec,jsc:jec), OSS%u_ocn_B(isc-1:iec+1,jsc-1:jec+1), &
            OSS%v_ocn_B(isc-1:iec+1,jsc-1:jec+1), IST%u_ice_B(isc-1:iec+1,jsc-1:jec+1), &
            IST%v_ice_B(isc-1:iec+1,jsc-1:jec+1), windstr_x, windstr_y, &
            OSS%sea_lev(isc-1:iec+1,jsc-1:jec+1), OSS%SST_C(isc:iec,jsc:jec),  &
            FIA%calving_hflx(isc:iec,jsc:jec), FIA%ice_cover(isc-1:iec+1,jsc-1:jec+1), &
            hi_avg(isc-1:iec+1,jsc-1:jec+1), stagger=BGRID_NE, &
            stress_stagger=stress_stagger, sss=OSS%s_surf(isc:iec,jsc:jec), &
            mass_berg=IOF%mass_berg, ustar_berg=IOF%ustar_berg, &
            area_berg=IOF%area_berg )
  endif

  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)
  if (IOF%id_ustar_berg>0 .and. associated(IOF%ustar_berg)) then
    call post_data(IOF%id_ustar_berg, IOF%ustar_berg, CS%diag)
  endif
  if (IOF%id_area_berg>0 .and. associated(IOF%area_berg)) then
    call post_data(IOF%id_area_berg, IOF%area_berg, CS%diag)
  endif
  if (IOF%id_mass_berg>0 .and. associated(IOF%mass_berg)) then
    call post_data(IOF%id_mass_berg, IOF%mass_berg, CS%diag)
  endif
  call disable_SIS_averaging(CS%diag)

end subroutine update_icebergs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dynamics_trans makes the calls to do ice dynamics and mass and tracer transport
subroutine SIS_dynamics_trans(IST, OSS, FIA, IOF, dt_slow, CS, icebergs_CS, G, IG, tracer_CSp)

  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow !< The slow ice dynamics timestep [s].
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module
  type(icebergs),             pointer       :: icebergs_CS !< A control structure for the iceberg model.
  type(SIS_tracer_flow_control_CS), pointer :: tracer_CSp !< The structure for controlling calls to
                                                   !! auxiliary ice tracer packages

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    mi_sum, &           ! Masses of ice per unit total area [kg m-2].
    misp_sum, &         ! Combined mass of snow, ice and melt pond water per unit total area [kg m-2].
    ice_free, &         ! The fractional open water [nondim], between 0 & 1.
    ice_cover           ! The fractional ice coverage, summed across all
                        ! thickness categories [nondim], between 0 & 1.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    WindStr_x_B, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_B, &      ! averaged over the ice categories on a B-grid [Pa].
    WindStr_x_ocn_B, &  ! Zonal wind stress on the ice-free ocean on a B-grid [Pa].
    WindStr_y_ocn_B, &  ! Meridional wind stress on the ice-free ocean on a B-grid [Pa].
    str_x_ice_ocn_B, &  ! Zonal ice-ocean stress on a B-grid [Pa].
    str_y_ice_ocn_B     ! Meridional ice-ocean stress on a B-grid [Pa].
  real, dimension(SZIB_(G),SZJ_(G))  :: &
    WindStr_x_Cu, &   ! Zonal wind stress averaged over the ice categores on C-grid u-points [Pa].
    WindStr_x_ocn_Cu, & ! Zonal wind stress on the ice-free ocean on C-grid u-points [Pa].
    str_x_ice_ocn_Cu   ! Zonal ice-ocean stress on C-grid u-points [Pa].
  real, dimension(SZI_(G),SZJB_(G))  :: &
    WindStr_y_Cv, &   ! Meridional wind stress averaged over the ice categores on C-grid v-points [Pa].
    WindStr_y_ocn_Cv, & ! Meridional wind stress on the ice-free ocean on C-grid v-points [Pa].
    str_y_ice_ocn_Cv  ! Meridional ice-ocean stress on C-grid v-points [Pa].

  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBx ! An temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBy ! An temporary array for diagnostics.
  real :: ps_vel   ! The fractional thickness catetory coverage at a velocity point.

  real :: dt_slow_dyn  ! The slow dynamics timestep [s].
  real :: dt_adv       ! The advective timestep [s].
  real :: Idt_slow     ! The inverse of dt_slow [s-1].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rdg_rate  ! A ridging rate [s-1], this will be calculated from the strain rates in the dynamics.
  integer :: i, j, k, n, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed
  integer :: ndyn_steps, nds ! The number of dynamic steps.
  integer :: nts_last ! The number of tracer advection steps before updating IST.

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  CS%n_calls = CS%n_calls + 1
  IOF%stress_count = 0

  ndyn_steps = 1
  if ((CS%dt_ice_dyn > 0.0) .and. (CS%dt_ice_dyn < dt_slow)) &
    ndyn_steps = max(CEILING(dt_slow/CS%dt_ice_dyn - 0.000001), 1)
  dt_slow_dyn = dt_slow / ndyn_steps
  if (CS%merged_cont .and. CS%adv_substeps > 0) dt_adv = dt_slow_dyn / real(CS%adv_substeps)
  nts_last = min(ndyn_steps*CS%adv_substeps, CS%adv_substeps*(CS%max_nts/CS%adv_substeps))

  do nds=1,ndyn_steps

    call enable_SIS_averaging(dt_slow_dyn, CS%Time - real_to_time((ndyn_steps-nds)*dt_slow_dyn), CS%diag)

    call mpp_clock_begin(iceClock4)

    ! Convert the category-resolved ice state into the simplified 2-d ice state.
    ! This should be called after a thermodynamic step or if ice_transport was called.
    if (CS%nts == 0) then
      misp_sum(:,:) = 0.0 ; mi_sum(:,:) = 0.0 ; ice_cover(:,:) = 0.0
      !$OMP parallel do default(shared)
      do j=jsd,jed ; do k=1,ncat ; do i=isd,ied
        misp_sum(i,j) = misp_sum(i,j) + IST%part_size(i,j,k) * &
                        (IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k)))
        mi_sum(i,j) = mi_sum(i,j) + (IG%H_to_kg_m2 * IST%mH_ice(i,j,k))  * IST%part_size(i,j,k)
        ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
      enddo ; enddo ; enddo
      do j=jsd,jed ; do i=isd,ied
        !### This can be merged in above, but it would change answers.
        misp_sum(i,j) = misp_sum(i,j) + mi_sum(i,j)
        ice_free(i,j) = IST%part_size(i,j,0)
      enddo ; enddo

      !  Determine the whole-cell averaged mass of snow and ice.
      call ice_state_to_cell_ave_state(IST, G, IG, CS%SIS_transport_CSp, CS%CAS)
    endif
    if (CS%merged_cont .and. (CS%nts == 0)) then
      do j=jsd,jed ; do i=isd,ied ; CS%mca_step(i,j,1) = misp_sum(i,j) ; enddo ; enddo
    endif
    call mpp_clock_end(iceClock4)

    !
    ! Dynamics - update ice velocities.
    !

    ! In the dynamics code, only the ice velocities are changed, and the ice-ocean
    ! stresses are calculated.  The gravity wave dynamics (i.e. the continuity
    ! equation) are not included in the dynamics.  All of the thickness categories
    ! are merged together.
    if (CS%Cgrid_dyn) then

      call mpp_clock_begin(iceClock4)
      ! Correct the wind stresses for changes in the fractional ice-coverage and set
      ! the wind stresses on the ice and the open ocean for a C-grid staggering.
      ! This block of code must be executed if ice_cover and ice_free or the various wind
      ! stresses were updated.
      call set_wind_stresses_C(FIA, ice_cover, ice_free, WindStr_x_Cu, WindStr_y_Cv, &
                               WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, G, ncat)


      if (CS%debug) then
        call uvchksum("Before SIS_C_dynamics [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G)
        call hchksum(ice_free, "ice_free before SIS_C_dynamics", G%HI)
        call hchksum(misp_sum, "misp_sum before SIS_C_dynamics", G%HI)
        call hchksum(mi_sum, "mi_sum before SIS_C_dynamics", G%HI)
        call hchksum(OSS%sea_lev, "sea_lev before SIS_C_dynamics", G%HI, haloshift=1)
        call hchksum(ice_cover, "ice_cover before SIS_C_dynamics", G%HI, haloshift=1)
        call uvchksum("[uv]_ocn before SIS_C_dynamics", OSS%u_ocn_C, OSS%v_ocn_C, G, halos=1)
        call uvchksum("WindStr_[xy] before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G, halos=1)
!        call hchksum_pair("WindStr_[xy]_A before SIS_C_dynamics", WindStr_x_A, WindStr_y_A, G, halos=1)
      endif

      !### Ridging needs to be added with C-grid dynamics.
      rdg_rate(:,:) = 0.0
      call mpp_clock_begin(iceClocka)
      if (CS%Warsaw_sum_order) then
        call SIS_C_dynamics(1.0-ice_free(:,:), misp_sum, mi_sum, IST%u_ice_C, IST%v_ice_C, &
                            OSS%u_ocn_C, OSS%v_ocn_C, WindStr_x_Cu, WindStr_y_Cv, OSS%sea_lev, &
                            str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, dt_slow_dyn, G, CS%SIS_C_dyn_CSp)
      else
        call SIS_C_dynamics(ice_cover, misp_sum, mi_sum, IST%u_ice_C, IST%v_ice_C, &
                            OSS%u_ocn_C, OSS%v_ocn_C, WindStr_x_Cu, WindStr_y_Cv, OSS%sea_lev, &
                            str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, dt_slow_dyn, G, CS%SIS_C_dyn_CSp)
      endif
      call mpp_clock_end(iceClocka)

      if (CS%debug) call uvchksum("After ice_dynamics [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G)

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
      call pass_vector(str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, G%Domain, stagger=CGRID_NE)
      call mpp_clock_end(iceClockb)

      ! Dynamics diagnostics
      call mpp_clock_begin(iceClockc)
      if (CS%id_fax>0) call post_data(CS%id_fax, WindStr_x_Cu, CS%diag)
      if (CS%id_fay>0) call post_data(CS%id_fay, WindStr_y_Cv, CS%diag)

      if (CS%debug) call uvchksum("Before set_ocean_top_stress_Cgrid [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G)

      ! Store all mechanical ocean forcing.
      if (CS%Warsaw_sum_order) then
        call set_ocean_top_stress_Cgrid(IOF, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, &
                                        str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, IST%part_size, G, IG)
      else
        call set_ocean_top_stress_C2(IOF, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, &
                                     str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, ice_free, ice_cover, G)
      endif

      call mpp_clock_end(iceClockc)

      call mpp_clock_end(iceClock4)

    else ! B-grid dynamics.

      ! Correct the wind stresses for changes in the fractional ice-coverage and set
      ! the wind stresses on the ice and the open ocean for a C-grid staggering.
      ! This block of code must be executed if ice_cover and ice_free or the various wind
      ! stresses were updated.

      call mpp_clock_begin(iceClock4)
      call set_wind_stresses_B(FIA, ice_cover, ice_free, WindStr_x_B, WindStr_y_B, &
                               WindStr_x_ocn_B, WindStr_y_ocn_B, G, ncat)

      if (CS%debug) then
        call Bchksum_pair("[uv]_ice_B before dynamics", IST%u_ice_B, IST%v_ice_B, G)
        call hchksum(ice_free, "ice_free before ice_dynamics", G%HI)
        call hchksum(misp_sum, "misp_sum before ice_dynamics", G%HI)
        call hchksum(mi_sum, "mi_sum before ice_dynamics", G%HI)
        call hchksum(OSS%sea_lev, "sea_lev before ice_dynamics", G%HI, haloshift=1)
        call Bchksum_pair("[uv]_ocn before ice_dynamics", OSS%u_ocn_B, OSS%v_ocn_B, G)
        call Bchksum_pair("WindStr_[xy]_B before ice_dynamics", WindStr_x_B, WindStr_y_B, G, halos=1)
      endif

      rdg_rate(:,:) = 0.0
      call mpp_clock_begin(iceClocka)
      if (CS%Warsaw_sum_order) then
        call SIS_B_dynamics(1.0-ice_free(:,:), misp_sum, mi_sum, IST%u_ice_B, IST%v_ice_B, &
                            OSS%u_ocn_B, OSS%v_ocn_B, WindStr_x_B, WindStr_y_B, OSS%sea_lev, &
                            str_x_ice_ocn_B, str_y_ice_ocn_B, CS%do_ridging, &
                            rdg_rate(isc:iec,jsc:jec), dt_slow_dyn, G, CS%SIS_B_dyn_CSp)
      else
        call SIS_B_dynamics(ice_cover, misp_sum, mi_sum, IST%u_ice_B, IST%v_ice_B, &
                            OSS%u_ocn_B, OSS%v_ocn_B, WindStr_x_B, WindStr_y_B, OSS%sea_lev, &
                            str_x_ice_ocn_B, str_y_ice_ocn_B, CS%do_ridging, &
                            rdg_rate(isc:iec,jsc:jec), dt_slow_dyn, G, CS%SIS_B_dyn_CSp)
      endif
      call mpp_clock_end(iceClocka)

      if (CS%debug) call Bchksum_pair("After dynamics [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G)

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
      call mpp_clock_end(iceClockb)

      ! Dynamics diagnostics
      call mpp_clock_begin(iceClockc)
      if ((CS%id_fax>0) .or. (CS%id_fay>0)) then
        !$OMP parallel do default(shared) private(ps_vel)
        do J=jsc-1,jec ; do I=isc-1,iec
          ps_vel = (1.0 - G%mask2dBu(I,J)) + 0.25*G%mask2dBu(I,J) * &
                ((ice_free(i+1,j+1) + ice_free(i,j)) + &
                 (ice_free(i+1,j) + ice_free(i,j+1)) )
          diagVarBx(I,J) = ps_vel * WindStr_x_ocn_B(I,J) + (1.0-ps_vel) * WindStr_x_B(I,J)
          diagVarBy(I,J) = ps_vel * WindStr_y_ocn_B(I,J) + (1.0-ps_vel) * WindStr_y_B(I,J)
        enddo ; enddo

        if (CS%id_fax>0) call post_data(CS%id_fax, diagVarBx, CS%diag)
        if (CS%id_fay>0) call post_data(CS%id_fay, diagVarBy, CS%diag)
      endif

      if (CS%debug) call Bchksum_pair("Before set_ocean_top_stress_Bgrid [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G)
      ! Store all mechanical ocean forcing.
      if (CS%Warsaw_sum_order) then
        call set_ocean_top_stress_Bgrid(IOF, WindStr_x_ocn_B, WindStr_y_ocn_B, &
                                        str_x_ice_ocn_B, str_y_ice_ocn_B, IST%part_size, G, IG)
      else
        call set_ocean_top_stress_B2(IOF, WindStr_x_ocn_B, WindStr_y_ocn_B, &
                                        str_x_ice_ocn_B, str_y_ice_ocn_B, ice_free, ice_cover, G)
      endif
      call mpp_clock_end(iceClockc)

      ! Convert the velocities to C-grid points for use in transport.
      do j=jsc,jec ; do I=isc-1,iec
        IST%u_ice_C(I,j) = 0.5 * ( IST%u_ice_B(I,J-1) + IST%u_ice_B(I,J) )
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        IST%v_ice_C(i,J) = 0.5 * ( IST%v_ice_B(I-1,J) + IST%v_ice_B(I,J) )
      enddo ; enddo
      call mpp_clock_end(iceClock4)
    endif ! End of B-grid dynamics

    ! Do ice mass transport and related tracer transport.  This updates the category-decomposed ice state.
    call mpp_clock_begin(iceClock8)
    if (CS%debug) call uvchksum("Before ice_transport [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G)
    call enable_SIS_averaging(dt_slow_dyn, CS%Time - real_to_time((ndyn_steps-nds)*dt_slow_dyn), CS%diag)

    if (CS%merged_cont) then
      if (CS%nts+CS%adv_substeps > CS%max_nts) call SIS_error(FATAL, &
        "Attempting to store more advective substeps than allocated space allows.  Increase MAX_NTS.")
      do n = CS%nts+1, CS%nts+CS%adv_substeps
        if (n < nts_last) then
          ! Some of the work is not needed for the last step before cat_ice_transport.
          call summed_continuity(IST%u_ice_C, IST%v_ice_C, CS%mca_step(:,:,n), CS%mca_step(:,:,n+1), &
                                 CS%uh_step(:,:,n), CS%vh_step(:,:,n), dt_adv, G, IG, CS%continuity_CSp)
          call pass_var(CS%mca_step(:,:,n+1), G%Domain)
        else
          call summed_continuity(IST%u_ice_C, IST%v_ice_C, CS%mca_step(:,:,n), CS%mca_step(:,:,n+1), &
                                 CS%uh_step(:,:,n), CS%vh_step(:,:,n), dt_adv, G, IG, CS%continuity_CSp)
        endif
      enddo
      CS%nts = CS%nts + CS%adv_substeps
    else
      call ice_cat_transport(CS%CAS, IST%TrReg, dt_slow_dyn, CS%adv_substeps, G, IG, CS%SIS_transport_CSp, &
                             uc=IST%u_ice_C, vc=IST%v_ice_C)
    endif

    if (CS%merged_cont .and. ((CS%nts + CS%adv_substeps > CS%max_nts) .or. (nds==ndyn_steps))) then
      if (CS%nts /= nts_last) call SIS_error(FATAL, "Bad logic in calculating nts_last.")
      n = CS%nts
      call ice_cat_transport(CS%CAS, IST%TrReg, dt_slow_dyn, CS%nts, G, IG, &
                             CS%SIS_transport_CSp, mca_tot=CS%mca_step(:,:,1:n+1), &
                             uh_tot=CS%uh_step(:,:,1:n), vh_tot=CS%vh_step(:,:,1:n))
      CS%nts = 0
    endif

    if (CS%nts==0) &
      call finish_ice_transport(CS%CAS, IST, IST%TrReg, G, IG, CS%SIS_transport_CSp, &
                                rdg_rate=rdg_rate)

    call mpp_clock_end(iceClock8)

    if (CS%column_check .and. (CS%nts==0)) &
      call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
                                message="      Post_transport")! , check_column=.true.)

  enddo ! nds=1,ndyn_steps
  call finish_ocean_top_stresses(IOF, G)

  ! Set appropriate surface quantities in categories with no ice.
  if (allocated(IST%t_surf)) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)<=0.0) &
      IST%t_surf(i,j,k) = T_0degC + OSS%T_fr_ocn(i,j)
    enddo ; enddo ; enddo
  endif

  ! Calculate and output various diagnostics of the ice state.
  call mpp_clock_begin(iceClock9)

  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)
  call post_ice_state_diagnostics(CS, IST, OSS, IOF, dt_slow, CS%Time, G, IG, CS%diag)
  call disable_SIS_averaging(CS%diag)

  if (CS%verbose) call ice_line(CS%Time, IST%part_size(isc:iec,jsc:jec,0), OSS%SST_C(:,:), G)
  if (CS%debug) call IST_chksum("End SIS_dynamics_trans", IST, G, IG)
  if (CS%bounds_check) call IST_bounds_check(IST, G, IG, "End of SIS_dynamics_trans", OSS=OSS)

  if (CS%Time + real_to_time(0.5*dt_slow) > CS%write_ice_stats_time) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
                              tracer_CSp=tracer_CSp)
    CS%write_ice_stats_time = CS%write_ice_stats_time + CS%ice_stats_interval
  elseif (CS%column_check) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp)
  endif

  call mpp_clock_end(iceClock9)

end subroutine SIS_dynamics_trans


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> slab_ice_dynamics_trans makes the calls to do the slab ice version of dynamics and mass and tracer transport
subroutine slab_ice_dyn_trans(IST, OSS, FIA, IOF, dt_slow, CS, G, IG, tracer_CSp)

  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow !< The slow ice dynamics timestep [s].
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module
  type(SIS_tracer_flow_control_CS), pointer :: tracer_CSp !< The structure for controlling calls to
                                                   !! auxiliary ice tracer packages

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    mi_sum, &           ! Masses of ice per unit total area [kg m-2].
    misp_sum            ! Combined mass of snow, ice and melt pond water per unit total area [kg m-2].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    WindStr_x_B, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_B, &      ! averaged over the ice categories on a B-grid [Pa].
    WindStr_x_ocn_B, &  ! Zonal wind stress on the ice-free ocean on a B-grid [Pa].
    WindStr_y_ocn_B, &  ! Meridional wind stress on the ice-free ocean on a B-grid [Pa].
    str_x_ice_ocn_B, &  ! Zonal ice-ocean stress on a B-grid [Pa].
    str_y_ice_ocn_B     ! Meridional ice-ocean stress on a B-grid [Pa].
  real, dimension(SZIB_(G),SZJ_(G))  :: &
    WindStr_x_Cu, &   ! Zonal wind stress averaged over the ice categores on C-grid u-points [Pa].
    WindStr_x_ocn_Cu, & ! Zonal wind stress on the ice-free ocean on C-grid u-points [Pa].
    str_x_ice_ocn_Cu   ! Zonal ice-ocean stress on C-grid u-points [Pa].
  real, dimension(SZI_(G),SZJB_(G))  :: &
    WindStr_y_Cv, &   ! Meridional wind stress averaged over the ice categores on C-grid v-points [Pa].
    WindStr_y_ocn_Cv, & ! Meridional wind stress on the ice-free ocean on C-grid v-points [Pa].
    str_y_ice_ocn_Cv  ! Meridional ice-ocean stress on C-grid v-points [Pa].

  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBx ! An temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBy ! An temporary array for diagnostics.
  real :: ps_vel   ! The fractional thickness catetory coverage at a velocity point.

  real :: dt_slow_dyn  ! The slow dynamics timestep [s].
  real :: Idt_slow     ! The inverse of dt_slow [s-1].
  integer :: i, j, k, n, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed
  integer :: ndyn_steps, nds ! The number of dynamic steps.

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = 1
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  CS%n_calls = CS%n_calls + 1
  IOF%stress_count = 0

  ndyn_steps = 1
  if ((CS%dt_ice_dyn > 0.0) .and. (CS%dt_ice_dyn < dt_slow)) &
    ndyn_steps = max(CEILING(dt_slow/CS%dt_ice_dyn - 0.000001), 1)
  dt_slow_dyn = dt_slow / ndyn_steps

  do nds=1,ndyn_steps

    call enable_SIS_averaging(dt_slow_dyn, CS%Time - real_to_time((ndyn_steps-nds)*dt_slow_dyn), CS%diag)

    call mpp_clock_begin(iceClock4)
    !$OMP parallel do default(shared)
    do j=jsd,jed ; do i=isd,ied
      mi_sum(i,j) = (IG%H_to_kg_m2 * IST%mH_ice(i,j,1))  * IST%part_size(i,j,1)
      misp_sum(i,j) = mi_sum(i,j) + IST%part_size(i,j,1) * &
                      (IG%H_to_kg_m2 * (IST%mH_snow(i,j,1) + IST%mH_pond(i,j,1)))
    enddo ; enddo
    call mpp_clock_end(iceClock4)

    !
    ! Dynamics - update ice velocities.
    !

    ! In the dynamics code, only the ice velocities are changed, and the ice-ocean
    ! stresses are calculated.  The gravity wave dynamics (i.e. the continuity
    ! equation) are not included in the dynamics.  All of the thickness categories
    ! are merged together.
    if (CS%Cgrid_dyn) then

      call mpp_clock_begin(iceClock4)
      ! Correct the wind stresses for changes in the fractional ice-coverage and set
      ! the wind stresses on the ice and the open ocean for a C-grid staggering.
      ! This block of code must be executed if ice_cover and ice_free or the various wind
      ! stresses were updated.
      call set_wind_stresses_C(FIA, IST%part_size(:,:,1), IST%part_size(:,:,0), WindStr_x_Cu, WindStr_y_Cv, &
                               WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, G, 1)

      if (CS%debug) then
        call uvchksum("Before SIS_C_dynamics [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G)
        call hchksum(IST%part_size(:,:,0), "ice_free before SIS_C_dynamics", G%HI)
        call hchksum(misp_sum, "misp_sum before SIS_C_dynamics", G%HI)
        call hchksum(mi_sum, "mi_sum before SIS_C_dynamics", G%HI)
        call hchksum(OSS%sea_lev, "sea_lev before SIS_C_dynamics", G%HI, haloshift=1)
        call hchksum(IST%part_size(:,:,1), "ice_cover before SIS_C_dynamics", G%HI, haloshift=1)
        call uvchksum("[uv]_ocn before SIS_C_dynamics", OSS%u_ocn_C, OSS%v_ocn_C, G, halos=1)
        call uvchksum("WindStr_[xy] before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G, halos=1)
!        call hchksum_pair("WindStr_[xy]_A before SIS_C_dynamics", WindStr_x_A, WindStr_y_A, G, halos=1)
      endif

      call mpp_clock_begin(iceClocka)
      call slab_ice_dynamics(IST%u_ice_C, IST%v_ice_C, OSS%u_ocn_C, OSS%v_ocn_C, &
                             WindStr_x_Cu, WindStr_y_Cv, str_x_ice_ocn_Cu, str_y_ice_ocn_Cv)
      call mpp_clock_end(iceClocka)

      if (CS%debug) call uvchksum("After ice_dynamics [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G)

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
      call pass_vector(str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, G%Domain, stagger=CGRID_NE)
      call mpp_clock_end(iceClockb)

      ! Dynamics diagnostics
      call mpp_clock_begin(iceClockc)
      if (CS%id_fax>0) call post_data(CS%id_fax, WindStr_x_Cu, CS%diag)
      if (CS%id_fay>0) call post_data(CS%id_fay, WindStr_y_Cv, CS%diag)

      if (CS%debug) call uvchksum("Before set_ocean_top_stress_Cgrid [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G)

      ! Store all mechanical ocean forcing.
      call set_ocean_top_stress_C2(IOF, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, &
                                   IST%part_size(:,:,0), IST%part_size(:,:,1), G)
      call mpp_clock_end(iceClockc)

      call mpp_clock_end(iceClock4)

    else ! B-grid dynamics.

      call mpp_clock_begin(iceClock4)
      ! Correct the wind stresses for changes in the fractional ice-coverage and set
      ! the wind stresses on the ice and the open ocean for a C-grid staggering.
      ! This block of code must be executed if ice_cover and ice_free or the various wind
      ! stresses were updated.

      call set_wind_stresses_B(FIA, IST%part_size(:,:,1), IST%part_size(:,:,0), WindStr_x_B, WindStr_y_B, &
                               WindStr_x_ocn_B, WindStr_y_ocn_B, G, 1)

      if (CS%debug) then
        call Bchksum_pair("[uv]_ice_B before dynamics", IST%u_ice_B, IST%v_ice_B, G)
        call hchksum(IST%part_size(:,:,0), "ice_free before ice_dynamics", G%HI)
        call hchksum(misp_sum, "misp_sum before ice_dynamics", G%HI)
        call hchksum(mi_sum, "mi_sum before ice_dynamics", G%HI)
        call hchksum(OSS%sea_lev, "sea_lev before ice_dynamics", G%HI, haloshift=1)
        call Bchksum_pair("[uv]_ocn before ice_dynamics", OSS%u_ocn_B, OSS%v_ocn_B, G)
        call Bchksum_pair("WindStr_[xy]_B before ice_dynamics", WindStr_x_B, WindStr_y_B, G, halos=1)
      endif

      call mpp_clock_begin(iceClocka)
      call slab_ice_dynamics(IST%u_ice_B, IST%v_ice_B, OSS%u_ocn_B, OSS%v_ocn_B, &
                             WindStr_x_B, WindStr_y_B, str_x_ice_ocn_B, str_y_ice_ocn_B)
      call mpp_clock_end(iceClocka)

      if (CS%debug) call Bchksum_pair("After dynamics [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G)

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
      call mpp_clock_end(iceClockb)

      ! Dynamics diagnostics
      call mpp_clock_begin(iceClockc)
      if ((CS%id_fax>0) .or. (CS%id_fay>0)) then
        !$OMP parallel do default(shared) private(ps_vel)
        do J=jsc-1,jec ; do I=isc-1,iec
          ps_vel = (1.0 - G%mask2dBu(I,J)) + 0.25*G%mask2dBu(I,J) * &
                ((IST%part_size(i+1,j+1,0) + IST%part_size(i,j,0)) + &
                 (IST%part_size(i+1,j,0) + IST%part_size(i,j+1,0)) )
          diagVarBx(I,J) = ps_vel * WindStr_x_ocn_B(I,J) + (1.0-ps_vel) * WindStr_x_B(I,J)
          diagVarBy(I,J) = ps_vel * WindStr_y_ocn_B(I,J) + (1.0-ps_vel) * WindStr_y_B(I,J)
        enddo ; enddo

        if (CS%id_fax>0) call post_data(CS%id_fax, diagVarBx, CS%diag)
        if (CS%id_fay>0) call post_data(CS%id_fay, diagVarBy, CS%diag)
      endif

      if (CS%debug) call Bchksum_pair("Before set_ocean_top_stress_Bgrid [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G)
      ! Store all mechanical ocean forcing.
      call set_ocean_top_stress_B2(IOF, WindStr_x_ocn_B, WindStr_y_ocn_B, str_x_ice_ocn_B, str_y_ice_ocn_B, &
                                   IST%part_size(:,:,0), IST%part_size(:,:,1), G)
      call mpp_clock_end(iceClockc)

       ! Convert the B-grid velocities to C-grid points for transport.
      if (CS%debug) call Bchksum_pair("Before ice_transport [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G)
      do j=jsc,jec ; do I=isc-1,iec
        IST%u_ice_C(I,j) = 0.5 * ( IST%u_ice_B(I,J-1) + IST%u_ice_B(I,J) )
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        IST%v_ice_C(i,J) = 0.5 * ( IST%v_ice_B(I-1,J) + IST%v_ice_B(I,J) )
      enddo ; enddo

      call mpp_clock_end(iceClock4)

    endif ! End of B-grid dynamics

    ! Do ice mass transport and related tracer transport.  This updates the category-decomposed ice state.
    call mpp_clock_begin(iceClock8)
    if (CS%debug) call uvchksum("Before ice_transport [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G)
    call enable_SIS_averaging(dt_slow_dyn, CS%Time - real_to_time((ndyn_steps-nds)*dt_slow_dyn), CS%diag)

    call slab_ice_advect(IST%u_ice_C, IST%v_ice_C, IST%mH_ice(:,:,1), 4.0*IG%kg_m2_to_H, &
                         dt_slow_dyn, G, IST%part_size(:,:,1), nsteps=CS%adv_substeps)
    call mpp_clock_end(iceClock8)

    if (CS%column_check .and. (CS%nts==0)) &
      call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
                                message="      Post_transport")! , check_column=.true.)

  enddo ! nds=1,ndyn_steps
  call finish_ocean_top_stresses(IOF, G)

  ! Set appropriate surface quantities in categories with no ice.
  if (allocated(IST%t_surf)) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do i=isc,iec ; if (IST%part_size(i,j,1)<=0.0) &
      IST%t_surf(i,j,1) = T_0degC + OSS%T_fr_ocn(i,j)
    enddo ; enddo
  endif

  ! Calculate and output various diagnostics of the ice state.
  call mpp_clock_begin(iceClock9)

  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)
  call post_ice_state_diagnostics(CS, IST, OSS, IOF, dt_slow, CS%Time, G, IG, CS%diag)
  call disable_SIS_averaging(CS%diag)

  if (CS%verbose) call ice_line(CS%Time, IST%part_size(isc:iec,jsc:jec,0), OSS%SST_C(:,:), G)
  if (CS%debug) call IST_chksum("End slab_ice_dyn_trans", IST, G, IG)
  if (CS%bounds_check) call IST_bounds_check(IST, G, IG, "End of slab_ice_dyn_trans", OSS=OSS)

  if (CS%Time + real_to_time(0.5*dt_slow) > CS%write_ice_stats_time) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
                              tracer_CSp=tracer_CSp)
    CS%write_ice_stats_time = CS%write_ice_stats_time + CS%ice_stats_interval
  elseif (CS%column_check) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp)
  endif

  call mpp_clock_end(iceClock9)

end subroutine slab_ice_dyn_trans


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> specified_ice_dynamics does an update of ice dynamic quantities with specified ice.
subroutine specified_ice_dynamics(IST, OSS, FIA, IOF, dt_slow, CS, G, IG)

  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow !< The slow ice dynamics timestep [s].
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module

  ! Local variables
  integer :: i, j, k, isc, iec, jsc, jec, ncat

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  CS%n_calls = CS%n_calls + 1

  IOF%stress_count = 0
  call set_ocean_top_stress_FIA(FIA, IOF, G)
  call finish_ocean_top_stresses(IOF, G)

  ! Set appropriate surface quantities in categories with no ice.
  if (allocated(IST%t_surf)) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)<=0.0) &
      IST%t_surf(i,j,k) = T_0degC + OSS%T_fr_ocn(i,j)
    enddo ; enddo ; enddo
  endif

  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)
  call post_ice_state_diagnostics(CS, IST, OSS, IOF, dt_slow, CS%Time, G, IG, CS%diag)
  call disable_SIS_averaging(CS%diag)

  if (CS%debug) call IST_chksum("End SIS_dynamics_trans", IST, G, IG)
  if (CS%bounds_check) call IST_bounds_check(IST, G, IG, "End of SIS_dynamics_trans", OSS=OSS)

  if (CS%Time + real_to_time(0.5*dt_slow) > CS%write_ice_stats_time) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp)
    CS%write_ice_stats_time = CS%write_ice_stats_time + CS%ice_stats_interval
  endif

end subroutine specified_ice_dynamics


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Offer diagnostics of the slowly evolving sea ice state.
subroutine post_ice_state_diagnostics(CS, IST, OSS, IOF, dt_slow, Time, G, IG, diag)
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
! type(fast_ice_avg_type),   intent (inout) :: FIA ! A type containing averages of fields
                                                   ! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(in)    :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow  !< The time interval of these diagnostics
  type(time_type),            intent(in)    :: Time     !< The ending time of these diagnostics
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module
  type(SIS_diag_ctrl),        pointer       :: diag !< A structure that is used to regulate diagnostic output

  ! Local variables
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: mass, mass_ice, mass_snow, tmp2d
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,IG%NkIce) :: &
    temp_ice    ! A diagnostic array with the ice temperature [degC].
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    temp_snow   ! A diagnostic array with the snow temperature [degC].
  ! ### This diagnostic does not exist yet.
  ! real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
  !   rdg_frac    ! fraction of ridged ice per category
  real, dimension(SZI_(G),SZJ_(G))   :: diagVar ! An temporary array for diagnostics.
  real, dimension(IG%NkIce) :: S_col ! Specified thermodynamic salinity of each
                                     ! ice layer if spec_thermo_sal is true.
  real :: rho_ice  ! The nominal density of sea ice [kg m-3].
  real :: rho_snow ! The nominal density of snow [kg m-3].
  real :: enth_units, I_enth_units
  real :: tmp_mca  ! A temporary cell averaged mass [H ~> kg m-2].
  real :: I_Nk        ! The inverse of the number of layers in the ice.
  real :: Idt_slow ! The inverse of the thermodynamic step [s-1].
  logical :: spec_thermo_sal
  logical :: do_temp_diags
  integer :: i, j, k, l, m, isc, iec, jsc, jec, ncat, NkIce

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  NkIce = IG%NkIce
  I_Nk = 1.0 / NkIce
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  ! Sum the concentration weighted mass for diagnostics.
  if (CS%id_mi>0 .or. CS%id_mib>0) then
    mass_ice(:,:) = 0.0
    mass_snow(:,:) = 0.0
    mass(:,:) = 0.0
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      mass_ice(i,j) = mass_ice(i,j) + IG%H_to_kg_m2*IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
      mass_snow(i,j) = mass_snow(i,j) + IG%H_to_kg_m2*IST%mH_snow(i,j,k)*IST%part_size(i,j,k)
      mass(i,j) = mass_ice(i,j) + mass_snow(i,j)
    enddo ; enddo ; enddo

    if (CS%id_simass>0) call post_data(CS%id_simass, mass_ice(isc:iec,jsc:jec), diag)
    if (CS%id_sisnmass>0) call post_data(CS%id_sisnmass, mass_snow(isc:iec,jsc:jec), diag)
    if (CS%id_mi>0) call post_data(CS%id_mi, mass(isc:iec,jsc:jec), diag)

    if (CS%id_mib>0) then
      if (associated(IOF%mass_berg)) then
        do j=jsc,jec ; do i=isc,iec
          mass(i,j) = (mass(i,j) + IOF%mass_berg(i,j)) ! Add icebergs mass [kg m-2]
        enddo ; enddo
      endif
      call post_data(CS%id_mib, mass(isc:iec,jsc:jec), diag)
    endif
  endif

  !
  ! Thermodynamic state diagnostics
  !
  if (CS%id_cn>0) call post_data(CS%id_cn, IST%part_size(:,:,1:ncat), diag)
  if (CS%id_siconc>0) call post_data(CS%id_siconc, sum(IST%part_size(:,:,1:ncat),3), diag)

  ! TK Mod: 10/18/02
  !  if (CS%id_obs_cn>0) call post_data(CS%id_obs_cn, Obs_cn_ice(:,:,2), diag)
  ! TK Mod: 10/18/02: (commented out...does not compile yet... add later)
  !  if (CS%id_obs_hi>0) &
  !    call post_avg(CS%id_obs_hi, Obs_h_ice(isc:iec,jsc:jec), IST%part_size(isc:iec,jsc:jec,1:), &
  !                  diag, G=G, wtd=.true.)

  !   Convert from ice and snow enthalpy back to temperature for diagnostic purposes.
  do_temp_diags = (CS%id_tsn > 0)
  do m=1,NkIce ; if (CS%id_t(m)>0) do_temp_diags = .true. ; enddo
  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             rho_ice=rho_ice, rho_snow=rho_snow, &
                             specified_thermo_salinity=spec_thermo_sal)
  I_enth_units = 1.0 / enth_units

  if (do_temp_diags) then
    !$OMP parallel do default(shared)
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

  if (CS%id_ext>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (IST%part_size(i,j,0) < 0.85) diagVar(i,j) = 1.0
    enddo ; enddo
    call post_data(CS%id_ext, diagVar, diag)
  endif
  if (CS%id_hp>0) call post_avg(CS%id_hp, IST%mH_pond, IST%part_size(:,:,1:), & ! mw/new
                                 diag, G=G, &
                                 scale=IG%H_to_kg_m2/1e3, wtd=.true.) ! rho_water=1e3
  if (CS%id_hs>0) call post_avg(CS%id_hs, IST%mH_snow, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_snow, wtd=.true.)
  if (CS%id_sisnthick>0) call post_avg(CS%id_sisnthick, IST%mH_snow, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_snow, wtd=.true.)
  if (CS%id_hi>0) call post_avg(CS%id_hi, IST%mH_ice, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_ice, wtd=.true.)
  if (CS%id_sithick>0) call post_avg(CS%id_sithick, IST%mH_ice, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_ice, wtd=.true.)
  if (CS%id_sivol>0) call post_avg(CS%id_sivol, IST%mH_ice, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_ice, wtd=.true.)
  if (CS%id_tsn>0) call post_avg(CS%id_tsn, temp_snow, IST%part_size(:,:,1:), &
                                 diag, G=G, wtd=.true.)
  if (CS%id_sitimefrac>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (IST%part_size(i,j,0) < 1.0) diagVar(i,j) = 1.0
    enddo ; enddo
    call post_data(CS%id_sitimefrac, diagVar, diag)
  endif
  if (CS%id_sisnconc>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec; do k=1,ncat
      if (IST%part_size(i,j,k) > 0.0 .and. IST%mH_snow(i,j,k) > 0.0) then
        diagVar(i,j) = diagVar(i,j) + IST%part_size(i,j,k)
      endif
    enddo ; enddo ; enddo
    call post_data(CS%id_sisnconc, diagVar, diag)
  endif

  do m=1,NkIce
    if (CS%id_t(m)>0) call post_avg(CS%id_t(m), temp_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   diag, G=G, wtd=.true.)
    if (CS%id_sal(m)>0) call post_avg(CS%id_sal(m), IST%sal_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   diag, G=G, wtd=.true.)
  enddo
  if (CS%id_t_iceav>0) call post_avg(CS%id_t_iceav, temp_ice, IST%part_size(:,:,1:), &
                                    diag, G=G, wtd=.true.)
  if (CS%id_S_iceav>0) call post_avg(CS%id_S_iceav, IST%sal_ice, IST%part_size(:,:,1:), &
                                    diag, G=G, wtd=.true.)

  ! Write out diagnostics of the ocean surface state, as seen by the slow sea ice.
  ! These fields do not change over the course of the sea-ice time stepping.
  call post_ocean_sfc_diagnostics(OSS, dt_slow, Time, G, diag)

  if (CS%id_e2m>0) then
    tmp2d(:,:) = 0.0
    !$OMP parallel do default(shared)
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
    call post_data(CS%id_e2m,  tmp2d(:,:), diag)
  endif

  if (CS%do_ridging) then
  !TOM> preparing output field fraction of ridged ice rdg_frac = (ridged ice volume) / (total ice volume)
  !     in each category; IST%rdg_mice is ridged ice mass per unit total area throughout the code.
!     if (CS%id_rdgf>0) then
!       !$OMP parallel do default(shared) private(tmp_mca)
!       do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
!         tmp_mca = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
!         if (tmp_mca > Rho_Ice*1.e-5*IG%kg_m2_to_H) then  ! 1 mm ice thickness x 1% ice concentration
!           rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp_mca
!         else
!           rdg_frac(i,j,k) = 0.0
!         endif
!       enddo ; enddo ; enddo
!       call post_data(CS%id_rdgf, rdg_frac(isc:iec,jsc:jec), diag)
!     endif
  endif

end subroutine post_ice_state_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Offer diagnostics of the ocean surface field, as seen by the sea ice.
subroutine post_ocean_sfc_diagnostics(OSS, dt_slow, Time, G, diag)
  type(ocean_sfc_state_type), intent(in)    :: OSS  !< A structure containing the arrays that describe
                                                    !! the ocean's surface state for the ice model.
  real,                       intent(in)    :: dt_slow  !< The time interval of these diagnostics
  type(time_type),            intent(in)    :: Time     !< The ending time of these diagnostics
  type(SIS_hor_grid_type),    intent(inout) :: G    !< The horizontal grid type
  type(SIS_diag_ctrl),        pointer       :: diag !< A structure that is used to regulate diagnostic output

  real :: Idt_slow ! The inverse of the thermodynamic step [s-1].
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  ! Write out diagnostics of the ocean surface state, as seen by the slow sea ice.
  ! These fields do not change over the course of the sea-ice time stepping.
  if (OSS%id_sst>0) call post_data(OSS%id_sst, OSS%SST_C, diag)
  if (OSS%id_sss>0) call post_data(OSS%id_sss, OSS%s_surf, diag)
  if (OSS%id_ssh>0) call post_data(OSS%id_ssh, OSS%sea_lev, diag)
  if (allocated(OSS%u_ocn_C)) then
    if (OSS%id_uo>0) call post_data(OSS%id_uo, OSS%u_ocn_C, diag)
    if (OSS%id_vo>0) call post_data(OSS%id_vo, OSS%v_ocn_C, diag)
  else
    if (OSS%id_uo>0) call post_data(OSS%id_uo, OSS%u_ocn_B, diag)
    if (OSS%id_vo>0) call post_data(OSS%id_vo, OSS%v_ocn_B, diag)
  endif
  if (OSS%id_frazil>0) &
    call post_data(OSS%id_frazil, OSS%frazil*Idt_slow, diag)

  if (coupler_type_initialized(OSS%tr_fields)) &
    call coupler_type_send_data(OSS%tr_fields, Time)

end subroutine post_ocean_sfc_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Finish setting the ice-ocean stresses by dividing the running sums of the
!! stresses by the number of times they have been augmented.
subroutine finish_ocean_top_stresses(IOF, G)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(in)    :: G   !< The horizontal grid type

  real :: taux2, tauy2  ! squared wind stresses (Pa^2)
  real :: I_count ! The number of times IOF has been incremented.

  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (IOF%stress_count > 1) then
    I_count = 1.0 / IOF%stress_count
    if (IOF%flux_uv_stagger == AGRID) then
      do j=jsc,jec ; do i=isc,iec
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) * I_count
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) * I_count
      enddo ; enddo
    elseif (IOF%flux_uv_stagger == BGRID_NE) then
      do J=jsc-1,jec ; do I=isc-1,iec
        IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) * I_count
        IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) * I_count
      enddo ; enddo
    elseif (IOF%flux_uv_stagger == CGRID_NE) then
      do j=jsc,jec ; do I=isc-1,iec
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) * I_count
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) * I_count
      enddo ; enddo
    else
      call SIS_error(FATAL, "finish_ocean_top_stresses: Unrecognized flux_uv_stagger.")
    endif
  endif

  if (allocated(IOF%stress_mag)) then
    ! if (IOF%simple_mag) then
    ! Determine the magnitude of the time and area mean stresses.
    call stresses_to_stress_mag(G, IOF%flux_u_ocn, IOF%flux_v_ocn, IOF%flux_uv_stagger, IOF%stress_mag)
    ! elseif (IOF%stress_count > 1) then
    !   ! Find the time mean magnitude of the stresses on the ocean.
    !   do j=jsc,jec ; do i=isc,iec
    !     IOF%stress_mag(i,j) = IOF%stress_mag(i,j) * I_count
    !   enddo ; enddo
    ! endif
  endif

end subroutine finish_ocean_top_stresses

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Translate the ice-ocean stresses with various staggering options into the
!! magnitude of the stresses averaged to tracer points.  The stresses must already
!! be set at all (symmetric) edge points unless stagger is AGRID.
subroutine stresses_to_stress_mag(G, str_x, str_y, stagger, stress_mag)
  type(SIS_hor_grid_type),   intent(in)    :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: str_x !< The x-direction ice to ocean stress [Pa].
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: str_y !< The y-direction ice to ocean stress [Pa].
  integer,                   intent(in)    :: stagger !< The staggering relative to the tracer points of the
                                                  !! two wind stress components. Valid entries include AGRID,
                                                  !! BGRID_NE, and CGRID_NE, following the Arakawa
                                                  !! grid-staggering  notation.  BGRID_SW and CGRID_SW are
                                                  !! possibilties that have not been implemented yet.
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(inout) :: stress_mag !< The magnitude of the stress at tracer points [Pa].

  ! Local variables
  real :: taux2, tauy2  ! squared wind stress components (Pa^2)
  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (stagger == AGRID) then
    do j=jsc,jec ; do i=isc,iec
      stress_mag(i,j) = sqrt(str_x(i,j)**2 + str_y(i,j)**2)
    enddo ; enddo
  elseif (stagger == BGRID_NE) then
    do j=jsc,jec ; do i=isc,iec
      if (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
           (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0) then
        stress_mag(i,j) = sqrt( &
          ((G%mask2dBu(I,J)*(str_x(I,J)**2 + str_y(I,J)**2) + &
            G%mask2dBu(I-1,J-1)*(str_x(I-1,J-1)**2 + str_y(I-1,J-1)**2)) + &
           (G%mask2dBu(I,J-1)*(str_x(I,J-1)**2 + str_y(I,J-1)**2) + &
            G%mask2dBu(I-1,J)*(str_x(I-1,J)**2 + str_y(I-1,J)**2)) ) / &
          ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) )
      else
        stress_mag(i,j) = 0.0
      endif
    enddo ; enddo
  elseif (stagger == CGRID_NE) then
    do j=jsc,jec ; do i=isc,iec
      taux2 = 0.0 ; tauy2 = 0.0
      if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0) &
        taux2 = (G%mask2dCu(I-1,j)*str_x(I-1,j)**2 + &
                 G%mask2dCu(I,j)*str_x(I,j)**2) / (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0) &
        tauy2 = (G%mask2dCv(i,J-1)*str_y(i,J-1)**2 + &
                 G%mask2dCv(i,J)*str_y(i,J)**2) / (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))
      stress_mag(i,j) = sqrt(taux2 + tauy2)
    enddo ; enddo
  else
    call SIS_error(FATAL, "stresses_to_stress_mag: Unrecognized stagger.")
  endif

end subroutine stresses_to_stress_mag

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Calculate the stresses on the ocean integrated across all the thickness categories with
!! the appropriate staggering, and store them in the public ice data type for use by the
!! ocean model.  This version of the routine uses wind and ice-ocean stresses on a B-grid.
subroutine set_ocean_top_stress_Bgrid(IOF, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, part_size, G, IG)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),       intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over open water [Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over open water [Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_x   !< The x-direction ice to ocean stress [Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y   !< The y-direction ice to ocean stress [Pa].
  real, dimension(SZI_(G),SZJ_(G),0:IG%CatIce), &
                             intent(in)    :: part_size !< The fractional area coverage of the ice
                                                  !! thickness categories [nondim], 0-1

  real    :: ps_vel ! part_size interpolated to a velocity point [nondim].
  integer :: i, j, k, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce


  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Bgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.
  if (IOF%flux_uv_stagger == AGRID) then
    !$OMP parallel do default(shared) private(ps_vel)
    do j=jsc,jec
      do i=isc,iec
        ps_vel = G%mask2dT(i,j) * part_size(i,j,0)
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * 0.25 * &
            ((windstr_x_water(I,J) + windstr_x_water(I-1,J-1)) + &
             (windstr_x_water(I-1,J) + windstr_x_water(I,J-1)))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * 0.25 * &
            ((windstr_y_water(I,J) + windstr_y_water(I-1,J-1)) + &
             (windstr_y_water(I-1,J) + windstr_y_water(I,J-1)))
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + part_size(i,j,k) * 0.25 * &
            ((str_ice_oce_x(I,J) + str_ice_oce_x(I-1,J-1)) + &
             (str_ice_oce_x(I-1,J) + str_ice_oce_x(I,J-1)))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + part_size(i,j,k) * 0.25 * &
            ((str_ice_oce_y(I,J) + str_ice_oce_y(I-1,J-1)) + &
             (str_ice_oce_y(I-1,J) + str_ice_oce_y(I,J-1)))
      endif ; enddo ; enddo
    enddo
  elseif (IOF%flux_uv_stagger == BGRID_NE) then
    !$OMP parallel do default(shared) private(ps_vel)
    do J=jsc-1,jec
      do I=isc-1,iec
        ps_vel = 1.0 ; if (G%mask2dBu(I,J)>0.5) ps_vel = &
                           0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                                 (part_size(i+1,j,0) + part_size(i,j+1,0)) )
        IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + windstr_x_water(I,J) * ps_vel
        IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + windstr_y_water(I,J) * ps_vel
      enddo
      do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dBu(I,J)>0.5) then
        ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                         (part_size(i+1,j,k) + part_size(i,j+1,k)) )
        IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + str_ice_oce_x(I,J) * ps_vel
        IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + str_ice_oce_y(I,J) * ps_vel
      endif ; enddo ; enddo
    enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
    !$OMP parallel do default(shared) private(ps_vel)
    do j=jsc,jec
      do I=isc-1,iec
        ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.5) ps_vel = &
                           0.5*(part_size(i+1,j,0) + part_size(i,j,0))
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * &
                0.5 * (windstr_x_water(I,J) + windstr_x_water(I,J-1))
      enddo
      do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dCu(I,j)>0.5) then
        ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * &
            0.5 * (str_ice_oce_x(I,J) + str_ice_oce_x(I,J-1))
      endif ; enddo ; enddo
    enddo
    !$OMP parallel do default(shared) private(ps_vel)
    do J=jsc-1,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.5) ps_vel = &
                           0.5*(part_size(i,j+1,0) + part_size(i,j,0))
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * &
                0.5 * (windstr_y_water(I,J) + windstr_y_water(I-1,J))
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dCv(i,J)>0.5) then
        ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * &
                0.5 * (str_ice_oce_y(I,J) + str_ice_oce_y(I-1,J))
      endif ; enddo ; enddo
    enddo
  else
    call SIS_error(FATAL, "set_ocean_top_stress_Bgrid: Unrecognized flux_uv_stagger.")
  endif
  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_Bgrid

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Calculate the stresses on the ocean integrated across all the thickness categories with the
!! appropriate staggering, and store them in the public ice data type for use by the ocean
!! model.  This version of the routine uses wind and ice-ocean stresses on a C-grid.
subroutine set_ocean_top_stress_Cgrid(IOF, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, part_size, G, IG)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),       intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over open water [Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over open water [Pa].
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: str_ice_oce_x   !< The x-direction ice to ocean stress [Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y   !< The y-direction ice to ocean stress [Pa].
  real, dimension(SZI_(G),SZJ_(G),0:IG%CatIce), &
                             intent(in)    :: part_size !< The fractional area coverage of the ice
                                                  !! thickness categories [nondim], 0-1

  real    :: ps_vel ! part_size interpolated to a velocity point [nondim].
  integer :: i, j, k, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Cgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.

  if (IOF%flux_uv_stagger == AGRID) then
    !$OMP parallel do default(shared) private(ps_vel)
    do j=jsc,jec
      do i=isc,iec
        ps_vel = G%mask2dT(i,j) * part_size(i,j,0)
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * 0.5 * &
                            (windstr_x_water(I,j) + windstr_x_water(I-1,j))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * 0.5 * &
                            (windstr_y_water(i,J) + windstr_y_water(i,J-1))
      enddo
      !### SIMPLIFY THIS TO USE THAT sum(part_size(i,j,1:ncat)) = 1.0-part_size(i,j,0) ?
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) +  part_size(i,j,k) * 0.5 * &
                            (str_ice_oce_x(I,j) + str_ice_oce_x(I-1,j))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + part_size(i,j,k) * 0.5 * &
                            (str_ice_oce_y(i,J) + str_ice_oce_y(i,J-1))
      endif ; enddo ; enddo
    enddo
  elseif (IOF%flux_uv_stagger == BGRID_NE) then
    !$OMP parallel do default(shared) private(ps_vel)
    do J=jsc-1,jec
      do I=isc-1,iec
        ps_vel = 1.0 ; if (G%mask2dBu(I,J)>0.5) ps_vel = &
                           0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                                 (part_size(i+1,j,0) + part_size(i,j+1,0)) )
        !### Consider deleting the masks here?  They probably do not change answers.
        IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
                (windstr_x_water(I,j) + windstr_x_water(I,j+1))
        IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
                (windstr_y_water(i,J) + windstr_y_water(i+1,J))
      enddo
      do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dBu(I,J)>0.5) then
        ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                         (part_size(i+1,j,k) + part_size(i,j+1,k)) )
        IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + ps_vel * 0.5 * &
                            (str_ice_oce_x(I,j) + str_ice_oce_x(I,j+1))
        IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + ps_vel * 0.5 * &
                            (str_ice_oce_y(i,J) + str_ice_oce_y(i+1,J))
      endif ; enddo ; enddo
    enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
    !$OMP parallel do default(shared) private(ps_vel)
    do j=jsc,jec
      do I=Isc-1,iec
        ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.5) ps_vel = &
                           0.5*(part_size(i+1,j,0) + part_size(i,j,0))
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * windstr_x_water(I,j)
      enddo
      do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dCu(I,j)>0.5) then
        ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * str_ice_oce_x(I,j)
      endif ; enddo ; enddo
    enddo
    !$OMP parallel do default(shared) private(ps_vel)
    do J=jsc-1,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.5) ps_vel = &
                           0.5*(part_size(i,j+1,0) + part_size(i,j,0))
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * windstr_y_water(i,J)
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dCv(i,J)>0.5) then
        ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * str_ice_oce_y(i,J)
      endif ; enddo ; enddo
    enddo
  else
    call SIS_error(FATAL, "set_ocean_top_stress_Cgrid: Unrecognized flux_uv_stagger.")
  endif

  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_Cgrid

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Calculate the stresses on the ocean integrated across all the thickness categories with
!! the appropriate staggering, and store them in the public ice data type for use by the
!! ocean model.  This version of the routine uses wind and ice-ocean stresses on a B-grid.
subroutine set_ocean_top_stress_B2(IOF, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, ice_free, ice_cover, G)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over open water [Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over open water [Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_x   !< The x-direction ice to ocean stress [Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y   !< The y-direction ice to ocean stress [Pa].
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: ice_free  !< The fractional open water area coverage [nondim], 0-1
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: ice_cover !< The fractional ice area coverage [nondim], 0-1

  real    :: ps_ice, ps_ocn ! ice_free and ice_cover interpolated to a velocity point [nondim].
  integer :: i, j, k, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Bgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.
  if (IOF%flux_uv_stagger == AGRID) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do i=isc,iec
      ps_ocn = G%mask2dT(i,j) * ice_free(i,j)
      ps_ice = G%mask2dT(i,j) * ice_cover(i,j)
      IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + 0.25 * &
            (ps_ocn * ((windstr_x_water(I,J) + windstr_x_water(I-1,J-1)) + &
                       (windstr_x_water(I-1,J) + windstr_x_water(I,J-1))) + &
             ps_ice * ((str_ice_oce_x(I,J) + str_ice_oce_x(I-1,J-1)) + &
                       (str_ice_oce_x(I-1,J) + str_ice_oce_x(I,J-1))) )
      IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + 0.25 * &
            (ps_ocn * ((windstr_y_water(I,J) + windstr_y_water(I-1,J-1)) + &
                       (windstr_y_water(I-1,J) + windstr_y_water(I,J-1))) + &
             ps_ice * ((str_ice_oce_y(I,J) + str_ice_oce_y(I-1,J-1)) + &
                       (str_ice_oce_y(I-1,J) + str_ice_oce_y(I,J-1))) )
    enddo ; enddo
  elseif (IOF%flux_uv_stagger == BGRID_NE) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do J=jsc-1,jec ; do I=isc-1,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dBu(I,J)>0.5) then
        ps_ocn = 0.25 * ((ice_free(i+1,j+1) + ice_free(i,j)) + &
                         (ice_free(i+1,j) + ice_free(i,j+1)) )
        ps_ice = 0.25 * ((ice_cover(i+1,j+1) + ice_cover(i,j)) + &
                         (ice_cover(i+1,j) + ice_cover(i,j+1)) )
      endif
      IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + (ps_ocn * windstr_x_water(I,J) + ps_ice * str_ice_oce_x(I,J))
      IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + (ps_ocn * windstr_y_water(I,J) + ps_ice * str_ice_oce_y(I,J))
    enddo ; enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do I=isc-1,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dCu(I,j)>0.5) then
        ps_ocn = 0.5*(ice_free(i+1,j) + ice_free(i,j))
        ps_ice = 0.5*(ice_cover(i+1,j) + ice_cover(i,j))
      endif
      IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + 0.5 * &
            (ps_ocn * (windstr_x_water(I,J) + windstr_x_water(I,J-1)) + &
             ps_ice * (str_ice_oce_x(I,J) + str_ice_oce_x(I,J-1)) )
    enddo ; enddo
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do J=jsc-1,jec ; do i=isc,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dCv(i,J)>0.5) then
        ps_ocn = 0.5*(ice_free(i,j+1) + ice_free(i,j))
        ps_ice = 0.5*(ice_cover(i,j+1) + ice_cover(i,j))
      endif
      IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + 0.5 * &
            (ps_ocn * (windstr_y_water(I,J) + windstr_y_water(I-1,J)) + &
             ps_ice * (str_ice_oce_y(I,J) + str_ice_oce_y(I-1,J)) )
    enddo ; enddo
  else
    call SIS_error(FATAL, "set_ocean_top_stress_B2: Unrecognized flux_uv_stagger.")
  endif
  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_B2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Calculate the stresses on the ocean integrated across all the thickness categories with the
!! appropriate staggering, and store them in the public ice data type for use by the ocean
!! model.  This version of the routine uses wind and ice-ocean stresses on a C-grid.
subroutine set_ocean_top_stress_C2(IOF, windstr_x_water, windstr_y_water, &
                                   str_ice_oce_x, str_ice_oce_y, ice_free, ice_cover, G)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over open water [Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over open water [Pa].
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: str_ice_oce_x   !< The x-direction ice to ocean stress [Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y   !< The y-direction ice to ocean stress [Pa].
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: ice_free  !< The fractional open water area coverage [nondim], 0-1
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: ice_cover !< The fractional ice area coverage [nondim], 0-1

  real    :: ps_ice, ps_ocn ! ice_free and ice_cover interpolated to a velocity point [nondim].
  integer :: i, j, k, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Cgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.

  if (IOF%flux_uv_stagger == AGRID) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do i=isc,iec
      ps_ocn = G%mask2dT(i,j) * ice_free(i,j)
      ps_ice = G%mask2dT(i,j) * ice_cover(i,j)
      IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + &
           (ps_ocn * 0.5 * (windstr_x_water(I,j) + windstr_x_water(I-1,j)) + &
            ps_ice * 0.5 * (str_ice_oce_x(I,j) + str_ice_oce_x(I-1,j)) )
      IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + &
           (ps_ocn * 0.5 * (windstr_y_water(i,J) + windstr_y_water(i,J-1)) + &
            ps_ice * 0.5 * (str_ice_oce_y(i,J) + str_ice_oce_y(i,J-1)) )
    enddo ; enddo
  elseif (IOF%flux_uv_stagger == BGRID_NE) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do J=jsc-1,jec ; do I=isc-1,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dBu(I,J)>0.5) then
        ps_ocn = 0.25 * ((ice_free(i+1,j+1) + ice_free(i,j)) + &
                         (ice_free(i+1,j) + ice_free(i,j+1)) )
        ps_ice = 0.25 * ((ice_cover(i+1,j+1) + ice_cover(i,j)) + &
                         (ice_cover(i+1,j) + ice_cover(i,j+1)) )
      endif
      IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + &
          (ps_ocn * 0.5 * (windstr_x_water(I,j) + windstr_x_water(I,j+1)) + &
           ps_ice * 0.5 * (str_ice_oce_x(I,j) + str_ice_oce_x(I,j+1)) )
      IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + &
          (ps_ocn * 0.5 * (windstr_y_water(i,J) + windstr_y_water(i+1,J)) + &
           ps_ice * 0.5 * (str_ice_oce_y(i,J) + str_ice_oce_y(i+1,J)) )
    enddo ; enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do I=Isc-1,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dCu(I,j)>0.5) then
        ps_ocn = 0.5*(ice_free(i+1,j) + ice_free(i,j))
        ps_ice = 0.5*(ice_cover(i+1,j) + ice_cover(i,j))
      endif
      IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + &
          (ps_ocn * windstr_x_water(I,j) + ps_ice * str_ice_oce_x(I,j))
    enddo ; enddo
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do J=jsc-1,jec ; do i=isc,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dCv(i,J)>0.5) then
        ps_ocn = 0.5*(ice_free(i,j+1) + ice_free(i,j))
        ps_ice = 0.5*(ice_cover(i,j+1) + ice_cover(i,j))
      endif
      IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + &
          (ps_ocn * windstr_y_water(i,J) + ps_ice * str_ice_oce_y(i,J))
    enddo ; enddo
  else
    call SIS_error(FATAL, "set_ocean_top_stress_C2: Unrecognized flux_uv_stagger.")
  endif

  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_C2



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Calculate the stresses on the ocean integrated across all the thickness categories
!! with the appropriate staggering, based on the information in a fast_ice_avg_type.
subroutine set_ocean_top_stress_FIA(FIA, IOF, G)
  type(fast_ice_avg_type),   intent(inout) :: FIA !< A type containing averages of fields
                                                  !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type

  real    :: ps_ice, ps_ocn ! ice_free and ice_cover interpolated to a velocity point [nondim].
  integer :: i, j, k, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Cgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.

  if (IOF%flux_uv_stagger == AGRID) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do i=isc,iec
      ps_ocn = G%mask2dT(i,j) * FIA%ice_free(i,j)
      ps_ice = G%mask2dT(i,j) * FIA%ice_cover(i,j)
      IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + &
           (ps_ocn * FIA%WindStr_ocn_x(i,j) + ps_ice * FIA%WindStr_x(i,j))
      IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + &
           (ps_ocn * FIA%WindStr_ocn_y(i,j) + ps_ice * FIA%WindStr_y(i,j))
    enddo ; enddo
  elseif (IOF%flux_uv_stagger == BGRID_NE) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do J=jsc-1,jec ; do I=isc-1,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dBu(I,J)>0.5) then
        ps_ocn = 0.25 * ((FIA%ice_free(i+1,j+1) + FIA%ice_free(i,j)) + &
                         (FIA%ice_free(i+1,j) + FIA%ice_free(i,j+1)) )
        ps_ice = 0.25 * ((FIA%ice_cover(i+1,j+1) + FIA%ice_cover(i,j)) + &
                         (FIA%ice_cover(i+1,j) + FIA%ice_cover(i,j+1)) )
      endif
      IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + &
          (ps_ocn * 0.25 * ((FIA%WindStr_ocn_x(i,j) + FIA%WindStr_ocn_x(i+1,j+1)) + &
                            (FIA%WindStr_ocn_x(i,j+1) + FIA%WindStr_ocn_x(i+1,j))) + &
           ps_ice * 0.25 * ((FIA%WindStr_x(i,j) + FIA%WindStr_x(i+1,j+1)) + &
                            (FIA%WindStr_x(i,j+1) + FIA%WindStr_x(i+1,J))) )
      IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + &
          (ps_ocn * 0.25 * ((FIA%WindStr_ocn_y(i,j) + FIA%WindStr_ocn_y(i+1,j+1)) + &
                            (FIA%WindStr_ocn_y(i,j+1) + FIA%WindStr_ocn_y(i+1,j))) + &
           ps_ice * 0.25 * ((FIA%WindStr_y(i,j) + FIA%WindStr_y(i+1,j+1)) + &
                            (FIA%WindStr_y(i,j+1) + FIA%WindStr_y(i+1,J))) )
    enddo ; enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do I=Isc-1,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dCu(I,j)>0.5) then
        ps_ocn = 0.5*(FIA%ice_free(i+1,j) + FIA%ice_free(i,j))
        ps_ice = 0.5*(FIA%ice_cover(i+1,j) + FIA%ice_cover(i,j))
      endif
      IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + &
           (ps_ocn * 0.5 * (FIA%WindStr_ocn_x(i+1,j) + FIA%WindStr_ocn_x(i,j)) + &
            ps_ice * 0.5 * (FIA%WindStr_x(i+1,j) + FIA%WindStr_x(i,j)) )
    enddo ; enddo
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do J=jsc-1,jec ; do i=isc,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dCv(i,J)>0.5) then
        ps_ocn = 0.5*(FIA%ice_free(i,j+1) + FIA%ice_free(i,j))
        ps_ice = 0.5*(FIA%ice_cover(i,j+1) + FIA%ice_cover(i,j))
      endif
      IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + &
          (ps_ocn * 0.5 * (FIA%WindStr_ocn_y(i,j+1) + FIA%WindStr_ocn_y(i,j)) + &
           ps_ice * 0.5 * (FIA%WindStr_y(i,j+1) + FIA%WindStr_y(i,j)) )
    enddo ; enddo
  else
    call SIS_error(FATAL, "set_ocean_top_stress_C2: Unrecognized flux_uv_stagger.")
  endif

  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_FIA




!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_wind_stresses_C determines the wind stresses on the ice and open ocean with
!!   a C-grid staggering of the points.
subroutine set_wind_stresses_C(FIA, ice_cover, ice_free, WindStr_x_Cu, WindStr_y_Cv, &
                               WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, G, ncat)
  type(fast_ice_avg_type),          intent(in)    :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),           intent(in)   :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)   :: &
    ice_cover, &        !< The fractional ice coverage, summed across all
                        !! thickness categories [nondim], between 0 & 1.
    ice_free            !< The fractional open water [nondim], between 0 & 1.
  real, dimension(SZIB_(G),SZJ_(G)), intent(out)  :: &
    WindStr_x_Cu, &   !< Zonal wind stress averaged over the ice categores on C-grid u-points [Pa].
    WindStr_x_ocn_Cu  !< Zonal wind stress on the ice-free ocean on C-grid u-points [Pa].
  real, dimension(SZI_(G),SZJB_(G)), intent(out)  :: &
    WindStr_y_Cv, &   !< Meridional wind stress averaged over the ice categores on C-grid v-points [Pa].
    WindStr_y_ocn_Cv  !< Meridional wind stress on the ice-free ocean on C-grid v-points [Pa].
  integer,                           intent(in)   :: ncat !< The number of ice thickness categories.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    WindStr_x_A, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_A, &      ! averaged over the ice categories on an A-grid [Pa].
    WindStr_x_ocn_A, &  ! Zonal (_x_) and meridional (_y_) wind stresses on the
    WindStr_y_ocn_A     ! ice-free ocean on an A-grid [Pa].
  real :: weights  ! A sum of the weights around a point.
  real :: I_wts    ! 1.0 / wts or 0 if wts is 0 [nondim].
  real :: max_ice_cover, FIA_ice_cover, ice_cover_now
  integer :: i, j, isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  max_ice_cover = 1.0 - 2.0*ncat*epsilon(max_ice_cover)
  !$OMP parallel do default(shared) private(FIA_ice_cover, ice_cover_now)
  do j=jsd,jed ; do i=isd,ied
    ! The use of these limits prevents the use of the ocean wind stresses if
    ! there is actually no open ocean and hence there may be no valid ocean
    ! stresses.  This can occur when ice_cover ~= 1 for both states, but
    ! they are not exactly 1.0 due to roundoff in the sum across categories above.
    ice_cover_now = min(ice_cover(i,j), max_ice_cover)
    FIA_ice_cover = min(FIA%ice_cover(i,j), max_ice_cover)

    if (ice_cover_now > FIA_ice_cover) then
      WindStr_x_A(i,j) = ((ice_cover_now-FIA_ice_cover)*FIA%WindStr_ocn_x(i,j) + &
                          FIA_ice_cover*FIA%WindStr_x(i,j)) / ice_cover_now
      WindStr_y_A(i,j) = ((ice_cover_now-FIA_ice_cover)*FIA%WindStr_ocn_y(i,j) + &
                          FIA_ice_cover*FIA%WindStr_y(i,j)) / ice_cover_now
    else
      WindStr_x_A(i,j) = FIA%WindStr_x(i,j)
      WindStr_y_A(i,j) = FIA%WindStr_y(i,j)
    endif

    if (ice_free(i,j) <= FIA%ice_free(i,j)) then
      WindStr_x_ocn_A(i,j) = FIA%WindStr_ocn_x(i,j)
      WindStr_y_ocn_A(i,j) = FIA%WindStr_ocn_y(i,j)
    else
      WindStr_x_ocn_A(i,j) = ((ice_free(i,j)-FIA%ice_free(i,j))*FIA%WindStr_x(i,j) + &
                          FIA%ice_free(i,j)*FIA%WindStr_ocn_x(i,j)) / ice_free(i,j)
      WindStr_y_ocn_A(i,j) = ((ice_free(i,j)-FIA%ice_free(i,j))*FIA%WindStr_y(i,j) + &
                          FIA%ice_free(i,j)*FIA%WindStr_ocn_y(i,j)) / ice_free(i,j)
    endif
  enddo ;  enddo

  !   The j-loop extents here are larger than they would normally be in case
  ! the stresses are being passed to the ocean on a B-grid.
  !$OMP parallel default(shared) private(weights,I_wts)
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

end subroutine set_wind_stresses_C


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_wind_stresses_B determines the wind stresses on the ice and open ocean with
!!   a B-grid staggering of the points.
subroutine set_wind_stresses_B(FIA, ice_cover, ice_free, WindStr_x_B, WindStr_y_B, &
                               WindStr_x_ocn_B, WindStr_y_ocn_B, G, ncat)
  type(fast_ice_avg_type),            intent(in)   :: FIA !< A type containing averages of fields
                                                          !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),            intent(in)   :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)   :: &
    ice_cover, &        !< The fractional ice coverage, summed across all
                        !! thickness categories [nondim], between 0 & 1.
    ice_free            !< The fractional open water [nondim], between 0 & 1.
  real, dimension(SZIB_(G),SZJB_(G)), intent(out) :: &
    WindStr_x_B, &      !< Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_B, &      !< averaged over the ice categories on a B-grid [Pa].
    WindStr_x_ocn_B, &  !< Zonal wind stress on the ice-free ocean on a B-grid [Pa].
    WindStr_y_ocn_B     !< Meridional wind stress on the ice-free ocean on a B-grid [Pa].
  integer,                            intent(in)   :: ncat !< The number of ice thickness categories.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    WindStr_x_A, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_A, &      ! averaged over the ice categories on an A-grid [Pa].
    WindStr_x_ocn_A, &  ! Zonal (_x_) and meridional (_y_) wind stresses on the
    WindStr_y_ocn_A     ! ice-free ocean on an A-grid [Pa].
  real :: weights  ! A sum of the weights around a point.
  real :: I_wts    ! 1.0 / wts or 0 if wts is 0 [nondim].
  real :: max_ice_cover, FIA_ice_cover, ice_cover_now
  integer :: i, j, isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  max_ice_cover = 1.0 - 2.0*ncat*epsilon(max_ice_cover)
  !$OMP parallel do default(shared) private(FIA_ice_cover, ice_cover_now)
  do j=jsd,jed ; do i=isd,ied
    ! The use of these limits prevents the use of the ocean wind stresses if
    ! there is actually no open ocean and hence there may be no valid ocean
    ! stresses.  This can occur when ice_cover ~= 1 for both states, but
    ! they are not exactly 1.0 due to roundoff in the sum across categories above.
    ice_cover_now = min(ice_cover(i,j), max_ice_cover)
    FIA_ice_cover = min(FIA%ice_cover(i,j), max_ice_cover)

    if (ice_cover_now > FIA_ice_cover) then
      WindStr_x_A(i,j) = ((ice_cover_now-FIA_ice_cover)*FIA%WindStr_ocn_x(i,j) + &
                          FIA_ice_cover*FIA%WindStr_x(i,j)) / ice_cover_now
      WindStr_y_A(i,j) = ((ice_cover_now-FIA_ice_cover)*FIA%WindStr_ocn_y(i,j) + &
                          FIA_ice_cover*FIA%WindStr_y(i,j)) / ice_cover_now
    else
      WindStr_x_A(i,j) = FIA%WindStr_x(i,j)
      WindStr_y_A(i,j) = FIA%WindStr_y(i,j)
    endif

    if (ice_free(i,j) <= FIA%ice_free(i,j)) then
      WindStr_x_ocn_A(i,j) = FIA%WindStr_ocn_x(i,j)
      WindStr_y_ocn_A(i,j) = FIA%WindStr_ocn_y(i,j)
    else
      WindStr_x_ocn_A(i,j) = ((ice_free(i,j)-FIA%ice_free(i,j))*FIA%WindStr_x(i,j) + &
                          FIA%ice_free(i,j)*FIA%WindStr_ocn_x(i,j)) / ice_free(i,j)
      WindStr_y_ocn_A(i,j) = ((ice_free(i,j)-FIA%ice_free(i,j))*FIA%WindStr_y(i,j) + &
                          FIA%ice_free(i,j)*FIA%WindStr_ocn_y(i,j)) / ice_free(i,j)
    endif
  enddo ;  enddo

  !$OMP parallel do default(shared) private(weights,I_wts)
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

end subroutine set_wind_stresses_B


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_register_restarts allocates and registers any variables associated
!!      slow ice dynamics and transport that need to be included in the restart files.
subroutine SIS_dyn_trans_register_restarts(mpp_domain, HI, IG, param_file, CS, &
                                      Ice_restart, restart_file)
  type(domain2d),          intent(in) :: mpp_domain !< The ice models' FMS domain type
  type(hor_index_type),    intent(in) :: HI     !< The horizontal index type describing the domain
  type(ice_grid_type),     intent(in) :: IG     !< The sea-ice grid type
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(dyn_trans_CS),      pointer    :: CS     !< The control structure for the SIS_dyn_trans module
  type(restart_file_type), pointer    :: Ice_restart !< The sea ice restart control structure
  character(len=*),        intent(in) :: restart_file !< The ice restart file name

!   This subroutine registers the restart variables associated with the
! the slow ice dynamics and thermodynamics.

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_dyn_trans_register_restarts called with an "//&
                            "associated control structure.")
    return
  endif
  allocate(CS)

  CS%Cgrid_dyn = .false. ; call read_param(param_file, "CGRID_ICE_DYNAMICS", CS%Cgrid_dyn)

  if (CS%Cgrid_dyn) then
    call SIS_C_dyn_register_restarts(mpp_domain, HI, param_file, &
                 CS%SIS_C_dyn_CSp, Ice_restart, restart_file)
  else
    call SIS_B_dyn_register_restarts(mpp_domain, HI, param_file, &
                 CS%SIS_B_dyn_CSp, Ice_restart, restart_file)
  endif
!  call SIS_transport_register_restarts(G, param_file, CS%SIS_transport_CSp, &
!                                       Ice_restart, restart_file)

end subroutine SIS_dyn_trans_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_register_restarts allocates and registers any variables associated
!!      slow ice dynamics and transport that need to be included in the restart files.
subroutine SIS_dyn_trans_read_alt_restarts(CS, G, Ice_restart, &
                                           restart_file, restart_dir)
  type(dyn_trans_CS),      pointer    :: CS  !< The control structure for the SIS_dyn_trans module
  type(SIS_hor_grid_type), intent(in) :: G   !< The horizontal grid type
  type(restart_file_type), pointer    :: Ice_restart !< The sea ice restart control structure
  character(len=*),        intent(in) :: restart_file !< The ice restart file name
  character(len=*),        intent(in) :: restart_dir !< The directory in which to find the restart files

  if (CS%Cgrid_dyn) then
    call SIS_C_dyn_read_alt_restarts(CS%SIS_C_dyn_CSp, G, Ice_restart, &
                                     restart_file, restart_dir)
  endif

end subroutine SIS_dyn_trans_read_alt_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_init initializes ice model data, parameters and diagnostics
!!   associated with the SIS2 dynamics and transport modules.
subroutine SIS_dyn_trans_init(Time, G, IG, param_file, diag, CS, output_dir, Time_init, &
                              slab_ice, specified_ice)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid structure
  type(ice_grid_type),         intent(in)    :: IG   !< The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(dyn_trans_CS),          pointer       :: CS   !< The control structure for the SIS_dyn_trans module
  character(len=*),            intent(in)    :: output_dir !< The directory to use for writing output
  type(time_type),             intent(in)    :: Time_Init !< Starting time of the model integration
  logical,           optional, intent(in)    :: specified_ice !< If present and true, use specified ice.
  logical,           optional, intent(in)    :: slab_ice  !< If true, use the archaic GFDL slab ice dynamics
                                                     !!  and transport.

  ! This include declares and sets the variable "version".
#  include "version_variable.h"
  character(len=40) :: mdl = "SIS_dyn_trans" ! This module's name.
  real :: Time_unit      ! The time unit for ICE_STATS_INTERVAL [s].
  character(len=8) :: nstr
  integer :: n, nLay
  logical :: spec_ice, do_slab_ice
  logical :: debug
  real, parameter :: missing = -1e34

  nLay = IG%NkIce
  spec_ice = .false. ; if (present(specified_ice)) spec_ice = specified_ice
  do_slab_ice = .false. ; if (present(slab_ice)) do_slab_ice = slab_ice

  call callTree_enter("SIS_dyn_trans_init(), SIS_dyn_trans.F90")

  if (.not.associated(CS)) call SIS_error(FATAL, &
      "SIS_dyn_trans_init called with an unassociated control structure.")
  if (CS%module_is_initialized) then
    call SIS_error(WARNING, "SIS_dyn_trans_init called with a control "// &
                            "structure that has already been initialized.")
    return
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, &
     "This module updates the ice momentum and does ice transport.")
  if ( .not.spec_ice ) then
    call get_param(param_file, mdl, "CGRID_ICE_DYNAMICS", CS%Cgrid_dyn, &
                 "If true, use a C-grid discretization of the sea-ice \n"//&
                 "dynamics; if false use a B-grid discretization.", &
                 default=.false.)
    call get_param(param_file, mdl, "DT_ICE_DYNAMICS", CS%dt_ice_dyn, &
                 "The time step used for the slow ice dynamics, including \n"//&
                 "stepping the continuity equation and interactions \n"//&
                 "between the ice mass field and velocities.  If 0 or \n"//&
                 "negative the coupling time step will be used.", &
                 units="seconds", default=-1.0)
    call get_param(param_file, mdl, "MERGED_CONTINUITY", CS%merged_cont, &
                 "If true, update the continuity equations for the ice, snow, \n"//&
                 "and melt pond water together summed across categories, with \n"//&
                 "proportionate fluxes for each part. Otherwise the media are \n"//&
                 "updated separately.", default=.false.)
    call get_param(param_file, mdl, "DO_RIDGING", CS%do_ridging, &
                 "If true, apply a ridging scheme to the convergent ice. \n"//&
                 "Otherwise, ice is compressed proportionately if the \n"//&
                 "concentration exceeds 1.  The original SIS2 implementation \n"//&
                 "is based on work by Torge Martin.", default=.false.)
    call get_param(param_file, mdl, "NSTEPS_ADV", CS%adv_substeps, &
                 "The number of advective iterations for each slow time \n"//&
                 "step.", default=1)

    call get_param(param_file, mdl, "MAX_TRACER_SUBSTEPS", CS%max_nts, &
                 "The maximum number of ice tracer transport steps that \n"//&
                 "can be stored before they are carried out.", default=CS%adv_substeps)

    call get_param(param_file, mdl, "ICEBERG_WINDSTRESS_BUG", CS%berg_windstress_bug, &
                 "If true, use older code that applied an old ice-ocean \n"//&
                 "stress to the icebergs in place of the current air-ocean \n"//&
                 "stress.  This option is here for backward compatibility, \n"//&
                 "but should be avoided.", default=.false.)
    call get_param(param_file, mdl, "WARSAW_SUM_ORDER", CS%Warsaw_sum_order, &
                 "If true, use the order of sums in the Warsaw version of SIS2. \n"//&
                 "The default is the opposite of MERGED_CONTINUITY. \n"//&
                 "This option exists for backward compatibilty but may \n"//&
                 "eventually be obsoleted.", default=.not.CS%merged_cont)
  endif

  call get_param(param_file, mdl, "TIMEUNIT", Time_unit, &
                 "The time unit for ICE_STATS_INTERVAL.", &
                 units="s", default=86400.0)
  call get_param(param_file, mdl, "ICE_STATS_INTERVAL", CS%ice_stats_interval, &
                 "The interval in units of TIMEUNIT between writes of the \n"//&
                 "globally summed ice statistics and conservation checks.", &
                 default=real_to_time(86400.0), timeunit=Time_unit)

  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_SLOW_ICE", CS%debug, &
                 "If true, write out verbose debugging data on the slow ice PEs.", &
                 default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "COLUMN_CHECK", CS%column_check, &
                 "If true, add code to allow debugging of conservation \n"//&
                 "column-by-column.  This does not change answers, but \n"//&
                 "can increase model run time.", default=.false., &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "IMBALANCE_TOLERANCE", CS%imb_tol, &
                 "The tolerance for imbalances to be flagged by COLUMN_CHECK.", &
                 units="nondim", default=1.0e-9, debuggingParam=.true.)
  call get_param(param_file, mdl, "ICE_BOUNDS_CHECK", CS%bounds_check, &
                 "If true, periodically check the values of ice and snow \n"//&
                 "temperatures and thicknesses to ensure that they are \n"//&
                 "sensible, and issue warnings if they are not.  This \n"//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)
  call get_param(param_file, mdl, "VERBOSE", CS%verbose, &
                 "If true, write out verbose diagnostics.", default=.false., &
                 debuggingParam=.true.)

  if (.not.(do_slab_ice .or. spec_ice)) then
    if (CS%Cgrid_dyn) then
      call SIS_C_dyn_init(CS%Time, G, param_file, CS%diag, CS%SIS_C_dyn_CSp, CS%ntrunc)
    else
      call SIS_B_dyn_init(CS%Time, G, param_file, CS%diag, CS%SIS_B_dyn_CSp)
    endif
    call SIS_transport_init(CS%Time, G, param_file, CS%diag, CS%SIS_transport_CSp, &
                            continuity_CSp=CS%continuity_CSp)

    call alloc_cell_average_state_type(CS%CAS, G%HI, IG, CS%SIS_transport_CSp)

    if (CS%merged_cont) then
      call safe_alloc_alloc(CS%mca_step, G%isd, G%ied, G%jsd, G%jed, max(CS%max_nts+1,1))
      call safe_alloc_alloc(CS%uh_step, G%isdB, G%IedB, G%jsd, G%jed, max(CS%max_nts,1))
      call safe_alloc_alloc(CS%vh_step, G%isd, G%ied, G%JsdB, G%JedB, max(CS%max_nts,1))
    endif

  endif

  call SIS_sum_output_init(G, param_file, output_dir, Time_Init, &
                           CS%sum_output_CSp, CS%ntrunc)

  CS%write_ice_stats_time = Time_Init + CS%ice_stats_interval * &
      (1 + (Time - Time_init) / CS%ice_stats_interval)


  ! Ice state diagnostics.
  CS%id_ext = register_diag_field('ice_model', 'EXT', diag%axesT1, Time, &
               'ice modeled', '0 or 1', missing_value=missing)
  CS%id_cn       = register_diag_field('ice_model', 'CN', diag%axesTc, Time, &
               'ice concentration', '0-1', missing_value=missing)
  CS%id_hp       = register_diag_field('ice_model', 'HP', diag%axesT1, Time, &
               'pond thickness', 'm-pond', missing_value=missing) ! mw/new
  CS%id_hs       = register_diag_field('ice_model', 'HS', diag%axesT1, Time, &
               'snow thickness', 'm-snow', missing_value=missing)
  CS%id_tsn      = register_diag_field('ice_model', 'TSN', diag%axesT1, Time, &
               'snow layer temperature', 'C',  missing_value=missing)
  CS%id_hi       = register_diag_field('ice_model', 'HI', diag%axesT1, Time, &
               'ice thickness', 'm-ice', missing_value=missing)
  CS%id_sitimefrac = register_diag_field('ice_model', 'sitimefrac', diag%axesT1, Time, &
               'time fraction of ice cover', '0-1', missing_value=missing)
  CS%id_siconc = register_diag_field('ice_model', 'siconc', diag%axesT1, Time, &
               'ice concentration', '0-1', missing_value=missing)
  CS%id_sithick  = register_diag_field('ice_model', 'sithick', diag%axesT1, Time, &
               'ice thickness', 'm-ice', missing_value=missing)
  CS%id_sivol  = register_diag_field('ice_model', 'sivol', diag%axesT1, Time, &
               'ice volume', 'm-ice', missing_value=missing)
  CS%id_sisnconc = register_diag_field('ice_model', 'sisnconc', diag%axesT1, Time, &
               'snow concentration', '0-1', missing_value=missing)
  CS%id_sisnthick= register_diag_field('ice_model', 'sisnthick', diag%axesT1, Time, &
               'snow thickness', 'm-snow', missing_value=missing)

  CS%id_t_iceav = register_diag_field('ice_model', 'T_bulkice', diag%axesT1, Time, &
               'Volume-averaged ice temperature', 'C', missing_value=missing)
  CS%id_s_iceav = register_diag_field('ice_model', 'S_bulkice', diag%axesT1, Time, &
               'Volume-averaged ice salinity', 'g/kg', missing_value=missing)
  call safe_alloc_ids_1d(CS%id_t, nLay)
  call safe_alloc_ids_1d(CS%id_sal, nLay)
  do n=1,nLay
    write(nstr, '(I4)') n ; nstr = adjustl(nstr)
    CS%id_t(n)   = register_diag_field('ice_model', 'T'//trim(nstr), &
                 diag%axesT1, Time, 'ice layer '//trim(nstr)//' temperature', &
                 'C',  missing_value=missing)
    CS%id_sal(n)   = register_diag_field('ice_model', 'Sal'//trim(nstr), &
               diag%axesT1, Time, 'ice layer '//trim(nstr)//' salinity', &
               'g/kg',  missing_value=missing)
  enddo

  ! Diagnostics that are specific to the C-grid or B-grid dynamics of the ice model
  if (.not.spec_ice) then ; if (CS%Cgrid_dyn) then
    CS%id_fax = register_diag_field('ice_model', 'FA_X', diag%axesCu1, Time, &
               'Air stress on ice on C-grid - x component', 'Pa', &
                missing_value=missing, interp_method='none')
    CS%id_fay = register_diag_field('ice_model', 'FA_Y', diag%axesCv1, Time, &
               'Air stress on ice on C-grid - y component', 'Pa', &
               missing_value=missing, interp_method='none')
  else
    CS%id_fax = register_diag_field('ice_model', 'FA_X', diag%axesB1, Time, &
               'air stress on ice - x component', 'Pa', &
               missing_value=missing, interp_method='none')
    CS%id_fay = register_diag_field('ice_model', 'FA_Y', diag%axesB1, Time, &
               'air stress on ice - y component', 'Pa', &
               missing_value=missing, interp_method='none')
  endif ; endif
  CS%id_mi   = register_diag_field('ice_model', 'MI', diag%axesT1, Time, &
               'ice + snow mass', 'kg/m^2', missing_value=missing)
  CS%id_simass = register_diag_field('ice_model', 'simass', diag%axesT1, Time, &
               'ice mass', 'kg/m^2', missing_value=missing)
  CS%id_sisnmass = register_diag_field('ice_model', 'sisnmass', diag%axesT1, Time, &
               'snow mass', 'kg/m^2', missing_value=missing)
  CS%id_mib  = register_diag_field('ice_model', 'MIB', diag%axesT1, Time, &
               'ice + snow + bergs mass', 'kg/m^2', missing_value=missing)
  CS%id_e2m  = register_diag_field('ice_model','E2MELT' ,diag%axesT1, Time, &
               'heat needed to melt ice', 'J/m^2', missing_value=missing)
!### THESE DIAGNOSTICS DO NOT EXIST YET.
!  CS%id_rdgf    = register_diag_field('ice_model','RDG_FRAC' ,diag%axesT1, Time, &
!               'ridged ice fraction', '0-1', missing_value=missing)
!### THIS DIAGNOSTIC IS MISSING.
!  CS%id_ta    = register_diag_field('ice_model', 'TA', diag%axesT1, Time, &
!            'surface air temperature', 'C', missing_value=missing)

  iceClock4 = mpp_clock_id( '  Ice: slow: dynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClocka = mpp_clock_id( '       slow: ice_dynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClockb = mpp_clock_id( '       slow: comm/cut check ', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClockc = mpp_clock_id( '       slow: diags', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock8 = mpp_clock_id( '  Ice: slow: transport', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock9 = mpp_clock_id( '  Ice: slow: thermodyn diags', flags=clock_flag_default, grain=CLOCK_LOOP )

  call callTree_leave("SIS_dyn_trans_init()")

end subroutine SIS_dyn_trans_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Allocate an array of integer diagnostic arrays and set them to -1, if they are not already allocated
subroutine safe_alloc_ids_1d(ids, nids)
  integer, allocatable, intent(inout) :: ids(:) !< An array of diagnostic IDs to allocate
  integer,              intent(in)    :: nids   !< The number of IDs to allocate

  if (.not.ALLOCATED(ids)) then
    allocate(ids(nids)) ; ids(:) = -1
  endif;
end subroutine safe_alloc_ids_1d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_transport_CS returns a pointer to the SIS_transport_CS type that
!!  the dyn_trans_CS points to.
function SIS_dyn_trans_transport_CS(CS) result(transport_CSp)
  type(dyn_trans_CS),     pointer :: CS    !< The control structure for the SIS_dyn_trans module
  type(SIS_transport_CS), pointer :: transport_CSp !< The SIS_transport_CS type used by SIS_dyn_trans

  transport_CSp => CS%SIS_transport_CSp
end function SIS_dyn_trans_transport_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_transport_CS returns a pointer to the sum_out_CS type that
!! the dyn_trans_CS points to.
function SIS_dyn_trans_sum_output_CS(CS) result(sum_out_CSp)
  type(dyn_trans_CS),   pointer :: CS    !< The control structure for the SIS_dyn_trans module
  type(SIS_sum_out_CS), pointer :: sum_out_CSp !< The SIS_sum_out_CS type used by SIS_dyn_trans

  sum_out_CSp => CS%sum_output_CSp
end function SIS_dyn_trans_sum_output_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_end deallocates memory associated with the dyn_trans_CS type
!! and calls similar routines for subsidiary modules.
subroutine SIS_dyn_trans_end(CS)
  type(dyn_trans_CS), pointer :: CS  !< The control structure for the SIS_dyn_trans module that
                                     !! is dellocated here

  if (CS%Cgrid_dyn) then
    call SIS_C_dyn_end(CS%SIS_C_dyn_CSp)
  else
    call SIS_B_dyn_end(CS%SIS_B_dyn_CSp)
  endif
  call SIS_transport_end(CS%SIS_transport_CSp)
  call dealloc_cell_average_state_type(CS%CAS)

  deallocate(CS)

end subroutine SIS_dyn_trans_end

end module SIS_dyn_trans
