!> Handles the main updates of the ice states at the slower time-scales of the coupling or
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
! scales of the coupling or the interactions with the ocean due to ice dynamics !
! and lateral transport.                                                       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP, CLOCK_ROUTINE
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, read_param, log_param, log_version, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_open_boundary, only : OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_time_manager,  only : time_type, time_type_to_real, real_to_time
use MOM_time_manager,  only : operator(+), operator(-)
use MOM_time_manager,  only : operator(>), operator(*), operator(/), operator(/=)
use MOM_unit_scaling,  only : unit_scale_type
use MOM_EOS,           only : EOS_type, calculate_density_derivs


use SIS_continuity,    only : SIS_continuity_CS, summed_continuity, ice_cover_transport
use SIS_debugging,     only : chksum, Bchksum, hchksum
use SIS_debugging,     only : hchksum_pair, Bchksum_pair, uvchksum
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_dyn_bgrid,     only : SIS_B_dyn_CS, SIS_B_dynamics, SIS_B_dyn_init
use SIS_dyn_bgrid,     only : SIS_B_dyn_register_restarts, SIS_B_dyn_end
use SIS_dyn_cgrid,     only : SIS_C_dyn_CS, SIS_C_dynamics, SIS_C_dyn_init
use SIS_dyn_cgrid,     only : SIS_C_dyn_register_restarts, SIS_C_dyn_end
use SIS_dyn_cgrid,     only : SIS_C_dyn_read_alt_restarts, basal_stress_coeff_C
use SIS_dyn_cgrid,     only : basal_stress_coeff_itd
use SIS_restart,       only : SIS_restart_CS
use SIS_framework,     only : coupler_type_initialized, coupler_type_send_data, safe_alloc
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_ice_diags,     only : ice_state_diags_type, register_ice_state_diagnostics
use SIS_ice_diags,     only : post_ocean_sfc_diagnostics, post_ice_state_diagnostics
use SIS_open_boundary, only : ice_OBC_type, OBC_segment_type
use SIS_sum_output,    only : write_ice_statistics, SIS_sum_output_init, SIS_sum_out_CS
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS
use SIS_transport,     only : SIS_transport_init, SIS_transport_end
use SIS_transport,     only : SIS_transport_CS, adjust_ice_categories, cell_average_state_type
use SIS_transport,     only : alloc_cell_average_state_type, dealloc_cell_average_state_type
use SIS_transport,     only : cell_ave_state_to_ice_state, ice_state_to_cell_ave_state
use SIS_transport,     only : ice_cat_transport, finish_ice_transport
use SIS_types,         only : ocean_sfc_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_types,         only : ice_state_type, IST_chksum, IST_bounds_check
use SIS_utils,         only : get_avg, post_avg, ice_line !, ice_grid_chksum
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs
use slab_ice,          only : slab_ice_advect, slab_ice_dynamics
use ice_bergs,         only : icebergs, icebergs_run, icebergs_init, icebergs_end
use ice_grid,          only : ice_grid_type

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_dynamics_trans, SIS_multi_dyn_trans, update_icebergs, dyn_trans_CS
public :: slab_ice_dyn_trans
public :: SIS_dyn_trans_register_restarts, SIS_dyn_trans_init, SIS_dyn_trans_end
public :: SIS_dyn_trans_read_alt_restarts, stresses_to_stress_mag
public :: SIS_dyn_trans_transport_CS, SIS_dyn_trans_sum_output_CS

!> The control structure for the SIS_dyn_trans module
type dyn_trans_CS ; private
  logical :: Cgrid_dyn    !< If true use a C-grid discretization of the sea-ice dynamics.
  real    :: dt_ice_dyn   !< The time step used for the slow ice dynamics, including
                          !! stepping the continuity equation and interactions
                          !! between the ice mass field and velocities [T ~> s]. If
                          !! 0 or negative, the coupling time step will be used.
  logical :: merged_cont  !< If true, update the continuity equations for the ice, snow,
                          !! and melt pond water together with proportionate fluxes.
                          !! Otherwise the three media are updated separately.
  real    :: dt_advect    !< The time step used for the advecting tracers and masses as
                          !! partitioned by thickness categories when merged_cont it true [T ~> s].
                          !! If 0 or negative, the coupling time step will be used.
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
                          !! of SIS2.  This option exists for backward compatibility
                          !! but may eventually be obsoleted.
  real :: complete_ice_cover !< The fractional ice coverage that is close enough to 1 to be
                          !! complete for the purpose of calculating wind stresses [nondim].

  logical :: debug        !< If true, write verbose checksums for debugging purposes.
  logical :: column_check !< If true, enable the heat check column by column.
  real    :: imb_tol      !< The tolerance for imbalances to be flagged by
                          !! column_check [nondim].
  logical :: bounds_check !< If true, check for sensible values of thicknesses
                          !! temperatures, fluxes, etc.
  logical :: verbose      !< A flag to control the printing of an ice-diagnostic
                          !! message.  When true, this will slow the model down.

  integer :: ntrunc = 0   !< The number of times the velocity has been truncated
                          !! since the last call to write_ice_statistics.

  integer :: n_calls = 0  !< The number of times SIS_dynamics_trans has been called.
  type(time_type) :: ice_stats_interval !< The interval between writes of the
                          !! globally summed ice statistics and conservation checks.
  type(time_type) :: write_ice_stats_time !< The next time to write out the ice statistics.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                                   !! timing of diagnostic output.
  logical :: lemieux_landfast !< If true, use the lemieux landfast ice parameterization.
  logical :: itd_landfast     !< If true, use the probabilistic landfast ice parameterization.

  !>@{ Diagnostic IDs
  integer :: id_fax=-1, id_fay=-1
  !!@}

  type(dyn_state_2d), pointer :: DS2d => NULL()
      !< A simplified 2-d description of the ice state integrated across thickness categories and layers.
  type(cell_average_state_type), pointer :: CAS => NULL()
      !< A structure with ocean-cell averaged masses.
  type(ice_state_diags_type), pointer :: IDs => NULL()
      !< A structure for regulating sea ice state diagnostics
  type(SIS_B_dyn_CS), pointer     :: SIS_B_dyn_CSp => NULL()
      !< Pointer to the control structure for the B-grid dynamics module
  type(SIS_C_dyn_CS), pointer     :: SIS_C_dyn_CSp => NULL()
      !< Pointer to the control structure for the C-grid dynamics module
  type(SIS_transport_CS), pointer :: SIS_transport_CSp => NULL()
      !< Pointer to the control structure for the ice transport module
  type(SIS_continuity_CS),    pointer :: continuity_CSp => NULL()
      !< The control structure for the SIS continuity module
  type(SIS_continuity_CS),    pointer :: cover_trans_CSp => NULL()
      !< The control structure for ice cover transport by the SIS continuity module
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()
     !< Pointer to the control structure for the summed diagnostics module
  logical :: module_is_initialized = .false. !< If true, this module has been initialized.
end type dyn_trans_CS

!> A simplified 2-d description of the ice state integrated across thickness categories and layers.
type, public :: dyn_state_2d ; private
  integer :: max_nts      !< The maximum number of transport steps that can be stored
                          !! before they are carried out.
  integer :: nts = 0      !< The number of accumulated transport steps since the last update.
  real :: ridge_rate_count !< The number of contributions to avg_ridge_rate

  real, allocatable, dimension(:,:) :: avg_ridge_rate !< The time average ridging rate in [T-1 ~> s-1].

  real, allocatable, dimension(:,:) :: mi_sum !< The total mass of ice per unit total area [R Z ~> kg m-2].
  real, allocatable, dimension(:,:) :: ice_cover !< The fractional ice coverage, summed across all
                          !! thickness categories [nondim], between 0 & 1.
  real, allocatable, dimension(:,:) :: u_ice_B  !< The pseudo-zonal ice velocity along the
                !! along the grid directions on a B-grid [L T-1 ~> m s-1].
                !! All thickness categories are assumed to have the same velocities.
  real, allocatable, dimension(:,:) :: v_ice_B  !< The pseudo-meridional ice velocity along the
                !! along the grid directions on a B-grid [L T-1 ~> m s-1].
  real, allocatable, dimension(:,:) :: u_ice_C  !< The pseudo-zonal ice velocity along the
                !! along the grid directions on a C-grid [L T-1 ~> m s-1].
                !! All thickness categories are assumed to have the same velocities.
  real, allocatable, dimension(:,:) :: v_ice_C  !< The pseudo-meridional ice velocity along the
                !! along the grid directions on a C-grid [L T-1 ~> m s-1].
  real, allocatable, dimension(:,:,:) :: mca_step !< The total mass per unit total area of snow, ice
                          !! and pond water summed across thickness categories in a cell, after each
                          !! transportation substep, with a 0 starting 3rd index [R Z ~> kg m-2].
  real, allocatable, dimension(:,:,:) :: uh_step !< The total zonal mass fluxes during each
                          !! transportation substep [R Z L2 T-1 ~> kg s-1].
  real, allocatable, dimension(:,:,:) :: vh_step !< The total meridional mass fluxes during each
                          !! transportation substep [R Z L2 T-1 ~> kg s-1].

end type dyn_state_2d

!>@{ CPU time clock IDs
integer :: iceClock4, iceClock8, iceClock9, iceClocka, iceClockb, iceClockc
!!@}

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> update_icebergs calls icebergs_run and offers diagnostics of some of the
!! iceberg fields that might drive the sea ice or ocean
subroutine update_icebergs(IST, OSS, IOF, FIA, icebergs_CS, dt_slow, G, US, IG, CS)
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
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module

  real, dimension(SZI_(G),SZJ_(G))   :: &
    hi_avg            ! The area-weighted average ice thickness in the units used for icebergs [m].
  real, dimension(G%isc:G%iec, G%jsc:G%jec)   :: &
    Temp_sfc, &       ! A local copy of the surface temperature in the units used for icebergs [degC]
    Saln_sfc, &       ! A local copy of the surface salinity in the units used for icebergs [gSalt kg-1]
    calving, &        ! A local copy of the calving rate in the units used for icebergs [kg m-2 s-1]
    calving_hflx, &   ! A local copy of the calving heat flux in the units used for icebergs [W m-2]
    windstr_x, &      ! The area-weighted average ice thickness in the units used for icebergs [Pa].
    windstr_y         ! The area-weighted average ice thickness in the units used for icebergs [Pa].
  real, dimension(G%isc-2:G%iec+1, G%jsc-1:G%jec+1)   :: &
    u_ice_C, &        ! The C-grid zonal ice velocity in the units used for icebergs [m s-1].
    u_ocn_C           ! The C-grid zonal ocean velocity in the units used for icebergs [m s-1].
  real, dimension(G%isc-1:G%iec+1, G%jsc-2:G%jec+1)   :: &
    v_ice_C, &        ! The C-grid meridional ice velocity in the units used for icebergs [m s-1].
    v_ocn_C           ! The C-grid meridional ocean velocity in the units used for icebergs [m s-1].
  real, dimension(G%isc-1:G%iec+1, G%jsc-1:G%jec+1)   :: &
    sea_lev, &        ! Sea level anomalies in the units used for icebergs [m].
    u_ice_B, &        ! The B-grid zonal ice velocity in the units used for icebergs [m s-1].
    u_ocn_B, &        ! The B-grid zonal ocean velocity in the units used for icebergs [m s-1].
    v_ice_B, &        ! The B-grid meridional ice velocity in the units used for icebergs [m s-1].
    v_ocn_B           ! The B-grid meridional ocean velocity in the units used for icebergs [m s-1].
  real :: rho_ice     ! The nominal density of sea ice [R ~> kg m-3].
  real :: H_to_m_ice  ! The specific volume of ice times the conversion factor from the SIS2
                      ! thickness units to those used by the icebergs [m R-1 Z-1 ~> m3 kg-1].
  integer :: stress_stagger
  integer :: i, j, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  H_to_m_ice = US%Z_to_m / rho_ice
  call get_avg(IST%mH_ice, IST%part_size(:,:,1:), hi_avg, wtd=.true.)
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    hi_avg(i,j) = hi_avg(i,j) * H_to_m_Ice
  enddo ; enddo
  do j=jsc,jec ; do i=isc,iec
    calving(i,j) = US%RZ_T_to_kg_m2s*FIA%calving(i,j)
    calving_hflx(i,j) = US%QRZ_T_to_W_m2*FIA%calving_hflx(i,j)
    Saln_sfc(i,j) = US%S_to_ppt*OSS%s_surf(i,j)
    Temp_sfc(i,j) = US%C_to_degC*OSS%SST_C(i,j)
  enddo ; enddo
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    sea_lev(i,j) = US%Z_to_m*OSS%sea_lev(i,j)
  enddo ; enddo

  if (CS%berg_windstress_bug) then
    ! This code reproduces a long-standing bug, in that the old ice-ocean
    ! stresses are being passed in place of the wind stresses on the icebergs.
    do j=jsc,jec ; do i=isc,iec
      windstr_x(i,j) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*IOF%flux_u_ocn(i,j)
      windstr_y(i,j) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*IOF%flux_v_ocn(i,j)
    enddo ; enddo
    stress_stagger = IOF%flux_uv_stagger
  else
    do j=jsc,jec ; do i=isc,iec
      windstr_x(i,j) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*FIA%WindStr_ocn_x(i,j)
      windstr_y(i,j) = US%RZ_T_to_kg_m2s*US%L_T_to_m_s*FIA%WindStr_ocn_y(i,j)
    enddo ; enddo
    stress_stagger = AGRID
  endif

  if (IST%Cgrid_dyn) then
    do j=jsc-1,jec+1 ; do I=isc-2,iec+1
      u_ice_C(I,j) = US%L_T_to_m_s*IST%u_ice_C(I,j) ; u_ocn_C(I,j) = US%L_T_to_m_s*OSS%u_ocn_C(I,j)
    enddo ; enddo
    do J=jsc-2,jec+1 ; do i=isc-1,iec+1
      v_ice_C(i,J) = US%L_T_to_m_s*IST%v_ice_C(i,J) ; v_ocn_C(i,J) = US%L_T_to_m_s*OSS%v_ocn_C(i,J)
    enddo ; enddo
    call icebergs_run( icebergs_CS, CS%Time, calving, &
            u_ocn_C, v_ocn_C, u_ice_C, v_ice_C, windstr_x, windstr_y, &
            sea_lev, Temp_sfc, calving_hflx, FIA%ice_cover(isc-1:iec+1,jsc-1:jec+1), &
            hi_avg(isc-1:iec+1,jsc-1:jec+1), stagger=CGRID_NE, &
            stress_stagger=stress_stagger, sss=Saln_sfc, &
            mass_berg=IOF%mass_berg, ustar_berg=IOF%ustar_berg, &
            area_berg=IOF%area_berg )
  else
    do J=jsc-1,jec+1 ; do I=isc-1,iec+1
      u_ice_B(I,J) = US%L_T_to_m_s*IST%u_ice_B(I,J) ; u_ocn_B(I,J) = US%L_T_to_m_s*OSS%u_ocn_B(I,J)
      v_ice_B(I,J) = US%L_T_to_m_s*IST%v_ice_B(I,J) ; v_ocn_B(I,J) = US%L_T_to_m_s*OSS%v_ocn_B(I,J)
    enddo ; enddo
    call icebergs_run( icebergs_CS, CS%Time, calving, &
            u_ocn_B, v_ocn_B, u_ice_B, v_ice_B, windstr_x, windstr_y, &
            sea_lev, Temp_sfc, calving_hflx, FIA%ice_cover(isc-1:iec+1,jsc-1:jec+1), &
            hi_avg(isc-1:iec+1,jsc-1:jec+1), stagger=BGRID_NE, &
            stress_stagger=stress_stagger, sss=Saln_sfc, &
            mass_berg=IOF%mass_berg, ustar_berg=IOF%ustar_berg, &
            area_berg=IOF%area_berg )
  endif

  do j=jsc,jec ; do i=isc,iec
    FIA%calving(i,j) = US%kg_m2s_to_RZ_T*calving(i,j)
    FIA%calving_hflx(i,j) = US%W_m2_to_QRZ_T*calving_hflx(i,j)
  enddo ; enddo
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
subroutine SIS_dynamics_trans(IST, OSS, FIA, IOF, dt_slow, CS, icebergs_CS, G, US, IG, tracer_CSp, OBC)
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow !< The slow ice dynamics timestep [T ~> s].
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module
  type(icebergs),             pointer       :: icebergs_CS !< A control structure for the iceberg model.
  type(SIS_tracer_flow_control_CS), pointer :: tracer_CSp !< The structure for controlling calls to
                                                   !! auxiliary ice tracer packages
  type(ice_OBC_type),         pointer       :: OBC  !< Open boundary structure.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    mi_sum, &           ! Masses of ice per unit total area [R Z ~> kg m-2].
    misp_sum, &         ! Combined mass of snow, ice and melt pond water per unit total area [R Z ~> kg m-2].
    ice_free, &         ! The fractional open water [nondim], between 0 & 1.
    ice_cover           ! The fractional ice coverage, summed across all
                        ! thickness categories [nondim], between 0 & 1.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    WindStr_x_B, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_B, &      ! averaged over the ice categories on a B-grid [R Z L T-2 ~> Pa].
    WindStr_x_ocn_B, &  ! Zonal wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    WindStr_y_ocn_B, &  ! Meridional wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    str_x_ice_ocn_B, &  ! Zonal ice-ocean stress on a B-grid [R Z L T-2 ~> Pa].
    str_y_ice_ocn_B     ! Meridional ice-ocean stress on a B-grid [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G))  :: &
    WindStr_x_Cu, &   ! Zonal wind stress averaged over the ice categories on C-grid u-points [R Z L T-2 ~> Pa].
    WindStr_x_ocn_Cu, & ! Zonal wind stress on the ice-free ocean on C-grid u-points [R Z L T-2 ~> Pa].
    str_x_ice_ocn_Cu   ! Zonal ice-ocean stress on C-grid u-points [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G))  :: &
    WindStr_y_Cv, &   ! Meridional wind stress averaged over the ice categories on C-grid v-points [R Z L T-2 ~> Pa].
    WindStr_y_ocn_Cv, & ! Meridional wind stress on the ice-free ocean on C-grid v-points [R Z L T-2 ~> Pa].
    str_y_ice_ocn_Cv  ! Meridional ice-ocean stress on C-grid v-points [R Z L T-2 ~> Pa].

  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBx ! A temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBy ! A temporary array for diagnostics.
  real :: ps_vel   ! The fractional thickness category coverage at a velocity point.

  type(time_type) :: Time_cycle_start ! The model's time at the start of an advective cycle.
  real :: dt_slow_dyn  ! The slow dynamics timestep [T ~> s].
  real :: dt_slow_dyn_sec ! The slow dynamics timestep [s].
  real :: dt_adv_cycle ! The length of the advective cycle timestep [T ~> s].
  real :: wt_new, wt_prev ! Weights in an average.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rdg_rate  ! A ridging rate [T-1 ~> s-1], calculated from the strain rates in the dynamics.
  type(dyn_state_2d), pointer :: DS2d => NULL()  ! A simplified 2-d description of the ice state
                                                 ! integrated across thickness categories and layers.
  integer :: i, j, k, n, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed
  integer :: ndyn_steps, nds ! The number of dynamic steps.
  integer :: nadv_cycle, nac ! The number of tracer advective cycles in this call.

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%n_calls = CS%n_calls + 1
  IOF%stress_count = 0

  DS2d => CS%DS2d

  ndyn_steps = 1 ; nadv_cycle = 1
  if ((CS%dt_advect > 0.0) .and. (CS%dt_advect < dt_slow)) &
    nadv_cycle = max(CEILING(dt_slow/CS%dt_advect - 1e-9), 1)
  dt_adv_cycle = dt_slow / real(nadv_cycle)

  if ((CS%dt_ice_dyn > 0.0) .and. (CS%dt_ice_dyn < dt_adv_cycle)) &
    ndyn_steps = max(CEILING(dt_adv_cycle/CS%dt_ice_dyn - 1e-6), 1)
  dt_slow_dyn = dt_adv_cycle / real(ndyn_steps)
  dt_slow_dyn_sec = US%T_to_s*dt_slow_dyn

  do nac=1,nadv_cycle
    Time_cycle_start = CS%Time - real_to_time((nadv_cycle-(nac-1))*US%T_to_s*dt_adv_cycle)

    if (CS%merged_cont) then
      ! Convert the category-resolved ice state into the simplified 2-d ice state.
      ! This should be called after a thermodynamic step or if ice_transport was called.
      call convert_IST_to_simple_state(IST, CS%DS2d, CS%CAS, G, US, IG, CS)

      ! Update the category-merged dynamics and use the merged continuity equation.
      call SIS_merged_dyn_cont(OSS, FIA, IOF, CS%DS2d, IST, dt_adv_cycle, Time_cycle_start, &
                               G, US, IG, CS, OBC)

      ! Complete the category-resolved mass and tracer transport and update the ice state type.
      call complete_IST_transport(CS%DS2d, CS%CAS, IST, dt_adv_cycle, G, US, IG, CS)

    else !  (.not.CS%merged_cont)

      do nds=1,ndyn_steps

        call cpu_clock_begin(iceClock4)
        ! The code timed by iceClock4 is the non-merged-cont equivalent of convert_IST_to_simple_state.

        ! Convert the category-resolved ice state into the simplified 2-d ice state.
        ! This should be called after a thermodynamic step or if ice_transport was called.
        if (DS2d%nts == 0) then  ! (This is always true.)
          misp_sum(:,:) = 0.0 ; mi_sum(:,:) = 0.0 ; ice_cover(:,:) = 0.0
          !$OMP parallel do default(shared)
          do j=jsd,jed ; do k=1,ncat ; do i=isd,ied
            misp_sum(i,j) = misp_sum(i,j) + IST%part_size(i,j,k) * &
                            (IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k))
            mi_sum(i,j) = mi_sum(i,j) + IST%mH_ice(i,j,k) * IST%part_size(i,j,k)
            ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
          enddo ; enddo ; enddo
          do j=jsd,jed ; do i=isd,ied
            misp_sum(i,j) = misp_sum(i,j) + mi_sum(i,j)
            ice_free(i,j) = IST%part_size(i,j,0)
          enddo ; enddo

          !  Determine the whole-cell averaged mass of snow and ice.
          call ice_state_to_cell_ave_state(IST, G, US, IG, CS%SIS_transport_CSp, CS%CAS)
        endif
        if (.not.CS%Warsaw_sum_order) then
          do j=jsd,jed ; do i=isd,ied ; ice_free(i,j) = max(1.0 - ice_cover(i,j), 0.0) ; enddo ; enddo
        endif
        call cpu_clock_end(iceClock4)

        !
        ! Dynamics - update ice velocities.
        !

        call enable_SIS_averaging(dt_slow_dyn_sec, Time_cycle_start + real_to_time(nds*dt_slow_dyn_sec), CS%diag)

        ! In the dynamics code, only the ice velocities are changed, and the ice-ocean
        ! stresses are calculated.  The gravity wave dynamics (i.e. the continuity
        ! equation) are not included in the dynamics.  All of the thickness categories
        ! are merged together.

        call cpu_clock_begin(iceClock4)
        ! The code timed by iceClock4 is the non-merged-cont equivalent of SIS_merged_dyn_cont.
        if (CS%Cgrid_dyn) then

          ! Correct the wind stresses for changes in the fractional ice-coverage and set
          ! the wind stresses on the ice and the open ocean for a C-grid staggering.
          ! This block of code must be executed if ice_cover and ice_free or the various wind
          ! stresses were updated.
          call set_wind_stresses_C(FIA, ice_cover, ice_free, WindStr_x_Cu, WindStr_y_Cv, &
                                   WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, G, US, CS%complete_ice_cover, OBC)

          if (CS%lemieux_landfast) then
            call basal_stress_coeff_C(G, mi_sum, ice_cover, OSS%sea_lev, CS%SIS_C_dyn_CSp)
          elseif (CS%itd_landfast) then
            call basal_stress_coeff_itd(G, IG, IST, OSS%sea_lev, CS%SIS_C_dyn_CSp)
          endif

          if (CS%debug) then
            call uvchksum("Before SIS_C_dynamics [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G, scale=US%L_T_to_m_s)
            call hchksum(ice_free, "ice_free before SIS_C_dynamics", G%HI)
            call hchksum(misp_sum, "misp_sum before SIS_C_dynamics", G%HI, scale=US%RZ_to_kg_m2)
            call hchksum(mi_sum, "mi_sum before SIS_C_dynamics", G%HI, scale=US%RZ_to_kg_m2)
            call hchksum(OSS%sea_lev, "sea_lev before SIS_C_dynamics", G%HI, haloshift=1, scale=US%Z_to_m)
            call hchksum(ice_cover, "ice_cover before SIS_C_dynamics", G%HI, haloshift=1)
            call uvchksum("[uv]_ocn before SIS_C_dynamics", OSS%u_ocn_C, OSS%v_ocn_C, G, halos=1, scale=US%L_T_to_m_s)
            call uvchksum("WindStr_[xy] before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G, &
                          halos=1, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    !        call hchksum_pair("WindStr_[xy]_A before SIS_C_dynamics", WindStr_x_A, WindStr_y_A, G, halos=1)
          endif

          !### Ridging needs to be added with C-grid dynamics.
          call cpu_clock_begin(iceClocka)
          if (CS%do_ridging) rdg_rate(:,:) = 0.0
          if (CS%Warsaw_sum_order) then
            call SIS_C_dynamics(1.0-ice_free(:,:), misp_sum, mi_sum, IST%u_ice_C, IST%v_ice_C, &
                                OSS%u_ocn_C, OSS%v_ocn_C, WindStr_x_Cu, WindStr_y_Cv, OSS%sea_lev, &
                                str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, dt_slow_dyn, G, US, CS%SIS_C_dyn_CSp)
          else
            call SIS_C_dynamics(ice_cover, misp_sum, mi_sum, IST%u_ice_C, IST%v_ice_C, &
                                OSS%u_ocn_C, OSS%v_ocn_C, WindStr_x_Cu, WindStr_y_Cv, OSS%sea_lev, &
                                str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, dt_slow_dyn, G, US, CS%SIS_C_dyn_CSp)
          endif
          call cpu_clock_end(iceClocka)

          if (CS%debug) call uvchksum("After ice_dynamics [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G, scale=US%L_T_to_m_s)

          call cpu_clock_begin(iceClockb)
          call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
          call pass_vector(str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, G%Domain, stagger=CGRID_NE)
          call cpu_clock_end(iceClockb)

          ! Dynamics diagnostics
          call cpu_clock_begin(iceClockc)
          if (CS%id_fax>0) call post_data(CS%id_fax, WindStr_x_Cu, CS%diag)
          if (CS%id_fay>0) call post_data(CS%id_fay, WindStr_y_Cv, CS%diag)

          if (CS%debug) call uvchksum("Before set_ocean_top_stress_Cgrid [uv]_ice_C", &
                                      IST%u_ice_C, IST%v_ice_C, G, scale=US%L_T_to_m_s)

          ! Store all mechanical ocean forcing.
          if (CS%Warsaw_sum_order) then
            call set_ocean_top_stress_Cgrid(IOF, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, &
                                            str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, IST%part_size, G, US, IG, OBC)
          else
            call set_ocean_top_stress_C2(IOF, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, &
                                         str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, ice_free, ice_cover, G, US, OBC)
          endif

          call cpu_clock_end(iceClockc)

        else ! B-grid dynamics.

          ! Correct the wind stresses for changes in the fractional ice-coverage and set
          ! the wind stresses on the ice and the open ocean for a C-grid staggering.
          ! This block of code must be executed if ice_cover and ice_free or the various wind
          ! stresses were updated.
          call set_wind_stresses_B(FIA, ice_cover, ice_free, WindStr_x_B, WindStr_y_B, &
                                   WindStr_x_ocn_B, WindStr_y_ocn_B, G, US, CS%complete_ice_cover)

          if (CS%debug) then
            call Bchksum_pair("[uv]_ice_B before dynamics", IST%u_ice_B, IST%v_ice_B, G, scale=US%L_T_to_m_s)
            call hchksum(ice_free, "ice_free before ice_dynamics", G%HI)
            call hchksum(misp_sum, "misp_sum before ice_dynamics", G%HI, scale=US%RZ_to_kg_m2)
            call hchksum(mi_sum, "mi_sum before ice_dynamics", G%HI, scale=US%RZ_to_kg_m2)
            call hchksum(OSS%sea_lev, "sea_lev before ice_dynamics", G%HI, haloshift=1, scale=US%Z_to_m)
            call Bchksum_pair("[uv]_ocn before ice_dynamics", OSS%u_ocn_B, OSS%v_ocn_B, G, scale=US%L_T_to_m_s)
            call Bchksum_pair("WindStr_[xy]_B before ice_dynamics", WindStr_x_B, WindStr_y_B, G, halos=1, &
                              scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
          endif

          call cpu_clock_begin(iceClocka)
          if (CS%do_ridging) rdg_rate(:,:) = 0.0
          if (CS%Warsaw_sum_order) then
            call SIS_B_dynamics(1.0-ice_free(:,:), misp_sum, mi_sum, IST%u_ice_B, IST%v_ice_B, &
                                OSS%u_ocn_B, OSS%v_ocn_B, WindStr_x_B, WindStr_y_B, OSS%sea_lev, &
                                str_x_ice_ocn_B, str_y_ice_ocn_B, CS%do_ridging, &
                                rdg_rate, dt_slow_dyn, G, US, CS%SIS_B_dyn_CSp)
          else
            call SIS_B_dynamics(ice_cover, misp_sum, mi_sum, IST%u_ice_B, IST%v_ice_B, &
                                OSS%u_ocn_B, OSS%v_ocn_B, WindStr_x_B, WindStr_y_B, OSS%sea_lev, &
                                str_x_ice_ocn_B, str_y_ice_ocn_B, CS%do_ridging, &
                                rdg_rate, dt_slow_dyn, G, US, CS%SIS_B_dyn_CSp)
          endif
          call cpu_clock_end(iceClocka)

          if (CS%debug) call Bchksum_pair("After dynamics [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G, scale=US%L_T_to_m_s)

          call cpu_clock_begin(iceClockb)
          call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
          call cpu_clock_end(iceClockb)

          ! Dynamics diagnostics
          call cpu_clock_begin(iceClockc)
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

          if (CS%debug) call Bchksum_pair("Before set_ocean_top_stress_Bgrid [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, &
                                          G, scale=US%L_T_to_m_s)
          ! Store all mechanical ocean forcing.
          if (CS%Warsaw_sum_order) then
            call set_ocean_top_stress_Bgrid(IOF, WindStr_x_ocn_B, WindStr_y_ocn_B, &
                                            str_x_ice_ocn_B, str_y_ice_ocn_B, IST%part_size, G, US, IG)
          else
            call set_ocean_top_stress_B2(IOF, WindStr_x_ocn_B, WindStr_y_ocn_B, &
                                            str_x_ice_ocn_B, str_y_ice_ocn_B, ice_free, ice_cover, G, US)
          endif
          call cpu_clock_end(iceClockc)

          ! Convert the velocities to C-grid points for use in transport.
          do j=jsc,jec ; do I=isc-1,iec
            IST%u_ice_C(I,j) = 0.5 * ( IST%u_ice_B(I,J-1) + IST%u_ice_B(I,J) )
          enddo ; enddo
          do J=jsc-1,jec ; do i=isc,iec
            IST%v_ice_C(i,J) = 0.5 * ( IST%v_ice_B(I-1,J) + IST%v_ice_B(I,J) )
          enddo ; enddo
        endif ! End of B-grid dynamics

        if (CS%do_ridging) then ! Accumulate the time-average ridging rate.
          DS2d%ridge_rate_count = DS2d%ridge_rate_count + 1.
          wt_new = 1.0 / DS2d%ridge_rate_count ; wt_prev = 1.0 - wt_new
          do j=jsc,jec ; do i=isc,iec
            DS2d%avg_ridge_rate(i,j) = wt_new * rdg_rate(i,j) + wt_prev * DS2d%avg_ridge_rate(i,j)
          enddo ; enddo
        endif

        call cpu_clock_end(iceClock4)

      enddo ! nds=1,ndyn_steps

      ! Do ice mass transport and related tracer transport.  This updates the category-decomposed ice state.
      call cpu_clock_begin(iceClock8)
      ! The code timed by iceClock8 is the non-merged_cont equivalent to complete_IST_transport.
      if (CS%debug) call uvchksum("Before ice_transport [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G, scale=US%L_T_to_m_s)
      call enable_SIS_averaging(dt_slow_dyn_sec, Time_cycle_start + real_to_time(nds*dt_slow_dyn_sec), CS%diag)

      call ice_cat_transport(CS%CAS, IST%TrReg, dt_slow_dyn, CS%adv_substeps, G, US, IG, CS%SIS_transport_CSp, &
                             uc=IST%u_ice_C, vc=IST%v_ice_C)

      if (DS2d%nts==0) then
        if (CS%do_ridging) then
          call finish_ice_transport(CS%CAS, IST, IST%TrReg, G, US, IG, dt_slow, CS%SIS_transport_CSp, &
               !                                    rdg_rate=DS2d%avg_ridge_rate)
                                    rdg_rate=IST%rdg_rate)
          DS2d%ridge_rate_count = 0. ; DS2d%avg_ridge_rate(:,:) = 0.0
        else
          call finish_ice_transport(CS%CAS, IST, IST%TrReg, G, US, IG, dt_slow,CS%SIS_transport_CSp)
        endif
      endif
      call cpu_clock_end(iceClock8)

    endif ! (.not.CS%merged_cont)

    if (CS%column_check .and. (DS2d%nts==0)) &
      call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
                                message="      Post_transport")! , check_column=.true.)

  enddo ! nac = 1,nadv_cycle

  ! Finalized the streses for use by the ocean.
  call finish_ocean_top_stresses(IOF, G)

  ! Do diagnostics and update some information for the atmosphere.
  call ice_state_cleanup(IST, OSS, IOF, dt_slow, G, US, IG, CS, tracer_CSp)

end subroutine SIS_dynamics_trans


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_multi_dyn_trans makes the calls to do ice dynamics and mass and tracer transport as
!! appropriate for a dynamic and advective update cycle with multiple calls.
subroutine SIS_multi_dyn_trans(IST, OSS, FIA, IOF, dt_slow, CS, icebergs_CS, G, US, IG, tracer_CSp, &
                               OBC, start_cycle, end_cycle, cycle_length)
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow !< The slow ice dynamics timestep [T ~> s].
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US    !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module
  type(icebergs),             pointer       :: icebergs_CS !< A control structure for the iceberg model.
  type(SIS_tracer_flow_control_CS), pointer :: tracer_CSp !< The structure for controlling calls to
                                                   !! auxiliary ice tracer packages
  type(ice_OBC_type),         pointer       :: OBC !< Open boundary structure.
  logical,          optional, intent(in)    :: start_cycle !< This indicates whether this call is to be
                                                   !! treated as the first call to SIS_multi_dyn_trans
                                                   !! in a time-stepping cycle; missing is like true.
  logical,          optional, intent(in)    :: end_cycle   !< This indicates whether this call is to be
                                                   !! treated as the last call to SIS_multi_dyn_trans
                                                   !! in a time-stepping cycle; missing is like true.
  real,             optional, intent(in)    :: cycle_length !< The duration of a coupled time stepping cycle [s].

  ! Local variables
  real :: dt_adv_cycle ! The length of the advective cycle timestep [T ~> s].
  real :: dt_diags     ! The length of time over which the diagnostics are valid [T ~> s].
  type(time_type) :: Time_cycle_start ! The model's time at the start of an advective cycle.
  integer :: nadv_cycle, nac ! The number of tracer advective cycles within this call.
  logical :: cycle_start, cycle_end, end_of_cycle

  CS%n_calls = CS%n_calls + 1
  IOF%stress_count = 0

  cycle_start = .true. ; if (present(start_cycle)) cycle_start = start_cycle
  cycle_end = .true. ; if (present(end_cycle)) cycle_end = end_cycle
  dt_diags = dt_slow ; if (present(cycle_length)) dt_diags = US%s_to_T*cycle_length

  if (.not.CS%merged_cont) call SIS_error(FATAL, &
          "SIS_multi_dyn_trans should not be called unless MERGED_CONTINUITY=True.")

  nadv_cycle = 1
  if ((CS%dt_advect > 0.0) .and. (CS%dt_advect < dt_slow)) &
    nadv_cycle = max(CEILING(dt_slow/CS%dt_advect - 1e-9), 1)
  dt_adv_cycle = dt_slow / real(nadv_cycle)

  do nac=1,nadv_cycle
    ! Convert the category-resolved ice state into the simplified 2-d ice state.
    ! This should be called after a thermodynamic step or if ice_transport was called.
    if ((nac > 1) .or. cycle_start) &
      call convert_IST_to_simple_state(IST, CS%DS2d, CS%CAS, G, US, IG, CS)

    ! Update the category-merged dynamics and use the merged continuity equation.
    ! This could be called as many times as necessary.
    Time_cycle_start = CS%Time - real_to_time((nadv_cycle-(nac-1))*US%T_to_s*dt_adv_cycle)
    end_of_cycle = (nac < nadv_cycle) .or. cycle_end
    call SIS_merged_dyn_cont(OSS, FIA, IOF, CS%DS2d, IST, dt_adv_cycle, Time_cycle_start, G, US, IG, CS, &
                             OBC, end_call=end_of_cycle)

    ! Complete the category-resolved mass and tracer transport and update the ice state type.
    ! This must be done before the next thermodynamic step.
    if (end_of_cycle) &
      call complete_IST_transport(CS%DS2d, CS%CAS, IST, dt_adv_cycle, G, US, IG, CS)

    if (CS%column_check .and. IST%valid_IST) &  ! This is just here from early debugging exercises,
      call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
                                message="      Post_transport")! , check_column=.true.)

  enddo ! nac=0,nadv_cycle-1
  ! This must be done before returning control to the ocean, but it does not require
  ! that complete_IST_transport be called.
  call finish_ocean_top_stresses(IOF, G, CS%DS2d, IG)

  ! This must be done before returning control to the atmosphere and before writing any diagnostics.
  if (cycle_end) &
    call ice_state_cleanup(IST, OSS, IOF, dt_diags, G, US, IG, CS, tracer_CSp)

end subroutine SIS_multi_dyn_trans

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Complete the category-resolved mass and tracer transport and update the ice state type.
subroutine complete_IST_transport(DS2d, CAS, IST, dt_adv_cycle, G, US, IG, CS)
  type(ice_state_type),          intent(inout) :: IST !< A type describing the state of the sea ice
  type(dyn_state_2d),            intent(inout) :: DS2d !< A simplified 2-d description of the ice state
                                                   !! integrated across thickness categories and layers.
  type(cell_average_state_type), intent(inout) :: CAS !< A structure with ocean-cell averaged masses.
  real,                          intent(in)    :: dt_adv_cycle !< The time since the last IST transport [T ~> s].
  type(SIS_hor_grid_type),       intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),         intent(in)    :: US    !< A structure with unit conversion factors
  type(ice_grid_type),           intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),            pointer       :: CS  !< The control structure for the SIS_dyn_trans module

  integer :: i, j, k, isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call cpu_clock_begin(iceClock8)
  ! Do the transport of mass and tracers by category and vertical layer.
  call ice_cat_transport(CS%CAS, IST%TrReg, dt_adv_cycle, DS2d%nts, G, US, IG, &
                         CS%SIS_transport_CSp, mca_tot=DS2d%mca_step(:,:,0:DS2d%nts), &
                         uh_tot=DS2d%uh_step(:,:,1:DS2d%nts), vh_tot=DS2d%vh_step(:,:,1:DS2d%nts))
  ! Convert the cell-averaged state back to the ice-state type, adjusting the
  ! category mass distributions, doing ridging, and updating the partition sizes.
  if (CS%do_ridging) then
    call finish_ice_transport(CS%CAS, IST, IST%TrReg, G, US, IG, dt_adv_cycle, CS%SIS_transport_CSp, &
         !                              rdg_rate=DS2d%avg_ridge_rate)
                              rdg_rate=IST%rdg_rate)
    DS2d%ridge_rate_count = 0. ; DS2d%avg_ridge_rate(:,:) = 0.0
  else
    call finish_ice_transport(CS%CAS, IST, IST%TrReg, G, US, IG, dt_adv_cycle, CS%SIS_transport_CSp)
  endif
  DS2d%nts = 0 ! There is no outstanding transport to be done and IST is up-to-date.

  ! Copy the velocities back to the ice state type
  if (CS%Cgrid_dyn) then
    do j=jsd,jed ; do I=IsdB,IedB ; IST%u_ice_C(I,j) = DS2d%u_ice_C(I,j) ; enddo ; enddo
    do J=JsdB,JedB ; do i=isd,ied ; IST%v_ice_C(i,J) = DS2d%v_ice_C(i,J) ; enddo ; enddo
  else
    do J=JsdB,JedB ; do I=IsdB,IedB
      IST%u_ice_B(I,J) = DS2d%u_ice_B(I,J) ; IST%v_ice_B(I,J) = DS2d%v_ice_B(I,J)
    enddo ; enddo
  endif

  IST%valid_IST = .true.
  call cpu_clock_end(iceClock8)

end subroutine complete_IST_transport

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Do final checks to set a consistent ice state and write diagnostics as appropriate.
subroutine ice_state_cleanup(IST, OSS, IOF, dt_slow, G, US, IG, CS, tracer_CSp)
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow !< The slow ice dynamics timestep [T ~> s].
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module
  type(SIS_tracer_flow_control_CS), optional, pointer :: tracer_CSp !< The structure for controlling
                                                   !! calls to auxiliary ice tracer packages

  ! Local variables
  integer :: i, j, k, n, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  ! Set appropriate surface quantities in categories with no ice.
  if (allocated(IST%t_surf)) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)<=0.0) &
      IST%t_surf(i,j,k) = IST%T_0degC + OSS%T_fr_ocn(i,j)
    enddo ; enddo ; enddo
  endif

  ! Calculate and output various diagnostics of the ice state.
  call cpu_clock_begin(iceClock9)

  call enable_SIS_averaging(US%T_to_s*dt_slow, CS%Time, CS%diag)
  call post_ice_state_diagnostics(CS%IDs, IST, OSS, IOF, dt_slow, CS%Time, G, US, IG, CS%diag)
  call disable_SIS_averaging(CS%diag)

  if (CS%verbose) call ice_line(CS%Time, IST%part_size(:,:,0), OSS%SST_C(:,:), G)
  if (CS%debug) call IST_chksum("End ice_state_cleanup", IST, G, US, IG)
  if (CS%bounds_check) call IST_bounds_check(IST, G, US, IG, "End of ice_state_cleanup", OSS=OSS)

  if (CS%Time + real_to_time(0.5*US%T_to_s*dt_slow) > CS%write_ice_stats_time) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
                              tracer_CSp=tracer_CSp)
    CS%write_ice_stats_time = CS%write_ice_stats_time + CS%ice_stats_interval
  elseif (CS%column_check) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp)
  endif

  call cpu_clock_end(iceClock9)

end subroutine ice_state_cleanup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Convert the category-resolved ice state into the simplified 2-d ice state and a cell averaged state.
subroutine convert_IST_to_simple_state(IST, DS2d, CAS, G, US, IG, CS)
  type(ice_state_type),          intent(inout) :: IST !< A type describing the state of the sea ice
  type(dyn_state_2d),            intent(inout) :: DS2d !< A simplified 2-d description of the ice state
                                                      !! integrated across thickness categories and layers.
  type(cell_average_state_type), intent(inout) :: CAS !< A structure with ocean-cell averaged masses.
  type(SIS_hor_grid_type),       intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),         intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),           intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),            pointer       :: CS  !< The control structure for the SIS_dyn_trans module

  ! Local variables
  integer :: i, j, k, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  if (DS2d%nts /= 0) then
    call SIS_error(WARNING, "convert_IST_to_simple_state called with incomplete transport.")
    return
  endif

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call cpu_clock_begin(iceClock4)

  DS2d%mca_step(:,:,0) = 0.0 ; DS2d%mi_sum(:,:) = 0.0 ; DS2d%ice_cover(:,:) = 0.0
  !$OMP parallel do default(shared)
  do j=jsd,jed ; do k=1,IG%CatIce ; do i=isd,ied
    DS2d%mca_step(i,j,0) = DS2d%mca_step(i,j,0) + IST%part_size(i,j,k) * &
                    ((IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k)))
    DS2d%mi_sum(i,j) = DS2d%mi_sum(i,j) + (IST%mH_ice(i,j,k))  * IST%part_size(i,j,k)
    DS2d%ice_cover(i,j) = DS2d%ice_cover(i,j) + IST%part_size(i,j,k)
  enddo ; enddo ; enddo
  do j=jsd,jed ; do i=isd,ied
    DS2d%mca_step(i,j,0) = DS2d%mca_step(i,j,0) + DS2d%mi_sum(i,j)
!    if ((abs(max(1.0-DS2d%ice_cover(i,j),0.0) - IST%part_size(i,j,0)) > 5.0e-15) .and. (G%mask2dT(i,j)>0.0)) then
!      write(mesg, '(3(ES13.5))') max(1.0 - DS2d%ice_cover(i,j), 0.0) - IST%part_size(i,j,0), &
!         max(1.0 - DS2d%ice_cover(i,j), 0.0), IST%part_size(i,j,0)
!      call SIS_error(FATAL, "Mismatch in ice_free values exceeding roundoff: "//trim(mesg))
!    endif
  enddo ; enddo
  if (CS%Cgrid_dyn) then
    do j=jsd,jed ; do I=IsdB,IedB ; DS2d%u_ice_C(I,j) = IST%u_ice_C(I,j) ; enddo ; enddo
    do J=JsdB,JedB ; do i=isd,ied ; DS2d%v_ice_C(i,J) = IST%v_ice_C(i,J) ; enddo ; enddo
  else
    do J=JsdB,JedB ; do I=IsdB,IedB
      DS2d%u_ice_B(I,J) = IST%u_ice_B(I,J) ; DS2d%v_ice_B(I,J) = IST%v_ice_B(I,J)
    enddo ; enddo
  endif

  !  Determine the whole-cell averaged mass of snow and ice.
  call ice_state_to_cell_ave_state(IST, G, US, IG, CS%SIS_transport_CSp, CAS)

  IST%valid_IST = .false.

  call cpu_clock_end(iceClock4)

end subroutine convert_IST_to_simple_state


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Update the category-merged ice state and call the merged continuity update.
subroutine SIS_merged_dyn_cont(OSS, FIA, IOF, DS2d, IST, dt_cycle, Time_start, G, US, IG, CS, OBC, end_call)
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(in)    :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  type(dyn_state_2d),         intent(inout) :: DS2d !< A simplified 2-d description of the ice state
                                                   !! integrated across thickness categories and layers.
  type(ice_state_type),       intent(in)    :: IST !< A type describing the state of the sea ice.
  real,                       intent(in)    :: dt_cycle !< The slow ice dynamics timestep [T ~> s].
  type(time_type),            intent(in)    :: TIme_start !< The starting time for this update cycle.
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module
  type(ice_OBC_type),         pointer       :: OBC !< Open boundary structure.
  logical,          optional, intent(in)    :: end_call !< If present and false, this call is
                                                   !! the last in the series of advective updates.

  ! This subroutine updates the 2-d sea-ice dynamics.
  !   Variables updated here: DS2d%ice_cover, DS2d%[uv]_ice_[BC], DS2d%mca_step, DS2d%mi_sum,
  !       CS%[uv]h_step, DS2d%nts, CS%SIS_[BC]_dyn_CSp,  IOF (stresses)

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    ice_free, &         ! The fractional open water [nondim], between 0 & 1.
    rdg_rate            ! A ridging rate [T-1 ~> s-1], this will be calculated from the strain rates
                        ! in the dynamics.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    WindStr_x_B, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_B, &      ! averaged over the ice categories on a B-grid [R Z L T-2 ~> Pa].
    WindStr_x_ocn_B, &  ! Zonal wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    WindStr_y_ocn_B, &  ! Meridional wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    str_x_ice_ocn_B, &  ! Zonal ice-ocean stress on a B-grid [R Z L T-2 ~> Pa].
    str_y_ice_ocn_B     ! Meridional ice-ocean stress on a B-grid [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G))  :: &
    WindStr_x_Cu, &   ! Zonal wind stress averaged over the ice categories on C-grid u-points [R Z L T-2 ~> Pa].
    WindStr_x_ocn_Cu, & ! Zonal wind stress on the ice-free ocean on C-grid u-points [R Z L T-2 ~> Pa].
    str_x_ice_ocn_Cu   ! Zonal ice-ocean stress on C-grid u-points [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G))  :: &
    WindStr_y_Cv, &   ! Meridional wind stress averaged over the ice categories on C-grid v-points [R Z L T-2 ~> Pa].
    WindStr_y_ocn_Cv, & ! Meridional wind stress on the ice-free ocean on C-grid v-points [R Z L T-2 ~> Pa].
    str_y_ice_ocn_Cv  ! Meridional ice-ocean stress on C-grid v-points [R Z L T-2 ~> Pa].

  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBx ! A temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBy ! A temporary array for diagnostics.

  real :: ps_vel   ! The fractional thickness category coverage at a velocity point.
  real :: wt_new, wt_prev ! Weights in an average.
  real :: dt_slow_dyn  ! The slow dynamics timestep [T ~> s].
  real :: dt_slow_dyn_sec ! The slow dynamics timestep [s].
  real :: dt_adv       ! The advective subcycle timestep [T ~> s].
  logical :: continuing_call ! If true, there are more in the series of advective updates
                             ! after this call.
  integer :: ndyn_steps, nds ! The number of dynamic steps in this call.
  integer :: i, j, k, n, isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed !, IsdB, IedB, JsdB, JedB

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ndyn_steps = 1
  if ((CS%dt_ice_dyn > 0.0) .and. (CS%dt_ice_dyn < dt_cycle)) &
    ndyn_steps = max(CEILING(dt_cycle/CS%dt_ice_dyn - 1e-6), 1)
  dt_slow_dyn = dt_cycle / ndyn_steps
  dt_slow_dyn_sec = US%T_to_s*dt_slow_dyn
  dt_adv = dt_slow_dyn / real(CS%adv_substeps)
  if (ndyn_steps*CS%adv_substeps > DS2d%max_nts) &
    call increase_max_tracer_step_memory(DS2d, G, ndyn_steps*CS%adv_substeps)
  continuing_call = .false. ; if (present(end_call)) continuing_call = .not.end_call

  do nds=1,ndyn_steps
    call cpu_clock_begin(iceClock4)
    call enable_SIS_averaging(dt_slow_dyn_sec, Time_start + real_to_time(nds*dt_slow_dyn_sec), CS%diag)
    do j=jsd,jed ; do i=isd,ied ; ice_free(i,j) = max(1.0 - DS2d%ice_cover(i,j), 0.0) ; enddo ; enddo

    ! In the dynamics code, only the ice velocities are changed, and the ice-ocean
    ! stresses are calculated.  The gravity wave dynamics (i.e. the continuity
    ! equation) are not included in the dynamics (yet).  All of the thickness categories
    ! are merged together.
    if (CS%Cgrid_dyn) then

      ! Correct the wind stresses for changes in the fractional ice-coverage and set
      ! the wind stresses on the ice and the open ocean for a C-grid staggering.
      ! This block of code must be executed if ice_cover and ice_free or the various wind
      ! stresses were updated.
      call set_wind_stresses_C(FIA, DS2d%ice_cover, ice_free, WindStr_x_Cu, WindStr_y_Cv, &
                               WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, G, US, CS%complete_ice_cover, OBC)

      if (CS%lemieux_landfast) then
        call basal_stress_coeff_C(G, DS2d%mi_sum, DS2d%ice_cover, OSS%sea_lev, CS%SIS_C_dyn_CSp)
      elseif (CS%itd_landfast) then
        call basal_stress_coeff_itd(G, IG, IST, OSS%sea_lev, CS%SIS_C_dyn_CSp)
      endif

      if (CS%debug) then
        call uvchksum("Before SIS_C_dynamics [uv]_ice_C", DS2d%u_ice_C, DS2d%v_ice_C, G, scale=US%L_T_to_m_s)
        call hchksum(ice_free, "ice_free before SIS_C_dynamics", G%HI)
        call hchksum(DS2d%mca_step(:,:,DS2d%nts), "misp_sum before SIS_C_dynamics", G%HI, scale=US%RZ_to_kg_m2)
        call hchksum(DS2d%mi_sum, "mi_sum before SIS_C_dynamics", G%HI, scale=US%RZ_to_kg_m2)
        call hchksum(OSS%sea_lev, "sea_lev before SIS_C_dynamics", G%HI, haloshift=1, scale=US%Z_to_m)
        call hchksum(DS2d%ice_cover, "ice_cover before SIS_C_dynamics", G%HI, haloshift=1)
        call uvchksum("[uv]_ocn before SIS_C_dynamics", OSS%u_ocn_C, OSS%v_ocn_C, G, halos=1, scale=US%L_T_to_m_s)
        call uvchksum("WindStr_[xy] before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G, &
                      halos=1, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
!        call hchksum_pair("WindStr_[xy]_A before SIS_C_dynamics", WindStr_x_A, WindStr_y_A, G, halos=1)
      endif

      call cpu_clock_begin(iceClocka)
      !### Ridging needs to be added with C-grid dynamics.
      if (CS%do_ridging) rdg_rate(:,:) = 0.0
      call SIS_C_dynamics(DS2d%ice_cover, DS2d%mca_step(:,:,DS2d%nts), DS2d%mi_sum, &
                          DS2d%u_ice_C, DS2d%v_ice_C, &
                          OSS%u_ocn_C, OSS%v_ocn_C, WindStr_x_Cu, WindStr_y_Cv, OSS%sea_lev, &
                          str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, dt_slow_dyn, G, US, CS%SIS_C_dyn_CSp)

      call cpu_clock_end(iceClocka)

      if (CS%debug) call uvchksum("After ice_dynamics [uv]_ice_C", DS2d%u_ice_C, DS2d%v_ice_C, G, scale=US%L_T_to_m_s)

      call cpu_clock_begin(iceClockb)
      call pass_vector(DS2d%u_ice_C, DS2d%v_ice_C, G%Domain, stagger=CGRID_NE)
      call pass_vector(str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, G%Domain, stagger=CGRID_NE)
      call cpu_clock_end(iceClockb)

      ! Dynamics diagnostics
      call cpu_clock_begin(iceClockc)
      if (CS%id_fax>0) call post_data(CS%id_fax, WindStr_x_Cu, CS%diag)
      if (CS%id_fay>0) call post_data(CS%id_fay, WindStr_y_Cv, CS%diag)
      if (CS%debug) call uvchksum("Before set_ocean_top_stress_Cgrid [uv]_ice_C", &
                                  DS2d%u_ice_C, DS2d%v_ice_C, G, scale=US%L_T_to_m_s)

      ! Store all mechanical ocean forcing.
      call set_ocean_top_stress_C2(IOF, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, &
                                   str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, ice_free, DS2d%ice_cover, G, US, OBC)
      call cpu_clock_end(iceClockc)

    else ! B-grid dynamics.

      ! Correct the wind stresses for changes in the fractional ice-coverage and set
      ! the wind stresses on the ice and the open ocean for a C-grid staggering.
      ! This block of code must be executed if ice_cover and ice_free or the various wind
      ! stresses were updated.

      call set_wind_stresses_B(FIA, DS2d%ice_cover, ice_free, WindStr_x_B, WindStr_y_B, &
                               WindStr_x_ocn_B, WindStr_y_ocn_B, G, US, CS%complete_ice_cover)

      if (CS%debug) then
        call Bchksum_pair("[uv]_ice_B before dynamics", DS2d%u_ice_B, DS2d%v_ice_B, G, scale=US%L_T_to_m_s)
        call hchksum(ice_free, "ice_free before ice_dynamics", G%HI)
        call hchksum(DS2d%mca_step(:,:,DS2d%nts), "misp_sum before ice_dynamics", G%HI, scale=US%RZ_to_kg_m2)
        call hchksum(DS2d%mi_sum, "mi_sum before ice_dynamics", G%HI, scale=US%RZ_to_kg_m2)
        call hchksum(OSS%sea_lev, "sea_lev before ice_dynamics", G%HI, haloshift=1, scale=US%Z_to_m)
        call Bchksum_pair("[uv]_ocn before ice_dynamics", OSS%u_ocn_B, OSS%v_ocn_B, G, scale=US%L_T_to_m_s)
        call Bchksum_pair("WindStr_[xy]_B before ice_dynamics", WindStr_x_B, WindStr_y_B, G, halos=1, &
                          scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      endif

      call cpu_clock_begin(iceClocka)
      if (CS%do_ridging) rdg_rate(:,:) = 0.0
      call SIS_B_dynamics(DS2d%ice_cover, DS2d%mca_step(:,:,DS2d%nts), DS2d%mi_sum, &
                          DS2d%u_ice_B, DS2d%v_ice_B, &
                          OSS%u_ocn_B, OSS%v_ocn_B, WindStr_x_B, WindStr_y_B, OSS%sea_lev, &
                          str_x_ice_ocn_B, str_y_ice_ocn_B, CS%do_ridging, &
                          rdg_rate, dt_slow_dyn, G, US, CS%SIS_B_dyn_CSp)
      call cpu_clock_end(iceClocka)

      if (CS%debug) call Bchksum_pair("After dynamics [uv]_ice_B", DS2d%u_ice_B, DS2d%v_ice_B, G, scale=US%L_T_to_m_s)

      call cpu_clock_begin(iceClockb)
      call pass_vector(DS2d%u_ice_B, DS2d%v_ice_B, G%Domain, stagger=BGRID_NE)
      call cpu_clock_end(iceClockb)

      ! Dynamics diagnostics
      call cpu_clock_begin(iceClockc)
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

      if (CS%debug) call Bchksum_pair("Before set_ocean_top_stress_Bgrid [uv]_ice_B", DS2d%u_ice_B, DS2d%v_ice_B, &
                                      G, scale=US%L_T_to_m_s)
      ! Store all mechanical ocean forcing.
      call set_ocean_top_stress_B2(IOF, WindStr_x_ocn_B, WindStr_y_ocn_B, &
                                   str_x_ice_ocn_B, str_y_ice_ocn_B, ice_free, DS2d%ice_cover, G, US)
      call cpu_clock_end(iceClockc)

      ! Convert the velocities to C-grid points for use in transport.
      do j=jsc,jec ; do I=isc-1,iec
        DS2d%u_ice_C(I,j) = 0.5 * ( DS2d%u_ice_B(I,J-1) + DS2d%u_ice_B(I,J) )
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        DS2d%v_ice_C(i,J) = 0.5 * ( DS2d%v_ice_B(I-1,J) + DS2d%v_ice_B(I,J) )
      enddo ; enddo
    endif ! End of B-grid dynamics

    if (CS%do_ridging) then ! Accumulate the time-average ridging rate.
      DS2d%ridge_rate_count = DS2d%ridge_rate_count + 1.
      wt_new = 1.0 / DS2d%ridge_rate_count ; wt_prev = 1.0 - wt_new
      do j=jsc,jec ; do i=isc,iec
        DS2d%avg_ridge_rate(i,j) = wt_new * rdg_rate(i,j) + wt_prev * DS2d%avg_ridge_rate(i,j)
      enddo ; enddo
    endif

    if (CS%debug) call uvchksum("Before ice_transport [uv]_ice_C", DS2d%u_ice_C, DS2d%v_ice_C, G, scale=US%L_T_to_m_s)
    call enable_SIS_averaging(dt_slow_dyn_sec, Time_start + real_to_time(nds*dt_slow_dyn_sec), CS%diag)

    ! Update the integrated ice mass and store the transports in each step.
    if (DS2d%nts+CS%adv_substeps > DS2d%max_nts) &
      call increase_max_tracer_step_memory(DS2d, G, DS2d%nts+CS%adv_substeps)

    do n = DS2d%nts+1, DS2d%nts+CS%adv_substeps
      if ((n < ndyn_steps*CS%adv_substeps) .or. continuing_call) then
        ! Some of the work is not needed for the last step before cat_ice_transport.
        call summed_continuity(DS2d%u_ice_C, DS2d%v_ice_C, DS2d%mca_step(:,:,n-1), DS2d%mca_step(:,:,n), &
                               DS2d%uh_step(:,:,n), DS2d%vh_step(:,:,n), dt_adv, G, US, IG, &
                               CS%continuity_CSp, h_ice=DS2d%mi_sum)
        call ice_cover_transport(DS2d%u_ice_C, DS2d%v_ice_C, DS2d%ice_cover, dt_adv, G, US, IG, CS%cover_trans_CSp, &
                                 masking_uhtot=DS2d%uh_step(:,:,n), masking_vhtot=DS2d%vh_step(:,:,n))
        call pass_var(DS2d%mi_sum, G%Domain, complete=.false.)
        call pass_var(DS2d%ice_cover, G%Domain, complete=.false.)
        call pass_var(DS2d%mca_step(:,:,n), G%Domain, complete=.true.)
      else
        call summed_continuity(DS2d%u_ice_C, DS2d%v_ice_C, DS2d%mca_step(:,:,n-1), DS2d%mca_step(:,:,n), &
                               DS2d%uh_step(:,:,n), DS2d%vh_step(:,:,n), dt_adv, G, US, IG, CS%continuity_CSp)
      endif
    enddo
    DS2d%nts = DS2d%nts + CS%adv_substeps
    call cpu_clock_end(iceClock4)

  enddo ! nds=1,ndyn_steps

end subroutine SIS_merged_dyn_cont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> slab_ice_dynamics_trans makes the calls to do the slab ice version of dynamics and mass and tracer transport
subroutine slab_ice_dyn_trans(IST, OSS, FIA, IOF, dt_slow, CS, G, US, IG, tracer_CSp, OBC)

  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(inout) :: FIA !< A type containing averages of fields
                                                   !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow !< The slow ice dynamics timestep [T ~> s].
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(dyn_trans_CS),         pointer       :: CS  !< The control structure for the SIS_dyn_trans module
  type(SIS_tracer_flow_control_CS), pointer :: tracer_CSp !< The structure for controlling calls to
                                                   !! auxiliary ice tracer packages
  type(ice_OBC_type),         pointer       :: OBC !< Open boundary structure.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    mi_sum, &           ! Masses of ice per unit total area [R Z ~> kg m-2].
    misp_sum            ! Combined mass of snow, ice and melt pond water per unit total area [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    WindStr_x_B, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_B, &      ! averaged over the ice categories on a B-grid [R Z L T-2 ~> Pa].
    WindStr_x_ocn_B, &  ! Zonal wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    WindStr_y_ocn_B, &  ! Meridional wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    str_x_ice_ocn_B, &  ! Zonal ice-ocean stress on a B-grid [R Z L T-2 ~> Pa].
    str_y_ice_ocn_B     ! Meridional ice-ocean stress on a B-grid [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G))  :: &
    WindStr_x_Cu, &   ! Zonal wind stress averaged over the ice categories on C-grid u-points [R Z L T-2 ~> Pa].
    WindStr_x_ocn_Cu, & ! Zonal wind stress on the ice-free ocean on C-grid u-points [R Z L T-2 ~> Pa].
    str_x_ice_ocn_Cu   ! Zonal ice-ocean stress on C-grid u-points [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G))  :: &
    WindStr_y_Cv, &   ! Meridional wind stress averaged over the ice categories on C-grid v-points [R Z L T-2 ~> Pa].
    WindStr_y_ocn_Cv, & ! Meridional wind stress on the ice-free ocean on C-grid v-points [R Z L T-2 ~> Pa].
    str_y_ice_ocn_Cv  ! Meridional ice-ocean stress on C-grid v-points [R Z L T-2 ~> Pa].

  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBx ! A temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBy ! A temporary array for diagnostics.
  real :: ps_vel   ! The fractional thickness category coverage at a velocity point.
  real :: dt_slow_dyn  ! The slow dynamics timestep [T ~> s].
  real :: dt_slow_dyn_sec ! The slow dynamics timestep [s].
  integer :: i, j, k, n, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed
  integer :: ndyn_steps, nds ! The number of dynamic steps.

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = 1
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%n_calls = CS%n_calls + 1
  IOF%stress_count = 0

  ndyn_steps = 1
  if ((CS%dt_ice_dyn > 0.0) .and. (CS%dt_ice_dyn < dt_slow)) &
    ndyn_steps = max(CEILING(dt_slow/CS%dt_ice_dyn - 0.000001), 1)
  dt_slow_dyn = dt_slow / ndyn_steps
  dt_slow_dyn_sec = US%T_to_s*dt_slow_dyn

  do nds=1,ndyn_steps

    call enable_SIS_averaging(dt_slow_dyn_sec, CS%Time - real_to_time((ndyn_steps-nds)*dt_slow_dyn_sec), CS%diag)

    call cpu_clock_begin(iceClock4)
    !$OMP parallel do default(shared)
    do j=jsd,jed ; do i=isd,ied
      mi_sum(i,j) = IST%mH_ice(i,j,1) * IST%part_size(i,j,1)
      misp_sum(i,j) = mi_sum(i,j) + IST%part_size(i,j,1) * &
                      (IST%mH_snow(i,j,1) + IST%mH_pond(i,j,1))
    enddo ; enddo
    call cpu_clock_end(iceClock4)

    !
    ! Dynamics - update ice velocities.
    !

    ! In the dynamics code, only the ice velocities are changed, and the ice-ocean
    ! stresses are calculated.  The gravity wave dynamics (i.e. the continuity
    ! equation) are not included in the dynamics.  All of the thickness categories
    ! are merged together.
    if (CS%Cgrid_dyn) then

      call cpu_clock_begin(iceClock4)
      ! Correct the wind stresses for changes in the fractional ice-coverage and set
      ! the wind stresses on the ice and the open ocean for a C-grid staggering.
      ! This block of code must be executed if ice_cover and ice_free or the various wind
      ! stresses were updated.
      call set_wind_stresses_C(FIA, IST%part_size(:,:,1), IST%part_size(:,:,0), WindStr_x_Cu, WindStr_y_Cv, &
                               WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, G, US, CS%complete_ice_cover, OBC)

      if (CS%debug) then
        call uvchksum("Before SIS_C_dynamics [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G, scale=US%L_T_to_m_s)
        call hchksum(IST%part_size(:,:,0), "ice_free before SIS_C_dynamics", G%HI)
        call hchksum(misp_sum, "misp_sum before SIS_C_dynamics", G%HI, scale=US%RZ_to_kg_m2)
        call hchksum(mi_sum, "mi_sum before SIS_C_dynamics", G%HI, scale=US%RZ_to_kg_m2)
        call hchksum(OSS%sea_lev, "sea_lev before SIS_C_dynamics", G%HI, haloshift=1, scale=US%Z_to_m)
        call hchksum(IST%part_size(:,:,1), "ice_cover before SIS_C_dynamics", G%HI, haloshift=1)
        call uvchksum("[uv]_ocn before SIS_C_dynamics", OSS%u_ocn_C, OSS%v_ocn_C, G, halos=1, scale=US%L_T_to_m_s)
        call uvchksum("WindStr_[xy] before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G, &
                      halos=1, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
!        call hchksum_pair("WindStr_[xy]_A before SIS_C_dynamics", WindStr_x_A, WindStr_y_A, G, halos=1)
      endif

      call cpu_clock_begin(iceClocka)
      call slab_ice_dynamics(IST%u_ice_C, IST%v_ice_C, OSS%u_ocn_C, OSS%v_ocn_C, &
                             WindStr_x_Cu, WindStr_y_Cv, str_x_ice_ocn_Cu, str_y_ice_ocn_Cv)
      call cpu_clock_end(iceClocka)

      if (CS%debug) call uvchksum("After ice_dynamics [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G, scale=US%L_T_to_m_s)

      call cpu_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_C, IST%v_ice_C, G%Domain, stagger=CGRID_NE)
      call pass_vector(str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, G%Domain, stagger=CGRID_NE)
      call cpu_clock_end(iceClockb)

      ! Dynamics diagnostics
      call cpu_clock_begin(iceClockc)
      if (CS%id_fax>0) call post_data(CS%id_fax, WindStr_x_Cu, CS%diag)
      if (CS%id_fay>0) call post_data(CS%id_fay, WindStr_y_Cv, CS%diag)

      if (CS%debug) call uvchksum("Before set_ocean_top_stress_Cgrid [uv]_ice_C", &
                                  IST%u_ice_C, IST%v_ice_C, G, scale=US%L_T_to_m_s)

      ! Store all mechanical ocean forcing.
      call set_ocean_top_stress_C2(IOF, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, &
                                   IST%part_size(:,:,0), IST%part_size(:,:,1), G, US, OBC)
      call cpu_clock_end(iceClockc)

      call cpu_clock_end(iceClock4)

    else ! B-grid dynamics.

      call cpu_clock_begin(iceClock4)
      ! Correct the wind stresses for changes in the fractional ice-coverage and set
      ! the wind stresses on the ice and the open ocean for a C-grid staggering.
      ! This block of code must be executed if ice_cover and ice_free or the various wind
      ! stresses were updated.

      call set_wind_stresses_B(FIA, IST%part_size(:,:,1), IST%part_size(:,:,0), WindStr_x_B, WindStr_y_B, &
                               WindStr_x_ocn_B, WindStr_y_ocn_B, G, US, CS%complete_ice_cover)

      if (CS%debug) then
        call Bchksum_pair("[uv]_ice_B before dynamics", IST%u_ice_B, IST%v_ice_B, G, scale=US%L_T_to_m_s)
        call hchksum(IST%part_size(:,:,0), "ice_free before ice_dynamics", G%HI)
        call hchksum(misp_sum, "misp_sum before ice_dynamics", G%HI, scale=US%RZ_to_kg_m2)
        call hchksum(mi_sum, "mi_sum before ice_dynamics", G%HI, scale=US%RZ_to_kg_m2)
        call hchksum(OSS%sea_lev, "sea_lev before ice_dynamics", G%HI, haloshift=1, scale=US%Z_to_m)
        call Bchksum_pair("[uv]_ocn before ice_dynamics", OSS%u_ocn_B, OSS%v_ocn_B, G, scale=US%L_T_to_m_s)
        call Bchksum_pair("WindStr_[xy]_B before ice_dynamics", WindStr_x_B, WindStr_y_B, G, halos=1, &
                          scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      endif

      call cpu_clock_begin(iceClocka)
      call slab_ice_dynamics(IST%u_ice_B, IST%v_ice_B, OSS%u_ocn_B, OSS%v_ocn_B, &
                             WindStr_x_B, WindStr_y_B, str_x_ice_ocn_B, str_y_ice_ocn_B)
      call cpu_clock_end(iceClocka)

      if (CS%debug) call Bchksum_pair("After dynamics [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G, scale=US%L_T_to_m_s)

      call cpu_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
      call cpu_clock_end(iceClockb)

      ! Dynamics diagnostics
      call cpu_clock_begin(iceClockc)
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

      if (CS%debug) call Bchksum_pair("Before set_ocean_top_stress_Bgrid [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G, &
                                      scale=US%L_T_to_m_s)
      ! Store all mechanical ocean forcing.
      call set_ocean_top_stress_B2(IOF, WindStr_x_ocn_B, WindStr_y_ocn_B, str_x_ice_ocn_B, str_y_ice_ocn_B, &
                                   IST%part_size(:,:,0), IST%part_size(:,:,1), G, US)
      call cpu_clock_end(iceClockc)

       ! Convert the B-grid velocities to C-grid points for transport.
      if (CS%debug) call Bchksum_pair("Before ice_transport [uv]_ice_B", IST%u_ice_B, IST%v_ice_B, G, &
                                      scale=US%L_T_to_m_s)
      do j=jsc,jec ; do I=isc-1,iec
        IST%u_ice_C(I,j) = 0.5 * ( IST%u_ice_B(I,J-1) + IST%u_ice_B(I,J) )
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        IST%v_ice_C(i,J) = 0.5 * ( IST%v_ice_B(I-1,J) + IST%v_ice_B(I,J) )
      enddo ; enddo

      call cpu_clock_end(iceClock4)

    endif ! End of B-grid dynamics

    ! Do ice mass transport and related tracer transport.  This updates the category-decomposed ice state.
    call cpu_clock_begin(iceClock8)
    if (CS%debug) call uvchksum("Before ice_transport [uv]_ice_C", IST%u_ice_C, IST%v_ice_C, G, scale=US%L_T_to_m_s)
    call enable_SIS_averaging(dt_slow_dyn_sec, CS%Time - real_to_time((ndyn_steps-nds)*dt_slow_dyn_sec), CS%diag)

    call slab_ice_advect(IST%u_ice_C, IST%v_ice_C, IST%mH_ice(:,:,1), 4.0*US%kg_m3_to_R*US%m_to_Z, &
                         dt_slow_dyn, G, US, IST%part_size(:,:,1), nsteps=CS%adv_substeps)
    call cpu_clock_end(iceClock8)

    if (CS%column_check) &
      call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
                                message="      Post_transport")! , check_column=.true.)

  enddo ! nds=1,ndyn_steps
  call finish_ocean_top_stresses(IOF, G)

  call ice_state_cleanup(IST, OSS, IOF, dt_slow, G, US, IG, CS, tracer_CSp)

end subroutine slab_ice_dyn_trans


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Finish setting the ice-ocean stresses by dividing the running sums of the
!! stresses by the number of times they have been augmented.  It may also record
!! the current ocean-cell averaged ice, snow and pond mass.
subroutine finish_ocean_top_stresses(IOF, G, DS2d, IG)
  type(ice_ocean_flux_type),    intent(inout) :: IOF  !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),      intent(in)    :: G    !< The horizontal grid type
  type(dyn_state_2d), optional, intent(in)    :: DS2d !< A simplified 2-d description of the ice state
                                                  !! integrated across thickness categories and layers.
  type(ice_grid_type), optional, intent(in)   :: IG  !< The sea-ice specific grid type

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

  if (present(DS2d) .and. present(IG)) then ; if (DS2d%nts > 0) then ; do j=jsc,jec ; do i=isc,iec
    IOF%mass_ice_sn_p(i,j) = DS2d%mca_step(i, j, DS2d%nts)
  enddo ; enddo ; endif ; endif

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
!! be set at all (symmetric) edge points unless stagger is AGRID.  The units of stress_mag
!! are the same as str_x and str_y.
subroutine stresses_to_stress_mag(G, str_x, str_y, stagger, stress_mag)
  type(SIS_hor_grid_type),   intent(in)    :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: str_x !< The x-direction ice to ocean stress [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: str_y !< The y-direction ice to ocean stress [R Z L T-2 ~> Pa].
  integer,                   intent(in)    :: stagger !< The staggering relative to the tracer points of the
                                                  !! two wind stress components. Valid entries include AGRID,
                                                  !! BGRID_NE, and CGRID_NE, following the Arakawa
                                                  !! grid-staggering  notation.  BGRID_SW and CGRID_SW are
                                                  !! possibilities that have not been implemented yet.
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(inout) :: stress_mag !< The magnitude of the stress at tracer points
                                                  !! in the same units as str_x and str_y [R Z L T-2 ~> Pa].

  ! Local variables
  real :: taux2, tauy2  ! squared wind stress components [R2 Z2 L2 T-4 ~> Pa2]
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
        taux2 = (G%mask2dCu(I-1,j)*str_x(I-1,j)**2 + G%mask2dCu(I,j)*str_x(I,j)**2) / &
                (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0) &
        tauy2 = (G%mask2dCv(i,J-1)*str_y(i,J-1)**2 + G%mask2dCv(i,J)*str_y(i,J)**2) / &
                (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))
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
                                      str_ice_oce_x, str_ice_oce_y, part_size, G, US, IG)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),       intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over
                                                  !! open water [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over
                                                  !! open water [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_x   !< The x-direction ice to ocean
                                                  !! stress [R Z L T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y   !< The y-direction ice to ocean
                                                  !! stress [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G),0:IG%CatIce), &
                             intent(in)    :: part_size !< The fractional area coverage of the ice
                                                  !! thickness categories [nondim], 0-1
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors

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
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then
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
        ps_vel = 1.0
        if (G%mask2dBu(I,J)>0.0) ps_vel = 0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                                                (part_size(i+1,j,0) + part_size(i,j+1,0)) )
        IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + windstr_x_water(I,J) * ps_vel
        IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + windstr_y_water(I,J) * ps_vel
      enddo
      do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dBu(I,J)>0.0) then
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
        ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.0) ps_vel = 0.5*(part_size(i+1,j,0) + part_size(i,j,0))
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * &
                0.5 * (windstr_x_water(I,J) + windstr_x_water(I,J-1))
      enddo
      do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dCu(I,j)>0.0) then
        ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * &
                0.5 * (str_ice_oce_x(I,J) + str_ice_oce_x(I,J-1))
      endif ; enddo ; enddo
    enddo
    !$OMP parallel do default(shared) private(ps_vel)
    do J=jsc-1,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.0) ps_vel = 0.5*(part_size(i,j+1,0) + part_size(i,j,0))
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * &
                0.5 * (windstr_y_water(I,J) + windstr_y_water(I-1,J))
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dCv(i,J)>0.0) then
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
                                      str_ice_oce_x, str_ice_oce_y, part_size, G, US, IG, OBC)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),       intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over open water
                                                  !! [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over open water
                                                  !! [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: str_ice_oce_x !< The x-direction ice to ocean stress [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y !< The y-direction ice to ocean stress [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G),0:IG%CatIce), &
                             intent(in)    :: part_size !< The fractional area coverage of the ice
                                                  !! thickness categories [nondim], 0-1
  type(ice_OBC_type),        pointer       :: OBC  !< Open boundary structure.

  real    :: ps_vel ! part_size interpolated to a velocity point [nondim].
  integer :: i, j, k, isc, iec, jsc, jec, ncat
  integer :: l_seg
  logical :: local_open_u_BC, local_open_v_BC
  type(OBC_segment_type), pointer :: segment => NULL()

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  local_open_u_BC = .false. ; local_open_v_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_open_u_BC = OBC%open_u_BCs_exist_globally
    local_open_v_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  if ((local_open_u_BC .or. local_open_v_BC) .and. &
      (IOF%flux_uv_stagger == AGRID) .or. (IOF%flux_uv_stagger == BGRID_NE)) &
        call SIS_error(FATAL, "No open boundaries for given flux staggering")

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
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.0) then
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
        ps_vel = 1.0
        if (G%mask2dBu(I,J)>0.0) ps_vel = 0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                                                (part_size(i+1,j,0) + part_size(i,j+1,0)) )
        IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
                (windstr_x_water(I,j) + windstr_x_water(I,j+1))
        IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
                (windstr_y_water(i,J) + windstr_y_water(i+1,J))
      enddo
      do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dBu(I,J)>0.0) then
        ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                         (part_size(i+1,j,k) + part_size(i,j+1,k)) )
        IOF%flux_u_ocn(I,J) = IOF%flux_u_ocn(I,J) + ps_vel * 0.5 * &
                            (str_ice_oce_x(I,j) + str_ice_oce_x(I,j+1))
        IOF%flux_v_ocn(I,J) = IOF%flux_v_ocn(I,J) + ps_vel * 0.5 * &
                            (str_ice_oce_y(i,J) + str_ice_oce_y(i+1,J))
      endif ; enddo ; enddo
    enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
    if (local_open_u_BC) then
      !$OMP parallel do default(shared) private(ps_vel)
      do j=jsc,jec
        do I=Isc-1,iec
          l_seg = OBC%segnum_u(I,j)
          ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.0) ps_vel = 0.5*(part_size(i+1,j,0) + part_size(i,j,0))
          if (l_seg /= OBC_NONE) then
            if (OBC%segment(l_seg)%open) then
              if (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
                ps_vel = part_size(i,j,0)
              else
                ps_vel = part_size(i+1,j,0)
              endif
          endif ; endif
          IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * windstr_x_water(I,j)
        enddo
        do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dCu(I,j)>0.0) then
          l_seg = OBC%segnum_u(I,j)
          ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
          if (l_seg /= OBC_NONE) then
            if (OBC%segment(l_seg)%open) then
              if (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
                ps_vel = part_size(i,j,k)
              else
                ps_vel = part_size(i+1,j,k)
              endif
          endif ; endif
          IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * str_ice_oce_x(I,j)
        endif ; enddo ; enddo
      enddo
    else
      !$OMP parallel do default(shared) private(ps_vel)
      do j=jsc,jec
        do I=Isc-1,iec
          ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.0) ps_vel = 0.5*(part_size(i+1,j,0) + part_size(i,j,0))
          IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * windstr_x_water(I,j)
        enddo
        do k=1,ncat ; do I=isc-1,iec ; if (G%mask2dCu(I,j)>0.0) then
          ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
          IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + ps_vel * str_ice_oce_x(I,j)
        endif ; enddo ; enddo
      enddo
    endif
    if (local_open_v_BC) then
      !$OMP parallel do default(shared) private(ps_vel)
      do J=jsc-1,jec
        do i=isc,iec
          l_seg = OBC%segnum_v(i,J)
          ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.0) ps_vel = 0.5*(part_size(i,j+1,0) + part_size(i,j,0))
          if (l_seg /= OBC_NONE) then
            if (OBC%segment(l_seg)%open) then
              if (OBC%segment(l_seg)%direction == OBC_DIRECTION_N) then
                ps_vel = part_size(i,j,0)
              else
                ps_vel = part_size(i,j+1,0)
              endif
          endif ; endif
          IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * windstr_y_water(i,J)
        enddo
        do k=1,ncat ; do i=isc,iec ; if (G%mask2dCv(i,J)>0.0) then
          l_seg = OBC%segnum_v(i,J)
          ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
          if (l_seg /= OBC_NONE) then
            if (OBC%segment(l_seg)%open) then
              if (OBC%segment(l_seg)%direction == OBC_DIRECTION_N) then
                ps_vel = part_size(i,j,k)
              else
                ps_vel = part_size(i,j+1,k)
              endif
          endif ; endif
          IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * str_ice_oce_y(i,J)
        endif ; enddo ; enddo
      enddo
    else
      !$OMP parallel do default(shared) private(ps_vel)
      do J=jsc-1,jec
        do i=isc,iec
          ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.0) ps_vel = 0.5*(part_size(i,j+1,0) + part_size(i,j,0))
          IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * windstr_y_water(i,J)
        enddo
        do k=1,ncat ; do i=isc,iec ; if (G%mask2dCv(i,J)>0.0) then
          ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
          IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + ps_vel * str_ice_oce_y(i,J)
        endif ; enddo ; enddo
      enddo
    endif
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
                                      str_ice_oce_x, str_ice_oce_y, ice_free, ice_cover, G, Us)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over
                                                  !! open water [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over
                                                  !! open water [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_x   !< The x-direction ice to ocean
                                                  !! stress [R Z L T-2 ~> Pa]
  real, dimension(SZIB_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y   !< The y-direction ice to ocean
                                                  !! stress [R Z L T-2 ~> Pa]
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
      if (G%mask2dBu(I,J)>0.0) then
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
      if (G%mask2dCu(I,j)>0.0) then
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
      if (G%mask2dCv(i,J)>0.0) then
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
                                   str_ice_oce_x, str_ice_oce_y, ice_free, ice_cover, G, US, OBC)
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: windstr_x_water !< The x-direction wind stress over
                                                  !! open water [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: windstr_y_water !< The y-direction wind stress over
                                                  !! open water [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G)), &
                             intent(in)    :: str_ice_oce_x !< The x-direction ice to ocean stress [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), &
                             intent(in)    :: str_ice_oce_y !< The y-direction ice to ocean stress [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: ice_free  !< The fractional open water area coverage [nondim], 0-1
  real, dimension(SZI_(G),SZJ_(G)), &
                             intent(in)    :: ice_cover !< The fractional ice area coverage [nondim], 0-1
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_OBC_type),        pointer       :: OBC  !< Open boundary structure.

  real    :: ps_ice, ps_ocn ! ice_free and ice_cover interpolated to a velocity point [nondim].
  integer :: i, j, k, isc, iec, jsc, jec
  integer :: l_seg
  logical :: local_open_u_BC, local_open_v_BC
  type(OBC_segment_type), pointer :: segment => NULL()

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  local_open_u_BC = .false. ; local_open_v_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_open_u_BC = OBC%open_u_BCs_exist_globally
    local_open_v_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  if ((local_open_u_BC .or. local_open_v_BC) .and. &
      (IOF%flux_uv_stagger == AGRID) .or. (IOF%flux_uv_stagger == BGRID_NE)) &
        call SIS_error(FATAL, "No open boundaries for given flux staggering")

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
      if (G%mask2dBu(I,J)>0.0) then
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
    if (local_open_u_BC) then
      !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
      do j=jsc,jec ; do I=Isc-1,iec
        ps_ocn = 1.0 ; ps_ice = 0.0
        l_seg = OBC%segnum_u(I,j)
        if (G%mask2dCu(I,j)>0.0) then
          ps_ocn = 0.5*(ice_free(i+1,j) + ice_free(i,j))
          ps_ice = 0.5*(ice_cover(i+1,j) + ice_cover(i,j))
        endif
        if (l_seg /= OBC_NONE) then
          if (OBC%segment(l_seg)%open) then
            if (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
              ps_ocn = ice_free(i,j)
              ps_ice = ice_cover(i,j)
            else
              ps_ocn = ice_free(i+1,j)
              ps_ice = ice_cover(i+1,j)
            endif
        endif ; endif
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + &
            (ps_ocn * windstr_x_water(I,j) + ps_ice * str_ice_oce_x(I,j))
      enddo ; enddo
    else
      !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
      do j=jsc,jec ; do I=Isc-1,iec
        ps_ocn = 1.0 ; ps_ice = 0.0
        if (G%mask2dCu(I,j)>0.0) then
          ps_ocn = 0.5*(ice_free(i+1,j) + ice_free(i,j))
          ps_ice = 0.5*(ice_cover(i+1,j) + ice_cover(i,j))
        endif
        IOF%flux_u_ocn(I,j) = IOF%flux_u_ocn(I,j) + &
            (ps_ocn * windstr_x_water(I,j) + ps_ice * str_ice_oce_x(I,j))
      enddo ; enddo
    endif
    if (local_open_v_BC) then
      !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
      do J=jsc-1,jec ; do i=isc,iec
        l_seg = OBC%segnum_v(i,J)
        ps_ocn = 1.0 ; ps_ice = 0.0
        if (G%mask2dCv(i,J)>0.0) then
          ps_ocn = 0.5*(ice_free(i,j+1) + ice_free(i,j))
          ps_ice = 0.5*(ice_cover(i,j+1) + ice_cover(i,j))
        endif
        if (l_seg /= OBC_NONE) then
          if (OBC%segment(l_seg)%open) then
            if (OBC%segment(l_seg)%direction == OBC_DIRECTION_N) then
              ps_ocn = ice_free(i,j)
              ps_ice = ice_cover(i,j)
            else
              ps_ocn = ice_free(i,j+1)
              ps_ice = ice_cover(i,j+1)
            endif
        endif ; endif
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + &
            (ps_ocn * windstr_y_water(i,J) + ps_ice * str_ice_oce_y(i,J))
      enddo ; enddo
    else
      !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
      do J=jsc-1,jec ; do i=isc,iec
        ps_ocn = 1.0 ; ps_ice = 0.0
        if (G%mask2dCv(i,J)>0.0) then
          ps_ocn = 0.5*(ice_free(i,j+1) + ice_free(i,j))
          ps_ice = 0.5*(ice_cover(i,j+1) + ice_cover(i,j))
        endif
        IOF%flux_v_ocn(i,J) = IOF%flux_v_ocn(i,J) + &
            (ps_ocn * windstr_y_water(i,J) + ps_ice * str_ice_oce_y(i,J))
      enddo ; enddo
    endif
  else
    call SIS_error(FATAL, "set_ocean_top_stress_C2: Unrecognized flux_uv_stagger.")
  endif

  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_C2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_wind_stresses_C determines the wind stresses on the ice and open ocean with
!!   a C-grid staggering of the points.
subroutine set_wind_stresses_C(FIA, ice_cover, ice_free, WindStr_x_Cu, WindStr_y_Cv, &
                               WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, G, US, max_ice_cover, OBC)
  type(fast_ice_avg_type),           intent(in)   :: FIA !< A type containing averages of fields
                                                         !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),           intent(in)   :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)   :: &
    ice_cover, &        !< The fractional ice coverage, summed across all
                        !! thickness categories [nondim], between 0 & 1.
    ice_free            !< The fractional open water [nondim], between 0 & 1.
  real, dimension(SZIB_(G),SZJ_(G)), intent(out)  :: &
    WindStr_x_Cu, &     !< Zonal wind stress averaged over the ice categories on C-grid u-points
                        !! [R Z L T-2 ~> Pa].
    WindStr_x_ocn_Cu    !< Zonal wind stress on the ice-free ocean on C-grid u-points [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G)), intent(out)  :: &
    WindStr_y_Cv, &     !< Meridional wind stress averaged over the ice categories on C-grid v-points
                        !! [R Z L T-2 ~> Pa].
    WindStr_y_ocn_Cv    !< Meridional wind stress on the ice-free ocean on C-grid v-points [R Z L T-2 ~> Pa].
  type(unit_scale_type),             intent(in)   :: US    !< A structure with unit conversion factors
  real,                              intent(in)   :: max_ice_cover !< The fractional ice coverage
                        !! that is close enough to 1 to be complete for the purpose of calculating
                        !! wind stresses [nondim].
  type(ice_OBC_type),                pointer      :: OBC  !< Open boundary structure.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    WindStr_x_A, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_A, &      ! averaged over the ice categories on an A-grid [R Z L T-2 ~> Pa].
    WindStr_x_ocn_A, &  ! Zonal (_x_) and meridional (_y_) wind stresses on the
    WindStr_y_ocn_A     ! ice-free ocean on an A-grid [R Z L T-2 ~> Pa].
  real :: weights       ! A sum of the weights around a point.
  real :: I_wts         ! 1.0 / wts or 0 if wts is 0 [nondim].
  real :: FIA_ice_cover, ice_cover_now
  integer :: i, j, isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  logical :: local_open_u_BC, local_open_v_BC

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  local_open_u_BC = .false. ; local_open_v_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_open_u_BC = OBC%open_u_BCs_exist_globally
    local_open_v_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  if (local_open_u_BC .or. local_open_v_BC) &
      call SIS_error(FATAL, "No OBCs coded yet in set_wind_stresses_C")

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
                               WindStr_x_ocn_B, WindStr_y_ocn_B, G, US, max_ice_cover)
  type(fast_ice_avg_type),            intent(in)   :: FIA !< A type containing averages of fields
                                                          !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),            intent(in)   :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)   :: &
    ice_cover, &        !< The fractional ice coverage, summed across all
                        !! thickness categories [nondim], between 0 & 1.
    ice_free            !< The fractional open water [nondim], between 0 & 1.
  real, dimension(SZIB_(G),SZJB_(G)), intent(out)  :: &
    WindStr_x_B, &      !< Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_B, &      !< averaged over the ice categories on a B-grid [R Z L T-2 ~> Pa].
    WindStr_x_ocn_B, &  !< Zonal wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    WindStr_y_ocn_B     !< Meridional wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
  type(unit_scale_type),              intent(in)   :: US    !< A structure with unit conversion factors
  real,                               intent(in)   :: max_ice_cover !< The fractional ice coverage
                        !! that is close enough to 1 to be complete for the purpose of calculating
                        !! wind stresses [nondim].

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))   :: &
    WindStr_x_A, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_A, &      ! averaged over the ice categories on an A-grid [R Z L T-2 ~> Pa].
    WindStr_x_ocn_A, &  ! Zonal (_x_) and meridional (_y_) wind stresses on the
    WindStr_y_ocn_A     ! ice-free ocean on an A-grid [R Z L T-2 ~> Pa].
  real :: weights       ! A sum of the weights around a point.
  real :: I_wts         ! 1.0 / wts or 0 if wts is 0 [nondim].
  real :: FIA_ice_cover, ice_cover_now
  integer :: i, j, isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

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
subroutine SIS_dyn_trans_register_restarts(HI, IG, param_file, CS, US, Ice_restart)
  type(hor_index_type),    intent(in) :: HI     !< The horizontal index type describing the domain
  type(ice_grid_type),     intent(in) :: IG     !< The sea-ice grid type
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(dyn_trans_CS),      pointer    :: CS     !< The control structure for the SIS_dyn_trans module
  type(unit_scale_type),   intent(in) :: US     !< A structure with unit conversion factors
  type(SIS_restart_CS),    pointer    :: Ice_restart !< The control structure for the ice restarts

!   This subroutine registers the restart variables associated with the
! the slow ice dynamics and thermodynamics.

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_dyn_trans_register_restarts called with an "//&
                            "associated control structure.")
    return
  endif
  allocate(CS)

  CS%Cgrid_dyn = .true. ; call read_param(param_file, "CGRID_ICE_DYNAMICS", CS%Cgrid_dyn)

  if (CS%Cgrid_dyn) then
    call SIS_C_dyn_register_restarts(HI, param_file, CS%SIS_C_dyn_CSp, US, Ice_restart)
  else
    call SIS_B_dyn_register_restarts(HI, param_file, CS%SIS_B_dyn_CSp, US, Ice_restart)
  endif
!  call SIS_transport_register_restarts(G, param_file, CS%SIS_transport_CSp, Ice_restart)

end subroutine SIS_dyn_trans_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_register_restarts allocates and registers any variables associated
!!      slow ice dynamics and transport that need to be included in the restart files.
subroutine SIS_dyn_trans_read_alt_restarts(CS, G, US, Ice_restart, restart_dir)
  type(dyn_trans_CS),      pointer    :: CS  !< The control structure for the SIS_dyn_trans module
  type(unit_scale_type),   intent(in) :: US  !< A structure with unit conversion factors
  type(SIS_hor_grid_type), intent(in) :: G   !< The horizontal grid type
  type(SIS_restart_CS),    pointer    :: Ice_restart !< The control structure for the ice restarts
  character(len=*),        intent(in) :: restart_dir !< The directory in which to find the restart files

  if (CS%Cgrid_dyn) then
    call SIS_C_dyn_read_alt_restarts(CS%SIS_C_dyn_CSp, G, US, Ice_restart, restart_dir)
  endif

end subroutine SIS_dyn_trans_read_alt_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_init initializes ice model data, parameters and diagnostics
!!   associated with the SIS2 dynamics and transport modules.
subroutine SIS_dyn_trans_init(Time, G, US, IG, param_file, diag, CS, output_dir, Time_init, &
                              slab_ice)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model time.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid structure
  type(unit_scale_type),       intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),         intent(in)    :: IG   !< The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(dyn_trans_CS),          pointer       :: CS   !< The control structure for the SIS_dyn_trans module
  character(len=*),            intent(in)    :: output_dir !< The directory to use for writing output
  type(time_type),             intent(in)    :: Time_Init !< Starting time of the model integration
  logical,           optional, intent(in)    :: slab_ice  !< If true, use the archaic GFDL slab ice dynamics
                                                     !!  and transport.

  ! This include declares and sets the variable "version".
#  include "version_variable.h"
  character(len=40) :: mdl = "SIS_dyn_trans" ! This module's name.
  real :: Time_unit      ! The time unit for ICE_STATS_INTERVAL [s].
  integer :: max_nts
  logical :: do_slab_ice
  logical :: debug
  real, parameter :: missing = -1e34

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
  call get_param(param_file, mdl, "CGRID_ICE_DYNAMICS", CS%Cgrid_dyn, &
                 "If true, use a C-grid discretization of the sea-ice "//&
                 "dynamics; if false use a B-grid discretization.", &
                 default=.true.)
  call get_param(param_file, mdl, "DT_ICE_DYNAMICS", CS%dt_ice_dyn, &
                 "The time step used for the slow ice dynamics, including "//&
                 "stepping the continuity equation and interactions "//&
                 "between the ice mass field and velocities.  If 0 or "//&
                 "negative the coupling time step will be used.", &
                 units="seconds", scale=US%s_to_T, default=-1.0)
  call get_param(param_file, mdl, "MERGED_CONTINUITY", CS%merged_cont, &
                 "If true, update the continuity equations for the ice, snow, "//&
                 "and melt pond water together summed across categories, with "//&
                 "proportionate fluxes for each part. Otherwise the media are "//&
                 "updated separately.", default=.false.)
  call get_param(param_file, mdl, "DT_ICE_ADVECT", CS%dt_advect, &
                 "The time step used for the advecting tracers and masses as "//&
                 "partitioned by thickness categories when merged_cont it true. "//&
                 "If 0 or negative, the coupling time step will be used.", &
                 units="seconds", scale=US%s_to_T, default=-1.0, do_not_log=.not.CS%merged_cont)
  if (.not.CS%merged_cont) CS%dt_advect = CS%dt_ice_dyn
  call get_param(param_file, mdl, "DO_RIDGING", CS%do_ridging, &
                 "If true, apply a ridging scheme to the convergent ice. "//&
                 "Otherwise, ice is compressed proportionately if the "//&
                 "concentration exceeds 1.  The original SIS2 implementation "//&
                 "is based on work by Torge Martin.", default=.false.)
  call get_param(param_file, mdl, "NSTEPS_ADV", CS%adv_substeps, &
                 "The number of advective iterations for each slow dynamics "//&
                 "time step.", default=1)
  if (CS%adv_substeps < 1) CS%adv_substeps = 1

  call get_param(param_file, mdl, "ICEBERG_WINDSTRESS_BUG", CS%berg_windstress_bug, &
                 "If true, use older code that applied an old ice-ocean "//&
                 "stress to the icebergs in place of the current air-ocean "//&
                 "stress.  This option is here for backward compatibility, "//&
                 "but should be avoided.", default=.false.)
  call get_param(param_file, mdl, "WARSAW_SUM_ORDER", CS%Warsaw_sum_order, &
                 "If true, use the order of sums in the Warsaw version of SIS2. "//&
                 "The default is the opposite of MERGED_CONTINUITY. "//&
                 "This option exists for backward compatibilty but may "//&
                 "eventually be obsoleted.", &
                 default=.not.CS%merged_cont, do_not_log=CS%merged_cont)
  if (CS%merged_cont .and. CS%Warsaw_sum_order) &
    call SIS_error(FATAL, "WARSAW_SUM_ORDER can not be true if MERGED_CONTINUITY=True.")
  call get_param(param_file, mdl, "LEMIEUX_LANDFAST", CS%lemieux_landfast, &
                   "If true, turn on Lemieux landfast ice parameterization.", default=.false., &
                   do_not_log=.true.)
  call get_param(param_file, mdl, "ITD_LANDFAST", CS%itd_landfast, &
                   "If true, turn on probabilistic landfast ice parameterization.", default=.false., &
                   do_not_log=.true.)
! if (CS%merged_cont .and. CS%itd_landfast) &
!   call SIS_error(FATAL, "ITD_LANDFAST can not be true if MERGED_CONTINUITY=True.")

  call get_param(param_file, mdl, "TIMEUNIT", Time_unit, &
                 "The time unit for ICE_STATS_INTERVAL.", &
                 units="s", default=86400.0)
  call get_param(param_file, mdl, "ICE_STATS_INTERVAL", CS%ice_stats_interval, &
                 "The interval in units of TIMEUNIT between writes of the "//&
                 "globally summed ice statistics and conservation checks.", &
                 default=real_to_time(86400.0), timeunit=Time_unit)

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
  call get_param(param_file, mdl, "VERBOSE", CS%verbose, &
                 "If true, write out verbose diagnostics.", default=.false., &
                 debuggingParam=.true.)

  CS%complete_ice_cover = 1.0 - 2.0*epsilon(CS%complete_ice_cover)

  if (.not.(do_slab_ice)) then
    CS%complete_ice_cover = 1.0 - 2.0*max(1,IG%CatIce)*epsilon(CS%complete_ice_cover)
    if (CS%Cgrid_dyn) then
      call SIS_C_dyn_init(CS%Time, G, US, param_file, CS%diag, CS%SIS_C_dyn_CSp, CS%ntrunc)
    else
      call SIS_B_dyn_init(CS%Time, G, US, param_file, CS%diag, CS%SIS_B_dyn_CSp)
    endif
    if (CS%merged_cont) then
      call SIS_transport_init(CS%Time, G, IG, US, param_file, CS%diag, CS%SIS_transport_CSp, &
                              continuity_CSp=CS%continuity_CSp, cover_trans_CSp=CS%cover_trans_CSp)
    else
      call SIS_transport_init(CS%Time, G, IG, US, param_file, CS%diag, CS%SIS_transport_CSp, &
                              continuity_CSp=CS%continuity_CSp)
    endif

    call alloc_cell_average_state_type(CS%CAS, G%HI, IG, CS%SIS_transport_CSp)

    if (.not.associated(CS%DS2d)) allocate(CS%DS2d)
    CS%DS2d%ridge_rate_count = 0.
    if (CS%do_ridging) call safe_alloc(CS%DS2d%avg_ridge_rate, G%isd, G%ied, G%jsd, G%jed)

    if (CS%merged_cont) then
      CS%DS2d%nts = 0 ; CS%DS2d%max_nts = 0
      call safe_alloc(CS%DS2d%mi_sum, G%isd, G%ied, G%jsd, G%jed)
      call safe_alloc(CS%DS2d%ice_cover, G%isd, G%ied, G%jsd, G%jed)
      max_nts = CS%adv_substeps
      if ((CS%dt_ice_dyn > 0.0) .and. (CS%dt_advect > CS%dt_ice_dyn)) &
        max_nts = CS%adv_substeps * max(CEILING(CS%dt_advect/CS%dt_ice_dyn - 1e-6), 1)
      call increase_max_tracer_step_memory(CS%DS2d, G, max_nts)

      call safe_alloc(CS%DS2d%u_ice_C, G%IsdB, G%IedB, G%jsd, G%jed)
      call safe_alloc(CS%DS2d%v_ice_C, G%isd, G%ied, G%JsdB, G%JedB)
      if (.not.CS%Cgrid_dyn) then
        call safe_alloc(CS%DS2d%u_ice_B, G%IsdB, G%IedB, G%JsdB, G%JedB)
        call safe_alloc(CS%DS2d%v_ice_B, G%IsdB, G%IedB, G%JsdB, G%JedB)
      endif
    endif

  endif

  call SIS_sum_output_init(G, param_file, output_dir, Time_Init, US, &
                           CS%sum_output_CSp, CS%ntrunc)

  CS%write_ice_stats_time = Time_Init + CS%ice_stats_interval * &
      (1 + (Time - Time_init) / CS%ice_stats_interval)

  ! Stress diagnostics that are specific to the C-grid or B-grid dynamics of the ice model
  if (CS%Cgrid_dyn) then
    CS%id_fax = register_diag_field('ice_model', 'FA_X', diag%axesCu1, Time, &
               'Air stress on ice on C-grid - x component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
                missing_value=missing, interp_method='none')
    CS%id_fay = register_diag_field('ice_model', 'FA_Y', diag%axesCv1, Time, &
               'Air stress on ice on C-grid - y component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
               missing_value=missing, interp_method='none')
  else
    CS%id_fax = register_diag_field('ice_model', 'FA_X', diag%axesB1, Time, &
               'air stress on ice - x component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
               missing_value=missing, interp_method='none')
    CS%id_fay = register_diag_field('ice_model', 'FA_Y', diag%axesB1, Time, &
               'air stress on ice - y component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
               missing_value=missing, interp_method='none')
  endif

  call register_ice_state_diagnostics(Time, IG, US, param_file, diag, CS%IDs)

  iceClock4 = cpu_clock_id( '  Ice: slow: dynamics', grain=CLOCK_ROUTINE )
  iceClocka = cpu_clock_id( '       slow: ice_dynamics', grain=CLOCK_LOOP )
  iceClockb = cpu_clock_id( '       slow: comm/cut check ', grain=CLOCK_LOOP )
  iceClockc = cpu_clock_id( '       slow: diags', grain=CLOCK_LOOP )
  iceClock8 = cpu_clock_id( '  Ice: slow: transport', grain=CLOCK_ROUTINE )
  iceClock9 = cpu_clock_id( '  Ice: slow: thermodyn diags', grain=CLOCK_ROUTINE )

  call callTree_leave("SIS_dyn_trans_init()")

end subroutine SIS_dyn_trans_init


!> Increase the memory available to store total ice and snow masses and mass fluxes for tracer advection.
!! Any data already stored in the fluxes is copied over to the new arrays.
subroutine increase_max_tracer_step_memory(DS2d, G, max_nts)
  type(dyn_state_2d),      intent(inout) :: DS2d   !< The control structure for the SIS_dyn_trans module
  type(SIS_hor_grid_type), intent(in)    :: G    !< The horizontal grid structure
  integer,                 intent(in)    :: max_nts !< The new maximum number of masses and mass fluxes
                                                 !! that can be stored for tracer advection.

  real, allocatable :: tmp_array(:,:,:)
  integer :: nts_prev

  if (DS2d%max_nts >= max(max_nts,1)) return

  nts_prev = DS2d%nts
  DS2d%max_nts = max(max_nts,1)

  if (allocated(DS2d%mca_step)) then
    allocate(tmp_array(G%isd:G%ied, G%jsd:G%jed, 0:nts_prev))
    tmp_array(:,:,0:nts_prev) = DS2d%mca_step(:,:,0:nts_prev)
    deallocate(DS2d%mca_step)
    allocate(DS2d%mca_step(G%isd:G%ied, G%jsd:G%jed, 0:DS2d%max_nts), source=0.0)
    ! Copy over the data that had been set before.
    DS2d%mca_step(:,:,0:nts_prev) = tmp_array(:,:,0:nts_prev)
    deallocate(tmp_array)
  else
    allocate(DS2d%mca_step(G%isd:G%ied, G%jsd:G%jed, 0:DS2d%max_nts), source=0.0)
  !  This is the equivalent for when the 6 argument version of safe_alloc is available.
  !      call safe_alloc(DS2d%mca_step, G%isd, G%ied, G%jsd, G%jed, 0, DS2d%max_nts)
  endif

  if (allocated(DS2d%uh_step)) then
    if (nts_prev > 0) then
      allocate(tmp_array(G%IsdB:G%IedB, G%jsd:G%jed, nts_prev))
      tmp_array(:,:,1:nts_prev) = DS2d%uh_step(:,:,1:nts_prev)
    endif
    deallocate(DS2d%uh_step)
    call safe_alloc(DS2d%uh_step, G%IsdB, G%IedB, G%jsd, G%jed, DS2d%max_nts)
    if (nts_prev > 0) then ! Copy over the data that had been set before.
      DS2d%uh_step(:,:,1:nts_prev) = tmp_array(:,:,1:nts_prev)
      deallocate(tmp_array)
    endif
  else
    call safe_alloc(DS2d%uh_step, G%IsdB, G%IedB, G%jsd, G%jed, DS2d%max_nts)
  endif

  if (allocated(DS2d%vh_step)) then
    if (nts_prev > 0) then
      allocate(tmp_array(G%isd:G%ied, G%JsdB:G%JedB, nts_prev))
      tmp_array(:,:,1:nts_prev) = DS2d%vh_step(:,:,1:nts_prev)
    endif
    deallocate(DS2d%vh_step)
    call safe_alloc(DS2d%vh_step, G%isd, G%ied, G%JsdB, G%JedB, DS2d%max_nts)
    if (nts_prev > 0) then ! Copy over the data that had been set before.
      DS2d%vh_step(:,:,1:nts_prev) = tmp_array(:,:,1:nts_prev)
      deallocate(tmp_array)
    endif
  else
    call safe_alloc(DS2d%vh_step, G%isd, G%ied, G%JsdB, G%JedB, DS2d%max_nts)
  endif

end subroutine increase_max_tracer_step_memory

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_transport_CS returns a pointer to the SIS_transport_CS type that
!!  the dyn_trans_CS points to.
function SIS_dyn_trans_transport_CS(CS) result(transport_CSp)
  type(dyn_trans_CS),     pointer :: CS    !< The control structure for the SIS_dyn_trans module
  type(SIS_transport_CS), pointer :: transport_CSp !< The SIS_transport_CS type used by SIS_dyn_trans

  transport_CSp => CS%SIS_transport_CSp
end function SIS_dyn_trans_transport_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_sum_output_CS returns a pointer to the sum_out_CS type that
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
                                     !! is deallocated here

  if (associated(CS%DS2d)) then
    if (allocated(CS%DS2d%mi_sum)) deallocate(CS%DS2d%mi_sum)
    if (allocated(CS%DS2d%ice_cover)) deallocate(CS%DS2d%ice_cover)
    if (allocated(CS%DS2d%mca_step)) deallocate(CS%DS2d%mca_step)
    if (allocated(CS%DS2d%uh_step)) deallocate(CS%DS2d%uh_step)
    if (allocated(CS%DS2d%vh_step)) deallocate(CS%DS2d%vh_step)
    if (allocated(CS%DS2d%u_ice_B)) deallocate(CS%DS2d%u_ice_B)
    if (allocated(CS%DS2d%v_ice_B)) deallocate(CS%DS2d%v_ice_B)
    if (allocated(CS%DS2d%u_ice_C)) deallocate(CS%DS2d%u_ice_C)
    if (allocated(CS%DS2d%v_ice_C)) deallocate(CS%DS2d%v_ice_C)
    if (allocated(CS%DS2d%avg_ridge_rate)) deallocate(CS%DS2d%avg_ridge_rate)
    deallocate(CS%DS2d)
  endif

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
