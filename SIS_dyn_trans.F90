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
! scales of the couplng or the interactions with the ocean due to ice dynamics !
! and lateral transport.                                                       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_dyn_trans

use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_debugging,     only : chksum, Bchksum, hchksum
use SIS_debugging,     only : vec_chksum_A, vec_chksum_B, vec_chksum_C
use SIS_sum_output, only : write_ice_statistics, SIS_sum_output_init, SIS_sum_out_CS

use mpp_domains_mod,  only  : domain2D
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : fill_symmetric_edges !, MOM_domains_init, clone_MOM_domain
! use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, log_version, param_file_type
use MOM_hor_index, only : hor_index_type ! , hor_index_init
! use MOM_string_functions, only : uppercase
use MOM_EOS, only : EOS_type, calculate_density_derivs

use fms_mod, only : clock_flag_default
! use fms_io_mod, only : restore_state, query_initialized
use fms_io_mod, only : register_restart_field, restart_file_type
use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
use mpp_mod, only : CLOCK_COMPONENT, CLOCK_LOOP, CLOCK_ROUTINE

use time_manager_mod, only : time_type, time_type_to_real, get_date, get_time
use time_manager_mod, only : set_date, set_time, operator(+), operator(-)
use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)

use SIS_types, only : ice_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_types, only : ocean_sfc_state_type
use SIS_types, only : IST_chksum,  IST_bounds_check
use SIS_utils, only : get_avg, post_avg, ice_line !, ice_grid_chksum
use SIS_hor_grid, only : SIS_hor_grid_type

use ice_grid, only : ice_grid_type

use SIS2_ice_thm,  only: get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,  only: enth_from_TS, Temp_from_En_S
use SIS_dyn_bgrid, only: SIS_B_dyn_CS, SIS_B_dynamics, SIS_B_dyn_init
use SIS_dyn_bgrid, only: SIS_B_dyn_register_restarts, SIS_B_dyn_end
use SIS_dyn_cgrid, only: SIS_C_dyn_CS, SIS_C_dynamics, SIS_C_dyn_init
use SIS_dyn_cgrid, only: SIS_C_dyn_register_restarts, SIS_C_dyn_end
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS
use SIS_transport, only : ice_transport, SIS_transport_init, SIS_transport_end
use SIS_transport, only : SIS_transport_CS

use ice_bergs,     only: icebergs, icebergs_run, icebergs_init, icebergs_end

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_dynamics_trans, update_icebergs, dyn_trans_CS
public :: SIS_dyn_trans_register_restarts, SIS_dyn_trans_init, SIS_dyn_trans_end
public :: SIS_dyn_trans_transport_CS, SIS_dyn_trans_sum_output_CS
public :: post_ocean_sfc_diagnostics, post_ice_state_diagnostics

type dyn_trans_CS ; private
  logical :: Cgrid_dyn ! If true use a C-grid discretization of the
                       ! sea-ice dynamics.
  logical :: specified_ice  ! If true, the sea ice is specified and there is
                            ! no need for ice dynamics.
  real    :: dt_ice_dyn   ! The time step used for the slow ice dynamics, including
                          ! stepping the continuity equation and interactions
                          ! between the ice mass field and velocities, in s. If
                          ! 0 or negative, the coupling time step will be used.
  logical :: do_ridging   !   If true, apply a ridging scheme to the convergent
                          ! ice.  The original SIS2 implementation is based on
                          ! work by Torge Martin.  Otherwise, ice is compressed
                          ! proportionately if the concentration exceeds 1.
  logical :: berg_windstress_bug  ! If true, use older code that applied an old
                          ! ice-ocean stress to the icebergs in place of the
                          ! current air-ice stress.  This option is here for
                          ! backward compatibility, but should be avoided.

  logical :: debug        ! If true, write verbose checksums for debugging purposes.
  logical :: column_check ! If true, enable the heat check column by column.
  real    :: imb_tol      ! The tolerance for imbalances to be flagged by
                          ! column_check, nondim.
  logical :: bounds_check ! If true, check for sensible values of thicknesses
                          ! temperatures, fluxes, etc.
  logical :: verbose      ! A flag to control the printing of an ice-diagnostic
                          ! message.  When true, this will slow the model down.

  integer :: ntrunc = 0   ! The number of times the velocity has been truncated
                          ! since the last call to write_ice_statistics.

  integer :: n_calls = 0  ! The number of times SIS_dynamics_trans has been called.
  type(time_type) :: ice_stats_interval ! The interval between writes of the
                          ! globally summed ice statistics and conservation checks.
  type(time_type) :: write_ice_stats_time ! The next time to write out the ice statistics.

  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                                   ! timing of diagnostic output.

  integer :: id_fax=-1, id_fay=-1, id_xprt=-1, id_mib=-1, id_mi=-1

  ! These are the diagnostic ids for describing the ice state.
  integer, dimension(:), allocatable :: id_t, id_sal
  integer :: id_cn=-1, id_hi=-1, id_hp=-1, id_hs=-1, id_tsn=-1, id_ext=-1 ! id_hp mw/new
  integer :: id_t_iceav=-1, id_s_iceav=-1, id_e2m=-1
  integer :: id_rdgr=-1 ! These do not exist yet: id_rdgf=-1, id_rdgo=-1, id_rdgv=-1

  integer :: id_simass=-1, id_sisnmass=-1, id_sivol=-1
  integer :: id_siconc=-1, id_sithick=-1, id_sisnconc=-1, id_sisnthick=-1
  integer :: id_siu=-1, id_siv=-1, id_sispeed=-1, id_sitimefrac=-1

  type(SIS_B_dyn_CS), pointer     :: SIS_B_dyn_CSp => NULL()
  type(SIS_C_dyn_CS), pointer     :: SIS_C_dyn_CSp => NULL()
  type(SIS_transport_CS), pointer :: SIS_transport_CSp => NULL()
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()
  logical :: module_is_initialized = .false.
end type dyn_trans_CS

integer :: iceClock4, iceClock8, iceClock9, iceClocka, iceClockb, iceClockc

contains

subroutine update_icebergs(IST, OSS, IOF, FIA, icebergs_CS, dt_slow, G, IG, CS)
  type(ice_state_type),       intent(inout) :: IST
  type(ocean_sfc_state_type), intent(in)    :: OSS
  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(ice_ocean_flux_type),  intent(inout) :: IOF
  real,                       intent(in)    :: dt_slow
  type(icebergs),             pointer       :: icebergs_CS
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(inout) :: IG
  type(dyn_trans_CS),         pointer       :: CS

  real, dimension(SZI_(G),SZJ_(G))   :: &
    hi_avg            ! The area-weighted average ice thickness, in m.
  real, dimension(G%isc:G%iec, G%jsc:G%jec)   :: &
    windstr_x, &      ! The area-weighted average ice thickness, in Pa.
    windstr_y         ! The area-weighted average ice thickness, in Pa.
  real :: rho_ice     ! The nominal density of sea ice in kg m-3.
  real :: H_to_m_ice  ! The specific volume of ice times the conversion factor
                      ! from thickness units, in m H-1.
  integer :: stress_stagger
  integer :: i, j, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  H_to_m_ice = IG%H_to_kg_m2 / rho_ice
  call get_avg(IST%mH_ice, IST%part_size(:,:,1:), hi_avg, wtd=.true.)
  hi_avg(:,:) = hi_avg(:,:) * H_to_m_Ice

  if (CS%berg_windstress_bug) then
    !  This code reproduces a long-standing bug, in that the old ice-ocean
    !   stresses are being passed in place of the wind stresses on the icebergs.
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
! SIS_dynamics_trans - do ice dynamics and mass and tracer transport           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_dynamics_trans(IST, OSS, FIA, IOF, dt_slow, CS, icebergs_CS, G, IG, tracer_CSp)

  type(ice_state_type),       intent(inout) :: IST
  type(ocean_sfc_state_type), intent(in)    :: OSS
  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(ice_ocean_flux_type),  intent(inout) :: IOF
  real,                       intent(in)    :: dt_slow
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(inout) :: IG
  type(dyn_trans_CS),         pointer       :: CS
  type(icebergs),             pointer       :: icebergs_CS
  type(SIS_tracer_flow_control_CS), pointer :: tracer_CSp

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: h2o_chg_xprt, mass, mass_ice, mass_snow, tmp2d
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,IG%NkIce) :: &
    temp_ice    ! A diagnostic array with the ice temperature in degC.
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    temp_snow   ! A diagnostic array with the snow temperature in degC.
  real, dimension(SZI_(G),SZJ_(G))   :: &
    ms_sum, mi_sum, &   ! Masses of snow and ice per unit total area, in kg m-2.
    ice_free, &         ! The fractional open water; nondimensional, between 0 & 1.
    ice_cover, &        ! The fractional ice coverage, summed across all
                        ! thickness categories; nondimensional, between 0 & 1.
    WindStr_x_A, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_A, &      ! averaged over the ice categories on an A-grid, in Pa.
    WindStr_x_ocn_A, &  ! Zonal (_x_) and meridional (_y_) wind stresses on the
    WindStr_y_ocn_A     ! ice-free ocean on an A-grid, in Pa.
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

  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBx ! An temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBy ! An temporary array for diagnostics.

  real :: weights  ! A sum of the weights around a point.
  real :: I_wts    ! 1.0 / wts or 0 if wts is 0, nondim.
  real :: ps_vel   ! The fractional thickness catetory coverage at a velocity point.

  real :: dt_slow_dyn
  real :: max_ice_cover, FIA_ice_cover, ice_cover_now
  integer :: ndyn_steps
  real :: Idt_slow
  integer :: i, j, k, l, m, isc, iec, jsc, jec, ncat, NkIce, nds
  integer :: isd, ied, jsd, jed
  integer :: iyr, imon, iday, ihr, imin, isec

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    rdg_frac    ! fraction of ridged ice per category
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rdg_open, & ! formation rate of open water due to ridging
    rdg_vosh, & ! rate of ice mass shifted from level to ridged ice
!!   rdg_s2o, &  ! snow mass [kg m-2] dumped into ocean during ridging
    rdg_rate, & ! Niki: Where should this come from?
    snow2ocn
  real    :: tmp3  ! This is a bad name - make it more descriptive!

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; NkIce = IG%NkIce
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  if (CS%specified_ice) then
    ndyn_steps = 0 ; dt_slow_dyn = 0.0
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,WindStr_x_A,WindStr_y_A,  &
!$OMP                                  ice_cover,ice_free,WindStr_x_ocn_A,       &
!$OMP                                  WindStr_y_ocn_A,FIA)
    do j=jsd,jed
      do i=isd,ied
        WindStr_x_ocn_A(i,j) = FIA%WindStr_ocn_x(i,j)
        WindStr_y_ocn_A(i,j) = FIA%WindStr_ocn_y(i,j)
        ice_cover(i,j) = FIA%ice_cover(i,j) ; ice_free(i,j) = FIA%ice_free(i,j)
        WindStr_x_A(i,j) = FIA%WindStr_x(i,j) ; WindStr_y_A(i,j) = FIA%WindStr_y(i,j)
      enddo
    enddo
  else
    ndyn_steps = 1
    if ((CS%dt_ice_dyn > 0.0) .and. (CS%dt_ice_dyn < dt_slow)) &
      ndyn_steps = max(CEILING(dt_slow/CS%dt_ice_dyn - 0.000001), 1)
    dt_slow_dyn = dt_slow / ndyn_steps
  endif
  IOF%stress_count = 0

  CS%n_calls = CS%n_calls + 1

  if (CS%id_xprt>0) then
    ! Store values to determine the ice and snow mass change due to transport.
    h2o_chg_xprt(:,:) = 0.0
  endif

  do nds=1,ndyn_steps

    call enable_SIS_averaging(dt_slow_dyn, CS%Time - set_time(int((ndyn_steps-nds)*dt_slow_dyn)), CS%diag)

    ! Correct the wind stresses for changes in the fractional ice-coverage.
    ice_cover(:,:) = 0.0
    max_ice_cover = 1.0 - 2.0*ncat*epsilon(max_ice_cover)
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,ncat,ice_cover,IST,FIA,ice_free, &
!$OMP                                  WindStr_x_A,WindStr_y_A,WindStr_x_ocn_A, &
!$OMP                                  max_ice_cover, WindStr_y_ocn_A) &
!$OMP                           private(FIA_ice_cover, ice_cover_now)
    do j=jsd,jed
      do k=1,ncat ; do i=isd,ied
        ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
      enddo ; enddo
      do i=isd,ied
        ! The use of these limits prevents the use of the ocean wind stresses
        ! there is actually no open ocean and hence there may be no valid ocean
        ! stresses.  This can occur when ice_cover ~= 1 for both states, but
        ! they are not exactly 1.0 due to roundoff in the sum above.
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

        ice_free(i,j) = IST%part_size(i,j,0)
        if (ice_free(i,j) <= FIA%ice_free(i,j)) then
          WindStr_x_ocn_A(i,j) = FIA%WindStr_ocn_x(i,j)
          WindStr_y_ocn_A(i,j) = FIA%WindStr_ocn_y(i,j)
        else
          WindStr_x_ocn_A(i,j) = ((ice_free(i,j)-FIA%ice_free(i,j))*FIA%WindStr_x(i,j) + &
                              FIA%ice_free(i,j)*FIA%WindStr_ocn_x(i,j)) / ice_free(i,j)
          WindStr_y_ocn_A(i,j) = ((ice_free(i,j)-FIA%ice_free(i,j))*FIA%WindStr_y(i,j) + &
                              FIA%ice_free(i,j)*FIA%WindStr_ocn_y(i,j)) / ice_free(i,j)
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
    if (CS%Cgrid_dyn) then
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

      if (CS%debug) then
        call IST_chksum("Before SIS_C_dynamics", IST, G, IG)
        call hchksum(IST%part_size(:,:,0), "ps(0) before SIS_C_dynamics", G%HI)
        call hchksum(ms_sum, "ms_sum before SIS_C_dynamics", G%HI)
        call hchksum(mi_sum, "mi_sum before SIS_C_dynamics", G%HI)
        call hchksum(OSS%sea_lev, "sea_lev before SIS_C_dynamics", G%HI, haloshift=1)
        call hchksum(ice_cover, "ice_cover before SIS_C_dynamics", G%HI, haloshift=1)
        call vec_chksum_C("[uv]_ocn before SIS_C_dynamics", OSS%u_ocn_C, OSS%v_ocn_C, G, halos=1)
        call vec_chksum_C("WindStr_[xy] before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G, halos=1)
        call vec_chksum_A("WindStr_[xy]_A before SIS_C_dynamics", WindStr_x_A, WindStr_y_A, G, halos=1)
     endif

      call mpp_clock_begin(iceClocka)
      !### Ridging needs to be added with C-grid dynamics.
      call SIS_C_dynamics(1.0-IST%part_size(:,:,0), ms_sum, mi_sum, IST%u_ice_C, IST%v_ice_C, &
                        OSS%u_ocn_C, OSS%v_ocn_C, &
                        WindStr_x_Cu, WindStr_y_Cv, OSS%sea_lev, str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, &
                        dt_slow_dyn, G, CS%SIS_C_dyn_CSp)
      call mpp_clock_end(iceClocka)

      if (CS%debug) then
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
      if (CS%id_fax>0) call post_data(CS%id_fax, WindStr_x_Cu, CS%diag)
      if (CS%id_fay>0) call post_data(CS%id_fay, WindStr_y_Cv, CS%diag)

      if (CS%debug) call IST_chksum("Before set_ocean_top_stress_Cgrid", IST, G, IG)

      call set_ocean_top_stress_Cgrid(IOF, WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, &
                                      str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, IST%part_size, G, IG)
      if (CS%debug) call IST_chksum("After set_ocean_top_stress_Cgrid", IST, G, IG)
      call mpp_clock_end(iceClockc)

    else ! B-grid dynamics.

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

      if (CS%debug) then
        call IST_chksum("Before ice_dynamics", IST, G, IG)
        call hchksum(IST%part_size(:,:,0), "ps(0) before ice_dynamics", G%HI)
        call hchksum(ms_sum, "ms_sum before ice_dynamics", G%HI)
        call hchksum(mi_sum, "mi_sum before ice_dynamics", G%HI)
        call hchksum(OSS%sea_lev, "sea_lev before ice_dynamics", G%HI, haloshift=1)
        call vec_chksum_B("uv_ocn before ice_dynamics", OSS%u_ocn_B, OSS%v_ocn_B, G)
        call vec_chksum_B("WindStr_x_B before ice_dynamics", WindStr_x_B, WindStr_y_B, G, halos=1)
      endif

      rdg_rate(:,:) = 0.0
      call mpp_clock_begin(iceClocka)
      call SIS_B_dynamics(1.0-IST%part_size(:,:,0), ms_sum, mi_sum, IST%u_ice_B, IST%v_ice_B, &
                        OSS%u_ocn_B, OSS%v_ocn_B, WindStr_x_B, WindStr_y_B, OSS%sea_lev, &
                        str_x_ice_ocn_B, str_y_ice_ocn_B, CS%do_ridging, &
                        rdg_rate(isc:iec,jsc:jec), dt_slow_dyn, G, CS%SIS_B_dyn_CSp)
      call mpp_clock_end(iceClocka)

      if (CS%debug) then
        call IST_chksum("After ice_dynamics", IST, G, IG)
      endif

      call mpp_clock_begin(iceClockb)
      call pass_vector(IST%u_ice_B, IST%v_ice_B, G%Domain, stagger=BGRID_NE)
      call mpp_clock_end(iceClockb)

      call mpp_clock_begin(iceClockc)
      !
      ! Dynamics diagnostics
      !
      if ((CS%id_fax>0) .or. (CS%id_fay>0)) then
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

        if (CS%id_fax>0) call post_data(CS%id_fax, diagVarBx, CS%diag)
        if (CS%id_fay>0) call post_data(CS%id_fay, diagVarBy, CS%diag)
      endif

      if (CS%debug) call IST_chksum("Before set_ocean_top_stress_Bgrid", IST, G, IG)
      call set_ocean_top_stress_Bgrid(IOF, WindStr_x_ocn_B, WindStr_y_ocn_B, &
                                      str_x_ice_ocn_B, str_y_ice_ocn_B, IST%part_size, G, IG)
      if (CS%debug) call IST_chksum("After set_ocean_top_stress_Bgrid", IST, G, IG)
      call mpp_clock_end(iceClockc)
    endif ! End of B-grid dynamics


    call mpp_clock_end(iceClock4)

    call enable_SIS_averaging(dt_slow_dyn, CS%Time - set_time(int((ndyn_steps-nds)*dt_slow_dyn)), CS%diag)
    !
    ! Do ice transport ... all ocean fluxes have been calculated by now.
    !
    call mpp_clock_begin(iceClock8)

    if (CS%id_xprt>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,h2o_chg_xprt,IST,G,IG)
      do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
        h2o_chg_xprt(i,j) = h2o_chg_xprt(i,j) - IST%part_size(i,j,k) * &
                          IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
      enddo ; enddo ; enddo
    endif

    if (CS%debug) then
      call IST_chksum("Before ice_transport", IST, G, IG)
    endif

    if (CS%Cgrid_dyn) then
      call ice_transport(IST%part_size, IST%mH_ice, IST%mH_snow, IST%mH_pond, &
                         IST%u_ice_C, IST%v_ice_C, IST%TrReg, &
                         dt_slow_dyn, G, IG, CS%SIS_transport_CSp,&
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

      call ice_transport(IST%part_size, IST%mH_ice, IST%mH_snow, IST%mH_pond, &
                         uc, vc, IST%TrReg, &
                         dt_slow_dyn, G, IG, CS%SIS_transport_CSp, &
                         IST%rdg_mice, snow2ocn, rdg_rate, rdg_open, rdg_vosh)
    endif
    if (CS%column_check) &
      call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
                                message="      Post_transport")! , check_column=.true.)

    if (CS%id_xprt>0) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,h2o_chg_xprt,IST,G,IG)
      do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      h2o_chg_xprt(i,j) = h2o_chg_xprt(i,j) + IST%part_size(i,j,k) * &
                        IG%H_to_kg_m2 * (IST%mH_snow(i,j,k) + IST%mH_ice(i,j,k))
    enddo ; enddo ; enddo ; endif

    call mpp_clock_end(iceClock8)

  enddo ! nds=1,ndyn_steps
  call finish_ocean_top_stresses(IOF, G%HI)

  ! Add snow mass dumped into ocean to flux of frozen precipitation:
  !### WARNING - rdg_s2o is never calculated!!!
!  if (CS%do_ridging) then ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
!    FIA%fprec_top(i,j,k) = FIA%fprec_top(i,j,k) + rdg_s2o(i,j)/dt_slow
!  enddo ; enddo ; enddo ; endif

  call mpp_clock_begin(iceClock9)

  ! Set appropriate surface quantities in categories with no ice.
  if (allocated(IST%t_surf)) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,OSS)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)<=0.0) &
      IST%t_surf(i,j,k) = T_0degC + OSS%T_fr_ocn(i,j)
    enddo ; enddo ; enddo
  endif

  if (CS%bounds_check) call IST_bounds_check(IST, G, IG, "After ice_transport", OSS=OSS)
  if (CS%debug) call IST_chksum("After ice_transport", IST, G, IG)

  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)

  call post_ice_state_diagnostics(CS, IST, OSS, IOF, dt_slow, G, IG, CS%diag, &
                                  h2o_chg_xprt=h2o_chg_xprt, rdg_rate=rdg_rate)

  call disable_SIS_averaging(CS%diag)

  if (CS%verbose) then
    call get_date(CS%Time, iyr, imon, iday, ihr, imin, isec)
    call get_time(CS%Time-set_date(iyr,1,1,0,0,0),isec,iday)
    call ice_line(iyr, iday+1, isec, IST%part_size(isc:iec,jsc:jec,0), &
                  OSS%SST_C(:,:), G)
  endif

  call mpp_clock_end(iceClock9)

  if (CS%debug) then
    call IST_chksum("End SIS_dynamics_trans", IST, G, IG)
  endif

  if (CS%bounds_check) then
    call IST_bounds_check(IST, G, IG, "End of SIS_dynamics_trans", OSS=OSS)
  endif

  if (CS%Time + set_time(int(floor(0.5*dt_slow+0.5))) > CS%write_ice_stats_time) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp, &
        tracer_CSp = tracer_CSp)
    CS%write_ice_stats_time = CS%write_ice_stats_time + CS%ice_stats_interval
  elseif (CS%column_check) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp)
  endif

end subroutine SIS_dynamics_trans

subroutine post_ice_state_diagnostics(CS, IST, OSS, IOF, dt_slow, G, IG, diag, &
                                       h2o_chg_xprt, rdg_rate)
  type(ice_state_type),       intent(inout) :: IST
  type(ocean_sfc_state_type), intent(in)    :: OSS
!  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(ice_ocean_flux_type),  intent(in)    :: IOF
  real,                       intent(in)    :: dt_slow
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(inout) :: IG
  type(dyn_trans_CS),         pointer       :: CS
  type(SIS_diag_ctrl),        pointer       :: diag
  real, dimension(G%isc:G%iec,G%jsc:G%jec), optional, intent(in) :: h2o_chg_xprt
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: rdg_rate

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: mass, mass_ice, mass_snow, tmp2d
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,IG%NkIce) :: &
    temp_ice    ! A diagnostic array with the ice temperature in degC.
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    temp_snow   ! A diagnostic array with the snow temperature in degC.
  real, dimension(SZI_(G),SZJ_(G))   :: diagVar ! An temporary array for diagnostics.
  real, dimension(IG%NkIce) :: S_col ! Specified thermodynamic salinity of each
                                     ! ice layer if spec_thermo_sal is true.
  real :: rho_ice  ! The nominal density of sea ice in kg m-3.
  real :: rho_snow ! The nominal density of snow in kg m-3.
  real :: enth_units, I_enth_units
  real :: I_Nk        ! The inverse of the number of layers in the ice.
  real :: Idt_slow ! The inverse of the thermodynamic step, in s-1.
  logical :: spec_thermo_sal
  logical :: do_temp_diags
  integer :: i, j, k, l, m, isc, iec, jsc, jec, ncat, NkIce ! , nds
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
!  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ;
  NkIce = IG%NkIce
  I_Nk = 1.0 / NkIce
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  ! Sum the concentration weighted mass for diagnostics.
  if (CS%id_mi>0 .or. CS%id_mib>0) then
    mass_ice(:,:) = 0.0
    mass_snow(:,:) = 0.0
    mass(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,mass,mass_ice,mass_snow,G,IST,IG)
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
          mass(i,j) = (mass(i,j) + IOF%mass_berg(i,j)) ! Add icebergs mass in kg/m^2
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
   call post_ocean_sfc_diagnostics(OSS, dt_slow, G, diag)

  if (CS%id_xprt>0 .and. present(h2o_chg_xprt)) then
    call post_data(CS%id_xprt, h2o_chg_xprt(isc:iec,jsc:jec)*864e2*365/dt_slow, &
                   diag)
  endif
  if (CS%id_e2m>0) then
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
    call post_data(CS%id_e2m,  tmp2d(:,:), diag)
  endif

  if (CS%do_ridging) then
  !TOM> preparing output field fraction of ridged ice rdg_frac = (ridged ice volume) / (total ice volume)
  !     in each category; IST%rdg_mice is ridged ice mass per unit total area throughout the code.
!     if (CS%id_rdgf>0) then
! !$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,G,rdg_frac,IG) &
! !$OMP                          private(tmp3)
!       do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
!         tmp3 = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
!         if (tmp3*IG%H_to_kg_m2 > Rho_Ice*1.e-5) then   ! 1 mm ice thickness x 1% ice concentration
!           rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp3
!         else
!           rdg_frac(i,j,k) = 0.0
!         endif
!       enddo ; enddo ; enddo
!       call post_data(CS%id_rdgf, rdg_frac(isc:iec,jsc:jec), diag)
!     endif

    if (CS%id_rdgr>0 .and. present(rdg_rate)) &
      call post_data(CS%id_rdgr, rdg_rate(isc:iec,jsc:jec), diag)
!    if (CS%id_rdgo>0) call post_data(CS%id_rdgo, rdg_open(isc:iec,jsc:jec), diag)
!    if (CS%id_rdgv>0) then
!      do j=jsc,jec ; do i=isc,iec
!        tmp2d(i,j) = rdg_vosh(i,j) * G%areaT(i,j) * G%mask2dT(i,j)
!      enddo ; enddo
!      call post_data(CS%id_rdgv, tmp2d, diag)
!    endif
  endif

end subroutine post_ice_state_diagnostics

subroutine post_ocean_sfc_diagnostics(OSS, dt_slow, G, diag)
  type(ocean_sfc_state_type), intent(in)    :: OSS
  real,                       intent(in)    :: dt_slow
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(SIS_diag_ctrl),        pointer       :: diag

  real :: Idt_slow ! The inverse of the thermodynamic step, in s-1.
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

end subroutine post_ocean_sfc_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! finish_ocean_top_stresses - Finish setting the ice-ocean stresses by dividing!
!   them through the stresses by the number of times they have been augmented. !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine finish_ocean_top_stresses(IOF, HI)
  type(hor_index_type), intent(in)    :: HI
  type(ice_ocean_flux_type), intent(inout) :: IOF

  real :: I_count
  integer :: i, j, isc, iec, jsc, jec
  isc = HI%isc ; iec = HI%iec ; jsc = HI%jsc ; jec = HI%jec

  if (IOF%stress_count > 1) then
    I_count = 1.0 / IOF%stress_count
    do j=jsc,jec ; do i=isc,iec
      IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) * I_count
      IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) * I_count
    enddo ; enddo
  endif

end subroutine finish_ocean_top_stresses

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ocean_top_stress_Bgrid - Calculate the stresses on the ocean integrated  !
!   across all the thickness categories with the appropriate staggering, and   !
!   store them in the public ice data type for use by the ocean model.  This   !
!   version of the routine uses wind and ice-ocean stresses on a B-grid.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ocean_top_stress_Bgrid(IOF, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, part_size, G, IG)
  type(ice_ocean_flux_type), intent(inout) :: IOF
  type(SIS_hor_grid_type),   intent(inout) :: G
  type(ice_grid_type),       intent(inout) :: IG
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: windstr_x_water, str_ice_oce_x
  real, dimension(SZIB_(G),SZJB_(G)), intent(in) :: windstr_y_water, str_ice_oce_y
  real, dimension (SZI_(G),SZJ_(G),0:IG%CatIce), intent(in) :: part_size

  real    :: ps_vel ! part_size interpolated to a velocity point, nondim.
  integer :: i, j, k, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce


  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Bgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.
!$OMP parallel default(none) shared(isc,iec,jsc,jec,ncat,G,IOF,                &
!$OMP                               part_size,windstr_x_water,windstr_y_water, &
!$OMP                               str_ice_oce_x,str_ice_oce_y)               &
!$OMP                       private(ps_vel)
  if (IOF%flux_uv_stagger == AGRID) then
!$OMP do
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
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dBu(I,J)>0.5) ps_vel = &
                           0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                                 (part_size(i+1,j,0) + part_size(i,j+1,0)) )
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + windstr_x_water(I,J) * ps_vel
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + windstr_y_water(I,J) * ps_vel
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dBu(I,J)>0.5) then
        ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                         (part_size(i+1,j,k) + part_size(i,j+1,k)) )
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + str_ice_oce_x(I,J) * ps_vel
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + str_ice_oce_y(I,J) * ps_vel
      endif ; enddo ; enddo
    enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.5) ps_vel = &
                           0.5*(part_size(i+1,j,0) + part_size(i,j,0))
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * &
                0.5 * (windstr_x_water(I,J) + windstr_x_water(I,J-1))
        ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.5) ps_vel = &
                           0.5*(part_size(i,j+1,0) + part_size(i,j,0))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * &
                0.5 * (windstr_y_water(I,J) + windstr_y_water(I-1,J))
      enddo
      do k=1,ncat ; do i=isc,iec
        if (G%mask2dCu(I,j)>0.5) then
          ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
          IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * &
              0.5 * (str_ice_oce_x(I,J) + str_ice_oce_x(I,J-1))
        endif
        if (G%mask2dCv(i,J)>0.5) then
          ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
          IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * &
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
  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_Bgrid

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ocean_top_stress_Cgrid - Calculate the stresses on the ocean integrated  !
!   across all the thickness categories with the appropriate staggering, and   !
!   store them in the public ice data type for use by the ocean model.  This   !
!   version of the routine uses wind and ice-ocean stresses on a C-grid.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ocean_top_stress_Cgrid(IOF, windstr_x_water, windstr_y_water, &
                                      str_ice_oce_x, str_ice_oce_y, part_size, G, IG)
  type(ice_ocean_flux_type), intent(inout) :: IOF
  type(SIS_hor_grid_type),   intent(inout) :: G
  type(ice_grid_type),       intent(inout) :: IG
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: windstr_x_water, str_ice_oce_x
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: windstr_y_water, str_ice_oce_y
  real, dimension (SZI_(G),SZJ_(G),0:IG%CatIce), intent(in) :: part_size

  real    :: ps_vel ! part_size interpolated to a velocity point, nondim.
  integer :: i, j, k, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
  endif

  !   Copy and interpolate the ice-ocean stress_Cgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.
!$OMP parallel default(none) shared(isc,iec,jsc,jec,ncat,G,IOF,    &
!$OMP                               part_size,windstr_x_water,windstr_y_water, &
!$OMP                               str_ice_oce_x,str_ice_oce_y)               &
!$OMP                       private(ps_vel)
  if (IOF%flux_uv_stagger == AGRID) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = G%mask2dT(i,j) * part_size(i,j,0)
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * 0.5 * &
                            (windstr_x_water(I,j) + windstr_x_water(I-1,j))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * 0.5 * &
                            (windstr_y_water(I,j) + windstr_y_water(i,J-1))
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dT(i,j)>0.5) then
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) +  part_size(i,j,k) * 0.5 * &
                            (str_ice_oce_x(I,j) + str_ice_oce_x(I-1,j))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + part_size(i,j,k) * 0.5 * &
                            (str_ice_oce_y(I,j) + str_ice_oce_y(i,J-1))
      endif ; enddo ; enddo
    enddo
  elseif (IOF%flux_uv_stagger == BGRID_NE) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dBu(I,J)>0.5) ps_vel = &
                           0.25*((part_size(i+1,j+1,0) + part_size(i,j,0)) + &
                                 (part_size(i+1,j,0) + part_size(i,j+1,0)) )
        ! Consider deleting the masks here?
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
                (windstr_x_water(I,j) + windstr_x_water(I,j+1))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * G%mask2dBu(I,J) * 0.5 * &
                (windstr_y_water(I,j) + windstr_y_water(i+1,J))
      enddo
      do k=1,ncat ; do i=isc,iec ; if (G%mask2dBu(I,J)>0.5) then
        ps_vel = 0.25 * ((part_size(i+1,j+1,k) + part_size(i,j,k)) + &
                         (part_size(i+1,j,k) + part_size(i,j+1,k)) )
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * 0.5 * &
                            (str_ice_oce_x(I,j) + str_ice_oce_x(I,j+1))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * 0.5 * &
                            (str_ice_oce_y(I,j) + str_ice_oce_y(i+1,J))
      endif ; enddo ; enddo
    enddo
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
!$OMP do
    do j=jsc,jec
      do i=isc,iec
        ps_vel = 1.0 ; if (G%mask2dCu(I,j)>0.5) ps_vel = &
                           0.5*(part_size(i+1,j,0) + part_size(i,j,0))
        IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * windstr_x_water(I,j)
        ps_vel = 1.0 ; if (G%mask2dCv(i,J)>0.5) ps_vel = &
                           0.5*(part_size(i,j+1,0) + part_size(i,j,0))
        IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * windstr_y_water(i,J)
      enddo
      do k=1,ncat ; do i=isc,iec
        if (G%mask2dCu(I,j)>0.5) then
          ps_vel = 0.5 * (part_size(i+1,j,k) + part_size(i,j,k))
          IOF%flux_u_ocn(i,j) = IOF%flux_u_ocn(i,j) + ps_vel * str_ice_oce_x(I,j)
        endif
        if (G%mask2dCv(i,J)>0.5) then
          ps_vel = 0.5 * (part_size(i,j+1,k) + part_size(i,j,k))
          IOF%flux_v_ocn(i,j) = IOF%flux_v_ocn(i,j) + ps_vel * str_ice_oce_y(I,j)
        endif
      enddo ; enddo
    enddo
  else
!$OMP single
    call SIS_error(FATAL, "set_ocean_top_stress_Cgrid: Unrecognized flux_uv_stagger.")
!$OMP end single
  endif
!$OMP end parallel

  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_Cgrid

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_register_restarts allocates and registers any variables for this
!!      module that need to be included in the restart files.
subroutine SIS_dyn_trans_register_restarts(mpp_domain, HI, IG, param_file, CS, &
                                      Ice_restart, restart_file)
  type(domain2d),          intent(in) :: mpp_domain
  type(hor_index_type),    intent(in) :: HI
  type(ice_grid_type),     intent(in) :: IG     ! The sea-ice grid type
  type(param_file_type),   intent(in) :: param_file
  type(dyn_trans_CS),      pointer    :: CS
  type(restart_file_type), pointer    :: Ice_restart
  character(len=*),        intent(in) :: restart_file

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
!

!   This subroutine registers the restart variables associated with the
! the slow ice dynamics and thermodynamics.

  logical :: Cgrid_dyn
!  integer :: isd, ied, jsd, jed, id
!  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_dyn_trans_register_restarts called with an "//&
                            "associated control structure.")
    return
  endif
  allocate(CS)

  Cgrid_dyn = .false. ; call read_param(param_file, "CGRID_ICE_DYNAMICS", Cgrid_dyn)

  if (Cgrid_dyn) then
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
! SIS_dyn_trans_init - initializes ice model data, parameters and diagnostics       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_dyn_trans_init(Time, G, IG, param_file, diag, CS, output_dir, Time_init)
  type(time_type),     target, intent(in)    :: Time   ! current time
  type(SIS_hor_grid_type),     intent(in)    :: G      ! The horizontal grid structure
  type(ice_grid_type),         intent(in)    :: IG     ! The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(dyn_trans_CS),          pointer       :: CS
  character(len=*),            intent(in)    :: output_dir
  type(time_type),             intent(in)    :: Time_Init  ! starting time of model integration

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_dyn_trans" ! This module's name.
  real :: Time_unit      ! The time unit in seconds for ICE_STATS_INTERVAL.
  character(len=8) :: nstr
  integer :: n, nLay
  real, parameter    :: missing = -1e34

  nLay = IG%NkIce

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
  call log_version(param_file, mod, version, &
     "This module updates the ice momentum and does ice transport.")
  call get_param(param_file, mod, "SPECIFIED_ICE", CS%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  call get_param(param_file, mod, "CGRID_ICE_DYNAMICS", CS%Cgrid_dyn, &
                 "If true, use a C-grid discretization of the sea-ice \n"//&
                 "dynamics; if false use a B-grid discretization.", &
                 default=.false.)
  call get_param(param_file, mod, "DT_ICE_DYNAMICS", CS%dt_ice_dyn, &
                 "The time step used for the slow ice dynamics, including \n"//&
                 "stepping the continuity equation and interactions \n"//&
                 "between the ice mass field and velocities.  If 0 or \n"//&
                 "negative the coupling time step will be used.", &
                 units="seconds", default=-1.0)
  call get_param(param_file, mod, "DO_RIDGING", CS%do_ridging, &
                 "If true, apply a ridging scheme to the convergent ice. \n"//&
                 "Otherwise, ice is compressed proportionately if the \n"//&
                 "concentration exceeds 1.  The original SIS2 implementation \n"//&
                 "is based on work by Torge Martin.", default=.false.)

  call get_param(param_file, mod, "ICEBERG_WINDSTRESS_BUG", CS%berg_windstress_bug, &
                 "If true, use older code that applied an old ice-ocean \n"//&
                 "stress to the icebergs in place of the current air-ocean \n"//&
                 "stress.  This option is here for backward compatibility, \n"//&
                 "but should be avoided.", default=.false.)

  call get_param(param_file, mod, "TIMEUNIT", Time_unit, &
                 "The time unit for ICE_STATS_INTERVAL.", &
                 units="s", default=86400.0)
  call get_param(param_file, mod, "ICE_STATS_INTERVAL", CS%ice_stats_interval, &
                 "The interval in units of TIMEUNIT between writes of the \n"//&
                 "globally summed ice statistics and conservation checks.", &
                 default=set_time(0,1), timeunit=Time_unit)

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
  call get_param(param_file, mod, "VERBOSE", CS%verbose, &
                 "If true, write out verbose diagnostics.", default=.false.)

  if (CS%Cgrid_dyn) then
    call SIS_C_dyn_init(CS%Time, G, param_file, CS%diag, CS%SIS_C_dyn_CSp, CS%ntrunc)
  else
    call SIS_B_dyn_init(CS%Time, G, param_file, CS%diag, CS%SIS_B_dyn_CSp)
  endif
  call SIS_transport_init(CS%Time, G, param_file, CS%diag, CS%SIS_transport_CSp)

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

  ! Diagnostics that are specific to C-grid dynamics of the ice model
  if (CS%Cgrid_dyn) then
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
  endif
  CS%id_xprt = register_diag_field('ice_model','XPRT',diag%axesT1, Time, &
               'frozen water transport convergence', 'kg/(m^2*yr)', missing_value=missing)
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
  CS%id_rdgr    = register_diag_field('ice_model','RDG_RATE' ,diag%axesT1, Time, &
               'ice ridging rate', '1/sec', missing_value=missing)
!### THESE DIAGNOSTICS DO NOT EXIST YET.
!  CS%id_rdgf    = register_diag_field('ice_model','RDG_FRAC' ,diag%axesT1, Time, &
!               'ridged ice fraction', '0-1', missing_value=missing)
!  CS%id_rdgo    = register_diag_field('ice_model','RDG_OPEN' ,diag%axesT1, Time, &
!               'opening due to ridging', '1/s', missing_value=missing)
!  CS%id_rdgv    = register_diag_field('ice_model','RDG_VOSH' ,diag%axesT1, Time, &
!               'volume shifted from level to ridged ice', 'm^3/s', missing_value=missing)
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

subroutine safe_alloc_ids_1d(ids, nids)
  integer, allocatable :: ids(:)
  integer, intent(in)  :: nids

  if (.not.ALLOCATED(ids)) then
    allocate(ids(nids)) ; ids(:) = -1
  endif
end subroutine safe_alloc_ids_1d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_transport_CS returns a pointer to the SIS_transport_CS type that
!!  the dyn_trans_CS points to.
function SIS_dyn_trans_transport_CS(CS) result(transport_CSp)
  type(dyn_trans_CS), pointer :: CS
  type(SIS_transport_CS), pointer :: transport_CSp

  transport_CSp => CS%SIS_transport_CSp
end function SIS_dyn_trans_transport_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_transport_CS returns a pointer to the sum_out_CS type that
!!  the dyn_trans_CS points to.
function SIS_dyn_trans_sum_output_CS(CS) result(sum_out_CSp)
  type(dyn_trans_CS), pointer :: CS
  type(SIS_sum_out_CS), pointer :: sum_out_CSp

  sum_out_CSp => CS%sum_output_CSp
end function SIS_dyn_trans_sum_output_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_dyn_trans_end deallocates memory associated with the dyn_trans_CS type
!! and calls similar routines for subsidiary modules.
subroutine SIS_dyn_trans_end(CS)
  type(dyn_trans_CS), pointer :: CS

  if (CS%Cgrid_dyn) then
    call SIS_C_dyn_end(CS%SIS_C_dyn_CSp)
  else
    call SIS_B_dyn_end(CS%SIS_B_dyn_CSp)
  endif
  call SIS_transport_end(CS%SIS_transport_CSp)

  deallocate(CS)

end subroutine SIS_dyn_trans_end

end module SIS_dyn_trans
