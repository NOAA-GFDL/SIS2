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
module SIS_slow_mod

use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use MOM_checksums,     only :  chksum, Bchksum, hchksum, uchksum, vchksum
use SIS_error_checking, only : check_redundant_B, check_redundant_C
use SIS_sum_output, only : write_ice_statistics, SIS_sum_output_init, SIS_sum_out_CS

use mpp_domains_mod,  only  : domain2D !, mpp_get_compute_domain, CORNER, EAST, NORTH
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_domains,       only : fill_symmetric_edges !, MOM_domains_init, clone_MOM_domain
! use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, log_version, param_file_type
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

use ice_type_mod, only : ice_state_type, ice_ocean_flux_type, fast_ice_avg_type
use ice_type_mod, only : ocean_sfc_state_type
use ice_type_mod, only : dyn_trans_CS
! use ice_type_mod, only : dealloc_IST_arrays, ice_state_register_restarts
! use ice_type_mod, only : ice_diagnostics_init
use ice_type_mod, only : IST_chksum,  IST_bounds_check
use ice_utils_mod, only : get_avg, post_avg, ice_line !, ice_grid_chksum
use SIS_hor_grid, only : SIS_hor_grid_type

use ice_grid, only : ice_grid_type

! use SIS_tracer_registry, only : register_SIS_tracer, register_SIS_tracer_pair
! use SIS_tracer_flow_control, only : SIS_call_tracer_register, SIS_tracer_flow_control_init
! use SIS_tracer_flow_control, only : SIS_call_tracer_column_fns
! use SIS_tracer_flow_control, only : SIS_tracer_flow_control_end

use SIS2_ice_thm,  only: get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,  only: enth_from_TS, Temp_from_En_S, T_freeze
use SIS_dyn_bgrid, only: SIS_B_dynamics, SIS_B_dyn_init, SIS_B_dyn_register_restarts, SIS_B_dyn_end
use SIS_dyn_cgrid, only: SIS_C_dynamics, SIS_C_dyn_init, SIS_C_dyn_register_restarts, SIS_C_dyn_end
use ice_transport_mod, only : ice_transport, ice_transport_init, ice_transport_end
use ice_transport_mod, only : ice_transport_CS

use ice_bergs,        only: icebergs, icebergs_run, icebergs_init, icebergs_end, icebergs_incr_mass

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_dynamics_trans, update_icebergs
public :: SIS_slow_register_restarts, SIS_slow_init, SIS_slow_end
public :: SIS_slow_transport_CS, SIS_slow_sum_output_CS

integer :: iceClock4, iceClock8, iceClock9, iceClocka, iceClockb, iceClockc

contains

subroutine update_icebergs(IST, OSS, IOF, FIA, icebergs_CS, G, IG)
  type(ice_state_type),       intent(inout) :: IST
  type(ocean_sfc_state_type), intent(in)    :: OSS
  type(fast_ice_avg_type),    intent(in)    :: FIA
  type(ice_ocean_flux_type),  intent(inout) :: IOF
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(inout) :: IG
  type(icebergs),             pointer       :: icebergs_CS

  real, dimension(SZI_(G),SZJ_(G))   :: &
    hi_avg            ! The area-weighted average ice thickness, in m.
  real :: H_to_m_ice  ! The specific volume of ice times the conversion factor
                      ! from thickness units, in m H-1.

  integer :: isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  H_to_m_ice = IG%H_to_kg_m2 / IST%Rho_ice
  call get_avg(IST%mH_ice, IST%part_size(:,:,1:), hi_avg, wtd=.true.)
  hi_avg(:,:) = hi_avg(:,:) * H_to_m_Ice

  !### I think that there is long-standing bug here, in that the old ice-ocean
  !###  stresses are being passed in place of the wind stresses on the icebergs. -RWH
  if (IST%Cgrid_dyn) then
    call icebergs_run( icebergs_CS, IST%Time, &
            IOF%calving(isc:iec,jsc:jec), OSS%u_ocn_C(isc-2:iec+1,jsc-1:jec+1), &
            OSS%v_ocn_C(isc-1:iec+1,jsc-2:jec+1), IST%u_ice_C(isc-2:iec+1,jsc-1:jec+1), &
            IST%v_ice_C(isc-1:iec+1,jsc-2:jec+1), &
            IOF%flux_u_ocn(isc:iec,jsc:jec), IOF%flux_v_ocn(isc:iec,jsc:jec), &
            OSS%sea_lev(isc-1:iec+1,jsc-1:jec+1), IST%t_surf(isc:iec,jsc:jec,0),  &
            IOF%calving_hflx(isc:iec,jsc:jec), FIA%ice_cover(isc-1:iec+1,jsc-1:jec+1), &
            hi_avg(isc-1:iec+1,jsc-1:jec+1), stagger=CGRID_NE, &
            stress_stagger=IOF%flux_uv_stagger)
  else
    call icebergs_run( icebergs_CS, IST%Time, &
            IOF%calving(isc:iec,jsc:jec), OSS%u_ocn_B(isc-1:iec+1,jsc-1:jec+1), &
            OSS%v_ocn_B(isc-1:iec+1,jsc-1:jec+1), IST%u_ice_B(isc-1:iec+1,jsc-1:jec+1), &
            IST%v_ice_B(isc-1:iec+1,jsc-1:jec+1), &
            IOF%flux_u_ocn(isc:iec,jsc:jec), IOF%flux_v_ocn(isc:iec,jsc:jec), &
            OSS%sea_lev(isc-1:iec+1,jsc-1:jec+1), IST%t_surf(isc:iec,jsc:jec,0),  &
            IOF%calving_hflx(isc:iec,jsc:jec), FIA%ice_cover(isc-1:iec+1,jsc-1:jec+1), &
            hi_avg(isc-1:iec+1,jsc-1:jec+1), stagger=BGRID_NE, &
            stress_stagger=IOF%flux_uv_stagger)
  endif

end subroutine update_icebergs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_dynamics_trans - do ice dynamics and mass and tracer transport           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_dynamics_trans(IST, OSS, FIA, IOF, dt_slow, CS, icebergs_CS, G, IG)

  type(ice_state_type),       intent(inout) :: IST
  type(ocean_sfc_state_type), intent(in)    :: OSS
  type(fast_ice_avg_type),    intent(inout) :: FIA
  type(ice_ocean_flux_type),  intent(inout) :: IOF
  real,                       intent(in)    :: dt_slow
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(inout) :: IG
  type(dyn_trans_CS),         pointer       :: CS
  type(icebergs),             pointer       :: icebergs_CS

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

  real, dimension(SZI_(G),SZJ_(G))   :: diagVar ! An temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBx ! An temporary array for diagnostics.
  real, dimension(SZIB_(G),SZJB_(G)) :: diagVarBy ! An temporary array for diagnostics.

  real :: weights  ! A sum of the weights around a point.
  real :: I_wts    ! 1.0 / wts or 0 if wts is 0, nondim.
  real :: ps_vel   ! The fractional thickness catetory coverage at a velocity point.

  real :: dt_slow_dyn
  integer :: ndyn_steps
  real :: Idt_slow
  real :: I_Nk        ! The inverse of the number of layers in the ice.
  integer :: i, j, k, l, m, isc, iec, jsc, jec, ncat, NkIce, nds
  integer :: isd, ied, jsd, jed
  integer ::iyr, imon, iday, ihr, imin, isec

  real, dimension(IG%NkIce) :: S_col ! Specified thermodynamic salinity of each
                                     ! ice layer if spec_thermo_sal is true.
  logical :: spec_thermo_sal
  logical :: do_temp_diags
  real :: enth_units, I_enth_units
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    rdg_frac    ! fraction of ridged ice per category
  real, dimension(SZI_(G),SZJ_(G)) :: &
    rdg_open, & ! formation rate of open water due to ridging
    rdg_vosh, & ! rate of ice volume shifted from level to ridged ice
!!   rdg_s2o, &  ! snow volume [m] dumped into ocean during ridging
    rdg_rate, & ! Niki: Where should this come from?
    snow2ocn
  real    :: tmp3  ! This is a bad name - make it more descriptive!

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; NkIce = IG%NkIce
  I_Nk = 1.0 / NkIce
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  if (IST%specified_ice) then
    ndyn_steps = 0.0 ; dt_slow_dyn = 0.0
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
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,WindStr_x_A,WindStr_y_A,  &
!$OMP                                  ice_cover,ice_free,WindStr_x_ocn_A,       &
!$OMP                                  WindStr_y_ocn_A)
  do j=jsd,jed
    do i=isd,ied
      WindStr_x_ocn_A(i,j) = FIA%flux_u_top(i,j,0)
      WindStr_y_ocn_A(i,j) = FIA%flux_v_top(i,j,0)

      ice_cover(i,j) = FIA%ice_cover(i,j) ; ice_free(i,j) = FIA%ice_free(i,j)
      WindStr_x_A(i,j) = FIA%WindStr_x(i,j) ; WindStr_y_A(i,j) = FIA%WindStr_y(i,j)
    enddo
  enddo

  do nds=1,ndyn_steps

    call enable_SIS_averaging(dt_slow_dyn, CS%Time - set_time(int((ndyn_steps-nds)*dt_slow_dyn)), CS%diag)

    ! Correct the wind stresses for changes in the fractional ice-coverage.
    ice_cover(:,:) = 0.0
!$OMP parallel do default(none) shared(isd,ied,jsd,jed,ncat,ice_cover,IST,FIA,ice_free, &
!$OMP                                  WindStr_x_A,WindStr_y_A,WindStr_x_ocn_A,WindStr_y_ocn_A)
    do j=jsd,jed
      do k=1,ncat ; do i=isd,ied
        ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
      enddo ; enddo
      do i=isd,ied
        ice_free(i,j) = IST%part_size(i,j,0)

        if (ice_cover(i,j) > FIA%ice_cover(i,j)) then
          WindStr_x_A(i,j) = ((ice_cover(i,j)-FIA%ice_cover(i,j))*FIA%flux_u_top(i,j,0) + &
                              FIA%ice_cover(i,j)*FIA%WindStr_x(i,j)) / ice_cover(i,j)
          WindStr_y_A(i,j) = ((ice_cover(i,j)-FIA%ice_cover(i,j))*FIA%flux_v_top(i,j,0) + &
                              FIA%ice_cover(i,j)*FIA%WindStr_y(i,j)) / ice_cover(i,j)
!        else
!          WindStr_x_A(i,j) = FIA%WindStr_x(i,j)
!          WindStr_y_A(i,j) = FIA%WindStr_y(i,j)
!        endif
!        if (ice_free(i,j) <= FIA%ice_free(i,j)) then
!          WindStr_x_ocn_A(i,j) = FIA%flux_u_top(i,j,0)
!          WindStr_y_ocn_A(i,j) = FIA%flux_v_top(i,j,0)
!        else
        elseif (ice_free(i,j) > FIA%ice_free(i,j)) then
          WindStr_x_ocn_A(i,j) = ((ice_free(i,j)-FIA%ice_free(i,j))*FIA%WindStr_x(i,j) + &
                              FIA%ice_free(i,j)*FIA%flux_u_top(i,j,0)) / ice_free(i,j)
          WindStr_y_ocn_A(i,j) = ((ice_free(i,j)-FIA%ice_free(i,j))*FIA%WindStr_y(i,j) + &
                              FIA%ice_free(i,j)*FIA%flux_v_top(i,j,0)) / ice_free(i,j)
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
        call uchksum(OSS%u_ocn_C, "u_ocn_C before SIS_C_dynamics", G%HI)
        call vchksum(OSS%v_ocn_C, "v_ocn_C before SIS_C_dynamics", G%HI)
        call uchksum(WindStr_x_Cu, "WindStr_x_Cu before SIS_C_dynamics", G%HI)
        call vchksum(WindStr_y_Cv, "WindStr_y_Cv before SIS_C_dynamics", G%HI)
        call check_redundant_C("WindStr before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, G)
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
        call Bchksum(OSS%u_ocn_B, "u_ocn before ice_dynamics", G%HI, symmetric=.true.)
        call Bchksum(OSS%v_ocn_B, "v_ocn before ice_dynamics", G%HI, symmetric=.true.)
        call Bchksum(WindStr_x_B, "WindStr_x_B before ice_dynamics", G%HI, symmetric=.true.)
        call Bchksum(WindStr_y_B, "WindStr_y_B before ice_dynamics", G%HI, symmetric=.true.)
        call check_redundant_B("WindStr before ice_dynamics",WindStr_x_B, WindStr_y_B, G)
      endif

      rdg_rate(:,:) = 0.0
      call mpp_clock_begin(iceClocka)
      call SIS_B_dynamics(1.0-IST%part_size(:,:,0), ms_sum, mi_sum, IST%u_ice_B, IST%v_ice_B, &
                        OSS%u_ocn_B, OSS%v_ocn_B, WindStr_x_B, WindStr_y_B, OSS%sea_lev, &
                        str_x_ice_ocn_B, str_y_ice_ocn_B, IST%do_ridging, &
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
      call ice_transport(IST%part_size, IST%mH_ice, IST%mH_snow, IST%u_ice_C, IST%v_ice_C, &
                         IST%TrReg, OSS%sea_lev, dt_slow_dyn, G, IG, CS%ice_transport_CSp, &
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
                         IST%TrReg, OSS%sea_lev, dt_slow_dyn, G, IG, CS%ice_transport_CSp, &
                         IST%rdg_mice, snow2ocn, rdg_rate, &
                         rdg_open, rdg_vosh)
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

  ! Add snow volume dumped into ocean to flux of frozen precipitation:
  !### WARNING - rdg_s2o is never calculated!!!
!  if (CS%do_ridging) then ; do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
!    FIA%fprec_top(i,j,k) = FIA%fprec_top(i,j,k) + rdg_s2o(i,j)*(IST%Rho_snow/dt_slow)
!  enddo ; enddo ; enddo ; endif

  call mpp_clock_begin(iceClock8)

  call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)

  ! Set appropriate surface quantities in categories with no ice.  Change <1e-10 to == 0?
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)<1e-10) &
    IST%t_surf(i,j,k) = T_0degC + T_Freeze(OSS%s_surf(i,j),IST%ITV)
  enddo ; enddo ; enddo

  if (CS%bounds_check) call IST_bounds_check(IST, G, IG, "After ice_transport", OSS=OSS)
  if (CS%debug) call IST_chksum("After ice_transport", IST, G, IG)

  ! Sum the concentration weighted mass for diagnostics.
  if (CS%id_mi>0 .or. CS%id_mib>0) then
    mass_ice(:,:) = 0.0
    mass_snow(:,:) = 0.0
    mass(:,:) = 0.0
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,mass,G,IST,IG)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      mass_ice(i,j) = mass_ice(i,j) + IG%H_to_kg_m2*IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
      mass_snow(i,j) = mass_snow(i,j) + IG%H_to_kg_m2*IST%mH_snow(i,j,k)*IST%part_size(i,j,k)
      mass(i,j) = mass_ice(i,j) + mass_snow(i,j)
    enddo ; enddo ; enddo
    if (CS%id_simass>0) call post_data(CS%id_simass, mass_ice(isc:iec,jsc:jec), CS%diag)
    if (CS%id_sisnmass>0) call post_data(CS%id_sisnmass, mass_snow(isc:iec,jsc:jec), CS%diag)
    if (CS%id_mi>0) call post_data(CS%id_mi, mass(isc:iec,jsc:jec), CS%diag)

    if (CS%id_mib>0) then
      if (IST%do_icebergs) call icebergs_incr_mass(icebergs_CS, mass(isc:iec,jsc:jec)) ! Add icebergs mass in kg/m^2
      call post_data(CS%id_mib, mass(isc:iec,jsc:jec), CS%diag)
    endif
  endif

  call mpp_clock_end(iceClock8)

  !
  ! Thermodynamic state diagnostics
  !
  call mpp_clock_begin(iceClock9)
  if (IST%id_cn>0) call post_data(IST%id_cn, IST%part_size(:,:,1:ncat), IST%diag)
  if (IST%id_siconc>0) call post_data(IST%id_siconc, sum(IST%part_size(:,:,1:ncat),3), IST%diag)
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
    call post_data(IST%id_ext, diagVar, CS%diag)
  endif
  if (IST%id_sitimefrac>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (IST%part_size(i,j,0) < 1.0) diagVar(i,j) = 1.0
    enddo ; enddo
    call post_data(IST%id_sitimefrac, diagVar, CS%diag)
  endif
  if (IST%id_sisnconc>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec; do k=1,ncat
      if (IST%part_size(i,j,k) > 0.0 .and. IST%mH_snow(i,j,k) > 0.0) then
        diagVar(i,j) = diagVar(i,j) + IST%part_size(i,j,k)
      endif
    enddo ; enddo ; enddo
    call post_data(IST%id_sisnconc, diagVar, CS%diag)
  endif
  if (IST%id_hs>0) call post_avg(IST%id_hs, IST%mH_snow, IST%part_size(:,:,1:), &
                                 CS%diag, G=G, &
                                 scale=IG%H_to_kg_m2/IST%Rho_snow, wtd=.true.)
  if (IST%id_sisnthick>0) call post_avg(IST%id_sisnthick, IST%mH_snow, IST%part_size(:,:,1:), &
                                 CS%diag, G=G, &
                                 scale=IG%H_to_kg_m2/IST%Rho_snow, wtd=.true.)
  if (IST%id_hi>0) call post_avg(IST%id_hi, IST%mH_ice, IST%part_size(:,:,1:), &
                                 CS%diag, G=G, &
                                 scale=IG%H_to_kg_m2/IST%Rho_ice, wtd=.true.)
  if (IST%id_sithick>0) call post_avg(IST%id_sithick, IST%mH_ice, IST%part_size(:,:,1:), &
                                 CS%diag, G=G, &
                                 scale=IG%H_to_kg_m2/IST%Rho_ice, wtd=.true.)
  if (IST%id_sivol>0) call post_avg(IST%id_sivol, IST%mH_ice, IST%part_size(:,:,1:), &
                                 CS%diag, G=G, &
                                 scale=IG%H_to_kg_m2/IST%Rho_ice, wtd=.true.)
  if (IST%id_tsfc>0) call post_avg(IST%id_tsfc, IST%t_surf(:,:,1:), IST%part_size(:,:,1:), &
                                 CS%diag, G=G, offset=-T_0degC, wtd=.true.)
  if (IST%id_sitemptop>0) call post_avg(IST%id_sitemptop, IST%t_surf(:,:,1:), IST%part_size(:,:,1:), &
                                 CS%diag, G=G, offset=-T_0degC, wtd=.true.)
  if (IST%id_tsn>0) call post_avg(IST%id_tsn, temp_snow, IST%part_size(:,:,1:), &
                                 CS%diag, G=G, wtd=.true.)
  do m=1,NkIce
    if (IST%id_t(m)>0) call post_avg(IST%id_t(m), temp_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   CS%diag, G=G, wtd=.true.)
    if (IST%id_sal(m)>0) call post_avg(IST%id_sal(m), IST%sal_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   CS%diag, G=G, wtd=.true.)
  enddo
  if (IST%id_t_iceav>0) call post_avg(IST%id_t_iceav, temp_ice, IST%part_size(:,:,1:), &
                                    CS%diag, G=G, wtd=.true.)
  if (IST%id_S_iceav>0) call post_avg(IST%id_S_iceav, IST%sal_ice, IST%part_size(:,:,1:), &
                                    CS%diag, G=G, wtd=.true.)


  if (CS%id_xprt>0) then
    call post_data(CS%id_xprt, h2o_chg_xprt(isc:iec,jsc:jec)*864e2*365/dt_slow, &
                   CS%diag)
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
    call post_data(IST%id_e2m,  tmp2d(:,:), CS%diag)
  endif
  call disable_SIS_averaging(CS%diag)

  !
  ! Ridging diagnostics
  !
  !TOM> preparing output field fraction of ridged ice rdg_frac = (ridged ice volume) / (total ice volume)
  !     in each category; IST%rdg_mice is ridged ice mass per unit total
  !     area throughout the code.
  if (IST%do_ridging) then
    call enable_SIS_averaging(dt_slow, CS%Time, CS%diag)

!     if (IST%id_rdgf>0) then
! !$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,G,rdg_frac,IG) &
! !$OMP                          private(tmp3)
!       do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
!         tmp3 = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
!         if (tmp3*IG%H_to_kg_m2 > IST%Rho_Ice*1.e-5) then   ! 1 mm ice thickness x 1% ice concentration
!           rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp3
!         else
!           rdg_frac(i,j,k) = 0.0
!         endif
!       enddo ; enddo ; enddo
!       call post_data(IST%id_rdgf, rdg_frac(isc:iec,jsc:jec), CS%diag)
!     endif

    if (IST%id_rdgr>0) call post_data(IST%id_rdgr, rdg_rate(isc:iec,jsc:jec), CS%diag)
!    if (IST%id_rdgo>0) call post_data(IST%id_rdgo, rdg_open(isc:iec,jsc:jec), CS%diag)
!    if (IST%id_rdgv>0) then
!      do j=jsc,jec ; do i=isc,iec
!        tmp2d(i,j) = rdg_vosh(i,j) * G%areaT(i,j) * G%mask2dT(i,j)
!      enddo ; enddo
!      call post_data(IST%id_rdgv, tmp2d, CS%diag)
!    endif
  endif

  if (CS%verbose) then
    call get_date(CS%Time, iyr, imon, iday, ihr, imin, isec)
    call get_time(CS%Time-set_date(iyr,1,1,0,0,0),isec,iday)
    call ice_line(iyr, iday+1, isec, IST%part_size(isc:iec,jsc:jec,0), &
                              IST%t_surf(:,:,0)-T_0degC, G)
  endif

  call mpp_clock_end(iceClock9)

  if (CS%debug) then
    call IST_chksum("End SIS_dynamics_trans", IST, G, IG)
  endif

  if (CS%bounds_check) then
    call IST_bounds_check(IST, G, IG, "End of SIS_dynamics_trans", OSS=OSS)
  endif

  if (CS%Time + (IST%Time_step_slow/2) > CS%write_ice_stats_time) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp)
    CS%write_ice_stats_time = CS%write_ice_stats_time + CS%ice_stats_interval
  elseif (CS%column_check) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, IG, CS%sum_output_CSp)
  endif

end subroutine SIS_dynamics_trans

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
!$OMP parallel default(none) shared(isc,iec,jsc,jec,ncat,G,                    &
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
!$OMP parallel default(none) shared(isc,iec,jsc,jec,ncat,G,    &
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
!> SIS_slow_register_restarts allocates and registers any variables for this 
!!      module that need to be included in the restart files.
subroutine SIS_slow_register_restarts(mpp_domain, HI, IG, param_file, CS, &
                                      Ice_restart, restart_file)
  type(domain2d),          intent(in)    :: mpp_domain
  type(hor_index_type),    intent(in)    :: HI
  type(ice_grid_type),     intent(in)    :: IG     ! The sea-ice grid type
  type(param_file_type),   intent(in)    :: param_file
  type(dyn_trans_CS),      pointer       :: CS
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

  logical :: Cgrid_dyn
!  integer :: isd, ied, jsd, jed, id
!  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_slow_register_restarts called with an "//&
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
!  call ice_transport_register_restarts(G, param_file, CS%ice_transport_CSp, &
!                                       Ice_restart, restart_file)

end subroutine SIS_slow_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_slow_init - initializes ice model data, parameters and diagnostics       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_slow_init(Time, G, IG, param_file, diag, CS, output_dir, Time_init)
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
  character(len=40)  :: mod = "SIS_slow" ! This module's name.
  real :: Time_unit      ! The time unit in seconds for ICE_STATS_INTERVAL.
  real, parameter    :: missing = -1e34

  call callTree_enter("SIS_slow_init(), SIS_slow.F90")

  if (.not.associated(CS)) call SIS_error(FATAL, &
      "SIS_slow_init called with an unassociated control structure.")
  if (CS%module_is_initialized) then
    call SIS_error(WARNING, "SIS_slow_init called with a control "// &
                            "structure that has already been initialized.")
    return
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, &
     "This module updates the ice momentum and does ice transport.")
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
  call ice_transport_init(CS%Time, G, param_file, CS%diag, CS%ice_transport_CSp)

  call SIS_sum_output_init(G, param_file, output_dir, Time_Init, &
                           CS%sum_output_CSp, CS%ntrunc)

  CS%write_ice_stats_time = Time_Init + CS%ice_stats_interval * &
      (1 + (Time - Time_init) / CS%ice_stats_interval)

  !
  ! diagnostics that are specific to C-grid dynamics of the ice model
  !
  if (CS%Cgrid_dyn) then
    CS%id_fax = register_diag_field('ice_model', 'FA_X', diag%axesCu1, Time, &
               'Air stress on ice on C-grid - x component', 'Pa', missing_value=missing)
    CS%id_fay = register_diag_field('ice_model', 'FA_Y', diag%axesCv1, Time, &
               'Air stress on ice on C-grid - y component', 'Pa', missing_value=missing)
  else
    CS%id_fax = register_diag_field('ice_model', 'FA_X', diag%axesB1, Time, &
               'air stress on ice - x component', 'Pa', missing_value=missing)
    CS%id_fay = register_diag_field('ice_model', 'FA_Y', diag%axesB1, Time, &
               'air stress on ice - y component', 'Pa', missing_value=missing)
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

  iceClock4 = mpp_clock_id( '  Ice: slow: dynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClocka = mpp_clock_id( '       slow: ice_dynamics', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClockb = mpp_clock_id( '       slow: comm/cut check ', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClockc = mpp_clock_id( '       slow: diags', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock8 = mpp_clock_id( '  Ice: slow: transport', flags=clock_flag_default, grain=CLOCK_LOOP )
  iceClock9 = mpp_clock_id( '  Ice: slow: thermodyn diags', flags=clock_flag_default, grain=CLOCK_LOOP )

  call callTree_leave("SIS_slow_init()")

end subroutine SIS_slow_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_slow_transport_CS returns a pointer to the ice_transport_CS type that
!!  the dyn_trans_CS points to.
function SIS_slow_transport_CS(CS) result(transport_CSp)
  type(dyn_trans_CS), pointer :: CS
  type(ice_transport_CS), pointer :: transport_CSp

  transport_CSp => CS%ice_transport_CSp
end function SIS_slow_transport_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_slow_transport_CS returns a pointer to the ice_transport_CS type that
!!  the dyn_trans_CS points to.
function SIS_slow_sum_output_CS(CS) result(sum_out_CSp)
  type(dyn_trans_CS), pointer :: CS
  type(SIS_sum_out_CS), pointer :: sum_out_CSp

  sum_out_CSp => CS%sum_output_CSp
end function SIS_slow_sum_output_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_slow_end - deallocates memory                                            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_slow_end(CS)
  type(dyn_trans_CS), pointer :: CS

  if (CS%Cgrid_dyn) then
    call SIS_C_dyn_end(CS%SIS_C_dyn_CSp)
  else
    call SIS_B_dyn_end(CS%SIS_B_dyn_CSp)
  endif
  call ice_transport_end(CS%ice_transport_CSp)

end subroutine SIS_slow_end

end module SIS_slow_mod
