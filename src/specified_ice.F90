!> Handles the stresses for specified ice.
module specified_ice

! This file is part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   This module handles the specified ice using SIS2 types and interfaces.     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_domains,       only : AGRID, BGRID_NE, CGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,  only : get_param, read_param, log_param, log_version, param_file_type
use MOM_time_manager, only : time_type, time_type_to_real, real_to_time
use MOM_time_manager, only : operator(+), operator(-)
use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)
use MOM_unit_scaling, only : unit_scale_type

use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_hor_grid,  only : SIS_hor_grid_type
use SIS_ice_diags, only : ice_state_diags_type, register_ice_state_diagnostics
use SIS_ice_diags, only : post_ice_state_diagnostics
use SIS_sum_output, only : write_ice_statistics, SIS_sum_output_init, SIS_sum_out_CS
use SIS_types,     only : ocean_sfc_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_types,     only : ice_state_type, IST_chksum, IST_bounds_check
use ice_grid,      only : ice_grid_type

implicit none ; private

#include <SIS2_memory.h>

public :: specified_ice_dynamics, specified_ice_init, specified_ice_end, specified_ice_sum_output_CS

!> The control structure for the specified_ice module
type, public :: specified_ice_CS ; private
  logical :: debug        !< If true, write verbose checksums for debugging purposes.
  logical :: bounds_check !< If true, check for sensible values of thicknesses
                          !! temperatures, fluxes, etc.
  integer :: ntrunc = 0   !< The number of times the velocity has been truncated
                          !! since the last call to write_ice_statistics.
  integer :: n_calls = 0  !< The number of times specified_ice_dynamics has been called.
  type(time_type) :: ice_stats_interval !< The interval between writes of the
                          !! globally summed ice statistics and conservation checks.
  type(time_type) :: write_ice_stats_time !< The next time to write out the ice statistics.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  type(ice_state_diags_type), pointer :: IDs => NULL()
      !< A structure for regulating sea ice state diagnostics
  type(SIS_sum_out_CS), pointer   :: sum_output_CSp => NULL()
     !< Pointer to the control structure for the summed diagnostics module
end type specified_ice_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> specified_ice_dynamics does an update of ice dynamic quantities with specified ice.
subroutine specified_ice_dynamics(IST, OSS, FIA, IOF, dt_slow, CS, G, US, IG)

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
  type(specified_ice_CS),     pointer       :: CS  !< The control structure for the specified_ice module

  ! Local variables
  integer :: i, j, k, isc, iec, jsc, jec, ncat

  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  CS%n_calls = CS%n_calls + 1

  IOF%stress_count = 0
  call set_ocean_top_stress_FIA(FIA, IOF, G, US)

  ! Set appropriate surface quantities in categories with no ice.
  if (allocated(IST%t_surf)) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)<=0.0) &
      IST%t_surf(i,j,k) = T_0degC + OSS%T_fr_ocn(i,j)
    enddo ; enddo ; enddo
  endif

  call enable_SIS_averaging(US%T_to_s*dt_slow, CS%Time, CS%diag)
  call post_ice_state_diagnostics(CS%IDs, IST, OSS, IOF, dt_slow, CS%Time, G, US, IG, CS%diag)
  call disable_SIS_averaging(CS%diag)

  if (CS%debug) call IST_chksum("End specified_ice_dynamics", IST, G, US, IG)
  if (CS%bounds_check) call IST_bounds_check(IST, G, US, IG, "End of specified_ice_dynamics", OSS=OSS)

  if (CS%Time + real_to_time(0.5*US%T_to_s*dt_slow) > CS%write_ice_stats_time) then
    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp)
    CS%write_ice_stats_time = CS%write_ice_stats_time + CS%ice_stats_interval
  endif

end subroutine specified_ice_dynamics


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Calculate the stresses on the ocean integrated across all the thickness categories
!! with the appropriate staggering, based on the information in a fast_ice_avg_type.
subroutine set_ocean_top_stress_FIA(FIA, IOF, G, US)
  type(fast_ice_avg_type),   intent(inout) :: FIA !< A type containing averages of fields
                                                  !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type), intent(inout) :: IOF !< A structure containing fluxes from the ice to
                                                  !! the ocean that are calculated by the ice model.
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors

  real :: ps_ice, ps_ocn  ! ice_free and ice_cover interpolated to a velocity point [nondim].
  real :: wt_prev, wt_now ! Relative weights of the previous average and the current step [nondim].
  real :: taux2, tauy2    ! Squared wind stresses [R2 Z2 L2 T-4 ~> Pa2]
  integer :: i, j, k, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (IOF%stress_count == 0) then
    IOF%flux_u_ocn(:,:) = 0.0 ; IOF%flux_v_ocn(:,:) = 0.0
    if (allocated(IOF%stress_mag)) IOF%stress_mag(:,:) = 0.0
  endif

  wt_now = 1.0 / (real(IOF%stress_count) + 1.0) ; wt_prev = 1.0 - wt_now

  !   Copy and interpolate the ice-ocean stress_Cgrid.  This code is slightly
  ! complicated because there are 3 different staggering options supported.

  if (IOF%flux_uv_stagger == AGRID) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do i=isc,iec
      ps_ocn = G%mask2dT(i,j) * FIA%ice_free(i,j)
      ps_ice = G%mask2dT(i,j) * FIA%ice_cover(i,j)
      IOF%flux_u_ocn(i,j) = wt_prev * IOF%flux_u_ocn(i,j) + wt_now * &
           (ps_ocn * FIA%WindStr_ocn_x(i,j) + ps_ice * FIA%WindStr_x(i,j))
      IOF%flux_v_ocn(i,j) = wt_prev * IOF%flux_v_ocn(i,j) + wt_now * &
           (ps_ocn * FIA%WindStr_ocn_y(i,j) + ps_ice * FIA%WindStr_y(i,j))
      if (allocated(IOF%stress_mag)) &
        IOF%stress_mag(i,j) = wt_prev * IOF%stress_mag(i,j) + wt_now * &
           sqrt(IOF%flux_u_ocn(i,j)**2 + IOF%flux_v_ocn(i,j)**2)
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
      IOF%flux_u_ocn(I,J) = wt_prev * IOF%flux_u_ocn(I,J) + wt_now * &
          (ps_ocn * 0.25 * ((FIA%WindStr_ocn_x(i,j) + FIA%WindStr_ocn_x(i+1,j+1)) + &
                            (FIA%WindStr_ocn_x(i,j+1) + FIA%WindStr_ocn_x(i+1,j))) + &
           ps_ice * 0.25 * ((FIA%WindStr_x(i,j) + FIA%WindStr_x(i+1,j+1)) + &
                            (FIA%WindStr_x(i,j+1) + FIA%WindStr_x(i+1,J))) )
      IOF%flux_v_ocn(I,J) = wt_prev * IOF%flux_v_ocn(I,J) + wt_now * &
          (ps_ocn * 0.25 * ((FIA%WindStr_ocn_y(i,j) + FIA%WindStr_ocn_y(i+1,j+1)) + &
                            (FIA%WindStr_ocn_y(i,j+1) + FIA%WindStr_ocn_y(i+1,j))) + &
           ps_ice * 0.25 * ((FIA%WindStr_y(i,j) + FIA%WindStr_y(i+1,j+1)) + &
                            (FIA%WindStr_y(i,j+1) + FIA%WindStr_y(i+1,J))) )
    enddo ; enddo
    if (allocated(IOF%stress_mag)) then ; do j=jsc,jec ; do i=isc,iec
      if (((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + &
           (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) > 0) then
        IOF%stress_mag(i,j) = wt_prev * IOF%stress_mag(i,j) + wt_now * sqrt( &
          ((G%mask2dBu(I,J)*(IOF%flux_u_ocn(I,J)**2 + IOF%flux_v_ocn(I,J)**2) + &
            G%mask2dBu(I-1,J-1)*(IOF%flux_u_ocn(I-1,J-1)**2 + IOF%flux_v_ocn(I-1,J-1)**2)) + &
           (G%mask2dBu(I,J-1)*(IOF%flux_u_ocn(I,J-1)**2 + IOF%flux_v_ocn(I,J-1)**2) + &
            G%mask2dBu(I-1,J)*(IOF%flux_u_ocn(I-1,J)**2 + IOF%flux_v_ocn(I-1,J)**2)) ) / &
          ((G%mask2dBu(I,J) + G%mask2dBu(I-1,J-1)) + (G%mask2dBu(I,J-1) + G%mask2dBu(I-1,J))) )
      else
        IOF%stress_mag(i,j) = 0.0
      endif
    enddo ; enddo ; endif
  elseif (IOF%flux_uv_stagger == CGRID_NE) then
    !$OMP parallel do default(shared) private(ps_ocn, ps_ice)
    do j=jsc,jec ; do I=Isc-1,iec
      ps_ocn = 1.0 ; ps_ice = 0.0
      if (G%mask2dCu(I,j)>0.5) then
        ps_ocn = 0.5*(FIA%ice_free(i+1,j) + FIA%ice_free(i,j))
        ps_ice = 0.5*(FIA%ice_cover(i+1,j) + FIA%ice_cover(i,j))
      endif
      IOF%flux_u_ocn(I,j) = wt_prev * IOF%flux_u_ocn(I,j) + wt_now * &
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
      IOF%flux_v_ocn(i,J) = wt_prev * IOF%flux_v_ocn(i,J) + wt_now * &
          (ps_ocn * 0.5 * (FIA%WindStr_ocn_y(i,j+1) + FIA%WindStr_ocn_y(i,j)) + &
           ps_ice * 0.5 * (FIA%WindStr_y(i,j+1) + FIA%WindStr_y(i,j)) )
    enddo ; enddo
    if (allocated(IOF%stress_mag)) then ; do j=jsc,jec ; do i=isc,iec
      taux2 = 0.0 ; tauy2 = 0.0
      if ((G%mask2dCu(I-1,j) + G%mask2dCu(I,j)) > 0) &
        taux2 = (G%mask2dCu(I-1,j)*IOF%flux_u_ocn(I-1,j)**2 + &
                 G%mask2dCu(I,j)*IOF%flux_u_ocn(I,j)**2) / (G%mask2dCu(I-1,j) + G%mask2dCu(I,j))
      if ((G%mask2dCv(i,J-1) + G%mask2dCv(i,J)) > 0) &
        tauy2 = (G%mask2dCv(i,J-1)*IOF%flux_v_ocn(i,J-1)**2 + &
                 G%mask2dCv(i,J)*IOF%flux_v_ocn(i,J)**2) / (G%mask2dCv(i,J-1) + G%mask2dCv(i,J))
      IOF%stress_mag(i,j) = wt_prev * IOF%stress_mag(i,j) + wt_now * sqrt(taux2 + tauy2)
    enddo ; enddo ; endif
  else
    call SIS_error(FATAL, "set_ocean_top_stress_FIA: Unrecognized flux_uv_stagger.")
  endif

  IOF%stress_count = IOF%stress_count + 1

end subroutine set_ocean_top_stress_FIA

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> specified_ice_init initializes ice model data, parameters and diagnostics
!!   associated with the SIS2 dynamics and transport modules.
subroutine specified_ice_init(Time, G, IG, param_file, diag, CS, output_dir, Time_init)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model time.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid structure
  type(ice_grid_type),         intent(in)    :: IG   !< The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(specified_ice_CS),      pointer       :: CS   !< The control structure for the specified_ice module
  character(len=*),            intent(in)    :: output_dir !< The directory to use for writing output
  type(time_type),             intent(in)    :: Time_Init !< Starting time of the model integration

  ! This include declares and sets the variable "version".
#  include "version_variable.h"
  character(len=40) :: mdl = "specified_ice" ! This module's name.
  real :: Time_unit      ! The time unit for ICE_STATS_INTERVAL [s].
  logical :: debug

  call callTree_enter("specified_ice_init(), specified_ice.F90")

  if (associated(CS)) then
    call SIS_error(WARNING, "specified_ice_init called with an "//&
                            "associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, &
     "This module updates the ice momentum and does ice transport.")

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
  call get_param(param_file, mdl, "ICE_BOUNDS_CHECK", CS%bounds_check, &
                 "If true, periodically check the values of ice and snow "//&
                 "temperatures and thicknesses to ensure that they are "//&
                 "sensible, and issue warnings if they are not.  This "//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)

  call SIS_sum_output_init(G, param_file, output_dir, Time_Init, G%US, &
                           CS%sum_output_CSp, CS%ntrunc)

  CS%write_ice_stats_time = Time_Init + CS%ice_stats_interval * &
      (1 + (Time - Time_init) / CS%ice_stats_interval)

  call register_ice_state_diagnostics(Time, IG, G%US, param_file, diag, CS%IDs)

  call callTree_leave("specified_ice_init()")

end subroutine specified_ice_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> specified_ice_sum_output_CS returns a pointer to the sum_out_CS type that
!! the specified_ice_CS points to.
function specified_ice_sum_output_CS(CS) result(sum_out_CSp)
  type(specified_ice_CS),   pointer :: CS    !< The control structure for the specified_ice module
  type(SIS_sum_out_CS), pointer :: sum_out_CSp !< The SIS_sum_out_CS type used by specified_ice

  sum_out_CSp => CS%sum_output_CSp
end function specified_ice_sum_output_CS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> specified_ice_end deallocates memory associated with the specified_ice_CS type.
subroutine specified_ice_end(CS)
  type(specified_ice_CS), pointer :: CS  !< The control structure for the specified_ice module that
                                         !! is dellocated here

  deallocate(CS)

end subroutine specified_ice_end

end module specified_ice
