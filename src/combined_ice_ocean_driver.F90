module combined_ice_ocean_driver
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM and SIS2.                                *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!-----------------------------------------------------------------------
!
! This module provides a common interface for jointly stepping SIS2 and
! MOM6, and will evolve as a platform for tightly integrating the ocean
! and sea ice models.
!
! <CONTACT EMAIL="Robert.Hallberg@noaa.gov"> Robert Hallberg
! </CONTACT>
!

use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave
! use MOM_file_parser, only : get_param, log_version, close_param_file, param_file_type
use MOM_time_manager, only : time_type, get_time !, set_time, operator(>)

use ice_model_mod,   only: ice_data_type, ice_model_end
use ice_model_mod,   only: update_ice_model_slow

use ocean_model_mod, only: update_ocean_model,  ocean_model_end! , ocean_model_init
use ocean_model_mod, only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
use flux_exchange_mod,  only: flux_ice_to_ocean

implicit none ; private

public update_slow_ice_and_ocean, ice_ocean_driver_init, ice_ocean_driver_end

type, public :: ice_ocean_driver_type ; private
  logical :: CS_is_initialized = .false.
end type ice_ocean_driver_type

contains

!=======================================================================
! <SUBROUTINE NAME="ice_ocean_driver_init">
!
! <DESCRIPTION>
! Initialize the coupling between the slow-ice and ocean models.
! </DESCRIPTION>
!
!>   This subroutine initializes the combined ice ocean coupling control type.
subroutine ice_ocean_driver_init(CS, Time_init, Time_in)
  type(ice_ocean_driver_type), pointer       :: CS        !< The control structure for combined ice-ocean driver
  type(time_type),             intent(in)    :: Time_init !< The start time for the coupled model's calendar.
  type(time_type),             intent(in)    :: Time_in   !< The time at which to initialize the coupled model.

!     real :: Time_unit   ! The time unit in seconds for ENERGYSAVEDAYS.
!   ! This include declares and sets the variable "version".
!   #include "version_variable.h"
!     character(len=40)  :: mod = "ice_ocean_driver_init"  ! This module's name.
!     character(len=48)  :: stagger
!     integer :: secs, days
!  type(param_file_type) :: param_file !< A structure to parse for run-time parameters

  call callTree_enter("ice_ocean_driver_init(), combined_ice_ocean_driver.F90")
  if (associated(CS)) then
    call MOM_error(WARNING, "ice_ocean_driver_init called with an associated "// &
                    "ice_ocean_driver_CS structure. Model is already initialized.")
    return
  endif
  allocate(CS)

!     OS%is_ocean_pe = Ocean_sfc%is_ocean_pe
!     if (.not.OS%is_ocean_pe) return

!     ! Read all relevant parameters and write them to the model log.
!     call log_version(param_file, mod, version, "")
!     call get_param(param_file, mod, "RESTART_CONTROL", OS%Restart_control, &
!                    "An integer whose bits encode which restart files are \n"//&
!                    "written. Add 2 (bit 1) for a time-stamped file, and odd \n"//&
!                    "(bit 0) for a non-time-stamped file.  A restart file \n"//&
!                    "will be saved at the end of the run segment for any \n"//&
!                    "non-negative value.", default=1)
!     call get_param(param_file, mod, "TIMEUNIT", Time_unit, &
!                    "The time unit for ENERGYSAVEDAYS.", &
!                    units="s", default=86400.0)

!     call get_param(param_file, mod, "OCEAN_SURFACE_STAGGER", stagger, &
!                    "A case-insensitive character string to indicate the \n"//&
!                    "staggering of the surface velocity field that is \n"//&
!                    "returned to the coupler.  Valid values include \n"//&
!                    "'A', 'B', or 'C'.", default="C")
!     if (uppercase(stagger(1:1)) == 'A') then ; Ocean_sfc%stagger = AGRID
!     elseif (uppercase(stagger(1:1)) == 'B') then ; Ocean_sfc%stagger = BGRID_NE
!     elseif (uppercase(stagger(1:1)) == 'C') then ; Ocean_sfc%stagger = CGRID_NE
!     else ; call MOM_error(FATAL,"ice_ocean_driver_init: OCEAN_SURFACE_STAGGER = "// &
!                           trim(stagger)//" is invalid.") ; endif

!     call close_param_file(param_file)
  CS%CS_is_initialized = .true.

  call callTree_leave("ice_ocean_driver_init(")
end subroutine ice_ocean_driver_init
! </SUBROUTINE> NAME="ice_ocean_driver_init"


!=======================================================================
! <SUBROUTINE NAME="update_slow_ice_and_ocean">
!
! <DESCRIPTION>
! Advance the slow portions of the sea-ice and the ocean in tandem.
! </DESCRIPTION>
!

!>   The subroutine update_slow_ice_and_ocean uses the forcing already stored in
!! the ice_data_type to advance both the sea-ice (and icebergs) and ocean states
!! for a time interval coupling_time_step.
subroutine update_slow_ice_and_ocean(CS, Ice, Ocn, Ocean_sfc, Ice_ocean_boundary, &
                                     time_start_update, coupling_time_step)
  type(ice_ocean_driver_type), &
                           pointer       :: CS   !< The control structure for this driver 
  type(ice_data_type),     intent(inout) :: Ice  !< The publicly visible ice data type
  type(ocean_state_type),  pointer       :: Ocn  !< The internal ocean state and control structures
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< The publicly visible ocean surface state type
  type(ice_ocean_boundary_type), &
                           intent(inout) :: Ice_ocean_boundary !< A structure containing the various forcing
                                                               !! fields going from the ice to the ocean
                                                               !! The arrays of this type are intent out.
  type(time_type),         intent(in)    :: time_start_update  !< The time at the beginning of the update step
  type(time_type),         intent(in)    :: coupling_time_step !< The amount of time over which to advance
                                                               !! the ocean and ice

  real :: time_step         ! The time step of a call to step_MOM in seconds.
  integer :: secs, days

  call callTree_enter("update_ice_and_ocean(), combined_ice_ocean_driver.F90")
  call get_time(coupling_time_step, secs, days)
  time_step = 86400.0*real(days) + real(secs)

!  if (time_start_update /= CS%Time) then
!    call MOM_error(WARNING, "update_ice_and_ocean: internal clock does not "//&
!                            "agree with time_start_update argument.")
!  endif
  if (.not.associated(CS)) then
    call MOM_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ice_ocean_driver_type. ice_ocean_driver_init must be "//  &
                    "called first to allocate this structure.")
    return
  endif
  if (.not.associated(Ocn)) then
    call MOM_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ocean_state_type structure. ocean_model_init must be "//  &
                    "called first to allocate this structure.")
    return
  endif

  ! Throw an error if the ice and ocean do not use the same grid and
  ! domain decomposition.

  ! Add clocks after there is an init call.

  call update_ice_model_slow(Ice)

  call flux_ice_to_ocean( time_start_update, Ice, Ocean_sfc, Ice_ocean_boundary )

  call update_ocean_model( Ice_ocean_boundary, Ocn, Ocean_sfc, &
                           time_start_update, coupling_time_step )

  call callTree_leave("update_ice_and_ocean()")
end subroutine update_slow_ice_and_ocean
! </SUBROUTINE> NAME="update_slow_ice_and_ocean"


!=======================================================================
! <SUBROUTINE NAME="ice_ocean_driver_end">
!
! <DESCRIPTION>
! Close down the sea-ice and ocean models
! </DESCRIPTION>
!
!>   The subroutine ice_ocean_driver_end terminates the model run, saving
!! the ocean and slow ice states in restart files and deallocating any data
!! associated with the ocean and slow ice.
subroutine ice_ocean_driver_end(CS, Ice, Ocean_sfc, Ocn, Time)
  type(ice_ocean_driver_type), pointer   :: CS   !< The control structure for combined ice-ocean driver
  type(ice_data_type),     intent(inout) :: Ice  !< The publicly visible ice data type
  type(ocean_state_type),  pointer       :: Ocn  !< The internal ocean state and control structures
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< The publicly visible ocean surface state type
  type(time_type),         intent(in)    :: Time !< The model time, used for writing restarts

  call ice_model_end (Ice)

  call ocean_model_end(Ocean_sfc, Ocn, Time)

  if (associated(CS)) deallocate(CS)

end subroutine ice_ocean_driver_end
! </SUBROUTINE> NAME="ice_ocean_driver_end"

end module combined_ice_ocean_driver
