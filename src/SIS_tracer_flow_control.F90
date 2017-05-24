module SIS_tracer_flow_control
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

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Andrew Shao, April 2016                                         *
!*      Adapted from MOM6 code MOM_tracer_flow_control.F90             *
!*                                                                     *
!*    This module contains subroutines into which calls to the tracer  *
!*    specific functions should be called. To add a new tracer the     *
!*    calls should be added to each of the following subroutines.      *
!*    Use ice_age_tracer.F90 as a model on which to base new types of  *
!*    ice tracers. Generally, most tracer packages will want to        *
!*    define their own tracer control structures which will get        *
!*    pointed to by the main tracer registry in                        *
!*    SIS_tracer_registry.F90.                                         *
!*                                                                     *
!*    SIS_call_tracer_register:                                        *
!*      Allocates the tracer flow control structure and reads the      *
!*      parameter file SIS_input to identify which tracers will be     *
!*      included in the current run. Tracer packages should allocate   *
!*      their tracer arrays here, establish whether the tracer will    *
!*      be initialized from a restart file, and register the tracer    *
!*      for the advection module.                                      *
!*                                                                     *
!*    SIS_tracer_flow_control_init:                                    *
!*      Sets the initial conditions for the tracer. Also, register     *
!*      tracer fields for output with the diag_manager.                *
!*                                                                     *
!*    SIS_call_tracer_column_fns:                                      *
!*      Apply any special treatment of the tracer that may occur in    *
!*      a vertical column. For example, source and sink terms should   *
!*      be calculated and applied here                                 *
!*                                                                     *
!*    SIS_call_tracer_stocks:                                          *
!*      NOTE: CURRENTLY NOT CALLED. This subroutine is intended to     *
!*      at some point calculate the inventory of each tracer and       *
!*      ready it for either output to a file or STDOUT                 *
!*                                                                     *
!*    SIS_tracer_flow_control_end:                                     *
!*      Deallocate arrays and prepare for model finalization           *
!*                                                                     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use SIS_diag_mediator, only : time_type, SIS_diag_ctrl
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING
use MOM_coms, only : sum_across_PEs
use ice_grid, only : ice_grid_type
use SIS_tracer_registry, only : SIS_tracer_registry_type
use SIS_tracer_registry, only : register_SIS_tracer, register_SIS_tracer_pair
use SIS_hor_grid, only : SIS_hor_grid_type

use fms_io_mod,      only : restart_file_type
use MOM_file_parser, only : get_param, log_version, param_file_type

#include <SIS2_memory.h>

! Add references to other user-provided tracer modules here.
use ice_age_tracer, only : register_ice_age_tracer, initialize_ice_age_tracer
use ice_age_tracer, only : ice_age_tracer_column_physics
use ice_age_tracer, only : ice_age_stock, ice_age_end
use ice_age_tracer, only : ice_age_tracer_CS

implicit none ; private

public SIS_call_tracer_register, SIS_tracer_flow_control_init
public SIS_call_tracer_column_fns, SIS_call_tracer_stocks, SIS_tracer_flow_control_end

type, public :: SIS_tracer_flow_control_CS ; private
    logical :: use_ice_age = .false.
    type(ice_age_tracer_CS), pointer :: ice_age_tracer_CSp => NULL()
end type SIS_tracer_flow_control_CS

contains

! The following subroutines and associated definitions provide the
! machinery to register and call the subroutines that initialize
! tracers and apply vertical column processes to tracers.

subroutine SIS_call_tracer_register(G, IG, param_file, CS, diag, TrReg, &
    Ice_restart, restart_file)
  type(SIS_hor_grid_type),                intent(in) :: G
  type(ice_grid_type),                    intent(in) :: IG
  type(param_file_type),                  intent(in) :: param_file
  type(SIS_tracer_flow_control_CS),       pointer    :: CS
  type(SIS_diag_ctrl),                    target     :: diag
  type(SIS_tracer_registry_type),         pointer    :: TrReg
  type(restart_file_type),                intent(inout) :: Ice_restart
  character(len=*),                       intent(in) :: restart_file

  ! Argument:  G - The ice model's horizontal grid structure.
  !  (in)      IG - The ice model's grid structure.
  !  (in)      param_file - A structure indicating the open file to parse for
  !                         model parameter values.
  !
  !  (in/out)  CS - A pointer that is set to point to the control structure
  !                 for the tracer flow control
  !  (in)      diag - A structure that is used to regulate diagnostic output.
  !  (in/out)  TrReg - A pointer that is set to point to the control structure
  !                  for the tracer advection and diffusion module.
  !  (in/out)  Ice model restart file to be written to
  !  (in)      Path to the restart file
  !
  ! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_tracer_flow_control" ! This module's name.

  if (associated(CS)) then
      call SIS_error(WARNING, "SIS_call_tracer_register called with an associated "// &
          "control structure.")
      return
  else ; allocate(CS) ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "USE_ICE_AGE_TRACER", CS%use_ice_age, &
      "If true, use the concentration based age tracer package.", &
      default=.false.)

  !    Add other user-provided calls to register tracers for restarting here. Each
  !  tracer package registration call returns a logical false if it cannot be run
  !  for some reason.  This then overrides the run-time selection from above.
  if (CS%use_ice_age) then
    CS%use_ice_age = register_ice_age_tracer(G, IG, param_file, CS%ice_age_tracer_CSp, &
        diag, TrReg, Ice_restart, restart_file)
  endif


end subroutine SIS_call_tracer_register

subroutine SIS_tracer_flow_control_init(day, G, IG, param_file, CS, is_restart)
  type(time_type), target,                    intent(in) :: day
  type(SIS_hor_grid_type),                    intent(inout) :: G
  type(ice_grid_type),                        intent(in) :: IG
  type(param_file_type),                      intent(in) :: param_file
  type(SIS_tracer_flow_control_CS),           pointer    :: CS
  logical,                                    intent(in) :: is_restart
  !   This subroutine calls all registered tracer initialization
  ! subroutines.

  ! Arguments:
  !  (in)      day - Time of the start of the run.
  !  (in)      G - The ice model's horizontal grid structure.
  !  (in)      IG - The ice model's vertical grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 call_tracer_register.
  !  (in)      is_restart - flag for whether tracer should be initialized from restart
  if (.not. associated(CS)) call SIS_error(FATAL, "tracer_flow_control_init: "// &
      "Module must be initialized via call_tracer_register before it is used.")

  !  Add other user-provided calls here.
  if (CS%use_ice_age) &
      call initialize_ice_age_tracer(day, G, IG, CS%ice_age_tracer_CSp, is_restart)

end subroutine SIS_tracer_flow_control_init

subroutine SIS_call_tracer_column_fns(dt, G, IG, CS, mi, mi_old)

  real,                                           intent(in) :: dt
  type(SIS_hor_grid_type),                        intent(in) :: G
  type(ice_grid_type),                            intent(in) :: IG
  type(SIS_tracer_flow_control_CS),               pointer    :: CS
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),    intent(in) :: mi
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),    intent(in) :: mi_old
  !   This subroutine calls all registered tracer column physics
  ! subroutines.

  ! Arguments:
  !  (in)      dt - The amount of time covered by this call, in s.
  !  (in)      G - The ice model's grid structure.
  !  (in)      IG - The ice model's vertical grid structure.

  if (.not. associated(CS)) call SIS_error(FATAL, "SIS_call_tracer_column_fns: "// &
      "Module must be initialized via call_tracer_register before it is used.")
  ! Add calls to tracer column functions here.
  if (CS%use_ice_age) &
      call ice_age_tracer_column_physics(dt, G, IG,  CS%ice_age_tracer_CSp, mi, mi_old)

end subroutine SIS_call_tracer_column_fns

subroutine SIS_call_tracer_stocks(G, IG, CS, mi, stock_values, stock_names, &
                                  stock_units, num_stocks)
  type(SIS_hor_grid_type),                        intent(in) :: G
  type(ice_grid_type),                            intent(in) :: IG
  type(SIS_tracer_flow_control_CS),               pointer    :: CS
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),    intent(in) :: mi
  real, dimension(:),                             intent(  out) :: stock_values
  character(len=*), dimension(:),      optional,  intent(  out) :: stock_names
  character(len=*), dimension(:),      optional,  intent(  out) :: stock_units
  integer,                             optional,  intent(  out) :: num_stocks

  ! Arguments:
  !  (in)      G - The ocean's grid structure.
  !  (in)      IG - The ocean's vertical grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 call_tracer_register.
  !  (in)      mi - mass of ice in a given category, used for summing
  !  (out)     stocks - Global integral of tracer
  !  (out)     nstocks - Number of passive tracer stocks

  character(len=200), dimension(MAX_FIELDS_) :: names, units
  character(len=200) :: set_pkg_name
  real, dimension(MAX_FIELDS_) :: values
  integer :: ns_tot, ns, index, m

  if (.not. associated(CS)) call SIS_error(FATAL, "SIS_call_tracer_stocks: "// &
      "Module must be initialized via call_tracer_register before it is used.")

  ns_tot = 0
  values(:) = 0.0
  stock_values(:) = 0.0

  !  Add other user-provided calls here.
  if (CS%use_ice_age) then
    ns = ice_age_stock(mi, values, G, IG, CS%ice_age_tracer_CSp, &
                         names, units)
    call SIS_store_stocks("ice_age_tracer", ns, names, units, values, stock_values, &
        ns_tot, stock_names, stock_units)

  endif
  num_stocks = ns_tot

end subroutine SIS_call_tracer_stocks

subroutine SIS_store_stocks(pkg_name, ns, names, units, values, stock_values, &
                        ns_tot, stock_names, stock_units)
  character(len=*),                         intent(in)    :: pkg_name
  integer,                                  intent(in)    :: ns
  character(len=*), dimension(:),           intent(in)    :: names, units
  real, dimension(:),                       intent(in)    :: values
  real, dimension(:),                       intent(inout) :: stock_values
  integer,                                  intent(inout) :: ns_tot
  character(len=*), dimension(:), optional, intent(inout) :: stock_names, stock_units

! This routine stores the stocks for SIS_call_tracer_stocks.
  integer :: n

  do n=1,ns
    stock_values(ns_tot+n) = values(n)
    if (present(stock_names)) stock_names(ns_tot+n) = names(n)
    if (present(stock_units)) stock_units(ns_tot+n) = units(n)
  enddo
  ns_tot = ns_tot + ns

end subroutine SIS_store_stocks

subroutine SIS_tracer_flow_control_end(CS)
  type(SIS_tracer_flow_control_CS), pointer :: CS

  if (CS%use_ice_age) call ice_age_end(CS%ice_age_tracer_CSp)

  if (associated(CS)) deallocate(CS)
end subroutine SIS_tracer_flow_control_end

end module SIS_tracer_flow_control
