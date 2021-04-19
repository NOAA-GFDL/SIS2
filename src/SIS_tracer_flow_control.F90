!> Contains subroutines into which calls to the tracer specific functions should be called.
module SIS_tracer_flow_control

! This file is a part of SIS2.  See LICENSE.md for the license.

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

use ice_grid,            only : ice_grid_type
use MOM_error_handler,   only : SIS_error=>MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_version, param_file_type
use SIS_diag_mediator,   only : time_type, SIS_diag_ctrl
use SIS_restart,         only : SIS_restart_CS
use SIS_hor_grid,        only : SIS_hor_grid_type
use SIS_tracer_registry, only : SIS_tracer_registry_type
use SIS_tracer_registry, only : register_SIS_tracer, register_SIS_tracer_pair

#include <SIS2_memory.h>

! Add references to other user-provided tracer modules here.
use ice_age_tracer, only : register_ice_age_tracer, initialize_ice_age_tracer
use ice_age_tracer, only : ice_age_tracer_column_physics
use ice_age_tracer, only : ice_age_stock, ice_age_end
use ice_age_tracer, only : ice_age_tracer_CS

implicit none ; private

public SIS_call_tracer_register, SIS_tracer_flow_control_init
public SIS_call_tracer_column_fns, SIS_call_tracer_stocks, SIS_tracer_flow_control_end

!> The control structure for orchestrating calls to the tracer packages
type, public :: SIS_tracer_flow_control_CS ; private
    logical :: use_ice_age = .false.  !< If true, use the ice age tracer packages.
    type(ice_age_tracer_CS), pointer :: ice_age_tracer_CSp => NULL()
                                      !< The control structure for the ice age tracer packages.
end type SIS_tracer_flow_control_CS

contains

! The following subroutines and associated definitions provide the
! machinery to register and call the subroutines that initialize
! tracers and apply vertical column processes to tracers.

!> Call the routines that register all of tracers in the tracer packages
subroutine SIS_call_tracer_register(G, IG, param_file, CS, diag, TrReg, Ice_restart)
  type(SIS_hor_grid_type),          intent(in) :: G   !< The horizontal grid type
  type(ice_grid_type),              intent(in) :: IG  !< The sea-ice specific grid type
  type(param_file_type),            intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS_tracer_flow_control_CS), pointer    :: CS  !< A pointer that is set to point to the
                                                      !! control structure for the tracer flow control
  type(SIS_diag_ctrl),              target     :: diag !< A structure that is used to regulate diagnostic output
  type(SIS_tracer_registry_type),   pointer    :: TrReg !< A pointer to the SIS tracer registry
  type(SIS_restart_CS),             pointer    :: Ice_restart !< The control structure for the ice restarts

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "SIS_tracer_flow_control" ! This module's name.

  if (associated(CS)) then
      call SIS_error(WARNING, "SIS_call_tracer_register called with an associated "// &
          "control structure.")
      return
  else ; allocate(CS) ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "USE_ICE_AGE_TRACER", CS%use_ice_age, &
      "If true, use the concentration based age tracer package.", &
      default=.false.)

  !    Add other user-provided calls to register tracers for restarting here. Each
  !  tracer package registration call returns a logical false if it cannot be run
  !  for some reason.  This then overrides the run-time selection from above.
  if (CS%use_ice_age) then
    CS%use_ice_age = register_ice_age_tracer(G, IG, param_file, CS%ice_age_tracer_CSp, &
                         diag, TrReg, Ice_restart)
  endif


end subroutine SIS_call_tracer_register

!> Call all registered tracer initialization subroutines.
subroutine SIS_tracer_flow_control_init(day, G, IG, param_file, CS, is_restart)
  type(time_type),          target, intent(in)    :: day !< The current model time
  type(SIS_hor_grid_type),          intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),              intent(in)    :: IG  !< The sea-ice specific grid type
  type(param_file_type),            intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_tracer_flow_control_CS), pointer       :: CS  !< The control structure returned by a
                                                         !! previous call to SIS_call_tracer_register.
  logical,                          intent(in)    :: is_restart !< A flag indicating whether this run
                                                         !! segment is being initialized from a restart file
  ! This subroutine calls all registered tracer initialization subroutines.

  if (.not. associated(CS)) call SIS_error(FATAL, "tracer_flow_control_init: "// &
      "Module must be initialized via call_tracer_register before it is used.")

  !  Add other user-provided calls here.
  if (CS%use_ice_age) &
      call initialize_ice_age_tracer(day, G, IG, CS%ice_age_tracer_CSp, is_restart)

end subroutine SIS_tracer_flow_control_init

!> Call all registered ice-tracer column physics subroutines
subroutine SIS_call_tracer_column_fns(dt, G, IG, CS, mi, mi_old)
  real,                    intent(in) :: dt  !< The amount of time covered by this call [T ~> s].
  type(SIS_hor_grid_type), intent(in) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG  !< The sea-ice specific grid type
  type(SIS_tracer_flow_control_CS), &
                           pointer    :: CS  !< The control structure returned by a
                                             !! previous call to SIS_call_tracer_register.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(in) :: mi  !< Mass of ice in a given category [R Z ~> kg m-2] at
                                             !! the end of the timestep
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(in) :: mi_old !< Mass of ice in a given category [R Z ~> kg m-2]
                                             !! at the beginning of the timestep

  ! This subroutine calls all registered ice-tracer column physics subroutines.

  if (.not. associated(CS)) call SIS_error(FATAL, "SIS_call_tracer_column_fns: "// &
      "Module must be initialized via call_tracer_register before it is used.")
  ! Add calls to tracer column functions here.
  if (CS%use_ice_age) &
      call ice_age_tracer_column_physics(dt, G, IG,  CS%ice_age_tracer_CSp, mi, mi_old)

end subroutine SIS_call_tracer_column_fns

!> Determine integrated stocks for the SIS tracer packages
subroutine SIS_call_tracer_stocks(G, IG, CS, mi, stock_values, stock_names, &
                                  stock_units, num_stocks)
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  type(ice_grid_type),              intent(in)  :: IG  !< The sea-ice specific grid type
  type(SIS_tracer_flow_control_CS), pointer     :: CS  !< The control structure returned by a
                                                       !! previous call to SIS_call_tracer_register.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                                    intent(in)  :: mi  !< Mass of ice in a given category used for
                                                       !! summing [R Z ~> kg m-2]
  real, dimension(:),               intent(out) :: stock_values !< The values of the summed tracer stocks.
  character(len=*), dimension(:), &
                         optional,  intent(out) :: stock_names !< The names of the summed tracer stocks.
  character(len=*), dimension(:), &
                         optional,  intent(out) :: stock_units !< The units of the tracer stocks
  integer,               optional,  intent(out) :: num_stocks  !< The number of summed tracer stocks.

  ! Local variables
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
    ns = ice_age_stock(mi, values, G, IG, CS%ice_age_tracer_CSp, names, units)
    call SIS_store_stocks("ice_age_tracer", ns, names, units, values, stock_values, &
        ns_tot, stock_names, stock_units)

  endif
  num_stocks = ns_tot

end subroutine SIS_call_tracer_stocks

!> Store the stocks for for a single tracer package for SIS_call_tracer_stocks.
subroutine SIS_store_stocks(pkg_name, ns, names, units, values, stock_values, &
                        ns_tot, stock_names, stock_units)
  character(len=*),   intent(in)    :: pkg_name !< The name of this tracer package
  integer,            intent(in)    :: ns    !< The number of stocks with this tracer package
  character(len=*), dimension(:), &
                      intent(in)    :: names !< The tracer names in this package
  character(len=*), dimension(:), &
                      intent(in)    :: units !< The tracer units
  real, dimension(:), intent(in)    :: values !< The values of the tracer stocks
  real, dimension(:), intent(inout) :: stock_values !< The stored tracer stock values
  integer,            intent(inout) :: ns_tot !< The total number of tracer stocks across all packages
  character(len=*), dimension(:), &
            optional, intent(inout) :: stock_names !< The tracer names for use in the stock reporting
  character(len=*), dimension(:), &
            optional, intent(inout) :: stock_units !< The units for use in stock reporting

! This routine stores the stocks for SIS_call_tracer_stocks.
  integer :: n

  do n=1,ns
    stock_values(ns_tot+n) = values(n)
    if (present(stock_names)) stock_names(ns_tot+n) = names(n)
    if (present(stock_units)) stock_units(ns_tot+n) = units(n)
  enddo
  ns_tot = ns_tot + ns

end subroutine SIS_store_stocks

!> Call all of the tracer package end routines and deallocate tracer memory
subroutine SIS_tracer_flow_control_end(CS)
  type(SIS_tracer_flow_control_CS), pointer :: CS  !< The control structure returned by a previous
                                                   !! call to SIS_call_tracer_register that is
                                                   !! deallocated here

  if (CS%use_ice_age) call ice_age_end(CS%ice_age_tracer_CSp)

  if (associated(CS)) deallocate(CS)
end subroutine SIS_tracer_flow_control_end

end module SIS_tracer_flow_control
