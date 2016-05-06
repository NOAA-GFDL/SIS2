module SIS_tracer_flow_control
    !***********************************************************************
    !*                   GNU General Public License                        *
    !* This file is a part of MOM.                                         *
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

    !********+*********+*********+*********+*********+*********+*********+**
    !*                                                                     *
    !*  By Andrew Shao, April 2016                                         *
    !*      Adapter from MOM6 code MOM_tracer_flow_control.F90             *
    !*                                                                     *
    !*    This module contains two subroutines into which calls to other   *
    !*  tracer initialization (call_tracer_init_fns) and column physics    *
    !*  routines (call_tracer_column_fns) can be inserted.                 *
    !*                                                                     *
    !********+*********+*********+*********+*********+*********+*********+**

    use SIS_diag_mediator, only : time_type, SIS_diag_ctrl
    use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING
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
    use ice_age_tracer_type, only : ice_age_tracer_CS

    implicit none ; private

    public SIS_call_tracer_register, SIS_tracer_flow_control_init
    public SIS_call_tracer_column_fns, SIS_call_tracer_stocks, SIS_tracer_flow_control_end

    type, public :: SIS_tracer_flow_control_CS ; private
        logical :: use_ice_age = .false.
        logical :: conc_age = .false.
        logical :: mass_age = .false.
        type(ice_age_tracer_CS), pointer :: ice_age_tracer_CSp => NULL()
    end type SIS_tracer_flow_control_CS

contains

    ! The following subroutines and associated definitions provide the
    ! machinery to register and call the subroutines that initialize
    ! tracers and apply vertical column processes to tracers.

    subroutine SIS_call_tracer_register(G, IG, param_file, CS, diag, TrReg, &
        Ice_restart, restart_file)
        type(SIS_hor_grid_type),        intent(in) :: G
        type(ice_grid_type),            intent(in) :: IG
        type(param_file_type),          intent(in) :: param_file
        type(SIS_tracer_flow_control_CS), pointer, intent(inout) :: CS
        type(SIS_diag_ctrl), target,      intent(in) :: diag
        type(SIS_tracer_registry_type),   pointer, intent(inout) :: TrReg
        type(restart_file_type),        intent(inout) :: Ice_restart
        character(len=*),               intent(in) :: restart_file

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
            if (CS%use_ice_age) CS%use_ice_age = &
                register_ice_age_tracer(G, IG, param_file, CS%ice_age_tracer_CSp, &
                diag, TrReg, Ice_restart, restart_file)


        end subroutine SIS_call_tracer_register

        subroutine SIS_tracer_flow_control_init(restart, day, G, IG, param_file, CS)
            logical,                               intent(in) :: restart
            type(time_type), target,               intent(in) :: day
            type(SIS_hor_grid_type),               intent(inout) :: G
            type(ice_grid_type),                   intent(in) :: IG
            type(param_file_type),                 intent(in) :: param_file
            type(SIS_tracer_flow_control_CS),      pointer    :: CS
            !   This subroutine calls all registered tracer initialization
            ! subroutines.

            ! Arguments: restart - 1 if the fields have already been read from
            !                     a restart file.
            !  (in)      day - Time of the start of the run.
            !  (in)      G - The ice model's horizontal grid structure.
            !  (in)      IG - The ice model's vertical grid structure.
            !  (in)      CS - The control structure returned by a previous call to
            !                 call_tracer_register.
            if (.not. associated(CS)) call SIS_error(FATAL, "tracer_flow_control_init: "// &
                "Module must be initialized via call_tracer_register before it is used.")

            !  Add other user-provided calls here.
!            if (CS%use_ice_age) &
!                call initialize_ice_age_tracer(restart, day, G, IG, CS%ice_age_tracer_CSp)
        end subroutine SIS_tracer_flow_control_init

        subroutine SIS_call_tracer_column_fns(dt, G, IG, CS)
            real,                                       intent(in) :: dt
            type(SIS_hor_grid_type),                    intent(inout) :: G
            type(ice_grid_type),                        intent(in) :: IG
            type(SIS_tracer_flow_control_CS), pointer,  intent(in) :: CS
            !   This subroutine calls all registered tracer column physics
            ! subroutines.

            ! Arguments:
            !  (in)      dt - The amount of time covered by this call, in s.
            !  (in)      G - The ice model's grid structure.
            !  (in)      IG - The ice model's vertical grid structure.

            if (.not. associated(CS)) call SIS_error(FATAL, "SIS_call_tracer_column_fns: "// &
                "Module must be initialized via call_tracer_register before it is used.")
            ! Add calls to tracer column functions here.
!             if (CS%use_ice_age) &
!                call ice_age_tracer_column_physics(dt, G, IG, &
!                    SIS_tracer_flow_CSp%ice_age_tracer_CSp)

        end subroutine SIS_call_tracer_column_fns

        subroutine SIS_call_tracer_stocks(stock_values, G, IG, CS, stock_names, stock_units, &
            num_stocks, stock_index, got_min_max, global_min, global_max,xgmin, ygmin, zgmin, xgmax, ygmax, zgmax)

            type(SIS_hor_grid_type),                  intent(in)  :: G
            real, dimension(:),                       intent(out) :: stock_values
            type(ice_grid_type),                      intent(in)  :: IG
            type(SIS_tracer_flow_control_CS),         pointer     :: CS
            character(len=*), dimension(:), optional, intent(out) :: stock_names
            character(len=*), dimension(:), optional, intent(out) :: stock_units
            integer,                        optional, intent(out) :: num_stocks
            integer,                        optional, intent(in)  :: stock_index
            logical,  dimension(:),         optional, intent(inout) :: got_min_max
            real, dimension(:),             optional, intent(out) :: global_min,  global_max
            real, dimension(:),             optional, intent(out) :: xgmin, ygmin, zgmin, xgmax, ygmax, zgmax
            !   This subroutine calls all registered tracer packages to enable them to
            ! add to the surface state returned to the coupler. These routines are optional.

            ! Arguments:
            !  (out)     stock_values - The integrated amounts of a tracer on the current
            !                           PE, usually in kg x concentration.
            !  (in)      G - The ocean's grid structure.
            !  (in)      IG - The ocean's vertical grid structure.
            !  (in)      CS - The control structure returned by a previous call to
            !                 call_tracer_register.
            !  (out,opt) stock_names - Diagnostic names to use for each stock.
            !  (out,opt) stock_units - Units to use in the metadata for each stock.
            !  (out,opt) num_stocks - The number of tracer stocks being returned.
            !  (in,opt)  stock_index - The integer stock index from stocks_constans_mod of
            !                          the stock to be returned.  If this is present and
            !                          greater than 0, only a single stock can be returned.
            character(len=200), dimension(MAX_FIELDS_) :: names, units
            character(len=200) :: set_pkg_name
            real, dimension(MAX_FIELDS_) :: values
            integer :: max_ns, ns_tot, ns, index, pkg, max_pkgs, nn

            if (.not. associated(CS)) call SIS_error(FATAL, "SIS_call_tracer_stocks: "// &
                "Module must be initialized via call_tracer_register before it is used.")

            index = -1 ; if (present(stock_index)) index = stock_index
            ns_tot = 0
            max_ns = size(stock_values)
            if (present(stock_names)) max_ns = min(max_ns,size(stock_names))
            if (present(stock_units)) max_ns = min(max_ns,size(stock_units))

            !  Add other user-provided calls here.
!            if (CS%use_ice_age) then
!                ns = ice_age_stock(h, values, G, IG, CS%ice_age_tracer_CSp, &
!                    names, units, stock_index)
!                call store_stocks("ice_age", ns, names, units, values, index, &
!                    stock_values, set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
!            endif

            if (ns_tot == 0) stock_values(1) = 0.0

            if (present(num_stocks)) num_stocks = ns_tot

        end subroutine SIS_call_tracer_stocks
  
        subroutine SIS_store_stocks(pkg_name, ns, names, units, values, index, stock_values, &
            set_pkg_name, max_ns, ns_tot, stock_names, stock_units)
            character(len=*),                         intent(in)    :: pkg_name
            integer,                                  intent(in)    :: ns
            character(len=*), dimension(:),           intent(in)    :: names, units
            real, dimension(:),                       intent(in)    :: values
            integer,                                  intent(in)    :: index
            real, dimension(:),                       intent(inout) :: stock_values
            character(len=*),                         intent(inout) :: set_pkg_name
            integer,                                  intent(in)    :: max_ns
            integer,                                  intent(inout) :: ns_tot
            character(len=*), dimension(:), optional, intent(inout) :: stock_names, stock_units

            ! This routine stores the stocks and does error handling for call_tracer_stocks.
            character(len=16) :: ind_text, ns_text, max_text
            integer :: n

            if ((index > 0) .and. (ns > 0)) then
                write(ind_text,'(i8)') index
                if (ns > 1) then
                    call SIS_error(FATAL,"Tracer package "//trim(pkg_name)//&
                        " is not permitted to return more than one value when queried"//&
                        " for specific stock index "//trim(adjustl(ind_text))//".")
                elseif (ns+ns_tot > 1) then
                    call SIS_error(FATAL,"Tracer packages "//trim(pkg_name)//" and "//&
                        trim(set_pkg_name)//" both attempted to set values for"//&
                        " specific stock index "//trim(adjustl(ind_text))//".")
                else
                    set_pkg_name = pkg_name
                endif
            endif

            if (ns_tot+ns > max_ns) then
                write(ns_text,'(i8)') ns_tot+ns ; write(max_text,'(i8)') max_ns
                call SIS_error(FATAL,"Attempted to return more tracer stock values (at least "//&
                    trim(adjustl(ns_text))//") than the size "//trim(adjustl(max_text))//&
                    "of the smallest value, name, or units array.")
            endif

            do n=1,ns
                stock_values(ns_tot+n) = values(n)
                if (present(stock_names)) stock_names(ns_tot+n) = names(n)
                if (present(stock_units)) stock_units(ns_tot+n) = units(n)
            enddo
            ns_tot = ns_tot + ns

        end subroutine SIS_store_stocks

        subroutine SIS_tracer_flow_control_end(CS)
            type(SIS_tracer_flow_control_CS), pointer :: CS

!            if (CS%use_ice_age) call ice_age_example_end(CS%ice_age_tracer_CSp)

            if (associated(CS)) deallocate(CS)
        end subroutine SIS_tracer_flow_control_end

    end module SIS_tracer_flow_control
