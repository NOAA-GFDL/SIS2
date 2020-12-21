!> Provides wrapped interfaces for infrastructure calls used by SIS2 that can not be
!! found in the MOM6 framework directory.
module SIS_framework

! This file is part of SIS2. See LICENSE.md for the license.

use fms_io_mod,        only : set_domain, nullify_domain
use fms_io_mod,        only : restart_file_type, FMS1_register_restart=>register_restart_field
use fms_io_mod,        only : save_restart_FMS1=>save_restart, FMS1_restore_state=>restore_state
use fms_io_mod,        only : FMS1_query_initialized=>query_initialized
! use fms2_io_mod,       only : query_initialized=>is_registered_to_restart
use mpp_mod,           only : SIS_chksum=>mpp_chksum
use mpp_domains_mod,   only : domain2D, CENTER, CORNER, EAST, NORTH
use mpp_domains_mod,   only : get_layout=>mpp_get_layout, get_compute_domain=>mpp_get_compute_domain
use mpp_domains_mod,   only : redistribute_data=>mpp_redistribute
use mpp_domains_mod,   only : broadcast_domain=>mpp_broadcast_domain

use MOM_domains,       only : MOM_domain_type
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, NOTE
use MOM_safe_alloc,    only : safe_alloc=>safe_alloc_alloc, safe_alloc_ptr
use MOM_file_parser,   only : get_param, read_param, log_param, log_version, param_file_type

implicit none ; private

public :: SIS_chksum, redistribute_data, domain2D, CENTER, CORNER, EAST, NORTH, axis_names_from_pos
public :: set_domain, nullify_domain, get_layout, get_compute_domain, broadcast_domain
public :: restart_file_type, restore_SIS_state, register_restart_field, save_restart, query_inited
public :: query_initialized, SIS_restart_init, SIS_restart_end, only_read_from_restarts
public :: SIS_initialize_framework, safe_alloc, safe_alloc_ptr

!> A restart registry and the control structure for restarts
type, public :: SIS_restart_CS ! ; private
  type(restart_file_type), pointer :: fms_restart => NULL() !< The FMS restart file type to use
  type(domain2d), pointer :: mpp_domain => NULL() !< The mpp domain to use for read of decomposed fields.
  character(len=240) :: restart_file !< The name or name root for MOM restart files.
  logical :: use_FMS2 = .false.      !< If true use the FMS2 interfaces for restarts.
end type SIS_restart_CS

!> Register fields for restarts
interface register_restart_field
  module procedure register_restart_field_4d
  module procedure register_restart_field_3d
  module procedure register_restart_field_2d
  module procedure register_restart_field_1d
  module procedure register_restart_field_0d
end interface

!> Register fields for restarts
interface only_read_from_restarts
  module procedure only_read_restart_field_4d
  module procedure only_read_restart_field_3d
  module procedure only_read_restart_field_2d
!  module procedure only_read_restart_field_1d
!  module procedure only_read_restart_field_0d
end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_initialize_framework is a template that might be used later to initialize structrures
!! and types that are used for the SIS2 infrastructure.
subroutine SIS_initialize_framework(PF)
  type(param_file_type),   intent(in)    :: PF  !< A structure indicating the open file
                                                !! to parse for model parameter values.

  character(len=200) :: inputdir   ! The directory where NetCDF input files are.
  character(len=200) :: config
  character(len=40)  :: mdl = "SIS_framework" ! This module's name.
  logical :: debug
  ! This include declares and sets the variable "version".
# include "version_variable.h"

  call callTree_enter("SIS_initialize_framework(), SIS_framework.F90")
  call log_version(PF, mdl, version, "")
  call get_param(PF, mdl, "DEBUG", debug, default=.false.)

  call callTree_leave('SIS_initialize_framework()')

end subroutine SIS_initialize_framework


!> Initialize and set up a restart control structure.
subroutine SIS_restart_init(CS, filename, domain, use_FMS2)
  type(SIS_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object that is allocated here
  character(len=*),      intent(in) :: filename !< The path to the restart file.
  type(MOM_domain_type), intent(in) :: domain   !< The MOM domain descriptor being used
  logical,     optional, intent(in) :: use_FMS2 !< If true use the FMS2 variant of calls.

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_restart_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  allocate(CS%fms_restart)
  CS%restart_file = trim(filename)
  CS%mpp_domain => domain%mpp_domain
  CS%use_FMS2 = .false. ; if (present(use_FMS2)) CS%use_FMS2 = use_FMS2

end subroutine SIS_restart_init

!======================= register_restart_field variants =======================

!> Register a 4-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_4d(CS, name, f_ptr, longname, units, position, dim_3, dim_4, &
                                     mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  integer,          optional, intent(in) :: position  !< A coded integer indicating the horizontal
                                                      !! position of this variable
  character(len=*), optional, intent(in) :: dim_3     !< The name of the 3rd dimension axis
  character(len=*), optional, intent(in) :: dim_4     !< The name of the 4th dimension axis
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=24), dimension(5) :: dim_names ! Dimension names to use with FMS2.
  logical :: is_optional
  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_framework: register_restart_field_4d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  if (.not.CS%use_FMS2) then
    ! This is the FMS1 variant of this call.
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, longname=longname, &
                                units=units, mandatory=mandatory, &
                                domain=CS%mpp_domain, position=position)
  else
    ! This is the FMS2 variant of this call.
    dim_names(1:5) = (/ "ih  ", "jh  ", "cat ", "zl  ", "Time" /)
    if (present(position)) call axis_names_from_pos(dim_names(1:2), position, varname=name)
    if (present(dim_3)) dim_names(3) = trim(dim_3)
    if (present(dim_4)) dim_names(4) = trim(dim_4)
    is_optional = .false. ; if (present(mandatory)) is_optional = .not.mandatory
    ! call FMS2_register_restart(CS%FMS2_restfile, name, f_ptr, dim_names, is_optional)
    ! if (present(units)) call FMS2_set_attribute(CS%FMS2_restfile, name, "units", units)
    ! if (present(longname)) call FMS2_set_attribute(CS%FMS2_restfile, name, "longname", longname)
    call SIS_error(FATAL, "register_restart_field_4d: The SIS_framework code does not work with FMS2 yet.")
  endif
end subroutine register_restart_field_4d

!> Register a 3-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_3d(CS, name, f_ptr, longname, units, position, dim_3, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  integer,          optional, intent(in) :: position  !< A coded integer indicating the horizontal
                                                      !! position of this variable
  character(len=*), optional, intent(in) :: dim_3     !< The name of the 3rd dimension axis
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=24), dimension(4) :: dim_names ! Dimension names to use with FMS2.
  logical :: is_optional
  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_framework: register_restart_field_3d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  if (.not.CS%use_FMS2) then
    ! This is the FMS1 variant of this call.
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, longname=longname, &
                                units=units, mandatory=mandatory, &
                                domain=CS%mpp_domain, position=position)
  else
    ! This is the FMS2 variant of this call.
    dim_names(1:4) = (/ "ih  ", "jh  ", "cat ", "Time" /)
    if (present(position)) call axis_names_from_pos(dim_names(1:2), position, varname=name)
    if (present(dim_3)) dim_names(3) = trim(dim_3)
    is_optional = .false. ; if (present(mandatory)) is_optional = .not.mandatory
    ! call FMS2_register_restart(CS%FMS2_restfile, name, f_ptr, dim_names, is_optional)
    ! if (present(units)) call FMS2_set_attribute(CS%FMS2_restfile, name, "units", units)
    ! if (present(longname)) call FMS2_set_attribute(CS%FMS2_restfile, name, "longname", longname)
    call SIS_error(FATAL, "register_restart_field_3d: The SIS_framework code does not work with FMS2 yet.")
  endif

end subroutine register_restart_field_3d

!> Register a 2-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_2d(CS, name, f_ptr, units, longname, position, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  integer,          optional, intent(in) :: position  !< A coded integer indicating the horizontal
                                                      !! position of this variable
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=24), dimension(3) :: dim_names ! Dimension names to use with FMS2.
  logical :: is_optional
  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_framework: register_restart_field_2d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  if (.not.CS%use_FMS2) then
    ! This is the FMS1 variant of this call.
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, longname=longname, &
                                units=units, mandatory=mandatory, &
                                domain=CS%mpp_domain, position=position)

  else
    ! This is the FMS2 variant of this call.
    dim_names(1:3) = (/ "ih  ", "jh  ", "Time" /)
    if (present(position)) call axis_names_from_pos(dim_names(1:2), position, varname=name)
    is_optional = .false. ; if (present(mandatory)) is_optional = .not.mandatory
    ! call FMS2_register_restart(CS%FMS2_restfile, name, f_ptr, dim_names, is_optional)
    ! if (present(units)) call FMS2_set_attribute(CS%FMS2_restfile, name, "units", units)
    ! if (present(longname)) call FMS2_set_attribute(CS%FMS2_restfile, name, "longname", longname)
    call SIS_error(FATAL, "register_restart_field_2d: The SIS_framework code does not work with FMS2 yet.")
  endif

end subroutine register_restart_field_2d

!> Register a 1-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_1d(CS, name, f_ptr, longname, units, dim_name, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real, dimension(:), target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  character(len=*), optional, intent(in) :: dim_name  !< The name of the axis
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=24), dimension(1) :: dim_names ! Dimension names to use with FMS2.
  logical :: is_optional
  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_framework: register_restart_field_1d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  if (.not.CS%use_FMS2) then
    ! This is the FMS1 variant of this call.
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, longname=longname, &
                                units=units, mandatory=mandatory, no_domain=.true.)
  else
    ! This is the FMS2 variant of this call.
    dim_names(1) = "cat" ; if (present(dim_name)) dim_names(1) = dim_name
    is_optional = .false. ; if (present(mandatory)) is_optional = .not.mandatory
    ! call FMS2_register_restart(CS%FMS2_restfile, name, f_ptr, dim_names, is_optional)
    ! if (present(units)) call FMS2_set_attribute(CS%FMS2_restfile, name, "units", units)
    ! if (present(longname)) call FMS2_set_attribute(CS%FMS2_restfile, name, "longname", longname)
    call SIS_error(FATAL, "register_restart_field_1d: The SIS_framework code does not work with FMS2 yet.")
  endif

end subroutine register_restart_field_1d

!> Register a 0-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_0d(CS, name, f_ptr, longname, units, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  logical :: is_optional
  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_framework: register_restart_field_0d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  if (.not.CS%use_FMS2) then
    ! This is the FMS1 variant of this call.
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, longname=longname, &
                                units=units, mandatory=mandatory, no_domain=.true.)
  else
    ! This is the FMS2 variant of this call.
    is_optional = .false. ; if (present(mandatory)) is_optional = .not.mandatory
    ! call FMS2_register_restart(CS%FMS2_restfile, name, f_ptr, dim_names, is_optional)
    ! if (present(units)) call FMS2_set_attribute(CS%FMS2_restfile, name, "units", units)
    ! if (present(longname)) call FMS2_set_attribute(CS%FMS2_restfile, name, "longname", longname)
    call SIS_error(FATAL, "register_restart_field_0d: The SIS_framework code does not work with FMS2 yet.")
  endif

end subroutine register_restart_field_0d

!====================== only_read_from_restarts variants =======================

!> Try to read a named 4-d field from the restart files
subroutine only_read_restart_field_4d(CS, name, f_ptr, position, directory, domain, success)
  type(SIS_restart_CS),            pointer     :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),                intent(in)  :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:,:),target, intent(in)  :: f_ptr     !< A pointer to the field to be read
  integer,               optional, intent(in)  :: position  !< A coded integer indicating the horizontal
                                                            !! position of this variable
  character(len=*),      optional, intent(in)  :: directory !< The directory in which to seek restart files.
  type(MOM_domain_type), optional, intent(in)  :: domain    !< The MOM domain descriptor being used
  logical,               optional, intent(out) :: success   !< True if the field was read successfully

  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_framework: only_read_restart_field_4d: "//&
      "Restart_CS must be initialized before it is used to read "//trim(name))

  if (present(domain)) then
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, read_only=.true., &
                                domain=domain%mpp_domain, mandatory=.false., position=position)
  else
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, read_only=.true., &
                                domain=CS%mpp_domain, mandatory=.false., position=position)
  endif
  call FMS1_restore_state(CS%fms_restart, idr, directory=directory)

  if (present(success)) success = query_initialized(CS, name)

end subroutine only_read_restart_field_4d

!> Try to read a named 3-d field from the restart files
subroutine only_read_restart_field_3d(CS, name, f_ptr, position, directory, domain, success)
  type(SIS_restart_CS),            pointer     :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),                intent(in)  :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:),  target, intent(in)  :: f_ptr     !< A pointer to the field to be read
  integer,               optional, intent(in)  :: position  !< A coded integer indicating the horizontal
                                                            !! position of this variable
  character(len=*),      optional, intent(in)  :: directory !< The directory in which to seek restart files.
  type(MOM_domain_type), optional, intent(in)  :: domain    !< The MOM domain descriptor being used
  logical,               optional, intent(out) :: success   !< True if the field was read successfully

  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_framework: only_read_restart_field_3d: "//&
      "Restart_CS must be initialized before it is used to read "//trim(name))

  if (present(domain)) then
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, read_only=.true., &
                                domain=domain%mpp_domain, mandatory=.false., position=position)
  else
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, read_only=.true., &
                                domain=CS%mpp_domain, mandatory=.false., position=position)
  endif
  call FMS1_restore_state(CS%fms_restart, idr, directory=directory)

  if (present(success)) success = query_initialized(CS, name)

end subroutine only_read_restart_field_3d

!> Try to read a named 2-d field from the restart files
subroutine only_read_restart_field_2d(CS, name, f_ptr, position, directory, domain, success)
  type(SIS_restart_CS),            pointer     :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),                intent(in)  :: name      !< The variable name to be used in the restart file
  real, dimension(:,:),    target, intent(in)  :: f_ptr     !< A pointer to the field to be read
  integer,               optional, intent(in)  :: position  !< A coded integer indicating the horizontal
                                                            !! position of this variable
  character(len=*),      optional, intent(in)  :: directory !< The directory in which to seek restart files.
  type(MOM_domain_type), optional, intent(in)  :: domain    !< The MOM domain descriptor being used
  logical,               optional, intent(out) :: success   !< True if the field was read successfully

  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_framework: only_read_restart_field_2d: "//&
      "Restart_CS must be initialized before it is used to read "//trim(name))

  if (present(domain)) then
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, read_only=.true., &
                                domain=domain%mpp_domain, mandatory=.false., position=position)
  else
    idr = FMS1_register_restart(CS%fms_restart, CS%restart_file, name, f_ptr, read_only=.true., &
                                domain=CS%mpp_domain, mandatory=.false., position=position)
  endif
  call FMS1_restore_state(CS%fms_restart, idr, directory=directory)

  if (present(success)) success = query_initialized(CS, name)

end subroutine only_read_restart_field_2d

!> query_initialized_name determines whether a named field has been successfully
!! read from a restart file yet.
function query_initialized(CS, name) result(query_init)
  type(SIS_restart_CS), pointer    :: CS !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  logical :: query_init        !< The returned value, set to true if the named field has been
                               !! read from a restart file
  !   This subroutine returns .true. if the field referred to by name has
  ! initialized from a restart file, and .false. otherwise.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  query_init = FMS1_query_initialized(CS%fms_restart, name)

end function query_initialized

!> query_inited determines whether a named field has been successfully read from a restart file yet.
!! It is identical to query_initialized, but has a separate name to deal with an unexplained
!! problem that the pgi compiler has with reused function names between modules.
function query_inited(CS, name) result(query_initialized)
  type(SIS_restart_CS), pointer    :: CS !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  logical :: query_initialized !< The returned value, set to true if the named field has been
                               !! read from a restart file
  !   This subroutine returns .true. if the field referred to by name has
  ! initialized from a restart file, and .false. otherwise.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  query_initialized = FMS1_query_initialized(CS%fms_restart, name)

end function query_inited

!> This subroutine sets standard horizontal axis names from a coded position integer
subroutine axis_names_from_pos(dim_names, position, varname)
  character(len=*), dimension(2), intent(out) :: dim_names !< The names of the i- and j-dimensions
  integer,                        intent(in)  :: position  !< A coded integer indicating the horizontal
                                                           !! position of this variable
  character(len=*), optional,     intent(in)  :: varname   !< The variable name to be used in the restart file

  if (position == CENTER) then
    dim_names(1:2) = (/ "ih", "jh" /)
  elseif (position == CORNER) then
    dim_names(1:2) = (/ "iq", "jq" /)
  elseif (position == NORTH) then
    dim_names(1:2) = (/ "iq", "jq" /)
  elseif (position == EAST) then
    dim_names(1:2) = (/ "iq", "jq" /)
  elseif (present(varname)) then
    call SIS_error(FATAL, "set_axis_names_from_pos: Unrecognized position setting for "//trim(varname))
  else
    call SIS_error(FATAL, "set_axis_names_from_pos: Unrecognized position setting")
  endif
end subroutine axis_names_from_pos


!> save_restart saves all registered variables to restart files.
subroutine save_restart(CS, time_stamp)
  type(SIS_restart_CS),        pointer    :: CS !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*), optional , intent(in) :: time_stamp !< A date stamp to use in the restart file name

  call save_restart_FMS1(CS%fms_restart, time_stamp)

end subroutine save_restart

!> Restore the entire state of the sea ice or a single varable using a SIS_restart control structure.
subroutine restore_SIS_state(CS, directory)
  type(SIS_restart_CS),       pointer    :: CS !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*), optional, intent(in) :: directory !< The directory in which to seek restart files.

  call FMS1_restore_state(CS%fms_restart, directory=directory)

end subroutine restore_SIS_state

!> Deallocate memory associated with a MOM_restart_CS variable.
subroutine SIS_restart_end(CS)
  type(SIS_restart_CS),  pointer    :: CS !< A pointer to a MOM_restart_CS object

  if (associated(CS%fms_restart)) deallocate(CS%fms_restart)
  deallocate(CS)

end subroutine SIS_restart_end


end module SIS_framework
