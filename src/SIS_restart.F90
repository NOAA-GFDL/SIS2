!> The SIS2 facility for reading and writing restart files, and querying what has been read.
module SIS_restart

! This file is part of SIS2. See LICENSE.md for the license.

use MOM_checksums,     only : chksum => rotated_field_chksum
use MOM_coupler_types, only : coupler_2d_bc_type, coupler_3d_bc_type
use MOM_domains,       only : PE_here, num_PEs, MOM_domain_type, clone_MOM_domain
use MOM_dyn_horgrid,   only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, is_root_pe, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : MOM_infra_file, MOM_field, create_MOM_file
use MOM_io,            only : MOM_read_data, MOM_write_field
use MOM_io,            only : file_exists, field_exists
use MOM_io,            only : axis_info, set_axis_info, delete_axis_info
use MOM_io,            only : vardesc, var_desc, query_vardesc, modify_vardesc, get_filename_appendix
use MOM_io,            only : MULTIPLE, READONLY_FILE, SINGLE_FILE
use MOM_io,            only : CENTER, CORNER, NORTH_FACE, EAST_FACE
use MOM_string_functions, only : lowercase
use MOM_time_manager,  only : time_type, time_type_to_real, real_to_time
use MOM_time_manager,  only : get_date, set_date
use ice_grid,          only : ice_grid_type
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_transcribe_grid, only : copy_SIS_horgrid_to_dyngrid

implicit none ; private

public :: SIS_restart_init, SIS_restart_end, restore_SIS_state, register_restart_field, save_restart
public :: query_initialized, query_inited, only_read_from_restarts, determine_is_new_run

! These are not used?
public :: restart_files_exist, is_new_run
public :: register_restart_field_as_obsolete, SIS_restart_init_end

!> A type for making arrays of pointers to 4-d arrays
type p4d
  real, dimension(:,:,:,:), pointer :: p => NULL() !< A pointer to a 4d array
end type p4d

!> A type for making arrays of pointers to 3-d arrays
type p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3d array
end type p3d

!> A type for making arrays of pointers to 2-d arrays
type p2d
  real, dimension(:,:), pointer :: p => NULL() !< A pointer to a 2d array
end type p2d

!> A type for making arrays of pointers to 1-d arrays
type p1d
  real, dimension(:), pointer :: p => NULL() !< A pointer to a 1d array
end type p1d

!> A type for making arrays of pointers to scalars
type p0d
  real, pointer :: p => NULL() !< A pointer to a scalar
end type p0d

!> A structure with information about a single restart field
type field_restart
  type(vardesc) :: vars         !< Description of a field that is to be read from or written
                                !! to the restart file.
  logical :: mand_var           !< If .true. the run will abort if this field is not successfully
                                !! read from the restart file.
  logical :: initialized        !< .true. if this field has been read from the restart file.
  character(len=32) :: var_name !< A name by which a variable may be queried.
  real    :: conv = 1.0         !< A factor by which a restart field should be multiplied before it
                                !! is written to a restart file, usually to convert it to MKS or
                                !! other standard units.  When read, the restart field is multiplied
                                !! by the Adcroft reciprocal of this factor.
end type field_restart

!> A structure to store information about restart fields that are no longer used
type obsolete_restart
  character(len=32) :: field_name       !< Name of restart field that is no longer in use
  character(len=32) :: replacement_name !< Name of replacement restart field, if applicable
end type obsolete_restart

!> A restart registry and the control structure for restarts
type, public :: SIS_restart_CS ; private
  logical :: restart                !< restart is set to .true. if the run has been started from a full restart
                                    !! file.  Otherwise some fields must be initialized approximately.
  integer :: novars = 0             !< The number of restart fields that have been registered.
  integer :: num_obsolete_vars = 0  !< The number of obsolete restart fields that have been registered.
  logical :: parallel_restartfiles  !< If true, each PE writes its own restart file,
                                    !! otherwise they are combined internally.
  logical :: new_run                !< If true, the input filenames and restart file existence will
                                    !! result in a new run that is not initialized from restart files.
  logical :: new_run_set = .false.  !< If true, new_run has been determined for this restart_CS.
  logical :: checksum_required      !< If true, require the restart checksums to match and error out otherwise.
                                    !! Users may want to avoid this comparison if for example the restarts are
                                    !! made from a run with a different mask_table than the current run,
                                    !! in which case the checksums will not match and cause crash.
  character(len=240) :: restartfile !< The name or name root for SIS restart files.

  !> An array of descriptions of the registered fields
  type(field_restart), pointer :: restart_field(:) => NULL()

  !> An array of obsolete restart fields
  type(obsolete_restart), pointer :: restart_obsolete(:) => NULL()

  !>@{ Pointers to the fields that have been registered for restarts
  type(p0d), pointer :: var_ptr0d(:) => NULL()
  type(p1d), pointer :: var_ptr1d(:) => NULL()
  type(p2d), pointer :: var_ptr2d(:) => NULL()
  type(p3d), pointer :: var_ptr3d(:) => NULL()
  type(p4d), pointer :: var_ptr4d(:) => NULL()
  !>@}
  integer :: max_fields !< The maximum number of restart fields
end type SIS_restart_CS

!> Register fields for restarts
interface register_restart_field
  module procedure register_restart_field_4d
  module procedure register_restart_field_3d
  module procedure register_restart_field_2d
  module procedure register_restart_field_1d
  module procedure register_restart_field_0d
  module procedure register_restart_coupler_type_3d
  module procedure register_restart_coupler_type_2d
end interface

!> Read optional variables from restart files.
interface only_read_from_restarts
  module procedure only_read_restart_field_4d
  module procedure only_read_restart_field_3d
  module procedure only_read_restart_field_2d
!  module procedure only_read_restart_field_1d
!  module procedure only_read_restart_field_0d
end interface

contains

!> Register a restart field as obsolete
subroutine register_restart_field_as_obsolete(field_name, replacement_name, CS)
  character(*), intent(in) :: field_name       !< Name of restart field that is no longer in use
  character(*), intent(in) :: replacement_name !< Name of replacement restart field, if applicable
  type(SIS_restart_CS), pointer :: CS          !< A pointer to a SIS_restart_CS object (intent in/out)

  CS%num_obsolete_vars = CS%num_obsolete_vars+1
  CS%restart_obsolete(CS%num_obsolete_vars)%field_name = field_name
  CS%restart_obsolete(CS%num_obsolete_vars)%replacement_name = replacement_name
end subroutine register_restart_field_as_obsolete

!> restart_files_exist determines whether any restart files exist.
logical function restart_files_exist(filename, directory, domain, CS)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(MOM_domain_type), intent(in)  :: domain    !< The MOM domain descriptor being used
  type(SIS_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to SIS_restart_init
  integer :: num_files

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "restart_files_exist: Module must be initialized before it is used.")

  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
    num_files = get_num_restart_files('r', directory, domain, CS)
  else
    num_files = get_num_restart_files(filename, directory, domain, CS)
  endif
  restart_files_exist = (num_files > 0)

end function restart_files_exist

!> determine_is_new_run determines from the value of filename and the existence of
!! automatically named restart files in directory whether this would be a new,
!! and as a side effect stores this information in CS.
function determine_is_new_run(filename, directory, G, CS) result(is_new_run)
  character(len=*),        intent(in)  :: filename  !< The list of restart file names or a single
                                                    !! character 'r' to read automatically named files
  character(len=*),        intent(in)  :: directory !< The directory in which to find restart files
  type(SIS_hor_grid_type), intent(in)  :: G         !< The horizontal grid type
  type(SIS_restart_CS),    pointer     :: CS        !< The control structure returned by a previous
                                                    !! call to SIS_restart_init
  logical :: is_new_run                             !< The function result, which indicates whether
                                                    !! this is a new run, based on the value of
                                                    !! filename and whether restart files exist

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "determine_is_new_run: Module must be initialized before it is used.")
  if (LEN_TRIM(filename) > 1) then
    CS%new_run = .false.
  elseif (LEN_TRIM(filename) == 0) then
    CS%new_run = .true.
  elseif (filename(1:1) == 'n') then
    CS%new_run = .true.
  elseif (filename(1:1) == 'F') then
    CS%new_run = (get_num_restart_files('r', directory, G%domain, CS) == 0)
  else
    CS%new_run = .false.
  endif

  CS%new_run_set = .true.
  is_new_run = CS%new_run
end function determine_is_new_run

!> is_new_run returns whether this is going to be a new run based on the
!! information stored in CS by a previous call to determine_is_new_run.
function is_new_run(CS)
  type(SIS_restart_CS),  pointer :: CS !< The control structure returned by a previous
                                       !! call to SIS_restart_init
  logical :: is_new_run                !< The function result, which had been stored in CS during
                                       !! a previous call to determine_is_new_run

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "is_new_run: Module must be initialized before it is used.")
  if (.not.CS%new_run_set) call SIS_error(FATAL, "SIS_restart " // &
      "determine_is_new_run must be called for a restart file before is_new_run.")

  is_new_run = CS%new_run
end function is_new_run

!> open_restart_units determines the number of existing restart files and optionally opens
!! them and returns unit ids and paths to the files, while the result is the number of files
!! that have been opened.
function open_restart_units(filename, directory, domain, CS, IO_handles, file_paths, &
                            global_files) result(num_files)
  character(len=*),      intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(MOM_domain_type), intent(in)  :: domain    !< The MOM domain descriptor being used
  type(SIS_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to SIS_restart_init
  type(MOM_infra_file), dimension(:), &
               optional, intent(out) :: IO_handles !< The I/O handles of all opened files
  character(len=*), dimension(:), &
               optional, intent(out) :: file_paths !< The full paths to the restart files
  logical, dimension(:), &
               optional, intent(out) :: global_files !< True if a file is global

  integer :: num_files  !< The number of files (both automatically named restart files
                        !! and others explicitly in filename) that have been opened.

  ! Local variables
  character(len=256) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=256) :: fname     ! The name of the current file.
  character(len=8)   :: suffix    ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
  integer :: num_restart     ! The number of restart files that have already
                             ! been opened using their numbered suffix.
  integer :: start_char      ! The location of the starting character in the
                             ! current file name.
  integer :: nf              ! The number of files that have been found so far
  integer :: m, length
  logical :: still_looking   ! If true, the code is still looking for automatically named files
  logical :: fexists         ! True if a file has been found
  character(len=32) :: filename_appendix = '' ! Filename appendix for ensemble runs
  character(len=80) :: restartname
  character(len=32) :: count_msg

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "open_restart_units: Module must be initialized before it is used.")

  ! Get NetCDF ids for all of the restart files.
  num_restart = 0 ; nf = 0 ; start_char = 1
  do while (start_char <= len_trim(filename) )
    do m=start_char,len_trim(filename)
      if (filename(m:m) == ' ') exit
    enddo
    fname = filename(start_char:m-1)
    start_char = m
    do while (start_char <= len_trim(filename))
      if (filename(start_char:start_char) == ' ') then
        start_char = start_char + 1
      else
        exit
      endif
    enddo

    if (((fname(1:1)=='r') .or. (fname(1:1)=='F')) .and. ( len_trim(fname) == 1)) then
      still_looking = (num_restart <= 0) ! Avoid going through the file list twice.
      do while (still_looking)
        restartname = trim(CS%restartfile)

        ! Determine if there is a filename_appendix (used for ensemble runs).
        call get_filename_appendix(filename_appendix)
        if (len_trim(filename_appendix) > 0) then
          length = len_trim(restartname)
          if (restartname(length-2:length) == ".nc") then
            restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//".nc"
          else
            restartname = restartname(1:length)  //'.'//trim(filename_appendix)
          endif
        endif
        filepath = trim(directory) // trim(restartname)

        if (num_restart < 10) then
          write(suffix,'("_",I1)') num_restart
        else
          write(suffix,'("_",I2)') num_restart
        endif
        length = len_trim(filepath)
        if (length < 3) then
          if (num_restart > 0) filepath = trim(filepath) // suffix
          filepath = trim(filepath)//".nc"
        elseif (filepath(length-2:length) == ".nc") then
          if (num_restart > 0) filepath = filepath(1:length-3)//trim(suffix)//".nc"
        else
          if (num_restart > 0) filepath = trim(filepath) // suffix
          filepath = trim(filepath)//".nc"
        endif

        num_restart = num_restart + 1
        ! Look for a global netCDF file.
        inquire(file=filepath, exist=fexists)
        if (fexists) then
          nf = nf + 1
          if (present(IO_handles)) &
            call IO_handles(nf)%open(trim(filepath), READONLY_FILE, &
                MOM_domain=domain, threading=MULTIPLE, fileset=SINGLE_FILE)
          if (present(global_files)) global_files(nf) = .true.
          if (present(file_paths)) file_paths(nf) = filepath
        elseif (CS%parallel_restartfiles) then
          ! Look for decomposed files using the I/O Layout.
          fexists = file_exists(filepath, Domain)
          if (fexists) then
            nf = nf + 1
            if (present(IO_handles)) &
              call IO_handles(nf)%open(trim(filepath), READONLY_FILE, &
                  MOM_domain=domain, threading=MULTIPLE, fileset=MULTIPLE)
            if (present(global_files)) global_files(nf) = .false.
            if (present(file_paths)) file_paths(nf) = filepath
          endif
        endif

        if (fexists) then
          if (is_root_pe() .and. (present(IO_handles))) &
            call SIS_mesg("SIS_restart: SIS run restarted using : "//trim(filepath))
        else
          still_looking = .false. ; exit
        endif
      enddo ! while (still_looking) loop
    else
      filepath = trim(directory)//trim(fname)
      inquire(file=filepath, exist=fexists)
      if (.not. fexists) then
        filepath = trim(filepath)//".nc"
        inquire(file=filepath, exist=fexists)
      endif
      if (fexists) then
        nf = nf + 1
        if (present(IO_handles)) &
          call IO_handles(nf)%open(trim(filepath), READONLY_FILE, &
              threading=MULTIPLE, fileset=SINGLE_FILE)
        if (present(global_files)) global_files(nf) = .true.
        if (present(file_paths)) file_paths(nf) = filepath
        if (is_root_pe() .and. (present(IO_handles))) &
          call SIS_mesg("SIS_restart: SIS run restarted using : "//trim(filepath))
      else
        if (present(IO_handles)) &
          call SIS_error(WARNING,"SIS_restart: Unable to find restart file : "//trim(filepath))
      endif

    endif
  enddo ! while (start_char < len_trim(filename)) loop

  write(count_msg, '(I0)') nf
  call SIS_mesg("SIS2: open_restart_units found "//trim(count_msg)//" files using "//&
                trim(filename)//" in directory "//trim(directory), 9)

  num_files = nf

end function open_restart_units

!> get_num_restart_files returns the number of existing restart files that match the provided
!! directory structure and other information stored in the control structure and optionally
!! also provides the full paths to these files.
function get_num_restart_files(filenames, directory, domain, CS) result(num_files)
  character(len=*),      intent(in)  :: filenames !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),      intent(in)  :: directory !< The directory in which to find restart files
  type(MOM_domain_type), intent(in)  :: domain    !< The MOM domain descriptor being used
  type(SIS_restart_CS),  pointer     :: CS        !< The control structure returned by a previous
                                                  !! call to SIS_restart_init
  integer :: num_files  !< The function result, the number of files (both automatically named
                        !! restart files and others explicitly in filename) that have been opened

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "get_num_restart_files: Module must be initialized before it is used.")

  ! This call uses open_restart_units without the optional arguments needed to actually
  ! open the files to determine the number of restart files.
  num_files = open_restart_units(filenames, directory, domain, CS)

end function get_num_restart_files


!> Initialize this module and set up a restart control structure.
subroutine SIS_restart_init(CS, filename, domain, param_file)
  type(SIS_restart_CS),  pointer    :: CS !< A pointer to a SIS_restart_CS object that is allocated here
  character(len=*),      intent(in) :: filename !< The path to the restart file.
  type(MOM_domain_type), intent(in) :: domain   !< The MOM domain descriptor being used
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters

!  logical :: rotate_index

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "SIS_restart"   ! This module's name.
  logical :: all_default   ! If true, all parameters are using their default values.

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_restart_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  ! Determine whether all parameters are set to their default values.
  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", CS%parallel_restartfiles, &
                 default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "MAX_FIELDS", CS%max_fields, default=100, do_not_log=.true.)
  call get_param(param_file, mdl, "RESTART_CHECKSUMS_REQUIRED", CS%checksum_required, &
                 default=.true., do_not_log=.true.)
  all_default = ((.not.CS%parallel_restartfiles) .and. (CS%max_fields == 100) .and. (CS%checksum_required))
  CS%restartfile = trim(filename)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "", all_default=all_default)
  call get_param(param_file, mdl, "PARALLEL_RESTARTFILES", CS%parallel_restartfiles, &
                 "If true, each processor writes its own restart file, "//&
                 "otherwise a single restart file is generated", &
                 default=.false.)
  call get_param(param_file, mdl, "MAX_FIELDS", CS%max_fields, &
                 "The maximum number of restart fields that can be used.", &
                 default=100)
  call get_param(param_file, mdl, "RESTART_CHECKSUMS_REQUIRED", CS%checksum_required, &
                 "If true, require the restart checksums to match and error out otherwise. "//&
                 "Users may want to avoid this comparison if for example the restarts are "//&
                 "made from a run with a different mask_table than the current run, "//&
                 "in which case the checksums will not match.", default=.true.)

  allocate(CS%restart_field(CS%max_fields))
  allocate(CS%restart_obsolete(CS%max_fields))
  allocate(CS%var_ptr0d(CS%max_fields))
  allocate(CS%var_ptr1d(CS%max_fields))
  allocate(CS%var_ptr2d(CS%max_fields))
  allocate(CS%var_ptr3d(CS%max_fields))
  allocate(CS%var_ptr4d(CS%max_fields))

end subroutine SIS_restart_init

!======================= register_restart_field variants =======================

!> Register a 4-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_4d(CS, name, f_ptr, longname, units, conversion, &
                                     position, dim_3, dim_4, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written, 1 by default.
  integer,          optional, intent(in) :: position  !< A coded integer indicating the horizontal
                                                      !! position of this variable
  character(len=*), optional, intent(in) :: dim_3     !< The name of the 3rd dimension axis
  character(len=*), optional, intent(in) :: dim_4     !< The name of the 4th dimension axis
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=32), dimension(4) :: dim_names ! Non-time dimension names to use with this variable
  character(len=:), allocatable :: long_name   ! The long name, or the name if the longname is not provided.
  character(len=:), allocatable :: var_units   ! The variable units, or 'none' if they are not provided.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_4d: " // &
      "Module must be initialized before it is used to register "//trim(name))

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                        ! once the total number of fields is known.

  CS%restart_field(CS%novars)%var_name = trim(name)
  dim_names(:) = "" ; dim_names(3) = "cat" ; dim_names(4) = "z_ice"
  if (present(dim_3)) dim_names(3) = trim(dim_3)
  if (present(dim_4)) dim_names(4) = trim(dim_4)
  long_name = trim(name) ; if (present(longname)) long_name = trim(longname)
  var_units = "none" ; if (present(units)) var_units = trim(units)
  CS%restart_field(CS%novars)%vars = var_desc(name, units=var_units, longname=long_name, &
                                              position=position, dim_names=dim_names)
  CS%restart_field(CS%novars)%mand_var = set_from_optional(.true., mandatory)

  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion

  CS%var_ptr4d(CS%novars)%p => f_ptr
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_4d

!> Register a 3-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_3d(CS, name, f_ptr, longname, units, conversion, position, dim_3, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written, 1 by default.
  integer,          optional, intent(in) :: position  !< A coded integer indicating the horizontal
                                                      !! position of this variable
  character(len=*), optional, intent(in) :: dim_3     !< The name of the 3rd dimension axis
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=32), dimension(3) :: dim_names ! Non-time dimension names to use with this variable
  character(len=:), allocatable :: long_name   ! The long name, or the name if the longname is not provided.
  character(len=:), allocatable :: var_units   ! The variable units, or 'none' if they are not provided.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_3d: " // &
      "Module must be initialized before it is used to register "//trim(name))

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                        ! once the total number of fields is known.

  CS%restart_field(CS%novars)%var_name = trim(name)
  dim_names(:) = "" ; dim_names(3) = "cat"
  if (present(dim_3)) dim_names(3) = trim(dim_3)
  long_name = trim(name) ; if (present(longname)) long_name = trim(longname)
  var_units = "none" ; if (present(units)) var_units = trim(units)
  CS%restart_field(CS%novars)%vars = var_desc(name, units=var_units, longname=long_name, &
                                              position=position, dim_names=dim_names)
  CS%restart_field(CS%novars)%mand_var = set_from_optional(.true., mandatory)

  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion

  CS%var_ptr3d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_3d

!> Register a 2-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_2d(CS, name, f_ptr, units, conversion, longname, position, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real, dimension(:,:), &
                      target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written, 1 by default.
  integer,          optional, intent(in) :: position  !< A coded integer indicating the horizontal
                                                      !! position of this variable
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=:), allocatable :: long_name   ! The long name, or the name if the longname is not provided.
  character(len=:), allocatable :: var_units   ! The variable units, or 'none' if they are not provided.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_2d: " // &
      "Module must be initialized before it is used to register "//trim(name))

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                        ! once the total number of fields is known.

  CS%restart_field(CS%novars)%var_name = trim(name)
  long_name = trim(name) ; if (present(longname)) long_name = trim(longname)
  var_units = "none" ; if (present(units)) var_units = trim(units)
  CS%restart_field(CS%novars)%vars = var_desc(name, units=var_units, longname=long_name, &
                                              z_grid='1', position=position)
  CS%restart_field(CS%novars)%mand_var = set_from_optional(.true., mandatory)
  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion

  CS%var_ptr2d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_2d

!> Register a 1-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_1d(CS, name, f_ptr, longname, units, conversion, dim_name, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real, dimension(:), target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written, 1 by default.
  character(len=*), optional, intent(in) :: dim_name  !< The name of the axis
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=:), allocatable :: long_name   ! The long name, or the name if the longname is not provided.
  character(len=:), allocatable :: var_units   ! The variable units, or 'none' if they are not provided.
  character(len=32), dimension(1) :: dim_names ! Non-time dimension names to use with this variable

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_1d: " // &
      "Module must be initialized before it is used to register "//trim(name))

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                        ! once the total number of fields is known.

  CS%restart_field(CS%novars)%var_name = trim(name)
  dim_names(:) = (/"cat"/)
  if (present(dim_name)) dim_names(1) = trim(dim_name)
  long_name = trim(name) ; if (present(longname)) long_name = trim(longname)
  var_units = "none" ; if (present(units)) var_units = trim(units)
  CS%restart_field(CS%novars)%vars = var_desc(name, units=var_units, longname=long_name, &
                                              position=0, dim_names=dim_names)
  CS%restart_field(CS%novars)%mand_var = set_from_optional(.true., mandatory)
  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion

  CS%var_ptr1d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr0d(CS%novars)%p => NULL()

end subroutine register_restart_field_1d

!> Register a 0-d field for restarts, providing the metadata as individual arguments
subroutine register_restart_field_0d(CS, name, f_ptr, longname, units, conversion, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),           intent(in) :: name      !< The variable name to be used in the restart file
  real,               target, intent(in) :: f_ptr     !< A pointer to the field to be read or written
  character(len=*), optional, intent(in) :: longname  !< variable long name
  character(len=*), optional, intent(in) :: units     !< variable units
  real,             optional, intent(in) :: conversion !< A factor to multiply a restart field by
                                                      !! before it is written, 1 by default.
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is not
                                                      !! successfully read from the restart file.

  ! Local variables
  character(len=:), allocatable :: long_name   ! The long name, or the name if the longname is not provided.
  character(len=:), allocatable :: var_units   ! The variable units, or 'none' if they are not provided.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_0d: " // &
      "Module must be initialized before it is used to register "//trim(name))

  CS%novars = CS%novars+1
  if (CS%novars > CS%max_fields) return ! This is an error that will be reported
                                        ! once the total number of fields is known.

  CS%restart_field(CS%novars)%var_name = trim(name)
  long_name = trim(name) ; if (present(longname)) long_name = trim(longname)
  var_units = "none" ; if (present(units)) var_units = trim(units)
  CS%restart_field(CS%novars)%vars = var_desc(name, units=var_units, longname=long_name, &
                                              position=0, z_grid='1')
  CS%restart_field(CS%novars)%mand_var = set_from_optional(.true., mandatory)
  CS%restart_field(CS%novars)%initialized = .false.
  CS%restart_field(CS%novars)%conv = 1.0
  if (present(conversion)) CS%restart_field(CS%novars)%conv = conversion

  CS%var_ptr0d(CS%novars)%p => f_ptr
  CS%var_ptr4d(CS%novars)%p => NULL()
  CS%var_ptr3d(CS%novars)%p => NULL()
  CS%var_ptr2d(CS%novars)%p => NULL()
  CS%var_ptr1d(CS%novars)%p => NULL()

end subroutine register_restart_field_0d


!> Register a 3-d coupler-type field for restarts.  The coupler type provides the relevant metadata.
subroutine register_restart_coupler_type_3d(CS, bc_ptr, varname_prefix, dim_3, mandatory)
  type(SIS_restart_CS),       pointer       :: CS     !< A pointer to a SIS_restart_CS object (intent in/out)
  type(coupler_3d_bc_type),   intent(inout) :: bc_ptr !< The coupler type field to be read or written
  character(len=*), optional, intent(in)    :: varname_prefix !< A prefix for the variable name in the
                                                      !! restart file, intended to allow multiple
                                                      !! BC_type variables to use the same restart
                                                      !! files.
  character(len=*), optional, intent(in) :: dim_3     !< The name of the 3rd dimension axis
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is
                                                      !! not successfully read from the restart file.

  ! Local variables
  character(len=32), dimension(3) :: dim_names ! Non-time dimension names to use with this variable
  character(len=:), allocatable :: var_name
  integer :: n, m

  if ((.not.associated(CS)) .and. (present(varname_prefix))) then
    call SIS_error(FATAL, "SIS_framework: register_restart_coupler_type_2d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(varname_prefix))
  elseif (.not.associated(CS)) then
    call SIS_error(FATAL, "SIS_framework: register_restart_coupler_type_2d: "//&
      "Restart_CS must be initialized before it is used to register an unspecified bc_type")
  endif

  dim_names(:) = "" ; dim_names(3) = "cat"
  if (present(dim_3)) dim_names(3) = trim(dim_3)

  do n = 1, bc_ptr%num_bcs
    do m = 1, bc_ptr%bc(n)%num_fields
      CS%novars = CS%novars+1
      if (CS%novars > CS%max_fields) cycle  ! Report this error once the total number of fields is known.

      var_name = trim(bc_ptr%bc(n)%field(m)%name)
      if (present(varname_prefix)) var_name = trim(varname_prefix)//trim(var_name)
      CS%restart_field(CS%novars)%var_name = trim(var_name)
      CS%restart_field(CS%novars)%vars = var_desc(var_name, units=bc_ptr%bc(n)%field(m)%units, &
                                                  longname=bc_ptr%bc(n)%field(m)%long_name, &
                                                  position=CENTER, dim_names=dim_names)
      CS%restart_field(CS%novars)%mand_var = set_from_optional(.true., mandatory)
      CS%restart_field(CS%novars)%initialized = .false.
      CS%restart_field(CS%novars)%conv = 1.0

      CS%var_ptr3d(CS%novars)%p => bc_ptr%bc(n)%field(m)%values
      CS%var_ptr4d(CS%novars)%p => NULL() ; CS%var_ptr2d(CS%novars)%p => NULL()
      CS%var_ptr1d(CS%novars)%p => NULL() ; CS%var_ptr0d(CS%novars)%p => NULL()
    enddo
  enddo

end subroutine register_restart_coupler_type_3d


!> Register a 2-d coupler-type field for restarts.  The coupler type provides the relevant metadata.
subroutine register_restart_coupler_type_2d(CS, bc_ptr, varname_prefix, mandatory)
  type(SIS_restart_CS),       pointer    :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  type(coupler_2d_bc_type),   intent(in) :: bc_ptr    !< The coupler type field to be written
  character(len=*), optional, intent(in) :: varname_prefix !< A prefix for the variable name in the
                                                      !! restart file, intended to allow multiple
                                                      !! BC_type variables to use the same restart
                                                      !! files.
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is
                                                      !! not successfully read from the restart file.

  ! Local variables
  character(len=:), allocatable :: var_name
  integer :: n, m

  if ((.not.associated(CS)) .and. (present(varname_prefix))) then
    call SIS_error(FATAL, "SIS_framework: register_restart_coupler_type_2d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(varname_prefix))
  elseif (.not.associated(CS)) then
    call SIS_error(FATAL, "SIS_framework: register_restart_coupler_type_2d: "//&
      "Restart_CS must be initialized before it is used to register an unspecified bc_type")
  endif

  do n = 1, bc_ptr%num_bcs
    do m = 1, bc_ptr%bc(n)%num_fields
      CS%novars = CS%novars+1
      if (CS%novars > CS%max_fields) cycle  ! Report this error once the total number of fields is known.

      var_name = trim(bc_ptr%bc(n)%field(m)%name)
      if (present(varname_prefix)) var_name = trim(varname_prefix)//trim(var_name)
      CS%restart_field(CS%novars)%var_name = trim(var_name)
      CS%restart_field(CS%novars)%vars = var_desc(var_name, units=bc_ptr%bc(n)%field(m)%units, &
                                                  longname=bc_ptr%bc(n)%field(m)%long_name, &
                                                  z_grid='1', position=CENTER)
      CS%restart_field(CS%novars)%mand_var = set_from_optional(.true., mandatory)
      CS%restart_field(CS%novars)%initialized = .false.
      CS%restart_field(CS%novars)%conv = 1.0

      CS%var_ptr2d(CS%novars)%p => bc_ptr%bc(n)%field(m)%values
      CS%var_ptr4d(CS%novars)%p => NULL() ; CS%var_ptr3d(CS%novars)%p => NULL()
      CS%var_ptr1d(CS%novars)%p => NULL() ; CS%var_ptr0d(CS%novars)%p => NULL()
    enddo
  enddo

end subroutine register_restart_coupler_type_2d

!====================== only_read_from_restarts variants =======================

!> Try to read a named 4-d field from the restart files
subroutine only_read_restart_field_4d(CS, name, f_ptr, domain, position, directory, success, scale)
  type(SIS_restart_CS),            pointer       :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),                intent(in)    :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:,:),        intent(inout) :: f_ptr     !< The array for the field to be read
  type(MOM_domain_type),           intent(in)    :: domain    !< The MOM domain descriptor being used
  integer,               optional, intent(in)    :: position  !< A coded integer indicating the horizontal
                                                              !! position of this variable
  character(len=*),      optional, intent(in)    :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out)   :: success   !< True if the field was read successfully
  real,                  optional, intent(in)    :: scale     !< A factor by which the field will be scaled

  ! Local variables
  character(len=240), allocatable, dimension(:) :: file_paths ! The possible file names.
  character(len=240) :: dir ! The directory to read from.
  logical, allocatable, dimension(:) :: global_file  ! True if the file is global
  integer :: n, num_files

  dir = "./INPUT/" ; if (present(directory)) dir = trim(directory)

  if (present(success)) success=.false.

  num_files = get_num_restart_files('r', directory, domain, CS)
  if (num_files == 0) return
  allocate(file_paths(num_files), global_file(num_files))
  num_files = open_restart_units('r', directory, domain, CS, file_paths=file_paths, global_files=global_file)

  do n=1,num_files ; if (field_exists(file_paths(n), name, MOM_Domain=domain)) then
    call MOM_read_data(file_paths(n), name, f_ptr, domain, timelevel=1, position=position, &
                       scale=scale, global_file=global_file(n))
    if (present(success)) success=.true.
    exit
  endif ; enddo

  deallocate(file_paths, global_file)

end subroutine only_read_restart_field_4d

!> Try to read a named 3-d field from the restart files
subroutine only_read_restart_field_3d(CS, name, f_ptr, domain, position, directory, success, scale)
  type(SIS_restart_CS),            pointer       :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),                intent(in)    :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:),          intent(inout) :: f_ptr     !< The array for the field to be read
  type(MOM_domain_type),           intent(in)    :: domain    !< The MOM domain descriptor being used
  integer,               optional, intent(in)    :: position  !< A coded integer indicating the horizontal
                                                              !! position of this variable
  character(len=*),      optional, intent(in)    :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out)   :: success   !< True if the field was read successfully
  real,                  optional, intent(in)    :: scale     !< A factor by which the field will be scaled

  ! Local variables
  character(len=240), allocatable, dimension(:) :: file_paths ! The possible file names
  character(len=240) :: dir ! The directory to read from.
  logical, allocatable, dimension(:) :: global_file  ! True if the file is global
  integer :: n, num_files

  dir = "./INPUT/" ; if (present(directory)) dir = trim(directory)

  if (present(success)) success=.false.

  num_files = get_num_restart_files('r', directory, domain, CS)
  if (num_files == 0) return
  allocate(file_paths(num_files), global_file(num_files))
  num_files = open_restart_units('r', directory, domain, CS, file_paths=file_paths, global_files=global_file)

  do n=1,num_files ; if (field_exists(file_paths(n), name, MOM_Domain=domain)) then
    call MOM_read_data(file_paths(n), name, f_ptr, domain, timelevel=1, position=position, &
                       scale=scale, global_file=global_file(n), file_may_be_4d=.true.)
    if (present(success)) success=.true.
    exit
  endif ; enddo

  deallocate(file_paths, global_file)

end subroutine only_read_restart_field_3d

!> Try to read a named 2-d field from the restart files
subroutine only_read_restart_field_2d(CS, name, f_ptr, domain, position, directory, success, scale)
  type(SIS_restart_CS),            pointer       :: CS        !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*),                intent(in)    :: name      !< The variable name to be used in the restart file
  real, dimension(:,:),            intent(inout) :: f_ptr     !< The array for the field to be read
  type(MOM_domain_type),           intent(in)    :: domain    !< The MOM domain descriptor being used
  integer,               optional, intent(in)    :: position  !< A coded integer indicating the horizontal
                                                              !! position of this variable
  character(len=*),      optional, intent(in)    :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out)   :: success   !< True if the field was read successfully
  real,                  optional, intent(in)    :: scale     !< A factor by which the field will be scaled

  ! Local variables
  character(len=240), allocatable, dimension(:) :: file_paths ! The possible file names.
  character(len=240) :: dir ! The directory to read from.
  logical, allocatable, dimension(:) :: global_file  ! True if the file is global
  integer :: n, num_files

  dir = "./INPUT/" ; if (present(directory)) dir = trim(directory)

  if (present(success)) success=.false.

  num_files = get_num_restart_files('r', directory, domain, CS)
  if (num_files == 0) return
  allocate(file_paths(num_files), global_file(num_files))
  num_files = open_restart_units('r', directory, domain, CS, file_paths=file_paths, global_files=global_file)

  do n=1,num_files ; if (field_exists(file_paths(n), name, MOM_Domain=domain)) then
    call MOM_read_data(file_paths(n), name, f_ptr, domain, timelevel=1, position=position, &
                       scale=scale, global_file=global_file(n), file_may_be_4d=.true.)
    if (present(success)) success=.true.
    exit
  endif ; enddo

  deallocate(file_paths, global_file)

end subroutine only_read_restart_field_2d

!> query_initialized determines whether a named field has been successfully read from a restart file yet.
logical function query_initialized(CS, name, will_initialize)
  type(SIS_restart_CS), pointer    :: CS   !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  logical,    optional, intent(in) :: will_initialize  !< If present and true, report this variable
                                           !! as having been initialized in future calls.
  integer :: m, n

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "query_initialized: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  query_initialized = .false.
  n = CS%novars+1
  do m=1,CS%novars
    if (trim(name) == CS%restart_field(m)%var_name) then
      if (CS%restart_field(m)%initialized) query_initialized = .true.
      n = m ; exit
    endif
  enddo
  ! Assume that you are going to initialize it now, so set flag to initialized if queried again.
  if ((n<=CS%novars) .and. present(will_initialize)) then
    if (will_initialize) CS%restart_field(n)%initialized = .true.
  endif

  if ((n==CS%novars+1) .and. (is_root_pe())) &
    call SIS_mesg("SIS_restart: Unknown restart variable "//name//" queried for initialization.")

  if ((is_root_pe()) .and. query_initialized) &
    call SIS_mesg("SIS_restart: "//name//" initialization confirmed by name.")

end function query_initialized

!> query_inited determines whether a named field has been successfully read from a restart file yet.
!! It is identical to query_initialized, but has a separate name to deal with an unexplained
!! problem that the pgi compiler has with reused function names between modules.
logical function query_inited(CS, name, will_initialize)
  type(SIS_restart_CS), pointer    :: CS   !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  logical,    optional, intent(in) :: will_initialize  !< If present and true, report this variable
                                           !! as having been initialized in future calls.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "query_inited: Module must be initialized before it is used.")

  query_inited = query_initialized(CS, name, will_initialize)

end function query_inited

!> Indicate that all variables have now been registered.
subroutine SIS_restart_init_end(CS)
  type(SIS_restart_CS),  pointer    :: CS !< A pointer to a SIS_restart_CS object

  if (associated(CS)) then
    if (CS%novars == 0) call SIS_restart_end(CS)
  endif

end subroutine SIS_restart_init_end

!> Set a value from an optional argument or a default value.
logical function set_from_optional(default, opt_val)
  logical,           intent(in) :: default !< The default value if opt_val is not present.
  logical, optional, intent(in) :: opt_val !< The optional value to use if present.

  set_from_optional = default
  if (present(opt_val)) set_from_optional = opt_val
end function set_from_optional

!> Issue a fatal error message as necessary.
subroutine restart_error(CS)
  type(SIS_restart_CS),  pointer    :: CS !< A pointer to a SIS_restart_CS object

  character(len=16)  :: num  ! String for error messages

  if (CS%novars > CS%max_fields) then
    write(num,'(I0)') CS%novars
    call SIS_error(FATAL,"SIS_restart: Too many fields registered for " // &
           "restart.  Set MAX_FIELDS to be at least " // &
           trim(adjustl(num)) // " in the MOM input file.")
  else
    call SIS_error(FATAL,"SIS_restart: Unspecified fatal error.")
  endif
end subroutine restart_error

!> Return bounds for computing checksums to store in restart files, or to read from them.
!! Note that some sea-ice restart variables are allocated only on the computational grid and
!! use the coupler's indexing convention, not that of the sea-ice, so offsets are needed.
!! The upper bounds are aligned to deal with the haloless symmetric memory case.
subroutine get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL, lbnd, ubnd)
  type(SIS_hor_grid_type), intent(in)  :: G   !< The horizontal grid type
  integer,                 intent(in)  :: pos !< An integer indicating staggering of variable
  integer,                 intent(out) :: isL !< i-start for checksum
  integer,                 intent(out) :: ieL !< i-end for checksum
  integer,                 intent(out) :: jsL !< j-start for checksum
  integer,                 intent(out) :: jeL !< j-end for checksum
  integer, dimension(:),   intent(in)  :: lbnd !< The lower index bounds of the array
  integer, dimension(:),   intent(in)  :: ubnd !< The upper index bounds of the array


  if ((size(lbnd) < 2) .or. (size(ubnd) < 2)) call SIS_error(FATAL, &
    "get_checksum_loop_ranges called with lbnd or ubnd arrays that are too small.")

  ! Regular non-symmetric compute domain
  isL = G%isc-G%isd+1
  ieL = G%iec-G%isd+1
  jsL = G%jsc-G%jsd+1
  jeL = G%jec-G%jsd+1

  ! Expand range east or south for symmetric arrays
  if (G%symmetric) then
    if ((pos == EAST_FACE) .or. (pos == CORNER)) then ! For u-, q-points only
      if (G%idg_offset == 0) isL = isL - 1 ! include western edge in checksums only for western PEs
    endif
    if ((pos == NORTH_FACE) .or. (pos == CORNER)) then ! For v-, q-points only
      if (G%jdg_offset == 0) jsL = jsL - 1 ! include western edge in checksums only for southern PEs
    endif
  endif

  ! Adjust the loop ranges for arrays that are allocated only on the computational domain and may
  ! use different indexing conventions than those in G.
  if ((ubnd(1)-lbnd(1) == ieL-isL) .or. (ubnd(1)-lbnd(1) == 1+ieL-isL)) then
    isL = ubnd(1) + isL - ieL
    ieL = ubnd(1)
  endif
  if ((ubnd(2)-lbnd(2) == jeL-jsL) .or. (ubnd(2)-lbnd(2) == 1+jeL-jsL)) then
    jsL = ubnd(2) + jsL - jeL
    jeL = ubnd(2)
  endif

end subroutine get_checksum_loop_ranges

!> save_restart saves all registered variables to restart files.
subroutine save_restart(directory, time, G, CS, IG, time_stamp)
  character(len=*),           intent(in)    :: directory !< The directory where the restart files
                                                     !! are to be written
  type(time_type),            intent(in)    :: time  !< The current model time
  type(SIS_hor_grid_type),    intent(inout) :: G     !< The horizontal grid type
  type(SIS_restart_CS),       pointer       :: CS    !< The control structure returned by a previous
                                                     !! call to SIS_restart_init
  type(ice_grid_type), &
                    optional, intent(in)    :: IG    !< The sea-ice grid type
  character(len=*), optional, intent(in)    :: time_stamp !< A date stamp to include in the restart file name

  ! Local variables
  type(vardesc) :: vars(CS%max_fields)  ! Descriptions of the fields that
                                        ! are to be read from the restart file.
  type(MOM_field) :: fields(CS%max_fields) ! Opaque types containing metadata describing
                                        ! each variable that will be written.
  type(axis_info)    :: extra_axes(5)   ! Descriptors for extra axes that might be used
  type(dyn_horgrid_type), pointer :: dG => NULL() ! Common horizontal grid type between SIS2 and MOM6

  character(len=512) :: restartpath     ! The restart file path (dir/file).
  character(len=256) :: restartname     ! The restart file name (no dir).
  character(len=8)   :: suffix          ! A suffix (like _2) that is appended
                                        ! to the name of files after the first.
  integer(kind=8) :: var_sz, size_in_file ! The size in bytes of each variable
                                        ! and the variables already in a file.
  integer(kind=8) :: max_file_size = 4294967292_8 ! The maximum size in bytes
                                        ! for any one file.
  integer :: start_var, next_var        ! The starting variables of the
                                        ! current and next files.
  type(MOM_infra_file) :: IO_handle     ! The I/O handle of the open fileset
  integer :: file_thread                ! A flag indicating whether to use parallel restart files
  integer :: m, n, nz
  integer :: num_files                  ! The number of restart files that will be used.
  character(len=8) :: hor_grid, z_grid, t_grid ! Variable grid info.
  character(len=64) :: var_name         ! A variable's name.
  real :: conv                          ! Shorthand for the conversion factor
  real :: restart_time
  character(len=32), dimension(5) :: dim_names ! Non-time dimension names to use with this variable
  character(len=32) :: filename_appendix = '' ! Appendix to filename for ensemble runs
  integer :: length                     ! The length of a text string.
  integer(kind=8) :: check_val(CS%max_fields,1)
  integer :: isL, ieL, jsL, jeL, pos

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "save_restart: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  num_files = 0
  next_var = 0
  nz = 1 ; if (present(IG)) nz = IG%NkIce

  restart_time = time_type_to_real(time) / 86400.0

  restartname = trim(CS%restartfile)
  if (PRESENT(time_stamp)) restartname = trim(CS%restartfile)//trim(time_stamp)

  if (present(IG)) then
    call set_axis_info(extra_axes(1), "cat", longname="Ice thickness categories", ax_size=IG%CatIce)
    call set_axis_info(extra_axes(2), "cat0", longname="Ice thickness categories with open water", ax_size=IG%CatIce+1)
    call set_axis_info(extra_axes(3), "z_ice", longname="Ice vertical layers", &
                       cartesian='Z', sense=-1, ax_size=IG%NkIce)
    call set_axis_info(extra_axes(4), "z_snow", longname="Snow vertical layers", &
                       cartesian='Z', sense=-1, ax_size=IG%NkSnow)
    call set_axis_info(extra_axes(5), "band", longname="Frequency and angular band of shortwave radiation", &
                       cartesian='Z', sense=-1, ax_size=4) !### This size is hard-coded for now.
  endif

  call create_dyn_horgrid(dG, G%HI)
  call clone_MOM_domain(G%Domain, dG%Domain)
  call copy_SIS_horgrid_to_dyngrid(G, dG)

  next_var = 1
  do while (next_var <= CS%novars )
    start_var = next_var
    size_in_file = 8*(2*G%Domain%niglobal+2*G%Domain%njglobal+2*nz+1000)

    do m=start_var,CS%novars
      call query_vardesc(CS%restart_field(m)%vars, position=pos, dim_names=dim_names, &
                         caller="SIS_save_restart")
      var_sz = 8
      if (pos /= 0) var_sz = 8*(G%Domain%niglobal+1)*(G%Domain%njglobal+1)
      do n=1,size(dim_names) ; if (len_trim(dim_names(n)) > 0) then
        if (trim(dim_names(n)) == "cat") then ; var_sz = var_sz*IG%CatIce
        elseif (trim(dim_names(n)) == "cat0") then ; var_sz = var_sz*(IG%CatIce+1)
        elseif (trim(dim_names(n)) == "z_ice") then ; var_sz = var_sz*IG%NkIce
        elseif (trim(dim_names(n)) == "z_snow") then ; var_sz = var_sz*IG%NkSnow
        elseif (trim(dim_names(n)) == "band") then ; var_sz = var_sz*4
        endif
      endif ; enddo

      if ((m==start_var) .OR. (size_in_file < max_file_size-var_sz)) then
        size_in_file = size_in_file + var_sz
      else ; exit
      endif

    enddo
    next_var = m

    ! Determine if there is a filename_appendix (used for ensemble runs).
    call get_filename_appendix(filename_appendix)
    if (len_trim(filename_appendix) > 0) then
      length = len_trim(restartname)
      if (restartname(length-2:length) == '.nc') then
        restartname = restartname(1:length-3)//'.'//trim(filename_appendix)//'.nc'
      else
        restartname = restartname(1:length)  //'.'//trim(filename_appendix)
      endif
    endif

    restartpath = trim(directory) // trim(restartname)

    if (num_files < 10) then
      write(suffix,'("_",I1)') num_files
    else
      write(suffix,'("_",I2)') num_files
    endif

    length = len_trim(restartpath)
    if (length < 3) then  ! This case is very uncommon but this test avoids segmentation-faults.
      if (num_files > 0) restartpath = trim(restartpath) // suffix
      restartpath = trim(restartpath)//".nc"
    elseif (restartpath(length-2:length) == ".nc") then
      if (num_files > 0) restartpath = restartpath(1:length-3)//trim(suffix)//".nc"
    else
      if (num_files > 0) restartpath = trim(restartpath) // suffix
      restartpath = trim(restartpath)//".nc"
    endif

    do m=start_var,next_var-1
      vars(m-start_var+1) = CS%restart_field(m)%vars
      call query_vardesc(vars(m), t_grid=t_grid, caller="save_restart")
      t_grid = adjustl(t_grid)
      if (t_grid(1:1) /= 'p') call modify_vardesc(vars(m), t_grid='s', caller="save_restart")
    enddo

    do m=start_var,next_var-1
      conv = CS%restart_field(m)%conv
      ! This gets the full range to do checksums on for arrays with halos.
      ! Note that some sea-ice restart variables are allocated only on the computational grid and
      ! use the coupler's indexing convention, not that of the sea-ice, so offsets are needed.  The
      ! upper bounds are aligned to deal with the haloless symmetric memory case.
      call query_vardesc(vars(m), position=pos, caller="save_restart")
      if (associated(CS%var_ptr4d(m)%p)) then
        call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL, lbound(CS%var_ptr4d(m)%p), ubound(CS%var_ptr4d(m)%p))
        check_val(m-start_var+1,1) = chksum(conv*CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
      elseif (associated(CS%var_ptr3d(m)%p)) then
        call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL, lbound(CS%var_ptr3d(m)%p), ubound(CS%var_ptr3d(m)%p))
        check_val(m-start_var+1,1) = chksum(conv*CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
      elseif (associated(CS%var_ptr2d(m)%p)) then
        call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL, lbound(CS%var_ptr2d(m)%p), ubound(CS%var_ptr2d(m)%p))
        check_val(m-start_var+1,1) = chksum(conv*CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
      elseif (associated(CS%var_ptr1d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(conv*CS%var_ptr1d(m)%p(:)) !?, pelist=(/PE_here()/))
      elseif (associated(CS%var_ptr0d(m)%p)) then
        check_val(m-start_var+1,1) = chksum(conv*CS%var_ptr0d(m)%p, pelist=(/PE_here()/))
      endif
    enddo

    file_thread = SINGLE_FILE ; if (CS%parallel_restartfiles) file_thread = MULTIPLE
    if (present(IG)) then
      call create_MOM_file(IO_handle, trim(restartpath), vars, (next_var-start_var), &
          fields, file_thread, dG=dG, checksums=check_val, extra_axes=extra_axes)
    else
      call create_MOM_file(IO_handle, trim(restartpath), vars, (next_var-start_var), &
          fields, file_thread, dG=dG, checksums=check_val)
    endif

    do m=start_var,next_var-1
      if (associated(CS%var_ptr3d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, CS%var_ptr3d(m)%p, &
                             restart_time, scale=CS%restart_field(m)%conv)
      elseif (associated(CS%var_ptr2d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, CS%var_ptr2d(m)%p, &
                             restart_time, scale=CS%restart_field(m)%conv)
      elseif (associated(CS%var_ptr4d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), G%Domain, CS%var_ptr4d(m)%p, &
                             restart_time, scale=CS%restart_field(m)%conv)
      elseif (associated(CS%var_ptr1d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), CS%var_ptr1d(m)%p, &
                             restart_time, scale=CS%restart_field(m)%conv)
      elseif (associated(CS%var_ptr0d(m)%p)) then
        call MOM_write_field(IO_handle, fields(m-start_var+1), CS%var_ptr0d(m)%p, &
                             restart_time, scale=CS%restart_field(m)%conv)
      endif
    enddo

    call IO_handle%close()

    num_files = num_files+1

  enddo

  call destroy_dyn_horgrid(dG)

end subroutine save_restart

!> restore_SIS_state reads the model state from previously generated files.  All restart
!! variables are read from the first file in the input list in which they are found.
subroutine restore_SIS_state(CS, directory, filelist, G)
  type(SIS_restart_CS),    pointer     :: CS        !< The control structure returned by a previous
                                                    !! call to SIS_restart_init with restart fields
                                                    !! already set by calls to register_restart_field
  character(len=*),        intent(in)  :: directory !< The directory in which to find restart files
  character(len=*),        intent(in)  :: filelist  !< The list of restart file names or just
                                                    !! 'r' or 'F' to read automatically named files
  type(SIS_hor_grid_type), intent(in)  :: G         !< The horizontal grid type

  ! Local variables
  real :: scale  ! A scaling factor for reading a field
  real :: conv   ! The output conversion factor for writing a field
  character(len=200) :: filepath  ! The path (dir/file) to the file being opened.
  character(len=80)  :: fname     ! The name of the current file.
  character(len=8)   :: suffix    ! A suffix (like "_2") that is added to any
                                  ! additional restart files.
  character(len=512) :: mesg      ! A message for warnings.
  character(len=80)  :: varname   ! A variable's name.
  integer :: num_file        ! The number of files (restart files and others
                             ! explicitly in filelist) that are open.
  integer :: i, n, m, missing_fields
  integer :: isL, ieL, jsL, jeL, is0, js0
  integer :: arr_sz(4), ubnd(4)         ! The sizes and upper index bounds of an array
  integer :: nvar ! The number of fields in a file
  integer :: pos  ! The staggering position of a variable

  type(MOM_infra_file) :: IO_handles(CS%max_fields) ! The I/O units of all open files.
  character(len=200) :: unit_path(CS%max_fields) ! The file names.

  character(len=8)   :: hor_grid ! String indicating horizontal grid position
  type(MOM_field), allocatable, dimension(:) :: fields ! Structures with information about each variable
  logical, allocatable, dimension(:) :: global_file  ! True if the file is global
  logical            :: is_there_a_checksum ! Is there a valid checksum that should be checked.
  integer(kind=8)    :: checksum_file  ! The checksum value recorded in the input file.
  integer(kind=8)    :: checksum_data  ! The checksum value for the data that was read in.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "restore_state: Module must be initialized before it is used.")
  if (CS%novars > CS%max_fields) call restart_error(CS)

  ! Get NetCDF ids for all of the restart files.
  if ((len_trim(filelist) == 1) .and. (filelist(1:1) == 'F')) then
    num_file = get_num_restart_files('r', directory, G%domain, CS)
    allocate(global_file(num_file))
    num_file = open_restart_units('r', directory, G%domain, CS, IO_handles=IO_handles, &
                                  file_paths=unit_path, global_files=global_file)
  else
    num_file = get_num_restart_files(filelist, directory, G%domain, CS)
    allocate(global_file(num_file))
    num_file = open_restart_units(filelist, directory, G%domain, CS, IO_handles=IO_handles, &
                                  file_paths=unit_path, global_files=global_file)
  endif

  if (num_file == 0) then
    call SIS_error(FATAL, "SIS_restart: Unable to find any restart files specified by "//&
                   trim(filelist)//"  in directory "//trim(directory))
  endif

  ! Read each variable from the first file in which it is found.
  do n=1,num_file
    call IO_handles(n)%get_file_info(nvar=nvar)

    allocate(fields(nvar))
    call IO_handles(n)%get_file_fields(fields(1:nvar))

    do m=1, nvar
      call IO_handles(n)%get_field_atts(fields(m), name=varname)
      do i=1,CS%num_obsolete_vars
        if (adjustl(lowercase(trim(varname))) == adjustl(lowercase(trim(CS%restart_obsolete(i)%field_name)))) then
            call SIS_error(FATAL, "SIS_restart restore_state: Attempting to use obsolete restart field "//&
                           trim(varname)//" - the new corresponding restart field is "//&
                           trim(CS%restart_obsolete(i)%replacement_name))
        endif
      enddo
    enddo

    missing_fields = 0

    do m=1,CS%novars
      if (CS%restart_field(m)%initialized) cycle
      call query_vardesc(CS%restart_field(m)%vars, hor_grid=hor_grid, caller="restore_state")
      select case (hor_grid)
        case ('q') ; pos = CORNER
        case ('h') ; pos = CENTER
        case ('u') ; pos = EAST_FACE
        case ('v') ; pos = NORTH_FACE
        case ('Bu') ; pos = CORNER
        case ('T')  ; pos = CENTER
        case ('Cu') ; pos = EAST_FACE
        case ('Cv') ; pos = NORTH_FACE
        case ('1') ; pos = 0
        case default ; pos = 0
      end select
      conv = CS%restart_field(m)%conv
      if (conv == 0.0) then ; scale = 1.0 ; else ; scale = 1.0 / conv ; endif

      call SIS_mesg("Attempting to read "//trim(CS%restart_field(m)%var_name), 9)

      do i=1,nvar  ! Loop through the fields that are in the file.
        call IO_handles(n)%get_field_atts(fields(i), name=varname)
        if (lowercase(trim(varname)) == lowercase(trim(CS%restart_field(m)%var_name))) then
          checksum_data = -1
          if (CS%checksum_required) then
            call IO_handles(n)%read_field_chksum(fields(i), checksum_file, is_there_a_checksum)
          else
            checksum_file = -1
            is_there_a_checksum = .false. ! Do not need to do data checksumming.
          endif

          if (associated(CS%var_ptr1d(m)%p))  then
            ! Read a 1d array, which should be invariant to domain decomposition.
            call MOM_read_data(unit_path(n), varname, CS%var_ptr1d(m)%p, timelevel=1, scale=scale, &
                               MOM_Domain=G%Domain, global_file=global_file(n), file_may_be_4d=.true.)
            if (is_there_a_checksum) checksum_data = chksum(conv*CS%var_ptr1d(m)%p(:)) !?, pelist=(/PE_here()/))
          elseif (associated(CS%var_ptr0d(m)%p)) then ! Read a scalar...
            call MOM_read_data(unit_path(n), varname, CS%var_ptr0d(m)%p, timelevel=1, scale=scale, &
                               MOM_Domain=G%Domain, global_file=global_file(n), file_may_be_4d=.true.)
            if (is_there_a_checksum) checksum_data = chksum(conv*CS%var_ptr0d(m)%p, pelist=(/PE_here()/))
          elseif (associated(CS%var_ptr2d(m)%p)) then  ! Read a 2d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr2d(m)%p, G%Domain, timelevel=1, &
                                 position=pos, scale=scale, global_file=global_file(n), file_may_be_4d=.true.)
            else ! This array is not domain-decomposed.  This variant is not yet implemented.
              call SIS_error(FATAL, &
                        "SIS_restart does not support 2-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr2d(m)%p, no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) then
              call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL, &
                                            lbound(CS%var_ptr2d(m)%p), ubound(CS%var_ptr2d(m)%p))
              checksum_data = chksum(conv*CS%var_ptr2d(m)%p(isL:ieL,jsL:jeL))
            endif
          elseif (associated(CS%var_ptr3d(m)%p)) then  ! Read a 3d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, G%Domain, timelevel=1, &
                                 position=pos, scale=scale, global_file=global_file(n), file_may_be_4d=.true.)
            else ! This array is not domain-decomposed.  This variant is not yet implemented.
              call SIS_error(FATAL, &
                        "SIS_restart does not support 3-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr3d(m)%p, no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) then
              call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL, &
                                            lbound(CS%var_ptr3d(m)%p), ubound(CS%var_ptr3d(m)%p))
              checksum_data = chksum(conv*CS%var_ptr3d(m)%p(isL:ieL,jsL:jeL,:))
            endif
          elseif (associated(CS%var_ptr4d(m)%p)) then  ! Read a 4d array.
            if (pos /= 0) then
              call MOM_read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, G%Domain, timelevel=1, &
                                 position=pos, scale=scale, global_file=global_file(n))
            else ! This array is not domain-decomposed.  This variant may be under-tested.
              call SIS_error(FATAL, &
                        "SIS_restart does not support 4-d arrays without domain decomposition.")
              ! call read_data(unit_path(n), varname, CS%var_ptr4d(m)%p, no_domain=.true., timelevel=1)
            endif
            if (is_there_a_checksum) then
              call get_checksum_loop_ranges(G, pos, isL, ieL, jsL, jeL, &
                                            lbound(CS%var_ptr4d(m)%p), ubound(CS%var_ptr4d(m)%p))
              checksum_data = chksum(conv*CS%var_ptr4d(m)%p(isL:ieL,jsL:jeL,:,:))
            endif
          else
            call SIS_error(FATAL, "SIS_restart restore_state: No pointers set for "//trim(varname))
          endif

          if (is_root_pe() .and. is_there_a_checksum .and. (checksum_file /= checksum_data)) then
             write (mesg,'(a,Z16,a,Z16,a)') "Checksum of input field "// trim(varname)//" ",checksum_data,&
                                          " does not match value ", checksum_file, &
                                          " stored in "//trim(unit_path(n)//"." )
             call SIS_error(FATAL, "SIS_restart(restore_state): "//trim(mesg) )
          endif

          CS%restart_field(m)%initialized = .true.
          call SIS_mesg("Done reading "//trim(varname), 9)

          exit ! Start search for next restart variable.
        endif
      enddo
      if (i>nvar) missing_fields = missing_fields+1
    enddo

    deallocate(fields)
    if (missing_fields == 0) exit
  enddo

  do n=1,num_file
    call IO_handles(n)%close()
  enddo

  deallocate(global_file)

  ! Check whether any mandatory fields have not been found.
  CS%restart = .true.
  do m=1,CS%novars
    if (.not.(CS%restart_field(m)%initialized)) then
      CS%restart = .false.
      if (CS%restart_field(m)%mand_var) then
        call SIS_error(FATAL,"SIS_restart: Unable to find mandatory variable " &
                       //trim(CS%restart_field(m)%var_name)//" in restart files.")
      endif
    endif
  enddo

end subroutine restore_SIS_state

!> Deallocate memory associated with a SIS_restart_CS variable.
subroutine SIS_restart_end(CS)
  type(SIS_restart_CS),  pointer    :: CS !< A pointer to a SIS_restart_CS object

  if (associated(CS%restart_field)) deallocate(CS%restart_field)
  if (associated(CS%restart_obsolete)) deallocate(CS%restart_obsolete)
  if (associated(CS%var_ptr0d)) deallocate(CS%var_ptr0d)
  if (associated(CS%var_ptr1d)) deallocate(CS%var_ptr1d)
  if (associated(CS%var_ptr2d)) deallocate(CS%var_ptr2d)
  if (associated(CS%var_ptr3d)) deallocate(CS%var_ptr3d)
  if (associated(CS%var_ptr4d)) deallocate(CS%var_ptr4d)
  deallocate(CS)

end subroutine SIS_restart_end

end module SIS_restart
