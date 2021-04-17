!> The SIS2 facility for reading and writing restart files, and querying what has been read.
module SIS_restart

! This file is part of SIS2. See LICENSE.md for the license.

use coupler_types_mod, only : coupler_type_register_restarts
use fms_io_mod,        only : restart_file_type, FMS1_register_restart=>register_restart_field
use fms_io_mod,        only : save_restart_FMS1=>save_restart, FMS1_restore_state=>restore_state
use fms_io_mod,        only : FMS1_query_initialized=>query_initialized

use ice_grid,          only : ice_grid_type
use SIS_hor_grid,      only : SIS_hor_grid_type

use MOM_coms_infra,    only : SIS_chksum=>field_chksum
use MOM_coupler_types, only : coupler_2d_bc_type, coupler_3d_bc_type
use MOM_domain_infra,  only : MOM_domain_type, domain2D
use MOM_domain_infra,  only : CENTER, CORNER, EAST_FACE, NORTH_FACE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, NOTE, is_root_pe, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, read_param, log_param, log_version, param_file_type
use MOM_io,            only : create_file, file_type, fieldtype, file_exists, open_file, close_file
use MOM_io,            only : get_filename_appendix
use MOM_io,            only : MULTIPLE, READONLY_FILE, SINGLE_FILE
use MOM_string_functions, only : slasher
use MOM_time_manager,  only : time_type

implicit none ; private

public :: restart_file_type, restore_SIS_state, register_restart_field, save_restart, query_inited
public :: query_initialized, SIS_restart_init, SIS_restart_end, only_read_from_restarts
public :: determine_is_new_run, is_new_run

!> A restart registry and the control structure for restarts
type, public :: SIS_restart_CS ; private
  type(restart_file_type), pointer :: fms_restart => NULL() !< The FMS restart file type to use
  type(domain2d), pointer :: mpp_domain => NULL() !< The mpp domain to use for read of decomposed fields.
  character(len=240) :: restartfile !< The name or name root for MOM restart files.
  logical :: parallel_restartfiles = .true. !< If true, each PE writes its own restart file,
                                    !! otherwise they are combined internally.
  logical :: new_run                !< If true, the input filenames and restart file existence will
                                    !! result in a new run that is not initialized from restart files.
  logical :: new_run_set = .false.  !< If true, new_run has been determined for this restart_CS.
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

!> restart_files_exist determines whether any restart files exist.
function restart_files_exist(filename, directory, G, CS)
  character(len=*),        intent(in)  :: filename  !< The list of restart file names or a single
                                                    !! character 'r' to read automatically named files
  character(len=*),        intent(in)  :: directory !< The directory in which to find restart files
  type(SIS_hor_grid_type), intent(in)  :: G         !< The horizontal grid type
  type(SIS_restart_CS),    pointer     :: CS        !< The control structure returned by a previous
                                                    !! call to SIS_restart_init
  logical :: restart_files_exist                    !< The function result, which indicates whether
                                                    !! any of the explicitly or automatically named
                                                    !! restart files exist in directory
  integer :: num_files

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "restart_files_exist: Module must be initialized before it is used.")

  if ((LEN_TRIM(filename) == 1) .and. (filename(1:1) == 'F')) then
    num_files = get_num_restart_files('r', directory, G, CS)
  else
    num_files = get_num_restart_files(filename, directory, G, CS)
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
    CS%new_run = (get_num_restart_files('r', directory, G, CS) == 0)
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
!! them and returns unit ids, paths and whether the files are global or spatially decomposed.
function open_restart_units(filename, directory, G, CS, IO_handles, file_paths, &
                            global_files) result(num_files)
  character(len=*),        intent(in)  :: filename  !< The list of restart file names or a single
                                                  !! character 'r' to read automatically named files
  character(len=*),        intent(in)  :: directory !< The directory in which to find restart files
  type(SIS_hor_grid_type), intent(in)  :: G       !< The horizontal grid type
  type(SIS_restart_CS),    pointer     :: CS      !< The control structure returned by a previous
                                                  !! call to SIS_restart_init
  type(file_type), dimension(:), &
                 optional, intent(out) :: IO_handles !< The I/O handles of all opened files
  character(len=*), dimension(:), &
                 optional, intent(out) :: file_paths !< The full paths to the restart files
  logical, dimension(:), &
               optional, intent(out) :: global_files !< True if a file is global
  integer :: num_files  !< The number of files (both automatically named restart
                        !! files and others explicitly in filename) that have been opened.

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
            call open_file(IO_handles(nf), trim(filepath), READONLY_FILE, &
                           threading=MULTIPLE, fileset=SINGLE_FILE)
          if (present(global_files)) global_files(nf) = .true.
          if (present(file_paths)) file_paths(nf) = filepath
        elseif (CS%parallel_restartfiles) then
          ! Look for decomposed files using the I/O Layout.
          fexists = file_exists(filepath, G%Domain)
          if (fexists) then
            nf = nf + 1
            if (present(IO_handles)) &
              call open_file(IO_handles(nf), trim(filepath), READONLY_FILE, MOM_domain=G%Domain)
            if (present(global_files)) global_files(nf) = .false.
            if (present(file_paths)) file_paths(nf) = filepath
          endif
        endif

        if (fexists) then
          if (is_root_pe() .and. (present(IO_handles))) &
            call SIS_mesg("SIS_restart: MOM run restarted using : "//trim(filepath))
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
          call open_file(IO_handles(nf), trim(filepath), READONLY_FILE, &
                       threading=MULTIPLE, fileset=SINGLE_FILE)
        if (present(global_files)) global_files(nf) = .true.
        if (present(file_paths)) file_paths(nf) = filepath
        if (is_root_pe() .and. (present(IO_handles))) &
          call SIS_mesg("SIS_restart: MOM run restarted using : "//trim(filepath))
      else
        if (present(IO_handles)) &
          call SIS_error(WARNING,"SIS_restart: Unable to find restart file : "//trim(filepath))
      endif

    endif
  enddo ! while (start_char < len_trim(filename)) loop

  write(count_msg, '(I)') nf
  call SIS_mesg("SIS2: open_restart_units found "//trim(count_msg)//" files using "//&
                trim(filename)//" in directory "//trim(directory))

  num_files = nf

end function open_restart_units

!> get_num_restart_files returns the number of existing restart files that match the provided
!! directory structure and other information stored in the control structure and optionally
!! also provides the full paths to these files.
function get_num_restart_files(filenames, directory, G, CS) result(num_files)
  character(len=*),        intent(in)  :: filenames !< The list of restart file names or a single
                                                    !! character 'r' to read automatically named files
  character(len=*),        intent(in)  :: directory !< The directory in which to find restart files
  type(SIS_hor_grid_type), intent(in)  :: G         !< The horizontal grid type
  type(SIS_restart_CS),    pointer     :: CS        !< The control structure returned by a previous
                                                    !! call to SIS_restart_init
  integer :: num_files  !< The function result, the number of files (both automatically named
                        !! restart files and others explicitly in filename) that have been opened

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "get_num_restart_files: Module must be initialized before it is used.")

  ! This call uses open_restart_units without the optional arguments needed to actually
  ! open the files to determine the number of restart files.
  num_files = open_restart_units(filenames, directory, G, CS)

end function get_num_restart_files


!> Initialize and set up a restart control structure.
subroutine SIS_restart_init(CS, filename, domain, param_file)
  type(SIS_restart_CS),  pointer    :: CS !< A pointer to a SIS_restart_CS object that is allocated here
  character(len=*),      intent(in) :: filename !< The path to the restart file.
  type(MOM_domain_type), intent(in) :: domain   !< The MOM domain descriptor being used
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_restart_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  allocate(CS%fms_restart)
  CS%restartfile = trim(filename)
  CS%mpp_domain => domain%mpp_domain

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

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_4d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  idr = FMS1_register_restart(CS%fms_restart, CS%restartfile, name, f_ptr, longname=longname, &
                              units=units, mandatory=mandatory, &
                              domain=CS%mpp_domain, position=position)

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

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_3d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  idr = FMS1_register_restart(CS%fms_restart, CS%restartfile, name, f_ptr, longname=longname, &
                              units=units, mandatory=mandatory, &
                              domain=CS%mpp_domain, position=position)

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

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_2d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  idr = FMS1_register_restart(CS%fms_restart, CS%restartfile, name, f_ptr, longname=longname, &
                              units=units, mandatory=mandatory, &
                              domain=CS%mpp_domain, position=position)

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

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_1d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  idr = FMS1_register_restart(CS%fms_restart, CS%restartfile, name, f_ptr, longname=longname, &
                              units=units, mandatory=mandatory, no_domain=.true.)

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

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_field_0d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(name))

  idr = FMS1_register_restart(CS%fms_restart, CS%restartfile, name, f_ptr, longname=longname, &
                              units=units, mandatory=mandatory, no_domain=.true.)

end subroutine register_restart_field_0d


!> Register a 3-d coupler-type field for restarts.  The coupler type provides the relevant metadata.
subroutine register_restart_coupler_type_3d(CS, bc_ptr, varname_prefix, mandatory)
  type(SIS_restart_CS),       pointer       :: CS     !< A pointer to a SIS_restart_CS object (intent in/out)
  type(coupler_3d_bc_type),   intent(inout) :: bc_ptr !< The coupler type field to be read or written
  character(len=*), optional, intent(in)    :: varname_prefix !< A prefix for the variable name in the
                                                      !! restart file, intended to allow multiple
                                                      !! BC_type variables to use the same restart
                                                      !! files.
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is
                                                      !! not successfully read from the restart file.

  ! Local variables
  character(len=:), allocatable :: bc_name

  bc_name = "unspecified" ; if (present(varname_prefix)) bc_name = trim(varname_prefix)

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_coupler_type_3d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(bc_name))

  call coupler_type_register_restarts(bc_ptr, CS%restartfile, CS%fms_restart, mpp_domain=CS%mpp_domain, &
                                      varname_prefix=varname_prefix)

end subroutine register_restart_coupler_type_3d


!> Register a 2-d coupler-type field for restarts.  The coupler type provides the relevant metadata.
subroutine register_restart_coupler_type_2d(CS, bc_ptr, varname_prefix, mandatory)
  type(SIS_restart_CS),       pointer       :: CS     !< A pointer to a SIS_restart_CS object (intent in/out)
  type(coupler_2d_bc_type),   intent(inout) :: bc_ptr !< The coupler type field to be read or written
  character(len=*), optional, intent(in)    :: varname_prefix !< A prefix for the variable name in the
                                                      !! restart file, intended to allow multiple
                                                      !! BC_type variables to use the same restart
                                                      !! files.
  logical,          optional, intent(in) :: mandatory !< If true, the run will abort if this field is
                                                      !! not successfully read from the restart file.

  ! Local variables
  character(len=:), allocatable :: bc_name

  bc_name = "unspecified" ; if (present(varname_prefix)) bc_name = trim(varname_prefix)

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: register_restart_coupler_type_2d: "//&
      "Restart_CS must be initialized before it is used to register "//trim(bc_name))

  call coupler_type_register_restarts(bc_ptr, CS%restartfile, CS%fms_restart, mpp_domain=CS%mpp_domain, &
                                      varname_prefix=varname_prefix)

end subroutine register_restart_coupler_type_2d

!====================== only_read_from_restarts variants =======================

!> Try to read a named 4-d field from the restart files
subroutine only_read_restart_field_4d(CS, name, f_ptr, domain, position, directory, success)
  type(SIS_restart_CS),            pointer     :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),                intent(in)  :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:,:),target, intent(in)  :: f_ptr     !< A pointer to the field to be read
  type(MOM_domain_type),           intent(in)  :: domain    !< The MOM domain descriptor being used
  integer,               optional, intent(in)  :: position  !< A coded integer indicating the horizontal
                                                            !! position of this variable
  character(len=*),      optional, intent(in)  :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out) :: success   !< True if the field was read successfully

  !For FMS2_io: type(FmsNetcdfDomainFile_t) :: file_ptr
  character(len=:), allocatable :: full_path
  logical :: opened
  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: only_read_restart_field_4d: "//&
      "Restart_CS must be initialized before it is used to read "//trim(name))

  idr = FMS1_register_restart(CS%fms_restart, CS%restartfile, name, f_ptr, read_only=.true., &
                              domain=domain%mpp_domain, mandatory=.false., position=position)
  call FMS1_restore_state(CS%fms_restart, idr, directory=directory)

  if (present(success)) success = query_initialized(CS, name)

end subroutine only_read_restart_field_4d

!> Try to read a named 3-d field from the restart files
subroutine only_read_restart_field_3d(CS, name, f_ptr, domain, position, directory, success)
  type(SIS_restart_CS),            pointer     :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),                intent(in)  :: name      !< The variable name to be used in the restart file
  real, dimension(:,:,:),  target, intent(in)  :: f_ptr     !< A pointer to the field to be read
  type(MOM_domain_type),           intent(in)  :: domain    !< The MOM domain descriptor being used
  integer,               optional, intent(in)  :: position  !< A coded integer indicating the horizontal
                                                            !! position of this variable
  character(len=*),      optional, intent(in)  :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out) :: success   !< True if the field was read successfully

  !For FMS2_io: type(FmsNetcdfDomainFile_t) :: file_ptr
  character(len=:), allocatable :: full_path
  logical :: opened
  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: only_read_restart_field_3d: "//&
      "Restart_CS must be initialized before it is used to read "//trim(name))

  idr = FMS1_register_restart(CS%fms_restart, CS%restartfile, name, f_ptr, read_only=.true., &
                              domain=domain%mpp_domain, mandatory=.false., position=position)
  call FMS1_restore_state(CS%fms_restart, idr, directory=directory)

  if (present(success)) success = query_initialized(CS, name)

end subroutine only_read_restart_field_3d

!> Try to read a named 2-d field from the restart files
subroutine only_read_restart_field_2d(CS, name, f_ptr, domain, position, directory, success)
  type(SIS_restart_CS),            pointer     :: CS        !< A pointer to a SIS_restart_CS object (intent in/out)
  character(len=*),                intent(in)  :: name      !< The variable name to be used in the restart file
  real, dimension(:,:),    target, intent(in)  :: f_ptr     !< A pointer to the field to be read
  type(MOM_domain_type),           intent(in)  :: domain    !< The MOM domain descriptor being used
  integer,               optional, intent(in)  :: position  !< A coded integer indicating the horizontal
                                                            !! position of this variable
  character(len=*),      optional, intent(in)  :: directory !< The directory in which to seek restart files.
  logical,               optional, intent(out) :: success   !< True if the field was read successfully

  !For FMS2_io: type(FmsNetcdfDomainFile_t) :: file_ptr
  character(len=:), allocatable :: full_path
  logical :: opened
  integer :: idr

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart: only_read_restart_field_2d: "//&
      "Restart_CS must be initialized before it is used to read "//trim(name))

  idr = FMS1_register_restart(CS%fms_restart, CS%restartfile, name, f_ptr, read_only=.true., &
                              domain=domain%mpp_domain, mandatory=.false., position=position)
  call FMS1_restore_state(CS%fms_restart, idr, directory=directory)

  if (present(success)) success = query_initialized(CS, name)

end subroutine only_read_restart_field_2d

!> query_initialized determines whether a named field has been successfully read from a restart file yet.
logical function query_initialized(CS, name)
  type(SIS_restart_CS), pointer    :: CS !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried
  !   This subroutine returns .true. if the field referred to by name has
  ! initialized from a restart file, and .false. otherwise.

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "query_initialized: Module must be initialized before it is used.")

  query_initialized = FMS1_query_initialized(CS%fms_restart, name)

end function query_initialized

!> query_inited determines whether a named field has been successfully read from a restart file yet.
!! It is identical to query_initialized, but has a separate name to deal with an unexplained
!! problem that the pgi compiler has with reused function names between modules.
logical function query_inited(CS, name)
  type(SIS_restart_CS), pointer    :: CS !< A pointer to a SIS_restart_CS object (intent in)
  character(len=*),     intent(in) :: name !< The name of the field that is being queried

  if (.not.associated(CS)) call SIS_error(FATAL, "SIS_restart " // &
      "query_inited: Module must be initialized before it is used.")

  query_inited = FMS1_query_initialized(CS%fms_restart, name)

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
  elseif (position == NORTH_FACE) then
    dim_names(1:2) = (/ "iq", "jq" /)
  elseif (position == EAST_FACE) then
    dim_names(1:2) = (/ "iq", "jq" /)
  elseif (present(varname)) then
    call SIS_error(FATAL, "set_axis_names_from_pos: Unrecognized position setting for "//trim(varname))
  else
    call SIS_error(FATAL, "set_axis_names_from_pos: Unrecognized position setting")
  endif
end subroutine axis_names_from_pos


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

  ! Several of the arguments are not needed here, but they will be needed with the new stand-alone SIS2
  ! restart capabilities.
  call save_restart_FMS1(CS%fms_restart, time_stamp=time_stamp, directory=directory)

end subroutine save_restart

!> Restore the entire state of the sea ice or a single varable using a SIS_restart control structure.
subroutine restore_SIS_state(CS, directory, filelist, G)
  type(SIS_restart_CS),    pointer     :: CS        !< The control structure returned by a previous
                                                    !! call to SIS_restart_init with restart fields
                                                    !! already set by calls to register_restart_field
  character(len=*),        intent(in)  :: directory !< The directory in which to find restart files
  character(len=*),        intent(in)  :: filelist  !< The list of restart file names or just
                                                    !! 'r' or 'F' to read automatically named files
  type(SIS_hor_grid_type), intent(in)  :: G         !< The horizontal grid type

  if ((len_trim(filelist) == 1) .and. (filelist(1:1) == 'n')) &
    call SIS_error(FATAL, "restore_SIS_state called for a new run.")

  call FMS1_restore_state(CS%fms_restart, directory=directory)

end subroutine restore_SIS_state

!> Deallocate memory associated with a SIS_restart_CS variable.
subroutine SIS_restart_end(CS)
  type(SIS_restart_CS),  pointer    :: CS !< A pointer to a SIS_restart_CS object

  if (associated(CS%fms_restart)) deallocate(CS%fms_restart)
  deallocate(CS)

end subroutine SIS_restart_end

end module SIS_restart
