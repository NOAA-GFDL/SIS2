!> Reads the SIS2 namelist input, which indicates which directories to use for certain types of
!! input and output, and where to look for the full parsable input file(s).
module SIS_get_input

! This file is a part of SIS2. See LICENSE.md for the license.

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : open_param_file, param_file_type, read_param
use MOM_io, only : file_exists, close_file, slasher, ensembler
use MOM_io, only : open_namelist_file, check_nml_error

implicit none ; private

public Get_SIS_Input

! This structure is to simplify communication with the calling code.

!> Container for paths and parameter file names.
type, public :: directories
  character(len=240) :: &
    restart_input_dir = ' ',& !< The directory to read restart and input files.
    restart_output_dir = ' ',&!< The directory into which to write restart files.
    output_directory = ' ', & !< The directory to use to write the model output.
    input_filename  = ' '     !< A string that indicates the input files or how
                              !! the run segment should be started.
end type directories

contains

!> Get_SIS_input reads the SIS namelist entries to see if the run is to be started from
!! a saved restart file, and get the names of the parameter files, I/O directories.
subroutine Get_SIS_Input(param_file, dirs, check_params, component)
  type(param_file_type), optional, intent(out) :: param_file !< A structure to parse for run-time parameters
  type(directories),     optional, intent(out) :: dirs         !< Container for paths and parameter file names.
  logical,               optional, intent(in)  :: check_params !< If present and False will stop error checking for
                                                               !! run-time parameters.
  character(len=*),      optional, intent(in)  :: component    !< An alternate component name, the default is "SIS"

!    See if the run is to be started from saved conditions, and get  !
!  the names of the I/O directories and initialization file.  This   !
!  subroutine also calls the subroutine that allows run-time changes !
!  in parameters.                                                    !
  integer, parameter :: npf = 5 ! Maximum number of parameter files
  character(len=240) :: &
    parameter_filename(npf) = ' ', & ! List of files containing parameters.
    output_directory = ' ', &   ! Directory to use to write the model output.
    restart_input_dir = ' ', &  ! Directory for reading restart and input files.
    restart_output_dir = ' ', & ! Directory into which to write restart files.
    input_filename  = ' '       ! A string that indicates the input files or how
                                ! the run segment should be started.
  character(len=240) :: output_dir
  character(len=24)  :: comp
  integer :: unit, io, ierr, valid_param_files

  namelist /SIS_input_nml/ output_directory, input_filename, parameter_filename, &
                           restart_input_dir, restart_output_dir

  if (file_exists('input.nml')) then
    unit = open_namelist_file(file='input.nml')
  else
    call SIS_error(FATAL,'Required namelist file input.nml does not exist.')
  endif

  ierr=1 ; do while (ierr /= 0)
    read(unit, nml=SIS_input_nml, iostat=io, end=10)
    ierr = check_nml_error(io, 'SIS_input_nml')
  enddo
10 call close_file(unit)
  if (present(dirs)) then
    dirs%output_directory = trim(slasher(ensembler(output_directory)))
    dirs%restart_output_dir = trim(slasher(ensembler(restart_output_dir)))
    dirs%restart_input_dir = trim(slasher(ensembler(restart_input_dir)))
    dirs%input_filename = trim(ensembler(input_filename))
  endif

  comp = "SIS" ; if (present(component)) comp = trim(adjustl(component))

  if (present(param_file)) then
    output_dir = trim(slasher(ensembler(output_directory)))
    valid_param_files = 0
    do io = 1, npf
      if (len_trim(trim(parameter_filename(io))) > 0) then
        call open_param_file(trim(parameter_filename(io)), param_file, &
                             check_params, component=comp, &
                             doc_file_dir=output_dir)
        valid_param_files = valid_param_files + 1
      endif
    enddo
    if (valid_param_files == 0) call SIS_error(FATAL, "There must be at "//&
         "least 1 valid entry in input_filename in SIS_input_nml in input.nml.")
  endif

end subroutine Get_SIS_Input

end module SIS_get_input
