module SIS_get_input
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
!*  By Robert Hallberg, July 2013                                      *
!*                                                                     *
!*    The subroutine in this file reads the MOM6 namelist input, which *
!*  indicates which directories to use for certain types of input and  *
!*  output, and where to look for the full parsable input file(s).     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : open_param_file, param_file_type, read_param
use MOM_io, only : file_exists, close_file, slasher
use MOM_io, only : open_namelist_file, check_nml_error

implicit none ; private

public Get_SIS_Input, archaic_nml_check

! This structure is to simplify communication with the calling code.

type, public :: directories
  character(len=120) :: &
    restart_input_dir = ' ',& ! The directory to read restart and input files.
    restart_output_dir = ' ',&! The directory into which to write restart files.
    output_directory = ' ', & ! The directory to use to write the model output.
    input_filename  = ' '     ! A string that indicates the input files or how
                              ! the run segment should be started.
end type directories

interface archaic_nml_check
  module procedure archaic_check_real, archaic_check_int, archaic_check_logical
end interface archaic_nml_check

contains

subroutine Get_SIS_Input(param_file, dirs, check_params)
  type(param_file_type), optional, intent(out) :: param_file
  type(directories),     optional, intent(out) :: dirs
  logical,               optional, intent(in)  :: check_params

!    See if the run is to be started from saved conditions, and get  !
!  the names of the I/O directories and initialization file.  This   !
!  subroutine also calls the subroutine that allows run-time changes !
!  in parameters.                                                    !
  integer, parameter :: npf = 5 ! Maximum number of parameter files
  character(len=120) :: &
    parameter_filename(npf) = ' ', & ! List of files containing parameters.
    output_directory = ' ', &   ! Directory to use to write the model output.
    restart_input_dir = ' ', &  ! Directory for reading restart and input files.
    restart_output_dir = ' ', & ! Directory into which to write restart files.
    input_filename  = ' '       ! A string that indicates the input files or how
                                ! the run segment should be started.
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
    dirs%output_directory = slasher(output_directory)
    dirs%restart_output_dir = slasher(restart_output_dir)
    dirs%restart_input_dir = slasher(restart_input_dir)
    dirs%input_filename = input_filename
  endif

  if (present(param_file)) then
    valid_param_files = 0
    do io = 1, npf
      if (len_trim(trim(parameter_filename(io))) > 0) then
        call open_param_file(trim(parameter_filename(io)), param_file, &
                             check_params, component="SIS")
        valid_param_files = valid_param_files + 1
      endif
    enddo
    if (valid_param_files == 0) call SIS_error(FATAL, "There must be at "//&
         "least 1 valid entry in input_filename in SIS_input_nml in input.nml.")
  endif

end subroutine Get_SIS_Input

subroutine archaic_check_real(param_file, pf_name, nl_name, var, missing, default)
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: pf_name, nl_name
  real,                  intent(in) :: var, missing
  real, optional,        intent(in) :: default
  
  real :: pf_val, def_val
  
  if (var == missing) return
  if (is_root_pe()) then
    def_val = missing ; if (present(default)) def_val = default
    pf_val = def_val ; call read_param(param_file, pf_name, pf_val)

    if (var == pf_val) then
      call SIS_error(WARNING, "Archaic SIS namelist variable "//trim(nl_name)//&
          " appears to be used, \n but is set the same as the SIS_input variable "//&
          trim(pf_name))
    elseif (pf_val == default) then
      call SIS_error(FATAL, "Archaic SIS namelist variable "//trim(nl_name)//&
          " appears to be used. \n Instead, use the SIS_input variable "//&
          trim(pf_name))
    else
      call SIS_error(FATAL, "Archaic SIS namelist variable "//trim(nl_name)//&
          " appears to be used, \n and is set differently from the SIS_input variable "//&
          trim(pf_name))
    endif
  endif
  
end subroutine archaic_check_real

subroutine archaic_check_int(param_file, pf_name, nl_name, var, missing, default)
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: pf_name, nl_name
  integer,               intent(in) :: var, missing
  integer, optional,     intent(in) :: default
  
  integer :: pf_val, def_val
  
  if (var == missing) return
  if (is_root_pe()) then
    def_val = missing ; if (present(default)) def_val = default
    pf_val = def_val ; call read_param(param_file, pf_name, pf_val)

    if (var == pf_val) then
      call SIS_error(WARNING, "Archaic SIS namelist variable "//trim(nl_name)//&
          " appears to be used, \n but is set the same as the SIS_input variable "//&
          trim(pf_name))
    elseif (pf_val == default) then
      call SIS_error(FATAL, "Archaic SIS namelist variable "//trim(nl_name)//&
          " appears to be used. \n Instead, use the SIS_input variable "//&
          trim(pf_name))
    else
      call SIS_error(FATAL, "Archaic SIS namelist variable "//trim(nl_name)//&
          " appears to be used, \n and is set differently from the SIS_input variable "//&
          trim(pf_name))
    endif
  endif
  
end subroutine archaic_check_int

subroutine archaic_check_logical(param_file, pf_name, nl_name, var, missing, default)
  type(param_file_type), intent(in) :: param_file
  character(len=*),      intent(in) :: pf_name, nl_name
  logical,               intent(in) :: var, missing
  logical, optional,     intent(in) :: default
  
  logical :: pf_val, def_val
  
  if (var .eqv. missing) return

  if (is_root_pe()) then
    def_val = missing ; if (present(default)) def_val = default
    pf_val = def_val ; call read_param(param_file, pf_name, pf_val)

    if (var .eqv. pf_val) then
      call SIS_error(WARNING, "Archaic SIS namelist variable "//trim(nl_name)//&
          " appears to be used, \n but is set the same as the SIS_input variable "//&
          trim(pf_name))
    else
      call SIS_error(FATAL, "Archaic SIS namelist variable "//trim(nl_name)//&
          " appears to be used. \n Instead, use the SIS_input variable "//&
          trim(pf_name))
    endif
  endif
  
end subroutine archaic_check_logical

end module SIS_get_input
