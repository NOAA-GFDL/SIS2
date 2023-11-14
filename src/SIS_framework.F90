!> Provides wrapped interfaces for infrastructure calls used by SIS2 that can not be
!! found in the MOM6 framework directory.
module SIS_framework

! This file is part of SIS2. See LICENSE.md for the license.

use MOM_coms_infra,    only : SIS_chksum=>field_chksum
use MOM_coupler_types, only : coupler_1d_bc_type, coupler_2d_bc_type, coupler_3d_bc_type
use MOM_coupler_types, only : coupler_type_spawn, coupler_type_initialized, coupler_type_send_data
use MOM_coupler_types, only : coupler_type_copy_data,  coupler_type_redistribute_data
use MOM_coupler_types, only : coupler_type_increment_data, coupler_type_rescale_data
use MOM_coupler_types, only : coupler_type_set_diags, coupler_type_write_chksums
use MOM_domain_infra,  only : MOM_domain_type, domain2D, get_domain_extent
use MOM_domain_infra,  only : global_field, redistribute_data=>redistribute_array, broadcast_domain
use MOM_domain_infra,  only : CENTER, CORNER, EAST=>EAST_FACE, NORTH=>NORTH_FACE, EAST_FACE, NORTH_FACE
use MOM_domain_infra,  only : set_domain, nullify_domain
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, read_param, log_param, log_version, param_file_type
! use MOM_io,            only : MULTIPLE, READONLY_FILE, SINGLE_FILE
use MOM_safe_alloc,    only : safe_alloc=>safe_alloc_alloc, safe_alloc_ptr
use MOM_time_manager,  only : time_type

implicit none ; private

public :: SIS_chksum, redistribute_data, domain2D
public :: set_domain, nullify_domain, get_domain_extent, broadcast_domain
public :: coupler_1d_bc_type, coupler_2d_bc_type, coupler_3d_bc_type
public :: coupler_type_spawn, coupler_type_initialized, coupler_type_send_data
public :: coupler_type_redistribute_data, coupler_type_copy_data, coupler_type_rescale_data
public :: coupler_type_increment_data, coupler_type_write_chksums, coupler_type_set_diags
public :: SIS_initialize_framework, safe_alloc, safe_alloc_ptr
! These encoding constants are used to indicate the discretization position of a variable
public :: CENTER, CORNER, EAST, NORTH, EAST_FACE, NORTH_FACE

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_initialize_framework is a template that might be used later to initialize structures
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

end module SIS_framework
