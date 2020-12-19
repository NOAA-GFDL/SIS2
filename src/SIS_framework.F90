!> Provides wrapped interfaces for infrastructure calls used by SIS2 that can not be
!! found in the MOM6 framework directory.
module SIS_framework

! This file is part of SIS2. See LICENSE.md for the license.

use fms_io_mod,        only : set_domain, nullify_domain
use fms_io_mod,        only : restart_file_type, register_restart_field
use fms_io_mod,        only : save_restart, restore_state, query_initialized
! use fms2_io_mod,       only : query_initialized=>is_registered_to_restart
use mpp_mod,           only : SIS_chksum=>mpp_chksum
use mpp_domains_mod,   only : domain2D, CORNER, EAST, NORTH
use mpp_domains_mod,   only : get_layout=>mpp_get_layout, get_compute_domain=>mpp_get_compute_domain
use mpp_domains_mod,   only : redistribute_data=>mpp_redistribute
use mpp_domains_mod,   only : broadcast_domain=>mpp_broadcast_domain

use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_safe_alloc,    only : safe_alloc=>safe_alloc_alloc, safe_alloc_ptr
use MOM_file_parser,   only : get_param, read_param, log_param, log_version, param_file_type

implicit none ; private

public :: SIS_chksum, redistribute_data, domain2D, CORNER, EAST, NORTH
public :: set_domain, nullify_domain, get_layout, get_compute_domain, broadcast_domain
public :: restart_file_type, register_restart_field, save_restart, restore_state, query_initialized
public :: SIS_initialize_framework, safe_alloc, safe_alloc_ptr

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

end module SIS_framework
