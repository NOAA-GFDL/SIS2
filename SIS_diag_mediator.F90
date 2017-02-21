module SIS_diag_mediator
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
!*    The subroutines here provide convenient wrappers to the fms      *
!*  diag_manager interfaces with additional diagnostic capabilies.     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type

use MOM_coms, only : PE_here
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_safe_alloc, only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions, only : lowercase, uppercase, slasher
use MOM_time_manager, only : time_type

use diag_manager_mod, only : diag_manager_init
use diag_manager_mod, only : send_data, diag_axis_init
use diag_manager_mod, only : register_diag_field_fms=>register_diag_field
use diag_manager_mod, only : register_static_field_fms=>register_static_field

implicit none ; private

public set_SIS_axes_info, post_SIS_data, register_SIS_diag_field, time_type
public safe_alloc_ptr, safe_alloc_alloc
public enable_SIS_averaging, disable_SIS_averaging, query_SIS_averaging_enabled
public SIS_diag_mediator_init, SIS_diag_mediator_end, set_SIS_diag_mediator_grid
public SIS_diag_mediator_close_registration, get_SIS_diag_time_end
public diag_axis_init, register_static_field

interface post_SIS_data
  module procedure post_data_2d, post_data_3d
end interface post_SIS_data

! 2D/3D axes type to contain 1D axes handles and pointers to masks
type, public :: axesType
  character(len=15) :: id ! This is the id string for this particular combination of handles
  integer :: rank ! The number of dimensions in the list of axes
  integer, dimension(:), allocatable :: handles ! Handles to 1D axes
  type(SIS_diag_ctrl), pointer :: diag_cs => null()
end type axesType

! This type is used to represent a diagnostic at the diag_mediator level.
type, private :: diag_type
  logical :: in_use
  integer :: fms_diag_id         ! underlying FMS diag id
  character(len=24) :: name
  real, pointer, dimension(:,:)   :: mask2d => null()
  real, pointer, dimension(:,:)   :: mask2d_comp => null()
  real, pointer, dimension(:,:,:) :: mask3d => null()
end type diag_type

!   The following data type contains pointers to diagnostic fields that might
! be shared between modules, and also to the variables that control the handling
! of model output.
type, public :: SIS_diag_ctrl
  integer :: doc_unit = -1 ! The unit number of a diagnostic documentation file.
                           ! This file is open if doc_unit is > 0.

! The following fields are used for the output of the data.
! These give the computational-domain sizes, and are relative to a start value
! of 1 in memory for the tracer-point arrays.
  integer :: is, ie, js, je
! These give the memory-domain sizes, and can be start at any value on each PE.
  integer :: isd, ied, jsd, jed
  real :: time_int              ! The time interval in s for any fields
                                ! that are offered for averaging.
  type(time_type) :: time_end   ! The end time of the valid
                                ! interval for any offered field.
  logical :: ave_enabled = .false. ! .true. if averaging is enabled.

  ! The following are axis types defined for output.
  type(axesType) :: axesBL, axesTL, axesCuL, axesCvL
  type(axesType) :: axesBi, axesTi, axesCui, axesCvi
  type(axesType) :: axesBc, axesTc, axesCuc, axesCvc
  type(axesType) :: axesBc0, axesTc0, axesCuc0, axesCvc0
  type(axesType) :: axesB1, axesT1, axesCu1, axesCv1
  type(axesType) :: axeszi, axeszL
  ! Mask arrays for diagnostics
  real, dimension(:,:),   pointer :: mask2dT   => null()
  real, dimension(:,:),   pointer :: mask2dBu  => null()
  real, dimension(:,:),   pointer :: mask2dCu  => null()
  real, dimension(:,:),   pointer :: mask2dCv  => null()
  real, dimension(:,:,:), pointer :: mask3dTL  => null()
  real, dimension(:,:,:), pointer :: mask3dBuL => null()
  real, dimension(:,:,:), pointer :: mask3dCuL => null()
  real, dimension(:,:,:), pointer :: mask3dCvL => null()
  real, dimension(:,:,:), pointer :: mask3dTi  => null()
  real, dimension(:,:,:), pointer :: mask3dBui => null()
  real, dimension(:,:,:), pointer :: mask3dCui => null()
  real, dimension(:,:,:), pointer :: mask3dCvi => null()
  real, dimension(:,:,:), pointer :: mask3dTC  => null()
  real, dimension(:,:,:), pointer :: mask3dBuC => null()
  real, dimension(:,:,:), pointer :: mask3dCuC => null()
  real, dimension(:,:,:), pointer :: mask3dCvC => null()
  ! Computational domain mask arrays for diagnostics.
  real, dimension(:,:),   pointer :: mask2dT_comp => null()

#define DIAG_ALLOC_CHUNK_SIZE 15
  type(diag_type), dimension(:), allocatable :: diags
  integer :: next_free_diag_id

  !default missing value to be sent to ALL diagnostics registerations
  real :: missing_value = -1.0e34

end type SIS_diag_ctrl

contains

subroutine set_SIS_axes_info(G, IG, param_file, diag_cs, set_vertical, axes_set_name)
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(inout) :: IG
  type(param_file_type),   intent(in)    :: param_file
  type(SIS_diag_ctrl),     intent(inout) :: diag_cs
  logical,          optional, intent(in)    :: set_vertical
  character(len=*), optional, intent(in)    :: axes_set_name
!   This subroutine sets up the grid and axis information for use by SIS.
!
! Arguments: G - The ocean's grid structure.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (inout)   diag_cs - A structure that is used to regulate diagnostic output.
!  (in,opt)  set_vertical - If true (or missing), set up the vertical axes.
!  (in,opt)  axes_set_name - A name to use for this set of axes.  The default is "ice".
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh, id_ct, id_ct0
  integer :: id_xhe, id_yhe
  integer :: k
  real :: zlev_ice(IG%NkIce), zinter_ice(IG%NkIce+1)
  logical :: set_vert, Cartesian_grid
  character(len=80) :: grid_config, units_temp, set_name
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "SIS_diag_mediator" ! This module's name.

  set_vert = .true. ; if (present(set_vertical)) set_vert = set_vertical
  set_name = "ice" ; if (present(axes_set_name)) set_name = trim(axes_set_name)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "GRID_CONFIG", grid_config, &
                 "The method for defining the horizontal grid.  Valid \n"//&
                 "entries include:\n"//&
                 "\t file - read the grid from GRID_FILE \n"//&
                 "\t mosaic - read the grid from a mosaic grid file \n"//&
                 "\t cartesian - a Cartesian grid \n"//&
                 "\t spherical - a spherical grid \n"//&
                 "\t mercator  - a Mercator grid", fail_if_missing=.true.)

  G%x_axis_units = "degrees_E"
  G%y_axis_units = "degrees_N"
  if (index(lowercase(trim(grid_config)),"cartesian") > 0) then
    ! This is a cartesian grid, and may have different axis units.
    Cartesian_grid = .true.
    call get_param(param_file, mod, "AXIS_UNITS", units_temp, &
                 "The units for the x- and y- axis labels.  AXIS_UNITS \n"//&
                 "should be defined as 'k' for km, 'm' for m, or 'd' \n"//&
                 "for degrees of latitude and longitude (the default). \n"//&
                 "Except on a Cartesian grid, only degrees are currently \n"//&
                 "implemented.", default='degrees')
    if (units_temp(1:1) == 'k') then
      G%x_axis_units = "kilometers" ; G%y_axis_units = "kilometers"
    elseif (units_temp(1:1) == 'm') then
      G%x_axis_units = "meters" ; G%y_axis_units = "meters"
    endif
    call log_param(param_file, mod, "explicit AXIS_UNITS", G%x_axis_units)
  else
    Cartesian_grid = .false.
  endif

  id_xq = diag_axis_init('xB', G%gridLonB(G%isgB:G%iegB), G%x_axis_units, 'x', &
            'Boundary point nominal longitude',set_name=set_name, &
             Domain2=G%Domain%mpp_domain)
  id_yq = diag_axis_init('yB', G%gridLatB(G%jsgB:G%jegB), G%y_axis_units, 'y', &
            'Boundary point nominal latitude', set_name=set_name, &
             Domain2=G%Domain%mpp_domain)

  id_xhe = diag_axis_init('xTe', G%gridLonB(G%isg-1:G%ieg), G%x_axis_units, 'x', &
            'T-cell edge nominal longitude', set_name=set_name, &
             Domain2=G%Domain%mpp_domain)
  id_yhe = diag_axis_init('yTe', G%gridLatB(G%jsg-1:G%jeg), G%y_axis_units, 'y', &
            'T-cell edge nominal latitude', set_name=set_name, &
            Domain2=G%Domain%mpp_domain)
  id_xh = diag_axis_init('xT', G%gridLonT(G%isg:G%ieg), G%x_axis_units, 'x', &
              'T point nominal longitude', set_name=set_name, edges=id_xhe, &
              Domain2=G%Domain%mpp_domain)
  id_yh = diag_axis_init('yT', G%gridLatT(G%jsg:G%jeg), G%y_axis_units, 'y', &
              'T point nominal latitude', set_name=set_name, edges=id_yhe, &
              Domain2=G%Domain%mpp_domain)

  if (set_vert) then
    do k=1,IG%NkIce+1 ; zinter_ice(k) = real(k-1) / real(IG%NkIce) ; enddo
    do k=1,IG%NkIce ; zlev_ice(k) = (k-0.5) / real(IG%NkIce) ; enddo
    id_zl = diag_axis_init('zl', zlev_ice, 'layer', 'z', 'Cell depth', &
                           set_name=set_name)
    id_zi = diag_axis_init('zi', zinter_ice, 'interface', 'z', &
                           'Cell interface depth', set_name=set_name)
  else
    id_zl = -1 ; id_zi = -1
  endif

  id_ct = diag_axis_init('ct', IG%cat_thick_lim(1:IG%CatIce), 'meters', 'n', & ! 'z',?
                         'Ice thickness category bounds', set_name=set_name)
  id_ct0 = diag_axis_init('ctu', IG%cat_thick_lim(1:IG%CatIce+1), 'meters', 'n', & ! 'z',?
                         'Ice thickness category upper bounds', set_name=set_name)

  ! Note that there are no 4-d spatial axis groupings yet.  Ferret only started
  ! allowing for 5-d data with version 6.8, which is later than the default for
  ! GFDL.  Once more recent versions come into widespread use, 4-d spatial grids
  ! should be reconsidered.  (R. Hallberg, 8/27/2013)

  ! Vertical axes for the interfaces and layers.
  call defineAxes(diag_cs, (/ id_zi /), diag_cs%axesZi)
  call defineAxes(diag_cs, (/ id_zL /), diag_cs%axesZL)

  ! Axis groupings for the model layers.
  call defineAxes(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%axesTL)
  call defineAxes(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%axesBL)
  call defineAxes(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%axesCuL)
  call defineAxes(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%axesCvL)

  ! Axis groupings for the model interfaces.
  call defineAxes(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%axesTi)
  call defineAxes(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%axesCui)
  call defineAxes(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%axesCvi)
  call defineAxes(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%axesBi)

  ! Axis groupings for the ice thickness categories.
  call defineAxes(diag_cs, (/ id_xh, id_yh, id_ct /), diag_cs%axesTc)
  call defineAxes(diag_cs, (/ id_xq, id_yh, id_ct /), diag_cs%axesCuc)
  call defineAxes(diag_cs, (/ id_xh, id_yq, id_ct /), diag_cs%axesCvc)
  call defineAxes(diag_cs, (/ id_xq, id_yq, id_ct /), diag_cs%axesBc)

  ! Axis groupings for the ocean and ice thickness categories.
  call defineAxes(diag_cs, (/ id_xh, id_yh, id_ct0 /), diag_cs%axesTc0)
  call defineAxes(diag_cs, (/ id_xq, id_yh, id_ct0 /), diag_cs%axesCuc0)
  call defineAxes(diag_cs, (/ id_xh, id_yq, id_ct0 /), diag_cs%axesCvc0)
  call defineAxes(diag_cs, (/ id_xq, id_yq, id_ct0 /), diag_cs%axesBc0)

  ! Axis groupings for 2-D arrays.
  call defineAxes(diag_cs, (/ id_xh, id_yh /), diag_cs%axesT1)
  call defineAxes(diag_cs, (/ id_xq, id_yq /), diag_cs%axesB1)
  call defineAxes(diag_cs, (/ id_xq, id_yh /), diag_cs%axesCu1)
  call defineAxes(diag_cs, (/ id_xh, id_yq /), diag_cs%axesCv1)

end subroutine set_SIS_axes_info

subroutine defineAxes(diag_cs, handles, axes)
  ! Defines "axes" from list of handle and associates mask
  type(SIS_diag_ctrl), target, intent(in)  :: diag_cs
  integer, dimension(:),       intent(in)  :: handles
  type(axesType),              intent(out) :: axes
  ! Local variables
  integer :: n
  n = size(handles)
  if (n<1 .or. n>3) call SIS_error(FATAL,"defineAxes: wrong size for list of handles!")
  allocate( axes%handles(n) )
  axes%id = i2s(handles, n) ! Identifying string
  axes%rank = n
  axes%handles(:) = handles(:)
  axes%diag_cs => diag_cs ! A [circular] link back to the SIS_diag_ctrl structure
end subroutine defineAxes

subroutine set_SIS_diag_mediator_grid(G, diag_cs)
  type(SIS_hor_grid_type), intent(inout) :: G
  type(SIS_diag_ctrl),     intent(inout) :: diag_cs
! Arguments: G - The ocean's grid structure.
!  (inout)   diag_cs - A structure that is used to regulate diagnostic output.
  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied ; diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed
end subroutine set_SIS_diag_mediator_grid

subroutine post_data_2d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id
  real,              intent(in) :: field(:,:)
  type(SIS_diag_ctrl), target, intent(in) :: diag_cs
  logical, optional, intent(in) :: is_static
  logical, optional, intent(in) :: mask(:,:)
! Arguments: diag_field_id - the id for an output variable returned by a
!                            previous call to register_diag_field.
!  (in)      field - The 2-d array being offered for output or averaging.
!  (inout)   diag_cs - A structure that is used to regulate diagnostic output.
!  (in,opt)  is_static - If true, this is a static field that is always offered.
!  (in,opt)  mask - If present, use this real array as the data mask.
  logical :: used, is_stat
  logical :: i_data, j_data
  integer :: isv, iev, jsv, jev
  integer :: fms_diag_id
  type(diag_type), pointer :: diag

  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Get a pointer to the SIS diag type for this id, and the FMS-level diag id.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_2d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  fms_diag_id = diag%fms_diag_id

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  if ( size(field,1) == diag_cs%ied-diag_cs%isd +1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie ; i_data = .true.   ! Data domain
  elseif ( size(field,1) == diag_cs%ied-diag_cs%isd +2 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1 ; i_data = .true. ! Symmetric data domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +1 ) then
    isv = 1 ; iev = diag_cs%ie + 1-diag_cs%is ; i_data = .false. ! Computational domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +2 ) then
    isv = 1 ; iev = diag_cs%ie + 2-diag_cs%is ; i_data = .false. ! Symmetric computational domain
  else
    call SIS_error(FATAL,"post_SIS_data_2d: peculiar size in i-direction of "//trim(diag%name))
  endif
  if ( size(field,2) == diag_cs%jed-diag_cs%jsd +1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je ; j_data = .true.   ! Data domain
  elseif ( size(field,2) == diag_cs%jed-diag_cs%jsd +2 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1 ; j_data = .true. ! Symmetric data domain
  elseif ( size(field,2) == diag_cs%je-diag_cs%js +1 ) then
    jsv = 1 ; jev = diag_cs%je + 1-diag_cs%js ; j_data = .false. ! Computational domain
  elseif ( size(field,1) == diag_cs%je-diag_cs%js +2 ) then
    jsv = 1 ; jev = diag_cs%je + 2-diag_cs%js ; j_data = .false. ! Symmetric computational domain
  else
    call SIS_error(FATAL,"post_SIS_data_2d: peculiar size in j-direction "//trim(diag%name))
  endif

  ! Handle cases where the data and computational domain are the same size.
  if (diag_cs%ied-diag_cs%isd == diag_cs%ie-diag_cs%is) i_data = j_data
  if (diag_cs%jed-diag_cs%jsd == diag_cs%je-diag_cs%js) j_data = i_data

  if (present(mask)) then
    if ((size(field,1) /= size(mask,1)) .or. &
        (size(field,2) /= size(mask,2))) then
      call SIS_error(FATAL, "post_SIS_data_2d: post_SIS_data called with a mask "//&
                            "that does not match the size of field "//trim(diag%name))
    endif
  elseif ( i_data .NEQV. j_data ) then
    call SIS_error(FATAL, "post_SIS_data_2d: post_SIS_data called for "//&
                   trim(diag%name)//" with mixed computational and data domain array sizes.")
  endif

  if (is_stat) then
    if (present(mask)) then
      used = send_data(fms_diag_id, field, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, mask=mask)
    elseif(i_data .and. associated(diag%mask2d)) then
      used = send_data(fms_diag_id, field, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%mask2d)
    elseif((.not.i_data) .and. associated(diag%mask2d_comp)) then
      used = send_data(fms_diag_id, field, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%mask2d_comp)
    else
      used = send_data(fms_diag_id, field, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev)
    endif
  elseif (diag_cs%ave_enabled) then
    if (present(mask)) then
      used = send_data(fms_diag_id, field, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, mask=mask)
    elseif(i_data .and. associated(diag%mask2d)) then
      used = send_data(fms_diag_id, field, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, rmask=diag%mask2d)
    elseif((.not.i_data) .and. associated(diag%mask2d_comp)) then
      used = send_data(fms_diag_id, field, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, rmask=diag%mask2d_comp)
    else
      used = send_data(fms_diag_id, field, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int)
    endif
  endif

end subroutine post_data_2d

subroutine post_data_3d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id
  real,              intent(in) :: field(:,:,:)
  type(SIS_diag_ctrl), target,  intent(in) :: diag_cs
  logical, optional, intent(in) :: is_static
  logical, optional, intent(in) :: mask(:,:,:)
! Arguments: diag_field_id - the id for an output variable returned by a
!                            previous call to register_diag_field.
!  (in)      field - The 3-d array being offered for output or averaging.
!  (inout)   diag_cs - A structure that is used to regulate diagnostic output.
!  (in)      static - If true, this is a static field that is always offered.
!  (in,opt)  mask - If present, use this real array as the data mask.
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  integer :: isv, iev, jsv, jev
  integer :: fms_diag_id
  type(diag_type), pointer :: diag

  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Get a pointer to the SIS diag type for this id, and the FMS-level diag id.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_3d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  fms_diag_id = diag%fms_diag_id

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  if ( size(field,1) == diag_cs%ied-diag_cs%isd +1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie        ! Data domain
  elseif ( size(field,1) == diag_cs%ied-diag_cs%isd +2 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1      ! Symmetric data domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +1 ) then
    isv = 1 ; iev = diag_cs%ie + 1-diag_cs%is  ! Computational domain
  elseif ( size(field,1) == diag_cs%ie-diag_cs%is +2 ) then
    isv = 1 ; iev = diag_cs%ie + 2-diag_cs%is  ! Symmetric computational domain
  else
    call SIS_error(FATAL,"post_SIS_data_3d: peculiar size in i-direction")
  endif
  if ( size(field,2) == diag_cs%jed-diag_cs%jsd +1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je        ! Data domain
  elseif ( size(field,2) == diag_cs%jed-diag_cs%jsd +2 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1      ! Symmetric data domain
  elseif ( size(field,2) == diag_cs%je-diag_cs%js +1 ) then
    jsv = 1 ; jev = diag_cs%je + 1-diag_cs%js  ! Computational domain
  elseif ( size(field,1) == diag_cs%je-diag_cs%js +2 ) then
    jsv = 1 ; jev = diag_cs%je + 2-diag_cs%js  ! Symmetric computational domain
  else
    call SIS_error(FATAL,"post_SIS_data_3d: peculiar size in j-direction")
  endif

  if (present(mask)) then
    if ((size(field,1) /= size(mask,1)) .or. &
        (size(field,2) /= size(mask,2)) .or. &
        (size(field,3) /= size(mask,3))) then
      call SIS_error(FATAL, "post_SIS_data_3d: post_SIS_data called with a mask "//&
                             "that does not match the size of field.")
    endif
  endif

  if (is_stat) then
    if (present(mask)) then
      used = send_data(fms_diag_id, field, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, mask=mask)
    elseif(associated(diag%mask3d)) then
      used = send_data(fms_diag_id, field, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, rmask=diag%mask3d)
    else
      used = send_data(fms_diag_id, field, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev)
    endif
  elseif (diag_cs%ave_enabled) then
    if (present(mask)) then
      used = send_data(fms_diag_id, field, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, mask=mask)
    elseif(associated(diag%mask3d)) then
      used = send_data(fms_diag_id, field, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int, rmask=diag%mask3d)
    else
      used = send_data(fms_diag_id, field, diag_cs%time_end, &
                       is_in=isv, js_in=jsv, ie_in=iev, je_in=jev, &
                       weight=diag_cs%time_int)
    endif
  endif

end subroutine post_data_3d


subroutine enable_SIS_averaging(time_int_in, time_end_in, diag_cs)
  real, intent(in) :: time_int_in
  type(time_type), intent(in) :: time_end_in
  type(SIS_diag_ctrl), intent(inout) :: diag_cs
! This subroutine enables the accumulation of time averages over the
! specified time interval.

! Arguments: time_int_in - the time interval in s over which any
!                          values that are offered are valid.
!  (in)      time_end_in - the end time in s of the valid interval.
!  (inout)   diag_cs - A structure that is used to regulate diagnostic output.

!  if (num_file==0) return
  diag_cs%time_int = time_int_in
  diag_cs%time_end = time_end_in
  diag_cs%ave_enabled = .true.
end subroutine enable_SIS_averaging

! Call this subroutine to avoid averaging any offered fields.
subroutine disable_SIS_averaging(diag_cs)
  type(SIS_diag_ctrl), intent(inout) :: diag_cs
! Argument: diag_cs - A structure that is used to regulate diagnostic output.

  diag_cs%time_int = 0.0
  diag_cs%ave_enabled = .false.

end subroutine disable_SIS_averaging

! Call this subroutine to determine whether the averaging is
! currently enabled.  .true. is returned if it is.
function query_SIS_averaging_enabled(diag_cs, time_int, time_end)
  type(SIS_diag_ctrl),           intent(in)  :: diag_cs
  real,            optional, intent(out) :: time_int
  type(time_type), optional, intent(out) :: time_end
  logical :: query_SIS_averaging_enabled
! Arguments: diag - A structure that is used to regulate diagnostic output.
!  (out,opt) time_int - The current setting of diag_cs%time_int, in s.
!  (out,opt) time_end - The current setting of diag_cs%time_end.

  if (present(time_int)) time_int = diag_cs%time_int
  if (present(time_end)) time_end = diag_cs%time_end
  query_SIS_averaging_enabled = diag_cs%ave_enabled
end function query_SIS_averaging_enabled

function get_SIS_diag_time_end(diag_cs)
  type(SIS_diag_ctrl),           intent(in)  :: diag_cs
  type(time_type) :: get_SIS_diag_time_end
! Argument: diag_cs - A structure that is used to regulate diagnostic output.

!   This function returns the valid end time for diagnostics that are handled
! outside of the MOM6 infrastructure, such as via the generic tracer code.

  get_SIS_diag_time_end = diag_cs%time_end
end function get_SIS_diag_time_end

function register_SIS_diag_field(module_name, field_name, axes, init_time, &
     long_name, units, missing_value, range, mask_variant, standard_name, &
     verbose, do_not_log, err_msg, interp_method, tile_count) result (register_diag_field)
  integer :: register_diag_field
  character(len=*), intent(in) :: module_name, field_name
  type(axesType),   intent(in) :: axes
  type(time_type),  intent(in) :: init_time
  character(len=*), optional, intent(in) :: long_name, units, standard_name
  real,             optional, intent(in) :: missing_value, range(2)
  logical,          optional, intent(in) :: mask_variant, verbose, do_not_log
  character(len=*), optional, intent(out):: err_msg
  character(len=*), optional, intent(in) :: interp_method
  integer,          optional, intent(in) :: tile_count
! Output:    An integer handle for a diagnostic array.
! Arguments: module_name - The name of this module, usually "ocean_model" or "ice_shelf_model".
!  (in)      field_name - The name of the diagnostic field.
!  (in)      axes - A type that indicates the axes for this field.
!  (in)      init_time - The time at which a field is first available?
!  (in,opt)  long_name - The long name of a field.
!  (in,opt)  units - The units of a field.
!  (in,opt)  standard_name - The standardized name associated with a field. (Not yet used in MOM.)
!  (in,opt)  missing_value - A value that indicates missing values.
!  (in,opt)  range - The valid range of a variable. (Not used in MOM.)
!  (in,opt)  mask_variant - If true a logical mask must be provided with post_data calls.  (Not used in MOM.)
!  (in,opt)  verbose - If true, FMS is verbosed. (Not used in MOM.)
!  (in,opt)  do_not_log - If true, do not log something. (Not used in MOM.)
!  (out,opt) err_msg - An character string into which an error message might be placed. (Not used in MOM.)
!  (in,opt)  interp_method - No clue. (Not used in MOM.)
!  (in,opt)  tile_count - No clue. (Not used in MOM.)
  character(len=240) :: mesg
  real :: MOM_missing_value
  integer :: primary_id, fms_id
  type(SIS_diag_ctrl), pointer :: diag_cs
  type(diag_type), pointer :: diag

  MOM_missing_value = axes%diag_cs%missing_value
  if(present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  primary_id = -1

  fms_id = register_diag_field_fms(module_name, field_name, axes%handles, &
         init_time, long_name=long_name, units=units, missing_value=MOM_missing_value, &
         range=range, mask_variant=mask_variant, standard_name=standard_name, &
         verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
         interp_method=interp_method, tile_count=tile_count)
  if (fms_id > 0) then
    primary_id = get_new_diag_id(diag_cs)
    diag => diag_cs%diags(primary_id)
    diag%fms_diag_id = fms_id
    if (len(field_name) > len(diag%name)) then
      diag%name = field_name(1:len(diag%name))
    else ; diag%name = field_name ; endif
  endif

  if (is_root_pe() .and. diag_CS%doc_unit > 0) then
    if (primary_id > 0) then
       mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Used]'
    else
       mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Unused]'
    endif
    write(diag_CS%doc_unit, '(a)') trim(mesg)
    if (present(long_name)) call describe_option("long_name", long_name, diag_CS)
    if (present(units)) call describe_option("units", units, diag_CS)
    if (present(standard_name)) &
      call describe_option("standard_name", standard_name, diag_CS)
  endif

  !Decide what mask to use based on the axes info
  if (primary_id > 0) then
  !3d masks
    if (axes%rank == 3) then
      diag%mask2d => null() ; diag%mask2d_comp => null() ; diag%mask3d => null()
      if (axes%id == diag_cs%axesTL%id) then
        diag%mask3d =>  diag_cs%mask3dTL
      elseif (axes%id == diag_cs%axesBL%id) then
        diag%mask3d =>  diag_cs%mask3dBuL
      elseif (axes%id == diag_cs%axesCuL%id ) then
        diag%mask3d =>  diag_cs%mask3dCuL
      elseif (axes%id == diag_cs%axesCvL%id) then
        diag%mask3d =>  diag_cs%mask3dCvL
      elseif (axes%id == diag_cs%axesTi%id) then
        diag%mask3d =>  diag_cs%mask3dTi
      elseif (axes%id == diag_cs%axesBi%id) then
        diag%mask3d =>  diag_cs%mask3dBui
      elseif (axes%id == diag_cs%axesCui%id ) then
        diag%mask3d =>  diag_cs%mask3dCui
      elseif (axes%id == diag_cs%axesCvi%id) then
        diag%mask3d =>  diag_cs%mask3dCvi
      elseif (axes%id == diag_cs%axesTc%id) then
        diag%mask3d =>  diag_cs%mask3dTC(:,:,1:)
      elseif (axes%id == diag_cs%axesBc%id) then
        diag%mask3d =>  diag_cs%mask3dBuC(:,:,1:)
      elseif (axes%id == diag_cs%axesCuc%id ) then
        diag%mask3d =>  diag_cs%mask3dCuC(:,:,1:)
      elseif (axes%id == diag_cs%axesCvc%id) then
        diag%mask3d =>  diag_cs%mask3dCvC(:,:,1:)
      elseif (axes%id == diag_cs%axesTc0%id) then
        diag%mask3d =>  diag_cs%mask3dTC(:,:,0:)
      elseif (axes%id == diag_cs%axesBc0%id) then
        diag%mask3d =>  diag_cs%mask3dBuC(:,:,0:)
      elseif (axes%id == diag_cs%axesCuc0%id ) then
        diag%mask3d =>  diag_cs%mask3dCuC(:,:,0:)
      elseif (axes%id == diag_cs%axesCvc0%id) then
        diag%mask3d =>  diag_cs%mask3dCvC(:,:,0:)
  !   else
  !       call SIS_error(FATAL, "SIS_diag_mediator:register_diag_field: " // &
  !            "unknown axes for diagnostic variable "//trim(field_name))
      endif
    !2d masks
    elseif (axes%rank == 2) then
      diag%mask2d => null() ; diag%mask2d_comp => null() ; diag%mask3d => null()
      if (axes%id == diag_cs%axesT1%id) then
        diag%mask2d =>  diag_cs%mask2dT
        diag%mask2d_comp => diag_cs%mask2dT_comp
      elseif (axes%id == diag_cs%axesB1%id) then
        diag%mask2d =>  diag_cs%mask2dBu
      elseif (axes%id == diag_cs%axesCu1%id) then
        diag%mask2d =>  diag_cs%mask2dCu
      elseif (axes%id == diag_cs%axesCv1%id) then
        diag%mask2d =>  diag_cs%mask2dCv
  !   else
  !       call SIS_error(FATAL, "SIS_diag_mediator:register_diag_field: " // &
  !            "unknown axes for diagnostic variable "//trim(field_name))
      endif
    else
      call SIS_error(FATAL, "SIS_diag_mediator:register_diag_field: " // &
           "unknown axes for diagnostic variable "//trim(field_name))
    endif
  endif ! if (primary_id>-1)

  register_diag_field = primary_id

end function register_SIS_diag_field

function register_static_field(module_name, field_name, axes, &
     long_name, units, missing_value, range, mask_variant, standard_name, &
     do_not_log, interp_method, tile_count)
  integer :: register_static_field
  character(len=*), intent(in) :: module_name, field_name
  type(axesType),   intent(in) :: axes
  character(len=*), optional, intent(in) :: long_name, units, standard_name
  real,             optional, intent(in) :: missing_value, range(2)
  logical,          optional, intent(in) :: mask_variant, do_not_log
  character(len=*), optional, intent(in) :: interp_method
  integer,          optional, intent(in) :: tile_count
  ! Output:    An integer handle for a diagnostic array.
  ! Arguments: module_name - The name of this module, usually "ocean_model" or "ice_shelf_model".
  !  (in)      field_name - The name of the diagnostic field.
  !  (in)      axes - A container with up to 3 integer handles that indicates the axes for this field.
  !  (in,opt)  long_name - The long name of a field.
  !  (in,opt)  units - The units of a field.
  !  (in,opt)  standard_name - The standardized name associated with a field. (Not yet used in MOM.)
  !  (in,opt)  missing_value - A value that indicates missing values.
  !  (in,opt)  range - The valid range of a variable. (Not used in MOM.)
  !  (in,opt)  mask_variant - If true a logical mask must be provided with post_data calls.  (Not used in MOM.)
  !  (in,opt)  do_not_log - If true, do not log something. (Not used in MOM.)
  !  (in,opt)  interp_method - No clue. (Not used in MOM.)
  !  (in,opt)  tile_count - No clue. (Not used in MOM.)
  character(len=240) :: mesg
  real :: MOM_missing_value
  integer :: primary_id, fms_id
  type(SIS_diag_ctrl), pointer :: diag_cs

  MOM_missing_value = axes%diag_cs%missing_value
  if(present(missing_value)) MOM_missing_value = missing_value

  diag_cs => axes%diag_cs
  primary_id = -1

  fms_id = register_static_field_fms(module_name, field_name, axes%handles, &
       long_name=long_name, units=units, missing_value=MOM_missing_value, &
       range=range, mask_variant=mask_variant, standard_name=standard_name, &
       do_not_log=do_not_log, &
       interp_method=interp_method, tile_count=tile_count)
  if (fms_id > 0) then
    primary_id = get_new_diag_id(diag_cs)
    diag_cs%diags(primary_id)%fms_diag_id = fms_id
  endif

  register_static_field = primary_id

end function register_static_field

subroutine describe_option(opt_name, value, diag_CS)
  character(len=*),    intent(in) :: opt_name, value
  type(SIS_diag_ctrl), intent(in) :: diag_CS

  character(len=240) :: mesg
  integer :: start_ind = 1, end_ind, len_ind

  len_ind = len_trim(value)

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(diag_CS%doc_unit, '(a)') trim(mesg)
end subroutine describe_option

function i2s(a,n_in)
!   "Convert the first n elements of an integer array to a string."
    integer, dimension(:), intent(in) :: a
    integer, optional    , intent(in) :: n_in
    character(len=15) :: i2s

    character(len=15) :: i2s_temp
    integer :: i,n

    n=size(a)
    if(present(n_in)) n = n_in

    i2s = ''
    do i=1,n
       write (i2s_temp, '(I4.4)') a(i)
       i2s = trim(i2s) //'_'// trim(i2s_temp)
    enddo
    i2s = adjustl(i2s)
end function i2s

subroutine SIS_diag_mediator_init(G, IG, param_file, diag_cs, component, err_msg, &
                                  doc_file_dir)
  type(SIS_hor_grid_type),    intent(inout) :: G
  type(ice_grid_type),        intent(inout) :: IG
  type(param_file_type),      intent(in)    :: param_file
  type(SIS_diag_ctrl),        intent(inout) :: diag_cs
  character(len=*), optional, intent(in)    :: component
  character(len=*), optional, intent(out)   :: err_msg
  character(len=*), optional, intent(in)    :: doc_file_dir

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.
  integer :: ios, new_unit
  logical :: opened, new_file
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt, doc_path
  character(len=40)  :: doc_file_param
  character(len=40)  :: mod  = "SIS_diag_mediator" ! This module's name.

  call diag_manager_init(err_msg=err_msg)

  ! Allocate list of all diagnostics
  allocate(diag_cs%diags(DIAG_ALLOC_CHUNK_SIZE))
  diag_cs%next_free_diag_id = 1
  diag_cs%diags(:)%in_use = .false.

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied ; diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

  if (is_root_pe() .and. (diag_CS%doc_unit < 0)) then
    if (present(component)) then
      doc_file_dflt = trim(component)//".available_diags"
      doc_file_param = trim(uppercase(component))//"_AVAILABLE_DIAGS_FILE"
    else
      write(this_pe,'(i6.6)') PE_here()
      doc_file_dflt = "available_diags."//this_pe
      doc_file_param = "AVAILABLE_DIAGS_FILE"
    endif
    call get_param(param_file, mod, trim(doc_file_param), doc_file, &
                 "A file into which to write a list of all available \n"//&
                 "sea ice diagnostics that can be included in a diag_table.", &
                 default=doc_file_dflt)
    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%doc_unit /= -1) new_file = .false.
    ! Find an unused unit number.
      do new_unit=512,42,-1
        inquire( new_unit, opened=opened)
        if (.not.opened) exit
      enddo

      if (opened) call SIS_error(FATAL, &
          "diag_mediator_init failed to find an unused unit number.")

      doc_path = doc_file
      if (present(doc_file_dir)) then ; if (len_trim(doc_file_dir) > 0) then
        doc_path = trim(slasher(doc_file_dir))//trim(doc_file)
      endif ; endif

      diag_CS%doc_unit = new_unit

      if (new_file) then
        open(diag_CS%doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%doc_unit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call SIS_error(FATAL, "Failed to open available diags file "//trim(doc_path)//".")
      endif
    endif
  endif

  call diag_masks_set(G, IG, -1.0e34, diag_cs)

end subroutine SIS_diag_mediator_init

subroutine diag_masks_set(G, IG, missing_value, diag_cs)
! Setup the 2d masks for diagnostics
  type(SIS_hor_grid_type), target, intent(in)    :: G
  type(ice_grid_type),             intent(inout) :: IG
  real,                            intent(in)    :: missing_value
  type(SIS_diag_ctrl),             intent(inout) :: diag_cs
  ! Local variables
  integer :: i, j, k, NkIce, CatIce

  NkIce = IG%NkIce ; CatIce = IG%CatIce

  diag_cs%mask2dT  => G%mask2dT
  diag_cs%mask2dBu => G%mask2dBu
  diag_cs%mask2dCu => G%mask2dCu
  diag_cs%mask2dCv => G%mask2dCv

  allocate(diag_cs%mask2dT_comp(G%isc:G%iec,G%jsc:G%jec))
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    diag_cs%mask2dT_comp(i,j) = diag_cs%mask2dT(i,j)
  enddo ; enddo

  allocate(diag_cs%mask3dTL(G%isd:G%ied,G%jsd:G%jed,1:NkIce))
  allocate(diag_cs%mask3dBuL(G%IsdB:G%IedB,G%JsdB:G%JedB,1:NkIce))
  allocate(diag_cs%mask3dCuL(G%IsdB:G%IedB,G%jsd:G%jed,1:NkIce))
  allocate(diag_cs%mask3dCvL(G%isd:G%ied,G%JsdB:G%JedB,1:NkIce))
  do k=1,NkIce
    diag_cs%mask3dTL(:,:,k)  = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBuL(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCuL(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvL(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo

  allocate(diag_cs%mask3dTi(G%isd:G%ied,G%jsd:G%jed,1:NkIce+1))
  allocate(diag_cs%mask3dBui(G%IsdB:G%IedB,G%JsdB:G%JedB,1:NkIce+1))
  allocate(diag_cs%mask3dCui(G%IsdB:G%IedB,G%jsd:G%jed,1:NkIce+1))
  allocate(diag_cs%mask3dCvi(G%isd:G%ied,G%JsdB:G%JedB,1:NkIce+1))
  do k=1,NkIce+1
    diag_cs%mask3dTi(:,:,k)  = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBui(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCui(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvi(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo

  allocate(diag_cs%mask3dTC(G%isd:G%ied,G%jsd:G%jed,0:CatIce))
  allocate(diag_cs%mask3dBuC(G%IsdB:G%IedB,G%JsdB:G%JedB,0:CatIce))
  allocate(diag_cs%mask3dCuC(G%IsdB:G%IedB,G%jsd:G%jed,0:CatIce))
  allocate(diag_cs%mask3dCvC(G%isd:G%ied,G%JsdB:G%JedB,0:CatIce))
  do k=0,CatIce
    diag_cs%mask3dTC(:,:,k)  = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBuC(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCuC(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvC(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo

  diag_cs%missing_value = missing_value

end subroutine diag_masks_set

subroutine SIS_diag_mediator_close_registration(diag_CS)
  type(SIS_diag_ctrl), intent(inout) :: diag_CS

  if (diag_CS%doc_unit > -1) then
    close(diag_CS%doc_unit) ; diag_CS%doc_unit = -2
  endif

end subroutine SIS_diag_mediator_close_registration

subroutine SIS_diag_mediator_end(time, diag_CS)
  type(time_type), intent(in) :: time
  type(SIS_diag_ctrl), intent(inout) :: diag_CS

  if (diag_CS%doc_unit > -1) then
    close(diag_CS%doc_unit) ; diag_CS%doc_unit = -3
  endif

end subroutine SIS_diag_mediator_end

! Allocate a new diagnostic id, it may be necessary to expand the diagnostics
! array.
function get_new_diag_id(diag_cs)

  integer :: get_new_diag_id
  type(SIS_diag_ctrl), intent(inout) :: diag_cs
  ! Arguments:
  !  (inout)   diag_cs  - diagnostics control structure

  type(diag_type), dimension(:), allocatable :: tmp
  integer :: i

  if (diag_cs%next_free_diag_id > size(diag_cs%diags)) then
    call assert(diag_cs%next_free_diag_id - size(diag_cs%diags) == 1, &
                'get_new_diag_id: inconsistent diag id')

    ! Increase the size of diag_cs%diags and copy data over.
    ! Do not use move_alloc() because it is not supported by Fortran 90
    allocate(tmp(size(diag_cs%diags)))
    tmp(:) = diag_cs%diags(:)
    deallocate(diag_cs%diags)
    allocate(diag_cs%diags(size(tmp) + DIAG_ALLOC_CHUNK_SIZE))
    diag_cs%diags(1:size(tmp)) = tmp(:)
    deallocate(tmp)

    ! Initialise new part of the diag array.
    do i=diag_cs%next_free_diag_id, size(diag_cs%diags)
      diag_cs%diags(i)%in_use = .false.
    enddo
  endif

  get_new_diag_id = diag_cs%next_free_diag_id
  diag_cs%next_free_diag_id = diag_cs%next_free_diag_id + 1

end function get_new_diag_id

subroutine assert(logical_arg, msg)

  logical, intent(in) :: logical_arg
  character(len=*), intent(in) :: msg

  if (.not. logical_arg) then
    call SIS_error(FATAL, 'Assert failed: '//msg)
  endif

end subroutine assert

end module SIS_diag_mediator
