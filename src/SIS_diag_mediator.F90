!> Convenient wrappers to the FMS diag_manager interfaces with additional diagnostic capabilities.
module SIS_diag_mediator

! This file is a part of SIS2. See LICENSE.md for the license.

use ice_grid,               only : ice_grid_type
use MOM_checksums,          only : chksum0, hchksum, uchksum, vchksum, Bchksum
use MOM_coms,               only : PE_here
use MOM_cpu_clock,          only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,          only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_manager_infra, only : diag_manager_init=>MOM_diag_manager_init, MOM_diag_manager_end
use MOM_diag_manager_infra, only : diag_axis_init=>MOM_diag_axis_init, get_MOM_diag_axis_name
use MOM_diag_manager_infra, only : send_data_infra, MOM_diag_field_add_attribute, EAST, NORTH
use MOM_diag_manager_infra, only : register_diag_field_infra, register_static_field_infra
use MOM_diag_manager_infra, only : get_MOM_diag_field_id, DIAG_FIELD_NOT_FOUND
use MOM_diag_manager_infra, only : diag_send_complete_infra
use MOM_error_handler,      only : SIS_error=>MOM_error, FATAL, is_root_pe, assert, callTree_showQuery
use MOM_error_handler,      only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,        only : get_param, log_param, log_version, param_file_type
use MOM_io,                 only : get_filename_appendix
use MOM_safe_alloc,         only : safe_alloc_ptr, safe_alloc_alloc
use MOM_string_functions,   only : lowercase, uppercase, slasher
use MOM_time_manager,       only : time_type, get_time
use MOM_unit_scaling,       only : unit_scale_type
use SIS_hor_grid,           only : SIS_hor_grid_type

implicit none ; private

public set_SIS_axes_info, post_SIS_data, register_SIS_diag_field, time_type
public safe_alloc_ptr, safe_alloc_alloc
public enable_SIS_averaging, enable_SIS_averages, disable_SIS_averaging, query_SIS_averaging_enabled
public SIS_diag_mediator_init, SIS_diag_mediator_end, set_SIS_diag_mediator_grid
public SIS_diag_mediator_close_registration, get_SIS_diag_time_end
public diag_axis_init, register_static_field, SIS_diag_send_complete
public register_scalar_field
public define_axes_group, diag_masks_set
public diag_register_area_ids
public found_in_diagtable

!> Make a diagnostic available for averaging or output.
interface post_SIS_data
  module procedure post_data_3d, post_data_2d, post_data_0d
end interface post_SIS_data

!> A group of 1D axes that comprise a 1D/2D/3D mesh
type, public :: axes_grp
  character(len=15) :: id   !< The id string for this particular combination of handles.
  integer           :: rank !< Number of dimensions in the list of axes.
  integer, dimension(:), allocatable :: handles !< Handles to 1D axes.
  type(SIS_diag_ctrl), pointer :: diag_cs => null() !< Circular link back to the main diagnostics control structure
                                                !! (Used to avoid passing said structure into every possible call).
  ! ID's for cell_methods
  character(len=9) :: x_cell_method = '' !< Default nature of data representation, if axes group
                                         !! includes x-direction.
  character(len=9) :: y_cell_method = '' !< Default nature of data representation, if axes group
                                         !! includes y-direction.
  ! For detecting position on the grid
  logical :: is_h_point = .false. !< If true, indicates that this axes group is for an h-point located field.
  logical :: is_q_point = .false. !< If true, indicates that this axes group is for a q-point located field.
  logical :: is_u_point = .false. !< If true, indicates that this axes group is for a u-point located field.
  logical :: is_v_point = .false. !< If true, indicates that this axes group is for a v-point located field.
  logical :: is_layer = .false. !< If true, indicates that this axes group is for a layer vertically-located field.
  logical :: is_interface = .false. !< If true, indicates that this axes group is for an interface
                                    !! vertically-located field.
  logical :: is_category = .false. !< If true, indicates that this axes group is for the ice thickness
                                   !! categories, exclusive of the ice-free category.
  logical :: is_cat_open = .false. !< If true, indicates that this axes group is for the ice thickness
                                   !! categories, inclusive of the ice-free open-water category.

  ! ID's for cell_measures
  integer :: id_area = -1 !< The diag_manager id for area to be used for cell_measure of variables with this axes_grp.
  ! For masking
  real, pointer, dimension(:,:)   :: mask2d => null() !< Mask for 2d (x-y) axes [nondim]
  real, pointer, dimension(:,:)   :: mask2d_comp => null() !< Mask for 2-d axes on the computational
                                                      !! domain for this diagnostic [nondim]
  real, pointer, dimension(:,:,:) :: mask3d => null() !< Mask for 3d axes [nondim]
end type axes_grp

!> This type is used to represent a diagnostic at the diag_mediator level.
!!
!! There can be both 'primary' and 'secondary' diagnostics. The primaries
!! reside in the diag_cs%diags array. They have an id which is an index
!! into this array. The secondaries are 'variations' on the primary diagnostic.
!! For example the CMOR diagnostics are secondary. The secondary diagnostics
!! are kept in a list with the primary diagnostic as the head.
type, private :: diag_type
  logical :: in_use              !< True if this entry is being used.
  integer :: fms_diag_id         !< Underlying FMS diag_manager id.
  character(len=64) :: debug_str = '' !< The diagnostic name and module or FATAL errors and debugging.
  type(axes_grp), pointer :: axes => null() !< The axis group for this diagnostic
  type(diag_type), pointer :: next => null() !< Pointer to the next diagnostic
  real :: conversion_factor = 0. !< If non-zero, a factor to multiply data by before posting to FMS,
                                 !! often including factors to undo internal scaling in units of [a A-1 ~> 1]
end type diag_type

!>   The SIS_diag_ctrl data type contains times to regulate diagnostics along with masks and
!! axes to use with diagnostics, and a list of structures with data about each diagnostic.
type, public :: SIS_diag_ctrl
  integer :: available_diag_doc_unit = -1 !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  integer :: chksum_iounit = -1           !< The unit number of a diagnostic documentation file.
                                          !! This file is open if available_diag_doc_unit is > 0.
  logical :: diag_as_chksum  !< If true, log chksums in a text file instead of posting diagnostics
  logical :: show_call_tree  !< Display the call tree while running. Set by VERBOSITY level.
  logical :: index_space_axes !< If true, diagnostic horizontal coordinates axes are in index space.

  ! The following fields are used for the output of the data.
  ! These give the computational-domain sizes, and are relative to a start value
  ! of 1 in memory for the tracer-point arrays.
  integer :: is  !< The start i-index of cell centers within the computational domain
  integer :: ie  !< The end i-index of cell centers within the computational domain
  integer :: js  !< The start j-index of cell centers within the computational domain
  integer :: je  !< The end j-index of cell centers within the computational domain
  ! These give the memory-domain sizes, and can be start at any value on each PE.
  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain
  real :: time_int              !< The time interval for any fields
                                !! that are offered for averaging [s].
  type(time_type) :: time_end   !< The end time of the valid interval for any offered field.
  logical :: ave_enabled = .false. !< True if averaging is enabled.

  !>@{ The following are 3D and 2D axis groups defined for output.  The names indicate
  !! the horizontal locations (B, T, Cu, or Cv) and vertical locations (L, i, or 1) or
  !! ice thickness categories (C or C0).
  type(axes_grp) :: axesBL, axesTL, axesCuL, axesCvL
  type(axes_grp) :: axesBi, axesTi, axesCui, axesCvi
  type(axes_grp) :: axesB1, axesT1, axesCu1, axesCv1
  type(axes_grp) :: axesBC0, axesTC0, axesCuC0, axesCvC0
  type(axes_grp) :: axesBC, axesTC, axesCuC, axesCvC
  !>@}
  type(axes_grp) :: axesZi !< A 1-D z-space axis at interfaces
  type(axes_grp) :: axesZL !< A 1-D z-space axis at layer centers
  type(axes_grp) :: axesNull !< An axis group for scalars

  ! Mask arrays for 2D diagnostics
  real, dimension(:,:),   pointer :: mask2dT   => null() !< 2D mask array for cell-center points [nondim]
  real, dimension(:,:),   pointer :: mask2dBu  => null() !< 2D mask array for cell-corner points [nondim]
  real, dimension(:,:),   pointer :: mask2dCu  => null() !< 2D mask array for east-face points [nondim]
  real, dimension(:,:),   pointer :: mask2dCv  => null() !< 2D mask array for north-face points [nondim]
  !> Computational domain mask arrays for 2D diagnostics [nondim]
  real, dimension(:,:),   pointer :: mask2dT_comp => null()

  ! 3D mask arrays for diagnostics at layers (mask...L) and interfaces (mask...i) all [nondim]
  real, dimension(:,:,:), pointer :: mask3dTL  => null() !< 3D mask array for layer cell-center points [nondim]
  real, dimension(:,:,:), pointer :: mask3dBL  => null() !< 3D mask array for layer cell-corner points [nondim]
  real, dimension(:,:,:), pointer :: mask3dCuL => null() !< 3D mask array for layer east-face points [nondim]
  real, dimension(:,:,:), pointer :: mask3dCvL => null() !< 3D mask array for layer north-faces [nondim]
  real, dimension(:,:,:), pointer :: mask3dTi  => null() !< 3D mask array for interface cell-centers [nondim]
  real, dimension(:,:,:), pointer :: mask3dBi  => null() !< 3D mask array for interface cell-corners [nondim]
  real, dimension(:,:,:), pointer :: mask3dCui => null() !< 3D mask array for interface east-faces [nondim]
  real, dimension(:,:,:), pointer :: mask3dCvi => null() !< 3D mask array for interface north-faces [nondim]

  ! 3D mask arrays for diagnostics by ice thickness category
  real, dimension(:,:,:), pointer :: mask3dTC0  => null() !< 3D mask array at cell-centers by ice thickness category
                                                          !! inclusive of the open water category [nondim]
  real, dimension(:,:,:), pointer :: mask3dBuC0 => null() !< 3D mask array at cell-corners by ice thickness category
                                                          !! inclusive of the open water category [nondim]
  real, dimension(:,:,:), pointer :: mask3dCuC0 => null() !< 3D mask array at east-faces by ice thickness category
                                                          !! inclusive of the open water category [nondim]
  real, dimension(:,:,:), pointer :: mask3dCvC0 => null() !< 3D mask array at north-faces by ice thickness category
                                                          !! inclusive of the open water category [nondim]
  real, dimension(:,:,:), pointer :: mask3dTC  => null() !< 3D mask array at cell-centers by ice thickness category
                                                         !! exclusive of the open water category [nondim]
  real, dimension(:,:,:), pointer :: mask3dBuC => null() !< 3D mask array at cell-corners by ice thickness category
                                                         !! exclusive of the open water category [nondim]
  real, dimension(:,:,:), pointer :: mask3dCuC => null() !< 3D mask array at east-faces by ice thickness category
                                                         !! exclusive of the open water category [nondim]
  real, dimension(:,:,:), pointer :: mask3dCvC => null() !< 3D mask array at north-faces by ice thickness category
                                                         !! exclusive of the open water category [nondim]

! Space for diagnostics is dynamically allocated as it is needed.
! The chunk size is how much the array should grow on each new allocation.
#define DIAG_ALLOC_CHUNK_SIZE 15
  type(diag_type), dimension(:), allocatable :: diags !< The list of diagnostics
  integer :: next_free_diag_id !< The next unused diagnostic ID

  !> default missing value to be sent to ALL diagnostics registrations [various]
  real :: missing_value = -1.0e34

  type(SIS_hor_grid_type), pointer :: G => null()  !< The ocean grid type
  type(unit_scale_type),   pointer :: US => null() !< A dimensional unit scaling type

  !> Number of checksum-only diagnostics
  integer :: num_chksum_diags

end type SIS_diag_ctrl

!>@{ CPU clocks
integer :: id_clock_diag_mediator
!>@}

contains

!> Set up the grid and axis information for use by SIS.
subroutine set_SIS_axes_info(G, IG, param_file, diag_cs, set_vertical, axes_set_name)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file structure
  type(SIS_diag_ctrl),     intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
  logical,          optional, intent(in) :: set_vertical !< If true (or missing), set up the vertical axes
  character(len=*), optional, intent(in) :: axes_set_name !<  A name to use for this set of axes.
                                                !! The default is "ice".

  ! Local variables
  integer :: id_xq, id_yq, id_zl, id_zi, id_xh, id_yh, id_null
  integer :: id_xhe, id_yhe
  integer :: id_ct, id_ct0 !, id_ct_bnd, id_ct0_bnd
  integer :: i, j, k, ncat
  logical :: set_vert
  real :: zlev_ice(IG%NkIce)   ! Fractional position labels for the vertical layers in the ice [nondim]
  real :: zinter_ice(IG%NkIce+1) ! Fractional position labels for the vertical interfaces in the ice [nondim]
  character(len=80) :: grid_config, units_temp, set_name
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "SIS_diag_mediator" ! This module's name.
  real, allocatable, dimension(:) :: IaxB, iax ! Index-based integer and half-integer i-axis labels [nondim]
  real, allocatable, dimension(:) :: JaxB, jax ! Index-based integer and half-integer j-axis labels [nondim]
  ! real, allocatable, dimension(:) :: cat_vals  ! Nominal ice thicknesses of the thickness categories [m]
  ! real, allocatable, dimension(:) :: cat_edge_vals ! Ice thicknesses bounds of the categories [m]

  set_vert = .true. ; if (present(set_vertical)) set_vert = set_vertical
  set_name = "ice" ; if (present(axes_set_name)) set_name = trim(axes_set_name)

  ! This is inconsistent with the labeling of axis units from MOM6, and it will be corrected in a subsequent commit.

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "GRID_CONFIG", grid_config, &
                 "The method for defining the horizontal grid.  Valid "//&
                 "entries include:\n"//&
                 "\t file - read the grid from GRID_FILE \n"//&
                 "\t mosaic - read the grid from a mosaic grid file \n"//&
                 "\t cartesian - a Cartesian grid \n"//&
                 "\t spherical - a spherical grid \n"//&
                 "\t mercator  - a Mercator grid", fail_if_missing=.true.)

  G%x_axis_units = "degrees_E"
  G%y_axis_units = "degrees_N"
  if (index(lowercase(trim(grid_config)),"cartesian") > 0) then
    ! This is a Cartesian grid, and may have different axis units.
    call get_param(param_file, mdl, "AXIS_UNITS", units_temp, &
                 "The units for the x- and y- axis labels.  AXIS_UNITS "//&
                 "should be defined as 'k' for km, 'm' for m, or 'd' "//&
                 "for degrees of latitude and longitude (the default). "//&
                 "Except on a Cartesian grid, only degrees are currently "//&
                 "implemented.", default='degrees')
    if (units_temp(1:1) == 'k') then
      G%x_axis_units = "kilometers" ; G%y_axis_units = "kilometers"
    elseif (units_temp(1:1) == 'm') then
      G%x_axis_units = "meters" ; G%y_axis_units = "meters"
    endif
    call log_param(param_file, mdl, "explicit AXIS_UNITS", G%x_axis_units)
  endif


  if (diag_cs%index_space_axes) then
    allocate(IaxB(G%IsgB:G%IegB))
    do i=G%IsgB, G%IegB
      Iaxb(i)=real(i)
    enddo
    allocate(iax(G%isg:G%ieg))
    do i=G%isg, G%ieg
      iax(i)=real(i)-0.5
    enddo
    allocate(JaxB(G%JsgB:G%JegB))
    do j=G%JsgB, G%JegB
      JaxB(j)=real(j)
    enddo
    allocate(jax(G%jsg:G%jeg))
    do j=G%jsg, G%jeg
      jax(j)=real(j)-0.5
    enddo
  endif

  ! Horizontal axes for the native grids, noting that SIS2 always uses symmetric memory.
  if (diag_cs%index_space_axes) then
    id_xq = diag_axis_init('xB', IaxB(G%isgB:G%iegB), 'none', 'x', &
        'Boundary point grid-space longitude', G%Domain, set_name=set_name, position=EAST)
    id_yq = diag_axis_init('yB', JaxB(G%jsgB:G%jegB), 'none', 'y', &
        'Boundary point grid-space latitude', G%Domain, set_name=set_name, position=NORTH)
    id_xhe = diag_axis_init('xTe', IaxB(G%isg-1:G%ieg), 'none', 'x', &
        'T-cell edge grid-space longitude', G%Domain, set_name=set_name, position=EAST)
    id_yhe = diag_axis_init('yTe', JaxB(G%jsg-1:G%jeg), 'none', 'y', &
        'T-cell edge grid-space latitude', G%Domain, set_name=set_name, position=NORTH)
    id_xh = diag_axis_init('xT', iax(G%isg:G%ieg), 'none', 'x', &
        'T point grid-space longitude', G%Domain, set_name=set_name, edges=id_xhe)
    id_yh = diag_axis_init('yT', jax(G%jsg:G%jeg), 'none', 'y', &
        'T point grid-space latitude', G%Domain, set_name=set_name, edges=id_yhe)
  else
    id_xq = diag_axis_init('xB', G%gridLonB(G%isgB:G%iegB), G%x_axis_units, 'x', &
        'Boundary point nominal longitude', G%Domain, set_name=set_name, position=EAST)
    id_yq = diag_axis_init('yB', G%gridLatB(G%jsgB:G%jegB), G%y_axis_units, 'y', &
        'Boundary point nominal latitude', G%Domain, set_name=set_name, position=NORTH)
    id_xhe = diag_axis_init('xTe', G%gridLonB(G%isg-1:G%ieg), G%x_axis_units, 'x', &
        'T-cell edge nominal longitude', G%Domain, set_name=set_name, position=EAST)
    id_yhe = diag_axis_init('yTe', G%gridLatB(G%jsg-1:G%jeg), G%y_axis_units, 'y', &
        'T-cell edge nominal latitude', G%Domain, set_name=set_name, position=NORTH)
    id_xh = diag_axis_init('xT', G%gridLonT(G%isg:G%ieg), G%x_axis_units, 'x', &
        'T point nominal longitude', G%Domain, set_name=set_name, edges=id_xhe)
    id_yh = diag_axis_init('yT', G%gridLatT(G%jsg:G%jeg), G%y_axis_units, 'y', &
        'T point nominal latitude', G%Domain, set_name=set_name, edges=id_yhe)
  endif

  if (set_vert) then
    do k=1,IG%NkIce+1 ; zinter_ice(k) = real(k-1) / real(IG%NkIce) ; enddo
    do k=1,IG%NkIce ; zlev_ice(k) = (k-0.5) / real(IG%NkIce) ; enddo
    id_zl = diag_axis_init('zl', zlev_ice, 'layer', 'z', 'Cell depth', set_name=set_name)
    id_zi = diag_axis_init('zi', zinter_ice, 'interface', 'z', &
                           'Cell interface depth', set_name=set_name)
  else
    id_zl = -1 ; id_zi = -1
  endif

  ! Ice thickness category axes

  ! The commented out version would set central values and bounds as separate axes.
  ! ncat = IG%CatIce
  ! allocate(cat_edge_vals(0:ncat+1), source=0.0)
  ! allocate(cat_vals(0:ncat), source=0.0)
  ! cat_edge_vals(0) = 0.0 ; cat_edge_vals(1:ncat+1) = IG%cat_thick_lim(1:ncat+1)
  ! cat_vals(0) = 0.0
  ! do k=1,ncat ; cat_vals(k) = 0.5*(cat_edge_vals(K+1) + cat_edge_vals(K)) ; enddo
  ! id_ct_bnd = diag_axis_init('cat_bnd', cat_vals(1:ncat+1), units="m", cart_name="n", &
  !                        long_name="Ice thickness category bounds")
  ! id_ct0_bnd = diag_axis_init('cat0_bnd', cat_edge_vals(0:ncat+1), units="m", cart_name="n", &
  !                         long_name="Ice thickness category bounds including open water")
  ! id_ct = diag_axis_init('cat', cat_vals(1:ncat), units="m", cart_name="n", &
  !                        long_name="Ice thickness category", edges=id_ct_bnd)
  ! id_ct0 = diag_axis_init('cat0', cat_vals(0:ncat), units="m", cart_name="n", &
  !                         long_name="Ice thickness category including open water", edges=id_ct0_bnd)
  ! deallocate(cat_edge_vals, cat_vals)

  id_ct = diag_axis_init('ct', IG%cat_thick_lim(1:IG%CatIce), 'meters', 'n', &
                         'Ice thickness category bounds', set_name=set_name)
  id_ct0 = diag_axis_init('ctu', IG%cat_thick_lim(1:IG%CatIce+1), 'meters', 'n', &
                         'Ice thickness category upper bounds', set_name=set_name)

  ! ToDo: Consider adding 4-d axis groupings.

  ! Vertical axes for the interfaces and layers
  call define_axes_group(diag_cs, (/ id_zi /), diag_cs%axesZi, is_interface=.true.)
  call define_axes_group(diag_cs, (/ id_zL /), diag_cs%axesZL, is_layer=.true.)

  ! Axis groupings for the model layers
  call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zL /), diag_cs%axesTL, &
       x_cell_method='mean', y_cell_method='mean', &
       is_h_point=.true., is_layer=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zL /), diag_cs%axesBL, &
       x_cell_method='point', y_cell_method='point', &
       is_q_point=.true., is_layer=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zL /), diag_cs%axesCuL, &
       x_cell_method='point', y_cell_method='mean', &
       is_u_point=.true., is_layer=.true.)
  call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zL /), diag_cs%axesCvL, &
       x_cell_method='mean', y_cell_method='point', &
       is_v_point=.true., is_layer=.true.)

  ! Axis groupings for the model interfaces
  call define_axes_group(diag_cs, (/ id_xh, id_yh, id_zi /), diag_cs%axesTi, &
       x_cell_method='mean', y_cell_method='mean', &
       is_h_point=.true., is_interface=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yq, id_zi /), diag_cs%axesBi, &
       x_cell_method='point', y_cell_method='point', &
       is_q_point=.true., is_interface=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh, id_zi /), diag_cs%axesCui, &
       x_cell_method='point', y_cell_method='mean', &
       is_u_point=.true., is_interface=.true.)
  call define_axes_group(diag_cs, (/ id_xh, id_yq, id_zi /), diag_cs%axesCvi, &
       x_cell_method='mean', y_cell_method='point', &
       is_v_point=.true., is_interface=.true.)

  ! Axis groupings for the ice thickness categories, exclusive of the open water category
  call define_axes_group(diag_cs, (/ id_xh, id_yh, id_ct /), diag_cs%axesTC, &
       x_cell_method='mean', y_cell_method='mean', &
       is_h_point=.true., is_category=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yq, id_ct /), diag_cs%axesBC, &
       x_cell_method='point', y_cell_method='point', &
       is_q_point=.true., is_category=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh, id_ct /), diag_cs%axesCuC, &
       x_cell_method='point', y_cell_method='mean', &
       is_u_point=.true., is_category=.true.)
  call define_axes_group(diag_cs, (/ id_xh, id_yq, id_ct /), diag_cs%axesCvC, &
       x_cell_method='mean', y_cell_method='point', &
       is_v_point=.true., is_category=.true.)

  ! Axis groupings for the ice thickness categories, inclusive of the open water category
  call define_axes_group(diag_cs, (/ id_xh, id_yh, id_ct0 /), diag_cs%axesTC0, &
       x_cell_method='mean', y_cell_method='mean', &
       is_h_point=.true., is_cat_open=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yq, id_ct0 /), diag_cs%axesBC0, &
       x_cell_method='point', y_cell_method='point', &
       is_q_point=.true., is_cat_open=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh, id_ct0 /), diag_cs%axesCuC0, &
       x_cell_method='point', y_cell_method='mean', &
       is_u_point=.true., is_cat_open=.true.)
  call define_axes_group(diag_cs, (/ id_xh, id_yq, id_ct0 /), diag_cs%axesCvC0, &
       x_cell_method='mean', y_cell_method='point', &
       is_v_point=.true., is_cat_open=.true.)

  ! Axis groupings for 2-D arrays
  call define_axes_group(diag_cs, (/ id_xh, id_yh /), diag_cs%axesT1, &
       x_cell_method='mean', y_cell_method='mean', is_h_point=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yq /), diag_cs%axesB1, &
       x_cell_method='point', y_cell_method='point', is_q_point=.true.)
  call define_axes_group(diag_cs, (/ id_xq, id_yh /), diag_cs%axesCu1, &
       x_cell_method='point', y_cell_method='mean', is_u_point=.true.)
  call define_axes_group(diag_cs, (/ id_xh, id_yq /), diag_cs%axesCv1, &
       x_cell_method='mean', y_cell_method='point', is_v_point=.true.)

  ! Axis group for special null axis for scalars from diag manager.
  id_null = diag_axis_init('scalar_axis', (/0./), 'none', 'N', 'none', null_axis=.true.)
  call define_axes_group(diag_cs, (/ id_null /), diag_cs%axesNull)

  if (diag_cs%index_space_axes) then
    deallocate(IaxB, iax, JaxB, jax)
  endif

end subroutine set_SIS_axes_info

!> Attaches the id of cell areas to axes groups for use with cell_measures
subroutine diag_register_area_ids(diag_cs, id_area_t, id_area_q)
  type(SIS_diag_ctrl), intent(inout) :: diag_cs   !< Diagnostics control structure
  integer,   optional, intent(in)    :: id_area_t !< Diag_mediator id for area of h-cells
  integer,   optional, intent(in)    :: id_area_q !< Diag_mediator id for area of q-cells
  ! Local variables
  integer :: fms_id, i
  if (present(id_area_t)) then
    fms_id = diag_cs%diags(id_area_t)%fms_diag_id
    diag_cs%axesT1%id_area = fms_id
    diag_cs%axesTi%id_area = fms_id
    diag_cs%axesTL%id_area = fms_id
    diag_cs%axesTC%id_area = fms_id
    diag_cs%axesTC0%id_area = fms_id
  endif
  if (present(id_area_q)) then
    fms_id = diag_cs%diags(id_area_q)%fms_diag_id
    diag_cs%axesB1%id_area = fms_id
    diag_cs%axesBi%id_area = fms_id
    diag_cs%axesBL%id_area = fms_id
    diag_cs%axesBC%id_area = fms_id
    diag_cs%axesBC0%id_area = fms_id
  endif
end subroutine diag_register_area_ids

!> Defines a group of "axes" from list of handles
subroutine define_axes_group(diag_cs, handles, axes, &
                             x_cell_method, y_cell_method, &
                             is_h_point, is_q_point, is_u_point, is_v_point, &
                             is_layer, is_interface, is_category, is_cat_open)
  type(SIS_diag_ctrl), target,    intent(in)  :: diag_cs !< Diagnostics control structure
  integer, dimension(:),      intent(in)  :: handles !< A list of 1D axis handles
  type(axes_grp),             intent(out) :: axes    !< The group of 1D axes
  character(len=*), optional, intent(in)  :: x_cell_method !< A x-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  character(len=*), optional, intent(in)  :: y_cell_method !< A y-direction cell method used to construct the
                                                           !! "cell_methods" attribute in CF convention
  logical,          optional, intent(in)  :: is_h_point !< If true, indicates this axes group for h-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_q_point !< If true, indicates this axes group for q-point
                                                        !! located fields
  logical,          optional, intent(in)  :: is_u_point !< If true, indicates this axes group for
                                                        !! u-point located fields
  logical,          optional, intent(in)  :: is_v_point !< If true, indicates this axes group for
                                                        !! v-point located fields
  logical,          optional, intent(in)  :: is_layer   !< If true, indicates that this axes group is
                                                        !! for a layer vertically-located field.
  logical,          optional, intent(in)  :: is_interface !< If true, indicates that this axes group
                                                        !! is for an interface vertically-located field.
  logical,          optional, intent(in)  :: is_category !< If true, indicates that this axes group is
                                                        !! for the ice thickness categories, exclusive
                                                        !! of the ice-free category.
  logical,          optional, intent(in)  :: is_cat_open !< If true, indicates that this axes group is
                                                        !! for the ice thickness categories, inclusive
                                                        !! of the ice-free open-water category.
  ! Local variables
  integer :: n

  n = size(handles)
  if (n<1 .or. n>3) call SIS_error(FATAL, "define_axes_group: wrong size for list of handles!")
  allocate( axes%handles(n) )
  axes%id = i2s(handles, n) ! Identifying string
  axes%rank = n
  axes%handles(:) = handles(:)
  axes%diag_cs => diag_cs ! A (circular) link back to the SIS_diag_ctrl structure

  if ((axes%rank<2) .and. (present(x_cell_method) .or. present(x_cell_method))) &
    call SIS_error(FATAL, 'define_axes_group: Can not set x_cell_method or y_cell_method for rank<2.')
  axes%x_cell_method = '' ; if (present(x_cell_method)) axes%x_cell_method = trim(x_cell_method)
  axes%y_cell_method = '' ; if (present(y_cell_method)) axes%y_cell_method = trim(y_cell_method)

  if (present(is_h_point)) axes%is_h_point = is_h_point
  if (present(is_q_point)) axes%is_q_point = is_q_point
  if (present(is_u_point)) axes%is_u_point = is_u_point
  if (present(is_v_point)) axes%is_v_point = is_v_point
  if (present(is_layer)) axes%is_layer = is_layer
  if (present(is_interface)) axes%is_interface = is_interface
  if (present(is_category)) axes%is_category = is_category
  if (present(is_cat_open)) axes%is_cat_open = is_cat_open

  ! Setup masks for this axes group
  axes%mask2d => null()
  if (axes%rank==2) then
    if (axes%is_h_point) axes%mask2d => diag_cs%mask2dT
    if (axes%is_h_point) axes%mask2d_comp => diag_cs%mask2dT_comp
    if (axes%is_u_point) axes%mask2d => diag_cs%mask2dCu
    if (axes%is_v_point) axes%mask2d => diag_cs%mask2dCv
    if (axes%is_q_point) axes%mask2d => diag_cs%mask2dBu
  endif

  axes%mask3d => null()
  if (axes%rank==3) then
    ! Native variables can/should use the native masks copied into diag_cs
    if (axes%is_layer) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTL
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCuL
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvL
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBL
    elseif (axes%is_interface) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTi
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCui
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvi
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBi
    elseif (axes%is_category) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTC
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCuC
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvC
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBuC
    elseif (axes%is_cat_open) then
      if (axes%is_h_point) axes%mask3d => diag_cs%mask3dTC0
      if (axes%is_u_point) axes%mask3d => diag_cs%mask3dCuC0
      if (axes%is_v_point) axes%mask3d => diag_cs%mask3dCvC0
      if (axes%is_q_point) axes%mask3d => diag_cs%mask3dBuC0
    endif
  endif

end subroutine define_axes_group

!> Set up the array extents for doing diagnostics
subroutine set_SIS_diag_mediator_grid(G, diag_cs)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(SIS_diag_ctrl),     intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

end subroutine set_SIS_diag_mediator_grid

!> Make a real scalar diagnostic available for averaging or output
subroutine post_data_0d(diag_field_id, field, diag_cs, is_static)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_SIS_diag_field.
  real,              intent(in) :: field         !< real value being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(SIS_diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.

  ! Local variables
  real :: locfield ! The field being offered in arbitrary unscaled units [a]
  logical :: used, is_stat
  type(diag_type), pointer :: diag => null()

  integer :: time_days
  integer :: time_seconds
  character(len=300) :: debug_mesg

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Iterate over list of diag 'variants', e.g. CMOR aliases, call send_data
  ! for each one.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_0d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)

  do while (associated(diag))
    locfield = field
    if (diag%conversion_factor /= 0.) &
      locfield = locfield * diag%conversion_factor

    if (diag_cs%diag_as_chksum) then
      ! Append timestep to mesg
      call get_time(diag_cs%time_end, time_seconds, days=time_days)
      write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
          trim(diag%debug_str), time_days, time_seconds

      call chksum0(locfield, debug_mesg, logunit=diag_cs%chksum_iounit)
    elseif (is_stat) then
      used = send_data_infra(diag%fms_diag_id, locfield)
    elseif (diag_cs%ave_enabled) then
      used = send_data_infra(diag%fms_diag_id, locfield, diag_cs%time_end)
    endif
    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end subroutine post_data_0d


!> Make a real 2-d array diagnostic available for averaging or output
subroutine post_data_2d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_SIS_diag_field.
  real,      target, intent(in) :: field(:,:)    !< 2-d array being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(SIS_diag_ctrl), target, intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables
  type(diag_type), pointer :: diag => null()

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! Iterate over list of diag 'variants' (e.g. CMOR aliases) and post each.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_2d: Unregistered diagnostic id')
  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    call post_data_2d_low(diag, field, diag_cs, is_static, mask)
    diag => diag%next
  enddo

  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)
end subroutine post_data_2d

!> Make a real 2-d array diagnostic available for averaging or output
!! using a diag_type instead of an integer id.
subroutine post_data_2d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag       !< A structure describing the diagnostic to post
  real,    target,   intent(in) :: field(:,:) !< 2-d array being offered for output or averaging
                                              !! in internally scaled arbitrary units [A ~> a]
  type(SIS_diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real, optional, target, intent(in) :: mask(:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables
  real, dimension(:,:), pointer :: locfield ! The field being offered in arbitrary unscaled units [a]
  real, dimension(:,:), pointer :: locmask  ! A pointer to the data mask to use [nondim]
  character(len=300) :: mesg
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  logical :: i_data, j_data ! True if the field is on the data domain in the i or j directions.
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, i, j
  integer :: time_days, time_seconds
  character(len=300) :: debug_mesg

  locfield => NULL()
  locmask => NULL()
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  cszi = (diag_cs%ie-diag_cs%is) +1 ; dszi = (diag_cs%ied-diag_cs%isd) +1
  cszj = (diag_cs%je-diag_cs%js) +1 ; dszj = (diag_cs%jed-diag_cs%jsd) +1
  if ( size(field,1) == dszi ) then
    isv = diag_cs%is ; iev = diag_cs%ie ; i_data = .true.   ! Data domain
  elseif ( size(field,1) == dszi + 1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1 ; i_data = .true. ! Symmetric data domain
  elseif ( size(field,1) == cszi ) then
    isv = 1 ; iev = cszi ; i_data = .false. ! Computational domain
  elseif ( size(field,1) == cszi + 1 ) then
    isv = 1 ; iev = cszi+1 ; i_data = .false. ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,1)," in i-direction\n"//&
       "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call SIS_error(FATAL,"post_SIS_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ( size(field,2) == dszj ) then
    jsv = diag_cs%js ; jev = diag_cs%je ; j_data = .true.   ! Data domain
  elseif ( size(field,2) == dszj + 1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1 ; j_data = .true. ! Symmetric data domain
  elseif ( size(field,2) == cszj ) then
    jsv = 1 ; jev = cszj ; j_data = .false. ! Computational domain
  elseif ( size(field,2) == cszj + 1 ) then
    jsv = 1 ; jev = cszj+1 ; j_data = .false. ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,2)," in j-direction\n"//&
       "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call SIS_error(FATAL,"post_SIS_data_2d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2) ) )
    do j=jsv,jev ; do i=isv,iev
      if (field(i,j) == diag_cs%missing_value) then
        locfield(i,j) = diag_cs%missing_value
      else
        locfield(i,j) = field(i,j) * diag%conversion_factor
      endif
    enddo ; enddo
    locfield(isv:iev,jsv:jev) = field(isv:iev,jsv:jev) * diag%conversion_factor
  else
    locfield => field
  endif

  ! Handle cases where the data and computational domain are the same size.
  if (diag_cs%ied-diag_cs%isd == diag_cs%ie-diag_cs%is) i_data = j_data
  if (diag_cs%jed-diag_cs%jsd == diag_cs%je-diag_cs%js) j_data = i_data
  if ( i_data .NEQV. j_data ) then
    call SIS_error(FATAL, "post_SIS_data_2d: post_SIS_data called for "//&
                   trim(diag%debug_str)//" with mixed computational and data domain array sizes.")
  endif

  if (present(mask)) then
    locmask => mask
  elseif (.not.is_stat) then  ! Static fields do not have assigned axes.
    if (i_data .and. associated(diag%axes%mask2d)) then
      locmask => diag%axes%mask2d
    elseif ((.not.i_data) .and. associated(diag%axes%mask2d_comp)) then
      locmask => diag%axes%mask2d_comp
    endif
  endif
  if (associated(locmask)) call assert(size(locfield) == size(locmask), &
        'post_data_2d_low: mask size mismatch: '//trim(diag%debug_str))

  if (diag_cs%diag_as_chksum) then
    ! Append timestep to mesg
    call get_time(diag_cs%time_end, time_seconds, days=time_days)
    write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
        trim(diag%debug_str), time_days, time_seconds

    if (diag%axes%is_h_point) then
      call hchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_u_point) then
      call uchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_v_point) then
      call vchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    elseif (diag%axes%is_q_point) then
      call Bchksum(locfield, debug_mesg, diag_cs%G%HI, &
                   logunit=diag_cs%chksum_iounit)
    else
      call SIS_error(FATAL, "post_data_2d_low: unknown axis type.")
    endif
  else
    if (is_stat) then
      if (associated(locmask)) then
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, rmask=locmask)
      else
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev)
      endif
    elseif (diag_cs%ave_enabled) then
      if (associated(locmask)) then
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, &
                         time=diag_cs%time_end, weight=diag_cs%time_int, rmask=locmask)
      else
        used = send_data_infra(diag%fms_diag_id, locfield, &
                         is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, &
                         time=diag_cs%time_end, weight=diag_cs%time_int)
      endif
    endif
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) &
    deallocate( locfield )
end subroutine post_data_2d_low

!> Make a real 3-d array diagnostic available for averaging or output.
subroutine post_data_3d(diag_field_id, field, diag_cs, is_static, mask)
  integer,           intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                 !! previous call to register_SIS_diag_field.
  real,      target, intent(in) :: field(:,:,:)  !< 3-d array being offered for output or averaging
                                                 !! in internally scaled arbitrary units [A ~> a]
  type(SIS_diag_ctrl), target, intent(in) :: diag_cs !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables
  type(diag_type), pointer :: diag => null()

  if (id_clock_diag_mediator>0) call cpu_clock_begin(id_clock_diag_mediator)

  ! Iterate over list of diag 'variants', e.g. CMOR aliases, different vertical
  ! grids, and post each.
  call assert(diag_field_id < diag_cs%next_free_diag_id, &
              'post_data_3d: Unregistered diagnostic id')

  if (diag_cs%show_call_tree) &
    call callTree_enter("post_data_3d("//trim(diag_cs%diags(diag_field_id)%debug_str)//")")

  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    call assert(associated(diag%axes), 'post_data_3d: axes is not associated')

    call post_data_3d_low(diag, field, diag_cs, is_static, mask)

    diag => diag%next
  enddo
  if (id_clock_diag_mediator>0) call cpu_clock_end(id_clock_diag_mediator)

  if (diag_cs%show_call_tree) &
    call callTree_leave("post_data_3d("//trim(diag_cs%diags(diag_field_id)%debug_str)//")")

end subroutine post_data_3d

!> Make a real 3-d array diagnostic available for averaging or output
!! using a diag_type instead of an integer id.
subroutine post_data_3d_low(diag, field, diag_cs, is_static, mask)
  type(diag_type),   intent(in) :: diag       !< A structure describing the diagnostic to post
  real,    target,   intent(in) :: field(:,:,:) !< 3-d array being offered for output or averaging
                                                !! in internally scaled arbitrary units [A ~> a]
  type(SIS_diag_ctrl),   intent(in) :: diag_CS !< Structure used to regulate diagnostic output
  logical, optional, intent(in) :: is_static !< If true, this is a static field that is always offered.
  real,    optional,target, intent(in) :: mask(:,:,:) !< If present, use this real array as the data mask [nondim]

  ! Local variables
  real, dimension(:,:,:), pointer :: locfield ! The field being offered in arbitrary unscaled units [a]
  real, dimension(:,:,:), pointer :: locmask  ! A pointer to the data mask to use [nondim]
  character(len=300) :: mesg
  logical :: used  ! The return value of send_data is not used for anything.
  logical :: is_stat
  integer :: cszi, cszj, dszi, dszj
  integer :: isv, iev, jsv, jev, ks, ke, i, j, k

  integer :: time_days
  integer :: time_seconds
  character(len=300) :: debug_mesg

  locfield => NULL()
  locmask => NULL()
  is_stat = .false. ; if (present(is_static)) is_stat = is_static

  ! Determine the proper array indices, noting that because of the (:,:)
  ! declaration of field, symmetric arrays are using a SW-grid indexing,
  ! but non-symmetric arrays are using a NE-grid indexing.  Send_data
  ! actually only uses the difference between ie and is to determine
  ! the output data size and assumes that halos are symmetric.
  isv = diag_cs%is ; iev = diag_cs%ie ; jsv = diag_cs%js ; jev = diag_cs%je

  cszi = (diag_cs%ie-diag_cs%is) +1 ; dszi = (diag_cs%ied-diag_cs%isd) +1
  cszj = (diag_cs%je-diag_cs%js) +1 ; dszj = (diag_cs%jed-diag_cs%jsd) +1
  if ( size(field,1) == dszi ) then
    isv = diag_cs%is ; iev = diag_cs%ie     ! Data domain
  elseif ( size(field,1) == dszi + 1 ) then
    isv = diag_cs%is ; iev = diag_cs%ie+1   ! Symmetric data domain
  elseif ( size(field,1) == cszi ) then
    isv = 1 ; iev = cszi                    ! Computational domain
  elseif ( size(field,1) == cszi + 1 ) then
    isv = 1 ; iev = cszi+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,1)," in i-direction\n"//&
       "does not match one of ", cszi, cszi+1, dszi, dszi+1
    call SIS_error(FATAL,"post_data_3d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  if ( size(field,2) == dszj ) then
    jsv = diag_cs%js ; jev = diag_cs%je     ! Data domain
  elseif ( size(field,2) == dszj + 1 ) then
    jsv = diag_cs%js ; jev = diag_cs%je+1   ! Symmetric data domain
  elseif ( size(field,2) == cszj ) then
    jsv = 1 ; jev = cszj                    ! Computational domain
  elseif ( size(field,2) == cszj + 1 ) then
    jsv = 1 ; jev = cszj+1                  ! Symmetric computational domain
  else
    write (mesg,*) " peculiar size ",size(field,2)," in j-direction\n"//&
       "does not match one of ", cszj, cszj+1, dszj, dszj+1
    call SIS_error(FATAL,"post_data_3d_low: "//trim(diag%debug_str)//trim(mesg))
  endif

  ks = lbound(field,3) ; ke = ubound(field,3)
  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) then
    allocate( locfield( lbound(field,1):ubound(field,1), lbound(field,2):ubound(field,2), ks:ke ), source=0.0 )
    do k=ks,ke ; do j=jsv,jev ; do i=isv,iev
      if (field(i,j,k) == diag_cs%missing_value) then
        locfield(i,j,k) = diag_cs%missing_value
      else
        locfield(i,j,k) = field(i,j,k) * diag%conversion_factor
      endif
    enddo ; enddo ; enddo
  else
    locfield => field
  endif

  if (present(mask)) then
    locmask => mask
  elseif (.not.is_stat) then  ! Static fields do not have assigned axes.
    if (associated(diag%axes%mask3d)) then
      locmask => diag%axes%mask3d
    endif
  endif
  if (associated(locmask)) call assert(size(locfield) == size(locmask), &
        'post_data_3d_low: mask size mismatch: '//diag%debug_str)

  if (diag%fms_diag_id>0) then
    if (diag_cs%diag_as_chksum) then
      ! Append timestep to mesg
      call get_time(diag_cs%time_end, time_seconds, days=time_days)
      write(debug_mesg, '(a, 1x, i0, 1x, i0)') &
          trim(diag%debug_str), time_days, time_seconds

      if (diag%axes%is_h_point) then
        call hchksum(locfield, debug_mesg, diag_cs%G%HI, &
                     logunit=diag_cs%chksum_iounit)
      elseif (diag%axes%is_u_point) then
        call uchksum(locfield, debug_mesg, diag_cs%G%HI, &
                     logunit=diag_cs%chksum_iounit)
      elseif (diag%axes%is_v_point) then
        call vchksum(locfield, debug_mesg, diag_cs%G%HI, &
                     logunit=diag_cs%chksum_iounit)
      elseif (diag%axes%is_q_point) then
        call Bchksum(locfield, debug_mesg, diag_cs%G%HI, &
                     logunit=diag_cs%chksum_iounit)
      else
        call SIS_error(FATAL, "post_data_3d_low: unknown axis type.")
      endif
    else
      if (is_stat) then
        if (associated(locmask)) then
          used = send_data_infra(diag%fms_diag_id, locfield, &
                           is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, rmask=locmask)
        else
          used = send_data_infra(diag%fms_diag_id, locfield, &
                           is_in=isv, ie_in=iev, js_in=jsv, je_in=jev)
        endif
      elseif (diag_cs%ave_enabled) then
        if (associated(locmask)) then
          used = send_data_infra(diag%fms_diag_id, locfield, &
                           is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, &
                           time=diag_cs%time_end, weight=diag_cs%time_int, rmask=locmask)
        else
          used = send_data_infra(diag%fms_diag_id, locfield, &
                           is_in=isv, ie_in=iev, js_in=jsv, je_in=jev, &
                           time=diag_cs%time_end, weight=diag_cs%time_int)
        endif
      endif
    endif
  endif

  if ((diag%conversion_factor /= 0.) .and. (diag%conversion_factor /= 1.)) deallocate( locfield )

end subroutine post_data_3d_low

!> Enable the accumulation of time averages over the specified time interval.
subroutine enable_SIS_averaging(time_int_in, time_end_in, diag_cs)
  real,                intent(in)    :: time_int_in !< The time interval [s] over which any
                                                    !! values that are offered are valid.
  type(time_type),     intent(in)    :: time_end_in !< The end time of the valid interval
  type(SIS_diag_ctrl), intent(inout) :: diag_CS     !< Structure used to regulate diagnostic output
! This subroutine enables the accumulation of time averages over the specified time interval.

!  if (num_file==0) return
  diag_cs%time_int = time_int_in
  diag_cs%time_end = time_end_in
  diag_cs%ave_enabled = .true.
end subroutine enable_SIS_averaging

!> Enable the accumulation of time averages over the specified time interval in time units.
subroutine enable_SIS_averages(time_int, time_end, diag_CS, T_to_s)
  real,            intent(in)    :: time_int !< The time interval over which any values
                                             !! that are offered are valid [T ~> s].
  type(time_type), intent(in)    :: time_end !< The end time of the valid interval.
  type(SIS_diag_ctrl), intent(inout) :: diag_CS  !< A structure that is used to regulate diagnostic output
  real,  optional, intent(in)    :: T_to_s   !< A conversion factor for time_int to [s].
! This subroutine enables the accumulation of time averages over the specified time interval.

  if (present(T_to_s)) then
    diag_cs%time_int = time_int*T_to_s
  elseif (associated(diag_CS%US)) then
    diag_cs%time_int = time_int*diag_CS%US%T_to_s
  else
    diag_cs%time_int = time_int
  endif
  diag_cs%time_end = time_end
  diag_cs%ave_enabled = .true.
end subroutine enable_SIS_averages

!> Call this subroutine to avoid averaging any offered fields.
subroutine disable_SIS_averaging(diag_cs)
  type(SIS_diag_ctrl), intent(inout) :: diag_cs !< Structure used to regulate diagnostic output

  diag_cs%time_int = 0.0
  diag_cs%ave_enabled = .false.
end subroutine disable_SIS_averaging

!> Indicate whether averaging diagnostics is currently enabled
logical function query_SIS_averaging_enabled(diag_cs, time_int, time_end)
  type(SIS_diag_ctrl),       intent(in)  :: diag_cs  !< Structure used to regulate diagnostic output
  real,            optional, intent(out) :: time_int !< Current setting of diag_cs%time_int [s]
  type(time_type), optional, intent(out) :: time_end !< Current setting of diag_cs%time_end

  if (present(time_int)) time_int = diag_cs%time_int
  if (present(time_end)) time_end = diag_cs%time_end
  query_SIS_averaging_enabled = diag_cs%ave_enabled
end function query_SIS_averaging_enabled

!> This function returns the valid end time for use with diagnostics that are
!! handled outside of the MOM6 diagnostics infrastructure.
function get_SIS_diag_time_end(diag_cs)
  type(SIS_diag_ctrl), intent(in)  :: diag_CS !< Structure used to regulate diagnostic output
  type(time_type) :: get_SIS_diag_time_end
  !   This function returns the valid end time for diagnostics that are handled
  ! outside of the MOM6 infrastructure, such as via the generic tracer code.

  get_SIS_diag_time_end = diag_cs%time_end
end function get_SIS_diag_time_end

!> Returns the "diag_mediator" handle for a group (native, CMOR, ...) of diagnostics
!! derived from one field.
integer function register_SIS_diag_field(module_name, field_name, axes_in, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name, &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, conversion)
  character(len=*),           intent(in) :: module_name !< Name of this module, usually "ice_model"
  character(len=*),           intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp),     target, intent(in) :: axes_in   !< Container with up to 3 (or 4?) integer handles that
                                                      !! indicates axes for this field
  type(time_type),            intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_data calls (not used in SIS?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in SIS?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in SIS?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in SIS?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count   !< no clue (not used in SIS?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute. Use '' to
                                                         !! have no attribute.  If present, this overrides the
                                                         !! default constructed from the default for
                                                         !! each individual axis direction.
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables
  real :: SIS_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(SIS_diag_ctrl), pointer :: diag_cs => NULL() ! A structure that is used
                                               ! to regulate diagnostic output
  type(axes_grp), pointer :: axes
  integer :: dm_id
  character(len=256) :: msg, cm_string
  character(len=256) :: new_module_name
  character(len=480) :: module_list, var_list
  character(len=24)  :: dimensions
  integer :: num_modnm, num_varnm
  logical :: active

  diag_cs => axes_in%diag_cs

  ! Check if the axes match a standard grid axis.
  ! If not, allocate the new axis and copy the contents.
  if (axes_in%id == diag_cs%axesTL%id) then
    axes => diag_cs%axesTL
  elseif (axes_in%id == diag_cs%axesBL%id) then
    axes => diag_cs%axesBL
  elseif (axes_in%id == diag_cs%axesCuL%id) then
    axes => diag_cs%axesCuL
  elseif (axes_in%id == diag_cs%axesCvL%id) then
    axes => diag_cs%axesCvL
  elseif (axes_in%id == diag_cs%axesTi%id) then
    axes => diag_cs%axesTi
  elseif (axes_in%id == diag_cs%axesBi%id) then
    axes => diag_cs%axesBi
  elseif (axes_in%id == diag_cs%axesCui%id) then
    axes => diag_cs%axesCui
  elseif (axes_in%id == diag_cs%axesCvi%id) then
    axes => diag_cs%axesCvi
  elseif (axes_in%id == diag_cs%axesTC%id) then
    axes => diag_cs%axesTC
  elseif (axes_in%id == diag_cs%axesCuC%id) then
    axes => diag_cs%axesCuC
  elseif (axes_in%id == diag_cs%axesCvC%id) then
    axes => diag_cs%axesCvC
  elseif (axes_in%id == diag_cs%axesBC%id) then
    axes => diag_cs%axesBC
  elseif (axes_in%id == diag_cs%axesTC0%id) then
    axes => diag_cs%axesTC0
  elseif (axes_in%id == diag_cs%axesCuC0%id) then
    axes => diag_cs%axesCuC0
  elseif (axes_in%id == diag_cs%axesCvC0%id) then
    axes => diag_cs%axesCvC0
  elseif (axes_in%id == diag_cs%axesBC0%id) then
    axes => diag_cs%axesBC0
  else
    allocate(axes)
    axes = axes_in
  endif

  SIS_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) SIS_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1

  module_list = "{"//trim(module_name)
  num_modnm = 1

  ! Register the native diagnostic
  active = register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, &
             init_time, long_name=long_name, units=units, missing_value=SIS_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count, &
             cmor_field_name=cmor_field_name, cmor_long_name=cmor_long_name, &
             cmor_units=cmor_units, cmor_standard_name=cmor_standard_name, &
             cell_methods=cell_methods, x_cell_method=x_cell_method, y_cell_method=y_cell_method, &
             conversion=conversion)
  num_varnm = 1 ; var_list = "{"//trim(field_name)
  if (present(cmor_field_name)) then
    num_varnm = num_varnm + 1
    var_list = trim(var_list)//","//trim(cmor_field_name)
  endif
  var_list = trim(var_list)//"}"

  dimensions = ""
  if (axes_in%is_h_point)   dimensions = trim(dimensions)//" xh, yh,"
  if (axes_in%is_q_point)   dimensions = trim(dimensions)//" xq, yq,"
  if (axes_in%is_u_point)   dimensions = trim(dimensions)//" xq, yh,"
  if (axes_in%is_v_point)   dimensions = trim(dimensions)//" xh, yq,"
  if (axes_in%is_category)  dimensions = trim(dimensions)//" ct,"
  if (axes_in%is_cat_open)  dimensions = trim(dimensions)//" ctu,"
  if (axes_in%is_layer)     dimensions = trim(dimensions)//" zl,"
  if (axes_in%is_interface) dimensions = trim(dimensions)//" zi,"

  if (len_trim(dimensions) > 0) then
    dimensions = trim(adjustl(dimensions))
    if (dimensions(len_trim(dimensions):len_trim(dimensions)) == ",") then
        dimensions = dimensions(1:len_trim(dimensions) - 1)
    endif
    dimensions = trim(dimensions)
  endif

  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit > 0)) then
    msg = ''
    if (present(cmor_field_name)) msg = 'CMOR equivalent is "'//trim(cmor_field_name)//'"'
    cm_string = ''
    !### Uncoment this to add cell methods:
    ! call attach_cell_methods(-1, axes, cm_string, cell_methods, x_cell_method, y_cell_method)
    module_list = trim(module_list)//"}"
    if (num_modnm <= 1) module_list = module_name
    if (num_varnm <= 1) var_list = ''

    call log_available_diag(dm_id>0, module_list, field_name, cm_string, msg, diag_CS, &
                            long_name, units, standard_name, variants=var_list, dimensions=dimensions)
  endif

  register_SIS_diag_field = dm_id

end function register_SIS_diag_field

!> Returns True if either the native or CMOR version of the diagnostic were registered. Updates 'dm_id'
!! after calling register_diag_field_expand_axes() for both native and CMOR variants of the field.
logical function register_diag_field_expand_cmor(dm_id, module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,      &
            verbose, do_not_log, err_msg, interp_method, tile_count, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, cell_methods, &
            x_cell_method, y_cell_method, conversion)
  integer,          intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ice_model" or "ice_model_fast"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp),   intent(in) :: axes !< Container w/ up to 3 integer handles that indicates axes
                                             !! for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided
                                                         !! with post_data calls (not used in SIS?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in SIS?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in SIS?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in SIS?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in SIS?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  character(len=*), optional, intent(in) :: cell_methods !< String to append as cell_methods attribute.
                                                         !! Use '' to have no attribute. If present, this
                                                         !! overrides the default constructed from the default
                                                         !! for each individual axis direction.
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]
  ! Local variables
  real :: SIS_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(SIS_diag_ctrl), pointer :: diag_cs => null()
  type(diag_type), pointer :: this_diag => null()
  integer :: fms_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name, cm_string

  SIS_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) SIS_missing_value = missing_value

  register_diag_field_expand_cmor = .false.
  diag_cs => axes%diag_cs

  ! Set up the 'primary' diagnostic, first get an underlying FMS id
  fms_id = register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
             long_name=long_name, units=units, missing_value=SIS_missing_value, &
             range=range, mask_variant=mask_variant, standard_name=standard_name, &
             verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
             interp_method=interp_method, tile_count=tile_count)
  if (.not. diag_cs%diag_as_chksum) &
    cm_string = ''
    !### Uncoment this to add cell methods:
    ! call attach_cell_methods(fms_id, axes, cm_string, cell_methods, x_cell_method, y_cell_method)

  this_diag => null()
  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
    if (present(conversion)) this_diag%conversion_factor = conversion
    register_diag_field_expand_cmor = .true.
  endif

  ! For the CMOR variation of the above diagnostic
  if (present(cmor_field_name) .and. .not. diag_cs%diag_as_chksum) then
    ! Fallback values for strings set to "NULL"
    posted_cmor_units = "not provided"         !
    posted_cmor_standard_name = "not provided" ! Values might be able to be replaced with a CS%missing field?
    posted_cmor_long_name = "not provided"     !

    ! If attributes are present for MOM variable names, use them first for the register_SIS_diag_field
    ! call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_SIS_diag_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_diag_field_expand_axes(module_name, cmor_field_name, axes, init_time,    &
               long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units),                  &
               missing_value=SIS_missing_value, range=range, mask_variant=mask_variant,               &
               standard_name=trim(posted_cmor_standard_name), verbose=verbose, do_not_log=do_not_log, &
               err_msg=err_msg, interp_method=interp_method, tile_count=tile_count)
    cm_string = ''
    !### Uncoment this to add cell methods:
    ! call attach_cell_methods(fms_id, axes, cm_string, cell_methods, x_cell_method, y_cell_method)

    this_diag => null()
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      call add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
      if (present(conversion)) this_diag%conversion_factor = conversion
      register_diag_field_expand_cmor = .true.
    endif
  endif

end function register_diag_field_expand_cmor

!> Returns an FMS id from register_diag_field_fms (the diag_manager routine) after expanding axes
!! (axes-group) into handles and conditionally adding an FMS area_id for cell_measures.
integer function register_diag_field_expand_axes(module_name, field_name, axes, init_time, &
            long_name, units, missing_value, range, mask_variant, standard_name,  &
            verbose, do_not_log, err_msg, interp_method, tile_count)
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ice_model"
                                              !! or "ice_model_fast"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container w/ up to 3 integer handles that indicates
                                             !! axes for this field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided
                                                         !! with post_data calls (not used in SIS?)
  logical,          optional, intent(in) :: verbose !< If true, FMS is verbose (not used in SIS?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something
                                                       !! (not used in SIS?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in SIS?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should
                                                         !! not be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count !< no clue (not used in SIS?)
  ! Local variables
  integer :: fms_id, area_id

  ! This gets the cell area associated with the grid location of this variable
  area_id = axes%id_area

  ! Get the FMS diagnostic id
  if (axes%diag_cs%diag_as_chksum) then
    fms_id = axes%diag_cs%num_chksum_diags + 1
    axes%diag_cs%num_chksum_diags = fms_id
  elseif (present(interp_method) .or. axes%is_h_point) then
    ! If interp_method is provided we must use it
    if (area_id>0) then
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method=interp_method, tile_count=tile_count, area=area_id)
    else
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method=interp_method, tile_count=tile_count)
    endif
  else
    ! If interp_method is not provided and the field is not at an h-point then interp_method='none'
    if (area_id>0) then
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method='none', tile_count=tile_count, area=area_id)
    else
      fms_id = register_diag_field_infra(module_name, field_name, axes%handles, &
               init_time, long_name=long_name, units=units, missing_value=missing_value, &
               range=range, mask_variant=mask_variant, standard_name=standard_name, &
               verbose=verbose, do_not_log=do_not_log, err_msg=err_msg, &
               interp_method='none', tile_count=tile_count)
    endif
  endif

  register_diag_field_expand_axes = fms_id

end function register_diag_field_expand_axes

!> Create a diagnostic type and attached to list
subroutine add_diag_to_list(diag_cs, dm_id, fms_id, this_diag, axes, module_name, field_name)
  type(SIS_diag_ctrl),        pointer       :: diag_cs !< Diagnostics mediator control structure
  integer,                intent(inout) :: dm_id !< The diag_mediator ID for this diagnostic group
  integer,                intent(in)    :: fms_id !< The FMS diag_manager ID for this diagnostic
  type(diag_type),        pointer       :: this_diag !< This diagnostic
  type(axes_grp), target, intent(in)    :: axes !< Container w/ up to 3 integer handles that
                                                !! indicates axes for this field
  character(len=*),       intent(in)    :: module_name !< Name of this module, usually
                                                       !! "ice_model" or "ice_model_fast"
  character(len=*),       intent(in)    :: field_name !< Name of diagnostic

  ! If the diagnostic is needed obtain a diag_mediator ID (if needed)
  if (dm_id == -1) dm_id = get_new_diag_id(diag_cs)
  ! Create a new diag_type to store links in
  call alloc_diag_with_id(dm_id, diag_cs, this_diag)
  call assert(associated(this_diag), 'add_diag_to_list: allocation failed for '//trim(field_name))
  ! Record FMS id, masks and conversion factor, in diag_type
  this_diag%fms_diag_id = fms_id
  this_diag%debug_str = trim(module_name)//"-"//trim(field_name)
  this_diag%axes => axes

end subroutine add_diag_to_list


!> Attaches "cell_methods" attribute to a variable based on defaults for axes_grp or optional arguments.
subroutine attach_cell_methods(id, axes, ostring, cell_methods, x_cell_method, y_cell_method)
  integer,                    intent(in)  :: id !< Handle to diagnostic
  type(axes_grp),             intent(in)  :: axes !< Container w/ up to 3 integer handles that indicates
                                                  !! axes for this field
  character(len=*),           intent(out) :: ostring !< The cell_methods strings that would appear in the file
  character(len=*), optional, intent(in)  :: cell_methods !< String to append as cell_methods attribute.
                                                         !! Use '' to have no attribute. If present, this
                                                         !! overrides the default constructed from the default
                                                         !! for each individual axis direction.
  character(len=*), optional, intent(in)  :: x_cell_method !< Specifies the cell method for the x-direction.
                                                         !! Use '' have no method.
  character(len=*), optional, intent(in)  :: y_cell_method !< Specifies the cell method for the y-direction.
                                                         !! Use '' have no method.
  ! Local variables
  character(len=9) :: axis_name
  logical :: x_mean, y_mean, x_sum, y_sum

  x_mean = .false.
  y_mean = .false.
  x_sum = .false.
  y_sum = .false.

  ostring = ''
  if (present(cell_methods)) then
    if (present(x_cell_method) .or. present(y_cell_method)) then
      call SIS_error(FATAL, "attach_cell_methods: " // &
           'Individual direction cell method was specified along with a "cell_methods" string.')
    endif
    if (len(trim(cell_methods))>0) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', trim(cell_methods))
      ostring = trim(cell_methods)
    endif
  else
    if (present(x_cell_method)) then
      if (len(trim(x_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(x_cell_method)
        if (trim(x_cell_method)=='mean') x_mean=.true.
        if (trim(x_cell_method)=='sum') x_sum=.true.
      endif
    else
      if (len(trim(axes%x_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%x_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%x_cell_method)
        if (trim(axes%x_cell_method)=='mean') x_mean=.true.
        if (trim(axes%x_cell_method)=='sum') x_sum=.true.
      endif
    endif
    if (present(y_cell_method)) then
      if (len(trim(y_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(y_cell_method)
        if (trim(y_cell_method)=='mean') y_mean=.true.
        if (trim(y_cell_method)=='sum') y_sum=.true.
      endif
    else
      if (len(trim(axes%y_cell_method))>0) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(id, 'cell_methods', trim(axis_name)//':'//trim(axes%y_cell_method))
        ostring = trim(adjustl(ostring))//' '//trim(axis_name)//':'//trim(axes%y_cell_method)
        if (trim(axes%y_cell_method)=='mean') y_mean=.true.
        if (trim(axes%y_cell_method)=='sum') y_sum=.true.
      endif
    endif
    if (x_mean .and. y_mean) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', 'area:mean')
      ostring = trim(adjustl(ostring))//' area:mean'
    elseif (x_sum .and. y_sum) then
      call MOM_diag_field_add_attribute(id, 'cell_methods', 'area:sum')
      ostring = trim(adjustl(ostring))//' area:sum'
    endif
  endif
  ostring = adjustl(ostring)
end subroutine attach_cell_methods

!> Registers a non-array scalar diagnostic, returning an integer handle
function register_scalar_field(module_name, field_name, init_time, diag_cs, &
            long_name, units, missing_value, range, standard_name, &
            do_not_log, err_msg, interp_method, cmor_field_name, &
            cmor_long_name, cmor_units, cmor_standard_name, conversion)
  integer :: register_scalar_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ice_model"
                                              !! or "ice_model_fast"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(time_type),  intent(in) :: init_time !< Time at which a field is first available?
  type(SIS_diag_ctrl),  intent(inout) :: diag_CS !< Structure used to regulate diagnostic output
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable in arbitrary units [a]
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in SIS?)
  character(len=*), optional, intent(out):: err_msg !< String into which an error message might be
                                                         !! placed (not used in SIS?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables
  real :: SIS_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  integer :: dm_id, fms_id
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=16)  :: dimensions

  SIS_missing_value = diag_cs%missing_value
  if (present(missing_value)) SIS_missing_value = missing_value

  dm_id = -1
  diag => null()
  cmor_diag => null()

  if (diag_cs%diag_as_chksum) then
    fms_id = diag_cs%num_chksum_diags + 1
    diag_cs%num_chksum_diags = fms_id
  else
    fms_id = register_diag_field_infra(module_name, field_name, init_time, &
                long_name=long_name, units=units, missing_value=SIS_missing_value, &
                range=range, standard_name=standard_name, do_not_log=do_not_log, &
                err_msg=err_msg)
  endif

  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    dm_id = get_new_diag_id(diag_cs)
    call alloc_diag_with_id(dm_id, diag_cs, diag)
    call assert(associated(diag), 'register_scalar_field: diag allocation failed')
    diag%fms_diag_id = fms_id
    diag%debug_str = trim(module_name)//"-"//trim(field_name)
    if (present(conversion)) diag%conversion_factor = conversion
  endif

  if (present(cmor_field_name)) then
    ! Fallback values for strings set to "not provided"
    posted_cmor_units = "not provided"
    posted_cmor_standard_name = "not provided"
    posted_cmor_long_name = "not provided"

    ! If attributes are present for MOM variable names, use them first for the register_static_field
    ! call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_static_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_diag_field_infra(module_name, cmor_field_name, init_time, &
           long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units), &
           missing_value=SIS_missing_value, range=range, &
           standard_name=trim(posted_cmor_standard_name), do_not_log=do_not_log, err_msg=err_msg)
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      if (dm_id == -1) then
        dm_id = get_new_diag_id(diag_cs)
      endif
      call alloc_diag_with_id(dm_id, diag_cs, cmor_diag)
      cmor_diag%fms_diag_id = fms_id
      cmor_diag%debug_str = trim(module_name)//"-"//trim(cmor_field_name)
      if (present(conversion)) cmor_diag%conversion_factor = conversion
    endif
  endif

  dimensions = "scalar"

  ! Document diagnostics in list of available diagnostics
  if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
    if (present(cmor_field_name)) then
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, &
                              variants="{"//trim(field_name)//","//trim(cmor_field_name)//"}", &
                              dimensions=dimensions)
    else
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, dimensions=dimensions)
    endif
  endif

  register_scalar_field = dm_id

end function register_scalar_field

!> Registers a static diagnostic, returning an integer handle
function register_static_field(module_name, field_name, axes, &
            long_name, units, missing_value, range, mask_variant, standard_name, &
            do_not_log, interp_method, tile_count, &
            cmor_field_name, cmor_long_name, cmor_units, cmor_standard_name, area, &
            x_cell_method, y_cell_method, area_cell_method, conversion)
  integer :: register_static_field !< An integer handle for a diagnostic array.
  character(len=*), intent(in) :: module_name !< Name of this module, usually "ice_model"
  character(len=*), intent(in) :: field_name !< Name of the diagnostic field
  type(axes_grp), target, intent(in) :: axes !< Container with up to 3 integer handles that
                                             !! indicates axes for this field
  character(len=*), optional, intent(in) :: long_name !< Long name of a field.
  character(len=*), optional, intent(in) :: units !< Units of a field.
  character(len=*), optional, intent(in) :: standard_name !< Standardized name associated with a field
  real,             optional, intent(in) :: missing_value !< A value that indicates missing values in
                                                          !! output files, in unscaled arbitrary units [a]
  real,             optional, intent(in) :: range(2) !< Valid range of a variable in arbitrary units [a]
  logical,          optional, intent(in) :: mask_variant !< If true a logical mask must be provided with
                                                         !! post_data calls (not used in SIS?)
  logical,          optional, intent(in) :: do_not_log !< If true, do not log something (not used in SIS?)
  character(len=*), optional, intent(in) :: interp_method !< If 'none' indicates the field should not
                                                         !! be interpolated as a scalar
  integer,          optional, intent(in) :: tile_count   !< no clue (not used in SIS?)
  character(len=*), optional, intent(in) :: cmor_field_name !< CMOR name of a field
  character(len=*), optional, intent(in) :: cmor_long_name !< CMOR long name of a field
  character(len=*), optional, intent(in) :: cmor_units !< CMOR units of a field
  character(len=*), optional, intent(in) :: cmor_standard_name !< CMOR standardized name associated with a field
  integer,          optional, intent(in) :: area !< fms_id for area_t
  character(len=*), optional, intent(in) :: x_cell_method !< Specifies the cell method for the x-direction.
  character(len=*), optional, intent(in) :: y_cell_method !< Specifies the cell method for the y-direction.
  character(len=*), optional, intent(in) :: area_cell_method !< Specifies the cell method for area
  real,             optional, intent(in) :: conversion !< A value to multiply data by before writing to files,
                                                       !! often including factors to undo internal scaling and
                                                       !! in units of [a A-1 ~> 1]

  ! Local variables
  real :: SIS_missing_value ! A value used to indicate missing values in output files, in arbitrary units [a]
  type(SIS_diag_ctrl), pointer :: diag_cs => null() !< A structure that is used to regulate diagnostic output
  type(diag_type), pointer :: diag => null(), cmor_diag => null()
  integer :: dm_id, fms_id
  character(len=256) :: posted_cmor_units, posted_cmor_standard_name, posted_cmor_long_name
  character(len=9) :: axis_name
  character(len=24) :: dimensions

  SIS_missing_value = axes%diag_cs%missing_value
  if (present(missing_value)) SIS_missing_value = missing_value

  diag_cs => axes%diag_cs
  dm_id = -1
  diag => null()
  cmor_diag => null()

  if (diag_cs%diag_as_chksum) then
    fms_id = diag_cs%num_chksum_diags + 1
    diag_cs%num_chksum_diags = fms_id
  else
    fms_id = register_static_field_infra(module_name, field_name, axes%handles, &
           long_name=long_name, units=units, missing_value=SIS_missing_value, &
           range=range, mask_variant=mask_variant, standard_name=standard_name, &
           do_not_log=do_not_log, &
           interp_method=interp_method, tile_count=tile_count, area=area)
  endif

  if (fms_id /= DIAG_FIELD_NOT_FOUND) then
    dm_id = get_new_diag_id(diag_cs)
    call alloc_diag_with_id(dm_id, diag_cs, diag)
    call assert(associated(diag), 'register_static_field: diag allocation failed')
    diag%fms_diag_id = fms_id
    diag%debug_str = trim(module_name)//"-"//trim(field_name)
    if (present(conversion)) diag%conversion_factor = conversion

    if (diag_cs%diag_as_chksum) then
      diag%axes => axes
    else
      if (present(x_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            trim(axis_name)//':'//trim(x_cell_method))
      endif
      if (present(y_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            trim(axis_name)//':'//trim(y_cell_method))
      endif
      if (present(area_cell_method)) then
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', &
            'area:'//trim(area_cell_method))
      endif
    endif
  endif

  if (present(cmor_field_name) .and. .not. diag_cs%diag_as_chksum) then
    ! Fallback values for strings set to "not provided"
    posted_cmor_units = "not provided"
    posted_cmor_standard_name = "not provided"
    posted_cmor_long_name = "not provided"

    ! If attributes are present for MOM variable names, use them first for the register_static_field
    ! call for CMOR verison of the variable
    if (present(units)) posted_cmor_units = units
    if (present(standard_name)) posted_cmor_standard_name = standard_name
    if (present(long_name)) posted_cmor_long_name = long_name

    ! If specified in the call to register_static_field, override attributes with the CMOR versions
    if (present(cmor_units)) posted_cmor_units = cmor_units
    if (present(cmor_standard_name)) posted_cmor_standard_name = cmor_standard_name
    if (present(cmor_long_name)) posted_cmor_long_name = cmor_long_name

    fms_id = register_static_field_infra(module_name, cmor_field_name, axes%handles, &
                long_name=trim(posted_cmor_long_name), units=trim(posted_cmor_units), &
                missing_value=SIS_missing_value, range=range, mask_variant=mask_variant, &
                standard_name=trim(posted_cmor_standard_name), do_not_log=do_not_log, &
                interp_method=interp_method, tile_count=tile_count, area=area)
    if (fms_id /= DIAG_FIELD_NOT_FOUND) then
      if (dm_id == -1) then
        dm_id = get_new_diag_id(diag_cs)
      endif
      call alloc_diag_with_id(dm_id, diag_cs, cmor_diag)
      cmor_diag%fms_diag_id = fms_id
      cmor_diag%debug_str = trim(module_name)//"-"//trim(cmor_field_name)
      if (present(conversion)) cmor_diag%conversion_factor = conversion
      if (present(x_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(1), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(x_cell_method))
      endif
      if (present(y_cell_method)) then
        call get_MOM_diag_axis_name(axes%handles(2), axis_name)
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', trim(axis_name)//':'//trim(y_cell_method))
      endif
      if (present(area_cell_method)) then
        call MOM_diag_field_add_attribute(fms_id, 'cell_methods', 'area:'//trim(area_cell_method))
      endif
    endif
  endif

  dimensions = ""
  if (axes%is_h_point)   dimensions = trim(dimensions)//" xh, yh,"
  if (axes%is_q_point)   dimensions = trim(dimensions)//" xq, yq,"
  if (axes%is_u_point)   dimensions = trim(dimensions)//" xq, yh,"
  if (axes%is_v_point)   dimensions = trim(dimensions)//" xh, yq,"
  if (axes%is_category)  dimensions = trim(dimensions)//" ct,"
  if (axes%is_cat_open)  dimensions = trim(dimensions)//" ctu,"
  if (axes%is_layer)     dimensions = trim(dimensions)//" zl,"
  if (axes%is_interface) dimensions = trim(dimensions)//" zi,"

  if (len_trim(dimensions) > 0) then
    dimensions = trim(adjustl(dimensions))
    if (dimensions(len_trim(dimensions):len_trim(dimensions)) == ",") then
        dimensions = dimensions(1:len_trim(dimensions) - 1)
    endif
    dimensions = trim(dimensions)
  endif

  ! Document diagnostics in list of available diagnostics
  !### if (is_root_pe() .and. diag_CS%available_diag_doc_unit > 0) then
  if (is_root_pe() .and. .false.) then  ! Replace this to work like MOM6.
    if (present(cmor_field_name)) then
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, &
                              variants="{"//trim(field_name)//","//trim(cmor_field_name)//"}", &
                              dimensions=dimensions)
    else
      call log_available_diag(associated(diag), module_name, field_name, '', '', diag_CS, &
                              long_name, units, standard_name, dimensions=dimensions)
    endif
  endif

  register_static_field = dm_id

end function register_static_field

!> Add a description of an option to the documentation file
subroutine describe_option(opt_name, value, diag_CS)
  character(len=*),    intent(in) :: opt_name !< The name of the option
  character(len=*),    intent(in) :: value    !< The value of the option
  type(SIS_diag_ctrl), intent(in) :: diag_CS  !< Structure used to regulate diagnostic output

  ! Local variables
  character(len=480) :: mesg
  integer :: len_ind

  len_ind = len_trim(value)

  mesg = "    ! "//trim(opt_name)//": "//trim(value)
  write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
end subroutine describe_option

!> Initialize the SIS diag_mediator and opens the available diagnostics file, if appropriate.
subroutine SIS_diag_mediator_init(G, IG, US, param_file, diag_cs, component, err_msg, &
                                  doc_file_dir)
  type(SIS_hor_grid_type), target, intent(inout) :: G  !< The horizontal grid type
  type(ice_grid_type),             intent(in)    :: IG !< The sea-ice specific grid type
  type(unit_scale_type),   target, intent(in)    :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl),        intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output
  character(len=*), optional, intent(in)    :: component !< An optional component name
  character(len=*), optional, intent(out)   :: err_msg !< A string for a returned error message
  character(len=*), optional, intent(in)    :: doc_file_dir !< A directory in which to create the file

  ! This subroutine initializes the diag_mediator and the diag_manager.
  ! The grid type should have its dimensions set by this point, but it
  ! is not necessary that the metrics and axis labels be set up yet.

  ! Local variables
  integer :: ios, i, new_unit
  logical :: opened, new_file
  character(len=8)   :: this_pe
  character(len=240) :: doc_file, doc_file_dflt, doc_path
  character(len=40)  :: doc_file_param
  character(len=240), allocatable :: diag_coords(:)
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40) :: mdl = "SIS_diag_mediator" ! This module's name.
  character(len=32) :: filename_appendix = '' !fms appendix to filename for ensemble runs

  call diag_manager_init(err_msg=err_msg)

  id_clock_diag_mediator = cpu_clock_id('(Ocean diagnostics framework)', grain=CLOCK_MODULE)

  ! Allocate and initialize list of all diagnostics (and variants)
  allocate(diag_cs%diags(DIAG_ALLOC_CHUNK_SIZE))
  diag_cs%next_free_diag_id = 1
  do i=1, DIAG_ALLOC_CHUNK_SIZE
    call initialize_diag_type(diag_cs%diags(i))
  enddo

  diag_cs%show_call_tree = callTree_showQuery()

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, 'USE_INDEX_DIAGNOSTIC_AXES', diag_cs%index_space_axes, &
                 'If true, use a grid index coordinate convention for diagnostic axes. ',&
                 default=.false.)

  call get_param(param_file, mdl, 'DIAG_MISVAL', diag_cs%missing_value, &
                 'Set the default missing value to use for diagnostics.', &
                 units="various", default=1.e20)
  call get_param(param_file, mdl, 'DIAG_AS_CHKSUM', diag_cs%diag_as_chksum, &
                 'Instead of writing diagnostics to the diag manager, write '//&
                 'a text file containing the checksum (bitcount) of the array.',  &
                 default=.false.)

  if (diag_cs%diag_as_chksum) &
    diag_cs%num_chksum_diags = 0

  ! Keep pointers to the grid for diagnostic checksums
  diag_cs%G => G
  diag_cs%US => US

  diag_cs%is = G%isc - (G%isd-1) ; diag_cs%ie = G%iec - (G%isd-1)
  diag_cs%js = G%jsc - (G%jsd-1) ; diag_cs%je = G%jec - (G%jsd-1)
  diag_cs%isd = G%isd ; diag_cs%ied = G%ied
  diag_cs%jsd = G%jsd ; diag_cs%jed = G%jed

  ! Initialize available diagnostic log file
  if (is_root_pe() .and. (diag_CS%available_diag_doc_unit < 0)) then
    if (present(component)) then
      doc_file_dflt = trim(component)//".available_diags"
      doc_file_param = trim(uppercase(component))//"_AVAILABLE_DIAGS_FILE"
    else
      write(this_pe,'(i6.6)') PE_here()
      doc_file_dflt = "available_diags."//this_pe
      doc_file_param = "AVAILABLE_DIAGS_FILE"
    endif
    call get_param(param_file, mdl, trim(doc_file_param), doc_file, &
                 "A file into which to write a list of all available "//&
                 "sea ice diagnostics that can be included in a diag_table.", &
                 default=doc_file_dflt, do_not_log=(diag_CS%available_diag_doc_unit/=-1))
    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%available_diag_doc_unit /= -1) new_file = .false.
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

      diag_CS%available_diag_doc_unit = new_unit

      if (new_file) then
        open(diag_CS%available_diag_doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%available_diag_doc_unit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%available_diag_doc_unit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call SIS_error(FATAL, "Failed to open available diags file "//trim(doc_path)//".")
      endif
    endif
  endif

  if (is_root_pe() .and. (diag_CS%chksum_iounit < 0) .and. diag_CS%diag_as_chksum) then
    !write(this_pe,'(i6.6)') PE_here()
    !doc_file_dflt = "chksum_diag."//this_pe
    doc_file_dflt = "chksum_diag"
    call get_param(param_file, mdl, "CHKSUM_DIAG_FILE", doc_file, &
                 "A file into which to write all checksums of the "//&
                 "diagnostics listed in the diag_table.", &
                 default=doc_file_dflt, do_not_log=(diag_CS%chksum_iounit/=-1))

    call get_filename_appendix(filename_appendix)
    if (len_trim(filename_appendix) > 0) then
      doc_file = trim(doc_file) //'.'//trim(filename_appendix)
    endif
#ifdef STATSLABEL
    doc_file = trim(doc_file)//"."//trim(adjustl(STATSLABEL))
#endif

    if (len_trim(doc_file) > 0) then
      new_file = .true. ; if (diag_CS%chksum_iounit /= -1) new_file = .false.
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

      diag_CS%chksum_iounit = new_unit

      if (new_file) then
        open(diag_CS%chksum_iounit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='REPLACE', iostat=ios)
      else ! This file is being reopened, and should be appended.
        open(diag_CS%chksum_iounit, file=trim(doc_path), access='SEQUENTIAL', form='FORMATTED', &
             action='WRITE', status='OLD', position='APPEND', iostat=ios)
      endif
      inquire(diag_CS%chksum_iounit, opened=opened)
      if ((.not.opened) .or. (ios /= 0)) then
        call SIS_error(FATAL, "Failed to open checksum diags file "//trim(doc_path)//".")
      endif
    endif
  endif

  call diag_masks_set(G, IG, -1.0e34, diag_cs)

end subroutine SIS_diag_mediator_init

!> Sets up the 2d and 3d masks for native diagnostics
subroutine diag_masks_set(G, IG, missing_value, diag_cs)
  type(SIS_hor_grid_type), target, intent(in)    :: G   !< The horizontal grid type
  type(ice_grid_type),             intent(in)    :: IG  !< The sea-ice specific grid type
  real,                            intent(in)    :: missing_value !< A fill value for missing points
  type(SIS_diag_ctrl),             intent(inout) :: diag_cs !< A structure that is used to regulate diagnostic output

  ! Local variables
  integer :: i, j, k, NkIce, CatIce

  NkIce = IG%NkIce ; CatIce = IG%CatIce

  ! 2d masks point to the model masks since they are identical
  diag_cs%mask2dT  => G%mask2dT
  diag_cs%mask2dBu => G%mask2dBu
  diag_cs%mask2dCu => G%mask2dCu
  diag_cs%mask2dCv => G%mask2dCv

  allocate(diag_cs%mask2dT_comp(G%isc:G%iec,G%jsc:G%jec))
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    diag_cs%mask2dT_comp(i,j) = diag_cs%mask2dT(i,j)
  enddo ; enddo

  allocate(diag_cs%mask3dTL(G%isd:G%ied,G%jsd:G%jed,1:NkIce))
  allocate(diag_cs%mask3dBL(G%IsdB:G%IedB,G%JsdB:G%JedB,1:NkIce))
  allocate(diag_cs%mask3dCuL(G%IsdB:G%IedB,G%jsd:G%jed,1:NkIce))
  allocate(diag_cs%mask3dCvL(G%isd:G%ied,G%JsdB:G%JedB,1:NkIce))
  do k=1,NkIce
    diag_cs%mask3dTL(:,:,k) = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBL(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCuL(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvL(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo
  allocate(diag_cs%mask3dTi(G%isd:G%ied,G%jsd:G%jed,1:NkIce+1))
  allocate(diag_cs%mask3dBi(G%IsdB:G%IedB,G%JsdB:G%JedB,1:NkIce+1))
  allocate(diag_cs%mask3dCui(G%IsdB:G%IedB,G%jsd:G%jed,1:NkIce+1))
  allocate(diag_cs%mask3dCvi(G%isd:G%ied,G%JsdB:G%JedB,1:NkIce+1))
  do k=1,NkIce+1
    diag_cs%mask3dTi(:,:,k) = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBi(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCui(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvi(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo

  allocate(diag_cs%mask3dTC0(G%isd:G%ied,G%jsd:G%jed,0:CatIce))
  allocate(diag_cs%mask3dBuC0(G%IsdB:G%IedB,G%JsdB:G%JedB,0:CatIce))
  allocate(diag_cs%mask3dCuC0(G%IsdB:G%IedB,G%jsd:G%jed,0:CatIce))
  allocate(diag_cs%mask3dCvC0(G%isd:G%ied,G%JsdB:G%JedB,0:CatIce))
  do k=0,CatIce
    diag_cs%mask3dTC0(:,:,k)  = diag_cs%mask2dT(:,:)
    diag_cs%mask3dBuC0(:,:,k) = diag_cs%mask2dBu(:,:)
    diag_cs%mask3dCuC0(:,:,k) = diag_cs%mask2dCu(:,:)
    diag_cs%mask3dCvC0(:,:,k) = diag_cs%mask2dCv(:,:)
  enddo
  diag_cs%mask3dTC => diag_cs%mask3dTC0(:,:,1:CatIce)
  diag_cs%mask3dBuC => diag_cs%mask3dBuC0(:,:,1:CatIce)
  diag_cs%mask3dCuC => diag_cs%mask3dCuC0(:,:,1:CatIce)
  diag_cs%mask3dCvC => diag_cs%mask3dCvC0(:,:,1:CatIce)

  diag_cs%missing_value = missing_value

end subroutine diag_masks_set

!> Prevent the registration of additional diagnostics, so that the creation of files can occur
subroutine SIS_diag_mediator_close_registration(diag_CS)
  type(SIS_diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output

  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -2
  endif

end subroutine SIS_diag_mediator_close_registration

!> Deallocate memory associated with the SIS diag mediator
subroutine SIS_diag_mediator_end(time, diag_CS, end_diag_manager)
  type(time_type),     intent(in)    :: time !< The current model time
  type(SIS_diag_ctrl), intent(inout) :: diag_CS !< Structure used to regulate diagnostic output
  logical,   optional, intent(in)    :: end_diag_manager !< If true, call diag_manager_end()

  ! Local variables
  type(diag_type), pointer :: diag, next_diag
  integer :: i

  if (diag_CS%available_diag_doc_unit > -1) then
    close(diag_CS%available_diag_doc_unit) ; diag_CS%available_diag_doc_unit = -3
  endif
  if (diag_CS%chksum_iounit > -1) then
    close(diag_CS%chksum_iounit) ; diag_CS%chksum_iounit = -3
  endif

  do i=1, diag_cs%next_free_diag_id - 1
    if (associated(diag_cs%diags(i)%next)) then
      next_diag => diag_cs%diags(i)%next
      do while (associated(next_diag))
        diag => next_diag
        next_diag => diag%next
        deallocate(diag)
      enddo
    endif
  enddo

  deallocate(diag_cs%diags)

  ! These points to arrays in the grid type, so they can not be deallocated here.
  if (associated(diag_cs%mask2dT))  diag_cs%mask2dT => NULL()
  if (associated(diag_cs%mask2dBu)) diag_cs%mask2dBu => NULL()
  if (associated(diag_cs%mask2dCu)) diag_cs%mask2dCu => NULL()
  if (associated(diag_cs%mask2dCv)) diag_cs%mask2dCv => NULL()
  if (associated(diag_cs%mask2dT_comp)) deallocate(diag_cs%mask2dT_comp)

  if (associated(diag_cs%mask3dTL))  deallocate(diag_cs%mask3dTL)
  if (associated(diag_cs%mask3dBL))  deallocate(diag_cs%mask3dBL)
  if (associated(diag_cs%mask3dCuL)) deallocate(diag_cs%mask3dCuL)
  if (associated(diag_cs%mask3dCvL)) deallocate(diag_cs%mask3dCvL)
  if (associated(diag_cs%mask3dTi))  deallocate(diag_cs%mask3dTi)
  if (associated(diag_cs%mask3dBi))  deallocate(diag_cs%mask3dBi)
  if (associated(diag_cs%mask3dCui)) deallocate(diag_cs%mask3dCui)
  if (associated(diag_cs%mask3dCvi)) deallocate(diag_cs%mask3dCvi)
  if (associated(diag_cs%mask3dTC0))  deallocate(diag_cs%mask3dTC0)
  if (associated(diag_cs%mask3dBuC0)) deallocate(diag_cs%mask3dBuC0)
  if (associated(diag_cs%mask3dCuC0)) deallocate(diag_cs%mask3dCuC0)
  if (associated(diag_cs%mask3dCvC0)) deallocate(diag_cs%mask3dCvC0)
  ! These were pointers into their ...C0 counterpart, so they are nullified instead.
  diag_cs%mask3dTC => NULL() ; diag_cs%mask3dBuC => NULL()
  diag_cs%mask3dCuC => NULL() ; diag_cs%mask3dCvC => NULL()

  if (present(end_diag_manager)) then
    if (end_diag_manager) call MOM_diag_manager_end(time)
  endif

end subroutine SIS_diag_mediator_end

!> Convert the first n elements (up to 3) of an integer array to an underscore delimited string.
function i2s(a, n_in)
  integer, dimension(:), intent(in) :: a    !< The array of integers to translate
  integer, optional    , intent(in) :: n_in !< The number of elements to translate, by default all
  character(len=15) :: i2s !< The returned string

  ! Local variables
  character(len=15) :: i2s_temp
  integer :: i,n

  n = size(a)
  if (present(n_in)) n = n_in

  i2s = ''
  do i=1,min(n,3)
    write (i2s_temp, '(I4.4)') a(i)
    i2s = trim(i2s) //'_'// trim(i2s_temp)
  enddo
  i2s = adjustl(i2s)
end function i2s

!> Returns a new diagnostic id, it may be necessary to expand the diagnostics array.
integer function get_new_diag_id(diag_cs)
  type(SIS_diag_ctrl), intent(inout) :: diag_cs !< Diagnostics control structure
  ! Local variables
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

    ! Initialize new part of the diag array.
    do i=diag_cs%next_free_diag_id, size(diag_cs%diags)
      call initialize_diag_type(diag_cs%diags(i))
    enddo
  endif

  get_new_diag_id = diag_cs%next_free_diag_id
  diag_cs%next_free_diag_id = diag_cs%next_free_diag_id + 1

end function get_new_diag_id

!> Initializes a diag_type (used after allocating new memory)
subroutine initialize_diag_type(diag)
  type(diag_type), intent(inout) :: diag !< diag_type to be initialized

  diag%in_use = .false.
  diag%fms_diag_id = -1
  diag%axes => null()
  diag%next => null()
  diag%conversion_factor = 0.

end subroutine initialize_diag_type

!> Make a new diagnostic. Either use memory which is in the array of 'primary'
!! diagnostics, or if that is in use, insert it to the list of secondary diags.
subroutine alloc_diag_with_id(diag_id, diag_cs, diag)
  integer,                 intent(in   ) :: diag_id !< id for the diagnostic
  type(SIS_diag_ctrl), target, intent(inout) :: diag_cs !< structure used to regulate diagnostic output
  type(diag_type),         pointer       :: diag    !< structure representing a diagnostic (inout)

  type(diag_type), pointer :: tmp => NULL()

  if (.not. diag_cs%diags(diag_id)%in_use) then
    diag => diag_cs%diags(diag_id)
  else
    allocate(diag)
    tmp => diag_cs%diags(diag_id)%next
    diag_cs%diags(diag_id)%next => diag
    diag%next => tmp
  endif
  diag%in_use = .true.

end subroutine alloc_diag_with_id

!> Log a diagnostic to the available diagnostics file.
subroutine log_available_diag(used, module_name, field_name, cell_methods_string, comment, &
                              diag_CS, long_name, units, standard_name, variants, dimensions)
  logical,          intent(in) :: used !< Whether this diagnostic was in the diag_table or not
  character(len=*), intent(in) :: module_name !< Name of the diagnostic module
  character(len=*), intent(in) :: field_name !< Name of this diagnostic field
  character(len=*), intent(in) :: cell_methods_string !< The spatial component of the CF cell_methods attribute
  character(len=*), intent(in) :: comment !< A comment to append after [Used|Unused]
  type(SIS_diag_ctrl), intent(in) :: diag_CS  !< The diagnotics control structure
  character(len=*), optional, intent(in) :: dimensions !< Descriptor of the horizontal and vertical dimensions
  character(len=*), optional, intent(in) :: long_name !< CF long name of diagnostic
  character(len=*), optional, intent(in) :: units !< Units for diagnostic
  character(len=*), optional, intent(in) :: standard_name !< CF standardized name of diagnostic
  character(len=*), optional, intent(in) :: variants !< Alternate modules and variable names for
                                                     !! this diagnostic and derived diagnostics
  ! Local variables
  character(len=240) :: mesg

  if (used) then
    mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Used]'
  else
    mesg = '"'//trim(module_name)//'", "'//trim(field_name)//'"  [Unused]'
  endif
  !### This form of output agrees with MOM6:
  ! if (used) then
  !   mesg = '"'//trim(field_name)//'"  [Used]'
  ! else
  !   mesg = '"'//trim(field_name)//'"  [Unused]'
  ! endif
  if (len(trim((comment)))>0) then
    write(diag_CS%available_diag_doc_unit, '(a,1x,"(",a,")")') trim(mesg),trim(comment)
  else
    write(diag_CS%available_diag_doc_unit, '(a)') trim(mesg)
  endif
  !### Thes should be uncommented later to align with MOM6.
! call describe_option("modules", module_name, diag_CS)
! if (present(dimensions)) then
!   if (len(trim(dimensions)) > 0) then
!     call describe_option("dimensions", dimensions, diag_CS)
!   endif
! endif
  if (present(long_name)) call describe_option("long_name", long_name, diag_CS)
  if (present(units)) call describe_option("units", units, diag_CS)
  if (present(standard_name)) &
    call describe_option("standard_name", standard_name, diag_CS)
! if (len(trim((cell_methods_string)))>0) &
!   call describe_option("cell_methods", trim(cell_methods_string), diag_CS)
  if (present(variants)) then ; if (len(trim(variants)) > 0) then
    call describe_option("variants", variants, diag_CS)
  endif ; endif
end subroutine log_available_diag

!> Log the diagnostic chksum to the chksum diag file
subroutine log_chksum_diag(docunit, description, chksum)
  integer,          intent(in) :: docunit     !< Handle of the log file
  character(len=*), intent(in) :: description !< Name of the diagnostic module
  integer,          intent(in) :: chksum      !< chksum of the diagnostic

  write(docunit, '(a,1x,i9.8)') description, chksum
  flush(docunit)

end subroutine log_chksum_diag

!> Fakes a register of a diagnostic to find out if an obsolete
!! parameter appears in the diag_table.
logical function found_in_diagtable(diag, varName)
  type(SIS_diag_ctrl), intent(in) :: diag     !< A structure used to control diagnostics.
  character(len=*),    intent(in) :: varName  !< The obsolete diagnostic name
  ! Local
  integer :: handle ! Integer handle returned from diag_manager

  ! We use register_static_field_fms() instead of register_static_field() so
  ! that the diagnostic does not appear in the available diagnostics list.
  handle = register_static_field_infra('ice_model', varName, diag%axesT1%handles)

  found_in_diagtable = (handle>0)

end function found_in_diagtable

!> Finishes the diag manager reduction methods as needed for the time_step
subroutine SIS_diag_send_complete()
  call diag_send_complete_infra()
end subroutine SIS_diag_send_complete

end module SIS_diag_mediator
