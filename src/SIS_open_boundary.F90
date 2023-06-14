!> Controls where open boundary conditions are applied
module SIS_open_boundary

! This file is part of SIS2. See LICENSE.md for the license.

use MOM_domains,              only : pass_var, pass_vector
use MOM_domains,              only : To_All, EAST_FACE, NORTH_FACE, SCALAR_PAIR, CGRID_NE, CORNER
use MOM_error_handler,        only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,          only : get_param, log_version, param_file_type, log_param
use MOM_grid,                 only : ocean_grid_type, hor_index_type
use MOM_dyn_horgrid,          only : dyn_horgrid_type
use MOM_interpolate,          only : init_external_field, time_interp_external, time_interp_external_init
use MOM_open_boundary,        only : OBC_NONE
use MOM_open_boundary,        only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary,        only : parse_segment_str, flood_fill, flood_fill2
use MOM_string_functions,     only : remove_spaces
use MOM_time_manager,         only : set_date, time_type, time_type_to_real, operator(-)
use MOM_unit_scaling,         only : unit_scale_type

implicit none ; private

#include <SIS2_memory.h>

public open_boundary_config
!public open_boundary_init
public open_boundary_impose_land_mask

integer, parameter         :: MAX_OBC_FIELDS = 100  !< Maximum number of data fields needed for OBC segments

!> Open boundary segment data from files (mostly).
type, public :: OBC_segment_data_type
  integer :: fid                            !< handle from FMS associated with segment data on disk
  integer :: fid_dz                         !< handle from FMS associated with segment thicknesses on disk
  character(len=8)                :: name   !< a name identifier for the segment data
  real, allocatable :: buffer_src(:,:,:)    !< buffer for segment data located at cell faces
  integer                         :: nk_src !< Number of vertical levels in the source data
  real, allocatable :: dz_src(:,:,:)        !< vertical grid cell spacing of the incoming segment
                                            !! data, set in [Z ~> m]
  real, allocatable :: buffer_dst(:,:,:)    !< buffer src data remapped to the target vertical grid
  real              :: value                !< constant value if fid is equal to -1
end type OBC_segment_data_type


!> Tracer on OBC segment data structure, for putting into a segment tracer registry.
type, public :: OBC_segment_tracer_type
  real, allocatable          :: t(:,:,:)              !< tracer concentration array
  real                       :: OBC_inflow_conc = 0.0 !< tracer concentration for generic inflows
  character(len=32)          :: name                  !< tracer name used for error messages
! type(tracer_type), pointer :: Tr => NULL()          !< metadata describing the tracer
  real, allocatable          :: tres(:,:,:)           !< tracer reservoir array
  logical                    :: is_initialized        !< reservoir values have been set when True
end type OBC_segment_tracer_type

!> Registry type for tracers on segments
type, public :: segment_tracer_registry_type
  integer                       :: ntseg = 0         !< number of registered tracer segments
  type(OBC_segment_tracer_type) :: Tr(MAX_FIELDS_)   !< array of registered tracers
  logical                       :: locked = .false.  !< New tracers may be registered if locked=.false.
                                                     !! When locked=.true.,no more tracers can be registered.
                                                     !! Not sure who should lock it or when...
end type segment_tracer_registry_type

!> Open boundary segment data structure.
type, public :: OBC_segment_type
  logical :: Flather        !< If true, applies Flather + Chapman radiation of barotropic gravity waves.
  logical :: radiation      !< If true, 1D Orlanksi radiation boundary conditions are applied.
                            !! If False, a gradient condition is applied.
  logical :: radiation_tan  !< If true, 1D Orlanksi radiation boundary conditions are applied to
                            !! tangential flows.
  logical :: radiation_grad !< If true, 1D Orlanksi radiation boundary conditions are applied to
                            !! dudv and dvdx.
  logical :: nudged         !< Optional supplement to radiation boundary.
  logical :: nudged_tan     !< Optional supplement to nudge tangential velocity.
  logical :: nudged_grad    !< Optional supplement to nudge normal gradient of tangential velocity.
  logical :: specified      !< Boundary normal velocity fixed to external value.
  logical :: specified_tan  !< Boundary tangential velocity fixed to external value.
  logical :: specified_grad !< Boundary gradient of tangential velocity fixed to external value.
  logical :: open           !< Boundary is open for continuity solver.
  logical :: gradient       !< Zero gradient at boundary.
  logical :: values_needed  !< Whether or not any external OBC fields are needed.
  logical :: u_values_needed      !< Whether or not external u OBC fields are needed.
  logical :: v_values_needed      !< Whether or not external v OBC fields are needed.
  logical :: vamp_values_needed   !< Whether or not external v amplitude OBC fields are needed.
  logical :: g_values_needed!< Whether or not external gradient OBC fields are needed.
  integer :: direction      !< Boundary faces one of the four directions.
  logical :: is_N_or_S      !< True if the OB is facing North or South and exists on this PE.
  logical :: is_E_or_W      !< True if the OB is facing East or West and exists on this PE.
  logical :: is_E_or_W_2    !< True if the OB is facing East or West anywhere.
  type(OBC_segment_data_type), pointer :: field(:) => NULL()  !< OBC data
  integer :: num_fields     !< number of OBC data fields (e.g. u_normal,u_parallel and eta for Flather)
  integer :: Is_obc         !< i-indices of boundary segment.
  integer :: Ie_obc         !< i-indices of boundary segment.
  integer :: Js_obc         !< j-indices of boundary segment.
  integer :: Je_obc         !< j-indices of boundary segment.
  real :: Velocity_nudging_timescale_in  !< Nudging timescale on inflow [T ~> s].
  real :: Velocity_nudging_timescale_out !< Nudging timescale on outflow [T ~> s].
  logical :: on_pe          !< true if any portion of the segment is located in this PE's data domain
  real, allocatable :: Cg(:,:)  !< The external gravity wave speed [L T-1 ~> m s-1]
                                !! at OBC-points.
  real, allocatable :: normal_vel(:,:,:)      !< The layer velocity normal to the OB
                                              !! segment [L T-1 ~> m s-1].
  real, allocatable :: tangential_vel(:,:,:)  !< The layer velocity tangential to the
                                              !! OB segment [L T-1 ~> m s-1].
  real, allocatable :: normal_trans(:,:,:)    !< The layer transport normal to the OB
                                              !! segment [H L2 T-1 ~> m3 s-1].
  type(segment_tracer_registry_type), pointer  :: tr_Reg=> NULL()!< A pointer to the tracer registry for the segment.
  type(hor_index_type) :: HI !< Horizontal index ranges
  real :: Tr_InvLscale_out                                  !< An effective inverse length scale for restoring
                                                            !! the tracer concentration in a fictitious
                                                            !! reservoir towards interior values when flow
                                                            !! is exiting the domain [L-1 ~> m-1]
  real :: Tr_InvLscale_in                                   !< An effective inverse length scale for restoring
                                                            !! the tracer concentration towards an externally
                                                            !! imposed value when flow is entering [L-1 ~> m-1]
end type OBC_segment_type

!> Open-boundary data
type, public :: ice_OBC_type
  integer :: number_of_segments = 0                   !< The number of open-boundary segments.
  integer :: ke = 0                                   !< The number of model layers
  logical :: open_u_BCs_exist_globally = .false.      !< True if any zonal velocity points
                                                      !! in the global domain use open BCs.
  logical :: open_v_BCs_exist_globally = .false.      !< True if any meridional velocity points
                                                      !! in the global domain use open BCs.
  logical :: Flather_u_BCs_exist_globally = .false.   !< True if any zonal velocity points
                                                      !! in the global domain use Flather BCs.
  logical :: Flather_v_BCs_exist_globally = .false.   !< True if any meridional velocity points
                                                      !! in the global domain use Flather BCs.
  logical :: nudged_u_BCs_exist_globally = .false.    !< True if any velocity points in the
                                                      !! global domain use nudged BCs.
  logical :: nudged_v_BCs_exist_globally = .false.    !< True if any velocity points in the
                                                      !! global domain use nudged BCs.
  logical :: specified_u_BCs_exist_globally = .false. !< True if any zonal velocity points
                                                      !! in the global domain use specified BCs.
  logical :: specified_v_BCs_exist_globally = .false. !< True if any meridional velocity points
                                                      !! in the global domain use specified BCs.
  logical :: radiation_BCs_exist_globally = .false.   !< True if radiations BCs are in use anywhere.
  logical :: user_BCs_set_globally = .false.          !< True if any OBC_USER_CONFIG is set
                                                      !! for input from user directory.
  logical :: update_OBC = .false.                     !< Is OBC data time-dependent
  logical :: needs_IO_for_data = .false.              !< Is any i/o needed for OBCs
  logical :: zero_vorticity = .false.                 !< If True, sets relative vorticity to zero on open boundaries.
  logical :: freeslip_vorticity = .false.             !< If True, sets normal gradient of tangential velocity to zero
                                                      !! in the relative vorticity on open boundaries.
  logical :: computed_vorticity = .false.             !< If True, uses external data for tangential velocity
                                                      !! in the relative vorticity on open boundaries.
  logical :: specified_vorticity = .false.            !< If True, uses external data for tangential velocity
                                                      !! gradients in the relative vorticity on open boundaries.
  logical :: zero_strain = .false.                    !< If True, sets strain to zero on open boundaries.
  logical :: freeslip_strain = .false.                !< If True, sets normal gradient of tangential velocity to zero
                                                      !! in the strain on open boundaries.
  logical :: computed_strain = .false.                !< If True, uses external data for tangential velocity to compute
                                                      !! normal gradient in the strain on open boundaries.
  logical :: specified_strain = .false.               !< If True, uses external data for tangential velocity gradients
                                                      !! to compute strain on open boundaries.
  logical :: zero_biharmonic = .false.                !< If True, zeros the Laplacian of flow on open boundaries for
                                                      !! use in the biharmonic viscosity term.
  logical :: brushcutter_mode = .false.               !< If True, read data on supergrid.
  logical, allocatable :: tracer_x_reservoirs_used(:) !< Dimensioned by the number of tracers, set globally,
                                                      !! true for those with x reservoirs (needed for restarts).
  logical, allocatable :: tracer_y_reservoirs_used(:) !< Dimensioned by the number of tracers, set globally,
                                                      !! true for those with y reservoirs (needed for restarts).
  integer                       :: ntr = 0            !< number of tracers

  ! Properties of the segments used.
  type(OBC_segment_type), allocatable :: segment(:)   !< List of segment objects.

  ! Which segment object describes the current point.
  integer, allocatable :: segnum_u(:,:) !< Segment number of u-points.
  integer, allocatable :: segnum_v(:,:) !< Segment number of v-points.
  logical :: OBC_pe                     !< Is there an open boundary on this tile?
  type(OBC_registry_type), pointer :: OBC_Reg => NULL()  !< Registry type for boundaries
  real, allocatable :: tres_x(:,:,:,:)   !< Array storage of tracer reservoirs for restarts [conc L ~> conc m]
  real, allocatable :: tres_y(:,:,:,:)   !< Array storage of tracer reservoirs for restarts [conc L ~> conc m]
  real :: silly_h  !< A silly value of thickness outside of the domain that can be used to test
                   !! the independence of the OBCs to this external data [H ~> m or kg m-2].
  real :: silly_u  !< A silly value of velocity outside of the domain that can be used to test
                   !! the independence of the OBCs to this external data [L T-1 ~> m s-1].
end type ice_OBC_type

!> Control structure for open boundaries that read from files.
!! Probably lots to update here.
type, public :: file_OBC_CS ; private
  real :: tide_flow = 3.0e6         !< Placeholder for now...
end type file_OBC_CS

!> Type to carry something (what??) for the OBC registry.
type, public :: OBC_struct_type
  character(len=32)               :: name             !< OBC name used for error messages
end type OBC_struct_type

!> Type to carry basic OBC information needed for updating values.
type, public :: OBC_registry_type
  integer               :: nobc = 0          !< number of registered open boundary types.
  type(OBC_struct_type) :: OB(MAX_FIELDS_)   !< array of registered boundary types.
  logical               :: locked = .false.  !< New OBC types may be registered if locked=.false.
                                             !! When locked=.true.,no more boundaries can be registered.
end type OBC_registry_type

integer :: id_clock_pass !< A CPU time clock

character(len=40)  :: mdl = "SIS_open_boundary" !< This module's name.

contains

!> Enables OBC module and reads configuration parameters
!! This routine is called from SIS_initialize_fixed which
!! occurs before the initialization of the vertical coordinate.
!! Therefore segment data are not fully initialized
!! here. The remainder of the segment data are initialized in a
!! later call to update_open_boundary_data
subroutine open_boundary_config(G, US, param_file, OBC)
  type(dyn_horgrid_type),  intent(inout) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< Parameter file handle
  type(ice_OBC_type),      pointer       :: OBC !< Open boundary control structure

  ! Local variables
  integer :: l ! For looping over segments
  logical :: debug_OBC, debug, mask_outside, reentrant_x, reentrant_y
  character(len=15) :: segment_param_str ! The run-time parameter name for each segment
  character(len=1024) :: segment_str      ! The contents (rhs) for parameter "segment_param_str"
  character(len=200) :: config1          ! String for OBC_USER_CONFIG
  real               :: Lscale_in, Lscale_out ! parameters controlling tracer values at the boundaries [L ~> m]
  character(len=128) :: inputdir
  logical :: check_reconstruction, check_remapping, force_bounds_in_subcell
  character(len=64)  :: remappingScheme
  logical :: default_2018_answers
  ! This include declares and sets the variable "version".
# include "version_variable.h"

  allocate(OBC)

  call get_param(param_file, mdl, "OBC_NUMBER_OF_SEGMENTS", OBC%number_of_segments, &
                 default=0, do_not_log=.true.)
  call log_version(param_file, mdl, version, &
                 "Controls where open boundaries are located, what kind of boundary condition "// &
                 "to impose, and what data to apply, if any.", &
                 all_default=(OBC%number_of_segments<=0))
  call get_param(param_file, mdl, "OBC_NUMBER_OF_SEGMENTS", OBC%number_of_segments, &
                 "The number of open boundary segments.", &
                 default=0)
  call get_param(param_file, mdl, "OBC_USER_CONFIG", config1, &
                 "A string that sets how the open boundary conditions are "// &
                 " configured: \n", default="none", do_not_log=.true.)
  call get_param(param_file, mdl, "NK", OBC%ke, &
                 "The number of model layers", default=0, do_not_log=.true.)

  if (config1 /= "none" .and. config1 /= "dyed_obcs") OBC%user_BCs_set_globally = .true.

  if (OBC%number_of_segments > 0) then
    call get_param(param_file, mdl, "OBC_ZERO_VORTICITY", OBC%zero_vorticity, &
         "If true, sets relative vorticity to zero on open boundaries.", &
         default=.false.)
    call get_param(param_file, mdl, "OBC_FREESLIP_VORTICITY", OBC%freeslip_vorticity, &
         "If true, sets the normal gradient of tangential velocity to "// &
         "zero in the relative vorticity on open boundaries. This cannot "// &
         "be true if another OBC_XXX_VORTICITY option is True.", default=.true.)
    call get_param(param_file, mdl, "OBC_COMPUTED_VORTICITY", OBC%computed_vorticity, &
         "If true, uses the external values of tangential velocity "// &
         "in the relative vorticity on open boundaries. This cannot "// &
         "be true if another OBC_XXX_VORTICITY option is True.", default=.false.)
    call get_param(param_file, mdl, "OBC_SPECIFIED_VORTICITY", OBC%specified_vorticity, &
         "If true, uses the external values of tangential velocity "// &
         "in the relative vorticity on open boundaries. This cannot "// &
         "be true if another OBC_XXX_VORTICITY option is True.", default=.false.)
    if ((OBC%zero_vorticity .and. OBC%freeslip_vorticity) .or.  &
        (OBC%zero_vorticity .and. OBC%computed_vorticity) .or.  &
        (OBC%zero_vorticity .and. OBC%specified_vorticity) .or.  &
        (OBC%freeslip_vorticity .and. OBC%computed_vorticity) .or.  &
        (OBC%freeslip_vorticity .and. OBC%specified_vorticity) .or.  &
        (OBC%computed_vorticity .and. OBC%specified_vorticity))  &
         call MOM_error(FATAL, "SIS_open_boundary.F90, open_boundary_config:\n"// &
         "Only one of OBC_ZERO_VORTICITY, OBC_FREESLIP_VORTICITY, OBC_COMPUTED_VORTICITY\n"// &
         "and OBC_IMPORTED_VORTICITY can be True at once.")
    call get_param(param_file, mdl, "OBC_ZERO_STRAIN", OBC%zero_strain, &
         "If true, sets the strain used in the stress tensor to zero on open boundaries.", &
         default=.false.)
    call get_param(param_file, mdl, "OBC_FREESLIP_STRAIN", OBC%freeslip_strain, &
         "If true, sets the normal gradient of tangential velocity to "// &
         "zero in the strain use in the stress tensor on open boundaries. This cannot "//&
         "be true if another OBC_XXX_STRAIN option is True.", default=.true.)
    call get_param(param_file, mdl, "OBC_COMPUTED_STRAIN", OBC%computed_strain, &
         "If true, sets the normal gradient of tangential velocity to "// &
         "zero in the strain use in the stress tensor on open boundaries. This cannot "// &
         "be true if another OBC_XXX_STRAIN option is True.", default=.false.)
    call get_param(param_file, mdl, "OBC_SPECIFIED_STRAIN", OBC%specified_strain, &
         "If true, sets the normal gradient of tangential velocity to "// &
         "zero in the strain use in the stress tensor on open boundaries. This cannot "// &
         "be true if another OBC_XXX_STRAIN option is True.", default=.false.)
    if ((OBC%zero_strain .and. OBC%freeslip_strain) .or.  &
        (OBC%zero_strain .and. OBC%computed_strain) .or.  &
        (OBC%zero_strain .and. OBC%specified_strain) .or.  &
        (OBC%freeslip_strain .and. OBC%computed_strain) .or.  &
        (OBC%freeslip_strain .and. OBC%specified_strain) .or.  &
        (OBC%computed_strain .and. OBC%specified_strain))  &
         call MOM_error(FATAL, "SIS_open_boundary.F90, open_boundary_config: \n"//&
         "Only one of OBC_ZERO_STRAIN, OBC_FREESLIP_STRAIN, OBC_COMPUTED_STRAIN \n"//&
         "and OBC_IMPORTED_STRAIN can be True at once.")
    call get_param(param_file, mdl, "OBC_ZERO_BIHARMONIC", OBC%zero_biharmonic, &
         "If true, zeros the Laplacian of flow on open boundaries in the biharmonic "// &
         "viscosity term.", default=.false.)
    call get_param(param_file, mdl, "MASK_OUTSIDE_OBCS", mask_outside, &
         "If true, set the areas outside open boundaries to be land.", &
         default=.false.)

    call get_param(param_file, mdl, "DEBUG", debug, default=.false.)
    call get_param(param_file, mdl, "DEBUG_OBC", debug_OBC, default=.false.)
    if (debug_OBC .or. debug) &
      call log_param(param_file, mdl, "DEBUG_OBC", debug_OBC, &
                 "If true, do additional calls to help debug the performance "//&
                 "of the open boundary condition code.", default=.false., &
                 debuggingParam=.true.)
    call get_param(param_file, mdl, "OBC_SILLY_THICK", OBC%silly_h, &
                 "A silly value of thicknesses used outside of open boundary "//&
                 "conditions for debugging.", units="m", default=0.0, scale=US%m_to_Z, &
                 do_not_log=.not.debug_OBC, debuggingParam=.true.)
    call get_param(param_file, mdl, "OBC_SILLY_VEL", OBC%silly_u, &
                 "A silly value of velocities used outside of open boundary "//&
                 "conditions for debugging.", units="m/s", default=0.0, scale=US%m_s_to_L_T, &
                 do_not_log=.not.debug_OBC, debuggingParam=.true.)
    reentrant_x = .false.
    call get_param(param_file, mdl, "REENTRANT_X", reentrant_x, default=.true.)
    reentrant_y = .false.
    call get_param(param_file, mdl, "REENTRANT_Y", reentrant_y, default=.false.)

    ! Allocate everything
    allocate(OBC%segment(1:OBC%number_of_segments))
    do l=1,OBC%number_of_segments
      OBC%segment(l)%Flather = .false.
      OBC%segment(l)%radiation = .false.
      OBC%segment(l)%radiation_tan = .false.
      OBC%segment(l)%radiation_grad = .false.
      OBC%segment(l)%nudged = .false.
      OBC%segment(l)%nudged_tan = .false.
      OBC%segment(l)%nudged_grad = .false.
      OBC%segment(l)%specified = .false.
      OBC%segment(l)%specified_tan = .false.
      OBC%segment(l)%specified_grad = .false.
      OBC%segment(l)%open = .false.
      OBC%segment(l)%gradient = .false.
      OBC%segment(l)%values_needed = .false.
      OBC%segment(l)%u_values_needed = .false.
      OBC%segment(l)%v_values_needed = .false.
      OBC%segment(l)%direction = OBC_NONE
      OBC%segment(l)%is_N_or_S = .false.
      OBC%segment(l)%is_E_or_W = .false.
      OBC%segment(l)%is_E_or_W_2 = .false.
      OBC%segment(l)%Velocity_nudging_timescale_in = 0.0
      OBC%segment(l)%Velocity_nudging_timescale_out = 0.0
      OBC%segment(l)%num_fields = 0
    enddo
    allocate(OBC%segnum_u(G%IsdB:G%IedB,G%jsd:G%jed), source=OBC_NONE)
    allocate(OBC%segnum_v(G%isd:G%ied,G%JsdB:G%JedB), source=OBC_NONE)

    do l = 1, OBC%number_of_segments
      write(segment_param_str(1:15),"('OBC_SEGMENT_',i3.3)") l
      call get_param(param_file, mdl, segment_param_str, segment_str, &
           "Documentation needs to be dynamic?????", &
           fail_if_missing=.true.)
      segment_str = remove_spaces(segment_str)
      if (segment_str(1:2) == 'I=') then
        call setup_u_point_obc(OBC, G, US, segment_str, l, param_file, reentrant_y)
      elseif (segment_str(1:2) == 'J=') then
        call setup_v_point_obc(OBC, G, US, segment_str, l, param_file, reentrant_x)
      else
        call MOM_error(FATAL, "SIS_open_boundary.F90, open_boundary_config: "//&
             "Unable to interpret "//segment_param_str//" = "//trim(segment_str))
      endif
    enddo

    ! Moved this earlier because time_interp_external_init needs to be called
    ! before anything that uses time_interp_external (such as initialize_segment_data)
    if (OBC%specified_u_BCs_exist_globally .or.  OBC%specified_v_BCs_exist_globally .or. &
      OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) then
      ! Need this for boundary interpolation.
      call time_interp_external_init()
    endif

    Lscale_in = 0.
    Lscale_out = 0.
    if (open_boundary_query(OBC, apply_open_OBC=.true.)) then
      call get_param(param_file, mdl, "OBC_TRACER_RESERVOIR_LENGTH_SCALE_OUT ", Lscale_out, &
                 "An effective length scale for restoring the tracer concentration "//&
                 "at the boundaries to externally imposed values when the flow "//&
                 "is exiting the domain.", units="m", default=0.0, scale=US%m_to_L)

      call get_param(param_file, mdl, "OBC_TRACER_RESERVOIR_LENGTH_SCALE_IN ", Lscale_in, &
                 "An effective length scale for restoring the tracer concentration "//&
                 "at the boundaries to values from the interior when the flow "//&
                 "is entering the domain.", units="m", default=0.0, scale=US%m_to_L)
    endif

    if (mask_outside) call mask_outside_OBCs(G, US, param_file, OBC)

    ! All tracers are using the same restoring length scale for now, but we may want to make this
    ! tracer-specific in the future for example, in cases where certain tracers are poorly constrained
    ! by data while others are well constrained - MJH.
    do l = 1, OBC%number_of_segments
      OBC%segment(l)%Tr_InvLscale_in = 0.0
      if (Lscale_in>0.) OBC%segment(l)%Tr_InvLscale_in =  1.0/Lscale_in
      OBC%segment(l)%Tr_InvLscale_out = 0.0
      if (Lscale_out>0.) OBC%segment(l)%Tr_InvLscale_out =  1.0/Lscale_out
    enddo

    call get_param(param_file, mdl, "BRUSHCUTTER_MODE", OBC%brushcutter_mode, &
         "If true, read external OBC data on the supergrid.", &
         default=.false.)

  endif ! OBC%number_of_segments > 0

  ! Safety check
  if ((OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) .and. &
       .not.G%symmetric ) call MOM_error(FATAL, &
       "SIS_open_boundary, open_boundary_config: "//&
       "Symmetric memory must be used when using Flather OBCs.")

  if (.not.(OBC%specified_u_BCs_exist_globally .or. OBC%specified_v_BCs_exist_globally .or. &
            OBC%open_u_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally)) then
    ! No open boundaries have been requested
    call open_boundary_dealloc(OBC)
  endif

end subroutine open_boundary_config


!> helper function for finding out about OBCs
logical function open_boundary_query(OBC, apply_open_OBC, apply_specified_OBC, apply_Flather_OBC, &
                                     apply_nudged_OBC, needs_ext_seg_data)
  type(ice_OBC_type),   pointer    :: OBC !< Open boundary control structure
  logical, optional,    intent(in) :: apply_open_OBC      !< Returns True if open_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_specified_OBC !< Returns True if specified_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_Flather_OBC   !< Returns True if Flather_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: apply_nudged_OBC    !< Returns True if nudged_*_BCs_exist_globally is true
  logical, optional,    intent(in) :: needs_ext_seg_data  !< Returns True if external segment data needed
  open_boundary_query = .false.
  if (.not. associated(OBC)) return
  if (present(apply_open_OBC)) open_boundary_query = OBC%open_u_BCs_exist_globally .or. &
                                                     OBC%open_v_BCs_exist_globally
  if (present(apply_specified_OBC)) open_boundary_query = OBC%specified_u_BCs_exist_globally .or. &
                                                          OBC%specified_v_BCs_exist_globally
  if (present(apply_Flather_OBC)) open_boundary_query = OBC%Flather_u_BCs_exist_globally .or. &
                                                        OBC%Flather_v_BCs_exist_globally
  if (present(apply_nudged_OBC)) open_boundary_query = OBC%nudged_u_BCs_exist_globally .or. &
                                                       OBC%nudged_v_BCs_exist_globally
  if (present(needs_ext_seg_data)) open_boundary_query = OBC%needs_IO_for_data

end function open_boundary_query

!> Define indices for segment and store in hor_index_type
!! using global segment bounds corresponding to q-points.
!! Copied from MOM6 for this segment type.
subroutine setup_segment_indices(G, seg, Is_obc, Ie_obc, Js_obc, Je_obc)
  type(dyn_horgrid_type), intent(in) :: G !< grid type
  type(OBC_segment_type), intent(inout) :: seg  !< Open boundary segment
  integer, intent(in) :: Is_obc !< Q-point global i-index of start of segment
  integer, intent(in) :: Ie_obc !< Q-point global i-index of end of segment
  integer, intent(in) :: Js_obc !< Q-point global j-index of start of segment
  integer, intent(in) :: Je_obc !< Q-point global j-index of end of segment
  ! Local variables
  integer :: IsgB, IegB, JsgB, JegB
  integer :: isg, ieg, jsg, jeg

  ! Isg, Ieg will be I*_obc in global space
  if (Ie_obc < Is_obc) then
    IsgB = Ie_obc
    IegB = Is_obc
  else
    IsgB = Is_obc
    IegB = Ie_obc
  endif

  if (Je_obc < Js_obc) then
    JsgB = Je_obc
    JegB = Js_obc
  else
    JsgB = Js_obc
    JegB = Je_obc
  endif

  ! NOTE: h-points are defined along the interior of the segment q-points.
  !   For a given segment and its start and end index pairs, [IJ][se]gB, the
  !   h-cell corresponding to this pair are shown in the figure below.
  !
  ! x-x----------------x-x
  ! | |        N       | |
  ! x-x   W         E  x-x
  !   |        S         |
  ! x-x----------------x-x
  ! | |                | |
  ! x-x                x-x
  !
  ! For segment points on the west and south, h-point indices are incremented
  ! in order to move to the interior cell.

  if (Is_obc > Ie_obc) then
    ! Northern boundary
    isg = IsgB + 1
    jsg = JsgB
    ieg = IegB
    jeg = JegB
  endif

  if (Is_obc < Ie_obc) then
    ! Southern boundary
    isg = IsgB + 1
    jsg = JsgB + 1
    ieg = IegB
    jeg = JegB + 1
  endif

  if (Js_obc < Je_obc) then
    ! Eastern boundary
    isg = IsgB
    jsg = JsgB + 1
    ieg = IegB
    jeg = JegB
  endif

  if (Js_obc > Je_obc) then
    ! Western boundary
    isg = IsgB + 1
    jsg = JsgB + 1
    ieg = IegB + 1
    jeg = JegB
  endif

  ! Global space I*_obc but sorted
  seg%HI%IsgB = IsgB
  seg%HI%JegB = JegB
  seg%HI%IegB = IegB
  seg%HI%JsgB = JsgB

  seg%HI%isg = isg
  seg%HI%jsg = jsg
  seg%HI%ieg = ieg
  seg%HI%jeg = jeg

  ! Move into local index space
  IsgB = IsgB - G%idg_offset
  JsgB = JsgB - G%jdg_offset
  IegB = IegB - G%idg_offset
  JegB = JegB - G%jdg_offset

  isg = isg - G%idg_offset
  jsg = jsg - G%jdg_offset
  ieg = ieg - G%idg_offset
  jeg = jeg - G%jdg_offset

  ! This is the i-extent of the segment on this PE.
  ! The values are nonsense if the segment is not on this PE.
  seg%HI%IsdB = min(max(IsgB, G%HI%IsdB), G%HI%IedB)
  seg%HI%IedB = min(max(IegB, G%HI%IsdB), G%HI%IedB)
  seg%HI%isd = min(max(isg, G%HI%isd), G%HI%ied)
  seg%HI%ied = min(max(ieg, G%HI%isd), G%HI%ied)
  seg%HI%IscB = min(max(IsgB, G%HI%IscB), G%HI%IecB)
  seg%HI%IecB = min(max(IegB, G%HI%IscB), G%HI%IecB)
  seg%HI%isc = min(max(isg, G%HI%isc), G%HI%iec)
  seg%HI%iec = min(max(ieg, G%HI%isc), G%HI%iec)

  ! This is the j-extent of the segment on this PE.
  ! The values are nonsense if the segment is not on this PE.
  seg%HI%JsdB = min(max(JsgB, G%HI%JsdB), G%HI%JedB)
  seg%HI%JedB = min(max(JegB, G%HI%JsdB), G%HI%JedB)
  seg%HI%jsd = min(max(jsg, G%HI%jsd), G%HI%jed)
  seg%HI%jed = min(max(jeg, G%HI%jsd), G%HI%jed)
  seg%HI%JscB = min(max(JsgB, G%HI%JscB), G%HI%JecB)
  seg%HI%JecB = min(max(JegB, G%HI%JscB), G%HI%JecB)
  seg%HI%jsc = min(max(jsg, G%HI%jsc), G%HI%jec)
  seg%HI%jec = min(max(jeg, G%HI%jsc), G%HI%jec)

end subroutine setup_segment_indices

!> Parse an OBC_SEGMENT_%%% string starting with "I=" and configure placement
!! and type of OBC accordingly
subroutine setup_u_point_obc(OBC, G, US, segment_str, l_seg, PF, reentrant_y)
  type(ice_OBC_type),   intent(inout) :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),  intent(in) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  character(len=*),        intent(in) :: segment_str !< A string in form of "I=%,J=%:%,string"
  integer,                 intent(in) :: l_seg !< which segment is this?
  type(param_file_type), intent(in)   :: PF  !< Parameter file handle
  logical, intent(in)                 :: reentrant_y !< is the domain reentrant in y?
  ! Local variables
  integer :: I_obc, Js_obc, Je_obc ! Position of segment in global index space
  integer :: j, a_loop
  character(len=32) :: action_str(8)
  character(len=128) :: segment_param_str
  real, allocatable, dimension(:)  :: tnudge ! Nudging timescales [T ~> s]
  ! This returns the global indices for the segment
  call parse_segment_str(G%ieg, G%jeg, segment_str, I_obc, Js_obc, Je_obc, action_str, reentrant_y)

  call setup_segment_indices(G, OBC%segment(l_seg),I_obc,I_obc,Js_obc,Je_obc)

  I_obc = I_obc - G%idg_offset ! Convert to local tile indices on this tile
  Js_obc = Js_obc - G%jdg_offset ! Convert to local tile indices on this tile
  Je_obc = Je_obc - G%jdg_offset ! Convert to local tile indices on this tile

  if (Je_obc>Js_obc) then
    OBC%segment(l_seg)%direction = OBC_DIRECTION_E
  elseif (Je_obc<Js_obc) then
    OBC%segment(l_seg)%direction = OBC_DIRECTION_W
    j=js_obc;js_obc=je_obc;je_obc=j
  endif

  OBC%segment(l_seg)%on_pe = .false.

  do a_loop = 1,8 ! up to 8 options available
    if (len_trim(action_str(a_loop)) == 0) then
      cycle
    elseif (trim(action_str(a_loop)) == 'FLATHER') then
      OBC%segment(l_seg)%Flather = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%Flather_u_BCs_exist_globally = .true.
      OBC%open_u_BCs_exist_globally = .true.
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_u_BCs_exist_globally = .true.
      OBC%radiation_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI_TAN') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%radiation_tan = .true.
      OBC%radiation_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI_GRAD') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%radiation_grad = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED') then
      OBC%segment(l_seg)%nudged = .true.
      OBC%nudged_u_BCs_exist_globally = .true.
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_TAN') then
      OBC%segment(l_seg)%nudged_tan = .true.
      OBC%nudged_u_BCs_exist_globally = .true.
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_GRAD') then
      OBC%segment(l_seg)%nudged_grad = .true.
      OBC%segment(l_seg)%g_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_u_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_u_BCs_exist_globally = .true. ! This avoids deallocation
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_TAN') then
      OBC%segment(l_seg)%specified_tan = .true.
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_GRAD') then
      OBC%segment(l_seg)%specified_grad = .true.
      OBC%segment(l_seg)%g_values_needed = .true.
    else
      call MOM_error(FATAL, "SIS_open_boundary.F90, setup_u_point_obc: "//&
                     "String '"//trim(action_str(a_loop))//"' not understood.")
    endif
    if (OBC%segment(l_seg)%nudged .or. OBC%segment(l_seg)%nudged_tan) then
      write(segment_param_str(1:43),"('OBC_SEGMENT_',i3.3,'_VELOCITY_NUDGING_TIMESCALES')") l_seg
      allocate(tnudge(2))
      call get_param(PF, mdl, segment_param_str(1:43), tnudge, &
                     "Timescales in days for nudging along a segment, "//&
                     "for inflow, then outflow. Setting both to zero should "//&
                     "behave like SIMPLE obcs for the baroclinic velocities.", &
                     fail_if_missing=.true., default=0., units="days", scale=86400.0*US%s_to_T)
      OBC%segment(l_seg)%Velocity_nudging_timescale_in = tnudge(1)
      OBC%segment(l_seg)%Velocity_nudging_timescale_out = tnudge(2)
      deallocate(tnudge)
    endif

  enddo ! a_loop

  OBC%segment(l_seg)%is_E_or_W_2 = .true.

  if (I_obc<=G%HI%IsdB+1 .or. I_obc>=G%HI%IedB-1) return ! Boundary is not on tile
  if (Je_obc<=G%HI%JsdB .or. Js_obc>=G%HI%JedB) return ! Segment is not on tile

  OBC%segment(l_seg)%on_pe = .true.
  OBC%segment(l_seg)%is_E_or_W = .true.

  do j=G%HI%jsd, G%HI%jed
    if (j>Js_obc .and. j<=Je_obc) then
      OBC%segnum_u(I_obc,j) = l_seg
    endif
  enddo
  OBC%segment(l_seg)%Is_obc = I_obc
  OBC%segment(l_seg)%Ie_obc = I_obc
  OBC%segment(l_seg)%Js_obc = Js_obc
  OBC%segment(l_seg)%Je_obc = Je_obc
! call allocate_OBC_segment_data(OBC, OBC%segment(l_seg))

  if (OBC%segment(l_seg)%u_values_needed .or.  OBC%segment(l_seg)%v_values_needed) &
    OBC%segment(l_seg)%values_needed = .true.
end subroutine setup_u_point_obc

!> Parse an OBC_SEGMENT_%%% string starting with "J=" and configure placement
!and type of OBC accordingly
subroutine setup_v_point_obc(OBC, G, US, segment_str, l_seg, PF, reentrant_x)
  type(ice_OBC_type),   intent(inout) :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),  intent(in) :: G   !< Ocean grid structure
  type(unit_scale_type),   intent(in) :: US  !< A dimensional unit scaling type
  character(len=*),        intent(in) :: segment_str !< A string in form of "J=%,I=%:%,string"
  integer,                 intent(in) :: l_seg !< which segment is this?
  type(param_file_type),   intent(in) :: PF  !< Parameter file handle
  logical, intent(in)                 :: reentrant_x !< is the domain reentrant in x?
  ! Local variables
  integer :: J_obc, Is_obc, Ie_obc ! Position of segment in global index space
  integer :: i, a_loop
  character(len=32) :: action_str(8)
  character(len=128) :: segment_param_str
  real, allocatable, dimension(:)  :: tnudge ! Nudging timescales [T ~> s]

  ! This returns the global indices for the segment
  call parse_segment_str(G%ieg, G%jeg, segment_str, J_obc, Is_obc, Ie_obc, action_str, reentrant_x)

  call setup_segment_indices(G, OBC%segment(l_seg),Is_obc,Ie_obc,J_obc,J_obc)
  J_obc = J_obc - G%jdg_offset ! Convert to local tile indices on this tile
  Is_obc = Is_obc - G%idg_offset ! Convert to local tile indices on this tile
  Ie_obc = Ie_obc - G%idg_offset ! Convert to local tile indices on this tile

  if (Ie_obc>Is_obc) then
    OBC%segment(l_seg)%direction = OBC_DIRECTION_S
  elseif (Ie_obc<Is_obc) then
    OBC%segment(l_seg)%direction = OBC_DIRECTION_N
    i = Is_obc ; Is_obc = Ie_obc ; Ie_obc = i
  endif

  OBC%segment(l_seg)%on_pe = .false.

  do a_loop = 1,8
    if (len_trim(action_str(a_loop)) == 0) then
      cycle
    elseif (trim(action_str(a_loop)) == 'FLATHER') then
      OBC%segment(l_seg)%Flather = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%Flather_v_BCs_exist_globally = .true.
      OBC%open_v_BCs_exist_globally = .true.
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_v_BCs_exist_globally = .true.
      OBC%radiation_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI_TAN') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%radiation_tan = .true.
      OBC%radiation_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'ORLANSKI_GRAD') then
      OBC%segment(l_seg)%radiation = .true.
      OBC%segment(l_seg)%radiation_grad = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED') then
      OBC%segment(l_seg)%nudged = .true.
      OBC%nudged_v_BCs_exist_globally = .true.
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_TAN') then
      OBC%segment(l_seg)%nudged_tan = .true.
      OBC%nudged_v_BCs_exist_globally = .true.
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'NUDGED_GRAD') then
      OBC%segment(l_seg)%nudged_grad = .true.
      OBC%segment(l_seg)%g_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'GRADIENT') then
      OBC%segment(l_seg)%gradient = .true.
      OBC%segment(l_seg)%open = .true.
      OBC%open_v_BCs_exist_globally = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE') then
      OBC%segment(l_seg)%specified = .true.
      OBC%specified_v_BCs_exist_globally = .true. ! This avoids deallocation
      OBC%segment(l_seg)%v_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_TAN') then
      OBC%segment(l_seg)%specified_tan = .true.
      OBC%segment(l_seg)%u_values_needed = .true.
    elseif (trim(action_str(a_loop)) == 'SIMPLE_GRAD') then
      OBC%segment(l_seg)%specified_grad = .true.
      OBC%segment(l_seg)%g_values_needed = .true.
    else
      call MOM_error(FATAL, "SIS_open_boundary.F90, setup_v_point_obc: "//&
                     "String '"//trim(action_str(a_loop))//"' not understood.")
    endif
    if (OBC%segment(l_seg)%nudged .or. OBC%segment(l_seg)%nudged_tan) then
      write(segment_param_str(1:43),"('OBC_SEGMENT_',i3.3,'_VELOCITY_NUDGING_TIMESCALES')") l_seg
      allocate(tnudge(2))
      call get_param(PF, mdl, segment_param_str(1:43), tnudge, &
                     "Timescales in days for nudging along a segment, "//&
                     "for inflow, then outflow. Setting both to zero should "//&
                     "behave like SIMPLE obcs for the baroclinic velocities.", &
                     fail_if_missing=.true., default=0., units="days", scale=86400.0*US%s_to_T)
      OBC%segment(l_seg)%Velocity_nudging_timescale_in = tnudge(1)
      OBC%segment(l_seg)%Velocity_nudging_timescale_out = tnudge(2)
      deallocate(tnudge)
    endif

  enddo ! a_loop

  if (J_obc<=G%HI%JsdB+1 .or. J_obc>=G%HI%JedB-1) return ! Boundary is not on tile
  if (Ie_obc<=G%HI%IsdB .or. Is_obc>=G%HI%IedB) return ! Segment is not on tile

  OBC%segment(l_seg)%on_pe = .true.
  OBC%segment(l_seg)%is_N_or_S = .true.

  do i=G%HI%isd, G%HI%ied
    if (i>Is_obc .and. i<=Ie_obc) then
      OBC%segnum_v(i,J_obc) = l_seg
    endif
  enddo
  OBC%segment(l_seg)%Is_obc = Is_obc
  OBC%segment(l_seg)%Ie_obc = Ie_obc
  OBC%segment(l_seg)%Js_obc = J_obc
  OBC%segment(l_seg)%Je_obc = J_obc
! call allocate_OBC_segment_data(OBC, OBC%segment(l_seg))

  if (OBC%segment(l_seg)%u_values_needed .or.  OBC%segment(l_seg)%v_values_needed) &
    OBC%segment(l_seg)%values_needed = .true.
end subroutine setup_v_point_obc

!> Reconcile masks and open boundaries, deallocate OBC on PEs where it is not needed.
!! Also adjust u- and v-point cell area on specified open boundaries and mask all
!! points outside open boundaries.
subroutine open_boundary_impose_land_mask(OBC, G, areaCu, areaCv, US)
  type(ice_OBC_type),                pointer       :: OBC !< Open boundary control structure
  type(dyn_horgrid_type),            intent(inout) :: G   !< Ocean grid structure
  type(unit_scale_type),             intent(in)    :: US  !< A dimensional unit scaling type
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: areaCu !< Area of a u-cell [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: areaCv !< Area of a u-cell [L2 ~> m2]
  ! Local variables
  integer :: i, j, n
  type(OBC_segment_type), pointer :: segment => NULL()
  logical :: any_U, any_V

  if (.not.associated(OBC)) return

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%is_E_or_W) then
      ! Sweep along u-segments and delete the OBC for blocked points.
      ! Also, mask all points outside.
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        if (G%mask2dCu(I,j) == 0) OBC%segnum_u(I,j) = OBC_NONE
        if (segment%direction == OBC_DIRECTION_W) then
          G%mask2dT(i,j) = 0
        else
          G%mask2dT(i+1,j) = 0
        endif
      enddo
      do J=segment%HI%JsdB+1,segment%HI%JedB-1
        if (segment%direction == OBC_DIRECTION_W) then
          G%mask2dCv(i,J) = 0
        else
          G%mask2dCv(i+1,J) = 0
        endif
      enddo
    else
      ! Sweep along v-segments and delete the OBC for blocked points.
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (G%mask2dCv(i,J) == 0) OBC%segnum_v(i,J) = OBC_NONE
        if (segment%direction == OBC_DIRECTION_S) then
          G%mask2dT(i,j) = 0
        else
          G%mask2dT(i,j+1) = 0
        endif
      enddo
      do I=segment%HI%IsdB+1,segment%HI%IedB-1
        if (segment%direction == OBC_DIRECTION_S) then
          G%mask2dCu(I,j) = 0
        else
          G%mask2dCu(I,j+1) = 0
        endif
      enddo
    endif
  enddo

  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe .or. .not. segment%specified) cycle
    if (segment%is_E_or_W) then
      ! Sweep along u-segments and for %specified BC points reset the u-point
      ! area which was masked out
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        if (segment%direction == OBC_DIRECTION_E) then
          areaCu(I,j) = G%areaT(i,j)   ! Both of these are in [L2 ~> m2]
        else   ! West
          areaCu(I,j) = G%areaT(i+1,j) ! Both of these are in [L2 ~> m2]
        endif
      enddo
    else
      ! Sweep along v-segments and for %specified BC points reset the v-point
      ! area which was masked out
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (segment%direction == OBC_DIRECTION_S) then
          areaCv(i,J) = G%areaT(i,j+1) ! Both of these are in [L2 ~> m2]
        else      ! North
          areaCu(i,J) = G%areaT(i,j)   ! Both of these are in [L2 ~> m2]
        endif
      enddo
    endif
  enddo

  ! G%mask2du will be open wherever bathymetry allows it.
  ! Bathymetry outside of the open boundary was adjusted to match
  ! the bathymetry inside so these points will be open unless the
  ! bathymetry inside the boundary was too shallow and flagged as land.
  any_U = .false.
  any_V = .false.
  do n=1,OBC%number_of_segments
    segment=>OBC%segment(n)
    if (.not. segment%on_pe) cycle
    if (segment%is_E_or_W) then
      I=segment%HI%IsdB
      do j=segment%HI%jsd,segment%HI%jed
        if (OBC%segnum_u(I,j) /= OBC_NONE) any_U = .true.
      enddo
    else
      J=segment%HI%JsdB
      do i=segment%HI%isd,segment%HI%ied
        if (OBC%segnum_v(i,J) /= OBC_NONE) any_V = .true.
      enddo
    endif
  enddo

  OBC%OBC_pe = .true.
  if (.not.(any_U .or. any_V)) OBC%OBC_pe = .false.

end subroutine open_boundary_impose_land_mask

!> Find the region outside of all open boundary segments and
!! make sure it is set to land mask. Gonna need to know global land
!! mask as well to get it right...
subroutine mask_outside_OBCs(G, US, param_file, OBC)
  type(dyn_horgrid_type),       intent(inout) :: G          !< Ocean grid structure
  type(param_file_type),        intent(in)    :: param_file !< Parameter file handle
  type(ice_OBC_type),           pointer       :: OBC        !< Open boundary structure
  type(unit_scale_type),        intent(in)    :: US         !< A dimensional unit scaling type

  ! Local variables
  integer :: isd, ied, IsdB, IedB, jsd, jed, JsdB, JedB, n
  integer :: i, j
  integer :: l_seg
  logical :: fatal_error = .False.
  real    :: min_depth ! The minimum depth for ocean points [Z ~> m]
  integer, parameter :: cin = 3, cout = 4, cland = -1, cedge = -2
  character(len=256) :: mesg    ! Message for error messages.
  type(OBC_segment_type), pointer :: segment => NULL() ! pointer to segment type list
  real, allocatable, dimension(:,:) :: color, color2  ! For sorting inside from outside,
                                                      ! two different ways

  if (.not. associated(OBC)) return

  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 units="m", default=0.0, scale=US%m_to_Z, do_not_log=.true.)
  ! The reference depth on a dyn_horgrid is 0, otherwise would need:
  !   min_depth = min_depth - G%Z_ref

  allocate(color(G%isd:G%ied, G%jsd:G%jed), source=0.0)
  allocate(color2(G%isd:G%ied, G%jsd:G%jed), source=0.0)

  ! Paint a frame around the outside.
  do j=G%jsd,G%jed
    color(G%isd,j) = cedge
    color(G%ied,j) = cedge
    color2(G%isd,j) = cedge
    color2(G%ied,j) = cedge
  enddo
  do i=G%isd,G%ied
    color(i,G%jsd) = cedge
    color(i,G%jed) = cedge
    color2(i,G%jsd) = cedge
    color2(i,G%jed) = cedge
  enddo

  ! Set color to cland in the land. Note that this is before the land
  ! mask has been initialized, set mask values based on depth.
  do j=G%jsd,G%jed
    do i=G%isd,G%ied
      if (G%bathyT(i,j) <= min_depth) color(i,j) = cland
      if (G%bathyT(i,j) <= min_depth) color2(i,j) = cland
    enddo
  enddo

  do j=G%jsd,G%jed ; do i=G%IsdB+1,G%IedB-1
    l_seg = OBC%segnum_u(I,j)
    if (l_seg == OBC_NONE) cycle

    if (OBC%segment(l_seg)%direction == OBC_DIRECTION_W) then
      if (color(i,j) == 0.0) color(i,j) = cout
      if (color(i+1,j) == 0.0) color(i+1,j) = cin
    elseif (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
      if (color(i,j) == 0.0) color(i,j) = cin
      if (color(i+1,j) == 0.0) color(i+1,j) = cout
    endif
  enddo ; enddo
  do J=G%JsdB+1,G%JedB-1 ; do i=G%isd,G%ied
    l_seg = OBC%segnum_v(i,J)
    if (l_seg == OBC_NONE) cycle

    if (OBC%segment(l_seg)%direction == OBC_DIRECTION_S) then
      if (color(i,j) == 0.0) color(i,j) = cout
      if (color(i,j+1) == 0.0) color(i,j+1) = cin
    elseif (OBC%segment(l_seg)%direction == OBC_DIRECTION_N) then
      if (color(i,j) == 0.0) color(i,j) = cin
      if (color(i,j+1) == 0.0) color(i,j+1) = cout
    endif
  enddo ; enddo

  do J=G%JsdB+1,G%JedB-1 ; do i=G%isd,G%ied
    l_seg = OBC%segnum_v(i,J)
    if (l_seg == OBC_NONE) cycle

    if (OBC%segment(l_seg)%direction == OBC_DIRECTION_S) then
      if (color2(i,j) == 0.0) color2(i,j) = cout
      if (color2(i,j+1) == 0.0) color2(i,j+1) = cin
    elseif (OBC%segment(l_seg)%direction == OBC_DIRECTION_N) then
      if (color2(i,j) == 0.0) color2(i,j) = cin
      if (color2(i,j+1) == 0.0) color2(i,j+1) = cout
    endif
  enddo ; enddo
  do j=G%jsd,G%jed ; do i=G%IsdB+1,G%IedB-1
    l_seg = OBC%segnum_u(I,j)
    if (l_seg == OBC_NONE) cycle

    if (OBC%segment(l_seg)%direction == OBC_DIRECTION_W) then
      if (color2(i,j) == 0.0) color2(i,j) = cout
      if (color2(i+1,j) == 0.0) color2(i+1,j) = cin
    elseif (OBC%segment(l_seg)%direction == OBC_DIRECTION_E) then
      if (color2(i,j) == 0.0) color2(i,j) = cin
      if (color2(i+1,j) == 0.0) color2(i+1,j) = cout
    endif
  enddo ; enddo

  ! Do the flood fill until there are no more uncolored cells.
  call flood_fill(G, color, cin, cout, cland)
  call flood_fill2(G, color2, cin, cout, cland)

  ! Use the color to set outside to min_depth on this process.
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    if (color(i,j) /= color2(i,j)) then
      fatal_error = .True.
      write(mesg,'("SIS_open_boundary: problem with OBC segments specification at ",I5,",",I5," during\n", &
          "the masking of the outside grid points.")') i, j
      call MOM_error(WARNING,"MOM register_tracer: "//mesg, all_print=.true.)
    endif
    if (color(i,j) == cout) G%bathyT(i,j) = min_depth
  enddo ; enddo
  if (fatal_error) call MOM_error(FATAL, &
      "SIS_open_boundary: inconsistent OBC segments.")

  deallocate(color)
  deallocate(color2)
end subroutine mask_outside_OBCs

!> Deallocate open boundary data
subroutine open_boundary_dealloc(OBC)
  type(ice_OBC_type), pointer :: OBC !< Open boundary control structure
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: n

  if (.not. associated(OBC)) return

  do n=1, OBC%number_of_segments
    segment => OBC%segment(n)
!   call deallocate_OBC_segment_data(segment)
  enddo
  if (allocated(OBC%segment)) deallocate(OBC%segment)
  if (allocated(OBC%segnum_u)) deallocate(OBC%segnum_u)
  if (allocated(OBC%segnum_v)) deallocate(OBC%segnum_v)
  if (allocated(OBC%tres_x)) deallocate(OBC%tres_x)
  if (allocated(OBC%tres_y)) deallocate(OBC%tres_y)
  deallocate(OBC)
end subroutine open_boundary_dealloc

end module SIS_open_boundary
