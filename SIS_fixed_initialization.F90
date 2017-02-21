!> Initializes fixed aspects of the model, such as horizontal grid metrics,
!! topography and Coriolis, in a way that is very similar to MOM6.
module SIS_fixed_initialization

! This file is part of SIS2. See LICENSE.md for the license.

use SIS_debugging, only : hchksum, qchksum, uchksum, vchksum, chksum
use MOM_domains, only : pass_var
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, log_version, param_file_type
use MOM_grid_initialize, only : initialize_masks, set_grid_metrics
use MOM_io, only : slasher
! use MOM_shared_initialization, only : MOM_shared_init_init
use MOM_shared_initialization, only : MOM_initialize_rotation, MOM_calculate_grad_Coriolis
use MOM_shared_initialization, only : initialize_topography_from_file, apply_topography_edits_from_file
use MOM_shared_initialization, only : initialize_topography_named, limit_topography, diagnoseMaximumDepth
use MOM_shared_initialization, only : set_rotation_planetary, set_rotation_beta_plane, initialize_grid_rotation_angle
use MOM_shared_initialization, only : reset_face_lengths_named, reset_face_lengths_file, reset_face_lengths_list
use MOM_shared_initialization, only : read_face_length_list, set_velocity_depth_max, set_velocity_depth_min
use MOM_shared_initialization, only : compute_global_grid_integrals, write_ocean_geometry_file

use netcdf

implicit none ; private

public :: SIS_initialize_fixed

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_initialize_fixed sets up time-invariant quantities related to SIS's
!!   horizontal grid, bathymetry, restricted channel widths and the Coriolis parameter.
subroutine SIS_initialize_fixed(G, PF, write_geom, output_dir)
  type(dyn_horgrid_type),  intent(inout) :: G    !< The ocean's grid structure.
  type(param_file_type),   intent(in)    :: PF  !< A structure indicating the open file
                                                !! to parse for model parameter values.
  logical,                 intent(in)    :: write_geom !< If true, write grid geometry files.
  character(len=*),        intent(in)    :: output_dir !< The directory into which to write files.

  real :: pi ! pi = 3.1415926... calculated as 4*atan(1)

  character(len=200) :: inputdir   ! The directory where NetCDF input files are.
  character(len=200) :: config
  character(len=40)  :: mod = "SIS_initialize_fixed" ! This module's name.
  logical :: debug
! This include declares and sets the variable "version".
#include "version_variable.h"

  call callTree_enter("SIS_initialize_fixed(), SIS_fixed_initialization.F90")
  call log_version(PF, mod, version, "")
  call get_param(PF, mod, "DEBUG", debug, default=.false.)

  call get_param(PF, mod, "INPUTDIR", inputdir, &
         "The directory in which input files are found.", default=".")
  inputdir = slasher(inputdir)

! Set up the parameters of the physical domain (i.e. the grid), G
  call set_grid_metrics(G, PF)

! Set up the bottom depth, G%bathyT, either analytically or from a file
  call SIS_initialize_topography(G%bathyT, G%max_depth, G, PF)

  ! To initialize masks, the bathymetry in halo regions must be filled in
  call pass_var(G%bathyT, G%Domain)

! Initialize the various masks and any masked metrics.
  call initialize_masks(G, PF)

  if (debug) then
    call hchksum(G%bathyT, 'SIS_initialize_fixed: depth ', G%HI, &
                 haloshift=min(1, G%ied-G%iec, G%jed-G%jec))
    call hchksum(G%mask2dT, 'SIS_initialize_fixed: mask2dT ', G%HI)
    call uchksum(G%mask2dCu, 'SIS_initialize_fixed: mask2dCu ', G%HI)
    call vchksum(G%mask2dCv, 'SIS_initialize_fixed: mask2dCv ', G%HI)
    call qchksum(G%mask2dBu, 'SIS_initialize_fixed: mask2dBu ', G%HI)
  endif

! Modulate geometric scales according to geography.
  call get_param(PF, mod, "CHANNEL_CONFIG", config, &
                 "A parameter that determines which set of channels are \n"//&
                 "restricted to specific  widths.  Options are:\n"//&
                 " \t none - All channels have the grid width.\n"//&
                 " \t global_1deg - Sets 16 specific channels appropriate \n"//&
                 " \t\t for a 1-degree model, as used in CM2G.\n"//&
                 " \t list - Read the channel locations and widths from a \n"//&
                 " \t\t text file, like MOM_channel_list in the MOM_SIS \n"//&
                 " \t\t test case.\n"//&
                 " \t file - Read open face widths everywhere from a \n"//&
                 " \t\t NetCDF file on the model grid.", &
                 default="none")
  select case ( trim(config) )
    case ("none")
    case ("list") ; call reset_face_lengths_list(G, PF)
    case ("file") ; call reset_face_lengths_file(G, PF)
    case ("global_1deg") ; call reset_face_lengths_named(G, PF, trim(config))
    case default ; call MOM_error(FATAL, "SIS_initialize_fixed: "// &
      "Unrecognized channel configuration "//trim(config))
  end select


!    Calculate the value of the Coriolis parameter at the latitude   !
!  of the q grid points, in s-1.
  call MOM_initialize_rotation(G%CoriolisBu, G, PF)
!   Calculate the components of grad f (beta)
  call MOM_calculate_grad_Coriolis(G%dF_dx, G%dF_dy, G)
  if (debug) then
    call qchksum(G%CoriolisBu, "SIS_initialize_fixed: f ", G%HI)
    call hchksum(G%dF_dx, "SIS_initialize_fixed: dF_dx ", G%HI)
    call hchksum(G%dF_dy, "SIS_initialize_fixed: dF_dy ", G%HI)
  endif

  call initialize_grid_rotation_angle(G, PF)

! Write out all of the grid data used by this run.
  if (write_geom) call write_ocean_geometry_file(G, PF, output_dir, &
                                                 geom_file="sea_ice_geometry")

  call callTree_leave('SIS_initialize_fixed()')

end subroutine SIS_initialize_fixed

!> SIS_initialize_topography makes the appropriate call to set up the bathymetry.
!! It is very similar to MOM_initialize_topography, but with fewer options.
subroutine SIS_initialize_topography(D, max_depth, G, PF)
  type(dyn_horgrid_type),           intent(in)  :: G  !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                    intent(out) :: D  !< Ocean bottom depth in m
  type(param_file_type),            intent(in)  :: PF !< Parameter file structure
  real,                             intent(out) :: max_depth !< Maximum depth of model in m

!  This subroutine makes the appropriate call to set up the bottom depth.
!  This is a separate subroutine so that it can be made public and shared with
!  the ice-sheet code or other components.
! Set up the bottom depth, G%bathyT either analytically or from file
  character(len=40)  :: mod = "SIS_initialize_topography" ! This subroutine's name.
  character(len=200) :: config

  call get_param(PF, mod, "TOPO_CONFIG", config, &
                 "This specifies how bathymetry is specified: \n"//&
                 " \t file - read bathymetric information from the file \n"//&
                 " \t\t specified by (TOPO_FILE).\n"//&
                 " \t flat - flat bottom set to MAXIMUM_DEPTH. \n"//&
                 " \t bowl - an analytically specified bowl-shaped basin \n"//&
                 " \t\t ranging between MAXIMUM_DEPTH and MINIMUM_DEPTH. \n"//&
                 " \t spoon - a similar shape to 'bowl', but with an vertical \n"//&
                 " \t\t wall at the southern face. \n"//&
                 " \t halfpipe - a zonally uniform channel with a half-sine \n"//&
                 " \t\t profile in the meridional direction.", &
                 default="file")
! The following are options that are available with MOM6 but not SIS2.
!                 " \t benchmark - use the benchmark test case topography. \n"//&
!                 " \t DOME - use a slope and channel configuration for the \n"//&
!                 " \t\t DOME sill-overflow test case. \n"//&
!                 " \t DOME2D - use a shelf and slope configuration for the \n"//&
!                 " \t\t DOME2D gravity current/overflow test case. \n"//&
!                 " \t seamount - Gaussian bump for spontaneous motion test case.\n"//&
!                 " \t Phillips - ACC-like idealized topography used in the Phillips config.\n"//&
!                 " \t USER - call a user modified routine.", &
!                 fail_if_missing=.true.)
  max_depth = -1.e9; call read_param(PF, "MAXIMUM_DEPTH", max_depth)
  select case ( trim(config) )
    case ("file");      call initialize_topography_from_file(D, G, PF)
    case ("flat");      call initialize_topography_named(D, G, PF, config, max_depth)
    case ("spoon");     call initialize_topography_named(D, G, PF, config, max_depth)
    case ("bowl");      call initialize_topography_named(D, G, PF, config, max_depth)
    case ("halfpipe");  call initialize_topography_named(D, G, PF, config, max_depth)
    case default ;      call MOM_error(FATAL,"SIS_initialize_topography: "// &
      "Unrecognized topography setup '"//trim(config)//"'")
  end select
  if (max_depth>0.) then
    call log_param(PF, mod, "MAXIMUM_DEPTH", max_depth, &
                   "The maximum depth of the ocean.", units="m")
  else
    max_depth = diagnoseMaximumDepth(D,G)
    call log_param(PF, mod, "!MAXIMUM_DEPTH", max_depth, &
                   "The (diagnosed) maximum depth of the ocean.", units="m")
  endif
  if (trim(config) .ne. "DOME") then
    call limit_topography(D, G, PF, max_depth)
  endif

end subroutine SIS_initialize_topography

end module SIS_fixed_initialization
