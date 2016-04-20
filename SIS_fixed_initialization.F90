!> Initializes fixed aspects of the model, such as horizontal grid metrics,
!! topography and Coriolis, in a way that is very similar to MOM6.
module SIS_fixed_initialization

  use constants_mod, only : omega

use SIS_hor_grid, only : SIS_hor_grid_type
use SIS_grid_initialize, only : initialize_SIS_masks

use MOM_domains, only : pass_var, pass_vector, sum_across_PEs, broadcast
use MOM_domains, only : root_PE, To_All, SCALAR_PAIR, CGRID_NE, AGRID
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_io, only : read_data, slasher, file_exists
use MOM_io, only : CORNER, NORTH_FACE, EAST_FACE

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_initialize_fixed

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_initialize_fixed sets up time-invariant quantities related to SIS's
!!   horizontal grid, bathymetry, restricted channel widths and the Coriolis parameter.
subroutine SIS_initialize_fixed(G, PF)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The model's horizontal grid structure.
  type(param_file_type),   intent(in)    :: PF  !< A structure indicating the open file
                                                !! to parse for model parameter values.

  real :: pi ! pi = 3.1415926... calculated as 4*atan(1)
  real    :: angle, lon_scale
  integer :: i, j

  character(len=200) :: mesg
  character(len=200) :: filename, topo_file, inputdir ! Strings for file/path
  character(len=200) :: topo_varname                  ! Variable name in file
  character(len=40)  :: mod_nm  = "SIS_fixed_init" ! This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

! Initialize the topography.

  ! Replace these with properly parsed input parameter calls.
!  call get_param(PF, mod, "INPUTDIR", inputdir, default=".")
!  call get_param(PF, mod, "TOPO_FILE", topo_file, &
!                 "The file from which the bathymetry is read.", &
!                 default="topog.nc")
!  call get_param(PF, mod, "TOPO_VARNAME", topo_varname, &
!                 "The name of the bathymetry variable in TOPO_FILE.", &
!                 default="depth")

  inputdir = "INPUT" ; topo_file = "topog.nc"
  topo_varname = "depth"
  
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(topo_file)

!  call log_param(PF, mod, "INPUTDIR/TOPO_FILE", filename)
!  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
!       " initialize_topography_from_file: Unable to open "//trim(filename))

  call read_data(filename,trim(topo_varname), G%bathyT, &
                 domain=G%Domain%mpp_domain)
!  call apply_topography_edits_from_file(D, G, PF)

  ! Initialize the various masks and any masked metrics.
  call initialize_SIS_masks(G, PF)

  ! This is where any channel information might be applied.

  pi = 4.0*atan(1.0)
  ! Set up the Coriolis parameter.
  do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
    G%CoriolisBu(I,J) = 2*omega*sin(G%geoLatBu(I,J)*pi/180)
  enddo ; enddo

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    lon_scale    = cos((G%geoLatBu(I-1,J-1) + G%geoLatBu(I,J-1  ) + &
                        G%geoLatBu(I-1,J) + G%geoLatBu(I,J)) * atan(1.0)/180)
    angle        = atan2((G%geoLonBu(I-1,J) + G%geoLonBu(I,J) - &
                          G%geoLonBu(I-1,J-1) - G%geoLonBu(I,J-1))*lon_scale, &
                          G%geoLatBu(I-1,J) + G%geoLatBu(I,J) - &
                          G%geoLatBu(I-1,J-1) - G%geoLatBu(I,J-1) )
    G%sin_rot(i,j) = sin(angle) ! angle is the clockwise angle from lat/lon to ocean
    G%cos_rot(i,j) = cos(angle) ! grid (e.g. angle of ocean "north" from true north)
  enddo ; enddo

  ! ### THIS DOESN'T SEEM RIGHT AT A CUBED-SPHERE FOLD -RWH
  call pass_var(G%cos_rot, G%Domain)
  call pass_var(G%sin_rot, G%Domain)

end subroutine SIS_initialize_fixed

end module SIS_fixed_initialization
