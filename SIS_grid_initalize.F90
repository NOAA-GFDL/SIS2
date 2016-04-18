!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_hor_grid - sets up grid and processor domains and a wide variety of  !
!   metric terms in a way that is very similar to MOM6. - Robert Hallberg      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_grid_initialize

  use constants_mod, only : omega, pi

use SIS_hor_grid, only : SIS_hor_grid_type, set_hor_grid, SIS_hor_grid_end

use MOM_domains, only : pass_var, pass_vector, pe_here, root_PE, broadcast
use MOM_domains, only : AGRID, BGRID_NE, CGRID_NE, To_All, Scalar_Pair
use MOM_domains, only : To_North, To_South, To_East, To_West
use MOM_domains, only : MOM_define_domain, MOM_define_IO_domain
use MOM_domains, only : MOM_domain_type
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_io, only : read_data, slasher, file_exists
use MOM_io, only : CORNER, NORTH_FACE, EAST_FACE

use mpp_domains_mod, only : mpp_get_domain_extents, mpp_deallocate_domain

implicit none ; private

#include <SIS2_memory.h>

public :: initialize_fixed_SIS_grid, Adcroft_reciprocal

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_fixed_SIS_grid is used to set the primary values in the model's horizontal
!!   grid.  The bathymetry, land-sea mask and any restricted channel widths are
!!   not know yet, so these are set later.
subroutine initialize_fixed_SIS_grid(G, param_file)
  type(SIS_hor_grid_type), intent(inout) :: G
  type(param_file_type)  , intent(in)    :: param_file
!   This subroutine sets up the necessary domain types and the sea-ice grid.

! Arguments: G - The sea-ice's horizontal grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This include declares and sets the variable "version".
#include "version_variable.h"

  real    :: angle, lon_scale
  integer :: i, j

  character(len=200) :: mesg
  character(len=200) :: filename, topo_file, inputdir ! Strings for file/path
  character(len=200) :: topo_varname                  ! Variable name in file
  character(len=40)  :: mod_nm  = "SIS_fixed_init" ! This module's name.

  ! Read all relevant parameters and write them to the model log.
!  call log_version(param_file, mod_nm, version)


  call set_grid_metrics_from_mosaic(G, param_file)

! Calculate derived metrics (i.e. reciprocals and products)
  call set_grid_derived_metrics(G, param_file)


! Initialize the topography.

  ! Replace these with properly parsed input parameter calls.
!  call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
!  call get_param(param_file, mod, "TOPO_FILE", topo_file, &
!                 "The file from which the bathymetry is read.", &
!                 default="topog.nc")
!  call get_param(param_file, mod, "TOPO_VARNAME", topo_varname, &
!                 "The name of the bathymetry variable in TOPO_FILE.", &
!                 default="depth")

  inputdir = "INPUT" ; topo_file = "topog.nc"
  topo_varname = "depth"
  
  inputdir = slasher(inputdir)
  filename = trim(inputdir)//trim(topo_file)

!  call log_param(param_file, mod, "INPUTDIR/TOPO_FILE", filename)
!  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
!       " initialize_topography_from_file: Unable to open "//trim(filename))

  call read_data(filename,trim(topo_varname), G%bathyT, &
                 domain=G%Domain%mpp_domain)
!  call apply_topography_edits_from_file(D, G, param_file)

  ! Initialize the various masks and any masked metrics.
  call initialize_SIS_masks(G, param_file)

  ! This is where any channel information might be applied.

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

end subroutine initialize_fixed_SIS_grid


!> set_grid_derived_metrics is sets additional grid metrics that can be derived
!!   from the basic grid lengths and areas.
subroutine set_grid_derived_metrics(G, param_file)
  type(SIS_hor_grid_type), intent(inout) :: G           !< The horizontal grid structure
  type(param_file_type),   intent(in)    :: param_file  !< Parameter file structure

!    Calculate the values of the metric terms that might be used
!  and save them in arrays.  This should be identical to the corresponding
!  subroutine in MOM_grid_initialize.F90.
!    Within this subroutine, the x- and y- grid spacings and their
!  inverses and the cell areas centered on T, Bu, Cu, and Cv points are
!  calculated, as are the geographic locations of each of these 4
!  sets of points.
  character( len = 128) :: warnmesg
  integer :: i, j, isd, ied, jsd, jed
  integer :: IsdB, IedB, JsdB, JedB

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call MOM_mesg("  MOM_grid_init.F90, set_grid_derived_metrics: deriving metrics", 5)

  do j=jsd,jed ; do i=isd,ied
    if (G%dxT(i,j) < 0.0) then
      write(warnmesg,68)  pe_here(),"dxT",i,j,G%dxT(i,j),0.0
      call MOM_mesg(warnmesg, all_print=.true.)
      G%dxT(i,j) = 0.0
    endif
    if (G%dyT(i,j) < 0.0) then
      write(warnmesg,68)  pe_here(),"dyT",i,j,G%dyT(i,j),0.0
      call MOM_mesg(warnmesg, all_print=.true.)
      G%dyT(i,j) = 0.0
    endif
    G%IdxT(i,j) = Adcroft_reciprocal(G%dxT(i,j))
    G%IdyT(i,j) = Adcroft_reciprocal(G%dyT(i,j))
    G%IareaT(i,j) = Adcroft_reciprocal(G%areaT(i,j))
  enddo ; enddo

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dxCu(I,j) < 0.0) then
      write(warnmesg,68)  pe_here(),"dxCu",I,j,G%dxCu(I,j),0.0
      call MOM_mesg(warnmesg, all_print=.true.)
      G%dxCu(I,j) = 0.0
    endif
    if (G%dyCu(I,j) < 0.0) then
      write(warnmesg,68)  pe_here(),"dyCu",I,j,G%dyCu(I,j),0.0
      call MOM_mesg(warnmesg, all_print=.true.)
      G%dyCu(I,j) = 0.0
    endif
    G%IdxCu(I,j) = Adcroft_reciprocal(G%dxCu(I,j))
    G%IdyCu(I,j) = Adcroft_reciprocal(G%dyCu(I,j))
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dxCv(i,J) < 0.0) then
      write(warnmesg,68)  pe_here(),"dxCv",i,J,G%dxCv(i,J),0.0
      call MOM_mesg(warnmesg, all_print=.true.)
      G%dxCv(i,J) = 0.0
    endif
    if (G%dyCv(i,J) < 0.0) then
      write(warnmesg,68)  pe_here(),"dyCv",i,J,G%dyCv(i,J),0.0
      call MOM_mesg(warnmesg, all_print=.true.)
      G%dyCv(i,J) = 0.0
    endif
    G%IdxCv(i,J) = Adcroft_reciprocal(G%dxCv(i,J))
    G%IdyCv(i,J) = Adcroft_reciprocal(G%dyCv(i,J))
  enddo ; enddo

  do J=JsdB,JedB ; do I=IsdB,IedB
    if (G%dxBu(I,J) < 0.0) then
      write(warnmesg,68)  pe_here(),"dxBu",I,J,G%dxBu(I,J),0.0
      call MOM_mesg(warnmesg, all_print=.true.)
      G%dxBu(I,J) = 0.0
    endif
    if (G%dyBu(I,J) < 0.0) then
      write(warnmesg,68)  pe_here(),"dyBu",I,J,G%dyBu(I,J),0.0
      call MOM_mesg(warnmesg, all_print=.true.)
      G%dyBu(I,J) = 0.0
    endif

    G%IdxBu(I,J) = Adcroft_reciprocal(G%dxBu(I,J))
    G%IdyBu(I,J) = Adcroft_reciprocal(G%dyBu(I,J))
    ! areaBu has usually been set to a positive area elsewhere.
    if (G%areaBu(I,J) <= 0.0) G%areaBu(I,J) = G%dxBu(I,J) * G%dyBu(I,J)
    G%IareaBu(I,J) = Adcroft_reciprocal(G%areaBu(I,J))
  enddo ; enddo

68 FORMAT ("WARNING: PE ",I4," ",a4,"(",I4,",",I4,") = ",ES12.4, &
           " is being changed to ",ES12.4,".")

end subroutine set_grid_derived_metrics


subroutine set_grid_metrics_from_mosaic(G,param_file)
  type(SIS_hor_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
!   This subroutine sets the grid metrics from a mosaic file.  This should be
! identical to the corresponding subroutine in MOM_grid_initialize.F90.
!
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

  real, dimension(G%isd :G%ied ,G%jsd :G%jed ) :: tempH1, tempH2, tempH3, tempH4
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB) :: tempQ1, tempQ2, tempQ3, tempQ4
  real, dimension(G%IsdB:G%IedB,G%jsd :G%jed ) :: tempE1, tempE2
  real, dimension(G%isd :G%ied ,G%JsdB:G%JedB) :: tempN1, tempN2
  ! These arrays are a holdover from earlier code in which the arrays in G were
  ! macros and may have had reduced dimensions.
  real, dimension(G%isd :G%ied ,G%jsd :G%jed ) :: dxT, dyT, areaT
  real, dimension(G%IsdB:G%IedB,G%jsd :G%jed ) :: dxCu, dyCu
  real, dimension(G%isd :G%ied ,G%JsdB:G%JedB) :: dxCv, dyCv
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB) :: dxBu, dyBu, areaBu
  ! This are symmetric arrays, corresponding to the data in the mosaic file
  real, dimension(2*G%isd-2:2*G%ied+1,2*G%jsd-2:2*G%jed+1) :: tmpT
  real, dimension(2*G%isd-3:2*G%ied+1,2*G%jsd-2:2*G%jed+1) :: tmpU
  real, dimension(2*G%isd-2:2*G%ied+1,2*G%jsd-3:2*G%jed+1) :: tmpV
  real, dimension(2*G%isd-3:2*G%ied+1,2*G%jsd-3:2*G%jed+1) :: tmpZ
  real, dimension(:,:), allocatable :: tmpGlbl
  character(len=200) :: filename, grid_file, inputdir
  character(len=64)  :: mod="MOM_grid_init set_grid_metrics_from_mosaic"
  integer :: err=0, ni, nj, global_indices(4)
  type(MOM_domain_type) :: SGdom ! Supergrid domain
  integer :: i, j, i2, j2
  integer :: npei,npej
  integer, dimension(:), allocatable :: exni,exnj
  integer        :: start(4), nread(4)

  call MOM_mesg("   MOM_grid_init.F90, set_grid_metrics_from_mosaic: reading grid", 5)

  call get_param(param_file, mod, "GRID_FILE", grid_file, &
                 "Name of the file from which to read horizontal grid data.", &
                 fail_if_missing=.true.)
  call get_param(param_file,  mod, "INPUTDIR", inputdir, &
         "The directory in which input files are found.", default=".")
  inputdir = slasher(inputdir)
  filename = trim(adjustl(inputdir)) // trim(adjustl(grid_file))
  call log_param(param_file, mod, "INPUTDIR/GRID_FILE", filename)
  if (.not.file_exists(filename)) &
    call MOM_error(FATAL," set_grid_metrics_from_mosaic: Unable to open "//&
                           trim(filename))

! Initialize everything to 0.
  dxCu(:,:) = 0.0 ; dyCu(:,:) = 0.0
  dxCv(:,:) = 0.0 ; dyCv(:,:) = 0.0
  dxBu(:,:) = 0.0 ; dyBu(:,:) = 0.0 ; areaBu(:,:) = 0.0

!<MISSING CODE TO READ REFINEMENT LEVEL>
  ni = 2*(G%iec-G%isc+1) ! i size of supergrid
  nj = 2*(G%jec-G%jsc+1) ! j size of supergrid

! Define a domain for the supergrid (SGdom)
  npei = G%domain%layout(1) ; npej = G%domain%layout(2)
  allocate(exni(npei)) ; allocate(exnj(npej))
  call mpp_get_domain_extents(G%domain%mpp_domain, exni, exnj)
  allocate(SGdom%mpp_domain)
  SGdom%nihalo = 2*G%domain%nihalo+1
  SGdom%njhalo = 2*G%domain%njhalo+1
  SGdom%niglobal = 2*G%domain%niglobal
  SGdom%njglobal = 2*G%domain%njglobal
  SGdom%layout(:) = G%domain%layout(:)
  SGdom%use_io_layout = G%domain%use_io_layout
  SGdom%io_layout(:) = G%domain%io_layout(:)
  global_indices(1) = 1+SGdom%nihalo
  global_indices(2) = SGdom%niglobal+SGdom%nihalo
  global_indices(3) = 1+SGdom%njhalo
  global_indices(4) = SGdom%njglobal+SGdom%njhalo
  exni(:) = 2*exni(:) ; exnj(:) = 2*exnj(:)
  if(ASSOCIATED(G%domain%maskmap)) then
     call MOM_define_domain(global_indices, SGdom%layout, SGdom%mpp_domain, &
            xflags=G%domain%X_FLAGS, yflags=G%domain%Y_FLAGS, &
            xhalo=SGdom%nihalo, yhalo=SGdom%njhalo, &
            xextent=exni,yextent=exnj, &
            symmetry=.true., name="MOM_MOSAIC", maskmap=G%domain%maskmap)
  else
     call MOM_define_domain(global_indices, SGdom%layout, SGdom%mpp_domain, &
            xflags=G%domain%X_FLAGS, yflags=G%domain%Y_FLAGS, &
            xhalo=SGdom%nihalo, yhalo=SGdom%njhalo, &
            xextent=exni,yextent=exnj, &
            symmetry=.true., name="MOM_MOSAIC")
  endif

  if (SGdom%use_io_layout) &
    call MOM_define_IO_domain(SGdom%mpp_domain, SGdom%io_layout)
  deallocate(exni)
  deallocate(exnj)

! Read X from the supergrid
  tmpZ(:,:) = 999.
  call read_data(filename, 'x', tmpZ, domain=SGdom%mpp_domain, position=CORNER)

  call pass_var(tmpZ, SGdom, position=CORNER)
  call extrapolate_metric(tmpZ, 2*(G%jsc-G%jsd)+2, missing=999.)
  do j=G%jsd,G%jed ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
    G%geoLonT(i,j) = tmpZ(i2-1,j2-1)
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB ; i2 = 2*I ; j2 = 2*J
    G%geoLonBu(I,J) = tmpZ(i2,j2)
  enddo ; enddo
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB ; i2 = 2*i ; j2 = 2*j
    G%geoLonCu(I,j) = tmpZ(i2,j2-1)
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*J
    G%geoLonCv(i,J) = tmpZ(i2-1,j2)
  enddo ; enddo
 ! For some reason, this messes up the solution...
 !   call pass_var(G%geoLonBu, G%domain, position=CORNER)

! Read Y from the supergrid
  tmpZ(:,:) = 999.
  call read_data(filename, 'y', tmpZ, domain=SGdom%mpp_domain, position=CORNER)

  call pass_var(tmpZ, SGdom, position=CORNER)
  call extrapolate_metric(tmpZ, 2*(G%jsc-G%jsd)+2, missing=999.)
  do j=G%jsd,G%jed ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
    G%geoLatT(i,j) = tmpZ(i2-1,j2-1)
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB ; i2 = 2*I ; j2 = 2*J
    G%geoLatBu(I,J) = tmpZ(i2,j2)
  enddo ; enddo
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB ; i2 = 2*i ; j2 = 2*j
    G%geoLatCu(I,j) = tmpZ(i2,j2-1)
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*J
    G%geoLatCv(i,J) = tmpZ(i2-1,j2)
  enddo ; enddo

! Read DX,DY from the supergrid
  tmpU(:,:) = 0. ; tmpV(:,:) = 0.
  call read_data(filename,'dx',tmpV,domain=SGdom%mpp_domain,position=NORTH_FACE)
  call read_data(filename,'dy',tmpU,domain=SGdom%mpp_domain,position=EAST_FACE)
  call pass_vector(tmpU, tmpV, SGdom, To_All+Scalar_Pair, CGRID_NE)
  call extrapolate_metric(tmpV, 2*(G%jsc-G%jsd)+2, missing=0.)
  call extrapolate_metric(tmpU, 2*(G%jsc-G%jsd)+2, missing=0.)

  do j=G%jsd,G%jed ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
    dxT(i,j) = tmpV(i2-1,j2-1) + tmpV(i2,j2-1)
    dyT(i,j) = tmpU(i2-1,j2-1) + tmpU(i2-1,j2)
  enddo ; enddo

  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB ; i2 = 2*i ; j2 = 2*j
    dxCu(I,j) = tmpV(i2,j2-1) + tmpV(i2+1,j2-1)
    dyCu(I,j) = tmpU(i2,j2-1) + tmpU(i2,j2)
  enddo ; enddo

  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
    dxCv(i,J) = tmpV(i2-1,j2) + tmpV(i2,j2)
    dyCv(i,J) = tmpU(i2-1,j2) + tmpU(i2-1,j2+1)
  enddo ; enddo

  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB ; i2 = 2*i ; j2 = 2*j
    dxBu(I,J) = tmpV(i2,j2) + tmpV(i2+1,j2)
    dyBu(I,J) = tmpU(i2,j2) + tmpU(i2,j2+1)
  enddo ; enddo

! Read AREA from the supergrid
  tmpT(:,:) = 0.
  call read_data(filename, 'area', tmpT, domain=SGdom%mpp_domain)
  call pass_var(tmpT, SGdom)
  call extrapolate_metric(tmpT, 2*(G%jsc-G%jsd)+2, missing=0.)

  do j=G%jsd,G%jed ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
    areaT(i,j) = (tmpT(i2-1,j2-1) + tmpT(i2,j2)) + &
                 (tmpT(i2-1,j2) + tmpT(i2,j2-1))
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB ; i2 = 2*i ; j2 = 2*j
    areaBu(i,j) = (tmpT(i2,j2) + tmpT(i2+1,j2+1)) + &
                  (tmpT(i2,j2+1) + tmpT(i2+1,j2))
  enddo ; enddo

  ni=SGdom%niglobal
  nj=SGdom%njglobal
  call mpp_deallocate_domain(SGdom%mpp_domain)
  deallocate(SGdom%mpp_domain)

  call pass_vector(dyCu, dxCv, G%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dxCu, dyCv, G%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dxBu, dyBu, G%Domain, To_All+Scalar_Pair, BGRID_NE)
  call pass_var(areaT, G%Domain)
  call pass_var(areaBu, G%Domain, position=CORNER)

  do i=G%isd,G%ied ; do j=G%jsd,G%jed
    G%dxT(i,j) = dxT(i,j) ; G%dyT(i,j) = dyT(i,j) ; G%areaT(i,j) = areaT(i,j)
  enddo ; enddo
  do I=G%IsdB,G%IedB ; do j=G%jsd,G%jed
    G%dxCu(I,j) = dxCu(I,j) ; G%dyCu(I,j) = dyCu(I,j)
  enddo ; enddo
  do i=G%isd,G%ied ; do J=G%JsdB,G%JedB
    G%dxCv(i,J) = dxCv(i,J) ; G%dyCv(i,J) = dyCv(i,J)
  enddo ; enddo
  do I=G%IsdB,G%IedB ; do J=G%JsdB,G%JedB
    G%dxBu(I,J) = dxBu(I,J) ; G%dyBu(I,J) = dyBu(I,J) ; G%areaBu(I,J) = areaBu(I,J)
  enddo ; enddo

  ! Construct axes for diagnostic output (only necessary because "ferret" uses
  ! broken convention for interpretting netCDF files).
  start(:) = 1 ; nread(:) = 1
  start(2) = 2 ; nread(1) = ni+1 ; nread(2) = 2
  allocate( tmpGlbl(ni+1,2) )
  if (is_root_PE()) &
    call read_data(filename, "x", tmpGlbl, start, nread, no_domain=.TRUE.)
  call broadcast(tmpGlbl, 2*(ni+1), root_PE())

  ! I don't know why the second axis is 1 or 2 here. -RWH
  do i=G%isg,G%ieg
    G%gridLonT(i) = tmpGlbl(2*(i-G%isg)+2,2)
  enddo
  do I=G%IsgB,G%IegB
    G%gridLonB(I) = tmpGlbl(2*(I-G%isg)+3,1)
  enddo
  deallocate( tmpGlbl )

  allocate( tmpGlbl(1, nj+1) )
  start(:) = 1 ; nread(:) = 1
  start(1) = int(ni/4)+1 ; nread(2) = nj+1
  if (is_root_PE()) &
    call read_data(filename, "y", tmpGlbl, start, nread, no_domain=.TRUE.)
  call broadcast(tmpGlbl, nj+1, root_PE())

  do j=G%jsg,G%jeg
    G%gridLatT(j) = tmpGlbl(1,2*(j-G%jsg)+2)
  enddo
  do J=G%JsgB,G%JegB
    G%gridLatB(J) = tmpGlbl(1,2*(j-G%jsg)+3)
  enddo
  deallocate( tmpGlbl )

end subroutine set_grid_metrics_from_mosaic

! ------------------------------------------------------------------------------

!> extrapolate_metric extrapolates missing metric data into all the halo regions.
subroutine extrapolate_metric(var, jh, missing)
  real, dimension(:,:), intent(inout) :: var     !< The array in which to fill in halos
  integer,              intent(in)    :: jh      !< The size of the halos to be filled
  real,       optional, intent(in)    :: missing !< The missing data fill value, 0 by default.
  real :: badval
  integer :: i,j

  badval = 0.0 ; if (present(missing)) badval = missing

  ! Fill in southern halo by extrapolating from the computational domain
  do j=lbound(var,2)+jh,lbound(var,2),-1 ; do i=lbound(var,1),ubound(var,1)
    if (var(i,j)==badval) var(i,j) = 2.0*var(i,j+1)-var(i,j+2)
  enddo ; enddo

  ! Fill in northern halo by extrapolating from the computational domain
  do j=ubound(var,2)-jh,ubound(var,2) ; do i=lbound(var,1),ubound(var,1)
    if (var(i,j)==badval) var(i,j) = 2.0*var(i,j-1)-var(i,j-2)
  enddo ; enddo

  ! Fill in western halo by extrapolating from the computational domain
  do j=lbound(var,2),ubound(var,2) ; do i=lbound(var,1)+jh,lbound(var,1),-1
    if (var(i,j)==badval) var(i,j) = 2.0*var(i+1,j)-var(i+2,j)
  enddo ; enddo

  ! Fill in eastern halo by extrapolating from the computational domain
  do j=lbound(var,2),ubound(var,2) ; do i=ubound(var,1)-jh,ubound(var,1)
    if (var(i,j)==badval) var(i,j) = 2.0*var(i-1,j)-var(i-2,j)
  enddo ; enddo

end subroutine extrapolate_metric

!> This function implements Adcroft's rule for reciprocals, namely that
!!   Adcroft_Inv(x) = 1/x for |x|>0 or 0 for x=0.
function Adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val  !< The value being inverted.
  real :: I_val            !< The Adcroft reciprocal of val.

  I_val = 0.0
  if (val /= 0.0) I_val = 1.0/val
end function Adcroft_reciprocal

!> initialize_SIS_masks initializes the grid masks and any metrics that come
!!    with masks already applied.
subroutine initialize_SIS_masks(G, PF)
  type(SIS_hor_grid_type), intent(inout) :: G  !< The horizontal grid structure
  type(param_file_type),   intent(in)    :: PF !< The param_file handle

!    Initialize_masks sets mask2dT, mask2dCu, mask2dCv, and mask2dBu to mask out
! flow over any points which are shallower than Dmin and permit an
! appropriate treatment of the boundary conditions.  mask2dCu and mask2dCv
! are 0.0 at any points adjacent to a land point.  mask2dBu is 0.0 at
! any land or boundary point.  For points in the interior, mask2dCu,
! mask2dCv, and mask2dBu are all 1.0.

  real :: Dmin, min_depth, mask_depth
  character(len=40)  :: mod = "initialize_SIS_masks"
  integer :: i, j


  call get_param(PF, mod, "MINIMUM_DEPTH", min_depth, &
                 "If MASKING_DEPTH is unspecified, then anything shallower than\n"//&
                 "MINIMUM_DEPTH is assumed to be land and all fluxes are masked out.\n"//&
                 "If MASKING_DEPTH is specified, then all depths shallower than\n"//&
                 "MINIMUM_DEPTH but deeper than MASKING_DEPTH are rounded to MINIMUM_DEPTH.", &
                 units="m", default=0.0)
  call get_param(PF, mod, "MASKING_DEPTH", mask_depth, &
                 "The depth below which to mask points as land points, for which all\n"//&
                 "fluxes are zeroed out. MASKING_DEPTH is ignored if negative.", &
                 units="m", default=-9999.0)
  Dmin = min_depth
  if (mask_depth>=0.) Dmin = mask_depth

  call pass_var(G%bathyT, G%Domain)
  G%mask2dCu(:,:) = 0.0 ; G%mask2dCv(:,:) = 0.0 ; G%mask2dBu(:,:) = 0.0

  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    if (G%bathyT(i,j) <= Dmin) then
      G%mask2dT(i,j) = 0.0
    else
      G%mask2dT(i,j) = 1.0
    endif
  enddo ; enddo

  do j=G%jsd,G%jed ; do I=G%isd,G%ied-1
    if ((G%bathyT(i,j) <= Dmin) .or. (G%bathyT(i+1,j) <= Dmin)) then
      G%mask2dCu(I,j) = 0.0
    else
      G%mask2dCu(I,j) = 1.0
    endif
  enddo ; enddo

  do J=G%jsd,G%jed-1 ; do i=G%isd,G%ied
    if ((G%bathyT(i,j) <= Dmin) .or. (G%bathyT(i,j+1) <= Dmin)) then
      G%mask2dCv(i,J) = 0.0
    else
      G%mask2dCv(i,J) = 1.0
    endif
  enddo ; enddo

  do J=G%jsd,G%jed-1 ; do I=G%isd,G%ied-1
    if ((G%bathyT(i+1,j) <= Dmin) .or. (G%bathyT(i+1,j+1) <= Dmin) .or. &
        (G%bathyT(i,j) <= Dmin) .or. (G%bathyT(i,j+1) <= Dmin)) then
      G%mask2dBu(I,J) = 0.0
    else
      G%mask2dBu(I,J) = 1.0
    endif
  enddo ; enddo

  call pass_var(G%mask2dBu, G%Domain, position=CORNER)
  call pass_vector(G%mask2dCu, G%mask2dCv, G%Domain, To_All+Scalar_Pair, CGRID_NE)

  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
    G%dy_Cu(I,j) = G%mask2dCu(I,j) * G%dyCu(I,j)
    G%areaCu(I,j) = G%dxCu(I,j) * G%dy_Cu(I,j)
    G%IareaCu(I,j) = G%mask2dCu(I,j) * Adcroft_reciprocal(G%areaCu(I,j))
  enddo ; enddo

  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
    G%dx_Cv(i,J) = G%mask2dCv(i,J) * G%dxCv(i,J)
    G%areaCv(i,J) = G%dyCv(i,J) * G%dx_Cv(i,J)
    G%IareaCv(i,J) = G%mask2dCv(i,J) * Adcroft_reciprocal(G%areaCv(i,J))
  enddo ; enddo


end subroutine initialize_SIS_masks

end module SIS_grid_initialize
