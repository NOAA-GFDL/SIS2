!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_hor_grid_mod - sets up grid and processor domains and a wide variety of  !
!   metric terms in a way that is very similar to MOM6. - Robert Hallberg      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_hor_grid_mod

  use constants_mod, only : omega, pi, grav

use mpp_domains_mod, only : mpp_define_domains, FOLD_NORTH_EDGE
use mpp_domains_mod, only : domain2D, mpp_global_field, YUPDATE, XUPDATE, CORNER
use mpp_domains_mod, only : CENTER, NORTH_FACE=>NORTH, EAST_FACE=>EAST
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only : mpp_define_io_domain, mpp_copy_domain, mpp_get_global_domain
use mpp_domains_mod, only : mpp_deallocate_domain, mpp_get_pelist, mpp_get_compute_domains
use mpp_domains_mod, only : domain1D, mpp_get_domain_components

use MOM_domains, only : MOM_domain_type, pass_var, pass_vector
use MOM_domains, only : PE_here, root_PE, broadcast, MOM_domains_init, clone_MOM_domain
use MOM_domains, only : num_PEs, SCALAR_PAIR, CGRID_NE, BGRID_NE, To_All
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_obsolete_params, only : obsolete_logical
use MOM_string_functions, only : slasher

use fms_io_mod, only : file_exist
use fms_mod,    only : field_exist, field_size, read_data
use mosaic_mod, only : get_mosaic_ntiles, get_mosaic_ncontacts, get_mosaic_contact

implicit none ; private

include 'netcdf.inc'
#include <SIS2_memory.h>

public :: set_hor_grid, SIS_hor_grid_end, isPointInCell

type, public :: SIS_hor_grid_type
  type(MOM_domain_type), pointer :: Domain => NULL()
  type(MOM_domain_type), pointer :: Domain_aux => NULL() ! A non-symmetric auxiliary domain type.
  integer :: isc, iec, jsc, jec ! The range of the computational domain indices
  integer :: isd, ied, jsd, jed ! and data domain indices at tracer cell centers.
  integer :: isg, ieg, jsg, jeg ! The range of the global domain tracer cell indices.
  integer :: IscB, IecB, JscB, JecB ! The range of the computational domain indices
  integer :: IsdB, IedB, JsdB, JedB ! and data domain indices at tracer cell vertices.
  integer :: IsgB, IegB, JsgB, JegB ! The range of the global domain vertex indices.
  integer :: isd_global         ! The values of isd and jsd in the global
  integer :: jsd_global         ! (decomposition invariant) index space.

  logical :: symmetric          ! True if symmetric memory is used.
  logical :: nonblocking_updates  ! If true, non-blocking halo updates are
                                  ! allowed.  The default is .false. (for now).
  integer :: first_direction ! An integer that indicates which direction is
                             ! to be updated first in directionally split
                             ! parts of the calculation.  This can be altered
                             ! during the course of the run via calls to
                             ! set_first_direction.

  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    mask2dT, &   ! 0 for land points and 1 for ocean points on the h-grid. Nd.
    geoLatT, & ! The geographic latitude at q points in degrees of latitude or m.
    geoLonT, & ! The geographic longitude at q points in degrees of longitude or m.
    dxT, IdxT, & ! dxT is delta x at h points, in m, and IdxT is 1/dxT in m-1.
    dyT, IdyT, & ! dyT is delta y at h points, in m, and IdyT is 1/dyT in m-1.
    areaT, &     ! areaT is the area of an h-cell, in m2.
    IareaT, &    ! IareaT = 1/areaT, in m-2.
    sin_rot, &   ! The sine and cosine of the angular rotation between the local
    cos_rot      ! model grid's northward and the true northward directions.

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    mask2dCu, &  ! 0 for boundary points and 1 for ocean points on the u grid.  Nondim.
    geoLatCu, &  ! The geographic latitude at u points in degrees of latitude or m.
    geoLonCu, &  ! The geographic longitude at u points in degrees of longitude or m.
    dxCu, IdxCu, & ! dxCu is delta x at u points, in m, and IdxCu is 1/dxCu in m-1.
    dyCu, IdyCu, & ! dyCu is delta y at u points, in m, and IdyCu is 1/dyCu in m-1.
    dy_Cu, &     ! The unblocked lengths of the u-faces of the h-cell in m.
    dy_Cu_obc, & ! The unblocked lengths of the u-faces of the h-cell in m for OBC.
    IareaCu, &   ! The masked inverse areas of u-grid cells in m2.
    areaCu       ! The areas of the u-grid cells in m2.

  real ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    mask2dCv, &  ! 0 for boundary points and 1 for ocean points on the v grid.  Nondim.
    geoLatCv, &  ! The geographic latitude at v points in degrees of latitude or m.
    geoLonCv, &  !  The geographic longitude at v points in degrees of longitude or m.
    dxCv, IdxCv, & ! dxCv is delta x at v points, in m, and IdxCv is 1/dxCv in m-1.
    dyCv, IdyCv, & ! dyCv is delta y at v points, in m, and IdyCv is 1/dyCv in m-1.
    dx_Cv, &     ! The unblocked lengths of the v-faces of the h-cell in m.
    dx_Cv_obc, & ! The unblocked lengths of the v-faces of the h-cell in m for OBC.
    IareaCv, &   ! The masked inverse areas of v-grid cells in m2.
    areaCv       ! The areas of the v-grid cells in m2.

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    mask2dBu, &  ! 0 for boundary points and 1 for ocean points on the q grid.  Nondim.
    geoLatBu, &  ! The geographic latitude at q points in degrees of latitude or m.
    geoLonBu, &  ! The geographic longitude at q points in degrees of longitude or m.
    dxBu, IdxBu, & ! dxBu is delta x at q points, in m, and IdxBu is 1/dxBu in m-1.
    dyBu, IdyBu, & ! dyBu is delta y at q points, in m, and IdyBu is 1/dyBu in m-1.
    areaBu, &    ! areaBu is the area of a q-cell, in m2
    IareaBu      ! IareaBu = 1/areaBu in m-2.

  real, pointer, dimension(:) :: &
    gridLatT => NULL(), gridLatB => NULL() ! The latitude of T or B points for
                        ! the purpose of labeling the output axes.
                        ! On many grids these are the same as geoLatT & geoLatBu.
  real, pointer, dimension(:) :: &
    gridLonT => NULL(), gridLonB => NULL() ! The longitude of T or B points for
                        ! the purpose of labeling the output axes.
                        ! On many grids these are the same as geoLonT & geoLonBu.
  character(len=40) :: &
    x_axis_units, &     !   The units that are used in labeling the coordinate
    y_axis_units        ! axes.

!  character(len=40) :: axis_units = ' '! Units for the horizontal coordinates.

  real :: g_Earth !   The gravitational acceleration in m s-2.
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    bathyT        ! Ocean bottom depth at tracer points, in m.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    CoriolisBu    ! The Coriolis parameter at corner points, in s-1.

end type SIS_hor_grid_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_hor_grid initializes the sea ice grid parameters.
subroutine set_hor_grid(G, param_file)
  type(SIS_hor_grid_type), intent(inout) :: G
  type(param_file_type)  , intent(in)    :: param_file
!   This subroutine sets up the necessary domain types and the sea-ice grid.

! Arguments: G - The sea-ice's horizontal grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

! This include declares and sets the variable "version".
#include "version_variable.h"

  real, allocatable, dimension(:,:)   :: depth, tmpx, tmpy, tmp_2d
  real, allocatable, dimension(:) :: xb1d, yb1d ! 1d global grid for diag_mgr
  real    :: angle, lon_scale

  integer, allocatable, dimension(:)  :: pelist, islist, ielist, jslist, jelist
  integer :: i, j, m, pe, ntiles, ncontacts
  integer :: isg, ieg, jsg, jeg
  integer :: is, ie, js, je, i_off, j_off
  integer :: ni, nj, dims(4)
  integer :: isca, ieca, jsca, jeca, isda, ieda, jsda, jeda
  integer :: npes

  logical :: symmetric       ! If true, use symmetric memory allocation.
  logical :: global_indexing
  character(len=256) :: grid_file, ocean_topog
  character(len=256) :: ocean_hgrid, ocean_mosaic
  character(len=200) :: mesg
  character(len=200) :: filename, topo_file, inputdir ! Strings for file/path
  character(len=200) :: topo_varname                  ! Variable name in file
  character(len=40)  :: mod_nm  = "hor_grid" ! This module's name.
  type(domain2d)     :: domain2
  type(domain2d), pointer :: Domain => NULL()

  grid_file = 'INPUT/grid_spec.nc'

  ! Set up the MOM_domain_type.  This will later occur via a call to MOM_domains_init.
  ! call MOM_domains_init(G%Domain, param_file, 1, dynamic=.true.)
  if (.not.associated(G%Domain)) then
    allocate(G%Domain)
    allocate(G%Domain%mpp_domain)
    allocate(G%Domain_aux)
    allocate(G%Domain_aux%mpp_domain)
  endif

!  pe = PE_here()
  npes = num_PEs()

#ifdef SYMMETRIC_MEMORY_
  symmetric = .true.
#else
  symmetric = .false.
#endif
#ifdef STATIC_MEMORY_
  call MOM_domains_init(G%domain, param_file, symmetric=symmetric, &
            static_memory=.true., NIHALO=NIHALO_, NJHALO=NJHALO_, &
            NIGLOBAL=NIGLOBAL_, NJGLOBAL=NJGLOBAL_, NIPROC=NIPROC_, &
            NJPROC=NJPROC_, domain_name="ice model", include_name="SIS2_memory.h")
#else
  call MOM_domains_init(G%domain, param_file, symmetric=symmetric, &
           domain_name="ice model", include_name="SIS2_memory.h")
#endif
  call clone_MOM_domain(G%domain, G%domain_aux, symmetric=.false., &
                        domain_name="ice model aux")

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod_nm, version)
  call get_param(param_file, mod_nm, "GLOBAL_INDEXING", global_indexing, &
                 "If true, use a global lateral indexing convention, so \n"//&
                 "that corresponding points on different processors have \n"//&
                 "the same index. This does not work with static memory.", &
                 default=.false.)
#ifdef STATIC_MEMORY_
  if (global_indexing) cal SIS_error(FATAL, "set_hor_grid : "//&
       "GLOBAL_INDEXING can not be true with STATIC_MEMORY.")
#endif

  call get_param(param_file, mod_nm, "FIRST_DIRECTION", G%first_direction, &
                 "An integer that indicates which direction goes first \n"//&
                 "in parts of the code that use directionally split \n"//&
                 "updates, with even numbers (or 0) used for x- first \n"//&
                 "and odd numbers used for y-first.", default=0)

  !--- first determine the if the grid file is using the correct format
  if (.not.(field_exist(grid_file, 'ocn_mosaic_file') .or. &
            field_exist(grid_file, 'gridfiles')) ) call SIS_error(FATAL, &
    'Error from ice_grid_mod(set_hor_grid): '//&
    'ocn_mosaic_file or gridfiles does not exist in file ' //trim(grid_file)//&
    '\nSIS2 only works with a mosaic format grid file.')

  call SIS_mesg("   Note from ice_grid_mod(set_hor_grid): "//&
                 "read grid from mosaic version grid", 5)

  if( field_exist(grid_file, "ocn_mosaic_file") ) then ! coupler mosaic
    call read_data(grid_file, "ocn_mosaic_file", ocean_mosaic)
    ocean_mosaic = "INPUT/"//trim(ocean_mosaic)
  else
    ocean_mosaic = trim(grid_file)
  end if
  ntiles = get_mosaic_ntiles(ocean_mosaic)
  if (ntiles /= 1) call SIS_error(FATAL, "Error from ice_grid_mod(set_hor_grid): "//&
      "ntiles should be 1 for ocean mosaic.")
  call read_data(ocean_mosaic, "gridfiles", ocean_hgrid)
  ocean_hgrid = 'INPUT/'//trim(ocean_hgrid)

  ! This code should be moved to MOM_domains_init once we start using a cubed-sphere grid.
  ! if (field_exist(ocean_mosaic, "contacts") ) then
  !   ncontacts = get_mosaic_ncontacts(ocean_mosaic)
  !   if (ncontacts < 1) call SIS_error(FATAL,'==>Error from ice_grid_mod(set_hor_grid): '//&
  !        'number of contacts should be larger than 0 when field contacts exist in file '//&
  !        trim(ocean_mosaic) )
  !   if (ncontacts > 2) call SIS_error(FATAL,'==>Error from ice_grid_mod(set_hor_grid): '//&
  !        'number of contacts should be no larger than 2')
  !   call get_mosaic_contact( ocean_mosaic, tile1(1:ncontacts), tile2(1:ncontacts),           &
  !        istart1(1:ncontacts), iend1(1:ncontacts), jstart1(1:ncontacts), jend1(1:ncontacts), &
  !        istart2(1:ncontacts), iend2(1:ncontacts), jstart2(1:ncontacts), jend2(1:ncontacts)  )
  !   do m = 1, ncontacts
  !     if (istart1(m) == iend1(m) ) then  ! x-direction contact, only cyclic condition
  !       if (istart2(m) /= iend2(m) ) call SIS_error(FATAL,  &
  !            "==>Error from ice_grid_mod(set_hor_grid): only cyclic condition is allowed for x-boundary")
  !       x_cyclic = .true.
  !     elseif ( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic or folded-north
  !       if ( jstart1(m) == jstart2(m) ) then ! folded north
  !          tripolar_grid=.true.
  !       else
  !          call SIS_error(FATAL, "==>Error from ice_grid_mod(set_hor_grid): "//&
  !            "only folded-north condition is allowed for y-boundary")
  !       endif
  !     else
  !       call SIS_error(FATAL,  &
  !            "==>Error from ice_grid_mod(set_hor_grid): invalid boundary contact")
  !     endif
  !   enddo
  ! endif

  !--- get grid size from the input file hgrid file.
  call field_size(ocean_hgrid, 'x', dims)
  if(mod(dims(1),2) /= 1) call SIS_error(FATAL, '==>Error from ice_grid_mod(set_hor_grid): '//&
      'x-size of x in file '//trim(ocean_hgrid)//' should be 2*niglobal+1')
  if(mod(dims(2),2) /= 1) call SIS_error(FATAL, '==>Error from ice_grid_mod(set_hor_grid): '//&
      'y-size of x in file '//trim(ocean_hgrid)//' should be 2*njglobal+1')
  ni = dims(1)/2
  nj = dims(2)/2

  if (ni /= G%Domain%niglobal) call SIS_error(FATAL, "set_hor_grid: "//&
    "The total i-grid size from file "//trim(ocean_hgrid)//" is inconsistent with SIS_input.")
  if (nj /= G%Domain%njglobal) call SIS_error(FATAL, "set_hor_grid: "//&
    "The total j-grid size from file "//trim(ocean_hgrid)//" is inconsistent with SIS_input.")

  call mpp_get_compute_domain(G%Domain%mpp_domain, isca, ieca, jsca, jeca )
  call mpp_get_data_domain(G%Domain%mpp_domain, isda, ieda, jsda, jeda )
  call mpp_get_global_domain(G%Domain%mpp_domain, isg, ieg, jsg, jeg )

  ! Allocate and fill in default values for elements of the sea ice grid type.
  if (global_indexing) then
    i_off = 0 ; j_off = 0
  else
    i_off = isda-1 ; j_off = jsda-1
    ! i_off = 1000 ; j_off = 1000 ! Use this for debugging.
  endif

  G%isc = isca-i_off ; G%iec = ieca-i_off ; G%jsc = jsca-j_off ; G%jec = jeca-j_off
  G%isd = isda-i_off ; G%ied = ieda-i_off ; G%jsd = jsda-j_off ; G%jed = jeda-j_off
  G%isg = isg ; G%ieg = ieg ; G%jsg = jsg ; G%jeg = jeg
!  G%ks = 0 ; G%ke = 0  ! Change this for shared ocean / ice grids.

  G%symmetric = G%Domain%symmetric
  G%nonblocking_updates = G%Domain%nonblocking_updates

  G%IscB = G%isc ; G%JscB = G%jsc
  G%IsdB = G%isd ; G%JsdB = G%jsd
  G%IsgB = G%isg ; G%JsgB = G%jsg
  if (G%symmetric) then
    G%IscB = G%isc-1 ; G%JscB = G%jsc-1
    G%IsdB = G%isd-1 ; G%JsdB = G%jsd-1
    G%IsgB = G%isg-1 ; G%JsgB = G%jsg-1
  endif
  G%IecB = G%iec ; G%JecB = G%jec
  G%IedB = G%ied ; G%JedB = G%jed
  G%IegB = G%ieg ; G%JegB = G%jeg

  i_off = isca - G%isc ; j_off = jsca - G%jsc

  call allocate_metrics(G)

  G%g_Earth = grav

  !--- z1l: loop through the pelist to find the symmetry processor.
  !--- This is needed to address the possibility that some of the all-land processor
  !--- regions are masked out. This is only needed for tripolar grid.
  if (G%Domain%Y_flags == FOLD_NORTH_EDGE) then
    allocate(pelist(npes), islist(npes), ielist(npes), jslist(npes), jelist(npes))
    call mpp_get_pelist(G%Domain%mpp_domain, pelist)
    call mpp_get_compute_domains(G%Domain%mpp_domain, &
             xbegin=islist, xend=ielist, ybegin=jslist, yend=jelist)

    do pe=1,npes ; if ( jslist(pe) == jsca .and. islist(pe) + ieca == ni+1 ) then
      if ( jelist(pe) /= jeca ) call SIS_error(FATAL, &
              "ice_model: jelist(p) /= jec but jslist(p) == jsc")
      if ( ielist(pe) + isca /= ni+1) call SIS_error(FATAL, &
              "ice_model: ielist(p) + isc /= ni+1 but islist(p) + iec == ni+1")
      exit
    endif ; enddo
    deallocate(pelist, islist, ielist, jslist, jelist)
  endif

end subroutine set_hor_grid


function Adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val
  real :: I_val
  ! This function implements Adcroft's rule for division by 0.

  I_val = 0.0
  if (val /= 0.0) I_val = 1.0/val
end function Adcroft_reciprocal

!---------------------------------------------------------------------

!> Returns true if the coordinates (x,y) are within the h-cell (i,j)
logical function isPointInCell(G, i, j, x, y)
  type(SIS_hor_grid_type),   intent(in) :: G    !< Grid type
  integer,                   intent(in) :: i, j !< i,j indices of cell to test
  real,                      intent(in) :: x, y !< x,y coordinates of point
! This is a crude calculation that assume a geographic coordinate system
  real :: xNE, xNW, xSE, xSW, yNE, yNW, ySE, ySW
  real :: p0, p1, p2, p3, l0, l1, l2, l3
  isPointInCell = .false.
  xNE = G%geoLonBu(i  ,j  ) ; yNE = G%geoLatBu(i  ,j  )
  xNW = G%geoLonBu(i-1,j  ) ; yNW = G%geoLatBu(i-1,j  )
  xSE = G%geoLonBu(i  ,j-1) ; ySE = G%geoLatBu(i  ,j-1)
  xSW = G%geoLonBu(i-1,j-1) ; ySW = G%geoLatBu(i-1,j-1)
  if (x<min(xNE,xNW,xSE,xSW) .or. x>max(xNE,xNW,xSE,xSW) .or. &
      y<min(yNE,yNW,ySE,ySW) .or. y>max(yNE,yNW,ySE,ySW) ) then
    return ! Avoid the more complicated calculation
  endif
  l0 = (x-xSW)*(ySE-ySW) - (y-ySW)*(xSE-xSW)
  l1 = (x-xSE)*(yNE-ySE) - (y-ySE)*(xNE-xSE)
  l2 = (x-xNE)*(yNW-yNE) - (y-yNE)*(xNW-xNE)
  l3 = (x-xNW)*(ySW-yNW) - (y-yNW)*(xSW-xNW)

  p0 = sign(1., l0) ; if (l0 == 0.) p0=0.
  p1 = sign(1., l1) ; if (l1 == 0.) p1=0.
  p2 = sign(1., l2) ; if (l2 == 0.) p2=0.
  p3 = sign(1., l3) ; if (l3 == 0.) p3=0.

  if ( (abs(p0)+abs(p2)) + (abs(p1)+abs(p3)) == abs((p0+p2) + (p1+p3)) ) then
    isPointInCell=.true.
  endif
end function isPointInCell

subroutine set_first_direction(G, y_first)
  type(SIS_hor_grid_type), intent(inout) :: G
  integer,               intent(in) :: y_first

  G%first_direction = y_first
end subroutine set_first_direction

!---------------------------------------------------------------------

subroutine allocate_metrics(G)
  type(SIS_hor_grid_type), intent(inout) :: G
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isg, ieg, jsg, jeg

  ! This subroutine allocates the lateral elements of the SIS_hor_grid_type that
  ! are always used and zeros them out.

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg

  ALLOC_(G%dxT(isd:ied,jsd:jed))       ; G%dxT(:,:) = 0.0
  ALLOC_(G%dxCu(IsdB:IedB,jsd:jed))    ; G%dxCu(:,:) = 0.0
  ALLOC_(G%dxCv(isd:ied,JsdB:JedB))    ; G%dxCv(:,:) = 0.0
  ALLOC_(G%dxBu(IsdB:IedB,JsdB:JedB))  ; G%dxBu(:,:) = 0.0
  ALLOC_(G%IdxT(isd:ied,jsd:jed))      ; G%IdxT(:,:) = 0.0
  ALLOC_(G%IdxCu(IsdB:IedB,jsd:jed))   ; G%IdxCu(:,:) = 0.0
  ALLOC_(G%IdxCv(isd:ied,JsdB:JedB))   ; G%IdxCv(:,:) = 0.0
  ALLOC_(G%IdxBu(IsdB:IedB,JsdB:JedB)) ; G%IdxBu(:,:) = 0.0

  ALLOC_(G%dyT(isd:ied,jsd:jed))       ; G%dyT(:,:) = 0.0
  ALLOC_(G%dyCu(IsdB:IedB,jsd:jed))    ; G%dyCu(:,:) = 0.0
  ALLOC_(G%dyCv(isd:ied,JsdB:JedB))    ; G%dyCv(:,:) = 0.0
  ALLOC_(G%dyBu(IsdB:IedB,JsdB:JedB))  ; G%dyBu(:,:) = 0.0
  ALLOC_(G%IdyT(isd:ied,jsd:jed))      ; G%IdyT(:,:) = 0.0
  ALLOC_(G%IdyCu(IsdB:IedB,jsd:jed))   ; G%IdyCu(:,:) = 0.0
  ALLOC_(G%IdyCv(isd:ied,JsdB:JedB))   ; G%IdyCv(:,:) = 0.0
  ALLOC_(G%IdyBu(IsdB:IedB,JsdB:JedB)) ; G%IdyBu(:,:) = 0.0

  ALLOC_(G%areaT(isd:ied,jsd:jed))       ; G%areaT(:,:) = 0.0
  ALLOC_(G%IareaT(isd:ied,jsd:jed))      ; G%IareaT(:,:) = 0.0
  ALLOC_(G%areaBu(IsdB:IedB,JsdB:JedB))  ; G%areaBu(:,:) = 0.0
  ALLOC_(G%IareaBu(IsdB:IedB,JsdB:JedB)) ; G%IareaBu(:,:) = 0.0

  ALLOC_(G%mask2dT(isd:ied,jsd:jed))      ; G%mask2dT(:,:) = 0.0
  ALLOC_(G%mask2dCu(IsdB:IedB,jsd:jed))   ; G%mask2dCu(:,:) = 0.0
  ALLOC_(G%mask2dCv(isd:ied,JsdB:JedB))   ; G%mask2dCv(:,:) = 0.0
  ALLOC_(G%mask2dBu(IsdB:IedB,JsdB:JedB)) ; G%mask2dBu(:,:) = 0.0
  ALLOC_(G%geoLatT(isd:ied,jsd:jed))      ; G%geoLatT(:,:) = 0.0
  ALLOC_(G%geoLatCu(IsdB:IedB,jsd:jed))   ; G%geoLatCu(:,:) = 0.0
  ALLOC_(G%geoLatCv(isd:ied,JsdB:JedB))   ; G%geoLatCv(:,:) = 0.0
  ALLOC_(G%geoLatBu(IsdB:IedB,JsdB:JedB)) ; G%geoLatBu(:,:) = 0.0
  ALLOC_(G%geoLonT(isd:ied,jsd:jed))      ; G%geoLonT(:,:) = 0.0
  ALLOC_(G%geoLonCu(IsdB:IedB,jsd:jed))   ; G%geoLonCu(:,:) = 0.0
  ALLOC_(G%geoLonCv(isd:ied,JsdB:JedB))   ; G%geoLonCv(:,:) = 0.0
  ALLOC_(G%geoLonBu(IsdB:IedB,JsdB:JedB)) ; G%geoLonBu(:,:) = 0.0

  ALLOC_(G%dx_Cv(isd:ied,JsdB:JedB))     ; G%dx_Cv(:,:) = 0.0
  ALLOC_(G%dy_Cu(IsdB:IedB,jsd:jed))     ; G%dy_Cu(:,:) = 0.0
  ALLOC_(G%dx_Cv_obc(isd:ied,JsdB:JedB)) ; G%dx_Cv_obc(:,:) = 0.0
  ALLOC_(G%dy_Cu_obc(IsdB:IedB,jsd:jed)) ; G%dy_Cu_obc(:,:) = 0.0

  ALLOC_(G%areaCu(IsdB:IedB,jsd:jed))  ; G%areaCu(:,:) = 0.0
  ALLOC_(G%areaCv(isd:ied,JsdB:JedB))  ; G%areaCv(:,:) = 0.0
  ALLOC_(G%IareaCu(IsdB:IedB,jsd:jed)) ; G%IareaCu(:,:) = 0.0
  ALLOC_(G%IareaCv(isd:ied,JsdB:JedB)) ; G%IareaCv(:,:) = 0.0

  ALLOC_(G%bathyT(isd:ied, jsd:jed)) ; G%bathyT(:,:) = 0.0
  ALLOC_(G%CoriolisBu(IsdB:IedB, JsdB:JedB)) ; G%CoriolisBu(:,:) = 0.0

  allocate(G%sin_rot(isd:ied,jsd:jed)) ; G%sin_rot(:,:) = 0.0
  allocate(G%cos_rot(isd:ied,jsd:jed)) ; G%cos_rot(:,:) = 1.0

  allocate(G%gridLonT(isg:ieg))   ; G%gridLonT(:) = 0.0
  allocate(G%gridLonB(isg-1:ieg)) ; G%gridLonB(:) = 0.0
  allocate(G%gridLatT(jsg:jeg))   ; G%gridLatT(:) = 0.0
  allocate(G%gridLatB(jsg-1:jeg)) ; G%gridLatB(:) = 0.0

end subroutine allocate_metrics

!---------------------------------------------------------------------
!> Release memory used by the SIS_hor_grid_type and related structures.
subroutine SIS_hor_grid_end(G)
  type(SIS_hor_grid_type), intent(inout) :: G

  DEALLOC_(G%dxT)  ; DEALLOC_(G%dxCu)  ; DEALLOC_(G%dxCv)  ; DEALLOC_(G%dxBu)
  DEALLOC_(G%IdxT) ; DEALLOC_(G%IdxCu) ; DEALLOC_(G%IdxCv) ; DEALLOC_(G%IdxBu)

  DEALLOC_(G%dyT)  ; DEALLOC_(G%dyCu)  ; DEALLOC_(G%dyCv)  ; DEALLOC_(G%dyBu)
  DEALLOC_(G%IdyT) ; DEALLOC_(G%IdyCu) ; DEALLOC_(G%IdyCv) ; DEALLOC_(G%IdyBu)

  DEALLOC_(G%areaT)  ; DEALLOC_(G%IareaT)
  DEALLOC_(G%areaBu) ; DEALLOC_(G%IareaBu)
  DEALLOC_(G%areaCu) ; DEALLOC_(G%IareaCu)
  DEALLOC_(G%areaCv)  ; DEALLOC_(G%IareaCv)

  DEALLOC_(G%mask2dT)  ; DEALLOC_(G%mask2dCu)
  DEALLOC_(G%mask2dCv) ; DEALLOC_(G%mask2dBu)

  DEALLOC_(G%geoLatT)  ; DEALLOC_(G%geoLatCu)
  DEALLOC_(G%geoLatCv) ; DEALLOC_(G%geoLatBu)
  DEALLOC_(G%geoLonT)  ; DEALLOC_(G%geoLonCu)
  DEALLOC_(G%geoLonCv) ; DEALLOC_(G%geoLonBu)

  DEALLOC_(G%dx_Cv) ; DEALLOC_(G%dy_Cu)
  DEALLOC_(G%dx_Cv_obc) ; DEALLOC_(G%dy_Cu_obc)

  DEALLOC_(G%bathyT)  ; DEALLOC_(G%CoriolisBu)
  DEALLOC_(G%sin_rot) ; DEALLOC_(G%cos_rot)

  deallocate(G%gridLonT) ; deallocate(G%gridLatT)
  deallocate(G%gridLonB) ; deallocate(G%gridLatB)

  deallocate(G%Domain%mpp_domain)
  deallocate(G%Domain)

end subroutine SIS_hor_grid_end

end module SIS_hor_grid_mod
