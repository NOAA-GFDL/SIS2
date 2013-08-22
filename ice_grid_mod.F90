!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_grid_mod - sets up grid and processor domain - Michael.Winton@noaa.gov   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_grid_mod

  use constants_mod,   only: radius, omega, pi, grav
  use mpp_mod,         only: mpp_pe, mpp_npes, mpp_root_pe, mpp_chksum
  use mpp_mod,         only: mpp_sync_self, mpp_send, mpp_recv, stdout, EVENT_RECV, COMM_TAG_1, NULL_PE
  use mpp_domains_mod, only: mpp_define_domains, CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
  use mpp_domains_mod, only: mpp_update_domains, domain2D, mpp_global_field, YUPDATE, XUPDATE, CORNER
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain, mpp_set_domain_symmetry
  use mpp_domains_mod, only: mpp_define_io_domain, mpp_copy_domain, mpp_get_global_domain
  use mpp_domains_mod, only: mpp_set_global_domain, mpp_set_data_domain, mpp_set_compute_domain
  use mpp_domains_mod, only: mpp_deallocate_domain, mpp_get_pelist, mpp_get_compute_domains
  use mpp_domains_mod, only: SCALAR_PAIR, CGRID_NE, BGRID_NE
  use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
  use MOM_file_parser, only : get_param, log_version, param_file_type
  use MOM_domains,     only : SIS_domain_type=>MOM_domain_type, pass_var, pass_vector
  use fms_mod,         only: field_exist, field_size, read_data
  use fms_mod,         only: get_global_att_value, stderr
  use mosaic_mod,      only: get_mosaic_ntiles, get_mosaic_ncontacts
  use mosaic_mod,      only: calc_mosaic_grid_area, get_mosaic_contact
  use grid_mod,        only: get_grid_cell_vertices

  implicit none ; private
  include 'netcdf.inc'
#include <SIS2_memory.h>

  public :: set_ice_grid, ice_grid_end, g_sum, get_avg
  public :: ice_line

  public :: Domain, im, jm
  public :: xb1d, yb1d
  public :: cell_area
  public :: grid_x_t,grid_y_t
  public :: tripolar_grid, x_cyclic

  type(domain2D), target, save :: Domain

type, public :: sea_ice_grid_type
  type(SIS_domain_type), pointer :: Domain => NULL()
  type(SIS_domain_type), pointer :: Domain_aux => NULL()
  integer :: isc, iec, jsc, jec ! The range of the computational domain indicies
  integer :: isd, ied, jsd, jed ! and data domain indicies at tracer cell centers.
  integer :: isg, ieg, jsg, jeg ! The range of the global domain tracer cell indicies.
  integer :: IscB, IecB, JscB, JecB ! The range of the computational domain indicies
  integer :: IsdB, IedB, JsdB, JedB ! and data domain indicies at tracer cell vertices.
  integer :: IsgB, IegB, JsgB, JegB ! The range of the global domain vertex indicies.
  integer :: isd_global         ! The values of isd and jsd in the global
  integer :: jsd_global         ! (decomposition invariant) index space.
  integer :: ks, ke             ! The range of ocean layer's vertical indicies.
  integer :: CatIce             ! The number of sea ice categories.
  integer :: NkIce              ! The number of vertical partitions within the
                                ! sea ice.
  integer :: NkSnow             ! The number of vertical partitions within the
                                ! snow atop the sea ice.
  
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
  logical ALLOCABLE_, dimension(NIMEM_,NJMEM_) :: &
    Lmask2dT     ! .true. for ocean points, .false. for land on the h-grid.

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
  logical ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEM_) :: &
    Lmask2dCu    ! .true. for ocean points, .false. for land on the v-grid.

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
  logical ALLOCABLE_, dimension(NIMEM_,NJMEMB_PTR_) :: &
    Lmask2dCv    ! .true. for ocean points, .false. for land on the v-grid.

  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    mask2dBu, &  ! 0 for boundary points and 1 for ocean points on the q grid.  Nondim.
    geoLatBu, &  ! The geographic latitude at q points in degrees of latitude or m.
    geoLonBu, &  ! The geographic longitude at q points in degrees of longitude or m.
    dxBu, IdxBu, & ! dxBu is delta x at q points, in m, and IdxBu is 1/dxBu in m-1.
    dyBu, IdyBu, & ! dyBu is delta y at q points, in m, and IdyBu is 1/dyBu in m-1.
    areaBu, &    ! areaBu is the area of a q-cell, in m2
    IareaBu      ! IareaBu = 1/areaBu in m-2.
  logical ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    Lmask2dBu    ! .true. for ocean points, .false. for land on the q-grid.

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

  character(len=40) :: axis_units = ' '! Units for the horizontal coordinates.

  real :: g_Earth !   The gravitational acceleration in m s-2.
  real ALLOCABLE_, dimension(NIMEMB_PTR_,NJMEMB_PTR_) :: &
    CoriolisBu    ! The Coriolis parameter at corner points, in s-1.

  ! The following are axis types defined for output.
  integer, dimension(3) :: axesBL, axesTL, axesCuL, axesCvL
  integer, dimension(3) :: axesBi, axesTi, axesCui, axesCvi
  integer, dimension(2) :: axesB1, axesT1, axesCu1, axesCv1
  integer, dimension(1) :: axeszi, axeszL

end type sea_ice_grid_type

type, public :: SIS2_domain_type
  type(domain2D), pointer :: mpp_domain => NULL() ! The domain with halos on
                                        ! this processor, centered at h points.
  integer :: niglobal, njglobal         ! The total horizontal domain sizes.
  integer :: nihalo, njhalo             ! The X- and Y- halo sizes in memory.
  logical :: symmetric                  ! True if symmetric memory is used with
                                        ! this domain.
  logical :: nonblocking_updates        ! If true, non-blocking halo updates are
                                        ! allowed.  The default is .false. (for now).
  integer :: layout(2), io_layout(2)    ! Saved data for sake of constructing
  integer :: X_FLAGS, Y_FLAGS           ! new domains of different resolution.
  logical :: use_io_layout              ! True if an I/O layout is available.
  logical, pointer :: maskmap(:,:)=> NULL() !option to mpp_define_domains
end type SIS2_domain_type


  integer                           :: im, jm ! , km         ! global domain and vertical size
  !
  ! grid geometry
  !
  logical                           ::  x_cyclic           ! x boundary condition
  logical                           ::  tripolar_grid      ! y boundary condition
  real, allocatable, dimension(:  ) ::  xb1d, yb1d         ! 1d global grid for diag_mgr
  real, allocatable, dimension(:  ) ::  grid_x_t,grid_y_t  ! 1d global grid for diag_mgr
  real, allocatable, dimension(:,:) ::  cell_area          ! grid cell area; sphere frac.
  
  integer            :: comm_pe                      ! pe to be communicated with

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_avg - take area weighted average over all partitions                     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine get_avg(x, cn, avg, wtd)
  real, dimension(:,:,:), intent(in)  :: x
  real, dimension(:,:,:), intent(in)  :: cn
  real, dimension(:,:),   intent(out) :: avg
  logical,      optional, intent(in)  :: wtd

  real, dimension(size(x,1),size(x,2)) :: wts
  logical :: do_wt
  integer :: i, j, k, ni, nj, nk

  do_wt = .false. ; if (present(wtd)) do_wt = wtd

  ni = size(x,1) ; nj = size(x,2); nk = size(x,3)
  if ((size(cn,1) /= ni) .or. (size(cn,2) /= nj) .or. (size(cn,3) /= nk)) &
    call SIS_error(FATAL, "Mismatched i- or j- sizes of x and cn in get_avg.")
  if ((size(avg,1) /= ni) .or. (size(avg,2) /= nj)) &
    call SIS_error(FATAL, "Mismatched i- or j- sizes of x and avg in get_avg.")
  if (size(cn,3) /= nk) &
    call SIS_error(FATAL, "Mismatched category sizes of x and cn in get_avg.")

  if (do_wt) then
    avg(:,:) = 0.0 ; wts(:,:) = 0.0
    do k=1,nk ; do j=1,nj ; do i=1,ni
      avg(i,j) = avg(i,j) + cn(i,j,k)*x(i,j,k)
      wts(i,j) = wts(i,j) + cn(i,j,k)
    enddo ; enddo ; enddo
     do j=1,nj ; do i=1,ni
      if (wts(i,j) > 0.) then
        avg(i,j) = avg(i,j) / wts(i,j)
      else
        avg(i,j) = 0.0
      endif
    enddo ; enddo
  else
    avg(:,:) = 0.0
    do k=1,nk ; do j=1,nj ; do i=1,ni
      avg(i,j) = avg(i,j) + cn(i,j,k)*x(i,j,k)
    enddo ; enddo ; enddo
  endif

end subroutine get_avg


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! g_sum - returns the global sum of a real array                               !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  real function g_sum(x)
    real, dimension(:,:) :: x

    real, dimension(1 :im, 1 :jm) :: g_x

    call mpp_global_field(Domain, x, g_x)
    g_sum = sum(g_x)

    return
  end function g_sum


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ice_grid - initialize sea ice grid for dynamics and transport            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ice_grid(G, param_file, ice_domain, km_in, layout, io_layout, maskmap )
  type(sea_ice_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file

  type(domain2D),        intent(inout) :: ice_domain
  integer,               intent(in)    :: km_in
  integer, dimension(2), intent(inout) :: layout
  integer, dimension(2), intent(inout) :: io_layout
  logical, optional,     intent(in)    :: maskmap(:,:)
! Arguments: G - The sea-ice's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
! This include declares and sets the variable "version".
#include "version_variable.h"

    real                                :: angle, lon_scale
    integer                             :: i, j, m, ntiles, ncontacts
    integer                             :: dims(4)
    real, allocatable, dimension(:,:,:) :: x_vert_t, y_vert_t
    real, allocatable,   dimension(:,:) :: depth, tmpx, tmpy, tmp_2d
    integer, dimension(2)               :: tile1, tile2
    integer, dimension(2)               :: istart1, iend1, jstart1, jend1
    integer, dimension(2)               :: istart2, iend2, jstart2, jend2
  integer :: x_flags, y_flags
    character(len=128)                  :: grid_file, ocean_topog
    character(len=256)                  :: ocean_hgrid, ocean_mosaic, attvalue
    type(domain2d)                      :: domain2
    integer                             :: isg, ieg, jsg, jeg
    integer                             :: is,  ie,  js,  je
    integer                             :: ni,nj,siz(4)
    integer, allocatable, dimension(:)  :: pelist, islist, ielist, jslist, jelist
    integer                             :: npes, p
    logical                             :: symmetrize, ndivx_is_even, im_is_even
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, km

  grid_file = 'INPUT/grid_spec.nc'
  ocean_topog = 'INPUT/topog.nc'

  !--- first determine the if the grid file is using the correct format
  if (.not.(field_exist(grid_file, 'ocn_mosaic_file') .or. &
            field_exist(grid_file, 'gridfiles')) ) call SIS_error(FATAL, &
    'Error from ice_grid_mod(set_ice_grid): '//&
    'ocn_mosaic_file or gridfiles does not exist in file ' //trim(grid_file)//&
    '\nSIS2 only works with a mosaic format grid file.')

  call SIS_mesg("   Note from ice_grid_mod(set_ice_grid): "//&
                 "read grid from mosaic version grid", 5)

  if( field_exist(grid_file, "ocn_mosaic_file") ) then ! coupler mosaic
    call read_data(grid_file, "ocn_mosaic_file", ocean_mosaic)
    ocean_mosaic = "INPUT/"//trim(ocean_mosaic)
  else
    ocean_mosaic = trim(grid_file)
  end if
  ntiles = get_mosaic_ntiles(ocean_mosaic)
  if (ntiles /= 1) call SIS_error(FATAL, "Error from ice_grid_mod(set_ice_grid): "//&
      "ntiles should be 1 for ocean mosaic.")
  call read_data(ocean_mosaic, "gridfiles", ocean_hgrid)
  ocean_hgrid = 'INPUT/'//trim(ocean_hgrid)
  call field_size(ocean_hgrid, 'x', siz)
  ni = siz(1)/2
  nj = siz(2)/2

  x_cyclic = .false.; tripolar_grid = .false.

  if (field_exist(ocean_mosaic, "contacts") ) then
    ncontacts = get_mosaic_ncontacts(ocean_mosaic)
    if (ncontacts < 1) call SIS_error(FATAL,'==>Error from ice_grid_mod(set_ice_grid): '//&
         'number of contacts should be larger than 0 when field contacts exist in file '//&
         trim(ocean_mosaic) )
    if (ncontacts > 2) call SIS_error(FATAL,'==>Error from ice_grid_mod(set_ice_grid): '//&
         'number of contacts should be no larger than 2')
    call get_mosaic_contact( ocean_mosaic, tile1(1:ncontacts), tile2(1:ncontacts),           &
         istart1(1:ncontacts), iend1(1:ncontacts), jstart1(1:ncontacts), jend1(1:ncontacts), &
         istart2(1:ncontacts), iend2(1:ncontacts), jstart2(1:ncontacts), jend2(1:ncontacts)  )
    do m = 1, ncontacts
      if (istart1(m) == iend1(m) ) then  ! x-direction contact, only cyclic condition
        if (istart2(m) /= iend2(m) ) call SIS_error(FATAL,  &
             "==>Error from ice_grid_mod(set_ice_grid): only cyclic condition is allowed for x-boundary")
        x_cyclic = .true.
      elseif ( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic or folded-north
        if ( jstart1(m) == jstart2(m) ) then ! folded north
           tripolar_grid=.true.
        else 
           call SIS_error(FATAL, "==>Error from ice_grid_mod(set_ice_grid): "//&
             "only folded-north condition is allowed for y-boundary")
        endif
      else 
        call SIS_error(FATAL,  &
             "==>Error from ice_grid_mod(set_ice_grid): invalid boundary contact")
      endif
    enddo
  endif

  !--- get grid size
  call field_size(ocean_hgrid, 'x', dims)
  if(mod(dims(1),2) .NE. 1) call SIS_error(FATAL, '==>Error from ice_grid_mod(set_ice_grid): '//&
      'x-size of x in file '//trim(ocean_hgrid)//' should be 2*ni+1')
  if(mod(dims(2),2) .NE. 1) call SIS_error(FATAL, '==>Error from ice_grid_mod(set_ice_grid): '//&
      'y-size of x in file '//trim(ocean_hgrid)//' should be 2*nj+1')
  im = dims(1)/2 
  jm = dims(2)/2 

  if (x_cyclic) then
    call SIS_mesg("==>Note from ice_grid_mod: x_boundary_type is cyclic")
  else
    call SIS_mesg("==>Note from ice_grid_mod: x_boundary_type is solid_walls")
  endif
  if (tripolar_grid) then
    call SIS_mesg("==>Note from ice_grid_mod: y_boundary_type is fold_north_edge")
  else
    call SIS_mesg("==>Note from ice_grid_mod: y_boundary_type is solid_walls")
  endif

  ! default is merdional domain decomp. to load balance xgrid
  if( layout(1)==0 .and. layout(2)==0 ) layout=(/ mpp_npes(), 1 /)
  if( layout(1)/=0 .and. layout(2)==0 ) layout(2) = mpp_npes()/layout(1)
  if( layout(1)==0 .and. layout(2)/=0 ) layout(1) = mpp_npes()/layout(2)
  x_flags = 0 ; y_flags = 0
  if (tripolar_grid) then    
    !z1l: Tripolar grid requires symmetry in i-direction domain decomposition
    ndivx_is_even = (mod(layout(1),2) == 0)
    im_is_even    = (mod(im,2) == 0)
    symmetrize = ( ndivx_is_even .AND. im_is_even ) .OR. &
            (  (.NOT.ndivx_is_even) .AND.  (.NOT.im_is_even) ) .OR. &
            (  (.NOT.ndivx_is_even) .AND. im_is_even .AND. layout(1) .LT. im/2 )

    if( .not. symmetrize) then
      call SIS_error(FATAL, "ice_model(set_ice_grid): tripolar regrid requires symmetry in i-direction domain decomposition")
    endif
    x_flags=CYCLIC_GLOBAL_DOMAIN ; y_flags=FOLD_NORTH_EDGE
  else if(x_cyclic) then
    x_flags=CYCLIC_GLOBAL_DOMAIN ; y_flags = 0
  else
    x_flags = 0 ; y_flags = 0
  endif
  call mpp_define_domains( (/1,im,1,jm/), layout, Domain, maskmap=maskmap, &
                          xflags=x_flags, xhalo=1, yflags=y_flags, yhalo=1, &
                          name='ice model' )
  call mpp_define_io_domain(Domain, io_layout)

  call mpp_get_compute_domain( Domain, isc, iec, jsc, jec )
  call mpp_get_data_domain( Domain, isd, ied, jsd, jed )
  call mpp_get_global_domain( Domain, isg, ieg, jsg, jeg )

  call mpp_define_domains( (/1,im,1,jm/), layout, ice_domain, maskmap=maskmap, name='ice_nohalo') ! domain without halo
  call mpp_define_io_domain(ice_domain, io_layout)

  ! Set up the SIS_domain_type.  This will later occur via a call to MOM_domains_init.
  ! call MOM_domains_init(G%Domain, param_file, 1, dynamic=.true.)
  if (.not.associated(G%Domain)) allocate(G%Domain)
  G%Domain%mpp_domain => Domain
  G%Domain%niglobal = im ; G%Domain%njglobal = jm
  G%Domain%nihalo = 1 ; G%Domain%njhalo = 1
  G%Domain%symmetric = .false.
  G%Domain%layout(:) = layout(:) ; G%Domain%io_layout(:) = io_layout(:)
  G%Domain%X_FLAGS = X_FLAGS ; G%Domain%Y_FLAGS = Y_FLAGS
  G%Domain%nonblocking_updates = .false.
  ! call get_domain_extent(G%Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, &
  !                          isg, ieg, jsg, jeg, i_offset, j_offset, G%Domain%Symmetric, .false.)

  ! Allocate and fill in default values for elements of the sea ice grid type.
  G%isc = isc ; G%iec = iec ; G%jsc = jsc ; G%jec = jec
  G%isd = isd ; G%ied = ied ; G%jsd = jsd ; G%jed = jed
  G%isg = isg ; G%ieg = ieg ; G%jsg = jsg ; G%jeg = jeg

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

  call allocate_metrics(G)

  allocate(G%sin_rot(isd:ied,jsd:jed)) ; G%sin_rot(:,:) = 0.0
  allocate(G%cos_rot(isd:ied,jsd:jed)) ; G%cos_rot(:,:) = 1.0

  !--- read data from grid_spec.nc
  allocate(depth(isc:iec,jsc:jec))
  call read_data(ocean_topog, 'depth', depth(isc:iec,jsc:jec), Domain)
  do j=jsc,jec ; do i=isc,iec ; if (depth(i,j) > 0) then
    G%mask2dT(i,j) = 1.0
  endif ; enddo ; enddo
  deallocate(depth)
  call mpp_update_domains(G%mask2dT, Domain)

  do J=jsc-1,jec ; do I=isc-1,iec
    if( G%mask2dT(i,j)>0.5 .and. G%mask2dT(i,j+1)>0.5 .and. &
        G%mask2dT(i+1,j)>0.5 .and. G%mask2dT(i+1,j+1)>0.5 ) then
       G%mask2dBu(I,J) = 1.0 ; G%Lmask2dBu(I,J) = .true.
    else
       G%mask2dBu(I,J) = 0.0 ; G%Lmask2dBu(I,J) = .false.
    endif
  enddo ; enddo

  do j=jsc,jec ; do I=isc-1,iec
    if( G%mask2dT(i,j)>0.5 .and. G%mask2dT(i+1,j)>0.5 ) then
       G%mask2dCu(I,j) = 1.0 ; G%Lmask2dCu(I,j) = .true.
    else
       G%mask2dCu(I,j) = 0.0 ; G%Lmask2dCu(I,j) = .false.
    endif
  enddo ; enddo

  do J=jsc-1,jec ; do i=isc,iec
    if( G%mask2dT(i,j)>0.5 .and. G%mask2dT(i,j+1)>0.5 ) then
       G%mask2dCv(i,J) = 1.0 ; G%Lmask2dCv(i,J) = .true.
    else
       G%mask2dCv(i,J) = 0.0 ; G%Lmask2dCv(i,J) = .false.
    endif
  enddo ; enddo

  do j=jsd,jed ; do i=isd,ied
    G%Lmask2dT(i,j) = (G%mask2dT(i,j) > 0.5) 
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do i=isd,ied
    G%Lmask2dCv(i,J) = (G%mask2dCv(i,J) > 0.5) 
  enddo ; enddo
  do j=jsd,jed ; do I=G%IsdB,G%IedB
    G%Lmask2dCu(I,j) = (G%mask2dCu(I,j) > 0.5) 
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB
    G%Lmask2dBu(I,J) = (G%mask2dBu(I,J) > 0.5) 
  enddo ; enddo

    if(tripolar_grid) then
       if (jsc==1.and.any(G%mask2dT(:,jsc)>0.5)) call SIS_error(FATAL, &
          'ice_model_mod: ice model requires southernmost row of land', all_print=.true.);
    endif

#ifdef STATIC_MEMORY_
  call get_param(param_file, "SIS_grid", "NCAT_ICE", G%CatIce, &
                 "The number of sea ice thickness categories.", units="nondim", &
                 default=km_in-1)
  if (G%CatIce /= NCAT_ICE_) call MOM_error(FATAL, "MOM_grid_init: " // &
       "Mismatched number of layers NK_ICE between SIS_memory.h and param_file "//&
       "or the input namelist file.")
  call get_param(param_file, "SIS_grid", "NK_ICE", G%NkIce, &
                 "The number of layers within the sea ice.", units="nondim", &
                 default=NK_ICE_)
  if (G%NkIce /= NK_ICE_) call MOM_error(FATAL, "MOM_grid_init: " // &
       "Mismatched number of layers NK_ICE between SIS_memory.h and param_file")

  call get_param(param_file, "SIS_grid", "NK_SNOW", G%NkSnow, &
                 "The number of layers within the snow atop the sea ice.", &
                 units="nondim", default=NK_SNOW_)
  if (G%NkSnow /= NK_SNOW_) call MOM_error(FATAL, "MOM_grid_init: " // &
       "Mismatched number of layers NK_SNOW between SIS_memory.h and param_file")

#else
  call get_param(param_file, "SIS_grid", "NCAT_ICE", G%CatIce, &
                 "The number of sea ice thickness categories.", units="nondim", &
                 default=km_in-1)
  call get_param(param_file, "SIS_grid", "NK_ICE", G%NkIce, &
                 "The number of layers within the sea ice.", units="nondim", &
                 default=4) ! Valid for SIS5L; Perhaps this should be ..., fail_if_missing=.true.
  call get_param(param_file, "SIS_grid", "NK_SNOW", G%NkSnow, &
                 "The number of layers within the snow atop the sea ice.", &
                 units="nondim", default=1) ! Perhaps this should be ..., fail_if_missing=.true.
#endif

    km = km_in
    km = G%CatIce + 1  ! Add 1 for the ice-free category.

    allocate ( cell_area(isc:iec,jsc:jec) )
    
    allocate (grid_x_t(ni))
    allocate (grid_y_t(nj))    
    grid_x_t = 0.0;grid_y_t = 0.0 

  call mpp_copy_domain(domain, domain2)
  call mpp_set_compute_domain(domain2, 2*isc-1, 2*iec+1, 2*jsc-1, 2*jec+1, 2*(iec-isc)+3, 2*(jec-jsc)+3 )
  call mpp_set_data_domain   (domain2, 2*isd-1, 2*ied+1, 2*jsd-1, 2*jed+1, 2*(ied-isd)+3, 2*(jed-jsd)+3 )   
  call mpp_set_global_domain (domain2, 2*isg-1, 2*ieg+1, 2*jsg-1, 2*jeg+1, 2*(ieg-isg)+3, 2*(jeg-jsg)+3 )   
  call mpp_get_compute_domain(domain2, is, ie, js, je)
  if(is .NE. 2*isc-1 .OR. ie .NE. 2*iec+1 .OR. js .NE. 2*jsc-1 .OR. je .NE. 2*jec+1) then
    call SIS_error(FATAL, 'ice_grid_mod: supergrid domain is not set properly')
  endif
  allocate(tmpx(is:ie, js:je), tmpy(is:ie, js:je) )
  call read_data(ocean_hgrid, 'x', tmpx, domain2)
  call read_data(ocean_hgrid, 'y', tmpy, domain2)     
  do J=jsc-1,jec ; do I=isc-1,iec
    G%geoLonBu(I,J) = tmpx(2*i+1,2*j+1)
    G%geoLatBu(I,J) = tmpy(2*i+1,2*j+1)
  enddo ; enddo
  deallocate(tmpx, tmpy)
  call calc_mosaic_grid_area(G%geoLonBu(isc-1:iec,jsc-1:jec)*pi/180, &
                             G%geoLatBu(isc-1:iec,jsc-1:jec)*pi/180, G%areaT(isc:iec,jsc:jec))
  do j=jsc,jec ; do i=isc,iec
    cell_area(i,j) = G%mask2dT(i,j) * G%areaT(i,j)/(4*PI*RADIUS*RADIUS)
  enddo ; enddo       
  call mpp_deallocate_domain(domain2)


    ! Determine the domain-averaged latitudes and longitudes on the grid for
    ! diagnostic purposes. This is a hold-over and will be changed to match
    ! the MOM6 convention of using the "nominal" latitudes and longitudes,
    ! which are the actual lat & lon on spherical portions of the grid.
    allocate ( xb1d (im+1), yb1d (jm+1) )
    if(PRESENT(maskmap)) then    
       allocate(tmpx(im+1,jm+1), tmpy(im+1,jm+1))
       call get_grid_cell_vertices('OCN', 1, tmpx, tmpy)
       xb1d(:) = sum(tmpx,2)/(jm+1)
       yb1d(:) = sum(tmpy,1)/(im+1)
       deallocate(tmpx, tmpy)
    else
       allocate ( tmpx(isc:iec+1, jm+1) )
       call mpp_set_domain_symmetry(Domain, .TRUE.)
       call mpp_global_field(Domain, G%geoLonBu(isc-1:iec,jsc-1:jec), tmpx, flags=YUPDATE, position=CORNER)
       allocate ( tmp_2d(isc:iec+1, jsc:jec+1) )
       tmp_2d = 0
       tmp_2d(isc:iec+1,jsc) = sum(tmpx,2)/(jm+1);
       deallocate(tmpx)
       allocate ( tmpx(im+1, jsc:jec+1) )

       call mpp_global_field(Domain, G%geoLatBu(isc-1:iec,jsc-1:jec), tmpx, flags=XUPDATE, position=CORNER)
       xb1d(:) = tmpx(:,jsc)
       deallocate(tmpx, tmp_2d)

       allocate ( tmpy(im+1, jsc:jec+1) )
       call mpp_global_field(Domain, G%geoLatBu(isc-1:iec,jsc-1:jec), tmpy, flags=XUPDATE, position=CORNER)
       allocate ( tmp_2d(isc:iec+1, jsc:jec+1) )
       tmp_2d = 0
       tmp_2d(isc,jsc:jec+1) = sum(tmpy,1)/(im+1);
       deallocate(tmpy)
       allocate ( tmpy(isc:iec+1, jm+1) )
       call mpp_global_field(Domain, tmp_2d, tmpy, flags=YUPDATE, position=CORNER)
       yb1d(:) = tmpy(isc,:)
       deallocate(tmpy, tmp_2d)
       call mpp_set_domain_symmetry(Domain, .FALSE.)
    endif


  G%gridLatB(jsg-1:jeg) = yb1d(1:jm+1)
  G%gridLonB(isg-1:ieg) = xb1d(1:im+1)
  do i=isg,ieg ; G%gridLonT(i) = 0.5*(G%gridLonB(i-1)+G%gridLonB(i)) ; enddo
  do j=jsg,jeg ; G%gridLatT(j) = 0.5*(G%gridLatB(j-1)+G%gridLatB(j)) ; enddo

  do j=jsc,jec ; do I=isc-1,iec
    G%dyCu(I,j) = edge_length(G%geoLonBu(I,J-1),G%geoLatBu(I,J-1), &
                                 G%geoLonBu(I,J),G%geoLatBu(I,J))
  enddo ; enddo
  do J=jsc-1,jec ; do i=isc,iec
    G%dxCv(i,J) = edge_length(G%geoLonBu(I-1,J), G%geoLatBu(I-1,J), &
                                 G%geoLonBu(I,J), G%geoLatBu(I,J))
  enddo ; enddo
  do j=jsc,jec ; do i=isc,iec
    lon_scale    = cos((G%geoLatBu(I-1,J-1) + G%geoLatBu(I,J-1  ) + &
                        G%geoLatBu(I-1,J) + G%geoLatBu(I,J)) * atan(1.0)/180)
    angle        = atan2((G%geoLonBu(I-1,J) + G%geoLonBu(I,J) - &
                          G%geoLonBu(I-1,J-1) - G%geoLonBu(I,J-1))*lon_scale, &
                          G%geoLatBu(I-1,J) + G%geoLatBu(I,J) - &
                          G%geoLatBu(I-1,J-1) - G%geoLatBu(I,J-1) )
    G%sin_rot(i,j) = sin(angle) ! angle is the clockwise angle from lat/lon to ocean
    G%cos_rot(i,j) = cos(angle) ! grid (e.g. angle of ocean "north" from true north)
  enddo ; enddo

  call mpp_update_domains(G%dyCu, G%dxCv, Domain, gridtype=CGRID_NE, flags = SCALAR_PAIR)

  ! ### THIS DOESN'T SEEM RIGHT AT THE TRIPOLAR FOLD. -RWH
    call mpp_update_domains(G%cos_rot, Domain)
    call mpp_update_domains(G%sin_rot, Domain)

    do j = jsc, jec
       do i = isc, iec
          G%dxT(i,j) = (G%dxCv(i,J-1) + G%dxCv(I,j) )/2
          if (G%mask2dT(i,j) > 0.0) then
!             G%dyT(i,j) = G%areaT(i,j)/G%dxT(i,j)  !### ANSWERS CHANGE?
            G%dyT(i,j) = cell_area(i,j)*4*pi*radius*radius/G%dxT(i,j)
          else
            G%dyT(i,j) = (G%dyCu(I-1,j) + G%dyCu(I,j) )/2
          endif
       enddo
    enddo

  ! ### THIS SHOULD BE A SCALAR PAIR FOR CUBED SPHERE, ETC. -RWH
    call mpp_update_domains(G%dxT, Domain )
    call mpp_update_domains(G%dyT, Domain )

    G%dxBu(:,:) = 1.0 ; G%IdxBu(:,:) = 1.0
    G%dyBu(:,:) = 1.0 ; G%IdyBu(:,:) = 1.0

  !### ADD PARENTHESIS FOR REPRODUCIBILITY.
  do J=jsc-1,jec ; do I=isc-1,iec
    G%dxBu(I,J) = 0.25*(G%dxT(i+1,j+1)+G%dxT(i+1,j)+G%dxT(i,j+1)+G%dxT(i,j) )
    G%IdxBu(I,J) = 1.0 / G%dxBu(I,j)
    G%dyBu(I,J) = 0.25*(G%dyT(i+1,j+1)+G%dyT(i+1,j)+G%dyT(i,j+1)+G%dyT(i,j) )
    G%IdyBu(I,J) = 1.0 / G%dyBu(I,j)
  enddo ; enddo

  call mpp_update_domains(G%dxBu, G%dyBu, Domain, gridtype=BGRID_NE, flags=SCALAR_PAIR )
  call mpp_update_domains(G%IdxBu, G%IdyBu, Domain, gridtype=BGRID_NE, flags=SCALAR_PAIR )

  do j=jsc,jec ; do i=isc,iec
    !### REGROUP FOR ROTATIONAL REPRODUCIBILITY
    G%geoLonT(i,j) = lon_avg( (/ G%geoLonBu(I-1,J-1), G%geoLonBu(I,J-1), &
                                 G%geoLonBu(I-1,J), G%geoLonBu(I,J) /) )
    G%geoLatT(i,j) = (G%geoLatBu(I-1,J-1) + G%geoLatBu(I,J-1) + &
                      G%geoLatBu(I-1,J)   + G%geoLatBu(I,J)) / 4
  enddo ; enddo

  do J=jsc-1,jec ; do I=isc-1,iec
    G%CoriolisBu(I,J) = 2*omega*sin(G%geoLatBu(I,J)*pi/180)
  enddo ; enddo

  G%g_Earth = grav

    !--- z1l: loop through the pelist to find the symmetry processor.
    !--- This is needed to address the possibility that some of the all-land processor 
    !--- regions are masked out. This is only needed for tripolar grid.
    if(tripolar_grid) then
       npes = mpp_npes()
       allocate(pelist(npes), islist(npes), ielist(npes), jslist(npes), jelist(npes))
       call mpp_get_pelist(Domain, pelist)
       call mpp_get_compute_domains(Domain, xbegin=islist, xend=ielist, ybegin=jslist, yend=jelist)

       comm_pe = NULL_PE
 
       do p = 1, npes
          if( jslist(p) == jsc .AND. islist(p) + iec == im+1 ) then
             if( jelist(p) .NE. jec ) then
                call SIS_error(FATAL, "ice_model: jelist(p) .NE. jec but jslist(p) == jsc")
             endif
             if( ielist(p) + isc .NE. im+1) then
                call SIS_error(FATAL, "ice_model: ielist(p) + isc .NE. im+1 but islist(p) + iec == im+1")
             endif
             comm_pe = pelist(p)
             exit 
          endif
       enddo
       deallocate(pelist, islist, ielist, jslist, jelist)
    endif
!    comm_pe = mpp_pe() + layout(1) - 2*mod(mpp_pe()-mpp_root_pe(),layout(1)) - 1

  end subroutine set_ice_grid

!#####################################################################
!--- release memory
subroutine ice_grid_end(G)
  type(sea_ice_grid_type), intent(inout) :: G

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

  DEALLOC_(G%CoriolisBu)
  DEALLOC_(G%sin_rot) ; DEALLOC_(G%cos_rot)

  deallocate(G%gridLatT) ; deallocate(G%gridLatB)
  deallocate(G%gridLonT) ; deallocate(G%gridLonB)

  deallocate( xb1d, yb1d, cell_area )

end subroutine ice_grid_end

!#####################################################################

real function lon_avg(lons)
  real, dimension(:), intent(in) :: lons

  real, dimension(size(lons(:))) :: lons2 ! lons relative to lon(1)
  integer                        :: i

  lons2(1) = 0.0
  do i=2,size(lons(:))
    lons2(i) = lons(i)-lons(1)
    if (lons2(i) >  180) lons2(i) = lons2(i) - 360;
    if (lons2(i) < -180) lons2(i) = lons2(i) + 360;
  end do
  lon_avg = lons(1)+sum(lons2)/size(lons(:))
end function lon_avg

!#####################################################################
function edge_length(x1, y1, x2, y2)
  real, intent(in) :: x1, x2, y1, y2 ! end-point coordinates in degrees
  real             :: edge_length
  real             :: dx, dy

  dx = (x2-x1)*cos((atan(1.0)/45)*(y2+y1)/2)
  dy = y2-y1
  edge_length = radius*(atan(1.0)/45)*(dx*dx+dy*dy)**0.5
end function edge_length

!#####################################################################
subroutine ice_line(year, day, second, cn_ocn, sst, G)
  integer,                         intent(in) :: year, day, second
  type(sea_ice_grid_type),         intent(in) :: G
  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: cn_ocn
  real, dimension(G%isc:G%iec,G%jsc:G%jec),   intent(in) :: sst

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: x
  real :: gx(3)
  integer :: n, i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

!    if (.not.(second==0 .and. mod(day,5)==0) ) return  ! safe?

  do n=-1,1,2
    do j=jsc,jec ; do i=isc,iec
      x(i,j) = 0.0
      if (cn_ocn(i,j)<0.85 .and. n*G%geoLatT(i,j)>0.0) &
        x(i,j) = G%mask2dT(i,j)*G%areaT(i,j)
    enddo ; enddo
    gx((n+3)/2) = g_sum(x(isc:iec,jsc:jec))/1e12
  enddo
  gx(3) = g_sum(sst(isc:iec,jsc:jec)*G%mask2dT(isc:iec,jsc:jec)*G%areaT(isc:iec,jsc:jec)) / &
         (g_sum(G%mask2dT(isc:iec,jsc:jec)*G%areaT(isc:iec,jsc:jec)) + 1e-10)
  !
  ! print info every 5 days
  !
  if ( mpp_pe()==0 .and. second==0 .and. mod(day,5)==0 ) &
    print '(a,2I4,3F10.5)','ICE y/d (SH_ext NH_ext SST):', year, day, gx
end subroutine ice_line


subroutine allocate_metrics(G)
  type(sea_ice_grid_type), intent(inout) :: G
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isg, ieg, jsg, jeg

  ! This subroutine allocates the lateral elements of the sea_ice_grid_type that
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
  ALLOC_(G%Lmask2dT(isd:ied,jsd:jed))      ; G%Lmask2dT(:,:) = .false.
  ALLOC_(G%Lmask2dCu(IsdB:IedB,jsd:jed))   ; G%Lmask2dCu(:,:) = .false.
  ALLOC_(G%Lmask2dCv(isd:ied,JsdB:JedB))   ; G%Lmask2dCv(:,:) = .false.
  ALLOC_(G%Lmask2dBu(IsdB:IedB,JsdB:JedB)) ; G%Lmask2dBu(:,:) = .false.
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

  ALLOC_(G%CoriolisBu(IsdB:IedB, JsdB:JedB)) ; G%CoriolisBu(:,:) = 0.0

  allocate(G%gridLonT(isg:ieg))   ; G%gridLonT(:) = 0.0
  allocate(G%gridLonB(isg-1:ieg)) ; G%gridLonB(:) = 0.0
  allocate(G%gridLatT(jsg:jeg))   ; G%gridLatT(:) = 0.0
  allocate(G%gridLatB(jsg-1:jeg)) ; G%gridLatB(:) = 0.0

end subroutine allocate_metrics

end module ice_grid_mod
