!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_grid_mod - sets up grid, processor domain, and does advection-Michael.Winton@noaa.gov!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_grid_mod

  use constants_mod,   only: radius, omega, pi
  use mpp_mod,         only: mpp_pe, mpp_npes, mpp_root_pe, mpp_error, NOTE, mpp_chksum
  use mpp_mod,         only: mpp_sync_self, mpp_send, mpp_recv, stdout, EVENT_RECV, COMM_TAG_1, NULL_PE
  use mpp_domains_mod, only: mpp_define_domains, CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
  use mpp_domains_mod, only: mpp_update_domains, domain2D, mpp_global_field, YUPDATE, XUPDATE, CORNER
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain, mpp_set_domain_symmetry
  use mpp_domains_mod, only: mpp_define_io_domain, mpp_copy_domain, mpp_get_global_domain
  use mpp_domains_mod, only: mpp_set_global_domain, mpp_set_data_domain, mpp_set_compute_domain
  use mpp_domains_mod, only: mpp_deallocate_domain, mpp_get_pelist, mpp_get_compute_domains
  use mpp_domains_mod, only: SCALAR_PAIR, CGRID_NE, BGRID_NE
  use fms_mod,         only: error_mesg, FATAL, field_exist, field_size, read_data
  use fms_mod,         only: get_global_att_value, stderr
  use mosaic_mod,      only: get_mosaic_ntiles, get_mosaic_ncontacts
  use mosaic_mod,      only: calc_mosaic_grid_area, get_mosaic_contact
  use grid_mod,        only: get_grid_cell_vertices

  implicit none ; private
  include 'netcdf.inc'
#include <SIS2_memory.h>

  public :: set_ice_grid, dt_evp, evp_sub_steps, g_sum, ice_avg, all_avg
  public :: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, im, jm, km
  public :: dtw, dte, dts, dtn, dxt, dxv, dyt, dyv, cor, xb1d, yb1d
  public :: geo_lon, geo_lat, sin_rot, cos_rot, cell_area, wett, wetv
  public :: geo_lonv_ib, geo_latv_ib
  public :: grid_x_t,grid_y_t
  public :: dTdx, dTdy, t_on_uv, t_to_uv, uv_to_t, ice_advect
  public :: ice_line, vel_t_to_uv, cut_check, latitude, slab_ice_advect
  public :: dxdy, dydx, ice_grid_end
  public :: tripolar_grid, x_cyclic, dt_adv
  public :: reproduce_siena_201303

  type(domain2D), save :: Domain

type, public :: sea_ice_grid_type
  type(SIS2_domain_type), pointer :: Domain => NULL()
  type(SIS2_domain_type), pointer :: Domain_aux => NULL()
  integer :: isc, iec, jsc, jec ! The range of the computational domain indicies
  integer :: isd, ied, jsd, jed ! and data domain indicies at tracer cell centers.
  integer :: isg, ieg, jsg, jeg ! The range of the global domain tracer cell indicies.
  integer :: IscB, IecB, JscB, JecB ! The range of the computational domain indicies
  integer :: IsdB, IedB, JsdB, JedB ! and data domain indicies at tracer cell vertices.
  integer :: IsgB, IegB, JsgB, JegB ! The range of the global domain vertex indicies.
  integer :: isd_global         ! The values of isd and jsd in the global
  integer :: jsd_global         ! (decomposition invariant) index space.
  integer :: ks, ke             ! The range of layer's vertical indicies.
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
    IareaT       ! IareaT = 1/areaT, in m-2.

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
  integer :: X_FLAGS,Y_FLAGS            ! new domains of different resolution.
  logical :: use_io_layout              ! True if an I/O layout is available.
  logical, pointer :: maskmap(:,:)=> NULL() !option to mpp_define_domains
end type SIS2_domain_type

  integer                           :: isc, iec, jsc, jec ! compute domain
  integer                           :: isd, ied, jsd, jed ! data domain
  integer                           :: im, jm, km         ! global domain and vertical size
  real, allocatable, dimension(:,:) :: wett, wetv         ! t and v cell masks
  !
  ! grid geometry
  !
  logical                           ::  x_cyclic           ! x boundary condition
  logical                           ::  tripolar_grid      ! y boundary condition
  real, allocatable, dimension(:,:) ::  dtw, dte, dts, dtn ! size of t cell sides
  real, allocatable, dimension(:,:) ::  dxt, dxv           ! x-extent of t and v cells
  real, allocatable, dimension(:,:) ::  dyt, dyv           ! y-extent of t and v cells
  real, allocatable, dimension(:,:) ::  dxdy, dydx         
  real, allocatable, dimension(:,:) ::  latitude           ! latitude of t cells
  real, allocatable, dimension(:,:) ::  cor                ! coriolis on v cells
  real, allocatable, dimension(:,:) ::  geo_lat            ! true latitude (rotated grid)
  real, allocatable, dimension(:,:) ::  geo_lon            ! true longitude              
  real, allocatable, dimension(:,:) ::  geo_latv_ib        ! true latitude at the grid corners
  real, allocatable, dimension(:,:) ::  geo_lonv_ib        ! true longitude at the grid corners
  real, allocatable, dimension(:  ) ::  xb1d, yb1d         ! 1d global grid for diag_mgr
  real, allocatable, dimension(:  ) ::  grid_x_t,grid_y_t  ! 1d global grid for diag_mgr
  real, allocatable, dimension(:,:) ::  sin_rot, cos_rot   ! sin/cos of vector rotation angle
  real, allocatable, dimension(:,:) ::  cell_area          ! grid cell area; sphere frac.
  
  !
  ! timestep parameters
  !
  integer            :: evp_sub_steps                ! evp subcycles / timestep
  real               :: dt_evp = 0.0                 ! evp timestep (sec)
  integer            :: adv_sub_steps                ! advection steps / timestep
  real               :: dt_adv = 0.0                 ! advection timestep (sec)
  integer            :: comm_pe                      ! pe to be communicated with

  logical            :: reproduce_siena_201303 = .TRUE.

contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_avg - take area weighted average over ice partiions                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function ice_avg(x,part)
    real, dimension(:,:,:),    intent(in) :: x
    real, dimension(:,:,:),    intent(in) :: part
    real, dimension(size(x,1),size(x,2) ) :: ice_avg, ice_part

    integer :: i, j, k

    ice_part = 0.0
    ice_avg  = 0.0
    if (size(x,3)==km) then
       do k=2,km
          do j = 1, size(x,2)
             do i = 1, size(x,1)
                ice_avg(i,j) = ice_avg(i,j) + part(i,j,k)*x(i,j,k)
             enddo
          enddo
       enddo
    else if (size(x,3)==km-1) then
       do k=2,km
          do j = 1, size(x,2)
             do i = 1, size(x,1)
                ice_avg(i,j) = ice_avg(i,j) + part(i,j,k)*x(i,j,k-1)
             enddo
          enddo
       enddo
    end if

    do k=2,km
       ice_part = ice_part + part(:,:,k)
    enddo

    do j = 1, size(x,2)
       do i = 1, size(x,1)
          if(ice_part(i,j) > 0 ) then
             ice_avg(i,j) = ice_avg(i,j)/ice_part(i,j)
          else
             ice_avg(i,j) = 0
          endif
       enddo
    enddo
    return
  end function ice_avg

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! all_avg - take area weighted average over all partiions                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function all_avg(x,part)
    real, dimension(:,:, :),             intent(in) :: x
    real, dimension(isc:,jsc:,:), intent(in) :: part
    real, dimension(isc:iec,jsc:jec)         :: all_avg

    integer :: k

    all_avg = 0.0
    if (size(x,3)==km) then
       do k=1,km
          all_avg = all_avg + part(:,:,k)*x(:,:,k)
       end do
    else
       do k=2,km
          all_avg = all_avg + part(:,:,k)*x(:,:,k-1)
       end do
    end if
    return
  end function all_avg

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
  ! uv_to_t - average v-values to t-points and apply mask                        !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine uv_to_t(uv, t)
    real,    intent(in   ), dimension(isd:,jsd:) :: uv
    real,    intent(  out), dimension(isc:,jsc:) :: t

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          if(wett(i,j) > 0.5 ) then
             t(i,j) = 0.25*(uv(i,j) + uv(i,j-1) + uv(i-1,j) + uv(i-1,j-1) )   
          else
             t(i,j) = 0.0
          endif
       enddo
    enddo

  end subroutine uv_to_t

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! dTdx - eastward difference of tracer, result on uv cells                     !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function dTdx(T)
    real, intent(in   ), dimension(isd:,jsd:) :: T 
    real, dimension(isc:iec,jsc:jec)          :: dTdx

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          DTdx(i,j) = 0.5*(T(i+1,j+1) - T(i,j+1) + T(i+1,j) - T(i,j) )
       enddo
    enddo

    return
  end function dTdx

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! dTdy - northward difference of tracer, result on uv cells                    !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function dTdy(T)
    real, intent(in   ), dimension(isd:,jsd:) :: T 
    real, dimension(isc:iec,jsc:jec)          :: dTdy

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          DTdy(i,j) = 0.5*(T(i+1,j+1) - T(i+1,j) + T(i,j+1) - T(i,j) )
       enddo
    enddo

    return
  end function dTdy

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! t_on_uv - average tracer to uv cells                                      !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function t_on_uv(T)
    real, intent(in   ), dimension(isd:,jsd:) :: T 
    real, dimension(isc:iec,jsc:jec)          :: t_on_uv

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          t_on_uv(i,j) = 0.25*(T(i+1,j+1)+T(i+1,j)+T(i,j+1)+T(i,j) )
       enddo
    enddo

    return
  end function t_on_uv

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! t_to_uv - average t-values to v-points and apply mask                        !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine t_to_uv(t, uv)
    real,    intent(in ), dimension(isd:,jsd:) :: t
    real,    intent(out), dimension(isc:,jsc:) :: uv

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          if(wetv(i,j) > 0.5 ) then
             uv(i,j) = 0.25*( t(i+1, j+1)+t(i+1, j) + t(i,j+1)+t(i,j) )
          else
             uv(i,j) = 0.0
          endif
       enddo
    enddo

  end subroutine t_to_uv

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! vel_t_to_uv - average vel component on t-points to v-points and apply mask   !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine vel_t_to_uv(tx, ty, uvx, uvy)
    real,    intent(in   ), dimension(isd:,jsd:) :: tx, ty
    real,    intent(  out), dimension(isc:,jsc:) :: uvx, uvy

    integer :: i, j

    do j = jsc, jec
       do i = isc, iec
          if( wetv(i,j) > 0.5 ) then
             uvx(i,j) = 0.25*( tx(i+1, j+1)+tx(i+1, j) + tx(i,j+1)+tx(i,j) )
             uvy(i,j) = 0.25*( ty(i+1, j+1)+ty(i+1, j) + ty(i,j+1)+ty(i,j) )
          else
             uvx(i,j) = 0.0
             uvy(i,j) = 0.0
          endif
       enddo
    enddo

  end subroutine vel_t_to_uv

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! set_ice_grid - initialize sea ice grid for dynamics and transport            !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine set_ice_grid(grid, ice_domain, dt_slow, dyn_sub_steps_in, &
                          adv_sub_steps_in, km_in, layout, io_layout, maskmap )
    type(sea_ice_grid_type), intent(inout) :: grid
    type(domain2D),        intent(inout) :: ice_domain
    real,                  intent(in)    :: dt_slow
    integer,               intent(in)    :: dyn_sub_steps_in
    integer,               intent(in)    :: adv_sub_steps_in
    integer,               intent(in)    :: km_in
    integer, dimension(2), intent(inout) :: layout
    integer, dimension(2), intent(inout) :: io_layout
    logical, optional,     intent(in)    :: maskmap(:,:)

    real                                :: angle, lon_scale
    integer                             :: i, j, m, ntiles, ncontacts, outunit
    integer                             :: dims(4)
    real, allocatable, dimension(:,:,:) :: x_vert_t, y_vert_t
    real, allocatable,   dimension(:,:) :: geo_latv   
    real, allocatable,   dimension(:,:) :: geo_lonv  
    real, allocatable,   dimension(:,:) :: depth, tmpx, tmpy, tmp_2d
    integer, dimension(2)               :: tile1, tile2
    integer, dimension(2)               :: istart1, iend1, jstart1, jend1
    integer, dimension(2)               :: istart2, iend2, jstart2, jend2
    character(len=80)                   :: domainname
    character(len=128)                  :: grid_file, ocean_topog
    character(len=256)                  :: ocean_hgrid, ocean_mosaic, attvalue
    type(domain2d)                      :: domain2
    integer                             :: isg, ieg, jsg, jeg
    integer                             :: is,  ie,  js,  je
    integer                             :: grid_version
    integer, parameter                  :: VERSION_0 = 0
    integer, parameter                  :: VERSION_1 = 1
    integer, parameter                  :: VERSION_2 = 2
    integer                             :: ni,nj,siz(4)
    integer, allocatable, dimension(:)  :: pelist, islist, ielist, jslist, jelist
    integer                             :: npes, p
    logical                             :: symmetrize, ndivx_is_even, im_is_even

    grid_file = 'INPUT/grid_spec.nc'
    ocean_topog = 'INPUT/topog.nc'
    outunit = stdout()
    !--- first determine the grid version
    if(field_exist(grid_file, 'ocn_mosaic_file') .or. field_exist(grid_file, 'gridfiles') ) then ! read from mosaic file
       write(outunit,*) '==>Note from ice_grid_mod(set_ice_grid): read grid from mosaic version grid'
       grid_version = VERSION_2
       if( field_exist(grid_file, 'ocn_mosaic_file') ) then ! coupler mosaic
          call read_data(grid_file, "ocn_mosaic_file", ocean_mosaic)
          ocean_mosaic = "INPUT/"//trim(ocean_mosaic)
       else
          ocean_mosaic = trim(grid_file)
       end if
       ntiles = get_mosaic_ntiles(ocean_mosaic)
       if(ntiles .NE. 1) call mpp_error(FATAL, '==>Error from ice_grid_mod(set_ice_grid): '//&
            'ntiles should be 1 for ocean mosaic, contact developer')
       call read_data(ocean_mosaic, "gridfiles", ocean_hgrid)
       ocean_hgrid = 'INPUT/'//trim(ocean_hgrid)
       call field_size(ocean_hgrid, 'x', siz)
       ni = siz(1)/2
       nj = siz(2)/2
    else  if(field_exist(grid_file, 'x_T')) then
       ocean_hgrid = grid_file
       write(outunit,*) '==>Note from ice_grid_mod(set_ice_grid): read grid from new version grid'
       grid_version = VERSION_1
       call field_size( ocean_hgrid, 'x_T', siz)
       ni = siz(1)
       nj = siz(2)
    else if(field_exist(grid_file, 'geolon_t')) then
       ocean_hgrid = grid_file
       write(outunit,*) '==>Note from ice_grid_mod(set_ice_grid): read grid from old version grid'
       grid_version = VERSION_0 
       call field_size( ocean_hgrid, 'geolon_t', siz)  
       ni = siz(1)
       nj = siz(2)
    else
       call mpp_error(FATAL, '==>Error from ice_grid_mod(set_ice_grid): '//&
            'x_T, geolon_t, ocn_mosaic_file, gridfiles does not exist in file ' //trim(grid_file))
    end if

    x_cyclic = .false.; tripolar_grid = .false.
    if(grid_version == VERSION_2) then
       if(field_exist(ocean_mosaic, "contacts") ) then
          ncontacts = get_mosaic_ncontacts(ocean_mosaic)
          if(ncontacts < 1) call mpp_error(FATAL,'==>Error from ice_grid_mod(set_ice_grid): '//&
               'number of contacts should be larger than 0 when field contacts exist in file '//trim(ocean_mosaic) )
          if(ncontacts > 2) call mpp_error(FATAL,'==>Error from ice_grid_mod(set_ice_grid): '//&
               'number of contacts should be no larger than 2')
          call get_mosaic_contact( ocean_mosaic, tile1(1:ncontacts), tile2(1:ncontacts),           &
               istart1(1:ncontacts), iend1(1:ncontacts), jstart1(1:ncontacts), jend1(1:ncontacts), &
               istart2(1:ncontacts), iend2(1:ncontacts), jstart2(1:ncontacts), jend2(1:ncontacts)  )
          do m = 1, ncontacts
             if(istart1(m) == iend1(m) ) then  ! x-direction contact, only cyclic condition
                if(istart2(m) .NE. iend2(m) ) call mpp_error(FATAL,  &
                     "==>Error from ice_grid_mod(set_ice_grid): only cyclic condition is allowed for x-boundary")
                x_cyclic = .true.
             else if( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic or folded-north
                if( jstart1(m) == jstart2(m) ) then ! folded north
                   tripolar_grid=.true.
                else 
                   call mpp_error(FATAL, "==>Error from ice_grid_mod(set_ice_grid): "//&
                     "only folded-north condition is allowed for y-boundary")
                end if
             else 
                call mpp_error(FATAL,  &
                     "==>Error from ice_grid_mod(set_ice_grid): invalid boundary contact")
             end if
          end do
       end if
       !--- get grid size
       call field_size(ocean_hgrid, 'x', dims)
       if(mod(dims(1),2) .NE. 1) call mpp_error(FATAL, '==>Error from ice_grid_mod(set_ice_grid): '//&
            'x-size of x in file '//trim(ocean_hgrid)//' should be 2*ni+1')
       if(mod(dims(2),2) .NE. 1) call mpp_error(FATAL, '==>Error from ice_grid_mod(set_ice_grid): '//&
            'y-size of x in file '//trim(ocean_hgrid)//' should be 2*nj+1')
       im = dims(1)/2 
       jm = dims(2)/2 
    else
       if( get_global_att_value(ocean_hgrid, "x_boundary_type", attvalue) ) then
          if(attvalue == 'cyclic') x_cyclic = .true.
       end if
       if( get_global_att_value(ocean_hgrid, "y_boundary_type", attvalue) ) then
          if(attvalue == 'cyclic') then
             call mpp_error(FATAL, "==>Error from ice_grid_mod(set_ice_grid): "//&
                  "y-cyclic boundary condition is not implemented yet.")
          else if(attvalue == 'fold_north_edge') then
             tripolar_grid = .true.
          end if
       end if
       !--- get the grid size
       call field_size(ocean_hgrid, 'wet', dims)
       im = dims(1)
       jm = dims(2)
    end if

    if(x_cyclic) then
       call mpp_error(NOTE, "==>Note from ice_grid_mod: x_boundary_type is cyclic")
    else
       call mpp_error(NOTE, "==>Note from ice_grid_mod: x_boundary_type is solid_walls")
    end if
   if(tripolar_grid) then
       call mpp_error(NOTE, "==>Note from ice_grid_mod: y_boundary_type is fold_north_edge")
    else
       call mpp_error(NOTE, "==>Note from ice_grid_mod: y_boundary_type is solid_walls")
    end if

    ! default is merdional domain decomp. to load balance xgrid
    if( layout(1)==0 .and. layout(2)==0 ) layout=(/ mpp_npes(), 1 /)
    if( layout(1)/=0 .and. layout(2)==0 ) layout(2) = mpp_npes()/layout(1)
    if( layout(1)==0 .and. layout(2)/=0 ) layout(1) = mpp_npes()/layout(2)
    domainname = 'ice model'
    if(tripolar_grid) then    
       !z1l: Tripolar grid requires symmetry in i-direction domain decomposition
       ndivx_is_even = (mod(layout(1),2) == 0)
       im_is_even    = (mod(im,2) == 0)
       symmetrize = ( ndivx_is_even .AND. im_is_even ) .OR. &
               (  (.NOT.ndivx_is_even) .AND.  (.NOT.im_is_even) ) .OR. &
               (  (.NOT.ndivx_is_even) .AND. im_is_even .AND. layout(1) .LT. im/2 )

       if( .not. symmetrize) then
          call mpp_error(FATAL, "ice_model(set_ice_grid): tripolar regrid requires symmetry in i-direction domain decomposition")
       endif
       call mpp_define_domains( (/1,im,1,jm/), layout, Domain, maskmap=maskmap,   &
                                xflags=CYCLIC_GLOBAL_DOMAIN, xhalo=1,             &
                                yflags=FOLD_NORTH_EDGE, yhalo=1, name=domainname )
    else if(x_cyclic) then
       call mpp_define_domains( (/1,im,1,jm/), layout, Domain, maskmap=maskmap,   &
                                xflags=CYCLIC_GLOBAL_DOMAIN, xhalo=1, yhalo=1, name=domainname )
    else
       call mpp_define_domains( (/1,im,1,jm/), layout, Domain, maskmap=maskmap,   &
                                xhalo=1, yhalo=1, name=domainname )
    endif
    call mpp_define_io_domain(Domain, io_layout)

    call mpp_get_compute_domain( Domain, isc, iec, jsc, jec )
    call mpp_get_data_domain( Domain, isd, ied, jsd, jed )
    call mpp_get_global_domain( Domain, isg, ieg, jsg, jeg )

    call mpp_define_domains( (/1,im,1,jm/), layout, ice_domain, maskmap=maskmap, name='ice_nohalo') ! domain without halo
    call mpp_define_io_domain(ice_domain, io_layout)

  ! Allocate and fill in default values for elements of the sea ice grid type.
  
  grid%isc = isc ; grid%iec = iec ; grid%jsc = jsc ; grid%jec = jec
  grid%isd = isd ; grid%ied = ied ; grid%jsd = jsd ; grid%jed = jed
  grid%isg = isg ; grid%ieg = ieg ; grid%jsg = jsg ; grid%jeg = jeg

  grid%symmetric = .false.          ! ###This will be made run-time selectable later.
  grid%nonblocking_updates = .true. ! ###This will be made run-time selectable later.

  grid%IscB = grid%isc ; grid%JscB = grid%jsc
  grid%IsdB = grid%isd ; grid%JsdB = grid%jsd
  grid%IsgB = grid%isg ; grid%JsgB = grid%jsg
  if (grid%symmetric) then
    grid%IscB = grid%isc-1 ; grid%JscB = grid%jsc-1
    grid%IsdB = grid%isd-1 ; grid%JsdB = grid%jsd-1
    grid%IsgB = grid%isg-1 ; grid%JsgB = grid%jsg-1
  endif
  grid%IecB = grid%iec ; grid%JecB = grid%jec
  grid%IedB = grid%ied ; grid%JedB = grid%jed
  grid%IegB = grid%ieg ; grid%JegB = grid%jeg

  call allocate_metrics(grid)

    allocate ( geo_lonv(isc:iec+1,jsc:jec+1), geo_latv(isc:iec+1,jsc:jec+1) )
    
    allocate ( wett   (isd:ied,jsd:jed), wetv   ( isc:iec,jsc:jec),  &
               dxt    (isd:ied,jsd:jed), dxv     (isd:ied,jsd:jed),  &
               dyt    (isd:ied,jsd:jed), dyv     (isd:ied,jsd:jed),  &
               dtw    (isc:iec,jsc:jec), dte     (isd:ied,jsd:jed),  &
               dts    (isc:iec,jsc:jec), dtn     (isd:ied,jsd:jed),  &
               dxdy   (isc:iec,jsc:jec), dydx    (isc:iec,jsc:jec),  &
               geo_lat(isc:iec,jsc:jec), geo_lon (isc:iec,jsc:jec),  &
               cor    (isc:iec,jsc:jec), sin_rot (isd:ied,jsd:jed),  &
               cos_rot(isd:ied,jsd:jed), latitude(isc:iec,jsc:jec) )
    allocate ( geo_latv_ib(isc:iec,jsc:jec), geo_lonv_ib(isc:iec,jsc:jec) )
    !--- read data from grid_spec.nc
    wett = 0.0
    if(grid_version == VERSION_2) then
       allocate(depth(isc:iec,jsc:jec))
       call read_data(ocean_topog, 'depth', depth(isc:iec,jsc:jec), Domain)
       do j = jsc, jec
          do i = isc, iec
             if(depth(i,j) > 0) wett(i,j) = 1.0
          end do
       end do
       deallocate(depth)
    else
       call read_data(grid_file, 'wet', wett(isc:iec,jsc:jec), Domain)
    end if
    call mpp_update_domains(wett, Domain)

    do j = jsc, jec
       do i = isc, iec
          if( wett(i,j)>0.5 .and. wett(i,j+1)>0.5 .and. wett(i+1,j)>0.5 .and. wett(i+1,j+1)>0.5 ) then
             wetv(i,j) = 1.0
          else
             wetv(i,j) = 0.0
          endif
       enddo
    enddo

    if(tripolar_grid) then
       if (jsc==1.and.any(wett(:,jsc)>0.5)) call error_mesg ('ice_model_mod', &
            'ice model requires southernmost row of land', FATAL);
    endif

    evp_sub_steps = dyn_sub_steps_in
    if (evp_sub_steps>0) then
       dt_evp = dt_slow/evp_sub_steps
    else                            ! evp_sub_steps==0 means no dynamics
       dt_evp = dt_slow              ! but set dt>0 to avoid divide by zero
    end if

    adv_sub_steps = adv_sub_steps_in
    if (adv_sub_steps>0) then
       dt_adv = dt_slow/adv_sub_steps
    else                            ! adv_sub_steps==0 means no advection
       dt_adv = dt_slow              ! but set dt>0 to avoid divide by zero
    end if

    km = km_in

    allocate ( cell_area(isc:iec,jsc:jec) )
    
    allocate (grid_x_t(ni))
    allocate (grid_y_t(nj))    
    grid_x_t = 0.0;grid_y_t = 0.0 

    select case(grid_version)
    case(VERSION_0)
       call mpp_copy_domain(domain, domain2)
       call mpp_set_compute_domain(domain2, isc, iec+1, jsc, jec+1, iec-isc+2, jec-jsc+2 )
       call mpp_set_data_domain   (domain2, isd, ied+1, jsd, jed+1, ied-isd+2, jed-jsd+2 )   
       call mpp_set_global_domain (domain2, isg, ieg+1, jsg, jeg+1, ieg-isg+2, jeg-jsg+2 )
       call read_data(grid_file, 'geolon_vert_t', geo_lonv, domain2)
       call read_data(grid_file, 'geolat_vert_t', geo_latv, domain2)  
       call read_data(grid_file, 'AREA_OCN', cell_area, Domain)    
       !Read similar to ocean
       call read_data(grid_file, "gridlon_t", grid_x_t, no_domain = .true.)
       call read_data(grid_file, "gridlat_t", grid_y_t, no_domain = .true.)
    case(VERSION_1)
       allocate ( x_vert_t(isc:iec,jsc:jec,4), y_vert_t(isc:iec,jsc:jec,4) )
       call read_data(grid_file, 'x_vert_T', x_vert_t,  Domain)
       call read_data(grid_file, 'y_vert_T', y_vert_t,  Domain)
       call read_data(grid_file, 'AREA_OCN', cell_area, Domain)
       geo_lonv(isc:iec,jsc:jec) = x_vert_t(isc:iec,jsc:jec,1)
       geo_lonv(iec+1,  jsc:jec) = x_vert_t(iec,jsc:jec,    2)
       geo_lonv(isc:iec,jec+1  ) = x_vert_t(isc:iec,jec,    4)
       geo_lonv(iec+1,  jec+1  ) = x_vert_t(iec,jec,        3)
       geo_latv(isc:iec,jsc:jec) = y_vert_t(isc:iec,jsc:jec,1)
       geo_latv(iec+1,  jsc:jec) = y_vert_t(iec,jsc:jec,    2)
       geo_latv(isc:iec,jec+1  ) = y_vert_t(isc:iec,jec,    4)
       geo_latv(iec+1,  jec+1  ) = y_vert_t(iec,    jec,    3)
       deallocate(x_vert_t, y_vert_t)
    case(VERSION_2)
       call mpp_copy_domain(domain, domain2)
       call mpp_set_compute_domain(domain2, 2*isc-1, 2*iec+1, 2*jsc-1, 2*jec+1, 2*(iec-isc)+3, 2*(jec-jsc)+3 )
       call mpp_set_data_domain   (domain2, 2*isd-1, 2*ied+1, 2*jsd-1, 2*jed+1, 2*(ied-isd)+3, 2*(jed-jsd)+3 )   
       call mpp_set_global_domain (domain2, 2*isg-1, 2*ieg+1, 2*jsg-1, 2*jeg+1, 2*(ieg-isg)+3, 2*(jeg-jsg)+3 )   
       call mpp_get_compute_domain(domain2, is, ie, js, je)
       if(is .NE. 2*isc-1 .OR. ie .NE. 2*iec+1 .OR. js .NE. 2*jsc-1 .OR. je .NE. 2*jec+1) then
          call mpp_error(FATAL, 'ice_grid_mod: supergrid domain is not set properly')
       endif
       allocate(tmpx(is:ie, js:je), tmpy(is:ie, js:je) )
       call read_data(ocean_hgrid, 'x', tmpx, domain2)
       call read_data(ocean_hgrid, 'y', tmpy, domain2)     
       do j = jsc, jec+1
          do i = isc, iec+1
             geo_lonv(i,j) = tmpx(2*i-1,2*j-1)
             geo_latv(i,j) = tmpy(2*i-1,2*j-1)
          end do
       end do
       call calc_mosaic_grid_area(geo_lonv(isc:iec+1, jsc:jec+1)*pi/180, geo_latv(isc:iec+1, jsc:jec+1)*pi/180, cell_area)
       cell_area = cell_area/(4*PI*RADIUS*RADIUS)
       deallocate(tmpx, tmpy)
       do j = jsc, jec
          do i = isc, iec
             if(wett(i,j) ==0) cell_area(i,j) = 0.0
          end do
       end do       
       call mpp_deallocate_domain(domain2)
    end select

    ! Determine the domain-averaged latitudes and longitudes on the grid for
    ! diagnostic purposes. This is a hold-over and will be changed to match
    ! the MOM6 convention of using the "nominal" latitudes and longitudes,
    ! which are the actual lat & lon on spherical portions of the grid.
    allocate ( xb1d (im+1), yb1d (jm+1) )
    if(PRESENT(maskmap)) then    
       allocate(tmpx(im+1,jm+1), tmpy(im+1,jm+1))
       call get_grid_cell_vertices('OCN', 1, tmpx, tmpy)
       xb1d = sum(tmpx,2)/(jm+1)
       yb1d = sum(tmpy,1)/(im+1)
       deallocate(tmpx, tmpy)
    else
       allocate ( tmpx(isc:iec+1, jm+1) )
       call mpp_set_domain_symmetry(Domain, .TRUE.)
       call mpp_global_field(Domain, geo_lonv, tmpx, flags=YUPDATE, position=CORNER)
       allocate ( tmp_2d(isc:iec+1, jsc:jec+1) )
       tmp_2d = 0
       tmp_2d(isc:iec+1,jsc) = sum(tmpx,2)/(jm+1);
       deallocate(tmpx)
       allocate ( tmpx(im+1, jsc:jec+1) )

       call mpp_global_field(Domain, tmp_2d, tmpx, flags=XUPDATE, position=CORNER)
       xb1d = tmpx(:,jsc)
       deallocate(tmpx, tmp_2d)

       allocate ( tmpy(im+1, jsc:jec+1) )
       call mpp_global_field(Domain, geo_latv, tmpy, flags=XUPDATE, position=CORNER)
       allocate ( tmp_2d(isc:iec+1, jsc:jec+1) )
       tmp_2d = 0
       tmp_2d(isc,jsc:jec+1) = sum(tmpy,1)/(im+1);
       deallocate(tmpy)
       allocate ( tmpy(isc:iec+1, jm+1) )
       call mpp_global_field(Domain, tmp_2d, tmpy, flags=YUPDATE, position=CORNER)
       yb1d = tmpy(isc,:)
       deallocate(tmpy, tmp_2d)
       call mpp_set_domain_symmetry(Domain, .FALSE.)
    endif

  ! Note here that geo_lonv & geo_latv are allocated on a SW grid, whereas the
  ! indexing convention everywhere else in SIS2 is to use a NE grid.
  do J=jsc-1,jec ; do I=isc-1,iec
    grid%geoLatBu(I,J) = geo_latv(I+1,J+1)
    grid%geoLonBu(I,J) = geo_lonv(I+1,J+1)
  enddo ; enddo

  dte(:,:) = 0.0 ! ; dtw(:,:) = 0.0
  dtn(:,:) = 0.0 ! ; dts(:,:) = 0.0
  cos_rot(:,:) = 0.0
  sin_rot(:,:) = 0.0

  do j=jsc,jec ; do i=isc,iec
    dts(i,j)     = edge_length(grid%geoLonBu(I-1,J-1), grid%geoLatBu(I-1,J-1), &
                               grid%geoLonBu(I,J-1), grid%geoLatBu(I,J-1))
    dtn(i,j)     = edge_length(grid%geoLonBu(I-1,J), grid%geoLatBu(I-1,J), &
                               grid%geoLonBu(I,J), grid%geoLatBu(I,J))
    dtw(i,j)     = edge_length(grid%geoLonBu(I-1,J-1),  grid%geoLatBu(I-1,J-1), &
                               grid%geoLonBu(I-1,J), grid%geoLatBu(I-1,J))
    dte(i,j)     = edge_length(grid%geoLonBu(I,J-1),grid%geoLatBu(I,J-1), &
                               grid%geoLonBu(I,J),grid%geoLatBu(I,J))
    lon_scale    = cos((grid%geoLatBu(I-1,J-1) + grid%geoLatBu(I,J-1  ) + &
                        grid%geoLatBu(I-1,J) + grid%geoLatBu(I,J)) * atan(1.0)/180)
    angle        = atan2((grid%geoLonBu(I-1,J) + grid%geoLonBu(I,J) - &
                          grid%geoLonBu(I-1,J-1) - grid%geoLonBu(I,J-1))*lon_scale, &
                          grid%geoLatBu(I-1,J) + grid%geoLatBu(I,J) - &
                          grid%geoLatBu(I-1,J-1) - grid%geoLatBu(I,J-1) )
    sin_rot(i,j) = sin(angle) ! angle is the clockwise angle from lat/lon to ocean
    cos_rot(i,j) = cos(angle) ! grid (e.g. angle of ocean "north" from true north)
  enddo ; enddo

    if(reproduce_siena_201303) then
       call mpp_update_domains(dte, Domain)
       call mpp_update_domains(dtn, Domain)
    else
       call mpp_update_domains(dte, dtn, Domain, gridtype=CGRID_NE, flags = SCALAR_PAIR)
    endif
    call mpp_update_domains(cos_rot, Domain)
    call mpp_update_domains(sin_rot, Domain)

    dxt(:,:) = 0.0
    dyt(:,:) = 0.0

    do j = jsc, jec
       do i = isc, iec
          dxt(i,j) = (dts(i,j) + dtn(i,j) )/2
          if(cell_area(i,j) > 0.0) then
             dyt(i,j) = cell_area(i,j)*4*pi*radius*radius/dxt(i,j)
          else
             dyt(i,j) = (dtw(i,j) + dte(i,j) )/2
          endif
       enddo
    enddo

    call mpp_update_domains(dxt, Domain )
    call mpp_update_domains(dyt, Domain )

    dxv(:,:) = 1.0
    dyv(:,:) = 1.0

    dxv(isc:iec,jsc:jec) = t_on_uv(dxt)
    dyv(isc:iec,jsc:jec) = t_on_uv(dyt)

    if(reproduce_siena_201303) then
       call mpp_update_domains(dxv, Domain )
       call mpp_update_domains(dyv, Domain )
    else
       call mpp_update_domains(dxv, dyv, Domain, gridtype=BGRID_NE, flags=SCALAR_PAIR )
    endif

    !--- dxdy and dydx to be used by ice_dyn_mod.
    dydx(:,:) = dTdx(dyt(:,:))
    dxdy(:,:) = dTdy(dxt(:,:))

  do j=jsc,jec ; do i=isc,iec
    !### REGROUP FOR ROTATIONAL REPRODUCIBILITY
    geo_lon(i,j) = lon_avg( (/ grid%geoLonBu(I-1,J-1), grid%geoLonBu(I,J-1), &
                               grid%geoLonBu(I-1,J), grid%geoLonBu(I,J) /) )
    geo_lat(i,j) = (grid%geoLatBu(I-1,J-1)   + grid%geoLatBu(I,J-1) + &
                    grid%geoLatBu(I-1,J) + grid%geoLatBu(I,J)) / 4
  enddo ; enddo

  do J=jsc,jec ; do I=isc,iec
    cor(I,J) = 2*omega*sin(grid%geoLatBu(I,J)*pi/180)
  enddo ; enddo
  latitude(:,:) = geo_lat(isc:iec,jsc:jec)*pi/180

    deallocate (geo_lonv, geo_latv)

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
                call mpp_error(FATAL, "ice_model: jelist(p) .NE. jec but jslist(p) == jsc")
             endif
             if( ielist(p) + isc .NE. im+1) then
                call mpp_error(FATAL, "ice_model: ielist(p) + isc .NE. im+1 but islist(p) + iec == im+1")
             endif
             comm_pe = pelist(p)
             exit 
          endif
       enddo
       deallocate(pelist, islist, ielist, jslist, jelist)
    endif
!    comm_pe = mpp_pe() + layout(1) - 2*mod(mpp_pe()-mpp_root_pe(),layout(1)) - 1

    return
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

  deallocate(wett, wetv, dtw, dte, dts, dtn, dxt, dxv, dyt, dyv, latitude )
  deallocate(cor, geo_lat, geo_lon, xb1d, yb1d, sin_rot, cos_rot, cell_area )


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
  subroutine ice_line(year, day, second, cn, sst)
    integer,                               intent(in) :: year, day, second
    real, dimension(isc:,jsc:,1:),   intent(in) :: cn
    real, dimension(isc:,jsc:),      intent(in) :: sst

    real, dimension(isc:iec,jsc:jec) :: x
    real                             :: gx(3)
    integer                          :: i

    do i=-1,1,2
       x = 0.0
       where (cn(:,:,1)<0.85 .and. i*geo_lat>0.0) x=cell_area
       gx((i+3)/2) = g_sum(x)*4*pi*radius*radius/1e12
    end do
    gx(3) = g_sum(sst*cell_area)/(g_sum(cell_area)+1e-10)
    !
    ! print info every 5 days
    !
    if ( mpp_pe()==0 .and. second==0 .and. mod(day,5)==0 ) &
       print '(a,2I4,3F10.5)','ICE y/d (SH_ext NH_ext SST):', year, day, gx
  end subroutine ice_line

  !#####################################################################
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! ice_advect - take adv_sub_steps upstream advection timesteps                 !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine ice_advect(uc, vc, trc, uf, vf)
    real, intent(in   ),           dimension(isd:,jsd:) :: uc, vc  ! advecting velocity on C-grid
    real, intent(inout),           dimension(isd:,jsd:) :: trc     ! tracer to advect
    real, optional, intent(inout), dimension(isc:,jsc:) :: uf, vf

    integer                          :: l, i, j
    real, dimension(isd:ied,jsd:jed) ::  uflx, vflx

    if (adv_sub_steps==0) return;

    if (present(uf)) uf = 0.0
    if (present(vf)) vf = 0.0

    uflx = 0.0
    vflx = 0.0

    do l=1,adv_sub_steps
       do j = jsd, jec
          do i = isd, iec
             if( uc(i,j) > 0.0 ) then
                uflx(i,j) = uc(i,j) * trc(i,j) * dte(i,j)
             else
                uflx(i,j) = uc(i,j) * trc(i+1,j) * dte(i,j)
             endif

             if( vc(i,j) > 0.0 ) then
                vflx(i,j) = vc(i,j) * trc(i,j) * dtn(i,j)
             else
                vflx(i,j) = vc(i,j) * trc(i,j+1) * dtn(i,j)
             endif
          enddo
       enddo

       do j = jsc, jec
          do i = isc, iec
             trc(i,j) = trc(i,j) + dt_adv * ( uflx(i-1,j) - uflx(i,j) + &
                        vflx(i,j-1) - vflx(i,j) )/ ( dxt(i,j) * dyt(i,j) ) 
          enddo
       enddo

       call mpp_update_domains(trc, Domain)

       if (present(uf)) then
          do j = jsc, jec
             do i = isc, iec       
                uf(i,j) = uf(i,j) + uflx(i,j)
             enddo
          enddo
       endif

       if (present(vf)) then
          do j = jsc, jec
             do i = isc, iec       
                vf(i,j) = vf(i,j) + vflx(i,j)
             enddo
          enddo
       endif

    end do

    if (present(uf)) uf = uf/adv_sub_steps;
    if (present(vf)) vf = vf/adv_sub_steps;

  end subroutine ice_advect

  !#####################################################################
  subroutine slab_ice_advect(ui, vi, trc, stop_lim)
    real, intent(in   ), dimension(isd:,jsd:) :: ui, vi       ! advecting velocity
    real, intent(inout), dimension(isd:,jsd:) :: trc          ! tracer to advect
    real, intent(in   )                             :: stop_lim

    integer                          :: l, i, j
    real, dimension(isd:ied,jsd:jed) :: ue, vn, uflx, vflx
    real                             :: avg, dif

    if (adv_sub_steps==0) return;

    ue(:,jsd) = 0.0; ue(:,jed) = 0.0
    vn(:,jsd) = 0.0; vn(:,jed) = 0.0

    do j = jsc, jec
       do i = isc, iec
          ue(i,j) = 0.5 * ( ui(i,j-1) + ui(i,j) )
          vn(i,j) = 0.5 * ( vi(i-1,j) + vi(i,j) )
       enddo
    enddo

    if(reproduce_siena_201303) then
       call mpp_update_domains(ue, Domain)
       call mpp_update_domains(vn, Domain)
    else
       call mpp_update_domains(ue, vn, Domain, gridtype=CGRID_NE)
    endif

    do l=1,adv_sub_steps
       do j = jsd, jec
          do i = isd, iec
             avg = ( trc(i,j) + trc(i+1,j) )/2
             dif = trc(i+1,j) - trc(i,j)
             if( avg > stop_lim .and. ue(i,j) * dif > 0.0) then
                uflx(i,j) = 0.0
             else if( ue(i,j) > 0.0 ) then
                uflx(i,j) = ue(i,j) * trc(i,j) * dte(i,j)
             else
                uflx(i,j) = ue(i,j) * trc(i+1,j) * dte(i,j)
             endif

             avg = ( trc(i,j) + trc(i,j+1) )/2
             dif = trc(i,j+1) - trc(i,j)
             if( avg > stop_lim .and. vn(i,j) * dif > 0.0) then
                vflx(i,j) = 0.0
             else if( vn(i,j) > 0.0 ) then
                vflx(i,j) = vn(i,j) * trc(i,j) * dtn(i,j)
             else
                vflx(i,j) = vn(i,j) * trc(i,j+1) * dtn(i,j)
             endif
          enddo
       enddo

       do j = jsc, jec
          do i = isc, iec
             trc(i,j) = trc(i,j) + dt_adv * ( uflx(i-1,j)-uflx(i,j) + vflx(i,j-1)-vflx(i,j) ) / (dxt(i,j)*dyt(i,j))
          enddo
       enddo

       call mpp_update_domains(trc, Domain)

    end do
  end subroutine slab_ice_advect

  !#####################################################################
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! cut_check - check for mirror symmetry of uv field along bipolar grid cut     !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine cut_check(mesg, uv)
    character(len=*), intent(in) :: mesg
    real, intent(inout), dimension(isd:,jsd:) :: uv

!    real, dimension(1:im, 1:jm)  :: global_uv

    integer :: i, l
    real, dimension(iec-isc+1) :: send_buffer, recv_buffer

    if(comm_pe == NULL_PE) return

    call mpp_recv(recv_buffer(1), glen=(iec-isc+1), from_pe=comm_pe, block=.false., tag=COMM_TAG_1)

    l = 0
    do i = iec-1, isd, -1
       l = l+1
       send_buffer(l) = uv(i,jec)
    enddo

    call mpp_send(send_buffer(1), plen=(iec-isc+1), to_pe=comm_pe, tag=COMM_TAG_1)
    call mpp_sync_self(check=EVENT_RECV)

    ! cut_check only at the north pole.
    if( jec == jm) then
       l = 0
       do i = isc, iec
          l=l+1
          ! excludes the point at im/2, and im.
          if(i == im/2 .or. i == im) cycle
          uv(i,jec) = (uv(i,jec) - recv_buffer(l))/2
       enddo
    endif

    call mpp_sync_self()


    !    call mpp_global_field ( Domain, uv, global_uv )

    !    do i=1,im/2-1
    !       if (global_uv(i,jm)/=-global_uv(im-i,jm)) then
    !      cut_error =.true.
    !     if (mpp_pe()==0) &
         !       print *, mesg, i, im-i, global_uv(i,jm), global_uv(im-i,jm)
    !          fix = (global_uv(i,jm)-global_uv(im-i,jm))/2
    !          global_uv(i   ,jm) =  fix
    !          global_uv(im-i,jm) = -fix
    !       end if
    !    end do

    !    uv(isc:iec,jsc:jec) = global_uv(isc:iec,jsc:jec)

    ! if (cut_error) call error_mesg ('ice_model_mod', &
         !                     'cut check of mirror anti-symmetry failed', FATAL)
  end subroutine cut_check
  !#####################################################################

subroutine allocate_metrics(G)
  type(sea_ice_grid_type), intent(inout) :: G
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  ! This subroutine allocates the lateral elements of the sea_ice_grid_type that
  ! are always used and zeros them out.

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

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

  ALLOC_(G%CoriolisBu(IsdB:IedB, JsdB:JedB)) ; G%CoriolisBu(:,:) = 0.0

end subroutine allocate_metrics

end module ice_grid_mod
