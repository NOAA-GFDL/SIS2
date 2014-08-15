!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_grid_mod - sets up grid and processor domain - Michael.Winton@noaa.gov   !
!   This module is in the process of extensive revision to harmonize it with   !
! MOM6.  -Robert Hallberg                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_grid_mod

  use constants_mod,   only: radius, omega, pi, grav

use mpp_domains_mod, only: mpp_define_domains, FOLD_NORTH_EDGE
use mpp_domains_mod, only: domain2D, mpp_global_field, YUPDATE, XUPDATE, CORNER
use mpp_domains_mod, only : CENTER, NORTH_FACE=>NORTH, EAST_FACE=>EAST
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain, mpp_set_domain_symmetry
use mpp_domains_mod, only: mpp_define_io_domain, mpp_copy_domain, mpp_get_global_domain
use mpp_domains_mod, only: mpp_set_global_domain, mpp_set_data_domain, mpp_set_compute_domain
use mpp_domains_mod, only: mpp_deallocate_domain, mpp_get_pelist, mpp_get_compute_domains
use mpp_domains_mod, only : domain1D, mpp_get_domain_components

use MOM_domains, only : SIS_domain_type=>MOM_domain_type, pass_var, pass_vector
use MOM_domains, only : PE_here, root_PE, broadcast, MOM_domains_init, clone_MOM_domain
use MOM_domains, only : num_PEs, SCALAR_PAIR, CGRID_NE, BGRID_NE, To_All, AGRID
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_string_functions, only : slasher

use fms_io_mod,      only : file_exist, parse_mask_table
use fms_mod,         only: field_exist, field_size, read_data
use fms_mod,         only: get_global_att_value, stderr
use mosaic_mod,      only: get_mosaic_ntiles, get_mosaic_ncontacts
use mosaic_mod,      only: calc_mosaic_grid_area, get_mosaic_contact
use grid_mod,        only: get_grid_cell_vertices

implicit none ; private

include 'netcdf.inc'
#include <SIS2_memory.h>

public :: set_ice_grid, ice_grid_end
public :: isPointInCell
public :: cell_area

type, public :: sea_ice_grid_type
  type(SIS_domain_type), pointer :: Domain => NULL()
  type(SIS_domain_type), pointer :: Domain_aux => NULL() ! A non-symmetric auxiliary domain type.
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

  real, allocatable, dimension(:) :: &
    H_cat_lim, &  ! The lower thickness limits for each ice category, in m.
    M_cat_lim     ! The lower mass-per-unit area limits for each ice category,
                  ! in kg m-2.
    
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
  logical, pointer :: maskmap(:,:)=>NULL() ! A pointer to an array indicating
                                ! which logical processors are actually used for
                                ! the ocean code. The other logical processors
                                ! would be all land points and are not assigned
                                ! to actual processors. This need not be
                                ! assigned if all logical processors are used.
end type SIS2_domain_type

! This is still here as an artefact of an older public interface and should go.
real, allocatable, dimension(:,:) ::  cell_area  ! grid cell area; sphere frac.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! set_ice_grid - initialize sea ice grid for dynamics and transport            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine set_ice_grid(G, param_file, ice_domain, NCat_dflt)
  type(sea_ice_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
  type(domain2D),        intent(inout) :: ice_domain
  integer,               intent(in)    :: NCat_dflt
!   This subroutine sets up the necessary domain types and the sea-ice grid.

! Arguments: G - The sea-ice's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (inout)   ice_domain - A domain with no halos that can be shared publicly.
!  (in)      NCat_dflt - The default number of ice categories.

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
  logical :: set_grid_like_SIS1
  character(len=256) :: grid_file, ocean_topog
  character(len=256) :: ocean_hgrid, ocean_mosaic
  character(len=200) :: mesg
  character(len=40)  :: mod_nm  = "ice_grid" ! This module's name.
  type(domain2d)     :: domain2
  type(domain2d), pointer :: Domain => NULL()

  grid_file = 'INPUT/grid_spec.nc'
  ocean_topog = 'INPUT/topog.nc'

  ! Set up the SIS_domain_type.  This will later occur via a call to MOM_domains_init.
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
  call clone_MOM_domain(G%domain, ice_domain, halo_size=0, symmetric=.false., &
                        domain_name="ice_nohalo")

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod_nm, version)
  call get_param(param_file, mod_nm, "GLOBAL_INDEXING", global_indexing, &
                 "If true, use a global lateral indexing convention, so \n"//&
                 "that corresponding points on different processors have \n"//&
                 "the same index. This does not work with static memory.", &
                 default=.false.)
#ifdef STATIC_MEMORY_
  call get_param(param_file, mod_nm, "NCAT_ICE", G%CatIce, &
                 "The number of sea ice thickness categories.", units="nondim", &
                 default=NCat_dflt)
  if (G%CatIce /= NCAT_ICE_) call SIS_error(FATAL, "set_ice_grid: " // &
       "Mismatched number of categories NCAT_ICE between SIS_memory.h and "//&
       "param_file or the input namelist file.")
  call get_param(param_file, mod_nm, "NK_ICE", G%NkIce, &
                 "The number of layers within the sea ice.", units="nondim", &
                 default=NK_ICE_)
  if (G%NkIce /= NK_ICE_) call SIS_error(FATAL, "set_ice_grid: " // &
       "Mismatched number of layers NK_ICE between SIS_memory.h and param_file")

  call get_param(param_file, mod_nm, "NK_SNOW", G%NkSnow, &
                 "The number of layers within the snow atop the sea ice.", &
                 units="nondim", default=NK_SNOW_)
  if (G%NkSnow /= NK_SNOW_) call SIS_error(FATAL, "set_ice_grid: " // &
       "Mismatched number of layers NK_SNOW between SIS_memory.h and param_file")
  if (global_indexing) cal SIS_error(FATAL, "set_ice_grid : "//&
       "GLOBAL_INDEXING can not be true with STATIC_MEMORY.")
#else
  call get_param(param_file, mod_nm, "NCAT_ICE", G%CatIce, &
                 "The number of sea ice thickness categories.", units="nondim", &
                 default=NCat_dflt)
  call get_param(param_file, mod_nm, "NK_ICE", G%NkIce, &
                 "The number of layers within the sea ice.", units="nondim", &
                 default=4) ! Valid for SIS5L; Perhaps this should be ..., fail_if_missing=.true.
  call get_param(param_file, mod_nm, "NK_SNOW", G%NkSnow, &
                 "The number of layers within the snow atop the sea ice.", &
                 units="nondim", default=1) ! Perhaps this should be ..., fail_if_missing=.true.
#endif

  call get_param(param_file, mod_nm, "SET_GRID_LIKE_SIS1", set_grid_like_SIS1, &
                 "If true, use SIS1 code to set the grid values.  Otherwise \n"//&
                 "use code derived from MOM6.", default=.false.)

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

  ! This code should be moved to MOM_domains_init once we start using a cubed-sphere grid.
  ! if (field_exist(ocean_mosaic, "contacts") ) then
  !   ncontacts = get_mosaic_ncontacts(ocean_mosaic)
  !   if (ncontacts < 1) call SIS_error(FATAL,'==>Error from ice_grid_mod(set_ice_grid): '//&
  !        'number of contacts should be larger than 0 when field contacts exist in file '//&
  !        trim(ocean_mosaic) )
  !   if (ncontacts > 2) call SIS_error(FATAL,'==>Error from ice_grid_mod(set_ice_grid): '//&
  !        'number of contacts should be no larger than 2')
  !   call get_mosaic_contact( ocean_mosaic, tile1(1:ncontacts), tile2(1:ncontacts),           &
  !        istart1(1:ncontacts), iend1(1:ncontacts), jstart1(1:ncontacts), jend1(1:ncontacts), &
  !        istart2(1:ncontacts), iend2(1:ncontacts), jstart2(1:ncontacts), jend2(1:ncontacts)  )
  !   do m = 1, ncontacts
  !     if (istart1(m) == iend1(m) ) then  ! x-direction contact, only cyclic condition
  !       if (istart2(m) /= iend2(m) ) call SIS_error(FATAL,  &
  !            "==>Error from ice_grid_mod(set_ice_grid): only cyclic condition is allowed for x-boundary")
  !       x_cyclic = .true.
  !     elseif ( jstart1(m) == jend1(m) ) then  ! y-direction contact, cyclic or folded-north
  !       if ( jstart1(m) == jstart2(m) ) then ! folded north
  !          tripolar_grid=.true.
  !       else 
  !          call SIS_error(FATAL, "==>Error from ice_grid_mod(set_ice_grid): "//&
  !            "only folded-north condition is allowed for y-boundary")
  !       endif
  !     else 
  !       call SIS_error(FATAL,  &
  !            "==>Error from ice_grid_mod(set_ice_grid): invalid boundary contact")
  !     endif
  !   enddo
  ! endif

  !--- get grid size from the input file hgrid file.
  call field_size(ocean_hgrid, 'x', dims)
  if(mod(dims(1),2) /= 1) call SIS_error(FATAL, '==>Error from ice_grid_mod(set_ice_grid): '//&
      'x-size of x in file '//trim(ocean_hgrid)//' should be 2*niglobal+1')
  if(mod(dims(2),2) /= 1) call SIS_error(FATAL, '==>Error from ice_grid_mod(set_ice_grid): '//&
      'y-size of x in file '//trim(ocean_hgrid)//' should be 2*njglobal+1')
  ni = dims(1)/2 
  nj = dims(2)/2 
  
  if (ni /= G%Domain%niglobal) call SIS_error(FATAL, "set_ice_grid: "//&
    "The total i-grid size from file "//trim(ocean_hgrid)//" is inconsistent with SIS_input.")
  if (nj /= G%Domain%njglobal) call SIS_error(FATAL, "set_ice_grid: "//&
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
  G%ks = 0 ; G%ke = 0  ! Change this for shared ocean / ice grids.

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

  !--- read data from grid_spec.nc
  allocate(depth(G%isc:G%iec,G%jsc:G%jec))
  call read_data(ocean_topog, 'depth', depth(G%isc:G%iec,G%jsc:G%jec), G%Domain%mpp_domain)
  do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (depth(i,j) > 0) then
    G%mask2dT(i,j) = 1.0
  endif ; enddo ; enddo
  deallocate(depth)
  call pass_var(G%mask2dT, G%Domain)

  do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
    if( G%mask2dT(i,j)>0.5 .and. G%mask2dT(i,j+1)>0.5 .and. &
        G%mask2dT(i+1,j)>0.5 .and. G%mask2dT(i+1,j+1)>0.5 ) then
       G%mask2dBu(I,J) = 1.0 ; G%Lmask2dBu(I,J) = .true.
    else
       G%mask2dBu(I,J) = 0.0 ; G%Lmask2dBu(I,J) = .false.
    endif
  enddo ; enddo

  do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
    if( G%mask2dT(i,j)>0.5 .and. G%mask2dT(i+1,j)>0.5 ) then
       G%mask2dCu(I,j) = 1.0 ; G%Lmask2dCu(I,j) = .true.
    else
       G%mask2dCu(I,j) = 0.0 ; G%Lmask2dCu(I,j) = .false.
    endif
  enddo ; enddo

  do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
    if( G%mask2dT(i,j)>0.5 .and. G%mask2dT(i,j+1)>0.5 ) then
       G%mask2dCv(i,J) = 1.0 ; G%Lmask2dCv(i,J) = .true.
    else
       G%mask2dCv(i,J) = 0.0 ; G%Lmask2dCv(i,J) = .false.
    endif
  enddo ; enddo
  call pass_var(G%mask2dBu, G%Domain, position=CORNER)
  call pass_vector(G%mask2dCu, G%mask2dCv, G%Domain, To_All+Scalar_pair)

  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    G%Lmask2dT(i,j) = (G%mask2dT(i,j) > 0.5) 
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied
    G%Lmask2dCv(i,J) = (G%mask2dCv(i,J) > 0.5) 
  enddo ; enddo
  do j=G%jsd,G%jed ; do I=G%IsdB,G%IedB
    G%Lmask2dCu(I,j) = (G%mask2dCu(I,j) > 0.5) 
  enddo ; enddo
  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB
    G%Lmask2dBu(I,J) = (G%mask2dBu(I,J) > 0.5) 
  enddo ; enddo


  if (set_grid_like_SIS1) then
    call mpp_copy_domain(G%Domain%mpp_domain, domain2)
    call mpp_set_compute_domain(domain2, 2*isca-1, 2*ieca+1, 2*jsca-1, 2*jeca+1, 2*(ieca-isca)+3, 2*(jeca-jsca)+3 )
    call mpp_set_data_domain   (domain2, 2*isda-1, 2*ieda+1, 2*jsda-1, 2*jeda+1, 2*(ieda-isda)+3, 2*(jeda-jsda)+3 )   
    call mpp_set_global_domain (domain2, 2*isg-1, 2*ieg+1, 2*jsg-1, 2*jeg+1, 2*(ieg-isg)+3, 2*(jeg-jsg)+3 )   
    call mpp_get_compute_domain(domain2, is, ie, js, je)
    if(is /= 2*isca-1 .or. ie /= 2*ieca+1 .or. js /= 2*jsca-1 .or. je /= 2*jeca+1) then
      call SIS_error(FATAL, 'ice_grid_mod: supergrid domain is not set properly')
    endif
    allocate(tmpx(is:ie, js:je), tmpy(is:ie, js:je) )
    call read_data(ocean_hgrid, 'x', tmpx, domain2)
    call read_data(ocean_hgrid, 'y', tmpy, domain2)     
    do J=jsca-1,jeca ; do I=isca-1,ieca
      G%geoLonBu(I-i_off,J-j_off) = tmpx(2*i+1,2*j+1)
      G%geoLatBu(I-i_off,J-j_off) = tmpy(2*i+1,2*j+1)
    enddo ; enddo
    deallocate(tmpx, tmpy)
    call calc_mosaic_grid_area(G%geoLonBu(G%isc-1:G%iec,G%jsc-1:G%jec)*pi/180, &
                               G%geoLatBu(G%isc-1:G%iec,G%jsc-1:G%jec)*pi/180, &
                               G%areaT(G%isc:G%iec,G%jsc:G%jec))
    call mpp_deallocate_domain(domain2)

!    call pass_var(G%geoLonBu, G%Domain, position=CORNER)
!    call pass_var(G%geoLatBu, G%Domain, position=CORNER)

    ! Determine the domain-averaged latitudes and longitudes on the grid for
    ! diagnostic purposes. This is a hold-over and will be changed to match
    ! the MOM6 convention of using the "nominal" latitudes and longitudes,
    ! which are the actual lat & lon on spherical portions of the grid.
    allocate ( xb1d (ni+1), yb1d (nj+1) )
    if (associated(G%Domain%maskmap)) then    
      allocate(tmpx(ni+1,nj+1), tmpy(ni+1,nj+1))
      call get_grid_cell_vertices('OCN', 1, tmpx, tmpy)
      xb1d(:) = sum(tmpx,2)/(nj+1)
      yb1d(:) = sum(tmpy,1)/(ni+1)
      deallocate(tmpx, tmpy)
    else
      Domain => G%Domain%mpp_domain
      allocate ( tmpx(isca:ieca+1, nj+1) )
      call mpp_set_domain_symmetry(Domain, .TRUE.)
      call mpp_global_field(Domain, G%geoLonBu(G%isc-1:G%iec,G%jsc-1:G%jec), &
                           tmpx, flags=YUPDATE, position=CORNER)
      allocate ( tmp_2d(isca:ieca+1, jsca:jeca+1) )
      tmp_2d = 0
      tmp_2d(isca:ieca+1,jsca) = sum(tmpx,2)/(nj+1);
      deallocate(tmpx)
      allocate ( tmpx(ni+1, jsca:jeca+1) )

      call mpp_global_field(Domain, G%geoLatBu(G%isc-1:G%iec,G%jsc-1:G%jec), &
                           tmpx, flags=XUPDATE, position=CORNER)
      xb1d(:) = tmpx(:,jsca)
      deallocate(tmpx, tmp_2d)

      allocate ( tmpy(ni+1, jsca:jeca+1) )
      call mpp_global_field(Domain, G%geoLatBu(G%isc-1:G%iec,G%jsc-1:G%jec), tmpy, &
                           flags=XUPDATE, position=CORNER)
      allocate ( tmp_2d(isca:ieca+1, jsca:jeca+1) )
      tmp_2d = 0
      tmp_2d(isca,jsca:jeca+1) = sum(tmpy,1)/(ni+1);
      deallocate(tmpy)
      allocate ( tmpy(isca:ieca+1, nj+1) )
      call mpp_global_field(Domain, tmp_2d, tmpy, flags=YUPDATE, position=CORNER)
      yb1d(:) = tmpy(isca,:)
      deallocate(tmpy, tmp_2d)
      call mpp_set_domain_symmetry(Domain, G%domain%symmetric)
    endif

    G%gridLatB(jsg-1:jeg) = yb1d(1:nj+1)
    G%gridLonB(isg-1:ieg) = xb1d(1:ni+1)
    do i=isg,ieg ; G%gridLonT(i) = 0.5*(G%gridLonB(i-1)+G%gridLonB(i)) ; enddo
    do j=jsg,jeg ; G%gridLatT(j) = 0.5*(G%gridLatB(j-1)+G%gridLatB(j)) ; enddo

    deallocate( xb1d, yb1d )

    do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
      G%dyCu(I,j) = edge_length(G%geoLonBu(I,J-1), G%geoLatBu(I,J-1), &
                                G%geoLonBu(I,J), G%geoLatBu(I,J))
    enddo ; enddo
    do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
      G%dxCv(i,J) = edge_length(G%geoLonBu(I-1,J), G%geoLatBu(I-1,J), &
                                G%geoLonBu(I,J), G%geoLatBu(I,J))
    enddo ; enddo
    call pass_vector(G%dyCu, G%dxCv, G%Domain, To_All+Scalar_Pair, CGRID_NE)

    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      G%dxT(i,j) = (G%dxCv(i,J-1) + G%dxCv(I,j) )/2
      if (G%mask2dT(i,j) > 0.0) then
        G%dyT(i,j) = G%areaT(i,j)/G%dxT(i,j)
      else
        G%dyT(i,j) = (G%dyCu(I-1,j) + G%dyCu(I,j) )/2
      endif
    enddo ; enddo
    call pass_vector(G%dxT, G%dyT, G%Domain, To_All+Scalar_Pair, AGRID)

    G%dxBu(:,:) = 1.0 ; G%IdxBu(:,:) = 1.0
    G%dyBu(:,:) = 1.0 ; G%IdyBu(:,:) = 1.0

    do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
      G%dxBu(I,J) = 0.25*((G%dxT(i+1,j+1) + G%dxT(i,j)) + &
                          (G%dxT(i+1,j) + G%dxT(i,j+1)) )
      G%IdxBu(I,J) = 1.0 / G%dxBu(I,j)
      G%dyBu(I,J) = 0.25*((G%dyT(i+1,j+1) + G%dyT(i,j)) + &
                          (G%dyT(i+1,j) + G%dyT(i,j+1)) )
      G%IdyBu(I,J) = 1.0 / G%dyBu(I,j)
      G%areaBu(I,J) = G%dxBu(I,J) * G%dyBu(I,J)
    enddo ; enddo
    call pass_vector(G%dxBu, G%dyBu, G%Domain, To_All+Scalar_Pair, BGRID_NE)
    call pass_vector(G%IdxBu, G%IdyBu, G%Domain, To_All+Scalar_Pair, BGRID_NE)

    do j=G%jsc,G%jec ; do I=G%isc-1,G%iec
      G%dxCu(I,j) = 0.5*(G%dxBu(I,J) + G%dxBu(I,J-1))
    enddo ; enddo
    do J=G%jsc-1,G%jec ; do i=G%isc,G%iec
      G%dyCv(i,J) = 0.5*(G%dyBu(I,J) + G%dyBu(I-1,J))
    enddo ; enddo
    call pass_vector(G%dxCu, G%dyCv, G%Domain, To_All+Scalar_Pair, CGRID_NE)

    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      G%geoLonT(i,j) = 0.5 * (lon_avg( (/ G%geoLonBu(I-1,J-1), G%geoLonBu(I,J) /) ) + &
                              lon_avg( (/ G%geoLonBu(I,J-1), G%geoLonBu(I-1,J) /) ) )
      G%geoLatT(i,j) = ((G%geoLatBu(I-1,J-1) + G%geoLatBu(I,J)) + &
                        (G%geoLatBu(I,J-1) + G%geoLatBu(I-1,J)) ) * 0.25
    enddo ; enddo
  else
    call set_grid_metrics_from_mosaic(G, param_file)
  endif
  call set_grid_derived_metrics(G, param_file)

  ! cell_area is unfortunately used outside of the ice model for various
  ! things, so it has to be set, but it should be eliminated. -RWH
  allocate ( cell_area(isca:ieca,jsca:jeca) )
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    cell_area(i+i_off,j+j_off) = G%mask2dT(i,j) * G%areaT(i,j)/(4*PI*RADIUS**2)
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

  do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
    G%CoriolisBu(I,J) = 2*omega*sin(G%geoLatBu(I,J)*pi/180)
  enddo ; enddo

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

end subroutine set_ice_grid

!---------------------------------------------------------------------

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

!---------------------------------------------------------------------
function edge_length(x1, y1, x2, y2)
  real, intent(in) :: x1, x2, y1, y2 ! end-point coordinates in degrees
  real             :: edge_length
  real             :: dx, dy

  dx = (x2-x1)*cos((atan(1.0)/45)*(y2+y1)/2)
  dy = y2-y1
  edge_length = radius*(atan(1.0)/45)*(dx*dx+dy*dy)**0.5
end function edge_length

! ------------------------------------------------------------------------------

subroutine set_grid_derived_metrics(G, param_file)
  type(sea_ice_grid_type), intent(inout) :: G
  type(param_file_type), intent(in)    :: param_file
! Arguments:
!  (inout)   G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.

!    Calculate the values of the metric terms that might be used
!  and save them in arrays.  This should be identical to the corresponding
!  subroutine in MOM_grid_initialize.F90.
!    Within this subroutine, the x- and y- grid spacings and their
!  inverses and the cell areas centered on T, Bu, Cu, and Cv points are
!  calculated, as are the geographic locations of each of these 4
!  sets of points.
  character( len = 128) :: warnmesg
  integer :: i,j, isd, ied, jsd, jed
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, IsdB, IedB, JsdB, JedB

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  call SIS_mesg("  MOM_grid_init.F90, set_grid_derived_metrics: deriving metrics", 5)
 
  do j=jsd,jed ; do i=isd,ied
    if (G%dxT(i,j) < 0.0) then
      write(warnmesg,68)  pe_here(),"dxT",i,j,G%dxT(i,j),0.0
      call SIS_mesg(warnmesg, all_print=.true.)
      G%dxT(i,j) = 0.0
    endif
    if (G%dyT(i,j) < 0.0) then
      write(warnmesg,68)  pe_here(),"dyT",i,j,G%dyT(i,j),0.0
      call SIS_mesg(warnmesg, all_print=.true.)
      G%dyT(i,j) = 0.0
    endif
    G%IdxT(i,j) = Adcroft_reciprocal(G%dxT(i,j))
    G%IdyT(i,j) = Adcroft_reciprocal(G%dyT(i,j))
    G%IareaT(i,j) = Adcroft_reciprocal(G%areaT(i,j))
  enddo ; enddo

  do j=jsd,jed ; do I=IsdB,IedB
    if (G%dxCu(I,j) < 0.0) then
      write(warnmesg,68)  pe_here(),"dxCu",I,j,G%dxCu(I,j),0.0
      call SIS_mesg(warnmesg, all_print=.true.)
      G%dxCu(I,j) = 0.0
    endif
    if (G%dyCu(I,j) < 0.0) then
      write(warnmesg,68)  pe_here(),"dyCu",I,j,G%dyCu(I,j),0.0
      call SIS_mesg(warnmesg, all_print=.true.)
      G%dyCu(I,j) = 0.0
    endif
    G%IdxCu(I,j) = Adcroft_reciprocal(G%dxCu(I,j))
    G%IdyCu(I,j) = Adcroft_reciprocal(G%dyCu(I,j))
    G%areaCu(I,j) = G%dxCu(I,j)*G%dyCu(I,j)  !### Replace with * G%dy_Cu(I,j)?
    G%IareaCu(I,j) = G%mask2dCu(I,j) * Adcroft_reciprocal(G%areaCu(I,j))
  enddo ; enddo

  do J=JsdB,JedB ; do i=isd,ied
    if (G%dxCv(i,J) < 0.0) then
      write(warnmesg,68)  pe_here(),"dxCv",i,j,G%dxCv(i,j),0.0
      call SIS_mesg(warnmesg, all_print=.true.)
      G%dxCv(i,j) = 0.0
    endif
    if (G%dyCv(i,J) <= 0.0) then
      write(warnmesg,68)  pe_here(),"dyCv",i,j,G%dyCv(i,j),0.0
      call SIS_mesg(warnmesg, all_print=.true.)
      G%dyCv(i,J) = 0.0
    endif
    G%IdxCv(i,J) = Adcroft_reciprocal(G%dxCv(i,J))
    G%IdyCv(i,J) = Adcroft_reciprocal(G%dyCv(i,J))
    G%areaCv(i,J) = G%dyCv(i,J)*G%dxCv(i,J)  !### Replace with * G%dx_Cv(i,J)?
    G%IareaCv(i,J) = G%mask2dCv(i,J) * Adcroft_reciprocal(G%areaCv(i,J))
  enddo ; enddo

  do J=JsdB,JedB ; do I=IsdB,IedB
    if (G%dxBu(I,J) < 0.0) then
      write(warnmesg,68)  pe_here(),"dxBu",I,J,G%dxBu(I,J),0.0
      call SIS_mesg(warnmesg, all_print=.true.)
      G%dxBu(I,J) = 0.0
    endif
    if (G%dyBu(I,J) < 0.0) then
      write(warnmesg,68)  pe_here(),"dyBu",I,J,G%dyBu(I,J),0.0
      call SIS_mesg(warnmesg, all_print=.true.)
      G%dyBu(I,J) = 0.0
    endif

    G%IdxBu(I,J) = Adcroft_reciprocal(G%dxBu(I,J))
    G%IdyBu(I,J) = Adcroft_reciprocal(G%dyBu(I,J))
    ! G%areaBu(I,J) = G%dxBu(I,J) * G%dyBu(I,J)
    G%IareaBu(I,J) = Adcroft_reciprocal(G%areaBu(I,J))
  enddo ; enddo

68 FORMAT ("WARNING: PE ",I4," ",a3,"(",I4,",",I4,") = ",ES10.4, &
           " is being changed to ",ES10.4,".")

end subroutine set_grid_derived_metrics

function Adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val
  real :: I_val
  ! This function implements Adcroft's rule for division by 0.  

  I_val = 0.0
  if (val /= 0.0) I_val = 1.0/val
end function Adcroft_reciprocal

subroutine set_grid_metrics_from_mosaic(G,param_file)
  type(sea_ice_grid_type), intent(inout) :: G
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
  type(SIS_domain_type) :: SGdom ! Supergrid domain
  integer :: i, j, i2, j2
  integer :: npei,npej
  integer, dimension(:), allocatable :: exni,exnj
  type(domain1D) :: domx, domy
  integer        :: start(4), nread(4)
 
  call SIS_mesg("   MOM_grid_init.F90, set_grid_metrics_from_mosaic: reading grid", 5)

  call get_param(param_file, mod, "GRID_FILE", grid_file, &
                 "Name of the file from which to read horizontal grid data.", &
                 fail_if_missing=.true.)
  call get_param(param_file,  mod, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  filename = trim(adjustl(inputdir)) // trim(adjustl(grid_file))
  call log_param(param_file, mod, "INPUTDIR/GRID_FILE", filename)
  if (.not.file_exist(filename)) &
    call SIS_error(FATAL," set_grid_metrics_from_mosaic: Unable to open "//&
                           trim(filename))

! Initialize everything to a small number
  dxCu(:,:) = 0.0 ; dyCu(:,:) = 0.0
  dxCv(:,:) = 0.0 ; dyCv(:,:) = 0.0
  dxBu(:,:) = 0.0 ; dyBu(:,:) = 0.0 ; areaBu(:,:) = 0.0

!<MISSING CODE TO READ REFINEMENT LEVEL>
  ni = 2*(G%iec-G%isc+1) ! i size of supergrid
  nj = 2*(G%jec-G%jsc+1) ! j size of supergrid

! Define a domain for the supergrid (SGdom)
  npei = G%domain%layout(1) ; npej = G%domain%layout(2)
  allocate(exni(npei)) ; allocate(exnj(npej))
  call mpp_get_domain_components(G%domain%mpp_domain, domx, domy)
  call mpp_get_compute_domains(domx, size=exni)
  call mpp_get_compute_domains(domy, size=exnj)
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
     call mpp_define_domains(global_indices, SGdom%layout, SGdom%mpp_domain, &
            xflags=G%domain%X_FLAGS, yflags=G%domain%Y_FLAGS, &
            xhalo=SGdom%nihalo, yhalo=SGdom%njhalo, &
            xextent=exni,yextent=exnj, &
            symmetry=.true., name="MOM_MOSAIC", maskmap=G%domain%maskmap)
  else
     call mpp_define_domains(global_indices, SGdom%layout, SGdom%mpp_domain, &
            xflags=G%domain%X_FLAGS, yflags=G%domain%Y_FLAGS, &
            xhalo=SGdom%nihalo, yhalo=SGdom%njhalo, &
            xextent=exni,yextent=exnj, &
            symmetry=.true., name="MOM_MOSAIC")
  endif

  if (SGdom%use_io_layout) &
    call mpp_define_IO_domain(SGdom%mpp_domain, SGdom%io_layout)
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
  
  G%gridLonT(:) = tmpGlbl(2:ni:2,2)
  G%gridLonB(:) = tmpGlbl(1:ni+1:2,1)
  deallocate( tmpGlbl )

  allocate  ( tmpGlbl(1, nj+1) )  
  start(:) = 1 ; nread(:) = 1
  start(1) = int(ni/4)+1 ; nread(2) = nj+1  
  if (is_root_PE()) &
    call read_data(filename, "y", tmpGlbl, start, nread, no_domain=.TRUE.)
  call broadcast(tmpGlbl, nj+1, root_PE())

  G%gridLatT(:) = tmpGlbl(1,2:nj:2)
  G%gridLatB(:) = tmpGlbl(1,1:nj+1:2)
  deallocate( tmpGlbl )

end subroutine set_grid_metrics_from_mosaic

! ------------------------------------------------------------------------------

subroutine extrapolate_metric(var, jh, missing)
  real, dimension(:,:), intent(inout) ::  var
  integer, intent(in) :: jh
  real, optional, intent(in) :: missing
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



!---------------------------------------------------------------------

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

  allocate(G%sin_rot(isd:ied,jsd:jed)) ; G%sin_rot(:,:) = 0.0
  allocate(G%cos_rot(isd:ied,jsd:jed)) ; G%cos_rot(:,:) = 1.0

  allocate(G%H_cat_lim(1:G%CatIce+1)) ; G%H_cat_lim(:) = 0.0
  allocate(G%M_cat_lim(1:G%CatIce+1)) ; G%M_cat_lim(:) = 0.0

  allocate(G%gridLonT(isg:ieg))   ; G%gridLonT(:) = 0.0
  allocate(G%gridLonB(isg-1:ieg)) ; G%gridLonB(:) = 0.0
  allocate(G%gridLatT(jsg:jeg))   ; G%gridLatT(:) = 0.0
  allocate(G%gridLatB(jsg-1:jeg)) ; G%gridLatB(:) = 0.0

end subroutine allocate_metrics

!> Returns true if the coordinates (x,y) are within the h-cell (i,j)
logical function isPointInCell(G, i, j, x, y)
  type(sea_ice_grid_type),   intent(in) :: G    !< Grid type
  integer,                   intent(in) :: i, j !< i,j indices of cell to test
  real,                      intent(in) :: x, y !< x,y coordinates of point
! This is a crude calculation that assume a geographic coordinate system
  real :: xNE, xNW, xSE, xSW, yNE, yNW, ySE, ySW
  real :: p0, p1, p2, p3, l0, l1, l2, l3
  isPointInCell = .false.
  xNE = G%geoLonBu(i  ,j  ); yNE = G%geoLatBu(i  ,j  )
  xNW = G%geoLonBu(i-1,j  ); yNW = G%geoLatBu(i-1,j  )
  xSE = G%geoLonBu(i  ,j-1); ySE = G%geoLatBu(i  ,j-1)
  xSW = G%geoLonBu(i-1,j-1); ySW = G%geoLatBu(i-1,j-1)
  if (x<min(xNE,xNW,xSE,xSW) .or. x>max(xNE,xNW,xSE,xSW) .or. &
      y<min(yNE,yNW,ySE,ySW) .or. y>max(yNE,yNW,ySE,ySW) ) then
    return ! Avoid the more complicated calculation
  endif
  l0=(x-xSW)*(ySE-ySW)-(y-ySW)*(xSE-xSW)
  l1=(x-xSE)*(yNE-ySE)-(y-ySE)*(xNE-xSE)
  l2=(x-xNE)*(yNW-yNE)-(y-yNE)*(xNW-xNE)
  l3=(x-xNW)*(ySW-yNW)-(y-yNW)*(xSW-xNW)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.

  if ( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) .eq. abs((p0+p2)+(p1+p3)) ) then
    isPointInCell=.true.
  endif
end function isPointInCell

!---------------------------------------------------------------------
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

  deallocate(cell_area)
  
  deallocate(G%Domain%mpp_domain)
  deallocate(G%Domain)

end subroutine ice_grid_end

end module ice_grid_mod
