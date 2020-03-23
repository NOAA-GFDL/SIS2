!>  This module contains the SIS_tracer_registry_type and subroutines that handle the registration
!! of tracers and related subroutines.
module SIS_tracer_registry

! This file is a part of SIS2.  See LICENSE.md for the license.

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!* By Robert Hallberg, May 2013                                        *
!*                                                                     *
!*   This module contains the SIS_tracer_registry_type and subroutines *
!* that handle the registration of tracers and related subroutines.    *
!* The primary subroutine, register_SIS_tracer, is called to indicate  *
!* the tracers that will be advected around with the sea-ice.  This    *
!* code was derived from its MOM6 counterpart, MOM_tracer_registry.F90 *
!*                                                                     *
!* Note by ashao (2016): This is a relatively low level module which   *
!* may not need to be modified for those seeking to add a tracer to    *
!* SIS2. Most users should look at SIS_tracer_flow_control.F90 and     *
!* ice_age_tracer.F90 for examples on which to model their own         *
!* tracer implementation.                                              *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use SIS_debugging,     only : hchksum
use SIS_diag_mediator, only : SIS_diag_ctrl
use MOM_domains,       only : pass_var, pe_here
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING
use MOM_error_handler, only : SIS_mesg=>MOM_mesg
use MOM_file_parser, only : get_param, log_version, param_file_type
use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type

implicit none ; private

#include <SIS2_memory.h>

public register_SIS_tracer, register_SIS_tracer_pair, get_SIS_tracer_pointer
public SIS_unpack_passive_ice_tr, SIS_repack_passive_ice_tr, SIS_count_passive_tracers
public SIS_tracer_chksum, add_SIS_tracer_diagnostics, add_SIS_tracer_OBC_values
public update_SIS_tracer_halos, set_massless_SIS_tracers, check_SIS_tracer_bounds
public SIS_tracer_registry_init, SIS_tracer_registry_end

!> A type that describes a SIS tracer field
type, public :: SIS_tracer_type
  real, dimension(:,:,:,:), &
    pointer :: t => NULL()   !< The array containing the tracer concentration, with dimensions
                             !! of x-, y-, category, and layer.
  integer :: nL = 0          !< The number of vertical layers for this tracer.
  real :: massless_val = 0.0 !< A value to use in massless layers.
  real, dimension(:,:), &
    pointer :: ad2d_x => NULL() !< The x-direction advective flux summed vertically and across
                                !! ice category [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:), &
    pointer :: ad2d_y => NULL() !< The y-direction advective flux summed vertically and across
                                !! ice category [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:), &
    pointer :: ad3d_x => NULL() !< The vertically summed x-direction advective flux [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:), &
    pointer :: ad3d_y => NULL() !< The vertically summed y-direction advective flux [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:,:), &
    pointer :: ad4d_x => NULL() !< The x-direction advective flux by ice category and layer in
                                !! units of [Conc R Z L2 T-1 ~> CONC kg s-1].
  real, dimension(:,:,:,:), &
    pointer :: ad4d_y => NULL() !< The y-direction advective flux by ice category and layer in
                                !! units of [Conc R Z L2 T-1 ~> CONC kg s-1].
!  real, dimension(:,:), &
!    pointer :: snow_flux_tr => NULL() !< Concentration of the tracer in snow (for salinity = 0.0)
  real, dimension(:,:,:), &
    pointer :: ocean_BC => NULL() !< Value of the tracer at the ice-ocean boundary
  real, dimension(:,:,:), &
    pointer :: snow_BC => NULL() !< Value of the tracer at the snow-ice boundary

  ! @ashao: OBC NOT IMPLEMENTED YET
  real :: OBC_inflow_conc = 0.0 !< A tracer concentration for generic inflows.
  real, dimension(:,:,:), &
    pointer :: OBC_in_u => NULL()  !< Structured values for flow into the domain through u-face open
                                   !! boundaries.
  real, dimension(:,:,:), &
    pointer :: OBC_in_v => NULL()  !< Structured values for flow into the domain through v-face open
                                   !! boundaries.
  character(len=32) :: name        !< A tracer name for error messages.
  logical :: nonnegative = .false. !< If true, this tracer can not be negative
  logical :: is_passive  = .false. !< True if this ice tracer is passive
end type SIS_tracer_type

!> A type that can be used to create arrays of pointers
type, public :: p3d
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3-d array
end type p3d

!> The SIS tracer registry
type, public :: SIS_tracer_registry_type
  integer :: ntr = 0         !< The number of registered tracers.
  type(SIS_tracer_type) :: Tr_snow(MAX_FIELDS_) !< The array of registered snow tracers.
  type(SIS_tracer_type) :: Tr_ice(MAX_FIELDS_)  !< The array of registered ice tracers.
  type(SIS_diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                             !! timing of diagnostic output.
end type SIS_tracer_registry_type

contains

!> Register an ice or snow tracer to be advected
subroutine register_SIS_tracer(tr1, G, IG, nLtr, name, param_file, TrReg, snow_tracer, &
                             massless_val, ad_2d_x, ad_2d_y, ad_3d_x, ad_3d_y, &
                             ad_4d_x, ad_4d_y, OBC_inflow, OBC_in_u, OBC_in_v, &
                             nonnegative, ocean_BC, snow_BC, is_passive)
  integer,                 intent(in) :: nLtr !<  The number of vertical levels for this tracer
  type(SIS_hor_grid_type), intent(in) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG),nLtr), &
                   target, intent(in) :: tr1 !< The pointer to the tracer, in arbitrary
                                             !! concentration units [Conc], and dimensions of
                                             !! i-, j-, category, and layer.
  character(len=*),        intent(in) :: name  !< The name to be used in messages about the tracer.
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS_tracer_registry_type), &
                           pointer    :: TrReg !< A pointer to the SIS tracer registry
  logical,       optional, intent(in) :: snow_tracer !< If present and true, this is a snow tracer.
  real,          optional, intent(in) :: massless_val !< The value to use to fill in massless categories.
  real, dimension(:,:), &
                 optional, pointer    :: ad_2d_x !< An array for the x-direction advective flux summed
                                             !! vertically and across ice category [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:), &
                 optional, pointer    :: ad_2d_y !< An array for the Y-direction advective flux summed
                                             !! vertically and across ice category [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:), &
                 optional, pointer    :: ad_3d_x !< An array for the vertically summed x-direction
                                             !! advective flux [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:), &
                 optional, pointer    :: ad_3d_y !< An array for the vertically summed y-direction
                                             !! advective flux [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:,:), &
                 optional, pointer    :: ad_4d_x !< An array for the x-direction advective flux by
                                             !! ice category and layer [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:,:), &
                 optional, pointer    :: ad_4d_y  !< An array for the x-direction advective flux by
                                             !! ice category and layer [Conc R Z L2 T-1 ~> Conc kg s-1].
  real,          optional, intent(in) :: OBC_inflow !<  The value of the tracer for all inflows via
                                             !! the open boundary conditions for which OBC_in_u or
                                             !! OBC_in_v are not specified, in the same units as tr [Conc].
  real, dimension(:,:,:), &
                 optional, pointer    :: OBC_in_u !< The value of the tracer at inflows through u-faces
                                             !! of tracer cells, in the same units as tr [Conc].
  real, dimension(:,:,:), &
                 optional, pointer    :: OBC_in_v !< The value of the tracer at inflows through v-faces
                                             !! of tracer cells, in the same units as tr [Conc].
  logical,       optional, intent(in) :: nonnegative !< If true, this tracer should never be negative.
  real, dimension(:,:,:), &
                 optional, pointer    :: ocean_BC !< Value of the tracer at the ice-ocean boundary
  real, dimension(:,:,:), &
                 optional, pointer    :: snow_BC  !< Value of the tracer at the snow-ice boundary
  logical,       optional, intent(in) :: is_passive  !< True if this ice tracer is passive
  ! This subroutine registers a tracer to be advected.

  ! Local variables
  logical :: snow_tr
  type(SIS_tracer_type), pointer :: Tr_here=>NULL()
  character(len=256) :: mesg    ! Message for error messages.

  if (.not. associated(TrReg)) call SIS_tracer_registry_init(param_file, TrReg)

  snow_tr = .false. ; if (present(snow_tracer)) snow_tr = snow_tracer

  if (TrReg%ntr>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in SIS_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_SIS_tracer.")') TrReg%ntr+1
    call SIS_error(FATAL,"MOM register_SIS_tracer: "//mesg)
  endif

  TrReg%ntr = TrReg%ntr + 1
  if (snow_tr) then
    Tr_here => TrReg%Tr_snow(TrReg%ntr)
  else
    Tr_here => TrReg%Tr_ice(TrReg%ntr)
  endif

  Tr_here%name = trim(name)
  Tr_here%t => tr1(:,:,:,1:nLtr)
  Tr_here%nL = nLtr
  if(present(ocean_BC)) then
    TrReg%Tr_ice(TrReg%ntr)%ocean_BC => ocean_BC
  endif
  if(present(snow_BC)) then
    TrReg%Tr_ice(TrReg%ntr)%snow_BC => snow_BC
  endif

  if (present(massless_val)) Tr_here%massless_val = massless_val

  if (present(ad_2d_x)) then ; if (associated(ad_2d_x)) Tr_here%ad2d_x => ad_2d_x ; endif
  if (present(ad_2d_y)) then ; if (associated(ad_2d_y)) Tr_here%ad2d_y => ad_2d_y ; endif
  if (present(ad_3d_x)) then ; if (associated(ad_3d_x)) Tr_here%ad3d_x => ad_3d_x ; endif
  if (present(ad_3d_y)) then ; if (associated(ad_3d_y)) Tr_here%ad3d_y => ad_3d_y ; endif
  if (present(ad_4d_x)) then ; if (associated(ad_4d_x)) then
    if (size(ad_4d_x,4) /= Tr_here%nL) call SIS_error(FATAL, &
         "Mismatch in register_SIS_tracer of the number of vertical levels "//&
         "in ad_4d_x for "//trim(Tr_here%name))
    Tr_here%ad4d_x => ad_4d_x
  endif ; endif
  if (present(ad_4d_y)) then ; if (associated(ad_4d_y)) then
    if (size(ad_4d_y,4) /= Tr_here%nL) call SIS_error(FATAL, &
         "Mismatch in register_SIS_tracer of the number of vertical levels "//&
         "in ad_4d_4 for "//trim(Tr_here%name))
    Tr_here%ad4d_y => ad_4d_y
  endif ; endif
  if (present(OBC_inflow)) Tr_here%OBC_inflow_conc = OBC_inflow
  if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
                                    Tr_here%OBC_in_u => OBC_in_u ; endif
  if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
                                    Tr_here%OBC_in_v => OBC_in_v ; endif
  Tr_here%nonnegative = .false.
  if (present(nonnegative)) Tr_here%nonnegative = nonnegative

  ! This is set as the default to guard against operations being performed on
  ! active tracers twice
  TrReg%Tr_ice(TrReg%ntr)%is_passive = .false.
  if(present(is_passive)) TrReg%Tr_ice(TrReg%ntr)%is_passive = is_passive

end subroutine register_SIS_tracer

!> Register a pair of ice and snow tracers to be advected
subroutine register_SIS_tracer_pair(ice_tr, nL_ice, name_ice, snow_tr, nL_snow, &
                                    name_snow, G, IG, param_file, TrReg, &
                                    massless_iceval, massless_snowval, nonnegative, is_passive)
  integer,                 intent(in) :: nL_ice !<  The number of vertical levels for the ice tracer.
  integer,                 intent(in) :: nL_snow !< The number of vertical levels for the snow tracer.
  type(SIS_hor_grid_type), intent(in) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG),nL_ice), &
                   target, intent(in) :: ice_tr !<  The pointer to the ice tracer, in arbitrary
                                               !! concentration units [Conc], and dimensions
                                               !! of i-, j-, category, and layer.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG),nL_snow), &
                   target, intent(in) :: snow_tr !< The pointer to the snow tracer, in arbitrary
                                               !! concentration units [Conc], and dimensions
                                               !! of i-, j-, category, and layer.
  character(len=*),        intent(in) :: name_ice !< The name to be used in messages about the tracer.
  character(len=*),        intent(in) :: name_snow !< The name to be used in messages about the tracer.
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS_tracer_registry_type), &
                           pointer    :: TrReg !< A pointer to the SIS tracer registry
  real,          optional, intent(in) :: massless_iceval !< The values to use to fill in
                                               !! massless ice categories.
  real,          optional, intent(in) :: massless_snowval !< The values to use to fill in
                                               !! massless snow categories.
  logical,       optional, intent(in) :: nonnegative !< If true, this tracer should never be negative.
  logical,       optional, intent(in) :: is_passive  !< True if this ice tracer is passive
! This subroutine registers a pair of ice and snow tracers to be advected.

  character(len=256) :: mesg    ! Message for error messages.

  if (.not. associated(TrReg)) call SIS_tracer_registry_init(param_file, TrReg)

  if (TrReg%ntr>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in SIS_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_SIS_tracer.")') TrReg%ntr+1
    call SIS_error(FATAL,"MOM register_SIS_tracer: "//mesg)
  endif

  TrReg%ntr = TrReg%ntr + 1
  TrReg%Tr_ice(TrReg%ntr)%name = trim(name_ice)
  TrReg%Tr_ice(TrReg%ntr)%t => ice_tr(:,:,:,1:nL_ice)
  TrReg%Tr_ice(TrReg%ntr)%nL = nL_ice

  if (present(massless_iceval)) &
    TrReg%Tr_ice(TrReg%ntr)%massless_val = massless_iceval

  TrReg%Tr_snow(TrReg%ntr)%name = trim(name_snow)
  TrReg%Tr_snow(TrReg%ntr)%t => snow_tr(:,:,:,1:nL_snow)
  TrReg%Tr_snow(TrReg%ntr)%nL = nL_snow

  if (present(massless_snowval)) &
    TrReg%Tr_snow(TrReg%ntr)%massless_val = massless_snowval

  TrReg%Tr_ice(TrReg%ntr)%nonnegative = .false.
  if (present(nonnegative)) TrReg%Tr_ice(TrReg%ntr)%nonnegative = nonnegative

  TrReg%Tr_ice(TrReg%ntr)%is_passive = .false.
  if (present(is_passive)) TrReg%Tr_ice(TrReg%ntr)%is_passive = is_passive

end subroutine register_SIS_tracer_pair

!> Unpacks only the passive tracer arrays into TrLay which is a slice of the
!! vertical layers within an ice thickness category
subroutine SIS_unpack_passive_ice_tr(i, j, cat, nkice, TrReg, TrLay)
  integer,                        intent(in)    :: i     !< Horizontal index i-direction
  integer,                        intent(in)    :: j     !< Horizontal index j-direction
  integer,                        intent(in)    :: cat   !< Index of ice thickness category
  integer,                        intent(in)    :: nkice !< Number of levels in the ice
  type(SIS_tracer_registry_type), intent(in)    :: TrReg !< Main tracer registry
  real, dimension(0:,:),          intent(inout) :: TrLay !< Array to hold vertical slice
                                                         !! of passive tracers

  integer :: m, tr, pass_idx

  pass_idx = 0
  do tr=1,TrReg%ntr
    if(TrReg%Tr_ice(tr)%is_passive) then
      pass_idx = pass_idx + 1
      ! Copy from main tracer array
      do m=1,nkice ; TrLay(m,pass_idx) = TrReg%Tr_ice(tr)%t(i,j,cat,m) ; enddo

      ! Set snow and ice boundary conditions (if they exist)
      if(associated(TrReg%Tr_ice(tr)%ocean_BC)) then
        TrLay(NkIce+1,pass_idx) = TrReg%Tr_ice(tr)%ocean_BC(i,j,cat)
      else
        TrLay(NkIce+1,pass_idx) = 0.0
      endif

      if(associated(TrReg%Tr_ice(tr)%snow_BC)) then
        TrLay(0,pass_idx) = TrReg%Tr_ice(tr)%snow_BC(i,j,cat)
      else
        TrLay(0,pass_idx) = 0.0
      endif
    endif
  enddo

end subroutine SIS_unpack_passive_ice_tr

!> Copy the vertical slice of passive tracers back into their 4D arrays
subroutine SIS_repack_passive_ice_tr(i, j, cat, nkice, TrReg, TrLay)
  integer,                        intent(in   ) :: i     !< Horizontal index i-direction
  integer,                        intent(in   ) :: j     !< Horizontal index j-direction
  integer,                        intent(in   ) :: cat   !< Index of ice thickness category
  integer,                        intent(in   ) :: nkice !< Number of levels in the ice
  type(SIS_tracer_registry_type), intent(inout) :: TrReg !< Main tracer registry
  real, dimension(0:,:),          intent(inout) :: TrLay !< Array to hold vertical slice
                                                         !! of passive tracers

  integer :: m, tr, pass_idx

  pass_idx = 0
  do tr=1,TrReg%ntr
    if(TrReg%Tr_ice(tr)%is_passive) then
      pass_idx = pass_idx + 1
      ! Copy values back to tracer array
      do m=1,nkice ; TrReg%Tr_ice(tr)%t(i,j,cat,m) = TrLay(m,pass_idx) ; enddo
    endif
  enddo

end subroutine SIS_repack_passive_ice_tr

!> Returns the total number of passive ice tracers
integer function SIS_count_passive_tracers(TrReg)
  type(SIS_tracer_registry_type), intent(in) :: TrReg !< A pointer to the SIS tracer registry

  integer :: tr

  SIS_count_passive_tracers = 0
  do tr=1,TrReg%ntr
    if( TrReg%Tr_ice(tr)%is_passive ) then
      SIS_count_passive_tracers = SIS_count_passive_tracers + 1
    endif
  enddo

end function SIS_count_passive_tracers

!> Set a pointer to a named SIS tracer
subroutine get_SIS_tracer_pointer(name, TrReg, Tr_ptr, nLayer)
  character(len=*),               intent(in) :: name  !< The tracer name being queried for
  type(SIS_tracer_registry_type), intent(in) :: TrReg !< A pointer to the SIS tracer registry
  real, dimension(:,:,:,:),       pointer    :: Tr_ptr !< A pointer to the named tracer
  integer,                        intent(out):: nLayer !< The number of layers for this tracer

  integer :: m

  do m=1,TrReg%ntr ; if (TrReg%Tr_ice(m)%name == trim(name)) exit ; enddo
  if (m <= TrReg%ntr) then ! This is an ice tracer.
    Tr_ptr => TrReg%Tr_ice(m)%t
    nLayer = TrReg%Tr_ice(m)%nL
  else  ! See whether this is a snow tracer.
    do m=1,TrReg%ntr ; if (TrReg%Tr_snow(m)%name == trim(name)) exit ; enddo
    if (m <= TrReg%ntr) then ! This is a snow tracer.
      Tr_ptr => TrReg%Tr_snow(m)%t
      nLayer = TrReg%Tr_snow(m)%nL
    else
      call SIS_error(FATAL, "SIS_tracer: register_SIS_tracer must be called for "//&
               trim(name)//" before get_registered_tracer_pointer is called for it.")
    endif
  endif

end subroutine get_SIS_tracer_pointer

!> Update the halos of the SIS tracers
subroutine update_SIS_tracer_halos(TrReg, G, complete)
  type(SIS_tracer_registry_type), intent(inout) :: TrReg !< A pointer to the SIS tracer registry
  type(SIS_hor_grid_type),        intent(inout) :: G   !< The horizontal grid type
  logical,              optional, intent(in)    :: complete !< If true, complete this update

  logical :: do_complete, comp_here
  integer :: m, n, last_ice_tr

  do_complete = .true. ; if (present(complete)) do_complete=complete
  if (TrReg%ntr==0) return

  last_ice_tr = -1
  do n=TrReg%ntr,1,-1 ; if (TrReg%Tr_ice(n)%nL > 0) then
    last_ice_tr = n ; exit
  endif ; enddo
  if (last_ice_tr < 0) &
    call SIS_error(FATAL, "At least 1 ice tracer needs to be registered for "//&
                   "update_SIS_tracer_halos to work as currently written.")

  do n=1,TrReg%ntr ; do m=1,TrReg%Tr_snow(n)%nL
    call pass_var(TrReg%Tr_snow(n)%t(:,:,:,m), G%Domain, complete=.false.)
  enddo ; enddo
  do n=1,TrReg%ntr ; do m=1,TrReg%Tr_ice(n)%nL
    comp_here = (do_complete .and. (n==last_ice_tr) .and. (m==TrReg%Tr_ice(n)%nL))
    call pass_var(TrReg%Tr_ice(n)%t(:,:,:,m), G%Domain, complete=comp_here)
  enddo ; enddo

end subroutine update_SIS_tracer_halos

!> Set the properties of the ice tracers where there is no ice mass
subroutine set_massless_SIS_tracers(mass, TrReg, G, IG, compute_domain, do_snow, do_ice)
  type(SIS_hor_grid_type),        intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),            intent(inout) :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                                  intent(in)    :: mass !< The ice or snow mass [R Z ~> kg m-2].
  type(SIS_tracer_registry_type), intent(inout) :: TrReg !< A pointer to the SIS tracer registry
  logical,              optional, intent(in) :: compute_domain !< If true, work over the whole data domain
  logical,              optional, intent(in) :: do_snow !< If true, work on snow tracers; the default is true.
  logical,              optional, intent(in) :: do_ice  !< If true, work on ice tracers; the default is true.

  integer :: i, j, k, m, n, is, ie, js, je, nCat
  logical :: do_snow_tr, do_ice_tr

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nCat = IG%CatIce
  if (present(compute_domain)) then ; if (compute_domain) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif

  do_snow_tr = .true. ; do_ice_tr = .true.
  if (present(do_snow)) do_snow_tr = do_snow
  if (present(do_ice))  do_ice_tr = do_ice

  if (do_snow_tr) then
    do n=1,TrReg%ntr ; do m=1,TrReg%Tr_snow(n)%nL
      do k=1,nCat ; do j=js,je ; do i=is,ie ; if (mass(i,j,k)<=0.0) &
        TrReg%Tr_snow(n)%t(i,j,k,m) = TrReg%Tr_snow(n)%massless_val
      enddo ; enddo ; enddo
    enddo ; enddo
  endif
  if (do_ice_tr) then
    do n=1,TrReg%ntr ; do m=1,TrReg%Tr_ice(n)%nL
      do k=1,nCat ; do j=js,je ; do i=is,ie ; if (mass(i,j,k)<=0.0) &
        TrReg%Tr_ice(n)%t(i,j,k,m) = TrReg%Tr_ice(n)%massless_val
      enddo ; enddo ; enddo
    enddo ; enddo
  endif

end subroutine set_massless_SIS_tracers

!> Check the ice tracers for unphysical values
subroutine check_SIS_tracer_bounds(TrReg, G, IG, msg, data_domain)
  type(SIS_hor_grid_type),        intent(in)    :: G   !< The horizontal grid type
  type(ice_grid_type),            intent(in)    :: IG  !< The sea-ice specific grid type
  character(len=*),               intent(in)    :: msg !< An identifying message fragment
  type(SIS_tracer_registry_type), intent(inout) :: TrReg !< A pointer to the SIS tracer registry
  logical,              optional, intent(in)    :: data_domain !< If true, work over the whole data domain

  character(len=512) :: mesg1, mesg2
  integer :: i, j, k, m, n, is, ie, js, je, nCat

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nCat = IG%CatIce
  if (present(data_domain)) then ; if (data_domain) then
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif ; endif

  do n=1,TrReg%ntr ; if (TrReg%Tr_ice(n)%nonnegative) then
    do m=1,TrReg%Tr_ice(n)%nL ; do k=1,nCat ; do j=js,je ; do i=is,ie
      if (TrReg%Tr_ice(n)%t(i,j,k,m) < 0.0) then
        write(mesg1,'(" = ",1pe12.4," at ", 2(F7.1)," or i,j,k,m = ",4i4,"; on pe ",i4)') &
           TrReg%Tr_ice(n)%t(i,j,k,m), G%geolonT(i,j), G%geolatT(i,j), i, j, k, m, pe_here()
        call SIS_error(WARNING, trim(msg)//": Negative ice tracer "//&
                       trim(TrReg%Tr_ice(n)%name)// trim(mesg1), all_print=.true.)
        TrReg%Tr_ice(n)%t(i,j,k,m) = 0.0
      endif
    enddo ; enddo ; enddo ; enddo
  endif ; enddo
  do n=1,TrReg%ntr ; if (TrReg%Tr_snow(n)%nonnegative) then
    do m=1,TrReg%Tr_snow(n)%nL ; do k=1,nCat ; do j=js,je ; do i=is,ie
      if (TrReg%Tr_snow(n)%t(i,j,k,m) < 0.0) then
        write(mesg1,'(" = ",1pe12.4," at ", 2(F7.1)," or i,j,k,m = ",4i4,"; on pe ",i4)') &
           TrReg%Tr_snow(n)%t(i,j,k,m), G%geolonT(i,j), G%geolatT(i,j), i, j, k, m, pe_here()
        call SIS_error(WARNING, trim(msg)//": Negative snow tracer "//&
                       trim(TrReg%Tr_snow(n)%name)// trim(mesg1), all_print=.true.)
        TrReg%Tr_snow(n)%t(i,j,k,m) = 0.0
      endif
    enddo ; enddo ; enddo ; enddo
  endif ; enddo

end subroutine check_SIS_tracer_bounds

!> Add open boundary condition concentrations for a tracer that has previously
!! been registered by a call to register_SIS_tracer.
subroutine add_SIS_tracer_OBC_values(name, TrReg, OBC_inflow, OBC_in_u, OBC_in_v)
  character(len=*),        intent(in) :: name  !< The tracer name whose properties are being set
  type(SIS_tracer_registry_type), &
                           pointer    :: TrReg !< A pointer to the SIS tracer registry
  real,          optional, intent(in) :: OBC_inflow !<  The value of the tracer for all inflows via
                                             !! the open boundary conditions for which OBC_in_u or
                                             !! OBC_in_v are not specified, in the same units as tr [Conc].
  real, dimension(:,:,:), &
                 optional, pointer    :: OBC_in_u !< The value of the tracer at inflows through u-faces
                                             !! of tracer cells, in the same units as tr [Conc].
  real, dimension(:,:,:), &
                 optional, pointer    :: OBC_in_v !< The value of the tracer at inflows through v-faces
                                             !! of tracer cells, in the same units as tr [Conc].
! This subroutine adds open boundary condition concentrations for a tracer that
! has previously been registered by a call to register_SIS_tracer.

  type(SIS_tracer_type), pointer :: Tr_here=>NULL()
  integer :: m

  if (.not. associated(TrReg)) call SIS_error(FATAL, "add_SIS_tracer_OBC_values :"// &
       "register_SIS_tracer must be called before add_SIS_tracer_OBC_values")

  do m=1,TrReg%ntr ; if (TrReg%Tr_ice(m)%name == trim(name)) exit ; enddo
  if (m <= TrReg%ntr) then ; Tr_here => TrReg%Tr_ice(m)
  else  ! See whether this is a snow tracer.
    do m=1,TrReg%ntr ; if (TrReg%Tr_snow(m)%name == trim(name)) exit ; enddo
    if (m <= TrReg%ntr) then ; Tr_here => TrReg%Tr_snow(m)
    else
      call SIS_error(FATAL, "SIS_tracer: register_SIS_tracer must be called for "//&
               trim(name)//" before add_SIS_tracer_OBC_values is called for it.")
    endif
  endif

  if (present(OBC_inflow)) Tr_here%OBC_inflow_conc = OBC_inflow
  if (present(OBC_in_u)) then ; if (associated(OBC_in_u)) &
                                    Tr_here%OBC_in_u => OBC_in_u ; endif
  if (present(OBC_in_v)) then ; if (associated(OBC_in_v)) &
                                    Tr_here%OBC_in_v => OBC_in_v ; endif

end subroutine add_SIS_tracer_OBC_values

!> Add diagnostic arrays for a tracer that has previously been registered via a
!! call to register_SIS_tracer.
subroutine add_SIS_tracer_diagnostics(name, TrReg, ad_2d_x, ad_2d_y, ad_3d_x, &
                                      ad_3d_y, ad_4d_x, ad_4d_y)
  character(len=*),        intent(in) :: name  !< The tracer name whose properties are being set
  type(SIS_tracer_registry_type), &
                           pointer    :: TrReg !< A pointer to the SIS tracer registry
  real, dimension(:,:), &
                 optional, pointer    :: ad_2d_x !< An array for the x-direction advective flux summed
                                             !! vertically and across ice category [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:), &
                 optional, pointer    :: ad_2d_y !< An array for the Y-direction advective flux summed
                                             !! vertically and across ice category [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:), &
                 optional, pointer    :: ad_3d_x !< An array for the vertically summed x-direction
                                             !! advective flux [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:), &
                 optional, pointer    :: ad_3d_y !< An array for the vertically summed y-direction
                                             !! advective flux [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:,:), &
                 optional, pointer    :: ad_4d_x !< An array for the x-direction advective flux by
                                             !! ice category and layer [Conc R Z L2 T-1 ~> Conc kg s-1].
  real, dimension(:,:,:,:), &
                 optional, pointer    :: ad_4d_y  !< An array for the x-direction advective flux by
                                             !! ice category and layer [Conc R Z L2 T-1 ~> Conc kg s-1].

! This subroutine adds diagnostic arrays for a tracer that has previously been
! registered by a call to register_SIS_tracer.

  type(SIS_tracer_type), pointer :: Tr_here=>NULL()
  integer :: m

  if (.not. associated(TrReg)) call SIS_error(FATAL, "add_SIS_tracer_diagnostics: "// &
       "register_SIS_tracer must be called before add_SIS_tracer_diagnostics")

  do m=1,TrReg%ntr ; if (TrReg%Tr_ice(m)%name == trim(name)) exit ; enddo
  if (m <= TrReg%ntr) then ; Tr_here => TrReg%Tr_ice(m)
  else  ! See whether this is a snow tracer.
    do m=1,TrReg%ntr ; if (TrReg%Tr_snow(m)%name == trim(name)) exit ; enddo
    if (m <= TrReg%ntr) then ; Tr_here => TrReg%Tr_snow(m)
    else
      call SIS_error(FATAL, "SIS_tracer: register_SIS_tracer must be called for "//&
               trim(name)//" before add_SIS_tracer_OBC_values is called for it.")
    endif
  endif

  if (present(ad_2d_x)) then ; if (associated(ad_2d_x)) Tr_here%ad2d_x => ad_2d_x ; endif
  if (present(ad_2d_y)) then ; if (associated(ad_2d_y)) Tr_here%ad2d_y => ad_2d_y ; endif
  if (present(ad_3d_x)) then ; if (associated(ad_3d_x)) Tr_here%ad3d_x => ad_3d_x ; endif
  if (present(ad_3d_y)) then ; if (associated(ad_3d_y)) Tr_here%ad3d_y => ad_3d_y ; endif
  if (present(ad_4d_x)) then ; if (associated(ad_4d_x)) then
    if (size(ad_4d_x,4) /= Tr_here%nL) call SIS_error(FATAL, &
         "Mismatch in register_SIS_tracer of the number of vertical levels "//&
         "in ad_4d_x for "//trim(Tr_here%name))
    Tr_here%ad4d_x => ad_4d_x
  endif ; endif
  if (present(ad_4d_y)) then ; if (associated(ad_4d_y)) then
    if (size(ad_4d_y,4) /= Tr_here%nL) call SIS_error(FATAL, &
         "Mismatch in register_SIS_tracer of the number of vertical levels "//&
         "in ad_4d_4 for "//trim(Tr_here%name))
    Tr_here%ad4d_y => ad_4d_y
  endif ; endif

end subroutine add_SIS_tracer_diagnostics

!> Write out chksums for the sea ice tracers
subroutine SIS_tracer_chksum(mesg, TrReg, G)
  character(len=*),               intent(in) :: mesg  !< A message that appears on the chksum lines
  type(SIS_tracer_registry_type), pointer    :: TrReg !< A pointer to the SIS tracer registry
  type(SIS_hor_grid_type),        intent(in) :: G     !< The horizontal grid type.
!   This subroutine writes out chksums for the sea ice tracers

  integer :: is, ie, js, je, l, m
  character(len=8) :: mesg_l
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do m=1,TrReg%ntr ; do l=1,TrReg%Tr_ice(m)%nL
    write(mesg_l,'("i")') l
    call hchksum(TrReg%Tr_ice(m)%t(:,:,:,l), mesg//trim(TrReg%Tr_ice(m)%name)//" "//&
                                 trim(adjustl(mesg_l)), G%HI)
  enddo ; enddo
  do m=1,TrReg%ntr ; do l=1,TrReg%Tr_snow(m)%nL
    write(mesg_l,'("i")') l
    call hchksum(TrReg%Tr_snow(m)%t(:,:,:,l), mesg//trim(TrReg%Tr_snow(m)%name)//" "//&
                                 trim(adjustl(mesg_l)), G%HI)
  enddo ; enddo
end subroutine SIS_tracer_chksum

!> Allocate the SIS tracer registry
subroutine SIS_tracer_registry_init(param_file, TrReg)
  type(param_file_type),          intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS_tracer_registry_type), pointer    :: TrReg !< A pointer to the SIS tracer registry that is
                                                      !! set up here.

  integer, save :: init_calls = 0
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "SIS_tracer_registry" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(TrReg)) then ; allocate(TrReg)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  init_calls = init_calls + 1
  if (init_calls > 1) then
    write(mesg,'("SIS_tracer_registry_init called ",I3, &
      &" times with different registry pointers.")') init_calls
    call SIS_error(WARNING,"SIS_tracer"//mesg)
  endif

end subroutine SIS_tracer_registry_init

!> Deallocate memory associated with the SIS tracer registry
subroutine SIS_tracer_registry_end(TrReg)
  type(SIS_tracer_registry_type), pointer :: TrReg !< A pointer to the SIS tracer registry that is
                                                   !! deallocated here.
  if (associated(TrReg)) deallocate(TrReg)
end subroutine SIS_tracer_registry_end

end module SIS_tracer_registry
