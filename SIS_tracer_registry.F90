module SIS_tracer_registry
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of SIS2.                                        *
!*                                                                     *
!* SIS2 is free software; you can redistribute it and/or modify it and *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* SIS2 is distributed in the hope that it will be useful, but WITHOUT *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

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
!********+*********+*********+*********+*********+*********+*********+**

use SIS_diag_mediator, only : SIS_diag_ctrl
use SIS_error_checking, only : hchksum
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING
use MOM_error_handler, only : SIS_mesg=>MOM_mesg
use MOM_file_parser, only : get_param, log_version, param_file_type
use ice_grid_mod, only : sea_ice_grid_type

implicit none ; private

#include <SIS2_memory.h>
#ifdef STATIC
#  define NINTMEM_(i) i
#else
#  define NINTMEM_(i) :
#endif

public register_SIS_tracer, get_SIS_tracer_pointer, SIS_tracer_chksum
public add_SIS_tracer_diagnostics, add_SIS_tracer_OBC_values
public SIS_tracer_registry_init, SIS_tracer_registry_end

type, public :: SIS_tracer_type
  real, dimension(:,:,:,:), pointer :: t => NULL()
             ! The array containing the tracer concentration, with dimensions
             ! of x-, y-, category, and layer.
  integer :: nL = 0 ! The number of vertical layers for this tracer.
  real, dimension(:,:), pointer :: ad2d_x => NULL(), ad2d_y => NULL()
             ! The arrays in which x- & y- advective fluxes summed
             ! vertically and across ice category are stored in units of
             ! CONC m3 s-1.
  real, dimension(:,:,:), pointer :: ad3d_x => NULL(), ad3d_y => NULL()
             ! The arrays in which vertically summed x- & y- advective fluxes
             ! are stored in units of CONC m3 s-1.
  real, dimension(:,:,:,:), pointer :: ad4d_x => NULL(), ad4d_y => NULL()
             ! The arrays in which x- & y- advective fluxes by ice category and
             ! layer are stored in units of CONC m3 s-1.
  real :: OBC_inflow_conc = 0.0  ! A tracer concentration for generic inflows.
  real, dimension(:,:,:), pointer :: OBC_in_u => NULL(), OBC_in_v => NULL()
             ! These arrays contain structured values for flow into the domain
             ! that are specified in open boundary conditions through u- and
             ! v- faces of the tracer cell.
  character(len=32) :: name  ! A tracer name for error messages.
end type SIS_tracer_type

type, public :: p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: SIS_tracer_registry_type
  integer :: ntr_snow = 0        ! The number of registered tracers.
  integer :: ntr_ice  = 0        ! The number of registered tracers.
  type(SIS_tracer_type) :: Tr_snow(MAX_FIELDS_)  ! The array of registered tracers.
  type(SIS_tracer_type) :: Tr_ice(MAX_FIELDS_)  ! The array of registered tracers.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
end type SIS_tracer_registry_type

contains

subroutine register_SIS_tracer(tr1, G, nLtr, name, param_file, Reg, snow_tracer, &
                 ad_2d_x, ad_2d_y, ad_3d_x, ad_3d_y, ad_4d_x, ad_4d_y, &
                 OBC_inflow, OBC_in_u, OBC_in_v)
  integer,                                  intent(in) :: nLtr              
  type(sea_ice_grid_type),     intent(in)    :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G),nLtr), target :: tr1
  character(len=*), intent(in)                :: name
  type(param_file_type), intent(in)           :: param_file
  type(SIS_tracer_registry_type), pointer     :: Reg
  logical,               intent(in), optional :: snow_tracer
  real, dimension(:,:),     pointer, optional :: ad_2d_x, ad_2d_y
  real, dimension(:,:,:),   pointer, optional :: ad_3d_x, ad_3d_y
  real, dimension(:,:,:,:), pointer, optional :: ad_4d_x, ad_4d_y
  real, intent(in), optional                  :: OBC_inflow
  real, pointer, dimension(:,:,:), optional   :: OBC_in_u, OBC_in_v
! This subroutine registers a tracer to be advected and horizontally
! diffused.

! Arguments: tr1 - The pointer to the tracer, in arbitrary concentration units
!                  (CONC), and dimensions of i-, j-, category, and layer.
!  (in)      nLtr - The number of vertical levels for this tracer.
!  (in)      name - The name to be used in messages about the tracer.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  Reg - A pointer to the tracer registry.
!  (in,opt)  snow_tracer - If present and true, this is a snow tracer.
!  (in,opt)  ad_2d_x - An array in which the zonal advective fluxes summed
!                      vertically and across ice category are stored in units of
!                      CONC m3 s-1.
!  (in,opt)  ad_2d_y - An array in which the meridional advective fluxes summed
!                      vertically and across ice category are stored in units of
!                      CONC m3 s-1.
!  (in,opt)  ad_3d_x - An array in which the zonal advective fluxes summed
!                      vertically are stored in units of CONC m3 s-1.
!  (in,opt)  ad_3d_y - An array in which the meridional advective fluxes summed
!                      vertically are stored in units of CONC m3 s-1.
!  (in,opt)  ad_4d_x - An array in which the zonal advective fluxes by ice
!                      category and layer are stored in units of CONC m3 s-1.
!  (in,opt)  ad_4d_y - An array in which the meridional fluxes by ice
!                      category and layer are stored in units of CONC m3 s-1.
!  (in)      OBC_inflow - The value of the tracer for all inflows via the open
!                         boundary conditions for which OBC_in_u or OBC_in_v are
!                         not specified, in the same units as tr (CONC).
!  (in)      OBC_in_u - The value of the tracer at inflows through u-faces of
!                       tracer cells, in the same units as tr (CONC).
!  (in)      OBC_in_v - The value of the tracer at inflows through v-faces of
!                       tracer cells, in the same units as tr (CONC).

  integer :: ntr
  logical :: snow_tr
  type(SIS_tracer_type), pointer :: Tr_here=>NULL()
  character(len=256) :: mesg    ! Message for error messages.

  if (.not. associated(Reg)) call SIS_tracer_registry_init(param_file, Reg)

  snow_tr = .false. ; if (present(snow_tracer)) snow_tr = snow_tracer

  ntr = Reg%ntr_ice ; if (snow_tr) ntr = Reg%ntr_snow
  if (ntr>=MAX_FIELDS_) then
    write(mesg,'("Increase MAX_FIELDS_ in SIS_memory.h to at least ",I3," to allow for &
        &all the tracers being registered via register_SIS_tracer.")') ntr+1
    call SIS_error(FATAL,"MOM register_SIS_tracer: "//mesg)
  endif
  
  if (snow_tr) then
    Reg%ntr_snow = Reg%ntr_snow + 1
    Tr_here => Reg%Tr_snow(Reg%ntr_snow)
  else
    Reg%ntr_ice = Reg%ntr_ice + 1
    Tr_here => Reg%Tr_ice(Reg%ntr_ice)
  endif

  Tr_here%name = trim(name)
  Tr_here%t => tr1(:,:,:,1:nLtr)
  Tr_here%nL = nLtr

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

end subroutine register_SIS_tracer

subroutine get_SIS_tracer_pointer(name, Reg, Tr_ptr, nLayer)
  character(len=*),          intent(in)      :: name
  type(SIS_tracer_registry_type), intent(in) :: Reg
  real, dimension(:,:,:,:),       pointer    :: Tr_ptr
  integer,                    intent(out)    :: nLayer

  integer :: m

  do m=1,Reg%ntr_ice ; if (Reg%Tr_ice(m)%name == trim(name)) exit ; enddo
  if (m <= Reg%ntr_ice) then ! This is an ice tracer.
    Tr_ptr => Reg%Tr_ice(m)%t
    nLayer = Reg%Tr_ice(m)%nL
  else  ! See whether this is a snow tracer.
    do m=1,Reg%ntr_snow ; if (Reg%Tr_snow(m)%name == trim(name)) exit ; enddo
    if (m <= Reg%ntr_snow) then ! This is a snow tracer.
      Tr_ptr => Reg%Tr_snow(m)%t
      nLayer = Reg%Tr_snow(m)%nL
    else
      call SIS_error(FATAL, "SIS_tracer: register_SIS_tracer must be called for "//&
               trim(name)//" before get_registered_tracer_pointer is called for it.")
    endif
  endif

end subroutine get_SIS_tracer_pointer

subroutine add_SIS_tracer_OBC_values(name, Reg, OBC_inflow, OBC_in_u, OBC_in_v)
  character(len=*), intent(in)               :: name
  type(SIS_tracer_registry_type), pointer    :: Reg
  real, intent(in), optional                 :: OBC_inflow
  real, pointer, dimension(:,:,:), optional  :: OBC_in_u, OBC_in_v
! This subroutine adds open boundary condition concentrations for a tracer that
! has previously been registered by a call to register_SIS_tracer.

! Arguments: name - The name of the tracer for which the diagnostic pointers.
!  (in/out)  Reg - A pointer to the tracer registry.
!  (in)      OBC_inflow - The value of the tracer for all inflows via the open
!                         boundary conditions for which OBC_in_u or OBC_in_v are
!                         not specified, in the same units as tr (CONC).
!  (in)      OBC_in_u - The value of the tracer at inflows through u-faces of
!                       tracer cells, in the same units as tr (CONC).
!  (in)      OBC_in_v - The value of the tracer at inflows through v-faces of
  type(SIS_tracer_type), pointer :: Tr_here=>NULL()
  integer :: m

  if (.not. associated(Reg)) call SIS_error(FATAL, "add_SIS_tracer_OBC_values :"// &
       "register_SIS_tracer must be called before add_SIS_tracer_OBC_values")

  do m=1,Reg%ntr_ice ; if (Reg%Tr_ice(m)%name == trim(name)) exit ; enddo
  if (m <= Reg%ntr_ice) then ; Tr_here => Reg%Tr_ice(m)
  else  ! See whether this is a snow tracer.
    do m=1,Reg%ntr_snow ; if (Reg%Tr_snow(m)%name == trim(name)) exit ; enddo
    if (m <= Reg%ntr_snow) then ; Tr_here => Reg%Tr_snow(m)
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

subroutine add_SIS_tracer_diagnostics(name, Reg, ad_2d_x, ad_2d_y, ad_3d_x, &
                                      ad_3d_y, ad_4d_x, ad_4d_y)
  character(len=*), intent(in)                :: name
  type(SIS_tracer_registry_type), pointer     :: Reg
  real, dimension(:,:),     pointer, optional :: ad_2d_x, ad_2d_y
  real, dimension(:,:,:),   pointer, optional :: ad_3d_x, ad_3d_y
  real, dimension(:,:,:,:), pointer, optional :: ad_4d_x, ad_4d_y
! This subroutine adds diagnostic arrays for a tracer that has previously been
! registered by a call to register_SIS_tracer.

! Arguments: name - The name of the tracer for which the diagnostic pointers.
!  (in/out)  Reg - A pointer to the tracer registry.
!  (in,opt)  ad_2d_x - An array in which the zonal advective fluxes summed
!                      vertically and across ice category are stored in units of
!                      CONC m3 s-1.
!  (in,opt)  ad_2d_y - An array in which the meridional advective fluxes summed
!                      vertically and across ice category are stored in units of
!                      CONC m3 s-1.
!  (in,opt)  ad_3d_x - An array in which the zonal advective fluxes summed
!                      vertically are stored in units of CONC m3 s-1.
!  (in,opt)  ad_3d_y - An array in which the meridional advective fluxes summed
!                      vertically are stored in units of CONC m3 s-1.
!  (in,opt)  ad_4d_x - An array in which the zonal advective fluxes by ice
!                      category and layer are stored in units of CONC m3 s-1.
!  (in,opt)  ad_4d_y - An array in which the meridional fluxes by ice
!                      category and layer are stored in units of CONC m3 s-1.
  type(SIS_tracer_type), pointer :: Tr_here=>NULL()
  integer :: m

  if (.not. associated(Reg)) call SIS_error(FATAL, "add_SIS_tracer_diagnostics: "// &
       "register_SIS_tracer must be called before add_SIS_tracer_diagnostics")

  do m=1,Reg%ntr_ice ; if (Reg%Tr_ice(m)%name == trim(name)) exit ; enddo
  if (m <= Reg%ntr_ice) then ; Tr_here => Reg%Tr_ice(m)
  else  ! See whether this is a snow tracer.
    do m=1,Reg%ntr_snow ; if (Reg%Tr_snow(m)%name == trim(name)) exit ; enddo
    if (m <= Reg%ntr_snow) then ; Tr_here => Reg%Tr_snow(m)
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

subroutine SIS_tracer_chksum(mesg, Reg, G)
  character(len=*),               intent(in) :: mesg
  type(SIS_tracer_registry_type), pointer    :: Reg
  type(sea_ice_grid_type),        intent(in) :: G
!   This subroutine writes out chksums for the model's thermodynamic state
! variables.
! Arguments: mesg - A message that appears on the chksum lines.
!  (in)      Reg - A pointer to the tracer registry.
!  (in)      G - The ocean's grid structure.
  integer :: is, ie, js, je, nz, l, m
  character(len=8) :: mesg_l
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  do m=1,Reg%ntr_ice ; do l=1,Reg%Tr_ice(m)%nL
    write(mesg_l,'("i")') l
    call hchksum(Reg%Tr_ice(m)%t(:,:,:,l), mesg//trim(Reg%Tr_ice(m)%name)//" "//&
                                 trim(adjustl(mesg_l)), G)
  enddo ; enddo
  do m=1,Reg%ntr_snow ; do l=1,Reg%Tr_snow(m)%nL
    write(mesg_l,'("i")') l
    call hchksum(Reg%Tr_snow(m)%t(:,:,:,l), mesg//trim(Reg%Tr_snow(m)%name)//" "//&
                                 trim(adjustl(mesg_l)), G)
  enddo ; enddo
end subroutine SIS_tracer_chksum

subroutine SIS_tracer_registry_init(param_file, Reg)
  type(param_file_type),      intent(in) :: param_file
  type(SIS_tracer_registry_type), pointer    :: Reg
! Arguments: param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  Reg - A pointer that is set to point to the tracer registry.
  integer, save :: init_calls = 0
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_tracer_registry" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (.not.associated(Reg)) then ; allocate(Reg)
  else ; return ; endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")

  init_calls = init_calls + 1
  if (init_calls > 1) then
    write(mesg,'("SIS_tracer_registry_init called ",I3, &
      &" times with different registry pointers.")') init_calls
    call SIS_error(WARNING,"SIS_tracer"//mesg)
  endif

end subroutine SIS_tracer_registry_init

subroutine SIS_tracer_registry_end(Reg)
  type(SIS_tracer_registry_type), pointer :: Reg
  if (associated(Reg)) deallocate(Reg)
end subroutine SIS_tracer_registry_end

end module SIS_tracer_registry
