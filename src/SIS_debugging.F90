!> Routines that perform various error checking and debugging functions for SIS2
module SIS_debugging

! This file is a part of SIS2.  See LICENSE.md for the license file.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   This module contains subroutines that perform various error checking and   !
! debugging functions for SIS2.  This routine is similar to it counterpart in  !
! the MOM6 code, except for the use of the SIS_hor_grid_type and by keeping it !
! separate we retain the ability to set up MOM6 and SIS2 debugging separately. !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_checksums, only : hchksum, Bchksum
use MOM_checksums, only : mom_hchksum_pair => hchksum_pair
use MOM_checksums, only : mom_Bchksum_pair => Bchksum_pair
use MOM_checksums, only : mom_uvchksum => uvchksum

use MOM_checksums, only : is_NaN, chksum, MOM_checksums_init
use MOM_coms, only : PE_here, root_PE, num_PEs, sum_across_PEs
use MOM_coms, only : min_across_PEs, max_across_PEs, reproducing_sum
use MOM_domains, only : pass_vector, pass_var, pe_here
use MOM_domains, only : BGRID_NE, AGRID, To_All, Scalar_Pair
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : log_version, param_file_type, get_param
use SIS_hor_grid, only : SIS_hor_grid_type
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_hor_index, only : hor_index_type

implicit none ; private

public :: check_redundant_C, check_redundant_B, check_redundant_T
public :: SIS_debugging_init

! These interfaces come from MOM_checksums.
public :: hchksum, Bchksum, is_NaN, chksum
public :: hchksum_pair, Bchksum_pair, uvchksum

!> do checksums on a pair of fields
interface hchksum_pair
  module procedure hchksum_pair_3d, hchksum_pair_2d
end interface hchksum_pair

!> do checksums on a pair of fields at corner points
interface Bchksum_pair
  module procedure Bchksum_pair_3d, Bchksum_pair_2d
end interface Bchksum_pair

!> do checksums on the components of a vector
interface uvchksum
  module procedure uvchksum_3d, uvchksum_2d
  module procedure uvchksum_3d_dG, uvchksum_2d_dG
end interface uvchksum

!> Check the duplicate points of the components of a C-grid vector
interface check_redundant_C
  module procedure check_redundant_vC3d, check_redundant_vC2d
end interface check_redundant_C
!> Check the duplicate points of the components of a B-grid vector or scalar
interface check_redundant_B
  module procedure check_redundant_vB3d, check_redundant_vB2d
  module procedure check_redundant_sB3d, check_redundant_sB2d
end interface check_redundant_B
!> Check the duplicate points of the components of an A-grid vector or scalar
interface check_redundant_T
  module procedure check_redundant_sT3d, check_redundant_sT2d
  module procedure check_redundant_vT3d, check_redundant_vT2d
end interface check_redundant_T

integer :: max_redundant_prints = 100 !< The maximum number of error messages to print
integer :: redundant_prints(3) = 0  !< The maximum number of error messages to print
logical :: debug = .false.          !< If true, write out verbose debugging data.
logical :: debug_chksums = .true.   !< If true, checksums are performed on arrays in the various vec_chksum routines.
logical :: debug_redundant = .true. !< If true, debug redundant data points during calls to the
                                    !! various vec_chksum routines.

contains

! =====================================================================

!> SIS_debugging_init initializes the SIS_debugging module, and sets
!! the parameters that control which checks are active for SIS2.
subroutine SIS_debugging_init(param_file)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "SIS_debugging" ! This module's name.

  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_CHKSUMS", debug_chksums, &
                 "If true, checksums are performed on arrays in the "//&
                 "various vec_chksum routines.", default=debug, &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_REDUNDANT", debug_redundant, &
                 "If true, debug redundant data points during calls to "//&
                 "the various vec_chksum routines.", default=debug, &
                 debuggingParam=.true.)

  call MOM_checksums_init(param_file)

end subroutine SIS_debugging_init

!> Check duplicated points of 3-d C-grid vectors for consistency and
!! document offending points
subroutine check_redundant_vC3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                    intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),             intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%IsdB:,G%jsd:,:),   intent(in)    :: u_comp !< The u-component of the vector being checked
  real, dimension(G%isd:,G%JsdB:,:),   intent(in)    :: v_comp !< The v-component of the vector being checked
  integer,                   optional, intent(in)    :: is   !< The starting i-index to work on
  integer,                   optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,                   optional, intent(in)    :: js   !< The starting j-index to work on
  integer,                   optional, intent(in)    :: je   !< The ending j-index to work on
  integer,                   optional, intent(in)    :: direction !< The direction flag to pass to pass_vector

  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(u_comp,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_vC2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction)
  enddo
end subroutine  check_redundant_vC3d

!> Check duplicated points of 2-d C-grid vectors for consistency and
!! document offending points
subroutine check_redundant_vC2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),         intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector being checked
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector being checked
  integer,               optional, intent(in)    :: is   !< The starting i-index to work on
  integer,               optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,               optional, intent(in)    :: js   !< The starting j-index to work on
  integer,               optional, intent(in)    :: je   !< The ending j-index to work on
  integer,               optional, intent(in)    :: direction !< The direction flag to pass to pass_vector

  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: u_resym(G%IsdB:G%IedB,G%jsd:G%jed)
  real :: v_resym(G%isd:G%ied,G%JsdB:G%JedB)
  character(len=128) :: mesg2

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL," check_redundant"//&
    " called with a non-associated auxiliary domain the grid type.")
  call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, direction)

  do I=IsdB,IedB ; do j=jsd,jed ; u_resym(I,j) = u_comp(I,j) ; enddo ; enddo
  do i=isd,ied ; do J=JsdB,JedB ; v_resym(i,J) = v_comp(i,J) ; enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    u_resym(i,j) = u_nonsym(i,j) ; v_resym(i,j) = v_nonsym(i,j)
  enddo ; enddo
  call pass_vector(u_resym, v_resym, G%Domain, direction)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch+1,je_ch
    if (u_resym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(3) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_resym(i,j),u_comp(i,j)-u_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(3) = redundant_prints(3) + 1
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(3) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_resym(i,j),v_comp(i,j)-v_resym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(3) = redundant_prints(3) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_vC2d

!> Check duplicated points of 3-d B-grid scalars for consistency and
!! document offending points
subroutine check_redundant_sB3d(mesg, array, G, is, ie, js, je)
  character(len=*),                     intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),              intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%IsdB:,G%JsdB:,:),   intent(in)    :: array !< The array being checked
  integer,                    optional, intent(in)    :: is   !< The starting i-index to work on
  integer,                    optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,                    optional, intent(in)    :: js   !< The starting j-index to work on
  integer,                    optional, intent(in)    :: je   !< The ending j-index to work on

  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(array,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_sB2d(trim(mesg)//trim(mesg_k), array(:,:,k), &
                              G, is, ie, js, je)
  enddo
end subroutine  check_redundant_sB3d

!> Check duplicated points of 2-d B-grid scalars for consistency and
!! document offending points
subroutine check_redundant_sB2d(mesg, array, G, is, ie, js, je)
  character(len=*),                 intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),          intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: array !< The array being checked
  integer,                optional, intent(in)    :: is   !< The starting i-index to work on
  integer,                optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,                optional, intent(in)    :: js   !< The starting j-index to work on
  integer,                optional, intent(in)    :: je   !< The ending j-index to work on

  real :: a_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: a_resym(G%IsdB:G%IedB,G%JsdB:G%JedB)
  character(len=128) :: mesg2

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  do i=isd,ied ; do j=jsd,jed
    a_nonsym(i,j) = array(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL," check_redundant"//&
    " called with a non-associated auxiliary domain the grid type.")
  call pass_vector(a_nonsym, a_nonsym, G%Domain_aux, &
                   direction=To_All+Scalar_Pair, stagger=BGRID_NE)

  do I=IsdB,IedB ; do J=JsdB,JedB ; a_resym(I,J) = array(I,J) ; enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    a_resym(i,j) = a_nonsym(i,j)
  enddo ; enddo
  call pass_vector(a_resym, a_resym, G%Domain, direction=To_All+Scalar_Pair, &
                   stagger=BGRID_NE)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (a_resym(i,j) /= array(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           array(i,j), a_resym(i,j),array(i,j)-a_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_sB2d

!> Check duplicated points of 3-d B-grid vectors for consistency and
!! document offending points
subroutine check_redundant_vB3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                    intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),             intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%IsdB:,G%JsdB:,:),  intent(in)    :: u_comp !< The u-component of the vector being checked
  real, dimension(G%IsdB:,G%JsdB:,:),  intent(in)    :: v_comp !< The v-component of the vector being checked
  integer,                   optional, intent(in)    :: is   !< The starting i-index to work on
  integer,                   optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,                   optional, intent(in)    :: js   !< The starting j-index to work on
  integer,                   optional, intent(in)    :: je   !< The ending j-index to work on
  integer,                   optional, intent(in)    :: direction !< The direction flag to pass to pass_vector

  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(u_comp,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_vB2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction)
  enddo
end subroutine  check_redundant_vB3d

!> Check duplicated points of 2-d B-grid vectors for consistency and
!! document offending points
subroutine check_redundant_vB2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),         intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%IsdB:,G%JsdB:), intent(in)   :: u_comp !< The u-component of the vector being checked
  real, dimension(G%IsdB:,G%JsdB:), intent(in)   :: v_comp !< The v-component of the vector being checked
  integer,               optional, intent(in)    :: is   !< The starting i-index to work on
  integer,               optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,               optional, intent(in)    :: js   !< The starting j-index to work on
  integer,               optional, intent(in)    :: je   !< The ending j-index to work on
  integer,               optional, intent(in)    :: direction !< The direction flag to pass to pass_vector

  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: u_resym(G%IsdB:G%IedB,G%JsdB:G%JedB)
  real :: v_resym(G%IsdB:G%IedB,G%JsdB:G%JedB)
  character(len=128) :: mesg2

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL," check_redundant"//&
    " called with a non-associated auxiliary domain the grid type.")
  call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, direction, stagger=BGRID_NE)

  do I=IsdB,IedB ; do J=JsdB,JedB
    u_resym(I,J) = u_comp(I,J) ; v_resym(I,J) = v_comp(I,J)
  enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    u_resym(i,j) = u_nonsym(i,j) ; v_resym(i,j) = v_nonsym(i,j)
  enddo ; enddo
  call pass_vector(u_resym, v_resym, G%Domain, direction, stagger=BGRID_NE)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (u_resym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_resym(i,j),u_comp(i,j)-u_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo
  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_resym(i,j),v_comp(i,j)-v_resym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_vB2d

!> Check duplicated points of 3-d cell-centered scalars for consistency and
!! document offending points
subroutine check_redundant_sT3d(mesg, array, G, is, ie, js, je)
  character(len=*),                     intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),              intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%isd:,G%jsd:,:),     intent(in)    :: array !< The array being checked
  integer,                    optional, intent(in)    :: is   !< The starting i-index to work on
  integer,                    optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,                    optional, intent(in)    :: js   !< The starting j-index to work on
  integer,                    optional, intent(in)    :: je   !< The ending j-index to work on

  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(array,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_sT2d(trim(mesg)//trim(mesg_k), array(:,:,k), &
                              G, is, ie, js, je)
  enddo
end subroutine  check_redundant_sT3d

!> Check duplicated points of 2-d cell-centered scalars for consistency and
!! document offending points
subroutine check_redundant_sT2d(mesg, array, G, is, ie, js, je)
  character(len=*),                 intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),          intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%isd:,G%jsd:),   intent(in)    :: array !< The array being checked
  integer,                optional, intent(in)    :: is   !< The starting i-index to work on
  integer,                optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,                optional, intent(in)    :: js   !< The starting j-index to work on
  integer,                optional, intent(in)    :: je   !< The ending j-index to work on

  real :: a_nonsym(G%isd:G%ied,G%jsd:G%jed)
  character(len=128) :: mesg2

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  is_ch = G%isc ; ie_ch = G%iec ; js_ch = G%jsc ; je_ch = G%jec
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  ! This only works on points outside of the standard computational domain.
  if ((is_ch == G%isc) .and. (ie_ch == G%iec) .and. &
      (js_ch == G%jsc) .and. (je_ch == G%jec)) return

  do i=isd,ied ; do j=jsd,jed
    a_nonsym(i,j) = array(i,j)
  enddo ; enddo

  call pass_var(a_nonsym, G%Domain)

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (a_nonsym(i,j) /= array(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           array(i,j), a_nonsym(i,j),array(i,j)-a_nonsym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
     redundant_prints(1) = redundant_prints(1) + 1
   endif
  enddo ; enddo

end subroutine  check_redundant_sT2d

!> Check duplicated points of 3-d A-grid vectors for consistency and
!! document offending points
subroutine check_redundant_vT3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction)
  character(len=*),                    intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),             intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%isd:,G%jsd:,:),    intent(in)    :: u_comp !< The u-component of the vector being checked
  real, dimension(G%isd:,G%jsd:,:),    intent(in)    :: v_comp !< The v-component of the vector being checked
  integer,                   optional, intent(in)    :: is   !< The starting i-index to work on
  integer,                   optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,                   optional, intent(in)    :: js   !< The starting j-index to work on
  integer,                   optional, intent(in)    :: je   !< The ending j-index to work on
  integer,                   optional, intent(in)    :: direction !< The direction flag to pass to pass_vector

  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(u_comp,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_vT2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction)
  enddo
end subroutine  check_redundant_vT3d

!> Check duplicated points of 2-d A-grid vectors for consistency and
!! document offending points
subroutine check_redundant_vT2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction)
  character(len=*),                intent(in)    :: mesg !< An identifying message
  type(SIS_hor_grid_type),         intent(inout) :: G    !< The horizontal grid type
  real, dimension(G%isd:,G%jsd:),  intent(in)    :: u_comp !< The u-component of the vector being checked
  real, dimension(G%isd:,G%jsd:),  intent(in)    :: v_comp !< The v-component of the vector being checked
  integer,               optional, intent(in)    :: is   !< The starting i-index to work on
  integer,               optional, intent(in)    :: ie   !< The ending i-index to work on
  integer,               optional, intent(in)    :: js   !< The starting j-index to work on
  integer,               optional, intent(in)    :: je   !< The ending j-index to work on
  integer,               optional, intent(in)    :: direction !< The direction flag to pass to pass_vector

  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed)
  character(len=128) :: mesg2

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  is_ch = G%isc ; ie_ch = G%iec ; js_ch = G%jsc ; je_ch = G%jec
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  ! This only works on points outside of the standard computational domain.
  if ((is_ch == G%isc) .and. (ie_ch == G%iec) .and. &
      (js_ch == G%jsc) .and. (je_ch == G%jec)) return

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  call pass_vector(u_nonsym, v_nonsym, G%Domain, direction, stagger=AGRID)

  do i=is_ch,ie_ch ; do j=js_ch+1,je_ch
    if (u_nonsym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_nonsym(i,j),u_comp(i,j)-u_nonsym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(1) = redundant_prints(1) + 1
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_nonsym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_nonsym(i,j),v_comp(i,j)-v_nonsym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(1) = redundant_prints(1) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_vT2d

! =====================================================================

!> This subroutine does a checksum and redundant point check on a 3d C-grid vector.
subroutine uvchksum_3d(mesg, u_comp, v_comp, G, halos, scalars, scale)
  character(len=*),                  intent(in)    :: mesg   !< An identifying message
  type(SIS_hor_grid_type),           intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                 optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                 optional, intent(in)    :: scalars !< If true this is a pair of
                                                              !! scalars that are being checked.
  real,                    optional, intent(in)    :: scale  !< A scaling factor for these arrays.

  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call mom_uvchksum(mesg, u_comp, v_comp, G%HI, halos, scale=scale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_C(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_C(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine uvchksum_3d

!> This subroutine does a checksum and redundant point check on a 2d C-grid vector.
subroutine uvchksum_2d(mesg, u_comp, v_comp, G, halos, scalars, scale)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(SIS_hor_grid_type),         intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector
  integer,               optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,               optional, intent(in)    :: scalars !< If true this is a pair of
                                                            !! scalars that are being checked.
  real,                  optional, intent(in)    :: scale  !< A scaling factor for these arrays.

  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call mom_uvchksum(mesg, u_comp, v_comp, G%HI, halos, scale=scale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_C(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_C(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine uvchksum_2d

!> This subroutine does a checksum and redundant point check on a 3d C-grid vector.
subroutine uvchksum_3d_dG(mesg, u_comp, v_comp, G, halos, scale)
  character(len=*),                  intent(in)    :: mesg   !< An identifying message
  type(dyn_horgrid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                 optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  real,                    optional, intent(in)    :: scale  !< A scaling factor for these arrays.

  if (debug_chksums) then
    call mom_uvchksum(mesg, u_comp, v_comp, G%HI, halos, scale=scale)
  endif

end subroutine uvchksum_3d_dG

!> This subroutine does a checksum and redundant point check on a 2d C-grid vector.
subroutine uvchksum_2d_dG(mesg, u_comp, v_comp, G, halos, scale)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(dyn_horgrid_type),          intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector
  integer,               optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  real,                  optional, intent(in)    :: scale  !< A scaling factor for these arrays.

  if (debug_chksums) then
    call mom_uvchksum(mesg, u_comp, v_comp, G%HI, halos, scale=scale)
  endif

end subroutine uvchksum_2d_dG

!> This subroutine does a checksum and redundant point check on a 3d B-grid vector.
subroutine Bchksum_pair_3d(mesg, u_comp, v_comp, G, halos, scalars, scale)
  character(len=*),                   intent(in)    :: mesg   !< An identifying message
  type(SIS_hor_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                  optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                  optional, intent(in)    :: scalars !< If true this is a pair of
                                                               !! scalars that are being checked.
  real,                     optional, intent(in)    :: scale  !< A scaling factor for these arrays.

  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call mom_Bchksum_pair(mesg, u_comp, v_comp, G%HI, halos, scale=scale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_B(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_B(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine Bchksum_pair_3d

!> This subroutine does a checksum and redundant point check on a 2d B-grid vector.
subroutine Bchksum_pair_2d(mesg, u_comp, v_comp, G, halos, scalars, symmetric, scale)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(SIS_hor_grid_type),          intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                optional, intent(in)    :: scalars !< If true this is a pair of
                                                             !! scalars that are being checked.
  logical,                optional, intent(in)    :: symmetric !< If true, do the checksums on the
                                                               !! full symmetric computational domain.
  real,                   optional, intent(in)    :: scale   !< A scaling factor for these arrays.

  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call mom_Bchksum_pair(mesg, u_comp, v_comp, G%HI, symmetric=symmetric, haloshift=halos, scale=scale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_B(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_B(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine Bchksum_pair_2d

!> This subroutine does a checksum and redundant point check on a 3d C-grid vector.
subroutine hchksum_pair_3d(mesg, u_comp, v_comp, G, halos, scalars, scale)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(SIS_hor_grid_type),          intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%jsd:,:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                optional, intent(in)    :: scalars !< If true this is a pair of
                                                             !! scalars that are being checked.
  real,                   optional, intent(in)    :: scale   !< A scaling factor for these arrays.

  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call mom_hchksum_pair(mesg, u_comp, v_comp, G%HI, halos, scale=scale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_T(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_T(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine hchksum_pair_3d

!> This subroutine does a checksum and redundant point check on a 2d C-grid vector.
subroutine hchksum_pair_2d(mesg, u_comp, v_comp, G, halos, scalars, scale)
  character(len=*),               intent(in)    :: mesg   !< An identifying message
  type(SIS_hor_grid_type),        intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%jsd:), intent(in)    :: v_comp !< The v-component of the vector
  integer,              optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,              optional, intent(in)    :: scalars !< If true this is a pair of
                                                           !! scalars that are being checked.
  real,                  optional, intent(in)   :: scale   !< A scaling factor for these arrays.

  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call mom_hchksum_pair(mesg, u_comp, v_comp, G%HI, halos, scale=scale)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_T(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_T(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine hchksum_pair_2d

end module SIS_debugging
