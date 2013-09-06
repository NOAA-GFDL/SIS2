module SIS_error_checking

!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

use MOM_coms, only : PE_here, root_PE, num_PEs, sum_across_PEs
use MOM_coms, only : min_across_PEs, max_across_PEs, reproducing_sum
use MOM_domains, only : pass_vector, pass_var, pe_here
use MOM_domains, only : BGRID_NE, AGRID, To_All, Scalar_Pair
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : log_version, param_file_type
use ice_grid_mod, only : sea_ice_grid_type

implicit none ; private

public :: hchksum, Bchksum, uchksum, vchksum, chksum, is_NaN
public :: check_redundant_C, check_redundant_B, check_redundant_T
public :: MOM_checksums_init

interface hchksum
  module procedure chksum_h_2d
  module procedure chksum_h_3d
end interface

interface Bchksum
  module procedure chksum_B_2d
  module procedure chksum_B_3d
end interface

interface uchksum
  module procedure chksum_u_2d
  module procedure chksum_u_3d
end interface

interface vchksum
  module procedure chksum_v_2d
  module procedure chksum_v_3d
end interface

interface chksum
  module procedure chksum1d
  module procedure chksum2d
  module procedure chksum3d
end interface

interface chk_sum_msg
  module procedure chk_sum_msg1
  module procedure chk_sum_msg2
  module procedure chk_sum_msg3
  module procedure chk_sum_msg5
end interface

interface is_NaN
  module procedure is_NaN_0d
  module procedure is_NaN_1d
  module procedure is_NaN_2d
  module procedure is_NaN_3d
end interface

interface check_redundant_C
  module procedure check_redundant_vC3d, check_redundant_vC2d
end interface check_redundant_C
interface check_redundant_B
  module procedure check_redundant_vB3d, check_redundant_vB2d
  module procedure check_redundant_sB3d, check_redundant_sB2d
end interface check_redundant_B
interface check_redundant_T
  module procedure check_redundant_sT3d, check_redundant_sT2d
  module procedure check_redundant_vT3d, check_redundant_vT2d
end interface check_redundant_T

integer, parameter :: default_shift=0
logical :: calculateStatistics=.false. ! If true, report min, max and mean
                               ! instead of the bitcount checksum
logical :: checkForNaNs=.true. ! If true, checks array for NaNs and cause
                               ! FATAL error is any are found

contains

! =====================================================================

subroutine chksum_h_2d(array, mesg, G, haloshift)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(G%isd:,G%jsd:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%isc:G%iec,G%jsc:G%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_h_2d: haloshift =',hshift
    write(0,*) 'chksum_h_2d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_h_2d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_h_2d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%jsd:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%jsd:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc)
    aMax = array(G%isc,G%jsc)
    n = 0
    do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j)
        aMin = min(aMin, array(i,j))
        aMax = max(aMax, array(i,j))
        n = n + 1
    enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_2d

! =====================================================================

subroutine chksum_B_2d(array, mesg, G, haloshift, symmetric)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(G%IsdB:,G%JsdB:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift
  logical, intent(in), optional :: symmetric

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift
  logical :: sym

  if (checkForNaNs) then
    if (is_NaN(array(G%IscB:G%IecB,G%JscB:G%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_B_2d: haloshift =',hshift
    write(0,*) 'chksum_B_2d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_B_2d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_B_2d '//trim(mesg))
  endif

  sym = .false. ; if (present(symmetric)) sym = symmetric

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
    if (sym) then
      bcSW=subchk(array, G, -1, -1)
      if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,mesg)
    else
      if (is_root_pe()) call chk_sum_msg("B-point:",bc0,mesg)
    endif
    return
  endif

  if (sym) then
    bcSW=subchk(array, G, -hshift-1, -hshift-1)
    bcSE=subchk(array, G, hshift, -hshift-1)
    bcNW=subchk(array, G, -hshift-1, hshift)
  else
    bcSW=subchk(array, G, -hshift, -hshift)
    bcSE=subchk(array, G, hshift, -hshift)
    bcNW=subchk(array, G, -hshift, hshift)
  endif
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%JsdB:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%JsdB:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc)
    aMax = array(G%isc,G%jsc)
    n = 0
    do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j)
        aMin = min(aMin, array(i,j))
        aMax = max(aMax, array(i,j))
        n = n + 1
    enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("B-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_B_2d

! =====================================================================

subroutine chksum_u_2d(array, mesg, G, haloshift)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(G%IsdB:,G%jsd:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%IscB:G%IecB,G%jsc:G%jec))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_u_2d: haloshift =',hshift
    write(0,*) 'chksum_u_2d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_u_2d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_u_2d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%jsd:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%jsd:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc)
    aMax = array(G%isc,G%jsc)
    n = 0
    do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j)
        aMin = min(aMin, array(i,j))
        aMax = max(aMax, array(i,j))
        n = n + 1
    enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_2d

! =====================================================================

subroutine chksum_v_2d(array, mesg, G, haloshift)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(G%isd:,G%JsdB:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%isc:G%iec,G%JscB:G%JecB))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_v_2d: haloshift =',hshift
    write(0,*) 'chksum_v_2d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_v_2d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_v_2d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%JsdB:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, bc
    subchk = 0
    do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j)))
        subchk = subchk + bc
    enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%JsdB:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc)
    aMax = array(G%isc,G%jsc)
    n = 0
    do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j)
        aMin = min(aMin, array(i,j))
        aMax = max(aMax, array(i,j))
        n = n + 1
    enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_v_2d

! =====================================================================

subroutine chksum_h_3d(array, mesg, G, haloshift)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(G%isd:,G%jsd:,:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%isc:G%iec,G%jsc:G%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_h_3d: haloshift =',hshift
    write(0,*) 'chksum_h_3d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_h_3d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_h_3d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("h-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("h-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%jsd:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j,k)))
        subchk = subchk + bc
    enddo; enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%jsd:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc,1)
    aMax = array(G%isc,G%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j,k)
        aMin = min(aMin, array(i,j,k))
        aMax = max(aMax, array(i,j,k))
        n = n + 1
    enddo; enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("h-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_h_3d

! =====================================================================

subroutine chksum_B_3d(array, mesg, G, haloshift)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%IscB:G%IecB,G%JscB:G%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_B_3d: haloshift =',hshift
    write(0,*) 'chksum_B_3d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_B_3d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_B_3d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("B-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("B-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%JsdB:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j,k)))
        subchk = subchk + bc
    enddo; enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%JsdB:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc,1)
    aMax = array(G%isc,G%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j,k)
        aMin = min(aMin, array(i,j,k))
        aMax = max(aMax, array(i,j,k))
        n = n + 1
    enddo; enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("B-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_B_3d

! =====================================================================

subroutine chksum_u_3d(array, mesg, G, haloshift)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(G%IsdB:,G%jsd:,:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%IscB:G%IecB,G%jsc:G%jec,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_u_3d: haloshift =',hshift
    write(0,*) 'chksum_u_3d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_u_3d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_u_3d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("u-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("u-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%jsd:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j,k)))
        subchk = subchk + bc
    enddo; enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%IsdB:,G%jsd:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc,1)
    aMax = array(G%isc,G%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j,k)
        aMin = min(aMin, array(i,j,k))
        aMax = max(aMax, array(i,j,k))
        n = n + 1
    enddo; enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("u-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_u_3d

! =====================================================================

subroutine chksum_v_3d(array, mesg, G, haloshift)
  type(sea_ice_grid_type), intent(in) :: G
  real, dimension(G%isd:,G%JsdB:,:), intent(in) :: array
  character(len=*), intent(in) :: mesg
  integer, intent(in), optional :: haloshift

  integer :: bc0,bcSW,bcSE,bcNW,bcNE,hshift

  if (checkForNaNs) then
    if (is_NaN(array(G%isc:G%iec,G%JscB:G%JecB,:))) &
      call chksum_error(FATAL, 'NaN detected: '//trim(mesg))
!   if (is_NaN(array)) &
!     call chksum_error(FATAL, 'NaN detected in halo: '//trim(mesg))
  endif

  if (calculateStatistics) then
    call subStats(G, array, mesg); return
  endif

  hshift=default_shift
  if (present(haloshift)) hshift=haloshift
  if (hshift<0) hshift=G%ied-G%iec

  if ( G%isc-hshift<G%isd .or. &
       G%iec+hshift>G%ied .or. &
       G%jsc-hshift<G%jsd .or. &
       G%jec+hshift>G%jed ) then
    write(0,*) 'chksum_v_3d: haloshift =',hshift
    write(0,*) 'chksum_v_3d: isd,isc,iec,ied=',G%isd,G%isc,G%iec,G%ied
    write(0,*) 'chksum_v_3d: jsd,jsc,jec,jed=',G%jsd,G%jsc,G%jec,G%jed
    call chksum_error(FATAL,'Error in chksum_v_3d '//trim(mesg))
  endif

  bc0=subchk(array, G, 0, 0)

  if (hshift==0) then
      if (is_root_pe()) call chk_sum_msg("v-point:",bc0,mesg)
      return
  endif

  bcSW=subchk(array, G, -hshift, -hshift)
  bcSE=subchk(array, G, hshift, -hshift)
  bcNW=subchk(array, G, -hshift, hshift)
  bcNE=subchk(array, G, hshift, hshift)

  if (is_root_pe()) call chk_sum_msg("v-point:",bc0,bcSW,bcSE,bcNW,bcNE,mesg)

  contains

  integer function subchk(array, G, di, dj)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%JsdB:,:), intent(in) :: array
    integer, intent(in) :: di, dj
    integer :: bitcount, i, j, k, bc
    subchk = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc+dj,G%jec+dj; do i=G%isc+di,G%iec+di
        bc = bitcount(abs(array(i,j,k)))
        subchk = subchk + bc
    enddo; enddo; enddo
    call sum_across_PEs(subchk)
    subchk=mod(subchk,1000000000)
  end function subchk

  subroutine subStats(G, array, mesg)
    type(sea_ice_grid_type), intent(in) :: G
    real, dimension(G%isd:,G%JsdB:,:), intent(in) :: array
    character(len=*), intent(in) :: mesg
    integer :: i, j, k, n
    real :: aMean, aMin, aMax
    aMean = 0.
    aMin = array(G%isc,G%jsc,1)
    aMax = array(G%isc,G%jsc,1)
    n = 0
    do k=LBOUND(array,3),UBOUND(array,3); do j=G%jsc,G%jec; do i=G%isc,G%iec
        aMean = aMean + array(i,j,k)
        aMin = min(aMin, array(i,j,k))
        aMax = max(aMax, array(i,j,k))
        n = n + 1
    enddo; enddo; enddo
    call sum_across_PEs(aMean)
    call sum_across_PEs(n)
    call min_across_PEs(aMin)
    call max_across_PEs(aMax)
    aMean = aMean / real(n)
    if (is_root_pe()) call chk_sum_msg("v-point:",aMean,aMin,aMax,mesg)
  end subroutine subStats

end subroutine chksum_v_3d


! =====================================================================
subroutine chksum1d(array, mesg, start_x, end_x)

  real, dimension(:) :: array
  character(len=*) :: mesg
  integer, optional :: start_x, end_x

  integer :: i, is, ie, bc, sum1
  integer :: bitcount
  real :: sum

  is = LBOUND(array,1) ; ie = UBOUND(array,1)
  if (present(start_x)) is = start_x
  if (present(end_x)) ie = end_x

  sum1 = 0 ; sum = 0
  do i=is,ie
    bc = bitcount(ABS(array(i)))
    sum1 = sum1 + bc
    sum = sum + array(i)
  enddo

  call sum_across_PEs(sum)
  call sum_across_PEs(sum1)

  if (is_root_pe()) &
    write(0,'(A40,1X,ES22.13,1X,I12)') mesg, sum, sum1

end subroutine chksum1d

! =====================================================================
subroutine chksum2d(array, mesg)

  real, dimension(:,:) :: array
  character(len=*) :: mesg

  integer :: i, j, bc, sum1
  integer :: bitcount
  real :: sum

  sum1 = 0
  do j=1,size(array,2); do i=1,size(array,1)
    bc = bitcount(ABS(array(i,j)))
    sum1 = sum1 + bc
  enddo ; enddo

  sum = reproducing_sum(array(:,:))
  call sum_across_PEs(sum1)

  if (is_root_pe()) &
    write(0,'(A40,1X,ES22.13,1X,I12)') mesg, sum, sum1

end subroutine chksum2d

subroutine chksum3d(array, mesg)

  real, dimension(:,:,:) :: array
  character(len=*) :: mesg

  integer :: i, j, k, bc, sum1
  integer :: bitcount
  real :: sum

  sum1 = 0
  do k=1,size(array,3); do j=1,size(array,2); do i=1,size(array,1)
    bc = bitcount(ABS(array(i,j,k)))
    sum1 = sum1 + bc
  enddo ; enddo ; enddo

  sum = reproducing_sum(array(:,:,:))
  call sum_across_PEs(sum1)

  if (is_root_pe()) &
    write(0,'(A40,1X,ES22.13,1X,I12)') mesg, sum, sum1

end subroutine chksum3d

! =====================================================================

function is_NaN_0d(x)
  real, intent(in) :: x
  logical :: is_NaN_0d
! This subroutine returns .true. if x is a NaN, and .false. otherwise.

 !is_NaN_0d = (((x < 0.0) .and. (x >= 0.0)) .or. &
 !          (.not.(x < 0.0) .and. .not.(x >= 0.0)))
  if (((x < 0.0) .and. (x >= 0.0)) .or. &
            (.not.(x < 0.0) .and. .not.(x >= 0.0))) then
    is_NaN_0d = .true.
  else
    is_NaN_0d = .false.
  endif

end function is_NaN_0d

! =====================================================================

function is_NaN_1d(x,skip_mpp)
  real, dimension(:), intent(in) :: x
  logical :: is_NaN_1d
  logical, optional :: skip_mpp
! This subroutine returns .true. if any x is a NaN, and .false. otherwise.
  integer :: i, n
  logical :: call_mpp

  n = 0
  do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i))) n = n + 1
  enddo
  call_mpp = .true.
  if (present(skip_mpp)) then
    if (skip_mpp) call_mpp = .false.
  endif
  if (call_mpp) call sum_across_PEs(n)
  is_NaN_1d = .false.
  if (n>0) is_NaN_1d = .true.

end function is_NaN_1d

! =====================================================================

function is_NaN_2d(x)
  real, dimension(:,:), intent(in) :: x
  logical :: is_NaN_2d
! This subroutine returns .true. if any x is a NaN, and .false. otherwise.
  integer :: i, j, n

  n = 0
  do j = LBOUND(x,2), UBOUND(x,2) ; do i = LBOUND(x,1), UBOUND(x,1)
    if (is_NaN_0d(x(i,j))) n = n + 1
  enddo ; enddo
  call sum_across_PEs(n)
  is_NaN_2d = .false.
  if (n>0) is_NaN_2d = .true.

end function is_NaN_2d

! =====================================================================

function is_NaN_3d(x)
  real, dimension(:,:,:), intent(in) :: x
  logical :: is_NaN_3d
! This subroutine returns .true. if any x is a NaN, and .false. otherwise.
  integer :: i, j, k, n

  n = 0
  do k = LBOUND(x,3), UBOUND(x,3)
    do j = LBOUND(x,2), UBOUND(x,2) ; do i = LBOUND(x,1), UBOUND(x,1)
      if (is_NaN_0d(x(i,j,k))) n = n + 1
    enddo ; enddo
  enddo
  call sum_across_PEs(n)
  is_NaN_3d = .false.
  if (n>0) is_NaN_3d = .true.

end function is_NaN_3d

! =====================================================================

subroutine chk_sum_msg1(fmsg,bc0,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0
  if (is_root_pe()) write(0,'(A,1(A,I9,X),A)') fmsg," c=",bc0,mesg
end subroutine chk_sum_msg1

! =====================================================================

subroutine chk_sum_msg5(fmsg,bc0,bcSW,bcSE,bcNW,bcNE,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0,bcSW,bcSE,bcNW,bcNE
  if (is_root_pe()) write(0,'(A,5(A,I9,1X),A)') &
     fmsg," c=",bc0,"sw=",bcSW,"se=",bcSE,"nw=",bcNW,"ne=",bcNE,mesg
end subroutine chk_sum_msg5

! =====================================================================

subroutine chk_sum_msg2(fmsg,bc0,bcSW,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  integer,          intent(in) :: bc0,bcSW
  if (is_root_pe()) write(0,'(A,2(A,I9,1X),A)') &
     fmsg," c=",bc0,"s/w=",bcSW,mesg
end subroutine chk_sum_msg2

! =====================================================================

subroutine chk_sum_msg3(fmsg,aMean,aMin,aMax,mesg)
  character(len=*), intent(in) :: fmsg, mesg
  real,             intent(in) :: aMean,aMin,aMax
  if (is_root_pe()) write(0,'(A,3(A,ES12.4,1X),A)') &
     fmsg," mean=",aMean,"min=",aMin,"max=",aMax,mesg
end subroutine chk_sum_msg3

! =====================================================================

subroutine MOM_checksums_init(param_file)
  type(param_file_type),   intent(in)    :: param_file
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_checksums" ! This module's name.

  call log_version(param_file, mod, version)

end subroutine MOM_checksums_init

! =====================================================================

subroutine chksum_error(signal, message)
  ! Wrapper for MOM_error to help place specific break points in
  ! debuggers
  integer, intent(in) :: signal
  character(len=*), intent(in) :: message
  call MOM_error(signal, message)
end subroutine chksum_error

! =====================================================================


subroutine check_redundant_vC3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                    intent(in)    :: mesg
  type(sea_ice_grid_type),             intent(inout) :: G
  real, dimension(G%IsdB:,G%jsd:,:),   intent(in)    :: u_comp
  real, dimension(G%isd:,G%JsdB:,:),   intent(in)    :: v_comp
  integer,                   optional, intent(in)    :: is, ie, js, je
  integer,                   optional, intent(in)    :: direction
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  direction - the direction flag to be passed to pass_vector.

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

subroutine check_redundant_vC2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                intent(in)    :: mesg
  type(sea_ice_grid_type),         intent(inout) :: G
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp
  integer,               optional, intent(in)    :: is, ie, js, je
  integer,               optional, intent(in)    :: direction
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  direction - the direction flag to be passed to pass_vector.

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
    if (u_resym(i,j) /= u_comp(i,j)) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_resym(i,j),u_comp(i,j)-u_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j)) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_resym(i,j),v_comp(i,j)-v_resym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo

end subroutine  check_redundant_vC2d

subroutine check_redundant_sB3d(mesg, array, G, is, ie, js, je)
  character(len=*),                     intent(in)    :: mesg
  type(sea_ice_grid_type),              intent(inout) :: G
  real, dimension(G%IsdB:,G%JsdB:,:),   intent(in) :: array
  integer,                    optional, intent(in)    :: is, ie, js, je
! Arguments: array - The array being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.

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


subroutine check_redundant_sB2d(mesg, array, G, is, ie, js, je)
  character(len=*),                intent(in)    :: mesg
  type(sea_ice_grid_type),         intent(inout) :: G
  real, dimension(G%IsdB:,G%JsdB:), intent(in)   :: array
  integer,               optional, intent(in)    :: is, ie, js, je
! Arguments: array - The array being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.

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
    if (a_resym(i,j) /= array(i,j)) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           array(i,j), a_resym(i,j),array(i,j)-a_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo

end subroutine  check_redundant_sB2d


subroutine check_redundant_vB3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction)
  character(len=*),                    intent(in)    :: mesg
  type(sea_ice_grid_type),             intent(inout) :: G
  real, dimension(G%IsdB:,G%JsdB:,:),  intent(in) :: u_comp
  real, dimension(G%IsdB:,G%JsdB:,:),  intent(in) :: v_comp
  integer,                   optional, intent(in)    :: is, ie, js, je
  integer,                   optional, intent(in)    :: direction
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  direction - the direction flag to be passed to pass_vector.

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

subroutine check_redundant_vB2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction)
  character(len=*),                intent(in)    :: mesg
  type(sea_ice_grid_type),         intent(inout) :: G
  real, dimension(G%IsdB:,G%JsdB:), intent(in)   :: u_comp
  real, dimension(G%IsdB:,G%JsdB:), intent(in)   :: v_comp
  integer,               optional, intent(in)    :: is, ie, js, je
  integer,               optional, intent(in)    :: direction
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  direction - the direction flag to be passed to pass_vector.

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
  
  do i=is_ch,ie_ch ; do j=js_ch+1,je_ch
    if (u_resym(i,j) /= u_comp(i,j)) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_resym(i,j),u_comp(i,j)-u_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j)) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_resym(i,j),v_comp(i,j)-v_resym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo

end subroutine  check_redundant_vB2d

subroutine check_redundant_sT3d(mesg, array, G, is, ie, js, je)
  character(len=*),                     intent(in)    :: mesg
  type(sea_ice_grid_type),              intent(inout) :: G
  real, dimension(G%IsdB:,G%JsdB:,:),   intent(in) :: array
  integer,                    optional, intent(in)    :: is, ie, js, je
! Arguments: array - The array being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.

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


subroutine check_redundant_sT2d(mesg, array, G, is, ie, js, je)
  character(len=*),                intent(in)    :: mesg
  type(sea_ice_grid_type),         intent(inout) :: G
  real, dimension(G%IsdB:,G%JsdB:), intent(in)   :: array
  integer,               optional, intent(in)    :: is, ie, js, je
! Arguments: array - The array being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.

  real :: a_nonsym(G%isd:G%ied,G%jsd:G%jed)
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
    a_nonsym(i,j) = array(i,j)
  enddo ; enddo

  call pass_var(a_nonsym, G%Domain)
  
  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (a_nonsym(i,j) /= array(i,j)) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           array(i,j), a_nonsym(i,j),array(i,j)-a_nonsym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo

end subroutine  check_redundant_sT2d


subroutine check_redundant_vT3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction)
  character(len=*),                    intent(in)    :: mesg
  type(sea_ice_grid_type),             intent(inout) :: G
  real, dimension(G%isd:,G%jsd:,:),    intent(in) :: u_comp
  real, dimension(G%isd:,G%jsd:,:),    intent(in) :: v_comp
  integer,                   optional, intent(in)    :: is, ie, js, je
  integer,                   optional, intent(in)    :: direction
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  direction - the direction flag to be passed to pass_vector.

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

subroutine check_redundant_vT2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction)
  character(len=*),                intent(in)    :: mesg
  type(sea_ice_grid_type),         intent(inout) :: G
  real, dimension(G%isd:,G%jsd:), intent(in)   :: u_comp
  real, dimension(G%isd:,G%jsd:), intent(in)   :: v_comp
  integer,               optional, intent(in)    :: is, ie, js, je
  integer,               optional, intent(in)    :: direction
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  direction - the direction flag to be passed to pass_vector.

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
    if (u_nonsym(i,j) /= u_comp(i,j)) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_nonsym(i,j),u_comp(i,j)-u_nonsym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_nonsym(i,j) /= v_comp(i,j)) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_nonsym(i,j),v_comp(i,j)-v_nonsym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo

end subroutine  check_redundant_vT2d

end module SIS_error_checking
