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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! This module contains convenient utilities for use by the SIS2 sea ice model. !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module SIS_utils

use MOM_coms,           only : g_sum=>reproducing_sum
use MOM_domains,        only : SCALAR_PAIR, CGRID_NE, BGRID_NE, To_All
use MOM_error_handler,  only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler,  only : is_root_pe
use SIS_diag_mediator,  only : post_SIS_data, SIS_diag_ctrl
use SIS_debugging,      only : hchksum, Bchksum, uchksum, vchksum
use SIS_debugging,      only : check_redundant_B, vec_chksum_A, vec_chksum_B, vec_chksum_C
use SIS_hor_grid,       only : SIS_hor_grid_type

implicit none ; private

public :: get_avg, post_avg, ice_line, is_NaN, g_sum, ice_grid_chksum

interface post_avg
  module procedure post_avg_3d, post_avg_4d
end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_avg - take area weighted average over some or all thickness categories.  !
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
!$OMP parallel do default(none) shared(ni,nj,nk,avg,cn,x,wts)
    do j=1,nj
      do k=1,nk ; do i=1,ni
        avg(i,j) = avg(i,j) + cn(i,j,k)*x(i,j,k)
        wts(i,j) = wts(i,j) + cn(i,j,k)
      enddo ; enddo
      do i=1,ni
        if (wts(i,j) > 0.) then
          avg(i,j) = avg(i,j) / wts(i,j)
        else
          avg(i,j) = 0.0
        endif
      enddo
    enddo
  else
    avg(:,:) = 0.0
!$OMP parallel do default(none) shared(ni,nj,nk,avg,cn,x)
    do j=1,nj
      do k=1,nk ; do i=1,ni
        avg(i,j) = avg(i,j) + cn(i,j,k)*x(i,j,k)
      enddo ; enddo
    enddo
  endif

end subroutine get_avg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_line - Write out a line with the northern and southern hemisphere ice    !
!            ice extents and global mean sea surface temperature.              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_line(year, day, second, cn_ocn, sst, G)
  integer,                         intent(in) :: year, day, second
  type(SIS_hor_grid_type),         intent(in) :: G
  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: cn_ocn
  real, dimension(G%isc:G%iec,G%jsc:G%jec), intent(in) :: sst

  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: x
  real :: gx(3)
  integer :: n, i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (.not.(second==0 .and. mod(day,5)==0) ) return

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
  if ( is_root_pe() .and. second==0 .and. mod(day,5)==0 ) &
    print '(a,2I4,3F10.5)','ICE y/d (SH_ext NH_ext SST):', year, day, gx
end subroutine ice_line

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! post_avg - take area weighted average over some or all thickness categories  !
!            and offer it for diagnostic output.                               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine post_avg_3d(id, val, part, diag, G, mask, scale, offset, wtd)
  integer,                 intent(in) :: id
  real, dimension(:,:,:),  intent(in) :: val, part
  type(SIS_diag_ctrl),     intent(in) :: diag
  type(SIS_hor_grid_type), optional, intent(in) :: G
  logical, dimension(:,:), optional, intent(in) :: mask
  real,                    optional, intent(in) :: scale, offset
  logical,                 optional, intent(in) :: wtd
  ! This subroutine determines the average of a quantity across thickness
  ! categories and does a send data on it.

  real :: avg(size(val,1),size(val,2)), wts(size(val,1),size(val,2))
  real :: scl, off
  logical :: do_wt
  integer :: i, j, k, ni, nj, nk, is, ie, js, je

  ni = size(val,1) ; nj = size(val,2) ; nk = size(val,3)
  if (size(part,1) /= ni) call SIS_error(FATAL, &
    "Mismatched i-sizes in post_avg.")
  if (size(part,2) /= nj) call SIS_error(FATAL, &
    "Mismatched j-sizes in post_avg.")
  if (size(part,3) /= nk) call SIS_error(FATAL, &
    "Mismatched k-sizes in post_avg.")

  if (present(G)) then  ! Account for the fact that arrays here start at 1.
    if ((ni == G%isc-G%iec + 1) .or. (ni == G%isc-G%iec + 2)) then
      is = 1; ie = ni  ! These arrays have no halos.
    elseif (ni == G%ied-(G%isd-1)) then  ! Arrays have halos.
      is = G%isc - (G%isd-1) ; ie = G%iec - (G%isd-1)
    elseif (ni == G%ied-(G%isd-1)+1) then ! Symmetric arrays with halos.
      is = G%isc - (G%isd-1) ; ie = G%iec - (G%isd-1) + 1
    else
      call SIS_error(FATAL,"post_avg: peculiar size in i-direction")
    endif

    if ((nj == G%jsc-G%jec + 1) .or. (nj == G%jsc-G%jec + 2)) then
      js = 1; je = nj  ! These arrays have no halos.
    elseif (nj == G%jed-(G%jsd-1)) then  ! Arrays have halos
      js = G%jsc - (G%jsd-1) ; je = G%jec - (G%jsd-1)
    elseif (nj == G%jed-(G%jsd-1)+1) then ! Symmetric arrays with halos.
      js = G%jsc - (G%jsd-1) ; je = G%jec - (G%jsd-1) + 1
    else
      call SIS_error(FATAL,"post_avg: peculiar size in j-direction")
    endif
  else
    is = 1; ie = ni ; js = 1 ; je = nj
  endif

  scl = 1.0 ; if (present(scale)) scl = scale
  off = 0.0 ; if (present(offset)) off = offset
  do_wt = .false. ; if (present(wtd)) do_wt = wtd

  if (do_wt) then
    avg(:,:) = 0.0 ; wts(:,:) = 0.0
    do k=1,nk ; do j=js,je ; do i=is,ie
      avg(i,j) = avg(i,j) + part(i,j,k)*(scl*val(i,j,k) + off)
      wts(i,j) = wts(i,j) + part(i,j,k)
    enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      if (wts(i,j) > 0.) then
        avg(i,j) = avg(i,j) / wts(i,j)
      else
        avg(i,j) = 0.0
      endif
    enddo ; enddo
  else
    avg(:,:) = 0.0
    do k=1,nk ; do j=js,je ; do i=is,ie
      avg(i,j) = avg(i,j) + part(i,j,k)*(scl*val(i,j,k) + off)
    enddo ; enddo ; enddo
  endif

  call post_SIS_data(id, avg, diag, mask=mask)

end subroutine post_avg_3d

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! post_avg - take area weighted average over some or all thickness categories  !
!            and offer it for diagnostic output.                               !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine post_avg_4d(id, val, part, diag, G, mask, scale, offset, wtd)
  integer, intent(in) :: id
  real, dimension(:,:,:,:), intent(in) :: val
  real, dimension(:,:,:),   intent(in) :: part
  type(SIS_diag_ctrl),      intent(in) :: diag
  type(SIS_hor_grid_type), optional, intent(in) :: G
  logical, dimension(:,:), optional, intent(in) :: mask
  real,                    optional, intent(in) :: scale, offset
  logical,                 optional, intent(in) :: wtd
  ! This subroutine determines the average of a quantity across thickness
  ! categories and does a send data on it.

  real :: avg(size(val,1),size(val,2)), wts(size(val,1),size(val,2))
  real :: scl, off, I_nLay
  logical :: do_wt
  integer :: i, j, k, L, ni, nj, nk, nLay, is, ie, js, je

  ni = size(val,1) ; nj = size(val,2) ; nk = size(val,3) ; nLay = size(val,4)
  if (size(part,1) /= ni) call SIS_error(FATAL, &
    "Mismatched i-sizes in post_avg.")
  if (size(part,2) /= nj) call SIS_error(FATAL, &
    "Mismatched j-sizes in post_avg.")
  if (size(part,3) /= nk) call SIS_error(FATAL, &
    "Mismatched k-sizes in post_avg.")

  if (present(G)) then  ! Account for the fact that arrays here start at 1.
    if ((ni == G%isc-G%iec + 1) .or. (ni == G%isc-G%iec + 2)) then
      is = 1; ie = ni  ! These arrays have no halos.
    elseif (ni == G%ied-(G%isd-1)) then  ! Arrays have halos.
      is = G%isc - (G%isd-1) ; ie = G%iec - (G%isd-1)
    elseif (ni == G%ied-(G%isd-1)+1) then ! Symmetric arrays with halos.
      is = G%isc - (G%isd-1) ; ie = G%iec - (G%isd-1) + 1
    else
      call SIS_error(FATAL,"post_avg: peculiar size in i-direction")
    endif

    if ((nj == G%jsc-G%jec + 1) .or. (nj == G%jsc-G%jec + 2)) then
      js = 1; je = nj  ! These arrays have no halos.
    elseif (nj == G%jed-(G%jsd-1)) then  ! Arrays have halos
      js = G%jsc - (G%jsd-1) ; je = G%jec - (G%jsd-1)
    elseif (nj == G%jed-(G%jsd-1)+1) then ! Symmetric arrays with halos.
      js = G%jsc - (G%jsd-1) ; je = G%jec - (G%jsd-1) + 1
    else
      call SIS_error(FATAL,"post_avg: peculiar size in j-direction")
    endif
  else
    is = 1; ie = ni ; js = 1 ; je = nj
  endif

  scl = 1.0 ; if (present(scale)) scl = scale
  off = 0.0 ; if (present(offset)) off = offset
  do_wt = .false. ; if (present(wtd)) do_wt = wtd

  if (do_wt) then
    avg(:,:) = 0.0 ; wts(:,:) = 0.0
    do L=1,nLay ; do k=1,nk ; do j=js,je ; do i=is,ie
      avg(i,j) = avg(i,j) + part(i,j,k)*(scl*val(i,j,k,L) + off)
    enddo ; enddo ; enddo ; enddo
    do k=1,nk ; do j=js,je ; do i=is,ie
      wts(i,j) = wts(i,j) + nLay * part(i,j,k)
    enddo ; enddo ; enddo
    do j=js,je ; do i=is,ie
      if (wts(i,j) > 0.) then
        avg(i,j) = avg(i,j) / wts(i,j)
      else
        avg(i,j) = 0.0
      endif
    enddo ; enddo
  else
    avg(:,:) = 0.0 ; I_nLay = 1.0/nLay
    do L=1,nLay ; do k=1,nk ; do j=js,je ; do i=is,ie
      avg(i,j) = avg(i,j) + (part(i,j,k)*I_nLay)*(scl*val(i,j,k,L) + off)
    enddo ; enddo ; enddo ; enddo
  endif

  call post_SIS_data(id, avg, diag, mask=mask)

end subroutine post_avg_4d

subroutine ice_grid_chksum(G, haloshift)
  type(SIS_hor_grid_type), optional, intent(inout) :: G
  integer, optional, intent(in) :: haloshift

  integer :: isc, iec, jsc, jec, hs
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  hs = 1 ; if (present(haloshift)) hs = haloshift

  call hchksum(G%mask2dT, "G%mask2dT", G%HI, haloshift=hs)
  call hchksum(G%geoLatT, "G%geoLatT", G%HI, haloshift=hs)
  call hchksum(G%geoLonT, "G%geoLonT", G%HI, haloshift=hs)
  call hchksum(G%dxT, "G%dxT", G%HI, haloshift=hs)
  call hchksum(G%IdxT, "G%IdxT", G%HI, haloshift=hs)
  call hchksum(G%IdyT, "G%IdyT", G%HI, haloshift=hs)
  call hchksum(G%dyT, "G%dyT", G%HI, haloshift=hs)
  call hchksum(G%areaT, "G%areaT", G%HI, haloshift=hs)
  call hchksum(G%IareaT, "G%IareaT", G%HI, haloshift=hs)
  call hchksum(G%mask2dT, "G%mask2dT", G%HI, haloshift=hs)
  call hchksum(G%cos_rot, "G%cos_rot", G%HI)
  call hchksum(G%sin_rot, "G%sin_rot", G%HI)

  call Bchksum(G%mask2dBu, "G%mask2dBu", G%HI, haloshift=hs)
  call Bchksum(G%geoLatBu, "G%geoLatBu", G%HI, haloshift=hs)
  call Bchksum(G%geoLonBu, "G%geoLonBu", G%HI, haloshift=hs)
  call vec_chksum_B("G%d[xy]Bu", G%dxBu, G%dyBu, G, halos=hs, scalars=.true.)
  call vec_chksum_B("G%Id[xy]Bu", G%IdxBu, G%IdyBu, G, halos=hs, scalars=.true.)
  call Bchksum(G%areaBu, "G%areaBu", G%HI, haloshift=hs)
  call Bchksum(G%IareaBu, "G%IareaBu", G%HI, haloshift=hs)

  call check_redundant_B("G%areaBu", G%areaBu, G, isc-1, iec+1, jsc-1, jec+1)
  call check_redundant_B("G%IareaBu", G%IareaBu, G, isc-1, iec+1, jsc-1, jec+1)

  call uchksum(G%mask2dCu, "G%mask2dCu", G%HI, haloshift=hs)
  call uchksum(G%geoLatCu, "G%geoLatCu", G%HI, haloshift=hs)
  call uchksum(G%geoLonCu, "G%geolonCu", G%HI, haloshift=hs)
  call vec_chksum_C("G%d[xy]C[uv]", G%dxCu, G%dyCv, G, halos=hs, scalars=.true.)
  call vec_chksum_C("G%d[yx]C[uv]", G%dyCu, G%dxCv, G, halos=hs, scalars=.true.)
  call vec_chksum_C("G%Id[xy]C[uv]", G%IdxCu, G%IdyCv, G, halos=hs, scalars=.true.)
  call vec_chksum_C("G%Id[yx]C[uv]", G%IdyCu, G%IdxCv, G, halos=hs, scalars=.true.)
  call uchksum(G%areaCu, "G%areaCu", G%HI, haloshift=hs)
  call uchksum(G%IareaCu, "G%IareaCu", G%HI, haloshift=hs)

  call vchksum(G%mask2dCv, "G%mask2dCv", G%HI, haloshift=hs)
  call vchksum(G%geoLatCv, "G%geoLatCv", G%HI, haloshift=hs)
  call vchksum(G%geoLonCv, "G%geoLonCv", G%HI, haloshift=hs)
  call uchksum(G%areaCu, "G%areaCv", G%HI, haloshift=hs)
  call uchksum(G%IareaCu, "G%IareaCv", G%HI, haloshift=hs)

  call hchksum(G%bathyT, "G%bathyT", G%HI, haloshift=hs)
  call Bchksum(G%CoriolisBu, "G%CoriolisBu", G%HI, haloshift=hs)
  call vec_chksum_A("G%dF_d[xy]", G%dF_dx, G%dF_dy, G, halos=hs)

end subroutine ice_grid_chksum


function is_NaN(x)
  real, intent(in) :: x
  logical :: is_nan
! This subroutine returns .true. if x is a NaN, and .false. otherwise.

  is_nan = (((x < 0.0) .and. (x >= 0.0)) .or. &
            (.not.(x < 0.0) .and. .not.(x >= 0.0)))

end function is_nan

end module SIS_utils
