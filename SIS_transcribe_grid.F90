module SIS_transcribe_grid

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

use MOM_domains, only : pass_var, pass_vector
use MOM_domains, only : To_All, SCALAR_PAIR, CGRID_NE, AGRID, BGRID_NE, CORNER
use MOM_dyn_horgrid, only : dyn_horgrid_type, set_derived_dyn_horgrid
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use SIS_hor_grid, only : SIS_hor_grid_type, set_derived_SIS_metrics

implicit none ; private

public copy_dyngrid_to_SIS_horgrid, copy_SIS_horgrid_to_dyngrid

contains

!> Copies information from a dynamic (shared) horizontal grid type into an
!! SIS_hor_grid_type.
subroutine copy_dyngrid_to_SIS_horgrid(dG, SG)
  type(dyn_horgrid_type),  intent(in)    :: dG  !< Common horizontal grid type
  type(SIS_hor_grid_type), intent(inout) :: SG  !< SIS2 horizontal grid type

  integer :: isd, ied, jsd, jed      ! Common data domains.
  integer :: IsdB, IedB, JsdB, JedB  ! Common data domains.
  integer :: ido, jdo, Ido2, Jdo2    ! Indexing offsets between the grids.
  integer :: Igst, Jgst              ! Global starting indices.
  integer :: i, j

  ! SIS_horgrid_init and create_dyn_horgrid are called outside of this routine.
  ! This routine copies over the fields that were set by MOM_initialized_fixed.

  ! Determine the indexing offsets between the grids.
  ido = dG%idg_offset - SG%idg_offset
  jdo = dG%jdg_offset - SG%jdg_offset

  isd = max(SG%isd, dG%isd+ido) ; jsd = max(SG%jsd, dG%jsd+jdo)
  ied = min(SG%ied, dG%ied+ido) ; jed = min(SG%jed, dG%jed+jdo)
  IsdB = max(SG%IsdB, dG%IsdB+ido) ; JsdB = max(SG%JsdB, dG%JsdB+jdo)
  IedB = min(SG%IedB, dG%IedB+ido) ; JedB = min(SG%JedB, dG%JedB+jdo)

  ! Check that the grids conform.
  if ((isd > SG%isc) .or. (ied < SG%ied) .or. (jsd > SG%jsc) .or. (jed > SG%jed)) &
    call MOM_error(FATAL, "copy_dyngrid_to_SIS_horgrid called with incompatible grids.")

  do i=isd,ied ; do j=jsd,jed
    SG%geoLonT(i,j) = dG%geoLonT(i+ido,j+jdo)
    SG%geoLatT(i,j) = dG%geoLatT(i+ido,j+jdo)
    SG%dxT(i,j) = dG%dxT(i+ido,j+jdo)
    SG%dyT(i,j) = dG%dyT(i+ido,j+jdo)
    SG%areaT(i,j) = dG%areaT(i+ido,j+jdo)
    SG%bathyT(i,j) = dG%bathyT(i+ido,j+jdo)

    SG%dF_dx(i,j) = dG%dF_dx(i+ido,j+jdo)
    SG%dF_dy(i,j) = dG%dF_dy(i+ido,j+jdo)
    SG%sin_rot(i,j) = dG%sin_rot(i+ido,j+jdo)
    SG%cos_rot(i,j) = dG%cos_rot(i+ido,j+jdo)
    SG%mask2dT(i,j) = dG%mask2dT(i+ido,j+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do j=jsd,jed
    SG%geoLonCu(I,j) = dG%geoLonCu(I+ido,j+jdo)
    SG%geoLatCu(I,j) = dG%geoLatCu(I+ido,j+jdo)
    SG%dxCu(I,j) = dG%dxCu(I+ido,j+jdo)
    SG%dyCu(I,j) = dG%dyCu(I+ido,j+jdo)
    SG%dy_Cu(I,j) = dG%dy_Cu(I+ido,j+jdo)
    SG%dy_Cu_obc(I,j) = dG%dy_Cu_obc(I+ido,j+jdo)

    SG%mask2dCu(I,j) = dG%mask2dCu(I+ido,j+jdo)
    SG%areaCu(I,j) = dG%areaCu(I+ido,j+jdo)
    SG%IareaCu(I,j) = dG%IareaCu(I+ido,j+jdo)
  enddo ; enddo

  do i=isd,ied ; do J=JsdB,JedB
    SG%geoLonCv(i,J) = dG%geoLonCv(i+ido,J+jdo)
    SG%geoLatCv(i,J) = dG%geoLatCv(i+ido,J+jdo)
    SG%dxCv(i,J) = dG%dxCv(i+ido,J+jdo)
    SG%dyCv(i,J) = dG%dyCv(i+ido,J+jdo)
    SG%dx_Cv(i,J) = dG%dx_Cv(i+ido,J+jdo)
    SG%dx_Cv_obc(i,J) = dG%dx_Cv_obc(i+ido,J+jdo)

    SG%mask2dCv(i,J) = dG%mask2dCv(i+ido,J+jdo)
    SG%areaCv(i,J) = dG%areaCv(i+ido,J+jdo)
    SG%IareaCv(i,J) = dG%IareaCv(i+ido,J+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do J=JsdB,JedB
    SG%geoLonBu(I,J) = dG%geoLonBu(I+ido,J+jdo)
    SG%geoLatBu(I,J) = dG%geoLatBu(I+ido,J+jdo)
    SG%dxBu(I,J) = dG%dxBu(I+ido,J+jdo)
    SG%dyBu(I,J) = dG%dyBu(I+ido,J+jdo)
    SG%areaBu(I,J) = dG%areaBu(I+ido,J+jdo)
    SG%CoriolisBu(I,J) = dG%CoriolisBu(I+ido,J+jdo)
    SG%mask2dBu(I,J) = dG%mask2dBu(I+ido,J+jdo)
  enddo ; enddo

  SG%gridLonT(SG%isg:SG%ieg) = dG%gridLonT(dG%isg:dG%ieg)
  SG%gridLatT(SG%jsg:SG%jeg) = dG%gridLatT(dG%jsg:dG%jeg)
  ! The both grids always use symmetric memory for gridLonB and gridLatB.
  SG%gridLonB(SG%isg-1:SG%ieg) = dG%gridLonB(dG%isg-1:dG%ieg)
  SG%gridLatB(SG%jsg-1:SG%jeg) = dG%gridLatB(dG%jsg-1:dG%jeg)

  ! The more complicated logic here avoids segmentation faults if one grid uses
  ! global symmetric memory while the other does not.  Because a northeast grid
  ! convention is being used, the upper bounds for each array correspond.
!  Ido2 = dG%IegB-SG%IegB ; Igst = max(SG%isg, dG%isg-Ido2)-1
!  Jdo2 = dG%JegB-SG%JegB ; Jgst = max(SG%jsg, dG%jsg-Jdo2)-1
!  do I=Igst,SG%IegB ; SG%gridLonB(I) = dG%gridLonB(I+Ido2) ; enddo
!  do J=Jgst,SG%JegB ; SG%gridLatB(J) = dG%gridLatB(J+Jdo2) ; enddo

  ! Copy various scalar variables and strings.
  SG%x_axis_units = dG%x_axis_units ; SG%y_axis_units = dG%y_axis_units
!   SG%areaT_global = dG%areaT_global ; SG%IareaT_global = dG%IareaT_global
  SG%south_lat = dG%south_lat ; SG%west_lon  = dG%west_lon
  SG%len_lat = dG%len_lat ; SG%len_lon = dG%len_lon
  SG%Rad_Earth = dG%Rad_Earth ; SG%max_depth = dG%max_depth

! Update the halos in case the dynamic grid has smaller halos than the ocean grid.
  call pass_var(SG%areaT, SG%Domain)
  call pass_var(SG%bathyT, SG%Domain)
  call pass_var(SG%geoLonT, SG%Domain)
  call pass_var(SG%geoLatT, SG%Domain)
  call pass_vector(SG%dxT, SG%dyT, SG%Domain, To_All+Scalar_Pair, AGRID)
  call pass_vector(SG%dF_dx, SG%dF_dy, SG%Domain, To_All, AGRID)
  call pass_vector(SG%cos_rot, SG%sin_rot, SG%Domain, To_All, AGRID)
  call pass_var(SG%mask2dT, SG%Domain)

  call pass_vector(SG%areaCu, SG%areaCv, SG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(SG%dyCu, SG%dxCv, SG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(SG%dxCu, SG%dyCv, SG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(SG%dy_Cu, SG%dx_Cv, SG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(SG%dy_Cu_obc, SG%dx_Cv_obc, SG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(SG%mask2dCu, SG%mask2dCv, SG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(SG%IareaCu, SG%IareaCv, SG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(SG%IareaCu, SG%IareaCv, SG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(SG%geoLatCu, SG%geoLatCv, SG%Domain, To_All+Scalar_Pair, CGRID_NE)

  call pass_var(SG%areaBu, SG%Domain, position=CORNER)
  call pass_var(SG%geoLonBu, SG%Domain, position=CORNER)
  call pass_var(SG%geoLatBu, SG%Domain, position=CORNER)
  call pass_vector(SG%dxBu, SG%dyBu, SG%Domain, To_All+Scalar_Pair, BGRID_NE)
  call pass_var(SG%CoriolisBu, SG%Domain, position=CORNER)
  call pass_var(SG%mask2dBu, SG%Domain, position=CORNER)

  call set_derived_SIS_metrics(SG)

end subroutine copy_dyngrid_to_SIS_horgrid


!> Copies information from an SIS_hor_grid_type into a dynamic (shared)
!! horizontal grid type.
subroutine copy_SIS_horgrid_to_dyngrid(SG, dG)
  type(SIS_hor_grid_type), intent(in)    :: SG  !< !< SIS2 horizontal grid type
  type(dyn_horgrid_type),  intent(inout) :: dG  !< Common horizontal grid type

  integer :: isd, ied, jsd, jed      ! Common data domains.
  integer :: IsdB, IedB, JsdB, JedB  ! Common data domains.
  integer :: ido, jdo, Ido2, Jdo2    ! Indexing offsets between the grids.
  integer :: Igst, Jgst              ! Global starting indices.
  integer :: i, j

  ! SIS_horgrid_init and create_dyn_horgrid are called outside of this routine.
  ! This routine copies over the fields that were set by MOM_initialized_fixed.

  ! Determine the indexing offsets between the grids.
  ido = SG%idG_offset - dG%idG_offset
  jdo = SG%jdG_offset - dG%jdG_offset

  isd = max(dG%isd, SG%isd+ido) ; jsd = max(dG%jsd, SG%jsd+jdo)
  ied = min(dG%ied, SG%ied+ido) ; jed = min(dG%jed, SG%jed+jdo)
  IsdB = max(dG%IsdB, SG%IsdB+ido) ; JsdB = max(dG%JsdB, SG%JsdB+jdo)
  IedB = min(dG%IedB, SG%IedB+ido) ; JedB = min(dG%JedB, SG%JedB+jdo)

  ! Check that the grids conform.
  if ((isd > dG%isc) .or. (ied < dG%ied) .or. (jsd > dG%jsc) .or. (jed > dG%jed)) &
    call MOM_error(FATAL, "copy_dyngrid_to_SIS_horgrid called with incompatible grids.")

  do i=isd,ied ; do j=jsd,jed
    dG%geoLonT(i,j) = SG%geoLonT(i+ido,j+jdo)
    dG%geoLatT(i,j) = SG%geoLatT(i+ido,j+jdo)
    dG%dxT(i,j) = SG%dxT(i+ido,j+jdo)
    dG%dyT(i,j) = SG%dyT(i+ido,j+jdo)
    dG%areaT(i,j) = SG%areaT(i+ido,j+jdo)
    dG%bathyT(i,j) = SG%bathyT(i+ido,j+jdo)

    dG%dF_dx(i,j) = SG%dF_dx(i+ido,j+jdo)
    dG%dF_dy(i,j) = SG%dF_dy(i+ido,j+jdo)
    dG%sin_rot(i,j) = SG%sin_rot(i+ido,j+jdo)
    dG%cos_rot(i,j) = SG%cos_rot(i+ido,j+jdo)
    dG%mask2dT(i,j) = SG%mask2dT(i+ido,j+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do j=jsd,jed
    dG%geoLonCu(I,j) = SG%geoLonCu(I+ido,j+jdo)
    dG%geoLatCu(I,j) = SG%geoLatCu(I+ido,j+jdo)
    dG%dxCu(I,j) = SG%dxCu(I+ido,j+jdo)
    dG%dyCu(I,j) = SG%dyCu(I+ido,j+jdo)
    dG%dy_Cu(I,j) = SG%dy_Cu(I+ido,j+jdo)
    dG%dy_Cu_obc(I,j) = SG%dy_Cu_obc(I+ido,j+jdo)

    dG%mask2dCu(I,j) = SG%mask2dCu(I+ido,j+jdo)
    dG%areaCu(I,j) = SG%areaCu(I+ido,j+jdo)
    dG%IareaCu(I,j) = SG%IareaCu(I+ido,j+jdo)
  enddo ; enddo

  do i=isd,ied ; do J=JsdB,JedB
    dG%geoLonCv(i,J) = SG%geoLonCv(i+ido,J+jdo)
    dG%geoLatCv(i,J) = SG%geoLatCv(i+ido,J+jdo)
    dG%dxCv(i,J) = SG%dxCv(i+ido,J+jdo)
    dG%dyCv(i,J) = SG%dyCv(i+ido,J+jdo)
    dG%dx_Cv(i,J) = SG%dx_Cv(i+ido,J+jdo)
    dG%dx_Cv_obc(i,J) = SG%dx_Cv_obc(i+ido,J+jdo)

    dG%mask2dCv(i,J) = SG%mask2dCv(i+ido,J+jdo)
    dG%areaCv(i,J) = SG%areaCv(i+ido,J+jdo)
    dG%IareaCv(i,J) = SG%IareaCv(i+ido,J+jdo)
  enddo ; enddo

  do I=IsdB,IedB ; do J=JsdB,JedB
    dG%geoLonBu(I,J) = SG%geoLonBu(I+ido,J+jdo)
    dG%geoLatBu(I,J) = SG%geoLatBu(I+ido,J+jdo)
    dG%dxBu(I,J) = SG%dxBu(I+ido,J+jdo)
    dG%dyBu(I,J) = SG%dyBu(I+ido,J+jdo)
    dG%areaBu(I,J) = SG%areaBu(I+ido,J+jdo)
    dG%CoriolisBu(I,J) = SG%CoriolisBu(I+ido,J+jdo)
    dG%mask2dBu(I,J) = SG%mask2dBu(I+ido,J+jdo)
  enddo ; enddo

  dG%gridLonT(dG%isg:dG%ieg) = SG%gridLonT(SG%isg:SG%ieg)
  dG%gridLatT(dG%jsg:dG%jeg) = SG%gridLatT(SG%jsg:SG%jeg)
  ! Both grids always use symmetric memory for gridLonB and gridLatB.
  dG%gridLonB(dG%isg-1:dG%ieg) = SG%gridLonB(SG%isg-1:SG%ieg)
  dG%gridLatB(dG%jsg-1:dG%jeg) = SG%gridLatB(SG%jsg-1:SG%jeg)

  ! The more complicated logic here avoids segmentation faults if one grid uses
  ! global symmetric memory while the other does not.  Because a northeast grid
  ! convention is being used, the upper bounds for each array correspond.
!  Ido2 = SG%IegB-dG%IegB ; Igst = max(dG%isg-1, (SG%isg-1)-Ido2)
!  Jdo2 = SG%JegB-dG%JegB ; Jgst = max(dG%jsg-1, (SG%jsg-1)-Jdo2)
!  do I=Igst,dG%IegB ; dG%gridLonB(I) = SG%gridLonB(I+Ido2) ; enddo
!  do J=Jgst,dG%JegB ; dG%gridLatB(J) = SG%gridLatB(J+Jdo2) ; enddo

  ! Copy various scalar variables and strings.
  dG%x_axis_units = SG%x_axis_units ; dG%y_axis_units = SG%y_axis_units
!   dG%areaT_global = SG%areaT_global ; dG%IareaT_global = SG%IareaT_global
  dG%south_lat = SG%south_lat ; dG%west_lon  = SG%west_lon
  dG%len_lat = SG%len_lat ; dG%len_lon = SG%len_lon
  dG%Rad_Earth = SG%Rad_Earth ; dG%max_depth = SG%max_depth

! Update the halos in case the dynamic grid has smaller halos than the ocean grid.
  call pass_var(dG%areaT, dG%Domain)
  call pass_var(dG%bathyT, dG%Domain)
  call pass_var(dG%geoLonT, dG%Domain)
  call pass_var(dG%geoLatT, dG%Domain)
  call pass_vector(dG%dxT, dG%dyT, dG%Domain, To_All+Scalar_Pair, AGRID)
  call pass_vector(dG%dF_dx, dG%dF_dy, dG%Domain, To_All, AGRID)
  call pass_vector(dG%cos_rot, dG%sin_rot, dG%Domain, To_All, AGRID)
  call pass_var(dG%mask2dT, dG%Domain)

  call pass_vector(dG%areaCu, dG%areaCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%dyCu, dG%dxCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%dxCu, dG%dyCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%dy_Cu, dG%dx_Cv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%dy_Cu_obc, dG%dx_Cv_obc, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%mask2dCu, dG%mask2dCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%IareaCu, dG%IareaCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%IareaCu, dG%IareaCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)
  call pass_vector(dG%geoLatCu, dG%geoLatCv, dG%Domain, To_All+Scalar_Pair, CGRID_NE)

  call pass_var(dG%areaBu, dG%Domain, position=CORNER)
  call pass_var(dG%geoLonBu, dG%Domain, position=CORNER)
  call pass_var(dG%geoLatBu, dG%Domain, position=CORNER)
  call pass_vector(dG%dxBu, dG%dyBu, dG%Domain, To_All+Scalar_Pair, BGRID_NE)
  call pass_var(dG%CoriolisBu, dG%Domain, position=CORNER)
  call pass_var(dG%mask2dBu, dG%Domain, position=CORNER)

  call  set_derived_dyn_horgrid(dG)

end subroutine copy_SIS_horgrid_to_dyngrid

end module SIS_transcribe_grid
