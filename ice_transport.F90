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
!   This module does the transport and redistribution between thickness        !
! categories for the SIS2 sea ice model.                                       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_transport_mod

use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use MOM_coms, only : reproducing_sum, EFP_type, EFP_to_real, EFP_real_diff
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING
use MOM_error_handler, only : SIS_mesg=>MOM_mesg, is_root_pe
use MOM_file_parser, only : get_param, log_param, read_param, log_version, param_file_type
use MOM_domains,     only : pass_var, pass_vector, BGRID_NE, CGRID_NE

use ice_grid_mod, only : sea_ice_grid_type
use ice_thm_mod,  only : get_thermo_coefs

implicit none ; private

#include <SIS2_memory.h>

public :: ice_transport_init, ice_transport, ice_transport_end !, ice_transport_register_restarts

type, public :: ice_transport_CS ; private

  ! parameters for doing advective and parameterized advection.
  logical :: SLAB_ICE = .false. ! should we do old style GFDL slab ice?
  real    :: chan_visc  = 0.     ! viscosity used in one-cell wide channels to parameterize transport (m^2/s)
  real    :: smag_ocn           = 0.15   ! Smagorinksy coefficient for viscosity (dimensionless)
  real    :: chan_cfl_limit     = 0.25   ! CFL limit for channel viscosity parameterization (dimensionless)
  real :: Rho_ocean = 1030.0  ! The nominal density of sea water, in kg m-3.
  real :: Rho_ice = 905.0     ! The nominal density of sea ice, in kg m-3.
  real :: Rho_snow = 330.0    ! The nominal density of snow on sea ice, in
                              ! kg m-3.
  logical :: specified_ice    ! If true, the sea ice is specified and there is
                              ! no need for ice dynamics.
  logical :: SIS1_transport   ! If true, use SIS1 code to solve the sea ice
                              ! continuity equation and transport tracers.
  logical :: check_conservation ! If true, write out verbose diagnostics of conservation.
  integer :: adv_sub_steps    ! The number of advective iterations for each slow
                              ! time step.
  type(time_type), pointer :: Time ! A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_ustar = -1, id_uocean = -1, id_uchan = -1
  integer :: id_vstar = -1, id_vocean = -1, id_vchan = -1
  integer :: id_ix_trans = -1, id_iy_trans = -1
end type ice_transport_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! transport - do ice transport and thickness class redistribution              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_transport (part_sz, h_ice, h_snow, uc, vc, t_ice, t_snow, &
                          sea_lev, hlim, dt_slow, G, CS)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT0_(G)), intent(inout) :: part_sz
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: h_ice, h_snow
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G),SZK_ICE_(G)), intent(inout) :: t_ice
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: t_snow
  real, dimension(SZIB_(G),SZJ_(G)),           intent(inout) :: uc
  real, dimension(SZI_(G),SZJB_(G)),           intent(inout) :: vc
  real, dimension(SZI_(G),SZJ_(G)),            intent(in)    :: sea_lev
  real, dimension(:),      intent(in) :: hlim  ! Move to grid type?
  real,                    intent(in) :: dt_slow
  type(ice_transport_CS), pointer :: CS
! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1, in/out.
!  (inout)   h_ice - The thickness of the ice in each category in m.
!  (inout)   h_snow - The thickness of the snow atop the ice in each category
!                     in m.
!  (in)      uc - The zonal ice velocity, in m s-1.
!  (in)      vc - The meridional ice velocity, in m s-1.
!  (inout)   t_ice - The temperature of the ice in each category and layer
!                    within the ice in degC.
!  (inout)   t_snow - The temperature of the snow atop the ice in each category
!                     in degC.
!  (in)      hlim - The lower thickness limit of each category, in m.
!  (in)      sea_lev - The height of the sea level, including contributions
!                      from non-levitating ice from an earlier time step, in m.
!  (in)      dt_slow - The amount of time over which the ice dynamics are to be
!                      advanced, in s.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.
  
  real, dimension(G%NkIce) :: pocket_enth ! Coefficients relating to the enthalpy
                                          ! contribution from melting in brine
                                          ! pockets, in degC2.
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(G)) :: &
    uh_ice, &  ! Zonal fluxes in m3 s-1 and kg s-1.
    uh_snow    ! Zonal fluxes in m3 s-1 and kg s-1.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    uf0, uf, & ! Zonal fluxes in m3 s-1 and kg s-1.
    ustar, ustaro, ustarv ! Local variables, transporting velocities
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(G)) :: &
    vh_ice, &  ! Meridional fluxes in m3 s-1 and kg s-1.
    vh_snow    ! Meridional fluxes in m3 s-1 and kg s-1.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    vf0, vf, & ! Meridional fluxes in m3 s-1 and kg s-1.
    vstar, vstaro, vstarv ! Local variables, transporting velocities
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)) :: &
    vol_ice, vol_snow, &  ! The total volume of snow and ice per unit area in a cell, in m.
    vol0_ice, vol0_snow   ! The initial total volume of snow and ice per unit area in a cell, in m.
  
!  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)) :: &
!    vol_ice_d1, h_ice_d1, vol_snow_d1, h_snow_d1, T_snow_d1, &
!    vol_ice_d2, h_ice_d2, vol_snow_d2, h_snow_d2, T_snow_d2
!  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G),SZK_ICE_(G)) :: &
!    T_ice_d1, T_ice_d2
  
  real, dimension(SZI_(G),SZJ_(G)) :: ice_cover ! The summed fractional ice concentration, ND.
  real :: u_visc, u_ocn, cnn, grad_eta ! Variables for channel parameterization
  type(EFP_type) :: tot_ice(2), tot_snow(2), enth_ice(2), enth_snow(2)
  real :: I_tot_ice, I_tot_snow

  real :: dt_adv
  character(len=200) :: mesg
  integer :: i, j, k, l, bad, isc, iec, jsc, jec, isd, ied, jsd, jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (CS%slab_ice) then
    call pass_vector(uc, vc, G%Domain, stagger=CGRID_NE)
    call slab_ice_advect(uc, vc, h_ice(:,:,1), 4.0, dt_slow, G, CS)
    call pass_var(h_ice(:,:,2), G%Domain)
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      if (h_ice(i,j,1) > 0.0) then
        part_sz(i,j,1) = 1.0
      else 
        part_sz(i,j,1) = 0.0
      endif
    enddo ; enddo
    return
  endif


  if (CS%chan_visc>0. .and. CS%adv_sub_steps>0) then
  ! This block of code is a parameterization of either (or both)
  ! i) the pressure driven oceanic flow in a narrow channel, or
  ! ii) pressure driven flow of ice itself.
  ! The latter is a speculative but both are missing due to
  ! masking of velocities to zero in a single-cell wide channel.
    dt_adv = dt_slow/CS%adv_sub_steps

    ice_cover(:,:) = min(sum(part_sz(:,:,1:G%CatIce),dim=3),1.0)
    ustar(:,:)=0.; vstar(:,:)=0.
    ustaro(:,:)=0.; vstaro(:,:)=0.
    ustarv(:,:)=0.; vstarv(:,:)=0.
    do j=jsc,jec ; do I=isc-1,iec
      if ((uc(I,j)==0.) .and. & ! this is a redundant test due to following line
          (G%mask2dBu(I,J)+G%mask2dBu(I,J-1)==0.) .and. &  ! =0 => no vels
          (G%mask2dT(i,j)*G%mask2dT(i+1,j)>0.)) then ! >0 => open for transport
        grad_eta=(sea_lev(i+1,j)-sea_lev(i,j)) * G%IdxCu(I,j) ! delta_i eta/dx
        u_visc=-G%g_Earth*((G%dyCu(I,j)*G%dyCu(I,j))/(12.*CS%chan_visc)) & ! -g*dy^2/(12*visc)
                 *grad_eta                                                  ! d/dx eta
        u_ocn=sqrt( G%g_Earth*G%dyCu(I,j)*abs(grad_eta)/(36.*CS%smag_ocn) ) ! Magnitude of ocean current
        u_ocn=sign(u_ocn, -grad_eta) ! Direct down the ssh gradient
        cnn=max(ice_cover(i,j),ice_cover(i+1,j))**2. ! Use the larger concentration
        uc(I,j)=cnn*u_visc+(1.-cnn)*u_ocn
        ! Limit flow to be stable for fully divergent flow
        if (uc(I,j)>0.) then
          uc(I,j)=min( uc(I,j), CS%chan_cfl_limit*G%dxT(i,j)/dt_adv)
        else
          uc(I,j)=max( uc(I,j),(-1*CS%chan_cfl_limit)*G%dxT(i+1,j)/dt_adv)
        endif
        if (CS%id_ustar>0) ustar(I,j)=uc(I,j)
        if (CS%id_uocean>0) ustaro(I,j)=u_ocn
        if (CS%id_uchan>0) ustarv(I,j)=u_visc
      endif
    enddo ; enddo
    do J=jsc-1,jec ; do i=isc,iec
      if ((vc(i,J)==0.) .and. & ! this is a redundant test due to following line
          (G%mask2dBu(I,J)+G%mask2dBu(I-1,J)==0.) .and. &  ! =0 => no vels
          (G%mask2dT(i,j)*G%mask2dT(i,j+1)>0.)) then ! >0 => open for transport
        grad_eta=(sea_lev(i,j+1)-sea_lev(i,j)) * G%IdyCv(i,J) ! delta_i eta/dy
        u_visc=-G%g_Earth*((G%dxCv(i,J)*G%dxCv(i,J))/(12.*CS%chan_visc)) & ! -g*dy^2/(12*visc)
                *grad_eta                                                  ! d/dx eta
        u_ocn=sqrt( G%g_Earth*G%dxCv(i,J)*abs(grad_eta)/(36.*CS%smag_ocn) ) ! Magnitude of ocean current
        u_ocn=sign(u_ocn, -grad_eta) ! Direct down the ssh gradient
        cnn=max(ice_cover(i,j),ice_cover(i,j+1))**2. ! Use the larger concentration
        vc(i,J)=cnn*u_visc+(1.-cnn)*u_ocn
        ! Limit flow to be stable for fully divergent flow
        if (vc(i,J)>0.) then
          vc(i,J)=min( vc(i,J), CS%chan_cfl_limit*G%dyT(i,j)/dt_adv)
        else
          vc(i,J)=max( vc(i,J),(-1*CS%chan_cfl_limit)*G%dyT(i,j+1)/dt_adv)
        endif
        if (CS%id_vstar>0) vstar(i,J)=vc(i,J)
        if (CS%id_vocean>0) vstaro(i,J)=u_ocn
        if (CS%id_vchan>0) vstarv(i,J)=u_visc
      endif
    enddo ; enddo
  endif

  call pass_vector(uc, vc, G%Domain, stagger=CGRID_NE)

  call get_thermo_coefs(layer_coefs=pocket_enth)

  if (CS%check_conservation) then
    call get_total_enthalpy(h_ice, h_snow, part_sz, t_ice, t_snow, &
                            pocket_enth, G, enth_ice(1), enth_snow(1))
  endif

  !   Convert from ice temperature (which is not conserved) to enthalpy, which
  ! includes the heat requirements for melting of brine pockets associated with 
  ! temperature changes.  These expressions stem from the assumptions that
  ! brine pockets will shrink or grow until their salinity gives a freezing
  ! point that matches the local temperature.
  ! This was prevously the subroutine thm_pack.
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    if (h_ice(i,j,k)>0.0) then
      vol_ice(i,j,k) = part_sz(i,j,k)*h_ice(i,j,k)
      vol_snow(i,j,k) = part_sz(i,j,k)*h_snow(i,j,k)
      if (CS%SIS1_transport) then
        do l=1,G%NkIce
          t_ice(i,j,k,l) = (t_ice(i,j,k,l) - pocket_enth(l)/t_ice(i,j,k,l)) * &
                            vol_ice(i,j,k)
        enddo
        t_snow(i,j,k) = t_snow(i,j,k)*vol_snow(i,j,k)
      else
        do l=1,G%NkIce
          t_ice(i,j,k,l) = (t_ice(i,j,k,l) - pocket_enth(l)/t_ice(i,j,k,l))
        enddo
      endif
    else
      if (part_sz(i,j,k)*h_snow(i,j,k) > 0.0) then
        call SIS_error(FATAL, "Input to ice_transport, non-zero snow mass rests atop no ice.")
      endif
      part_sz(i,j,k) = 0.0 ; vol_ice(i,j,k) = 0.0
      vol_snow(i,j,k) = 0.0

      do l=1,G%NkIce ; t_ice(i,j,k,l) = 0.0 ; enddo
      t_snow(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo

  if (CS%check_conservation) then
    call get_total_amounts(vol_ice, vol_snow, G, tot_ice(1), tot_snow(1))
  endif
  
  call pass_var(part_sz, G%Domain) ! cannot be combined with updates below
  do l=1,G%NkIce ! The do loop allows these to be combined with the following.
    call pass_var(t_ice(:,:,:,l), G%Domain, complete=.false.)
  enddo
  call pass_var(t_snow, G%Domain, complete=.false.)
  call pass_var(vol_ice,  G%Domain, complete=.false.)
  call pass_var(vol_snow, G%Domain, complete=.false.)
  call pass_var(h_ice, G%Domain, complete=.true.)

  uf(:,:) = 0.0; vf(:,:) = 0.0
  if (CS%SIS1_transport) then
    do k=1,G%CatIce
      call ice_advect(uc, vc, part_sz(:,:,k), dt_slow, G, CS)
      call ice_advect(uc, vc, vol_snow(:,:,k), dt_slow, G, CS, uf0, vf0)
      uf(:,:) = uf(:,:) + CS%Rho_snow*uf0(:,:)
      vf(:,:) = vf(:,:) + CS%Rho_snow*vf0(:,:)
      call ice_advect(uc, vc, vol_ice(:,:,k), dt_slow, G, CS, uf0, vf0)
      uf(:,:) = uf(:,:) + CS%Rho_ice*uf0(:,:)
      vf(:,:) = vf(:,:) + CS%Rho_ice*vf0(:,:)
      call ice_advect(uc, vc, t_snow(:,:,k), dt_slow, G, CS)
      do l=1,G%NkIce
        call ice_advect(uc, vc, t_ice(:,:,k,l), dt_slow, G, CS)
      enddo
    enddo

    do j=jsc,jec ; do i=isc,iec
      if (sum(vol_ice(i,j,:))>0) &
        call ice_redistribute(part_sz(i,j,1:G%CatIce), G, &
           vol_snow(i,j,:), t_snow(i,j,:), vol_ice (i,j,:), &
           t_ice(i,j,:,1), t_ice(i,j,:,2), &
           t_ice(i,j,:,3), t_ice(i,j,:,4), hlim)
    enddo ; enddo

    !   Convert from enthalpy back to ice temperature. These expressions stem from
    ! the assumptions that brine pockets will shrink or grow until their salinity
    ! gives a freezing point that matches the local temperature.
    ! This was prevously the subroutine thm_unpack.
    do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
      if (part_sz(i,j,k)<1e-10) vol_ice(i,j,k) = 0.0
      if (vol_ice(i,j,k)>0.0) then
        do l=1,G%NkIce
          t_ice(i,j,k,l) = t_ice(i,j,k,l) / vol_ice(i,j,k)
          t_ice(i,j,k,l) = 0.5*(t_ice(i,j,k,l) - &
                                sqrt(t_ice(i,j,k,l)*t_ice(i,j,k,l) + 4.0*pocket_enth(l)))
        enddo
        h_ice(i,j,k) = vol_ice(i,j,k)/part_sz(i,j,k)
        if (vol_snow(i,j,k)>0.0) then
          t_snow(i,j,k) = t_snow(i,j,k)/vol_snow(i,j,k)
        else
          t_snow(i,j,k) = 0.0
        endif
        h_snow(i,j,k) = vol_snow(i,j,k)/part_sz(i,j,k)
      else
        part_sz(i,j,k) = 0.0 ; h_ice(i,j,k) = 0.0
        do l=1,G%NkIce ; t_ice(i,j,k,l) = 0.0 ; enddo
        h_snow(i,j,k) = 0.0
        t_snow(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo
  else ! Not SIS1_transport.
    ! Do the transport via the continuity equations and tracer conservation
    ! equations for h_ice and tracers, inverting for the fractional size of
    ! each partition.
    

    ! Part-size no longer matters, but make sure that ice is in the right thickness
    ! category before advection.

    ! Copy the arrays for debugging purposes.
    ! ###
!    vol_ice_d1(:,:,:) = vol_ice(:,:,:)
!    h_ice_d1(:,:,:) = h_ice(:,:,:)
!    vol_snow_d1(:,:,:) = vol_snow(:,:,:)
!    h_snow_d1(:,:,:) = h_snow(:,:,:)
!    T_snow_d1(:,:,:) = T_snow(:,:,:)
!    T_ice_d1(:,:,:,:) = T_ice(:,:,:,:)

    call adjust_ice_categories(hlim, vol_ice, vol_snow, h_ice, t_ice, t_snow, G, CS)

    ! Copy the arrays for debugging purposes.
    ! ###
!    vol_ice_d2(:,:,:) = vol_ice(:,:,:)
!    h_ice_d2(:,:,:) = h_ice(:,:,:)
!    vol_snow_d2(:,:,:) = vol_snow(:,:,:)
!    h_snow_d2(:,:,:) = h_snow(:,:,:)
!    T_snow_d2(:,:,:) = T_snow(:,:,:)
!    T_ice_d2(:,:,:,:) = T_ice(:,:,:,:)

    do k=1,G%CatIce ; do j=jsd,jed ; do i=isd,ied
      vol0_ice(i,j,k) = vol_ice(i,j,k)
      vol0_snow(i,j,k) = vol_snow(i,j,k)
    enddo ; enddo ; enddo
    call ice_continuity(uc, vc, vol0_ice, vol_ice, uh_ice, vh_ice, dt_slow, G, CS)
    call ice_continuity(uc, vc, vol0_snow, vol_snow, uh_snow, vh_snow, dt_slow, G, CS)

    call advect_ice_tracer(vol0_ice, vol_ice, uh_ice, vh_ice, h_ice, dt_slow, G, CS)
    do l=1,G%NkIce
      call advect_ice_tracer(vol0_ice, vol_ice, uh_ice, vh_ice, t_ice(:,:,:,l), dt_slow, G, CS)
    enddo
    call advect_ice_tracer(vol0_snow, vol_snow, uh_snow, vh_snow, t_snow, dt_slow, G, CS)

    ice_cover(:,:) = 0.0
    do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
      if (vol_ice(i,j,k) > 0.0) then
        part_sz(i,j,k) = vol_ice(i,j,k) / h_ice(i,j,k)
        h_snow(i,j,k) = vol_snow(i,j,k) / part_sz(i,j,k)
!### Consider...        h_snow(i,j,k) = h_ice(i,j,k) * vol_snow(i,j,k) / vol_ice(i,j,k)
        ice_cover(i,j) = ice_cover(i,j) + part_sz(i,j,k)
      else
        part_sz(i,j,k) = 0.0 ; h_ice(i,j,k) = 0.0
        if (vol_snow(i,j,k) > 0.0) &
          call SIS_error(FATAL, &
            "Positive snow volumes should not exist without ice.")
        h_snow(i,j,k) = 0.0
      endif  
    enddo ; enddo ; enddo
    do j=jsc,jec ; do i=isc,iec
      part_sz(i,j,0) = 1.0-ice_cover(i,j)
    enddo ; enddo

    ! Compress the ice where the fractional coverage exceeds 1, starting with
    ! the thinnest categories.
    call compress_ice(part_sz, hlim, vol_ice, vol_snow, h_ice, h_snow, t_ice, t_snow, G, CS)

    if ((CS%id_ix_trans>0) .or. (CS%id_iy_trans>0)) then ; do k=1,G%CatIce 
      do j=jsc,jec ; do I=isc-1,iec
        uf(I,j) = uf(I,j) + (CS%Rho_snow*uh_snow(I,j,k) + CS%Rho_ice*uh_ice(I,j,k))
      enddo ; enddo
      do J=jsc-1,jec ; do i=isc,iec
        vf(i,J) = vf(i,J) + (CS%Rho_snow*vh_snow(i,J,k) + CS%Rho_ice*vh_ice(i,J,k))
      enddo ; enddo
    enddo ; endif

    !   Convert from enthalpy back to ice temperature. These expressions stem from
    ! the assumptions that brine pockets will shrink or grow until their salinity
    ! gives a freezing point that matches the local temperature.
    do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
      if (vol_ice(i,j,k)>0.0) then
        do l=1,G%NkIce
          t_ice(i,j,k,l) = 0.5*(t_ice(i,j,k,l) - &
                                sqrt(t_ice(i,j,k,l)*t_ice(i,j,k,l) + 4.0*pocket_enth(l)))
        enddo
      else
        part_sz(i,j,k) = 0.0 ; h_ice(i,j,k) = 0.0
        do l=1,G%NkIce ; t_ice(i,j,k,l) = 0.0 ; enddo
        h_snow(i,j,k) = 0.0
        t_snow(i,j,k) = 0.0
      endif
    enddo ; enddo ; enddo

  endif ! Not SIS1_transport.

  call pass_var(part_sz, G%Domain) ! cannot be combined with the two updates below
  call pass_var(h_snow, G%Domain, complete=.false.)
  call pass_var(h_ice, G%Domain, complete=.true.)

  ! Recalculate part_sz(:,:,0) to ensure that the sum of part_sz adds up to 1.
  part_sz(:,:,0) = 1.0
  do k=1,G%CatIce ; part_sz(:,:,0) = part_sz(:,:,0) - part_sz(:,:,k) ; enddo
!  ice_cover(:,:) = 0.0
!  do k=1,G%CatIce ; ice_cover(:,:) = ice_cover(:,:) + part_sz(:,:,k) ; enddo
!  part_sz(:,:,0) = 1.0 - ice_cover(:,:)

  if (CS%check_conservation) then
    do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
      vol_ice(i,j,k) = part_sz(i,j,k)*h_ice(i,j,k)
      vol_snow(i,j,k) = part_sz(i,j,k)*h_snow(i,j,k)
    enddo ; enddo ; enddo
    
    call get_total_amounts(vol_ice, vol_snow, G, tot_ice(2), tot_snow(2))

    call get_total_enthalpy(h_ice, h_snow, part_sz, t_ice, t_snow, &
                            pocket_enth, G, enth_ice(2), enth_snow(2))
    
    if (is_root_pe()) then
      I_tot_ice  = abs(EFP_to_real(tot_ice(1)))
      if (I_tot_ice > 0.0) I_tot_ice = 1.0 / I_tot_ice    ! Adcroft's rule inverse.
      I_tot_snow = abs(EFP_to_real(tot_snow(1)))
      if (I_tot_snow > 0.0) I_tot_snow = 1.0 / I_tot_snow ! Adcroft's rule inverse.
      write(*,'("  Total Ice vol:  ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(tot_ice(2)), EFP_real_diff(tot_ice(2),tot_ice(1)), &
        EFP_real_diff(tot_ice(2),tot_ice(1)) * I_tot_ice
      write(*,'("  Total Snow vol: ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(tot_snow(2)), EFP_real_diff(tot_snow(2),tot_snow(1)), &
        EFP_real_diff(tot_snow(2),tot_snow(1)) * I_tot_snow


      I_tot_ice  = abs(EFP_to_real(enth_ice(1)))
      if (I_tot_ice > 0.0) I_tot_ice = 1.0 / I_tot_ice    ! Adcroft's rule inverse.
      I_tot_snow = abs(EFP_to_real(enth_snow(1)))
      if (I_tot_snow > 0.0) I_tot_snow = 1.0 / I_tot_snow ! Adcroft's rule inverse.
      write(*,'("  Enthalpy Ice:  ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(enth_ice(2)), EFP_real_diff(enth_ice(2),enth_ice(1)), &
        EFP_real_diff(enth_ice(2),enth_ice(1)) * I_tot_ice
      write(*,'("  Enthalpy Snow: ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(enth_snow(2)), EFP_real_diff(enth_snow(2),enth_snow(1)), &
        EFP_real_diff(enth_snow(2),enth_snow(1)) * I_tot_snow
    endif
  endif

  ! Check for bad values of thickness and temperature.
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    ! Check for bad values.
    bad = 0
    if (h_snow(i,j,k) <0.0 .or. h_snow(i,j,k) > 1e3 ) &
      bad = bad+1
    if (h_ice(i,j,k) <0.0 .or. h_ice(i,j,k) > 1e3 ) &
      bad = bad+1
    if ((bad == 0) .and. (h_ice(i,j,k) == 0.0)) cycle
    if (t_snow(i,j,k)>0.0.or.t_snow(i,j,k)<-100.0) &
      bad = bad+1
    do l=1,G%NkIce
      if (t_ice(i,j,k,l) > 0.0 .or. t_ice(i,j,k,l) < -100.0) &
        bad = bad+1
    enddo      
!### USE BETTER ERROR HANDLING LATER.
    if (bad>0) &
      print *, 'BAD ICE AFTER UNPACK ', 'hs/hi=',h_snow(i,j,k),h_ice(i,j,k),'tsn/tice=', &
                        t_snow(i,j,k),t_ice(i,j,k,1),t_ice(i,j,k,2),t_ice(i,j,k,3), &
                        t_ice(i,j,k,4),'cn=',part_sz(i,j,k)
  enddo ; enddo ; enddo
  
  if (CS%id_ix_trans>0) call post_SIS_data(CS%id_ix_trans, uf, CS%diag)
  if (CS%id_iy_trans>0) call post_SIS_data(CS%id_iy_trans, vf, CS%diag)
  if (CS%id_ustar >0) call post_SIS_data(CS%id_ustar, ustar, CS%diag)
  if (CS%id_vstar >0) call post_SIS_data(CS%id_vstar , vstar, CS%diag)
  if (CS%id_vocean>0) call post_SIS_data(CS%id_vocean, vstaro, CS%diag)
  if (CS%id_uocean>0) call post_SIS_data(CS%id_uocean, ustaro, CS%diag)
  if (CS%id_vchan>0)  call post_SIS_data(CS%id_vchan,  vstarv, CS%diag)
  if (CS%id_uchan>0)  call post_SIS_data(CS%id_uchan,  ustarv, CS%diag)

end subroutine ice_transport

subroutine ice_continuity(u, v, h_in, h, uh, vh, dt, G, CS)
  type(sea_ice_grid_type),                     intent(inout) :: G
  real, dimension(SZIB_(G),SZJ_(G)),           intent(in)    :: u
  real, dimension(SZI_(G),SZJB_(G)),           intent(in)    :: v
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(in)    :: h_in
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(out)   :: h
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(G)), intent(out)   :: uh
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(G)), intent(out)   :: vh
  real,                                        intent(in)    :: dt
  type(ice_transport_CS),                      pointer       :: CS
!    This subroutine time steps the category thicknesses averaged over the whole
!  grid cell, initially using directionally unsplit using upwind advection.  But
!  in subsequent versions this will be replaced with a directionally split 
!  algorithm derived from MOM_continuity_PPM.F90.  In the following
!  documentation, H is used for the units of thickness (usually m or kg m-2.)

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h_in - Initial layer thickness, in H.
!  (out)     h - Final layer thickness, in H.
!  (out)     uh - Volume flux through zonal faces = u*h*dy, H m2 s-1.
!  (out)     vh - Volume flux through meridional faces = v*h*dx,
!                  in H m2 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 ice_transport_init.
  real h_up
  integer :: i, j, k, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do k=1,G%CatIce ; do j=js,je ; do I=is-1,ie
    if (u(I,j) >= 0.0) then ; h_up = h(i,j,k)
    else ; h_up = h(i+1,j,k) ; endif
!###   uh(I,j,k) = G%dy_Cu(I,j) * u(I,j) * h_up
    uh(I,j,k) = G%dyCu(I,j) * u(I,j) * h_up
  enddo ; enddo ; enddo

  do k=1,G%CatIce ; do J=js-1,je ; do i=is,ie
    if (v(i,J) >= 0.0) then ; h_up = h(i,j,k)
    else ; h_up = h(i,j+1,k) ; endif
!###   vh(i,J,k) = G%dy_Cv(i,J) * v(i,J) * h_up
    vh(i,J,k) = G%dyCv(i,J) * v(i,J) * h_up
  enddo ; enddo ; enddo

  do k=1,G%CatIce ; do j=js,je ; do i=is,ie
    h(i,j,k) = h_in(i,j,k) - dt* G%IareaT(i,j) * &
         ((uh(I,j,k) - uh(I-1,j,k)) + (vh(i,J,k) - vh(i,J-1,k)))

!   This line should be unnecessary.
    if (h(i,j,k) < 0.0) h(i,j,k) = 0.0
  enddo ; enddo ; enddo

end subroutine ice_continuity

subroutine advect_ice_tracer(h_prev, h_end, uhtr, vhtr, tr, dt, G, CS) !, Reg)
  type(sea_ice_grid_type),                     intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(in)    :: h_prev, h_end
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(G)), intent(in)    :: uhtr
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(G)), intent(in)    :: vhtr
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: tr
  real,                                        intent(in)    :: dt
  type(ice_transport_CS),                      pointer       :: CS
!  type(tracer_registry_type),                 pointer       :: Reg
!    This subroutine time steps the tracer concentrations.
!  A monotonic, conservative, weakly diffusive scheme will eventually used, but
!  for now a simple upwind scheme is used.

! Arguments: h_prev - Category thickness before advection, in m or kg m-2.
!  (in)      h_end - Layer thickness after advection, in m or kg m-2.
!  (in)      uhtr - Accumulated volume or mass fluxes through zonal faces,
!                   in m3 or kg.
!  (in)      vhtr - Accumulated volume or mass fluxes through meridional faces,
!                   in m3 or kg.
!    (inout) tr - The tracer concentration being worked on (###Move to registry.)
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 tracer_advect_init.
!  (in)      Reg - A pointer to the tracer registry.

  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(G)) :: flux_x
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(G)) :: flux_y
  real    :: tr_up
  real    :: Ihnew
  integer :: i, j, k, m, is, ie, js, je !, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
!  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Reconstruct the old value of h ???
  ! if (h_prev(i,j,k) > 0.0) then
  ! h_last(i,j,k) = h_end(i,j,k) + dt * G%IareaT(i,j) * &
  !        ((uh(I,j,k) - uh(I-1,j,k)) + (vh(i,J,k) - vh(i,J-1,k)))

  ! For now this is just non-directionally split upwind advection.
  do k=1,G%CatIce
    do j=js,je ; do I=is-1,ie
      if (uhtr(I,j,k) >= 0.0) then ; tr_up = tr(i,j,k)
      else ; tr_up = tr(i+1,j,k) ; endif
      flux_x(I,j,k) = uhtr(I,j,k) * tr_up
    enddo ; enddo

    do J=js-1,je ; do i=is,ie
      if (vhtr(i,J,k) >= 0.0) then ; tr_up = tr(i,j,k)
      else ; tr_up = tr(i,j+1,k) ; endif
      flux_y(i,J,k) = vhtr(i,J,k) * tr_up
    enddo ; enddo

    do j=js,je ; do i=is,ie
      Ihnew = 0.0
      if (h_end(i,j,k) > 0.0) Ihnew = 1.0 / h_end(i,j,k)
      tr(i,j,k) = ( h_prev(i,j,k)*tr(i,j,k) - dt* G%IareaT(i,j) * &
           ((flux_x(I,j,k) - flux_x(I-1,j,k)) + &
            (flux_y(i,J,k) - flux_y(i,J-1,k))) ) * Ihnew
    enddo ; enddo
  enddo

end subroutine advect_ice_tracer

subroutine adjust_ice_categories(hlim, vol_ice, vol_snow, h_ice, t_ice, t_snow, G, CS)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(:),                          intent(in)    :: hlim  ! Move to grid type?
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: vol_ice, vol_snow, h_ice
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G),SZK_ICE_(G)), intent(inout) :: t_ice
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: t_snow
  type(ice_transport_CS),                      pointer       :: CS

!   This subroutine moves mass between thickness categories if it is thinner or
! thicker than the bounding limits of each category.

! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      hlim - The lower thickness limit of each category, in m.
!  (inout)   vol_ice - The volume per unit grid-cell area of the ice in each
!                      category in m.
!  (inout)   vol_snow - The volume per unit grid-cell area of the snow atop the
!                       ice in each category in m.
!  (inout)   h_ice - The thickness of the ice in each category in m.
!  (inout)   t_ice - The temperature of the ice in each category and layer
!                    within the ice in degC.
!  (inout)   t_snow - The temperature of the snow atop the ice in each category
!                     in degC.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.
  real :: vol_trans
  real :: snow_trans
  real :: Ivol_snow
  real :: Ivol_new
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do k=1,G%CatIce-1 ; do j=jsd,jed ; do i=isd,ied
    if ((vol_ice(i,j,k) > 0.0) .and. (h_ice(i,j,k) > hlim(k+1))) then
      ! Move some or all of the ice to a thicker category.
      ! For now move all of it.
      vol_trans = vol_ice(i,j,k) ; Ivol_new = 1.0 / (vol_trans + vol_ice(i,j,k+1))
      snow_trans = vol_snow(i,j,k) ! * (vol_trans / vol_ice) = 1

      ! This is upwind advection across categories.  Improve it later.
      do m=1,G%NkIce
        T_ice(i,j,k+1,m) = (vol_trans*T_ice(i,j,k,m) + &
                            vol_ice(i,j,k+1)*T_ice(i,j,k+1,m)) * Ivol_new
      enddo
      h_ice(i,j,k+1) = (vol_trans*h_ice(i,j,k) + &
                        vol_ice(i,j,k+1)*h_ice(i,j,k+1)) * Ivol_new
      ! h should be the first thing to correct via a non-constant profile, and
      ! can be improved independent of T & S.
      h_ice(i,j,k) = hlim(k+1)

      vol_ice(i,j,k+1) = vol_ice(i,j,k+1) + vol_trans
      vol_ice(i,j,k) = vol_ice(i,j,k) - vol_trans 

      if (snow_trans > 0.0) then
        Ivol_snow = 1.0 / (snow_trans + vol_snow(i,j,k+1))
        t_snow(i,j,k+1) = (snow_trans*T_snow(i,j,k) + &
                           vol_snow(i,j,k+1)*T_snow(i,j,k+1)) * Ivol_snow
        vol_snow(i,j,k+1) = vol_snow(i,j,k+1) + snow_trans
        vol_snow(i,j,k) = vol_snow(i,j,k) - snow_trans
      endif
    endif
  enddo ; enddo ; enddo

  do k=G%CatIce,2,-1 ; do j=jsd,jed ; do i=isd,ied
    if ((vol_ice(i,j,k) > 0.0) .and. (h_ice(i,j,k) < hlim(k))) then
      ! Move some or all of the ice to a thinner category.
      ! For now move all of it.
      vol_trans = vol_ice(i,j,k) ; Ivol_new = 1.0 / (vol_trans + vol_ice(i,j,k-1))
      snow_trans = vol_snow(i,j,k) ! * (vol_trans / vol_ice) = 1

      ! This is upwind advection across categories.  Improve it later.
      do m=1,G%NkIce
        T_ice(i,j,k-1,m) = (vol_trans*T_ice(i,j,k,m) + &
                            vol_ice(i,j,k-1)*T_ice(i,j,k-1,m)) * Ivol_new
      enddo
      h_ice(i,j,k-1) = (vol_trans*h_ice(i,j,k) + &
                        vol_ice(i,j,k-1)*h_ice(i,j,k-1)) * Ivol_new
      ! h should be the first thing to correct via a non-constant profile, and
      ! can be improved independent of T & S.
      h_ice(i,j,k) = hlim(k)

      vol_ice(i,j,k-1) = vol_ice(i,j,k-1) + vol_trans
      vol_ice(i,j,k) = vol_ice(i,j,k) - vol_trans 

      if (snow_trans > 0.0) then
        Ivol_snow = 1.0 / (snow_trans + vol_snow(i,j,k-1))
        t_snow(i,j,k-1) = (snow_trans*T_snow(i,j,k) + &
                           vol_snow(i,j,k-1)*T_snow(i,j,k-1)) * Ivol_snow
        vol_snow(i,j,k-1) = vol_snow(i,j,k-1) + snow_trans
        vol_snow(i,j,k) = vol_snow(i,j,k) - snow_trans
      endif
    endif
  enddo ; enddo ; enddo

end subroutine adjust_ice_categories

subroutine compress_ice(part_sz, hlim, vol_ice, vol_snow, h_ice, h_snow, t_ice, t_snow, G, CS)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT0_(G)), intent(inout) :: part_sz
  real, dimension(:),                          intent(in)    :: hlim  ! Move to grid type?
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: vol_ice, vol_snow, h_ice, h_snow
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G),SZK_ICE_(G)), intent(inout) :: t_ice
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(inout) :: t_snow
  type(ice_transport_CS),                      pointer       :: CS
!   This subroutine compresses the ice, starting with the thinnest category, if
! the total fractional ice coverage exceeds 1.  It is assumed at the start that
! the sum over all categories (including ice free) of part_sz is 1, but that the
! part_sz of the ice free category may be negative to make this so.  In this
! routine, the volume (mass) is conserved, while the fractional coverage is
! solved for, while the new thicknesses are diagnosed.

! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      hlim - The lower thickness limit of each category, in m.
!  (inout)   vol_ice - The volume per unit grid-cell area of the ice in each
!                      category in m.
!  (inout)   vol_snow - The volume per unit grid-cell area of the snow atop the
!                       ice in each category in m.
!  (inout)   h_ice - The thickness of the ice in each category in m.
!  (inout)   h_snow - The thickness of the snow atop the ice in each category
!                     in m.
!  (inout)   t_ice - The temperature of the ice in each category and layer
!                    within the ice in degC.
!  (inout)   t_snow - The temperature of the snow atop the ice in each category
!                     in degC.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.

  real, dimension(SZI_(G),SZJ_(G)) :: excess_cover
  real :: compression_ratio
  real :: Icompress_here
  real :: vol_trans, vol_old
  real :: snow_trans, snow_old
  real :: Ivol_snow, Ivol_new
  logical :: do_j(SZJ_(G))
  character(len=200) :: mesg
  integer :: i, j, k, m, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  do_j(:) = .false.
  do j=jsc,jec ; do i=isc,iec
    if (part_sz(i,j,0) < 0.0) then
      excess_cover(i,j) = -part_sz(i,j,0) ; part_sz(i,j,0) = 0.0
      do_j(j) = .true.
    else
      excess_cover(i,j) = 0.0
    endif
  enddo ; enddo

  do j=jsc,jec ; if (do_j(j)) then ; do k=1,G%CatIce-1 ; do i=isc,iec
    if ((excess_cover(i,j) > 0.0) .and. (vol_ice(i,j,k) > 0.0)) then
      compression_ratio = h_ice(i,j,k) / hlim(k+1)
      if (part_sz(i,j,k)*(1.0-compression_ratio) >= excess_cover(i,j)) then
        ! This category is compacted, but not to the point that it needs to
        ! be transferred to a thicker layer.
        Icompress_here = part_sz(i,j,k) / (part_sz(i,j,k) - excess_cover(i,j))
        h_ice(i,j,k) = h_ice(i,j,k) * Icompress_here
        h_snow(i,j,k) = h_snow(i,j,k) * Icompress_here
        part_sz(i,j,k) = part_sz(i,j,k) - excess_cover(i,j)
        excess_cover(i,j) = 0.0
      else
        ! Mass from this category needs to be transfered to the next thicker
        ! category after being compacted to thickness hlim(k+1).
        excess_cover(i,j) = excess_cover(i,j) - part_sz(i,j,k)*(1.0-compression_ratio)
        part_sz(i,j,k+1) = part_sz(i,j,k+1) + part_sz(i,j,k)*compression_ratio

        vol_trans = vol_ice(i,j,k) ; vol_old = vol_ice(i,j,k+1)
        vol_ice(i,j,k+1) = vol_ice(i,j,k+1) + vol_trans
        Ivol_new = 1.0 / vol_ice(i,j,k+1)
      ! This is upwind advection across categories.  Improve it later.
        do m=1,G%NkIce
          T_ice(i,j,k+1,m) = (vol_trans*T_ice(i,j,k,m) + &
                              vol_old*T_ice(i,j,k+1,m)) * Ivol_new
        enddo
!        This is not quite right, or at least not consistent.
!        h_ice(i,j,k+1) = (vol_trans*hlim(k+1) + &
!                          vol_old*h_ice(i,j,k+1)) * Ivol_new
        h_ice(i,j,k+1) = vol_ice(i,j,k+1) / part_sz(i,j,k+1)

        if (vol_snow(i,j,k) > 0.0) then
          snow_trans = vol_snow(i,j,k) ; snow_old = vol_snow(i,j,k+1)
          vol_snow(i,j,k+1) = vol_snow(i,j,k+1) + vol_snow(i,j,k)
          Ivol_snow = 1.0 / vol_snow(i,j,k+1)

          t_snow(i,j,k+1) = (snow_trans*t_snow(i,j,k) + &
                             snow_old*t_snow(i,j,k+1)) * Ivol_snow
        endif
        h_snow(i,j,k+1) = vol_snow(i,j,k+1) / part_sz(i,j,k+1)

        vol_ice(i,j,k) = 0.0 ; vol_snow(i,j,k) = 0.0 ; part_sz(i,j,k) = 0.0
        h_ice(i,j,k) = 0.0 ; h_snow(i,j,k) = 0.0
      endif
    endif
  enddo ; enddo ; endif ; enddo

  k=G%CatIce
  do j=jsc,jec ; if (do_j(j)) then ; do i=isc,iec
    if (excess_cover(i,j) > 0.0) then
      if (part_sz(i,j,k) <= 1.0) &
        call SIS_error(FATAL, &
            "Category CatIce part_sz inconsistent with excess cover.")
      Icompress_here = part_sz(i,j,k) / (part_sz(i,j,k) - excess_cover(i,j))
      h_ice(i,j,k) = h_ice(i,j,k) * Icompress_here
      h_snow(i,j,k) = h_snow(i,j,k) * Icompress_here
      part_sz(i,j,k) = part_sz(i,j,k) - excess_cover(i,j)
      excess_cover(i,j) = 0.0
    endif
  enddo ; endif ; enddo

  if (CS%check_conservation) then
    ! Check for consistency between vol_ice, h_ice, and part_sz.
    do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
      if ((vol_ice(i,j,k) == 0.0) .and. (h_ice(i,j,k)*part_sz(i,j,k) /= 0.0)) then
        write(mesg,'("Compress mismatch at ",3(i8),": vol, h, part, hxp = zero, ",3(1pe15.6))') &
           i, j, k, h_ice(i,j,k), part_sz(i,j,k), h_ice(i,j,k)*part_sz(i,j,k)
        call SIS_error(WARNING, mesg, all_print=.true.)
      endif
      if (abs(vol_ice(i,j,k) - h_ice(i,j,k)*part_sz(i,j,k)) > 1e-12*vol_ice(i,j,k)) then
        write(mesg,'("Compress mismatch at ",3(i8),": vol, h, part, hxp = ",4(1pe15.6))') &
           i, j, k, vol_ice(i,j,k), h_ice(i,j,k), part_sz(i,j,k), h_ice(i,j,k)*part_sz(i,j,k)
        call SIS_error(WARNING, mesg, all_print=.true.)
      endif
    enddo ; enddo ; enddo
  endif

end subroutine compress_ice

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_advect - take adv_sub_steps upstream advection timesteps                 !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_advect(uc, vc, trc, dt_slow, G, CS, uf, vf)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: uc  ! x-face advecting velocity
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: vc  ! y-face advecting velocity
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: trc ! tracer to advect
  real,                              intent(in   ) :: dt_slow
  type(ice_transport_CS),            pointer       :: CS
  real, dimension(SZIB_(G),SZJ_(G)), optional, intent(  out) :: uf
  real, dimension(SZI_(G),SZJB_(G)), optional, intent(  out) :: vf
! Arguments: uc - The zonal ice velocity, in m s-1.
!  (in)      vc - The meridional ice velocity, in m s-1.
!  (inout)   trc - A tracer concentration times thickness, in m kg kg-1 or
!                  other units.
!  (in)      dt_slow - The amount of time over which the ice dynamics are to be
!                      advanced, in s.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.
!  (out)     vf - The averaged zonal tracer flux, in m3 kg kg-1 s-1.
!  (out)     vf - The averaged meridional tracer flux, in m3 kg kg-1 s-1.

  real, dimension(SZIB_(G),SZJ_(G)) :: uflx
  real, dimension(SZI_(G),SZJB_(G)) :: vflx
  real :: dt_adv, I_adv_steps
  integer :: l, i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (CS%adv_sub_steps==0) return;
  dt_adv = dt_slow/CS%adv_sub_steps

  if (present(uf)) uf(:,:) = 0.0
  if (present(vf)) vf(:,:) = 0.0

  uflx(:,:) = 0.0
  vflx(:,:) = 0.0

  do l=1,CS%adv_sub_steps
    do j=jsc,jec ; do I=isc-1,iec
      if( uc(I,j) > 0.0 ) then
         uflx(I,j) = uc(I,j) * trc(i,j) * G%dyCu(I,j)
      else
         uflx(I,j) = uc(I,j) * trc(i+1,j) * G%dyCu(I,j)
      endif
    enddo ; enddo

    do J=jsc-1,jec ; do i=isc,iec
      if( vc(i,J) > 0.0 ) then
         vflx(i,J) = vc(i,J) * trc(i,j) * G%dxCv(i,J)
      else
         vflx(i,J) = vc(i,J) * trc(i,j+1) * G%dxCv(i,J)
      endif
    enddo ; enddo

    do j=jsc,jec ; do i=isc,iec
      trc(i,j) = trc(i,j) + dt_adv * ( (uflx(I-1,j) - uflx(I,j)) + &
                 (vflx(i,J-1) - vflx(i,J)) ) * G%IareaT(i,j)
    enddo ; enddo

    call pass_var(trc, G%Domain)

    if (present(uf)) then ; do j=jsc,jec ; do I=isc-1,iec       
      uf(I,j) = uf(I,j) + uflx(I,j)
    enddo ; enddo ; endif

    if (present(vf)) then ;  do J=jsc-1,jec ; do i=isc,iec       
      vf(i,J) = vf(i,J) + vflx(i,J)
    enddo ; enddo ; endif
  enddo

  if (CS%adv_sub_steps>1) then
    I_adv_steps = 1.0/CS%adv_sub_steps
    if (present(uf)) then ; do j=jsc,jec ; do I=isc-1,iec       
      uf(I,j) = uf(I,j) * I_adv_steps
    enddo ; enddo ; endif
    if (present(vf)) then ;  do J=jsc-1,jec ; do i=isc,iec       
      vf(i,J) = vf(i,J) * I_adv_steps
    enddo ; enddo ; endif
  endif

end subroutine ice_advect

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine slab_ice_advect(uc, vc, trc, stop_lim, dt_slow, G, CS)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: uc  ! x-face advecting velocity
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: vc  ! y-face advecting velocity
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: trc ! tracer to advect
  real,                              intent(in   ) :: stop_lim
  real,                              intent(in   ) :: dt_slow
  type(ice_transport_CS),            pointer       :: CS
! Arguments: uc - The zonal ice velocity, in m s-1.
!  (in)      vc - The meridional ice velocity, in m s-1.
!  (inout)   trc - A tracer concentration times thickness, in m kg kg-1 or
!                  other units.
!  (in)      stop_lim - ?
!  (in)      dt_slow - The amount of time over which the ice dynamics are to be
!                      advanced, in s.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.

  real, dimension(SZIB_(G),SZJ_(G)) :: uflx
  real, dimension(SZI_(G),SZJB_(G)) :: vflx
  real :: avg, dif
  real :: dt_adv
  integer :: l, i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec


  if (CS%adv_sub_steps==0) return;
  dt_adv = dt_slow/CS%adv_sub_steps


  do l=1,CS%adv_sub_steps
    do j=jsc,jec ; do I=isc-1,iec
      avg = ( trc(i,j) + trc(i+1,j) )/2
      dif = trc(i+1,j) - trc(i,j)
      if( avg > stop_lim .and. uc(I,j) * dif > 0.0) then
        uflx(I,j) = 0.0
      else if( uc(i,j) > 0.0 ) then
        uflx(I,j) = uc(I,j) * trc(i,j) * G%dyCu(I,j)
      else
        uflx(I,j) = uc(I,j) * trc(i+1,j) * G%dyCu(I,j)
      endif
    enddo ; enddo

    do J=jsc-1,jec ; do i=isc,iec
      avg = ( trc(i,j) + trc(i,j+1) )/2
      dif = trc(i,j+1) - trc(i,j)
      if( avg > stop_lim .and. vc(i,J) * dif > 0.0) then
        vflx(i,J) = 0.0
      else if( vc(i,J) > 0.0 ) then
        vflx(i,J) = vc(i,J) * trc(i,j) * G%dxCv(i,J)
      else
        vflx(i,J) = vc(i,J) * trc(i,j+1) * G%dxCv(i,J)
      endif
    enddo ; enddo

    do j=jsc,jec ; do i=isc,iec
      trc(i,j) = trc(i,j) + dt_adv * ((uflx(I-1,j) - uflx(I,j)) + &
                                      (vflx(i,J-1) - vflx(i,J)) ) * G%IareaT(i,j)
    enddo ; enddo

    call pass_var(trc, G%Domain)
  enddo

end subroutine slab_ice_advect

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_redistribute - a simple ice redistribution scheme from Igor Polyakov     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_redistribute(cn, G, hs, tsn, hi, t1, t2, t3, t4, hlim)
!!$    real, intent(inout), dimension(1:G%CatIce)           :: cn, hs, hi, t1, t2
  type(sea_ice_grid_type), intent(inout) :: G
  real, intent(inout), dimension(:) :: cn, hs, tsn, hi, t1, t2, t3, t4
  real, dimension(:), intent(in) :: hlim

  real    :: cw                                 ! open water concentration
  integer :: k

  cw = 1-sum(cn)
  if (cw<0) cn(1) = cw + cn(1) ! open water has been eliminated by convergence

  do k=1,G%CatIce-1 ; if (cn(k)<0) then
    cn(k+1) = cn(k+1)+cn(k) ; cn(k) = 0 ! pass concentration deficit up to
    hs(k+1) = hs(k+1)+hs(k) ; hs(k) = 0 ! next thicker category
    tsn(k+1) = tsn(k+1)+tsn(k) ; tsn(k) = 0
    hi(k+1) = hi(k+1)+hi(k) ; hi(k) = 0
    t1(k+1) = t1(k+1)+t1(k) ; t1(k) = 0 ! NOTE:  here between the thm_pack and
    t2(k+1) = t2(k+1)+t2(k) ; t2(k) = 0 ! thm_unpack calls, all quantities are
    t3(k+1) = t3(k+1)+t3(k) ; t3(k) = 0 ! extensive, so we add instead of
    t4(k+1) = t4(k+1)+t4(k) ; t4(k) = 0 ! averaging
  endif ; enddo

  do k=1,G%CatIce-1 ; if (hi(k)>hlim(k+1)*cn(k)) then
    cn(k+1) = cn(k+1)+cn(k) ; cn(k) = 0 ! upper thickness limit exceeded
    hs(k+1) = hs(k+1)+hs(k) ; hs(k) = 0 ! move ice up to next thicker category
    tsn(k+1) = tsn(k+1)+tsn(k) ; tsn(k) = 0
    hi(k+1) = hi(k+1)+hi(k) ; hi(k) = 0
    t1(k+1) = t1(k+1)+t1(k) ; t1(k) = 0
    t2(k+1) = t2(k+1)+t2(k) ; t2(k) = 0
    t3(k+1) = t3(k+1)+t3(k) ; t3(k) = 0
    t4(k+1) = t4(k+1)+t4(k) ; t4(k) = 0
  endif ; enddo

  do k=G%CatIce,2,-1 ; if (hi(k)<hlim(k)*cn(k)) then
    cn(k-1) = cn(k-1)+cn(k) ; cn(k) = 0  ! lower thickness limit exceeded;
    hs(k-1) = hs(k-1)+hs(k) ; hs(k) = 0  ! move ice down to thinner category
    tsn(k-1) = tsn(k-1)+tsn(k) ; tsn(k) = 0
    hi(k-1) = hi(k-1)+hi(k) ; hi(k) = 0
    t1(k-1) = t1(k-1)+t1(k) ; t1(k) = 0
    t2(k-1) = t2(k-1)+t2(k) ; t2(k) = 0
    t3(k-1) = t3(k-1)+t3(k) ; t3(k) = 0
    t4(k-1) = t4(k-1)+t4(k) ; t4(k) = 0
  endif ; enddo

end subroutine ice_redistribute

subroutine get_total_amounts(vol_ice, vol_snow, G, tot_ice, tot_snow)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(in) :: vol_ice, vol_snow
  type(EFP_type), intent(out) :: tot_ice, tot_snow
! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      vol_ice - The volume per unit grid-cell area of the ice in each
!                      category in m.
!  (in)      vol_snow - The volume per unit grid-cell area of the snow atop the
!                       ice in each category in m.
!  (out)     tot_ice - The globally integrated total ice, in m3.
!  (out)     tot_snow - The globally integrated total snow, in m3.

  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: sum_vol_ice, sum_vol_snow
  real :: total
  integer :: i, j, k, m, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  sum_vol_ice(:,:) = 0.0
  sum_vol_snow(:,:) = 0.0
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_vol_ice(i,j) = sum_vol_ice(i,j) + G%areaT(i,j) * vol_ice(i,j,k)
    sum_vol_snow(i,j) = sum_vol_snow(i,j) + G%areaT(i,j) * vol_snow(i,j,k)
  enddo ; enddo ; enddo

  total = reproducing_sum(sum_vol_ice, EFP_sum=tot_ice)
  total = reproducing_sum(sum_vol_snow, EFP_sum=tot_snow)

end subroutine get_total_amounts

subroutine get_total_enthalpy(h_ice, h_snow, part_sz, t_ice, t_snow, &
                              pocket_enth, G, enth_ice, enth_snow)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(in) :: h_ice, h_snow, t_snow
  real, dimension(SZI_(G),SZJ_(G),SZCAT0_(G)), intent(in) :: part_sz
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G),SZK_ICE_(G)),  intent(in) :: t_ice
  real, dimension(SZK_ICE_(G)),  intent(in) :: pocket_enth
  type(EFP_type), intent(out) :: enth_ice, enth_snow
! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      vol_ice - The volume per unit grid-cell area of the ice in each
!                      category in m.
!  (in)      vol_snow - The volume per unit grid-cell area of the snow atop the
!                       ice in each category in m.
!  (out)     tot_ice - The globally integrated total ice, in m3.
!  (out)     tot_snow - The globally integrated total snow, in m3.

  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: sum_enth_ice, sum_enth_snow
  real :: total, LI_CI, enth_chg
  integer :: i, j, k, m, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call get_thermo_coefs(max_enthalpy_chg=LI_CI)

  sum_enth_ice(:,:) = 0.0
  sum_enth_snow(:,:) = 0.0
  do m=1,G%NkIce ; do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    ! Determine the enthalpy correction due to brine pockets.
    if ((t_ice(i,j,k,m) >= 0.0) .or. (pocket_enth(m) > LI_CI*(-t_ice(i,j,k,m)))) then
      enth_chg = LI_CI
    else
      enth_chg = -pocket_enth(m)/t_ice(i,j,k,m)  ! The enthalpy correction due to brine pockets.
    endif
    
    sum_enth_ice(i,j) = sum_enth_ice(i,j) + (G%areaT(i,j) * (h_ice(i,j,k)*part_sz(i,j,k))) * &
                           (t_ice(i,j,k,m) + enth_chg) 
  enddo ; enddo ; enddo ; enddo
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_enth_snow(i,j) = sum_enth_snow(i,j) + (G%areaT(i,j) * (h_snow(i,j,k)*part_sz(i,j,k))) * &
                        T_snow(i,j,k)
  enddo ; enddo ; enddo

  total = reproducing_sum(sum_enth_ice, EFP_sum=enth_ice)
  total = reproducing_sum(sum_enth_snow, EFP_sum=enth_snow)

end subroutine get_total_enthalpy

subroutine get_snow_heat(vol_snow, t_snow, G, enth_snow)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),  intent(in) :: vol_snow, t_snow
  type(EFP_type), intent(out) ::  enth_snow
! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      vol_ice - The volume per unit grid-cell area of the ice in each
!                      category in m.
!  (in)      vol_snow - The volume per unit grid-cell area of the snow atop the
!                       ice in each category in m.
!  (out)     tot_ice - The globally integrated total ice, in m3.
!  (out)     tot_snow - The globally integrated total snow, in m3.

  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: sum_enth_snow
  real :: total
  integer :: i, j, k, m, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  sum_enth_snow(:,:) = 0.0
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_enth_snow(i,j) = sum_enth_snow(i,j) + (G%areaT(i,j) * vol_snow(i,j,k)) * &
                        T_snow(i,j,k)
  enddo ; enddo ; enddo

  total = reproducing_sum(sum_enth_snow, EFP_sum=enth_snow)

end subroutine get_snow_heat

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_transport_init - initialize the ice transport and set parameters.        !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_transport_init(Time, G, param_file, diag, CS)
  type(time_type),     target, intent(in)    :: Time
  type(sea_ice_grid_type),     intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(ice_transport_CS),      pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.

!   This subroutine sets the parameters and registers the diagnostics associated
! with the ice dynamics.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "ice_transport" ! This module's name.
  real, parameter :: missing = -1e34

  if (associated(CS)) then
    call SIS_error(WARNING, "ice_transport_init called with an associated control structure.")
!    call SIS_error(FATAL, "ice_transport_init called with an unassociated control structure. \n"//&
!                    "ice_transport_register_restarts must be called before ice_transport_init.")
    return
  endif
  allocate(CS)

  CS%diag => diag
  CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "SPECIFIED_ICE", CS%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  if ( CS%specified_ice ) then
    CS%adv_sub_steps = 0
    call log_param(param_file, mod, "NSTEPS_ADV", CS%adv_sub_steps, &
                 "The number of advective iterations for each slow time \n"//&
                 "step.  With SPECIFIED_ICE this is always 0.")
  else
    call get_param(param_file, mod, "NSTEPS_ADV", CS%adv_sub_steps, &
                 "The number of advective iterations for each slow time \n"//&
                 "step.", default=1)
  endif

  call get_param(param_file, mod, "ICE_CHANNEL_VISCOSITY", CS%chan_visc, &
                 "A viscosity used in one-cell wide channels to \n"//&
                 "parameterize transport, especially with B-grid sea ice \n"//&
                 "coupled to a C-grid ocean model.", units="m2 s-1", default=0.0)
  call get_param(param_file, mod, "ICE_CHANNEL_SMAG_COEF", CS%smag_ocn, &
                 "A Smagorinsky coefficient for viscosity in channels.", &
                 units="Nondim", default=0.15)
  call get_param(param_file, mod, "ICE_CHANNEL_CFL_LIMIT", CS%chan_cfl_limit, &
                 "The CFL limit that is applied to the parameterized \n"//&
                 "viscous transport in single-point channels.", &
                 units="Nondim", default=0.25)

  call get_param(param_file, mod, "RHO_OCEAN", CS%Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0)
  call get_param(param_file, mod, "RHO_ICE", CS%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0)
  call get_param(param_file, mod, "RHO_SNOW", CS%Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0)

  call get_param(param_file, mod, "USE_SLAB_ICE", CS%SLAB_ICE, &
                 "If true, use the very old slab-style ice.", default=.false.)
  call get_param(param_file, mod, "SIS1_ICE_TRANSPORT", CS%SIS1_transport, &
                 "If true, use SIS1 code to solve the ice continuity \n"//&
                 "equation and transport tracers.", default=.true.)
  call get_param(param_file, mod, "CHECK_ICE_TRANSPORT_CONSERVATION", CS%check_conservation, &
                 "If true, use add multiple diagnostics of ice and snow \n"//&
                 "mass conservation in the sea-ice transport code.  This \n"//&
                 "is expensive and should be used sparingly.", default=.false.)

  CS%id_ustar = register_diag_field('ice_model', 'U_STAR', diag%axesCu1, Time, &
              'channel transport velocity - x component', 'm/s', missing_value=missing)
  CS%id_vstar = register_diag_field('ice_model', 'V_STAR', diag%axesCv1, Time, &
              'channel transport velocity - y component', 'm/s', missing_value=missing)
  CS%id_uocean = register_diag_field('ice_model', 'U_CHAN_OCN', diag%axesCu1, Time, &
              'ocean component of channel transport - x', 'm/s', missing_value=missing)
  CS%id_vocean = register_diag_field('ice_model', 'V_CHAN_OCN', diag%axesCv1, Time, &
              'ocean component of channel transport - y', 'm/s', missing_value=missing)
  CS%id_uchan = register_diag_field('ice_model', 'U_CHAN_VISC', diag%axesCu1, Time, &
              'viscous component of channel transport - x', 'm/s', missing_value=missing)
  CS%id_vchan = register_diag_field('ice_model', 'V_CHAN_VISC', diag%axesCv1, Time, &
              'viscous component of channel transport - y', 'm/s', missing_value=missing)
  CS%id_ix_trans = register_diag_field('ice_model', 'IX_TRANS', diag%axesCu1, Time, &
               'x-direction ice transport', 'kg/s', missing_value=missing)
  CS%id_iy_trans = register_diag_field('ice_model', 'IY_TRANS', diag%axesCv1, Time, &
               'y-direction ice transport', 'kg/s', missing_value=missing)

end subroutine ice_transport_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_transport_end - deallocate the memory associated with this module.       !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_transport_end(CS)
  type(ice_transport_CS), pointer :: CS

  deallocate(CS)
end subroutine ice_transport_end

end module ice_transport_mod
