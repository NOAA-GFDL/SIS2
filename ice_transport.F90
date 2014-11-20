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
use SIS_tracer_registry, only : SIS_tracer_registry_type, get_SIS_tracer_pointer
use SIS_tracer_registry, only : update_SIS_tracer_halos, set_massless_SIS_tracers
use SIS_tracer_advect, only : advect_tracers_thicker, SIS_tracer_advect_CS
use SIS_tracer_advect, only : advect_SIS_tracers, SIS_tracer_advect_init, SIS_tracer_advect_end
use SIS_tracer_advect, only : advect_scalar
use SIS_continuity, only :  SIS_continuity_init, SIS_continuity_end
use SIS_continuity, only :  continuity=>ice_continuity, SIS_continuity_CS

use ice_grid_mod, only : sea_ice_grid_type
use ice_ridging_mod, only : ice_ridging

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
  real :: Roll_factor         ! A factor by which the propensity of small
                              ! amounts of thick sea-ice to become thinner by
                              ! rolling is increased, or 0 to disable rolling.
                              ! Sensible values are 0 or larger than 1.
  logical :: specified_ice    ! If true, the sea ice is specified and there is
                              ! no need for ice dynamics.
  logical :: check_conservation ! If true, write out verbose diagnostics of conservation.
  integer :: adv_sub_steps    ! The number of advective iterations for each slow
                              ! time step.
  type(time_type), pointer :: Time ! A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  logical :: do_ridging
  logical :: use_SIS_continuity ! Use the more modern continuity solver from SIS.
  logical :: use_SIS_thickness_advection ! Use the SIS tracer advection schemes for thickness.
  type(SIS_continuity_CS),    pointer :: continuity_CSp => NULL()
  type(SIS_tracer_advect_CS), pointer :: SIS_tr_adv_CSp => NULL()
  integer :: id_ustar = -1, id_uocean = -1, id_uchan = -1
  integer :: id_vstar = -1, id_vocean = -1, id_vchan = -1
  integer :: id_ix_trans = -1, id_iy_trans = -1
end type ice_transport_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! transport - do ice transport and thickness class redistribution              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_transport(part_sz, mH_ice, mH_snow, uc, vc, TrReg, &
                         sea_lev, mH_lim, dt_slow, G, CS, &
                         rdg_hice, age_ice, snow2ocn, rdg_rate, rdg_open, &
                         rdg_vosh)
  type(sea_ice_grid_type),                      intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(G)), intent(inout) :: part_sz
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),   intent(inout) :: mH_ice, mH_snow
  type(SIS_tracer_registry_type),               pointer       :: TrReg
  real, dimension(SZIB_(G),SZJ_(G)),            intent(inout) :: uc
  real, dimension(SZI_(G),SZJB_(G)),            intent(inout) :: vc
  real, dimension(SZI_(G),SZJ_(G)),             intent(in)    :: sea_lev
  real, dimension(:),                           intent(in)    :: mH_lim  ! Move to grid type?
  real,                                         intent(in)    :: dt_slow
  type(ice_transport_CS),                       pointer       :: CS
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),   intent(inout) :: rdg_hice
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),   intent(inout) :: age_ice
  real, dimension(SZI_(G),SZJ_(G)),             intent(inout) :: snow2ocn ! snow volume [m] dumped into ocean during ridging
  real, dimension(SZI_(G),SZJ_(G)),             intent(inout) :: rdg_rate
  real, dimension(SZI_(G),SZJ_(G)),             intent(inout) :: rdg_open ! formation rate of open water due to ridging
  real, dimension(SZI_(G),SZJ_(G)),             intent(inout) :: rdg_vosh ! rate of ice volume shifted from level to ridged ice
! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1, in/out.
!  (inout)   mH_ice - The mass per unit area of the ice in each category in H (often kg m-2).
!  (inout)   mH_snow - The mass per unit area of the snow atop the ice in each
!                     category in H (often kg m-2).
!  (in)      uc - The zonal ice velocity, in m s-1.
!  (in)      vc - The meridional ice velocity, in m s-1.
!  (inout)   TrReg - The registry of registered SIS ice and snow tracers.
!  (in)      mH_lim - The lower ice-loading limit of each category, in H (often kg m-2).
!  (in)      sea_lev - The height of the sea level, including contributions
!                      from non-levitating ice from an earlier time step, in m.
!  (in)      dt_slow - The amount of time over which the ice dynamics are to be
!                      advanced, in s.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.
  real, dimension(:,:,:,:), pointer :: &
    heat_ice=>NULL(), & ! Pointers to the enth_ice and enth_snow arrays from the
    heat_snow=>NULL()   ! SIS tracer registry.  enth_ice is the enthalpy of the
                        ! ice in each category and layer, while enth_snow is the
                        ! enthalpy of the snow atop the ice in each category.
                        ! Both are in enth_units (J or rescaled).
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(G)) :: &
    uh_ice, &  ! Zonal fluxes in H m2 s-1.
    uh_snow    ! Zonal fluxes in H m2 s-1.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    uf, & ! Zonal fluxes in m3 s-1 and kg s-1.
    ustar, ustaro, ustarv ! Local variables, transporting velocities
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(G)) :: &
    vh_ice, &  ! Meridional fluxes in m3 s-1 and kg s-1.
    vh_snow    ! Meridional fluxes in m3 s-1 and kg s-1.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    vf, & ! Meridional fluxes in m3 s-1 and kg s-1.
    vstar, vstaro, vstarv ! Local variables, transporting velocities
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)) :: &
    mca_ice, mca_snow, &  ! The mass of snow and ice per unit total area in a
                          ! cell, in units of H (often kg m-2).  "mca" stands
                          ! for "mass cell averaged"
    mca0_ice, mca0_snow   ! The initial mass of snow and ice per unit total
                          ! area in a cell, in units of H (often kg m-2).
  real :: h_in_m          ! The ice thickness in m.
  real :: hca_in_m        ! The ice thickness averaged over the whole cell in m.
  real, dimension(SZI_(G),SZJ_(G)) :: opnwtr
  real, dimension(SZI_(G),SZJ_(G)) :: ice_cover ! The summed fractional ice concentration, ND.
  real :: u_visc, u_ocn, cnn, grad_eta ! Variables for channel parameterization

  real :: I_mca_ice

  type(EFP_type) :: tot_ice(2), tot_snow(2), enth_ice(2), enth_snow(2)
  real :: I_tot_ice, I_tot_snow
  real :: C1_3 = 1.0/3.0

  real :: dt_adv
  character(len=200) :: mesg
  integer :: i, j, k, m, bad, isc, iec, jsc, jec, isd, ied, jsd, jed, nL
  integer :: iTransportSubcycles ! For transport sub-cycling
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (CS%slab_ice) then
    call pass_vector(uc, vc, G%Domain, stagger=CGRID_NE)
    call slab_ice_advect(uc, vc, mH_ice(:,:,1), 4.0*G%kg_m2_to_H, dt_slow, G, CS)
    call pass_var(mH_ice(:,:,2), G%Domain)
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      if (mH_ice(i,j,1) > 0.0) then
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
    ustar(:,:) = 0. ; vstar(:,:) = 0.
    ustaro(:,:) = 0. ; vstaro(:,:) = 0.
    ustarv(:,:) = 0. ; vstarv(:,:) = 0.
    do j=jsc,jec ; do I=isc-1,iec
      if ((uc(I,j)==0.) .and. & ! this is a redundant test due to following line
          (G%mask2dBu(I,J)+G%mask2dBu(I,J-1)==0.) .and. &  ! =0 => no vels
          (G%mask2dT(i,j)*G%mask2dT(i+1,j)>0.)) then ! >0 => open for transport
        grad_eta=(sea_lev(i+1,j)-sea_lev(i,j)) * G%IdxCu(I,j) ! delta_i eta/dx
        u_visc=-G%g_Earth*((G%dy_Cu(I,j)*G%dy_Cu(I,j))/(12.*CS%chan_visc)) & ! -g*dy^2/(12*visc)
                 *grad_eta                                                  ! d/dx eta
        u_ocn=sqrt( G%g_Earth*G%dy_Cu(I,j)*abs(grad_eta)/(36.*CS%smag_ocn) ) ! Magnitude of ocean current
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
        u_visc=-G%g_Earth*((G%dx_Cv(i,J)*G%dx_Cv(i,J))/(12.*CS%chan_visc)) & ! -g*dy^2/(12*visc)
                *grad_eta                                                  ! d/dx eta
        u_ocn=sqrt( G%g_Earth*G%dx_Cv(i,J)*abs(grad_eta)/(36.*CS%smag_ocn) ) ! Magnitude of ocean current
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

  if (CS%check_conservation) then
    call get_total_amounts(mH_ice, mH_snow, part_sz, G, tot_ice(1), tot_snow(1))
    call get_total_enthalpy(mH_ice, mH_snow, part_sz, TrReg, G, enth_ice(1), &
                            enth_snow(1))
  endif

  ! Make sure that ice is in the right thickness category before advection.
  call adjust_ice_categories(mH_lim, mH_ice, mH_snow, part_sz, &
                             TrReg, G, CS) !Niki: add ridging?

  !   Determine the whole-cell averaged mass of snow and ice.
  mca_ice(:,:,:) = 0.0 ; mca_snow(:,:,:) = 0.0
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    if (mH_ice(i,j,k)>0.0) then
      mca_ice(i,j,k) = part_sz(i,j,k)*mH_ice(i,j,k)
      mca_snow(i,j,k) = part_sz(i,j,k)*mH_snow(i,j,k)
    else
      if (part_sz(i,j,k)*mH_snow(i,j,k) > 0.0) then
        call SIS_error(FATAL, "Input to ice_transport, non-zero snow mass rests atop no ice.")
      endif
      part_sz(i,j,k) = 0.0 ; mca_ice(i,j,k) = 0.0
      mca_snow(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo

  ! Do the transport via the continuity equations and tracer conservation
  ! equations for mH_ice and tracers, inverting for the fractional size of
  ! each partition.
  call pass_var(part_sz, G%Domain) ! cannot be combined with updates below
  call update_SIS_tracer_halos(TrReg, G, complete=.false.)
  call pass_var(mca_ice,  G%Domain, complete=.false.)
  call pass_var(mca_snow, G%Domain, complete=.false.)
  call pass_var(mH_ice, G%Domain, complete=.true.)


  dt_adv = dt_slow / real(CS%adv_sub_steps)
  do iTransportSubcycles = 1, CS%adv_sub_steps
    if (iTransportSubcycles>1) then ! Do not need to update on first iteration
      call update_SIS_tracer_halos(TrReg, G, complete=.false.)
      call pass_var(mca_ice,  G%Domain, complete=.false.)
      call pass_var(mca_snow, G%Domain, complete=.false.)
      call pass_var(mH_ice, G%Domain, complete=.true.)
    endif

    do k=1,G%CatIce ; do j=jsd,jed ; do i=isd,ied
      mca0_ice(i,j,k) = mca_ice(i,j,k)
      mca0_snow(i,j,k) = mca_snow(i,j,k)
    enddo ; enddo ; enddo
    if (CS%use_SIS_continuity) then
      call continuity(uc, vc, mca0_ice, mca_ice, uh_ice, vh_ice, dt_adv, G, CS%continuity_CSp)
      call continuity(uc, vc, mca0_snow, mca_snow, uh_snow, vh_snow, dt_adv, G, CS%continuity_CSp)
    else
      call ice_continuity(uc, vc, mca0_ice, mca_ice, uh_ice, vh_ice, dt_adv, G, CS)
      call ice_continuity(uc, vc, mca0_snow, mca_snow, uh_snow, vh_snow, dt_adv, G, CS)
    endif

    if (CS%use_SIS_thickness_advection) then
      call advect_scalar(mH_ice, mca0_ice, mca_ice, uh_ice, vh_ice, dt_adv, G, CS%SIS_tr_adv_CSp)
    else
      call advect_ice_tracer(mca0_ice, mca_ice, uh_ice, vh_ice, mH_ice, dt_adv, G, CS)
    endif

    call advect_SIS_tracers(mca0_ice, mca_ice, uh_ice, vh_ice, dt_adv, G, &
                            CS%SIS_tr_adv_CSp, TrReg, snow_tr=.false.)
    call advect_SIS_tracers(mca0_snow, mca_snow, uh_snow, vh_snow, dt_adv, G, &
                            CS%SIS_tr_adv_CSp, TrReg, snow_tr=.true.)
  enddo ! iTransportSubcycles

  ! Add code to make sure that mH_ice(i,j,1) > mH_lim(1).
  do j=jsc,jec ; do i=isc,iec
    if ((mca_ice(i,j,1) > 0.0) .and. (mH_ice(i,j,1) < mH_lim(1))) then
      mH_ice(i,j,1) = mH_lim(1)
    endif
  enddo ; enddo

  ! Convert mca_ice and mca_snow back to part_sz and mH_snow.
  ice_cover(:,:) = 0.0
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    if (mca_ice(i,j,k) > 0.0) then
      if (CS%roll_factor * (mH_ice(i,j,k)*G%H_to_kg_m2/CS%Rho_Ice)**3 > &
          (mca_ice(i,j,k)*G%H_to_kg_m2/CS%Rho_Ice)*G%areaT(i,j)) then
        ! This ice is thicker than it is wide even if all the ice in a grid
        ! cell is collected into a single cube, so it will roll.  Any snow on
        ! top will simply be redistributed into a thinner layer, although it
        ! should probably be dumped into the ocean.  Rolling makes the ice
        ! thinner so that it melts faster, but it should never be made thinner
        ! than mH_lim(1).
        mH_ice(i,j,k) = max((CS%Rho_ice*G%kg_m2_to_H) * &
             sqrt((mca_ice(i,j,k)*G%areaT(i,j)) / &
                  (CS%roll_factor * mH_ice(i,j,k)) ), mH_lim(1))
      endif

      part_sz(i,j,k) = mca_ice(i,j,k) / mH_ice(i,j,k)
      mH_snow(i,j,k) = mH_ice(i,j,k) * (mca_snow(i,j,k) / mca_ice(i,j,k))
      ice_cover(i,j) = ice_cover(i,j) + part_sz(i,j,k)
    else
      part_sz(i,j,k) = 0.0 ; mH_ice(i,j,k) = 0.0
      if (mca_snow(i,j,k) > 0.0) &
        call SIS_error(FATAL, &
          "Positive mca_snow values should not exist without ice.")
      mH_snow(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo
  do j=jsc,jec ; do i=isc,iec
    part_sz(i,j,0) = 1.0-ice_cover(i,j)
  enddo ; enddo

  ! Compress the ice where the fractional coverage exceeds 1, starting with
  ! the thinnest categories.  This is a minimalist version of a sea-ice
  ! ridging scheme.  A more complete ridging scheme would also compress
  ! thicker ice and allow the fractional ice coverage to drop below 1.
  call compress_ice(part_sz, mH_lim, mca_ice, mca_snow, mH_ice, mH_snow, &
                    TrReg, G, CS)

  !   Handle massless categories.
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    if (mca_ice(i,j,k)<=0.0) then
      part_sz(i,j,k) = 0.0 ; mH_ice(i,j,k) = 0.0
      mH_snow(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo
  call set_massless_SIS_tracers(mca_ice, TrReg, G, compute_domain=.true.)

    ! Is sum(part_sz) = 1 ?

!  Niki: TOM does the ridging after redistribute which would need IST%age_ice and IST%rdg_hice below.
!   !  ### THIS IS HARD-CODED ONLY TO WORK WITH 2 LAYERS.
!   !  ### heat_snow AND OTHER TRACERS ARE OMITTED.
!   if (CS%do_ridging) then
!     do j=jsc,jec ; do i=isc,iec
!       snow2ocn(i,j) = 0.0 !TOM> initializing snow2ocean
!       if (sum(mH_ice(i,j,:)) > 1.e-10*CS%Rho_ice .and. &
!           sum(part_sz(i,j,1:G%CatIce)) > 0.01) &
!         call ice_ridging(G%CatIce, part_sz(i,j,:), mH_ice(i,j,:), &
!             mH_snow(i,j,:), &
!             heat_ice(i,j,:,1), heat_ice(i,j,:,2), & !Niki: Is this correct? Bob: No, 2-layers hard-coded.
!             age_ice(i,j,:), snow2ocn(i,j), rdg_rate(i,j), rdg_hice(i,j,:), &
!             dt_slow, mH_lim, rdg_open(i,j), rdg_vosh(i,j))
!     enddo ; enddo
!   endif   ! do_ridging

  uf(:,:) = 0.0; vf(:,:) = 0.0
  if ((CS%id_ix_trans>0) .or. (CS%id_iy_trans>0)) then ; do k=1,G%CatIce
    do j=jsc,jec ; do I=isc-1,iec
      uf(I,j) = uf(I,j) + (uh_snow(I,j,k) + uh_ice(I,j,k))
    enddo ; enddo
    do J=jsc-1,jec ; do i=isc,iec
      vf(i,J) = vf(i,J) + (vh_snow(i,J,k) + vh_ice(i,J,k))
    enddo ; enddo
  enddo ; endif

  call pass_var(part_sz, G%Domain) ! cannot be combined with the two updates below
  call pass_var(mH_snow, G%Domain, complete=.false.)
  call pass_var(mH_ice, G%Domain, complete=.true.)

  !   Recalculate part_sz(:,:,0) to ensure that the sum of part_sz adds up to 1.
  ! Compress_ice should already have taken care of this within the computational
  ! domain, but with a slightly different order of arithmetic.  The max is here
  ! to avoid tiny negative values of order -1e-16 from round-off in the
  ! difference between ice_cover and 1.
  ice_cover(:,:) = 0.0
  do k=1,G%CatIce ; ice_cover(:,:) = ice_cover(:,:) + part_sz(:,:,k) ; enddo
  part_sz(:,:,0) = max(1.0 - ice_cover(:,:), 0.0)

  if (CS%check_conservation) then
    call get_total_amounts(mH_ice, mH_snow, part_sz, G, tot_ice(2), tot_snow(2))

    call get_total_enthalpy(mH_ice, mH_snow, part_sz, TrReg, G, enth_ice(2), &
                            enth_snow(2))

    if (is_root_pe()) then
      I_tot_ice  = abs(EFP_to_real(tot_ice(1)))
      if (I_tot_ice > 0.0) I_tot_ice = 1.0 / I_tot_ice    ! Adcroft's rule inverse.
      I_tot_snow = abs(EFP_to_real(tot_snow(1)))
      if (I_tot_snow > 0.0) I_tot_snow = 1.0 / I_tot_snow ! Adcroft's rule inverse.
      write(*,'("  Total Ice mass:  ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(tot_ice(2)), EFP_real_diff(tot_ice(2),tot_ice(1)), &
        EFP_real_diff(tot_ice(2),tot_ice(1)) * I_tot_ice
      write(*,'("  Total Snow mass: ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
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
!  (in)      h_in - Initial layer thickness, in H (kg m-2).
!  (out)     h - Final layer thickness, in H (kg m-2).
!  (out)     uh - Mass flux through zonal faces = u*h*dy, H m2 s-1.
!  (out)     vh - Mass flux through meridional faces = v*h*dx, in H m2 s-1.
!  (in)      dt - Time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 ice_transport_init.
  real h_up
  integer :: i, j, k, is, ie, js, je
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do k=1,G%CatIce ; do j=js,je ; do I=is-1,ie
    if (u(I,j) >= 0.0) then ; h_up = h_in(i,j,k)
    else ; h_up = h_in(i+1,j,k) ; endif
    uh(I,j,k) = G%dy_Cu(I,j) * u(I,j) * h_up
  enddo ; enddo ; enddo

  do k=1,G%CatIce ; do J=js-1,je ; do i=is,ie
    if (v(i,J) >= 0.0) then ; h_up = h_in(i,j,k)
    else ; h_up = h_in(i,j+1,k) ; endif
    vh(i,J,k) = G%dx_Cv(i,J) * v(i,J) * h_up
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

! Arguments: h_prev - Category thickness times fractional coverage before advection, in m or kg m-2.
!  (in)      h_end - Layer thickness times fractional coverage after advection, in m or kg m-2.
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
  real    :: vol_end, Ivol_end
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
      flux_x(I,j,k) = (dt*uhtr(I,j,k)) * tr_up
    enddo ; enddo

    do J=js-1,je ; do i=is,ie
      if (vhtr(i,J,k) >= 0.0) then ; tr_up = tr(i,j,k)
      else ; tr_up = tr(i,j+1,k) ; endif
      flux_y(i,J,k) = (dt*vhtr(i,J,k)) * tr_up
    enddo ; enddo

    do j=js,je ; do i=is,ie
      vol_end = (G%areaT(i,j) * h_end(i,j,k))
      Ivol_end = 0.0 ; if (vol_end > 0.0) Ivol_end = 1.0 / vol_end
      tr(i,j,k) = ( (G%areaT(i,j)*h_prev(i,j,k))*tr(i,j,k) - &
                   ((flux_x(I,j,k) - flux_x(I-1,j,k)) + &
                    (flux_y(i,J,k) - flux_y(i,J-1,k))) ) * Ivol_end
    enddo ; enddo
  enddo

end subroutine advect_ice_tracer

subroutine adjust_ice_categories(mH_lim, mH_ice, mH_snow, part_sz, &
                                 TrReg, G, CS)
  type(sea_ice_grid_type),                    intent(inout) :: G
  real, dimension(:),                         intent(in)    :: mH_lim  ! Move to grid type?
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)), intent(inout) :: mH_ice, mH_snow
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(G)), intent(inout) :: part_sz
  type(SIS_tracer_registry_type),             pointer       :: TrReg
  type(ice_transport_CS),                     pointer       :: CS

!   This subroutine moves mass between thickness categories if it is thinner or
! thicker than the bounding limits of each category.

! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      mH_lim - The lower thickness limit of each category, in H (often kg m-2).
!  (inout)   mca_ice - The mass per unit grid-cell area of the ice in each
!                      category in H (often kg m-2).
!  (inout)   mca_snow - The mass per unit grid-cell area of the snow atop the
!                       ice in each category in H (often kg m-2).
!  (inout)   mH_ice - The thickness of the ice in each category in H (often kg m-2).
!  (inout)   TrReg - The registry of registered SIS ice and snow tracers.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.
  real :: mca_trans  ! The cell-averaged ice mass transfered between categories, in kg m-2.
  real :: part_trans ! The fractional area transfered between categories, nondim.
  real :: snow_trans ! The cell-averaged snow transfered between categories, in kg m-2.
  real :: I_mH_lim1  ! The inverse of the lower thickness limit, in m2 kg-1.
  real, dimension(SZI_(G),SZCAT_(G)) :: &
    mca_ice, mca_snow, &  ! The mass of snow and ice per unit total area in a
                          ! cell, in units of H (often kg m-2).  "mca" stands
                          ! for "mass cell averaged"
    mca0_ice, mca0_snow, & ! Initial ice and snow masses per unit cell area, in kg m-2.
    trans_ice, trans_snow  ! Cross-catagory transfers of ice and snow mass, in kg m-2.
  logical :: do_any, do_j(SZJ_(G))
  integer :: i, j, k, m, is, ie, js, je, nLay
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
 
  I_mH_lim1 = 1.0 / mH_lim(1)

  ! Zero out the part_size of any massless categories.
  do k=1,G%CatIce ; do j=js,je ; do i=is,ie ; if (mH_ice(i,j,k) <= 0.0) then
    if (mH_ice(i,j,k) < 0.0) then
      call SIS_error(FATAL, "Input to adjust_ice_categories, negative ice mass.")
    endif
    if (mH_snow(i,j,k) > 0.0) then
      call SIS_error(FATAL, "Input to adjust_ice_categories, non-zero snow mass rests atop no ice.")
    endif
    part_sz(i,j,k) = 0.0
  endif ; enddo ; enddo ; enddo

  ! Only work on rows that have any sea ice at all.
  do j=js,je
    do_j(j) = .false.
    do k=1,G%CatIce
      do i=is,ie ; if (part_sz(i,j,k)*mH_ice(i,j,k) > 0.0) then
        do_j(j) = .true. ; exit
      endif ; enddo
      if (do_j(j)) exit
    enddo
  enddo

  do j=js,je ; if (do_j(j)) then
    do k=1,G%CatIce ; do i=is,ie
      mca_ice(i,k) = part_sz(i,j,k)*mH_ice(i,j,k)
      mca_snow(i,k) = part_sz(i,j,k)*mH_snow(i,j,k)
      mca0_ice(i,k) = mca_ice(i,k) ; mca0_snow(i,k) = mca_snow(i,k)
    enddo ; enddo
    trans_ice(:,:) = 0.0 ; trans_snow(:,:) = 0.0
    do_any = .false.

    do k=1,G%CatIce-1 ; do i=is,ie
      if ((mca_ice(i,k) > 0.0) .and. (mH_ice(i,j,k) > mH_lim(k+1))) then
        ! Move some or all of the ice to a thicker category.
        ! For now move all of it.
        mca_trans = mca_ice(i,k)
        part_trans = part_sz(i,j,k) ! * (mca_trans / mca_ice) * (mH_ice / h_trans)
        snow_trans = mca_snow(i,k) ! * (part_trans / part_sz) = 1

        trans_ice(i,K) = mca_trans
        trans_snow(i,K) = snow_trans
        do_any = .true.

        ! Use area-weighted remapped thicknesses so that the total ice area and
        ! mass are both conserved in the remapping operation.  Using a mass-
        ! weighted average thickness instead would cause the ice to systematically
        ! contract in area.
        mH_ice(i,j,k+1) = (part_trans*mH_ice(i,j,k) + &
                           part_sz(i,j,k+1)*mH_ice(i,j,k+1)) / &
                          (part_trans + part_sz(i,j,k+1))
        ! h should be the first thing to correct via a non-constant profile, and
        ! can be improved independent of T & S.
        mH_ice(i,j,k) = mH_lim(k+1)
        part_sz(i,j,k+1) = part_sz(i,j,k+1) + part_trans
        part_sz(i,j,k) = part_sz(i,j,k) - part_trans

        mca_ice(i,k+1) = mca_ice(i,k+1) + mca_trans
        mca_ice(i,k) = mca_ice(i,k) - mca_trans

        mca_snow(i,k+1) = mca_snow(i,k+1) + snow_trans
        mca_snow(i,k) = mca_snow(i,k) - snow_trans
        
        mH_snow(i,j,k) = 0.0 ; mH_snow(i,j,k+1) = 0.0
        if (part_sz(i,j,k)>0.0) mH_snow(i,j,k) = mca_snow(i,k) / part_sz(i,j,k)
        if (part_sz(i,j,k+1)>0.0) mH_snow(i,j,k+1) = mca_snow(i,k+1) / part_sz(i,j,k+1)
      endif
    enddo ; enddo

    if (do_any) then
      call advect_tracers_thicker(mca0_ice, trans_ice, G, CS%SIS_tr_adv_CSp, &
                                  TrReg, .false., j, is, ie)
      call advect_tracers_thicker(mca0_snow, trans_snow, G, CS%SIS_tr_adv_CSp, &
                                  TrReg, .true., j, is, ie)
    endif

    do k=1,G%CatIce ; do i=is,ie
      mca0_ice(i,k) = mca_ice(i,k) ; mca0_snow(i,k) = mca_snow(i,k)
    enddo ; enddo
    trans_ice(:,:) = 0.0 ; trans_snow(:,:) = 0.0
    do_any = .false.

    do k=G%CatIce,2,-1 ; do i=is,ie
      if ((mca_ice(i,k) > 0.0) .and. (mH_ice(i,j,k) < mH_lim(k))) then
        ! Move some or all of the ice to a thinner category.
        ! For now move all of it.
        mca_trans = mca_ice(i,k)
        part_trans = part_sz(i,j,k) ! * (mca_trans / mca_ice) * (mH_ice / h_trans)
        snow_trans = mca_snow(i,k) ! * (part_trans / part_sz) = 1

        do_any = .true.
        trans_ice(i,K-1) = -mca_trans  ! Note the shifted index conventions!
        trans_snow(i,K-1) = -snow_trans

        ! Use area-weighted remapped thicknesses so that the total ice area and
        ! mass are both conserved in the remapping operation.
        mH_ice(i,j,k-1) = (part_trans*mH_ice(i,j,k) + &
                           part_sz(i,j,k-1)*mH_ice(i,j,k-1)) / &
                          (part_trans + part_sz(i,j,k-1))

        ! h should be the first thing to correct via a non-constant profile, and
        ! can be improved independently from T & S.
        mH_ice(i,j,k) = mH_lim(k)
        part_sz(i,j,k-1) = part_sz(i,j,k-1) + part_trans
        part_sz(i,j,k) = part_sz(i,j,k) - part_trans

        mca_ice(i,k-1) = mca_ice(i,k-1) + mca_trans
        mca_ice(i,k) = mca_ice(i,k) - mca_trans

        mca_snow(i,k-1) = mca_snow(i,k-1) + snow_trans
        mca_snow(i,k) = mca_snow(i,k) - snow_trans

        mH_snow(i,j,k) = 0.0 ; mH_snow(i,j,k-1) = 0.0
        if (part_sz(i,j,k)>0.0) mH_snow(i,j,k) = mca_snow(i,k) / part_sz(i,j,k)
        if (part_sz(i,j,k-1)>0.0) mH_snow(i,j,k-1) = mca_snow(i,k-1) / part_sz(i,j,k-1)
      endif
    enddo ; enddo

    if (do_any) then
      call advect_tracers_thicker(mca0_ice, trans_ice, G, CS%SIS_tr_adv_CSp, &
                                  TrReg, .false., j, is, ie)
      call advect_tracers_thicker(mca0_snow, trans_snow, G, CS%SIS_tr_adv_CSp, &
                                  TrReg, .true., j, is, ie)
    endif

    ! Compress the ice in category 1 if it is thinner than the minimum.
    if (mH_lim(1) > 0.0) then
      do i=is,ie
        if ((mH_ice(i,j,1)*part_sz(i,j,1) > 0.0) .and. (mH_ice(i,j,1) < mH_lim(1))) then
          ! Compress the ice in this category to the minimum thickness.
          ! mH_snow(i,j,1) = mH_snow(i,j,1) * (mH_lim(1) / mH_ice(i,j,1))
          part_sz(i,j,1) = part_sz(i,j,1) * (mH_ice(i,j,1) * I_mH_lim1)
          mH_ice(i,j,1) = mH_lim(1)
          mH_snow(i,j,1) = mca_snow(i,1) / part_sz(i,j,1)
        endif
      enddo
    endif

    ! Recalculate _all_ snow thicknesses.
    !  This should be unnecessary but probably changes answers at the level of roundoff.
    do k=1,G%CatIce ; do i=is,ie
      if ((mca_snow(i,k)>0.0) .and. (part_sz(i,j,k) <= 0.0)) then
        call SIS_error(FATAL, "Nonzero snow mass exists on zero area.")
      endif
      if (mca_snow(i,k)>0.0) then
        mH_snow(i,j,k) = mca_snow(i,k) / part_sz(i,j,k)
      else
        mH_snow(i,j,k) = 0.0
      endif
    enddo ; enddo

  endif ; enddo  ! j-loop and do_j

end subroutine adjust_ice_categories

subroutine compress_ice(part_sz, mH_lim, mca_ice, mca_snow, mH_ice, mH_snow, &
                        TrReg, G, CS)
  type(sea_ice_grid_type),                      intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(G)), intent(inout) :: part_sz
  real, dimension(:),                           intent(in)    :: mH_lim  ! Move to grid type?
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),   intent(inout) :: mca_ice, mca_snow, mH_ice, mH_snow
  type(SIS_tracer_registry_type),               pointer       :: TrReg
  type(ice_transport_CS),                       pointer       :: CS
!   This subroutine compresses the ice, starting with the thinnest category, if
! the total fractional ice coverage exceeds 1.  It is assumed at the start that
! the sum over all categories (including ice free) of part_sz is 1, but that the
! part_sz of the ice free category may be negative to make this so.  In this
! routine, the mass (volume) is conserved, while the fractional coverage is
! solved for, while the new thicknesses are diagnosed.

!   This subroutine is effectively a minimalist version of a sea-ice ridging
! scheme.  A more complete ridging scheme would also compress thicker ice and
! allow the fractional ice coverage to drop below 1.

! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      mH_lim - The lower thickness limit of each category, in H (often kg m-2).
!  (inout)   mca_ice - The mass per unit grid-cell area of the ice in each
!                      category in H (often kg m-2).
!  (inout)   mca_snow - The mass per unit grid-cell area of the snow atop the
!                       ice in each category in H.
!  (inout)   mH_ice - The thickness of the ice in each category in H.
!  (inout)   mH_snow - The thickness of the snow atop the ice in each category
!                     in H.
!  (inout)   TrReg - The registry of registered SIS ice and snow tracers.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.

  real, dimension(SZI_(G),SZJ_(G)) :: excess_cover
  real :: compression_ratio
  real :: Icompress_here
  real :: mca_trans, mca_old
  real :: snow_trans, snow_old
  real :: Imca_new
  real, dimension(SZI_(G),SZCAT_(G)) :: &
    mca0_ice, mca0_snow, trans_ice, trans_snow
  logical :: do_any, do_j(SZJ_(G))
  character(len=200) :: mesg
  integer :: i, j, k, m, isc, iec, jsc, jec, nLay
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

  do j=jsc,jec ; if (do_j(j)) then
    do k=1,G%CatIce ; do i=isc,iec
      mca0_ice(i,k) = mca_ice(i,j,k) ; mca0_snow(i,k) = mca_snow(i,j,k)
    enddo ; enddo
    trans_ice(:,:) = 0.0 ; trans_snow(:,:) = 0.0
    do_any = .false.

    do k=1,G%CatIce-1 ; do i=isc,iec
      if ((excess_cover(i,j) > 0.0) .and. (mca_ice(i,j,k) > 0.0)) then
        compression_ratio = mH_ice(i,j,k) / mH_lim(k+1)
        if (part_sz(i,j,k)*(1.0-compression_ratio) >= excess_cover(i,j)) then
          ! This category is compacted, but not to the point that it needs to
          ! be transferred to a thicker layer.
          Icompress_here = part_sz(i,j,k) / (part_sz(i,j,k) - excess_cover(i,j))
          mH_ice(i,j,k) = mH_ice(i,j,k) * Icompress_here
          mH_snow(i,j,k) = mH_snow(i,j,k) * Icompress_here
          part_sz(i,j,k) = part_sz(i,j,k) - excess_cover(i,j)
          excess_cover(i,j) = 0.0
        else
          ! Mass from this category needs to be transfered to the next thicker
          ! category after being compacted to thickness mH_lim(k+1).
          excess_cover(i,j) = excess_cover(i,j) - part_sz(i,j,k)*(1.0-compression_ratio)
          part_sz(i,j,k+1) = part_sz(i,j,k+1) + part_sz(i,j,k)*compression_ratio

          mca_trans = mca_ice(i,j,k) ; mca_old = mca_ice(i,j,k+1)
          trans_ice(i,K) = mca_trans ; do_any = .true.
          mca_ice(i,j,k+1) = mca_ice(i,j,k+1) + mca_trans
          Imca_new = 1.0 / mca_ice(i,j,k+1)
  !        This is not quite right, or at least not consistent.
  !        mH_ice(i,j,k+1) = (mca_trans*mH_lim(k+1) + &
  !                          mca_old*mH_ice(i,j,k+1)) * Imca_new
          mH_ice(i,j,k+1) = mca_ice(i,j,k+1) / part_sz(i,j,k+1)

          if (mca_snow(i,j,k) > 0.0) then
            snow_trans = mca_snow(i,j,k) ; snow_old = mca_snow(i,j,k+1)
            trans_snow(i,K) = snow_trans
            mca_snow(i,j,k+1) = mca_snow(i,j,k+1) + mca_snow(i,j,k)
          endif
          mH_snow(i,j,k+1) = mca_snow(i,j,k+1) / part_sz(i,j,k+1)

          mca_ice(i,j,k) = 0.0 ; mca_snow(i,j,k) = 0.0 ; part_sz(i,j,k) = 0.0
          mH_ice(i,j,k) = 0.0 ; mH_snow(i,j,k) = 0.0
        endif
      endif
    enddo ; enddo

    if (do_any) then
      call advect_tracers_thicker(mca0_ice, trans_ice, G, CS%SIS_tr_adv_CSp, &
                                  TrReg, .false., j, isc, iec)
      call advect_tracers_thicker(mca0_snow, trans_snow, G, CS%SIS_tr_adv_CSp, &
                                  TrReg, .true., j, isc, iec)
    endif
  endif ; enddo ! j-loop.

  k=G%CatIce
  do j=jsc,jec ; if (do_j(j)) then ; do i=isc,iec
    if (excess_cover(i,j) > 0.0) then
      if (part_sz(i,j,k) <= 1.0) &
        call SIS_error(FATAL, &
            "Category CatIce part_sz inconsistent with excess cover.")
      Icompress_here = part_sz(i,j,k) / (part_sz(i,j,k) - excess_cover(i,j))
      mH_ice(i,j,k) = mH_ice(i,j,k) * Icompress_here
      mH_snow(i,j,k) = mH_snow(i,j,k) * Icompress_here
      part_sz(i,j,k) = part_sz(i,j,k) - excess_cover(i,j)
      excess_cover(i,j) = 0.0
    endif
  enddo ; endif ; enddo

  if (CS%check_conservation) then
    ! Check for consistency between mca_ice, mH_ice, and part_sz.
    do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
      if ((mca_ice(i,j,k) == 0.0) .and. (mH_ice(i,j,k)*part_sz(i,j,k) /= 0.0)) then
        write(mesg,'("Compress mismatch at ",3(i8),": mca, h, part, hxp = zero, ",3(1pe15.6))') &
           i, j, k, mH_ice(i,j,k), part_sz(i,j,k), mH_ice(i,j,k)*part_sz(i,j,k)
        call SIS_error(WARNING, mesg, all_print=.true.)
      endif
      if (abs(mca_ice(i,j,k) - mH_ice(i,j,k)*part_sz(i,j,k)) > 1e-12*mca_ice(i,j,k)) then
        write(mesg,'("Compress mismatch at ",3(i8),": mca, h, part, hxp = ",4(1pe15.6))') &
           i, j, k, mca_ice(i,j,k), mH_ice(i,j,k), part_sz(i,j,k), mH_ice(i,j,k)*part_sz(i,j,k)
        call SIS_error(WARNING, mesg, all_print=.true.)
      endif
    enddo ; enddo ; enddo
  endif

end subroutine compress_ice


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
        uflx(I,j) = uc(I,j) * trc(i,j) * G%dy_Cu(I,j)
      else
        uflx(I,j) = uc(I,j) * trc(i+1,j) * G%dy_Cu(I,j)
      endif
    enddo ; enddo

    do J=jsc-1,jec ; do i=isc,iec
      avg = ( trc(i,j) + trc(i,j+1) )/2
      dif = trc(i,j+1) - trc(i,j)
      if( avg > stop_lim .and. vc(i,J) * dif > 0.0) then
        vflx(i,J) = 0.0
      else if( vc(i,J) > 0.0 ) then
        vflx(i,J) = vc(i,J) * trc(i,j) * G%dx_Cv(i,J)
      else
        vflx(i,J) = vc(i,J) * trc(i,j+1) * G%dx_Cv(i,J)
      endif
    enddo ; enddo

    do j=jsc,jec ; do i=isc,iec
      trc(i,j) = trc(i,j) + dt_adv * ((uflx(I-1,j) - uflx(I,j)) + &
                                      (vflx(i,J-1) - vflx(i,J)) ) * G%IareaT(i,j)
    enddo ; enddo

    call pass_var(trc, G%Domain)
  enddo

end subroutine slab_ice_advect

subroutine get_total_amounts(mH_ice, mH_snow, part_sz, G, tot_ice, tot_snow)
  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),   intent(in)  :: mH_ice, mH_snow
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(G)), intent(in)  :: part_sz
  type(EFP_type), intent(out) :: tot_ice, tot_snow
! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      mH_ice - The mass per unit area of the ice in each
!                    category in units of H (often kg m-2).
!  (in)      mH_snow - The mass per unit area of the snow atop the
!                     ice in each category in units of H (often kg m-2).
!  (in)      G - The ocean's grid structure.
!  (out)     tot_ice - The globally integrated total ice, in kg.
!  (out)     tot_snow - The globally integrated total snow, in kg.

  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: sum_mca_ice, sum_mca_snow
  real :: total
  integer :: i, j, k, m, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  sum_mca_ice(:,:) = 0.0
  sum_mca_snow(:,:) = 0.0
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_mca_ice(i,j) = sum_mca_ice(i,j) + G%areaT(i,j) * (part_sz(i,j,k)*mH_ice(i,j,k))
    sum_mca_snow(i,j) = sum_mca_snow(i,j) + G%areaT(i,j) * (part_sz(i,j,k)*mH_snow(i,j,k))
  enddo ; enddo ; enddo

  total = reproducing_sum(sum_mca_ice, EFP_sum=tot_ice)
  total = reproducing_sum(sum_mca_snow, EFP_sum=tot_snow)

end subroutine get_total_amounts

subroutine get_total_enthalpy(mH_ice, mH_snow, part_sz, TrReg, &
                              G, enth_ice, enth_snow)
  type(sea_ice_grid_type),                      intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(G)),   intent(in)  :: mH_ice, mH_snow
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(G)), intent(in)  :: part_sz
  type(SIS_tracer_registry_type),               pointer     :: TrReg
  type(EFP_type),                               intent(out) :: enth_ice, enth_snow
! Arguments: part_sz - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1 at the end, in/out.
!  (in)      mH_ice - The mass per unit area of the ice in each
!                    category in units of H (often kg m-2).
!  (in)      mH_snow - The mass per unit area of the snow atop the
!                     ice in each category in units of H (often kg m-2).
!  (in)      TrReg - The registry of registered SIS ice and snow tracers.
!  (in)      G - The ocean's grid structure.
!  (out)     enth_ice - The globally integrated total ice enthalpy in J.
!  (out)     enth_snow - The globally integrated total snow enthalpy in J.

  real, dimension(:,:,:,:), pointer :: &
    heat_ice=>NULL(), & ! Pointers to the enth_ice and enth_snow arrays from the
    heat_snow=>NULL()   ! SIS tracer registry.  enth_ice is the enthalpy of the
                        ! ice in each category and layer, while enth_snow is the
                        ! enthalpy of the snow atop the ice in each category.
                        ! Both are in enth_units (J or rescaled).
  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: sum_enth_ice, sum_enth_snow
  real :: total, I_Nk
  integer :: i, j, k, m, isc, iec, jsc, jec, nLay
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call get_SIS_tracer_pointer("enth_ice", TrReg, heat_ice, nLay)
  call get_SIS_tracer_pointer("enth_snow", TrReg, heat_snow, nLay)
  sum_enth_ice(:,:) = 0.0 ; sum_enth_snow(:,:) = 0.0

  I_Nk = 1.0 / G%NkIce
  do m=1,G%NkIce ; do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_enth_ice(i,j) = sum_enth_ice(i,j) + (G%areaT(i,j) * &
              ((mH_ice(i,j,k)*part_sz(i,j,k))*I_Nk)) * heat_ice(i,j,k,m)
  enddo ; enddo ; enddo ; enddo
  do k=1,G%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_enth_snow(i,j) = sum_enth_snow(i,j) + (G%areaT(i,j) * &
              (mH_snow(i,j,k)*part_sz(i,j,k))) * heat_snow(i,j,k,1)
  enddo ; enddo ; enddo

  total = reproducing_sum(sum_enth_ice, EFP_sum=enth_ice)
  total = reproducing_sum(sum_enth_snow, EFP_sum=enth_snow)

end subroutine get_total_enthalpy

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

  call get_param(param_file, mod, "SEA_ICE_ROLL_FACTOR", CS%Roll_factor, &
                 "A factor by which the propensity of small amounts of \n"//&
                 "thick sea-ice to become thinner by rolling is increased \n"//&
                 "or 0 to disable rolling.  This can be thought of as the \n"//&
                 "minimum number of ice floes in a grid cell divided by \n"//&
                 "the horizontal floe aspect ratio.  Sensible values are \n"//&
                 "0 (no rolling) or larger than 1.", units="Nondim",default=1.0)

  call get_param(param_file, mod, "USE_SLAB_ICE", CS%SLAB_ICE, &
                 "If true, use the very old slab-style ice.", default=.false.)
  call get_param(param_file, mod, "CHECK_ICE_TRANSPORT_CONSERVATION", CS%check_conservation, &
                 "If true, use add multiple diagnostics of ice and snow \n"//&
                 "mass conservation in the sea-ice transport code.  This \n"//&
                 "is expensive and should be used sparingly.", default=.false.)
  call get_param(param_file, mod, "DO_RIDGING", CS%do_ridging, &
                 "Apply a ridging scheme as imported by Torge Martin.", default=.false.)
  call get_param(param_file, mod, "USE_SIS_CONTINUITY", CS%use_SIS_continuity, &
                 "If true, uses a continuity solver from SIS that has more options.", default=.true.)
  call get_param(param_file, mod, "USE_SIS_THICKNESS_ADVECTION", CS%use_SIS_thickness_advection, &
                 "If true, uses the SIS tracer transport scheme for thickness.", default=.false.)

  call SIS_continuity_init(Time, G, param_file, diag, CS%continuity_CSp)
  call SIS_tracer_advect_init(Time, G, param_file, diag, CS%SIS_tr_adv_CSp)

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

  call SIS_continuity_end(CS%continuity_CSp)
  call SIS_tracer_advect_end(CS%SIS_tr_adv_CSp)

  deallocate(CS)
end subroutine ice_transport_end

end module ice_transport_mod
