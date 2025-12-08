!> A full implementation of Icepack ponds parameterizations (to come)
module ice_ponds_mod

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use MOM_domains,       only : pass_var, pass_vector, BGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, log_param, read_param, log_version, param_file_type
use MOM_unit_scaling,  only : unit_scale_type
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_types,         only : ice_state_type, ist_chksum
use SIS_tracer_registry, only : SIS_tracer_registry_type, SIS_tracer_type, get_SIS_tracer_pointer
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs
use ice_grid,          only : ice_grid_type
!Icepack modules
use icepack_kinds
use icepack_itd, only: icepack_init_itd, cleanup_itd
use icepack_meltpond_lvl,  only: compute_ponds_lvl
use icepack_meltpond_sealvl,  only: compute_ponds_sealvl
use icepack_meltpond_topo, only: compute_ponds_topo
use icepack_warnings, only: icepack_warnings_flush, icepack_warnings_aborted, &
                            icepack_warnings_setabort
use icepack_tracers, only: icepack_init_tracer_indices, icepack_init_tracer_sizes
use icepack_tracers, only: icepack_query_tracer_sizes
use icepack_parameters, only : icepack_init_parameters
implicit none ; private

#include <SIS2_memory.h>

public :: ice_ponds, ice_ponds_init

!> Ice ponds control structure
type, public :: ice_ponds_CS ; private
  logical :: &
  level_pond = .false., &      !< .true. = preferred ponds on level ice
  sealevel_pond = .false., &   !< .true. = ponds not draining below sealevel
  topo_pond = .false.          !< .true. = topographic ponds
  real :: area_underflow = 0.0 !< a non-dimesional fractional area underflow limit for the sea-ice
                               !! ponding scheme. This is defaulted to zero, but a reasonable
                               !! value might be 10^-26 which for a km square grid cell
                               !! would equate to an Angstrom scale ice patch.
end type ice_ponds_CS

contains

!> Initialize the ice ponds
subroutine ice_ponds_init(G, IG, PF, CS, US)
  type(SIS_hor_grid_type),    intent(in) :: G      !<  G The ocean's grid structure.
  type(ice_grid_type),        intent(in) :: IG     !<   The sea-ice-specific grid structure.
  type(param_file_type),      intent(in) :: PF     !< A structure to parse for run-time parameters
  type(ice_ponds_CS),       pointer    :: CS     !< The ponds control structure.
  type(unit_scale_type),      intent(in) :: US     !< A structure with unit conversion factors.

  integer (kind=int_kind) :: ntrcr, ncat, nilyr, nslyr, nblyr, nfsd, n_iso, n_aero
  integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_alvl, nt_vlvl, nt_qsno
  character(len=40) :: mdl = "ice_ponds_init" ! This module's name.

  if (.not.associated(CS)) allocate(CS)
  call get_param(PF, mdl, "MELTPOND_LEVEL", CS%level_pond, &
                 "Use level melt ponds - not supported yet", default=.false.)
  call get_param(PF, mdl, "MELTPOND_SEALEVEL", CS%sealevel_pond, &
                 "Use sealevel melt ponds", default=.false.)
  call get_param(PF, mdl, "MELTPOND_TOPO", CS%topo_pond, &
                 "Use topographic melt ponds - not supported yet", default=.false.)

  ncat = IG%CatIce ! The number of sea-ice thickness categories
  nilyr = IG%NkIce ! The number of ice layers per category
  nslyr = IG%NkSnow ! The number if snow layers per category
  nblyr = 0 ! The number of bio/brine layers per category
  nfsd = 0 ! The number of floe size distribution layers
  n_iso = 0 ! The number of isotopes in use
  n_aero = 0 ! The number of aerosols in use
  nt_Tsfc = 1 ! Tracer index for ice/snow surface temperature
  nt_qice = 2 ! Starting index for ice enthalpy in layers
  nt_qsno = 2 + nilyr ! Starting index for snow enthalpy
  nt_sice = 2 + nilyr + nslyr ! Index for ice salinity
! nt_alvl=2+2*nilyr+nslyr ! Index for level ice fraction
! nt_vlvl=3+2*nilyr+nslyr ! Index for level ice volume fraction
  ntrcr = 2 + 2 * nilyr + nslyr ! number of tracers in use
  ! (1,2) snow/ice surface temperature +
  ! (3) ice salinity*nilyr  + (4) pond thickness

  call icepack_init_tracer_sizes(ntrcr_in=ntrcr, &
       ncat_in=ncat, nilyr_in=nilyr, nslyr_in=nslyr, nblyr_in=nblyr, &
       nfsd_in=nfsd, n_iso_in=n_iso, n_aero_in=n_aero)

  call icepack_init_tracer_indices(nt_Tsfc_in=nt_Tsfc, &
           nt_sice_in=nt_sice, nt_qice_in=nt_qice, nt_qsno_in=nt_qsno)
!          nt_alvl_in=nt_alvl, nt_vlvl_in=nt_vlvl )

! call icepack_init_parameters(mu_rdg_in=CS%mu_rdg, conserv_check_in=.true.)

end subroutine ice_ponds_init

!> ice_ponds is a wrapper for the Icepack pond routines
subroutine ice_ponds(IST, G, IG, mca_ice, mca_snow, mca_pond, TrReg, CS, US, dt, &
                       rdg_rate, rdg_height)
  type(ice_state_type),              intent(inout) :: IST !< A type describing the state of the sea ice.
  type(SIS_hor_grid_type),           intent(inout) :: G   !< G The ocean's grid structure.
  type(ice_grid_type),               intent(inout) :: IG  !< The sea-ice-specific grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), intent(inout) :: mca_ice  !< mass of ice?
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), intent(inout) :: mca_snow !< mass of snow?
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), intent(inout) :: mca_pond !< mass of pond water?
  type(SIS_tracer_registry_type),    pointer       :: TrReg  !< TrReg - The registry of registered SIS ice and
                                                          !! snow tracers.
  type(ice_ponds_CS),                intent(in)    :: CS  !< The ponds control structure.
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors.
  real,                              intent(in)    :: dt  !< The amount of time over which the ice dynamics are to be.
                                                          !!    advanced in seconds. [T ~> s]
  real, dimension(SZI_(G),SZJ_(G)), intent(out), optional :: rdg_rate !< Diagnostic of the rate of fractional
                                                              !! area loss-gain due to ridging (1/s)
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), intent(inout), optional :: rdg_height !< A diagnostic of the ridged ice
                                                              !! height [Z ~> m]

! logical,                            intent(in)    :: dyn_Cgrid !<  True if using C-grid velocities, B-grid if False.

  real :: dt_sec ! timestep in seconds

  integer :: i, j, k, n ! loop vars
  integer :: isc, iec, jsc, jec ! loop bounds

  integer :: &
       ncat  , & ! number of thickness categories
       nilyr , & ! number of ice layers
       nslyr     ! number of snow layers

  real, dimension(0:IG%CatIce) :: hin_max   ! category limits (m)

  real, dimension(IG%CatIce) :: &
       aicen, & ! concentration of ice
       vicen, & ! volume per unit area of ice          (m)
       vsnon, & ! volume per unit area of snow         (m)
       tr_tmp   ! for temporary storage
  ! ice tracers; ntr*(NkIce+NkSnow) guaranteed to be enough for all (intensive)
  real, dimension(4+2*IG%NkIce+IG%NkSnow,IG%CatIce) :: trcrn

  integer, dimension(4+2*IG%NkIce+IG%NkSnow) :: &
       trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon (weighting to use)
       n_trcr_strata  ! number of underlying tracer layers

  real, dimension(4+2*IG%NkIce+IG%NkSnow,3) :: &
       trcr_base      ! = 0 or 1 depending on tracer dependency
                    ! argument 2:  (1) aice, (2) vice, (3) vsno

  real :: meltt,  & ! top melt rate (m/s)
          melts,  & ! snow melt rate (m/s)
          frain,  & ! rainfall rate (kg/m2/s)
          Tair,   & ! air temperature (K)
          fsurfn, & ! atm-ice surface heat flux  (W/m2)
          Tsfcn,  & ! surface temperature (C)
          dhs,    & ! depth difference for snow on sea ice and pond ice
          ffrac,  & ! fraction of fsurfn over pond used to melt ipond
          meltsliqn,    & ! liquid contribution to meltponds in dt (kg/m^2)
          apnd, hpnd, ipnd, & ! pond tracers
          dpnd_freebdn, & ! pond drainage rate due to freeboard constraint (m/step)
          dpnd_dlidn,   & ! pond loss/gain due to ice lid (m/step)
          dpnd_flushn     ! pond flushing rate due to ice permeability (m/s)

  real, dimension(IG%NkIce) :: &
          qicen, &      ! ice layer enthalpy (J m-3)
          sicen         ! salinity (ppt)

  real :: rho_ice, rho_snow ! Density of ice and snow [R ~> kg m-3]
  real :: divu_adv
  integer :: m, ntrcr ! loop vars for tracer; n is tracer #; m is tracer layer
  integer :: nt_tsfc_in, nt_qice_in, nt_qsno_in, nt_sice_in
  integer :: nL_ice, nL_snow ! number of tracer levels
  integer :: ncat_out, ntrcr_out, nilyr_out, nslyr_out ! array sizes returned from Icepack query
  character(len=1256) :: mesg

  nSlyr = IG%NkSnow
  nIlyr = IG%NkIce
  nCat  = IG%CatIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  call get_SIS2_thermo_coefs(IST%ITV, rho_snow=rho_snow)
  dt_sec = dt*US%T_to_s

  call icepack_query_tracer_sizes(ncat_out=ncat_out, ntrcr_out=ntrcr_out, nilyr_out=nilyr_out, nslyr_out=nslyr_out)

  if (nIlyr .ne. nilyr_out .or. nSlyr .ne. nslyr_out ) &
    call SIS_error(FATAL, "Oops!! It looks like you are trying to use sea-ice ponds "//&
                          "but did not include the Icepack (https://github.com/CICE-Consortium/Icepack)"//&
                          "source code repository in your compilation procedure, and are instead using the default "//&
                          "stub routine contained in config_src/external. Adjust your compilation accordingly." )

  ! set category limits; Icepack has a max on the largest, unlimited, category (why?)

  hin_max(0)=0.0
  do k=1,nCat
    hin_max(k) = US%Z_to_m * IG%mH_cat_bound(k) / Rho_ice
  end do

! call get_SIS_tracer_pointer("level_area", TrReg, Tr_ice_alvl_ptr, 1)
! call get_SIS_tracer_pointer("level_mass", TrReg, Tr_ice_mlvl_ptr, 1)

!  call IST_chksum('before ice ponds ', IST, G, US, IG)

  do j=jsc,jec; do i=isc,iec
  if ((G%mask2dT(i,j) .gt. 0.0) .and. (sum(IST%part_size(i,j,1:nCat)) .gt. 0.0)) then
  ! feed locations to Icepack's ridge_ice

    ! start like we're putting ALL the snow and pond in the ocean
    IST%snow_to_ocn(i,j) = IST%snow_to_ocn(i,j) + sum(mca_snow(i,j,:))
    IST%enth_snow_to_ocn(i,j) = IST%enth_snow_to_ocn(i,j) + sum(mca_snow(i,j,:)*TrReg%Tr_snow(1)%t(i,j,:,1))
    IST%water_to_ocn(i,j) = IST%water_to_ocn(i,j) + sum(mca_pond(i,j,:))
    aicen(1:nCat) = IST%part_size(i,j,1:nCat)

    if (sum(aicen) .eq. 0.0) then ! no ice -> no ponds
      IST%part_size(i,j,0) = 1.0
    else
      ! set up ice and snow volumes
      vicen(1:nCat) = mca_ice(i,j,1:nCat) /Rho_ice * US%Z_to_m  ! volume per unit area of ice (m)
      vsnon(1:nCat) = mca_snow(i,j,1:nCat)/Rho_snow * US%Z_to_m ! volume per unit area of snow (m)

      ! call Icepack routine; how are ponds treated?
      do n=1,nCat
        if (CS%sealevel_pond) then
          ipnd = IST%mH_pond_ice(i,j,n)
          apnd = IST%mH_pond_ice(i,j,n)
          hpnd = IST%mH_pond(i,j,n)
          call compute_ponds_sealvl( dt_sec,                &
                                     meltt,  melts,  frain, &
                                     Tair,   fsurfn, Tsfcn, &
                                     dhs,    ffrac,         &
                                     aicen(n),  vicen(n),  vsnon(n), &
                                     qicen,  sicen,         &
                                     apnd,   hpnd,  ipnd,   &
                                     meltsliqn,             &
                                     dpnd_freebdn,          &
                                     dpnd_dlidn, dpnd_flushn)
!         call compute_ponds_lvl (dt=dt_sec,        &
!                               nilyr=nilyr,      &
!                               ktherm=ktherm,    &
!                               hi_min=hi_min,    &
!                               dpscale=dpscale,  &
!                               frzpnd=frzpnd,    &
!                               rfrac=rfrac,      &
!                               meltt=melttn (n), &
!                               melts=meltsn (n), &
!                               frain=frain,      &
!                               Tair=Tair,        &
!                               fsurfn=fsurfn(n), &
!                               dhs=dhsn     (n), &
!                               ffrac=ffracn (n), &
!                               aicen=aicen  (n), &
!                               vicen=vicen  (n), &
!                               vsnon=vsnon  (n), &
!                               qicen=zqin (:,n), &
!                               sicen=zSin (:,n), &
!                               Tsfcn=Tsfc   (n), &
!                               alvl=alvl    (n), &
!                               apnd=apnd    (n), &
!                               hpnd=hpnd    (n), &
!                               ipnd=ipnd    (n), &
!                               meltsliqn=l_meltsliqn(n))

!         call compute_ponds_topo(dt,       ncat,      nilyr,     &
!                               ktherm,                         &
!                               aice,     aicen,                &
!                               vice,     vicen,                &
!                               vsno,     vsnon,                &
!                               meltt,                &
!                               fsurf,    fpond,                &
!                               Tsfc,     Tf,                   &
!                               zqin,     zSin,                 &
!                               apnd,     hpnd,      ipnd       )
!       if (icepack_warnings_aborted(subname)) return
        endif
      enddo

      if ( icepack_warnings_aborted() ) then
        call icepack_warnings_flush(0)
        call icepack_warnings_setabort(.false.)
        call SIS_error(WARNING, 'icepack compute_ponds error')
      endif

      ! pop pond off top of stack
      tr_tmp(1:nCat)=trcrn(ntrcr,1:nCat)

      do k=1,nCat
        IST%mH_pond(i,j,k) = tr_tmp(k)
        mca_pond(i,j,k) = IST%mH_pond(i,j,k)*aicen(k)
      enddo

      ! ! output: snow/ice masses/thicknesses
      do k=1,nCat
        if (aicen(k) < CS%area_underflow) then
           aicen(k)=0.0
           vicen(k)=0.0
        endif
        if (aicen(k) > 0.0) then
          IST%part_size(i,j,k)  = aicen(k)
          mca_ice(i,j,k)  = vicen(k)*Rho_ice * US%m_to_Z
          IST%mH_ice(i,j,k)   = vicen(k)*Rho_ice/aicen(k) * US%m_to_Z
          mca_snow(i,j,k) = vsnon(k)*Rho_snow * US%m_to_Z
          IST%mH_snow(i,j,k)  = vsnon(k)*Rho_snow/aicen(k) * US%m_to_Z
        else
          IST%part_size(i,j,k) = 0.0
          mca_ice(i,j,k)  = 0.0
          IST%mH_ice(i,j,k) = 0.0
          mca_snow(i,j,k) = 0.0
          IST%mH_snow(i,j,k) = 0.0
       endif

     enddo

     IST%part_size(i,j,0) = 1.0 - sum(IST%part_size(i,j,1:nCat))

    endif
    ! subtract new snow/pond mass and energy on ice to sum net fluxes to ocean
    IST%snow_to_ocn(i,j) = IST%snow_to_ocn(i,j) - sum(mca_snow(i,j,:))
    IST%enth_snow_to_ocn(i,j) = IST%enth_snow_to_ocn(i,j) - sum(mca_snow(i,j,:)*TrReg%Tr_snow(1)%t(i,j,:,1))
    IST%water_to_ocn(i,j) = IST%water_to_ocn(i,j) - sum(mca_pond(i,j,:))

  endif; enddo; enddo ! part_sz, j, i

!  call IST_chksum('after ice ponds ', IST, G, US, IG)

end subroutine ice_ponds

!> ice_ponds_end deallocates the memory associated with this module.
subroutine ice_ponds_end()

end subroutine ice_ponds_end

end module ice_ponds_mod
