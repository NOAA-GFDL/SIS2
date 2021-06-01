!> A full implementation of Icepack ridging parameterizations
module ice_ridging_mod

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! replaced T. Martin code with wrapper for Icepack ridging function - mw 1/18  !
!                                                                              !
! Prioritized to do list as of 6/4/19 (mw):                                    !
!                                                                              !
! 1) implement new snow_to_ocean diagnostic to record this flux.               !
! 2) implement ridging_rate diagnostics: ridging_shear, ridging_conv           !
! 3) implement "do_j" style optimization as in "compress_ice" or               !
!    "adjust_ice_categories" (SIS_transport.F90) if deemed necessary           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use MOM_domains,       only : pass_var, pass_vector, BGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, log_param, read_param, log_version, param_file_type
use MOM_unit_scaling,  only : unit_scale_type
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_types,         only : ice_state_type, ist_chksum
use fms_io_mod,        only : register_restart_field, restart_file_type
use SIS_tracer_registry, only : SIS_tracer_registry_type, SIS_tracer_type, get_SIS_tracer_pointer
use SIS2_ice_thm,    only : get_SIS2_thermo_coefs
use ice_grid,          only : ice_grid_type
!Icepack modules
use icepack_kinds
use icepack_itd, only: icepack_init_itd, cleanup_itd
use icepack_mechred, only: ridge_ice
use icepack_warnings, only: icepack_warnings_flush, icepack_warnings_aborted, &
                            icepack_warnings_setabort
use icepack_tracers, only: icepack_init_tracer_indices, icepack_init_tracer_sizes
use icepack_tracers, only: icepack_query_tracer_sizes
use icepack_parameters, only : icepack_init_parameters
implicit none ; private

#include <SIS2_memory.h>



public :: ice_ridging, ridge_rate, ice_ridging_init

real, parameter :: hlim_unlim = 1.e8   !< Arbitrary huge number used in ice_ridging
real    :: s2o_frac       = 0.5        !< Fraction of snow dumped into ocean during ridging [nondim]
logical :: rdg_lipscomb = .true.       !< If true, use the Lipscomb ridging scheme
                                       !! TODO: These parameters belong in a control structure

! future namelist parameters?
integer (kind=int_kind), parameter :: &
     krdg_partic = 0  , & ! 1 = new participation, 0 = Thorndike et al 75
     krdg_redist = 0      ! 1 = new redistribution, 0 = Hibler 80

! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
real(kind=dbl_kind), parameter ::  mu_rdg = 3.0


contains



!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ridge_rate finds the deformation rate or total energy dissipation rate due
!!   to ridging (Flato and Hibler, 1995, JGR) or the net area loss in riding
!!   (CICE documentation) depending on the state of the ice drift
function ridge_rate(del2, div) result (rnet)
  real, intent(in)  :: del2 !< The magnitude squared of the shear rates [T-2 ~> s-2].
  real, intent(in)  :: div  !< The ice flow divergence [T-1 ~> s-1]

  ! Local variables
  real :: rnet ! The net rate of area loss or energy dissipation rate due to ridging [T-1 ~> s-1]
  real :: del, rconv, rshear ! Local variables with rates [T-1 ~> s-1]
  !TOM> cs is now set in namelist:
  !Niki: this was commented out
  real, parameter   :: cs=0.25 !(CICE documentation)

  del=sqrt(del2)

  rconv  = -min(div,0.0)           ! energy dissipated by convergence ...
  rshear = 0.5*(del-abs(div))      ! ... and by shear
  rnet   = rconv + cs*rshear       ! net energy contains only part of the
                                   !  shear energy as only a fraction is
           !  dissipated in the ridging process
end function ridge_rate


subroutine ice_ridging_init(G,IG,TrReg,US)
  type(SIS_hor_grid_type),                       intent(inout) :: G !<  G The ocean's grid structure.
  type(ice_grid_type),                           intent(inout) :: IG !<   The sea-ice-specific grid structure.
  type(SIS_tracer_registry_type),                pointer       :: TrReg  ! TrReg - The registry of registered SIS ice and snow tracers.
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors.

  integer (kind=int_kind) :: ntrcr, ncat, nilyr, nslyr, nblyr, nfsd, n_iso, n_aero
  integer (kind=int_kind) :: nt_Tsfc, nt_sice, nt_qice, nt_vlvl, nt_qsno

  ncat=IG%CatIce ! The number of sea-ice thickness categories
  nilyr=IG%NkIce ! The number of ice layers per category
  nslyr=IG%NkSnow ! The number if snow layers per category
  nblyr=0 ! The number of bio/brine layers per category
  nfsd=0 ! The number of floe size distribution layers
  n_iso=0 ! The number of isotopes in use
  n_aero=0 ! The number of aerosols in use
  nt_Tsfc=1 ! Tracer index for ice/snow surface temperatore
  nt_qice=2 ! Starting index for ice enthalpy in layers
  nt_qsno=2+nilyr ! Starting index for snow enthalpy
  nt_sice=2+nilyr+nslyr ! Index for ice salinity
  nt_vlvl=2+2*nilyr+nslyr ! Index for level ice volume fraction
  ntrcr=2+2*nilyr+nslyr ! number of tracers in use
  ! (1,2) snow/ice surface temperature +
  ! (3) ice salinity*nilyr  + (4) pond thickness

  call icepack_init_tracer_sizes(ntrcr_in=ntrcr, &
       ncat_in=ncat, nilyr_in=nilyr, nslyr_in=nslyr, nblyr_in=nblyr, &
       nfsd_in=nfsd, n_iso_in=n_iso, n_aero_in=n_aero)

  call icepack_init_tracer_indices(nt_Tsfc_in=nt_Tsfc, &
           nt_sice_in=nt_sice, nt_qice_in=nt_qice, &
           nt_qsno_in=nt_qsno, nt_vlvl_in=nt_vlvl )

  call icepack_init_parameters(mu_rdg_in=mu_rdg,conserv_check_in=.true.)

end subroutine ice_ridging_init
!
! ice_ridging is a wrapper for the icepack ridging routine ridge_ice
!
subroutine ice_ridging(IST, G, IG, mca_ice, mca_snow, mca_pond, TrReg, US, dt, rdg_rate)
  type(ice_state_type),              intent(inout) :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type),                       intent(inout) :: G !<  G The ocean's grid structure.
  type(ice_grid_type),                           intent(inout) :: IG !<   The sea-ice-specific grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),   intent(inout) :: mca_ice, mca_snow, mca_pond
  type(SIS_tracer_registry_type),                pointer       :: TrReg  ! TrReg - The registry of registered SIS ice and snow tracers.
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors.
  real (kind=dbl_kind),                          intent(in)    :: dt !<   The amount of time over which the ice dynamics are to be.
                                                                     !    advanced in seconds.
  real, dimension(SZI_(G),SZJ_(G)), intent(inout), optional :: rdg_rate !< Diagnostic of the rate of fractional
                                                                        !! area loss-gain due to ridging (1/s)

!  logical,                                       intent(in)    :: dyn_Cgrid !<  True if using C-gridd velocities, B-grid if False.



  ! these strain metrics are calculated here from the velocities used for advection
  real :: sh_Dt ! sh_Dt is the horizontal tension (du/dx - dv/dy) including
                ! all metric terms, in s-1.
  real :: sh_Dd ! sh_Dd is the flow divergence (du/dx + dv/dy) including all
                ! metric terms, in s-1.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    sh_Ds       ! sh_Ds is the horizontal shearing strain (du/dy + dv/dx)
                ! including all metric terms, in s-1.

  integer :: i, j, k ! loop vars
  integer isc, iec, jsc, jec ! loop bounds
  integer :: halo_sh_Ds  ! The halo size that can be used in calculating sh_Ds.

  integer (kind=int_kind) :: &
       ncat  , & ! number of thickness categories
       nilyr , & ! number of ice layers
       nslyr  ! number of snow layers


  real (kind=dbl_kind), dimension(0:IG%CatIce) :: hin_max   ! category limits (m)

  logical (kind=log_kind) :: &
       closing_flag, &! flag if closing is valid
       tr_brine       ! if .true., brine height differs from ice thickness

  ! optional history fields
  real (kind=dbl_kind) :: &
       dardg1dt   , & ! rate of fractional area loss by ridging ice (1/s)
       dardg2dt   , & ! rate of fractional area gain by new ridges (1/s)
       dvirdgdt   , & ! rate of ice volume ridged (m/s)
       opening    , & ! rate of opening due to divergence/shear (1/s)
       closing    , & ! rate of closing due to divergence/shear (1/s)
       fpond      , & ! fresh water flux to ponds (kg/m^2/s)
       fresh      , & ! fresh water flux to ocean (kg/m^2/s)
       fhocn          ! net heat flux to ocean (W/m^2)

  real (kind=dbl_kind), dimension(IG%CatIce) :: &
       dardg1ndt  , & ! rate of fractional area loss by ridging ice (1/s)
       dardg2ndt  , & ! rate of fractional area gain by new ridges (1/s)
       dvirdgndt  , & ! rate of ice volume ridged (m/s)
       aparticn   , & ! participation function
       krdgn      , & ! mean ridge thickness/thickness of ridging ice
       araftn     , & ! rafting ice area
       vraftn     , & ! rafting ice volume
       aredistn   , & ! redistribution function: fraction of new ridge area
       vredistn       ! redistribution function: fraction of new ridge volume

  real (kind=dbl_kind), dimension(IG%CatIce) :: &
       faero_ocn      ! aerosol flux to ocean (kg/m^2/s)

  real (kind=dbl_kind), dimension(IG%CatIce) :: &
       fiso_ocn       ! isotope flux to ocean (kg/m^2/s)

  integer (kind=int_kind) :: &
       ndtd = 1  , & ! number of dynamics subcycles
       n_aero = 0, & ! number of aerosol tracers
       ntrcr = 0     ! number of tracer level


  real(kind=dbl_kind) :: &
       del_sh        , & ! shear strain measure
       rdg_conv = 0.0, & ! normalized energy dissipation from convergence (1/s)
       rdg_shear= 0.0    ! normalized energy dissipation from shear (1/s)

  real(kind=dbl_kind), dimension(IG%CatIce) :: &
       aicen, & ! concentration of ice
       vicen, & ! volume per unit area of ice          (m)
       vsnon, &    ! volume per unit area of snow         (m)
       tr_tmp ! for temporary storade
  ! ice tracers; ntr*(NkIce+NkSnow) guaranteed to be enough for all (intensive)
  real(kind=dbl_kind), dimension(2+2*IG%NkIce+IG%NkSnow,IG%CatIce) :: trcrn

  real(kind=dbl_kind) :: aice0          ! concentration of open water

  integer (kind=int_kind), dimension(2+2*IG%NkIce+IG%NkSnow) :: &
       trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon (weighting to use)
       n_trcr_strata  ! number of underlying tracer layers

  real(kind=dbl_kind), dimension(2+2*IG%NkIce+IG%NkSnow,3) :: &
       trcr_base      ! = 0 or 1 depending on tracer dependency
                    ! argument 2:  (1) aice, (2) vice, (3) vsno

  integer(kind=int_kind), dimension(2+2*IG%NkIce+IG%NkSnow,IG%CatIce) :: &
       nt_strata      ! indices of underlying tracer layers

  type(SIS_tracer_type), dimension(:), pointer :: Tr=>NULL() ! SIS2 tracers
  real, dimension(:,:,:,:),       pointer    :: Tr_ice_enth_ptr=>NULL() !< A pointer to the named tracer
  real, dimension(:,:,:,:),       pointer    :: Tr_snow_enth_ptr=>NULL() !< A pointer to the named tracer
  real, dimension(:,:,:,:),       pointer    :: Tr_ice_salin_ptr=>NULL() !< A pointer to the named tracer

  real :: rho_ice, rho_snow, divu_adv
  integer :: m, n ! loop vars for tracer; n is tracer #; m is tracer layer
  integer(kind=int_kind) :: nt_tsfc_in, nt_qice_in, nt_qsno_in, nt_sice_in
  integer(kind=int_kind) :: nL_ice, nL_snow ! number of tracer levels
  integer(kind=int_kind) :: ncat_out, ntrcr_out, nilyr_out, nslyr_out ! array sizes returned from Icepack query

  nSlyr = IG%NkSnow
  nIlyr = IG%NkIce
  nCat  = IG%CatIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  call get_SIS2_thermo_coefs(IST%ITV, rho_snow=rho_snow)

  call icepack_query_tracer_sizes(ncat_out=ncat_out,ntrcr_out=ntrcr_out, nilyr_out=nilyr_out, nslyr_out=nslyr_out)


  if (nIlyr .ne. nilyr_out .or. nSlyr .ne. nslyr_out ) call SIS_error(FATAL,'nilyr or nslyr mismatch with Icepack')

  ! copy strain calculation code from SIS_C_dynamics; might be a more elegant way ...
  !
  halo_sh_Ds = min(isc-G%isd, jsc-G%jsd, 2)
  !  if (dyn_Cgrid) then
  do J=jsc-halo_sh_Ds,jec+halo_sh_Ds-1 ; do I=isc-halo_sh_Ds,iec+halo_sh_Ds-1
    ! This uses a no-slip boundary condition.
    sh_Ds(I,J) = (2.0-G%mask2dBu(I,J)) * &
         (G%dxBu(I,J)*G%IdyBu(I,J)*(IST%u_ice_C(I,j+1)*G%IdxCu(I,j+1) - IST%u_ice_C(I,j)*G%IdxCu(I,j)) + &
         G%dyBu(I,J)*G%IdxBu(I,J)*(IST%v_ice_C(i+1,J)*G%IdyCv(i+1,J) - IST%v_ice_C(i,J)*G%IdyCv(i,J)))
  enddo; enddo

  ! set category limits; Icepack has a max on the largest, unlimited, category (why?)

  hin_max(0)=0.0
  do k=1,nCat
    hin_max(k) = IG%mH_cat_bound(k)/(Rho_ice*IG%kg_m2_to_H)
  end do

  trcr_base = 0.0; n_trcr_strata = 0; nt_strata = 0; ! init some tracer vars
  ! When would we use icepack tracer "strata"?

  ! set icepack tracer index "nt_lvl" to (last) pond tracer so it gets dumped when
  ! ridging in ridge_ice (this is what happens to "level" ponds); first add up ntrcr;
  ! then set nt_lvl to ntrcr+1; could move this to an initializer - mw

  call get_SIS_tracer_pointer("enth_ice",TrReg,Tr_ice_enth_ptr,nL_ice)
  call get_SIS_tracer_pointer("enth_snow",TrReg,Tr_snow_enth_ptr,nL_snow)
  call get_SIS_tracer_pointer("salin_ice",TrReg,Tr_ice_salin_ptr,nL_ice)

!  call IST_chksum('before ice ridging ',IST,G,US,IG)

  if (present(rdg_rate)) rdg_rate(:,:)=0.0
  do j=jsc,jec; do i=isc,iec
  if ((G%mask2dT(i,j) .gt. 0.0) .and. (sum(IST%part_size(i,j,1:nCat)) .gt. 0.0)) then
  ! feed locations to Icepack's ridge_ice

    ! start like we're putting ALL the snow and pond in the ocean
    IST%snow_to_ocn(i,j) = IST%snow_to_ocn(i,j) + sum(mca_snow(i,j,:))
    IST%enth_snow_to_ocn(i,j) = IST%enth_snow_to_ocn(i,j) + sum(mca_snow(i,j,:)*TrReg%Tr_snow(1)%t(i,j,:,1));
    IST%water_to_ocn(i,j) = IST%water_to_ocn(i,j) + sum(mca_pond(i,j,:));
    aicen(1:nCat) = IST%part_size(i,j,1:nCat);


    if (sum(aicen) .eq. 0.0) then ! no ice -> no ridging
      IST%part_size(i,j,0) = 1.0;
    else
      ! set up ice and snow volumes
      vicen(1:nCat) = mca_ice(i,j,1:nCat) /Rho_ice ! volume per unit area of ice (m)
      vsnon(1:nCat) = mca_snow(i,j,1:nCat)/Rho_snow ! volume per unit area of snow (m)

      sh_Dt = (G%dyT(i,j)*G%IdxT(i,j)*(G%IdyCu(I,j) * IST%u_ice_C(I,j) - &
                                       G%IdyCu(I-1,j)*IST%u_ice_C(I-1,j)) - &
               G%dxT(i,j)*G%IdyT(i,j)*(G%IdxCv(i,J) * IST%v_ice_C(i,J) - &
                                       G%IdxCv(i,J-1)*IST%v_ice_C(i,J-1)))
      sh_Dd = (G%IareaT(i,j)*(G%dyCu(I,j) * IST%u_ice_C(I,j) - &
                              G%dyCu(I-1,j)*IST%u_ice_C(I-1,j)) + &
               G%IareaT(i,j)*(G%dxCv(i,J) * IST%v_ice_C(i,J) - &
                              G%dxCv(i,J-1)*IST%v_ice_C(i,J-1)))

      del_sh = sqrt(sh_Dd**2 + 0.25 * (sh_Dt**2 + &
                   (0.25 * ((sh_Ds(I-1,J-1) + sh_Ds(I,J)) + &
                            (sh_Ds(I-1,J) + sh_Ds(I,J-1))))**2 ) ) ! H&D eqn 9
      rdg_conv  = -min(sh_Dd,0.0)              ! energy dissipated by convergence ...
      rdg_shear = 0.5*(del_sh-abs(sh_Dd))      ! ... and by shear

      aice0 = IST%part_size(i,j,0)
      if (aice0<0.) then
         call SIS_error(WARNING,'aice0<0. before call to ridge ice.')
         aice0=0.
      endif

      ntrcr = 0
!      Tr_ptr=>NULL()
      if (TrReg%ntr>0) then ! load tracer array
        ntrcr=ntrcr+1
        trcrn(ntrcr,1) = Tr_ice_enth_ptr(i,j,1,1) ! surface temperature? this is taken from the thinnest ice category
        trcr_depend(ntrcr) = 1; ! 1 means ice-based tracer
        trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,2) = 1.0; ! 2nd index for ice
        do k=1,nL_ice
          ntrcr=ntrcr+1
          trcrn(ntrcr,1:ncat)=Tr_ice_enth_ptr(i,j,1:nCat,k)
          trcr_depend(ntrcr) = 1; ! 1 means ice-based tracer
          trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,2) = 1.0; ! 2nd index for ice
        enddo
        do k=1,nL_snow
          ntrcr=ntrcr+1
          trcrn(ntrcr,1:nCat)=Tr_snow_enth_ptr(i,j,1:nCat,k)
          trcr_depend(ntrcr) = 2; ! 2 means snow-based tracer
          trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,3) = 1.0; ! 3rd index for snow
        enddo
        do k=1,nL_ice
          ntrcr=ntrcr+1
          trcrn(ntrcr,1:nCat)=Tr_ice_salin_ptr(i,j,1:nCat,k)
          trcr_depend(ntrcr) = 1; ! 1 means ice-based tracer
          trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,2) = 1.0; ! 2nd index for ice
        enddo
      endif ! have tracers to load

      ! load pond on top of stack
      ntrcr = ntrcr + 1
      trcrn(ntrcr,1:nCat) = IST%mH_pond(i,j,1:nCat)
      trcr_depend(ntrcr) = 0; ! 0 means ice area-based tracer
      trcr_base(ntrcr,:) = 0.0; trcr_base(ntrcr,1) = 1.0; ! 1st index for ice area

      if (ntrcr .ne. ntrcr_out ) call SIS_error(FATAL,'ntrcr mismatch with Icepack')

      tr_brine = .false.
      dardg1dt=0.0
      dardg2dt=0.0
      opening=0.0
      fpond=0.0
      fresh=0.0
      fhocn=0.0
      faero_ocn=0.0
      fiso_ocn=0.0
      aparticn=0.0
      krdgn=0.0
      aredistn=0.0
      vredistn=0.0
      dardg1ndt=0.0
      dardg2ndt=0.0
      dvirdgndt=0.0
      araftn=0.0
      vraftn=0.0
      closing_flag=.false.

      ! call Icepack routine; how are ponds treated?
      call ridge_ice (dt,           ndtd,           &
                      ncat,         n_aero,         &
                      nilyr,        nslyr,          &
                      ntrcr,        hin_max,        &
                      rdg_conv,     rdg_shear,      &
                      aicen,                        &
                      trcrn,                        &
                      vicen,        vsnon,          &
                      aice0,                        &
                      trcr_depend,                  &
                      trcr_base,                    &
                      n_trcr_strata,                &
                      nt_strata,                    &
                      krdg_partic,  krdg_redist,    &
                      mu_rdg,       tr_brine,       &
                      dardg1dt=dardg1dt,     dardg2dt=dardg2dt,       &
                      dvirdgdt=dvirdgdt,     opening=opening,        &
                      fpond=fpond,                        &
                      fresh=fresh,        fhocn=fhocn,          &
                      faero_ocn=faero_ocn,   fiso_ocn=fiso_ocn,   &
                      aparticn=aparticn,  &
                      krdgn=krdgn,             &
                      aredistn=aredistn,       &
                      vredistn=vredistn,       &
                      dardg1ndt=dardg1ndt,     &
                      dardg2ndt=dardg2ndt,      &
                      dvirdgndt=dvirdgndt,                    &
                      araftn=araftn,       &
                      vraftn=vraftn,         &
                      closing_flag=closing_flag ,closing=closing)

      rdg_rate(i,j) = dardg1dt-dardg2dt

      if ( icepack_warnings_aborted() ) then
        call icepack_warnings_flush(0);
        call icepack_warnings_setabort(.false.)
        call SIS_error(WARNING,'icepack ridge_ice error');
      endif

      ! pop pond off top of stack
      tr_tmp(1:nCat)=trcrn(ntrcr,1:nCat)

      do k=1,nCat
        IST%mH_pond(i,j,k) = tr_tmp(k)
        mca_pond(i,j,k) = IST%mH_pond(i,j,k)*aicen(k)
      enddo

      if (TrReg%ntr>0) then
        ! unload tracer array reversing order of load -- stack-like fashion

         do k=1,nL_ice
           tr_tmp(1:nCat)=trcrn(1+k,1:nCat)
           Tr_ice_enth_ptr(i,j,1:nCat,k) = tr_tmp(1:nCat)
         enddo

         do k=1,nL_snow
           tr_tmp(1:nCat)=trcrn(1+nl_ice+k,1:ncat)
           Tr_snow_enth_ptr(i,j,1:nCat,k) = tr_tmp(1:nCat)
         enddo

         do k=1,nL_ice
           tr_tmp(1:nCat)=trcrn(1+k+nl_ice+nl_snow,1:nCat)
           Tr_ice_salin_ptr(i,j,1:nCat,k) =  tr_tmp(1:nCat)
         enddo

      endif ! have tracers to unload

      ! ! output: snow/ice masses/thicknesses
      do k=1,nCat

        if (aicen(k) > 0.0) then
          IST%part_size(i,j,k)  = aicen(k)
          mca_ice(i,j,k)  = vicen(k)*Rho_ice
          IST%mH_ice(i,j,k)   = vicen(k)*Rho_ice/aicen(k)
          mca_snow(i,j,k) = vsnon(k)*Rho_snow
          IST%mH_snow(i,j,k)  = vsnon(k)*Rho_snow/aicen(k)
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
    IST%snow_to_ocn(i,j) = IST%snow_to_ocn(i,j) - sum(mca_snow(i,j,:));
    IST%enth_snow_to_ocn(i,j) = IST%enth_snow_to_ocn(i,j) - sum(mca_snow(i,j,:)*TrReg%Tr_snow(1)%t(i,j,:,1));
    IST%water_to_ocn(i,j) = IST%water_to_ocn(i,j) - sum(mca_pond(i,j,:));

  endif; enddo; enddo ! part_sz, j, i

!  call IST_chksum('after ice ridging ',IST,G,US,IG)

end subroutine ice_ridging


!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_ridging_end deallocates the memory associated with this module.
subroutine ice_ridging_end()

end subroutine ice_ridging_end


end module ice_ridging_mod
