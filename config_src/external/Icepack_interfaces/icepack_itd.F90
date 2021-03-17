      module icepack_itd

      use icepack_kinds
!      use icepack_parameters, only: c0, c1, c2, c3, c15, c25, c100, p1, p01, p001, p5, puny
!      use icepack_parameters, only: Lfresh, rhos, ice_ref_salinity, hs_min, cp_ice, Tocnfrz, rhoi
!      use icepack_parameters, only: rhosi, sk_l, hs_ssl, min_salin
!      use icepack_tracers,    only: nt_Tsfc, nt_qice, nt_qsno, nt_aero, nt_isosno, nt_isoice
!      use icepack_tracers,    only: nt_apnd, nt_hpnd, nt_fbri, tr_brine, nt_bgc_S, bio_index
!      use icepack_tracers,    only: n_iso
!      use icepack_tracers,    only: tr_iso
!      use icepack_tracers,    only: icepack_compute_tracers
!      use icepack_parameters, only: solve_zsal, skl_bgc, z_tracers
!      use icepack_parameters, only: kcatbound, kitd
!      use icepack_therm_shared, only: Tmin, hi_min
!      use icepack_warnings,   only: warnstr, icepack_warnings_add
!      use icepack_warnings,   only: icepack_warnings_setabort, icepack_warnings_aborted
!      use icepack_zbgc_shared,only: zap_small_bgc

      implicit none

      private
      public :: cleanup_itd ,icepack_init_itd



!=======================================================================

      contains

!=======================================================================

      subroutine icepack_init_itd(ncat, hin_max)

      integer (kind=int_kind), intent(in) :: &
           ncat ! number of thickness categories

      real (kind=dbl_kind), intent(out) :: &
           hin_max(0:ncat)  ! category limits (m)


      end subroutine icepack_init_itd

      subroutine cleanup_itd (dt,          ntrcr,      &
                              nilyr,       nslyr,      &
                              ncat,        hin_max,    &
                              aicen,       trcrn,      &
                              vicen,       vsnon,      &
                              aice0,       aice,       &
                              n_aero,                  &
                              nbtrcr,      nblyr,      &
                              tr_aero,                 &
                              tr_pond_topo,            &
                              heat_capacity,           &
                              first_ice,               &
                              trcr_depend, trcr_base,  &
                              n_trcr_strata,nt_strata, &
                              fpond,       fresh,      &
                              fsalt,       fhocn,      &
                              faero_ocn,   fiso_ocn,   &
                              fzsal,                   &
                              flux_bio,    limit_aice_in)

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         nslyr , & ! number of snow layers
         ntrcr , & ! number of tracers in use
         nbtrcr, & ! number of bio tracers in use
         n_aero    ! number of aerosol tracers

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step

      real (kind=dbl_kind), dimension(0:ncat), intent(in) :: &
         hin_max   ! category boundaries (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! ice tracers

      real (kind=dbl_kind), intent(inout) :: &
         aice  , & ! total ice concentration
         aice0     ! concentration of open water

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      logical (kind=log_kind), intent(in) :: &
         tr_aero,      & ! aerosol flag
         tr_pond_topo, & ! topo pond flag
         heat_capacity   ! if false, ice and snow have zero heat capacity

      logical (kind=log_kind), dimension(ncat), intent(inout) :: &
         first_ice   ! For bgc and S tracers. set to true if zapping ice.

      ! ice-ocean fluxes (required for strict conservation)

      real (kind=dbl_kind), intent(inout), optional :: &
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean        (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean     (W/m^2)
         fzsal        ! net salt flux to ocean from zsalinity (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout), optional :: &
         flux_bio     ! net tracer flux to ocean from biology (mmol/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout), optional :: &
         faero_ocn    ! aerosol flux to ocean     (kg/m^2/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         fiso_ocn     ! isotope flux to ocean     (kg/m^2/s)

      logical (kind=log_kind), intent(in), optional ::   &
         limit_aice_in      ! if false, allow aice to be out of bounds
                            ! may want to allow this for unit tests


      end subroutine cleanup_itd

      end module icepack_itd

!=======================================================================
