      module icepack_mechred

      use icepack_kinds
      use icepack_tracers, only : n_iso, n_aero

      implicit none

      private
      public :: ridge_ice
      contains

      subroutine ridge_ice (dt,          ndtd,       &
                            ncat,        n_aero,     &
                            nilyr,       nslyr,      &
                            ntrcr,       hin_max,    &
                            rdg_conv,    rdg_shear,  &
                            aicen,       trcrn,      &
                            vicen,       vsnon,      &
                            aice0,                   &
                            trcr_depend, trcr_base,  &
                            n_trcr_strata,           &
                            nt_strata,               &
                            krdg_partic, krdg_redist,&
                            mu_rdg,      tr_brine,   &
                            dardg1dt,    dardg2dt,   &
                            dvirdgdt,    opening,    &
                            fpond,                   &
                            fresh,       fhocn,      &
                            faero_ocn,   fiso_ocn,   &
                            aparticn,    krdgn,      &
                            aredistn,    vredistn,   &
                            dardg1ndt,   dardg2ndt,  &
                            dvirdgndt,               &
                            araftn,      vraftn,     &
                            closing_flag,closing )

      integer (kind=int_kind), intent(in) :: &
         ndtd       , & ! number of dynamics subcycles
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nslyr , & ! number of snow layers
         n_aero, & ! number of aerosol tracers
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), intent(in) :: &
         mu_rdg , & ! gives e-folding scale of ridged ice (m^.5)
         dt             ! time step

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max   ! category limits (m)

      real (kind=dbl_kind), intent(in) :: &
         rdg_conv   , & ! normalized energy dissipation due to convergence (1/s)
         rdg_shear      ! normalized energy dissipation due to shear (1/s)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen      , & ! concentration of ice
         vicen      , & ! volume per unit area of ice          (m)
         vsnon          ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn          ! ice tracers

      real (kind=dbl_kind), intent(inout) :: &
         aice0          ! concentration of open water

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      integer (kind=int_kind), intent(in) :: &
         krdg_partic, & ! selects participation function
         krdg_redist    ! selects redistribution function

      logical (kind=log_kind), intent(in) :: &
         closing_flag, &! flag if closing is valid
         tr_brine       ! if .true., brine height differs from ice thickness

      ! optional history fields
      real (kind=dbl_kind), intent(inout), optional :: &
         dardg1dt   , & ! rate of fractional area loss by ridging ice (1/s)
         dardg2dt   , & ! rate of fractional area gain by new ridges (1/s)
         dvirdgdt   , & ! rate of ice volume ridged (m/s)
         opening    , & ! rate of opening due to divergence/shear (1/s)
         closing    , & ! rate of closing due to divergence/shear (1/s)
         fpond      , & ! fresh water flux to ponds (kg/m^2/s)
         fresh      , & ! fresh water flux to ocean (kg/m^2/s)
         fhocn          ! net heat flux to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         dardg1ndt  , & ! rate of fractional area loss by ridging ice (1/s)
         dardg2ndt  , & ! rate of fractional area gain by new ridges (1/s)
         dvirdgndt  , & ! rate of ice volume ridged (m/s)
         aparticn   , & ! participation function
         krdgn      , & ! mean ridge thickness/thickness of ridging ice
         araftn     , & ! rafting ice area
         vraftn     , & ! rafting ice volume
         aredistn   , & ! redistribution function: fraction of new ridge area
         vredistn       ! redistribution function: fraction of new ridge volume

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         faero_ocn      ! aerosol flux to ocean (kg/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         fiso_ocn       ! isotope flux to ocean (kg/m^2/s)

      ! local variables

      real (kind=dbl_kind), dimension (ncat) :: &
         eicen          ! energy of melting for each ice layer (J/m^2)

      real (kind=dbl_kind), dimension (ncat) :: &
         esnon, &       ! energy of melting for each snow layer (J/m^2)
         vbrin, &       ! ice volume with defined by brine height (m)
         sicen          ! Bulk salt in h ice (ppt*m)

      real (kind=dbl_kind) :: &
         asum       , & ! sum of ice and open water area
         aksum      , & ! ratio of area removed to area ridged
         msnow_mlt  , & ! mass of snow added to ocean (kg m-2)
         esnow_mlt  , & ! energy needed to melt snow in ocean (J m-2)
         mpond      , & ! mass of pond added to ocean (kg m-2)
         closing_net, & ! net rate at which area is removed    (1/s)
                        ! (ridging ice area - area of new ridges) / dt
         divu_adv   , & ! divu as implied by transport scheme  (1/s)
         opning     , & ! rate of opening due to divergence/shear
                        ! opning is a local variable;
                        ! opening is the history diagnostic variable
         ardg1      , & ! fractional area loss by ridging ice
         ardg2      , & ! fractional area gain by new ridges
         virdg      , & ! ice volume ridged
         aopen          ! area opening due to divergence/shear

      real (kind=dbl_kind), dimension (n_aero) :: &
         maero          ! aerosol mass added to ocean (kg m-2)

      real (kind=dbl_kind), dimension (n_iso) :: &
         miso          ! isotope mass added to ocean (kg m-2)

      real (kind=dbl_kind), dimension (0:ncat) :: &
         apartic        ! participation function; fraction of ridging
                        ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (ncat) :: &
         hrmin      , & ! minimum ridge thickness
         hrmax      , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp      , & ! ridge e-folding thickness (krdg_redist = 1)
         krdg       , & ! mean ridge thickness/thickness of ridging ice
         ardg1n     , & ! area of ice ridged
         ardg2n     , & ! area of new ridges
         virdgn     , & ! ridging ice volume
         mraftn         ! rafting ice mask

      real (kind=dbl_kind) :: &
         vice_init, vice_final, & ! ice volume summed over categories
         vsno_init, vsno_final, & ! snow volume summed over categories
         eice_init, eice_final, & ! ice energy summed over layers
         vbri_init, vbri_final, & ! ice volume in fbri*vicen summed over categories
         sice_init ,sice_final, & ! ice bulk salinity summed over categories
         esno_init, esno_final    ! snow energy summed over layers

      integer (kind=int_kind), parameter :: &
         nitermax = 20  ! max number of ridging iterations

      integer (kind=int_kind) :: &
         n          , & ! thickness category index
         niter      , & ! iteration counter
         k          , & ! vertical index
         it             ! tracer index

      real (kind=dbl_kind) :: &
         dti            ! 1 / dt

      logical (kind=log_kind) :: &
         iterate_ridging ! if true, repeat the ridging

      character (len=char_len) :: &
         fieldid        ! field identifier



      end subroutine ridge_ice


      end module icepack_mechred

!=======================================================================
