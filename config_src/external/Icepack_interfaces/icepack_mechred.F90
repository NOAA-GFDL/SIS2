      module icepack_mechred

      use icepack_kinds
      use icepack_tracers, only : n_iso, n_aero

      implicit none

      private
      public :: ridge_ice
      contains

!-----------------------------------------------------------------------
!> Interface for updating the sea-ice state due to ice ridging processes
!! using Icepack.
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
         ndtd       , & !< Thenumber of dynamics subcycles
         ncat  , & !< The number of thickness categories
         nilyr , & !< The number of ice layers
         nslyr , & !< The number of snow layers
         n_aero, & !< The number of aerosol tracers
         ntrcr     !< The number of tracers in use

      real (kind=dbl_kind), intent(in) :: &
         mu_rdg , & !< The e-folding scale of ridged ice [m^0.5]
         dt         !< The time step over which ridging occurs [s]

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max   !< category limits [m]

      real (kind=dbl_kind), intent(in) :: &
         rdg_conv   , & !< Normalized energy dissipation due to convergence [s-1]
         rdg_shear      !< Normalized energy dissipation due to shear [s-1]

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         aicen      , & !< The concentration of ice [nondim]
         vicen      , & !< The volume per unit area of ice [m]
         vsnon          !< The volume per unit area of snow [m]

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn          !< Ice tracer concentrations [kg m-3]

      real (kind=dbl_kind), intent(inout) :: &
         aice0          !< The concentration of open water [nondim]

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & !<  = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  !< The number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      !< = 0 or 1 depending on tracer dependency [nondim]
                        !! argument 2:  (1) aice, (2) vice, (3) vsno [nondim]

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      !< Tndices of underlying tracer layers

      integer (kind=int_kind), intent(in) :: &
         krdg_partic, & !< Selects participation function
         krdg_redist    !< Selects redistribution function

      logical (kind=log_kind), intent(in) :: &
         closing_flag, &!< flag if closing is valid
         tr_brine       !< if .true., brine height differs from ice thickness

      ! optional history fields
      real (kind=dbl_kind), intent(inout), optional :: &
         dardg1dt   , & !< The rate of fractional area loss by ridging ice [s-1]
         dardg2dt   , & !< The rate of fractional area gain by new ridges [s-1]
         dvirdgdt   , & !< The rate of ice volume ridged [m s-1]
         opening    , & !< The rate of opening due to divergence/shear [s-1]
         closing    , & !< The rate of closing due to divergence/shear [s-1]
         fpond      , & !< The fresh water flux to ponds [kg m2 s-1]
         fresh      , & !< The fresh water flux to ocean [kg m2 s-1]
         fhocn          !< The net heat flux to ocean [W m-2]

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         dardg1ndt  , & !< The rate of fractional area loss by ridging ice [s-1]
         dardg2ndt  , & !< The rate of fractional area gain by new ridges [s-1]
         dvirdgndt  , & !< The rate of ice volume ridged [m s-1]
         aparticn   , & !< The participation function [nondim]
         krdgn      , & !< The mean ridge thickness/thickness of ridging ice [m]
         araftn     , & !< The rafting ice area [m2]
         vraftn     , & !< The rafting ice volume [m3]
         aredistn   , & !< The redistribution function: fraction of new ridge area [nondim]
         vredistn       !< The redistribution function: fraction of new ridge volume [nondim]

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         faero_ocn      !< The aerosol flux to ocean [kg m-2 s-1]

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         fiso_ocn       !< isotope flux to ocean [kg m-2 s-1]




      end subroutine ridge_ice


      end module icepack_mechred

!=======================================================================
