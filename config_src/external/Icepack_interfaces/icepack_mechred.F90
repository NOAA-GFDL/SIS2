      module icepack_mechred

      use icepack_kinds
      use icepack_tracers, only : n_iso, n_aero

      implicit none

      private
      public :: icepack_step_ridge
      contains

!-----------------------------------------------------------------------
!> Interface for updating the sea-ice state due to ice ridging processes
!! using Icepack.
      subroutine icepack_step_ridge(dt,           ndtd,          &
                                    nilyr,        nslyr,         &
                                    nblyr,                       &
                                    ncat,         hin_max,       &
                                    rdg_conv,     rdg_shear,     &
                                    aicen,                       &
                                    trcrn,                       &
                                    vicen,        vsnon,         &
                                    aice0,        trcr_depend,   &
                                    trcr_base,    n_trcr_strata, &
                                    nt_strata,                   &
                                    dardg1dt,     dardg2dt,      &
                                    dvirdgdt,     opening,       &
                                    fpond,                       &
                                    fresh,        fhocn,         &
                                    n_aero,                      &
                                    faero_ocn,    fiso_ocn,      &
                                    aparticn,     krdgn,         &
                                    aredistn,     vredistn,      &
                                    dardg1ndt,    dardg2ndt,     &
                                    dvirdgndt,                   &
                                    araftn,       vraftn,        &
                                    aice,         fsalt,         &
                                    first_ice,    fzsal,         &
                                    flux_bio,     closing,       &
                                    Tf)

      real (kind=dbl_kind), intent(in) :: &
         dt        !< The time step over which ridging occurs [s]

      integer (kind=int_kind), intent(in) :: &
         ncat  , & !< The number of thickness categories
         ndtd  , & !< Thenumber of dynamics subcycles
         nblyr , & !< The number of bio layers
         nilyr , & !< The number of ice layers
         nslyr , & !< The number of snow layers
         n_aero    !< The number of aerosol tracers

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max   !< category limits [m]

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & !< = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  !< The number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      !< = 0 or 1 depending on tracer dependency
                        !! argument 2:  (1) aice, (2) vice, (3) vsno [nondim]

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      !< Indices of underlying tracer layers

      real (kind=dbl_kind), intent(inout) :: &
         aice     , & !< sea ice concentration
         aice0    , & !< The concentration of open water [nondim]
         rdg_conv , & !< Normalized energy dissipation due to convergence [s-1]
         rdg_shear, & !< Normalized energy dissipation due to shear [s-1]
         dardg1dt , & !< rate of area loss by ridging ice [s-1]
         dardg2dt , & !< rate of area gain by new ridges [s-1]
         dvirdgdt , & !< rate of ice volume ridged [m s-1]
         opening  , & !< rate of opening due to divergence/shear [s-1]
         fpond    , & !< fresh water flux to ponds [kg m2 s-1]
         fresh    , & !< fresh water flux to ocean [kg m2 s-1]
         fsalt    , & !< salt flux to ocean [kg m2 s-1]
         fhocn        !< net heat flux to ocean [W m-2]

      real (kind=dbl_kind), intent(inout), optional :: &
         fzsal        !< zsalinity flux to ocean [kg m2 s-1] (deprecated)

      real (kind=dbl_kind), intent(inout), optional :: &
         closing      !< rate of closing due to divergence/shear [s-1]

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen    , & !< The concentration of ice [nondim]
         vicen    , & !< The volume per unit area of ice [m]
         vsnon    , & !< The volume per unit area of snow [m]
         dardg1ndt, & !< The rate of area loss by ridging ice [s-1]
         dardg2ndt, & !< The rate of area gain by new ridges [s-1]
         dvirdgndt, & !< The rate of ice volume ridged [m s-1]
         aparticn , & !< The participation function [nondim]
         krdgn    , & !< The mean ridge thickness/thickness of ridging ice [m]
         araftn   , & !< The rafting ice area [m2]
         vraftn   , & !< The rafting ice volume [m3]
         aredistn , & !< The redistribution function: fraction of new ridge area [nondim]
         vredistn , & !< The redistribution function: fraction of new ridge volume [nondim]
         faero_ocn, & !< The aerosol flux to ocean [kg m-2 s-1]
         flux_bio     !< All bio fluxes to ocean

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
         fiso_ocn     !< isotope flux to ocean [kg m-2 s-1]

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn        !< Ice tracer concentrations [kg m-3]

      !logical (kind=log_kind), intent(in) :: &
         !tr_pond_topo,& ! If .true., use explicit topography-based ponds
         !tr_aero     ,& ! If .true., use aerosol tracers
         !tr_brine    ,& ! If .true., brine height differs from ice thickness

      logical (kind=log_kind), dimension(:), intent(inout) :: &
         first_ice    !< True until ice forms
      real (kind=dbl_kind), intent(in) :: &
         Tf           ! freezing temperature

      end subroutine icepack_step_ridge

      end module icepack_mechred

!=======================================================================
