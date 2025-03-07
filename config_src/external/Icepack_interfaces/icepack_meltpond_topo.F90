module icepack_meltpond_topo

  use icepack_kinds, only : int_kind, dbl_kind
  use icepack_tracers, only : n_iso, n_aero

  implicit none

  private
  public :: compute_ponds_topo
  contains

!> Interface for updating the melt ponds using Icepack.
  subroutine compute_ponds_topo(dt, ncat, nilyr,        &
                                    ktherm,             &
                                    aice,  aicen,       &
                                    vice,  vicen,       &
                                    vsno,  vsnon,       &
                                    meltt,              &
                                    fsurf, fpond,       &
                                    Tsfcn, Tf,          &
                                    qicen, sicen,       &
                                    apnd,  hpnd, ipnd   )

    integer (kind=int_kind), intent(in) :: &
         ncat , &   !< number of thickness categories
         nilyr, &   !< number of ice layers
         ktherm     !< type of thermodynamics (0 0-layer, 1 BL99, 2 mushy)

    real (kind=dbl_kind), intent(in) :: &
         dt !< time step (s)

    real (kind=dbl_kind), intent(in) :: &
         aice, &    !< total ice area fraction
         vsno, &    !< total snow volume (m)
         Tf   !< ocean freezing temperature [= ice bottom temperature] (degC)

    real (kind=dbl_kind), intent(inout) :: &
         vice, &    !< total ice volume (m)
         fpond      !< fresh water flux to ponds (m)

    real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen, &   !< ice area fraction, per category
         vsnon      !< snow volume, per category (m)

    real (kind=dbl_kind), dimension (:), intent(inout) :: &
         vicen      !< ice volume, per category (m)

    real (kind=dbl_kind), dimension (:), intent(in) :: &
         Tsfcn      !< surface temperature

    real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         qicen, &   !< ice enthalpy
         sicen      !< ice salinity

    real (kind=dbl_kind), dimension (:), intent(inout) :: &
         apnd, &    !< pond fractional area
         hpnd, &    !< pond depth (m)
         ipnd       !< pond ice thickness (m)

    real (kind=dbl_kind), intent(in) :: &
         meltt, &   !< total surface meltwater flux
         fsurf      !< thermodynamic heat flux at ice/snow surface (W/m^2)

  end subroutine compute_ponds_topo

end module icepack_meltpond_topo
