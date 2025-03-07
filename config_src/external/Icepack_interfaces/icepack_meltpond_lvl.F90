module icepack_meltpond_lvl

  use icepack_kinds, only : int_kind, dbl_kind, char_len
  use icepack_tracers, only : n_iso, n_aero

  implicit none

  private
  public :: compute_ponds_lvl
  contains

!> Interface for updating the melt ponds using Icepack.
  subroutine compute_ponds_lvl(dt,     nilyr,        &
                                   ktherm,               &
                                   hi_min, dpscale,      &
                                   frzpnd,               &
                                   rfrac,  meltt, melts, &
                                   frain,  Tair,  fsurfn,&
                                   dhs,    ffrac,        &
                                   aicen,  vicen, vsnon, &
                                   qicen,  sicen,        &
                                   Tsfcn,  alvl,         &
                                   apnd,   hpnd,  ipnd,  &
                                   meltsliqn)

    integer (kind=int_kind), intent(in) :: &
         nilyr, &    !< number of ice layers
         ktherm      !< type of thermodynamics (-1 none, 1 BL99, 2 mushy)

    real (kind=dbl_kind), intent(in) :: &
         dt,       & !< time step (s)
         hi_min,   & !< minimum ice thickness allowed for thermo (m)
         dpscale     !< alter e-folding time scale for flushing

    character (len=char_len), intent(in) :: &
         frzpnd      !< pond refreezing parameterization

    real (kind=dbl_kind), intent(in) :: &
         Tsfcn, &    !< surface temperature (C)
         alvl,  &    !< fraction of level ice
         rfrac, &    !< water fraction retained for melt ponds
         meltt, &    !< top melt rate (m/s)
         melts, &    !< snow melt rate (m/s)
         frain, &    !< rainfall rate (kg/m2/s)
         Tair,  &    !< air temperature (K)
         fsurfn,&    !< atm-ice surface heat flux  (W/m2)
         aicen, &    !< ice area fraction
         vicen, &    !< ice volume (m)
         vsnon, &    !< snow volume (m)
         meltsliqn   !< liquid contribution to meltponds in dt (kg/m^2)

    real (kind=dbl_kind), intent(inout) :: &
         apnd, &     !< pond fraction
         hpnd, &     !< pond depth (m)
         ipnd        !< ice pond thickness (m)

    real (kind=dbl_kind), dimension (:), intent(in) :: &
         qicen, &  !< ice layer enthalpy (J m-3)
         sicen     !< salinity (ppt)

    real (kind=dbl_kind), intent(in) :: &
         dhs       !< depth difference for snow on sea ice and pond ice (m)

    real (kind=dbl_kind), intent(out) :: &
         ffrac     !< fraction of fsurfn over pond used to melt ipond

  end subroutine compute_ponds_lvl

end module icepack_meltpond_lvl
