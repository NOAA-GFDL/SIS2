module icepack_meltpond_sealvl

  use icepack_kinds, only : int_kind, dbl_kind, char_len
  use icepack_tracers, only : n_iso, n_aero

  implicit none

  private
  public ::   compute_ponds_sealvl,   &
              pond_hypsometry,        &
              pond_height
  contains

!> Interface for updating the melt ponds using Icepack.
  subroutine compute_ponds_sealvl( dt,                    &
                                   meltt,  melts,  frain, &
                                   Tair,   fsurfn, Tsfcn, &
                                   dhs,    ffrac,         &
                                   aicen,  vicen,  vsnon, &
                                   qicen,  sicen,         &
                                   apnd,   hpnd,  ipnd,   &
                                   meltsliqn,             &
                                   dpnd_freebdn,          &
                                   dpnd_dlidn, dpnd_flushn)

  real(kind=dbl_kind), intent(in) :: &
       dt          !< time step (s)

  real(kind=dbl_kind), intent(in) :: &
       Tsfcn, &    !< surface temperature (C)
       meltt, &    !< top melt rate (m/s)
       melts, &    !< snow melt rate (m/s)
       frain, &    !< rainfall rate (kg/m2/s)
       Tair,  &    !< air temperature (K)
       fsurfn,&    !< atm-ice surface heat flux  (W/m2)
       aicen, &    !< ice area fraction
       vicen, &    !< ice volume (m)
       vsnon, &    !< snow volume (m)
       meltsliqn   !< liquid contribution to meltponds in dt (kg/m^2)

  real(kind=dbl_kind), intent(inout) :: &
       apnd, hpnd, ipnd, & !< pond tracers
       dpnd_freebdn,     & !< pond drainage rate due to freeboard constraint (m/step)
       dpnd_dlidn,       & !< pond loss/gain due to ice lid (m/step)
       dpnd_flushn         !< pond flushing rate due to ice permeability (m/s)

  real(kind=dbl_kind), dimension (:), intent(in) :: &
       qicen, &    !< ice layer enthalpy (J m-3)
       sicen       !< salinity (ppt)

  real(kind=dbl_kind), intent(in) :: &
       dhs         !< depth difference for snow on sea ice and pond ice

  real(kind=dbl_kind), intent(out) :: &
       ffrac       !< fraction of fsurfn over pond used to melt ipond

  end subroutine compute_ponds_sealvl

  subroutine pond_hypsometry(hpnd, apnd, dhpond, dvpond, hin)

  real(kind=dbl_kind), intent(inout) :: &
       hpnd, &   !< pond depth of ponded area tracer [m]
       apnd      !< pond fractional area of category tracer

  real(kind=dbl_kind), intent(in), optional :: &
       dvpond, & !< incoming change in pond volume per category area
       dhpond, & !< incoming change in pond depth [m]
       hin       !< category ice thickness [m]

  end subroutine pond_hypsometry

  subroutine pond_height(apond, hpnd, hin, hpsurf)

  real(kind=dbl_kind), intent(in) :: &
       hin, &   !< category mean ice thickness [m]
       apond, & !< pond area fraction of the category
       hpnd     !< mean pond depth [m]

  real(kind=dbl_kind), intent(out) :: &
       hpsurf   !< height of pond surface above base of the ice [m]

  end subroutine pond_height

end module icepack_meltpond_sealvl
