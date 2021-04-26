!> Routines to do SIS2 ice thermodynamic calcualtions
module SIS2_ice_thm

! This file is part of SIS2. See LICENSE.md for the license.

use ice_thm_mod,         only : get_thermo_coefs
use MOM_EOS,             only : EOS_type, EOS_init, EOS_end
use MOM_error_handler,   only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,     only : get_param, log_param, read_param, log_version, param_file_type
use MOM_obsolete_params, only : obsolete_logical, obsolete_real
use MOM_unit_scaling,    only : unit_scale_type

implicit none ; private

public :: SIS2_ice_thm_init, SIS2_ice_thm_end, ice_temp_SIS2, estimate_tsurf
public :: ice_resize_SIS2, add_frazil_SIS2, rebalance_ice_layers

public :: get_SIS2_thermo_coefs, ice_thermo_init, ice_thermo_end
public :: Temp_from_Enth_S, Temp_from_En_S, enth_from_TS, enthalpy_from_TS
public :: enthalpy_liquid_freeze, T_Freeze, calculate_T_Freeze, enthalpy_liquid
public :: e_to_melt_TS, energy_melt_enthS, latent_sublimation

!> This type contains the parameters regulating sea-ice thermodyanmics
type, public :: ice_thermo_type ; private
  real :: Cp_ice            !< The heat capacity of ice [Q degC-1 ~> J kg-1 degC-1].
  real :: Cp_water          !< The heat capacity of liquid water in the ice model,
                            !! but not in the brine pockets [Q degC-1 ~> J kg-1 degC-1].
  real :: Cp_brine          !< The heat capacity of liquid water in the brine
                            !! pockets within the ice [Q degC-1 ~> J kg-1 degC-1].  Cp_brine
                            !! should be set equal to Cp_Water, but for
                            !! algorithmic convenience can be set equal to Cp_ice.
  real :: rho_ice           !< The nominal density of ice [R ~> kg m-3].
  real :: rho_snow          !< The nominal density of snow [R ~> kg m-3].
  real :: rho_water         !< The nominal density of water [R ~> kg m-3].
  real :: Lat_Fus           !< The latent heat of fusion [Q ~> J kg-1].
  real :: Lat_Vapor         !< The latent heat of vaporization [Q ~> J kg-1].
  real :: dTf_dS            !< The derivative of the freezing point with salinity,
                            !! [degC kg gSalt-1].  (dTf_dS is negative.)

  real :: enth_liq_0 = 0.0  !< The value of enthalpy for liquid fresh water at 0 degC [Q ~> J kg-1].
  real :: Q_to_J_kg = 1.0   !< A conversion factor for enthalpy back to Joules kg-1 [J kg-1 Q-1 ~> 1].
  logical :: slab_ice = .false. !< If true use the very old slab ice thermodynamics,
                            !! with effectively zero heat capacity of ice and snow.
  logical :: sublimation_bug = .false.  !< If true use an older calculation that omits the
                            !! latent heat of fusion from the latent heat of sublimation.
                            !! This variable should be obsoleted.
  type(EOS_type), pointer :: EOS=>NULL() !< A pointer to the shared MOM6/SIS2
                            !! equation-of-state type. This is here to encourage the use of
                            !! common and consistent thermodynamics between the ice and ocean.
end type ice_thermo_type


!> The control structure for the SIS2 ice thermodynamics
type, public :: SIS2_ice_thm_CS ; private
  ! properties of ice, snow, and seawater (NCAR CSM values)
  real :: KS   !< Thermal conductivity of snow, often 0.31 [Q R Z2 T-1 ~> W m-1 degC-1]
  real :: KI   !< Thermal conductivity of ice, often 2.03 [Q R Z2 T-1 ~> W m-1 degC-1]

  real :: temp_ice_freeze !< The freezing temperature of the top ice layer [degC].
  real :: temp_range_est  !< An estimate of the range of snow and ice temperatures
                          !! that is used to evaluate whether an explicit
                          !! diffusive form of the heat fluxes or an inversion
                          !! based on the layer heat budget is more likely to
                          !! be the most accurate [degC].

  real :: H_lo_lim        !< The lower limit on snow or ice thickness for temperature calculations [Z ~> m]
  real :: frazil_temp_offset = 0.5 !< An offset between the temperature with
                          !! which frazil ice forms and the freezing point of
                          !! each sublayer of the ice [degC].  This functionality could
                          !! later be accounted for using liq_lim instead.

  real :: liq_lim = .99   !< A limit on the liquid fraction of sea-ice before it can no longer
                          !! be considered sea-ice, as it will simply contract due to the rising
                          !! ice shards [nondim].
  ! In the ice temperature calculation we place a limit to below (salinity
  ! dependent) freezing point on the prognosed temperatures.  For ice_resize
  ! it is better to make a slightly more restrictive limit that requires the
  ! temperature to be such that the brine content is less than "liq_lim" of
  ! the total mass.  That is T_f/T < liq_lim implying T<T_f/liq_lim

  logical :: do_pond = .false. !< If true, activate melt pond scheme
  ! mw/new - these melt pond control data are temporarily placed here
  real :: tdrain = -0.8 !< If average ice temp. > tdrain, drain pond [degC]
  real :: r_min_pond = 0.15 !< Pond retention of meltwater
  real :: r_max_pond = 0.9  !< See CICE5 doc
  real :: max_pond_frac = 0.5  !< Pond water beyond this is dumped
  real :: min_pond_frac = 0.2  !< Ponds below sea level don't drain
  ! mw/new - end of melt pond control data
end type SIS2_ice_thm_CS


contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS2_ice_thm_init initializes the control structure for the ice thermodynamic
!! update code.
subroutine SIS2_ice_thm_init(US, param_file, CS)

  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors
  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS2_ice_thm_CS), pointer    :: CS !< The SIS2 ice thermodynamics control structure

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "SIS2_ice_thm (updates)" ! This module's name.

  if (.not.associated(CS)) allocate(CS)

  call log_version(param_file, mdl, version, &
     "This sub-module does updates of the sea-ice due to thermodynamic changes.")

  call get_param(param_file, mdl, "SNOW_CONDUCTIVITY", CS%Ks, &
                 "The conductivity of heat in snow.", &
                 units="W m-1 K-1", default=0.31, scale=US%W_m2_to_QRZ_T*US%m_to_Z)
  call get_param(param_file, mdl, "ICE_CONDUCTIVITY", CS%Ki, &
                 "The conductivity of heat in ice.", &
                 units="W m-1 K-1", default=2.03, scale=US%W_m2_to_QRZ_T*US%m_to_Z)
  call get_param(param_file, mdl, "MIN_H_FOR_TEMP_CALC", CS%H_lo_lim, &
                 "The minimum ice thickness at which to do temperature "//&
                 "calculations.", units="m", default=0.0, scale=US%m_to_Z)

  call get_param(param_file, mdl, "DO_POND", CS%do_pond, &
                 "If true, calculate melt ponds and use them for "//&
                 "shortwave radiation calculation.", default=.false.)
  call get_param(param_file, mdl, "TDRAIN", CS%tdrain, &
                 "Melt ponds drain to sea level when ice average temp. "//&
                 "exceeds TDRAIN (stand-in for mushy layer thermo)", default=-0.8, units="degC")
  call get_param(param_file, mdl, "R_MIN_POND", CS%r_min_pond, &
                 "Minimum retention rate of surface water sources in melt pond "//&
                 "(retention scales linearly with ice cover)", default=0.15, units="nondim")
  call get_param(param_file, mdl, "R_MAX_POND", CS%r_max_pond, &
                 "Maximum retention rate of surface water sources in melt pond "//&
                 "(retention scales linearly with ice cover)", default=0.9, units="nondim")
  call get_param(param_file, mdl, "MIN_POND_FRAC", CS%min_pond_frac, &
                 "Minimum melt pond cover (by ponds at sea level) "//&
                 "pond drains to this when ice is porous.", default=0.2, units="nondim")
  call get_param(param_file, mdl, "MAX_POND_FRAC", CS%max_pond_frac, &
                 "Maximum melt pond cover - associated with pond volume "//&
                 "that suppresses ice top to waterline", default=0.5, units="nondim")
  call get_param(param_file, mdl, "ICE_TEMP_RANGE_ESTIMATE", CS%temp_range_est,&
                 "An estimate of the range of snow and ice temperatures "//&
                 "that is used to evaluate whether an explicit diffusive "//&
                 "form of the heat fluxes or an inversion based on the "//&
                 "layer heat budget is more likely to be more accurate. "//&
                 "Setting this to 0 causes the explicit diffusive form. "//&
                 "to always be used.", units="degC", default=40.0)
  CS%temp_range_est = abs(CS%temp_range_est)

  call obsolete_logical(param_file, "OLD_ICE_HEAT_CAPACITY", warning_val=.false.)

end subroutine SIS2_ice_thm_init



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_temp_SIS2 calculates the updated snow and ice enthalpy and new skin
!!    temperature after due to surface forcing and vertical diffusion of heat.
subroutine ice_temp_SIS2(m_pond, m_snow, m_ice, enthalpy, sice, SF_0, dSF_dT, sol, tfw, fb, &
                         tsurf, dtt, NkIce, tmelt, bmelt, CS, US, ITV, check_conserve)

  real, intent(in   ) :: m_pond  !< pond mass per unit area [R Z ~> kg m-2]
  real, intent(in   ) :: m_snow  !< snow mass per unit area [R Z ~> kg m-2]
  real, intent(in   ) :: m_ice   !< ice mass per unit area [R Z ~> kg m-2]
  real, dimension(0:NkIce) , &
        intent(inout) :: enthalpy !< The enthalpy of each layer in a column of
                                  !! snow and ice [Q ~> J kg-1].
  real, dimension(NkIce), &
        intent(in)    :: Sice  !< ice salinity by layer [gSalt kg-1]
  real, intent(in   ) :: SF_0  !< net upward surface heat flux at ts=0 [Q R Z T-1 ~> W m-2]
  real, intent(in   ) :: dSF_dT !< d(sfc heat flux)/d(ts) [Q R Z T-1 degC-1 ~> W m-2 degC-1]
  real, dimension(0:NkIce), &
        intent(in   ) :: sol   !< Solar heating of the snow and ice layers [Q R Z T-1 ~> W m-2]
  real, intent(in   ) :: tfw   !< seawater freezing temperature [degC]
  real, intent(in   ) :: fb    !< heat flux upward from ocean to ice bottom [Q R Z T-1 ~> W m-2]
  real, intent(  out) :: tsurf !< surface temperature [degC]
  real, intent(in   ) :: dtt   !< timestep [T ~> s]
  integer, intent(in   ) :: NkIce !< The number of ice layers.
  real, intent(inout) :: tmelt !< accumulated top melting energy  [Q R Z ~> J m-2]
  real, intent(inout) :: bmelt !< accumulated bottom melting energy [Q R Z ~>  J m-2]
  type(SIS2_ice_thm_CS), intent(in) :: CS  !< The SIS2 ice thermodynamics control structure
  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  logical, optional, intent(in) :: check_conserve !< If true, check for local heat conservation.

  ! Local variables for temperature calculation [see Winton (1999) section II.A.]
  ! note:  here equations are multiplied by hi to improve thin ice accuracy
  real, dimension(0:NkIce) :: temp_est   ! An estimated snow and ice temperature [degC].
  real, dimension(0:NkIce) :: temp_IC    ! The temperatures of the snow and ice based on the initial
                                         ! enthalpy [degC].
  real, dimension(0:NkIce) :: temp_new   ! The updated temperatures [degC].
  real, dimension(NkIce) :: tfi  ! The ice freezing temperatures [degC].
  real :: mL_ice   ! The mass-per-unit-area of each ice layer [R Z ~> kg m-2]
  real :: mL_snow  ! The mass-per-unit-area of each snow layer [R Z ~> kg m-2]
  real :: e_extra  ! [Q R Z ~> J m-2]
  real, dimension(0:NkIce) :: m_lay     ! Masses of all layers [R Z ~> kg m-2].
  real :: enth_fp  ! The enthalpy at the freezing point (solid for fresh ice).
  real :: kk        ! The conductive thermal coupling coefficient between adjacent
                    ! ice sublayers [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: k10       ! The conductive thermal coupling coefficient between the
                    ! snow and the topmost ice sublayer [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: k0a       ! The implicit conductive thermal coupling coefficient between
                    ! the snow and the skin temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: k0skin    ! The conductive thermal coupling coefficient between the
                    ! snow and the skin temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: I_bb      ! [degC Q-1 R-1 Z-1 ~> degC m2 J-1].
  real :: b_denom_1 ! [Q R Z degC-1 ~> J m-2 degC-1].
  real :: comp_rat  ! The complement of cc_bb, going from 0 to 1 [nondim].
  real :: tsf       ! The surface freezing temperature [degC].
  real :: k0a_x_ta  ! K0a times the atmospheric temperature [Q R Z T-1 ~> W m-2]
  ! real :: tsno_est, salt_part ! , rat
  real :: tsurf_est  ! An estimate of the surface temperature [degC].
  real, dimension(0:NkIce+1) :: cc ! Interfacial coupling coefficients [Q R Z degC-1 ~> J m-2 degC-1].
  real, dimension(0:NkIce) :: bb   ! Effective layer heat capacities [Q R Z degC-1 ~> J m-2 degC-1].
  real, dimension(0:NkIce) :: cc_bb ! Remaining coupling ratios [nondim].
  real, dimension(-1:NkIce) :: heat_flux_int ! The downward heat fluxes at the
                                ! interfaces between layers [Q R Z T-1 ~> W m-2].
                                ! heat_flux_int uses the index convention from
                                ! MOM6 that interface K is below layer k.
  real :: I_liq_lim     ! The inverse of CS%liq_lim [nondim].
  real :: heat_flux_err_rat ! A factor used to estimate which of two approaches is more
                            ! accurate [Q T-1 degC-1 ~> W kg-1 degC-1]
  real :: col_enth1, col_enth2, col_enth2b, col_enth3  ! [Q R Z ~> J m-2]
  real :: d_e_extra   ! [Q R Z ~> J m-2]
  real :: e_extra_sum ! [Q R Z ~> J m-2]
  real :: tflux_bot   ! Heat flux in the ice at the ice-ocean interface [Q R Z ~> J m-2]
  real :: tflux_sfc   ! [Q R Z ~> J m-2]
  ! These are used in diagnostic calculations that could be uncommented for debugging.
  ! real :: sum_sol     ! The time integrated shortwave heating [Q R Z ~> J m-2]
  ! real :: tfb_diff_err  ! A diagostic estimate of errors from roundoff [Q R Z ~> J m-2]
  ! real :: tfb_resid_err ! A diagnostic estimate of the residual due to roundoff [Q R Z ~> J m-2]
  ! real :: d_tflux_bot, tflux_bot_diff, tflux_bot_resid ! Diagnostic terms [Q R Z ~> J m-2]
  real :: hsnow_eff    ! The effective thickness of the snow layer [Z ~> m]
  real :: snow_temp_new, snow_temp_max ! Snow temperatures [degC]
  real :: hL_ice_eff   ! The effective thickness of each ice layer [Z ~> m]
  real :: enth_liq_lim ! The enthalpy at the point where the ice or snow stops
                       ! acting as a solid, and all extra heat goes into
                       ! melting [Q ~> J kg-1].
  real :: rho_ice  ! The nominal density of sea ice [R ~> kg m-3].
  real :: rho_snow ! The nominal density of snow [R ~> kg m-3].
  real :: Cp_ice   ! The heat capacity of ice [Q degC-1 ~> J kg-1 degC-1].
  real :: Cp_brine ! The heat capacity of liquid water in the brine pockets,
                   ! [Q degC-1 ~> J kg-1 degC-1].
  real :: Lat_fus  ! The latent heat of fusion [Q ~> J kg-1].
  logical :: col_check ! If true, check for local heat conservation.
  integer :: k

  col_check = .false. ; if (present(check_conserve)) col_check = check_conserve

  temp_IC(0) = temp_from_En_S(enthalpy(0), 0.0, ITV)
  call temp_from_Enth_S(enthalpy(1:), sice(1:), temp_IC(1:), ITV)

  call get_SIS2_thermo_coefs(ITV, rho_ice=rho_ice, rho_snow=rho_snow, &
                             Cp_Ice=Cp_ice, Latent_Fusion=Lat_fus, Cp_Brine=Cp_Brine)

  mL_ice = m_ice / NkIce   ! ice mass per unit area of each layer
  mL_snow = m_snow         ! snow mass per unit area [R Z ~> kg m-2].
  call calculate_T_Freeze(sice, tfi, ITV)    ! freezing temperature of ice layers

  ! Set the effective thickness of each ice and snow layer, limited to avoid
  ! instabilities for thin layers.
  hL_ice_eff = max(mL_ice / rho_ice, CS%H_lo_lim)
  hsnow_eff = mL_snow / rho_snow + max(1.0e-35*US%m_to_Z, 1.0e-20*CS%H_lo_lim)

  ! Set the surface freezing temperature
  tsf = tfi(1) ; if (mL_snow>0.0) tsf = 0.0

  ! Set some diffusive coupling coefficients between the various layers.
  kk = CS%Ki / hL_ice_eff     ! full ice layer conductivity
  k10 = 2.0*(CS%Ks*CS%Ki) / (hL_ice_eff*CS%Ks + hsnow_eff*CS%Ki) ! coupling ice layer 1 to snow
  k0a = (CS%Ks*dSF_dT) / (0.5*dSF_dT*hsnow_eff + CS%Ks)      ! coupling snow to "air"
  k0skin = 2.0*CS%Ks / hsnow_eff
  k0a_x_ta = (CS%Ks*SF_0) / (0.5*dSF_dT*hsnow_eff + CS%Ks) ! coupling times "air" temperature

  ! Determine the enthalpy for conservation checks.
  m_lay(0) = mL_snow ; do k=1,NkIce ; m_lay(k) = mL_ice ; enddo

  if (col_check) then
    col_enth1 = 0.0
    do k=0,NkIce ; col_enth1 = col_enth1 + m_lay(k)*enthalpy(k) ; enddo
  endif
  e_extra_sum = 0.0

  !   First get a non-conservative estimate with a linearized fully implicit
  ! treatment of layer coupling.

  ! Determine the effective layer heat capacities.
  !   bb = dheat/dTemp, as derived from a linearization of the enthalpy equation
  !     but using the value for ice near the melt point if temp_IC is above T_fr.
  bb(0) = mL_snow* Cp_ice
  do k=1,NkIce   ! Store the heat capacity term in bb.
    if (tfi(k) >= 0.0) then ! This is pure ice.
      bb(k) = mL_ice * Cp_ice
    elseif (temp_IC(k) < tfi(k)) then
      bb(k) = mL_ice * (Cp_ice - (tfi(k) / temp_IC(k)**2) * &
                        (Lat_Fus - (Cp_brine-Cp_ice) * temp_IC(k)) )
      ! Or more generally:
      !  S_Sf = sice(k) / calculate_S_freeze(temp_IC(k))
      !  bb(k) = mL_ice * (Cp_ice + dSf_dT*Lat_fus + S_Sf * (Cp_brine-Cp_ice))
    else ! Use the value when temp_IC = tfi.
      bb(k) = mL_ice * (Cp_brine - Lat_fus / tfi(k))
    endif
  enddo

  ! The following expressions could be modified to permit there to be
  ! variable diffusivities in the ice/brine mixtures.
  cc(0) = k0a*dtt  ! Atmosphere-snow coupling
  cc(1) = k10*dtt  ! Snow-ice coupling
  do k=2,NkIce ; cc(k) = kk*dtt ; enddo
  cc(NkIce+1) = 2*kk*dtt ! bottom of ice coupling with ocean at freezing point.

  !   This is a version of a tridiagonal solver where the relationship between
  ! the diagonal elements and the coupling between layers has been used to
  ! ensure that there are no differences ever taken, so there will be minimal
  ! truncation errors.  This closely follows schemes that have previously been
  ! found to work very well in GOLD and MOM6.

  ! Go UP the ice column.
  b_denom_1 = (bb(NkIce) + cc(NkIce+1))
  I_bb = 1.0 / (b_denom_1 + cc(NkIce))
  temp_est(NkIce) = ( (sol(NkIce)*dtt + bb(NkIce)*temp_IC(NkIce)) + cc(NkIce+1)*tfw ) * I_bb
  comp_rat = b_denom_1 * I_bb
  cc_bb(NkIce) = cc(NkIce) * I_bb

  do k=NkIce-1,1,-1
    b_denom_1 = bb(k) + comp_rat*cc(k+1)
    I_bb =  1.0 / (b_denom_1 + cc(k))
    temp_est(k) = ((sol(k)*dtt + bb(k)*temp_IC(k)) + cc(k+1)*temp_est(k+1)) * I_bb

    comp_rat = b_denom_1 * I_bb ! 1.0 >= comp_rat >= 0.0
    cc_bb(k) = cc(k) * I_bb
  end do

  b_denom_1 = bb(0) + comp_rat*cc(1)
  I_bb =  1.0 / (b_denom_1 + cc(0))
  !   This is a complete calculation of temp_est(0), assuming that the upward surface
  ! flux is given by SF_0 + dSF_dT*tsurf, with tsurf estimated as part of this calculation.
  temp_est(0) = (((sol(0)*dtt + bb(0)*temp_IC(0)) - k0a_x_ta*dtt) + cc(1)*temp_est(1)) * I_bb

  ! Diagnose the surface skin temperature by matching the diffusive fluxes in
  ! the snow with the atmospheric fluxes.  I.e. solve the following for tsurf_est:
  !  (SF_0 + dSF_dT*tsurf_est) = k0skin * (temp_est(0) - tsurf_est)
  tsurf_est = (k0skin*temp_est(0) - SF_0) / (dSF_dT + k0skin)

  if (tsurf_est > tsf .or. m_pond > 0.0 ) then ! mw/new - liq. h2o @ sfc
                                               ! also pins temp at freezing
    ! The surface is melting, set tsurf to melt temp. and recalculate I_bb.
    tsurf_est = tsf
    ! cc(0) = k0skin*dtt

    I_bb =  1.0 / (b_denom_1 + k0skin*dtt)
    temp_est(0) = min(tsf, &
      (((sol(0)*dtt + bb(0)*temp_IC(0)) + k0skin*dtt*tsf) + cc(1)*temp_est(1)) * I_bb)
  endif
  ! Go back DOWN the ice column to get the estimated temperatures, subject to the
  ! limitation that all temperatures are assumed to be at or below freezing.
  do k=1,NkIce
    temp_est(k) = min(temp_est(k) + cc_bb(k) * temp_est(k-1), tfi(k))
  end do

  !
  ! Quasi-conservative iterative pass going UP the ice column
  !
  if (NkIce == 1) then
    temp_est(1) = laytemp_SIS2(mL_ice, tfi(1), sol(1) + (2*kk*tfw + k10*temp_est(0)), &
                               2*kk + k10, temp_IC(1), enthalpy(1), sice(1), dtt, ITV, US)
  else
    temp_est(NkIce) = laytemp_SIS2(mL_ice, tfi(NkIce), sol(NkIce) + kk*(2*tfw+temp_est(NkIce-1)), &
                                   3*kk, temp_IC(NkIce), enthalpy(NkIce), sice(NkIce), dtt, ITV, US)
    do k=NkIce-1,2,-1
      temp_est(k) = laytemp_SIS2(mL_ice, tfi(k), sol(k) + kk*(temp_est(k-1)+temp_est(k+1)), &
                                 2*kk, temp_IC(k), enthalpy(k), sice(k), dtt, ITV, US)
    enddo
    temp_est(1) = laytemp_SIS2(mL_ice, tfi(1), sol(1) + (kk*temp_est(2) + k10*temp_est(0)), &
                               kk + k10, temp_IC(1), enthalpy(1), sice(1), dtt, ITV, US)
  endif

  ! Calculate the bulk snow temperature and surface skin temperature together.
  temp_est(0) = laytemp_SIS2(mL_snow, 0.0, sol(0) + (k10*temp_est(1)-k0a_x_ta), &
                             k10+k0a, temp_IC(0), enthalpy(0), 0.0, dtt, ITV, US)
  tsurf = (k0skin*temp_est(0) - SF_0) / (dSF_dT + k0skin)  ! diagnose surface skin temp.

  !
  !   The following conservative update pass going DOWN the ice column is where
  ! the layer enthalpies are actually updated and extra heat that should drive
  ! melting is accumulated and stored.
  !

  !   Pre-calculate a conversion factor that can be used to figure out whether it
  ! is more accurate to use the explicit form for the fluxes or to invert enthalpy
  ! conservation for the fluxes, depending on the relative layer masses and the
  ! strength of the coupling between layers.  The expression below is a simplified
  ! form that assume that temperatures are measured in Celsius and are less than
  ! about CS%temp_err_est, and that the heat budget has 4 terms of comparable
  ! magnitude to the enthalpy content of the layer.  If there were a larger
  ! tolerance for the iterative estimate of a new temperature, that would go into
  ! the numerator but not the denominator.  The present form is based on the
  ! assumption that floating-point roundoff dominates the errors.
  heat_flux_err_rat = 0.7071 * dtt * (CS%temp_range_est) / &
                      (CS%temp_range_est * Cp_ice + Lat_fus)

  e_extra = 0.0
  if (tsurf > tsf .or. m_pond > 0.0) then
                  ! The surface is melting: update enthalpy with the surface at melt temp.
                  ! mw/new - liq h2o at surface also pins temp at freezing
    tsurf = tsf
    ! Accumulate surface melt energy.
    if (mL_snow>0.0) then
      heat_flux_int(-1) = k0skin * tsf
      heat_flux_int(0) = -k10*temp_est(1)
      call update_lay_enth(mL_snow, 0.0, enthalpy(0), heat_flux_int(-1), &
                           sol(0), heat_flux_int(0), -k0skin, k10, dtt, &
                           heat_flux_err_rat, ITV, US, e_extra)

      tmelt = tmelt + e_extra - dtt*((SF_0 + dSF_dT*tsf) + heat_flux_int(-1))
      tflux_sfc = dtt*heat_flux_int(-1)
      e_extra_sum = e_extra_sum + e_extra
    else
      ! There is no snow mass, so just convert tsnow = tsf to enthalpy.
      enthalpy(0) = enth_from_TS(tsf, 0.0, ITV)

      heat_flux_int(0) = k10*(tsf - temp_est(1))
      heat_flux_int(-1) = heat_flux_int(0)

      ! Replace use tsurf in the calculation of tmelt
      tmelt = tmelt + dtt*((sol(0) - (SF_0 + dSF_dT*tsf)) - heat_flux_int(0))
      tflux_sfc = dtt*heat_flux_int(0)

    endif
  else
    heat_flux_int(-1) = -k0a_x_ta
    heat_flux_int(0) = -k10*temp_est(1)
    snow_temp_max = (tsf*(dSF_dT + k0skin) + SF_0) / k0skin
    call update_lay_enth(mL_snow, 0.0, enthalpy(0), heat_flux_int(-1), &
                         sol(0), heat_flux_int(0), -k0a, k10, dtt, &
                         heat_flux_err_rat, ITV, US, e_extra, &
                         temp_new=snow_temp_new, temp_max=snow_temp_max)
    tsurf = (k0skin*snow_temp_new - SF_0) / (dSF_dT + k0skin)  ! diagnose surface skin temp.
    ! This is equivalent to, but safer than, tsurf = -(SF_0 + heat_flux_int(-1)) / dSF_dT

    tflux_sfc = dtt*heat_flux_int(-1)
    e_extra_sum = e_extra_sum + e_extra
    tmelt = tmelt + e_extra
  endif

  do k=1,NkIce-1 ! The flux from above is fixed, only have downward feedback
    heat_flux_int(K) = -kk*temp_est(k+1)
    call update_lay_enth(mL_ice, Sice(k), enthalpy(k), heat_flux_int(K-1), &
                         sol(k), heat_flux_int(K), 0.0, kk, dtt, &
                         heat_flux_err_rat, ITV, US, e_extra)

    e_extra_sum = e_extra_sum + e_extra
    if (k <= NkIce/2) then ; tmelt = tmelt + e_extra
    else ; bmelt = bmelt + e_extra ; endif
  enddo

  heat_flux_int(NkIce) = -2.0*kk*tfw
  call update_lay_enth(mL_ice, Sice(NkIce), enthalpy(NkIce), heat_flux_int(NkIce-1), &
                       sol(NkIce), heat_flux_int(NkIce), 0.0, 2.0*kk, dtt, &
                       heat_flux_err_rat, ITV, US, e_extra)
  e_extra_sum = e_extra_sum + e_extra
  bmelt = bmelt + e_extra
  !
  ! END of the conservative update of enthalpy.
  !

  if (col_check) then  ! Extra diagnostic expressions are commented out here.
    col_enth2b = 0.0 ! ; col_enth2 = e_extra_sum ; sum_sol = 0.0
    do k=0,NkIce
      ! col_enth2 = col_enth2 + m_lay(k)*enthalpy(k)
      col_enth2b = col_enth2b + m_lay(k)*enthalpy(k)
      ! sum_sol = sum_sol + dtt*sol(k)
    enddo
    ! ! tflux_bot_resid and tflux_bot_diff are two mathematically equivalent
    ! ! estimates of the heat flux at the base of the ice.
    ! tflux_bot_resid = (col_enth2 - col_enth1) - (sum_sol + tflux_sfc)
    ! tflux_bot_diff = -heat_flux_int(NkIce)*dtt

    ! ! Estimate the errors with these two expressions from 64-bit roundoff.
    ! tfb_diff_err = 1e-15*2.0*kk*dtt * sqrt(tfw**2 + 10.0**2)  ! The -10 degC is arbitrary but good enough?
    ! tfb_resid_err = 1e-15*sqrt(col_enth2**2 + col_enth1**2 + sum_sol**2 + tflux_sfc**2)
    ! d_tflux_bot = tflux_bot_diff - tflux_bot_resid
  endif

  tflux_bot = -heat_flux_int(NkIce)*dtt

  !   Accumulate the difference between the ocean's heat flux to the ice-ocean
  ! interface and the sea-ice heat flux to evaluate bottom melting/freezing.
  bmelt = bmelt + (dtt*fb - tflux_bot)

  e_extra_sum = 0.0

  enth_liq_lim = Enth_from_TS(0.0, 0.0, ITV)
  if (enthalpy(0) > enth_liq_lim) then ! put excess snow energy into top melt.
    e_extra = (enthalpy(0) - enth_liq_lim) * mL_snow
    tmelt = tmelt + e_extra
    e_extra_sum = e_extra_sum + e_extra
    enthalpy(0) = enth_liq_lim
  endif

  I_liq_lim = 1.0 / CS%liq_lim
  do k=1,NkIce
    enth_liq_lim = enth_from_TS(tfi(k)*I_liq_lim, sice(k), ITV)
    if (enthalpy(k) > enth_liq_lim) then
      ! Put excess energy into the closer of top or bottom melt.
      e_extra = (enthalpy(k) - enth_liq_lim) * mL_ice
      e_extra_sum = e_extra_sum + e_extra
      enthalpy(k) = enth_liq_lim
      if (k<=NkIce/2) then
        tmelt = tmelt + e_extra
      else
        bmelt = bmelt + e_extra
      endif
    endif
  enddo

  if (col_check) then
    col_enth3 = 0.0
    do k=0,NkIce ; col_enth3 = col_enth3 + m_lay(k)*enthalpy(k) ; enddo

    d_e_extra = col_enth3 - (col_enth2b - e_extra_sum)
    if (abs(d_e_extra) > 1.0e-12*(abs(col_enth3) + abs(col_enth2b) + &
                                  abs(e_extra_sum))) then
      d_e_extra = (col_enth3 - col_enth2b) - e_extra_sum
    endif
  endif

  call ice_check(mL_snow, NkIce*mL_ice, enthalpy, sice, NkIce, &
           "at end of ice_temp_SIS2", ITV, bmelt=bmelt, tmelt=tmelt, t_sfc=tsurf, US=US)

end subroutine ice_temp_SIS2


!### I think that estimate_tsurf is never used.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> estimate_tsurf estimates the surface skin temperature using the first (nonconservative)
!!    part of the calculations in ice_temp_SIS2.
subroutine estimate_tsurf(m_pond, m_snow, m_ice, enthalpy, sice, SF_0, dSF_dT, &
                          sol, tfw, tsurf, dtt, NkIce, CS, US, ITV)

  real, intent(in   ) :: m_pond  !< pond mass per unit area [R Z ~> kg m-2]
  real, intent(in   ) :: m_snow  !< snow mass per unit area [R Z ~> kg m-2]
  real, intent(in   ) :: m_ice   !< ice mass per unit area [R Z ~> kg m-2]
  real, dimension(0:NkIce) , &
        intent(in   ) :: enthalpy !< The enthalpy of each layer in a column of
                                  !! snow and ice, in enthalpy units [Q ~> J kg-1].
  real, dimension(NkIce), &
        intent(in)    :: Sice  !< ice salinity by layer [gSalt kg-1]
  real, intent(in   ) :: SF_0  !< net upward surface heat flux when Tsurf=0 [Q R Z T-1 ~> W m-2]
  real, intent(in   ) :: dSF_dT !< d(sfc heat flux)/d(ts) [Q R Z T-1 degC-1 ~> W m-2 degC-1]
  real, dimension(0:NkIce), &
        intent(in)    :: sol   !< Solar heating of the snow and ice layers [Q R Z T-1 ~> W m-2]
  real, intent(in   ) :: tfw   !< seawater freezing temperature [degC]
  real, intent(  out) :: tsurf !< surface temperature [degC]
  real, intent(in   ) :: dtt   !< timestep [T ~> s]
  integer, intent(in   ) :: NkIce !< The number of ice layers.
  type(SIS2_ice_thm_CS), intent(in) :: CS  !< The SIS2 ice thermodynamics control structure
  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.

!
! variables for temperature calculation [see Winton (1999) section II.A.]
! note:  here equations are multiplied by hi to improve thin ice accuracy
!
  real, dimension(0:NkIce) :: temp_est   ! An estimated snow and ice temperature [degC].
  real, dimension(0:NkIce) :: temp_IC    ! The temperatures of the snow and ice based on the initial
                                         ! enthalpy [degC].
  real, dimension(NkIce) :: tfi  ! The ice freezing temperatures [degC].
  real :: mL_ice    ! The mass-per-unit-area of each ice layer [R Z ~> kg m-2]
  real :: mL_snow   ! The mass-per-unit-area of each snow layer [R Z ~> kg m-2]
  real :: kk        ! The conductive thermal coupling coefficient between adjacent
                    ! ice sublayers [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: k10       ! The conductive thermal coupling coefficient between the
                    ! snow and the topmost ice sublayer [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: k0a       ! The implicit conductive thermal coupling coefficient between
                    ! the snow and the skin temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: k0skin    ! The conductive thermal coupling coefficient between the
                    ! snow and the skin temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: b_denom_1 ! A temporary variable [Q R Z degC-1 ~> J m-2 degC-1]
  real :: I_bb      ! The inverse of the tridiagonal pivot [degC Q-1 R-1 Z-1 ~> degC m2 J-1]
  real :: comp_rat  ! The complement of cc_bb, going from 0 to 1.
  real :: tsf       ! The surface freezing temperature [degC].
  real :: k0a_x_ta  ! The surface heat flux times normalized by 1 + the ratio
                    ! of the temperature feedback on surface fluces to the
                    ! skin-snow conductive sensitivity [Q R Z T-1 ~> W m-2].
  real :: tsurf_est ! An estimate of the surface temperature [degC].
  real, dimension(0:NkIce+1) :: cc ! Interfacial coupling coefficients [Q R Z degC-1 ~> J m-2 degC-1].
  real, dimension(0:NkIce) :: bb   ! Effective layer heat capacities [Q R Z degC-1 ~> J m-2 degC-1].
  real :: hsnow_eff  ! The effective thickness of the snow layer [Z ~> m].
  real :: hL_ice_eff ! The effective thickness of each ice sub-layer [Z ~> m].
  real :: rho_ice  ! The nominal density of sea ice [R ~> kg m-3].
  real :: rho_snow ! The nominal density of snow [R ~> kg m-3].
  real :: Cp_ice   ! The heat capacity of ice [Q degC-1 ~> J kg-1 degC-1].
  real :: Cp_brine ! The heat capacity of liquid water in the brine pockets,
                   ! [Q degC-1 ~> J kg-1 degC-1].
  real :: Lat_fus  ! The latent heat of fusion [Q degC-1 ~> J kg-1].
  integer :: k

  temp_IC(0) = temp_from_En_S(enthalpy(0), 0.0, ITV)
  call temp_from_Enth_S(enthalpy(1:), sice(1:), temp_IC(1:), ITV)

  call get_SIS2_thermo_coefs(ITV, rho_ice=rho_ice, rho_snow=rho_snow, &
                             Cp_Ice=Cp_ice, Latent_Fusion=Lat_fus, Cp_Brine=Cp_Brine)

  mL_ice = m_ice / NkIce   ! ice mass per unit area of each layer
  mL_snow = m_snow         ! snow mass per unit area (in kg m-2).
  call calculate_T_Freeze(sice, tfi, ITV)    ! freezing temperature of ice layers

  ! Set the effective thickness of each ice and snow layer, limited to avoid
  ! instabilities for thin layers.
  hL_ice_eff = max(mL_ice / rho_ice, CS%H_lo_lim)
  hsnow_eff = mL_snow / rho_snow + max(1.0e-35*US%m_to_Z, 1.0e-20*CS%H_lo_lim)

  ! Set the surface freezing temperature
  tsf = tfi(1) ; if (mL_snow>0.0) tsf = 0.0

  if ((dSF_dT*hsnow_eff + 2.0*CS%Ks) <= 0.0) then
    ! This is an unstable forcing situation, with the skin temperature effectively
    ! decoupled from the temperature of the snow layer, so the skin temperature
    ! will either go cold or warm.
    if (SF_0 <= 0.0) then ! The surface is being heated from above.
      tsurf = tsf
    else
      tsurf = tsf - 40.0 ! ### CS%dT_cold
    endif
    return
  endif

  kk = CS%Ki / hL_ice_eff       ! full ice layer conductivity
  k10 = 2.0*(CS%Ks*CS%Ki) / (hL_ice_eff*CS%Ks + hsnow_eff*CS%Ki) ! coupling ice layer 1 to snow
  k0a = (CS%Ks*dSF_dT) / (0.5*dSF_dT*hsnow_eff + CS%Ks)   ! coupling snow to "air"
  k0skin = 2.0*CS%Ks / hsnow_eff
  k0a_x_ta = (CS%Ks*SF_0) / (0.5*dSF_dT*hsnow_eff + CS%Ks) ! coupling times "air" temperature

  !   Estimate the surface temperature with a linearized fully implicit
  ! treatment of layer coupling.

  ! Determine the effective layer heat capacities.
  !   bb = dheat/dTemp, as derived from a linearization of the enthalpy equation
  !     but using the value for ice near the melt point if temp_IC is above T_fr.
  bb(0) = mL_snow* Cp_ice
  do k=1,NkIce   ! Store the heat capacity term in bb.
    if (tfi(k) >= 0.0) then ! This is pure ice.
      bb(k) = mL_ice * Cp_ice
    elseif (temp_IC(k) < tfi(k)) then
      bb(k) = mL_ice * (Cp_ice - (tfi(k) / temp_IC(k)**2) * &
                        (Lat_Fus - (Cp_brine-Cp_ice) * temp_IC(k)) )
      ! Or more generally:
      !  S_Sf = sice(k) / calculate_S_freeze(temp_IC(k))
      ! bb(k) = mL_ice * (Cp_ice + dSf_dT*Lat_fus + S_Sf * &
      !                   ((Cp_brine-Cp_ice))
    else ! Use the value when temp_IC = tfi.
      bb(k) = mL_ice * (Cp_brine - Lat_fus / tfi(k))
    endif
  enddo

  ! The following expressions could be modified to permit there to be
  ! variable diffusivities in the ice/brine mixtures.
  cc(0) = k0a*dtt  ! Atmosphere-snow coupling
  cc(1) = k10*dtt  ! Snow-ice coupling
  do k=2,NkIce ; cc(k) = kk*dtt ; enddo
  cc(NkIce+1) = 2*kk*dtt ! bottom of ice coupling with ocean at freezing point.

  !   This is a version of a tridiagonal solver where the relationship between
  ! the diagonal elements and the coupling between layers has been used to
  ! ensure that there are no differences ever taken, so there will be minimal
  ! truncation errors.  This closely follows schemes that have previously been
  ! found to work very well in GOLD and MOM6.

  ! Go up the ice column; there will be no need for the downward pass of the
  ! tridiagonal solver, since this routine only outputs the skin temperature.
  b_denom_1 = (bb(NkIce) + cc(NkIce+1))
  I_bb = 1.0 / (b_denom_1 + cc(NkIce))
  temp_est(NkIce) = ( (sol(NkIce)*dtt + bb(NkIce)*temp_IC(NkIce)) + &
                      cc(NkIce+1)*tfw ) * I_bb
  comp_rat = b_denom_1 * I_bb

  do k=NkIce-1,1,-1
    b_denom_1 = bb(k) + comp_rat*cc(k+1)
    I_bb =  1.0 / (b_denom_1 + cc(k))
    temp_est(k) = ((sol(k)*dtt + bb(k)*temp_IC(k)) + cc(k+1)*temp_est(k+1)) * I_bb

    comp_rat = b_denom_1 * I_bb ! 1.0 >= comp_rat >= 0.0
  end do

  I_bb =  1.0 / ((bb(0) + comp_rat*cc(1)) + cc(0))
  !   This is a complete calculation of temp_est(0), assuming that the upward surface
  ! flux is given by SF_0 + dSF_dT*tsurf, with tsurf estimated in the following line
  ! as part of this calculation.
  temp_est(0) = (((sol(0)*dtt + bb(0)*temp_IC(0)) - k0a_x_ta*dtt) + &
                   cc(1)*temp_est(1)) * I_bb

  ! Diagnose the surface skin temperature by matching the diffusive fluxes in
  ! the snow with the atmospheric fluxes.  I.e. solve the following for tsurf_est:
  !  (SF_0 + dSF_dT*tsurf_est) = k0skin * (temp_est(0) - tsurf_est)
  tsurf = min(tsf, (k0skin*temp_est(0) - SF_0) / (dSF_dT + k0skin))

end subroutine estimate_tsurf

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> laytemp_SIS2 does an implicit calculation of new layer temperature
function laytemp_SIS2(m_ice, T_fr, Qf, bf, tp, enth, salin, dtt, ITV, US) result (new_temp)
  real :: new_temp         !< The new temperature of the ice layer [degC]
  real, intent(in) :: m_ice !< mass of ice [R Z ~> kg m-2]
  real, intent(in) :: T_fr !< ice freezing temp [degC] (determined by salinity)
  real, intent(in) :: Qf   !< Inward forcing [Q R Z T-1 ~> W m-2]
  real, intent(in) :: bf   !< response of outward heat flux to local temperature [Q R Z T-1 ~> W m-2 degC]
  real, intent(in) :: tp   !< prior step temperature [degC]
  real, intent(in) :: enth !< prior step enthalpy
  real, intent(in) :: salin !< ice salinity [gSalt kg-1].
  real, intent(in) :: dtt  !< timestep [T ~> s]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors

  real :: T_g    ! The latest best guess at Temp [degC].
  real :: T_deriv  ! The value of Temp at which to evaluate dErr_dT [degC].
  real :: T_max, T_min ! Bracketing temperatures [degC].
  real :: Err    ! The error in the heat budget at T_guess [Q R Z ~> J m-2].
  real :: Err_Tmin, Err_Tmax ! The errors at T_min and T_max [Q R Z ~> J m-2].
  real :: T_prev  ! The previous value of T_g [degC].
  real :: dErr_dT ! The partial derivative of Err with T_g [Q R Z degC-1 ~> J m-2 degC-1].
  real :: Enth_tol = 1.0e-15 ! The fractional Enthalpy difference tolerance for convergence.
  real :: TfmxdCp_BI ! The vertically integrated freezing temperature times the difference
                   ! between the freezing capacities of brine and ice.  [Q R Z ~> J m-2]

  real :: E0   ! Starting heat relative to salinity dependent freezing [Q ~> J kg-1].
  real :: AA, BB, CC !  [Q R Z ~> J m-2]
  real :: Cp_ice   ! The heat capacity of ice [Q degC-1 ~> J kg-1 degC-1].
  real :: Cp_brine ! The heat capacity of liquid water in the brine pockets,
                   ! [Q degC-1 ~> J kg-1 degC-1].
  real :: Cp_water ! The heat capacity of liquid water in the ice model,
                   ! but not in the brine pockets [Q degC-1 ~> J kg-1 degC-1].
  real :: LI       ! The latent heat of fusion [Q ~> J kg-1].

  integer :: itt
!  real :: T_itt(20), dTemp(20), Err_itt(20)

  call get_SIS2_thermo_coefs(ITV, Cp_Ice=Cp_ice, Latent_Fusion=LI, &
                             Cp_Brine=Cp_Brine, Cp_Water=Cp_Water)

  if ( T_fr == 0.0 ) then
    ! For fresh water, avoid the degeneracy of the enthalpy-temperature
    ! relationship by extending the linear expression for frozen water, then
    ! limit it later to be at or below freezing.
    !
    !   m * {Cp_Ice} * (tn-tp) = dtt * (Qf - bf*tn)
    !
    new_temp = (m_ice*Cp_Ice*tp + Qf*dtt) / (m_ice*Cp_Ice + bf*dtt) ! = -BB/AA

  else
    if (tp >= T_fr) then
      E0 = Cp_Water*(tp - T_fr)  ! >= 0
    else
      E0 = Cp_Ice*(tp - T_fr) - LI*(1 - T_fr/tp)  ! < 0
      if (Cp_ice /= Cp_brine) E0 = E0 + (Cp_brine - Cp_ice) * T_fr*log(tp/T_fr)
    endif
    ! Determine whether the new solution will be above or below freezing.

    if (m_ice*E0 + dtt * (Qf - bf*T_fr) >= 0) then
      ! This layer will be completely melted, so return the freezing value.
      ! new_temp = T_fr + (m_ice*E0 + dtt*(Qf - bf*T_fr)) / (Cp_water*m_ice + dtt*bf)
      new_temp = T_fr
    else
      ! This layer will be partly melted.
      ! Solve a quadratic equation for the new layer temperature, tn:
      !
      !   m_ice * {Cp_Ice-LI*T_fr/(tn*tp)} * (tn-tp) = dtt * (Qf - bf*tn)
      !   m_ice * {{Cp_Ice*(tn - T_fr) - LI*(1 - T_fr/tn)} - E0} = dtt * (Qf - bf*tn)
      !        En0(tn) = Cp_Ice*(tn - T_fr) - LI*(1 - T_fr/tn)
      !
      AA = m_ice*Cp_Ice + bf*dtt
      BB = -(m_ice*((E0 + LI) + Cp_Ice*T_fr) + Qf*dtt)
      CC = m_ice*LI*T_fr
      ! This form avoids round-off errors.
      if (BB >= 0) then
        new_temp = -(BB + sqrt(BB*BB - 4*AA*CC)) / (2*AA)
      else
        new_temp = (2*CC) / (-BB + sqrt(BB*BB - 4*AA*CC))
      endif

      if (Cp_ice /= Cp_brine) then
      ! At this point, new_temp is just a good starting guess that needs to be iterated to convergence.

      ! Solve the following expression for the new layer temperature, tn:
      !
        TfmxdCp_BI = T_fr*m_ice*(Cp_Brine-Cp_Ice)
      !   Err = m_ice * ((Cp_Ice*(tn - T_fr) - LI*(1 - T_fr/tn)) - &
      !               TfmxdCp_BI*(log(T_fr/tn)) - E0) - dtt * (Qf - bf*tn)
      !   Err = -(m_ice*((E0 + LI) + Cp_Ice*T_fr) + f*dtt) + &
      !         m_ice * (Cp_Ice*tn + LI*(T_fr/tn)) - TfmxdCp_BI*(log(T_fr/tn))) + dtt*bf*tn

        ! This might be a good enough first guess that bracketing is unnecessary
        ! but it is better to play it safe...
        T_g = new_temp
        Err = BB + ((m_ice * (Cp_Ice*t_g + LI*(T_fr/t_g)) - &
                     TfmxdCp_BI*(log(T_fr/t_g))) + dtt*bf*t_g)

        if (Err <= 0.0) then
          T_min = T_g ; Err_Tmin = Err
          T_max = T_fr ; Err_Tmax = BB + (m_ice * (Cp_Ice*T_fr + LI) + dtt*bf*T_fr)
        else
          T_max = T_g ; Err_Tmax = Err
          T_min = -273.15
          Err_Tmin = BB + ((m_ice * (Cp_Ice*T_min + LI*(T_fr/T_min)) - &
                           TfmxdCp_BI*(log(T_fr/T_min))) + dtt*bf*T_min)
        endif

!        T_itt(:) = 0.0 ; dTemp(:) = 0.0 ; Err_itt(:) = 0.0
        do itt=1,20 ! Note that 3 or 4 iterations usually are enough.
          Err = BB + ((m_ice * (Cp_Ice*t_g + LI*(T_fr/t_g)) - &
                       TfmxdCp_BI*(log(T_fr/t_g))) + dtt*bf*t_g)
!          T_itt(itt) = T_g ; Err_itt(itt) = Err

          if (abs(Err) <= Enth_tol*(abs(BB) + m_ice*LI + abs(dtt*bf*T_fr))) then
            new_temp = T_g ; exit
          elseif (Err < 0.0) then
            T_min = T_g ; Err_Tmin = Err
          else
            T_max = T_g ; Err_Tmax = Err
          endif

        ! Use the more efficient Newton's method of McDougall & Witherspoon (2014),
        ! Appl. Math. Lett., 29, 20-25.
          T_deriv = T_g
          if (itt > 1) then ! Reuse the estimate of dT_dEn from the last iteration.
            if ((dErr_dT*T_g - 0.5*Err > dErr_dT*T_min) .and. &
                (dErr_dT*T_g - 0.5*Err < dErr_dT*T_max)) &
              T_deriv = T_g - 0.5*Err / dErr_dT
          endif

          dErr_dt = m_ice * (Cp_Ice - LI*T_fr/(t_deriv**2)) + (TfmxdCp_BI / t_deriv + bf*dtt) ! >= 0.0
          T_prev = T_g
          if ((dErr_dT*T_g - Err > dErr_dT*T_min) .and. &
              (dErr_dT*T_g - Err < dErr_dT*T_max)) then
            T_g = T_g - Err / dErr_dT
          else
            T_g = (Err_Tmax * T_min - Err_Tmin * T_max) / &
                  (Err_Tmax - Err_Tmin)
          endif
!          dTemp(itt) = T_g - T_prev

        enddo
        new_temp = T_g ! Use the best guess.

!        write (*,'("T_itt = ",8F14.8)') T_itt(1:8)
!        write (*,'("Err_itt = ",8(1PE14.6))') Err_itt(1:8)
!        write (*,'("dTemp = ",8(1Pe12.4))') dTemp(1:8)
      endif

    endif
  endif

  ! Only return temperatures that are at or below the freezing point.
  new_temp = min(new_temp, T_fr)

end function laytemp_SIS2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> update_lay_enth does an implicit calculation of new layer enthalpy
subroutine update_lay_enth(m_lay, sice, enth, ftop, ht_body, fbot, dftop_dT, &
                           dfbot_dT, dtt, hf_err_rat, ITV, US, extra_heat, temp_new, temp_max)
  real, intent(in) :: m_lay    !< This layers mass of ice [R Z ~> kg m-2]
  real, intent(in) :: sice     !< ice salinity [gSalt kg-1]
  real, intent(inout) :: enth  !< ice enthalpy in enthaly units [Q ~> J kg-1].
  real, intent(inout) :: ftop  !< Downward heat flux atop the layer at T = 0 degC, or
                               !! the prescribed heat flux if dftop_dT = 0 [Q R Z T-1 ~> W m-2].
  real, intent(in) :: ht_body  !< Body heating to layer [Q R Z T-1 ~> W m-2]
  real, intent(inout) :: fbot  !< Downward heat below the layer at T = 0 degC [Q R Z T-1 ~> W m-2].
  real, intent(in) :: dftop_dT !< The linearization of ftop with layer temperature [Q R Z T-1 ~> W m-2 degC-1].
  real, intent(in) :: dfbot_dT !< The linearization of fbot with layer temperature [Q R Z T-1 ~> W m-2 degC-1].
  real, intent(in) :: dtt      !< The timestep [T ~> s]
  real, intent(in) :: hf_err_rat  !< A conversion factor for comparing the errors
                               !! in explicit and implicit estimates of the updated
                               !! heat fluxes [degC T Q-1 ~> kg degC W-1].
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors
  real, intent(out) :: extra_heat !< The heat above the melt point [Q R Z ~> J m-2].
  real, optional, intent(out) :: temp_new !< The new temperature [degC].
  real, optional, intent(in)  :: temp_max !< The maximum new temperature [degC].

  real :: htg      ! The rate of heating of the layer [Q R Z T-1 ~> W m-2].
  real :: new_temp ! The new layer temperature [degC].
  real :: max_temp ! The maximum new layer temperature [degC].
  real :: max_enth ! The maximum new layer enthalpy [degC].
  real :: fb       ! The negative of the dependence of layer heating on
                   ! temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1]. fb > 0.
  real :: extra_enth ! Excess enthalpy above the melt point [Q R Z ~> J m-2].
  real :: enth_in  ! The initial enthalpy [Q ~> J kg-1].
  real :: enth_fp  ! The enthalpy at the freezing point [Q ~> J kg-1].
  real :: AA, BB, CC ! Temporary variables used to solve a quadratic equation [Q R Z ~> J m-2].
  real :: dT_dEnth ! The partial derivative of temperature with enthalpy
                   ! [degC Q-1 ~> degC kg J-1].
  real :: En_J     ! The enthalpy with 0 offset for liquid at 0 degC [Q ~> J kg-1]
  real :: T_fr     ! Ice freezing temperature (determined by bulk salinity) [degC].
  real :: fbot_in, ftop_in ! Input values of fbot and ftop [Q R Z T-1 ~> W m-2].
  real :: dflux_dtot_dT  ! A temporary work array [degC].
  real :: T_g      ! The latest best guess at Temp [degC].
  real :: T_deriv  ! The value of Temp at which to evaluate dErr_dT [degC].
  real :: T_max, T_min ! Bracketing temperatures [degC].
  real :: Err      ! The heat content error at T_guess [Q R Z ~> J m-2].
  real :: Err_Tmin, Err_Tmax ! The errors at T_max and T_min [Q R Z ~> J m-2].
  real :: Err_Tind ! The error at some temperature [Q R Z ~> J m-2]
  real :: T_prev   ! The previous value of T_g [degC].
  real :: dErr_dT  ! The partial derivative of Err with T_g [Q R Z degC-1 ~> J m-2 degC-1].
  real :: Enth_tol = 1.0e-15 ! The fractional Enthalpy difference tolerance for convergence [nondim].
  real :: TfxdCp_WI ! The freezing temperature times the heat capacity from water to ice [degC Q ~> degC J kg-1]
  real :: TfxdCp_BI ! The freezing temperature times the heat capacity from brine to ice [degC Q ~> degC J kg-1]
  real :: Cp_ice   ! The heat capacity of ice [Q degC-1 ~> J kg-1 degC-1].
  real :: Cp_brine ! The heat capacity of liquid water in the brine pockets,
                   ! [Q degC-1 ~> J kg-1 degC-1].
  real :: Cp_water ! The heat capacity of liquid water in the ice model,
                   ! but not in the brine pockets [Q degC-1 ~> J kg-1 degC-1].
  real :: LI       ! The latent heat of fusion [Q degC-1 ~> J kg-1].
  integer :: itt
  ! real :: T_itt(20), dTemp(20), Err_itt(20)

  call get_SIS2_thermo_coefs(ITV, Latent_Fusion=LI, &
                             Cp_Ice=Cp_ice, Cp_Brine=Cp_brine, Cp_Water=Cp_water)

  ! Solve m_lay*(enth - enth_in) + extra_heat = dt * (ht_body + ftop - fbot)
  !  ftop = ftop_in + temp*dftop_dT
  !  fbot = fbot_in + temp*dfbot_dT
  ! subject to  enth <= enth_fp and extra_heat >= 0

  ftop_in = ftop ; fbot_in = fbot
  htg = (ht_body + ftop_in) - fbot_in
  fb = -(dftop_dT - dfbot_dT)   !  = -dhgt_dt > 0

  extra_heat = 0.0 ; extra_enth = 0.0
  if (sice > 0.0) then
    T_fr = T_freeze(sice, ITV)
    enth_fp = enthalpy_liquid_freeze(sice, ITV)
  else
    T_fr = 0.0
    enth_fp = enth_from_TS(0.0, 0.0, ITV)
  endif
  max_temp = T_fr ; max_enth = enth_fp
  if (present(temp_max)) then ; if (temp_max < T_fr) then
    max_temp = temp_max ; max_enth = enth_from_TS(temp_max, sice, ITV)
  endif ; endif
  enth_in = enth

  ! Solve m_lay * (enth_new - enth) = dtt * (htg - fb*t_new)
  !       t_new = Temp_from_En_S(enth_new, Sice, ITV)
  if (m_lay == 0.0) then
    new_temp = min(htg / fb, max_temp)
    enth = enth_from_TS(new_temp, sice, ITV)
  elseif (dtt * (htg - fb*max_temp) >= m_lay*(max_enth - enth_in)) then
    ! There is enough heat being applied here that the ice would be above the
    ! maximum temperature (often the freezing point).  The ice should be set to
    ! the maximum temperature and enthalpy, and the extra heat stored for later
    ! use in melting.
    extra_enth = m_lay*(enth_in - max_enth) + dtt * (htg - fb*max_temp)
    extra_heat = extra_enth
    new_temp = max_temp
    enth = max_enth
  elseif ( sice == 0.0 ) then  ! Note that T_fr = 0.
    ! dT_dEnth is 0 for enth > enth_fp.
    !   dT_dEnth = dTemp_dEnth_EnS(enth_in, Sice, ITV)
    dT_dEnth = 1.0 / Cp_Ice

    ! Solve for enth:  m_lay  * (enth - enth_in) =
    !       dtt * (htg - fb*T_fr - fb*dT_dEnth*(enth - enth_fp))
    !  enth = enth_in + dtt * (htg - fb*(0.0 - dT_dEnth*(enth_fp-enth_in))) / &
    !                          (m_lay + dtt*b*dT_dEnth)
    ! Or equivalently...  (noting that T_fr = 0.0)
    enth = enth_fp + (dtt * (htg - fb*0.0) + m_lay * (enth_in-enth_fp)) / &
                     (m_lay  + dtt*(fb*dT_dEnth))
    ! The following is equivalent to new_temp = Temp_from_En_S(enth, 0.0, ITV)
    !     or  new_temp = dT_dEnth * (enth - enth_fp) ! + T_fr==0.
    ! but it avoids serious roundoff issues later on when b is large.
    new_temp = dT_dEnth * ((dtt * (htg - fb*0.0) + m_lay * (enth_in-enth_fp)) / &
                           (m_lay  + dtt*(fb*dT_dEnth)))
  else
    En_J = (enth_in - enthalpy_liquid(0.0, 0.0, ITV))
    ! For a first guess solve the case if (Cp_ice == Cp_brine) by solving a
    ! quadratic equation for the new layer temperature, tn:
    !
    !   m * (En_J - (Cp_Water-Cp_Ice)*T_fr + LI + htg*dt/m) =
    !        (m*Cp_Ice + b*dt) *tn + m*LI*T_fr/tn
    !
    AA = m_lay *Cp_Ice + fb*dtt
    BB = -(m_lay*((En_J - (Cp_Water-Cp_Ice)*T_fr) + LI) + htg*dtt)
    CC = m_lay * LI*T_fr
    ! This form avoids round-off errors.
    if (BB >= 0) then
      new_temp = -(BB + sqrt(BB*BB - 4*AA*CC)) / (2*AA)
    else
      new_temp = (2*CC) / (-BB + sqrt(BB*BB - 4*AA*CC))
    endif
    if (Cp_ice == Cp_brine) then  ! Keep this solution.
      enth = enth_from_TS(new_temp, sice, ITV)
    else ! (Cp_ice /= Cp_brine)
      ! Correct the new temperature estimate.
      ! Solve for enth & -273.15 < T_g < T_fr < 0
      ! m_lay*(enth - En_J) = dtt * (htg - fb*T_g)
      ! enth = (-LI * (1.0 - T_fr/T_g)) + &
      !         ((Cp_Ice*Tg + TfxdCp_WI) - TfxdCp_BI*log(T_fr/T_g))

      ! Err = m_lay*(enth - En_J) + dtt * (fb*T_g - htg)
      ! Err = m_lay*((-LI * (1.0 - T_fr/T_g)) + &
      !       ((Cp_Ice*Tg + TfxdCp_WI) - TfxdCp_BI*log(T_fr/T_g)) - En_J) + dt * (fb*T_g - htg)

      ! En_J = enth_in - enth_liq_0
      TfxdCp_WI = T_fr*(Cp_Water-Cp_Ice)
      TfxdCp_BI = T_fr*(Cp_Brine-Cp_Ice)
      Err_Tind = (m_lay*(-LI + TfxdCp_WI - En_J) - dtt*htg)

      T_min = -273.15
      Err_Tmin = m_lay*(LI * (T_fr/T_min) + (Cp_Ice*T_min - TfxdCp_BI*log(T_fr/T_min))) + &
            (dtt * fb * T_min + Err_Tind)
      ! Approximate this as
      ! Err_Tmin = m_lay*((Cp_Ice*T_min - TfxdCp_BI*log(T_fr/T_min))) + &
      !      (dtt * fb * T_min + Err_Tind) ?
      T_max = T_fr
      Err_Tmax = m_lay*(T_fr*Cp_Water - En_J) + dtt * (fb * T_fr - htg)

      T_g = new_temp
      ! Using a false position method first-guess instead adds about 2 iterations.
      !  T_g = (Err_Tmax * T_min - Err_Tmin * T_max) / (Err_Tmax - Err_Tmin)

  !   T_itt(:) = 0.0 ; dTemp(:) = 0.0 ; Err_itt(:) = 0.0
      do itt=1,20 ! Note that 3 or 4 iterations usually are enough.
        Err = m_lay*(LI * (T_fr/T_g) + (Cp_Ice*T_g - TfxdCp_BI*log(T_fr/T_g))) + &
              (dtt * fb * T_g + Err_Tind)
  !     T_itt(itt) = T_g ; Err_itt(itt) = Err

        if (abs(Err) <= Enth_tol*(abs(Err_Tind) + m_lay*LI + abs(dtt*fb*T_fr))) then
          new_temp = T_g ; exit
        elseif (Err < 0.0) then
          T_min = T_g ; Err_Tmin = Err
        else
          T_max = T_g ; Err_Tmax = Err
        endif

      ! Use the more efficient Newton's method of McDougall & Witherspoon (2014),
      ! Appl. Math. Lett., 29, 20-25.
        T_deriv = T_g
        if (itt > 1) then ! Reuse the estimate of dT_dEn from the last iteration.
          if ((dErr_dT*T_g - 0.5*Err > dErr_dT*T_min) .and. &
              (dErr_dT*T_g - 0.5*Err < dErr_dT*T_max)) &
            T_deriv = T_g - 0.5*Err / dErr_dT
        endif

        dErr_dT = m_lay*( -LI * (T_fr / T_deriv**2) + &
                         (Cp_Ice + TfxdCp_BI / T_deriv)) + dtt*fb ! >= 0.0
  !     T_prev = T_g
        if ((dErr_dT*T_g - Err > dErr_dT*T_min) .and. &
            (dErr_dT*T_g - Err < dErr_dT*T_max)) then
          T_g = T_g - Err / dErr_dT
        else
          T_g = (Err_Tmax * T_min - Err_Tmin * T_max) / &
                (Err_Tmax - Err_Tmin)
        endif
  !     dTemp(itt) = T_g - T_prev

      enddo
      new_temp = T_g ! Use the best guess.

      enth = enth_from_TS(new_temp, sice, ITV)
    endif
  endif


!    Figure out whether it is more accurate to use the explicit form for the
!  fluxes or to invert enthalpy conservation for the fluxes.
  if (abs(hf_err_rat*dftop_dT) <= m_lay) then
    ! This branch currently applies in all interior ice layers.
    ftop = ftop_in + dftop_dT*new_temp
    ! An explicit estimate (or no update) is used for the top flux.
    if (hf_err_rat*dfbot_dT <= m_lay) then  ! Use the explicit expression for fbot.
      fbot = fbot_in + dfbot_dT*new_temp
    else ! Use conservation to invert for fbot.
      fbot = (ht_body + ftop) - (m_lay*(enth - enth_in) + extra_enth)/dtt
    endif
  elseif (hf_err_rat*dfbot_dT <= m_lay) then
    ! Use the explicit expression for fbot and invert for ftop.
    fbot = fbot_in + dfbot_dT*new_temp
    ftop = (fbot - ht_body) + (m_lay*(enth - enth_in) + extra_enth)/dtt
  else
    ! Conservation is used to invert for both fbot and ftop, partitioning
    ! the changes in proportion to their sensitivities.
    !   dflux = (htg - (m*(enth - enth_in) + extra_enth)/dtt)
    !   dflux = dfbot - dftop ; dftop / dfbot = dftop_dT / dfbot_dT

    if (dfbot_dT - dftop_dT > 0.0) then
      dflux_dtot_dT = (htg - (m_lay*(enth - enth_in) + extra_enth)/dtt) / &
                      (dfbot_dT - dftop_dT)
    else
      dflux_dtot_dT = 0.0 ! This should never occur.
    endif

    ftop = ftop_in + dftop_dT * dflux_dtot_dT
    fbot = fbot_in + dfbot_dT * dflux_dtot_dT
  endif

  if (present(temp_new)) temp_new = new_temp

end subroutine update_lay_enth

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_check evaluates whether the ice has any extreme temperatures or masses
!! and writes messages about any offending columns.
subroutine ice_check(ms, mi, enthalpy, s_ice, NkIce, msg_part, ITV, &
                      bmelt, tmelt, t_sfc, US)
  real, intent(in) :: ms !< The mass of snow [R Z ~> kg m-2] or [kg m-2]
  real, intent(in) :: mi !< The mass of ice [R Z ~> kg m-2] or [kg m-2]
  real, dimension(0:NkIce), &
        intent(in) :: enthalpy !< The ice enthalpy, in enthalpy units [Q ~> J kg-1]
  real, dimension(NkIce), intent(in) :: s_ice !< The ice bulk salinity [gSalt kg-1]
  integer, intent(in) :: NkIce !< The number of vertical temperature layers in the ice
  character(len=*), intent(in) :: msg_part !< An identifying message
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real, optional, intent(in) :: bmelt  !< The heat flux assocated with bottom melt [Q R Z ~> J m-2] or [J m-2]
  real, optional, intent(in) :: tmelt  !< The heat flux assocated with top melt [Q R Z ~> J m-2] or [J m-2]
  real, optional, intent(in) :: t_sfc  !< The ice surface temperature [degC]
  type(unit_scale_type), optional, intent(in) :: US  !< A structure with unit conversion factors

  ! Local variables
  real :: RZ_scale  ! A column mass conversion factor
  real :: QRZ_scale  ! A column mass conversion factor
  character(len=300) :: mesg
  character(len=80) :: msg2
  integer :: k, bad
  real, dimension(0:NkIce) :: t_col ! Temperatures of the ice and slow layers [degC]

  RZ_scale = 1.0 ; if (present(US)) RZ_scale = US%RZ_to_kg_m2
  QRZ_scale = 1.0 ; if (present(US)) QRZ_scale = US%Q_to_J_kg*US%RZ_to_kg_m2

  bad = 0
  if (present(t_sfc)) then
    if ((t_sfc > 1.e-14) .or. (t_sfc < -100.0)) bad = bad+1
  endif
  if (RZ_scale*ms <0.0 .or. RZ_scale*ms > 3.30e5  ) bad = bad+1
  if (RZ_scale*mi <0.0 .or. RZ_scale*mi > 1.0e6  ) bad = bad+1

  t_col(0) = temp_from_En_S(enthalpy(0), 0.0, ITV)
  call temp_from_Enth_S(enthalpy(1:), s_ice(:), t_col(1:), ITV)
  do k=0,NkIce ; if ((t_col(k) > 0.0) .or. (t_col(k) < -100.0)) bad = bad+1 ; enddo

  if (bad>0) then
    mesg = "BAD ICE "//trim(msg_part)
    write (msg2,'(" ms,mi=",2(ES11.3))') ms*RZ_scale*1e-3,mi*RZ_scale*1e-3 ; mesg = trim(mesg)//trim(msg2)
    if (present(t_sfc)) then
      write (msg2,'(" t_sfc,tsn,tice=",ES11.3)') t_sfc ; mesg = trim(mesg)//trim(msg2)
    else
      mesg = trim(mesg)//"tsn,tice="
    endif
    do k=0,NkIce ; write (msg2,'(ES11.3)') t_col(k) ; mesg = trim(mesg)//trim(msg2) ; enddo
    if (present(bmelt)) then
      write (msg2,'(" bmelt=",ES11.3)') QRZ_scale*bmelt ; mesg = trim(mesg)//trim(msg2)
    endif
    if (present(tmelt)) then
      write (msg2,'(" tmelt=",ES11.3)') QRZ_scale*tmelt ; mesg = trim(mesg)//trim(msg2)
    endif
    call SIS_error(WARNING, mesg, all_print=.true.)
  endif

end subroutine ice_check

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_resize_SIS2 makes snow and ice thickness and temperature changes in ice
!! columns due to thermodynamic forcing.
subroutine ice_resize_SIS2(a_ice, m_pond, m_lay, Enthalpy, Sice_therm, Salin, &
                           snow, rain, evap, tmlt, bmlt, NkIce, npassive, TrLay, &
                           heat_to_ocn, h2o_ice_to_ocn, h2o_ocn_to_ice, evap_from_ocn, &
                           snow_to_ice, salt_to_ice, ITV, US, CS, ablation, &
                           enthalpy_evap, enthalpy_melt, enthalpy_freeze)
  ! mw/new - melt pond - added first two arguments & rain
  real, intent(in   ) :: a_ice       !< area of ice (1-open_water_frac) for pond retention [nondim]
  real, intent(inout) :: m_pond      !< melt pond mass [R Z ~> kg m-2]
  real, dimension(0:NkIce), &
        intent(inout) :: m_lay       !< Snow and ice mass per unit area by layer [R Z ~> kg m-2].
  real, dimension(0:NkIce+1), &
        intent(inout) :: Enthalpy    !< Snow, ice, and ocean enthalpy by layer [Q ~> J kg-1].
  real, dimension(NkIce), &
        intent(in)    :: Sice_therm  !< ice salinity by layer, as used for thermodynamics [gSalt kg-1]
  real, dimension(NkIce+1), &
        intent(inout) :: Salin       !< Conserved ice bulk salinity by layer [gSalt kg-1]
  real, intent(in   ) :: snow        !< new snow [R Z ~> kg m-2]
  real, intent(in   ) :: rain        !< rain for pond source [R Z ~> kg m-2] - not yet active
  real, intent(in   ) :: evap        !< ice evaporation/sublimation [R Z ~> kg m-2]
  real, intent(in   ) :: tmlt        !< top melting energy [Q R Z ~> J m-2]
  real, intent(in   ) :: bmlt        !< bottom melting energy [Q R Z ~> J m-2]
  integer, intent(in) :: NkIce       !< The number of ice layers.
  integer, intent(in) :: npassive    !< Number of passive tracers
  real, dimension(0:NkIce+1,npassive), &
        intent(inout) :: TrLay       !< Passive tracer slice
  real, intent(  out) :: heat_to_ocn !< energy left after ice all melted [Q R Z ~> J m-2]
  real, intent(  out) :: h2o_ice_to_ocn !< liquid water flux to ocean [R Z ~> kg m-2]
  real, intent(  out) :: h2o_ocn_to_ice !< liquid water flux from ocean [R Z ~> kg m-2]
  real, intent(  out) :: evap_from_ocn  !< evaporation flux from ocean [R Z ~> kg m-2]
  real, intent(  out) :: snow_to_ice !< snow below waterline becomes ice [R Z ~> kg m-2]
  real, intent(  out) :: salt_to_ice !< Net flux of salt to the ice [R Z gSalt kg-1 ~> gSalt m-2].
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors
  type(SIS2_ice_thm_CS), intent(in) :: CS  !< The SIS2_ice_thm control structure.

  real, intent(  out) :: ablation      !< The mass loss from bottom melt [R Z ~> kg m-2].
  real, intent(  out) :: enthalpy_evap !< The enthalpy loss due to the mass loss
                                       !! by evaporation / sublimation. [Q kg m-2 ~> J m-2]
  real, intent(  out) :: enthalpy_melt !< The enthalpy loss due to the mass loss
                                       !! by melting [Q R Z ~> J m-2].
  real, intent(  out) :: enthalpy_freeze !< The enthalpy gain due to the mass gain
                                       !! by freezing [Q kg m-2 ~> J m-2].

  real :: top_melt, bot_melt, melt_left ! Heating amounts, all in [Q R Z ~> J m-2]
  real :: mtot_ice    ! The summed ice mass [R Z ~> kg m-2].
  real :: enth_freeze ! The enthalpy of newly formed congelation ice [Q ~> J kg-1].
  real, dimension(0:NkIce) :: enth_fr ! The snow and ice layers' freezing point
                                      ! enthalpy [Q ~> J kg-1].
  real :: min_dEnth_freeze    ! The minimum enthalpy change that must occur when freezing water,
                              ! usually enough to account for the latent heat of fusion
                              ! in a small fraction of the water [Q ~> J kg-2].
  real :: m_freeze            ! The newly formed ice from freezing [R Z ~> kg m-2].
  real :: M_melt              ! The ice mass lost to melting [R Z ~> kg m-2].
  real :: evap_left           ! The remaining evaporation [R Z ~> kg m-2].
  real :: evap_here           ! The evaporation from the current layer [R Z ~> kg m-2].
  real :: m_submerged         ! The submerged mass of ice [R Z ~> kg m-2].
  real :: salin_freeze        ! The salinity of newly frozen ice [gSalt kg-1].
  real :: enthM_evap, enthM_melt, enthM_freezing, enthM_snowfall ! [Q R Z ~> J m-2]
  real :: LI          ! The latent heat of fusion [Q ~> J kg-1].
  real :: Lat_vapor   ! The latent heat of vaporization [Q ~> J kg-1].
  real :: rho_ice     ! The nominal density of sea ice [R ~> kg m-3].
  real :: rho_water   ! The nominal density of seawater [R ~> kg m-3].
  real :: h2o_to_ocn, h2o_orig, h2o_imb  ! Water budget terms [R Z ~> kg m-2]
  real :: pond_rate   ! The fraction of the meltwater that is retained in ponds [nondim]
  real :: h2o_to_pond, h2o_from_pond  ! Water budget terms [R Z ~> kg m-2]
  real :: tavg        ! Mass-weighted vertically averaged ice temperature [degC]
  real :: mp_min, mp_max ! Melt pond mass load limits [R Z ~> kg m-2]
  integer :: k, tr
  logical :: debug = .false.

  call get_SIS2_thermo_coefs(ITV, Latent_Fusion=LI, &
                             Latent_vapor=Lat_vapor, Rho_water=rho_water, Rho_ice=rho_ice)
  min_dEnth_freeze = LI * (1.0-CS%liq_lim)

  pond_rate = CS%r_min_pond + (CS%r_max_pond-CS%r_min_pond)*a_ice

  top_melt = tmlt ; bot_melt = bmlt

  ! set mass mark; will subtract the mass at end to find the melt flux to ocean
  if (debug) then
    h2o_orig = 0.0 ; do k=0,NkIce ; h2o_orig = h2o_orig + m_lay(k) ; enddo
  endif

  enth_fr(0) = enthalpy_liquid_freeze(0.0, ITV)
  do k=1,NkIce ! break out individual layers
    enth_fr(k) = enthalpy_liquid_freeze(sice_therm(k), ITV)
  enddo

  heat_to_ocn = 0.0   ! for excess melt energy

  evap_from_ocn = 0.0 ! for excess evap-melt
  h2o_ocn_to_ice = 0.0 ; h2o_ice_to_ocn = 0.0 ; snow_to_ice  = 0.0
  h2o_to_pond = 0.0
  h2o_from_pond = 0.0
  salt_to_ice = 0.0
  enthM_freezing = 0.0 ; enthM_melt = 0.0 ; enthM_evap = 0.0 ; enthM_snowfall = 0.0

  ! raining on cold ice led to unphysical temperature oscillations
  ! in single column test; need to pass all rain through to ocean for now - mw
  ! m_pond = m_pond + pond_rate*rain ! mw/new pond intercepts rain
  ! h2o_ice_to_ocn = h2o_ice_to_ocn + (1-pond_rate)*rain
  ! h2o_ice_to_ocn = h2o_ice_to_ocn + rain

  ! Delete this later, since it should not happen.
  mtot_ice = 0.0 ; do k=1,NkIce ; mtot_ice = mtot_ice + m_lay(k) ; enddo
  if (mtot_ice == 0.0) m_lay(0) = 0.0  ! This should already be true! Trap an error?  Convert the snow to ice?

  m_lay(0) = m_lay(0) + snow ! add snow
  enthM_snowfall = snow*enthalpy(0)
  if (evap < 0.0) then
    m_lay(0) = m_lay(0) - evap ! Treat frost formation like snow.
    enthM_snowfall = enthM_snowfall - evap*enthalpy(0)
  endif

  if (top_melt < 0.0 .and. CS%do_pond) then ! mw/new: add fresh/0C ice to top layer
    ! enth_freeze = -LI   ! this is right for prognostic salinity (i think)
    enth_freeze = Enthalpy(1) ! this is right for fixed salinity
    m_freeze = top_melt/enth_freeze;
    if (m_freeze > m_pond) then
      ! Excess freezing energy goes to cooling
      Enthalpy(1) = (m_lay(1)*Enthalpy(1) + m_freeze*enth_freeze) / (m_lay(1) + m_pond)
      Salin(1) = (m_lay(1)*Salin(1)) / (m_lay(1) + m_pond) ! for bulk salinity
      ! Passive tracer array
      do tr = 1,npassive
        TrLay(1,tr) = (m_lay(1)*TrLay(1,tr)) / (m_lay(1) + m_pond)
      enddo

      m_lay(1) = m_lay(1) + m_pond
      m_pond = 0.0
    else
      Enthalpy(1) = (m_lay(1)*Enthalpy(1) + m_freeze*enth_freeze) / (m_lay(1) + m_freeze)
      Salin(1) = (m_lay(1)*Salin(1)) / (m_lay(1) + m_freeze) ! for bulk salinity
      do tr = 1,npassive
        TrLay(1,tr) = (m_lay(1)*TrLay(1,tr)) / (m_lay(1) + m_freeze) ! for passive tracers
      enddo
      m_lay(1) = m_lay(1) + m_freeze
      m_pond = m_pond - m_freeze
    endif
    top_melt = 0.0
  endif

  if (top_melt < 0.0 .and. .not. CS%do_pond) then ! shouldn't happen but can handle
    bot_melt = bot_melt + top_melt
    top_melt = 0.0
  endif

  ! Assume that Salin(NkIce+1) already is the bottom freezing salinity.
  salin_freeze = Salin(NkIce+1)
  if (bot_melt < 0.0 ) then ! add freezing to bottom layer at tice and salin_freeze.
    ! Enth_freeze is based on the colder of the properties of the existing ice
    ! or the heat content of ocean water after a small fraction has frozen.
    enth_freeze = min(Enthalpy(NkIce), enthalpy(NkIce+1) - min_dEnth_freeze)
    m_freeze = -bot_melt / (Enthalpy(NkIce+1) - enth_freeze)

    Enthalpy(NkIce) = (m_lay(NkIce)*Enthalpy(NkIce) + m_freeze*enth_freeze) / (m_lay(NkIce) + m_freeze)
    Salin(NkIce) = (m_lay(NkIce)*Salin(NkIce) + m_freeze*salin_freeze) / (m_lay(NkIce) + m_freeze)

    do tr=1,npassive
      ! Change passive tracer, assuming that surface boundary value is stored in NkIce+1
      TrLay(NkIce,tr) = (m_lay(NkIce)*TrLay(NkIce,tr) + m_freeze*TrLay(NkIce+1,tr)) / &
                        (m_lay(NkIce) + m_freeze)
    enddo

    salt_to_ice = salt_to_ice + m_freeze*salin_freeze

    m_lay(NkIce) = m_lay(NkIce) + m_freeze
    h2o_ocn_to_ice = h2o_ocn_to_ice + m_freeze
    enthM_freezing = enthM_freezing + m_freeze*enthalpy(NkIce+1)

    bot_melt = 0.0
  endif

  ! Apply mass losses from evaporation.

  if (evap > 0.0) then ! apply evaporation mass flux
    evap_left = evap
    do k=0,NkIce
      evap_here = min(evap_left, m_lay(k))
      evap_left = evap_left - evap_here
      m_lay(k) = m_lay(k) - evap_here
      ! Assume that evaporation does not make ice salty?
      if (k>0) salt_to_ice = salt_to_ice - Salin(k) * evap_here
      enthM_evap = enthM_evap + evap_here * enthalpy(k)

      if (evap_left <= 0.0) exit
    enddo

    evap_from_ocn = evap_left
    !   The energy required to evaporate was already taken into account in
    ! ice_thm, but there is excess energy that has not been used here that needs
    ! to be passed on to the ocean.
    heat_to_ocn = heat_to_ocn + evap_left*(Lat_vapor+LI)
  endif

  if (top_melt > 0.0 ) then ! apply top melt heat flux
    melt_left = top_melt
    do k=0,NkIce
      if (melt_left < m_lay(k) * (enth_fr(k) - Enthalpy(k))) then
        M_melt = melt_left / (enth_fr(k) - Enthalpy(k)) ! melt part of this layer
        melt_left = 0.0
      else
        M_melt = m_lay(k) ! melt this whole layer
        melt_left = melt_left - M_melt*(enth_fr(k) - Enthalpy(k))
      endif
      m_lay(k) = m_lay(k) - M_melt
      if (k>0) salt_to_ice = salt_to_ice - Salin(k) * M_melt
      if (CS%do_pond) then
        h2o_to_pond = h2o_to_pond + pond_rate*M_melt
        h2o_ice_to_ocn = h2o_ice_to_ocn + (1-pond_rate)*M_melt
      else
        h2o_ice_to_ocn = h2o_ice_to_ocn + M_melt
      endif
      enthM_melt = enthM_melt + M_melt*enth_fr(k)

      if (melt_left <= 0.0) exit ! All melt energy has been used.
    enddo

    m_pond = m_pond + h2o_to_pond ! mw/new - add to pond from rain and surface melt
    heat_to_ocn = heat_to_ocn + melt_left ! melt heat left after snow & ice gone
  endif

  ! apply bottom melt heat flux

  melt_left = bot_melt ; ablation = 0.0
  if (melt_left > 0.0) then ! melt ice and snow from below
    do k=NkIce,0,-1
      if (melt_left < m_lay(k) * (enth_fr(k) - Enthalpy(k))) then
        M_melt = melt_left / (enth_fr(k) - Enthalpy(k)) ! melt part of this layer
        melt_left = 0.0
      else
        M_melt = m_lay(k) ! melt this whole layer
        melt_left = melt_left - M_melt*(enth_fr(k) - Enthalpy(k))
      endif
      m_lay(k) = m_lay(k) - M_melt
      if (k>0) salt_to_ice = salt_to_ice - Salin(k) * M_melt
      h2o_ice_to_ocn = h2o_ice_to_ocn + M_melt
      enthM_melt = enthM_melt + M_melt*enth_fr(k)
      ablation = ablation + M_melt

      if (melt_left <= 0.0) exit ! All melt energy has been used.
    enddo

    heat_to_ocn = heat_to_ocn + melt_left
  endif

  ! There are no further heat or mass losses or gains by the ice+snow.
  enthalpy_evap = enthM_evap
  enthalpy_melt = enthM_melt
  enthalpy_freeze = enthM_freezing

  ! calculate total ice for pond drainage and waterline adjustments below
  mtot_ice = 0.0 ; do k=1,NkIce ; mtot_ice = mtot_ice + m_lay(k) ; enddo

  if ( m_pond > 0.0 ) then ! consider pond drainage through ice
    tavg = 0.0
    do k=1,NkIce ! calculate liquid fraction as in Hunke et al 2013
      tavg = tavg + Temp_from_En_S(Enthalpy(k), Sice_therm(k), ITV)*m_lay(k)
    end do
    tavg = tavg/mtot_ice  ! average ice temperature
    if (tavg > CS%tdrain) then ! drain pond based on tunable ice temp. criterion
      mp_max = mtot_ice*((rho_water/rho_ice) - 1.0)
      mp_min = mp_max*(CS%min_pond_frac/CS%max_pond_frac)**2
      h2o_ice_to_ocn = h2o_ice_to_ocn + max(m_pond - mp_min, 0.0)
      m_pond = m_pond - max(m_pond - mp_min, 0.0)
    endif
  endif

  ! calculate mass to take from pond to bring ice top to waterline - mw/new
  h2o_from_pond = m_pond + m_lay(0) - ((rho_water/rho_ice) - 1.0)*mtot_ice
  if (h2o_from_pond > 0.0) then ! reduce pond to raise ice top toward waterline
    h2o_from_pond = min(h2o_from_pond, m_pond) ! keep m_pond >= 0.0
    m_pond = m_pond - h2o_from_pond
    h2o_ice_to_ocn = h2o_ice_to_ocn + h2o_from_pond ! pass dumped h2o to ocean
  endif

  ! The mass of ice that must be submerged (when floating according to
  ! Archimedes principle) for ice top to be at waterline.
  m_submerged = (mtot_ice+m_lay(0)+m_pond)* (rho_ice/rho_water)
  if (m_submerged > mtot_ice) then
    snow_to_ice = min(m_submerged - mtot_ice, m_lay(0)) ! need ice from snow

    m_lay(0) = m_lay(0) - snow_to_ice

    ! Add ice to the topmost layer and dilute its salinity.
    Enthalpy(1) = (m_lay(1)*Enthalpy(1) + snow_to_ice*Enthalpy(0)) / &
                  (m_lay(1) + snow_to_ice)
    Salin(1) = Salin(1) * m_lay(1) / (m_lay(1) + snow_to_ice)
    do tr = 1,npassive
      TrLay(1,tr) = TrLay(1,tr) * m_lay(1) / (m_lay(1) + snow_to_ice)
    enddo
    m_lay(1) = m_lay(1) + snow_to_ice
  else
    snow_to_ice = 0.0
  endif

  if (debug) then
    mtot_ice = 0.0 ; do k=1,NkIce ; mtot_ice = mtot_ice + m_lay(k) ; enddo
    h2o_to_ocn = h2o_orig + snow - (evap-evap_from_ocn) - (m_lay(0) + mtot_ice)
    h2o_imb = h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)
    if (abs(h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)) > &
        max(1.0e-10*US%kg_m3_to_R*US%m_to_Z, 1.0e-12*h2o_orig, &
            1.0e-12*(abs(h2o_ice_to_ocn)+abs(h2o_ocn_to_ice)))) then
      h2o_imb = h2o_to_ocn - (h2o_ice_to_ocn - h2o_ocn_to_ice)
    endif

    call ice_check(m_lay(0), mtot_ice, enthalpy, Sice_therm, &
                      NkIce, "at end of ice_resize_SIS2", ITV, &
                      bmelt=bot_melt, tmelt=top_melt, US=US)
  endif

end subroutine ice_resize_SIS2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> add_frazil_SIS2 increases the mass in ice columns due to the accretion of frazil ice.
subroutine add_frazil_SIS2(m_lay, Enthalpy, Sice_therm, Salin, npassive, TrLay, &
                           frazil, tfw, NkIce, h2o_ocn_to_ice, &
                           salt_to_ice, ITV, US, CS, enthalpy_freeze)
  real, dimension(0:NkIce), &
        intent(inout) :: m_lay       !< Snow and ice mass per unit area by layer [R Z ~> kg m-2].
  real, dimension(0:NkIce+1), &
        intent(inout) :: Enthalpy    !< Snow, ice, and ocean enthalpy by layer [Q ~> J kg-1].
  real, dimension(NkIce), &
        intent(in)    :: Sice_therm  !< ice salinity by layer, as used for thermodynamics [gSalt kg-1]
  real, dimension(NkIce+1), &
        intent(inout) :: Salin       !< Conserved ice bulk salinity by layer [gSalt kg-1]
  integer, intent(in) :: npassive    !< Number of passive tracers
  real, dimension(NkIce+1,npassive), &
        intent(inout) :: TrLay       !< Passive tracer in the column layer [Conc]
  real, intent(in   ) :: frazil      !< frazil in energy units [Q R Z ~> J m-2]
  real, intent(in   ) :: tfw         !< seawater freezing temperature [degC]
  integer, intent(in) :: NkIce       !< The number of ice layers.
  real, intent(  out) :: h2o_ocn_to_ice !< liquid water flux from ocean [R Z ~> kg m-2]
  real, intent(  out) :: salt_to_ice !< Net flux of salt to the ice [R Z gSalt kg-1 ~> gSalt m-2].
  type(SIS2_ice_thm_CS), intent(in) :: CS  !< The SIS2_ice_thm control structure.
  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real, intent(  out) :: enthalpy_freeze !< The enthalpy gain due to the
                                     !! mass gain by freezing [Q R Z ~> J m-2].

  ! Local variables
  real :: enth_frazil       ! The enthalpy of newly formed frazil ice [Q ~> J kg-1].
  real :: frazil_per_layer  ! The frazil heat sink from each of the sublayers of
                            ! of the ice [Q R Z ~> J kg-1].
  real :: t_frazil          ! The temperature which with the frazil-ice is created [degC].
  real :: m_frazil          ! The newly-formed mass per unit area of frazil ice [R Z ~> kg m-2].
  real :: min_dEnth_freeze  ! The minimum enthalpy change that must occur when freezing water,
                            ! usually enough to account for the latent heat of fusion
                            ! in a small fraction of the water [Q ~> J kg-2].
  real :: salin_freeze      ! The salinity of newly frozen ice [gSalt kg-1].
  real :: enthM_freezing    ! The enthalpy gain due to the mass gain by freezing [Q R Z ~> J m-2].
  real :: LI                ! The latent heat of fusion [Q ~> J kg-1].
  ! These variables are used only for debugging.
  real :: mtot_ice    ! The summed ice mass [R Z ~> kg m-2].
  real :: h2o_to_ocn, h2o_orig, h2o_imb ! Water budget terms [R Z ~> kg m-2].
  integer :: k, tr
  logical :: debug = .false.

  call get_SIS2_thermo_coefs(ITV, Latent_Fusion=LI)
  min_dEnth_freeze = LI * (1.0-CS%liq_lim)

  ! set mass mark; will subtract the mass at end to find the melt flux to ocean
  if (debug) then
    h2o_orig = 0.0 ; do k=0,NkIce ; h2o_orig = h2o_orig + m_lay(k) ; enddo
  endif

  h2o_ocn_to_ice = 0.0 ; salt_to_ice = 0.0 ; enthM_freezing = 0.0

  ! Assume that Salin(NkIce+1) already is the freezing salinity.
  salin_freeze = Salin(NkIce+1)

  ! Add frazil:
  !   Frazil mostly forms in leads, so add its heat uniformly over all of the
  ! layers rather than just adding it to the ice bottom.
  if (frazil > 0.0) then
    frazil_per_layer = frazil / NkIce
    do k=1,NkIce
    ! ### t_frazil and enth_frazil are calculated in a kludgey way here; revisit this?
      t_frazil = min(tfw, T_Freeze(sice_therm(k), ITV) - CS%Frazil_temp_offset)
      enth_frazil = min(enth_from_TS(t_frazil, sice_therm(k), ITV), &
                        enthalpy(NkIce+1) - min_dEnth_freeze)
      m_frazil = frazil_per_layer / (enthalpy(NkIce+1) - enth_frazil)

      Enthalpy(k) = (m_lay(k)*Enthalpy(k) + m_frazil*enth_frazil) / (m_lay(k) + m_frazil)
      Salin(k) = (m_lay(k)*Salin(k) + m_frazil*salin_freeze) / (m_lay(k) + m_frazil)
      do tr = 1,npassive
        ! Passive tracer assume that the flux of tracer from frazil is in NkIce+1
        TrLay(k,tr) = (m_lay(k)*TrLay(k,tr) + m_frazil*TrLay(NkIce+1,tr)) / (m_lay(k) + m_frazil)
      enddo
      Salt_to_ice = Salt_to_ice + m_frazil*salin_freeze

      m_lay(k) = m_lay(k) + m_frazil
      h2o_ocn_to_ice = h2o_ocn_to_ice + m_frazil

      enthM_freezing = enthM_freezing + m_frazil*enthalpy(NkIce+1)
    enddo
  endif

  ! There are no further heat or mass losses or gains by the ice+snow.
  enthalpy_freeze = enthM_freezing

  ! With the addition of frazil only, there is no need to make the snow below
  ! waterline adjustment, and no ice is converted to seawater.

  if (debug) then
    mtot_ice = 0.0 ; do k=1,NkIce ; mtot_ice = mtot_ice + m_lay(k) ; enddo
    h2o_to_ocn = h2o_orig - (m_lay(0) + mtot_ice)
    h2o_imb = h2o_to_ocn + h2o_ocn_to_ice
    if (abs(h2o_to_ocn + h2o_ocn_to_ice) > &
        max(1.0e-10*US%kg_m3_to_R*US%m_to_Z, 1.0e-12*h2o_orig, 1.0e-12*(abs(h2o_ocn_to_ice)))) then
      h2o_imb = h2o_to_ocn + h2o_ocn_to_ice
    endif

    call ice_check(m_lay(0), mtot_ice, enthalpy, Sice_therm, &
                   NkIce, "at end of add_frazil_SIS2", ITV, US=US)
  endif

end subroutine add_frazil_SIS2


!> Adjust the mass of the various ice layers to give the prescribed relative thicknesses.
subroutine rebalance_ice_layers(m_lay, mtot_ice, Enthalpy, Salin, NkIce, npassive, TrLay)
  real, dimension(0:NkIce),   intent(inout) :: m_lay !< The ice mass by layer [R Z ~> kg m-2].
  real,                       intent(out)   :: mtot_ice !< The summed ice mass [R Z ~> kg m-2].
  real, dimension(0:NkIce+1), intent(inout) :: Enthalpy !< Snow, ice, and ocean enthalpy by layer [Q ~> J kg-1].
  real, dimension(NkIce+1),   intent(inout) :: Salin  !< Conserved ice bulk salinity by layer [gSalt kg-1]
  integer,                    intent(in)    :: NkIce  !< The number of ice layers.
  integer,                    intent(in)    :: npassive !< Number of passive tracers
  real, dimension(0:NkIce+1,npassive), &
                              intent(inout) :: TrLay !< Passive tracer in the column layer [Conc]

  real, dimension(NkIce) :: mlay_new     ! New mass in a layer [R Z ~> kg m-2]
  real, dimension(NkIce) :: enth_ice_new ! New integrated enthalpy [R Z Q ~> kg m-2 Enth]
  real, dimension(NkIce) :: sal_ice_new  ! New integrated salinity [R Z gSalt kg-1 ~> gSalt m-2]
  real, dimension(NkIce,npassive) :: tr_ice_new ! New integrated tracer amount [R Z Conc ~> kg Conc m-2]
  real :: m_k1_to_k2   ! The being transferred from layer k1 to layer k2 [R Z ~> kg m-2]
  real :: m_ice_avg    ! The average mass per layer [R Z ~> kg m-2]
  integer :: k, k1, k2, kold, knew, tr

  mtot_ice = 0.0 ; do k=1,NkIce ; mtot_ice = mtot_ice + m_lay(k) ; enddo
  if (mtot_ice == 0.0) then ! There is no ice, so quit.
    return
  else
    do k=1,NkIce ; mlay_new(k) = 0.0 ; enth_ice_new(k) = 0.0 ; enddo
    do k=1,NkIce ; sal_ice_new(k) = 0.0 ; enddo
    do k=1,NkIce ; do tr = 1,npassive ; tr_ice_new(k,tr) = 0.0 ; enddo ; enddo
    m_ice_avg = mtot_ice / NkIce
    k1 = 1 ; k2 = 1
    do  ! Add ice from k1 to k2 to even up layer thicknesses.
      ! m_k1_to_k2 = min(m_ice_avg - mlay_new(k2), m_lay(k1))
      if ((mlay_new(k2) >= m_ice_avg) .and. (k2 < NkIce)) then ; k2 = k2+1
      elseif (m_lay(k1) <= 0.0) then ; k1 = k1+1
      elseif ((m_ice_avg - mlay_new(k2) > m_lay(k1)) .or. (k2 == NkIce)) then
        ! Move all remaining ice from k1 to k2.
        m_k1_to_k2 = m_lay(k1)
        enth_ice_new(k2) = enth_ice_new(k2) + m_k1_to_k2 * Enthalpy(k1)
        sal_ice_new(k2) = sal_ice_new(k2) + m_k1_to_k2 * Salin(k1)
        do tr = 1,npassive
          tr_ice_new(k2,tr) = tr_ice_new(k2,tr) + m_k1_to_k2 * TrLay(k1,tr)
        enddo
        mlay_new(k2) = mlay_new(k2) + m_k1_to_k2
        m_lay(k1) = 0.0 ! = m_lay(k1) - m_k1_to_k2
        k1 = k1+1
      else
        ! Move some of the ice from k1 to k2.
        m_k1_to_k2 = m_ice_avg - mlay_new(k2)
        enth_ice_new(k2) = enth_ice_new(k2) + m_k1_to_k2 * Enthalpy(k1)
        sal_ice_new(k2) = sal_ice_new(k2) + m_k1_to_k2 * Salin(k1)
        do tr = 1,npassive
          tr_ice_new(k2,tr) = tr_ice_new(k2,tr) + m_k1_to_k2 * TrLay(k1,tr)
        enddo
        mlay_new(k2) = m_ice_avg ! = mlay_new(k2) + m_k1_to_k2
        m_lay(k1) = m_lay(k1) - m_k1_to_k2
        k2 = k2+1
      endif
      if (k1 > NkIce) exit
    enddo
    do k=1,NkIce ; if (mlay_new(k) > 0.0) then
      Enthalpy(k) = enth_ice_new(k)/mlay_new(k)
      Salin(k) = sal_ice_new(k)/mlay_new(k)
      do tr = 1,npassive ; TrLay(k,tr) = tr_ice_new(k,tr)/mlay_new(k) ; enddo
    endif ; enddo
    do k=1,NkIce ; m_lay(k) = mlay_new(k) ; enddo
  endif

end subroutine rebalance_ice_layers

!> SIS2_ice_thm_end deallocates a SIS2_ice_thm_CS type and any sub-types.
subroutine SIS2_ice_thm_end(CS)
  type(SIS2_ice_thm_CS), pointer :: CS !< A pointer to the SIS2_ice_thm control structure.

  deallocate(CS)
end subroutine SIS2_ice_thm_end



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Everything above this point pertains to the updating of the ice state due to
!  themrmodyanmic forcing.  Everything after this point relates to more basic
!  calculations of sea ice thermodynamics.  They could be separated into two
!  modules, but by keeping them together compliers stand a better chance of
!  inlining various small subroutines and achieving better performance.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_thermo_init initializes the sea-ice ice thermodynamics parameter structure.
subroutine ice_thermo_init(param_file, ITV, US, init_EOS )

  type(param_file_type), intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ice_thermo_type), pointer    :: ITV !< A pointer to the ice thermodynamic parameter structure.
  type(unit_scale_type), intent(in) :: US  !< A structure with unit conversion factors
  logical,     optional, intent(in) :: init_EOS !< If true initialize the equation of state for seawater

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "SIS2_ice_thm (thermo)" ! This module's name.
  logical :: specified_ice

  if (.not.associated(ITV)) allocate(ITV)

  call log_version(param_file, mdl, version, &
     "This sub-module calculates ice thermodynamic quantities.")
  ITV%Q_to_J_kg = US%Q_to_J_kg
  call get_param(param_file, mdl, "LATENT_HEAT_FUSION", ITV%Lat_Fus, &
                 "The latent heat of fusion as used by SIS.", &
                 units="J kg-1", default=3.34e5, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "LATENT_HEAT_VAPOR", ITV%Lat_Vapor, &
                 "The latent heat of vaporization of water at 0C as used by SIS.", &
                 units="J kg-1", default=2.500e6, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "RHO_OCEAN", ITV%Rho_water, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RHO_ICE", ITV%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RHO_SNOW", ITV%Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "CP_ICE", ITV%Cp_ice, &
                 "The heat capacity of fresh ice, approximated as a constant.", &
                 units="J kg-1 K-1", default=2100.0, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "CP_SEAWATER", ITV%Cp_Water, &
                 "The heat capacity of sea water, approximated as a constant.", &
                 units="J kg-1 K-1", default=4200.0, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "CP_BRINE", ITV%Cp_brine, &
                 "The heat capacity of water in brine pockets within the "//&
                 "sea-ice, approximated as a constant.  CP_BRINE and "//&
                 "CP_SEAWATER should be equal, but for computational "//&
                 "convenience CP_BRINE can be set equal to CP_ICE.", &
                 units="J kg-1 K-1", default=US%Q_to_J_kg*ITV%Cp_Water, scale=US%J_kg_to_Q)
  call get_param(param_file, mdl, "DTFREEZE_DS", ITV%dTf_dS, &
                 "The derivative of the freezing temperature with salinity.", &
                 units="deg C PSU-1", default=-0.054)

  call get_param(param_file, mdl, "ENTHALPY_LIQUID_0", ITV%enth_liq_0, &
                 "The enthalpy of liquid fresh water at 0 C.  The solutions "//&
                 "should be physically consistent when this is adjusted, "//&
                 "because only the relative value is of physical meaning, "//&
                 "but roundoff errors can change the solution.", &
                 units="J kg-1", default=0.0, scale=US%J_kg_to_Q)
  call obsolete_real(param_file, "ENTHALPY_UNITS", warning_val=1.0)
  call get_param(param_file, mdl, "SUBLIMATION_BUG", ITV%sublimation_bug, &
                 "If true use an older calculation that omits the latent "//&
                 "heat of fusion from the latent heat of sublimation. "//&
                 "This variable should be obsoleted as soon as possible.", &
                 default=.false.)

  call get_param(param_file, mdl, "SPECIFIED_ICE", specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  if (specified_ice) then
    ITV%slab_ice = .true.
    call log_param(param_file, mdl, "USE_SLAB_ICE", ITV%slab_ice, &
                 "Use the very old slab-style ice.  With SPECIFIED_ICE, "//&
                 "USE_SLAB_ICE is always true.")
  else
    call get_param(param_file, mdl, "USE_SLAB_ICE", ITV%slab_ice, &
                 "If true, use the very old slab-style ice.", default=.false.)
  endif

  if (present(init_EOS)) then ; if (init_EOS) then
    if (.not.associated(ITV%EOS)) call EOS_init(param_file, ITV%EOS)
  endif ; endif

end subroutine ice_thermo_init

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> T_Freeze returns the freezing temperature as a function of salinity.
!### (and possibly later pressure).
function T_Freeze(S, ITV)
  real, intent(in) :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: T_Freeze !< The freezing point temperature [degC]

  T_Freeze = 0.0 + ITV%dTf_dS * S

end function T_Freeze

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> calculate_T_Freeze calculates an array of freezing temperatures for an
!!            an array of salinities (and maybe later pressures).
subroutine calculate_T_Freeze(S, T_Freeze, ITV)
  real, dimension(:), intent(in) :: S  !< The ice bulk salinity [gSalt kg-1]
  real, dimension(:), intent(out) :: T_Freeze !< The freezing point temperature [degC]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.

  integer :: k, nk_ice
  nk_ice = size(S)

  do k=1,nk_ice ; T_Freeze(k) = 0.0 + ITV%dTf_dS * S(k) ; enddo

end subroutine calculate_T_Freeze

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> enthalpy_from_TS sets a column of enthalpies from temperature and salinity.
subroutine enthalpy_from_TS(T, S, enthalpy, ITV)
  real, dimension(:), intent(in) :: T  !< The ice temperature [degC]
  real, dimension(:), intent(in) :: S  !< The ice bulk salinity [gSalt kg-1]
  real, dimension(:), intent(out) :: enthalpy !< The ice enthalpy, in enthalpy units [Q ~> J kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.

  integer :: k, nk_ice
  nk_ice = size(T)

  do k=1,nk_ice ; enthalpy(k) = enth_from_TS(T(k), S(k), ITV) ; enddo

end subroutine enthalpy_from_TS

!> enth_from_TS returns an ice enthalpy given temperature and salinity.
function enth_from_TS(T, S, ITV) result(enthalpy)
  real, intent(in)  :: T  !< The ice temperature [degC]
  real, intent(in)  :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: enthalpy !< The ice enthalpy, in enthalpy units [Q ~> J kg-1]

  real :: T_fr   ! The freezing temperature [degC].

  T_fr = ITV%dTf_dS*max(0.0,S)

  if ((S == 0.0) .and. (T <= 0.0)) then
    ! Note that at the freezing point, fresh water is assumed to be all ice,
    ! due to the degeneracy in inverting temperature for enthalpy.
    enthalpy = ((ITV%enth_liq_0 - ITV%Lat_Fus) + ITV%Cp_Ice*T)
  elseif (T >= T_fr) then ! This layer is already melted, so the enthalpy is
    ! just what is required to warm or cool it to 0 degC.
    enthalpy = (ITV%enth_liq_0 + ITV%Cp_Water*T)
  elseif (ITV%Cp_Ice == ITV%Cp_Brine) then
    enthalpy = ((ITV%enth_liq_0 - ITV%Lat_Fus * (1.0 - T_fr/T)) + &
                  (ITV%Cp_Ice*T + (ITV%Cp_Water - ITV%Cp_Ice)*T_fr))
  else  ! The derivation of this expression can be found in the SIS2 manual;
    ! it assumes that the freezing temperature varies linearly with salinity.
    enthalpy = ((ITV%enth_liq_0 - ITV%Lat_Fus * (1.0 - T_fr/T)) + &
                  ((ITV%Cp_Ice*T + (ITV%Cp_Water - ITV%Cp_Ice)*T_fr) - &
                   (ITV%Cp_Brine - ITV%Cp_Ice) * T_fr*log(T_fr/T))) ! Note that log(1/a) = -log(a).
  endif

end function enth_from_TS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> enthalpy_liquid_freeze returns the enthalpy of liquid water at the freezing
!!    point for a given salinity [Q ~> J kg-1].
function enthalpy_liquid_freeze(S, ITV)
  real, intent(in)  :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: enthalpy_liquid_freeze

  enthalpy_liquid_freeze = ITV%Cp_water*(ITV%dTf_dS*S) + ITV%enth_liq_0

end function enthalpy_liquid_freeze

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> enthalpy_liquid returns the enthalpy of liquid water at the given
!!     temperature and salinity, in enthalpy units [Q ~> J kg-1]
function enthalpy_liquid(T, S, ITV)
  real, intent(in)  :: T  !< The ice temperature [degC]
  real, intent(in)  :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: enthalpy_liquid !< The enthalpy of liquid water at the given temperature
                          !! and salinity [Q ~> J kg-1]

  enthalpy_liquid = ITV%enth_liq_0 + ITV%Cp_Water*T

end function enthalpy_liquid

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> enth_melt returns the enthalpy change associated with melting water of
!! a given temperature (T) and salinity (S), in enthalpy units [Q ~> J kg-1].
function enth_melt(T, S, ITV) result (emelt)
  real, intent(in)  :: T  !< The ice temperature [degC]
  real, intent(in)  :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real             :: emelt !< The enthaply change needed to melt the frozen water [Q ~> J kg-1]

  emelt = enthalpy_liquid_freeze(S, ITV) - enth_from_TS(T, S, ITV)
end function enth_melt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> latent_sublimation returns the latent heat of sublimation as vapor at 0 degC
!!   for a given ice and snow enthalpy (in [Q ~> J kg-1]), and a weighting factor
!!   for the sublimation between the snow and ice.
function latent_sublimation(enth_snow, enth_ice, wt_snow, ITV) result (latent)
  real, intent(in) :: enth_snow  !< The enthalpy of the snow [Q ~> J kg-1].
  real, intent(in) :: enth_ice   !< The enthalpy of the ice surface [Q ~> J kg-1].
  real, intent(in) :: wt_snow    !< A weighting factor (0-1) for the snow areal
                                 !! coverage [nondim]; the complement is for the ice.
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: latent !< The latent heat of sublimation [Q ~> J kg-1].

  real :: enth_liq_0 ! The value of enthalpy for liquid fresh water at 0 degC [Q ~> J kg-1].
                     ! This should become ITV%Enth_liq_0, but it is not due to
                     ! a bug in how this is calculated.

  if (ITV%sublimation_bug) then  ! This is false by default and should be obsoleted.
    enth_liq_0 = Enth_from_TS(0.0, 0.0, ITV)
  else
    enth_liq_0 = ITV%enth_liq_0
  endif

  if (ITV%slab_ice) then
    latent = ITV%Lat_Vapor + ITV%Lat_Fus
  else
    latent = ITV%Lat_Vapor + (enth_liq_0 - ((1.0 - wt_snow) * enth_ice + wt_snow * enth_snow))
  endif

end function latent_sublimation


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Temp_from_Enth_S sets a column of temperatures from enthalpy and salinity.
subroutine Temp_from_Enth_S(En, S, Temp, ITV)
  real, dimension(:), intent(in) :: En !< The ice enthalpy, in enthalpy units [Q ~> J kg-1]
  real, dimension(:), intent(in) :: S  !< The ice bulk salinity [gSalt kg-1]
  real, dimension(:), intent(out) :: Temp !< The ice temperature [degC]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.

  integer :: k, nk_ice

  nk_ice = size(Temp)

  do k=1,nk_ice
    Temp(k) = Temp_from_En_S(En(k), S(k), ITV)
  enddo

end subroutine Temp_from_Enth_S

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dTemp_dEnth_EnS returns the partial deriative of sea ice temperature with enthalpy
!! for ice of a given enthalpy and salinity.
function dTemp_dEnth_EnS(En, S, ITV) result(dT_dE)
  real, intent(in)  :: En !< The ice enthalpy, in enthalpy units [Q ~> J kg-1]
  real, intent(in)  :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: dT_dE  !< Partial derivative of temperature with enthalpy [degC Q-1 ~> degC kg J-1].

  real :: I_CpI_Eu, I_CpW_Eu ! Inverse heat capacities [degC Q-1 ~> degC kg J-1]
  real :: BB    ! A temporary variable [Q ~> J kg-1]
  real :: T_fr  ! The freezing temperature [degC].
  real :: En_J  ! Enthalpy with 0 offset at 0 degC [Q ~> J kg-1].

  I_CpI_Eu = 1.0 / (ITV%Cp_Ice)
  I_CpW_Eu = 1.0 / (ITV%Cp_Water)

  T_fr = ITV%dTf_dS*S

  En_J = (En - ITV%enth_liq_0)
  if (S <= 0.0) then ! There is a step function for fresh water.
    if (En_J >= 0.0) then ; dT_dE = I_CpW_Eu
    elseif (En_J > -ITV%Lat_Fus) then ; dT_dE = 0.0
    else ; dT_dE = I_CpI_Eu ; endif
  elseif (ITV%Cp_Ice == ITV%Cp_Brine) then
    ! This makes the assumption that all water in the ice and snow categories,
    ! both fluid and in pockets, has the same heat capacity.
    if (En_J < T_fr * ITV%Cp_water) then
      BB = 0.5*((En_J - T_fr*(ITV%Cp_water-ITV%Cp_ice)) + ITV%Lat_Fus)
      dT_dE = I_CpI_Eu * 0.5 * (1.0 - BB / sqrt(BB**2 - T_fr*ITV%Cp_Ice*ITV%Lat_Fus))
    else  ! This layer is already melted, so just warm it to 0 degC.
      dT_dE = I_CpW_Eu
    endif
  else
    call SIS_error(FATAL, "Write dTemp_dEnth_Enth for Cp_ice /= Cp_brine.")
  endif

end function dTemp_dEnth_EnS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> dTemp_dEnth_TS returns the partial deriative of sea ice temperature with enthalpy
!! for ice of a given temperature and salinity.
function dTemp_dEnth_TS(Temp, S, ITV) result(dT_dE)
  real, intent(in)  :: Temp !< The ice temperature [degC]
  real, intent(in)  :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: dT_dE  !< Partial derivative of temperature with enthalpy [degC Q-1 ~> degC kg J-1].

  real :: I_CpI_Eu, I_CpW_Eu ! Inverse heat capacities [degC Q-1 ~> degC kg J-1]
  real :: T_fr  ! The freezing temperature [degC].

  I_CpI_Eu = 1.0 / (ITV%Cp_Ice)
  I_CpW_Eu = 1.0 / (ITV%Cp_Water)

  T_fr = ITV%dTf_dS*S

  if (S <= 0.0) then ! There is a step function for fresh water.
    if (Temp > T_fr) then ; dT_dE = I_CpW_Eu
    elseif (Temp == T_Fr) then ; dT_dE = 0.0
    else ; dT_dE = I_CpI_Eu ; endif
  else
    ! This makes the assumption that all water in the ice and snow categories,
    ! both fluid and in pockets, has the same heat capacity.
    if (Temp < T_fr) then
        ! These are equivalent expressions.
        !  dEn_dT = ( -ITV%Lat_Fus * (T_fr / Temp**2)) + &
        !             ITV%Cp_Ice + (ITV%Cp_Brine - ITV%Cp_Ice) * (T_fr/Temp)
        dT_dE = (-Temp) / ((ITV%Lat_Fus * (T_fr / Temp) + &
                           (ITV%Cp_Ice*(-Temp) + (ITV%Cp_Brine - ITV%Cp_Ice) * (-T_Fr))))
    else  ! This layer is already melted, so just warm it to 0 degC.
      dT_dE = I_CpW_Eu
    endif
  endif

end function dTemp_dEnth_TS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Temp_from_En_S returns the sea ice temperature for ice of a given enthalpy
!! and salinity.
function Temp_from_En_S(En, S, ITV) result(Temp)
  real, intent(in)  :: En !< The ice enthalpy, in enthalpy units [Q ~> J kg-1]
  real, intent(in)  :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: Temp  !< Temperature [degC].

  real :: I_Cp_Ice, I_Cp_Water  ! Inverse heat capacities [degC Q ~> kg degC J-1].
  real :: BB    ! A temporary variable [Q ~> J kg-1]
  real :: T_fr  ! The freezing temperature [degC].
  real :: En_J  ! Enthalpy with 0 offset [Q ~> J kg-1].
  real :: T_guess  ! The latest best guess at Temp [degC].
  real :: T_deriv  ! The value of Temp at which to evaluate dT_dEn [degC].
  real :: T_next   ! The tentative next value for T_guess [degC].
  real :: T_max, T_min ! Bracketing temperatures [degC].
  real :: En_Tg    ! The enthalpy at T_guess [Q ~> J kg-1].
  real :: En_Tmin, En_Tmax ! The enthalpies at T_max and T_min [Q ~> J kg-1].
  real :: dT_dEn   ! The partial derivative of temperature with enthalpy [degC Q-1 ~> degC kg J-1].
  real :: Enth_tol = 1.0e-15 ! The fractional Enthalpy difference tolerance for convergence [nondim].

!  real :: dTemp(20), T_itt(20)
  integer :: itt

  I_Cp_Ice = 1.0 / ITV%Cp_Ice
  I_Cp_Water = 1.0 / ITV%CP_Water

  T_fr = ITV%dTf_dS*S
  En_J = En - ITV%enth_liq_0

  if (S <= 0.0) then ! There is a step function for fresh water.
    if (En_J >= 0.0) then ; Temp = En_J * I_Cp_Water
    elseif (En_J >= -ITV%Lat_Fus) then ; Temp = 0.0
    else ; Temp = I_Cp_Ice * (En_J + ITV%Lat_Fus) ; endif
  elseif (En_J >= T_fr*ITV%Cp_Water) then  ! This layer is completely melted.
    Temp = En_J * I_Cp_Water
  else
    !   This makes the assumption that all water in the ice and snow categories,
    ! both fluid and in pockets, has the same heat capacity. This may be the
    ! final solution or it may be a good first guess to start the iterations.
    !   ITV%Lat_Fus * (T_fr/T) + ITV%Cp_Ice*T = En_J - (ITV%Cp_Water-ITV%Cp_Ice)*T_fr + ITV%Lat_Fus
    !   ITV%Lat_Fus * (T_fr/T) + ITV%Cp_Ice*T = 2.0*BB
    !   ITV%Cp_Ice*T**2 - 2.0*BB*T + ITV%Lat_Fus * T_fr = 0.0

    ! Note that (En_J < T_fr*ITV%Cp_Water)
    BB = 0.5*((En_J - T_fr*(ITV%Cp_water-ITV%Cp_ice)) + ITV%Lat_Fus)
    Temp = I_Cp_Ice * (BB - sqrt(BB**2 - T_fr*ITV%Cp_Ice*ITV%Lat_Fus))

    if (ITV%Cp_Ice /= ITV%Cp_brine) then
      T_min = -273.15
      En_Tmin = ( - ITV%Lat_Fus * (1.0 - T_fr/T_min) ) + &
              ((ITV%Cp_Ice*T_min + (ITV%Cp_Water-ITV%Cp_Ice)*T_fr) - &
               (ITV%Cp_Brine - ITV%Cp_Ice) * (T_fr * log(T_fr/T_min)))
      ! Could En_Tmin be approximated as -ITV%Lat_Fus + (ITV%Cp_Ice*T_min + (ITV%Cp_Water-ITV%Cp_Ice)*T_fr) ?
      T_max = T_fr ; En_Tmax = (ITV%Cp_Water*T_fr)

      ! This might be a good enough first guess that bracketing is unnecessary.
      T_guess = Temp
!      dTemp(:) = 0.0 ; T_itt(:) = 0.0
      do itt=1,20 ! Note that 3 or 4 iterations usually are enough.
!        T_itt(itt) = T_guess
        ! This expression uses the fact that the freezing point varies linearly
        ! with salinity.
        En_Tg = ( - ITV%Lat_Fus * (1.0 - T_fr/T_guess)) + &
                ((ITV%Cp_Ice*T_guess + (ITV%Cp_Water-ITV%Cp_Ice)*T_fr) - &
                 (ITV%Cp_Brine - ITV%Cp_Ice) * (T_fr * log(T_fr/T_guess)))

        if (abs(En_Tg - En_J) <= 1.0e-15*(abs(En_J) + abs(En_Tg))) then
          Temp = T_guess ; exit ! (The exit could be a return?)
        elseif (En_Tg < En_J) then
          T_min = T_guess ; En_Tmin = En_Tg
        else
          T_max = T_guess ; En_Tmax = En_Tg
        endif

        ! Use the more efficient Newton's method of McDougall & Witherspoon (2014),
        ! Appl. Math. Lett., 29, 20-25.
        T_deriv = T_guess
        if (itt > 1) then ! Reuse the estimate of dT_dEn from the last iteration.
          T_deriv = T_guess + 0.5*dT_dEn * (En_J - En_Tg)
          if ((T_deriv < T_min) .or. (T_deriv > T_max)) T_deriv = T_guess
        endif

        ! These are equivalent expressions.
        !  dEn_dT = ( -ITV%Lat_Fus * (T_fr / T_deriv**2)) + &
        !             ITV%Cp_Ice + (ITV%Cp_Brine - ITV%Cp_Ice) * (T_fr/T_deriv)
        dT_dEn = (-T_deriv) / (ITV%Lat_Fus * (T_fr / T_deriv) + &
                  (ITV%Cp_Ice*(-T_deriv) + (ITV%Cp_Brine - ITV%Cp_Ice) * (-T_fr)))
        T_next = T_guess + dT_dEn * (En_J - En_Tg)

!        dTemp(itt) = T_next - T_guess
        if ((T_next > T_max) .or. (T_next < T_min)) then ! Use the false position method.
          T_guess = ((En_Tmax - En_J) * T_min + (En_J - En_Tmin) * T_max) / &
                     (En_Tmax - En_Tmin)
        else
          T_guess = T_next
        endif
      enddo
      Temp = T_guess

!   These were used to debug this routine.
!      if (itt > 5) then
!        write (*,'("dTemp = ",8E12.4)') dTemp(1:8)
!        write (*,'("T_itt = ",8E14.6)') T_itt(1:8)
!      endif
    endif
  endif

end function Temp_from_En_S

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> e_to_melt_TS - return the energy needed to melt a given snow/ice
!!      configuration [J kg-1].
function e_to_melt_TS(T, S, ITV)
  real, intent(in) :: T  !< The ice temperature [degC]
  real, intent(in) :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real :: e_to_melt_TS   !< The energy required to melt this mixture of ice and brine
                         !! and warm it to its bulk freezing temperature [J kg-1].

  real :: e_to_melt  ! The energy required to melt this mixture of ice and brine
                     ! and warm it to its bulk freezing temperature [Q ~> J kg-1].
  real :: T_fr       ! The freezing temperature [degC].
  T_fr = ITV%dTf_dS*S

  if (T >= T_Fr) then ! This layer is already melted and has excess heat.
    e_to_melt = ITV%Cp_Water * (-T + T_Fr)
  elseif (ITV%Cp_Ice == ITV%Cp_brine) then
    e_to_melt = (ITV%Lat_Fus - ITV%Cp_ice*T) * (1.0 - T_Fr/T)
            ! =  ITV%Lat_Fus * (1.0 - T_Fr/T) + ITV%Cp_ice * (T_Fr - T)
  else
    e_to_melt = ITV%Lat_Fus * (1.0 - T_fr/T) + (ITV%Cp_Ice * (T_Fr - T) + &
                   (ITV%Cp_Brine - ITV%Cp_Ice) * (T_fr*log(T_fr/T)))
  endif
  e_to_melt_TS = ITV%Q_to_J_kg*e_to_melt

end function e_to_melt_TS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> energy_melt_enthS returns the energy needed to melt a given snow/ice
!!      configuration [J kg-1].
function energy_melt_enthS(En, S, ITV) result(e_to_melt)
  real, intent(in) :: En !< The ice enthalpy, in enthalpy units [Q ~> J kg-1]
  real, intent(in) :: S  !< The ice bulk salinity [gSalt kg-1]
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.

  real :: e_to_melt  !< The energy required to melt this mixture of ice and brine
                     !! and warm it to its bulk freezing temperature [J kg-1].

  e_to_melt = ITV%Q_to_J_kg * (enthalpy_liquid_freeze(S, ITV) - En)

end function energy_melt_enthS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> get_SIS2_thermo_coefs returns various thermodynamic coefficients, rescaling the units
!! appropriately if an optional unit_scale_type argument is provided.
subroutine get_SIS2_thermo_coefs(ITV, ice_salinity, Cp_Ice, Cp_brine, Cp_water, &
                                 rho_ice, rho_snow, rho_water, Latent_fusion, Latent_vapor, &
                                 EOS, spec_thermo_salin, slab_ice)
  type(ice_thermo_type), intent(in) :: ITV !< The ice thermodynamic parameter structure.
  real, dimension(:), &
           optional, intent(out) :: ice_salinity  !< The specified salinity of each layer when the
                    !! thermodynamic salinities are pre-specified [gSalt kg-1].
  real,    optional, intent(out) :: Cp_Ice    !< The heat capacity of ice [Q degC-1 ~> J kg-1 degC-1].
  real,    optional, intent(out) :: Cp_Water  !< The heat capacity of seawater [Q degC-1 ~> J kg-1 degC-1].
  real,    optional, intent(out) :: Cp_Brine  !< The heat capacity of liquid water in brine pockets
                    !! within the sea-ice [Q degC-1 ~> J kg-1 degC-1].  Cp_Brine and Cp_Water should be equal,
                    !! but for computational convenience Cp_Brine has often been set equal to Cp_Ice instead.
  real,    optional, intent(out) :: rho_ice   !< A nominal density of ice [R ~> kg m-3] or [kg m-3].
  real,    optional, intent(out) :: rho_snow  !< A nominal density of snow [R ~> kg m-3] or [kg m-3].
  real,    optional, intent(out) :: rho_water !< A nominal density of water [R ~> kg m-3] or [kg m-3].
  real,    optional, intent(out) :: Latent_fusion !< The latent heat of fusion [Q ~> J kg-1].
  real,    optional, intent(out) :: Latent_vapor !< The latent heat of vaporization [Q ~> J kg-1].
  type(EOS_type), &
           optional, pointer     :: EOS       !< A pointer to the MOM6/SIS2 ocean equation-of-state type.
  logical, optional, intent(out) ::  spec_thermo_salin !< If true, all thermodynamic calculations
                    !! are done with a specified salinity profile that may be independent of the ice
                    !! bulk salinity.
  logical, optional, intent(out) ::  slab_ice !< If true, use the very old slab ice thermodynamics,
                    !! with effectively zero heat capacity of ice and snow.

  call get_thermo_coefs(ice_salinity=ice_salinity)

  if (present(Cp_Ice)) Cp_Ice = ITV%Cp_Ice
  if (present(Cp_Water)) Cp_Water = ITV%Cp_Water
  if (present(Cp_Brine)) Cp_Brine = ITV%Cp_Brine
  if (present(rho_ice)) rho_ice = ITV%rho_ice
  if (present(rho_snow)) rho_snow = ITV%rho_snow
  if (present(rho_water)) rho_water = ITV%rho_water
  if (present(Latent_fusion)) Latent_fusion = ITV%Lat_Fus
  if (present(Latent_vapor)) Latent_vapor = ITV%Lat_Vapor

  if (present(spec_thermo_salin)) spec_thermo_salin = .true.
  if (present(slab_ice)) slab_ice = ITV%slab_ice
  if (present(EOS)) then
    if (.not.associated(ITV%EOS)) call SIS_error(FATAL, &
      "An EOS pointer was requested via get_SIS2_thermo_coefs, but ITV%EOS "//&
      "has not yet been allocated.")
    EOS => ITV%EOS
  endif

end subroutine get_SIS2_thermo_coefs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_thermo_end deallocates an ice_thermo_type and any sub-types.
subroutine ice_thermo_end(ITV)
  type(ice_thermo_type), pointer :: ITV !< A pointer to the ice thermodynamic parameter structure.

  if (associated(ITV%EOS)) call EOS_end(ITV%EOS)
  deallocate(ITV)

end subroutine ice_thermo_end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!                       N-LAYER VERTICAL THERMODYNAMICS                        !
!                                                                              !
! References:                                                                  !
!   Hallberg, R., and M. Winton, 2016:  The SIS2.0 sea-ice model, in prep.     !
!                                                                              !
!   Winton, M., 2011:  A conservative non-iterative n-layer sea ice            !
!     temperature solver, in prep.                                             !
!                                                                              !
!                                                                              !
!         ->+---------+ <- ts - diagnostic surface temperature ( <= 0C )       !
!        /  |         |                                                        !
!      hs   |  snow   | <- tsn   One snow layer with heat capacity             !
!        \  |         |                                                        !
!         =>+---------+                                                        !
!        /  |         |                                                        !
!       /   |         | <- t1    N salty ice layers with heat capacity         !
!      /    |         |                                                        !
!     /     |         | <- t2                                                  !
!   hi      |...ice...|                                                        !
!     \     |         | <- tN-1                                                !
!      \    |         |                                                        !
!       \   |         | <- tN                                                  !
!        \  |         |                                                        !
!         ->+---------+ <- base of ice fixed at seawater freezing temp.        !
!                                                                              !
!                                       Bob Hallberg (Robert.Hallberg@noaa.gov)!
!                                       Mike Winton  (Michael.Winton@noaa.gov) !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Enhancement for melt ponds cohabiting ice top with snow layer:               !
!                                                                              !
!    ->+--------+ <- ts==0C when pond (net heat here -> pond freeze/melt)      !
!   /  |        |--------|<-+                                                  !
! hs   |  snow  |  pond  |   hp - layer has no heat capacity; "surface" energy !
!      |        |        |  .     balance takes place at ice/snow top          !
!   \  |        |        |  .                                                  !
!    =>+--------+--------|<-+                                                  !
!   /  |                 |                                                     !
! hi   |    ...ice...    |                   do_pond = true activates scheme   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
!             SIS2 "getting-started" melt pond scheme
!
! This melt pond scheme adds a single layer melt pond to each ice thickness
! category.  The layer does not have heat capacity.  It is assumed to be at
! freshwater freezing temperature and well mixed.  All pond surface fluxes
! are communicated directly to its bottom where surface energy balance is
! calculated.  The pond layer is advected and redistributed between ice
! thickness categories similarly to the snow and ice layers.  As is the case
! with snow, pond cannot exist without an ice layer below.  Basic descriptions
! of pond sources and sinks, and pond fraction for radiation treatments follow:
!
! Source water:  Surface melting of snow and ice are the source of pond water.
! The runoff scheme keys on the total ice covered area.  This is the only
! interaction between category variables in the scheme.  The scheme uses r, the
! fraction of melt (r)etained given by the CICE5 scheme (p. 42):
! r = r_min + (r_max-r_min)*ice_area.  The controlling namelist
! parameters for r_min and r_max are pond_r_min and pond_r_max.
!
! Freezing sink:  The surface energy budget continues to be performed at the
! top of the snow or ice when pond is present but the interface is fixed at
! freezing temperature and the residual between the fluxes on either side of the
! interface is made up with melt or freezing.  Freezing in this calculation is a
! sink of pond water.
!
! Freeboard sink:  if the pond and snow are sufficiently massive to push the top
! of the sea ice below sea level, pond is dumped into the ocean until the ice
! top is brought back to sea level, or the pond is completely depleted.  This
! adjustment follows the CICE5 Hunke "level ice" pond scheme.
!
! Porous-ice sink:  No through-ice drainage occurs until the ice average
! temperature exceeds a specified value.  The namelist parameter for this is
! pond_porous_temp.  Once this limit is exceeded the pond drains to a minimum
! value intended to represent coverage by ponds at sea level.  This scheme is
! a placeholder for the mushy-layer thermodynamics to be implemented later.
!
! Pond fraction for radiation:  In the snow-free case we assume that pond
! fraction ranges between a specified minimum, where surface cavities below
! sea level are filled with pond water, and a specified maximum value.
! The pond depth and pond fraction are considered to be proportional as was
! found at the SHEBA site and incorporated into the Bailey melt pond
! parameterization.  This means that the pond fraction is proportional to
! the square root of the pond volume.  The pond volume has a maximum
! determined by the non-negative freeboard requirement.  This pond volume
! is associated with the maximum pond fraction, so less pond water is needed
! to cover thinner ice.
!
! New namelist parameters: do_pond, pond_r_max, pond_r_min, pond_porous_temp,
! pond_frac_max, pond_frac_min
!
! New diagnostics: hp, fp, pond_source, pond_sink_freeboard, pond_sink_porous,
! pond_sink_tot <only hp implemented so far; add pond transport diagnostics?>
!
! M. Winton (6/16)
!
end module SIS2_ice_thm
