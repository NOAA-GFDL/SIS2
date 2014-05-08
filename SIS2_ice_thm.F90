!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!                       N-LAYER VERTICAL THERMODYNAMICS                        !
!                                                                              !
! Reference:                                                                   !
!   Winton, M., 2011:  A conservative non-iterative n-layer sea ice            !
!   temperature solver, in prep.                                               !
!                                                                              !
!                                                                              !
!         ->+---------+ <- ts - diagnostic surface temperature ( <= 0C )       !
!        /  |         |                                                        !
!      hs   |  snow   | <- tsn   One snow layer with heat capacity             !
!        \  |         |                                                        !
!         =>+---------+                                                        !
!        /  |         |                                                        !
!       /   |         | <- t1    Four salty ice layers with heat capacity      !
!      /    |         |                                                        !
!     /     |         | <- t2                                                  !
!   hi      |...ice...|                                                        !
!     \     |         | <- tN-1                                                !
!      \    |         |                                                        !
!       \   |         | <- tN                                                  !
!        \  |         |                                                        !
!         ->+---------+ <- base of ice fixed at seawater freezing temp.        !
!                                                                              !
!                                         Mike Winton (Michael.Winton@noaa.gov)!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module SIS2_ice_thm

! for calling delta-Eddington shortwave from ice_optics
use ice_shortwave_dEdd, only : shortwave_dEdd0_set_snow, shortwave_dEdd0_set_pond, &
                              shortwave_dEdd0, dbl_kind, int_kind, nilyr, nslyr
use ice_thm_mod, only : DS, DI, DW, TFI, MU_TS, CI, get_thermo_coefs

use constants_mod, only : LI => hlf ! latent heat of fusion - 334e3 J/(kg-ice)

implicit none ; private

public :: DS, DI, DW, MU_TS, CI, get_thermo_coefs
public :: SIS2_ice_thm_param, e_to_melt_TS, ice_optics_SIS2
public :: ice_temp_SIS2, ice_resize_SIS2, ecolumn_SIS2


!
! properties of ice, snow, and seawater (NCAR CSM values)
!

real            :: KS    = 0.31      ! conductivity of snow - 0.31 W/(mK)
real, parameter :: KI    = 2.03      ! conductivity of ice  - 2.03 W/(mK)
real, parameter :: CW    = 4.2e3     ! heat capacity of seawater 4200 J/(kg K)

! albedos are from CSIM4 assumming 0.53 visible and 0.47 near-ir insolation
real            :: ALB_SNO = 0.85       ! albedo of snow (not melting)
real            :: ALB_ICE = 0.5826     ! albedo of ice (not melting)
real            :: PEN_ICE = 0.3        ! ice surface penetrating solar fraction
real            :: OPT_DEP_ICE = 0.67   ! ice optical depth (m)
real            :: T_RANGE_MELT = 1.0   ! melt albedos scaled in below melting T

real            :: H_LO_LIM = 0.0       ! hi/hs lower limit for temp. calc.
real            :: H_SUBROUNDOFF = 1e-35 ! A miniscule value compared with H_LO_LIM
real            :: FRAZIL_TEMP_OFFSET = 0.5 ! A temperature offset between the 

integer, parameter :: NN = 4
real               :: DT = 1800.0

!
! In the ice temperature calculation we place a limit to below (salinity
! dependent) freezing point on the prognosed temperatures.  For ice_resize
! it is better to make a slightly more restrictive limit that requires the
! temperature to be such that the brine content is less than "liq_lim" of
! the total mass.  That is T_f/T < liq_lim implying T<T_f/liq_lim
real, parameter :: liq_lim = .99

! for calling delta-Eddington shortwave from ice_optics
logical :: do_deltaEdd = .true.
integer (kind=int_kind) :: &
   nx_block, ny_block, & ! block dimensions
   icells                ! number of ice-covered grid cells

integer (kind=int_kind), dimension (1) :: &
   indxi   , & ! compressed indices for ice-covered cells
   indxj

! inputs
real (kind=dbl_kind), dimension (1,1) :: &
   aice   , & ! concentration of ice
   vice   , & ! volume of ice
   vsno   , & ! volume of snow
   Tsfc   , & ! surface temperature
   coszen , & ! cosine of solar zenith angle
   tarea  , & ! cell area - not used
   swvdr  , & ! sw down, visible, direct  (W/m^2)
   swvdf  , & ! sw down, visible, diffuse (W/m^2)
   swidr  , & ! sw down, near IR, direct  (W/m^2)
   swidf      ! sw down, near IR, diffuse (W/m^2)


! outputs
real (kind=dbl_kind), dimension (1,1) :: &
   fs     , & ! horizontal coverage of snow
   fp     , & ! pond fractional coverage (0 to 1)
   hp         ! pond depth (m)

real (kind=dbl_kind), dimension (1,1,1) :: &
   rhosnw , & ! density in snow layer (kg/m3)
   rsnw       ! grain radius in snow layer (micro-meters)

real (kind=dbl_kind), dimension (1,1,18) :: &
         trcr        ! aerosol tracers


real (kind=dbl_kind), dimension (1,1) :: &
   alvdr   , & ! visible, direct, albedo (fraction) 
   alvdf   , & ! visible, diffuse, albedo (fraction) 
   alidr   , & ! near-ir, direct, albedo (fraction) 
   alidf   , & ! near-ir, diffuse, albedo (fraction) 
   fswsfc  , & ! SW absorbed at snow/bare ice/pondedi ice surface (W m-2)
   fswint  , & ! SW interior absorption (below surface, above ocean,W m-2)
   fswthru     ! SW through snow/bare ice/ponded ice into ocean (W m-2)

real (kind=dbl_kind), dimension (1,1,1) :: &
   Sswabs      ! SW absorbed in snow layer (W m-2)

real (kind=dbl_kind), dimension (1,1,nilyr) :: &
   Iswabs      ! SW absorbed in ice layer (W m-2)

real (kind=dbl_kind), dimension (1,1) :: &
   albice  , & ! bare ice albedo, for history  
   albsno  , & ! snow albedo, for history  
   albpnd      ! pond albedo, for history  

contains

!
! energy needed to melt kg of ice with given temp/salinity
!
function emelt(temp, salt) result (retval)
real :: retval
real :: temp, salt

  real :: tfi

  tfi = -MU_TS*salt
  if (tfi == 0.0) then
    retval = CI*(tfi-temp)+LI
  else
    retval = CI*(tfi-temp)+LI*(1-tfi/temp)
  endif
!
!Niki: The above formulation is not valid when temp>tfi and is causing problems,e.g. hlay(k)<0
!      The following is a trial to fix the issue by changing the enthalpy reference point.
!
!  if (tfi <= temp) then
!    retval = -CW*temp     !all ice already melted. bring saline water to 0
!  elseif (tfi == 0.0) then
!    retval = CI*(tfi-temp) + LI
!  else
!    retval = CI*(tfi-temp) + LI*(1-tfi/temp) - CW*tfi  !bring ice to tfi, melt the unmelted portion then bring saline water to 0
!  endif

end function emelt

!
! convert enthalpy back to temp. (inverse of emelt)
!
function emelt2temp(emelt, salt) result (retval)
real :: retval
real :: emelt, salt

  real :: tfi, A, B, C

  tfi = -MU_TS*salt
  A = CI
  B = emelt - LI - CI*tfi !!!Niki: for modified emelt +CW*tfi
  C = LI*tfi

!!!!  if (emelt < -CW*tfi) then ; retval = -emelt / CW
  if ( tfi == 0.0 ) then
    retval = -B/A
  else
    retval = -(B+sqrt(B*B-4*A*C))/(2*A)
  endif

end function emelt2temp

!
! snow/ice column energy (freezing temp. liquid water is "0")
!
function ecolumn_SIS2(hsno, hice, tice, sice, NkIce) result (retval)
  real :: retval
  real,    intent(in) :: hsno, hice
  integer, intent(in) :: NkIce
  real, dimension(NkIce), intent(in) :: tice, sice

  integer :: k

  retval = 0.0
  if ( hsno>0.0) then
    retval = -DS*hsno*LI
  endif

  do k=1,NkIce
    retval = retval-DI*(hice/NN)*emelt(tice(k), sice(k))
  enddo

end function ecolumn_SIS2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_thm_param - set ice thermodynamic parameters                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS2_ice_thm_param(alb_sno_in, alb_ice_in, pen_ice_in, opt_dep_ice_in, &
                         t_range_melt_in, ks_in, h_lo_lim_in, deltaEdd )
  real, intent(in)    :: alb_sno_in, alb_ice_in, pen_ice_in 
  real, intent(in)    :: opt_dep_ice_in, t_range_melt_in
  real, intent(in)    :: ks_in
  real, intent(in)    :: h_lo_lim_in
  logical, intent(in) :: deltaEdd 

  ALB_SNO     = alb_sno_in
  ALB_ICE     = alb_ice_in
  PEN_ICE     = pen_ice_in
  OPT_DEP_ICE = opt_dep_ice_in
  T_RANGE_MELT = t_range_melt_in

  KS          = ks_in
  H_LO_LIM    = h_lo_lim_in
  H_SUBROUNDOFF = max(1e-35, 1e-20*H_LO_LIM)
  do_deltaEdd = deltaEdd

end subroutine SIS2_ice_thm_param

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_optics - set albedo, penetrating solar, and ice/snow transmissivity      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_optics_SIS2(hs, hi, ts, tfw, alb_vis_dir, alb_vis_dif, alb_nir_dir, &
     alb_nir_dif, abs_sfc, abs_snow, abs_ice1, abs_ice2, &
     abs_ice3, abs_ice4, abs_ocn, abs_int, pen, trn, coszen_in)
  real, intent(in   ) :: hs  ! snow thickness (m-snow)
  real, intent(in   ) :: hi  ! ice thickness (m-ice)
  real, intent(in   ) :: ts  ! surface temperature
  real, intent(in   ) :: tfw ! seawater freezing temperature
  real, intent(  out) :: alb_vis_dir ! ice surface albedo (0-1)
  real, intent(  out) :: alb_vis_dif ! ice surface albedo (0-1)
  real, intent(  out) :: alb_nir_dir ! ice surface albedo (0-1)
  real, intent(  out) :: alb_nir_dif ! ice surface albedo (0-1)
  real, intent(  out) :: abs_sfc  ! frac abs sw abs at surface
  real, intent(  out) :: abs_snow ! frac abs sw abs in snow
  real, intent(  out) :: abs_ice1 ! frac abs sw abs at ice layer 1
  real, intent(  out) :: abs_ice2 ! frac abs sw abs at ice layer 2
  real, intent(  out) :: abs_ice3 ! frac abs sw abs at ice layer 3
  real, intent(  out) :: abs_ice4 ! frac abs sw abs at ice layer 4
  real, intent(  out) :: abs_ocn  ! frac abs sw abs in ocean
  real, intent(  out) :: abs_int  ! frac abs sw abs in ice interior
  real, intent(  out) :: pen      ! frac     sw passed below the surface (frac 1-pen absrobed at the surface)
  real, intent(  out) :: trn      ! frac     sw passed below the bottom  (frac 1-trn absorbed at the interior)
  real, intent(in),optional :: coszen_in
  real :: alb, as, ai, cs
  real :: thick_ice_alb, tcrit, fh

  if(do_deltaEdd) then

     ! temporary for delta-Eddington shortwave call
     nx_block = 1
     ny_block = 1
     icells = 1
     indxi(1) = 1
     indxj(1) = 1
     aice(1,1) = 1.0
     tarea(1,1) = 1.0 ! not used

     ! stuff that matters
     coszen(1,1) = cos(3.14*67.0/180.0) ! NP summer solstice
     if(present(coszen_in))  coszen(1,1) = max(0.01,coszen_in)
     Tsfc(1,1) = ts
     vsno(1,1) = hs
     vice(1,1) = hi
     swvdr(1,1) = 0.25
     swvdf(1,1) = 0.25
     swidr(1,1) = 0.25
     swidf(1,1) = 0.25

     call shortwave_dEdd0_set_snow(nx_block, ny_block, &
          icells,             &
          indxi,    indxj,    &
          aice,     vsno,     &
          Tsfc,     fs,       &
          rhosnw,   rsnw) ! out: fs, rhosnw, rsnw

     call shortwave_dEdd0_set_pond(nx_block, ny_block, &
          icells,             &
          indxi,    indxj,    &
          aice,     Tsfc,     &
          fs,       fp,       &
          hp) ! out: fp, hp
     call shortwave_dEdd0  (nx_block, ny_block,    &
          icells,   indxi,       &
          indxj,    coszen,      &
          aice,     vice,        &
          vsno,     fs,          &
          rhosnw,   rsnw,        &
          fp,       hp,          &
          swvdr,    swvdf,       &
          swidr,    swidf,       &
          alvdf,    alvdr,       & ! out: these and below
          alidr,    alidf,       &
          fswsfc,   fswint,      &
          fswthru,  Sswabs,      &
          Iswabs,   albice,      &
          albsno,   albpnd)

     ! ### ADD PARENTHESES AND MULTIPLY BY A RECIPROCAL.
     alb = 1-fswsfc(1,1)-fswint(1,1)-fswthru(1,1)
     abs_sfc  = fswsfc(1,1)  /(1-alb)
     abs_snow = Sswabs(1,1,1)/(1-alb)
     abs_ice1 = Iswabs(1,1,1)/(1-alb)
     abs_ice2 = Iswabs(1,1,2)/(1-alb)
     abs_ice3 = Iswabs(1,1,3)/(1-alb)
     abs_ice4 = Iswabs(1,1,4)/(1-alb)
     abs_ocn  = fswthru(1,1) /(1-alb)

     alb_vis_dir = alvdr(1,1)
     alb_vis_dif = alvdf(1,1)
     alb_nir_dir = alidr(1,1)
     alb_nir_dif = alidf(1,1)

     ! pen = (fswint(1,1)+fswthru(1,1))/(fswsfc(1,1)+fswint(1,1)+fswthru(1,1))
     pen = abs_snow + abs_ice1 + abs_ice2 + abs_ice3 + abs_ice4 + abs_ocn
     ! trn = fswthru(1,1)/(fswint(1,1)+fswthru(1,1))
     trn = 0.0
     if(pen > 0.0) trn = abs_ocn/pen
     abs_int = 1.0 - trn

  else
     as = ALB_SNO; ai = ALB_ICE
     cs = hs/(hs+0.02)                        ! thin snow partially covers ice

     fh = min(atan(5.0*hi)/atan(5.0*0.5),1.0) ! use this form from CSIM4 to
     ! reduce albedo for thin ice
      if (ts+T_RANGE_MELT > TFI) then        ! reduce albedo for melting as in
         ! CSIM4 assuming 0.53/0.47 vis/ir
         as = as-0.1235*min((ts+T_RANGE_MELT-TFI)/T_RANGE_MELT,1.0)
         ai = ai-0.075 *min((ts+T_RANGE_MELT-TFI)/T_RANGE_MELT,1.0)
      endif
      ai = fh*ai+(1-fh)*0.06                 ! reduce albedo for thin ice

     alb = cs*as+(1-cs)*ai
     pen = (1-cs)*PEN_ICE
     trn = exp(-hi/OPT_DEP_ICE);
     alb_vis_dir = alb
     alb_vis_dif = alb
     alb_nir_dir = alb
     alb_nir_dif = alb

     !! check for ice albdeos out of range (0 to 1)
     ! if (alb.lt.0.0 .or. alb.gt.1.0) then
     !    print *,'ice_optics: albedo out of range, alb=',alb
     !    print *,'cs=',cs,  'as=',as, 'ai=',ai
     !    print *,'ts=',ts,  'fh=',fh, 'hs=',hs, 'hi=',hi, 'tfw=',tfw
     !    print *,'ALB_SNO=',ALB_SNO,  'ALB_ICE=',ALB_ICE, 'T_RANGE_MELT,=',T_RANGE_MELT, 'TFI=',TFI
     !    stop
     ! end if
  endif

end subroutine ice_optics_SIS2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice5lay_temp - A subroutine that calculates the snow and ice enthalpy        !
!    changes due to surface forcing.                                           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_temp_SIS2(hsno, tsn, hice, tice, sice, sh_T0, B, &
                         sol, tfw, fb, tsurf, dtt, NkIce, tmelt, bmelt)

  real, intent(in   ) :: hsno  ! snow thickness (m)
  real, intent(inout) :: tsn   ! snow temperature (deg-C)
  real, intent(in   ) :: hice  ! ice thickness (m)
  real, dimension(NkIce), &
        intent(inout) :: tice ! ice temperature by layer (deg-C)
  real, dimension(NkIce), &
        intent(in)    :: Sice ! ice salinity by layer (g/kg)
  real, intent(in   ) :: sh_T0 ! net surface heat flux (+ up) at ts=0 (W/m^2)
  real, intent(in   ) :: B     ! d(sfc heat flux)/d(ts) [W/(m^2 deg-C)]
  real, dimension(0:NkIce), &
        intent(in)    :: sol   ! Solar heating of the snow and ice layers (W m-2)
  real, intent(in   ) :: tfw   ! seawater freezing temperature (deg-C)
  real, intent(in   ) :: fb    ! heat flux from ocean to ice bottom (W/m^2)
  real, intent(  out) :: tsurf ! surface temperature (deg-C)
  real, intent(in   ) :: dtt   ! timestep (sec)
  integer, intent(in   ) :: NkIce ! The number of ice layers.
  real, intent(inout) :: tmelt ! accumulated top melting energy  (J/m^2)
  real, intent(inout) :: bmelt ! accumulated bottom melting energy (J/m^2)
!
! variables for temperature calculation [see Winton (1999) section II.A.]
! note:  here equations are multiplied by hi to improve thin ice accuracy
!
  real :: A ! Net downward surface heat flux from the atmosphere at 0C (W/m^2)
  real, dimension(0:NkIce) :: temp_est
  real, dimension(NkIce) :: tfi, tice_est ! estimated new ice temperatures
  real :: mi, ms, e_extra
  real :: kk, k10, k0a, k0skin
  real :: I_bb, b_denom_1
  real :: comp_rat ! The complement of rat, going from 0 to 1.
  real :: tsf, k0a_x_ta, tsno_est, hie, salt_part, rat, tsurf_est
  real, dimension(0:NkIce+1) :: cc ! Interfacial coupling coefficients.
  real, dimension(0:NkIce) :: bb   ! Effective layer heat capacities.
  real, dimension(0:NkIce) :: cc_bb ! Remaining coupling ratios.
  real :: hsnow_eff ! , Ks_h
  integer :: k

  A = -sh_T0
  DT = dtt ! set timestep from argument - awkward, remove later

!   call ice_temp(A, B, sol, tsurf, hsno, tsn, hice, tice, sice, tfw, fb, tmelt, bmelt)
  mi = DI*hice/NkIce           ! full ice layer mass
  ms = DS*hsno              ! full snow layer mass
  tfi(:) = -MU_TS*sice(:)   ! freezing temperature of ice layers
  hie = max(hice, H_LO_LIM); ! prevent thin ice inaccuracy (mw)
  kk = NkIce*KI/hie                     ! full ice layer conductivity

  tsf = tfi(1)                          ! surface freezing temperature
  if (hsno>0.0) tsf = 0.0
  hsnow_eff = hsno + H_SUBROUNDOFF
  
  k10 = 2.0*(KS*(NkIce*KI)) / (hice*KS + hsnow_eff*(NkIce*KI)) ! coupling ice layer 1 to snow
  k0a = (KS*B) / (0.5*B*hsnow_eff + KS)      ! coupling snow to "air"
  k0skin = 2.0*KS / hsnow_eff
  k0a_x_ta = (KS*A) / (0.5*B*hsnow_eff + KS) ! coupling times "air" temperture

  ! First get non-conservative estimate with implicit treatment of layer coupling.
  
  ! This is the start of what was ice_temp_est.

  ! Determine the effective layer heat capacities.
  bb(0) = hsno*DS*CI
  do k=1,NkIce   ! load bb with heat capacity term.
    salt_part = 0.0
    if (sice(k)>0.0) salt_part = -MU_TS*sice(k)*LI/(tice(k)*tice(k))
    bb(k) = (hice/NkIce)*DI*(CI-salt_part) ! add coupling to this later
  enddo

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
  temp_est(NkIce) = ( (sol(NkIce)*dtt + bb(NkIce)*tice(NkIce)) + &
                      cc(NkIce+1)*tfw ) * I_bb
  comp_rat = b_denom_1 * I_bb
  cc_bb(NkIce) = cc(NkIce) * I_bb

  do k=NkIce-1,1,-1
    b_denom_1 = bb(k) + comp_rat*cc(k+1)
    I_bb =  1.0 / (b_denom_1 + cc(k))
    temp_est(k) = ((sol(k)*dtt + bb(k)*tice(k)) + cc(k+1)*temp_est(k+1)) * I_bb

    comp_rat = b_denom_1 * I_bb ! 1.0 >= comp_rat >= 0.0
    cc_bb(k) = cc(k) * I_bb
  end do

  b_denom_1 = bb(0) + comp_rat*cc(1)
  I_bb =  1.0 / (b_denom_1 + cc(0))
  !   This is a complete calculation of temp_est(0), assuming that the surface
  ! flux is given by A - B*tsurf, with tsurf as estimated
  temp_est(0) = ((bb(0)*tsn + k0a_x_ta*dtt) + cc(1)*temp_est(1)) * I_bb
  ! Note that the SIS1 code omits sol(0) at its equivalent of this point.
  temp_est(0) = temp_est(0) + (sol(0)*dtt) * I_bb
  
  ! Diagnose the surface skin temperature by matching the diffusive fluxes in
  ! the snow with the atmospheric fluxes.  I.e. solve the following for tsurf_est:
  !  (A - B*tsurf_est) = k0skin * (tsurf_est - temp_est(0))
  tsurf_est = (A + k0skin*temp_est(0)) / (B + k0skin)

  if (tsurf_est > tsf) then
    ! The surface is melting, set tsurf to melt temp. and recalculate I_bb.
    tsurf_est = tsf
    ! cc(0) = k0skin*dtt

    I_bb =  1.0 / (b_denom_1 + (k0skin*dtt))
    temp_est(0) = ((bb(0)*tsn + k0skin*dtt*tsf) + cc(1)*temp_est(1)) * I_bb
    ! Note that the SIS1 code omits sol(0) at its equivalent of this point.
    temp_est(0) = temp_est(0) + (sol(0)*dtt) * I_bb
  endif
  ! Go back DOWN the ice column to get the final temperatures.
  do k=1,NkIce
    temp_est(k) = temp_est(k) + cc_bb(k) * temp_est(k-1)
  end do
  ! This is the end of what was ice_temp_est.

  tice_est(:) = temp_est(1:NkIce)
  tsno_est = temp_est(0)
!
! The following two lines disable implicit coupling estimate;  in a 10 year CORE
! test the skin temperature without implicit estimate was mostly within .05K
! but was cooler around topography in NH (up to 1K - Baffin Island).  Thickness
! looked very similar
!
! tice_est(:) = tice(:) ; tsno_est=tsn

  !
  ! conservative pass going UP the ice column
  !
  tice_est(NkIce) = laytemp_SIS2(mi, tfi(NkIce), sol(NkIce) + kk*(2*tfw+tice_est(NkIce-1)), &
                                 3*kk, tice(NkIce), dtt)
  do k=NkIce-1,2,-1
    tice_est(k) = laytemp_SIS2(mi, tfi(k), sol(k) + kk*(tice_est(k-1)+tice_est(k+1)), &
                               2*kk, tice(k), dtt)
  enddo
  tice_est(1) = laytemp_SIS2(mi, tfi(1), sol(1) + (kk*tice_est(2) + k10*tsno_est), &
                             kk + k10, tice(1), dtt)
  tsno_est = laytemp_SIS2(ms, 0.0, sol(0) + (k10*tice_est(1)+k0a_x_ta), k10+k0a, tsn, dtt)
  tsurf = (A + k0skin*tsno_est) / (B + k0skin)  ! diagnose surface skin temp.

  if (tsurf > tsf) then ! surface is melting, redo with surf. at melt temp.
    tsurf = tsf
    tsno_est = laytemp_SIS2(ms, 0.0, sol(0) + k10*tice_est(1) + k0skin*tsf,&
                            k10+k0skin, tsn, dtt)
    ! add in surf. melt
    if (hsno>0.0) then
      tmelt = tmelt+DT*((A-B*tsurf)-2*KS*(tsurf-tsno_est)/hsno)
    else
      tmelt = tmelt+DT*((sol(0)+A-B*tsurf)-k10*(tsurf-tice_est(1))) ! tsno = tsurf
    endif
  endif
  tsn = tsno_est ! finalize snow temperature

  !
  ! conservative pass going DOWN the ice column
  !
  tice(1)  = laytemp_SIS2(mi, tfi(1), sol(1)+kk*tice_est(1+1) &
                          + k10*(tsno_est-tice_est(1)), kk, tice(1), dtt)
  do k=2,NkIce-1 ! flux from above is fixed, only have downward feedback
    tice(k) = laytemp_SIS2(mi, tfi(k), sol(k)+kk*tice_est(k+1) &
                                       +kk*(tice(k-1)-tice_est(k)), kk, tice(k), dtt)
  enddo
  tice(NkIce) = laytemp_SIS2(mi, tfi(NkIce), sol(NkIce)+2*kk*tfw &
                             +kk*(tice(NkIce-1)-tice_est(NkIce)), 2*kk, tice(NkIce), dtt)
  !
  ! END conservative update
  !
  
  ! The following is a dangerous calculation if kk is too large, in part because
  ! all of the other fluxes are calculated implicitly for the layer above, so
  ! the temperatures are well bounded.  Can this be rearranged or calculated
  ! as a residual of the heat changes in the ice and snow?
  bmelt = bmelt + DT*(fb - 2*kk*(tfw-tice(NkIce))) ! add in bottom melting/freezing

  if (tsn > 0.0) then ! put excess snow energy into top melt
    e_extra = CI*DS*hsno*tsn
    tmelt = tmelt + e_extra
    tsn = 0.0
  endif

  do k=1,NkIce ; if (tice(k)>tfi(k)/liq_lim) then ! push excess energy to closer of top or bottom melt
    e_extra = (emelt(tfi(k)/liq_lim,sice(k))-emelt(tice(k),sice(k)))*DI*hie/NkIce
    tice(k) = tfi(k)/liq_lim
    if (k<=NkIce/2) then
      tmelt = tmelt+e_extra
    else
      bmelt = bmelt+e_extra
    endif
  endif ; enddo

  call temp_check(tsurf, hsno, tsn, hice, tice, NkIce, bmelt, tmelt)

end subroutine ice_temp_SIS2

!
! laytemp_SIS2 - implicit calculation of new layer temperature
!
function laytemp_SIS2(m, tfi, f, b, tp, dtt) result (retval)
  real retval
  real, intent(in) :: m    ! mass of ice - kg/m2
  real, intent(in) :: tfi  ! ice freezing temp. (determined by salinity)
  real, intent(in) :: f    ! Inward forcing - W/m2
  real, intent(in) :: b    ! response of outward heat flux to local temperature - W/m2/K
  real, intent(in) :: tp   ! prior step temperature
  real, intent(in) :: dtt  ! timestep in s.

  real :: AA, BB, CC
  ! solve quadratic equation for layer temperature, t:
  !
  !   m * {CI-LI*tfi/(t*tp)} * (t-tp) = dtt * (f - b*t)
  !

  if ( tfi == 0.0 ) then
    retval = (m*CI*tp + f*dtt) / (m*CI + b*dtt) ! = -BB/AA
  else
    AA = m*CI + b*dtt
    BB = -(m*LI*tfi/tp + m*CI*tp + f*dtt)
    CC = m*LI*tfi
    ! This form avoids round-off errors.
    if (BB >= 0) then
      retval = -(BB + sqrt(BB*BB - 4*AA*CC)) / (2*AA)
    else
      retval = (2*CC) / (-BB + sqrt(BB*BB - 4*AA*CC))
    endif
  endif

end function laytemp_SIS2

subroutine temp_check(ts, hs, tsn, hi, t_ice, NkIce, bmelt, tmelt)
  real, intent(in) :: ts, hs, tsn, hi, bmelt, tmelt
  real, dimension(NkIce), intent(in) :: t_ice
  integer, intent(in) :: NkIce
  integer :: k, bad

  bad = 0
  if (ts >0.0.or.ts <-100.0) bad = bad+1
  if (tsn>0.0.or.tsn<-100.0) bad = bad+1
  do k=1,NkIce ; if (t_ice(k) >0.0 .or. t_ice(k) < -100.0) bad = bad+1 ; enddo

  if (bad>0) then
    print *, 'BAD ICE AFTER TEMP ', 'hs/hi=',hs,hi,'ts/tsn/tice=',ts, &
                      tsn,t_ice(:),'tmelt/bmelt=',tmelt,bmelt
  endif
end subroutine temp_check

subroutine resize_check(hs, tsn, hi, t_ice, NkIce, bmelt, tmelt)
  real, intent(in) :: hs, tsn, hi, bmelt, tmelt
  real, dimension(NkIce), intent(in) :: t_ice
  integer, intent(in) :: NkIce
  integer :: k, bad

  bad = 0
  if (hs <0.0.or.hs > 1e3  ) bad = bad+1
  if (hi <0.0.or.hi > 1e3  ) bad = bad+1
  if (tsn>0.0.or.tsn<-100.0) bad = bad+1
  do k=1,NkIce ; if (t_ice(k) >0.0 .or. t_ice(k) < -100.0) bad = bad+1 ; enddo

  if (bad>0) then
    print *, 'BAD ICE AFTER RESIZE ', 'hs/hi=',hs,hi,'tsn/tice=',&
                      tsn,t_ice(:),'tmelt/bmelt=',tmelt,bmelt
  endif
end subroutine resize_check

subroutine unpack_check(hs, tsn, hi, t_ice, NkIce, cn)
  real, intent(in) :: hs, tsn, hi, cn
  real, dimension(NkIce), intent(in) :: t_ice
  integer, intent(in) :: NkIce
  integer :: k, bad

  bad = 0
  if (hs <0.0.or.hs > 1e3  ) bad = bad+1
  if (hi <0.0.or.hi > 1e3  ) bad = bad+1
  if (tsn>0.0.or.tsn<-100.0) bad = bad+1
  do k=1,NkIce ; if (t_ice(k) >0.0 .or. t_ice(k) < -100.0) bad = bad+1 ; enddo

  if (bad>0) print *, 'BAD ICE AFTER UNPACK ', 'hs/hi=',hs,hi,'tsn/tice=', &
                      tsn,t_ice(:),'cn=',cn
end subroutine unpack_check


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! e_to_melt_TS - energy needed to melt a given snow/ice configuration.         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function e_to_melt_TS(hs, tsn, hi, T, S)
  real, intent(in) :: hs, tsn, hi
  real, dimension(:), intent(in) :: T, S
  real             :: e_to_melt_TS

  integer :: k, nk_ice
  real :: I_nk_ice
  
  nk_ice = size(T)
  I_nk_ice = 1.0 / real(nk_ice)

  e_to_melt_TS = (DS*hs) * (LI - CI*tsn)
  do k=1,nk_ice
    ! The commented out lines are correct, but will cause a change in answers
    ! from the previous solution.
    if (T(k) < -MU_TS*S(k)) then
      e_to_melt_TS = e_to_melt_TS + ((DI*hi)*I_nk_ice) * &
                     (LI - CI*T(k)) * (1.0 + MU_TS*S(k)/T(k))
    else  ! This layer is already melted and has excess heat.
      e_to_melt_TS = e_to_melt_TS + ((DI*hi)*I_nk_ice) * &
                     CW * (-T(k) - MU_TS*S(k))
    endif
  enddo

end function e_to_melt_TS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_resize_SIS2 - An n-layer code for applying snow and ice thickness and    !
!    temperature changes due to thermodynamic forcing.                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_resize_SIS2(hsno, tsn, hice, tice, Sice, snow, frazil, evap, &
                          tmlt, bmlt, tfw, NkIce, heat_to_ocn, h2o_to_ocn, &
                          h2o_from_ocn, snow_to_ice, bablt )
  real, intent(inout) :: hsno        ! snow thickness (m-snow)
  real, intent(inout) :: tsn         ! snow temperature (deg-C)
  real, intent(inout) :: hice        ! ice thickness (m-ice)
  real, dimension(NkIce), &
        intent(inout) :: tice ! ice temperature by layer (deg-C)
  real, dimension(NkIce), &
        intent(in)    :: Sice ! ice salinity by layer (g/kg)
  real, intent(in   ) :: snow        ! new snow (kg/m^2-snow)
  real, intent(in   ) :: frazil      ! frazil in energy units
  real, intent(in   ) :: evap        ! ice evaporation (kg/m^2)
  real, intent(in   ) :: tmlt       ! top melting energy (J/m^2)
  real, intent(in   ) :: bmlt       ! bottom melting energy (J/m^2)
  real, intent(in   ) :: tfw         ! seawater freezing temperature (deg-C)
  integer, intent(in   ) :: NkIce ! The number of ice layers.
  real, intent(  out) :: heat_to_ocn ! energy left after ice all melted (J/m^2)
  real, intent(  out) :: h2o_to_ocn  ! liquid water flux to ocean (kg/m^2)
  real, intent(  out) :: h2o_from_ocn! evaporation flux from ocean (kg/m^2)
  real, intent(  out) :: snow_to_ice ! snow below waterline becomes ice
  real, intent(  out), optional :: bablt ! bottom ablation (kg/m^2)

  real :: tmelt, bmelt
  real, dimension(NkIce) :: hlay ! temporary ice layer thicknesses
  real, dimension(NkIce) :: lo_old, lo_new, hi_old, hi_new ! layer bounds
  real, dimension(NkIce) :: ntice
  real :: hw                  ! waterline height above ice base
  real :: evap_left, melt_left
  real :: mass, etot
  real :: overlap
  integer :: k, kold, knew

  heat_to_ocn  = 0.0
  h2o_to_ocn   = DS*hsno + DI*hice + snow-evap ! - from ice at end gives ocean h2o flux
  h2o_from_ocn = 0.0
  snow_to_ice  = 0.0

  tmelt = tmlt ; bmelt = bmlt

!  call ice_resize(hsno, tsn, hice, tice, sice, snow, frazil, evap, tmlt, bmlt, &
!             tfw, heat_to_ocn, h2o_to_ocn, h2o_from_ocn, snow_to_ice, bablt)
  do k=1,NkIce ! break out individual layers
    hlay(k) = hice/NkIce
  enddo
  ! set mass mark; will subtract mass at end for melt flux to ocean
  h2o_to_ocn = DS*hsno + DI*hice + snow - evap 
  h2o_from_ocn = 0.0 ! for excess evap-melt
  heat_to_ocn = 0.0  ! for excess melt energy

  if (hice == 0.0) hsno = 0.0
  if (hsno == 0.0) tsn = 0.0
  hsno = hsno+snow/DS ! add snow

  ! add frazil
  if (frazil > 0.0 .and. hice == 0.0) then
    do k=1,NkIce
      tice(k) = min(tfw,-MU_TS*sice(k)-Frazil_temp_offset) !was tfw  ! Why the 0.5?
      hlay(k) = hlay(k) + ((frazil/NkIce)/emelt(tice(k), sice(k)))/DI
    enddo
  endif

  if (tmelt < 0.0) then  ! this shouldn't happen
    bmelt = bmelt+tmelt
    tmelt = 0.0
  endif
  if (bmelt < 0.0 ) then ! add freezing to bottom layer
    hlay(NkIce) = hlay(NkIce) + (-bmelt/emelt(tice(NkIce),sice(NkIce)))/DI
    bmelt = 0.0
  endif

  if (evap > 0.0) then ! apply evaporation mass flux
    if (evap < hsno*DS) then
      hsno = hsno-evap/DS ! evaporate part of snow
    else ! evaporate ice
      evap_left = evap - hsno*DS
      hsno = 0.0
      do k=1,NkIce
        if (evap_left < hlay(k)*DI) then
          hlay(k) = hlay(k)-evap_left/DI
          evap_left = 0.0
          exit
        endif
        evap_left = evap_left-hlay(k)*DI
        hlay(k) = 0.0
      enddo
      h2o_from_ocn = evap_left
    endif
  endif

  if (tmelt > 0.0 ) then ! apply top melt heat flux
    if (tmelt < hsno*DS*(LI-CI*tsn)) then
      hsno = hsno-(tmelt/(LI-CI*tsn))/DS ! melt part of snow
    else ! melt ice
      melt_left = tmelt - hsno*DS*(LI-CI*tsn)
      hsno = 0.0
      do k=1,NkIce
        if (melt_left < hlay(k)*DI*emelt(tice(k),sice(k))) then
          hlay(k) = hlay(k) - (melt_left/emelt(tice(k), sice(k)))/DI ! melt part layer
          melt_left = 0.0
          exit
        endif
        melt_left = melt_left - hlay(k)*DI*emelt(tice(k), sice(k)) ! melt whole layer
        hlay(k) = 0.0
      enddo
      heat_to_ocn = heat_to_ocn + melt_left ! melt heat left after snow & ice gone
    endif
  endif

  ! apply bottom melt heat flux
  if (present(bablt)) then ! set mark to capture bottom melt
    bablt = DS*hsno
    do k=1,NkIce
      bablt = bablt+DI*hlay(k)
    enddo
  endif

  melt_left = bmelt
  if (melt_left>0.0) then ! melt ice from below
    do k=NkIce,1,-1
      if (melt_left < hlay(k)*DI*emelt(tice(k),sice(k))) then
        hlay(k) = hlay(k) - (melt_left/emelt(tice(k), sice(k)))/DI ! melt part layer
        melt_left = 0.0
        exit
      endif
      melt_left = melt_left - hlay(k)*DI*emelt(tice(k), sice(k)) ! melt whole layer
      hlay(k) = 0.0
    enddo
  
  endif
  if (melt_left > 0.0 ) then ! melt snow from below
    if (melt_left < hsno*DS*(LI-CI*tsn)) then
      hsno = hsno - (melt_left/(LI-CI*tsn))/DS ! melt part of snow
    else
      heat_to_ocn = heat_to_ocn + melt_left - DS*hsno*(LI-CI*tsn)
      hsno = 0.0              ! melt heat left after snow & ice gone
    endif
  endif

  if (present(bablt)) then ! dif mark to capture bottom melt
    bablt = bablt-DS*hsno
    do k=1,NkIce
      bablt = bablt-DI*hlay(k)
    enddo
  endif

  ! snow below waterline adjustment
  hice = 0.0
  do k=1,NkIce
    hice = hice+hlay(k) ! sum ice thickness
  enddo
  hw = (DI*hice+DS*hsno)/DW ! water height above ice base
  if (hw>hice) then ! convert snow to ice to maintain ice top at waterline
    snow_to_ice = (hw - hice)*DI ! need this much ice mass from snow
    if (snow_to_ice <= hsno*DS) then
      hsno = hsno - snow_to_ice/DS
    else
      ! snow_to_ice = (hw-hice)*DI = (((DI-DW)*hice+DS*hsno)/DW)*DI =
      !             = (DI/DW)*(DS*hsno) - (DW/DI-1)*(DI*hice)
      !             <= (DI/DW)*(DS*hsno)  (since DW > DI > 0 & hice >= 0)
      !             <= (DS*hsno)         (since DW > DI > 0 & hsno >= 0) 
      ! So this else branch can never happen.
      hsno = 0.0
    endif    

! add ice to a layer (currently, salinity of target layer does not change).
!    call add_ice(hlay, tice, sice, 1, snow_to_ice, &
!                 emelt2temp(LI-CI*tsn, sice(1)), sice(1))
    mass = hlay(1)*DI + snow_to_ice
    etot = emelt(tice(1), sice(1)) * hlay(1)*DI + &
           emelt(emelt2temp(LI-CI*tsn, sice(1)), sice(1)) * snow_to_ice
    hlay(1) = mass/DI
    tice(1) = emelt2temp(etot/mass, sice(1))
  else
    snow_to_ice = 0.0;
  endif

  ! even up ice layers
  ! call even_up(hice, hlay, tice, sice)
  hice = 0.0
  do k=1,NkIce ! set old layer bounds while summing thickness
    lo_old(k) = hice
    hice = hice+hlay(k)
    hi_old(k) = hice
    ntice(k) = 0.0
  enddo

  if (hice==0.0) then ! no ice - don't bother with rest of evening process
    tice(:) = -MU_TS*sice(:)
  else
    do k=1,NkIce ! new layer bounds - evenly spaced
      lo_new(k) = (k-1)*hice/NkIce
      hi_new(k) = lo_new(k)+hice/NkIce
    enddo

    do kold=1,NkIce
      do knew=1,NkIce
        overlap = min(hi_old(kold),hi_new(knew))-max(lo_old(kold),lo_new(knew))
        if (overlap > 0.0) then
          ntice(knew) = ntice(knew)+overlap*emelt(tice(kold),sice(kold)) ! add in h*emelt
        endif
      enddo
    enddo

    do k=1,NkIce
      tice(k) = emelt2temp(ntice(k)*NkIce/hice, sice(k))  ! retrieve temp.
      if (tice(k)>-MU_TS*sice(k)+1.0d-10) print *, 'ERROR: ICE TEMP=',tice(k), &
                                           '>' ,-MU_TS*sice(k)
    enddo
  endif

  h2o_to_ocn = h2o_to_ocn+h2o_from_ocn ! correct mark for leftover evap thru ice
  h2o_to_ocn = h2o_to_ocn-DS*hsno-DI*hice

  call resize_check(hsno, tsn, hice, tice, NkIce, bmelt, tmelt)

end subroutine ice_resize_SIS2

end module SIS2_ice_thm
