!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!                       FIVE-LAYER VERTICAL THERMODYNAMICS                     !
!                                                                              !
! Reference:                                                                   !
!   Winton, M., 2011:  A conservative non-iterative n-layer sea ice            !
!   temperature solver, in prep.                                               !
!                                                                              !
!                                                                              !
!         ->+---------+ <- ts - diagnostic surface temperature ( <= 0C )       !
!        /  |         |                                                        !
!      hs   |  snow   | <- tsn   One snow layer with heat capapcity
!        \  |         |                                                        !
!         =>+---------+                                                        !
!        /  |         |                                                        !
!       /   |         | <- t1    Four salty ice layers with heat capacity      !
!      /    |         |                                                        !
!     /     |         | <- t2                                                  !
!   hi      |...ice...|                                                        !
!     \     |         | <- t3                                                  !
!      \    |         |                                                        !
!       \   |         | <- t4                                                  !
!        \  |         |                                                        !
!         ->+---------+ <- base of ice fixed at seawater freezing temp.        !
!                                                                              !
!                                         Mike Winton (Michael.Winton@noaa.gov)!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

! TK mod:  SLAB_ICE treatment modified to follow supersource
!          (after Bryan 1969).  The conductive heat
!           flux from ice to atmosphere is computed based on
!           an effective ice thickness which ensures a minimum
!           thickness of 1.7cm for the calculation.

module ice_thm_mod


use constants_mod, only : LI => hlf ! latent heat of fusion - 334e3 J/(kg-ice)

implicit none ; private

public :: DS, DI, DW, MU_TS, TFI, CI, slab_ice_optics, get_thermo_coefs
public :: ice5lay_temp, ice5lay_resize, ice_thm_param, e_to_melt
          ! test driver needs line below
          !,LI, KS, KI, CI, DT, SI1, SI2, SI3, SI4

interface e_to_melt
  module procedure e_to_melt_4, e_to_melt_TS
end interface

!
! properties of ice, snow, and seawater (NCAR CSM values)
!

!real, parameter :: LI = 334e3 ! temporary replacement for constants_mod ref

real            :: KS    = 0.31      ! conductivity of snow - 0.31 W/(mK)
real, parameter :: DS    = 330.0     ! density of snow - 330 kg/(m^3)
real, parameter :: KI    = 2.03      ! conductivity of ice  - 2.03 W/(mK)
real, parameter :: DI    = 905.0     ! density of ice  - 905 kg/(m^3)
real, parameter :: CI    = 2.1e3      ! heat cap. of fresh ice - 2100 J/(kg K)
! salinities from S(z) = 0.5*3.2*(1.0-cos(3.1416*z**(0.407/(z+0.573))))
! z=[1 3 5 7]/8 ; ref: Hunke et al: CICE V. 4.0, 2008, p. 26
real, parameter :: SI1   = 0.65      ! salinity of sea ice top layer
real, parameter :: SI2   = 2.35      ! salinity of sea ice second layer
real, parameter :: SI3   = 3.03      ! salinity of sea ice third layer
real, parameter :: SI4   = 3.19      ! salinity of sea ice bottom layer

real, parameter :: MU_TS = 0.054     ! relates freezing temp. to salinity
real, parameter :: TFI   = -MU_TS*SI1! top ice freezing temp. = -mu*salinity
real, parameter :: CW    = 4.2e3     ! heat capacity of seawater 4200 J/(kg K)
real, parameter :: DW    = 1030.0    ! density of water for waterline - kg/(m^3)

! albedos are from CSIM4 assumming 0.53 visible and 0.47 near-ir insolation

real            :: H_LO_LIM = 0.0       ! hi/hs lower limit for temp. calc.
real            :: H_SUBROUNDOFF = 1e-35 ! A miniscule value compared with H_LO_LIM
real            :: FRAZIL_TEMP_OFFSET = 0.5 ! A temperature offset between the
logical         :: SLAB_ICE = .false.   ! should we do old style GFDL slab ice?
!
! slab ice specific parameters
!
real, parameter :: CRIT_THICKNESS       = 1.00
real, parameter :: T_RANGE              = 10.0
real, parameter :: MIN_ICE_ALB          = 0.55   ! coupled model uses 0.55
real, parameter :: MAX_ICE_ALB          = 0.80
real, parameter :: ALB_OCEAN            = 0.10

integer, parameter :: NN = 4
real               :: DT = 1800.0

!
! In the ice temperature calculation we place a limit to below (salinity
! dependent) freezing point on the prognosed temperatures.  For ice_resize
! it is better to make a slightly more restrictive limit that requires the
! temperature to be such that the brine content is less than "liq_lim" of
! the total mass.  That is T_f/T < liq_lim implying T<T_f/liq_lim
real, parameter :: liq_lim = .99

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
function ecolumn(hsno, hice, tice, sice) result (retval)
real retval
real :: hsno, hice
real, dimension(NN) :: tice, sice

  integer :: k

  retval = 0.0
  if ( hsno>0.0) then
    retval = -DS*hsno*LI
  endif

  do k=1,NN
    retval = retval-DI*(hice/NN)*emelt(tice(k), sice(k))
  enddo

end function ecolumn

function ice_temp_est(A, B, sol, hsno, tsn, hice, tice, sice, tfw) &
result (temp_est)
real, dimension(0:NN) :: temp_est ! snow (0) and ice (1:NN) estimated temps
real, intent(in) :: A        ! net down heat flux from above at 0C (W/m^2)
real, intent(in) :: B        ! d(up sfc heat flux)/d(ts) [W/(m^2 deg-C)]
real, intent(in), dimension(0:NN) :: sol  ! solar absorbed by ice layers (W/m^2)
real, intent(in) :: hsno                  ! snow thickness (m)
real, intent(in) :: tsn                   ! snow temperature
real, intent(in) :: hice                  ! ice thickness (m)
real, intent(in), dimension(NN) :: tice   ! ice temperature (deg-C)
real, intent(in), dimension(NN) :: sice   ! ice salinity (ppt)
real, intent(in) :: tfw                   ! seawater freezing temperature (deg-C)

  real :: kk, k10, k0a, tsf, k0a_x_ta, salt_part, rat, tsurf
  real, dimension(0:NN) :: aa, bb, bbb, cc, ff ! tridiagonal coefficients
  real, dimension(0:NN) :: bb_new, ff_new ! modified by tridiag. algorithm
  integer :: k

  kk = NN*KI/hice                        ! full ice layer conductivity
!  k10 = 1/((hice/(2*NN*KI))+hsno/(2*KS)) ! coupling ice layer 1 to snow
!  k0a = 1/(hsno/(2*KS)+1/B)            ! coupling snow to "air"
!  ta  = A/B                            ! "air" temperature

  k10 = 2.0*(KS*(NN*KI)) / (hice*KS + hsno*(NN*KI)) ! coupling ice layer 1 to snow
  k0a = (KS*B) / (0.5*B*hsno + KS)       ! coupling snow to "air"
  k0a_x_ta = (KS*A) / (0.5*B*hsno + KS)  ! coupling times "air" temperture
  tsf = -MU_TS*sice(1)                   ! surface freezing temperature
  if (hsno>0.0) tsf = 0.0

  ! initialize tridiagonal matrix coefficients
  do k=1,NN   ! load bb with heat capacity term (also used in ff)
    salt_part = 0.0
    if (sice(k)>0.0) salt_part = -MU_TS*sice(k)*LI/(tice(k)*tice(k))
    bb(k) = (hice/NN)*(DI/DT)*(CI-salt_part) ! add coupling to this later
  enddo

  aa(NN) = 0.0              ! bottom of ice is at top of matrix, aa=0
  cc(NN) = -kk
  ff(NN) = sol(NN)+bb(NN)*tice(NN)+2*kk*tfw
  bbb(NN) = bb(NN) + 3*kk    ! add coupling

  do k=2,NN-1
    aa(k) = -kk; cc(k) = -kk;
    ff(k) = sol(k)+bb(k)*tice(k)
    bbb(k) = bb(k) + 2*kk  ! = bb(k) - aa(k) - cc(k)
  enddo

  aa(1) = -kk
  cc(1) = -k10
  ff(1) = sol(1)+bb(1)*tice(1)
  bbb(1) = bb(1) + kk + k10  ! = bb(1) - aa(1) - cc(1)

  ! if melting, change coefficients for fixed surf. temp. @ melting
  aa(0) = -k10
  bb(0) = hsno*(DS/DT)*CI+k0a
  bbb(0) = hsno*(DS/DT)*CI+k10+k0a    ! if melting, change this
  cc(0) = 0.0              ! snow is at bottom of matrix, cc=0
  ff(0) = hsno*(DS/DT)*CI*tsn+k0a_x_ta ! if melting, change this

  ! Notes:  aa(k) = cc(k+1) <= 0 ; bbb(k) = bb(k) - aa(k) - cc(k) ; bb(k) >= 0.0

  !
  ! going UP the ice column (down the tridiagonal matrix)
  !
  bb_new(NN) = bbb(NN); ff_new(NN) = ff(NN)
  do k=NN-1,0,-1
    rat = aa(k)/bb_new(k+1) ! -1 <= rat <= 0
    bb_new(k) = bbb(k) - rat*cc(k+1)
    ff_new(k) = ff(k) - rat*ff_new(k+1)
  end do

  temp_est(0) = ff_new(0)/bb_new(0)
  tsurf = (A*hsno+2*KS*temp_est(0))/(B*hsno+2*KS) ! diagnose surface skin temp.

  if (tsurf > tsf) then ! surface is melting, redo with surf. at melt temp.
    tsurf = tsf
    aa(0) = -k10*hsno                      ! mult. thru by hsno for hsno==0 case
    bbb(0) = hsno*(hsno*(DS/DT)*CI + k10) + 2*KS ! swap coupling to atmos out and
    ff(0) = hsno*hsno*(DS/DT)*CI*tsn +2*KS*tsf  ! coupling to melting surface in
    rat = aa(0)/bb_new(0+1)
    bb_new(0) = bbb(0) - rat*cc(0+1)   ! reset, rat is left over from end of loop
    ff_new(0) = ff(0) - rat*ff_new(0+1) ! rat=aa(1)/bb(2), so it's unchanged
    temp_est(0) = ff_new(0)/bb_new(0)
  endif

  !
  ! going DOWN the ice column (up the tridiagonal matrix)
  !
  do k=1,NN
    temp_est(k) = (ff_new(k) - cc(k)*temp_est(k-1)) / bb_new(k)
  end do

  return
end function ice_temp_est

!
! laytemp - implicit calculation of new layer temperature
!
function laytemp(m, tfi, f, b, tp) result (retval)
real retval
real :: m    ! mass of ice - kg/m2
real :: tfi  ! ice freezing temp. (determined by salinity)
real :: f    ! forcing - W/m2
real :: b    ! response of outward heat flux to local temperature - W/m2/K
real :: tp   ! prior step temperature

  real :: mdt, AA, BB, CC
  ! solve quadratic equation for layer temperature, t:
  !
  !   (m/DT) * {CI-LI*tfi/(t*tp)} * (t-tp) = f - b*t
  !
  !Niki: There could be a problem if tp==0 but tfi/=0 in which case the above equation is ill-defined
  !      but this code returns -f/b
  !
  !if( tp == 0.0 .and. tfi /= 0.0 ) print*, 'BAD laytemp INPUT ',tp,tfi
  !

  mdt  = m/DT                   ! mass divided by timestep - kg/m2/s

  if ( tfi == 0.0 ) then
    retval = (mdt*CI*tp+f)/(mdt*CI+b) ! = -BB/AA
  else
    AA = mdt*CI+b
    BB = -(mdt*LI*tfi/tp+mdt*CI*tp+f)
    CC = mdt*LI*tfi
    ! This form avoids round-off errors.
    !  if (BB >= 0) then
    retval = -(BB + sqrt(BB*BB-4*AA*CC))/(2*AA)
    !  else
    !    retval = (2*CC) / (-BB + sqrt(BB*BB - 4*AA*CC))
    !  endif
  endif
  return
end function laytemp

subroutine ice_temp(A, B, sol, tsurf, hsno, tsn, hice, tice, sice, tfw, fb, &
                                                                tmelt, bmelt)
real, intent(in   ) :: A        ! net down heat flux from above at 0C (W/m^2)
real, intent(in   ) :: B        ! d(up sfc heat flux)/d(ts) [W/(m^2 deg-C)]
! sol - solar absorbed by snow/ice layers (W/m^2)
real, intent(in   ), dimension(0:NN) :: sol
real, intent(inout) :: tsurf    ! surface temperature (deg-C)
real, intent(in   ) :: hsno     ! snow thickness (m)
real, intent(inout) :: tsn      ! snow temperature (deg-C)
real, intent(in   ) :: hice     ! ice thickness (m)
real, intent(inout), dimension(NN) :: tice ! ice temperature (deg-C)
real, intent(in   ), dimension(NN) :: sice ! ice salinity (ppt)
real, intent(in   ) :: tfw      ! seawater freezing temperature (deg-C)
real, intent(in   ) :: fb       ! heat flux from ocean to ice bottom (W/m^2)
real, intent(inout) :: tmelt    ! accumulated top melting energy  (J/m^2)
real, intent(inout) :: bmelt    ! accumulated bottom melting energy (J/m^2)

  real, dimension(0:NN) :: temp_est
  real, dimension(NN) :: tfi, tice_est ! estimated new ice temperatures
  real :: mi, ms, e_extra
  real :: kk, k10, k0a, tsf, k0a_x_ta, tsno_est,hie
  integer :: k

  mi = DI*hice/NN           ! full ice layer mass
  ms = DS*hsno              ! full snow layer mass
  tfi(:) = -MU_TS*sice(:)   ! freezing temperature of ice layers
  hie = max(hice, H_LO_LIM); ! prevent thin ice inaccuracy (mw)
  kk = NN*KI/hie                        ! full ice layer conductivity
!  k10 = 1/(hie/(2*NN*KI)+hsno/(2*KS))   ! coupling ice layer 1 to snow
!  k0a = 1/(hsno/(2*KS)+1/B)            ! coupling snow to "air"
!  ta  = A/B                            ! "air" temperature

  k10 = 2.0*(KS*(NN*KI)) / (hice*KS + hsno*(NN*KI)) ! coupling ice layer 1 to snow
  k0a = (KS*B) / (0.5*B*hsno + KS)      ! coupling snow to "air"
  k0a_x_ta = (KS*A) / (0.5*B*hsno + KS) ! coupling times "air" temperture
  tsf = tfi(1)                          ! surface freezing temperature
  if (hsno>0.0) tsf = 0.0

  ! 1st get non-conservative estimate with implicit treatment of layer coupling
 temp_est=ice_temp_est(A, B, sol, hsno, tsn, hie, tice, sice, tfw)
 tice_est=temp_est(1:NN)
 tsno_est=temp_est(0)
!
! following two lines disable implicit coupling estimate;  in a 10 year CORE
! test the skin temperature without implicit estimate was mostly within .05K
! but was cooler around topography in NH (up to 1K - Baffin Island).  Thickness
! looked very similar
!
! tice_est=tice
! tsno_est=tsn

  !
  ! conservative pass going UP the ice column
  !
  tice_est(NN) = laytemp(mi, tfi(NN), sol(NN)+kk*(2*tfw+tice_est(NN-1)), &
                                                             3*kk, tice(NN))
  do k=NN-1,2,-1
    tice_est(k) = laytemp(mi, tfi(k), &
                       sol(k)+kk*(tice_est(k-1)+tice_est(k+1)), 2*kk, tice(k))
  enddo
  tice_est(1) = laytemp(mi, tfi(1), sol(1)+kk*tice_est(2)+k10*tsno_est, &
                                                              kk+k10, tice(1))
  tsno_est = laytemp(ms, 0.0, sol(0)+k10*tice_est(1)+k0a_x_ta, k10+k0a, tsn)
  tsurf = (A*hsno+2*KS*tsno_est)/(B*hsno+2*KS)  ! diagnose surface skin temp.

  if (tsurf > tsf) then ! surface is melting, redo with surf. at melt temp.
    tsurf = tsf
    tsno_est = laytemp(hsno*ms, 0.0, hsno*sol(0)+hsno*k10*tice_est(1)+2*KS*tsf,&
                       hsno*k10+2*KS, tsn) ! note: mult. thru by hsno
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
  tice(1)  = laytemp(mi, tfi(1), sol(1)+kk*tice_est(1+1) &
                                       +k10*(tsno_est-tice_est(1)), kk, tice(1))
  do k=2,NN-1 ! flux from above is fixed, only have downward feedback
    tice(k) = laytemp(mi, tfi(k), sol(k)+kk*tice_est(k+1) &
                                       +kk*(tice(k-1)-tice_est(k)), kk, tice(k))
  enddo
  tice(NN) = laytemp(mi, tfi(NN), sol(NN)+2*kk*tfw &
                              +kk*(tice(NN-1)-tice_est(NN)), 2*kk, tice(NN))
  !
  ! END conservative update
  !

  bmelt = bmelt+DT*(fb-2*kk*(tfw-tice(NN))) ! add in bottom melting/freezing

  if (tsn > 0.0) then ! put excess snow energy into top melt
    e_extra = CI*DS*hsno*tsn
    tmelt = tmelt+e_extra
    tsn = 0.0
  endif

  do k=1,NN
    if (tice(k)>tfi(k)/liq_lim) then ! push excess energy to closer of top or bottom melt
      e_extra = (emelt(tfi(k)/liq_lim,sice(k))-emelt(tice(k),sice(k)))*DI*hie/NN
      tice(k) = tfi(k)/liq_lim
      if (k<=NN/2) then
        tmelt = tmelt+e_extra
      else
        bmelt = bmelt+e_extra
      endif
    endif
  enddo
  return
end subroutine ice_temp

!
! even ice layers after growth/melt; returns new ice thickness
!
subroutine even_up(hice, hlay,tice,sice)
real, intent(inout) :: hice
real, intent(in   ), dimension(NN) :: hlay
real, intent(inout), dimension(NN) :: tice
real, intent(in   ), dimension(NN) :: sice

  real, dimension(NN) :: lo_old, lo_new, hi_old, hi_new ! layer bounds
  real, dimension(NN) :: ntice
  integer :: k, kold, knew
  real :: overlap

  hice = 0.0
  do k=1,NN ! set old layer bounds while summing thickness
    lo_old(k) = hice
    hice = hice+hlay(k)
    hi_old(k) = hice
    ntice(k) = 0.0
  enddo

  if (hice==0.0) then ! no ice - don't bother with rest of evening process
    tice = -MU_TS*sice
    return
  endif

  do k=1,NN ! new layer bounds - evenly spaced
    lo_new(k) = (k-1)*hice/NN
    hi_new(k) = lo_new(k)+hice/NN
  enddo

  do kold=1,NN
    do knew=1,NN
      overlap = min(hi_old(kold),hi_new(knew))-max(lo_old(kold),lo_new(knew))
      if (overlap > 0.0) then
        ntice(knew) = ntice(knew)+overlap*emelt(tice(kold),sice(kold)) ! add in h*emelt
      endif
    enddo
  enddo

  do k=1,NN
    tice(k) = emelt2temp(ntice(k)*NN/hice, sice(k))  ! retrieve temp.
    if (tice(k)>-MU_TS*sice(k)+1.0d-10) print *, 'ERROR: ICE TEMP=',tice(k), &
                                         '>' ,-MU_TS*sice(k)
  enddo

  return
end subroutine even_up

!
! add ice to a layer (currently, salinity of target layer does not change)
!
subroutine add_ice(hlay, tice, sice, lay, mi, ti, si)
real, intent(inout), dimension(NN) :: hlay ! ice layer thicknesses
real, intent(inout), dimension(NN) :: tice ! ice layer temperatures
real, intent(in   ), dimension(NN) :: sice ! ice layer salinities
integer, intent(in   ) :: lay      ! target layer
real, intent(in   ) ::    mi       ! mass of added ice
real, intent(in   ) ::    ti       ! temperature of added ice
real, intent(in   ) ::    si       ! salinity of added ice

  real mass, etot

  mass = hlay(lay)*DI+mi
  etot = emelt(tice(lay), sice(lay))*hlay(lay)*DI + emelt(ti, si)*mi
  hlay(lay) = mass/DI
  tice(lay) = emelt2temp(etot/mass, sice(lay))
  return
end subroutine add_ice

!
!  ice_resize - ice & snow thickness change
!
subroutine ice_resize(hsno, tsn, hice, tice, sice, snow, frazil, evap, tmelt, &
                      bmelt, tfw, heat_to_ocn, h2o_to_ocn, h2o_from_ocn, &
                      snow_to_ice, bablt)
real, intent(inout) :: hsno         ! snow thickness (m-snow)
real, intent(inout) :: tsn          ! snow temperature (deg-C)
real, intent(inout) :: hice         ! ice thickness (m-ice)
real, intent(inout), dimension(NN) ::  tice ! ice layer temperature (deg-C)
real, intent(in   ), dimension(NN) ::  sice ! ice layer salinity (g/kg)
real, intent(in   )  :: snow        ! new snow (kg/m^2-snow)
real, intent(in   )  :: frazil      ! frazil in energy units (J/m^2 to be extracted)
real, intent(in   )  :: evap        ! ice evaporation (kg/m^2)
real, intent(inout)  :: tmelt       ! top melting energy (J/m^2)
real, intent(inout)  :: bmelt       ! bottom melting energy (J/m^2)
real, intent(in   )  :: tfw         ! frazil temperature - seawater freezing (deg-C)
real, intent(inout) :: heat_to_ocn  ! energy left after ice all melted (J/m^2)
real, intent(inout) :: h2o_to_ocn   ! liquid water flux to ocean (kg/m^2)
real, intent(inout) :: h2o_from_ocn ! evaporation flux from ocean (kg/m^2)
real, intent(inout) :: snow_to_ice  ! snow below waterline becomes ice
real, intent(inout), optional :: bablt ! bottom ablation (kg/m^2)

  real, dimension(NN) :: hlay ! temporary ice layer thicknesses
  real :: hw                  ! waterline height above ice base
  real :: evap_left, melt_left
  real, dimension(NN) :: t_frazil  ! The temperature which with the frazil-ice is created, in C.
  real, dimension(NN) :: h_frazil  ! The newly-formed thickness of frazil ice, in m.
  integer :: k

  do k=1,NN ! break out individual layers
    hlay(k) = hice/NN
  enddo
  ! set mass mark; will subtract mass at end for melt flux to ocean
  h2o_to_ocn = DS*hsno+DI*hice+snow-evap
  h2o_from_ocn = 0.0 ! for excess evap-melt
  heat_to_ocn = 0.0  ! for excess melt energy

  if (hice == 0.0) hsno = 0.0
  if (hsno == 0.0) tsn = 0.0
  hsno = hsno+snow/DS ! add snow

  ! add frazil
  if (frazil > 0.0) then
    do k=1,NN
      t_frazil(k) = min(tfw,-MU_TS*sice(k)-Frazil_temp_offset) !was tfw  ! Why the 0.5?
      h_frazil(k) = ((frazil/NN)/emelt(t_frazil(k), sice(k)))/DI
    enddo
    if (hice == 0.0) then
      do k=1,NN
        tice(k) = T_frazil(k)
        hlay(k) = hlay(k) + h_frazil(k)
      enddo
    else
      do k=1,NN
        tice(k) = (hlay(k)*tice(k) + h_frazil(k)*T_frazil(k)) / (hlay(k) + h_frazil(k))
        hlay(k) = hlay(k) + h_frazil(k)
      enddo
    endif
  endif

  if (tmelt < 0.0) then  ! this shouldn't happen
    bmelt = bmelt+tmelt
    tmelt = 0.0
  endif
  if (bmelt < 0.0 ) then ! add freezing to bottom layer
    hlay(NN) = hlay(NN) + (-bmelt/emelt(tice(NN),sice(NN)))/DI
    bmelt = 0.0
  endif

  if (evap > 0.0) then ! apply evaporation mass flux
    if (evap < hsno*DS) then
      hsno = hsno-evap/DS ! evaporate part of snow
    else ! evaporate ice
      evap_left = evap - hsno*DS
      hsno = 0.0
      do k=1,NN
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
      do k=1,NN
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
    do k=1,NN
      bablt = bablt+DI*hlay(k)
    enddo
  endif

  melt_left = bmelt
  if (melt_left>0.0) then ! melt ice from below
    do k=NN,1,-1
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
    do k=1,NN
      bablt = bablt-DI*hlay(k)
    enddo
  endif

  ! snow below waterline adjustment
  hice = 0.0
  do k=1,NN
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

       !Niki: How to balance this?
!       print*, 'UNBALANCED hsno conversion ',hsno,snow_to_ice/DS
      hsno = 0.0
    endif
    call add_ice(hlay, tice, sice, 1, snow_to_ice, &
                 emelt2temp(LI-CI*tsn, sice(1)), sice(1))
  else
    snow_to_ice = 0.0;
  endif

  ! even up ice layers
  call even_up(hice, hlay, tice, sice)

  h2o_to_ocn = h2o_to_ocn+h2o_from_ocn ! correct mark for leftover evap thru ice
  h2o_to_ocn = h2o_to_ocn-DS*hsno-DI*hice

  return
end subroutine ice_resize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_thm_param - set ice thermodynamic parameters                             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_thm_param(slab_ice_in, ks_in, h_lo_lim_in )
  logical, intent(in) :: slab_ice_in
  real, intent(in)    :: ks_in
  real, intent(in)    :: h_lo_lim_in


  SLAB_ICE    = slab_ice_in

  KS          = ks_in
  H_LO_LIM    = h_lo_lim_in
  H_SUBROUNDOFF = max(1e-35, 1e-20*H_LO_LIM)

end subroutine ice_thm_param

subroutine slab_ice_optics(hs, hi, ts, tfw, albedo)
  real, intent(in   ) :: hs  ! snow thickness (m-snow)
  real, intent(in   ) :: hi  ! ice thickness (m-ice)
  real, intent(in   ) :: ts  ! surface temperature
  real, intent(in   ) :: tfw ! seawater freezing temperature
  real, intent(  out) :: albedo ! ice surface albedo (0-1)
  real :: alb, as, ai, cs
  real :: thick_ice_alb, tcrit, fh

  tcrit = tfw - T_RANGE
  if (ts <= tcrit) then
    thick_ice_alb = MAX_ICE_ALB
  else if (ts >= tfw) then
    thick_ice_alb = MIN_ICE_ALB
  else
    thick_ice_alb = MAX_ICE_ALB + (MIN_ICE_ALB-MAX_ICE_ALB)*(ts-tcrit)/T_RANGE
  endif

  if (hi >= crit_thickness) then
    albedo = thick_ice_alb
  else
    albedo = ALB_OCEAN + (thick_ice_alb-ALB_OCEAN)*sqrt(hi/CRIT_THICKNESS)
  endif

end subroutine slab_ice_optics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice5lay_temp - ice & snow temp. change [Winton (2000) section 2.a]           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice5lay_temp(hs, tsn, hi, t1, t2, t3, t4, ts, A, B, &
                        IS, I1, I2, I3, I4, tfw, fb, dtt, tmelt, bmelt)
real, intent(in   ) :: hs    ! snow thickness (m)
real, intent(inout) :: tsn   ! snow temperature (deg-C)
real, intent(in   ) :: hi    ! ice thickness (m)
real, intent(inout) :: t1    ! top ice temperature (deg-C)
real, intent(inout) :: t2    ! second layer ice temperature (deg-C)
real, intent(inout) :: t3    ! third layer ice temperature (deg-C)
real, intent(inout) :: t4    ! bottom ice temperature (deg-C)
real, intent(  out) :: ts    ! surface temperature (deg-C)
real, intent(in   ) :: A     ! net surface heat flux (+ up) at ts=0 (W/m^2)
real, intent(in   ) :: B     ! d(sfc heat flux)/d(ts) [W/(m^2 deg-C)]
real, intent(in   ) :: IS    ! solar absorbed by snow layer (W/m^2)
real, intent(in   ) :: I1    ! solar absorbed by ice layer 1 (W/m^2)
real, intent(in   ) :: I2    ! solar absorbed by ice layer 1 (W/m^2)
real, intent(in   ) :: I3    ! solar absorbed by ice layer 1 (W/m^2)
real, intent(in   ) :: I4    ! solar absorbed by ice layer 1 (W/m^2)
real, intent(in   ) :: tfw   ! seawater freezing temperature (deg-C)
real, intent(in   ) :: fb    ! heat flux from ocean to ice bottom (W/m^2)
real, intent(in   ) :: dtt   ! timestep (sec)
real, intent(inout) :: tmelt ! accumulated top melting energy  (J/m^2)
real, intent(inout) :: bmelt ! accumulated bottom melting energy (J/m^2)
!
! variables for temperature calculation [see Winton (1999) section II.A.]
! note:  here equations are multiplied by hi to improve thin ice accuracy
!
    real, dimension(NN) :: tice, sice
    real, dimension(0:NN) :: sol ! layer 0 for snow, 1:NN for ice

! TK Mods:
real :: hi_effective
real :: KI_over_eps = 1.7065e-2     ! 5/2.93 from Bryan (1969);
!                                 Value used in SS tsc.F (1.7065 cm)
!                                  converted to meters...

  DT = dtt ! set timestep from argument - awkward, remove later

  if (SLAB_ICE) then
    hi_effective = hi + KI_over_eps     ! TK added
    ts = (KI*tfw-A*hi_effective) / (KI+B*hi_effective)     ! TK mod
    if (ts > 0.0) then       ! surface melting conditions
       ts = 0.0
       if (hi>0.0) tmelt = tmelt + (KI*tfw/hi_effective - A)*dt ! TK mod
    endif
    if (hi>0.0) then
       bmelt = bmelt + (fb-KI*(tfw-ts)/hi_effective)*dt     ! TK mod
    else
       bmelt = bmelt + (fb-A-B*tfw)*dt
    endif
    return
  endif

  sol(0) = IS
  sol(1) = I1 ;  sol(2) = I2 ;  sol(3) = I3 ;  sol(4) = I4
  tice(1) = t1; tice(2) = t2; tice(3) = t3; tice(4) = t4
  sice(1) = SI1; sice(2) = SI2; sice(3) = SI3; sice(4) = SI4

  call ice_temp(-A, B, sol, ts, hs, tsn, hi, tice, sice, tfw, fb, tmelt, bmelt)
  t1 = tice(1); t2 = tice(2); t3 = tice(3); t4 = tice(4);
  call temp_check(ts, hs, tsn, hi, tice(:), 4, bmelt, tmelt)

end subroutine ice5lay_temp


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
! e_to_melt - energy needed to melt a given 4-layer snow/ice configuration.    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function e_to_melt_4(hs, tsn, hi, t1, t2, t3, t4)
    real, intent(in) :: hs, tsn, hi, t1, t2, t3, t4
    real             :: e_to_melt_4

  e_to_melt_4 = DS*hs*(LI-CI*tsn) &
                      +DI*hi*(CI-LI/t1)*(-MU_TS*SI1-t1)/4 &
                      +DI*hi*(CI-LI/t2)*(-MU_TS*SI2-t2)/4 &
                      +DI*hi*(CI-LI/t3)*(-MU_TS*SI3-t3)/4 &
                      +DI*hi*(CI-LI/t4)*(-MU_TS*SI4-t4)/4
end function e_to_melt_4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! e_to_melt - energy needed to melt a given snow/ice configuration             !
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
! ice5lay_resize - ice & snow thickness change [Winton (1998) section II.B.]   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice5lay_resize(hs, tsn, hi, t1, t2, t3, t4, snow, frazil, evap, &
                          tmelt, bmelt, tfw, heat_to_ocn, h2o_to_ocn,      &
                          h2o_from_ocn, snow_to_ice, bablt                 )
real, intent(inout) :: hs          ! snow thickness (m-snow)
real, intent(inout) :: tsn         ! snow temperature (deg-C)
real, intent(inout) :: hi          ! ice thickness (m-ice)
real, intent(inout) :: t1          ! temperature of top ice (deg-C)
real, intent(inout) :: t2          ! temperature of second ice (deg-C)
real, intent(inout) :: t3          ! temperature of third ice (deg-C)
real, intent(inout) :: t4          ! temperature of bottom ice (deg-C)
real, intent(in   ) :: snow        ! new snow (kg/m^2-snow)
real, intent(in   ) :: frazil      ! frazil in energy units
real, intent(in   ) :: evap        ! ice evaporation (kg/m^2)
real, intent(in   ) :: tmelt       ! top melting energy (J/m^2)
real, intent(in   ) :: bmelt       ! bottom melting energy (J/m^2)
real, intent(in   ) :: tfw         ! seawater freezing temperature (deg-C)
real, intent(  out) :: heat_to_ocn ! energy left after ice all melted (J/m^2)
real, intent(  out) :: h2o_to_ocn  ! liquid water flux to ocean (kg/m^2)
real, intent(  out) :: h2o_from_ocn! evaporation flux from ocean (kg/m^2)
real, intent(  out) :: snow_to_ice ! snow below waterline becomes ice
real, intent(  out), optional :: bablt ! bottom ablation (kg/m^2)

real :: h1, h2, dh, hw
real, dimension(NN) :: sice, tice
real :: tmlt, bmlt

  heat_to_ocn  = 0.0
  h2o_to_ocn   = DS*hs+DI*hi+snow-evap ! - from ice at end gives ocean h2o flux
  h2o_from_ocn = 0.0
  snow_to_ice  = 0.0

  if (SLAB_ICE) then
    !
    ! add snow and frazil
    !
    hi = hi + snow/DI + frazil/(DI*LI)
    t1 = tfw; t2 = tfw;
    !
    ! atmospheric evaporation
    !
    if (evap <= hi*DI) then
      hi = hi - evap/DI
    else
      h2o_from_ocn = evap-hi*DI
      hi = 0.0
    end if
    !
    ! ... melting
    !
    hi = hi - (tmelt+bmelt)/(DI*LI)
    if (hi<0.0) then
       heat_to_ocn = -hi*DI*LI
       hi = 0.0
    else
       heat_to_ocn = 0.0
    endif

    h2o_to_ocn = h2o_to_ocn+h2o_from_ocn ! reset mark for leftover evap thru ice
    h2o_to_ocn = h2o_to_ocn-DI*hi-DS*hs  ! hs should be zero
    if (present(bablt)) bablt = bmelt/LI
    return
  endif


  tice(1) = t1; tice(2) = t2; tice(3) = t3; tice(4) = t4;
  sice(1) = SI1; sice(2) = SI2; sice(3) = SI3; sice(4) = SI4;
  tmlt = tmelt; bmlt = bmelt;
  call ice_resize(hs, tsn, hi, tice, sice, snow, frazil, evap, tmlt, bmlt, &
             tfw, heat_to_ocn, h2o_to_ocn, h2o_from_ocn, snow_to_ice, bablt)
  t1 = tice(1); t2 = tice(2); t3 = tice(3); t4 = tice(4);
  call resize_check(hs, tsn, hi, tice, 4, bmelt, tmelt)
  return
end subroutine ice5lay_resize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! get_thermo_coefs - return various thermodynamic coefficients.                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine get_thermo_coefs(ice_salinity)
  real, dimension(:), optional, intent(out) :: ice_salinity
! Arguments: ice_salinity - The specified salinity of each layer when the
!                           thermodynamic salinities are pre-specified.
  integer k, nk

  if (present(ice_salinity)) then
    nk = size(ice_salinity)
    if (nk >= 1) ice_salinity(1) = SI1
    if (nk >= 2) ice_salinity(2) = SI2
    if (nk >= 3) ice_salinity(3) = SI3
    if (nk >= 4) ice_salinity(4) = SI4
    do k=5,nk ;  ice_salinity(k) = SI4 ; enddo
  endif

end subroutine get_thermo_coefs

end module ice_thm_mod
