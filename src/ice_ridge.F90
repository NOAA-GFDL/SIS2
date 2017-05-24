module ice_ridging_mod
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of SIS2.                                        *
!*                                                                     *
!* SIS2 is free software; you can redistribute it and/or modify it and *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* SIS2 is distributed in the hope that it will be useful, but WITHOUT *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser, only : get_param, log_param, read_param, log_version, param_file_type
use MOM_domains,     only : pass_var, pass_vector, BGRID_NE
use SIS_hor_grid, only : SIS_hor_grid_type
use fms_io_mod,      only : register_restart_field, restart_file_type

implicit none ; private

#include <SIS2_memory.h>

public :: ice_ridging, ridge_rate, ice_ridging_init

real, parameter :: hlim_unlim = 1.e8   ! arbitrary huge number used in ice_ridging
real    :: s2o_frac       = 0.5        ! fraction of snow dumped into ocean during ridging
logical :: rdg_lipscomb = .true.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_ridging_init - initialize the ice ridging and set parameters.            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_ridging_init - some preparations for the ridging routine,                !
!                    called within subroutine ridging                          !
!                    T. Martin, 2008                                           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_ridging_init(km,cn, hi, part_undef, part_undef_sum, &
                            hmin, hmax, efold, rdg_ratio)
  integer, intent(in)  :: km
  real, dimension(0:  ), intent(in)  :: cn             ! including open water fraction
  real, dimension(1:  ), intent(in)  :: hi             ! CAUTION: hi represents ice volume here
  real, dimension(1:km), intent(out) :: hmin, hmax     ! minimum and maximum ice thickness involved in Hibler's ridged ice distribution
  real, dimension(1:km), intent(out) :: efold          ! e-folding scale lambda of Lipscomb's ridged ice distribution
  real, dimension(1:km), intent(out) :: rdg_ratio      ! ratio of ; k_n in CICE
  real, dimension(0:km), intent(out) :: part_undef     ! fraction of undeformed ice or open water participating in ridging
  real,                  intent(out) :: part_undef_sum

  real, parameter       :: raft_max = 1.0 ! maximum thickness [m] of ice that rafts
  real, dimension(-1:km) :: ccn            ! cumulative ITD (or G in CICE documentation)
  real, dimension(1:km) :: thick
! now set in namelist (see ice_dyn_param):
!Niki: The following two were commented out
  real                  :: part_par       ! participation function shape parameter in parent ice categories, G* or a* in CICE
  real                  :: dist_par       ! distribution function shape parameter in receiving categories, H* or mu in CICE
  real                  :: part_pari
  real, dimension(0:km) :: rdgtmp
  integer :: k
  logical :: rdg_lipscomb = .true.

  ! actual ice thickness
  do k=2,km
    ! hi/cn is actual ice thickness; hi is volume here
    thick(k) = 0.0 ; if (cn(k)>1.e-10) thick(k) = hi(k)/cn(k)
  enddo

  ! cumulative ice thickness distribution
  ccn(0) = 0.0     ! helps to include calculations for open water in k-loops
  do k=1,km
    if (cn(k) > 1.e-10) then
      ccn(k) = ccn(k-1)+cn(k)
    else
      ccn(k) = ccn(k-1)
    endif
  enddo
  ! normalize ccn
  do k=1,km ; ccn(k) = ccn(k) / ccn(km) ; enddo

  !----------------------------------------------------------------------------------------
  ! i)  participation function for undeformed, thin ice:
  !     amount of ice per thin ice category participating in ridging
  ! ii) parameters defining the range of the thick ice categories
  !     to which the newly ridged ice is redistributed
  !
  ! A) scheme following Lipscomb et al., 2007, JGR
  !   i)  negative exponential participation function for level ice
  !   ii) negative exponential distribution of newly ridged ice
  !
  ! B) scheme following Thorndike et al., 1975, JGR
  !   i)  linear participation function with negative slope and root at (part_par,0.0)
  !    and Hibler, 1980, Mon.Wea.Rev.
  !   ii) uniform distribution of newly ridged ice in receiving categories
  !----------------------------------------------------------------------------------------
  if (rdg_lipscomb) then
    ! ************
    ! *   Ai     *
    ! ************
    part_par = -0.05    ! this is -a* for practical reasons,
                        ! part_par(lipscomb)=part_par(thorn-hib)/3 for best comparability of schemes
    !do k=1,km
    !part_undef(k) = (exp(ccn(k-1)/part_par)-exp(ccn(k)/part_par)) / (1.-exp(1./part_par))
    !enddo
    ! set in namelist: part_par  = -20.   ! this is -1/a*; CICE standard is a* = 0.05
    part_pari = 1./(1.-exp(part_par))
    do k=0,km ; rdgtmp(k) = exp(ccn(k)*part_par) * part_pari ; enddo
    do k=1,km ; part_undef(k) = rdgtmp(k-1) - rdgtmp(k) ; enddo
    ! ************
    ! *   Aii    *
    ! ************
    dist_par = 4.0
    ! set in namelist: dist_par = 4.0   ! unit [m**0.5], e-folding scale of ridged ice,
                    ! for comparable results of Lipscomb and Thorn-Hib schemes choose
                    ! 3 & 25, 4 & 50, 5 & 75 or 6 & 100
    do k=2,km
      if (thick(k)>0.0) then
        efold(k) = dist_par * sqrt(thick(k))
        hmin(k)  = min(2.*thick(k), thick(k)+raft_max)
        rdg_ratio(k) = (hmin(k)+efold(k)) / thick(k)
      else
        efold(k)=0.0; hmin(k)=0.0; rdg_ratio(k)=1.0
      endif
    enddo
  !----------------------------------------------------------------------------------------
  else
    ! ************
    ! *   Bi     *
    ! ************
    part_par = 0.15
    ! set in namelist: part_par = 0.15   ! CICE standard is 0.15
    do k=1,km
      if (ccn(k) < part_par) then
        part_undef(k) = (ccn(k)  -ccn(k-1)) * (2.-( (ccn(k-1)+ccn(k))  /part_par )) / part_par
      else if (ccn(k) >= part_par .and. ccn(k-1) < part_par) then
        part_undef(k) = (part_par-ccn(k-1)) * (2.-( (ccn(k-1)+part_par)/part_par )) / part_par
      else
        part_undef(k) = 0.0
      endif
    enddo
       ! ************
       ! *   Bii    *
       ! ************
  ! set in namelist: dist_par = 50.0   ! 25m suggested by CICE and is appropriate for multi-year ridges [Flato and Hibler, 1995, JGR],
                               ! 50m gives better fit to first-year ridges [Amundrud et al., 2004, JGR]
    do k=2,km
      if (thick(k)>0.0) then
        hmin(k) = min(2.*thick(k), thick(k)+raft_max)             ! minimum and maximum defining range in ice thickness space
        hmax(k) = max(2.*sqrt(dist_par*thick(k)), hmin(k)+1.e-10) !  that receives ice in the ridging process
        rdg_ratio(k) = 0.5*(hmin(k)+hmax(k))/thick(k)             ! ratio of mean ridge thickness to thickness of parent level ice
      else
        hmin(k)=0.0; hmax(k)=0.0; rdg_ratio(k)=1.0
      endif
    enddo

  endif
  !----------------------------------------------------------------------------------------

  ! ratio of net ice area removed / total area participating
  part_undef_sum = ccn(1)
  do k=2,km
    if (rdg_ratio(k)>0.0) part_undef_sum = part_undef_sum + part_undef(k)*(1.-1./rdg_ratio(k))
  enddo

end subroutine ice_ridging_init


!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ridge_rate - deformation rate or                                             !
!              total energy dissipation rate due to ridging                    !
!              (Flato and Hibler, 1995, JGR) or                                !
!              net area loss in riding (CICE documentation)                    !
!              depending on the state of the ice drift                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
function ridge_rate(del2, div) result (rnet)
  real, intent(in)  :: del2, div
  real              :: del, rnet, rconv, rshear
  !TOM> cs is now set in namelist:
!Niki: this was commented out
  real, parameter   :: cs=0.25 !(CICE documentation)

  del=sqrt(del2)

  rconv  = -min(div,0.0)           ! energy dissipated by convergence ...
  rshear = 0.5*(del-abs(div))      ! ... and by shear
  rnet   = rconv + cs*rshear       ! net energy contains only part of the
                                   !  shear energy as only a fraction is
           !  dissipated in the ridging process
  return
end function ridge_rate

!TOM>~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_ridging - mechanical redistribution of thin (undeformed) ice into        !
!               thicker (deformed/ridged) ice categories                       !
!               T. Martin, 2008                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_ridging(km, cn, hi, hs, t1, t2, age, snow_to_ocn, rdg_rate, hi_rdg, &
                       dt, hlim_in, rdg_open, vlev)

  integer, intent(in) :: km
  real, dimension(0:), intent(inout) :: cn                    ! including open water fraction
  real, dimension(1:), intent(inout) :: hi, t1, t2, hs, age   ! CAUTION: these quantities are extensive here,
                    !          i.e. hi represents ice volume
  real,                intent(out)   :: snow_to_ocn           ! total snow volume dumped into ocean during ridging
  real,                intent(in)    :: rdg_rate              ! ridging rate from subroutine ridge_rate
  real, dimension(1:), intent(inout) :: hi_rdg                ! A diagnostic of the ridged ice volume in each category.
  real,                intent(in)    :: dt                    ! time step dt has units seconds
  real, dimension(1:), intent(in)    :: hlim_in               ! ice thickness category limits
  real,                intent(out)   :: rdg_open              ! change in open water area due to newly formed ridges
  real,                intent(out)   :: vlev                  ! volume of level ice participating in ridging

  integer :: k, kd, kr, n_iterate
  integer, parameter :: n_itermax = 10 ! maximum number of iterations for redistribution
!    real, parameter :: frac_hs_rdg = 0.5 ! fraction of snow that remains on ridged ice;
  real             :: frac_hs_rdg      ! fraction of snow that remains on ridged ice;
         !  (1.0-frac_hs_rdg) falls into ocean
  real, dimension(0:km) :: part_undef  ! fraction of undeformed ice or open water participating in ridging
  real                  :: area_undef  ! fractional area of parent and ...
  real, dimension(1:km) :: area_def    ! ... newly ridged ice, respectively
  real, dimension(1:km) :: vol_def     ! fractional volume of newly ridged ice
  real, dimension(1:km) :: cn_old, hi_old   ! state of quantities at beginning of iteration loop
  real, dimension(1:km) :: rdg_frac    ! ratio of ridged and total ice volume
  real                  :: alev        ! area of level ice participating in ridging
  real                  :: ardg, vrdg  ! area and volume of newly formed rdiged (vlev=vrdg!!!)
  real, dimension(1:km) :: hmin, hmax, efold, rdg_ratio, hlim
  real                  :: hl, hr
  real                  :: cn_tot, part_undef_sum
  real                  :: div_adv, Rnet, Rdiv, Rtot, rdg_area, rdgtmp, hlimtmp
  real                  :: area_frac
  real, dimension(1:km) :: area_rdg
  real, dimension(1:km) :: frac_hi, frac_hs, frac_t1, frac_t2, frac_age
  logical               :: rdg_iterate
  !-------------------------------------------------------------------
  ! some preparations
  !-------------------------------------------------------------------
  hlimtmp = hlim_in(km)
  hlim(km) = hlim_unlim   ! ensures all ridged ice is smaller than thickest ice allowed
  frac_hs_rdg = 1.0-s2o_frac
  !snow_to_ocn = 0.0 -> done in subroutine transport
  alev=0.0; ardg=0.0; vlev=0.0; vrdg=0.0
  !
  call ice_ridging_init(km, cn, hi, part_undef, part_undef_sum, &
                        hmin, hmax, efold, rdg_ratio)

  !-------------------------------------------------------------------
  ! opening and closing rates of the ice cover
  !-------------------------------------------------------------------

  ! update total area fraction as this may exceed 1 after transportation/advection
  ! (cn_tot <= 1 after ridging!)
  cn_tot = sum(cn(0:km))

  ! dissipated energy in ridging from state of ice drift
  !  after Flato and Hibler (1995, JGR)
  !  (see subroutine ridge_rate in ice_dyn_mod),
  !  equals net closing rate times ice strength
  ! by passing to new, local variable rdg_rate is saved for diagnostic output
  Rnet = rdg_rate
  ! the divergence rate given by the advection scheme ...
  div_adv = (1.-cn_tot) / dt
  ! ... may exceed the rate derived from the drift state (Rnet)
  if (div_adv < 0.) Rnet = max(Rnet, -div_adv)
  ! non-negative opening rate that ensures cn_tot <=1 after ridging
  Rdiv = Rnet + div_adv
  ! total closing rate
  Rtot = Rnet / part_undef_sum

  !-------------------------------------------------------------------
  ! iteration of ridging redistribution
  do n_iterate=1, n_itermax
  !-------------------------------------------------------------------

    ! save initial state of ice concentration, total and ridged ice volume
    !  at beginning of each iteration loop
    do k=1,km
      cn_old(k) = cn(k)
      hi_old(k) = hi(k)

      rdg_frac(k) = 0.0 ; if (hi(k)>0.0) rdg_frac(k) = hi_rdg(k)/hi(k)
    enddo

    ! reduce rates in case more than 100% of any category would be removed
    do k=1,km ; if (cn(k)>1.e-10 .and. part_undef(k)>0.0) then
      rdg_area = part_undef(k) * Rtot * dt   ! area ridged in category k
      if (rdg_area > cn(k)) then
        rdgtmp = cn(k)/rdg_area
        Rtot = Rtot * rdgtmp
        Rdiv = Rdiv * rdgtmp
      endif
    endif ; enddo

    !-------------------------------------------------------------------
    ! redistribution of ice
    !-------------------------------------------------------------------

    ! changes in open water area
    cn(0) = max(cn(0) + (Rdiv - part_undef(1)*Rtot) * dt, 0.0)

    if (Rtot>0.0) then

      ! area, volume and energy changes in each category
      do kd=1,km   ! donating category
        area_undef = min(part_undef(kd)*Rtot*dt, cn_old(kd))   ! area that experiences ridging in category k,
                                                               ! make sure that not more than 100% are used
        if (cn_old(kd) > 1.e-10) then
          area_frac    = area_undef / cn_old(kd)              ! fraction of level ice area involved in ridging
          area_rdg(kd) = area_undef / rdg_ratio(kd)           ! area of new ridges in category k
        else
          area_frac    = 0.0
          area_rdg(kd) = 0.0
        endif
        !if (rdg_ratio(kd) > 0.0) then     ! distinguish between level and ridged ice in
        !else                              !  each category: let only level ice ridge;
        !endif                             !  sea also change of hi_rdg below

        ! reduce area, volume and energy of snow and ice in source category
        frac_hi(kd)  = hi(kd)   * area_frac
        frac_hs(kd)  = hs(kd)   * area_frac
        frac_t1(kd)  = t1(kd)   * area_frac
        frac_t2(kd)  = t2(kd)   * area_frac
        frac_age(kd) = age(kd)  * area_frac

        cn(kd)  = cn(kd)  - area_undef
        hi(kd)  = hi(kd)  - frac_hi(kd)
        hs(kd)  = hs(kd)  - frac_hs(kd)
        t1(kd)  = t1(kd)  - frac_t1(kd)
        t2(kd)  = t2(kd)  - frac_t2(kd)
        age(kd) = age(kd) - frac_age(kd)

        alev = alev + area_undef   ! diagnosing area of level ice participating in ridging
        vlev = vlev + frac_hi(kd)  ! diagnosing total ice volume moved due to ridging
                                   !  (here donating categories)

        !    Here it is assumed that level and ridged ice
        !  of a category participate in ridging in equal
        !  measure; this also means that ridged ice may be ridged again
        hi_rdg(kd) = hi_rdg(kd) - rdg_frac(kd)*frac_hi(kd)
        hi_rdg(kd) = max(hi_rdg(kd),0.0)      ! ensure hi_rdg >= 0

        ! dump part of the snow in ocean (here, sum volume, transformed to flux in update_ice_model_slow)
        snow_to_ocn = snow_to_ocn + frac_hs(kd)*(1.0-frac_hs_rdg)

      enddo

      ! split loop in order to derive frac_... variables with initial status (before ridging redistribution)
      do kd=1,km

        !----------------------------------------------------------------------------------------
        ! add area, volume and energy in receiving category :
        ! A) after Lipscomb, 2007 (negative exponential distribution)
        ! B) after Hibler, 1980, Mon. Weather Rev. (uniform distribution)
        !----------------------------------------------------------------------------------------
        if (rdg_lipscomb) then
          ! ************
          ! *   A      *
          ! ************
          if (efold(kd)>0.0) then
            do kr=1,km-1   ! receiving categories
              if (hmin(kd) >= hlim(kr)) then
                area_def(kr) = 0.0
                vol_def(kr)  = 0.0
              else
                hl = max(hmin(kd), hlim(kr-1))
                hr = hlim(kr)
                area_def(kr) = exp((hmin(kd)-hl)/efold(kd))   &
                 -             exp((hmin(kd)-hr)/efold(kd))
                vol_def(kr)  = ( (hl+efold(kd))*exp((hmin(kd)-hl)/efold(kd))   &
                 -   (hr+efold(kd))*exp((hmin(kd)-hr)/efold(kd)) ) &
                 / (hmin(kd)+efold(kd))
              endif
            enddo   ! k receiving
            ! thickest categery is a special case:
            hl = max(hmin(kd), hlim(km-1))
            area_def(km) =                  exp((hmin(kd)-hl)/efold(kd))
            vol_def(km)  = ( (hl+efold(kd))*exp((hmin(kd)-hl)/efold(kd)) ) &
             / (hmin(kd)+efold(kd))
          else
            do kr=1,km
              area_def(kr) = 0.0
              vol_def(kr)  = 0.0
            enddo
          endif
        !----------------------------------------------------------------------------------------
        else ! not rdg_lipscomb
          ! ************
          ! *   B      *
          ! ************
          if (hmax(kd)==hmin(kd)) then
            do kr=1,km ; area_def(kr) = 0.0 ; vol_def(kr)  = 0.0 ; enddo
          else
            do kr=1,km   ! receiving categories
              if (hmin(kd) >= hlim(kr) .or. hmax(kd) <= hlim(kr-1)) then
                hl = 0.0
                hr = 0.0
              else
                hl = max(hmin(kd), hlim(kr-1))
                hr = min(hmax(kd), hlim(kr)  )
              endif
              area_def(kr) = (hr   -hl   ) / (hmax(kd)   -hmin(kd)   )
              !vol_def(kr) = (hr**2-hl**2) / (hmax(kd)**2-hmin(kd)**2)
              vol_def(kr)  = area_def(kr) * (hr+hl) / (hmax(kd)+hmin(kd))
            enddo   ! k receiving
          endif

        endif
        !----------------------------------------------------------------------------------------

        ! update ice/snow area, volume, energy for receiving categories
        do kr=1,km   ! receiving categories
          cn(kr)  = cn(kr)  + area_def(kr) * area_rdg(kd)
          rdgtmp  = vol_def(kr)  * frac_hi(kd)
          hi(kr)  = hi(kr)  + rdgtmp
          hs(kr)  = hs(kr)  + vol_def(kr)  * frac_hs(kd) * frac_hs_rdg
          t1(kr)  = t1(kr)  + vol_def(kr)  * frac_t1(kd)
          t2(kr)  = t2(kr)  + vol_def(kr)  * frac_t2(kd)
          age(kr) = age(kr) + vol_def(kr)  * frac_age(kd)

          ardg = ardg + area_def(kr) * area_rdg(kd) ! diagnosing area of newly ridged ice
          vrdg = vrdg + rdgtmp                      ! diagnosing total ice volume moved due to ridging
                                                    !  (here receiving categories, cross check with vlev)

          ! add newly ridged ice volume to total ridged ice in each category
          hi_rdg(kr) = hi_rdg(kr) + rdgtmp
        enddo

      enddo   ! kd loop over donating categories

    endif ! Rtot>0.0

    ! update total area fraction and check if this is now <= 1
    ! and rerun the ice redistribution when necessary
    cn_tot = sum(cn(0:km))
    rdg_iterate = .false.
    if (abs(cn_tot-1.) > 1.e-10) then
       rdg_iterate = .true.
       div_adv = (1.-cn_tot) / dt
       Rnet    = max(0.0, -div_adv)
       Rdiv    = max(0.0,  div_adv)
       call ice_ridging_init(km, cn, hi, part_undef, part_undef_sum, &
                             hmin, hmax, efold, rdg_ratio)
       Rtot    = Rnet / part_undef_sum
    endif

    !-------------------------------------------------------------------
    if (.not. rdg_iterate) exit
  enddo   ! n_iterate
  !-------------------------------------------------------------------

  ! check ridged ice volume for natural limits
  do k=1,km
    hi_rdg(k) = max(hi_rdg(k),0.0)   ! ridged ice volume positive
    hi_rdg(k) = min(hi_rdg(k),hi(k)) ! ridged ice volume not to exceed total ice volume
  enddo

  ! calculate opening rate of ridging
  rdg_open = (alev - ardg) / dt

  ! cross check ice volume transferred from level to ridged ice
  if (abs(vlev - vrdg) > 1e-10) then
    print *,'WARNING: subroutine ice_ridging: parent ice volume does not equal ridge volume', vlev, vrdg
  endif
  ! turn vlev into a rate for diagnostics
  vlev = vlev / dt

  ! return to true upper most ice thickness category limit
  !hlim(km) = hlimtmp

end subroutine ice_ridging


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_ridging_end - deallocate the memory associated with this module.             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_ridging_end()

end subroutine ice_ridging_end

end module ice_ridging_mod
