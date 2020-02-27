!> Handles the diagnostics of the ice state.
module SIS_ice_diags

! This file is part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   SIS2 is a SEA ICE MODEL for coupling through the GFDL exchange grid. SIS2  !
! is a revision of the original SIS with have extended capabilities, including !
! the option of using a B-grid or C-grid spatial discretization.  The SIS2     !
! software has been extensively reformulated from SIS for greater consistency  !
! with the Modular Ocean Model, version 6 (MOM6), and to permit might tighter  !
! dynamical coupling between the ocean and sea-ice.                            !
!   This module handles diagnostics of the sea-ice state.                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser, only : get_param, read_param, log_param, log_version, param_file_type
use MOM_time_manager, only : time_type, time_type_to_real, real_to_time
! use MOM_time_manager, only : operator(+), operator(-)
! use MOM_time_manager, only : operator(>), operator(*), operator(/), operator(/=)
use coupler_types_mod, only: coupler_type_initialized, coupler_type_send_data

use SIS_continuity,    only : SIS_continuity_CS, summed_continuity, ice_cover_transport
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl, safe_alloc_alloc
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_hor_grid,  only : SIS_hor_grid_type
use SIS_types,     only : ocean_sfc_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_types,     only : ice_state_type, IST_chksum, IST_bounds_check
use SIS_utils,     only : get_avg, post_avg
use SIS2_ice_thm,  only : get_SIS2_thermo_coefs, enthalpy_liquid_freeze, Temp_from_En_S
use ice_grid,      only : ice_grid_type

implicit none ; private

#include <SIS2_memory.h>

public :: post_ocean_sfc_diagnostics, post_ice_state_diagnostics, register_ice_state_diagnostics

!> This structure has the IDs used for sea-ice state diagnostics.
type, public :: ice_state_diags_type ; private
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                                   !! timing of diagnostic output.

  !>@{ Diagnostic IDs
  integer :: id_fax=-1, id_fay=-1

  ! These are the diagnostic ids for describing the ice state.
  integer :: id_mib=-1, id_mi=-1
  integer, dimension(:), allocatable :: id_t, id_sal
  integer :: id_cn=-1, id_hi=-1, id_hp=-1, id_hs=-1, id_tsn=-1, id_ext=-1 ! id_hp mw/new
  integer :: id_t_iceav=-1, id_s_iceav=-1, id_e2m=-1

  integer :: id_simass=-1, id_sisnmass=-1, id_sivol=-1
  integer :: id_siconc=-1, id_sithick=-1, id_sisnconc=-1, id_sisnthick=-1
  integer :: id_siu=-1, id_siv=-1, id_sispeed=-1, id_sitimefrac=-1
  !!@}
end type ice_state_diags_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Offer diagnostics of the slowly evolving sea ice state.
subroutine post_ice_state_diagnostics(IDs, IST, OSS, IOF, dt_slow, Time, G, IG, diag)
  type(ice_state_diags_type), pointer       :: IDs !< The control structure for the SIS_dyn_trans module
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in)    :: OSS !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(ice_ocean_flux_type),  intent(in)    :: IOF !< A structure containing fluxes from the ice to
                                                   !! the ocean that are calculated by the ice model.
  real,                       intent(in)    :: dt_slow  !< The time interval of these diagnostics
  type(time_type),            intent(in)    :: Time     !< The ending time of these diagnostics
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(SIS_diag_ctrl),        pointer       :: diag !< A structure that is used to regulate diagnostic output

  ! Local variables
  real, dimension(G%isc:G%iec,G%jsc:G%jec) :: mass, mass_ice, mass_snow, tmp2d
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,IG%NkIce) :: &
    temp_ice    ! A diagnostic array with the ice temperature [degC].
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
    temp_snow   ! A diagnostic array with the snow temperature [degC].
  ! ### This diagnostic does not exist yet.
  ! real, dimension(SZI_(G),SZJ_(G),IG%CatIce) :: &
  !   rdg_frac    ! fraction of ridged ice per category
  real, dimension(SZI_(G),SZJ_(G)) :: diagVar ! A temporary array for diagnostics.
  real, dimension(IG%NkIce) :: S_col ! Specified thermodynamic salinity of each
                                     ! ice layer if spec_thermo_sal is true.
  real :: rho_ice  ! The nominal density of sea ice [kg m-3].
  real :: rho_snow ! The nominal density of snow [kg m-3].
  real :: enth_units, I_enth_units
  real :: tmp_mca  ! A temporary cell averaged mass [H ~> kg m-2].
  real :: I_Nk        ! The inverse of the number of layers in the ice.
  logical :: spec_thermo_sal
  logical :: do_temp_diags
  integer :: i, j, k, l, m, isc, iec, jsc, jec, ncat, NkIce

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  NkIce = IG%NkIce
  I_Nk = 1.0 / NkIce

  ! Sum the concentration weighted mass for diagnostics.
  if ((IDs%id_mi>0) .or. (IDs%id_mib>0) .or. (IDs%id_simass>0) .or. (IDs%id_sisnmass>0)) then
    mass_ice(:,:) = 0.0
    mass_snow(:,:) = 0.0
    mass(:,:) = 0.0
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      mass_ice(i,j) = mass_ice(i,j) + IG%H_to_kg_m2*IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
      mass_snow(i,j) = mass_snow(i,j) + IG%H_to_kg_m2*IST%mH_snow(i,j,k)*IST%part_size(i,j,k)
      mass(i,j) = mass_ice(i,j) + mass_snow(i,j)
    enddo ; enddo ; enddo

    if (IDs%id_simass>0) call post_data(IDs%id_simass, mass_ice(isc:iec,jsc:jec), diag)
    if (IDs%id_sisnmass>0) call post_data(IDs%id_sisnmass, mass_snow(isc:iec,jsc:jec), diag)
    if (IDs%id_mi>0) call post_data(IDs%id_mi, mass(isc:iec,jsc:jec), diag)

    if (IDs%id_mib>0) then
      if (associated(IOF%mass_berg)) then ; do j=jsc,jec ; do i=isc,iec
        mass(i,j) = (mass(i,j) + IOF%mass_berg(i,j)) ! Add icebergs mass [kg m-2]
      enddo ; enddo ; endif
      call post_data(IDs%id_mib, mass(isc:iec,jsc:jec), diag)
    endif
  endif

  !
  ! Thermodynamic state diagnostics
  !
  if (IDs%id_cn>0) call post_data(IDs%id_cn, IST%part_size(:,:,1:ncat), diag)
  if (IDs%id_siconc>0) call post_data(IDs%id_siconc, sum(IST%part_size(:,:,1:ncat),3), diag)

  ! TK Mod: 10/18/02
  !  if (IDs%id_obs_cn>0) call post_data(IDs%id_obs_cn, Obs_cn_ice(:,:,2), diag)
  ! TK Mod: 10/18/02: (commented out...does not compile yet... add later)
  !  if (IDs%id_obs_hi>0) &
  !    call post_avg(IDs%id_obs_hi, Obs_h_ice(isc:iec,jsc:jec), IST%part_size(isc:iec,jsc:jec,1:), &
  !                  diag, G=G, wtd=.true.)

  !   Convert from ice and snow enthalpy back to temperature for diagnostic purposes.
  do_temp_diags = (IDs%id_tsn > 0)
  do m=1,NkIce ; if (IDs%id_t(m)>0) do_temp_diags = .true. ; enddo
  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, enthalpy_units=enth_units, &
                             rho_ice=rho_ice, rho_snow=rho_snow, &
                             specified_thermo_salinity=spec_thermo_sal)
  I_enth_units = 1.0 / enth_units

  if (do_temp_diags) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
      if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k) > 0.0) then
        if (spec_thermo_sal) then ; do m=1,NkIce
          temp_ice(i,j,k,m) = temp_from_En_S(IST%enth_ice(i,j,k,m), S_col(m), IST%ITV)
        enddo ; else ; do m=1,NkIce
          temp_ice(i,j,k,m) = temp_from_En_S(IST%enth_ice(i,j,k,m), &
                                              IST%sal_ice(i,j,k,m), IST%ITV)
        enddo ; endif
      else
        do m=1,NkIce ; temp_ice(i,j,k,m) = 0.0 ; enddo
      endif
      if (IST%part_size(i,j,k)*IST%mH_snow(i,j,k) > 0.0) then
        temp_snow(i,j,k) = temp_from_En_S(IST%enth_snow(i,j,k,1), 0.0, IST%ITV)
      else
        temp_snow(i,j,k) = 0.0 ! ### Should this be = temp_ice(i,j,k,1)?
      endif
    enddo ; enddo ; enddo
  endif

  if (IDs%id_ext>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (IST%part_size(i,j,0) < 0.85) diagVar(i,j) = 1.0
    enddo ; enddo
    call post_data(IDs%id_ext, diagVar, diag)
  endif
  if (IDs%id_hp>0) call post_avg(IDs%id_hp, IST%mH_pond, IST%part_size(:,:,1:), & ! mw/new
                                 diag, G=G, &
                                 scale=IG%H_to_kg_m2/1e3, wtd=.true.) ! rho_water=1e3
  if (IDs%id_hs>0) call post_avg(IDs%id_hs, IST%mH_snow, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_snow, wtd=.true.)
  if (IDs%id_sisnthick>0) call post_avg(IDs%id_sisnthick, IST%mH_snow, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_snow, wtd=.true.)
  if (IDs%id_hi>0) call post_avg(IDs%id_hi, IST%mH_ice, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_ice, wtd=.true.)
  if (IDs%id_sithick>0) call post_avg(IDs%id_sithick, IST%mH_ice, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_ice, wtd=.true.)
  if (IDs%id_sivol>0) call post_avg(IDs%id_sivol, IST%mH_ice, IST%part_size(:,:,1:), &
                                 diag, G=G, scale=IG%H_to_kg_m2/Rho_ice, wtd=.true.)
  if (IDs%id_tsn>0) call post_avg(IDs%id_tsn, temp_snow, IST%part_size(:,:,1:), &
                                 diag, G=G, wtd=.true.)
  if (IDs%id_sitimefrac>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (IST%part_size(i,j,0) < 1.0) diagVar(i,j) = 1.0
    enddo ; enddo
    call post_data(IDs%id_sitimefrac, diagVar, diag)
  endif
  if (IDs%id_sisnconc>0) then
    diagVar(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec; do k=1,ncat
      if (IST%part_size(i,j,k) > 0.0 .and. IST%mH_snow(i,j,k) > 0.0) then
        diagVar(i,j) = diagVar(i,j) + IST%part_size(i,j,k)
      endif
    enddo ; enddo ; enddo
    call post_data(IDs%id_sisnconc, diagVar, diag)
  endif

  do m=1,NkIce
    if (IDs%id_t(m)>0) call post_avg(IDs%id_t(m), temp_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   diag, G=G, wtd=.true.)
    if (IDs%id_sal(m)>0) call post_avg(IDs%id_sal(m), IST%sal_ice(:,:,:,m), IST%part_size(:,:,1:), &
                                   diag, G=G, wtd=.true.)
  enddo
  if (IDs%id_t_iceav>0) call post_avg(IDs%id_t_iceav, temp_ice, IST%part_size(:,:,1:), &
                                    diag, G=G, wtd=.true.)
  if (IDs%id_S_iceav>0) call post_avg(IDs%id_S_iceav, IST%sal_ice, IST%part_size(:,:,1:), &
                                    diag, G=G, wtd=.true.)

  ! Write out diagnostics of the ocean surface state, as seen by the slow sea ice.
  ! These fields do not change over the course of the sea-ice time stepping.
  call post_ocean_sfc_diagnostics(OSS, dt_slow, Time, G, diag)

  if (IDs%id_e2m>0) then
    tmp2d(:,:) = 0.0
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)>0.0) then
      tmp2d(i,j) = tmp2d(i,j) + IST%part_size(i,j,k)*IST%mH_snow(i,j,k)*IG%H_to_kg_m2 * &
                       ((enthalpy_liquid_freeze(0.0, IST%ITV) - &
                         IST%enth_snow(i,j,k,1)) * I_enth_units)
      if (spec_thermo_sal) then ; do m=1,NkIce
        tmp2d(i,j) = tmp2d(i,j) + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*IG%H_to_kg_m2*I_Nk) * &
                       ((enthalpy_liquid_freeze(S_col(m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
      enddo ; else ; do m=1,NkIce
        tmp2d(i,j) = tmp2d(i,j) + (IST%part_size(i,j,k)*IST%mH_ice(i,j,k)*IG%H_to_kg_m2*I_Nk) * &
                       ((enthalpy_liquid_freeze(IST%sal_ice(i,j,k,m), IST%ITV) - &
                         IST%enth_ice(i,j,k,m)) * I_enth_units)
      enddo ; endif
    endif ; enddo ; enddo ; enddo
    call post_data(IDs%id_e2m,  tmp2d(:,:), diag)
  endif

  !TOM> preparing output field fraction of ridged ice rdg_frac = (ridged ice volume) / (total ice volume)
  !     in each category; IST%rdg_mice is ridged ice mass per unit total area throughout the code.
!     if (IDs%id_rdgf>0) then
!       !$OMP parallel do default(shared) private(tmp_mca)
!       do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
!         tmp_mca = IST%mH_ice(i,j,k)*IST%part_size(i,j,k)
!         if (tmp_mca > US%kg_m3_to_R*Rho_Ice*1.e-5*US%m_to_Z) then  ! 1 mm ice thickness x 1% ice concentration
!           rdg_frac(i,j,k) = IST%rdg_mice(i,j,k) / tmp_mca
!         else
!           rdg_frac(i,j,k) = 0.0
!         endif
!       enddo ; enddo ; enddo
!       call post_data(IDs%id_rdgf, rdg_frac(isc:iec,jsc:jec), diag)
!     endif

end subroutine post_ice_state_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Offer diagnostics of the ocean surface field, as seen by the sea ice.
subroutine post_ocean_sfc_diagnostics(OSS, dt_slow, Time, G, diag)
  type(ocean_sfc_state_type), intent(in)    :: OSS  !< A structure containing the arrays that describe
                                                    !! the ocean's surface state for the ice model.
  real,                       intent(in)    :: dt_slow  !< The time interval of these diagnostics
  type(time_type),            intent(in)    :: Time     !< The ending time of these diagnostics
  type(SIS_hor_grid_type),    intent(inout) :: G    !< The horizontal grid type
  type(SIS_diag_ctrl),        pointer       :: diag !< A structure that is used to regulate diagnostic output

  real :: Idt_slow ! The inverse of the thermodynamic step [s-1].
  Idt_slow = 0.0 ; if (dt_slow > 0.0) Idt_slow = 1.0/dt_slow

  ! Write out diagnostics of the ocean surface state, as seen by the slow sea ice.
  ! These fields do not change over the course of the sea-ice time stepping.
  if (OSS%id_sst>0) call post_data(OSS%id_sst, OSS%SST_C, diag)
  if (OSS%id_sss>0) call post_data(OSS%id_sss, OSS%s_surf, diag)
  if (OSS%id_ssh>0) call post_data(OSS%id_ssh, OSS%sea_lev, diag)
  if (allocated(OSS%u_ocn_C)) then
    if (OSS%id_uo>0) call post_data(OSS%id_uo, OSS%u_ocn_C, diag)
    if (OSS%id_vo>0) call post_data(OSS%id_vo, OSS%v_ocn_C, diag)
  else
    if (OSS%id_uo>0) call post_data(OSS%id_uo, OSS%u_ocn_B, diag)
    if (OSS%id_vo>0) call post_data(OSS%id_vo, OSS%v_ocn_B, diag)
  endif
  if (OSS%id_frazil>0) &
    call post_data(OSS%id_frazil, OSS%frazil*Idt_slow, diag)

  if (coupler_type_initialized(OSS%tr_fields)) &
    call coupler_type_send_data(OSS%tr_fields, Time)

end subroutine post_ocean_sfc_diagnostics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Do the registration calls for diagnostics of the ice state.
subroutine register_ice_state_diagnostics(Time, IG, param_file, diag, IDs)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model time.
  type(ice_grid_type),         intent(in)    :: IG   !< The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(ice_state_diags_type),  pointer       :: IDs  !< A structure for regulating sea ice state diagnostics.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=8) :: nstr
  integer :: n, nLay
  real, parameter :: missing = -1e34

  nLay = IG%NkIce

  if (.not.associated(IDs)) allocate(IDs)
  call log_version(param_file, "SIS_ice_diagnostics", version, &
     "This module handles sea-ice state diagnostics.")

  IDs%diag => diag

  ! Ice state diagnostics.
  IDs%id_ext = register_diag_field('ice_model', 'EXT', diag%axesT1, Time, &
               'ice modeled', '0 or 1', missing_value=missing)
  IDs%id_cn       = register_diag_field('ice_model', 'CN', diag%axesTc, Time, &
               'ice concentration', '0-1', missing_value=missing)
  IDs%id_hp       = register_diag_field('ice_model', 'HP', diag%axesT1, Time, &
               'pond thickness', 'm-pond', missing_value=missing) ! mw/new
  IDs%id_hs       = register_diag_field('ice_model', 'HS', diag%axesT1, Time, &
               'snow thickness', 'm-snow', missing_value=missing)
  IDs%id_tsn      = register_diag_field('ice_model', 'TSN', diag%axesT1, Time, &
               'snow layer temperature', 'C',  missing_value=missing)
  IDs%id_hi       = register_diag_field('ice_model', 'HI', diag%axesT1, Time, &
               'ice thickness', 'm-ice', missing_value=missing)
  IDs%id_sitimefrac = register_diag_field('ice_model', 'sitimefrac', diag%axesT1, Time, &
               'time fraction of ice cover', '0-1', missing_value=missing)
  IDs%id_siconc = register_diag_field('ice_model', 'siconc', diag%axesT1, Time, &
               'ice concentration', '0-1', missing_value=missing)
  IDs%id_sithick  = register_diag_field('ice_model', 'sithick', diag%axesT1, Time, &
               'ice thickness', 'm-ice', missing_value=missing)
  IDs%id_sivol  = register_diag_field('ice_model', 'sivol', diag%axesT1, Time, &
               'ice volume', 'm-ice', missing_value=missing)
  IDs%id_sisnconc = register_diag_field('ice_model', 'sisnconc', diag%axesT1, Time, &
               'snow concentration', '0-1', missing_value=missing)
  IDs%id_sisnthick= register_diag_field('ice_model', 'sisnthick', diag%axesT1, Time, &
               'snow thickness', 'm-snow', missing_value=missing)

  IDs%id_t_iceav = register_diag_field('ice_model', 'T_bulkice', diag%axesT1, Time, &
               'Volume-averaged ice temperature', 'C', missing_value=missing)
  IDs%id_s_iceav = register_diag_field('ice_model', 'S_bulkice', diag%axesT1, Time, &
               'Volume-averaged ice salinity', 'g/kg', missing_value=missing)
  call safe_alloc_ids_1d(IDs%id_t, nLay)
  call safe_alloc_ids_1d(IDs%id_sal, nLay)
  do n=1,nLay
    write(nstr, '(I4)') n ; nstr = adjustl(nstr)
    IDs%id_t(n)   = register_diag_field('ice_model', 'T'//trim(nstr), &
                 diag%axesT1, Time, 'ice layer '//trim(nstr)//' temperature', &
                 'C',  missing_value=missing)
    IDs%id_sal(n)   = register_diag_field('ice_model', 'Sal'//trim(nstr), &
               diag%axesT1, Time, 'ice layer '//trim(nstr)//' salinity', &
               'g/kg',  missing_value=missing)
  enddo

  IDs%id_mi   = register_diag_field('ice_model', 'MI', diag%axesT1, Time, &
               'ice + snow mass', 'kg/m^2', missing_value=missing)
  IDs%id_simass = register_diag_field('ice_model', 'simass', diag%axesT1, Time, &
               'ice mass', 'kg/m^2', missing_value=missing)
  IDs%id_sisnmass = register_diag_field('ice_model', 'sisnmass', diag%axesT1, Time, &
               'snow mass', 'kg/m^2', missing_value=missing)
  IDs%id_mib  = register_diag_field('ice_model', 'MIB', diag%axesT1, Time, &
               'ice + snow + bergs mass', 'kg/m^2', missing_value=missing)
  IDs%id_e2m  = register_diag_field('ice_model','E2MELT' ,diag%axesT1, Time, &
               'heat needed to melt ice', 'J/m^2', missing_value=missing)

!### THIS DIAGNOSTIC IS MISSING.
!  IDs%id_ta    = register_diag_field('ice_model', 'TA', diag%axesT1, Time, &
!            'surface air temperature', 'C', missing_value=missing)
!### THESE DIAGNOSTICS DO NOT EXIST YET.
!  IDs%id_rdgf    = register_diag_field('ice_model','RDG_FRAC' ,diag%axesT1, Time, &
!               'ridged ice fraction', '0-1', missing_value=missing)
end subroutine register_ice_state_diagnostics


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Allocate an array of integer diagnostic arrays and set them to -1, if they are not already allocated
subroutine safe_alloc_ids_1d(ids, nids)
  integer, allocatable, intent(inout) :: ids(:) !< An array of diagnostic IDs to allocate
  integer,              intent(in)    :: nids   !< The number of IDs to allocate

  if (.not.ALLOCATED(ids)) then
    allocate(ids(nids)) ; ids(:) = -1
  endif;
end subroutine safe_alloc_ids_1d


end module SIS_ice_diags
