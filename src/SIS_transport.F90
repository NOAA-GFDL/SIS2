!> Does the transport and redistribution between thickness categories for the SIS2 sea ice model.
module SIS_transport

! This file is a part of SIS2.  See LICENSE.md for the license.

use MOM_coms,          only : reproducing_sum, EFP_type, EFP_to_real, EFP_real_diff
use MOM_domains,       only : pass_var, pass_vector, BGRID_NE, CGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg, is_root_pe
use MOM_file_parser,   only : get_param, log_param, read_param, log_version, param_file_type
use MOM_hor_index,     only : hor_index_type
use MOM_obsolete_params, only : obsolete_logical, obsolete_real
use MOM_unit_scaling,  only : unit_scale_type
use SIS_continuity,    only : SIS_continuity_init, SIS_continuity_end
use SIS_continuity,    only : continuity=>ice_continuity, SIS_continuity_CS
use SIS_continuity,    only : summed_continuity, proportionate_continuity
use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use SIS_framework,     only : safe_alloc
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_tracer_advect, only : advect_tracers_thicker, SIS_tracer_advect_CS
use SIS_tracer_advect, only : advect_SIS_tracers, SIS_tracer_advect_init, SIS_tracer_advect_end
use SIS_tracer_advect, only : advect_scalar
use SIS_tracer_registry, only : SIS_tracer_registry_type, get_SIS_tracer_pointer
use SIS_tracer_registry, only : update_SIS_tracer_halos, set_massless_SIS_tracers
use SIS_tracer_registry, only : check_SIS_tracer_bounds
use SIS_types,         only : ice_state_type
use ice_grid,          only : ice_grid_type
use ice_ridging_mod,   only : ice_ridging

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_transport_init, SIS_transport_end, adjust_ice_categories
public :: alloc_cell_average_state_type, dealloc_cell_average_state_type
public :: cell_ave_state_to_ice_state, ice_state_to_cell_ave_state, cell_mass_from_CAS
public :: ice_cat_transport, finish_ice_transport

!> The SIS_transport_CS contains parameters for doing advective and parameterized advection.
type, public :: SIS_transport_CS ; private

  real :: Rho_ice             !< The nominal density of sea ice [R ~> kg m-3], used here only in
                              !! rolling and setting ridging parameters
  real :: roll_factor         !< A factor by which the propensity of small amounts of thick sea-ice
                              !! to become thinner by rolling is increased, or 0 to disable rolling.
                              !! Sensible values are 0 or larger than 1.
  real :: ice_cover_discard   !< A tiny fractional ice coverage which if positive causes the mass
                              !! in categories with less than this coverage to be discarded.

  logical :: readjust_categories !< If true, readjust the distribution into
                              !! ice thickness categories after advection.
  logical :: check_conservation !< If true, write out verbose diagnostics of conservation.
  logical :: bounds_check     !< If true, check for sensible values of thicknesses,
                              !! temperatures, salinities, tracers, etc.
  logical :: inconsistent_cover_bug !< If true, omit a recalculation of the fractional ice-free
                              !! areal coverage after the adjustment of the ice categories.
                              !! The default should be changed to false and then this option obsoleted.
  type(time_type), pointer :: Time !< A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                              !! timing of diagnostic output.
  logical :: do_ridging       !< If true, the ridging scheme is enabled.
  type(SIS_continuity_CS),    pointer :: continuity_CSp => NULL()
          !< The control structure for the SIS continuity module
  type(SIS_tracer_advect_CS), pointer :: SIS_tr_adv_CSp => NULL()
          !< The control structure for the SIS tracer advection module
  type(SIS_tracer_advect_CS), pointer :: SIS_thick_adv_CSp => NULL()
          !< The control structure for the SIS thickness advection module

  !>@{ Diagnostic IDs
  integer :: id_ix_trans = -1, id_iy_trans = -1, id_xprt = -1, id_rdgr = -1
  ! integer :: id_rdgo=-1, id_rdgv=-1 ! These do not exist yet
  !!@}

end type SIS_transport_CS

!> This structure contains a variation of the ice model state where the variables are in
!! mass per unit ocean cell area (not per unit ice area).  These are useful for conservative
!! advection, but not so useful for diagnosing ice thickness.
type, public :: cell_average_state_type ; private
  real, allocatable, dimension(:,:,:) :: m_ice  !< The mass of ice in each thickness category
                                                !! per unit total area in a cell [R Z ~> kg m-2].
  real, allocatable, dimension(:,:,:) :: m_snow !< The mass of ice in each thickness category
                                                !! per unit total area in a cell [R Z ~> kg m-2].
  real, allocatable, dimension(:,:,:) :: m_pond !< The mass of melt pond water in each thickness
                                                !! category per unit total area in a cell [R Z ~> kg m-2].
  real, allocatable, dimension(:,:,:) :: mH_ice !< The mass of ice in each thickness category
                                                !! per unit of ice area in a cell [R Z ~> kg m-2].  The
                                                !! ratio of m_ice / mH_ice gives the fractional
                                                !! ice coverage of each category.  Massless cells
                                                !! still are given plausible values of mH_ice.

  ! The following fields are used for diagnostics.
  real :: dt_sum = 0.0 !< The accumulated time since the fields were populated from an ice state type [T ~> s].
  real, allocatable, dimension(:,:) :: mass0    !< The total mass of ice, snow and melt pond water
                                                !! when the fields were populated [R Z ~> kg m-2].
  real, allocatable, dimension(:,:) :: uh_sum   !< The accumulated zonal mass fluxes of ice, snow
                                                !! and melt pond water, summed across categories,
                                                !! since the fields were populated [R Z L2 ~> kg].
  real, allocatable, dimension(:,:) :: vh_sum   !< The accumulated meridional mass fluxes of ice, snow
                                                !! and melt pond water, summed across categories,
                                                !! since the fields were populated [R Z L2 ~> kg].
  type(EFP_type) :: tot_ice                     !< The globally integrated mass of sea ice [kg].
  type(EFP_type) :: tot_snow                    !< The globally integrated mass of snow [kg].
  type(EFP_type) :: enth_ice                    !< The globally integrated sea ice enthalpy [J].
  type(EFP_type) :: enth_snow                   !< The globally integrated snow enthalpy [J].
end type cell_average_state_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_cat_transport does ice transport of mass and tracers by thickness category
subroutine ice_cat_transport(CAS, TrReg, dt_slow, nsteps, G, US, IG, CS, uc, vc, mca_tot, uh_tot, vh_tot)
  type(cell_average_state_type),     intent(inout) :: CAS !< A structure with ocean-cell averaged masses.
  type(SIS_hor_grid_type),           intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),               intent(inout) :: IG  !< The sea-ice specific grid type
  type(SIS_tracer_registry_type),    pointer       :: TrReg !< The registry of SIS ice and snow tracers.
  real,                              intent(in)    :: dt_slow !< The amount of time over which the
                                                          !! ice dynamics are to be advanced [T ~> s].
  integer,                           intent(in)    :: nsteps  !< The number of advective iterations
                                                          !! to use within this time step.
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_transport_CS),            pointer       :: CS  !< A pointer to the control structure for this module
  real, dimension(SZIB_(G),SZJ_(G)), optional, intent(in)    :: uc  !< The zonal ice velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), optional, intent(in)    :: vc  !< The meridional ice velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),0:max(nsteps,1)), optional, intent(in) :: &
    mca_tot    !< The total mass per unit total area of snow and ice summed across thickness
               !! categories in a cell, after each substep [R Z ~> kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),max(nsteps,1)), optional, intent(in) :: &
    uh_tot     !< Total zonal fluxes during each substep [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G),max(nsteps,1)), optional, intent(in) :: &
    vh_tot     !< Total meridional fluxes during each substep [R Z L2 T-1 ~> kg s-1].

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),SZCAT_(IG)) :: &
    uh_ice, &  ! Zonal fluxes of ice [R Z L2 T-1 ~> kg s-1].
    uh_snow, & ! Zonal fluxes of snow [R Z L2 T-1 ~> kg s-1].
    uh_pond    ! Zonal fluxes of melt pond water [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZCAT_(IG)) :: &
    vh_ice, &  ! Meridional fluxes of ice [R Z L2 T-1 ~> kg s-1].
    vh_snow, & ! Meridional fluxes of snow [R Z L2 T-1 ~> kg s-1].
    vh_pond    ! Meridional fluxes of melt pond water [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)) :: &
    mca0_ice, &  ! The initial mass of ice per unit ocean area in a cell [R Z ~> kg m-2].
    mca0_snow, & ! The initial mass of snow per unit ocean area in a cell [R Z ~> kg m-2].
    mca0_pond    ! The initial mass of melt pond water per unit ocean area
                 ! in a cell [R Z ~> kg m-2].
  real :: dt_adv ! An advective timestep [T ~> s]
  logical :: merged_cont
  character(len=200) :: mesg
  integer :: i, j, k, n, isc, iec, jsc, jec, isd, ied, jsd, jed, nCat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nCat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (CAS%dt_sum <= 0.0) then
    call set_massless_SIS_tracers(CAS%m_snow, TrReg, G, IG, compute_domain=.true., do_ice=.false.)
    call set_massless_SIS_tracers(CAS%m_ice, TrReg, G, IG, compute_domain=.true., do_snow=.false.)

    if (CS%bounds_check) call check_SIS_tracer_bounds(TrReg, G, IG, "SIS_transport set massless 1")
  endif

  merged_cont = (present(mca_tot) .and. present(uh_tot) .and. present(vh_tot))
  if (merged_cont .and. (present(uc) .or. present(vc))) call SIS_error(WARNING, &
    "Velocities should not be provided to ice_cat_transport when mass fluxes are provided.")
  if ((.not. merged_cont) .and. .not.(present(uc) .and. present(vc))) call SIS_error(FATAL, &
    "Either velocities or masses and mass fluxes must appear in a call to ice_cat_transport.")

  ! Do the transport via the continuity equations and tracer conservation equations
  ! for CAS%mH_ice and tracers, inverting for the fractional size of each partition.
  if (nsteps > 0) dt_adv = dt_slow / real(nsteps)
  do n = 1, nsteps
    call update_SIS_tracer_halos(TrReg, G, complete=.false.)
    call pass_var(CAS%m_ice,  G%Domain, complete=.false.)
    call pass_var(CAS%m_snow, G%Domain, complete=.false.)
    call pass_var(CAS%m_pond, G%Domain, complete=.false.)
    call pass_var(CAS%mH_ice, G%Domain, complete=.true.)

    do k=1,nCat ; do j=jsd,jed ; do i=isd,ied
      mca0_ice(i,j,k) = CAS%m_ice(i,j,k)
      mca0_snow(i,j,k) = CAS%m_snow(i,j,k)
      mca0_pond(i,j,k) = CAS%m_pond(i,j,k)
    enddo ; enddo ; enddo

    if (merged_cont) then
      call proportionate_continuity(mca_tot(:,:,n-1), uh_tot(:,:,n), vh_tot(:,:,n), &
                                    dt_adv, G, US, IG, CS%continuity_CSp, &
                                    h1=CAS%m_ice,  uh1=uh_ice,  vh1=vh_ice, &
                                    h2=CAS%m_snow, uh2=uh_snow, vh2=vh_snow, &
                                    h3=CAS%m_pond, uh3=uh_pond, vh3=vh_pond)
    else
      call continuity(uc, vc, mca0_ice, CAS%m_ice, uh_ice, vh_ice, dt_adv, &
                      G, US, IG, CS%continuity_CSp, use_h_neg=.true.)
      call continuity(uc, vc, mca0_snow, CAS%m_snow, uh_snow, vh_snow, dt_adv, &
                      G, US, IG, CS%continuity_CSp, masking_uh=uh_ice, masking_vh=vh_ice)
      call continuity(uc, vc, mca0_pond, CAS%m_pond, uh_pond, vh_pond, dt_adv, &
                      G, US, IG, CS%continuity_CSp, masking_uh=uh_ice, masking_vh=vh_ice)
    endif

    call advect_scalar(CAS%mH_ice, mca0_ice, CAS%m_ice, uh_ice, vh_ice, &
                            dt_adv, G, US, IG, CS%SIS_thick_adv_CSp)
    call advect_SIS_tracers(mca0_ice, CAS%m_ice, uh_ice, vh_ice, &
                            dt_adv, G, US, IG, CS%SIS_tr_adv_CSp, TrReg, snow_tr=.false.)
    call advect_SIS_tracers(mca0_snow, CAS%m_snow, uh_snow, vh_snow, &
                            dt_adv, G, US, IG, CS%SIS_tr_adv_CSp, TrReg, snow_tr=.true.)

    ! Accumulated diagnostics
    CAS%dt_sum = CAS%dt_sum + dt_adv
    if (allocated(CAS%uh_sum)) then ; do k=1,nCat ; do j=jsc,jec ; do I=isc-1,iec
      CAS%uh_sum(I,j) = CAS%uh_sum(I,j) + dt_adv * ((uh_pond(I,j,k) + uh_snow(I,j,k)) + uh_ice(I,j,k))
    enddo ; enddo ; enddo ; endif
    if (allocated(CAS%vh_sum)) then ; do k=1,nCat ; do J=jsc-1,jec ; do i=isc,iec
      CAS%vh_sum(i,J) = CAS%vh_sum(i,J) + dt_adv * ((vh_pond(i,J,k) + vh_snow(i,J,k)) + vh_ice(i,J,k))
    enddo ; enddo ; enddo ; endif

    if (CS%bounds_check) then
      write(mesg,'(i4)') n
      call check_SIS_tracer_bounds(TrReg, G, IG, "After advect_SIS_tracers "//trim(mesg))
    endif
  enddo

end subroutine ice_cat_transport

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> finish_ice_transport completes the ice transport and thickness class redistribution
subroutine finish_ice_transport(CAS, IST, TrReg, G, US, IG, dt, CS, rdg_rate)
  type(cell_average_state_type),     intent(inout) :: CAS !< A structure with ocean-cell averaged masses.
  type(ice_state_type),              intent(inout) :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type),           intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),               intent(inout) :: IG  !< The sea-ice specific grid type
  real,                              intent(in)    :: dt  !< The timestep used for ridging [T -> s].
  type(SIS_tracer_registry_type),    pointer       :: TrReg !< The registry of SIS ice and snow tracers.
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_transport_CS),            pointer       :: CS  !< A pointer to the control structure for this module
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: rdg_rate !< The ice ridging rate [T-1 ~> s-1].

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    uf           ! Total zonal fluxes [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    vf           ! Total meridional fluxes [R Z L2 T-1 ~> kg s-1].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)) :: &
    mca0_ice, &  ! The initial mass of ice per unit ocean area in a cell [R Z ~> kg m-2].
    mca0_snow    ! The initial mass of snow per unit ocean area in a cell [R Z ~> kg m-2].
!### These will be needed when the ice ridging is properly implemented.
!  real :: snow2ocn !< Snow dumped into ocean during ridging [R Z ~> kg m-2]
!  real :: enth_snow2ocn !< Mass-averaged enthalpy of the now dumped into ocean during ridging [Q ~> J kg-1]
!  real, dimension(SZI_(G),SZJ_(G)) :: &
!    rdg_open, & ! formation rate of open water due to ridging [T-1 ~> s-1]
!    rdg_vosh    ! rate of ice mass shifted from level to ridged ice [R Z T-1 ~> kg m-2 s-1]
  real :: yr_dt           ! Tne number of timesteps in a year [nondim].
  real, dimension(SZI_(G),SZJ_(G)) :: trans_conv ! The convergence of frozen water transport [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G)) :: ice_cover ! The summed fractional ice concentration [nondim].
  type(EFP_type) :: tot_ice, tot_snow, enth_ice, enth_snow
  real :: I_tot_ice, I_tot_snow
  real :: Idt  ! The reciprocal of the accumulated time [T-1 ~> s-1]
  integer :: i, j, k, isc, iec, jsc, jec, nCat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nCat = IG%CatIce

  !  Convert the ocean-cell averaged properties back into the ice_state_type.
  call cell_ave_state_to_ice_state(CAS, G, US, IG, CS, IST, TrReg)

  if (CS%do_ridging) then
    ! Compress the ice using the ridging scheme taken from the CICE-Icepack module
    call ice_ridging(IST, G, IG, CAS%m_ice, CAS%m_snow, CAS%m_pond, TrReg, US, dt, IST%rdg_rate)
    ! Clean up any residuals
    call compress_ice(IST%part_size, IST%mH_ice, IST%mH_snow, IST%mH_pond, TrReg, G, US, IG, CS, CAS)
  else

  ! Compress the ice where the fractional coverage exceeds 1, starting with the
  ! thinnest category, in what amounts to a minimalist version of a sea-ice
  ! ridging scheme.  A more complete ridging scheme would also compress
  ! thicker ice and allow the fractional ice coverage to drop below 1.
    call compress_ice(IST%part_size, IST%mH_ice, IST%mH_snow, IST%mH_pond, TrReg, G, US, IG, CS, CAS)
  endif
  if (CS%bounds_check) call check_SIS_tracer_bounds(TrReg, G, IG, "After compress_ice")

  if (CS%readjust_categories) then
    call adjust_ice_categories(IST%mH_ice, IST%mH_snow, IST%mH_pond, IST%part_size, &
                               TrReg, G, IG, CS)
    if (CS%bounds_check) call check_SIS_tracer_bounds(TrReg, G, IG, "After adjust_ice_categories")
  endif

  ! Recalculating m_ice and m_snow for consistency when handling tracer
  ! concentrations in massless categories.
  do k=1,nCat ; do j=jsc,jec ; do i=isc,iec
    mca0_ice(i,j,k) = IST%part_size(i,j,k)*IST%mH_ice(i,j,k)
    mca0_snow(i,j,k) = IST%part_size(i,j,k)*IST%mH_snow(i,j,k)
  enddo ; enddo ; enddo
  call set_massless_SIS_tracers(mca0_snow, TrReg, G, IG, compute_domain=.true., do_ice=.false.)
  call set_massless_SIS_tracers(mca0_ice, TrReg, G, IG, compute_domain=.true., do_snow=.false.)

  if (CS%bounds_check) call check_SIS_tracer_bounds(TrReg, G, IG, "SIS_transport set massless 2")

! Niki: TOM does the ridging after redistribute which would need age_ice and IST%rgd_mice below.
!  !  ### THIS IS HARD-CODED ONLY TO WORK WITH 2 LAYERS.
!  !  ### heat_snow AND OTHER TRACERS ARE OMITTED.
!  if (CS%do_ridging) then
!    do j=jsc,jec ; do i=isc,iec
!      if (sum(IST%mH_ice(i,j,:)) > 1.e-10*US%m_to_Z*CS%Rho_ice .and. &
!          sum(IST%part_size(i,j,1:nCat)) > 0.01) then
!        call ice_ridging(nCat, IST%part_size(i,j,:), IST%mH_ice(i,j,:), &
!            IST%mH_snow(i,j,:), CS%Rho_ice, &
!            heat_ice(i,j,:,1), heat_ice(i,j,:,2), & !Niki: Is this correct? Bob: No, 2-layers hard-coded.
!            age_ice(i,j,:), snow2ocn, enth_snow2ocn, rdg_rate(i,j), IST%rgd_mice(i,j,:), &
!            CAS%dt_sum, IG%mH_cat_bound, rdg_open(i,j), rdg_vosh(i,j), US)
!        ! Store the snow mass (and related properties?) that will be passed to the ocean at the
!        ! next opportunity.
!        if (snow2ocn > 0.0) then
!          IST%enth_snow_to_ocn(i,j) = (IST%enth_snow_to_ocn(i,j) * IST%snow_to_ocn(i,j) + &
!                                       enth_snow2ocn * snow2ocn) / &
!                                      (IST%snow_to_ocn(i,j) + snow2ocn)
!          IST%snow_to_ocn(i,j) = IST%snow_to_ocn(i,j) + snow2ocn
!        endif
!      endif
!    enddo ; enddo
!  endif   ! do_ridging

  !   Recalculate IST%part_size(:,:,0) to ensure that the sum of IST%part_size adds up to 1.
  ! Compress_ice should already have taken care of this within the computational
  ! domain, but with a slightly different order of arithmetic.  The max is here
  ! to avoid tiny negative values of order -1e-16 from round-off in the
  ! difference between ice_cover and 1, or to set the fractional open ocean area
  ! to a miniscule positive value so that the ocean-air fluxes are always
  ! calculated.
  ice_cover(:,:) = 0.0
  do k=1,nCat ; do j=jsc,jec ; do i=isc,iec
    ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
  enddo ; enddo ; enddo
  do j=jsc,jec ; do i=isc,iec
    IST%part_size(i,j,0) = max(1.0 - ice_cover(i,j), IG%ocean_part_min)
  enddo ; enddo

  call pass_var(IST%part_size, G%Domain) ! cannot be combined with the three updates below
  call pass_var(IST%mH_pond, G%Domain, complete=.false.)
  call pass_var(IST%mH_snow, G%Domain, complete=.false.)
  call pass_var(IST%mH_ice, G%Domain, complete=.true.)

  if (CS%check_conservation) then
    call get_total_mass(IST, G, US, IG, tot_ice, tot_snow, scale=US%RZ_to_kg_m2)
    call get_total_enthalpy(IST, G, US, IG, enth_ice, enth_snow, scale=US%RZ_to_kg_m2)

    if (is_root_pe()) then
      I_tot_ice  = abs(EFP_to_real(CAS%tot_ice))
      if (I_tot_ice > 0.0) I_tot_ice = 1.0 / I_tot_ice    ! Adcroft's rule inverse.
      I_tot_snow = abs(EFP_to_real(CAS%tot_snow))
      if (I_tot_snow > 0.0) I_tot_snow = 1.0 / I_tot_snow ! Adcroft's rule inverse.
      write(*,'("  Total Ice mass:  ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(tot_ice), EFP_real_diff(tot_ice, CAS%tot_ice), &
        EFP_real_diff(tot_ice, CAS%tot_ice) * I_tot_ice
      write(*,'("  Total Snow mass: ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(tot_snow), EFP_real_diff(tot_snow, CAS%tot_snow), &
        EFP_real_diff(tot_snow, CAS%tot_snow) * I_tot_snow

      I_tot_ice  = abs(EFP_to_real(CAS%enth_ice))
      if (I_tot_ice > 0.0) I_tot_ice = 1.0 / I_tot_ice    ! Adcroft's rule inverse.
      I_tot_snow = abs(EFP_to_real(CAS%enth_snow))
      if (I_tot_snow > 0.0) I_tot_snow = 1.0 / I_tot_snow ! Adcroft's rule inverse.
      write(*,'("  Enthalpy Ice:  ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(enth_ice), EFP_real_diff(enth_ice, CAS%enth_ice), &
        EFP_real_diff(enth_ice, CAS%enth_ice) * I_tot_ice
      write(*,'("  Enthalpy Snow: ",ES24.16,", Error: ",ES12.5," (",ES8.1,")")') &
        EFP_to_real(enth_snow), EFP_real_diff(enth_snow, CAS%enth_snow), &
        EFP_real_diff(enth_snow, CAS%enth_snow) * I_tot_snow
    endif
  endif

  ! Calculate and send transport-related diagnostics.
  Idt = 0.0 ; if (CAS%dt_sum > 0.0) Idt = 1.0 / CAS%dt_sum
  if (CS%id_xprt>0) then
    yr_dt = (8.64e4 * 365.0) * US%s_to_T * Idt
    call get_cell_mass(IST, G, IG, trans_conv)
    do j=jsc,jec ; do i=isc,iec
      trans_conv(i,j) = (trans_conv(i,j) - CAS%mass0(i,j)) * yr_dt
    enddo ; enddo
    call post_SIS_data(CS%id_xprt, trans_conv, CS%diag)
  endif
  if (CS%id_ix_trans>0) then
    do j=jsc,jec ; do I=isc-1,iec ; uf(I,j) = Idt * CAS%uh_sum(I,j) ; enddo ; enddo
    call post_SIS_data(CS%id_ix_trans, uf, CS%diag)
  endif
  if (CS%id_iy_trans>0) then
    do J=jsc-1,jec ; do i=isc,iec ; vf(i,J) = Idt * CAS%vh_sum(i,J) ; enddo ; enddo
    call post_SIS_data(CS%id_iy_trans, vf, CS%diag)
  endif
  if (CS%do_ridging) then
    if (CS%id_rdgr>0 .and. present(rdg_rate)) &
      call post_SIS_data(CS%id_rdgr, rdg_rate, CS%diag)
!    if (CS%id_rdgo>0) call post_SIS_data(CS%id_rdgo, rdg_open, diag)
!    if (CS%id_rdgv>0) then
!      do j=jsc,jec ; do i=isc,iec
!        tmp2d(i,j) = rdg_vosh(i,j) * G%areaT(i,j) * G%mask2dT(i,j)
!      enddo ; enddo
!      call post_SIS_data(CS%id_rdgv, tmp2d, diag)
!    endif
  endif

  if (CS%bounds_check) call check_SIS_tracer_bounds(TrReg, G, IG, "At end of SIS_transport")

end subroutine finish_ice_transport


!>  Determine the whole-cell averaged mass of snow and ice by thickness category based
!! on the information in the ice state type.
subroutine ice_state_to_cell_ave_state(IST, G, US, IG, CS, CAS)
  type(ice_state_type),          intent(in)    :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type),       intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),         intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),           intent(in)    :: IG  !< The sea-ice specific grid type
  type(SIS_transport_CS),        pointer       :: CS  !< A pointer to the control structure for this module
  type(cell_average_state_type), intent(inout) :: CAS !< A structure with ocean-cell averaged masses.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: ice_cover ! The summed fractional ice concentration [nondim].
  real, dimension(SZI_(G),SZJ_(G)) :: mHi_avg   ! The average ice mass-thickness [R Z ~> kg m-2].
  integer :: i, j, k, isc, iec, jsc, jec, nCat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nCat = IG%CatIce

  CAS%m_ice(:,:,:) = 0.0 ; CAS%m_snow(:,:,:) = 0.0 ; CAS%m_pond(:,:,:) = 0.0 ; CAS%mH_ice(:,:,:) = 0.0
  ice_cover(:,:) = 0.0 ; mHi_avg(:,:) = 0.0
  !$OMP parallel do default(shared)
  do j=jsc,jec
    do k=1,nCat ; do i=isc,iec
      if (IST%mH_ice(i,j,k)>0.0) then
        CAS%m_ice(i,j,k)  = IST%part_size(i,j,k) * IST%mH_ice(i,j,k)
        CAS%m_snow(i,j,k) = IST%part_size(i,j,k) * IST%mH_snow(i,j,k)
        CAS%m_pond(i,j,k) = IST%part_size(i,j,k) * IST%mH_pond(i,j,k)
        CAS%mH_ice(i,j,k) = IST%mH_ice(i,j,k)
        ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
        mHi_avg(i,j) = mHi_avg(i,j) + CAS%m_ice(i,j,k)
      else
        if (IST%part_size(i,j,k)*IST%mH_snow(i,j,k) > 0.0) then
          call SIS_error(FATAL, "Input to SIS_transport, non-zero snow mass rests atop no ice.")
        endif
        if (IST%part_size(i,j,k)*IST%mH_pond(i,j,k) > 0.0) then
          call SIS_error(FATAL, "Input to SIS_transport, non-zero pond mass rests atop no ice.")
        endif
        CAS%m_ice(i,j,k) = 0.0 ; CAS%m_snow(i,j,k) = 0.0 ; CAS%m_pond(i,j,k) = 0.0
      endif
    enddo ; enddo
    do i=isc,iec ; if (ice_cover(i,j) > 0.0) then
      mHi_avg(i,j) = mHi_avg(i,j) / ice_cover(i,j)
    endif ; enddo

    !   Handle massless categories.
    do k=1,nCat ; do i=isc,iec
      if (CAS%m_ice(i,j,k)<=0.0 .and. (G%mask2dT(i,j) > 0.0)) then
        if (mHi_avg(i,j) <= IG%mH_cat_bound(k)) then
          CAS%mH_ice(i,j,k) = IG%mH_cat_bound(k)
        elseif (mHi_avg(i,j) >= IG%mH_cat_bound(k+1)) then
          CAS%mH_ice(i,j,k) = IG%mH_cat_bound(k+1)
        else
          CAS%mH_ice(i,j,k) = mHi_avg(i,j)
        endif
      endif
    enddo ; enddo
  enddo

  ! Handle diagnostics
  CAS%dt_sum = 0.0
  if (allocated(CAS%mass0))  call get_cell_mass(IST, G, IG, CAS%mass0)
  if (allocated(CAS%uh_sum)) CAS%uh_sum(:,:) = 0.0
  if (allocated(CAS%vh_sum)) CAS%vh_sum(:,:) = 0.0

  if (CS%check_conservation) then ! mw/new - need to update this for pond ?
    call get_total_mass(IST, G, US, IG, CAS%tot_ice, CAS%tot_snow, scale=US%RZ_to_kg_m2)
    call get_total_enthalpy(IST, G, US, IG, CAS%enth_ice, CAS%enth_snow, scale=US%RZ_to_kg_m2)
  endif

end subroutine ice_state_to_cell_ave_state

!> Convert the ocean-cell averaged properties back into the ice_state_type.
subroutine cell_ave_state_to_ice_state(CAS, G, US, IG, CS, IST, TrReg)
  type(cell_average_state_type),  intent(inout) :: CAS !< A structure with ocean-cell averaged masses.
  type(SIS_hor_grid_type),        intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),          intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),            intent(in)    :: IG  !< The sea-ice specific grid type
  type(SIS_transport_CS),         pointer       :: CS  !< A pointer to the control structure for this module
  type(ice_state_type),           intent(inout) :: IST !< A type describing the state of the sea ice
  type(SIS_tracer_registry_type), pointer       :: TrReg !< The registry of SIS ice and snow tracers.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: ice_cover ! The summed fractional ice concentration [nondim].
  real :: mass_neglect    ! A negligible mass per unit area [R Z ~> kg m-2].
  real :: L_to_H
  integer :: i, j, k, isc, iec, jsc, jec, nCat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nCat = IG%CatIce
  mass_neglect = US%kg_m3_to_R*US%m_to_Z*1.0e-60

  ! Ensure that CAS%mH_ice(i,j,1) >= IG%mH_cat_bound(1).
  do j=jsc,jec ; do i=isc,iec
    if ((CAS%m_ice(i,j,1) > 0.0) .and. (CAS%mH_ice(i,j,1) < IG%mH_cat_bound(1))) &
      CAS%mH_ice(i,j,1) = IG%mH_cat_bound(1)
  enddo ; enddo

  ! Convert CAS%m_ice and CAS%m_snow back to IST%part_size and IST%mH_snow.
  ice_cover(:,:) = 0.0
  L_to_H = US%L_to_Z * CS%Rho_ice
  !$OMP parallel do default(shared)
  do j=jsc,jec ; do k=1,nCat ; do i=isc,iec
    if (CAS%m_ice(i,j,k) > 0.0) then
      !### This is a simplified version of the test, but it could rarely change answers at roundoff.
      ! if (CS%roll_factor * CAS%mH_ice(i,j,k)**3 > L_to_H**2 * (CAS%m_ice(i,j,k)*G%areaT(i,j))) then
      if (CS%roll_factor * (CAS%mH_ice(i,j,k)/(US%L_to_Z*CS%Rho_Ice))**3 > &
          (CAS%m_ice(i,j,k)/(US%L_to_Z*CS%Rho_Ice))*G%areaT(i,j)) then
        ! This ice is thicker than it is wide even if all the ice in a grid cell is collected
        ! into a single cube, so it will roll.  Any snow on top will simply be redistributed
        ! into a thinner layer, although it should probably be dumped into the ocean.  Rolling
        ! makes the ice thinner so that it melts faster, but it should never be made thinner
        ! than IG%mH_cat_bound(1).
        CAS%mH_ice(i,j,k) = max(IG%mH_cat_bound(1), L_to_H * &
             sqrt((CAS%m_ice(i,j,k)*G%areaT(i,j)) / (CS%roll_factor * CAS%mH_ice(i,j,k)) ))
      endif

      ! Make sure that CAS%mH_ice(i,j,k) > IG%mH_cat_bound(1).
      if (CAS%mH_ice(i,j,k) < IG%mH_cat_bound(1)) CAS%mH_ice(i,j,k) = IG%mH_cat_bound(1)

      IST%part_size(i,j,k) = CAS%m_ice(i,j,k) / CAS%mH_ice(i,j,k)
      IST%mH_snow(i,j,k) = CAS%mH_ice(i,j,k) * (CAS%m_snow(i,j,k) / CAS%m_ice(i,j,k))
      IST%mH_pond(i,j,k) = CAS%mH_ice(i,j,k) * (CAS%m_pond(i,j,k) / CAS%m_ice(i,j,k))
      IST%mH_ice(i,j,k) = CAS%mH_ice(i,j,k)
      ice_cover(i,j) = ice_cover(i,j) + IST%part_size(i,j,k)
    else
      IST%part_size(i,j,k) = 0.0 ; IST%mH_ice(i,j,k) = 0.0
      if (CAS%m_snow(i,j,k) > mass_neglect) &
        call SIS_error(FATAL, &
          "Positive CAS%m_snow values should not exist without ice.")
      if (CAS%m_pond(i,j,k) > mass_neglect ) &
        call SIS_error(FATAL, &
          "Something needs to be done with positive CAS%m_pond values without ice.")
      IST%mH_snow(i,j,k) = 0.0 ; IST%mH_pond(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo
  do j=jsc,jec ; do i=isc,iec
    IST%part_size(i,j,0) = 1.0-ice_cover(i,j)
  enddo ; enddo

end subroutine cell_ave_state_to_ice_state

!> adjust_ice_categories moves mass between thickness categories if it is thinner or
!! thicker than the bounding limits of each category.
subroutine adjust_ice_categories(mH_ice, mH_snow, mH_pond, part_sz, TrReg, G, IG, CS)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(inout) :: mH_ice  !< The mass per unit area of the ice
                                                !! in each category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(inout) :: mH_snow !< The mass per unit area of the snow
                                                !! atop the ice in each category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                           intent(inout) :: mH_pond !< The mass per unit area of the pond
                                                !! on the ice in each category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(IG)), &
                           intent(inout) :: part_sz !< The fractional ice concentration
                                                !! within a cell in each thickness
                                                !! camcategory [nondim], 0-1.
  type(SIS_tracer_registry_type), &
                           pointer       :: TrReg !< The registry of SIS ice and snow tracers.
  type(SIS_transport_CS),  pointer       :: CS  !< A pointer to the control structure for this module

!   This subroutine moves mass between thickness categories if it is thinner or
! thicker than the bounding limits of each category.

  ! Local variables
  real :: mca_trans  ! The cell-averaged ice mass transferred between categories [R Z ~> kg m-2].
  real :: part_trans ! The fractional area transferred between categories [nondim].
  real :: snow_trans ! The cell-averaged snow transferred between categories [R Z ~> kg m-2].
  real :: pond_trans ! The cell-averaged pond transferred between categories [R Z ~> kg m-2].
  real :: I_mH_lim1  ! The inverse of the lower thickness limit [R-1 Z-1 ~> m2 kg-1].
  real, dimension(SZI_(G),SZCAT_(IG)) :: &
    ! The mass of snow, pond and ice per unit total area in a cell [R Z ~> kg m-2].
    ! "mca" stands for "mass cell averaged"
    mca_ice, mca_snow, mca_pond, &
    ! Initial ice, snow and pond masses per unit cell area [R Z ~> kg m-2].
    mca0_ice, mca0_snow, mca0_pond, &
    ! Cross-catagory transfers of ice, snow and pond mass [R Z ~> kg m-2].
    trans_ice, trans_snow, trans_pond
  real, dimension(SZI_(G)) :: ice_cover ! The summed fractional ice coverage [nondim].
  logical :: do_any, do_j(SZJ_(G)), resum_cat(SZI_(G), SZJ_(G))
  integer :: i, j, k, m, is, ie, js, je, nCat
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  nCat = IG%CatIce

  I_mH_lim1 = 1.0 / IG%mH_cat_bound(1)

  resum_cat(:,:) = .false.

  ! Zero out the part_size of any massless categories.
  do k=1,nCat ; do j=js,je ; do i=is,ie ; if (mH_ice(i,j,k) <= 0.0) then
    if (mH_ice(i,j,k) < 0.0) then
      call SIS_error(FATAL, "Input to adjust_ice_categories, negative ice mass.")
    endif
    if (mH_snow(i,j,k) > 0.0) then
      call SIS_error(FATAL, "Input to adjust_ice_categories, non-zero snow mass rests atop no ice.")
    endif
    if (mH_pond(i,j,k) > 0.0) then
      call SIS_error(FATAL, "Input to adjust_ice_categories, non-zero pond mass rests atop no ice.")
    endif
    if (part_sz(i,j,k) > 0.0) resum_cat(i,j) = .true.
    part_sz(i,j,k) = 0.0
  endif ; enddo ; enddo ; enddo

  ! Only work on rows that have any sea ice at all.
  do j=js,je
    do_j(j) = .false.
    do k=1,nCat
      do i=is,ie ; if (part_sz(i,j,k)*mH_ice(i,j,k) > 0.0) then
        do_j(j) = .true. ; exit
      endif ; enddo
      if (do_j(j)) exit
    enddo
  enddo

  do j=js,je ; if (do_j(j)) then
    do k=1,nCat ; do i=is,ie
      mca_ice(i,k) = part_sz(i,j,k)*mH_ice(i,j,k)
      mca_snow(i,k) = part_sz(i,j,k)*mH_snow(i,j,k)
      mca_pond(i,k) = part_sz(i,j,k)*mH_pond(i,j,k)

      mca0_ice(i,k) = mca_ice(i,k)
      mca0_snow(i,k) = mca_snow(i,k)
      mca0_pond(i,k) = mca_pond(i,k)
    enddo ; enddo
    trans_ice(:,:) = 0.0 ; trans_snow(:,:) = 0.0 ; trans_pond(:,:) = 0.0
    do_any = .false.

    do k=1,nCat-1 ; do i=is,ie
      if ((mca_ice(i,k) > 0.0) .and. (mH_ice(i,j,k) > IG%mH_cat_bound(k+1))) then
        ! Move some or all of the ice to a thicker category.
        ! For now move all of it.
        mca_trans = mca_ice(i,k)
        part_trans = part_sz(i,j,k) ! * (mca_trans / mca_ice) * (mH_ice / h_trans)
        snow_trans = mca_snow(i,k) ! * (part_trans / part_sz) = 1
        pond_trans = mca_pond(i,k) ! * (part_trans / part_sz) = 1

        trans_ice(i,K) = mca_trans
        trans_snow(i,K) = snow_trans
        trans_pond(i,K) = pond_trans
        do_any = .true.

        ! Use area-weighted remapped thicknesses so that the total ice area and
        ! mass are both conserved in the remapping operation.  Using a mass-
        ! weighted average thickness instead would cause the ice to systematically
        ! contract in area.
        mH_ice(i,j,k+1) = (part_trans*mH_ice(i,j,k) + &
                           part_sz(i,j,k+1)*mH_ice(i,j,k+1)) / &
                          (part_trans + part_sz(i,j,k+1))
        ! h should be the first thing to correct via a non-constant profile, and
        ! can be improved independent of T & S.
        mH_ice(i,j,k) = IG%mH_cat_bound(k+1)
        part_sz(i,j,k+1) = part_sz(i,j,k+1) + part_trans
        part_sz(i,j,k) = part_sz(i,j,k) - part_trans

        mca_ice(i,k+1) = mca_ice(i,k+1) + mca_trans
        mca_ice(i,k) = mca_ice(i,k) - mca_trans

        mca_snow(i,k+1) = mca_snow(i,k+1) + snow_trans
        mca_snow(i,k) = mca_snow(i,k) - snow_trans

        mH_snow(i,j,k) = 0.0 ; mH_snow(i,j,k+1) = 0.0
        if (part_sz(i,j,k)>0.0) mH_snow(i,j,k) = mca_snow(i,k) / part_sz(i,j,k)
        if (part_sz(i,j,k+1)>0.0) mH_snow(i,j,k+1) = mca_snow(i,k+1) / part_sz(i,j,k+1)

        mca_pond(i,k+1) = mca_pond(i,k+1) + pond_trans
        mca_pond(i,k) = mca_pond(i,k) - pond_trans

        mH_pond(i,j,k) = 0.0 ; mH_pond(i,j,k+1) = 0.0
        if (part_sz(i,j,k)>0.0) mH_pond(i,j,k) = mca_pond(i,k) / part_sz(i,j,k)
        if (part_sz(i,j,k+1)>0.0) mH_pond(i,j,k+1) = mca_pond(i,k+1) / part_sz(i,j,k+1)
      endif
    enddo ; enddo

    if (do_any) then ! no pond tracers yet - mw
      call advect_tracers_thicker(mca0_ice, trans_ice, G, IG, CS%SIS_tr_adv_CSp, &
                                  TrReg, .false., j, is, ie)
      call advect_tracers_thicker(mca0_snow, trans_snow, G, IG, CS%SIS_tr_adv_CSp, &
                                  TrReg, .true., j, is, ie)
    endif

    do k=1,nCat ; do i=is,ie
      mca0_ice(i,k) = mca_ice(i,k)
      mca0_snow(i,k) = mca_snow(i,k)
      mca0_pond(i,k) = mca_pond(i,k)
    enddo ; enddo
    trans_ice(:,:) = 0.0 ; trans_snow(:,:) = 0.0 ; trans_pond(:,:) = 0.0
    do_any = .false.

    do k=nCat,2,-1 ; do i=is,ie
      if ((mca_ice(i,k) > 0.0) .and. (mH_ice(i,j,k) < IG%mH_cat_bound(k))) then
        ! Move some or all of the ice to a thinner category.
        ! For now move all of it.
        mca_trans = mca_ice(i,k)
        part_trans = part_sz(i,j,k) ! * (mca_trans / mca_ice) * (mH_ice / h_trans)
        snow_trans = mca_snow(i,k) ! * (part_trans / part_sz) = 1
        pond_trans = mca_pond(i,k) ! * (part_trans / part_sz) = 1

        do_any = .true.
        trans_ice(i,K-1) = -mca_trans  ! Note the shifted index conventions!
        trans_snow(i,K-1) = -snow_trans
        trans_pond(i,K-1) = -pond_trans

        ! Use area-weighted remapped thicknesses so that the total ice area and
        ! mass are both conserved in the remapping operation.
        mH_ice(i,j,k-1) = (part_trans*mH_ice(i,j,k) + &
                           part_sz(i,j,k-1)*mH_ice(i,j,k-1)) / &
                          (part_trans + part_sz(i,j,k-1))

        ! h should be the first thing to correct via a non-constant profile, and
        ! can be improved independently from T & S.
        mH_ice(i,j,k) = IG%mH_cat_bound(k)

        part_sz(i,j,k-1) = part_sz(i,j,k-1) + part_trans
        part_sz(i,j,k) = part_sz(i,j,k) - part_trans

        mca_ice(i,k-1) = mca_ice(i,k-1) + mca_trans
        mca_ice(i,k) = mca_ice(i,k) - mca_trans

        mca_snow(i,k-1) = mca_snow(i,k-1) + snow_trans
        mca_snow(i,k) = mca_snow(i,k) - snow_trans

        mH_snow(i,j,k) = 0.0 ; mH_snow(i,j,k-1) = 0.0
        if (part_sz(i,j,k)>0.0) mH_snow(i,j,k) = mca_snow(i,k) / part_sz(i,j,k)
        if (part_sz(i,j,k-1)>0.0) mH_snow(i,j,k-1) = mca_snow(i,k-1) / part_sz(i,j,k-1)

        mca_pond(i,k-1) = mca_pond(i,k-1) + pond_trans
        mca_pond(i,k) = mca_pond(i,k) - pond_trans

        mH_pond(i,j,k) = 0.0 ; mH_pond(i,j,k-1) = 0.0
        if (part_sz(i,j,k)>0.0) mH_pond(i,j,k) = mca_pond(i,k) / part_sz(i,j,k)
        if (part_sz(i,j,k-1)>0.0) mH_pond(i,j,k-1) = mca_pond(i,k-1) / part_sz(i,j,k-1)
      endif
    enddo ; enddo

    if (do_any) then
      call advect_tracers_thicker(mca0_ice, trans_ice, G, IG, CS%SIS_tr_adv_CSp, &
                                  TrReg, .false., j, is, ie)
      call advect_tracers_thicker(mca0_snow, trans_snow, G, IG, CS%SIS_tr_adv_CSp, &
                                  TrReg, .true., j, is, ie)
    endif

    ! Compress the ice in category 1 if it is thinner than the minimum.  This
    ! does not affect any tracer concentrations.
    if (IG%mH_cat_bound(1) > 0.0) then
      do i=is,ie
        if ((mH_ice(i,j,1)*part_sz(i,j,1) > 0.0) .and. &
            (mH_ice(i,j,1) < IG%mH_cat_bound(1))) then
          ! Compress the ice in this category to the minimum thickness.
          part_sz(i,j,1) = part_sz(i,j,1) * (mH_ice(i,j,1) * I_mH_lim1)
          mH_snow(i,j,1) = mH_snow(i,j,1) * (IG%mH_cat_bound(1) / mH_ice(i,j,1))
          mH_pond(i,j,1) = mH_pond(i,j,1) * (IG%mH_cat_bound(1) / mH_ice(i,j,1))
          ! This is equivalent to mH_snow(i,j,1) = mca_snow(i,1) / part_sz(i,j,1)
          mH_ice(i,j,1) = IG%mH_cat_bound(1)
          resum_cat(i,j) = .true.
        endif
      enddo
    endif
    if (CS%ice_cover_discard > 0.0) then
      do i=is,ie ; if ((part_sz(i,j,1) > 0.0) .and. (part_sz(i,j,1) < CS%ice_cover_discard)) then
        part_sz(i,j,1) = 0.0
        resum_cat(i,j) = .true.
      endif ; enddo
    endif

  endif ; enddo  ! j-loop and do_j

  if (.not.CS%inconsistent_cover_bug) then
    do j=js,je
      do_any = .false.
      do i=is,ie
        ice_cover(i) = 0.0
        if (resum_cat(i,j)) do_any = .true.
      enddo
      if (.not.do_any) cycle
      do k=1,nCat ; do i=is,ie
        ice_cover(i) = ice_cover(i) + part_sz(i,j,k)
      enddo ; enddo
      do i=is,ie ; if (resum_cat(i,j)) then
        part_sz(i,j,0) = max(1.0 - ice_cover(i), IG%ocean_part_min)
      endif ; enddo
    enddo
  endif

end subroutine adjust_ice_categories

!> compress_ice compresses the ice, starting with the thinnest category, if the total fractional
!! ice coverage exceeds 1.  It is assumed at the start that the sum over all categories (including
!! ice free) of part_sz is 1, but that the part_sz of the ice free category may be negative to make
!! this so.  In this routine, the mass (volume) is conserved, while the fractional coverage is
!! solved for, while the new thicknesses are diagnosed.
subroutine compress_ice(part_sz, mH_ice, mH_snow, mH_pond, TrReg, G, US, IG, CS, CAS)
  type(SIS_hor_grid_type),           intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),               intent(in)    :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G),0:SZCAT_(IG)), &
                                     intent(inout) :: part_sz !< The fractional ice concentration
                                                          !! within a cell in each thickness
                                                          !! category [nondim], 0-1.
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                                     intent(inout) :: mH_ice  !< The mass per unit area of the ice
                                                          !! in each category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                                     intent(inout) :: mH_snow !< The mass per unit area of the snow
                                                          !! atop the ice in each category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), &
                                     intent(inout) :: mH_pond !< The mass per unit area of the pond
                                                          !! on the ice in each category [R Z ~> kg m-2].
  type(SIS_tracer_registry_type),    pointer       :: TrReg !< The registry of SIS ice and snow tracers.
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors
  type(SIS_transport_CS),            pointer       :: CS  !< A pointer to the control structure for this module
  type(cell_average_state_type), optional, intent(in) :: CAS !< A structure with ocean-cell averaged masses.
!   This subroutine compresses the ice, starting with the thinnest category, if
! the total fractional ice coverage exceeds 1.  It is assumed at the start that
! the sum over all categories (including ice free) of part_sz is 1, but that the
! part_sz of the ice free category may be negative to make this so.  In this
! routine, the mass (volume) is conserved, while the fractional coverage is
! solved for, while the new thicknesses are diagnosed.

!   This subroutine is effectively a minimalist version of a sea-ice ridging
! scheme.  A more complete ridging scheme would also compress thicker ice and
! allow the fractional ice coverage to drop below 1.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: excess_cover
  real :: compression_ratio
  real :: Icompress_here
  real :: mca_old
!  real :: Imca_new
  real :: mass_neglect
  real :: part_trans ! The fractional area transferred into a thicker category [nondim].
  real, dimension(SZI_(G),SZCAT_(IG)) :: &
    m0_ice, &  ! The initial mass per unit grid-cell area of ice in each category [R Z ~> kg m-2].
    m0_snow, & ! The initial mass per unit grid-cell area of snow in each category [R Z ~> kg m-2].
    m0_pond    ! The initial mass per unit grid-cell pond melt water in each category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZCAT_(IG)) :: &
    trans_ice, trans_snow, trans_pond ! The masses transferred into the next thicker category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZCAT_(IG)) :: mca_ice  ! The mass per unit grid-cell area
                                                  ! of the ice in each category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZCAT_(IG)) :: mca_snow ! The mass per unit grid-cell area
                                                  ! of the snow atop the ice in each category [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZCAT_(IG)) :: mca_pond ! The mass per unit grid-cell area of the melt
                                                  ! ponds atop the ice in each category [R Z ~> kg m-2].
  logical :: do_any, do_j(SZJ_(G))
  character(len=200) :: mesg
  integer :: i, j, k, m, isc, iec, jsc, jec, nCat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  nCat = IG%CatIce

  !### Consider recalculating mca_ice and mca_snow here, as it is not reused again outside.

  ! 1.0e-40 kg/m2 is roughly the mass of one molecule of water divided by the surface area of the Earth.
  mass_neglect = US%kg_m3_to_R*US%m_to_Z*1.0e-60

  do_j(:) = .false.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,do_j,G,IG,part_sz,excess_cover, &
!$OMP                                  mH_ice,mH_snow,mH_pond,&
!$OMP                                  mass_neglect,CS,CAS,TrReg,nCat) &
!$OMP                          private(m0_ice,do_any,m0_snow,trans_ice,trans_snow, &
!$OMP                                  m0_pond,trans_pond,compression_ratio,Icompress_here, &
!$OMP                                  mca_ice,mca_snow,mca_pond,mca_old,part_trans)
  do j=jsc,jec
    do i=isc,iec
      if (part_sz(i,j,0) < 0.0) then
        excess_cover(i,j) = -part_sz(i,j,0) ; part_sz(i,j,0) = 0.0
        do_j(j) = .true.
      else
        excess_cover(i,j) = 0.0
      endif
    enddo

    if (do_j(j)) then
      if (present(CAS)) then
        do k=1,nCat ; do i=isc,iec
          m0_ice(i,k) = CAS%m_ice(i,j,k)
          m0_snow(i,k) = CAS%m_snow(i,j,k)
          m0_pond(i,k) = CAS%m_pond(i,j,k)
          mca_ice(i,k) = m0_ice(i,k) ; mca_snow(i,k) = m0_snow(i,k) ; mca_pond(i,k) = m0_pond(i,k)
        enddo ; enddo
      else  ! This is mathematically equivalent ot the code above, but can differ at roundoff.
        do k=1,nCat ; do i=isc,iec
          m0_ice(i,k) = part_sz(i,j,k) * mH_ice(i,j,k)
          m0_snow(i,k) = part_sz(i,j,k) * mH_snow(i,j,k)
          m0_pond(i,k) = part_sz(i,j,k) * mH_pond(i,j,k)
          mca_ice(i,k) = m0_ice(i,k) ; mca_snow(i,k) = m0_snow(i,k) ; mca_pond(i,k) = m0_pond(i,k)
        enddo ; enddo
      endif

      trans_ice(:,:) = 0.0 ; trans_snow(:,:) = 0.0 ; trans_pond(:,:) = 0.0
      do_any = .false.
      do k=1,nCat-1 ; do i=isc,iec
        if ((excess_cover(i,j) > 0.0) .and. (mca_ice(i,k) > 0.0)) then
          compression_ratio = mH_ice(i,j,k) / IG%mH_cat_bound(k+1)
          if (part_sz(i,j,k)*(1.0-compression_ratio) >= excess_cover(i,j)) then
            ! This category is compacted, but not to the point that it needs to
            ! be transferred to a thicker layer.
            Icompress_here = part_sz(i,j,k) / (part_sz(i,j,k) - excess_cover(i,j))
            mH_ice(i,j,k) = mH_ice(i,j,k) * Icompress_here
            mH_snow(i,j,k) = mH_snow(i,j,k) * Icompress_here
            mH_pond(i,j,k) = mH_pond(i,j,k) * Icompress_here
            part_sz(i,j,k) = part_sz(i,j,k) - excess_cover(i,j)
            excess_cover(i,j) = 0.0
          else
            ! Mass from this category needs to be transferred to the next thicker
            ! category after being compacted to thickness IG%mH_cat_bound(k+1).
            excess_cover(i,j) = excess_cover(i,j) - part_sz(i,j,k)*(1.0-compression_ratio)

            if (mca_ice(i,k) > mass_neglect) then
              part_sz(i,j,k+1) = part_sz(i,j,k+1) + part_sz(i,j,k)*compression_ratio

              mca_old = mca_ice(i,k+1)
              trans_ice(i,K) = mca_ice(i,k) ; do_any = .true.
              mca_ice(i,k+1) = mca_ice(i,k+1) + mca_ice(i,k)

              if (part_sz(i,j,k+1) > 1.0e-60) then ! For 32-bit reals this should be 1.0e-30.
                ! This is the usual case, and underflow is no problem.
                mH_ice(i,j,k+1) = mca_ice(i,k+1) / part_sz(i,j,k+1)
              elseif (trans_ice(i,K) > mca_old) then
                ! Set the ice category's thickness to its lower bound.
                part_sz(i,j,k+1) = mca_ice(i,k+1) / IG%mH_cat_bound(k+1)
                mH_ice(i,j,k+1) = IG%mH_cat_bound(k+1)
              else  ! Keep the ice category's thickness at its previous value.
                part_sz(i,j,k+1) = mca_ice(i,k+1) / mH_ice(i,j,k+1)
              endif

              if (mca_snow(i,k) > 0.0) then
                trans_snow(i,K) = mca_snow(i,k)
                mca_snow(i,k+1) = mca_snow(i,k+1) + mca_snow(i,k)
              endif
              mH_snow(i,j,k+1) = mca_snow(i,k+1) / part_sz(i,j,k+1)

              if (mca_pond(i,k) > 0.0) then
                trans_pond(i,K) = mca_pond(i,k)
                mca_pond(i,k+1) = mca_pond(i,k+1) + mca_pond(i,k)
              endif
              mH_pond(i,j,k+1) = mca_pond(i,k+1) / part_sz(i,j,k+1)
            endif

            mca_ice(i,k) = 0.0 ; mca_snow(i,k) = 0.0 ; mca_pond(i,k) = 0.0
            mH_ice(i,j,k) = 0.0 ; mH_snow(i,j,k) = 0.0 ; mH_pond(i,j,k) = 0.0
            part_sz(i,j,k) = 0.0
          endif
        endif
      enddo ; enddo

      if (do_any) then
        ! The following subroutine calls are not thread-safe. There is a pointer in the subroutine
        ! (Tr) that could be redirected from underneath a thread when another goes in.
        !$OMP CRITICAL (safepointer)
        call advect_tracers_thicker(m0_ice, trans_ice, G, IG, CS%SIS_tr_adv_CSp, &
                                    TrReg, .false., j, isc, iec)
        call advect_tracers_thicker(m0_snow, trans_snow, G, IG, CS%SIS_tr_adv_CSp, &
                                    TrReg, .true., j, isc, iec)
        !$OMP END CRITICAL (safepointer)
      endif

      k=nCat
      do i=isc,iec
        if (excess_cover(i,j) > 0.0) then
          if ((part_sz(i,j,k) <= 1.0) .and. &
              (excess_cover(i,j) > 2.0*nCat*epsilon(Icompress_here))) then
            call SIS_error(FATAL, &
                "Category CatIce part_sz inconsistent with excess cover.")
          endif
          Icompress_here = part_sz(i,j,k) / (part_sz(i,j,k) - excess_cover(i,j))
          mH_ice(i,j,k) = mH_ice(i,j,k) * Icompress_here
          mH_snow(i,j,k) = mH_snow(i,j,k) * Icompress_here
          mH_pond(i,j,k) = mH_pond(i,j,k) * Icompress_here
          part_sz(i,j,k) = part_sz(i,j,k) - excess_cover(i,j)
          excess_cover(i,j) = 0.0
        endif
      enddo

    ! if (CS%check_conservation) then
    !   ! Check for consistency between mca_ice, mH_ice, and part_sz.
    !   do k=1,nCat ; do i=isc,iec
    !     if ((mca_ice(i,k) == 0.0) .and. (mH_ice(i,j,k)*part_sz(i,j,k) /= 0.0)) then
    !       write(mesg,'("Compress mismatch at ",3(i8),": mca, h, part, hxp = zero, ",3(1pe15.6))') &
    !          i, j, k, mH_ice(i,j,k), part_sz(i,j,k), mH_ice(i,j,k)*part_sz(i,j,k)
    !       call SIS_error(WARNING, mesg, all_print=.true.)
    !     endif
    !     if (abs(mca_ice(i,k) - mH_ice(i,j,k)*part_sz(i,j,k)) > 1e-12*mca_ice(i,k)) then
    !       write(mesg,'("Compress mismatch at ",3(i8),": mca, h, part, hxp = ",4(1pe15.6))') &
    !          i, j, k, mca_ice(i,k), mH_ice(i,j,k), part_sz(i,j,k), mH_ice(i,j,k)*part_sz(i,j,k)
    !       call SIS_error(WARNING, mesg, all_print=.true.)
    !     endif
    !   enddo ; enddo
    ! endif

    endif  ! Any compression occurs in this j-loop.
  enddo ! j-loop


end subroutine compress_ice

!> get_total_mass determines the globally integrated mass of snow and ice
subroutine get_total_mass(IST, G, US, IG, tot_ice, tot_snow, tot_pond, scale)
  type(ice_state_type),    intent(in)    :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(in)    :: G   !< The horizontal grid type
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  type(EFP_type),          intent(out)   :: tot_ice  !< The globally integrated total ice [kg].
  type(EFP_type),          intent(out)   :: tot_snow !< The globally integrated total snow [kg].
  type(EFP_type),optional, intent(out)   :: tot_pond !< The globally integrated total snow [kg].
  real,          optional, intent(in)    :: scale !< A scaling factor from H to the desired units.

  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: sum_ice, sum_snow, sum_pond
  real :: H_to_units ! A conversion factor from H to the desired output units.
  real :: total
  integer :: i, j, k, m, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  H_to_units = US%RZ_to_kg_m2*US%L_to_m**2 ; if (present(scale)) H_to_units = scale*US%L_to_m**2

  sum_ice(:,:) = 0.0
  sum_snow(:,:) = 0.0
  do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_ice(i,j) = sum_ice(i,j) + G%areaT(i,j) * &
                       (IST%part_size(i,j,k) * (H_to_units*IST%mH_ice(i,j,k)))
    sum_snow(i,j) = sum_snow(i,j) + G%areaT(i,j) * &
                       (IST%part_size(i,j,k) * (H_to_units*IST%mH_snow(i,j,k)))
    if (present(tot_pond)) &
      sum_pond(i,j) = sum_pond(i,j) + G%areaT(i,j) * &
                       (IST%part_size(i,j,k) * (H_to_units*IST%mH_pond(i,j,k)))
  enddo ; enddo ; enddo

  total = reproducing_sum(sum_ice, EFP_sum=tot_ice)
  total = reproducing_sum(sum_snow, EFP_sum=tot_snow)
  if (present(tot_pond)) total = reproducing_sum(sum_pond, EFP_sum=tot_pond)

end subroutine get_total_mass

!> get_cell_mass determines the integrated mass of snow and ice in each cell
subroutine get_cell_mass(IST, G, IG, cell_mass, scale)
  type(ice_state_type),             intent(in)  :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  type(ice_grid_type),              intent(in)  :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: cell_mass !< The total amount of ice and snow [R Z ~> kg m-2].
  real,                   optional, intent(in)  :: scale !< A scaling factor from H to the desired units.

  real :: H_to_units ! A conversion factor from H to the desired output units.
  integer :: i, j, k, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  H_to_units = 1.0 ; if (present(scale)) H_to_units = scale

  cell_mass(:,:) = 0.0
  do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
    cell_mass(i,j) = cell_mass(i,j) + IST%part_size(i,j,k) * H_to_units * &
                          ((IST%mH_snow(i,j,k) + IST%mH_pond(i,j,k)) + IST%mH_ice(i,j,k))
  enddo ; enddo ; enddo

end subroutine get_cell_mass

subroutine cell_mass_from_CAS(CAS, G, IG, mca, scale)
  type(cell_average_state_type),    intent(in)  :: CAS !< A structure with ocean-cell averaged masses by
                                                       !! category and phase of water.
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  type(ice_grid_type),              intent(in)  :: IG  !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: mca !< The combined mass of ice, snow, and
                                                       !! melt pond water in each cell [R Z ~> kg m-2].
  real,                   optional, intent(in)  :: scale !< A scaling factor from H to the desired units.

  real :: H_to_units ! A conversion factor from H to the desired output units.
  integer :: i, j, k, isc, iec, jsc, jec, nCat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; nCat = IG%CatIce

  H_to_units = 1.0 ; if (present(scale)) H_to_units = scale

  do j=jsc,jec ; do i=isc,iec ; mca(i,j) = 0.0 ; enddo ; enddo
  do k=1,nCat ; do j=jsc,jec ; do i=isc,iec
    mca(i,j) = mca(i,j) + H_to_units * (CAS%m_ice(i,j,k) + (CAS%m_snow(i,j,k) + CAS%m_pond(i,j,k)))
  enddo ; enddo ; enddo

end subroutine cell_mass_from_CAS

!> get_total_enthalpy determines the globally integrated enthalpy of snow and ice
subroutine get_total_enthalpy(IST, G, US, IG, enth_ice, enth_snow, scale)
  type(ice_state_type),    intent(in)    :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(in)    :: G   !< The horizontal grid type
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type
  type(EFP_type),          intent(out)   :: enth_ice !< The globally integrated total ice enthalpy [J].
  type(EFP_type),          intent(out)   :: enth_snow !< The globally integrated total snow enthalpy [J].
  real,          optional, intent(in)    :: scale !< A scaling factor from H to the desired units.

  ! Local variables
  real, dimension(:,:,:,:), &
    pointer :: heat_ice=>NULL() ! Pointer to the enth_ice array from the SIS tracer registry.
                        ! Enth_ice is the enthalpy of the ice in each category and layer [Q ~> J kg-1].
  real, dimension(:,:,:,:), &
    pointer :: heat_snow=>NULL() ! Pointer to the enth_snow array from the SIS tracer registry.
                        ! Enth_snow is the enthalpy of the snow atop the ice in each category [Q ~> J kg-1].
  real, dimension(G%isc:G%iec, G%jsc:G%jec) :: sum_enth_ice, sum_enth_snow
  real :: H_to_units ! A conversion factor from H to the desired output units.
  real :: total      ! The returned total enthalpy [J]
  real :: I_Nk       ! The inverse of the number of internal ice layers [nondim].
  integer :: i, j, k, m, isc, iec, jsc, jec, nLay
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  H_to_units = US%RZ_to_kg_m2*US%L_to_m**2 ; if (present(scale)) H_to_units = scale*US%L_to_m**2

  call get_SIS_tracer_pointer("enth_ice", IST%TrReg, heat_ice, nLay)
  call get_SIS_tracer_pointer("enth_snow", IST%TrReg, heat_snow, nLay)
  sum_enth_ice(:,:) = 0.0 ; sum_enth_snow(:,:) = 0.0

  I_Nk = 1.0 / IG%NkIce
  do m=1,IG%NkIce ; do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_enth_ice(i,j) = sum_enth_ice(i,j) + (G%areaT(i,j) * &
              (((H_to_units*IST%mH_ice(i,j,k))*IST%part_size(i,j,k))*I_Nk)) * heat_ice(i,j,k,m)
  enddo ; enddo ; enddo ; enddo
  do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
    sum_enth_snow(i,j) = sum_enth_snow(i,j) + (G%areaT(i,j) * &
              ((H_to_units*IST%mH_snow(i,j,k))*IST%part_size(i,j,k))) * heat_snow(i,j,k,1)
  enddo ; enddo ; enddo
  !### What about sum_enth_pond?

  total = reproducing_sum(sum_enth_ice, EFP_sum=enth_ice)
  total = reproducing_sum(sum_enth_snow, EFP_sum=enth_snow)

end subroutine get_total_enthalpy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_transport_init initializes the ice transport and sets parameters.
subroutine SIS_transport_init(Time, G, US, param_file, diag, CS, continuity_CSp, cover_trans_CSp)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model time.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid type
  type(unit_scale_type),       intent(in)    :: US  !< A structure with unit conversion factors
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(SIS_transport_CS),      pointer       :: CS   !< The control structure for this module
                                                     !! that is allocated and populated here
  type(SIS_continuity_CS), optional, pointer :: continuity_CSp !< The control structure for the
                                                     !!  SIS continuity module
  type(SIS_continuity_CS), optional, pointer :: cover_trans_CSp !< The control structure for ice cover
                                                     !!  transport by the SIS continuity module
!   This subroutine sets the parameters and registers the diagnostics associated
! with the ice dynamics.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "SIS_transport" ! This module's name.
  character(len=80)  :: scheme   ! A string for identifying an advection scheme.
  logical :: merged_cont
  real, parameter :: missing = -1e34

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_transport_init called with an associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag
  CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "RECATEGORIZE_ICE", CS%readjust_categories, &
                 "If true, readjust the distribution into ice thickness "//&
                 "categories after advection.", default=.true.)
  call get_param(param_file, mdl, "RHO_ICE", CS%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "SEA_ICE_ROLL_FACTOR", CS%roll_factor, &
                 "A factor by which the propensity of small amounts of "//&
                 "thick sea-ice to become thinner by rolling is increased "//&
                 "or 0 to disable rolling.  This can be thought of as the "//&
                 "minimum number of ice floes in a grid cell divided by "//&
                 "the horizontal floe aspect ratio.  Sensible values are "//&
                 "0 (no rolling) or larger than 1.", units="Nondim", default=1.0)
  call get_param(param_file, mdl, "ICE_COVER_DISCARD", CS%ice_cover_discard, &
                 "A tiny fractional ice coverage which if positive causes the mass "//&
                 "in categories with less than this coverage to be discarded.", &
                 units="nondim", default=-1.0)


  call get_param(param_file, mdl, "CHECK_ICE_TRANSPORT_CONSERVATION", CS%check_conservation, &
                 "If true, use add multiple diagnostics of ice and snow "//&
                 "mass conservation in the sea-ice transport code.  This "//&
                 "is expensive and should be used sparingly.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DO_RIDGING", CS%do_ridging, &
                 "If true, apply a ridging scheme to the convergent ice. "//&
                 "Otherwise, ice is compressed proportionately if the "//&
                 "concentration exceeds 1.  The original SIS2 implementation "//&
                 "is based on work by Torge Martin.", default=.false.)
  call get_param(param_file, mdl, "SIS_THICKNESS_ADVECTION_SCHEME", scheme, &
          desc="The horizontal transport scheme for thickness:\n"//&
          "  UPWIND_2D - Non-directionally split upwind\n"//&
          "  PCM    - Directionally split piecewise constant\n"//&
          "  PLM    - Piecewise Linear Method\n"//&
          "  PPM:H3 - Piecewise Parabolic Method (Huyhn 3rd order)", &
          default='UPWIND_2D')
  call get_param(param_file, mdl, "ICE_BOUNDS_CHECK", CS%bounds_check, &
                 "If true, periodically check the values of ice and snow "//&
                 "temperatures and thicknesses to ensure that they are "//&
                 "sensible, and issue warnings if they are not.  This "//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)
  call get_param(param_file, mdl, "MERGED_CONTINUITY", merged_cont, &
               "If true, update the continuity equations for the ice, snow, "//&
               "and melt pond water together summed across categories, with "//&
               "proportionate fluxes for each part. Otherwise the media are "//&
               "updated separately.", default=.false., do_not_log=.true.)
  call get_param(param_file, mdl, "INCONSISTENT_COVER_BUG", CS%inconsistent_cover_bug, &
                 "If true, omit a recalculation of the fractional ice-free "//&
                 "areal coverage after the adjustment of the ice categories.", &
                 default=.false., do_not_log=merged_cont)
  if (merged_cont) CS%inconsistent_cover_bug = .false.
  call obsolete_logical(param_file, "ADVECT_TSURF", warning_val=.false.)
  call obsolete_real(param_file, "ICE_CHANNEL_VISCOSITY", warning_val=0.0)
  call obsolete_real(param_file, "ICE_CHANNEL_CFL_LIMIT", warning_val=0.25)
  call obsolete_logical(param_file, "USE_SIS_CONTINUITY", .true.)
  call obsolete_logical(param_file, "USE_SIS_THICKNESS_ADVECTION", .true.)

  call SIS_continuity_init(Time, G, US, param_file, diag, CS%continuity_CSp, &
                           CS_cvr=cover_trans_CSp)
  call SIS_tracer_advect_init(Time, G, param_file, diag, CS%SIS_tr_adv_CSp)

  if (present(continuity_CSp)) continuity_CSp => CS%continuity_CSp

  call SIS_tracer_advect_init(Time, G, param_file, diag, CS%SIS_thick_adv_CSp, scheme=scheme)

  CS%id_ix_trans = register_diag_field('ice_model', 'IX_TRANS', diag%axesCu1, Time, &
               'x-direction ice transport', 'kg/s', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2, &
               missing_value=missing, interp_method='none')
  CS%id_iy_trans = register_diag_field('ice_model', 'IY_TRANS', diag%axesCv1, Time, &
               'y-direction ice transport', 'kg/s', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2, &
               missing_value=missing, interp_method='none')
  CS%id_xprt = register_diag_field('ice_model', 'XPRT', diag%axesT1, Time, &
               'frozen water transport convergence', 'kg/(m^2*yr)', conversion=US%RZ_to_kg_m2, &
               missing_value=missing)
  CS%id_rdgr = register_diag_field('ice_model', 'RDG_RATE', diag%axesT1, Time, &
               'ice ridging rate', '1/sec', conversion=US%s_to_T, missing_value=missing)
!### THESE DIAGNOSTICS DO NOT EXIST YET.
!  CS%id_rdgo = register_diag_field('ice_model', 'RDG_OPEN', diag%axesT1, Time, &
!               'rate of opening due to ridging', '1/s', conversion=US%s_to_T, missing_value=missing)
!  CS%id_rdgv = register_diag_field('ice_model', 'RDG_VOSH', diag%axesT1, Time, &
!               'volume shifted from level to ridged ice', 'm^3/s', conversion=US%RZ_T_to_kg_m2s*US%L_to_m**2, &
!                missing_value=missing)

end subroutine SIS_transport_init

!> Allocate a cell_average_state_type and its elements
subroutine alloc_cell_average_state_type(CAS, HI, IG, CS)
  type(cell_average_state_type), pointer    :: CAS !< A structure with ocean-cell averaged masses
                                                   !! that is being allocated here.
  type(hor_index_type),          intent(in) :: HI  !< The horizontal index type describing the domain
  type(ice_grid_type),           intent(in) :: IG  !< The sea-ice specific grid type
  type(SIS_transport_CS), optional, pointer :: CS  !< A pointer to the control structure for this module

  integer :: isd, ied, jsd, jed, nCat
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nCat = IG%CatIce

  if (.not.associated(CAS)) allocate(CAS)
  call safe_alloc(CAS%m_ice, isd, ied, jsd, jed, ncat)
  call safe_alloc(CAS%m_snow, isd, ied, jsd, jed, ncat)
  call safe_alloc(CAS%m_pond, isd, ied, jsd, jed, ncat)
  call safe_alloc(CAS%mH_ice, isd, ied, jsd, jed, ncat)

  if (present(CS)) then
    if (CS%id_xprt>0) &
      call safe_alloc(CAS%mass0, isd, ied, jsd, jed)
    if (CS%id_ix_trans>0) &
      call safe_alloc(CAS%uh_sum, HI%IsdB, HI%IedB, jsd, jed)
    if (CS%id_iy_trans>0) &
      call safe_alloc(CAS%vh_sum, isd, ied, HI%JsdB, HI%JedB)
  endif
end subroutine alloc_cell_average_state_type

!> Allocate a cell_average_state_type and its elements
subroutine dealloc_cell_average_state_type(CAS)
  type(cell_average_state_type), pointer    :: CAS !< A structure with ocean-cell averaged masses
                                                   !! that is being allocated here.
  if (.not.associated(CAS)) return
  deallocate(CAS%m_ice, CAS%m_snow, CAS%m_pond, CAS%mH_ice)
  if (allocated(CAS%mass0)) deallocate(CAS%mass0)
  if (allocated(CAS%uh_sum)) deallocate(CAS%uh_sum)
  if (allocated(CAS%vh_sum)) deallocate(CAS%vh_sum)
  deallocate(CAS)

end subroutine dealloc_cell_average_state_type

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_transport_end deallocates the memory associated with this module.
subroutine SIS_transport_end(CS)
  type(SIS_transport_CS), pointer :: CS  !< The control structure for this module that
                                         !! is deallocated here

  call SIS_continuity_end(CS%continuity_CSp)
  call SIS_tracer_advect_end(CS%SIS_tr_adv_CSp)

  deallocate(CS)
end subroutine SIS_transport_end

end module SIS_transport
