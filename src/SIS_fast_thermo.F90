!> Handles the rapid thermodynamic interactions between the ice and the atmosphere, including
!! heating and the accumulation of fluxes, but not changes to the ice or snow mass.
module SIS_fast_thermo

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   SIS2 is a SEA ICE MODEL for coupling through the GFDL exchange grid. SIS2  !
! is a revision of the original SIS with have extended capabilities, including !
! the option of using a B-grid or C-grid spatial discretization.  The SIS2     !
! software has been extensively reformulated from SIS for greater consistency  !
! with the Modular Ocean Model, version 6 (MOM6), and to permit might tighter  !
! dynamical coupling between the ocean and sea-ice.                            !
!   This module handles the rapid thermodynamic interactions between the ice   !
! and the atmosphere, including heating and the accumulation of fluxes, but    !
! not changes to the ice or snow mass.                                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use ice_boundary_types, only : atmos_ice_boundary_type ! , land_ice_boundary_type
use ice_grid,           only : ice_grid_type

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_time_manager,  only : time_type, time_type_to_real, operator(+), operator(-)
use MOM_time_manager,  only : operator(>), operator(*), operator(/), operator(/=)
use MOM_unit_scaling,  only : unit_scale_type

use SIS_debugging,     only : hchksum
use SIS_diag_mediator, only : SIS_diag_ctrl
use SIS_framework,     only : coupler_3d_bc_type, coupler_type_spawn
use SIS_framework,     only : coupler_type_increment_data, coupler_type_rescale_data
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_optics,        only : ice_optics_SIS2, bright_ice_temp, SIS_optics_CS
use SIS_optics,        only : VIS_DIR, VIS_DIF, NIR_DIR, NIR_DIF
use SIS_types,         only : ice_state_type, IST_chksum, IST_bounds_check, ice_rad_type
use SIS_types,         only : fast_ice_avg_type, simple_OSS_type, total_sfc_flux_type, FIA_chksum
use SIS2_ice_thm,      only : SIS2_ice_thm_CS, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,      only : ice_temp_SIS2, latent_sublimation
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs, enth_from_TS, Temp_from_En_S

implicit none ; private

public :: do_update_ice_model_fast, SIS_fast_thermo_init, SIS_fast_thermo_end
public :: accumulate_deposition_fluxes, convert_frost_to_snow
public :: fast_thermo_CS, avg_top_quantities, total_top_quantities, infill_array
public :: redo_update_ice_model_fast, find_excess_fluxes

!> The control structure for the SIS fast thermodynamics module
type fast_thermo_CS ; private
  ! These two arrays are used with column_check when evaluating the enthalpy
  ! conservation with the fast thermodynamics code.
  real, pointer, dimension(:,:,:) :: enth_prev => NULL() !< The previous enthalpy [Q R Z ~> J m-2], used with
                                     !! column_check when evaluating the enthalpy conservation
                                     !! with the fast thermodynamics code
  real, pointer, dimension(:,:,:) :: heat_in => NULL() !< The heat input [Q R Z ~> J m-2],  used with
                                     !! column_check when evaluating the enthalpy conservation
                                     !! with the fast thermodynamics code

  logical :: debug_fast   !< If true, write verbose checksums of code that is
                          !! executed on fast ice PEs for debugging purposes.
  logical :: debug_slow   !< If true, write verbose checksums of code that is
                          !! executed on slow ice PEs for debugging purposes.
  logical :: column_check !< If true, enable the heat check column by column.
  real    :: imb_tol      !< The tolerance for imbalances to be flagged by
                          !! column_check [nondim].
  logical :: bounds_check !< If true, check for sensible values of thicknesses
                          !! temperatures, fluxes, etc.

  integer :: n_fast = 0   !< The number of times update_ice_model_fast has been called.
  logical :: Reorder_0C_heatflux !< If true, rearrange the calculation of the heat fluxes projected
                          !! back to 0C to work on each contribution separately, so that they can
                          !! be identically replicated if there is a single fast timestep per
                          !!  coupled timestep and REDO_FAST_ICE_UPDATE=True
  integer :: max_tskin_itt !< The maximum number of iterations of the skin temperature and
                          !! optical properties during redo_update_ice_model_fast.

  !> A pointer to the control structures for subsidiary modules.
  type(SIS2_ice_thm_CS), pointer  :: ice_thm_CSp => NULL()
end type fast_thermo_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> sum_top_quantities does a running sum of fluxes for later use by the slow ice
!!   physics and the ocean.  Nothing here will be exposed to other modules until
!!   after it has passed through avg_top_quantities.
subroutine sum_top_quantities (FIA, ABT, flux_u, flux_v, flux_sh, evap, &
       flux_sw, flux_lw, lprec, fprec, flux_lh, t_skin, SST, &
       sh_T0, evap_T0, lw_T0, dshdt, devapdt, dlwdt, G, US, IG)
  type(fast_ice_avg_type),       intent(inout) :: FIA !< A type containing averages of fields
                                                      !! (mostly fluxes) over the fast updates
  type(atmos_ice_boundary_type), intent(in)    :: ABT !< A type containing atmospheric boundary
                                                      !! forcing fields that are used to drive the ice
  type(SIS_hor_grid_type),       intent(in)    :: G   !< The horizontal grid type
  type(ice_grid_type),           intent(in)    :: IG  !< The sea-ice specific grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: flux_u   !< The grid-wise quasi-zonal wind stress on the ice [R L Z T-2 ~> Pa].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: flux_v   !< The grid-wise quasi-meridional wind stress on the ice [R L Z T-2 ~> Pa].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: flux_sh  !< The upward sensible heat flux from the top of the ice into
                           !! the atmosphere [Q R Z T-1 ~> W m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: evap     !< The upward flux of water due to sublimation or evaporation
                           !! from the top of the ice to the atmosphere [R Z T-1 ~> kg m-2 s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: flux_lw  !< The net longwave heat flux from the atmosphere into the
                           !! ice or ocean [Q R Z T-1 ~> W m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: lprec    !< The liquid precipitation onto the ice [R Z T-1 ~> kg m-2 s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: fprec    !< The frozen precipitation onto the ice [R Z T-1 ~> kg m-2 s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: flux_lh  !< The upward latent heat flux associated with sublimation or
                           !! evaporation [Q R Z T-1 ~> W m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: sh_T0    !< The upward sensible heat flux from the top of the ice into
                           !! the atmosphere when the skin temperature is 0 degC [Q R Z T-1 ~> W m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: evap_T0  !< The sublimation rate when the skin temperature is 0 degC,
                           !! [R Z T-1 ~> kg m-2 s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: lw_T0    !< The downward longwave heat flux from the atmosphere into the
                           !! ice or ocean when the skin temperature is 0 degC [Q R Z T-1 ~> W m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: dshdt    !< The derivative of the upward sensible heat flux from the
                           !! the top of the ice into the atmosphere with ice skin
                           !! temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: devapdt  !< The derivative of the sublimation rate with the surface
                           !! temperature [R Z T-1 degC-1 ~> kg m-2 s-1 degC-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
    intent(in) :: dlwdt    !< The derivative of the longwave heat flux from the atmosphere
                           !! into the ice or ocean with ice skin temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,IG%CatIce), &
    intent(in) :: t_skin   !< The sea ice surface skin temperature [degC].
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
    intent(in) :: SST      !< The sea surface temperature [degC].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce,size(FIA%flux_sw_top,4)), &
    intent(in) :: flux_sw  !< The downward shortwave heat fluxes [Q R Z T-1 ~> W m-2]. The 4th
                           !! dimension is a combination of angular orientation & frequency.
  type(unit_scale_type),         intent(in)    :: US        !< A structure with unit conversion factors

  real :: t_sfc
  integer :: i, j, k, m, n, b, nb, i2, j2, k2, isc, iec, jsc, jec, i_off, j_off, ncat
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  nb = size(FIA%flux_sw_top,4)

  i_off = LBOUND(ABT%t_flux,1) - G%isc
  j_off = LBOUND(ABT%t_flux,2) - G%jsc

  call coupler_type_spawn(ABT%fluxes, FIA%tr_flux, (/isd, isc, iec, ied/), &
                          (/jsd, jsc, jec, jed/), (/0, ncat/), as_needed=.true.)

  if (FIA%avg_count == 0) then
    ! zero_top_quantities - zero fluxes to begin summing in ice fast physics.
    FIA%flux_u_top(:,:,:) = 0.0 ; FIA%flux_v_top(:,:,:) = 0.0
    FIA%flux_sh_top(:,:,:) = 0.0 ; FIA%evap_top(:,:,:) = 0.0
    FIA%flux_lw_top(:,:,:) = 0.0 ; FIA%flux_lh_top(:,:,:) = 0.0
    FIA%flux_sw_top(:,:,:,:) = 0.0
    FIA%lprec_top(:,:,:) = 0.0 ; FIA%fprec_top(:,:,:) = 0.0
    FIA%flux_sw_dn(:,:,:) = 0.0 ; FIA%Tskin_avg(:,:) = 0.0

    if (allocated(FIA%flux_sh0)) then
      FIA%dshdt(:,:,:) = 0.0 ; FIA%devapdt(:,:,:) = 0.0 ; FIA%dlwdt(:,:,:) = 0.0
      FIA%flux_sh0(:,:,:) = 0.0 ; FIA%evap0(:,:,:) = 0.0
      FIA%flux_lw0(:,:,:) = 0.0 ; FIA%Tskin_cat(:,:,:) = 0.0
    endif

    call coupler_type_rescale_data(FIA%tr_flux, 0.0)
  endif

  !$OMP parallel do default(shared)
  do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
    FIA%flux_u_top(i,j,k)  = FIA%flux_u_top(i,j,k)  + flux_u(i,j,k)
    FIA%flux_v_top(i,j,k)  = FIA%flux_v_top(i,j,k)  + flux_v(i,j,k)
    FIA%flux_sh_top(i,j,k)  = FIA%flux_sh_top(i,j,k) + flux_sh(i,j,k)
    FIA%evap_top(i,j,k)  = FIA%evap_top(i,j,k)  + evap(i,j,k)
    do b=1,nb ; FIA%flux_sw_top(i,j,k,b) = FIA%flux_sw_top(i,j,k,b) + flux_sw(i,j,k,b) ; enddo
    FIA%flux_lw_top(i,j,k) = FIA%flux_lw_top(i,j,k) + flux_lw(i,j,k)
    FIA%lprec_top(i,j,k)   = FIA%lprec_top(i,j,k)   + lprec(i,j,k)
    FIA%fprec_top(i,j,k)   = FIA%fprec_top(i,j,k)   + fprec(i,j,k)
    FIA%flux_lh_top(i,j,k) = FIA%flux_lh_top(i,j,k) + flux_lh(i,j,k)
  enddo ; enddo ; enddo
  ! FIA%flux_sw_dn is accumulated where the fast radiation diagnostics are output
  ! because it depends on arrays that are stored in the public ice_data_type.

  !Do not handle air_sea_deposition fluxes here, they need to be handled after atmos_down
  call coupler_type_increment_data( ABT%fluxes, FIA%tr_flux, exclude_flux_type='air_sea_deposition')

  if (allocated(FIA%flux_sh0)) then
    !$OMP parallel do default(shared) private(t_sfc)
    do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
      FIA%dshdt(i,j,k) = FIA%dshdt(i,j,k) + dshdt(i,j,k)
      FIA%devapdt(i,j,k) = FIA%devapdt(i,j,k) + devapdt(i,j,k)
      FIA%dlwdt(i,j,k) = FIA%dlwdt(i,j,k) + dlwdt(i,j,k)
      FIA%flux_sh0(i,j,k) = FIA%flux_sh0(i,j,k) + sh_T0(i,j,k)
      FIA%evap0(i,j,k) = FIA%evap0(i,j,k) + evap_T0(i,j,k)
      FIA%flux_lw0(i,j,k) = FIA%flux_lw0(i,j,k) + lw_T0(i,j,k)
      t_sfc = SST(i,j) ; if (k>0) t_sfc = t_skin(i,j,k)
      FIA%Tskin_cat(i,j,k) = FIA%Tskin_cat(i,j,k) + t_sfc
    enddo ; enddo ; enddo
  endif

  FIA%avg_count = FIA%avg_count + 1

end subroutine sum_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> avg_top_quantities determines time average fluxes for later use by the
!!    slow ice physics and by the ocean.
subroutine avg_top_quantities(FIA, Rad, IST, G, US, IG)
  type(fast_ice_avg_type), intent(inout) :: FIA !< A type containing averages of fields
                                                !! (mostly fluxes) over the fast updates
  type(ice_rad_type),      intent(in)    :: Rad !< A structure with fields related to the absorption,
                                                !! reflection and transmission of shortwave radiation.
  type(ice_state_type),    intent(in)    :: IST !< A type describing the state of the sea ice
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),   intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),     intent(in)    :: IG  !< The sea-ice specific grid type

  real    :: u, v, divid, sign
  real    :: I_avc    ! The inverse of the number of contributions.
  real    :: I_wts    ! 1.0 / ice_cover or 0 if ice_cover is 0 [nondim].
  integer :: i, j, k, m, n, b, nb, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  nb = size(FIA%flux_sw_top,4)

  !
  ! compute average fluxes
  !
  if (FIA%avg_count == 0) call SIS_error(FATAL,'avg_top_quantities: '//&
       'no ocean model fluxes have been averaged')

  ! Rotate the stress from lat/lon to ocean coordinates and possibly change the
  ! sign to positive for downward fluxes of positive momentum.
  sign = 1.0 ; if (FIA%atmos_winds) sign = -1.0
  I_avc = 1.0/real(FIA%avg_count)

  !$OMP parallel do default(shared) private(u,v)
  do j=jsc,jec
    do k=0,ncat ; do i=isc,iec
      u = FIA%flux_u_top(i,j,k) * (sign*I_avc)
      v = FIA%flux_v_top(i,j,k) * (sign*I_avc)
      FIA%flux_u_top(i,j,k) = u*G%cos_rot(i,j)-v*G%sin_rot(i,j) ! rotate stress from lat/lon
      FIA%flux_v_top(i,j,k) = v*G%cos_rot(i,j)+u*G%sin_rot(i,j) ! to ocean coordinates
      FIA%flux_sh_top(i,j,k)  = FIA%flux_sh_top(i,j,k)  * I_avc
      FIA%evap_top(i,j,k)  = FIA%evap_top(i,j,k)  * I_avc
      do b=1,nb ; FIA%flux_sw_top(i,j,k,b) = FIA%flux_sw_top(i,j,k,b) * I_avc ; enddo
      FIA%flux_lw_top(i,j,k) = FIA%flux_lw_top(i,j,k) * I_avc
      FIA%fprec_top(i,j,k)   = FIA%fprec_top(i,j,k)   * I_avc
      FIA%lprec_top(i,j,k)   = FIA%lprec_top(i,j,k)   * I_avc
      FIA%flux_lh_top(i,j,k) = FIA%flux_lh_top(i,j,k) * I_avc

      ! Copy radiation fields from the fast to the slow states.
      if (k>0) FIA%sw_abs_ocn(i,j,k) = Rad%sw_abs_ocn(i,j,k)
      ! Convert frost forming atop sea ice into frozen precip.
!      if ((k>0) .and. (FIA%evap_top(i,j,k) < 0.0)) then
!        FIA%fprec_top(i,j,k) = FIA%fprec_top(i,j,k) - FIA%evap_top(i,j,k)
!        FIA%evap_top(i,j,k) = 0.0
!      endif
    enddo ; enddo

    do b=1,nb ; do i=isc,iec
      FIA%flux_sw_dn(i,j,b) = FIA%flux_sw_dn(i,j,b)*I_avc
    enddo ; enddo
    do i=isc,iec
      FIA%Tskin_avg(i,j) = FIA%Tskin_avg(i,j) * I_avc
    enddo
  enddo
  call coupler_type_rescale_data(FIA%tr_flux, I_avc)

  ! Determine the fractional ice coverage and the wind stresses averaged
  ! across all the ice thickness categories on an A-grid.
  FIA%WindStr_x(:,:) = 0.0 ; FIA%WindStr_y(:,:) = 0.0 ; FIA%ice_cover(:,:) = 0.0

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,FIA,Rad,IST) &
!$OMP                           private(I_wts)
  do j=jsc,jec
    do k=1,ncat ; do i=isc,iec
      FIA%WindStr_x(i,j) = FIA%WindStr_x(i,j) + IST%part_size(i,j,k) * FIA%flux_u_top(i,j,k)
      FIA%WindStr_y(i,j) = FIA%WindStr_y(i,j) + IST%part_size(i,j,k) * FIA%flux_v_top(i,j,k)
      FIA%ice_cover(i,j) = FIA%ice_cover(i,j) + IST%part_size(i,j,k)
    enddo ; enddo
    do i=isc,iec
      FIA%WindStr_ocn_x(i,j) = FIA%flux_u_top(i,j,0)
      FIA%WindStr_ocn_y(i,j) = FIA%flux_v_top(i,j,0)
      if (FIA%ice_cover(i,j) > 0.0) then
        I_wts = 1.0 / FIA%ice_cover(i,j)
        FIA%WindStr_x(i,j) = FIA%WindStr_x(i,j) * I_wts
        FIA%WindStr_y(i,j) = FIA%WindStr_y(i,j) * I_wts
        if (FIA%ice_cover(i,j) > 1.0) FIA%ice_cover(i,j) = 1.0

        ! The max with 0 in the following line is here for safety; the only known
        ! instance where it has been required is when reading a SIS-1-derived
        ! restart file with tiny negative concentrations. SIS2 should not need it.
        FIA%ice_free(i,j) = max(IST%part_size(i,j,0), 0.0)

    !    Rescale to add up to 1?
    !    I_wts = 1.0 / (FIA%ice_free(i,j) + FIA%ice_cover(i,j))
    !    FIA%ice_free(i,j) = FIA%ice_free(i,j) * I_wts
    !    FIA%ice_cover(i,j) = FIA%ice_cover(i,j) * I_wts
      else
        FIA%ice_free(i,j) = 1.0 ; FIA%ice_cover(i,j) = 0.0
        FIA%WindStr_x(i,j) = FIA%flux_u_top(i,j,0)
        FIA%WindStr_y(i,j) = FIA%flux_v_top(i,j,0)
      endif
    enddo
  enddo

  if (allocated(FIA%flux_sh0)) then
    !$OMP parallel do default(shared)
    do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
      FIA%dshdt(i,j,k) = FIA%dshdt(i,j,k) * I_avc
      FIA%devapdt(i,j,k) = FIA%devapdt(i,j,k) * I_avc
      FIA%dlwdt(i,j,k) = FIA%dlwdt(i,j,k) * I_avc
      FIA%flux_sh0(i,j,k) = FIA%flux_sh0(i,j,k) * I_avc
      FIA%evap0(i,j,k) = FIA%evap0(i,j,k) * I_avc
      FIA%flux_lw0(i,j,k) = FIA%flux_lw0(i,j,k) * I_avc
      FIA%Tskin_cat(i,j,k) = FIA%Tskin_cat(i,j,k) * I_avc
    enddo ; enddo ; enddo

    ! Fill in the information to reconstruct the fluxes for any area-less categories.
    ! The open-ocean category must always be calculated for this to work properly.
    call infill_array(IST, FIA%dshdt(:,:,0), FIA%dshdt(:,:,1:), G, IG)
    call infill_array(IST, FIA%devapdt(:,:,0), FIA%devapdt(:,:,1:), G, IG)
    call infill_array(IST, FIA%dlwdt(:,:,0), FIA%dlwdt(:,:,1:), G, IG)
    call infill_array(IST, FIA%flux_sh0(:,:,0), FIA%flux_sh0(:,:,1:), G, IG)
    call infill_array(IST, FIA%evap0(:,:,0), FIA%evap0(:,:,1:), G, IG)
    call infill_array(IST, FIA%flux_lw0(:,:,0), FIA%flux_lw0(:,:,1:), G, IG)
    call infill_array(IST, FIA%Tskin_cat(:,:,0), FIA%Tskin_cat(:,:,1:), G, IG)
  endif

  ! set count to zero and fluxes will be zeroed before the next sum
  FIA%avg_count = 0

end subroutine avg_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> total_top_quantities determines the sum across partitions of various fluxes
!! for later use on a potentially different ice state on the slow side.
subroutine total_top_quantities(FIA, TSF, part_size, G, US, IG)
  type(fast_ice_avg_type),   intent(in)    :: FIA !< A type containing averages of fields
                                                  !! (mostly fluxes) over the fast updates
  type(total_sfc_flux_type), intent(inout) :: TSF !< A type with fluxes that are averaged across
                                                  !! the fast updates and integrated across thickness
                                                  !! categories from the fast ice update
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),       intent(in)    :: IG  !< The sea-ice specific grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
                             intent(in)    :: part_size !< The fractional area coverage of the ice
                                                  !! thickness categories [nondim], 0-1

  integer :: i, j, k, m, n, b, nb, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  nb = size(FIA%flux_sw_top,4)

  call coupler_type_spawn(FIA%tr_flux, TSF%tr_flux, (/isd, isc, iec, ied/), &
                          (/jsd, jsc, jec, jed/), as_needed=.true. )

  TSF%flux_u(:,:) = 0.0 ; TSF%flux_v(:,:) = 0.0
  TSF%flux_sh(:,:) = 0.0 ; TSF%flux_lw(:,:) = 0.0 ; TSF%flux_lh(:,:) = 0.0
  TSF%flux_sw(:,:,:) = 0.0
  TSF%evap(:,:) = 0.0 ; TSF%fprec(:,:) = 0.0 ; TSF%lprec(:,:) = 0.0
  call coupler_type_rescale_data(TSF%tr_flux, 0.0)

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    TSF%flux_u(i,j) = TSF%flux_u(i,j) + part_size(i,j,k) * FIA%flux_u_top(i,j,k)
    TSF%flux_v(i,j) = TSF%flux_v(i,j) + part_size(i,j,k) * FIA%flux_v_top(i,j,k)
    TSF%flux_sh(i,j) = TSF%flux_sh(i,j) + part_size(i,j,k) * FIA%flux_sh_top(i,j,k)
    TSF%flux_lw(i,j) = TSF%flux_lw(i,j) + part_size(i,j,k) * FIA%flux_lw_top(i,j,k)
    TSF%flux_lh(i,j) = TSF%flux_lh(i,j) + part_size(i,j,k) * FIA%flux_lh_top(i,j,k)
    do b=1,nb
      TSF%flux_sw(i,j,b) = TSF%flux_sw(i,j,b) + part_size(i,j,k) * FIA%flux_sw_top(i,j,k,b)
    enddo
    TSF%evap(i,j) = TSF%evap(i,j) + part_size(i,j,k) * FIA%evap_top(i,j,k)
    TSF%fprec(i,j) = TSF%fprec(i,j) + part_size(i,j,k) * FIA%fprec_top(i,j,k)
    TSF%lprec(i,j) = TSF%lprec(i,j) + part_size(i,j,k) * FIA%lprec_top(i,j,k)
  enddo ; enddo ; enddo
  call coupler_type_increment_data(FIA%tr_flux, part_size, TSF%tr_flux)

  !   If the sum of part_size across all the ice and ocean categories is not
  ! exactly 1, rescaling might be advisable, but for now it is assumed that
  ! part_size is properly scaled.

end subroutine total_top_quantities

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> find_excess_fluxes determines the difference between the sum across
!! partitions of various fluxes and the sum previously found by total_top_quantities.
subroutine find_excess_fluxes(FIA, TSF, XSF, part_size, G, US, IG)
  type(fast_ice_avg_type),   intent(in)    :: FIA !< A type containing averages of fields
                                                  !! (mostly fluxes) over the fast updates
  type(total_sfc_flux_type), intent(in)    :: TSF !< A type with fluxes that are averaged across
                                                  !! the fast updates and integrated across thickness
                                                  !! categories from the fast ice update
  type(total_sfc_flux_type), intent(inout) :: XSF !< A structure of the excess fluxes between the
                                                  !! atmosphere and the ice or ocean relative to
                                                  !! those stored in TSF
  type(SIS_hor_grid_type),   intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),     intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),       intent(in)    :: IG  !< The sea-ice specific grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce), &
                             intent(in)    :: part_size !< The fractional area coverage of the ice
                                                  !! thickness categories [nondim], 0-1
  integer :: i, j, k, m, n, b, nb, isc, iec, jsc, jec, ncat
  integer :: isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  nb = size(FIA%flux_sw_top,4)

  ! Check whether this is the first call when the number of tracer fluxes are
  ! known, and the XSF tracer flux arrays need to be allocated now.  This call
  ! only does anything the first or second time that total_top_quantities is called.
  call coupler_type_spawn(FIA%tr_flux, XSF%tr_flux, (/isd, isc, iec, ied/), &
                          (/jsd, jsc, jec, jed/), as_needed=.true. )

  ! Note that XSF%flux_u and XSF%flux_v are not necessary as the changing ice cover is
  ! already being taken into account.

  ! This could be done more succinctly by initializing XSF = -TSF, but the
  ! answers would then be more accutely impacted by roundoff.
  XSF%flux_sh(:,:) = 0.0 ; XSF%flux_lw(:,:) = 0.0 ; XSF%flux_lh(:,:) = 0.0
  XSF%flux_sw(:,:,:) = 0.0
  XSF%evap(:,:) = 0.0 ; XSF%fprec(:,:) = 0.0 ; XSF%lprec(:,:) = 0.0
  call coupler_type_rescale_data(XSF%tr_flux, 0.0)

  do k=0,ncat ; do j=jsc,jec ; do i=isc,iec
    XSF%flux_sh(i,j) = XSF%flux_sh(i,j) + part_size(i,j,k) * FIA%flux_sh_top(i,j,k)
    XSF%flux_lw(i,j) = XSF%flux_lw(i,j) + part_size(i,j,k) * FIA%flux_lw_top(i,j,k)
    XSF%flux_lh(i,j) = XSF%flux_lh(i,j) + part_size(i,j,k) * FIA%flux_lh_top(i,j,k)
    do b=1,nb
      XSF%flux_sw(i,j,b) = XSF%flux_sw(i,j,b) + part_size(i,j,k) * FIA%flux_sw_top(i,j,k,b)
    enddo
    XSF%evap(i,j) = XSF%evap(i,j) + part_size(i,j,k) * FIA%evap_top(i,j,k)
    XSF%fprec(i,j) = XSF%fprec(i,j) + part_size(i,j,k) * FIA%fprec_top(i,j,k)
    XSF%lprec(i,j) = XSF%lprec(i,j) + part_size(i,j,k) * FIA%lprec_top(i,j,k)
  enddo ; enddo ; enddo
  call coupler_type_increment_data(FIA%tr_flux, part_size, XSF%tr_flux)

  do j=jsc,jec ; do i=isc,iec
    XSF%flux_sh(i,j) = XSF%flux_sh(i,j) - TSF%flux_sh(i,j)
    XSF%flux_lw(i,j) = XSF%flux_lw(i,j) - TSF%flux_lw(i,j)
    XSF%flux_lh(i,j) = XSF%flux_lh(i,j) - TSF%flux_lh(i,j)
    do b=1,nb ; XSF%flux_sw(i,j,b) = XSF%flux_sw(i,j,b) - TSF%flux_sw(i,j,b) ; enddo
    XSF%evap(i,j) = XSF%evap(i,j)   - TSF%evap(i,j)
    XSF%fprec(i,j) = XSF%fprec(i,j) - TSF%fprec(i,j)
    XSF%lprec(i,j) = XSF%lprec(i,j) - TSF%lprec(i,j)
  enddo ; enddo
  call coupler_type_increment_data(TSF%tr_flux, XSF%tr_flux, scale_factor=-1.0)

end subroutine find_excess_fluxes


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> infill_array fills in an array with actual, interpolated or plausible values
!! from other ice categories.
subroutine infill_array(IST, ta_ocn, ta, G, IG)
  type(ice_state_type),    intent(in)    :: IST  !< The ice state type whose values of part_size and
                                                 !! mH_ice are being used to control the infilling.
  type(SIS_hor_grid_type), intent(in)    :: G    !< The sea-ice lateral grid type.
  type(ice_grid_type),     intent(in)    :: IG   !< The sea-ice grid type.
  real, dimension(G%isd:G%ied, G%jsd:G%jed), &
                           intent(in)    :: ta_ocn !< The value of the array for ocean partitions.
  real, dimension(G%isd:G%ied, G%jsd:G%jed, IG%CatIce), &
                           intent(inout) :: ta   !< The array that is being infilled.

  integer, dimension(G%isd:G%ied) :: k_thick, k_thin
  real, dimension(G%isd:G%ied) :: &
    ta_thick, ta_thin, mH_thin
  real :: wt2, mH_cat_tgt
  logical, dimension(G%isd:G%ied) :: infill
  logical :: any_ice
  integer :: i, j, k, k2, k3, isc, iec, jsc, jec, ncat, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; NkIce = IG%NkIce

  do j=jsc,jec
    ! Determine which categories exist.
    do i=isc,iec ; k_thick(i) = -1 ; k_thin(i) = ncat+1 ; mH_thin(i) = 0.0 ; enddo
    any_ice = .false.
    do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k) > 0.0) then
      k_thick(i) = k
      ta_thick(i) = ta(i,j,k)
      any_ice = .true.
    endif ; enddo ; enddo

    do i=isc,iec
      if (k_thick(i) < 1) ta_thick(i) = G%mask2dT(i,j) * ta_ocn(i,j)
    enddo

    if (.not.any_ice) then
      ! This entire row is ice free.
      do k=1,ncat ; do i=isc,iec
        ta(i,j,k) = G%mask2dT(i,j)*ta_ocn(i,j)
      enddo ; enddo
    else
      do k=ncat,1,-1 ; do i=isc,iec ; if (IST%part_size(i,j,k) > 0.0) then
        k_thin(i) = k
        ta_thin(i) = ta(i,j,k)
        mH_thin(i) = IST%mH_ice(i,j,k)
      endif ; enddo ; enddo

      ! The logic for the k=1 (thinest) category is particularly simple.
      do i=isc,iec ; if (IST%part_size(i,j,1) <= 0.0) then
        ta(i,j,1) = G%mask2dT(i,j)*ta_ocn(i,j)
      endif ; enddo
      do k=2,ncat ; do i=isc,iec
        if (G%mask2dT(i,j) > 0.0) then
          if (k > k_thick(i)) then
            ! This is the most common case, since it applies to ice-free water.
            ta(i,j,k) = ta_thick(i)
          elseif (IST%part_size(i,j,k) > 0.0) then
            ! This is a valid value - no interpolation is needed.
            cycle
          elseif (k < k_thin(i)) then
            ! Linearly interpolate to the ocean's freezing temperature.
            mH_cat_tgt = 0.5*(IG%mH_cat_bound(K) + IG%mH_cat_bound(K+1))
            wt2 = 0.0  ! If in doubt, use ta_thin.  This should not happen.
            if (mH_thin(i) > mH_cat_tgt) wt2 = (mH_thin(i) - mH_cat_tgt) / mH_thin(i)

            ta(i,j,k) = wt2*ta_ocn(i,j) + (1.0-wt2)*ta_thin(i)
          else
            ! Linearly interpolate between bracketing neighboring
            ! categories with valid values.
            do k2=k+1,ncat ; if (IST%part_size(i,j,k2) > 0.0) exit ; enddo
            do k3=k-1,1,-1 ; if (IST%part_size(i,j,k3) > 0.0) exit ; enddo
            mH_cat_tgt = 0.5*(IG%mH_cat_bound(K) + IG%mH_cat_bound(K+1))
            wt2 = 0.5
            if ((IST%mH_ice(i,j,k3) > mH_cat_tgt) .and. &
                (IST%mH_ice(i,j,k2) < mH_cat_tgt)) &
              wt2 = (IST%mH_ice(i,j,k3) - mH_cat_tgt) / &
                    (IST%mH_ice(i,j,k3) - IST%mH_ice(i,j,k2))

            ta(i,j,k) = wt2*ta(i,j,k2) + (1.0-wt2) * ta(i,j,k3)
          endif
        else ! This is a land point.
          ta(i,j,k) = 0.0
        endif
      enddo ; enddo

    endif ! .not.any_ice
  enddo

end subroutine infill_array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> do_update_ice_model_fast applies the surface heat fluxes, shortwave radiation
!!   diffusion of heat to the sea-ice to implicitly determine a new temperature
!!   profile, subject to the constraint that ice and snow temperatures are never
!!   above freezing.  Melting and freezing occur elsewhere.
subroutine do_update_ice_model_fast(Atmos_boundary, IST, sOSS, Rad, FIA, &
                                    Time_step, CS, G, US, IG )

  type(atmos_ice_boundary_type), intent(in) :: Atmos_boundary !< A type containing atmospheric boundary
                                                        !! forcing fields that are used to drive the ice
  type(ice_state_type),      intent(inout) :: IST       !< A type describing the state of the sea ice
  type(simple_OSS_type),     intent(in)    :: sOSS      !< A type describing the ocean surface state
  type(ice_rad_type),        intent(inout) :: Rad       !< A type containing fields related to shortwave
                                                        !! radiation in the ice, including the skin temperature
  type(fast_ice_avg_type),   intent(inout) :: FIA       !< A type that is used to accumulate averages of
                                                        !! fields (mostly fluxes) over the fast updates
  type(time_type),           intent(in)    :: Time_step !< The amount of time over which to advance the ice
  type(fast_thermo_CS),      pointer       :: CS        !< The control structure for this module
  type(SIS_hor_grid_type),   intent(inout) :: G         !< The horizontal grid type
  type(unit_scale_type),     intent(in)    :: US        !< A structure with unit conversion factors
  type(ice_grid_type),       intent(in)    :: IG        !< The ice vertical grid type


  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce) :: &
    flux_sh, &  ! The upward sensible heat flux from the ice to the atmosphere
                ! at the surface of the ice [Q R Z T-1 ~> W m-2].
    evap, &     ! The upward flux of water due to sublimation or evaporation
                ! from the top of the ice to the atmosphere [R Z T-1 ~> kg m-2 s-1].
    flux_lh, &  ! The upward latent heat flux associated with sublimation or
                ! evaporation [Q R Z T-1 ~> W m-2].
    flux_lw, &  ! The net downward longwave heat flux into the ice [Q R Z T-1 ~> W m-2].
    flux_u, &   ! The grid-aligned quasi-zonal wind stress on the ice [R L Z T-2 ~> Pa].
    flux_v, &   ! The grid-aligned quasi-meridional wind stress on the ice [R L Z T-2 ~> Pa].
    lprec, &    ! The liquid precipitation onto the ice [R Z T-1 ~> kg m-2 s-1].
    fprec, &    ! The frozen precipitation onto the ice [R Z T-1 ~> kg m-2 s-1].
    sh_T0, &    ! The upward sensible heat flux from the top of the ice into
                ! the atmosphere when the skin temperature is 0 degC [Q R Z T-1 ~> W m-2].
    evap_T0, &  ! The sublimation rate  when the skin temperature is 0 degC,
                ! [R Z T-1 ~> kg m-2 s-1].
    lw_T0, &    ! The downward longwave heat flux from the atmosphere into the
                ! ice or ocean when the skin temperature is 0 degC [Q R Z T-1 ~> W m-2].
    dshdt, &    ! The derivative of the upward sensible heat flux with the surface
                ! temperature [Q R Z T-1 degC-1 ~> W m-2 degC-1].
    devapdt, &  ! The derivative of the sublimation rate with the surface
                ! temperature [R Z T-1 degC-1 ~> kg m-2 s-1 degC-1].
    dlwdt       ! The derivative of the downward radiative heat flux with surface
                ! temperature (i.e. d(flux_lw)/d(surf_temp)) [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,0:IG%CatIce,size(FIA%flux_sw_top,4)) :: &
    flux_sw     ! The downward shortwave heat fluxes [Q R Z T-1 ~> W m-2].  The fourth
                ! dimension is a combination of angular orientation and frequency.
  real, dimension(0:IG%NkIce) :: T_col ! The temperature of a column of ice and snow [degC].
  real, dimension(IG%NkIce)   :: S_col ! The thermodynamic salinity of a column of ice [gSalt kg-1].
  real, dimension(0:IG%NkIce) :: enth_col   ! The enthalpy of a column of snow and ice [Q ~> J kg-1].
  real, dimension(0:IG%NkIce) :: SW_abs_col ! The shortwave absorption within a column of snow and
                  ! ice [Q R Z T-1 ~> W m-2].
  real :: dt_fast ! The fast thermodynamic time step [T ~> s].
  real :: Tskin   ! The new skin temperature [degC].
  real :: dTskin  ! The change in the skin temperature [degC].
  real :: latent  ! The latent heat of sublimation of ice or snow [Q ~> J kg-1].
  real :: hf_0    ! The positive upward surface heat flux when T_sfc = 0 degC [Q R Z T-1 ~> W m-2].
  real :: dhf_dt  ! The derivative of the upward surface heat flux with Ts [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: sw_tot  ! sum over all shortwave (dir/dif and vis/nir) components [Q R Z T-1 ~> W m-2].
  real :: snow_wt ! A fractional weighting of snow in the category surface area [nondim].
  real :: LatHtVap       ! The latent heat of vaporization of water at 0C [Q ~> J kg-1].
  integer :: i, j, k, m, i2, j2, k2, isc, iec, jsc, jec, ncat, i_off, j_off, NkIce, b, nb
  character(len=8) :: nstr

  real :: tot_heat_in, enth_here, enth_imb, norm_enth_imb
  real :: SW_absorbed ! Absorbed shortwave heating [Q R Z T-1 ~> W m-2]
  real :: I_Nk     ! The inverse of the number of internal ice layers [nondim].

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_fast_thermo: Module must be initialized before it is used.")

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  i_off = LBOUND(Atmos_boundary%t_flux,1) - G%isc
  j_off = LBOUND(Atmos_boundary%t_flux,2) - G%jsc
  NkIce = IG%NkIce ; I_Nk = 1.0 / NkIce
  nb = size(FIA%flux_sw_top,4)

  CS%n_fast = CS%n_fast + 1

  if (CS%column_check .and. (FIA%avg_count==0)) then
    CS%heat_in(:,:,:) = 0.0
    CS%enth_prev(:,:,:) = 0.0
    do k=1,ncat ; do j=jsc,jec ; do i=isc,iec ; if (IST%mH_ice(i,j,k)>0.0) then
      CS%enth_prev(i,j,k) = IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1)
      do m=1,IG%NkIce
        CS%enth_prev(i,j,k) = CS%enth_prev(i,j,k) + &
                               (IST%mH_ice(i,j,k) * I_Nk) * IST%enth_ice(i,j,k,m)
      enddo
    endif ; enddo ; enddo ; enddo
  endif

  if (CS%debug_fast) &
    call IST_chksum("Start do_update_ice_model_fast", IST, G, US, IG)

  !$OMP parallel do default(shared) private(i2,j2,k2)
  do j=jsc,jec
    !   Set up local copies of fluxes.  The Atmos_boundary arrays may have
    ! different index conventions than are used internally in this component.
    do k=0,ncat ; do i=isc,iec
      i2 = i+i_off ; j2 = j+j_off ; k2 = k+1
      flux_u(i,j,k)  = US%kg_m2s_to_RZ_T*US%m_s_to_L_T*Atmos_boundary%u_flux(i2,j2,k2)
      flux_v(i,j,k)  = US%kg_m2s_to_RZ_T*US%m_s_to_L_T*Atmos_boundary%v_flux(i2,j2,k2)
      flux_sh(i,j,k)  = US%W_m2_to_QRZ_T*Atmos_boundary%t_flux(i2,j2,k2)
      evap(i,j,k)  = US%kg_m2s_to_RZ_T*Atmos_boundary%q_flux(i2,j2,k2)
      flux_lw(i,j,k) = US%W_m2_to_QRZ_T*Atmos_boundary%lw_flux(i2,j2,k2)
      flux_sw(i,j,k,NIR_DIR) = US%W_m2_to_QRZ_T*Atmos_boundary%sw_flux_nir_dir(i2,j2,k2)
      flux_sw(i,j,k,NIR_DIF) = US%W_m2_to_QRZ_T*Atmos_boundary%sw_flux_nir_dif(i2,j2,k2)
      flux_sw(i,j,k,VIS_DIR) = US%W_m2_to_QRZ_T*Atmos_boundary%sw_flux_vis_dir(i2,j2,k2)
      flux_sw(i,j,k,VIS_DIF) = US%W_m2_to_QRZ_T*Atmos_boundary%sw_flux_vis_dif(i2,j2,k2)
      lprec(i,j,k)   = US%kg_m2s_to_RZ_T*Atmos_boundary%lprec(i2,j2,k2)
      fprec(i,j,k)   = US%kg_m2s_to_RZ_T*Atmos_boundary%fprec(i2,j2,k2)
      dshdt(i,j,k) = US%W_m2_to_QRZ_T*Atmos_boundary%dhdt(i2,j2,k2)
      devapdt(i,j,k) = US%kg_m2s_to_RZ_T*Atmos_boundary%dedt(i2,j2,k2)
      dlwdt(i,j,k) = -1.*US%W_m2_to_QRZ_T*Atmos_boundary%drdt(i2,j2,k2)
    enddo ; enddo
    do i=isc,iec
      sh_T0(i,j,0) = flux_sh(i,j,0) - dshdt(i,j,0) * sOSS%SST_C(i,j)
      evap_T0(i,j,0) = evap(i,j,0) - devapdt(i,j,0) * sOSS%SST_C(i,j)
      lw_T0(i,j,0) = flux_lw(i,j,0) - dlwdt(i,j,0) * sOSS%SST_C(i,j)
    enddo
    do k=1,ncat ; do i=isc,iec
      sh_T0(i,j,k) = flux_sh(i,j,k) - dshdt(i,j,k) * Rad%t_skin(i,j,k)
      evap_T0(i,j,k) = evap(i,j,k) - devapdt(i,j,k) * Rad%t_skin(i,j,k)
      lw_T0(i,j,k) = flux_lw(i,j,k) - dlwdt(i,j,k) * Rad%t_skin(i,j,k)
    enddo ; enddo
  enddo

  if (CS%debug_fast) then
    call hchksum(flux_u(:,:,1:), "Mid do_fast flux_u", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    call hchksum(flux_v(:,:,1:), "Mid do_fast flux_v", G%HI, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    call hchksum(flux_sh(:,:,1:), "Mid do_fast flux_sh", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(evap(:,:,1:), "Mid do_fast evap", G%HI, scale=US%RZ_T_to_kg_m2s)
    call hchksum(flux_lw(:,:,1:), "Mid do_fast flux_lw", G%HI, scale=US%QRZ_T_to_W_m2)
    do b=1,size(flux_sw,4)
      write(nstr, '(I4)') b ; nstr = adjustl(nstr)
      call hchksum(flux_sw(:,:,1:,b), "Mid do_fast flux_sw("//trim(nstr)//")", G%HI, scale=US%QRZ_T_to_W_m2)
    enddo
    call hchksum(lprec(:,:,1:), "Mid do_fast lprec", G%HI, scale=US%RZ_T_to_kg_m2s)
    call hchksum(fprec(:,:,1:), "Mid do_fast fprec", G%HI, scale=US%RZ_T_to_kg_m2s)
    call hchksum(dshdt(:,:,1:), "Mid do_fast dshdt", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(devapdt(:,:,1:), "Mid do_fast devapdt", G%HI, scale=US%RZ_T_to_kg_m2s)
    call hchksum(dlwdt(:,:,1:), "Mid do_fast dlwdt", G%HI, scale=US%QRZ_T_to_W_m2)
  endif

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, Latent_vapor=LatHtVap)

  do j=jsc,jec ; do i=isc,iec
    flux_lh(i,j,0) = LatHtVap * evap(i,j,0)
  enddo ; enddo

  !
  ! implicit update of ice surface temperature
  !
  dt_fast = US%s_to_T*time_type_to_real(Time_step)

  !$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,NkIce,nb,IST,dshdt,devapdt,dlwdt, &
  !$OMP                                  flux_sw,flux_sh,evap,flux_lw,dt_fast,flux_lh,G,US,&
  !$OMP                                  S_col,I_Nk,LatHtVap,IG,sOSS,FIA,Rad,CS) &
  !$OMP                          private(latent,enth_col,sw_tot,dhf_dt,snow_wt,hf_0,Tskin,dTskin, &
  !$OMP                                  SW_abs_col,SW_absorbed,enth_here,tot_heat_in,enth_imb,&
  !$OMP                                  norm_enth_imb)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    if (IST%part_size(i,j,k) > 0.0) then
      enth_col(0) = IST%enth_snow(i,j,k,1)
      do m=1,NkIce ; enth_col(m) = IST%enth_ice(i,j,k,m) ; enddo

      ! This is for sublimation into water vapor at 0 degC; if the vapor should be
      ! at a different temperature, a correction would be made here.
      snow_wt = 0.0 ; if (IST%mH_snow(i,j,k)>0.0) snow_wt = 1.0
      latent = latent_sublimation(IST%enth_snow(i,j,k,1), IST%enth_ice(i,j,k,1), snow_wt, IST%ITV)

      sw_tot = 0.0
      do b=2,nb,2 ! This sum combines direct and diffuse fluxes to preserve answers.
        sw_tot = sw_tot + (flux_sw(i,j,k,b-1) + flux_sw(i,j,k,b))
      enddo
      sw_tot = ( (flux_sw(i,j,k,VIS_DIR) + flux_sw(i,j,k,VIS_DIF)) + &
               (flux_sw(i,j,k,NIR_DIR) + flux_sw(i,j,k,NIR_DIF)) )

      dhf_dt = (dshdt(i,j,k) + devapdt(i,j,k)*latent) - dlwdt(i,j,k)
      if (CS%Reorder_0C_heatflux) then
        ! This form seperately projects each contribution to the surface heat
        ! flux back to its value when T=0, so that a bitwise identical result
        ! can be obtained if the heating is redone.  The two forms are
        ! mathematically equivalent.
        hf_0 = ( ((flux_sh(i,j,k)  - dshdt(i,j,k) * Rad%t_skin(i,j,k)) + &
                (evap(i,j,k) - devapdt(i,j,k) * Rad%t_skin(i,j,k)) * latent) - &
               ((flux_lw(i,j,k) - dlwdt(i,j,k) * Rad%t_skin(i,j,k)) + &
                Rad%sw_abs_sfc(i,j,k) * sw_tot) )
      else
        hf_0 = ((flux_sh(i,j,k) + evap(i,j,k)*latent) - &
                (flux_lw(i,j,k) + Rad%sw_abs_sfc(i,j,k)* sw_tot)) - &
               dhf_dt * Rad%t_skin(i,j,k)
      endif

      SW_abs_col(0) = Rad%sw_abs_snow(i,j,k)*sw_tot
      do m=1,NkIce ; SW_abs_col(m) = Rad%sw_abs_ice(i,j,k,m)*sw_tot ; enddo

      !   This call updates the snow and ice temperatures and accumulates the
      ! surface and bottom melting/freezing energy.  The ice and snow do not
      ! actually lose or gain any mass from freezing or melting.
      ! mw/new - pass melt pond (surface temp fixed at freezing when present)
      call ice_temp_SIS2(IST%mH_pond(i,j,k), IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), &
                         enth_col, S_col, hf_0, dhf_dt, SW_abs_col, &
                         sOSS%T_fr_ocn(i,j), sOSS%bheat(i,j), Tskin, &
                         dt_fast, NkIce, FIA%tmelt(i,j,k), FIA%bmelt(i,j,k), &
                         CS%ice_thm_CSp, US, IST%ITV, CS%column_check)
      IST%enth_snow(i,j,k,1) = enth_col(0)
      do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_col(m) ; enddo

      if (CS%Reorder_0C_heatflux) then
        ! These extended expressions for the new fluxes will reproduce answers
        ! with redo_update_ice_model_fast.  They are mathematically equivalent
        ! to the next set of expressions.
        flux_sh(i,j,k) = (flux_sh(i,j,k) - dshdt(i,j,k) * Rad%t_skin(i,j,k)) + Tskin * dshdt(i,j,k)
        evap(i,j,k)    = (evap(i,j,k) - devapdt(i,j,k) * Rad%t_skin(i,j,k)) + Tskin * devapdt(i,j,k)
        flux_lw(i,j,k) = (flux_lw(i,j,k) - dlwdt(i,j,k) * Rad%t_skin(i,j,k)) + Tskin * dlwdt(i,j,k)
      else
        dTskin = Tskin - Rad%t_skin(i,j,k)
        flux_sh(i,j,k) = flux_sh(i,j,k) + dTskin * dshdt(i,j,k)
        evap(i,j,k)  = evap(i,j,k) + dTskin * devapdt(i,j,k)
        flux_lw(i,j,k) = flux_lw(i,j,k) + dTskin * dlwdt(i,j,k)
      endif
      flux_lh(i,j,k) = latent * evap(i,j,k)
      Rad%t_skin(i,j,k) = Tskin

      if (CS%column_check) then
        SW_absorbed = SW_abs_col(0)
        do m=1,NkIce ; SW_absorbed = SW_absorbed + SW_abs_col(m) ; enddo
        CS%heat_in(i,j,k) = CS%heat_in(i,j,k) + dt_fast * &
          ((flux_lw(i,j,k) + Rad%sw_abs_sfc(i,j,k)*sw_tot) + SW_absorbed + &
           sOSS%bheat(i,j) - (flux_sh(i,j,k) + flux_lh(i,j,k)))

        enth_here = (IST%mH_snow(i,j,k)) * enth_col(0)
        do m=1,NkIce
          enth_here = enth_here + (IST%mH_ice(i,j,k) * I_Nk) * enth_col(m)
        enddo
        tot_heat_in = CS%heat_in(i,j,k) - (FIA%bmelt(i,j,k) + FIA%tmelt(i,j,k))
        enth_imb = enth_here - (CS%enth_prev(i,j,k) + tot_heat_in)
        if (abs(enth_imb) > CS%imb_tol * (abs(enth_here) + &
                  abs(CS%enth_prev(i,j,k)) + abs(tot_heat_in)) ) then
          norm_enth_imb = enth_imb / (abs(enth_here) + abs(CS%enth_prev(i,j,k)) + abs(tot_heat_in))
          enth_imb = enth_here - (CS%enth_prev(i,j,k) + tot_heat_in)
        endif
      endif
    else ! IST%mH_ice <= 0
      flux_lh(i,j,k) = LatHtVap * evap(i,j,k)
    endif
  enddo ; enddo ; enddo

  call sum_top_quantities(FIA, Atmos_boundary, flux_u, flux_v, flux_sh, evap, &
                          flux_sw, flux_lw, lprec, fprec, flux_lh, Rad%t_skin, sOSS%SST_C, &
                          sh_T0, evap_T0, lw_T0, dshdt, devapdt, dlwdt, G, US, IG )

  if (CS%debug_fast) &
    call IST_chksum("End do_update_ice_model_fast", IST, G, US, IG)

  if (CS%bounds_check) &
    call IST_bounds_check(IST, G, US, IG, "End of update_ice_fast", Rad=Rad) !, OSS=sOSS)

end subroutine do_update_ice_model_fast

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> accumulate_deposition_fluxes accumulates any fluxes that are labeled as having
!! the type 'air_sea_deposition'.  These fluxes were not available when
!! update_ice_model_fast was called and the other fluxes were accumulated.
subroutine accumulate_deposition_fluxes(ABT, FIA, G, IG)
  type(atmos_ice_boundary_type), intent(in)    :: ABT !< A type containing atmospheric boundary
                                                      !! forcing fields that are used to drive the ice
  type(fast_ice_avg_type),       intent(inout) :: FIA !< A type containing averages of fields
                                                      !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),       intent(in)    :: G   !< The horizontal grid type
  type(ice_grid_type),           intent(in)    :: IG  !< The sea-ice specific grid type

  call coupler_type_increment_data( ABT%fluxes, FIA%tr_flux, only_flux_type='air_sea_deposition')

end subroutine accumulate_deposition_fluxes

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> redo_update_ice_model_fast applies the surface heat fluxes, shortwave radiation
!!   diffusion of heat to the sea-ice to implicitly determine a new temperature
!!   profile, subject to the constraint that ice and snow temperatures are never
!!   above freezing, using fluxes that have been determined during previous calls
!!   to do_update_ice_model_fast and stored in the fast_ice_avg_type.
subroutine redo_update_ice_model_fast(IST, sOSS, Rad, FIA, TSF, optics_CSp, &
                                      Time_step, CS, G, US, IG )
  type(ice_state_type),      intent(inout) :: IST        !< A type describing the state of the sea ice
  type(simple_OSS_type),     intent(in)    :: sOSS       !< A type describing the ocean surface state
  type(ice_rad_type),        intent(inout) :: Rad        !< A type containing fields related to
                                                         !! shortwave radiation in the ice
  type(fast_ice_avg_type),   intent(inout) :: FIA        !< A type containing averages of fields
                                                         !! (mostly fluxes) over the fast updates
  type(total_sfc_flux_type), intent(in)    :: TSF        !< A type with fluxes that are averaged across
                                                         !! the fast updates and integrated across thickness
                                                         !! categories from the fast ice update
  type(SIS_optics_CS),       intent(in)    :: optics_CSp !< The control structure for the sea ice optics module
  type(time_type),           intent(in)    :: Time_step  !< The amount of time over which to advance the ice
  type(fast_thermo_CS),      pointer       :: CS         !< The control structure for this module
  type(SIS_hor_grid_type),   intent(inout) :: G          !< The horizontal grid type
  type(unit_scale_type),     intent(in)    :: US         !< A structure with unit conversion factors
  type(ice_grid_type),       intent(in)    :: IG         !< The ice vertical grid type

  real, dimension(IG%NkIce)   :: &
    S_col         ! The thermodynamic salinity of a column of ice [gSalt kg-1].
  real, dimension(0:IG%NkIce) :: &
    T_col, &      ! The temperature of a column of ice and snow [degC].
    SW_abs_col, & ! The shortwave absorption within a column of snow and ice [Q R Z T-1 ~> W m-2].
    enth_col, &   ! The enthalpy of a column of snow and ice [Q ~> J kg-1].
    enth_col_in   ! The initial enthalpy of a column of snow and ice [Q ~> J kg-1].

  real :: dt_here ! The time step here [T ~> s].
  real :: Tskin   ! The new skin temperature [degC].
  real :: latent  ! The latent heat of sublimation of ice or snow [Q ~> J kg-1].
  real :: hf_0    ! The positive upward surface heat flux when T_sfc = 0 degC [Q R Z T-1 ~> W m-2].
  real :: dhf_dt  ! The derivative of the upward surface heat flux with Ts [Q R Z T-1 degC-1 ~> W m-2 degC-1].
  real :: sw_tot  ! sum over dir/dif vis/nir components [Q R Z T-1 ~> W m-2]
  real, dimension(size(FIA%flux_sw_top,4)) :: &
    albedos             ! The ice albedos by directional and wavelength band.
  real, dimension(IG%NkIce) :: &
    sw_abs_lay          ! The fractional shortwave absorption by each ice layer.
  real :: snow_wt       ! A fractional weighting of snow in the category surface area.
  real, dimension(G%isd:G%ied,size(FIA%flux_sw_top,4)) :: &
    sw_tot_ice_band     !   The total shortwave radiation by band, integrated
                        ! across the ice thickness partitions, but not the open
                        ! ocean partition [Q R Z T-1 ~> W m-2].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,IG%CatIce,size(FIA%flux_sw_top,4)) :: &
    sw_top_chg          !   The change in the shortwave down due to the new albedos [Q R Z T-1 ~> W m-2].
  real    :: flux_sw_prev  ! The previous value of flux_sw_top [Q R Z T-1 ~> W m-2].
  real    :: rescale    ! A rescaling factor between 0 and 1 [nondim].
  real    :: bmelt_tmp, tmelt_tmp ! Temporary arrays [Q R Z ~> J m-2].
  real    :: dSWt_dt    ! The derivative of SW_tot with skin temperature [Q R Z T-1 ~> W m-2 degC-1].
  real    :: Tskin_prev ! The previous value of Tskin [degC]
  real    :: T_bright   ! A skin temperature below which the snow and ice attain
                        ! their greatest brightness and albedo no longer varies [degC].
!  real    :: Tskin_itt(0:max(1,CS%max_tskin_itt))
!  real    :: SW_tot_itt(max(1,CS%max_tskin_itt))
  logical :: do_optics(G%isd:G%ied,G%jsd:G%jed)
  logical :: any_sw, any_ice
  logical :: do_any_j(G%jsd:G%jed)
  logical :: use_new_albedos
  integer :: i, j, k, m, i2, j2, k2, isc, iec, jsc, jec, ncat, NkIce
  integer :: b, b2, nb, nbmerge, itt, max_itt

  real :: ice_sw_tot ! The sum of shortwave fluxes into the ice and snow, but
                     ! excluding the fluxes transmitted to the ocean [Q R Z T-1 ~> W m-2].
  real :: TSF_sw_tot ! The total of all shortwave fluxes into the snow, ice,
                     ! and ocean that were previously stored in TSF [Q R Z T-1 ~> W m-2].
  real :: I_Nk       ! The inverse of the number of internal ice layers [nondim].

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_fast_thermo: Module must be initialized before it is used.")

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  nb = size(FIA%flux_sw_top,4)
  NkIce = IG%NkIce ; I_Nk = 1.0 / NkIce

  T_bright = bright_ice_temp(optics_CSp, IST%ITV)

  if (CS%debug_slow) then
    call IST_chksum("Start redo_update_ice_model_fast", IST, G, US, IG)
  endif

  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col)

  !
  ! implicit update of ice surface temperature
  !
  dt_here = US%s_to_T*time_type_to_real(Time_step)

  ! Ice can scatter direct shortwave into diffuse without loss of energy
  ! conservation, so it only the total of the shortwave in each frequency
  ! band that needs to be considered.  The effect of setting, nbmerge=2 is to
  ! combine direct and diffuse radiation bands of the same wavelength.
  nbmerge = 2
  ! Setting nbmerge=nb only considers the total of all shortwave bands, which
  ! is appropriate considering that all heat fluxes to the ocean are currently
  ! treated as diffuse visible light by SIS2.
  nbmerge = nb

  ! If there are multiple calls to redo_update_ice_model_fast between calls to
  ! slow_thermodynamics, this call would only occur during the first such call.
  FIA%tmelt(:,:,:) = 0.0 ; FIA%bmelt(:,:,:) = 0.0
  sw_top_chg(:,:,:,:) = 0.0
  use_new_albedos = .true.  ! Changing this value changes answers at the level of roundoff.
!  max_itt = 1

  ! This is here to reset the whole array after a restart.  It might be unnecessary.
  do k=1,ncat ; do j=jsc,jec ; do i=isc,iec
    Rad%sw_abs_ocn(i,j,k) = FIA%sw_abs_ocn(i,j,k)
  enddo ; enddo ; enddo

  !$OMP parallel do default(none) &
  !$OMP    shared( isc,iec,jsc,jec,nb,ncat,FIA,IST,do_any_j,do_optics) &
  !$OMP    private(any_sw,any_ice)
  do j=jsc,jec
    do_any_j(j) = .false.
    do i=isc,iec
      any_sw = .false.
      do b=1,nb ; if (FIA%flux_sw_dn(i,j,b) > 0.0) then
        any_sw = .true. ; exit
      endif ; enddo
      any_ice = .false.
      do k=1,ncat ; if (IST%part_size(i,j,k) > 0.0) then
        any_ice = .true. ; exit
      endif ; enddo
      do_optics(i,j) = (any_sw .and. any_ice)
      if (any_ice) do_any_j(j) = .true.
    enddo
  enddo

  if (CS%debug_slow) then
    call hchksum(Rad%coszen_lastrad, "Redo optics coszen_lastrad", G%HI)
    call hchksum(Rad%Tskin_rad, "Redo optics Tskin_rad", G%HI)
    call hchksum(FIA%flux_sw_dn, "Redo optics FIA%flux_sw_dn", G%HI, scale=US%QRZ_T_to_W_m2)
  endif

  !$OMP parallel do default(none) &
  !$OMP    shared( isc,iec,jsc,jec,nb,ncat,NkIce,FIA,IST,sOSS,Rad,US,IG,CS,optics_CSp,dt_here, &
  !$OMP            use_new_albedos,sw_top_chg,S_col,T_bright,max_itt,do_any_j,do_optics) &
  !$OMP    private(albedos,sw_abs_lay,flux_sw_prev,latent,enth_col,sw_tot,dSWt_dt,dhf_dt,hf_0, &
  !$OMP            Tskin,SW_abs_col,snow_wt,enth_col_in,tmelt_tmp,bmelt_tmp,Tskin_prev)
   !,Tskin_itt,SW_tot_itt)
  do j=jsc,jec ; if (do_any_j(j)) then
    ! Only work on j-rows with some ice in them.

    ! Determine the optical properties of the ice.  Because the optical properties
    ! can depend on the surface skin temperature, which is not yet well known,
    ! there may be some iteration for self-consistency.
    do k=1,ncat ; do i=isc,iec ; if (do_optics(i,j) .and. IST%part_size(i,j,k) > 0.0) then
      call ice_optics_SIS2(IST%mH_pond(i,j,k), IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), &
               Rad%Tskin_Rad(i,j,k), sOSS%T_fr_ocn(i,j), IG%NkIce, &
               albedos, Rad%sw_abs_sfc(i,j,k), Rad%sw_abs_snow(i,j,k), &
               sw_abs_lay, Rad%sw_abs_ocn(i,j,k), Rad%sw_abs_int(i,j,k), &
               US, optics_CSp, IST%ITV, coszen_in=Rad%coszen_lastrad(i,j))

      if (CS%max_Tskin_itt > 0) then
        ! Determine a new skin temperature that is consistent with the updated
        ! optics, and then recalculate the surface optical properties.  Ultimately
        ! this should be reiterated to convergence, perhaps using Newton's method.
        enth_col_in(0) = IST%enth_snow(i,j,k,1)
        do m=1,NkIce ; enth_col_in(m) = IST%enth_ice(i,j,k,m) ; enddo

        ! This is for sublimation into water vapor at 0 degC; if the vapor should be
        ! at a different temperature, a correction would be made here.
        snow_wt = 0.0 ; if (IST%mH_snow(i,j,k)>0.0) snow_wt = 1.0
        latent = latent_sublimation(IST%enth_snow(i,j,k,1), IST%enth_ice(i,j,k,1), snow_wt, IST%ITV)
        Tskin = Rad%Tskin_Rad(i,j,k)

!        ! These are here for debugging.
!        Tskin_itt(1:) = 10.0 ; SW_tot_itt(:) = 0.0
!        Tskin_itt(0) = Tskin

        do itt=1,CS%max_Tskin_itt
          sw_tot = 0.0 ; dSWt_dt = 0.0
          do b=1,nb
            sw_tot = sw_tot + (1.0 - albedos(b))*FIA%flux_sw_dn(i,j,b)
  !           dSWt_dt = dSWt_dt - dAlb_dt(b)*FIA%flux_sw_dn(i,j,b)
          enddo

          dhf_dt = (FIA%dshdt(i,j,k) + FIA%devapdt(i,j,k)*latent) - &
                 (FIA%dlwdt(i,j,k) + Rad%sw_abs_sfc(i,j,k)*dSWt_dt)
          hf_0 = (FIA%flux_sh0(i,j,k) + FIA%evap0(i,j,k)*latent) - &
                 (FIA%flux_lw0(i,j,k) + Rad%sw_abs_sfc(i,j,k)*sw_tot)

          SW_abs_col(0) = Rad%sw_abs_snow(i,j,k)*sw_tot
          do m=1,NkIce ; SW_abs_col(m) = sw_abs_lay(m)*sw_tot ; enddo

          do m=0,NkIce ; enth_col(m) = enth_col_in(m) ; enddo
          tmelt_tmp = 0.0 ; bmelt_tmp = 0.0 ; Tskin_prev = Tskin
          !   This call only estimates an updated skin temperature for
          ! calculating the ice optical properties.
          call ice_temp_SIS2(IST%mH_pond(i,j,k), IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), &
                   enth_col, S_col, hf_0, dhf_dt, SW_abs_col, &
                   sOSS%T_fr_ocn(i,j), sOSS%bheat(i,j), Tskin, &
                   0.5*dt_here, NkIce, tmelt_tmp, bmelt_tmp, CS%ice_thm_CSp, US, IST%ITV)

!         ! These are here to debug the iterations.
!         Tskin_itt(itt) = Tskin
!         SW_tot_itt(itt) = SW_tot

          call ice_optics_SIS2(IST%mH_pond(i,j,k), IST%mH_snow(i,j,k), &
                  IST%mH_ice(i,j,k), Tskin, sOSS%T_fr_ocn(i,j), IG%NkIce, &
                  albedos, Rad%sw_abs_sfc(i,j,k), Rad%sw_abs_snow(i,j,k), &
                  sw_abs_lay, Rad%sw_abs_ocn(i,j,k), Rad%sw_abs_int(i,j,k), &
                  US, optics_CSp, IST%ITV, coszen_in=Rad%coszen_lastrad(i,j))
          ! The feedbacks on the shortwave radiation are destabilizing, but only
          ! over a limited temperature range, so stop iterating (1) if the skin
          ! temperature is changing by a small enough amount, (2) there is
          ! surface melting, or (3) the skin is cold enough to be at a
          ! a temperature where the ice and snow have achieved maximal
          ! brightness.  Eventually this skin temperature dependence should
          ! be replaced by code that evolves the snow and ice grain sizes.
          if ((tmelt_tmp > 0.0) .or. (Tskin <= T_bright) .or. &
              (abs(Tskin_prev - Tskin) < 1e-5*abs(T_bright))) then
            exit
          endif
        enddo
!        if (itt > max_itt) then
!          max_itt = itt
!        endif
      endif

      if (use_new_albedos) then ; do b=1,nb ; if (FIA%flux_sw_dn(i,j,b) > 0.0) then
        flux_sw_prev = FIA%flux_sw_top(i,j,k,b)
        FIA%flux_sw_top(i,j,k,b) = (1.0 - albedos(b))*FIA%flux_sw_dn(i,j,b)
        sw_top_chg(i,j,k,b) = FIA%flux_sw_top(i,j,k,b) - flux_sw_prev
      endif ; enddo ; endif

      do m=1,IG%NkIce ; Rad%sw_abs_ice(i,j,k,m) = sw_abs_lay(m) ; enddo
    endif ; enddo ; enddo
  endif ; enddo

  ! The j-loops above and below could be combined, but they have been split to
  ! allow the following intermediate diagnostics to be added.
  if (CS%debug_slow) then
    call flux_redo_chksum("Middle redo_update_fast", IST, Rad, FIA, TSF, G, US, IG)
  endif

  !$OMP parallel do default(none) &
  !$OMP    shared( isc,iec,jsc,jec,nb,ncat,NkIce,FIA,IST,TSF,sOSS,Rad,IG,US,CS,dt_here, &
  !$OMP            nbmerge,S_col,do_any_j,do_optics) &
  !$OMP    private(rescale,sw_tot_ice_band,ice_sw_tot,TSF_sw_tot,latent,enth_col,sw_tot, &
  !$OMP            dhf_dt,hf_0,Tskin,SW_abs_col,snow_wt,enth_col_in)
  do j=jsc,jec ; if (do_any_j(j)) then
    ! Only work on j-rows with some ice in them.

    !    Determine whether the shortwave fluxes absorbed by the ice and snow in
    ! any shortwave frequency bands (or groups of bands) exceed the total
    ! shortwave absorption during the atmospheric steps, and if so scale them
    ! back for energy conservation.  If the shortwave absorption by the ice in
    ! any bands have decreased or increased only slightly, the difference will
    ! later be applied to the ocean.

    sw_tot_ice_band(:,:) = 0.0
    ! Note that the flux to the ocean is deliberately omitted here.
    ! Properly the radiative properties should be treated separately for each band.
    do k=1,ncat ; do b=1,nb ; do i=isc,iec
      sw_tot_ice_band(i,b) = sw_tot_ice_band(i,b) + IST%part_size(i,j,k) * &
               ((1.0 - Rad%sw_abs_ocn(i,j,k)) * FIA%flux_sw_top(i,j,k,b))
    enddo ; enddo ; enddo

    do i=isc,iec ; if (do_optics(i,j)) then ;  do b=1,nb,nbmerge
      ice_sw_tot = 0.0 ; TSF_sw_tot = 0.0
      do b2=0,nbmerge-1
        ice_sw_tot = ice_sw_tot + sw_tot_ice_band(i,b+b2)
        TSF_sw_tot = TSF_sw_tot + TSF%flux_sw(i,j,b+b2)
      enddo
      if (abs(ice_sw_tot) > abs(TSF_sw_tot)) then
        rescale = abs(TSF_sw_tot) /  abs(ice_sw_tot)
        ! Note that the shortwave flux to the ocean will be adjusted later
        ! so that the total shortwave heat fluxes agree with the initial
        ! calculation as passed from the atmosphere; that is where conservation
        ! is achieved.  This is the least aggressive rescaling that will avoid
        ! having negative shortwave fluxes into the ocean.
        do b2=0,nbmerge-1 ; do k=0,ncat
          FIA%flux_sw_top(i,j,k,b+b2) = rescale * FIA%flux_sw_top(i,j,k,b+b2)
        enddo ; enddo
      endif
    enddo ; endif ; enddo

    ! Now actually change the ice state due to the fast thermodynamic processes.
    do k=1,ncat ; do i=isc,iec ; if (IST%part_size(i,j,k) > 0.0) then
      enth_col(0) = IST%enth_snow(i,j,k,1)
      do m=1,NkIce ; enth_col(m) = IST%enth_ice(i,j,k,m) ; enddo

      ! This is for sublimation into water vapor at 0 degC; if the vapor should be
      ! at a different temperature, a correction would be made here.
      snow_wt = 0.0 ; if (IST%mH_snow(i,j,k)>0.0) snow_wt = 1.0
      latent = latent_sublimation(IST%enth_snow(i,j,k,1), IST%enth_ice(i,j,k,1), snow_wt, IST%ITV)

      sw_tot = (FIA%flux_sw_top(i,j,k,VIS_DIR) + FIA%flux_sw_top(i,j,k,VIS_DIF)) + &
               (FIA%flux_sw_top(i,j,k,NIR_DIR) + FIA%flux_sw_top(i,j,k,NIR_DIF))

      dhf_dt = (FIA%dshdt(i,j,k) + FIA%devapdt(i,j,k)*latent) - FIA%dlwdt(i,j,k)
      hf_0 = ( (FIA%flux_sh0(i,j,k) + FIA%evap0(i,j,k)*latent) - &
             (FIA%flux_lw0(i,j,k) + Rad%sw_abs_sfc(i,j,k)*sw_tot) )

      SW_abs_col(0) = Rad%sw_abs_snow(i,j,k)*sw_tot
      do m=1,NkIce ; SW_abs_col(m) = Rad%sw_abs_ice(i,j,k,m)*sw_tot ; enddo

      !   This call updates the snow and ice temperatures and accumulates the
      ! surface and bottom melting/freezing energy.  The ice and snow do not
      ! actually lose or gain any mass from freezing or melting.
      ! mw/new - pass melt pond (surface temp fixed at freezing when present)
      call ice_temp_SIS2(IST%mH_pond(i,j,k), IST%mH_snow(i,j,k), IST%mH_ice(i,j,k), &
                         enth_col, S_col, hf_0, dhf_dt, SW_abs_col, &
                         sOSS%T_fr_ocn(i,j), sOSS%bheat(i,j), Tskin, &
                         dt_here, NkIce, FIA%tmelt(i,j,k), FIA%bmelt(i,j,k), &
                         CS%ice_thm_CSp, US, IST%ITV, CS%column_check)
      IST%enth_snow(i,j,k,1) = enth_col(0)
      do m=1,NkIce ; IST%enth_ice(i,j,k,m) = enth_col(m) ; enddo

!      Rad%t_skin(i,j,k) = Tskin
      FIA%flux_sh_top(i,j,k)  = FIA%flux_sh0(i,j,k)  + Tskin * FIA%dshdt(i,j,k)
      FIA%evap_top(i,j,k)  = FIA%evap0(i,j,k) + Tskin * FIA%devapdt(i,j,k)
      FIA%flux_lw_top(i,j,k) = FIA%flux_lw0(i,j,k) + Tskin * FIA%dlwdt(i,j,k)
      FIA%flux_lh_top(i,j,k) = latent * FIA%evap_top(i,j,k)

      ! Copy radiation fields from the fast to the slow states.
      FIA%sw_abs_ocn(i,j,k) = Rad%sw_abs_ocn(i,j,k)
    endif ; enddo ; enddo
  endif ; enddo

  if (CS%debug_slow) &
    call IST_chksum("End redo_update_ice_model_fast", IST, G, US, IG)

  if (CS%bounds_check) &
    call IST_bounds_check(IST, G, US, IG, "End of redo_update_ice_fast", Rad=Rad) !, OSS=sOSS)

end subroutine redo_update_ice_model_fast


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Perform checksums on various arrays used in redoing the fast ice update.
subroutine flux_redo_chksum(mesg, IST, Rad, FIA, TSF, G, US, IG)
  character(len=*),          intent(in) :: mesg  !< A message that appears on the chksum lines.
  type(ice_state_type),      intent(in) :: IST   !< A type describing the state of the sea ice
  type(ice_rad_type),        intent(in) :: Rad   !< A type containing fields related to
                                                 !! shortwave radiation in the ice
  type(fast_ice_avg_type),   intent(in) :: FIA   !< A type containing averages of fields
                                                 !! (mostly fluxes) over the fast updates
  type(total_sfc_flux_type), intent(in) :: TSF   !< A type with fluxes that are averaged across
                                                 !! the fast updates and integrated across thickness
                                                 !! categories from the fast ice update
  type(SIS_hor_grid_type),   intent(inout) :: G  !< The ice-model's horizontal grid type.
  type(unit_scale_type),     intent(in)    :: US !< A structure with unit conversion factors
  type(ice_grid_type),       intent(in)    :: IG !< The ice vertical grid type

  real, dimension(G%isd:G%ied,G%jsd:G%jed,IG%CatIce) :: tmp_diag  ! A temporary diagnostic array
  character(len=8) :: nstr
  integer :: b, i, j, k, m

  call hchksum(IST%part_size(:,:,1:), trim(mesg)//" IST%part_size", G%HI)
  call hchksum(Rad%sw_abs_ocn(:,:,1:), trim(mesg)//" Rad%sw_abs_ocn", G%HI)
  do b=1,size(FIA%flux_sw_top,4)
    write(nstr, '(I4)') b ; nstr = adjustl(nstr)
    tmp_diag(:,:,:) = 0.0
    do k=1,IG%CatIce ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (IST%part_size(i,j,k) > 0.0) tmp_diag(i,j,k) = FIA%flux_sw_top(i,j,k,b)
    enddo ; enddo ; enddo
    call hchksum(tmp_diag, & ! similar to FIA%flux_sw_top(:,:,1:,b), &
                 trim(mesg)//" FIA%flux_sw_top("//trim(nstr)//")", G%HI, scale=US%QRZ_T_to_W_m2)
    call hchksum(TSF%flux_sw(:,:,b), &
                 trim(mesg)//" TSF%flux_sw("//trim(nstr)//")", G%HI, scale=US%QRZ_T_to_W_m2)
  enddo
  call hchksum(FIA%flux_sh0(:,:,1:), trim(mesg)//" FIA%flux_sh0", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%dshdt(:,:,1:), trim(mesg)//" FIA%dshdt", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%flux_lw0(:,:,1:), trim(mesg)//" FIA%flux_lw0", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%dlwdt(:,:,1:), trim(mesg)//" FIA%dlwdt", G%HI, scale=US%QRZ_T_to_W_m2)
  call hchksum(FIA%evap0(:,:,1:), trim(mesg)//" FIA%evap0", G%HI, scale=US%RZ_T_to_kg_m2s)
  call hchksum(FIA%devapdt(:,:,1:), trim(mesg)//" FIA%devapdt", G%HI, scale=US%RZ_T_to_kg_m2s)
  do m=1,size(Rad%sw_abs_ice,4)
    write(nstr, '(I4)') m ; nstr = adjustl(nstr)
    tmp_diag(:,:,:) = 0.0
    do k=1,IG%CatIce ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (IST%part_size(i,j,k) > 0.0) tmp_diag(i,j,k) = Rad%sw_abs_ice(i,j,k,m)
    enddo ; enddo ; enddo
    call hchksum(tmp_diag, & ! similar to Rad%sw_abs_ice(:,:,:,m), &
                 trim(mesg)//" Rad%sw_abs_ice("//trim(nstr)//")", G%HI)
  enddo
  call hchksum(Rad%sw_abs_snow(:,:,:), trim(mesg)//" Rad%sw_abs_snow", G%HI)
  call hchksum(Rad%sw_abs_ocn(:,:,:), trim(mesg)//" Rad%sw_abs_ocn", G%HI)

end subroutine flux_redo_chksum


!> Convert negative evaporation over ice (i.e. frost formation) into snow.
subroutine convert_frost_to_snow(FIA, G, US, IG)
  type(fast_ice_avg_type),       intent(inout) :: FIA !< A type containing averages of fields
                                                      !! (mostly fluxes) over the fast updates
  type(SIS_hor_grid_type),       intent(in)    :: G   !< The horizontal grid type
  type(unit_scale_type),         intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),           intent(in)    :: IG  !< The sea-ice specific grid type

  integer :: i, j, k, isc, iec, jsc, jec, ncat
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec ; if (FIA%evap_top(i,j,k) < 0.0) then
    FIA%fprec_top(i,j,k) = FIA%fprec_top(i,j,k) - FIA%evap_top(i,j,k)
    FIA%evap_top(i,j,k) = 0.0
  endif ; enddo ; enddo ; enddo

end subroutine convert_frost_to_snow

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_fast_thermo_init - initializes the parameters and diagnostics associated
!!    with the SIS_fast_thermo module.
subroutine SIS_fast_thermo_init(Time, G, IG, param_file, diag, CS)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid structure
  type(ice_grid_type),         intent(in)    :: IG   !< The sea-ice grid type
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(fast_thermo_CS),        pointer       :: CS   !< The control structure for the SIS_fast_thermo
                                                     !! module that is initialized here


! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40) :: mdl = "SIS_fast_thermo" ! This module's name.
  logical           :: debug

  call callTree_enter("SIS_fast_thermo_init(), SIS_fast_thermo.F90")

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_fast_thermo_init called with associated control structure.")
!    return
  else
    allocate(CS)
  endif

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, &
     "This module applies rapidly varying heat fluxes to the ice and does an "//&
     "implicit surface temperature calculation.")

  call get_param(param_file, mdl, "REORDER_0C_HEATFLUX", CS%Reorder_0C_heatflux, &
                 "If true, rearrange the calculation of the heat fluxes "//&
                 "projected back to 0C to work on each contribution "//&
                 "separately, so that they can be indentically replicated "//&
                 "if there is a single fast timestep per coupled timestep "//&
                 "and REDO_FAST_ICE_UPDATE=True.", default=.false.)
  call get_param(param_file, mdl, "MAX_TSKIN_ITT", CS%max_tskin_itt, &
                 "The maximum number of iterations of the skin temperature "//&
                 "and optical properties during redo_update_ice_model_fast.", &
                 default=10)
  call get_param(param_file, mdl, "COLUMN_CHECK", CS%column_check, &
                 "If true, add code to allow debugging of conservation "//&
                 "column-by-column.  This does not change answers, but "//&
                 "can increase model run time.", default=.false., &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "IMBALANCE_TOLERANCE", CS%imb_tol, &
                 "The tolerance for imbalances to be flagged by COLUMN_CHECK.", &
                 units="nondim", default=1.0e-9, debuggingParam=.true.)
  call get_param(param_file, mdl, "ICE_BOUNDS_CHECK", CS%bounds_check, &
                 "If true, periodically check the values of ice and snow "//&
                 "temperatures and thicknesses to ensure that they are "//&
                 "sensible, and issue warnings if they are not.  This "//&
                 "does not change answers, but can increase model run time.", &
                 default=.true.)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_SLOW_ICE", CS%debug_slow, &
                 "If true, write out verbose debugging data on the slow ice PEs.", &
                 default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_FAST_ICE", CS%debug_fast, &
                 "If true, write out verbose debugging data on the fast ice PEs.", &
                 default=debug, debuggingParam=.true.)

  call SIS2_ice_thm_init(G%US, param_file, CS%ice_thm_CSp)

  if (CS%column_check) then
    allocate(CS%enth_prev(G%HI%isd:G%HI%ied, G%HI%jsd:G%HI%jed, IG%CatIce)) ; CS%enth_prev(:,:,:) = 0.0
    allocate(CS%heat_in(G%HI%isd:G%HI%ied, G%HI%jsd:G%HI%jed, IG%CatIce)) ; CS%heat_in(:,:,:) = 0.0
  endif

  call callTree_leave("SIS_fast_thermo_init()")

end subroutine SIS_fast_thermo_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_fast_thermo_end deallocates any memory associated with this module.
subroutine SIS_fast_thermo_end(CS)
  type(fast_thermo_CS), pointer :: CS   !< The control structure for the SIS_slow_thermo
                                        !! module that is deallocated here

  call SIS2_ice_thm_end(CS%ice_thm_CSp)

  if (associated(CS)) deallocate(CS)

end subroutine SIS_fast_thermo_end

end module SIS_fast_thermo
