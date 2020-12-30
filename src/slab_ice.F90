!> Does the transport and redistribution between thickness categories for the SIS2 sea ice model.
module slab_ice

! This file is a part of SIS2.  See LICENSE.md for the license.

! use MOM_coms, only : reproducing_sum, EFP_type, EFP_to_real, EFP_real_diff
use MOM_domains,       only : pass_var, pass_vector, BGRID_NE, CGRID_NE
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
! use MOM_file_parser, only : get_param, log_param, read_param, log_version, param_file_type
use MOM_hor_index,   only : hor_index_type
use MOM_obsolete_params, only : obsolete_logical, obsolete_real
use MOM_unit_scaling,   only : unit_scale_type
! use SIS_diag_mediator, only : post_SIS_data, query_SIS_averaging_enabled, SIS_diag_ctrl
! use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field, time_type
use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type

implicit none ; private

#include <SIS2_memory.h>

public :: slab_ice_advect, slab_ice_dynamics

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> Advect an ice tracer or the thickness using a very old slab-ice algorithm
!! dating back to the Manabe model.
subroutine slab_ice_advect(uc, vc, trc, stop_lim, dt_slow, G, US, part_sz, nsteps)
  type(SIS_hor_grid_type),           intent(inout) :: G   !< The horizontal grid type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: uc  !< x-face advecting velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: vc  !< y-face advecting velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: trc !< Depth integrated amount of the tracer to advect,
                                                          !! in [Conc R Z ~> Conc kg m-2] or other units, or
                                                          !! the total ice mass [R Z ~> kg m-2].
  real,                              intent(in   ) :: stop_lim !< A tracer amount below which to
                                                          !! stop advection, in the same units as tr [Conc]
  real,                              intent(in   ) :: dt_slow !< The time covered by this call [T ~> s].
  type(unit_scale_type),             intent(in)    :: US  !< A structure with unit conversion factors
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(out) :: part_sz !< A part size that is set based on
                                                          !! whether trc (which may be mass) exceeds 0.
  integer,                          optional, intent(in ) :: nsteps !< The number of advective substeps.

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: uflx ! Zonal tracer fluxes [Conc R Z L2 T-1 ~> Conc kg s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: vflx ! Meridional tracer fluxes [Conc R Z L2 T-1 ~> Conc kg s-1]
  real :: avg, dif ! Average and forward difference of integrated tracer concentrations [Conc R Z ~> Conc kg m-2]
  real :: dt_adv  ! The advective timestep [T ~> s]
  integer :: i, j, n, isc, iec, jsc, jec, n_substeps
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  n_substeps = 1 ; if (present(nsteps)) n_substeps = nsteps
  if (n_substeps==0) return
  dt_adv = dt_slow / n_substeps

  do n=1,n_substeps
    do j=jsc,jec ; do I=isc-1,iec
      avg = 0.5*( trc(i,j) + trc(i+1,j) )
      dif = trc(i+1,j) - trc(i,j)
      if ( avg > stop_lim .and. uc(I,j) * dif > 0.0) then
        uflx(I,j) = 0.0
      elseif ( uc(i,j) > 0.0 ) then
        uflx(I,j) = uc(I,j) * trc(i,j) * G%dy_Cu(I,j)
      else
        uflx(I,j) = uc(I,j) * trc(i+1,j) * G%dy_Cu(I,j)
      endif
    enddo ; enddo

    do J=jsc-1,jec ; do i=isc,iec
      avg = 0.5*( trc(i,j) + trc(i,j+1) )
      dif = trc(i,j+1) - trc(i,j)
      if (avg > stop_lim .and. vc(i,J) * dif > 0.0) then
        vflx(i,J) = 0.0
      elseif ( vc(i,J) > 0.0 ) then
        vflx(i,J) = vc(i,J) * trc(i,j) * G%dx_Cv(i,J)
      else
        vflx(i,J) = vc(i,J) * trc(i,j+1) * G%dx_Cv(i,J)
      endif
    enddo ; enddo

    do j=jsc,jec ; do i=isc,iec
      trc(i,j) = trc(i,j) + dt_adv * ((uflx(I-1,j) - uflx(I,j)) + &
                                      (vflx(i,J-1) - vflx(i,J)) ) * G%IareaT(i,j)
    enddo ; enddo

    call pass_var(trc, G%Domain)
  enddo

  if (present(part_sz)) then ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
    part_sz(i,j) = 0.0 ; if (trc(i,j) > 0.0) part_sz(i,j) = 1.0
  enddo ; enddo ; endif

end subroutine slab_ice_advect

!> slab_ice_dynamics updates the B-grid or C-grid ice velocities and ice-ocean stresses as in the
!! very old slab-ice algorithm dating back to the Manabe model.  This code works for either
!! B-grid or C-grid discretiztions, but the velocity and stress variables must have consistent
!! array sizes and units.
subroutine slab_ice_dynamics(ui, vi, uo, vo, fxat, fyat, fxoc, fyoc)
  real, dimension(:,:), intent(inout) :: ui    !< Zonal ice velocity [L T-1 ~> m s-1]
  real, dimension(:,:), intent(inout) :: vi    !< Meridional ice velocity [L T-1 ~> m s-1]
  real, dimension(:,:), intent(in   ) :: uo    !< Zonal ocean velocity [L T-1 ~> m s-1]
  real, dimension(:,:), intent(in   ) :: vo    !< Meridional ocean velocity [L T-1 ~> m s-1]
  real, dimension(:,:), intent(in   ) :: fxat  !< Zonal air stress on ice [kg m-2 L T-2 ~> Pa]
  real, dimension(:,:), intent(in   ) :: fyat  !< Meridional air stress on ice [kg m-2 L T-2 ~> Pa]
  real, dimension(:,:), intent(  out) :: fxoc  !< Zonal ice stress on ocean [kg m-2 L T-2 ~> Pa]
  real, dimension(:,:), intent(  out) :: fyoc  !< Meridional ice stress on ocean [kg m-2 L T-2 ~> Pa]

  ui(:,:) = uo(:,:) ; vi(:,:) = vo(:,:)
  fxoc(:,:) = fxat(:,:) ; fyoc(:,:) = fyat(:,:)

end subroutine slab_ice_dynamics

end module slab_ice
