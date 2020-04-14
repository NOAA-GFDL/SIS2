!> sets up sea-ice specific grid information, including the category thicknesses and
!! vertical structure of the ice.
module ice_grid

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_grid - sets up sea-ice specific grid information, including the          !
!   category thicknesses and vertical structure of the ice. -Robert Hallberg   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_obsolete_params, only : obsolete_real
use MOM_unit_scaling, only : unit_scale_type

implicit none ; private

#include <SIS2_memory.h>
include 'netcdf.inc'

public :: set_ice_grid, ice_grid_end

!> Contains sea-ice specific grid elements, including information about the category descriptions
!! and the vertical layers in the snow and ice.
type, public :: ice_grid_type
  integer :: CatIce     !< The number of sea ice categories.
  integer :: NkIce      !< The number of vertical partitions within the sea ice.
  integer :: NkSnow     !< The number of vertical partitions within the snow atop the sea ice.
  real :: H_to_kg_m2    !< A constant that translates thicknesses from the internal units
                        !! of thickness to kg m-2 [kg m-2 R-1 Z-1 ~> 1].
  real :: kg_m2_to_H    !< A constant that translates thicknesses from kg m-2 to the internal
                        !! units of thickness [R Z m2 kg-1 ~> 1].
  real :: H_subroundoff !<   A thickness that is so small that it can be added to
                        !! any physically meaningful positive thickness without
                        !! changing it at the bit level [R Z ~> kg m-2].
  real :: ocean_part_min !< The minimum value for the fractional open-ocean area [nondim].  This can
                        !! be 0, but for some purposes it may be useful to set this to a miniscule
                        !! value (like 1e-40) that will be lost to roundoff during any sums so that
                        !! the open ocean fluxes can be used in interpolation across categories.
  real, allocatable, dimension(:) :: &
    cat_thick_lim, &  !< The lower thickness limits for each ice category [m].
    mH_cat_bound      !< The lower mass-per-unit area limits for each ice category [R Z ~> kg m-2].

end type ice_grid_type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> set_ice_grid initializes sea ice specific grid parameters
subroutine set_ice_grid(IG, US, param_file, NCat_dflt, ocean_part_min_dflt)
  type(ice_grid_type),   intent(inout) :: IG  !< The sea-ice specific grid type
  type(unit_scale_type), intent(in)    :: US  !< A structure with unit conversion factors
  type(param_file_type), intent(in)    :: param_file !< A structure to parse for run-time parameters
  integer,               intent(in)    :: NCat_dflt !< The default number of ice categories
  real, optional,        intent(in)    :: ocean_part_min_dflt !< The default value for the minimum open water area
!   This subroutine sets up the necessary domain types and the sea-ice grid.

! This include declares and sets the variable "version".
#include "version_variable.h"

! character(len=200) :: mesg
  character(len=40)  :: mod_nm  = "ice_grid" ! This module's name.
  real :: opm_dflt

  opm_dflt = 0.0 ; if (present(ocean_part_min_dflt)) opm_dflt = ocean_part_min_dflt

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod_nm, version)
#ifdef STATIC_MEMORY_
  call get_param(param_file, mod_nm, "NCAT_ICE", IG%CatIce, &
                 "The number of sea ice thickness categories.", units="nondim", &
                 default=NCat_dflt)
  if (IG%CatIce /= NCAT_ICE_) call SIS_error(FATAL, "set_ice_grid: " // &
       "Mismatched number of categories NCAT_ICE between SIS_memory.h and "//&
       "param_file or the input namelist file.")
  call get_param(param_file, mod_nm, "NK_ICE", IG%NkIce, &
                 "The number of layers within the sea ice.", units="nondim", &
                 default=NK_ICE_)
  if (IG%NkIce /= NK_ICE_) call SIS_error(FATAL, "set_ice_grid: " // &
       "Mismatched number of layers NK_ICE between SIS_memory.h and param_file")

  call get_param(param_file, mod_nm, "NK_SNOW", IG%NkSnow, &
                 "The number of layers within the snow atop the sea ice.", &
                 units="nondim", default=NK_SNOW_)
  if (IG%NkSnow /= NK_SNOW_) call SIS_error(FATAL, "set_ice_grid: " // &
       "Mismatched number of layers NK_SNOW between SIS_memory.h and param_file")
  if (global_indexing) cal SIS_error(FATAL, "set_ice_grid : "//&
       "GLOBAL_INDEXING can not be true with STATIC_MEMORY.")
#else
  call get_param(param_file, mod_nm, "NCAT_ICE", IG%CatIce, &
                 "The number of sea ice thickness categories.", units="nondim", &
                 default=NCat_dflt)
  call get_param(param_file, mod_nm, "NK_ICE", IG%NkIce, &
                 "The number of layers within the sea ice.", units="nondim", &
                 default=4) ! Valid for SIS5L; Perhaps this should be ..., fail_if_missing=.true.
  call get_param(param_file, mod_nm, "NK_SNOW", IG%NkSnow, &
                 "The number of layers within the snow atop the sea ice.", &
                 units="nondim", default=1) ! Perhaps this should be ..., fail_if_missing=.true.
#endif

  call obsolete_real(param_file, "H_TO_KG_M2", warning_val=1.0)
  IG%H_to_kg_m2 = US%RZ_to_kg_m2
  IG%kg_m2_to_H = US%kg_m3_to_R * US%m_to_Z
  IG%H_subroundoff = 1.0e-30*US%kg_m3_to_R*US%m_to_Z

  call get_param(param_file, mod_nm, "MIN_OCEAN_PARTSIZE", IG%ocean_part_min, &
                 "The minimum value for the fractional open-ocean area. "//&
                 "This can be 0, but for some purposes it may be useful "//&
                 "to set this to a miniscule value (like 1e-40) that will "//&
                 "be lost to roundoff during any sums so that the open "//&
                 "ocean fluxes can be used in with new categories.", &
                 units="nondim", default=opm_dflt)
  if ((IG%ocean_part_min < 0.0) .or. (IG%ocean_part_min > 1.0e-20)) &
    call SIS_error(FATAL, "MIN_OCEAN_PARTSIZE has been set outside of the valid"//&
                   "range of 0 to 1e-20.")

  call allocate_ice_metrics(IG)

end subroutine set_ice_grid


!> Allocate any required arrays in the ice_grid_type.
subroutine allocate_ice_metrics(IG)
  type(ice_grid_type), intent(inout) :: IG  !< The sea-ice specific grid type
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, isg, ieg, jsg, jeg

  ! This subroutine allocates any extensive elements of the ice_grid_type
  ! and zeros them out.
  allocate(IG%cat_thick_lim(1:IG%CatIce+1)) ; IG%cat_thick_lim(:) = 0.0
  allocate(IG%mH_cat_bound(1:IG%CatIce+1)) ; IG%mH_cat_bound(:) = 0.0
end subroutine allocate_ice_metrics

!---------------------------------------------------------------------
!> Release memory used by the ice_grid_type and related structures.
subroutine ice_grid_end(IG)
  type(ice_grid_type), intent(inout) :: IG  !< The sea-ice specific grid type

  deallocate(IG%cat_thick_lim)
  deallocate(IG%mH_cat_bound)

end subroutine ice_grid_end

end module ice_grid
