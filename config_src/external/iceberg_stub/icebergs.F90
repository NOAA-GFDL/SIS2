!> Top-level/entry functions that step forward the governing equations
module ice_bergs

! This file is part of NOAA-GFDL/icebergs. See LICENSE.md for the license.

use time_manager_mod, only: time_type
use ice_bergs_framework, only: icebergs

implicit none ; private

!ice_model.F90:49:use ice_bergs,          only : icebergs, icebergs_run, icebergs_init, icebergs_end
!ice_type.F90:4:use ice_bergs,         only : icebergs, icebergs_stock_pe, icebergs_save_restart
!SIS_dyn_trans.F90:66:use ice_bergs,         only : icebergs, icebergs_run, icebergs_init, icebergs_end

public icebergs_init, icebergs_end, icebergs_run, icebergs_stock_pe, icebergs
public icebergs_save_restart

#ifndef _FILE_VERSION
! Version of file provided can be set to git hash via a CPP macro but if not set we use 'unknown'
#define _FILE_VERSION 'unknown'
#endif
character(len=128) :: version = _FILE_VERSION !< Version of file

contains

!> Initializes icebergs container "bergs"
subroutine icebergs_init(bergs, &
             gni, gnj, layout, io_layout, axes, dom_x_flags, dom_y_flags, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot, ocean_depth, maskmap, fractional_area)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  integer, intent(in) :: gni !< Number of global points in i-direction
  integer, intent(in) :: gnj !< Number of global points in j-direction
  integer, intent(in) :: layout(2) !< Parallel decomposition of computational processors in i/j direction
  integer, intent(in) :: io_layout(2) !< Parallel decomposition of i/o processors in i/j direction
  integer, intent(in) :: axes(2) !< Diagnostic axes
  integer, intent(in) :: dom_x_flags !< Decomposition flags for i-direction
  integer, intent(in) :: dom_y_flags !< Decomposition flags for j-direction
  real, intent(in) :: dt !< Time step (s)
  type (time_type), intent(in) :: Time !< Model time
  real, dimension(:,:), intent(in) :: ice_lon !< Longitude of cell corners using NE convention (degree E)
  real, dimension(:,:), intent(in) :: ice_lat !< Latitude of cell corners using NE conventino (degree N)
  real, dimension(:,:), intent(in) :: ice_wet !< Wet/dry mask (1 is wet, 0 is dry) of cell centers
  real, dimension(:,:), intent(in) :: ice_dx !< Zonal length of cell on northern side (m)
  real, dimension(:,:), intent(in) :: ice_dy !< Meridional length of cell on eastern side (m)
  real, dimension(:,:), intent(in) :: ice_area !< Area of cells (m^2, or non-dim is fractional_area=True)
  real, dimension(:,:), intent(in) :: cos_rot !< Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), intent(in) :: sin_rot !< Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), intent(in),optional :: ocean_depth !< Depth of ocean bottom (m)
  logical, intent(in), optional :: maskmap(:,:) !< Masks out parallel cores
  logical, intent(in), optional :: fractional_area !< If true, ice_area contains cell area as fraction of entire spherical surface

end subroutine icebergs_init

!> The main driver the steps updates icebergs
subroutine icebergs_run(bergs, time, calving, uo, vo, ui, vi, tauxa, tauya, ssh, sst, calving_hflx, cn, hi, &
                        stagger, stress_stagger, sss, mass_berg, ustar_berg, area_berg)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  type(time_type), intent(in) :: time !< Model time
  real, dimension(:,:), intent(inout) :: calving !< Calving (kg/s). This field is updated with melt by bergs.
  real, dimension(:,:), intent(inout) :: calving_hflx !< Calving heat flux (W/m2)
  real, dimension(:,:), intent(in) :: uo !< Ocean zonal velocity (m/s)
  real, dimension(:,:), intent(in) :: vo !< Ocean meridional velocity (m/s)
  real, dimension(:,:), intent(in) :: ui !< Ice zonal velocity (m/s)
  real, dimension(:,:), intent(in) :: vi !< Ice meridional velocity (m/s)
  real, dimension(:,:), intent(in) :: tauxa !< Zonal wind stress (Pa)
  real, dimension(:,:), intent(in) :: tauya !< Meridional wind stress (Pa)
  real, dimension(:,:), intent(in) :: ssh !< Effective sea-surface height (m)
  real, dimension(:,:), intent(in) :: sst !< Sea-surface temperature (C or K)
  real, dimension(:,:), intent(in) :: cn !< Sea-ice concentration (nondim)
  real, dimension(:,:), intent(in) :: hi !< Sea-ice thickness (m)
  integer, optional, intent(in) :: stagger !< Enumerated value indicating staggering of ocean/ice u,v variables
  integer, optional, intent(in) :: stress_stagger !< Enumerated value indicating staggering of stress variables
  real, dimension(:,:), optional, intent(in) :: sss !< Sea-surface salinity (1e-3)
  real, dimension(:,:), optional, pointer :: mass_berg !< Mass of bergs (kg)
  real, dimension(:,:), optional, pointer :: ustar_berg !< Friction velocity on base of bergs (m/s)
  real, dimension(:,:), optional, pointer :: area_berg !< Area of bergs (m2)

end subroutine icebergs_run

!> Calculate stocks of water and heat
subroutine icebergs_stock_pe(bergs, index, value)
  ! Modules
  !use stock_constants_mod, only : ISTOCK_WATER, ISTOCK_HEAT
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  integer, intent(in) :: index !< =ISTOCK_WATER or ISTOCK_HEAT
  real, intent(out) :: value !< Amount of ice or water

end subroutine icebergs_stock_pe

!> Write restart files
subroutine icebergs_save_restart(bergs, time_stamp)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory
  character(len=*),    intent(in), optional :: time_stamp !< Timestamp for restart file

end subroutine icebergs_save_restart

!> Deallocate all memory and disassociated pointer
subroutine icebergs_end(bergs)
  ! Arguments
  type(icebergs), pointer :: bergs !< Container for all types and memory

end subroutine icebergs_end

end module
