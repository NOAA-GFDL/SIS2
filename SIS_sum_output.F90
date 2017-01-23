module SIS_sum_output
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
!   This file contains the subroutines that calculate globally integrated      !
! sea-ice quantities for SIS2, and writes them to a netcdf file and and an     !
! ASCII output file.  This code was originally adapted from MOM_sum_output.F90 !
! by Robert Hallberg in May 2014.                                              !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_coms, only : sum_across_PEs, PE_here, root_PE, num_PEs, max_across_PEs
use MOM_coms, only : reproducing_sum
use MOM_coms, only : EFP_type, operator(+), operator(-), assignment(=), EFP_to_real, real_to_EFP
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
! use MOM_io, only : create_file, fieldtype, flush_file, reopen_file, vardesc, write_field
use MOM_io, only : open_file
use MOM_io, only : APPEND_FILE, ASCII_FILE, SINGLE_FILE, WRITEONLY_FILE
use MOM_string_functions, only : slasher
use MOM_time_manager, only : time_type, get_time, set_time, operator(>), operator(-)
use MOM_time_manager, only : get_date, get_calendar_type, NO_CALENDAR
! use MOM_tracer_flow_control, only : tracer_flow_control_CS, call_tracer_stocks

use SIS_types, only : ice_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_types, only : ocean_sfc_state_type
use SIS_hor_grid, only : SIS_hor_grid_type
use ice_grid, only : ice_grid_type
use SIS2_ice_thm, only : enth_from_TS, get_SIS2_thermo_coefs, ice_thermo_type
use SIS_sum_output_type, only : SIS_sum_out_CS
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS, SIS_call_tracer_stocks

use netcdf

implicit none ; private

#include <SIS2_memory.h>

public write_ice_statistics, accumulate_bottom_input
public SIS_sum_output_init, SIS_sum_output_end, SIS_sum_out_CS
public accumulate_input_1, accumulate_input_2

!-----------------------------------------------------------------------

! integer, parameter :: NUM_FIELDS = 17

! type, public :: SIS_sum_out_CS ; private

!   real    :: fresh_water_input  !   The total mass of fresh water added by
!                                 ! surface fluxes since the last time that
!   real    :: mass_prev          !   The total sea ice mass the last time that
!                                 ! write_ice_statistics was called, in kg.
!   real    :: salt_prev          !   The total amount of salt in the sea ice the last
!                                 ! time that write_ice_statistics was called, in PSU kg.
!   real    :: net_salt_input     !   The total salt added by surface fluxes since
!                                 ! the last time that write_ice_statistics was called,
!                                 ! in PSU kg.
!   real    :: heat_prev          !   The total amount of heat in the sea ice the last
!                                 ! time that write_ice_statistics was called, in Joules.
!   real    :: net_heat_input     !   The total heat added by surface fluxes since
!                                 ! the last time that write_ice_statistics was called,
!                                 ! in Joules.
!   type(EFP_type) :: &
!     fresh_water_in_EFP, &       ! These are extended fixed point versions of the
!     net_salt_in_EFP, &          ! correspondingly named variables above.
!     net_heat_in_EFP, heat_prev_EFP, salt_prev_EFP, mass_prev_EFP
!   real    :: dt                 ! The baroclinic dynamics time step, in s.
!   real    :: timeunit           !   The length of the units for the time
!                                 ! axis, in s.
!   type(time_type) :: Start_time ! The start time of the simulation.
!                                 ! Start_time is set in MOM_initialization.F90
!   logical :: write_stdout       ! If true, periodically write sea ice statistics
!                                 ! to stdout to allow the progress to be seen.
!   logical :: write_stocks       ! If true, write the integrated tracer amounts
!                                 ! to stdout when the statistics files are written.
!   integer :: previous_calls = 0 ! The number of times write_ice_statistics has been called.
!   integer :: prev_n = 0         ! The value of n from the last call.
! !  integer :: statsfile_nc       ! NetCDF id of the statistics file.
!   integer :: statsfile_ascii    ! The unit number of the ascii version of the statistics file.
! !  type(fieldtype), dimension(NUM_FIELDS+MAX_FIELDS_) :: &
! !             fields             ! fieldtype variables for the output fields.
!   character(len=200) :: statsfile  ! The name of the statistics file with path.
! end type SIS_sum_out_CS

contains

subroutine SIS_sum_output_init(G, param_file, directory, Input_start_time, CS, &
                               ntrunc)
  type(SIS_hor_grid_type),  intent(in)    :: G
  type(param_file_type),    intent(in)    :: param_file
  character(len=*),         intent(in)    :: directory
  type(time_type),          intent(in)    :: Input_start_time
  type(SIS_sum_out_CS),     pointer       :: CS
  integer, target, optional,intent(inout) :: ntrunc
! Arguments: G - The sea ice model's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      directory - The directory where the statistics file goes.
!  (in)      Input_start_time - The start time of the simulation.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
!  (in/out,opt)  ntrunc - The integer that stores the number of times the velocity
!                     has been truncated since the last call to write_ice_statistics.
  real :: Rho_0, maxvel
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_sum_output" ! This module's name.
  character(len=200) :: statsfile  ! The name of the statistics file.

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_sum_output_init called with associated control structure.")
    return
  endif
  allocate(CS)

  if (present(ntrunc)) then ; CS%ntrunc => ntrunc ; else ; allocate(CS%ntrunc) ; endif
  CS%ntrunc = 0

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "WRITE_STOCKS", CS%write_stocks, &
                 "If true, write the integrated tracer amounts to stdout \n"//&
                 "when the statistics files are written.", default=.true.)
  call get_param(param_file, mod, "STDOUT_HEARTBEAT", CS%write_stdout, &
                 "If true, periodically write sea ice statistics to \n"//&
                 "stdout to allow the progress to be seen.", default=.true.)
  call get_param(param_file, mod, "DT_ICE_DYNAMICS", CS%dt, &
                 "The time step used for the slow ice dynamics, including "//&
                 "stepping the continuity equation and interactions between "//&
                 "the ice mass field and velocities.", units="s", &
                 default=-1.0, do_not_log=.true.)
  call get_param(param_file, mod, "MAXTRUNC", CS%maxtrunc, &
                 "The run will be stopped, and the day set to a very \n"//&
                 "large value if the velocity is truncated more than \n"//&
                 "MAXTRUNC times between  writing ice statistics. \n"//&
                 "Set MAXTRUNC to 0 to stop if there is any truncation \n"//&
                 "of sea ice velocities.", units="truncations save_interval-1", default=0)

  call get_param(param_file, mod, "STATISTICS_FILE", statsfile, &
                 "The file to use to write the globally integrated \n"//&
                 "statistics.", default="seaice.stats")

  CS%statsfile = trim(slasher(directory))//trim(statsfile)
  call log_param(param_file, mod, "output_path/STATISTICS_FILE", CS%statsfile)
#ifdef STATSLABEL
  CS%statsfile = trim(CS%statsfile)//"."//trim(adjustl(STATSLABEL))
#endif

  call get_param(param_file, mod, "TIMEUNIT", CS%Timeunit, &
                 "The time unit in seconds a number of input fields", &
                 units="s", default=86400.0)
  if (CS%Timeunit < 0.0) CS%Timeunit = 86400.0
  call get_param(param_file, mod, "COLUMN_CHECK", CS%column_check, &
                 "If true, add code to allow debugging of conservation \n"//&
                 "column-by-column.  This does not change answers, but \n"//&
                 "can increase model run time.", default=.false.)
  call get_param(param_file, mod, "IMBALANCE_TOLERANCE", CS%imb_tol, &
                 "The tolerance for imbalances to be flagged by COLUMN_CHECK.", &
                 units="nondim", default=1.0e-9)

  CS%Start_time = Input_start_time

  allocate(CS%water_in_col(G%isd:G%ied, G%jsd:G%jed)) ; CS%water_in_col(:,:) = 0.0
  allocate(CS%heat_in_col(G%isd:G%ied, G%jsd:G%jed))  ; CS%heat_in_col(:,:) = 0.0
  allocate(CS%salt_in_col(G%isd:G%ied, G%jsd:G%jed))  ; CS%salt_in_col(:,:) = 0.0
  if (CS%column_check) then
    allocate(CS%water_col_prev(G%isd:G%ied, G%jsd:G%jed)) ; CS%water_col_prev(:,:) = 0.0
    allocate(CS%heat_col_prev(G%isd:G%ied, G%jsd:G%jed))  ; CS%heat_col_prev(:,:) = 0.0
    allocate(CS%salt_col_prev(G%isd:G%ied, G%jsd:G%jed))  ; CS%salt_col_prev(:,:) = 0.0
  endif

end subroutine SIS_sum_output_init

subroutine SIS_sum_output_end(CS)
  type(SIS_sum_out_CS), pointer :: CS
!   This subroutine deallocates the memory owned by this module.
! Argument: CS - The control structure returned by a previous call to
!                SIS_sum_output_init.

  if (associated(CS)) then
    deallocate(CS)
  endif
end subroutine SIS_sum_output_end

subroutine write_ice_statistics(IST, day, n, G, IG, CS, message, check_column, tracer_CSp)
  type(ice_state_type),    intent(inout) :: IST

  type(time_type),         intent(inout) :: day
  integer,                 intent(in)    :: n
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(inout) :: IG
  type(SIS_sum_out_CS),    pointer       :: CS
  character(len=*), optional, intent(in) :: message
  logical,          optional, intent(in) :: check_column
  type(SIS_tracer_flow_control_CS), optional, pointer       :: tracer_CSp

!  This subroutine calculates and writes the total sea-ice mass by
! hemisphere, heat, salt, and other globally integrated quantities.

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields, including potential temperature and
!                 salinity or mixed layer density. Absent fields have NULL ptrs.
!  (in/out)  day - The current model time.
!  (in)      n - The time step number of the current execution.
!  (in)      G - The sea ice model's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 SIS_sum_output_init.
!  (in,opt)  message - A text message to use with this output.
!  (in,opt)  check_column - If true, check for column-wise heat and mass conservation.

  real, dimension(SZI_(G),SZJ_(G), 2) :: &
    ice_area, &    ! The area of ice in each cell and hemisphere, in m2.
    ice_extent, &  ! The extent (cells with >10% coverage) of ice in each
                   ! cell and hemisphere, in m2.
    col_mass, &    ! The column integrated ice and snow mass in each cell and
                   ! hemisphere, in kg.
    col_heat, &    ! The column integrated ice and snow heat in each cell and
                   ! hemisphere, in J.
    col_salt       ! The column integrated salt in the ice in each cell and
                   ! hemisphere in kg.

  real, dimension(2) :: &
    Area_NS, &     ! The total sea-ice area in the two hemispheres, in m2.
    Extent_NS, &   ! The total sea-ice extent in the two hemispheres, in m2.
    Heat_NS, &     ! The total sea-ice enthalpy in the two hemispheres, in J.
    mass_NS, &     ! The total sea-ice mass in the two hemispheres, in kg.
    salt_NS, &     ! The total sea-ice salt in the two hemispheres, in kg.
    salinity_NS    ! The average sea-ice salinity in the two hemispheres, in g/kg.

  real :: Mass         ! The total mass of the sea ice and snow atop it in kg.
  real :: mass_chg     ! The change in total sea ice mass of fresh water since
                       ! the last call to this subroutine, in kg.
  real :: mass_anom    ! The change in fresh water that cannot be accounted for
                       ! by the surface fluxes, in kg.
  real :: I_Mass       ! Adcroft's rule reciprocal of mass: 1/Mass or 0, in kg-1.
  real :: Salt         ! The total amount of salt in the ocean, in PSU kg.
  real :: Salt_chg     ! The change in total sea ice salt since the last call
                       ! to this subroutine, in PSU kg.
  real :: Salt_anom    ! The change in salt that cannot be accounted for by
                       ! the surface fluxes, in PSU kg.
  real :: Salt_anom_norm ! The salt anomaly normalized by salt (if it is nonzero).
  real :: salin        ! The mean salinity of the ocean, in PSU.
  real :: salin_chg    ! The change in total salt since the last call
                       ! to this subroutine divided by total mass, in PSU.
  real :: salin_anom   ! The change in total salt that cannot be accounted for by
                       ! the surface fluxes divided by total mass in PSU.
  real :: salin_mass_in ! The mass of salt input since the last call, kg.
  real :: Heat         ! The total amount of Heat in the ocean, in Joules.
  real :: Heat_chg     ! The change in total sea ice heat since the last call
                       ! to this subroutine, in Joules.
  real :: Heat_anom    ! The change in heat that cannot be accounted for by
                       ! the surface fluxes, in Joules.
  real :: Heat_anom_norm ! The heat anomaly normalized by heat (if it is nonzero).
  real :: temp         ! The mean potential temperature of the ocean, in C.
  real :: temp_anom    ! The change in total heat that cannot be accounted for
                       ! by the surface fluxes, divided by the total heat
                       ! capacity of the ocean, in C.
  real :: Area         ! The total area of the sea ice in m2.
  real :: Extent       ! The total extent of the sea ice in m2.
  real :: heat_imb     ! The column integrated heat imbalance in enth_unit kg m-2.
  real :: mass_imb     ! The column integrated mass imbalance in kg.
  real :: enth_liq_0   ! The enthalpy of liquid water at the freezing point, in enth_unit.
  real :: I_nlay, kg_H_nlay, area_pt
  real :: area_h       ! The masked area of a column.
  type(EFP_type) :: &
    mass_EFP, &        ! Extended fixed point sums of total mass, etc.
    salt_EFP, heat_EFP, salt_chg_EFP, heat_chg_EFP, mass_chg_EFP, &
    mass_anom_EFP, salt_anom_EFP, heat_anom_EFP

  real :: CFL_trans    ! A transport-based definition of the CFL number, nondim.
  real :: CFL_u, CFL_v ! Simple CFL numbers for u- and v- advection, nondim.
  real :: dt_CFL       ! The timestep for calculating the CFL number, in s.
  real :: max_CFL      ! The maximum of the CFL numbers, nondim.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    Temp_int, Salt_int
  logical :: check_col
  integer :: num_nc_fields  ! The number of fields that will actually go into
                            ! the NetCDF file.
  integer :: i, j, k, is, ie, js, je, L, m, nlay, ncat, hem
  integer :: start_of_day, num_days
  integer :: iyear, imonth, iday, ihour, iminute, isecond, itick ! For call to get_date()
  real    :: reday, var
  character(len=120) :: statspath_nc
  character(len=300) :: mesg
  character(len=48)  :: msg_start
  character(len=32)  :: mesg_intro, time_units, day_str, n_str, trunc_str

  integer :: isc, iec, jsc, jec

  real :: Tr_stocks(MAX_FIELDS_)
  character(len=40), dimension(MAX_FIELDS_) :: &
      Tr_names, Tr_units
  integer :: nTr_stocks


! real :: Tr_stocks(MAX_FIELDS_)
! real :: Tr_min(MAX_FIELDS_),Tr_max(MAX_FIELDS_)
! real :: Tr_min_x(MAX_FIELDS_), Tr_min_y(MAX_FIELDS_), Tr_min_z(MAX_FIELDS_)
! real :: Tr_max_x(MAX_FIELDS_), Tr_max_y(MAX_FIELDS_), Tr_max_z(MAX_FIELDS_)
! logical :: Tr_minmax_got(MAX_FIELDS_) = .false.
! character(len=40), dimension(MAX_FIELDS_) :: &
!   Tr_names, Tr_units
! integer :: nTr_stocks

! A description for output of each of the fields.
!  type(vardesc) :: vars(NUM_FIELDS+MAX_FIELDS_)

!  num_nc_fields = 17
! vars(6) = vardesc("Mass_cat","Total Ice Mass by Category",'1','L','s',"kg")
! vars(7) = vardesc("Mass","Total Mass",'1','1','s',"kg")
! vars(8) = vardesc("Mass_chg","Total Mass Change between Entries",'1','1','s',"kg")
! vars(9) = vardesc("Mass_anom","Anomalous Total Mass Change",'1','1','s',"kg")

! vars(12) = vardesc("Salt","Total Salt",'1','1','s',"kg")
! vars(13) = vardesc("Salt_chg","Total Salt Change between Entries",'1','1','s',"kg")
! vars(14) = vardesc("Salt_anom","Anomalous Total Salt Change",'1','1','s',"kg")
! vars(15) = vardesc("Heat","Total Heat",'1','1','s',"Joules")
! vars(16) = vardesc("Heat_chg","Total Heat Change between Entries",'1','1','s',"Joules")
! vars(17) = vardesc("Heat_anom","Anomalous Total Heat Change",'1','1','s',"Joules")

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
!  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce ; nlay = IG%NkIce
  check_col = .false. ; if (present(check_column) .and. CS%column_check) check_col = check_column

  I_nlay = 1.0 / (1.0*nlay)
  kg_H_nlay = IG%H_to_kg_m2 * I_nlay

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "write_ice_statistics: Module must be initialized before it is used.")

  nTr_stocks = 0
  if (present(tracer_CSp)) then
    call SIS_call_tracer_stocks(G, IG, tracer_CSp, IST%mH_ice, Tr_stocks, &
                                stock_names=Tr_names, stock_units=Tr_units, num_stocks=nTr_stocks)
  endif

! nTr_stocks = 0
! if (present(tracer_CSp)) then
!   call call_tracer_stocks(h, Tr_stocks, G, tracer_CSp, stock_names=Tr_names, stock_units=Tr_units, num_stocks=nTr_stocks,&
!                              got_min_max=Tr_minmax_got, global_min=Tr_min, global_max=Tr_max, &
!                              xgmin=Tr_min_x, ygmin=Tr_min_y, zgmin=Tr_min_z,&
!                              xgmax=Tr_max_x, ygmax=Tr_max_y, zgmax=Tr_max_z)
!   if (nTr_stocks > 0) then
!     do m=1,nTr_stocks
!       vars(num_nc_fields+m) = &
!         vardesc(Tr_names(m), Tr_names(m),'1','1','s',Tr_units(m))
!     enddo
!     num_nc_fields = num_nc_fields + nTr_stocks
!   endif
! endif

  if (CS%previous_calls == 0) then

    !  Reopen or create a text output file, with an explanatory header line.
    if (is_root_pe()) then
      if (day > CS%Start_time) then
        call open_file(CS%statsfile_ascii, trim(CS%statsfile), &
                       action=APPEND_FILE, form=ASCII_FILE, nohdrs=.true.)
      else
        call open_file(CS%statsfile_ascii, trim(CS%statsfile), &
                       action=WRITEONLY_FILE, form=ASCII_FILE, nohdrs=.true.)
        if (abs(CS%timeunit - 86400.0) < 1.0) then
          write(CS%statsfile_ascii,'("  Step,",7x,"Day,",20x,"Area(N/S),",22x,"Extent(N/S),",17x,&
              &"Mass(N/S),",22x,"Heat(N/S),",14x,"Salinty(N/S),   Frac Mass Err,   Temp Err,   Salin Err")')
          write(CS%statsfile_ascii,'(12x,"[days]",23x,"[m2]",28x,"[m2]",24x,"[kg]",29x,&
              &"[J]",21x,"[g/kg]",10x,"[Nondim]",6x,"[Nondim]",6x,"[Nondim]")')
        else
          if ((CS%timeunit >= 0.99) .and. (CS%timeunit < 1.01)) then
            time_units = "           [seconds]     "
          else if ((CS%timeunit >= 3599.0) .and. (CS%timeunit < 3601.0)) then
            time_units = "            [hours]      "
          else if ((CS%timeunit >= 86399.0) .and. (CS%timeunit < 86401.0)) then
            time_units = "             [days]      "
          else if ((CS%timeunit >= 3.0e7) .and. (CS%timeunit < 3.2e7)) then
            time_units = "            [years]      "
          else
            write(time_units,'(9x,"[",es8.2," s]    ")') CS%timeunit
          endif

          write(CS%statsfile_ascii,'("  Step,",7x,"Time,  Area(N/S),  Extent(N/S),   &
              &Mass, Heat,  Salt,  Frac Mass Err,    Heat Err,   Salin Err")')
          write(CS%statsfile_ascii,'(A25,10x,"[m2]",11x,"[m2]",7x,"[kg]",13x,&
              &"[J]",9x,"[kg]",6x,"[Nondim]",8x,"[J]",8x,"[kg]")') time_units
        endif
      endif
    endif

    statspath_nc = trim(CS%statsfile) // ".nc"
!   if (day > CS%Start_time) then
!     call reopen_file(CS%statsfile_nc, trim(statspath_nc), vars, &
!                      num_nc_fields, G, CS%fields, SINGLE_FILE, CS%timeunit)
!   else
!     call create_file(CS%statsfile_nc, trim(statspath_nc), vars, &
!                      num_nc_fields, G, CS%fields, SINGLE_FILE, CS%timeunit)
!   endif
  endif

! The following quantities are to be written by hemisphere:
!   Ice area, ice extent, Ice+snow mass, enthalpy, salt
! Error analysis on mass, enthalpy, salt

  ice_area(:,:,:) = 0.0
  ice_extent(:,:,:) = 0.0
  col_mass(:,:,:) = 0.0
  col_heat(:,:,:) = 0.0
  col_salt(:,:,:) = 0.0

  enth_liq_0 = Enth_from_TS(0.0, 0.0, IST%ITV)
  do j=js,je ; do i=is,ie
    hem = 1 ; if (G%geolatT(i,j) < 0.0) hem = 2
    do k=1,ncat ; if (G%mask2dT(i,j) * IST%part_size(i,j,k) > 0.0) then
      area_pt = G%areaT(i,j) * G%mask2dT(i,j) * IST%part_size(i,j,k)

      ice_area(i,j,hem) = ice_area(i,j,hem) + area_pt
      col_mass(i,j,hem) = col_mass(i,j,hem) + area_pt * IG%H_to_kg_m2 * &
                          (IST%mH_ice(i,j,k) + (IST%mH_snow(i,j,k) + &
                           IST%mH_pond(i,j,k))) ! mw/new - assumed pond heat/salt = 0

      col_heat(i,j,hem) = col_heat(i,j,hem) + area_pt * IG%H_to_kg_m2 * &
                          (IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1) + &
                           IST%mH_pond(i,j,k) * enth_liq_0)
      do L=1,nlay
        col_heat(i,j,hem) = col_heat(i,j,hem) + area_pt * &
                            ((IST%mH_ice(i,j,k)*kg_H_nlay) * IST%enth_ice(i,j,k,L))
        col_salt(i,j,hem) = col_salt(i,j,hem) + area_pt * &
                  ((0.001*IST%mH_ice(i,j,k)*kg_H_nlay) * IST%sal_ice(i,j,k,L))
      enddo
    endif ; enddo
    if (ice_area(i,j,hem) > 0.1*G%AreaT(i,j)) ice_extent(i,j,hem) = G%AreaT(i,j)

  enddo ; enddo
  Area = reproducing_sum(ice_area, sums=Area_NS)
  Extent = reproducing_sum(ice_extent, sums=Extent_NS)
  Heat = reproducing_sum(col_heat, sums=Heat_NS, EFP_sum=heat_EFP)
  Mass = reproducing_sum(col_mass, sums=Mass_NS, EFP_sum=mass_EFP)
  Salt = reproducing_sum(col_salt, sums=Salt_NS, EFP_sum=salt_EFP)

  salinity_NS(:) = 0.0
  do hem=1,2
    if (mass_NS(hem) > 0.0) salinity_NS(hem) = salt_NS(hem) / mass_NS(hem)
  enddo


! Calculate the maximum CFL numbers.
  max_CFL = 0.0
  dt_CFL = max(CS%dt, 0.)
  if (allocated(IST%u_ice_C)) then ; do j=js,je ; do I=is-1,ie
    if (IST%u_ice_C(I,j) < 0.0) then
      CFL_trans = (-IST%u_ice_C(I,j) * dt_CFL) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
    else
      CFL_trans = (IST%u_ice_C(I,j) * dt_CFL) * (G%dy_Cu(I,j) * G%IareaT(i,j))
    endif
    max_CFL = max(max_CFL, CFL_trans)
  enddo ; enddo ; endif
  if (allocated(IST%v_ice_C)) then ; do J=js-1,je ; do i=is,ie
    if (IST%v_ice_C(i,J) < 0.0) then
      CFL_trans = (-IST%v_ice_C(i,J) * dt_CFL) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
    else
      CFL_trans = (IST%v_ice_C(i,J) * dt_CFL) * (G%dx_Cv(i,J) * G%IareaT(i,j))
    endif
    max_CFL = max(max_CFL, CFL_trans)
  enddo ; enddo ; endif
  if ( .not.(allocated(IST%u_ice_C) .or. allocated(IST%v_ice_C)) .and. &
       (allocated(IST%u_ice_B) .and. allocated(IST%v_ice_B)) ) then
    do J=js-1,je ; do I=is-1,ie
      CFL_u = abs(IST%u_ice_B(I,J)) * dt_CFL * G%IdxBu(I,J)
      CFL_v = abs(IST%v_ice_B(I,J)) * dt_CFL * G%IdyBu(I,J)
      max_CFL = max(max_CFL, CFL_u, CFL_v)
    enddo ; enddo
  endif

  call sum_across_PEs(CS%ntrunc)
  if (nTr_stocks > 0) call sum_across_PEs(Tr_stocks,nTr_stocks)
  call max_across_PEs(max_CFL)

  if (CS%previous_calls == 0) then
    CS%mass_prev = Mass ; CS%fresh_water_input = 0.0
    CS%salt_prev = Salt ; CS%net_salt_input = 0.0
    CS%heat_prev = Heat ; CS%net_heat_input = 0.0

    CS%mass_prev_EFP = mass_EFP ; CS%fresh_water_in_EFP = real_to_EFP(0.0)
    CS%salt_prev_EFP = salt_EFP ; CS%net_salt_in_EFP = real_to_EFP(0.0)
    CS%heat_prev_EFP = heat_EFP ; CS%net_heat_in_EFP = real_to_EFP(0.0)
  else
    do j=js,je ; do i=is,ie
      area_h = G%areaT(i,j) * G%mask2dT(i,j)
      CS%water_in_col(i,j) = area_h * CS%water_in_col(i,j)
      CS%heat_in_col(i,j) = area_h * CS%heat_in_col(i,j)
      CS%salt_in_col(i,j) = area_h * CS%salt_in_col(i,j)
    enddo ; enddo

    CS%fresh_water_input = reproducing_sum(CS%water_in_col, EFP_sum=CS%fresh_water_in_EFP)
    CS%net_salt_input = reproducing_sum(CS%salt_in_col, EFP_sum=CS%net_salt_in_EFP)
    CS%net_heat_input = reproducing_sum(CS%heat_in_col, EFP_sum=CS%net_heat_in_EFP)
  endif

  Salt_chg_EFP = Salt_EFP - CS%salt_prev_EFP
  Salt_anom_EFP = Salt_chg_EFP - CS%net_salt_in_EFP
  Salt_chg = EFP_to_real(Salt_chg_EFP) ; Salt_anom = EFP_to_real(Salt_anom_EFP)
  Heat_chg_EFP = Heat_EFP - CS%heat_prev_EFP
  Heat_anom_EFP = Heat_chg_EFP - CS%net_heat_in_EFP
  Heat_chg = EFP_to_real(Heat_chg_EFP) ; Heat_anom = EFP_to_real(Heat_anom_EFP)

  mass_chg_EFP = mass_EFP - CS%mass_prev_EFP
  salin_mass_in = 0.0
!  if (G%Boussinesq) then
    mass_anom_EFP = mass_chg_EFP - CS%fresh_water_in_EFP
!  else
    ! net_salt_input needs to be converted from psu m s-1 to kg m-2 s-1.
!    mass_anom_EFP = mass_chg_EFP - CS%fresh_water_in_EFP
!    salin_mass_in = 0.001*EFP_to_real(CS%net_salt_in_EFP)
!  endif
  mass_chg = EFP_to_real(mass_chg_EFP)
  mass_anom = EFP_to_real(mass_anom_EFP) - salin_mass_in

  I_Mass = 0.0 ; if (Mass > 0.0) I_Mass = 1.0/Mass
  salin = Salt * I_Mass ; salin_anom = Salt_anom * I_Mass

 ! salin_chg = Salt_chg / Mass
!  temp = heat / (Mass* CI) ; temp_anom = Heat_anom / (Mass* CI)

  call get_time(day, start_of_day, num_days)
  if (abs(CS%timeunit - 86400.0) < 1.0) then
    reday = REAL(num_days)+ (REAL(start_of_day)/86400.0)
    mesg_intro = "SIS Day "
  else
    reday = REAL(num_days)*(86400.0/CS%timeunit) + &
            REAL(start_of_day)/abs(CS%timeunit)
    mesg_intro = "SIS Time "
  endif
  if (reday < 1.0e8) then ;      write(day_str, '(F12.3)') reday
  elseif (reday < 1.0e11) then ; write(day_str, '(F15.3)') reday
  else ;                         write(day_str, '(ES15.9)') reday ; endif

  if     (n < 1000000)   then ; write(n_str, '(I6)')  n
  elseif (n < 10000000)  then ; write(n_str, '(I7)')  n
  elseif (n < 100000000) then ; write(n_str, '(I8)')  n
  else                        ; write(n_str, '(I10)') n ; endif

  if     (CS%ntrunc < 1000000)   then ; write(trunc_str, '(I6)')  CS%ntrunc
  elseif (CS%ntrunc < 10000000)  then ; write(trunc_str, '(I7)')  CS%ntrunc
  elseif (CS%ntrunc < 100000000) then ; write(trunc_str, '(I8)')  CS%ntrunc
  else                                ; write(trunc_str, '(I10)') CS%ntrunc ; endif

  msg_start = trim(n_str)//","//trim(day_str)
  if (present(message)) msg_start = trim(message)
  msg_start = trim(msg_start)//", "//trim(trunc_str)

  if (is_root_pe()) then
    Heat_anom_norm = 0.0 ; if (Heat /= 0.0) Heat_anom_norm = Heat_anom/Heat
    Salt_anom_norm = 0.0 ; if (Salt /= 0.0) Salt_anom_norm = Salt_anom/Salt
    write(CS%statsfile_ascii,'(A,", Area", 2(ES19.12), ", Ext", 2(es11.4), ", CFL", F6.3, &
                &", M",2(ES12.5),", Enth",2(ES13.5),", S ",2(f8.4),", Me ",ES9.2,&
                &", Te ",ES9.2,", Se ",ES9.2)') &
          trim(msg_start), Area_NS(1:2), Extent_NS(1:2), max_CFL, mass_NS(1:2), &
          heat_NS(1:2), 1000.*salinity_NS(1:2), mass_anom * I_Mass, &
          Heat_anom_norm, salt_anom_norm
  endif

  if (is_root_pe() .and. CS%write_stdout) then
    if (get_calendar_type() == NO_CALENDAR) then
      write(*,'(A,A," ",A,": Area", 2(ES19.12), ", Mass ", 2(ES18.11))') &
        trim(mesg_intro), trim(day_str(1:3))//trim(day_str(4:)), trim(n_str), &
        Area_NS(1:2), mass_NS(1:2)
    else
      call get_date(day, iyear, imonth, iday, ihour, iminute, isecond, itick)
      write(*,'("SIS Date",i7,2("/",i2.2)," ",i2.2,2(":",i2.2)," ",A, &
        &": Area", 2(ES19.12), ", Mass ", 2(ES18.11))') &
        iyear, imonth, iday, ihour, iminute, isecond, trim(n_str), Area_NS(1:2), mass_NS(1:2)
    endif

    if (CS%ntrunc > 0) then
      write(*,'(A," Sea Ice Truncations ",I0)') &
        trim(mesg_intro)//trim(day_str), CS%ntrunc
    endif

    if (CS%write_stocks) then
      msg_start = " Total"
      if (present(message)) msg_start = trim(message)
      write(*,'(A," Ice Mass: ",ES24.16,", Change: ",ES12.5," Error: ",ES12.5," (",ES8.1,")")') &
            trim(msg_start), Mass, mass_chg, mass_anom, mass_anom * I_Mass
      if (Salt == 0.) then
        write(*,'(A," Ice Salt: ",ES24.16,", Change: ",ES12.5," Error: ",ES12.5)') &
            trim(msg_start), Salt*0.001, Salt_chg*0.001, Salt_anom*0.001
      else
        write(*,'(A," Ice Salt: ",ES24.16,", Change: ",ES12.5," Error: ",ES12.5," (",ES8.1,")")') &
            trim(msg_start), Salt*0.001, Salt_chg*0.001, Salt_anom*0.001, Salt_anom/Salt
      endif
      if (Heat == 0.) then
        write(*,'(A," Ice Heat: ",ES24.16,", Change: ",ES12.5," Error: ",ES12.5)') &
            trim(msg_start), Heat, Heat_chg, Heat_anom
      else
        write(*,'(A," Ice Heat: ",ES24.16,", Change: ",ES12.5," Error: ",ES12.5," (",ES8.1,")")') &
            trim(msg_start), Heat, Heat_chg, Heat_anom, Heat_anom/Heat
      endif

      if (present(tracer_CSp)) then
        do m=1,nTr_stocks
          write(*,'("      Total ",a,": ",ES24.16,X,a)') &
             trim(Tr_names(m)), Tr_stocks(m), trim(Tr_units(m))
        enddo
      endif

!     do m=1,nTr_stocks

!        write(*,'("      Total ",a,": ",ES24.16,X,a)') &
!             trim(Tr_names(m)), Tr_stocks(m), trim(Tr_units(m))
!
!        if(Tr_minmax_got(m)) then
!          write(*,'(64X,"Global Min:",ES24.16,X,"at: (", f7.2,","f7.2,","f8.2,")"  )') &
!               Tr_min(m),Tr_min_x(m),Tr_min_y(m),Tr_min_z(m)
!          write(*,'(64X,"Global Max:",ES24.16,X,"at: (", f7.2,","f7.2,","f8.2,")"  )') &
!               Tr_max(m),Tr_max_x(m),Tr_max_y(m),Tr_max_z(m)
!       endif

!     enddo
    endif ! write_stocks
  endif ! write_stdout

  if (check_col .and. (CS%previous_calls > 0)) then ; do j=js,je ; do i=is,ie
    hem = 1 ; if (G%geolatT(i,j) < 0.0) hem = 2
    heat_imb = (col_heat(i,j,hem) - CS%heat_col_prev(i,j)) - CS%heat_in_col(i,j)
    mass_imb = (col_mass(i,j,hem) - CS%water_col_prev(i,j)) - CS%water_in_col(i,j)
    if (abs(mass_imb) > CS%imb_tol*abs(Mass) .and. (abs(Mass) > 0.0)) then
      write(mesg,'("Mass imbalance of ",ES11.4," (",ES8.1,") detected at i,j=",2(i4), &
                  &" Lon/Lat = ",2(f8.2))') &
                  mass_imb, mass_imb/max(abs(mass),abs(mass_imb)), &
                  i, j, G%geolonT(i,j), G%geolatT(i,j)
      call SIS_error(WARNING, mesg, all_print=.true.)
    endif
    if (abs(heat_imb) > CS%imb_tol*abs(Heat) .and. (abs(Heat) > 0.0)) then
      write(mesg,'("Heat imbalance of ",ES11.4," (",ES8.1,") detected at i,j=",2(i4), &
                  &" Lon/Lat = ",2(f8.2))') &
                  heat_imb, heat_imb/max(abs(heat),abs(heat_imb)), i, j, &
                  G%geolonT(i,j), G%geolatT(i,j)
      call SIS_error(WARNING, mesg, all_print=.true.)
    endif
  enddo ; enddo ; endif

! var = real(CS%ntrunc)
! call write_field(CS%statsfile_nc, CS%fields(1), var, reday)
! call write_field(CS%statsfile_nc, CS%fields(2), toten, reday)
! call write_field(CS%statsfile_nc, CS%fields(3), PE, reday)
! call write_field(CS%statsfile_nc, CS%fields(4), KE, reday)
! call write_field(CS%statsfile_nc, CS%fields(5), H_0APE, reday)
! call write_field(CS%statsfile_nc, CS%fields(6), mass_lay, reday)

! call write_field(CS%statsfile_nc, CS%fields(7), Mass, reday)
! call write_field(CS%statsfile_nc, CS%fields(8), mass_chg, reday)
! call write_field(CS%statsfile_nc, CS%fields(9), mass_anom, reday)
! call write_field(CS%statsfile_nc, CS%fields(10), max_CFL(1), reday)
! call write_field(CS%statsfile_nc, CS%fields(11), max_CFL(1), reday)

! call write_field(CS%statsfile_nc, CS%fields(12), 0.001*Salt, reday)
! call write_field(CS%statsfile_nc, CS%fields(13), 0.001*salt_chg, reday)
! call write_field(CS%statsfile_nc, CS%fields(14), 0.001*salt_anom, reday)
! call write_field(CS%statsfile_nc, CS%fields(15), Heat, reday)
! call write_field(CS%statsfile_nc, CS%fields(16), heat_chg, reday)
! call write_field(CS%statsfile_nc, CS%fields(17), heat_anom, reday)
! do m=1,nTr_stocks
!   call write_field(CS%statsfile_nc, CS%fields(17+m), Tr_stocks(m), reday)
! enddo

! call flush_file(CS%statsfile_nc)

  if (is_root_pe() .and. (CS%ntrunc>CS%maxtrunc)) then
    call SIS_error(FATAL, "write_ice_statistics: Sea ice velocity has been "//&
                          "truncated too many times.")
  endif
  CS%ntrunc = 0
  CS%previous_calls = CS%previous_calls + 1
  if (CS%column_check) then ; do j=js,je ; do i=is,ie
    CS%water_col_prev(i,j) = col_mass(i,j,1) + col_mass(i,j,2)
    CS%heat_col_prev(i,j) = col_heat(i,j,1) + col_heat(i,j,2)
    CS%salt_col_prev(i,j) = col_salt(i,j,1) + col_salt(i,j,2)
  enddo ; enddo ; endif
  CS%mass_prev = Mass ; CS%fresh_water_input = 0.0
  CS%salt_prev = Salt ; CS%net_salt_input = 0.0
  CS%heat_prev = Heat ; CS%net_heat_input = 0.0

  CS%water_in_col(:,:) = 0.0
  CS%heat_in_col(:,:) = 0.0
  CS%salt_in_col(:,:) = 0.0

  CS%mass_prev_EFP = mass_EFP ; CS%fresh_water_in_EFP = real_to_EFP(0.0)
  CS%salt_prev_EFP = Salt_EFP ; CS%net_salt_in_EFP = real_to_EFP(0.0)
  CS%heat_prev_EFP = Heat_EFP ; CS%net_heat_in_EFP = real_to_EFP(0.0)
end subroutine write_ice_statistics


subroutine accumulate_bottom_input(IST, OSS, FIA, IOF, dt, G, IG, CS)
!   This subroutine accumulates the net input of fresh water and heat through
! the bottom of the sea-ice for conservation checks.
! Arguments: IST - The internal sea ice state type.
!  (in)      dt - The amount of time over which to average.
!  (in)      G - The sea ice model's grid structure.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 SIS_sum_output_init.
  type(SIS_hor_grid_type),    intent(in) :: G
  type(ice_grid_type),        intent(in) :: IG
  type(ice_state_type),       intent(in) :: IST
  type(ocean_sfc_state_type), intent(in) :: OSS
  type(fast_ice_avg_type),    intent(in) :: FIA
  type(ice_ocean_flux_type),  intent(in) :: IOF
  real,                       intent(in) :: dt
  type(SIS_sum_out_CS),       pointer    :: CS

  real :: Flux_SW, enth_units, LI

  integer :: i, j, k, isc, iec, jsc, jec, ncat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  call get_SIS2_thermo_coefs(IST%ITV, enthalpy_units=enth_units, Latent_fusion=LI)

  if (CS%dt < 0.0) CS%dt = dt

  do j=jsc,jec ; do i=isc,iec
    CS%water_in_col(i,j) = CS%water_in_col(i,j) - dt * &
           ( ((FIA%runoff(i,j) + FIA%calving(i,j)) + &
              (IOF%lprec_ocn_top(i,j) + IOF%fprec_ocn_top(i,j))) - IOF%flux_q_ocn_top(i,j) )
    Flux_SW = (IOF%flux_sw_vis_dir_ocn(i,j) + IOF%flux_sw_vis_dif_ocn(i,j)) + &
              (IOF%flux_sw_nir_dir_ocn(i,j) + IOF%flux_sw_nir_dif_ocn(i,j))
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) - (dt * enth_units) * &
          ( Flux_SW + &
           ((IOF%flux_lw_ocn_top(i,j) - IOF%flux_lh_ocn_top(i,j)) - IOF%flux_t_ocn_top(i,j)) + &
            (-LI)*(IOF%fprec_ocn_top(i,j) + FIA%calving(i,j)) )
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) - enth_units * &
           (OSS%frazil(i,j)-FIA%frazil_left(i,j))
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + &
           ((IOF%Enth_Mass_in_atm(i,j) + IOF%Enth_Mass_in_ocn(i,j)) + &
            (IOF%Enth_Mass_out_atm(i,j) + IOF%Enth_Mass_out_ocn(i,j)) )
    CS%salt_in_col(i,j) = CS%salt_in_col(i,j) + dt * IOF%flux_salt(i,j)
  enddo ; enddo

end subroutine accumulate_bottom_input

subroutine accumulate_input_1(IST, FIA, dt, G, IG, CS)
!   This subroutine accumulates the net input of fresh water and heat through
! the top of the sea-ice for conservation checks.

! Arguments: IST - The internal sea ice state type.
!  (in)      dt - The amount of time over which to average.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      G - The sea ice model's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 SIS_sum_output_init.
  type(ice_state_type),    intent(in) :: IST
  type(fast_ice_avg_type), intent(in) :: FIA
  real,                    intent(in) :: dt
  type(SIS_hor_grid_type), intent(in) :: G
  type(ice_grid_type),     intent(in) :: IG
  type(SIS_sum_out_CS),    pointer    :: CS

  real, dimension(SZI_(G),SZJ_(G)) :: &
    FW_in, &   ! The net fresh water input, integrated over a timestep in kg.
    salt_in, & ! The total salt added by surface fluxes, integrated
               ! over a time step in PSU kg.
    heat_in    ! The total heat added by surface fluxes, integrated
               ! over a time step in Joules.
  real :: FW_input   ! The net fresh water input, integrated over a timestep
                  ! and summed over space, in kg.
  real :: salt_input ! The total salt added by surface fluxes, integrated
                  ! over a time step and summed over space, in kg.
  real :: heat_input ! The total heat added by surface fluxes, integrated
                  ! over a time step and summed over space, in Joules.
  real :: area_h, area_pt, Flux_SW
  real :: enth_units
  type(EFP_type) :: &
    FW_in_EFP, &   ! Extended fixed point versions of FW_input, salt_input, and
    salt_in_EFP, & ! heat_input, in kg, PSU kg, and Joules.
    heat_in_EFP    !
  integer :: i, j, k, isc, iec, jsc, jec, ncat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  call get_SIS2_thermo_coefs(IST%ITV, enthalpy_units=enth_units)

  FW_in(:,:) = 0.0 ; salt_in(:,:) = 0.0 ; heat_in(:,:) = 0.0

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,IST,CS,enth_units,dt,FIA) &
!$OMP                          private(area_pt,Flux_SW)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    area_pt = IST%part_size(i,j,k)
    Flux_SW = (FIA%flux_sw_vis_dir_top(i,j,k) + FIA%flux_sw_vis_dif_top(i,j,k)) + &
              (FIA%flux_sw_nir_dir_top(i,j,k) + FIA%flux_sw_nir_dif_top(i,j,k))
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + ((dt * area_pt) * enth_units) * &
        ( Flux_SW * (1.0 - FIA%sw_abs_ocn(i,j,k)) + &
          ((FIA%flux_lw_top(i,j,k) - FIA%flux_t_top(i,j,k)) )  + &
           (-FIA%flux_lh_top(i,j,k)) + FIA%bheat(i,j))
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) - (enth_units * area_pt) * &
                   (FIA%bmelt(i,j,k) + FIA%tmelt(i,j,k))
  enddo ; enddo ; enddo

end subroutine accumulate_input_1

subroutine accumulate_input_2(IST, FIA, IOF, part_size, dt, G, IG, CS)
!   This subroutine accumulates the net input of fresh water and heat through
! the top of the sea-ice for conservation checks.

! Arguments: IST - The internal sea ice state type.
!  (in)      part_size - The fractional ice concentration within a cell in each
!                      thickness category, nondimensional, 0-1.
!  (in)      dt - The amount of time over which to average.
!  (in)      G - The sea ice model's grid structure.
!  (in)      IG - The sea-ice-specific grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 SIS_sum_output_init.
  type(SIS_hor_grid_type), intent(inout) :: G
  type(ice_grid_type),     intent(inout) :: IG
  type(ice_state_type),    intent(inout) :: IST
  type(fast_ice_avg_type),    intent(in) :: FIA
  type(ice_ocean_flux_type),  intent(in) :: IOF
  real, dimension(SZI_(G),SZJ_(G),SZCAT0_(IG)), intent(in) :: part_size
  real,                    intent(in) :: dt
  type(SIS_sum_out_CS),    pointer    :: CS

  real :: area_pt, Flux_SW, pen_frac
  real :: enth_units, LI
  integer :: i, j, k, m, isc, iec, jsc, jec, ncat

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  ! This subroutine includes the accumulation of mass fluxes and heat fluxes
  ! into the ice that are known before SIS#_thermodynamics, as well the
  ! ice-top fluxes that will be passed on directly to the ocean.  It does
  ! not include the enthalpy changes due to net mass changes in the ice,
  ! as these are not yet known.

  call get_SIS2_thermo_coefs(IST%ITV, enthalpy_units=enth_units, Latent_fusion=LI)
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,CS,dt,IST,FIA,IOF,&
!$OMP                                  enth_units, LI) &
!$OMP                          private(area_pt)
  do j=jsc,jec ; do i=isc,iec
    ! Runoff and calving are passed directly on to the ocean.
    CS%water_in_col(i,j) = CS%water_in_col(i,j) + dt * &
          (FIA%runoff(i,j) + FIA%calving(i,j))

    area_pt = IST%part_size(i,j,0)
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + ((dt * area_pt) * enth_units) * &
          ((FIA%flux_lw_top(i,j,0) - FIA%flux_lh_top(i,j,0)) - FIA%flux_t_top(i,j,0))

    ! These are mass fluxes that are simply passed through to the ocean.
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + (dt * enth_units) * (-LI) * &
                      (area_pt * FIA%fprec_top(i,j,0) + FIA%calving(i,j))

  enddo ; enddo

  ! The terms that are added here include surface fluxes that will be passed
  ! directly on into the ocean.
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ncat,part_size,IST,CS,dt,enth_units,FIA)&
!$OMP                          private(area_pt,pen_frac,Flux_SW)
    do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
      area_pt = part_size(i,j,k)
      pen_frac = 1.0 ; if (k>0) pen_frac = FIA%sw_abs_ocn(i,j,k)
      Flux_SW = (FIA%flux_sw_vis_dir_top(i,j,k) + FIA%flux_sw_vis_dif_top(i,j,k)) + &
                (FIA%flux_sw_nir_dir_top(i,j,k) + FIA%flux_sw_nir_dif_top(i,j,k))

      CS%water_in_col(i,j) = CS%water_in_col(i,j) + (dt * area_pt) * &
          ( (FIA%lprec_top(i,j,k) + FIA%fprec_top(i,j,k)) - FIA%flux_q_top(i,j,k) )
      CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + ((dt * area_pt) * enth_units) * &
           ( pen_frac*Flux_SW )

      if (k>0) &
        CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + (area_pt * enth_units) * &
           ((FIA%bmelt(i,j,k) + FIA%tmelt(i,j,k)) - dt*FIA%bheat(i,j))
    enddo ; enddo ; enddo

 ! Runoff and calving do not bring in salt, so salt_in(i,j) = 0.0

end subroutine accumulate_input_2

end module SIS_sum_output
