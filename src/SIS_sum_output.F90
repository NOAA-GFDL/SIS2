!>  Subroutines that calculate globally integrated sea-ice quantities for SIS2, and writes them
!! to a netcdf file and an ASCII output file.
module SIS_sum_output

! This file is a part of SIS2.  See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!   This file contains the subroutines that calculate globally integrated      !
! sea-ice quantities for SIS2, and writes them to a netcdf file and an ASCII   !
! output file.  This code was originally adapted from MOM_sum_output.F90       !
! by Robert Hallberg in May 2014.                                              !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_coms,          only : EFP_type, operator(+), operator(-), assignment(=)
use MOM_coms,          only : reproducing_sum, reproducing_sum_EFP, EFP_to_real, real_to_EFP
use MOM_coms,          only : EFP_sum_across_PEs, max_across_PEs
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
! use MOM_io,          only : create_file, fieldtype, flush_file, reopen_file, vardesc, write_field
use MOM_io,            only : open_ASCII_file, APPEND_FILE, ASCII_FILE, SINGLE_FILE, WRITEONLY_FILE
use MOM_string_functions, only : slasher
use MOM_time_manager,  only : time_type, get_time, operator(>), operator(-)
use MOM_time_manager,  only : get_date, get_calendar_type, NO_CALENDAR
use MOM_unit_scaling,  only : unit_scale_type
use SIS_types,         only : ice_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS_types,         only : ocean_sfc_state_type
use SIS_hor_grid,      only : SIS_hor_grid_type
use ice_grid,          only : ice_grid_type
use SIS2_ice_thm,      only : enthalpy_liquid_freeze, get_SIS2_thermo_coefs, ice_thermo_type
use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS, SIS_call_tracer_stocks

implicit none ; private

#include <SIS2_memory.h>

public write_ice_statistics, accumulate_bottom_input
public SIS_sum_output_init, SIS_sum_output_end
public accumulate_input_1, accumulate_input_2

!-----------------------------------------------------------------------

! integer, parameter :: NUM_FIELDS = 17

!> This structure contains the parameters that regulate the summed output.
type, public :: SIS_sum_out_CS ; private
  real, dimension(:,:), allocatable :: &
    water_in_col, &             !< The water that has been input to the ice and snow in a column since
                                !! the last time that write_ice_statistics was called [R Z ~> kg m-2]
    heat_in_col, &              !< The heat that has been input to the ice and snow in a column since
                                !! the last time that write_ice_statistics was called [Q R Z ~> J m-2]
    salt_in_col, &              !< The salt that has been input to the ice and snow in a column since
                                !! the last time that write_ice_statistics was called [1e3 S R Z ~> kgSalt m-2]
    ! These three arrays are only allocated and used for monitoring column-wise conservation.
    water_cell_prev, &          !< The cell integrated water that was in the ice and snow the last
                                !! time that write_ice_statistics was called [kg].
    heat_cell_prev, &           !< The cell integrated heat that was in the ice and snow the last
                                !! time that write_ice_statistics was called [J].
    salt_cell_prev              !< The cell integrated salt that was in the ice and snow the last
                                !! time that write_ice_statistics was called [kgSalt].

  type(EFP_type) :: heat_prev_EFP !<   The total amount of heat in the sea ice the last
                                !! time that write_ice_statistics was called [J], in EFP form.
  type(EFP_type) :: salt_prev_EFP !<   The total amount of salt in the sea ice the last
                                !! time that write_ice_statistics was called [kgSalt], in EFP form.
  type(EFP_type) :: mass_prev_EFP !<   The total sea ice mass the last time that
                                !! write_ice_statistics was called [kg], in EFP form.
  real    :: dt                 !< The ice dynamics time step [T ~> s].
  real    :: timeunit           !<   The length of the units for the time axis [s].
  type(time_type) :: Start_time !< The start time of the simulation.
                                !< Start_time is set in SIS_initialization.F90
  logical :: column_check       !< If true, enable the column by column heat and
                                !! mass conservation check
  real    :: imb_tol            !< The tolerance for imbalances to be flagged by
                                !! column_check [nondim].
  integer :: maxtrunc           !< The number of truncations per ice statistics
                                !! save interval at which the run is stopped.
  logical :: write_stdout       !< If true, periodically write sea ice statistics
                                !! to stdout to allow the progress to be seen.
  logical :: write_stocks       !< If true, write the integrated tracer amounts
                                !! to stdout when the statistics files are written.
  integer :: previous_calls = 0 !< The number of times write_ice_statistics has been called.
  integer :: prev_n = 0         !< The value of n from the last call.
  integer, pointer :: ntrunc => NULL() !< The number of times the velocity has been truncated
                                !! since the last call to write_ice_statistics.
! integer :: statsfile_nc       !< NetCDF id of the statistics file.
  integer :: statsfile_ascii    !< The unit number of the ascii version of the statistics file.
! type(fieldtype), dimension(NUM_FIELDS+MAX_FIELDS_) :: &
!            fields             !< fieldtype variables for the output fields.
  character(len=200) :: statsfile  !< The name of the statistics file with path.

end type SIS_sum_out_CS

contains

!> Initialize the SIS_sum_output control structure, allocate memory and store runtime parameters.
subroutine SIS_sum_output_init(G, param_file, directory, Input_start_time, US, CS, ntrunc)
  type(SIS_hor_grid_type),  intent(in)    :: G      !< The horizontal grid type
  type(param_file_type),    intent(in)    :: param_file !< A structure to parse for run-time parameters
  character(len=*),         intent(in)    :: directory  !<  The directory where the statistics file goes
  type(time_type),          intent(in)    :: Input_start_time !< The start time of the simulation
  type(unit_scale_type),    intent(in)    :: US     !< A structure with unit conversion factors
  type(SIS_sum_out_CS),     pointer       :: CS     !< A pointer that is set to point to the control
                                                    !! structure for this module
  integer, target, optional,intent(inout) :: ntrunc !< The integer that stores the number of times
                                                    !! the velocity has been truncated since the
                                                    !! last call to write_ice_statistics

  ! Local variables
  character(len=40)  :: mdl = "SIS_sum_output" ! This module's name.
  character(len=200) :: statsfile  ! The name of the statistics file.
  ! This include declares and sets the variable "version".
# include "version_variable.h"

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_sum_output_init called with associated control structure.")
    return
  endif
  allocate(CS)

  if (present(ntrunc)) then ; CS%ntrunc => ntrunc ; else ; allocate(CS%ntrunc) ; endif
  CS%ntrunc = 0

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "WRITE_STOCKS", CS%write_stocks, &
                 "If true, write the integrated tracer amounts to stdout "//&
                 "when the statistics files are written.", default=.true.)
  call get_param(param_file, mdl, "STDOUT_HEARTBEAT", CS%write_stdout, &
                 "If true, periodically write sea ice statistics to "//&
                 "stdout to allow the progress to be seen.", default=.true.)
  call get_param(param_file, mdl, "DT_ICE_DYNAMICS", CS%dt, &
                 "The time step used for the slow ice dynamics, including "//&
                 "stepping the continuity equation and interactions between "//&
                 "the ice mass field and velocities.", units="s", scale=US%s_to_T, &
                 default=-1.0, do_not_log=.true.)
  call get_param(param_file, mdl, "MAXTRUNC", CS%maxtrunc, &
                 "The run will be stopped, and the day set to a very "//&
                 "large value if the velocity is truncated more than "//&
                 "MAXTRUNC times between  writing ice statistics. "//&
                 "Set MAXTRUNC to 0 to stop if there is any truncation "//&
                 "of sea ice velocities.", units="truncations save_interval-1", default=0)

  call get_param(param_file, mdl, "STATISTICS_FILE", statsfile, &
                 "The file to use to write the globally integrated "//&
                 "statistics.", default="seaice.stats")

  CS%statsfile = trim(slasher(directory))//trim(statsfile)
  call log_param(param_file, mdl, "output_path/STATISTICS_FILE", CS%statsfile)
#ifdef STATSLABEL
  CS%statsfile = trim(CS%statsfile)//"."//trim(adjustl(STATSLABEL))
#endif

  call get_param(param_file, mdl, "TIMEUNIT", CS%Timeunit, &
                 "The time unit in seconds a number of input fields", &
                 units="s", default=86400.0)
  if (CS%Timeunit < 0.0) CS%Timeunit = 86400.0
  call get_param(param_file, mdl, "COLUMN_CHECK", CS%column_check, &
                 "If true, add code to allow debugging of conservation "//&
                 "column-by-column.  This does not change answers, but "//&
                 "can increase model run time.", default=.false., &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "IMBALANCE_TOLERANCE", CS%imb_tol, &
                 "The tolerance for imbalances to be flagged by COLUMN_CHECK.", &
                 units="nondim", default=1.0e-9, debuggingParam=.true.)

  CS%Start_time = Input_start_time

  allocate(CS%water_in_col(G%isd:G%ied, G%jsd:G%jed), source=0.0)
  allocate(CS%heat_in_col(G%isd:G%ied, G%jsd:G%jed), source=0.0)
  allocate(CS%salt_in_col(G%isd:G%ied, G%jsd:G%jed), source=0.0)
  if (CS%column_check) then
    allocate(CS%water_cell_prev(G%isd:G%ied, G%jsd:G%jed), source=0.0)
    allocate(CS%heat_cell_prev(G%isd:G%ied, G%jsd:G%jed), source=0.0)
    allocate(CS%salt_cell_prev(G%isd:G%ied, G%jsd:G%jed), source=0.0)
  endif

end subroutine SIS_sum_output_init

!> Deallocate memory associated with the SIS_sum_out control structure
subroutine SIS_sum_output_end(CS)
  type(SIS_sum_out_CS), pointer :: CS !< The control structure returned by a previous call to
                                      !! SIS_sum_output_init that is deallocated here
!   This subroutine deallocates the memory owned by this module.

  if (associated(CS)) deallocate(CS)

end subroutine SIS_sum_output_end

!> Write out the sea ice statistics of the total sea-ice mass, heat and salt by
!! hemisphere and other globally integrated quantities.
subroutine write_ice_statistics(IST, day, n, G, US, IG, CS, message, check_column, tracer_CSp)
  type(ice_state_type),       intent(inout) :: IST !< A type describing the state of the sea ice
  type(time_type),            intent(inout) :: day !< The current model time.
  integer,                    intent(in)    :: n   !< The time step number of the current execution
  type(SIS_hor_grid_type),    intent(inout) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in)    :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(inout) :: IG  !< The sea-ice specific grid type
  type(SIS_sum_out_CS),       pointer       :: CS  !< The control structure returned by a previous
                                                   !! call to SIS_sum_output_init
  character(len=*), optional, intent(in) :: message !< A text message to use with this output
  logical,          optional, intent(in) :: check_column !< If true, check for column-wise heat and
                                                   !! mass conservation.
  type(SIS_tracer_flow_control_CS), &
                    optional, pointer    :: tracer_CSp !< A pointer to the tracer package flow
                                                   !! control module to enable the writing of
                                                   !! the passive tracer statistics.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G), 2) :: &
    ice_area, &    ! The area of ice in each cell and hemisphere [m2].
    ice_extent, &  ! The extent (cells with >10% coverage) of ice in each cell and hemisphere [m2].
    col_mass, &    ! The column integrated ice and snow mass in each cell and hemisphere [kg].
    col_heat, &    ! The column integrated ice and snow heat in each cell and hemisphere [J].
    col_salt       ! The column integrated salt in the ice in each cell and hemisphere [kgSalt].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    water_into_cell, & ! The area integral of the water that has been input to the ice and snow
                   ! in a grid cell since the last time that write_ice_statistics was called [kg]
    heat_into_cell, &  ! The area integral of the heat that has been input to the ice and snow
                   ! in a grid cell since the last time that write_ice_statistics was called [J]
    salt_into_cell ! The area integral of the salt that has been input to the ice and snow
                   ! in a grid cell since the last time that write_ice_statistics was called [kgSalt]

  real, dimension(2) :: &
    Area_NS, &     ! The total sea-ice area in the two hemispheres [m2].
    Extent_NS, &   ! The total sea-ice extent in the two hemispheres [m2].
    heat_NS, &     ! The total sea-ice enthalpy in the two hemispheres [J].
    mass_NS, &     ! The total sea-ice mass in the two hemispheres [kg].
    salt_NS, &     ! The total sea-ice salt in the two hemispheres [kgSalt].
    salinity_NS    ! The average sea-ice salinity in the two hemispheres [gSalt kg-1].

  real :: Mass         ! The total mass of the sea ice and snow atop it [kg].
  real :: mass_chg     ! The change in total sea ice mass of fresh water since
                       ! the last call to this subroutine [kg].
  real :: mass_anom    ! The change in fresh water that cannot be accounted for
                       ! by the surface fluxes [kg].
  real :: I_Mass       ! Adcroft's rule reciprocal of mass: 1/Mass or 0 [kg-1].
  real :: Salt         ! The total amount of salt in the brine pockets in sea ice [kgSalt].
  real :: Salt_chg     ! The change in total sea ice salt since the last call
                       ! to this subroutine [kgSalt].
  real :: Salt_anom    ! The change in salt that cannot be accounted for by the surface fluxes [kgSalt].
  real :: Salt_anom_norm ! The salt anomaly normalized by salt (if it is nonzero) [nondim].
  real :: Heat         ! The total amount of enthalpy in the sea ice, brine pockets and melt ponds [J].
  real :: Heat_chg     ! The change in total sea ice heat since the last call to this subroutine [J].
  real :: Heat_anom    ! The change in heat that cannot be accounted for bythe surface fluxes [J].
  real :: Heat_anom_norm ! The heat anomaly normalized by heat (if it is nonzero) [nondim].
  real :: heat_imb     ! The column integrated heat imbalance [J].
  real :: mass_imb     ! The column integrated mass imbalance [kg].
  real :: enth_liq_0   ! The enthalpy of liquid water at the freezing point [Q ~> J kg-1].
  real :: I_nlay       ! The inverse of the number of layers [nondim]
  real :: mass_scale   ! A mass unit conversion factor for area-integrated mass [kg R-1 Z-1 L-2 ~> 1]
  real :: heat_scale   ! A mass unit conversion factor for area-integrated heat [J Q-1 R-1 Z-1 L-2 ~> 1]
  real :: salt_scale   ! A mass unit conversion factor for area-integrated salt [gSalt S-1 R-1 Z-1 L-2 ~> 1]
  real :: area_scale   ! A mass unit conversion factor for cell area [m2 L-2 ~> 1]
  real :: area_pt      ! The area of a thickness category in a cell [L2 ~> m2].
  real :: area_h       ! The masked area of a column [L2 ~> m2].
  type(EFP_type) :: &
    fresh_water_in_EFP, & ! The total mass of fresh water added by surface fluxes
                       ! since the last time that write_ice_statistics was called [kg],
                       ! in extended fixed point (EFP) form.
    net_salt_in_EFP, & !   The total salt added by surface fluxes since the last
                       ! time that write_ice_statistics was called [kgSalt], in EFP form.
    net_heat_in_EFP, & !   The total heat added by surface fluxes since the last
                       ! time that write_ice_statistics was called [J], in EFP form.
    mass_EFP, &        ! The total water mass of the sea ice and snow and ponds atop it [kg].
    mass_chg_EFP, &    ! The change in total sea ice mass of fresh water since
                       ! the last call to this subroutine [kg].
    mass_anom_EFP, &   ! The change in fresh water that cannot be accounted for
                       ! by the surface fluxes [kg].
    salt_EFP, &        ! The total amount of salt in the brine pockets in sea ice [kgSalt].
    salt_chg_EFP, &    ! The change in total sea ice salt since the last call
                       ! to this subroutine [kgSalt].
    salt_anom_EFP, &   ! The change in total salt that cannot be accounted for by
                       ! the surface fluxes divided by total mass [kgSalt].
    heat_EFP, &        ! The total amount of enthalpy in the sea ice, snow, brine pockets
                       ! and melt ponds [J].
    heat_chg_EFP, &    ! The change in total sea ice heat since the last call
                       ! to this subroutine [J].
    heat_anom_EFP      ! The change in heat that cannot be accounted for by the surface fluxes [J].
  type(EFP_type), dimension(14+MAX_FIELDS_) :: EFP_list ! An array of EFP types for joint global sums.

  ! These real versions of other variables are here for debugging.
  real    :: fresh_water_input  !   The total mass of fresh water added by surface fluxes since
                                ! the last time that write_ice_statistics was called [kg].
  real    :: net_salt_input     !   The total salt added by surface fluxes since the last
                                ! time that write_ice_statistics was called [kgSalt].
  real    :: net_heat_input     !   The total heat added by surface fluxes since the last
                                ! time that write_ice_statistics was called [J].

  real :: CFL_trans    ! A transport-based definition of the CFL number [nondim].
  real :: CFL_u, CFL_v ! Simple CFL numbers for u- and v- advection [nondim].
  real :: dt_CFL       ! The timestep for calculating the CFL number [T ~> s].
  real :: max_CFL      ! The maximum of the CFL numbers [nondim].
  logical :: check_col
  integer :: num_nc_fields  ! The number of fields that will actually go into the NetCDF file.
  integer :: i, j, k, isc, iec, jsc, jec, isr, ier, jsr, jer, L, m, nlay, ncat, hem
  integer :: start_of_day, num_days
  integer :: iyear, imonth, iday, ihour, iminute, isecond, itick ! For call to get_date()
  real    :: reday     ! A real representation of the output time.
  character(len=120) :: statspath_nc
  character(len=300) :: mesg
  character(len=48)  :: msg_start
  character(len=32)  :: mesg_intro, time_units, day_str, n_str, trunc_str

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

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isr = isc - (G%isd-1) ; ier = iec - (G%isd-1) ; jsr = jsc - (G%jsd-1) ; jer = jec - (G%jsd-1)
  ncat = IG%CatIce ; nlay = IG%NkIce
  check_col = .false. ; if (present(check_column) .and. CS%column_check) check_col = check_column

  I_nlay = 1.0 / nlay

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "write_ice_statistics: Module must be initialized before it is used.")

  nTr_stocks = 0
  if (present(tracer_CSp)) then
    call SIS_call_tracer_stocks(G, IG, tracer_CSp, IST%mH_ice, Tr_stocks, &
                                stock_names=Tr_names, stock_units=Tr_units, num_stocks=nTr_stocks)
  endif

! nTr_stocks = 0
! if (present(tracer_CSp)) then
!   call call_tracer_stocks(h, Tr_stocks, G, tracer_CSp, stock_names=Tr_names, stock_units=Tr_units, &
!                           num_stocks=nTr_stocks, got_min_max=Tr_minmax_got, global_min=Tr_min, global_max=Tr_max, &
!                           xgmin=Tr_min_x, ygmin=Tr_min_y, zgmin=Tr_min_z,&
!                           xgmax=Tr_max_x, ygmax=Tr_max_y, zgmax=Tr_max_z)
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
        call open_ASCII_file(CS%statsfile_ascii, trim(CS%statsfile), &
            action=APPEND_FILE)
      else
        call open_ASCII_file(CS%statsfile_ascii, trim(CS%statsfile), &
            action=WRITEONLY_FILE)
        if (abs(CS%timeunit - 86400.0) < 1.0) then
          write(CS%statsfile_ascii,'("  Step,",7x,"Day,",28x,"Area(N/S),",22x,"Extent(N/S),",27x,&
              &"Mass(N/S),",22x,"Heat(N/S),",14x,"Salinty(N/S),   Frac Mass Err,   Temp Err,   Salin Err")')
          write(CS%statsfile_ascii,'(12x,"[days]",31x,"[m2]",28x,"[m2]",34x,"[kg]",29x,&
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

  ! Calculate the global maximum advective CFL numbers.
  max_CFL = 0.0
  dt_CFL = max(CS%dt, 0.)
  if (IST%Cgrid_dyn) then
    if (allocated(IST%u_ice_C)) then ; do j=jsc,jec ; do I=isc-1,iec
      if (IST%u_ice_C(I,j) < 0.0) then
        CFL_trans = (-IST%u_ice_C(I,j) * dt_CFL) * (G%dy_Cu(I,j) * G%IareaT(i+1,j))
      else
        CFL_trans = (IST%u_ice_C(I,j) * dt_CFL) * (G%dy_Cu(I,j) * G%IareaT(i,j))
      endif
      max_CFL = max(max_CFL, CFL_trans)
    enddo ; enddo ; endif
    if (allocated(IST%v_ice_C)) then ; do J=jsc-1,jec ; do i=isc,iec
      if (IST%v_ice_C(i,J) < 0.0) then
        CFL_trans = (-IST%v_ice_C(i,J) * dt_CFL) * (G%dx_Cv(i,J) * G%IareaT(i,j+1))
      else
        CFL_trans = (IST%v_ice_C(i,J) * dt_CFL) * (G%dx_Cv(i,J) * G%IareaT(i,j))
      endif
      max_CFL = max(max_CFL, CFL_trans)
    enddo ; enddo ; endif
  elseif (allocated(IST%u_ice_B) .and. allocated(IST%v_ice_B)) then
    do J=jsc-1,jec ; do I=isc-1,iec
      CFL_u = abs(IST%u_ice_B(I,J)) * dt_CFL * G%IdxBu(I,J)
      CFL_v = abs(IST%v_ice_B(I,J)) * dt_CFL * G%IdyBu(I,J)
      max_CFL = max(max_CFL, CFL_u, CFL_v)
    enddo ; enddo
  endif
  call max_across_PEs(max_CFL)

  ! Set combinations of scalign factors that rescale back to MKS units for output
  area_scale = US%L_to_m**2
  mass_scale = US%L_to_m**2 * US%RZ_to_kg_m2
  salt_scale = US%L_to_m**2 * US%RZ_to_kg_m2 * US%S_to_ppt
  heat_scale = US%L_to_m**2 * US%RZ_to_kg_m2 * US%Q_to_J_kg

  ! The following quantities are to be written by hemisphere:
  !   Ice area, ice extent, Ice+snow mass, enthalpy, salt
  ice_area(:,:,:) = 0.0
  ice_extent(:,:,:) = 0.0
  col_mass(:,:,:) = 0.0
  col_heat(:,:,:) = 0.0
  col_salt(:,:,:) = 0.0

  enth_liq_0 = enthalpy_liquid_freeze(0.0, IST%ITV)
  do j=jsc,jec ; do i=isc,iec
    hem = 1 ; if (G%geolatT(i,j) < 0.0) hem = 2
    do k=1,ncat ; if (G%mask2dT(i,j) * IST%part_size(i,j,k) > 0.0) then
      area_pt = G%areaT(i,j) * G%mask2dT(i,j) * IST%part_size(i,j,k)

      ice_area(i,j,hem) = ice_area(i,j,hem) + area_pt * area_scale
      col_mass(i,j,hem) = col_mass(i,j,hem) + area_pt * mass_scale * &
                          (IST%mH_ice(i,j,k) + (IST%mH_snow(i,j,k) + &
                           IST%mH_pond(i,j,k))) ! mw/new - assumed pond heat/salt = 0

      col_heat(i,j,hem) = col_heat(i,j,hem) + area_pt * heat_scale * &
                          (IST%mH_snow(i,j,k) * IST%enth_snow(i,j,k,1) + &
                           IST%mH_pond(i,j,k) * enth_liq_0)
      do L=1,nlay
        col_heat(i,j,hem) = col_heat(i,j,hem) + area_pt * heat_scale * &
                            ((IST%mH_ice(i,j,k) * I_nlay) * IST%enth_ice(i,j,k,L))
        col_salt(i,j,hem) = col_salt(i,j,hem) + area_pt * salt_scale * &
                  ((0.001*IST%mH_ice(i,j,k) * I_nlay) * IST%sal_ice(i,j,k,L))
      enddo
    endif ; enddo
    if (allocated(IST%snow_to_ocn)) then ; if (IST%snow_to_ocn(i,j) > 0.0) then
      area_pt = G%areaT(i,j) * G%mask2dT(i,j)
      col_mass(i,j,hem) = col_mass(i,j,hem) + area_pt * mass_scale * IST%snow_to_ocn(i,j)
      col_heat(i,j,hem) = col_heat(i,j,hem) + area_pt * heat_scale * &
                              (IST%snow_to_ocn(i,j) * IST%enth_snow_to_ocn(i,j))
    endif ; endif
    if (ice_area(i,j,hem) > 0.1*area_scale*G%areaT(i,j)) ice_extent(i,j,hem) = area_scale*G%areaT(i,j)

  enddo ; enddo

  ! Combining the sums adds code complexity, but avoids multiple blocking all-PE updates.
  do hem=1,2
    EFP_list(0+hem) = reproducing_sum_EFP(ice_area(:,:,hem), isr, ier, jsr, jer, only_on_PE=.true.)
    EFP_list(2+hem) = reproducing_sum_EFP(ice_extent(:,:,hem), isr, ier, jsr, jer, only_on_PE=.true.)
    EFP_list(4+hem) = reproducing_sum_EFP(col_heat(:,:,hem), isr, ier, jsr, jer, only_on_PE=.true.)
    EFP_list(6+hem) = reproducing_sum_EFP(col_mass(:,:,hem), isr, ier, jsr, jer, only_on_PE=.true.)
    EFP_list(8+hem) = reproducing_sum_EFP(col_salt(:,:,hem), isr, ier, jsr, jer, only_on_PE=.true.)
  enddo
  EFP_list(11) = real_to_EFP(real(CS%ntrunc))
  do m=1,nTr_stocks ; EFP_list(11+m) = real_to_EFP(Tr_stocks(11+m)) ; enddo

  if (CS%previous_calls > 0) then
    do j=jsc,jec ; do i=isc,iec
      ! Convert the mass, heat and salt input per unit area into cell integrals in
      ! units of [kg], [J] or [kgSalt].
      area_h = G%areaT(i,j) * G%mask2dT(i,j)
      water_into_cell(i,j) = area_h * mass_scale * CS%water_in_col(i,j)
      heat_into_cell(i,j) = area_h * heat_scale * CS%heat_in_col(i,j)
      salt_into_cell(i,j) = area_h * salt_scale * CS%salt_in_col(i,j)
    enddo ; enddo
    EFP_list(12+nTr_stocks) = reproducing_sum_EFP(water_into_cell, isr, ier, jsr, jer, only_on_PE=.true.)
    EFP_list(13+nTr_stocks) = reproducing_sum_EFP(salt_into_cell, isr, ier, jsr, jer, only_on_PE=.true.)
    EFP_list(14+nTr_stocks) = reproducing_sum_EFP(heat_into_cell, isr, ier, jsr, jer, only_on_PE=.true.)

    call EFP_sum_across_PEs(EFP_list, 14+nTr_stocks)

    fresh_water_in_EFP = EFP_list(12+nTr_stocks)
    net_salt_in_EFP = EFP_list(13+nTr_stocks)
    net_heat_in_EFP = EFP_list(14+nTr_stocks)
  else
    call EFP_sum_across_PEs(EFP_list, 11+nTr_stocks)
  endif

  ! Unpack the list of sums
  do hem=1,2
    Area_NS(hem) = EFP_to_real(EFP_list(0+hem))
    Extent_NS(hem) = EFP_to_real(EFP_list(2+hem))
    heat_NS(hem) = EFP_to_real(EFP_list(4+hem))
    mass_NS(hem) = EFP_to_real(EFP_list(6+hem))
    salt_NS(hem) = EFP_to_real(EFP_list(8+hem))
  enddo
  heat_EFP = EFP_list(5) + EFP_list(6)  ; Heat = heat_NS(1) + heat_NS(2)
  mass_EFP = EFP_list(7) + EFP_list(8)  ; Mass = mass_NS(1) + mass_NS(2)
  salt_EFP = EFP_list(9) + EFP_list(10) ; Salt = salt_NS(1) + salt_NS(2)
  CS%ntrunc = EFP_to_real(EFP_list(11))
  do m=1,nTr_stocks ; Tr_stocks(11+m) = EFP_to_real(EFP_list(11+m)) ; enddo

  ! Find the changes and anomalies for error analysis in salt, enthalpy and mass.
  if (CS%previous_calls == 0) then
    Salt_chg = 0.0 ; Salt_anom = 0.0
    Heat_chg = 0.0 ; Heat_anom = 0.0
    mass_chg = 0.0 ; mass_anom = 0.0
  else
    Salt_chg_EFP = Salt_EFP - CS%salt_prev_EFP
    Salt_anom_EFP = Salt_chg_EFP - net_salt_in_EFP
    Salt_chg = EFP_to_real(Salt_chg_EFP) ; Salt_anom = EFP_to_real(Salt_anom_EFP)
    Heat_chg_EFP = Heat_EFP - CS%heat_prev_EFP
    Heat_anom_EFP = Heat_chg_EFP - net_heat_in_EFP
    Heat_chg = EFP_to_real(Heat_chg_EFP) ; Heat_anom = EFP_to_real(Heat_anom_EFP)

    mass_chg_EFP = mass_EFP - CS%mass_prev_EFP
    mass_anom_EFP = mass_chg_EFP - fresh_water_in_EFP
    mass_chg = EFP_to_real(mass_chg_EFP) ; mass_anom = EFP_to_real(mass_anom_EFP)

    ! net_salt_input needs to be accounted for if mass includes salt.
    ! mass_anom = mass_anom - EFP_to_real(net_salt_in_EFP)
  endif

  I_Mass = 0.0 ; if (Mass > 0.0) I_Mass = 1.0/Mass
  salinity_NS(:) = 0.0
  do hem=1,2 ; if (mass_NS(hem) > 0.0) salinity_NS(hem) = 1000.*(salt_NS(hem) / mass_NS(hem)) ; enddo

  ! All quantities have been calculated at this point.  Output messages are prepared next.

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
    write(CS%statsfile_ascii,'(A,", Area", 2(ES23.16), ", Ext", 2(es11.4), ", CFL", F6.3, &
                &", M",2(ES12.5),", Enth",2(ES13.5),", S ",2(f8.4),", Me ",ES9.2,&
                &", Te ",ES9.2,", Se ",ES9.2)') &
          trim(msg_start), Area_NS(1:2), Extent_NS(1:2), max_CFL, mass_NS(1:2), &
          heat_NS(1:2), salinity_NS(1:2), mass_anom * I_Mass, &
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
            trim(msg_start), Salt, Salt_chg, Salt_anom
      else
        write(*,'(A," Ice Salt: ",ES24.16,", Change: ",ES12.5," Error: ",ES12.5," (",ES8.1,")")') &
            trim(msg_start), Salt, Salt_chg, Salt_anom, Salt_anom/Salt
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

  if (check_col .and. (CS%previous_calls > 0)) then ; do j=jsc,jec ; do i=isc,iec
    hem = 1 ; if (G%geolatT(i,j) < 0.0) hem = 2
    heat_imb = (col_heat(i,j,hem) - CS%heat_cell_prev(i,j)) - heat_into_cell(i,j)
    mass_imb = (col_mass(i,j,hem) - CS%water_cell_prev(i,j)) - water_into_cell(i,j)
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

! call write_field(CS%statsfile_nc, CS%fields(1), real(CS%ntrunc), reday)
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

  ! Reset the cumulative fluxes and store the net properties.
  CS%ntrunc = 0
  CS%previous_calls = CS%previous_calls + 1
  if (CS%column_check) then ; do j=jsc,jec ; do i=isc,iec
    CS%water_cell_prev(i,j) = col_mass(i,j,1) + col_mass(i,j,2)
    CS%heat_cell_prev(i,j) = col_heat(i,j,1) + col_heat(i,j,2)
    CS%salt_cell_prev(i,j) = col_salt(i,j,1) + col_salt(i,j,2)
  enddo ; enddo ; endif

  CS%water_in_col(:,:) = 0.0
  CS%heat_in_col(:,:) = 0.0
  CS%salt_in_col(:,:) = 0.0

  CS%mass_prev_EFP = mass_EFP
  CS%salt_prev_EFP = Salt_EFP
  CS%heat_prev_EFP = Heat_EFP

end subroutine write_ice_statistics


!> Accumulate the net input of fresh water and heat through the bottom of the
!! sea-ice for conservation checks.
subroutine accumulate_bottom_input(IST, OSS, FIA, IOF, dt, G, US, IG, CS)
  type(SIS_hor_grid_type),    intent(in) :: G   !< The horizontal grid type
  type(ice_grid_type),        intent(in) :: IG  !< The sea-ice specific grid type
  type(ice_state_type),       intent(in) :: IST !< A type describing the state of the sea ice
  type(ocean_sfc_state_type), intent(in) :: OSS !< A structure containing the arrays that describe
                                                !! the ocean's surface state for the ice model.
  type(fast_ice_avg_type),    intent(in) :: FIA !< A type containing averages of fields
                                                !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(in) :: IOF !< A structure containing fluxes from the ice to
                                                !! the ocean that are calculated by the ice model.
  real,                       intent(in) :: dt  !< The amount of time over which to average [T ~> s].
  type(unit_scale_type),      intent(in) :: US  !< A structure with unit conversion factors
  type(SIS_sum_out_CS),       pointer    :: CS  !< The control structure returned by a previous call
                                                !! to SIS_sum_output_init.

  ! Local variables
  real :: Flux_SW ! Total shortwave flux [Q R Z T-1 ~> W m-2]
  real :: LI      ! Latent heat of fusion [Q ~> J kg-1]

  integer :: i, j, k, isc, iec, jsc, jec, ncat, b, nb

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  nb = size(IOF%flux_sw_ocn, 3)

  call get_SIS2_thermo_coefs(IST%ITV, Latent_fusion=LI)

  if (CS%dt < 0.0) CS%dt = dt

  do j=jsc,jec ; do i=isc,iec
    CS%water_in_col(i,j) = CS%water_in_col(i,j) - dt * &
           ( ((FIA%runoff(i,j) + FIA%calving(i,j)) + &
              (IOF%lprec_ocn_top(i,j) + IOF%fprec_ocn_top(i,j))) - IOF%evap_ocn_top(i,j) )

    Flux_SW = 0.0
    do b=2,nb,2 ! This sum combines direct and diffuse fluxes to preserve answers.
      Flux_SW = Flux_SW + (IOF%flux_sw_ocn(i,j,b-1) + IOF%flux_sw_ocn(i,j,b))
    enddo
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) - dt * &
          (Flux_SW + ((IOF%flux_lw_ocn_top(i,j) - IOF%flux_lh_ocn_top(i,j)) - IOF%flux_sh_ocn_top(i,j)) + &
           (-LI)*(IOF%fprec_ocn_top(i,j) + FIA%calving(i,j)) )
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) - (OSS%frazil(i,j) - FIA%frazil_left(i,j))
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + &
           ((IOF%Enth_Mass_in_atm(i,j) + IOF%Enth_Mass_in_ocn(i,j)) + &
            (IOF%Enth_Mass_out_atm(i,j) + IOF%Enth_Mass_out_ocn(i,j)) )
    if (allocated(IOF%transmutation_enth)) &
      CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + IOF%transmutation_enth(i,j)

    CS%salt_in_col(i,j) = CS%salt_in_col(i,j) + dt * IOF%flux_salt(i,j)
    if (allocated(IOF%transmutation_salt_flux)) &
      CS%salt_in_col(i,j) = CS%salt_in_col(i,j) + dt * IOF%transmutation_salt_flux(i,j)
  enddo ; enddo

end subroutine accumulate_bottom_input

!> Accumulate the net input of fresh water and heat through the top of the
!! sea-ice for conservation checks with the first phase of the updates
subroutine accumulate_input_1(IST, FIA, OSS, dt, G, US, IG, CS)
  type(ice_state_type),       intent(in) :: IST !< A type describing the state of the sea ice
  type(fast_ice_avg_type),    intent(in) :: FIA !< A type containing averages of fields
                                                !! (mostly fluxes) over the fast updates
  type(ocean_sfc_state_type), intent(in) :: OSS !< A structure containing the arrays that describe
                                                !! the ocean's surface state for the ice model.
  real,                       intent(in) :: dt  !< The amount of time over which to average [T ~> s].
  type(SIS_hor_grid_type),    intent(in) :: G   !< The horizontal grid type
  type(unit_scale_type),      intent(in) :: US  !< A structure with unit conversion factors
  type(ice_grid_type),        intent(in) :: IG  !< The sea-ice specific grid type
  type(SIS_sum_out_CS),       pointer    :: CS  !< The control structure returned by a previous call
                                                !! to SIS_sum_output_init.

  ! Local variables
  real :: area_pt  ! The fractional area of a thickness partition in a cell [nondim]
  real :: Flux_SW  ! Total shortwave flux [Q R Z T-1 ~> W m-2]
  integer :: i, j, k, isc, iec, jsc, jec, ncat, b, nb

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  nb = size(FIA%flux_sw_top, 4)

  !$OMP parallel do default(shared) private(area_pt,Flux_SW)
  do j=jsc,jec ; do k=1,ncat ; do i=isc,iec
    area_pt = IST%part_size(i,j,k)
    Flux_SW = 0.0
    do b=2,nb,2 ! This sum combines direct and diffuse fluxes to preserve answers.
      Flux_SW = Flux_SW + (FIA%flux_sw_top(i,j,k,b-1) + FIA%flux_sw_top(i,j,k,b))
    enddo
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + (dt * area_pt) * &
        ( Flux_SW * (1.0 - FIA%sw_abs_ocn(i,j,k)) + &
           (FIA%flux_lw_top(i,j,k) - FIA%flux_sh_top(i,j,k))  + &
           (-FIA%flux_lh_top(i,j,k)) + OSS%bheat(i,j))
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) - area_pt * (FIA%bmelt(i,j,k) + FIA%tmelt(i,j,k))
  enddo ; enddo ; enddo

end subroutine accumulate_input_1

!> Accumulate the net input of fresh water and heat through the top of the
!! sea-ice for conservation checks, with a second phase of the updates
subroutine accumulate_input_2(IST, FIA, IOF, OSS, part_size, dt, G, US, IG, CS)
  type(SIS_hor_grid_type), intent(inout) :: G   !< The horizontal grid type
  type(ice_grid_type),     intent(inout) :: IG  !< The sea-ice specific grid type
  type(ice_state_type),    intent(inout) :: IST !< A type describing the state of the sea ice
  type(fast_ice_avg_type),    intent(in) :: FIA !< A type containing averages of fields
                                                !! (mostly fluxes) over the fast updates
  type(ice_ocean_flux_type),  intent(in) :: IOF !< A structure containing fluxes from the ice to
                                                !! the ocean that are calculated by the ice model.
  type(ocean_sfc_state_type), intent(in) :: OSS !< A structure containing the arrays that describe
                                                !! the ocean's surface state for the ice model.
  real, dimension(SZI_(G),SZJ_(G),SZCAT0_(IG)), &
                              intent(in) :: part_size !< The fractional ice concentration within a
                                                !! cell in each thickness category [nondim], 0-1.
  real,                       intent(in) :: dt  !< The amount of time over which to average [T ~> s].
  type(unit_scale_type),      intent(in) :: US  !< A structure with unit conversion factors
  type(SIS_sum_out_CS),       pointer    :: CS  !< The control structure returned by a previous call
                                                !! to SIS_sum_output_init.

  ! Local variables
  real :: area_pt  ! The fractional area of a thickness partition in a cell [nondim]
  real :: Flux_SW  ! The total shortwave flux, summed across frequency and angular bands [Q R Z T-1 ~> W m-2]
  real :: pen_frac ! The fraction of the shortwave that is absorbed by the ocean [nondim]
  real :: LI       ! Latent heat of fusion [Q ~> J kg-1]
  integer :: i, j, k, m, isc, iec, jsc, jec, ncat, b, nb

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  nb = size(FIA%flux_sw_top, 4)

  ! This subroutine includes the accumulation of mass fluxes and heat fluxes
  ! into the ice that are known before SIS#_thermodynamics, as well the
  ! ice-top fluxes that will be passed on directly to the ocean.  It does
  ! not include the enthalpy changes due to net mass changes in the ice,
  ! as these are not yet known.

  call get_SIS2_thermo_coefs(IST%ITV, Latent_fusion=LI)
  !$OMP parallel do default(shared) private(area_pt)
  do j=jsc,jec ; do i=isc,iec
    ! Runoff and calving are passed directly on to the ocean.
    CS%water_in_col(i,j) = CS%water_in_col(i,j) + dt * (FIA%runoff(i,j) + FIA%calving(i,j))

    area_pt = IST%part_size(i,j,0)
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + (dt * area_pt) * &
          ((FIA%flux_lw_top(i,j,0) - FIA%flux_lh_top(i,j,0)) - FIA%flux_sh_top(i,j,0))

    ! These are mass fluxes that are simply passed through to the ocean.
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + dt * (-LI) * &
                      (area_pt * FIA%fprec_top(i,j,0) + FIA%calving(i,j))

  enddo ; enddo

  ! The terms that are added here include surface fluxes that will be passed
  ! directly on into the ocean.
  !$OMP parallel do default(shared) private(area_pt,pen_frac,Flux_SW)
  do j=jsc,jec ; do k=0,ncat ; do i=isc,iec
    area_pt = part_size(i,j,k)
    pen_frac = 1.0 ; if (k>0) pen_frac = FIA%sw_abs_ocn(i,j,k)
    Flux_SW = 0.0
    do b=2,nb,2 ! This sum combines direct and diffuse fluxes to preserve answers.
      Flux_SW = Flux_SW + (FIA%flux_sw_top(i,j,k,b-1) + FIA%flux_sw_top(i,j,k,b))
    enddo

    CS%water_in_col(i,j) = CS%water_in_col(i,j) + (dt * area_pt) * &
        ( (FIA%lprec_top(i,j,k) + FIA%fprec_top(i,j,k)) - FIA%evap_top(i,j,k) )
    CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + (dt * area_pt) * ( pen_frac*Flux_SW )

    if (k>0) &
      CS%heat_in_col(i,j) = CS%heat_in_col(i,j) + area_pt * &
         ((FIA%bmelt(i,j,k) + FIA%tmelt(i,j,k)) - dt*OSS%bheat(i,j))
  enddo ; enddo ; enddo

  ! Runoff and calving do not bring in salt, so there is no modification of salt_in_col

end subroutine accumulate_input_2

end module SIS_sum_output
