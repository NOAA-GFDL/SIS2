module SIS_sum_output_type
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
!   This file contains the control structure for SIS_sum_output.F90.  It has   !
! been moved temporarily to its own file to avoid a predecessor cycle.         !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_coms, only : EFP_type
use MOM_time_manager, only : time_type

implicit none ; private


type, public :: SIS_sum_out_CS ! ; private

  real    :: mass_prev          !   The total sea ice mass the last time that
                                ! write_ice_statistics was called, in kg.
  real    :: fresh_water_input  !   The total mass of fresh water added by
                                ! surface fluxes since the last time that
  real    :: salt_prev          !   The total amount of salt in the sea ice the last
                                ! time that write_ice_statistics was called, in PSU kg.
  real    :: net_salt_input     !   The total salt added by surface fluxes since
                                ! the last time that write_ice_statistics was called,
                                ! in PSU kg.
  real    :: heat_prev          !   The total amount of heat in the sea ice the last
                                ! time that write_ice_statistics was called, in Joules.
  real    :: net_heat_input     !   The total heat added by surface fluxes since
                                ! the last time that write_ice_statistics was called,
                                ! in Joules.
  real, dimension(:,:), pointer :: &
    water_in_col, &             !   The water, heat, and salt that have been
    heat_in_col, &              ! input to the ice and snow in a column since
    salt_in_col, &              ! the last time that write_ice_statistics was
                                ! called, in kg m-2, J m-2, and kg m-2.
    water_col_prev, &           !   The column integrated water, heat, and salt
    heat_col_prev, &            ! that were in the ice and snow the last time
    salt_col_prev               ! that write_ice_statistics was called,
                                ! in kg m-2, J m-2, and kg m-2.

  type(EFP_type) :: &
    fresh_water_in_EFP, &       ! These are extended fixed point versions of the
    net_salt_in_EFP, &          ! correspondingly named variables above.
    net_heat_in_EFP, heat_prev_EFP, salt_prev_EFP, mass_prev_EFP
  real    :: dt                 ! The baroclinic dynamics time step, in s.
  real    :: timeunit           !   The length of the units for the time
                                ! axis, in s.
  type(time_type) :: Start_time ! The start time of the simulation.
                                ! Start_time is set in MOM_initialization.F90
  logical :: column_check       ! If true, enable the column by column heat and
                                ! mass conservation check
  real    :: imb_tol            ! The tolerance for imbalances to be flagged by
                                ! column_check, nondim.
  integer :: maxtrunc           ! The number of truncations per ice statistics
                                ! save interval at which the run is stopped.
  logical :: write_stdout       ! If true, periodically write sea ice statistics
                                ! to stdout to allow the progress to be seen.
  logical :: write_stocks       ! If true, write the integrated tracer amounts
                                ! to stdout when the statistics files are written.
  integer :: previous_calls = 0 ! The number of times write_ice_statistics has been called.
  integer :: prev_n = 0         ! The value of n from the last call.
  integer, pointer :: ntrunc    ! The number of times the velocity has been truncated
                                ! since the last call to write_ice_statistics.
!  integer :: statsfile_nc       ! NetCDF id of the statistics file.
  integer :: statsfile_ascii    ! The unit number of the ascii version of the statistics file.
!  type(fieldtype), dimension(NUM_FIELDS+MAX_FIELDS_) :: &
!             fields             ! fieldtype variables for the output fields.
  character(len=200) :: statsfile  ! The name of the statistics file with path.

end type SIS_sum_out_CS


end module SIS_sum_output_type
