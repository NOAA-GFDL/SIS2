!> Provides a common interface for jointly stepping SIS2 and MOM6, and will
!! evolve as a platform for tightly integrating the ocean and sea ice models.
module combined_ice_ocean_driver

! This file is a part of SIS2. See LICENSE.md for the license.

!-----------------------------------------------------------------------
!
! This module provides a common interface for jointly stepping SIS2 and
! MOM6, and will evolve as a platform for tightly integrating the ocean
! and sea ice models.

use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_error_handler, only : callTree_enter, callTree_leave
use MOM_file_parser,   only : param_file_type, open_param_file, close_param_file
use MOM_file_parser,   only : read_param, get_param, log_param, log_version
use MOM_io,            only : file_exists, close_file, slasher, ensembler
use MOM_io,            only : open_namelist_file, check_nml_error
use MOM_time_manager,  only : time_type, time_type_to_real !, operator(>)

use ice_model_mod,   only : ice_data_type, ice_model_end
use ice_model_mod,   only : update_ice_slow_thermo, update_ice_dynamics_trans
use ocean_model_mod, only : update_ocean_model, ocean_model_end
use ocean_model_mod, only : ocean_public_type, ocean_state_type, ice_ocean_boundary_type

use coupler_types_mod, only : coupler_type_send_data, coupler_type_data_override
use coupler_types_mod, only : coupler_type_copy_data
use data_override_mod, only : data_override
use mpp_domains_mod,   only : domain2D, mpp_get_layout, mpp_get_compute_domain

implicit none ; private

public :: update_slow_ice_and_ocean, ice_ocean_driver_init, ice_ocean_driver_end

!> The control structure for the combined ice-ocean driver
type, public :: ice_ocean_driver_type ; private
  logical :: CS_is_initialized = .false. !< If true, this module has been initialized
  logical :: single_MOM_call  !< If true, advance the state of MOM with a single
                              !! step including both dynamics and thermodynamics.
                              !! If false, the two phases are advanced with
                              !! separate calls. The default is true.
end type ice_ocean_driver_type

contains

!=======================================================================
!>  This subroutine initializes the combined ice ocean coupling control type.
subroutine ice_ocean_driver_init(CS, Time_init, Time_in)
  type(ice_ocean_driver_type), pointer       :: CS        !< The control structure for combined ice-ocean driver
  type(time_type),             intent(in)    :: Time_init !< The start time for the coupled model's calendar.
  type(time_type),             intent(in)    :: Time_in   !< The time at which to initialize the coupled model.

  ! Local variables
  integer, parameter :: npf = 5 ! Maximum number of parameter files
  character(len=240) :: &
    output_directory = ' ', &      ! Directory to use to write the model output.
    parameter_filename(npf) = ' '  ! List of files containing parameters.
  character(len=240) :: output_dir ! A directory for logging output.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "ice_ocean_driver_init"  ! This module's name.
!     real :: Time_unit   ! The time unit in seconds for ENERGYSAVEDAYS.
!     character(len=48)  :: stagger
  type(param_file_type) :: param_file !< A structure to parse for run-time parameters
  integer :: unit, io, ierr, valid_param_files

  namelist /ice_ocean_driver_nml/ output_directory, parameter_filename

  call callTree_enter("ice_ocean_driver_init(), combined_ice_ocean_driver.F90")
  if (associated(CS)) then
    call MOM_error(WARNING, "ice_ocean_driver_init called with an associated "// &
                    "ice_ocean_driver_CS structure. Model is already initialized.")
    return
  endif
  allocate(CS)

  ! Read the relevant input parameters.
  if (file_exists('input.nml')) then
    unit = open_namelist_file(file='input.nml')
  else
    call MOM_error(FATAL, 'Required namelist file input.nml does not exist.')
  endif

  ierr=1 ; do while (ierr /= 0)
    read(unit, nml=ice_ocean_driver_nml, iostat=io, end=10)
    ierr = check_nml_error(io, 'ice_ocean_driver_nml')
  enddo
10 call close_file(unit)

  output_dir = trim(slasher(ensembler(output_directory)))
  valid_param_files = 0
  do io=1,npf ; if (len_trim(trim(parameter_filename(io))) > 0) then
    call open_param_file(trim(parameter_filename(io)), param_file, &
                         component="Ice_Ocean_driver", doc_file_dir=output_dir)
    valid_param_files = valid_param_files + 1
  endif ; enddo
  if (valid_param_files == 0) call MOM_error(FATAL, "There must be at least "//&
       "1 valid entry in input_filename in ice_ocean_driver_nml in input.nml.")

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "SINGLE_MOM_CALL", CS%single_MOM_call, &
                 "If true, advance the state of MOM with a single step \n"//&
                 "including both dynamics and thermodynamics.  If false, \n"//&
                 "the two phases are advanced with separate calls.", default=.true.)

!  OS%is_ocean_pe = Ocean_sfc%is_ocean_pe
!  if (.not.OS%is_ocean_pe) return

!  call get_param(param_file, mdl, "TIMEUNIT", Time_unit, &
!                 "The time unit for ENERGYSAVEDAYS.", units="s", default=86400.0)

  call close_param_file(param_file)
  CS%CS_is_initialized = .true.

  call callTree_leave("ice_ocean_driver_init(")
end subroutine ice_ocean_driver_init


!>   The subroutine update_slow_ice_and_ocean uses the forcing already stored in
!! the ice_data_type to advance both the sea-ice (and icebergs) and ocean states
!! for a time interval coupling_time_step.
subroutine update_slow_ice_and_ocean(CS, Ice, Ocn, Ocean_sfc, Ice_ocean_boundary, &
                                     time_start_update, coupling_time_step)
  type(ice_ocean_driver_type), &
                           pointer       :: CS   !< The control structure for this driver
  type(ice_data_type),     intent(inout) :: Ice  !< The publicly visible ice data type
  type(ocean_state_type),  pointer       :: Ocn  !< The internal ocean state and control structures
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< The publicly visible ocean surface state type
  type(ice_ocean_boundary_type), &
                           intent(inout) :: Ice_ocean_boundary !< A structure containing the various forcing
                                                               !! fields going from the ice to the ocean
                                                               !! The arrays of this type are intent out.
  type(time_type),         intent(in)    :: time_start_update  !< The time at the beginning of the update step
  type(time_type),         intent(in)    :: coupling_time_step !< The amount of time over which to advance
                                                               !! the ocean and ice

  real :: time_step         ! The time step of a call to step_MOM in seconds.

  call callTree_enter("update_ice_and_ocean(), combined_ice_ocean_driver.F90")
  time_step = time_type_to_real(coupling_time_step)

!  if (time_start_update /= CS%Time) then
!    call MOM_error(WARNING, "update_ice_and_ocean: internal clock does not "//&
!                            "agree with time_start_update argument.")
!  endif
  if (.not.associated(CS)) then
    call MOM_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ice_ocean_driver_type. ice_ocean_driver_init must be "//  &
                    "called first to allocate this structure.")
  endif
  if (.not.associated(Ocn)) then
    call MOM_error(FATAL, "update_ocean_model called with an unassociated "// &
                    "ocean_state_type structure. ocean_model_init must be "//  &
                    "called first to allocate this structure.")
  endif

  if (.not.(Ocean_sfc%is_ocean_pe .and. Ice%slow_ice_pe)) call MOM_error(FATAL, &
        "update_slow_ice_and_ocean can only be called from PEs that handle both "//&
        "the ocean and the slow ice processes.")

  if (.not.same_domain(Ocean_sfc%domain, Ice%slow_Domain_NH)) &
    call MOM_error(FATAL, "update_slow_ice_and_ocean can only be used if the "//&
        "ocean and slow ice layouts and domain sizes are identical.")

  !### Add clocks of the various calls.

  call update_ice_slow_thermo(Ice)

  call update_ice_dynamics_trans(Ice)

!    call mpp_clock_begin(fluxIceOceanClock)
  call direct_flux_ice_to_IOB( time_start_update, Ice, Ice_ocean_boundary )
!    call mpp_clock_end(fluxIceOceanClock)

  if (CS%single_MOM_call) then
    call update_ocean_model(Ice_ocean_boundary, Ocn, Ocean_sfc, &
                            time_start_update, coupling_time_step )
  else
    !### This is here as a temporary measure to avoid using newer arguments
    !### to update_ocean_model.
    call update_ocean_model(Ice_ocean_boundary, Ocn, Ocean_sfc, &
                            time_start_update, coupling_time_step )
!### This pair of calls works properly with MOM6 in place of the single call above.
!    call update_ocean_model(Ice_ocean_boundary, Ocn, Ocean_sfc, time_start_update, &
!                            coupling_time_step,  update_dyn=.true., update_thermo=.false., &
!                            start_cycle=.true., end_cycle=.false., cycle_length=time_step)
!    call update_ocean_model(Ice_ocean_boundary, Ocn, Ocean_sfc, time_start_update, &
!                            coupling_time_step,  update_dyn=.false., update_thermo=.true., &
!                            start_cycle=.false., end_cycle=.true., cycle_length=time_step)
  endif

  call callTree_leave("update_ice_and_ocean()")
end subroutine update_slow_ice_and_ocean

!> same_domain returns true if two domains use the same list of PEs and have
!! the same size computational domains.
function same_domain(a, b)
  type(domain2D), intent(in) :: a !< The first domain in the comparison
  type(domain2D), intent(in) :: b !< The second domain in the comparison
  integer :: isize_a, jsize_a, isize_b, jsize_b
  integer :: layout_a(2), layout_b(2)
  logical :: same_domain

  ! This does a limited number of checks for consistent domain sizes.

  call mpp_get_layout(a, layout_a)
  call mpp_get_layout(b, layout_b)
  same_domain = ((layout_a(1) == layout_b(1)) .and. (layout_a(2) == layout_b(2)))

  call mpp_get_compute_domain(a, xsize=isize_a, ysize=jsize_a)
  call mpp_get_compute_domain(b, xsize=isize_b, ysize=jsize_b)

  same_domain = same_domain .and. ((layout_a(1) == layout_b(1)) .and. &
                                   (layout_a(2) == layout_b(2)))

end function same_domain

!> This subroutine does a direct copy of the fluxes from the ice data type into
!! a ice-ocean boundary type on the same grid.
subroutine direct_flux_ice_to_IOB( Time, Ice, IOB )
  type(time_type),    intent(in)    :: Time !< Current time
  type(ice_data_type),intent(in)    :: Ice  !< A derived data type to specify ice boundary data
  type(ice_ocean_boundary_type), &
                      intent(inout) :: IOB !< A derived data type to specify
                                    !!  properties and fluxes passed from ice to ocean

  integer :: i, j, is, ie, js, je, i_off, j_off, n, m
  logical :: used

  ! Do a direct copy of the ice surface fluxes into the Ice-ocean-boundary type.

  if (ASSOCIATED(IOB%u_flux)) IOB%u_flux(:,:) = Ice%flux_u(:,:)
  if (ASSOCIATED(IOB%v_flux)) IOB%v_flux(:,:) = Ice%flux_v(:,:)
  if (ASSOCIATED(IOB%p     )) IOB%p(:,:) = Ice%p_surf(:,:)
  if (ASSOCIATED(IOB%mi    )) IOB%mi(:,:) = Ice%mi(:,:)
  if (ASSOCIATED(IOB%t_flux)) IOB%t_flux(:,:) = Ice%flux_t(:,:)
  if (ASSOCIATED(IOB%salt_flux)) IOB%salt_flux(:,:) = Ice%flux_salt(:,:)
  if (ASSOCIATED(IOB%sw_flux_nir_dir)) IOB%sw_flux_nir_dir(:,:) = Ice%flux_sw_nir_dir(:,:)
  if (ASSOCIATED(IOB%sw_flux_nir_dif)) IOB%sw_flux_nir_dif (:,:) = Ice%flux_sw_nir_dif (:,:)
  if (ASSOCIATED(IOB%sw_flux_vis_dir)) IOB%sw_flux_vis_dir(:,:) = Ice%flux_sw_vis_dir(:,:)
  if (ASSOCIATED(IOB%sw_flux_vis_dif)) IOB%sw_flux_vis_dif (:,:) = Ice%flux_sw_vis_dif (:,:)
  if (ASSOCIATED(IOB%lw_flux)) IOB%lw_flux(:,:) = Ice%flux_lw(:,:)
  if (ASSOCIATED(IOB%lprec)) IOB%lprec(:,:) = Ice%lprec(:,:)
  if (ASSOCIATED(IOB%fprec)) IOB%fprec(:,:) = Ice%fprec(:,:)
  if (ASSOCIATED(IOB%runoff)) IOB%runoff(:,:) = Ice%runoff(:,:)
  if (ASSOCIATED(IOB%calving)) IOB%calving(:,:) = Ice%calving
  if (ASSOCIATED(IOB%stress_mag)) IOB%stress_mag(:,:) = Ice%stress_mag(:,:)
  if (ASSOCIATED(IOB%ustar_berg)) IOB%ustar_berg(:,:) = Ice%ustar_berg(:,:)
  if (ASSOCIATED(IOB%area_berg)) IOB%area_berg(:,:) = Ice%area_berg(:,:)
  if (ASSOCIATED(IOB%mass_berg)) IOB%mass_berg(:,:) = Ice%mass_berg(:,:)
  if (ASSOCIATED(IOB%runoff_hflx)) IOB%runoff_hflx(:,:) = Ice%runoff_hflx(:,:)
  if (ASSOCIATED(IOB%calving_hflx)) IOB%calving_hflx(:,:) = Ice%calving_hflx(:,:)
  if (ASSOCIATED(IOB%q_flux)) IOB%q_flux(:,:) = Ice%flux_q(:,:)

  ! Extra fluxes
  call coupler_type_copy_data(Ice%ocean_fluxes, IOB%fluxes)

  ! These lines allow the data override code to reset the fluxes to the ocean.
  call data_override('OCN', 'u_flux',    IOB%u_flux   , Time )
  call data_override('OCN', 'v_flux',    IOB%v_flux   , Time)
  call data_override('OCN', 't_flux',    IOB%t_flux   , Time)
  call data_override('OCN', 'q_flux',    IOB%q_flux   , Time)
  call data_override('OCN', 'salt_flux', IOB%salt_flux, Time)
  call data_override('OCN', 'lw_flux',   IOB%lw_flux  , Time)
  call data_override('OCN', 'sw_flux_nir_dir', IOB%sw_flux_nir_dir, Time)
  call data_override('OCN', 'sw_flux_nir_dif', IOB%sw_flux_nir_dif, Time)
  call data_override('OCN', 'sw_flux_vis_dir', IOB%sw_flux_vis_dir, Time)
  call data_override('OCN', 'sw_flux_vis_dif', IOB%sw_flux_vis_dif, Time)
  call data_override('OCN', 'lprec',     IOB%lprec    , Time)
  call data_override('OCN', 'fprec',     IOB%fprec    , Time)
  call data_override('OCN', 'runoff',    IOB%runoff   , Time)
  call data_override('OCN', 'calving',   IOB%calving  , Time)
  call data_override('OCN', 'runoff_hflx',  IOB%runoff_hflx   , Time)
  call data_override('OCN', 'calving_hflx', IOB%calving_hflx  , Time)
  call data_override('OCN', 'p',         IOB%p        , Time)
  call data_override('OCN', 'mi',        IOB%mi       , Time)
  if (ASSOCIATED(IOB%stress_mag) ) &
    call data_override('OCN', 'stress_mag', IOB%stress_mag, Time )
  !Are these if statements needed, or does data_override routine check if variable is associated?
  if (ASSOCIATED(IOB%ustar_berg)) &
    call data_override('OCN', 'ustar_berg', IOB%ustar_berg, Time)
  if (ASSOCIATED(IOB%area_berg) ) &
    call data_override('OCN', 'area_berg',  IOB%area_berg , Time)
  if (ASSOCIATED(IOB%mass_berg) ) &
    call data_override('OCN', 'mass_berg',  IOB%mass_berg , Time)

  ! Override and output extra fluxes of tracers or gasses
  call coupler_type_data_override('OCN', IOB%fluxes, Time )

  call coupler_type_send_data(IOB%fluxes, Time )

end subroutine direct_flux_ice_to_IOB

!=======================================================================
!>   The subroutine ice_ocean_driver_end terminates the model run, saving
!! the ocean and slow ice states in restart files and deallocating any data
!! associated with the ocean and slow ice.
subroutine ice_ocean_driver_end(CS, Ice, Ocean_sfc, Ocn, Time)
  type(ice_ocean_driver_type), pointer   :: CS   !< The control structure for combined ice-ocean driver
  type(ice_data_type),     intent(inout) :: Ice  !< The publicly visible ice data type
  type(ocean_state_type),  pointer       :: Ocn  !< The internal ocean state and control structures
  type(ocean_public_type), intent(inout) :: Ocean_sfc !< The publicly visible ocean surface state type
  type(time_type),         intent(in)    :: Time !< The model time, used for writing restarts

  call ice_model_end(Ice)

  call ocean_model_end(Ocean_sfc, Ocn, Time)

  if (associated(CS)) deallocate(CS)

end subroutine ice_ocean_driver_end

end module combined_ice_ocean_driver
