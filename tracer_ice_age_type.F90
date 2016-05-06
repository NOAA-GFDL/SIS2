module ice_age_tracer_type

use MOM_time_manager, only : time_type
use SIS_tracer_registry, only : SIS_tracer_registry_type
use SIS_diag_mediator, only: SIS_diag_ctrl
use MOM_io, only : vardesc

integer, parameter :: NTR_MAX = 2

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: ice_age_tracer_CS
  integer :: ntr    ! The number of tracers that are actually used.
  character(len = 200) :: IC_file ! The file in which the age-tracer initial values
                    ! can be found, or an empty string for internal initialization.
  type(time_type), pointer  :: Time ! A pointer to the ocean model's clock.
  type(SIS_tracer_registry_type), pointer :: TrReg => NULL()
  real, pointer :: tr(:,:,:,:) => NULL()   ! The array of tracers used in this
                                           ! subroutine, in g m-3?
  real, pointer :: tr_aux(:,:,:,:) => NULL() ! The masked tracer concentration
                                             ! for output, in g m-3.
  type(p3d), dimension(NTR_MAX) :: &
    tr_adx, &! Tracer zonal advective fluxes in g m-3 m3 s-1.
    tr_ady! Tracer meridional advective fluxes in g m-3 m3 s-1.

  real, dimension(NTR_MAX) :: &
    IC_val = 0.0, &    ! The (uniform) initial condition value.
    young_val = 0.0, & ! The value assigned to tr at the surface.
    land_val = -1.0, & ! The value of tr used where land is masked out.
    tracer_start_year  ! The year in which tracers start aging, or at which the
                       ! surface value equals young_val, in years.
  logical :: mask_tracers  ! If true, tracers are masked out in massless layers.
  logical :: tracers_may_reinit  ! If true, tracers may go through the
                           ! initialization code if they are not found in the
                           ! restart files.
  logical :: tracer_ages(NTR_MAX)
  logical :: new_ice_is_sink(NTR_MAX)
  logical :: do_areal_ice_age
  logical :: do_mass_ice_age
  real :: min_thick_age, min_conc_age

  integer, dimension(NTR_MAX) :: &
    id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1 ! Indices used for the diagnostic
                                                   ! mediator
  integer, dimension(NTR_MAX) :: nlevels

  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.

  type(vardesc) :: tr_desc(NTR_MAX)
end type ice_age_tracer_CS

end module ice_age_tracer_type
