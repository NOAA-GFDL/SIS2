!> Implements relaxation regions in SIS2 model
!> for sea ice concentration (partial area) and thickness by categories
!> The algorithm can be used to impose open boundary conditions
!> for sea ice thickness and partial area in regional SIS2 applications,
!> or for correcting the fields within the domain.
!>
!> Dmitry Dukhovskoy NOAA OAR PSL 2025
!>
module SIS_sponge

!! Module to read in Time of the ice fields, ice concentration, ice thickness, relaxation time scale
use MOM_coms,          only : sum_across_PEs, max_across_PEs
use MOM_coms,          only : PE_here   !! debugging
use MOM_time_manager,  only : time_type, set_date, get_time, get_date
use MOM_unit_scaling,  only : unit_scale_type
use ice_grid,          only : ice_grid_type

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : MOM_get_verbosity
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : file_exists, MOM_read_data, slasher
use MOM_io,            only : axis_info
use MOM_interpolate,   only : init_external_field, get_external_field_info, time_interp_external_init
use MOM_interpolate,   only : time_interp_external
use MOM_interpolate,   only : external_field

use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_types,         only : ice_state_type, IST_chksum, IST_bounds_check, total_sfc_flux_type
use SIS_types,         only : ocean_sfc_state_type, ice_ocean_flux_type, fast_ice_avg_type
use SIS2_ice_thm,      only : SIS2_ice_thm_CS, SIS2_ice_thm_init, SIS2_ice_thm_end
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs, enthalpy_liquid_freeze
use SIS2_ice_thm,      only : enth_from_TS, Temp_from_En_S, enthalpy_liquid, calculate_T_freeze
use SIS_optics,        only : VIS_DIR, VIS_DIF, NIR_DIR, NIR_DIF    ! debugging only, delete later
use SIS_utils,         only : is_NaN

implicit none; private

#include <SIS2_memory.h>

public initialize_icerelax_file, apply_isponge, set_up_isponge_field, SIS_sponge_end

!> A structure for creating arrays of pointers to 3D arrays
type, public :: p3d; private
  integer :: nz_data                             !< The number of vertical levels in the input field.
  integer :: num_tlevs                           !< The number of time records contained in the file
  real, dimension(:,:,:), pointer :: p => NULL() !< A pointer to a 3D array [various]
  character(len=15)               :: fld_name    !< Name of the ice field being relaxed
end type p3d
!> A structure for creating arrays of pointers to 2D arrays
type, public :: p2d; private
  type(external_field) :: field !< Time interpolator field handle
  integer :: ncat_data          !< The number of sea ice categories
  integer :: num_tlevs          !< The number of time records contained in the file
  real :: scale = 1.0           !< A multiplicative factor by which to rescale input data [various]
  real, dimension(:,:), pointer :: p => NULL()   !< A pointer to a 2D array [various]
  character(len=:), allocatable  :: name         !< The name of the input field
  character(len=:), allocatable  :: long_name    !< The long name of the input field
  character(len=:), allocatable  :: unit         !< The unit of the input field
  type(axis_info),  allocatable  :: axes_data(:) !< Axis types for the input field
                                                 !! name, longname, cartesian("X","Y",...) ax_size,...
end type p2d
!
!> A structure for 2D arrays
type, public :: f2d
  real, allocatable, dimension(:,:) :: fld
end type f2d
!> A structure for 3D arrays
type, public :: f3d
  real, allocatable, dimension(:,:,:) :: fld3
end type f3d

!> This control structure holds memory and parameters for the SIS_sponge module
type, public :: isponge_CS ; private
  logical, public :: use_isponge = .false.  !< If true, ice tracer fields may be relaxed somewhere in the domain
  integer         :: num_col         !< The number of relaxation points within the computational domain.
  integer, public :: fldno = 0       !< The number of fields which have already been
                                     !! registered by calls to set_up_isponge_field

  integer, pointer :: col_i(:) => NULL()         !< Array of the i-indcs of each of the columns being relaxed.
  integer, pointer :: col_j(:) => NULL()         !< Array of the j-indcs of each of the columns being relaxed.
  real, pointer    :: Iresttime_col(:) => NULL() !< The inverse restoring time of each column [T-1 ~> s-1].

  type(p3d) :: var(MAX_FIELDS_RLX_)     !< Pointers to the fields that will be relaxed
  type(p2d) :: Ref_val(MAX_FIELDS_RLX_) !< Relaxation values - The values to which the fields are
                                        ! relaxed distributed by ice cats (linear_index, ice_cat)
  type(f2d) :: Old_val(MAX_FIELDS_RLX_) !< Keep old values of relaxed fields prior to relaxation
                                        !! debug only, will need to get rid off later
  type(f2d) :: Ref_orig(MAX_FIELDS_RLX_) !< Relaxation values original input fields, i.e.
                                         !! on 2d grid not distributed by ice cats.
  logical :: time_varying_sponges       !< True if using newer sponge code
  logical :: spongeDataOngrid           !< True if the sponge data are on the model horizontal grid
end type isponge_CS

contains

!> This subroutine sets the inverse restoration time (Idamp) for sea ice fields and
!! the values towards which an arbitrary
!! number of tracers should be restored within the relaxation zone.
subroutine initialize_icerelax_file(param_file, G, IG, CS, US, IST, Time)
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  type(isponge_CS),        pointer    :: CS         !< A pointer to the SIS_isponge control structure
                                                    !! for this module
  type(unit_scale_type),   intent(in) :: US         !< A structure with unit conversion factors
  type(ice_state_type),    intent(in) :: IST        !< A type describing the state of the sea ice
  type(time_type),         intent(in) :: Time       !< The sea-ice model's clock,

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G))  :: Irelax !< The inverse of the restoring time [T-1 ~> s-1].
  real, allocatable, dimension(:,:) :: rlx_H  !< A temporary array for reading relax target ice thickness
                                              !! mean grid cell value  kg/m [R Z L ~> kg m-1]
  real, allocatable, dimension(:,:) :: rlx_C  !< A temporary array for reading relax target ice partial area

  integer, parameter :: verb_msg = 8 !< verbosity level for messages
  integer :: i, j, k, is, ie, js, je, ncat
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: year     !< The current model year
  integer :: day      !< The current model year-day
  integer :: second   !< The second of the day
  integer :: mon, hr, minute, itick  !< Time variables
  integer :: verbosity   !< MOM verbosity level
  real    :: max_rlxrate !< The strongest relaxation over the domain [T ~> s]
  real    :: rho_ice     !< The nominal density of sea ice [R ~> kg m-3]
  integer, dimension(4) :: siz
  character(len=40) :: ithck_var, iarea_var, rlxrate_var, rlx_unit
  character(len=40) :: mdl = "initialize_icerelax_file"
  character(len=50) :: rlx_long_name
  character(len=200) :: relaxrate_file, state_file  !< relax filenames: inverse time, target fields
  character(len=200) :: filename, inputdir          !< Strings for file/path and path.
  character(len=256) :: mesg

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  Irelax = 0.0

  verbosity = MOM_get_verbosity()

  call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(param_file, mdl, "ISPONGE_RELAX_FILE", relaxrate_file, &
                 "The name of the file with the sponge relaxation rates.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "ISPONGE_STATE_FILE", state_file, &
                 "The name of the file with the state to relax toward.", &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "ISPONGE_ITHCK_VAR", ithck_var, &
                 "The name of the ice thickness variable in "//&
                 "ISPONGE_STATE_FILE.", default="ithkn")
  call get_param(param_file, mdl, "ISPONGE_IAREA_VAR", iarea_var, &
                 "The name of the ice partial area variable in "//&
                 "ISPONGE_STATE_FILE.", default="iarea")
  call get_param(param_file, mdl, "ISPONGE_RLXRATE_VAR", rlxrate_var, &
                 "The name of the relaxation rate variable in "//&
                 "ISPONGE_RELAX_FILE.", default="relax_rate")

  ! Read in relaxation rate, s-1, for ice thickness and partial area
  filename = trim(inputdir)//trim(relaxrate_file)
  call log_param(param_file, mdl, "INPUTDIR/ISPONGE_RELAX_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call SIS_error(FATAL, " initialize_icerelax_file: Unable to open "//trim(filename))

  call MOM_read_data(filename, rlxrate_var, Irelax(:,:), G%Domain, scale=US%s_to_T)

  ! Check overall ice relax. rate only if verbosity allows printing the diagnostics
  if (verbosity > verb_msg) then
    max_rlxrate =  maxval(Irelax*US%T_to_s)
    call max_across_PEs(max_rlxrate)
    call get_date(Time, year, mon, day, hr, minute, second, itick)
    write(mesg,'("SIS Time:",i6,2("/",i2.2),1x,3(":",i2.2),"; max(Irelax)=",D13.4," s-1")') &
        year, mon, day, hr, minute, second, max_rlxrate
    call SIS_mesg(mesg, verb_msg)
  endif

  call initialize_isponge(param_file, Irelax, G, IG, CS)

  ! Now register all of the fields which are nudged in the relaxation region.
  filename = trim(inputdir)//trim(state_file)
  call log_param(param_file, mdl, "INPUTDIR/ISPONGE_STATE_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call SIS_error(FATAL, " initialize_icerelax_files: Unable to open "//trim(filename))
!

  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  call SIS_mesg('initialize_icerelax_file: Calling set_up_isponge_field: mH_ice', verb_msg)
  call set_up_isponge_field(filename, ithck_var, Time, 1, IG%CatIce, G, IG, US, IST%mH_ice, CS, &
       'mH_ice', rlx_long_name='ice_thickness', rlx_unit='kg m-2', scale=US%m_to_Z * rho_ice)
  call SIS_mesg('initialize_icerelax_file: Calling set_up_isponge_field: part_size', verb_msg)
  call set_up_isponge_field(filename, iarea_var, Time, 0, IG%CatIce, G, IG, US, IST%part_size, CS, &
         'part_size', rlx_long_name='partial_area', rlx_unit='none')

end subroutine initialize_icerelax_file

!> This subroutine determines the number of points which are within ice relaxation region in
!! this computational domain.  Only points that have positive values of
!! Iresttime and which mask2dT indicates are ocean points are included as the
!! relaxation points.
subroutine initialize_isponge(param_file, Iresttime, G, IG, CS, time_var_rlx, sponge_ongrid)
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: Iresttime  !< The inverse of the restoring time [T-1 ~> s-1].
  type(isponge_CS),        pointer    :: CS         !< A pointer to the SIS_isponge control structure
                                                    !! for this module
  logical, optional, intent(in) :: time_var_rlx, sponge_ongrid !< place-holders, currently both true

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  ! Local variables
  integer, parameter :: verb_msg = 8 !< verbosity level for messages
  character(len=40)  :: mdl = "initialize_isponge"  ! This module's name.
  character(len=256) :: mesg
  logical :: use_isponge
  integer :: i, j, k, m, n, b, nb, isc, iec, jsc, jec, ncat
  integer :: col, total_isponge_cols

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_sponge: initialize_isponge called with "// &
                            "an associated control structure.")
    return
  endif

  ! Set default, read and log parameters
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "SIS_SPONGE", use_isponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  if (.not.use_isponge) return
  allocate(CS)

  write(mesg,'(A,": SIS_SPONGE IS ON")') trim(mdl)
  call SIS_mesg(trim(mesg))
  CS%time_varying_sponges = .true.  ! TODO: add option to SIS_input for not time varying rlx fields
  CS%spongeDataOngrid = .true.
  if (present(time_var_rlx)) CS%time_varying_sponges = time_var_rlx
  if (present(sponge_ongrid)) CS%spongeDataOngrid = sponge_ongrid

  CS%use_isponge = use_isponge

  CS%num_col = 0 ; CS%fldno = 0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((Iresttime(i,j) > 0.0) .and. (G%mask2dT(i,j) > 0.0)) &
      CS%num_col = CS%num_col + 1
  enddo ; enddo

  if (CS%num_col > 0) then
    allocate(CS%Iresttime_col(CS%num_col), source=0.0)
    allocate(CS%col_i(CS%num_col), source=0)
    allocate(CS%col_j(CS%num_col), source=0)
    col = 1
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((Iresttime(i,j) > 0.0) .and. (G%mask2dT(i,j) > 0.0)) then
        CS%col_i(col) = i ; CS%col_j(col) = j
        CS%Iresttime_col(col) = Iresttime(i,j)
        col = col +1
      endif
    enddo ; enddo
  endif

  total_isponge_cols = CS%num_col
  call sum_across_PEs(total_isponge_cols)

  write(mesg,'(A,": total isponge cols=",i8)') trim(mdl), total_isponge_cols
  call SIS_mesg(mesg, verb_msg)
  call log_param(param_file, mdl, "!Total isponge columns", total_isponge_cols, &
                 "The total number of ice columns where relaxation is applied.")

end subroutine initialize_isponge

!> This subroutine stores the reference (target) profile for the SIS variable whose
!! address is given by f_ptr. Reference profile = values towards which the
!! SIS field is being relaxed to.
!! Current version assumes 2D input fields.
subroutine set_up_isponge_field(filename, fieldname, Time, kdS, kdE, G, IG, US, f_ptr, CS, &
                                rlxfld_name, rlx_long_name, rlx_unit, scale)
  character(len=*),        intent(in) :: filename   !< The name of the file with the
                                                    !! time varying field data
  character(len=*),        intent(in) :: fieldname  !< The name of the field in the file
                                                    !! with the time varying field data
  type(time_type),         intent(in) :: Time       !< The current model time
  integer,                 intent(in) :: kdS, kdE   !< start/end indices for Cats, conc(0)= open water area
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in) :: US         !< A structure with unit conversion factors
  real, dimension(SZI_(G), SZJ_(G), kdS:kdE), &
                   target, intent(in) :: f_ptr      !< a pointer to the field which will be relaxed [various]
                                                    !! note: IST%part_size(isd:ied, jsd:jed, 0:CatIce)
  type(isponge_CS),     pointer       :: CS         !< A pointer to the control structure for this module that
                                                    !! is set by a previous call to initialize_sponge.
  character(*),            intent(in) :: rlxfld_name   !< Name of the relaxed field
  character(len=*),        optional,  &
                           intent(in) :: rlx_long_name !< The long name of the tracer field
                                                       !! if not given, use the sp_name
  character(len=*),        optional,  &
                           intent(in) :: rlx_unit      !< The unit of the tracer field
                                                       !! if not given, use 'none'
  real,          optional, intent(in) :: scale !< A factor by which to rescale the input data, including any
                                               !! contributions due to dimensional rescaling [various ~> 1].

  ! Local variables
  integer, parameter :: verb_msg = 8 !< verbosity level for messages
  integer, dimension(4) :: fld_sz    !< Dimensions of the input ice target fields
  integer :: isd, ied, jsd, jed      !< Start/end indices of the domain
  integer :: i, j, k, col            !< Dummy indices
  integer :: CatIce                  !< The number of ice categories
  character(len=256) :: mesg         !< String for error messages
  character(len=256) :: long_name    !< The long name of the tracer field
  character(len=256) :: unit         !< The unit of the tracer field
  character(len=40)  :: mdl          !< This module name

  long_name = rlxfld_name; if (present(rlx_long_name)) long_name = rlx_long_name
  unit = 'none'; if (present(rlx_unit)) unit = rlx_unit

  CatIce = IG%CatIce
  mdl = 'set_up_isponge_field'

  if (.not.associated(CS)) return
  ! initialize time interpolator module
  call time_interp_external_init()
  isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
  CS%fldno = CS%fldno + 1
  write(mesg,'("set_up_isponge: fldno=",I0)') CS%fldno
  call SIS_mesg(mesg, verb_msg)
  if (CS%fldno > MAX_FIELDS_RLX_) then
    write(mesg,'("Increase MAX_FIELDS_RLX_ to at least ",I3," in SIS_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_sponge." )') CS%fldno
    call SIS_error(FATAL,"set_up_isponge_field: "//mesg)
  endif
  ! get a unique time interp id for this field. Ice relax target fields are on-grid
  if (CS%spongeDataOngrid) then
    call SIS_mesg("set_up_isponge_field: calling init_external_field", verb_msg)
    CS%Ref_val(CS%fldno)%field = init_external_field(filename, fieldname, MOM_domain=G%Domain, &
               verbose=.true.)
  else
    call SIS_error(FATAL,"set_up_isponge_field: SIS2 relaxation fields on a not-native grid not implemented")
  endif
  CS%Ref_val(CS%fldno)%name = rlxfld_name
  CS%Ref_val(CS%fldno)%long_name = long_name
  CS%Ref_val(CS%fldno)%unit = unit
  fld_sz(1:4) = -1
  call get_external_field_info(CS%Ref_val(CS%fldno)%field, size=fld_sz, axes=CS%Ref_val(CS%fldno)%axes_data)
  CS%Ref_val(CS%fldno)%ncat_data = CatIce ! individual relax fields should have same # of categories
  CS%Ref_val(CS%fldno)%num_tlevs = fld_sz(4)
  CS%Ref_val(CS%fldno)%scale = 1.0 ; if (present(scale)) CS%Ref_val(CS%fldno)%scale = scale

  ! initializes the target profile array for this field
  ! for all columns which will be masked
  select case(trim(rlxfld_name))
    case('part_size')
      allocate(CS%Ref_val(CS%fldno)%p(CS%num_col,0:CatIce), source=0.0)
    case('mH_ice')
      allocate(CS%Ref_val(CS%fldno)%p(CS%num_col,CatIce), source=0.0)
    case default
      write(mesg,'(A," Unknown relaxation field: ",A," setting default pointer dimensions =",2(i5,1x))') &
           trim(mdl),trim(rlxfld_name),CS%num_col,CatIce
      call SIS_error(WARNING,trim(mesg))
      allocate(CS%Ref_val(CS%fldno)%p(CS%num_col,CatIce), source=0.0)
  end select
  allocate(CS%Old_val(CS%fldno)%fld(CS%num_col,CatIce), source=0.0)
  allocate(CS%Ref_orig(CS%fldno)%fld(isd:ied,jsd:jed), source=0.0)

  CS%var(CS%fldno)%p => f_ptr    ! points to the actual ice fields that will be relaxed
  CS%var(CS%fldno)%fld_name = rlxfld_name

  write(mesg,'("set_up_isponge_field: ",A," fld_sz(1:4)=",4(I5,1x)," scale=",f14.6)') &
               rlxfld_name, fld_sz(1:4), CS%Ref_val(CS%fldno)%scale
  call SIS_mesg(mesg, verb_msg)

end subroutine set_up_isponge_field

!> This subroutine applies relaxation (nudging) to ice thickness (by categories) and ice concentration
!! tracers for every column where the relaxation time scale > 0.
subroutine apply_isponge(dt_slow, CS, G, IG, IST, US, OSS, Time)
  real,                      intent(in)  :: dt_slow   !< The amount of time covered by this call [T ~> s].
  type(isponge_CS),          pointer     :: CS     !< A pointer that is set to point to the ice sponge control
                                                   !! structure for this module
  type(ice_grid_type),       intent(in)  :: IG     !< The sea-ice specific grid type
  type(SIS_hor_grid_type),   intent(in)  :: G      !< The horizontal grid type
  type(ice_state_type),   intent(inout)  :: IST    !< A type describing the state of the sea ice
  type(unit_scale_type),     intent(in)  :: US     !< A structure with unit conversion factors
  type(ocean_sfc_state_type), intent(in) :: OSS    !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(time_type),           intent(in)  :: Time   !< The current model date

  ! Local variables
  real :: damp         !< The timestep times the local damping coefficient [nondim].
  real :: I1pdamp      !< I1pdamp is 1/(1 + damp). [nondim]
  real :: dt           !< time step [s]
  real :: s_ice_bulk   !< ice bulk S for filling S values in the newly created ice
  real, allocatable :: sice(:), tfi(:)  ! Arrays for ice salinity and freezing temperature
  real :: enth_ice     !< The enthalpy of ice [Q ~> J kg-1]
  real :: Tfrz         !< The freezeing temperature of sea water [C ~> degC]
  real :: coeff        !< conversion coefficient from unscaled to scaled units
  real :: enth_Tfrz    !< Ice enthalpy at the freezing temperature for a given ice salinity
  character(len=40)  :: mdl = "apply_isponge"  ! This subroutine's name.
  character(len=256) :: mesg
  character(len=15)  :: fld_name
  real    :: Idt_slow                   !< The inverse of the thermodynamic step [T-1 ~> s-1].
  real    :: iconc_old, ithk_old, iconc_new, ithk_new !< old and updated ice thkn and conc
  real    :: iconc_tot, iconc_tot_old   !< aggregated old and updated ice conc
  real    :: ithk_tot_new, ithk_tot_old !< aggregated old and updated ice thickness
  real    :: ice_salin                  !< average ice column S [gSalt kg-1]
  real, dimension(:,:), allocatable :: data_in  !< A buffer for storing the full 2-d time-interpolated array
  real, dimension(:,:), allocatable :: mask_in  !< A 2-d mask for extended input grid [nondim]
  real    :: I_Nk                 !< The inverse of the number of vertical ice layers
  real    :: part_water           !< partial area of open water
  integer :: id, jd, kd, jdp      !< Input dataset data sizes
  type(axis_info), dimension(4) :: axes_data
  integer, parameter :: verb_msg = 8 !< verbosity level for messages
  integer :: i, j, k, l, m, col
  integer :: ii, jj, iiG, jjG
  integer :: CatIce             !< The number of sea ice categories.
  integer :: NkIce              !< The number of vertical layers within the sea ice.
  integer :: verbosity   !< MOM verbosity level
  integer :: is, ie, js, je     !< compute domain indices
  integer :: isg, ieg, jsg, jeg !< global extent
  integer :: isd, ied, jsd, jed !< data domain indices
  integer :: nid, njd, isdG, iedG, jsdG, jedG
  integer, dimension(4) :: fld_sz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isg = G%isg ; ieg = G%ieg ; jsg = G%jsg ; jeg = G%jeg
  nid  = G%ied - G%isd + 1
  njd  = G%jed - G%jsd + 1
  isdG = G%isd_global; iedG = isdG + nid
  jsdG = G%jsd_global; jedG = jsdG + njd

  verbosity = MOM_get_verbosity()

  CatIce = IG%CatIce
  NkIce  = IG%NkIce
  I_Nk  = 1. / NkIce
  dt = dt_slow*US%T_to_s
  s_ice_bulk = 3.0*US%ppt_to_S

  if (CS%num_col == 0) return

  ! First get relax fields and interp. in time:
  allocate(data_in(isd:ied,jsd:jed))
  allocate(sice(NkIce), tfi(NkIce), source=-999.)
  do m=1,CS%fldno
    if (verbosity > verb_msg) then
      call time_interp_external(CS%Ref_val(m)%field, Time, data_in, verbose=.true.)
    else
      call time_interp_external(CS%Ref_val(m)%field, Time, data_in, verbose=.false.)
    endif
    CS%Ref_orig(m)%fld(:,:) = data_in(:,:)
  enddo

  ! Convert input 2D fields --> 3D ice thicknesses and concentration by categories
  ! Input hice is aggregated ice volume per m2, i.e. hice=voli=sum(hice(k)*cice(k))
  ! In each category:
  ! scale ice thickness m --> kg m-2 and unscale US%m_to_Z
  call distribute_ice2cats(CS, IG, G)

  do col=1,CS%num_col
    i = CS%col_i(col) ; j = CS%col_j(col)
    damp = dt * CS%Iresttime_col(col); I1pdamp = 1.0 / (1.0 + damp)
    do k=1,IG%CatIce
      do m=1,CS%fldno
        CS%Old_val(m)%fld(col,k) = CS%var(m)%p(i,j,k)
        CS%var(m)%p(i,j,k) = I1pdamp * &
           (CS%var(m)%p(i,j,k) + CS%Ref_val(m)%p(col,k)*damp)
      enddo
      ! Adjust enth and S in the newly formed ice if needed:
      ! Note ice enthalpy < 0
      do l=1,NkIce
        if (IST%sal_ice(i,j,k,l) < s_ice_bulk) &
          IST%sal_ice(i,j,k,l) = s_ice_bulk
        sice(l) = IST%sal_ice(i,j,k,l)
      enddo

      ! Enth should be at least enth(T_freez)
      ! Make ice T below T frz and/or keep at ocean SST if it is < ice Tfrz
      ! to prevent rapid ice melt in the relaxation zone
      call calculate_T_Freeze(sice, tfi, IST%ITV)
      tfi = min(tfi-0.1*US%degC_to_C, OSS%SST_C(i,j)*US%degC_to_C)

      do l=1,NkIce
        enth_ice = IST%enth_ice(i,j,k,l)
        enth_Tfrz = enth_from_TS(tfi(l), sice(l), IST%ITV)

        if (enth_ice > enth_Tfrz) then
          IST%enth_ice(i,j,k,l) = enth_Tfrz
        endif
      enddo
    enddo  ! CatIce

    ! Adjust open water partial area:
    do m=1,CS%fldno
      if (CS%var(m)%fld_name(1:9)=='part_size') then
        part_water = 1.0 - sum(CS%var(m)%p(i,j,1:CatIce))
        part_water = max(0.0, part_water)
        part_water = min(1.0, part_water)
        CS%var(m)%p(i,j,0) = part_water
      endif
    enddo

    ! Remove all snow if ice conc or thickness = 0
    do k=1,IG%CatIce
      do m=1,CS%fldno
        if (CS%var(m)%p(i,j,k) < 1.e-10) IST%mH_snow(i,j,k)=0.0
      enddo
    enddo
  enddo

  if (allocated(sice)) deallocate(sice)
  if (allocated(tfi)) deallocate(tfi)
  if (allocated(data_in)) deallocate(data_in)

end subroutine apply_isponge

!> Map local indices to global
subroutine local_to_global_indx(G, i, j, iiG, jjG)
  type(SIS_hor_grid_type),   intent(in)  :: G          !< The horizontal grid type
  integer, intent(in)                    :: i, j       !< Indices on the current tile
  integer, intent(out)                   :: iiG, jjG   !< Global indices

  iiG = G%isd_global + (i-1)  ; jjG = G%jsd_global + (j-1)

end subroutine local_to_global_indx

!> Redistribute input target 2D hice and iconc into ice thickness categories
!! Thick-to-thin algorithm:
!! Assign all ice initially to the thickest ice category
!! based on original ice thickness (hice).
!! Redistribute a small volume of ice from the thickest
!! category into the thinner categories.
!! This avoids 0s ice concentration and thickness.
!! Check that the total ice volume and concentration are conserved
!! during the distribution.
subroutine distribute_ice2cats(CS, IG, G, scaled, eps_err)
  type(isponge_CS),        pointer     :: CS       !< A pointer that is set to point to the ice sponge control
                                                   !! structure for this module
  type(ice_grid_type),     intent(in)  :: IG       !< The sea-ice specific grid type
  type(SIS_hor_grid_type), intent(in)  :: G        !< The horizontal grid type
  logical, optional,       intent(in)  :: scaled   !< true if input hice (m) converted to kg/m2 and scaled
                                                   !! default = .false. hice units: (m) = m3/m2 vol/m2
  real, optional,          intent(in)  :: eps_err  !< error allowed for hice, cice after redistribution

  ! local variables
  integer :: isd, ied, jsd, jed      !< Data domain indices
  integer :: isdG, jsdG              !< Global indx, start pnts. data domain
  integer :: m, i, j, k, col         !< Dummy indices
  integer :: iiG, jjG                !< Global indices
  integer :: CatIce                  !< The number of sea ice categories.
  integer :: icat0                   !< Ice thickness category where to assign original
                                     !! hice and cice to begin distrbution

  real, allocatable, dimension(:,:) :: cice2d, hice2d !< 2D arrays for input target ice area
                                                      !! and volume per unit area [m3*m-2]
  real, allocatable, dimension(:) :: hLim_vals    !< Ice thickness categories [m]
  real, allocatable, dimension(:) :: hcat, ccat   !< 1D arrays for thkn and conc in each ice category
  real, allocatable, dimension(:) :: volcat       !< Ice volume per unit area in each category [m3*m-2]
  real :: scale_cf            !< Scaling factor applied to the input ice target relaxation fields
  real :: Iscale              !< Inverse scale to "unscale" the data
  real :: hice, cice          !< Target relax. aggregated ice volume/area [m3*m-2] and partial area
  real :: hice_k              !< Ice thkn in category k [m]
  real :: hice_tot, cice_tot  !< Aggregated ice volume per unit area and partial area
  real :: eps0                !< Error allowed for hice, cice after redistribution
  real :: ck_min              !< The minimum ice concentration used in the
                              !! lower categories (wrt icat0) to distribute hice/cice
  real :: ccat_k, hcat_k, dch_k !< Ice concentration, volume, volume change in category k
  real :: cnew                !< Updated partial area in the thickest category
  real :: rmm                 !< Dummy variable
  real :: htot_min            !< The minimum ice volume [m3*m-2]  required
                              !! to distribute over the ice cats (1,icat0)
  real :: part_water          !< Partial area of open water
  character(len=40)  :: mdl = "distribute_ice2cats"  ! This module's name`
  character(len=256) :: mesg
  logical :: scaled_hice         !< Flag indicating if the target variable has been scaled
  logical :: err_hice, err_cice  !< Checks for conserved total ice thkn and conc

  scaled_hice = .false.
  if (present(scaled)) scaled_hice = scaled
  eps0 = 1.e-10  !< 0-checking
  if (present(eps_err)) eps0 = eps_err

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isdG = G%isd_global;   jsdG = G%jsd_global
  CatIce = IG%CatIce

  allocate(cice2d(isd:ied,jsd:jed), hice2d(isd:ied,jsd:jed))
  allocate(hLim_vals(CatIce+1), ccat(1:CatIce), hcat(1:CatIce), volcat(CatIce))

  hLim_vals(:) = IG%cat_thick_lim(:)  !< ice thkn cats are not scaled [m]

  cice2d=0.0 ; hice2d=0.0
  do m=1,CS%fldno
    scale_cf = CS%Ref_val(m)%scale
    Iscale = 1.0
    if (abs(1.-scale_cf) > 1.e-10 .and. scale_cf > 0.) &
      Iscale = 1.0/scale_cf
    select case (trim(CS%var(m)%fld_name))
      case('mH_ice')
        hice2d = CS%Ref_orig(m)%fld
        if (scaled_hice .and. abs(1.-scale_cf) > 1.e-10) &
            hice2d = CS%Ref_orig(m)%fld*Iscale  !< unscale input hice to original units (m) to find ice cat
      case('part_size')
        cice2d = CS%Ref_orig(m)%fld             !< partial area (conc) is not scaled
    end select
  enddo

  do col=1,CS%num_col
    i = CS%col_i(col) ; j = CS%col_j(col)
    hice = hice2d(i,j)    !< in input units [m3*m-2], also equivallent to cell-mean ice thickness
    cice = cice2d(i,j)    !< total partial area (conc)
    if (hice < 1.e-10 .or. cice < 1.e-10) then
      hice = 0.0 ; cice = 0.0
      do m=1,CS%fldno ; do k=1,CatIce
        select case (trim(CS%var(m)%fld_name))
          case('mH_ice')
            CS%Ref_val(m)%p(col,k) = 0.0
          case('part_size')
            CS%Ref_val(m)%p(col,k) = 0.0
            if (k == 1) CS%Ref_val(m)%p(col,0) = 1.0     !< open water partial area
        end select
      enddo ; enddo
      cycle
    endif

    ck_min = 1.e-2    !< conc in the lower cats, some small value but >> eps0
    icat0 = 1000
    hice_k = hice/cice   ! ice thickness in category k from ice volume (m3/m2)
    icat0 = find_icat(hice_k, hLim_vals)

    ! Check for ice cat. error:
    if (icat0 < 1 .or. icat0 > CatIce) then
      iiG = isdG + (i-1) ; jjG = jsdG + (j-1)
      write(mesg,'(A,"ERROR: iG,jG=",2(i4,1x)," hice=",D16.4," cice=",D16.4," hice_k=",D16.4)') &
             trim(mdl), iiG, jjG, hice, cice, hice_k
      call SIS_mesg(mesg, all_print=.true.)
      write(mesg,'(A," error: icat0 ",i2," hice=",f12.4," cice=",f12.4," hice_k=",f12.6)') &
            trim(mdl), hice, cice, hice_k
      call SIS_error(FATAL, trim(mesg))
    endif

    ! Adjust min conc in cats for low partial areas
    rmm = 1./float(icat0)
    if (cice < ck_min*float(icat0)) ck_min=cice*rmm

    ! Check if there is enough ice for distribution:
    htot_min = sum(ck_min*hLim_vals(1:icat0))
    if (hice < htot_min .or. cice < (ck_min*icat0) .or. (hice*cice) < 1.e-10) then
      do m=1,CS%fldno ; do k=1,CatIce
        select case (trim(CS%var(m)%fld_name))
          case('mH_ice')
            CS%Ref_val(m)%p(col,k) = 0.0
          case('part_size')
            CS%Ref_val(m)%p(col,k) = 0.0
            if (k == 1) CS%Ref_val(m)%p(col,0) = 1.0  ! open water fraction size
        end select
      enddo ; enddo
      cycle   ! not enough ice, skip the following lines
    endif

    ! Distribute ice thicknesses and conc by cats
    ! initial step: place all ice in icat0
    hcat = 0.0 ; ccat = 0.0 ; volcat = 0.0
    hcat(icat0) = hice_k ; ccat(icat0) = cice ; volcat(icat0) = hice
    iiG = isdG + (i-1) ; jjG = jsdG + (j-1)

    ! Sanity check: Ensure that the initial sea ice distribution conserves
    ! total ice volume and concentration.
    call check_hcice(hcat, ccat, hice, cice, iiG, jjG, str=' Initial ice distr. :')

    do k=1,icat0-1
      ccat_k = ck_min
      hcat_k = hLim_vals(k) + eps0
      dch_k = ccat_k*hcat_k   !< ice vol moved into this cat
      if (dch_k > volcat(icat0)) exit !< not enough ice left in the thickest cat
      ccat(k) = ccat_k
      hcat(k) = hcat_k
      volcat(k) = dch_k
      ! Update initial ice conc & vol in the thickest cat:
      cnew = ccat(icat0)-ccat(k)
      cnew = max(cnew, ck_min)
      ccat(icat0) = cnew
      volcat(icat0) = volcat(icat0)-dch_k
      hcat(icat0) = volcat(icat0)/ccat(icat0)
    enddo

    ! Register fields into control structure:
    do m=1,CS%fldno ; do k=1,CatIce
      select case (trim(CS%var(m)%fld_name))
        case('mH_ice')
          CS%Ref_val(m)%p(col,k) = hcat(k)*CS%Ref_val(m)%scale  !< scaled
        case('part_size')
          CS%Ref_val(m)%p(col,k) = ccat(k)
      end select
    enddo ; enddo

    ! Partial area of open water:
    do m=1,CS%fldno
      select case (trim(CS%var(m)%fld_name))
        case('part_size')
          part_water = 1.0 - sum(CS%Ref_val(m)%p(col,1:CatIce))
          part_water = max(0.0, part_water)
          part_water = min(1.0, part_water)
          CS%Ref_val(m)%p(col,0) = part_water
      end select
    enddo

    ! Check that the mean ice thickn and total conc. are conserved:
    hcat = 0.0 ; ccat = 0.0
    do m=1,CS%fldno ; do k=1,CatIce
      select case (trim(CS%var(m)%fld_name))
        case('mH_ice')
          if (abs(1.-CS%Ref_val(m)%scale) > 1.e-10 .and. &
              CS%Ref_val(m)%scale > 0.) then
            Iscale = 1./CS%Ref_val(m)%scale
            hcat(k) = CS%Ref_val(m)%p(col,k)*Iscale
          else
            hcat(k) = CS%Ref_val(m)%p(col,k)
          endif
        case('part_size')
          ccat(k) = CS%Ref_val(m)%p(col,k)
      end select
    enddo ; enddo

    call check_hcice(hcat, ccat, hice, cice, iiG, jjG, str=' After ice distr. :')

  enddo   !< do col

  if (allocated(hice2d)) deallocate(hice2d)
  if (allocated(cice2d)) deallocate(cice2d)
  if (allocated(hLim_vals)) deallocate(hLim_vals)
  if (allocated(ccat)) deallocate(ccat)
  if (allocated(hcat)) deallocate(hcat)
  if (allocated(volcat)) deallocate(volcat)

end subroutine distribute_ice2cats

!> Check if total ice thkn*conc and conc are conserved (i.e. equal to the original hice, cice)
subroutine check_hcice(hcat, ccat, hice, cice, iG, jG, str)
  real, intent(in) :: hcat(:)                    !< Ice thicknesses by categories
  real, intent(in) :: ccat(:)                    !< Ice concentrations by categories
  real, intent(in) :: hice                       !< Original aggregated ice volume per unit area [m3*m-2]
  real, intent(in) :: cice                       !< Original aggregated ice concentration
  character(len=*),  optional, intent(in) :: str !< A string for an error message

  ! Local variables
  real    :: hice_tot, cice_tot  !< Aggregated ice volume [m3*m-2] derived from hcat and ccat
  logical :: err_hice, err_cice  !< Errors of the total ice volume and concentration
  real    :: eps_err             !< Error tolerance margin
  integer :: iG, jG              !< Global ice indices for warning message
  integer :: CatIce
  character(len=40) :: msg_info
  character(len=200) :: mesg

  err_hice = .false. ; err_cice = .false.
  eps_err = 1.e-6      !< allow small error during distribution
  CatIce = size(hcat)
  call partial_area_total(ccat(1:CatIce), cice_tot)
  call ice_thkn_total(ccat(1:CatIce), hcat, hice_tot)

  if (abs(hice_tot-hice) > eps_err) err_hice=.true.
  if (abs(cice_tot-cice) > eps_err) err_cice=.true.

  ! Provide Error information
  msg_info = "check_hcice "
  if (present(str)) msg_info = trim(msg_info)//str
  if (err_hice) then
    write(mesg,'(A," iG, jG: ",2(i0,1x)," hice not conserved: ",f12.6," hice=",f12.6," cice=",f12.6)') &
          trim(msg_info), iG, jG, hice_tot, hice, cice
    call SIS_error(WARNING, trim(mesg))
  elseif (err_cice) then
    write(mesg,'(A," iG, jG: ",2(i0,1x)," cice not conserved: ",f12.6," cice=",f12.6," hice=",f12.6)') &
          trim(msg_info), iG, jG, cice_tot, cice, hice
    call SIS_error(WARNING, trim(mesg))
  endif

end subroutine check_hcice

!> Find ice thickness category for given ice thickness (m)
function find_icat(hice_k, hLim_vals) result (icat0)
  integer          :: icat0         !< The ice thkn category where hice_k belongs
  real, intent(in) :: hLim_vals(:)  !< ice thkn cats, not scaled [m]
  real, intent(in) :: hice_k        !< ice thkn in a category [m]

  ! Local variable
  integer :: k
  integer :: CatIce     !< The number of ice thkn cats

  CatIce = size(hLim_vals)-1
  icat0 = 1e6
  if (hice_k  >=  hLim_vals(CatIce)) then
    icat0 = CatIce
  elseif (hice_k  <  hLim_vals(1)) then
    icat0=1
  else
    do k=1,CatIce
      if (hice_k  >=  hLim_vals(k) .and. hice_k  <  hLim_vals(k+1)) then
        icat0 = k
        exit
      endif
    enddo
  endif

end function find_icat

!> The subroutine computes aggregated partial area
!! given a 1D array of cice(1:CatIce) partial areas by cats.
subroutine partial_area_total(cice_cat, cice_tot)
  real, intent(in)    :: cice_cat(:)   !< 1D array of partial areas by cats.
  real, intent(inout) :: cice_tot      !< Aggregated ice partial area of the grid cell

  ! Local variable
  integer :: k

  cice_tot = 0.0
  do k=1,size(cice_cat)
    cice_tot = cice_tot + cice_cat(k)
  enddo

end subroutine partial_area_total

!> The subroutine computes the aggregated ice volume per unit area [m3*m-2],
!! which is also equal to the grid cell mean ice thickness,
!! for 1D arrays of thiknesses and partial area by cats.
subroutine ice_thkn_total(cice_cat, hice_cat, hice_tot)
  real, intent(in)    :: cice_cat(:)  !< 1D array of partial areas by cats.
  real, intent(in)    :: hice_cat(:)  !< 1D array of ice thickn. by cats.
  real, intent(inout) :: hice_tot     !< Aggregated ice volume per unit area [m3*m-2]

  ! Local variables
  integer :: k

  hice_tot = 0.0
  do k=1,size(cice_cat)
    hice_tot = hice_tot + cice_cat(k)*hice_cat(k)
  enddo

end subroutine ice_thkn_total

!> Deallocate memory associated with the SIS_sponge module
subroutine SIS_sponge_end(CS)
  type(isponge_CS), pointer :: CS !< The ice sponge control structure that is deallocated here

  deallocate(CS)
end subroutine SIS_sponge_end

end module SIS_sponge

