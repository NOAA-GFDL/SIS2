!> Implements relaxation regions in SIS2 model
!> mainly for sea ice concentration (partial area) and thickness by categories
!> This is a quick fix for the problem caused by SIS closed boundaries 
!> unrealistic ice concentrations and thicknesses are simulated when sea ice is pulled off / piled up along the 
!> closed boundaries
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
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : file_exists, MOM_read_data, slasher
use MOM_io,            only : axis_info
use MOM_interpolate,   only : init_external_field, get_external_field_info, time_interp_external_init
use MOM_interpolate,   only : time_interp_external
use MOM_interpolate,   only : external_field     
                            
use SIS_diag_mediator, only : post_SIS_data, post_data=>post_SIS_data
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_sum_output,    only : SIS_sum_out_CS, write_ice_statistics! , SIS_sum_output_init
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
public global_to_local_ij, print_ice_thkn_conc

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
  integer, public :: itest, jtest    ! debugging, output at idices on PE
  integer         :: num_col         !< The number of relaxation points within the computational domain.
  integer, public :: fldno = 0       !< The number of fields which have already been
                                     !! registered by calls to set_up_sponge_field

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
!! the values towards which the interface heights and an arbitrary
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

  real, dimension(SZI_(G),SZJ_(G))  :: Irelax !< The inverse of the restoring time [T-1 ~> s-1].
  real, allocatable, dimension(:,:) :: rlx_H  !< A temporary array for reading relax target ice thickness
                                              !! mean grid cell value  kg/m [R Z L ~> kg m-1]
  real, allocatable, dimension(:,:) :: rlx_C  !< A temporary array for reading relax target ice partial area

  integer :: i, j, k, is, ie, js, je, ncat
  integer :: isd, ied, jsd, jed
  integer :: isc, iec, jsc, jec
  integer :: itestG, jtestG, itest, jtest
  integer :: year     !< The current model year
  integer :: day      !< The current model year-day
  integer :: second   !< The second of the day
  integer :: mon, hr, minute, itick
  integer :: start_of_day, num_days
  real :: max_rlxrate, rho_ice
  integer, dimension(4) :: siz
  character(len=40) :: ithck_var, iarea_var, rlxrate_var, rlx_unit
  character(len=40) :: mdl = "initialize_icerelax_file"
  character(len=50) :: rlx_long_name
  character(len=200) :: relaxrate_file, state_file  ! relax filenames: inverse time, target fields
  character(len=200) :: filename, inputdir ! Strings for file/path and path.
  character(len=256) :: mesg

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; ncat = IG%CatIce
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  Irelax = 0.0 ; itestG = 0 ; jtestG = 0 ; itest = 0 ; jtest = 0

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
                 "The name of the ice partial area/concentration variable in "//&
                 "ISPONGE_STATE_FILE.", default="ithkn")
  call get_param(param_file, mdl, "ISPONGE_RLXRATE_VAR", rlxrate_var, &
                 "The name of the relaxation rate variable in "//&
                 "ISPONGE_RELAX_FILE.", default="relax_rate")
  call get_param(param_file, mdl, "ISPONGE_ITEST", itestG, &
                 "I index of a test point to check ice relaxation dumped to log file")
  call get_param(param_file, mdl, "ISPONGE_JTEST", jtestG, &
                 "J index of a test point to check ice relaxation dumped to log file")

  ! Read in relaxation rate, s-1, for ice thickness and partial area
  filename = trim(inputdir)//trim(relaxrate_file)
  call log_param(param_file, mdl, "INPUTDIR/ISPONGE_RELAX_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call SIS_error(FATAL, " initialize_icerelax_file: Unable to open "//trim(filename))

  call MOM_read_data(filename, rlxrate_var, Irelax(:,:), G%Domain, scale=US%s_to_T) 
  max_rlxrate =  maxval(Irelax*US%T_to_s)
  call max_across_PEs(max_rlxrate)
  call get_date(Time, year, mon, day, hr, minute, second, itick)
  write(mesg,'("SIS Time:",i6,2("/",i2.2),1x,3(":",i2.2),"; max(Irelax)=",D13.4," s-1")') &
        year, mon, day, hr, minute, second, max_rlxrate
  call SIS_mesg(mesg) 
  call get_time(Time, start_of_day, num_days)

  if (itestG.gt.0 .and. jtestG.gt.0) &
    call global_to_local_ij(G, itestG, jtestG, itest, jtest)

  call SIS_mesg('initialize_icerelax_file: Calling initialize_isponge') 
  if (itest.gt.0 .and. jtest.gt.0) then
    call initialize_isponge(param_file, Irelax, G, IG, CS, itest=itest, jtest=jtest)
  else
    call initialize_isponge(param_file, Irelax, G, IG, CS)
  endif

  ! Now register all of the fields which are nudged in the relaxation region.
  filename = trim(inputdir)//trim(state_file)
  call log_param(param_file, mdl, "INPUTDIR/ISPONGE_STATE_FILE", filename)
  if (.not.file_exists(filename, G%Domain)) &
    call SIS_error(FATAL, " initialize_icerelax_files: Unable to open "//trim(filename))
!
  
  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  call SIS_mesg('initialize_icerelax_file: Calling set_up_isponge_field: mH_ice') 
  call set_up_isponge_field(filename, ithck_var, Time, 1, IG%CatIce, G, IG, US, IST%mH_ice, CS, &
       'mH_ice', rlx_long_name='ice_thickness', rlx_unit='kg m-2', scale=US%m_to_Z * rho_ice)
  call SIS_mesg('initialize_icerelax_file: Calling set_up_isponge_field: part_size') 
  call set_up_isponge_field(filename, iarea_var, Time, 0, IG%CatIce, G, IG, US, IST%part_size, CS, &
         'part_size', rlx_long_name='partial_area', rlx_unit='none')

end subroutine initialize_icerelax_file

!> This subroutine determines the number of points which are within ice relaxation region in
!! this computational domain.  Only points that have positive values of
!! Iresttime and which mask2dT indicates are ocean points are included as the
!! relaxation points.  
subroutine initialize_isponge(param_file, Iresttime, G, IG, CS, itest, jtest, time_var_rlx, sponge_ongrid)
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(ice_grid_type),     intent(in) :: IG         !< The sea-ice specific grid type
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: Iresttime  !< The inverse of the restoring time [T-1 ~> s-1].
  type(isponge_CS),        pointer    :: CS         !< A pointer to the SIS_isponge control structure
                                                    !! for this module
  integer, optional, intent(in) :: itest, jtest     !< test grid indices for debugging
  logical, optional, intent(in) :: time_var_rlx, sponge_ongrid !< place-holders, currently both true

  ! This include declares and sets the variable "version".
# include "version_variable.h"
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
! get_param (procedure --> get_param_logical) - checks if variable is set true in the param_file
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "SIS_SPONGE", use_isponge, &
                 "If true, sponges may be applied anywhere in the domain. "//&
                 "The exact location and properties of those sponges are "//&
                 "specified from MOM_initialization.F90.", default=.false.)

  if (use_isponge) then
    call SIS_mesg("initialize_isponge: use_isponge=True")
  else
    call SIS_mesg("initialize_isponge: use_isponge=False")
  endif
  if (.not.use_isponge) return
  allocate(CS)

  CS%time_varying_sponges = .true.  ! TODO: add option to SIS_input for not time varying rlx fields
  CS%spongeDataOngrid = .true.
  if (present(time_var_rlx)) CS%time_varying_sponges = time_var_rlx
  if (present(sponge_ongrid)) CS%spongeDataOngrid = sponge_ongrid

  CS%use_isponge = use_isponge
  if (present(itest) .and. present(jtest)) then
    write(mesg,'(A," itest/jtest =",2(i5,1x))') trim(mdl), itest, jtest
    write(*,'(A)') trim(mesg)
    CS%itest = itest
    CS%jtest = jtest
  else
    CS%itest = -999
    CS%jtest = -999
  endif

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
  call SIS_mesg(mesg)
  
  call log_param(param_file, mdl, "!Total isponge columns", total_isponge_cols, &
                 "The total number of ice columns where relaxation is applied.")

end subroutine initialize_isponge

!> This subroutine stores the reference profile for the SIS variable whose
!! address is given by f_ptr. Reference profile = values towards which the
!! SIS field is being relaxed to. 
!! Current version assumes 2D fields only (+ 1D for categories) such as Hice
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
  integer :: isd, ied, jsd, jed
  integer, dimension(4) :: fld_sz
  integer :: i, j, k, col, CatIce
  character(len=256) :: mesg      ! String for error messages
  character(len=256) :: long_name ! The long name of the tracer field
  character(len=256) :: unit      ! The unit of the tracer field
  character(len=40)  :: mdl       ! this module name

  long_name = rlxfld_name; if (present(rlx_long_name)) long_name = rlx_long_name
  unit = 'none'; if (present(rlx_unit)) unit = rlx_unit

  CatIce = IG%CatIce
  mdl = 'set_up_isponge_field'

  if (.not.associated(CS)) return
  ! initialize time interpolator module
  call time_interp_external_init()
  isd = G%isd; ied = G%ied; jsd = G%jsd; jed = G%jed
  CS%fldno = CS%fldno + 1
  write(mesg,'("set_up_isponge: fldno=",I)') CS%fldno
  call SIS_mesg(mesg)
  if (CS%fldno > MAX_FIELDS_RLX_) then
    write(mesg,'("Increase MAX_FIELDS_RLX_ to at least ",I3," in SIS_memory.h or decrease &
           &the number of fields to be damped in the call to &
           &initialize_sponge." )') CS%fldno
    call SIS_error(FATAL,"set_up_isponge_field: "//mesg)
  endif
  ! get a unique time interp id for this field. Ice relax target fields are on-grid
  if (CS%spongeDataOngrid) then
    call SIS_mesg("set_up_isponge_field: calling init_external_field")
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
  CS%Ref_val(CS%fldno)%ncat_data = CatIce !< individual relax fields should have same # of categories
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
      call SIS_error(WARNING,"print_ice_thkn_conc: "//mesg)
      allocate(CS%Ref_val(CS%fldno)%p(CS%num_col,CatIce), source=0.0)
  end select
  allocate(CS%Old_val(CS%fldno)%fld(CS%num_col,CatIce), source=0.0)
  allocate(CS%Ref_orig(CS%fldno)%fld(isd:ied,jsd:jed), source=0.0)

  CS%var(CS%fldno)%p => f_ptr    ! points to the actual ice fields that will be relaxed
  CS%var(CS%fldno)%fld_name = rlxfld_name

  write(mesg,'("set_up_isponge_field: ",A," fld_sz(1:4)=",4(I5,1x)," scale=",f14.6)') &
               rlxfld_name, fld_sz(1:4), CS%Ref_val(CS%fldno)%scale
  call SIS_mesg(mesg)

end subroutine set_up_isponge_field

!> This subroutine applies relaxation ("damping") to ice thickness (by categories) and ice concentration
!! tracers for every column where the relaxation time scale > 0.
subroutine apply_isponge(dt_slow, CS, G, IG, IST, US, OSS, Time)
  real,                      intent(in)  :: dt_slow   !< The amount of time covered by this call [T ~> s].
  type(isponge_CS),          pointer     :: CS     !< A pointer that is set to point to the ice sponge control
                                                   !! structure for this module
  type(ice_grid_type),       intent(in)  :: IG     !< The sea-ice specific grid type
  type(SIS_hor_grid_type),   intent(in)  :: G          !< The horizontal grid type
  type(ice_state_type),   intent(inout)  :: IST    !< A type describing the state of the sea ice
  type(unit_scale_type),     intent(in)  :: US     !< A structure with unit conversion factors
  type(ocean_sfc_state_type), intent(in) :: OSS    !< A structure containing the arrays that describe
                                                   !! the ocean's surface state for the ice model.
  type(time_type),           intent(in)  :: Time   !< The current model date

  ! Local variables
!  real :: H_rescale_ice, H_rescale_snow
  real :: damp         !< The timestep times the local damping coefficient [nondim].
  real :: I1pdamp      !< I1pdamp is 1/(1 + damp). [nondim]
  real :: p_old        !< debugging  
  real :: dt           !< time step in s
  real :: s_ice_bulk   !< ice bulk S for filling S values in the newly ceated ice 
  real, allocatable :: sice(:), tfi(:)
  real :: enth_ice, Tfrz, coeff, enth_Tfrz
  character(len=40)  :: mdl = "apply_isponge"  ! This subroutine's name.
  character(len=256) :: mesg
  character(len=15)  :: fld_name
  real    :: Idt_slow    !< The inverse of the thermodynamic step [T-1 ~> s-1].
  real    :: iconc_old, ithk_old, iconc_new, ithk_new
  real    :: iconc_tot, iconc_tot_old
  real    :: ithk_tot_new, ithk_tot_old
  real    :: dlt_iconc, dlt_ithk, enthalpy_ocn, enthalpy_ocn_tfrz
  real    :: dlt_salt, dlt_heat, dlt_water, dlt_snow
  real    :: dlt_ice           !< total change of ice due to conc and thickness relaxation
  real    :: ice_salin         !< average ice column S gSalt kg-1 
  real    :: water_ice_ocn, heat_ice_ocn, salt_ice_ocn
  real    :: enthalpy_ocn0
  real    :: dlt_enth
  logical :: f_debug
!
  real, dimension(:,:), allocatable :: data_in  !< A buffer for storing the full 2-d time-interpolated array
  real, dimension(:,:), allocatable :: mask_in  !< A 2-d mask for extended input grid [nondim]

  real    :: I_Nk
  real    :: part_water
  integer :: id, jd, kd, jdp      !< Input dataset data sizes
  type(axis_info), dimension(4) :: axes_data
  integer :: i, j, k, l, m, col
  integer :: ii, jj, iiG, jjG
  integer :: CatIce             !< The number of sea ice categories.
  integer :: NkIce              !< The number of vertical layers within the sea ice.
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
!
  f_debug = .true.
  CatIce = IG%CatIce
  NkIce  = IG%NkIce
  I_Nk  = 1. / NkIce
  dt = dt_slow*US%T_to_s
  s_ice_bulk = 3.0*US%ppt_to_S
  !dgr2rad = atan(1.0)/45.  

  if (CS%num_col == 0) return
  
! First get relax fields and interp. in time:
  allocate(data_in(isd:ied,jsd:jed))  
  allocate(sice(NkIce), tfi(NkIce), source=-999.)
  do m=1,CS%fldno
    call time_interp_external(CS%Ref_val(m)%field, Time, data_in, verbose=.true.)
    CS%Ref_orig(m)%fld(:,:) = data_in(:,:)
    ! Debug
    do col=1,CS%num_col
      i = CS%col_i(col) ; j = CS%col_j(col)
      if (CS%itest.eq.i .and. CS%jtest.eq.j) then
        iiG = isdG + (i-1)  ; jjG = jsdG + (j-1)
        write(mesg,'("apply_isponge: test i/j=",2(i5,1x)," time_iterp data_in=",f8.4)') &
        iiG, jjG, data_in(i,j)
        write(*,'(A)') trim(mesg)
      endif
    enddo
  enddo

  ! Convert input 2D fields --> 3D ice thicknesses and concentration by categories
  ! Input hice is "ice volume" per m2, i.e. hice=voli=sum(hice(k)*cice(k))
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
  !
  ! Diagnostics at test point if specified in SIS_input
        if (i==CS%itest .and. j==CS%jtest) then
          select case (trim(CS%var(m)%fld_name))
            case('mH_ice')    ; coeff = US%RZ_to_kg_m2
            case('part_size') ; coeff = 1.0
            case default
              write(mesg,'("SIS_sponge: Unknown relaxation field: ",A)') trim(CS%var(m)%fld_name)
              call SIS_error(FATAL,"apply_isponge: "//mesg)          
          end select
          write(mesg,'(A8," k=",I2," old:=",D12.4," new=",D12.4,&
                      " tau=",D12.4," refval=",D12.4," dt=",f7.1)') &
            CS%var(m)%fld_name(1:8), k, CS%Old_val(m)%fld(col,k)*coeff, &
            CS%var(m)%p(i,j,k)*coeff, &
            CS%Iresttime_col(col), CS%Ref_val(m)%p(col,k)*coeff, dt
          write(*,'(A)') trim(mesg)
        endif
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

        dlt_enth = 0.0
        if (enth_ice > enth_Tfrz) then
!          dlt_enth = enth_Tfrz - IST%enth_ice(i,j,k,l)
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
        if (CS%var(m)%p(i,j,k).lt.1.e-10) IST%mH_snow(i,j,k)=0.0
      enddo
    enddo
    
! Diagnostics at test point
    if (i.eq.CS%itest .and. j.eq.CS%jtest) then  
      iconc_tot = 0.0 ; iconc_tot_old = 0.0 ; ithk_tot_new = 0.0 ; ithk_tot_old = 0.0
      do k=1,IG%CatIce
        do m=1,CS%fldno
          fld_name = CS%var(m)%fld_name
          select case (trim(fld_name))
            case('mH_ice')
              ithk_old = CS%Old_val(m)%fld(col,k)/CS%Ref_val(m)%scale
              ithk_new = CS%var(m)%p(i,j,k)/CS%Ref_val(m)%scale
            case('part_size')
              iconc_old = CS%Old_val(m)%fld(col,k)
              iconc_new = CS%var(m)%p(i,j,k)
          end select
        enddo
        iconc_tot_old = iconc_tot_old + iconc_old
        iconc_tot    = iconc_tot + iconc_new
        ithk_tot_old = ithk_tot_old + ithk_old*iconc_old
        ithk_tot_new = ithk_tot_new + ithk_new*iconc_new
      enddo
!
      write(mesg, '("conc old=",f6.4," new=",f6.4," thick (m) old=",f8.4," new=",f8.4)') &
            iconc_tot_old, iconc_tot, ithk_tot_old, ithk_tot_new
      write(*,'(A)') trim(mesg)
    endif
  
  enddo

  if (allocated(sice)) deallocate(sice)
  if (allocated(tfi)) deallocate(tfi)
  if (allocated(data_in)) deallocate(data_in)

end subroutine apply_isponge

!> Convert global indices (itestG,jtestG) to indices on current tile 
subroutine global_to_local_ij(G, itestG, jtestG, itest, jtest)
  type(SIS_hor_grid_type), intent(in) :: G          !< The horizontal grid type
  integer, intent(in) :: itestG, jtestG
  integer, intent(out) :: itest, jtest

  integer :: current_pe, nihalo, njhalo, iscG, iecG, jscG, jecG
  integer :: isdG, jsdG, iedG, jedG
  integer :: nic, njc, nid, njd

  character(len=50) :: mdl  ! subroutine name
  character(len=256) :: mesg

  mdl = 'global_to_local_ij'
  nihalo = G%Domain%nihalo
  njhalo = G%Domain%njhalo

  current_pe = PE_here()

  ! Exclude halo points, computational domain:
  nic = G%iec - G%isc + 1
  njc = G%jec - G%jsc + 1
  iscG = G%isd_global + nihalo; iecG = iscG + nic
  jscG = G%jsd_global + njhalo; jecG = jscG + njc

  ! Data domain:
  nid  = G%ied - G%isd + 1
  njd  = G%jed - G%jsd + 1
  isdG = G%isd_global; iedG = isdG + nid
  jsdG = G%jsd_global; jedG = jsdG + njd

  ! Find test point:
  itest = 0; jtest = 0
  if (iscG <= itestG .and. itestG <= iecG .and. jscG <= jtestG .and. jtestG <= jecG) then
    itest = itestG - isdG + 1; jtest = jtestG - jsdG + 1
  endif
                 
  if (itest > 0 .and. jtest > 0) then
    write(mesg, '(A," PE=",i5," Test pnt Global i, j=", 2(i5,1x)," local i, j=", 2(i5,1x), &
                  "isdG/iedG=", 2(i5,1x), "jsdG/jedG=", 2(i5,1x))') &
         trim(mdl), current_pe, itestG, jtestG, itest, jtest, isdG, iedG, jsdG, jedG 
    write(*, '(A)') trim(mesg)
  endif        

end subroutine global_to_local_ij

! Map local indices to global
subroutine local_to_global_indx(G, i, j, iiG, jjG)
  type(SIS_hor_grid_type),   intent(in)  :: G          !< The horizontal grid type
  integer, intent(in)                    :: i, j 
  integer, intent(out)                   :: iiG, jjG

  iiG = G%isd_global + (i-1)  ; jjG = G%jsd_global + (j-1)

end subroutine local_to_global_indx

!> Redistribute input target 2D hice and iconc into ice thickness categories
!! place all ice into 1 category based on original ice thickness (hice)
!! Fill the lower cats with "some" ice to avoid 0s ice allowing better convergence
!! of the Icepack ITD iteration algorithm
subroutine distribute_ice2cats(CS, IG, G, scaled, eps_err)
  type(isponge_CS),        pointer     :: CS       !< A pointer that is set to point to the ice sponge control
                                                   !! structure for this module
  type(ice_grid_type),     intent(in)  :: IG       !< The sea-ice specific grid type
  type(SIS_hor_grid_type), intent(in)  :: G        !< The horizontal grid type
  logical, optional,       intent(in)  :: scaled   !< true if input hice (m) converted to kg/m2 and scaled
                                                   !! default = .false. hice units: (m) = m3/m2 vol/m2
  real, optional,          intent(in)  :: eps_err  !< error allowed for hice, cice after redistribution 

  ! local variables
  integer :: isd, ied, jsd, jed      !< data domain indices
  integer :: isdG, jsdG              !< Global indx, start pnts. data domain
  integer :: m, i, j, k, col
  integer :: iiG, jjG
  integer :: CatIce                  !< number of sea ice categories.
  integer :: icat0                   !< cat where to assign original hice and cice to begin distr. 

  real, allocatable, dimension(:,:) :: cice2d, hice2d
  real, allocatable, dimension(:) :: hLim_vals
  real, allocatable, dimension(:) :: hcat, ccat   !< 1D arrays for thkn and conc distributed by cats
  real, allocatable, dimension(:) :: volcat       !< ice volume (m3/m2) by cats = hcat*ccat, sum(volcat)=hice
  real :: Iscale              !< inverse scale to "unscale" the data
  real :: scale_cf
  real :: hice, cice          !< relax total ice thkn (vol m3/m2) and conc at a grid pnt
  real :: hice_k              !< hice_cat ice thkn in a category k, i.e. scaled by cice(cat=k)
  real :: hice_tot, cice_tot  !< total (sum over cats) ice thickn (volume, i.e. thkn*conc) and conc 
  real :: eps0                !< error allowed for hice, cice after redistribution
  real :: ck_min              !< min ice conc used in the lower cats to distribute hice/cice
  real :: ccat_k, hcat_k, dch_k !< conc, ivol, ice vol change in a cat=k
  real :: cnew 
  real :: htot_min
  real :: part_water          !< partial area of open water
  character(len=40)  :: mdl = "distribute_ice2cats"  ! This module's name.
  character(len=256) :: mesg
  logical :: scaled_hice
  logical :: err_hice, err_cice  !< checks for conserved total ice thkn and conc

  scaled_hice = .false.
  if (present(scaled)) scaled_hice = scaled
  eps0 = 1.e-10  !< 0-checking
  if (present(eps_err)) eps0 = eps_err

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isdG = G%isd_global;   jsdG = G%jsd_global
  CatIce = IG%CatIce

  allocate(cice2d(isd:ied,jsd:jed), hice2d(isd:ied,jsd:jed))
  allocate(hLim_vals(CatIce+1), ccat(CatIce), hcat(CatIce), volcat(CatIce))

  hLim_vals(:) = IG%cat_thick_lim(:)  !< ice thkn cats are not scaled (in m)

  cice2d=0.0 ; hice2d=0.0
  do m=1,CS%fldno
    scale_cf = CS%Ref_val(m)%scale
    Iscale = 1.0
    if (abs(1.-scale_cf).gt.1.e-10 .and. scale_cf.gt.0.) &
      Iscale = 1.0/scale_cf
    select case (trim(CS%var(m)%fld_name))
      case('mH_ice')
        hice2d = CS%Ref_orig(m)%fld
        if (scaled_hice .and. abs(1.-scale_cf).gt.1.e-10) &
            hice2d = CS%Ref_orig(m)%fld*Iscale     !< unscale input hice to original units (m) to find ice cat 
      case('part_size')
        cice2d = CS%Ref_orig(m)%fld                !< partial area (conc) is not scaled
    end select
  enddo

  do col=1,CS%num_col
    i = CS%col_i(col) ; j = CS%col_j(col)
    hice = hice2d(i,j)      !< in input units (m), note this is "volume/m2", i.e. (ice thkn)*(ice conc)
    cice = cice2d(i,j)      !< total partial area (conc)
    !if (hice.lt.1.e-10 .or. cice.lt.1.e-10) cycle  ! input ice fields are 0s
    if (hice.lt.1.e-10 .or. cice.lt.1e-10) then
      hice = 0.0 ; cice = 0.0
      do m=1,CS%fldno ; do k=1,CatIce
        select case (trim(CS%var(m)%fld_name))
          case('mH_ice')
            CS%Ref_val(m)%p(col,k) = IG%mH_cat_bound(k+1) !< The lower mass/unit area limits for ice cat [R Z ~> kg m-2].
          case('part_size')
            CS%Ref_val(m)%p(col,k) = 0.0
            if (k.eq.1) CS%Ref_val(m)%p(col,0) = 1.0     !< open water partial area
        end select
      enddo ; enddo
      cycle
    endif    

    ck_min = 1.e-2    !< conc in the lower cats, some small value but >> eps0 
    icat0 = 1e6
    hice_k = hice/cice   ! ice thkn (vol/part. area)  from ice volume (m3/m2) in cat=k
    icat0 = find_icat(hice_k, CatIce, hLim_vals)
! Check for ice cat. error:
    if (icat0<1 .or. icat0>CatIce) then
      iiG = isdG + (i-1) ; jjG = jsdG + (j-1)
      write(mesg,'(A,"ERROR: iG,jG=",2(i4,1x)," hice=",D16.4," cice=",D16.4," hice_k=",D16.4)') &
             trim(mdl), iiG, jjG, hice, cice, hice_k
      write(*,'(A)') trim(mesg)
      print*,"isnan cice=",is_nan(cice), "isnan hice=", is_nan(hice)
      print*,"hice.lt.1.e-10:",(hice.lt.1.e-10)," cice.lt.1.e-10:",(cice.lt.1.e-10)
      print*,"cice=cice",(cice.eq.cice)
      write(mesg,'(A," error: icat0 ",i2," hice=",f12.4," cice=",f12.4," hice_k=",f12.6)') &
            trim(mdl), hice, cice, hice_k
      call SIS_error(FATAL, trim(mesg)) 
    endif

    ! Adjust min conc in cats for low partial areas
    if (cice < ck_min*float(icat0)) ck_min=cice/(float(icat0))

! Check if there is enough ice for distribution:
    htot_min = sum(ck_min*hLim_vals(1:icat0))
    if (hice.lt.htot_min .or. cice.lt.(ck_min*icat0) .or. (hice*cice).lt.1.e-10) then
      do m=1,CS%fldno ; do k=1,CatIce
        select case (trim(CS%var(m)%fld_name))
          case('mH_ice')
            CS%Ref_val(m)%p(col,k) = 0.0 
          case('part_size')
            CS%Ref_val(m)%p(col,k) = 0.0
            if (k.eq.1) CS%Ref_val(m)%p(col,0) = 1.0  ! open water fraction size
        end select
      enddo ; enddo
      cycle
    endif
! distribute ice thicknesses and conc by cats
    hcat = 0.0 ; ccat = 0.0 ; volcat = 0.0
    hcat(icat0) = hice_k ; ccat(icat0) = cice ; volcat(icat0) = hice
    iiG = isdG + (i-1) ; jjG = jsdG + (j-1)

    call check_hcice(CatIce, hcat, ccat, hice, cice, err_hice, err_cice, &
                     hice_tot, cice_tot, iiG, jjG, str='1.', verb=.true.)

    do k=1,icat0-1
      ccat_k = ck_min
      hcat_k = hLim_vals(k) + eps0  
      dch_k = ccat_k*hcat_k   !< ice vol moved into this cat
      if (dch_k.gt.volcat(icat0)) exit !< not enough ice left in the thickest cat         
      ccat(k) = ccat_k  
      hcat(k) = hcat_k
      volcat(k) = dch_k

      ! Update thickest cat with initial ice:
      cnew = ccat(icat0)-ccat(k)
      cnew = max(cnew, ck_min)
      ccat(icat0) = cnew
      volcat(icat0) = volcat(icat0)-dch_k
      hcat(icat0) = volcat(icat0)/ccat(icat0)

    enddo
 
    call check_hcice(CatIce, hcat, ccat, hice, cice, err_hice, err_cice, &
                     hice_tot, cice_tot, iiG, jjG, str='2.', verb=.true.)

! Diagnostics at the test point:
    if (i.eq.CS%itest .and. j.eq.CS%jtest) then
      write(mesg,'(A," test pnt: hice=",f7.3," cice=",f6.3," htot=",f7.3,&
                  " ctot=",f7.3)') &
            trim(mdl), hice, cice, hice_tot, cice_tot
      write(*,'(A)') trim(mesg)
    endif

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
    do m=1,CS%fldno
      select case (trim(CS%var(m)%fld_name))
        case('mH_ice')               
          hcat(1:CatIce) = CS%Ref_val(m)%p(col,1:CatIce)*Iscale
        case('part_size')
          ccat(1:CatIce) = CS%Ref_val(m)%p(col,1:CatIce)
      end select
    enddo

    call partial_area_total(CatIce, ccat(1:CatIce), cice_tot)
    call ice_thkn_total(CatIce, ccat(1:CatIce), hcat, hice_tot)

    if (abs(cice_tot - cice).gt.eps0) then
      write(mesg,'(A,"ice conc. not conserved: init=",f6.3," after redistr.=",f6.3," err=",d14.4)') &
            cice, cice_tot, abs(cice_tot - cice)
      call SIS_error(WARNING, trim(mesg))
    endif

    if (abs(hice_tot - hice).gt.eps0) then
      write(mesg,'(A,"ice thkn  not conserved: init=",f6.3," after redistr.=",f6.3," err=",d14.4)') &
            hice, hice_tot, abs(hice_tot - hice)
      call SIS_error(WARNING, trim(mesg))
    endif

  enddo   !< do col

  if (allocated(hice2d)) deallocate(hice2d)
  if (allocated(cice2d)) deallocate(cice2d)
  if (allocated(hLim_vals)) deallocate(hLim_vals)
  if (allocated(ccat)) deallocate(ccat)
  if (allocated(hcat)) deallocate(hcat)
  if (allocated(volcat)) deallocate(volcat)

end subroutine distribute_ice2cats

!< Check if total ice thkn*conc and conc are conserved (i.e. equal original hice, cice)
subroutine check_hcice(CatIce, hcat, ccat, hice, cice, err_hice, err_cice, &
                       hice_tot, cice_tot, iG, jG, str, verb)
  integer, intent(in) :: CatIce
  real, dimension(CatIce), intent(in) :: hcat
  real, dimension(CatIce), intent(in) :: ccat
  real, intent(in) :: hice
  real, intent(in) :: cice
  real, intent(inout)    :: hice_tot, cice_tot
  logical, intent(inout) :: err_hice, err_cice
  logical,           optional, intent(in) :: verb
  character(len=*),  optional, intent(in) :: str

  real :: eps_err
  integer :: iG, jG
  logical :: verbose
  character(len=40) :: msg_info
  character(len=200) :: mesg

  err_hice = .false. ; err_cice = .false.
  eps_err = 1.e-6      !< allow small error during distribution
  hice_tot = sum(hcat*ccat)  !< Total ice volume m3/m2
  cice_tot = sum(ccat)
 
  if (abs(hice_tot-hice).gt.eps_err) err_hice=.true. 
  if (abs(cice_tot-cice).gt.eps_err) err_cice=.true.

  verbose=.false.
  if (present(verb)) verbose=verb

  if (.not.verbose) return

! Error information:
  msg_info = "check hice cice "
  if (present(str)) msg_info = trim(msg_info)//str
  if (err_hice) then
    write(mesg,'(A," iG, jG=",2(i4,1x)," hice not conserved: ",f12.6," hice=",f12.6)') &
          trim(msg_info), iG, jG, hice_tot, hice  
    write(*,'(A)') trim(mesg)
  elseif (err_cice) then
    write(mesg,'(A," iG, jG=",2(i4,1x)," cice not conserved: ",f12.6," cice=",f12.6)') &
          trim(msg_info), iG, jG, cice_tot, cice  
    write(*,'(A)') trim(mesg)
  endif

end subroutine check_hcice

function find_icat(hice_k, CatIce, hLim_vals) result (icat0)
  integer :: icat0                           !< The ice thkn category where hice_k belongs
  integer, intent(in)          :: CatIce     !< The number of ice thkn cats
  real, dimension(CatIce+1), &
                    intent(in) :: hLim_vals  !< ice thkn cats, not scaled (m)
  real, intent(in)             :: hice_k     !< ice thkn in a category (i.e. voli/cice, voli=hice*cice m3/m2

  integer :: k 

  icat0 = 1e6
  if (hice_k .ge. hLim_vals(CatIce)) then
    icat0 = CatIce
  elseif (hice_k .lt. hLim_vals(1)) then
    icat0=1
  else
    do k=1,CatIce
      if (hice_k .ge. hLim_vals(k) .and. hice_k .lt. hLim_vals(k+1)) then
        icat0 = k
        exit
      endif
    enddo
  endif

end function find_icat

!< Subroutine computes total partial area for 1D array of cice(1:CatIce) partial areas by cats.
subroutine partial_area_total(CatIce, cice_cat, cice_tot)
  integer, intent(in) :: CatIce
  real, dimension(CatIce), intent(in) :: cice_cat  !< 1D array of partial areas by cats.
  real, intent(inout) :: cice_tot                  !< total ice partial area of the grid cell

  integer :: k

  cice_tot = 0.0
  do k=1,CatIce
    cice_tot = cice_tot + cice_cat(k)
  enddo

end subroutine partial_area_total

!< Subroutine computes grid cell mean ice thickness for 1D arrays of thikn and partial area by cats.
subroutine ice_thkn_total(CatIce, cice_cat, hice_cat, hice_tot)
  integer, intent(in) :: CatIce
  real, dimension(CatIce), intent(in) :: cice_cat  !< 1D array of partial areas by cats.
  real, dimension(CatIce), intent(in) :: hice_cat  !< 1D array of ice thickn. by cats.
  real, intent(inout) :: hice_tot                  !< grid cell mean ice thickness

  integer :: k
  real :: ithkn_tot, iconc_tot

  hice_tot = 0.0
  do k=1,CatIce
    hice_tot = hice_tot + cice_cat(k)*hice_cat(k)
  enddo

end subroutine ice_thkn_total

!
! For debugging: print ice thkn by cats and total at 1 point
! derive the fields from isponge CS pointer structure or 
! Directly from the control structure IST (target for isponge CS)
subroutine print_ice_thkn_conc(IST, CS, G, IG, US, mesg_in, use_IST)
  type(isponge_CS),        pointer     :: CS      !< A pointer that is set to point to the ice sponge control
                                                  !! structure for this module
  type(SIS_hor_grid_type), intent(in)  :: G       !< The horizontal grid type
  type(ice_grid_type),     intent(in)  :: IG      !< The sea-ice specific grid type
  type(ice_state_type),    intent(in)  :: IST     !< A type describing the state of the sea ice
  type(unit_scale_type),   intent(in)  :: US      !< A structure with unit conversion factors
  character(len=*), optional, &
                           intent(in)  :: mesg_in
  logical, optional, intent(in) :: use_IST
  
  integer            :: CatIce, col, k, m, iiG, jjG, i, j
  real               :: rho_ice, mHice_k, iconc_k, ithkn_k, ithkn_tot, iconc_tot
  character(len=40)  :: str_in
  character(len=256) :: mesg
  character(len=15)  :: cs_name
  logical :: fld_IST

  fld_IST = .false.
  if (present(use_IST)) fld_IST = use_IST

  if (fld_IST) then
    cs_name = 'IST'
  else
    cs_name = 'ispongeCS'
  endif

  str_in = ' '
  if (present(mesg_in)) str_in = trim(mesg_in)
  str_in = trim(str_in) // trim(cs_name)

  CatIce = IG%CatIce
  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice)
  !H_to_m_ice = US%Z_to_m / rho_ice  ! convert kg/m2 scaled to m of ice (unscaled)

  do col=1,CS%num_col
    i = CS%col_i(col) ; j = CS%col_j(col)
    if (i.ne.CS%itest .or. j.ne.CS%jtest) cycle

    iiG = G%isd_global + (i-1)  ; jjG = G%jsd_global + (j-1)
    do k=1,IG%CatIce
      do m=1,CS%fldno
        select case (trim(CS%var(m)%fld_name))
          case('mH_ice') 
            if (fld_IST) then
              mHice_k = IST%mH_ice(i,j,k) * US%RZ_to_kg_m2  ! thkn in kg/m2
            else
              mHice_k = CS%var(m)%p(i,j,k) * US%RZ_to_kg_m2  ! thkn in kg/m2
            endif
          case('part_size') 
            if (fld_IST) then
              iconc_k = IST%part_size(i,j,k)
            else
              iconc_k = CS%var(m)%p(i,j,k)  ! partial area in cat
            endif
          case default
            write(mesg,'("SIS_sponge: Unknown relaxation field: ",A)') trim(CS%var(m)%fld_name)
            call SIS_error(FATAL,"print_ice_thkn_conc: "//mesg)
        end select
      enddo
      ithkn_k = mHice_k * iconc_k / rho_ice  ! vol/m2 = m3/m2 = [m], thickness
      ithkn_tot = ithkn_tot + ithkn_k
      iconc_tot = iconc_tot + iconc_k

      write(mesg,'(A," Test iG/jG=",2(i4,1x)," k=",i2," mHice=",f12.5," thkn=",f12.5," conc=",f12.5)') &
          trim(str_in), iiG, jjG, k, mHice_k, ithkn_k, iconc_k
      write(*,'(A)') trim(mesg)
    enddo

    write(mesg,'(A," ====  TOTAL: thkn=",f12.5," conc=",f12.5)') trim(str_in), ithkn_tot, iconc_tot
    write(*,'(A)') trim(mesg) 
    exit

  enddo

end subroutine print_ice_thkn_conc

!> Deallocate memory associated with the SIS_sponge module
subroutine SIS_sponge_end(CS)
  type(isponge_CS), pointer :: CS !< The ice sponge control structure that is deallocated here
  
  deallocate(CS) 
                 
end subroutine SIS_sponge_end

end module SIS_sponge 

