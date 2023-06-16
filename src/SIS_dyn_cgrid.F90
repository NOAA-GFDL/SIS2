!> Update sea-ice dynamics using elastic-viscous-plastic rheology with a C-grid discretization
module SIS_dyn_cgrid

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!                                                                              !
! C-grid SEA ICE DYNAMICS using ELASTIC-VISCOUS-PLASTIC RHEOLOGY adapted from  !
! Hunke and Dukowicz (JPO 1997, H&D hereafter) with some derivation from SIS1  !
! and with guidance from the C-grid implementations of sea-ice in MITgcm as    !
! documented in MITgcm user notes by Martin Losch and in LIM3 by S. Bouillon   !
! et al. (Ocean Modelling, 2009 & 2013). This code initially written by        !
! Robert Hallberg in 2013.                                                     !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use ice_grid,          only : ice_grid_type

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, NOTE, SIS_mesg=>MOM_mesg
use MOM_file_parser,   only : get_param, log_param, read_param, log_version, param_file_type
use MOM_domains,       only : pass_var, pass_vector, CGRID_NE, CORNER, pe_here
use MOM_domains,       only : MOM_domain_type, clone_MOM_domain
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : open_ASCII_file, APPEND_FILE, ASCII_FILE, MULTIPLE, SINGLE_FILE
use MOM_io,            only : MOM_read_data
use MOM_time_manager,  only : time_type, real_to_time, operator(+), operator(-)
use MOM_time_manager,  only : set_date, get_time, get_date
use MOM_unit_scaling,  only : unit_scale_type

use SIS_diag_mediator, only : post_SIS_data, SIS_diag_ctrl
use SIS_diag_mediator, only : query_SIS_averaging_enabled, enable_SIS_averaging
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_debugging,     only : chksum, Bchksum, hchksum, uvchksum
use SIS_debugging,     only : check_redundant_B, check_redundant_C
use SIS_restart,       only : register_restart_field, only_read_from_restarts, SIS_restart_CS
use SIS_restart,       only : query_initialized=>query_inited
use SIS_framework,     only : safe_alloc
use SIS_hor_grid,      only : SIS_hor_grid_type
use SIS_types,         only : ice_state_type
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_C_dyn_init, SIS_C_dynamics, SIS_C_dyn_end
public :: SIS_C_dyn_register_restarts, SIS_C_dyn_read_alt_restarts
public :: basal_stress_coeff_C, basal_stress_coeff_itd

!> The control structure with parameters regulating C-grid ice dynamics
type, public :: SIS_C_dyn_CS ; private
  real, allocatable, dimension(:,:) :: &
    str_t, &  !< The tension stress tensor component [R Z L2 T-2 ~> Pa m].
    str_d, &  !< The divergence stress tensor component [R Z L2 T-2 ~> Pa m].
    str_s     !< The shearing stress tensor component (cross term) [R Z L2 T-2 ~> Pa m].

  ! parameters for calculating water drag and internal ice stresses
  real :: p0                  !< Pressure constant in the Hibler rheology [R L2 T-2 ~> Pa]
  real :: p0_rho              !< The pressure constant divided by ice density [L2 T-2 ~> N m kg-1].
  real :: c0                  !< another pressure constant, c* in Hunke & Dukowicz 1997 [nondim]
  real :: cdw                 !< The drag coefficient between the sea ice and water. [nondim]
  real :: EC = 2.0            !< yield curve axis ratio [nondim]
  real :: Rho_ocean           !< The nominal density of sea water [R ~> kg m-3].
  real :: Rho_ice             !< The nominal density of sea ice [R ~> kg m-3].
  real :: drag_bg_vel2 = 0.0  !< A background (subgridscale) velocity for drag with the ocean
                              !< squared [L2 T-2 ~> m2 s-2].  This is always 0 for now.
  real :: min_ocn_inertial_h  !< A minimum ocean thickness used to limit the viscous coupling
                              !! rate implied for the ocean by the ice-ocean stress [Z ~> m].
  real :: Tdamp               !< The damping timescale of the stress tensor components toward their
                              !! equilibrium solution due to the elastic terms [T ~> s] or [nondim].
  real :: del_sh_min_scale    !< A scaling factor for the minimum permitted value of minimum
                              !! shears used in the denominator of the stress equations [nondim].
                              !  I suspect that this needs to be greater than 1. -RWH
  real    :: vel_underflow    !< Velocity components smaller than vel_underflow
                              !! are set to 0 [L T-1 ~> m s-1].
  real    :: str_underflow    !< Stress tensor components smaller than str_underflow
                              !! are set to 0 [R Z L2 T-2 ~> Pa m].
  real    :: CFL_trunc        !< Velocity components will be truncated when they are large enough
                              !! that the corresponding CFL number exceeds this value [nondim].
  logical :: CFL_check_its    !< If true, check the CFL number for every iteration
                              !! of the rheology solver; otherwise only check the
                              !! final velocities that are used for transport.
  logical :: debug            !< If true, write verbose checksums for debugging purposes.
  logical :: debug_EVP        !< If true, write out verbose debugging data for each of
                              !! the steps within the EVP solver.
  logical :: debug_redundant  !< If true, debug redundant points.
  logical :: project_drag_vel !< If true, project forward the ice velocity used in the drag
                              !! calculation to avoid an instability that can occur when an finite
                              !! stress is applied to thin ice moving with the velocity of the ocean.
  logical :: project_ci       !< If true, project the ice concentration and related ice strength
                              !! changes due to the convergent or divergent ice flow.
  logical :: weak_coast_stress = .false. !< If true, do not use land masks in determining the area
                              !! for stress convergence, which acts to weaken the stress-driven
                              !! acceleration in coastal points.
  logical :: weak_low_shear = .false. !< If true, the divergent stresses go toward 0 in the C-grid
                              !! dynamics when the shear magnitudes are very weak.
                              !! Otherwise they go to -P_ice.  This setting is temporary.
  integer :: evp_sub_steps    !< The number of iterations in the EVP dynamics
                              !! for each slow time step.
  real    :: dt_Rheo          !< The maximum sub-cycling time step for the EVP dynamics [T ~> s].
  type(time_type), pointer :: Time => NULL() !< A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the
                              !! timing of diagnostic output.
  integer, pointer :: ntrunc => NULL() !< The number of times the velocity has been truncated
                              !! since the last call to write_ice_statistics.
  character(len = 200) :: u_trunc_file !< The complete path to the file in which a column's worth
                              !! of u-accelerations are written if velocity truncations occur.
  character(len = 200) :: v_trunc_file !< The complete path to the file in which a column's worth
                              !! of v-accelerations are written if velocity truncations occur.
  integer :: u_file = -1      !< The unit number for an opened u-truncation file, or -1 if it has
                              !! not been opened.
  integer :: v_file = -1      !< The unit number for an opened v-truncation file, or -1 if it has
                              !! not been opened.
  integer :: cols_written     !< The number of columns whose output has been
                              !! written by this PE during the current run.
  integer :: max_writes       !< The maximum number of times any PE can write out
                              !! a column's worth of accelerations during a run.
  logical :: lemieux_landfast !< If true, use the Lsemieux landfast ice parameterization.
  real :: lemieux_k1          !< 1st free parameter for landfast parameterization [nondim]
  real :: lemieux_k2          !< second free parameter (N/m^3) for landfast parametrization [R L T-2 ~> N m-3]
  real :: lemieux_alphab      !< Cb factor in Lemieux et al 2015 [nondim]
  real :: lemieux_threshold_hw !< max water depth for grounding [Z ~> m]
                              !! see keel data from Amundrud et al. 2004 (JGR)
  real :: lemieux_u0          !< residual velocity for basal stress [L T-1 ~> m s-1]
  logical :: itd_landfast     !< If true, use the probabilistic landfast ice parameterization.
  real :: basal_stress_min_thick !< min ice thickness for grounding [Z ~> m]
  real :: basal_stress_max_depth !< max water depth for grounding [Z ~> m]
  real :: basal_stress_mu_s   !< bottom drag parameter [L Z-1 ~> nondim]
  real :: bathy_roughness_min !< minimum bathymetric roughness [z ~> m]
  real :: bathy_roughness_max !< maximum bathymetric roughness [z ~> m]
  real :: puny                !< small number [nondim]
  real :: onemeter            !< make the units work out (hopefully) [Z ~> m]
  real :: basal_stress_cutoff !< tunable parameter for the bottom drag [nondim]
  integer :: ncat_b           ! number of bathymetry categories
  integer :: ncat_i           ! number of ice thickness categories (log-normal)

  real, pointer, dimension(:,:) :: Tb_u=>NULL() !< Basal stress component at u-points
                                                !! [R Z L T-2 -> kg m-1 s-2]
  real, pointer, dimension(:,:) :: Tb_v=>NULL() !< Basal stress component at v-points
                                                !! [R Z L T-2 -> kg m-1 s-2]
  real, pointer, dimension(:,:) :: sigma_b=>NULL()   !< !< Bottom depth variance [Z ~> m].

  logical :: FirstCall = .true. !< If true, this module has not been called before
  !>@{ Diagnostic IDs
  integer :: id_fix = -1, id_fiy = -1, id_fcx = -1, id_fcy = -1
  integer :: id_fwx = -1, id_fwy = -1, id_sigi = -1, id_sigii = -1
  integer :: id_flfx = -1, id_flfy = -1, id_stren = -1, id_stren0 = -1
  integer :: id_ui = -1, id_vi = -1, id_Coru = -1, id_Corv = -1
  integer :: id_PFu = -1, id_PFv = -1, id_fpx = -1, id_fpy = -1
  integer :: id_fix_d = -1, id_fix_t = -1, id_fix_s = -1
  integer :: id_fiy_d = -1, id_fiy_t = -1, id_fiy_s = -1
  integer :: id_str_d = -1, id_str_t = -1, id_str_s = -1
  integer :: id_sh_d = -1, id_sh_t = -1, id_sh_s = -1
  integer :: id_del_sh = -1, id_del_sh_min = -1
  integer :: id_mis = -1, id_ci = -1, id_ci0 = -1, id_miu = -1, id_miv = -1
  integer :: id_ui_hifreq = -1, id_vi_hifreq = -1
  integer :: id_str_d_hifreq = -1, id_str_t_hifreq = -1, id_str_s_hifreq = -1
  integer :: id_sh_d_hifreq = -1, id_sh_t_hifreq = -1, id_sh_s_hifreq = -1
  integer :: id_sigi_hifreq = -1, id_sigii_hifreq = -1
  integer :: id_stren_hifreq = -1, id_ci_hifreq = -1
  integer :: id_siu = -1, id_siv = -1, id_sispeed = -1 ! SIMIP diagnostics
  !!@}
end type SIS_C_dyn_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_C_dyn_init initializes the ice dynamics, sets parameters, and registers diagnostics
subroutine SIS_C_dyn_init(Time, G, US, param_file, diag, CS, ntrunc)
  type(time_type),     target, intent(in)    :: Time !< The sea-ice model's clock,
                                                     !! set with the current model time.
  type(SIS_hor_grid_type),     intent(in)    :: G    !< The horizontal grid type
  type(unit_scale_type),       intent(in)    :: US    !< A structure with unit conversion factors
  type(param_file_type),       intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(SIS_diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate diagnostic output
  type(SIS_C_dyn_CS),          pointer       :: CS   !< The control structure for this module
  integer,   target, optional, intent(inout) :: ntrunc !< The integer that stores the number of times
                                                     !! the velocity has been truncated since the
                                                     !! last call to write_ice_statistics.

!   This subroutine sets the parameters and registers the diagnostics associated
! with the ice dynamics.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40) :: mdl = "SIS_C_dyn" ! This module's name.
  character(len=200) :: filename, h2_file, inputdir
  logical           :: debug
  real, parameter   :: missing = -1e34

  if (.not.associated(CS)) then
    call SIS_error(FATAL, "SIS_C_dyn_init called with an unassociated control structure. \n"//&
                    "SIS_C_dyn_register_restarts must be called before SIS_C_dyn_init.")
    return
  endif

  CS%diag => diag
  CS%Time => Time
  if (present(ntrunc)) then ; CS%ntrunc => ntrunc ; else ; allocate(CS%ntrunc) ; endif
  CS%ntrunc = 0

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version)
  call get_param(param_file, mdl, "DT_RHEOLOGY", CS%dt_Rheo, &
                 "The sub-cycling time step for iterating the rheology "//&
                 "and ice momentum equations. If DT_RHEOLOGY is negative, "//&
                 "the time step is set via NSTEPS_DYN.", units="seconds", &
                 default=-1.0, scale=US%s_to_T)
  CS%evp_sub_steps = -1
  if (CS%dt_Rheo <= 0.0) &
    call get_param(param_file, mdl, "NSTEPS_DYN", CS%evp_sub_steps, &
                 "The number of iterations of the rheology and ice "//&
                 "momentum equations for each slow ice time step.", default=432)
  call get_param(param_file, mdl, "ICE_TDAMP_ELASTIC", CS%Tdamp, &
                 "The damping timescale associated with the elastic terms "//&
                 "in the sea-ice dynamics equations (if positive) or the "//&
                 "fraction of DT_ICE_DYNAMICS (if negative).", &
                 units="s or nondim", default=-0.2)
  if (CS%Tdamp > 0.0) CS%Tdamp = CS%Tdamp*US%s_to_T
  call get_param(param_file, mdl, "WEAK_LOW_SHEAR_ICE", CS%weak_low_shear, &
                 "If true, the divergent stresses go toward 0 in the C-grid "//&
                 "dynamics when the shear magnitudes are very weak. "//&
                 "Otherwise they go to -P_ice.  This setting is temporary.", &
                 default=.false.)

  call get_param(param_file, mdl, "PROJECT_ICE_DRAG_VEL", CS%project_drag_vel, &
                 "If true, project forward the ice velocity used in the "//&
                 "drag calculation to avoid an instability that can occur "//&
                 "when a finite stress is applied to thin ice moving with "//&
                 "the velocity of the ocean.", default=.true.)
  call get_param(param_file, mdl, "ICE_YIELD_ELLIPTICITY", CS%EC, &
                 "The ellipticity coefficient for the plastic yield curve "//&
                 "in the sea-ice rheology.  For an infinite ellipticity "//&
                 "(i.e., a cavitating fluid rheology), use 0.", &
                 units="Nondim", default=2.0)

  call get_param(param_file, mdl, "ICE_STRENGTH_PSTAR", CS%p0, &
                 "A constant in the expression for the ice strength, "//&
                 "P* in Hunke & Dukowicz 1997.", &
                 units="Pa", default=2.75e4, scale=US%kg_m3_to_R*US%m_s_to_L_T**2)
  call get_param(param_file, mdl, "ICE_STRENGTH_CSTAR", CS%c0, &
                 "A constant in the exponent of the expression for the "//&
                 "ice strength, c* in Hunke & Dukowicz 1997.", &
                 units="nondim", default=20.)
  call get_param(param_file, mdl, "ICE_CDRAG_WATER", CS%cdw, &
                 "The drag coefficient between the sea ice and water.", &
                 units="nondim", default=3.24e-3)
  call get_param(param_file, mdl, "MIN_OCN_INTERTIAL_H", CS%min_ocn_inertial_h, &
                 "A minimum ocean thickness used to limit the viscous coupling rate "//&
                 "implied for the ocean by the ice-ocean stress. Only used if positive.", &
                 units="m", default=0.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "ICE_DEL_SH_MIN_SCALE", CS%del_sh_min_scale, &
                 "A scaling factor for the lower bound on the shear rates "//&
                 "used in the denominator of the stress calculation. This "//&
                 "probably needs to be greater than 1.", units="nondim", default=2.0)
  call get_param(param_file, mdl, "PROJECT_ICE_CONCENTRATION", CS%project_ci, &
                 "If true, project the evolution of the ice concentration "//&
                 "due to the convergence or divergence of the ice flow.", default=.true.)

  call get_param(param_file, mdl, "RHO_OCEAN", CS%Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "RHO_ICE", CS%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0, scale=US%kg_m3_to_R)
  CS%p0_rho = CS%p0 / CS%Rho_ice

  call get_param(param_file, mdl, "VEL_UNDERFLOW", CS%vel_underflow, &
                 "A negligibly small velocity magnitude below which velocity "//&
                 "components are set to 0.  A reasonable value might be 1e-30 m/s, "//&
                 "which is less than an Angstrom divided by the age of the universe.", &
                 units="m s-1", default=0.0, scale=US%m_s_to_L_T)
  call get_param(param_file, mdl, "STRESS_UNDERFLOW", CS%str_underflow, &
                 "A negligibly small magnitude below which ice stress tensor "//&
                 "components are set to 0.  A reasonable value might be "//&
                 "1e-15 kg m-1 s-1 times vel_underflow.", &
                 units="Pa m", default=0.0, scale=US%m_s_to_L_T**2*US%kg_m3_to_R*US%m_to_Z)
  call get_param(param_file, mdl, "CFL_TRUNCATE", CS%CFL_trunc, &
                 "The value of the CFL number that will cause ice velocity "//&
                 "components to be truncated; instability can occur past 0.5.", &
                 units="nondim", default=0.5)
  call get_param(param_file, mdl, "CFL_TRUNC_DYN_ITS", CS%CFL_check_its, &
                 "If true, check the CFL number for every iteration of the "//&
                 "rheology solver; otherwise only the final velocities that "//&
                 "are used for transport are checked.", &
                 default=.false.)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_SLOW_ICE", CS%debug, &
                 "If true, write out verbose debugging data on the slow ice PEs.", &
                 default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_EVP_SUBSTEPS", CS%debug_EVP, &
                 "If true, write out verbose debugging data for each of the "//&
                 "steps within the EVP solver.", default=debug, debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_REDUNDANT", CS%debug_redundant, &
                 "If true, debug redundant data points.", default=CS%debug, &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to the file where the accelerations "//&
                 "leading to zonal velocity truncations are written. "//&
                 "Leave this empty for efficiency if this diagnostic is "//&
                 "not needed.", default="", debuggingParam=.true.)
  call get_param(param_file, mdl, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to the file where the accelerations "//&
                 "leading to meridional velocity truncations are written. "//&
                 "Leave this empty for efficiency if this diagnostic is not needed.", &
                 default="", debuggingParam=.true.)
  call get_param(param_file, mdl, "MAX_TRUNC_FILE_SIZE_PER_PE", CS%max_writes, &
                 "The maximum number of colums of truncations that any PE "//&
                 "will write out during a run.", default=50, debuggingParam=.true.)
  call get_param(param_file, mdl, "LEMIEUX_LANDFAST", CS%lemieux_landfast, &
                 "If true, turn on Lemieux landfast ice parameterization.", default=.false.)
  if (CS%lemieux_landfast) then
    call get_param(param_file, mdl, "LEMIEUX_K1", CS%lemieux_k1, &
                   "The value of the first Lemieux landfast ice tuneable parameter.", &
                   units="Nondim", default=8.0)
    call get_param(param_file, mdl, "LEMIEUX_K2", CS%lemieux_k2, &
                   "The value of the second Lemieux landfast ice tuneable parameter.", &
                   units="N m-3", default=15.0, scale=US%kg_m3_to_R*US%m_s_to_L_T*US%T_to_s)
    call get_param(param_file, mdl, "LEMIEUX_THRESHOLD_HW", CS%lemieux_threshold_hw, &
                   "Maximum water depth for grounding in Lemieux landfast ice.", &
                   units="m", default=30.0, scale=US%m_to_Z)
  endif
  call get_param(param_file, mdl, "ITD_LANDFAST", CS%itd_landfast, &
                 "If true, turn on probabilistic landfast ice parameterization.", default=.false.)

  if (CS%lemieux_landfast .or. CS%itd_landfast) then
    call get_param(param_file, mdl, "LEMIEUX_ALPHA_B", CS%lemieux_alphab, &
                   "The value of a third Lemieux landfast ice tuneable parameter.", &
                   units="Nondim", default=20.0)
    call get_param(param_file, mdl, "LEMIEUX_U0", CS%lemieux_u0, &
                   "Velocity for Lemieux landfast ice.", &
                   units="m s-1", default=5.e-5, scale=US%m_s_to_L_T)

    allocate(CS%Tb_u(G%IsdB:G%IedB,G%jsd:G%jed), source=0.0)
    allocate(CS%Tb_v(G%isd:G%ied,G%JsdB:G%JedB), source=0.0)
  endif
  if (CS%itd_landfast) then
    call get_param(param_file, mdl, "BASAL_STRESS_MIN_THICK", CS%basal_stress_min_thick, &
                   "Minimum ice thickness for grounding in ITD landfast ice.", &
                   units="m", default=0.01, scale=US%m_to_Z)
    call get_param(param_file, mdl, "BASAL_STRESS_MAX_DEPTH", CS%basal_stress_max_depth, &
                   "Maximum water depth for grounding in ITD landfast ice.", &
                   units="m", default=50.0, scale=US%m_to_Z)
    call get_param(param_file, mdl, "BASAL_STRESS_MU_S", CS%basal_stress_mu_s, &
                   "Drag coefficient in ITD landfast ice.", &
                   units="nondim", default=0.1, scale=US%L_to_Z)
    call get_param(param_file, mdl, "BASAL_STRESS_PUNY", CS%puny, &
                   "Small number in ITD landfast ice.", &
                   units="nondim", default=1.0e-20)
    call get_param(param_file, mdl, "BASAL_STRESS_SCALE", CS%onemeter, &
                   "Scale factor in ITD landfast ice.", &
                   units="m", default=1.0, scale=US%m_to_Z)
    call get_param(param_file, mdl, "BASAL_STRESS_CUTOFF", CS%basal_stress_cutoff, &
                   "Scale factor in ITD landfast ice.", &
                   units="nondim", default=1.9430)
    call get_param(param_file, mdl, "H2_FILE", h2_file, &
                 "The path to the file containing the sub-grid-scale "//&
                 "topographic roughness amplitude with ITD_LANDFAST.", &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    filename = trim(inputdir) // "/" // trim(h2_file)
    allocate(CS%sigma_b(G%isd:G%ied,G%jsd:G%jed), source=0.0)
    call MOM_read_data(filename, 'h2', CS%sigma_b, G%domain, scale=US%m_to_Z**2)
    call get_param(param_file, mdl, "BATHY_ROUGHNESS_MIN", CS%bathy_roughness_min, &
                   "Minimum bathymetric roughness.", &
                   units="m", default=2.5, scale=US%m_to_Z)
    call get_param(param_file, mdl, "BATHY_ROUGHNESS_MAX", CS%bathy_roughness_max, &
                   "Maximum bathymetric roughness.", &
                   units="m", default=2.5, scale=US%m_to_Z)
    CS%sigma_b(:,:) = min(max(sqrt(CS%sigma_b(:,:)), CS%bathy_roughness_min), CS%bathy_roughness_max)
    call pass_var(CS%sigma_b, G%Domain)
    call get_param(param_file, mdl, "BASAL_STRESS_NCAT_B", CS%ncat_b, &
                   "Number of bathymetric depth categories in landfast ice computation.", &
                   default=100)
    call get_param(param_file, mdl, "BASAL_STRESS_NCAT_I", CS%ncat_i, &
                   "Number of ice thickness categories in landfast ice computation.", &
                   default=100)
  endif

!  if (len_trim(dirs%output_directory) > 0) then
!    if (len_trim(CS%u_trunc_file) > 0) &
!      CS%u_trunc_file = trim(dirs%output_directory)//trim(CS%u_trunc_file)
!    if (len_trim(CS%v_trunc_file) > 0) &
!      CS%v_trunc_file = trim(dirs%output_directory)//trim(CS%v_trunc_file)
!    call log_param(param_file, mdl, "output_dir/U_TRUNC_FILE", CS%u_trunc_file)
!    call log_param(param_file, mdl, "output_dir/V_TRUNC_FILE", CS%v_trunc_file)
!  endif
  CS%u_file = -1 ; CS%v_file = -1 ; CS%cols_written = 0

  CS%id_sigi  = register_diag_field('ice_model','SIGI' ,diag%axesT1, Time,     &
            'first stress invariant', 'none', missing_value=missing)
  CS%id_sigii = register_diag_field('ice_model','SIGII' ,diag%axesT1, Time,    &
            'second stress invariant', 'none', missing_value=missing)
  CS%id_stren = register_diag_field('ice_model','STRENGTH' ,diag%axesT1, Time, &
            'ice strength', 'Pa*m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, missing_value=missing)
  CS%id_stren0 = register_diag_field('ice_model','STREN_0' ,diag%axesT1, Time, &
            'ice strength at start of rheology', &
            'Pa*m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, missing_value=missing)
  CS%id_fix   = register_diag_field('ice_model', 'FI_X', diag%axesCu1, Time,   &
            'ice internal stress - x component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fiy   = register_diag_field('ice_model', 'FI_Y', diag%axesCv1, Time,   &
            'ice internal stress - y component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fcx   = register_diag_field('ice_model', 'FC_X', diag%axesCu1, Time,   &
            'Coriolis force - x component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fcy   = register_diag_field('ice_model', 'FC_Y', diag%axesCv1, Time,   &
            'Coriolis force - y component', 'Pa', conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_Coru   = register_diag_field('ice_model', 'Cor_ui', diag%axesCu1, Time,&
            'Coriolis ice acceleration - x component', &
            'm s-2', conversion=US%L_T_to_m_s*US%s_to_T, &
            missing_value=missing, interp_method='none')
  CS%id_Corv   = register_diag_field('ice_model', 'Cor_vi', diag%axesCv1, Time,&
            'Coriolis ice acceleration - y component', &
            'm s-2', conversion=US%L_T_to_m_s*US%s_to_T, &
            missing_value=missing, interp_method='none')
  CS%id_fpx   = register_diag_field('ice_model', 'FP_X', diag%axesCu1, Time,   &
            'Pressure force - x component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fpy   = register_diag_field('ice_model', 'FP_Y', diag%axesCv1, Time,   &
            'Pressure force - y component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_PFu   = register_diag_field('ice_model', 'Pfa_ui', diag%axesCu1, Time, &
            'Pressure-force ice acceleration - x component', &
            'm s-2',  conversion=US%L_T_to_m_s*US%s_to_T, &
            missing_value=missing, interp_method='none')
  CS%id_PFv   = register_diag_field('ice_model', 'Pfa_vi', diag%axesCv1, Time, &
            'Pressure-force ice acceleration - y component', &
            'm s-2',  conversion=US%L_T_to_m_s*US%s_to_T, &
            missing_value=missing,  interp_method='none')
  CS%id_fwx   = register_diag_field('ice_model', 'FW_X', diag%axesCu1, Time,   &
            'water stress on ice - x component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fwy   = register_diag_field('ice_model', 'FW_Y', diag%axesCv1, Time,   &
            'water stress on ice - y component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_flfx  = register_diag_field('ice_model', 'FLF_X', diag%axesCu1, Time,   &
            'land-fast bottom stress on ice - x component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_flfy  = register_diag_field('ice_model', 'FLF_Y', diag%axesCv1, Time,   &
            'land-fast bottom stress on ice - y component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_ui    = register_diag_field('ice_model', 'UI', diag%axesCu1, Time,     &
            'ice velocity - x component', 'm/s', missing_value=missing,        &
            interp_method='none', conversion=US%L_T_to_m_s)
  CS%id_vi    = register_diag_field('ice_model', 'VI', diag%axesCv1, Time,     &
            'ice velocity - y component', 'm/s', missing_value=missing,        &
            interp_method='none', conversion=US%L_T_to_m_s)
  CS%id_mis  = register_diag_field('ice_model', 'MIS_tot', diag%axesT1, Time,  &
            'Mass of ice and snow at t-points', 'kg m-2', conversion=US%RZ_to_kg_m2, missing_value=missing)
  CS%id_ci0  = register_diag_field('ice_model', 'CI_tot', diag%axesT1, Time,   &
            'Initial summed concentration of ice at t-points', 'nondim',       &
            missing_value=missing)
  CS%id_ci  = register_diag_field('ice_model', 'CI_proj', diag%axesT1, Time,   &
            'Projected summed concentration of ice at t-points', 'nondim',     &
            missing_value=missing)
  CS%id_miu = register_diag_field('ice_model', 'MI_U', diag%axesCu1, Time,   &
            'Mass of ice and snow at u-points', 'kg m-2', conversion=US%RZ_to_kg_m2, &
            missing_value=missing, interp_method='none')
  CS%id_miv = register_diag_field('ice_model', 'MI_V', diag%axesCv1, Time,   &
            'Mass of ice and snow at v-points', 'kg m-2', conversion=US%RZ_to_kg_m2, &
            missing_value=missing, interp_method='none')

  CS%id_fix_d   = register_diag_field('ice_model', 'FI_d_X', diag%axesCu1, Time,         &
            'ice divergence internal stress - x component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fiy_d   = register_diag_field('ice_model', 'FI_d_Y', diag%axesCv1, Time,         &
            'ice divergence internal stress - y component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fix_t   = register_diag_field('ice_model', 'FI_t_X', diag%axesCu1, Time,        &
            'ice tension internal stress - x component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fiy_t   = register_diag_field('ice_model', 'FI_t_Y', diag%axesCv1, Time,        &
            'ice tension internal stress - y component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fix_s   = register_diag_field('ice_model', 'FI_s_X', diag%axesCu1, Time,        &
            'ice shearing internal stress - x component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')
  CS%id_fiy_s   = register_diag_field('ice_model', 'FI_s_Y', diag%axesCv1, Time,        &
            'ice shearing internal stress - y component', &
            'Pa',  conversion=US%RZ_T_to_kg_m2s*US%L_T_to_m_s, &
            missing_value=missing, interp_method='none')

  CS%id_str_d   = register_diag_field('ice_model', 'str_d', diag%axesT1, Time, &
            'ice divergence internal stress', 'Pa m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, &
            missing_value=missing)
  CS%id_str_t   = register_diag_field('ice_model', 'str_t', diag%axesT1, Time, &
            'ice tension internal stress', 'Pa m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, &
            missing_value=missing)
  CS%id_str_s   = register_diag_field('ice_model', 'str_s', diag%axesB1, Time, &
            'ice shearing internal stress', 'Pa m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, &
            missing_value=missing)
  CS%id_sh_d   = register_diag_field('ice_model', 'sh_d', diag%axesT1, Time,   &
            'ice divergence strain rate', 's-1', conversion=US%s_to_T, missing_value=missing)
  CS%id_sh_t   = register_diag_field('ice_model', 'sh_t', diag%axesT1, Time,   &
            'ice tension strain rate', 's-1', conversion=US%s_to_T, missing_value=missing)
  CS%id_sh_s   = register_diag_field('ice_model', 'sh_s', diag%axesB1, Time,   &
            'ice shearing strain rate', 's-1', conversion=US%s_to_T, missing_value=missing)
  CS%id_del_sh = register_diag_field('ice_model', 'del_sh', diag%axesT1, Time, &
            'ice strain rate magnitude', 's-1', conversion=US%s_to_T, missing_value=missing)
  CS%id_del_sh_min = register_diag_field('ice_model', 'del_sh_min', diag%axesT1, Time, &
            'minimum ice strain rate magnitude', 's-1', conversion=US%s_to_T, missing_value=missing)

  CS%id_ui_hifreq = register_diag_field('ice_model', 'ui_hf', diag%axesCu1, Time, &
            'ice velocity - x component', 'm/s', missing_value=missing,        &
            interp_method='none', conversion=US%L_T_to_m_s)
  CS%id_vi_hifreq = register_diag_field('ice_model', 'vi_hf', diag%axesCv1, Time, &
            'ice velocity - y component', 'm/s', missing_value=missing,        &
            interp_method='none, conversion=US%L_T_to_m_s')
  CS%id_str_d_hifreq = register_diag_field('ice_model', 'str_d_hf', diag%axesT1, Time, &
            'ice divergence internal stress', 'Pa m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, &
            missing_value=missing)
  CS%id_str_t_hifreq = register_diag_field('ice_model', 'str_t_hf', diag%axesT1, Time, &
            'ice tension internal stress', 'Pa m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, &
            missing_value=missing)
  CS%id_str_s_hifreq = register_diag_field('ice_model', 'str_s_hf', diag%axesB1, Time, &
            'ice shearing internal stress', 'Pa m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, &
            missing_value=missing)
  CS%id_sh_d_hifreq = register_diag_field('ice_model', 'sh_d_hf', diag%axesT1, Time, &
            'ice divergence rate', 's-1', conversion=US%s_to_T, missing_value=missing)
  CS%id_sh_t_hifreq = register_diag_field('ice_model', 'sh_t_hf', diag%axesT1, Time, &
            'ice tension rate', 's-1', conversion=US%s_to_T, missing_value=missing)
  CS%id_sh_s_hifreq = register_diag_field('ice_model', 'sh_s_hf', diag%axesB1, Time, &
            'ice shearing rate', 's-1', conversion=US%s_to_T, missing_value=missing)
  CS%id_sigi_hifreq  = register_diag_field('ice_model','sigI_hf' ,diag%axesT1, Time, &
            'first stress invariant', 'none', missing_value=missing)
  CS%id_sigii_hifreq = register_diag_field('ice_model','sigII_hf' ,diag%axesT1, Time, &
            'second stress invariant', 'none', missing_value=missing)
  CS%id_ci_hifreq  = register_diag_field('ice_model', 'CI_hf', diag%axesT1, Time, &
            'Summed concentration of ice at t-points', 'nondim', missing_value=missing)
  CS%id_stren_hifreq = register_diag_field('ice_model','STRENGTH_hf' ,diag%axesT1, Time, &
            'ice strength', 'Pa*m', conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2, missing_value=missing)

  CS%id_siu = register_diag_field('ice_model', 'siu', diag%axesT1, Time, &
            'ice velocity - x component', 'm/s', missing_value=missing,  &
            interp_method='none', conversion=US%L_T_to_m_s)
  CS%id_siv = register_diag_field('ice_model', 'siv', diag%axesT1, Time, &
            'ice velocity - y component', 'm/s', missing_value=missing,  &
            interp_method='none', conversion=US%L_T_to_m_s)
  CS%id_sispeed = register_diag_field('ice_model', 'sispeed', diag%axesT1, Time, &
            'ice speed', 'm/s', missing_value=missing, conversion=US%L_T_to_m_s)

end subroutine SIS_C_dyn_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> find_ice_strength returns the magnitude of force on ice in plastic deformation
subroutine find_ice_strength(mi, ci, ice_strength, G, US, CS, halo_sz) ! ??? may change to do loop
  type(SIS_hor_grid_type),          intent(in)  :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: mi  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: ice_strength !< The ice strength [R Z L2 T-2 ~> Pa m].
  type(SIS_C_dyn_CS),               pointer     :: CS  !< The control structure for this module
  type(unit_scale_type),            intent(in)  :: US  !< A structure with unit conversion factors
  integer,                optional, intent(in)  :: halo_sz !< The halo size to work on
  integer :: i, j, isc, iec, jsc, jec, halo
  halo = 0 ; if (present(halo_sz)) halo = halo_sz
  isc = G%isc-halo ; iec = G%iec+halo ; jsc = G%jsc-halo ; jec = G%jec+halo

  do j=jsc,jec ; do i=isc,iec
    ice_strength(i,j) = CS%p0_rho*mi(i,j)*exp(-CS%c0*max(1.0-ci(i,j),0.0))
  enddo ; enddo

end subroutine find_ice_strength

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_C_dynamics takes a single dynamics timestep with EVP subcycles
subroutine SIS_C_dynamics(ci, mis, mice, ui, vi, uo, vo, fxat, fyat, &
                          sea_lev, fxoc, fyoc, dt_slow, G, US, CS)

  type(SIS_hor_grid_type),           intent(inout) :: G   !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: mis   !< Mass per unit ocean area of sea ice,
                                                            !! snow and melt pond water [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: mice  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: ui    !< Zonal ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: vi    !< Meridional ice velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: uo    !< Zonal ocean velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: vo    !< Meridional ocean velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: fxat  !< Zonal air stress on ice [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: fyat  !< Meridional air stress on ice [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: sea_lev !< The height of the sea level, including
                                                            !! contributions from non-levitating ice converted
                                                            !! to sea-water equivalents, as determined
                                                            !! by the ocean [Z ~> m].

  real, dimension(SZIB_(G),SZJ_(G)), intent(  out) :: fxoc  !< Zonal ice stress on ocean [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)), intent(  out) :: fyoc  !< Meridional ice stress on ocean [R Z L T-2 ~> Pa]
  real,                              intent(in   ) :: dt_slow !< The amount of time over which the ice
                                                            !! dynamics are to be advanced [T ~> s].
  type(unit_scale_type),             intent(in)    :: US    !< A structure with unit conversion factors
  type(SIS_C_dyn_CS),                pointer       :: CS    !< The control structure for this module

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    sh_Dt, &    ! sh_Dt is the horizontal tension (du/dx - dv/dy) including
                ! all metric terms [T-1 ~> s-1].
    sh_Dd       ! sh_Dd is the flow divergence (du/dx + dv/dy) including all
                ! metric terms [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    sh_Ds       ! sh_Ds is the horizontal shearing strain (du/dy + dv/dx)
                ! including all metric terms [T-1 ~> s-1].


  real, dimension(SZI_(G),SZJ_(G)) :: &
    pres_mice, & ! The ice internal pressure per unit column mass [L2 T-2 ~> N m kg-1].
    ci_proj, &  ! The projected ice concentration [nondim].
    zeta, &     ! The ice bulk viscosity [R Z L2 T-1 ~> Pa m s] (i.e., [N s m-1]).
    del_sh, &   ! The magnitude of the shear rates [T-1 ~> s-1].
    diag_val, & ! A temporary diagnostic array.
    del_sh_min_pr, &  ! When multiplied by pres_mice, this gives the minimum
                ! value of del_sh that is used in the calculation of zeta [T-1 ~> s-1].
                ! This is set based on considerations of numerical stability,
                ! and varies with the grid spacing.
    dx2T, dy2T, &   ! dx^2 or dy^2 at T points [L2 ~> m2].
    dx_dyT, dy_dxT, &  ! dx/dy or dy_dx at T points [nondim].
    siu, siv, sispeed  ! diagnostics on T points [L T-1 ~> m s-1].

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    fxic, &   ! Zonal force due to internal stresses [R Z L T-2 ~> Pa].
    fxic_d, & ! Zonal force due to divergence internal stress [R Z L T-2 ~> Pa].
    fxic_t, & ! Zonal force due to tension internal stress [R Z L T-2 ~> Pa].
    fxic_s, & ! Zonal force due to shearing internal stress [R Z L T-2 ~> Pa].
    fxlf, &   ! Zonal landfast ice stress [R Z L T-2 ~> Pa]
    ui_min_trunc, &  ! The range of v-velocities beyond which the velocities
    ui_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells.
    Cor_u, & ! Zonal Coriolis acceleration [L T-2 ~> m s-2].
    PFu, &   ! Zonal hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    diag_val_u, & ! A temporary diagnostic array.
    u_tmp, & ! A temporary copy of the old values of ui [L T-1 ~> m s-1].
    u_IC, &  ! The initial zonal ice velocities [L T-1 ~> m s-1].
    mi_u, &  ! The total ice and snow mass interpolated to u points [R Z ~> kg m-2].
    f2dt_u, &! The squared effective Coriolis parameter at u-points times a
             ! time step [T-1 ~> s-1].
    I1_f2dt2_u  ! 1 / ( 1 + f^2 dt^2) at u-points [nondim].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    fyic, &   ! Meridional force due to internal stresses [R Z L T-2 ~> Pa].
    fyic_d, & ! Meridional force due to divergence internal stress [R Z L T-2 ~> Pa].
    fyic_t, & ! Meridional force due to tension internal stress [R Z L T-2 ~> Pa].
    fyic_s, & ! Meridional force due to shearing internal stress [R Z L T-2 ~> Pa].
    fylf, &   ! Meridional landfast ice stress [R Z L T-2 ~> Pa]
    vi_min_trunc, &  ! The range of v-velocities beyond which the velocities
    vi_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells.
    Cor_v, &  ! Meridional Coriolis acceleration [L T-2 ~> m s-2].
    PFv, &   ! Meridional hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    diag_val_v, & ! A temporary diagnostic array.
    v_IC, &  ! The initial meridional ice velocities [L T-1 ~> m s-1].
    mi_v, &  ! The total ice and snow mass interpolated to v points [R Z ~> kg m-2].
    f2dt_v, &! The squared effective Coriolis parameter at v-points times a
             ! time step [T-1 ~> s-1].
    I1_f2dt2_v  ! 1 / ( 1 + f^2 dt^2) at v-points [nondim].

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    mi_ratio_A_q, & ! A ratio of the masses interpolated to the faces around a
             ! vorticity point that ranges between (4 mi_min/mi_max) and 1,
             ! divided by the sum of the ocean areas around a point [L-2 ~> m-2].
    q, &     ! A potential-vorticity-like field for the ice, the Coriolis parameter
             ! divided by a spatially averaged mass per unit area [T-1 R-1 Z-1 ~> s-1 m2 kg-1].
    dx2B, dy2B, &   ! dx^2 or dy^2 at B points [L2 ~> m2].
    dx_dyB, dy_dxB  ! dx/dy or dy_dx at B points [nondim].
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    azon, bzon, & !  _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vi & ui,
    amer, bmer, & ! respectively to get the barotropic inertial rotation,
    cmer, dmer    ! in units of [T-1 ~> s-1].  azon and amer couple the same pair of
                  ! velocities, but with the influence going in opposite
                  ! directions.

  real :: Cor       ! A Coriolis acceleration [L T-2 ~> m s-2].
  real :: fxic_now  ! Zonal ice internal stress convergence [R Z L T-2 ~> Pa].
  real :: fyic_now  ! Meridional ice internal stress convergence [R Z L T-2 ~> Pa].
  real :: drag_u, drag_v ! Drag rates with the ocean at u & v points [R Z T-1 ~> kg m-2 s-1].
  real :: drag_LFu  ! Drag rates to the land for landfast ice at u points [R Z T-1 ~> kg m-2 s-1].
  real :: drag_LFv  ! Drag rates to the land for landfast ice at v points [R Z T-1 ~> kg m-2 s-1].
  real :: drag_max  ! A maximum drag rate allowed in the ocean [R Z T-1 ~> kg m-2 s-1].
  real :: tot_area  ! The sum of the area of the four neighboring cells [L2 ~> m2].
  real :: dxharm    ! The harmonic mean of the x- and y- grid spacings [L ~> m].
  real :: muq2, mvq2  ! The product of the u- and v-face masses per unit cell
                      ! area surrounding a vorticity point [R2 Z2 ~> kg2 m-4].
  real :: muq, mvq    ! The u- and v-face masses per unit cell area extrapolated
                      ! to a vorticity point on the coast [R Z ~> kg m-2].
  real :: min_rescale ! The smallest of the 4 surrounding values of rescale [nondim].
  real :: I_1pdt_T    ! 1.0 / (1.0 + dt_2Tdamp) [nondim].
  real :: I_1pE2dt_T  ! 1.0 / (1.0 + EC^2 * dt_2Tdamp) [nondim].

  real :: v2_at_u     ! The squared v-velocity interpolated to u points [L2 T-2 ~> m2 s-2].
  real :: u2_at_v     ! The squared u-velocity interpolated to v points [L2 T-2 ~> m2 s-2].
  real :: v2_at_u_min ! The squared v-velocity interpolated to u points [L2 T-2 ~> m2 s-2].
  real :: u2_at_v_min ! The squared u-velocity interpolated to v points [L2 T-2 ~> m2 s-2].
  real :: uio_init    ! Ice-ocean velocity differences [L T-1 ~> m s-1]
  real :: vio_init    ! Ice-ocean velocity differences [L T-1 ~> m s-1]
  real :: m_uio_explicit ! Ice-ocean x-velocity differences times the ice mass [R Z L T-1 ~> kg m-1 s-1]
  real :: m_vio_explicit ! Ice-ocean y-velocity differences times the ice mass [R Z L T-1 ~> kg m-1 s-1]
  real :: uio_pred    ! Ice-ocean x-velocity differences [L T-1 ~> m s-1]
  real :: vio_pred    ! Ice-ocean y-velocity differences [L T-1 ~> m s-1]
  real :: I_cdRhoDt   ! The inverse of the product of the drag coefficient, ocean density and
                      ! timestep [L Z-1 R-1 T-1 ~> m3 kg-1 s-1].
  real :: cdRho       ! The ice density times the drag coefficient and rescaling factors [R Z L-1 ~> kg m-3]
  real :: b_vel0      ! The initial difference between the velocity magnitude
                      ! and the absolute value of the u- or v- component, plus
                      ! the ice thickness divided by the time step and the drag
                      ! coefficient [L T-1 ~> m s-1].
  real :: uio_C   ! A u-velocity difference between the ocean and ice [L T-1 ~> m s-1].
  real :: vio_C   ! A v-velocity difference between the ocean and ice [L T-1 ~> m s-1].

  real :: Tdamp   ! The damping timescale of the stress tensor components
                  ! toward their equilibrium solution due to the elastic terms [T ~> s].
  real :: dt      ! The short timestep associated with the EVP dynamics [T ~> s].
  real :: dt_2Tdamp ! The ratio of the timestep to the elastic damping timescale [nondim].
  real :: dt_cumulative ! The elapsed time within this call to EVP dynamics [T ~> s].
  integer :: EVP_steps ! The number of EVP sub-steps that will actually be taken.
  real :: I_sub_steps  ! The number inverse of the number of EVP time steps per
                  ! slow time step.
  real :: EC2     ! EC^2, where EC is the yield curve axis ratio.
  real :: I_EC2   ! 1/EC^2, where EC is the yield curve axis ratio.
  real :: I_EC    ! 1/EC, where EC is the yield curve axis ratio.
  real :: I_2EC   ! 1/(2*EC), where EC is the yield curve axis ratio.
  real, parameter :: H_subroundoff = 1e-30 ! A negligible ice thickness [m].
  real :: m_neglect  ! A tiny mass per unit area [R Z ~> kg m-2].
  real :: m_neglect2 ! A tiny mass per unit area squared [R2 Z2 ~> kg2 m-4].
  real :: m_neglect4 ! A tiny mass per unit area to the 4th power [R4 Z4 ~> kg4 m-8].
  real :: sum_area   ! The sum of ocean areas around a vorticity point [L2 ~> m2].

  type(time_type) :: &
    time_it_start, &  ! The starting time of the iterative steps.
    time_step_end, &  ! The end time of an iterative step.
    time_end_in       ! The end time for diagnostics when this routine started.
  real :: time_int_in ! The diagnostics' time interval when this routine started.
  logical :: do_hifreq_output  ! If true, output occurs every iterative step.
  logical :: do_trunc_its  ! If true, overly large velocities in the iterations are truncated.
  integer :: halo_sh_Ds  ! The halo size that can be used in calculating sh_Ds.
  integer :: i, j, isc, iec, jsc, jec, n
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "SIS_C_dynamics: Module must be initialized before it is used.")

  if ((isc - G%isdB < 2) .or. (jsc - G%jsdB < 2)) call SIS_error(FATAL, &
         "SIS_C_dynamics is written to require a 2-point halo or 1-point and symmetric memory.")

  halo_sh_Ds = min(isc-G%isd, jsc-G%jsd, 2)

  ! Zero these arrays to accumulate sums.
  fxoc(:,:) = 0.0 ; fyoc(:,:) = 0.0
  fxlf(:,:) = 0.0 ; fylf(:,:) = 0.0
  fxic(:,:) = 0.0 ; fyic(:,:) = 0.0
  Cor_u(:,:) = 0.0 ; Cor_v(:,:) = 0.0
  fxic_d(:,:) = 0.0 ; fyic_d(:,:) = 0.0 ; fxic_t(:,:) = 0.0 ; fyic_t(:,:) = 0.0
  fxic_s(:,:) = 0.0 ; fyic_s(:,:) = 0.0

  if ((CS%evp_sub_steps<=0) .and. (CS%dt_Rheo<=0.0)) return

  if (CS%FirstCall) then
    !   These halo updates are only needed if the str_... arrays have just
    ! been read from a restart file.  Otherwise the halos contain valid data.
    call pass_var(CS%str_d, G%Domain) ; call pass_var(CS%str_t, G%Domain)
    call pass_var(CS%str_s, G%Domain, position=CORNER)
    CS%FirstCall = .false.
  endif

  if (CS%dt_Rheo > 0.0) then
    EVP_steps = max(CEILING(dt_slow/CS%dt_Rheo - 0.0001), 1)
  else
    EVP_steps = CS%evp_sub_steps
  endif
  dt = dt_slow/EVP_steps

  drag_max = CS%Rho_ocean * CS%min_ocn_inertial_h / dt_slow
  I_cdRhoDt = 1.0 / (CS%cdw * US%L_to_Z*CS%Rho_ocean * dt)
  do_trunc_its = (CS%CFL_check_its .and. (CS%CFL_trunc > 0.0) .and. (dt_slow > 0.0))

  EC2 = CS%EC**2
  I_EC = 0.0 ; if (CS%EC > 0.0) I_EC = 1.0 / CS%EC
  I_2EC = 0.0 ; if (CS%EC > 0.0) I_2EC = 0.5 / CS%EC
  I_EC2 = 0.0 ; if (EC2 > 0.0) I_EC2 = 1.0 / EC2

  do_hifreq_output = .false.
  if ((CS%id_ui_hifreq > 0) .or. (CS%id_vi_hifreq > 0) .or. &
      (CS%id_str_d_hifreq > 0) .or. (CS%id_str_t_hifreq > 0) .or. &
      (CS%id_str_s_hifreq > 0) .or. (CS%id_sh_d_hifreq > 0) .or. &
      (CS%id_sh_t_hifreq > 0) .or. (CS%id_sh_s_hifreq > 0) .or. &
      (CS%id_ci_hifreq > 0) .or. (CS%id_stren_hifreq > 0)) then
    do_hifreq_output = query_SIS_averaging_enabled(CS%diag, time_int_in, time_end_in)
    if (do_hifreq_output) &
      time_it_start = time_end_in - real_to_time(US%T_to_s*dt_slow)
  endif

  Tdamp = CS%Tdamp
  if (CS%Tdamp == 0.0) then
    ! Hunke (2001) chooses a specified multiple (0.36) of dt_slow for Tdamp, and shows that
    ! stability requires Tdamp > 2*dt.  Here 0.2 is used instead for greater stability.
    Tdamp = max(0.2*dt_slow, 3.0*dt)
  elseif (CS%Tdamp < 0.0) then
    Tdamp = max(-CS%Tdamp*dt_slow, 3.0*dt)
  endif
  dt_2Tdamp = dt / (2.0 * Tdamp)

  ui_min_trunc(:,:) = 0.0 ; ui_max_trunc(:,:) = 0.0
  vi_min_trunc(:,:) = 0.0 ; vi_max_trunc(:,:) = 0.0

  m_neglect = H_subroundoff*US%m_to_Z*CS%Rho_ice
  m_neglect2 = m_neglect**2 ; m_neglect4 = m_neglect**4
!$OMP parallel default(none) shared(isc,iec,jsc,jec,G,US,CS,dt_slow,ui_min_trunc,u_IC,ui,   &
!$OMP                               ui_max_trunc,vi_min_trunc,vi_max_trunc,v_IC,vi,mice, &
!$OMP                               mis,ci,dt,Tdamp,I_2EC,ci_proj,pres_mice,       &
!$OMP                               dx2B,dy2B,dx_dyB,dy_dxB,dx2T,dy2T,dx_dyT,dy_dxT,     &
!$OMP                               mi_ratio_A_q,m_neglect4,m_neglect2,mi_u,mi_v,q,      &
!$OMP                               m_neglect,azon,bzon,czon,dzon,f2dt_u,I1_f2dt2_u,PFu, &
!$OMP                               sea_lev,amer,bmer,cmer,dmer,f2dt_v,I1_f2dt2_v,PFv,   &
!$OMP                               del_sh_min_pr                         )              &
!$OMP                       private(dxharm,sum_area,muq2,mvq2,muq,mvq,tot_area)
  if ((CS%CFL_trunc > 0.0) .and. (dt_slow > 0.0)) then
!$OMP do
    do j=jsc,jec
      do I=isc-1,iec ; if (G%dy_Cu(I,j) > 0.0) then
        ui_min_trunc(I,j) = (-CS%CFL_trunc) * G%areaT(i+1,j) / (dt_slow*G%dy_Cu(I,j))
        ui_max_trunc(I,j) = CS%CFL_trunc * G%areaT(i,j) / (dt_slow*G%dy_Cu(I,j))
      endif ; enddo
      do I=isc-1,iec ; u_IC(I,j) = ui(I,j) ; enddo
    enddo
!$OMP end do nowait
!$OMP do
    do J=jsc-1,jec
      do i=isc,iec ; if (G%dx_Cv(i,J) > 0.0) then
        vi_min_trunc(i,J) = (-CS%CFL_trunc) * G%areaT(i,j+1) / (dt_slow*G%dx_Cv(i,J))
        vi_max_trunc(i,J) = CS%CFL_trunc * G%areaT(i,j) / (dt_slow*G%dx_Cv(i,J))
      endif ; enddo
      do i=isc,iec ; v_IC(i,J) = vi(i,j) ; enddo
    enddo
!$OMP end do nowait
  endif
!$OMP do
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    ci_proj(i,j) = ci(i,j)

    ! Precompute pres_mice and the minimum value of del_sh for stability.
    pres_mice(i,j) = CS%p0_rho*exp(-CS%c0*max(1.0-ci(i,j),0.0))

    dxharm = 2.0*G%dxT(i,j)*G%dyT(i,j) / (G%dxT(i,j) + G%dyT(i,j))
    !   Setting a minimum value of del_sh is sufficient to guarantee numerical
    ! stability of the overall time-stepping for the velocities and stresses.
    ! Setting a minimum value of the shear magnitudes is equivalent to setting
    ! a maximum value of the effective lateral viscosities.
    ! I think that this is stable when CS%del_sh_min_scale >= 1.  -RWH
    if (dxharm > 0.) then
      del_sh_min_pr(i,j) = (2.0*CS%del_sh_min_scale * dt**2) / (Tdamp * dxharm**2)
    else
      del_sh_min_pr(i,j) = 0.
    endif
  enddo ; enddo

  ! Ensure that the input stresses are not larger than could be justified by
  ! the ice pressure now, as the ice might have melted or been advected away
  ! during the thermodynamic and transport phases.
  call limit_stresses(pres_mice, mice(:,:), CS%str_d, CS%str_t, CS%str_s, G, US, CS)

  ! Zero out ice velocities with no mass.
!$OMP do
  do j=jsc,jec ; do I=isc-1,iec
    if (G%mask2dCu(I,j) * (mis(i,j)+mis(i+1,j)) == 0.0) ui(I,j) = 0.0
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-1,jec ; do i=isc,iec
    if (G%mask2dCv(i,J) * (mis(i,j)+mis(i,j+1)) == 0.0) vi(I,j) = 0.0
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-1,jec ; do I=isc-1,iec
    dx2B(I,J) = G%dxBu(I,J)*G%dxBu(I,J) ; dy2B(I,J) = G%dyBu(I,J)*G%dyBu(I,J)
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-2,jec+1 ; do I=isc-2,iec+1
    dx_dyB(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; dy_dxB(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    dx2T(i,j) = G%dxT(i,j)*G%dxT(i,j) ; dy2T(i,j) = G%dyT(i,j)*G%dyT(i,j)
    dx_dyT(i,j) = G%dxT(i,j)*G%IdyT(i,j) ; dy_dxT(i,j) = G%dyT(i,j)*G%IdxT(i,j)
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-1,jec ; do I=isc-1,iec
    if (CS%weak_coast_stress) then
      sum_area = (G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i,j+1) + G%areaT(i+1,j))
    else
      sum_area = (G%mask2dT(i,j)*G%areaT(i,j) + G%mask2dT(i+1,j+1)*G%areaT(i+1,j+1)) + &
                 (G%mask2dT(i,j+1)*G%areaT(i,j+1) + G%mask2dT(i+1,j)*G%areaT(i+1,j))
    endif
    if (sum_area <= 0.0) then
      ! This is a land point.
      mi_ratio_A_q(I,J) = 0.0
    elseif (G%mask2dBu(I,J)>0.0) then
      ! This is an interior ocean point.
      !   Determine an appropriately averaged mass on q-points. The following
      ! expression for mi_q is mi when the masses are all equal, and goes to 4
      ! times the smallest mass averaged onto the 4 adjacent velocity points.  It
      ! comes from taking the harmonic means of the harmonic means of the
      ! arithmetic mean masses at the velocity points.  mi_ratio goes from 4 times
      ! the ratio of the smallest mass at velocity points over the largest mass
      ! at velocity points up to 1.
      muq2 = 0.25 * (mis(i,j) + mis(i+1,j)) * (mis(i,j+1) + mis(i+1,j+1))
      mvq2 = 0.25 * (mis(i,j) + mis(i,j+1)) * (mis(i+1,j) + mis(i+1,j+1))
      mi_ratio_A_q(I,J) = 32.0 * muq2 * mvq2 / ((m_neglect4 + (muq2 + mvq2) * &
                 ((mis(i,j) + mis(i+1,j+1)) + (mis(i,j+1) + mis(i+1,j)))**2) * sum_area)
    elseif ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
            (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) > 1.5) then
      !   This is a corner point, and there are 1 unmasked u-point and 1 v-point.
      ! The ratio below goes from 4 times the ratio of the smaller of the two
      ! masses at velocity points over the larger up to 1.
      muq = 0.5 * (G%mask2dCu(I,j) * (mis(i,j) + mis(i+1,j)) + &
                   G%mask2dCu(I,j+1) * (mis(i,j+1) + mis(i+1,j+1)) )
      mvq = 0.5 * (G%mask2dCv(i,J) * (mis(i,j) + mis(i,j+1)) + &
                   G%mask2dCv(i+1,J) * (mis(i+1,j) + mis(i+1,j+1)) )
      mi_ratio_A_q(I,J) = 4.0 * muq * mvq / ((m_neglect2 + (muq + mvq)**2) * sum_area)
    else
      ! This is a straight coastline or all neighboring velocity points are
      ! masked out.  In any case, with just 1 point, the ratio is always 1.
      mi_ratio_A_q(I,J) = 1.0 / sum_area
    endif
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do j=jsc-1,jec+1 ; do I=isc-1,iec
    mi_u(I,j) = 0.5*(mis(i+1,j) + mis(i,j))
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-1,jec ; do i=isc-1,iec+1
    mi_v(i,J) = 0.5*(mis(i,j+1) + mis(i,j))
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-1,jec ; do I=isc-1,iec
    tot_area = (G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i+1,j) + G%areaT(i,j+1))
    q(I,J) = G%CoriolisBu(I,J) * tot_area / &
         (((G%areaT(i,j) * mis(i,j) + G%areaT(i+1,j+1) * mis(i+1,j+1)) + &
           (G%areaT(i+1,j) * mis(i+1,j) + G%areaT(i,j+1) * mis(i,j+1))) + tot_area * m_neglect)
  enddo ; enddo
!$OMP do
  do j=jsc,jec ; do I=isc-1,iec
    ! Calculate terms related to the Coriolis force on the zonal velocity.
    azon(I,j) = 0.25 * mi_v(i+1,J) * q(I,J)
    bzon(I,j) = 0.25 * mi_v(i,J) * q(I,J)
    czon(I,j) = 0.25 * mi_v(i,J-1) * q(I,J-1)
    dzon(I,j) = 0.25 * mi_v(i+1,J-1) * q(I,J-1)

    f2dt_u(I,j) = dt * 4.0 * ((azon(I,j)**2 + czon(I,j)**2) + &
                              (bzon(I,j)**2 + dzon(I,j)**2))
    I1_f2dt2_u(I,j) = 1.0 / ( 1.0 + dt * f2dt_u(I,j) )

    ! Calculate the zonal acceleration due to the sea level slope.
    PFu(I,j) = -G%g_Earth*(sea_lev(i+1,j)-sea_lev(i,j)) * G%IdxCu(I,j)
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-1,jec ; do i=isc,iec
    ! Calculate terms related to the Coriolis force on the meridional velocity.
    amer(I-1,j) = 0.25 * mi_u(I-1,j) * q(I-1,J)
    bmer(I,j) = 0.25 * mi_u(I,j) * q(I,J)
    cmer(I,j+1) = 0.25 * mi_u(I,j+1) * q(I,J)
    dmer(I-1,j+1) = 0.25 * mi_u(I-1,j+1) * q(I-1,J)

    f2dt_v(i,J) = dt * 4.0 * ((amer(I-1,j)**2 + cmer(I,j+1)**2) + &
                              (bmer(I,j)**2 + dmer(I-1,j+1)**2))
    I1_f2dt2_v(i,J) = 1.0 / ( 1.0 + dt * f2dt_v(i,J) )

    ! Calculate the meridional acceleration due to the sea level slope.
    PFv(i,J) = -G%g_Earth*(sea_lev(i,j+1)-sea_lev(i,j)) * G%IdyCv(i,J)
  enddo ; enddo
!$OMP end parallel

  if (CS%debug .or. CS%debug_redundant) then
    call uvchksum("PF[uv] in SIS_C_dynamics", PFu, PFv, G, scale=US%L_T_to_m_s*US%s_to_T)
    call uvchksum("f[xy]at in SIS_C_dynamics", fxat, fyat, G, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    call uvchksum("[uv]i pre-steps SIS_C_dynamics", ui, vi, G, scale=US%L_T_to_m_s)
    call uvchksum("[uv]o in SIS_C_dynamics", uo, vo, G, scale=US%L_T_to_m_s)
  endif

  dt_cumulative = 0.0

  ! Do the iterative time steps.
  do n=1,EVP_steps

    dt_cumulative = dt_cumulative + dt
    ! If there is a 2-point wide halo and symmetric memory, this is the only
    ! halo update that is needed per iteration.  With a 1-point wide halo and
    ! symmetric memory, an update is also needed for sh_Ds.
    call pass_vector(ui, vi, G%Domain, stagger=CGRID_NE)

    !    Calculate the strain tensor for viscosities and forcing elastic eqn.
    !  The following are the forms of the horizontal tension and horizontal
    !  shearing strain advocated by Smagorinsky (1993) and discussed in
    !  Griffies and Hallberg (MWR, 2000).  Similar forms are used in the sea
    !  ice model of Bouillon et al. (Ocean Modelling, 2009).

    !   The calculation of sh_Ds has the widest halo. The logic below avoids
    ! a halo update when possible.
    !   With a halo of >= 2 this is:  do J=jsc-2,jec+1 ; do I=isc-2,iec+1
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,halo_sh_Ds,sh_Ds,G, &
!$OMP                                  dx_dyB,dy_dxB,ui,vi)
    do J=jsc-halo_sh_Ds,jec+halo_sh_Ds-1 ; do I=isc-halo_sh_Ds,iec+halo_sh_Ds-1
      ! This uses a no-slip boundary condition.
      sh_Ds(I,J) = (2.0-G%mask2dBu(I,J)) * &
          (dx_dyB(I,J)*(ui(I,j+1)*G%IdxCu(I,j+1) - ui(I,j)*G%IdxCu(I,j)) + &
           dy_dxB(I,J)*(vi(i+1,J)*G%IdyCv(i+1,J) - vi(i,J)*G%IdyCv(i,J)))
    enddo ; enddo
    if (halo_sh_Ds < 2) call pass_var(sh_Ds, G%Domain, position=CORNER)
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,sh_Dt,sh_Dd,dy_dxT,dx_dyT,G,ui,vi)
    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      sh_Dt(i,j) = (dy_dxT(i,j)*(G%IdyCu(I,j) * ui(I,j) - &
                                 G%IdyCu(I-1,j)*ui(I-1,j)) - &
                    dx_dyT(i,j)*(G%IdxCv(i,J) * vi(i,J) - &
                                 G%IdxCv(i,J-1)*vi(i,J-1)))
      sh_Dd(i,j) = (G%IareaT(i,j)*(G%dyCu(I,j) * ui(I,j) - &
                                   G%dyCu(I-1,j)*ui(I-1,j)) + &
                    G%IareaT(i,j)*(G%dxCv(i,J) * vi(i,J) - &
                                   G%dxCv(i,J-1)*vi(i,J-1)))
    enddo ; enddo

   if (CS%project_ci) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,ci_proj,ci,dt_cumulative, &
!$OMP                                  sh_Dd,pres_mice,CS)
     do j=jsc-1,jec+1 ; do i=isc-1,iec+1
       ! Estimate future ice concentrations from the approximate expression
       !   d ci / dt = - ci * sh_Dt
       ! The choice to base this on the final velocity, the initial concentration
       ! and the elapsed time is because it is that final velocity that will drive
       ! ice convergence.
       ci_proj(i,j) = ci(i,j) * exp(-dt_cumulative*sh_Dd(i,j))
       ! Recompute pres_mice.
       pres_mice(i,j) = CS%p0_rho*exp(-CS%c0*max(1.0-ci_proj(i,j),0.0))
     enddo ; enddo
   endif

   ! calculate viscosities - how often should we do this ?
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,del_sh,zeta,sh_Dd,sh_Dt, &
!$OMP                                  I_EC2,sh_Ds,pres_mice,mice,del_sh_min_pr)
    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      ! Averaging the squared shearing strain is larger than squaring
      ! the averaged strain.  I don't know what is better. -RWH
      del_sh(i,j) = sqrt(sh_Dd(i,j)**2 + I_EC2 * (sh_Dt(i,j)**2 + &
                   (0.25 * ((sh_Ds(I-1,J-1) + sh_Ds(I,J)) + &
                            (sh_Ds(I-1,J) + sh_Ds(I,J-1))))**2 ) ) ! H&D eqn 9

      if (max(del_sh(i,j), del_sh_min_pr(i,j)*pres_mice(i,j)) /= 0.) then
        zeta(i,j) = 0.5*pres_mice(i,j)*mice(i,j) / &
           max(del_sh(i,j), del_sh_min_pr(i,j)*pres_mice(i,j))
      else
        zeta(i,j) = 0.
      endif
    enddo ; enddo

    ! Step the stress component equations semi-implicitly.
    I_1pdt_T = 1.0 / (1.0 + dt_2Tdamp)
    I_1pE2dt_T = 1.0 / (1.0 + EC2*dt_2Tdamp)
    if (CS%weak_low_shear) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,CS,I_1pdt_T,dt_2Tdamp,zeta, &
!$OMP                                  sh_Dd,del_sh,I_EC2,sh_Dt)
      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        ! This expression uses that Pres=2*del_sh*zeta with an elliptic yield curve.
        CS%str_d(i,j) = I_1pdt_T * ( CS%str_d(i,j) + dt_2Tdamp * &
                    ( zeta(i,j) * (sh_Dd(i,j) - del_sh(i,j)) ) )
        CS%str_t(i,j) = I_1pdt_T * ( CS%str_t(i,j) + (I_EC2 * dt_2Tdamp) * &
                    ( zeta(i,j) * sh_Dt(i,j) ) )
      enddo ; enddo
    else
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,CS,I_1pdt_T,dt_2Tdamp,zeta, &
!$OMP                                  sh_Dd,I_EC2,sh_Dt,pres_mice,mice)
      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        ! This expression uses that Pres=2*del_sh*zeta with an elliptic yield curve.
        CS%str_d(i,j) = I_1pdt_T * ( CS%str_d(i,j) + dt_2Tdamp * &
                    ( zeta(i,j) * sh_Dd(i,j) - 0.5*pres_mice(i,j)*mice(i,j) ) )
        CS%str_t(i,j) = I_1pdt_T * ( CS%str_t(i,j) + (I_EC2 * dt_2Tdamp) * &
                    ( zeta(i,j) * sh_Dt(i,j) ) )
      enddo ; enddo
    endif
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,CS,I_1pdt_T,I_EC2,dt_2Tdamp, &
!$OMP                                  G,zeta,mi_ratio_A_q,sh_Ds)
    do J=jsc-1,jec ; do I=isc-1,iec
      ! zeta is already set to 0 over land.
      CS%str_s(I,J) = I_1pdt_T * ( CS%str_s(I,J) + (I_EC2 * dt_2Tdamp) * &
                  ( ((G%areaT(i,j)*zeta(i,j) + G%areaT(i+1,j+1)*zeta(i+1,j+1)) + &
                     (G%areaT(i+1,j)*zeta(i+1,j) + G%areaT(i,j+1)*zeta(i,j+1))) * &
                   mi_ratio_A_q(I,J) * sh_Ds(I,J) ) )
    enddo ; enddo

    if (CS%str_underflow > 0.0) then
      !$OMP parallel do default(shared)
      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        if (abs(CS%str_d(i,j)) < CS%str_underflow) CS%str_d(i,j) = 0.0
        if (abs(CS%str_t(i,j)) < CS%str_underflow) CS%str_t(i,j) = 0.0
      enddo ; enddo
      !$OMP parallel do default(shared)
      do J=jsc-1,jec ; do I=isc-1,iec
        if (abs(CS%str_s(I,J)) < CS%str_underflow) CS%str_s(I,J) = 0.0
      enddo ; enddo
    endif

    cdRho = CS%cdw * US%L_to_Z*CS%Rho_ocean
    ! Save the current values of u for later use in updating v.
    do I=isc-1,iec
      u_tmp(I,jsc-1) = ui(I,jsc-1) ; u_tmp(I,jec+1) = ui(I,jec+1) ;
    enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,u_tmp,ui,vi,azon,bzon,czon,dzon, &
!$OMP                                  G,CS,dy2T,dx2B,vo,uo,Cor_u,f2dt_u,I1_f2dt2_u,    &
!$OMP                                  mi_u,dt,PFu,fxat,I_cdRhoDt,cdRho,m_neglect,fxoc, &
!$OMP                                  fxlf,fxic,fxic_d,fxic_t,fxic_s,do_trunc_its,drag_max) &
!$OMP                          private(Cor,fxic_now,v2_at_u,v2_at_u_min,uio_init,drag_u,drag_LFu,b_vel0, &
!$OMP                                  m_uio_explicit,uio_pred,uio_C)
    do j=jsc,jec ; do I=isc-1,iec
      ! Save the current values of u for later use in updating v.
      u_tmp(I,j) = ui(I,j)

      Cor = ((azon(I,j) * vi(i+1,J) + czon(I,j) * vi(i,J-1)) + &
             (bzon(I,j) * vi(i,J) + dzon(I,j) * vi(i+1,J-1))) ! - Cor_ref_u(I,j)
      !  Evaluate 1/m x.Div(m strain).  This expressions include all metric terms
      !  for an orthogonal grid.  The str_d term integrates out to no curl, while
      !  str_s & str_t terms impose no divergence and do not act on solid body rotation.
      fxic_now = G%IdxCu(I,j) * (CS%str_d(i+1,j) - CS%str_d(i,j)) + &
            (G%IdyCu(I,j)*(dy2T(i+1,j)*CS%str_t(i+1,j) - &
                           dy2T(i,j)  *CS%str_t(i,j)) + &
             G%IdxCu(I,j)*(dx2B(I,J)  *CS%str_s(I,J) - &
                           dx2B(I,J-1)*CS%str_s(I,J-1)) ) * G%IareaCu(I,j)
      v2_at_u =  CS%drag_bg_vel2 + 0.25 * &
                     (((vi(i,J)-vo(i,J))**2 + (vi(i+1,J-1)-vo(i+1,J-1))**2) + &
                      ((vi(i+1,J)-vo(i+1,J))**2 + (vi(i,J-1)-vo(i,J-1))**2))
      if (CS%lemieux_landfast .or. CS%itd_landfast) &
               v2_at_u_min = min(abs(vi(I,j)), abs(vi(i+1,J-1)), &
                                 abs(vi(i+1,J)), abs(vi(i,J-1)))**2

      uio_init = (ui(I,j)-uo(I,j))

      ! Determine the Coriolis acceleration and sum for averages...
      Cor_u(I,j) = Cor_u(I,j) + (Cor - f2dt_u(I,j) * ui(I,j)) * I1_f2dt2_u(I,j)

      if (CS%project_drag_vel) then
      ! Project the new u-velocity using a quasi-analytic implicit treatment for
      ! drag, but explicit treatments for everything else, to estimate the drag
      ! coefficient, then take the larger of the two estimates of
      ! the ice-ocean drag.
        drag_u = 0.0
        if (G%mask2dCu(I,j) > 0.0) then
          m_uio_explicit = uio_init*mi_u(I,j) + dt * &
               ((Cor + PFu(I,j))*mi_u(I,j) + (fxic_now + fxat(I,j)))
          b_vel0 = mi_u(I,j) * I_cdRhoDt + &
                   ( sqrt(uio_init**2 + v2_at_u) - abs(uio_init) )
          if (b_vel0**2 > 1e8*I_cdRhoDt*abs(m_uio_explicit)) then
            uio_pred = m_uio_explicit * I_cdRhoDt / b_vel0
          else
            uio_pred = 0.5 * (sqrt(b_vel0**2 + 4.0*I_cdRhoDt*abs(m_uio_explicit)) - b_vel0)
          endif
          drag_u = cdRho * sqrt(max(uio_init**2, uio_pred**2) + v2_at_u )
        endif
      else
        drag_u = cdRho * sqrt(uio_init**2 + v2_at_u )
      endif
      if (drag_max>0.) drag_u = min( drag_u, drag_max )
      drag_LFu = 0.0
      if (CS%lemieux_landfast .or. CS%itd_landfast) then
        drag_LFu = CS%Tb_u(I,j) / (sqrt(ui(I,j)**2 + v2_at_u_min ) + CS%lemieux_u0)
      endif

      !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
      ! update of the other terms and an implicit bottom drag calculation.
      uio_C =  G%mask2dCu(I,j) * ( mi_u(I,j) * &
               ((ui(I,j) + dt * Cor) * I1_f2dt2_u(I,j) - uo(I,j)) + &
                dt * ((mi_u(I,j) * PFu(I,j) + (fxic_now + fxat(I,j))) - drag_LFu*uo(I,j)) ) / &
               (mi_u(I,j) + m_neglect + dt * (drag_u + drag_LFu))

      ui(I,j) = (uio_C + uo(I,j)) * G%mask2dCu(I,j)

      ! Note that fxoc is the stress felt by the ocean.
      fxoc(I,j) = fxoc(I,j) + drag_u*uio_C

      ! Here fxlf is the stress felt by the landfast ice.
      fxlf(I,j) = fxlf(I,j) - drag_LFu*ui(I,j)

      ! sum accelerations to take averages.
      fxic(I,j) = fxic(I,j) + fxic_now

      if (CS%id_fix_d>0) fxic_d(I,j) = fxic_d(I,j) + G%mask2dCu(I,j) * &
                 G%IdxCu(I,j) * (CS%str_d(i+1,j) - CS%str_d(i,j))
      if (CS%id_fix_t>0) fxic_t(I,j) = fxic_t(I,j) + G%mask2dCu(I,j) * &
                  G%IdyCu(I,j)*(dy2T(i+1,j)* CS%str_t(i+1,j) - &
                                dy2T(i,j)  * CS%str_t(i,j) ) * G%IareaCu(I,j)
      if (CS%id_fix_s>0) fxic_s(I,j) = fxic_s(I,j) + G%mask2dCu(I,j) * &
                  G%IdxCu(I,j)*(dx2B(I,J)  *CS%str_s(I,J) - &
                                dx2B(I,J-1)*CS%str_s(I,J-1)) * G%IareaCu(I,j)

    enddo ; enddo

!$OMP parallel do default(none) shared(isc,iec,jsc,jec,amer,bmer,cmer,dmer,u_tmp,G,CS, &
!$OMP                                  dx2T,dy2B,uo,vo,vi,Cor_v,f2dt_v,I1_f2dt2_v,mi_v, &
!$OMP                                  dt,PFv,fyat,I_cdRhoDt,cdRho,m_neglect,fyoc,fyic, &
!$OMP                                  fylf,fyic_d,fyic_t,fyic_s,do_trunc_its,vi_min_trunc,  &
!$OMP                                  vi_max_trunc,drag_max) &
!$OMP                          private(Cor,fyic_now,u2_at_v,vio_init,drag_v,drag_LFv,u2_at_v_min, &
!$OMP                                  m_vio_explicit,b_vel0,vio_pred,vio_C)
    do J=jsc-1,jec ; do i=isc,iec
      Cor = -1.0*((amer(I-1,j) * u_tmp(I-1,j) + cmer(I,j+1) * u_tmp(I,j+1)) + &
                  (bmer(I,j) * u_tmp(I,j) + dmer(I-1,j+1) * u_tmp(I-1,j+1)))
      !  Evaluate 1/m y.Div(m strain).  This expressions include all metric terms
      !  for an orthogonal grid.  The str_d term integrates out to no curl, while
      !  str_s & str_t terms impose no divergence and do not act on solid body rotation.
      fyic_now = G%IdyCv(i,J) * (CS%str_d(i,j+1)-CS%str_d(i,j)) + &
            (-G%IdxCv(i,J)*(dx2T(i,j+1)*CS%str_t(i,j+1) - &
                            dx2T(i,j)  *CS%str_t(i,j)) + &
              G%IdyCv(i,J)*(dy2B(I,J)  *CS%str_s(I,J) - &
                            dy2B(I-1,J)*CS%str_s(I-1,J)) )*G%IareaCv(i,J)
      u2_at_v = CS%drag_bg_vel2 + 0.25 * &
                (((u_tmp(I,j)-uo(I,j))**2 + (u_tmp(I-1,j+1)-uo(I-1,j+1))**2) + &
                 ((u_tmp(I,j+1)-uo(I,j+1))**2 + (u_tmp(I-1,j)-uo(I-1,j))**2))
      if (CS%lemieux_landfast .or. CS%itd_landfast) &
                u2_at_v_min = min(abs(u_tmp(i,J)), abs(u_tmp(I-1,j+1)), &
                                  abs(u_tmp(I,j+1)), abs(u_tmp(I-1,j)))**2

      vio_init = (vi(i,J)-vo(i,J))

      ! Determine the Coriolis acceleration and sum for averages...
      Cor_v(I,J) = Cor_v(I,J) + (Cor - f2dt_v(i,J) * vi(i,J)) * I1_f2dt2_v(i,J)

      if (CS%project_drag_vel) then
      ! Project the new v-velocity using a quasi-analytic implicit treatment for
      ! drag, but explicit treatments for everything else, to estimate the drag
      ! coefficient, then take the larger of the two estimates of
      ! the ice-ocean drag.

        drag_v = 0.0
        if (G%mask2dCv(i,J) > 0.0) then
          m_vio_explicit = vio_init*mi_v(i,J) + dt * &
               ((Cor + PFv(i,J))*mi_v(i,J) + (fyic_now + fyat(i,J)))
          b_vel0 = mi_v(i,J) * I_cdRhoDt + (sqrt(vio_init**2 + u2_at_v) - abs(vio_init))
          if (b_vel0**2 > 1e8*I_cdRhoDt*abs(m_vio_explicit)) then
            vio_pred = m_vio_explicit * I_cdRhoDt / b_vel0
          else
            vio_pred = 0.5 * (sqrt(b_vel0**2 + 4.0*I_cdRhoDt*abs(m_vio_explicit)) - b_vel0)
          endif
          drag_v = cdRho * sqrt(max(vio_init**2, vio_pred**2) + u2_at_v )
        endif
      else
        drag_v = cdRho * sqrt(vio_init**2 + u2_at_v )
      endif
      if (drag_max>0.) drag_v = min( drag_v, drag_max )
      drag_LFv = 0.0
      if (CS%lemieux_landfast .or. CS%itd_landfast) then
        drag_LFv = CS%Tb_v(i,J) / (sqrt(vi(i,J)**2 + u2_at_v_min ) + CS%lemieux_u0)
      endif

      !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
      ! update of the other terms and an implicit bottom drag calculation.
      vio_C =  G%mask2dCv(i,J) * ( mi_v(i,J) * &
               ((vi(i,J) + dt * Cor) * I1_f2dt2_v(i,J) - vo(i,J)) + &
                dt * ((mi_v(i,J) * PFv(i,J) + (fyic_now + fyat(i,J))) - drag_LFv*vo(i,J)) ) / &
               (mi_v(i,J) + m_neglect + dt * (drag_v + drag_LFv))

      vi(i,J) = (vio_C + vo(i,J)) * G%mask2dCv(i,J)

      ! Note that fyoc is the stress felt by the ocean.
      fyoc(i,J) = fyoc(i,J) + drag_v*vio_C

      ! Here fylf is the stress felt by the landfast ice.
      fylf(I,j) = fylf(I,j) - drag_LFv*vi(I,j)

      ! sum accelerations to take averages.
      fyic(i,J) = fyic(i,J) + fyic_now

      if (CS%id_fiy_d>0) fyic_d(i,J) = fyic_d(i,J) + G%mask2dCv(i,J) * &
                 G%IdyCv(i,J) * (CS%str_d(i,j+1)-CS%str_d(i,j))
      if (CS%id_fiy_t>0) fyic_t(i,J) = fyic_t(i,J) + G%mask2dCv(i,J) * &
                 (G%IdxCv(i,J)*(dx2T(i,j+1)*(-CS%str_t(i,j+1)) - &
                                dx2T(i,j)  *(-CS%str_t(i,j))) ) * G%IareaCv(i,J)
      if (CS%id_fiy_s>0) fyic_s(i,J) = fyic_s(i,J) + G%mask2dCv(i,J) * &
                 (G%IdyCv(i,J)*(dy2B(I,J)  *CS%str_s(I,J) - &
                                dy2B(I-1,J)*CS%str_s(I-1,J)) ) * G%IareaCv(i,J)

    enddo ; enddo

    ! Apply appropriate limits on the magnitude of the velocies, both to handle
    ! underflow and to keep failing runs going so that they can be diagnosed.
    if (do_trunc_its .or. (CS%vel_underflow > 0.0)) then
      !$OMP parallel do default(shared)
      do j=jsc,jec ; do I=isc-1,iec
        if (abs(ui(I,j)) < CS%vel_underflow) ui(I,j) = 0.0
        if (do_trunc_its) then
          if (ui(I,j) < ui_min_trunc(I,j)) then
            ui(I,j) = ui_min_trunc(I,j)
          elseif (ui(I,j) > ui_max_trunc(I,j)) then
            ui(I,j) = ui_max_trunc(I,j)
          endif
        endif
      enddo ; enddo

      !$OMP parallel do default(shared)
      do J=jsc-1,jec ; do i=isc,iec
        if (abs(vi(i,J)) < CS%vel_underflow) vi(i,J) = 0.0
        if (do_trunc_its) then
          if (vi(i,J) < vi_min_trunc(i,J)) then
            vi(i,J) = vi_min_trunc(i,J)
          elseif (vi(i,J) > vi_max_trunc(i,J)) then
            vi(i,J) = vi_max_trunc(i,J)
          endif
        endif
      enddo ; enddo
    endif

    if (do_hifreq_output) then
      time_step_end = time_it_start + real_to_time(n*US%T_to_s*dt)
      call enable_SIS_averaging(US%T_to_s*dt, time_step_end, CS%diag)
      if (CS%id_ui_hifreq > 0) call post_SIS_data(CS%id_ui_hifreq, ui, CS%diag)
      if (CS%id_vi_hifreq > 0) call post_SIS_data(CS%id_vi_hifreq, vi, CS%diag)
      if (CS%id_str_d_hifreq > 0) call post_SIS_data(CS%id_str_d_hifreq, CS%str_d, CS%diag)
      if (CS%id_str_t_hifreq > 0) call post_SIS_data(CS%id_str_t_hifreq, CS%str_t, CS%diag)
      if (CS%id_str_s_hifreq > 0) call post_SIS_data(CS%id_str_s_hifreq, CS%str_s, CS%diag)
      if (CS%id_sh_d_hifreq > 0) call post_SIS_data(CS%id_sh_d_hifreq, sh_Dd, CS%diag)
      if (CS%id_sh_t_hifreq > 0) call post_SIS_data(CS%id_sh_t_hifreq, sh_Dt, CS%diag)
      if (CS%id_sh_s_hifreq > 0) call post_SIS_data(CS%id_sh_s_hifreq, sh_Ds, CS%diag)
      if (CS%id_sigi_hifreq>0) then
        call find_sigI(mice, ci_proj, CS%str_d, diag_val, G, US, CS)
        call post_SIS_data(CS%id_sigi_hifreq, diag_val, CS%diag)
      endif
      if (CS%id_sigii_hifreq>0) then
        call find_sigII(mice, ci_proj, CS%str_t, CS%str_s, diag_val, G, US, CS)
        call post_SIS_data(CS%id_sigii_hifreq, diag_val, CS%diag)
      endif
      if (CS%id_ci_hifreq>0) call post_SIS_data(CS%id_ci_hifreq, ci_proj, CS%diag)
      if (CS%id_stren_hifreq>0) then
        do j=jsc,jec ; do i=isc,iec
          diag_val(i,j) = pres_mice(i,j)*mice(i,j)
        enddo ; enddo
        call post_SIS_data(CS%id_stren_hifreq, diag_val, CS%diag)
      endif
    endif

    if (CS%debug_EVP .and. CS%debug) then
      call hchksum(CS%str_d, "str_d in SIS_C_dynamics", G%HI, haloshift=1, scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
      call hchksum(CS%str_t, "str_t in SIS_C_dynamics", G%HI, haloshift=1, scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
      call Bchksum(CS%str_s, "str_s in SIS_C_dynamics", G%HI, &
                   haloshift=0, symmetric=.true., scale=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
    endif
    if (CS%debug_EVP .and. (CS%debug .or. CS%debug_redundant)) then
      call uvchksum("f[xy]ic in SIS_C_dynamics", fxic, fyic, G, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call uvchksum("f[xy]oc in SIS_C_dynamics", fxoc, fyoc, G, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call uvchksum("f[xy]lf in SIS_C_dynamics", fxlf, fylf, G, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
      call uvchksum("Cor_[uv] in SIS_C_dynamics", Cor_u, Cor_v, G, scale=US%L_T_to_m_s*US%s_to_T)
      call uvchksum("[uv]i in SIS_C_dynamics", ui, vi, G, scale=US%L_T_to_m_s)
    endif

  enddo ! l=1,EVP_steps

  if (CS%debug .or. CS%debug_redundant) &
    call uvchksum("[uv]i end SIS_C_dynamics", ui, vi, G, scale=US%L_T_to_m_s)

  ! Reset the time information in the diag type.
  if (do_hifreq_output) call enable_SIS_averaging(time_int_in, time_end_in, CS%diag)

  ! make averages
  I_sub_steps = 1.0/EVP_steps
!$OMP parallel default(none) shared(isc,iec,jsc,jec,G,fxoc,fxlf,fxic,Cor_u,fxic_d,fxic_t, &
!$OMP                               fxic_s,I_sub_steps,fyoc,fylf,fyic,Cor_v,fyic_d,       &
!$OMP                               fyic_t,fyic_s)
!$OMP do
  do j=jsc,jec ; do I=isc-1,iec
    fxoc(I,j) = fxoc(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxlf(I,j) = fxlf(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic(I,j) = fxic(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    Cor_u(I,j) = Cor_u(I,j) * (G%mask2dCu(I,j) * I_sub_steps)

    fxic_d(I,j) = fxic_d(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic_t(I,j) = fxic_t(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic_s(I,j) = fxic_s(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-1,jec ; do i=isc,iec
    fyoc(i,J) = fyoc(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fylf(i,J) = fylf(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic(i,J) = fyic(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    Cor_v(i,J) = Cor_v(i,J) * (G%mask2dCv(i,J) * I_sub_steps)

    fyic_d(i,J) = fyic_d(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic_t(i,J) = fyic_t(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic_s(i,J) = fyic_s(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
  enddo ; enddo
!$OMP end parallel
  !   Truncate any overly large velocity components.  These final velocities
  ! are the ones that are used for continuity and transport, and hence have
  ! CFL limitations that must be satisfied for numerical stability.
  if ((CS%CFL_trunc > 0.0) .and. (dt_slow > 0.0)) then
    if (len_trim(CS%u_trunc_file) > 0) then
      do j=jsc,jec ; do I=isc-1,iec
        if ((ui(I,j) < ui_min_trunc(I,j)) .or. (ui(I,j) > ui_max_trunc(I,j))) then
          if (mi_u(I,j) > m_neglect) then
            CS%ntrunc = CS%ntrunc + 1
            call write_u_trunc(I, j, ui, u_IC, uo, mis, fxoc, fxic, Cor_u, &
                               PFu, fxat, dt_slow, G, US, CS)
          endif
          if (ui(I,j) < ui_min_trunc(I,j)) then
            ui(I,j) = 0.95 * ui_min_trunc(I,j)
          else
            ui(I,j) = 0.95 * ui_max_trunc(I,j)
          endif
        endif
      enddo ; enddo
    else
      do j=jsc,jec ; do I=isc-1,iec
        if (ui(I,j) < ui_min_trunc(I,j)) then
          ui(I,j) = 0.95 * ui_min_trunc(I,j)
          if (mi_u(I,j) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        elseif (ui(I,j) > ui_max_trunc(I,j)) then
          ui(I,j) = 0.95 * ui_max_trunc(I,j)
          if (mi_u(I,j) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo
    endif
    if (len_trim(CS%v_trunc_file) > 0) then
      do J=jsc-1,jec ; do i=isc,iec
        if ((vi(i,J) < vi_min_trunc(i,J)) .or. (vi(i,J) > vi_max_trunc(i,J))) then
          if (mi_v(i,J) > m_neglect) then
            CS%ntrunc = CS%ntrunc + 1
            call write_v_trunc(i, J, vi, v_IC, vo, mis, fyoc, fyic, Cor_v, &
                               PFv, fyat, dt_slow, G, US, CS)
          endif
          if (vi(i,J) < vi_min_trunc(i,J)) then
            vi(i,J) = 0.95 * vi_min_trunc(i,J)
          else
            vi(i,J) = 0.95*vi_max_trunc(i,J)
          endif
        endif
      enddo ; enddo
    else
      do J=jsc-1,jec ; do i=isc,iec
        if (vi(i,J) < vi_min_trunc(i,J)) then
          vi(i,J) = 0.95 * vi_min_trunc(i,J)
          if (mi_v(i,J) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        elseif (vi(i,J) > vi_max_trunc(i,J)) then
          vi(i,J) = 0.95*vi_max_trunc(i,J)
          if (mi_v(i,J) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo
    endif
  endif

  ! Write out diagnostics associated with the ice dynamics.
  if (query_SIS_averaging_enabled(CS%diag)) then
    if (CS%id_fix>0) call post_SIS_data(CS%id_fix, fxic, CS%diag)
    if (CS%id_fiy>0) call post_SIS_data(CS%id_fiy, fyic, CS%diag)
    if (CS%id_fcx>0) then
      do j=jsc,jec ; do I=isc-1,iec ; diag_val_u(I,j) = Cor_u(I,j)*mi_u(I,j) ; enddo ; enddo
      call post_SIS_data(CS%id_fcx, diag_val_u, CS%diag)
    endif
    if (CS%id_fcy>0) then
      do J=jsc-1,jec ; do i=isc,iec ; diag_val_v(i,J) = Cor_v(i,J)*mi_v(i,J) ; enddo ; enddo
      call post_SIS_data(CS%id_fcy, diag_val_v, CS%diag)
    endif
    if (CS%id_Coru>0) call post_SIS_data(CS%id_Coru, Cor_u, CS%diag)
    if (CS%id_Corv>0) call post_SIS_data(CS%id_Corv, Cor_v, CS%diag)
    if (CS%id_PFu>0) call post_SIS_data(CS%id_PFu, PFu, CS%diag)
    if (CS%id_PFv>0) call post_SIS_data(CS%id_PFv, PFv, CS%diag)
    if (CS%id_fpx>0) then
      do j=jsc,jec ; do I=isc-1,iec ; diag_val_u(I,j) = PFu(I,j)*mi_u(I,j) ; enddo ; enddo
      call post_SIS_data(CS%id_fpx, diag_val_u, CS%diag)
    endif
    if (CS%id_fpy>0) then
      do J=jsc-1,jec ; do i=isc,iec ; diag_val_v(i,J) = PFv(i,J)*mi_v(i,J) ; enddo ; enddo
      call post_SIS_data(CS%id_fpy, diag_val_v, CS%diag)
    endif
    if (CS%id_fwx>0) call post_SIS_data(CS%id_fwx, -fxoc, CS%diag) ! water force on ice
    if (CS%id_fwy>0) call post_SIS_data(CS%id_fwy, -fyoc, CS%diag) ! ...= -ice on water
    if (CS%id_flfx>0) call post_SIS_data(CS%id_flfx, fxlf, CS%diag) ! water force on ice
    if (CS%id_flfy>0) call post_SIS_data(CS%id_flfy, fylf, CS%diag) ! ...= -ice on water
!  The diagnostics of fxat and fyat are supposed to be taken over all partitions
!  (ocean & ice), whereas fxat and fyat here are only averaged over the ice.

    if (CS%id_fix_d>0) call post_SIS_data(CS%id_fix_d, fxic_d, CS%diag)
    if (CS%id_fiy_d>0) call post_SIS_data(CS%id_fiy_d, fyic_d, CS%diag)
    if (CS%id_fix_t>0) call post_SIS_data(CS%id_fix_t, fxic_t, CS%diag)
    if (CS%id_fiy_t>0) call post_SIS_data(CS%id_fiy_t, fyic_t, CS%diag)
    if (CS%id_fix_s>0) call post_SIS_data(CS%id_fix_s, fxic_s, CS%diag)
    if (CS%id_fiy_s>0) call post_SIS_data(CS%id_fiy_s, fyic_s, CS%diag)

    if (CS%id_sigi>0) then
      call find_sigI(mice, ci, CS%str_d, diag_val, G, US, CS)
      call post_SIS_data(CS%id_sigi, diag_val, CS%diag)
    endif
    if (CS%id_sigii>0) then
      call find_sigII(mice, ci, CS%str_t, CS%str_s, diag_val, G, US, CS)
      call post_SIS_data(CS%id_sigii, diag_val, CS%diag)
    endif
    if (CS%id_stren>0) then
      if (CS%project_ci) then
        call find_ice_strength(mice, ci_proj, diag_val, G, US, CS)
      else
        call find_ice_strength(mice, ci, diag_val, G, US, CS)
      endif
      call post_SIS_data(CS%id_stren, diag_val, CS%diag)
    endif
    if (CS%id_stren0>0) then
      call find_ice_strength(mice, ci, diag_val, G, US, CS)
      call post_SIS_data(CS%id_stren0, diag_val, CS%diag)
    endif

    if (CS%id_ui>0) call post_SIS_data(CS%id_ui, ui, CS%diag)
    if (CS%id_vi>0) call post_SIS_data(CS%id_vi, vi, CS%diag)
    if (CS%id_miu>0) call post_SIS_data(CS%id_miu, mi_u, CS%diag)
    if (CS%id_miv>0) call post_SIS_data(CS%id_miv, mi_v, CS%diag)
    if (CS%id_mis>0) call post_SIS_data(CS%id_mis, mis, CS%diag)
    if (CS%id_ci0>0) call post_SIS_data(CS%id_ci0, ci, CS%diag)
    if (CS%id_ci>0)  call post_SIS_data(CS%id_ci, ci_proj, CS%diag)

    if (CS%id_str_d>0) call post_SIS_data(CS%id_str_d, CS%str_d, CS%diag)
    if (CS%id_str_t>0) call post_SIS_data(CS%id_str_t, CS%str_t, CS%diag)
    if (CS%id_str_s>0) call post_SIS_data(CS%id_str_s, CS%str_s, CS%diag)

    if (CS%id_sh_d>0) call post_SIS_data(CS%id_sh_d, sh_Dd, CS%diag)
    if (CS%id_sh_t>0) call post_SIS_data(CS%id_sh_t, sh_Dt, CS%diag)
    if (CS%id_sh_s>0) call post_SIS_data(CS%id_sh_s, sh_Ds, CS%diag)

    if (CS%id_del_sh>0) call post_SIS_data(CS%id_del_sh, del_sh, CS%diag)
    if (CS%id_del_sh_min>0) then
      do j=jsc,jec ; do i=isc,iec
        diag_val(i,j) = del_sh_min_pr(i,j)*pres_mice(i,j)
      enddo ; enddo
      call post_SIS_data(CS%id_del_sh_min, diag_val, CS%diag)
    endif
    if (CS%id_siu>0 .or. CS%id_siv>0 .or. CS%id_sispeed>0) then

      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        if (mis(i,j) > 0.0) then
          siu(i,j) = (ui(I-1,j) + ui(I,j))/2
          siv(i,j) = (vi(i,J-1) + vi(i,J))/2
          sispeed(i,j) = (siu(i,j)*siu(i,j)+siv(i,j)*siv(i,j))**0.5
        else
          siu(i,j) = 0.0; siv(i,j) = 0.0; sispeed(i,j) = 0.0;
        endif
      enddo ; enddo
      if (CS%id_siu>0) call post_SIS_data(CS%id_siu, siu, CS%diag)
      if (CS%id_siv>0) call post_SIS_data(CS%id_siv, siv, CS%diag)
      if (CS%id_sispeed>0) call post_SIS_data(CS%id_sispeed, sispeed, CS%diag)
    endif

  endif

end subroutine SIS_C_dynamics

!> limit_stresses ensures that the input stresses are not larger than could be justified by the ice
!! pressure now, as the ice might have melted or been advected away during the thermodynamic and
!! transport phases, or the ice flow convergence or divergence may have altered the ice concentration.
subroutine limit_stresses(pres_mice, mice, str_d, str_t, str_s, G, US, CS, limit)
  type(SIS_hor_grid_type),            intent(in)    :: G     !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)    :: pres_mice !< The ice internal pressure per
                                                             !! unit column mass [L2 T-2 ~> N m kg-1].
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)    :: mice  !< The mass per unit total area (ice covered
                                                             !! and ice free) of the ice [R Z ~> kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),   intent(inout) :: str_d !< The divergence stress tensor component
                                                             !! [R Z L2 T-2 ~> Pa m].
  real, dimension(SZI_(G),SZJ_(G)),   intent(inout) :: str_t !< The tension stress tensor component
                                                             !! [R Z L2 T-2 ~> Pa m].
  real, dimension(SZIB_(G),SZJB_(G)), intent(inout) :: str_s !< The shearing stress tensor component
                                                             !! [R Z L2 T-2 ~> Pa m].
  type(unit_scale_type),              intent(in)    :: US    !< A structure with unit conversion factors
  type(SIS_C_dyn_CS),                 pointer       :: CS    !< The control structure for this module
  real, optional,                     intent(in)    :: limit !< A factor by which the strength limits
                                                             !! are changed [nondim]

!   This subroutine ensures that the input stresses are not larger than could
! be justified by the ice pressure now, as the ice might have melted or been
! advected away during the thermodynamic and transport phases, or the
! ice flow convergence or divergence may have altered the ice concentration.

  ! Local variables
  real :: pressure  ! The integrated internal ice pressure at a point [R Z L2 T-2 ~> Pa m].
  real :: pres_avg  ! The average of the internal ice pressures around a point [R Z L2 T-2 ~> Pa m].
  real :: sum_area  ! The sum of ocean areas around a vorticity point [L2 ~> m2].
  real :: I_2EC     ! 1/(2*EC), where EC is the yield curve axis ratio [nondim].
  real :: lim       ! A local copy of the factor by which the limits are changed [nondim].
  real :: lim_2     ! The limit divided by 2 [nondim].
!  real :: EC2       ! EC^2, where EC is the yield curve axis ratio [nondim].
!  real :: rescale_str ! A factor by which to rescale the internal stresses [nondim].
!  real :: stress_mag  ! The magnitude of the stress at a point [R Z L2 T-2 ~> Pa m].
!  real :: str_d_q     ! CS%str_d interpolated to a vorticity point [R Z L2 T-2 ~> Pa m].
!  real :: str_t_q     ! CS%str_t interpolated to a vorticity point [R Z L2 T-2 ~> Pa m].

  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  lim = 1.0 ; if (present(limit)) lim = limit
  I_2EC = 0.0 ; if (CS%EC > 0.0) I_2EC = (0.5*lim) / CS%EC
  lim_2 = 0.5 * lim

  ! The rescaling here is done separately for each component.
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    pressure = pres_mice(i,j)*mice(i,j)
    if (str_d(i,j) < -pressure) str_d(i,j) = -pressure
    if (CS%EC*str_t(i,j) > lim_2*pressure) str_t(i,j) = I_2EC*pressure
    if (CS%EC*str_t(i,j) < -lim_2*pressure) str_t(i,j) = -I_2EC*pressure
  enddo ; enddo
  do J=jsc-1,jec ; do I=isc-1,iec
    if (CS%weak_coast_stress) then
      sum_area = (G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i,j+1) + G%areaT(i+1,j))
    else
      sum_area = (G%mask2dT(i,j)*G%areaT(i,j) + G%mask2dT(i+1,j+1)*G%areaT(i+1,j+1)) + &
                 (G%mask2dT(i,j+1)*G%areaT(i,j+1) + G%mask2dT(i+1,j)*G%areaT(i+1,j))
    endif
    pres_avg = 0.0
    if (sum_area > 0.0) &
      pres_avg = ((G%areaT(i,j)   * (pres_mice(i,j)*mice(i,j)) + &
                   G%areaT(i+1,j+1)*(pres_mice(i+1,j+1)*mice(i+1,j+1))) + &
                  (G%areaT(i+1,j) * (pres_mice(i+1,j)*mice(i+1,j)) + &
                   G%areaT(i,j+1) * (pres_mice(i,j+1)*mice(i,j+1)))) / sum_area

    if (CS%EC*str_s(I,J) > lim_2*pres_avg) str_s(I,J) = I_2EC*pres_avg
    if (CS%EC*str_s(I,J) < -lim_2*pres_avg) str_s(I,J) = -I_2EC*pres_avg
  enddo ; enddo

!    This commented out version seems to work, but is not obviously better than
! treating each component separately, and the later is simpler.
!  EC2 = CS%EC**2
!  do J=jsc-1,jec ; do I=isc-1,iec
!    ! Rescale str_s based on interpolated values of str_d and str_t, which works
!    ! because the ice strengths are also interpolated.
!    if (CS%weak_coast_stress) then
!      sum_area = (G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i,j+1) + G%areaT(i+1,j))
!    else
!      sum_area = (G%mask2dT(i,j)*G%areaT(i,j) + G%mask2dT(i+1,j+1)*G%areaT(i+1,j+1)) + &
!                 (G%mask2dT(i,j+1)*G%areaT(i,j+1) + G%mask2dT(i+1,j)*G%areaT(i+1,j))
!    endif
!    if (sum_area > 0.0) then
!      pres_avg = ((G%areaT(i,j)   * (pres_mice(i,j)*mice(i,j)) + &
!                   G%areaT(i+1,j+1)*(pres_mice(i+1,j+1)*mice(i+1,j+1))) + &
!                  (G%areaT(i+1,j) * (pres_mice(i+1,j)*mice(i+1,j)) + &
!                   G%areaT(i,j+1) * (pres_mice(i,j+1)*mice(i,j+1)))) / sum_area
!      if (pres_avg <= 0.0) then
!        str_s(I,J) = 0.0
!      else
!        str_d_q = 0.25 * ((str_d(i,j) + str_d(i+1,j+1)) + &
!                          (str_d(i+1,j) + str_d(i,j+1)))
!        str_t_q = 0.25 * ((str_t(i,j) + str_t(i+1,j+1)) + &
!                          (str_t(i+1,j) + str_t(i,j+1)))
!        ! The factor of 2 here arises because of the definitions of the str_#.
!        stress_mag = 2.0 * sqrt(min(0.0, str_d_q + 0.5*pres_avg)**2 + &
!                                EC2 * (str_s(I,J)**2 + str_t_q**2))
!        if ((stress_mag > pres_avg) .and. (G%mask2dBu(I,j)>0.0)) &
!          str_s(I,J) = str_s(I,J) * (pres_avg / stress_mag)
!      endif
!    endif
!  enddo ; enddo

!  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
!    ! Rescale str_d and str_t without regard to the values of str_s.
!    pressure = pres_mice(i,j)*mice(i,j)
!    if (pressure <= 0.0) then
!      str_d(i,j) = 0.0 ; str_t(i,j) = 0.0
!    elseif (str_d(i,j) < -0.5*pressure) then
!    ! The factor of 2 here arises because of the definitions of str_d and str_t.
!      stress_mag = 2.0 * sqrt((str_d(i,j) + 0.5*pressure)**2 + &
!                              EC2 * str_t(i,j)**2)
!      if (stress_mag > pressure) then
!        rescale_str = pressure / stress_mag
!        str_d(i,j) = rescale_str * (str_d(i,j) + 0.5*pressure) - &
!                        0.5*pressure
!        str_t(i,j) = str_t(i,j) * rescale_str
!      endif
!    elseif (CS%EC*abs(str_t(i,j)) > 0.5*pressure) then
!      ! Only reduce excessively large values of str_t, but do not increase the
!      ! magnitude of str_d.
!      str_t(i,j) = sign(I_2EC*pressure, str_t(i,j))
!    endif
!  enddo ; enddo

end subroutine limit_stresses

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> find_sigI finds the first stress invariant
subroutine find_sigI(mi, ci, str_d, sigI, G, US, CS)
  type(SIS_hor_grid_type),          intent(in)  :: G     !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: mi    !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: ci    !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: str_d !< The divergence stress tensor component
                                                         !! [R Z L2 T-2 ~> Pa m].
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: sigI  !< The first stress invariant [nondim]
  type(unit_scale_type),            intent(in)  :: US    !< A structure with unit conversion factors
  type(SIS_C_dyn_CS),               pointer     :: CS    !< The control structure for this module

  real, dimension(SZI_(G),SZJ_(G)) :: &
    strength ! The ice strength [R Z L2 T-2 ~> Pa m].
  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call find_ice_strength(mi, ci, strength, G, US, CS)

  do j=jsc,jec ; do i=isc,iec
    sigI(i,j) = 0.0
    if (strength(i,j) > 0.0) sigI(i,j) = 2.0 * str_d(i,j) / strength(i,j)
  enddo ; enddo

end subroutine find_sigI

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> find_sigII finds the second stress invariant
subroutine find_sigII(mi, ci, str_t, str_s, sigII, G, US, CS)
  type(SIS_hor_grid_type),            intent(in)  :: G     !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: mi    !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: ci    !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: str_t !< The tension stress tensor component
                                                           !! [R Z L2 T-2 ~> Pa m].
  real, dimension(SZIB_(G),SZJB_(G)), intent(in)  :: str_s !< The shearing stress tensor component
                                                           !! [R Z L2 T-2 ~> Pa m].
  real, dimension(SZI_(G),SZJ_(G)),   intent(out) :: sigII !< The second stress invariant [nondim].
  type(unit_scale_type),              intent(in)  :: US    !< A structure with unit conversion factors
  type(SIS_C_dyn_CS),                 pointer     :: CS    !< The control structure for this module

  real, dimension(SZI_(G),SZJ_(G)) :: &
    strength ! The ice strength [R Z L2 T-2 ~> Pa m].
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    str_s_ss ! Str_s divided by the sum of the neighboring ice strengths [nondim].
  real :: strength_sum  ! The sum of the 4 neighboring strengths times areas [R Z L4 T-2 ~> Pa m3].
  real :: sum_area   ! The sum of ocean areas around a vorticity point [L2 ~> m2].
  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    strength(i,j) = mi(i,j) * CS%p0_rho*exp(-CS%c0*max(1.0-ci(i,j),0.0))
  enddo ; enddo

  do J=jsc-1,jec ; do I=isc-1,iec
    str_s_ss(I,J) = 0.0
    sum_area = (G%mask2dT(i,j)*G%areaT(i,j) + G%mask2dT(i+1,j+1)*G%areaT(i+1,j+1)) + &
               (G%mask2dT(i,j+1)*G%areaT(i,j+1) + G%mask2dT(i+1,j)*G%areaT(i+1,j))
    if (sum_area > 0.0) then
      strength_sum = (G%areaT(i,j)*strength(i,j) + G%areaT(i+1,j+1)*strength(i+1,j+1)) + &
                     (G%areaT(i+1,j)*strength(i+1,j) + G%areaT(i,j+1)*strength(i,j+1))
      if (strength_sum > 0.0) str_s_ss(I,J) = str_s(I,J) * sum_area / strength_sum
    endif
  enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    sigII(i,j) = 0.0

    ! This distributes str_s according to the strength of the neighboring cells.
    if (strength(i,j) > 0.0) &
      sigII(i,j) = sqrt((2.0*str_t(i,j)/strength(i,j))**2 + &
                        (0.5*((str_s_ss(I,J) + str_s_ss(I-1,J-1)) + &
                              (str_s_ss(I-1,J) + str_s_ss(I,J-1))))**2 )
  enddo ; enddo

end subroutine find_sigII

!> Computes basal stress Tbu coefficients (landfast ice)
!!
!! Lemieux, J. F., B. Tremblay, F. Dupont, M. Plante, G.C. Smith, D. Dumont (2015).
!! A basal stress parameterization form modeling landfast ice, J. Geophys. Res.
!! Oceans, 120, 3157-3173.
!!
!! Lemieux, J. F., F. Dupont, P. Blain, F. Roy, G.C. Smith, G.M. Flato (2016).
!! Improving the simulation of landfast ice by combining tensile strength and a
!! parameterization for grounded ridges, J. Geophys. Res. Oceans, 121.
!!
!! author: JF Lemieux, Philippe Blain (ECCC)
!!
!! note: Tb_u and Tb_v are parts of the Cb as defined in Lemieux et al. 2015 and 2016.
!!
subroutine basal_stress_coeff_C(G, mi, ci, sea_lev, CS)

  type(SIS_hor_grid_type),            intent(in)  :: G     !< The horizontal grid type
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: mi    !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: ci    !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: sea_lev !< Sea level [Z ~> m]
  type(SIS_C_dyn_CS),                 pointer     :: CS    !< The control structure for this module

  real :: &
         h_u, & ! volume per unit area of ice at u location (mean thickness) [Z ~> m]
         h_v, & ! volume per unit area of ice at v location (mean thickness) [Z ~> m]
         hw_u, & ! water depth at u location [Z ~> m]
         hw_v, & ! water depth at v location [Z ~> m]
         hc_u, & ! critical thickness at u location [Z ~> m]
         hc_v    ! critical thickness at v location [Z ~> m]

  integer :: i, j, isc, iec, jsc, jec
  real :: ci_u ! Concentration at u-points [nondim]
  real :: ci_v ! Concentration at u-points [nondim]
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  ! Compute the term (h_u - h_cu) * exp(-C(1-A_u)) (at u points)
  do j=jsc,jec
    do I=isc-1,iec
      hw_u = min(G%bathyT(i,j) + sea_lev(i,j), G%bathyT(i+1,j) + sea_lev(i+1,j))
      ci_u = max(ci(i,j), ci(i+1,j))
      if (ci_u > 0.01 .and. G%mask2dCu(I,j) > 0.0 .and. hw_u < CS%lemieux_threshold_hw) then
        h_u = max(mi(i,j), mi(i+1,j))/CS%Rho_ice

        ! 1- calculate critical thickness
        hc_u = hw_u * ci_u / CS%lemieux_k1

        ! 2- calculate basal stress factor
        CS%Tb_u(I,j) = CS%lemieux_k2 * MAX(0.0,                                   &
                 (h_u - hc_u)) * exp(-CS%lemieux_alphab*(1.0-ci_u))
      else
        CS%Tb_u(I,j) = 0.0
      endif
    enddo
  enddo
!
! Compute the term (h_v - h_cv) * exp(-C(1-A_v)) (at v points)
!
  do J=jsc-1,jec
    do i=isc,iec
      hw_v = min(G%bathyT(i,j) + sea_lev(i,j), G%bathyT(i,j+1) + sea_lev(i,j+1))
      ci_v = max(ci(i,j), ci(i,j+1))
      if (ci_v > 0.01 .and. G%mask2dCv(i,J) > 0.0 .and. hw_v < CS%lemieux_threshold_hw) then
        h_v = max(mi(i,j), mi(i,j+1))/CS%Rho_ice

        ! 1- calculate critical thickness
        hc_v = hw_v * ci_v / CS%lemieux_k1

        ! 2- calculate basal stress factor
        CS%Tb_v(i,J) = CS%lemieux_k2 * MAX(0.0,                                   &
                  (h_v - hc_v)) * exp(-CS%lemieux_alphab*(1.0-ci_v))
      else
        CS%Tb_v(i,J) = 0.0
      endif
    enddo
  enddo
!          call uvchksum("Tb_[uv] before SIS_C_dynamics", CS%Tb_u, CS%Tb_v, G, &
!                         halos=1, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)

end subroutine basal_stress_coeff_C

!> Computes basal stress Tbu coefficients (landfast ice)
!!
!! Computes seabed (basal) stress factor Tbu (landfast ice) based on
!! probability of contact between the ITD and the seabed. The water depth
!! could take into account variations of the SSH. In the simplest
!! formulation, hwater is simply the value of the bathymetry. To calculate
!! the probability of contact, it is assumed that the bathymetry follows
!! a normal distribution with sigma_b = 2.5d0. An improvement would
!! be to provide the distribution based on high resolution data.
!!
!! Dupont, F., D. Dumont, J.F. Lemieux, E. Dumas-Lefebvre, A. Caya (2022).
!! A probabilistic seabed-ice keel interaction model, The Cryosphere, 16,
!! 1963-1977.
!!
!! authors: D. Dumont, J.F. Lemieux, E. Dumas-Lefebvre, F. Dupont
!!
subroutine basal_stress_coeff_itd(G, IG, IST, sea_lev, CS)

  type(SIS_hor_grid_type),            intent(in)  :: G     !< The horizontal grid type
  type(ice_grid_type),                intent(in)  :: IG    !< The sea-ice specific grid type
  type(ice_state_type),               intent(in)  :: IST   !< A type describing the state of the sea ice
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: sea_lev !< Sea level [Z ~> m]
  type(SIS_C_dyn_CS),                 pointer     :: CS    !< The control structure for this module

  real, dimension(CS%ncat_i) :: & ! log-normal for ice thickness
           x_k, & ! center of thickness categories (m)
           g_k, & ! probability density function (thickness, 1/m)
           P_x    ! probability for each thickness category

  real, dimension(CS%ncat_b) :: & ! normal dist for bathymetry
           y_n, & ! center of bathymetry categories (m)
           b_n, & ! probability density function (bathymetry, 1/m)
           P_y    ! probability for each bathymetry category

  integer, dimension(CS%ncat_b) :: tmp
  real, dimension(0:IG%CatIce) :: hin_max   ! category limits (m)

  logical, dimension (CS%ncat_b) :: gt    ! result of comparison between x_k and y_n

  real, dimension(IG%CatIce) :: vcat   ! category limits (m)
  real, dimension(IG%CatIce) :: acat   ! category limits (m)

  !  parameters for PDFs
  real :: wid_i     ! (m)
  real :: wid_b     ! (m)
  real :: mu_i      ! [nondim]
  real :: sigma_i   ! [nondim]
  real :: mu_b      ! (m)
  real :: m_i       ! (m)
  real :: v_i       ! (m2)

  real, dimension(CS%ncat_i):: tb_tmp
  real, dimension (SZI_(G),SZJ_(G)):: Tbt ! seabed stress factor at t point [R Z L T-2 -> kg m-1 s-2]
  real :: atot       ! [nondim]
  real :: x_kmax     ! thickest ice? [m]
  real :: cut        ! thinnest ice thickness [m]
  real :: rho_ice    ! ice density [R ~> kg m-3]
  real :: rho_water  ! water density [R ~> kg m-3]
  real :: pi         ! [nondim]
  integer :: i, ii, j, isc, iec, jsc, jec, k, n, ncat
  real :: ci_u       ! Concentration at u-points [nondim]
  real :: ci_v       ! Concentration at u-points [nondim]

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ncat = IG%CatIce
  call get_SIS2_thermo_coefs(IST%ITV, rho_ice=rho_ice, rho_water=rho_water)
  pi = 4.0*atan(1.0)

! a note regarding hi_min and hin_max(0):
! both represent a minimum ice thickness.  hin_max(0) is
! intended to be used for particular numerical implementations
! of category conversions in the ice thickness distribution.
! hi_min is a more general purpose parameter, but is specifically
! for maintaining stability in the thermodynamics.
! hin_max(0) = 0.1 m for the delta function itd
! hin_max(0) = 0.0 m for linear remapping
!
! Also note that the upper limit on the thickest category
! is only used for the linear remapping scheme
! and it is not a true upper limit on the thickness
  hin_max(0)=0.0
  do k=1,nCat
    hin_max(k) = IG%mH_cat_bound(k)/rho_ice
  end do

  Tbt=0.0

  do j=jsc-1,jec+1
    do i=isc-1,iec+1

      acat(1:ncat) = IST%part_size(i,j,1:ncat)
      atot = sum(acat)
      vcat(1:ncat) = IST%mH_ice(i,j,1:ncat) / rho_ice
      m_i = sum(vcat)

      if (atot > 0.05 .and. m_i > CS%basal_stress_min_thick .and. &
          sea_lev(i,j) < CS%basal_stress_max_depth) then

        mu_b = G%bathyT(i,j) + sea_lev(i,j)         ! (hwater) mean of PDF (normal dist) bathymetry
        wid_i = CS%basal_stress_max_depth/CS%ncat_i    ! width of ice categories
        wid_b = 6.0*CS%sigma_b(i,j)/CS%ncat_b          ! width of bathymetry categories (6 sigma_b = 2x3 sigma_b)

        x_k(:) = (/( wid_i*( real(ii) - 0.5 ), ii=1, CS%ncat_i )/)
        y_n(:) = (/( ( mu_b - 3.0*CS%sigma_b(i,j) ) + (real(ii) - 0.5) * (6.0*CS%sigma_b(i,j)/CS%ncat_b), &
                         ii=1, CS%ncat_b )/)

        v_i=0.0
        do n =1, ncat
          v_i = v_i + vcat(n)**2 / (max(acat(n), CS%puny))
        enddo
        v_i = v_i - m_i**2

        ! parameters for the log-normal
        mu_i    = log(m_i/(CS%onemeter * sqrt(1.0 + v_i/m_i**2)))
        sigma_i = sqrt(log(1.0 + v_i/m_i**2))

        ! max thickness associated with percentile of log-normal PDF
        ! x_kmax=x997 was obtained from an optimization procedure (Dupont et al. 2022)

        if (sigma_i > 0) then
          x_kmax = CS%onemeter * exp(mu_i + sqrt(2.0*sigma_i)*CS%basal_stress_cutoff)

          ! Set x_kmax to hlev of the last category where there is ice
          ! when there is no ice in the last category
          cut = x_k(CS%ncat_i)
          do n = ncat,-1,1
            if (acat(n) < CS%puny) then
              cut = hin_max(n-1)
            else
              exit
            endif
          enddo
          x_kmax = min(cut, x_kmax)

          g_k(:) = exp(-(log(x_k(:)/CS%onemeter) - mu_i) ** 2 / (2.0 * sigma_i ** 2)) / &
                   (x_k(:) * sigma_i * sqrt(2.0 * pi))

          b_n(:)  = exp(-(y_n(:) - mu_b) ** 2 / (2.0 * CS%sigma_b(i,j) ** 2)) / (CS%sigma_b(i,j) * sqrt(2.0*pi))

          P_x(:) = g_k(:) * wid_i
          P_y(:) = b_n(:) * wid_b

          do n =1, CS%ncat_i
            if (x_k(n) > x_kmax) P_x(n)=0.0
          enddo

          ! calculate Tb factor at t-location
          do n=1, CS%ncat_i
            gt(:) = (y_n(:) <= rho_ice*x_k(n)/rho_water)
            tmp(:) = merge(1,0,gt(:))
            ii = sum(tmp)
            if (ii == 0) then
              tb_tmp(n) = 0.0
            else
              tb_tmp(n) = max(CS%basal_stress_mu_s * G%g_Earth * P_x(n) * &
                          sum(P_y(1:ii)*(rho_ice*x_k(n) - rho_water*y_n(1:ii))), 0.0)
            endif
          enddo
          Tbt(i,j) = sum(tb_tmp) * exp(-CS%lemieux_alphab * (1.0 - atot))
        else
          Tbt(i,j) = 0.0
        endif
      endif
    enddo
  enddo

  do j=jsc,jec
    do I=isc-1,iec
      ! convert quantities to u-location
      CS%Tb_u(I,j)  = max(Tbt(i,j),Tbt(i+1,j))
    enddo
  enddo
  do J=jsc-1,jec
    do i=isc,iec
      ! convert quantities to v-location
      CS%Tb_v(i,J)  = max(Tbt(i,j),Tbt(i,j+1))
    enddo
  enddo

end subroutine basal_stress_coeff_itd

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_C_dyn_register_restarts allocates and registers any variables for the
!!   SIS C-grid dynamics module that need to be included in the restart files.
subroutine SIS_C_dyn_register_restarts(HI, param_file, CS, US, Ice_restart)
  type(hor_index_type),    intent(in) :: HI    !< The horizontal index type describing the domain
  type(param_file_type),   intent(in) :: param_file !< A structure to parse for run-time parameters
  type(SIS_C_dyn_CS),      pointer    :: CS    !< The control structure for this module
  type(unit_scale_type),   intent(in) :: US    !< A structure with unit conversion factors
  type(SIS_restart_CS),    pointer    :: Ice_restart  !< The control structure for the ice restarts

!   This subroutine registers the restart variables associated with the
! the ice dynamics.

  integer :: isd, ied, jsd, jed, id
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed

  if (associated(CS)) then
    call SIS_error(WARNING, "SIS_C_dyn_register_restarts called with an "//&
                            "associated control structure.")
    return
  endif
  allocate(CS)

  call safe_alloc(CS%str_d, isd, ied, jsd, jed)
  call safe_alloc(CS%str_t, isd, ied, jsd, jed)
  call safe_alloc(CS%str_s, HI%IsdB, HI%IedB, HI%JsdB, HI%JedB)

  if (associated(Ice_restart)) then
    call register_restart_field(Ice_restart, 'str_d', CS%str_d, mandatory=.false., &
                                units="Pa m", conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
    call register_restart_field(Ice_restart, 'str_t', CS%str_t, mandatory=.false., &
                                units="Pa m", conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
    if (HI%symmetric) then
      call register_restart_field(Ice_restart, 'sym_str_s', CS%str_s, position=CORNER, mandatory=.false., &
                                  units="Pa m", conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
    else
      call register_restart_field(Ice_restart, 'str_s', CS%str_s, position=CORNER, mandatory=.false., &
                                  units="Pa m", conversion=US%RZ_to_kg_m2*US%L_T_to_m_s**2)
    endif
  endif
end subroutine SIS_C_dyn_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_C_dyn_read_alt_restarts reads in alternative variables for the SIS C-grid dynamics module
!!   that might have been in the restart file, specifically dealing with changing between symmetric
!!   and non-symmetric memory restart files.  It also handles any changes in dimensional rescaling
!!   of these variables between what is stored in the restart file and what is done for the current
!!   run segment.
subroutine SIS_C_dyn_read_alt_restarts(CS, G, US, Ice_restart, restart_dir)
  type(SIS_C_dyn_CS),      pointer    :: CS    !< The control structure for this module
  type(SIS_hor_grid_type), intent(in) :: G   !< The horizontal grid type
  type(unit_scale_type),   intent(in) :: US  !< A structure with unit conversion factors
  type(SIS_restart_CS),    pointer    :: Ice_restart !< The control structure for the ice restarts
  character(len=*),        intent(in) :: restart_dir !< The directory in which to find the restart files

  ! These are temporary variables that will be used only here for reading and
  ! then discarded.
  real, allocatable, target, dimension(:,:) :: str_tmp
  type(MOM_domain_type),   pointer :: domain_tmp => NULL()
  real :: stress_rescale
  logical :: read_values
  integer :: i, j, id

  if (.not.associated(Ice_restart)) return
  if (G%symmetric .and. (.not.query_initialized(Ice_restart, 'sym_str_s'))) then

    call clone_MOM_domain(G%domain, domain_tmp, symmetric=.false., &
                          domain_name="ice temporary domain")
    allocate(str_tmp(G%isd:G%ied, G%jsd:G%jed), source=0.0)

    call only_read_from_restarts(Ice_restart, 'str_s', str_tmp, domain_tmp, position=CORNER, &
                                 directory=restart_dir, success=read_values)
    if (read_values) then
      ! The non-symmetric variant of this variable has been successfully read.
      call pass_var(str_tmp, domain_tmp, position=CORNER)
      do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
        CS%str_s(I,J) = str_tmp(I,J)
      enddo ; enddo
    endif

  elseif ((.not.G%symmetric) .and. (.not.query_initialized(Ice_restart, 'str_s'))) then

    call clone_MOM_domain(G%domain, domain_tmp, symmetric=.true., &
                          domain_name="ice temporary domain")
    allocate(str_tmp(G%isd-1:G%ied, G%jsd-1:G%jed), source=0.0)

    call only_read_from_restarts(Ice_restart, 'sym_str_s', str_tmp, domain_tmp, position=CORNER, &
                                 directory=restart_dir, success=read_values)
    if (read_values) then
      ! The symmetric variant of this variable has been successfully read.
      do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
        CS%str_s(I,J) = str_tmp(I,J)
      enddo ; enddo
    endif
  endif

  if (allocated(str_tmp)) deallocate(str_tmp)
  if (associated(domain_tmp)) then ; deallocate(domain_tmp%mpp_domain) ; deallocate(domain_tmp) ; endif

  ! Now redo the dimensional rescaling of the stresses if necessary.
  if (US%s_to_T_restart*US%m_to_L_restart*US%kg_m3_to_R_restart*US%m_to_Z_restart /= 0.0) then
    stress_rescale = US%s_to_T_restart**2 / &
                     (US%kg_m3_to_R_restart * US%m_to_Z_restart * US%m_to_L_restart**2)
    do J=G%jsc-1,G%jec ; do I=G%isc-1,G%iec
      CS%str_s(I,J) = stress_rescale * CS%str_s(I,J)
      if (abs(CS%str_s(I,J)) < CS%str_underflow) CS%str_s(I,J) = 0.0
    enddo ; enddo
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      CS%str_d(i,j) = stress_rescale * CS%str_d(i,j)
      CS%str_t(i,j) = stress_rescale * CS%str_t(i,j)
      if (abs(CS%str_d(i,j)) < CS%str_underflow) CS%str_d(i,j) = 0.0
      if (abs(CS%str_t(i,j)) < CS%str_underflow) CS%str_t(i,j) = 0.0
    enddo ; enddo
  endif

end subroutine SIS_C_dyn_read_alt_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> write_u_trunc is used to record the location of any pseudo-zonal velocity
!! truncations and related fields.
subroutine write_u_trunc(I, j, ui, u_IC, uo, mis, fxoc, fxic, Cor_u, PFu, fxat, &
                         dt_slow, G, US, CS)
  integer,                           intent(in) :: I    !< The i-index of the column to report on
  integer,                           intent(in) :: j    !< The j-index of the column to report on
  type(SIS_hor_grid_type),           intent(in) :: G    !< The horizontal grid type
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: ui   !< The zonal ice velocity [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: u_IC !< The initial zonal ice velocity [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: uo   !< The zonal ocean velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: mis  !< The mass of ice and snow per unit ocean area [R Z ~> kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: fxoc !< The zonal ocean-to-ice force [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: fxic !< The ice internal force [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: Cor_u !< The zonal Coriolis acceleration [L T-2 ~> m s-2].
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: PFu  !< The zonal Pressure force acceleration [L T-2 ~> m s-2].
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: fxat !< The zonal wind stress [R Z L T-2 ~> Pa].
  real,                              intent(in) :: dt_slow !< The slow ice dynamics timestep [T ~> s].
  type(unit_scale_type),             intent(in) :: US   !< A structure with unit conversion factors
  type(SIS_C_dyn_CS),                pointer    :: CS   !< The control structure for this module

  real :: dt_mi, dt_usc, u_scale, CFL
  real :: m_neglect  ! A tiny mass per unit area [R Z ~> kg m-2].
  real, parameter :: H_subroundoff = 1e-30 ! A negligible ice thickness [m]
  integer :: file
  integer :: yr, mo, day, hr, minute, sec, yearday

  m_neglect = H_subroundoff*US%m_to_Z*CS%Rho_ice

  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

  ! Open up the file for output if this is the first call.
    if (CS%u_file < 0) then
      if (len_trim(CS%u_trunc_file) < 1) return
      call open_ASCII_file(CS%u_file, trim(CS%u_trunc_file), &
          action=APPEND_FILE)
      if (CS%u_file < 0) then
        call SIS_error(NOTE, 'Unable to open file '//trim(CS%u_trunc_file)//'.')
        return
      endif
    endif
    file = CS%u_file

    if (ui(I,j) > 0.0) then
      CFL = (ui(I,j) * (dt_slow*G%dy_Cu(I,j))) / G%areaT(i,j)
    else
      CFL = (ui(I,j) * (dt_slow*G%dy_Cu(I,j))) / G%areaT(i+1,j)
    endif

    u_scale = US%L_T_to_m_s
    dt_usc = dt_slow * u_scale

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'("Time ",i5,i4,F6.2," U-trunc at ",I4,": ",2(I3), &
        & " (",F7.2," E "F7.2," N) u = ",ES10.3," (CFL ",ES9.2,") was ",ES10.3," dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), I, j, &
        G%geoLonCu(I,j), G%geoLatCu(I,j), u_scale*ui(I,j), CFL, u_scale*u_IC(I,j), US%T_to_s*dt_slow

    dt_mi = dt_usc / (0.5*(mis(i,j) + mis(i+1,j)) + m_neglect)

    write (file, '("ui, uo, dui = ", 3ES11.3, " ;  mice+snow = ",2ES11.3)') &
           u_scale*ui(I,j), u_scale*uo(I,j), u_scale*(ui(I,j) - u_IC(I,j)), &
           US%RZ_to_kg_m2*mis(i,j), US%RZ_to_kg_m2*mis(i+1,j)

    write (file, '("U change due to fxat, fxoc, fxic, Cor_u, PFu = ", 5ES11.3, " sum = ",ES11.3)') &
      fxat(I,j)*dt_mi, -fxoc(I,j)*dt_mi, fxic(I,j)*dt_mi, Cor_u(I,j)*dt_usc, PFu(I,j)*dt_usc, &
      (fxat(I,j) - fxoc(I,j) + fxic(I,j))*dt_mi + (Cor_u(I,j) + PFu(I,j))*dt_usc

    call flush(file)
  endif

end subroutine write_u_trunc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> write_v_trunc is used to record the location of any pseudo-meridional velocity
!! truncations and related fields.
subroutine write_v_trunc(i, J, vi, v_IC, vo, mis, fyoc, fyic, Cor_v, PFv, fyat, &
                         dt_slow, G, US, CS)
  integer,                           intent(in) :: i    !< The i-index of the column to report on
  integer,                           intent(in) :: J    !< The j-index of the column to report on
  type(SIS_hor_grid_type),           intent(in) :: G    !< The horizontal grid type
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: vi   !< The meridional ice velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: v_IC !< The initial meridional ice velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: vo   !< The meridional ocean velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: mis  !< The mass of ice and snow per unit ocean area [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: fyoc !< The meridional ocean-to-ice force [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: fyic !< The ice internal force [R Z L T-2 ~> Pa].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: Cor_v !< The meridional Coriolis acceleration [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: PFv  !< The meridional pressure force acceleration [L T-2 ~> m s-2].
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: fyat !< The meridional wind stress [R Z L T-2 ~> Pa].
  real,                              intent(in) :: dt_slow !< The slow ice dynamics timestep [T ~> s].
  type(unit_scale_type),             intent(in) :: US   !< A structure with unit conversion factors
  type(SIS_C_dyn_CS),                pointer    :: CS   !< The control structure for this module

  real :: dt_mi, dt_usc, u_scale, CFL
  real :: m_neglect  ! A tiny mass per unit area [R Z ~> kg m-2].
  real, parameter :: H_subroundoff = 1e-30 ! A negligible ice thickness [m].
  integer :: file
  integer :: yr, mo, day, hr, minute, sec, yearday

  m_neglect = H_subroundoff*US%m_to_Z*CS%Rho_ice

  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

  ! Open up the file for output if this is the first call.
    if (CS%v_file < 0) then
      if (len_trim(CS%v_trunc_file) < 1) return
      call open_ASCII_file(CS%v_file, trim(CS%v_trunc_file), &
          action=APPEND_FILE)
      if (CS%v_file < 0) then
        call SIS_error(NOTE, 'Unable to open file '//trim(CS%v_trunc_file)//'.')
        return
      endif
    endif
    file = CS%v_file

    if (vi(i,J) > 0.0) then
      CFL = (vi(i,J) * (dt_slow*G%dx_Cv(i,J))) / G%areaT(i,j)
    else
      CFL = (vi(i,J) * (dt_slow*G%dx_Cv(i,J))) / G%areaT(i,j+1)
    endif

    u_scale = US%L_T_to_m_s
    dt_usc = dt_slow * u_scale

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'("Time ",i5,i4,F6.2," V-trunc at ",I4,": ",2(I3), &
        & " (",F7.2," E ",F7.2," N) v = ",ES10.3," (CFL ",ES9.2,") was ",ES10.3," dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), i, J, &
        G%geoLonCv(i,J), G%geoLatCv(i,J), u_scale*vi(i,J), CFL, u_scale*v_IC(i,J), US%T_to_s*dt_slow

    dt_mi = dt_usc / (0.5*(mis(i,j) + mis(i,j+1)) + m_neglect)

    write (file, '("vi, vo, dvi = ", 3ES11.3, " ;  mice+snow = ",2ES11.3)') &
           u_scale*vi(i,J), u_scale*vo(i,J), u_scale*(vi(i,J) - v_IC(i,J)), &
           US%RZ_to_kg_m2*mis(i,j), US%RZ_to_kg_m2*mis(i,j+1)

    write (file, '("V change due to fyat, fyoc, fyic, Cor_v, PFv = ", 5ES11.3, " sum = ",ES11.3)') &
      fyat(i,J)*dt_mi, -fyoc(i,J)*dt_mi, fyic(i,J)*dt_mi, Cor_v(i,J)*dt_usc, PFv(i,J)*dt_usc, &
      (fyat(i,J) - fyoc(i,J) + fyic(i,J))*dt_mi + (Cor_v(i,J) + PFv(i,J))*dt_usc

    call flush(file)
  endif

end subroutine write_v_trunc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_C_dyn_end deallocates the memory associated with this module.
subroutine SIS_C_dyn_end(CS)
  type(SIS_C_dyn_CS), pointer :: CS    !< The control structure for this module

  deallocate(CS%str_d) ; deallocate(CS%str_t) ; deallocate(CS%str_s)
  if (associated(CS%Tb_u)) deallocate(CS%Tb_u)
  if (associated(CS%Tb_v)) deallocate(CS%Tb_v)
  if (associated(CS%sigma_b)) deallocate(CS%sigma_b)

  deallocate(CS)
end subroutine SIS_C_dyn_end

end module SIS_dyn_cgrid
