module SIS_dyn_cgrid
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
! C-grid SEA ICE DYNAMICS using ELASTIC-VISCOUS-PLASTIC RHEOLOGY adapted from  !
! Hunke and Dukowicz (JPO 1997, H&D hereafter) with some derivation from SIS1  !
! and with guidance from the C-grid implementations of sea-ice in MITgcm as    !
! documented in MITgcm user notes by Martin Losch and in LIM3 by S. Bouillon   !
! et al. (Ocean Modelling, 2009 & 2013). This code initially written by        !
! Robert Hallberg in 2013.                                                     !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_diag_mediator, only : post_SIS_data, SIS_diag_ctrl
use SIS_diag_mediator, only : query_SIS_averaging_enabled, enable_SIS_averaging
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_debugging,     only : chksum, Bchksum, hchksum, vec_chksum_C
use SIS_debugging,     only : check_redundant_B, check_redundant_C, vec_chksum_C
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, NOTE, SIS_mesg=>MOM_mesg
use MOM_file_parser,  only : get_param, log_param, read_param, log_version, param_file_type
use MOM_domains,      only : pass_var, pass_vector, CGRID_NE, CORNER, pe_here
use MOM_hor_index,    only : hor_index_type
use MOM_io, only : open_file
use MOM_io, only : APPEND_FILE, ASCII_FILE, MULTIPLE, SINGLE_FILE
use SIS_hor_grid, only : SIS_hor_grid_type
use fms_io_mod,       only : register_restart_field, restart_file_type
use mpp_domains_mod,  only : domain2D
use time_manager_mod, only : time_type, set_time, operator(+), operator(-)
use time_manager_mod, only : set_date, get_time, get_date

implicit none ; private

#include <SIS2_memory.h>

public :: SIS_C_dyn_init, SIS_C_dynamics, SIS_C_dyn_end, SIS_C_dyn_register_restarts

type, public :: SIS_C_dyn_CS ; private
  real, allocatable, dimension(:,:) :: &
    str_t, &  ! The tension stress tensor component, in Pa m.
    str_d, &  ! The divergence stress tensor component, in Pa m.
    str_s     ! The shearing stress tensor component (cross term), in Pa m.

  ! parameters for calculating water drag and internal ice stresses
  logical :: SLAB_ICE = .false. ! should we do old style GFDL slab ice?
  real :: p0 = 2.75e4         ! pressure constant (Pa)
  real :: p0_rho              ! The pressure constant divided by ice density, N m kg-1.
  real :: c0 = 20.0           ! another pressure constant
  real :: cdw = 3.24e-3       ! ice/water drag coef. (nondim)
  real :: EC = 2.0            ! yield curve axis ratio
  real :: Rho_ocean = 1030.0  ! The nominal density of sea water, in kg m-3.
  real :: Rho_ice = 905.0     ! The nominal density of sea ice, in kg m-3.
  real :: drag_bg_vel2 = 0.0  ! A background (subgridscale) velocity for drag
                              ! with the ocean squared, in m2 s-2.
  real :: min_ocn_inertial_h = 0. ! A minimum ocean thickness used to limit the viscous coupling
                              ! rate implied for the ocean by the ice-ocean stress.
  real :: Tdamp   ! The damping timescale of the stress tensor components
                  ! toward their equilibrium solution due to the elastic terms,
                  ! in s.
  real :: del_sh_min_scale = 2.0 ! A scaling factor for the minimum permitted
                              ! value of minimum shears used in the denominator
                              ! of the stress equations, nondim.  I suspect that
                              ! this needs to be greater than 1.
  real    :: CFL_trunc        ! Velocity components will be truncated when they
                              ! are large enough that the corresponding CFL number
                              ! exceeds this value, nondim.
  logical :: CFL_check_its    ! If true, check the CFL number for every iteration
                              ! of the rheology solver; otherwise only check the
                              ! final velocities that are used for transport.
  logical :: specified_ice    ! If true, the sea ice is specified and there is
                              ! no need for ice dynamics.
  logical :: debug            ! If true, write verbose checksums for debugging purposes.
  logical :: debug_redundant  ! If true, debug redundant points.
  logical :: project_drag_vel ! If true, project forward the ice velocity used
                              ! in the drag calculation to avoid an instability
                              ! that can occur when an finite stress is applied
                              ! to thin ice moving with the velocity of the ocean.
  logical :: project_ci       ! If true, project the ice concentration and
                              ! related ice strength changes due to the convergent
                              ! or divergent ice flow.
  logical :: weak_coast_stress = .false.
  logical :: weak_low_shear = .false.
  integer :: evp_sub_steps    ! The number of iterations in the EVP dynamics
                              ! for each slow time step.
  real    :: dt_Rheo           ! The maximum sub-cycling time step for the EVP dynamics.
  type(time_type), pointer :: Time ! A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer, pointer :: ntrunc  ! The number of times the velocity has been truncated
                              ! since the last call to write_ice_statistics.
  character(len = 200) :: u_trunc_file ! The complete path to files in which a
  character(len = 200) :: v_trunc_file ! column's worth of accelerations are
                                       ! written if velocity truncations occur.
  integer :: u_file, v_file ! The unit numbers for opened u- or v- truncation
                            ! files, or -1 if they have not yet been opened.
  integer :: cols_written   ! The number of columns whose output has been
                            ! written by this PE during the current run.
  integer :: max_writes     ! The maximum number of times any PE can write out
                            ! a column's worth of accelerations during a run.

  logical :: FirstCall = .true.
  integer :: id_fix = -1, id_fiy = -1, id_fcx = -1, id_fcy = -1
  integer :: id_fwx = -1, id_fwy = -1, id_sigi = -1, id_sigii = -1
  integer :: id_stren = -1, id_stren0 = -1
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
end type SIS_C_dyn_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_C_dyn_init - initialize the ice dynamics and set parameters.             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_C_dyn_init(Time, G, param_file, diag, CS, ntrunc)
  type(time_type),     target, intent(in)    :: Time
  type(SIS_hor_grid_type),     intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(SIS_C_dyn_CS),          pointer       :: CS
  integer, target, optional,   intent(inout) :: ntrunc
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
!  (in/out,opt)  ntrunc - The integer that stores the number of times the velocity
!                     has been truncated since the last call to write_ice_statistics.

!   This subroutine sets the parameters and registers the diagnostics associated
! with the ice dynamics.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "SIS_C_dyn" ! This module's name.

  real, parameter       :: missing = -1e34

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
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "SPECIFIED_ICE", CS%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  if ( CS%specified_ice ) then
    CS%evp_sub_steps = 0 ; CS%dt_Rheo = -1.0
    call log_param(param_file, mod, "NSTEPS_DYN", CS%evp_sub_steps, &
                 "The number of iterations in the EVP dynamics for each \n"//&
                 "slow time step.  With SPECIFIED_ICE this is always 0.")
  else
    call get_param(param_file, mod, "DT_RHEOLOGY", CS%dt_Rheo, &
                 "The sub-cycling time step for iterating the rheology \n"//&
                 "and ice momentum equations. If DT_RHEOLOGY is negative, \n"//&
                 "the time step is set via NSTEPS_DYN.", units="seconds", &
                 default=-1.0)
    CS%evp_sub_steps = -1
    if (CS%dt_Rheo <= 0.0) &
      call get_param(param_file, mod, "NSTEPS_DYN", CS%evp_sub_steps, &
                 "The number of iterations of the rheology and ice \n"//&
                 "momentum equations for each slow ice time step.", default=432)
  endif
  call get_param(param_file, mod, "ICE_TDAMP_ELASTIC", CS%Tdamp, &
                 "The damping timescale associated with the elastic terms \n"//&
                 "in the sea-ice dynamics equations (if positive) or the \n"//&
                 "fraction of DT_ICE_DYNAMICS (if negative).", &
                 units = "s or nondim", default=-0.2)
  call get_param(param_file, mod, "WEAK_LOW_SHEAR_ICE", CS%weak_low_shear, &
                 "If true, the divergent stresses go toward 0 in the C-grid \n"//&
                 "dynamics when the shear magnitudes are very weak. \n"//&
                 "Otherwise they go to -P_ice.  This setting is temporary.", &
                 default=.false.)

  call get_param(param_file, mod, "PROJECT_ICE_DRAG_VEL", CS%project_drag_vel, &
                 "If true, project forward the ice velocity used in the \n"//&
                 "drag calculation to avoid an instability that can occur \n"//&
                 "when an finite stress is applied to thin ice moving with \n"//&
                 "the velocity of the ocean.", default=.true.)
  call get_param(param_file, mod, "ICE_YIELD_ELLIPTICITY", CS%EC, &
                 "The ellipticity coefficient for the plastic yield curve \n"//&
                 "in the sea-ice rheology.  For an infinite ellipticity \n"//&
                 "(i.e., a cavitating fluid rheology), use 0.", &
                 units = "Nondim", default=2.0)

  call get_param(param_file, mod, "ICE_STRENGTH_PSTAR", CS%p0, &
                 "A constant in the expression for the ice strength, \n"//&
                 "P* in Hunke & Dukowicz 1997.", units="Pa", default=2.75e4)
  call get_param(param_file, mod, "ICE_STRENGTH_CSTAR", CS%c0, &
                 "A constant in the exponent of the expression for the \n"//&
                 "ice strength, c* in Hunke & Dukowicz 1997.", &
                 units="nondim", default=20.)
  call get_param(param_file, mod, "ICE_CDRAG_WATER", CS%cdw, &
                 "The drag coefficient between the sea ice and water.", &
                 units="nondim", default=3.24e-3)
  call get_param(param_file, mod, "MIN_OCN_INTERTIAL_H", CS%min_ocn_inertial_h, &
                 "A minimum ocean thickness used to limit the viscous coupling rate\n"//&
                 "implied for the ocean by the ice-ocean stress. Only used if positive.", &
                 units="m", default=0.0)
  call get_param(param_file, mod, "ICE_DEL_SH_MIN_SCALE", CS%del_sh_min_scale, &
                 "A scaling factor for the lower bound on the shear rates \n"//&
                 "used in the denominator of the stress calculation. This \n"//&
                 "probably needs to be greater than 1.", units="nondim", default=2.0)
  call get_param(param_file, mod, "PROJECT_ICE_CONCENTRATION", CS%project_ci, &
                 "If true, project the evolution of the ice concentration \n"//&
                 "due to the convergence or divergence of the ice flow.", default=.true.)

  call get_param(param_file, mod, "RHO_OCEAN", CS%Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0)
  call get_param(param_file, mod, "RHO_ICE", CS%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0)
  CS%p0_rho = CS%p0 / CS%Rho_ice

  call get_param(param_file, mod, "CFL_TRUNCATE", CS%CFL_trunc, &
                 "The value of the CFL number that will cause ice velocity \n"//&
                 "components to be truncated; instability can occur past 0.5.", &
                 units="nondim", default=0.5)
  call get_param(param_file, mod, "CFL_TRUNC_DYN_ITS", CS%CFL_check_its, &
                 "If true, check the CFL number for every iteration of the \n"//&
                 "rheology solver; otherwise only the final velocities that \n"//&
                 "are used for transport are checked.", &
                 default=.false.)
  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call get_param(param_file, mod, "DEBUG_REDUNDANT", CS%debug_redundant, &
                 "If true, debug redundant data points.", default=CS%debug)
  if ( CS%specified_ice ) then
    CS%slab_ice = .true.
    call log_param(param_file, mod, "USE_SLAB_ICE", CS%slab_ice, &
                 "Use the very old slab-style ice.  With SPECIFIED_ICE, \n"//&
                 "USE_SLAB_ICE is always true.")
  else
    call get_param(param_file, mod, "USE_SLAB_ICE", CS%slab_ice, &
                 "If true, use the very old slab-style ice.", default=.false.)
  endif
  call get_param(param_file, mod, "U_TRUNC_FILE", CS%u_trunc_file, &
                 "The absolute path to the file where the accelerations \n"//&
                 "leading to zonal velocity truncations are written. \n"//&
                 "Leave this empty for efficiency if this diagnostic is \n"//&
                 "not needed.", default="")
  call get_param(param_file, mod, "V_TRUNC_FILE", CS%v_trunc_file, &
                 "The absolute path to the file where the accelerations \n"//&
                 "leading to meridional velocity truncations are written. \n"//&
                 "Leave this empty for efficiency if this diagnostic is \n"//&
                 "not needed.", default="")
  call get_param(param_file, mod, "MAX_TRUNC_FILE_SIZE_PER_PE", CS%max_writes, &
                 "The maximum number of colums of truncations that any PE \n"//&
                 "will write out during a run.", default=50)

!  if (len_trim(dirs%output_directory) > 0) then
!    if (len_trim(CS%u_trunc_file) > 0) &
!      CS%u_trunc_file = trim(dirs%output_directory)//trim(CS%u_trunc_file)
!    if (len_trim(CS%v_trunc_file) > 0) &
!      CS%v_trunc_file = trim(dirs%output_directory)//trim(CS%v_trunc_file)
!    call log_param(param_file, mod, "output_dir/U_TRUNC_FILE", CS%u_trunc_file)
!    call log_param(param_file, mod, "output_dir/V_TRUNC_FILE", CS%v_trunc_file)
!  endif
  CS%u_file = -1 ; CS%v_file = -1 ; CS%cols_written = 0

  CS%id_sigi  = register_diag_field('ice_model','SIGI' ,diag%axesT1, Time,     &
            'first stress invariant', 'none', missing_value=missing)
  CS%id_sigii = register_diag_field('ice_model','SIGII' ,diag%axesT1, Time,    &
            'second stress invariant', 'none', missing_value=missing)
  CS%id_stren = register_diag_field('ice_model','STRENGTH' ,diag%axesT1, Time, &
            'ice strength', 'Pa*m', missing_value=missing)
  CS%id_stren0 = register_diag_field('ice_model','STREN_0' ,diag%axesT1, Time, &
            'ice strength at start of rheology', 'Pa*m', missing_value=missing)
  CS%id_fix   = register_diag_field('ice_model', 'FI_X', diag%axesCu1, Time,   &
            'ice internal stress - x component', 'Pa', missing_value=missing,  &
            interp_method='none')
  CS%id_fiy   = register_diag_field('ice_model', 'FI_Y', diag%axesCv1, Time,   &
            'ice internal stress - y component', 'Pa', missing_value=missing,  &
            interp_method='none')
  CS%id_fcx   = register_diag_field('ice_model', 'FC_X', diag%axesCu1, Time,   &
            'Coriolis force - x component', 'Pa', missing_value=missing,       &
            interp_method='none')
  CS%id_fcy   = register_diag_field('ice_model', 'FC_Y', diag%axesCv1, Time,   &
            'Coriolis force - y component', 'Pa', missing_value=missing,       &
            interp_method='none')
  CS%id_Coru   = register_diag_field('ice_model', 'Cor_ui', diag%axesCu1, Time,&
            'Coriolis ice acceleration - x component', 'm s-2',                &
            missing_value=missing, interp_method='none')
  CS%id_Corv   = register_diag_field('ice_model', 'Cor_vi', diag%axesCv1, Time,&
            'Coriolis ice acceleration - y component', 'm s-2',                &
            missing_value=missing, interp_method='none')
  CS%id_fpx   = register_diag_field('ice_model', 'FP_X', diag%axesCu1, Time,   &
            'Pressure force - x component', 'Pa', missing_value=missing,       &
            interp_method='none')
  CS%id_fpy   = register_diag_field('ice_model', 'FP_Y', diag%axesCv1, Time,   &
            'Pressure force - y component', 'Pa', missing_value=missing,       &
            interp_method='none')
  CS%id_PFu   = register_diag_field('ice_model', 'Pfa_ui', diag%axesCu1, Time, &
            'Pressure-force ice acceleration - x component', 'm s-2',          &
            missing_value=missing, interp_method='none')
  CS%id_PFv   = register_diag_field('ice_model', 'Pfa_vi', diag%axesCv1, Time, &
            'Pressure-force ice acceleration - y component', 'm s-2',          &
            missing_value=missing,  interp_method='none')
  CS%id_fwx   = register_diag_field('ice_model', 'FW_X', diag%axesCu1, Time,   &
            'water stress on ice - x component', 'Pa', missing_value=missing,  &
            interp_method='none')
  CS%id_fwy   = register_diag_field('ice_model', 'FW_Y', diag%axesCv1, Time,   &
            'water stress on ice - y component', 'Pa', missing_value=missing,  &
            interp_method='none')
  CS%id_ui    = register_diag_field('ice_model', 'UI', diag%axesCu1, Time,     &
            'ice velocity - x component', 'm/s', missing_value=missing,        &
            interp_method='none')
  CS%id_vi    = register_diag_field('ice_model', 'VI', diag%axesCv1, Time,     &
            'ice velocity - y component', 'm/s', missing_value=missing,        &
            interp_method='none')
  CS%id_mis  = register_diag_field('ice_model', 'MIS_tot', diag%axesT1, Time,  &
            'Mass of ice and snow at t-points', 'kg m-2', missing_value=missing)
  CS%id_ci0  = register_diag_field('ice_model', 'CI_tot', diag%axesT1, Time,   &
            'Initial summed concentration of ice at t-points', 'nondim',       &
            missing_value=missing)
  CS%id_ci  = register_diag_field('ice_model', 'CI_proj', diag%axesT1, Time,   &
            'Projected summed concentration of ice at t-points', 'nondim',     &
            missing_value=missing)
  CS%id_miu   = register_diag_field('ice_model', 'MI_U', diag%axesCu1, Time,   &
            'Mass of ice and snow at u-points', 'kg m-2',                      &
            missing_value=missing, interp_method='none')
  CS%id_miv   = register_diag_field('ice_model', 'MI_V', diag%axesCv1, Time,   &
            'Mass of ice and snow at v-points', 'kg m-2',                      &
            missing_value=missing, interp_method='none')

  CS%id_fix_d   = register_diag_field('ice_model', 'FI_d_X', diag%axesCu1, Time,         &
            'ice divergence internal stress - x component', 'Pa', missing_value=missing, &
            interp_method='none')
  CS%id_fiy_d   = register_diag_field('ice_model', 'FI_d_Y', diag%axesCv1, Time,         &
            'ice divergence internal stress - y component', 'Pa', missing_value=missing, &
            interp_method='none')
  CS%id_fix_t   = register_diag_field('ice_model', 'FI_t_X', diag%axesCu1, Time,        &
            'ice tension internal stress - x component', 'Pa', missing_value=missing,   &
            interp_method='none')
  CS%id_fiy_t   = register_diag_field('ice_model', 'FI_t_Y', diag%axesCv1, Time,        &
            'ice tension internal stress - y component', 'Pa', missing_value=missing,   &
            interp_method='none')
  CS%id_fix_s   = register_diag_field('ice_model', 'FI_s_X', diag%axesCu1, Time,        &
            'ice shearing internal stress - x component', 'Pa', missing_value=missing,  &
            interp_method='none')
  CS%id_fiy_s   = register_diag_field('ice_model', 'FI_s_Y', diag%axesCv1, Time,        &
            'ice shearing internal stress - y component', 'Pa', missing_value=missing,  &
            interp_method='none')

  CS%id_str_d   = register_diag_field('ice_model', 'str_d', diag%axesT1, Time, &
            'ice divergence internal stress', 'Pa', missing_value=missing)
  CS%id_str_t   = register_diag_field('ice_model', 'str_t', diag%axesT1, Time, &
            'ice tension internal stress', 'Pa', missing_value=missing)
  CS%id_str_s   = register_diag_field('ice_model', 'str_s', diag%axesB1, Time, &
            'ice shearing internal stress', 'Pa', missing_value=missing)
  CS%id_sh_d   = register_diag_field('ice_model', 'sh_d', diag%axesT1, Time,   &
            'ice divergence strain rate', 's-1', missing_value=missing)
  CS%id_sh_t   = register_diag_field('ice_model', 'sh_t', diag%axesT1, Time,   &
            'ice tension strain rate', 's-1', missing_value=missing)
  CS%id_sh_s   = register_diag_field('ice_model', 'sh_s', diag%axesB1, Time,   &
            'ice shearing strain rate', 's-1', missing_value=missing)
  CS%id_del_sh = register_diag_field('ice_model', 'del_sh', diag%axesT1, Time, &
            'ice strain rate magnitude', 's-1', missing_value=missing)
  CS%id_del_sh_min = register_diag_field('ice_model', 'del_sh_min', diag%axesT1, Time, &
            'minimum ice strain rate magnitude', 's-1', missing_value=missing)

  CS%id_ui_hifreq = register_diag_field('ice_model', 'ui_hf', diag%axesCu1, Time, &
            'ice velocity - x component', 'm/s', missing_value=missing,        &
            interp_method='none')
  CS%id_vi_hifreq = register_diag_field('ice_model', 'vi_hf', diag%axesCv1, Time, &
            'ice velocity - y component', 'm/s', missing_value=missing,        &
            interp_method='none')
  CS%id_str_d_hifreq = register_diag_field('ice_model', 'str_d_hf', diag%axesT1, Time, &
            'ice divergence internal stress', 'Pa', missing_value=missing)
  CS%id_str_t_hifreq = register_diag_field('ice_model', 'str_t_hf', diag%axesT1, Time, &
            'ice tension internal stress', 'Pa', missing_value=missing)
  CS%id_str_s_hifreq = register_diag_field('ice_model', 'str_s_hf', diag%axesB1, Time, &
            'ice shearing internal stress', 'Pa', missing_value=missing)
  CS%id_sh_d_hifreq = register_diag_field('ice_model', 'sh_d_hf', diag%axesT1, Time, &
            'ice divergence rate', 's-1', missing_value=missing)
  CS%id_sh_t_hifreq = register_diag_field('ice_model', 'sh_t_hf', diag%axesT1, Time, &
            'ice tension rate', 's-1', missing_value=missing)
  CS%id_sh_s_hifreq = register_diag_field('ice_model', 'sh_s_hf', diag%axesB1, Time, &
            'ice shearing rate', 's-1', missing_value=missing)
  CS%id_sigi_hifreq  = register_diag_field('ice_model','sigI_hf' ,diag%axesT1, Time, &
            'first stress invariant', 'none', missing_value=missing)
  CS%id_sigii_hifreq = register_diag_field('ice_model','sigII_hf' ,diag%axesT1, Time, &
            'second stress invariant', 'none', missing_value=missing)
  CS%id_ci_hifreq  = register_diag_field('ice_model', 'CI_hf', diag%axesT1, Time, &
            'Summed concentration of ice at t-points', 'nondim', missing_value=missing)
  CS%id_stren_hifreq = register_diag_field('ice_model','STRENGTH_hf' ,diag%axesT1, Time, &
            'ice strength', 'Pa*m', missing_value=missing)

  CS%id_siu = register_diag_field('ice_model', 'siu', diag%axesT1, Time, &
            'ice velocity - x component', 'm/s', missing_value=missing,  &
            interp_method='none')
  CS%id_siv = register_diag_field('ice_model', 'siv', diag%axesT1, Time, &
            'ice velocity - y component', 'm/s', missing_value=missing,  &
            interp_method='none')
  CS%id_sispeed = register_diag_field('ice_model', 'sispeed', diag%axesT1, Time, &
            'ice speed', 'm/s', missing_value=missing)

end subroutine SIS_C_dyn_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! find_ice_strength - magnitude of force on ice in plastic deformation         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine find_ice_strength(mi, ci, ice_strength, G, CS, halo_sz) ! ??? may change to do loop
  type(SIS_hor_grid_type),          intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: mi, ci
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: ice_strength
  type(SIS_C_dyn_CS),               pointer     :: CS
  integer,                optional, intent(in)  :: halo_sz
  integer :: i, j, isc, iec, jsc, jec, halo
  halo = 0 ; if (present(halo_sz)) halo = halo_sz
  isc = G%isc-halo ; iec = G%iec+halo ; jsc = G%jsc-halo ; jec = G%jec+halo

  do j=jsc,jec ; do i=isc,iec
    ice_strength(i,j) = CS%p0_rho*mi(i,j)*exp(-CS%c0*max(1.0-ci(i,j),0.0))
  enddo ; enddo

end subroutine find_ice_strength

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_C_dynamics - take a single dynamics timestep with EVP subcycles          !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_C_dynamics(ci, msnow, mice, ui, vi, uo, vo, &
                          fxat, fyat, sea_lev, fxoc, fyoc, dt_slow, G, CS)

  type(SIS_hor_grid_type),           intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: ci, msnow, mice  ! ice properties
  real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: ui      ! ice velocity
  real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: vi      ! ice velocity
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: uo      ! ocean velocity
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: vo      ! ocean velocity
  real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: fxat  ! air stress on ice
  real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: fyat  ! air stress on ice
  real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: sea_lev     ! sea level
  real, dimension(SZIB_(G),SZJ_(G)), intent(  out) :: fxoc  ! ice stress on ocean
  real, dimension(SZI_(G),SZJB_(G)), intent(  out) :: fyoc  ! ice stress on ocean
  real,                              intent(in   ) :: dt_slow
  type(SIS_C_dyn_CS),                pointer       :: CS
! Arguments: ci - The sea ice concentration, nondim.
!  (in)      msnow - The mass per unit total area (ice covered and ice free)
!                    of the snow, in kg m-2.
!  (in)      mice - The mass per unit total area (ice covered and ice free)
!                   of the ice, in kg m-2.
!  (inout)   ui - The zonal ice velocity, in m s-1.
!  (inout)   vi - The meridional ice velocity, in m s-1.
!  (in)      uo - The zonal ocean velocity, in m s-1.
!  (in)      vo - The meridional ocean velocity, in m s-1.
!  (in)      sea_lev - The height of the sea level, including contributions
!                      from non-levitating ice from an earlier time step, in m.
!  (in)      dt_slow - The amount of time over which the ice dynamics are to be
!                      advanced, in s.
!  (in)      G - The ocean's grid structure.
!  (in/out)  CS - A pointer to the control structure for this module.

  real, dimension(SZI_(G),SZJ_(G)) :: &
    sh_Dt, &    ! sh_Dt is the horizontal tension (du/dx - dv/dy) including
                ! all metric terms, in s-1.
    sh_Dd       ! sh_Dd is the flow divergence (du/dx + dv/dy) including all
                ! metric terms, in s-1.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    sh_Ds       ! sh_Ds is the horizontal shearing strain (du/dy + dv/dx)
                ! including all metric terms, in s-1.


  real, dimension(SZI_(G),SZJ_(G)) :: &
    mis, &      ! Total snow and ice mass per unit area, in kg m-2.
    pres_mice, & ! The ice internal pressure per unit column mass, in N m / kg.
    ci_proj, &  ! The projected ice concentration, nondim.
    zeta, &     ! The ice bulk viscosity, in Pa m s or N s / m.
    del_sh, &   ! The magnitude of the shear rates, in s-1.
    diag_val, & ! A temporary diagnostic array.
    del_sh_min_pr, &  ! When multiplied by pres_mice, this gives the minimum
                ! value of del_sh that is used in the calculation of zeta,
                ! in s-1.  This is set based on considerations of numerical
                ! stability, and varies with the grid spacing.
    dx2T, dy2T, &   ! dx^2 or dy^2 at T points, in m2.
    dx_dyT, dy_dxT, &  ! dx/dy or dy_dx at T points, nondim.
    siu, siv, sispeed  ! diagnostics on T points, m/s

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    fxic, &  ! Zonal force due to internal stresses, in Pa.
    fxic_d, fxic_t, fxic_s, &
    ui_min_trunc, &  ! The range of v-velocities beyond which the velocities
    ui_max_trunc, &  ! are truncated, in m s-1, or 0 for land cells.
    Cor_u, & ! Zonal Coriolis acceleration, in m s-2.
    PFu, &   ! Zonal hydrostatic pressure driven acceleration, in m s-2.
    u_tmp, & ! A temporary copy of the old values of ui, in m s-1.
    u_IC, &  ! The initial zonal ice velocities, in m s-1.
    mi_u, &  ! The total ice and snow mass interpolated to u points, in kg m-2.
    f2dt_u, &! The squared effective Coriolis parameter at u-points times a
             ! time step, in s-1.
    I1_f2dt2_u  ! 1 / ( 1 + f^2 dt^2) at u-points, nondimensional.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    fyic, &  ! Meridional force due to internal stresses, in Pa.
    fyic_d, fyic_t, fyic_s, &
    vi_min_trunc, &  ! The range of v-velocities beyond which the velocities
    vi_max_trunc, &  ! are truncated, in m s-1, or 0 for land cells.
    Cor_v, &  ! Meridional Coriolis acceleration, in m s-2.
    PFv, &   ! Meridional hydrostatic pressure driven acceleration, in m s-2.
    v_IC, &  ! The initial meridional ice velocities, in m s-1.
    mi_v, &  ! The total ice and snow mass interpolated to v points, in kg m-2.
    f2dt_v, &! The squared effective Coriolis parameter at v-points times a
             ! time step, in s-1.
    I1_f2dt2_v  ! 1 / ( 1 + f^2 dt^2) at v-points, nondimensional.

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    mi_ratio_A_q, & ! A ratio of the masses interpolated to the faces around a
             ! vorticity point that ranges between (4 mi_min/mi_max) and 1,
             ! divided by the sum of the ocean areas around a point, in m-2.
    q, &     ! A potential-vorticity-like field for the ice, the Coriolis
             ! parameter divided by a spatially averaged mass per unit area,
             ! in s-1 m2 kg-1.
    dx2B, dy2B, &   ! dx^2 or dy^2 at B points, in m2.
    dx_dyB, dy_dxB  ! dx/dy or dy_dx at B points, nondim.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    azon, bzon, & !  _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vi & ui,
    amer, bmer, & ! respectively to get the barotropic inertial rotation,
    cmer, dmer    ! in units of s-1.  azon and amer couple the same pair of
                  ! velocities, but with the influence going in opposite
                  ! directions.

  real :: Cor       ! A Coriolis accleration, in m s-2.
  real :: fxic_now, fyic_now  ! ice internal stress convergence, in kg m-1 s-2.
  real :: drag_u, drag_v      ! Drag rates with the ocean at u & v points, in kg m-2 s-1.
  real :: drag_max  ! A maximum drag rate allowed in the ocean, in kg m-2 s-1.
  real :: tot_area  ! The sum of the area of the four neighboring cells, in m2.
  real :: dxharm    ! The harmonic mean of the x- and y- grid spacings, in m.
  real :: muq2, mvq2  ! The product of the u- and v-face masses per unit cell
                      ! area surrounding a vorticity point, in kg2 m-4.
  real :: muq, mvq    ! The u- and v-face masses per unit cell area extrapolated
                      ! to a vorticity point on the coast, in kg m-2.
  real :: pres_sum    ! The sum of the internal ice pressures aroung a point, in Pa.
  real :: min_rescale ! The smallest of the 4 surrounding values of rescale, ND.
  real :: I_1pdt_T    ! 1.0 / (1.0 + dt_2Tdamp)
  real :: I_1pE2dt_T  ! 1.0 / (1.0 + EC^2 * dt_2Tdamp)

  real :: v2_at_u     ! The squared v-velocity interpolated to u points, in m s-1.
  real :: u2_at_v     ! The squared u-velocity interpolated to v points, in m s-1.
  real :: uio_init, m_uio_explicit, uio_pred ! , uio
  real :: vio_init, m_vio_explicit, vio_pred ! , vio
  real :: I_cdRhoDt, cdRho
  real :: b_vel0      ! The initial difference between the velocity magnitude
                      ! and the absolute value of the u- or v- component, plus
                      ! the ice thickness divided by the time step and the drag
                      ! coefficient, all in m s-1.
  real :: uio_C   ! A u-velocity difference between the ocean and ice, in m s-1.
  real :: vio_C   ! A v-velocity difference between the ocean and ice, in m s-1.

  real :: Tdamp   ! The damping timescale of the stress tensor components
                  ! toward their equilibrium solution due to the elastic terms,
                  ! in s.
  real :: dt      ! The short timestep associated with the EVP dynamics, in s.
  real :: dt_2Tdamp ! The ratio of the timestep to the elastic damping timescale.
  real :: dt_cumulative ! The elapsed time within this call to EVP dynamics, in s.
  integer :: EVP_steps ! The number of EVP sub-steps that will actually be taken.
  real :: I_sub_steps  ! The number inverse of the number of EVP time steps per
                  ! slow time step.
  real :: EC2     ! EC^2, where EC is the yield curve axis ratio.
  real :: I_EC2   ! 1/EC^2, where EC is the yield curve axis ratio.
  real :: I_EC    ! 1/EC, where EC is the yield curve axis ratio.
  real :: I_2EC   ! 1/(2*EC), where EC is the yield curve axis ratio.
  real, parameter :: H_subroundoff = 1e-30 ! A negligible thickness, in m, that
                                           ! can be cubed without underflow.
  real :: m_neglect  ! A tiny mass per unit area, in kg m-2.
  real :: m_neglect2 ! A tiny mass per unit area squared, in kg2 m-4.
  real :: m_neglect4 ! A tiny mass per unit area to the 4th power, in kg4 m-8.
  real :: sum_area   ! The sum of ocean areas around a vorticity point, in m2.

  type(time_type) :: &
    time_it_start, &  ! The starting time of the iteratve steps.
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
  fxic(:,:) = 0.0 ; fyic(:,:) = 0.0
  Cor_u(:,:) = 0.0 ; Cor_v(:,:) = 0.0
  fxic_d(:,:) = 0.0 ; fyic_d(:,:) = 0.0 ; fxic_t(:,:) = 0.0 ; fyic_t(:,:) = 0.0
  fxic_s(:,:) = 0.0 ; fyic_s(:,:) = 0.0

  if (CS%SLAB_ICE) then
    ui(:,:) = uo(:,:) ; vi(:,:) = vo(:,:)
    fxoc(:,:) = fxat(:,:) ; fyoc(:,:) = fyat(:,:)
    return
  end if

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
  I_cdRhoDt = 1.0 / (CS%cdw * CS%Rho_ocean * dt)
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
      time_it_start = time_end_in - set_time(int(floor(dt_slow+0.5)))
  endif

  Tdamp = CS%Tdamp
  if (CS%Tdamp == 0.0) then
    ! Hunke (2001) chooses a specified multiple (0.36) of dt_slow for Tdamp,
    ! and shows that stability requires Tdamp > 2*dt.  Here 0.2 is used instead
    ! for greater stability.
    Tdamp = max(0.2*dt_slow, 3.0*dt)
  elseif (CS%Tdamp < 0.0) then
    Tdamp = max(-CS%Tdamp*dt_slow, 3.0*dt)
  endif
  dt_2Tdamp = dt / (2.0 * Tdamp)

  ui_min_trunc(:,:) = 0.0 ; ui_max_trunc(:,:) = 0.0
  vi_min_trunc(:,:) = 0.0 ; vi_max_trunc(:,:) = 0.0

  m_neglect = H_subroundoff*CS%Rho_ice
  m_neglect2 = m_neglect**2 ; m_neglect4 = m_neglect**4
!$OMP parallel default(none) shared(isc,iec,jsc,jec,G,CS,dt_slow,ui_min_trunc,u_IC,ui,   &
!$OMP                               ui_max_trunc,vi_min_trunc,vi_max_trunc,v_IC,vi,mice, &
!$OMP                               msnow,ci,dt,Tdamp,I_2EC,mis,ci_proj,pres_mice,       &
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
    ! Store the total snow and ice mass.
    mis(i,j) = mice(i,j) + msnow(i,j)
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
  call limit_stresses(pres_mice, mice, CS%str_d, CS%str_t, CS%str_s, G, CS)

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
    elseif (G%mask2dBu(I,J)>0.5) then
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
                 ((mis(i,j) + mis(i+1,j+1)) + (mis(i,j+1) + mis(i+1,j)))**2) * &
                 sum_area)
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
    call vec_chksum_C("PF[uv] in SIS_C_dynamics", PFu, PFv, G)
    call vec_chksum_C("f[xy]at in SIS_C_dynamics", fxat, fyat, G)
    call vec_chksum_C("[uv]i pre-steps SIS_C_dynamics", ui, vi, G)
    call vec_chksum_C("[uv]o in SIS_C_dynamics", uo, vo, G)
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
    !  The following are the forms of the horizontal tension and hori-
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

      if (max(del_sh(i,j), del_sh_min_pr(i,j)*pres_mice(i,j)) /=0.) then
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


    cdRho = CS%cdw * CS%Rho_ocean
    ! Save the current values of u for later use in updating v.
    do I=isc-1,iec
      u_tmp(I,jsc-1) = ui(I,jsc-1) ; u_tmp(I,jec+1) = ui(I,jec+1) ;
    enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,u_tmp,ui,vi,azon,bzon,czon,dzon, &
!$OMP                                  G,CS,dy2T,dx2B,vo,uo,Cor_u,f2dt_u,I1_f2dt2_u,    &
!$OMP                                  mi_u,dt,PFu,fxat,I_cdRhoDt,cdRho,m_neglect,fxoc, &
!$OMP                                  fxic,fxic_d,fxic_t,fxic_s,do_trunc_its,          &
!$OMP                                  ui_min_trunc,ui_max_trunc,drag_max) &
!$OMP                          private(Cor,fxic_now,v2_at_u,uio_init,drag_u,b_vel0, &
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

      uio_init = (ui(I,j)-uo(I,j))

      ! Determine the Coriolis acceleration and sum for averages...
      Cor_u(I,j) = Cor_u(I,j) + (Cor  - f2dt_u(I,j) * ui(I,j)) * I1_f2dt2_u(I,j)

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
            uio_pred = 0.5 * (sqrt(b_vel0**2 + 4.0*I_cdRhoDt*abs(m_uio_explicit)) - &
                              b_vel0)
          endif
          drag_u = cdRho * sqrt(max(uio_init**2, uio_pred**2) + v2_at_u )
        endif
      else
        drag_u = cdRho * sqrt(uio_init**2 + v2_at_u )
      endif
      if (drag_max>0.) drag_u = min( drag_u, drag_max )

      !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
      ! update of the other terms and an implicit bottom drag calculation.
      uio_C =  G%mask2dCu(I,j) * ( mi_u(I,j) * &
                 ((ui(I,j) + dt * Cor) * I1_f2dt2_u(I,j) - uo(I,j)) + &
                  dt * (mi_u(I,j) * PFu(I,j) + (fxic_now + fxat(I,j))) ) / &
                 (mi_u(I,j) + m_neglect + dt * drag_u)

      ui(I,j) = (uio_C + uo(I,j)) * G%mask2dCu(I,j)
      ! Note that fxoc is the stress felt by the ocean.
      fxoc(I,j) = fxoc(I,j) + drag_u*uio_C

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

      if (do_trunc_its) then
        if (ui(I,j) < ui_min_trunc(I,j)) then
          ui(I,j) = ui_min_trunc(I,j)
        elseif (ui(I,j) > ui_max_trunc(I,j)) then
          ui(I,j) = ui_max_trunc(I,j)
        endif
      endif
    enddo ; enddo
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,amer,bmer,cmer,dmer,u_tmp,G,CS, &
!$OMP                                  dx2T,dy2B,uo,vo,vi,Cor_v,f2dt_v,I1_f2dt2_v,mi_v, &
!$OMP                                  dt,PFv,fyat,I_cdRhoDt,cdRho,m_neglect,fyoc,fyic, &
!$OMP                                  fyic_d,fyic_t,fyic_s,do_trunc_its,vi_min_trunc,  &
!$OMP                                  vi_max_trunc,drag_max) &
!$OMP                          private(Cor,fyic_now,u2_at_v,vio_init,drag_v,    &
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
          b_vel0 = mi_v(i,J) * I_cdRhoDt + &
                   (sqrt(vio_init**2 + u2_at_v) - abs(vio_init))
          if (b_vel0**2 > 1e8*I_cdRhoDt*abs(m_vio_explicit)) then
            vio_pred = m_vio_explicit * I_cdRhoDt / b_vel0
          else
            vio_pred = 0.5 * (sqrt(b_vel0**2 + 4.0*I_cdRhoDt*abs(m_vio_explicit)) - &
                              b_vel0)
          endif
          drag_v = cdRho * sqrt(max(vio_init**2, vio_pred**2) + u2_at_v )
        endif
      else
        drag_v = cdRho * sqrt(vio_init**2 + u2_at_v )
      endif
      if (drag_max>0.) drag_v = min( drag_v, drag_max )

      !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
      ! update of the other terms and an implicit bottom drag calculation.
      vio_C =  G%mask2dCv(i,J) * ( mi_v(i,J) * &
                 ((vi(i,J) + dt * Cor) * I1_f2dt2_v(i,J) - vo(i,J)) + &
                  dt * (mi_v(i,J) * PFv(i,J) + (fyic_now + fyat(i,J))) ) / &
                 (mi_v(i,J) + m_neglect + dt * drag_v)

      vi(i,J) = (vio_C + vo(i,J)) * G%mask2dCv(i,J)
      ! Note that fyoc is the stress felt by the ocean.
      fyoc(i,J) = fyoc(i,J) + drag_v*vio_C

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

      if (do_trunc_its) then
        if (vi(i,J) < vi_min_trunc(i,J)) then
          vi(i,J) = vi_min_trunc(i,J)
        elseif (vi(i,J) > vi_max_trunc(i,J)) then
          vi(i,J) = vi_max_trunc(i,J)
        endif
      endif
    enddo ; enddo

    if (do_hifreq_output) then
      time_step_end = time_it_start + set_time(int(floor(n*dt+0.5)))
      call enable_SIS_averaging(dt, time_step_end, CS%diag)
      if (CS%id_ui_hifreq > 0) call post_SIS_data(CS%id_ui_hifreq, ui, CS%diag)
      if (CS%id_vi_hifreq > 0) call post_SIS_data(CS%id_vi_hifreq, vi, CS%diag)
      if (CS%id_str_d_hifreq > 0) call post_SIS_data(CS%id_str_d_hifreq, CS%str_d, CS%diag)
      if (CS%id_str_t_hifreq > 0) call post_SIS_data(CS%id_str_t_hifreq, CS%str_t, CS%diag)
      if (CS%id_str_s_hifreq > 0) call post_SIS_data(CS%id_str_s_hifreq, CS%str_s, CS%diag)
      if (CS%id_sh_d_hifreq > 0) call post_SIS_data(CS%id_sh_d_hifreq, sh_Dd, CS%diag)
      if (CS%id_sh_t_hifreq > 0) call post_SIS_data(CS%id_sh_t_hifreq, sh_Dt, CS%diag)
      if (CS%id_sh_s_hifreq > 0) call post_SIS_data(CS%id_sh_s_hifreq, sh_Ds, CS%diag)
      if (CS%id_sigi_hifreq>0) then
        call find_sigI(mice, ci_proj, CS%str_d, diag_val, G, CS)
        call post_SIS_data(CS%id_sigi_hifreq, diag_val, CS%diag)
      endif
      if (CS%id_sigii_hifreq>0) then
        call find_sigII(mice, ci_proj, CS%str_t, CS%str_s, diag_val, G, CS)
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

    if (CS%debug) then
      call hchksum(CS%str_d, "str_d in SIS_C_dynamics", G%HI, haloshift=1)
      call hchksum(CS%str_t, "str_t in SIS_C_dynamics", G%HI, haloshift=1)
      call Bchksum(CS%str_s, "str_s in SIS_C_dynamics", G%HI, &
                   haloshift=0, symmetric=.true.)
    endif
    if (CS%debug .or. CS%debug_redundant) then
      call vec_chksum_C("f[xy]ic in SIS_C_dynamics", fxic, fyic, G)
      call vec_chksum_C("f[xy]oc in SIS_C_dynamics", fxoc, fyoc, G)
      call vec_chksum_C("Cor_[uv] in SIS_C_dynamics", Cor_u, Cor_v, G)
      call vec_chksum_C("[uv]i in SIS_C_dynamics", ui, vi, G)
    endif

  enddo ! l=1,EVP_steps

  if (CS%debug .or. CS%debug_redundant) &
    call vec_chksum_C("[uv]i end SIS_C_dynamics", ui, vi, G)

  ! Reset the time information in the diag type.
  if (do_hifreq_output) call enable_SIS_averaging(time_int_in, time_end_in, CS%diag)

  ! make averages
  I_sub_steps = 1.0/EVP_steps
!$OMP parallel default(none) shared(isc,iec,jsc,jec,G,fxoc,fxic,Cor_u,fxic_d,fxic_t, &
!$OMP                               fxic_s,I_sub_steps,fyoc,fyic,Cor_v,fyic_d,       &
!$OMP                               fyic_t,fyic_s)
!$OMP do
  do j=jsc,jec ; do I=isc-1,iec
    fxoc(I,j) = fxoc(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
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
                               PFu, fxat, dt_slow, G, CS)
          endif
          if (ui(I,j) < ui_min_trunc(I,j)) then
            ui(I,j) = 0.95 * ui_min_trunc(I,j)
          else
            ui(I,j) = 0.95*ui_max_trunc(I,j)
          endif
        endif
      enddo ; enddo
    else
      do j=jsc,jec ; do I=isc-1,iec
        if (ui(I,j) < ui_min_trunc(I,j)) then
          ui(I,j) = 0.95 * ui_min_trunc(I,j)
          if (mi_u(I,j) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        elseif (ui(I,j) > ui_max_trunc(I,j)) then
          ui(I,j) = 0.95*ui_max_trunc(I,j)
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
                               PFv, fyat, dt_slow, G, CS)
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
    if (CS%id_fcx>0) call post_SIS_data(CS%id_fcx, Cor_u(:,:)*mi_u(:,:), CS%diag)
    if (CS%id_fcy>0) call post_SIS_data(CS%id_fcy, Cor_v(:,:)*mi_v(:,:), CS%diag)
    if (CS%id_Coru>0) call post_SIS_data(CS%id_fcx, Cor_u, CS%diag)
    if (CS%id_Corv>0) call post_SIS_data(CS%id_fcy, Cor_v, CS%diag)
    if (CS%id_PFu>0) call post_SIS_data(CS%id_PFu, PFu, CS%diag)
    if (CS%id_PFv>0) call post_SIS_data(CS%id_PFv, PFv, CS%diag)
    if (CS%id_fpx>0) call post_SIS_data(CS%id_fpx, PFu(:,:)*mi_u(:,:), CS%diag)
    if (CS%id_fpy>0) call post_SIS_data(CS%id_fpy, PFv(:,:)*mi_v(:,:), CS%diag)
    if (CS%id_fwx>0) call post_SIS_data(CS%id_fwx, -fxoc, CS%diag) ! water force on ice
    if (CS%id_fwy>0) call post_SIS_data(CS%id_fwy, -fyoc, CS%diag) ! ...= -ice on water
!  The diagnostics of fxat and fyat are supposed to be taken over all partitions
!  (ocean & ice), whereas fxat and fyat here are only averaged over the ice.

    if (CS%id_fix_d>0) call post_SIS_data(CS%id_fix_d, fxic_d, CS%diag)
    if (CS%id_fiy_d>0) call post_SIS_data(CS%id_fiy_d, fyic_d, CS%diag)
    if (CS%id_fix_t>0) call post_SIS_data(CS%id_fix_t, fxic_t, CS%diag)
    if (CS%id_fiy_t>0) call post_SIS_data(CS%id_fiy_t, fyic_t, CS%diag)
    if (CS%id_fix_s>0) call post_SIS_data(CS%id_fix_s, fxic_s, CS%diag)
    if (CS%id_fiy_s>0) call post_SIS_data(CS%id_fiy_s, fyic_s, CS%diag)

    if (CS%id_sigi>0) then
      call find_sigI(mice, ci, CS%str_d, diag_val, G, CS)
      call post_SIS_data(CS%id_sigi, diag_val, CS%diag)
    endif
    if (CS%id_sigii>0) then
      call find_sigII(mice, ci, CS%str_t, CS%str_s, diag_val, G, CS)
      call post_SIS_data(CS%id_sigii, diag_val, CS%diag)
    endif
    if (CS%id_stren>0) then
      if (CS%project_ci) then
        call find_ice_strength(mice, ci_proj, diag_val, G, CS)
      else
        call find_ice_strength(mice, ci, diag_val, G, CS)
      endif
      call post_SIS_data(CS%id_stren, diag_val, CS%diag)
    endif
    if (CS%id_stren0>0) then
      call find_ice_strength(mice, ci, diag_val, G, CS)
      call post_SIS_data(CS%id_stren0, diag_val, CS%diag)
    endif

    if (CS%id_ui>0) call post_SIS_data(CS%id_ui, ui, CS%diag)
    if (CS%id_vi>0) call post_SIS_data(CS%id_vi, vi, CS%diag)
    if (CS%id_miu>0) call post_SIS_data(CS%id_miu, mi_u, CS%diag)
    if (CS%id_miv>0) call post_SIS_data(CS%id_miv, mi_v, CS%diag)
    if (CS%id_mis>0) call post_SIS_data(CS%id_mis, mice, CS%diag)
    if (CS%id_ci0>0) call post_SIS_data(CS%id_ci0, ci, CS%diag)
    if (CS%id_ci>0)  call post_SIS_data(CS%id_ci, ci_proj, CS%diag)

    if (CS%id_str_d>0) call post_SIS_data(CS%id_str_d, CS%str_d, CS%diag)
    if (CS%id_str_t>0) call post_SIS_data(CS%id_str_t, CS%str_t, CS%diag)
    if (CS%id_str_s>0) call post_SIS_data(CS%id_str_s, CS%str_s, CS%diag)

    if (CS%id_sh_d>0) call post_SIS_data(CS%id_sh_d, sh_Dd, CS%diag)
    if (CS%id_sh_t>0) call post_SIS_data(CS%id_sh_t, sh_Dt, CS%diag)
    if (CS%id_sh_s>0) call post_SIS_data(CS%id_sh_s, sh_Ds, CS%diag)

    if (CS%id_del_sh>0) call post_SIS_data(CS%id_del_sh, del_sh, CS%diag)
    if (CS%id_del_sh_min>0) call post_SIS_data(CS%id_del_sh_min, &
                                    (del_sh_min_pr(:,:)*pres_mice(:,:)), CS%diag)
    if (Cs%id_siu>0 .or. Cs%id_siv>0 .or. Cs%id_sispeed>0) then

      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        if (mis(i,j) > 0.0) then
          siu(i,j) = (ui(I-1,j) + ui(I,j))/2
          siv(i,j) = (vi(i,J-1) + vi(i,J))/2
          sispeed(i,j) = (siu(i,j)*siu(i,j)+siv(i,j)*siv(i,j))**0.5
        else
          siu(i,j) = 0.0; siv(i,j) = 0.0; sispeed(i,j) = 0.0;
        endif
      enddo ; enddo
      if (Cs%id_siu>0) call post_SIS_data(CS%id_siu, siu, CS%diag)
      if (Cs%id_siv>0) call post_SIS_data(CS%id_siv, siv, CS%diag)
      if (Cs%id_sispeed>0) call post_SIS_data(CS%id_sispeed, sispeed, CS%diag)
    endif

  endif

end subroutine SIS_C_dynamics

subroutine limit_stresses(pres_mice, mice, str_d, str_t, str_s, G, CS, limit)
  type(SIS_hor_grid_type),            intent(in)    :: G
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)    :: pres_mice, mice
  real, dimension(SZI_(G),SZJ_(G)),   intent(inout) :: str_d, str_t
  real, dimension(SZIB_(G),SZJB_(G)), intent(inout) :: str_s
  type(SIS_C_dyn_CS),                 pointer       :: CS
  real, optional,                     intent(in)    :: limit
! Arguments: pres_mice - The ice internal pressure per unit column mass, in N m / kg.
!  (in)      mice - The mass per unit total area (ice covered and ice free)
!                   of the ice, in kg m-2.
!  (in/out)  str_t - The tension stress tensor component, in Pa m.
!  (in/out)  str_d - The divergence stress tensor component, in Pa m.
!  (in/out)  str_s - The shearing stress tensor component (cross term), in Pa m.
!  (in)      G - The ocean's grid structure.
!  (in)      CS - A pointer to the control structure for this module.
!  (in,opt)  limit - a factor by which the strength limits are changed.

!   This subroutine ensures that the input stresses are not larger than could
! be justified by the ice pressure now, as the ice might have melted or been
! advected away during the thermodynamic and transport phases, or the
! ice flow convergence or divergence may have altered the ice concentration.

  real :: pressure  ! The internal ice pressure at a point, in Pa.
  real :: pres_avg  ! The average of the internal ice pressures around a point, in Pa.
  real :: sum_area  ! The sum of ocean areas around a vorticity point, in m2.
  real :: I_2EC     ! 1/(2*EC), where EC is the yield curve axis ratio.
  real :: lim       ! A local copy of the factor by which the limits are changed.
  real :: lim_2     ! The limit divided by 2.
!  real :: EC2       ! EC^2, where EC is the yield curve axis ratio.
!  real :: rescale_str ! A factor by which to rescale the internal stresses, ND.
!  real :: stress_mag  ! The magnitude of the stress at a point.
!  real :: str_d_q     ! CS%str_d interpolated to a vorticity point, in Pa m.
!  real :: str_t_q     ! CS%str_t interpolated to a vorticity point, in Pa m.

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
!        if ((stress_mag > pres_avg) .and. (G%mask2dBu(I,j)>0.5)) &
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
! sigI - first stress invariant                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine find_sigI(mi, ci, str_d, sigI, G, CS)
  type(SIS_hor_grid_type),          intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: mi, ci, str_d
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: sigI
  type(SIS_C_dyn_CS),               pointer     :: CS

  real, dimension(SZI_(G),SZJ_(G)) :: &
    strength ! The ice strength, in Pa.
  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call find_ice_strength(mi, ci, strength, G, CS)

  do j=jsc,jec ; do i=isc,iec
    sigI(i,j) = 0.0
    if (strength(i,j) > 0.0) sigI(i,j) = 2.0 * str_d(i,j) / strength(i,j)
  enddo ; enddo

end subroutine find_sigI

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! sigII - second stress invariant                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine find_sigII(mi, ci, str_t, str_s, sigII, G, CS)
  type(SIS_hor_grid_type),            intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: mi, ci, str_t
  real, dimension(SZIB_(G),SZJB_(G)), intent(in)  :: str_s
  real, dimension(SZI_(G),SZJ_(G)),   intent(out) :: sigII
  type(SIS_C_dyn_CS),                 pointer     :: CS

  real, dimension(SZI_(G),SZJ_(G)) :: &
    strength ! The ice strength, in Pa.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    str_s_ss ! Str_s divided by the sum of the neighboring ice strengths.
  real :: strength_sum  ! The sum of the 4 neighboring strengths, in Pa.
  real :: sum_area   ! The sum of ocean areas around a vorticity point, in m2.
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_C_dyn_register_restarts - allocate and register any variables for this   !
!      module that need to be included in the restart files.                   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_C_dyn_register_restarts(mpp_domain, HI, param_file, CS, &
                                       Ice_restart, restart_file)
  type(domain2d),          intent(in) :: mpp_domain
  type(hor_index_type),    intent(in) :: HI
  type(param_file_type),   intent(in) :: param_file
  type(SIS_C_dyn_CS),      pointer    :: CS
  type(restart_file_type), pointer    :: Ice_restart
  character(len=*),        intent(in) :: restart_file

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
!

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

  allocate(CS%str_d(isd:ied, jsd:jed)) ; CS%str_d(:,:) = 0.0
  allocate(CS%str_t(isd:ied, jsd:jed)) ; CS%str_t(:,:) = 0.0
  allocate(CS%str_s(HI%IsdB:HI%IedB, HI%JsdB:HI%JedB)) ; CS%str_s(:,:) = 0.0
  if (associated(Ice_restart)) then
    id = register_restart_field(Ice_restart, restart_file, 'str_d', CS%str_d, &
                                domain=mpp_domain, mandatory=.false.)
    id = register_restart_field(Ice_restart, restart_file, 'str_t', CS%str_t, &
                                domain=mpp_domain, mandatory=.false.)
    id = register_restart_field(Ice_restart, restart_file, 'str_s', CS%str_s, &
                     domain=mpp_domain, position=CORNER, mandatory=.false.)
   endif
end subroutine SIS_C_dyn_register_restarts


! The following two subroutines are used to record the location of velocity
! truncations and related fields.

subroutine write_u_trunc(I, j, ui, u_IC, uo, mis, fxoc, fxic, Cor_u, PFu, fxat, &
                         dt_slow, G, CS)
  integer, intent(in) :: I, j
  type(SIS_hor_grid_type),           intent(in) :: G
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: ui, u_IC, uo
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: mis
  real, dimension(SZIB_(G),SZJ_(G)), intent(in) :: fxoc, fxic, Cor_u, PFu, fxat
  real,                              intent(in) :: dt_slow
  type(SIS_C_dyn_CS),                pointer    :: CS

  real :: dt_mi, CFL
  real, parameter :: H_subroundoff = 1e-30 ! A negligible thickness, in m, that
                                           ! can be cubed without underflow.
  integer :: file
  integer :: yr, mo, day, hr, minute, sec, yearday

  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

  ! Open up the file for output if this is the first call.
    if (CS%u_file < 0) then
      if (len_trim(CS%u_trunc_file) < 1) return
      call open_file(CS%u_file, trim(CS%u_trunc_file), action=APPEND_FILE, &
                     form=ASCII_FILE, threading=MULTIPLE, fileset=SINGLE_FILE)
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


    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'("Time ",i5,i4,F6.2," U-trunc at ",I4,": ",2(I3), &
        & " (",F7.2," E "F7.2," N) u = ",ES10.3," (CFL ",ES9.2,") was ",ES10.3," dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), I, j, &
        G%geoLonCu(I,j), G%geoLatCu(I,j), ui(I,j), CFL, u_IC(I,j), dt_slow

    dt_mi = dt_slow / (0.5*(mis(i,j) + mis(i+1,j)) + H_subroundoff*CS%Rho_ice)

    write (file, '("ui, uo, dui = ", 3ES11.3, " ;  mice+snow = ",2ES11.3)') &
      ui(I,j), uo(I,j), ui(I,j) - u_IC(I,j), mis(i,j), mis(i+1,j)

    write (file, '("U change due to fxat, fxoc, fxic, Cor_u, PFu = ", 5ES11.3, " sum = ",ES11.3)') &
      fxat(I,j)*dt_mi, -fxoc(I,j)*dt_mi, fxic(I,j)*dt_mi, Cor_u(I,j)*dt_slow, PFu(I,j)*dt_slow, &
      (fxat(I,j) - fxoc(I,j) + fxic(I,j))*dt_mi + (Cor_u(I,j) + PFu(I,j))*dt_slow

    call flush(file)
  endif

end subroutine write_u_trunc

subroutine write_v_trunc(i, J, vi, v_IC, vo, mis, fyoc, fyic, Cor_v, PFv, fyat, &
                         dt_slow, G, CS)
  integer, intent(in) :: i, j
  type(SIS_hor_grid_type),           intent(in) :: G
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: vi, v_IC, vo
  real, dimension(SZI_(G),SZJ_(G)),  intent(in) :: mis
  real, dimension(SZI_(G),SZJB_(G)), intent(in) :: fyoc, fyic, Cor_v, PFv, fyat
  real,                              intent(in) :: dt_slow
  type(SIS_C_dyn_CS),                pointer    :: CS

  real :: dt_mi, CFL
  real, parameter :: H_subroundoff = 1e-30 ! A negligible thickness, in m, that
                                           ! can be cubed without underflow.
  integer :: file
  integer :: yr, mo, day, hr, minute, sec, yearday


  if (CS%cols_written < CS%max_writes) then
    CS%cols_written = CS%cols_written + 1

  ! Open up the file for output if this is the first call.
    if (CS%v_file < 0) then
      if (len_trim(CS%v_trunc_file) < 1) return
      call open_file(CS%v_file, trim(CS%v_trunc_file), action=APPEND_FILE, &
                     form=ASCII_FILE, threading=MULTIPLE, fileset=SINGLE_FILE)
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

    call get_date(CS%Time, yr, mo, day, hr, minute, sec)
    call get_time((CS%Time - set_date(yr, 1, 1, 0, 0, 0)), sec, yearday)
    write (file,'(/,"--------------------------")')
    write (file,'("Time ",i5,i4,F6.2," V-trunc at ",I4,": ",2(I3), &
        & " (",F7.2," E ",F7.2," N) v = ",ES10.3," (CFL ",ES9.2,") was ",ES10.3," dt = ",1PG10.4)') &
        yr, yearday, (REAL(sec)/3600.0), pe_here(), i, J, &
        G%geoLonCv(i,J), G%geoLatCv(i,J), vi(i,J), CFL, v_IC(i,J), dt_slow

    dt_mi = dt_slow / (0.5*(mis(i,j) + mis(i,j+1)) + H_subroundoff*CS%Rho_ice)

    write (file, '("vi, vo, dvi = ", 3ES11.3, " ;  mice+snow = ",2ES11.3)') &
      vi(i,J), vo(i,J), vi(i,J) - v_IC(i,J), mis(i,j), mis(i,j+1)

    write (file, '("V change due to fyat, fyoc, fyic, Cor_v, PFv = ", 5ES11.3, " sum = ",ES11.3)') &
      fyat(i,J)*dt_mi, -fyoc(i,J)*dt_mi, fyic(i,J)*dt_mi, Cor_v(i,J)*dt_slow, PFv(i,J)*dt_slow, &
      (fyat(i,J) - fyoc(i,J) + fyic(i,J))*dt_mi + (Cor_v(i,J) + PFv(i,J))*dt_slow

    call flush(file)
  endif

end subroutine write_v_trunc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! SIS_C_dyn_end - deallocate the memory associated with this module.           !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine SIS_C_dyn_end(CS)
  type(SIS_C_dyn_CS), pointer :: CS

  deallocate(CS%str_d) ; deallocate(CS%str_t) ; deallocate(CS%str_s)

  deallocate(CS)
end subroutine SIS_C_dyn_end

end module SIS_dyn_cgrid
