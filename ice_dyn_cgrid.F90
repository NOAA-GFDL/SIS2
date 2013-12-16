module ice_dyn_cgrid
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
! Robert Hallberg in 2013   .                                                  !
!                                                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use SIS_diag_mediator, only : post_SIS_data, SIS_diag_ctrl
use SIS_diag_mediator, only : query_SIS_averaging_enabled, enable_SIS_averaging
use SIS_diag_mediator, only : register_diag_field=>register_SIS_diag_field
use SIS_error_checking, only : chksum, Bchksum, uchksum, vchksum, hchksum
use SIS_error_checking, only : check_redundant_B, check_redundant_C
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_file_parser,  only : get_param, log_param, read_param, log_version, param_file_type
use MOM_domains,      only : pass_var, pass_vector, BGRID_NE, CGRID_NE, CORNER
use ice_grid_mod,     only : sea_ice_grid_type
use fms_io_mod,       only : register_restart_field, restart_file_type
use time_manager_mod, only : time_type, set_time, operator(+), operator(-)

implicit none ; private

#include <SIS2_memory.h>

public :: ice_C_dyn_init, ice_C_dynamics, ice_C_dyn_end, ice_C_dyn_register_restarts

type, public :: ice_C_dyn_CS ; private
  real, allocatable, dimension(:,:) :: &
    str_t, &  ! The tension stress tensor component, in Pa m.
    str_d, &  ! The divergence stress tensor component, in Pa m.
    str_s     ! The shearing stress tensor component (cross term), in Pa m.

  ! parameters for calculating water drag and internal ice stresses
  logical :: SLAB_ICE = .false. ! should we do old style GFDL slab ice?
  real :: p0 = 2.75e4         ! pressure constant (Pa)
  real :: c0 = 20.0           ! another pressure constant
  real :: cdw = 3.24e-3       ! ice/water drag coef. (nondim)
  real :: EC = 2.0            ! yield curve axis ratio
  real :: Rho_ocean = 1030.0  ! The nominal density of sea water, in kg m-3.
  real :: Rho_ice = 905.0     ! The nominal density of sea ice, in kg m-3.
  real :: Rho_snow = 330.0    ! The nominal density of snow on sea ice, in kg m-3.
  real :: drag_bg_vel2 = 0.0  ! A background (subgridscale) velocity for drag
                              ! with the ocean squared, in m2 s-2.
  real :: Tdamp   ! The damping timescale of the stress tensor components
                  ! toward their equilibrium solution due to the elastic terms,
                  ! in s.
  real :: del_sh_min_scale = 2.0 ! A scaling factor for the minimum permitted
                              ! value of minimum shears used in the denominator
                              ! of the stress equations, nondim.  I suspect that
                              ! this needs to be greater than 1.
  logical :: specified_ice    ! If true, the sea ice is specified and there is
                              ! no need for ice dynamics.
  logical :: debug            ! If true, write verbose checksums for debugging purposes.
  logical :: debug_redundant  ! If true, debug redundant points.
  logical :: weak_coast_stress = .false.
  integer :: evp_sub_steps    ! The number of iterations in the EVP dynamics
                              ! for each slow time step.
  type(time_type), pointer :: Time ! A pointer to the ice model's clock.
  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_fix = -1, id_fiy = -1, id_fcx = -1, id_fcy = -1
  integer :: id_fwx = -1, id_fwy = -1, id_sigi = -1, id_sigii = -1
  integer :: id_stren = -1, id_ui = -1, id_vi = -1, id_Coru = -1, id_Corv = -1
  integer :: id_PFu = -1, id_PFv = -1, id_fpx = -1, id_fpy = -1
  integer :: id_fix_d = -1, id_fix_t = -1, id_fix_s = -1
  integer :: id_fiy_d = -1, id_fiy_t = -1, id_fiy_s = -1
  integer :: id_str_d = -1, id_str_t = -1, id_str_s = -1
  integer :: id_sh_d = -1, id_sh_t = -1, id_sh_s = -1
  integer :: id_del_sh = -1, id_del_sh_min = -1
  integer :: id_ui_hifreq = -1, id_vi_hifreq = -1
  integer :: id_str_d_hifreq = -1, id_str_t_hifreq = -1, id_str_s_hifreq = -1
  integer :: id_sh_d_hifreq = -1, id_sh_t_hifreq = -1, id_sh_s_hifreq = -1
  integer :: id_sigi_hifreq = -1, id_sigii_hifreq = -1
end type ice_C_dyn_CS

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_C_dyn_init - initialize the ice dynamics and set parameters.             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_C_dyn_init(Time, G, param_file, diag, CS)
  type(time_type),     target, intent(in)    :: Time
  type(sea_ice_grid_type),     intent(in)    :: G
  type(param_file_type),       intent(in)    :: param_file
  type(SIS_diag_ctrl), target, intent(inout) :: diag
  type(ice_C_dyn_CS),          pointer       :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.

!   This subroutine sets the parameters and registers the diagnostics associated
! with the ice dynamics.

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "ice_C_dyn" ! This module's name.

  real, parameter       :: missing = -1e34

  if (.not.associated(CS)) then
    call SIS_error(FATAL, "ice_dyn_init called with an unassociated control structure. \n"//&
                    "ice_dyn_register_restarts must be called before ice_dyn_init.")
    return
  endif

  CS%diag => diag
  CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version)
  call get_param(param_file, mod, "SPECIFIED_ICE", CS%specified_ice, &
                 "If true, the ice is specified and there is no dynamics.", &
                 default=.false.)
  if ( CS%specified_ice ) then
    CS%evp_sub_steps = 0
    call log_param(param_file, mod, "NSTEPS_DYN", CS%evp_sub_steps, &
                 "The number of iterations in the EVP dynamics for each \n"//&
                 "slow time step.  With SPECIFIED_ICE this is always 0.")
  else
    call get_param(param_file, mod, "NSTEPS_DYN", CS%evp_sub_steps, &
                 "The number of iterations in the EVP dynamics for each \n"//&
                 "slow time step.", default=432)
  endif
  call get_param(param_file, mod, "ICE_TDAMP_ELASTIC", CS%Tdamp, &
                 "The damping timescale associated with the elastic terms \n"//&
                 "in the sea-ice dynamics equations.", units = "s", default=0.0)

  call get_param(param_file, mod, "ICE_YEILD_ELLIPTICITY", CS%EC, &
                 "The ellipticity coefficient for the plastic yeild curve \n"//&
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
  call get_param(param_file, mod, "ICE_DEL_SH_MIN_SCALE", CS%del_sh_min_scale, &
                 "A scaling factor for the lower bound on the shear rates \n"//&
                 "used in the denominator of the stress calculation. This \n"//&
                 "probably needs to be greater than 1.", units="nondim", default=2.0)

  call get_param(param_file, mod, "RHO_OCEAN", CS%Rho_ocean, &
                 "The nominal density of sea water as used by SIS.", &
                 units="kg m-3", default=1030.0)
  call get_param(param_file, mod, "RHO_ICE", CS%Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0)
  call get_param(param_file, mod, "RHO_SNOW", CS%Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0)

  call get_param(param_file, mod, "DEBUG", CS%debug, &
                 "If true, write out verbose debugging data.", default=.false.)
  call get_param(param_file, mod, "DEBUG_REDUNDANT", CS%debug_redundant, &
                 "If true, debug redundant data points.", default=CS%debug)
  call get_param(param_file, mod, "USE_SLAB_ICE", CS%SLAB_ICE, &
                 "If true, use the very old slab-style ice.", default=.false.)

  CS%id_sigi  = register_diag_field('ice_model','SIGI' ,diag%axesT1, Time,         &
            'first stress invariant', 'none', missing_value=missing)
  CS%id_sigii = register_diag_field('ice_model','SIGII' ,diag%axesT1, Time,        &
            'second stress invariant', 'none', missing_value=missing)
  CS%id_stren = register_diag_field('ice_model','STRENGTH' ,diag%axesT1, Time,     &
            'ice strength', 'Pa*m', missing_value=missing)
  CS%id_fix   = register_diag_field('ice_model', 'FI_X', diag%axesCu1, Time,        &
            'ice internal stress - x component', 'Pa', missing_value=missing)
  CS%id_fiy   = register_diag_field('ice_model', 'FI_Y', diag%axesCv1, Time,        &
            'ice internal stress - y component', 'Pa', missing_value=missing)
  CS%id_fcx   = register_diag_field('ice_model', 'FC_X', diag%axesCu1, Time,        &
            'Coriolis force - x component', 'Pa', missing_value=missing)
  CS%id_fcy   = register_diag_field('ice_model', 'FC_Y', diag%axesCv1, Time,        &
            'Coriolis force - y component', 'Pa', missing_value=missing)
  CS%id_Coru   = register_diag_field('ice_model', 'Cor_ui', diag%axesCu1, Time,     &
            'Coriolis ice acceleration - x component', 'm s-2', missing_value=missing)
  CS%id_Corv   = register_diag_field('ice_model', 'Cor_vi', diag%axesCv1, Time,     &
            'Coriolis ice acceleration - y component', 'm s-2', missing_value=missing)
  CS%id_fpx   = register_diag_field('ice_model', 'FP_X', diag%axesCu1, Time,        &
            'Pressure force - x component', 'Pa', missing_value=missing)
  CS%id_fpy   = register_diag_field('ice_model', 'FP_Y', diag%axesCv1, Time,        &
            'Pressure force - y component', 'Pa', missing_value=missing)
  CS%id_PFu   = register_diag_field('ice_model', 'Pfa_ui', diag%axesCu1, Time,     &
            'Pressure-force ice acceleration - x component', 'm s-2', missing_value=missing)
  CS%id_PFv   = register_diag_field('ice_model', 'Pfa_vi', diag%axesCv1, Time,     &
            'Pressure-force ice acceleration - y component', 'm s-2', missing_value=missing)
  CS%id_fwx   = register_diag_field('ice_model', 'FW_X', diag%axesCu1, Time,        &
            'water stress on ice - x component', 'Pa', missing_value=missing)
  CS%id_fwy   = register_diag_field('ice_model', 'FW_Y', diag%axesCv1, Time,        &
            'water stress on ice - y component', 'Pa', missing_value=missing)
  CS%id_ui    = register_diag_field('ice_model', 'UI', diag%axesCu1, Time,          &
            'ice velocity - x component', 'm/s', missing_value=missing)
  CS%id_vi    = register_diag_field('ice_model', 'VI', diag%axesCv1, Time,          &
            'ice velocity - y component', 'm/s', missing_value=missing)

  CS%id_fix_d   = register_diag_field('ice_model', 'FI_d_X', diag%axesCu1, Time,        &
            'ice divergence internal stress - x component', 'Pa', missing_value=missing)
  CS%id_fiy_d   = register_diag_field('ice_model', 'FI_d_Y', diag%axesCv1, Time,        &
            'ice divergence internal stress - y component', 'Pa', missing_value=missing)
  CS%id_fix_t   = register_diag_field('ice_model', 'FI_t_X', diag%axesCu1, Time,        &
            'ice tension internal stress - x component', 'Pa', missing_value=missing)
  CS%id_fiy_t   = register_diag_field('ice_model', 'FI_t_Y', diag%axesCv1, Time,        &
            'ice tension internal stress - y component', 'Pa', missing_value=missing)
  CS%id_fix_s   = register_diag_field('ice_model', 'FI_s_X', diag%axesCu1, Time,        &
            'ice shearing internal stress - x component', 'Pa', missing_value=missing)
  CS%id_fiy_s   = register_diag_field('ice_model', 'FI_s_Y', diag%axesCv1, Time,        &
            'ice shearing internal stress - y component', 'Pa', missing_value=missing)

  CS%id_str_d   = register_diag_field('ice_model', 'str_d', diag%axesT1, Time, &
            'ice divergence internal stress', 'Pa', missing_value=missing)
  CS%id_str_t   = register_diag_field('ice_model', 'str_t', diag%axesT1, Time, &
            'ice tension internal stress', 'Pa', missing_value=missing)
  CS%id_str_s   = register_diag_field('ice_model', 'str_s', diag%axesB1, Time, &
            'ice shearing internal stress', 'Pa', missing_value=missing)
  CS%id_sh_d   = register_diag_field('ice_model', 'sh_d', diag%axesT1, Time, &
            'ice divergence strain rate', 's-1', missing_value=missing)
  CS%id_sh_t   = register_diag_field('ice_model', 'sh_t', diag%axesT1, Time, &
            'ice tension strain rate', 's-1', missing_value=missing)
  CS%id_sh_s   = register_diag_field('ice_model', 'sh_s', diag%axesB1, Time, &
            'ice shearing strain rate', 's-1', missing_value=missing)
  CS%id_del_sh = register_diag_field('ice_model', 'del_sh', diag%axesT1, Time, &
            'ice strain rate magnitude', 's-1', missing_value=missing)
  CS%id_del_sh_min = register_diag_field('ice_model', 'del_sh_min', diag%axesT1, Time, &
            'minimum ice strain rate magnitude', 's-1', missing_value=missing)

  CS%id_ui_hifreq = register_diag_field('ice_model', 'ui_hf', diag%axesCu1, Time, &
            'ice velocity - x component', 'm/s', missing_value=missing)
  CS%id_vi_hifreq = register_diag_field('ice_model', 'vi_hf', diag%axesCv1, Time, &
            'ice velocity - y component', 'm/s', missing_value=missing)
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
  CS%id_sigi_hifreq  = register_diag_field('ice_model','sigI_hf' ,diag%axesT1, Time,         &
            'first stress invariant', 'none', missing_value=missing)
  CS%id_sigii_hifreq = register_diag_field('ice_model','sigII_hf' ,diag%axesT1, Time,        &
            'second stress invariant', 'none', missing_value=missing)

end subroutine ice_C_dyn_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! find_ice_strength - magnitude of force on ice in plastic deformation         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine find_ice_strength(hi, ci, ice_strength, G, CS, halo_sz) ! ??? may change to do loop
  type(sea_ice_grid_type),          intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: hi, ci
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: ice_strength
  type(ice_C_dyn_CS),               pointer     :: CS
  integer,                optional, intent(in)  :: halo_sz

  integer :: i, j, isc, iec, jsc, jec, halo
  halo = 0 ; if (present(halo_sz)) halo = halo_sz
  isc = G%isc-halo ; iec = G%iec+halo ; jsc = G%jsc-halo ; jec = G%jec+halo

  do j=jsc,jec ; do i=isc,iec
    ice_strength(i,j) = CS%p0*hi(i,j)*ci(i,j)*exp(-CS%c0*(1-ci(i,j)))
  enddo ; enddo

end subroutine find_ice_strength

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_C_dynamics - take a single dynamics timestep with EVP subcycles          !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_C_dynamics(ci, hs, hi, ui, vi, uo, vo, &
                          fxat, fyat, sea_lev, fxoc, fyoc, dt_slow, G, CS)

  type(sea_ice_grid_type), intent(inout) :: G
  real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: ci, hs, hi  ! ice properties
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
  type(ice_C_dyn_CS),                pointer       :: CS
! Arguments: ci - The sea ice concentration, nondim.
!  (in)      hs - The thickness of the snow, in m.
!  (in)      hi - The thickness of the ice, in m.
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
    mi, &       ! Total snow and ice mass per unit area, in kg m-2.
    pres_ice, & ! The ice internal pressure, in N/kg.
    zeta, &     ! The ice bulk viscosity, in ???
    del_sh, &   ! The magnitude of the shear rates, in s-1.
    diag_val, & ! A temporary diagnostic array.
    del_sh_min, &   ! The minimum value of del_sh that is used in the calculation
                    ! of zeta, in s-1.  This is set based on considerations of
                    ! numerical stability, and varies with the grid spacing.
    dx2T, dy2T, &   ! dx^2 or dy^2 at T points, in m2.
    dx_dyT, dy_dxT  ! dx/dy or dy_dx at T points, nondim.

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    fxic, &  ! Zonal force due to internal stresses, in Pa.
    fxic_d, fxic_t, fxic_s, &
    Cor_u, & ! Zonal Coriolis acceleration, in m s-2.
    PFu, &   ! Zonal hydrostatic pressure driven acceleration, in m s-2.
    u_tmp, & ! A temporary copy of the old values of ui, in m s-1.
    mi_u, &  ! The total ice and snow mass interpolated to u points, in kg m-2.
    I_mi_u, &! The inverse of mi_u (plus a tiny mass to avoid NaNs), in m2 kg-1.
    f2dt_u, &! The squared effective Coriolis parameter at u-points times a
             ! time step, in s-1.
    I1_f2dt2_u  ! 1 / ( 1 + f^2 dt^2) at u-points, nondimensional.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    fyic, &  ! Meridional force due to internal stresses, in Pa.
    fyic_d, fyic_t, fyic_s, &
    Cor_v, &  ! Meridional Coriolis acceleration, in m s-2.
    PFv, &   ! Meridional hydrostatic pressure driven acceleration, in m s-2.
    mi_v, &  ! The total ice and snow mass interpolated to v points, in kg m-2.
    I_mi_v, &! The inverse of mi_v (plus a tiny mass to avoid NaNs), in m2 kg-1.
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
  real :: fxic_now, fyic_now  ! ice internal stress accelerations, in m s-2.
  real :: drag_u, drag_v      ! Drag rates with the ocean at u & v points, in s-1.
  real :: tot_area  ! The sum of the area of the four neighboring cells, in m2.
  real :: dxharm    ! The harmonic mean of the x- and y- grid spacings, in m.
  real :: muq2, mvq2  ! The product of the u- and v-face masses per unit cell
                      ! area surrounding a vorticity point, in kg2 m-4.
  real :: muq, mvq    ! The u- and v-face masses per unit cell area extrapolated
                      ! to a vorticity point on the coast, in kg m-2.
  real :: stress_mag  ! The magnitude of the stress at a point.
  real :: str_d_q     ! CS%str_d interpolated to a vorticity point, in Pa m.
  real :: str_t_q     ! CS%str_t interpolated to a vorticity point, in Pa m.
  real :: rescale_str ! A factor by which to rescale the internal stresses, ND.
  real :: pr_vol      ! The internal ice pressure per unit volume, in Pa.
  real :: pres_sum    ! The sum of the internal ice pressures aroung a point, in Pa.
  real :: pres_avg    ! The average of the internal ice pressures around a point, in Pa.
  real :: min_rescale ! The smallest of the 4 surrounding values of rescale, ND.
  real :: I_1pdt_T    ! 1.0 / (1.0 + dt_2Tdamp)
  real :: I_1pE2dt_T  ! 1.0 / (1.0 + EC^2 * dt_2Tdamp)

  real :: Tdamp   ! The damping timescale of the stress tensor components
                  ! toward their equilibrium solution due to the elastic terms,
                  ! in s.
  real :: dt      ! The short timestep associated with the EVP dynamics, in s.
  real :: dt_2Tdamp ! The ratio of the timestep to the elastic damping timescale.
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
  real :: m_neglect3 ! A tiny mass per unit area cubed, in kg3 m-6.
  real :: m_neglect4 ! A tiny mass per unit area to the 4th power, in kg4 m-8.
  real :: sum_area   ! The sum of ocean areas around a vorticity point, in m2.

  type(time_type) :: &
    time_it_start, &  ! The starting time of the iteratve steps.
    time_step_end, &  ! The end time of an iterative step.
    time_end_in       ! The end time for diagnostics when this routine started.
  real :: time_int_in ! The diagnostics' time interval when this routine started.
  logical :: do_hifreq_output  ! If true, output occurs every iterative step.

  integer :: halo_sh_Ds  ! The halo size that can be used in calculating sh_Ds.
  integer :: i, j, isc, iec, jsc, jec, n
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (.not.associated(CS)) call SIS_error(FATAL, &
         "ice_C_dynamics: Module must be initialized before it is used.")

  if ((isc - G%isdB < 2) .or. (jsc - G%jsdB < 2)) call SIS_error(FATAL, &
         "ice_C_dynamics is written to require a 2-point halo or 1-point and symmetric memory.")

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

  if (CS%evp_sub_steps==0) return

  dt = dt_slow/CS%evp_sub_steps
  EC2 = CS%EC**2
  I_EC = 0.0 ; if (CS%EC > 0.0) I_EC = 1.0 / CS%EC
  I_2EC = 0.0 ; if (CS%EC > 0.0) I_2EC = 0.5 / CS%EC
  I_EC2 = 0.0 ; if (EC2 > 0.0) I_EC2 = 1.0 / EC2

  do_hifreq_output = .false.
  if ((CS%id_ui_hifreq > 0) .or. (CS%id_vi_hifreq > 0) .or. &
      (CS%id_str_d_hifreq > 0) .or. (CS%id_str_t_hifreq > 0) .or. &
      (CS%id_str_s_hifreq > 0) .or. (CS%id_sh_d_hifreq > 0) .or. &
      (CS%id_sh_t_hifreq > 0) .or. (CS%id_sh_s_hifreq > 0)) then
    do_hifreq_output = query_SIS_averaging_enabled(CS%diag, time_int_in, time_end_in)
    if (do_hifreq_output) &
      time_it_start = time_end_in - set_time(int(floor(dt_slow+0.5)))
  endif

  Tdamp = CS%Tdamp
  if (CS%Tdamp <= 0.0) then
    ! Hunke (2001) chooses a specified multiple of dt_slow for Tdamp, and shows
    ! that stability requires Tdamp > 2*dt.
    Tdamp = max(0.36*dt_slow, 3.0*dt)
  endif
  dt_2Tdamp = dt / (2.0 * Tdamp)

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    ! Store the total snow and ice mass.
    mi(i,j) = ci(i,j)*(hi(i,j)*CS%Rho_ice + hs(i,j)*CS%Rho_snow)

    ! Precompute pres_ice and the minimum value of del_sh for stability.
    pr_vol = CS%p0*exp(-CS%c0*(1-ci(i,j)))
    pres_ice(i,j) = pr_vol * (hi(i,j)*ci(i,j))

    dxharm = 2.0*G%dxT(i,j)*G%dyT(i,j) / (G%dxT(i,j) + G%dyT(i,j))
    !   Setting a minimum value of del_sh is sufficient to guarantee numerical
    ! stability of the overall time-stepping for the velocities and stresses.
    ! Setting a minimum value of the shear magnitudes is equivalent to setting
    ! a maximum value of the effective lateral viscosities.
    ! I think that this is stable when CS%del_sh_min_scale >= 1.  -RWH
    del_sh_min(i,j) = (2.0 * CS%del_sh_min_scale * pr_vol * dt**2) / &
                      (Tdamp * CS%Rho_ice * dxharm**2)
  enddo ; enddo

  ! Ensure that the input stresses are not larger than could be justified by
  ! the ice pressure now, as the ice might have melted or been advected away
  ! during the thermodynamic and transport phases.
 
!    This commented out version seems to work, but is not obviously better than
! treating each component separately, and the later is simpler.
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
!      pres_avg = ((G%areaT(i,j)*pres_ice(i,j) + G%areaT(i+1,j+1)*pres_ice(i+1,j+1)) + &
!                  (G%areaT(i+1,j)*pres_ice(i+1,j) + G%areaT(i,j+1)*pres_ice(i,j+1))) / &
!                  sum_area
!      if (pres_avg <= 0.0) then
!        CS%str_s(I,J) = 0.0
!      else
!        str_d_q = 0.25 * ((CS%str_d(i,j) + CS%str_d(i+1,j+1)) + &
!                          (CS%str_d(i+1,j) + CS%str_d(i,j+1)))
!        str_t_q = 0.25 * ((CS%str_t(i,j) + CS%str_t(i+1,j+1)) + &
!                          (CS%str_t(i+1,j) + CS%str_t(i,j+1)))
!        ! The factor of 2 here arises because of the definitions of the str_#.
!        stress_mag = 2.0 * sqrt(min(0.0, str_d_q + 0.5*pres_avg)**2 + &
!                                EC2 * (CS%str_s(I,J)**2 + str_t_q**2))
!        if ((stress_mag > pres_avg) .and. G%Lmask2dBu(I,j)) &
!          CS%str_s(I,J) = CS%str_s(I,J) * (pres_avg / stress_mag)
!      endif
!    endif   
!  enddo ; enddo

!  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
!    ! Rescale str_d and str_t without regard to the values of str_s.
!    if (pres_ice(i,j) <= 0.0) then
!      CS%str_d(i,j) = 0.0 ; CS%str_t(i,j) = 0.0
!    elseif (CS%str_d(i,j) < -0.5*pres_ice(i,j)) then
!    ! The factor of 2 here arises because of the definitions of str_d and str_t.
!      stress_mag = 2.0 * sqrt((CS%str_d(i,j) + 0.5*pres_ice(i,j))**2 + &
!                              EC2 * CS%str_t(i,j)**2)
!      if (stress_mag > pres_ice(i,j)) then
!        rescale_str = pres_ice(i,j) / stress_mag
!        CS%str_d(i,j) = rescale_str * (CS%str_d(i,j) + 0.5*pres_ice(i,j)) - &
!                        0.5*pres_ice(i,j)
!        CS%str_t(i,j) = CS%str_t(i,j) * rescale_str
!      endif
!    elseif (CS%EC*abs(CS%str_t(i,j)) > 0.5*pres_ice(i,j)) then
!      ! Only reduce excessively large values of str_t, but do not increase the
!      ! magnitude of str_d.
!      CS%str_t(i,j) = sign(I_2EC*pres_ice(i,j), CS%str_t(i,j))
!    endif
!  enddo ; enddo

  ! The rescaling here is done separately for each component.
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if (CS%str_d(i,j) < -pres_ice(i,j)) CS%str_d(i,j) = -pres_ice(i,j)
    if (CS%EC*CS%str_t(i,j) > 0.5*pres_ice(i,j)) CS%str_t(i,j) = I_2EC*pres_ice(i,j)
    if (CS%EC*CS%str_t(i,j) < -0.5*pres_ice(i,j)) CS%str_t(i,j) = -I_2EC*pres_ice(i,j)
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
      pres_avg = ((G%areaT(i,j)*pres_ice(i,j) + G%areaT(i+1,j+1)*pres_ice(i+1,j+1)) + &
                  (G%areaT(i+1,j)*pres_ice(i+1,j) + G%areaT(i,j+1)*pres_ice(i,j+1))) / &
                  sum_area
    if (CS%EC*CS%str_s(I,J) > 0.5*pres_avg) CS%str_s(I,J) = I_2EC*pres_avg
    if (CS%EC*CS%str_s(I,J) < -0.5*pres_avg) CS%str_s(I,J) = -I_2EC*pres_avg
  enddo ; enddo

  ! Zero out ice velocities with no mass.
  do j=jsc,jec ; do I=isc-1,iec
    if (G%mask2dCu(I,j) * (mi(i,j)+mi(i+1,j)) == 0.0) ui(I,j) = 0.0
  enddo ; enddo
  do J=jsc-1,jec ; do i=isc,iec
    if (G%mask2dCv(i,J) * (mi(i,j)+mi(i,j+1)) == 0.0) vi(I,j) = 0.0
  enddo ; enddo

  do J=jsc-1,jec ; do I=isc-1,iec
    dx2B(I,J) = G%dxBu(I,J)*G%dxBu(I,J) ; dy2B(I,J) = G%dyBu(I,J)*G%dyBu(I,J)
  enddo ; enddo
  do J=jsc-2,jec+1 ; do I=isc-2,iec+1
    dx_dyB(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; dy_dxB(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
  enddo ; enddo
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    dx2T(i,j) = G%dxT(i,j)*G%dxT(i,j) ; dy2T(i,j) = G%dyT(i,j)*G%dyT(i,j)
    dx_dyT(i,j) = G%dxT(i,j)*G%IdyT(i,j) ; dy_dxT(i,j) = G%dyT(i,j)*G%IdxT(i,j)
  enddo ; enddo

  m_neglect = H_subroundoff*CS%Rho_ice
  m_neglect2 = m_neglect**2 ; m_neglect3 = m_neglect**3 ; m_neglect4 = m_neglect**4
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
    elseif (G%Lmask2dBu(I,J)) then
      ! This is an interior ocean point.
      !   Determine an appropriately averaged mass on q-points. The following
      ! expression for mi_q is mi when the masses are all equal, and goes to 4
      ! times the smallest mass averaged onto the 4 adjacent velocity points.  It
      ! comes from taking the harmonic means of the harmonic means of the
      ! arithmetic mean masses at the velocity points.  mi_ratio goes from 4 times
      ! the ratio of the smallest mass at velocity points over the largest mass
      ! at velocity points up to 1.
      muq2 = 0.25 * (mi(i,j) + mi(i+1,j)) * (mi(i,j+1) + mi(i+1,j+1))
      mvq2 = 0.25 * (mi(i,j) + mi(i,j+1)) * (mi(i+1,j) + mi(i+1,j+1))
      mi_ratio_A_q(I,J) = 32.0 * muq2 * mvq2 / ((m_neglect4 + (muq2 + mvq2) * &
                 ((mi(i,j) + mi(i+1,j+1)) + (mi(i,j+1) + mi(i+1,j)))**2) * &
                 sum_area)
    elseif ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
            (G%mask2dCv(i,J) + G%mask2dCu(i+1,J)) > 1.5) then
      !   This is a corner point, and there are 1 unmasked u-point and 1 v-point.
      ! The ratio below goes from 4 times the ratio of the smaller of the two
      ! masses at velocity points over the larger up to 1.
      muq = 0.5 * (G%mask2dCu(I,j) * (mi(i,j) + mi(i+1,j)) + &
                   G%mask2dCu(I,j+1) * (mi(i,j+1) + mi(i+1,j+1)) )
      mvq = 0.5 * (G%mask2dCv(i,J) * (mi(i,j) + mi(i,j+1)) + &
                   G%mask2dCv(i+1,J) * (mi(i+1,j) + mi(i+1,j+1)) )
      mi_ratio_A_q(I,J) = 4.0 * muq * mvq / ((m_neglect2 + (muq + mvq)**2) * sum_area)
    else
      ! This is a straight coastline or all neighboring velocity points are
      ! masked out.  In any case, with just 1 point, the ratio is always 1.
      mi_ratio_A_q(I,J) = 1.0 / sum_area
    endif
  enddo ; enddo

  do j=jsc-1,jec+1 ; do I=isc-1,iec
    mi_u(I,j) = 0.5*(mi(i+1,j) + mi(i,j))
    I_mi_u(I,j) = G%mask2dCu(I,j) / (mi_u(I,j) + m_neglect)
  enddo ; enddo

  do J=jsc-1,jec ; do i=isc-1,iec+1
    mi_v(i,J) = 0.5*(mi(i,j+1) + mi(i,j))
    I_mi_v(i,J) = G%mask2dCv(i,J) / (mi_v(i,J) + m_neglect)
  enddo ; enddo

  do J=jsc-1,jec ; do I=isc-1,iec
    tot_area = (G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i+1,j) + G%areaT(i,j+1))
    q(I,J) = G%CoriolisBu(I,J) * tot_area / &
         (((G%areaT(i,j) * mi(i,j) + G%areaT(i+1,j+1) * mi(i+1,j+1)) + &
           (G%areaT(i+1,j) * mi(i+1,j) + G%areaT(i,j+1) * mi(i,j+1))) + tot_area * m_neglect)
  enddo ; enddo

  ! sea level slope force
  do j=jsc,jec ; do I=isc-1,iec
  enddo ; enddo
  do J=jsc-1,jec ; do i=isc,iec
  enddo ; enddo

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

  if (CS%debug) then
    call uchksum(PFu, "PFu in ice_C_dynamics", G)
    call vchksum(PFv, "PFv in ice_C_dynamics", G)

    call uchksum(ui, "ui pre-steps ice_C_dynamics", G)
    call vchksum(vi, "vi pre-steps ice_C_dynamics", G)
  endif
  if (CS%debug_redundant) then
    call check_redundant_C("PFu/PFv in ice_C_dynamics", PFu, PFv, G)
    call check_redundant_C("fxat/fyat in ice_C_dynamics", fxat, fyat, G)
    call check_redundant_C("uo/vo in ice_C_dynamics",uo, vo, G)
    call check_redundant_C("ui/vi pre-steps ice_C_dynamics",ui, vi, G)
  endif

  ! Do the iterative time steps.
  do n=1,CS%evp_sub_steps

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
    do J=jsc-halo_sh_Ds,jec+halo_sh_Ds-1 ; do I=isc-halo_sh_Ds,iec+halo_sh_Ds-1
      ! This uses a no-slip boundary condition.
      sh_Ds(I,J) = (2.0-G%mask2dBu(I,J)) * &
          (dx_dyB(I,J)*(ui(I,j+1)*G%IdxCu(I,j+1) - ui(I,j)*G%IdxCu(I,j)) + &
           dy_dxB(I,J)*(vi(i+1,J)*G%IdyCv(i+1,J) - vi(i,J)*G%IdyCv(i,J)))
    enddo ; enddo
    if (halo_sh_Ds < 2) call pass_var(sh_Ds, G%Domain, position=BGRID_NE)

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

   ! calculate viscosities - how often should we do this ?

    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      ! Averaging the squared shearing strain is larger than squaring
      ! the averaged strain.  I don't know what is better. -RWH
      del_sh(i,j) = sqrt(sh_Dd(i,j)**2 + I_EC2 * (sh_Dt(i,j)**2 + &
                   (0.25 * ((sh_Ds(I-1,J-1) + sh_Ds(I,j)) + &
                            (sh_Ds(I-1,J) + sh_Ds(I-1,j))))**2 ) ) ! H&D eqn 9

      zeta(i,j) = 0.5*pres_ice(i,j) / max(del_sh(i,j), del_sh_min(i,j))
    enddo ; enddo

    ! Step the stress component equations semi-implicitly.
    I_1pdt_T = 1.0 / (1.0 + dt_2Tdamp)
    I_1pE2dt_T = 1.0 / (1.0 + EC2*dt_2Tdamp)
    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      ! This expression uses that Pres=2*del_sh*zeta with an elliptic yield curve.
      CS%str_d(i,j) = I_1pdt_T * ( CS%str_d(i,j) + dt_2Tdamp * &
                  ( zeta(i,j) * (sh_Dd(i,j) - del_sh(i,j)) ) )
      CS%str_t(i,j) = I_1pdt_T * ( CS%str_t(i,j) + (I_EC2 * dt_2Tdamp) * &
                  ( zeta(i,j) * sh_Dt(i,j) ) )
    enddo ; enddo

    do J=jsc-1,jec ; do I=isc-1,iec
      ! zeta is already set to 0 over land.
      CS%str_s(I,J) = I_1pdt_T * ( CS%str_s(I,J) + (I_EC2 * dt_2Tdamp) * &
                  ( ((G%areaT(i,j)*zeta(i,j) + G%areaT(i+1,j+1)*zeta(i+1,j+1)) + &
                     (G%areaT(i+1,j)*zeta(i+1,j) + G%areaT(i,j+1)*zeta(i,j+1))) * &
                   mi_ratio_A_q(I,J) * sh_Ds(I,J) ) )
    enddo ; enddo


    ! Save the current values of u for later use in updating v.
    do I=isc-1,iec
      u_tmp(I,jsc-1) = ui(I,jsc-1) ; u_tmp(I,jec+1) = ui(I,jec+1) ;
    enddo
    do j=jsc,jec ; do I=isc-1,iec
      ! Save the current values of u for later use in updating v.
      u_tmp(I,j) = ui(I,j)

      Cor = ((azon(I,j) * vi(i+1,J) + czon(I,j) * vi(i,J-1)) + &
             (bzon(I,j) * vi(i,J) + dzon(I,j) * vi(i+1,J-1))) ! - Cor_ref_u(I,j)
      !  Evaluate 1/m x.Div(m strain).  This expressions include all metric terms
      !  for an orthogonal grid.  The str_d term integrates out to no curl, while
      !  str_s & str_t terms impose no divergence and do not act on solid body rotation.
      fxic_now = (G%IdxCu(I,j) * (CS%str_d(i+1,j) - CS%str_d(i,j)) + &
            (G%IdyCu(I,j)*(dy2T(i+1,j)*CS%str_t(i+1,j) - &
                           dy2T(i,j)  *CS%str_t(i,j)) + &
             G%IdxCu(I,j)*(dx2B(I,J)  *CS%str_s(I,J) - &
                           dx2B(I,J-1)*CS%str_s(I,J-1)) ) * G%IareaCu(I,j) ) * &
                  I_mi_u(i,j)
      drag_u = CS%cdw * CS%Rho_ocean * sqrt((ui(I,j)-uo(I,j))**2 + 0.25 * &
                  (((vi(i,J)-vo(i,J))**2 + (vi(i+1,J-1)-vo(i+1,J-1))**2) + &
                   ((vi(i+1,J)-vo(i+1,J))**2 + (vi(i,J-1)-vo(i,J-1))**2)) + &
                  CS%drag_bg_vel2 ) * I_mi_u(I,j)

      ! Determine the Coriolis acceleration and sum for averages...
      Cor_u(I,j) = Cor_u(I,j) + (Cor  - f2dt_u(I,j) * ui(I,j)) * I1_f2dt2_u(I,j)

      !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
      ! update of the other terms and an implicit bottom drag calculation.
      ui(I,j) = G%mask2dCu(I,j) * ( (ui(I,j) + dt * Cor) * I1_f2dt2_u(I,j) + &
                 dt * ((PFu(I,j) + fxat(I,j)*I_mi_u(I,j)) + &
                       (fxic_now + drag_u * uo(I,j))) ) / &
                (1.0 + dt * drag_u)

      ! sum accelerations to take averages.
      fxic(I,j) = fxic(I,j) + fxic_now*mi_u(I,j)
      ! Note that fxoc is the stress felt by the ocean.
      fxoc(I,j) = fxoc(I,j) - drag_u*(uo(I,j) - ui(I,j))*mi_u(I,j)
      if (CS%id_fix_d>0) fxic_d(I,j) = fxic_d(I,j) + mi_u(I,j) * &
               ( (G%IdyCu(I,j)*(dy2T(i+1,j)*(CS%str_d(i+1,j)) - &
                                dy2T(i,j)  *(CS%str_d(i,j))) ) * &
                 (G%IareaCu(I,j) * I_mi_u(i,j)) )
      if (CS%id_fix_t>0) fxic_t(I,j) = fxic_t(I,j) + mi_u(I,j) * &
               ( (G%IdyCu(I,j)*(dy2T(i+1,j)*(CS%str_t(i+1,j)) - &
                                dy2T(i,j)  *(CS%str_t(i,j))) ) * &
                 (G%IareaCu(I,j) * I_mi_u(i,j)) )
      if (CS%id_fix_s>0) fxic_s(I,j) = fxic_s(I,j) + mi_u(I,j) * &
               ( (G%IdxCu(I,j)*(dx2B(I,J)  *CS%str_s(I,J) - &
                                dx2B(I,J-1)*CS%str_s(I,J-1)) ) * &
                 (G%IareaCu(I,j) * I_mi_u(i,j)) )
    enddo ; enddo
    do J=jsc-1,jec ; do i=isc,iec
      Cor = -1.0*((amer(I-1,j) * u_tmp(I-1,j) + cmer(I,j+1) * u_tmp(I,j+1)) + &
                  (bmer(I,j) * u_tmp(I,j) + dmer(I-1,j+1) * u_tmp(I-1,j+1)))
      !  Evaluate 1/m y.Div(m strain).  This expressions include all metric terms
      !  for an orthogonal grid.  The str_d term integrates out to no curl, while
      !  str_s & str_t terms impose no divergence and do not act on solid body rotation.
      fyic_now = (G%IdyCv(i,J) * (CS%str_d(i,j+1)-CS%str_d(i,j)) + &
            (-G%IdxCv(i,J)*(dx2T(i,j+1)*CS%str_t(i,j+1) - &
                            dx2T(i,j)  *CS%str_t(i,j)) + &
              G%IdyCv(i,J)*(dy2B(I,J)  *CS%str_s(I,J) - &
                            dy2B(I-1,J)*CS%str_s(I-1,J)) )*G%IareaCv(i,J) ) * &
                  I_mi_v(i,j)
      drag_v = CS%cdw*CS%Rho_ocean * sqrt((vi(i,J)-vo(i,J))**2 + 0.25 * &
                  (((u_tmp(I,j)-uo(I,j))**2 + (u_tmp(I-1,j+1)-uo(I-1,j+1))**2) + &
                   ((u_tmp(I,j+1)-uo(I,j+1))**2 + (u_tmp(I-1,j)-uo(I-1,j))**2)) + &
                  CS%drag_bg_vel2 )  * I_mi_v(i,J)

      ! Determine the Coriolis acceleration and sum for averages...
      Cor_v(I,J) = Cor_v(I,J) + (Cor - f2dt_v(i,J) * vi(i,J)) * I1_f2dt2_v(i,J)

      !   This is a quasi-implicit timestep of Coriolis, followed by an explicit
      ! update of the other terms and an implicit bottom drag calculation.
      vi(i,J) = G%mask2dCv(i,J) * ((vi(i,J) + dt * Cor) * I1_f2dt2_v(i,J) + &
                 dt * ((PFv(i,J) + fyat(i,J)*I_mi_v(i,J)) + &
                       (fyic_now + drag_v * vo(i,J))) ) / &
                (1.0 + dt * drag_v)

      ! sum accelerations to take averages.
      fyic(i,J) = fyic(i,J) + fyic_now*mi_v(i,J)
      ! Note that fyoc is the stress felt by the ocean.
      fyoc(i,J) = fyoc(i,J) - drag_v*(vo(i,J) - vi(i,J))*mi_v(i,J)

      if (CS%id_fiy_d>0) fyic_d(i,J) = fyic_d(i,J) + mi_v(i,J) * &
               ( (G%IdxCv(i,J)*(dx2T(i,j+1)*(CS%str_d(i,j+1)) - &
                                dx2T(i,j)  *(CS%str_d(i,j))) ) * &
                 (G%IareaCv(i,J) * I_mi_v(i,j)) )
      if (CS%id_fiy_t>0) fyic_t(i,J) = fyic_t(i,J) + mi_v(i,J) * &
               ( (G%IdxCv(i,J)*(dx2T(i,j+1)*(-CS%str_t(i,j+1)) - &
                                dx2T(i,j)  *(-CS%str_t(i,j))) ) * &
                 (G%IareaCv(i,J) * I_mi_v(i,j)) )
      if (CS%id_fiy_s>0) fyic_s(i,J) = fyic_s(i,J) + mi_v(i,J) * &
               ( (G%IdyCv(i,J)*(dy2B(I,J)  *CS%str_s(I,J) - &
                                dy2B(I-1,J)*CS%str_s(I-1,J)) ) * &
                 (G%IareaCv(i,J) * I_mi_v(i,j)) )
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
        call find_sigI(hi, ci, CS%str_d, diag_val, G, CS)
        call post_SIS_data(CS%id_sigi_hifreq, diag_val, CS%diag, mask=G%Lmask2dT)
      endif
      if (CS%id_sigii_hifreq>0) then
        call find_sigII(hi, ci, CS%str_t, CS%str_s, diag_val, G, CS)
        call post_SIS_data(CS%id_sigii_hifreq, diag_val, CS%diag, mask=G%Lmask2dT)
      endif
    endif

    if (CS%debug) then
      call hchksum(CS%str_d, "str_d in ice_C_dynamics", G, haloshift=1)
      call hchksum(CS%str_t, "str_t in ice_C_dynamics", G, haloshift=1)
      call Bchksum(CS%str_s, "str_s in ice_C_dynamics", G, haloshift=1)

      call uchksum(fxic, "fxic in ice_C_dynamics", G)
      call vchksum(fyic, "fyic in ice_C_dynamics", G)
      call uchksum(fxoc, "fxoc in ice_C_dynamics", G)
      call vchksum(fyoc, "fyoc in ice_C_dynamics", G)
      call uchksum(Cor_u, "Cor_u in ice_C_dynamics", G)
      call vchksum(Cor_v, "Cor_v in ice_C_dynamics", G)
      call uchksum(ui, "ui in ice_C_dynamics", G)
      call vchksum(vi, "vi in ice_C_dynamics", G)
    endif
    if (CS%debug_redundant) then
      call check_redundant_C("fxic/fyic in ice_C_dynamics steps",fxic, fyic, G)
      call check_redundant_C("Cor_u/Cor_v in ice_C_dynamics steps", Cor_u, Cor_v, G)
      call check_redundant_C("fxoc in ice_C_dynamics steps", fxoc, fyoc, G)
      call check_redundant_C("ui/vi in ice_C_dynamics steps", ui, vi, G)
    endif

  enddo ! l=1,CS%evp_sub_steps

  if (CS%debug) then
    call uchksum(ui, "ui end ice_C_dynamics", G)
    call vchksum(vi, "vi end ice_C_dynamics", G)
  endif
  if (CS%debug_redundant) &
    call check_redundant_C("ui/vi end ice_C_dynamics", ui, vi, G)

  ! Reset the time information in the diag type.
  if (do_hifreq_output) call enable_SIS_averaging(time_int_in, time_end_in, CS%diag)

  ! make averages
  I_sub_steps = 1.0/CS%evp_sub_steps
  do j=jsc,jec ; do I=isc-1,iec
    fxoc(I,j) = fxoc(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic(I,j) = fxic(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    Cor_u(I,j) = Cor_u(I,j) * (G%mask2dCu(I,j) * I_sub_steps)

    fxic_d(I,j) = fxic_d(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic_t(I,j) = fxic_t(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic_s(I,j) = fxic_s(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
  enddo ; enddo
  do J=jsc-1,jec ; do i=isc,iec
    fyoc(i,J) = fyoc(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic(i,J) = fyic(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    Cor_v(i,J) = Cor_v(i,J) * (G%mask2dCv(i,J) * I_sub_steps)

    fyic_d(i,J) = fyic_d(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic_t(i,J) = fyic_t(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic_s(i,J) = fyic_s(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
  enddo ; enddo

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
      call find_sigI(hi, ci, CS%str_d, diag_val, G, CS)
      call post_SIS_data(CS%id_sigi, diag_val, CS%diag, mask=G%Lmask2dT)
    endif
    if (CS%id_sigii>0) then
      call find_sigII(hi, ci, CS%str_t, CS%str_s, diag_val, G, CS)
      call post_SIS_data(CS%id_sigii, diag_val, CS%diag, mask=G%Lmask2dT)
    endif
    if (CS%id_stren>0) then
      call find_ice_strength(hi, ci, diag_val, G, CS)
      call post_SIS_data(CS%id_stren, diag_val, CS%diag, mask=G%Lmask2dT)
    endif

    if (CS%id_ui>0) call post_SIS_data(CS%id_ui, ui, CS%diag)
    if (CS%id_vi>0) call post_SIS_data(CS%id_vi, vi, CS%diag)

    if (CS%id_str_d>0) call post_SIS_data(CS%id_str_d, CS%str_d, CS%diag)
    if (CS%id_str_t>0) call post_SIS_data(CS%id_str_t, CS%str_t, CS%diag)
    if (CS%id_str_s>0) call post_SIS_data(CS%id_str_s, CS%str_s, CS%diag)

    if (CS%id_sh_d>0) call post_SIS_data(CS%id_sh_d, sh_Dd, CS%diag)
    if (CS%id_sh_t>0) call post_SIS_data(CS%id_sh_t, sh_Dt, CS%diag)
    if (CS%id_sh_s>0) call post_SIS_data(CS%id_sh_s, sh_Ds, CS%diag)

    if (CS%id_del_sh>0) call post_SIS_data(CS%id_del_sh, del_sh, CS%diag)
    if (CS%id_del_sh_min>0) call post_SIS_data(CS%id_del_sh_min, del_sh_min, CS%diag)

  endif

end subroutine ice_C_dynamics

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! sigI - first stress invariant                                                !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine find_sigI(hi, ci, str_d, sigI, G, CS)
  type(sea_ice_grid_type), intent(in)    :: G
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: hi, ci, str_d
  real, dimension(SZI_(G),SZJ_(G)), intent(out) :: sigI
  type(ice_C_dyn_CS),               pointer     :: CS

  real, dimension(SZI_(G),SZJ_(G)) :: &
    strength ! The ice strength, in Pa.
  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  call find_ice_strength(hi, ci, strength, G, CS)

  do j=jsc,jec ; do i=isc,iec
    sigI(i,j) = 0.0
    if (strength(i,j) > 0.0) sigI(i,j) = 2.0 * str_d(i,j) / strength(i,j)
  enddo ; enddo

end subroutine find_sigI

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! sigII - second stress invariant                                              !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine find_sigII(hi, ci, str_t, str_s, sigII, G, CS)
  type(sea_ice_grid_type),            intent(in)  :: G
  real, dimension(SZI_(G),SZJ_(G)),   intent(in)  :: hi, ci, str_t
  real, dimension(SZIB_(G),SZJB_(G)), intent(in)  :: str_s
  real, dimension(SZI_(G),SZJ_(G)),   intent(out) :: sigII
  type(ice_C_dyn_CS),                  pointer    :: CS

  real, dimension(SZI_(G),SZJ_(G)) :: &
    strength ! The ice strength, in Pa.
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    str_s_ss ! Str_s divided by the sum of the neighboring ice strengths.
  real :: strength_sum  ! The sum of the 4 neighboring strengths, in Pa.
  real :: sum_area   ! The sum of ocean areas around a vorticity point, in m2.
  real, parameter :: H_subroundoff = 1e-30 ! A thickness that is so small it is
                    ! usually lost in roundoff and can be neglected, but can be
                    ! cubed without being lost to underflow, in H.
  integer :: i, j, isc, iec, jsc, jec
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    strength(i,j) = (hi(i,j)*ci(i,j)) * CS%p0*exp(-CS%c0*(1-ci(i,j)))
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
! ice_dyn_register_restarts - allocate and register any variables for this     !
!      module that need to be included in the restart files.                   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_C_dyn_register_restarts(G, param_file, CS, Ice_restart, restart_file)
  type(sea_ice_grid_type), intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(ice_C_dyn_CS),      pointer       :: CS
  type(restart_file_type), intent(inout) :: Ice_restart
  character(len=*),        intent(in)    :: restart_file

! Arguments: G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module.
!

!   This subroutine registers the restart variables associated with the
! the ice dynamics.

  integer :: isd, ied, jsd, jed, id
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call SIS_error(WARNING, "ice_dyn_register_restarts called with an "//&
                            "associated control structure.")
    return
  endif
  allocate(CS)

  allocate(CS%str_d(isd:ied, jsd:jed)) ; CS%str_d(:,:) = 0.0
  allocate(CS%str_t(isd:ied, jsd:jed)) ; CS%str_t(:,:) = 0.0
  allocate(CS%str_s(G%IsdB:G%IedB, G%JsdB:G%JedB)) ; CS%str_s(:,:) = 0.0
  id = register_restart_field(Ice_restart, restart_file, 'str_d', CS%str_d, &
                              domain=G%Domain%mpp_domain)
  id = register_restart_field(Ice_restart, restart_file, 'str_t', CS%str_t, &
                              domain=G%Domain%mpp_domain)
  id = register_restart_field(Ice_restart, restart_file, 'str_s', CS%str_s, &
                              domain=G%Domain%mpp_domain, position=CORNER)
end subroutine ice_C_dyn_register_restarts

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_dyn_end - deallocate the memory associated with this module.             !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine ice_C_dyn_end(CS)
  type(ice_C_dyn_CS), pointer :: CS

  deallocate(CS%str_d) ; deallocate(CS%str_t) ; deallocate(CS%str_s)

  deallocate(CS)
end subroutine ice_C_dyn_end

end module ice_dyn_cgrid
