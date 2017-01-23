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

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Andrew Shao 2016                                                *
!*                                                                     *
!*    This file contains an example of the code that is needed to set  *
!*  up and use two passive tracers that represent sea-ice age.         *
!*                                                                     *
!*  Areal age tracer:                                                  *
!*      For this tracer, the age of the ice in each thickness category *
!*  that exceeds a given threshold is increased by the length of the   *
!*  timestep. The maximum of the age tracer may be expected to be      *
!*  similar to observational estimates of ice age.                     *
!*  Reference: Hunke and Bitz [2009], JGR                              *
!*                                                                     *
!*  Mass balance age tracer:                                           *
!*      As before, all ice that exceeds a given proscribed area is     *
!*  increased by the length of the timestep. However, any ice that is  *
!*  formed from seawater and added to existing ice has an age equal    *
!*  to the length of the timestep                                      *
!*  Reference: Lietaer et al. [2011], Ocean Modeling                   *
!*                                                                     *
!*                                                                     *
!*  This file is adapted from ideal_age_example.F90 included as        *
!*  part of the MOM6 repository. This code replaces and incorporates   *
!*  some of the ice age tracer work implemented in SIS by Torge Martin *
!*                                                                     *
!*  Macros written all in capital letters are defined in SIS2_memory.h *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**
module ice_age_tracer
! ashao: Get all the dependencies from other modules (check against
!   ideal_age_example.F90)

use SIS_diag_mediator, only     : register_SIS_diag_field, &
                                  safe_alloc_ptr
use SIS_diag_mediator, only     : SIS_diag_ctrl, post_data=>post_SIS_data
use SIS_tracer_registry, only   : register_SIS_tracer, SIS_tracer_registry_type
use SIS_hor_grid, only          : sis_hor_grid_type
use SIS_utils, only             : post_avg

use MOM_file_parser, only       : get_param, log_param, log_version, param_file_type
use MOM_restart, only           : query_initialized, MOM_restart_CS
use MOM_io, only                : vardesc, var_desc, query_vardesc, file_exists
use MOM_time_manager, only      : time_type, get_time
use MOM_sponge, only            : set_up_sponge_field, sponge_CS
use MOM_error_handler, only     : SIS_error=>MOM_error, FATAL, WARNING
use MOM_error_handler, only     : SIS_mesg=>MOM_mesg
use MOM_string_functions, only  : slasher

use fms_mod, only               : read_data
use fms_io_mod, only            : register_restart_field, restore_state
use fms_io_mod, only            : restart_file_type

use ice_grid, only              : ice_grid_type

implicit none ; private
#include "SIS2_memory.h"

public register_ice_age_tracer, initialize_ice_age_tracer
public ice_age_tracer_column_physics
public ice_age_stock, ice_age_end

integer, parameter :: NTR_MAX = 2

type p3d
  real, dimension(:,:,:), pointer :: p => NULL()
end type p3d

type, public :: ice_age_tracer_CS
  integer :: ntr                              ! The number of tracers that are actually used.
  character(len = 200) :: IC_file             ! The file in which the age-tracer initial values
                                              ! can be found, or an empty string for internal initialization.
  type(time_type), pointer :: Time            ! A pointer to the ocean model's clock.
  type(SIS_tracer_registry_type), pointer :: TrReg => NULL()
  real, pointer :: tr(:,:,:,:,:) => NULL()    ! The array of tracers used in this
                                              ! subroutine, in g m-3?
  real, pointer :: tr_aux(:,:,:,:,:) => NULL()! The masked tracer concentration
                                              ! for output, in g m-3.
  type(p3d), dimension(NTR_MAX) :: &
      tr_adx, &                               ! Tracer zonal advective fluxes in g m-3 m3 s-1.
      tr_ady                                  ! Tracer meridional advective fluxes in g m-3 m3 s-1.

  real, pointer :: ocean_BC(:,:,:,:)=>NULL()  ! Ocean boundary value of the tracer by category
  real, pointer :: snow_BC(:,:,:,:)=>NULL()   ! Snow boundary value of the tracer by category

  real, dimension(NTR_MAX) :: &
      IC_val = 0.0, &                         ! The (uniform) initial condition value.
      young_val = 0.0, &                      ! The value assigned to tr at the surface.
      land_val = -1.0, &                      ! The value of tr used where land is masked out.
      tracer_start_year = 0.0                 ! The year in which tracers start aging, or at which the
                                              ! surface value equals young_val, in years.
  logical :: mask_tracers                     ! If true, tracers are masked out in massless layers.
  logical :: tracers_may_reinit               ! If true, tracers may go through the
                                              ! initialization code if they are not found in the
                                              ! restart files.
  logical :: tracer_ages(NTR_MAX)             ! Flag if the tracer ages as a function of time
  logical :: new_ice_is_sink(NTR_MAX)         ! Flag for whether the formation of new ice acts
                                              ! as a sink for age
  logical :: advect_vertical(NTR_MAX)         ! Whether the tracer should be advected "vertically"
                                              ! to thicker ice
  logical :: uniform_vertical(NTR_MAX)        ! Whether the tracer should uniform across ice thickness
                                              ! categories if so, "mix" the tracer by setting it to
                                              ! the maximum age at the grid pont
  logical :: do_ice_age_areal
  logical :: do_ice_age_mass
  real    :: min_thick_age, min_conc_age

  integer, dimension(NTR_MAX) :: &
      id_tracer = -1, id_tr_adx = -1, id_tr_ady = -1, id_avg =-1 ! Indices used for the diagnostic
                                                                 ! mediator
  integer, dimension(NTR_MAX) :: nlevels

  type(SIS_diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.

  real :: min_mass                            ! Minimum mass of ice in thickness category for it
                                              ! to 'exist'

  type(vardesc) :: tr_desc(NTR_MAX)
end type ice_age_tracer_CS

contains

logical function register_ice_age_tracer(G, IG, param_file, CS, diag, TrReg, &
    Ice_restart, restart_file)
  type(sis_hor_grid_type),                intent(in) :: G
  type(ice_grid_type),                    intent(in) :: IG
  type(param_file_type),                  intent(in) :: param_file
  type(ice_age_tracer_CS),                pointer    :: CS
  type(SIS_diag_ctrl),                    target     :: diag
  type(SIS_tracer_registry_type),         pointer    :: TrReg
  type(restart_file_type),                intent(inout) :: Ice_restart
  character(len=*)                                   :: restart_file
  ! This subroutine is used to age register tracer fields and subroutines
  ! to be used with SIS.
  ! Arguments:
  !  (in)      Ice - The ice data type
  !  (in)      G - The ocean's grid structure.
  !  (in)      IG - Ice model's grid structure
  !  (in)      param_file - A structure indicating the open file to parse for
  !                         model parameter values.
  !  (in/out)  CS - A pointer that is set to point to the control structure
  !                 for the ice age tracer
  !  (in)      diag - A structure that is used to regulate diagnostic output.
  !  (in/out)  TrReg - A pointer that is set to point to the control structure
  !                  for the tracer advection and diffusion module.
  !  (in)      restart_file - Name of the restart file for the ice model.

  ! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "ice_age_tracer" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  real, dimension(:,:,:), pointer :: ocean_BC_ptr, snow_BC_ptr
  integer :: isc, iec, jsc, jec, k, m, tr
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (associated(CS)) then
      call SIS_error(WARNING, "register_ice_age_tracer called with an "// &
          "associated control structure.")
      register_ice_age_tracer = .false.
      return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "DO_ICE_AGE_AREAL", CS%do_ice_age_areal, &
      "If true, use an ice age tracer that ages at \n"//&
      "a unit rate if the ice exists.", &
      default=.false.)
  call get_param(param_file, mod, "DO_ICE_AGE_MASS", CS%do_ice_age_mass, &
      "If true, use an ice age tracer that is set to 0 age \n"//&
      "when ice is first formed, ages at unit rate in the ice pack \n"// &
      "and any new ice added to the ice pack represents a decrease in age", &
      default=.false.)
  call get_param(param_file, mod, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
      "If true, tracers may go through the initialization code \n"//&
      "if they are not found in the restart files.  Otherwise \n"//&
      "it is a fatal error if the tracers are not found in the \n"//&
      "restart files of a restarted run.", default=.false.)

  CS%ntr = 0
  if (CS%do_ice_age_areal) then
      CS%ntr = CS%ntr + 1 ; m = CS%ntr
      CS%tr_desc(m) = var_desc("ice_age_areal", "years", "Area based ice age tracer", caller=mod)
      CS%tracer_ages(m) = .true.
      CS%new_ice_is_sink(m) = .false.
      CS%uniform_vertical(m) = .true.
      CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = -1.0
      CS%land_val(m) = 0.0
      CS%nlevels(m) = IG%NkIce
  endif
  if (CS%do_ice_age_mass) then
      CS%ntr = CS%ntr + 1 ; m = CS%ntr
      CS%tr_desc(m) = var_desc("ice_age_mass", "years", "Mass-based ice age tracer", caller=mod)
      CS%tracer_ages(m) = .true.
      CS%new_ice_is_sink(m) = .true.
      CS%uniform_vertical(m) = .false.
      CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = -1.0
      CS%land_val(m) = 0.0
      CS%nlevels(m) = IG%NkIce
  endif
  CS%min_mass = 1.0e-7*IG%H_to_kg_m2

  ! Allocate the main tracer arrays
  allocate(CS%tr(SZI_(G), SZJ_(G),IG%CatIce,IG%NkIce,CS%ntr)) ; CS%tr(:,:,:,:,:) = 0.0
  ! Boundary condition arrays
  allocate(CS%ocean_BC(SZI_(G), SZJ_(G),IG%CatIce,CS%ntr)) ; CS%ocean_BC(:,:,:,:) = 0.0
  allocate(CS%snow_BC(SZI_(G), SZJ_(G),IG%CatIce,CS%ntr)) ; CS%snow_BC(:,:,:,:)  = 0.0

  ! Make sure that diag manager is assigned
  CS%diag => diag

  do m=1,CS%ntr

    call query_vardesc(CS%tr_desc(m), name=var_name, &
        caller="register_ice_age_tracer")

    ! Register the tracer for the restart file.
    CS%id_tracer(m) = register_restart_field(Ice_restart, restart_file, var_name, &
        CS%tr(:,:,:,1,m), domain=G%domain%mpp_domain, &
        mandatory=.false.)

    ocean_BC_ptr => CS%ocean_BC(:,:,:,m)
    snow_BC_ptr  => CS%snow_BC(:,:,:,m)
    ! Register the tracer for horizontal advection & diffusion. Note that the argument
    ! of nLTr is IG%NkIce
    call register_SIS_tracer(CS%tr(:,:,:,:,m), G, IG, IG%NkIce, var_name, param_file, &
        TrReg, snow_tracer=.false.,ad_3d_x=CS%tr_adx(m)%p, &
        ad_3d_y=CS%tr_ady(m)%p, ocean_BC=ocean_BC_ptr, snow_BC=snow_BC_ptr, is_passive=.true.)

    ! Set boundary conditions
    if(CS%new_ice_is_sink(m)) then
      ! Newly formed ice has zero age
      CS%ocean_BC(:,:,:,m) = 0.0
    else
      ! New ice does not decrease the age (i.e. enters the ice at the same concentration as
      ! For now, all levels within a thickness category should have the same tracer concentration
      ! Equivalently, this statement would mean that new ice has the same age as the bottom-most
      ! level of the ice
      CS%ocean_BC(:,:,:,m) = CS%tr(:,:,:,IG%NkIce,m)
    endif
    CS%snow_BC(:,:,:,m) = CS%tr(:,:,:,1,m)
  enddo ! do m=1, CS%ntr

  CS%TrReg => TrReg
  register_ice_age_tracer = .true.

end function register_ice_age_tracer


subroutine initialize_ice_age_tracer( day, G, IG, CS, is_restart )
  type(time_type), target,                intent(in) :: day
  type(sis_hor_grid_type),                intent(in) :: G
  type(ice_grid_type),                    intent(in) :: IG
  type(ice_age_tracer_CS),                pointer    :: CS
  logical,                                intent(in) :: is_restart

  !   This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
  ! and it sets up the tracer output.

  ! Arguments: restart - .true. if the fields have already been read from
  !                     a restart file.
  !  (in)      day - Time of the start of the run.
  !  (in)      G - The ocean's grid structure.
  !  (in)      IG - The ice model's grid structure.
  !  (in/out)  CS - The control structure returned by a previous call to
  !                 register_ideal_age_tracer.
  real, parameter   :: missing = -1.0e34
  character(len=24) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  integer :: i, j, k, isc, iec, jsc, jec, nz, m, tr
  integer :: IscB, IecB, JscB, JecB


  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  IscB = G%IscB ; IecB = G%IecB ; JscB = G%JscB ; JecB = G%JecB

  CS%Time => day

  do tr=1,CS%ntr
    do m = 1,IG%NkIce ; do k=1,CS%nlevels(tr) ; do j=jsc,jec ; do i=isc,iec
      if (G%mask2dT(i,j) < 0.5) then
        CS%tr(i,j,k,m,tr) = CS%land_val(tr)
      elseif(.not. is_restart) then
        CS%tr(i,j,k,m,tr) = CS%IC_val(tr)
      endif
    enddo ; enddo ; enddo ; enddo
  enddo ! Tracer loop

  ! This needs to be changed if the units of tracer are changed above.
  flux_units = "years kg s-1"

  do tr=1,CS%ntr

    call query_vardesc(CS%tr_desc(tr), name, units=units, longname=longname, &
        caller="initialize_ice_age_tracer")
    CS%id_tracer(tr) = register_SIS_diag_field("ice_model", trim(name), CS%diag%axesTc, &
        CS%Time, trim(longname) , trim(units),missing_value = missing)
    CS%id_tr_adx(tr) = register_SIS_diag_field("ice_model", trim(name)//"_adx", &
        CS%diag%axesCuc, CS%Time, trim(longname)//" advective zonal flux" , &
        trim(flux_units))
    CS%id_tr_ady(tr) = register_SIS_diag_field("ice_model", trim(name)//"_ady", &
        CS%diag%axesCvc, CS%Time, trim(longname)//" advective meridional flux" , &
        trim(flux_units))
    CS%id_avg(tr) = register_SIS_diag_field("ice_model", trim(name)//"_avg", &
        CS%diag%axesT1, CS%Time, trim(longname)//" mass-averaged over all thickness categories" , &
        trim(units))


    if (CS%id_tr_adx(tr) > 0) call safe_alloc_ptr(CS%tr_adx(tr)%p,IscB,IecB,jsc,jec,CS%nlevels(tr))
    if (CS%id_tr_ady(tr) > 0) call safe_alloc_ptr(CS%tr_ady(tr)%p,isc,iec,JscB,JecB,CS%nlevels(tr))

  enddo

end subroutine initialize_ice_age_tracer

subroutine ice_age_tracer_column_physics(dt, G, IG, CS,  mi, mi_old)
  real,                                       intent(in) :: dt
  type(sis_hor_grid_type),                    intent(in) :: G
  type(ice_grid_type),                        intent(in) :: IG
  type(ice_age_tracer_CS)                                :: CS
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),intent(in) :: mi
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),intent(in) :: mi_old

  ! Arguments:
  !  (in)      dt - The amount of time covered by this call, in s.
  !  (in)      mi_old - Mass of ice at the beginning of the ice model timestep
  !  (in)      G - The ocean model's grid structure.
  !  (in)      IG - The ice model's grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 register_ideal_age_tracer.

  real :: Isecs_per_year  ! The number of seconds in a year.
  real :: year            ! The time in years.
  real :: dt_year         ! Timestep in units of years
  real :: min_age         ! Minimum age of ice to avoid being set to 0
  real :: mi_min          ! Minimum mass in ice category
  real :: max_age         ! Maximum age at a grid point
  real, dimension(SZI_(G),SZJ_(G)) :: vertsum_mi, vertsum_mi_old
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)) :: tr_avg
  integer :: secs, days   ! Integer components of the time type.
  integer :: i, j, k, m, tr
  integer :: isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ;

  if (CS%ntr < 1) return

  Isecs_per_year = 1.0 / (365.0*86400.0)
  dt_year = dt * Isecs_per_year

  min_age = 0.0

  call get_time(CS%Time, secs, days)
  year = (86400.0*days + real(secs)) * Isecs_per_year

  do tr=1,CS%ntr
    ! Increment the age of the ice if it exists
    if (CS%tracer_ages(tr) .and. &
        (year>=CS%tracer_start_year(tr))) then
        do j=jsc,jec ; do i=isc,iec
          if(G%mask2dT(i,j)>0.0) then
            do k=1,IG%CatIce; do m=1,CS%nlevels(tr)
              if(CS%tr(i,j,k,m,tr)<min_age) CS%tr(i,j,k,m,tr) = 0.0

              if(mi(i,j,k)>CS%min_mass) then
                CS%tr(i,j,k,m,tr) = CS%tr(i,j,k,m,tr) + dt_year
              else
                CS%tr(i,j,k,m,tr) = 0.0
              endif
            enddo ; enddo
          else
            do m=1,CS%nlevels(tr) ; do k=1,IG%CatIce
              CS%tr(i,j,k,m,tr) = 0.0
            enddo ; enddo
          endif

        enddo ; enddo
    endif

    if (CS%uniform_vertical(tr)) then
      do k=1,IG%CatIce; do j=jsc,jec ; do i=isc,iec
        max_age = maxval(CS%tr(i,j,k,:,tr))
        do m=1,CS%nlevels(tr)
          CS%tr(i,j,k,m,tr) = max_age
        enddo
      enddo ; enddo ; enddo
    endif

    ! Update boundary conditions
    CS%ocean_BC(:,:,:,tr) = CS%tr(:,:,:,IG%NkIce,tr)
    CS%snow_BC(:,:,:,tr) = CS%tr(:,:,:,1,tr)

    ! If levels with different thicknesses are implemented, this averaging
    ! will need to be updated
    tr_avg(:,:,:) = 0.0
    do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
      do m=1,IG%NkIce
        tr_avg(i,j,k) = tr_avg(i,j,k) + CS%tr(i,j,k,m,tr)/IG%NkIce
      enddo
    enddo ; enddo ; enddo

    if (CS%id_tracer(tr)>0) &
        call post_data(CS%id_tracer(tr),CS%tr(:,:,:,1,tr),CS%diag)
    if (CS%id_tr_adx(tr)>0) &
        call post_data(CS%id_tr_adx(tr),CS%tr_adx(tr)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(tr)>0) &
        call post_data(CS%id_tr_ady(tr),CS%tr_ady(tr)%p(:,:,:),CS%diag)
    if (CS%id_avg(tr)>0) then
      call post_avg(CS%id_avg(tr), tr_avg, mi, &
          CS%diag, G=G, wtd=.true.)
    endif
  enddo

end subroutine ice_age_tracer_column_physics

function ice_age_stock(mi, stocks, G, IG, CS, names, units)
  real, dimension(:),                          intent(out)   :: stocks
  type(sis_hor_grid_type),                     intent(in)    :: G
  type(ice_grid_type),                         intent(in)    :: IG
  real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), intent(in)    :: mi
  type(ice_age_tracer_CS), pointer                           :: CS
  character(len=*), dimension(:),              intent(out)   :: names
  character(len=*), dimension(:),              intent(out)   :: units

  integer ice_age_stock
  ! This function calculates the mass-weighted integral of all tracer stocks,
  ! returning the number of stocks it has calculated.  If the stock_index
  ! is present, only the stock corresponding to that coded index is returned.

  ! Arguments: stocks - the mass-weighted integrated amount of each tracer,
  !                     in kg times concentration units.
  !  (in)      G - The ocean's grid structure.
  !  (in)      G - The ice model's grid structure.
  !  (in)      CS - The control structure returned by a previous call to
  !                 register_ideal_age_tracer.
  !  (out)     names - the names of the stocks calculated.
  !  (out)     units - the units of the stocks calculated.
  ! Return value: the number of stocks calculated here.

  integer :: i, j, k, m, tr, nstocks
  integer :: isc, iec, jsc, jec
  real :: avg_tr
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (CS%ntr < 1) return

  do tr=1,CS%ntr
    call query_vardesc(CS%tr_desc(tr), name=names(tr), units=units(tr), caller="ice_age_stock")
    units(tr) = trim(units(tr))//" kg"
    stocks(tr) = 0.0
    do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
      avg_tr = 0.0
      ! For now an equally weighted average over the layers is used to calculate the stocks
      ! In the future, if ice levels have different thicknesses, this will need to be updated
      do m=1,IG%NkIce ; avg_tr = avg_tr + CS%tr(i,j,k,m,tr) ; enddo
      avg_tr = avg_tr/IG%NkIce

      stocks(tr) = stocks(tr) + avg_tr * &
          (G%mask2dT(i,j) * G%areaT(i,j) * mi(i,j,k))
    enddo ; enddo ; enddo

  enddo

  ice_age_stock = CS%ntr

end function ice_age_stock

subroutine ice_age_end(CS)
  type(ice_age_tracer_CS), pointer :: CS
  integer :: m

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    if (associated(CS%tr_aux)) deallocate(CS%tr_aux)
    do m=1,CS%ntr
      if (associated(CS%tr_adx(m)%p)) deallocate(CS%tr_adx(m)%p)
      if (associated(CS%tr_ady(m)%p)) deallocate(CS%tr_ady(m)%p)
    enddo
    deallocate(CS%ocean_BC)
    deallocate(CS%snow_BC)

    deallocate(CS)
  endif
end subroutine ice_age_end

end module ice_age_tracer
