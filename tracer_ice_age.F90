!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the impliec warranty of MERCHANTABILITY  *
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
!*  By Andrew Shao and Robert Hallberg, 2016                           *
!*                                                                     *
!*    This file contains an example of the code that is needed to set  *
!*  up and use two passive tracers that represent sea-ice age.         *
!*                                                                     *
!*  Areal age tracer:                                                  *
!*      For this tracer, the age of the ice that exceeds a given area  *
!*  threshold is increased by the length of the timestep. The entire   *
!*  ice layer in a given grid point is thus uniform in age. This       *
!*  age tracer admits ready comparison to observationally derived      *
!*  estimates of ice age.                                              *
!*  Reference: Hunke and Bitz [2009], JGR                              *
!*                                                                     *
!*  Mass balance age tracer:                                           *
!*      As before, all ice that exceeds a given proscribed aread is    *
!*  increased by the length of the timestep. However, any ice that is  *
!*  formed from seawater and added to existing ice has zero age.       *
!*  Reference: Lietaer et al. [2011], Ocean Modeling                   *
!*                                                                     *
!*    A single subroutine is called from within each file to register  *
!*  each of the tracers for reinitialization and advection and to      *
!*  register the subroutine that initializes the tracers and set up    *
!*  their output and the subroutine that does any tracer physics or    *
!*  chemistry along with diapycnal mixing (included here because some  *
!*  tracers may float or swim vertically or dye diapycnal processes).  *
!*                                                                     *
!*  	This file is adapted from ideal_age_example.F90 included as    *
!*  part of the MOM6 repository.                                       *
!*                                                                     *
!*  For now, because SIS2 lacks an equivalent to tracer_flow_control   *
!*  a separate module exists for the age control structure and the     *
!*  age-related subroutines                                            *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**



module ice_age_tracer
! ashao: Get all the dependencies from other modules (check against
!   ideal_age_example.F90)
!   Also, need to check whether using the MOM6 or SIS2 versions of these modules
!   04082016: shouldn't matter too much right now, except for
!    future changes in static memory cases
use SIS_diag_mediator, only : post_SIS_data, register_SIS_diag_field, &
    safe_alloc_ptr
use SIS_diag_mediator, only: SIS_diag_ctrl, post_data=>post_SIS_data
use SIS_tracer_registry, only : register_SIS_tracer, SIS_tracer_registry_type
use SIS_hor_grid, only : sis_hor_grid_type

use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_io, only : vardesc, var_desc, query_vardesc, file_exists
use MOM_time_manager, only : time_type, get_time
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING
use MOM_error_handler, only : SIS_mesg=>MOM_mesg
use MOM_string_functions, only : slasher

use fms_mod,          only : read_data
use fms_io_mod, only : register_restart_field
use fms_io_mod, only : restart_file_type
use ice_age_tracer_type, only : ice_age_tracer_CS

use ice_grid, only : ice_grid_type

implicit none ; private
#include "SIS2_memory.h"

public register_ice_age_tracer, initialize_ice_age_tracer
public ice_age_tracer_column_physics
public ice_age_stock, ice_age_end

contains

logical function register_ice_age_tracer(G, IG, param_file, CS, diag, TrReg, &
                                   Ice_restart, restart_file)
  type(sis_hor_grid_type),  intent(in) :: G
  type(ice_grid_type),      intent(in) :: IG
  type(param_file_type),    intent(in) :: param_file
  type(ice_age_tracer_CS), pointer     :: CS
  type(SIS_diag_ctrl), target,   intent(in) :: diag
  type(SIS_tracer_registry_type),    pointer    :: TrReg
  type(restart_file_type), intent(inout) :: Ice_restart
  character(len=*)       :: restart_file
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
  real, pointer :: tr_ptr(:,:,:) => NULL()
  integer :: isc, iec, jsc, jec, m, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; NkIce = IG%NkIce

  if (associated(CS)) then
    call SIS_error(WARNING, "register_ice_age_tracer called with an "// &
                             "associated control structure.")
    register_ice_age_tracer = .false.
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "DO_AREAL_ICE_AGE", CS%do_areal_ice_age, &
                 "If true, use an ice age tracer that is set to 0 age \n"//&
                 "when ice is first formed and ages at unit rate in the ice pack.", &
                 default=.true.)
  call get_param(param_file, mod, "DO_MASS_ICE_AGE", CS%do_mass_ice_age, &
                 "If true, use an ice age tracer that is set to 0 age \n"//&
                 "when ice is first formed, ages at unit rate in the ice pack \n"// &
                 "and any new ice added to the ice pack represents a decrease in age", &
                 default=.true.)
  call get_param(param_file, mod, "MINIMUM_THICKNESS_AGE", CS%min_thick_age, &
                 "Minimum thickness of ice that should be considered as first \n"//&
                 "year ice (in m). According to the WMO, this should be 30cm \n.",  &
                 default=0.3)
  call get_param(param_file, mod, "MINIMUM_CONCENTRATION_AGE", CS%min_conc_age, &
                 "Minimum concentration of ice that should be considered to exist. \n"//&
                 "According to Rigor and Wallace this should 0.15 \n.",  &
                 default=0.15)
  call get_param(param_file, mod, "ICE_AGE_IC_FILE", CS%IC_file, &
                 "The file in which the age-tracer initial values can be \n"//&
                 "found, or an empty string for internal initialization.", &
                 default=" ")
  if ((len_trim(CS%IC_file) > 0) .and. (scan(CS%IC_file,'/') == 0)) then
    ! Add the directory if CS%IC_file is not already a complete path.
    call get_param(param_file, mod, "INPUTDIR", inputdir, default=".")
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mod, "INPUTDIR/ICE_AGE_IC_FILE", CS%IC_file)
  endif
  call get_param(param_file, mod, "MASK_MASSLESS_TRACERS", CS%mask_tracers, &
                 "If true, the tracers are masked out in massless layer. \n"//&
                 "This can be a problem with time-averages.", default=.false.)
  call get_param(param_file, mod, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
                 "If true, tracers may go through the initialization code \n"//&
                 "if they are not found in the restart files.  Otherwise \n"//&
                 "it is a fatal error if the tracers are not found in the \n"//&
                 "restart files of a restarted run.", default=.false.)

  CS%ntr = 0
  if (CS%do_areal_ice_age) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("ice_age_areal", "years", "Areal Ice Age Tracer", caller=mod)
    CS%tracer_ages(m) = .true. ;
    CS%new_ice_is_sink = .false. ;
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
    CS%nlevels(m) = 1;
  endif
  if (CS%do_mass_ice_age) then
    CS%ntr = CS%ntr + 1 ; m = CS%ntr
    CS%tr_desc(m) = var_desc("ice_age_mass", "years", "Mass Ice Age Tracer", caller=mod)
    CS%tracer_ages(m) = .true. ;
    CS%new_ice_is_sink = .true. ;
    CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = 0.0
    CS%nlevels(m) = 1;
  endif


  allocate(CS%tr(SZI_(G), SZJ_(G),CS%nlevels(m),CS%ntr)) ; CS%tr(:,:,:,:) = 0.0
  if (CS%mask_tracers) then
    allocate(CS%tr_aux(SZI_(G), SZJ_(G),CS%nlevels(m),CS%ntr)) ; CS%tr_aux(:,:,:,:) = 0.0
  endif
    write(*,*) "Registering for restart"
  do m=1,CS%ntr
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)

    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_ice_age_tracer")

    ! Register the tracer for the restart file.
    CS%id_tracer(m) = register_restart_field(Ice_restart, restart_file, var_name, &
                                     tr_ptr, domain=G%domain%mpp_domain, &
                                     mandatory=.false., read_only=.true.)

    ! Register the tracer for horizontal advection & diffusion.
    call register_SIS_tracer(tr_ptr, G, IG, CS%nlevels(m), var_name, param_file, &
                             TrReg, snow_tracer=.false.)
                             

  enddo ! do m=1, CS%ntr

  CS%TrReg => TrReg
  register_ice_age_tracer = .true.
  write(*,*) "Ice age tracer registered"
end function register_ice_age_tracer


subroutine initialize_ice_age_tracer(restart, day, G, IG, CS )
  logical,                            intent(in) :: restart
  type(time_type), target,               intent(in) :: day
  type(sis_hor_grid_type),              intent(in) :: G
  type(ice_grid_type),              intent(in) :: IG
  type(ice_age_tracer_CS),          pointer    :: CS

!   This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
! and it sets up the tracer output.

! Arguments: restart - .true. if the fields have already been read from
!                     a restart file.
!  (in)      day - Time of the start of the run.
!  (in)      G - The ocean's grid structure.
!  (in)      IG - The ice model's grid structure.
!  (in/out)  CS - The control structure returned by a previous call to
!                 register_ideal_age_tracer.

  character(len=24) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  integer :: i, j, k, isc, iec, jsc, jec, nz, m
  integer :: IscB, IecB, JscB, JecB
  real, pointer :: tr_ptr(:,:,:) => NULL()


  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  IscB = G%IscB ; IecB = G%IecB ; JscB = G%JscB ; JecB = G%JecB

  CS%Time => day

  do m=1,CS%ntr
      call query_vardesc(CS%tr_desc(m), name=name, &
          caller="initialize_ice_age_tracer")
      if ((.not.restart) .or. (CS%tracers_may_reinit)) then

          if (len_trim(CS%IC_file) > 0) then
              !  Read the tracer concentrations from a netcdf file.
              if (.not.file_exists(CS%IC_file, G%Domain)) &
                  call SIS_error(FATAL, "initialize_ice_age_tracer: "// &
                  "Unable to open "//CS%IC_file)
          else
              call read_data(CS%IC_file, trim(name), CS%tr(:,:,:,m), &
                  domain=G%Domain%mpp_domain)
          endif
      else
          do k=1,CS%nlevels(m) ; do j=jsc,jec ; do i=isc,iec
              if (G%mask2dT(i,j) < 0.5) then
                  CS%tr(i,j,k,m) = CS%land_val(m)
              else
                  CS%tr(i,j,k,m) = CS%IC_val(m)
              endif
          enddo ; enddo ; enddo
      endif ! restart
  enddo ! Tracer loop

  ! This needs to be changed if the units of tracer are changed above.
  flux_units = "years kg s-1" ;

    write(*,*) "Registering with diag manager"

  do m=1,CS%ntr
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for the restart file.
    call query_vardesc(CS%tr_desc(m), name, units=units, longname=longname, &
                       caller="initialize_ice_age_tracer")

    CS%id_tracer(m) = register_SIS_diag_field("ice_model", trim(name), CS%diag%axesTL, &
        CS%Time, trim(longname) , trim(units))
    CS%id_tr_adx(m) = register_SIS_diag_field("ice_model", trim(name)//"_adx", &
        CS%diag%axesCuc, CS%Time, trim(longname)//" advective zonal flux" , &
        trim(flux_units))
    CS%id_tr_ady(m) = register_SIS_diag_field("ice_model", trim(name)//"_ady", &
        CS%diag%axesCvc, CS%Time, trim(longname)//" advective meridional flux" , &
        trim(flux_units))

    if (CS%id_tr_adx(m) > 0) call safe_alloc_ptr(CS%tr_adx(m)%p,IscB,IecB,jsc,jec,CS%nlevels(m))
    if (CS%id_tr_ady(m) > 0) call safe_alloc_ptr(CS%tr_ady(m)%p,isc,iec,JscB,JecB,CS%nlevels(m))

  enddo

  write(*,*) "Ice age tracer initialized"

end subroutine initialize_ice_age_tracer

subroutine ice_age_tracer_column_physics(dt, G, IG,  mi, mi_old, CS)
  real,                               intent(in) :: dt
  type(sis_hor_grid_type),                intent(in) :: G
  type(ice_grid_type),                intent(in) :: IG
  type(ice_age_tracer_CS),          pointer    :: CS
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce), intent(in) :: mi, mi_old

! Arguments:
!  (in)      dt - The amount of time covered by this call, in s.
!  (in)      mi_old - Mass of ice at the beginning of the ice model timestep
!  (in)      G - The ocean model's grid structure.
!  (in)      IG - The ice model's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_ideal_age_tracer.

  real :: Isecs_per_year  ! The number of seconds in a year.
  real :: year            ! The time in years.
  real,dimension(SZI_(G),SZJ_(G)) :: vertsum_mi, vertsum_mi_old
  integer :: secs, days   ! Integer components of the time type.
  integer :: i, j, k, m
  integer :: isc, iec, jsc, jec, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; NkIce = IG%NkIce

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  Isecs_per_year = 1.0 / (365.0*86400.0)
  call get_time(CS%Time, secs, days)
  year = (86400.0*days + real(secs)) * Isecs_per_year

    vertsum_mi = 0.0
    vertsum_mi_old = 0.0;

    do j=jsc,jec; do i=isc,iec
        do k=1,NkIce
            vertsum_mi(i,j) = vertsum_mi(i,j) + mi(i,j,k)
            vertsum_mi_old(i,j) = vertsum_mi_old(i,j) + mi_old(i,j,k)
        enddo
    enddo; enddo

  ! Increment the age of the ice if it exists
  do m=1,CS%ntr ; if (CS%tracer_ages(m) .and. &
      (year>=CS%tracer_start_year(m))) then
      !$OMP parallel do default(none) shared(is,ie,js,je,CS,IG,dt,Isecs_per_year,m)
      do k=1,CS%nlevels(m) ; do j=jsc,jec ; do i=isc,iec

          if(vertsum_mi(i,j) > 0.0) then
              CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j)*dt*Isecs_per_year
          else
              CS%tr(i,j,k,m) = 0.0;
          endif
      enddo ; enddo ; enddo
  endif ; enddo



! If newly formed ice reduces the age, then apply the net sink term
  do m=1,CS%ntr ; if (CS%new_ice_is_sink(m) .and. &
                      (year>=CS%tracer_start_year(m))) then
!$OMP parallel do default(none) shared(isc,iec,jsc,jec,CS,nkIce,IG,dt,Isecs_per_year,m)
    do k=1,CS%nlevels(m) ; do j=jsc,jec ; do i=isc,iec

        ! ashao: Need to think about how and where the sink associated
        ! with newly formed sea ice gets appliec

    enddo ; enddo ; enddo
  endif ; enddo


! Apply masks if requested
  do m=1,CS%ntr
    if (CS%mask_tracers) then
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr_aux(:,:,:,m),CS%diag)
    else
      if (CS%id_tracer(m)>0) &
        call post_data(CS%id_tracer(m),CS%tr(:,:,:,m),CS%diag)
    endif
    if (CS%id_tr_adx(m)>0) &
      call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
    if (CS%id_tr_ady(m)>0) &
      call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
  enddo

end subroutine ice_age_tracer_column_physics

function ice_age_stock(h, stocks, G, IG, CS, names, units, stock_index)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in)    :: h
  real, dimension(:),                 intent(out)   :: stocks
  type(sis_hor_grid_type),              intent(in)    :: G
  type(ice_grid_type),              intent(in)    :: IG
  type(ice_age_tracer_CS),          pointer       :: CS
  character(len=*), dimension(:),     intent(out)   :: names
  character(len=*), dimension(:),     intent(out)   :: units
  integer, optional,                  intent(in)    :: stock_index
  integer                                           :: ice_age_stock
! This function calculates the mass-weighted integral of all tracer stocks,
! returning the number of stocks it has calculated.  If the stock_index
! is present, only the stock corresponding to that coded index is returned.

! Arguments: h - Layer thickness, in m or kg m-2.
!  (out)     stocks - the mass-weighted integrated amount of each tracer,
!                     in kg times concentration units.
!  (in)      G - The ocean's grid structure.
!  (in)      G - The ice model's grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 register_ideal_age_tracer.
!  (out)     names - the names of the stocks calculated.
!  (out)     units - the units of the stocks calculated.
!  (in,opt)  stock_index - the coded index of a specific stock being sought.
! Return value: the number of stocks calculated here.

  integer :: i, j, k, m
  integer :: isc, iec, jsc, jec, NkIce
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; NkIce = IG%NkIce

  ice_age_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="ice_age_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = 0.0
    do k=1,NkIce ; do j=jsc,jec ; do i=isc,iec
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                             (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo 
    stocks(m) = IG%H_to_kg_m2 * stocks(m)
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

    deallocate(CS)
  endif
end subroutine ice_age_end

end module ice_age_tracer
