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
    use SIS_diag_mediator, only : register_SIS_diag_field, &
        safe_alloc_ptr
    use SIS_diag_mediator, only: SIS_diag_ctrl, post_data=>post_SIS_data
    use SIS_tracer_registry, only : register_SIS_tracer, SIS_tracer_registry_type
    use SIS_hor_grid_mod, only : sis_hor_grid_type

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

    use ice_grid_mod, only : ice_grid_type

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
            tracer_start_year = 0.0  ! The year in which tracers start aging, or at which the
                             ! surface value equals young_val, in years.
        logical :: mask_tracers  ! If true, tracers are masked out in massless layers.
        logical :: tracers_may_reinit  ! If true, tracers may go through the
                                 ! initialization code if they are not found in the
                                 ! restart files.
        logical :: tracer_ages(NTR_MAX) ! Flag to if the tracer ages as a function of time
        logical :: new_ice_is_sink(NTR_MAX) ! Flag for whether the formation of new ice acts
                                            ! as a sink for age
        logical :: advect_vertical(NTR_MAX) ! Whether the tracer should be advected "vertically"
                                            ! to thicker ice
        logical :: uniform_vertical(NTR_MAX)! Whether the tracer should uniform across ice thickness
                                            ! categories if so, "mix" the tracer by setting it to
                                            ! the maximum age at the grid pont
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



contains

    logical function register_ice_age_tracer(G, IG, param_file, CS, diag, TrReg, &
        Ice_restart, restart_file)
        type(sis_hor_grid_type),  intent(in) :: G
        type(ice_grid_type),      intent(in) :: IG
        type(param_file_type),    intent(in) :: param_file
        type(ice_age_tracer_CS), pointer     :: CS
        type(SIS_diag_ctrl), target,   intent(inout) :: diag
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
        integer :: isc, iec, jsc, jec, m
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
        call get_param(param_file, mod, "DO_AREAL_ICE_AGE", CS%do_areal_ice_age, &
            "If true, use an ice age tracer that is set to 0 age \n"//&
            "when ice is first formed and ages at unit rate in the ice pack.", &
            default=.false.)
        call get_param(param_file, mod, "DO_MASS_ICE_AGE", CS%do_mass_ice_age, &
            "If true, use an ice age tracer that is set to 0 age \n"//&
            "when ice is first formed, ages at unit rate in the ice pack \n"// &
            "and any new ice added to the ice pack represents a decrease in age", &
            default=.false.)

        call get_param(param_file, mod, "ICE_AGE_IC_FILE", CS%IC_file, &
            "The file in which the age-tracer initial values can be \n"//&
            "found, or an empty string for internal initialization.", &
            default="")
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
            CS%tracer_ages(m) = .true.
            CS%new_ice_is_sink(m) = .false.
            CS%uniform_vertical(m) = .true.
            CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = -1.0
            CS%land_val(m) = 0.0
            CS%nlevels(m) = IG%CatIce
        endif
        if (CS%do_mass_ice_age) then
            CS%ntr = CS%ntr + 1 ; m = CS%ntr
            CS%tr_desc(m) = var_desc("ice_age_mass", "years", "Mass Ice Age Tracer", caller=mod)
            CS%tracer_ages(m) = .true.
            CS%new_ice_is_sink(m) = .true.
            CS%uniform_vertical(m) = .false.
            CS%IC_val(m) = 0.0 ; CS%young_val(m) = 0.0 ; CS%tracer_start_year(m) = -1.0
            CS%land_val(m) = 0.0
            CS%nlevels(m) = IG%CatIce
        endif

        allocate(CS%tr(SZI_(G), SZJ_(G),IG%CatIce,CS%ntr)) ; CS%tr(:,:,:,:) = 0.0

        CS%diag => diag
        do m=1,CS%ntr

            call query_vardesc(CS%tr_desc(m), name=var_name, &
                caller="register_ice_age_tracer")

            ! Register the tracer for the restart file.
            CS%id_tracer(m) = register_restart_field(Ice_restart, restart_file, var_name, &
                CS%tr(:,:,:,m), domain=G%domain%mpp_domain, &
                mandatory=.false.)

            ! Register the tracer for horizontal advection & diffusion.
            call register_SIS_tracer(CS%tr(:,:,:,m), G, IG, 1, var_name, param_file, &
                TrReg, snow_tracer=.false.)
                             

        enddo ! do m=1, CS%ntr

        CS%TrReg => TrReg
        register_ice_age_tracer = .true.

    end function register_ice_age_tracer


    subroutine initialize_ice_age_tracer( day, G, IG, CS )
        type(time_type), target,               intent(in) :: day
        type(sis_hor_grid_type),              intent(in) :: G
        type(ice_grid_type),              intent(in) :: IG
        type(ice_age_tracer_CS),          pointer    :: CS
        real, parameter                              :: missing = -1.0e34
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


        if (.not.associated(CS)) return
        if (CS%ntr < 1) return

        isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
        IscB = G%IscB ; IecB = G%IecB ; JscB = G%JscB ; JecB = G%JecB

        CS%Time => day
        CS%tr(:,:,:,:) = 0.0
        do m=1,CS%ntr
            do k=1,CS%nlevels(m) ; do j=jsc,jec ; do i=isc,iec
                if (G%mask2dT(i,j) < 0.5) then
                    CS%tr(i,j,k,m) = CS%land_val(m)
                else
                    CS%tr(i,j,k,m) = CS%IC_val(m)
                endif
            enddo ; enddo ; enddo
        enddo ! Tracer loop

        ! This needs to be changed if the units of tracer are changed above.
        flux_units = "years kg s-1"

        do m=1,CS%ntr

            call query_vardesc(CS%tr_desc(m), name, units=units, longname=longname, &
                caller="initialize_ice_age_tracer")

            CS%id_tracer(m) = register_SIS_diag_field("ice_model", trim(name), CS%diag%axesTC, &
                CS%Time, trim(longname) , trim(units),missing_value = missing)
            CS%id_tr_adx(m) = register_SIS_diag_field("ice_model", trim(name)//"_adx", &
                CS%diag%axesCuc, CS%Time, trim(longname)//" advective zonal flux" , &
                trim(flux_units))
            CS%id_tr_ady(m) = register_SIS_diag_field("ice_model", trim(name)//"_ady", &
                CS%diag%axesCvc, CS%Time, trim(longname)//" advective meridional flux" , &
                trim(flux_units))

            if (CS%id_tr_adx(m) > 0) call safe_alloc_ptr(CS%tr_adx(m)%p,IscB,IecB,jsc,jec,CS%nlevels(m))
            if (CS%id_tr_ady(m) > 0) call safe_alloc_ptr(CS%tr_ady(m)%p,isc,iec,JscB,JecB,CS%nlevels(m))

        enddo

    end subroutine initialize_ice_age_tracer

    subroutine ice_age_tracer_column_physics(dt, G, IG, CS,  mi, mi_old)
        real,                               intent(in) :: dt
        type(sis_hor_grid_type),                intent(in) :: G
        type(ice_grid_type),                intent(in) :: IG
        type(ice_age_tracer_CS),          pointer    :: CS
        real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)), optional :: mi, mi_old

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
        real,dimension(SZI_(G),SZJ_(G)) :: vertsum_mi, vertsum_mi_old
        integer :: secs, days   ! Integer components of the time type.
        integer :: i, j, k, m
        integer :: isc, iec, jsc, jec
        isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ;

        if (.not.associated(CS)) return
        if (CS%ntr < 1) return

        mi_min = 1.0e-7*IG%kg_m2_to_H

        Isecs_per_year = 1.0 / (365.0*86400.0)
        dt_year = dt * Isecs_per_year
        min_age = dt_year * 0.1

        call get_time(CS%Time, secs, days)
        year = (86400.0*days + real(secs)) * Isecs_per_year


        vertsum_mi(:,:) = 0.0
        vertsum_mi_old(:,:) = 0.0

        do j=jsc,jec; do i=isc,iec
            do k=1,IG%CatIce
                vertsum_mi(i,j) = vertsum_mi(i,j) + mi(i,j,k)
                vertsum_mi_old(i,j) = vertsum_mi_old(i,j) + mi_old(i,j,k)
            enddo
        enddo; enddo

        ! Increment the age of the ice if it exists
        do m=1,CS%ntr ; if (CS%tracer_ages(m) .and. &
            (year>=CS%tracer_start_year(m))) then
            do k=1,CS%nlevels(m) ; do j=jsc,jec ; do i=isc,iec

                if(CS%tr(i,j,k,m)<min_age) CS%tr(i,j,k,m) = 0.0

                if(mi(i,j,k)>mi_min) then
                    CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + dt_year
                else
                    CS%tr(i,j,k,m) = 0.0
                endif

            enddo ; enddo ; enddo

        endif ; enddo

        ! If newly formed ice reduces the age, then apply the net sink term
        do m=1,CS%ntr ; if (CS%new_ice_is_sink(m)) then
            do k=1,CS%nlevels(m) ; do j=jsc,jec ; do i=isc,iec
                ! @ashao: Need to think about how and where the sink associated
                ! with newly formed sea ice gets applied. For now, apply to
                ! the entire ice pack assuming that new ice is formed equally
                ! on every ice thickness category

                if(mi(i,j,k)>mi_old(i,j,k)) then
                    CS%tr(i,j,k,m) = CS%tr(i,j,k,m) * (mi_old(i,j,k) / mi(i,j,k))
                endif

            enddo ; enddo ; enddo
        endif ; enddo

        ! If the tracer should be uniform, set age at every grid point to the maximum
        do m=1,CS%ntr ; if (CS%uniform_vertical(m)) then
            do j=jsc,jec ; do i=isc,iec
                if(vertsum_mi(i,j) > 0.0) then
                    max_age = maxval(CS%tr(i,j,:,m))
                    do k=1,CS%nlevels(m)
                        CS%tr(i,j,k,m) = max_age
                    enddo

                endif
                print *, "Max ice age ", max_age
            enddo ; enddo
        endif ; enddo

        do m=1,CS%ntr
            if (CS%id_tracer(m)>0) &
                call post_data(CS%id_tracer(m),CS%tr(:,:,:,m),CS%diag)
            if (CS%id_tr_adx(m)>0) &
                call post_data(CS%id_tr_adx(m),CS%tr_adx(m)%p(:,:,:),CS%diag)
            if (CS%id_tr_ady(m)>0) &
                call post_data(CS%id_tr_ady(m),CS%tr_ady(m)%p(:,:,:),CS%diag)
        enddo

    end subroutine ice_age_tracer_column_physics

    subroutine ice_age_stock(nstocks, stocks, names, units, G, IG, CS, mi)
        integer                                           :: nstocks
        real, dimension(:)                                :: stocks
        character(len=*), dimension(:)                    :: names
        character(len=*), dimension(:)                    :: units
        type(sis_hor_grid_type),               intent(in) :: G
        type(ice_grid_type),                   intent(in) :: IG
        type(ice_age_tracer_CS),               pointer    :: CS
        real, dimension(SZI_(G),SZJ_(G),SZCAT_(IG)),intent(in) :: mi

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

        integer :: i, j, k, m
        integer :: isc, iec, jsc, jec
        isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

        if (.not.associated(CS)) return
        if (CS%ntr < 1) return

        do m=1,CS%ntr
            nstocks = nstocks + 1
            call query_vardesc(CS%tr_desc(m), name=names(nstocks), units=units(nstocks), caller="ice_age_stock")
            units(nstocks) = trim(units(m))//" kg"
            stocks(nstocks) = 0.0
            do k=1,IG%CatIce ; do j=jsc,jec ; do i=isc,iec
                stocks(nstocks) = stocks(nstocks) + CS%tr(i,j,k,m) * &
                    (G%mask2dT(i,j) * G%areaT(i,j) * mi(i,j,k))
            enddo ; enddo ; enddo
!            stocks(nstocks) = IG%H_to_kg_m2 * stocks(nstocks)
        enddo


    end subroutine ice_age_stock

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
