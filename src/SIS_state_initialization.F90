!> This module is used to initalize the sea ice state for SIS2
module SIS_state_initialization

! This file is a part of SIS2. See LICENSE.md for the license.

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!   This module has code that initializes the sea ice state at the start of a  !
! run.  If a run segment is initialzed from a restart file, some of these      !
! routines have options that just read and log their input parameters.         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use ice_grid,          only : ice_grid_type
use ice_type_mod,      only : ice_data_type, dealloc_ice_arrays
use ice_type_mod,      only : ice_type_slow_reg_restarts
use MOM_data_override, only : data_override, data_override_init, data_override_unset_domains
use MOM_domains,       only : MOM_domain_type
use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_error_handler, only : callTree_enter, callTree_leave, callTree_waypoint
use MOM_file_parser,   only : get_param, log_param, log_version, read_param, param_file_type
use MOM_hor_index,     only : hor_index_type, hor_index_init
use MOM_io,            only : file_exists, MOM_read_data, slasher
use MOM_time_manager,  only : time_type, time_type_to_real, real_to_time
use MOM_unit_scaling,  only : unit_scale_type
use SIS_restart,       only : restore_SIS_state, query_initialized=>query_inited
use SIS_restart,       only : register_restart_field, only_read_from_restarts
use SIS_get_input,     only : directories
use SIS_types,         only : ice_state_type
use SIS_hor_grid,      only : SIS_hor_grid_type, set_hor_grid, SIS_hor_grid_end
use SIS2_ice_thm,      only : get_SIS2_thermo_coefs, enth_from_TS, Temp_from_En_S, T_freeze, ice_thermo_type

implicit none ; private

#include <SIS2_memory.h>

public :: ice_state_mass_init, ice_state_thermo_init, initialize_ice_categories
public :: read_archaic_thermo_restarts

contains


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_ice_categories sets the bounds of the ice thickness categories.
subroutine initialize_ice_categories(IG, Rho_ice, US, PF, hLim_vals)
  type(ice_grid_type),          intent(inout) :: IG  !< The sea-ice specific grid type
  real,                         intent(in)    :: Rho_ice !< The nominal ice density [R ~> kg m-3].
  type(unit_scale_type),        intent(in)    :: US  !< A structure with unit conversion factors
  type(param_file_type),        intent(in)    :: PF  !< A structure to parse for run-time parameters
  real, dimension(:), optional, intent(in)    :: hLim_vals !< The ice category thickness limits [m].

  ! Initialize IG%cat_thick_lim and IG%mH_cat_bound here.
  !  ###This subroutine should be extended to add more options.

  real :: hlim_dflt(8) = (/ 1.0e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5 /) ! lower thickness limits 1...CatIce
  integer :: k, CatIce, list_size

  CatIce = IG%CatIce
  list_size = -1
  if (present(hLim_vals)) then ; if (size(hLim_vals(:)) > 1) then
    list_size = size(hlim_vals(:))
    do k=1,min(CatIce+1,list_size) ; IG%cat_thick_lim(k) = hlim_vals(k) ; enddo
  endif ; endif
  if (list_size < 2) then  ! Use the default categories.
    list_size = size(hlim_dflt(:))
    do k=1,min(CatIce+1,list_size) ; IG%cat_thick_lim(k) = hlim_dflt(k) ; enddo
  endif

  if ((CatIce+1 > list_size) .and. (list_size > 1)) then
    do k=list_size+1, CatIce+1
      IG%cat_thick_lim(k) =  2.0*IG%cat_thick_lim(k-1) - IG%cat_thick_lim(k-2)
    enddo
  endif

  do k=1,IG%CatIce+1
    IG%mH_cat_bound(k) = IG%cat_thick_lim(k)*US%m_to_Z * Rho_ice
  enddo
end subroutine initialize_ice_categories


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_state_mass_init initializes the sea ice concentration and ice and snow mass per unit area.
!! They will be distributed to the correct thickness category later in the initialization.
subroutine ice_state_mass_init(IST, Ice, G, IG, US, PF, init_Time, just_read_params)

  type(ice_state_type),    intent(inout) :: IST  !< The sea ice state type being modified
  type(ice_data_type),     intent(inout) :: Ice  !< The ice data type that is being initialized.
  type(SIS_hor_grid_type), intent(in)    :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)    :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)    :: US   !< A structure with unit conversion factors
  type(param_file_type),   intent(in)    :: PF   !< A structure to parse for run-time parameters
  type(time_type),         intent(in)    :: init_Time  !< The initialization time of the run segment
  logical,       optional, intent(in)    :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing ice properties.

  real :: rho_ice         ! The nominal density of sea ice [R ~> kg m-3].
  real :: rho_snow        ! The nominal density of snow [R ~> kg m-3].
  real :: mH_ice_uniform  ! A uniform initial sea ice mass per unit area [R Z ~> kg m-2]
  real :: mH_snow_uniform ! A uniform initial snow mass per unit area [R Z ~> kg m-2]
  real, dimension(SZI_(G),SZJ_(G)) :: h_input    ! Temporary ice thickness array [m]
  real, dimension(SZI_(G),SZJ_(G)) :: total_conc ! Summed ice concentration [nondim]
# include "version_variable.h"
  character(len=40)  :: mdl = "SIS_state_initialization" ! This module's name.
  character(len=200) :: conc_config, hIce_config, hSnow_config
  logical :: just_read, any_data_override
  integer :: i, j, k, n
  integer :: isc, iec, jsc, jec, CatIce

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec ; CatIce = IG%CatIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(PF, mdl, "RHO_ICE", Rho_ice, &
                 "The nominal density of sea ice as used by SIS.", &
                 units="kg m-3", default=905.0, scale=US%kg_m3_to_R, do_not_log=.true.)
  call get_param(PF, mdl, "RHO_SNOW", Rho_snow, &
                 "The nominal density of snow as used by SIS.", &
                 units="kg m-3", default=330.0, scale=US%kg_m3_to_R, do_not_log=.true.)

  call get_param(PF, mdl, "CONCENTRATION_INIT_CONFIG", conc_config, &
            "A string that determines how the initial total sea ice concentration "//&
            "is initialized for a new run: \n"//&
            "    file - read sea ice concentrations from a specified file \n"//&
            "    data_override - use the data_override capability or zero everywhere \n"//&
            "    zero - there is no sea ice anywhere \n"//&
            "    latitudes - initial sea ice concentration is a function of latitude.", &
            default='data_override', do_not_log=just_read)

  call get_param(PF, mdl, "ICE_THICKNESS_INIT_CONFIG", hIce_config, &
            "A string that determines how the initial sea ice thickness "//&
            "is initialized for a new run: \n"//&
            "    file - read sea ice thickesses from a specified file \n"//&
            "    data_override - use the data_override capability or zero everywhere \n"//&
            "    uniform - sea ice has uniform thickness where the concentration is nonzero.", &
            default='data_override', do_not_log=just_read)

  call get_param(PF, mdl, "SNOW_THICKNESS_INIT_CONFIG", hSnow_config, &
            "A string that determines how the initial total snow thickness "//&
            "is initialized for a new run: \n"//&
            "    file - read sea ice concentrations from a specified file \n"//&
            "    data_override - use the data_override capability or zero everywhere \n"//&
            "    uniform - snow has uniform thickness where the concentration is nonzero.", &
            default='data_override', do_not_log=just_read)

  any_data_override = (((trim(conc_config)=="data_override") .or. &
                        (trim(hIce_config)=="data_override") .or. &
                        (trim(hSnow_config)=="data_override")) .and. .not.just_read)

  if (.not.just_read) then
    IST%part_size(:,:,:) = 0.0
    IST%part_size(:,:,0) = 1.0
  endif

  if (any_data_override) then
    call data_override_init()
    call data_override_unset_domains(unset_Ice=.true., must_be_set=.false.)
    call data_override_init(Ice_domain_in=Ice%slow_domain_NH)
  endif

  select case (trim(conc_config))
    case ("data_override")
      if (.not.just_read) &
        call data_override('ICE', 'sic_obs', IST%part_size(isc:iec,jsc:jec,1), init_Time)
    case ("file")
      call initialize_concentration_from_file(IST%part_size, G, IG, US, PF, just_read_params=just_read)
    case ("zero")
      if (.not.just_read) &
        IST%part_size(:,:,1:CatIce) = 0.0
    case ("latitudes")
      call initialize_concentration_from_latitudes(IST%part_size, G, IG, US, PF, just_read_params=just_read)
  end select

  select case (trim(hIce_config))
    case ("data_override")
      if (.not.just_read) then
        h_input(:,:) = 0.0
        call data_override('ICE', 'sit_obs', h_input(isc:iec,jsc:jec), init_Time)
        do j=jsc,jec ; do i=isc,iec
          IST%mH_ice(i,j,1) = h_input(i,j)*US%m_to_Z * Rho_ice
        enddo ; enddo
      endif
    case ("file")
      call initialize_thickness_from_file(IST%mH_ice, G, IG, US, PF, just_read_params=just_read)
    case ("uniform")
      call get_param(PF, mdl, "ICE_INIT_MASS", mH_ice_uniform, &
                 "A uniform initial sea ice mass per unit area where there is sea ice.", &
                 units="kg m-2", default=0.0, scale=US%kg_m3_to_R*US%m_to_Z, do_not_log=just_read)
      if (.not.just_read) then ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
        if (IST%part_size(i,j,k) > 0.0) IST%mH_ice(i,j,k) = mH_ice_uniform
      enddo ; enddo ; enddo ; endif
  end select

  select case (trim(hSnow_config))
    case ("data_override")
      if (.not.just_read) then
        h_input(:,:) = 0.0
        call data_override('ICE', 'si_snow_obs', h_input(isc:iec,jsc:jec), init_Time)
        do j=jsc,jec ; do i=isc,iec
          IST%mH_snow(i,j,1) = h_input(i,j)*US%m_to_Z * Rho_ice
        enddo ; enddo
      endif
    case ("file")
      call initialize_snow_thick_from_file(IST%mH_snow, G, IG, US, PF, just_read_params=just_read)
    case ("uniform")
      call get_param(PF, mdl, "SNOW_INIT_MASS", mH_snow_uniform, &
                 "A uniform initial snow mass per unit area where there is sea ice.", &
                 units="kg m-2", default=0.0, scale=US%kg_m3_to_R*US%m_to_Z, do_not_log=just_read)
      if (.not.just_read) then ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
        if (IST%part_size(i,j,k) > 0.0) IST%mH_snow(i,j,k) = mH_snow_uniform
      enddo ; enddo ; enddo ; endif
  end select

  if (just_read) return

  total_conc(:,:) = 0.0
  do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
    total_conc(i,j) = total_conc(i,j) + IST%part_size(i,j,k)
    if (IST%part_size(i,j,k) == 0.0) then
      IST%mH_ice(i,j,k) = IG%mH_cat_bound(k+1)
      IST%mH_snow(i,j,k) = 0.0
    endif
  enddo ; enddo ; enddo

  do j=jsc,jec ; do i=isc,iec
    IST%part_size(i,j,0) = 1.0 - total_conc(i,j)
  enddo ; enddo

  if (any_data_override) call data_override_unset_domains(unset_Ice=.true.)

end subroutine ice_state_mass_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> ice_state_thermo_init initializes the sea ice and snow enthalpies and salinities.  The ice
!! concentration and ice and snow mass must already be set, but they might not be in the correct
!! thickness category yet.  That redistribution will occur later in the initialization.
subroutine ice_state_thermo_init(IST, Ice, G, IG, US, PF, init_Time, just_read_params)
  type(ice_state_type),    intent(inout) :: IST  !< The sea ice state type being modified
  type(ice_data_type),     intent(inout) :: Ice  !< The ice data type that is being initialized.
  type(SIS_hor_grid_type), intent(in)    :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)    :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)    :: US   !< A structure with unit conversion factors
  type(param_file_type),   intent(in)    :: PF   !< A structure to parse for run-time parameters
   type(time_type),        intent(in)    :: init_Time      !< The initialization time of the run segment
  logical,       optional, intent(in)    :: just_read_params !< If true only read parameters, but
                                                 !! do not log them or change the ice state.
!
  ! Local variables
  real :: enth_spec_snow, enth_spec_ice  ! Specified enthalpy of snow and ice [Q ~> J kg-1]
  real :: ice_bulk_salin ! The globally constant sea ice bulk salinity [gSalt kg-1] = [ppt]
                         ! that is used to calculate the ocean salt flux.
  real :: ice_rel_salin  ! The initial bulk salinity of sea-ice relative to the
                         ! salinity of the water from which it formed [nondim].
  real :: ice_salin_IC   ! The initial ice bulk salinity [gSalt kg-1] = [ppt]
  real :: ice_temp_IC    ! The initial ice temperature [degC]
  real :: ice_rel_temp_IC ! The initial ice temperature relative to the freezing point [degC]
  real :: S_col(IG%NkIce) ! Specified ice column salinity used for ice thermodynamics [gSalt kg-1]
  real, dimension(SZI_(G),SZJ_(G)) :: salin_input  ! Temporary ice salinity [gSalt kg-1]
  real, dimension(SZI_(G),SZJ_(G)) :: temp_input   ! Temporary ice temperature [degC]
# include "version_variable.h"
  character(len=40)  :: mdl = "SIS_state_initialization" ! This module's name.
  character(len=200) :: salin_config, enth_ice_config, enth_snow_config
  logical :: just_read    ! If true, just read parameters but set nothing.
  logical :: any_data_override
  logical :: spec_thermo_sal ! If true used a specified ice profile for thermodynamics instead
                             ! of the actual ice salinity profile.
  integer :: i, j, k, n
  integer :: isc, iec, jsc, jec, CatIce, NkIce

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  CatIce = IG%CatIce ; NkIce = IG%NkIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(PF, mdl, "ICE_SALINITY_INIT_CONFIG", salin_config, &
            "A string that determines how the sea ice salinity is initialized for a new run: \n"//&
            " \t uniform -  Use a constant ice salinity initial condition \n"//&
            " \t file - Read sea ice salinities from a specified file \n"//&
            " \t data_override - use the data_override capability or zero everywhere.", &
            default='uniform', do_not_log=just_read)

  call get_param(PF, mdl, "ICE_ENTHALPY_INIT_CONFIG", enth_ice_config, &
            "A string that determines how the sea ice enthalpy is initialized for a new run: \n"//&
            " \t uniform_temp - Use a constant ice temperature initial condition \n"//&
            " \t relative_temp - Use an ice temperature initial condition with a \n"//&
            " \t\t specified depression below the bulk ice freezing point \n"//&
            " \t file - Read sea ice temperatures or enthalpies from a specified file \n"//&
            " \t data_override - use the data_override capability or freezing enthalpy everywhere.", &
            default='uniform_temp', do_not_log=just_read)

  call get_param(PF, mdl, "SNOW_ENTHALPY_INIT_CONFIG", enth_snow_config, &
            "A string that determines how the snow enthalpy is initialized for a new run: \n"//&
            " \t uniform_temp - Use a constant ice temperature initial condition \n"//&
            " \t relative_temp - Use an ice temperature initial condition with a \n"//&
            " \t\t specified depression below the bulk ice freezing point \n"//&
            " \t file - Read sea ice temperatures or enthalpies from a specified file \n"//&
            " \t data_override - use the data_override capability or freezing enthalpy everywhere.", &
            default='uniform_temp', do_not_log=just_read)

  call get_param(PF, mdl, "ICE_BULK_SALINITY", ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "g/kg", &
                 default=4.0, do_not_log=.true.)
  call get_param(PF, mdl, "ICE_RELATIVE_SALINITY", ice_rel_salin, &
                 "The initial salinity of sea ice as a fraction of the "//&
                 "salinity of the seawater from which it formed.", &
                 units = "nondim", default=0.0, do_not_log=.true.)
  if ((ice_bulk_salin < 0.0) .or. (ice_rel_salin > 0.0)) ice_bulk_salin = 0.0

  S_col(:) = 0.0
  call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, spec_thermo_salin=spec_thermo_sal)

  any_data_override = (((trim(salin_config)=="data_override") .or. &
                        (trim(enth_ice_config)=="data_override") .or. &
                        (trim(enth_snow_config)=="data_override")) .and. .not.just_read)

  if (any_data_override) then
    call data_override_init()
    call data_override_unset_domains(unset_Ice=.true., must_be_set=.false.)
    call data_override_init(Ice_domain_in=Ice%slow_domain_NH)
  endif

  ! Initialize the sea ice salinity.  All ice layers and categories are set here.
  select case (trim(salin_config))
    case ("uniform")
      call get_param(PF, mdl, "ICE_SALINITY_IC", ice_salin_IC, &
                 "The uniform sea ice salinity used for the initial condition", &
                 units="g kg-1", default=ice_bulk_salin, do_not_log=just_read)
      if (.not.just_read) then ; do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
        !### This should not change answers but does:  if (IST%part_size(i,j,k) > 0.0)
        IST%sal_ice(i,j,k,n) = ice_salin_IC
      enddo ; enddo ; enddo ; enddo ; endif
    case ("data_override")
      if (.not.just_read) then
        salin_input(:,:) = 0.0
        call data_override('ICE', 'si_salin_obs', salin_input(isc:iec,jsc:jec), init_Time)
        do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%sal_ice(i,j,k,n) = salin_input(i,j)
        enddo ; enddo ; enddo ; enddo
      endif
    case ("file")
      call initialize_salinity_from_file(IST%sal_ice, G, IG, US, PF, just_read_params=just_read)
  end select

  ! Initialize the sea ice enthalpy.  All ice layers and categories are set here.
  select case (trim(enth_ice_config))
    case ("uniform_temp")
      call get_param(PF, mdl, "ICE_TEMPERATURE_IC", ice_temp_IC, &
                 "The uniform sea ice and snow temperature used for the initial condition", &
                 units="degC", default=-4.0, do_not_log=just_read)
      if (spec_thermo_sal .and. (.not.just_read)) then
        do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,n) = Enth_from_TS(ice_temp_IC, S_col(n), IST%ITV)
        enddo ; enddo ; enddo ; enddo
      elseif (.not.just_read) then
        do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,n) = Enth_from_TS(ice_temp_IC, IST%sal_ice(i,j,k,n), IST%ITV)
        enddo ; enddo ; enddo ; enddo
      endif
    case ("relative_temp")
      call get_param(PF, mdl, "ICE_RELATIVE_TEMP_IC", ice_rel_temp_IC, &
                 "The sea ice and snow temperature relative to the local bulk freezing point "//&
                 "used for the initial condition", units="degC", default=-4.0, do_not_log=just_read)
      if (.not.just_read) then ; do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
        IST%enth_ice(i,j,k,n) = Enth_from_TS(ice_rel_temp_IC, 0.0, IST%ITV)
      enddo ; enddo ; enddo ; enddo ; endif
      if (spec_thermo_sal .and. (.not.just_read)) then
        do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,n) = Enth_from_TS(T_Freeze(S_col(n), IST%ITV) + ice_rel_temp_IC, &
                                                S_col(n), IST%ITV)
        enddo ; enddo ; enddo ; enddo
      elseif (.not.just_read) then
        do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,n) = Enth_from_TS(T_Freeze(IST%sal_ice(i,j,k,n), IST%ITV) + &
                                                ice_rel_temp_IC, IST%sal_ice(i,j,k,n), IST%ITV)
        enddo ; enddo ; enddo ; enddo
      endif
    case ("data_override")
      if (.not.just_read) then
        temp_input(:,:) = 0.0
        call data_override('ICE', 'si_temp_obs', temp_input(isc:iec,jsc:jec), init_Time)
        if (spec_thermo_sal .and. (.not.just_read)) then
          do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
            IST%enth_ice(i,j,k,n) = Enth_from_TS(temp_input(i,j), S_col(n), IST%ITV)
          enddo ; enddo ; enddo ; enddo
        elseif (.not.just_read) then
          do n=1,NkIce ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
            IST%enth_ice(i,j,k,n) = Enth_from_TS(temp_input(i,j), IST%sal_ice(i,j,k,n), IST%ITV)
          enddo ; enddo ; enddo ; enddo
        endif
      endif
    case ("file")
      call initialize_ice_enthalpy_from_file(IST%enth_ice, IST%sal_ice, G, IG, US, &
                                             IST%ITV, PF, just_read_params=just_read)
  end select


  ! Initialize the snow enthalpy.  All ice layers and categories are set here.
  select case (trim(enth_snow_config))
    case ("uniform_temp")
      call get_param(PF, mdl, "ICE_TEMPERATURE_IC", ice_temp_IC, &
                 "The uniform sea ice and snow temperature used for the initial condition", &
                 units="degC", default=-4.0, do_not_log=just_read)
      if (.not.just_read) then ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
        IST%enth_snow(i,j,k,1) = Enth_from_TS(ice_temp_IC, 0.0, IST%ITV)
      enddo ; enddo ; enddo ; endif
    case ("relative_temp")
      call get_param(PF, mdl, "ICE_RELATIVE_TEMP_IC", ice_rel_temp_IC, &
                 "The sea ice and snow temperature relative to the local bulk freezing point "//&
                 "used for the initial condition", units="degC", default=-4.0, do_not_log=just_read)
      if (.not.just_read) then ; do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
        IST%enth_snow(i,j,k,1) = Enth_from_TS(ice_rel_temp_IC, 0.0, IST%ITV)
      enddo ; enddo ; enddo ; endif
    case ("data_override")
      if (.not.just_read) then
        temp_input(:,:) = 0.0
        call data_override('ICE', 'si_temp_obs', temp_input(isc:iec,jsc:jec), init_Time)
        do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_snow(i,j,k,1) = Enth_from_TS(temp_input(i,j), 0.0, IST%ITV)
        enddo ; enddo ; enddo
      endif
    case ("file")
      call initialize_snow_enthalpy_from_file(IST%enth_snow, G, IG, US, IST%ITV, PF, just_read_params=just_read)
  end select

  if (just_read) return

  if (any_data_override) call data_override_unset_domains(unset_Ice=.true.)

end subroutine ice_state_thermo_init


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_concentration_from_file reads the ice concentration from a 2-d file.  This concentration
!! is placed in category 1, and will be distributed to the correct thickness category later.
subroutine initialize_concentration_from_file(part_size, G, IG, US, PF, just_read_params)
  type(SIS_hor_grid_type), intent(in)  :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)  :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),0:IG%CatIce), &
                           intent(out) :: part_size !< The ice concentration that is being initialized [nondim]
  type(param_file_type),   intent(in)  :: PF   !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing concentrations.

  ! Local variables
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "initialize_concentration_from_file" ! This subroutine's name.
  character(len=64)  :: conc_var ! Ice concentration variable name in files
  character(len=200) :: filename, concentration_file, inputdir ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, CatIce

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; CatIce = IG%CatIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), SIS_state_initialization.F90")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=just_read)
  call get_param(PF, mdl, "ICE_CONCENTRATION_FILE", concentration_file, &
                 "The name of the sea ice concentration file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)

  filename = trim(slasher(inputdir))//trim(concentration_file)
  if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/ICE_CONCENTRATION_FILE", filename)
  call get_param(PF, mdl, "ICE_CONCENTRATION_IC_VAR", conc_var, &
                 "The initial condition variable for ice mass per unit area.", &
                 default="conc_ice", do_not_log=just_read)
  if (just_read) return ! All run-time parameters have been read, so return.

  if (.not.file_exists(filename, G%Domain)) call SIS_error(FATAL, &
         " initialize_concentration_from_file: Unable to open "//trim(filename))

  do k=2,CatIce ; part_size(:,:,k) = 0.0 ; enddo
  call MOM_read_data(filename, conc_var, part_size(:,:,1), G%Domain, scale=US%kg_m3_to_R*US%m_to_Z)
  do j=js,je ; do i=is,ie ; part_size(i,j,0) = 1.0 - part_size(i,j,1) ; enddo ; enddo

  call callTree_leave(trim(mdl)//'()')

end subroutine initialize_concentration_from_file


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_concentration_from_file reads the ice concentration based on latitude.  This concentration
!! is placed in category 1, and will be distributed to the correct thickness category later.
subroutine initialize_concentration_from_latitudes(part_size, G, IG, US, PF, just_read_params)
  type(SIS_hor_grid_type), intent(in)  :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)  :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),0:IG%CatIce), &
                           intent(out) :: part_size !< The ice concentration that is being initialized [nondim]
  type(param_file_type),   intent(in)  :: PF   !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing concentrations.

  ! Local variables
  real :: Arctic_ice_edge    ! The southern latitude of Arctic ice in an initial condition [degrees of latitude]
  real :: Antarctic_ice_edge ! The nouthern latitude of Antarctic ice in an initial condition [degrees of latitude]
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "initialize_concentration_from_latitudes" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, CatIce

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; CatIce = IG%CatIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), SIS_state_initialization.F90")

  call get_param(PF, mdl, "ARCTIC_ICE_EDGE_IC", Arctic_ice_edge, &
                 "The southern latitude of Arctic ice in an initial condition.", &
                 default=91.0, units="degrees of latitude", do_not_log=just_read)
  call get_param(PF, mdl, "ANTARCTIC_ICE_EDGE_IC", Antarctic_ice_edge, &
                 "The northern latitude of Antarctic ice in an initial condition.", &
                 default=-91.0, units="degrees of latitude", do_not_log=just_read)
  if (just_read) return ! All run-time parameters have been read, so return.

  do j=js,je ; do i=is,ie
    part_size(i,j,1) = 0.0
    if ((G%geolatT(i,j) > Arctic_ice_edge) .or. (G%geolatT(i,j) < Antarctic_ice_edge)) &
      part_size(i,j,1) = 1.0
    part_size(i,j,0) = 1.0 - part_size(i,j,1)
  enddo ; enddo
  do k=2,CatIce ; part_size(:,:,k) = 0.0 ; enddo

  call callTree_leave(trim(mdl)//'()')

end subroutine initialize_concentration_from_latitudes


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_thickness_from_file reads the ice thickness from a 2-d file.  This thickness is
!! placed in category 1, and will be distributed to the correct thickness category later.
subroutine initialize_thickness_from_file(mH_ice, G, IG, US, PF, just_read_params)
  type(SIS_hor_grid_type), intent(in)  :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)  :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce), &
                           intent(out) :: mH_ice !< The ice thickness that is being initialized [R Z ~> kg m-2]
  type(param_file_type),   intent(in)  :: PF   !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing ice thickness.

  ! Local variables
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "initialize_thickness_from_file" ! This subroutine's name.
  character(len=64)  :: thick_var ! Ice thickness variable name in files
  character(len=200) :: filename, thickness_file, inputdir ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, CatIce

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; CatIce = IG%CatIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), SIS_state_initialization.F90")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=just_read)
  call get_param(PF, mdl, "ICE_THICKNESS_FILE", thickness_file, &
                 "The name of the sea ice thickness file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)

  filename = trim(slasher(inputdir))//trim(thickness_file)
  if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/ICE_THICKNESS_FILE", filename)
  call get_param(PF, mdl, "ICE_THICKNESS_IC_VAR", thick_var, &
                 "The initial condition variable for ice mass per unit area.", &
                 default="mH_ice", do_not_log=just_read)
  if (just_read) return ! All run-time parameters have been read, so return.

  if (.not.file_exists(filename, G%Domain)) call SIS_error(FATAL, &
         " initialize_thickness_from_file: Unable to open "//trim(filename))

  do k=2,CatIce ; mH_ice(:,:,k) = 0.0 ; enddo
  call MOM_read_data(filename, thick_var, mH_ice(:,:,1), G%Domain, scale=US%kg_m3_to_R*US%m_to_Z)

  call callTree_leave(trim(mdl)//'()')

end subroutine initialize_thickness_from_file


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_snow_thick_from_file reads the ice thickness from a 2-d file.  This thickness is
!! placed atop category 1, and will be distributed to the correct thickness category later.
subroutine initialize_snow_thick_from_file(mH_snow, G, IG, US, PF, just_read_params)
  type(SIS_hor_grid_type), intent(in)  :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)  :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce), &
                           intent(out) :: mH_snow !< The snow thickness that is being initialized [R Z ~> kg m-2]
  type(param_file_type),   intent(in)  :: PF   !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing snow thickness.

  ! Local variables
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40)  :: mdl = "initialize_snow_thick_from_file" ! This subroutine's name.
  character(len=64)  :: thick_var ! Snow thickness variable name in files
  character(len=200) :: filename, thickness_file, inputdir ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, CatIce

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; CatIce = IG%CatIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), SIS_state_initialization.F90")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=just_read)
  call get_param(PF, mdl, "SNOW_THICKNESS_FILE", thickness_file, &
                 "The name of the snow thickness on sea ice file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)

  filename = trim(slasher(inputdir))//trim(thickness_file)
  if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/SNOW_THICKNESS_FILE", filename)
  call get_param(PF, mdl, "SNOW_THICKNESS_IC_VAR", thick_var, &
                 "The initial condition variable for snow mass per unit area.", &
                 default="mH_ice", do_not_log=just_read)
  if (just_read) return ! All run-time parameters have been read, so return.

  if (.not.file_exists(filename, G%Domain)) call SIS_error(FATAL, &
         " initialize_snow_thick_from_file: Unable to open "//trim(filename))

  do k=2,CatIce ; mH_snow(:,:,k) = 0.0 ; enddo
  call MOM_read_data(filename, thick_var, mH_snow(:,:,1), G%Domain, scale=US%kg_m3_to_R*US%m_to_Z)

  call callTree_leave(trim(mdl)//'()')

end subroutine initialize_snow_thick_from_file


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_salinity_from_file reads the ice salinity from a file with 2-d or 3-d
!! (horizontal position and depth).  This salinity is used for all ice thickness categories
subroutine initialize_salinity_from_file(salin, G, IG, US, PF, just_read_params)
  type(SIS_hor_grid_type), intent(in)  :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)  :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,IG%NkIce), &
                           intent(out) :: salin !< The ice salinity that is being initialized [gSalt kg-1]
  type(param_file_type),   intent(in)  :: PF   !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing salinty.

  ! Local variables
  logical :: just_read    ! If true, just read parameters but set nothing.
  logical :: file_is_2d   ! If true, the salinity file has 2-d data.  Otherwise it includes a depth profile.
  real :: salin_scale     ! A scaling factor to use when reading salinity.
  character(len=40)  :: mdl = "initialize_salinity_from_file" ! This subroutine's name.
  character(len=64)  :: salin_var ! Ice salinity variable name in files
  character(len=200) :: filename, salinity_file, inputdir ! Strings for file/path
  integer :: i, j, k, n, is, ie, js, je, CatIce, NkIce

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; CatIce = IG%CatIce ; NkIce = IG%NkIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), SIS_state_initialization.F90")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=just_read)
  call get_param(PF, mdl, "ICE_SALINITY_FILE", salinity_file, &
                 "The name of the sea ice salinity file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)

  filename = trim(slasher(inputdir))//trim(salinity_file)
  if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/ICE_SALINITY_FILE", filename)
  call get_param(PF, mdl, "ICE_SALINITY_IC_VAR", salin_var, &
                 "The initial condition variable for sea ice salinity.", &
                 default="salinity", do_not_log=just_read)
  call get_param(PF, mdl, "ICE_SALINITY_IC_RESCALE", salin_scale, &
                 "A rescaling factor to use when reading ice salnity.", &
                 default=1.0, units="nondim", do_not_log=just_read)
  call get_param(PF, mdl, "ICE_SALINITY_FILE_IS_2D", file_is_2d, &
                 "If true, the salinity file has 2-data; otherwise it includes salinity profiles", &
                 default=.true., do_not_log=just_read)
  if (just_read) return ! All run-time parameters have been read, so return.

  if (.not.file_exists(filename, G%Domain)) call SIS_error(FATAL, &
         " initialize_salinity_from_file: Unable to open "//trim(filename))

  if (file_is_2d) then
    call MOM_read_data(filename, salin_var, salin(:,:,1,1), G%Domain, scale=salin_scale)
    do n=1,NkIce ; do k=1,CatIce ; do j=js,je ; do i=is,ie
      salin(i,j,k,n) = salin(i,j,1,1)
    enddo ; enddo ; enddo ; enddo
  else
    call MOM_read_data(filename, salin_var, salin(:,:,1,:), G%Domain, scale=salin_scale)
    do n=1,NkIce ; do k=1,CatIce ; do j=js,je ; do i=is,ie
      salin(i,j,k,n) = salin(i,j,1,n)
    enddo ; enddo ; enddo ; enddo
  endif

  call callTree_leave(trim(mdl)//'()')

end subroutine initialize_salinity_from_file


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_ice_enthalpy_from_file reads the ice enthalpy from a file with 2-d or 3-d
!! (horizontal position and depth) enthalpy or temperature.  This enthalpy or temperature is
!! used along with the pre-set salinity to set the enthalpy for all ice thickness categories.
subroutine initialize_ice_enthalpy_from_file(enth_ice, sal_ice, G, IG, US, ITV, PF, just_read_params)
  type(SIS_hor_grid_type), intent(in)  :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)  :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,IG%NkIce), &
                           intent(out) :: enth_ice !< The ice enthalpy that is being initialized [Q ~> J kg-1]
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,IG%NkIce), &
                           intent(in)  :: sal_ice !< The ice salinity [gSalt kg-1]
  type(ice_thermo_type),   intent(in)  :: ITV  !< The ice themodynamics parameter structure.
  type(param_file_type),   intent(in)  :: PF   !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing ice enthalpy.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp_input_2d         ! Temporary 2-d (horizontal position) ice temperature array [degC]
  real, dimension(SZI_(G),SZJ_(G),IG%NkIce) :: &
    temp_input_3d         ! Temporary 3-d (horizontal position and depth) ice temperature array [degC]
  real :: S_col(IG%NkIce) ! Specified ice column salinity used for ice thermodynamics [gSalt kg-1]
  logical :: just_read    ! If true, just read parameters but set nothing.
  logical :: file_is_2d   ! If true, the ice_enthalpy file has 2-d data.  Otherwise it includes a depth profile.
  logical :: enthalpy_file   ! If true, the file has enthalpy data in [J kg-1]; otherwise it has
                             ! temperatures in [degC].
  logical :: spec_thermo_sal ! If true used a specified ice profile for thermodynamics instead
                             ! of the actual ice salinity profile.
  character(len=40)  :: mdl = "initialize_ice_enthalpy_from_file" ! This subroutine's name.
  character(len=64)  :: varname ! Ice enthalpy or temperature variable name in files
  character(len=200) :: filename, ice_file, inputdir ! Strings for file/path
  integer :: i, j, k, n, is, ie, js, je, CatIce, NkIce

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; CatIce = IG%CatIce ; NkIce = IG%NkIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), SIS_state_initialization.F90")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=just_read)
  call get_param(PF, mdl, "READ_ENTHALPY_FILE", enthalpy_file, &
                 "If true, the file being read has ice enthalpy; otherwise it has temperature", &
                 default=.false., do_not_log=just_read)
  if (enthalpy_file) then
    call get_param(PF, mdl, "ICE_ENTHALPY_FILE", ice_file, &
                 "The name of the sea ice enthalpy file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
    filename = trim(slasher(inputdir))//trim(ice_file)
    if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/ICE_ENTHALPY_FILE", filename)
    call get_param(PF, mdl, "ICE_ENTHALPY_IC_VAR", varname, &
                 "The initial condition variable for sea ice enthalpy.", &
                 default="ice_enthalpy", do_not_log=just_read)
    call get_param(PF, mdl, "ICE_ENTHALPY_FILE_IS_2D", file_is_2d, &
                 "If true, the ice enthalpy file has 2-data; otherwise it includes enthalpy profiles", &
                 default=.true., do_not_log=just_read)
  else
    call get_param(PF, mdl, "ICE_TEMPERATURE_FILE", ice_file, &
                 "The name of the sea ice temperature file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
    filename = trim(slasher(inputdir))//trim(ice_file)
    if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/ICE_TEMPERATURE_FILE", filename)
    call get_param(PF, mdl, "ICE_TEMPERATURE_IC_VAR", varname, &
                 "The initial condition variable for sea ice temperature.", &
                 default="ice_temperature", do_not_log=just_read)
    call get_param(PF, mdl, "ICE_TEMPERATURE_FILE_IS_2D", file_is_2d, &
                 "If true, the ice temperature file has 2-data; otherwise it includes temperature profiles", &
                 default=.true., do_not_log=just_read)
  endif
  if (just_read) return ! All run-time parameters have been read, so return.

  if (.not.file_exists(filename, G%Domain)) call SIS_error(FATAL, &
         " initialize_ice_enthalpy_from_file: Unable to open "//trim(filename))

  S_col(:) = 0.0
  call get_SIS2_thermo_coefs(ITV, ice_salinity=S_col, spec_thermo_salin=spec_thermo_sal)

  if (file_is_2d) then
    if (enthalpy_file) then
      call MOM_read_data(filename, varname, enth_ice(:,:,1,1), G%Domain, scale=US%J_kg_to_Q)
      do n=1,NkIce ; do k=1,CatIce ; do j=js,je ; do i=is,ie
        enth_ice(i,j,k,n) = enth_ice(i,j,1,1)
      enddo ; enddo ; enddo ; enddo
    else
      call MOM_read_data(filename, varname, temp_input_2d, G%Domain)
      if (spec_thermo_sal) then
        do n=1,NkIce ; do k=1,CatIce ; do j=js,je ; do i=is,ie
          enth_ice(i,j,k,n) = Enth_from_TS(temp_input_2d(i,j), S_col(n), ITV)
        enddo ; enddo ; enddo ; enddo
      else
        do n=1,NkIce ; do k=1,CatIce ; do j=js,je ; do i=is,ie
          enth_ice(i,j,k,n) = Enth_from_TS(temp_input_2d(i,j), sal_ice(i,j,k,n), ITV)
        enddo ; enddo ; enddo ; enddo
      endif
    endif
  else
    if (enthalpy_file) then
      call MOM_read_data(filename, varname, enth_ice(:,:,1,:), G%Domain, scale=US%J_kg_to_Q)
      do n=1,NkIce ; do k=1,CatIce ; do j=js,je ; do i=is,ie
        enth_ice(i,j,k,n) = enth_ice(i,j,1,n)
      enddo ; enddo ; enddo ; enddo
    else
      call MOM_read_data(filename, varname, temp_input_3d, G%Domain)
      if (spec_thermo_sal) then
        do n=1,NkIce ; do k=1,CatIce ; do j=js,je ; do i=is,ie
          enth_ice(i,j,k,n) = Enth_from_TS(temp_input_3d(i,j,n), S_col(n), ITV)
        enddo ; enddo ; enddo ; enddo
      else
        do n=1,NkIce ; do k=1,CatIce ; do j=js,je ; do i=is,ie
          enth_ice(i,j,k,n) = Enth_from_TS(temp_input_3d(i,j,n), sal_ice(i,j,k,n), ITV)
        enddo ; enddo ; enddo ; enddo
      endif
    endif
  endif

  call callTree_leave(trim(mdl)//'()')

end subroutine initialize_ice_enthalpy_from_file


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> initialize_snow_enthalpy_from_file reads the snow enthalpy from a file with 2-d enthalpy or
!! temperature.  This enthalpy or temperature is used along with the pre-set salinity to set the
!! enthalpy for all snow thickness categories.
subroutine initialize_snow_enthalpy_from_file(enth_snow, G, IG, US, ITV, PF, just_read_params)
  type(SIS_hor_grid_type), intent(in)  :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)  :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)  :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),IG%CatIce,1), &
                           intent(out) :: enth_snow !< The snow enthalpy that is being initialized [Q ~> J kg-1]
  type(ice_thermo_type),   intent(in)  :: ITV  !< The ice themodynamics parameter structure.
  type(param_file_type),   intent(in)  :: PF   !< A structure indicating the open file
                                               !! to parse for model parameter values.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                               !! only read parameters without changing snow enthalpy.

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    temp_input_2d         ! Temporary 2-d (horizontal position) snow temperature array [degC]
  logical :: just_read    ! If true, just read parameters but set nothing.
  logical :: enthalpy_file  ! If true, the file has enthalpy data in [J kg-1]; otherwise it has
                            ! temperatures in [degC].
  character(len=40)  :: mdl = "initialize_snow_enthalpy_from_file" ! This subroutine's name.
  character(len=64)  :: varname ! Ice enthalpy or temperature variable name in files
  character(len=200) :: filename, snow_file, inputdir ! Strings for file/path
  integer :: i, j, k, is, ie, js, je, CatIce

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; CatIce = IG%CatIce

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) call callTree_enter(trim(mdl)//"(), SIS_state_initialization.F90")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=just_read)
  call get_param(PF, mdl, "READ_ENTHALPY_FILE", enthalpy_file, &
                 "If true, the file being read has snow enthalpy; otherwise it has temperature", &
                 default=.false., do_not_log=just_read)
  if (enthalpy_file) then
    call get_param(PF, mdl, "SNOW_ENTHALPY_FILE", snow_file, &
                 "The name of the snow enthalpy file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
    filename = trim(slasher(inputdir))//trim(snow_file)
    if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/SNOW_ENTHALPY_FILE", filename)
    call get_param(PF, mdl, "SNOW_ENTHALPY_IC_VAR", varname, &
                 "The initial condition variable for snow enthalpy.", &
                 default="snow_enthalpy", do_not_log=just_read)
  else
    call get_param(PF, mdl, "SNOW_TEMPERATURE_FILE", snow_file, &
                 "The name of the snow temperature file.", &
                 fail_if_missing=.not.just_read, do_not_log=just_read)
    filename = trim(slasher(inputdir))//trim(snow_file)
    if (.not.just_read) call log_param(PF, mdl, "INPUTDIR/SNOW_TEMPERATURE_FILE", filename)
    call get_param(PF, mdl, "SNOW_TEMPERATURE_IC_VAR", varname, &
                 "The initial condition variable for snow temperature.", &
                 default="snow_temperature", do_not_log=just_read)
  endif
  if (just_read) return ! All run-time parameters have been read, so return.

  if (.not.file_exists(filename, G%Domain)) call SIS_error(FATAL, &
         " initialize_snow_enthalpy_from_file: Unable to open "//trim(filename))

  if (enthalpy_file) then
    call MOM_read_data(filename, varname, enth_snow(:,:,1,1), G%Domain, scale=US%J_kg_to_Q)
    do k=1,CatIce ; do j=js,je ; do i=is,ie
      enth_snow(i,j,k,1) = enth_snow(i,j,1,1)
    enddo ; enddo ; enddo
  else
    call MOM_read_data(filename, varname, temp_input_2d, G%Domain)
    do k=1,CatIce ; do j=js,je ; do i=is,ie
      enth_snow(i,j,k,1) = Enth_from_TS(temp_input_2d(i,j), 0.0, ITV)
    enddo ; enddo ; enddo
  endif

  call callTree_leave(trim(mdl)//'()')

end subroutine initialize_snow_enthalpy_from_file


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> read_archaic_thermo_restarts is called if the restart file exists, but some parts of the sea ice state
!! have not been initialized using standard variable names, but might still exist.
subroutine read_archaic_thermo_restarts(Ice, IST, G, IG, US, PF, dirs, restart_file)
  type(ice_data_type),     intent(inout) :: Ice  !< The ice data type that is being initialized.
  type(ice_state_type),    intent(inout) :: IST  !< The sea ice state type being modified
  type(SIS_hor_grid_type), intent(in)    :: G    !< The horizontal grid type
  type(ice_grid_type),     intent(in)    :: IG   !< The sea-ice specific grid type
  type(unit_scale_type),   intent(in)    :: US   !< A structure with unit conversion factors
  type(param_file_type),   intent(in)    :: PF   !< A structure to parse for run-time parameters
  type(directories),       intent(in)    :: dirs          !< A structure containing several relevant directory paths.
  character(len=*),        intent(in)    :: restart_file  !< The restart file name to read;
                                                          !! the directory comes from dirs.

  ! Local variables
  real :: S_col(IG%NkIce) ! Specified ice column salinity used for ice thermodynamics [gSalt kg-1]
  logical :: spec_thermo_sal
  character(len=240) :: restart_path
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "SIS_state_initialization" ! This module's name.
  character(len=8)   :: nstr
  real, allocatable, target, dimension(:,:,:,:) :: t_ice_tmp
  real, allocatable, target, dimension(:,:,:) :: t_snow_tmp, sal_ice_tmp

  real :: ice_bulk_salin ! The globally constant sea ice bulk salinity [gSalt kg-1] = [ppt]
                         ! that is used to calculate the ocean salt flux.
  real :: ice_rel_salin  ! The initial bulk salinity of sea-ice relative to the
                         ! salinity of the water from which it formed [nondim].
  ! This pointers is used only for coding convenience.
  type(MOM_domain_type),   pointer :: sGD => NULL()

  integer :: i, j, k, l, i2, j2, k2, n
  integer :: isc, iec, jsc, jec, CatIce, NkIce
  logical :: read_values, read_t_ice(IG%NkIce)

  if (associated(Ice%sCS)) then ; if (.not.associated(Ice%sCS%IST)) then
    call SIS_error(FATAL, "read_archaic_thermo_restarts called with an unassociated Ice%sCS%Ice_state structure.")
  endif ; endif

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  CatIce = IG%CatIce ; NkIce = IG%NkIce

  ! Set some pointers for convenience.
  sGD => G%Domain

  call callTree_enter("read_archaic_thermo_restarts(), SIS_state_initialization.F90")

  ! Read all relevant parameters and write them to the model log.
!  call log_version(PF, mdl, version, "")

  call get_param(PF, mdl, "ICE_BULK_SALINITY", ice_bulk_salin, &
                 "The fixed bulk salinity of sea ice.", units = "g/kg", &
                 default=4.0, do_not_log=.true.)
  call get_param(PF, mdl, "ICE_RELATIVE_SALINITY", ice_rel_salin, &
                 "The initial salinity of sea ice as a fraction of the "//&
                 "salinity of the seawater from which it formed.", &
                 units = "nondim", default=0.0, do_not_log=.true.)
  if ((ice_bulk_salin < 0.0) .or. (ice_rel_salin > 0.0)) ice_bulk_salin = 0.0

  restart_path = trim(dirs%restart_input_dir)//trim(restart_file)

  ! Approximately initialize state fields that are not present in the restart files that were read.

  if (.not.query_initialized(Ice%Ice_restart, 'sal_ice')) then
    ! Initialize the ice salinity from separate variables for each layer, perhaps from a SIS1 restart.
    allocate(sal_ice_tmp(SZI_(G), SZJ_(G), CatIce)) ; sal_ice_tmp(:,:,:) = 0.0
    do n=1,NkIce
      write(nstr, '(I4)') n ; nstr = adjustl(nstr)
      call only_read_from_restarts(Ice%Ice_restart, 'sal_ice'//trim(nstr), sal_ice_tmp(:,:,:), &
                                   G%domain, directory=dirs%restart_input_dir, success=read_values)
      if (read_values) then
        do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%sal_ice(i,j,k,n) = sal_ice_tmp(i,j,k)
        enddo ; enddo ; enddo
      elseif (n==1) then
        IST%sal_ice(:,:,:,1) = ice_bulk_salin
      else
        IST%sal_ice(:,:,:,n) = IST%sal_ice(:,:,:,n-1)
      endif
    enddo

    deallocate(sal_ice_tmp)
  endif

  read_t_ice(:) = .false.
  if ((.not.query_initialized(Ice%Ice_restart, 'enth_ice')) .or. &
      (.not.query_initialized(Ice%Ice_restart, 'enth_snow'))) then
    allocate(t_ice_tmp(SZI_(G), SZJ_(G), CatIce, NkIce)) ; t_ice_tmp(:,:,:,:) = 0.0

    do n=1,NkIce
      write(nstr, '(I4)') n ; nstr = adjustl(nstr)
      call only_read_from_restarts(Ice%Ice_restart, 't_ice'//trim(nstr), t_ice_tmp(:,:,:,n), &
                                   G%domain, directory=dirs%restart_input_dir, success=read_values)
      read_t_ice(n) = read_values
    enddo
  endif

  ! Initialize the ice enthalpy.
  if (.not.query_initialized(Ice%Ice_restart, 'enth_ice')) then
    ! Try to initialize the ice enthalpy from separate temperature variables for each layer,
    ! perhaps from a SIS1 restart.
    if (.not.read_t_ice(1)) then
      call SIS_error(FATAL, "Either t_ice1 or enth_ice must be present in the SIS2 restart file "//restart_path)
    endif

    S_col(:) = 0.0
    call get_SIS2_thermo_coefs(IST%ITV, ice_salinity=S_col, spec_thermo_salin=spec_thermo_sal)

    do n=1,NkIce
      if ((n > 1) .and. (.not.read_t_ice(n))) &
        t_ice_tmp(:,:,:,n) = t_ice_tmp(:,:,:,n-1)

      if (spec_thermo_sal) then
        do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,n) = Enth_from_TS(t_ice_tmp(i,j,k,n), S_col(n), IST%ITV)
        enddo ; enddo ; enddo
      else
        do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          IST%enth_ice(i,j,k,n) = Enth_from_TS(t_ice_tmp(i,j,k,n), IST%sal_ice(i,j,k,n), IST%ITV)
        enddo ; enddo ; enddo
      endif
    enddo
  endif

  ! Initialize the snow enthalpy.
  if (.not.query_initialized(Ice%Ice_restart, 'enth_snow')) then
    ! Try to initialize the snow enthalpy from separate temperature variables for each layer,
    ! perhaps from a SIS1 restart.
    allocate(t_snow_tmp(SZI_(G), SZJ_(G), CatIce)) ; t_snow_tmp(:,:,:) = 0.0
    call only_read_from_restarts(Ice%Ice_restart, 't_snow', t_snow_tmp, G%domain, &
                                 directory=dirs%restart_input_dir, success=read_values)
    if (.not.read_values) then ! Try reading the ice temperature if snow is not available.
      if (read_t_ice(1)) then
        t_snow_tmp(:,:,:) = t_ice_tmp(:,:,:,1)
      else
        do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
          t_snow_tmp(i,j,k) = Temp_from_En_S(IST%enth_ice(i,j,k,1), IST%sal_ice(i,j,k,1), IST%ITV)
        enddo ; enddo ; enddo
      endif
    endif
    do k=1,CatIce ; do j=jsc,jec ; do i=isc,iec
      IST%enth_snow(i,j,k,1) = Enth_from_TS(t_snow_tmp(i,j,k), 0.0, IST%ITV)
    enddo ; enddo ; enddo
    deallocate(t_snow_tmp)
  endif

  if (allocated(t_ice_tmp)) deallocate(t_ice_tmp)

  call callTree_leave("read_archaic_thermo_restarts(), SIS_state_initialization.F90")

end subroutine read_archaic_thermo_restarts

end module SIS_state_initialization
