!> sea ice and SST specified from data as per GFDL climate group
module ice_spec_mod

use fms_mod, only : write_version_number
use mpp_mod, only : input_nml_file

use MOM_data_override, only : data_override, data_override_init, data_override_unset_domains
use MOM_error_handler, only : stdlog, stdout
use MOM_hor_index,     only : hor_index_type
use MOM_io,            only : open_namelist_file, check_nml_error, close_file
use MOM_time_manager,  only : time_type, get_date, set_date
use SIS_framework,     only : domain2d

implicit none ;  private

include 'netcdf.inc'

public :: get_sea_surface

logical :: module_is_initialized = .false. !< If true, this module has been called before.

logical :: mcm_ice = .false. !< When mcm_ice=.true., ice is handled as in supersource
real    :: sst_pert = 0.     !< global temperature perturbation used for sensitivity experiments [degC]

real    :: minimum_ice_concentration = 0.2 !< A minimum ice concentration [nondim]
real    :: minimum_ice_thickness     = 1.0 !< A minimum ice thickness [m]
logical :: do_leads = .true.   !< when do_leads=false there is no fractional ice concentration
                               !! also you should set the minimum_ice_concentration = 0.5
logical :: sst_degk = .false.  !< when sst_degk=true the input sst data is in degrees Kelvin
                               !! otherwise it is assumed to be in degrees Celsius

 integer :: repeat_date(3)=(/-1,-1,-1/) !< amip date for repeating single day (rsd) option

namelist / ice_spec_nml / mcm_ice, do_leads, minimum_ice_concentration, &
                          minimum_ice_thickness, sst_degk, sst_pert, repeat_date

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> get_sea_surface obtains some combination of SST, ice concentration and ice thickness, from data override files
subroutine get_sea_surface(Time, HI, SST, ice_conc, ice_thick, ice_domain, ice_domain_end)
  type (time_type),         intent(in)  :: Time !< The current model time
  type(hor_index_type),     intent(in)  :: HI  !< The horizontal index type describing the domain
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                  optional, intent(out) :: SST   !< The surface temperature [degC]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                  optional, intent(out) :: ice_conc   !< The fractional ice concentration [nondim]
  real, dimension(HI%isd:HI%ied,HI%jsd:HI%jed), &
                  optional, intent(out) :: ice_thick !< The ice thickness [m]
  type(domain2d), optional, intent(in)  :: ice_domain !< The domain used to read this data
  type(domain2d), optional, intent(in)  :: ice_domain_end !< If present reset the data override ice
                                                   !! domain back to this one at the end of this routine

  ! These local variables do not need halos, so they are declared without them.
  real, dimension(HI%isc:HI%iec,HI%jsc:HI%jec) :: sst_obs ! Observed sea surface temperature [degC] or [degK]
  real, dimension(HI%isc:HI%iec,HI%jsc:HI%jec) :: icec ! Observed sea ice concentration [nondim]
  real, dimension(HI%isc:HI%iec,HI%jsc:HI%jec) :: iceh ! Observed sea ice thickness [m]

  character(len=128), parameter :: version = '$Id: ice_spec.F90,v 1.1.2.1.6.1.2.1.2.1 2013/06/18 22:24:14 nnz Exp $'
  character(len=128), parameter :: tagname = '$Name: siena_201305_ice_sis2_5layer_dEdd_nnz $'
  real ::  t_sw_freeze0 = -1.8
  real ::  t_sw_freeze
  real :: SST_offset ! An offset between the observed SST and the output value, either to correct
                     ! for differences in temperature units or to apply a perturbation [degC]
  real, parameter :: T_0degC = 273.15 ! 0 degrees C in Kelvin [degK]
  integer :: i, j
  integer :: ierr, io, unit
  type(time_type) :: Spec_Time
  integer :: tod(3), dum1, dum2, dum3

  if (.not.module_is_initialized) then
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=ice_spec_nml, iostat=io)
#else
    unit = open_namelist_file()
    read  (unit, ice_spec_nml,iostat=io)
    call close_file (unit)
#endif
    ierr = check_nml_error(io,'ice_spec_nml')
    write (stdout(),'(/)')
    write (stdout(), ice_spec_nml)
    write (stdlog(), ice_spec_nml)

    call write_version_number(version, tagname)
    module_is_initialized = .true.
  endif

  if (present(ice_domain)) then
    call data_override_unset_domains(unset_Ice=.true., must_be_set=.false.)
    call data_override_init(Ice_domain_in = Ice_domain)
  endif

! modify time repeating single day option
  if (all(repeat_date>0)) then
    call get_date(Time,dum1,dum2,dum3,tod(1),tod(2),tod(3))
    Spec_Time = set_date(repeat_date(1),repeat_date(2),repeat_date(3),tod(1),tod(2),tod(3))
  else
    Spec_Time = Time
  endif

  icec(:,:) = 0.0 ; iceh(:,:) = 0.0
  call data_override('ICE', 'sic_obs', icec, Spec_Time)
  call data_override('ICE', 'sit_obs', iceh, Spec_Time)

  if (present(SST)) then
    t_sw_freeze = t_sw_freeze0
    if (sst_degk) then
      t_sw_freeze = t_sw_freeze0 + T_0degC ! convert sea water freeze point to degK
    endif

    sst_obs(:,:) = t_sw_freeze
    call data_override('ICE', 'sst_obs', sst_obs,  Spec_Time)
  endif

  if (present(ice_domain)) then
    call data_override_unset_domains(unset_Ice=.true.)

    ! Reset the data_override ice domain back to the one expected by the coupler.
    if (present(ice_domain_end)) call data_override_init(Ice_domain_in = Ice_domain_end)
  endif

  if (mcm_ice) then
    icec = 0.0
!   TK Mod: Limit minimum non-zero sea ice thickness to 0.01m.
!           This is to eliminate some very thin but non-zero
!           sea ice thickness values, where they really should be zero
!           but have become nonzero due to spatial interpolation
!           where the input grid and model grid are not
!           EXACTLY the same.  0.01 was obtained by trial and
!           error to roughly match supersource behavior.
!           5/22/01; 8/23/01

    where (iceh < 0.01) iceh=0.0
    where (iceh > 0.0)
      icec = 1.0
    end where
  else
    where (icec >= minimum_ice_concentration)
      iceh = max(iceh, minimum_ice_thickness)
    elsewhere
      icec = 0.0
      iceh = 0.0
    end where
    if (.not.do_leads) then
      where (icec >= minimum_ice_concentration) icec = 1.0
    endif
  endif

  if (present(ice_thick)) then
    do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      ice_thick(i,j) = iceh(i,j)
    enddo ; enddo
  endif

  if (present(ice_conc)) then
    do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      ice_conc(i,j) = icec(i,j)
    enddo ; enddo
  endif

  if (present(SST)) then
    ! SST is in Celsius, but sst_obs may be in Kelvin.
    SST_offset = 0.0 ; if (sst_degk) SST_offset = -T_0degC

    ! Add on non-zero sea surface temperature perturbation (namelist option)
    ! this perturbation may be useful in accessing model sensitivities
    if ( abs(sst_pert) > 0.0001 ) SST_offset = SST_offset + sst_pert

    do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
      if (icec(i,j) > 0.0) sst_obs(i,j) = t_sw_freeze
      if ((icec(i,j) == 0.0) .and. (sst_obs(i,j) <= t_sw_freeze)) sst_obs(i,j) = t_sw_freeze + 1e-10

      SST(i,j) = sst_obs(i,j) + SST_offset
    enddo ; enddo
  endif

end subroutine get_sea_surface

end module ice_spec_mod
