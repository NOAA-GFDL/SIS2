!> Contains the types that are used for exchanging information  with the atmosphere, land and ocean
!! components via the FMS coupler.
module ice_boundary_types

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_exchange_types contains the types that are used for exchanging information
!   with the atmosphere, land and ocean components via the FMS coupler.  These
!   types should be altered only in close coordination with the entire FMS
!   development effort.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

use MOM_error_handler, only : stdout, is_root_pe
use MOM_domains,       only : CGRID_NE, BGRID_NE, AGRID
use SIS_framework,     only : coupler_2d_bc_type, coupler_3d_bc_type
use SIS_framework,     only : SIS_chksum, coupler_type_write_chksums
use iso_fortran_env,   only : int64

implicit none ; private

public :: ocean_ice_boundary_type, atmos_ice_boundary_type
public :: land_ice_boundary_type
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum

!   The following three types are for data exchange with the FMS coupler
! they are defined here but declared in coupler_main and allocated in flux_init.

!> A type for exchange between the ocean and the sea ice
type ocean_ice_boundary_type
  real, dimension(:,:), pointer :: &
    u      => NULL(), &  !< The x-direction ocean velocity at a position
                         !! determined by stagger [m s-1].
    v      => NULL(), &  !< The y-direction ocean velocity at a position
                         !! determined by stagger [m s-1].
    t      => NULL(), &  !< The ocean's surface temperature [Kelvin].
    s      => NULL(), &  !< The ocean's surface salinity [gSalt kg-1].
    frazil => NULL(), &  !< The frazil heat rejected by the ocean [J m-2].
    sea_level => NULL()  !< The sea level after adjustment for any surface
                         !! pressure that the ocean allows to be expressed [m].
  real, dimension(:,:,:), pointer :: data =>NULL() !< S collective field for "named" fields above
  integer   :: stagger = BGRID_NE  !< A flag indicating how the velocities are staggered.
  integer   :: xtype     !< A flag indicating the exchange type, which may be set to
                         !! REGRID, REDIST or DIRECT and is used by coupler
  type(coupler_2d_bc_type) :: fields !< An array of fields used for additional tracers
end type ocean_ice_boundary_type

!> A type for exchange between the atmosphere and the sea ice
type atmos_ice_boundary_type
  real, dimension(:,:,:), pointer :: &
    u_flux  => NULL(), & !< The true-eastward stresses (momentum fluxes) from the atmosphere
                         !! to the ocean or ice in each category, discretized on an A-grid,
                         !! and _not_ rotated to align with the model grid [Pa].
    v_flux  => NULL(), & !< The true-northward stresses (momentum fluxes) from the atmosphere
                         !! to the ocean or ice in each category, discretized on an A-grid,
                         !! and _not_ rotated to align with the model grid [Pa].
    u_star  => NULL(), & !< The atmospheric friction velocity on an A-grid [Pa].
    t_flux  => NULL(), & !< The net sensible heat flux flux from the ocean or ice into the
                         !! atmosphere at the surface [W m-2].
    q_flux  => NULL(), & !< The flux of moisture from the ice or ocean to the
                         !! atmosphere due to evaporation or sublimation [kg m-2 s-1].
    lw_flux => NULL(), & !< The net flux of longwave radiation from the atmosphere into the
                         !! ice or ocean [W m-2].
    !! sw_flux_tot_down => NULL(), & !< The total downward flux of shortwave radiation
    !!                      !! at the surface of the ice or ocean [W m-2].
    sw_flux_vis_dir => NULL(), & !< The visible (_vis) or near-infrared (_nir),
    sw_flux_vis_dif => NULL(), & !< direct (_dir) or diffuse (_dif) net shortwave
    sw_flux_nir_dir => NULL(), & !< radiation fluxes from the atmosphere into
    sw_flux_nir_dif => NULL(), & !< the ice or ocean [W m-2].
    sw_down_vis_dir => NULL(), & !< The visible (_vis) or near-infrared (_nir),
    sw_down_vis_dif => NULL(), & !< direct (_dir) or diffuse (_dif) downward
    sw_down_nir_dir => NULL(), & !< shortwave radiation fluxes from the atmosphere
    sw_down_nir_dif => NULL(), & !< into the ice or ocean [W m-2].

    lprec   => NULL(), & !< The liquid precipitation from the atmosphere onto the
                         !! atmosphere or ice in each thickness category [kg m-2 s-1].
                         !! Rain falling on snow is currently assumed to pass or drain
                         !! directly through the ice into the ocean; this should be
                         !! revisited!
    fprec   => NULL(), & !< The frozen precipitation (snowfall) from the atmosphere
                         !! to the ice or ocean [kg m-2 s-1].  Currently in SIS2
                         !! all frozen precipitation, including snow, sleet, hail
                         !! and graupel, are all treated as snow.
    dhdt    => NULL(), & !< The derivative of the upward sensible heat flux with the
                         !! surface temperature [W m-2 degC-1].
    dedt    => NULL(), & !< The derivative of the sublimation and evaporation rate
                         !! with the surface temperature [kg m-2 s-1 degC-1].
    drdt    => NULL(), & !< The derivative of the net UPWARD longwave radiative
                         !! heat flux (-lw_flux) with surface temperature [W m-2 degC-1].
    coszen  => NULL(), & !< The cosine of the solar zenith angle averaged over the
                         !! next radiation timestep (not the one that was used to
                         !! calculate the sw_flux fields), nondim and <=1.
    p       => NULL()    !< The atmospheric surface pressure [Pa], often ~1e5 Pa.
!    data    => NULL() ! This can probably be removed.
  integer   :: xtype     !< A flag indicating the exchange type, which may be set to
                         !! REGRID, REDIST or DIRECT and is used by coupler
  type(coupler_3d_bc_type)  :: fluxes !< An array of fluxes used for additional tracers
end type atmos_ice_boundary_type

!> A type for exchange between the land and the sea ice
type land_ice_boundary_type
  real, dimension(:,:),   pointer :: &
    runoff  =>NULL(), &  !< The liquid runoff into the ocean [kg m-2].
    calving =>NULL(), &  !< The frozen runoff into each cell, that is offered
                         !! first to the icebergs (if any), where it might be
                         !! used or modified before being passed to the ocean [kg m-2].
    runoff_hflx  =>NULL(), & !< The heat flux associated with the temperature of
                         !! of the liquid runoff, relative to liquid water
                         !! at 0 degC [W m-2].
    calving_hflx =>NULL() !< The heat flux associated with the temperature of
                         !! of the frozen runoff, relative to liquid? (or frozen?) water
                         !! at 0 degC [W m-2].
  real, dimension(:,:,:), pointer :: data => NULL() !< A collective field for "named" fields above
  integer   :: xtype     !< A flag indicating the exchange type, which may be set to
                         !! REGRID, REDIST or DIRECT and is used by coupler
end type land_ice_boundary_type

contains

!> Write checksums of the fields in an ocean_ice_boundary_type
subroutine ocn_ice_bnd_type_chksum(id, timestep, bnd_type)

  character(len=*), intent(in) :: id !< An identifying message fragment
  integer         , intent(in) :: timestep !< The timestep number
  type(ocean_ice_boundary_type), intent(in) :: bnd_type !< The structure whose elements are to be checksummed

  ! Local variables
  integer(kind=int64) :: chks ! A checksum for the field
  logical :: root    ! True only on the root PE
  integer :: outunit ! The output unit to write to

  outunit = stdout()
  root = is_root_pe()

  if (root) write(outunit,*) 'BEGIN CHECKSUM(ocean_ice_boundary_type):: ', id, timestep
  chks = SIS_chksum(bnd_type%u        ) ; if (root) write(outunit,100) 'ocn_ice_bnd_type%u        ', chks
  chks = SIS_chksum(bnd_type%v        ) ; if (root) write(outunit,100) 'ocn_ice_bnd_type%v        ', chks
  chks = SIS_chksum(bnd_type%t        ) ; if (root) write(outunit,100) 'ocn_ice_bnd_type%t        ', chks
  chks = SIS_chksum(bnd_type%s        ) ; if (root) write(outunit,100) 'ocn_ice_bnd_type%s        ', chks
  chks = SIS_chksum(bnd_type%frazil   ) ; if (root) write(outunit,100) 'ocn_ice_bnd_type%frazil   ', chks
  chks = SIS_chksum(bnd_type%sea_level) ; if (root) write(outunit,100) 'ocn_ice_bnd_type%sea_level', chks
  !    write(outunit,100) 'ocn_ice_bnd_type%data     ', SIS_chksum(bnd_type%data     )
  100 FORMAT("CHECKSUM::",A32," = ",Z20)

  call coupler_type_write_chksums(bnd_type%fields, outunit, 'oibt%')

end subroutine ocn_ice_bnd_type_chksum

!> Write checksums of the fields in an atmos_ice_boundary_type
subroutine atm_ice_bnd_type_chksum(id, timestep, bnd_type)
  character(len=*), intent(in) :: id !< An identifying message fragment
  integer         , intent(in) :: timestep !< The timestep number
  type(atmos_ice_boundary_type), intent(in) :: bnd_type !< The structure whose elements are to be checksummed

  ! Local variables
  integer(kind=int64) :: chks ! A checksum for the field
  logical :: root    ! True only on the root PE
  integer :: outunit ! The output unit to write to

  outunit = stdout()
  root = is_root_pe()

  if (root) write(outunit,*) 'BEGIN CHECKSUM(atmos_ice_boundary_type):: ', id, timestep
  chks = SIS_chksum(bnd_type%u_flux)  ; if (root) write(outunit,100) 'atm_ice_bnd_type%u_flux  ', chks
  chks = SIS_chksum(bnd_type%v_flux)  ; if (root) write(outunit,100) 'atm_ice_bnd_type%v_flux  ', chks
  chks = SIS_chksum(bnd_type%u_star)  ; if (root) write(outunit,100) 'atm_ice_bnd_type%u_star  ', chks
  chks = SIS_chksum(bnd_type%t_flux)  ; if (root) write(outunit,100) 'atm_ice_bnd_type%t_flux  ', chks
  chks = SIS_chksum(bnd_type%q_flux)  ; if (root) write(outunit,100) 'atm_ice_bnd_type%q_flux  ', chks
  chks = SIS_chksum(bnd_type%lw_flux) ; if (root) write(outunit,100) 'atm_ice_bnd_type%lw_flux ', chks
  chks = SIS_chksum(bnd_type%sw_flux_vis_dir) ; if (root) write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dir ', chks
  chks = SIS_chksum(bnd_type%sw_flux_vis_dif) ; if (root) write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dif ', chks
  chks = SIS_chksum(bnd_type%sw_flux_nir_dir) ; if (root) write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dir ', chks
  chks = SIS_chksum(bnd_type%sw_flux_nir_dif) ; if (root) write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dif ', chks
  if (associated(bnd_type%sw_down_vis_dir)) then
    chks = SIS_chksum(bnd_type%sw_down_vis_dir) ; if (root) write(outunit,100) 'atm_ice_bnd_type%sw_down_vis_dir ', chks
  endif
  if (associated(bnd_type%sw_down_vis_dif)) then
    chks = SIS_chksum(bnd_type%sw_down_vis_dif) ; if (root) write(outunit,100) 'atm_ice_bnd_type%sw_down_vis_dif ', chks
  endif
  if (associated(bnd_type%sw_down_nir_dir)) then
    chks = SIS_chksum(bnd_type%sw_down_nir_dir) ; if (root) write(outunit,100) 'atm_ice_bnd_type%sw_down_nir_dir ', chks
  endif
  if (associated(bnd_type%sw_down_nir_dif)) then
    chks = SIS_chksum(bnd_type%sw_down_nir_dif) ; if (root) write(outunit,100) 'atm_ice_bnd_type%sw_down_nir_dif ', chks
  endif
  chks = SIS_chksum(bnd_type%lprec)  ; if (root) write(outunit,100) 'atm_ice_bnd_type%lprec  ', chks
  chks = SIS_chksum(bnd_type%fprec)  ; if (root) write(outunit,100) 'atm_ice_bnd_type%fprec  ', chks
  chks = SIS_chksum(bnd_type%dhdt)   ; if (root) write(outunit,100) 'atm_ice_bnd_type%dhdt   ', chks
  chks = SIS_chksum(bnd_type%dedt)   ; if (root) write(outunit,100) 'atm_ice_bnd_type%dedt   ', chks
  chks = SIS_chksum(bnd_type%drdt)   ; if (root) write(outunit,100) 'atm_ice_bnd_type%drdt   ', chks
  chks = SIS_chksum(bnd_type%coszen) ; if (root) write(outunit,100) 'atm_ice_bnd_type%coszen ', chks
  chks = SIS_chksum(bnd_type%p)      ; if (root) write(outunit,100) 'atm_ice_bnd_type%p      ', chks
  ! chks = SIS_chksum(bnd_type%data) ; if (root) write(outunit,100) 'atm_ice_bnd_type%data   ', chks
100 FORMAT("CHECKSUM::",A32," = ",Z20)

end subroutine atm_ice_bnd_type_chksum

!> Write checksums of the fields in a land_ice_boundary_type
subroutine lnd_ice_bnd_type_chksum(id, timestep, bnd_type)
  character(len=*), intent(in) :: id !< An identifying message fragment
  integer         , intent(in) :: timestep !< The timestep number
  type(land_ice_boundary_type), intent(in) :: bnd_type !< The structure whose elements are to be checksummed

  ! Local variables
  integer(kind=int64) :: chks ! A checksum for the field
  logical :: root    ! True only on the root PE
  integer :: outunit ! The output unit to write to

  outunit = stdout()
  root = is_root_pe()

  if (root) write(outunit,*) 'BEGIN CHECKSUM(land_ice_boundary_type):: ', id, timestep
  chks = SIS_chksum(bnd_type%runoff)       ; if (root) write(outunit,100) 'lnd_ice_bnd_type%runoff  ', chks
  chks = SIS_chksum(bnd_type%calving)      ; if (root) write(outunit,100) 'lnd_ice_bnd_type%calving ', chks
  chks = SIS_chksum(bnd_type%runoff_hflx)  ; if (root) write(outunit,100) 'lnd_ice_bnd_type%runoff_hflx ', chks
  chks = SIS_chksum(bnd_type%calving_hflx) ; if (root) write(outunit,100) 'lnd_ice_bnd_type%calving_hflx', chks
  ! chks = SIS_chksum(bnd_type%data) ; if (root) write(outunit,100) 'lnd_ice_bnd_type%data    ', chks
  100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine lnd_ice_bnd_type_chksum

end module ice_boundary_types
