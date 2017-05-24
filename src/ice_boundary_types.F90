!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! ice_exchange_types contains the types that are used for exchanging information
!   with the atmosphere, land and ocean components via the FMS coupler.  These
!   types should be altered only in close coordination with the entire FMS
!   develoment effort.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
module ice_boundary_types

use coupler_types_mod, only : coupler_2d_bc_type, coupler_3d_bc_type
use fms_mod,           only : stdout
use mpp_mod,           only : mpp_chksum
use mpp_parameter_mod, only : CGRID_NE, BGRID_NE, AGRID

implicit none ; private

public :: ocean_ice_boundary_type, atmos_ice_boundary_type
public :: land_ice_boundary_type
public :: ocn_ice_bnd_type_chksum, atm_ice_bnd_type_chksum
public :: lnd_ice_bnd_type_chksum

!   The following three types are for data exchange with the FMS coupler
! they are defined here but declared in coupler_main and allocated in flux_init.

type ocean_ice_boundary_type
  real, dimension(:,:),   pointer :: &
    u      => NULL(), &  ! The x-direction ocean velocity at a position
                         ! determined by stagger, in m s-1.
    v      => NULL(), &  ! The y-direction ocean velocity at a position
                         ! determined by stagger, in m s-1.
    t      => NULL(), &  ! The ocean's surface temperature in Kelvin.
    s      => NULL(), &  ! The ocean's surface temperature in g/kg.
    frazil => NULL(), &  ! The frazil heat rejected by the ocean, in J m-2.
    sea_level => NULL()  ! The sea level after adjustment for any surface
                         ! pressure that the ocean allows to be expressed, in m.
  real, dimension(:,:,:), pointer :: data      =>NULL() ! collective field for "named" fields above
  integer                         :: stagger = BGRID_NE
  integer                         :: xtype              ! REGRID, REDIST or DIRECT used by coupler
  type(coupler_2d_bc_type)        :: fields     ! array of fields used for additional tracers
end type ocean_ice_boundary_type

type atmos_ice_boundary_type
  real, dimension(:,:,:), pointer :: &
    u_flux  => NULL(), & ! The true-eastward stresses (momentum fluxes) from the atmosphere
                         ! to the ocean or ice in each category, discretized on an A-grid,
                         ! and _not_ rotated to align with the model grid, in Pa.
    v_flux  => NULL(), & ! The true-northward stresses (momentum fluxes) from the atmosphere
                         ! to the ocean or ice in each category, discretized on an A-grid,
                         ! and _not_ rotated to align with the model grid, in Pa.
    u_star  => NULL(), & ! The atmospheric friction velocity on an A-grid, in Pa.
    t_flux  => NULL(), & ! The net sensible heat flux flux from the ocean or ice into the
                         ! atmosphere at the surface, in W m-2.
    q_flux  => NULL(), & ! The flux of moisture from the ice or ocean to the
                         ! atmosphere due to evaporation or sublimation, in kg m-2 s-1.
    lw_flux => NULL(), & ! The net flux of longwave radiation from the atmosphere into the
                         ! ice or ocean, in W m-2.
    ! sw_flux_tot_down => NULL(), & ! The total downward flux of shortwave radiation
    !                      ! at the surface of the ice or ocean, in W m-2.
    sw_flux_vis_dir => NULL(), & ! The visible (_vis) or near-infrared (_nir),
    sw_flux_vis_dif => NULL(), & ! direct (_dir) or diffuse (_dif) net shortwave
    sw_flux_nir_dir => NULL(), & ! radiation fluxes from the atmosphere into
    sw_flux_nir_dif => NULL(), & ! the ice or ocean, in W m-2.
    sw_down_vis_dir => NULL(), & ! The visible (_vis) or near-infrared (_nir),
    sw_down_vis_dif => NULL(), & ! direct (_dir) or diffuse (_dif) downward
    sw_down_nir_dir => NULL(), & ! shortwave radiation fluxes from the atmosphere
    sw_down_nir_dif => NULL(), & ! into the ice or ocean, in W m-2.

    lprec   => NULL(), & ! The liquid precipitation from the atmosphere onto the
                         ! atmosphere or ice in each thickness category, in kg m-2 s-1.
                         ! Rain falling on snow is currently assumed to pass or drain
                         ! directly through the ice into the ocean; this should be
                         ! revisited!
    fprec   => NULL(), & ! The frozen precipitation (snowfall) from the atmosphere
                         ! to the ice or ocean, in kg m-2 s-1.  Currently in SIS2
                         ! all frozen precipitation, including snow, sleet, hail
                         ! and graupel, are all treated as snow.
    dhdt    => NULL(), & ! The derivative of the upward sensible heat flux with the
                         ! surface temperature in W m-2 K-1.
    dedt    => NULL(), & ! The derivative of the sublimation and evaporation rate
                         ! with the surface temperature, in kg m-2 s-1 K-1.
    drdt    => NULL(), & ! The derivative of the net UPWARD longwave radiative
                         ! heat flux (-lw_flux) with surface temperature, in W m-2 K-1.
    coszen  => NULL(), & ! The cosine of the solar zenith angle averged over the
                         ! next radiation timestep (not the one that was used to
                         ! calculate the sw_flux fields), nondim and <=1.
    p       => NULL(), & ! The atmospheric surface pressure, in Pa, often ~1e5 Pa.
    data    => NULL()
  integer                   :: xtype  ! DIRECT or REDIST - used by coupler.
  type(coupler_3d_bc_type)  :: fluxes ! array of fluxes used for additional tracers
end type atmos_ice_boundary_type

type land_ice_boundary_type
  real, dimension(:,:),   pointer :: &
    runoff  =>NULL(), &  ! The liquid runoff into the ocean, in kg m-2.
    calving =>NULL(), &  ! The frozen runoff into each cell, that is offered
                         ! first to the icebergs (if any), where it might be
                         ! used or modified before being passed to the ocean,
                         ! in kg m-2.
    runoff_hflx  =>NULL(), & ! The heat flux associated with the temperature of
                             ! of the liquid runoff, relative to liquid water
                             ! at 0 deg C, in W m-2.
    calving_hflx =>NULL()    ! The heat flux associated with the temperature of
                             ! of the frozen runoff, relative to liquid? (or frozen?) water
                             ! at 0 deg C, in W m-2.
  real, dimension(:,:,:), pointer :: data    =>NULL() ! collective field for "named" fields above
  integer                         :: xtype  ! REGRID, REDIST or DIRECT - used by coupler.
end type land_ice_boundary_type

contains

subroutine ocn_ice_bnd_type_chksum(id, timestep, bnd_type)

  character(len=*), intent(in) :: id
  integer         , intent(in) :: timestep
  type(ocean_ice_boundary_type), intent(in) :: bnd_type
  integer ::   n, m, outunit

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(ocean_ice_boundary_type):: ', id, timestep
  write(outunit,100) 'ocn_ice_bnd_type%u        ',mpp_chksum(bnd_type%u        )
  write(outunit,100) 'ocn_ice_bnd_type%v        ',mpp_chksum(bnd_type%v        )
  write(outunit,100) 'ocn_ice_bnd_type%t        ',mpp_chksum(bnd_type%t        )
  write(outunit,100) 'ocn_ice_bnd_type%s        ',mpp_chksum(bnd_type%s        )
  write(outunit,100) 'ocn_ice_bnd_type%frazil   ',mpp_chksum(bnd_type%frazil   )
  write(outunit,100) 'ocn_ice_bnd_type%sea_level',mpp_chksum(bnd_type%sea_level)
  !    write(outunit,100) 'ocn_ice_bnd_type%data     ',mpp_chksum(bnd_type%data     )
  100 FORMAT("CHECKSUM::",A32," = ",Z20)

  do n = 1, bnd_type%fields%num_bcs  !{
    do m = 1, bnd_type%fields%bc(n)%num_fields  !{
        write(outunit,101) 'oibt%',trim(bnd_type%fields%bc(n)%name), &
             trim(bnd_type%fields%bc(n)%field(m)%name), &
             mpp_chksum(bnd_type%fields%bc(n)%field(m)%values)
    enddo  !} m
  enddo  !} n
  101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

end subroutine ocn_ice_bnd_type_chksum

subroutine atm_ice_bnd_type_chksum(id, timestep, bnd_type)
  character(len=*), intent(in) :: id
  integer         , intent(in) :: timestep
  type(atmos_ice_boundary_type), intent(in) :: bnd_type
  integer ::   n, outunit

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(atmos_ice_boundary_type):: ', id, timestep
  write(outunit,100) 'atm_ice_bnd_type%u_flux          ',mpp_chksum(bnd_type%u_flux)
  write(outunit,100) 'atm_ice_bnd_type%v_flux          ',mpp_chksum(bnd_type%v_flux)
  write(outunit,100) 'atm_ice_bnd_type%u_star          ',mpp_chksum(bnd_type%u_star)
  write(outunit,100) 'atm_ice_bnd_type%t_flux          ',mpp_chksum(bnd_type%t_flux)
  write(outunit,100) 'atm_ice_bnd_type%q_flux          ',mpp_chksum(bnd_type%q_flux)
  write(outunit,100) 'atm_ice_bnd_type%lw_flux         ',mpp_chksum(bnd_type%lw_flux)
  write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dir ',mpp_chksum(bnd_type%sw_flux_vis_dir)
  write(outunit,100) 'atm_ice_bnd_type%sw_flux_vis_dif ',mpp_chksum(bnd_type%sw_flux_vis_dif)
  write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dir ',mpp_chksum(bnd_type%sw_flux_nir_dir)
  write(outunit,100) 'atm_ice_bnd_type%sw_flux_nir_dif ',mpp_chksum(bnd_type%sw_flux_nir_dif)
  if (associated(bnd_type%sw_down_vis_dir)) &
    write(outunit,100) 'atm_ice_bnd_type%sw_down_vis_dir ',mpp_chksum(bnd_type%sw_down_vis_dir)
  if (associated(bnd_type%sw_down_vis_dif)) &
    write(outunit,100) 'atm_ice_bnd_type%sw_down_vis_dif ',mpp_chksum(bnd_type%sw_down_vis_dif)
  if (associated(bnd_type%sw_down_nir_dir)) &
    write(outunit,100) 'atm_ice_bnd_type%sw_down_nir_dir ',mpp_chksum(bnd_type%sw_down_nir_dir)
  if (associated(bnd_type%sw_down_nir_dif)) &
    write(outunit,100) 'atm_ice_bnd_type%sw_down_nir_dif ',mpp_chksum(bnd_type%sw_down_nir_dif)
  write(outunit,100) 'atm_ice_bnd_type%lprec           ',mpp_chksum(bnd_type%lprec)
  write(outunit,100) 'atm_ice_bnd_type%fprec           ',mpp_chksum(bnd_type%fprec)
  write(outunit,100) 'atm_ice_bnd_type%dhdt            ',mpp_chksum(bnd_type%dhdt)
  write(outunit,100) 'atm_ice_bnd_type%dedt            ',mpp_chksum(bnd_type%dedt)
  write(outunit,100) 'atm_ice_bnd_type%drdt            ',mpp_chksum(bnd_type%drdt)
  write(outunit,100) 'atm_ice_bnd_type%coszen          ',mpp_chksum(bnd_type%coszen)
  write(outunit,100) 'atm_ice_bnd_type%p               ',mpp_chksum(bnd_type%p)
!    write(outunit,100) 'atm_ice_bnd_type%data            ',mpp_chksum(bnd_type%data)
100 FORMAT("CHECKSUM::",A32," = ",Z20)

end subroutine atm_ice_bnd_type_chksum

subroutine lnd_ice_bnd_type_chksum(id, timestep, bnd_type)
  character(len=*), intent(in) :: id
  integer         , intent(in) :: timestep
  type(land_ice_boundary_type), intent(in) :: bnd_type
  integer ::   n, outunit

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(land_ice_boundary_type):: ', id, timestep
  write(outunit,100) 'lnd_ice_bnd_type%runoff  ',mpp_chksum(bnd_type%runoff)
  write(outunit,100) 'lnd_ice_bnd_type%calving ',mpp_chksum(bnd_type%calving)
  write(outunit,100) 'lnd_ice_bnd_type%runoff_hflx ',mpp_chksum(bnd_type%runoff_hflx)
  write(outunit,100) 'lnd_ice_bnd_type%calving_hflx',mpp_chksum(bnd_type%calving_hflx)
  !    write(outunit,100) 'lnd_ice_bnd_type%data    ',mpp_chksum(bnd_type%data)
  100 FORMAT("CHECKSUM::",A32," = ",Z20)
end subroutine lnd_ice_bnd_type_chksum

end module ice_boundary_types
