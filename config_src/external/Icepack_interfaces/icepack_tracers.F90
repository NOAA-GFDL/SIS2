!=======================================================================
! Indices and flags associated with the tracer infrastructure.
! Grid-dependent and max_trcr-dependent arrays are declared in ice_state.F90.
!
! author Elizabeth C. Hunke, LANL

      module icepack_tracers

      use icepack_kinds

      implicit none

      private
      public :: icepack_init_tracer_indices

      real(kind=dbl_kind), parameter, public :: n_iso=0,n_aero=0

      contains


!=======================================================================
!autodocument_start icepack_init_tracer_indices
! set the number of column tracer indices

      subroutine icepack_init_tracer_indices(&
           nt_Tsfc_in, nt_qice_in, nt_qsno_in, nt_sice_in, &
           nt_fbri_in, nt_iage_in, nt_FY_in, &
           nt_alvl_in, nt_vlvl_in, nt_apnd_in, nt_hpnd_in, nt_ipnd_in, &
           nt_fsd_in, nt_isosno_in, nt_isoice_in, &
           nt_aero_in, nt_zaero_in, nt_bgc_C_in, &
           nt_bgc_N_in, nt_bgc_chl_in, nt_bgc_DOC_in, nt_bgc_DON_in, &
           nt_bgc_DIC_in, nt_bgc_Fed_in, nt_bgc_Fep_in, nt_bgc_Nit_in, nt_bgc_Am_in, &
           nt_bgc_Sil_in, nt_bgc_DMSPp_in, nt_bgc_DMSPd_in, nt_bgc_DMS_in, nt_bgc_hum_in, &
           nt_bgc_PON_in, nlt_zaero_in, nlt_bgc_C_in, nlt_bgc_N_in, nlt_bgc_chl_in, &
           nlt_bgc_DOC_in, nlt_bgc_DON_in, nlt_bgc_DIC_in, nlt_bgc_Fed_in, &
           nlt_bgc_Fep_in, nlt_bgc_Nit_in, nlt_bgc_Am_in, nlt_bgc_Sil_in, &
           nlt_bgc_DMSPp_in, nlt_bgc_DMSPd_in, nlt_bgc_DMS_in, nlt_bgc_hum_in, &
           nlt_bgc_PON_in, nt_zbgc_frac_in, nt_bgc_S_in, nlt_chl_sw_in, &
           nlt_zaero_sw_in, &
           bio_index_o_in, bio_index_in)

        integer, intent(in), optional :: &
             nt_Tsfc_in, & ! ice/snow temperature
             nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
             nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage_in, & ! volume-weighted ice age
             nt_FY_in, & ! area-weighted first-year ice area
             nt_alvl_in, & ! level ice area fraction
             nt_vlvl_in, & ! level ice volume fraction
             nt_apnd_in, & ! melt pond area fraction
             nt_hpnd_in, & ! melt pond depth
             nt_ipnd_in, & ! melt pond refrozen lid thickness
             nt_fsd_in,  & ! floe size distribution
             nt_isosno_in,  & ! starting index for isotopes in snow
             nt_isoice_in,  & ! starting index for isotopes in ice
             nt_aero_in,    & ! starting index for aerosols in ice
             nt_bgc_Nit_in, & ! nutrients
             nt_bgc_Am_in,  & !
             nt_bgc_Sil_in, & !
             nt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nt_bgc_DMSPd_in,&!
             nt_bgc_DMS_in, & !
             nt_bgc_hum_in, & !
             nt_bgc_PON_in, & ! zooplankton and detritus
             nlt_bgc_Nit_in,& ! nutrients
             nlt_bgc_Am_in, & !
             nlt_bgc_Sil_in,& !
             nlt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd_in,&!
             nlt_bgc_DMS_in,& !
             nlt_bgc_hum_in,& !
             nlt_bgc_PON_in,& ! zooplankton and detritus
             nt_zbgc_frac_in,&! fraction of tracer in the mobile phase
             nt_bgc_S_in,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw_in    ! points to total chla in trcrn_sw

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             bio_index_o_in, &
             bio_index_in

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small
             nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small
             nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small
             nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small
             nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small
             nlt_bgc_chl_in   ! diatoms, phaeocystis, pico/small

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DOC_in, & !  dissolved organic carbon
             nlt_bgc_DOC_in   !  dissolved organic carbon

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DON_in, & !  dissolved organic nitrogen
             nlt_bgc_DON_in   !  dissolved organic nitrogen

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_DIC_in, & ! dissolved inorganic carbon
             nlt_bgc_DIC_in   !  dissolved inorganic carbon

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_bgc_Fed_in, & !  dissolved iron
             nt_bgc_Fep_in, & !  particulate iron
             nlt_bgc_Fed_in,& !  dissolved iron
             nlt_bgc_Fep_in   !  particulate iron

        integer (kind=int_kind), dimension(:), intent(in), optional :: &
             nt_zaero_in,   & !  black carbon and other aerosols
             nlt_zaero_in,  & !  black carbon and other aerosols
             nlt_zaero_sw_in  ! black carbon and dust in trcrn_sw

      end subroutine icepack_init_tracer_indices


      end module icepack_tracers

!=======================================================================
