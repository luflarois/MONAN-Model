module monan_chemistry_vars
    !! ## memory module for chemistry
    !!
    !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
    !! ## MONAN
    !!
    !! Author: Rodrigues, L.F.
    !!
    !! E-mail: luflarois@gmail.com
    !!
    !! Date: 2026-03-27
    !!
    !! #####Version: <>
    !!
    !! —
    !! **Full description**:
    !!
    !! Module container for chemistry variables. This module is used to declare the variables that are needed for the 
    !! chemistry parameterization. It also contains a subroutine to initialize the variables, which is currently a dummy 
    !! subroutine that does not do anything.
    !!
    !! ** History**:
    !!
    !! - 2026-03-27: Modified to eliminate j dimension and invert indices to (k,i) for consistency with chemistry interface.
    !!—
    !! ** Licence **:
    !!
    !! <img src='https://www.gnu.org/graphics/gplv3-127x51.png' width='63'>
    !!
    use mpas_kind_types
    
    implicit none
    public
    save

    !=================================================================================================================
    !wrf-variables: these variables are needed to keep calls to different physics parameterizations
    !as in wrf model.
    !=================================================================================================================

    logical:: l_radtlw                   !controls call to longwave radiation parameterization.
    logical:: l_radtsw                   !controls call to shortwave radiation parameterization.
    logical:: l_conv                     !controls call to convective parameterization.
    logical:: l_camlw                    !controls when to save local CAM LW abs and ems arrays.
    logical:: l_diags                    !controls when to calculate physics diagnostics.
    logical:: l_acrain                   !when .true., limit to accumulated rain is applied.
    logical:: l_acradt                   !when .true., limit to lw and sw radiation is applied.
    logical:: l_mp_tables                !when .true., read look-up tables for Thompson cloud microphysics scheme.

    integer,public:: ids,ide,jds,jde,kds,kde
    integer,public:: ims,ime,jms,jme,kms,kme
    integer,public:: its,ite,jts,jte,kts,kte
    integer,public:: iall
    integer,public:: n_microp

    integer,public:: num_months          !number of months                                         [-]

    real(kind=RKIND),public:: dt_dyn     !time-step for dynamics
    real(kind=RKIND),public:: dt_microp  !time-step for cloud microphysics parameterization.
    real(kind=RKIND),public:: dt_radtlw  !time-step for longwave radiation parameterization      [mns]
    real(kind=RKIND),public:: dt_radtsw  !time-step for shortwave radiation parameterization     [mns]

    real(kind=RKIND),public:: xice_threshold

    !... arrays related to surface: now 1D (i)
    real(kind=RKIND),dimension(:),allocatable:: &
        xlon_p,           &!longitude, west is negative                                      [degrees]
        xlat_p,           &!latitude, south is negative                                      [degrees]
        lwupb_p,          &!all-sky upwelling longwave flux at bottom-of-atmosphere          [J m-2]
        sfc_albedo_p,     &!surface albedo                                                   [-]
        psfc_p,           &!surface pressure                                                 [Pa]
        ptop_p,           &!model-top pressure                                               [Pa]
        coszr_p            !cosine of the solar zenith angle                                 [-]

    !... arrays at model levels (k,i)
    real(kind=RKIND),dimension(:,:),allocatable:: &
        o3_p,             &!ozone mixing ratio                                               [-]
        u_p,              &!u-velocity interpolated to theta points                          [m/s]
        v_p,              &!v-velocity interpolated to theta points                          [m/s]
        fzm_p,            &!weight for interpolation to w points                             [-]
        fzp_p,            &!weight for interpolation to w points                             [-]
        zz_p,             &!height (or something)                                            [m]
        pres_p,           &!pressure                                                         [Pa]
        pi_p,             &!(p_phy/p0)**(r_d/cp)                                             [-]
        z_p,              &!height of layer                                                  [m]
        zmid_p,           &!height of middle of layer                                        [m]
        dz_p,             &!layer thickness                                                  [m]
        t_p,              &!temperature                                                      [K]
        th_p,             &!potential temperature                                            [K]
        al_p,             &!inverse of air density                                           [m3/kg]
        rho_p,            &!air density                                                      [kg/m3]
        rh_p,             &!relative humidity                                                [-]
        znu_p,            &! (pres_hyd_p / P0) needed in the Tiedtke convection scheme       [Pa]
        w_p,              &!vertical velocity at w-points                                    [m/s]
        pres2_p,          &!pressure at w-points                                             [Pa]
        t2_p,             &!temperature at w-points                                          [K]
        qv_p,             &!water vapor mixing ratio                                         [kg/kg]
        qc_p,             &!cloud water mixing ratio                                         [kg/kg]
        qr_p,             &!rain mixing ratio                                                [kg/kg]
        qi_p,             &!cloud ice mixing ratio                                           [kg/kg]
        qs_p,             &!snow mixing ratio                                                [kg/kg]
        qg_p               !graupel mixing ratio                                             [kg/kg]

    real(kind=RKIND),dimension(:,:,:),allocatable:: &
        chem_conc_p,      &!chemical species mixing ratio (ppbm) (k,i,j) (after chemistry integration, before transport)
        chem_tend_p,      &!chemical species tendency (ppbm) (k,i,j) (after chemistry integration, before transport)
        chem_tend_dyn_p    !chemical species tendency (ppbm) (k,i,j) (after chemistry integration, before transport)
        

    !... arrays for hydrostatic pressure and exner (1D and 2D)
    real(kind=RKIND),dimension(:),allocatable:: &
        psfc_hyd_p,       &!surface pressure                                                 [Pa]
        psfc_hydd_p        !"dry" surface pressure                                           [Pa]
    real(kind=RKIND),dimension(:,:),allocatable:: &
        pres_hyd_p,       &!pressure located at theta levels                                 [Pa]
        pres_hydd_p,      &!"dry" pressure located at theta levels                           [Pa]
        pres2_hyd_p,      &!pressure located at w-velocity levels                            [Pa]
        pres2_hydd_p,     &!"dry" pressure located at w-velocity levels                      [Pa]
        znu_hyd_p          !(pres_hyd_p / P0) needed in the Tiedtke convection scheme        [Pa]

    !=================================================================================================================
    !... variables related to ozone climatology:
    !=================================================================================================================
    real(kind=RKIND),dimension(:,:),allocatable:: &
        o3clim_p           !climatological ozone volume mixing ratio                         [???] (k, i)

    !=================================================================================================================
    !... variables and arrays related to parameterization of cloud microphysics:
    !=================================================================================================================

    logical,parameter:: &
        warm_rain = .false. !warm-phase cloud microphysics only (used in WRF).

    logical:: &
        f_qc,             &!parameter set to true to include the cloud water mixing ratio.
        f_qr,             &!parameter set to true to include the rain mixing ratio.
        f_qi,             &!parameter set to true to include the cloud ice mixing ratio.
        f_qs,             &!parameter set to true to include the snow mixing ratio.
        f_qg,             &!parameter set to true to include the graupel mixing ratio.
        f_qoz              !parameter set to true to include the ozone mixing ratio.

    logical:: &
        f_nc,             &!parameter set to true to include the cloud water number concentration.
        f_ni,             &!parameter set to true to include the cloud ice number concentration.
        f_nifa,           &!parameter set to true to include the number concentration of hygroscopic aerosols.
        f_nwfa,           &!parameter set to true to include the number concentration of hydrophobic aerosols.
        f_nbca             !parameter set to true to include the number concentration of black carbon.

    ! (outras variáveis de outras físicas não modificadas permanecem como estavam)
    ! ... aqui você pode incluir as demais declarações que não foram alteradas, ou omitir se não são usadas.

contains

    subroutine chemistry_vars_init()
        ! dummy
    end subroutine chemistry_vars_init

end module monan_chemistry_vars