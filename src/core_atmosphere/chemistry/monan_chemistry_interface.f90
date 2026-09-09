module monan_chemistry_interface
    use mpas_kind_types
    use mpas_pool_routines
    use mpas_atmphys_constants
    use monan_chemistry_vars
    use chem_list, only: nSpecies,spc_name
    use modMemoryChem, only: chem_g

    implicit none
    private
    public:: allocate_forall_chemistry,   &
             deallocate_forall_chemistry, &
             MPAS_to_chemistry,           &
             chemistry_to_MPAS, &
             zero_chemistry_scalars

    integer :: idxChem(nSpecies)

    contains

    !=================================================================================================================
    ! Converte uma string ASCII para minusculas. As dimensoes de indice
    ! geradas pelo Registry (index_o3, index_no2, ...) sao sempre em
    ! minusculo, mas spc_name (vindo de chem_list.f90) esta em maiusculo -
    ! por isso essa conversao e necessaria antes de montar a chave de busca.
    !=================================================================================================================
    pure function to_lower(str_in) result(str_out)
       character(len=*), intent(in) :: str_in
       character(len=len(str_in)) :: str_out
       integer :: i, code

       do i = 1, len(str_in)
          code = iachar(str_in(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             str_out(i:i) = achar(code + (iachar('a') - iachar('A')))
          else
             str_out(i:i) = str_in(i:i)
          end if
       end do
    end function to_lower

    !=================================================================================================================
    subroutine allocate_forall_chemistry(nCells, nVertLevels, nChemSpecies)
    !=================================================================================================================
    integer, intent(in) :: nCells, nVertLevels, nChemSpecies

!print *,'LFR-DBG: Block 1'
    ! Alocações inalteradas (já estão corretas)
    if(.not.allocated(xlon_p) ) allocate(xlon_p(nCells)               )
    if(.not.allocated(xlat_p) ) allocate(xlat_p(nCells)               )
    if(.not.allocated(lwupb_p)) allocate(lwupb_p(nCells)              )
    if(.not.allocated(sfc_albedo_p) ) allocate(sfc_albedo_p(nCells) )
    if(.not.allocated(psfc_p) ) allocate(psfc_p(nCells)         )
    if(.not.allocated(ptop_p) ) allocate(ptop_p(nCells)         )
    if(.not.allocated(coszr_p)      ) allocate(coszr_p(nCells)              )
!print *,'LFR-DBG: Block 2'
    if(.not.allocated(o3_p)   ) allocate(o3_p(nVertlevels, nCells)    )
    if(.not.allocated(u_p)    ) allocate(u_p(nVertlevels, nCells)    )
    if(.not.allocated(v_p)    ) allocate(v_p(nVertlevels, nCells)    )
    if(.not.allocated(fzm_p)  ) allocate(fzm_p(nVertlevels, nCells)  )
    if(.not.allocated(fzp_p)  ) allocate(fzp_p(nVertlevels, nCells)  )
    if(.not.allocated(zz_p)   ) allocate(zz_p(nVertlevels, nCells)   )
    if(.not.allocated(pres_p) ) allocate(pres_p(nVertlevels, nCells) )
    if(.not.allocated(pi_p)   ) allocate(pi_p(nVertlevels, nCells)   )
    if(.not.allocated(z_p)    ) allocate(z_p(nVertlevels, nCells)    )
    if(.not.allocated(zmid_p) ) allocate(zmid_p(nVertlevels, nCells) )
    if(.not.allocated(dz_p)   ) allocate(dz_p(nVertlevels, nCells)   )
    if(.not.allocated(t_p)    ) allocate(t_p(nVertlevels, nCells)    )
    if(.not.allocated(th_p)   ) allocate(th_p(nVertlevels, nCells)   )
    if(.not.allocated(al_p)   ) allocate(al_p(nVertlevels, nCells)   )
    if(.not.allocated(rho_p)  ) allocate(rho_p(nVertlevels, nCells)  )
    if(.not.allocated(rh_p)   ) allocate(rh_p(nVertlevels, nCells)   )
    if(.not.allocated(znu_p)  ) allocate(znu_p(nVertlevels, nCells)  )
!print *,'LFR-DBG: Block 3'
    if(.not. allocated(chem_conc_p)) allocate(chem_conc_p(nVertlevels, nCells, nChemSpecies) )
    if(.not. allocated(chem_tend_p)) allocate(chem_tend_p(nVertlevels, nCells, nChemSpecies) )
    if(.not. allocated(chem_tend_dyn_p)) allocate(chem_tend_dyn_p(nVertlevels, nCells, nChemSpecies) )
!print *,'LFR-DBG: Block 4'
    if(.not.allocated(w_p)    ) allocate(w_p(nVertlevels+1, nCells)    )
    if(.not.allocated(pres2_p)) allocate(pres2_p(nVertlevels+1, nCells))
    if(.not.allocated(t2_p)   ) allocate(t2_p(nVertlevels+1, nCells)   )
!print *,'LFR-DBG: Block 5'
    if(.not.allocated(qv_p)   ) allocate(qv_p(nVertlevels, nCells)   )
    if(.not.allocated(qc_p)   ) allocate(qc_p(nVertlevels, nCells)   )
    if(.not.allocated(qr_p)   ) allocate(qr_p(nVertlevels, nCells)   )
    if(.not.allocated(qi_p)   ) allocate(qi_p(nVertlevels, nCells)   )
    if(.not.allocated(qs_p)   ) allocate(qs_p(nVertlevels, nCells)   )
    if(.not.allocated(qg_p)   ) allocate(qg_p(nVertlevels, nCells)   )
!print *,'LFR-DBG: Block 6'
    if(.not.allocated(psfc_hyd_p)  ) allocate(psfc_hyd_p(nCells)          )
    if(.not.allocated(psfc_hydd_p) ) allocate(psfc_hydd_p(nCells)         )
    if(.not.allocated(pres_hyd_p)  ) allocate(pres_hyd_p(nVertlevels, nCells)  )
    if(.not.allocated(pres_hydd_p) ) allocate(pres_hydd_p(nVertlevels, nCells) )
    if(.not.allocated(pres2_hyd_p) ) allocate(pres2_hyd_p(nVertlevels+1, nCells) )
    if(.not.allocated(pres2_hydd_p)) allocate(pres2_hydd_p(nVertlevels+1, nCells))
    if(.not.allocated(znu_hyd_p)   ) allocate(znu_hyd_p(nVertlevels, nCells)   )
!print *,'LFR-DBG: Block 7'
    end subroutine allocate_forall_chemistry

    !=================================================================================================================
    subroutine deallocate_forall_chemistry(configs)
    !=================================================================================================================
    type(mpas_pool_type),intent(in):: configs

        if(allocated(xlon_p) ) deallocate(xlon_p)
        if(allocated(xlat_p) ) deallocate(xlat_p)
        if(allocated(lwupb_p)) deallocate(lwupb_p)
        if(allocated(sfc_albedo_p) ) deallocate(sfc_albedo_p)
        if(allocated(psfc_p) ) deallocate(psfc_p)
        if(allocated(ptop_p) ) deallocate(ptop_p)
        if(allocated(coszr_p)      ) deallocate(coszr_p)

        if(allocated(o3_p)   ) deallocate(o3_p)
        if(allocated(u_p)    ) deallocate(u_p)
        if(allocated(v_p)    ) deallocate(v_p)
        if(allocated(fzm_p)  ) deallocate(fzm_p)
        if(allocated(fzp_p)  ) deallocate(fzp_p)
        if(allocated(zz_p)   ) deallocate(zz_p)
        if(allocated(pres_p) ) deallocate(pres_p)
        if(allocated(pi_p)   ) deallocate(pi_p)
        if(allocated(z_p)    ) deallocate(z_p)
        if(allocated(zmid_p) ) deallocate(zmid_p)
        if(allocated(dz_p)   ) deallocate(dz_p)
        if(allocated(t_p)    ) deallocate(t_p)
        if(allocated(th_p)   ) deallocate(th_p)
        if(allocated(al_p)   ) deallocate(al_p)
        if(allocated(rho_p)  ) deallocate(rho_p)
        if(allocated(rh_p)   ) deallocate(rh_p)
        if(allocated(znu_hyd_p)  ) deallocate(znu_hyd_p)

        if(allocated(chem_conc_p)) deallocate(chem_conc_p)
        if(allocated(chem_tend_p)) deallocate(chem_tend_p)
        if(allocated(chem_tend_dyn_p)) deallocate(chem_tend_dyn_p)

        if(allocated(w_p)    ) deallocate(w_p)
        if(allocated(pres2_p)) deallocate(pres2_p)
        if(allocated(t2_p)   ) deallocate(t2_p)

        if(allocated(qv_p)   ) deallocate(qv_p)
        if(allocated(qc_p)   ) deallocate(qc_p)
        if(allocated(qr_p)   ) deallocate(qr_p)
        if(allocated(qi_p)   ) deallocate(qi_p)
        if(allocated(qs_p)   ) deallocate(qs_p)
        if(allocated(qg_p)   ) deallocate(qg_p)

        if(allocated(psfc_hyd_p)  ) deallocate(psfc_hyd_p)
        if(allocated(psfc_hydd_p) ) deallocate(psfc_hydd_p)
        if(allocated(pres_hyd_p)  ) deallocate(pres_hyd_p)
        if(allocated(pres_hydd_p) ) deallocate(pres_hydd_p)
        if(allocated(pres2_hyd_p) ) deallocate(pres2_hyd_p)
        if(allocated(pres2_hydd_p)) deallocate(pres2_hydd_p)
        if(allocated(znu_hyd_p)   ) deallocate(znu_hyd_p)
    ! (inalterado, omitido por brevidade)
    end subroutine deallocate_forall_chemistry

    !=================================================================================================================
    subroutine MPAS_to_chemistry(configs,mesh,state,time_lev,diag,diag_physics,nCells, nVertLevels, nChemSpecies)
    !=================================================================================================================
    type(mpas_pool_type),intent(in):: configs
    type(mpas_pool_type),intent(in):: mesh
    type(mpas_pool_type),intent(in):: state
    type(mpas_pool_type),intent(in):: diag
    integer,intent(in):: nCells, nVertLevels, nChemSpecies
    integer,intent(in):: time_lev
    type(mpas_pool_type),intent(inout):: diag_physics

    integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
    logical,pointer:: config_o3climatology
    real(kind=RKIND),dimension(:),pointer    :: latCell,coszr,loncell
    real(kind=RKIND),dimension(:),pointer    :: fzm,fzp,rdzw, lwupb, sfc_albedo
    real(kind=RKIND),dimension(:),pointer    :: surface_pressure,plrad,plradm1
    real(kind=RKIND),dimension(:,:),pointer  :: zgrid, o3clim, o3
    real(kind=RKIND),dimension(:,:),pointer  :: zz,exner,pressure_b,rtheta_p,rtheta_b
    real(kind=RKIND),dimension(:,:),pointer  :: rho_zz,theta_m,pressure_p,u,v,w
    real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qr,qi,qs,qg
    real(kind=RKIND),dimension(:,:,:),pointer:: scalars, chem_conc, chem_tend, chem_tend_dyn

    integer:: i,k,n,spc
    real(kind=RKIND):: z0,z1,z2,w1,w2
    real(kind=RKIND):: rho_a,rho1,rho2,tem1,tem2
    character(len=14) :: idxChemStr
    integer, pointer :: idxTmp

    !initialization:
    call mpas_pool_get_config(configs,'config_o3climatology' ,config_o3climatology)

    call mpas_pool_get_array(mesh,'latCell',latCell)
    call mpas_pool_get_array(mesh,'lonCell',lonCell)
    call mpas_pool_get_array(mesh,'fzm'    ,fzm    )
    call mpas_pool_get_array(mesh,'fzp'    ,fzp    )
    call mpas_pool_get_array(mesh,'rdzw'   ,rdzw   )
    call mpas_pool_get_array(mesh,'zgrid'  ,zgrid  )
    call mpas_pool_get_array(mesh,'zz'     ,zz     )
    call mpas_pool_get_array(diag,'surface_pressure'      ,surface_pressure)
    call mpas_pool_get_array(diag,'exner'                 ,exner           )
    call mpas_pool_get_array(diag,'pressure_base'         ,pressure_b      )
    call mpas_pool_get_array(diag,'pressure_p'            ,pressure_p      )
    call mpas_pool_get_array(diag,'rtheta_base'           ,rtheta_b        )
    call mpas_pool_get_array(diag,'rtheta_p'              ,rtheta_p        )
    call mpas_pool_get_array(diag,'uReconstructZonal'     ,u               )
    call mpas_pool_get_array(diag,'uReconstructMeridional',v               )
    call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
    call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
    call mpas_pool_get_array(state,'w'      ,w      ,time_lev)
    call mpas_pool_get_array(mesh,'latCell',latCell)
    call mpas_pool_get_array(diag_physics,'lwupb' ,lwupb )
    call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
    call mpas_pool_get_array(diag_physics,'coszr'     ,coszr     )
    call mpas_pool_get_dimension(state,'index_qv',index_qv)

    call mpas_pool_get_array(state,'scalars',scalars,time_lev)
    qv => scalars(index_qv,:,:)

    call mpas_pool_get_array(diag_physics,'plrad',plrad)
    call mpas_pool_get_array(diag_physics,'o3clim'    ,o3clim    )
    call mpas_pool_get_array(diag,'o3'    ,o3    )
    
    do spc = 1, nSpecies
       write(idxChemStr, fmt='("index_",A)') trim(to_lower(spc_name(spc)))
       !print *,'idxChemStr = ',trim(idxChemStr)
       nullify(idxTmp)
       call mpas_pool_get_dimension(state, trim(idxChemStr), idxTmp)
    
       if (associated(idxTmp)) then
          idxChem(spc) = idxTmp
       else
          call mpas_log_write('indice nao encontrado para '//trim(spc_name(spc)), &
               messageType=MPAS_LOG_WARN)
          idxChem(spc) = -1   ! sentinela, pra facilitar debug depois
       end if
        end do

    do i = 1,nCells
        do k = 1, nVertLevels
            do spc = 1, nspecies
                chem_g(spc)%sc_p(k,i) = scalars(idxChem(spc), k, i)
            end do
        end do
    end do

    do i = 1,nCells
        do k = 1, nVertLevels
            qv_p(k,i) = max(0.,qv(k,i))

            u_p(k,i) = u(k,i)
            v_p(k,i) = v(k,i)

            zz_p(k,i)  = zz(k,i)
            rho_p(k,i) = zz(k,i) * rho_zz(k,i)
            rho_p(k,i) = rho_p(k,i)*(1._RKIND + qv_p(k,i))
            th_p(k,i)  = theta_m(k,i) / (1._RKIND + R_v/R_d * qv_p(k,i))
            t_p(k,i)   = th_p(k,i)*exner(k,i)

            pi_p(k,i)   = exner(k,i)
            pres_p(k,i) = pressure_p(k,i) + pressure_b(k,i)

            zmid_p(k,i) = 0.5*(zgrid(k+1,i)+zgrid(k,i))
            dz_p(k,i)   = zgrid(k+1,i)-zgrid(k,i)

        end do
    end do
 !print *, 'LFR-DBG: MPAS_to_chemistry: after filling _p arrays'; call flush(6)
    do i = 1,nCells
        xlat_p(i)       = latCell(i) / degrad
        xlon_p(i)       = lonCell(i) / degrad
        lwupb_p(i)      = lwupb(i)
        sfc_albedo_p(i) = sfc_albedo(i)
        coszr_p(i)      = coszr(i)
    end do

    !print *, 'LFR-DBG: MPAS_to_chemistry: calculating surface pressure (hydrostatic)'; call flush(6)
    do i = 1,nCells
        tem1 = zgrid(2,i)-zgrid(1,i)
        tem2 = zgrid(3,i)-zgrid(2,i)
        rho1 = rho_zz(1,i) * zz(1,i) * (1. + qv_p(1,i))
        rho2 = rho_zz(2,i) * zz(2,i) * (1. + qv_p(2,i))
        surface_pressure(i) = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
                            * (rho1 + 0.5*(rho2-rho1)*tem1/(tem1+tem2))
        surface_pressure(i) = 0.5*gravity*(zgrid(2,i)-zgrid(1,i)) &
                            * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
        surface_pressure(i) = surface_pressure(i) + pressure_p(1,i) + pressure_b(1,i)
    end do

    do i = 1,nCells
        do k = 1, nVertLevels
            znu_p(k,i) = pres_p(k,i) / surface_pressure(i)
        end do
    end do

    !print *, 'LFR-DBG: MPAS_to_chemistry: Arrays em níveis w (kts:kte+1)'; call flush(6)
    do i = 1,nCells
        do k = 1, nVertLevels  !+1
            w_p(k,i) = w(k,i)
            z_p(k,i) = zgrid(k,i)
        end do
    end do

    !print *, 'LFR-DBG: MPAS_to_chemistry: Interpolação de pressão e temperatura para níveis w'; call flush(6)
    do i = 1,nCells
        do k = 2, nVertLevels
            tem1 = 1./(zgrid(k+1,i)-zgrid(k-1,i))
            fzm_p(k,i) = (zgrid(k,i)-zgrid(k-1,i)) * tem1
            fzp_p(k,i) = (zgrid(k+1,i)-zgrid(k,i)) * tem1
            t2_p(k,i)    = fzm_p(k,i)*t_p(k,i) + fzp_p(k,i)*t_p(k-1,i)
            pres2_p(k,i) = fzm_p(k,i)*pres_p(k,i) + fzp_p(k,i)*pres_p(k-1,i)
        end do
    end do

    !print *, 'LFR-DBG: MPAS_to_chemistry: Topo do modelo'; call flush(6)
    k = nVertLevels+1
    do i = 1,nCells
        z0 = zgrid(k,i)
        z1 = 0.5*(zgrid(k,i)+zgrid(k-1,i)) 
        z2 = 0.5*(zgrid(k-1,i)+zgrid(k-2,i))
        w1 = (z0-z2)/(z1-z2)
        w2 = 1.-w1
        t2_p(k,i) = w1*t_p(k-1,i) + w2*t_p(k-2,i)
        pres2_p(k,i) = exp(w1*log(pres_p(k-1,i))+w2*log(pres_p(k-2,i)))
    end do

    !print *, 'LFR-DBG: MPAS_to_chemistry: Extrapolação para superfície'; call flush(6)
    k = nVertLevels-1
    do i = 1,nCells
        z0 = zgrid(k,i)
        z1 = 0.5*(zgrid(k,i)+zgrid(k+1,i)) 
        z2 = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
        w1 = (z0-z2)/(z1-z2)
        w2 = 1.-w1
        t2_p(k,i)    = w1*t_p(k,i)+w2*t_p(k+1,i)
        pres2_p(k,i) = w1*pres_p(k,i)+w2*pres_p(k+1,i)
        psfc_p(i) = pres2_p(k,i)
        psfc_p(i) = surface_pressure(i)
    end do

    !print *, 'LFR-DBG: MPAS_to_chemistry: Calculating hydrostatic pressure'; call flush(6)
    do i = 1,nCells
        k = nVertLevels+1
        pres2_hyd_p(k,i)  = pres2_p(k,i)
        pres2_hydd_p(k,i) = pres2_p(k,i)
        do k = nVertLevels, 1, -1
            rho_a = rho_p(k,i) / (1.+qv_p(k,i))
            pres2_hyd_p(k,i)  = pres2_hyd_p(k+1,i)  + gravity*rho_p(k,i)*dz_p(k,i)
            pres2_hydd_p(k,i) = pres2_hydd_p(k+1,i) + gravity*rho_a*dz_p(k,i)
        end do
        do k = nVertLevels, 1, -1
            pres_hyd_p(k,i)  = 0.5*(pres2_hyd_p(k+1,i)+pres2_hyd_p(k,i))
            pres_hydd_p(k,i) = 0.5*(pres2_hydd_p(k+1,i)+pres2_hydd_p(k,i))
        end do
        psfc_hyd_p(i) = pres2_hyd_p(1,i)
        psfc_hydd_p(i) = pres2_hydd_p(1,i)
        do k = nVertLevels, 1, -1
            znu_hyd_p(k,i) = pres_hyd_p(k,i) / psfc_hyd_p(i) 
        end do
    end do

    ! Salvar pressão no topo
    !print *, 'LFR-DBG: MPAS_to_chemistry: Saving pressure at the top'; call flush(6)
    do i = 1,nCells
        plrad(i) = pres2_p(nVertLevels+1,i) 
    end do

    end subroutine MPAS_to_chemistry


     !=================================================================================================================
     subroutine chemistry_to_MPAS(configs,state,time_lev,nCells, nVertLevels)
     !=================================================================================================================
     integer,intent(in):: nVertLevels, nCells
     type(mpas_pool_type),intent(in):: configs
     type(mpas_pool_type),intent(in):: state
     integer,intent(in):: time_lev

     integer :: i,k,spc
     real(kind=RKIND),dimension(:,:,:),pointer:: scalars

     call mpas_pool_get_array(state,'scalars',scalars,time_lev)

     do i = 1,nCells
        do k = 1, nVertLevels
            do spc = 1, nspecies
                scalars(idxChem(spc), k, i) = chem_g(spc)%sc_p(k,i)
            end do
        end do
     end do

   end subroutine chemistry_to_MPAS

    !===========================================================================
    ! zero_chemistry_scalars.f90
    !
    ! Gerado automaticamente por gen_registry_chem.py a partir de chem_list.f90.
    ! Zera todas as 47 especies quimicas (SPACK/BRAMS) dentro do
    ! var_array "scalars", usando os indices index_<especie> registrados pelo
    ! Registry (Registry_chem_scalars.inc.xml). Nao mexe em qv/qc/qr e demais
    ! constituintes que nao sao de quimica.
    !
    ! Uso tipico: chamar uma unica vez, no setup do bloco (antes do primeiro
    ! timestep), para garantir que os campos nasçam com 0.0, ja que o Registry
    ! nao suporta default_value para constituintes de var_array.
    !===========================================================================
    subroutine zero_chemistry_scalars(state)
    
       use mpas_derived_types, only : mpas_pool_type
       use mpas_pool_routines, only : mpas_pool_get_array, mpas_pool_get_dimension
       use mpas_kind_types, only : RKIND
       use mpas_log, only : mpas_log_write
       use mpas_derived_types, only : MPAS_LOG_WARN
    
       implicit none
    
       type (mpas_pool_type), intent(inout) :: state
    
       real (kind=RKIND), dimension(:,:,:), pointer :: scalars
       integer, pointer :: idx
       integer :: n, timeLevel
       integer, parameter :: num_species = 47
    
       character(len=8), dimension(num_species) :: species_names = [ &
            'o3      ', 'h2o2    ', 'no      ', 'no2     ', &
            'no3     ', 'n2o5    ', 'hono    ', 'hno3    ', &
            'hno4    ', 'so2     ', 'sulf    ', 'co      ', &
            'co2     ', 'n2      ', 'o2      ', 'h2o     ', &
            'h2      ', 'o3p     ', 'o1d     ', 'ho      ', &
            'ho2     ', 'ch4     ', 'eth     ', 'alka    ', &
            'alke    ', 'bio     ', 'aro     ', 'hcho    ', &
            'ald     ', 'ket     ', 'crbo    ', 'onit    ', &
            'pan     ', 'op1     ', 'op2     ', 'ora1    ', &
            'ora2    ', 'mo2     ', 'akap    ', 'akep    ', &
            'biop    ', 'pho     ', 'add     ', 'arop    ', &
            'cbop    ', 'oln     ', 'xo2     ' ]
    
       ! Zera para os time levels 1 e 2 (os dois niveis de tempo do "state" no
       ! nucleo dinamico MPAS-A padrao). Ajuste se o seu build usar outro numero.
       do timeLevel = 1, 2
    
          nullify(scalars)
          call mpas_pool_get_array(state, 'scalars', scalars, timeLevel)
    
          if (.not. associated(scalars)) then
             call mpas_log_write('zero_chemistry_scalars: scalars nao associado para timeLevel=$i', &
                  messageType=MPAS_LOG_WARN, intArgs=[timeLevel])
             cycle
          end if
    
          do n = 1, num_species
             nullify(idx)
             call mpas_pool_get_dimension(state, 'index_'//trim(species_names(n)), idx)
    
             if (associated(idx)) then
                scalars(idx, :, :) = 0.0_RKIND
             else
                call mpas_log_write('zero_chemistry_scalars: indice nao encontrado para especie ('// &
                     trim(species_names(n))//'). Pacote chemistry ativo?', &
                     messageType=MPAS_LOG_WARN)
             end if
          end do
    
       end do
    
       call mpas_log_write('zero_chemistry_scalars: $i especies quimicas zeradas.', &
            intArgs=[num_species])
    
    end subroutine zero_chemistry_scalars


end module monan_chemistry_interface