module modChem
  use mpas_derived_types
  use ModTuv_driver, only : &
      Tuv_driver
  use monan_chemistry_vars
  use monan_chemistry_interface
  use mpas_pool_routines
  use mpas_timekeeping, only: mpas_get_time, mpas_get_clock_time
  use chem_list, only: nspecies, weight, spc_name, nr_photo
  use mod_chem_spack_rodas3_dyndt, only: chem_rodas3_dyndt !, test_rodas3_dynt
  use modMemoryChem, only: alloc_chem &
      , n_dyn_chem &
      , chem_timestep &
      , chemistry &
      , chem_g &
      , last_accepted_dt  &
      , maxblock_size &
      , split_method &
      , nspecies_chem_transported &
      , nspecies_chem_no_transported &
      , transp_chem_index &
      , no_transp_chem_index &
      , block_end
  use chem_list, only: &
      nspecies  &
      ,spc_name  &
      ,spc_alloc &
      ,transport &
      ,on
  use mem_spack, only: spack, spack_alloc, alloc_spack

  !use modFilesChem, only: read_Cams_Chem, mass_frac_to_molec_cm3, read_static_test, readBramsOut

  implicit none
  private
  public :: chemistry_driver

contains        

    subroutine chemistry_driver(domain, iTimestep, currTime)
        !! ## driver para química atmosférica
        !!
        !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
        !! ## MONAN
        !!
        !! Author: rodrigues, L.F.
        !!
        !! E-mail: luflarois@gmail.com
        !!
        !! Date: 2026-03-03
        !!
        !! #####Version: 0.1.0
        !!
        !! —
        !! **Full description**:
        !!
        !! driver para a química atmosférica.  
        !! Este módulo é responsável por orquestrar a execução dos processos químicos, 
        !! incluindo a chamada de rotinas específicas para cada processo (fotólise, química gasosa, 
        !! química de aerossóis, etc.).  Ele também pode ser responsável por gerenciar o acoplamento 
        !! entre os processos químicos e os outros componentes do modelo (dinâmica, física, etc.).
        !!
        !! ** History**:
        !!
        !! - Itenizado_as_alterações_ao_longo_do_tempo (genérica)
        !!—
        !! ** Licence **:
        !!
        !! CC-GPL 3.0 License (https://creativecommons.org/licenses/GPL/3.0/)
        !!
        implicit none

        type(domain_type),intent(inout):: domain
        integer, intent(in) :: iTimestep
        !! Current timestep
        type(MPAS_Time_Type):: currTime
        !!

        !local pointers:
         type(mpas_pool_type),pointer::  configs,      &
                                         mesh,         &
                                         state,        &
                                         diag,         &
                                         diag_physics, &
                                         tend_physics, &
                                         atm_input,    &
                                         sfc_input
        type(block_type),pointer:: block
        
        integer :: mynum, nBlocks, thread, time_lev, i,k,n, ierr
        integer, pointer:: nThreads, nCells
        integer, pointer :: nVertLevels
        real (kind=RKIND), pointer :: config_dt
        real (kind=RKIND) :: conc

        integer,dimension(:),pointer:: cellSolveThreadStart, cellSolveThreadEnd
        real(kind=RKIND),dimension(:,:),pointer  :: o3
        character(len=2) ::  ctime
        character(len=StrKIND) :: timeStamp
        real(kind=RKIND), allocatable :: jphoto(:,:,:)


        block => domain % blocklist
        myNum = domain%dminfo%my_proc_id
print *,'LFR-DBG: Starting chemistry_driver - iTimestep : ', iTimestep
        call mpas_pool_get_subpool(block%structs, 'mesh', mesh)
print *,'LFR-DBG: 01'
        call mpas_pool_get_subpool(block%structs,'diag_physics',diag_physics)
print *,'LFR-DBG: 02'
        call mpas_pool_get_subpool(block%structs, 'state', state)
print *,'LFR-DBG: 03'
        call mpas_pool_get_subpool(block%structs, 'diag', diag)
print *,'LFR-DBG: 04'
        call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
print *,'LFR-DBG: 05'
!        call mpas_pool_get_dimension(mesh,'nSpecies', nSpecies)
print *,'LFR-DBG: 06'
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
print *,'LFR-DBG: 07'
!print *,'LFR-DBG: Retrieved dimensions - nVertLevels = ', nVertLevels, ' nSpecies = ', nSpecies, ' nCells = ', nCells
        call mpas_pool_get_config(block % configs, 'config_dt', config_dt)
print *,'LFR-DBG: 08'

        call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
print *,'LFR-DBG: 09'
!print *,'LFR-DBG: Retrieved config_dt = ', config_dt       
        chem_timestep = config_dt/4.0 !LFR-DBG TBD: Chemistry timestep in seconds (to be set in config file later)
print *,'LFR-DBG: Set chem_timestep = ', chem_timestep

         !- set the number of dynamics cycles inside each chemistry cycle:
         !- observe that 'config_dt' (timestep of grid) is used.
        !- set the number of dynamics cycles inside each chemistry cycle:
        !- observe that 'config_dt' (timestep of grid) is used.
        n_dyn_chem = max(1,nint(chem_timestep/config_dt))
print *,'LFR-DBG: Set n_dyn_chem = ', n_dyn_chem

         !- chemistry is called every 'n_dyn_chem' steps:
         !- observe that 'iTimestep' is the current step of the grid.
        do while(associated(block))
print *,'LFR-DBG: 10'
            call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)
print *,'LFR-DBG: 11'      
            allocate(jphoto(nVertlevels, nCells, nr_photo))
print *,'LFR-DBG: 12'

            if(mod(iTimestep, 4) == 0 .or. iTimestep == 1) then
                !print *, 'Hello from chemistry_driver', iTimestep
print *,'LFR-DBG: 12.1'

                if(iTimestep == 1) then
                    !print *, 'Allocating chemistry arrays...'
                    maxblock_size = nCells
                    call alloc_chem(nVertLevels,nCells,nVertLevels)
                    call alloc_spack(maxblock_size = maxblock_size)
                end if
print *,'LFR-DBG: 12.2'
                !print *, 'Allocating chemistry arrays...'
                call allocate_forall_chemistry(nCells = nCells, nVertLevels = nVertLevels, nChemSpecies = nSpecies)
print *,'LFR-DBG: 12.3'
                !chemistry prep step:
                time_lev = 1
print *,'LFR-DBG: 13'
!$OMP PARALLEL DO
                do thread=1,nThreads
                    !print *,'Fazendo thread ', thread !, cellSolveThreadStart(thread), cellSolveThreadEnd(thread)
                    call MPAS_to_chemistry(block%configs,mesh,state,time_lev,diag,diag_physics, nCells, nVertLevels, nSpecies)
                end do
!$OMP END PARALLEL DO     
print *,'LFR-DBG: 14'
                !call test_rodas3_dynt() !Only to test if works

                !if(iTimestep == 1) then
                !    chem_conc_p = 0.0
                    !call read_static_test(filename = 'data_chem_tst.csv', chem_out = chem_conc_p, nVertLevels = nVertLevels, nCells = nCells, xlat_p = xlat_p, xlon_p = xlon_p, nspecies = nSpecies)
                    !call read_Cams_Chem(filename = 'data_plev.nc',zmid = zmid_p,xlat = xlat_p, xlon = xlon_p &
                    !                   , chem_out = chem_conc_p,nVertLevels = nVertLevels, nCells = nCells, nSpecies = nSpecies)
                !    call readBramsOut(chem_out = chem_conc_p,xlat = xlat_p,xlon = xlon_p,nVertLevels = nVertLevels,nCells = nCells,nspecies_in = nSpecies)

                !end if
print *,'LFR-DBG: 15'
!                print *, 'Calling chemistry driver...'
                call Tuv_driver(           &
                 domain     = domain       &
                ,iTimestep  = iTimestep    &             
                ,press      = pres_hyd_p   &
                ,temp       = t_p          &
                ,zt_        = z_p          &
                ,zm_        = zmid_p       &
                ,dzp        = dz_p         &
                ,rho        = rho_p        &
                ,pp         = pres_p       &
                ,coszr      = coszr_p      &
                ,rlongup    = lwupb_p      &
                ,glat       = xlat_p       &
                ,glon       = xlon_p       &
                ,qv         = qv_p         &
                ,sfc_albedo = sfc_albedo_p &
                ,nCells = nCells           &
                ,nVertLevels = nVertLevels &
                ,mynum = mynum             &
                ,jphoto = jphoto           &
                )

write(ctime,fmt='(I2.2)') iTimestep     
!open(unit=22, file='valJ_'//ctime//'.dat', status='replace', action='write')
!do i=1,nCells
!    write(22,fmt='(F10.5,1X,F10.5,1X,E10.5)') xlat_p(i),xlon_p(i),jphoto(1,i,1)
!end do
!close(unit=22)

                ! call Rodas3_dynt(                &
                !  press         = pres_hyd_p      &
                ! ,temp          = t_p             &
                ! ,rho           = rho_p           &
                ! ,pp            = pres_p          &
                ! ,coszr         = coszr_p         &
                ! ,glat          = xlat_p          &
                ! ,glon          = xlon_p          &
                ! ,qv            = qv_p            &
                ! ,qc            = qc_p            &
                ! ,chem_conc     = chem_conc_p     &
                ! ,chem_tend     = chem_tend_p     &
                ! ,chem_tend_dyn = chem_tend_dyn_p &
                ! ,jphoto        = jphoto          &
                ! ,dtlt          = config_dt       &
                ! ,itimestep     = iTimestep       &
                ! ,nCells        = nCells          &
                ! ,nVertLevels   = nVertLevels     &
                ! ,timestamp     = timeStamp       &
                ! ,mynum         = mynum           &
                ! )


                call chem_rodas3_dyndt( &
                    domain     = domain       &
                  , nob = nVertLevels &
                  , block_end = block_end &
                  , dtlt = config_dt &
                  , press = pres_hyd_p &
                  , temp = t_p &
                  , vapp = qv_p&
                  , last_accepted_dt  = last_accepted_dt &
                  , n_dyn_chem = n_dyn_chem &
                  , split_method = split_method &
                  , jphoto = jphoto &
                  , chem_g = chem_g &
                  , nspecies_chem_transported = nspecies_chem_transported &
                  , nspecies_chem_no_transported = nspecies_chem_no_transported &
                  , transp_chem_index = transp_chem_index &
                  , no_transp_chem_index = no_transp_chem_index &
                  , chemistry = chemistry &
                  , maxblock_size = maxblock_size &
                  )

                call chemistry_to_MPAS(block%configs,diag,nSpecies,nVertLevels,nCells)

                call deallocate_forall_chemistry(block%configs)
    
            end if !End of valid
            
            block => block % next
        end do

    end subroutine chemistry_driver


    subroutine chem_accum(nCells,nVertLevels,ntask,nspecies_chem_transported, & !chem1_g, &
                        transp_chem_index,n_dyn_chem)

        integer         , intent(in) :: ntask 
        !- 1: accumulate tendencies, 2: reset tendencies
        integer         , intent(in) :: nCells
        integer         , intent(in) :: nVertLevels

        ! mem_chem1
        integer          , intent(in)    :: nspecies_chem_transported
!        type (chem1_vars), intent(inout) :: chem1_g(nspecies_chem_transported)
        integer          , intent(inout) :: transp_chem_index(nspecies_chem_transported)
        integer          , intent(in)    :: n_dyn_chem


        !- local var
        integer :: ispc,n,ixyz,ntps,i,j,k
        real n_dyn_chem_i

        if(ntask == 1) then
            n_dyn_chem_i= real(1./n_dyn_chem)
            do ispc=1,nspecies_chem_transported

                !- map the species to transported ones
                n=transp_chem_index(ispc)

                !- calculate the mean dynamic tendency for the entire chemistry timestep
 	            do i= 1,nCells
                    do k=1,nVertLevels
!TBD                        chem1_g(n)%sc_t_dyn(ixyz) = chem1_g(n)%sc_t_dyn(ixyz) + &
!TBD                                              n_dyn_chem_i * chem1_g(n)%sc_t(ixyz)
                    end do
 	            end do
            end do
        else
            do ispc=1,nspecies_chem_transported
                n=transp_chem_index(ispc) !- map the species to transported ones
                !- set to zero arrays for  accumulation over the next chem timestep
!TBD                chem1_g(n)%sc_t_dyn(1:ntps) = 0.
            end do
        endif

    end  subroutine chem_accum


end module modchem