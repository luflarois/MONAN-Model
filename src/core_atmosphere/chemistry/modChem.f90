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
  use modTest, onlY: readMergedChemFile, mapChemToGridNearest
  use mod_chem_output, only: write_chem_netcdf

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
        logical, pointer :: config_chemistry
        real (kind=RKIND), pointer :: config_chem_timestep

        block => domain % blocklist
        myNum = domain%dminfo%my_proc_id


        configs => domain % configs
        call mpas_pool_get_config(configs, 'config_chemistry', config_chemistry)
        if(.not. config_chemistry) return
        call mpas_pool_get_config(configs, 'config_chem_timestep', config_chem_timestep)

        if (myNum == 0) call mpas_log_write(message='Chemistry beggining ...')

        call mpas_pool_get_subpool(block%structs, 'mesh', mesh)
        call mpas_pool_get_subpool(block%structs,'diag_physics',diag_physics)
        call mpas_pool_get_subpool(block%structs, 'state', state)
        call mpas_pool_get_subpool(block%structs, 'diag', diag)
        call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_config(block % configs, 'config_dt', config_dt)

        call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)

         !- set the number of dynamics cycles inside each chemistry cycle:
         !- observe that 'config_dt' (timestep of grid) is used.
        !- set the number of dynamics cycles inside each chemistry cycle:
        !- observe that 'config_dt' (timestep of grid) is used.
        n_dyn_chem = max(1,nint(config_chem_timestep/config_dt))

         !- chemistry is called every 'n_dyn_chem' steps:
         !- observe that 'iTimestep' is the current step of the grid.
        do while(associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nThreads', nThreads)
            allocate(jphoto(nVertlevels, nCells, nr_photo))

            if(mod(iTimestep, n_dyn_chem) == 0 .or. iTimestep == 1) then
                if(iTimestep == 1) then
                    !print *, 'Allocating chemistry arrays...'
                    maxblock_size = nCells
                    call alloc_chem(nVertLevels,nCells,nVertLevels)
                    call alloc_spack(maxblock_size = maxblock_size)
                    !For the first time set last_accept_dt to a fixed value
                    last_accepted_dt = config_dt*n_dyn_chem
                end if
                !print *, 'Allocating chemistry arrays...'
                call allocate_forall_chemistry(nCells = nCells, nVertLevels = nVertLevels, nChemSpecies = nSpecies)
                !chemistry prep step:
                time_lev = 1
!$OMP PARALLEL DO
                do thread=1,nThreads
                    !print *,'Fazendo thread ', thread !, cellSolveThreadStart(thread), cellSolveThreadEnd(thread)
                    call MPAS_to_chemistry(block%configs,mesh,state,time_lev,diag,diag_physics, nCells, nVertLevels, nSpecies)
                end do
!$OMP END PARALLEL DO     

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

!call writeJphoto(iTimestep,nCells,nVertLevels,mynum &
!                 ,xlat_p,xlon_p,coszr_p,jphoto &
!                 ,pres_hyd_p,t_p,z_p,zmid_p,dz_p &
!                 ,rho_p,pres_p,qv_p,sfc_albedo_p,lwupb_p)

!write(ctime,fmt='(I2.2)') iTimestep  
if (iTimestep == 1) call readMergedChemFile("./chem_merged_t60.bin", nSpecies)
if (iTimestep == 1) call mapChemToGridNearest(xlat_p, xlon_p, nCells, nVertLevels &
                                              ,t_p, pres_hyd_p, qv_p, zmid_p, nSpecies)

                call chem_rodas3_dyndt( &
                    nob = nVertLevels &
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

!call write_chem_netcdf( &
!    chem_g = chem_g      &
!  , xlat_p = xlat_p      &
!  , xlon_p = xlon_p      &
!  , nVertLevels = nVertLevels&
!  , nCells = nCells      &
!  , iTimestep = iTimestep   &
!  , dt = 300.          &
!  , start_date = "20220714000000"  &
!    )

!                call chemistry_to_MPAS(block%configs,diag,nSpecies,nVertLevels,nCells)

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


    !==============================================================================
    subroutine writeJphoto(iTimestep,nCells,nVertLevels,mynum &
                 ,glat,glon,coszr,jphoto &
                 ,press,temp,zt_,zm_,dzp &
                 ,rho,pp,qv,sfc_albedo,rlongup)

        !> @brief Brief description
        !>
        !> @project MONAN-Model
        !>
        !> @author luflarois <luflarois@gmail.com>
        !> @date 2026-09-01
        !> @version 0.1.0
        !> @details 
        !> Brief description
        !> 
        !> @license This code is under GPLv3 License, see it in <https://www.gnu.org/licenses/gpl-3.0.en.html>
        !==============================================================================
        implicit none
    
        !Parameters:
        character(len=*), parameter :: source_file    = 'modChem.f90'
        !> name of source file - used for debug
        character(len=*), parameter :: procedure_name = 'writeJphoto'
        !> name of this procedure - used for debug
    
        !Inputs:
        integer, intent(in) :: iTimestep,nCells,nVertLevels,mynum
        real,dimension(nVertLevels,nCells),intent(in) :: press
        !! Pressure field [Pa]
        real,dimension(nVertLevels,nCells),intent(in) :: temp
        !! Temperature field [K]
        real,dimension(nVertLevels,nCells),intent(in) :: zt_
        !! Geometric height of layer top [m]
        real,dimension(nVertLevels,nCells),intent(in) :: zm_
        !! Geometric height of layer midpoint [m]
        real,dimension(nVertLevels,nCells),intent(in) :: dzp
        !! Thickness of each layer [m]
        real,dimension(nVertLevels,nCells),intent(in) :: rho
        !! Air density at layer midpoint [kg/m3]
        real,dimension(nVertLevels,nCells),intent(in) :: pp
        !! Pressure at layer midpoint [Pa]  - SOundings
        real,dimension(nCells),intent(in)         :: coszr
        !! Cosine of solar zenith angle at each cell [unitless]
        real,dimension(nCells),intent(in)         :: rlongup
        !! !all-sky upwelling longwave flux at bottom-of-atmosphere
        real,dimension(nCells),intent(in)         :: glat
        !! latitude, south is negative [degrees]
        real,dimension(nCells),intent(in)         :: glon
        !! longitude, west is negative [degrees]
        real,dimension(nVertLevels,nCells),intent(in) :: qv
        !! water vapor mixing ratio   [kg/kg]
        real,dimension(nCells),intent(in) :: sfc_albedo
        !! Surface albedo
        !!
        real(kind=RKIND), intent(in) :: jphoto(nVertlevels, nCells, nr_photo)
 
    
        integer :: iunit, nl, nc
        character(len=17) :: filename

        !Code area
        if(iTimestep == 1 .or. iTimestep == 144 ) then
            write(filename,fmt='("Jphoto_",I2.2,"_",I3.3,".csv")') mynum,iTimestep
            open(newunit=iunit,file=filename,status='replace',action='write')
            write(iunit,*) "nc,glat,glon,zt,coszr," &
                     //"press, temp," &
                     //"dzp, rho, pp," & 
                     //"rlongup, qv, sfc_albedo," &
                     //"jphoto(min),jphoto(max)"
            do nc = 1,nCells
                do nl = 1,nVertLevels
                    write(iunit,fmt='(I4.4,",",4(F10.4,","),8(E15.6,","),1(E18.6,","),E18.6)') &
                     nc,glat(nc),glon(nc),zt_(nl,nc),coszr(nc), &
                     press(nl,nc), temp(nl,nc), dzp(nl,nc), rho(nl,nc), pp(nl,nc), & 
                     rlongup(nc), qv(nl,nc), sfc_albedo(nc), &
                     minval(jphoto(nl,nc,1:17)),maxval(jphoto(nl,nc,1:17))
                end do
            end do
            close(unit=iunit)
        end if
    
    end subroutine writeJphoto

    !---
end module modchem