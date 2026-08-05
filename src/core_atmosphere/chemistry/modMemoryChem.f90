module modMemoryChem
    use chem_list, only: &
         nspecies  &
        ,spc_name  &
        ,spc_alloc &
        ,transport &
        ,on

    implicit none

   type chem_vars   
        real, allocatable :: sc_p(:,:)
        real, allocatable :: sc_t(:,:) 
        real, allocatable :: sc_t_dyn(:,:)
    end type chem_vars
    type (chem_vars)  :: chem_g(nspecies)

    real,allocatable :: last_accepted_dt(:)
    integer, allocatable :: block_end(:)
    
    integer :: nspecies_chem_transported
    integer :: nspecies_chem_no_transported
    integer :: transp_chem_index(nSpecies)
    integer :: no_transp_chem_index(nSpecies)
    integer :: n_dyn_chem
    real    :: chem_timestep 
    character(len=8) :: split_method
    integer :: chemistry
    integer :: maxblock_size

contains

    subroutine alloc_chem(dim1,dim2,dim3)
        implicit none

        integer, intent(in) :: dim1
        integer, intent(in) :: dim2
        integer, intent(in) :: dim3
        integer :: ispc

        do ispc = 1,nSpecies
            if(allocated(chem_g(iSpc)%sc_p)) then
                print *,'Error: specie ',iSpc," "//trim(spc_name(ispc)),', used to sc_p, already allocated!'
                print *,'Please, check it!'
                stop 'ERROR!'
            end if
            allocate(chem_g(iSpc)%sc_p(dim1,dim2))
            chem_g(iSpc)%sc_p = 0.0

            !-srf: only tendencies arrays for transported species are allocated (save memory)
            if (spc_alloc(transport,ispc) == on) then  
                if(allocated(chem_g(iSpc)%sc_t)) then
                    print *,'Error: specie ',iSpc," "//trim(spc_name(ispc)),' used to sc_t, already allocated!'
                    print *,'Please, check it!'
                    stop 'ERROR!'
                end if
                allocate(chem_g(iSpc)%sc_t(dim1,dim2))
            else
            !- for non-transported species, arrays are allocated with one-dimension
                allocate(chem_g(iSpc)%sc_t(1,1))
            end if
            chem_g(iSpc)%sc_t = 0.0

            !- allocate memory for the tendency tracer mixing ratio if parallel spliting 
            !- operator will be used
            if (spc_alloc(transport,ispc) == on .and. trim(adjustl(split_method)) == 'PARALLEL') then  
                if(allocated(chem_g(iSpc)%sc_t_dyn)) then
                    print *,'Error: specie ',iSpc," "//trim(spc_name(ispc)),' used to sc_t_dyn, already allocated!'
                    print *,'Please, check it!'
                    stop 'ERROR!'
                end if
                allocate(chem_g(iSpc)%sc_t_dyn(dim1,dim2))  
            else
                allocate(chem_g(iSpc)%sc_t_dyn(1,1))  
            end if
            chem_g(iSpc)%sc_t_dyn = 0.0                                
        end do

        if(allocated(last_accepted_dt)) then
            print *,'Error: last_accepted_dt, already allocated!'
            print *,'Please, check it!'
            stop 'ERROR!'
        end if
        allocate(last_accepted_dt(dim3), block_end(dim3))

        nspecies_chem_transported = 0
        transp_chem_index   (:)   = 0

        nspecies_chem_no_transported = 0
        no_transp_chem_index(:)      = 0

        do ispc=1,nspecies
           !- Fill pointers to scalar arrays into scalar tables 
           !-srf - only for the "transported" species
           if (spc_alloc(transport,ispc) == on .and. allocated(chem_g(ispc)%sc_t)) then
              !- number of chem transported species
              nspecies_chem_transported = nspecies_chem_transported + 1

              !- mapping between ispc and transported chem species
              transp_chem_index(nspecies_chem_transported) = ispc    
           else
              nspecies_chem_no_transported = nspecies_chem_no_transported + 1
              no_transp_chem_index(nspecies_chem_no_transported) = ispc
           endif
        enddo

    end subroutine alloc_chem

end module modMemoryChem