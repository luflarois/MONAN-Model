module mem_spack

    use chem_list, only: &
    nspecies, & ! parameter
    nr, & ! parameter
    nr_photo, & ! parameter
    photojmethod       ! parameter

    implicit none

    type, public :: spack_type

        !3d real
        double precision, allocatable, dimension(:, :, :) :: dldrdc

        !2d real
        double precision, allocatable, dimension(:, :) :: sc_p_new
        double precision, allocatable, dimension(:, :) :: sc_p_4
        double precision, allocatable, dimension(:, :) :: dlr
        double precision, allocatable, dimension(:, :) :: dlr3


        double precision, allocatable, dimension(:, :) :: jphoto
        double precision, allocatable, dimension(:, :) :: rk
        double precision, allocatable, dimension(:, :) :: w
        double precision, allocatable, dimension(:, :) :: sc_p

        !1d real
        double precision, allocatable, dimension(:) :: temp
        double precision, allocatable, dimension(:) :: press
        double precision, allocatable, dimension(:) :: cosz
        double precision, allocatable, dimension(:) :: att
        double precision, allocatable, dimension(:) :: vapp
        double precision, allocatable, dimension(:) :: volmol
        double precision, allocatable, dimension(:) :: volmol_i
        double precision, allocatable, dimension(:) :: xlw
        double precision, allocatable, dimension(:) :: err

    end type spack_type

    type, public :: spack_type_2d
        double precision, allocatable, dimension(:, :) :: dlmat
        double precision, allocatable, dimension(:) :: dlb1
        double precision, allocatable, dimension(:) :: dlb2
        double precision, allocatable, dimension(:) :: dlb3
        double precision, allocatable, dimension(:) :: dlb4
        double precision, allocatable, dimension(:) :: dlk1
        double precision, allocatable, dimension(:) :: dlk2
        double precision, allocatable, dimension(:) :: dlk3
        double precision, allocatable, dimension(:) :: dlk4

    end type spack_type_2d

    private

    double precision, parameter :: rtols = 1.d-3 ! 1e-2 means two digits
    double precision, parameter :: atols = 1.d+7 ! jacobson (1998, smvgear) range 1.e3-1.e7! 1.d0

    double precision, public, dimension(nspecies) :: atol, rtol
    type(spack_type), public, allocatable, dimension(:) :: spack
    type(spack_type_2d), public, allocatable, dimension(:, :) :: spack_2d
    logical, public :: spack_alloc = .false.

    public :: alloc_spack


    contains
    !========================================================================

    subroutine alloc_spack(chemistry, maxblock_size)

        implicit none
        !integer, intent(in) :: nob_mem
        integer i, ii, nob, n
        integer, intent(in) :: chemistry
        integer, intent(in) :: maxblock_size

        if (spack_alloc) then
            print *, 'error: spack_alloc already allocated'
            print *, 'routine: spack_alloc file: mem_spack.f90'
            stop
        end if

        !  if(nob_mem==0) nob=nob_real
        !  if(nob_mem==1) nob=1


        do n = 1, nspecies
            atol(n) = atols
            rtol(n) = rtols
        end do

        !- allocating spaces to copy structure

        nob = 1 ! to save memory, the scratch will be re-used
        allocate(spack(nob))

        !- maxblock_size is the maximum block size including all grids
        !- all grids/blocks will share the same array/memory area
    

        do i = 1, nob

            !- 3d variables
            allocate(spack(i)%dldrdc (1:maxblock_size, nspecies, nspecies)) ;
            spack(i)%dldrdc = 0.0d0

            !- 2d variables
            allocate(spack(i)%jphoto (1:maxblock_size, nr_photo)) ;
            spack(i)%jphoto = 0.0d0

            allocate(spack(i)%rk (1:maxblock_size, nr)) ;
            spack(i)%rk = 0.0d0
            allocate(spack(i)%w (1:maxblock_size, nr)) ;
            spack(i)%w = 0.0d0
            allocate(spack(i)%sc_p (1:maxblock_size, nspecies)) ;
            spack(i)%sc_p = 0.0d0
            allocate(spack(i)%sc_p_new(1:maxblock_size, nspecies)) ;
            spack(i)%sc_p_new = 0.0d0

            allocate(spack(i)%dlr (1:maxblock_size, nspecies)) ;
            spack(i)%dlr = 0.0d0

            !- for rodas 3 only for version 1
            !if( chemistry == 4) then
            !  allocate(spack(i)%dlr3  (1:maxblock_size,nspecies))    ;spack(i)%dlr3    = 0.0d0
            !  allocate(spack(i)%sc_p_4 (1:maxblock_size,nspecies))   ;spack(i)%sc_p_4  = 0.0d0
            !endif

            !- 1d variables
            allocate(spack(i)%temp (1:maxblock_size)) ;
            spack(i)%temp = 0.0d0
            allocate(spack(i)%press (1:maxblock_size)) ;
            spack(i)%press = 0.0d0
            allocate(spack(i)%cosz (1:maxblock_size)) ;
            spack(i)%cosz = 0.0d0
            allocate(spack(i)%att (1:maxblock_size)) ;
            spack(i)%att = 0.0d0
            allocate(spack(i)%vapp (1:maxblock_size)) ;
            spack(i)%vapp = 0.0d0
            allocate(spack(i)%volmol (1:maxblock_size)) ;
            spack(i)%volmol = 0.0d0
            allocate(spack(i)%volmol_i (1:maxblock_size)) ;
            spack(i)%volmol_i = 0.0d0
            allocate(spack(i)%xlw (1:maxblock_size)) ;
            spack(i)%xlw = 0.0d0
            allocate(spack(i)%err (1:maxblock_size)) ;
            spack(i)%err = 0.0d0

        enddo

        spack_alloc = .true.

    end subroutine alloc_spack
    !-----------------------------------------------------------------
    !subroutine dealloc_spack(nob_mem)
    !implicit none
    !integer, intent(in) :: nob_mem
    !integer i,ii,nob
    !  !if(nob_mem==0) nob=nob_real
    !  !if(nob_mem==1) nob=1
    !
    !  do i=1,nob_real
    !
    !    !if (associated(spack(i)%dlmat   )) deallocate(spack(i)%dlmat   )
    !    !if (associated(spack(i)%dlmatlu )) deallocate(spack(i)%dlmatlu )
    !    if (associated(spack(i)%dldrdc  )) deallocate(spack(i)%dldrdc  )
    !
    !    !2d variables
    !
    !    if  (associated(spack(i)%jphoto))  deallocate (spack(i)%jphoto)
    !    if  (associated(spack(i)%rk    ))  deallocate (spack(i)%rk    )
    !    if  (associated(spack(i)%w     ))  deallocate (spack(i)%w     )
    !    if  (associated(spack(i)%sc_p  ))  deallocate (spack(i)%sc_p  )
    !
    !    if  (associated(spack(i)%sc_p_new))deallocate (spack(i)%sc_p_new )
    !    if  (associated(spack(i)%dlr     ))deallocate (spack(i)%dlr      )
    !    !if  (associated(spack(i)%dlk1    ))deallocate (spack(i)%dlk1     )
    !    !if  (associated(spack(i)%dlk2    ))deallocate (spack(i)%dlk2     )
    !    !if (associated(spack(i)%dlb1    ))deallocate (spack(i)%dlb1     )
    !    !if (associated(spack(i)%dlb2    ))deallocate (spack(i)%dlb2     )
    !
    !    !1d variables
    !    if  (associated(spack(i)%temp   ))deallocate (spack(i)%temp   )
    !    if  (associated(spack(i)%press  ))deallocate (spack(i)%press  )
    !    if  (associated(spack(i)%cosz   ))deallocate (spack(i)%cosz   )
    !    if  (associated(spack(i)%vapp   ))deallocate (spack(i)%vapp   )
    !    if  (associated(spack(i)%volmol ))deallocate (spack(i)%volmol )
    !    if  (associated(spack(i)%xlw    ))deallocate (spack(i)%xlw    )
    !
    !   do ii=1,maxblock_size
    !    if  (associated(spack_2d(ii,i)%dlmat)) deallocate (spack_2d(ii,i)%dlmat)
    !    if  (associated(spack_2d(ii,i)% dlb1)) deallocate (spack_2d(ii,i)% dlb1)
    !    if  (associated(spack_2d(ii,i)% dlb2)) deallocate (spack_2d(ii,i)% dlb2)
    !   enddo
    !
    !  enddo
    !
    ! !deallocating  spaces to copy structure
    !  if(allocated(spack)) deallocate(spack)
    !
    !  if(allocated(spack_2d)) deallocate(spack_2d)
    !
    !end subroutine dealloc_spack

end module mem_spack
