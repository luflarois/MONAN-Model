module mem_spack

    ! use spack_utils,only: nob_real=>nob,maxblock_size

    !lfr-monan use spack_utils,only: &
    !lfr-monan      maxblock_size      ! in

    use chem_list, only: &
    nspecies, & ! parameter
    nr, & ! parameter
    nr_photo, & ! parameter
    photojmethod       ! parameter

    use, intrinsic :: iso_c_binding, only: c_double

    implicit none

    type, public :: spack_type

        !3d real
        real(kind=c_double), allocatable, dimension(:, :, :) :: dldrdc

        !2d real
        real(kind=c_double), allocatable, dimension(:, :) :: sc_p_new
        real(kind=c_double), allocatable, dimension(:, :) :: sc_p_4
        real(kind=c_double), allocatable, dimension(:, :) :: dlr
        real(kind=c_double), allocatable, dimension(:, :) :: dlr3


        real(kind=c_double), allocatable, dimension(:, :) :: jphoto
        real(kind=c_double), allocatable, dimension(:, :) :: rk
        real(kind=c_double), allocatable, dimension(:, :) :: w
        real(kind=c_double), allocatable, dimension(:, :) :: sc_p

        !1d real
        real(kind=c_double), allocatable, dimension(:) :: temp
        real(kind=c_double), allocatable, dimension(:) :: press
        real(kind=c_double), allocatable, dimension(:) :: cosz
        real(kind=c_double), allocatable, dimension(:) :: att
        real(kind=c_double), allocatable, dimension(:) :: vapp
        real(kind=c_double), allocatable, dimension(:) :: volmol
        real(kind=c_double), allocatable, dimension(:) :: volmol_i
        real(kind=c_double), allocatable, dimension(:) :: xlw
        real(kind=c_double), allocatable, dimension(:) :: err

    end type spack_type

    private

    real(kind=c_double), parameter :: rtols = 1.e-3_c_double ! 1e-2 means two digits
    real(kind=c_double), parameter :: atols = 1.e+7_c_double ! jacobson (1998, smvgear) range 1.e3-1.e7! 1.d0

    real(kind=c_double), public, dimension(nspecies) :: atol, rtol
    type(spack_type), public :: spack
    logical, public :: spack_alloc = .false.

    public :: alloc_spack


    contains
    !========================================================================

    subroutine alloc_spack(maxblock_size)

        implicit none
        integer n
        integer, intent(in) :: maxblock_size

        if (spack_alloc) then
            print *, 'error: spack_alloc already allocated'
            print *, 'routine: spack_alloc file: mem_spack.f90'
            stop
        end if

        do n = 1, nspecies
            atol(n) = atols
            rtol(n) = rtols
        end do

        !- allocating spaces to copy structure
        !- maxblock_size is the maximum block size including all grids
        !- all grids/blocks will share the same array/memory area
        !- spack is a single scratch structure (not an array): there is only
        !- one block being processed at a time, so no outer dimension is needed

        !- 3d variables
        allocate(spack%dldrdc (1:maxblock_size, nspecies, nspecies)) ;
        spack%dldrdc = 0.0_c_double

        !- 2d variables
        allocate(spack%jphoto (1:maxblock_size, nr_photo)) ;
        spack%jphoto = 0.0_c_double

        allocate(spack%rk (1:maxblock_size, nr)) ;
        spack%rk = 0.0_c_double
        allocate(spack%w (1:maxblock_size, nr)) ;
        spack%w = 0.0_c_double
        allocate(spack%sc_p (1:maxblock_size, nspecies)) ;
        spack%sc_p = 0.0_c_double
        allocate(spack%sc_p_new(1:maxblock_size, nspecies)) ;
        spack%sc_p_new = 0.0_c_double

        allocate(spack%dlr (1:maxblock_size, nspecies)) ;
        spack%dlr = 0.0_c_double

        !- for rodas 3 only for version 1
        !if( chemistry == 4) then
        !  allocate(spack%dlr3  (1:maxblock_size,nspecies))    ;spack%dlr3    = 0.0d0
        !  allocate(spack%sc_p_4 (1:maxblock_size,nspecies))   ;spack%sc_p_4  = 0.0d0
        !endif

        !- 1d variables
        allocate(spack%temp (1:maxblock_size)) ;
        spack%temp = 0.0_c_double
        allocate(spack%press (1:maxblock_size)) ;
        spack%press = 0.0_c_double
        allocate(spack%cosz (1:maxblock_size)) ;
        spack%cosz = 0.0_c_double
        allocate(spack%att (1:maxblock_size)) ;
        spack%att = 0.0_c_double
        allocate(spack%vapp (1:maxblock_size)) ;
        spack%vapp = 0.0_c_double
        allocate(spack%volmol (1:maxblock_size)) ;
        spack%volmol = 0.0_c_double
        allocate(spack%volmol_i (1:maxblock_size)) ;
        spack%volmol_i = 0.0_c_double
        allocate(spack%xlw (1:maxblock_size)) ;
        spack%xlw = 0.0_c_double
        allocate(spack%err (1:maxblock_size)) ;
        spack%err = 0.0_c_double

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
