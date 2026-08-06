module mod_chem_spack_rodas3_dyndt

    use modMemoryChem, only: &
    chem_vars        ! Type

    use mem_spack, only: &
    spack_type, & ! Type
    spack_type_2d, &     ! Type
    spack

    use mod_chem_spack_jacdchemdc, only: &
    jacdchemdc        ! Subroutine

    use mod_chem_spack_kinetic, only: &
    kinetic           ! Subroutine

    use mod_chem_spack_fexchem, only: &
    fexchem           ! Subroutine


    use chem_list, only: &
    nspecies &
    ,nr &
    ,PhotojMethod &
    ,nr_photo &
    ,weight &
    ,maxnspecies


    implicit none

    private


    public chem_rodas3_dyndt


contains

    !========================================================================================
    subroutine chem_rodas3_dyndt( &
        nob &
      , block_end &
      , dtlt &
      , press &
      , temp &
      , vapp &
      , last_accepted_dt &
      , n_dyn_chem &
      , split_method &
      , jphoto &
      , chem_g &
      , nspecies_chem_transported &
      , nspecies_chem_no_transported &
      , transp_chem_index &
      , no_transp_chem_index &
      , chemistry &
      , maxblock_size &
    )

        use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_ptr, c_loc

        implicit none

        integer,             intent(in) :: nob
        integer,             intent(in) :: block_end(:)
        integer,             intent(in) :: n_dyn_chem
        integer,             intent(in) :: chemistry
        integer,             intent(in) :: nspecies_chem_transported
        integer,             intent(in) :: nspecies_chem_no_transported
        integer,             intent(in) :: transp_chem_index(:)
        integer,             intent(in) :: no_transp_chem_index(:)
        integer,             intent(in) :: maxblock_size
        real,                intent(in) :: dtlt
        real,                intent(in) :: press(:,:)
        real,                intent(in) :: temp(:,:)
        real,                intent(in) :: vapp(:,:)
        character(len = 20), intent(in) :: split_method
        real,                intent(in) :: jphoto(:,:,:)
        !
        type(chem_vars),     intent(inout) :: chem_g(:)
        real,                intent(inout) :: last_accepted_dt(:)

        interface
            ! Cria a estrutura da matriz esparsa
            function sfcreate_solve(n, complex_flag, error) result(ptr)
                integer, intent(in) :: n, complex_flag
                integer, intent(out) :: error
                integer(kind=8) :: ptr
            end function sfcreate_solve

            ! Obtém o identificador de um elemento da matriz
            function sfgetelement(mat_id, irow, icol) result(elem)
                integer(kind=8), intent(in) :: mat_id
                integer, intent(in) :: irow, icol
                integer(kind=8) :: elem
            end function sfgetelement

            ! Adiciona um valor real a um elemento da matriz
            subroutine sfadd1real(elem, value)
                integer(kind=8), intent(in) :: elem
                double precision, intent(in) :: value
            end subroutine sfadd1real

            ! Zera todos os elementos da matriz
            subroutine sfzero(mat_id)
                integer(kind=8), intent(in) :: mat_id
            end subroutine sfzero

            ! Fatora a matriz (LU, etc.)
            function sffactor(mat_id) result(ierr)
                integer(kind=8), intent(in) :: mat_id
                integer :: ierr
            end function sffactor

            subroutine sfsolve_c(mat_id, rhs, sol) bind(C, name="spSolve")
                import :: c_int64_t, c_double, c_ptr
                integer(c_int64_t), value, intent(in) :: mat_id
                type(c_ptr), value, intent(in) :: rhs
                type(c_ptr), value, intent(in) :: sol
            end subroutine sfsolve_c

        end interface

        real, parameter :: cp = 1004.
        real, parameter :: rgas = 287.
        real, parameter :: cpor = cp / rgas
        real, parameter :: p00 = 1.0e5
        double precision, parameter :: pmar = 28.96d0
        double precision, parameter :: threshold1 = 0.0d0
        double precision, parameter :: threshold2 = -1.d-1
        double precision, parameter :: igamma = 0.5d0
        double precision, parameter :: c43 = 4.d0 / 3.d0
        double precision, parameter :: c83 = 8.d0 / 3.d0
        double precision, parameter :: c56 = 5.d0 / 6.d0
        double precision, parameter :: c16 = 1.d0 / 6.d0
        double precision, parameter :: c112 = 1.d0 / 12.d0
        !- default number of allocatable blocks (re-use scratch arrays spack()% )
        integer, parameter :: inob = 1
        !- parameters for dynamic timestep control
        double precision, parameter :: facmin = 0.2d0 ! lower bound on step decrease factor (default=0.2)
        double precision, parameter :: facmax = 6.0d0 ! upper bound on step increase factor (default=6)
        double precision, parameter :: facrej = 0.1d0 ! step decrease factor after multiple rejections
        double precision, parameter :: facsafe = 0.9d0 ! step by which the new step is slightly smaller
        ! than the predicted value  (default=0.9)
        double precision, parameter :: uround = 1.d-15, elo = 3.0d0, roundoff = 1.d-8
        integer, parameter :: complex_t = 0
        double precision, parameter ::  rtols=1.D-3 ! 1e-2 means two digits
        double precision, parameter ::  atols=1.D+7 ! Jacobson (1998, SMVGEAR) range 1.e3-1.e7! 1.D0


        double precision, allocatable :: dlmat(:, :, :)
        integer, allocatable :: ipos(:)
        integer, allocatable :: jpos(:)

        integer :: numberofnonzeros, nz
        integer :: offset, offsetdlmat, offsetnz
        integer :: blocksize, sizeofmatrix, maxnonzeros, blocknonzeros
        real :: start, finish
        real :: elapsed_time_solver, elapsed_time, elapsed_time_alloc, elapsed_time_dealloc, elapsed_time_copy

        integer(kind = 8) :: matrix_id
        integer :: error
        integer(kind = 8), allocatable, dimension(:) :: element

        integer :: all_accepted
        double precision :: dt_min, dt_max, dt_actual, dt_new
        double precision :: time_f, time_c
        double precision fac, tol, err1, max_err1
        double precision dble_dtlt_i, fxc, igamma_dtstep, dt_chem, dt_chem_i
        integer :: i, ijk, n, j, k, ispc, ji, jj, k_, i_, j_, kij_, kij, ii
        double precision :: atol(nspecies)
        double precision :: rtol(nspecies)

        !integer :: maxblock_size
        double precision, allocatable, target :: rhs_tmp(:), sol_tmp(:)

        type(spack_type_2d), allocatable, dimension(:, :) :: spack_2d
        

        call get_number_nonzeros(nr_photo, nr, nspecies, maxnonzeros)

        sizeofmatrix = nspecies
        !maxblock_size = nvertlevels
        atol = atols
        rtol = rtols

        allocate(ipos(maxnonzeros)) ;
        ipos = 1
        allocate(jpos(maxnonzeros)) ;
        jpos = 1
print *,'LFR-DBG: maxnonzeros, nob, maxblock_size: ',maxnonzeros, nob, maxblock_size
        allocate(spack_2d(maxblock_size, nob))
        do i = 1, nob
            do ii = 1, maxblock_size
            !print *,'*DBG: ',i,ii,nspecies
                allocate(spack_2d(ii, i)%dlmat(nspecies, nspecies)) ;
                spack_2d(ii, i)%dlmat = 0.0d0
                allocate(spack_2d(ii, i)% dlb1(nspecies)) ;
                spack_2d(ii, i)% dlb1 = 0.0d0
                allocate(spack_2d(ii, i)% dlb2(nspecies)) ;
                spack_2d(ii, i)% dlb2 = 0.0d0
                allocate(spack_2d(ii, i)% dlk1(nspecies)) ;
                spack_2d(ii, i)% dlk1 = 0.0d0
                allocate(spack_2d(ii, i)% dlk2(nspecies)) ;
                spack_2d(ii, i)% dlk2 = 0.0d0
                !- for rodas 3 only
                if (chemistry == 4) then
                    allocate(spack_2d(ii, i)% dlb3(nspecies)) ;
                    spack_2d(ii, i)% dlb3 = 0.0d0
                    allocate(spack_2d(ii, i)% dlb4(nspecies)) ;
                    spack_2d(ii, i)% dlb4 = 0.0d0
                    allocate(spack_2d(ii, i)% dlk3(nspecies)) ;
                    spack_2d(ii, i)% dlk3 = 0.0d0
                    allocate(spack_2d(ii, i)% dlk4(nspecies)) ;
                    spack_2d(ii, i)% dlk4 = 0.0d0
                endif

            enddo
        end do


        do i = 1, nob !- loop over all blocks i) = nVertLevels
            !- copying structure from input to internal
            do ijk = 1, block_end(i) !index_g%block_end(i) - MONAN: the block_end is the number of levels
                spack(inob)%press(ijk) = press(ijk,i)
                spack(inob)%temp(ijk) = temp(ijk,i)
                spack(inob)%vapp(ijk) = vapp(ijk,i)
                spack(inob)%volmol(ijk) = (6.02d23 * 1d-15 * pmar) * (press(ijk,i)) / (8.314d0 * temp(ijk,i))
                spack(inob)%volmol_i(ijk) = 1.0d0 / spack(inob)%volmol(ijk)

                !- no transported species section
                do ispc = 1, nspecies_chem_no_transported
                    !- map the species to NO transported ones
                    n = no_transp_chem_index(ispc)
                    !- initialize no-transported species (don't need to convert, because these
                    !- species are already saved using molecule/cm^3 units)
                    spack(inob)%sc_p(ijk, n) = chem_g(n)%sc_p(ijk, i)
                end do

            end do
            !- convert from brams chem (ppbm) arrays to spack (molec/cm3)
            !- transported species section
            if (split_method == 'PARALLEL' .and. n_dyn_chem > 1) then
                do ijk = 1, block_end(i) !index_g%block_end(i)

                    do ispc = 1, nspecies_chem_transported

                        !- map the species to transported ones
                        n = transp_chem_index(ispc)

                        spack(inob)%sc_p(ijk, n) = (chem_g(n)%sc_p(ijk, i) - &  ! updated mixing ratio
                        chem_g(n)%sc_t_dyn(ijk,i) * n_dyn_chem * dtlt)&  ! accumulated tendency
                        * spack(inob)%volmol(ijk) / weight(n)

                    end do
                end do
            else
                do ijk = 1, block_end(i) !index_g%block_end(i)

                    do ispc = 1, nspecies_chem_transported

                        !- map the species to transported ones
                        n = transp_chem_index(ispc)

                        !- conversion from ppbm to molecule/cm^3
                        spack(inob)%sc_p(ijk, n) = chem_g(n)%sc_p(ijk, i) * spack(inob)%volmol(ijk) / weight(n)
                        spack(inob)%sc_p(ijk, n) = max(0.d0, spack(inob)%sc_p(ijk, n))
                    end do

                end do
            endif

            !- Photolysis section
            if (trim(photojmethod) == 'FAST-JX' .or. trim(photojmethod) == 'FAST-TUV') then

                do ijk = 1, block_end(i) !index_g%block_end(i)
                    do n = 1, nr_photo
                        !LFR-MONAN spack(inob)%jphoto(ijk, n) = jphoto(n, k_, i_, j_) !fast_JX_g%jphoto(n,k_,i_,j_)
                        spack(inob)%jphoto(ijk, n) = jphoto(n, ijk, i)
                    end do
                enddo

            elseif(trim(photojmethod) == 'LUT') then

                do ijk = 1, block_end(i) !index_g%block_end(i)

                    !- UV  attenuation un function AOT
                    !spack(inob)%att(ijk) = att(i)

                    !- get zenital angle (for LUT PhotojMethod)
                    !spack(inob)%cosz(ijk) = cosz(i)

                enddo
            endif


            !- compute kinetical and photochemical reactions (array: spack(inob)%rk)
            call kinetic(nr_photo, spack(inob)%jphoto &
            , spack(inob)%rk &
            , spack(inob)%temp &
            , spack(inob)%vapp &
            , spack(inob)%press &
            , 1 &
            , block_end(i) & !index_g%block_end(i), &
            , maxblock_size, nr)

            !- ROSENBROCK METHOD ----------------------------------------------------------------------
            !- kml: Inicio do equivalente a roschem
            !- srf: extending to RODAS 3 (4 stages, order 3)

            dt_chem = last_accepted_dt(i) !index_g%last_accepted_dt(i) ! dble(dtlt)

            dt_min = max(1.0d0, 1.d-2 * dble(dtlt * n_dyn_chem))
            dt_max = dble(dtlt * n_dyn_chem)
            dt_new = 0.0d0
            !print *,'LFR-DBG: initial dt_chem,dt_min,dt_max,dt_new, time_c,time_f: ',i,dt_chem,dt_min, dt_max, dt_new, time_c, time_f
            time_c = 0.0d0
            time_f = dble(dtlt * n_dyn_chem)

            run_until_integr_ends: do while (time_c + roundoff < time_f)

                !-   Compute the Jacobian (DLRDC).
                call jacdchemdc (spack(inob)%sc_p &
                !     CALL jacdchemdc (spack(inob)%sc_p_4  &
                , spack(inob)%rk &
                , spack(inob)%dldrdc& ! Jacobian matrix
                , nspecies, 1 &
                , block_end(i) & !index_g%block_end(i), &
                , maxblock_size, nr)

                !-   Compute chemical net production terms (actually, P-L) at initial time (array spack(inob)%DLr)
                call fexchem (spack(inob)%sc_p &
                !     CALL fexchem (spack(inob)%sc_p_4 &
                , spack(inob)%rk &
                , spack(inob)%dlr & !production term
                , nspecies, 1 &
                , block_end(i) & !index_g%block_end(i), &
                , maxblock_size, nr)

                do ji = 1, nspecies
                    do ijk = 1, block_end(i) !index_g%block_end(i)
                        !                  PRINT *,'LFR-DBG: spack(inob)%dlr(ijk,ji)',ijk,ji,spack(inob)%dlr(ijk, ji)
                        spack_2d(ijk, inob)%dlb1(ji) = spack(inob)%dlr(ijk, ji)
                    enddo
                enddo

                !-   compute matrix (1/(Igamma*dt)  - Jacobian)  where Jacobian = DLdrdc
                !-   fill DLMAT with non-diagonal (fixed in the timestep) Jacobian
                do jj = 1, nspecies
                    do ji = 1, nspecies
                        do ijk = 1, block_end(i) !index_g%block_end(i)

                            spack_2d(ijk, inob)%dlmat(ji, jj) = - spack(inob)%dldrdc(ijk, ji, jj)

                        enddo
                    enddo
                enddo

                untilaccepted: do

                    !- fill DLMAT with diagonal Jacobian (which changes in time, because of dt_chem)
                    igamma_dtstep = 1.0d0 / (igamma * dt_chem)

                    do jj = 1, nspecies
                        do ijk = 1, block_end(i) !index_g%block_end(i)
                            spack_2d(ijk, inob)%dlmat(jj, jj) = igamma_dtstep - spack(inob)%dldrdc(ijk, jj, jj)
                        enddo
                    enddo

                    if ((i .eq. 1)) then
                        blocknonzeros = 0
                        do ji = 1, nspecies
                            do jj = 1, nspecies
                                !IF(spack_2d(1,inob)%DLmat(Ji,Jj)<=0 .AND. spack_2d(1,inob)%DLmat(Ji,Jj)>=(-1)*0) CYCLE
                                if ((spack_2d(1, inob)%dlmat(ji, jj) .ne. 0.0d+00)) then
                                    blocknonzeros = blocknonzeros + 1
                                    ipos(blocknonzeros) = ji
                                    jpos(blocknonzeros) = jj
                                endif
                            enddo
                        enddo
                        numberofnonzeros = blocknonzeros

                        !create matrix
                        matrix_id = sfcreate_solve(sizeofmatrix, complex_t, error)

                        if (allocated(element)) deallocate (element)
                        allocate(element(numberofnonzeros))

                        do nz = 1, numberofnonzeros
                            ji = ipos(nz)
                            jj = jpos(nz)
                            element(nz) = sfgetelement(matrix_id, ji, jj)
                            !element(nz)=sfGetElement(matrix_Id,INT(DLmat(1,nz,1)),INT(DLmat(1,nz,2)))
                        end do

                        call sfzero(matrix_id)
                        do nz = 1, numberofnonzeros
                            ji = ipos(nz)
                            jj = jpos(nz)
                            call sfadd1real(element(nz), spack_2d(1, inob)%dlmat(ji, jj))
                            !CALL sfAdd1Real(element(nz),DLmat(1,nz,3))
                        end do

                        error = sffactor(matrix_id)
                    endif

                    !@LNCC: begin points loop
                    do ijk = 1, block_end(i) !index_g%block_end(i)
                        call sfzero(matrix_id)
                        do nz = 1, numberofnonzeros
                            ji = ipos(nz)
                            jj = jpos(nz)
                            call sfadd1real(element(nz), spack_2d(ijk, inob)%dlmat(ji, jj))
                            !CALL sfAdd1Real(element(nz),DLmat(ijk,nz,3))
                        end do
                        error = sffactor(matrix_id)

                        !---------------------------------------------------------------------------------------------------------------------
                        !    1- First step
                        !-   Compute DLk1 by Solving (1/(Igamma*dt) - DLRDC) DLk1=DLR
                        !- solver sparse
                        allocate(rhs_tmp(nspecies), sol_tmp(nspecies))
                        rhs_tmp = spack_2d(ijk, inob)%dlb1
                        call sfsolve_c(matrix_id, c_loc(rhs_tmp), c_loc(sol_tmp))
                        spack_2d(ijk, inob)%dlk1 = sol_tmp
                        deallocate(rhs_tmp, sol_tmp)
                        !
                        !---------------------------------------------------------------------------------------------------------------------
                        !    2- Second step
                        !    compute   K2 by solving (1/0.5 h JAC)K2 = 4/h * K1 +  F(Yn)
                        !    compute DLK2 by solving (1/Igama*dt - DLRDC)DLK2 = 4/h * DLK1 +  DLb1
                        dt_chem_i = 1.0d0 / dt_chem

                        do ji = 1, nspecies
                            spack_2d(ijk, inob)%dlb2(ji) = (4.0d0 * dt_chem_i) * spack_2d(ijk, inob)%dlk1(ji) + &
                            spack_2d(ijk, inob)%dlb1(ji)
                        enddo
                        allocate(rhs_tmp(nspecies), sol_tmp(nspecies))
                        rhs_tmp = spack_2d(ijk, inob)%dlb2
                        call sfsolve_c(matrix_id, c_loc(rhs_tmp), c_loc(sol_tmp))
                        spack_2d(ijk, inob)%dlk2 = sol_tmp
                        deallocate(rhs_tmp, sol_tmp)
                        !do n=1,nspecies
                        !   PRINT *,'LFR-DBG: DLk2',ijk,n,spack_2d(ijk,inob)%dlk2(n)
                        !enddo
                        !---------------------------------------------------------------------------------------------------------------------
                        !    3- Third step
                        !    a) update concentrations

                        !dt_chem_i = 1.0d0/dt_chem
                        do ji = 1, nspecies
                            spack(inob)%sc_p_new(ijk, ji) = spack(inob)%sc_p(ijk, ji) + 2.0d0 * spack_2d(ijk, inob)%dlk1(ji)

                            if (spack(inob)%sc_p_new(ijk, ji) .lt. threshold1) then
                                spack(inob)%sc_p_new(ijk, ji) = threshold1
                                spack_2d(ijk, inob)%dlk1(ji) = 0.5d0 * (spack(inob)%sc_p_new(ijk, ji) - spack(inob)%sc_p(ijk, ji))
                            endif
                        enddo
                        !
                        !    b) update the net production term (DLr= P-L = F(Y3)) at this stage with the first-order
                        !       approximation with the new concentration
                        !
                        call fexchem (spack(inob)%sc_p_new &
                        , spack(inob)%rk &
                        , spack(inob)%dlr &
                        , nspecies, ijk &
                        , ijk & !index_g%block_end(i), &
                        , maxblock_size, nr)

                        !    c) compute   K3 by solving (1 /(0.5 h) - JAC)K3 =   F(Y3) + 0.5 (K1-K2)
                        do ji = 1, nspecies
                            spack_2d(ijk, inob)%dlb3(ji) = spack(inob)%dlr(ijk, ji) + dt_chem_i * &
                            (spack_2d(ijk, inob)%dlk1(ji) - spack_2d(ijk, inob)%dlk2(ji))
                        enddo
                        allocate(rhs_tmp(nspecies), sol_tmp(nspecies))
                        rhs_tmp = spack_2d(ijk, inob)%dlb3
                        call sfsolve_c(matrix_id, c_loc(rhs_tmp), c_loc(sol_tmp))
                        spack_2d(ijk, inob)%dlk3 = sol_tmp
                        deallocate(rhs_tmp, sol_tmp)
                        !do n=1,nspecies
                        !   PRINT *,'LFR-DBG: DLk3',ijk,n,spack_2d(ijk,inob)%dlk3(n)
                        !enddo
                        !---------------------------------------------------------------------------------------------------------------------
                        !    4- Fourth step
                        !    a) update concentrations
                        !       Y4 = Yn + 2 * k1 +  K3
                        dt_chem_i = 1.0d0 / dt_chem
                        do ji = 1, nspecies
                            spack(inob)%sc_p_new(ijk, ji) = spack(inob)%sc_p(ijk, ji) + &
                            2.0d0 * spack_2d(ijk, inob)%dlk1(ji) + &
                            spack_2d(ijk, inob)%dlk3(ji)

                            if (spack(inob)%sc_p_new(ijk, ji) .lt. threshold1) then
                                spack(inob)%sc_p_new(ijk, ji) = threshold1
                                spack_2d(ijk, inob)%dlk3(ji) = (spack(inob)%sc_p_new(ijk, ji) - spack(inob)%sc_p(ijk, ji)) &
                                -2.0d0 * spack_2d(ijk, inob)%dlk1(ji)
                            endif

                        enddo
                        !    b) update the net production term (DLr= P-L = F(Y4) ) at this stage with the 3rd-order
                        !       approximation with the new concentration
                        !
                        call fexchem (spack(inob)%sc_p_new & ! Y4
                        , spack(inob)%rk &
                        , spack(inob)%dlr & ! F(Y4)
                        , nspecies, ijk &
                        , ijk & !index_g%block_end(i)
                        , maxblock_size, nr)

                        !    c) compute   K4 by solving (1/(0.5 h)- JAC)K4 =    F(Y3) + K1/h -K2/h -8/3 K3/h
                        do ji = 1, nspecies
                            spack_2d(ijk, inob)%dlb4(ji) = spack(inob)%dlr (ijk, ji) + dt_chem_i * & ! F(Y4)
                            (spack_2d(ijk, inob)%dlk1(ji) &
                            - spack_2d(ijk, inob)%dlk2(ji) &
                            - c83 * spack_2d(ijk, inob)%dlk3(ji))
                        enddo
                        allocate(rhs_tmp(nspecies), sol_tmp(nspecies))
                        rhs_tmp = spack_2d(ijk, inob)%dlb4
                        call sfsolve_c(matrix_id, c_loc(rhs_tmp), c_loc(sol_tmp))
                        spack_2d(ijk, inob)%dlk4 = sol_tmp
                        deallocate(rhs_tmp, sol_tmp)

                    enddo
                    !@LNCC: end points loop


                    !---------------------------------------------------------------------------------------------------------------------
                    !   - the solution
                    dt_chem_i = 1.0d0 / dt_chem
                    do ji = 1, nspecies
                        do ijk = 1, block_end(i) !index_g%block_end(i)

                            spack(inob)%sc_p_new(ijk, ji) = spack(inob)%sc_p(ijk, ji) + &
                            2.0d0 * spack_2d(ijk, inob)%dlk1(ji) &
                            + spack_2d(ijk, inob)%dlk3(ji) &
                            + spack_2d(ijk, inob)%dlk4(ji)
                            spack(inob)%sc_p_new(ijk, ji) = max (spack(inob)%sc_p_new(ijk, ji), threshold1)

                        enddo
                    enddo

                    !
                    !-  Compute the error estimation : spack(inob)%err
                    do ijk = 1, block_end(i) !index_g%block_end(i)

                        spack(inob)%err(ijk) = 0.0d0
                        do ji = 1, nspecies
                            if (spack(inob)%sc_p_new(ijk, ji) .gt. 1.0e+14) then
                                print *, 'scp e scp_new = ', ji, ijk, spack(inob)%sc_p(ijk, ji), spack(inob)%sc_p_new(ijk, ji)
                            endif
                            tol = atol(ji) + rtol(ji) * dmax1(dabs(spack(inob)%sc_p(ijk, ji)), dabs(spack(inob)%sc_p_new(ijk, ji)))

                            err1 = spack_2d(ijk, inob)%dlk4(ji)

                            spack(inob)%err(ijk) = spack(inob)%err(ijk) + (err1 / tol)**2.0d0

                        enddo

                        spack(inob)%err(ijk) = dmax1(uround, dsqrt(spack(inob)%err(ijk) / nspecies))
                    enddo

                    all_accepted = 1
                    do ijk = 1, block_end(i) !index_g%block_end(i)
                        if (spack(inob)%err(ijk) - roundoff > 1.0d0) then
                            all_accepted = 0 ;
                            exit
                        endif
                    enddo

                    !- find the maximum error occurred
                    max_err1 = maxval(spack(inob)%err(1:block_end(i))) !index_g%block_end(i) ) )

                    !- use it to determine the new time step for all block elements
                    !- new step size is bounded by FacMin <= Hnew/H <= FacMax
                    fac = min(facmax, max(facmin, facsafe / max_err1**(1.0d0 / elo)))

                    !- possible new timestep
                    dt_new = dt_chem * fac
                    dt_new = max(dt_min, min(dt_max, dt_new))

                    !- to reset the timestep resizing in function of the error estimation, use the statements below:
                    ! all_accepted = 1; dt_new=dt_max

                    if (all_accepted == 0 .and. dt_new > dt_min) then  ! current solution is not accepted

                        !- resize the timestep and try again
                        dt_chem = dt_new

                    else    ! current solution is     accepted
                        !count_blocks_accept = count_blocks_accept+1
                        !- go ahead, updating spack(inob)%sc_p with the solution (spack(inob)%sc_p_new)
                        !- next time
                        time_c = time_c + dt_chem

                        !- next timestep (dt_new but limited by the time_f-time_c, the rest of time integration interval)
                        dt_chem = min(dt_new, time_f - time_c)
                        !- save the accepted timestep for the next integration interval
                        if (time_c < time_f) last_accepted_dt(i) = dt_new !index_g%last_accepted_dt(i) =    dt_new

                        !- pointer (does not work yet)
                        ! spack(inob)%sc_p=>spack(inob)%sc_p_new     ! POINTER
                        !- copy
                        do ji = 1, nspecies
                            do ijk = 1, block_end(i) !index_g%block_end(i)
                                spack(inob)%sc_p(ijk, ji) = spack(inob)%sc_p_new(ijk, ji)
                                ! spack(inob)%sc_p_4(ijk,Ji) = spack(inob)%sc_p_new(ijk,Ji)
                            enddo
                        enddo

                        exit untilaccepted

                    endif

                end do untilaccepted


            enddo run_until_integr_ends ! time-spliting


            !--------------------------------------------------------------------------------------------
            !- Restoring species tendencies OR updated mixing ratios from internal to brams structure

            !- transported species section

            if (split_method == 'PARALLEL' .and. n_dyn_chem > 1) then

                do ijk = 1, block_end(i) !index_g%block_end(i)

                    do ispc = 1, nspecies_chem_transported

                        !- map the species to transported ones
                        n = transp_chem_index(ispc)

                        !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
                        chem_g(n)%sc_p(ijk, i) = chem_g(n)%sc_t_dyn(ijk,i) * n_dyn_chem * dtlt + &
                        spack(inob)%sc_p(ijk, n) * weight(n) * spack(inob)%volmol_i(ijk)
                        chem_g(n)%sc_p(ijk, i) = max(0., chem_g(n)%sc_p(ijk, i))

                    end do
                enddo

            elseif(split_method == 'PARALLEL' .and. n_dyn_chem == 1) then

                dble_dtlt_i = 1.0d0 / dble(dtlt)
                do ijk = 1, block_end(i) !index_g%block_end(i)

                    do ispc = 1, nspecies_chem_transported

                        !- map the species to transported ones
                        n = transp_chem_index(ispc)

                        !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
                        !           chem_g(n)%sc_t(kij_) =          +  &! use this for update only chemistry (No dyn/emissions
                        chem_g(n)%sc_t(ijk,i) = chem_g(n)%sc_t(ijk,i) + &! previous tendency
                        (spack(inob)%sc_p(ijk, n) * weight(n) * spack(inob)%volmol_i(ijk) - &! new mixing ratio
                        chem_g(n)%sc_p(ijk, i)) &! old mixing ratio
                         * dble_dtlt_i                     ! inverse of timestep
                    end do
                end do

            else

                do ijk = 1, block_end(i) !index_g%block_end(i)

                    do ispc = 1, nspecies_chem_transported

                        !- map the species to transported ones
                        n = transp_chem_index(ispc)

                        chem_g(n)%sc_p(ijk, i) = spack(inob)%sc_p(ijk, n) * weight(n) * spack(inob)%volmol_i(ijk)
                        chem_g(n)%sc_p(ijk, i) = max(0., chem_g(n)%sc_p(ijk, i))
                    end do
                end do

            endif


            !- no transported species section
            do ijk = 1, block_end(i) !index_g%block_end(i)

                do ispc = 1, nspecies_chem_no_transported

                    !- map the species to no transported ones
                    n = no_transp_chem_index(ispc)

                    !- save no-transported species (keep current unit : molec/cm3)
                    !LFR-MONAN chem_g(n)%sc_p(k_, i_, j_) = max(0., real (spack(inob)%sc_p(ijk, n)))
                    chem_g(n)%sc_p(ijk, i) = max(0., real (spack(inob)%sc_p(ijk, n)))
                end do

            end do

        end do ! enddo loop over all blocks


    end subroutine chem_rodas3_dyndt


!--------------------------------------------------------------------------  
  subroutine get_number_nonzeros(jppj,nr,nspecies,nonzeros)

    integer          , intent(in)    :: jppj
    integer          , intent(in)    :: nr
    integer          , intent(in)    :: nspecies
    integer          , intent(out) :: nonzeros

    double precision ,dimension(nspecies,nspecies) :: def_non_zeros 
    double precision ,dimension(nspecies) :: sc_p
    double precision  :: xlw,vapp(1),cosz(1),temp(1),press(1),att(1)
    double precision  :: rk(nr)
    double precision  :: jphoto(jppj)
    integer :: i,ji,jj

    jphoto(:) = 2.3333331d0
    xlw	  = 1.d0
    vapp(:) = 1.d15
    cosz(:) = 1.d0
    att(:)  = 1.0d0
    temp(:) = 273.15d0
    press(:)= 1.0d5
    sc_p    = 1.d15 ! dummy concentration to get the maximum number
                    ! of possible non zero elements

    call kinetic(jppj,jphoto	 &
     		     ,rk(1:nr)   &
     		     ,temp	 &
     		     ,vapp	 &
     		     ,press	 &
     		     ,1,1,1,nr  )

    def_non_zeros = 0.d0

    call jacdchemdc(sc_p,rk(1:nr),def_non_zeros, & ! jacobian matrix
                    nspecies,1,1,1,nr)

    nonzeros = 0
    do jj=1,nspecies
       def_non_zeros(jj,jj)=1.d0+ def_non_zeros(jj,jj)
       do ji=1,nspecies
          if (def_non_zeros(ji,jj) .ne. 0.d0) then
             nonzeros = nonzeros + 1
          endif
       enddo
    enddo

    print *,'LFR-DBG: get_number_nonzeros: ',nonzeros

  end subroutine get_number_nonzeros


end module mod_chem_spack_rodas3_dyndt
