
module mod_chem_spack_rodas3_dyndt
    use mpas_pool_routines
    use mpas_derived_types
    use modMemoryChem, only: &
    chem_vars        ! Type

    use mem_spack, only: &
    spack_type, & ! Type
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
    ,maxnspecies &
    , spc_name

    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_ptr, c_loc

    implicit none

    !private


    !public chem_rodas3_dyndt

    !- interface do solver de matriz esparsa (Sangiovanni-Vincentelli / netlib sparse),
    !- movida para o nível do módulo para ser compartilhada por chem_rodas3_dyndt e rosenbrok
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
            import :: c_double
            integer(kind=8), intent(in) :: elem
            real(kind=c_double), intent(in) :: value
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


contains

    !========================================================================================
    !>  Ponto de entrada da integração química do RELACS via RODAS3. Para cada
    !!  bloco de pontos de grade `i = 1..nob`. No MONAN nob=nVertLevels
    !!
    !! @author Rodrigues, L.F. 
    !! @date 2026-08-26
    !! @version 0.1.0
    !!
    !!  1. converte as concentrações de ppbm (BRAMS) para molec/cm³ (`spack%sc_p`);
    !!  2. calcula as constantes de taxa das 128 reações do RELACS (`kinetic`),
    !!     incluindo as 17 fotolíticas (`jphoto`);
    !!  3. chama `rosenbrok`, que integra a química daquele bloco do tempo
    !!     \( 0 \) até \( dtlt\times n\_dyn\_chem \);
    !!  4. restaura `sc_p`/`sc_t` na estrutura do BRAMS, conforme `split_method`.
    !!
    !!  A matemática do solver (estágios do RODAS3, controle de passo, fatoração
    !!  esparsa) está documentada em `rosenbrok`.
    !! 
    !! @license This code is under GPLv3 License, see it in <https://www.gnu.org/licenses/gpl-3.0.en.html>
    !==============================================================================
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

        real, parameter :: cp = 1004.
        real, parameter :: rgas = 287.
        real, parameter :: cpor = cp / rgas
        real, parameter :: p00 = 1.0e5
        real(kind=c_double), parameter :: pmar = 28.96_c_double
        real(kind=c_double), parameter :: threshold1 = 0.0_c_double
        real(kind=c_double), parameter :: threshold2 = -1.e-1_c_double
        real(kind=c_double), parameter :: igamma = 0.5_c_double
        real(kind=c_double), parameter :: c43 = 4._c_double / 3._c_double
        real(kind=c_double), parameter :: c83 = 8._c_double / 3._c_double
        real(kind=c_double), parameter :: c56 = 5._c_double / 6._c_double
        real(kind=c_double), parameter :: c16 = 1._c_double / 6._c_double
        real(kind=c_double), parameter :: c112 = 1._c_double / 12._c_double
        !- parameters for dynamic timestep control
        real(kind=c_double), parameter :: facmin = 0.2_c_double ! lower bound on step decrease factor (default=0.2)
        real(kind=c_double), parameter :: facmax = 6.0_c_double ! upper bound on step increase factor (default=6)
        real(kind=c_double), parameter :: facrej = 0.1_c_double ! step decrease factor after multiple rejections
        real(kind=c_double), parameter :: facsafe = 0.9_c_double ! step by which the new step is slightly smaller
        ! than the predicted value  (default=0.9)
        real(kind=c_double), parameter :: uround = 1.e-15_c_double, elo = 3.0_c_double, roundoff = 1.e-8_c_double
        integer, parameter :: complex_t = 0
        real(kind=c_double), parameter ::  rtols=1.e-3_c_double ! 1e-2 means two digits
        real(kind=c_double), parameter ::  atols=1.e+7_c_double ! Jacobson (1998, SMVGEAR) range 1.e3-1.e7! 1.D0


        real(kind=c_double), allocatable :: dlmat(:, :, :)
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
        real(kind=c_double) :: dt_min, dt_max, dt_actual, dt_new
        real(kind=c_double) :: time_f, time_c
        real(kind=c_double) :: fac, tol, err1, max_err1
        real(kind=c_double) :: fxc, igamma_dtstep, dt_chem, dt_chem_i
        integer :: i, ijk, n, j, k, ispc, ji, jj, k_, i_, j_, kij_, kij
        real(kind=c_double) :: atol(nspecies)
        real(kind=c_double) :: rtol(nspecies)

        !integer :: maxblock_size
        real(kind=c_double), allocatable, target :: rhs_tmp(:), sol_tmp(:)

        call get_number_nonzeros(nr_photo, nr, nspecies, spack%rk(1, 1:nr), &
        spack%jphoto(1, 1:nr_photo), p00, maxnonzeros)

        sizeofmatrix = nspecies
        !maxblock_size = nvertlevels
        atol = atols
        rtol = rtols

        allocate(ipos(maxnonzeros)) ;
        ipos = 1
        allocate(jpos(maxnonzeros)) ;
        jpos = 1


        do i = 1, nob !- loop over all blocks i) = nVertLevels

            !> Converte as concentrações do bloco `i` de ppbm para molec/cm³, preenche
            !! `spack%jphoto` e calcula as constantes de taxa `spack%rk` (via `kinetic`).
            call prepare_spack_block( &
                i                            &
              , block_end                    &
              , dtlt                         &
              , press                        &
              , temp                         &
              , vapp                         &
              , n_dyn_chem                   &
              , split_method                 &
              , jphoto                       &
              , chem_g                      &
              , nspecies_chem_transported    &
              , nspecies_chem_no_transported &
              , transp_chem_index            &
              , no_transp_chem_index         &
              , maxblock_size                &
            )

            dt_chem = last_accepted_dt(i) !index_g%last_accepted_dt(i) ! dble(dtlt)
            dt_min = max(1.0_c_double, 1.e-2_c_double * real(dtlt * n_dyn_chem, kind=c_double))
            dt_max = real(dtlt * n_dyn_chem, kind=c_double)
            dt_new = 0.0_c_double
            time_c = 0.0_c_double
            time_f = real(dtlt * n_dyn_chem, kind=c_double)

            !- ROSENBROCK METHOD ----------------------------------------------------------------------
            !- integra a química deste bloco (i) do tempo 0 até time_f, com sub-passos adaptativos,
            !- dentro da subroutine rosenbrok (antigo laço run_until_integr_ends)
            call rosenbrok( &
                i                  &
              , block_end          &
              , maxblock_size      &
              , sizeofmatrix       &
              , dt_min             &
              , dt_max             &
              , time_f             &
              , time_c             &
              , dt_chem            &
              , last_accepted_dt   &
              , atol               &
              , rtol               &
              , matrix_id          &
              , numberofnonzeros   &
              , ipos               &
              , jpos               &
              , element            &
            )


            !> Restaura, na estrutura do modelo (`chem_g`), a concentração ou a
            !! tendência resultante da integração química deste bloco, conforme
            !! `split_method`/`n_dyn_chem`.
            call restore_chem1g_block( &
                i                            &
              , block_end                    &
              , dtlt                         &
              , n_dyn_chem                   &
              , split_method                 &
              , chem_g                      &
              , nspecies_chem_transported    &
              , nspecies_chem_no_transported &
              , transp_chem_index            &
              , no_transp_chem_index         &
            )

        end do ! enddo loop over all blocks


    end subroutine chem_rodas3_dyndt

    !========================================================================================
    !>  Prepara, em `spack`, o estado de entrada da química para **um único bloco** `i`
    !!  de pontos de grade, antes da integração propriamente dita (feita por
    !!  `rosenbrok`). Três tarefas:
    !!
    !!  1. converte a concentração de cada espécie de ppbm (mixing ratio de massa,
    !!     como armazenado em `chem_g`) para molec/cm³ (`spack%sc_p`), via
    !!     \( c = c_{ppbm}\cdot volmol / weight \), com
    !!     \( volmol = \dfrac{N_A\cdot 10^{-15}\cdot pmar\cdot press}{8{,}314\cdot temp} \);
    !!  2. copia as taxas de fotólise pré-calculadas (`jphoto`) para `spack%jphoto`;
    !!  3. calcula as constantes de taxa das 128 reações do RELACS (`spack%rk`),
    !!     chamando `kinetic` com `press`/`temp`/`vapp` convertidos para double
    !!     precision só no momento da chamada (`dble(...)`), sem copiá-los antes.
    subroutine prepare_spack_block( &
        i                            &
      , block_end                    &
      , dtlt                         &
      , press                        &
      , temp                         &
      , vapp                         &
      , n_dyn_chem                   &
      , split_method                 &
      , jphoto                       &
      , chem_g                       &
      , nspecies_chem_transported    &
      , nspecies_chem_no_transported &
      , transp_chem_index            &
      , no_transp_chem_index         &
      , maxblock_size                &
    )

        implicit none

        !> índice do bloco sendo preparado nesta chamada
        integer,             intent(in) :: i
        !> número de pontos de grade válidos em cada bloco
        integer,             intent(in) :: block_end(:)
        !> passo de tempo do modelo (s)
        real,                intent(in) :: dtlt
        !> pressão, temperatura e razão de mistura de vapor d'água (dimensões originais, por bloco)
        real,                intent(in) :: press(:,:), temp(:,:), vapp(:,:)
        !> número de sub-passos de química dentro de um passo de tempo do modelo
        integer,             intent(in) :: n_dyn_chem
        !> esquema de *splitting* de tendências ('PARALLEL' ou sequencial)
        character(len = 20), intent(in) :: split_method
        !> taxas de fotólise pré-calculadas (ex.: pelo Fast-TUV), por nível/reação/ponto
        real,                intent(in) :: jphoto(:,:,:)
        !> estrutura de espécies do BRAMS; só lida aqui (`sc_p`, `sc_t_dyn`)
        type(chem_vars),    intent(in) :: chem_g(:)
        !> número de espécies transportadas / não-transportadas e seus índices de mapeamento
        integer,             intent(in) :: nspecies_chem_transported
        integer,             intent(in) :: nspecies_chem_no_transported
        integer,             intent(in) :: transp_chem_index(:)
        integer,             intent(in) :: no_transp_chem_index(:)
        !> tamanho máximo de bloco
        integer,             intent(in) :: maxblock_size

        !> massa molar do ar seco (g/mol), usada no cálculo de `volmol`
        real(kind=c_double), parameter :: pmar = 28.96_c_double

        integer :: ijk,ispc, n

        !- copying structure from input to internal
        do ijk = 1, block_end(i) !index_g%block_end(i) - MONAN: the block_end is the number of levels
            spack%volmol(ijk) = (6.02e23_c_double * 1e-15_c_double * pmar) * (press(i,ijk)) / (8.314_c_double * temp(i,ijk))
            spack%volmol_i(ijk) = 1.0_c_double / spack%volmol(ijk)

            !- no transported species section
            do ispc = 1, nspecies_chem_no_transported
                !- map the species to NO transported ones
                n = no_transp_chem_index(ispc)
                !- initialize no-transported species (don't need to convert, because these
                !- species are already saved using molecule/cm^3 units)
                spack%sc_p(ijk, n) = chem_g(n)%sc_p(i,ijk)
            end do

        end do
        !- convert from brams chem (ppbm) arrays to spack (molec/cm3)
        !- transported species section
        if (split_method == 'PARALLEL' .and. n_dyn_chem > 1) then
            do ijk = 1, block_end(i) !index_g%block_end(i)

                do ispc = 1, nspecies_chem_transported

                    !- map the species to transported ones
                    n = transp_chem_index(ispc)

                    spack%sc_p(ijk, n) = (chem_g(n)%sc_p(i,ijk) - &  ! updated mixing ratio
                    chem_g(n)%sc_t_dyn(i,ijk) * n_dyn_chem * dtlt)&  ! accumulated tendency
                    * spack%volmol(ijk) / weight(n)

                end do
            end do
        else
            do ijk = 1, block_end(i) !index_g%block_end(i)

                do ispc = 1, nspecies_chem_transported

                    !- map the species to transported ones
                    n = transp_chem_index(ispc)

                    !- conversion from ppbm to molecule/cm^3
                    spack%sc_p(ijk, n) = chem_g(n)%sc_p(i,ijk) * spack%volmol(ijk) / weight(n)
                    spack%sc_p(ijk, n) = max(0._c_double, spack%sc_p(ijk, n))
                end do

            end do
        endif

        !- Photolysis section
        if (trim(photojmethod) == 'FAST-JX' .or. trim(photojmethod) == 'FAST-TUV') then

            do ijk = 1, block_end(i) !index_g%block_end(i)
                do n = 1, nr_photo 
                    spack%jphoto(ijk, n) = jphoto(i,ijk,n)
                end do
            enddo

        elseif(trim(photojmethod) == 'LUT') then

            do ijk = 1, block_end(i) !index_g%block_end(i)

                !- UV  attenuation un function AOT
                !spack%att(ijk) = att(i)

                !- get zenital angle (for LUT PhotojMethod)
                !spack%cosz(ijk) = cosz(i)

            enddo
        endif

        !- compute kinetical and photochemical reactions (array: spack%rk)
        !- press/temp/vapp usados diretamente nas dimensões originais (i,ijk), sem
        !- cópia prévia para spack%press/%temp/%vapp; dble() converte de real
        !- (precisão simples, como chegam em chem_rodas3_dyndt) para real(kind=c_double)
        !- (o que kinetic() espera), só no momento da chamada
        call kinetic(nr_photo, spack%jphoto &
        , spack%rk &
        , real(temp(i,:), kind=c_double) &
        , real(vapp(i, :), kind=c_double) &
        , real(press(i,:), kind=c_double) &
        , 1 &
        , block_end(i) & !index_g%block_end(i), &
        , maxblock_size, nr)

    end subroutine prepare_spack_block

    !========================================================================================
    !>  Devolve, para a estrutura do modelo (`chem_g`), o resultado da integração
    !!  química deste bloco `i` — a contraparte de saída de `prepare_spack_block`.
    !!  O que é escrito depende de `split_method`/`n_dyn_chem`:
    !!
    !!  - `'PARALLEL'` com \( n\_dyn\_chem>1 \): soma a tendência dinâmica já
    !!    acumulada (`chem_g(n)%sc_t_dyn`) à concentração pós-química, convertida
    !!    de volta para ppbm, e grava direto em `chem_g(n)%sc_p`;
    !!  - `'PARALLEL'` com \( n\_dyn\_chem=1 \): não sobrescreve `sc_p`, só
    !!    **acrescenta** a `chem_g(n)%sc_t` a tendência devida à química,
    !!    \( \big(c^{\text{novo}}_{ppbm}-c^{\text{antigo}}_{ppbm}\big)/dtlt \);
    !!  - caso contrário (*splitting* sequencial): sobrescreve `chem_g(n)%sc_p`
    !!    direto com a concentração pós-química, convertida para ppbm.
    !!
    !!  Espécies não-transportadas sempre sobrescrevem `sc_p` direto, sem
    !!  tendência, mantendo a unidade interna (molec/cm³).
    subroutine restore_chem1g_block( &
        i                            &
      , block_end                    &
      , dtlt                         &
      , n_dyn_chem                   &
      , split_method                 &
      , chem_g                      &
      , nspecies_chem_transported    &
      , nspecies_chem_no_transported &
      , transp_chem_index            &
      , no_transp_chem_index         &
    )

        implicit none

        !> índice do bloco cujo resultado está sendo devolvido nesta chamada
        integer,             intent(in)    :: i
        !> número de pontos de grade válidos em cada bloco
        integer,             intent(in)    :: block_end(:)
        !> passo de tempo do modelo (s)
        real,                intent(in)    :: dtlt
        !> número de sub-passos de química dentro de um passo de tempo do modelo
        integer,             intent(in)    :: n_dyn_chem
        !> esquema de *splitting* de tendências ('PARALLEL' ou sequencial)
        character(len = 20), intent(in)    :: split_method
        !> estrutura de espécies do BRAMS; `sc_p` e/ou `sc_t` são atualizados aqui
        type(chem_vars),    intent(inout) :: chem_g(:)
        !> número de espécies transportadas / não-transportadas e seus índices de mapeamento
        integer,             intent(in)    :: nspecies_chem_transported
        integer,             intent(in)    :: nspecies_chem_no_transported
        integer,             intent(in)    :: transp_chem_index(:)
        integer,             intent(in)    :: no_transp_chem_index(:)

        integer :: ijk,n,ispc
        real(kind=c_double) :: dble_dtlt_i

        !--------------------------------------------------------------------------------------------
        !- Restoring species tendencies OR updated mixing ratios from internal to brams structure

        !- transported species section

        if (split_method == 'PARALLEL' .and. n_dyn_chem > 1) then

            do ijk = 1, block_end(i) !index_g%block_end(i)

                do ispc = 1, nspecies_chem_transported

                    !- map the species to transported ones
                    n = transp_chem_index(ispc)

                    !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
                    chem_g(n)%sc_p(i,ijk) = chem_g(n)%sc_t_dyn(i,ijk) * n_dyn_chem * dtlt + &
                    spack%sc_p(ijk, n) * weight(n) * spack%volmol_i(ijk)
                    chem_g(n)%sc_p(i,ijk) = max(0., chem_g(n)%sc_p(i,ijk))

                end do
            enddo

        elseif(split_method == 'PARALLEL' .and. n_dyn_chem == 1) then

            dble_dtlt_i = 1.0_c_double / real(dtlt, kind=c_double)
            do ijk = 1, block_end(i) !index_g%block_end(i)

                do ispc = 1, nspecies_chem_transported

                    !- map the species to transported ones
                    n = transp_chem_index(ispc)

                    !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
                    !           chem_g(n)%sc_t(kij_) =          +  &! use this for update only chemistry (No dyn/emissions
                    chem_g(n)%sc_t(i,ijk) = chem_g(n)%sc_t(i,ijk) + &! previous tendency
                    (spack%sc_p(ijk, n) * weight(n) * spack%volmol_i(ijk) - &! new mixing ratio
                    chem_g(n)%sc_p(i,ijk)) &! old mixing ratio
                     * dble_dtlt_i                     ! inverse of timestep
                end do
            end do

        else

            do ijk = 1, block_end(i) !index_g%block_end(i)

                do ispc = 1, nspecies_chem_transported

                    !- map the species to transported ones
                    n = transp_chem_index(ispc)

                    chem_g(n)%sc_p(i,ijk) = spack%sc_p(ijk, n) * weight(n) * spack%volmol_i(ijk)
                    chem_g(n)%sc_p(i,ijk) = max(0., chem_g(n)%sc_p(i,ijk))
                end do
            end do

        endif


        !- no transported species section
        do ijk = 1, block_end(i) !index_g%block_end(i)

            do ispc = 1, nspecies_chem_no_transported

                !- map the species to no transported ones
                n = no_transp_chem_index(ispc)

                !- save no-transported species (keep current unit : molec/cm3)
                !LFR-MONAN chem_g(n)%sc_p(k_, i_, j_) = max(0., real (spack%sc_p(ijk, n)))
                chem_g(n)%sc_p(i,ijk) = max(0., real (spack%sc_p(ijk, n)))
            end do

        end do

    end subroutine restore_chem1g_block

    !========================================================================================
    !>  Integra a química de **um único bloco** de pontos de grade, do tempo
    !!  \( t_c = 0 \) até \( t_c = t_f \) (`time_f`, tipicamente
    !!  \( dtlt \times n\_dyn\_chem \)), usando o método de Rosenbrock **RODAS3**:
    !!  um Runge-Kutta linearmente implícito, de 4 estágios e 3ª ordem, L-estável,
    !!  adequado à rigidez (*stiffness*) do sistema de EDOs da química atmosférica.
    !!
    !!  Extraída do antigo laço `run_until_integr_ends` de [[chem_rodas3_dyndt]].
    !!
    !!  ### O sistema integrado
    !!
    !!  Para cada espécie \( j \), a EDO química tem a forma
    !!  \[
    !!  \frac{dc_j}{dt} = F_j(\mathbf{c}) = P_j(\mathbf{c}) - L_j(\mathbf{c})\,c_j
    !!  \]
    !!  onde \( P_j \) é a produção e \( L_j\,c_j \) a perda química, calculadas por
    !!  `fexchem` a partir das taxas de reação `rk` (já calculadas por `kinetic`,
    !!  fora desta rotina, em [[chem_rodas3_dyndt]]).
    !!
    !!  ### Os 4 estágios do RODAS3
    !!
    !!  A cada tentativa de passo \( h \) (`dt_chem`), com \( \gamma=1/2 \)
    !!  (`igamma`), resolve-se sempre o **mesmo** sistema linear \( A\,k=b \), com
    !!  \[
    !!  A = \frac{1}{\gamma h}\,\mathbf{I} - \mathbf{J}, \qquad
    !!  \mathbf{J} = \left.\frac{\partial F}{\partial\mathbf{c}}\right|_{\mathbf{c}=\mathbf{y}_0}
    !!  \]
    !!  (a Jacobiana `spack%dldrdc`, calculada uma única vez por sub-passo, fora do
    !!  laço de tentativas), variando apenas o lado direito \( b \) a cada estágio:
    !!
    !!  | Estágio | Lado direito \( b \)                                                    | Nova avaliação de \(F\)? |
    !!  |---------|--------------------------------------------------------------------------|---------------------------|
    !!  | \(k_1\) | \( F(\mathbf{y}_0) \)                                                     | não (reaproveitado)      |
    !!  | \(k_2\) | \( \dfrac{4}{h}\,k_1 + F(\mathbf{y}_0) \)                                 | não                       |
    !!  | \(k_3\) | \( F(\mathbf{y}_3) + \dfrac{1}{h}(k_1-k_2) \)                             | sim, em \(\mathbf{y}_3=\mathbf{y}_0+2k_1\) |
    !!  | \(k_4\) | \( F(\mathbf{y}_4) + \dfrac{1}{h}\Big(k_1-k_2-\dfrac{8}{3}k_3\Big) \)     | sim, em \(\mathbf{y}_4=\mathbf{y}_0+2k_1+k_3\) |
    !!
    !!  A solução de 3ª ordem do sub-passo é
    !!  \[
    !!  \mathbf{y}_{n+1} = \mathbf{y}_0 + 2k_1 + k_3 + k_4
    !!  \]
    !!
    !!  ### Controle adaptativo do passo \( h \)
    !!
    !!  O erro local é estimado direto do último estágio \( k_4 \), sem custo extra
    !!  (a diferença entre a solução de 3ª ordem e uma solução embutida de ordem
    !!  inferior colapsa nele):
    !!  \[
    !!  \text{err} = \max\!\left(\text{uround},\ \sqrt{\frac{1}{n_{esp}}
    !!  \sum_{j=1}^{n_{esp}}\left(\frac{k_{4,j}}
    !!  {\,atol_j+rtol_j\cdot\max(|c_j^{\text{antigo}}|,|c_j^{\text{novo}}|)\,}
    !!  \right)^{2}}\right)
    !!  \]
    !!  Se \( \text{err}>1 \) em qualquer ponto do bloco, o passo é **rejeitado** e
    !!  refeito com um novo \( h \):
    !!  \[
    !!  fac=\text{clamp}\!\left(\frac{facsafe}{\text{max\_err1}^{1/elo}},
    !!  [facmin,facmax]\right), \qquad
    !!  h_{novo}=\text{clamp}(h\cdot fac,\,[dt_{min},dt_{max}])
    !!  \]
    !!  Se aceito, \( t_c \mathrel{+}= h \) e `last_accepted_dt` guarda \( h_{novo} \)
    !!  como chute inicial da próxima chamada deste mesmo bloco.
    !!
    !!  ### Fatoração esparsa
    !!
    !!  A matriz \( A \) é esparsa (topologia do mecanismo RELACS); sua fatoração LU
    !!  é feita pelo pacote `sparse` (Sangiovanni-Vincentelli/Kundert, via as
    !!  rotinas `sf*`, `bind(C)`). O **padrão** de esparsidade (`ipos`/`jpos`) é
    !!  determinado uma única vez, na primeira chamada (`i==1`), e reaproveitado
    !!  em todas as chamadas seguintes — só a fatoração **numérica** é refeita a
    !!  cada ponto de grade e a cada tentativa de passo.
    subroutine rosenbrok( &
        i                  &
      , block_end          &
      , maxblock_size      &
      , sizeofmatrix       &
      , dt_min             &
      , dt_max             &
      , time_f             &
      , time_c             &
      , dt_chem             &
      , last_accepted_dt   &
      , atol               &
      , rtol               &
      , matrix_id          &
      , numberofnonzeros   &
      , ipos               &
      , jpos               &
      , element            &
    )

        implicit none

        !> índice do bloco sendo integrado nesta chamada
        integer,             intent(in)    :: i
        !> número de pontos de grade válidos em cada bloco
        integer,             intent(in)    :: block_end(:)
        !> tamanho máximo de bloco (dimensão dos arrays automáticos por ponto)
        integer,             intent(in)    :: maxblock_size
        !> tamanho da matriz química (=`nspecies`), usado por `sfcreate_solve`
        integer,             intent(in)    :: sizeofmatrix
        !> limites do passo interno de integração \( h \): \( dt_{min} \le h \le dt_{max} \)
        real(kind=c_double),    intent(in)    :: dt_min, dt_max
        !> tempo total do sub-ciclo químico a cobrir, \( t_f = dtlt\times n\_dyn\_chem \)
        real(kind=c_double),    intent(in)    :: time_f
        !> tempo interno já integrado dentro do sub-ciclo \( t_c \); avança a cada passo aceito
        real(kind=c_double),    intent(inout) :: time_c
        !> passo de integração interno \( h \); entra com o chute inicial e sai com o último \( h \) tentado
        real(kind=c_double),    intent(inout) :: dt_chem
        !> memória do último \( h \) aceito por bloco, usada como chute inicial na próxima chamada
        real,                intent(inout) :: last_accepted_dt(:)
        !> tolerâncias absoluta (`atol`) e relativa (`rtol`) por espécie, usadas em \( \text{err} \)
        real(kind=c_double),    intent(in)    :: atol(:), rtol(:)
        !> identificador da matriz esparsa \( A \) (criado só na 1ª chamada, `i==1`, reaproveitado depois)
        integer(kind = 8),   intent(inout) :: matrix_id
        !> número de elementos não-nulos da estrutura esparsa (fixado na 1ª chamada)
        integer,             intent(inout) :: numberofnonzeros
        !> posições (linha,coluna) dos elementos não-nulos de \( A \) (fixadas na 1ª chamada)
        integer,             intent(inout) :: ipos(:), jpos(:)
        !> handles dos elementos da matriz esparsa (alocados/realocados na 1ª chamada)
        integer(kind = 8), allocatable, intent(inout) :: element(:)

        !> \( \mathbf{J} \) reconstruída, na forma \( A=\frac{1}{\gamma h}\mathbf{I}-\mathbf{J} \):
        !! preenchida para todos os pontos antes do laço de pontos (a parte fora da
        !! diagonal muda só uma vez por sub-passo; a diagonal, a cada tentativa,
        !! pois depende de \( h \)); lida ponto a ponto dentro do laço de pontos.
        real(kind=c_double) :: dlmat(maxblock_size, nspecies, nspecies)
        !> lado direito do estágio \( k_1 \), \( b_1=F(\mathbf{y}_0) \); calculado uma
        !! vez por sub-passo (fora do laço de tentativas) e reaproveitado em toda
        !! rejeição/nova tentativa daquele sub-passo
        real(kind=c_double) :: dlb1(maxblock_size, nspecies)
        !> estágio \( k_1 \); escrito dentro do laço de pontos, relido depois dele
        !! (na combinação da solução \( \mathbf{y}_{n+1} \))
        real(kind=c_double) :: dlk1(maxblock_size, nspecies)
        !> estágio \( k_3 \); mesmo padrão de `dlk1`
        real(kind=c_double) :: dlk3(maxblock_size, nspecies)
        !> estágio \( k_4 \); mesmo padrão de `dlk1`, também relido na estimativa de erro
        real(kind=c_double) :: dlk4(maxblock_size, nspecies)

        !> lado direito do estágio \( k_2 \); escrátio de um ponto por vez (escrito e
        !! consumido na mesma passada do laço de pontos, sem dimensão de ponto)
        real(kind=c_double) :: dlb2(nspecies)
        !> lado direito do estágio \( k_3 \); mesmo padrão de `dlb2`
        real(kind=c_double) :: dlb3(nspecies)
        !> lado direito do estágio \( k_4 \); mesmo padrão de `dlb2`
        real(kind=c_double) :: dlb4(nspecies)
        !> estágio \( k_2 \); usado só dentro da mesma passada do laço de pontos
        !! (para montar \( b_3 \) e \( b_4 \)), nunca lido fora dela — mesmo padrão de `dlb2`
        real(kind=c_double) :: dlk2(nspecies)

        !- parâmetros do método RODAS3 (idênticos aos de chem_rodas3_dyndt)
        real(kind=c_double), parameter :: threshold1 = 0.0_c_double
        !> \( \gamma \), parâmetro do tableau de Rosenbrock: define \( A=\frac{1}{\gamma h}\mathbf{I}-\mathbf{J} \)
        real(kind=c_double), parameter :: igamma = 0.5_c_double
        !> coeficiente \( 8/3 \) usado no lado direito \( b_4 \) do estágio \( k_4 \)
        real(kind=c_double), parameter :: c83 = 8._c_double / 3._c_double
        !- parâmetros de controle adaptativo do passo (h)
        real(kind=c_double), parameter :: facmin = 0.2_c_double ! lower bound on step decrease factor (default=0.2)
        real(kind=c_double), parameter :: facmax = 6.0_c_double ! upper bound on step increase factor (default=6)
        real(kind=c_double), parameter :: facsafe = 0.9_c_double ! step by which the new step is slightly smaller
        ! than the predicted value  (default=0.9)
        !> `elo`: ordem do estimador de erro embutido (\(1/elo\) no expoente de \(fac\));
        !! `uround`: piso numérico para \( \text{err} \); `roundoff`: folga na aceitação (\( \text{err}-roundoff>1 \))
        real(kind=c_double), parameter :: uround = 1.e-15_c_double, elo = 3.0_c_double, roundoff = 1.e-8_c_double
        integer, parameter :: complex_t = 0

        !- variáveis locais (índices de laço e escalares de uma tentativa de passo)
        integer :: ji, jj, ijk, nz, blocknonzeros, error, all_accepted
        real(kind=c_double) :: igamma_dtstep, dt_chem_i, dt_new, fac, tol, err1, max_err1
        real(kind=c_double), allocatable, target :: rhs_tmp(:), sol_tmp(:)

            !> **Laço externo de sub-passos**: repete a integração em sub-intervalos de
            !! tamanho \( h \) até cobrir todo o intervalo \( [0,t_f] \), acumulando o
            !! tempo integrado em \( t_c \).
            run_until_integr_ends: do while (time_c + roundoff < time_f)

                !> Jacobiana \( \mathbf{J}=\partial F/\partial\mathbf{c} \) em
                !! \( \mathbf{c}=\mathbf{y}_0 \) (`spack%sc_p`), calculada **uma única vez**
                !! por sub-passo — reaproveitada em todas as tentativas de \( h \) deste
                !! sub-passo, mesmo que rejeitadas.
                call jacdchemdc (spack%sc_p &
                !     CALL jacdchemdc (spack%sc_p_4  &
                , spack%rk &
                , spack%dldrdc& ! Jacobian matrix
                , nspecies, 1 &
                , block_end(i) & !index_g%block_end(i), &
                , maxblock_size, nr)

                !> \( F(\mathbf{y}_0)=P-L \) em \( \mathbf{y}_0 \), também calculado uma
                !! única vez por sub-passo — vira \( b_1 \) do estágio \( k_1 \).
                call fexchem (spack%sc_p &
                !     CALL fexchem (spack%sc_p_4 &
                , spack%rk &
                , spack%dlr & !production term
                , nspecies, 1 &
                , block_end(i) & !index_g%block_end(i), &
                , maxblock_size, nr)

                !> \( b_1(ijk,j) = F_j(\mathbf{y}_0)_{ijk} \) — copiado para **todos** os
                !! pontos do bloco de uma vez, antes do laço de pontos começar.
                do ji = 1, nspecies
                    do ijk = 1, block_end(i) !index_g%block_end(i)
                        !                  PRINT *,'LFR-DBG: spack%dlr(ijk,ji)',ijk,ji,spack%dlr(ijk, ji)
                        dlb1(ijk, ji) = spack%dlr(ijk, ji)
                    enddo
                enddo

                !> Parte fora da diagonal de \( A=\frac{1}{\gamma h}\mathbf{I}-\mathbf{J} \):
                !! \( dlmat(ijk,j,k) = -\mathbf{J}_{jk} \), a mesma para todas as tentativas
                !! de \( h \) deste sub-passo (só a diagonal depende de \( h \), preenchida
                !! logo abaixo, dentro de `untilaccepted`).
                do jj = 1, nspecies
                    do ji = 1, nspecies
                        do ijk = 1, block_end(i) !index_g%block_end(i)

                            dlmat(ijk, ji, jj) = - spack%dldrdc(ijk, ji, jj)

                        enddo
                    enddo
                enddo

                !> **Laço de aceitação/rejeição do passo**: tenta integrar o sub-passo
                !! com o \( h \) corrente; se o erro estimado ultrapassar a tolerância,
                !! reduz \( h \) (via `fac`) e repete — sem recalcular \( \mathbf{J} \)
                !! nem \( F(\mathbf{y}_0) \), só a diagonal de \( A \) e sua fatoração LU.
                untilaccepted: do

                    !> Preenche a diagonal de \( A \): \( dlmat(ijk,j,j) = \frac{1}{\gamma h} - \mathbf{J}_{jj} \).
                    !! Único termo de \( A \) que muda a cada tentativa, pois depende de \( h \).
                    igamma_dtstep = 1.0_c_double / (igamma * dt_chem)

                    do jj = 1, nspecies
                        do ijk = 1, block_end(i) !index_g%block_end(i)
                            dlmat(ijk, jj, jj) = igamma_dtstep - spack%dldrdc(ijk, jj, jj)
                        enddo
                    enddo

                    !> **Padrão de esparsidade, determinado uma única vez** (primeira
                    !! chamada de `rosenbrok` no bloco `i==1`): varre \( dlmat(1,\cdot,\cdot) \)
                    !! (a matriz do primeiro ponto de grade) marcando as posições
                    !! \( (j,k) \) não-nulas em `ipos`/`jpos`, cria a estrutura esparsa
                    !! (`sfcreate_solve`), obtém os *handles* de cada elemento
                    !! (`sfgetelement`) e faz a primeira fatoração LU (`sffactor`) — tudo
                    !! isso é reaproveitado, sem refazer, em todas as chamadas seguintes.
                    if ((i .eq. 1)) then
                        blocknonzeros = 0
                        do ji = 1, nspecies
                            do jj = 1, nspecies
                                !IF(dlmat(1,Ji,Jj)<=0 .AND. dlmat(1,Ji,Jj)>=(-1)*0) CYCLE
                                if ((dlmat(1, ji, jj) .ne. 0.0_c_double)) then
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
                        end do

                        call sfzero(matrix_id)
                        do nz = 1, numberofnonzeros
                            ji = ipos(nz)
                            jj = jpos(nz)
                            call sfadd1real(element(nz), dlmat(1, ji, jj))
                        end do

                        error = sffactor(matrix_id)
                    endif

                    !> **Laço de pontos**: para cada ponto `ijk` do bloco, refatora
                    !! numericamente \( A \) com os valores de `dlmat(ijk,:,:)` (reaproveitando
                    !! o padrão/estrutura simbólica já fixada acima) e resolve os 4 estágios
                    !! do RODAS3, um após o outro, todos com a mesma matriz \( A \) já fatorada.
                    !@LNCC: begin points loop
                    do ijk = 1, block_end(i) !index_g%block_end(i)
                        call sfzero(matrix_id)
                        do nz = 1, numberofnonzeros
                            ji = ipos(nz)
                            jj = jpos(nz)
                            call sfadd1real(element(nz), dlmat(ijk, ji, jj))
                        end do
                        error = sffactor(matrix_id)

                        !---------------------------------------------------------------------------------------------------------------------
                        !> **Estágio 1**: resolve \( A\,k_1=b_1 \), com \( b_1=F(\mathbf{y}_0) \)
                        !! (já calculado, em `dlb1`). Estágio puramente linear — não avalia \(F\).
                        !-   Compute DLk1 by Solving (1/(Igamma*dt) - DLRDC) DLk1=DLR
                        !- solver sparse
                        allocate(rhs_tmp(nspecies), sol_tmp(nspecies))
                        rhs_tmp = dlb1(ijk, :)
                        call sfsolve_c(matrix_id, c_loc(rhs_tmp), c_loc(sol_tmp))
                        dlk1(ijk, :) = sol_tmp
                        deallocate(rhs_tmp, sol_tmp)
                        !
                        !---------------------------------------------------------------------------------------------------------------------
                        !> **Estágio 2**: resolve \( A\,k_2=b_2 \), com
                        !! \( b_2=\dfrac{4}{h}k_1+F(\mathbf{y}_0) \). Também linear — ainda não
                        !! avalia \(F\) de novo, só reaproveita \( k_1 \) e \( b_1 \).
                        !    2- Second step
                        !    compute   K2 by solving (1/0.5 h JAC)K2 = 4/h * K1 +  F(Yn)
                        !    compute DLK2 by solving (1/Igama*dt - DLRDC)DLK2 = 4/h * DLK1 +  DLb1
                        dt_chem_i = 1.0_c_double / dt_chem

                        do ji = 1, nspecies
                            dlb2(ji) = (4.0_c_double * dt_chem_i) * dlk1(ijk, ji) + &
                            dlb1(ijk, ji)
                        enddo
                        allocate(rhs_tmp(nspecies), sol_tmp(nspecies))
                        rhs_tmp = dlb2
                        call sfsolve_c(matrix_id, c_loc(rhs_tmp), c_loc(sol_tmp))
                        dlk2 = sol_tmp
                        deallocate(rhs_tmp, sol_tmp)
                        !---------------------------------------------------------------------------------------------------------------------
                        !> **Estágio 3**: reavalia \( F \) (não-linear). Prediz
                        !! \( \mathbf{y}_3=\mathbf{y}_0+2k_1 \) (com piso em `threshold1`,
                        !! corrigindo \( k_1 \) retroativamente se houve corte), calcula
                        !! \( F(\mathbf{y}_3) \) por `fexchem` e resolve
                        !! \( A\,k_3=b_3=F(\mathbf{y}_3)+\dfrac{1}{h}(k_1-k_2) \).
                        !    3- Third step
                        !    a) update concentrations

                        !dt_chem_i = 1.0d0/dt_chem
                        do ji = 1, nspecies
                            spack%sc_p_new(ijk, ji) = spack%sc_p(ijk, ji) + 2.0_c_double * dlk1(ijk, ji)

                            if (spack%sc_p_new(ijk, ji) .lt. threshold1) then
                                spack%sc_p_new(ijk, ji) = threshold1
                                dlk1(ijk, ji) = 0.5_c_double * (spack%sc_p_new(ijk, ji) - spack%sc_p(ijk, ji))
                            endif
                        enddo
                        !
                        !    b) update the net production term (DLr= P-L = F(Y3)) at this stage with the first-order
                        !       approximation with the new concentration
                        !
                        call fexchem (spack%sc_p_new &
                        , spack%rk &
                        , spack%dlr &
                        , nspecies, ijk &
                        , ijk & !index_g%block_end(i), &
                        , maxblock_size, nr)

                        !    c) compute   K3 by solving (1 /(0.5 h) - JAC)K3 =   F(Y3) + 0.5 (K1-K2)
                        do ji = 1, nspecies
                            dlb3(ji) = spack%dlr(ijk, ji) + dt_chem_i * &
                            (dlk1(ijk, ji) - dlk2(ji))
                        enddo
                        allocate(rhs_tmp(nspecies), sol_tmp(nspecies))
                        rhs_tmp = dlb3
                        call sfsolve_c(matrix_id, c_loc(rhs_tmp), c_loc(sol_tmp))
                        dlk3(ijk, :) = sol_tmp
                        deallocate(rhs_tmp, sol_tmp)

                        !---------------------------------------------------------------------------------------------------------------------
                        !> **Estágio 4**: reavalia \( F \) de novo. Prediz
                        !! \( \mathbf{y}_4=\mathbf{y}_0+2k_1+k_3 \) (mesmo piso/correção
                        !! retroativa de `threshold1`, agora sobre \( k_3 \)), calcula
                        !! \( F(\mathbf{y}_4) \) e resolve
                        !! \( A\,k_4=b_4=F(\mathbf{y}_4)+\dfrac{1}{h}\big(k_1-k_2-\tfrac{8}{3}k_3\big) \).
                        !    4- Fourth step
                        !    a) update concentrations
                        !       Y4 = Yn + 2 * k1 +  K3
                        dt_chem_i = 1.0_c_double / dt_chem
                        do ji = 1, nspecies
                            spack%sc_p_new(ijk, ji) = spack%sc_p(ijk, ji) + &
                            2.0_c_double * dlk1(ijk, ji) + &
                            dlk3(ijk, ji)

                            if (spack%sc_p_new(ijk, ji) .lt. threshold1) then
                                spack%sc_p_new(ijk, ji) = threshold1
                                dlk3(ijk, ji) = (spack%sc_p_new(ijk, ji) - spack%sc_p(ijk, ji)) &
                                -2.0_c_double * dlk1(ijk, ji)
                            endif

                        enddo
                        !    b) update the net production term (DLr= P-L = F(Y4) ) at this stage with the 3rd-order
                        !       approximation with the new concentration
                        !
                        call fexchem (spack%sc_p_new & ! Y4
                        , spack%rk &
                        , spack%dlr & ! F(Y4)
                        , nspecies, ijk &
                        , ijk & !index_g%block_end(i)
                        , maxblock_size, nr)

                        !    c) compute   K4 by solving (1/(0.5 h)- JAC)K4 =    F(Y3) + K1/h -K2/h -8/3 K3/h
                        do ji = 1, nspecies
                            dlb4(ji) = spack%dlr (ijk, ji) + dt_chem_i * & ! F(Y4)
                            (dlk1(ijk, ji) &
                            - dlk2(ji) &
                            - c83 * dlk3(ijk, ji))
                        enddo
                        allocate(rhs_tmp(nspecies), sol_tmp(nspecies))
                        rhs_tmp = dlb4
                        call sfsolve_c(matrix_id, c_loc(rhs_tmp), c_loc(sol_tmp))
                        dlk4(ijk, :) = sol_tmp
                        deallocate(rhs_tmp, sol_tmp)

                    enddo
                    !@LNCC: end points loop


                    !---------------------------------------------------------------------------------------------------------------------
                    !> **Solução do sub-passo**, para todos os pontos:
                    !! \( \mathbf{y}_{n+1}=\mathbf{y}_0+2k_1+k_3+k_4 \) (com piso em `threshold1`).
                    !! Ainda **candidata** — só substitui `spack%sc_p` se o passo for aceito
                    !! mais abaixo.
                    dt_chem_i = 1.0_c_double / dt_chem
                    do ji = 1, nspecies
                        do ijk = 1, block_end(i) !index_g%block_end(i)

                            spack%sc_p_new(ijk, ji) = spack%sc_p(ijk, ji) + &
                            2.0_c_double * dlk1(ijk, ji) &
                            + dlk3(ijk, ji) &
                            + dlk4(ijk, ji)
                            spack%sc_p_new(ijk, ji) = max (spack%sc_p_new(ijk, ji), threshold1)

                        enddo
                    enddo

                    !> **Estimativa de erro do sub-passo**, ponto a ponto, direto de \( k_4 \)
                    !! (sem custo de uma solução de ordem inferior separada):
                    !! \[
                    !! tol_j = atol_j + rtol_j\cdot\max(|c_j^{\text{antigo}}|,|c_j^{\text{novo}}|),
                    !! \qquad
                    !! \text{err}(ijk) = \max\!\left(\text{uround},\ \sqrt{\tfrac{1}{n_{esp}}
                    !! \textstyle\sum_j\big(k_{4,j}/tol_j\big)^2}\right)
                    !! \]
                    do ijk = 1, block_end(i) !index_g%block_end(i)

                        spack%err(ijk) = 0.0_c_double
                        do ji = 1, nspecies
                            if (spack%sc_p_new(ijk, ji) .gt. 2.5e+19_c_double) then
                                print *, 'scp e scp_new = ', spc_name(ji), "Lev = ",i,"Cell=", ijk, spack%sc_p(ijk, ji), spack%sc_p_new(ijk, ji)
                                stop 555
                            endif
                            tol = atol(ji) + rtol(ji) * dmax1(dabs(spack%sc_p(ijk, ji)), dabs(spack%sc_p_new(ijk, ji)))

                            err1 = dlk4(ijk, ji)

                            spack%err(ijk) = spack%err(ijk) + (err1 / tol)**2.0_c_double

                        enddo

                        spack%err(ijk) = dmax1(uround, dsqrt(spack%err(ijk) / nspecies))
                    enddo

                    !> Passo **aceito** só se `err ≤ 1` (com folga `roundoff`) em **todos**
                    !! os pontos do bloco — um único ponto ruim rejeita o bloco inteiro.
                    all_accepted = 1
                    do ijk = 1, block_end(i) !index_g%block_end(i)
                        if (spack%err(ijk) - roundoff > 1.0_c_double) then
                            all_accepted = 0 ;
                            exit
                        endif
                    enddo

                    !- find the maximum error occurred
                    max_err1 = maxval(spack%err(1:block_end(i))) !index_g%block_end(i) ) )

                    !> Controle adaptativo do passo, a partir do pior erro do bloco:
                    !! \[
                    !! fac=\text{clamp}\!\left(\frac{facsafe}{\text{max\_err1}^{1/elo}},
                    !! [facmin,facmax]\right), \qquad
                    !! h_{novo}=\text{clamp}(h\cdot fac,\,[dt_{min},dt_{max}])
                    !! \]
                    !- use it to determine the new time step for all block elements
                    !- new step size is bounded by FacMin <= Hnew/H <= FacMax
                    fac = min(facmax, max(facmin, facsafe / max_err1**(1.0_c_double / elo)))

                    !- possible new timestep
                    dt_new = dt_chem * fac
                    dt_new = max(dt_min, min(dt_max, dt_new))

                    !- to reset the timestep resizing in function of the error estimation, use the statements below:
                    ! all_accepted = 1; dt_new=dt_max

                    !> **Rejeita**: reduz \( h \) e repete `untilaccepted` a partir da
                    !! mesma \( \mathbf{J} \)/\( F(\mathbf{y}_0) \) — só a diagonal de \( A \)
                    !! e sua fatoração mudam. Exceção: se \( h_{novo} \) já está no piso
                    !! `dt_min`, força a aceitação para não travar a integração.
                    if (all_accepted == 0 .and. dt_new > dt_min) then  ! current solution is not accepted

                        !- resize the timestep and try again
                        dt_chem = dt_new

                    else    ! current solution is     accepted
                        !count_blocks_accept = count_blocks_accept+1
                        !> **Aceita**: avança \( t_c \), guarda \( h_{novo} \) em
                        !! `last_accepted_dt(i)` (chute inicial da próxima chamada deste
                        !! bloco) e efetiva a solução (`spack%sc_p ← spack%sc_p_new`).
                        !- go ahead, updating spack%sc_p with the solution (spack%sc_p_new)
                        !- next time
                        time_c = time_c + dt_chem

                        !- next timestep (dt_new but limited by the time_f-time_c, the rest of time integration interval)
                        dt_chem = min(dt_new, time_f - time_c)
                        !- save the accepted timestep for the next integration interval
                        if (time_c < time_f) last_accepted_dt(i) = dt_new !index_g%last_accepted_dt(i) =    dt_new

                        !- pointer (does not work yet)
                        ! spack%sc_p=>spack%sc_p_new     ! POINTER
                        !- copy
                        do ji = 1, nspecies
                            do ijk = 1, block_end(i) !index_g%block_end(i)
                                spack%sc_p(ijk, ji) = spack%sc_p_new(ijk, ji)
                            enddo
                        enddo

                        exit untilaccepted

                    endif

                end do untilaccepted


            enddo run_until_integr_ends ! time-spliting

    end subroutine rosenbrok



!--------------------------------------------------------------------------  
  subroutine get_number_nonzeros(jppj,nr,nspecies,rk,jphoto,p00,nonzeros)

    integer          , intent(in)    :: jppj
    integer          , intent(in)    :: nr
    integer          , intent(in)    :: nspecies
    real(kind=c_double) , intent(inout) :: rk(nr)
    real(kind=c_double) , intent(inout) :: jphoto(jppj)
    real             , intent(in)    :: p00
    integer          , intent(inout) :: nonzeros

    real(kind=c_double) ,dimension(nspecies,nspecies) :: def_non_zeros 
    real(kind=c_double) ,dimension(nspecies) :: sc_p
    real(kind=c_double)  :: xlw,vapp(1),cosz(1),temp(1),press(1),att(1)
    integer :: i,ji,jj

    jphoto(:) = 2.3333331_c_double
    xlw	  = 1._c_double
    vapp(:) = 1.e15_c_double
    cosz(:) = 1._c_double
    att(:)  = 1.0_c_double
    temp(:) = 273.15_c_double
    press(:)= real(p00, kind=c_double)
    sc_p    = 1.e15_c_double ! dummy concentration to get the maximum number
                    ! of possible non zero elements

    call kinetic(jppj,jphoto	 &
     		     ,rk(1:nr)   &
     		     ,temp	 &
     		     ,vapp	 &
     		     ,press	 &
     		     ,1,1,1,nr  )

    def_non_zeros = 0._c_double

    call jacdchemdc(sc_p,rk(1:nr),def_non_zeros, & ! jacobian matrix
                    nspecies,1,1,1,nr)

    nonzeros = 0
    do jj=1,nspecies
       def_non_zeros(jj,jj)=1._c_double+ def_non_zeros(jj,jj)
       do ji=1,nspecies
          if (def_non_zeros(ji,jj) .ne. 0._c_double) then
             nonzeros = nonzeros + 1
          endif
       enddo
    enddo

  end subroutine get_number_nonzeros


end module mod_chem_spack_rodas3_dyndt
