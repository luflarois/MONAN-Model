!===============================================================================
! modTest.f90
!
! Modulo com duas subrotinas:
!
!   1) readMergedChemFile(filename, nSpecies)
!        Le o arquivo binario gerado por merge_chem_files.f90 (registros:
!        (time,nponts) e depois, por ponto, (k_,glat,glon,press,temp,vapp) |
!        jphoto(:) | scp(:) | sct_dyn(:)) e guarda k_, glat, glon, scp e
!        sct_dyn em variaveis globais ao modulo (g_k, g_glat, g_glon, g_scp,
!        g_sctDyn).
!
!   2) mapChemToGridNearest(xlat, xlon, nCells, nVertLevels, nSpecies)
!        Usa os dados lidos pela rotina (1) e, para cada celula do grid de
!        destino (xlat,xlon), encontra a coluna mais proxima ("nearest
!        point") entre os dados lidos e copia os valores verticais (scp,
!        sct_dyn) para chem_g(nSpecies)%sc_p(nVertLevels,nCells) e
!        chem_g(nSpecies)%sc_t_dyn(nVertLevels,nCells), sem nenhuma
!        interpolacao (nem horizontal nem vertical) - os valores de saida
!        sao EXATAMENTE os valores lidos na posicao/coluna escolhida.
!
! PREMISSAS IMPORTANTES (leia antes de usar):
!
!   - Os niveis verticais de cada coluna sao ordenados por k_ (ordem
!     CRESCENTE), que agora vem gravado no arquivo. Se no seu MONAN/BRAMS
!     k_=1 e a superficie e k_ cresce para cima, o nivel 1 da saida
!     (scp_out(1,:)) sera a superficie; se a convencao for invertida, troque
!     "res = (k1 < k2)" por "res = (k1 > k2)" em ltTriple mais abaixo.
!
!   - "Nearest point" e calculado por distancia euclidiana simples em
!     (glat,glon) (nao e haversine/grande-circulo). Para grades regionais
!     pequenas isso e uma aproximacao razoavel; para dominios muito
!     grandes/latitudes altas, troque ltTriple/dist2 por uma metrica
!     esferica se precisar de mais precisao.
!
!   - chem_g e um array alocavel de um tipo derivado ja declarado em
!     modMemoryChem (module procedure/variable "chem_g"), com componentes
!     alocaveis sc_p(:,:) e sc_t_dyn(:,:). Este modulo apenas usa
!     chem_g via "use", nao redefine o tipo.
!
!   - Numero de niveis de saida = MAX(nVertLevels, maior numero de niveis
!     encontrado em qualquer coluna de entrada), conforme pedido. Isso
!     garante que nenhum dado lido seja truncado. Colunas com menos
!     niveis do que nLevOut ficam com os niveis excedentes preenchidos
!     com 0.0 (ajuste se preferir outro valor de preenchimento).
!
! CORRECAO (2026-08-27): os valores lidos vinham como lixo (nponts gigante,
! time = denormal) porque os arquivos chem_in_*.bin/chem_merged_*.bin sao
! escritos em real(kind=RKIND) pelo MONAN/BRAMS (CONFIRMADO -DSINGLE_PRECISION
! => precisao simples, 4 bytes), mas este modulo lia tudo como "real" com o
! kind errado - isso desalinha os bytes de TODOS os campos seguintes do
! registro. Agora todo real ligado ao arquivo usa explicitamente kind=RKIND
! (via "use mpas_kind_types, only: RKIND"). Byte-order (endianness) tambem
! e fixado explicitamente como 'BIG_ENDIAN' no open(), para casar com
! OUTPUT_CONVERT do merge_chem_files.f90, independente da flag global de
! compilacao do MONAN.
!
! CORRECAO (2026-08-28): merge_chem_files.f90 passou a gravar k_ no primeiro
! registro de cada ponto (k_,glat,glon,press,temp,vapp) e a filtrar so a
! regiao sudeste do Brasil. Este modulo foi atualizado para ler k_ e usa-lo
! para ordenar os niveis verticais de cada coluna corretamente (antes a
! ordem vinha apenas da ordem de escrita no arquivo, sem garantia de bater
! com a ordem fisica dos niveis).
!
! CORRECAO (2026-08-29): xlat_p/xlon_p do MONAN vem com longitude em 0..360,
! mas g_glon (lido do arquivo fundido) ja esta em -180..+180 (o
! merge_chem_files.f90 normaliza antes de gravar). mapChemToGridNearest
! agora normaliza uma copia local de xlon para -180..+180 (funcao lonTo180)
! antes de calcular a distancia/achar a coluna mais proxima - sem isso, o
! nearest-neighbor buscava na posicao errada para qualquer xlon > 180.
!===============================================================================
module modTest

    use mpas_kind_types, only: RKIND  ! MESMO RKIND usado pelo resto do MONAN -
                                       ! confirmado -DSINGLE_PRECISION =>
                                       ! selected_real_kind(6) (4 bytes).
                                       ! Se o modulo que define RKIND no seu
                                       ! MONAN tiver outro nome, ajuste este USE.
    use modMemoryChem, only: chem_g   ! chem_g(:) ja declarado em modMemoryChem,
                                       ! com campos sc_p(:,:) e sc_t_dyn(:,:)
    use chem_list

    implicit none
    private

    public :: readMergedChemFile
    public :: mapChemToGridNearest

    !--------------------------- dados globais do modulo -----------------------
    integer(kind=8),      save :: g_nPts     = 0_8
    real(kind=RKIND),     save :: g_time     = -1.0_RKIND
    integer,              save :: g_nSpecies = 0
    integer,          allocatable, save :: g_k(:)         ! (g_nPts) indice vertical original
    real(kind=RKIND), allocatable, save :: g_glat(:)      ! (g_nPts)
    real(kind=RKIND), allocatable, save :: g_glon(:)      ! (g_nPts)
    real(kind=RKIND), allocatable, save :: g_scp(:,:)     ! (g_nSpecies, g_nPts)
    real(kind=RKIND), allocatable, save :: g_sctDyn(:,:)  ! (g_nSpecies, g_nPts)
    !-----------------------------------------------------------------------------

contains

    !===========================================================================
    ! Le o arquivo fundido e armazena k_, glat, glon, scp e sct_dyn nas
    ! variaveis globais do modulo. press, temp, vapp e jphoto sao descartados
    ! na leitura.
    !===========================================================================
    subroutine readMergedChemFile(filename, nSpecies)
        implicit none
        character(len=*), intent(in) :: filename
        integer,           intent(in) :: nSpecies

        integer :: unit_in, ios
        integer(kind=8) :: p
        integer :: k_r
        real(kind=RKIND) :: glat_r, glon_r

        unit_in = 40
        open(unit=unit_in, file=trim(filename), form='unformatted', &
             access='sequential', status='old', action='read', &
             iostat=ios)
        if (ios /= 0) then
            print *, 'ERRO (modTest/readMergedChemFile): nao consegui abrir ', &
                     trim(filename)
            stop 1
        end if

        ! registro 1: time, nponts (gravado assim pelo merge_chem_files.f90)
        read(unit_in) g_time, g_nPts
        g_nSpecies = nSpecies

        if (allocated(g_k))      deallocate(g_k)
        if (allocated(g_glat))   deallocate(g_glat)
        if (allocated(g_glon))   deallocate(g_glon)
        if (allocated(g_scp))    deallocate(g_scp)
        if (allocated(g_sctDyn)) deallocate(g_sctDyn)

        allocate(g_k(g_nPts))
        allocate(g_glat(g_nPts))
        allocate(g_glon(g_nPts))
        allocate(g_scp(g_nSpecies, g_nPts))
        allocate(g_sctDyn(g_nSpecies, g_nPts))

        do p = 1_8, g_nPts
            ! registro do ponto 1: (k_,glat,glon,press,temp,vapp) - le so os
            ! 3 primeiros valores; press,temp,vapp (resto do registro) sao
            ! descartados automaticamente ao ler menos itens do que foram
            ! gravados naquele write.
            read(unit_in) k_r, glat_r, glon_r
            g_k(p)    = k_r
            g_glat(p) = glat_r
            g_glon(p) = glon_r

            ! registro do ponto 2: jphoto(:) - nao interessa aqui, pula o
            ! registro inteiro (read sem lista de variaveis avanca 1 registro).
            read(unit_in)

            ! registro do ponto 3: scp(1:nSpecies)
            read(unit_in) g_scp(1:g_nSpecies, p)

            ! registro do ponto 4: sct_dyn(1:nSpecies)
            read(unit_in) g_sctDyn(1:g_nSpecies, p)
        end do

        close(unit_in)

        print *, 'modTest: li ', g_nPts, ' pontos de ', trim(filename), &
                 ' (time=', g_time, ', nSpecies=', g_nSpecies, ')'

    end subroutine readMergedChemFile


    !===========================================================================
    ! Mapeia (nearest point, sem interpolacao horizontal ou vertical) os dados
    ! lidos por readMergedChemFile para o grid de destino (xlat,xlon), gravando
    ! em chem_g(nSpecies)%sc_p(nVertLevels,nCells) e
    ! chem_g(nSpecies)%sc_t_dyn(nVertLevels,nCells).
    !===========================================================================
    subroutine mapChemToGridNearest(xlat, xlon, nCells, nVertLevels, temp, press, qv, z, nSpecies)
        implicit none
        integer, intent(in) :: nCells
        integer, intent(in) :: nVertLevels
        integer, intent(in) :: nSpecies
        real(kind=RKIND), intent(in) :: xlat(nCells)
        real(kind=RKIND), intent(in) :: xlon(nCells)
        real, intent(in) :: temp(nVertlevels,nCells)
        real, intent(in) :: press(nVertlevels,nCells)
        real, intent(in) :: qv(nVertlevels,nCells)
        real,dimension(nVertLevels,nCells),intent(in) :: z

        integer(kind=8), allocatable :: idx(:)        ! indices dos pontos, ordenados por (glat,glon)
        integer(kind=8), allocatable :: colStart(:)    ! colStart(c):colStart(c+1)-1 = pontos da coluna c em idx
        integer(kind=8), allocatable :: colFirstPt(:)  ! um ponto representante de cada coluna (para o nearest search)
        integer :: nCol, maxLevSrc, nLevOut, nLevHere
        integer :: iCell, iCol, isp, lev, bestCol
        integer(kind=8) :: p, iStart
        real(kind=RKIND) :: dist2, bestDist2
        real(kind=RKIND), allocatable :: xlonN(:)   ! xlon normalizado p/ -180..+180
        real, parameter :: moleculas = 2.5e+19
        real(kind=RKIND) :: volmol,factorH2O

        double precision, parameter :: pmar = 28.96d0

        if (g_nPts <= 0_8) then
            print *, 'ERRO (modTest/mapChemToGridNearest): nenhum dado lido.'
            print *, '   Chame readMergedChemFile antes desta rotina.'
            stop 1
        end if
        if (nSpecies /= g_nSpecies) then
            print *, 'ERRO (modTest/mapChemToGridNearest): nSpecies (', nSpecies, &
                     ') difere do nSpecies usado na leitura (', g_nSpecies, ').'
            stop 1
        end if

        !-----------------------------------------------------------------------
        ! 0) xlon (do modelo) vem em 0..360, mas g_glon (lido do arquivo
        !    fundido) ja esta em -180..+180 (o merge_chem_files.f90 normaliza
        !    antes de gravar). Normaliza uma copia local de xlon para o mesmo
        !    intervalo antes de qualquer comparacao/distancia.
        !-----------------------------------------------------------------------
        allocate(xlonN(nCells))
        do iCell = 1, nCells
            xlonN(iCell) = lonTo180(xlon(iCell))
        end do

        !-----------------------------------------------------------------------
        ! 1) Ordena os indices dos pontos por (glat,glon,k_) para agrupar, de
        !    forma eficiente (O(n log n)), os pontos que pertencem a mesma
        !    coluna (mesma posicao horizontal), ja em ordem crescente de k_
        !    (nivel vertical) dentro de cada coluna.
        !-----------------------------------------------------------------------
        allocate(idx(g_nPts))
        do p = 1_8, g_nPts
            idx(p) = p
        end do
        call quicksortIdx(idx, 1_8, g_nPts)

        allocate(colStart(g_nPts + 1_8))
        allocate(colFirstPt(g_nPts))
        nCol = 0
        maxLevSrc = 0
        p = 1_8
        do while (p <= g_nPts)
            nCol = nCol + 1
            colStart(nCol) = p
            colFirstPt(nCol) = idx(p)
            iStart = p
            do while (p + 1_8 <= g_nPts)
                if (g_glat(idx(p+1_8)) == g_glat(idx(iStart)) .and. &
                    g_glon(idx(p+1_8)) == g_glon(idx(iStart))) then
                    p = p + 1_8
                else
                    exit
                end if
            end do
            maxLevSrc = max(maxLevSrc, int(p - iStart + 1_8))
            p = p + 1_8
        end do
        colStart(nCol+1) = g_nPts + 1_8

        !-----------------------------------------------------------------------
        ! 2) Numero de niveis de saida = o MAIOR entre nVertLevels (pedido
        !    pelo chamador) e o maior numero de niveis encontrado na entrada,
        !    para nao truncar nenhum dado lido.
        !-----------------------------------------------------------------------
        nLevOut = max(nVertLevels, maxLevSrc)

        !if (.not. allocated(chem_g)) allocate(chem_g(nSpecies))

        do isp = 1, nSpecies
        !    if (allocated(chem_g(isp)%sc_p))     deallocate(chem_g(isp)%sc_p)
        !    if (allocated(chem_g(isp)%sc_t_dyn))  deallocate(chem_g(isp)%sc_t_dyn)
        !    allocate(chem_g(isp)%sc_p(nLevOut, nCells))
        !    allocate(chem_g(isp)%sc_t_dyn(nLevOut, nCells))
            chem_g(isp)%sc_p     = 0.0_RKIND
            chem_g(isp)%sc_t_dyn = 0.0_RKIND
        end do




        !-----------------------------------------------------------------------
        ! 3) Para cada celula do grid de destino, acha a coluna de entrada mais
        !    proxima (nearest neighbor em glat/glon, busca exaustiva sobre as
        !    nCol colunas unicas) e copia os valores lidos, nivel a nivel, sem
        !    nenhuma interpolacao.
        !-----------------------------------------------------------------------
        do iCell = 1, nCells
            bestCol   = 1
            bestDist2 = huge(1.0)
            do iCol = 1, nCol
                p = colFirstPt(iCol)
                dist2 = (g_glat(p) - xlat(iCell))**2 + (g_glon(p) - xlonN(iCell))**2
                if (dist2 < bestDist2) then
                    bestDist2 = dist2
                    bestCol   = iCol
                end if
            end do

            nLevHere = int(colStart(bestCol+1) - colStart(bestCol))

            do lev = 1, min(nLevHere, nLevOut)
                factorH2O = water_vapor_factor(z(lev,iCell))
                !if (iCell==1) print *,'Factor=',lev,z(lev,iCell),factorH2O
                p = idx(colStart(bestCol) + int(lev-1, kind=8))
                volmol = (6.02d23 * 1d-15 * pmar) * (press(lev,iCell)) / (8.314d0 * temp(lev,icell))
                do isp = 1, nSpecies
                    chem_g(isp)%sc_p(lev, iCell)     = g_scp(isp, p)
                    chem_g(isp)%sc_t_dyn(lev, iCell) = g_sctDyn(isp, p)
                end do
                !preenchendo os gases background 
                chem_g(O2)%sc_p = moleculas*0.21*weight(O2)/volmol
                chem_g(N2)%sc_p = moleculas*0.78*weight(N2)/volmol
                chem_g(CO2)%sc_p = moleculas*0.0043*weight(CO2)/volmol
                chem_g(H2O)%sc_p(lev,:) = (moleculas*0.02*weight(H2O)/volmol)*factorH2O
            end do
        end do

        deallocate(idx, colStart, colFirstPt, xlonN)

        print *, 'modTest: mapeamento concluido - nCol=', nCol, &
                 ' maxLevSrc=', maxLevSrc, ' nLevOut=', nLevOut
        print *,'Amostra dos gases interpolados:'
        do lev=1,min(nLevHere, nLevOut),5
            print *,'lev=',lev
            print *,'O2  = ',(chem_g(O2)%sc_p(lev, iCell),iCell=1,5)
            print *,'O3  = ',(chem_g(O3)%sc_p(lev, iCell),iCell=1,5)
            print *,'N2  = ',(chem_g(N2)%sc_p(lev, iCell),iCell=1,5)
            print *,'CO2 = ',(chem_g(CO2)%sc_p(lev, iCell),iCell=1,5)
            print *,'H2O = ',(chem_g(H2O)%sc_p(lev, iCell),iCell=1,5)
            print *,'NO2 = ',(chem_g(NO2)%sc_p(lev, iCell),iCell=1,5)
            print *,'--------------------------------------------------------------------------------'
        end do


    end subroutine mapChemToGridNearest


    !===========================================================================
    ! Normaliza uma longitude para a convencao -180 a +180. Usada para
    ! converter xlon (0..360, convencao do modelo) para o mesmo intervalo
    ! usado por g_glon (ja normalizado pelo merge_chem_files.f90).
    !===========================================================================
    pure function lonTo180(lonIn) result(lonOut)
        real(kind=RKIND), intent(in) :: lonIn
        real(kind=RKIND) :: lonOut

        lonOut = lonIn
        if (lonOut > 180.0_RKIND) then
            lonOut = lonOut - 360.0_RKIND
        else if (lonOut < -180.0_RKIND) then
            lonOut = lonOut + 360.0_RKIND
        end if
    end function lonTo180


    !===========================================================================
    ! Quicksort (indices) por (glat,glon,k_), usado apenas internamente para
    ! agrupar pontos da mesma coluna, ja ordenados por nivel vertical.
    !===========================================================================
    recursive subroutine quicksortIdx(idx, loIn, hiIn)
        implicit none
        integer(kind=8), intent(inout) :: idx(:)
        integer(kind=8), intent(in)    :: loIn, hiIn

        integer(kind=8) :: i, j, tmp, pivotPt
        real(kind=RKIND) :: pivotLat, pivotLon
        integer :: pivotK

        if (loIn >= hiIn) return

        i = loIn
        j = hiIn
        pivotPt  = idx((loIn + hiIn) / 2_8)
        pivotLat = g_glat(pivotPt)
        pivotLon = g_glon(pivotPt)
        pivotK   = g_k(pivotPt)

        do while (i <= j)
            do while (ltTriple(g_glat(idx(i)), g_glon(idx(i)), g_k(idx(i)), &
                                pivotLat, pivotLon, pivotK))
                i = i + 1_8
            end do
            do while (ltTriple(pivotLat, pivotLon, pivotK, &
                                g_glat(idx(j)), g_glon(idx(j)), g_k(idx(j))))
                j = j - 1_8
            end do
            if (i <= j) then
                tmp    = idx(i)
                idx(i) = idx(j)
                idx(j) = tmp
                i = i + 1_8
                j = j - 1_8
            end if
        end do

        if (loIn < j)  call quicksortIdx(idx, loIn, j)
        if (i < hiIn)  call quicksortIdx(idx, i, hiIn)

    end subroutine quicksortIdx


    !===========================================================================
    ! Comparacao lexicografica (lat,lon,k) - usada pelo quicksortIdx. k_ e o
    ! criterio de desempate DENTRO da mesma coluna (glat,glon), garantindo
    ! que os niveis saiam em ordem crescente de k_.
    !===========================================================================
    pure logical function ltTriple(lat1, lon1, k1, lat2, lon2, k2) result(res)
        implicit none
        real(kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2
        integer,           intent(in) :: k1, k2

        if (lat1 < lat2) then
            res = .true.
        else if (lat1 > lat2) then
            res = .false.
        else
            if (lon1 < lon2) then
                res = .true.
            else if (lon1 > lon2) then
                res = .false.
            else
                res = (k1 < k2)   ! <-- se k_=1 for o TOPO em vez da
                                  !     superficie no seu MONAN, troque para (k1 > k2)
            end if
        end if
    end function ltTriple


    real function water_vapor_factor(z) result(factor)
        implicit none
        ! z : altitude em metros
        ! Retorna o fator multiplicador em relação ao valor de superfície (camada 0 km = 1.0)
        ! Baseado na tabela do perfil vertical aproximado de vapor d'água
    
        real, intent(in) :: z
    
        integer, parameter :: n = 17
        real, dimension(n) :: alt_km, ppm, fator
        real :: z_km
        integer :: i
    
        ! Altitudes da tabela (km)
        alt_km = (/ 0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 6.0d0, 8.0d0, 10.0d0, &
                    12.0d0, 16.0d0, 20.0d0, 30.0d0, 50.0d0, 80.0d0, 90.0d0, &
                    100.0d0, 150.0d0 /)
    
        ! Concentrações da tabela (ppm)
        ppm = (/ 15000.0d0, 8500.0d0, 4500.0d0, 2200.0d0, 900.0d0, 300.0d0, &
                 90.0d0, 25.0d0, 5.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 9.0d0, &
                 6.0d0, 2.0d0, 0.3d0 /)
    
        ! Fator relativo à superfície (ppm(i) / ppm(1))
        fator = ppm / ppm(1)
    
        z_km = z / 1000.0d0
    
        ! Fora dos limites da tabela: mantém o valor da borda (clamp)
        if (z_km <= alt_km(1)) then
            factor = fator(1)
            return
        else if (z_km >= alt_km(n)) then
            factor = fator(n)
            return
        end if
    
        ! Busca o intervalo [alt_km(i), alt_km(i+1)] que contém z_km
        do i = 1, n - 1
            if (z_km >= alt_km(i) .and. z_km <= alt_km(i+1)) then
                ! Interpolação linear entre os dois pontos
                factor = fator(i) + (fator(i+1) - fator(i)) * &
                          (z_km - alt_km(i)) / (alt_km(i+1) - alt_km(i))
                return
            end if
        end do
    
    end function water_vapor_factor

end module modTest
