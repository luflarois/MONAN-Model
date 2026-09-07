!========================================================================================
!>  Escrita da saída química (concentrações por espécie) em NetCDF, regridada
!!  para um mapa lat/lon regular.
!!
!!  A malha do MPAS/MONAN é não-estruturada: `xlat_p(ijk)`/`xlon_p(ijk)` dão a
!!  posição real de cada célula, mas o índice `ijk` NÃO segue uma varredura
!!  lat/lon (não é sequencial). Por isso, antes de escrever, cada célula é
!!  jogada (via a lat/lon real dela, não via `ijk`) numa caixa de uma grade
!!  regular `(lon, lat)`; caixas com mais de uma célula recebem a média, e
!!  caixas sem nenhuma célula ficam com `_FillValue`.
!!
!!  Gera um arquivo por chamada, nomeado `chem_out_AAAAMMDDHHMMSS.nc`, onde
!!  `AAAAMMDDHHMMSS` é a data/hora deste timestep, calculada a partir da data
!!  inicial da rodada (`start_date`) e de `iTimestep * dt` segundos.
!!
!!  Uma variável NetCDF é criada para cada espécie de `chem_list` (`spc_name`),
!!  com dimensões (lon, lat, nVertLevels).
!!
!!  ATENÇÃO — ajustar antes de compilar:
!!   1) o `use` de `chem_vars` abaixo aponta para `modMemoryChem`, que foi a
!!      origem identificada para `chem_g` em `modChem.f90`; troque para o
!!      módulo/nome corretos se for outro;
!!   2) o `kind` dos reais (`RKIND`, aqui usando kind=8) deve bater com o
!!      usado no restante do MONAN/MPAS para `xlat_p`, `xlon_p` e `sc_p`;
!!   3) o atributo `units` de cada variável está como `'model_units'` porque
!!      não é possível saber, só a partir de `chem_list`, se `sc_p` está em
!!      ppbm (BRAMS) ou molec/cm3 (spack) neste ponto da chamada — ajuste
!!      para a unidade real de `chem_g(n)%sc_p` no momento da escrita;
!!   4) `chem_list` não tem um array de "descrição" por espécie (só
!!      `spc_name`, de 8 caracteres, e propriedades físico-químicas como
!!      `weight`); a variável NetCDF de cada espécie usa `spc_name` como
!!      nome/`long_name`, e o peso molecular (`weight`) como atributo extra;
!!   5) a resolução da grade de saída (`dlon`/`dlat`, em graus) tem default de
!!      1.0°; se sair com muitas caixas vazias (malha mais grosseira que a
!!      grade escolhida), aumente `dlon_in`/`dlat_in` na chamada. Este é um
!!      "box-average" simples (rápido, O(nCells)), não uma interpolação
!!      conservativa de verdade — para pós-processamento científico mais
!!      rigoroso, prefira uma ferramenta de remapeamento dedicada (ex.:
!!      ESMF_RegridWeightGen/xESMF, CDO remapbil/remapcon) sobre este mesmo
!!      `chem_g`/`xlat_p`/`xlon_p`;
!!   6) assume-se que o domínio não cruza a linha internacional de data
!!      (antimeridiano) na convenção de longitude usada em `xlon_p` — se
!!      cruzar, ajuste `xlon_p` para uma convenção contínua antes de chamar.
!!
!!  Requer linkar com a biblioteca netcdf-fortran (`-lnetcdff -lnetcdf`) e
!!  `use netcdf` disponível no include path.
!!
!! @date 2026-09-04
!========================================================================================
module mod_chem_output

    use netcdf

    use chem_list, only: &
        nspecies           &
      , spc_name            &
      , weight              &
      , chemical_mechanism

    !> tipo de chem_g — ver nota (1) no cabeçalho do módulo
    use modMemoryChem, only: chem_vars

    implicit none

    private
    public :: write_chem_netcdf

    integer, parameter :: RKIND = 8   !> ajustar para o kind real usado no MONAN/MPAS

contains

    !========================================================================================
    !>  Escreve `chem_g`, regridado para lat/lon, em `chem_out_AAAAMMDDHHMMSS.nc`.
    !!
    !! @param chem_g       estrutura de espécies: `chem_g(n)%sc_p(i,ijk)`,
    !!                     `n`=espécie, `i`=nível vertical, `ijk`=célula      [nspecies]
    !! @param xlat_p       latitude de cada célula, -90 a 90 (graus)          [nCells]
    !! @param xlon_p       longitude de cada célula, 0 a 360 (graus)         [nCells]
    !! @param nVertLevels  número de níveis verticais (dimensão `i` de sc_p)
    !! @param nCells       número de células horizontais (dimensão `ijk`)
    !! @param iTimestep    índice do passo de tempo do modelo (0, 1, 2, ...)
    !! @param dt           passo de tempo do modelo, em segundos (ex.: 300.)
    !! @param start_date   data/hora inicial da rodada, 'AAAAMMDDHHMMSS' (14 caracteres)
    !! @param dlon_in      [opcional] resolução da grade de saída em longitude (graus, default 1.0)
    !! @param dlat_in      [opcional] resolução da grade de saída em latitude  (graus, default 1.0)
    !========================================================================================
    subroutine write_chem_netcdf( &
        chem_g       &
      , xlat_p       &
      , xlon_p       &
      , nVertLevels  &
      , nCells       &
      , iTimestep    &
      , dt           &
      , start_date   &
      , dlon_in      &
      , dlat_in      &
    )

        implicit none

        type(chem_vars),    intent(in) :: chem_g(:)
        real,   intent(in) :: xlat_p(:)
        real,   intent(in) :: xlon_p(:)
        integer,            intent(in) :: nVertLevels
        integer,            intent(in) :: nCells
        integer,            intent(in) :: iTimestep
        real,               intent(in) :: dt
        character(len=14),  intent(in) :: start_date
        real, optional,     intent(in) :: dlon_in
        real, optional,     intent(in) :: dlat_in

        character(len=64)  :: fname
        character(len=14)  :: valid_date
        integer :: ncid
        integer :: dimid_lon, dimid_lat, dimid_lev
        integer :: varid_lon, varid_lat, varid_lev, varid_ncells
        integer :: varid_spc(nspecies)
        integer :: n, ijk, ilon, ilat

        real :: dlon, dlat
        real(kind=RKIND) :: lon0, lon1, lat0, lat1
        integer :: nlon, nlat

        real(kind=RKIND), allocatable :: lon_grid(:), lat_grid(:)
        integer,           allocatable :: ilon_cell(:), ilat_cell(:)
        integer,           allocatable :: ncells_box(:,:)
        real,              allocatable :: accum(:,:,:)   ! (nlon, nlat, nVertLevels)

        real, parameter :: fill_value = -9.99e33

        !- 0) resolução da grade de saída ------------------------------------------------
        dlon = 1.0
        dlat = 1.0
        if (present(dlon_in)) dlon = dlon_in
        if (present(dlat_in)) dlat = dlat_in

        !- 1) data/hora deste timestep = start_date + iTimestep*dt segundos ------------
        call compute_valid_date(start_date, iTimestep, dt, valid_date)

        !- 2) nome do arquivo -----------------------------------------------------------
        write(fname, '(A,A,A)') 'chem_out_', trim(valid_date), '.nc'

        !- 3) monta a grade lat/lon regular a partir da extensão real das células -------
        !-    (min/max de xlat_p/xlon_p, não valores fixos -90/90, 0/360: funciona tanto
        !-    para malha global quanto para domínio regional)
        lon0 = minval(xlon_p)
        lon1 = maxval(xlon_p)
        lat0 = minval(xlat_p)
        lat1 = maxval(xlat_p)

        nlon = max(1, nint(real(lon1 - lon0) / dlon)) + 1
        nlat = max(1, nint(real(lat1 - lat0) / dlat)) + 1

        allocate(lon_grid(nlon), lat_grid(nlat))
        do ilon = 1, nlon
            lon_grid(ilon) = lon0 + real(ilon - 1, kind=RKIND) * dlon
        end do
        do ilat = 1, nlat
            lat_grid(ilat) = lat0 + real(ilat - 1, kind=RKIND) * dlat
        end do

        !- 4) associa cada célula à caixa da grade regular mais próxima, USANDO a
        !-    lat/lon real da célula (xlat_p/xlon_p) — não o índice `ijk`, que não
        !-    tem relação nenhuma com posição espacial ------------------------------------
        allocate(ilon_cell(nCells), ilat_cell(nCells))
        allocate(ncells_box(nlon, nlat))
        ncells_box = 0

        do ijk = 1, nCells
            ilon = nint(real(xlon_p(ijk) - lon0) / dlon) + 1
            ilat = nint(real(xlat_p(ijk) - lat0) / dlat) + 1
            ilon = min(max(ilon, 1), nlon)
            ilat = min(max(ilat, 1), nlat)
            ilon_cell(ijk) = ilon
            ilat_cell(ijk) = ilat
            ncells_box(ilon, ilat) = ncells_box(ilon, ilat) + 1
        end do

        !- 5) cria o arquivo --------------------------------------------------------------
        call check( nf90_create(trim(fname), NF90_CLOBBER, ncid) )

        !- atributos globais
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'title',              'MONAN chemistry output (regridded to lat/lon)') )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'chemical_mechanism', trim(chemical_mechanism)) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'start_date',         start_date) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'valid_date',         valid_date) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'itimestep',          iTimestep) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'dt_seconds',         dt) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'regrid_method', &
            'nearest-box average from unstructured MPAS cells (see module header note 5)') )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'dlon_deg', dlon) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'dlat_deg', dlat) )

        !- 6) dimensões ---------------------------------------------------------------
        call check( nf90_def_dim(ncid, 'lon', nlon,        dimid_lon) )
        call check( nf90_def_dim(ncid, 'lat', nlat,        dimid_lat) )
        call check( nf90_def_dim(ncid, 'lev', nVertLevels, dimid_lev) )

        !- 7) coordenadas -----------------------------------------------------------------
        call check( nf90_def_var(ncid, 'lon', NF90_DOUBLE, (/ dimid_lon /), varid_lon) )
        call check( nf90_put_att(ncid, varid_lon, 'units',     'degrees_east') )
        call check( nf90_put_att(ncid, varid_lon, 'long_name', 'longitude') )

        call check( nf90_def_var(ncid, 'lat', NF90_DOUBLE, (/ dimid_lat /), varid_lat) )
        call check( nf90_put_att(ncid, varid_lat, 'units',     'degrees_north') )
        call check( nf90_put_att(ncid, varid_lat, 'long_name', 'latitude') )

        call check( nf90_def_var(ncid, 'lev', NF90_INT, (/ dimid_lev /), varid_lev) )
        call check( nf90_put_att(ncid, varid_lev, 'long_name', 'model vertical level index') )

        !- variável de diagnóstico: nº de células da malha original em cada caixa
        !- (ajuda a identificar caixas vazias/preenchidas com fill value)
        call check( nf90_def_var(ncid, 'ncells_per_box', NF90_INT, (/ dimid_lon, dimid_lat /), varid_ncells) )
        call check( nf90_put_att(ncid, varid_ncells, 'long_name', &
            'number of unstructured MPAS cells averaged into this lon/lat box') )

        !- 8) uma variável por espécie química, dimensionada (lon, lat, lev) -------------
        do n = 1, nspecies
            call check( nf90_def_var(ncid, trim(spc_name(n)), NF90_REAL, &
                                      (/ dimid_lon, dimid_lat, dimid_lev /), varid_spc(n)) )
            call check( nf90_put_att(ncid, varid_spc(n), 'long_name',        trim(spc_name(n))) )
            call check( nf90_put_att(ncid, varid_spc(n), 'molar_mass_g_mol', weight(n)) )
            call check( nf90_put_att(ncid, varid_spc(n), 'units',            'model_units') ) ! ver nota (3)
            call check( nf90_put_att(ncid, varid_spc(n), '_FillValue',       fill_value) )
        end do

        call check( nf90_enddef(ncid) )

        !- 9) escreve coordenadas e diagnóstico --------------------------------------------
        call check( nf90_put_var(ncid, varid_lon, lon_grid) )
        call check( nf90_put_var(ncid, varid_lat, lat_grid) )
        call check( nf90_put_var(ncid, varid_lev, (/ (n, n = 1, nVertLevels) /)) )
        call check( nf90_put_var(ncid, varid_ncells, ncells_box) )

        !- 10) regrid + escreve cada espécie -----------------------------------------------
        allocate(accum(nlon, nlat, nVertLevels))

        do n = 1, nspecies

            accum = 0.0

            !- acumula, célula a célula, na caixa correspondente (mapeada pela
            !- posição real xlat_p/xlon_p da célula, não pelo índice ijk)
            do ijk = 1, nCells
                accum(ilon_cell(ijk), ilat_cell(ijk), :) = accum(ilon_cell(ijk), ilat_cell(ijk), :) &
                                                          + real(chem_g(n)%sc_p(:, ijk))
            end do

            !- média nas caixas com pelo menos 1 célula; _FillValue nas vazias
            do ilat = 1, nlat
                do ilon = 1, nlon
                    if (ncells_box(ilon, ilat) > 0) then
                        accum(ilon, ilat, :) = accum(ilon, ilat, :) / real(ncells_box(ilon, ilat))
                    else
                        accum(ilon, ilat, :) = fill_value
                    endif
                end do
            end do

            call check( nf90_put_var(ncid, varid_spc(n), accum) )

        end do

        deallocate(accum)

        !- 11) fecha o arquivo -------------------------------------------------------------
        call check( nf90_close(ncid) )

        deallocate(lon_grid, lat_grid, ilon_cell, ilat_cell, ncells_box)

        print *, 'mod_chem_output: escrito ', trim(fname), ' (grade ', nlon, ' x ', nlat, ')'

    end subroutine write_chem_netcdf


    !========================================================================================
    !>  Calcula `start_date + iTimestep*dt` segundos, no formato AAAAMMDDHHMMSS.
    !========================================================================================
    subroutine compute_valid_date(start_date, iTimestep, dt, valid_date)

        implicit none

        character(len=14), intent(in)  :: start_date   ! 'AAAAMMDDHHMMSS'
        integer,            intent(in)  :: iTimestep
        real,               intent(in)  :: dt            ! segundos
        character(len=14), intent(out) :: valid_date

        integer :: yyyy, mm, dd, hh, mn, ss
        integer :: jdn0, jdn1
        integer :: total_seconds, total_days, rem_seconds, sec0

        !- decompõe a data inicial ------------------------------------------------------
        read(start_date( 1: 4), '(I4)') yyyy
        read(start_date( 5: 6), '(I2)') mm
        read(start_date( 7: 8), '(I2)') dd
        read(start_date( 9:10), '(I2)') hh
        read(start_date(11:12), '(I2)') mn
        read(start_date(13:14), '(I2)') ss

        !- data inicial -> dia juliano (inteiro) -----------------------------------------
        jdn0 = ymd_to_jdn(yyyy, mm, dd)

        !- segundos totais decorridos desde o início da rodada ---------------------------
        sec0          = hh*3600 + mn*60 + ss
        total_seconds = sec0 + nint( real(iTimestep) * dt )

        total_days  = total_seconds / 86400
        rem_seconds = mod(total_seconds, 86400)
        if (rem_seconds < 0) then
            rem_seconds = rem_seconds + 86400
            total_days  = total_days - 1
        endif

        jdn1 = jdn0 + total_days

        call jdn_to_ymd(jdn1, yyyy, mm, dd)

        hh = rem_seconds / 3600
        mn = mod(rem_seconds, 3600) / 60
        ss = mod(rem_seconds, 60)

        write(valid_date, '(I4.4,I2.2,I2.2,I2.2,I2.2,I2.2)') yyyy, mm, dd, hh, mn, ss

    end subroutine compute_valid_date


    !========================================================================================
    !>  Data gregoriana -> número de dia juliano (algoritmo de Fliegel & Van
    !!  Flandern), válido para o calendário gregoriano proléptico.
    !========================================================================================
    integer function ymd_to_jdn(yyyy, mm, dd) result(jdn)
        implicit none
        integer, intent(in) :: yyyy, mm, dd
        integer :: a, y, m

        a = (14 - mm) / 12
        y = yyyy + 4800 - a
        m = mm + 12*a - 3

        jdn = dd + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045

    end function ymd_to_jdn


    !========================================================================================
    !>  Número de dia juliano -> data gregoriana (inverso do algoritmo acima).
    !========================================================================================
    subroutine jdn_to_ymd(jdn, yyyy, mm, dd)
        implicit none
        integer, intent(in)  :: jdn
        integer, intent(out) :: yyyy, mm, dd
        integer :: a, b, c, d, e, m

        a = jdn + 32044
        b = (4*a + 3) / 146097
        c = a - (146097*b) / 4
        d = (4*c + 3) / 1461
        e = c - (1461*d) / 4
        m = (5*e + 2) / 153

        dd   = e - (153*m + 2)/5 + 1
        mm   = m + 3 - 12*(m/10)
        yyyy = 100*b + d - 4800 + m/10

    end subroutine jdn_to_ymd


    !========================================================================================
    !>  Checagem de erro padrão para chamadas da biblioteca NetCDF.
    !========================================================================================
    subroutine check(istat)
        implicit none
        integer, intent(in) :: istat
        if (istat /= nf90_noerr) then
            print *, 'ERRO NetCDF: ', trim(nf90_strerror(istat))
            stop 1
        endif
    end subroutine check

end module mod_chem_output
