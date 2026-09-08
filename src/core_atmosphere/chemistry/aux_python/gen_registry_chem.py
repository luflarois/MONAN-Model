#!/usr/bin/env python3
"""
gen_registry_chem.py

Lê um arquivo chem_list.f90 (estilo SPACK/BRAMS), extrai a lista de nomes de
especies quimicas declarada no array `spc_name` e gera dois arquivos de
include para o Registry.xml do MPAS/MONAN:

  - Registry_chem_scalars.inc.xml
  - Registry_chem_lbc.inc.xml

Uso:
    python3 gen_registry_chem.py chem_list.f90 [--outdir DIR]
"""

import argparse
import re
import sys
from pathlib import Path


def extract_species(f90_path: Path) -> list[str]:
    """Extrai os nomes de especies do array spc_name em um chem_list.f90.

    O array normalmente aparece como:

        CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies) :: spc_name=(/ &
         'O3  ' & !
           ,'H2O2' & !
           ...
        /)

    Retorna os nomes em maiusculas, sem espacos, na ordem original.
    """
    text = f90_path.read_text(encoding="utf-8", errors="ignore")

    # Remove comentarios de linha inteira e continuacoes '&' para simplificar
    # a busca do bloco do array spc_name.
    match = re.search(
        r"::\s*spc_name\s*=\s*\(/(.*?)/\)",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if not match:
        raise ValueError(
            "Nao foi possivel encontrar a declaracao do array 'spc_name' "
            f"em {f90_path}"
        )

    block = match.group(1)

    # Extrai todos os literais entre aspas simples (nomes das especies)
    names = re.findall(r"'([^']*)'", block)
    if not names:
        raise ValueError(f"Nenhum nome de especie encontrado em spc_name ({f90_path})")

    # Remove espacos em branco (o Fortran usa CHARACTER(LEN=8) com padding)
    species = [n.strip() for n in names if n.strip()]
    return species


def build_scalars_inc(species: list[str]) -> str:
    lines = ['<!-- Chemical Species Scalars - Spack/BRAMS -->']
    for sp in species:
        low = sp.lower()
        # NB: o parser do Registry ja prefixa "index_" sozinho ao gerar a
        # dimensao de indice do constituinte; name_in_code deve ser o nome
        # puro (sem "index_"), senao a dimensao gerada fica "index_index_<nome>".
        lines.append(f'<var name="{low}" array_group="chem" name_in_code="{low}"')
        lines.append(f'     units="ppmv" description="{sp.upper()} concentration"')
        lines.append('     packages="chemistry"')
        lines.append('     default_value="0.0"')
        lines.append('     streams="input;output"/>')
    lines.append('<!-- End of chemical species -->')
    return "\n".join(lines) + "\n"


def build_lbc_inc(species: list[str]) -> str:
    lines = ['<!-- Chemical Species - lbc- Spack/BRAMS -->']
    for sp in species:
        low = sp.lower()
        lines.append(f'<var name="lbc_{low}" name_in_code="{low}" array_group="chem"')
        lines.append('     packages="chemistry" units="ppmv"')
        lines.append('     default_value="0.0"')
        lines.append(f'     description="{sp.upper()} concentration (LBC)"/>')
    lines.append('<!-- End of chemical species -->')
    return "\n".join(lines) + "\n"


def build_scalars_tend_inc(species: list[str]) -> str:
    """Constituintes de tendencia (var_array scalars_tend). Sem isso,
    num_scalars_tend fica menor que num_scalars, causando estouro de
    array (out-of-bounds) na fisica/advecao — foi exatamente esse o bug
    identificado pelo Valgrind."""
    lines = ['<!-- Chemical Species Tendencies - Spack/BRAMS -->']
    for sp in species:
        low = sp.lower()
        lines.append(f'<var name="tend_{low}" array_group="chem" name_in_code="{low}"')
        lines.append(f'     units="ppmv s^{{-1}}" description="Tendency of {sp.upper()} concentration"')
        lines.append('     packages="chemistry"')
        lines.append('     default_value="0.0"/>')
    lines.append('<!-- End of chemical species tendencies -->')
    return "\n".join(lines) + "\n"


def build_scalars_amb_inc(species: list[str]) -> str:
    """Constituintes de 'ambient' (var_array scalars_amb), usado no
    blending de fronteira em area limitada. Mesmo motivo do tend: precisa
    ter os mesmos constituintes de 'scalars' para num_scalars_amb bater."""
    lines = ['<!-- Chemical Species Ambient (LBC blending) - Spack/BRAMS -->']
    for sp in species:
        low = sp.lower()
        lines.append(f'<var name="{low}_amb" name_in_code="{low}" array_group="chem"')
        lines.append('     packages="chemistry" units="ppmv"')
        lines.append('     default_value="0.0"')
        lines.append(f'     description="{sp.upper()} concentration increment (ambient/LBC blending)"/>')
    lines.append('<!-- End of chemical species ambient -->')
    return "\n".join(lines) + "\n"


def build_zero_fortran(species: list[str]) -> str:
    """Gera a subrotina Fortran zero_chemistry_scalars, com a lista de
    especies (indices index_<nome>) derivada diretamente do chem_list.f90,
    para nao precisar manter essa lista manualmente em dois lugares.
    """
    names_lower = [sp.lower() for sp in species]

    # Monta o array de nomes em blocos de 4 por linha, no estilo Fortran,
    # todos com o mesmo comprimento (len=8, como CHARACTER(LEN=8) do chem_list.f90)
    padded = [f"'{n:<8}'" for n in names_lower]
    rows = []
    for i in range(0, len(padded), 4):
        chunk = ", ".join(padded[i:i + 4])
        cont = " &" if i + 4 < len(padded) else " ]"
        rows.append("        " + chunk + ("," if i + 4 < len(padded) else "") + cont)

    species_array_lines = "\n".join(rows)

    return f"""!===========================================================================
! zero_chemistry_scalars.f90
!
! Gerado automaticamente por gen_registry_chem.py a partir de chem_list.f90.
! Zera todas as {len(names_lower)} especies quimicas (SPACK/BRAMS) dentro do
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
   integer, parameter :: num_species = {len(names_lower)}

   character(len=8), dimension(num_species) :: species_names = [ &
{species_array_lines}

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
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Gera os includes Registry_chem_scalars.inc.xml e "
        "Registry_chem_lbc.inc.xml a partir de um chem_list.f90"
    )
    parser.add_argument("chem_list", type=Path, help="Caminho para o chem_list.f90")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="Diretorio de saida para os arquivos .inc.xml (default: diretorio atual)",
    )
    args = parser.parse_args()

    if not args.chem_list.is_file():
        print(f"Erro: arquivo nao encontrado: {args.chem_list}", file=sys.stderr)
        return 1

    species = extract_species(args.chem_list)
    print(f"Encontradas {len(species)} especies em {args.chem_list.name}:")
    print(", ".join(species))

    args.outdir.mkdir(parents=True, exist_ok=True)

    scalars_path = args.outdir / "Registry_chem_scalars.inc.xml"
    lbc_path = args.outdir / "Registry_chem_lbc.inc.xml"
    tend_path = args.outdir / "Registry_chem_scalars_tend.inc.xml"
    amb_path = args.outdir / "Registry_chem_scalars_amb.inc.xml"
    zero_fortran_path = args.outdir / "zero_chemistry_scalars.f90"

    scalars_path.write_text(build_scalars_inc(species), encoding="utf-8")
    lbc_path.write_text(build_lbc_inc(species), encoding="utf-8")
    tend_path.write_text(build_scalars_tend_inc(species), encoding="utf-8")
    amb_path.write_text(build_scalars_amb_inc(species), encoding="utf-8")
    zero_fortran_path.write_text(build_zero_fortran(species), encoding="utf-8")

    print(f"\nGerado: {scalars_path}")
    print(f"Gerado: {lbc_path}")
    print(f"Gerado: {tend_path}")
    print(f"Gerado: {amb_path}")
    print(f"Gerado: {zero_fortran_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
