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
        lines.append(f'<var name="{low}" array_group="chem" name_in_code="index_{low}"')
        lines.append(f'     units="ppmv" description="{sp.upper()} concentration"')
        lines.append('     packages="chemistry"')
        lines.append('     streams="input;output"/>')
    lines.append('<!-- End of chemical species -->')
    return "\n".join(lines) + "\n"


def build_lbc_inc(species: list[str]) -> str:
    lines = ['<!-- Chemical Species - lbc- Spack/BRAMS -->']
    for sp in species:
        low = sp.lower()
        lines.append(f'<var name="lbc_{low}" name_in_code="{low}" array_group="chem"')
        lines.append('     packages="chemistry" units="ppmv"')
        lines.append(f'     description="{sp.upper()} concentration (LBC)"/>')
    lines.append('<!-- End of chemical species -->')
    return "\n".join(lines) + "\n"


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

    scalars_path.write_text(build_scalars_inc(species), encoding="utf-8")
    lbc_path.write_text(build_lbc_inc(species), encoding="utf-8")

    print(f"\nGerado: {scalars_path}")
    print(f"Gerado: {lbc_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
