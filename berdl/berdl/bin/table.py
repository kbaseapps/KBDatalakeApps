import os
import json
import argparse
from pathlib import Path
from berdl.tables.datalake_table import DatalakeTableBuilder
from berdl.genome_paths import GenomePaths
from berdl.pangenome.paths_pangenome import PathsPangenome


def main(input_params, selected_clade_member):
    scratch = input_params['_config']['scratch']

    root_pangenome = Path(scratch) / 'pangenome' / selected_clade_member
    paths_pangenome = PathsPangenome(root=root_pangenome)
    paths_root = GenomePaths(root=Path(input_params['_config']['scratch']).resolve())
    input_genomes = []
    builder = DatalakeTableBuilder(paths_root, paths_pangenome, input_genomes, True, True)
    builder.build()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run BERDL pangenome pipeline"
    )
    parser.add_argument(
        "input_params",
        help="Path to input params JSON file"
    )
    parser.add_argument(
        "selected_clade_member",
        help="clade member to build pangenome"
    )
    #  read input params
    args = parser.parse_args()
    filename_input_params = args.input_params
    _selected_clade_member = args.selected_clade_member

    if not os.path.exists(filename_input_params):
        raise FileNotFoundError(
            f"Input params file not found: {filename_input_params}"
        )

    with open(filename_input_params, 'r') as fh:
        _input_params = json.load(fh)

    main(_input_params, _selected_clade_member)
