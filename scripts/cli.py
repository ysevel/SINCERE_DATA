# sincere_data/cli.py

import sys

from . import Sincere-data_reads_filtering
from . import splashome_filter
from . import Sincere-data_outlier_predictor
from . import Sincere-data_meta_analysis
from . import Sincere-data_decontamination


COMMANDS = {
    "read_filter": Sincere-data_reads_filtering.main,
    "splashome_finder": splashome_filter.main,
    "outlier_predictor": Sincere-data_outlier_predictor.main,
    "meta_analysis": Sincere-data_meta_analysis.main,
    "decontamination": Sincere-data_decontamination.main,
}


def main():

    if len(sys.argv) < 2:
        print(
            "Usage:\n"
            "  sincere-data read_filter [options]\n"
            "  sincere-data splashome_finder [options]\n"
            "  sincere-data outlier_predictor [options]\n"
            "  sincere-data meta_analysis [options]\n"
            "  sincere-data decontamination [options]"
        )
        sys.exit(1)

    command = sys.argv[1]

    if command not in COMMANDS:
        print(f"Unknown command: {command}")
        sys.exit(1)

    # remove sub-commands name:
    sys.argv = [sys.argv[0]] + sys.argv[2:]

    COMMANDS[command]()
