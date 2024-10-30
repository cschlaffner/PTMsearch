import logging
import subprocess
from argparse import ArgumentParser
from pathlib import Path
from typing import FrozenSet, List, Union

import pandas as pd
from pyopenms import MSExperiment, MzMLFile

from src.config.config import Config
from src.diagnostic_ions.detection import DiagnosticIonDetector
from src.diagnostic_ions.summary import (
    get_detected_modifications_with_combinations,
    plot_detected_ion_combinations,
    plot_detected_ions,
)
from src.diagnostic_ions.utils import modification_unimod_format_to_dia_nn_varmod_format
from src.mzml_processing.extraction import ScanWindowExtractor
from src.mzml_processing.utils import get_diann_compatible_mzml_output_file
from src.split_processing.scan_window_splitting import ScanWindowSplitting
from src.split_processing.spectral_library_splitting import split_library_by_mods

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


# TODO: maybe move this to utils file
def make_dia_nn_var_mod_commands(
    mods: List[Union[str, FrozenSet[str]]], additional_mods: List[str]
) -> List[str]:
    # add all modifications from the combination (or single mod)
    if mods == "unmodified":
        var_mod_commands = []
    else:
        mods_list = [mods] if isinstance(mods, str) else mods
        var_mod_commands = [
            f"--var-mod {modification_unimod_format_to_dia_nn_varmod_format(mod)}"
            for mod in mods_list
        ]

    # add all mods that should be searched additionally
    var_mod_commands += [
        f"--var-mod {modification_unimod_format_to_dia_nn_varmod_format(mod)}"
        for mod in additional_mods
    ]

    return var_mod_commands


def get_mod_combination_str(mod_combination: Union[str, FrozenSet[str]]) -> str:
    if isinstance(mod_combination, str):
        return mod_combination
    return "_".join(sorted(mod_combination))


def main(config_path: Path):
    config = Config.from_path(config_path)
    result_path = Path(config.result_dir)
    # TODO. write logging to file
    logger.info("Loading spectra...")

    exp = MSExperiment()
    MzMLFile().load(str(config.mzml_file), exp)

    extractor = ScanWindowExtractor(
        config.lower_collision_energy, config.higher_collision_energy
    )

    # TODO: add some validation (including whether all paths exist etc)
    ms1_and_higher_energy_windows = extractor.extract_ms1_and_higher_energy_windows(exp)
    higher_energy_windows = extractor.extract_higher_energy_windows(
        ms1_and_higher_energy_windows
    )
    ms1_windows = extractor.extract_ms1_windows(ms1_and_higher_energy_windows)

    ms1_and_lower_energy_windows = extractor.extract_ms1_and_lower_energy_windows(exp)

    ions_file = result_path / "detected_ions.csv"

    if not ions_file.exists():
        logger.info("Searching for diagnostic ions in higher-energy spectra...")

        detected_ions_df = DiagnosticIonDetector(
            Path(config.known_diagnostic_ions_file),
            config.diagnostic_ions_mass_tolerance,
            config.diagnostic_ions_mass_tolerance_unit,
            config.snr_threshold,
            config.higher_collision_energy,
        ).extract_diagnostic_ions_for_spectra(higher_energy_windows.getSpectra())

        detected_ions_df.to_csv(ions_file, index=False)
        plot_detected_ions(detected_ions_df, result_path)
        plot_detected_ion_combinations(
            detected_ions_df, result_path, config.detection_count_percentile
        )

        logger.info("Saved diagnostic ion detection results in %s.", ions_file)
    else:
        detected_ions_df = pd.read_csv(ions_file)

    # TODO: add intermediate output of the results as plots and/or list the combinations
    # Have an option for running only until this point.

    if config.modifications_to_search != []:
        modifications_to_search = config.modifications_to_search
        modification_combinations = config.modification_combinations
        # TODO: add validation that config mod combinations are not single mods and all mods in comb. have to be listed in mod list
    else:
        modifications_to_search, modification_combinations = (
            get_detected_modifications_with_combinations(
                detected_ions_df,
                config.detection_count_percentile,
                config.detection_count_min,
            )
        )

    detected_ions_df = detected_ions_df[
        detected_ions_df["letter_and_unimod_format_mod"].isin(modifications_to_search)
    ]

    logger.info(
        "Considering only mods %s and combinations %s for window splitting and spectrum library handling.",
        modifications_to_search,
        [
            get_mod_combination_str(mod_combination)
            for mod_combination in modification_combinations
        ],
    )
    # TODO: add information about additional combinations

    if config.library_free:
        # TODO: validate that the path is there
        database_path = config.database_for_library_prediction
        mods_for_lib_creation = (
            modifications_to_search + modification_combinations + ["unmodified"]
        )
        spectral_library_files_by_mod = {}

        for mods in mods_for_lib_creation:
            mod_combination_str = get_mod_combination_str(mods)
            library_path_for_diann = (
                result_path / f"spectral_library_{mod_combination_str}"
            )
            predicted_library_path = Path(f"{library_path_for_diann}.predicted.speclib")

            if not predicted_library_path.exists():
                var_mod_commands_for_lib = make_dia_nn_var_mod_commands(
                    mods, config.modifications_additions
                )

                # TODO: make more configurable/use extra DIA-NN config file
                dia_nn_library_command_for_mod = [
                    f"{config.dia_nn_path}",
                    "--gen-spec-lib",
                    "--fasta-search",
                    "--predictor",
                    "--cut K*,R*",
                    "--threads 8",
                    "--min-pep-len 7",
                    "--max-pep-len 30",
                    "--max-pr-charge 4",
                    "--min-pr-mz 350",
                    "--max-pr-mz 1200",
                    "--min-fr-mz 200",
                    "--max-fr-mz 2000",
                    "--missed-cleavages 2",
                    "--strip-unknown-mods",
                    "--var-mods 3",
                    "--pg-level 2",
                    # DIA-NN does automatic extension removal for the library path,
                    # i.e. takes the path before the last '.', so if there is a '.'
                    # in the path to the result directory, the resulting library
                    # path is messed up. That is the reason
                    # for this unnecessary .predicted here.
                    f"--out-lib {library_path_for_diann}.predicted",
                    f"--fasta {database_path}",
                ] + var_mod_commands_for_lib

                logger.info(
                    "Running DIA-NN to create spectral library for mod(s) %s  with command %s...",
                    mod_combination_str,
                    dia_nn_library_command_for_mod,
                )
                subprocess.run(
                    dia_nn_library_command_for_mod,
                    check=True,
                )
            spectral_library_files_by_mod[mods] = predicted_library_path
    else:
        if config.spectral_library_files_by_mod:
            # TODO: add validation that the keys match the mods (and combinations) to search and contains one "unmodified"
            spectral_library_files_by_mod = config.spectral_library_files_by_mod
        else:
            library = pd.read_csv(
                config.spectral_library_for_filtering_path, delimiter="\t"
            )
            # TODO: add including configurable mods (to search and to ignore)
            # and throwing an error or at least a notification if the SL for one of the specified mods/combs is empty/contains only unmod spectra
            spectral_library_df_by_mod = split_library_by_mods(
                library, False, modification_combinations
            )
            spectral_library_files_by_mod = {}
            for mods, library_df in spectral_library_df_by_mod.items():
                library_path = (
                    result_path
                    / f"spectral_library_{get_mod_combination_str(mods)}.tsv"
                )
                library_df.to_csv(library_path, sep="\t")
                spectral_library_files_by_mod[mods] = library_path

    logger.info("Splitting scan windows by modifications ...")

    windows_by_mods = ScanWindowSplitting(
        config.lower_collision_energy, config.higher_collision_energy
    ).split_windows_by_mods(
        ms1_windows,
        ms1_and_lower_energy_windows,
        ms1_and_higher_energy_windows,
        detected_ions_df,
        modification_combinations,
    )

    logger.info(
        "Searching for the following splits: %s.",
        [
            get_mod_combination_str(mod_combination)
            for mod_combination in list(windows_by_mods.keys())
        ],
    )

    for mods, windows_for_mod in windows_by_mods.items():
        mod_combination_str = get_mod_combination_str(mods)
        logger.info("Preparing search for modification(s) %s...", mod_combination_str)

        logger.info(
            "Renaming spectrum IDs to satisfy condition "
            "of incremental IDs without missing numbers for DIA-NN."
        )
        windows_for_mod, spectrum_id_mapping = extractor.rename_spectrum_ids(
            windows_for_mod, return_id_mapping=True
        )

        # TODO: this should probably be some temporary directory later
        mzml_path_for_mod = (
            result_path / f"lower_energy_windows_{mod_combination_str}.mzML"
        )
        output_file = get_diann_compatible_mzml_output_file()
        output_file.store(str(mzml_path_for_mod), windows_for_mod)
        logger.info(
            "Saved windows to search (MS1 and lower-energy) in %s.", mzml_path_for_mod
        )

        spectrum_id_mapping_path = (
            result_path / f"spectrum_id_mapping_{mod_combination_str}.csv"
        )
        spectrum_id_mapping.to_csv(spectrum_id_mapping_path, index=False)
        logger.info(
            "Saved mapping from original to renamed spectrum IDs in %s.",
            spectrum_id_mapping_path,
        )

        spectral_library_file_for_mod = spectral_library_files_by_mod[mods]

        var_mod_commands_for_mods = make_dia_nn_var_mod_commands(
            mods, config.modifications_additions
        )

        # TODO: make more configurable/use extra DIA-NN config file
        dia_nn_command_for_mod = [
            f"{config.dia_nn_path}",
            f"--f {mzml_path_for_mod}",
            f"--lib {spectral_library_file_for_mod}",
            f"--out {result_path}/report_{mod_combination_str}.tsv",
            "--mass-acc 10",
            "--mass-acc-ms1 20",
            "--window 0",
            "--threads 8",
            "--qvalue 1",
            "--pg-level 2",
            "--decoy-report",
            "--no-prot-inf",
        ] + var_mod_commands_for_mods

        logger.info(
            "DIA-NN command to run for modification %s is %s. Starting DIA-NN run...",
            mod_combination_str,
            dia_nn_command_for_mod,
        )

        try:
            subprocess.run(
                dia_nn_command_for_mod,
                # TODO: print warning if it fails
                check=True,
            )
            logger.info(
                "DIA-NN run for modification %s has finished.", mod_combination_str
            )
        except subprocess.CalledProcessError as e:
            # TODO: capture error message
            logger.warning(
                "DIA-NN run for modification %s crashed.", mod_combination_str
            )

        # and do some aggregation


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "config_path",
        type=str,
        help="Path to config JSON file that defines values for higher and lower collision energy",
    )
    args = parser.parse_args()

    main(Path(args.config_path))
