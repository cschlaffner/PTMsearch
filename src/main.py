import logging
import subprocess
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
from pyopenms import MSExperiment, MzMLFile

from src.config.config import Config
from src.diagnostic_ions.detection import DiagnosticIonDetector
from src.diagnostic_ions.summary import (
    get_detected_modifications_with_combinations,
    plot_detected_ion_combinations,
    plot_detected_ions,
)
from src.diagnostic_ions.utils import get_mod_combination_str
from src.mzml_processing.extraction import ScanWindowExtractor
from src.mzml_processing.utils import get_diann_compatible_mzml_output_file
from src.split_processing.result_aggregation import ResultAggregation
from src.split_processing.scan_window_splitting import ScanWindowSplitting
from src.split_processing.spectral_library_splitting import split_library_by_mods
from src.utils import (
    make_dia_nn_additional_param_commands,
    make_dia_nn_var_mod_commands,
)
from src.validation import (
    validate_filepath_exists,
    validate_modification_combinations,
    validate_modifications,
    validate_number_lower_higher_energy_windows,
    validate_spectral_library_files_by_mod,
    validate_spectral_library_not_empty,
)


def main(config: Config, logger: logging.Logger):
    result_path = Path(config.result_dir)

    validate_filepath_exists(Path(config.mzml_file))
    validate_filepath_exists(Path(config.dia_nn_path))

    logger.info("Loading spectra...")

    exp = MSExperiment()
    MzMLFile().load(config.mzml_file, exp)

    extractor = ScanWindowExtractor(
        config.lower_collision_energy, config.higher_collision_energy
    )

    ms1_and_higher_energy_windows = extractor.extract_ms1_and_higher_energy_windows(exp)
    higher_energy_windows = extractor.extract_higher_energy_windows(
        ms1_and_higher_energy_windows
    )
    ms1_windows = extractor.extract_ms1_windows(ms1_and_higher_energy_windows)
    ms1_and_lower_energy_windows = extractor.extract_ms1_and_lower_energy_windows(exp)
    validate_number_lower_higher_energy_windows(
        ms1_and_lower_energy_windows, ms1_and_higher_energy_windows
    )

    # ---------------------- Diagnostic ion detection -------------------

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
        plot_detected_ion_combinations(detected_ions_df, result_path)

        logger.info("Saved diagnostic ion detection results in %s.", ions_file)
    else:
        logger.info("Loading already existing ion file...")
        detected_ions_df = pd.read_csv(ions_file)

    if config.run_only_ion_detection:
        logger.info(
            "The run was configured as ion detection only, so no search is done. Exiting."
        )
        return

    # ---------------------- Determining PTMs and combinations to search -------------------

    if config.modifications_to_search != []:
        modifications_to_search = config.modifications_to_search
        modification_combinations = config.modification_combinations

        validate_modifications(modifications_to_search)
        validate_modification_combinations(
            modification_combinations, len(config.modifications_additional)
        )
    else:
        modifications_to_search, modification_combinations = (
            get_detected_modifications_with_combinations(
                detected_ions_df,
                config.detection_count_min,
                len(config.modifications_additional),
            )
        )

    modification_combinations_set = (
        frozenset.union(*modification_combinations)
        if len(modification_combinations) > 0
        else set()
    )
    all_modifications_set = set(modifications_to_search).union(
        modification_combinations_set
    )
    logger.info(
        "Keeping ions for mods %s.",
        all_modifications_set,
    )
    detected_ions_df = detected_ions_df[
        detected_ions_df["letter_and_unimod_format_mod"].isin(all_modifications_set)
    ]

    logger.info(
        "Considering single mods %s and combinations %s "
        "for window splitting and spectral library preparation.",
        modifications_to_search,
        [
            get_mod_combination_str(mod_combination)
            for mod_combination in modification_combinations
        ],
    )

    if len(config.modifications_additional) > 0:
        validate_modifications(config.modifications_additional)

        if config.library_free:
            library_log_text = "and added to spectral library prediction"
        else:
            if config.spectral_library_files_by_mod:
                library_log_text = ""
            else:
                library_log_text = (
                    "and not discarded during spectral library filtering."
                )

        logger.info(
            "Additional modifications %s will be added to search %s.",
            config.modifications_additional,
            library_log_text,
        )

    # ---------------------- Spectral library preparation -------------------

    if config.library_free:
        logger.info(
            "Library-free mode was selected. Spectral libraries are predicted by DIA-NN."
        )

        database_path = config.database_for_library_prediction
        validate_filepath_exists(Path(database_path))
        mods_for_lib_creation = np.concatenate(
            [modifications_to_search, modification_combinations, ["unmodified"]]
        )
        spectral_library_files_by_mod = {}

        for mods in mods_for_lib_creation:
            mod_combination_str = get_mod_combination_str(mods)
            library_path_for_diann = (
                result_path / f"spectral_library_{mod_combination_str}"
            )
            predicted_library_path = Path(f"{library_path_for_diann}.predicted.speclib")

            if predicted_library_path.exists():
                logger.info(
                    "Spectral library for mod(s) %s already exists.",
                    mod_combination_str,
                )
                continue

            var_mod_commands_for_lib = make_dia_nn_var_mod_commands(
                mods, config.modifications_additional
            )

            additional_param_commands_lib = make_dia_nn_additional_param_commands(
                config.dia_nn_library_params
            )

            dia_nn_library_command_for_mod = (
                [
                    f"{config.dia_nn_path}",
                    "--gen-spec-lib",
                    "--fasta-search",
                    "--predictor",
                    "--strip-unknown-mods",
                    # DIA-NN does automatic extension removal for the library path,
                    # i.e. takes the path before the last '.', so if there is a '.'
                    # in the path to the result directory, the resulting library
                    # path is messed up. That is the reason
                    # for this unnecessary .predicted here.
                    f"--out-lib {library_path_for_diann}.predicted",
                    f"--fasta {database_path}",
                ]
                + additional_param_commands_lib
                + var_mod_commands_for_lib
            )

            logger.info(
                "Running DIA-NN to create spectral library for mod(s) %s  with command %s.",
                mod_combination_str,
                dia_nn_library_command_for_mod,
            )
            logger.info(
                "The logs of the run will be found in the respective DIA-NN logfile. "
                "Starting run..."
            )
            try:
                subprocess.run(
                    dia_nn_library_command_for_mod,
                    check=True,
                )
            except subprocess.CalledProcessError as e:
                logger.exception(
                    "DIA-NN spectra library creation for mod(s) %s crashed.",
                    mod_combination_str,
                )
                raise e

            logger.info(
                "Spectral library for mod(s) %s has been created.",
                mod_combination_str,
            )

            spectral_library_files_by_mod[mods] = predicted_library_path
        logger.info("Prediction of all required spectral libraries has finished.")
    else:
        if config.spectral_library_files_by_mod:
            logger.info(
                "A spectral library per split was provided. Using those libraries for search."
            )
            validate_spectral_library_files_by_mod(
                config.spectral_library_files_by_mod,
                modifications_to_search,
                modification_combinations,
            )
            spectral_library_files_by_mod = config.spectral_library_files_by_mod
        else:
            logger.info("Splitting provided spectral library...")
            validate_filepath_exists(Path(config.spectral_library_for_filtering_path))

            library = pd.read_csv(
                config.spectral_library_for_filtering_path, delimiter="\t"
            )
            spectral_library_df_by_mod = split_library_by_mods(
                library,
                modifications_to_search,
                mod_combinations_to_search=modification_combinations,
                additional_mods_to_search=config.modifications_additional,
                logger=logger,
            )
            spectral_library_files_by_mod = {}

            for mods, library_df in spectral_library_df_by_mod.items():
                validate_spectral_library_not_empty(
                    library_df, get_mod_combination_str(mods)
                )
                library_path = (
                    result_path
                    / f"spectral_library_{get_mod_combination_str(mods)}.tsv"
                )
                library_df.to_csv(library_path, sep="\t", index=False)
                spectral_library_files_by_mod[mods] = library_path
            logger.info("Spectral library splitting has finished.")

    # ---------------------- Scan window splitting -------------------

    logger.info("Splitting scan windows by modifications ...")

    windows_by_mods = ScanWindowSplitting(
        config.lower_collision_energy, config.higher_collision_energy
    ).split_windows_by_mods(
        ms1_windows,
        ms1_and_lower_energy_windows,
        ms1_and_higher_energy_windows,
        detected_ions_df,
        modifications_to_search,
        modification_combinations,
    )

    if config.modifications_to_search != []:
        for mod in config.modifications_to_search:
            if not mod in windows_by_mods:
                logger.warning(
                    "Ions for modification %s were not detected in any windows.",
                    mod,
                )
    if config.modification_combinations != []:
        for mod_combination in config.modification_combinations:
            if not mod_combination in windows_by_mods:
                logger.warning(
                    "Ion combination %s was not detected in any windows.",
                    get_mod_combination_str(mod_combination),
                )

    logger.info(
        "Searching for the following splits: %s.",
        [
            get_mod_combination_str(mod_combination)
            for mod_combination in list(windows_by_mods.keys())
        ],
    )

    # ---------------------- DIA-NN search -------------------

    file_paths_by_mods = {}

    for mods, windows_for_mod in windows_by_mods.items():
        mod_combination_str = get_mod_combination_str(mods)
        dia_nn_report_path = result_path / f"report_{mod_combination_str}.tsv"
        spectrum_id_mapping_path = (
            result_path / f"spectrum_id_mapping_{mod_combination_str}.csv"
        )
        mzml_path = result_path / f"lower_energy_windows_{mod_combination_str}.mzML"
        file_paths_by_mods[mods] = {
            "dia_nn_report_path": dia_nn_report_path,
            "mzml_path": mzml_path,
            "spectrum_id_mapping_path": spectrum_id_mapping_path,
        }

        logger.info("Preparing search for modification(s) %s...", mod_combination_str)

        if mzml_path.exists() and spectrum_id_mapping_path.exists():
            logger.info(
                "MzML with windows to search and ID mapping file already exist."
            )
        else:
            logger.info(
                "Renaming spectrum IDs to satisfy condition "
                "of incremental IDs without missing numbers for DIA-NN."
            )
            windows_for_mod, spectrum_id_mapping = extractor.rename_spectrum_ids(
                windows_for_mod, return_id_mapping=True
            )

            output_file = get_diann_compatible_mzml_output_file()
            output_file.store(str(mzml_path), windows_for_mod)
            logger.info(
                "Saved windows to search (MS1 and lower-energy) in %s.", mzml_path
            )

            spectrum_id_mapping.to_csv(spectrum_id_mapping_path, index=False)
            logger.info(
                "Saved mapping from original to renamed spectrum IDs in %s.",
                spectrum_id_mapping_path,
            )

        spectral_library_file_for_mod = spectral_library_files_by_mod[mods]

        if dia_nn_report_path.exists():
            logger.info("DIA-NN search report already exists.")
            continue

        var_mod_commands_for_mods = make_dia_nn_var_mod_commands(
            mods, config.modifications_additional
        )

        additional_param_commands_search = make_dia_nn_additional_param_commands(
            config.dia_nn_search_params
        )

        dia_nn_command_for_mod = (
            [
                f"{config.dia_nn_path}",
                f"--f {mzml_path}",
                f"--lib {spectral_library_file_for_mod}",
                f"--out {dia_nn_report_path}",
                "--qvalue 1",
                "--decoy-report",
                "--no-prot-inf",
            ]
            + additional_param_commands_search
            + var_mod_commands_for_mods
        )

        logger.info(
            "DIA-NN command to run for modification %s is %s.",
            mod_combination_str,
            dia_nn_command_for_mod,
        )

        logger.info(
            "The logs of the run will be found in the respective DIA-NN logfile. Starting run..."
        )

        try:
            subprocess.run(
                dia_nn_command_for_mod,
                check=True,
            )
            logger.info(
                "DIA-NN run for modification %s has finished.", mod_combination_str
            )
        except subprocess.CalledProcessError as e:
            logger.exception(
                "DIA-NN run for modification %s crashed.", mod_combination_str
            )
            raise e

    # ---------------------- Aggregation -------------------

    logger.info("Starting aggregation of DIA-NN reports...")

    aggregation = ResultAggregation(config.fdr_threshold, config.normalize_cscores)
    report_all_targets_with_decoys, report_fdr_filtered = aggregation.aggregate_results(
        result_path, file_paths_by_mods
    )

    logger.info("Aggregation has finished.")

    report_all_targets_with_decoys.to_csv(
        result_path / "report_aggregated_all_targets_with_decoys.csv", index=False
    )
    report_fdr_filtered.to_csv(
        result_path / "report_aggregated_fdr_filtered.csv", index=False
    )

    logger.info(
        "Result reports were written to the specified directory %s .", result_path
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "config_path",
        type=str,
        help="Path to config JSON file that defines values for higher and lower collision energy",
    )
    args = parser.parse_args()

    config = Config.from_path(Path(args.config_path))
    result_path = Path(config.result_dir)
    result_path.mkdir(exist_ok=True)

    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=result_path / "log.txt", level=logging.INFO)

    try:
        main(config, logger)
    except BaseException as e:
        logger.exception(e)
        raise e
