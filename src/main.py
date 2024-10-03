import logging
import subprocess
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from pyopenms import MSExperiment, MzMLFile

from src.config.config import Config
from src.diagnostic_ions.detection import DiagnosticIonDetector
from src.diagnostic_ions.utils import modification_unimod_format_to_dia_nn_varmod_format
from src.mzml_processing.extraction import ScanWindowExtractor
from src.mzml_processing.utils import get_diann_compatible_mzml_output_file
from src.split_processing.scan_window_splitting import ScanWindowSplitting
from src.split_processing.spectral_library_splitting import split_library_by_mods

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def main(config_path: Path):
    config = Config.from_path(config_path)
    result_path = Path(config.result_dir)

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

    logger.info("Searching for diagnostic ions in higher-energy spectra...")

    # TODO: what about the mods that are not in UniMod? handle and/or add some validation

    detected_ions_df = DiagnosticIonDetector(
        Path(config.known_diagnostic_ions_file),
        config.diagnostic_ions_mass_tolerance,
        config.diagnostic_ions_mass_tolerance_unit,
        config.snr_threshold,
        config.higher_collision_energy,
    ).extract_diagnostic_ions_for_spectra(higher_energy_windows.getSpectra())

    ions_file = result_path / "detected_ions.csv"

    detected_ions_df.to_csv(ions_file, index=False)

    logger.info("Saved diagnostic ion detection results in %s.", ions_file)

    detected_mods = set(detected_ions_df["letter_and_unimod_format_mod"].unique())

    if config.spectral_library_files_by_mod:
        spectral_library_files_by_mod = config.spectral_library_files_by_mod
    else:
        library = pd.read_csv(
            config.spectral_library_for_filtering_path, delimiter="\t"
        )
        spectral_library_df_by_mod = split_library_by_mods(library, False)
        spectral_library_files_by_mod = {}
        for mod, library_df in spectral_library_df_by_mod.items():
            library_path = result_path / f"spectral_library_{mod}.tsv"
            library_df.to_csv(library_path, sep="\t")
            spectral_library_files_by_mod[mod] = library_path

    spectral_library_mods = set(spectral_library_files_by_mod.keys()).difference(
        {"unmodified"}
    )

    matching_mods = detected_mods.intersection(spectral_library_mods)
    window_only_mods = detected_mods.difference(spectral_library_mods)
    spectral_library_only_mods = spectral_library_mods.difference(detected_mods)

    if len(window_only_mods) > 0:
        logger.info(
            "Diagnostic ions were detected for some mods that are not in your spectral library: %s",
            window_only_mods,
        )

    if len(spectral_library_only_mods) > 0:
        logger.info(
            "Diagnostic ions for some mods in your spectral library were not detected: %s",
            spectral_library_only_mods,
        )
    # TODO: add modification handling, validate that the list of search combinations does not contain single entries

    # Split only by modifications that are also in the spectral library so that all windows with
    # other modifications are searched as 'unmodified'.
    detected_ions_df = detected_ions_df[
        detected_ions_df["letter_and_unimod_format_mod"].isin(matching_mods)
    ]

    logger.info("Splitting scan windows by modifications: %s...", matching_mods)

    windows_by_mod = ScanWindowSplitting(
        config.lower_collision_energy, config.higher_collision_energy
    ).split_windows_by_mods(
        ms1_windows,
        ms1_and_lower_energy_windows,
        ms1_and_higher_energy_windows,
        detected_ions_df,
    )

    logger.info("Searching for the following splits: %s.", list(windows_by_mod.keys()))

    for mod, windows_for_mod in windows_by_mod.items():
        logger.info("Preparing search for modification %s...", mod)

        logger.info(
            "Renaming spectrum IDs to satisfy condition "
            "of incremental IDs without missing numbers for DIA-NN."
        )
        windows_for_mod, spectrum_id_mapping = extractor.rename_spectrum_ids(
            windows_for_mod, return_id_mapping=True
        )

        # TODO: this should probably be some temporary directory later
        mzml_path_for_mod = result_path / f"lower_energy_windows_{mod}.mzML"
        output_file = get_diann_compatible_mzml_output_file()
        output_file.store(str(mzml_path_for_mod), windows_for_mod)
        logger.info(
            "Saved windows to search (MS1 and lower-energy) in %s.", mzml_path_for_mod
        )

        spectrum_id_mapping_path = result_path / f"spectrum_id_mapping_{mod}.csv"
        spectrum_id_mapping.to_csv(spectrum_id_mapping_path, index=False)
        logger.info(
            "Saved mapping from original to renamed spectrum IDs in %s.",
            spectrum_id_mapping_path,
        )

        spectral_library_file_for_mod = spectral_library_files_by_mod[mod]

        var_mod_command_for_mod = (
            ""
            if mod == "unmodified"
            else f"--var-mod {modification_unimod_format_to_dia_nn_varmod_format(mod)}"
        )

        # TODO: make more configurable/use extra DIA-NN config file
        dia_nn_command_for_mod = [
            f"{config.dia_nn_path}",
            f"--f {mzml_path_for_mod}",
            f"--lib {spectral_library_file_for_mod}",
            f"--out {result_path}/report_{mod}.tsv",
            "--mass-acc 10",
            "--mass-acc-ms1 20",
            "--window 0",
            "--threads 8",
            "--qvalue 1",
            "--pg-level 2",
            "--decoy-report",
            var_mod_command_for_mod,
        ]

        logger.info(
            "DIA-NN command to run for modification %s is %s. Starting DIA-NN run...",
            mod,
            dia_nn_command_for_mod,
        )

        subprocess.run(
            dia_nn_command_for_mod,
            check=True,
        )

        logger.info("DIA-NN run for modification %s has finished.", mod)

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
