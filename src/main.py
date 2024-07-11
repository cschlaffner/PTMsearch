import logging
import subprocess
from pathlib import Path

from pyopenms import MSExperiment, MzMLFile

from src.config.config import Config
from src.diagnostic_ions.detection import DiagnosticIonDetector
from src.mzml_processing.extraction import ScanWindowExtractor
from src.mzml_processing.utils import get_diann_compatible_mzml_output_file
from src.split_processing.scan_window_splitting import ScanWindowSplitting

logger = logging.getLogger(__name__)


def main(config_path: Path):
    config = Config.from_path(config_path)
    tmp_path = Path(config.tmp_dir)

    logger.info("Loading spectra...")

    exp = MSExperiment()
    MzMLFile().load(str(config.mzml_file), exp)

    extractor = ScanWindowExtractor(
        config.lower_collision_energy, config.higher_collision_energy
    )

    # TODO: make a better function in the extractor
    # and add some validation
    ms1_and_higher_energy_windows = extractor.extract_ms1_and_higher_energy_windows(exp)
    higher_energy_windows = extractor.extract_higher_energy_windows(
        ms1_and_higher_energy_windows
    )
    ms1_windows = extractor.extract_ms1_windows(ms1_and_higher_energy_windows)

    ms1_and_lower_energy_windows = extractor.extract_ms1_and_lower_energy_windows(exp)

    logger.info("Searching for diagnostic ions in higher-energy spectra...")

    detected_ions_df = DiagnosticIonDetector(
        config.known_diagnostic_ions_file,
        config.diagnostic_ions_mass_tolerance,
        config.diagnostic_ions_mass_tolerance_unit,
        config.snr_threshold,
        config.higher_collision_energy,
    ).extract_diagnostic_ions_for_spectra(higher_energy_windows.getSpectra())

    ions_file = tmp_path / "detected_ions.csv"

    detected_ions_df.to_csv(ions_file, index=False)

    logger.info("Saved diagnostic ion detection results in %s.", ions_file)

    # HERE, BEFORE splitting: reduce the detected ions df to only include ions that are also in the SL! So that
    # the others are searched with "unmdifed!"

    windows_by_mod = ScanWindowSplitting(
        config.lower_collision_energy, config.higher_collision_energy
    ).split_windows_by_mods(
        ms1_windows,
        ms1_and_lower_energy_windows,
        ms1_and_higher_energy_windows,
        detected_ions_df,
    )

    window_mods = set(windows_by_mod.keys())
    spectral_library_mods = set(config.spectral_library_files_by_mod.keys())

    # this is trash! Redo, so that the mods from the windows that are not in the SL are searched as unmodified!
    matching_mods = window_mods.intersection(spectral_library_mods)
    window_only_mods = window_mods.difference(spectral_library_mods).difference(
        {"unmodified"}
    )
    spectral_library_only_mods = spectral_library_mods.difference(
        window_mods
    ).difference({"unmodified"})

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

    for mod in matching_mods:
        windows_for_mod = windows_by_mod[mod]
        # DIA-NN requires spectra to have incremental IDs without missing numbers
        windows_for_mod = extractor.rename_spectrum_ids(windows_for_mod)

        mzml_path_for_mod = tmp_path / f"lower_energy_windows_{mod}.mzML"
        output_file = get_diann_compatible_mzml_output_file()
        output_file.store(str(mzml_path_for_mod), windows_for_mod)

        spectral_library_file_for_mod = config.spectral_library_files_by_mod[mod]

        subprocess.run(
            f"{config.dia_nn_path} "
            f"--f {mzml_path_for_mod} "
            f"--lib {spectral_library_file_for_mod} "
            f"--out {tmp_path}/report_{mod}.tsv "
            "--mass-acc 5 --mass-acc-ms1 20 --window 0 --threads 8",
            check=True,
        )

        # and do some aggregation
