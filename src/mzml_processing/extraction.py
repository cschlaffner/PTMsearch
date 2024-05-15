from argparse import ArgumentParser
from pathlib import Path

from pyopenms import MSExperiment, MSSpectrum, MzMLFile

from src.config.config import Config
from src.mzml_processing.utils import (
    get_diann_compatible_mzml_output_file,
    get_spectrum_collision_energy,
)

OUTPUT_SUFFIX_LOWER_ENERGY = "_lower_energy.mzML"


def check_ms1_or_low_energy_spectrum(
    spectrum: MSSpectrum, lower_collision_energy: float
) -> bool:
    return (
        spectrum.getMSLevel() == 1
        or get_spectrum_collision_energy(spectrum) == lower_collision_energy
    )


def extract_lower_energy_windows(
    exp: MSExperiment, lower_collision_energy: float
) -> MSExperiment:
    """Extracts all lower-energy and MS1 scan windows to a second MSExperiment"""
    spectra = exp.getSpectra()

    extracted_spectra = [
        spectrum
        for spectrum in spectra
        if check_ms1_or_low_energy_spectrum(spectrum, lower_collision_energy)
    ]

    output_exp = MSExperiment()
    output_exp.setSpectra(extracted_spectra)
    return output_exp


def extract_and_store_lower_energy_windows(mzml_path: Path, config_path: Path) -> None:
    """Given an mzML file with data consisting of higher- and lower-energy scan windows,
    create an mzML file that contains only the lower-energy and MS1 windows.
    Implementation loosely based on
    https://pyopenms.readthedocs.io/en/release_2.3.0/data_manipulation.html
    """
    output_path = mzml_path.parent / (mzml_path.stem + OUTPUT_SUFFIX_LOWER_ENERGY)
    config = Config.from_path(config_path)

    exp = MSExperiment()
    MzMLFile().load(str(mzml_path), exp)

    output_exp = extract_lower_energy_windows(exp, config.lower_collision_energy)

    output_file = get_diann_compatible_mzml_output_file()
    output_file.store(str(output_path), output_exp)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "mzml_path",
        type=str,
        help="Path to mzML file with higher and lower collision energy windows",
    )
    parser.add_argument(
        "--config_path",
        type=str,
        help="Path to config JSON file that defines values for higher and lower collision energy",
    )
    args = parser.parse_args()

    extract_and_store_lower_energy_windows(Path(args.mzml_path), Path(args.config_path))
