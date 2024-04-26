from argparse import ArgumentParser
from pathlib import Path

from pyopenms import (
    MSExperiment,
    MSSpectrum,
    MzMLFile,
    OnDiscMSExperiment,
    PeakFileOptions,
    PlainMSDataWritingConsumer,
)
from src.config.config import Config

OUTPUT_SUFFIX_LOWER_ENERGY = "_lower_energy.mzML"


def check_ms1_or_low_energy_spectrum(
    spectrum: MSSpectrum, lower_collision_energy: int
) -> bool:
    return (
        spectrum.getMSLevel() == 1
        or spectrum.getPrecursors()[0].getMetaValue("collision energy")
        == lower_collision_energy
    )


def extract_lower_energy_windows(mzml_path: Path, config_path: Path) -> None:
    """Given an mzML file with data consisting of higher- and lower-energy scan windows,
    create an mzML file that contains only the lower-energy windows.
    Implementation based on
    https://pyopenms.readthedocs.io/en/latest/user_guide/memory_management.html
    """
    output_path = mzml_path.parent / (mzml_path.stem + OUTPUT_SUFFIX_LOWER_ENERGY)
    config = Config.from_path(config_path)

    exp = MSExperiment()
    MzMLFile().load(str(mzml_path), exp)

    spectra = exp.getSpectra()

    extracted_spectra = [
        spectrum
        for spectrum in spectra
        if check_ms1_or_low_energy_spectrum(spectrum, config.lower_collision_energy)
    ]

    output_exp = MSExperiment()
    output_exp.setSpectra(extracted_spectra)

    options = PeakFileOptions()
    options.setMz32Bit(True)
    options.setIntensity32Bit(True)

    output_file = MzMLFile()
    output_file.setOptions(options)
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

    extract_lower_energy_windows(Path(args.mzml_path), Path(args.config_path))
