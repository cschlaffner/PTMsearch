from argparse import ArgumentParser
from pathlib import Path

from pyopenms import MSSpectrum, OnDiscMSExperiment, PlainMSDataWritingConsumer
from src.mzml_processing.constants import LOWER_COLLISION_ENERGY

OUTPUT_SUFFIX_LOWER_ENERGY = "_lower_energy.mzML"


def check_ms1_or_low_energy_spectrum(spectrum: MSSpectrum) -> bool:
    return (
        spectrum.getMSLevel() == 1
        or spectrum.getPrecursors()[0].getMetaValue("collision energy")
        == LOWER_COLLISION_ENERGY
    )


def count_lower_energy_spectra(exp: OnDiscMSExperiment) -> int:
    number_spectra = 0

    for i in range(exp.getNrSpectra()):
        spectrum = exp.getSpectrum(i)
        if check_ms1_or_low_energy_spectrum(spectrum):
            number_spectra += 1

    return number_spectra


def extract_lower_energy_windows(mzml_path: Path) -> None:
    """Given an mzML file with data consisting of higher- and lower-energy scan windows,
    create an mzML file that contains only the lower-energy windows.
    Implementation based on
    https://pyopenms.readthedocs.io/en/latest/user_guide/memory_management.html
    """
    output_path = mzml_path.parent / (mzml_path.stem + OUTPUT_SUFFIX_LOWER_ENERGY)

    exp = OnDiscMSExperiment()
    exp.openFile(str(mzml_path))
    consumer = PlainMSDataWritingConsumer(str(output_path))

    number_expected_spectra = count_lower_energy_spectra(exp)
    consumer.setExpectedSize(number_expected_spectra, 0)

    # TODO: a try run with subset and whole number of spectra, keep experimental settings? Chromatograms?
    for i in range(exp.getNrSpectra()):
        spectrum = exp.getSpectrum(i)
        if check_ms1_or_low_energy_spectrum(spectrum):
            consumer.consumeSpectrum(spectrum)

    del consumer


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "mzml_path",
        type=str,
        help="Path to mzML file with higher and lower collision energy windows",
    )
    args = parser.parse_args()

    extract_lower_energy_windows(Path(args.mzml_path))
