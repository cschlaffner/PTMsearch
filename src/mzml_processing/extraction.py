import re
from argparse import ArgumentParser
from pathlib import Path
from typing import Callable, Tuple, Union

import pandas as pd
from pyopenms import MSExperiment, MSSpectrum, MzMLFile

from src.config.config import Config
from src.mzml_processing.utils import (
    check_collision_energy_ms2_spectrum,
    check_ms1_spectrum,
    get_diann_compatible_mzml_output_file,
)

OUTPUT_SUFFIX_LOWER_ENERGY = "_lower_energy.mzML"


class ScanWindowExtractor:
    def __init__(
        self,
        lower_collision_energy: Union[float, int],
        higher_collision_energy: Union[float, int],
    ) -> None:
        self.lower_collision_energy = lower_collision_energy
        self.higher_collision_energy = higher_collision_energy

    def _check_ms1_or_lower_energy_spectrum(self, spectrum: MSSpectrum) -> bool:
        return check_ms1_spectrum(spectrum) or check_collision_energy_ms2_spectrum(
            spectrum, self.lower_collision_energy
        )

    def _check_ms1_or_higher_energy_spectrum(self, spectrum: MSSpectrum) -> bool:
        return check_ms1_spectrum(spectrum) or check_collision_energy_ms2_spectrum(
            spectrum, self.higher_collision_energy
        )

    def _check_higher_energy_spectrum(self, spectrum: MSSpectrum) -> bool:
        return check_collision_energy_ms2_spectrum(
            spectrum, self.higher_collision_energy
        )

    def _extract_windows_for_criterion(
        self,
        exp: MSExperiment,
        filtering_criterion: Callable[[MSSpectrum], bool],
    ) -> MSExperiment:
        """Extracts all scan windows that fullfil the criterion to a second MSExperiment"""
        spectra = exp.getSpectra()

        extracted_spectra = [
            spectrum for spectrum in spectra if filtering_criterion(spectrum)
        ]

        output_exp = MSExperiment()
        output_exp.setSpectra(extracted_spectra)
        return output_exp

    def extract_ms1_and_lower_energy_windows(self, exp: MSExperiment) -> MSExperiment:
        """Extracts all lower-energy and MS1 scan windows to a second MSExperiment"""
        return self._extract_windows_for_criterion(
            exp, self._check_ms1_or_lower_energy_spectrum
        )

    def extract_ms1_and_higher_energy_windows(self, exp: MSExperiment) -> MSExperiment:
        """Extracts all higher-energy and MS1 scan windows to a second MSExperiment"""
        return self._extract_windows_for_criterion(
            exp, self._check_ms1_or_higher_energy_spectrum
        )

    def extract_higher_energy_windows(self, exp: MSExperiment) -> MSExperiment:
        """Extracts all higher-energy MS2 scan windows to a second MSExperiment"""
        return self._extract_windows_for_criterion(
            exp, self._check_higher_energy_spectrum
        )

    def extract_ms1_windows(self, exp: MSExperiment) -> MSExperiment:
        """Extracts all MS1 scan windows to a second MSExperiment"""
        return self._extract_windows_for_criterion(exp, check_ms1_spectrum)

    def rename_spectrum_ids(
        self, exp: MSExperiment, return_id_mapping: bool
    ) -> Union[MSExperiment, Tuple[MSExperiment, pd.DataFrame]]:
        # TODO: add test
        spectra = exp.getSpectra()
        original_ids = []
        renamed_ids = []

        scan_id_regex = re.compile("scan=[0-9]+")

        for i, spectrum in enumerate(spectra):
            original_id = spectrum.getNativeID()
            original_title = spectrum.getMetaValue("spectrum title")
            original_ids.append(original_id)
            original_id_number = re.findall(scan_id_regex, original_id)[0][5:]
            new_id_number = i + 1

            new_id = re.sub(scan_id_regex, f"scan={new_id_number}", original_id)
            new_title_part = re.sub(
                scan_id_regex, f"scan={new_id_number}", original_title
            )
            new_title = re.sub(
                f"{original_id_number}.{original_id_number}",
                f"{new_id_number}.{new_id_number}",
                new_title_part,
            )

            spectrum.setMetaValue("spectrum title", new_title)
            spectrum.setNativeID(new_id)
            renamed_ids.append(new_id)

        exp.setSpectra(spectra)
        if not return_id_mapping:
            return exp

        id_mapping_df = pd.DataFrame(
            {"original_id": original_ids, "renamed_id": renamed_ids}
        )
        return exp, id_mapping_df


def extract_and_store_ms1_and_lower_energy_windows(
    mzml_path: Path, config_path: Path
) -> None:
    """Given an mzML file with data consisting of higher- and lower-energy scan windows,
    create an mzML file that contains only the lower-energy and MS1 windows.
    Implementation loosely based on
    https://pyopenms.readthedocs.io/en/release_2.3.0/data_manipulation.html
    """
    output_path = mzml_path.parent / (mzml_path.stem + OUTPUT_SUFFIX_LOWER_ENERGY)
    config = Config.from_path(config_path)

    exp = MSExperiment()
    MzMLFile().load(str(mzml_path), exp)

    extractor = ScanWindowExtractor(
        config.lower_collision_energy, config.higher_collision_energy
    )

    output_exp = extractor.extract_ms1_and_lower_energy_windows(exp)
    output_exp, mapping = extractor.rename_spectrum_ids(output_exp, True)

    output_file = get_diann_compatible_mzml_output_file()
    output_file.store(str(output_path), output_exp)
    mapping.to_csv(mzml_path.parent / (mzml_path.stem + "_id_mapping.csv"), index=False)


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

    extract_and_store_ms1_and_lower_energy_windows(
        Path(args.mzml_path), Path(args.config_path)
    )
