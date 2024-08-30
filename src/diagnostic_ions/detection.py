from enum import Enum
from pathlib import Path
from typing import List, Tuple, Union, cast

import numpy as np
import pandas as pd
from pyopenms import MSSpectrum

from src.diagnostic_ions.utils import get_modification_unimod_format
from src.mzml_processing.utils import validate_collision_energy_ms2_spectrum

MassToleranceUnit = Enum("MassToleranceUnit", ["ppm", "Da"])


class DiagnosticIonDetector:
    """Extracts diagnostic ion information from higher-energy MS2 spectra."""

    def __init__(
        self,
        known_ions_file: Path,
        mass_tolerance: float,
        mass_tolerance_unit: str,
        snr_threshold: float,
        higher_collision_energy: Union[float, int],
    ) -> None:
        self.mass_tolerance = mass_tolerance
        self.mass_tolerance_unit = mass_tolerance_unit
        self.known_ions = pd.read_csv(known_ions_file, header=0, dtype={"mz": float})
        self.snr_threshold = snr_threshold
        self.higher_collision_energy = higher_collision_energy

        self._validate_config()

        # Add single-letter amino acid and UniMod accession format
        # for compatibility with mod format in spectral libraries
        letter_and_unimod_format_mod = self.known_ions[
            ["amino_acid", "mod_name"]
        ].apply(
            lambda mod: get_modification_unimod_format(mod.amino_acid, mod.mod_name),
            axis=1,
        )
        self.known_ions.insert(
            0, "letter_and_unimod_format_mod", letter_and_unimod_format_mod
        )

    def _validate_config(self):
        assert np.array_equal(
            self.known_ions.columns, ["amino_acid", "mod_name", "type", "mz"]
        ), (
            "Known ions CSV file has wrong format, should include columns: "
            "amino_acid (name of the amino acid that is modified), "
            "mod_name (name of the modification), type (immonium ion, "
            "neutral loss or others)) and mz (m/z value of the diagnostic ion).",
        )
        allowed_units = [unit.name for unit in MassToleranceUnit]
        assert (
            self.mass_tolerance_unit in allowed_units
        ), f"Mass tolerance unit is {self.mass_tolerance_unit}, allowed are {allowed_units}."

    def _remove_noise(
        self, spectrum_mz: np.ndarray, intensities: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """check if the spectrum should be analysed at all with regard to the
        signal to noise ratio (maximum/mean) threshold and select all peaks larger
        than the mean intensity."""
        if len(intensities) == 0:
            return spectrum_mz, intensities

        mean_intensity = np.mean(intensities)
        max_peak_intensity = np.max(intensities)
        signal_to_noise_ratio = max_peak_intensity / mean_intensity

        if signal_to_noise_ratio < self.snr_threshold:
            return np.array([]), np.array([])

        thresholding_mask = intensities >= mean_intensity
        intensities_filtered = intensities[thresholding_mask]
        mz_filtered = spectrum_mz[thresholding_mask]
        return mz_filtered, intensities_filtered

    def _get_absolute_mass_tolerance(self, ion_mz: float) -> float:
        """Compute the actual mass tolerance from the specified tolerance and unit."""
        return (
            self.mass_tolerance * ion_mz / 1e6
            if self.mass_tolerance_unit == MassToleranceUnit.ppm.name
            else self.mass_tolerance
        )

    def extract_diagnostic_ions_for_spectrum(
        self, spectrum: MSSpectrum
    ) -> pd.DataFrame:
        """Extract information about the modifications
        for which a diagnostic ion matches a spectrum peak considering the tolerance."""
        validate_collision_energy_ms2_spectrum(spectrum, self.higher_collision_energy)

        spectrum_mz, intensities = spectrum.get_peaks()

        spectrum_mz, intensities = self._remove_noise(spectrum_mz, intensities)

        empty_df = pd.DataFrame(
            columns=[
                "spectrum_id",
                "amino_acid",
                "mod_name",
                "letter_and_unimod_format_mod",
                "type",
                "theoretical_mz",
                "detected_mz",
                "detected_intensity",
            ]
        )

        if len(spectrum_mz) == 0:
            return empty_df

        detected_ions_dfs = []

        for known_ion in self.known_ions.itertuples():
            # Casting to please the typechecking
            known_ion_mz = cast(float, known_ion.mz)
            tolerance = self._get_absolute_mass_tolerance(known_ion_mz)
            lower_border_mz = known_ion_mz - tolerance
            higher_border_mz = known_ion_mz + tolerance

            lower_border = np.searchsorted(spectrum_mz, lower_border_mz, side="left")
            higher_border = np.searchsorted(spectrum_mz, higher_border_mz, side="right")

            detected_peaks_mz = spectrum_mz[lower_border:higher_border]
            detected_peaks_intensities = intensities[lower_border:higher_border]
            if len(detected_peaks_mz) == 0:
                continue

            max_peak_idx = np.argmax(detected_peaks_intensities)
            max_peak_mz = detected_peaks_mz[max_peak_idx]
            max_peak_intensity = detected_peaks_intensities[max_peak_idx]
            detected_ions_dfs.append(
                pd.DataFrame(
                    {
                        "spectrum_id": [spectrum.getNativeID()],
                        "amino_acid": [known_ion.amino_acid],
                        "mod_name": [known_ion.mod_name],
                        "letter_and_unimod_format_mod": [
                            known_ion.letter_and_unimod_format_mod
                        ],
                        "type": [known_ion.type],
                        "theoretical_mz": [known_ion_mz],
                        "detected_mz": [max_peak_mz],
                        "detected_intensity": [max_peak_intensity],
                    }
                ),
            )

        if len(detected_ions_dfs) == 0:
            return empty_df

        detected_ions_df = pd.concat(
            detected_ions_dfs,
            ignore_index=True,
        )

        return detected_ions_df

    def validate_diagnostic_ion_spectra(self, spectra: List[MSSpectrum]) -> None:
        """Validates that spectra are suitable for diagnostic ion extraction
        (MS2 and higher collision energy)."""
        for spectrum in spectra:
            validate_collision_energy_ms2_spectrum(
                spectrum, self.higher_collision_energy
            )

    def extract_diagnostic_ions_for_spectra(
        self, spectra: List[MSSpectrum]
    ) -> pd.DataFrame:
        """Extract modification names for the provided spectra."""
        assert len(spectra) > 0, "Empty list of spectra was given."
        results = [
            self.extract_diagnostic_ions_for_spectrum(spectrum) for spectrum in spectra
        ]
        result = pd.concat(
            results,
            ignore_index=True,
        )

        return result

    # TODO: update docs to the df using version
