from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from pyopenms import MSSpectrum
from scipy.stats import median_abs_deviation

from src.mzml_processing.utils import get_spectrum_collision_energy


class DiagnosticIonDetector:
    """Extracts diagnostic ion information from higher-energy MS2 spectra."""

    def __init__(
        self,
        known_ions_file: Path,
        mass_tolerance: float,
        has_noisy_spectra: bool,
        signal_to_noise_ratio_threshold: Optional[float] = None,
    ) -> None:
        self.mass_tolerance = mass_tolerance

        self.known_ions = pd.read_csv(
            known_ions_file, header=0, index_col=0, dtype={"mz": float}
        )
        assert self.known_ions.index.name == "mod_name" and self.known_ions.columns == [
            "mz"
        ], (
            "Known ions CSV file has wrong format, should include columns mod_name",
            "(name of the modification), and mz (m/z value of the diagnostic ion).",
        )

        self.has_noisy_spectra = has_noisy_spectra
        if has_noisy_spectra:
            assert signal_to_noise_ratio_threshold is not None, (
                "To remove noise from noisy spectra, a threshold for the",
                "signal-to-noise ratio (signal / median absolute deviation) is required.",
            )
            self.signal_to_noise_ratio_threshold = signal_to_noise_ratio_threshold

    def _remove_noise(
        self, spectrum_mz: np.ndarray, intensities: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Select peaks that pass the signal-to-noise ratio
        (signal / median absolute deviation) threshold.
        TODO: add source for MAD usage?"""
        signal_to_noise_ratio = intensities / median_abs_deviation(intensities)
        intensities_filtered = intensities[
            signal_to_noise_ratio >= self.signal_to_noise_ratio_threshold
        ]
        mz_filtered = spectrum_mz[
            signal_to_noise_ratio >= self.signal_to_noise_ratio_threshold
        ]
        return mz_filtered, intensities_filtered

    def _get_absolute_mass_tolerance(self, ion_mz: float) -> float:
        """Compute the actual mass tolerance from the specified tolerance and unit,
        currently only ppm."""
        # TODO: add Dalton
        return self.mass_tolerance * ion_mz / 10e6

    def extract_diagnostic_ions_for_spectrum(self, spectrum: MSSpectrum) -> List[str]:
        """Extract the names (TODO: add more information) of the modifications
        for which a diagnostic ion matches a spectrum peak considering the tolerance."""

        spectrum_mz, intensities = spectrum.get_peaks()

        if self.has_noisy_spectra:
            spectrum_mz, intensities = self._remove_noise(spectrum_mz, intensities)

        detected_ions_mask = [
            np.any(
                [  # TODO: add ppm/Dalton computation for mass tolerance
                    np.isclose(
                        known_ion_mz,
                        peak_mz,
                        rtol=0,
                        atol=self._get_absolute_mass_tolerance(known_ion_mz),
                    )
                    for peak_mz in spectrum_mz
                ]
            )
            for known_ion_mz in self.known_ions["mz"]
        ]

        return self.known_ions[detected_ions_mask].index.to_list()

    def validate_diagnostic_ion_spectra(
        self, spectra: List[MSSpectrum], higher_collision_energy: float
    ) -> None:
        """Validates that spectra are of MS level 2 and with higher collision energy in order
        to extract diagnostic ions from them."""
        for spectrum in spectra:
            assert (
                spectrum.getMSLevel() == 2
                and get_spectrum_collision_energy(spectrum) == higher_collision_energy
            )

    def extract_diagnostic_ions_for_spectra(
        self, spectra: List[MSSpectrum]
    ) -> List[List[str]]:
        """Extract modification names for the provided spectra."""
        return [
            self.extract_diagnostic_ions_for_spectrum(spectrum) for spectrum in spectra
        ]

    def extract_diagnostic_ions_for_spectra_with_validation(
        self, spectra: List[MSSpectrum], higher_collision_energy: float
    ) -> List[List[str]]:
        """Extract modification names for the provided spectra
        including validation (MS level 2 and higher collision energy)."""
        self.validate_diagnostic_ion_spectra(spectra, higher_collision_energy)
        return self.extract_diagnostic_ions_for_spectra(spectra)
