from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from pyopenms import MSSpectrum

from src.mzml_processing.utils import get_spectrum_collision_energy


class DiagnosticIonDetector:
    def __init__(
        self,
        known_ions_file: Path,
        mass_tolerance: float,
    ) -> None:
        self.mass_tolerance = mass_tolerance
        self.known_ions = pd.read_csv(
            known_ions_file, header=0, index_col=0, dtype={"mz": float}
        )
        assert self.known_ions.index.name == "mod_name" and self.known_ions.columns == [
            "mz"
        ], "Known ions CSV file has wrong format, should include columns mod_name (name of the modification), and mz (m/z value of the diagnostic ion).",
    

    def extract_diagnostic_ions_for_spectrum(self, spectrum: MSSpectrum) -> List[str]:
        """Extract the names of the modifications for which a diagnostic ion
        matches a spectrum peak (TODO: noise handling) considering the tolerance."""

        spectrum_mz, intensities = spectrum.get_peaks()
        # WIP, TODO: add signal-to-noise ratio
        detected_ions_mask = [
            np.any(
                [
                    np.isclose(known_ion_mz, peak_mz, rtol=0, atol=self.mass_tolerance)
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
        """Extract modification names for the provided spectra including validation (MS level 2 and higher collision energy"""
        self.validate_diagnostic_ion_spectra(spectra, higher_collision_energy)
        return self.extract_diagnostic_ions_for_spectra(spectra)
