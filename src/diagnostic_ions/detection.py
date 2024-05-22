import datetime
from enum import Enum
from pathlib import Path
from time import time
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
from pyopenms import MSSpectrum
from src.mzml_processing.utils import check_collision_energy_ms2_spectrum

MassToleranceUnit = Enum("MassToleranceUnit", ["ppm", "Da"])


class DiagnosticIonDetector:
    """Extracts diagnostic ion information from higher-energy MS2 spectra."""

    def __init__(
        self,
        known_ions_file: Path,
        mass_tolerance: float,
        mass_tolerance_unit: str,
        intensity_threshold: float,
        higher_collision_energy: Union[float, int],
    ) -> None:
        self.mass_tolerance = mass_tolerance
        self.mass_tolerance_unit = mass_tolerance_unit
        self.known_ions = pd.read_csv(known_ions_file, header=0, dtype={"mz": float})
        self.intensity_threshold = intensity_threshold
        self.higher_collision_energy = higher_collision_energy

        self._validate_config()

        # TODO: delete
        self.search_times = []
        self.concat_times = []
        self.final_concat_time = []

    def _validate_config(self):
        assert np.array_equal(
            self.known_ions.columns, ["amino_acid", "mod_name", "mz"]
        ), (
            "Known ions CSV file has wrong format, should include columns: amino_acid",
            "(name of the amino acid that is modified), mod_name (name of the modification),",
            "and mz (m/z value of the diagnostic ion).",
        )
        assert self.mass_tolerance_unit in [unit.name for unit in MassToleranceUnit]

    # TODO: clean up
    def _remove_noise(
        self, spectrum_mz: np.ndarray, intensities: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        # """Select peaks that pass the signal-to-noise ratio
        # (signal / median absolute deviation) threshold.
        # TODO: add source for MAD usage?"""
        # # TODO: redo
        # signal_to_noise_ratio = intensities / median_abs_deviation(intensities)
        """Select peaks that pass the intensity threshold."""
        thresholding_mask = intensities >= self.intensity_threshold
        intensities_filtered = intensities[thresholding_mask]
        mz_filtered = spectrum_mz[thresholding_mask]
        return mz_filtered, intensities_filtered

    def _get_absolute_mass_tolerance(self, ion_mz: float) -> float:
        """Compute the actual mass tolerance from the specified tolerance and unit."""
        # TODO: check Dalton
        return (
            self.mass_tolerance * ion_mz / 1e6
            if self.mass_tolerance_unit == MassToleranceUnit.ppm.name
            else self.mass_tolerance
        )

    def _validate_spectrum(self, spectrum: MSSpectrum) -> None:
        """Validates that a spectrum is of MS level 2 and with higher collision energy in order
        to extract diagnostic ions from it."""
        assert check_collision_energy_ms2_spectrum(
            spectrum, self.higher_collision_energy
        ), "Spectrum for diagnostic ion extraction must be a higher-energy MS2 scan."

    def extract_diagnostic_ions_for_spectrum(
        self, spectrum: MSSpectrum
    ) -> pd.DataFrame:
        """Extract information about the modifications
        for which a diagnostic ion matches a spectrum peak considering the tolerance."""
        self._validate_spectrum(spectrum)

        spectrum_mz, intensities = spectrum.get_peaks()

        spectrum_mz, intensities = self._remove_noise(spectrum_mz, intensities)

        detected_ions_dfs = [
            pd.DataFrame(
                columns=[
                    "spectrum_id",
                    "amino_acid",
                    "mod_name",
                    "theoretical_mz",
                    "detected_mz",
                    "detected_intensity",
                ]
            )
        ]
        a = time()
        for known_ion in self.known_ions.itertuples():
            tolerance = self._get_absolute_mass_tolerance(known_ion.mz)
            lower_border_mz = known_ion.mz - tolerance
            higher_border_mz = known_ion.mz + tolerance

            lower_border = np.searchsorted(spectrum_mz, lower_border_mz, side="left")
            higher_border = np.searchsorted(spectrum_mz, higher_border_mz, side="right")

            detected_peaks_mz = spectrum_mz[lower_border:higher_border]
            detected_peaks_intensities = intensities[lower_border:higher_border]
            num_hits = len(detected_peaks_mz)

            detected_ions_dfs.append(
                pd.DataFrame(
                    {
                        "spectrum_id": np.repeat(spectrum.getNativeID(), num_hits),
                        "amino_acid": np.repeat(known_ion.amino_acid, num_hits),
                        "mod_name": np.repeat(known_ion.mod_name, num_hits),
                        "theoretical_mz": np.repeat(known_ion.mz, num_hits),
                        "detected_mz": detected_peaks_mz,
                        "detected_intensity": detected_peaks_intensities,
                    }
                ),
            )

        b = time()
        self.search_times.append(b - a)

        c = time()
        detected_ions_df = pd.concat(
            detected_ions_dfs,
            ignore_index=True,
        )

        d = time()
        self.concat_times.append(d - c)

        return detected_ions_df

    def validate_diagnostic_ion_spectra(self, spectra: List[MSSpectrum]) -> None:
        """Validates that spectra are suitable for diagnostic ion extraction
        (MS2 and higher collision energy)."""
        for spectrum in spectra:
            self._validate_spectrum(spectrum)

    def extract_diagnostic_ions_for_spectra(
        self, spectra: List[MSSpectrum]
    ) -> pd.DataFrame:
        """Extract modification names for the provided spectra."""
        results = []
        a = time()
        for i, spectrum in enumerate(spectra):
            results.append(self.extract_diagnostic_ions_for_spectrum(spectrum))
            if i % 1000 == 0:
                print(f"{datetime.datetime.now()}: {i} spectra analyzed", flush=True)
        b = time()

        print(f"average search time: {np.mean(self.search_times)}", flush=True)
        print(f"average concat time: {np.mean(self.concat_times)}", flush=True)
        c = time()
        result = pd.concat(
            results,
            ignore_index=True,
        )
        d = time()

        print(f"total search time: {b - a}", flush=True)
        print(f"final concat time: {d - c}", flush=True)
        return result

    # TODO: update docs to the df using version, make sure to do duplicate handling afterwards!
