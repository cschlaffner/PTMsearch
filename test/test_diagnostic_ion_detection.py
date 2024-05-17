from test.spectra import (
    COLLISION_ENERGY_HIGHER,
    spectrum_ms1,
    spectrum_ms2_higher_energy,
    spectrum_ms2_lower_energy,
)

import numpy as np
import pandas as pd
import pytest
from pyopenms import MSSpectrum

from src.diagnostic_ions.detection import DiagnosticIonDetector

intensity_threshold = 1000
known_ions_file = "test/mock_ions.csv"
mass_tolerance = 5
unit_ppm = "ppm"

spectrum_native_id = "s_1"
peak = 2000.0
noise = 500.0
mz_mod_1 = 100.0
mz_mod_2 = 200.0
mz_mod_3 = 300.0

mz_intensities_empty = (np.array([]), np.array([]))
mz_intensities_only_ions = (
    np.array([mz_mod_1, mz_mod_2, mz_mod_3]),
    np.array([peak, peak, peak]),
)
mz_intensities_noise_at_ion_positions = (
    np.array([mz_mod_1, mz_mod_2, mz_mod_3]),
    np.array([noise, noise, noise]),
)
mz_intensities_irrelevant_peaks = (
    np.array([400.0, 500.0, 600.0]),
    np.array([peak, peak, peak]),
)
mz_intensities_ions_additional_noise = (
    np.array([mz_mod_1, 150.0, 170.0, mz_mod_2, 220.0, mz_mod_3, 400.0]),
    np.array([peak, noise, noise, peak, noise, peak, noise]),
)

mz_intensities_ions_additional_noise_irrelevant_peaks = (
    np.array([mz_mod_1, 150.0, 170.0, mz_mod_2, 220.0, mz_mod_3, 400.0]),
    np.array([peak, peak, noise, peak, noise, peak, peak]),
)


empty_results_df = pd.DataFrame(
    columns=[
        "spectrum_id",
        "amino_acid",
        "mod_name",
        "theoretical_mz",
        "detected_mz",
        "detected_intensity",
    ]
)

all_ions_from_spectrum_exact_matching_df = pd.DataFrame(
    {
        "spectrum_id": [spectrum_native_id, spectrum_native_id, spectrum_native_id],
        "amino_acid": ["a_1", "a_1", "a_2"],
        "mod_name": ["m_1", "m_2", "m_3"],
        "theoretical_mz": [mz_mod_1, mz_mod_2, mz_mod_3],
        "detected_mz": [mz_mod_1, mz_mod_2, mz_mod_3],
        "detected_intensity": [peak, peak, peak],
    }
)


@pytest.fixture(
    params=[
        mz_intensities_only_ions,
        mz_intensities_ions_additional_noise,
        mz_intensities_ions_additional_noise_irrelevant_peaks,
    ],
    ids=[
        "only_ions",
        "ions_additional_noise",
        "ions_additional_noise_irrelevant_peaks",
    ],
)
def higher_energy_spectrum_exact_matching(request) -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id, peaks_mz_intensities=request.param
    )


def test_detect_ions_spectrum_exact_matching_applied_intensity_threshold(
    higher_energy_spectrum_exact_matching,
) -> None:
    """The diagnostic ions should be detected in the spectrum, additional noise
    or peaks at irrelevant positions should be ignored."""
    detector = DiagnosticIonDetector(
        known_ions_file,
        mass_tolerance,
        unit_ppm,
        intensity_threshold,
        COLLISION_ENERGY_HIGHER,
    )
    detected_ions_df = detector.extract_diagnostic_ions_for_spectrum(
        higher_energy_spectrum_exact_matching
    )
    pd.testing.assert_frame_equal(
        # not comparing the dtype to avoid failing because of float64 VS float32
        detected_ions_df,
        all_ions_from_spectrum_exact_matching_df,
        check_dtype=False,
    )


@pytest.fixture(
    params=[
        mz_intensities_empty,
        mz_intensities_noise_at_ion_positions,
        mz_intensities_irrelevant_peaks,
    ],
    ids=[
        "empty",
        "noise_at_ion_positions",
        "irrelevant_peaks",
    ],
)
def higher_energy_spectrum_no_ions(request) -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id, peaks_mz_intensities=request.param
    )


def test_no_ions_in_spectrum_applied_intensity_threshold(
    higher_energy_spectrum_no_ions,
) -> None:
    """If there are no ions in the spectra, no ions should be detected. Further,
    empty m/z and intensity arrays should not lead to errors."""
    detector = DiagnosticIonDetector(
        known_ions_file,
        mass_tolerance,
        unit_ppm,
        intensity_threshold,
        COLLISION_ENERGY_HIGHER,
    )
    detected_ions_df = detector.extract_diagnostic_ions_for_spectrum(
        higher_energy_spectrum_no_ions
    )
    pd.testing.assert_frame_equal(
        detected_ions_df,
        empty_results_df,
    )


@pytest.fixture
def higher_energy_spectrum_noise_at_ion_position() -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id,
        # TODO: only one
        peaks_mz_intensities=mz_intensities_noise_at_ion_positions,
    )


def test_noise_detected_without_intensity_threshold(
    higher_energy_spectrum_noise_at_ion_positions,
) -> None:
    """When no intensity threshold is applied, the noise at the ion positions
    is detected as ions."""
    pass


"""
plus:
    - test mass tolerance
        -ppm
        -dalton (TBD)
    - test multiple peaks within mass tolerance threshold
    - should fail for MS1 spectrum
    - should fail for lower-energy spectrum
    
    - test spectrum list"""
