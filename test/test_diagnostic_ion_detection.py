from test.spectra import (
    COLLISION_ENERGY_HIGHER,
    spectrum_ms1,
    spectrum_ms2_higher_energy,
    spectrum_ms2_lower_energy,
)
from typing import List

import numpy as np
import pandas as pd
import pytest
from _pytest.fixtures import SubRequest
from pyopenms import MSSpectrum

from src.diagnostic_ions.detection import DiagnosticIonDetector

# ------------------------ Setup: values ---------------------------------------

intensity_threshold = 1000
known_ions_file = "test/mock_ions.csv"
mass_tolerance = 5
unit_ppm = "ppm"

spectrum_native_id = "s_1"
spectrum_2_native_id = "s_2"
peak = 2000.0
noise = 500.0
mz_mod_1 = 100.0
mz_mod_2 = 200.0
mz_mod_3 = 300.0

# not exact borders due to floating point errors, TODO: fix?
mz_mod_1_threshold_border_higher = 100.000499999
mz_mod_1_threshold_border_lower = 99.99950000001


# ------------------------ Setup: m/z and intensities ---------------------------------------


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

mz_intensities_below_threshold_peak = (
    np.array([mz_mod_1]),
    np.array([noise]),
)

mz_intensities_within_ppm_tolerance = (
    np.array([mz_mod_1_threshold_border_higher]),
    np.array([peak]),
)

mz_intensities_within_ppm_tolerance_multiple = (
    np.array([mz_mod_1_threshold_border_lower, mz_mod_1_threshold_border_higher]),
    np.array([peak, peak]),
)

mz_intensities_only_mod1 = (
    np.array([mz_mod_1]),
    np.array([peak]),
)

mz_intensities_only_mod12 = (
    np.array([mz_mod_1, mz_mod_2]),
    np.array([peak, peak]),
)


# ------------------------ Setup: correct result DataFrames ---------------------------------------


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

below_threshold_ion_df = pd.DataFrame(
    {
        "spectrum_id": [spectrum_native_id],
        "amino_acid": ["a_1"],
        "mod_name": ["m_1"],
        "theoretical_mz": [mz_mod_1],
        "detected_mz": [mz_mod_1],
        "detected_intensity": [noise],
    }
)

within_ppm_tolerance_df = pd.DataFrame(
    {
        "spectrum_id": [spectrum_native_id],
        "amino_acid": ["a_1"],
        "mod_name": ["m_1"],
        "theoretical_mz": [mz_mod_1],
        "detected_mz": [mz_mod_1_threshold_border_higher],
        "detected_intensity": [peak],
    }
)

within_ppm_tolerance_multiple_df = pd.DataFrame(
    {
        "spectrum_id": [spectrum_native_id, spectrum_native_id],
        "amino_acid": ["a_1", "a_1"],
        "mod_name": ["m_1", "m_1"],
        "theoretical_mz": [mz_mod_1, mz_mod_1],
        "detected_mz": [
            mz_mod_1_threshold_border_lower,
            mz_mod_1_threshold_border_higher,
        ],
        "detected_intensity": [peak, peak],
    }
)

multiple_spectra_df = pd.DataFrame(
    {
        "spectrum_id": [
            spectrum_native_id,
            spectrum_native_id,
            spectrum_2_native_id,
        ],
        "amino_acid": ["a_1", "a_1", "a_1"],
        "mod_name": ["m_1", "m_2", "m_1"],
        "theoretical_mz": [mz_mod_1, mz_mod_2, mz_mod_1],
        "detected_mz": [mz_mod_1, mz_mod_2, mz_mod_1],
        "detected_intensity": [peak, peak, peak],
    }
)


# ------------------------ Setup: diagnostic ion detectors ---------------------------------------


@pytest.fixture
def detector_exact_matching() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        0,
        unit_ppm,
        intensity_threshold,
        COLLISION_ENERGY_HIGHER,
    )


@pytest.fixture
def detector_with_ppm_tolerance() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        mass_tolerance,
        unit_ppm,
        intensity_threshold,
        COLLISION_ENERGY_HIGHER,
    )


@pytest.fixture
def detector_exact_matching_no_intensity_threshold() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        0,
        unit_ppm,
        0,
        COLLISION_ENERGY_HIGHER,
    )


def assert_detection_results_correct(
    detector: DiagnosticIonDetector, spectrum: MSSpectrum, expected_df: pd.DataFrame
) -> None:
    detected_ions_df = detector.extract_diagnostic_ions_for_spectrum(spectrum)
    pd.testing.assert_frame_equal(
        detected_ions_df,
        expected_df,
        # not comparing the dtype to avoid failing because of float64 VS float32
        check_dtype=False,
    )


# ------------------------ Tests: ion detection single spectrum ---------------------------------------


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
def higher_energy_spectrum_exact_matching(request: SubRequest) -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id, peaks_mz_intensities=request.param
    )


def test_detect_ions_spectrum_exact_matching_applied_intensity_threshold(
    higher_energy_spectrum_exact_matching: MSSpectrum,
    detector_exact_matching: DiagnosticIonDetector,
) -> None:
    """The diagnostic ions should be detected in the spectrum, additional noise
    or peaks at irrelevant positions should be ignored."""
    assert_detection_results_correct(
        detector_exact_matching,
        higher_energy_spectrum_exact_matching,
        all_ions_from_spectrum_exact_matching_df,
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
def higher_energy_spectrum_no_ions(request: SubRequest) -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id, peaks_mz_intensities=request.param
    )


def test_no_ions_in_spectrum_applied_intensity_threshold(
    higher_energy_spectrum_no_ions: MSSpectrum,
    detector_exact_matching: DiagnosticIonDetector,
) -> None:
    """If there are no ions in the spectra, no ions should be detected. Further,
    empty m/z and intensity arrays should not lead to errors."""
    assert_detection_results_correct(
        detector_exact_matching, higher_energy_spectrum_no_ions, empty_results_df
    )


@pytest.fixture
def higher_energy_spectrum_below_threshold() -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id,
        peaks_mz_intensities=mz_intensities_below_threshold_peak,
    )


def test_lower_intensity_peak_detected_without_intensity_threshold(
    higher_energy_spectrum_below_threshold: MSSpectrum,
    detector_exact_matching_no_intensity_threshold: DiagnosticIonDetector,
) -> None:
    """When no intensity threshold is applied, peaks with a lower intensity at the
    ion positions are detected."""
    assert_detection_results_correct(
        detector_exact_matching_no_intensity_threshold,
        higher_energy_spectrum_below_threshold,
        below_threshold_ion_df,
    )


@pytest.fixture
def higher_energy_spectrum_within_ppm_tolerance() -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id,
        peaks_mz_intensities=mz_intensities_within_ppm_tolerance,
    )


def test_within_ppm_tolerance_detected(
    higher_energy_spectrum_within_ppm_tolerance: MSSpectrum,
    detector_with_ppm_tolerance: DiagnosticIonDetector,
) -> None:
    """When a ppm tolerance is set, peaks within the tolerance range
    should be detected."""
    assert_detection_results_correct(
        detector_with_ppm_tolerance,
        higher_energy_spectrum_within_ppm_tolerance,
        within_ppm_tolerance_df,
    )


# TODO: add test for Da tolerance


@pytest.fixture
def higher_energy_spectrum_within_ppm_tolerance_multiple() -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id,
        peaks_mz_intensities=mz_intensities_within_ppm_tolerance_multiple,
    )


def test_within_ppm_tolerance_multiple_detected(
    higher_energy_spectrum_within_ppm_tolerance_multiple: MSSpectrum,
    detector_with_ppm_tolerance: DiagnosticIonDetector,
) -> None:
    """When a ppm tolerance is set, and multiple peaks are within the tolerance
    of a known ion m/z, all those peaks should be detected."""
    assert_detection_results_correct(
        detector_with_ppm_tolerance,
        higher_energy_spectrum_within_ppm_tolerance_multiple,
        within_ppm_tolerance_multiple_df,
    )


# ------------------------ Tests: error cases ---------------------------------------


@pytest.fixture
def ms1_spectrum() -> MSSpectrum:
    return spectrum_ms1(native_id=spectrum_native_id)


def test_ion_detection_ms1_fails(
    ms1_spectrum: MSSpectrum,
    detector_exact_matching: DiagnosticIonDetector,
) -> None:
    """Diagnostic ion detection for an MS1 scan should fail."""
    with pytest.raises(AssertionError, match="MS2"):
        detector_exact_matching.extract_diagnostic_ions_for_spectrum(ms1_spectrum)


@pytest.fixture
def lower_energy_spectrum() -> MSSpectrum:
    return spectrum_ms2_lower_energy(native_id=spectrum_native_id)


def test_ion_detection_lower_energy_fails(
    lower_energy_spectrum: MSSpectrum,
    detector_exact_matching: DiagnosticIonDetector,
) -> None:
    """Diagnostic ion detection for a lower-energy MS2 scan should fail."""
    with pytest.raises(AssertionError, match="higher-energy"):
        detector_exact_matching.extract_diagnostic_ions_for_spectrum(
            lower_energy_spectrum
        )


# ------------------------ Tests: multiple spectra ---------------------------------------


@pytest.fixture
def higher_energy_multiple_spectra_exact_matching() -> List[MSSpectrum]:
    return [
        spectrum_ms2_higher_energy(
            native_id=spectrum_native_id,
            peaks_mz_intensities=mz_intensities_only_mod12,
        ),
        spectrum_ms2_higher_energy(
            native_id=spectrum_2_native_id,
            peaks_mz_intensities=mz_intensities_only_mod1,
        ),
    ]


def test_detect_ions_multiple_spectra_exact_matching(
    higher_energy_multiple_spectra_exact_matching: List[MSSpectrum],
    detector_exact_matching: DiagnosticIonDetector,
) -> None:
    """The diagnostic ions should be detected in the spectra and stored with the corresponding
    spectrum id."""
    detected_ions_df = detector_exact_matching.extract_diagnostic_ions_for_spectra(
        higher_energy_multiple_spectra_exact_matching
    )
    pd.testing.assert_frame_equal(
        detected_ions_df,
        multiple_spectra_df,
        # not comparing the dtype to avoid failing because of float64 VS float32
        check_dtype=False,
    )


def test_validate_spectra_successful(
    higher_energy_multiple_spectra_exact_matching: List[MSSpectrum],
    detector_exact_matching: DiagnosticIonDetector,
) -> None:
    """Validation of a list of higher-energy MS2 spectra should succeed."""
    detector_exact_matching.validate_diagnostic_ion_spectra(
        higher_energy_multiple_spectra_exact_matching
    )


def test_validate_MS1_spectra_fails(
    ms1_spectrum: MSSpectrum,
    detector_exact_matching: DiagnosticIonDetector,
) -> None:
    """Validation of a list of MS1 spectra should fail."""
    with pytest.raises(AssertionError, match="MS2"):
        detector_exact_matching.validate_diagnostic_ion_spectra([ms1_spectrum])


def test_validate_lower_energy_spectra_fails(
    lower_energy_spectrum: MSSpectrum,
    detector_exact_matching: DiagnosticIonDetector,
) -> None:
    """Validation of a list of lower-energy spectra should fail."""
    with pytest.raises(AssertionError, match="higher-energy"):
        detector_exact_matching.validate_diagnostic_ion_spectra([lower_energy_spectrum])
