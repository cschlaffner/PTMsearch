from pathlib import Path
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

snr_threshold = 2
snr_threshold_higher = 3
known_ions_file = Path("test/mock_ions.csv")
mass_tolerance_ppm = 5
unit_ppm = "ppm"
mass_tolerance_da = 0.0005
unit_da = "Da"

spectrum_native_id = "s_1"
spectrum_2_native_id = "s_2"
peak = 10000.0
peak_higher = 1000000.0
noise = 500.0
mz_mod_1 = 100.0
mz_mod_2 = 200.0
mz_mod_3 = 300.0

# not exact borders due to floating point errors
mz_mod_1_threshold_border_higher = 100.000499999
mz_mod_1_threshold_border_lower = 99.99950000001


# ------------------------ Setup: m/z and intensities ---------------------------------------


mz_intensities_empty = (np.array([]), np.array([]))
mz_intensities_noise_at_ion_positions = (
    np.array([mz_mod_1, mz_mod_2, mz_mod_3]),
    np.array([noise, noise, noise]),
)
mz_intensities_ions_additional_noise = (
    np.array([mz_mod_1, 150.0, 170.0, mz_mod_2, 220.0, mz_mod_3, 400.0]),
    np.array([peak, noise, noise, peak, noise, peak, noise]),
)

mz_intensities_single_higher_peak_additional_noise = (
    np.array([mz_mod_1, 150.0, 170.0, mz_mod_2, 220.0, mz_mod_3, 400.0]),
    np.array([peak, noise, noise, peak_higher, noise, peak, noise]),
)

mz_intensities_ions_additional_noise_irrelevant_peaks = (
    np.array(
        [
            mz_mod_1,
            150.0,
            170.0,
            mz_mod_2,
            220.0,
            mz_mod_3,
            400.0,
            500.0,
            600.0,
            700.0,
            800.0,
        ]
    ),
    np.array([peak, peak, noise, peak, noise, peak, peak, noise, noise, noise, noise]),
)

mz_intensities_within_mass_tolerance_irrelevant_peaks = (
    np.array([mz_mod_1_threshold_border_higher, 150.0, 400.0]),
    np.array([peak, peak, peak]),
)

mz_intensities_within_mass_tolerance_multiple = (
    np.array([mz_mod_1_threshold_border_lower, mz_mod_1_threshold_border_higher]),
    np.array([peak, peak_higher]),
)

mz_intensities_within_mass_tolerance_additional_noise = (
    np.array([mz_mod_1_threshold_border_higher, 200.0, 300.0]),
    np.array([peak, noise, noise]),
)

mz_intensities_only_mod1_additional_noise = (
    np.array([mz_mod_1, 150.0, 20.0]),
    np.array([peak, noise, noise]),
)

mz_intensities_only_mod12_additional_noise = (
    np.array([mz_mod_1, 150.0, mz_mod_2, 250.0, 400.0]),
    np.array([peak, noise, peak, noise, noise]),
)


# ------------------------ Setup: correct result DataFrames ---------------------------------------


empty_results_df = pd.DataFrame(
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

all_ions_from_spectrum_exact_matching_df = pd.DataFrame(
    {
        "spectrum_id": [spectrum_native_id, spectrum_native_id, spectrum_native_id],
        "amino_acid": ["Lysine", "Lysine", "Tyrosine"],
        "mod_name": ["Acetyl", "Methyl", "Phospho"],
        "letter_and_unimod_format_mod": [
            "K(UniMod:1)",
            "K(UniMod:34)",
            "Y(UniMod:21)",
        ],
        "type": ["IM", "NL", "IM"],
        "theoretical_mz": [mz_mod_1, mz_mod_2, mz_mod_3],
        "detected_mz": [mz_mod_1, mz_mod_2, mz_mod_3],
        "detected_intensity": [peak, peak, peak],
    }
)

single_higher_peak_ion_df = pd.DataFrame(
    {
        "spectrum_id": [spectrum_native_id],
        "amino_acid": ["Lysine"],
        "mod_name": ["Methyl"],
        "letter_and_unimod_format_mod": [
            "K(UniMod:34)",
        ],
        "type": ["NL"],
        "theoretical_mz": [mz_mod_2],
        "detected_mz": [mz_mod_2],
        "detected_intensity": [peak_higher],
    }
)

within_mass_tolerance_df = pd.DataFrame(
    {
        "spectrum_id": [spectrum_native_id],
        "amino_acid": ["Lysine"],
        "mod_name": ["Acetyl"],
        "letter_and_unimod_format_mod": [
            "K(UniMod:1)",
        ],
        "type": ["IM"],
        "theoretical_mz": [mz_mod_1],
        "detected_mz": [mz_mod_1_threshold_border_higher],
        "detected_intensity": [peak],
    }
)

within_mass_tolerance_multiple_df = pd.DataFrame(
    {
        "spectrum_id": [spectrum_native_id],
        "amino_acid": ["Lysine"],
        "mod_name": ["Acetyl"],
        "letter_and_unimod_format_mod": [
            "K(UniMod:1)",
        ],
        "type": ["IM"],
        "theoretical_mz": [mz_mod_1],
        "detected_mz": [
            mz_mod_1_threshold_border_higher,
        ],
        "detected_intensity": [peak_higher],
    }
)

multiple_spectra_df = pd.DataFrame(
    {
        "spectrum_id": [
            spectrum_native_id,
            spectrum_native_id,
            spectrum_2_native_id,
        ],
        "amino_acid": ["Lysine", "Lysine", "Lysine"],
        "mod_name": ["Acetyl", "Methyl", "Acetyl"],
        "letter_and_unimod_format_mod": [
            "K(UniMod:1)",
            "K(UniMod:34)",
            "K(UniMod:1)",
        ],
        "type": ["IM", "NL", "IM"],
        "theoretical_mz": [mz_mod_1, mz_mod_2, mz_mod_1],
        "detected_mz": [mz_mod_1, mz_mod_2, mz_mod_1],
        "detected_intensity": [peak, peak, peak],
    }
)


# ------------------------ Setup: diagnostic ion detectors ---------------------------------------


@pytest.fixture
def detector_exact_matching_snr_threshold() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        0,
        unit_ppm,
        snr_threshold,
        COLLISION_ENERGY_HIGHER,
    )


@pytest.fixture
def detector_exact_matching_higher_snr_threshold() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        0,
        unit_ppm,
        snr_threshold_higher,
        COLLISION_ENERGY_HIGHER,
    )


@pytest.fixture
def detector_exact_matching_no_snr_threshold() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        0,
        unit_ppm,
        0,
        COLLISION_ENERGY_HIGHER,
    )


@pytest.fixture
def detector_with_ppm_tolerance_snr_threshold() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        mass_tolerance_ppm,
        unit_ppm,
        snr_threshold,
        COLLISION_ENERGY_HIGHER,
    )


@pytest.fixture
def detector_with_ppm_tolerance_no_snr_threshold() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        mass_tolerance_ppm,
        unit_ppm,
        0,
        COLLISION_ENERGY_HIGHER,
    )


@pytest.fixture
def detector_with_da_tolerance_no_snr_threshold() -> DiagnosticIonDetector:
    return DiagnosticIonDetector(
        known_ions_file,
        mass_tolerance_da,
        unit_da,
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


# ------------------------ Tests: ion detection SNR threshold ---------------------------------------


@pytest.fixture(
    params=[
        mz_intensities_ions_additional_noise,
        mz_intensities_ions_additional_noise_irrelevant_peaks,
    ],
    ids=[
        "ions_additional_noise",
        "ions_additional_noise_irrelevant_peaks",
    ],
)
def higher_energy_spectrum_exact_matching(request: SubRequest) -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id, peaks_mz_intensities=request.param
    )


def test_detect_ions_spectrum_exact_matching_applied_snr_threshold(
    higher_energy_spectrum_exact_matching: MSSpectrum,
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """The diagnostic ions should be detected in the spectrum, additional noise
    or peaks at irrelevant positions should be ignored."""
    assert_detection_results_correct(
        detector_exact_matching_snr_threshold,
        higher_energy_spectrum_exact_matching,
        all_ions_from_spectrum_exact_matching_df,
    )


def test_exact_matching_higher_snr_threshold(
    higher_energy_spectrum_exact_matching: MSSpectrum,
    detector_exact_matching_higher_snr_threshold: DiagnosticIonDetector,
) -> None:
    """Applying a higher SNR threshold should lead to no detections if the
    spectrum does not pass the SNR threshold anymore."""
    assert_detection_results_correct(
        detector_exact_matching_higher_snr_threshold,
        higher_energy_spectrum_exact_matching,
        empty_results_df,
    )


@pytest.fixture
def higher_energy_spectrum_single_higher_peak() -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id,
        peaks_mz_intensities=mz_intensities_single_higher_peak_additional_noise,
    )


def test_detect_single_higher_peak_applied_snr_threshold(
    higher_energy_spectrum_single_higher_peak: MSSpectrum,
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """When a peak is much higher than the others, the lower peaks should
    not be detected even if the are higher than some other peaks."""
    assert_detection_results_correct(
        detector_exact_matching_snr_threshold,
        higher_energy_spectrum_single_higher_peak,
        single_higher_peak_ion_df,
    )


@pytest.fixture(
    params=[
        mz_intensities_empty,
        mz_intensities_noise_at_ion_positions,
    ],
    ids=[
        "empty",
        "noise_at_ion_positions",
    ],
)
def higher_energy_spectrum_no_ions(request: SubRequest) -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id, peaks_mz_intensities=request.param
    )


def test_no_ions_in_spectrum_applied_snr_threshold(
    higher_energy_spectrum_no_ions: MSSpectrum,
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """If there are no ions in the spectra (only noise), no ions should be detected. Further,
    empty m/z and intensity arrays should not lead to errors."""
    assert_detection_results_correct(
        detector_exact_matching_snr_threshold,
        higher_energy_spectrum_no_ions,
        empty_results_df,
    )


# ------------------------ Tests: ion detection mass tolerance ---------------------------------------


@pytest.fixture
def higher_energy_spectrum_within_mass_tolerance() -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id,
        peaks_mz_intensities=mz_intensities_within_mass_tolerance_irrelevant_peaks,
    )


def test_within_ppm_tolerance_detected(
    higher_energy_spectrum_within_mass_tolerance: MSSpectrum,
    detector_with_ppm_tolerance_no_snr_threshold: DiagnosticIonDetector,
) -> None:
    """When a ppm tolerance is set, peaks within the tolerance range
    should be detected, other peaks should still be ignored."""
    assert_detection_results_correct(
        detector_with_ppm_tolerance_no_snr_threshold,
        higher_energy_spectrum_within_mass_tolerance,
        within_mass_tolerance_df,
    )


def test_within_da_tolerance_detected(
    higher_energy_spectrum_within_mass_tolerance: MSSpectrum,
    detector_with_da_tolerance_no_snr_threshold: DiagnosticIonDetector,
) -> None:
    """When a Da tolerance is set, peaks within the tolerance range
    should be detected."""
    assert_detection_results_correct(
        detector_with_da_tolerance_no_snr_threshold,
        higher_energy_spectrum_within_mass_tolerance,
        within_mass_tolerance_df,
    )


@pytest.fixture
def higher_energy_spectrum_within_mass_tolerance_multiple() -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id,
        peaks_mz_intensities=mz_intensities_within_mass_tolerance_multiple,
    )


def test_within_ppm_tolerance_multiple_detected(
    higher_energy_spectrum_within_mass_tolerance_multiple: MSSpectrum,
    detector_with_ppm_tolerance_no_snr_threshold: DiagnosticIonDetector,
) -> None:
    """When a ppm tolerance is set, and multiple peaks are within the tolerance
    of a known ion m/z, the highest peak should be detected."""
    assert_detection_results_correct(
        detector_with_ppm_tolerance_no_snr_threshold,
        higher_energy_spectrum_within_mass_tolerance_multiple,
        within_mass_tolerance_multiple_df,
    )


@pytest.fixture
def higher_energy_spectrum_within_ppm_tolerance_additional_noise() -> MSSpectrum:
    return spectrum_ms2_higher_energy(
        native_id=spectrum_native_id,
        peaks_mz_intensities=mz_intensities_within_mass_tolerance_additional_noise,
    )


def test_within_ppm_tolerance_detected_additional_noise(
    higher_energy_spectrum_within_ppm_tolerance_additional_noise: MSSpectrum,
    detector_with_ppm_tolerance_snr_threshold: DiagnosticIonDetector,
) -> None:
    """Tolerance should also work simulaneously with noise reduction."""
    assert_detection_results_correct(
        detector_with_ppm_tolerance_snr_threshold,
        higher_energy_spectrum_within_ppm_tolerance_additional_noise,
        within_mass_tolerance_df,
    )


# ------------------------ Tests: error cases ---------------------------------------


@pytest.fixture
def ms1_spectrum() -> MSSpectrum:
    return spectrum_ms1(native_id=spectrum_native_id)


def test_ion_detection_ms1_fails(
    ms1_spectrum: MSSpectrum,
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """Diagnostic ion detection for an MS1 scan should fail."""
    with pytest.raises(AssertionError, match="MS2"):
        detector_exact_matching_snr_threshold.extract_diagnostic_ions_for_spectrum(
            ms1_spectrum
        )


@pytest.fixture
def lower_energy_spectrum() -> MSSpectrum:
    return spectrum_ms2_lower_energy(native_id=spectrum_native_id)


def test_ion_detection_lower_energy_fails(
    lower_energy_spectrum: MSSpectrum,
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """Diagnostic ion detection for a lower-energy MS2 scan should fail."""
    with pytest.raises(AssertionError, match="collision energy"):
        detector_exact_matching_snr_threshold.extract_diagnostic_ions_for_spectrum(
            lower_energy_spectrum
        )


# ------------------------ Tests: multiple spectra ---------------------------------------


@pytest.fixture
def higher_energy_multiple_spectra_exact_matching() -> List[MSSpectrum]:
    return [
        spectrum_ms2_higher_energy(
            native_id=spectrum_native_id,
            peaks_mz_intensities=mz_intensities_only_mod12_additional_noise,
        ),
        spectrum_ms2_higher_energy(
            native_id=spectrum_2_native_id,
            peaks_mz_intensities=mz_intensities_only_mod1_additional_noise,
        ),
    ]


def test_detect_ions_multiple_spectra_exact_matching(
    higher_energy_multiple_spectra_exact_matching: List[MSSpectrum],
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """The diagnostic ions should be detected in the spectra and stored with the corresponding
    spectrum id."""
    detected_ions_df = (
        detector_exact_matching_snr_threshold.extract_diagnostic_ions_for_spectra(
            higher_energy_multiple_spectra_exact_matching
        )
    )
    pd.testing.assert_frame_equal(
        detected_ions_df,
        multiple_spectra_df,
        # not comparing the dtype to avoid failing because of float64 VS float32
        check_dtype=False,
    )


def test_validate_spectra_successful(
    higher_energy_multiple_spectra_exact_matching: List[MSSpectrum],
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """Validation of a list of higher-energy MS2 spectra should succeed."""
    detector_exact_matching_snr_threshold.validate_diagnostic_ion_spectra(
        higher_energy_multiple_spectra_exact_matching
    )


def test_validate_MS1_spectra_fails(
    ms1_spectrum: MSSpectrum,
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """Validation of a list of MS1 spectra should fail."""
    with pytest.raises(AssertionError, match="MS2"):
        detector_exact_matching_snr_threshold.validate_diagnostic_ion_spectra(
            [ms1_spectrum]
        )


def test_validate_lower_energy_spectra_fails(
    lower_energy_spectrum: MSSpectrum,
    detector_exact_matching_snr_threshold: DiagnosticIonDetector,
) -> None:
    """Validation of a list of lower-energy spectra should fail."""
    with pytest.raises(AssertionError, match="collision energy"):
        detector_exact_matching_snr_threshold.validate_diagnostic_ion_spectra(
            [lower_energy_spectrum]
        )
