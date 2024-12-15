from typing import FrozenSet, List, Tuple

import numpy as np
import pandas as pd
import pytest
from _pytest.fixtures import SubRequest

from src.diagnostic_ions.summary import get_detected_modifications_with_combinations

higher_energy_spectrum_1_native_id = "ms2_id_1"
higher_energy_spectrum_2_native_id = "ms2_id_2"
higher_energy_spectrum_3_native_id = "ms2_id_3"
higher_energy_spectrum_4_native_id = "ms2_id_4"
higher_energy_spectrum_5_native_id = "ms2_id_5"
higher_energy_spectrum_6_native_id = "ms2_id_6"

k_1 = "K(UniMod:1)"
k_34 = "K(UniMod:34)"
y_21 = "Y(UniMod:21)"

window_minimum_threshold = 2
percentile_threshold = 2 / 3
additional_mods_count = 4


@pytest.fixture
def detected_ions_df() -> pd.DataFrame:
    # K1
    # K1
    # Y21, K1
    # K1, L34, Y21
    # K34
    # K1, Y21
    return pd.DataFrame(
        {
            "spectrum_id": [
                higher_energy_spectrum_1_native_id,
                higher_energy_spectrum_2_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_5_native_id,
                higher_energy_spectrum_6_native_id,
                higher_energy_spectrum_6_native_id,
            ],
            "amino_acid": [
                "Lysine",
                "Lysine",
                "Tyrosine",
                "Lysine",
                "Lysine",
                "Lysine",
                "Tyrosine",
                "Lysine",
                "Lysine",
                "Tyrosine",
            ],
            "mod_name": [
                "Acetyl",
                "Acetyl",
                "Phospho",
                "Acetyl",
                "Acetyl",
                "Methyl",
                "Phospho",
                "Methyl",
                "Acetyl",
                "Phospho",
            ],
            "letter_and_unimod_format_mod": [
                k_1,
                k_1,
                y_21,
                k_1,
                k_1,
                k_34,
                y_21,
                k_34,
                k_1,
                y_21,
            ],
        }
    )


@pytest.fixture
def mods_combinations_no_threshold() -> Tuple[List[str], List[FrozenSet[str]]]:
    return [k_1, k_34], [frozenset([k_1, y_21]), frozenset([k_1, y_21, k_34])]


@pytest.fixture
def mods_combinations_no_threshold_additional_mods() -> (
    Tuple[List[str], List[FrozenSet[str]]]
):
    return [k_1, k_34, y_21], []


@pytest.fixture
def mods_combinations_window_threshold() -> Tuple[List[str], List[FrozenSet[str]]]:
    return [k_1, k_34], [frozenset([k_1, y_21])]


@pytest.fixture
def mods_combinations_percentile_threshold() -> Tuple[List[str], List[FrozenSet[str]]]:
    return [k_1, k_34, y_21], [frozenset([k_1, y_21])]


@pytest.mark.parametrize(
    "detection_count_percentile, detection_count_min, num_additional_modifications, mods_combinations_result_fixture",
    [
        (1.0, 0, 0, "mods_combinations_no_threshold"),
        (
            1.0,
            0,
            additional_mods_count,
            "mods_combinations_no_threshold_additional_mods",
        ),
        (1.0, window_minimum_threshold, 0, "mods_combinations_window_threshold"),
        (percentile_threshold, 0, 0, "mods_combinations_percentile_threshold"),
        (
            percentile_threshold,
            window_minimum_threshold,
            0,
            "mods_combinations_window_threshold",
        ),
    ],
)
def test_mod_selection(
    detection_count_percentile: float,
    detection_count_min: int,
    num_additional_modifications: int,
    mods_combinations_result_fixture: str,
    request: SubRequest,
    detected_ions_df: pd.DataFrame,
) -> None:
    expected_mods, expected_combinations = request.getfixturevalue(
        mods_combinations_result_fixture
    )
    mods, combinations = get_detected_modifications_with_combinations(
        detected_ions_df,
        detection_count_percentile,
        detection_count_min,
        num_additional_modifications,
    )
    assert np.array_equal(sorted(mods), sorted(expected_mods))
    assert np.array_equal(sorted(combinations), sorted(expected_combinations))
