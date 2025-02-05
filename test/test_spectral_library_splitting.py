from typing import Dict, FrozenSet, List, Tuple

import pandas as pd
import pytest
from _pytest.fixtures import SubRequest

from src.split_processing.spectral_library_splitting import split_library_by_mods

k_1 = "K(UniMod:1)"
k_1_mass_diff = "K[42.010565]"
k_3 = "K(UniMod:3)"
k_34 = "K(UniMod:34)"
y_21 = "Y(UniMod:21)"
y_21_mass_diff = "Y[79.966331]"
m_35 = "M(UniMod:35)"

k_1_y_21 = frozenset([k_1, y_21])
k_3_y_21 = frozenset([k_3, y_21])
k_1_y_21_k_34 = frozenset([k_1, y_21, k_34])


@pytest.fixture
def mods_single() -> Tuple[List[str], List[FrozenSet[str]], List[str]]:
    return [k_1, y_21, k_34], [], []


@pytest.fixture
def mods_combinations() -> Tuple[List[str], List[FrozenSet[str]], List[str]]:
    return [k_1, y_21, k_34], [k_1_y_21, k_3_y_21, k_1_y_21_k_34], []


@pytest.fixture
def mods_combinations_reduced() -> Tuple[List[str], List[FrozenSet[str]], List[str]]:
    return [k_1, y_21], [k_1_y_21, k_3_y_21, k_1_y_21_k_34], []


@pytest.fixture
def mods_combinations_additional() -> Tuple[List[str], List[FrozenSet[str]], List[str]]:
    return (
        [k_1, y_21, k_34],
        [k_1_y_21, k_3_y_21, k_1_y_21_k_34],
        [m_35],
    )


@pytest.fixture
def spectral_library_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "ModifiedPeptide": [
                "ABCDEF",
                f"ABC{k_1}A",
                f"ABC{k_1}A",
                f"ABC{k_1}B{k_1}",
                f"ABC{k_34}A",
                f"ABC{k_1}A{y_21}",
                f"ABC{y_21}A{k_1}",
                f"ABC{y_21}A{k_34}",
                f"ABC{y_21}A{k_34}{k_1}",
                f"ABC{k_1}B{m_35}",
                f"ABC{m_35}",
            ]
        }
    )


# pylint: disable=unnecessary-lambda-assignment
library_unmod = lambda: pd.DataFrame({"ModifiedPeptide": ["ABCDEF"]})

library_unmod_m_35 = lambda: pd.DataFrame({"ModifiedPeptide": ["ABCDEF", f"ABC{m_35}"]})

library_k_1 = lambda: pd.DataFrame(
    {"ModifiedPeptide": [f"ABC{k_1}A", f"ABC{k_1}A", f"ABC{k_1}B{k_1}", "ABCDEF"]}
)
library_k_1_m_35 = lambda: pd.DataFrame(
    {
        "ModifiedPeptide": [
            f"ABC{k_1}A",
            f"ABC{k_1}A",
            f"ABC{k_1}B{k_1}",
            f"ABC{k_1}B{m_35}",
            "ABCDEF",
            f"ABC{m_35}",
        ]
    }
)

library_y_21 = lambda: pd.DataFrame({"ModifiedPeptide": ["ABCDEF"]})
library_y_21_m_35 = lambda: pd.DataFrame({"ModifiedPeptide": ["ABCDEF", f"ABC{m_35}"]})

library_k_34 = lambda: pd.DataFrame({"ModifiedPeptide": [f"ABC{k_34}A", "ABCDEF"]})
library_k_34_m_35 = lambda: pd.DataFrame(
    {"ModifiedPeptide": [f"ABC{k_34}A", "ABCDEF", f"ABC{m_35}"]}
)

library_k_1_y_21 = lambda: pd.DataFrame(
    {
        "ModifiedPeptide": [
            f"ABC{k_1}A",
            f"ABC{k_1}A",
            f"ABC{k_1}B{k_1}",
            f"ABC{k_1}A{y_21}",
            f"ABC{y_21}A{k_1}",
            "ABCDEF",
        ]
    }
)
library_k_1_y_21_m_35 = lambda: pd.DataFrame(
    {
        "ModifiedPeptide": [
            f"ABC{k_1}A",
            f"ABC{k_1}A",
            f"ABC{k_1}B{k_1}",
            f"ABC{k_1}A{y_21}",
            f"ABC{y_21}A{k_1}",
            f"ABC{k_1}B{m_35}",
            "ABCDEF",
            f"ABC{m_35}",
        ]
    }
)

library_k_3_y_21 = lambda: pd.DataFrame({"ModifiedPeptide": ["ABCDEF"]})
library_k_3_y_21_m_35 = lambda: pd.DataFrame(
    {"ModifiedPeptide": ["ABCDEF", f"ABC{m_35}"]}
)

library_k_1_y_21_k_34 = lambda: pd.DataFrame(
    {
        "ModifiedPeptide": [
            f"ABC{k_1}A",
            f"ABC{k_1}A",
            f"ABC{k_1}B{k_1}",
            f"ABC{k_34}A",
            f"ABC{k_1}A{y_21}",
            f"ABC{y_21}A{k_1}",
            f"ABC{y_21}A{k_34}",
            f"ABC{y_21}A{k_34}{k_1}",
            "ABCDEF",
        ]
    }
)
library_k_1_y_21_k_34_m_35 = lambda: pd.DataFrame(
    {
        "ModifiedPeptide": [
            f"ABC{k_1}A",
            f"ABC{k_1}A",
            f"ABC{k_1}B{k_1}",
            f"ABC{k_34}A",
            f"ABC{k_1}A{y_21}",
            f"ABC{y_21}A{k_1}",
            f"ABC{y_21}A{k_34}",
            f"ABC{y_21}A{k_34}{k_1}",
            f"ABC{k_1}B{m_35}",
            "ABCDEF",
            f"ABC{m_35}",
        ]
    }
)


@pytest.fixture
def libraries_by_mods_single() -> Dict[str, pd.DataFrame]:
    return {
        k_1: library_k_1(),
        y_21: library_y_21(),
        k_34: library_k_34(),
        "unmodified": library_unmod(),
    }


@pytest.fixture
def libraries_by_mods_combinations() -> Dict[str, pd.DataFrame]:
    return {
        k_1: library_k_1(),
        y_21: library_y_21(),
        k_34: library_k_34(),
        k_1_y_21: library_k_1_y_21(),
        k_3_y_21: library_k_3_y_21(),
        k_1_y_21_k_34: library_k_1_y_21_k_34(),
        "unmodified": library_unmod(),
    }


@pytest.fixture
def libraries_by_mods_combinations_reduced() -> Dict[str, pd.DataFrame]:
    return {
        k_1: library_k_1(),
        y_21: library_y_21(),
        k_1_y_21: library_k_1_y_21(),
        k_3_y_21: library_k_3_y_21(),
        k_1_y_21_k_34: library_k_1_y_21_k_34(),
        "unmodified": library_unmod(),
    }


@pytest.fixture
def libraries_by_mods_combinations_additional() -> Dict[str, pd.DataFrame]:
    return {
        k_1: library_k_1_m_35(),
        y_21: library_y_21_m_35(),
        k_34: library_k_34_m_35(),
        k_1_y_21: library_k_1_y_21_m_35(),
        k_3_y_21: library_k_3_y_21_m_35(),
        k_1_y_21_k_34: library_k_1_y_21_k_34_m_35(),
        "unmodified": library_unmod_m_35(),
    }


@pytest.mark.parametrize(
    "mods_to_search_fixture, expected_result_libs_fixture",
    [
        ("mods_single", "libraries_by_mods_single"),
        ("mods_combinations", "libraries_by_mods_combinations"),
        ("mods_combinations_reduced", "libraries_by_mods_combinations_reduced"),
        ("mods_combinations_additional", "libraries_by_mods_combinations_additional"),
    ],
)
def test_library_splitting_unimod(
    mods_to_search_fixture: str,
    expected_result_libs_fixture: str,
    request: SubRequest,
    spectral_library_df: pd.DataFrame,
) -> None:
    mods, combinations, additional_mods = request.getfixturevalue(
        mods_to_search_fixture
    )
    expected_result_libs = request.getfixturevalue(expected_result_libs_fixture)

    result_libs = split_library_by_mods(
        spectral_library_df, mods, combinations, additional_mods
    )
    assert result_libs.keys() == expected_result_libs.keys()
    for mod_split, result_lib in result_libs.items():
        assert result_lib.equals(expected_result_libs[mod_split])


@pytest.fixture
def spectral_library_df_mass_diff_format() -> pd.DataFrame:
    return pd.DataFrame(
        {"ModifiedPeptide": [f"ABC{k_1_mass_diff}A", f"ABC{y_21_mass_diff}A"]}
    )


@pytest.fixture
def libraries_by_mods_mass_diff() -> Dict[str, pd.DataFrame]:
    return {
        k_1: pd.DataFrame({"ModifiedPeptide": [f"ABC{k_1_mass_diff}A"]}),
        y_21: pd.DataFrame({"ModifiedPeptide": [f"ABC{y_21_mass_diff}A"]}),
        "unmodified": pd.DataFrame(),
    }


def test_library_splitting_mass_diff(
    spectral_library_df_mass_diff_format: pd.DataFrame,
    libraries_by_mods_mass_diff: Dict[str, pd.DataFrame],
) -> None:
    result_libs = split_library_by_mods(
        spectral_library_df_mass_diff_format,
        [k_1, y_21],
        mod_combinations_to_search=[],
        additional_mods_to_search=[],
        has_unimod_format=False,
    )
    assert result_libs.keys() == libraries_by_mods_mass_diff.keys()
    for mod_split, result_lib in result_libs.items():
        assert result_lib.equals(libraries_by_mods_mass_diff[mod_split])
