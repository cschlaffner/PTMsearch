from typing import FrozenSet, List, Tuple

import pandas as pd
import pytest
from _pytest.fixtures import SubRequest

from src.split_processing.spectral_library_splitting import split_library_by_mods

k_1 = "K(UniMod:1)"
k_34 = "K(UniMod:34)"
y_21 = "Y(UniMod:21)"
m_35 = "M(UniMod:35)"

# TODO: test other PTM format


@pytest.fixture
def mods_single() -> Tuple[List[str], List[FrozenSet[str]], List[str]]:
    return [k_1, k_34], [frozenset([k_1, y_21]), frozenset([k_1, y_21, k_34])], []


@pytest.fixture
def mods_combinations() -> Tuple[List[str], List[FrozenSet[str]], List[str]]:
    return [k_1, k_34], [frozenset([k_1, y_21]), frozenset([k_1, y_21, k_34])], []


@pytest.fixture
def mods_combinations_additional() -> Tuple[List[str], List[FrozenSet[str]], List[str]]:
    return [k_1, k_34], [frozenset([k_1, y_21]), frozenset([k_1, y_21, k_34])], []


@pytest.fixture
def spectral_library_df() -> pd.DataFrame:
    return pd.DataFrame({"ModifiedPeptide": []})


@pytest.mark.parametrize(
    "mods_to_search_fixture, expected_result_libs_fixture",
    [("mods_single",), ("mods_combinations",), ("mods_combinations_additional",)],
)
def test_library_splitting_unimod(
    mods_to_search_fixture: str,
    expected_result_libs_fixture: str,
    request: SubRequest,
    spectral_library_df: pd.DataFrame,
) -> None:
    mods, mod_combinations, additional_mods = request.getfixturevalue(
        mods_to_search_fixture
    )
    expected_result_libs = request.getfixturevalue(expected_result_libs_fixture)

    result_libs = split_library_by_mods(
        spectral_library_df, True, mods, mods_combinations, additional_mods
    )

    assert result_libs.keys() == expected_result_libs.keys()
    for mod_split in result_libs:
        assert result_libs[mod_split].equals(expected_result_libs[mod_split])
