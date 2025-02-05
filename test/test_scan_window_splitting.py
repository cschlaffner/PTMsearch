from test.spectra import (
    COLLISION_ENERGY_HIGHER,
    COLLISION_ENERGY_LOWER,
    spectrum_ms1,
    spectrum_ms2_higher_energy,
    spectrum_ms2_lower_energy,
)
from typing import Dict, FrozenSet, List, Tuple, Union

import pandas as pd
import pytest
from _pytest.fixtures import SubRequest
from pyopenms import MSExperiment, MSSpectrum

from src.mzml_processing.utils import (
    get_ms2_spectrum_collision_energy,
    get_ms2_spectrum_precursors,
)
from src.split_processing.scan_window_splitting import ScanWindowSplitting

ms1_spectrum_1_native_id = "ms1_id_1"
ms1_spectrum_2_native_id = "ms1_id_2"
ms1_spectrum_3_native_id = "ms1_id_10"

higher_energy_spectrum_1_native_id = "ms2_id_1"
higher_energy_spectrum_2_native_id = "ms2_id_2"
higher_energy_spectrum_3_native_id = "ms2_id_3"
higher_energy_spectrum_4_native_id = "ms2_id_10"
higher_energy_spectrum_5_native_id = "ms2_id_5"

precursor_mz_1 = 1.1
precursor_mz_2 = 20.2


@pytest.fixture
def scan_window_splitting() -> ScanWindowSplitting:
    return ScanWindowSplitting(COLLISION_ENERGY_LOWER, COLLISION_ENERGY_HIGHER)


@pytest.fixture
def mods_empty() -> List[str]:
    return []


@pytest.fixture
def mods_single() -> List[str]:
    return ["K(UniMod:1)"]


@pytest.fixture
def mods_multiple() -> List[str]:
    return ["K(UniMod:1)", "K(UniMod:34)", "Y(UniMod:21)"]


@pytest.fixture
def mod_combinations_some() -> List[FrozenSet[str]]:
    return [{"K(UniMod:1)", "K(UniMod:34)", "Y(UniMod:21)"}]


@pytest.fixture
def mod_combinations_all() -> List[FrozenSet[str]]:
    return [
        {"K(UniMod:1)", "Y(UniMod:21)"},
        {"K(UniMod:1)", "K(UniMod:34)", "Y(UniMod:21)"},
    ]


def spectra_ms1() -> MSExperiment:
    exp_ms1_spectra = MSExperiment()
    ms1_spectra = [
        spectrum_ms1(native_id=ms1_spectrum_1_native_id),
        spectrum_ms1(native_id=ms1_spectrum_2_native_id),
        spectrum_ms1(native_id=ms1_spectrum_3_native_id),
    ]
    exp_ms1_spectra.setSpectra(ms1_spectra)
    return exp_ms1_spectra


def spectra_ms1_and_lower_energy_matching() -> MSExperiment:
    exp_ms1_spectra = spectra_ms1()
    ms1_spectra = exp_ms1_spectra.getSpectra()

    exp_ms1_and_lower_energy_spectra = MSExperiment()
    exp_ms1_and_lower_energy_spectra.setSpectra(
        [
            ms1_spectra[0],
            spectrum_ms2_lower_energy(precursor_mz=precursor_mz_1),
            spectrum_ms2_lower_energy(precursor_mz=precursor_mz_2),
            ms1_spectra[1],
            spectrum_ms2_lower_energy(precursor_mz=precursor_mz_1),
            spectrum_ms2_lower_energy(precursor_mz=precursor_mz_2),
            ms1_spectra[2],
            spectrum_ms2_lower_energy(precursor_mz=precursor_mz_1),
        ]
    )
    return exp_ms1_and_lower_energy_spectra


def spectra_ms1_and_higher_energy_matching() -> MSExperiment:
    exp_ms1_spectra = spectra_ms1()
    ms1_spectra = exp_ms1_spectra.getSpectra()

    exp_ms1_and_higher_energy_spectra = MSExperiment()
    exp_ms1_and_higher_energy_spectra.setSpectra(
        [
            ms1_spectra[0],
            spectrum_ms2_higher_energy(
                native_id=higher_energy_spectrum_1_native_id,
                precursor_mz=precursor_mz_1,
            ),
            spectrum_ms2_higher_energy(
                native_id=higher_energy_spectrum_2_native_id,
                precursor_mz=precursor_mz_2,
            ),
            ms1_spectra[1],
            spectrum_ms2_higher_energy(
                native_id=higher_energy_spectrum_3_native_id,
                precursor_mz=precursor_mz_1,
            ),
            spectrum_ms2_higher_energy(
                native_id=higher_energy_spectrum_4_native_id,
                precursor_mz=precursor_mz_2,
            ),
            ms1_spectra[2],
            spectrum_ms2_higher_energy(
                native_id=higher_energy_spectrum_5_native_id,
                precursor_mz=precursor_mz_1,
            ),
        ]
    )
    return exp_ms1_and_higher_energy_spectra


@pytest.fixture
def spectra_matching() -> Tuple[MSExperiment]:
    return (
        spectra_ms1(),
        spectra_ms1_and_lower_energy_matching(),
        spectra_ms1_and_higher_energy_matching(),
    )


@pytest.fixture
def detected_ions_df_single_mod() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "spectrum_id": [
                higher_energy_spectrum_1_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_4_native_id,
            ],
            "amino_acid": ["Lysine", "Lysine", "Lysine"],
            "mod_name": ["Acetyl", "Acetyl", "Acetyl"],
            "letter_and_unimod_format_mod": [
                "K(UniMod:1)",
                "K(UniMod:1)",
                "K(UniMod:1)",
            ],
        }
    )


@pytest.fixture
def spectra_matching_with_result_single_mod() -> (
    Tuple[Union[MSExperiment, Dict[str, pd.DataFrame]]]
):
    exp_ms1_and_lower_energy_spectra = spectra_ms1_and_lower_energy_matching()
    ms1_and_lower_energy_spectra = exp_ms1_and_lower_energy_spectra.getSpectra()
    unmod_exp = MSExperiment()
    unmod_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[2],
        ms1_and_lower_energy_spectra[6],
        ms1_and_lower_energy_spectra[7],
    ]
    unmod_exp.setSpectra(unmod_spectra)

    mod_exp = MSExperiment()
    mod_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[1],
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[4],
        ms1_and_lower_energy_spectra[5],
    ]
    mod_exp.setSpectra(mod_spectra)

    return (
        spectra_ms1(),
        exp_ms1_and_lower_energy_spectra,
        spectra_ms1_and_higher_energy_matching(),
        {"unmodified": unmod_exp, "K(UniMod:1)": mod_exp},
    )


@pytest.fixture
def detected_ions_df_multiple_mods() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "spectrum_id": [
                higher_energy_spectrum_1_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_5_native_id,
            ],
            "amino_acid": [
                "Lysine",
                "Tyrosine",
                "Lysine",
                "Lysine",
                "Lysine",
                "Tyrosine",
                "Tyrosine",
            ],
            "mod_name": [
                "Acetyl",
                "Phospho",
                "Acetyl",
                "Acetyl",
                "Methyl",
                "Phospho",
                "Phospho",
            ],
            "letter_and_unimod_format_mod": [
                "K(UniMod:1)",
                "Y(UniMod:21)",
                "K(UniMod:1)",
                "K(UniMod:1)",
                "K(UniMod:34)",
                "Y(UniMod:21)",
                "Y(UniMod:21)",
            ],
        }
    )


@pytest.fixture
def detected_ions_df_multiple_mods_duplicate_entries() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "spectrum_id": [
                higher_energy_spectrum_1_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_3_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_4_native_id,
                higher_energy_spectrum_5_native_id,
            ],
            "amino_acid": [
                "Lysine",
                "Tyrosine",
                "Tyrosine",
                "Tyrosine",
                "Lysine",
                "Lysine",
                "Lysine",
                "Lysine",
                "Tyrosine",
                "Tyrosine",
            ],
            "mod_name": [
                "Acetyl",
                "Phospho",
                "Phospho",
                "Phospho",
                "Acetyl",
                "Acetyl",
                "Methyl",
                "Methyl",
                "Phospho",
                "Phospho",
            ],
            "letter_and_unimod_format_mod": [
                "K(UniMod:1)",
                "Y(UniMod:21)",
                "Y(UniMod:21)",
                "Y(UniMod:21)",
                "K(UniMod:1)",
                "K(UniMod:1)",
                "K(UniMod:34)",
                "K(UniMod:34)",
                "Y(UniMod:21)",
                "Y(UniMod:21)",
            ],
        }
    )


@pytest.fixture
def spectra_matching_with_result_multiple_mods() -> (
    Tuple[Union[MSExperiment, Dict[str, pd.DataFrame]]]
):
    exp_ms1_and_lower_energy_spectra = spectra_ms1_and_lower_energy_matching()
    ms1_and_lower_energy_spectra = exp_ms1_and_lower_energy_spectra.getSpectra()
    unmod_exp = MSExperiment()
    unmod_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[2],
    ]
    unmod_exp.setSpectra(unmod_spectra)

    mod_1_exp = MSExperiment()
    mod_1_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[1],
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[4],
        ms1_and_lower_energy_spectra[5],
    ]
    mod_1_exp.setSpectra(mod_1_spectra)

    mod_2_exp = MSExperiment()
    mod_2_spectra = [
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[5],
    ]
    mod_2_exp.setSpectra(mod_2_spectra)

    mod_3_exp = MSExperiment()
    mod_3_spectra = [
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[4],
        ms1_and_lower_energy_spectra[5],
        ms1_and_lower_energy_spectra[6],
        ms1_and_lower_energy_spectra[7],
    ]
    mod_3_exp.setSpectra(mod_3_spectra)

    return (
        spectra_ms1(),
        exp_ms1_and_lower_energy_spectra,
        spectra_ms1_and_higher_energy_matching(),
        {
            "unmodified": unmod_exp,
            "K(UniMod:1)": mod_1_exp,
            "K(UniMod:34)": mod_2_exp,
            "Y(UniMod:21)": mod_3_exp,
        },
    )


@pytest.fixture
def spectra_matching_with_result_multiple_mods_combinations_and_single() -> (
    Tuple[Union[MSExperiment, Dict[str, pd.DataFrame]]]
):
    exp_ms1_and_lower_energy_spectra = spectra_ms1_and_lower_energy_matching()
    ms1_and_lower_energy_spectra = exp_ms1_and_lower_energy_spectra.getSpectra()
    unmod_exp = MSExperiment()
    unmod_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[2],
    ]
    unmod_exp.setSpectra(unmod_spectra)

    mod_1_exp = MSExperiment()
    mod_1_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[1],
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[4],
    ]
    mod_1_exp.setSpectra(mod_1_spectra)

    mod_2_exp = MSExperiment()
    mod_2_spectra = [
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[4],
        ms1_and_lower_energy_spectra[6],
        ms1_and_lower_energy_spectra[7],
    ]
    mod_2_exp.setSpectra(mod_2_spectra)

    mod_3_exp = MSExperiment()
    mod_3_spectra = [
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[5],
    ]
    mod_3_exp.setSpectra(mod_3_spectra)

    return (
        spectra_ms1(),
        exp_ms1_and_lower_energy_spectra,
        spectra_ms1_and_higher_energy_matching(),
        {
            "unmodified": unmod_exp,
            "K(UniMod:1)": mod_1_exp,
            "Y(UniMod:21)": mod_2_exp,
            frozenset({"K(UniMod:1)", "K(UniMod:34)", "Y(UniMod:21)"}): mod_3_exp,
        },
    )


@pytest.fixture
def spectra_matching_with_result_multiple_mods_combinations_and_single_reduced() -> (
    Tuple[Union[MSExperiment, Dict[str, pd.DataFrame]]]
):
    exp_ms1_and_lower_energy_spectra = spectra_ms1_and_lower_energy_matching()
    ms1_and_lower_energy_spectra = exp_ms1_and_lower_energy_spectra.getSpectra()
    unmod_exp = MSExperiment()
    unmod_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[2],
        ms1_and_lower_energy_spectra[6],
        ms1_and_lower_energy_spectra[7],
    ]
    unmod_exp.setSpectra(unmod_spectra)

    mod_1_exp = MSExperiment()
    mod_1_spectra = [  # 1, 3, 5
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[1],
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[4],
    ]
    mod_1_exp.setSpectra(mod_1_spectra)

    mod_2_exp = MSExperiment()  # 4
    mod_2_spectra = [
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[5],
    ]
    mod_2_exp.setSpectra(mod_2_spectra)

    return (
        spectra_ms1(),
        exp_ms1_and_lower_energy_spectra,
        spectra_ms1_and_higher_energy_matching(),
        {
            "unmodified": unmod_exp,
            "K(UniMod:1)": mod_1_exp,
            frozenset({"K(UniMod:1)", "K(UniMod:34)", "Y(UniMod:21)"}): mod_2_exp,
        },
    )


@pytest.fixture
def spectra_matching_with_result_multiple_mods_combinations() -> (
    Tuple[Union[MSExperiment, Dict[str, pd.DataFrame]]]
):
    exp_ms1_and_lower_energy_spectra = spectra_ms1_and_lower_energy_matching()
    ms1_and_lower_energy_spectra = exp_ms1_and_lower_energy_spectra.getSpectra()
    unmod_exp = MSExperiment()
    unmod_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[2],
    ]
    unmod_exp.setSpectra(unmod_spectra)

    mod_1_exp = MSExperiment()
    mod_1_spectra = [
        ms1_and_lower_energy_spectra[0],
        ms1_and_lower_energy_spectra[1],
    ]
    mod_1_exp.setSpectra(mod_1_spectra)

    mod_2_exp = MSExperiment()
    mod_2_spectra = [
        ms1_and_lower_energy_spectra[6],
        ms1_and_lower_energy_spectra[7],
    ]
    mod_2_exp.setSpectra(mod_2_spectra)

    mod_3_exp = MSExperiment()
    mod_3_spectra = [
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[4],
    ]
    mod_3_exp.setSpectra(mod_3_spectra)

    mod_4_exp = MSExperiment()
    mod_4_spectra = [
        ms1_and_lower_energy_spectra[3],
        ms1_and_lower_energy_spectra[5],
    ]
    mod_4_exp.setSpectra(mod_4_spectra)

    return (
        spectra_ms1(),
        exp_ms1_and_lower_energy_spectra,
        spectra_ms1_and_higher_energy_matching(),
        {
            "unmodified": unmod_exp,
            "K(UniMod:1)": mod_1_exp,
            "Y(UniMod:21)": mod_2_exp,
            frozenset({"K(UniMod:1)", "Y(UniMod:21)"}): mod_3_exp,
            frozenset({"K(UniMod:1)", "K(UniMod:34)", "Y(UniMod:21)"}): mod_4_exp,
        },
    )


@pytest.fixture
def spectra_mismatch_higher_lower_energy_by_precursor_mz() -> Tuple[MSExperiment]:
    exp_ms1_spectra = spectra_ms1()
    ms1_spectra = exp_ms1_spectra.getSpectra()

    exp_ms1_and_lower_energy_spectra = spectra_ms1_and_lower_energy_matching()

    exp_ms1_and_higher_energy_spectra = MSExperiment()
    exp_ms1_and_higher_energy_spectra.setSpectra(
        [
            ms1_spectra[0],
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_1),
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_2),
            ms1_spectra[1],
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_1),
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_2 + 1),
            ms1_spectra[2],
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_1),
        ]
    )
    return (
        exp_ms1_spectra,
        exp_ms1_and_lower_energy_spectra,
        exp_ms1_and_higher_energy_spectra,
    )


def test_splitting_mismatching_higher_lower_energy_by_precursor_mz_should_fail(
    spectra_mismatch_higher_lower_energy_by_precursor_mz: Tuple[MSExperiment],
    detected_ions_df_single_mod: pd.DataFrame,
    mods_single: List[str],
    scan_window_splitting: ScanWindowSplitting,
) -> None:
    with pytest.raises(AssertionError, match="precursor m/z"):
        scan_window_splitting.split_windows_by_mods(
            *spectra_mismatch_higher_lower_energy_by_precursor_mz,
            detected_ions_df_single_mod,
            mods_single
        )


@pytest.fixture
def spectra_mismatch_higher_lower_energy_by_ms1() -> Tuple[MSExperiment]:
    exp_ms1_spectra = spectra_ms1()
    ms1_spectra = exp_ms1_spectra.getSpectra()

    exp_ms1_and_lower_energy_spectra = spectra_ms1_and_lower_energy_matching()

    exp_ms1_and_higher_energy_spectra = MSExperiment()
    exp_ms1_and_higher_energy_spectra.setSpectra(
        [
            ms1_spectra[0],
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_1),
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_2),
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_1),
            ms1_spectra[1],
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_2),
            ms1_spectra[2],
            spectrum_ms2_higher_energy(precursor_mz=precursor_mz_1),
        ]
    )
    return (
        exp_ms1_spectra,
        exp_ms1_and_lower_energy_spectra,
        exp_ms1_and_higher_energy_spectra,
    )


def test_splitting_mismatching_higher_lower_energy_by_ms1_spectra_should_fail(
    spectra_mismatch_higher_lower_energy_by_ms1: Tuple[MSExperiment],
    detected_ions_df_single_mod: pd.DataFrame,
    mods_single: List[str],
    scan_window_splitting: ScanWindowSplitting,
) -> None:
    with pytest.raises(AssertionError, match="MS1"):
        scan_window_splitting.split_windows_by_mods(
            *spectra_mismatch_higher_lower_energy_by_ms1,
            detected_ions_df_single_mod,
            mods_single
        )


@pytest.fixture
def spectra_mismatching_ms1() -> Tuple[MSExperiment]:
    exp_ms1_spectra_mismatching = MSExperiment()
    exp_ms1_spectra_mismatching.setSpectra(spectra_ms1().getSpectra()[:-1])

    exp_ms1_and_lower_energy_spectra = spectra_ms1_and_lower_energy_matching()

    exp_ms1_and_higher_energy_spectra = spectra_ms1_and_higher_energy_matching()
    return (
        exp_ms1_spectra_mismatching,
        exp_ms1_and_lower_energy_spectra,
        exp_ms1_and_higher_energy_spectra,
    )


def test_splitting_mismatching_ms1_spectra_should_fail(
    spectra_mismatching_ms1: Tuple[MSExperiment],
    detected_ions_df_single_mod: pd.DataFrame,
    mods_single: List[str],
    scan_window_splitting: ScanWindowSplitting,
) -> None:
    with pytest.raises(AssertionError, match="MS1"):
        scan_window_splitting.split_windows_by_mods(
            *spectra_mismatching_ms1, detected_ions_df_single_mod, mods_single
        )


def test_splitting_wrong_mslevel_or_energy_spectra_should_fail(
    spectra_matching: Tuple[MSExperiment],
    detected_ions_df_single_mod: pd.DataFrame,
    mods_single: List[str],
    scan_window_splitting: ScanWindowSplitting,
) -> None:
    with pytest.raises(AssertionError, match="MSLevel"):
        scan_window_splitting.split_windows_by_mods(
            spectra_matching[1],
            spectra_matching[1],
            spectra_matching[2],
            detected_ions_df_single_mod,
            mods_single,
        )

    with pytest.raises(AssertionError, match="collision energy"):
        scan_window_splitting.split_windows_by_mods(
            spectra_matching[0],
            spectra_matching[1],
            spectra_matching[1],
            detected_ions_df_single_mod,
            mods_single,
        )

    with pytest.raises(AssertionError, match="collision energy"):
        scan_window_splitting.split_windows_by_mods(
            spectra_matching[0],
            spectra_matching[2],
            spectra_matching[2],
            detected_ions_df_single_mod,
            mods_single,
        )

    with pytest.raises(AssertionError, match="collision energy"):
        scan_window_splitting.split_windows_by_mods(
            spectra_matching[0],
            spectra_matching[2],
            spectra_matching[1],
            detected_ions_df_single_mod,
            mods_single,
        )


def test_splitting_no_mods(
    spectra_matching: Tuple[MSExperiment],
    mods_empty: List[str],
    scan_window_splitting: ScanWindowSplitting,
) -> None:
    windows_by_mod = scan_window_splitting.split_windows_by_mods(
        *spectra_matching, pd.DataFrame(), mods_empty
    )
    ms1_lower_energy_spectra = spectra_matching[1]
    assert windows_by_mod == {"unmodified": ms1_lower_energy_spectra}


@pytest.mark.parametrize(
    "spectra_matching_with_result_fixture, detected_ions_df_fixture, mods_to_search_fixture",
    [
        (
            "spectra_matching_with_result_single_mod",
            "detected_ions_df_single_mod",
            "mods_single",
        ),
        (
            "spectra_matching_with_result_single_mod",
            "detected_ions_df_multiple_mods",
            "mods_single",
        ),
        (
            "spectra_matching_with_result_multiple_mods",
            "detected_ions_df_multiple_mods",
            "mods_multiple",
        ),
        (
            "spectra_matching_with_result_multiple_mods",
            "detected_ions_df_multiple_mods_duplicate_entries",
            "mods_multiple",
        ),
    ],
)
def test_splitting_with_mods(
    spectra_matching_with_result_fixture: str,
    detected_ions_df_fixture: str,
    mods_to_search_fixture: str,
    request: SubRequest,
    scan_window_splitting: ScanWindowSplitting,
) -> None:
    spectra_matching_with_result = request.getfixturevalue(
        spectra_matching_with_result_fixture
    )
    detected_ions_df = request.getfixturevalue(detected_ions_df_fixture)
    mods_to_search = request.getfixturevalue(mods_to_search_fixture)

    windows_by_mod = scan_window_splitting.split_windows_by_mods(
        *spectra_matching_with_result[:-1], detected_ions_df, mods_to_search
    )
    expected_result = spectra_matching_with_result[-1]

    assert windows_by_mod == expected_result


@pytest.mark.parametrize(
    "spectra_matching_with_result_fixture, detected_ions_df_fixture, mods_to_search_fixture, mod_combinations_fixture",
    [
        (
            "spectra_matching_with_result_multiple_mods_combinations_and_single",
            "detected_ions_df_multiple_mods",
            "mods_multiple",
            "mod_combinations_some",
        ),
        (
            "spectra_matching_with_result_multiple_mods_combinations_and_single_reduced",
            "detected_ions_df_multiple_mods",
            "mods_single",
            "mod_combinations_some",
        ),
        (
            "spectra_matching_with_result_multiple_mods_combinations_and_single",
            "detected_ions_df_multiple_mods_duplicate_entries",
            "mods_multiple",
            "mod_combinations_some",
        ),
        (
            "spectra_matching_with_result_multiple_mods_combinations",
            "detected_ions_df_multiple_mods",
            "mods_multiple",
            "mod_combinations_all",
        ),
        (
            "spectra_matching_with_result_multiple_mods_combinations",
            "detected_ions_df_multiple_mods_duplicate_entries",
            "mods_multiple",
            "mod_combinations_all",
        ),
    ],
)
def test_splitting_with_mods_combinations(
    spectra_matching_with_result_fixture: str,
    detected_ions_df_fixture: str,
    mod_combinations_fixture: str,
    mods_to_search_fixture: str,
    request: SubRequest,
    scan_window_splitting: ScanWindowSplitting,
) -> None:
    spectra_matching_with_result = request.getfixturevalue(
        spectra_matching_with_result_fixture
    )
    detected_ions_df = request.getfixturevalue(detected_ions_df_fixture)
    mod_combinations = request.getfixturevalue(mod_combinations_fixture)
    mods_to_search = request.getfixturevalue(mods_to_search_fixture)

    windows_by_mod = scan_window_splitting.split_windows_by_mods(
        *spectra_matching_with_result[:-1],
        detected_ions_df,
        mods_to_search,
        mod_combinations_to_search=mod_combinations
    )
    expected_result = spectra_matching_with_result[-1]

    assert windows_by_mod == expected_result
