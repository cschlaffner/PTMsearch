import re
from pathlib import Path
from typing import Dict, FrozenSet, List, Union

import pandas as pd
from pyopenms import MSExperiment

from src.diagnostic_ions.utils import modifications_db, residue_db


def validate_number_lower_higher_energy_windows(
    lower_energy_windows_exp: MSExperiment, higher_energy_windows_exp: MSExperiment
) -> None:
    num_lower = len(lower_energy_windows_exp.getSpectra())
    num_higher = len(higher_energy_windows_exp.getSpectra())
    assert (
        num_lower == num_higher
    ), f"Numbers of higher- and lower-energy windows are not equal: {num_lower} VS {num_higher}."


def validate_modifications(modifications: Union[List[str], FrozenSet[str]]) -> None:
    assert not "unmodified" in modifications, (
        "Searching for a split without PTMs is included automatically, 'unmodified' "
        "should not be specified explicitly."
    )

    for mod in modifications:
        amino_acid = mod[0]
        unimod_id = re.findall("[0-9]+", mod)[0]
        assert (
            mod == f"{amino_acid}(UniMod:{unimod_id})"
        ), f"PTM {mod} has a wrong format."
        assert residue_db.hasResidue(amino_acid), f"PTM {mod} has an invalid residue."
        assert modifications_db.has(
            f"UniMod:{unimod_id}"
        ), f"PTM {mod} not found in UniMod database."


def validate_modification_combinations(
    modification_combinations: List[FrozenSet[str]],
) -> None:
    for combination in modification_combinations:
        assert not len(combination) > 1, (
            f"Combination {combination} contains only a single PTM and "
            "therefore has to be listed in modifications_to_search."
        )
        validate_modifications(combination)


def validate_path_exists(path: Path) -> None:
    assert path.exists(), f"Path {path} does not exist."


def validate_filepath_exists(path: Path) -> None:
    validate_path_exists(path)
    assert not path.is_dir(), f"Path {path} must be a file but is a directory."


def validate_spectral_library_files_by_mod(
    spectral_library_files_by_mod: Dict[Union[str, FrozenSet], str],
    modifications: List[str],
    modification_combinations: List[FrozenSet[str]],
) -> None:
    modifications_and_combinations = set(
        modifications + modification_combinations + ["unmodified"]
    )
    assert (
        "unmodified" in spectral_library_files_by_mod.keys()
    ), "A spectral library file for the split without PTMs has to be supplied."
    assert spectral_library_files_by_mod.keys() == modifications_and_combinations, (
        "The PTMs and combinations for which spectral library files are supplied "
        "do not match the PTMs and combinations to search for."
    )

    for library_path in spectral_library_files_by_mod.values():
        validate_filepath_exists(Path(library_path))


def validate_spectral_library_not_empty(library_df: pd.DataFrame, mod: str) -> None:
    assert len(library_df) > 0, (
        f"Filtered library for {mod} is empty. "
        f"Exclude {mod} from the selection or fix the library."
    )
