import logging
import re
from typing import Dict, FrozenSet, List, Optional, Tuple

import numpy as np
import pandas as pd

from src.diagnostic_ions.utils import (
    diff_mono_mass_to_unimod_format,
    get_mod_combination_str,
)


def split_library_by_mods(
    spectral_library: pd.DataFrame,
    mods_to_search: List[str],
    mod_combinations_to_search: Optional[List[FrozenSet[str]]] = None,
    additional_mods_to_search: Optional[List[str]] = None,
    has_unimod_format: bool = True,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, pd.DataFrame]:
    """Splits the spectral library according to PTMs, taking into account
    the combinations. Each single PTM and combination gets its own library.
    For a combination, all precursors that carry a subset of the set of
    PTMs are sorted into the corresponding library. Additional mods are
    ignored during search, i.e., the respective precursors are sorted into
    the libraries they would be sorted in without the additional mods.
    Precursors that carry PTMs that were not specified are discarded.
    Each resulting library contains all unmodified precursors.
    The input library is expected to be in TSV format and must have DIA-NN
    compatible columns. The PTMs within the library must be either given in
    UniMod format (e.g Y(UniMod:21)) or in mass difference format
    (e.g. Y[79.9663]). Mass diff format has not been tested with
    an actual library and for compatibility with the rest of the
    software, so it is currently not used in the main workflow."""

    if mod_combinations_to_search is None:
        mod_combinations_to_search = []
    if additional_mods_to_search is None:
        additional_mods_to_search = []

    library_entry_lists_by_mods: Dict[str, List] = {}
    unimod_regex = re.compile(r".\(UniMod:[0-9]+\)")
    mass_diff_regex = re.compile(r".\[[0-9]+.[0-9]+\]")

    for lib_entry in spectral_library.itertuples(index=False):
        # explicit string setting to please typecheck
        sequence = str(lib_entry.ModifiedPeptide)
        if has_unimod_format:
            mods = np.unique(re.findall(unimod_regex, sequence))
        else:
            mods = np.unique(re.findall(mass_diff_regex, sequence))
            mods = [
                diff_mono_mass_to_unimod_format(mod[0], float(mod[2:-1]))
                for mod in mods
            ]

        # take out additional mods that should be searched but not taken
        # into account for splitting
        if additional_mods_to_search != []:
            mods = [mod for mod in mods if mod not in additional_mods_to_search]

        if len(mods) == 0:
            mods = ["unmodified"]

        # Single-mod searches
        if len(mods) == 1:
            mod = mods[0]
            if mod in mods_to_search or mod == "unmodified":
                if mod not in library_entry_lists_by_mods:
                    library_entry_lists_by_mods[mod] = [lib_entry]
                else:
                    library_entry_lists_by_mods[mod].append(lib_entry)

        # Mod combination searches
        if mod_combinations_to_search is not None:
            for mod_combination in mod_combinations_to_search:
                if set(mods).issubset(mod_combination):
                    if mod_combination not in library_entry_lists_by_mods:
                        library_entry_lists_by_mods[mod_combination] = [lib_entry]
                    else:
                        library_entry_lists_by_mods[mod_combination].append(lib_entry)

    # If no precursors for a PTM or combination were detected: create empty libraries
    for mod in mods_to_search + mod_combinations_to_search + ["unmodified"]:
        if mod not in library_entry_lists_by_mods:
            library_entry_lists_by_mods[mod] = []
        if logger is not None:
            logger.info("%s: %s precursors", mod, len(library_entry_lists_by_mods[mod]))

    # Add unmodified entries to all other libraries
    for mod in library_entry_lists_by_mods:
        if mod == "unmodified":
            continue
        library_entry_lists_by_mods[mod] = (
            library_entry_lists_by_mods[mod] + library_entry_lists_by_mods["unmodified"]
        )

    libraries_by_mod = {}
    for mod, library_entry_list in library_entry_lists_by_mods.items():
        library_df = pd.DataFrame(library_entry_list)
        libraries_by_mod[mod] = library_df

        if logger is None:
            continue

        if len(library_df) == 0:
            logger.warning(
                f"No precursors for {mod} found in the library!"
                f" The split {mod} will not yield any results."
            )
            continue

        if mod == "unmodified":
            continue

        mod_indication_string = "UniMod" if has_unimod_format else "["
        num_mods_in_lib = (
            library_df["ModifiedPeptide"].str.contains(mod_indication_string).sum()
        )
        if num_mods_in_lib == 0:
            logger.warning(
                f"No modified precursors for split {get_mod_combination_str(mods)} were found in"
                " the library. Search for this split will only be run with unmodified precursors."
                " Consider adapting your spectral library."
            )

    return libraries_by_mod
