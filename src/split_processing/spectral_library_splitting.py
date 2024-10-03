import re
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from src.diagnostic_ions.utils import diff_mono_mass_to_unimod_format


# expecting library to be tsv
def split_library_by_mods(
    spectrum_library: pd.DataFrame,
    has_unimod_format: bool,
    mod_combinations_to_search: List[Tuple[str]],
) -> Dict[str, pd.DataFrame]:
    library_entry_lists_by_mods: Dict[str, List] = {}
    unimod_regex = re.compile(".(UniMod:[0-9]+)")
    mass_diff_regex = re.compile(r".\[[0-9]+.[0-9]+\]")

    for lib_entry in spectrum_library.itertuples(index=False):
        # explicit string setting to please typecheck
        sequence = str(lib_entry.ModifiedPeptide)
        if has_unimod_format:
            mods = np.unique(re.findall(unimod_regex, sequence))
        else:
            mods = np.unique(re.findall(mass_diff_regex, sequence))
            # TODO: check
            mods = [
                diff_mono_mass_to_unimod_format(mod[0], float(mod[2:-1]))
                for mod in mods
            ]

        lib_entry_df = pd.DataFrame(data=[lib_entry])
        if len(mods) == 0:
            mods = ["unmodified"]

        # Single-mod searches
        # TODO: how to handle it? When some mods only appear in combination with others?
        # for mod in mods:
        if len(mods) == 1:
            mod = mods[0]
            if mod not in library_entry_lists_by_mods:
                library_entry_lists_by_mods[mod] = [lib_entry_df]
            else:
                library_entry_lists_by_mods[mod].append(lib_entry_df)

        # Mod combination searches
        for mod_combination in mod_combinations_to_search:
            if set(mods).issubset(set(mod_combination)):
                if mod_combination not in library_entry_lists_by_mods:
                    library_entry_lists_by_mods[mod_combination] = [lib_entry_df]
                else:
                    library_entry_lists_by_mods[mod_combination].append(lib_entry_df)

    # Add unmodified entries to all other libraries
    for mod in library_entry_lists_by_mods:
        if mod == "unmodified":
            continue
        library_entry_lists_by_mods[mod] = np.concatenate(
            [
                library_entry_lists_by_mods[mod],
                library_entry_lists_by_mods["unmodified"],
            ]
        )

    print("Done filtering, concatenating...")
    libraries_by_mod = {}
    for mod, library_entry_list in library_entry_lists_by_mods.items():
        libraries_by_mod[mod] = pd.concat(library_entry_list, ignore_index=True)

    return libraries_by_mod
