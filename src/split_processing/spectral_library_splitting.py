import re
from typing import Dict, List

import numpy as np
import pandas as pd
from src.diagnostic_ions.utils import diff_mono_mass_to_unimod_format


# expecting library to be tsv
# only take those with the mods for the windows with mods, do a run with unmodified for all windows
def split_library_by_mods(
    spectrum_library: pd.DataFrame, has_unimod_format: bool
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

        # TODO: check how to handle those multiple-mod cases
        for mod in mods:
            if mod not in library_entry_lists_by_mods:
                library_entry_lists_by_mods[mod] = [lib_entry_df]
            else:
                library_entry_lists_by_mods[mod].append(lib_entry_df)
    print("Done filtering, concatenating...")
    libraries_by_mod = {}
    for mod, library_entry_list in library_entry_lists_by_mods.items():
        libraries_by_mod[mod] = pd.concat(library_entry_list, ignore_index=True)

    # TODO: all unmodified to all with mods?

    return libraries_by_mod
    # TODO: maybe make faster: take only the wanted mods?
