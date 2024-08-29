import re
from typing import Dict, List

import numpy as np
import pandas as pd
import pyopenms as oms
from pyopenms import MSSpectrum

from src.diagnostic_ions.utils import diff_mono_mass_to_unimod_format
from src.mzml_processing.utils import check_ms1_spectrum


# expecting library to be tsv
# only take those with the mods for the windows with mods, do a run with unmodified for all windows
def split_library_by_mods(
    spectrum_library: pd.DataFrame, has_unimod_format: bool
) -> Dict[str, pd.DataFrame]:
    library_entry_lists_by_mods = {
        "unmodified": [pd.DataFrame(columns=spectrum_library.columns)]
    }
    for lib_entry in spectrum_library.itertuples():
        # explicit string setting to please typecheck
        sequence = str(spectrum_library.transition_group_id)
        if has_unimod_format:
            mods = np.unique(re.findall(".(UniMod:[0-9]+)", sequence))
        else:
            mods = np.unique(re.findall(".\[[0-9]+.[0-9]+\])", sequence))
            # TODO: check
            mods = [diff_mono_mass_to_unimod_format(mod[0], float(mod[2:-1]))]

        if len(mods) == 0:
            library_entry_lists_by_mods["unmodified"].append(lib_entry)
            continue

        # TODO: check how to handle those multiple-mod cases
        for mod in mods:
            library_entry_lists_by_mods[mod].append(lib_entry)

    libraries_by_mod = {}
    for mod, library_entry_list in library_entry_lists_by_mods.items():
        libraries_by_mod[mod] = pd.concat(library_entry_list, ignore_index=True)

    # TODO: all unmodified to all with mods?

    return libraries_by_mod
    # maybe make faster: take only the wanted mods?
    # later: use pyopenms.ModificationsDB() and amino acid stuff to unify found ions and library

    # afterwards: compare list of keys (except unmodified) by intersection, get all that are in library but not in windows (give
    # log note or so), get all that are in windows but not in library, match matching and write SL and windows and start DIA-NN process
