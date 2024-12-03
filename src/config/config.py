import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Union


@dataclass
class Config:
    # TODO: add documentation
    result_dir: str
    """The directory in which all intermediate and end results will be stored."""

    mzml_file: str
    """The file that contains higher- and lower-energy windows.
    Currently, only one file (= one run) is supported.
    Must be compatible with DIA-NN, see
    https://github.com/vdemichev/DiaNN?tab=readme-ov-file#raw-data-formats"""

    dia_nn_path: str
    """Path to the DIA-NN executable. Currently, this code is tested with DIA-NN 1.8.1,
    so there might be issues with other versions."""

    library_free: bool
    """In library-free mode, the spectral libraries per split will be predicted by
    DIA-NN."""

    modifications_to_search: List[str]
    modification_combinations: List[FrozenSet[str]]
    modifications_additions: List[str]

    lower_collision_energy: Union[float, int]
    """Collision energy of the fragment ion windows, i.e., the MS2 spectra
    that should be used for the DIA-NN search."""

    higher_collision_energy: Union[float, int]
    """Collision energy of the immonium ion windows, i.e., the MS2 spectra
    that should be used for the ion/modification detection."""

    diagnostic_ions_mass_tolerance: Union[float, int]
    """Mass tolerance (numeric value only) for ion detection and DIA-NN search."""

    diagnostic_ions_mass_tolerance_unit: str
    """Unit of the mass tolerance (supported: ppm or Da)."""

    snr_threshold: Union[float, int]
    """Signal-to-noise ratio threshold for ion detection: immonium ion spectra with a lower SNR will
    be discarded. SNR of a spectrum is calculated as peak_intensity_max / peak_intensity_mean."""

    fdr_threshold: float

    known_diagnostic_ions_file: str = "src/diagnostic_ions/known_ions_only_unimod.csv"
    """Path to a CSV file containing information about diagnostic ions for the detection.
    Current one contains selected ions from (TODO: cite) that are listed in UniMod.
    If you provide your own ions file, make sure that it contains only ions for modifications that
    are listed in UniMod. Otherwise, it's not compatible with prediction and search with DIA-NN."""

    # TODO: change this to support also combinations
    spectral_library_files_by_mod: Dict[Union[str, FrozenSet[str]], str] = field(
        default_factory=lambda: {}
    )
    """If you already have spectral libraries: The spectral library
    for each modification. Should contain precursors for the modification
    and also the unmodified ones. Those files will not be validated/checked
    for correct modifications."""

    spectral_library_for_filtering_path: str = ""
    """If you already have one spectral library containing precursors
    for all modifications (and potentially combinations): this library
    will be filtered and split into one library per split (determined by the
    modifications/ions). Precursors containing mods that should not be searched
    for are discarded."""

    database_for_library_prediction: str = ""
    """If library-free mode is used, you must specify a database that will
    be used for the spectral library prediction by DIA-NN."""

    detection_count_percentile: float = 1.0
    """If modifications and combinations should be detected automatically: Consider only the
    combinations that contain all together (=in sum) this fraction of the total number
    of windows with modifications (combinations are sorted beforehand by number of windows 
    descending). This means that lowering this threshold will lead
    to less selected combinations. Applies only to combination selection, not to single
    modification selection."""

    detection_count_min: int = 0
    """If modifications and combinations should be detected automatically: Consider
    only modifications and combinations that occur in at least this many windows.
    For single modifications, the occurences in combinations that are selected will
    be excluded when calculating the number of occurences for the single modification.
    """

    normalize_cscores: bool = False
    """If set to true, min-max normalization to 0-1 range will be run for the c-scores
    of each split before aggregation, to avoid exceptionally large c-scores for some splits
    (which can occur if DIA-NN falls back to a linear classifier for c-score calculation
    due to low number of identifications)."""

    run_only_ion_detection: bool = False
    """If this is set to True, the workflow stops after ion detection
    and does not perform any search. Use this if you want to get an overview
    over the immonium ions first and want to do the search in another run."""

    @classmethod
    def from_path(cls, config_path: Path):
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)
            mods_combinations = [
                frozenset(combination)
                for combination in config_json["modification_combinations"]
            ]
            config_json["modification_combinations"] = mods_combinations
        return cls(**config_json)
