import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Union


@dataclass
class Config:
    """Config for the runs, to be specified in a JSON file."""

    result_dir: str
    """The directory in which all intermediate and end results will be stored."""

    mzml_file: str
    """The file that contains higher- and lower-energy windows. Assumes equal
    number of higher-and lower-energy windows, i.e., that each lower-energy window
    has exactly one corresponding higher-energy window. Assumes same number of MS2
    windows for each MS1 window.
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
    """PTMs that should be considered during window splitting and search.
    Must be in DIA-NN-compatible UniMod format, e.g. Y(UniMod:21). If specified,
    diagnostic ions found for other PTMs will be ignored. If not specified, PTMs and
    combinations will be selected automatically based on the detected ions. In that
    case, you should configure detection_count_percentile and detection_count_min for
    automatic detection."""

    modification_combinations: List[FrozenSet[str]]
    """PTMs that should be considered together, i.e., extra splits will be
    created for windows that contain those exact combinations. Will also be taken
    into account for library filtering or prediction for the respective splits.
    Specification in the config file should be in list form, e.g.
    [["K(UniMod:1)","P(UniMod:35)"], ["K(UniMod:1)", "P(UniMod:35)"]]
    The PTMs mentioned here must also be specified in modifications_to_search.
    If no PTMs are specified in modifications_to_search, PTMs and combinations
    will be selected automatically."""

    modifications_additions: List[str]
    """Additional PTMs that do not (or rarely) emit immonium ions but are wanted
    to be included in the search. They will be added to the search for each split
    and also to each library prediction. If library-based mode with a single library
    is used, they will be ignored during library filtering (i.e., precursors with
    those PTMs will not be discarded, but added to the splits based on the
    other PTMs they carry)."""

    lower_collision_energy: Union[float, int]
    """Collision energy of the fragment ion windows, i.e., the MS2 spectra
    that should be used for the DIA-NN search."""

    higher_collision_energy: Union[float, int]
    """Collision energy of the immonium ion windows, i.e., the MS2 spectra
    that should be used for the ion/PTM detection."""

    diagnostic_ions_mass_tolerance: Union[float, int]
    """Mass tolerance (numeric value only) for ion detection and DIA-NN search."""

    diagnostic_ions_mass_tolerance_unit: str
    """Unit of the mass tolerance (supported: ppm or Da)."""

    snr_threshold: Union[float, int]
    """Signal-to-noise ratio threshold for ion detection: immonium ion spectra with a lower SNR will
    be discarded. SNR of a spectrum is calculated as peak_intensity_max / peak_intensity_mean."""

    fdr_threshold: float = 0.01
    """The FDR/q-value at which the aggregated output report should be filtered.
    Default is 1%."""

    known_diagnostic_ions_file: str = "src/diagnostic_ions/known_ions_only_unimod.csv"
    """Path to a CSV file containing information about diagnostic ions for the detection.
    Current one contains selected ions from (TODO: cite) that are listed in UniMod.
    If you provide your own ions file, make sure that it contains only ions for PTMs that
    are listed in UniMod. Otherwise, it's not compatible with prediction and search with DIA-NN."""

    spectral_library_files_by_mod: Dict[Union[str, FrozenSet[str]], str] = field(
        default_factory=lambda: {}
    )
    """If you already have spectral libraries: The spectral library
    for each PTM. Should contain precursors for the PTM
    and also the unmodified ones. Those files will not be validated/checked
    for correct PTMs."""

    spectral_library_for_filtering_path: str = ""
    """If you already have one spectral library containing precursors
    for all PTMs (and potentially combinations): this library
    will be filtered and split into one library per split (determined by the
    PTMs/ions). Precursors containing mods that should not be searched
    for are discarded."""

    database_for_library_prediction: str = ""
    """If library-free mode is used, you must specify a database that will
    be used for the spectral library prediction by DIA-NN."""

    dia_nn_library_params: Dict[str, Union[str, int]] = field(
        default_factory=lambda: {}
    )
    """Additional parameters that should be added for the library prediction
    with DIA-NN (if library-free mode is used), such as cut, min-pr-mz etc.
    The mandatory parameters (gen-spec-lib, fasta-search, predictor, out-lib,
    fasta, strip-unknown-mods), additional parameters (pg-level(automatically
    set to 2)) and the respective PTMs are already added automatically."""

    dia_nn_search_params: Dict[str, Union[str, int]] = field(default_factory=lambda: {})
    """Additional parameters that should be added for the DIA-NN search, such as mass-acc,
    window etc. Ideally, mass-acc should be set to the same value as the tolerance for
    diagnostic ion detection (although DIA-NN always treats the value as ppm, not Da).
    The mandatory parameters (f (the input file), lib, out) further parameters
    (qvalue (set to 1 because aggregated q-value calculation is conducted afterwards),
    decoy-report, no-prot-inf, pg-level (automatically set to 2))
    and the respective PTMs are already added automatically."""

    detection_count_percentile: float = 1.0
    """If PTMs and combinations should be detected automatically: Consider only the
    combinations that contain all together (=in sum) this fraction of the total number
    of windows with PTMs (combinations are sorted beforehand by number of windows 
    descending). This means that lowering this threshold will lead
    to less selected combinations. Applies only to combination selection, not to single
    PTM selection. Thus, you should also specify a detection_count_min to avoid small
    splits for single PTMs."""

    detection_count_min: int = 0
    """If PTMs and combinations should be detected automatically: Consider
    only PTMs and combinations that occur in at least this many windows.
    For single PTMs, the occurences in combinations that are selected will
    be excluded when calculating the number of occurences for the single PTM.
    """

    normalize_cscores: bool = False
    """If set to true, min-max normalization to 0-1 range will be run for the c-scores
    of each split before aggregation, to avoid exceptionally large c-scores for some splits."""

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
