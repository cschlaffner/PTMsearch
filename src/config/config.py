import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, FrozenSet, List, Union


@dataclass
class Config:
    # TODO: set some defaults(?)
    result_dir: str
    # TODO: allow multiple files to allow multi-run
    mzml_file: str
    dia_nn_path: str
    # TODO: change later (this setting is WIP)
    spectral_library_files_by_mod: Dict[str, str]
    spectral_library_for_filtering_path: str
    library_free: bool
    database_for_library_prediction: str
    modifications_to_search: List[str]
    modification_combinations: List[FrozenSet[str]]
    modifications_additions: List[str]
    lower_collision_energy: Union[float, int]
    higher_collision_energy: Union[float, int]
    known_diagnostic_ions_file: str
    diagnostic_ions_mass_tolerance: Union[float, int]
    diagnostic_ions_mass_tolerance_unit: str
    snr_threshold: Union[float, int]

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
