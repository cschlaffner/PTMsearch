import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Union


@dataclass
class Config:
    tmp_dir: str
    # TODO: allow multiple files to allow multi-run
    mzml_file: str
    dia_nn_path: str
    # TODO: change later (this setting is WIP)
    spectral_library_files_by_mod: Dict[str, str]
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
        return cls(**config_json)
