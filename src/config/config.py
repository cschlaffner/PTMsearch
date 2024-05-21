import json
from dataclasses import dataclass
from pathlib import Path
from typing import Union


@dataclass
class Config:
    lower_collision_energy: Union[float, int]
    higher_collision_energy: Union[float, int]
    known_diagnostic_ions_file: str
    diagnostic_ions_mass_tolerance: Union[float, int]
    diagnostic_ions_mass_tolerance_unit: str
    intensity_threshold: Union[float, int]

    @classmethod
    def from_path(cls, config_path: Path):
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)
        return cls(**config_json)
