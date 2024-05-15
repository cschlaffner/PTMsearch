import json
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Config:
    lower_collision_energy: float
    higher_collision_energy: float
    known_ions_file: str

    @classmethod
    def from_path(cls, config_path: Path):
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)
        return cls(**config_json)
