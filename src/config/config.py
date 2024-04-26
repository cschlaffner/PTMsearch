import json
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Config:
    lower_collision_energy: int
    higher_collision_energy: int

    @classmethod
    def from_path(cls, config_path: Path):
        with open(config_path, "r") as config_file:
            config_json = json.load(config_file)
        return cls(**config_json)
