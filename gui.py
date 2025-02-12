import logging
import tkinter as tk
from pathlib import Path
from tkinter import font as tkFont
from typing import Any, Dict, FrozenSet, List

from src.config.config import Config
from src.diagnostic_ions.detection import MassToleranceUnit
from src.main import main


def get_ptm_unimod(ptm_input: str) -> str:
    amino_acid = ptm_input[0]
    unimod_id = ptm_input[1:]
    return f"{amino_acid}(UniMod:{unimod_id})"


def get_ptm_list(ptm_list_input: str) -> List[str]:
    if ptm_list_input == "":
        return []
    ptm_list = ptm_list_input.split(";")
    return [get_ptm_unimod(mod) for mod in ptm_list]


def get_ptm_combinations_list(ptm_combinations_list_input: str) -> List[FrozenSet[str]]:
    if ptm_combinations_list_input == "":
        return []
    ptm_combinations_list = ptm_combinations_list_input.split("|")
    ptm_combinations_sets_list = [
        frozenset(get_ptm_list(combination)) for combination in ptm_combinations_list
    ]
    return ptm_combinations_sets_list


def get_dictionary(dict_input: str) -> Dict[str, Any]:
    if dict_input == "":
        return []

    def try_to_number(value):
        try:
            value_float = float(value)
            if value_float % 1 == 0:
                return int(value_float)
            return value_float
        except (ValueError, TypeError):
            return value

    dict_list = dict_input.split(";")
    key_value_list = [key_value.split(":") for key_value in dict_list]
    dict_config = {
        key_value[0]: try_to_number(key_value[1]) for key_value in key_value_list
    }
    return dict_config


# TODO: handle splits for PTM combinations


def create_config_json() -> Config:
    config = Config(
        result_dir=result_dir.get(),
        mzml_file=mzml_file.get(),
        dia_nn_path=dia_nn_path.get(),
        library_free=library_free.get(),
        modifications_to_search=get_ptm_list(modifications_to_search.get()),
        modification_combinations=get_ptm_combinations_list(
            modification_combinations.get()
        ),
        modifications_additional=get_ptm_list(modifications_additional.get()),
        lower_collision_energy=float(lower_collision_energy.get()),
        higher_collision_energy=float(higher_collision_energy.get()),
        diagnostic_ions_mass_tolerance=float(diagnostic_ions_mass_tolerance.get()),
        diagnostic_ions_mass_tolerance_unit=diagnostic_ions_mass_tolerance_unit.get(),
        snr_threshold=float(snr_threshold.get()),
        fdr_threshold=float(fdr_threshold.get()),
        known_diagnostic_ions_file=known_diagnostic_ions_file.get(),
        spectral_library_files_by_mod=get_dictionary(
            spectral_library_files_by_mod.get()
        ),
        spectral_library_for_filtering_path=spectral_library_for_filtering_path.get(),
        database_for_library_prediction=database_for_library_prediction.get(),
        dia_nn_library_params=get_dictionary(dia_nn_library_params.get()),
        dia_nn_search_params=get_dictionary(dia_nn_search_params.get()),
        detection_count_percentile=float(detection_count_percentile.get()),
        detection_count_min=int(detection_count_min.get()),
        normalize_cscores=normalize_cscores.get(),
        run_only_ion_detection=run_only_ion_detection.get(),
    )
    config_path = Path(config_dir.get()) / f"{experiment_name.get()}.json"
    config.save(config_path)
    return config


def create_config_and_run() -> None:
    config = create_config_json()

    logger = logging.getLogger(__name__)
    logging.basicConfig(
        filename=Path(config.result_dir) / "log.txt", level=logging.INFO
    )

    try:
        main(config, logger)
    except BaseException as e:
        logger.exception(e)
        raise e


root = tk.Tk(screenName=None, baseName=None, className="Tk", useTk=1)
root.title("PTMsearch")
root.minsize(500, 700)

arial = tkFont.Font(family="Arial", size=12)
arial_headline = tkFont.Font(family="Arial", size=12, weight="bold")

# -------------------------- Filepaths -------------------
experiment_name = tk.StringVar()
experiment_name_entry = tk.Entry(root, textvariable=experiment_name)
experiment_name_entry.grid(row=0, column=1)
tk.Label(root, text="Experiment name", font=arial).grid(row=0, sticky=tk.W)

config_dir = tk.StringVar()
config_dir_entry = tk.Entry(root, textvariable=config_dir)
config_dir_entry.grid(row=1, column=1)
tk.Label(root, text="Config folder path", font=arial).grid(row=1, sticky=tk.W)

result_dir = tk.StringVar()
result_dir_entry = tk.Entry(root, textvariable=result_dir)
result_dir_entry.grid(row=2, column=1)
tk.Label(root, text="Result folder path", font=arial).grid(row=2, sticky=tk.W)

mzml_file = tk.StringVar()
mzml_file_entry = tk.Entry(root, textvariable=mzml_file)
mzml_file_entry.grid(row=3, column=1)
tk.Label(root, text="MzML input file path", font=arial).grid(row=3, sticky=tk.W)

dia_nn_path = tk.StringVar()
dia_nn_path_entry = tk.Entry(root, textvariable=dia_nn_path)
dia_nn_path_entry.grid(row=4, column=1)
tk.Label(root, text="DIA-NN executable path", font=arial).grid(row=4, sticky=tk.W)

known_diagnostic_ions_file = tk.StringVar()
known_diagnostic_ions_file_entry = tk.Entry(
    root, textvariable=known_diagnostic_ions_file
)
known_diagnostic_ions_file_entry.grid(row=5, column=1)
known_diagnostic_ions_file_entry.insert(
    0, "src/diagnostic_ions/known_ions_only_unimod.csv"
)
tk.Label(root, text="Diagnostic ions data path", font=arial).grid(row=5, sticky=tk.W)

# -------------------------- Ion detection -------------------

row_ion_detection = 6
tk.Label(root, text="Ion detection config", font=arial_headline).grid(
    row=row_ion_detection, sticky=tk.W
)

lower_collision_energy = tk.StringVar()
lower_collision_energy_entry = tk.Entry(root, textvariable=lower_collision_energy)
lower_collision_energy_entry.grid(row=row_ion_detection + 1, column=1)
tk.Label(root, text="Lower NCE value", font=arial).grid(
    row=row_ion_detection + 1, sticky=tk.W
)

higher_collision_energy = tk.StringVar()
higher_collision_energy_entry = tk.Entry(root, textvariable=higher_collision_energy)
higher_collision_energy_entry.grid(row=row_ion_detection + 2, column=1)
tk.Label(root, text="Higher NCE value", font=arial).grid(
    row=row_ion_detection + 2, sticky=tk.W
)

diagnostic_ions_mass_tolerance = tk.StringVar()
diagnostic_ions_mass_tolerance_entry = tk.Entry(
    root, textvariable=diagnostic_ions_mass_tolerance
)
diagnostic_ions_mass_tolerance_entry.grid(row=row_ion_detection + 3, column=1)
tk.Label(root, text="Ion detection mass tolerance", font=arial).grid(
    row=row_ion_detection + 3, sticky=tk.W
)

diagnostic_ions_mass_tolerance_unit = tk.StringVar(None, MassToleranceUnit.ppm.name)
tk.Radiobutton(
    root,
    text="ppm",
    variable=diagnostic_ions_mass_tolerance_unit,
    value=MassToleranceUnit.ppm.name,
    font=arial,
).grid(row=row_ion_detection + 4, sticky=tk.W)
tk.Radiobutton(
    root,
    text="Dalton",
    variable=diagnostic_ions_mass_tolerance_unit,
    value=MassToleranceUnit.Da.name,
    font=arial,
).grid(row=row_ion_detection + 5, sticky=tk.W)

snr_threshold = tk.StringVar()
snr_threshold_entry = tk.Entry(root, textvariable=snr_threshold)
snr_threshold_entry.grid(row=row_ion_detection + 6, column=1)
tk.Label(root, text="Ion detection SNR threshold", font=arial).grid(
    row=row_ion_detection + 6, sticky=tk.W
)

run_only_ion_detection = tk.BooleanVar()
tk.Checkbutton(
    root, text="Run only ion detection", variable=run_only_ion_detection, font=arial
).grid(row=row_ion_detection + 7, sticky=tk.W)

# -------------------------- PTM selection -------------------
row_ptm_selection = 15
tk.Label(root, text="PTM selection config", font=arial_headline).grid(
    row=row_ptm_selection, sticky=tk.W
)

tk.Label(
    root,
    text="PTMs to search alone",
    font=arial,
).grid(row=row_ptm_selection + 1, sticky=tk.W)
modifications_to_search = tk.StringVar()
modifications_to_search_entry = tk.Entry(root, textvariable=modifications_to_search)
modifications_to_search_entry.grid(row=row_ptm_selection + 1, column=1)

tk.Label(
    root,
    text="PTMs to search in combination",
    font=arial,
).grid(row=row_ptm_selection + 2, sticky=tk.W)
modification_combinations = tk.StringVar()
modification_combinations_entry = tk.Entry(root, textvariable=modification_combinations)
modification_combinations_entry.grid(row=row_ptm_selection + 2, column=1)

tk.Label(
    root,
    text="Additional PTMs to search",
    font=arial,
).grid(row=row_ptm_selection + 3, sticky=tk.W)
modifications_additional = tk.StringVar()
modifications_additional_entry = tk.Entry(root, textvariable=modifications_additional)
modifications_additional_entry.grid(row=row_ptm_selection + 3, column=1)

tk.Label(
    root,
    text="Automatic selection: percentile",
    font=arial,
).grid(row=row_ptm_selection + 4, sticky=tk.W)
detection_count_percentile = tk.StringVar()
detection_count_percentile_entry = tk.Entry(
    root, textvariable=detection_count_percentile
)
detection_count_percentile_entry.grid(row=row_ptm_selection + 4, column=1)
detection_count_percentile_entry.insert(0, 1.0)

tk.Label(
    root,
    text="Automatic selection: minimum number of windows",
    font=arial,
).grid(row=row_ptm_selection + 5, sticky=tk.W)
detection_count_min = tk.StringVar()
detection_count_min_entry = tk.Entry(root, textvariable=detection_count_min)
detection_count_min_entry.grid(row=row_ptm_selection + 5, column=1)
detection_count_min_entry.insert(0, 0)

# -------------------------- Library handling -------------------

row_library_handling = 21
tk.Label(root, text="Spectral library config", font=arial_headline).grid(
    row=row_library_handling, sticky=tk.W
)

tk.Label(
    root,
    text="Spectral library paths by PTMs",
    font=arial,
).grid(row=row_library_handling + 1, sticky=tk.W)
spectral_library_files_by_mod = tk.StringVar()
spectral_library_files_by_mod_entry = tk.Entry(
    root, textvariable=spectral_library_files_by_mod
)
spectral_library_files_by_mod_entry.grid(row=row_library_handling + 1, column=1)

tk.Label(
    root,
    text="Filtering: spectral library path",
    font=arial,
).grid(row=row_library_handling + 2, sticky=tk.W)
spectral_library_for_filtering_path = tk.StringVar()
spectral_library_for_filtering_path_entry = tk.Entry(
    root, textvariable=spectral_library_for_filtering_path
)
spectral_library_for_filtering_path_entry.grid(row=row_library_handling + 2, column=1)

library_free = tk.BooleanVar()
tk.Checkbutton(root, text="Library-free mode", variable=library_free, font=arial).grid(
    row=row_library_handling + 3, sticky=tk.W
)

tk.Label(
    root,
    text="Library-free: database path",
    font=arial,
).grid(row=row_library_handling + 4, sticky=tk.W)
database_for_library_prediction = tk.StringVar()
database_for_library_prediction_entry = tk.Entry(
    root, textvariable=database_for_library_prediction
)
database_for_library_prediction_entry.grid(row=row_library_handling + 4, column=1)

tk.Label(
    root,
    text="Library-free: DIA-NN prediction params",
    font=arial,
).grid(row=row_library_handling + 5, sticky=tk.W)
dia_nn_library_params = tk.StringVar()
dia_nn_library_params_entry = tk.Entry(root, textvariable=dia_nn_library_params)
dia_nn_library_params_entry.grid(row=row_library_handling + 5, column=1)


# -------------------------- DIA-NN search and aggregation -------------------

row_search_aggregation = 27
tk.Label(root, text="DIA-NN search and aggregation config", font=arial_headline).grid(
    row=row_search_aggregation, sticky=tk.W
)

tk.Label(
    root,
    text="DIA-NN search params",
    font=arial,
).grid(row=row_search_aggregation + 1, sticky=tk.W)
dia_nn_search_params = tk.StringVar()
dia_nn_search_params_entry = tk.Entry(root, textvariable=dia_nn_search_params)
dia_nn_search_params_entry.grid(row=row_search_aggregation + 1, column=1)

tk.Label(
    root,
    text="FDR threshold",
    font=arial,
).grid(row=row_search_aggregation + 2, sticky=tk.W)
fdr_threshold = tk.StringVar()
fdr_threshold_entry = tk.Entry(root, textvariable=fdr_threshold)
fdr_threshold_entry.grid(row=row_search_aggregation + 2, column=1)
fdr_threshold_entry.insert(0, 0.01)

normalize_cscores = tk.BooleanVar()
tk.Checkbutton(
    root, text="Normalize c-scores", variable=normalize_cscores, font=arial
).grid(row=row_search_aggregation + 3, sticky=tk.W)

"""
Buttons: generate config JSON, generate config JSON and run"""

# -------------------------- Buttons -------------------
row_buttons = 31
button_config = tk.Button(
    root, text="Create config JSON", width=25, command=create_config_json, font=arial
)
button_config.grid(row=row_buttons, column=0)
button_config_and_run = tk.Button(
    root,
    text="Create config JSON and run",
    width=25,
    command=create_config_and_run,
    font=arial,
)
button_config_and_run.grid(row=row_buttons, column=1)

root.mainloop()
