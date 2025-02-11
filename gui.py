import tkinter as tk
from tkinter import font as tkFont

from src.diagnostic_ions.detection import MassToleranceUnit

root = tk.Tk(screenName=None, baseName=None, className="Tk", useTk=1)
root.title("PTMsearch")
root.minsize(600, 500)

arial = tkFont.Font(family="Arial", size=12)
arial_headline = tkFont.Font(family="Arial", size=12, weight="bold")

# -------------------------- Filepaths -------------------
experiment_name = tk.StringVar()
experiment_name_entry = tk.Entry(root, textvariable=experiment_name)
experiment_name_entry.grid(row=0, column=1)
tk.Label(root, text="Experiment name", font=arial).grid(row=0, sticky=tk.W)

result_dir = tk.StringVar()
result_dir_entry = tk.Entry(root, textvariable=result_dir)
result_dir_entry.grid(row=1, column=1)
tk.Label(root, text="Result folder path", font=arial).grid(row=1, sticky=tk.W)

mzml_path = tk.StringVar()
mzml_path_entry = tk.Entry(root, textvariable=mzml_path)
mzml_path_entry.grid(row=2, column=1)
tk.Label(root, text="MzML input file path", font=arial).grid(row=2, sticky=tk.W)

dia_nn_path = tk.StringVar()
dia_nn_path_entry = tk.Entry(root, textvariable=dia_nn_path)
dia_nn_path_entry.grid(row=3, column=1)
tk.Label(root, text="DIA-NN executable path", font=arial).grid(row=3, sticky=tk.W)

known_diagnostic_ions_file = tk.StringVar()
known_diagnostic_ions_file_entry = tk.Entry(
    root, textvariable=known_diagnostic_ions_file
)
known_diagnostic_ions_file_entry.grid(row=4, column=1)
known_diagnostic_ions_file_entry.insert(
    0, "src/diagnostic_ions/known_ions_only_unimod.csv"
)
tk.Label(root, text="Diagnostic ions data path", font=arial).grid(row=4, sticky=tk.W)

# -------------------------- Ion detection -------------------

row_ion_detection = 5
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
row_ptm_selection = 14
tk.Label(root, text="PTM selection config", font=arial_headline).grid(
    row=row_ptm_selection, sticky=tk.W
)

tk.Label(
    root,
    text="PTMs to search alone (specify as semicolon-separated list,",
    font=arial,
).grid(row=row_ptm_selection + 1, sticky=tk.W)
tk.Label(
    root,
    text="leave empty for automatic PTM selection)",
    font=arial,
).grid(row=row_ptm_selection + 2, sticky=tk.W)
modifications_to_search = tk.StringVar()
modifications_to_search_entry = tk.Entry(root, textvariable=modifications_to_search)
modifications_to_search_entry.grid(row=row_ptm_selection + 2, column=1)

tk.Label(
    root,
    text="PTMs to search in combination (specify as vertical-bar-",
    font=arial,
).grid(row=row_ptm_selection + 3, sticky=tk.W)
tk.Label(
    root,
    text="separated list of semicolon-separated lists)",
    font=arial,
).grid(row=row_ptm_selection + 4, sticky=tk.W)
modification_combinations = tk.StringVar()
modification_combinations_entry = tk.Entry(root, textvariable=modification_combinations)
modification_combinations_entry.grid(row=row_ptm_selection + 4, column=1)

tk.Label(
    root,
    text="Additional PTMs to search (specify as semicolon-separated list)",
    font=arial,
).grid(row=row_ptm_selection + 5, sticky=tk.W)
modifications_additional = tk.StringVar()
modifications_additional_entry = tk.Entry(root, textvariable=modifications_additional)
modifications_additional_entry.grid(row=row_ptm_selection + 5, column=1)

tk.Label(
    root,
    text="Automatic selection: percentile",
    font=arial,
).grid(row=row_ptm_selection + 6, sticky=tk.W)
detection_count_percentile = tk.StringVar()
detection_count_percentile_entry = tk.Entry(
    root, textvariable=detection_count_percentile
)
detection_count_percentile_entry.grid(row=row_ptm_selection + 6, column=1)
detection_count_percentile_entry.insert(0, 1.0)

tk.Label(
    root,
    text="Automatic selection: minimum number of windows",
    font=arial,
).grid(row=row_ptm_selection + 7, sticky=tk.W)
detection_count_min = tk.StringVar()
detection_count_min_entry = tk.Entry(root, textvariable=detection_count_min)
detection_count_min_entry.grid(row=row_ptm_selection + 7, column=1)
detection_count_min_entry.insert(0, 0)

# -------------------------- Library handling -------------------

row_library_handling = 22
library_free = tk.BooleanVar()
tk.Checkbutton(root, text="Library-free mode", variable=library_free, font=arial).grid(
    row=row_library_handling, sticky=tk.W
)
""" - spectral_library_files_by_mod (dict)
 - spectral_library_for_filtering_path (path)
 - database_for_library_prediction (path)
  - dia_nn_library_params (dict)"""

# -------------------------- DIA-NN search and aggregation -------------------

row_search_aggregation = 23
"""
 - dia_nn_search_params (dict)
   - fdr_threshold (number)"""

normalize_cscores = tk.BooleanVar()
tk.Checkbutton(
    root, text="Normalize c-scores", variable=normalize_cscores, font=arial
).grid(row=row_search_aggregation + 1, sticky=tk.W)

"""
Buttons: generate config JSON, generate config JSON and run"""

root.mainloop()
