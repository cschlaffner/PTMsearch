# PTMsearch

This software leverages immonium ion information from stepped fragmentation DIA data for pre-filtering PTMs for peptide identification search. The higher-energy windows of the stepped fragmentation data are extracted and used for immonium ion detection. The information about which immonium ions were detected in which spectra is then used for splitting the corresponding lower-energy spectra accordingly with subsequent split-wise search by [DIA-NN](https://github.com/vdemichev/DiaNN) ([Demichev et al. (2020)](https://doi.org/10.1038/s41592-019-0638-x)). The split-wise DIA-NN reports are then aggregated into a single report.

## Getting started

The software uses Python 3.11. A conda environment is included in the `environment.yml` file and can be set up with `conda env create -f environment.yml`.

For running in command-line mode, create a JSON configuration for your experiment and run `python -m src.main <path/to/input/config>`.

To use the GUI, run `python gui.py` in the terminal and add the config parameters in the GUI. A config JSON will be automatically created when clicking the respective button.

![](https://github.com/cschlaffner/PTMsearch/blob/main/ScreenshotGUI.png)

## Configuration

In the GUI, the following config parameters have to be provided (lists should be specified without whitespaces):
- **Experiment name**: Name of the experiment that will be used for creating the config.
- **Config folder path**: Path where the config JSON is stored.
- **Result folder path**: Path where intermediate and end results will be stored. If the folder does not exist yet, it will be created.
- **MzML input file path**: An mzML file containing the stepped fragmentation DIA data. Must contain an equal number of higher-and lower-energy windows with matching precursor m/z ranges, so that for each lower-energy window, a higher-energy window can be matched by precursor m/z and spectrum ID of the preceding MS1 window. The mzML must be compatible for DIA-NN (see [DIA-NN documentation](https://github.com/vdemichev/DiaNN/tree/719be81544c70888f65a34a07f0643ae1be59570?tab=readme-ov-file#raw-data-formats)). Currently, only a single run (i.e., a single mzML file) is supported.
- **DIA-NN executable path**: Path to DIA-NN executable that is used for split-wise search and, if library-free mode is selected, also for library prediction. The software has been tested with [DIA-NN 1.8.1](https://github.com/vdemichev/DiaNN/tree/719be81544c70888f65a34a07f0643ae1be59570), other versions might lead to errors.
- **Diagnostic ions data path**: Path to a CSV file containing data about the immonium ions. The default, in `src/diagnostic_ions/known_ions_only_unimod.csv`, is a subset of the immonium ion data reported by [Hung et al. (2007)](https://doi.org/10.1007/s00216-007-1449-y) and [Zolg et al. (2018)](https://doi.org/10.1074/mcp.tir118.000783). It is advisable to use this default. If you want to search for ions that are not included in the default file, you can extend the file or use your own (which has to have the same format as the default). Only PTMs that are listed in UniMod are supported.
- **Ion detection config**:
    - **Lower NCE value**: NCE of the MS2 spectra containing the fragment ions.
    - **Higher NCE value**: NCE of the MS2 spectra containing the immonium ions.
    - **Ion detection mass tolerance**: Mass tolerance that is applied during ion detection, can be in ppm or Dalton.
    - **Ion detection SNR threshold**: Signal-to-noise ratio (spectrum intensity maximum divided by average) that is applied to discard noisy spectra, i.e., spectra with a lower SNR, during detection.
    - **Run only ion detection**: If selected, the software stops after ion detection and does not run subsequent searches.
- **PTM selection config**:
    - **PTMs to search alone**: PTMs that will be considered during spectral library handling, window splitting and searching. Will be searched alone, i.e., in a separate split per PTM. Specify as semicolon-separated list, with PTMs written as one-letter code and UniMod ID, e.g. `K1;Y21` to search for K acetylation and Y phosphorylation. If left blank, PTMs and combinations will be selected automatically from the detected diagnostic ions.
    - **PTM combinations to search**: PTM combinations to consider during spectral library handling, window splitting and searching. E.g. if searching for the combination of K acetylation, K methylation, and Y phosphorylation, all windows that contain exactly those three ions are added to the split for this combination, and will be searched with a library containing precursors that carry any single or multiple of those PTMs. A PTM can be listed as to search alone and as part of combinations, where searching in combinations takes precedence and the PTM is searched alone in all remaining windows in which its ion is detected. E.g., if K acetylation is listed as part of combinations, each combination will be searched in a separate split of all windows containing exactly that combination, and all remaining windows where K acetylation is detected are searched in another split for K acetylation. Specify as semicolon-separated list of combinations, in each combination the PTMs have to be separated by vertical bars, e.g. `K1|K34|Y21;K1|Y21` to search for the combination of K acetylation, K methylation, and Y phosphorylation, and the combination of K acetylation and Y phosphorylation. The number of PTMs allowed per combination is limited by DIA-NN, for DIA-NN 1.8.1 it was 5. The PTM combinations will be selected automatically if the field for PTMs to search alone is left blank.
    - **Additional PTMs to search**: PTMs that are searched in every split regardless of detected immonium ions, to include PTMs that rarely or never emit immonium ions. Specification format is the same as for PTMs to search alone. They count towards the maximum number of PTMs allowed per combination.
    - **Automatic selection: minimum number of windows**: If the field PTMs to search alone is left blank, PTMs and combinations will be selected automatically from the detected ions. Only the PTMs and combinations detected in at least this many windows are selected.
 - **Spectral library config**: Either split-wise libraries (1st field), a library for filtering (2nd field) or library-free mode with a sequence database (3rd and following fields) must be specified.
    - **Spectral library paths by PTMs**: A user-provided spectral library for each PTM (and possibly combination) to search for, listed by PTM (and possibly combination). In this case, the spectral library contents are your own responsibility regarding containing the respective PTMs (and unmodified precursors to search in the splits). The libraries (and PTM format) must be compatible with DIA-NN (see [DIA-NN documentation](https://github.com/vdemichev/DiaNN/tree/719be81544c70888f65a34a07f0643ae1be59570?tab=readme-ov-file#spectral-library-formats)). The PTMs and combinations must be exactly the same as the PTMs in the PTM selection config, with the addition of a library for "unmodified" containing only unmodified precursors for the split without PTMs. Specify as semicolon-separated list of entries `split:path/to/library`, e.g. `K1:path/example1.speclib;K1|Y21:path/example2.speclib`.
    - **Filtering: spectral library path**: A single library to be filtered and segmented into split-wise libraries for the respective PTMs. It must be a TSV file containing PTMs in the same format as the ones in DIA-NN-predicted libraries, e.g. "Y(UniMod:21)", and must be compatible with DIA-NN.
    - **Library-free mode**: If selected, the spectral library for each split is predicted by DIA-NN from a sequence database. The caveat: The DIA-NN predictor [is only trained for some PTMs](https://github.com/vdemichev/DiaNN?tab=readme-ov-file#creation-of-spectral-libraries), so prediction quality for other PTMs is questionable.
    - **Library-free: database path**: Path to the sequence database for spectral library prediction in library-free mode.
    - **Library-free: DIA-NN prediction params**: Additional parameters for DIA-NN for library prediction. Specify as semicolon-separated list of entries `parameter:value` (or `parameter:` for parameters without value), e.g., `cut:K*,R*;threads:8`. The mandatory parameters (`gen-spec-lib`, `fasta-search`, `predictor`, `out-lib`, `fasta`, `strip-unknown-mods`) are already added automatically.
-  **DIA-NN search and aggregation config**:
    - **DIA-NN search params**: Additional parameters for DIA-NN search. Specify as semicolon-separated list of entries `parameter:value` (or `parameter:` for parameters without value), e.g., `window:0;no-lib-filter:;threads:8`. The mandatory parameters (`f`, `lib`, `out`) and further parameters (`qvalue` (set to 1 because aggregated q-value calculation is conducted afterwards), `decoy-report`, `no-prot-inf`) are already added automatically.
    - **FDR threshold**: The FDR at which the aggregated output report should be filtered.
    - **Normalize c-scores**: If selected, c-scores are normalized to the 0-1 range for each split to migitate artifacts that are a result of the search with DIA-NN being conducted split-wise.

For manually creating the config JSON for command-line mode, have a look at `src/config/config.py` and `gui.py` for the corresponding names and formatting of the configurable parameters.


## Output

The software provides an aggregated report with all target and decoys precursors unfiltered (`report_aggregated_all_targets_with_decoys.csv`) and a report with targets filtered by the configured FDR threshold (`report_aggregated_fdr_filtered.csv`).
The report includes all columns from the split-wise DIA-NN reports. For decoys, a lot of columns are empty and some are shifted because the decoys output from DIA-NN does not assign all columns correctly, which is fixed for the most important columns in the aggregated report. Further, the protein information was calculated split-wise by DIA-NN. This software supports only precursor identification and no subsequent protein inference or quantification, thus, most protein-related columns might not be useful.

The following columns are added in the aggregated report:
- `q_value_aggregated`: Q-value computed after aggregation. This should be used instead of the split-wise calculated `Q.Value` from the DIA-NN report which is incorrect after aggregation.
- `remapped_id`: The renamed identifier of the spectrum in which the precursor was found. Because the spectra identifiers per split have to be renamed before search as DIA-NN requires continuous numbers.
- `original_id`: The corresponding spectra identifiers as they were in the original provided mzML file before renaming.
- `CScore_normalized`: This column is only added if c-score normalization is enabled. Then, the c-scores within each split are normalized to the 0-1 range to avoid skewed distributions if some splits lead to large c-score values. The normalized c-scores are then used for q-value calculation and FDR filtering, the `CScore` column calculated by DIA-NN should be ignored in this case.

Currently, the software stores all intermediate results (mzML files and spectrum id mappings, predicted or filtered libraries, DIA-NN reports for splits) in the specified result directory. This can require a larger amount of storage space, so you might want to delete the intermediate results afterwards.
