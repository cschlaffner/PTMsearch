# PTMsearch

This repository contains the software for my master thesis "A Software for Dynamic Post-Translational Modification
Searching in DIA Proteomics Data". It leverages immonium ion from stepped fragmentation DIA data for pre-filtering PTMs for search. The higher-energy windows of the stepped fragmentation data are extracted and used for immonium ion detection. The information about which immonium ions were detected in which spectra is then used for splitting the corresponding lower-energy spectra accordingly with subsequent split-wise search by [DIA-NN](https://github.com/vdemichev/DiaNN) (Demichev et al. (2020)). The split-wise DIA-NN reports are then aggregated into a single report.

So far, this software has only been run on Linux. Although it should probably also run on Windows if the file paths are adapted accordingly, but that has not been tested yet. Some fuctionality has not been optimized yet, e.g., DIA-NN runs are conducted sequentially and not in parallel. Currently, only precursor identification is included, but no protein inference after aggregation. The list of immonium ions that is used during detection, stored in `src/diagnostic_ions/known_ions_only_unimod.csv`, is a subset of the immonium ion data reported by Hung et al. (2007) and Zolg et al. (2018). It can be edited or an own provided file in the same format can be used.

The software is run via command line with Python, tested for Python version 3.11.8. The required packages with versions are listed in `enviroment.yml`. After setup, the software can be run with the command 
```
python -m src.main <path/to/input/config>
```

## Input

All input paths, settings etc. have to be provided via a JSON file. `src/config/config.py` lists all configurable parameters with their documentation. Configs for all experiments of my thesis are in the folder `tooling_configs`. Following is an overview over important required inputs for the config:
- `mzml_file`: An mzML file containing the stepped fragmentation DIA data. Must contain an equal number higher-and lower-energy windows with matching precursor m/z ranges, so that for each lower-eneergy window, a higher-energy window can be matched by precursor m/z and spectrum ID of the preceding MS1 window. The mzML must be compatible for DIA-NN (see [DIA-NN documentation](https://github.com/vdemichev/DiaNN/tree/719be81544c70888f65a34a07f0643ae1be59570?tab=readme-ov-file#raw-data-formats)). Currently, only a single run (which means a single mzML file) is supported.
- Input for spectral library can be one of the following:
    - `database_for_library_prediction`: a sequence databse for split-wise automatic spectral librar prediction with DIA-NN. The caveat: The DIA-NN predictor [is only trained for some PTMs](https://github.com/vdemichev/DiaNN?tab=readme-ov-file#creation-of-spectral-libraries), so prediction quality for other PTMs is questioable. When using this option, `library_free` must be set to true in the config.
    - `spectral_library_files_by_mod`: a user-provided spectral library for each PTM (and possibly combination) to search for, listed by PTM (and possibly combination). In this case, the spectral library contents are your own responsibility regarding containing the respective PTMs. The libraries must be compatible with DIA-NN (see [DIA-NN documentation](https://github.com/vdemichev/DiaNN/tree/719be81544c70888f65a34a07f0643ae1be59570?tab=readme-ov-file#spectral-library-formats)). When using this option, `library_free` must be set to false.
    - `spectral_library_for_filtering_path`: a single library to be filtered and segmented into split-wise libraries for the respective PTMs. It must be a TSV file containing PTMs in the same format as the ones in DIA-NN-predicted libraries, e.g. "Y(UniMod:21)", and must be compatible with DIA-NN. When using this option, `library_free` must be set to false. Currently, this feature has been tested (aside from unit tests) for a very small reduced library (10 precursors). A large library with about 4.4M precursors led to an out-of-memory error. Since this feature was not needed for experiments so far, it was left it in an experimental state requiring optimization for upscaling.
 - Input for PTMs to search for:
     - single PTMs: `modifications_to_search`
     - PTM combinations: `modification_combinations`
     - possible additional PTMs that do not/rarely emit immonium ions and should be searched for in every split: `modifications_additional`
     - They must be specified in the same format as in DIA-NN-predicted libraries (e.g. "Y(UniMod:21)"). See example configs in the `tooling_configs` folder for how exactly the specification structure for single PTMs and combinations looks like. 
     - If you do not want to specify any of those, specify them as empty list. If `modifications_to_search` is empty, PTMs and combinations are selected automatically based on detected immonium ions. The automatic detection can be configured with `detection_count_percentile` and `detection_count_min`.
- `dia_nn_path`: Path to DIA-NN executable that is used for split-wise search and, if library-free mode is selected, also for library prediction. The software has been tested with [DIA-NN 1.8.1](https://github.com/vdemichev/DiaNN/tree/719be81544c70888f65a34a07f0643ae1be59570), other versions might lead to errors.
- All other required inputs are listed in `src/config/config.py`.
    

- Results
    - one report with all target and decoys precursors unfiltered (report_aggregated_all_targets_with_decoys.csv), one with targets filtered by FDR threshold (report_aggregated_fdr_filtered.csv).
    - report columns same as in DIA-NN report (ref) For decoys a lot of columns are empty and some are shifted because the decoys output from DIA-NN does not assign all columns correctly. Further, the protein information columns were added split-wise. Since protein inference after aggregation has not been implemented yet, most protein-related columns are probably useless.
    - additional columns: q_value_aggregated: q-value computed after aggregation (use this instead of Q.Value, that one is within the split and useless).
    remapped_id: the spectra identifiers per split have to be renamed because DIA-NN needs continuous input. The renamed identifier of the spectrum in which the precursor was found. (MS2.Scan from DIA-NN only gives the index of the MS2 spectra when listing only the MS2 spectra -> less useful)
    original_id: the corresponding spectra identifiers as they were in the original mzML file before renaming.
    CScore_normalized: this column only occurs if c-score normalization is enabled. Then, the c-scores within each split are normalized to the 0-1 range to avoid skewed distributions if some splits lead to resulting large c-score values. The normalized c-scores are then used for q-value calculation and FDR filtering, the CScore column by DIA-NN should be ignored then.
    - software currently stores all intermediate results (mzML files and spectrum id mappings, predicted or filtered libraries, DIA-NN reports (add naming)) in the specified result directory -> can require a larger amount of storage space, you might want to delete the intermediate results afterwards

- Folder structure
    - src: the actual software implementation: core functionality and overall workflow. Includes a few TODOs that could be solved for better user-friendliness (mostly additional validation of inputs).
    - test: unit tests for core functions
    - notebooks: which I used for data post-processing and plot creation. Also database completion (from UniProt, contaminated and synthetic) and database filtering for a reduced one. The notebooks were for interediate/additional usage during the thesis, therefore are not documented as well and code quality is not emphasized there. Also, some notebooks were used on my laptop and some on the server, so file paths can be inconsistent.
    - tooling_configs: includes all software configs that I used for my experiments. The two configs for single mods (_single_mods and _single_mods_all_predicted) lead to the same results, the difference is only that one uses spectral libraries that were predicted in standalone runs before, the other includes library prediction in the workflow. Also includes all DIA-NN scripts in an extra folder.
