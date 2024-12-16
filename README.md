# PTMsearch

- software for my master thesis
- summary of algo
- has been tested on Linux. Although I think it should probably also run on Windows if the file paths are adapted accordingly, but that has not been tested yet.
- currently only command line tool
- currently not very optimized, e.g. DIA-NN runs are sequentially and not in parallel.
- Currently only for precursor/peptide identification. Protein inference after aggregation not included.

- input:
    - stepped fragmentation data: dia-nn compatible (link, command line command) stepped fragmentation mzML file. Currently, only a single run ("single mzML file) is supported.
    - library: spectral library for each PTM to search for (own responsibility, must match specified PTMs and combinations), or single library to be searched/filtered (must contain PTMs in UniMod or mass diff format, must be TSV and compatible for DIA-NN), or just set library-free mode -> will predict libs split-wise automatic with DIA-NN. Last option caveat: DIA-NN predictor only trained for some mods, rest prediction ambiguous/take with grain of salt. (DIA-NN doc VS my results)
    - PTMs to search for: single PTMs and combinations, possible additional PTMs (searched for in every split). PTMs and Combs can also be selected automatically based on detected ions.
    - A DIA-NN executable that is used for search and, if library-free mode is selected, also for library prediction. The software has been tested with DIA-NN 1.8.1, other versions might lead to errors.
    - All input paths, library settings etc. are provided via a JSON file. Config.py lists all configurable parameters with their documentation. Configs for all experiments of my thesis are in folder.

- Run software
    - command with config file
    

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
    - tooling_configs: includes all software configs that I used for my experiments. Also includes all DIA-NN scripts in an extra folder.
