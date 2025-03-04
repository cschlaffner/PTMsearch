#!/bin/bash

/hpi/fs00/home/andrea.nathansen/diann/diann-1.8.1 \
--f ptm-search-data/data/230928_JL_Immonium_ions_Unmodified_DIA_lower_energy.mzML \
--out ptm-search-data/results_thesis/dia_nn_results_sl_whole_database/report_whole_database_with_lib_gen_file_unmodifiedDIA_pglevel2_mzlimits_with_synthetic.tsv \
--mass-acc 10 --mass-acc-ms1 20 --window 0 --threads 8 --pg-level 2 \
--gen-spec-lib --fasta-search --predictor \
--cut K*,R* \
--min-pep-len 7 \
--max-pep-len 30 \
--max-pr-charge 4 \
--missed-cleavages 2 \
--min-pr-mz 350 --max-pr-mz 1200 --min-fr-mz 200 --max-fr-mz 2000 \
--fasta ptm-search-data/data/database_for_diann/uniprotkb_AND_reviewed_true_AND_model_o_2024_07_16_with_crap_and_synthetic.fasta