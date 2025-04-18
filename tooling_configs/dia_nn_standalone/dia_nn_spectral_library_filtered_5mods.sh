#!/bin/bash


/hpi/fs00/home/andrea.nathansen/diann/diann-1.8.1 --gen-spec-lib --fasta-search --predictor \
--threads 8 \
--cut K*,R* \
--min-pep-len 7 \
--max-pep-len 30 \
--max-pr-charge 4 \
--min-pr-mz 350 --max-pr-mz 1200 --min-fr-mz 200 --max-fr-mz 2000 \
--missed-cleavages 2 \
--var-mod UniMod:1,42.010565,K \
--var-mod UniMod:34,14.01565,K \
--var-mod UniMod:35,15.994915,P \
--var-mod UniMod:354,44.985078,Y \
--var-mod UniMod:21,79.966331,Y \
--strip-unknown-mods \
--var-mods 3 \
--pg-level 2 \
--out-lib ptm-search-data/data/library_filtered_database/in_silico_library_filtered_5mods.predicted \
--fasta ptm-search-data/data/database_for_diann/uniprotkb_AND_reviewed_true_AND_model_o_2024_07_16_with_crap_and_synthetic_filtered.fasta \
