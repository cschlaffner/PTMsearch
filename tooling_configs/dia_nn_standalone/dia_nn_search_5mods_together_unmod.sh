#!/bin/bash

/hpi/fs00/home/andrea.nathansen/diann/diann-1.8.1 \
--f ptm-search-data/data/230928_JL_Immonium_ions_Unmodified_DIA_lower_energy.mzML \
--lib ptm-search-data/data/library_filtered_database/in_silico_library_filtered_5mods.predicted.speclib \
--out ptm-search-data/results_thesis/dia_nn_results_filtered/report_filtered_5mods_together_unmod.tsv \
--mass-acc 10 --mass-acc-ms1 20 --window 0 --threads 8 --pg-level 2 \
--var-mod UniMod:1,42.010565,K \
--var-mod UniMod:34,14.01565,K \
--var-mod UniMod:35,15.994915,P \
--var-mod UniMod:354,44.985078,Y \
--var-mod UniMod:21,79.966331,Y 