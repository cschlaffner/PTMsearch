#!/bin/bash

/hpi/fs00/home/andrea.nathansen/diann/diann-1.8.1 \
--f /hpi/fs00/home/andrea.nathansen/proteomics_data/230928_JL_Immonium_ions_Modified_DIA_lower_energy.mzML \
--lib /hpi/fs00/home/andrea.nathansen/proteomics_data/library_filtered_database/in_silico_library_filtered_acetyl.predicted.speclib \
--out /hpi/fs00/home/andrea.nathansen/dia_nn_results_filtered/report_filtered_acetyl_all_spectra.tsv \
--mass-acc 10 --mass-acc-ms1 20 --window 0 --threads 8 --pg-level 2 \
--var-mod UniMod:1,42.010565,K