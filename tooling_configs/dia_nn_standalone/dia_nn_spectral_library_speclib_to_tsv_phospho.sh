#!/bin/bash

/hpi/fs00/home/andrea.nathansen/diann/diann-1.8.1 --gen-spec-lib \
--threads 8 \
--lib ptm-search-data/data/library_filtered_database/in_silico_library_filtered_phospho.predicted.speclib \
--out-lib ptm-search-data/data/library_filtered_database/in_silico_library_filtered_phospho.tsv
