#!/bin/bash

for energy in 30 35 40 45 50 55 60 :
do
python ion_extraction_tryout_snrs_10ppm.py ../proteomics_data/230228_Immonium_${energy}NCE_mod.mzML --higher_collision_energy=${energy}
done