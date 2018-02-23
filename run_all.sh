#!/bin/bash

RED='\033[0;31m'
NC='\033[0m' # No Color

/usr/local/shared/python/2.7.6-gcc/bin/python 1_tsv_to_fits.py ../data/data_cfhtlens.tsv
echo -e ${RED}"Finished 1_tsv_to_fits"${NC}

/usr/local/shared/python/2.7.6-gcc/bin/python 2_unify_fits.py -d ../data/data_cfhtlens__
echo -e ${RED}"Finished 2_unify_fits"${NC}

/usr/local/shared/python/2.7.6-gcc/bin/python 3_binned_photo-z.py ../data/data_cfhtlens.fits
echo -e ${RED}"Finished 3_binned_photo-z"${NC}

/usr/local/shared/python/2.7.6-gcc/bin/python 4_bins_to_cl.py -s ../data/data_cfhtlens_bins.fits
echo -e ${RED}"Finished 4_bins_to_cl"${NC}

#addqueue -c '5hours' -q cmb -m 40 -n 1 ./run_all.sh
