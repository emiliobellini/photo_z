#!/bin/bash

/usr/local/shared/python/2.7.6-gcc/bin/python 1_tsv_to_fits.py ../data/data_cfhtlens.tsv
echo "Finished 1_tsv_to_fits"

/usr/local/shared/python/2.7.6-gcc/bin/python 2_unify_fits.py -d ../data/data_cfhtlens__
echo "Finished 2_unify_fits"

/usr/local/shared/python/2.7.6-gcc/bin/python 3_binned_data.py ../data/data_cfhtlens.fits
echo "Finished 3_binned_photo-z"

/usr/local/shared/python/2.7.6-gcc/bin/python 4_bins_to_cl.py ../data/data_cfhtlens_bins.fits
echo "Finished 4_bins_to_cl"

/usr/local/shared/python/2.7.6-gcc/bin/python 5_kl_transform.py ../data/data_cfhtlens_bins.fits
echo "Finished 5_kl_transform"

#addqueue -c '5hours' -q cmb -m 40 -n 1 ./run_all.sh
