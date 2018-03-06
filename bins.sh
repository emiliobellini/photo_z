#!/bin/bash

/usr/local/shared/python/2.7.6-gcc/bin/python 1_binned_data.py ../data/data_cfhtlens.fits
echo "Finished 1_binned_photo-z"

/usr/local/shared/python/2.7.6-gcc/bin/python 2_bins_to_cl.py ../data/data_cfhtlens_bins.fits
echo "Finished 2_bins_to_cl"

/usr/local/shared/python/2.7.6-gcc/bin/python 3_kl_transform.py ../data/data_cfhtlens_bins.fits
echo "Finished 3_kl_transform"

#addqueue -c '3hours' -q cmb -m 40 -n 1 ./bins.sh
