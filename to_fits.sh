#!/bin/bash

/usr/local/shared/python/2.7.6-gcc/bin/python P1_tsv_to_fits.py ../data/data_cfhtlens.tsv
echo "Finished P1_tsv_to_fits"

/usr/local/shared/python/2.7.6-gcc/bin/python P2_unify_fits.py -d ../data/data_cfhtlens__
echo "Finished P2_unify_fits"

#addqueue -c '3hours' -q cmb -m 40 -n 1 ./to_tsv.sh
