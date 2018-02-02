import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import pyccl as ccl
import pylab as plt
import math



# Parse the given arguments
parser = argparse.ArgumentParser("Generate Cl's from binned Photo-z")
parser.add_argument("input_file", type=str, help="Input FITS file")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)


#Open files and store datas
#Read table
table = Table.read(paths['input'], hdu='Photo-z')

print table
# cosmo = ccl.Cosmology(Omega_c=0.27, Omega_b=0.045, h=0.67, A_s=2.1e-9, n_s=0.96,
#                       transfer_function='boltzmann_class')
#
# print ccl.core.transfer_function_types.keys()
#
# k = np.logspace(-4., 1., 100) # Wavenumber
# a = 1. # Scale factor
#
# pk_lin = ccl.linear_matter_power(cosmo, k, a)
#
# pk_nl = ccl.nonlin_matter_power(cosmo, k, a)
#
# plt.plot(k, pk_lin, 'b-')
# plt.plot(k, pk_nl, 'r-')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$k$',fontsize=22)
# plt.ylabel(r'$P(k)$',fontsize=22)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.show()
