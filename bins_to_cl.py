import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import pyccl as ccl
import matplotlib
matplotlib.use('Agg')
import pylab as plt



#Define l_min and l_max
L_MIN = 2
L_MAX = 2000


# Parse the given arguments
parser = argparse.ArgumentParser("Generate Cl's from binned Photo-z")
parser.add_argument("input_file", type=str, help="Input FITS file")
parser.add_argument("--output_plot", "-p", type=str, default = None, help="Output plot")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)
if args.output_plot:
    paths['output_plot'] = os.path.abspath(args.output_plot)
else:
    paths['output_plot'] = os.path.splitext(paths['input'])[0]
    paths['output_plot'] = paths['output_plot'] + '.pdf'


#Open files and store datas
#Read table
table = Table.read(paths['input'], hdu='Photo-z')


#Extract binned photo-z
z = table['z'].data
keys = table.colnames
del keys[0]
pz = {}
for key in keys:
    pz[key] = table[key].data


#Calculate cosmology
cosmo = ccl.Cosmology(Omega_c=0.27, Omega_b=0.045, h=0.67, A_s=1e-10, n_s=0.96)

#Calculate tracers
lens = {}
for key in keys:
    lens[key] = ccl.ClTracerLensing(cosmo, False, z=z, n=pz[key])

#Calculate Cl's
ell = np.arange(L_MIN, L_MAX)
angular_cl = {}
for key1 in keys:
    angular_cl[key1] = {}
    for key2 in keys:
        angular_cl[key1][key2] = ccl.angular_cl(cosmo, lens[key1], lens[key2], ell)


#Plot diagonal Cl's
for key in keys:
    plt.plot(ell, angular_cl[key][key], label=r'$'+key+'$')
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C^\ell_{ii}$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(L_MIN, L_MAX)
plt.legend()
plot_path = os.path.splitext(paths['output_plot'])
plt.savefig(plot_path[0]+'_cls_ii'+plot_path[1])
plt.close()

#Plot off-diagonal Cl's (a plot for each row of the matrix)
for key1 in keys:
    for key2 in keys:
        plt.plot(ell, angular_cl[key1][key2], label=r'$'+key2+'$')
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C^\ell_{ij}$')
    plt.title(r'$'+key+'$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(L_MIN, L_MAX)
    plt.legend()
    plot_path = os.path.splitext(paths['output_plot'])
    plt.savefig(plot_path[0]+'_cls_ij_'+key1+plot_path[1])
    plt.close()
