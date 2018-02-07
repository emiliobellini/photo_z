import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import pyccl as ccl
import matplotlib
matplotlib.use('Agg')
import pylab as plt



#Define ell_max.
#Ell_min is always 0, so it is easier to read the Cl's array.
#Ell_min is there just for plotting purposes.
ELL_MIN = 2
ELL_MAX = 2000
#Cosmology
Omega_c=0.27
Omega_b=0.045
h=0.67
A_s=1e-10
n_s=0.96


# Parse the given arguments
parser = argparse.ArgumentParser("Generate Cl's from binned Photo-z")
parser.add_argument("input_file", type=str, help="Input FITS file")
parser.add_argument("--output_plot", "-p", type=str, default = None, help="Output plot")
parser.add_argument("--skip_plots", "-s", help="Skip plots", action="store_true")
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
keys = sorted(keys)
pz = {}
for key in keys:
    pz[key] = table[key].data
print 'Read data from ' + os.path.relpath(paths['input'])
sys.stdout.flush()


#Calculate cosmology
cosmo = ccl.Cosmology(Omega_c=Omega_c, Omega_b=Omega_b, h=h, A_s=A_s, n_s=n_s)
print 'Calculated cosmology'
sys.stdout.flush()

#Calculate tracers
lens = {}
for key in keys:
    lens[key] = ccl.ClTracerLensing(cosmo, False, z=z, n=pz[key])
print 'Calculated tracers'
sys.stdout.flush()

#Calculate Cl's
ell = np.arange(0, ELL_MAX+1)
angular_cl = np.zeros((len(keys), len(keys), len(ell)))
for count1 in range(len(keys)):
    for count2 in range(len(keys)):
        lens1 = lens[keys[count1]]
        lens2 = lens[keys[count2]]
        angular_cl[count1][count2] = ccl.angular_cl(cosmo, lens1, lens2, ell)
print 'Calculated Cl\'s'
sys.stdout.flush()


#Update input file with a table containing Cl's
#Create table
columns = []
for count in range(len(keys)):
    l = str(ELL_MAX+1)
    columns.append(fits.Column(name=keys[count],array=angular_cl[count],format=l+'E'))
hdu = fits.BinTableHDU.from_columns(columns, name='Cls')
#Update existing file
with fits.open(paths['input'], mode='update') as hdul:
    hdul.append(hdu)
    hdul.flush()
print 'Updated input file at ' + os.path.relpath(paths['input'])
sys.stdout.flush()


#Plots
if not args.skip_plots:
    #Plot diagonal Cl's
    for count in range(len(keys)):
        plt.plot(ell, angular_cl[count][count], label=r'$'+keys[count]+'$')
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C^\ell_{ii}$')
    plt.title(r'Auto-correlations')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(ELL_MIN, ELL_MAX)
    plt.legend()
    plot_path = os.path.splitext(paths['output_plot'])
    plot_path = plot_path[0]+'_cls_ii'+plot_path[1]
    plt.savefig(plot_path)
    plt.close()
    print 'Saved auto-correlation plot at ' + os.path.relpath(plot_path)
    sys.stdout.flush()

    #Plot off-diagonal Cl's (a plot for each row of the matrix)
    for count1 in range(len(keys)):
        for count2 in range(len(keys)):
            plt.plot(ell, angular_cl[count1][count2], label=r'$'+keys[count2]+'$')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$C^\ell_{ij}$')
        plt.title(r'Cross-correlation with $'+keys[count1]+'$')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(ELL_MIN, ELL_MAX)
        plt.legend()
        plot_path = os.path.splitext(paths['output_plot'])
        plot_path = plot_path[0]+'_cls_ij_'+keys[count1]+plot_path[1]
        plt.savefig(plot_path)
        plt.close()
        print 'Saved cross-correlation plot for bin ' + keys[count1] + ' at ' + os.path.relpath(plot_path)
        sys.stdout.flush()

print 'Success!!'
