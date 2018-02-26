import os, sys, fnmatch
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from astropy.io import fits
from astropy.table import Table, vstack
import variables as vrs



# Parse the given arguments
parser = argparse.ArgumentParser("Generate binned Photo-z")
parser.add_argument("input_file", type=str, help="Input FITS file")
parser.add_argument("--output_file", "-o", type=str, default = None, help="Output file")
parser.add_argument("--output_plot", "-p", type=str, default = None, help="Output plot")
parser.add_argument("--skip_info", "-s", help="Skip info messages and save time", action="store_true")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)
if args.output_file:
    paths['output_file'] = os.path.abspath(args.output_file)
else:
    paths['output_file'] = os.path.splitext(paths['input'])[0]
    paths['output_file'] = paths['output_file'] + '_bins.fits'
if args.output_plot:
    paths['output_plot'] = os.path.abspath(args.output_plot)
else:
    paths['output_plot'] = os.path.splitext(paths['output_file'])[0]
    paths['output_plot'] = paths['output_plot'] + '.pdf'


#Open files and store datas
#Read table
table = Table.read(paths['input'])
#Read image
image = fits.open(paths['input'])['PZ_full'].data


#Print effective densities for each field
#TODO: insert correct values for A_eff
if not args.skip_info:
    sel_bins = vrs.select_data(table, vrs.Z_MIN, vrs.Z_MAX)
    w, w2 = {'TOT': 0.}, {'TOT': 0.}
    n_eff = {}
    n_gal = {'TOT': 0}
    A_eff = {'W1': 42.9, 'W2': 12.1, 'W3': 26.1, 'W4': 13.3}
    A_eff['TOT'] = np.array(A_eff.values()).sum()
    for key in ['W1', 'W2', 'W3', 'W4']:
        w[key] = np.array([x['weight'] for x in table[sel_bins] if key in x['field']])
        n_gal[key] = len(w[key])
        w2[key] = (w[key]**2).sum()
        w[key] = w[key].sum()
        n_eff[key] = w[key]**2/w2[key]/(A_eff[key]*60.*60.)
        n_gal['TOT'] = n_gal['TOT'] + n_gal[key]
        w['TOT'] = w['TOT'] + w[key]
        w2['TOT'] = w2['TOT'] + w2[key]
        print('Selected ' + str(n_gal[key]) + ' galaxies for field ' + key
            + ' (' + '{:2.4f}'.format(n_eff[key]) + ' per sq arcmin).')
        sys.stdout.flush()
    n_eff['TOT'] = w['TOT']**2/w2['TOT']/(A_eff['TOT']*60.*60.)
    print('Selected in total ' + str(n_gal['TOT']) + ' galaxies ('
        + '{:2.4f}'.format(n_eff['TOT']) + ' per sq arcmin).')
    sys.stdout.flush()


#Get arrays with selected data for each bin
sel_bins = {}
for key in vrs.Z_BINS.keys():
    sel_bins[key] = vrs.select_data(table, vrs.Z_BINS[key][0], vrs.Z_BINS[key][1])


#Generate binned Photo-z

#Step size (provided by CFHTlens to generate the Photo-z arrays)
z_step = 0.05
photoz_x = (np.arange(len(image[0]))+1./2.)*z_step

#Normalize photoz
image = (image.T / image.sum(axis=1)).T / z_step
print 'Normalized Photo-z'
sys.stdout.flush()

#Calculate binned Photo-z
photoz_y = {}
for key in vrs.Z_BINS.keys():
    w_tot = table['weight'][sel_bins[key]].sum()
    photoz_y[key] = np.dot(table['weight'][sel_bins[key]].data, image[sel_bins[key]])/w_tot
    print 'Calculated Photo-z for bin ' + key
    sys.stdout.flush()


#Save Photo-z data

#Create the first empty image
hdu = fits.PrimaryHDU()
hdul = fits.HDUList([hdu])

#Create a table with all the galaxies for each bin
del hdu
for key in sorted(vrs.Z_BINS.keys()):
    hdu_part = table[sel_bins[key]]
    new_col = np.full(len(hdu_part),key)
    hdu_part['BIN'] = new_col
    try:
        hdu = vstack([hdu, hdu_part], join_type='exact')
    except NameError:
        hdu = hdu_part
hdu = fits.table_to_hdu(hdu)
hdul.append(hdu)
print 'Created table with data'
sys.stdout.flush()

#Create image with redshift bins
hdu = fits.ImageHDU([vrs.Z_BINS[x] for x in sorted(vrs.Z_BINS.keys())], name='Z_bins')
hdul.append(hdu)

#Create table with binned Photo-z
columns = []
columns.append(fits.Column(name='z',array=photoz_x,format='E'))
for key in sorted(photoz_y.keys()):
    columns.append(fits.Column(name=key,array=photoz_y[key],format='E'))
hdu = fits.BinTableHDU.from_columns(columns, name='Photo_z')
hdul.append(hdu)
print 'Created table with binned Photo-z'
sys.stdout.flush()

#Write the output file
try:
    os.remove(paths['output_file'])
except:
    pass
hdul.writeto(paths['output_file'])
print hdul.info()
sys.stdout.flush()
print 'Written output file at ' + os.path.relpath(paths['output_file'])
sys.stdout.flush()


#Generate plot
for key in sorted(vrs.Z_BINS.keys()):
    plt.plot(photoz_x, photoz_y[key], label = key)
plt.legend(loc="best", frameon = False, fontsize=10)
if not args.skip_info:
    plt.title(str(n_gal['TOT']) + ' galaxies (' + '{:2.2f}'.format(n_eff['TOT']) + ' per sq arcmin)')
plt.xlabel(r'$z$')
plt.ylabel(r'$P$')
plt.xlim(0.,2.)
plt.savefig(paths['output_plot'])
plt.close()
print 'Saved plot at ' + os.path.relpath(paths['output_plot'])
sys.stdout.flush()

print 'Success!!'
