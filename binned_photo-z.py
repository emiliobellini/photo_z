import os, sys, fnmatch
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from astropy.io import fits
from astropy.table import Table, vstack



#Photo-z Bins
z_mins = [0.15,0.29,0.43,0.57,0.70,0.90,1.10,1.30]
z_bins = {}
for n in np.arange(len(z_mins)-1):
    key = str(z_mins[n])+'<z_B<'+str(z_mins[n+1])
    z_bins[key] = [z_mins[n], z_mins[n+1]]


#Criteria used to select the data
def select_data(data, z_min, z_max):

    sel = data['Z_B']>=z_min
    sel = (data['Z_B']<z_max)*sel
    sel = (data['MASK']<=1)*sel
    sel = ([x in good_fit_patterns for x in data['field']])*sel

    return sel


#Good fit patterns
good_fit_patterns = ['W1m0m0', 'W1m0m3', 'W1m0m4', 'W1m0p1', 'W1m0p2', 'W1m0p3', 'W1m1m0',
                     'W1m1m2', 'W1m1m3', 'W1m1m4', 'W1m1p3', 'W1m2m1', 'W1m2m2', 'W1m2m3',
                     'W1m2p1', 'W1m2p2', 'W1m3m0', 'W1m3m2', 'W1m3m4', 'W1m3p1', 'W1m3p3',
                     'W1m4m0', 'W1m4m1', 'W1m4m3', 'W1m4m4', 'W1m4p1', 'W1p1m1', 'W1p1m2',
                     'W1p1m3', 'W1p1m4', 'W1p1p1', 'W1p1p2', 'W1p1p3', 'W1p2m0', 'W1p2m2',
                     'W1p2m3', 'W1p2m4', 'W1p2p1', 'W1p2p2', 'W1p2p3', 'W1p3m1', 'W1p3m2',
                     'W1p3m3', 'W1p3m4', 'W1p3p1', 'W1p3p2', 'W1p3p3', 'W1p4m0', 'W1p4m1',
                     'W1p4m2', 'W1p4m3', 'W1p4m4', 'W1p4p1', 'W1p4p2', 'W1p4p3',

                     'W2m0m0', 'W2m0m1', 'W2m0p1', 'W2m0p2', 'W2m1m0', 'W2m1m1', 'W2m1p1',
                     'W2m1p3', 'W2p1m0', 'W2p1p1', 'W2p1p2', 'W2p2m0', 'W2p2m1', 'W2p2p1',
                     'W2p2p2', 'W2p3m0', 'W2p3m1', 'W2p3p1', 'W2p3p3',

                     'W3m0m1', 'W3m0m2', 'W3m0m3', 'W3m0p2', 'W3m0p3', 'W3m1m0', 'W3m1m2',
                     'W3m1m3', 'W3m1p1', 'W3m1p2', 'W3m1p3', 'W3m2m1', 'W3m2m2', 'W3m2m3',
                     'W3m2p1', 'W3m2p2', 'W3m3m0', 'W3m3m1', 'W3m3m2', 'W3m3m3', 'W3m3p1',
                     'W3m3p2', 'W3p1m0', 'W3p1m1', 'W3p1m2', 'W3p1m3', 'W3p1p2', 'W3p1p3',
                     'W3p2m0', 'W3p2m3', 'W3p2p3', 'W3p3m1', 'W3p3m3', 'W3p3p1', 'W3p3p2',
                     'W3p3p3',

                     'W4m0m2', 'W4m0p1', 'W4m1m0', 'W4m1m1', 'W4m1m2', 'W4m1p1', 'W4m2m0',
                     'W4m2p1', 'W4m2p3', 'W4m3m0', 'W4m3p1', 'W4m3p2', 'W4m3p3', 'W4p1m0',
                     'W4p1m1', 'W4p1m2', 'W4p2m0', 'W4p2m1', 'W4p2m2']


# Parse the given arguments
parser = argparse.ArgumentParser("Generate binned Photo-z")
parser.add_argument("input_file", type=str, help="Input FITS file")
parser.add_argument("--output_file", "-o", type=str, default = None, help="Output file")
parser.add_argument("--output_plot", "-p", type=str, default = None, help="Output plot")
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
sel_bins = select_data(table, z_mins[0], z_mins[-1])
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
for key in z_bins.keys():
    sel_bins[key] = select_data(table, z_bins[key][0], z_bins[key][1])


#Generate binned Photo-z

#Step size
z_step = 0.05
photoz_x = (np.arange(len(image[0]))+1./2.)*z_step

#Normalize photoz
image = (image.T / image.sum(axis=1)).T / z_step
print 'Normalized Photo-z'
sys.stdout.flush()

#Calculate binned Photo-z
photoz_y = {}
for key in z_bins.keys():
    w_tot = table['weight'][sel_bins[key]].sum()
    photoz_y[key] = np.dot(table['weight'][sel_bins[key]].data, image[sel_bins[key]])/w_tot
    print 'Calculated Photo-z for bin ' + key
    sys.stdout.flush()


#Save Photo-z data

#Create the first empty image
hdu = fits.PrimaryHDU()
hdul = fits.HDUList([hdu])

#Create a table with all the galaxies for each bin
for key in sorted(z_bins.keys()):
    columns = []
    hdu = table[sel_bins[key]]
    hdu.meta['EXTNAME'] = key
    hdu = fits.table_to_hdu(hdu)
    hdul.append(hdu)
    print 'Created table for bin ' + key
    sys.stdout.flush()

#Create table with binned Photo-z
columns = []
columns.append(fits.Column(name='z',array=photoz_x,format='E'))
for key in sorted(photoz_y.keys()):
    columns.append(fits.Column(name=key,array=photoz_y[key],format='E'))
hdu = fits.BinTableHDU.from_columns(columns, name='Photo-z')
hdul.append(hdu)
print 'Created table with binned Photo-z'
sys.stdout.flush()

#Write the output file
try:
    os.remove(paths['output_file'])
except:
    pass
hdul.writeto(paths['output_file'])
print 'Written output file at ' + os.path.relpath(paths['output_file'])
sys.stdout.flush()


#Generate plot
for key in sorted(z_bins.keys()):
    plt.plot(photoz_x, photoz_y[key], label = key)
plt.legend(loc="best", frameon = False, fontsize=10)
plt.title(str(n_gal['TOT']) + ' galaxies (' + '{:2.2f}'.format(n_eff['TOT']) + ' per sq arcmin)')
plt.xlabel(r'$z$')
plt.ylabel(r'$P$')
plt.xlim(0.,2.)
plt.savefig(paths['output_plot'])
plt.close()
print 'Saved plot at ' + os.path.relpath(paths['output_plot'])
sys.stdout.flush()

print 'Success!!'
