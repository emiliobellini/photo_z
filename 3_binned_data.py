import os, sys, fnmatch
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
import variables as vrs



# Parse the given arguments
parser = argparse.ArgumentParser("Generate binned Photo-z")
parser.add_argument("input_file", type=str, help="Input FITS file")
parser.add_argument("--output_file", "-o", type=str, default = None, help="Output file")
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


#Open files and store datas
#Read table
table = Table.read(paths['input'])
#Read image
image = fits.open(paths['input'])['PZ_full'].data


#Create redshift bins
z_bins = np.array([[vrs.Z_BINS[n], vrs.Z_BINS[n+1]] for n in np.arange(len(vrs.Z_BINS)-1)])


#Print effective densities for each field
if not args.skip_info:
    n_eff = {}
    n_gal = {'TOT' : 0}
    w_tot = 0.
    w2_tot = 0.
    sel_bins = vrs.select_data(table, np.min(z_bins), np.max(z_bins))
    for field in ['W1', 'W2', 'W3', 'W4', 'TOT']:
        if field != 'TOT':
            w = np.array([x['weight'] for x in table[sel_bins] if field in x['field']])
            n_gal[field] = len(w)
            n_eff[field] = (w.sum())**2./(w**2.).sum()/vrs.CFHTlens_A_eff[field]
            n_gal['TOT'] = n_gal['TOT'] + n_gal[field]
            w_tot = w_tot + w.sum()
            w2_tot = w2_tot + (w**2.).sum()
        n_eff['TOT'] = w_tot**2./w2_tot/vrs.CFHTlens_A_eff['TOT']
        print '----> Survey region ' + field + ':'
        sys.stdout.flush()
        print '--------> n_galaxies = ' + '{:2d}'.format(n_gal[field])
        sys.stdout.flush()
        print '--------> n_eff = ' + '{:2.2f}'.format(n_eff[field]) + ' arcmin^-2'
        sys.stdout.flush()
        print '--------> A_eff = ' + '{:2.2f}'.format(vrs.CFHTlens_A_eff[field]/(60.**2)) + ' deg^2'
        sys.stdout.flush()


#Get arrays with selected data for each bin
sel_bins = np.array([vrs.select_data(table, z_bins[n][0], z_bins[n][1]) for n in range(len(z_bins))])


#Generate binned Photo-z

#Step size (provided by CFHTlens to generate the Photo-z arrays)
photoz_x = (np.arange(len(image[0]))+1./2.)*vrs.CFHTlens_dZ

#Normalize photoz
image = (image.T / image.sum(axis=1)).T/vrs.CFHTlens_dZ
print 'Normalized Photo-z'
sys.stdout.flush()

#Calculate binned Photo-z, n_eff and sigma_g
photoz_y = np.zeros((len(z_bins),len(photoz_x)))
n_eff = np.zeros(len(z_bins))
sigma_g = np.zeros(len(z_bins))
for n in range(len(z_bins)):
    w_sum = table['weight'][sel_bins[n]].sum()
    w2_sum = (table['weight'][sel_bins[n]]**2.).sum()
    m = np.average(table['e1'][sel_bins[n]])
    e1 = table['e1'][sel_bins[n]]/(1+m)
    e2 = table['e2'][sel_bins[n]]-table['c2'][sel_bins[n]]/(1+m)
    photoz_y[n] = np.dot(table['weight'][sel_bins[n]].data, image[sel_bins[n]])/w_sum
    n_eff[n] = w_sum**2/w2_sum/vrs.CFHTlens_A_eff['TOT']
    sigma_g[n] = np.dot(table['weight'][sel_bins[n]].data**2., e1**2. + e2**2.)/w2_sum
    sigma_g[n] = sigma_g[n]**0.5
    print 'Calculated Photo-z, n_eff and sigma_g for bin ' + str(n+1)
    sys.stdout.flush()
print '----> n_eff for each bin is ' + str(['{:2.2e}'.format(i) for i in n_eff])
sys.stdout.flush()
print '----> sigma_g for each bin is ' + str(['{:0.3f}'.format(i) for i in sigma_g])
sys.stdout.flush()


#Save Photo-z data

#Create the first empty image
hdu = fits.PrimaryHDU()
hdul = fits.HDUList([hdu])

#Create image with redshift bins
hdu = fits.ImageHDU(z_bins, name='Z_bins')
hdul.append(hdu)

#Create images with binned n_eff
hdu = fits.ImageHDU(n_eff, name='n_eff')
hdul.append(hdu)

#Create images with binned n_eff
hdu = fits.ImageHDU(sigma_g, name='sigma_g')
hdul.append(hdu)

#Create images with binned Photo-z
hdu = fits.ImageHDU(photoz_x, name='Photoz_z')
hdul.append(hdu)
hdu = fits.ImageHDU(photoz_y, name='Photoz_p')
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


print 'Success!!'
