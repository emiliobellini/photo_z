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


# #Print effective densities for each field
# #TODO: insert correct values for A_eff
# if not args.skip_info:
#     sel_bins = vrs.select_data(table, vrs.Z_MIN, vrs.Z_MAX)
#     w, w2 = {'TOT': 0.}, {'TOT': 0.}
#     n_eff = {}
#     n_gal = {'TOT': 0}
#     A_eff = {'W1': 42.9, 'W2': 12.1, 'W3': 26.1, 'W4': 13.3}
#     A_eff['TOT'] = np.array(A_eff.values()).sum()
#     for key in ['W1', 'W2', 'W3', 'W4']:
#         w[key] = np.array([x['weight'] for x in table[sel_bins] if key in x['field']])
#         n_gal[key] = len(w[key])
#         w2[key] = (w[key]**2).sum()
#         w[key] = w[key].sum()
#         n_eff[key] = w[key]**2/w2[key]/(A_eff[key]*60.*60.)
#         n_gal['TOT'] = n_gal['TOT'] + n_gal[key]
#         w['TOT'] = w['TOT'] + w[key]
#         w2['TOT'] = w2['TOT'] + w2[key]
#         print('Selected ' + str(n_gal[key]) + ' galaxies for field ' + key
#             + ' (' + '{:2.4f}'.format(n_eff[key]) + ' per sq arcmin).')
#         sys.stdout.flush()
#     n_eff['TOT'] = w['TOT']**2/w2['TOT']/(A_eff['TOT']*60.*60.)
#     print('Selected in total ' + str(n_gal['TOT']) + ' galaxies ('
#         + '{:2.4f}'.format(n_eff['TOT']) + ' per sq arcmin).')
#     sys.stdout.flush()


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
    ellipticity = table['e1'][sel_bins[n]]**2. + table['e2'][sel_bins[n]]**2.
    photoz_y[n] = np.dot(table['weight'][sel_bins[n]].data, image[sel_bins[n]])/w_sum
    n_eff[n] = w_sum**2/w2_sum/vrs.CFHTlens_A_eff
    sigma_g[n] = np.dot(table['weight'][sel_bins[n]].data**2., ellipticity)/w2_sum
    sigma_g[n] = sigma_g[n]**0.5
    print 'Calculated Photo-z, n_eff and sigma_g for bin ' + str(n+1)
    sys.stdout.flush()
print '---> n_eff for each bin is ' + str(['{:2.2e}'.format(i) for i in n_eff])
print '---> sigma_g for each bin is ' + str(['{:0.3f}'.format(i) for i in sigma_g])


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
