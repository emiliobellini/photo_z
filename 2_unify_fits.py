import os, sys, fnmatch
import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack



# Parse the given arguments
parser = argparse.ArgumentParser("Unify several fits files into a single one")
parser.add_argument("input_path", type=str, help="Input FITS path")
parser.add_argument("--output_file", "-o", type=str, default = None, help="Output file")
parser.add_argument("--delete", "-d", help="Delete input files", action="store_true")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
split_input = os.path.split(os.path.abspath(args.input_path))
paths['input'] = []
for fn in os.listdir(split_input[0]):
    if fnmatch.fnmatch(fn, '*.fits'):
        if fnmatch.fnmatch(fn, split_input[1]+'*'):
            paths['input'].append(split_input[0] + '/' + fn)
if args.output_file:
    paths['output'] = os.path.abspath(args.output_file)
else:
    paths['output'] = os.path.abspath(args.input_path)
    paths['output'] = paths['output'].rstrip('_')
    paths['output'] = paths['output'] + '.fits'


#Open all the files and join them in a single one

#Create the first empty image
hdu = fits.PrimaryHDU()
hdul = fits.HDUList([hdu])

#Create table
del hdu
for fn in paths['input']:
    try:
        hdu = vstack([hdu, Table.read(fn)], join_type='exact')
    except NameError:
        hdu = Table.read(fn)
hdu = fits.table_to_hdu(hdu)
hdul.append(hdu)
print 'Created final table'
sys.stdout.flush()

#Get the list of images stored
with fits.open(paths['input'][0]) as imgs:
    images = [x.name for x in imgs if x.is_image and x.size>0]
#Create images
for key in images:
    del hdu
    for fn in paths['input']:
        try:
            hdu = np.vstack((hdu, fits.open(fn)[key].data))
        except NameError:
            hdu = fits.open(fn)[key].data
    hdu = fits.ImageHDU(hdu, name=key)
    hdul.append(hdu)
    print 'Created final image for ' + key
    sys.stdout.flush()

#Write the output file
try:
    os.remove(paths['output'])
except:
    pass
hdul.writeto(paths['output'])
print hdul.info()
sys.stdout.flush()
print 'Written output file at ' + os.path.relpath(paths['output'])
sys.stdout.flush()

#Remove partial files
if args.delete:
    for fn in paths['input']:
        try:
            os.remove(fn)
        except:
            pass

print 'Success!!'

sys.exit()
