import os, sys
import argparse
import csv
import numpy as np
from astropy.io import fits


#Fields to be stored in images (all the rest will be stored in a table)
#Chenge this if you want different fields to be stored in images
image_fields = ['PZ_full']


# Parse the given arguments
parser = argparse.ArgumentParser("Read a TSV file and convert it to FITS format")
parser.add_argument("input_file", type=str, help="Input TSV file")
parser.add_argument("--output_file", "-o", type=str, default = None, help="Output FITS file")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)
if args.output_file:
    paths['output'] = os.path.abspath(args.output_file)
else:
    paths['output'] = os.path.splitext(paths['input'])[0] + '.fits'


def floatify(x):
    try:
        return float(x)
    except:
        return x

def listify(data):
    try:
        array = data.split(',')
        array[1]
        return [floatify(x) for x in array]
    except:
        return data

#Read the TSV file and store the info in a dictionary
count = 1
data = {}
with open(paths['input']) as f:
    tsv_file = csv.DictReader(f, dialect='excel-tab')
    for key in tsv_file.fieldnames:
        data[key] = []
    for row in tsv_file:
        for key in tsv_file.fieldnames:
            value = row[key]
            value = listify(value)
            value = floatify(value)
            data[key].append(value)
        if np.mod(count,1e6)==0:
            print '----> Read first ' + '%e' % count +  ' rows'
            sys.stdout.flush()
        count = count + 1
print 'Finished reading data!'
sys.stdout.flush()


#Create the first empty image and store it as the first element of hdul
hdu = fits.PrimaryHDU()
hdul = fits.HDUList([hdu])

#Create the table and append it to hdul
table_fields = [x for x in data.keys() if x not in image_fields]
columns = []
for key in table_fields:
    if type(data[key][0]) is str:
        l = str(np.max([len(x) for x in data[key]]))
        columns.append(fits.Column(name=key,array=data[key],format=l+'A'))
    if type(data[key][0]) is list:
        l = str(np.max([len(x) for x in data[key]]))
        if type(data[key][0][0]) is float:
            columns.append(fits.Column(name=key,array=data[key],format=l+'E'))
    if type(data[key][0]) is float:
        columns.append(fits.Column(name=key,array=data[key],format='E'))
hdu = fits.BinTableHDU.from_columns(columns, name='Table')
hdul.append(hdu)
print 'Created table!'
sys.stdout.flush()

#Create the images and append them to hdul
for key in image_fields:
    hdu = fits.ImageHDU(data[key], name=key)
    hdul.append(hdu)
    print 'Created image for ' + key + '!'
    sys.stdout.flush()


#Write the output file
try:
    os.remove(paths['output'])
except:
    pass
hdul.writeto(paths['output'])
print 'Done! Written output file at ' + os.path.relpath(paths['output'])
