import os, sys
import argparse
import csv
import numpy as np
from astropy.io import fits

#Redshift bins width used to divide the files
Z_WIDTH = 0.2
#Threshold redshift after which all data are stored in the same file
Z_TH = 1.3
#Minimum and Maximum redshift of the data
Z_MIN = 0.0
Z_MAX = 3.5
#Fields to be stored in images (all the rest will be stored in a table)
#Chenge this if you want different fields to be stored in images
image_fields = ['PZ_full']


# Parse the given arguments
parser = argparse.ArgumentParser("Read a TSV file and convert it to FITS format")
parser.add_argument("input_file", type=str, help="Input TSV file")
parser.add_argument("--output_file", "-o", type=str, default = None, help="Output FITS file")
args = parser.parse_args()


#Calculate the redshift bin for each file
#z_mins = [Z_MIN+Z_WIDTH*n for n in np.arange(np.trunc((Z_TH-Z_MIN)/Z_WIDTH))]
num_bins = int((Z_TH-Z_MIN)/Z_WIDTH) + 2
z_mins = [round(Z_MIN+Z_WIDTH*n,4) for n in np.arange(num_bins)]
z_mins.append(round(Z_MAX,4))
z_bins = [[z_mins[n], z_mins[n+1]] for n in np.arange(num_bins)]


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)
paths['output'] = []
if args.output_file:
    tmp_output = os.path.splitext(os.path.abspath(args.output_file))
else:
    tmp_output = (os.path.splitext(paths['input'])[0], '.fits')
for count in np.arange(num_bins):
    fp = '__' + str(z_bins[count][0]) + 'z' + str(z_bins[count][1])
    paths['output'].append(tmp_output[0] + fp + tmp_output[1])


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


#Read the TSV file, store the info in a dictionary and write the output files
#The loop is to divide the data in redshift bins and save them in different
#files. This way we save memory.
for num_bin in np.arange(num_bins):
    bin_label = str(z_bins[num_bin][0]) + '<=z<' + str(z_bins[num_bin][1])
    #Read the TSV file and store the info in a dictionary
    count = 1
    data = {}
    with open(paths['input']) as f:
        tsv_file = csv.DictReader(f, dialect='excel-tab')
        for key in tsv_file.fieldnames:
            data[key] = []
        for row in tsv_file:
            if z_bins[num_bin][0]<=floatify(row['Z_B'])<z_bins[num_bin][1]:
                for key in tsv_file.fieldnames:
                    value = row[key]
                    value = listify(value)
                    value = floatify(value)
                    data[key].append(value)
                if np.mod(count,1e5)==0:
                    print('----> Read first ' + '%e' % count +  ' rows of bin '
                        + bin_label)
                    sys.stdout.flush()
                count = count + 1
    print('Read bin ' + bin_label)
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
    print('Created table for bin' + bin_label)
    sys.stdout.flush()

    #Create the images and append them to hdul
    for key in image_fields:
        hdu = fits.ImageHDU(data[key], name=key)
        hdul.append(hdu)
        print 'Created image for ' + key + '!'
        sys.stdout.flush()


    #Write the output file
    try:
        os.remove(paths['output'][num_bin])
    except:
        pass
    hdul.writeto(paths['output'][num_bin])
    print('Written output file for bin' + bin_label + 'at '
        + os.path.relpath(paths['output']))

print('Success!!')
