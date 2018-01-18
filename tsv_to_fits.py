import os
import argparse
import csv
import numpy as np
from astropy.io import fits
import astropy


#Fields to be stored in images (all the rest will be stored in a table)
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
            print 'Read first ' + '%e' % count +  ' rows'
        count = count + 1


#Initialise the hdu list
#hdu = []

#Create the first empty image, hdu[0]
#hdu.append(fits.PrimaryHDU())

#Create the table, hdu[1]
#table_fields = [x for x in data.keys() if x not in image_fields]
#columns = []
#for key in table_fields:
#    if type(data[key][0]) is str:
#        l = str(np.max([len(x) for x in data[key]]))
#        columns.append(fits.Column(name=key,array=data[key],format=l+'A'))
#    if type(data[key][0]) is list:
#        l = str(np.max([len(x) for x in data[key]]))
#        if type(data[key][0][0]) is float:
#            columns.append(fits.Column(name=key,array=data[key],format=l+'E'))
#    if type(data[key][0]) is float:
#        columns.append(fits.Column(name=key,array=data[key],format='E'))


counts = np.array([312, 334, 308, 317])
names = np.array(['NGC1', 'NGC2', 'NGC3', 'NGC4'])
col1 = fits.Column(name='target', format='10A', array=names)
col2 = fits.Column(name='counts', format='J', unit='DN', array=counts)
col3 = fits.Column(name='notes', format='A10')
col4 = fits.Column(name='spectrum', format='1000E')
col5 = fits.Column(name='flag', format='L', array=[True, False, True, True])
#hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])
#coldefs = fits.ColDefs([col1, col2, col3, col4, col5])
#hdu = fits.BinTableHDU.from_columns(coldefs)
print astropy.version.version

#print columns
#columns = fits.ColDefs(columns)
#print fits.BinTableHDU.from_columns(columns)
#hdu.append(fits.BinTableHDU.from_columns(columns))
#fits.BinTableHDU.from_columns
#    print key
#fits.Column(name='RA',array=ra_arr,format='E')
#
#print data.keys()
