import os, fnmatch
import argparse
import numpy as np
from astropy.io import fits


#Photo-z Bins
z_mins = [0.15,0.29,0.43,0.57,0.70,0.90,1.10,1.30]
z_bins = [[z_mins[n], z_mins[n+1]] for n in np.arange(len(z_mins)-1)]


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
parser = argparse.ArgumentParser("Generate plots of binned Photo-z")
parser.add_argument("input_path", type=str, help="Input FITS path")
parser.add_argument("--output_file", "-o", type=str, default = None, help="Output file")
parser.add_argument("--output_plot", "-p", type=str, default = None, help="Output plot")
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
    paths['output_file'] = os.path.abspath(args.output_file)
else:
    paths['output_file'] = os.path.abspath(args.input_path)
    paths['output_file'] = os.path.splitext(paths['output_file'])[0]
    paths['output_file'] = paths['output_file'].rstrip('_')
    paths['output_file'] = paths['output_file'] + '_bins.fits'
if args.output_plot:
    paths['output_plot'] = os.path.abspath(args.output_plot)
else:
    paths['output_plot'] = os.path.splitext(paths['output_file'])[0]
    paths['output_plot'] = paths['output_plot'] + '.pdf'


#Open files and store datas
table = []
image = []
for fn in paths['input']:
    table.append(fits.open(fn)[1].data)
    image.append(fits.open(fn)[2].data)


#Criteria used to select the data
def select_data_pos(data, z_min, z_max):
    pos = np.argwhere(data['Z_B']>=z_min).flatten()
    pos = np.intersect1d(pos, np.argwhere(data['Z_B']<z_max).flatten())
    pos = np.intersect1d(pos, np.argwhere(data['MASK']<=1).flatten())
#    pos = np.intersect1d(pos, np.argwhere(data['field'] in good_fit_patterns).flatten())

    return np.array(pos)


#Get an array with positions of data for each bin
pos_bins = {}
for [z_min, z_max] in z_bins:
    key = str(z_min)+'<z_B<'+str(z_max)
    pos_bins[key] = []
    for count in range(len(table)):
        pos_bins[key].append(select_data_pos(table[count], z_min, z_max))


#Generate binned Photo-z
z_step = 0.05
photoz_x = (np.arange(len(image[0][0]))+1./2.)*z_step
photoz_y = {}
for key in pos_bins.keys():
    for count in range(len(image)):
        print pos_bins[key][count]