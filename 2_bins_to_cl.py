import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import pyccl as ccl
import variables as vrs



# Parse the given arguments
parser = argparse.ArgumentParser("Generate Cl's from binned Photo-z")
parser.add_argument("input_file", type=str, help="Input FITS file")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)


#Open files and store datas
with fits.open(paths['input']) as fn:
    z = fn['Photoz_z'].data
    pz = fn['Photoz_p'].data
#Check that the input quantities have the right dimensions
if (pz.shape[1] != len(z)):
    raise IOError('Input arrays have disagreeing shapes. z.shape = ' + str(z.shape) +
        '. pz.shape = ' + str(pz.shape) + '.')


#Calculate cosmology
cosmo = ccl.Cosmology(Omega_c=vrs.Omega_c, Omega_b=vrs.Omega_b, h=vrs.h, sigma8=vrs.sigma8, n_s=vrs.n_s)#A_s=vrs.A_s
print 'Calculated cosmology'
sys.stdout.flush()

#Calculate tracers
lens = np.array([ccl.ClTracerLensing(cosmo, False, z=z.astype(np.float64), n=pz[x].astype(np.float64)) for x in range(len(pz))])
print 'Calculated tracers'
sys.stdout.flush()

#Calculate Cl's
angular_cl = np.zeros((len(pz), len(pz), vrs.L_MAX+1))
for count1 in range(len(pz)):
    for count2 in range(len(pz)):
        lens1 = lens[count1]
        lens2 = lens[count2]
        angular_cl[count1][count2] = ccl.angular_cl(cosmo, lens1, lens2, range(vrs.L_MAX+1))
angular_cl = np.swapaxes(angular_cl,2,0)
print 'Calculated Cl\'s'
sys.stdout.flush()


#Update input file with an image containing Cl's
#Create image
image_name = 'Cls'
hdu = fits.ImageHDU(angular_cl, name=image_name)
#Update existing file
with fits.open(paths['input'], mode='update') as hdul:
    try:
        hdul.__delitem__(image_name)
    except:
        pass
    hdul.append(hdu)
    hdul.flush()
    print hdul.info()
    sys.stdout.flush()
print 'Updated input file at ' + os.path.relpath(paths['input'])
sys.stdout.flush()


print 'Success!!'
