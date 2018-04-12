import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import pyccl as ccl
import variables as vrs



# Parse the given arguments
parser = argparse.ArgumentParser("Apply the KL transform to the correlation functions and the covariance matrix")
parser.add_argument("input_file", type=str, help="Input FITS file")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)


#Open files and store datas
with fits.open(paths['input']) as fn:
    z = fn['photoz_z'].data
    pz = fn['photoz_p'].data
    kl_t = fn['kl_t_avg'].data
    theta = fn['theta'].data
    xi_obs = fn['xi_obs'].data
    cov_mat = fn['cov_mat'].data
n_bins = len(kl_t)
n_theta = len(theta)
n_ells = vrs.L_MAX+1


#Calculate theory correlation function
#Cosmology
cosmo = ccl.Cosmology(Omega_c=vrs.Omega_c, Omega_b=vrs.Omega_b, h=vrs.h, sigma8=vrs.sigma8, n_s=vrs.n_s)
#Tracers
lens = np.array([ccl.ClTracerLensing(cosmo, False, z=z.astype(np.float64), n=pz[x].astype(np.float64)) for x in range(n_bins)])
#Cl's
ell = np.arange(n_ells)
cls = np.zeros((n_bins, n_bins, n_ells))
for count1 in range(n_bins):
    for count2 in range(n_bins):
        cls[count1,count2] = ccl.angular_cl(cosmo, lens[count1], lens[count2], ell)
cls = np.transpose(cls,axes=[2,0,1])
#Correlation function
xi_th = np.zeros((2, n_bins, n_bins, n_theta))
for count1 in range(n_bins):
    for count2 in range(n_bins):
        for count3 in range(n_theta):
            xi_th[0,count1,count2,count3] = ccl.correlation(cosmo, ell, cls[:,count1,count2], theta[count3], corr_type='L+', method='FFTLog')
            xi_th[1,count1,count2,count3] = ccl.correlation(cosmo, ell, cls[:,count1,count2], theta[count3], corr_type='L-', method='FFTLog')
xi_th = np.transpose(xi_th,axes=[0,3,1,2])
print 'Calculated correlation function'
sys.stdout.flush()


#KL transform (correlation function)
xi_th_kl = kl_t.dot(xi_th).dot(kl_t.T)
xi_th_kl = np.transpose(xi_th_kl, axes=[1, 2, 0, 3])
xi_th_kl = np.diagonal(xi_th_kl,  axis1=2, axis2=3)


#Update input file
hdu = {}
hdu['xi_th'] = fits.ImageHDU(xi_th, name='xi_th')
hdu['xi_th_kl'] = fits.ImageHDU(xi_th_kl, name='xi_th_kl')
with fits.open(paths['input'], mode='update') as hdul:
    for im in hdu.keys():
        try:
            hdul.__delitem__(im)
        except:
            pass
        hdul.append(hdu[im])
    #Update existing file
    hdul.flush()
    print hdul.info()
    sys.stdout.flush()
print 'Updated input file at ' + os.path.relpath(paths['input'])
sys.stdout.flush()


print 'Success!!'
