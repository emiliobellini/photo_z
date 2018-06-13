import os
import sys
import re
import tarfile
import numpy as np

from astropy.io import fits

import settings


def file_exists_or_error(fname):
    abspath = os.path.abspath(fname)
    if os.path.exists(abspath):
        return abspath
    raise IOError('File '+abspath+' not found!')

def write_to_fits(fname, array, name):
    if not os.path.exists(fname):
        hdul = fits.HDUList([fits.PrimaryHDU()])
        hdul.writeto(fname)
    with fits.open(fname, mode='update') as hdul:
        try:
            hdul.__delitem__(name)
        except:
            pass
        hdul.append(fits.ImageHDU(array, name=name))
    print('Appended ' + name.upper() + ' to ' + os.path.relpath(fname))
    sys.stdout.flush()
    return

def print_info_fits(fname):
    with fits.open(fname) as hdul:
        print(hdul.info())
        sys.stdout.flush()
    return

#Given the position in the file, find the corresponding position in the array
def position_xipm(n):
    n_bins = len(settings.Z_BINS)-1
    n_theta = len(settings.THETA_ARCMIN)
    n_theta_xip = np.array(settings.MASK_THETA[0]).astype(int).sum()
    n_theta_xim = np.array(settings.MASK_THETA[1]).astype(int).sum()
    p_max = (n_theta_xip+n_theta_xim)*n_bins*(n_bins+1)/2
    if n>=p_max:
        raise ValueError("The input number is larger than expected!")
    div, mod = np.divmod(n, n_theta_xip+n_theta_xim)
    if mod<n_theta_xip:
        p_pm = 0
        p_theta = mod
    else:
        p_pm = 1
        p_theta = 3+mod-n_theta_xip
    intervals = np.flip(np.array([np.arange(x,n_bins+1).sum() for x in np.arange(2,n_bins+2)]),0)
    p_bin_1 = np.where(intervals<=div)[0][-1]
    p_bin_2 = div - intervals[p_bin_1] + p_bin_1
    return p_pm, p_theta, p_bin_1, p_bin_2

def unflatten_xipm(array):
    n_bins = len(settings.Z_BINS)-1
    n_theta = len(settings.THETA_ARCMIN)
    xipm = np.zeros((2, n_theta, n_bins, n_bins))
    for count in range(len(array)):
        p_pm, p_theta, p_bin_1, p_bin_2 = position_xipm(count)
        xipm[p_pm,p_theta,p_bin_1,p_bin_2] = array[count]
        xipm[p_pm,p_theta,p_bin_2,p_bin_1] = array[count]
    return xipm

def unpack_and_stack(fname):
    n_bins = len(settings.Z_BINS)-1
    mask_theta = np.array(settings.MASK_THETA)
    n_theta_masked = sum(1 for x in mask_theta.flatten() if x)
    base_name = 'mockxipm/xipm_cfhtlens_sub2real0001_maskCLW1_blind1_z1_z1_athena.dat'
    tar = tarfile.open(fname, 'r')
    n_sims, mod = np.divmod(sum(1 for x in tar.getmembers() if x.isreg()), n_bins*(n_bins+1)/2)
    if mod != 0:
        raise IOError('The number of files in ' + fname + ' is not correct!')
    xipm_sims = np.zeros((n_sims, n_theta_masked*n_bins*(n_bins+1)/2))
    for n_sim in range(n_sims):
        for n_bin1 in range(n_bins):
            for n_bin2 in range(n_bin1, n_bins):
                pos = np.flip(np.arange(n_bins+1),0)[:n_bin1].sum()
                pos = (pos + n_bin2 - n_bin1)*n_theta_masked
                new_name = base_name.replace('real0001', 'real{0:04d}'.format(n_sim+1))
                new_name = new_name.replace('z1_athena', 'z{0:01d}_athena'.format(n_bin1+1))
                new_name = new_name.replace('blind1_z1', 'blind1_z{0:01d}'.format(n_bin2+1))
                f = tar.extractfile(new_name)
                if f:
                    xi = np.loadtxt(f)
                    xi = np.hstack((xi[:,1][mask_theta[0]], xi[:,2][mask_theta[1]]))
                    for i, xi_val in enumerate(xi):
                        xipm_sims[n_sim][pos+i] = xi_val
        if (n_sim+1)%100==0 or n_sim+1==n_sims:
            print('----> Unpacked {}/{} correlation functions'.format(n_sim+1, n_sims))
            sys.stdout.flush()
    return xipm_sims

def read_fits_data(fname):
    hdul = fits.open(fname, memmap=True)
    table = hdul['data'].data
    image = hdul['PZ_full'].data
    z_bins = np.array([[settings.Z_BINS[n], settings.Z_BINS[n+1]] for n in np.arange(len(settings.Z_BINS)-1)])
    sel_bins = np.array([settings.get_mask(table, z_bins[n][0], z_bins[n][1]) for n in range(len(z_bins))])
    photo_z = np.zeros((len(z_bins)+1,len(image[0])))
    n_eff = np.zeros(len(z_bins))
    sigma_g = np.zeros(len(z_bins))
    photo_z[0] = (np.arange(len(image[0]))+1./2.)*settings.CFHTlens_dZ
    for n in range(len(z_bins)):
        w_sum = table['weight'][sel_bins[n]].sum()
        w2_sum = (table['weight'][sel_bins[n]]**2.).sum()
        #TODO: Correct ellipticities
        m = np.average(table['e1'][sel_bins[n]])
        e1 = table['e1'][sel_bins[n]]/(1+m)
        e2 = table['e2'][sel_bins[n]]-table['c2'][sel_bins[n]]/(1+m)
        photo_z[n+1] = np.dot(table['weight'][sel_bins[n]], image[sel_bins[n]])/w_sum
        n_eff[n] = w_sum**2/w2_sum/settings.CFHTlens_A_eff
        sigma_g[n] = np.dot(table['weight'][sel_bins[n]]**2., (e1**2. + e2**2.)/2.)/w2_sum
        sigma_g[n] = sigma_g[n]**0.5
        print('----> Completed bin {}/{}'.format(n+1, len(z_bins)))
        sys.stdout.flush()
    return photo_z, n_eff, sigma_g
