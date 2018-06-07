import os
import sys
import numpy as np

from astropy.io import fits

import settings


def file_exist_or_error(fname):
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

def read_xipm(fname):
    xi_file = np.loadtxt(fname, dtype='float64', unpack=True)[1]
xi_obs = np.zeros((2, n_theta, n_bins, n_bins))
for count in range(len(xi_file)):
    p_pm, p_theta, p_bin_1, p_bin_2 = find_pos(count, n_theta, n_bins)
    xi_obs[p_pm,p_theta,p_bin_1,p_bin_2] = xi_file[count]
    xi_obs[p_pm,p_theta,p_bin_2,p_bin_1] = xi_file[count]
