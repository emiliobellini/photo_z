import os
import sys

from astropy.io import fits


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
def position_xipm_in_array(
    n,
    n_theta,
    n_bins,
    ):

    #Neglect firs/last thetas for pm
    neg_p = vrs.NEGLECT_THETA_PLUS
    neg_m = vrs.NEGLECT_THETA_MINUS

    p_max = (n_theta-len(neg_p))*n_bins*(n_bins+1)/2 + (n_theta-len(neg_m))*n_bins*(n_bins+1)/2
    if n>=p_max:
        raise ValueError("The input number is larger than expected!")

    div, mod = np.divmod(n, 2*n_theta-len(neg_p)-len(neg_m))

    if mod<n_theta-len(neg_p):
        p_pm = 0
        p_theta = np.delete(np.arange(n_theta), neg_p)[mod]
    else:
        p_pm = 1
        p_theta = np.delete(np.arange(n_theta), neg_m)[mod-(n_theta-len(neg_p))]

    intervals = np.flip(np.array([np.arange(x,n_theta+1).sum() for x in np.arange(2,n_theta+2)]),0)
    p_bin_1 = np.where(intervals<=div)[0][-1]
    p_bin_2 = div - intervals[p_bin_1] + p_bin_1

    return p_pm, p_theta, p_bin_1, p_bin_2
