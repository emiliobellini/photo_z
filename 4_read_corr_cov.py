import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import variables as vrs



# Parse the given arguments
parser = argparse.ArgumentParser("Read the correlation function and the covariance matrix, then apply the KL transform and store everything")
parser.add_argument("input_file", type=str, help="Input FITS file")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)


#Open files and store datas
with fits.open(paths['input']) as fn:
    kl_t = fn['kl_t_avg'].data
n_bins = len(kl_t)


#Angular separation
theta_deg = np.array(vrs.THETA_ARCMIN)/60. # Theta is in degrees
n_theta = len(theta_deg)


#Given the position in the file, find the corresponding position in the array
def find_pos(n, n_theta, n_bins):
    
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


#Function to check symmetries
def check_symm(arr, perm):
    if len(arr.shape) != len(perm):
        raise ValueError('The number of indices of the array is different from the permutation provided')
    test_n = abs(np.transpose(arr, perm) - arr)
    test_d = abs(np.transpose(arr, perm) + arr)/2.
    test_d[test_d < vrs.MP_THRESHOLD] = vrs.MP_THRESHOLD
    if (test_n/test_d).max() > vrs.MP_THRESHOLD:
        print 'WARNING: the array is not symmetric w.r.t. the permutation ' + str(perm) + ' !!!'

#Read observed correlation function
xi_file = np.loadtxt(os.path.abspath('./xipm.dat'), dtype='float64', unpack=True)[1]
xi_obs = np.zeros((2, n_theta, n_bins, n_bins))
for count in range(len(xi_file)):
    p_pm, p_theta, p_bin_1, p_bin_2 = find_pos(count, n_theta, n_bins)
    xi_obs[p_pm,p_theta,p_bin_1,p_bin_2] = xi_file[count]
    xi_obs[p_pm,p_theta,p_bin_2,p_bin_1] = xi_file[count]
#Check that the correlation function is symmetric
check_symm(xi_obs, (0, 1, 3, 2))
print 'Read correlation function'
sys.stdout.flush()


#Read covariance matrix
cov_file = np.loadtxt(os.path.abspath('./xipm_cov.dat'), dtype='float64', unpack=True)[2]
cov_mat = np.zeros((2, 2, n_theta, n_theta, n_bins, n_bins, n_bins, n_bins))
n_max = (n_theta-len(vrs.NEGLECT_THETA_PLUS))*n_bins*(n_bins+1)/2 + (n_theta-len(vrs.NEGLECT_THETA_MINUS))*n_bins*(n_bins+1)/2
for count in range(len(cov_file)):
    div, mod = np.divmod(count, n_max)
    p_pm_1, p_theta_1, p_i_1, p_j_1 = find_pos(div, n_theta, n_bins)
    p_pm_2, p_theta_2, p_i_2, p_j_2 = find_pos(mod, n_theta, n_bins)
    cov_mat[p_pm_1,p_pm_2,p_theta_1,p_theta_2,p_i_1,p_i_2,p_j_1,p_j_2] = cov_file[count]
    cov_mat[p_pm_1,p_pm_2,p_theta_1,p_theta_2,p_i_1,p_j_2,p_j_1,p_i_2] = cov_file[count]
    cov_mat[p_pm_1,p_pm_2,p_theta_1,p_theta_2,p_j_1,p_i_2,p_i_1,p_j_2] = cov_file[count]
    cov_mat[p_pm_1,p_pm_2,p_theta_1,p_theta_2,p_j_1,p_j_2,p_i_1,p_i_2] = cov_file[count]
#Check that the covariance matrix is symmetric
check_symm(cov_mat, (0, 1, 2, 3, 4, 7, 6, 5))
check_symm(cov_mat, (0, 1, 2, 3, 6, 5, 4, 7))
check_symm(cov_mat, (0, 1, 2, 3, 6, 7, 4, 5))
check_symm(cov_mat, (1, 0, 3, 2, 5, 4, 7, 6))
print 'Read covariance matrix'
sys.stdout.flush()


#KL transform (correlation function)
xi_obs_kl = kl_t.dot(xi_obs).dot(kl_t.T)
xi_obs_kl = np.transpose(xi_obs_kl, axes=[1, 2, 0, 3])
xi_obs_kl = np.diagonal(xi_obs_kl,  axis1=2, axis2=3)
#KL transform (covariance matrix)
cov_mat_kl = kl_t.dot(cov_mat).dot(kl_t.T)
cov_mat_kl = np.transpose(cov_mat_kl, axes=[1, 2, 3, 4, 0, 7, 5, 6])
cov_mat_kl = kl_t.dot(cov_mat_kl).dot(kl_t.T)
cov_mat_kl = np.transpose(cov_mat_kl, axes=[1, 2, 3, 4, 5, 6, 0, 7])
cov_mat_kl = np.diagonal(cov_mat_kl,  axis1=4, axis2=6)
cov_mat_kl = np.diagonal(cov_mat_kl,  axis1=4, axis2=5)


#Update input file
hdu = {}
hdu['theta'] = fits.ImageHDU(theta_deg, name='theta')
hdu['xi_obs'] = fits.ImageHDU(xi_obs, name='xi_obs')
hdu['cov_mat'] = fits.ImageHDU(cov_mat, name='cov_mat')
hdu['xi_obs_kl'] = fits.ImageHDU(xi_obs_kl, name='xi_obs_kl')
hdu['cov_mat_kl'] = fits.ImageHDU(cov_mat_kl, name='cov_mat_kl')
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
