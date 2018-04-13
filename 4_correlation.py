import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import pyccl as ccl
import variables as vrs



# Parse the given arguments
parser = argparse.ArgumentParser("Calculate the correlation function, apply the KL transform and calculate the chi^2")
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
n_bins = len(kl_t)
n_ells = vrs.L_MAX+1


#Angular separation
theta = np.array(vrs.THETA_ARCMIN)/60. # Theta is in degrees
n_theta = len(theta)


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
    if (test_n/test_d).max() > 1.e2*vrs.MP_THRESHOLD:
        print 'WARNING: the array is not symmetric w.r.t. the permutation ' + str(perm) + ' !!!'


#Read observed correlation function
xi_file = np.loadtxt(os.path.abspath(vrs.CORR_FILE), dtype='float64', unpack=True)[1]
xi_obs = np.zeros((2, n_theta, n_bins, n_bins))
for count in range(len(xi_file)):
    p_pm, p_theta, p_bin_1, p_bin_2 = find_pos(count, n_theta, n_bins)
    xi_obs[p_pm,p_theta,p_bin_1,p_bin_2] = xi_file[count]
    xi_obs[p_pm,p_theta,p_bin_2,p_bin_1] = xi_file[count]
#Check that the correlation function is symmetric
check_symm(xi_obs, (0, 1, 3, 2))
print 'Read observed correlation function'
sys.stdout.flush()


#Read covariance matrix
cov_file = np.loadtxt(os.path.abspath(vrs.COVMAT_FILE), dtype='float64', unpack=True)[2]
cov_mat = np.zeros((2, n_theta, n_bins, n_bins, 2, n_theta, n_bins, n_bins))
n_max = (n_theta-len(vrs.NEGLECT_THETA_PLUS))*n_bins*(n_bins+1)/2 + (n_theta-len(vrs.NEGLECT_THETA_MINUS))*n_bins*(n_bins+1)/2
for count in range(len(cov_file)):
    div, mod = np.divmod(count, n_max)
    p_pm_1, p_theta_1, p_i_1, p_j_1 = find_pos(div, n_theta, n_bins)
    p_pm_2, p_theta_2, p_i_2, p_j_2 = find_pos(mod, n_theta, n_bins)
    cov_mat[p_pm_1,p_theta_1,p_i_1,p_j_1,p_pm_2,p_theta_2,p_i_2,p_j_2] = cov_file[count]
    cov_mat[p_pm_1,p_theta_1,p_i_1,p_j_1,p_pm_2,p_theta_2,p_j_2,p_i_2] = cov_file[count]
    cov_mat[p_pm_1,p_theta_1,p_j_1,p_i_1,p_pm_2,p_theta_2,p_i_2,p_j_2] = cov_file[count]
    cov_mat[p_pm_1,p_theta_1,p_j_1,p_i_1,p_pm_2,p_theta_2,p_j_2,p_i_2] = cov_file[count]
#Check that the covariance matrix is symmetric
check_symm(cov_mat, (0, 1, 2, 3, 4, 5, 7, 6))
check_symm(cov_mat, (0, 1, 3, 2, 4, 5, 6, 7))
check_symm(cov_mat, (0, 1, 3, 2, 4, 5, 7, 6))
check_symm(cov_mat, (4, 5, 6, 7, 0, 1, 2, 3))
print 'Read covariance matrix'
sys.stdout.flush()


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
check_symm(xi_th, (0, 1, 3, 2))
print 'Calculated theory correlation function'
sys.stdout.flush()


#Mask
mask_p = np.array([x not in vrs.NEGLECT_THETA_PLUS for x in range(n_theta)])
mask_m = np.array([x not in vrs.NEGLECT_THETA_MINUS for x in range(n_theta)])
mask = np.hstack((mask_p,mask_m))
n_data = len(mask[mask])*n_bins*(n_bins+1)/2

#Compute the chi^2
#Reshape and Flatten
xi_th_f = xi_th.reshape((2*n_theta,n_bins,n_bins))
xi_th_f = np.triu(xi_th_f[mask]).flatten()
xi_th_f = xi_th_f[xi_th_f != 0]
xi_obs_f = xi_obs.reshape((2*n_theta,n_bins,n_bins))
xi_obs_f = np.triu(xi_obs_f[mask]).flatten()
xi_obs_f = xi_obs_f[xi_obs_f != 0]
cov_mat_f = cov_mat.reshape((2*n_theta,n_bins,n_bins,2*n_theta,n_bins,n_bins))
cov_mat_f = np.triu(cov_mat_f[:,:,:,mask])
cov_mat_f = np.transpose(cov_mat_f,axes=[3,4,5,0,1,2])
cov_mat_f = np.triu(cov_mat_f[:,:,:,mask]).flatten()
cov_mat_f = cov_mat_f[cov_mat_f != 0]
cov_mat_f = cov_mat_f.reshape((n_data,n_data))
#Inverse cov_mat
inv_cov_mat_f = (vrs.N_SIMULATIONS-n_data-2.)/(vrs.N_SIMULATIONS-1.)*np.linalg.inv(cov_mat_f)
#Chi^2
chi2 = (xi_obs_f-xi_th_f).dot(inv_cov_mat_f).dot(xi_obs_f-xi_th_f)
print 'chi^2 = ' + '{:.2f}'.format(chi2) + ', reduced = ' + '{:.2f}'.format(chi2/(n_data-5.)) + ' (before KL)'
sys.stdout.flush()




#KL transform (correlation function)
xi_th_kl = kl_t.dot(xi_th).dot(kl_t.T)
xi_th_kl = np.transpose(xi_th_kl, axes=[1, 2, 0, 3])
xi_obs_kl = kl_t.dot(xi_obs).dot(kl_t.T)
xi_obs_kl = np.transpose(xi_obs_kl, axes=[1, 2, 0, 3])
#KL transform (covariance matrix)
cov_mat_kl = kl_t.dot(cov_mat).dot(kl_t.T)
cov_mat_kl = np.transpose(cov_mat_kl, axes=[5, 6, 0, 7, 1, 2, 3, 4])
cov_mat_kl = kl_t.dot(cov_mat_kl).dot(kl_t.T)
cov_mat_kl = np.transpose(cov_mat_kl, axes=[1, 2, 3, 4, 5, 6, 0, 7])


#Compute the chi^2
#Reshape and Flatten
xi_obs_f = xi_obs_kl.reshape((2*n_theta,n_bins,n_bins))
xi_obs_f = np.triu(xi_obs_f[mask]).flatten()
xi_obs_f = xi_obs_f[xi_obs_f != 0]
xi_th_f = xi_th_kl.reshape((2*n_theta,n_bins,n_bins))
xi_th_f = np.triu(xi_th_f[mask]).flatten()
xi_th_f = xi_th_f[xi_th_f != 0]
cov_mat_f = cov_mat_kl.reshape((2*n_theta,n_bins,n_bins,2*n_theta,n_bins,n_bins))
cov_mat_f = np.triu(cov_mat_f[:,:,:,mask])
cov_mat_f = np.transpose(cov_mat_f,axes=[3,4,5,0,1,2])
cov_mat_f = np.triu(cov_mat_f[:,:,:,mask]).flatten()
cov_mat_f = cov_mat_f[cov_mat_f != 0]
cov_mat_f = cov_mat_f.reshape((n_data,n_data))
#Inverse cov_mat
inv_cov_mat_f = (vrs.N_SIMULATIONS-n_data-2.)/(vrs.N_SIMULATIONS-1.)*np.linalg.inv(cov_mat_f)
#Chi^2
chi2 = (xi_obs_f-xi_th_f).dot(inv_cov_mat_f).dot(xi_obs_f-xi_th_f)
print 'chi^2 = ' + '{:.2f}'.format(chi2) + ', reduced = ' + '{:.2f}'.format(chi2/(n_data-5.)) + ' (after KL, without diagonal)'
sys.stdout.flush()


#KL transform (correlation function)
xi_th_kl = np.diagonal(xi_th_kl,  axis1=2, axis2=3)
xi_obs_kl = np.diagonal(xi_obs_kl,  axis1=2, axis2=3)
#KL transform (covariance matrix)
cov_mat_kl = np.diagonal(cov_mat_kl,  axis1=6, axis2=7)
cov_mat_kl = np.diagonal(cov_mat_kl,  axis1=2, axis2=3)
cov_mat_kl = np.transpose(cov_mat_kl, axes=[0,1,5,2,3,4])


#Compute the chi^2
#Apply mask
xi_obs_f = xi_obs_kl.reshape((2*n_theta,n_bins))
xi_obs_f = xi_obs_f[mask]
xi_th_f = xi_th_kl.reshape((2*n_theta,n_bins))
xi_th_f = xi_th_f[mask]
cov_mat_f = cov_mat_kl.reshape((2*n_theta,n_bins,2*n_theta,n_bins))
cov_mat_f = cov_mat_f[:,:,mask]
cov_mat_f = cov_mat_f[mask]
for count in range(n_bins):
    n_data = len(mask[mask])*(n_bins-count)
    #Select bins
    xi_obs_f_bin = xi_obs_f[:,:n_bins-count].flatten()
    xi_th_f_bin = xi_th_f[:,:n_bins-count].flatten()
    cov_mat_f_bin = cov_mat_f[:,:n_bins-count,:,:n_bins-count].reshape((n_data,n_data))
    #Inverse cov_mat
    inv_cov_mat_f_bin = (vrs.N_SIMULATIONS-n_data-2.)/(vrs.N_SIMULATIONS-1.)*np.linalg.inv(cov_mat_f_bin)
    #Chi^2
    chi2 = (xi_obs_f_bin-xi_th_f_bin).dot(inv_cov_mat_f_bin).dot(xi_obs_f_bin-xi_th_f_bin)
    print 'chi^2 = ' + '{:.2f}'.format(chi2) + ', reduced = ' + '{:.2f}'.format(chi2/(n_data-5.)) + ' (after KL, diagonal with ' +str(n_bins-count)+ ' bins)'
    sys.stdout.flush()


#Update input file
hdu = {}
hdu['theta'] = fits.ImageHDU(theta, name='theta')
hdu['xi_th'] = fits.ImageHDU(xi_th, name='xi_th')
hdu['xi_obs'] = fits.ImageHDU(xi_obs, name='xi_obs')
hdu['cov_mat'] = fits.ImageHDU(cov_mat, name='cov_mat')
hdu['xi_th_kl'] = fits.ImageHDU(xi_th_kl, name='xi_th_kl')
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
