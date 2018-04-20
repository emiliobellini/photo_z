import os, sys, re
import argparse
from astropy.io import fits
import numpy as np
import pyccl as ccl
import emcee



# Parse the given arguments
parser = argparse.ArgumentParser("Calculate the correlation function, apply the KL transform and calculate the chi^2")
parser.add_argument("input_file", type=str, help="Input parameters file")
parser.add_argument("data_file", type=str, help="Input data file")
parser.add_argument("--output_dir", "-o", type=str, default = None, help="Output dir")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)
paths['data'] = os.path.abspath(args.data_file)
if args.output_dir:
    paths['output_dir'] = os.path.abspath(args.output_dir)
else:
    paths['output_dir'] = os.path.splitext(paths['input'])[0] + '.dat'


#Default parameters [h, Omega_c, Omega_b, sigma8, n_s]
cosmo_pars = np.array([
    [None, 0.61197750, None],
    [None, 0.31111823, None],
    [None, 0.08743234, None],
    [None, 0.60773680, None],
    [None, 1.25771300, None]
    ])
n_kl = 3
n_ells = 2001
n_walkers = 50
n_steps = 1000
n_threads = 1
n_sim = 1988


#Read data file
with fits.open(paths['data']) as fn:
    z = fn['photoz_z'].data
    pz = fn['photoz_p'].data
    kl_t = fn['kl_t_avg'].data
    theta = fn['theta'].data
    xi_obs = fn['xi_obs_kl'].data
    cov_mat = fn['cov_mat_kl'].data
    mask_theta = fn['mask_theta'].data.astype(bool)
n_bins = len(pz)
n_theta = len(theta)


#Reshape correlation function and covariance matrix
xi_obs = xi_obs.reshape((2*n_theta,n_bins))
xi_obs = xi_obs[mask_theta]
xi_obs = xi_obs[:,:n_kl].flatten()
cov_mat = cov_mat.reshape((2*n_theta,n_bins,2*n_theta,n_bins))
cov_mat = cov_mat[:,:,mask_theta]
cov_mat = cov_mat[mask_theta]
cov_mat = cov_mat[:,:n_kl,:,:n_kl]
cov_mat = cov_mat.reshape((len(mask_theta[mask_theta])*n_kl,len(mask_theta[mask_theta])*n_kl))
inv_cov_mat = (n_sim-len(mask_theta[mask_theta])*n_kl-2.)/(n_sim-1.)*np.linalg.inv(cov_mat)


#Read input file
def floatify(x):
    try:
        return float(x)
    except:
        return x

def listify(data):
    try:
        array = data.split(',')
        array[1]
        return np.array([floatify(x) for x in array])
    except:
        return data

def read_line(file_path, par):
    with open(file_path) as fn:
        for line in fn:
            if "=" in line and line[0] != '#':
                line = re.sub('#.+', '', line)
                name , value = line.split('=')
                name = name.strip()
                if name == par:
                    value = value.strip()
                    value = listify(value)
                    if type(value) != str:
                        if value[0]>=value[2]:
                            raise IOError('The left bound is larger than the right bound for parameter ' + par)
                        if (value[0]>=value[1]) or (value[1]>=value[2]):
                            raise IOError('The central value is out of the bounds for parameter ' + par)
                    if type(value) is str:
                        value = np.array([None, floatify(value), None])
                    if len(value) != 3:
                        raise IOError('Incorrect length of input for parameter ' + par)
                    return value
        raise ValueError()

try:
    cosmo_pars[0] = read_line(paths['input'], 'h')
except ValueError:
    pass
try:
    cosmo_pars[1] = read_line(paths['input'], 'Omega_c')
except ValueError:
    pass
try:
    cosmo_pars[2] = read_line(paths['input'], 'Omega_b')
except ValueError:
    pass
try:
    cosmo_pars[3] = read_line(paths['input'], 'sigma8')
except ValueError:
    pass
try:
    cosmo_pars[4] = read_line(paths['input'], 'n_s')
except ValueError:
    pass
try:
    n_kl = int(read_line(paths['input'], 'n_kl')[1])
except ValueError:
    pass
try:
    n_ells = int(read_line(paths['input'], 'ell_max')[1]+1)
except ValueError:
    pass
try:
    n_walkers = int(read_line(paths['input'], 'n_walkers')[1])
except ValueError:
    pass
try:
    n_steps = int(read_line(paths['input'], 'n_steps')[1])
except ValueError:
    pass
try:
    n_sim = int(read_line(paths['input'], 'n_sim')[1])
except ValueError:
    pass
try:
    n_threads = int(read_line(paths['input'], 'n_threads')[1])
except ValueError:
    pass


#Define mask for variables
mask_vars = np.array([type(x)==float for x in cosmo_pars[:,0]])
n_dim = len(mask_vars[mask_vars])


#Get random initial points
def get_random(pars, squeeze):
    rnd_pars = np.array([])
    for count in range(len(pars)):
        rnd = pars[count][1] + 2.*(np.random.rand()-.5)*min(pars[count][2]-pars[count][1], pars[count][1]-pars[count][0])/squeeze
        rnd_pars = np.append(rnd_pars, rnd)
    return rnd_pars


#Initial point
vars_0 = np.array([get_random(cosmo_pars[mask_vars], 1.e1) for x in range(n_walkers)])


#Construct full array of cosmological parameters
def get_cosmo(var, full=cosmo_pars, mask=mask_vars):
    pars = np.zeros(len(mask))
    count1 = 0
    for count2 in range(len(pars)):
        if not mask[count2]:
            pars[count2] = full[count2][1]
        else:
            pars[count2] = var[count1]
            count1 = count1+1
    return pars


#Get theory correlation function
def get_theory(var):
    #Get cosmological parameters
    var_tot = get_cosmo(var)
    #Cosmology
    cosmo = ccl.Cosmology(h=var_tot[0], Omega_c=var_tot[1], Omega_b=var_tot[2], sigma8=var_tot[3], n_s=var_tot[4])
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
    #KL transform
    xi_th = kl_t.dot(xi_th).dot(kl_t.T)
    xi_th = np.transpose(xi_th, axes=[1, 2, 0, 3])
    xi_th = np.diagonal(xi_th,  axis1=2, axis2=3)
    xi_th = xi_th.reshape((2*n_theta,n_bins))
    xi_th = xi_th[mask_theta,:n_kl].flatten()

    return xi_th


#Define priors
def lnprior(var, full=cosmo_pars, mask=mask_vars):
    is_in = (full[mask][:,0] <= var).all()
    is_in = is_in*(var <= full[mask][:,2]).all()
    if is_in:
        return 0.0
    return -np.inf


#Define likelihood
def lnlike(var, obs=xi_obs, icov=inv_cov_mat):
    #Get theory
    xi_th = get_theory(var)
    #Get chi2
    chi2 = (obs-xi_th).dot(icov).dot(obs-xi_th)
    return -chi2


#Define posterior
def lnprob(var):
    lp = lnprior(var)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(var)


#Initialize sampler
sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, threads=n_threads)


#Print usefule stuff
print 'Starting the chains!'
print 'Number of threads = ' + str(n_threads)
print 'Number of steps = ' + str(n_steps)
print 'Number of walkers = ' + str(n_walkers)
print 'Number of KL modes = ' + str(n_kl)
print 'Maximum ell = ' + str(n_ells-1)
sys.stdout.flush()

#Creat file
f = open(paths['output_dir'], 'w')
f.close()

for count, result in enumerate(sampler.sample(vars_0, iterations=n_steps, storechain=False)):
    pos = result[0]
    prob = result[1]
    f = open(paths['output_dir'], 'a')
    for k in range(pos.shape[0]):
        f.write(str(prob[k]) + "    " + "    ".join([str(x) for x in pos[k]]) + "\n")#"{0:4d} {1:s}\n".format(k, " ".join(position[k]))
    f.close()
    if (count+1) % 10 == 0:
        print '----> Computed ' + str(count+1) + ' over a total of ' + str(n_steps) + ' steps'
        sys.stdout.flush()


print 'Success!!'
