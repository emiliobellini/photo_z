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


#Read data file
with fits.open(paths['data']) as fn:
    z = fn['photoz_z'].data
    pz = fn['photoz_p'].data
    kl_t = fn['kl_t_avg'].data
    theta = fn['theta'].data
    xi_obs = fn['xi_obs_kl'].data
    cov_mat = fn['cov_mat_kl'].data
    mask = fn['mask_theta'].data.astype(bool)
n_bins = len(pz)
n_theta = len(theta)


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


#Default parameters
pars = {
    'h' : 0.61197750,
    'Omega_c' : 0.31111823454300624,
    'Omega_b' : 0.08743233863669807,
    'sigma8' : 0.60773680,
    'n_s' : 1.2577130,
    'ell_max' : 2000,
    'n_kl' : 3,
    'n_sim' : 1988,
    'n_walkers' : 50,
    'n_steps' : 1000
    }
cosmo_pars = ['h', 'Omega_c', 'Omega_b', 'sigma8', 'n_s']


#Read input file
with open(paths['input']) as fn:
    for line in fn:
        if "=" in line and line[0] != '#':
            line = re.sub('#.+', '', line)
            name , value = line.split('=')
            name = name.strip()
            if name in pars.keys():
                value = value.strip()
                value = floatify(listify(value))
                pars[name] = value
n_ells = int(pars['ell_max']+1)
n_walkers = int(pars['n_walkers'])
n_steps = int(pars['n_steps'])

#Build array with constant and variable parameters
cosmo_vars = np.array([])
for count, par in enumerate(cosmo_pars):
    if type(pars[par]) is float:
        cosmo_pars[count] = pars[par]
    else:
        cosmo_vars = np.append(cosmo_vars, pars[par][1])
n_dim = len(cosmo_vars)


#Reshape correlation function and covariance matrix
xi_obs = xi_obs.reshape((2*n_theta,n_bins))
xi_obs = xi_obs[mask]
xi_obs = xi_obs[:,:pars['n_kl']].flatten()
cov_mat = cov_mat.reshape((2*n_theta,n_bins,2*n_theta,n_bins))
cov_mat = cov_mat[:,:,mask]
cov_mat = cov_mat[mask]
cov_mat = cov_mat[:,:pars['n_kl'],:,:pars['n_kl']]
cov_mat = cov_mat.reshape((len(mask[mask])*pars['n_kl'],len(mask[mask])*pars['n_kl']))
inv_cov_mat = (pars['n_sim']-len(mask[mask])*pars['n_kl']-2.)/(pars['n_sim']-1.)*np.linalg.inv(cov_mat)


#Construct full array of cosmological parameters
def get_cosmo(cosmo_vars):
    final = cosmo_pars
    count1 = 0
    for count2, par in enumerate(final):
        if type(final[count2]) is str:
            final[count2] = cosmo_vars[count1]
            count1 = count1+1
    return final


#Get theory correlation function
def get_theory(cosmo_vars):
    #Get cosmological parameters
    cosmo_vars_tot = get_cosmo(cosmo_vars)
    print cosmo_vars_tot
    #Cosmology
    cosmo = ccl.Cosmology(h=cosmo_vars_tot[0], Omega_c=cosmo_vars_tot[1], Omega_b=cosmo_vars_tot[2], sigma8=cosmo_vars_tot[3], n_s=cosmo_vars_tot[4])
    #Tracers
    lens = np.array([ccl.ClTracerLensing(cosmo, False, z=z.astype(np.float64), n=pz[x].astype(np.float64)) for x in range(n_bins)])
    #Cl's
    ell = np.arange(n_ells)
    cls = np.zeros((n_bins, n_bins, n_ells))
    ccl.angular_cl(cosmo, lens[0], lens[0], ell)
    #for count1 in range(n_bins):
        #for count2 in range(n_bins):
            #cls[count1,count2] = ccl.angular_cl(cosmo, lens[count1], lens[count2], ell)
    #cls = np.transpose(cls,axes=[2,0,1])
    ##Correlation function
    #xi_th = np.zeros((2, n_bins, n_bins, n_theta))
    #for count1 in range(n_bins):
        #for count2 in range(n_bins):
            #for count3 in range(n_theta):
                #xi_th[0,count1,count2,count3] = ccl.correlation(cosmo, ell, cls[:,count1,count2], theta[count3], corr_type='L+', method='FFTLog')
                #xi_th[1,count1,count2,count3] = ccl.correlation(cosmo, ell, cls[:,count1,count2], theta[count3], corr_type='L-', method='FFTLog')
    #xi_th = np.transpose(xi_th,axes=[0,3,1,2])
    ##KL transform
    #xi_th = kl_t.dot(xi_th).dot(kl_t.T)
    #xi_th = np.transpose(xi_th, axes=[1, 2, 0, 3])
    #xi_th = np.diagonal(xi_th,  axis1=2, axis2=3)
    #xi_th = xi_th.reshape((2*n_theta,n_bins))
    #xi_th = xi_th[mask,:pars['n_kl']].flatten()

    #return xi_th


#Define priors
def lnprior(cosmo_vars):
    count1 = 0
    is_in = True
    for count2, par in enumerate(cosmo_pars):
        if type(cosmo_pars[count2]) is str:
            is_in = is_in & (pars[par][0] <= cosmo_vars[count1] <= pars[par][2])
            count1 = count1+1
    if is_in:
        return 0.0
    return -np.inf

#Define likelihood
def lnlike(cosmo_vars, xi_obs, inv_cov_mat):
    #Get theory
    xi_th = get_theory(cosmo_vars)
    #Get chi2
    chi2 = (xi_obs-xi_th).dot(inv_cov_mat).dot(xi_obs-xi_th)
    return -chi2

#Define posterior
def lnprob(cosmo_vars, xi_obs, inv_cov_mat):
    lp = lnprior(cosmo_vars)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(cosmo_vars, xi_obs, inv_cov_mat)


#Get random initial points
def get_random(mean, pars, cosmo_pars):
    count1 = 0
    rnd_pars = np.array([])
    for count2, par in enumerate(cosmo_pars):
        if type(par) is str:
            rnd = (np.random.random_sample()-0.5)/(pars[par][2]-pars[par][0])/1.e1
            rnd_pars = np.append(rnd_pars, mean[count1] + rnd)
            count1 = count1+1
    return rnd_pars

#Initial point
cosmo_vars_0 = np.array([get_random(cosmo_vars, pars, cosmo_pars) for x in range(n_walkers)])

#get_theory(np.array([0.54, 0.28, 0.59]))
get_theory(np.array([0.68, 0.34, 0.71]))

#print lnprob([0.53847335941, 0.280045328237, 0.585219830549], xi_obs, inv_cov_mat)
#print lnprob(np.array([0.694041927815, 0.295496912016, 0.599160329902]), xi_obs, inv_cov_mat)

##Initialize sampler
#sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, args=[xi_obs, inv_cov_mat])


##Creat file
#f = open(paths['output_dir'], 'w')
#f.close()

#for result in sampler.sample(cosmo_vars_0, iterations=n_steps, storechain=False):
    #pos = result[0]
    #prob = result[1]
    #f = open(paths['output_dir'], 'a')
    #for k in range(pos.shape[0]):
        #f.write(str(prob[k]) + "    " + "    ".join([str(x) for x in pos[k]]) + "\n")#"{0:4d} {1:s}\n".format(k, " ".join(position[k]))
    #f.close()

