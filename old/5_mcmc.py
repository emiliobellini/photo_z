import os, sys, re
import argparse
from astropy.io import fits
import numpy as np
import pyccl as ccl
import emcee



# Parse the given arguments
parser = argparse.ArgumentParser("Calculate the correlation function, apply the KL transform and calculate the chi^2")
parser.add_argument("--input_file", "-i", type=str, default = None, help="Input parameters file")
parser.add_argument("--data_file", "-d", type=str, default = None, help="Input data file")
parser.add_argument("--output_file", "-o", type=str, default = None, help="Output file")
parser.add_argument("--restart", "-r", help="Restart the chains from the last point of the output file", action="store_true")
args = parser.parse_args()
if not(args.input_file):
    raise IOError('You should specify an input file!')
if not(args.data_file):
    raise IOError('You should specify a data file!')


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)
paths['data'] = os.path.abspath(args.data_file)


#Default parameters [h, omega_c, omega_b, ln10_A_s, n_s]
cosmo_pars = np.array([
    [None, 0.61197750, None],
    [None, 0.11651890, None],
    [None, 0.03274485, None],
    [None, 2.47363700, None],
    [None, 1.25771300, None]
    ])
n_ells = 2001
n_walkers = 10
n_steps = 2
n_threads = 2
n_sim = 1988



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

def read_line(file_path, par, exp='floats'):
    with open(file_path) as fn:
        for line in fn:
            if '=' in line and line[0] != '#':
                line = re.sub('#.+', '', line)
                name , value = line.split('=')
                name = name.strip()
                if name == par:
                    value = value.strip()
                    if exp=='string':
                        return value
                    else:
                        value = listify(value)
                        if type(value) != str:
                            if value[0]>=value[2]:
                                raise IOError('The left bound is larger than the right bound for parameter ' + par)
                            if (value[0]>=value[1]) or (value[1]>=value[2]):
                                raise IOError('The central value is out of the bounds for parameter ' + par)
                        else:
                            value = np.array([None, floatify(value), None])
                        if len(value) != 3:
                            raise IOError('Incorrect length of input for parameter ' + par)
                        return value
        raise ValueError()

try:
    cosmo_pars[0] = read_line(paths['input'], 'h', exp='floats')
except ValueError:
    pass
try:
    cosmo_pars[1] = read_line(paths['input'], 'omega_c', exp='floats')
except ValueError:
    pass
try:
    cosmo_pars[2] = read_line(paths['input'], 'omega_b', exp='floats')
except ValueError:
    pass
try:
    cosmo_pars[3] = read_line(paths['input'], 'ln10_A_s', exp='floats')
except ValueError:
    pass
try:
    cosmo_pars[4] = read_line(paths['input'], 'n_s', exp='floats')
except ValueError:
    pass
try:
    n_ells = int(read_line(paths['input'], 'ell_max', exp='floats')[1]+1)
except ValueError:
    pass
try:
    n_walkers = int(read_line(paths['input'], 'n_walkers', exp='floats')[1])
except ValueError:
    pass
try:
    n_steps = int(read_line(paths['input'], 'n_steps', exp='floats')[1])
except ValueError:
    pass
try:
    n_sim = int(read_line(paths['input'], 'n_sim', exp='floats')[1])
except ValueError:
    pass
try:
    n_threads = int(read_line(paths['input'], 'n_threads', exp='floats')[1])
except ValueError:
    pass
try:
    method = read_line(paths['input'], 'method', exp='string')
    if method=='full':
        is_kl = False
        is_diag = False
    elif method=='kl_off_diag':
        is_kl = True
        is_diag = False
    elif method=='kl_diag':
        is_kl = True
        is_diag = True
    else:
        raise IOError('Method not recognized. Options are: full, kl_off_diag, kl_diag')
except ValueError:
    raise ValueError('Missing method in the ini file')
if is_kl:
    try:
        n_kl = int(read_line(paths['input'], 'n_kl', exp='floats')[1])
    except ValueError:
        raise ValueError('Missing n_kl in the ini file')



#Read data file
with fits.open(paths['data']) as fn:
    z = fn['photoz_z'].data
    pz = fn['photoz_p'].data
    theta = fn['theta'].data
    mask_theta = fn['mask_theta'].data.astype(bool)
    xi_obs = fn['xi_obs'].data
    cov_mat = fn['cov_mat'].data
    if is_kl:
        kl_t = fn['kl_t_avg'].data
n_bins = len(pz)
n_theta = len(theta)



#KL transform and reshape data
def kl_transform(data, datat='corr'):
    data_kl = kl_t.dot(data).dot(kl_t.T)
    if datat=='corr':
        return np.transpose(data_kl, axes=[1, 2, 0, 3])
    elif datat=='cov':
        data_kl = np.transpose(data_kl, axes=[5, 6, 0, 7, 1, 2, 3, 4])
        data_kl = kl_t.dot(data_kl).dot(kl_t.T)
        return np.transpose(data_kl, axes=[1, 2, 3, 4, 5, 6, 0, 7])


def reshape(data, datat='corr'):
    data_r = data
    if datat=='corr':
        data_r = data_r.reshape((2*n_theta,n_bins,n_bins))
        data_r = np.triu(data_r[mask_theta])
        if is_kl:
            data_r = data_r[:,:n_kl,:n_kl]
        if is_diag:
            data_r = np.diagonal(data_r,  axis1=1, axis2=2)
    elif datat=='cov':
        data_r = data_r.reshape((2*n_theta,n_bins,n_bins,2*n_theta,n_bins,n_bins))
        data_r = np.triu(data_r[:,:,:,mask_theta])
        data_r = np.transpose(data_r,axes=[3,4,5,0,1,2])
        data_r = np.triu(data_r[:,:,:,mask_theta])
        if is_kl:
            data_r = data_r[:,:n_kl,:n_kl,:,:n_kl,:n_kl]
        if is_diag:
            data_r = np.diagonal(data_r,  axis1=4, axis2=5)
            data_r = np.diagonal(data_r,  axis1=1, axis2=2)
            data_r = np.transpose(data_r, axes=[0,2,1,3])
    data_r = data_r.flatten()
    data_r = data_r[data_r != 0]
    if datat=='cov':
        global n_data
        n_data = int(np.sqrt(len(data_r)))
        data_r = data_r.reshape((n_data,n_data))
    return data_r

if is_kl:
    xi_obs = kl_transform(xi_obs, datat='corr')
    cov_mat = kl_transform(cov_mat, datat='cov')
xi_obs = reshape(xi_obs, datat='corr')
cov_mat = reshape(cov_mat, datat='cov')
inv_cov_mat = (n_sim-n_data-2.)/(n_sim-1.)*np.linalg.inv(cov_mat)



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
    cosmo = ccl.Cosmology(h=var_tot[0], Omega_c=var_tot[1]/var_tot[0]**2., Omega_b=var_tot[2]/var_tot[0]**2., A_s=(10.**(-10.))*np.exp(var_tot[3]), n_s=var_tot[4])
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
    #Reshape and eventually KL transform
    if is_kl:
        xi_th = kl_transform(xi_th, datat='corr')
    xi_th = reshape(xi_th, datat='corr')
    return xi_th


#Get sigma_8
def get_sigma_8(var):
    #Get cosmological parameters
    var_tot = get_cosmo(var)
    #Cosmology
    cosmo = ccl.Cosmology(h=var_tot[0], Omega_c=var_tot[1]/var_tot[0]**2., Omega_b=var_tot[2]/var_tot[0]**2., A_s=(10.**(-10.))*np.exp(var_tot[3]), n_s=var_tot[4])
    sigma8 = ccl.sigma8(cosmo)
    return sigma8


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
    try:
        xi_th = get_theory(var)
    except:
        print 'CCL failure with pars = ' + str(var)
        sys.stdout.flush()
        return -np.inf
    #Get chi2
    chi2 = (obs-xi_th).dot(icov).dot(obs-xi_th)
    return -chi2/2.


#Define posterior
def lnprob(var):
    lp = lnprior(var)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(var)


if n_dim==0:
    print 'The number of varying parameters is less than 2, i.e. ' + str(n_dim)
    print 'Only the likelihood at the initial point will be evaluated'
    
    print 'Cosmological parameters:'
    print '----> h             = ' + '{0:2.4e}'.format(cosmo_pars[0,1])
    print '----> Omega_c h^2   = ' + '{0:2.4e}'.format(cosmo_pars[1,1])
    print '----> Omega_b h^2   = ' + '{0:2.4e}'.format(cosmo_pars[2,1])
    print '----> ln(10^10 A_s) = ' + '{0:2.4e}'.format(cosmo_pars[3,1])
    print '----> n_s           = ' + '{0:2.4e}'.format(cosmo_pars[4,1])
    print 'Derived parameters:'
    print '----> sigma_8       = ' + '{0:2.4e}'.format(get_sigma_8([]))
    print 'Likelihood:'
    print '----> -ln(like)     = ' + '{0:4.4f}'.format(-lnprob([]))
    
elif n_dim==1:
    raise IOError('To run a chain you should specify at least two varying parameters!')
else:
    if not(args.output_file):
        raise IOError('You should specify an output file')
    paths['output'] = os.path.abspath(args.output_file)


    #Initialize sampler
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob, threads=n_threads)


    #Print useful stuff
    print 'Starting the chains!'
    print 'Number of threads = ' + str(n_threads)
    print 'Number of steps = ' + str(n_steps)
    print 'Number of walkers = ' + str(n_walkers)
    print 'Maximum ell = ' + str(n_ells-1)
    print 'Method = ' + method
    if is_kl:
        print 'Number of KL modes = ' + str(n_kl)
    sys.stdout.flush()


    if args.restart:
        #Initial point from data
        vars_0 = np.loadtxt(paths['output'],unpack=True)
        vars_0 = vars_0[2:2+n_dim]
        vars_0 = vars_0[:,-n_walkers:].T
    else:
        #Initial point
        vars_0 = np.array([get_random(cosmo_pars[mask_vars], 1.e1) for x in range(n_walkers)])
        #Create file
        f = open(paths['output'], 'w')
        f.close()

    for count, result in enumerate(sampler.sample(vars_0, iterations=n_steps, storechain=False)):
        pos = result[0]
        prob = result[1]
        f = open(paths['output'], 'a')
        for k in range(pos.shape[0]):
            out = np.append(np.array([1., -prob[k]]), pos[k])
            out = np.append(out, get_sigma_8(pos[k]))
            f.write('    '.join(['{0:.10e}'.format(x) for x in out]) + '\n')
        f.close()
        if (count+1) % 10 == 0:
            print '----> Computed ' + '{0:5.1%}'.format(float(count+1) / n_steps) + ' of the steps'
            sys.stdout.flush()



print 'Success!!'
