import os
import re
import numpy as np

from astropy.io import fits

import settings


def file_exists_or_error(fname):
    abspath = os.path.abspath(fname)
    if os.path.exists(abspath):
        return abspath
    raise IOError('File '+abspath+' not found!')

def listify(data):
    try:
        return np.array([float(data), float(data), float(data)])
    except:
        try:
            array = data.split(',')
            array = [x.strip() for x in array]
            return np.array([None if x=='None' else float(x) for x in array])
        except:
            return data

def read_param_from_file(fname, par):
    with open(fname) as fn:
        for line in fn:
            line = re.sub('#.+', '', line)
            if '=' in line:
                name , value = line.split('=')
                name = name.strip()
                value = value.strip()
                if name == par:
                    return value
    return None

def get_param(fname, par, type='string'):
    value = read_param_from_file(fname, par)
    if value is None:
        value = settings.default_params[par]
    if type=='cosmo':
        return listify(value)
    elif type=='float':
        return float(value)
    elif type=='int':
        return int(value)
    elif type=='path':
        return os.path.abspath(value)
    elif type=='bool':
        if value=='yes':
            return True
        elif value=='no':
            return False
        else:
            raise IOError('Parameter value for ' +par+ ' not recognized!')
    elif type=='string':
        return value
    else:
        raise IOError('Parameter type for ' +par+ ' not recognized!')

def get_cosmo_array(fname, pars):
    cosmo_pars = []
    for n, par in enumerate(pars):
        value = get_param(fname, par, type='cosmo')
        if value is None:
            cosmo_pars.append(settings.default_params[par])
        elif (len(value)==3) and (type(value) is not str):
            cosmo_pars.append(value)
        else:
            raise IOError('Check the syntax of the cosmological parameters!')
    return np.array(cosmo_pars)

def folder_exists_or_create(fname):
    abspath = os.path.abspath(fname)
    folder, _ = os.path.split(abspath)
    if not os.path.exists(folder):
        os.makedirs(folder)
    return abspath

def sanity_checks(cosmo_params, cosmo_list, settings, path):
    def check(test, msg):
        if test:
            raise IOError(msg)
    #Cosmo checks
    for n, par in enumerate(cosmo_params):
        check(par[1] is None,
            'Central value for ' +cosmo_list[n]+ ' is None!')
        check((par[0] is not None) and (par[0]>par[1]),
            'Central value for ' +cosmo_list[n]+ ' is lower than the left bound!')
        check((par[2] is not None) and (par[1]>par[2]),
            'Central value for ' +cosmo_list[n]+ ' is larger than the right bound!')
    #Check strings
    check(settings['method'] not in ['full', 'kl_off_diag', 'kl_diag'],
        'method not recognized! Options: full, kl_off_diag, kl_diag')
    check(settings['sampler'] not in ['emcee', 'fisher', 'single_point'],
        'sampler not recognized! Options: emcee, fisher, single_point')
    check(settings['space'] not in ['real', 'fourier'],
        'space not recognized! Options: real, fourier')
    if settings['n_sims'] not in ['auto', 'all']:
        try:
            if int(settings['n_sims'])<1:
                raise IOError('If n_sims is an integer, it has to be strictly positive!')
        except:
            raise IOError('n_sims not recognized! Options: auto, all, a positive integer')
    #Check numbers
    check(settings['n_threads']<1,
        'n_threads should be at least 1!')
    check(settings['ell_max']<2,
        'ell_max too little!')
    #Consistency between parameters
    mask = mask_cosmo_params(cosmo_params)
    check(settings['sampler'] in ['emcee', 'fisher'] and len(mask[mask]) < 2,
        'For emcee and fisher the minimum number of varied parameters is 2!')
    check(settings['sampler']=='emcee' and  settings['n_walkers'] < 1,
        'For emcee the minimum number of walkers is 1!')
    check(settings['sampler']=='emcee' and  settings['n_steps'] < 1,
        'For emcee the minimum number of steps is 1!')
    check(settings['method'] in ['kl_off_diag', 'kl_diag'] and
        settings['n_kl'] < 1,
        'n_kl should be at least 1!')
    #Check data existence
    with fits.open(path['data']) as hdul:
        imgs = [hdul[n].name for n in range(1,len(hdul))]
        for name in ['PHOTO_Z', 'N_EFF', 'SIGMA_G']:
            check(name not in imgs,
                name + ' was not found in data file!')
        if settings['space']=='real':
            for name in ['THETA', 'MASK_THETA', 'XIPM_OBS', 'XIPM_SIM']:
                check(name not in imgs,
                    name + ' was not found in data file!')
        if settings['space']=='fourier':
            raise ValueError("Implement conditions for Fourier space")
            for name in ['ELL', 'CLS_OBS', 'CLS_SIM']:
                check(name not in imgs,
                    name + ' was not found in data file!')
        #Check data dimensions
        n_bins = hdul['PHOTO_Z'].shape[0]-1
        check(len(hdul['PHOTO_Z'].shape)!=2,
            'PHOTO_Z has wrong dimensions!')
        check(hdul['N_EFF'].shape!=(n_bins,),
            'N_EFF has wrong dimensions!')
        check(hdul['SIGMA_G'].shape!=(n_bins,),
            'SIGMA_G has wrong dimensions!')
        if settings['space']=='real':
            n_theta = hdul['THETA'].shape[0]
            n_sims = hdul['XIPM_SIM'].shape[0]
            check(len(hdul['THETA'].shape)!=1,
                'THETA has wrong dimensions!')
            check(hdul['MASK_THETA'].shape!=(2,n_theta),
                'MASK_THETA has wrong dimensions!')
            check(hdul['XIPM_OBS'].shape!=(2,n_theta,n_bins,n_bins),
                'XIPM_OBS has wrong dimensions!')
            check(hdul['XIPM_SIM'].shape!=(n_sims,2,n_theta,n_bins,n_bins),
                'XIPM_SIM has wrong dimensions!')
        elif settings['space']=='fourier':
            raise ValueError("Implement conditions for Fourier space")
    return

def mask_cosmo_params(params):
    def mask_line(param):
        if (param[0] is None) or (param[2] is None):
            return True
        if (param[0]<param[1]) or (param[1]<param[2]):
            return True
        return False
    return np.array([mask_line(x) for x in params])

def how_many_sims(settings, data):
    all_sim = len(data['corr_sim'])
    n_x_var = data['mask_x_var'].flatten()
    n_x_var = len(n_x_var[n_x_var])
    n_bins = len(data['photo_z'])-1
    all_points = n_x_var*n_bins*(n_bins+1)/2
    ratio = (all_sim-all_points-2.)/(all_sim-1.)
    if settings['n_sims']=='all':
        return all_sim
    elif settings['n_sims']=='auto':
        n_kl = settings['n_kl']
        if settings['method']=='kl_off_diag':
            all_points = n_x_var*n_kl*(n_kl+1)/2
        elif settings['method']=='kl_diag':
            all_points = n_x_var*n_kl
        return int(round((2.+all_points-ratio)/(1.-ratio)))
    else:
        return int(settings['n_sims'])

def compute_covmat(data, settings):
    return None
