# Initialize
# - Parse arguments
# - Read parameters
# - Read variables
# - Sanity checks


# Preliminary calculations
# - Compute how many simulations
# - Compute KL transform (conditional)
# - Apply KL transform
# - Compute inverse cov_mat

# Run (input: array with cosmo_params, kl_t, corr_obs, inv_cov_mat)



# import argparse
#
# from astropy.io import fits
#
# import tools_run as tools
# import settings
#
#
# #Parse the given arguments
# parser = argparse.ArgumentParser('Run MCMC, Fisher matrix or single point likelihood evaluation')
# parser.add_argument("--params_file", "-p", type=str, default = None, help="Input parameters file")
# parser.add_argument("--restart", "-r", help="Restart the chains from the last point of the output file (only for emcee)", action="store_true")
# args = parser.parse_args()
# if not(args.params_file):
#     raise IOError('You should specify a params file!')
#
#
# #Define absolute paths and check the existence of each required file
# path = {
#     'params' : tools.file_exists_or_error(args.params_file),
#     }
# path['data'] = tools.get_param(path['params'], 'data_file', type='path')
# tools.file_exists_or_error(path['data'])
# path['output'] = tools.get_param(path['params'], 'output', type='path')
# tools.folder_exists_or_create(path['output'])
#
# #Create array with cosmo parameters [h, omega_c, omega_b, ln10_A_s, n_s]
# cosmo_params_list = ['h', 'omega_c', 'omega_b', 'ln10_A_s', 'n_s']
# cosmo_params = tools.get_cosmo_array(path['params'], cosmo_params_list)
#
#
# #Read and store the other parameters
# settings = {
#     'n_walkers' : tools.get_param(path['params'], 'n_walkers', type='int'),
#     'n_steps' : tools.get_param(path['params'], 'n_steps', type='int'),
#     'n_threads' : tools.get_param(path['params'], 'n_threads', type='int'),
#     'ell_max' : tools.get_param(path['params'], 'ell_max', type='int'),
#     'n_kl' : tools.get_param(path['params'], 'n_kl', type='int'),
#     'method' : tools.get_param(path['params'], 'method', type='string'),
#     'sampler' : tools.get_param(path['params'], 'sampler', type='string'),
#     'space' : tools.get_param(path['params'], 'space', type='string'),
#     'n_sims' : tools.get_param(path['params'], 'n_sims', type='string'),
#     'kl_scale_dependence' : tools.get_param(path['params'], 'kl_scale_dependence', type='bool')
# }
#
#
# #Sanity checks
# tools.sanity_checks(cosmo_params, cosmo_params_list, settings, path)
#
#
# #Read data file
# data = {}
# with fits.open(path['data']) as fn:
#     data['photo_z'] = fn['photo_z'].data
#     data['n_eff'] = fn['n_eff'].data
#     data['sigma_g'] = fn['sigma_g'].data
#     if settings['space']=='real':
#         data['x_var'] = fn['theta'].data
#         data['mask_x_var'] = fn['mask_theta'].data.astype(bool)
#         data['corr_obs'] = fn['xipm_obs'].data
#         data['corr_sim'] = fn['xipm_sim'].data
#     elif settings['space']=='fourier':
#         data['x_var'] = fn['ell'].data
#         data['mask_x_var'] = fn['mask_ell'].data.astype(bool)
#         data['corr_obs'] = fn['cls_obs'].data
#         data['corr_sim'] = fn['cls_sim'].data
#
#
# #Compute how many simulations are used
# settings['n_sims'] = tools.how_many_sims(settings, data)
#
#
# #Compute covariance matrix
# data['cov_mat'] = tools.compute_covmat(data, settings)
#
# print settings
