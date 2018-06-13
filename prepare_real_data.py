import argparse
import os
import numpy as np

import settings
import tools_prep as tools


#Parse the given arguments
parser = argparse.ArgumentParser('Create data FITS file in real space')
parser.add_argument('input_folder', type=str, help='Input folder')
args = parser.parse_args()


#Define absolute paths and check the existence of each required file
path = {
    'data' : tools.file_exists_or_error(args.input_folder + 'data.fits'),
    'xipm' : tools.file_exists_or_error(args.input_folder + 'xipm.dat'),
    'sims' : tools.file_exists_or_error(args.input_folder + 'mockxipm.tar.gz'),
    'output' : tools.file_exists_or_error(os.path.abspath('') + '/data') + '/data_real.fits'
    }


#Angular separation (theta)
theta = np.array(settings.THETA_ARCMIN)/60. # Theta is in degrees
tools.write_to_fits(fname=path['output'], array=theta, name='theta')
#Mask for theta
mask_theta = np.array(settings.MASK_THETA).astype(int)
tools.write_to_fits(fname=path['output'], array=mask_theta, name='mask_theta')


#Read and reshape xipm observed
xipm = np.loadtxt(path['xipm'], dtype='float64')
xipm = tools.unflatten_xipm(xipm[:,1])
tools.write_to_fits(fname=path['output'], array=xipm, name='xipm_obs')


#Read and reshape xipm from simulations
xipm = tools.unpack_and_stack(fname=path['sims'])
xipm = np.array([tools.unflatten_xipm(x) for x in xipm])
tools.write_to_fits(fname=path['output'], array=xipm, name='xipm_sim')


#Calculate photo-z sigma_g and n_eff
photo_z, n_eff, sigma_g = tools.read_fits_data(path['data'])
tools.write_to_fits(fname=path['output'], array=photo_z, name='photo_z')
tools.write_to_fits(fname=path['output'], array=n_eff, name='n_eff')
tools.write_to_fits(fname=path['output'], array=sigma_g, name='sigma_g')


#Print info about the fits file
tools.print_info_fits(fname=path['output'])


print('Success!!')
