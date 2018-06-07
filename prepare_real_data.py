import argparse
import os
import numpy as np

import settings
import tools


#Parse the given arguments
parser = argparse.ArgumentParser('Create data FITS file in real space')
parser.add_argument('input_folder', type=str, help='Input folder')
args = parser.parse_args()


#Define absolute paths and check the existence of each required file
path = {
    # 'data' : tools.file_exist_or_error(args.input_folder + 'data.fits'),
    'xipm' : tools.file_exist_or_error(args.input_folder + 'xipm.dat'),
    'sims' : tools.file_exist_or_error(args.input_folder + 'mockxipm.tgz'),
    'output' : tools.file_exist_or_error(os.path.abspath('') + '/output') + '/data_real.fits'
    }


#Angular separation (theta)
theta = np.array(settings.THETA_ARCMIN)/60. # Theta is in degrees
tools.write_to_fits(fname=path['output'], array=theta, name='theta')
#Mask for theta
mask_theta = np.array(settings.MASK_THETA).astype(int)
tools.write_to_fits(fname=path['output'], array=mask_theta, name='mask_theta')

print(tools.position_xipm(6))
#Print info about the fits file
tools.print_info_fits(fname=path['output'])
