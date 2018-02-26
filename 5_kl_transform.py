import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import variables as vrs



# Parse the given arguments
parser = argparse.ArgumentParser("Calculate the KL transform")
parser.add_argument("input_file", type=str, help="Input FITS file")
parser.add_argument("--output_plot", "-p", type=str, default = None, help="Output plot")
parser.add_argument("--skip_plots", "-s", help="Skip plots", action="store_true")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)
if args.output_plot:
    paths['output_plot'] = os.path.abspath(args.output_plot)
else:
    paths['output_plot'] = os.path.splitext(paths['input'])[0]
    paths['output_plot'] = paths['output_plot'] + '.pdf'


#Open files and store datas
#Read image with the Cl's (S, signal)
with fits.open(paths['input']) as fn:
    S = fn['Cls'].data
    z_bins = fn['z_bins'].data


#Calculate n_eff and sigma_g for each bin
#TODO: insert correct values for n_eff (A_eff) and sigma_g
n_eff = np.array([10., 11., 9., 12., 8., 13., 7.]) #n_eff in arcmin^-2
n_eff = n_eff*3282.8*60.*60. #n_eff in srad^-1
sigma_g = 0.28
diag_noise = sigma_g**2/n_eff

#Calculate noise (N), and Cholesky decompose it (N=LL^+)
N = np.array([np.diag(diag_noise) for x in vrs.L_RANGE])
L = np.linalg.cholesky(N)
inv_L = np.linalg.inv(L)

#Calculate matrix for which we want to calculate eigenvalues and eigenvectors
M = np.array([np.dot(inv_L[x],S[x]+N[x]) for x in vrs.L_RANGE])
M = np.array([np.dot(M[x],inv_L[x].T) for x in vrs.L_RANGE])

#Calculate eigenvalues and eigenvectors
eigval, eigvec = np.linalg.eigh(M)
#Re-order eigenvalues and eigenvectors sorted from smallest to largest eigenvalue
new_ord = np.array([np.argsort(eigval[x])[::-1] for x in vrs.L_RANGE])
eigval = np.array([eigval[x][new_ord[x]] for x in vrs.L_RANGE])
eigvec = np.array([eigvec[x][:,new_ord[x]] for x in vrs.L_RANGE])

#Calculate transformation matrix (E) from eigenvectors and L^-1
E = np.array([np.dot(eigvec[x].T,inv_L[x]) for x in vrs.L_RANGE])
#Change sign to eigenvectors according to the first element
signs = np.array([[np.sign(E[ell][x][0]/E[2][x][0]) for x in vrs.BINS_RANGE] for ell in vrs.L_RANGE])
E = np.array([(E[x].T*signs[x]).T for x in vrs.L_RANGE])


#Calculate the KL transformed Cl's
angular_cl = np.array([np.diag(eigval[x]) for x in vrs.L_RANGE])

#Test if the transformation matrix gives the correct new Cl's
test = np.array([np.dot(E[x],S[x]+N[x]) for x in vrs.L_RANGE])
test = np.array([np.dot(test[x],E[x].T) for x in vrs.L_RANGE])
test = np.array([abs(angular_cl[x]-test[x]) for x in vrs.L_RANGE])
test = test[1:].max()/abs(angular_cl[1:]).max()
if test>vrs.MP_THRESHOLD:
    print('WARNING: the transformation matrix does not reproduce the correct Cl\'s.'
        + ' The relative difference is ' + '{:1.2e}'.format(test)
        + ' and the maximum accepted is ' + '{:1.2e}'.format(vrs.MP_THRESHOLD) + '.')
    sys.stdout.flush()
else:
    print 'Calculated new Cl\'s and KL transformation matrix'
    sys.stdout.flush()



#Update input file with an image containing Cl's after KL and the transformation matrix
#Create image for Cl's
hdu = {}
hdu['Cls_KL'] = fits.ImageHDU(angular_cl, name='Cls_KL')
hdu['KL_T'] = fits.ImageHDU(E, name='KL_T')
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


#Plots
if not args.skip_plots:
    #Plot KL Cl's
    cm = plt.get_cmap('Blues')
    for count in vrs.BINS_RANGE:
        y = [angular_cl[x][count][count] for x in vrs.L_RANGE]
        plt.plot(vrs.L_RANGE, y, label=r'$\alpha_'+str(count+1)+'$',
        color=cm(0.2+0.6*(1.-float(count)/float(len(angular_cl[0])))))
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$D_\ell^{\alpha\alpha}$')
    plt.title(r'Eigenvalues')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(vrs.L_MIN, vrs.L_MAX)
    plt.legend()
    plot_path = os.path.splitext(paths['output_plot'])
    plot_path = plot_path[0]+'_eigenvalues'+plot_path[1]
    plt.savefig(plot_path)
    plt.close()
    print 'Saved KL eigenvalues plot at ' + os.path.relpath(plot_path)
    sys.stdout.flush()

    #Plot eigenvectors
    for count in [0,1,2]:
        for ell in vrs.L_RANGE[2:5]:
            x = np.array([np.average(z_bins[x]) for x in vrs.BINS_RANGE])
            plt.plot(x, E[ell][count]/np.linalg.norm(E[ell][count]))
    plot_path = os.path.splitext(paths['output_plot'])
    plot_path = plot_path[0]+'_eigenvectors'+plot_path[1]
    plt.savefig(plot_path)
    plt.close()


print 'Success!!'
