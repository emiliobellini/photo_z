import os, sys
import argparse
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import math
import variables as vrs



# Parse the given arguments
parser = argparse.ArgumentParser("Calculate the KL transform")
parser.add_argument("input_file", type=str, help="Input FITS file")
args = parser.parse_args()


#Define absolute paths to input and output files
paths = {}
paths['input'] = os.path.abspath(args.input_file)


#Open files and store datas
#Read image with the Cl's (S, signal)
with fits.open(paths['input']) as fn:
    S = fn['Cls'].data
    n_eff = fn['n_eff'].data*(180.*60./math.pi)**2. #converted in stedrad^-1
    sigma_g = fn['sigma_g'].data
#Check that the input quantities have the right dimensions
if (len(sigma_g) != len(n_eff)) or (S.shape[1] != len(n_eff)) or (S.shape[2] != len(n_eff)):
    raise IOError('Input arrays have disagreeing shapes. n_eff.shape = ' + str(n_eff.shape) +
        '. sigma_g.shape = ' + str(sigma_g.shape) + '. CLS.shape = ' + str(S.shape) + '.')


#Calculate noise (N), and Cholesky decompose it (N=LL^+)
N = np.array([np.diag(sigma_g**2/n_eff) for x in range(len(S))])
L = np.linalg.cholesky(N)
inv_L = np.linalg.inv(L)

#Calculate matrix for which we want to calculate eigenvalues and eigenvectors
M = np.array([np.dot(inv_L[x],S[x]+N[x]) for x in range(len(S))])
M = np.array([np.dot(M[x],inv_L[x].T) for x in range(len(S))])

#Calculate eigenvalues and eigenvectors
eigval, eigvec = np.linalg.eigh(M)
#Re-order eigenvalues and eigenvectors sorted from smallest to largest eigenvalue
new_ord = np.array([np.argsort(eigval[x])[::-1] for x in range(len(S))])
eigval = np.array([eigval[x][new_ord[x]] for x in range(len(S))])
eigvec = np.array([eigvec[x][:,new_ord[x]] for x in range(len(S))])

#Calculate transformation matrix (E) from eigenvectors and L^-1
E = np.array([np.dot(eigvec[x].T,inv_L[x]) for x in range(len(S))])
#Change sign to eigenvectors according to the first element
signs = np.array([[np.sign(E[ell][x][0]/E[2][x][0]) for x in range(len(S[0]))] for ell in range(len(S))])
E = np.array([(E[x].T*signs[x]).T for x in range(len(S))])

#Calculate the KL transformed Cl's
angular_cl = np.array([np.diag(eigval[x]) for x in range(len(S))])

#Test if the transformation matrix gives the correct new Cl's
test1 = np.array([np.dot(E[x],S[x]+N[x]) for x in range(len(S))])
test1 = np.array([np.dot(test1[x],E[x].T) for x in range(len(S))])
test1 = np.array([abs(test1[x]-angular_cl[x]) for x in range(len(S))])
test1 = test1[1:].max()/abs(angular_cl[1:]).max()
test2 = np.array([np.dot(L[x].T,E[x].T) for x in range(len(S))])
test2 = np.array([np.dot(test2[x],E[x]) for x in range(len(S))])
test2 = np.array([np.dot(test2[x],L[x]) for x in range(len(S))])
test2 = np.array([abs(test2[x]-np.identity(len(range(len(S[0]))))) for x in range(len(S))])
test2 = test2[1:].max()
if test1>vrs.MP_THRESHOLD or test2>vrs.MP_THRESHOLD:
    print('WARNING: the transformation matrix does not reproduce the correct Cl\'s.'
        + ' The relative difference is ' + '{:1.2e}'.format(max(test1,test2))
        + ' and the maximum accepted is ' + '{:1.2e}'.format(vrs.MP_THRESHOLD) + '.')
    sys.stdout.flush()
else:
    print 'Calculated new Cl\'s and KL transformation matrix'
    sys.stdout.flush()


#Average the kl transform to be l independent
E_avg = np.zeros((len(E[0]),len(E[0])))
den = np.array([(2.*x+1) for x in range(2,len(E))]).sum()
for n in range(len(E[0])):
    for m in range(len(E[0])):
        num = np.array([(2.*x+1)*E[:,n][:,m][x] for x in range(2,len(E))]).sum()
        E_avg[n][m] = num/den



#Update input file
hdu = {}
hdu['NOISE'] = fits.ImageHDU(N, name='NOISE')
hdu['ClS_KL'] = fits.ImageHDU(angular_cl, name='ClS_KL')
hdu['KL_T'] = fits.ImageHDU(E, name='KL_T')
hdu['KL_T_AVG'] = fits.ImageHDU(E_avg, name='KL_T_AVG')
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
