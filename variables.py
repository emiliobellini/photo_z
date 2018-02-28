#List and values of all the variables and functions used in the code


#Photo-z Bins (minimum, maximum and intermediate bins)
Z_BINS = [0.15,0.29,0.43,0.57,0.70,0.90,1.10,1.30]


#CFHTlens specifications
CFHTlens_dZ = 0.05
CFHTlens_A_eff = {}
CFHTlens_A_eff['W1'] = (60.**2)*sum([0.76, 0.89, 0.85, 0.72, 0.79, 0.78, 0.85,
                                     0.81, 0.84, 0.89, 0.83, 0.84, 0.78, 0.79,
                                     0.83, 0.87, 0.87, 0.79, 0.84, 0.84, 0.83,
                                     0.82, 0.86, 0.82, 0.84, 0.80, 0.79, 0.83,
                                     0.80, 0.85, 0.75, 0.85, 0.81, 0.83, 0.79,
                                     0.85, 0.83, 0.76, 0.78, 0.82, 0.75, 0.90,
                                     0.81, 0.85, 0.84, 0.82, 0.83, 0.76, 0.76,
                                     0.86, 0.87, 0.85, 0.79, 0.83, 0.80]) #arcmin^2
CFHTlens_A_eff['W2'] = (60.**2)*sum([0.65, 0.63, 0.63, 0.58, 0.68, 0.61, 0.62,
                                     0.70, 0.72, 0.62, 0.66, 0.63, 0.65, 0.62,
                                     0.63, 0.58, 0.72, 0.70, 0.76]) #arcmin^2
CFHTlens_A_eff['W3'] = (60.**2)*sum([0.79, 0.79, 0.68, 0.72, 0.75, 0.78, 0.81,
                                     0.79, 0.84, 0.76, 0.81, 0.77, 0.81, 0.77,
                                     0.79, 0.66, 0.77, 0.75, 0.80, 0.79, 0.78,
                                     0.74, 0.79, 0.71, 0.79, 0.69, 0.76, 0.79,
                                     0.82, 0.84, 0.76, 0.79, 0.76, 0.75, 0.80,
                                     0.75]) #arcmin^2
CFHTlens_A_eff['W4'] = (60.**2)*sum([0.73, 0.76, 0.72, 0.70, 0.77, 0.77, 0.61,
                                     0.62, 0.71, 0.64, 0.66, 0.67, 0.72, 0.74,
                                     0.73, 0.72, 0.71, 0.65, 0.73]) #arcmin^2
CFHTlens_A_eff['TOT'] = sum(CFHTlens_A_eff.values()) #arcmin^2


#Define ell_min and ell_max
L_MIN = 2
L_MAX = 2000


#Cosmology
Omega_c=0.27
Omega_b=0.045
h=0.67
sigma8=0.8
#A_s=2.215e-9
n_s=0.96


#Define Machine Precision Threshold, above which a warning is raised
MP_THRESHOLD = 1.e-14


#Criteria used to select the data
def select_data(data, z_min, z_max):

    sel = data['Z_B']>=z_min
    sel = (data['Z_B']<z_max)*sel
    sel = (data['MASK']<=1)*sel
    sel = ([x in good_fit_patterns for x in data['field']])*sel

    return sel


#Good fit patterns
good_fit_patterns = ['W1m0m0', 'W1m0m3', 'W1m0m4', 'W1m0p1', 'W1m0p2', 'W1m0p3', 'W1m1m0',
                     'W1m1m2', 'W1m1m3', 'W1m1m4', 'W1m1p3', 'W1m2m1', 'W1m2m2', 'W1m2m3',
                     'W1m2p1', 'W1m2p2', 'W1m3m0', 'W1m3m2', 'W1m3m4', 'W1m3p1', 'W1m3p3',
                     'W1m4m0', 'W1m4m1', 'W1m4m3', 'W1m4m4', 'W1m4p1', 'W1p1m1', 'W1p1m2',
                     'W1p1m3', 'W1p1m4', 'W1p1p1', 'W1p1p2', 'W1p1p3', 'W1p2m0', 'W1p2m2',
                     'W1p2m3', 'W1p2m4', 'W1p2p1', 'W1p2p2', 'W1p2p3', 'W1p3m1', 'W1p3m2',
                     'W1p3m3', 'W1p3m4', 'W1p3p1', 'W1p3p2', 'W1p3p3', 'W1p4m0', 'W1p4m1',
                     'W1p4m2', 'W1p4m3', 'W1p4m4', 'W1p4p1', 'W1p4p2', 'W1p4p3',

                     'W2m0m0', 'W2m0m1', 'W2m0p1', 'W2m0p2', 'W2m1m0', 'W2m1m1', 'W2m1p1',
                     'W2m1p3', 'W2p1m0', 'W2p1p1', 'W2p1p2', 'W2p2m0', 'W2p2m1', 'W2p2p1',
                     'W2p2p2', 'W2p3m0', 'W2p3m1', 'W2p3p1', 'W2p3p3',

                     'W3m0m1', 'W3m0m2', 'W3m0m3', 'W3m0p2', 'W3m0p3', 'W3m1m0', 'W3m1m2',
                     'W3m1m3', 'W3m1p1', 'W3m1p2', 'W3m1p3', 'W3m2m1', 'W3m2m2', 'W3m2m3',
                     'W3m2p1', 'W3m2p2', 'W3m3m0', 'W3m3m1', 'W3m3m2', 'W3m3m3', 'W3m3p1',
                     'W3m3p2', 'W3p1m0', 'W3p1m1', 'W3p1m2', 'W3p1m3', 'W3p1p2', 'W3p1p3',
                     'W3p2m0', 'W3p2m3', 'W3p2p3', 'W3p3m1', 'W3p3m3', 'W3p3p1', 'W3p3p2',
                     'W3p3p3',

                     'W4m0m2', 'W4m0p1', 'W4m1m0', 'W4m1m1', 'W4m1m2', 'W4m1p1', 'W4m2m0',
                     'W4m2p1', 'W4m2p3', 'W4m3m0', 'W4m3p1', 'W4m3p2', 'W4m3p3', 'W4p1m0',
                     'W4p1m1', 'W4p1m2', 'W4p2m0', 'W4p2m1', 'W4p2m2']
