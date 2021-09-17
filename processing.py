import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import math
import scipy
import emcee
import corner
import lalsimulation as ls
import h5py
import snr_calculation as s
import populations as p


SMALL_SIZE = 12
MEDIUM_SIZE = 15
BIGGER_SIZE = 18
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

def get_quantiles(fname, index, quantiles, hist=False, folder='mTOV_convergence', columns=None):
    root = '/mnt/c/users/christine/gwcosmology/spin_results/' + folder
    # print(root)
    array = np.genfromtxt(root+fname)

    if columns is not None:
        array[:,0] = array[:, columns[0]] - array[:, columns[1]]
        index = 0
    array = array[:,index]
    if hist:
        plt.hist(array,bins=50)
    return np.quantile(array, quantiles)


def quantiles_from_data(array, quantiles, hist=False):
    if hist:
        plt.hist(array,bins=50)
    return np.quantile(array, quantiles)

def get_constraints(data, name, folder='real_data'):
    constraints = np.zeros(7)
    sigma_1 = quantiles_from_data(data, [0.159, 0.841])
    sigma_2 = quantiles_from_data(data, [0.025, 0.975])
    sigma_3 = quantiles_from_data(data, [0.0015, 0.997])

    med = quantiles_from_data(data, [0.5])

    constraints = [sigma_3[0], sigma_2[0], sigma_1[0], med[0], sigma_1[1], sigma_2[1], sigma_3[1]]

    #print('../spin_results/outputs/{}/{}'.format(folder, name))
    np.savetxt('../spin_results/outputs/{}/{}'.format(folder, name), constraints)
    return constraints

import os

def processing(folder, mode, plot=True):
    base = '../spin_results/{}/'.format(folder)


    for i in os.listdir(base):
        if mode == 'LMG_txt':
            if 'likes' not in i and os.path.isfile(os.path.join(base,i)):
                root = i.split('.txt')[0]
                real = np.genfromtxt(base+'{}.txt'.format(root))
                likes = np.genfromtxt(base+'{}_likes.txt'.format(root))

                if '_u' in root:
                    print(root)
                    massgap = real[:,3]-real[:,2]
                    get_constraints(massgap, '{}.txt'.format(root), folder=folder)
                    print(np.median(real[:,2]), np.median(real[:,3]), np.median(massgap), massgap[massgap>0].shape[0]/massgap.shape[0])
                else:

                    massgap = real[:,6]-real[:,5]
                    get_constraints(massgap, '{}.txt'.format(root), folder=folder)
                    print(np.median(real[:,6]), np.median(real[:,5]), np.median(massgap), massgap[massgap>0].shape[0]/massgap.shape[0])

        elif mode == 'slope_txt':
            if 'likes' not in i and os.path.isfile(os.path.join(base, i))  and 'slope' in i:
                root = i.split('.txt')[0]
                print(root)
                real = np.genfromtxt(base + '{}.txt'.format(root))
                likes = np.genfromtxt(base + '{}_likes.txt'.format(root))

                if '2c' in root:
                    constraints = get_constraints(real[:,8], 'slope/{}.txt'.format(root), folder=folder)
                    print(constraints)
                else:
                    constraints = get_constraints(real[:,5], 'slope/{}.txt'.format(root), folder=folder)
                    print(constraints)

        elif mode == 'mtov_txt':
            if 'likes' not in i and os.path.isfile(os.path.join(base, i))  and 'slope' not in i:
                root = i.split('.txt')[0]
                print(root)
                real = np.genfromtxt(base + '{}.txt'.format(root))
                likes = np.genfromtxt(base + '{}_likes.txt'.format(root))

                if '2c' in root:
                    constraints = get_constraints(real[:,5], 'mTOV/{}.txt'.format(root), folder=folder)
                    print(constraints)
                else:
                    constraints = get_constraints(real[:,2], 'mTOV/{}.txt'.format(root), folder=folder)
                    print(constraints)

        elif mode == 'real_data_99pc':
            if 'likes' not in i and os.path.isfile(os.path.join(base, i)):
                root = i.split('.')[0]
                print(root)
                real = np.genfromtxt(base + '{}.txt'.format(root))
                likes = np.genfromtxt(base + '{}_likes.txt'.format(root))
                indices = np.arange(real.shape[0])
                N = 1000
                cutoff = np.zeros(N)
                minbh = np.zeros(1000)

                if '2c' in root:
                    print('2c')
                elif '1c' in root:
                    for j in range(N):
                        point = int(np.random.choice(indices, 1))
                        point = real[point]
                        minbh[j] = point[3]
                        cutoff[j] = truncnorm.ppf(0.99, loc=point[0], scale=point[1], a=(1-point[0])/point[1], b=(point[2]-point[0])/point[1])
                        if plot:
                            if j % 10 == 0:
                                x = np.linspace(1, 3.5)
                                plt.plot(x, truncnorm.pdf(x, loc=point[0], scale=point[1], a=(1-point[0])/point[1], b=(point[2]-point[0])/point[1]), c='b', alpha=0.1)

                    constraints = get_constraints(minbh-cutoff, '99pc_massgap/{}.txt'.format(root), folder=folder)
                    constraints = get_constraints(cutoff, '99pc/{}.txt'.format(root), folder=folder)
                    print(np.sum(minbh-cutoff > 0)/cutoff.shape[0])

                elif 'u' in root:
                    for j in range(N):
                        point = int(np.random.choice(indices, 1))
                        point = real[point]
                        minbh[j] = point[3]
                        cutoff[j] = (point[2]-1)*0.99 + 1
                        if plot:
                            if j % 10 == 0:
                                x = np.linspace(1,3.5)
                                plt.plot(x, np.ones(50)/(point[2]-1), c='b', alpha=0.1)

                    constraints = get_constraints(minbh-cutoff, '99pc_massgap/{}.txt'.format(root), folder=folder)
                    constraints = get_constraints(cutoff, '99pc/{}.txt'.format(root), folder=folder)
                    print(np.sum(minbh-cutoff > 0)/cutoff.shape[0])

        elif mode == 'real_data_txt':
            if 'likes' not in i and os.path.isfile(os.path.join('../spin_results/real_data/',i)):
                root = i.split('.')[0]
                print(root)
                real = np.genfromtxt(base + '{}.txt'.format(root))
                likes = np.genfromtxt(base + '{}_likes.txt'.format(root))

                if '2c' in root:
                    get_constraints(real[:,5], 'mTOV/{}.txt'.format(root), folder=folder)
                    get_constraints(real[:,6], 'bhmin/{}.txt'.format(root), folder=folder)
                    get_constraints(real[:,7], 'bhslope/{}.txt'.format(root), folder=folder)
                    massgap = real[:,6]-real[:,5]
                    get_constraints(massgap, 'massgap/{}.txt'.format(root), folder=folder)
                    print(np.median(real[:,5]), np.median(real[:,6]),
                          np.median(real[:,7]), massgap[massgap>0].shape[0]/massgap.shape[0])
                else:
                    get_constraints(real[real[:,5]<=1][:,2], 'mTOV/{}.txt'.format(root), folder=folder)
                    get_constraints(real[real[:,5]<=1][:,3], 'bhmin/{}.txt'.format(root), folder=folder)
                    get_constraints(real[real[:,5]<=1][:,4], 'bhslope/{}.txt'.format(root), folder=folder)
                    massgap = real[real[:,5]<=1][:,3]-real[real[:,5]<=1][:,2]
                    get_constraints(massgap, 'massgap/{}.txt'.format(root), folder=folder)
                    print(np.median(real[real[:,5]<=1][:,2]), np.median(real[real[:,5]<=1][:,3]),
                          np.median(real[real[:,5]<=1][:,4]), massgap[massgap>0].shape[0]/massgap.shape[0])
