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
import pandas as pd
import os
from scipy.stats import truncnorm


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

def get_constraints(data, name, folder='real_data', boundaries=None):
    constraints = np.zeros(7)
    sigma_1 = quantiles_from_data(data, [0.159, 0.841])
    sigma_2 = quantiles_from_data(data, [0.025, 0.975])
    sigma_3 = quantiles_from_data(data, [0.0015, 0.997])

    med = quantiles_from_data(data, [0.5])

    constraints = [sigma_3[0], sigma_2[0], sigma_1[0], med[0], sigma_1[1], sigma_2[1], sigma_3[1]]
    if boundaries is not None:
        for i in boundaries:
            constraints.append(data[data>i].shape[0]/data.shape[0])

    #print('../spin_results/outputs/{}/{}'.format(folder, name))
    np.savetxt('../spin_results/outputs/{}/{}'.format(folder, name), constraints)
    return constraints

def plot_constraints(constraints, counts, color='b'):
    # plt.fill_between(counts, constraints[:,0], constraints[:,6], color='b', alpha=0.1, label = '$3\sigma$')
    plt.fill_between(counts, constraints[:,1], constraints[:,5], color=color, alpha=0.2, label = '$2\sigma$')
    plt.fill_between(counts, constraints[:,2], constraints[:,4], color=color, alpha=0.3, label = '$1\sigma$')
    plt.plot(counts, constraints[:,3], c='k')
    plt.scatter(counts, constraints[:,3], c='k', s=30)
    plt.legend()



def processing(folder, mode, dlist=None, plot=True, name='default'):
    base = '../spin_results/{}/'.format(folder)

    if dlist is None:
        dlist = os.listdir(base)

    for i in dlist:


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
                    print(root)
                    massgap = real[:,6]-real[:,5]
                    get_constraints(massgap, '{}.txt'.format(root), folder=folder)
                    print(np.median(real[:,6]), np.median(real[:,5]), np.median(massgap), massgap[massgap>0].shape[0]/massgap.shape[0])


        if mode == 'LMG_txt_99pc':
            if 'likes' not in i and os.path.isfile(os.path.join(base,i)):
                root = i.split('.txt')[0]
                real = np.genfromtxt(base+'{}.txt'.format(root))
                likes = np.genfromtxt(base+'{}_likes.txt'.format(root))


                indices = np.arange(real.shape[0])
                N = 1000
                cutoff = np.zeros(N)
                minbh = np.zeros(1000)

                if '1c' in root:
                    print(root)
                    massgap = real[:,3]-real[:,2]
                    get_constraints(massgap, '{}.txt'.format(root), folder=folder)
                    print(np.median(real[:,2]), np.median(real[:,3]), np.median(massgap), massgap[massgap>0].shape[0]/massgap.shape[0])

                    for j in range(N):
                        point = int(np.random.choice(indices, 1))
                        point = real[point]
                        minbh[j] = point[3]
                        cutoff[j] = truncnorm.ppf(0.99, loc=point[0], scale=point[1], a=(1-point[0])/point[1], b=(point[2]-point[0])/point[1])

                else:
                    print(root)
                    massgap = real[:,6]-real[:,5]
                    get_constraints(massgap, '{}.txt'.format(root), folder=folder)
                    print(np.median(real[:,6]), np.median(real[:,5]), np.median(massgap), massgap[massgap>0].shape[0]/massgap.shape[0])

                    for j in range(N):
                        point = int(np.random.choice(indices, 1))
                        point = real[point]
                        minbh[j] = point[6]
                        cutoff[j] = get_99pc_2c(1000, point[0], point[1], point[2], point[3], point[4], 1, point[5])

                constraints = get_constraints(minbh-cutoff, '99pc_massgap/{}.txt'.format(root), folder=folder, boundaries=[0,1])
                print(constraints)

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
            if 'likes' not in i and os.path.isfile(os.path.join(base, i))  and 'slope' not in i and 'bias' not in i:
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

        elif mode == 'bias_txt':
            if 'likes' not in i and os.path.isfile(os.path.join(base, i))  and 'bias' in i:
                root = i.split('.txt')[0]
                print(root)
                real = np.genfromtxt(base + '{}.txt'.format(root))
                likes = np.genfromtxt(base + '{}_likes.txt'.format(root))

                if '2c' in root:
                    constraints = get_constraints(real[:,5], 'bias/{}.txt'.format(root), folder=folder)
                    print(constraints)
                else:
                    constraints = get_constraints(real[:,2], 'bias/{}.txt'.format(root), folder=folder)
                    print(constraints)

                # 99pc bias

                indices = np.arange(real.shape[0])
                N = 1000
                cutoff = np.zeros(N)

                if '2c' in root:
                    for j in range(N):
                        point = int(np.random.choice(indices, 1))
                        point = real[point]
                        cutoff[j] = get_99pc_2c(1000, point[0], point[1], point[2], point[3], point[4], 1, point[5])

                elif '1c' in root:
                    for j in range(N):
                        point = int(np.random.choice(indices, 1))
                        point = real[point]
                        cutoff[j] = truncnorm.ppf(0.99, loc=point[0], scale=point[1], a=(1-point[0])/point[1], b=(point[2]-point[0])/point[1])

                elif '_u' in root:
                    for j in range(N):
                        point = int(np.random.choice(indices, 1))
                        point = real[point]
                        cutoff[j] = (point[2]-1)*0.99 + 1

                constraints = get_constraints(cutoff, 'bias_99pc/{}.txt'.format(root), folder=folder)
                print(constraints)

        elif mode == 'real_data_99pc':
            if 'likes' not in i and os.path.isfile(os.path.join(base, i)) and name in i:
                root = i.split('.')[0]
                print(root)
                real = np.genfromtxt(base + '{}.txt'.format(root))
                likes = np.genfromtxt(base + '{}_likes.txt'.format(root))
                indices = np.arange(real.shape[0])
                N = 1000
                cutoff = np.zeros(N)
                minbh = np.zeros(1000)

                if '2c' in root:
                    for j in range(N):
                        point = int(np.random.choice(indices, 1))
                        point = real[point]
                        minbh[j] = point[6]
                        cutoff[j] = get_99pc_2c(1000, point[0], point[1], point[2], point[3], point[4], 1, point[5])

                        if plot:
                            if j % 10 == 0:
                                x = np.linspace(1, 3.5)
                                plt.plot(x, p.two_truncnormal_like(x, point[0], point[1], point[2], point[3], point[4], 1, np.ones(50)*point[5]), c='b', alpha=0.1)
                            plt.savefig('../spin_results/outputs/real_data/{}_trace.png'.format(root))
                            plt.close()



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
                            plt.savefig('../spin_results/outputs/real_data/{}_trace.png'.format(root))
                            plt.close()

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
                            plt.savefig('../spin_results/outputs/real_data/{}_trace.png'.format(root))
                            plt.close()

                constraints = get_constraints(minbh-cutoff, '99pc_massgap/{}.txt'.format(root), folder=folder, boundaries=[0,1])
                constraints = get_constraints(cutoff, '99pc/{}.txt'.format(root), folder=folder)
                print(np.sum(minbh-cutoff > 0)/cutoff.shape[0])

        elif mode == 'real_data_txt':
            if 'likes' not in i and os.path.isfile(os.path.join('../spin_results/real_data/',i)) and name in i:
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
                    get_constraints(massgap, 'massgap/{}.txt'.format(root), folder=folder, boundaries=[0,1])
                    print(np.median(real[real[:,5]<=1][:,2]), np.median(real[real[:,5]<=1][:,3]),
                          np.median(real[real[:,5]<=1][:,4]), massgap[massgap>0].shape[0]/massgap.shape[0])

        elif mode == 'real_data_spin':
            dat_190814 = np.genfromtxt('./real_data/gw190814_pos_chieff_samples.txt')
            N = 1000
            spins = np.zeros(N)

            if 'likes' not in i and 'withgw190814' not in i and os.path.isfile(os.path.join('../spin_results/real_data/',i)) and name in i:
                root = i.split('.')[0]
                print(root)
                real = np.genfromtxt(base + '{}.txt'.format(root))
                likes = np.genfromtxt(base + '{}_likes.txt'.format(root))
                indices = np.arange(real.shape[0])

                for j in range(N):
                    m2 = np.random.choice(dat_190814[:,1])

                    if '2c' in i:
                        mTOV = np.random.choice(real[:,5])
                    elif '1c' in i:
                        mTOV = np.random.choice(real[:,2])
                    else:
                        mTOV = np.random.choice(real[:,2])

                    if m2 > p.m_crit(mTOV, 1):
                        spins[j] = 1
                    elif m2 > mTOV:
                        spins[j] = get_spin(m2, mTOV)
                    else:
                        spins[j] = 0

                    get_constraints(spins, 'minspin/{}.txt'.format(root), folder=folder)
                    if plot:
                        plt.hist(spins)
                        plt.savefig('../spin_results/outputs/{}/minspin/{}.png'.format(folder, i))
                        plt.close()


                print(np.quantile(spins, 0.3), np.quantile(spins, 0.5), np.quantile(spins, 0.7))

        elif mode == 'real_data_spin_99pc':
            dat_190814 = np.genfromtxt('./real_data/gw190814_pos_chieff_samples.txt')
            N = 500
            spins = np.zeros(N)

            if 'likes' not in i and 'withgw190814' not in i and os.path.isfile(os.path.join('../spin_results/real_data/',i)) and not os.path.exists('../spin_results/outputs/{}/minspin_99pc/{}.png'.format(folder, i)) and name in i:
                root = i.split('.')[0]
                print(root)
                real = np.genfromtxt(base + '{}.txt'.format(root))
                likes = np.genfromtxt(base + '{}_likes.txt'.format(root))
                indices = np.arange(real.shape[0])

                for j in range(N):
                    m2 = np.random.choice(dat_190814[:,1])

                    if '2c' in i:
                        point = int(np.random.choice(indices))
                        point = real[point]
                        mTOV = get_99pc_2c(1000, point[0], point[1], point[2], point[3], point[4], 1, point[5])
                    elif '1c' in i:
                        point = int(np.random.choice(indices))
                        point = real[point]
                        mTOV = truncnorm.ppf(0.99, loc=point[0], scale=point[1], a=(1-point[0])/point[1], b=(point[2]-point[0])/point[1])
                    else:
                        mTOV = (np.random.choice(real[:,2])-1)*0.99+1

                    if m2 > p.m_crit(mTOV, 1):
                        spins[j] = 1
                    elif m2 > mTOV:
                        spins[j] = get_spin(m2, mTOV)
                    else:
                        spins[j] = 0

                    get_constraints(spins, 'minspin_99pc/{}.txt'.format(root), folder=folder)

                    if plot:
                        plt.hist(spins)
                        plt.savefig('../spin_results/outputs/{}/minspin_99pc/{}.png'.format(folder, i))
                        plt.close()

                print(np.quantile(spins, 0.3), np.quantile(spins, 0.5), np.quantile(spins, 0.7))


def table(folder, mode, already_done = None, dlist=None, save=None, plot=True):
    rows_list = []
    if already_done is None:
        already_done = []
    base_dir = '../spin_results/outputs/{}/'.format(folder)
    if dlist is None:
        dlist = os.listdir(base_dir)

    for i in dlist:
        test = []
        test.append(i.split('_')[0]) # detector
        if test[0] == 'APlus':
            counts = np.linspace(30, 150, num=5, dtype='int16')
        elif test[0] == 'Design':
            counts = np.linspace(10, 50, num=5, dtype='int16')

        if 'png' not in i and os.path.isfile(os.path.join(base_dir, i)):
            # print(i)

            if mode == 'LMG_txt':
                if '_u' in i: # mass model
                    test.append('u')
                elif '_1c' in i:
                    test.append('1c')
                else:
                    test.append('2c')

                test.append(i.split('_')[2]) # mTOV
                test.append(i.split('_')[4]) # bhmin

                if [test[0], test[1], test[2], test[3]] not in already_done:
                    print(i)
                    widths = []
                    constraints = np.zeros((5,7))
                    index = 0
                    for count in counts:
                        mtov_path = base_dir + test[0] + '_mTOV_' + test[2] + '_bhmin_' + test[3] + '_' + str(count) + '.txt'
                        if test[1] is 'u':
                            mtov_path = base_dir + test[0] + '_mTOV_' + test[2] + '_bhmin_' + test[3] + '_' + str(count) + '_u.txt'
                        if test[1] is '1c':
                            mtov_path = base_dir + test[0] + '_mTOV_' + test[2] + '_bhmin_' + test[3] + '_' + str(count) + '_1c.txt'
                        mtov_dat = np.genfromtxt(mtov_path)

                        mtov_str = '${0:.2f}^{{+{1:.2f}}}_{{-{2:.2f}}}$'.format(mtov_dat[3], mtov_dat[4]-mtov_dat[3], mtov_dat[3]-mtov_dat[2])
                        widths.append(mtov_dat[4]-mtov_dat[2])
                        test.append(mtov_str)
                        constraints[index] = mtov_dat[:7]
                        index += 1


                    already_done.append([test[0], test[1], test[2], test[3]])

                    slope, intercept = np.polyfit(np.log10(counts), np.log10(widths), deg=1)
                    test.append(slope)
                    test.append(10**intercept)
                    rows_list.append(test)

                    if plot:
                        plot_constraints(constraints, counts)
                        plt.savefig('{}.png'.format(mtov_path.split('.txt')[0]))
                        print('{}.png'.format(mtov_path.split('.txt')[0]))
                        plt.close()


            elif mode == 'mTOV_txt' and 'bias' not in i:
                if '2c' in i:
                    test.append('2c')
                elif '_u' in i:
                    test.append('u')
                else:
                    test.append('1c')

                test.append(i.split('_')[2]) # mTOV

                if i.split('_run_')[0]+'_'+test[1] not in already_done:
                    print(i)
                    widths = []
                    constraints = np.zeros((5,7))
                    index = 0
                    for count in counts:
                        mtov_path = base_dir + i.split('_run_')[0] + '_run_' + str(count) + '.txt'
                        if test[1] is '2c':
                            mtov_path = base_dir + i.split('_run_')[0] + '_run_' + str(count) + '_2component.txt'
                        if test[1] is 'u':
                            mtov_path = base_dir + i.split('_run_')[0] + '_run_' + str(count) + '_uniform.txt'

                        mtov_dat = np.genfromtxt(mtov_path)

                        mtov_str = '${0:.2f}^{{+{1:.2f}}}_{{-{2:.2f}}}$'.format(mtov_dat[3], mtov_dat[4]-mtov_dat[3], mtov_dat[3]-mtov_dat[2])
                        widths.append(mtov_dat[4]-mtov_dat[2])
                        test.append(mtov_str)
                        constraints[index] = mtov_dat
                        index += 1

                    already_done.append(i.split('_run_')[0]+'_'+test[1])
                    slope, intercept = np.polyfit(np.log10(counts), np.log10(widths), deg=1)

                    if plot:
                        plot_constraints(constraints, counts)
                        plt.savefig('{}{}_{}.png'.format(base_dir, i.split('_run_')[0], test[1]))
                        plt.close()

                    test.append(slope)
                    test.append(10**intercept)
                    rows_list.append(test)

            elif mode == 'bias' and 'bias' in i:
                if '2c' in i:
                    test.append('2c')
                elif '_u' in i:
                    test.append('u')
                elif '1c' in i:
                    test.append('1c')
                test.append(i.split('_')[2])

                if i.split('_run_')[0]+'_'+test[1] not in already_done:
                    print(i)
                    widths = []
                    constraints = np.zeros((5,7))
                    index = 0
                    for count in counts:
                        mtov_path = base_dir.split('/bias')[0] + '/bias/'+ i.split('_run_')[0] + '_run_' + str(count) + '_1c_bias.txt'
                        if test[1] is '2c':
                            mtov_path = base_dir.split('/bias')[0] + '/bias/'+ i.split('_run_')[0] + '_run_' + str(count) + '_2c_bias.txt'
                        if test[1] is 'u':
                            mtov_path = base_dir.split('/bias')[0] + '/bias/'+ i.split('_run_')[0] + '_run_' + str(count) + '_u_bias.txt'

                        mtov_dat = np.genfromtxt(mtov_path)

                        mtov_str = '${0:.2f}^{{+{1:.2f}}}_{{-{2:.2f}}}$'.format(mtov_dat[3], mtov_dat[4]-mtov_dat[3], mtov_dat[3]-mtov_dat[2])
                        test.append(mtov_str)
                        widths.append(mtov_dat[4]-mtov_dat[2])
                        constraints[index] = mtov_dat



                        mtov_path = base_dir.split('/bias')[0] + '/bias_99pc/'+ i.split('_run_')[0] + '_run_' + str(count) + '_1c_bias.txt'
                        if test[1] is '2c':
                            mtov_path = base_dir.split('/bias')[0] + '/bias_99pc/'+ i.split('_run_')[0] + '_run_' + str(count) + '_2c_bias.txt'
                        if test[1] is 'u':
                            mtov_path = base_dir.split('/bias')[0] + '/bias_99pc/'+ i.split('_run_')[0] + '_run_' + str(count) + '_u_bias.txt'

                        mtov_dat = np.genfromtxt(mtov_path)

                        mtov_str = '${0:.2f}^{{+{1:.2f}}}_{{-{2:.2f}}}$'.format(mtov_dat[3], mtov_dat[4]-mtov_dat[3], mtov_dat[3]-mtov_dat[2])
                        test.append(mtov_str)

                        index += 1

                    already_done.append(i.split('_run_')[0]+'_'+test[1])
                    slope, intercept = np.polyfit(np.log10(counts), np.log10(widths), deg=1)

                    if plot:
                        plot_constraints(constraints, counts)
                        plt.savefig('{}{}_{}_{}.png'.format(base_dir, i.split('_run_')[0], 'bias', test[1]))
                        plt.close()
                    #plt.plot(np.log10(counts), np.log10(widths))
                    test.append(slope)
                    test.append(10**intercept)
                    rows_list.append(test)

            elif mode == 'slope':
                if '2c' in i:
                    test.append('2c')
                elif '_u' in i:
                    test.append('u')
                elif '1c' in i:
                    test.append('1c')

                test.append(i.split('_')[2])
                test.append(i.split('_')[-1].split('.txt')[0])
                #print(i, test)
                if 'slope' in i and i.split('_run_')[0]+'_'+test[1]+'_'+test[3] not in already_done and '0.' in i:
                    print(i.split('_run_')[0]+'_'+test[1]+'_'+test[3])
                    widths = []
                    constraints = np.zeros((5,7))
                    index = 0
                    for count in counts:
                        mtov_path = base_dir + i.split('_run_')[0] + '_run_' + str(count) + '_1c_slope_{}.txt'.format(test[3])
                        if test[1] is '2c':
                            mtov_path = base_dir + i.split('_run_')[0] + '_run_' + str(count) + '_2c_slope_{}.txt'.format(test[3])
                        if test[1] is 'u':
                            mtov_path = base_dir + i.split('_run_')[0] + '_run_' + str(count) + '_u_slope_{}.txt'.format(test[3])

                        mtov_dat = np.genfromtxt(mtov_path)

                        mtov_str = '${0:.2f}^{{+{1:.2f}}}_{{-{2:.2f}}}$'.format(mtov_dat[3], mtov_dat[4]-mtov_dat[3], mtov_dat[3]-mtov_dat[2])
                        widths.append(mtov_dat[4]-mtov_dat[2])
                        test.append(mtov_str)
                        constraints[index] = mtov_dat
                        index += 1

                    already_done.append(i.split('_run_')[0]+'_'+test[1]+'_'+test[3])
                    slope, intercept = np.polyfit(np.log10(counts), np.log10(widths), deg=1)
                    if plot:
                        plot_constraints(constraints, counts)
                        plt.savefig('{}{}_{}_{}.png'.format(base_dir, i.split('_run_')[0], 'slope', test[1]))
                        plt.close()
                    #plt.plot(np.log10(counts), np.log10(widths))
                    test.append(slope)
                    test.append(10**intercept)
                    rows_list.append(test)

        if save is not None:
            pd.DataFrame(rows_list).to_csv('../spin_results/{}.csv'.format(save))

def table_realdata(save=True, name='default'):
    rows_list = []
    for i in os.listdir('../spin_results/outputs/real_data/bhmin'):
        if name in i:
            #if '2c' in i and 'nogw190814' in i and 'nogw190426' in i and '0q' in i:
            print(i)
            test = []
            #test.append(index)
            test.append(i.split('_')[2]) #type

            if 'withgw190814' in i:
                test.append(True)
            else:
                test.append(False)
            if 'withgw190426' in i:
                test.append(True)
            else:
                test.append(False)

            test.append(i.split('_')[0]) # mass model

            if '0q' in i:
                test.append(0)
            else:
                test.append(3)

            mtov_path = '../spin_results/outputs/real_data/mTOV/' + i
            mtov_dat = np.genfromtxt(mtov_path)
            mtov_str = '${0:.1f}^{{+{1:.1f}}}_{{-{2:.1f}}}$'.format(mtov_dat[3], mtov_dat[4]-mtov_dat[3], mtov_dat[3]-mtov_dat[2])
            test.append(mtov_str)
            min_str = mtov_dat[1]
            test.append('{0:.1f}'.format(min_str))


            pc_path = '../spin_results/outputs/real_data/99pc/' + i
            pc_dat = np.genfromtxt(pc_path)
            pc_str = '${0:.1f}^{{+{1:.1f}}}_{{-{2:.1f}}}$'.format(pc_dat[3], pc_dat[4]-pc_dat[3], pc_dat[3]-pc_dat[2])
            test.append(pc_str)

            bhmin_path = '../spin_results/outputs/real_data/bhmin/' + i
            bhmin_dat = np.genfromtxt(bhmin_path)
            bhmin_str = '${0:.1f}^{{+{1:.1f}}}_{{-{2:.1f}}}$'.format(bhmin_dat[3], bhmin_dat[4]-bhmin_dat[3], bhmin_dat[3]-bhmin_dat[2])
            test.append(bhmin_str)

            abh_path = '../spin_results/outputs/real_data/bhslope/' + i
            abh_dat = np.genfromtxt(abh_path)
            abh_str = '${0:.1f}^{{+{1:.1f}}}_{{-{2:.1f}}}$'.format(abh_dat[3], abh_dat[4]-abh_dat[3], abh_dat[3]-abh_dat[2])
            test.append(abh_str)

            lmg_path = '../spin_results/outputs/real_data/massgap/' + i
            lmg_dat = np.genfromtxt(lmg_path)
            lmg_str = '${0:.1f}^{{+{1:.1f}}}_{{-{2:.1f}}}$'.format(lmg_dat[3], lmg_dat[4]-lmg_dat[3], lmg_dat[3]-lmg_dat[2])
            test.append(lmg_str)


            lmg_path = '../spin_results/outputs/real_data/99pc_massgap/' + i
            lmg_dat = np.genfromtxt(lmg_path)
            lmg_str = '${0:.1f}^{{+{1:.1f}}}_{{-{2:.1f}}}$'.format(lmg_dat[3], lmg_dat[4]-lmg_dat[3], lmg_dat[3]-lmg_dat[2])
            test.append(lmg_str)
            test.append(lmg_dat[7])
            test.append(lmg_dat[8])

            #if 'withgw190814' in i:
            #    spin_path = '../spin_results/outputs/real_data/minspin/' + i
            #    spin_dat = np.genfromtxt(spin_path)
            #    spin_str = '${0:.1f}^{{+{1:.1f}}}_{{-{2:.1f}}}$'.format(spin_dat[3], spin_dat[4]-spin_dat[3], spin_dat[3]-spin_dat[2])
            #    test.append(spin_str)
            #else:
            #    test.append('')

            #print(test)
            rows_list.append(test)
    if save:
        pd.DataFrame(rows_list).to_csv('../spin_results/real_data.csv')
    return rows_list

def get_99pc_2c(size, a, mu_1, sigma_1, mu_2, sigma_2, lower, upper):
    samples = generate_two_truncnormal(size, a, mu_1, sigma_1, mu_2, sigma_2, lower, upper)
    return np.quantile(samples, 0.99)

def generate_two_truncnormal(size, a, mu_1, sigma_1, mu_2, sigma_2, lower, upper):
    samples = np.zeros(size)
    rands = np.random.rand(size)
    samples[rands<a] = p.generate_truncnormal(samples[rands<a].shape, mu_1, sigma_1, lower, upper)
    samples[rands>=a] = p.generate_truncnormal(samples[rands>=a].shape, mu_2, sigma_2, lower, upper)

    return samples

def get_spin(mTrue, mTOV):
    spins = np.linspace(0,2)
    ms = p.m_crit(mTOV, spins)
    interpfunc = scipy.interpolate.interp1d(ms, spins)
    return interpfunc(mTrue)
