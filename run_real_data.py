import sys
sys.path.append('~/NSBH_spins')
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import math
import scipy
import emcee
import corner
import lalsimulation as ls
import snr_calculation as s
import astropy.cosmology as cosmo
import astropy.units as u
import populations as p
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--gw190814", action='store_true')
parser.add_argument("--gw190426", action='store_true')
parser.add_argument('--type', default='direct', type=str)
parser.add_argument('--nsamp', default=3000, type=int)
parser.add_argument('--name', type=str)

args = parser.parse_args()

gw190426 =  np.genfromtxt('./real_data/gw190426_{}_samples.txt'.format(args.type))
gw190426 = gw190426[np.random.choice(np.arange(gw190426.shape[0]), args.nsamp)]
gw190426 = np.expand_dims(gw190426, 0)

gw200115 =  np.genfromtxt('./real_data/gw200115_{}_samples.txt'.format(args.type))
gw200115 = gw200115[np.random.choice(np.arange(gw200115.shape[0]), args.nsamp)]
gw200115 = np.expand_dims(gw200115, 0)

gw200105 =  np.genfromtxt('./real_data/gw200105_{}_samples.txt'.format(args.type))
gw200105 = gw200105[np.random.choice(np.arange(gw200105.shape[0]), args.nsamp)]
gw200105 = np.expand_dims(gw200105, 0)

gw190814 =  np.genfromtxt('./real_data/gw190814_{}_samples.txt'.format(args.type))
gw190814 = gw190814[np.random.choice(np.arange(gw190814.shape[0]), args.nsamp)]
gw190814 = np.expand_dims(gw190814, 0)

all_samples = np.append(gw200115, gw200105, axis=0)
gw190814_str = 'nogw190814'
if args.gw190814:
    all_samples = np.append(all_samples, gw190814, axis=0)
    gw190814_str = 'withgw190814'

gw190426_str = 'nogw190426'
if args.gw190814:
    all_samples = np.append(all_samples, gw190426, axis=0)
    gw190426_str = 'withgw190426'


max_jjkep = 1
spin_slope = 0
bh_min = 5
bh_slope = 2
mtov_True = 2

p.set_detector("O3")

nsbh_population = p.Population([1.5, 100, 2, 1, 3, 5, 2], 'nsbh_one', vary_slope=False, selection=True, m1_nospin = True, spinning=True, spin_params=[max_jjkep, spin_slope])
nsbh_population.set_injection_spins(p.injection_set)
nsbh_population.samples = True

fixed = {"mu": 1.5, "sigma":100, "m_TOV":[mtov_True,1.8,4], "bh_min":[bh_min, 1, 10], "bh_slope": [bh_slope, 0, 10], "max_jjkep": [max_jjkep, 0, 1], "spin_slope": [spin_slope, -1, 5]}
samples, likes = nsbh_population.infer(all_samples, 10000, mult=True, save_to = None,fixed=fixed)

np.savetxt('../real_data/{}_{}_{}_{}.txt'.format(args.name, args.type, gw190814_str, gw190426_str), samples)
np.savetxt('../real_data/{}_{}_{}_{}_likes.txt'.format(args.name, args.type, gw190814_str, gw190426_str), likes)
