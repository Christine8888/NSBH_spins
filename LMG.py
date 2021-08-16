
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
parser.add_argument("--detector", type=str, default="Design")
parser.add_argument("--event_min", type=int, default=10)
parser.add_argument("--event_max", type=int, default=50)
parser.add_argument("--n_events", type=int, default=5)
parser.add_argument("--mtov_true", type=float, default=2.0)
parser.add_argument("--max_jjkep", type=float, default=1.0)
parser.add_argument("--spin_slope", type=float, default=0.0)
parser.add_argument("--bh_min", type=float, default=5.0)
parser.add_argument("--bh_slope", type=float, default=2.0)
parser.add_argument("--folder", type=str, default="LMG_convergence")

args = parser.parse_args()

event_counts = np.linspace(args.event_min, args.event_max, args.n_events, dtype='int')

mtov_True = args.mtov_true
detector = args.detector
p.set_detector(detector)
max_jjkep = args.max_jjkep
spin_slope = args.spin_slope
bh_min = args.bh_min
bh_slope = args.bh_slope
folder = args.folder

nsbh_population = p.Population([0.63, 1.35, 0.07, 1.85, 0.35, mtov_True, 1, 3, bh_min, bh_slope], 'nsbh', vary_slope=False, selection=True, m1_nospin = True, spinning=True, spin_params=[max_jjkep, spin_slope])
nsbh_population.set_injection_spins(p.injection_set)

pop_samples = nsbh_population.get_population(event_counts[0], True)

fixed = {"a":0.63, "mu_1": 1.35, "sigma_1":0.07, "mu_2": 1.85, "sigma_2":0.35, "m_TOV":[mtov_True,1.7,3.2], "bh_min":[bh_min, bh_min-2, bh_min+2], "bh_slope": bh_slope, "max_jjkep": max_jjkep, "spin_slope": spin_slope}
print(fixed)

for i in range(5):
    samples, likes = nsbh_population.infer(pop_samples, 4000, mult=True, save_to = None,fixed=fixed)
    np.savetxt('../{}/{}_mTOV_{}_bhmin_{}_{}.txt'.format(folder, detector, mtov_True, bh_min, str(event_counts[i])), samples)
    np.savetxt('../{}/{}_mTOV_{}_bhmin_{}_{}_likes.txt'.format(folder, detector, mtov_True, bh_min, str(event_counts[i])), likes)
    if i != 4:
    	new_samples = nsbh_population.get_population(event_counts[i+1]-event_counts[i], True)
    	pop_samples = np.vstack([pop_samples, new_samples])
    print(pop_samples.shape)
