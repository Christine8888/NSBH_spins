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
parser.add_argument("--folder", type=str, default="mTOV_convergence")
parser.add_argument("--free",  action="store_true")


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

nsbh_population = p.Population([1.5, 0.5, mtov_True, 1, 3, bh_min, bh_slope], 'nsbh_one', vary_slope=False, selection=True, m1_nospin = True, spinning=True, spin_params=[max_jjkep, spin_slope])
nsbh_population.set_injection_spins(p.injection_set)

pop_samples = nsbh_population.get_population(event_counts[0], True)

if args.free:
    N = 4000
    fixed = {"m_TOV":[mtov_True,1.7,3.2], "bh_min": bh_min, "bh_slope": bh_slope}
else:
    N = 2000
    fixed = {"mu": 1.5, "sigma":0.5, "m_TOV":[mtov_True,1.7,3.2], "bh_min": bh_min, "bh_slope": bh_slope, "max_jjkep": max_jjkep, "spin_slope": spin_slope}

for i in range(5):
    samples, likes = nsbh_population.infer(pop_samples, N, save_to=None, fixed=fixed, mult=True)
    np.savetxt('../{}/{}_mTOV_{}_run_{}.txt'.format(folder, detector, mtov_True, str(event_counts[i])), samples)
    np.savetxt('../{}/{}_mTOV_{}_run_{}_likes.txt'.format(folder, detector, mtov_True, str(event_counts[i])), likes)
    if i != 4:
        new_samples = nsbh_population.get_population(event_counts[i+1]-event_counts[i], True)
        pop_samples = np.vstack([pop_samples, new_samples])
    print(pop_samples.shape)

# pop_samples = pop_samples.reshape((180000,4))
# np.savetxt('./mTOV_{}_run_{}_popsamples'.format(mtov_True, str(event_counts[i])), pop_samples)
