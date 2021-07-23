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

event_counts = [30, 60, 90, 120, 150]

mtov_True = 2
detector = "APlus"
p.set_detector(detector)
max_jjkep = 1.0
spin_slope = 0.0
bh_min = 5
bh_slope = 2

nsbh_population = p.Population([1.4, 0.5, mtov_True, 1, 3, bh_min, bh_slope], 'nsbh_one', vary_slope=False, selection=True, m1_nospin = True, spinning=True, spin_params=[max_jjkep, spin_slope])
nsbh_population.set_injection_spins(p.injection_set)

pop_samples = nsbh_population.get_population(10, True)
fixed = {"mu": 1.4, "sigma":0.5, "m_TOV":[mtov_True,1.7,3.2], "bh_min": bh_min, "bh_slope": bh_slope, "max_jjkep": max_jjkep, "spin_slope": spin_slope}
for i in range(5):
    samples, likes = nsbh_population.infer(pop_samples, 2000, save_to=None, fixed=fixed, mult=True)
    np.savetxt('./mTOV_convergence/{}_mTOV_{}_run_{}.txt'.format(detector, mtov_True, str(event_counts[i])), samples)
    np.savetxt('./mTOV_convergence/{}_mTOV_{}_run_{}_likes.txt'.format(detector, mtov_True, str(event_counts[i])), likes)

    new_samples = nsbh_population.get_population(10, True)
    pop_samples = np.vstack([pop_samples, new_samples])
    print(pop_samples.shape)

# pop_samples = pop_samples.reshape((180000,4))
# np.savetxt('./mTOV_{}_run_{}_popsamples'.format(mtov_True, str(event_counts[i])), pop_samples)
