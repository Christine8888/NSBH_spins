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

event_counts = [10, 20, 30, 40, 50]

mtov_True = 2

p.set_detector("Design")
nsbh_population = p.Population([1.4, 0.5, mtov_True, 1, 3, 5, 2], 'nsbh_one', vary_slope=False, selection=True, m1_nospin = True, spinning=True, spin_params=[1.0, 0.0])
nsbh_population.set_injection_spins(p.injection_set)

pop_samples = nsbh_population.get_population(10, True)
fixed = {"mu":1.4, "sigma":0.5, "m_TOV":[mtov_True,1.7,3.2], "bh_min": 5, "bh_slope": 2, "max_jjkep":1, "spin_slope":0}
for i in range(5):
    samples, likes = nsbh_population.infer(pop_samples, 2000, save_to=None, fixed=fixed, mult=True)
    np.savetxt('./mTOV_convergence/mTOV_{}_run_{}.txt'.format(mtov_True, str(event_counts[i])), samples)
    np.savetxt('./mTOV_convergence/mTOV_{}_run_{}_likes.txt'.format(mtov_True, str(event_counts[i])), likes)

    new_samples = nsbh_population.get_population(10, True)
    pop_samples = np.vstack([pop_samples, new_samples])
    print(pop_samples.shape)

pop_samples = pop_samples.reshape((180000,4))
np.savetxt('./mTOV_{}_run_{}_popsamples'.format(mtov_True, str(event_counts[i])), pop_samples)
