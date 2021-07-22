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
p.set_detector("Design")
nsbh_population = p.Population([0.63, 1.35, 0.07, 1.85, 0.35, 2, 1, 3, 5, 2], 'nsbh', vary_slope=False, selection=True, m1_nospin = True, spinning=True, spin_params=[1.0, 0.0])
nsbh_population.set_injection_spins(p.injection_set)

pop_samples = nsbh_population.get_population(10, True)

fixed = {"a":0.63, "mu_1": 1.35, "sigma_1":0.07, "mu_2": 1.85, "sigma_2":0.35, "m_TOV":2, "bh_min": 5, "bh_slope": 2, "max_jjkep":1, "spin_slope":0}

for i in range(5):
    samples, likes = nsbh_population.infer(pop_samples, 2000, mult=True, fixed=fixed)
    np.savetxt('./mTOV_2_run_{}_2component.txt'.format(str(event_counts[i])), samples)
    np.savetxt('./mTOV_2_run_{}_2component.txt'.format(str(event_counts[i])), likes)

    new_samples = nsbh_population.get_population(10, True)
    pop_samples = np.vstack([pop_samples, new_samples])
    print(pop_samples.shape)
