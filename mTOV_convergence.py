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
nsbh_population = p.Population([1.4, 0.5, 2, 1, 3, 5, 2], 'nsbh_one', vary_slope=False, selection=True, m1_nospin = True, spinning=True, spin_params=[1.0, 0.0])
nsbh_population.set_injection_spins(p.injection_set)

pop_samples = nsbh_population.get_population(10, True)

for i in range(5):
    samples, likes = nsbh_population.infer(pop_samples, 2000, mult=True)
    np.savetxt('./mTOV_2_run_{}.txt'.format(str(event_counts[i])), samples)
    np.savetxt('./mTOV_2_run_{}.txt'.format(str(event_counts[i])), likes)

    new_samples = nsbh_population.get_population(10, True)
    pop_samples = np.vstack([pop_samples, new_samples])
    print(pop_samples.shape)
