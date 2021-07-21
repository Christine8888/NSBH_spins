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

p.set_detector("Design")
nsbh_population = p.Population([1.4, 0.5, 2, 1, 3, 5, 4], 'nsbh_one', vary_slope=False, selection=True, m1_nospin = True, spinning=True, spin_params=[1.0, 0.0])
nsbh_population.set_injection_spins(p.injection_set)
pop_samples = nsbh_population.get_population(50, True)
print(pop_samples.shape)
pop_samples = pop_samples.reshape((int(50*3000), 4))
print(pop_samples)
np.savetxt('./pop_samples.txt', pop_samples)
#samples, likes = nsbh_population.infer(pop_samples, 20000, save_to = './spinning_test_run.h5')
#np.savetxt('./full_test_run_samples.txt', samples)
#np.savetxt('./full_test_run_likes.txt', likes)

#np.savetxt('./full_test_run_popsamples.txt', pop_samples)
