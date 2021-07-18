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
nsbh_population = p.Population([1.4, 0.5, 2, 1, 3, 5, 4], 'nsbh_one', False, selection=True, m1_nospin = True)
pop_samples = nsbh_population.get_population(200, False)
samples, likes = nsbh_population.infer(pop_samples, 20000, save_to = './test_run.h5')
np.savetxt('./test_run_samples.txt', samples)
np.savetxt('./test_run_likes.txt', likes)
np.save_txt('./test_run_popsamples.txt', pop_samples)
