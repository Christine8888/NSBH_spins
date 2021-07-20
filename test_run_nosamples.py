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

p.set_detector("APlus")
nsbh_population = p.Population([1.4, 0.5, 2, 1, 3, 5, 4], 'nsbh_one', False, selection=True, m1_nospin = True)
nsbh_population.set_injection_spins(p.injection_set)
pop_samples = nsbh_population.get_population(100, False)
samples, likes = nsbh_population.infer(pop_samples, 10000, save_to = './nosample_test_run.h5')
