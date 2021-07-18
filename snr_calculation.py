import lalsimulation as ls
import astropy.cosmology as cosmo
import astropy.units as u
import lal
import lalsimulation as ls
import numpy as np

import scipy.interpolate as sip
import h5py


def draw_thetas(N):
    """Draw `N` random angular factors for the SNR.

    Theta is as defined in [Finn & Chernoff
    (1993)](https://ui.adsabs.harvard.edu/#abs/1993PhRvD..47.2198F/abstract).

    Author: Will Farr
    """

    cos_thetas = np.random.uniform(low=-1, high=1, size=N)
    cos_incs = np.random.uniform(low=-1, high=1, size=N)
    phis = np.random.uniform(low=0, high=2*np.pi, size=N)
    zetas = np.random.uniform(low=0, high=2*np.pi, size=N)

    Fps = 0.5*np.cos(2*zetas)*(1 + np.square(cos_thetas))*np.cos(2*phis) - np.sin(2*zetas)*cos_thetas*np.sin(2*phis)
    Fxs = 0.5*np.sin(2*zetas)*(1 + np.square(cos_thetas))*np.cos(2*phis) + np.cos(2*zetas)*cos_thetas*np.sin(2*phis)

    return np.sqrt(0.25*np.square(Fps)*np.square(1 + np.square(cos_incs)) + np.square(Fxs)*np.square(cos_incs))

def next_pow_two(x):
    """Return the next (integer) power of two above `x`.

    Author: Will Farr
    """

    x2 = 1
    while x2 < x:
        x2 = x2 << 1
    return x2

def optimal_snr_detframe(m1_detector, m2_detector, dL, psd_fn=ls.SimNoisePSDaLIGOEarlyHighSensitivityP1200087,psd_file = None, dist_unit = u.Mpc, a1x = 0.0, a1y = 0.0, a1z = 0.0, a2x = 0.0, a2y = 0.0, a2z = 0.0):
    """Return the optimal SNR of a signal.

        :param m1_detector: The detector-frame mass 1 in Msun.

        :param m2_detector: The detector-frame mass 2 in Msun.

        :param dL: The distance in units dist_unit

        :param psd_fn: A function that returns the detector PSD at a given
        frequency (default is early aLIGO high sensitivity, defined in
        [P1200087](https://dcc.ligo.org/LIGO-P1200087/public).

        :param psd_file: A txt file containing frequency in column 1, PSD in column 2. Optional and only used if psd_fn is unused.

        :param dis_unit: The distance unit used

        :a1x, a1y, a1z: cartesian spin coordinates of the primary, by default non-spinning

        :a2x, a2y, a2z: cartesian spin coordinates of the secondary, by default non-spinning

        :return: The SNR of a face-on, overhead source.

        Author: Maya Fishbach and Will Farr
        """

    # Get dL, Gpc
    dL = (dL*dist_unit).to(u.Gpc).value

    # Basic setup: min frequency for w.f., PSD start freq, etc.
    fmin = 19.0
    fref = 40.0
    psdstart = 20.0

    # This is a conservative estimate of the chirp time + MR time (2 seconds)
    tmax = ls.SimInspiralChirpTimeBound(fmin, m1_detector*lal.MSUN_SI, m2_detector*lal.MSUN_SI, a1z, a2z) + 2

    df = 1.0/next_pow_two(tmax)
    fmax = 2048.0 # Hz --- based on max freq of 5-5 inspiral

    # Generate the waveform, redshifted as we would see it in the detector, but with zero angles (i.e. phase = 0, inclination = 0)
    hp, hc = ls.SimInspiralChooseFDWaveform(m1_detector*lal.MSUN_SI, m2_detector*lal.MSUN_SI, a1x, a1y, a1z, a2x, a2y, a2z, dL*1e9*lal.PC_SI, 0.0, 0.0, 0.0, 0.0, 0.0, df, fmin, fmax, fref, None, ls.IMRPhenomPv2)

    Nf = int(round(fmax/df)) + 1
    fs = np.linspace(0, fmax, Nf)
    sel = fs > psdstart

    # PSD
    sffs = lal.CreateREAL8FrequencySeries("psds", 0, 0.0, df, lal.HertzUnit, fs.shape[0])
    if  psd_file is None:
        psd_fn(sffs, psdstart)
        return ls.MeasureSNRFD(hp, sffs, psdstart, -1.0) #http://software.ligo.org/docs/lalsuite/lalsimulation/group___l_a_l_sim_utils__h.html#gabecc316af1ce5186a981fae2a447a194
    else:
        sffs = lal.CreateREAL8FrequencySeries("psds", 0, 0.0, df, lal.HertzUnit, fs.shape[0])
        ls.SimNoisePSDFromFile(sffs,psdstart,psd_file)
        #sffs.data.data = sffs.data.data**0.5 #expects ASD file, so take square root
        return ls.MeasureSNRFD(hp, sffs, psdstart, -1.0)

def snr_sourceframe(m1_source, m2_source, z, psd_fn=ls.SimNoisePSDaLIGOEarlyHighSensitivityP1200087,psd_file = None, a1x = 0.0, a1y = 0.0, a1z = 0.0, a2x = 0.0, a2y = 0.0, a2z = 0.0):

    '''
    Calculate the SNR for a binary with source-frame masses m1_source, m2_source, and spins (a1x, a1y, a1z), (a2x, a2y, a2z)
    with a random sky position and inclination in a single detector with a noise curve described by psd_fn or psd_file
    '''

    #draw a random angular projection factor for the source
    theta = draw_thetas(1)

    #calculate the luminosity distance corresponding to redshift z
    dL = cosmo.Planck15.luminosity_distance(z).to(u.Mpc).value

    #redshift the masses, and calculate the optimal SNR (zero angles)
    opt_snr = optimal_snr_detframe(m1_source * (1+z), m2_source * (1+z), dL, psd_fn,psd_file, u.Mpc, a1x, a1y, a1z, a2x, a2y, a2z)

    #apply the angular projection factor to each source
    snr = opt_snr * theta

    return snr


def optimal_snr_grid(mmin, mmax, dL = 1000.0, dist_unit = u.Mpc, nm = 100, psd_fn=ls.SimNoisePSDaLIGOEarlyHighSensitivityP1200087,psd_file = None, a1x = 0.0, a1y = 0.0, a1z = 0.0, a2x = 0.0, a2y = 0.0, a2z = 0.0):
    ms = np.exp(np.linspace(np.log(mmin), np.log(mmax), nm))
    osnrs = np.zeros((nm, nm))
    for i, m1 in enumerate(ms):
        # print(m1)
        for j in range(i+1):
            #print(i, j)
            m2 = ms[j]
            osnrs[i,j] = optimal_snr_detframe(m1, m2, dL, psd_fn, psd_file, dist_unit, a1x, a1y, a1z, a2x, a2y, a2z)
            osnrs[j,i] = osnrs[i,j]
    return ms, osnrs


def snr_sourceframe_from_grid(m1_source, m2_source, z, ms, osnrs, ref_dL = 1000.0, dist_unit = u.Mpc):
    '''m1_source, m2_source and z are arrays.
    Takes the output of optimal_snr_grid and interpolates
    to speed up snr calculation, assuming all systems have the same spins
    '''

    #construct interpolant for the optimal SNR
    osnr_interp = sip.RectBivariateSpline(ms,ms,osnrs,s=0.0)

    #draw a random angular projection factor for each source
    thetas = draw_thetas(len(m1_source))

    #calculate the luminosity distances corresponding to redshifts z
    dL = cosmo.Planck15.luminosity_distance(z).to(dist_unit).value

    #redshift the masses, and calculate the optimal SNR using the interpolant, using the fact that the SNR scales inversely with the luminosity distance
    opt_snr = osnr_interp.ev(m1_source*(1+z), m2_source*(1+z))*ref_dL/dL

    #apply the angular projection factor to each source
    snr = opt_snr * thetas

    return snr, osnr_interp
