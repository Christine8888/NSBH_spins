import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import math
import scipy
import emcee
import corner
import snr_calculation as s
import astropy.cosmology as cosmo
import astropy.units as u
from multiprocessing.pool import ThreadPool as Pool
from scipy.integrate import cumtrapz
import scipy.interpolate as sip
from scipy.interpolate import interp1d

DET_THRESHOLD = 8

cosmology_1 = cosmo.wCDM(H0 = 70, Om0 = 0.3, Ode0 = 0.7, w0 = -1.0) #can change this

def Rz_MD(z): #this has arbitrary units
    #return (1+z)**2.7/(1 + ((1+z)/2.9)**5.6)
    return 1

def p_z(zs, cosmology): #here the zs have to be dense enough to permit numerical integration via trapz
    #print(zs)
    dVdz = cosmology.differential_comoving_volume(zs).to(u.Gpc**3/u.sr).value*4*math.pi
    pz_unnorm = Rz_MD(zs) * dVdz * (1+zs)**-1 #not normalized
    pz = pz_unnorm/np.trapz(pz_unnorm, zs) #numerically normalize
    return pz

def set_detector(instrument):
    global injection_set
    global injection_set_bns
    if instrument == "EarlyHigh":
        horizon = 478.4 * u.Mpc
        ms = np.genfromtxt('./injections/nospin_mgrid_EH.txt')
        osnrs = np.genfromtxt('./injections/nospin_snrgrid_EH.txt')

        injection_set = np.genfromtxt('./injections/threshold_8_injections_EH.txt')
        injection_set_bns = np.genfromtxt('./injections/threshold_8_bns_injections_EH.txt')
        max_z = cosmo.z_at_value(cosmo.Planck15.luminosity_distance, horizon, 0, 1)
        # print(type(max_z))
        zs = np.linspace(0.0,max_z, 10000)

    if instrument == "LateHigh":
        horizon = 1196.5 * u.Mpc # horizon for a 2+20 merger
        ms = np.genfromtxt('./injections/nospin_mgrid_LH.txt')
        osnrs = np.genfromtxt('./injections/nospin_snrgrid_LH.txt')
        injection_set = np.genfromtxt('./injections/threshold_8_injections_LH.txt')
        injection_set_bns = np.genfromtxt('./injections/threshold_8_bns_injections_LH.txt')

        max_z = cosmo.z_at_value(cosmo.Planck15.luminosity_distance, horizon, 0, 1)
        zs = np.linspace(0.0,max_z, 10000)

    if instrument == "Design":
        horizon = 1370.9 * u.Mpc # horizon for a 2+20 merger
        ms = np.genfromtxt('./injections/nospin_mgrid_D.txt')
        osnrs = np.genfromtxt('./injections/nospin_snrgrid_D.txt')
        injection_set = np.genfromtxt('./injections/threshold_8_injections_D.txt')
        injection_set_bns = np.genfromtxt('./injections/threshold_8_bns_injections_D.txt')

        max_z = cosmo.z_at_value(cosmo.Planck15.luminosity_distance, horizon, 0, 1)
        zs = np.linspace(0.0,max_z, 10000)

    if instrument == "APlus":
        horizon = 2744.7 * u.Mpc # horizon for a 2+20 merger
        ms = np.genfromtxt('./injections/nospin_mgrid_AP.txt')
        osnrs = np.genfromtxt('./injections/nospin_snrgrid_AP.txt')
        injection_set = np.genfromtxt('./injections/threshold_8_injections_AP.txt')
        injection_set_bns = np.genfromtxt('./injections/threshold_8_bns_injections_AP.txt')

        max_z = cosmo.z_at_value(cosmo.Planck15.luminosity_distance, horizon, 0, 1)
        zs = np.linspace(0.0,max_z, 10000)

    pz_vals = p_z(zs, cosmology_1)
    cumulative_pz = cumtrapz(pz_vals, zs, initial = 0.0)
    global inv_cumulative_pz
    inv_cumulative_pz = sip.interp1d(cumulative_pz, zs)
    snr_local, interp_local = s.snr_sourceframe_from_grid(np.array([10]), np.array([10]), np.array([0.01]), ms, osnrs)
    global snr
    snr = snr_local
    global interp
    interp = interp_local

set_detector("EarlyHigh")

def p_inject_bns(m1, m2, m_min=1, m_max=3, beta=3):
    p_m1 = np.ones(m1.shape[0])/(m_max - m_min)
    p_m2 = like_m2(m2, m1, m_min, beta=3)
    return p_m1 * p_m2

def p_inject_list(m1, m2, m1_min=3, m2_min=1, m2_max=3):
    p_m1 = like_plmin(m1, m1_min, alpha=4)
    p_m2 = np.ones(m2.shape[0])/(m2_max - m2_min)
    return p_m1 * p_m2

def p_inject_one(m1, m2, m1_min=3, m2_min=1, m2_max=3):
    p_m1 = like_plmin_one(m1, m1_min, alpha=4)
    p_m2 = 1/(m2_max - m2_min)
    return p_m1 * p_m2

def generate_distance(N, d_min=1, d_max=1.1): # doesn't work????
    # p_dist goes as D_L^2
    scale = d_max**3 -1
    rand = np.random.rand(N)
    return (rand*scale + 1)**(1/3) + d_min

def generate_z(N):
    rand = np.random.rand(N)
    return inv_cumulative_pz(rand)

def generate_injection(m1_min=3, alpha=4, m2_min=1, m2_max=3):
    m1 = float(generate_plmin(1, m1_min, alpha))
    m2 = np.random.rand()*(m2_max - m2_min) + m2_min
    return np.array([m1, m2])

def generate_injection_bns(m_min=1, m_max=3, beta=3):
    m_1 = np.random.rand()*(m_max-m_min) + m_min
    m_2 = float(generate_q(1, beta, m_1, m_min = m_min))
    return np.array([m_1, m_2])

def check_injection(m1, m2, interp, ref_dL = 1, threshold = DET_THRESHOLD): # ref is 1000 Mpc or 1 Gpc
    # dL = float(generate_distance(1))
    # z = cosmo.z_at_value(cosmo.Planck15.luminosity_distance, dL*u.Gpc)
    z = float(generate_z(1))
    dL = cosmo.Planck15.luminosity_distance(z).to(u.Gpc).value
    theta = float(s.draw_thetas(1))
    # print(interp.ev(m1*(1+z), m2*(1+z)))
    snr = interp.ev(m1*(1+z), m2*(1+z)) * ref_dL/dL
    # print(snr*theta)
    return (snr * theta) > threshold, snr*theta

def create_injection_set(N, generate_injection, interp, threshold = DET_THRESHOLD):
    i = 0
    injection_set = np.zeros((N, 2))
    while i < N:
        injection = generate_injection()
        recovered, value = check_injection(injection[0], injection[1], interp, threshold = threshold)
        if recovered:
            # print(value)
            injection_set[i] = injection
            i += 1
    return injection_set

def generate_normal(size, mu, sigma):
    return np.random.normal(0.0, sigma, size) + mu

def generate_two_normal(size, a, mu_1, sigma_1, mu_2, sigma_2):
    return np.hstack([np.random.normal(0.0, sigma_1, round(size*a)) + mu_1, np.random.normal(0.0, sigma_2, round(size*(1-a))) + mu_2])

def normal_cdf(x, mu, sigma):
    return 0.5*(1+scipy.special.erf((x-mu)/(sigma*np.sqrt(2))))

pi = float(math.pi)

def phi(x, mu, sigma):
    return (1/(sigma*np.sqrt(2*pi))) * np.exp(-(x-mu)**2/(2*sigma**2))

def truncnormal_like(x, mu, sigma, lower, upper):
    #return st.truncnorm.pdf(x, a=-mu/sigma, b=100, loc = mu, scale = sigma)
    result = np.zeros(x.shape[0])
    mask = np.logical_and(x<=upper, x>=lower)
    result[mask] = (phi(x,mu,sigma)/(normal_cdf(upper,mu,sigma)-normal_cdf(lower,mu,sigma)))[mask]
    result[np.logical_or(x>upper, x<lower)] = 0
    return result

def truncnormal_like_one(x, mu, sigma, lower, upper):
    #return st.truncnorm.pdf(x, a=-mu/sigma, b=100, loc = mu, scale = sigma)
    if x <= upper and x >= lower:
        return phi(x,mu,sigma)/(normal_cdf(upper,mu,sigma)-normal_cdf(lower,mu,sigma))
    return 0

def generate_truncnormal(N, mu, sigma, lower, upper):
    return st.truncnorm.rvs(a=(lower-mu)/sigma, b=(upper-mu)/sigma, loc = mu, scale = sigma, size = N)

def generate_two_truncnormal(size, a, mu_1, sigma_1, mu_2, sigma_2, lower, upper):
    samples = np.zeros(size)
    rands = np.random.rand(size)
    for i in range(size):
        if rands[i] < a:
            samples[i] = generate_truncnormal(1, mu_1, sigma_1, lower, upper)
        else:
            samples[i] = generate_truncnormal(1, mu_2, sigma_2, lower, upper)
    return samples

def two_truncnormal_like(x, a, mu_1, sigma_1, mu_2, sigma_2, lower, upper):
    #return st.truncnorm.pdf(x, a=-mu/sigma, b=100, loc = mu, scale = sigma)
    result = np.zeros(x.shape[0])
    mask = np.logical_and(x<=upper, x>=lower)
    # print(result[mask].shape, x[mask].shape)
    result[mask] = a*phi(x[mask],mu_1,sigma_1)/(normal_cdf(upper[mask],mu_1,sigma_1)-normal_cdf(lower,mu_1,sigma_1))+(1-a)*phi(x[mask],mu_2,sigma_2)/(normal_cdf(upper[mask],mu_2,sigma_2)-normal_cdf(lower,mu_2,sigma_2))
    result[~mask] = 0
    return result

def two_truncnormal_like_one(x, a, mu_1, sigma_1, mu_2, sigma_2, lower, upper):
    #return st.truncnorm.pdf(x, a=-mu/sigma, b=100, loc = mu, scale = sigma)
    if x<= upper and x>= lower:
        return a*phi(x,mu_1,sigma_1)/(normal_cdf(upper,mu_1,sigma_1)-normal_cdf(lower,mu_1,sigma_1))+(1-a)*phi(x,mu_2,sigma_2)/(normal_cdf(upper,mu_2,sigma_2)-normal_cdf(lower,mu_2,sigma_2))
    return 0

def m_crit(m_TOV, j_jkep):
    a2 = 1.316e-1
    a4 = 7.111e-2
    return (1 + a2*(j_jkep)**2 + a4*(j_jkep)**4)*m_TOV

def m_crit_slope(m_TOV, slope, j_jkep):
    return m_TOV + slope*j_jkep

def like_m2(m_2, m_1, m_min=1, beta=3):
    like = m_2**beta * (beta+1)
    like /= (m_1**(beta+1)-m_min**(beta+1))
    return like

def generate_q(N, beta, m_1, m_min = 1):
    draws = np.random.rand(N, 2)
    draws = draws * ((m_1**(beta+1))-(m_min**(beta+1))) # first scale to within (0,1)
    draws += m_min**(beta+1) # then shift so m_min is 0
    draws = draws**(1/(beta+1)) # then complete inverse transform
    return draws[:, 0]

def get_rho(N, index=4, rhomin = 10):
    r = np.random.rand(N)
    return rhomin * (1-r) ** (-1/(index-1))

def chirp_mass(m_1, m_2):
    return (m_1*m_2)**(3/5)/((m_1+m_2)**(1/5))

def chirp_mass_sigma(m_1, m_2, rho=10):
    alpha = 0
    v = mass_ratio(m_1, m_2)

    if v > 0.1:
        alpha = 0.01
    elif v < 0.1 and v > 0.05:
        alpha = 0.03
    elif v < 0.05:
        alpha = 0.1

    return chirp_mass(m_1, m_2)*alpha*12/rho

def mass_ratio(m_1, m_2):
    q=m_2/m_1
    v=q/((1+q)**2)
    return v

def mass_ratio_sigma(m_1, m_2,rho=10):
    return mass_ratio(m_1, m_2)*0.03*12/rho

# from https://arxiv.org/pdf/1711.09226.pdf LALInference simulations. does it need to be multipled by actual chieff?

def chieff_sigma(rho):
    return 0.0086 * ((rho/15)**(-0.41))

# assuming perfectly aligned
def chieff(m_1, m_2, chi_1, chi_2):
    val = m_1 * chi_1 + m_2 * chi_2
    val /= (m_1 + m_2)
    return val

def jacobian(m_1, m_2):
    pa = m_1**(13/5) * m_2**(3/5)
    pa /= (m_1 + m_2)**(21/5)

    pb = m_1**(8/5) * m_2**(8/5)
    pb /= (m_1 + m_2)**(21/5)

    return pa - pb

def pl_cdf(x, x_min, alpha):
    return 1-(x/x_min)**(-alpha+1)

def like_plminmax(x, m_min, m_max, alpha):
    result = np.zeros(x.shape[0])
    mask = np.logical_and(x<=m_max, x>=m_min)
    result[mask] = ((alpha-1)/m_min)*(x[mask]/m_min)**(-alpha)/pl_cdf(m_max, m_min, alpha)
    result[~mask] = 0
    return result

def like_plmin(x, m_min, alpha):
    result = np.zeros(x.shape[0])
    mask = x>=m_min
    result[mask] = ((alpha-1)/m_min)*(x[mask]/m_min)**(-alpha)
    result[~mask] = 0
    return result

def like_plminmax_one(x, m_min, m_max, alpha):
    if x < m_max and x > m_min:
        return ((alpha-1)/m_min)*(x/m_min)**(-alpha)/pl_cdf(m_max, m_min, alpha)
    return 0

def like_plmin_one(x, m_min, alpha):
    if x > m_min:
        return ((alpha-1)/m_min)*(x/m_min)**(-alpha)
    return 0

def like_beta(x, beta):
    # beta = beta+1
    result = np.zeros(x.shape[0])
    mask = np.logical_and(x<=1, x>0)
    result[mask] = x[mask]**beta/(beta+1)
    return result

def like_beta_one(x, beta):
    if x <=1 and x >= 0:
        return x**beta/(beta+1)
    return 0

def generate_plmin(N, x_min, alpha):
    rand = np.random.rand(N)
    pl_val = x_min * (1 - rand) ** (-1 / (alpha - 1))
    return pl_val

def generate_plminmax(N, x_min, x_max, alpha):
    rand = np.random.rand(N)
    rand *= pl_cdf(x_max, x_min, alpha)
    pl_val = x_min * (1 - rand) ** (-1 / (alpha - 1))
    return pl_val

def draw_NSBH(params, spin, vary_slope = False, a_bh = 0.5, pop_type = 'nsbh'):
    if pop_type == 'nsbh':
        a_ns2 = params[0]
        mu_1 = params[1]
        sigma_1 = params[2]
        mu_2 = params[3]
        sigma_2 = params[4]
        m_TOV = params[5]
        max_jjkep = params[6]
        beta = params[7]
        bh_min = params[8]
        bh_slope = params[9]
    else:
        mu = params[0]
        sigma = params[1]
        m_TOV = params[2]
        max_jjkep = params[3]
        beta = params[4]
        bh_min = params[5]
        bh_slope = params[6]

    if vary_slope:
        if a_bh < np.random.rand():
            # return float(generate_plminmax(1, bh_min, bh_max, bh_slope))
            return float(generate_plmin(1, bh_min, bh_slope))
        else:
            if pop_type == 'nsbh':
                return float(generate_two_truncnormal(1, a_ns2, mu_1, sigma_1, mu_2, sigma_2, lower=1, upper=m_crit_slope(m_TOV, params[10], spin)))
            elif pop_type == 'nsbh_one':
                return float(generate_truncnormal(1, mu, sigma, lower=1, upper=m_crit_slope(m_TOV, params[7], spin)))
    else:
        if a_bh < np.random.rand():
            # return float(generate_plminmax(1, bh_min, bh_max, bh_slope))
            return float(generate_plmin(1, bh_min, bh_slope))
        else:
            if pop_type == 'nsbh':

                return float(generate_two_truncnormal(1, a_ns2, mu_1, sigma_1, mu_2, sigma_2, lower=1, upper=m_crit(m_TOV, spin)))
            elif pop_type == 'nsbh_one':
                return float(generate_truncnormal(1, mu, sigma, lower=1, upper=m_crit(m_TOV, spin)))

def generate_NSBH(N, params, nsbh_only = False, vary_slope = False, spinning=False, a_bh = 0.5, pop_type = 'nsbh', spin_params = [1.0, 0.0]):
    total = 0
    pop = np.zeros((N, 4))

    if pop_type == 'nsbh':
        beta = params[7]
        max_jjkep = params[6]
    elif pop_type == 'nsbh_one':
        beta = params[4]
        max_jjkep = params[3]

    while total < N:
        if spinning:
            s1 = float(1-generate_q(1, spin_params[1], 1, 1-spin_params[0]))
            s2 = float(1-generate_q(1, spin_params[1], 1, 1-spin_params[0]))
        else:
            s1 = np.random.rand()*max_jjkep
            s2 = np.random.rand()*max_jjkep

        if nsbh_only:
            d1 = draw_NSBH(params, s1, vary_slope, a_bh = 1, pop_type = pop_type)
            d2 = draw_NSBH(params, s2, vary_slope, a_bh = 0, pop_type = pop_type)
        else:
            d1 = draw_NSBH(params, s1, vary_slope, a_bh = a_bh, pop_type = pop_type)
            d2 = draw_NSBH(params, s2, vary_slope, a_bh = a_bh, pop_type = pop_type)

        if d1 > d2:
            m_1 = d1
            m_2 = d2
            spin_1 = s1 # set to 0 for non-spinning black hole case
            spin_2 = s2
        elif d2 > d1:
            m_1 = d2
            m_2 = d1
            spin_1 = s2
            spin_2 = s1

        q = m_2/m_1

        if np.random.rand() < q**beta:
            pop[total, :] = np.array([m_1, m_2, spin_1, spin_2])
            total += 1

    return pop

class Population():
    def __init__(self, params, pop_type, vary_slope=False, selection=False, spinning = False, m1_nospin = False, spin_params = [1.0, 0.0]):
        """
        Population of compact objects.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal, "nsbh" for NSBH
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?
        spinning -- free spin parameters?
        spin_params -- maximum j/j_kep, slope (1-chi^b)
        """

        if pop_type == "one":
            self.mu_1 = params[0]
            self.sigma_1 = params[1]
            self.m_TOV = params[2]
            self.max_jjkep = params[3]
            self.beta = params[4]
            if vary_slope:
                self.slope = params[5]

        elif pop_type == "two":
            self.a = params[0]
            self.mu_1 = params[1]
            self.sigma_1 = params[2]
            self.mu_2 = params[3]
            # print(self.mu_2)
            self.sigma_2 = params[4]
            self.m_TOV = params[5]
            self.max_jjkep = params[6]
            self.beta = params[7]
            if vary_slope:
                self.slope = params[8]

        elif pop_type == "nsbh_one":
            self.mu = params[0]
            self.sigma = params[1]
            self.m_TOV = params[2]
            self.max_jjkep = params[3]
            self.beta = params[4]
            self.bh_min = params[5]
            self.bh_slope = params[6]
            self.a_bh = 0.5
            if vary_slope:
                self.slope = params[7]

        elif pop_type == "nsbh":
            self.a = params[0]
            self.mu_1 = params[1]
            self.sigma_1 = params[2]
            self.mu_2 = params[3]
            self.sigma_2 = params[4]
            self.m_TOV = params[5]
            self.max_jjkep = params[6]
            self.beta = params[7]
            self.bh_min = params[8]
            self.bh_slope = params[9]
            self.a_bh = 0.5
            if vary_slope:
                self.slope = params[10]

        self.pop_type = pop_type
        self.vary_slope = vary_slope
        print('vary slope: {}'.format(self.vary_slope))
        self.selection = selection
        print('selection: {}'.format(self.selection))
        self.spinning = spinning
        print('spinning: {}'.format(self.spinning))
        self.m1_nospin = m1_nospin
        print('no m1 spin: {}'.format(self.m1_nospin))
        if self.spinning:
            self.max_jjkep = spin_params[0]
            self.spin_slope = spin_params[1]

    def uniform_spin_one(self, spin, max_jjkep):
        if spin < max_jjkep:
            return 1/max_jjkep
        return 0

    def uniform_spin(self, spin, max_jjkep):
        result = np.zeros(spin.shape[0])
        result[spin < max_jjkep] = 1/max_jjkep
        result[spin >= max_jjkep] = 0
        return result

    def pl_spin_one(self, spin, max_jjkep, spin_slope):
        if spin < max_jjkep:
            return like_m2(1-spin, 1, 1-max_jjkep, spin_slope)
        return 0

    def pl_spin(self, spin, max_jjkep, spin_slope):
        result = np.zeros(spin.shape[0])
        mask = spin < max_jjkep
        # result[mask] = p.like_m2(1, m_1 = 1, m_min = 1-max_jjkep, beta = spin_slope)-p.like_m2(spin[mask], 1, 1-max_jjkep, spin_slope)
        result[mask] = like_m2(1-spin[mask], 1, 1-max_jjkep, spin_slope)
        result[~mask] = 0
        return result

    def selection_norm(self, params, injection_set):
        N = injection_set.shape[0]
        new_set = np.hstack([injection_set, np.zeros((N,2))])
        if self.pop_type == 'one':
            return (1/N) * np.sum(self.event_likelihood_one_one_samples(new_set, params, nomean=True)/p_inject_bns(new_set[:,0], new_set[:,1]))
        elif self.pop_type == 'two':
            # print([self.event_likelihood_two_single(np.array([i]), params) for i in new_set])
            return (1/N) * np.sum(self.event_likelihood_two_samples(new_set, params, nomean=True)/p_inject_bns(new_set[:,0], new_set[:,1]))
        elif self.pop_type == 'nsbh':
            return (1/N) * np.sum(self.event_likelihood_nsbh_samples(new_set, params, nomean=True)/p_inject_list(new_set[:,0], new_set[:,1]))
        elif self.pop_type == 'nsbh_one':
            return (1/N) * np.sum(self.event_likelihood_nsbh_one_samples(new_set, params, nomean=True)/p_inject_list(new_set[:,0], new_set[:,1]))


    def get_population(self, N, samples=True, N_samples=500):
        self.N = N
        self.samples = samples
        self.N_samples = N_samples

        n = self.N_samples

        if self.samples:
            population = np.zeros((N, int((n*0.75)*8), 4))
            i = 0
            while i < N:
                new_draw = self.get_samples(n, return_p0=False)
                if new_draw is not None:
                    population[i] = new_draw
                    i += 1
            return population
        else:
            population = np.zeros((N, 1, 4))
            i = 0
            while i < N:
                new_draw = self.get_samples(n, return_p0=True)
                if new_draw is not None:
                    population[i] = new_draw
                    i += 1
            return population



    def get_samples(self, n, return_p0=False):
        if self.m1_nospin:
            test_chi_1 = 0
        else:
            if self.spinning:
                test_chi_1 = float(1-generate_q(1, self.spin_slope, 1, 1-self.max_jjkep))
            else:
                test_chi_1 = np.random.rand() * self.max_jjkep

        if self.spinning:
            # print('here')
            test_chi_2 = float(1-generate_q(1, self.spin_slope, 1, 1-self.max_jjkep))
        else:
            test_chi_2 = np.random.rand() * self.max_jjkep


        if self.pop_type == "one":
            if self.vary_slope:
                test_m_1 = float(generate_truncnormal(1, self.mu_1, self.sigma_1, lower=1, \
                                                      upper=m_crit_slope(self.m_TOV, self.slope, test_chi_1)))

                test_m_2 = float(generate_q(1, self.beta, np.min([test_m_1, m_crit_slope(self.m_TOV, self.slope, test_chi_2)]), m_min=1))
            else:
                test_m_1 = float(generate_truncnormal(1, self.mu_1, self.sigma_1, lower=1, \
                                                      upper=m_crit(self.m_TOV, test_chi_1)))

                test_m_2 = float(generate_q(1, self.beta, np.min([test_m_1, m_crit(self.m_TOV, test_chi_2)]), m_min=1))


        elif self.pop_type == "two":
            if self.vary_slope:
                test_m_1 = float(generate_two_truncnormal(1, self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, \
                                                          lower=1, upper=m_crit_slope(self.m_TOV, self.slope, test_chi_1)))
                test_m_2 = float(generate_q(1, self.beta, np.min([test_m_1, m_crit_slope(self.m_TOV, self.slope, test_chi_2)]), m_min=1))
            else:
                test_m_1 = float(generate_two_truncnormal(1, self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, \
                                                          lower=1, upper=m_crit(self.m_TOV, test_chi_1)))
                test_m_2 = float(generate_q(1, self.beta, np.min([test_m_1, m_crit(self.m_TOV, test_chi_2)]), m_min=1))

        elif self.pop_type == "nsbh":
            if self.vary_slope:
                params = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.max_jjkep, self.beta, \
                      self.bh_min, self.bh_slope, self.slope]
            else:
                params = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.max_jjkep, self.beta, \
                      self.bh_min, self.bh_slope]

            if self.spinning:
                pop = generate_NSBH(1, params, nsbh_only = True, vary_slope = self.vary_slope, spinning = self.spinning, spin_params = [self.max_jjkep, self.spin_slope])
            else:
                pop = generate_NSBH(1, params, nsbh_only = True, vary_slope = self.vary_slope)
            # print(pop)
            test_m_1 = pop[0,0]
            test_m_2 = pop[0,1]
            test_chi_1 = pop[0,2]
            if self.m1_nospin:
                test_chi_1 = 0
            test_chi_2 = pop[0,3]

        elif self.pop_type == "nsbh_one":
            if self.vary_slope:
                params = [self.mu, self.sigma, self.m_TOV, self.max_jjkep, self.beta, \
                      self.bh_min, self.bh_slope, self.slope]
            else:
                params = [self.mu, self.sigma, self.m_TOV, self.max_jjkep, self.beta, \
                      self.bh_min, self.bh_slope]

            if self.spinning:
                pop = generate_NSBH(1, params, nsbh_only = True, vary_slope = self.vary_slope, pop_type = 'nsbh_one', spinning = self.spinning, spin_params = [self.max_jjkep, self.spin_slope])
            else:
                pop = generate_NSBH(1, params, nsbh_only = True, vary_slope = self.vary_slope, pop_type = 'nsbh_one',)

            # print(pop)
            test_m_1 = pop[0,0]
            test_m_2 = pop[0,1]
            test_chi_1 = pop[0,2]
            if self.m1_nospin:
                test_chi_1 = 0
            test_chi_2 = pop[0,3]

        if self.selection:
            detected, event_snr = check_injection(test_m_1, test_m_2, interp)
            # print(event_snr)
            if not detected:
                return None
            test_rho = event_snr

        else:
            test_rho = get_rho(1, 4, 12) # should replace w/ a distance
        # print(test_rho)
        test_mchirp = chirp_mass(test_m_1, test_m_2)
        test_mchirp_sigma = float(chirp_mass_sigma(test_m_1, test_m_2, test_rho))
        test_mchirp += float(np.random.randn()*test_mchirp_sigma)

        test_ratio = mass_ratio(test_m_1, test_m_2)
        test_ratio_sigma = float(mass_ratio_sigma(test_m_1, test_m_2,test_rho))
        test_ratio += float(np.random.randn()*test_ratio_sigma)

        test_chieff = chieff(test_m_1, test_m_2, test_chi_1, test_chi_2)
        test_chieff_sigma = float(chieff_sigma(test_rho))
        test_chieff += float(np.random.randn()*test_chieff_sigma)

        p0 = [test_m_1, test_m_2, test_chi_1, test_chi_2]
        # print(p0)
        if return_p0:
            return p0

        if self.m1_nospin:
            p0 = [test_m_1, test_m_2, test_chi_1, test_chi_2]
            pscale = [0.1, 0.1, 0.01, 0.05]
        else:
            pscale = [0.1, 0.1, 0.05, 0.05]

        pscale /= (test_rho/8)
        pos = p0 + pscale*np.random.randn(8, 4)
        pos = np.abs(pos)
        nwalkers, ndim = pos.shape

        def loglike_one(params, m_chirp, m_chirp_sigma, m_ratio, m_ratio_sigma, chi_eff, chi_eff_sigma):
            #print(params)
            # params: m_1, m_2, chi_1, chi_2
            like_mchirp = st.norm.pdf(chirp_mass(params[0], params[1]), loc=m_chirp, scale=m_chirp_sigma)
            like_ratio = st.norm.pdf(mass_ratio(params[0], params[1]), loc=m_ratio, scale=m_ratio_sigma)
            like_chieff = st.norm.pdf(chieff(params[0], params[1], params[2], params[3]), loc=chi_eff, scale=chi_eff_sigma)
            return np.log(like_mchirp) + np.log(like_ratio) + np.log(like_chieff) + np.log(jacobian(params[0], params[1]))


        def logpost_one(params, m_chirp, m_chirp_sigma, mass_ratio, mass_ratio_sigma, chi_eff, chi_eff_sigma):
            if self.pop_type == 'nsbh' or self.pop_type == 'nsbh_one':
                m1range = [3, 25]
                m2range = [0.5, 3.5]
            else:
                m1range = [0.5, 3.5]
                m2range = [0.5, 3.5]
            if params[0] > m1range[0] and params[0] < m1range[1]:
                #
                if params[1] > m2range[0] and params[1] < m2range[1]:

                    if params[0] > params[1]:

                        if self.m1_nospin:

                            if params[2] >= 0 and params[2] <= 0.01:
                                # print('here')
                                if params[3] >= 0 and params[3] <= 1:
                                    loglike = loglike_one(params, m_chirp, m_chirp_sigma, mass_ratio, mass_ratio_sigma, chi_eff, chi_eff_sigma)
                                    return loglike
                        else:
                            if params[2] >= 0 and params[2] <= 1:
                                if params[3] >= 0 and params[3] <= 1:
                                    loglike = loglike_one(params, m_chirp, m_chirp_sigma, mass_ratio, mass_ratio_sigma, chi_eff, chi_eff_sigma)
                                    return loglike
            return -np.inf

        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost_one, args=(test_mchirp, test_mchirp_sigma, test_ratio, test_ratio_sigma, test_chieff, test_chieff_sigma))
        sampler.run_mcmc(pos, n, progress=False)
        #t_auto = np.mean(sampler.get_autocorr_time())
        samples = sampler.get_chain(discard = int(n/4), flat=True, thin=1)
        log_prob_samples = sampler.get_log_prob(discard = 100, flat=True)
        # print(log_prob_samples)
        return samples #, log_prob_samples

    def event_likelihood_one_samples(self, samples, params, nomean=False):
        # params: mu, sigma, m_TOV, slope (optional), spin parameters
        if self.vary_slope:
            truncs = truncnormal_like(samples[:,0], params[0], params[1], lower = 1, upper = m_crit_slope(params[2], params[3], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:,0], m_crit_slope(params[2], params[3], samples[:,3])], axis=0))
            if self.spinning:
                # spin_likes = self.uniform_spin(samples[:,2], params[4]) * self.uniform_spin(samples[:,3], params[4])
                spin_likes = self.pl_spin(samples[:,3], params[4], params[5])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin(samples[:,2], params[4], params[5])
            else:
                spin_likes = 1
        else:
            truncs = truncnormal_like(samples[:,0], params[0], params[1], lower = 1, upper = m_crit(params[2], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:,0], m_crit(params[2], samples[:,3])], axis=0))
            if self.spinning:
                # spin_likes = self.uniform_spin(samples[:,2], params[3]) * self.uniform_spin(samples[:,3], params[3])
                spin_likes = self.pl_spin(samples[:,3], params[3], params[4])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin(samples[:,2], params[3], params[4])
            else:
                spin_likes = 1
        if nomean:
            return truncs*qlikes*spin_likes
        return np.mean(truncs*qlikes*spin_likes)

    def event_likelihood_two_samples(self, samples, params, nomean=False):
        # params: 0 a, 1 mu_1, 2 sigma_1, 3 mu_2, 4 sigma_2, 5 m_TOV, 6 slope (optional), spin parameters
        if self.vary_slope:
            truncs = two_truncnormal_like(samples[:,0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[6], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:, 0], m_crit_slope(params[5], params[6], samples[:, 3])], axis=0))
            if self.spinning:
                # spin_likes = self.uniform_spin(samples[:,2], params[7]) * self.uniform_spin(samples[:,3], params[7])
                spin_likes = self.pl_spin(samples[:,3], params[7], params[8])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin(samples[:,2], params[7], params[8])
            else:
                spin_likes = 1
        else:
            truncs = two_truncnormal_like(samples[:,0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:, 0], m_crit(params[5], samples[:, 3])], axis=0))
            if self.spinning:
                # spin_likes = self.uniform_spin(samples[:,2], params[6]) * self.uniform_spin(samples[:,3], params[6])
                spin_likes = self.pl_spin(samples[:,3], params[6], params[7])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin(samples[:,2], params[6], params[7])
            else:
                spin_likes = 1
        if nomean:
            return truncs*qlikes*spin_likes
        return np.mean(truncs*qlikes*spin_likes)

    def event_likelihood_nsbh_one_samples(self, samples, params, nomean=False):
        # params: mu, sigma, m_TOV, bh_min, bh_slope, slope (optional)
        # REPLACE LIKELIHOODS!

        p_m1 = like_plmin(samples[:,0], params[3], params[4])
        q = samples[:,1]/samples[:,0]
        p_q = like_beta(q, self.beta)

        if self.vary_slope:

            p_m2 = truncnormal_like(samples[:,1], params[0], params[1], lower = 1, upper = m_crit_slope(params[2], params[5], samples[:,3]))

            if self.spinning:
                # spin_likes = self.uniform_spin(samples[:,2], params[9]) * self.uniform_spin(samples[:,3], params[9])
                spin_likes = self.pl_spin(samples[:,3], params[6], params[7])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin(samples[:,2], params[6], params[7])
            else:
                spin_likes = 1
        else:
            p_m2 = truncnormal_like(samples[:,1], params[0], params[1], lower = 1, upper = m_crit(params[2], samples[:,3]))

            if self.spinning:
                # spin_likes = self.uniform_spin(samples[:,2], params[8]) * self.uniform_spin(samples[:,3], params[8])
                spin_likes = self.pl_spin(samples[:,3], params[5], params[6])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin(samples[:,2], params[5], params[6])
            else:
                spin_likes = 1
        if nomean:
            return p_m1*p_m2*p_q*spin_likes
        return np.mean(p_m1*p_m2*p_q*spin_likes)

    def event_likelihood_nsbh_samples(self, samples, params, nomean=False):
        # params: a, mu_1, sigma_1, mu_2, sigma_2, m_TOV, bh_min, bh_slope, slope (optional)
        # REPLACE LIKELIHOODS!
        if self.vary_slope:
            p_m1 = like_plmin(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[8], samples[:,3]))

            q = samples[:,1]/samples[:,0]
            p_q = like_beta(q, self.beta)

            if self.spinning:
                # spin_likes = self.uniform_spin(samples[:,2], params[9]) * self.uniform_spin(samples[:,3], params[9])
                spin_likes = self.pl_spin(samples[:,3], params[9], params[10])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin(samples[:,2], params[9], params[10])
            else:
                spin_likes = 1
        else:
            p_m1 = like_plmin(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], samples[:,3]))
            q = samples[:,1]/samples[:,0]
            p_q = like_beta(q, self.beta)

            if self.spinning:
                # spin_likes = self.uniform_spin(samples[:,2], params[8]) * self.uniform_spin(samples[:,3], params[8])
                spin_likes = self.pl_spin(samples[:,3], params[8], params[9])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin(samples[:,2], params[8], params[9])
            else:
                spin_likes = 1
        if nomean:
            return p_m1*p_m2*p_q*spin_likes
        return np.mean(p_m1*p_m2*p_q*spin_likes)

    def event_likelihood_one_single(self, samples, params):
        i = samples[0]

        if self.vary_slope:
            truncs = truncnormal_like_one(i[0], params[0], params[1], lower = 1, upper = m_crit_slope(params[2], params[3], i[2]))
            qlikes = like_m2(i[1], min([i[0], m_crit_slope(params[2], params[3], i[3])]))
            if self.spinning:
                spin_likes = self.pl_spin_one(i[3], params[4], params[5])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin_one(i[2], params[4], params[5])
            else:
                spin_likes = 1
        else:
            truncs = truncnormal_like_one(i[0], params[0], params[1], lower = 1, upper = m_crit(params[2], i[2]))
            qlikes = like_m2(i[1], min([i[0], m_crit(params[2], i[3])]))
            if self.spinning:
                spin_likes = self.pl_spin_one(i[3], params[3], params[4])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin_one(i[2], params[3], params[4])
            else:
                spin_likes = 1
        return np.mean(truncs*qlikes*spin_likes)

    def event_likelihood_two_single(self, samples, params):
        i = samples[0]
        if self.vary_slope:
            truncs = two_truncnormal_like_one(x = i[0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                            mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[6], i[2]))
            qlikes = like_m2(i[1], min([i[0], m_crit_slope(params[5], params[6], i[3])]))
            if self.spinning:
                spin_likes = self.pl_spin_one(i[3], params[7], params[8])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin_one(i[2], params[7], params[8])
            else:
                spin_likes = 1
        else:
            truncs = two_truncnormal_like_one(x = i[0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                            mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], i[2]))
            qlikes = like_m2(i[1], min([i[0], m_crit(params[5], i[3])]))
            if self.spinning:
                spin_likes = self.pl_spin_one(i[3], params[6], params[7])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin_one(i[2], params[6], params[7])
            else:
                spin_likes = 1
        return np.mean(truncs*qlikes*spin_likes)

    def event_likelihood_nsbh_single(self, samples, params):
        i = samples[0]
        # params: a, mu_1, sigma_1, mu_2, sigma_2, m_TOV, bh_min, bh_slope, slope (optional)
        # REPLACE LIKELIHOODS!
        if self.vary_slope:
            p_m1 = like_plmin_one(i[0], params[6], params[7])
            p_m2 = two_truncnormal_like_one(i[1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[8], i[3]))

            q = i[1]/i[0]
            p_q = like_beta_one(q, self.beta)
            if self.spinning:
                spin_likes = self.pl_spin_one(i[3], params[9], params[10])

                if not self.m1_nospin:
                    spin_likes *= self.pl_spin_one(i[2], params[9], params[10])
            else:
                spin_likes = 1

        else:
            p_m1 = like_plmin_one(i[0], params[6], params[7])
            # print(p_m1)
            p_m2 = two_truncnormal_like_one(i[1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], i[3]))
            q = i[1]/i[0]
            p_q = like_beta_one(q, self.beta)
            if self.spinning:
                spin_likes = self.pl_spin_one(i[3], params[8], params[9])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin_one(i[2], params[8], params[9])
            else:
                spin_likes = 1
            # print(spin_likes)
        return np.mean(p_m1*p_m2*p_q*spin_likes)

    def event_likelihood_nsbh_one_single(self, samples, params):
        i = samples[0]
        p_m1 = like_plmin_one(i[0], params[3], params[4])
        # #print(i)
        q = i[1]/i[0]
        p_q = like_beta_one(q, self.beta)

        if self.vary_slope:

            p_m2 = truncnormal_like_one(i[1], params[0], params[1], lower = 1, upper = m_crit_slope(params[2], params[5], i[3]))

            if self.spinning:
                spin_likes = self.pl_spin_one(i[3], params[6], params[7])

                if not self.m1_nospin:
                    spin_likes *= self.pl_spin_one(i[2], params[6], params[7])
            else:
                spin_likes = 1

        else:
            p_m2 = truncnormal_like_one(i[1], params[0], params[1], lower = 1, upper = m_crit(params[2], i[3]))
            if self.spinning:
                spin_likes = self.pl_spin_one(i[3], params[5], params[6])
                if not self.m1_nospin:
                    spin_likes *= self.pl_spin_one(i[2], params[5], params[6])
            else:
                spin_likes = 1
            # print(spin_likes)
        # print(i, p_m2)
        # print(np.mean(p_m1*p_m2*p_q*spin_likes))
        return np.mean(p_m1*p_m2*p_q*spin_likes)


    def pop_like(self, samples, params):
        # global single
        if self.pop_type == "one":
            if self.selection:
                mu = self.selection_norm(params, injection_set_bns)
                # print(mu)
            else:
                mu = 1

            if self.samples:
                result = np.sum([np.log(self.event_likelihood_one_samples(i, params))/mu for i in samples])
            else:
                result = np.sum([np.log(self.event_likelihood_one_single(i, params))/mu for i in samples])
        elif self.pop_type =="two":
            if self.selection:
                mu = self.selection_norm(params, injection_set_bns)
                # print(mu)
            else:
                mu = 1
            if self.samples:
                result = np.sum([np.log(self.event_likelihood_two_samples(i, params))/mu for i in samples])
            else:
                # x = np.linspace(1.01, 2, 200)
                # likes = two_truncnormal_like(x, params[0], params[1], params[2], params[3], params[4], 1, params[5])
                # print(likes)
                # amu = np.trapz(likes*(x**2.2), x=x)
                # print(amu)
                result = np.sum([np.log(self.event_likelihood_two_single(i, params)/mu) for i in samples])
        elif self.pop_type =="nsbh":
            if self.selection:
                mu = self.selection_norm(params, injection_set)
                # print(mu)
            else:
                mu = 1
            if self.samples:
                result = np.sum([np.log(self.event_likelihood_nsbh_samples(i, params)/mu) for i in samples])
            else:
                result = np.sum([np.log(self.event_likelihood_nsbh_single(i, params)/mu) for i in samples])

        elif self.pop_type =="nsbh_one":
            if self.selection:
                mu = self.selection_norm(params, injection_set)
                # print(mu)
            else:
                mu = 1 # could replace w/ threshold 0
            if self.samples:
                result = np.sum([np.log(self.event_likelihood_nsbh_one_samples(i, params)/mu) for i in samples])
            else:
                # print([np.log(self.event_likelihood_nsbh_one_single(i, params)) for i in samples])
                result = np.sum([np.log(self.event_likelihood_nsbh_one_single(i, params)/mu) for i in samples])

        if math.isnan(result):
            return -np.inf

        return result

    def infer(self, samples, steps, save_to='./default.h5', fixed = {}, mult=False):
        """
        Perform inference on samples.

        samples -- generated samples for this population
        steps -- number of steps
        save_to -- .h5 file backend
        fixed -- a dictionary including values
        """

        if self.pop_type == "one":
            ranges, pscale = fix_params_one(fixed, self.vary_slope, self.spinning)
            # print(ranges)
            # print(pscale)
            def logpost_one(params, data):
                if params[0] > ranges[0,0] and params[0] < ranges[0,1]: # mu
                    if params[1] > ranges[1,0] and params[1] < ranges[1,1]: # sigma
                        if params[2] > ranges[2, 0] and params[2] < ranges[2, 1]: # m_TOV
                            if self.vary_slope:
                                if params[3] > ranges[3, 0] and params[3] < ranges[3, 1]: # slope
                                    if self.spinning:
                                        if params[4] > ranges[4,0] and params[4] < ranges[4,1]: # max jjkep
                                            if params[5] >= ranges[5,0] and params[5] < ranges[5,1]:
                                                return self.pop_like(data, params)
                                    else:
                                        return self.pop_like(data, params)
                                else:
                                    return -np.inf
                            else:
                                if self.spinning:
                                    if params[3] > ranges[3,0] and params[3] < ranges[3,1]: # max jjkep
                                        if params[4] >= ranges[4,0] and params[4] < ranges[4,1]:
                                            return self.pop_like(data, params)
                                    else:
                                        return -np.inf
                                return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope and self.spinning:
                p0 = [self.mu_1, self.sigma_1, self.m_TOV, self.slope, self.max_jjkep, self.spin_slope]
                # pscale = [0.1, 0.05, 0.1, 0.05, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(12, 6)
            elif self.vary_slope:
                p0 = [self.mu_1, self.sigma_1, self.m_TOV, self.slope]
                # pscale = [0.1, 0.05, 0.1, 0.05]
                pos = p0 + pscale*np.random.randn(8, 4)
            elif self.spinning:
                p0 = [self.mu_1, self.sigma_1, self.m_TOV, self.max_jjkep, self.spin_slope]
                # pscale = [0.1, 0.05, 0.1, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(10, 5)
            else:
                p0 = [self.mu_1, self.sigma_1, self.m_TOV]
                # pscale = [0.1, 0.05, 0.1]
                pos = p0 + pscale*np.random.randn(6, 3)

        elif self.pop_type == "nsbh":
            ranges, pscale = fix_params_nsbh(fixed, self.vary_slope, self.spinning)
            print(ranges, pscale)
            def logpost_one(params, data):
                # params = a, mu_1, sigma_1, mu_2, sigma_2, m_TOV, bh_min, bh_slope, slope (optional)
                if params[0] > ranges[0,0] and params[0] < ranges[0,1]: # a
                     #must get diff peaks
                    if params[1] > ranges[1,0] and params[1] < ranges[1,1]: # mu_1
                        if params[2] > ranges[2,0] and params[2] < ranges[2,1]: #sigma_1
                            if params[1] < params[3]:
                                if params[2] < params[4]: # must get diff sdevs
                                    if params[3] > ranges[3,0] and params[3] < ranges[3,1]: # mu_2
                                        if params[4] > ranges[4,0] and params[4] < ranges[4,1]: # sigma_2
                                            if params[5] > ranges[5,0] and params[5] < ranges[5,1]: # m_TOV
                                                if params[6] > ranges[6,0] and params[6] < ranges[6,1]: # bh_min
                                                    if params[7] > ranges[7,0] and params[7] < ranges[7,1]: # bh_slope
                                                        if self.vary_slope:
                                                            if params[8] > ranges[8,0] and params[8] < ranges[8,1]: # slope
                                                                if self.spinning:
                                                                    if params[9] > ranges[9,0] and params[9] < ranges[9,1]: # max jjkep
                                                                        if params[10] >= ranges[10,0] and params[10] < ranges[10,1]:
                                                                            return self.pop_like(data, params)
                                                                else:
                                                                    return self.pop_like(data, params)
                                                            else:
                                                                return -np.inf
                                                        else:
                                                            if self.spinning:
                                                                # print('here')
                                                                if params[8] > ranges[8,0] and params[8] < ranges[8,1]: # max jjkep
                                                                    if params[9] >= ranges[9,0] and params[9] < ranges[9,1]:
                                                                        return self.pop_like(data, params)
                                                                else:
                                                                    return -np.inf
                                                            return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope and self.spinning:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.bh_min, self.bh_slope, self.slope, self.max_jjkep, self.spin_slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.2, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(22, 11)
            elif self.vary_slope:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.bh_min, self.bh_slope, self.slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.2, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(18, 9)
            elif self.spinning:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.bh_min, self.bh_slope, self.max_jjkep, self.spin_slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.2, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(20, 10)
            else:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.bh_min, self.bh_slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.2, 0.2]
                pos = p0 + pscale*np.random.randn(16, 8)
            # print(p0)

        elif self.pop_type == "nsbh_one":
            ranges, pscale = fix_params_nsbh_one(fixed, self.vary_slope, self.spinning)
            print(ranges, pscale)
            def logpost_one(params, data):
                if params[0] > ranges[0,0] and params[0] < ranges[0,1]: # mu
                     #must get diff peaks
                    if params[1] > ranges[1,0] and params[1] < ranges[1,1]: # sigma
                        if params[2] > ranges[2,0] and params[2] < ranges[2,1]:
                            if params[3] > ranges[3,0] and params[3] < ranges[3,1]: #
                                if params[4] > ranges[4,0] and params[4] < ranges[4,1]: # bh_slope
                                    if self.vary_slope:
                                        if params[5] > ranges[5,0] and params[5] < ranges[5,1]: # slope
                                            if self.spinning:
                                                if params[6] > ranges[6,0] and params[6] < ranges[6,1]: # max jjkep
                                                    if params[7] >= ranges[7,0] and params[7] < ranges[7,1]:
                                                        return self.pop_like(data, params)
                                            else:
                                                return self.pop_like(data, params)
                                        else:
                                            return -np.inf
                                    else:
                                        if self.spinning:
                                            # print('here')
                                            if params[5] > ranges[5,0] and params[5] < ranges[5,1]: # max jjkep
                                                if params[6] >= ranges[6,0] and params[6] < ranges[6,1]:
                                                    return self.pop_like(data, params)
                                            else:
                                                return -np.inf
                                        return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope and self.spinning:
                p0 = [self.mu, self.sigma, self.m_TOV, self.bh_min, self.bh_slope, self.slope, self.max_jjkep, self.spin_slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.2, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(16, 8)
            elif self.vary_slope:
                p0 = [self.mu, self.sigma, self.m_TOV, self.bh_min, self.bh_slope, self.slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.2, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(12, 6)
            elif self.spinning:
                p0 = [self.mu, self.sigma, self.m_TOV, self.bh_min, self.bh_slope, self.max_jjkep, self.spin_slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.2, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(14, 7)
            else:
                p0 = [self.mu, self.sigma, self.m_TOV, self.bh_min, self.bh_slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.2, 0.2]
                pos = p0 + pscale*np.random.randn(10, 5)
            # print(p0, pos)

        elif self.pop_type == "two":
            ranges, pscale = fix_params_two(fixed, self.vary_slope, self.spinning)
            print(ranges, pscale)
            def logpost_one(params, data):
                # params = a, mu1, sigma1, mu2, sigma2, m_TOV
                if params[0] > ranges[0,0] and params[0] < ranges[0,1]: # a
                     #must get diff peaks
                    if params[1] > ranges[1,0] and params[1] < ranges[1,1]: # mu_1
                        if params[2] > ranges[2,0] and params[2] < ranges[2,1]: #sigma_1
                            if params[1] < params[3]:
                                if params[2] < params[4]: # must get diff sdevs
                                    if params[3] > ranges[3,0] and params[3] < ranges[3,1]: # mu_2
                                        if params[4] > ranges[4,0] and params[4] < ranges[4,1]: # sigma_2
                                            if params[5] > ranges[5,0] and params[5] < ranges[5,1]: # m_TOV
                                                if self.vary_slope:
                                                    if params[6] > ranges[6,0] and params[6] < ranges[6,1]: # slope
                                                        if self.spinning:
                                                            if params[7] > ranges[7,0] and params[7] < ranges[7,1]: # max jjkep
                                                                if params[8] >= ranges[8,0] and params[8] < ranges[8,1]:
                                                                    return self.pop_like(data, params)
                                                        else:
                                                            return self.pop_like(data, params)
                                                    else:
                                                        return -np.inf
                                                else:
                                                    if self.spinning:
                                                        if params[6] > ranges[6,0] and params[6] < ranges[6,1]: # max jjkep
                                                            if params[7] >= ranges[7,0] and params[7] < ranges[7,1]:
                                                                return self.pop_like(data, params)
                                                        else:
                                                            return -np.inf
                                                    return self.pop_like(data, params)
                return -np.inf

            pscale /= np.sqrt(samples.shape[0])

            if self.vary_slope and self.spinning:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.slope, self.max_jjkep, self.spin_slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(18, 9)
            elif self.vary_slope:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05]
                pos = p0 + pscale*np.random.randn(14, 7)
            elif self.spinning:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.max_jjkep, self.spin_slope]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(16, 8)
            else:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV]
                # pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1]
                pos = p0 + pscale*np.random.randn(12, 6)


        nwalkers, ndim = pos.shape
        if save_to is not None:
            backend = emcee.backends.HDFBackend(save_to)
            backend.reset(nwalkers, ndim)

        if mult:
            with Pool() as pool:
                sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost_one, pool=pool, args=([samples]))
                sampler.run_mcmc(pos, steps, progress=True, skip_initial_state_check=True)
                posterior_samples = sampler.get_chain(discard = 100, flat=True)
                log_prob_samples = sampler.get_log_prob(discard = 100, flat=True)
                pool.close()
                pool.join()

        backend = emcee.backends.HDFBackend(save_to)
        backend.reset(nwalkers, ndim)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost_one, backend=backend, args=([samples]))

        sampler.run_mcmc(pos, steps, progress=True, skip_initial_state_check=True)
        posterior_samples = sampler.get_chain(discard = 100, flat=True)
        log_prob_samples = sampler.get_log_prob(discard = 100, flat=True)

        return posterior_samples, log_prob_samples




class Population_One(Population):
    def __init__(self, params, vary_slope=False, selection=False):
        """
        Population of compact objects, with varying spin.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?
        """
        self.pop_type = "one"
        super().__init__(params, pop_type = self.pop_type, vary_slope = vary_slope, selection = selection)
        print(self.vary_slope)

class Population_One_Spinning(Population):
    def __init__(self, params, vary_slope=False, selection=False, spin_params = [1.0]):
        """
        Population of compact objects, with varying spin.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?
        """
        self.pop_type = "one"
        super().__init__(params, pop_type = self.pop_type, vary_slope = vary_slope, selection = selection, spinning=True, spin_params = spin_params)
        print(self.vary_slope)

class Population_Two(Population):
    def __init__(self, params, vary_slope=False, selection=False):
        """
        Population of compact objects, with varying spin.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?
        """
        self.pop_type = "two"
        super().__init__(params, pop_type = self.pop_type, vary_slope = vary_slope, selection = selection)
        print(self.vary_slope)

class Population_Two_Spinning(Population):
    def __init__(self, params, vary_slope=False, selection=False, spin_params = [1.0]):
        """
        Population of compact objects, with varying spin.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?
        """
        self.pop_type = "two"
        super().__init__(params, pop_type = self.pop_type, vary_slope = vary_slope, selection = selection, spinning=True, spin_params = spin_params)
        print(self.vary_slope)

class Population_NSBH(Population):
    def __init__(self, params, vary_slope=False, selection=False):
        """
        Population of compact objects, with varying spin.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal, "nsbh" for NSBH
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?
        """
        self.pop_type = "nsbh"
        super().__init__(params, pop_type = self.pop_type, vary_slope = vary_slope, selection = selection)
        print(self.vary_slope)

class Population_NSBH_Spinning(Population):
    def __init__(self, params, vary_slope=False, selection=False, spin_params = [1.0]):
        """
        Population of compact objects, with varying spin.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?
        """
        self.pop_type = "nsbh"
        super().__init__(params, pop_type = self.pop_type, vary_slope = vary_slope, selection = selection, spinning=True, spin_params = spin_params)
        # print(self.vary_slope)

def fix_params_one(fixed, vary_slope, spinning):

    count = 3
    if vary_slope:
        count +=1
    if spinning:
        count +=2

    ranges = np.zeros((count, 2))
    pscale = np.zeros(count) + 0.00

    if "mu" in fixed:
        if isinstance(fixed["mu"], list):
            list_ = fixed["mu"]
            ranges[0] = [list_[1], list_[2]]
        else:
            ranges[0] = [fixed["mu"]-0.01, fixed["mu"]+0.01]
    else:
        ranges[0] = [1, 2]
        pscale[0] = 0.1

    if "sigma" in fixed:
        if isinstance(fixed["sigma"], list):
            list_ = fixed["sigma"]
            ranges[1] = [list_[1], list_[2]]
        else:
            ranges[1] = [fixed["sigma"]-0.01, fixed["sigma"]+0.01]
    else:
        ranges[1] = [0.01, 1]
        pscale[1] = 0.1

    if "m_TOV" in fixed:
        if isinstance(fixed["sigma"], list):
            list_ = fixed["sigma"]
            ranges[1] = [list_[1], list_[2]]
        else:
            ranges[2] = [fixed["m_TOV"]-0.01, fixed["m_TOV"]+0.01]
    else:
        ranges[2] = [1.5, 2.5]
        pscale[2] = 0.1

    if vary_slope:
        if "slope" in fixed:
            if isinstance(fixed["slope"], list):
                list_ = fixed["slope"]
                ranges[3] = [list_[1], list_[2]]
            else:
                ranges[3] = [fixed["slope"]-0.01, fixed["slope"]+0.01]
        else:
            ranges[3] = [0, 0.5]
            pscale[3] = 0.05

        if spinning:
            if "max_jjkep" in fixed:
                if isinstance(fixed["max_jjkep"], list):
                    list_ = fixed["max_jjkep"]
                    ranges[4] = [list_[1], list_[2]]
                else:
                    ranges[4] = [fixed["max_jjkep"]-0.01, fixed["max_jjkep"]+0.01]
            else:
                ranges[4] = [0, 1]
                pscale[4] = 0.05

            if "spin_slope" in fixed:
                if isinstance(fixed["spin_slope"], list):
                    list_ = fixed["spin_slope"]
                    ranges[5] = [list_[1], list_[2]]
                else:
                    ranges[5] = [fixed["spin_slope"]-0.01, fixed["spin_slope"]+0.01]
            else:
                ranges[5] = [-0.01, 6]
                pscale[5] = 0.2

    else:
        if spinning:
            if "max_jjkep" in fixed:
                if isinstance(fixed["max_jjkep"], list):
                    list_ = fixed["max_jjkep"]
                    ranges[3] = [list_[1], list_[2]]
                else:
                    ranges[3] = [fixed["max_jjkep"]-0.01, fixed["max_jjkep"]+0.01]
            else:
                ranges[3] = [0, 1]
                pscale[3] = 0.05

            if "spin_slope" in fixed:
                if isinstance(fixed["spin_slope"], list):
                    list_ = fixed["spin_slope"]
                    ranges[4] = [list_[1], list_[2]]
                else:
                    ranges[4] = [fixed["spin_slope"]-0.01, fixed["spin_slope"]+0.01]
            else:
                ranges[4] = [-0.01, 6]
                pscale[4] = 0.2

    return ranges, pscale

def fix_params_two(fixed, vary_slope, spinning):

    count = 6
    if vary_slope:
        count +=1
    if spinning:
        count +=2

    ranges = np.zeros((count, 2))
    pscale = np.zeros(count) + 0.00

    if "a" in fixed:
        if isinstance(fixed["a"], list):
            list_ = fixed["a"]
            ranges[0] = [list_[1], list_[2]]
        else:
            ranges[0] = [fixed["a"]-0.01, fixed["a"]+0.01]
    else:
        ranges[0] = [0, 1]
        pscale[0] = 0.05

    if "mu_1" in fixed:
        if isinstance(fixed["mu_1"], list):
            list_ = fixed["mu_1"]
            ranges[1] = [list_[1], list_[2]]
        else:
            ranges[1] = [fixed["mu_1"]-0.01, fixed["mu_1"]+0.01]
    else:
        ranges[1] = [1, 2]
        pscale[1] = 0.1

    if "sigma_1" in fixed:
        if isinstance(fixed["sigma_1"], list):
            list_ = fixed["sigma_1"]
            ranges[2] = [list_[1], list_[2]]
        else:
            ranges[2] = [fixed["sigma_1"]-0.01, fixed["sigma_1"]+0.01]
    else:
        ranges[2] = [0.01, 0.5]
        pscale[2] = 0.05

    if "mu_2" in fixed:
        if isinstance(fixed["mu_2"], list):
            list_ = fixed["mu_2"]
            ranges[3] = [list_[1], list_[2]]
        else:
            ranges[3] = [fixed["mu_2"]-0.01, fixed["mu_2"]+0.01]
    else:
        ranges[3] = [1, 2]
        pscale[3] = 0.1

    if "sigma_2" in fixed:
        if isinstance(fixed["sigma_2"], list):
            list_ = fixed["sigma_2"]
            ranges[4] = [list_[1], list_[2]]
        else:
            ranges[4] = [fixed["sigma_2"]-0.01, fixed["sigma_2"]+0.01]
    else:
        ranges[4] = [0.01, 0.5]
        pscale[4] = 0.05

    if "m_TOV" in fixed:
        if isinstance(fixed["m_TOV"], list):
            list_ = fixed["m_TOV"]
            ranges[5] = [list_[1], list_[2]]
        else:
            ranges[5] = [fixed["m_TOV"]-0.01, fixed["m_TOV"]+0.01]
    else:
        ranges[5] = [1.5, 2.5]
        pscale[5] = 0.1

    if vary_slope:
        if "slope" in fixed:
            if isinstance(fixed["slope"], list):
                list_ = fixed["slope"]
                ranges[6] = [list_[1], list_[2]]
            else:
                ranges[6] = [fixed["slope"]-0.01, fixed["slope"]+0.01]
        else:
            ranges[6] = [0, 0.5]
            pscale[6] = 0.05

        if spinning:
            if "max_jjkep" in fixed:
                if isinstance(fixed["max_jjkep"], list):
                    list_ = fixed["max_jjkep"]
                    ranges[7] = [list_[1], list_[2]]
                else:
                    ranges[7] = [fixed["max_jjkep"]-0.01, fixed["max_jjkep"]+0.01]
            else:
                ranges[7] = [0, 1]
                pscale[7] = 0.05

            if "spin_slope" in fixed:
                if isinstance(fixed["spin_slope"], list):
                    list_ = fixed["spin_slope"]
                    ranges[8] = [list_[1], list_[2]]
                else:
                    ranges[8] = [fixed["spin_slope"]-0.01, fixed["spin_slope"]+0.01]
            else:
                ranges[8] = [-0.01, 6]
                pscale[8] = 0.2

    else:
        if spinning:
            if "max_jjkep" in fixed:
                if isinstance(fixed["max_jjkep"], list):
                    list_ = fixed["max_jjkep"]
                    ranges[6] = [list_[1], list_[2]]
                else:
                    ranges[6] = [fixed["max_jjkep"]-0.01, fixed["max_jjkep"]+0.01]
            else:
                ranges[6] = [0, 1]
                pscale[6] = 0.05

            if "spin_slope" in fixed:
                if isinstance(fixed["spin_slope"], list):
                    list_ = fixed["spin_slope"]
                    ranges[7] = [list_[1], list_[2]]
                else:
                    ranges[7] = [fixed["spin_slope"]-0.01, fixed["spin_slope"]+0.01]
            else:
                ranges[7] = [-0.01, 6]
                pscale[7] = 0.2

    return ranges, pscale

def fix_params_nsbh(fixed, vary_slope, spinning):

    count = 8
    if vary_slope:
        count +=1
    if spinning:
        count +=2

    ranges = np.zeros((count, 2))
    pscale = np.zeros(count) + 0.00

    if "a" in fixed:
        if isinstance(fixed["a"], list):
            list_ = fixed["a"]
            ranges[0] = [list_[1], list_[2]]
        else:
            ranges[0] = [fixed["a"]-0.01, fixed["a"]+0.01]
    else:
        ranges[0] = [0, 1]
        pscale[0] = 0.05

    if "mu_1" in fixed:
        if isinstance(fixed["mu_1"], list):
            list_ = fixed["mu_1"]
            ranges[1] = [list_[1], list_[2]]
        else:
            ranges[1] = [fixed["mu_1"]-0.01, fixed["mu_1"]+0.01]
    else:
        ranges[1] = [1, 2]
        pscale[1] = 0.1

    if "sigma_1" in fixed:
        if isinstance(fixed["sigma_1"], list):
            list_ = fixed["sigma_1"]
            ranges[2] = [list_[1], list_[2]]
        else:
            ranges[2] = [fixed["sigma_1"]-0.01, fixed["sigma_1"]+0.01]
    else:
        ranges[2] = [0.01, 0.5]
        pscale[2] = 0.05

    if "mu_2" in fixed:
        if isinstance(fixed["mu_2"], list):
            list_ = fixed["mu_2"]
            ranges[3] = [list_[1], list_[2]]
        else:
            ranges[3] = [fixed["mu_2"]-0.01, fixed["mu_2"]+0.01]
    else:
        ranges[3] = [1, 2]
        pscale[3] = 0.1

    if "sigma_2" in fixed:
        if isinstance(fixed["sigma_2"], list):
            list_ = fixed["sigma_2"]
            ranges[4] = [list_[1], list_[2]]
        else:
            ranges[4] = [fixed["sigma_2"]-0.01, fixed["sigma_2"]+0.01]
    else:
        ranges[4] = [0.01, 0.5]
        pscale[4] = 0.05

    if "m_TOV" in fixed:
        if isinstance(fixed["m_TOV"], list):
            list_ = fixed["m_TOV"]
            ranges[5] = [list_[1], list_[2]]
        else:
            ranges[5] = [fixed["m_TOV"]-0.01, fixed["m_TOV"]+0.01]
    else:
        ranges[5] = [1.5, 2.5]
        pscale[5] = 0.1

    if "bh_min" in fixed:
        if isinstance(fixed["bh_min"], list):
            list_ = fixed["bh_min"]
            ranges[6] = [list_[1], list_[2]]
        else:
            ranges[6] = [fixed["bh_min"]-0.01, fixed["bh_min"]+0.01]
    else:
        ranges[6] = [3, 10]
        pscale[6] = 0.2

    if "bh_slope" in fixed:
        if isinstance(fixed["bh_slope"], list):
            list_ = fixed["bh_slope"]
            ranges[7] = [list_[1], list_[2]]
        else:
            ranges[7] = [fixed["bh_slope"]-0.01, fixed["bh_slope"]+0.01]
    else:
        ranges[7] = [1, 6]
        pscale[7] = 0.2

    if vary_slope:
        if "slope" in fixed:
            if isinstance(fixed["slope"], list):
                list_ = fixed["slope"]
                ranges[8] = [list_[1], list_[2]]
            else:
                ranges[8] = [fixed["slope"]-0.01, fixed["slope"]+0.01]
        else:
            ranges[8] = [0, 0.5]
            pscale[8] = 0.05

        if spinning:
            if "max_jjkep" in fixed:
                if isinstance(fixed["max_jjkep"], list):
                    list_ = fixed["max_jjkep"]
                    ranges[9] = [list_[1], list_[2]]
                else:
                    ranges[9] = [fixed["max_jjkep"]-0.01, fixed["max_jjkep"]+0.01]
            else:
                ranges[9] = [0, 1]
                pscale[9] = 0.05

            if "spin_slope" in fixed:
                if isinstance(fixed["spin_slope"], list):
                    list_ = fixed["spin_slope"]
                    ranges[10] = [list_[1], list_[2]]
                else:
                    ranges[10] = [fixed["spin_slope"]-0.01, fixed["spin_slope"]+0.01]
            else:
                ranges[10] = [-0.01, 6]
                pscale[10] = 0.2

    else:
        if spinning:
            if "max_jjkep" in fixed:
                if isinstance(fixed["max_jjkep"], list):
                    list_ = fixed["max_jjkep"]
                    ranges[8] = [list_[1], list_[2]]
                else:
                    ranges[8] = [fixed["max_jjkep"]-0.01, fixed["max_jjkep"]+0.01]
            else:
                ranges[8] = [0, 1]
                pscale[8] = 0.05

            if "spin_slope" in fixed:
                if isinstance(fixed["spin_slope"], list):
                    list_ = fixed["spin_slope"]
                    ranges[9] = [list_[1], list_[2]]
                else:
                    ranges[9] = [fixed["spin_slope"]-0.01, fixed["spin_slope"]+0.01]
            else:
                ranges[9] = [-0.01, 6]
                pscale[9] = 0.2

    return ranges, pscale


def fix_params_nsbh_one(fixed, vary_slope, spinning):

    count = 5
    if vary_slope:
        count +=1
    if spinning:
        count +=2

    ranges = np.zeros((count, 2))
    pscale = np.zeros(count) + 0.00

    if "mu" in fixed:
        if isinstance(fixed["mu"], list):
            list_ = fixed["mu"]
            ranges[0] = [list_[1], list_[2]]
        else:
            ranges[0] = [fixed["mu"]-0.01, fixed["mu"]+0.01]
    else:
        ranges[0] = [1, 2]
        pscale[0] = 0.1

    if "sigma" in fixed:
        if isinstance(fixed["sigma"], list):
            list_ = fixed["sigma"]
            ranges[1] = [list_[1], list_[2]]
        else:
            ranges[1] = [fixed["sigma"]-0.01, fixed["sigma"]+0.01]
    else:
        ranges[1] = [0.01, 1]
        pscale[1] = 0.1

    if "m_TOV" in fixed:
        if isinstance(fixed["sigma"], list):
            list_ = fixed["sigma"]
            ranges[1] = [list_[1], list_[2]]
        else:
            ranges[2] = [fixed["m_TOV"]-0.01, fixed["m_TOV"]+0.01]
    else:
        ranges[2] = [1.5, 2.5]
        pscale[2] = 0.1

    if "bh_min" in fixed:
        if isinstance(fixed["bh_min"], list):
            list_ = fixed["bh_min"]
            ranges[3] = [list_[1], list_[2]]
        else:
            ranges[3] = [fixed["bh_min"]-0.01, fixed["bh_min"]+0.01]
    else:
        ranges[3] = [3, 10]
        pscale[3] = 0.2

    if "bh_slope" in fixed:
        if isinstance(fixed["bh_slope"], list):
            list_ = fixed["bh_slope"]
            ranges[4] = [list_[1], list_[2]]
        else:
            ranges[4] = [fixed["bh_slope"]-0.01, fixed["bh_slope"]+0.01]
    else:
        ranges[4] = [1, 6]
        pscale[4] = 0.2

    if vary_slope:
        if "slope" in fixed:
            if isinstance(fixed["slope"], list):
                list_ = fixed["slope"]
                ranges[5] = [list_[1], list_[2]]
            else:
                ranges[5] = [fixed["slope"]-0.01, fixed["slope"]+0.01]
        else:
            ranges[5] = [0, 0.5]
            pscale[5] = 0.05
        if spinning:
            if "max_jjkep" in fixed:
                if isinstance(fixed["max_jjkep"], list):
                    list_ = fixed["max_jjkep"]
                    ranges[6] = [list_[1], list_[2]]
                else:
                    ranges[6] = [fixed["max_jjkep"]-0.01, fixed["max_jjkep"]+0.01]
            else:
                ranges[6] = [0, 1]
                pscale[6] = 0.05

            if "spin_slope" in fixed:
                if isinstance(fixed["spin_slope"], list):
                    list_ = fixed["spin_slope"]
                    ranges[7] = [list_[1], list_[2]]
                else:
                    ranges[7] = [fixed["spin_slope"]-0.01, fixed["spin_slope"]+0.01]
            else:
                ranges[7] = [-0.01, 6]
                pscale[7] = 0.2

    else:
        if spinning:
            if "max_jjkep" in fixed:
                if isinstance(fixed["max_jjkep"], list):
                    list_ = fixed["max_jjkep"]
                    ranges[5] = [list_[1], list_[2]]
                else:
                    ranges[5] = [fixed["max_jjkep"]-0.01, fixed["max_jjkep"]+0.01]
            else:
                ranges[5] = [0, 1]
                pscale[5] = 0.05

            if "spin_slope" in fixed:
                if isinstance(fixed["spin_slope"], list):
                    list_ = fixed["spin_slope"]
                    ranges[6] = [list_[1], list_[2]]
                else:
                    ranges[6] = [fixed["spin_slope"]-0.01, fixed["spin_slope"]+0.01]
            else:
                ranges[6] = [-0.01, 6]
                pscale[6] = 0.2

    return ranges, pscale
