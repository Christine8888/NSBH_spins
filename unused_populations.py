import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import math
import scipy
import emcee
import corner

# utility functions

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
    result[np.logical_and(x<=upper, x>=lower)] = phi(x[np.logical_and(x<=upper, x>=lower)],mu,sigma)/(normal_cdf(upper,mu,sigma)-normal_cdf(lower,mu,sigma))
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
    result[np.logical_and(x<=upper, x>=lower)] = a*phi(x[np.logical_and(x<=upper, x>=lower)],mu_1,sigma_1)/(normal_cdf(upper,mu_1,sigma_1)-normal_cdf(lower,mu_1,sigma_1))+(1-a)*phi(x[np.logical_and(x<=upper, x>=lower)],mu_2,sigma_2)/(normal_cdf(upper,mu_2,sigma_2)-normal_cdf(lower,mu_2,sigma_2))
    result[np.logical_or(x>upper, x<lower)] = 0
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
    draws = np.random.rand(N*100, 2)
    draws = draws * ((m_1**(beta+1))-(m_min**(beta+1))) # first scale to within (0,1)
    draws += m_min**(beta+1) # then shift so m_min is 0
    draws = draws**(1/(beta+1)) # then complete inverse transform
    return draws[:N, 0]

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

def generate_plminmax(N, x_min, x_max, alpha):
    rand = np.random.rand(N)
    rand *= pl_cdf(x_max, x_min, alpha)
    pl_val = x_min * (1 - rand) ** (-1 / (alpha - 1))
    return pl_val




class Population_SpinDist(Population):
    def __init__(self, params, pop_type, vary_slope=False, selection = False, spin_params=[1.0]):
        """
        Population of compact objects, with varying spin.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal, "nsbh" for NSBH
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?
        selection -- include selection effects
        spin_params -- maximum j/j_kep

        """

        super().__init__(params, pop_type = pop_type, vary_slope = vary_slope, selection = selection)
        self.max_jjkep = spin_params[0]

    ## CAN ADD MORE LIKELIHOODS FOR OTHER SPIN DISTRIBUTIONS!

    def uniform_spin_one(self, spin, max_jjkep):
        if spin < max_jjkep:
            return 1/max_jjkep
        return 0

    def uniform_spin(self, spin, max_jjkep):
        result = np.zeros(spin.shape[0])
        result[spin < max_jjkep] = 1/max_jjkep
        result[spin >= max_jjkep] = 0
        return result

    def event_likelihood_one_samples(self, samples, params):
        # params: mu, sigma, m_TOV, slope (optional), max_jjkep
        if self.vary_slope:
            truncs = truncnormal_like(samples[:,0], params[0], params[1], lower = 1, upper = m_crit_slope(params[2], params[3], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:,0], m_crit_slope(params[2], params[3], samples[:,3])], axis=0))
            spin_likes = self.uniform_spin(samples[:,2], params[4]) * self.uniform_spin(samples[:,3], params[4])
        else:
            truncs = truncnormal_like(samples[:,0], params[0], params[1], lower = 1, upper = m_crit(params[2], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:,0], m_crit(params[2], samples[:,3])], axis=0))
            spin_likes = self.uniform_spin(samples[:,2], params[3]) * self.uniform_spin(samples[:,3], params[3])
        return np.mean(truncs*qlikes*spin_likes)

    def event_likelihood_two_samples(self, samples, params):
        # params: a, mu_1, sigma_1, mu_2, sigma_2, m_TOV, slope (optional), max_jjkep
        if self.vary_slope:
            truncs = two_truncnormal_like(samples[:,0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[6], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:, 0], m_crit_slope(params[5], params[6], samples[:, 3])], axis=0))
            spin_likes = self.uniform_spin(samples[:,2], params[7]) * self.uniform_spin(samples[:,3], params[7])
        else:
            truncs = two_truncnormal_like(samples[:,0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:, 0], m_crit(params[5], samples[:, 3])], axis=0))
            spin_likes = self.uniform_spin(samples[:,2], params[6]) * self.uniform_spin(samples[:,3], params[6])
        return np.mean(truncs*qlikes*spin_likes)

    def event_likelihood_nsbh_samples(self, samples, params):
        # params: 0 a, 1 mu_1, 2 sigma_1, 3 mu_2, 4 sigma_2, 5 m_TOV, 6 bh_min, 7 bh_slope, 8 slope (optional), 8 or 9 max_jjkep
        if self.vary_slope:
            p_m1 = like_plmin(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[8], samples[:,3]))

            q = samples[:,1]/samples[:,0]
            p_q = like_beta(q, self.beta)
            spin_likes = self.uniform_spin(samples[:,2], params[9]) * self.uniform_spin(samples[:,3], params[9])
        else:
            p_m1 = like_plmin(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], samples[:,3]))
            q = samples[:,1]/samples[:,0]
            p_q = like_beta(q, self.beta)
            spin_likes = self.uniform_spin(samples[:,2], params[8]) * self.uniform_spin(samples[:,3], params[8])
        return np.mean(p_m1*p_m2*p_q*spin_likes)

    def event_likelihood_one_single(self, samples, params):
        i = samples[0]
        if self.vary_slope:
            truncs = truncnormal_like_one(i[0], params[0], params[1], lower = 1, upper = m_crit_slope(params[2], params[3], i[2]))
            qlikes = like_m2(i[1], min([i[0], m_crit_slope(params[2], params[3], i[3])]))
            spin_likes = self.uniform_spin_one(samples[:,2], params[4]) * self.uniform_spin_one(samples[:,3], params[4])
        else:
            truncs = np.array([truncnormal_like_one(i[0], params[0], params[1], lower = 1, upper = m_crit(params[2], i[2])) \
                           for i in samples])
            qlikes = np.array([like_m2(i[1], min([i[0], m_crit(params[2], i[3])])) for i in samples])
            spin_likes = self.uniform_spin_one(samples[:,2], params[3]) * self.uniform_spin_one(samples[:,3], params[3])
        return np.mean(truncs*qlikes*spin_likes)

    def event_likelihood_two_single(self, samples, params):
        i = samples[0]
        if self.vary_slope:
            truncs = two_truncnormal_like_one(x = i[0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                            mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[6], i[2]))
            qlikes = like_m2(i[1], min([i[0], m_crit_slope(params[5], params[6], i[3])]))
            spin_likes = self.uniform_spin_one(samples[:,2], params[7]) * self.uniform_spin_one(samples[:,3], params[7])

        else:
            truncs = two_truncnormal_like_one(x = i[0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                            mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], i[2]))
            qlikes = like_m2(i[1], min([i[0], m_crit(params[5], i[3])]))
            spin_likes = self.uniform_spin_one(samples[:,2], params[6]) * self.uniform_spin_one(samples[:,3], params[6])
        return np.mean(truncs*qlikes*spin_likes)

    def event_likelihood_nsbh_single(self, samples, params):
        print('here')
        i = samples[0]
        # params: a, mu_1, sigma_1, mu_2, sigma_2, m_TOV, bh_min, bh_slope, slope (optional)
        # REPLACE LIKELIHOODS!
        if self.vary_slope:
            p_m1 = like_plmin_one(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like_one(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[8], samples[:,3]))

            q = samples[:,1]/samples[:,0]
            p_q = like_beta_one(q, self.beta)
            spin_likes = self.uniform_spin(samples[:,2], params[9]) * self.uniform_spin(samples[:,3], params[9])
        else:
            p_m1 = like_plmin_one(samples[:,0], params[6], params[7])
            print(p_m1)
            p_m2 = two_truncnormal_like_one(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], samples[:,3]))
            q = samples[:,1]/samples[:,0]
            p_q = like_beta_one(q, self.beta)
            spin_likes = self.uniform_spin(samples[:,2], params[8]) * self.uniform_spin(samples[:,3], params[8])
        return np.mean(p_m1*p_m2*p_q*spin_likes)


    def infer(self, samples, steps, save_to=None):
        """
        Perform inference on samples.

        samples -- generated samples for this population
        steps -- number of steps
        save_to -- .h5 file backend
        """
        if self.pop_type == "one":

            def logpost_one(params, data):
                if params[0] > np.min(data) and params[0] < np.max(data): # mu
                    if params[1] < 1 and params[1] > 0.01: # sigma
                        if params[2] > 1.5 and params[2] < 2.5: # m_TOV
                            if self.vary_slope:
                                if params[3] > 0 and params[3] < 0.5: # slope
                                    if params[4] > 0 and params[4] < 1: # max jjkep
                                        return self.pop_like(data, params)
                            else:
                                if params[3] > 0 and params[3] < 1: # max jjkep
                                    return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope:
                p0 = [self.mu_1, self.sigma_1, self.m_TOV, self.slope, self.max_jjkep]
                pscale = [0.1, 0.05, 0.1, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(10, 5)
            else:
                p0 = [self.mu_1, self.sigma_1, self.m_TOV, self.max_jjkep]
                pscale = [0.1, 0.05, 0.1, 0.05]
                pos = p0 + pscale*np.random.randn(8, 4)

        elif self.pop_type == "two":

            def logpost_one(params, data):
                # params = a, mu1, sigma1, mu2, sigma2, m_TOV
                if params[0] < 1 and params[0] > 0: # a
                    if params[1] < params[3]: # must get diff peaks
                        if params[2] < params[4]: # must get diff sdevs
                            if params[1] < 2 and params[1] > 1: # mu 1
                                if params[2] > 0.01 and params[2] < 0.5: # sigma 1
                                    if params[3] > 1 and params[3] < 2: # mu 2
                                        if params[4] > 0.01 and params[4] < 0.5: # sigma 2
                                            if params[5] > 1.5 and params[5] < 2.5: # m_TOV
                                                if self.vary_slope:
                                                    if params[6] > 0 and params[6] < 0.5: # slope
                                                        if params[7] > 0 and params[7] < 1: # max_jjkep
                                                            return self.pop_like(data, params)
                                                else:
                                                    if params[6] > 0 and params[6] < 1: # max_jjkep
                                                        return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.slope+0.2, self.max_jjkep]
                pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(16, 8)
            else:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.max_jjkep]
                pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05]
                pos = p0 + pscale*np.random.randn(14, 7)

        elif self.pop_type == "nsbh":

            def logpost_one(params, data):
                # params: 0 a, 1 mu_1, 2 sigma_1, 3 mu_2, 4 sigma_2, 5 m_TOV, 6 bh_min, 7 bh_slope, 8 slope (optional), 8 or 9 max_jjkep
                # REPLACE LIKELIHOODS!

                if params[0] < 1 and params[0] > 0: # a
                    if params[1] < params[3]: # must get diff peaks
                        if params[2] < params[4]: # must get diff sdevs
                            if params[1] < 2 and params[1] > 1: # mu 1
                                if params[2] > 0.01 and params[2] < 0.5: # sigma 1
                                    if params[3] > 1 and params[3] < 2: # mu 2
                                        if params[4] > 0.01 and params[4] < 0.5: # sigma 2
                                            if params[5] > 1.5 and params[5] < 2.5: # m_TOV
                                                if params[6] > params[5] and params[6] < params[7]: # bh_min
                                                    if params[7] > 1 and params[7] < 6: # bh_slope
                                                        if self.vary_slope:
                                                            if params[8] > 0 and params[8] < 0.5: # slope
                                                                if params[9] > 0 and params[9] < 1: # max_jjkep
                                                                    return self.pop_like(data, params)
                                                        else:
                                                            if params[8] > 0 and params[8] < 1: # max_jjkep
                                                                return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.bh_min, self.bh_slope, self.max_jjkep, self.slope+0.2]
                pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.2, 0.2, 0.05, 0.05]
                pos = p0 + pscale*np.random.randn(20, 10)
            else:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.bh_min, self.bh_slope, self.max_jjkep]
                pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.2, 0.2, 0.05]
                pos = p0 + pscale*np.random.randn(18, 9)


        nwalkers, ndim = pos.shape
        if save_to is not None:
            backend = emcee.backends.HDFBackend(filename)
            backend.reset(nwalkers, ndim)
        print(pos)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost_one, args=([samples]))
        sampler.run_mcmc(pos, steps, progress=True)
        posterior_samples = sampler.get_chain(discard = 100, flat=True)
        log_prob_samples = sampler.get_log_prob(discard = 100, flat=True)

        return posterior_samples, log_prob_samples
