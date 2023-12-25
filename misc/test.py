class Population():
    def __init__(self, params, pop_type, vary_slope=False, selection=False):
        """
        Population of compact objects.

        Keyword arguments:
        N -- number of binaries to generate
        params -- list of population parameters AND hyperparameters
        pop_type -- "one" for single, "two" for bimodal, "nsbh" for NSBH
        samples -- boolean, return samples or points?
        N_samples -- number of samples, if applicable
        vary_slope -- measure slope of spin/mass relationship?

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
        self.selection = selection
        print('selection: {}'.format(self.selection))

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

        test_chi_1 = np.random.rand() * self.max_jjkep
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
            if self.selection:
                if not np.random.rand() < (test_m_1**2.2)/(m_crit(self.m_TOV, self.max_jjkep)**2.2):
                    return None


        elif self.pop_type == "two":
            if self.vary_slope:
                test_m_1 = float(generate_two_truncnormal(1, self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, \
                                                          lower=1, upper=m_crit_slope(self.m_TOV, self.slope, test_chi_1)))
                test_m_2 = float(generate_q(1, self.beta, np.min([test_m_1, m_crit_slope(self.m_TOV, self.slope, test_chi_2)]), m_min=1))
            else:
                test_m_1 = float(generate_two_truncnormal(1, self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, \
                                                          lower=1, upper=m_crit(self.m_TOV, test_chi_1)))
                test_m_2 = float(generate_q(1, self.beta, np.min([test_m_1, m_crit(self.m_TOV, test_chi_2)]), m_min=1))
            if self.selection:
                if not np.random.rand() < (test_m_1**2.2)/(m_crit(self.m_TOV, self.max_jjkep)**2.2):
                    return None

        elif self.pop_type == "nsbh":
            if self.vary_slope:
                params = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.max_jjkep, self.beta, \
                      self.bh_min, self.bh_slope, self.slope]
            else:
                params = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.max_jjkep, self.beta, \
                      self.bh_min, self.bh_slope]

            pop = generate_NSBH(1, params, nsbh_only = True, vary_slope = self.vary_slope)
            # print(pop)
            test_m_1 = pop[0,0]
            test_m_2 = pop[0,1]
            test_chi_1 = pop[0,2]
            test_chi_2 = pop[0,3]

            if self.selection:
                if not np.random.rand() < (test_m_1**2.2)/(10**2.2): # reminder that you're cutting it at 20!
                    return None

        test_rho = get_rho(1, 4, 10)

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

        if return_p0:
            return p0

        pscale = [0.1, 0.1, 0.1, 0.1]
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
            if params[0] > 0.5 and params[0] < 3.5:
                if params[1] > 0.5 and params[1] < 3.5:
                    if params[0] > params[1]:
                        if params[2] > 0 and params[2] < 1:
                            if params[3] > 0 and params[3] < 1:
                                loglike = loglike_one(params, m_chirp, m_chirp_sigma, mass_ratio, mass_ratio_sigma, chi_eff, chi_eff_sigma)
                                return loglike
            return -np.inf

        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost_one, args=(test_mchirp, test_mchirp_sigma, test_ratio, test_ratio_sigma, test_chieff, test_chieff_sigma))
        sampler.run_mcmc(pos, n, progress=False)
        #t_auto = np.mean(sampler.get_autocorr_time())
        samples = sampler.get_chain(discard = int(n/4), flat=True, thin=1)

        return samples

    def event_likelihood_one_samples(self, samples, params):
        # params: mu, sigma, m_TOV, slope (optional)
        if self.vary_slope:
            truncs = truncnormal_like(samples[:,0], params[0], params[1], lower = 1, upper = m_crit_slope(params[2], params[3], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:,0], m_crit_slope(params[2], params[3], samples[:,3])], axis=0))
        else:
            truncs = truncnormal_like(samples[:,0], params[0], params[1], lower = 1, upper = m_crit(params[2], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:,0], m_crit(params[2], samples[:,3])], axis=0))
        return np.mean(truncs*qlikes)

    def event_likelihood_two_samples(self, samples, params):
        # params: 0 a, 1 mu_1, 2 sigma_1, 3 mu_2, 4 sigma_2, 5 m_TOV, 6 slope (optional)
        if self.vary_slope:
            truncs = two_truncnormal_like(samples[:,0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[6], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:, 0], m_crit_slope(params[5], params[6], samples[:, 3])], axis=0))
        else:
            truncs = two_truncnormal_like(samples[:,0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], samples[:,2]))

            qlikes = like_m2(samples[:,1], np.min([samples[:, 0], m_crit(params[5], samples[:, 3])], axis=0))
        return np.mean(truncs*qlikes)

    def event_likelihood_nsbh_samples(self, samples, params):
        # params: a, mu_1, sigma_1, mu_2, sigma_2, m_TOV, bh_min, bh_slope, slope (optional)
        # REPLACE LIKELIHOODS!
        if self.vary_slope:
            p_m1 = like_plmin(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[8], samples[:,3]))

            q = samples[:,1]/samples[:,0]
            p_q = like_beta(q, self.beta)
        else:
            p_m1 = like_plmin(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], samples[:,3]))
            q = samples[:,1]/samples[:,0]
            p_q = like_beta(q, self.beta)
        return np.mean(p_m1*p_m2*p_q)

    def event_likelihood_one_single(self, samples, params):
        if self.vary_slope:
            truncs = np.array([truncnormal_like_one(i[0], params[0], params[1], lower = 1, upper = m_crit_slope(params[2], params[3], i[2])) \
                           for i in samples])
            qlikes = np.array([like_m2(i[1], min([i[0], m_crit_slope(params[2], params[3], i[3])])) for i in samples])
        else:
            truncs = np.array([truncnormal_like_one(i[0], params[0], params[1], lower = 1, upper = m_crit(params[2], i[2])) \
                           for i in samples])
            qlikes = np.array([like_m2(i[1], min([i[0], m_crit(params[2], i[3])])) for i in samples])
        return np.mean(truncs*qlikes)

    def event_likelihood_two_single(self, samples, params):
        if self.vary_slope:
            truncs = np.array([two_truncnormal_like_one(x = i[0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                            mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[6], i[2])) \
                           for i in samples])
            qlikes = np.array([like_m2(i[1], min([i[0], m_crit_slope(params[5], params[6], i[3])])) for i in samples])
        else:
            truncs = np.array([two_truncnormal_like_one(x = i[0], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                            mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], i[2])) \
                           for i in samples])
            qlikes = np.array([like_m2(i[1], min([i[0], m_crit(params[5], i[3])])) for i in samples])
        return np.mean(truncs*qlikes)

    def event_likelihood_nsbh_single(self, samples, params):
        # params: a, mu_1, sigma_1, mu_2, sigma_2, m_TOV, bh_min, bh_slope, slope (optional)
        # REPLACE LIKELIHOODS!
        if self.vary_slope:
            p_m1 = like_plmin_one(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like_one(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit_slope(params[5], params[8], samples[:,3]))

            q = samples[:,1]/samples[:,0]
            p_q = like_beta_one(q, self.beta)
        else:
            p_m1 = like_plmin_one(samples[:,0], params[6], params[7])
            p_m2 = two_truncnormal_like_one(samples[:,1], a = params[0], mu_1 = params[1], sigma_1 = params[2], \
                                      mu_2 = params[3], sigma_2 = params[4], lower = 1, upper = m_crit(params[5], samples[:,3]))
            q = samples[:,1]/samples[:,0]
            p_q = like_beta_one(q, self.beta)
        return np.mean(p_m1*p_m2*p_q)


    def pop_like(self, samples, params):
        if self.pop_type == "one":
            if self.samples:
                result = np.sum([np.log(self.event_likelihood_one_samples(i, params)) for i in samples])
            else:
                result = np.sum([np.log(self.event_likelihood_one_single(i, params)) for i in samples])
        elif self.pop_type =="two":

            if self.samples:
                result = np.sum([np.log(self.event_likelihood_two_samples(i, params)) for i in samples])
            else:
                x = np.linspace(1.01, 2, 200)
                likes = two_truncnormal_like(x, params[0], params[1], params[2], params[3], params[4], 1, params[5])
                # print(likes)
                amu = np.trapz(likes*(x**2.2), x=x)
                # print(amu)
                result = np.sum([np.log(self.event_likelihood_two_single(i, params)/amu) for i in samples])
        elif self.pop_type =="nsbh":
            if self.samples:
                result = np.sum([np.log(self.event_likelihood_nsbh_samples(i, params)) for i in samples])
            else:
                result = np.sum([np.log(self.event_likelihood_nsbh_single(i, params)) for i in samples])

        if math.isnan(result):
            return -np.inf

        return result

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
                                if params[3] > 0 and params[3] < 0.5: #
                                    return self.pop_like(data, params)
                                else:
                                    return -np.inf
                            else:
                                return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope:
                p0 = [self.mu_1, self.sigma_1, self.m_TOV, self.slope]
                pscale = [0.1, 0.05, 0.1, 0.05]
                pos = p0 + pscale*np.random.randn(8, 4)
            else:
                p0 = [self.mu_1, self.sigma_1, self.m_TOV]
                pscale = [0.1, 0.05, 0.1]
                pos = p0 + pscale*np.random.randn(6, 3)

        elif self.pop_type == "nsbh":

            def logpost_one(params, data):
                # params = a, mu_1, sigma_1, mu_2, sigma_2, m_TOV, bh_min, bh_slope, slope (optional)
                if params[0] < 1 and params[0] > 0: # a
                    if params[1] < params[3]: #must get diff peaks
                        if params[1] < 2 and params[1] > 1: # mu_1
                            if params[2] > 0.01 and params[2] < 0.5: #sigma_1
                                if params[2] < params[4]: # must get diff sdevs
                                    if params[3] > 1 and params[3] < 2: # mu_2
                                        if params[4] > 0.01 and params[4] < 0.5: # sigma_2
                                            if params[5] > 1.5 and params[5] < 2.5: # m_TOV
                                                if params[6] > params[5] and params[6] < params[7]: # bh_min
                                                    if params[7] > 1 and params[7] < 6: # bh_slope
                                                        if self.vary_slope:
                                                            if params[8] > 0 and params[8] < 0.5: # slope
                                                                return self.pop_like(data, params)
                                                            else:
                                                                return -np.inf
                                                        else:
                                                            return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.bh_min, self.bh_slope, self.slope]
                pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05, 0.2, 0.2]
                pos = p0 + pscale*np.random.randn(18, 9)
            else:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.bh_min, self.bh_slope]
                pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.2, 0.2]
                pos = p0 + pscale*np.random.randn(16, 8)

        elif self.pop_type == "two":

            def logpost_one(params, data):
                # params = a, mu1, sigma1, mu2, sigma2, m_TOV
                if params[0] < 1 and params[0] > 0: # a
                    if params[1] < params[3]: #must get diff peaks
                        if params[1] < 2 and params[1] > 1: # mu_1
                            if params[2] > 0.01 and params[2] < 0.5: #sigma_1
                                if params[2] < params[4]: # must get diff sdevs
                                    if params[3] > 1 and params[3] < 2: # mu_2
                                        if params[4] > 0.01 and params[4] < 0.5: # sigma_2
                                            if params[5] > 1.5 and params[5] < 2.5: # m_TOV
                                                if self.vary_slope:
                                                    if params[6] > 0 and params[6] < 0.5: # slope
                                                        return self.pop_like(data, params)
                                                    else:
                                                        return -np.inf
                                                else:
                                                    return self.pop_like(data, params)
                return -np.inf

            if self.vary_slope:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV, self.slope]
                pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1, 0.05]
                pos = p0 + pscale*np.random.randn(14, 7)
            else:
                p0 = [self.a, self.mu_1, self.sigma_1, self.mu_2, self.sigma_2, self.m_TOV]
                pscale = [0.05, 0.1, 0.01, 0.1, 0.05, 0.1]
                pos = p0 + pscale*np.random.randn(12, 6)


        nwalkers, ndim = pos.shape
        if save_to is not None:
            backend = emcee.backends.HDFBackend(filename)
            backend.reset(nwalkers, ndim)

        sampler = emcee.EnsembleSampler(nwalkers, ndim, logpost_one, args=([samples]))
        sampler.run_mcmc(pos, steps, progress=True)
        posterior_samples = sampler.get_chain(discard = 100, flat=True)
        log_prob_samples = sampler.get_log_prob(discard = 100, flat=True)

        return posterior_samples, log_prob_samples
