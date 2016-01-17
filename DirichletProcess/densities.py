__author__ = 'mateusz'

from math import log, lgamma
from collections import OrderedDict, namedtuple
from utils import log_sum_exp

def log_beta_pdf(x, a, b):
    if x == 0 or x == 1:
        return float('-inf')

    return -log_beta(a, b) + (a - 1) * log(x) + (b - 1) * log(1 - x)


def log_beta(a, b):
    if a <= 0 or b <= 0:
        return float('-inf')

    return lgamma(a) + lgamma(b) - lgamma(a + b)

def log_binomial_pdf(x, n, p):
    if p == 0:
        if x == 0:
            return 0
        else:
            return float('-inf')

    if p == 1:
        if x == n:
            return 0
        else:
            return float('-inf')

    return log_binomial_coefficient(n, x) + x * log(p) + (n - x) * log(1 - p)


def log_binomial_coefficient(n, x):
    return log_factorial(n) - log_factorial(x) - log_factorial(n - x)


def log_factorial(n):
    return lgamma(n + 1)


class Density(object):
    def __init__(self, params=None):
        self.params = params

        self.cache = OrderedDict()

        self.max_cache_size = 10000

    def log_p(self, data, params):
        '''
        Args:
            data : (nametuple) Data for density.

            params : (nametuple) Parameters in density.

        Kwargs:
            global_params: (namedtuple) Parameters which are shared across all atoms. If this is None it will use the
                                        current value.
        '''
        key = (data, params, self.params)

        if key not in self.cache:
            self.cache[key] = self._log_p(data, params)

            if len(self.cache) > self.max_cache_size:
                self.cache.popitem(last=False)

        return self.cache[key]

    def _log_p(self, data, params):
        raise NotImplemented

class PyCloneBinomialDensity(Density):
    def _log_p(self, data, params):
        ll = []

        for cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi  in zip(data.cn_n, data.cn_r, data.cn_v, data.mu_n, data.mu_r, data.mu_v, data.log_pi):
            temp = log_pi + self._log_binomial_likelihood(data.b,
                                                          data.d,
                                                          cn_n,
                                                          cn_r,
                                                          cn_v,
                                                          mu_n,
                                                          mu_r,
                                                          mu_v,
                                                          params.x)

            ll.append(temp)

        return log_sum_exp(ll)

    def _log_binomial_likelihood(self, b, d, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, cellular_frequency):
        f = cellular_frequency
        t = self.params.tumour_content

        p_n = (1 - t) * cn_n
        p_r = t * (1 - f) * cn_r
        p_v = t * f * cn_v

        norm_const = p_n + p_r + p_v

        p_n = p_n / norm_const
        p_r = p_r / norm_const
        p_v = p_v / norm_const

        mu = p_n * mu_n + p_r * mu_r + p_v * mu_v

        return log_binomial_pdf(b, d, mu)

class MultiSampleDensity(Density):
    '''
    Wraps a collection of univariate densities.
    '''
    def __init__(self, cluster_densities, shared_params=False):
        '''
        Args:
            cluster_densities: (dict) A collection of Density objects for each sample.
        '''
        self.cluster_densities = cluster_densities

        self.shared_params = shared_params

    @property
    def params(self):
        if self.shared_params:
            for cluster_id in self.cluster_densities:
                return self.cluster_densities[cluster_id].params

        else:
            params = OrderedDict()

            for cluster_id in self.cluster_densities:
                params[cluster_id] = self.cluster_densities[cluster_id].params

            return params

    @params.setter
    def params(self, value):
        if self.shared_params:
            for cluster_id in self.cluster_densities:
                self.cluster_densities[cluster_id].params = value

        elif isinstance(value, namedtuple):
            for cluster_id in self.cluster_densities:
                self.cluster_densities[cluster_id].params = value[cluster_id]

        else:
            raise Exception('Cannot set object type {0} as a density parameter'.format(type(value)))


    def log_p(self, data, params):
        log_p = 0

        for sample_id in self.cluster_densities:
            density = self.cluster_densities[sample_id]

            log_p += density.log_p(data[sample_id], params[sample_id])

        return log_p