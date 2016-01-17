__author__ = 'mateusz'

from collections import namedtuple, OrderedDict
from densities import log_beta_pdf
from random import betavariate

BetaData = namedtuple('BetaData', 'x')
BetaParameter = namedtuple('BetaPriorData', ['a', 'b'])

class BaseMeasure(object):
    '''
    Base class for base measures.
    '''
    def log_p(self, data):
        '''
        Return the log probability of the density.

        Args:
            data : An data object of the same type as returned by self.random()
        '''
        raise NotImplemented


    def random(self):
        '''
        Return a random sample from the base measure.
        '''
        raise NotImplemented


class BetaBaseMeasure:
    def __init__(self, a, b):
        self.params = BetaParameter(a, b)


    def log_p(self, data):
        return log_beta_pdf(data.x, self.params.a, self.params.b)


    def random(self):
        x = betavariate(self.params.a, self.params.b)

        return BetaData(x)


class MultiSampleBaseMeasure(BaseMeasure):
    def __init__(self, base_measures):
        '''
        Args:
            base_measures: (dict) Mapping of sample IDs to base measures.
        '''
        self.base_measures = base_measures

    def log_p(self, data):
        log_p = 0

        for sample_id in self.base_measures:
            log_p += self.base_measures[sample_id].log_p(data[sample_id])

        return log_p

    def random(self):
        random_sample = OrderedDict()

        for sample_id in self.base_measures:
            random_sample[sample_id] = self.base_measures[sample_id].random()

        return random_sample