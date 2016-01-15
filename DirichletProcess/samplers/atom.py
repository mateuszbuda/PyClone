__author__ = 'mateusz'

from math import log
from random import uniform
from collections import OrderedDict
from ..partition import PartitionCell

class AtomSampler(object):
    '''
    Base class for samplers to update the cell values in the partition (atoms of DP).
    '''
    def __init__(self, base_measure, cluster_density):
        '''
        Args:
            base_measure : (BaseMeasure) Base measure for DP process.

            cluster_density : (Density) Cluster density for DP process.
        '''
        self.base_measure = base_measure

        self.cluster_density = cluster_density

    def sample(self, data, partition):
        '''
        Sample a new value for atoms in the partition. The partition passed in will be updated in place.

        Args:
            data : (list) List of data points appropriate for cluster_density.

            partition : (Partition) Partition of dp.
        '''
        for cell in partition.cells:
            cell.value = self.sample_atom(data, cell)

    def sample_atom(self, data, cell):
        '''
        Sample a new value for the atom associated with the cell. Returns a suitable value for the cell.
        '''
        raise NotImplemented


class MetropolisHastingsAtomSampler(AtomSampler):
    '''
    Update the atom values using a Metropolis-Hastings steps with a user specified proposal function which takes
    the previous cell value as an argument.
    '''
    def __init__(self, base_measure, cluster_density, proposal_func):
        AtomSampler.__init__(self, base_measure, cluster_density)

        self.proposal_func = proposal_func

    def sample_atom(self, data, cell):
        old_param = cell.value
        new_param = self.proposal_func.random(old_param)

        old_ll = self.base_measure.log_p(old_param)
        new_ll = self.base_measure.log_p(new_param)

        for j in cell.items:
            old_ll += self.cluster_density.log_p(data[j], old_param)
            new_ll += self.cluster_density.log_p(data[j], new_param)

        forward_log_ratio = new_ll - self.proposal_func.log_p(new_param, old_param)
        reverse_log_ratio = old_ll - self.proposal_func.log_p(old_param, new_param)

        log_ratio = forward_log_ratio - reverse_log_ratio

        u = uniform(0, 1)

        if log_ratio >= log(u):
            return new_param
        else:
            return old_param

class BaseMeasureAtomSampler(MetropolisHastingsAtomSampler):
    '''
    Update the atom values using a Metropolis-Hastings steps with the base measure as a proposal density.
    '''
    def __init__(self, base_measure, cluster_density):
        proposal_func = BaseMeasureProposalFunction(base_measure)

        MetropolisHastingsAtomSampler.__init__(self, base_measure, cluster_density, proposal_func)

class BaseMeasureProposalFunction(object):
    def __init__(self, base_measure):
        self.base_measure = base_measure

    def log_p(self, data, params):
        return self.base_measure.log_p(data)

    def random(self, params):
        return self.base_measure.random()


class MultiSampleAtomSampler(AtomSampler):
    def __init__(self, base_measure, cluster_density, atom_samplers):
        AtomSampler.__init__(self, base_measure, cluster_density)

        self.atom_samplers = atom_samplers

    def sample_atom(self, data, cell):
        new_atom = OrderedDict()

        for sample_id in self.atom_samplers:
            sample_data = [x[sample_id] for x in data]

            sample_cell = PartitionCell(cell.value[sample_id])

            sample_cell._items = cell._items

            new_atom[sample_id] = self.atom_samplers[sample_id].sample_atom(sample_data, sample_cell)

        return new_atom