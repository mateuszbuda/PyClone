__author__ = 'mateusz'

from collections import OrderedDict, namedtuple
from trace import DiskTrace
from DirichletProcess.measures import BetaBaseMeasure, MultiSampleBaseMeasure
from DirichletProcess.densities import PyCloneBinomialDensity, MultiSampleDensity
from DirichletProcess.samplers.atom import BaseMeasureAtomSampler, MultiSampleAtomSampler
from DirichletProcess.samplers.partition import AuxillaryParameterPartitionSampler
from DirichletProcess.samplers.dp import DirichletProcessSampler

PyCloneBinomialParameter = namedtuple('PyCloneBinomialParameter', 'tumour_content')

def run_pyclone_binomial_analysis(data, sample_ids, tumour_content, trace_dir, num_iters, alpha, alpha_priors):

    sample_atom_samplers = OrderedDict()

    sample_base_measures = OrderedDict()

    sample_cluster_densities = OrderedDict()

    alpha = 1
    beta = 1

    for sample_id in sample_ids:
        sample_base_measures[sample_id] = BetaBaseMeasure(alpha, beta)

        sample_cluster_densities[sample_id] = PyCloneBinomialDensity(PyCloneBinomialParameter(tumour_content[sample_id]))

        sample_atom_samplers[sample_id] = BaseMeasureAtomSampler(sample_base_measures[sample_id],
                                                                 sample_cluster_densities[sample_id])

    base_measure = MultiSampleBaseMeasure(sample_base_measures)

    cluster_density = MultiSampleDensity(sample_cluster_densities)

    atom_sampler = MultiSampleAtomSampler(base_measure, cluster_density, sample_atom_samplers)

    partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)

    sampler = DirichletProcessSampler(atom_sampler, partition_sampler, alpha, alpha_priors)

    trace = DiskTrace(trace_dir, sample_ids, data.keys(), {'cellular_frequencies' : 'x'})

    trace.open()

    sampler.sample(data.values(), trace, num_iters)

    trace.close()
