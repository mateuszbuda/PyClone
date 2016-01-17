'''
  File name           : utility.py
  Author              : Fredrik Dahlin
  Date created        : 5/1/2016
  Date last modified  : 5/1/2016
  Python Version      : 3.4
'''
import csv 
import os
import yaml
from collections import OrderedDict, namedtuple
from mutation import Mutation, State
from math import log


def getSamplingVariantAlleleProbability():
  """ getSamplingVariantAlleleProbability | u
  Returns the probability of sampling a variant allele from a cell 
  with genotype g when b(g) is not 0 and b(g) does not equal c(g)
  
  Example: 
  getNumberOfVariantAlleles(AAB) = 1
  ------------------------------------------------------------
  """
  # c = getGenotypeCopyNumber()
  # b = getNumberOfVariantAlleles()
  # e = 0.1   # vet inte vad detta ska vara

  # if (b != 0) and (c != b):
  #   u = b/c
  # elif (b == 0):
  #   u = e
  # elif (b == c):
  #   u = 1 - e

  return "getSamplingVariantAlleleProbability"


def loadData():
  """ Load data
  mutation_id   : unique identifier for a mutation. In general 
                  specifying the gene for the mutation is a bad idea 
                  in case a gene contains multiple mutations. Usually 
                  some combination of gene name and genomic coordinates 
                  is a good choice. In this case I have used the case 
                  with the variant, the genotype of mutation in the 
                  variant case and the genomic coordinates.

  ref_counts    : the number of reads which contain the reference 
                  allele for the mutation.

  var_counts    : the number of reads which contain the variant (mutant) 
                  allele for the mutation.

  normal_cn     : the copy number of the mutant locus for the normal cells 
                  in the sample. In most cases this will be 2, with the 
                  following exceptions (there may be some others I haven't considered).
    i.  If the sample is from a male and the mutation is on a sex chromosomes (X or Y) 
        you would expect the normal cells to have copy number 1.
    ii. If the normal tissue has a germline copy number variant you would need 
        to set the copy number to the predicted value. The only way to get this 
        is to run a copy number analysis on normal tissue from the same donor.

  minor_cn      : the minor parental copy number predicted from the tumour sample.

  major_cn      : the major parental copy number predicted from the tumour sample.
  
  ------------------------------------------------------------
  Input   : None
  Output  : Array with each line as a dict
  """
  data = {}

  # Opens the data directory and reads each individual file
  for filename in os.listdir("./Data"):
    with open("./Data/"+filename,'r') as tsv:
      sample_id = filename.split('.')[0]
      data[sample_id] = []
      reader = csv.DictReader(tsv, dialect="excel-tab")
      for line in reader:
          data[sample_id].append(line)

  return data


def loadDataPyClone():
    sample_data = OrderedDict()

    for sample_id in ['SRR385938', 'SRR385939', 'SRR385940', 'SRR385941']:
        file_name = 'Data/' + sample_id + '.yaml'

        file_name = os.path.join('./', file_name)

        sample_data[sample_id] = _load_sample_data(file_name)


    sample_ids = sample_data.keys()

    common_mutations = set.intersection(*[set(x.keys()) for x in sample_data.values()])

    data = OrderedDict()

    for mutation_id in common_mutations:
        data[mutation_id] = OrderedDict()

        for sample_id in sample_ids:
            data[mutation_id][sample_id] = sample_data[sample_id][mutation_id]

    return data, sample_ids


def _load_sample_data(file_name):
    '''
    Load data from PyClone formatted input file.
    '''
    data = OrderedDict()

    fh = open(file_name)

    config = yaml.load(fh)

    fh.close()

    for mutation_dict in config['mutations']:
        mutation = load_mutation_from_dict(mutation_dict)

        data[mutation.id] = _get_pyclone_data(mutation, error_rate = 0.001)

    return data


def _get_pyclone_data(mutation, error_rate):
    a = mutation.ref_counts
    b = mutation.var_counts

    d = a + b

    cn_n = tuple([x.cn_n for x in mutation.states])
    cn_r = tuple([x.cn_r for x in mutation.states])
    cn_v = tuple([x.cn_v for x in mutation.states])

    mu_n = tuple([x.get_mu_n(error_rate) for x in mutation.states])
    mu_r = tuple([x.get_mu_r(error_rate) for x in mutation.states])
    mu_v = tuple([x.get_mu_v(error_rate) for x in mutation.states])

    prior_weights = tuple([x.prior_weight for x in mutation.states])

    log_pi = _get_log_pi(prior_weights)

    return PyCloneBinomialData(b, d, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi)


PyCloneBinomialData = namedtuple('PyCloneBinomialData',
                                 ['b', 'd', 'cn_n', 'cn_r', 'cn_v', 'mu_n', 'mu_r', 'mu_v', 'log_pi'])


def _get_log_pi(weights):
    pi = [x / sum(weights) for x in weights]

    return tuple([log(x) for x in pi])


def load_mutation_from_dict(d):
    mutation_id = d['id']

    ref_counts = int(d['ref_counts'])
    var_counts = int(d['var_counts'])

    mutation = Mutation(mutation_id, ref_counts, var_counts)

    for state_dict in d['states']:
        state = load_state_from_dict(state_dict)

        mutation.add_state(state)

    return mutation


def load_state_from_dict(d):
    g_n = d['g_n']
    g_r = d['g_r']
    g_v = d['g_v']

    prior_weight = float(d['prior_weight'])

    return State(g_n, g_r, g_v, prior_weight)


def make_directory(target_dir):
    '''
    Make target directory if it does not exist.
    '''
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)


def make_parent_directory(file_name):
    '''
    Given a file name, make the parent directory if it does not exist using make_directory.

    For example, given /some/where/foo.bar make the folder /some/where.
    '''
    file_name = os.path.abspath(file_name)

    parent_dir = os.path.dirname(file_name)

    make_directory(parent_dir)