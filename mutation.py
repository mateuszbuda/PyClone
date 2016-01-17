'''
  File name           : mutation.py
  Author              : Fredrik Dahlin
  Date created        : 11/1/2016
  Date last modified  : 11/1/2016
  Python Version      : 3.4
'''



class Mutation():
  """ Mutation
  Helper class that keeps track of each mutation and it's states
  """
  def __init__(self, mutation_id, ref_counts, var_counts):
    self.id = mutation_id
    self.ref_counts = ref_counts
    self.var_counts = var_counts
    self.states = []


  def setStates(self, states):
    self.states = states


  def getStates(self):
    return self.states

  def add_state(self, state):
        self.states.append(state)

  def to_dict(self):
        return {
                'id' : self.id,
                'ref_counts' : self.ref_counts,
                'var_counts' : self.var_counts,
                'states' : [x.to_dict() for x in self.states]
                }

class State(object):
    def __init__(self, g_n, g_r, g_v, prior_weight):
        self.g_n = g_n

        self.g_r = g_r

        self.g_v = g_v

        self.prior_weight = prior_weight

    @property
    def cn_n(self):
        return len(self.g_n)

    @property
    def cn_r(self):
        return len(self.g_r)

    @property
    def cn_v(self):
        return len(self.g_v)

    def get_mu_n(self, error_rate):
        return self._get_variant_allele_probability(self.g_n, error_rate)

    def get_mu_r(self, error_rate):
        return self._get_variant_allele_probability(self.g_r, error_rate)

    def get_mu_v(self, error_rate):
        return self._get_variant_allele_probability(self.g_v, error_rate)

    def to_dict(self):
        return {'g_n' : self.g_n, 'g_r' : self.g_r, 'g_v' : self.g_v, 'prior_weight' : self.prior_weight}

    def _get_copy_number(self, genotype):
        if genotype is None:
            return 0
        else:
            return len(genotype)

    def _get_variant_allele_probability(self, genotype, error_rate):
        if genotype is None:
            return error_rate

        num_ref_alleles = genotype.count("A")
        num_var_alleles = genotype.count("B")

        cn = len(genotype)

        if cn != num_ref_alleles + num_var_alleles:
            raise Exception("{0} is not a valid genotype. Only A or B are allowed as alleles.")

        if num_ref_alleles == 0:
            return 1 - error_rate
        elif num_var_alleles == 0:
            return error_rate
        else:
            return num_var_alleles / cn
