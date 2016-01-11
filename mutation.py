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
