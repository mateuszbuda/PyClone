'''
  File name           : priors.py
  Author              : Fredrik Dahlin
  Date created        : 5/1/2016
  Date last modified  : 11/1/2016
  Python Version      : 3.4
'''
import mutation



def getMutations(prior, data):
  ''' getMutations
  Creates a mutation object for each mutation and create it's corresponding 
  states. Save all mutations in an array and return it to main. 
  
  ------------------------------------------------------------------------------
  Input   : prior name 
          : data array
  Output  : mutation array
  '''  
  mutations = []

  for i in data:
    # Creates a mutation object
    m = mutation.Mutation(i["mutation_id"], i["ref_counts"], i["var_counts"])

    # Creates states corresponding to that mutation
    states = getPrior(prior, int(i["normal_cn"]), int(i["minor_cn"]), int(i["major_cn"]))

    # Save the states in the mutation object
    m.setStates(states)

    mutations.append(m)

  return mutations


def getPrior(prior, normal_cn, minor_cn, major_cn):
  ''' getPrior
  Flow control, returns one of theses priors (AB | BB | NoZygosity | TCN | PCN )
  
  ------------------------------------------------------------------------------
  Input   : prior name 
  Output  : prior array
  '''
  if prior == "AB":
    return ABPrior()
  elif prior == "BB":
    return BBPrior()
  elif prior == "NoZygosity":
    return NoZygosityPrior(normal_cn, minor_cn, major_cn)
  elif prior == "TCN":
    return TCNPrior(normal_cn, minor_cn, major_cn)
  elif prior == "PCN":
    return PCNPrior(normal_cn, minor_cn, major_cn)
  else:
    raise ValueError(prior,"is not a implemented prior. Available priors are AB, BB, NoZygosity, TCN and PCN.")


def ABPrior():
  """ ABPrior
  Assumes that g_n = AA, g_r = AA and g_v = AB. 
  Each mutation is assumed to be diploid and heterozygous.

  ------------------------------------------------------------
  Input   : None
  Output  : ("AA", "AA", "AB") - (g_n, g_r, g_v)
          : Prior weight
  """

  return [("AA", "AA", "AB", 1)]


def BBPrior():
  """ BBPrior
  Assumes that g_n = AA, g_r = AA and g_v = BB. 
  Each mutation is assumed to be diploid and homozygous.

  ------------------------------------------------------------
  Input   : None
  Output  : ("AA", "AA", "BB") - (g_n, g_r, g_v)  
          : Prior weight  
  """

  return [("AA", "AA", "BB", 1)]


def NoZygosityPrior(normal_cn, minor_cn, major_cn):
  """ NoZygosityPrior
  Assumes that g_n = AA, g_r = AA, c(g_v) = total_cn and b(g_v) = 1. 

  Assumes the genotype of the variant population has a copy number 
  which matches the predicted total value. It assumes the genotype has only 
  one B allele. Thus if the predicted total copy number was 4, 
  the variant genotype would be set to AAAB.

  ------------------------------------------------------------
  Input   : normal_cn - copy number of the mutant locus for the normal cells in the sample. 
            minor_cn  - minor parental copy number predicted from the tumour sample.
            major_cn  - major parental copy number predicted from the tumour sample.
  Output  : ("AA", "AA", g_v) - (g_n, g_r, g_v)  
          : Prior weight  
  """
  total_cn = minor_cn + major_cn
  g_n = "A" * normal_cn
  g_v = "A" * (total_cn - 1) + "B"

  return [(g_n, "AA", g_v, 1)]


def TCNPrior(normal_cn, minor_cn, major_cn):
  """ TCNPrior
  Assumes that g_n = AA, c(g_v) = total_cn and b(g_v) is an element of {1,...,total_cn}. 
  The genotype of the variant population has the predicted copy number and at least one 
  variant allele. 

  This means that g_v has the length total_cn and g_v consists of 0 or (total_cn - 1)
  'A' and 1 to total_cn 'B'. Thus if the total copy number was 4, the possible genotypes 
  for g_v would be AAAB, AABB, ABBB, BBBB; 

  Also assumes with equal probability that g_r = AA or c(g_r) = total_cn and b(g_r) = 0.
  The genotype of the variant population at the locus has the predicted total copy number 
  and we consider the possibility that any number of copies (> 0) of the locus contains the
  mutant allele.  

  This means that g_r are either 'AA' or 'A' * total_cn with equal probability. 

  ------------------------------------------------------------
  Input   : normal_cn - copy number of the mutant locus for the normal cells in the sample. 
            minor_cn  - minor parental copy number predicted from the tumour sample.
            major_cn  - major parental copy number predicted from the tumour sample.
  Output  : ("AA", "AA", g_v) - (g_n, g_r, g_v)    
          : Prior weight  
  """
  total_cn = minor_cn + major_cn
  g_n = ["A" * normal_cn] * 2 * total_cn
  g_r = ["AA", "A" * total_cn] * total_cn
  g_v = []
  prior_weight = [1] * 2 * total_cn
  
  for i in range(1, total_cn + 1):
    genotype = "A" * (total_cn - i) + "B" * i
    g_v.append(genotype)
    g_v.append(genotype)

  states = zip(g_n, g_r, g_v, prior_weight)

  return states


def PCNPrior(normal_cn, minor_cn, major_cn):
  """ PCNPrior
  Assumes that g_n = AA, c(g_v) = total_cn and b(g_v) is an element in {1, C1, C2}. 
  The genotype of the variant population has the genotype with the predicted copy 
  number and one variant allele, or as many variant alleles as one of the parental 
  copy numbers. 

  if b(g_v) is an element in {C1, C2}
    g_r = g_n 
    The mutation occurs before copy number events.

  if b(g_v) = 1
    c(g_r) = total_cn
    b(g_r) = 0
    The mutation occurs after the copy number events. 

  Assumes the genotype of the variant population has a copy number which matches 
  the predicted total value (the sum of minor_cn and major_cn in the input .tsv file). 
  It further assumes only genotypes with as many B alleles as the major or 
  minor copy number or one B allele are allowed. If the number of B allels in 
  the genotype of the variant population matches the major or minor copy number 
  then the genotype of the reference population is assumed to match the normal population. 
  If the genotype of the variant population has only one B allele the genotype of the 
  reference population is assumed to have the predicted total copy number but be all A alleles.

  ------------------------------------------------------------
  Input   : normal_cn - This is the copy number of the mutant locus for the normal cells in the sample. 
            minor_cn  - This is the minor parental copy number predicted from the tumour sample.
            major_cn  - This is the major parental copy number predicted from the tumour sample.
  Output  : ("AA", < g_n | g_r >, b(g_v)) - (g_n, g_r, g_v)  
          : Prior weight  
  """
  total_cn = minor_cn + major_cn
  g_n = "A" * normal_cn


  if normal_cn == total_cn:
    if minor_cn == 0:
      g_v = ["A" * minor_cn + "B" * major_cn, "A" * (total_cn - 1) + "B"]
      g_r = [g_n, g_n]
    else:
      g_v = ["A" * major_cn + "B" * minor_cn, "A" * minor_cn + "B" * major_cn, "A" * (total_cn - 1) + "B"]
      g_r = [g_n, g_n, g_n]
        
  else:
    if minor_cn == 0:
      g_v = ["A" * minor_cn + "B" * major_cn, "A" * (total_cn - 1) + "B", "A" * (total_cn - 1) + "B"]
      g_r = [g_n, g_n, "A" * total_cn]
    else:
      g_v = ["A" * major_cn + "B" * minor_cn, "A" * minor_cn + "B" * major_cn, "A" * (total_cn - 1) + "B", "A" * (total_cn - 1) + "B"]
      g_r = [g_n, g_n, g_n, "A" * total_cn]


  states = []
  for s in g_v:
    states.append((g_n, g_r[s], g_v[s], 1))

  return states