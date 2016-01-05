'''
    File name           : priors.py
    Author              : Fredrik Dahlin
    Date created        : 5/1/2016
    Date last modified  : 5/1/2016
    Python Version      : 3.4
'''


def getPrior(prior):
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
    return NoZygosityPrior()
  elif prior == "TCN":
    return TCNPrior()
  elif prior == "PCN":
    return PCNPrior()
  else:
    raise ValueError(prior,"is not a implemented prior. Available priors are AB, BB, NoZygosity, TCN and PCN.")


def ABPrior():
  """ ABPrior
  Assumes that gn = AA, gr = AA and gv = AB. 
  Each mutation is assumed to be diploid and heterozygous.
  ------------------------------------------------------------
  """
  gr = "AA"
  gv = "AB"
  gn = "AA"

  return "ABPrior"


def BBPrior():
  """ BBPrior
  Assumes that gn = AA, gr = AA and gv = BB. 
  Each mutation is assumed to be diploid and homozygous.
  ------------------------------------------------------------
  """
  gr = "AA"
  gv = "BB"
  gn = "AA"

  return "BBPrior"


def NoZygosityPrior():
  """ NoZygosityPrior
  Assumes that gn = AA, gr = AA, c(gv) = C and b(gv) = 1. 
  The genotype of the variant population has the predicted 
  copy number with exactly one mutant allele.
  ------------------------------------------------------------
  """
  gr = "AA"
  gn = "AA"

  return "NoZygosityPrior"


def TCNPrior():
  """ TCNPrior
  Assumes that gn = AA, c(gv) = C and b(gv) is an element of {1,...,C}. 
  The genotype of the variant population has the predicted 
  copy number and at least on variant allele. 

  Also assumes with equal probability that gr = AA or 
  c(gr) = C and b(gr) = 0.
  The genotype of the variant population at the locus has the 
  predicted total copy number and we consider the possibility 
  that any number of copies (> 0) of the locus contains the
  mutant allele.  
  ------------------------------------------------------------
  """
  gr = "AA"
  gn = "AA"

  return "TCNPrior"


def PCNPrior():
  """ PCNPrior
  Assumes that gn = AA, c(gv) = C and b(gv) is an element in {1, C1, C2}. 
  The genotype of the variant population has the genotype with
  the predicted copy number and one variant allele, or as many
  variant alleles as one of the parental copy numbers. 

  if b(gv) is an element in {C1, C2}
    gr = gn 
    The mutation occurs before copy number events.
  if b(gv) = 1
    c(gr) = C
    b(gr) = 0
    The mutation occurs after the copy number events. 
  ------------------------------------------------------------
  """
  gn = "AA"

  return "PCNPrior"


prior = input("Select prior: ")
print(getPrior(prior))