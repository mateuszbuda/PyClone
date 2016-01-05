'''
    File name           : utility.py
    Author              : Fredrik Dahlin
    Date created        : 5/1/2016
    Date last modified  : 5/1/2016
    Python Version      : 3.4
'''


def getGenotypeCopyNumber():
  """ getGenotypeCopyNumber | c
  Returns the copy number of the genotype.
  
  Example: 
  getGenotypeCopyNumber(AAB) = 3
  ------------------------------------------------------------
  """

  return "getGenotypeCopyNumber"


def getNumberOfVariantAlleles():
  """ getNumberOfVariantAlleles | b
  Returns the number of variant alleles in the genotype. 
  
  Example: 
  getNumberOfVariantAlleles(AAB) = 1
  ------------------------------------------------------------
  """

  return "getNumberOfVariantAlleles"


def getSamplingVariantAlleleProbability():
  """ getSamplingVariantAlleleProbability | u
  Returns the probability of sampling a variant allele from a cell 
  with genotype g when b(g) is not 0 and b(g) does not equal c(g)
  
  Example: 
  getNumberOfVariantAlleles(AAB) = 1
  ------------------------------------------------------------
  """
  c = getGenotypeCopyNumber()
  b = getNumberOfVariantAlleles()
  e = 0.1   # vet inte vad detta ska vara

  if (b != 0) and (c != b):
    u = b/c
  elif (b == 0):
    u = e
  elif (b == c):
    u = 1 - e

  return "getSamplingVariantAlleleProbability"