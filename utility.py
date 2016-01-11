'''
  File name           : utility.py
  Author              : Fredrik Dahlin
  Date created        : 5/1/2016
  Date last modified  : 5/1/2016
  Python Version      : 3.4
'''
import csv 
import os



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
  data = []

  # Opens the data directory and reads each individual file
  for filename in os.listdir("./Data"):
    with open("./Data/"+filename,'r') as tsv:
      reader = csv.DictReader(tsv, dialect="excel-tab")
      for line in reader:
        data.append(line)

  return data
