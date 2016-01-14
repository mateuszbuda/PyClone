'''
  File name           : main.py
  Author              : Fredrik Dahlin
  Date created        : 11/1/2016
  Date last modified  : 11/1/2016
  Python Version      : 3.4
'''
import priors
import utility
import pyclone_binomial



# Load data
raw_data = utility.loadData()

# Define prior < AB | BB | NoZygosity | TCN | PCN >
prior = "TCN"

# Get possible states for each mutation
mutations = priors.getMutations(prior, raw_data)

#data, sample_ids, tumour_content, trace_dir, num_iters, alpha, alpha_priors
