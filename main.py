'''
  File name           : main.py
  Author              : Fredrik Dahlin
  Date created        : 11/1/2016
  Date last modified  : 11/1/2016
  Python Version      : 3.4
'''
import priors
import utility
# import pyclone_binomial



# Load data
raw_data = utility.loadData()


# Define prior < AB | BB | NoZygosity | TCN | PCN >
prior = "TCN"


# Get possible states for each mutation
data = {}
for sample in raw_data:
  mutations = priors.getMutations(prior, raw_data[sample])
  data[sample] = mutations
  

#data, sample_ids, tumour_content, trace_dir, num_iters, alpha, alpha_priors

error_rate = 0.001

tumour_content = {}
for id in data:
  tumour_content[id] = 1.0


trace_dir = 'trace'

num_iters = 10000

alpha = 1

alpha_priors = {
	'shape': 1.0,
	'rate': 0.001
}

#pyclone_binomial.run_pyclone_binomial_analysis(mutations, sample_ids, tumour_content, trace_dir, num_iters, alpha, alpha_priors)

