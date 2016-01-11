'''
  File name           : main.py
  Author              : Fredrik Dahlin
  Date created        : 11/1/2016
  Date last modified  : 11/1/2016
  Python Version      : 3.4
'''
import priors
import utility



# Load data
data = utility.loadData()

# Define prior < AB | BB | NoZygosity | TCN | PCN >
prior = "TCN"

# Get possible states for each mutation
mutations = priors.getMutations(prior, data)