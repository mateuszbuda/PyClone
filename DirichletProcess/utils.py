__author__ = 'mateusz'

from math import isinf, exp, log

def log_sum_exp(log_X):
    '''
    Given a list of values in log space, log_X. Compute exp(log_X[0] + log_X[1] + ... log_X[n])

    Numerically safer than naive method.
    '''
    max_exp = max(log_X)

    if isinf(max_exp):
        return max_exp

    total = 0

    for x in log_X:
        total += exp(x - max_exp)

    return log(total) + max_exp
