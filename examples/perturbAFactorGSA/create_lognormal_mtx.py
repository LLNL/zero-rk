import numpy as np

# -----------------------------------------------------------------------------
# user defined inputs:
matrix_file = 'hydrogen_1000.mtx'
num_reactions = 21  # number of columns, must match mechanism count
num_samples = 1000   # number of rows
seed = 0

# The same lognormal distribution is used to sample all the A-Factor
# Multipliers. This section should be modified to generate different
# distributions for each reaction.
std_dev = 0.8 # standard deviation of the underlying normal distribution
              # used to generate the lognormal distribution 
afactor_std_dev = std_dev*np.ones(num_reactions)
# -----------------------------------------------------------------------------

# if seed is not set, then different samples as generated from run to run
np.random.seed(seed)

# open the matrix file
file_ptr = open(matrix_file,'w')

# write the matrix size to file (# rows, # columns)
file_ptr.write('{0:6d} {1:6d}\n'.format(num_samples, num_reactions))

# create the matrix of A-Factor multiplicative perturbations
for j in range(num_samples):
  for k in range(len(afactor_std_dev)):

    # sample A-Factor multiplier from a lognormal distribution with mean 0.0
    # and a standard deviation for log(afactor_mult) for reaction k.

    afactor_mult = np.random.lognormal(0.0, # mean of log(afactor_mult)
                                       afactor_std_dev[k])
    
    file_ptr.write('{0:24.16e}\n'.format(afactor_mult))
