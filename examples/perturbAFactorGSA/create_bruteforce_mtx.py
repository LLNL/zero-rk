import numpy as np

# -----------------------------------------------------------------------------
# user defined inputs:
matrix_file = 'hydrogen_brute.mtx'
num_reactions = 21  # number of columns, must match mechanism count

# creates a GSA matrix file where each reaction is perturbed by afactor_mult
# individually, it is expected to produce the same IDTs as the brute force
# A-Factor utility with a single perturbation (i.e. multiplication only).
# The GSA matrix is then a square matrix of all ones except for the diagonal
# which has a value of afactor_mult.
afactor_mult = 2.0

# -----------------------------------------------------------------------------

# open the matrix file
file_ptr = open(matrix_file,'w')

# write the matrix size to file (# rows, # columns)
file_ptr.write('{0:6d} {1:6d}\n'.format(num_reactions, num_reactions))

# create the matrix of A-Factor multiplicative perturbations
for j in range(num_reactions):
  for k in range(num_reactions):
    val = 1.0
    if (j == k):
      val = afactor_mult
    
    file_ptr.write('{0:24.16e}\n'.format(val))
