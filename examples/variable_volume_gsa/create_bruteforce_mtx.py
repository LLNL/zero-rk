
# -----------------------------------------------------------------------------
# user defined inputs:
matrix_file = 'hydrogen_bruteforce.mtx'
num_reactions = 21              # number of columns, must match mechanism count
num_samples = num_reactions     # number of rows
afactor_mult = 2.0

# open the matrix file
file_ptr = open(matrix_file,'w')

# write the matrix size to file (# rows, # columns)
file_ptr.write('{0:6d} {1:6d}\n'.format(num_samples, num_reactions))

# create the matrix of A-Factor multiplicative perturbations
for j in range(num_samples):
  multipliers = [1.0]*num_reactions
  multipliers[j] = afactor_mult
  for m in multipliers:
    file_ptr.write('{0:24.16e}\n'.format(m))

