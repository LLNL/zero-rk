import numpy as np
import struct
import matplotlib.pyplot as plt

datafile1= "datap100000.0phi1.0Tu400.0"
gridfile1= "gridp100000.0phi1.0Tu400.0"

### Read grid files
grid1 = np.genfromtxt(gridfile1, delimiter='')

# Read data file
data1 = open(datafile1, "rb").read()

(num_points1, num_vars, time) = struct.unpack_from("iid", data1)
print('Number of points and variables:',num_points1, num_vars)

index = 0
# Read variable names
name = []
for j in range(num_vars):
    name.append(struct.unpack_from("64c",data1,offset=4+4+8+index))
    index = index + 64

# Read variable values
Y1=np.zeros((num_vars+1,num_points1))
for j in range(num_vars+1):
    for i in range(num_points1):
        Y1[j,i] = struct.unpack_from("d",data1,offset=4+4+8+index)[0]
        index = index + 8

# Write to ascii file
fh = open("ascii_data","w")
# Write data
for i in range(num_points1):
    fh.write(str(grid1[i]) + '\t')
    for j in range(num_vars+1):
        fh.write(str(Y1[j,i]) + '\t')
    fh.write('\n')
fh.close()
