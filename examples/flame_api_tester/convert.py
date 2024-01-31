
import sys
import numpy as np
import struct

#input binary file (from unsteady solver)
datafile=sys.argv[1]
#output text file (used as input to flame_api_tester)
asciifile=sys.argv[2]

#grid length in meters
#length =  0.015
length = float(sys.argv[3])

# Read data file
data = open(datafile, "rb").read()

(num_points, num_vars, time) = struct.unpack_from("iid", data)
#N.B. we skip the first point
print('Number of points and variables:',num_points-1, num_vars)

grid = np.linspace(0,length,num_points)

index = 0
# Read variable names
names = []
for j in range(num_vars):
    n = struct.unpack_from("64s",data,offset=4+4+8+index)[0].split(b'\0')[0].decode()
    names.append(n)
    index = index + 64

# Read variable values
Y=np.zeros((num_vars,num_points))
for j in range(num_vars):
    for i in range(num_points):
        Y[j,i] = struct.unpack_from("d",data,offset=4+4+8+index)[0]
        index = index + 8

# Write to ascii file
fh = open(asciifile,"w")
# Write header
fh.write('# x(m)\t')
for j in range(num_vars):
    if(names[j] == "RelativeVolume"):
        continue
    fh.write(names[j].strip('MassFraction_')+'\t')
fh.write('\n')
# Write data
for i in range(num_points):
    fh.write(str(grid[i]) + '\t')
    for j in range(num_vars):
        if(names[j] == "RelativeVolume"):
            continue
        fh.write(str(Y[j,i]) + '\t')
    fh.write('\n')
fh.close()

