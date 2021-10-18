#!/usr/bin/python3
'''This code calculates the jackknife error
   First check/edit : September 12, 2015
   Checked further October 27, 2015 
   '''
import sys
import numpy as np
import itertools 
from math import *
data = []; data_tot = 0. ; Data = [] ; data_jack = []

if len( sys.argv ) > 4:
    filename = sys.argv[1]
    therm_cut = int(sys.argv[2])
    blocksize  = int(sys.argv[3])
    which_column  = int(sys.argv[4])

if len( sys.argv ) <= 4:
    print("NEED 4 ARGUMENTS : FILE THERM-CUT BLOCKSIZE COLUMN_TO_PARSE")
    sys.exit()

file = open(filename, "r")
for line in itertools.islice(file, therm_cut, None):
    
    line = line.split()
    if which_column > int(np.shape(line)[0])-1:
        print ("Column to average does not exist")
        sys.exit(1)
    data_i = float(line[which_column])
    data.append(data_i)
    data_tot += data_i
    n = len(data)

n_b = int(n/blocksize)
B = 0.

for k in range(n_b):
    for w in range((k*blocksize)+1,(k*blocksize)+blocksize+1):
        B += data[w-1]
    Data.insert(k,B)
    B = 0

''' Do the jackknife estimates '''

for i in range(n_b-1):
    data_jack.append((data_tot - Data[i]) / (n - blocksize))
    data_av = data_tot / n   # Do the overall averages
    data_av = data_av
    data_jack_av = 0.; data_jack_err = 0.
for i in range(n_b-1):
    dR = data_jack[i]
    data_jack_av += dR
    data_jack_err += dR**2

data_jack_av /= n_b-1
data_jack_err /= n_b-1

data_jack_err = sqrt((n_b - 2) * abs(data_jack_err - data_jack_av**2))
print(" %8.7f "  " %6.7f"   " %6.2f" % (data_jack_av, data_jack_err, n_b))
