import h5py
import numpy as np
import sys

# hack to overwrite parameters inside the h5 files
# useful for forcing a reset of statistics and sweeps in case the user decides that the thermalization is insufficient but the user wants to keep the configuration
# possible extensions : simulated annealing approaches when changing system parameters (not implemented in the code)

new_value = 0

base = sys.argv[1]
h5InOut = h5py.File(f'{base}.clone.h5','r+')

h5InOut['parameters']['reset_statistics'][()] = new_value
h5InOut.close()


