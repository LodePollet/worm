import h5py
import numpy as np
import sys

base = sys.argv[1]

# for MPI the filename below must be changed to see the statistics for rank > 0
# add a suffix  ".n" with n > 0 the rank 
h5In = h5py.File(f'{base}.clone.h5','r')


data = h5In['simulation']['realizations']['0']['clones']['0']['checkpoint']['update_statistics'][()]

names = ["INSERT WORM ", "MOVE WORM   ", "INSERT KINK ", "DELETE KINK ", "GLUE WORM   "]
print("# UPDATE STATISTICS")
print("# col 1 : all updates")
print("# col 2 : impossible updates")
print("# col 3 : rejected updates")
print("# col 4 : accepted updates")
print("# col 5 : acceptance factor with respect to Metropolis ratio only")
print("# col 6 : acceptance factor with respect to all attempts.\n")

s = 0
for i in np.arange(data.shape[1]):
  s += data[-1][i]
  print('# ', names[i], '%15d %15d %15d %15d %5.3f %5.3f' % (data[-1][i], data[0][i], data[1][i], data[2][i], data[2][i] / (data[1][i] + data[2][i]), data[2][i] / data[-1][i] ))

print()
print()
print("# Total number of steps : ", s , " or 10^", np.log10(s))
print()
print()

h5In.close()


