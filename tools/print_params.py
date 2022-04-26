import h5py
import numpy as np
import sys

base = sys.argv[1]
h5In = h5py.File(f'{base}.clone.h5','r')

for key in h5In['parameters'].keys():
  #print('{:20s} {:3s} {%s}'.format(key, " : ",  h5In['parameters'][key][()]))
  print('{:20s} {:3s}'.format(key, " : "), end='')
  val =  h5In['parameters'][key][()]
  if (isinstance(val, bytes)):
    print(val.decode())
  else:
    print(val)

h5In.close()


