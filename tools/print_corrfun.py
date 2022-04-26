import sys
import h5py
import numpy as np


lattice = sys.argv[1]
h5Out = h5py.File(f'./{lattice}.out.h5','r')
h5res = h5Out['simulation']['results']

print("Density_Distribution")
nval = np.array(h5res['Density_Distribution']['mean']['value'])
nerr = np.array(h5res['Density_Distribution']['mean']['error'])
print(nval)
print(nerr)
print()


print("DensDens_CorrFun")
dd = np.array(h5res['DensDens_CorrFun']['mean']['value'])
dderr = np.array(h5res['DensDens_CorrFun']['mean']['error'])
print(dd)
print(dderr)
print()

print("Density_Matrix")
dm = np.array(h5res['Density_Matrix']['mean']['value'])
dmerr = np.array(h5res['Density_Matrix']['mean']['error'])
##erase first element of dm and dmerr (distance 0)
#dm = dm[1:]
#dmerr = dmerr[1:]
print(dm)
print(dmerr)
print()

#print("Green's function")
#gft = np.array(h5res['Greenfun_p0_tau']['mean']['value'])
#gfterr = np.array(h5res['Greenfun_p0_tau']['mean']['error'])
#print(gft)
#print(gfterr)

