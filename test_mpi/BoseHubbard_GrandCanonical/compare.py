import json
import sys
import h5py
import numpy as np


lattice = sys.argv[1]
with open(f"./ED_GC/{lattice}.json",'r') as edjsonf:
	edOut = json.load(edjsonf)
h5Out = h5py.File(f'./run/{lattice}.out.h5','r')
h5res = h5Out['simulation']['results']


def close_enough(refval, testval, error):
	assert error >= 0
	within_range = testval - error <= refval and refval <= testval + error
	if not within_range:
      		within_range2 = testval - error*2 <= refval and refval <= testval + error*2
      		if not within_range2:
        		print(f"ERROR: refval {refval} is not within 2sigma of testval {testval} ± {error}")
      		else:
        		print(f"WARNING: refval {refval} is not within sigma (but within 2 sigma) of testval {testval} ± {error}")
	return within_range

def veccheck(ved,vmean,verr,tolfac=1):
	assert ved.size == vmean.size and vmean.size == verr.size
	assert tolfac > 0
	b = True
	for i in range(ved.size):
		if not close_enough(ved[i], vmean[i], tolfac*verr[i]):
			b = False
	return b


print("DensDens_CorrFun")
try: #might not exist in edOut (json file)
	dd_ed = np.array(edOut['dd'], dtype=float) 
	dd = np.array(h5res['DensDens_CorrFun']['mean']['value'])
	dderr = np.array(h5res['DensDens_CorrFun']['mean']['error'])
	if veccheck(dd_ed, dd, dderr):
		print("success")
	else:
		print(dd_ed)
		print(dd, dderr)
	print()
except:
	print("skip")
	pass
print()

print("Density_Matrix")
try: #might not exist in edOut (json file)
	dm_ed = np.array(edOut['dm'], dtype=float)
	dm = np.array(h5res['Density_Matrix']['mean']['value'])
	dmerr = np.array(h5res['Density_Matrix']['mean']['error'])
	##erase first element of dm and dmerr (distance 0)
	#dm = dm[1:]
	#dmerr = dmerr[1:]
	if veccheck(dm_ed, dm, dmerr):
		print("success")
	else:
		print(dm_ed)
		print(dm, dmerr)
except:
	print("skip")
	pass
print()

print("Kinetic Energy")

Ekin_ed = float(edOut['Ekin'])
Ekin = np.array(h5res['Kinetic_Energy']['mean']['value'])
Ekinerr = np.array(h5res['Kinetic_Energy']['mean']['error'])
if close_enough(Ekin_ed, Ekin, Ekinerr):
	print("Ekin success")

print("Potential Energy")
Epot_ed = float(edOut['Epot_gc'])
Epot = np.array(h5res['Potential_Energy']['mean']['value'])
Epoterr = np.array(h5res['Potential_Energy']['mean']['error'])
if close_enough(Epot_ed, Epot, Epoterr):
	print("Epot success")

print("Total Energy")
Etot_ed = float(edOut['Etot'])
Etot = np.array(h5res['Total_Energy']['mean']['value'])
Etoterr = np.array(h5res['Total_Energy']['mean']['error'])
if close_enough(Etot_ed, Etot, Etoterr):
	print("Etot success")



print()
print("done")

