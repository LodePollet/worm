# usage: bash test.sh <lattice> [case]; case (default: the lattice) names
# parameter_files/<case>.ini and ED_GC/<case>.json
case="${2:-$1}"
mkdir build
cd build
cmake .. -DLATTICE="$1"	-DCMAKE_BUILD_TYPE=Release
make -j 4 VERBOSE=1
cd ..
mkdir run
cd run
mpirun -n 4 ../build/qmc_worm_mpi ../parameter_files/"$case".ini
cd ..
