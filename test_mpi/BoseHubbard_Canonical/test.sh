mkdir build
cd build
cmake .. -DLATTICE="$1"	-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-11
make -j 4 VERBOSE=1
cd ..
mkdir run
cd run
mpirun -n 4 ../build/qmc_worm_mpi ../parameter_files/"$1".ini
cd ..
