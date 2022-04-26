mkdir build
cd build
cmake .. -DLATTICE="$1"	-DCMAKE_BUILD_TYPE=Release
make -j 4 VERBOSE=1
cd ..
mkdir run
cd run
../build/qmc_worm ../parameter_files/"$1".ini
cd ..
