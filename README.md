Efficient and scalable Path Integral Monte Carlo Simulations with worm-type updates for Bose-Hubbard and XXZ models
===================================================================================================================

This is an open source implementation of the worm algorithm for sign-positive models 
as explained in our preprint [arXiv:2204.XXX](https://arxiv.org/abs/12204.XXX).
Please refer to the main text for an in-depth discussion of the method,
particularly the algorithm and the updates.
This README focuses on the technical details in order to run the codes and
reproduce the data shown in the figures of the main text.

Structure
---------

There are 5 directories in this repo
  * paper : contains the pdf of the accompanying paper plus the data used to generate the figures in the paper.
  * src : the source code of the worm algorithm
  * parameter_files : examples of parameter files for Bose-Hubbard and spin-XXZ models. 
  * test_mpi : parameter files to test the code against ground state Lanczos results for small system sizes.
  * tools : simple helper python scripts to extract information from the output hdf5 files
    

Requirements
------------

These codes are based on the [ALPSCore](https://github.com/ALPSCore/ALPSCore)
library. Refer to [their website](http://alpscore.org/) for installation
instructions. At the time of this writing, ALPSCore imposed the following system
requirements:

  * C++ compiler, along with suitable MPI wrappers (e.g. OpenMPI)
  * CMake build system (version 3.1 or later),
  * Boost headers and `program_options` library (version 1.56 or later),
  * HDF5 library version 1.8.x (version 1.10 does _not_ work, see below).
  
Beyond these, our codes require

  * a C++11-capable compiler. We have only tested our code when ALPSCore is installed with modern compilers that have C++11 features turned on by default.





Building and Installation
-------------------------

### Building and installing ALPSCore

Detailed instructions on
[how to build ALPSCore](https://github.com/ALPSCore/ALPSCore/wiki/Installation)
can be fournd in the project's wiki. The procedure revolves around the following:

    $ cd alpscore
    $ mkdir build.tmp && cd build.tmp
    $ cmake ..
    $ make -jN
    $ make test
    $ make install
    
Replace `N` with the number of processors you want to use to build, e.g. `-j8`.
You may want to specify additional flags to `cmake`:

  * `-DCMAKE_INSTALL_PREFIX=$HOME/.local`, or another custom install location.
    This is required if you don't have permission to write at the default
    install prefix (`/usr/local`). Mind that ALPSCore installs a CMake script
    that has to be picked up by CMake when building our codes. Thus, any
    non-standard install location needs to be matched by a
    `-DCMAKE_PREFIX_PATH=<...>` flag when configuring the client code.
  * If a local version of boost has been installed
    ([see above](#regarding-boost-and-c11)), point CMake to it by specifying
    `-DBOOST_ROOT=/path/to/boost/install`. Otherwise your local version may not
    be found, or be shadowed by an incompatible version.
  * We strongly recommend using a recent compiler that has (at least) C++11 features turned on by default. 

### Building our client codes

Our codes also use CMake to configure the build environment. The procedure is
analogous to ALPSCore's, e.g.:


    $ mkdir build && cd build
    $ cmake ../src
    $ make -jN all

Again, provide `-DCMAKE_PREFIX_PATH=/path/to/alpscore/install` if ALPSCore has
been installed in a non-standard location. Refer to the READMEs in the
subdirectories of the individual codes for additional flags to customize their
behavior.

General Usage
-------------

### Running a simulation from a parameter file

Simulation parameters may be specified in an INI-style parameter file, e.g.

    sweeps = 1000000
    thermalization = 10000
    timelimit = 18000

followed by simulation-specific parameters. The parameter files which have been
use to produce some of the figures in the main text are provided in the `parameter_files`
subdirectory of the individual codes. Assuming this file is saved as `job.ini`, the simulation is started
by running :

    $ ./qmc_worm job.ini

which would run it on a single core. After an initial thermalization phase of
10000 sweeps, the code would run for another one million sweeps while measuring
the observables after each one. The code would run for at most 5 hours (=18000
seconds), which is useful when working on a cluster with wallclock constraints.
A `timelimit` of 0 means that no time limit is imposed at all.

When the desired number of sweeps has been performed or the job ran out of time,
the results are written to `job.out.h5`.

#### Using MPI parallelization

To run the simulation on multiple cores or even nodes, use the executable
without the `_single` suffix in combination with the MPI wrapper script:

    $ mpiexec -n $NUM_MPI ./qmc_worm_mpi job.ini
    
This will start `$NUM_MPI` independent simulations with identical parameters
(with the exception of the random seed: each MPI process uses its own RNG,
seeded with the `SEED` parameter plus its MPI rank).

Each MPI process will produce independent checkpoint files, enumerated by the
MPI rank, e.g. `job.clone.h5.<RANK>`, but only one output file `job.out.h5` will
be written which contains the collected results from all MPI processes.

Periodically, the total number of measurements among all the processes will be
accumulated. If it exceeds the value specified in the `sweeps` parameter (or the
`timelimit` is reached), the simulation will terminate. Since this check
requires synchronization of the processes, it is not done after each sweep but
rather at intervals between `Tmin` and `Tmax` seconds (which may be specified in
the parameters file). Thus, when working on a cluster with wallclock
constraints, one should reserve at least `timelimit+Tmax` seconds to avoid
premature forceful termination of the job.

Keep in mind that the thermalization phase has to be done for each MPI process
independently, i.e. `$NUM_MPI` × `thermalization` sweeps will be carried out in
total before measurement samples can be taken.

### Overriding parameters and defaults

The complete list of all configurable parameters can be obtained by calling the
executables with the `--help` flag:

    $ ./qmc_worm --help
    # Simulation of the Bose-Hubbard or the XXZ model with the worm algorithm
    [...]
    Available options:
       help (bool):                         Print help message
       SEED (long int):                     PRNG seed (default value: 42)
       timelimit (unsigned long int):       time limit for the simulation (default value: 0)
       outputfile (std::string):            name of the output file (default value: qmc_worm.out.h5)
       checkpoint (std::string):            name of the checkpoint file to save to (default value: qmc_worm.clone.h5)

    [...]

Most parameters provide default values, but some don't to force the user to
consciously specify them.

Any parameter may be overridden on the command line, e.g.

    $ ./qmc_worm job.ini --sweeps=10000
    
Parameters provided on the command line take precedence over those in the
parameter file or in the parameters stored in a checkpoint (see below).

### Resuming a simulation from a checkpoint

In case the simulation terminated because it reached the `timelimit`, or the job
terminated due to e.g. a node failure, one can resume the simulation from the
checkpoint:

    $ ./qmc_worm job.clone.h5

In case MPI is used, any one checkpoint file can be specified on the command
line and the individual MPI processes will find their respective checkpoint file
automatically:

    $ mpiexec -n $NUM_MPI ./qmc_worm.clone.h5
    
However, `$NUM_MPI` needs to match the amount used in the previous run.

When resuming, the `timelimit` is basically reset. ~~In case the simulation
terminated because it had taken `sweeps` measurements but the results turned out
unsatisfactory, one can override the `sweeps` parameter on the command line to
sample further:

    $ ./qmc_worm job.clone.h5 --sweeps=10000000

Not all parameters can be overridden when resuming from a checkpoint this way.

### Understanding the output

Simulation results are output in a HDF5 data file `job.out.h5` containing the
simulation parameters, measurements of the observables along with binning
analyses. This file follows the standard format used in the
example codes of the ALPSCore project; please refer to their documentation for
details.

Acknowledgements
----------------

These codes make use of the [ALPSCore library][4], based on the original
[ALPS project][5]. ALPSCore makes use of the [HDF5 data format][6], as well as
the [Boost C++ libraries][7]. We also rely on some functionality of the 
[Froehlich polaron lecture notes][8] and [code][9], and on the [TKSVM library][10].

  [4]: http://alpscore.org/
  [5]: http://alps.comp-phys.org/
  [6]: https://www.hdfgroup.org/
  [7]: http://www.boost.org/
  [8]: https://scipost.org/SciPostPhysLectNotes.2
  [9]: https://gitlab.lrz.de/Lode.Pollet/LecturesDiagrammaticMonteCarlo
  [10]: https://gitlab.physik.uni-muenchen.de/tk-svm/tksvm-op

License
-------

Copyright © 2022  Nicolas Sadoune and Lode Pollet

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available in the
file [LICENSE.txt](LICENSE).
