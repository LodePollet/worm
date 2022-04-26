#include "worm.hpp"
#include "copyright.hpp"
#include <iostream>
#include <alps/mc/mpiadapter.hpp>

int main(int argc, char** argv)
{
  alps::mpi::environment env(argc, argv);
  alps::mpi::communicator comm;
  const int rank=comm.rank();
  const bool is_master=(rank==0);
  
  try {
    if (is_master) {
      std::cout << "# " << worm::code_name() << std::endl;
      print_copyright(std::cout);

      // Creates the parameters for the simulation
      // If an hdf5 file is supplied, reads the parameters there
      std::cout << "# Initializing parameters..." << std::endl;
    }
    alps::parameters_type<alps::mcmpiadapter<worm>>::type parameters(argc, argv, comm);
    alps::mcmpiadapter<worm>::define_parameters(parameters);
    if (parameters.help_requested(std::cout)) {
      exit(0);
    }
    
    std::cout << "# Constructing worm on rank " << rank  << std::endl;
    alps::mcmpiadapter<worm> sim(parameters, comm);
    // If needed, restore the last checkpoint
    std::string checkpoint_file = parameters["checkpoint"].as<std::string>();
    if (!is_master) checkpoint_file+="."+std::to_string(rank);
    if (parameters.is_restored()) {
      std::cout << "# Restoring checkpoint from " << checkpoint_file << " on rank " << rank << std::endl;
      sim.load(checkpoint_file);
    }
    else {
      sim.initialize();
    }
    if (is_master) {
      std::cout << std::endl << "# Simulation parameters:" << std::endl;
      sim.print_params(std::cout);
    }
    // Run the simulation
    std::cout << "# Running simulation on rank " << rank << std::endl;

    sim.run(alps::stop_callback(size_t(parameters["runtimelimit"])));

    // Checkpoint the simulation
    std::cout << "# Checkpointing simulation after testing on rank " << rank <<  std::endl;
    try {
      sim.test_conf();
    }
    catch (const char* e) {
      cerr << e << endl;
      exit(1);
    }
    sim.save(checkpoint_file);
    alps::results_type<alps::mcmpiadapter<worm>>::type results;
    results = alps::collect_results(sim);

    if (is_master) {
      
      // Print results
      std::cout << "Simulation ran for " << results["Total_Energy"].count() << " steps." << std::endl;
      std::cout << "# Result:" << std::endl;
      std::cout << results << std::endl;

      // Saving to the output file
      std::string output_file = parameters["outputfile"];
      alps::hdf5::archive ar(output_file, "w");
      ar["/parameters"] << parameters;
      ar["/simulation/results"] << results;
      
      std::cout << "# Finished.\n";
    }
    return 0;
  } catch (const std::runtime_error& exc) {
      std::cout << "Exception caught: " << exc.what() << std::endl;
      env.abort(2);
      return 2;
  } catch (...) {
      std::cout << "Unknown exception caught." << std::endl;
      env.abort(2);
      return 2;
  }
}
