#include "worm.hpp"
#include "copyright.hpp"
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        std::cout << "# " << worm::code_name() << std::endl;
        print_copyright(std::cout);

        // Creates the parameters for the simulation
        std::cout << "# Initializing parameters..." << std::endl;
        alps::parameters_type<worm>::type parameters(argc, argv, "/parameters");
        
        // Define parameters BEFORE parsing
        worm::define_parameters(parameters);
        
        if (parameters.help_requested(std::cout)) {
            return 0;  // Use return instead of exit
        }

        std::cout << "# Constructing worm..." << std::endl;
        worm sim(parameters);

        // Get checkpoint and output files EXPLICITLY
        std::string checkpoint_file = parameters["checkpoint"].as<std::string>();
        std::string output_file = parameters["outputfile"].as<std::string>();

        // Restore from checkpoint if available
        if (parameters.is_restored()) {
            std::cout << "# Restoring checkpoint from " << checkpoint_file << std::endl;
            try {
                sim.load(checkpoint_file);
            }
            catch (const std::exception& e) {
                std::cerr << "# ERROR restoring checkpoint: " << e.what() << std::endl;
                return 1;
            }
        }
        else {
            sim.initialize();
        }

        std::cout << std::endl << "# Simulation parameters:" << std::endl;
        sim.print_params(std::cout);

        // Run the simulation
        std::cout << "# Running simulation..." << std::endl;
        sim.run(alps::stop_callback(size_t(parameters["runtimelimit"])));

        // Test configuration
        std::cout << "# Testing configuration..." << std::endl;
        try {
            sim.test_conf();
            std::cout << "# Configuration test passed." << std::endl;
        }
        catch (const std::exception& e) {
            std::cerr << "# ERROR in configuration test: " << e.what() << std::endl;
            return 1;
        }

        // Checkpoint the simulation
        std::cout << "# Saving checkpoint..." << std::endl;
        try {
            sim.save(checkpoint_file);
        }
        catch (const std::exception& e) {
            std::cerr << "# WARNING: Failed to save checkpoint: " << e.what() << std::endl;
            // Continue to save results anyway
        }

        // Collect and save results
        alps::results_type<worm>::type results = alps::collect_results(sim);

        std::cout << "# Result:" << std::endl;
        std::cout << results << std::endl;

        std::cout << "# Saving results to " << output_file << std::endl;
        try {
            alps::hdf5::archive ar(output_file, "w");
            ar["/parameters"] << parameters;
            ar["/simulation/results"] << results;
        }
        catch (const std::exception& e) {
            std::cerr << "# ERROR saving results: " << e.what() << std::endl;
            return 1;
        }

        std::cout << "# Simulation completed successfully.\n";
        return 0;

    }
    catch (const std::exception& e) {
        std::cerr << "# FATAL ERROR: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "# FATAL ERROR: Unknown exception" << std::endl;
        return 1;
    }
    return 0;
}

