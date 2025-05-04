/**
 * \ingroup heat_solver_1d
 * \file main.cpp
 * \brief Main file for testing 1D Heat equation.
 */

/** main.cpp
 * 
 * This is the main file for running a 1D heat equation simulation and writing the
 * output to a file. It uses the pseudo-spectral method to solve
 * 
 * The simulation steps an initial heat condition forward in time using RK4.
 * You can choose the finite difference method or spectral method and the initial
 * condition (a combination of sine waves or a Gaussian pulse).
 * 
 * The output is saved to file in the /data directory so they can be at regular
 * time intervals so the results can be visualized with visualize.py later
 * 
 * To use, run `./heat_sim <method> <initial_condition>`
 * A method and initial conditon are required:
 *   <method> is either `fd` for finite difference or `spectral` for spectral
 *   <initial_condition> is either `sine` or `gaussian`.
 * 
 * Example: `./heat_solver.exe spectral sine`
 * 
 * The default simulation uses a grid of 101 points over the domain [0,1] and runs for
 * 100000 steps with a timestep of 0.0001. If you want to use different parameters,
 * you can change them at the top of the main function.
 * 
 */
 
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
#include "heat_solver.h"

/**
 * @brief Main entry point for running 1D heat equation simulation.
 * 
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 * @return Exit status code.
 */
int main(int argc, char* argv[]) {
    const double L = 1.0;
    const size_t n = 101;
    const double alpha = 0.01;
    const double dt = 0.0001;
    const size_t steps = 100000;
    const size_t save_interval = 100;

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fd|spectral> <sine|gaussian>" << std::endl;
        return 1;
    }

    std::string method = argv[1];
    std::string ic_type = argv[2];

    bool useSpectral = false;
    if (method == "spectral") {
        useSpectral = true;
    } else if (method != "fd") {
        std::cerr << "Invalid method: " << method << ". Use either 'fd' or 'spectral'." << std::endl;
        return 1;
    }

    if (ic_type != "sine" && ic_type != "gaussian") {
        std::cerr << "Invalid initial condition: " << ic_type << ". Use either 'sine' or 'gaussian'." << std::endl;
        return 1;
    }

    const double dx = L / (n - 1);

    std::cout << "Method: " << (useSpectral ? "spectral (RK4 in Fourier space)" : "finite-difference (RK4 in real space)") << std::endl;
    std::cout << "Initial condition: " << ic_type << std::endl;

    mkdir("data", 0777);

    MDArray<double> u(n);

    if (ic_type == "sine") {
        for (size_t i = 0; i < n; i++) {
            double x = i * dx;
            u[i] = sin(2 * M_PI * x / L) + 0.5 * sin(4 * M_PI * x / L);
        }
    } else if (ic_type == "gaussian") {
        const double x0 = L / 2.0;
        const double sigma = 0.05;
        for (size_t i = 0; i < n; i++) {
            double x = i * dx;
            u[i] = exp(-((x - x0) * (x - x0)) / (2 * sigma * sigma));
        }
    }

    saveSolution(u, "solution_0.dat", dx, 0.0);

    for (size_t step = 1; step <= steps; step++) {
        if (useSpectral) {
            spectral_RK4_step(u, dt, alpha, L);
        } else {
            RK4_step(u, dt, dx, alpha);
        }

        if (step % save_interval == 0 || step == steps) {
            double current_time = step * dt;
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution(u, fname.str(), dx, current_time);
            std::cout << "Saved " << fname.str() << " at time t=" << current_time << std::endl;
        }
    }

    std::cout << "Simulation complete." << std::endl;
    return 0;
}
