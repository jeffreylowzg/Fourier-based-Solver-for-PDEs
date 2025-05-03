/**
 * \ingroup reaction_2d
 * \file main.cpp
 * \brief Main file for 2d reaction solver
 */

/** main.cpp
 * 
 * @brief main file for 2D reaction solver.
 * 
 * @details Can run with pseudo-spectral backward Euler or RK4. There are two
 * reaction models, three initial conditions, and five source terms.
 * 
 * The solver works by using Strang splitting.
 *
 * You can run with:
 *   ./reaction2D.exe <method> <model> <initial_condition> <source_term>
 *
 * All arguments are required
 * 
 * <method> is the method to use. The options are:
 *  - spectral_rk4 // pseudo-spectral
 *  - spectral_be  // pseudo-spectral
 * 
 * <model> is the reaction model to use:
 *  - logistic
 *  - grayscott
 * 
 * <initial_condition> is the initial condition type:
 *  - gaussian
 *  - gaussian16
 *  - zero
 * 
 * <source_term>
 *  - none
 *  - gaussian_pulse
 *  - standing_wave
 *  - traveling_gaussian
 *  - pulsed_ring
 *
 * When run, this file creates a data directory for file outputs. Frames are
 * saved every save_interval steps.
 *
 * You can vary several parameters at the top of the main function
 *  - L (domain size)
 *  - n (grid size - n x n)
 *  - D (diffusion coefficient)
 *  - r (logistic model coefficient)
 *  - dt (time step)
 *  - steps (number of timesteps to simulate)
 *  - save_interval (how often to save results)
 */

#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "reaction2D.h"

int main(int argc, char* argv[]) {

    // simulation parameters
    const double L             = 1.0;
    const size_t n             = 128;
    const double D             = 0.01;
    const double r             = 1.0;
    const double dt            = 0.0001;
    const size_t steps         = 50000;
    const size_t save_interval = 100;

    // usage
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <spectral_rk4|spectral_be>"
                  << " <logistic|grayscott>"
                  << " <gaussian|gaussian16|zero>"
                  << " <none|gaussian_pulse|standing_wave|traveling_gaussian|pulsed_ring>"
                  << std::endl;
        return 1;
    }

    std::string method      = argv[1];
    std::string model_type  = argv[2];
    std::string ic_type     = argv[3];
    std::string source_type = argv[4];
    bool include_source     = (source_type != "none");

    const double dx            = L / n;

    // Create output directory
    mkdir("data", 0777);

    // Allocate solution fields
    MDArray2D<double> u(n, n);
    MDArray2D<double> v(n, n); // for Gray-Scott only

    // --- Initial conditions ---
    if (model_type == "logistic") {
        if (ic_type == "gaussian") {
            double x0 = 0.5 * L, y0 = 0.5 * L, sigma = 0.05;
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j) {
                    double x = i * dx, y = j * dx;
                    double r2 = (x - x0)*(x - x0) + (y - y0)*(y - y0);
                    u(i,j) = std::exp(-r2 / (2.0 * sigma * sigma));
                }
        } else if (ic_type == "gaussian16") {
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j)
                    u(i,j) = 0.0;
            int gridSize = 4;
            double sigma = 0.05;
            for (int p = 0; p < gridSize; ++p) {
                for (int q = 0; q < gridSize; ++q) {
                    double x0 = (p + 1) * L / (gridSize + 1);
                    double y0 = (q + 1) * L / (gridSize + 1);
                    for (size_t i = 0; i < n; ++i) {
                        double x = i * dx;
                        for (size_t j = 0; j < n; ++j) {
                            double y = j * dx;
                            double r2 = (x - x0)*(x - x0) + (y - y0)*(y - y0);
                            u(i,j) += std::exp(-r2 / (2.0 * sigma * sigma));
                        }
                    }
                }
            }
        } else if (ic_type == "zero") {
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j)
                    u(i,j) = 0.0;
        } else {
            std::cerr << "Invalid initial condition: " << ic_type << std::endl;
            return 1;
        }

        // v is unused in logistic
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                v(i,j) = 0.0;

    } else if (model_type == "grayscott") {
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) {
                u(i,j) = 1.0;
                v(i,j) = 0.0;
            }
        // Central v blob
        double x0 = 0.5 * L, y0 = 0.5 * L, sigma = 0.05;
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j) {
                double x = i * dx, y = j * dx;
                double r2 = (x - x0)*(x - x0) + (y - y0)*(y - y0);
                v(i,j) = std::exp(-r2 / (2.0 * sigma * sigma));
            }
    } else {
        std::cerr << "Invalid model type: " << model_type << std::endl;
        return 1;
    }

    // Save initial state
    saveSolution2D(u, "solution_0.dat", dx, 0.0);

    // --- Time stepping loop ---
    for (size_t step = 1; step <= steps; ++step) {
        double current_time = step * dt;

        if (method == "spectral_rk4") {
            reactiondiffusion_RK4_step_2d_combined(
                u, v, dt, D, r, L, current_time,
                source_type, model_type
            );
        } else if (method == "spectral_be") {
            reactiondiffusion_BE_step_2d_combined(
                u, v, dt, D, r, L, current_time,
                source_type, model_type
            );
        } else {
            std::cerr << "Invalid method: " << method << std::endl;
            return 1;
        }

        if (step % save_interval == 0 || step == steps) {
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution2D(u, fname.str(), dx, current_time);
            std::cout << "Saved " << fname.str()
                      << " at time t = " << current_time << std::endl;
        }
    }

    std::cout << "Simulation complete." << std::endl;
    return 0;
}
