/** 
 * 
 * \ingroup convection_2d
 * \file main.cpp
 * \brief main file for convection simulation
 * 
 * main file for 2D convection solver with full and pseudo spectral methods.
 * It can run with backward Euler or RK4 and three different initial conditions.
 * Additionally, you can compare the full and pseudo RK4 results.
 *
 * You can run it with
 *   ./convection2d.exe <method> <initial_condition>
 *
 * Both arguments are required
 * 
 * <method> is the method to use. The options are:
 *  - spectral_rk4 // pseudo
 *  - spectral_be  // pseudo
 *  - full_spectral_rk4 // full
 *  - compare // runs both pseudo and full and saves results from both
 * 
 * <initial_condition> is the initial condition type:
 *  - sine
 *  - gaussian
 *  - disk
 *
 * When run, this file creates a data directory for file outputs. Frames are
 * saved every save_interval steps.
 *
 * You can vary several parameters at the top of the main function
 *  - L (domain size)
 *  - n (grid size - n x n)
 *  - D (diffusion coefficient)
 *  - vx (velocity in x direction)
 *  - vy (velocity in y direction)
 *  - dt (time step)
 *  - steps (number of timesteps to simulate)
 *  - save_interval (how often to save results)
 */

#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include "convection2d.h"

int main(int argc, char* argv[]) {

    // simulation parameters
    const double L = 1.0;        // Domain size
    const size_t n = 128;        // Grid size (n x n)
    const double D = 0.01;       // Diffusion coefficient
    const double vx = 1.0;       // Velocity in x
    const double vy = 0.5;       // Velocity in y
    const double dt = 0.0001;    // Time step
    const size_t steps = 50000;  // Number of steps
    const size_t save_interval = 100; // Save every N steps


    // Expect two command-line arguments: method and initial condition
    if (argc < 3) {
            std::cerr << "Usage: " << argv[0]
            << " <spectral_rk4|spectral_be|full_spectral_rk4|compare> <sine|gaussian|disk>"
            << std::endl;
        return 1;
    }

    std::string method = argv[1];
    std::string ic_type = argv[2];

    const double dx = L / n;     // Grid spacing

    // Create output directory
    mkdir("data", 0777);

    // Allocate 2D grid
    MDArray2D<double> u(n, n);

    // Initialize the grid with chosen initial condition
    if (ic_type == "sine") {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double x = i * dx;
                double y = j * dx;
                u(i, j) = std::sin(2.0 * M_PI * x / L) * std::sin(2.0 * M_PI * y / L);
            }
        }
    } else if (ic_type == "gaussian") {
        double x0 = 0.5 * L;
        double y0 = 0.5 * L;
        double sigma = 0.05;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double x = i * dx;
                double y = j * dx;
                double r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
                u(i, j) = std::exp(-r2 / (2.0 * sigma * sigma));
            }
        }
    } else if (ic_type == "disk") {
        double x0 = 0.5 * L;
        double y0 = 0.5 * L;
        double r = 0.15;
        double w = 0.05;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double x = i * dx;
                double y = j * dx;
                double d2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
                bool in_disk = (d2 < r * r);
                bool in_slot = (std::fabs(x - x0) < w / 2 && y < y0);
                u(i, j) = (in_disk && !in_slot) ? 1.0 : 0.0;
            }
        }
    } else {
        std::cerr << "Invalid initial condition: " << ic_type
                  << ". Use 'sine', 'gaussian', or 'disk'." << std::endl;
        return 1;
    }

    // Save the initial condition
    saveSolution2D(u, "solution_0.dat", dx, 0.0);

    // Time integration loop
    for (size_t step = 1; step <= steps; ++step) {
        double current_time = step * dt;

        if (method == "spectral_rk4") {
            spectral_RK4_step_2d_convection(u, dt, D, vx, vy, L);
        } else if (method == "spectral_be") {
            spectral_BE_step_2d_convection(u, dt, D, vx, vy, L);
        } else if (method == "full_spectral_rk4") {
            full_spectral_RK4_step_2d_convection(u, dt, D, vx, vy, L);
        } else if (method == "compare") {
            // Allocate two fields
            MDArray2D<double> u_pseudo = u;  // Pseudo-spectral RK4
            MDArray2D<double> u_full = u;    // Full-spectral RK4

            for (size_t step = 1; step <= steps; ++step) {
                double current_time = step * dt;

                // Advance separately
                spectral_RK4_step_2d_convection(u_pseudo, dt, D, vx, vy, L);
                full_spectral_RK4_step_2d_convection(u_full, dt, D, vx, vy, L);

                // Save both solutions periodically
                if (step % save_interval == 0 || step == steps) {
                    std::ostringstream pseudo_fname, full_fname;
                    pseudo_fname << "pseudo_solution_" << step << ".dat";
                    full_fname << "full_solution_" << step << ".dat";

                    saveSolution2D(u_pseudo, pseudo_fname.str(), dx, current_time);
                    saveSolution2D(u_full, full_fname.str(), dx, current_time);

                    std::cout << "Saved pseudo and full solutions at time t=" << current_time << std::endl;
                }
            }
            std::cout << "Comparison simulation complete." << std::endl;
            return 0;
        } else {
            std::cerr << "Invalid method: " << method
                      << ". Use 'spectral_rk4', 'spectral_be', or 'full_spectral_rk4'." << std::endl;
            return 1;
        }

        // Save the solution periodically
        if (step % save_interval == 0 || step == steps) {
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution2D(u, fname.str(), dx, current_time);
            std::cout << "Saved " << fname.str() << " at time t=" << current_time << std::endl;
        }
    }

    std::cout << "2D advection-diffusion simulation complete." << std::endl;
    return 0;
}
