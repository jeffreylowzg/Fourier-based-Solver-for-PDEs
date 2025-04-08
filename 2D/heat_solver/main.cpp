#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include "heat_solver2d.h"

int main(int argc, char* argv[]) {
    // Require two command-line arguments: method and initial condition.
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fd|spectral_rk4|spectral_be> <sine|gaussian|source>" << std::endl;
        return 1;
    }
    
    std::string method = argv[1];
    std::string ic_type = argv[2];
    
    // Simulation parameters.
    const double L = 1.0;
    const size_t n = 101;      // grid size (n x n)
    const double dx = L / (n - 1);
    const double alpha = 0.01;
    const double dt = 0.0001;
    const size_t steps = 50000;
    const size_t save_interval = 100; // save every 100 steps

    // Create output directory.
    mkdir("data", 0777);
    
    // Initialize the 2D grid.
    MDArray2D<double> u(n, n);
    
    // Set initial condition.
    if (ic_type == "sine") {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double x = i * dx;
                double y = j * dx;
                u(i, j) = sin(2 * M_PI * x / L) * sin(2 * M_PI * y / L);
            }
        }
    } else if (ic_type == "gaussian") {
        const double x0 = L / 2.0, y0 = L / 2.0;
        const double sigma = 0.05;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double x = i * dx;
                double y = j * dx;
                u(i, j) = exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (2 * sigma * sigma));
            }
        }
    } else if (ic_type == "source") {
        // For source mode, start with zero initial conditions.
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                u(i, j) = 0.0;
            }
        }
    } else {
        std::cerr << "Invalid initial condition: " << ic_type << ". Use 'sine', 'gaussian', or 'source'." << std::endl;
        return 1;
    }
    
    // Save the initial condition.
    saveSolution2D(u, "solution_0.dat", dx, 0.0);
    
    // Time integration loop.
    for (size_t step = 1; step <= steps; step++) {
        double current_time = step * dt;
        if (ic_type == "source") {
            // Use the solver with the source term.
            if (method == "fd") {
                RK4_step_2d_source(u, dt, dx, alpha, current_time, L);
            } else if (method == "spectral_rk4") {
                spectral_RK4_step_2d_source(u, dt, alpha, L, current_time);
            } else if (method == "spectral_be") {
                spectral_BE_step_2d_source(u, dt, alpha, L, current_time);
            } else {
                std::cerr << "Invalid method for source mode: " << method << ". Use 'fd', 'spectral_rk4', or 'spectral_be'." << std::endl;
                return 1;
            }
        } else {
            // Standard evolution for non-source initial conditions.
            if (method == "spectral_rk4") {
                spectral_RK4_step_2d(u, dt, alpha, L);
            } else if (method == "spectral_be") {
                spectral_BE_step_2d(u, dt, alpha, L);
            } else if (method == "fd") {
                RK4_step_2d(u, dt, dx, alpha);
            } else {
                std::cerr << "Invalid method: " << method << ". Use 'fd', 'spectral_rk4', or 'spectral_be'." << std::endl;
                return 1;
            }
        }
        
        if (step % save_interval == 0 || step == steps) {
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution2D(u, fname.str(), dx, current_time);
            std::cout << "Saved " << fname.str() << " at time t=" << current_time << std::endl;
        }
    }
    
    std::cout << "2D simulation complete." << std::endl;
    return 0;
}

