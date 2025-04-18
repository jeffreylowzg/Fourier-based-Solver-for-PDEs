#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>

// Include your spectral advection–diffusion header:
#include "convection2d.h"

int main(int argc, char* argv[]) {
    // Expect two command-line arguments: method and initial condition
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <spectral_rk4|spectral_be> <sine|gaussian>"
                  << std::endl;
        return 1;
    }
    
    std::string method = argv[1];
    std::string ic_type = argv[2];
    
    // Simulation parameters
    const double L = 1.0;      // domain size (assume [0,L] x [0,L])
    const size_t n = 128;      // grid size (n x n), must be consistent with FFT
    const double dx = L / n;   // grid spacing
    const double D = 0.01;     // diffusion coefficient
    const double vx = 1.0;     // velocity in x
    const double vy = 0.5;     // velocity in y
    const double dt = 0.0001;  // time step
    const size_t steps = 50000; // total number of steps
    const size_t save_interval = 100; // how often to save data

    // Create output directory
    mkdir("data", 0777);

    // Allocate 2D grid
    MDArray2D<double> u(n, n);

    // Initialize the grid with chosen initial condition
    if (ic_type == "sine") {
        // A simple sine wave in x and y
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double x = (double)i * dx;
                double y = (double)j * dx;
                // e.g. sine in x and y
                u(i, j) = std::sin(2.0 * M_PI * x / L)
                        * std::sin(2.0 * M_PI * y / L);
            }
        }
    } else if (ic_type == "gaussian") {
        // A Gaussian bump in the center
        double x0 = 0.5 * L;
        double y0 = 0.5 * L;
        double sigma = 0.05;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double x = (double)i * dx;
                double y = (double)j * dx;
                double r2 = (x - x0)*(x - x0) + (y - y0)*(y - y0);
                u(i, j) = std::exp(-r2 / (2.0 * sigma * sigma));
            }
        }
    } else {
        std::cerr << "Invalid initial condition: " << ic_type
                  << ". Use 'sine' or 'gaussian'." << std::endl;
        return 1;
    }

    // Save the initial condition
    saveSolution2D(u, "solution_0.dat", dx, 0.0);

    // Time integration loop
    for (size_t step = 1; step <= steps; step++) {
        double current_time = step * dt;

        if (method == "spectral_rk4") {
            // RK4 in Fourier space for advection–diffusion
            spectral_RK4_step_2d_convection(u, dt, D, vx, vy, L);
        } else if (method == "spectral_be") {
            // Backward Euler in Fourier space for advection–diffusion
            spectral_BE_step_2d_convection(u, dt, D, vx, vy, L);
        } else {
            std::cerr << "Invalid method: " << method
                      << ". Use 'spectral_rk4' or 'spectral_be'." << std::endl;
            return 1;
        }

        // Save the solution periodically
        if (step % save_interval == 0 || step == steps) {
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution2D(u, fname.str(), dx, current_time);
            std::cout << "Saved " << fname.str() << " at time t="
                      << current_time << std::endl;
        }
    }

    std::cout << "2D advection–diffusion simulation complete." << std::endl;
    return 0;
}
