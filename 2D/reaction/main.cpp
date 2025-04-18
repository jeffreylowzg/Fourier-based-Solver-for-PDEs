#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include "reaction2D.h"

// Main program for the reaction–diffusion simulation.
int main(int argc, char* argv[]) {
    // Command-line arguments:
    //   method: "spectral_rk4" or "spectral_be"
    //   ic_type: "gaussian", "gaussian16", etc.
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <spectral_rk4|spectral_be> <gaussian|gaussian16>" << std::endl;
        return 1;
    }
    
    std::string method = argv[1];
    std::string ic_type = argv[2];
    
    // Simulation parameters.
    const double L = 1.0;       // Domain size [0,L] x [0,L]
    const size_t n = 128;       // Number of grid points in each direction.
    const double dx = L / n;    // Grid spacing.
    const double D = 0.01;      // Diffusion coefficient.
    const double r = 1.0;       // Reaction rate.
    const double dt = 0.0001;   // Time step.
    const size_t steps = 50000; // Number of time steps.
    const size_t save_interval = 100; // How often to save data.
    
    // Create the output directory (if it does not exist)
    mkdir("data", 0777);

    // Allocate the grid.
    MDArray2D<double> u(n, n);
    
    // Set the initial condition.
    if (ic_type == "gaussian") {
        // A single Gaussian bump in the center.
        double x0 = 0.5 * L;
        double y0 = 0.5 * L;
        double sigma = 0.05;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                double x = i * dx;
                double y = j * dx;
                double r2 = (x - x0)*(x - x0) + (y - y0)*(y - y0);
                u(i, j) = std::exp(-r2 / (2.0 * sigma * sigma));
            }
        }
    } else if (ic_type == "gaussian16") {
        // Create 16 Gaussian spikes arranged in a 4x4 grid.
        // Clear the array (MDArray2D initializes elements to zero by default,
        // but we explicitly set them to 0 to ensure no residual values).
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                u(i, j) = 0.0;
            }
        }
        int gridSize = 4;         // 4 spikes per row => total 16 spikes.
        double sigma = 0.05;      // Standard deviation of each Gaussian.
        // Loop over the grid positions.
        // The centers are placed at positions (p+1)*L/(gridSize+1) so that they are
        // evenly spaced inside the domain.
        for (int p = 0; p < gridSize; p++) {
            for (int q = 0; q < gridSize; q++) {
                double x0 = (p + 1) * L / (gridSize + 1); // centers in x
                double y0 = (q + 1) * L / (gridSize + 1); // centers in y
                // Add this Gaussian spike to every grid point.
                for (size_t i = 0; i < n; i++) {
                    double x = i * dx;
                    for (size_t j = 0; j < n; j++) {
                        double y = j * dx;
                        double r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
                        // Sum contributions if Gaussians overlap.
                        u(i, j) += std::exp(-r2 / (2.0 * sigma * sigma));
                    }
                }
            }
        }
    } else {
        std::cerr << "Invalid initial condition: " << ic_type
                  << ". Use 'gaussian' or 'gaussian16'." << std::endl;
        return 1;
    }

    // Save the initial condition.
    saveSolution2D(u, "solution_0.dat", dx, 0.0);

    // Time integration loop.
    for (size_t step = 1; step <= steps; step++) {
        double current_time = step * dt;
        
        if (method == "spectral_rk4") {
            reactiondiffusion_RK4_step_2d(u, dt, D, r, L);
        } else if (method == "spectral_be") {
            reactiondiffusion_BE_step_2d(u, dt, D, r, L);
        } else {
            std::cerr << "Invalid method: " << method
                      << ". Use 'spectral_rk4' or 'spectral_be'." << std::endl;
            return 1;
        }
        
        // Save the solution periodically.
        if (step % save_interval == 0 || step == steps) {
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution2D(u, fname.str(), dx, current_time);
            std::cout << "Saved " << fname.str() << " at time t = "
                      << current_time << std::endl;
        }
    }
    
    std::cout << "2D reaction–diffusion simulation complete." << std::endl;
    return 0;
}
