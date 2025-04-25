#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "reaction2D.h"

int main(int argc, char* argv[]) {
    // Usage: <program> <spectral_rk4|spectral_be> <gaussian|gaussian16|zero>
    //                <none|gaussian_pulse|standing_wave|traveling_gaussian|pulsed_ring>
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <spectral_rk4|spectral_be>"
                  << " <gaussian|gaussian16|zero>"
                  << " <none|gaussian_pulse|standing_wave|traveling_gaussian|pulsed_ring>"
                  << std::endl;
        std::cerr << "Source options:\n"
                  << "  none               (no source term)\n"
                  << "  gaussian_pulse     (growing Gaussian pulse at center)\n"
                  << "  standing_wave      (sinusoidal standing wave)\n"
                  << "  traveling_gaussian (Gaussian moving across domain)\n"
                  << "  pulsed_ring        (ring-shaped pulse in time window)\n";
        return 1;
    }

    std::string method      = argv[1];
    std::string ic_type     = argv[2];
    std::string source_type = argv[3];
    bool include_source     = (source_type != "none");

    // Simulation parameters
    const double L             = 1.0;      // Domain size
    const size_t n             = 128;      // Grid points per side
    const double dx            = L / n;    // Grid spacing
    const double D             = 0.01;     // Diffusion coefficient
    const double r             = 1.0;      // Reaction rate
    const double dt            = 0.0001;   // Time step
    const size_t steps         = 50000;    // Number of time steps
    const size_t save_interval = 100;      // How often to save

    // Create output directory
    mkdir("data", 0777);

    // Allocate solution array
    MDArray2D<double> u(n, n);

    // --- Initial condition ---
    if (ic_type == "gaussian") {
        double x0 = 0.5*L, y0 = 0.5*L, sigma = 0.05;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                double x = i*dx, y = j*dx;
                double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);
                u(i,j) = std::exp(-r2 / (2.0*sigma*sigma));
            }
        }
    }
    else if (ic_type == "gaussian16") {
        // zero-out then add 16 Gaussians
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                u(i,j) = 0.0;

        int gridSize = 4;
        double sigma = 0.05;
        for (int p = 0; p < gridSize; ++p) {
            for (int q = 0; q < gridSize; ++q) {
                double x0 = (p+1)*L/(gridSize+1);
                double y0 = (q+1)*L/(gridSize+1);
                for (size_t i = 0; i < n; ++i) {
                    double x = i*dx;
                    for (size_t j = 0; j < n; ++j) {
                        double y = j*dx;
                        double r2 = (x-x0)*(x-x0) + (y-y0)*(y-y0);
                        u(i,j) += std::exp(-r2 / (2.0*sigma*sigma));
                    }
                }
            }
        }
    }
    else if (ic_type == "zero") {
        // zero initial condition
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                u(i,j) = 0.0;
    }
    else {
        std::cerr << "Invalid initial condition: " << ic_type
                  << ". Use 'gaussian', 'gaussian16' or 'zero'." << std::endl;
        return 1;
    }

    // Save initial state
    saveSolution2D(u, "solution_0.dat", dx, 0.0);

    // --- Time-integration loop ---
    for (size_t step = 1; step <= steps; ++step) {
        double current_time = step * dt;

        if (include_source) {
            if (method == "spectral_rk4") {
                reactiondiffusion_RK4_step_2d_source(
                    u, dt, D, r, L,
                    current_time,
                    source_type
                );
            }
            else if (method == "spectral_be") {
                reactiondiffusion_BE_step_2d_source(
                    u, dt, D, r, L,
                    current_time,
                    source_type
                );
            }
            else {
                std::cerr << "Invalid method: " << method << std::endl;
                return 1;
            }
        }
        else {
            if (method == "spectral_rk4") {
                reactiondiffusion_RK4_step_2d(u, dt, D, r, L);
            }
            else if (method == "spectral_be") {
                reactiondiffusion_BE_step_2d(u, dt, D, r, L);
            }
            else {
                std::cerr << "Invalid method: " << method << std::endl;
                return 1;
            }
        }

        // Save output
        if (step % save_interval == 0 || step == steps) {
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution2D(u, fname.str(), dx, current_time);
            std::cout << "Saved " << fname.str()
                      << " at time t = " << current_time << std::endl;
        }
    }

    std::cout << "2D reaction–diffusion simulation complete." << std::endl;
    return 0;
}
