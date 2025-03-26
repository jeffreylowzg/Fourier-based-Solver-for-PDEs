// main.cpp
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
#include "heat_solver.h"

int main(int argc, char* argv[]) {
    const double L = 1.0;
    const size_t n = 101;
    const double dx = L / (n - 1);
    const double alpha = 0.01;
    const double dt = 0.0001;
    const size_t steps = 50000;
    const size_t save_interval = 100; // Save every 100 steps

    // Determine method: default to finite-difference ("fd") if not specified.
    std::string method = "fd";
    if (argc >= 2) {
        method = argv[1];
    }
    bool useSpectral = (method == "spectral");

    // Determine initial condition type:
    std::string ic_type;
    if (argc >= 3) {
        ic_type = argv[2];
    } else {
        // Default: "sine" for spectral and "gaussian" for finite-difference.
        ic_type = useSpectral ? "sine" : "gaussian";
    }
    
    std::cout << "Method: " << (useSpectral ? "spectral (FFTW)" : "finite-difference") << std::endl;
    std::cout << "Initial condition: " << ic_type << std::endl;

    // Create output directory
    mkdir("data", 0777);

    MDArray<double> u(n);

    // Set initial condition based on the chosen type:
    if (ic_type == "sine") {
        // Periodic initial condition: combination of sine waves.
        for (size_t i = 0; i < n; i++) {
            double x = i * dx;
            u[i] = sin(2 * M_PI * x / L) + 0.5 * sin(4 * M_PI * x / L);
        }
    } else if (ic_type == "gaussian") {
        // Gaussian pulse centered at L/2.
        const double x0 = L / 2.0;
        const double sigma = 0.05;
        for (size_t i = 0; i < n; i++) {
            double x = i * dx;
            u[i] = exp(-((x - x0) * (x - x0)) / (2 * sigma * sigma));
        }
    } else {
        std::cerr << "Unknown initial condition type: " << ic_type << ". Using gaussian." << std::endl;
        const double x0 = L / 2.0;
        const double sigma = 0.05;
        for (size_t i = 0; i < n; i++) {
            double x = i * dx;
            u[i] = exp(-((x - x0) * (x - x0)) / (2 * sigma * sigma));
        }
    }

    // Save initial condition
    saveSolution(u, "solution_0.dat", dx, 0.0);

    // Time integration loop
    for (size_t step = 1; step <= steps; step++) {
        if (useSpectral) {
            // Evolve in Fourier space using the spectral method.
            spectral_step(u, dt, alpha, L);
        } else {
            // Use finite-difference method (RK4 integration).
            RK4_step(u, dt, dx, alpha);
        }

        // Save the solution at intervals.
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
