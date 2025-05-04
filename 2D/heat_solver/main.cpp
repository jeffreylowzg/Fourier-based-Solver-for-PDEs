/** 
 * \ingroup heat_solver_2d
 * \file main.cpp
 * \brief main file for 2D heat solver
 *
 * @details Solver can use either finite difference or pseudo-spectral methods.
 * It allows for various initial conditions and source terms.
 *
 * You can run it with
 *   ./heat_solver2d.exe <method> <initial_condition> <source_type>
 *
 * All three arguments are required
 * <method> is the numerical method to use:
 *  - fd
 *  - spectral_rk4
 *  - spectral_be
 *
 * <initial_condition> is the initial condition type:
 *  - sine
 *  - gaussian
 *  - zero
 *
 * <source_type> is the source term type:
 *  - none
 *  - gaussian_pulse
 *  - standing_wave
 *  - traveling_gaussian
 *  - pulsed_ring
 *
 * When run, this file creates a data directory for file outputs. Frames are saved
 * as "solution_<step>.dat" every save_interval steps.
 *
 * You can vary several parameters at the top of the main function
 *  - L (domain size)
 *  - n (grid points)
 *  - alpha (heat diffusion constant)
 *  - dt (time step)
 *  - steps (number of timesteps to simulate)
 *  - save_interval (how often to save results)
 */


#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "heat_solver2d.h"

int main(int argc, char* argv[]) {

    // simulation parameters
    const double L            = 1.0;
    const size_t n            = 101;
    const double alpha        = 0.01;
    const double dt           = 0.0001;
    const size_t steps        = 50000;
    const size_t save_interval= 100;

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <fd|spectral_rk4|spectral_be>"
                  << " <sine|gaussian|zero>"
                  << " <none|gaussian_pulse|standing_wave|"
                     "traveling_gaussian|pulsed_ring>"
                  << std::endl;
        return 1;
    }

    std::string method      = argv[1];
    std::string ic_type     = argv[2];
    std::string source_type = argv[3];
    bool include_source     = (source_type != "none");

    const double dx           = L / (n - 1);

    mkdir("data", 0777);
    MDArray2D<double> u(n,n);

    // Initial conditions
    if (ic_type == "sine") {
        for (size_t i=0; i<n; i++)
        for (size_t j=0; j<n; j++) {
            double x = i*dx, y = j*dx;
            u(i,j) = sin(2*M_PI*x/L)*sin(2*M_PI*y/L);
        }
    }
    else if (ic_type == "gaussian") {
        double x0 = 0.5*L, y0 = 0.5*L, sigma = 0.05;
        for (size_t i=0; i<n; i++)
        for (size_t j=0; j<n; j++) {
            double x = i*dx, y = j*dx;
            u(i,j) = exp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))
                        /(2*sigma*sigma));
        }
    }
    else if (ic_type == "zero") {
        // already zeroed by constructor, but do it explicitly:
        for (size_t i=0; i<n; i++)
        for (size_t j=0; j<n; j++)
            u(i,j) = 0.0;
    }
    else {
        std::cerr << "Invalid IC: " << ic_type
                  << ". Use 'sine','gaussian' or 'zero'.\n";
        return 1;
    }

    saveSolution2D(u, "solution_0.dat", dx, 0.0);

    // Time-integration
    for (size_t step=1; step<=steps; step++) {
        double t = step*dt;

        if (include_source) {
            if (method == "fd")
                RK4_step_2d_source(u, dt, dx, alpha, t, L, source_type);
            else if (method == "spectral_rk4")
                spectral_RK4_step_2d_source(u, dt, alpha, L, t, source_type);
            else if (method == "spectral_be")
                spectral_BE_step_2d_source(u, dt, alpha, L, t, source_type);
            else {
                std::cerr << "Invalid method: " << method << std::endl;
                return 1;
            }
        }
        else {
            if (method == "fd")
                RK4_step_2d(u, dt, dx, alpha);
            else if (method == "spectral_rk4")
                spectral_RK4_step_2d(u, dt, alpha, L);
            else if (method == "spectral_be")
                spectral_BE_step_2d(u, dt, alpha, L);
            else {
                std::cerr << "Invalid method: " << method << std::endl;
                return 1;
            }
        }

        if (step % save_interval == 0 || step == steps) {
            std::ostringstream fname;
            fname << "solution_" << step << ".dat";
            saveSolution2D(u, fname.str(), dx, t);
            std::cout << "Saved " << fname.str()
                      << " at time t=" << t << std::endl;
        }
    }

    std::cout << "2D simulation complete.\n";
    return 0;
}
