/** main.cpp
 * 
 * This is the main file for testing convergence.
 * This file tests convergence by stepping the simulation forward in time and
 * comparing the result to the analytical solution.
 * 
 * The initial condition is a sine wave with period 1 on the domain [0,1].
 * The analytical solution is a decaying exponential.
 * 
 * The convergence test will run multiple times and you can vary the grid
 * spacing and/or the timestep over these runs
 * 
 * To use, run `./convergence_test <test> [options]`
 * where <test> is either `be` or `rk4`
 * and the options are `vary_n` and/or `vary_dt`.
 * 
 * vary_n will halve the grid spacing after each test (finer grid)
 * vary_dt will halve the timestep after each test (smaller timestep)
 * 
 * If you want to change the starting conditions of the test, you can edit the
 * top of the main function to change the grid spacing, number of timesteps,
 * and the end time of the simulation.
 * 
 * If you want to change how many iterations of tests to run, you can change the
 * value of num_tests
 */

#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include "heat_solver.h"

int main(int argc, char* argv[]) {
  
    
    /**
     * You can change these variables if you want to run different tests
     *  1/n is the spacing of the grid
     *  steps is the total number of timesteps to run
     *  end_time is the final simulated time (timestep = end_time/steps)
     *  num_tests = the number of tests to run
     */
    int n = 64;
    int steps = 256;
    const double end_time = 1.;
    const int num_tests = 5;

    if (argc < 2) {
        printf("argc=%d\n", argc);
        std::cerr << "Usage: " << argv[0] << " <be|rk4> [vary_n] [vary_dt]" << std::endl;
        std::cerr << "vary_n tests convergence as the grid gets finer" << std::endl;
        std::cerr << "vary_dt tests convergence as the timestep decreases" << std::endl;
        return 1;
    }
    
    bool vary_n = false;
    bool vary_dt = false;
    // read in method
    std::string method = argv[1];
    if ((method != "rk4") && (method != "be")) {
        std::cerr << "Method must be either `be` or `rk4`" << std::endl;
        return 1;
    }
    // check for vary_n and vary_dt
    for (int i = 2; i < argc; i++) {
        std::string s = argv[i];
        if (s == "vary_n") {
            vary_n = true;
        } else if (s == "vary_dt") {
            vary_dt = true;
        }
    }

    for (int test = 0; test < num_tests; test++) {
        printf("n = %d, dt = 1/%d, steps = %d\n", n, steps, steps);

        // Simulation parameters
        const double L = 1.0;
        const double dx = L / n; // grid spacing
        const double alpha = 0.01;
        const double dt = end_time / (double)steps;
        const double t_final = dt * steps;

        // Initialize the 2D grid.
        MDArray<double> u(n+1);
        MDArray<double> analytical(n+1);
        
        // initial condition
        for (int i = 0; i <= n; i++) {
            double x = i * dx;
            u[i] = sin(2 * M_PI * x / L);
            // analytical solution = original * decaying exponential
            analytical[i] = u[i] * exp(-4. * M_PI * M_PI * alpha * t_final);
        }

        // Time integration.
        if (method == "rk4") {
            spectral_RK4_steps(u, dt, alpha, L, steps);
        } else { // method == "be"
            spectral_BE_steps(u, dt, alpha, L, steps);
        }

        // compare against analytical solution to determine error
        double error = (u - analytical).euclidian_norm();
        printf("error = %.4f * 1e-6\n", error * 1e6);

        // make grid finer - increase n
        if (vary_n) {
            n *= 2;
        }
        // decrease timestep and increase number of steps
        if (vary_dt) {
            steps *= 2;
        }
    }
    
    std::cout << "Convergence test complete." << std::endl;
    return 0;
}
