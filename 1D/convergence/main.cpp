#include <iostream>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include "heat_solver.h"

int main(int argc, char* argv[]) {
  
  if (argc < 2) {
    printf("argc=%d\n", argc);
    std::cerr << "Usage: " << argv[0] << " <be|rk4> [vary_n] [vary_dt]" << std::endl;
    std::cerr << "vary_n tests convergence as the grid gets finer" << std::endl;
    std::cerr << "vary_dt tests convergence as the timestep decreases" << std::endl;
    return 1;
  }
  
  bool vary_n = false;
  bool vary_dt = false;
  std::string method = argv[1];
  if ((method != "rk4") && (method != "be")) {
    std::cerr << "Method must be either be or rk4" << std::endl;
    return 1;
  }
  for (int i = 2; i < argc; i++) {
    std::string s = argv[i];
    if (s == "vary_n") {
      vary_n = true;
    } else if (s == "vary_dt") {
      vary_dt = true;
    }
  }

  // starting values for grid and timestep - the solver runs until t = 1
  int n = 512;
  int steps = 10;
  double end_time = 1;
  for (int test = 0; test < 5; test++) {
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
    for (int i = 0; i < n; i++) {
      double x = i * dx;
      u[i] = sin(2 * M_PI * x / L);
      // analytical solution = original * decaying exponential
      analytical[i] = u[i] * exp(-4. * alpha * M_PI * M_PI * t_final);
    }

    // Time integration loop.
    for (int step = 1; step <= steps; step++) {
      if (method == "rk4") {
        spectral_RK4_step(u, dt, alpha, L);
      } else { // method == "be"
        spectral_BE_step(u, dt, alpha, L);
      }
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
      
  //     if (step % save_interval == 0 || step == steps) {
  //         std::ostringstream fname;
  //         fname << "solution_" << step << ".dat";
  //         saveSolution2D(u, fname.str(), dx, current_time);
  //         std::cout << "Saved " << fname.str() << " at time t=" << current_time << std::endl;
  //     }
  // }
  
  // std::cout << "2D simulation complete." << std::endl;
  return 0;
}
