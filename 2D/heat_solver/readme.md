# 2D Heat Solver
This is a 2D heat solver. Solver can use either finite difference or pseudo-spectral methods. It allows for various initial conditions and source terms.

## Usage

You can run it with
```
  ./heat_solver2d.exe <method> <initial_condition> <source_type>
```
Where:

  `<method>` is the numerical method to use:
  - fd
  - spectral_rk4
  - spectral_be

 `<initial_condition>` is the initial condition type:
  - sine
  - gaussian
  - zero

 `<source_type>` is the source term type:
  - none
  - gaussian_pulse
  - standing_wave
  - traveling_gaussian
  - pulsed_ring

All three arguments are required
 
Example: `./heat_solver2d.exe spectral_rk4 sine gaussian_pulse`

When run, this file creates a data directory for file outputs. Frames are saved as "solution_<step>.dat" every save_interval steps.