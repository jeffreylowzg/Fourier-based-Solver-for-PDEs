# 2D Reaction Solver
Can run with pseudo-spectral backward Euler or RK4. There are two reaction models, three initial conditions, and five source terms.

The solver works by using Strang splitting.

## Usage

You can run it with
```
  ./reaction2D.exe <method> <model> <initial_condition> <source_term>
```
Where:

`<method>` is the method to use. The options are:
 - spectral_rk4 // pseudo-spectral
 - spectral_be  // pseudo-spectral

`<model>` is the reaction model to use:
 - logistic
 - grayscott

`<initial_condition>` is the initial condition type:
 - gaussian
 - gaussian16
 - zero

`<source_term>`
 - none
 - gaussian_pulse
 - standing_wave
 - traveling_gaussian
 - pulsed_ring

All arguments are required
 
Example: `./heat_solver2d.exe spectral_rk4 sine gaussian_pulse`

When run, this file creates a data directory for file outputs. Frames are saved as "solution_<step>.dat" every save_interval steps.