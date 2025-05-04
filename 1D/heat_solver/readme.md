# 1D Heat Solver

This is the main file for running a 1D heat equation simulation and writing the
output to a file. It uses the pseudo-spectral method to solve

The simulation steps an initial heat condition forward in time using RK4.
You can choose the finite difference method or spectral method and the initial
condition (a combination of sine waves or a Gaussian pulse).

The output is saved to file in the /data directory so they can be at regular
time intervals so the results can be visualized with visualize.py later

## Usage

```
./heat_sim <method> <initial_condition>
```

where:

 `<method>` is either `fd` for finite difference or `spectral` for spectral

`<initial_condition>` is either `sine` or `gaussian`.

Example: `./heat_solver.exe spectral sine`

The default simulation uses a grid of 101 points over the domain [0,1] and runs for
100000 steps with a timestep of 0.0001. If you want to use different parameters,
you can change them at the top of the main function.