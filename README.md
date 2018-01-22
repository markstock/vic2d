# vic2d

A command-line two-dimensional fluid simulator using a novel semi-Lagrangian advection technique

## Building vic2d

This should be easy, at least on Linux. This should be do it.

    make

## Running vic2d

To get a list of options, run

    vic2d -h

Note that all diffusivities are coefficients, so for Reynolds number 10000, you should use `-vd 0.0001` for momentum diffusivity.

For some sample command lines and output, including links to YouTube videos, check out the [old vic2d page](http://markjstock.org/vic2d/).

## Theory

This code uses a semi-Lagrangian vortex method, which means that pressure doesn't enter in the equations: the formulation is in velocity-vorticity coordinates. The velocity-vorticity inversion is accomplished with MUDPACK, a multigrid solution method. Interpolation to and from the grid is done with a 4th order kernel. Advection uses a semi-Lagrangian method, in which the vorticity at each grid point is drawn from a 4th order backward-looking interpolation.

The idea for a semi-Lagrangian vortex methods comes from a 1969 paper on atmospheric physics that I, unfortunately, cannot find any more.

## Thanks

This code uses libpng and mudpack, a multigrid solver. All the rest of this ugly, ugly code is mine.

