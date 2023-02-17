# vic2d

A command-line two-dimensional fluid simulator using a novel semi-Lagrangian advection technique

![vorticies in a flow](media/sampleA.png?raw=true "vic2d simulation of vorticity rings")

## Building vic2d

We are now using CMake to build the code, though you can still use the unmaintained unix `Makefile`.

    git clone https://github.com/markstock/vic2d.git
    cd vic2d
    mkdir build
    cd buil
    cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=ON ..
    make

Or, Windows users can run the cross-compiled binary included in this distribution (`vic2d.exe`).

To build this Windows binary on a Linux machine, the following should work.

    sudo dnf install mingw64-libpng-static.noarch mingw64-gcc-gfortran.x86_64 mingw64-gcc-c++.x86_64 mingw64-gcc.x86_64 mingw64-zlib-static.noarch mingw64-zlib.noarch mingw64-libpng.noarch mingw64-libgomp.x86_64
    mkdir build_win
    cd build_win
    cmake -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++ -DCMAKE_Fortran_COMPILER=x86_64-w64-mingw32-gfortran -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_OMP=ON ..
    make

## Running vic2d

To get a list of options, run

    vic2d -h

Note that all diffusivities are coefficients, so for Reynolds number 10000, you should set momentum diffusivity with `-vd 0.0001`

Also note that the actual resolution of the output images will be one higher than that given on the command-line. So `-x 1024 -y 1024` will generate square PNG images with 1025 pixels on a side. And any input images for such a run should have 1025 pixels on a side. I know, this is not intuitive. I created it this way because for the multigrid solver to function efficiently, you should try to run at resolutions that are small multiples of large factors of 2, so 2560 (5 * 128) is much better than 2544 (159 * 16), and those numbers are easier to remember.

For some sample command lines and output, including links to YouTube videos, check out the [old vic2d page](http://markjstock.org/vic2d/).

## Theory

This code uses a vortex method solver, which means that pressure doesn't enter in the equations: the formulation is in velocity-vorticity coordinates. The velocity-vorticity inversion is accomplished with MUDPACK, a multigrid solution method. Interpolation to and from the grid is done with a 4th order kernel (M4'). Advection uses a semi-Lagrangian method, in which the vorticity at each grid point is drawn from a 4th order (Runge-Kutta) backward-looking projection. All of these methods contribute to the code's remarkable lack of numerical diffusion - you can run problems with Reynolds numbers in the hundreds of millions - though the drawback is that the simulation is not unconditionally stable.

The first implementation of a semi-Lagrangian vorticity-based solver is in Sawyer & Portnall (1963) A semi-Lagrangian method of solving the vorticity advection equation, Tellus, 15:4, 336-342, DOI: 10.3402/tellusa.v15i4.8862 [link](https://doi.org/10.3402/tellusa.v15i4.8862).

## Notes

This program can create digital output very similar to that made by the traditional technique of ebru, or paper marbling. [Here](http://digital-marbling.de/) is an excellent web site about this art-making method.

## Thanks

This code uses [libpng](https://github.com/glennrp/libpng) and [mudpack](https://www2.cisl.ucar.edu/resources/legacy/mudpack), a multigrid solver. All the rest of this ugly, ugly code is mine.

