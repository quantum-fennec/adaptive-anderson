#ADAPTIVE ANDERSON SOLVER
This software package provides a fortran implementation of
Adaptive Anderson Mixing - a modification of the well known
Anderson algorithm for solving large non-linear root finding
problems.
This adaptive version automatically adjust the mixing parameter
to achieve a better convergence.

For more informations about the package and its use, please read the paper
in the doc directory. To make the program, modify Makefile.inc
and run

    make

If the build suceeds, the library to be linked will appear in the lib subdirectory.

For Python programmers, a Cython wrapper for the package is included.
