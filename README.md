# ADAPTIVE ANDERSON SOLVER

This software package provides a Fortran implementation of Adaptive Anderson
Mixing - a modification of the well known Anderson algorithm for solving large
non-linear root finding problems. This adaptive version automatically adjusts
the mixing parameter during iterations to achieve a better convergence.

The Anderson algorithm is commonly used e.g. in the iterative search for a
self-consistent charge density (or potential) whithin the DFT
(density-functional theory) electronic structure calculations. For more
information about the package and its use, consult the paper in the doc
directory.

To build the program, check/modify Makefile.inc and run

    make

If the build succeeds, the library to be linked will appear in the lib
subdirectory.

For Python programmers, a Cython wrapper for the package is included. It should
be already compiled in-place by the above "make" command. Run

    make python

to compile the wrapper in-place (will be placed in the src directory), or

    make python_install

to install the Python package. For the "editable installation" of the Python
package, i.e. the installation that uses the source directory, run

    make python_editable_install

## Example Applications

- DFT self-consistence loop in FENNEC [1](https://www.sciencedirect.com/science/article/abs/pii/S0378475416301173) [2](https://dspace5.zcu.cz/bitstream/11025/37042/2/Novak_Vackar_Cimrman_Evaluating%20Hellmann.pdf)
- DFT self-consistence loop in the development version of SPRKKR [3](https://www.ebert.cup.uni-muenchen.de/old/index.php?option=com_content&view=article&id=8&catid=4&Itemid=7&lang=en)
