#ADAPTIVE ANDERSON SOLVER
This software package provides a fortran implementation of
Adaptive Anderson Mixing - a modification of the well known
Anderson algorithm for solving large non-linear root finding
problems.
This adaptive version automatically adjust the mixing parameter
during iterations to achieve a better convergence.

For more informations about the package and its use, please read the paper
in the doc directory. To make the program, modify Makefile.inc
and run

    make

If the build suceeds, the library to be linked will appear in the lib subdirectory.

For Python programmers, a Cython wrapper for the package is included.
Run

		make python

to compile the wrapper inplace (will be placed in the src directory), or

		make python_install

to install the package. For the "editable installation" of the python package,
i.e. the installation that uses the source directory, you can run

		make python_editable_install

## Example Applications

- DFT self-consistence loop in FENNEC [1](https://www.sciencedirect.com/science/article/abs/pii/S0378475416301173) [2](https://dspace5.zcu.cz/bitstream/11025/37042/2/Novak_Vackar_Cimrman_Evaluating%20Hellmann.pdf)
- DFT self-consistence loop in the development version of SPRKKR [3](https://www.ebert.cup.uni-muenchen.de/old/index.php?option=com_content&view=article&id=8&catid=4&Itemid=7&lang=en)
