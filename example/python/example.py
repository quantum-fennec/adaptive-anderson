import sys
import os

""" Add path to the compiled module to the sys.path"""
path = os.path.realpath(__file__)
for i in range(3):
    path = os.path.dirname(path)
sys.path.insert(0, path)

import adaptive_anderson_solver
import numpy as np

def function(x, do_output=True):
    """ A function whose root is sought """
    function = lambda x,y: [x**2-10+y, y**2-5]
    if do_output:
        print(f"Trying {x}")
    return np.array(function(*x))

out=adaptive_anderson_solver.solve(function, np.array([4.,4]))

print(f"The sought root is in {out}, the residual norm is {np.linalg.norm(function(out, False))}")
