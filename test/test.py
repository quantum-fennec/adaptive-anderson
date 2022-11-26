import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import adaptive_anderson_solver
import numpy as np

def fce(x):
    fce = lambda x,y: [x**2-10+y, y**2-5]
    return np.array(fce(*x))

out=adaptive_anderson_solver.solve(fce, np.array([4.,4]))
assert np.linalg.norm(fce(out))< 1e-5
