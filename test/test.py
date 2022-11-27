import sys
import numpy as np
import os

def test():
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    import adaptive_anderson_solver



    def fce(x):
        fce = lambda x,y: [x**2-10+y, y**2-5]
        return np.array(fce(*x))

    out=adaptive_anderson_solver.solve(fce, np.array([4.,4]), tolerance=1e-5)
    assert np.linalg.norm(fce(out))< 1e-5
    out=adaptive_anderson_solver.solve(fce, out, weights=np.array([1.,2]))

    try:
        out=adaptive_anderson_solver.solve(fce, out, weights=np.array([1.,-2.]), tolerance=1e-5)
        ok=True
    except ValueError:
        ok=False

    assert ok==False

    x=np.array([4.,4.])
    aa=adaptive_anderson_solver.AdaptiveAndersonSolver(x, tolerance=1e-5)
    while not aa.step(fce(x),x):
        pass
    assert np.linalg.norm(fce(x))< 1e-5

if __name__ == "__main__":
    test()
