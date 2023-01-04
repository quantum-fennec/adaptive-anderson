import sys
import numpy as np
import os

def test():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src'))
    import adaptive_anderson_solver

    brk=0
    iterations = 0
    def fce(x):
        nonlocal iterations
        iterations += 1
        fce = lambda x,y: [x**2-10+y, y**2-5]
        if brk:
           print()
           print(x, ':', fce(*x))
           if brk==2:
              breakpoint()
        return np.array(fce(*x))

    iterations = 0
    out=adaptive_anderson_solver.solve(fce, np.array([4.,4]), tolerance=1e-5)
    assert np.linalg.norm(fce(out)) < 1e-5

    num_it = iterations
    iterations = 0
    out=adaptive_anderson_solver.solve(fce, np.array([4.,4]), tolerance=1e-3)
    assert np.linalg.norm(fce(out)) < 1e-3
    assert iterations < num_it

    iterations = 0
    out=adaptive_anderson_solver.solve(fce, np.array([4.,4]), adapt_from=10, tolerance=1e-3)
    assert iterations > num_it

    iterations = 0
    out=adaptive_anderson_solver.solve(fce, np.array([4.,4]), adaptive_alpha=False, tolerance=1e-3, iterations=num_it+1)
    assert iterations > num_it

    iterations = 0
    out=adaptive_anderson_solver.solve(fce, np.array([4.,4]), alpha=0.3, tolerance=1e-5)
    out=adaptive_anderson_solver.solve(fce, np.array([4., 4.]), weights=np.array([1.,2]), tolerance=1e-5)
    out=adaptive_anderson_solver.solve(fce, np.array([4., 4.]), weights=np.array([1.,2]), tolerance=1e-5, norm_tolerance=True)
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
