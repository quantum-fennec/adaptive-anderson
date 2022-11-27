cimport numpy as np
import numpy as np
import cython

cdef extern from "adaptive_anderson_solver.h":
    ctypedef void adaptive_anderson_solver_state;

    adaptive_anderson_solver_state* __adaptive_anderson_solver_MOD_adaptive_anderson_init(int* size, double* x0,
                               int* history, double* tolerance, double* alpha,
                               int* adaptive_alpha, double* delta, double* delta_per_vector,
                               double* weights, bint* norm_tolerance,  double* collinearity_threshold,
                               double* regularization_lambda, int* adapt_from,
                               double* restart_threshold, double* b_ii_switch_to_linear,
                               double* linear_if_cycling, int* debug_store_to_file,
                               int* verbosity)

    int __adaptive_anderson_solver_MOD_adaptive_anderson_step(adaptive_anderson_solver_state* state, double* residuum, double* x)
    void __adaptive_anderson_solver_MOD_adaptive_anderson_end(adaptive_anderson_solver_state** state)
    cdef double __adaptive_anderson_solver_MOD_adaptive_anderson_residual_norm(adaptive_anderson_solver_state* state)


cdef class AdaptiveAndersonSolver:

    cdef adaptive_anderson_solver_state* state;
    cdef double[::1] x0;

    def __cinit__(AdaptiveAndersonSolver self, double[::1] x0, double tolerance=1e-10,
          int history=6, double alpha=0.5, int adaptive_alpha=1, double delta=1.0,
          double delta_per_vector=0.04, double[::1] weights=None,
          bint norm_tolerance=True,  double collinearity_threshold=1e-6,
          double regularization_lambda=0.0, int adapt_from=0,
          double restart_threshold=0.0, double b_ii_switch_to_linear=0.0,
          double linear_if_cycling=0.0, int debug_store_to_file=0,
          int verbosity=0):

          cdef int size = x0.size

          self.x0 = x0
          self.state = __adaptive_anderson_solver_MOD_adaptive_anderson_init(&size, &x0[0],
                        history=&history, tolerance=&tolerance, alpha=&alpha,
                        adaptive_alpha=&adaptive_alpha, delta=&delta,
                        delta_per_vector = &delta_per_vector, weights=NULL if weights is None else &weights[0],
                        norm_tolerance=&norm_tolerance, collinearity_threshold=&collinearity_threshold,
                        regularization_lambda=&regularization_lambda, adapt_from=&adapt_from,
                        restart_threshold=&restart_threshold, b_ii_switch_to_linear=&b_ii_switch_to_linear,
                        linear_if_cycling=&linear_if_cycling, debug_store_to_file=&debug_store_to_file,
                        verbosity=&verbosity)

          if not self.state:
              raise ValueError('Illegal arguments passed to the Adaptive Anderson solver constructor')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef bint step(AdaptiveAndersonSolver self, double[::1] residuum, double[::1] x):
          return __adaptive_anderson_solver_MOD_adaptive_anderson_step(self.state, &residuum[0], &x[0])

    cpdef double norm(AdaptiveAndersonSolver self):
          return __adaptive_anderson_solver_MOD_adaptive_anderson_residual_norm(self.state)

    def __dealloc__(AdaptiveAndersonSolver self):
          if self.state is not NULL:
             __adaptive_anderson_solver_MOD_adaptive_anderson_end(&self.state)

    def solve(AdaptiveAndersonSolver self, fce):
          x=np.copy(self.x0)
          try:
            while True:
                res = fce(x)
                if self.step(res, x):
                    break
          except StopIteration:
            pass
          return x


def solve(fce, double[::1] x0, double tolerance=1e-10,
          int history=6, double alpha=0.5, int adaptive_alpha=1, double delta=1.0,
          double delta_per_vector=0.04, double[::1] weights=None,
          bint norm_tolerance=True,  double collinearity_threshold=1e-6,
          double regularization_lambda=0.0, int adapt_from=0,
          double restart_threshold=0.0, double b_ii_switch_to_linear=0.0,
          double linear_if_cycling=0.0, int debug_store_to_file=0,
          int verbosity=0):
    """
    Adaptive Anderson mixing algorithm

    This function solves the mixing problem in the form of the root finding problem.

    Parameters
    ----------
    fce: callable
        Function of the signature f(np.ndarray[double, ndim=1, size=x]) -> np.ndarray[double, ndim=1, size=x]
    """

    cdef AdaptiveAndersonSolver aa = AdaptiveAndersonSolver(x0,
                                history=history, tolerance=tolerance, alpha=alpha,
                                adaptive_alpha=adaptive_alpha, delta=delta,
                                delta_per_vector = delta_per_vector, weights=weights,
                                norm_tolerance=norm_tolerance, collinearity_threshold=collinearity_threshold,
                                regularization_lambda=regularization_lambda, adapt_from=adapt_from,
                                restart_threshold=restart_threshold, b_ii_switch_to_linear=b_ii_switch_to_linear,
                                linear_if_cycling=linear_if_cycling, debug_store_to_file=debug_store_to_file,
                                verbosity=verbosity)
    return aa.solve(fce)
