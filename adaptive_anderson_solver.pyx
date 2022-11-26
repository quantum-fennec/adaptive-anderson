cimport numpy as np
import cython

cdef extern from "adaptive_anderson_solver.h":
    void * __adaptive_anderson_solver_MOD_adaptive_anderson_init(int* size, double* x0,
                               int* history, double* precision, double* alpha,
                               int* adaptive_alpha, double* delta, double* delta_per_vector,
                               double* weights, bint* norm_threshold,  double* collinearity_threshold,
                               double* regularization_lambda, int* adapt_from,
                               double* restart_threshold, double* b_ii_switch_to_linear,
                               double* linear_if_cycling, int* debug_store_to_file,
                               int* verbosity)

    int __adaptive_anderson_solver_MOD_adaptive_anderson_step(int* size, void* state, double* residuum, double* x)
    void __adaptive_anderson_solver_MOD_adaptive_anderson_end(void* state)
    void __adaptive_anderson_solver_MOD_adaptive_anderson_residuum_norm(void* state)


@cython.boundscheck(False)
@cython.wraparound(False)
def solve(fce, np.ndarray[double, ndim=1] x0, double precision=1e-10,
          int history=6, double alpha=0.5, int adaptive_alpha=1, double delta=1.0,
          double delta_per_vector=0.04, np.ndarray[double, ndim=1] weights=None,
          bint norm_threshold=True,  double collinearity_threshold=1e-6,
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

    cdef void * state
    cdef np.ndarray[double, ndim=1] res
    cdef int size = x0.size
    cdef np.ndarray[double, ndim=1] x = x0.copy()

    state = __adaptive_anderson_solver_MOD_adaptive_anderson_init(&size, &x0[0],
                                history=&history, precision=&precision, alpha=&alpha,
                                adaptive_alpha=&adaptive_alpha, delta=&delta,
                                delta_per_vector = &delta_per_vector, weights=&weights[0],
                                norm_threshold=&norm_threshold, collinearity_threshold=&collinearity_threshold,
                                regularization_lambda=&regularization_lambda, adapt_from=&adapt_from,
                                restart_threshold=&restart_threshold, b_ii_switch_to_linear=&b_ii_switch_to_linear,
                                linear_if_cycling=&linear_if_cycling, debug_store_to_file=&debug_store_to_file,
                                verbosity=&verbosity)

    try:
        while True:
            res = fce(x)
            if __adaptive_anderson_solver_MOD_adaptive_anderson_step(&size, state, &res[0], &x[0]):
                break
    except StopIteration:
        pass
    finally:
        __adaptive_anderson_solver_MOD_adaptive_anderson_end(&state)
    return x
