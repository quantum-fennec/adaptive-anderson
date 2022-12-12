""" This module defines a wrapper for Fortran Adaptive Anderson Solver. """

cimport numpy as np
import numpy as np
import cython

cdef extern from "adaptive_anderson_solver.h":
    ctypedef void adaptive_anderson_solver_state;

    adaptive_anderson_solver_state* __adaptive_anderson_solver_MOD_adaptive_anderson_init(int* size, double* x0,
                               int* history, double* tolerance, double* alpha,
                               int* adaptive_alpha, double* delta, double* delta_per_vector,
                               double* weights, bint* norm_tolerance, int* adapt_from,
                               double* collinearity_threshold, double* regularization_lambda,
                               double* restart_threshold, double* b_ii_switch_to_linear,
                               double* linear_if_cycling, int* debug_store_to_file,
                               int* verbosity)

    int __adaptive_anderson_solver_MOD_adaptive_anderson_step(adaptive_anderson_solver_state* state, double* residuum, double* x)
    void __adaptive_anderson_solver_MOD_adaptive_anderson_end(adaptive_anderson_solver_state** state)
    cdef double __adaptive_anderson_solver_MOD_adaptive_anderson_residual_norm(adaptive_anderson_solver_state* state)

cdef class AdaptiveAndersonSolver:
    """
    Class that wraps Adaptive Anderson solver.
    """

    cdef adaptive_anderson_solver_state* state;
    cdef double[::1] x;

    def __init__(AdaptiveAndersonSolver self, double[::1] x0, double tolerance=1e-10,
          int history=6, double alpha=0.5, int adaptive_alpha=1, double delta=1.0,
          double delta_per_vector=0.04, double[::1] weights=None,
          bint norm_tolerance=True, int adapt_from=0,
          double collinearity_threshold=1e-6,
          double regularization_lambda=0.0,
          double restart_threshold=0.0, double b_ii_switch_to_linear=0.0,
          double linear_if_cycling=0.0, int debug_store_to_file=0,
          int verbosity=0):
          r"""
          x0: np.ndarray[dtype=double]
              The array containing the initial guess $\rho^{\textrm{in}}_0$
          history: int
              History length --- the number of remembered input density/residual pairs.
          tolerance: double
              The convergence threshold. Negative value means, that the algorithm never
              stop (however, you can stop it raising StopIteration exeption)
          alpha: double
              Mixing parameter $a_0$
          adaptive\_alpha: bool
              If it is false, use the standard (non-adaptive) Anderson mixing.
          delta: double
              Adaptation coeficient $d_{\textrm{base}}$.
          delta_per_vector:double
              Adaptation coeficient $d$.
          weights: np.ndarray[double]
              If the vector of weights ($w_i$) is given, the $L^2$ norms used thorough the algorithm are evaluated as follows: $\sqrt{\sum x_i^2 w_i} $.
              To be used e.g. if components of the $\rho$ vector correspond to spatial elements with varying volumes.
          norm_tolerance: double
              If it is true, the {\tt threshold} is adjusted to be ``{\tt weights} independent'', i.e. it is multiplied by $\sqrt{\nicefrac{n}{\sum_i w_i}}$
          collinearity_threshold: double
              Residuals, whose linearly-independent component has norm lower than $\textrm{\tt collinearity_threshold}*|\tau_i|$ are omitted from minimizing.
          adapt_from: int
              Do not adapt $a_0$ in the first $\textrm{\tt adapt_from}-1$ iterations.

          reqularization_lambda: double
              A simple regularization: this coefficient is added to the diagonal of the matrix of the scalar product of the residuals. Experimental.
          restart_threshold: double
              In each step, discard the residuals whose norm is larger than $|\tau_i| / \textrm{\tt restart_threshold}$. Experimental.
          b_ii_switch_to_linear: double
              If $b_{i,i} / b_{i-1,i-1} > \textrm{\tt b_ii_switch_to_linear}$ or $< 1 / \textrm{\tt b_ii_switch_to_linear}$, switch to linear mixing.
          linear_if_cycling: double
              If $\exists j \leq i: |\rho^{\textrm{in}}_{i+1} - \rho^{\textrm{in}}_{j}| / |\rho^{\textrm{in}}_{i}| < a_i * \textrm{\tt linear_if_cycling} $, switch to linear mixing.

          debug_store_to_file: int
              If it is nonzero, fortran file handlers {\tt debug_store_to_file} and $\textrm{{\tt debug_store_to_file}}+1$ are used for storing $\rho^{\textrm{in}}_i$, {\tt weights} and $\tau_i$ to files {\tt and\_(inputs|residuals|weights).data}.
          verbosity: int
              The amount of the printed information about the mixing. All lines are prepended by {\tt AAMIX}. Zero means no output.
          """

    def __cinit__(AdaptiveAndersonSolver self, double[::1] x0, double tolerance=1e-10,
          int history=6, double alpha=0.5, int adaptive_alpha=1, double delta=1.0,
          double delta_per_vector=0.04, double[::1] weights=None,
          bint norm_tolerance=True, int adapt_from=0,
          double collinearity_threshold=1e-6,
          double regularization_lambda=0.0,
          double restart_threshold=0.0, double b_ii_switch_to_linear=0.0,
          double linear_if_cycling=0.0, int debug_store_to_file=0,
          int verbosity=0):

          cdef int size = x0.size

          self.x = x0
          self.state = __adaptive_anderson_solver_MOD_adaptive_anderson_init(&size, &x0[0],
                        history=&history, tolerance=&tolerance, alpha=&alpha,
                        adaptive_alpha=&adaptive_alpha, delta=&delta,
                        delta_per_vector = &delta_per_vector, weights=NULL if weights is None else &weights[0],
                        norm_tolerance=&norm_tolerance, adapt_from=&adapt_from,
                        collinearity_threshold=&collinearity_threshold, regularization_lambda=&regularization_lambda,
                        restart_threshold=&restart_threshold, b_ii_switch_to_linear=&b_ii_switch_to_linear,
                        linear_if_cycling=&linear_if_cycling, debug_store_to_file=&debug_store_to_file,
                        verbosity=&verbosity)

          if not self.state:
              raise ValueError('Illegal arguments passed to the Adaptive Anderson solver constructor')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef bint step(AdaptiveAndersonSolver self, double[::1] residuum, double[::1] x):
          out = __adaptive_anderson_solver_MOD_adaptive_anderson_step(self.state, &residuum[0], &x[0])
          self.x = x
          return out

    cpdef double norm(AdaptiveAndersonSolver self):
          """
          Return the L2 norm of the last residual.
          """
          return __adaptive_anderson_solver_MOD_adaptive_anderson_residual_norm(self.state)

    def __dealloc__(AdaptiveAndersonSolver self):
          if self.state is not NULL:
             __adaptive_anderson_solver_MOD_adaptive_anderson_end(&self.state)

    def solve(AdaptiveAndersonSolver self, function, iterations=0):
        """
        Run (or complete) the "selfconsisten" cycle.
        Stop, when the desired precision is reached or if the StopIteration have been raised.

        Parameters
        ----------
        function: callable
            The function, whose root is sought

        """
        x=np.copy(self.x)
        cdef int iteration=0
        try:
          while True:
              iteration+=1
              res = function(x)
              if self.step(res, x):
                  break
              if iterations>0 and iteration==iterations:
                  return x
        except StopIteration:
          pass
        return x


def solve(function, double[::1] x0, double tolerance=1e-10,
          int history=6, double alpha=0.5, int adaptive_alpha=1, double delta=1.0,
          double delta_per_vector=0.04, double[::1] weights=None,
          bint norm_tolerance=True,  double collinearity_threshold=1e-6,
          double regularization_lambda=0.0, int adapt_from=0,
          double restart_threshold=0.0, double b_ii_switch_to_linear=0.0,
          double linear_if_cycling=0.0, int debug_store_to_file=0,
          int verbosity=0, int iterations=0):
    r"""
    Adaptive Anderson mixing algorithm

    This function solves the mixing problem in the form of the root finding problem.

    Parameters
    ----------
    function: callable
        Function of the signature f(np.ndarray[double, ndim=1, size=x]) -> np.ndarray[double, ndim=1, size=x]
    x0: np.ndarray[dtype=double]
        The array containing the initial guess $\rho^{\textrm{in}}_0$
    history: int
        History length --- the number of remembered input density/residual pairs.
    tolerance: double
        The convergence threshold. Negative value means, that the algorithm never
        stop (however, you can stop it raising StopIteration exeption)
    alpha: double
        Mixing parameter $a_0$
    adaptive\_alpha: bool
        If it is false, use the standard (non-adaptive) Anderson mixing.
    delta: double
        Adaptation coeficient $d_{\textrm{base}}$.
    delta_per_vector:double
        Adaptation coeficient $d$.
    weights: np.ndarray[double]
        If the vector of weights ($w_i$) is given, the $L^2$ norms used thorough the algorithm are evaluated as follows: $\sqrt{\sum x_i^2 w_i} $.
        To be used e.g. if components of the $\rho$ vector correspond to spatial elements with varying volumes.
    norm_tolerance: double
        If it is true, the {\tt threshold} is adjusted to be ``{\tt weights} independent'', i.e. it is multiplied by $\sqrt{\nicefrac{n}{\sum_i w_i}}$
    collinearity_threshold: double
        Residuals, whose linearly-independent component has norm lower than $\textrm{\tt collinearity_threshold}*|\tau_i|$ are omitted from minimizing.
    adapt_from: int
        Do not adapt $a_0$ in the first $\textrm{\tt adapt_from}-1$ iterations.

    reqularization_lambda: double
        A simple regularization: this coefficient is added to the diagonal of the matrix of the scalar product of the residuals. Experimental.
    restart_threshold: double
        In each step, discard the residuals whose norm is larger than $|\tau_i| / \textrm{\tt restart_threshold}$. Experimental.
    b_ii_switch_to_linear: double
        If $b_{i,i} / b_{i-1,i-1} > \textrm{\tt b_ii_switch_to_linear}$ or $< 1 / \textrm{\tt b_ii_switch_to_linear}$, switch to linear mixing.
    linear_if_cycling: double
        If $\exists j \leq i: |\rho^{\textrm{in}}_{i+1} - \rho^{\textrm{in}}_{j}| / |\rho^{\textrm{in}}_{i}| < a_i * \textrm{\tt linear_if_cycling} $, switch to linear mixing.

    debug_store_to_file: int
        If it is nonzero, fortran file handlers {\tt debug_store_to_file} and $\textrm{{\tt debug_store_to_file}}+1$ are used for storing $\rho^{\textrm{in}}_i$, {\tt weights} and $\tau_i$ to files {\tt and\_(inputs|residuals|weights).data}.
    verbosity: int
        The amount of the printed information about the mixing. All lines are prepended by {\tt AAMIX}. Zero means no output.
    """
    cdef AdaptiveAndersonSolver aa = AdaptiveAndersonSolver(x0,
                                history=history, tolerance=tolerance, alpha=alpha,
                                adaptive_alpha=adaptive_alpha, delta=delta,
                                delta_per_vector = delta_per_vector, weights=weights,
                                norm_tolerance=norm_tolerance, adapt_from=adapt_from,
                                collinearity_threshold=collinearity_threshold, regularization_lambda=regularization_lambda,
                                restart_threshold=restart_threshold, b_ii_switch_to_linear=b_ii_switch_to_linear,
                                linear_if_cycling=linear_if_cycling, debug_store_to_file=debug_store_to_file,
                                verbosity=verbosity)
    return aa.solve(function, iterations=iterations)
