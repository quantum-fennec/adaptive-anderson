typedef void adaptive_anderson_solver_state;

adaptive_anderson_solver_state* __adaptive_anderson_solver_MOD_adaptive_anderson_init(int* size, double* x0,
                                int* history, double* tolerance, double* alpha,
                               int* adaptive_alpha, double* delta, double* delta_per_vector,
                               double* weights, int* norm_tolerance, double* collinearity_threshold,
                               double* regularization_lambda, int* adapt_from,
                               double* restart_threshold, double* b_ii_switch_to_linear,
                               double* linear_if_cycling, int* debug_store_to_file,
                               int* verbosity);

int __adaptive_anderson_solver_MOD_adaptive_anderson_step(adaptive_anderson_solver_state* state, double* residuum, double* x);
void __adaptive_anderson_solver_MOD_adaptive_anderson_end(adaptive_anderson_solver_state** state);
double __adaptive_anderson_solver_MOD_adaptive_anderson_residual_norm(adaptive_anderson_solver_state* state);
