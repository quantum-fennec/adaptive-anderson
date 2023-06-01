typedef void adaptive_anderson_solver_state;

adaptive_anderson_solver_state* __adaptive_anderson_solver_MOD_adaptive_anderson_init(int* size, double* x0,
                                int* history, double* tolerance, double* alpha,
                               int* adaptive_alpha, double* delta, double* delta_per_vector, double* delta_gap,
                               double* weights, int* norm_tolerance, int* adapt_from,
                               double* collinearity_threshold, double* regularization_lambda,
                               double* restart_threshold, double* b_ii_switch_to_linear,
                               double* linear_if_cycling,
                               int *discard_first, int *forgot_first, int *forgot_from,
                               int *broyden_each, int *choose_worst,
                               int* debug_store_to_file, int* verbosity,
                               const char* read_from_file,
                               int* read_from_file_desc,
                               int read_from_file_length);

int __adaptive_anderson_solver_MOD_adaptive_anderson_step(adaptive_anderson_solver_state* state, double* residuum, double* x);
void __adaptive_anderson_solver_MOD_adaptive_anderson_end(adaptive_anderson_solver_state** state);
double __adaptive_anderson_solver_MOD_adaptive_anderson_residual_norm(adaptive_anderson_solver_state* state);
