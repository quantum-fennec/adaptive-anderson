module adaptive_anderson_solver

use, intrinsic :: iso_fortran_env, only : stdout=>output_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nonlinear root solver: solve nonlinear f(x) = 0
! using Anderson algorithm with adaptive_alpha alpha coefficient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


implicit none

private :: check_lapack, anderson_pulay_remove_residual, &
           read_find_position, &
           read_int, read_float, read_bool, &
           write_int, write_float, write_bool

!!! Holds the state of nonlinear solver
type, public :: adaptive_anderson_solver_state
   !dimension of solved problem
   integer :: n

  !OPTIONS
  !max length of history
  integer :: history
  !desired accuracy - threshold for the residuum norm
  real*8 :: tolerance
  !mixing coeficient
  real*8 :: alpha
  !adaptive change of alpha
  logical :: adaptive_alpha
  !adapt alpha to make the coefficient of the last residuum to be
  !delta + (k-1) * delta_per_vector
  !where k is the current number of used residuals
  real*8 :: delta
  real*8 :: delta_per_vector

  !If the alpha is between alpha and delta_gap, it is ok
  real*8 :: delta_gap

  !INNER STATE
  !number of used items in history
  integer :: used
  !current index in circular buffers
  integer :: current
  integer :: previous
  !circular buffer for inputs
  real*8, dimension (:,:), allocatable :: inputs
  !circular buffer for residuals
  real*8, dimension (:,:), allocatable :: residuals

  !circular buffer for residual differences (for broyden)
  real*8, dimension (:,:), allocatable :: resdiff
  !dot product of residual differences with themselves
  real*8, dimension (:), allocatable :: resdiff_dots
  !matrix for coefficients of broyden
  real*8, dimension (:,:), allocatable :: matrix

  !pulay matrix of residual scalar products for solving the optimal residuum
  real*8, dimension (:,:), allocatable :: broyden_matrix
  !the lapack qr decomposition of the matrix
  real*8, dimension (:,:), allocatable :: qr_matrix

  !use custom integration
  real*8, dimension (:), pointer :: weights

  !working place - order of the residuals this round,
  !collinear removed (zero)
  integer, allocatable :: order(:)
  !the number of non collinear vectors in the residuals matrix
  integer :: non_collinear

  !previous update of alpha was positive (+1) or negative (-1)
  integer :: adaptation_direction
  !number of changes in the direction
  integer :: adaptation_changed
  !iteration #
  integer  :: iteration

  !Do not use the vector in the current iteration, if its diagonal element in the
  !R matrix of QR factorization is lower than the threshold
  real*8 :: collinearity_threshold

  !regularization parameter, set to <= 0 to not to use
  real*8 :: regularization_lambda
  !adapt from iteration
  integer :: adapt_from

  !If greater than zero, discard all the residuals (but not the last one)
  !with norm > min(residual norms) / restart_threshold
  real*8 :: restart_threshold

  !If b_ii / b_{i-1}{i-1} is > b_ii_switch_to_linear or < 1 / b_ii_switch_to_linear
  !and the last turn was not linear, switch to the linear mixing (this turn only)
  real*8 :: b_ii_switch_to_linear
  !last turn was linear
  integer :: switched_to_linear

  !Switch to linear if norm the new input vector is
  !near an old one: the norm of the difference is < alpha * linear_if_cycling
  real*8 :: linear_if_cycling

  !integer - debug - store the mixing states to
  ! and_inputs.data
  ! and_residuals.data
  ! and_weights.data
  !using io handle debug_store_to_file and (debug_store_to_file + 1)
  integer :: debug_store_to_file

  !the amount of printed info during the mixing
  integer :: verbosity

  integer :: discard_first
  integer :: forgot_first
  integer :: forgot_from

  !do a broyden2 iteration each <broyden_each step>
  integer :: broyden_each

  !do not discard the oldest density, choose the worst from n oldest
  integer :: choose_worst

  !computed desired_alpha from delta and delta_per_vector
  real*8 :: desired_alpha
  !diagnostic  values
  real*8 :: last_bii
  real*8, allocatable :: solution(:)
  real*8, allocatable :: previous_solution(:)
  real*8 :: last_adapt_coef
  integer :: adaptation_delayed
  integer :: adaptation_count(2)
  integer :: adapted_in_iteration
  integer :: direction_changed
  integer :: direction_unchanged
  real*8 :: last_alpha
  real*8 :: minimal_alpha(2)
  real*8 :: new_minimal_alpha(2)
  integer :: minimal_alpha_from


  integer :: adjusted_at_iteration

end type adaptive_anderson_solver_state

  real*8, external :: ddot

contains

  !!! Init solver
  function adaptive_anderson_init(n, x0, history, tolerance, alpha, &
                               adaptive_alpha, delta, delta_per_vector, &
                               delta_gap, &
                               weights, norm_tolerance, adapt_from, &
                               collinearity_threshold, regularization_lambda, &
                               restart_threshold, &
                               b_ii_switch_to_linear, linear_if_cycling, &
                               discard_first, &
                               forgot_first, &
                               forgot_from, &
                               broyden_each, &
                               choose_worst, &
                               debug_store_to_file, verbosity, &
                               read_from_file, &
                               read_from_file_desc &
                               ) result(state)
     use, intrinsic :: iso_c_binding
     !number of unknowns
     integer, intent(in) :: n
     !initial (guess of the) solution
     real*8, dimension(n), intent(in) :: x0

     !number of remembered pairs of input and residuals, default 10
     integer, intent(in), optional :: history
     !mixing coefficient default 0.5
     real*8, intent(in), optional :: alpha
     !desired accuracy - threshold for residuum norm, default not applied.
     !If not applied, it is user responsibility to decide when to stop iterating.
     real*8, intent(in), optional :: tolerance
     !Adapt mixing coefficients? Default False.
     logical, intent(in), optional :: adaptive_alpha
     !Coefficients for mixing coefficient adaptation.
     !It will try to make the coefficient of the last residuum to be
     ! delta + (k-1) * delta_per_vector where k is number of residuals.
     real*8, intent(in), optional :: delta
     real*8, intent(in), optional :: delta_per_vector
     real*8, intent(in), optional :: delta_gap

     !(integral) weights (e.g. gauss quadrature weights) for residuum, used
     !for norm of the residuum.
     real*8, dimension(n), intent(in), optional, target :: weights
     !Adjust the threshold according to sum of the weights, default False
     logical, intent(in), optional :: norm_tolerance
     !Start the adaptation of alpha at the <adapt_from> iteration
     integer, intent(in), optional :: adapt_from
     !Agresivity for removing collinear vectors
     real*8, intent(in), optional :: collinearity_threshold
     !Regularization -  lambda > 0 slows convergence, but can solve instability
     real*8, intent(in), optional :: regularization_lambda

     !If greater than zero, discard all the residuals (but not the last one)
     !with norm > min(residual norms) / restart_threshold
     real*8, intent(in), optional :: restart_threshold

     !If b_ii / b_{i-1}{i-1} is > b_ii_switch_to_linear or < 1 / b_ii_switch_to_linear
     !and the last turn was not linear, switch to the linear mixing (this turn only)
     real*8, intent(in), optional :: b_ii_switch_to_linear

     !Switch to linear if norm the new input vector is
     !near an old one: the norm of the difference is < alpha * linear_if_cycling
     real*8, intent(in), optional :: linear_if_cycling

     !Discard the first n iterations (just do linear mixing, do not
     !remembering the results
     integer, intent(in), optional :: discard_first

     !If there is forgot_from iterations, start to discarding the forgot_first
     !iterations
     integer, intent(in), optional :: forgot_first
     integer, intent(in), optional :: forgot_from

     !Do broyden2 iteration each ith step
     integer, intent(in), optional :: broyden_each

     ! Forgotten density (throwed out from history) will be the one
     ! with the worst residuum from the N densities
     integer, intent(in), optional :: choose_worst
     !integer - debug - store the mixing states to files
     ! and_inputs.data
     ! and_residuals.data
     ! and_weights.data
     !using io handle debug_store_to_file and (debug_store_to_file + 1)
     integer, intent(in), optional :: debug_store_to_file

     ! > 0 : print some debug messages
!!     integer, intent(in), optional :: verbosity
     integer, intent(in), optional :: verbosity

     character(len=*, kind=C_char), intent(in), optional :: read_from_file
     integer, intent(in), optional :: read_from_file_desc



     !---- end of the arguments

     integer i
     logical ok
     real*8 r

     !returns NULL if invalid parameters, allocated solver state otherwise
     type(adaptive_anderson_solver_state), pointer :: state

     allocate(state)
     state%verbosity = merge(verbosity, 0, present(verbosity))

     if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Initialization start'

     ok = .True.
     state%n = n
     state%used = 0
     state%current = 0
     state%previous = 0


     state%discard_first = merge(discard_first, 0, present(discard_first))
     state%forgot_first = merge(forgot_first, 0, present(forgot_first))
     state%forgot_from = merge(forgot_from, 0, present(forgot_from))
     state%broyden_each = merge(broyden_each, 0, present(broyden_each))
     state%choose_worst = merge(choose_worst, 0, present(choose_worst))
     state%debug_store_to_file = merge(debug_store_to_file, 0, present(debug_store_to_file))

     if ( present(history) ) then
       if (history < 1) then
          WRITE (*,*) 'History should be at least 1'
          ok = .False.
       else
          state%history = history
       end if
     else
       state%history = 10
     end if

     if ( present( tolerance )) then
          !! negative threshold is allowed - it means no threshold (use it for your own custom
          !! stopping criterium
          if ( isnan(tolerance) ) then
            WRITE (*,*) 'Tolerance can not be NaN'
            ok = .False.
          end if
          state%tolerance = tolerance
     else
          state%tolerance = -1
     end if

     if ( present(weights)) then
       r = 0d0
       do i=1,n
        if (weights(i) <= 0d0 ) then
          WRITE (*,*) 'Weights entries should be positive, ',i,'th item is negative.'
          ok = .False.
          exit
        end if
        r = r + weights(i)
       end do
       state%weights =>  weights
       if (present(norm_tolerance)) then
          if (norm_tolerance) state%tolerance = state%tolerance * n / r
       end if
     else
       state%weights => null()
     end if

     if ( present(alpha) ) then
       if ( isnan(alpha) ) then
            WRITE (*,*) 'Alpha can not be NaN'
            ok = .False.
       end if

       state%alpha = alpha
     else
       state%alpha = 0.5D0
     end if

     if ( present(adaptive_alpha) ) then
       state%adaptive_alpha = adaptive_alpha
     else
       state%adaptive_alpha = .FALSE.
     end if

     if ( present(delta_gap) ) then
       state%delta_gap = delta_gap
     else
       state%delta_gap = 0d0
     end if


     if ( present(delta) ) then
       if ( isnan(delta) ) then
            WRITE (*,*) 'Delta can not be NaN'
            ok = .False.
       end if
       state%delta = delta
     else
       state%delta = 1.0
     end if

     if ( present(delta_per_vector) ) then
       if ( isnan(delta_per_vector) ) then
            WRITE (*,*) 'Delta_per_vector can not be NaN'
            ok = .False.
       end if
       state%delta_per_vector = delta_per_vector
     else
       state%delta_per_vector = 0.005
     end if

     if ( present(collinearity_threshold) ) then
       if ( isnan(collinearity_threshold) ) then
            WRITE (*,*) 'Collinearity_threshold can not be NaN'
            ok = .False.
       end if
       state%collinearity_threshold = collinearity_threshold
     else
       state%collinearity_threshold = 1d-10
     end if

     if ( present(regularization_lambda) ) then
       if ( isnan(regularization_lambda) ) then
            WRITE (*,*) 'Regularization_lambda can not be NaN'
            ok = .False.
       end if
       state%regularization_lambda = regularization_lambda
     else
       state%regularization_lambda = 0d0
     end if

     if ( present(adapt_from) ) then
       state%adapt_from = adapt_from
     else
       state%adapt_from = 0
     end if

     if (present(restart_threshold) ) then
       if ( isnan(restart_threshold) ) then
            WRITE (*,*) 'Restart_threshold can not be NaN'
            ok = .False.
       end if
       state%restart_threshold = restart_threshold
     else
       state%restart_threshold = 0d0
     end if

     if (.not. ok) then
        deallocate(state)
        state => null()
        return
     end if

     state%b_ii_switch_to_linear = merge(b_ii_switch_to_linear, 0d0, present(b_ii_switch_to_linear))
     state%linear_if_cycling = merge(linear_if_cycling, 0d0, present(linear_if_cycling))
     state%switched_to_linear = 0

     state%adaptation_direction = 0
     state%adaptation_changed = 0
     state%iteration = 0
     state%last_adapt_coef = 0d0
     state%last_alpha = 0d0
     state%minimal_alpha(:) = state%alpha
     state%new_minimal_alpha(:) = state%alpha
     state%minimal_alpha_from = 10

     if (state%debug_store_to_file > 0) then
       if(state%verbosity > 0) write (*,*) "AAMIX(DEBUG) - Opening file to store the mixing input in FD ", &
                    state%debug_store_to_file
        open(unit=state%debug_store_to_file, file='and_inputs.data')
        open(unit=state%debug_store_to_file+1, file='and_weights.data')
        if(present(weights)) write (state%debug_store_to_file+1,*) weights
        close( state%debug_store_to_file+1)
        open(unit=state%debug_store_to_file+1, file='and_residuals.data')
        write (state%debug_store_to_file,*) x0
     end if
     state%adaptation_delayed = 2
     state%direction_changed = 0
     state%direction_unchanged = 0
     state%adapted_in_iteration = 0
     state%adaptation_count = 0
     state%adjusted_at_iteration = 0

     state%verbosity=9
     state%delta = 0.95
     state%delta_per_vector = 0.0
     state%delta_gap=0.15
     write (*,*) state%delta, state%delta_gap

     !CUSTOM SETTINGS
     !REMOVE
     !state%verbosity=6
     !state%discard_first=0
     !state%history=40
     !state%choose_worst=20
     !state%forgot_from=6
     !state%forgot_first=1
     !state%broyden_each=2
     !state%adaptive_alpha=.true.

     !state%weights => null()
     !state%discard_first = 1
     !state%broyden_each = 1

     if (present(read_from_file)) then
        if (read_from_file .ne. "") then
          i = merge(read_from_file_desc, 1998, present(read_from_file_desc))
          open (i, file=read_from_file, action='read', status='old', form='unformatted', ACCESS="STREAM")
          call adaptive_anderson_read_config(state, i)
          close(i)
        end if
     end if

     if (state%verbosity > 0) then
        call adaptive_anderson_write_config(state, stdout)
        call adaptive_anderson_write_config(state, stdout, 'AAMIX(CONF) ')
     end if


     allocate( state%residuals( state%n, state%history ) )
     allocate( state%inputs(state%n, state%history ) )
     allocate( state%matrix(state%history, state%history) )
     allocate( state%qr_matrix(state%history, state%history) )
     allocate( state%solution(state%history) )
     allocate( state%previous_solution(state%history) )
     allocate( state%order(state%history) )
     if (state%broyden_each .ne. 0) then
         allocate( state%resdiff( state%n, state%history ) )
         allocate( state%resdiff_dots( state%history ) )
         allocate( state%broyden_matrix( state%history, state%history ) )
     end if

     state%previous_solution = 0
     call adaptive_anderson_shift_circular_buffer(state)
     call dcopy(state%n, x0(1), 1, state%inputs(1, state%current), 1)


     if( state%forgot_from > state%history) then
         !the remaining will be already forgotten
         state%forgot_first = max(0, state%forgot_first - (state%forgot_from - state%history))
     end if

     state%verbosity = 6
     if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Initialization finished'


  end function

  function adaptive_anderson_residual_norm(state)
    real*8 adaptive_anderson_residual_norm
    type(adaptive_anderson_solver_state), intent(in) :: state
    adaptive_anderson_residual_norm = sqrt(state%matrix(state%previous, state%previous))
  end function

  subroutine adaptive_anderson_adaptation_regularization(state)
    type(adaptive_anderson_solver_state) :: state
    integer i,j
    i=state%current
    j=state%previous

    if(j .eq. 0) return
    if(state%matrix(i,i) < state%matrix(j,j)/5) return

   !Regularization if the coef is lowered
   if ( &
       (state%last_adapt_coef .ne. 0) .and. (&
       (state%last_adapt_coef .le. 0.5 .and. state%alpha > state%minimal_alpha(2) / 3) .or. &
       (state%last_adapt_coef .ge. 2 .and. state%alpha < state%minimal_alpha(1) / 3) &
       )) then
      do j=2, min(state%used, state%iteration - state%adjusted_at_iteration +1)
         i = state%order(j)
         state%matrix(i,i) = state%matrix(i,i) * 1.25
      end do
      if( state%verbosity > 3) WRITE (*,*) 'AAMIX(NA R): Regularization'
      state%adjusted_at_iteration=state%iteration
   end if
  end subroutine

  !!! Given result of previous iteration, return pointer to new
  !!! input
  logical function adaptive_anderson_step(state, residuum, x)
     type(adaptive_anderson_solver_state), intent(inout), target :: state
     real*8, intent(inout), target :: residuum(state%n)
     real*8, intent(inout), optional :: x(state%n)
     !state%used can be shifted
     if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Step start'

     if (state%debug_store_to_file > 0) then
        write (state%debug_store_to_file+1,*) residuum
     end if

     if ( state%discard_first > 0 ) then
          if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Step linear'
          if(present(x)) then
            call dcopy(state%n, state%inputs(1,1), 1, x, 1)
            call daxpy(state%n, state%alpha, residuum(1), 1, x(1), 1)
            call dcopy(state%n, x, 1, state%inputs(1,1), 1)
          else
            residuum(:) = residuum * state%alpha
            residuum(:) = residuum + state%inputs(:,1)
            call dcopy(state%n, residuum, 1, state%inputs(1,1), 1)
          end if
          state%discard_first = state%discard_first - 1
          adaptive_anderson_step = .false.
          return
     end if

     if (state%verbosity > 0) then
         write(6,'(a,2i6,f12.8,2i6)') ' AAMIX(N:n,history,alpha,current):  n, history, alpha =', &
         state%n, state%history, state%alpha, state%current
     end if

     call adaptive_anderson_store_result(state, residuum)
     if (adaptive_anderson_check_convergence(state, state%tolerance)) then
        adaptive_anderson_step = .True.
        if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Step end - converged', &
            state%matrix(state%current, state%current), state%tolerance
        call adaptive_anderson_shift_circular_buffer(state)
        return
     end if

     if (state%restart_threshold > 0) then
         call adaptive_anderson_implicit_restart(state, state%tolerance)
     end if

     call adaptive_anderson_init_order(state)
     call adaptive_anderson_adaptation_regularization(state)

     !just a hack to avoid zero division in non-short-circuit evaluation
     if (state%iteration > 1 .and. state%broyden_each>0 .and. &
         mod(state%iteration, abs(state%broyden_each - 1) + 1)  == 0) then
        call adaptive_anderson_broyden_update(state, state%solution)
     else
        call adaptive_anderson_minimize_residual(state, state%solution, state%collinearity_threshold)
     end if
     state%solution(state%used+1:) = 0d0
     !Adaptation can solve the problem again, so remember the coefficients

     if(state%verbosity > 0) write (*,*) "AAMIX(O:order):", state%order(:state%non_collinear)
     if(state%verbosity > 0) write (*,*) "AAMIX(S:solution):", state%solution(:state%used)

     if (state%adaptive_alpha .and. state%iteration >= state%adapt_from) then
       call adaptive_anderson_adapt_alpha(state)
     end if

     !Set it here, since so far the value from the previous iteration have been used
     state%last_bii = state%solution(state%current)

     if(state%verbosity > 0) write (*,*) "AAMIX(D:iteration,bii,|residual|,alpha):", &
          state%iteration, state%solution(state%current), &
          sqrt(state%matrix(state%current, state%current)),state%alpha

     call adaptive_anderson_shift_circular_buffer(state)
     if (present(x)) then
        call adaptive_anderson_form_new_input(state, state%solution, x)
     else
        call adaptive_anderson_form_new_input(state, state%solution, residuum)
     end if

     adaptive_anderson_step = .False.
     if (state%debug_store_to_file > 0) then
        if (present(x)) then
          write (state%debug_store_to_file,*) x
        else 
          write (state%debug_store_to_file,*) residuum
        end if
     end if
     if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Step end'
  end function

  !!! Just return, whether the algorithm converged
  logical function adaptive_anderson_check_convergence(state, threshold)
     type(adaptive_anderson_solver_state), intent(inout) :: state
     real*8 :: norm
     real*8 :: threshold

     norm = sqrt(state%matrix(state%current, state%current))
     adaptive_anderson_check_convergence = &
        threshold > 0d0 .and.  norm < threshold
  end function

  !Just a simple formula for linear mixing
  subroutine adaptive_anderson_linear_mix(state, x)
     type(adaptive_anderson_solver_state), intent(inout) :: state
     real*8, intent(out):: x(:)
     x = state%inputs(:,state%current) + state%residuals(:, state%current) * state%alpha
     call dcopy(state%n, x, 1, state%inputs(1,state%current), 1)
  end subroutine


  !Working routine - form from the solution the new input value
  !Switching to linear mixing, if the critera are met

  subroutine adaptive_anderson_form_new_input(state, coefficients, x)
     type(adaptive_anderson_solver_state), intent(inout) :: state
     real*8, intent(in) :: coefficients(:)
     real*8, intent(out):: x(:)

     !Working variables
     real*8 :: norm, tmp
     integer :: i

     !switch to linear mixing, if b_i,i / b_i-1 is too large or too small
      if ( state%b_ii_switch_to_linear > 0 .and. state%iteration > 1 ) then
         !ratio of b_ii from the last two iterations
         norm = state%last_bii / coefficients(state%order(1))
         !if they differs too much !switch to linear
         if ( state%switched_to_linear == 0 .and. &
            (norm > state%b_ii_switch_to_linear  .or. 1./norm < state%b_ii_switch_to_linear) ) then
           if(state%verbosity > 0) write (*,*) "AAMIX(L1:Convergence => linear):", state%last_bii, 1/norm*state%last_bii, norm
           call adaptive_anderson_linear_mix(state, x)
           state%switched_to_linear = 1
           return
         end if
      end if

     !inputs = inputs @ coefficients        //@ means matrix multiply
     call dgemv('N', state%n, state%used, 1D0, state%inputs, state%n, &
           coefficients, 1, &
           0.d0, x, 1)

     !inputs += alpha * residuum @ coefficients        //@ means matrix multiply
     call dgemv('N', state%n, state%used, state%alpha, state%residuals, state%n, &
           coefficients, 1, &
           1.0d0, x, 1)

     !If the new input is too close to a previous one (norm of their difference
     !is less than norm(residual) * state%alpha * state%linear_if_cycling)
     !sitch to the linear mixing
     if (state%linear_if_cycling > 0) then
       tmp = 1d99
       norm = sqrt(norm)

       do i=1,state%used
           if (i == state%current ) cycle
           tmp = min(tmp, sqrt(sum(((state%inputs(:,i) - x)*state%weights)**2)))

           if (tmp / norm < state%alpha * state%linear_if_cycling) then
              if (state%verbosity >0 ) write(6,*) 'AAMIX(L2:Cycling => linear): ', i, tmp
              state%switched_to_linear = 2
              call adaptive_anderson_linear_mix(state, x)
              return
           end if
       end do
       if (state%verbosity > 2 ) write(6,*) 'AAMIX(LT:i,|rho_i-rho_last|,|rho_last|): ', i, tmp, norm
     end if
     state%switched_to_linear = 0
     call dcopy(state%n, x, 1, state%inputs(1,state%current), 1)
  end subroutine

  !Working routine - shift the internal counter to the new iteration
  !If a residuum (or more of them) have been removed in the middle of the
  !buffer, shift them accordingly. However, in this version of the code,
  !there is only soft removals (only for given iterations, not from the
  !history at all), so this possibility can not occure.

  subroutine adaptive_anderson_shift_circular_buffer(state)

      type(adaptive_anderson_solver_state), intent(inout) :: state

      !shift the circular buffer by (zero means replace current
      !tuple (input, residual))
      integer :: i,oi,ooi,mx, omx
      real*8 :: norm

      state%previous = state%current


      !discard the worst from choose_worst densities, instead the oldest one
      if (state%choose_worst > 0 .and. state%used == state%history) then
         norm = 0
         call adaptive_anderson_init_order(state)
         if ( state%non_collinear == state%history ) then

           do i=state%history,state%history-state%choose_worst+1, -1
              oi = state%order(i)
              if (state%matrix(oi,oi) > norm) then
                  mx = i
                  norm = state%matrix(oi,oi)
              end if
           end do

           if (mx .ne. state%order(state%used) ) then
               ooi = state%order(mx)
               do i=mx+1, state%history
                  oi = state%order(i)
                  state%matrix(ooi,:) = state%matrix(oi,:)
                  state%matrix(:,ooi) = state%matrix(:,oi)
                  call dcopy(state%n, state%residuals(1,oi), 1, state%residuals(1, ooi), 1)
                  call dcopy(state%n, state%inputs(1, oi), 1, state%inputs(1,ooi), 1)
                  if (state%broyden_each .ne. 0) then
                    call dcopy(state%n, state%resdiff(1,oi), 1, state%resdiff(1, ooi), 1)
                    state%broyden_matrix(ooi,:) = state%broyden_matrix(oi,:)
                    state%broyden_matrix(:,ooi) = state%broyden_matrix(:,oi)
                    state%resdiff_dots(ooi) = state%resdiff_dots(oi)
                  end if
                  ooi = oi
               end do
           end if
         end if
      end if


      state%current = state%current + 1
      if (state%current > state%history) then
           state%current = 1
      end if

      !some residuum(s) have been removed due to lost of the orthogonality, shift the arrays
      if (state%used < state%history .and. state%used >= state%current) then
         state%used = state%used + 1
         do i=state%used, state%current+1, -1
              state%matrix(i,:) = state%matrix(i-1,:)
              state%matrix(:,i) = state%matrix(:,i-1)
              call dcopy(state%n, state%residuals(1,i-1), 1, state%residuals(1, i), 1)
              call dcopy(state%n, state%inputs(1, i-1), 1, state%inputs(1,i), 1)
         end do
      !just not all the history has been filled yet
      elseif (state%used < state%current) then
         state%used = state%current
      end if

      state%iteration = state%iteration + 1
  end subroutine

  subroutine adaptive_anderson_store_result(state, residuum)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8, intent(in), target :: residuum(state%n)
      real*8, target :: tmp(state%n)
      real*8, pointer :: normed_residuum(:)
      integer last, curr
      integer i,j

      call dcopy(state%n, residuum, 1, state%residuals(1, state%current), 1)
      !TODO INTEGRACE
      !rs = residuals(current) * weights
      if(associated(state%weights)) then
        tmp = residuum * state%weights
        normed_residuum => tmp
      else
        normed_residuum => residuum
      end if
      state%previous_solution(:) = state%solution(:)

      !compute scalar products with residuals
!!      call dgemv('T', state%n, state%used, 1D0, state%residuals(:,:state%used), state%n,  &
!!                normed_residuum, 1,                                                        &
!!                0D0, state%matrix(2:state%used+1,state%current+1),1)
      call dgemv('T', state%n, state%used, 1D0, state%residuals(1,1), state%n,   &
                 normed_residuum, 1,                                             &
                 0D0, state%matrix(1,state%current),1)
      !u . v = v . u
      state%matrix(state%current, :) = state%matrix(:, state%current)

      if(state%broyden_each .ne. 0 .and. state%iteration > 1) then
          call dcopy(state%n, residuum, 1, state%resdiff(1, state%current), 1)
          call daxpy(state%n, -1D0, state%residuals(1, state%previous), 1, state%resdiff(1, state%current), 1)
          if(associated(state%weights)) then
              call dcopy(state%n, state%resdiff(1, state%current), 1, tmp, 1)
              state%resdiff(:, state%current) = state%resdiff(:, state%current) * state%weights
              state%resdiff_dots(state%current) = ddot(state%n, tmp, 1, state%resdiff(1, state%current), 1)
          else
              state%resdiff_dots(state%current) = &
                 ddot(state%n, state%resdiff(1, state%current), 1, state%resdiff(1, state%current), 1)
          end if

         do i=1, max(state%current, state%used)
              if ((i == 1 .and. state%used == state%current) .or. (i == state%current +1 )) cycle
              state%broyden_matrix(i, state%current) = ddot(state%n,state%resdiff(1,i), 1, residuum, 1) &
                                                     / state%resdiff_dots(i)
         end do
      end if
  end subroutine

  subroutine adaptive_anderson_broyden_update(state, solution)
      !Solves Anderson-like mixing coefficients for the Broyden2 update
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8, dimension(:) :: solution
      real*8, dimension(state%used, state%used) :: mat
      integer i,j,oio, io,jo
      integer nupdates
      real tmp

      !We use non_collinear, since some densities can be 
      !forgotten
      nupdates = state%non_collinear - 1 
      mat = 0

      oio = state%order(state%non_collinear)
      do i=nupdates,1,-1
         io = state%order(i)
         mat(io,io) = 1d0
         do j=state%non_collinear,i+1,-1
            jo = state%order(j)
            mat(:,io) = mat(:, io) - mat(:, jo) * &
                    (state%broyden_matrix(jo,io) - state%broyden_matrix(jo,oio))
         end do
         oio = io
      end do
      do i=1, state%used
      end do

      solution = 0
      jo = state%order(1)
      do i=nupdates,1,-1
         io=state%order(i)
         call daxpy(state%used, -state%broyden_matrix(io,jo), mat(1,io), 1, solution(1), 1)
      end do

      oio = state%order(state%non_collinear)
      do i=nupdates,1,-1
         io = state%order(i)
         solution(oio) = solution(oio) - solution(io)
         oio = io
      end do
      !solution(:state%n-1) = solution(2:) + solution(:state%n-1)
      solution(state%current) = solution(state%current) + 1
  end subroutine




  !! do an implicit restart - remove all 'too big' residuals:
  !! those, whose norm is > "last residual norm" / threshold
  subroutine adaptive_anderson_implicit_restart(state, threshold)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8 :: threshold
      integer :: i,shift

      threshold = sqrt(state%matrix(state%current, state%current)) / threshold

      shift = 0
      do i=1, state%used
        if ( .not. i == state%current .and. sqrt(state%matrix(i,i)) > threshold )  then
          if(state%verbosity > 0) write (*,*) "AAMIX(R:remove)", i, threshold, sqrt(state%matrix(i,i))
          shift = shift + 1
          cycle
        else if (shift > 0) then
           state%matrix(i-shift, :state%used) = state%matrix(i,:state%used)
           state%matrix(:state%used, i-shift) = state%matrix(:state%used, i)
           call dcopy(state%n, state%residuals(1,i), 1, state%residuals(1, i-shift), 1)
           call dcopy(state%n, state%inputs(1,i), 1, state%inputs(1, i-shift), 1)
           if (state%broyden_each .ne. 0) then
              call dcopy(state%n, state%resdiff(1,i), 1, state%resdiff(1, i-shift), 1)
              state%broyden_matrix(i-shift, :state%used) = state%broyden_matrix(i,:state%used)
              state%broyden_matrix(:state%used, i-shift) = state%broyden_matrix(:state%used, i)
              state%resdiff_dots(i-shift) = state%resdiff_dots(i)
           end if
           if (i == state%current) then
             state%current = state%current - shift
          end if
           if (i == state%previous) then
             state%previous = state%previous - shift
          end if
        end if
      end do
      state%used = state%used - shift
  end subroutine

  !! Return index of the previous input/residual pair
  integer function adaptive_anderson_previous_index(state)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      adaptive_anderson_previous_index = state%previous
  end function

  !!Work routine - check collinearity in QR decomposition of residual matrix
  !!If the funcion found a collinear vectors (with respect to the given threshold)
  !!it removes it from the state%order and return FALSE. Otherwise, just return TRUE
  logical function adaptive_anderson_check_collinearity_in_qr(state, threshold)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8 :: threshold
      integer i,j,oi,oj
      real*8 norm, mnorm

      if (state%verbosity > 3) write (*,*) "AAMIX(STATE): Collinearity check"

      norm = abs(state%qr_matrix(1,1))
      mnorm = max(1d-290, norm * threshold)
      !Check for the collinear vectors
      do j=2, state%non_collinear
         i=state%order(j)
         if(state%verbosity > 0) write (*,'(E16.8)', advance="no") state%qr_matrix(i,i)
         if ( abs(state%qr_matrix(i,i)) < mnorm ) then
            write (*,*) 'AAMIX(L):', state%qr_matrix(i,:)
            write (*,*) 'AAMIX(L):', state%qr_matrix(:,i)
            if(state%verbosity > 0)  then
              write (*,*)
              write (*,*) 'AAMIX(COL):', state%qr_matrix(i,i), mnorm,  &
                                                                state%qr_matrix(1,1), threshold
            end if
            call anderson_pulay_remove_residual(state, i)
            adaptive_anderson_check_collinearity_in_qr = .False.
            return
         end if
      end do
      if(state%verbosity > 0) write (*,*) 'AAMIX(CL) OK', threshold
      adaptive_anderson_check_collinearity_in_qr = .True.
  end function

  !! Work routine
  !! Set up the order array for this iteration
  subroutine adaptive_anderson_init_order(state)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      integer :: i
      integer :: forgot
      integer :: overlap

      overlap = state%used - state%current
      state%non_collinear = state%used
      !!!Order of vectors in householder orthogonalization
      if (state%iteration >= state%forgot_from) then
          forgot = min(state%iteration - state%forgot_from + 1, state%forgot_first)
          forgot = min(forgot, state%forgot_first - (state%iteration - state%history))
          if ( forgot > 0 ) then
             if (state%verbosity > 0) then
               write (*,*) "AAMIX(Forgot): First ", forgot, " densities are not considered."
             end if
             overlap = overlap - forgot
             state%non_collinear = state%non_collinear - forgot
          end if
      end if

      !!!If overlap is negative, it means, that items 1..-overlap should be forgotten
      do i=1, state%current + min(0, overlap)
         state%order(i) = state%current - i + 1
      end do
      do i=1, overlap
         state%order(i+state%current) = state%used - i + 1
      end do
  end subroutine


  !!! Just form the matrix from the state%order
  !!! (note: currently, the solving of the linear system is with pivoting
  !!! so after discarding the vectors, the order is not significant, just
  !!! it matters, whether the vector is discarded or not
  subroutine adaptive_anderson_form_matrix(state, matrix)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8, intent(inout), dimension(:,:) :: matrix

      integer :: i,j

      do i=1,state%non_collinear
         matrix(i,i) = state%matrix(state%order(i), state%order(i))
         do j=1,i-1
            matrix(i,j) = state%matrix(state%order(i),state%order(j))
            matrix(j,i) = matrix(i,j)
         end do
      end do
  end subroutine


  !! Solve problem defined by pulay matrix using modified gram-schmidt orthogonalization
  !! Discard the oldest vectors, that are (nearly) linear combination of other vectors
  !! So the MGS algorithm have to start with the newest ones
  subroutine adaptive_anderson_minimize_residual(state, solution, collinearity_threshold)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8 :: solution(state%used)
      real*8, intent(in)  :: collinearity_threshold

      real*8 :: rhs(state%used)
      real*8 :: matrix(state%used,state%used)

      integer :: lwork
      real*8  :: work(state%used * state%used + 10)
      integer :: pivot(state%used)
      real*8 :: tau(state%used)
      integer :: info
      real*8 :: norm
      integer :: i,ja

      lwork= state%used * state%used + 10

      !If a collinear vector is found in QR factorization, remove it and do the factorization again
      if (state%verbosity > 3) write (*,*) "AAMIX(STATE) Minimalization"
      do
        !Create the matrix (in the right order)
        call adaptive_anderson_form_matrix(state, state%qr_matrix)

        !regularization
        if (state%regularization_lambda > 0) then
          norm = state%qr_matrix(1,1)
          do i=1,state%non_collinear
            state%qr_matrix(i,i) = state%qr_matrix(i,i) + state%regularization_lambda * norm
          end do
        end if

        !Do QR factorization
        !call dgeqrf(state%non_collinear, state%non_collinear, state%qr_matrix, state%history, tau, work,  info)
        call dgeqrf(state%non_collinear, state%non_collinear, state%qr_matrix, state%history, tau, work, lwork, info)
        call check_lapack(info, "DGERQF")

        !If collinear vector is found, it is removed from the order and a new check is performed
        if (.not. adaptive_anderson_check_collinearity_in_qr(state, collinearity_threshold)) cycle

        rhs(:state%non_collinear) = 1D0
        !!! Do not use QR for solving - it is not accurate, even with pivoting
        !!! call DORMQR('L', 'T', state%non_collinear, 1, state%non_collinear, state%qr_matrix, state%history, tau, rhs, &
        !!!                      state%non_collinear, work, lwork, info)
        !!! call check_lapack(info, "DORMQR")

        !!! call DTRTRS('U', 'N', 'N', state%non_collinear, 1, state%qr_matrix, state%history, rhs, state%used, info)
        !!! call check_lapack(info, "DTRTRS")

        call adaptive_anderson_form_matrix(state, matrix)

        if (state%verbosity > 1) then
          do i=1,state%non_collinear
            write (*,'(A9,*(E26.18))') ' AAMIX(M)', matrix(i,:state%non_collinear)
          end do
        end if

        call dgesv(state%non_collinear, 1, matrix, state%used, pivot, rhs, state%used, info)
        if (info > 0 ) then
          if (info .eq. 1) info = state%non_collinear
          call anderson_pulay_remove_residual(state, info)
          cycle
        end if
        call check_lapack(info, "DGESV")

        if (state%verbosity > 1) write (*,*) 'AAMIX(R)', rhs(:state%non_collinear)
        !All is OK, leave the cycle
        rhs(:state%non_collinear) = rhs(:state%non_collinear) / sum(rhs(:state%non_collinear))
        if (state%verbosity > 0) write (*,*) 'AAMIX(RES)', rhs(:state%non_collinear)
        exit
      end do

      solution(:state%used) = 0D0
      do i=1,state%non_collinear
         solution(state%order(i)) = rhs(i)
      end do
   end subroutine


   ! Check, whether the adaptive coefficient shoud be accepted or not.
   function adaptive_anderson_adapt_alpha_check(state, solution, coef) result(ok)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8 :: solution(:)
      real*8, intent(inout) :: coef
      integer :: ok

      integer :: i,j
      real*8 :: r

      j = 0
      if (coef > 1.2d0) then
        do i=2, state%non_collinear
            if (1.2 * abs(solution(state%order(i))) >= abs(solution(state%order(1)))) then
               j = j + 1
               if (j >= 2) then
                  coef = sqrt(coef)
                  if (state%verbosity > 0) write (*,*) 'AAMIX(NA 5) - no adaptation', &
                                                 solution(state%order(i)), solution(state%order(1)), coef
                  !ok = 5
                  !return
                  exit
               end if
            end if
        end do
      end if

      if ( coef < state%desired_alpha  ) then
         r=0d0
         do i=2, state%non_collinear
               if (solution(state%order(i)) <  0 .and. abs(solution(state%order(i))) < coef ) then
                  r = r - solution(state%order(i)) * max(0d0, state%matrix(state%order(i), state%current)) / &
                  sqrt(state%matrix(state%order(i), state%order(i))*state%matrix(state%current, state%current))
                  !if ( r > coef .and.( state%last_adapt_coef <= 1.2 .or. state%adaptation_delayed == 1 ) ) then
                  !   if (state%verbosity > 0) write (*,*) 'AAMIX(NA 6) - no adaptation', coef, solution(state%order(i))
                  !   if (state%no_adaptation == 6) then
                  !       coef = coef ** 0.5
                  !   else
                  !       ok = 6
                  !       return
                  !   end if
                  !end if
              end if
         end do
         if (r > 0 .and. coef < state%desired_alpha) then
             write (*,*) 'AAMIX(NA 9) - correction', coef, MIN(coef + r/3, state%desired_alpha)
             coef = MIN(coef + r/2, state%desired_alpha)
         end if
      end if
      ok = 0
   end function

   function adaptive_anderson_convergence_ratio(state) result(ratio)
   !Return the ratio of the last and the current residual norms
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8 ratio
      if (state%previous >=0) then
         ratio = state%matrix(state%previous, state%previous) / state%matrix(state%current, state%current)
      end if
   end function

   !Change alpha accordingly the sollution of the minimal residual problem
   !A very empirical routine given by observation, that the anderson iteration
   !has the best convergence when the coeeficient by the last added solution
   !is cca one
   subroutine adaptive_anderson_adapt_alpha(state)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      !coefficient of last vector - to be "normed" and compared with desired value
      real*8 :: coef, last_adapt_coef, r, s, norm
      !target value for the coefficient
      real*8 :: desired

      !direction to modify the coefficient
      integer :: new_direction, res_direction
      logical :: direction_changed, aggressive
      integer :: changed, i, j, k, oi, oj, vectors
      real*8 :: histories(4)
      real*8 :: change_threshold
      real*8 :: solution(state%used)

      real*8 :: coefs
      integer :: coefs_count

      integer :: ok

      integer :: tmp(1)
      logical :: strong

      if (state%minimal_alpha_from .ne. 0 .and.  mod(state%iteration, state%minimal_alpha_from) == 0) then
         state%minimal_alpha = state%new_minimal_alpha
         state%new_minimal_alpha = state%alpha
      else
         state%new_minimal_alpha(1) = min(state%new_minimal_alpha(1), state%alpha)
         state%new_minimal_alpha(2) = max(state%new_minimal_alpha(2), state%alpha)
      end if

      histories = (/10,8,6,5/)
      solution = state%solution(:state%used)
      desired = state%delta + (state%non_collinear - 1) * min(state%delta_per_vector,6d0)
      state%desired_alpha = desired
      state%last_alpha = state%alpha

      if (state%non_collinear <= 1) then
        if (state%used > 3) then
            state%alpha = state%alpha * 3
            return
        end if
        return
      end if

      if (state%adaptation_changed == 0 .and. state%matrix(state%current,state%current) > state%matrix(1,1) &
                                        .and. state%iteration < state%history .and. state%iteration < 3) then
      !if(state%iteration < 3) then
        if (state%verbosity > 0) write (*,*) 'AAMIX(NA 1)', state%iteration,  &
                                            state%matrix(state%current,state%current), state%matrix(1,1)
        call adaptive_anderson_no_adaptation(state, 1)
        return
      end if
      coef = abs(solution(state%current))

      if (state%previous .ne. 0) then
        if ( &
           (state%previous_solution(state%previous) < 0.3) .and. &
          ( abs(state%previous_solution(state%previous)) / maxval(abs(state%previous_solution(:state%used) )) < 0.2) &
          ) then
         !if the  previous soulution is not present in the mix, the adaptation would not be related to the current alpha
          if (state%verbosity > 0) write (*,*) "AAMIX(NA 2)", state%matrix(state%current, state%current), &
                                    state%matrix(state%previous, state%previous), sqrt(coef)
          coef = sqrt(coef)
          !call adaptive_anderson_no_adaptation(state, 2)
          !return
        end if
      end if

      ok = adaptive_anderson_adapt_alpha_check(state, solution, coef)
      if (ok .ne.  0 ) then
        call adaptive_anderson_no_adaptation(state, ok)
        return
      end if


      if(state%verbosity > 0) write (*,*) "AAMIX(CO):", coef, desired, state%delta_gap

      if ( coef < desired ) then
         new_direction = -1
         coef = desired / coef
         if (coef < 0.01D0) coef = 0.01d0
      elseif ( coef > desired + state%delta_gap) then
         new_direction = 1
         coef = coef / (desired + state%delta_gap)
      else
         call adaptive_anderson_no_adaptation(state, 0)
         if(state%verbosity > 0) write (*,*) 'AAMIX(NA D) - coefficient == 1.0'
         return
      end if

     !set some debug informations
     if(state%adaptation_direction .eq. 0) then
           state%adaptation_changed=state%adaptation_changed+1
           direction_changed = .false.
           state%direction_changed = 0
           state%direction_unchanged = 0
     else
           direction_changed = state%adaptation_direction .ne. new_direction
           if (direction_changed) then
               state%direction_unchanged = 0
               state%direction_changed = state%direction_changed + 1
               state%adaptation_changed=state%adaptation_changed+1
               state%adaptation_direction = new_direction
           else
               state%direction_changed = 0
               state%direction_unchanged = state%direction_unchanged + 1
           end if
     end if
     !change_threshold = merge(4d0,2d0, state%adaptation_changed < 2)
     change_threshold = 2.0

     !do not too large changes of alpha
     if(coef > change_threshold) then
          if(state%verbosity > 0) write (*,*) "AAMIX(C: coef, reduced coef):", coef, change_threshold + log(coef / change_threshold)
          coef = change_threshold + log(coef / change_threshold)
     end if

     if (new_direction .eq. -1) coef = 1d0 / coef
     !state%alpha = state%alpha * coef

     if(state%verbosity > 1) write (*,*) "AAMIX(AD1: debug values):", &
                   solution(state%current), coef, new_direction, state%adaptation_changed, change_threshold, &
                   direction_changed, state%adaptation_direction

     if (coef > 1 .and. state%previous > 0)  then
         if (state%previous_solution(state%previous) > 0) then
            tmp = maxloc(state%previous_solution(:state%used))
         else
            tmp = minloc(state%previous_solution(:state%used))
         end if
         i = tmp(1)

         if (i .ne. state%previous .and. state%previous_solution(i) .ne. 0d0                                                   &
                                   .and. abs(state%previous_solution(state%previous) / state%previous_solution(i)) < 0.2 ) then
             if(state%verbosity > 0) write (*,*) "AAMIX(NA 3)", coef, i, state%previous_solution(i), &
                                                                state%previous_solution(state%previous)
             !call adaptive_anderson_no_adaptation(state, 3, coef)
             !return
             coef = sqrt(coef)
         end if

         if ( adaptive_anderson_convergence_ratio(state) < 0.2) then
            if(state%verbosity > 0) write (*,*) "AAMIX(NA 4)",state%matrix(state%current, state%current), &
                                                              state%matrix(state%previous, state%previous)
            coef = sqrt(coef)
            !call adaptive_anderson_no_adaptation(state, 4, coef)
           return
         end if
      elseif (coef <  1) then
         do i=2,state%non_collinear
           if (state%solution(i) > 1.2) then
            coef = sqrt(coef)
            if(state%verbosity > 0) write (*,*) "AAMIX(NA 8)", state%solution(i)
            !call adaptive_anderson_no_adaptation(state, 4, coef)
            !return
           end if
         end do
     end if

     !ZPRAC/
     if (state%last_alpha / state%alpha > 1.) then
        if (coef < 1.0 ) then
           coef = coef / max(sqrt(coef), state%alpha / state%last_alpha)
           if(state%verbosity > 0) write (*,*) "AAMIX(NA 11)", state%solution(i)
        end if
     end if

     last_adapt_coef = state%last_adapt_coef
     state%last_adapt_coef = coef
     state%adaptation_direction = new_direction


     strong = .FALSE.
     if (coef < 1.0 .and. state%previous > 0) then
       if(&
                                  state%matrix(state%current, state%current) > &
                                  state%matrix(state%previous, state%previous)/5 ) then
          !strong = .True.
       end if
    end if

     if(strong) then
       if (state%verbosity > 0) WRITE (*,*) 'AAMIX(CR) 3', coef, coef**(1d0/3d0)
       if (new_direction == 1) then
         coef = coef**(1d0/2d0)
       else
         coef = coef**(1d0/3d0)
       end if
     else
       if (state%verbosity > 0) WRITE (*,*) 'AAMIX(CR) 2', coef, coef**0.5
       if (new_direction == 1) then
         coef = coef**(1d0/1.5d0)
       else
         coef = coef**(1d0/2d0)
       end if
     end if


     if(coef < 1d0) then
          r = state%matrix(state%previous, state%previous)
          do i=2, state%used
             j = state%order(i)
             r = min(r, state%matrix(j,j))
          end do
          if(r / state%matrix(state%current, state%current) > 5 ) then
              if (state%verbosity > 0) WRITE (*,*) 'AAMIX(NA 7)', coef,r, &
                    state%matrix(state%current, state%current), &
                    r / state%matrix(state%current, state%current)
              !call adaptive_anderson_no_adaptation(state, 7, coef)
              !return
              coef = sqrt(coef)
           end if
     end if

     if (adaptive_anderson_check_raise(state, state%alpha*coef)) return

     if (coef .ne. 1d0 ) then
        state%adapted_in_iteration = state%iteration
        state%alpha = state%alpha*coef
        i = merge(1,2, coef > 1d0)
        state%adaptation_count(i) =  state%adaptation_count(i) + 1
        state%adaptation_delayed = 0
     end if


     if(state%verbosity > 0) write (*,*) "AAMIX(AD2:b_ii,coef,coef_i-1,changed,alpha,direction):", &
                          solution(state%current), coef, last_adapt_coef, &
                          changed, state%alpha, new_direction
     return

     contains
        function adaptive_anderson_check_raise(state, cap) result (corrected)
             type(adaptive_anderson_solver_state) :: state
             real*8 cap
             logical corrected
             real*8 new

             new = sqrt(state%alpha * state%last_alpha)
             write (*,*) "AAMIX(RCOR)", state%alpha, new, state%matrix(state%current, state%current) ,&
                                        state%matrix(state%previous, state%previous)

             if (state%non_collinear .le. state%used / 3 .and. state%used > state%history / 2) then
                     state%direction_changed = 1
                     state%direction_unchanged = 0
                     state%adaptation_delayed = 0
                     state%last_alpha = state%alpha
                     if(state%verbosity > 0) write (*,*) "AAMIX(NA C) CORRECTED", state%alpha, state%minimal_alpha
                     state%alpha = min(state%alpha * 10, state%minimal_alpha(1)*10)
                     corrected = .TRUE.
             else if (cap < state%minimal_alpha(2) / 10) then
                     state%direction_changed = 1
                     state%direction_unchanged = 0
                     state%adaptation_delayed = 0
                     state%last_alpha = state%alpha
                     if(state%verbosity > 0) write (*,*) "AAMIX(NA B) CORRECTED", state%alpha, state%minimal_alpha
                     state%alpha = state%minimal_alpha(2) / 10
                     corrected = .TRUE.
             else if (state%last_alpha .ne. 0d0 .and. state%alpha > state%last_alpha .and. &
                 state%matrix(state%current, state%current) > state%matrix(state%previous, state%previous) &
                 .and. new < cap) then
                     if(state%verbosity > 0) write (*,*) "AAMIX(NA A) CORRECTED", state%alpha, new
                     state%direction_changed = 1
                     state%direction_unchanged = 0
                     state%adaptation_delayed = 0
                     state%last_alpha = state%alpha
                     state%alpha = new
                     corrected = .TRUE.
             else
                     corrected = .False.
             end if
        end function

        subroutine adaptive_anderson_no_adaptation(state, why, coef)
             type(adaptive_anderson_solver_state) :: state
             integer :: why
             real*8, optional :: coef

             if (present(coef)) state%last_adapt_coef = coef
             if (adaptive_anderson_check_raise(state, state%alpha)) return
             state%direction_changed = 0
             state%direction_unchanged = 0
             state%adaptation_delayed = state%adaptation_delayed+1
        end subroutine

  end subroutine


  !Destroy solver object and do all necessary clean up
  subroutine adaptive_anderson_end(state)
       type(adaptive_anderson_solver_state), pointer :: state
       deallocate(state%residuals)
       deallocate(state%inputs)
       deallocate(state%matrix)
       deallocate(state%qr_matrix)
       deallocate(state%solution)
       deallocate(state%previous_solution)
       if (state%debug_store_to_file > 0) then
          close(state%debug_store_to_file)
          close(state%debug_store_to_file+1)
       end if
       if (state%broyden_each .ne. 0) then
           deallocate( state%resdiff )
           deallocate( state%resdiff_dots )
           deallocate( state%broyden_matrix )
       end if

       if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Cleanup'
       deallocate(state)
       state => null()
  end subroutine

  !!! Private routines !!!

  !Unified handling of lapack errors
  subroutine check_lapack(info, msg)
    integer :: info
    character(len=*) :: msg
    if (info .ne. 0) then
       write (*,*) "LAPACK CALL ", msg, " FAILED WITH CODE: ", INFO
       stop 255
    end if
  end subroutine

  !Remove the given residual from the residual matrix (for this iteration)
  subroutine anderson_pulay_remove_residual(state, i)
       type(adaptive_anderson_solver_state), intent(inout) :: state
       integer :: i
       state%non_collinear = state%non_collinear - 1
       if (i <= state%non_collinear) then
         state%order(i:state%non_collinear) = state%order(i+1:state%non_collinear+1)
       end if
  end subroutine

 function read_find_position(datas, name)
   integer :: read_find_position
   character(len=*) :: datas
   character(len=*), intent(in) :: name

   integer :: pos, npos, ln, endd
   integer :: ok

   ln = len(name)
   endd = len(datas)
   pos = 1

   do while (pos < endd)
     npos = index(datas(pos:), name)
     if ( npos == 0) return
     pos = npos + pos - 1
     if ( pos > 1 .and. ichar(datas(pos-1:)) > 32) then
        pos = pos + 1
        cycle
     end if

     pos=pos+ln+1
     if ( ichar(datas(pos-1:)) > 32) then
        pos=pos - ln
        cycle
     end if

     if ( pos >= len(datas) ) exit
     read_find_position = pos
     return
   end do
   read_find_position = 0
end function

subroutine read_int(datas, name, value)
   character(len=*) :: datas
   character(len=*), intent(in) :: name
   integer :: value
   integer :: pos, ok

   pos = read_find_position(datas, name )
   if ( pos > 0 ) read (datas(pos:) , *, iostat=ok) value
end subroutine


subroutine read_float(datas, name, value)
   character(len=*) :: datas
   character(len=*), intent(in) :: name
   real*8 :: value

   integer :: pos, ok

   pos = read_find_position(datas, name )
   if ( pos > 0 ) read (datas(pos:) , *, iostat=ok) value
end subroutine

subroutine read_bool(datas, name, value)
   character(len=*) :: datas
   character(len=*), intent(in) :: name
   logical :: value

   integer :: pos, ok, val

   pos = read_find_position(datas, name )
   if ( pos > 0 ) read (datas(pos:) , *, iostat=ok) val
   value = val .ne. 0
end subroutine

subroutine write_bool(unt, name, value)
  integer unt
  character(len=*), intent(in) :: name
  logical value

  write (unt,*) name, ' ', merge(0,1, value)
end subroutine

subroutine write_float(unt, name, value)
  integer unt
  character(len=*), intent(in) :: name
  real*8 value

  write (unt,*) name, ' ', value
end subroutine

subroutine write_int(unt, name, value)
  integer unt
  character(len=*), intent(in) :: name
  integer value

  write (unt,*) name, ' ', value
end subroutine

subroutine adaptive_anderson_read_config(state, unt)
  type(adaptive_anderson_solver_state) :: state
  integer :: unt
  integer :: file_size
  character(len=:), allocatable :: config

  inquire(unt, SIZE=file_size)
  allocate( character(len=file_size) :: config)
  read(unt) config

  call read_int(config, 'HISTORY', state%history)
  call read_int(config, 'DISCARD_FIRST', state%discard_first)
  call read_int(config, 'FORGOT_FROM', state%forgot_from)
  call read_int(config, 'FORGOT_FIRST', state%forgot_first)
  call read_int(config, 'CHOOSE_WORST', state%choose_worst)
  call read_int(config, 'BROYDEN_EACH', state%broyden_each)
  call read_float(config, 'TOLERANCE', state%tolerance)
  call read_float(config, 'ALPHA', state%alpha)
  call read_bool(config, 'ADAPTIVE_ALPHA', state%adaptive_alpha)
  call read_float(config, 'DELTA', state%delta)
  call read_float(config, 'DELTA_GAP', state%delta_gap)
  call read_float(config, 'DELTA_PER_VECTOR', state%delta_per_vector)
  call read_int(config, 'ADAPT_FROM', state%adapt_from)
  call read_float(config, 'REGULARIZATION_LAMBDA', state%regularization_lambda)
  call read_float(config, 'RESTART_TRESHOLD', state%restart_threshold)
  call read_float(config, 'B_II_SWITCH_TO_LINEAR', state%b_ii_switch_to_linear)
  call read_float(config, 'LINEAR_IF_CYCLING', state%linear_if_cycling)

  deallocate(config)

end subroutine

subroutine adaptive_anderson_write_config(state, unt , prefix)
  type(adaptive_anderson_solver_state) :: state
  integer :: unt
  character(len=*), optional :: prefix

  if( .not. present(prefix)) then
      prefix = ''
  end if

  call write_int(unt, prefix // 'HISTORY', state%history)
  call write_int(unt, prefix // 'DISCARD_FIRST', state%discard_first)
  call write_int(unt, prefix // 'FORGOT_FROM', state%forgot_from)
  call write_int(unt, prefix // 'FORGOT_FIRST', state%forgot_first)
  call write_int(unt, prefix // 'CHOOSE_WORST', state%choose_worst)
  call write_int(unt, prefix // 'BROYDEN_EACH', state%broyden_each)
  call write_float(unt, prefix // 'TOLERANCE', state%tolerance)
  call write_float(unt, prefix // 'ALPHA', state%alpha)
  call write_bool(unt, prefix // 'ADAPTIVE_ALPHA', state%adaptive_alpha)
  call write_float(unt, prefix // 'DELTA', state%delta)
  call write_float(unt, prefix // 'DELTA_GAP', state%delta_gap)
  call write_float(unt, prefix // 'DELTA_PER_VECTOR', state%delta_per_vector)
  call write_int(unt, prefix // 'ADAPT_FROM', state%adapt_from)
  call write_float(unt, prefix // 'REGULARIZATION_LAMBDA', state%regularization_lambda)
  call write_float(unt, prefix // 'RESTART_TRESHOLD', state%restart_threshold)
  call write_float(unt, prefix // 'B_II_SWITCH_TO_LINEAR', state%b_ii_switch_to_linear)
  call write_float(unt, prefix // 'LINEAR_IF_CYCLING', state%linear_if_cycling)

end subroutine

end module
