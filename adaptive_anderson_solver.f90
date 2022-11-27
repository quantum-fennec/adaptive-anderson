module adaptive_anderson_solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nonlinear root solver: solve nonlinear f(x) = 0
! using Anderson algorithm with adaptive_alpha alpha coefficient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


implicit none

private :: check_lapack

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

  !pulay matrix of residual scalar products for solving the optimal residuum
  real*8, dimension (:,:), allocatable :: matrix
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
  integer*8 :: switched_to_linear

  !Switch to linear if norm the new input vector is
  !near an old one: the norm of the difference is < alpha * linear_if_cycling
  real*8 :: linear_if_cycling

  !integer - debug - store the mixing states to
  ! and_inputs.data
  ! and_residuals.data
  ! and_weights.data
  !using io handle debug_store_to_file and (debug_store_to_file + 1)
  integer*8 :: debug_store_to_file

  !the amount of printed info during the mixing
  integer*8 :: verbosity

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

end type adaptive_anderson_solver_state

contains

  !!! Init solver
  function adaptive_anderson_init(n, x0, history, tolerance, alpha, &
                               adaptive_alpha, delta, delta_per_vector, &
                               weights, norm_tolerance, adapt_from, &
                               collinearity_threshold, regularization_lambda, &
                               restart_threshold, &
                               b_ii_switch_to_linear, linear_if_cycling, &
                               debug_store_to_file, verbosity &
                               ) result(state)
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

     !integer - debug - store the mixing states to
     ! and_inputs.data
     ! and_residuals.data
     ! and_weights.data
     !using io handle debug_store_to_file and (debug_store_to_file + 1)
     integer, intent(in), optional :: debug_store_to_file

     ! > 0 : print some debug messages
     integer, intent(in), optional :: verbosity

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

     state%debug_store_to_file = merge(debug_store_to_file, 0,present(debug_store_to_file))

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
     state%history = 10

     if ( present( tolerance )) then
          !! negative threshold is allowed - it means no threshold (use it for your own custom
          !! stopping criterium
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
       if (alpha <= 0d0 ) then
          WRITE (*,*) 'Alpha should be positive, it is actually:', alpha
          ok = .False.
       else
          state%alpha = alpha
       end if
     else
       state%alpha = 0.5D0
     end if

     if ( present(adaptive_alpha) ) then
       state%adaptive_alpha = adaptive_alpha
     else
       state%adaptive_alpha = .FALSE.
     end if

     if ( present(delta) ) then
       state%delta = delta
     else
       state%delta = 1.0
     end if

     if ( present(delta_per_vector) ) then
       state%delta_per_vector = delta_per_vector
     else
       state%delta_per_vector = 0.005
     end if

     if ( present(collinearity_threshold) ) then
       state%collinearity_threshold = collinearity_threshold
     else
       state%collinearity_threshold = 1d-10
     end if
     state%collinearity_threshold = 1d-10

     if ( present(regularization_lambda) ) then
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

     allocate( state%residuals( state%n, state%history ) )
     allocate( state%inputs(state%n, state%history ) )
     allocate( state%matrix(state%history, state%history) )
     allocate( state%qr_matrix(state%history, state%history) )
     allocate( state%solution(state%history) )
     allocate( state%previous_solution(state%history) )
     allocate( state%order(state%history) )

     state%adaptation_direction = 0
     state%adaptation_changed = 0
     state%iteration = 0
     state%last_adapt_coef = 0d0
     state%last_alpha = 0d0

     call adaptive_anderson_shift_circular_buffer(state)
     call dcopy(state%n, x0(1), 1, state%inputs(1, state%current), 1)

     if (state%debug_store_to_file > 0) then
       if(state%verbosity > 0) write (*,*) "AAMIX(DEBUG) - Opening file to store the mixing input in FD ", &
                    state%debug_store_to_file
        open(unit=state%debug_store_to_file, file='and_inputs.data')
        open(unit=state%debug_store_to_file+1, file='and_weights.data')
        if(associated(state%weights)) write (state%debug_store_to_file+1,*) state%weights
        close( state%debug_store_to_file+1)
        open(unit=state%debug_store_to_file+1, file='and_residuals.data')
        write (state%debug_store_to_file,*) x0
     end if
     state%adaptation_delayed = 2
     state%direction_changed = 0
     state%direction_unchanged = 0
     state%adapted_in_iteration = 0
     state%previous_solution = 0
     state%adaptation_count = 0
     if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Initialization finished'
  end function

  function adaptive_anderson_residual_norm(state)
    real*8 adaptive_anderson_residual_norm
    type(adaptive_anderson_solver_state), intent(in) :: state
    adaptive_anderson_residual_norm = sqrt(state%matrix(state%previous, state%previous))
  end function

  !!! Given result of previous iteration, return pointer to new
  !!! input
  logical function adaptive_anderson_step(state, residuum, x)
     type(adaptive_anderson_solver_state), intent(inout), target :: state
     real*8, intent(in), target :: residuum(state%n)
     real*8, intent(out) :: x(state%n)
     !state%used can be shifted
     if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Step start'

     if (state%debug_store_to_file > 0) then
        write (state%debug_store_to_file+1,*) residuum
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
     call adaptive_anderson_minimize_residual(state, state%solution, state%collinearity_threshold)
     state%solution(state%used+1:) = 0d0
     !Adaptation can solve the problem again, so remember the coefficients

     if(state%verbosity > 0) write (*,*) "AAMIX(O:order):", state%order(:state%non_collinear)
     if(state%verbosity > 0) write (*,*) "AAMIX(S:solution):", state%solution(:state%used)

     if (state%adaptive_alpha .and. state%iteration >= state%adapt_from) then
       call adaptive_anderson_adapt_alpha(state)
     end if
     call adaptive_anderson_form_new_input(state, state%solution, x)

     !Set it here, since so far the value from the previous iteration have been used
     state%last_bii = state%solution(state%current)

     if(state%verbosity > 0) write (*,*) "AAMIX(D:iteration,bii,|residual|,alpha):", &
          state%iteration, state%solution(state%current), &
          sqrt(state%matrix(state%current, state%current)),state%alpha

     call adaptive_anderson_shift_circular_buffer(state)
     call dcopy(state%n, x, 1, state%inputs(1,state%current), 1)

     adaptive_anderson_step = .False.
     if (state%debug_store_to_file > 0) then
        write (state%debug_store_to_file,*) x
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
      integer :: increment
      integer :: j

      state%previous = state%current
      state%current = state%current + 1

      if (state%current > state%history) then
           state%current = 1
      end if

      !some residuum(s) have been removed due to lost of the orthogonality, shift the arrays
      if (state%used < state%history .and. state%used >= state%current) then
         state%used = state%used + 1
         do j=state%used, state%current+1, -1
              state%matrix(j,:) = state%matrix(j-1,:)
              state%matrix(:,j) = state%matrix(:,j-1)
              call dcopy(state%n, state%residuals(1,j-1), 1, state%residuals(1, j), 1)
              call dcopy(state%n, state%inputs(1, j-1), 1, state%inputs(1,j), 1)
         end do
      !just not all the history has been filled yet
      elseif (state%used < state%current) then
         state%used = state%current
      end if

      !TODO
      !if ((increment > 0) .and. (state%used > state%current) .and. (state%replace)):
      !replace the worst residuum instead of the last

      state%iteration = state%iteration + 1
  end subroutine

  subroutine adaptive_anderson_store_result(state, residuum)
      type(adaptive_anderson_solver_state), intent(inout) :: state
      real*8, intent(in), target :: residuum(state%n)
      real*8, target :: tmp(state%n)
      real*8, pointer :: normed_residuum(:)
      integer last, curr

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
      last = adaptive_anderson_previous_index(state)
      curr = state%current
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
      do i=2, state%non_collinear
         if(state%verbosity > 0) write (*,'(E16.8)', advance="no") state%qr_matrix(i,i)
         if ( abs(state%qr_matrix(i,i)) < mnorm ) then
            if(state%verbosity > 0)  write (*,*) 'AAMIX(COL):', state%qr_matrix(i,i), mnorm,  &
                                                                state%qr_matrix(1,1), threshold
            state%non_collinear = state%non_collinear - 1
            state%order(i:state%non_collinear) = state%order(i+1:state%non_collinear+1)
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

      !!!Order of vectors in householder orthogonalization
      do i=1,state%current
         state%order(i) = state%current - i + 1
      end do
      do i=1, state%used - state%current
         state%order(i+state%current) = state%used - i + 1
      end do
      state%non_collinear = state%used
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
      integer :: i,j

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
      logical :: ok

      integer :: i,j
      real*8 :: r

      j = 0
      if (coef > 1.2d0) then
        do i=2, state%non_collinear
            if (abs(solution(state%order(i))) >= abs(solution(state%order(1)))) then
               j = j + 1
               if (j >= 2) then
                  if (state%verbosity > 0) write (*,*) 'AAMIX(NA 5) - no adaptation', &
                                                 solution(state%order(i)), solution(state%order(j))
                  ok = .False.
                  return
               end if
            end if
        end do
      end if

      if ( coef < 1.2 ) then
         do i=2, state%non_collinear
               if (solution(state%order(i)) < -coef * 0.5) then
                  if (state%verbosity > 0) write (*,*) 'AAMIX(NA 6) - no adaptation', coef, solution(state%order(i))
                  ok = .False.
                  return
                end if
         end do
      end if
      ok = .True.
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
      integer*8 :: new_direction, res_direction
      logical :: direction_changed, aggressive
      integer*8 :: changed, i, j, k, oi, oj, vectors
      integer*8 :: directions
      real*8 :: histories(4)
      real*8 :: change_threshold
      real*8 :: solution(state%used)

      real*8 :: coefs
      integer*8 :: coefs_count

      logical :: ok

      integer :: tmp(1)

      histories = (/10,8,6,5/)
      solution = state%solution(:state%used)

      if (state%non_collinear <= 1) return

      if (state%adaptation_changed == 0 .and. state%matrix(state%current,state%current) > state%matrix(1,1) &
                                        .and. state%iteration < state%history .and. state%iteration < 3) then
        if (state%verbosity > 0) write (*,*) 'AAMIX(NA 1)', state%iteration,  &
                                             state%matrix(state%current,state%current), state%matrix(1,1)
        call adaptive_anderson_no_adaptation(state)
        return
      end if

      if (state%previous .ne. 0) then
        if ( &
          abs(state%previous_solution(state%previous)) / maxval(abs(state%previous_solution(:state%used) )) < 0.2 &
          .or. abs(state%previous_solution(state%previous)) < 0.2 ) then
         !if the  previous soulution is not present in the mix, the adaptation would not be related to the current alpha
          if (state%verbosity > 0) write (*,*) "AAMIX(NA 2)", state%matrix(state%current, state%current), &
                                    state%matrix(state%previous, state%previous)
          call adaptive_anderson_no_adaptation(state)
          return
        end if
      end if

      coef = abs(solution(state%current))
      ok = adaptive_anderson_adapt_alpha_check(state, solution, coef)
      if (.not. ok) then
        if (state%verbosity > 0) write (*,*) 'AAMIX(NA 3) - not ok(1)'
        call adaptive_anderson_no_adaptation(state)
        return
      end if

      desired = state%delta + (state%non_collinear - 1) * min(state%delta_per_vector,6d0)
      coef = coef / desired
      directions = merge(1,-1, coef > 0)


      if(state%verbosity > 0) write (*,*) "AAMIX(CO):", coef, desired

      if (coef > 1.D0) then
          new_direction = 1
      else if (coef < 1.D0) then
          if (coef < 0.01D0) coef = 0.01d0
          new_direction = -1
          coef = 1d0/coef
     else
          call adaptive_anderson_no_adaptation(state)
          if(state%verbosity > 0) write (*,*) 'AAMIX(NA 4) - coefficient == 1.0'
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

     if (coef > 1) then
         if (state%previous_solution(state%previous) > 0) then
            tmp = maxloc(state%previous_solution(:state%used))
         else
            tmp = minloc(state%previous_solution(:state%used))
         end if
         i = tmp(1)

         if (i .ne. state%previous .and. state%previous_solution(i) .ne. 0d0                                                   &
                                   .and. abs(state%previous_solution(state%previous) / state%previous_solution(i)) < 0.75 ) then
             if(state%verbosity > 0) write (*,*) "AAMIX(NA 3)", coef, i, state%previous_solution(i), &
                                                                state%previous_solution(state%previous)
             call adaptive_anderson_no_adaptation(state, coef)
             return
         end if

         if ( adaptive_anderson_convergence_ratio(state) < 0.2) then
            if(state%verbosity > 0) write (*,*) "AAMIX(NA 4)",state%matrix(state%current, state%current), &
                                                              state%matrix(state%previous, state%previous)
            call adaptive_anderson_no_adaptation(state, coef)
           return
         end if

     end if


     last_adapt_coef = state%last_adapt_coef
     state%last_adapt_coef = coef
     state%adaptation_direction = new_direction

     if (direction_changed .and. state%matrix(state%current, state%current) > &
                                  state%matrix(state%previous, state%previous)/10 ) then
       if (state%verbosity > 0) WRITE (*,*) 'AAMIX(CR) 3', coef, coef**(1d0/3d0)
       coef = coef**(1d0/3d0)
     else
       if (state%verbosity > 0) WRITE (*,*) 'AAMIX(CR) 2', coef, coef**0.5
       coef = coef**0.5
     end if

     if (coef .ne. 1d0 ) then
        state%last_alpha = state%alpha
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

        subroutine adaptive_anderson_no_adaptation(state, coef)
             type(adaptive_anderson_solver_state) :: state
             real*8, optional :: coef

             state%direction_changed = 0
             state%direction_unchanged = 0
             state%adaptation_delayed = state%adaptation_delayed+1
             if (present(coef)) state%last_adapt_coef = coef
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
       if( state%verbosity > 3) WRITE (*,*) 'AAMIX(STATE): Cleanup'
       deallocate(state)
       state => null()
  end subroutine

  !Unified handling of lapack errors
  subroutine check_lapack(info, msg)
    integer :: info
    character(len=*) :: msg
    if (info .ne. 0) then
       write (*,*) INFO
       write (*,*) msg
       stop 255
    end if
  end subroutine


end module
