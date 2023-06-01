      !A WEIRD FUNCTION, WHOSE RESIDUUM SHOULD BE FOUND
      subroutine SFC(n, x, residuum)
        integer n
        double precision :: x(n)
        double precision :: residuum(n)
        integer i
        do i=1,n
           residuum(i) = 10 - x(i)**2 + x(mod(i+5, n)+1) + &
                         x(mod(i+3, n)+1)
        end do

         residuum(1) = residuum(1) + x(2)**2 / 100 -15
         residuum(3) = residuum(3) + sqrt(abs(x(4))) / 10 - 15
         residuum(7) = residuum(7) + sin(residuum(6))
      end subroutine

      program example
        use adaptive_anderson_solver
        !In this object, the inner state of the solver will be hold
        type(adaptive_anderson_solver_state), pointer :: state
        !The size of the problem (the size of input/residuum vector)
        integer, parameter :: n = 10
        !The current (on the beggining the initial) input vector
        double precision :: x(n)
        !The residuum corresponding to the input vector
        double precision :: residuum(n)

        !Set the initial input
        x = 0
        !Init the solver, set the parameters, etc...
        state => adaptive_anderson_init(n, x, tolerance=1d-10, history=5)

        !Run the selfconsisten cycle
        do
          !Compute residuum for the current input
          CALL SFC(n, x, residuum)
          !Supply the residuum to the solver, obtaining a new input vector.
          !If the convergence criterion is met, the function returns .TRUE..
          IF ( adaptive_anderson_step(state, residuum, x)) EXIT
          !Just some unnecessary degug output to see the convergence
          WRITE (*,*) "Residual norm:", &
                      adaptive_anderson_residual_norm(state)
        end do

        !Output the result
        WRITE (*,*) "The solution:", x
        WRITE (*,*) "Final residual norm:", &
                      adaptive_anderson_residual_norm(state)

        !Final cleanup (deallocation of the memory used by solver)
        call adaptive_anderson_end(state)
      end
