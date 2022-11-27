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
        type(adaptive_anderson_solver_state), pointer :: state
        integer, parameter :: n = 10
        double precision :: x(n)
        double precision :: residuum(n)

        x = 0
        state => adaptive_anderson_init(n, x, threshold=1d-10, history=5)
        do
          CALL SFC(n, x, residuum)
          IF ( adaptive_anderson_step(n, state, residuum, x)) EXIT
          WRITE (*,*) "Residual norm:", &
                      adaptive_anderson_residual_norm(state)
        end do
        WRITE (*,*) "Final residual norm:", &
                      adaptive_anderson_residual_norm(state)
        call adaptive_anderson_end(state)
      end
