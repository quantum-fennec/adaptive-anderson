      subroutine assert(val, msg)
        logical val
        character (len=* ) msg
        if (.not. val) then
          write (*,*) msg
          stop 255
        end if
      end subroutine

      subroutine assert_arr(a,b,msg)
        double precision :: a(4),b(4)
        character (len=* ) msg
        if ( .not. all(abs(a-b) < 1e-7) ) then
           write (*,*) msg
           write (*,*) a
           write (*,*) b
           stop 255
        end if
      end subroutine

      program test
        use adaptive_anderson_solver
        type(adaptive_anderson_solver_state), pointer :: state
        double precision :: x0(4)
        double precision :: x(4)

        double precision, parameter :: x1(4) = (/1.5,0.,0.,0. /)
        double precision, parameter :: x2(4) = (/2.0,0.,0.,0. /)
        double precision, parameter :: x3(4) = (/1.8,0.2,0.,0. /)

        logical :: res
        integer :: i

        x0 = 0.
        x0(1) = 1

        state => adaptive_anderson_init(4, x0, alpha = 0.5d0)
        state%inputs(:,2:)=1d100
        state%residuals(:,1:)=1d100
        res=adaptive_anderson_step(4,state, x0, x)
        call assert(.not. res, "Not converged yet")
        call assert(state%used .eq. 2, "used must be 2")
        call assert_arr(x,x1, "first x should be 1.5" )
        !1.0 0 0 0  => 1.0 0 0 0
        !1.5 0 0 0  => ?
        res=adaptive_anderson_step(4,state, x0, x)
        call assert(.not. res, "Not converged yet 2")
        call assert(state%used .eq. 3, "used must be 3")
        !the space for the result is allready "allocated"
        write(*,*) state%non_collinear
        call assert(state%non_collinear .eq. 1, &
                                "non_collinear must be 1")
        !1.5 0 0 0  => 1.0
        !2.0 0 0 0  => ?
        call assert_arr(x,x2, "first x should be 2" )
        x0(1) = -1.0
        x0(2) = 1

        state%adaptive_alpha = .false.
        res=adaptive_anderson_step(4,state, x0,x )
        !1.5 0 0 0 => 1.0 0 0 0
        !2.0 0 0 0 => 0.5 0.5 0 0
        !2.3 0.2 0 0
        do i=1,state%used
          WRITE (*,*) state%matrix(i, :state%used)
        end do
        WRITE (*,*) state%solution

        call assert_arr(x,x3, "assert x3" )
        call assert(.not. res, "Not yet converged 3")
        state%adaptive_alpha = .true.

        x0(1) = 0.3
        x0(2) = -0.2
        x0(3) = 0.2
        res = adaptive_anderson_step(4,state, x0,x )
        x0(3) = 0.1
        res = adaptive_anderson_step(4,state, x0,x )
        x0(2) = -0.1
        res = adaptive_anderson_step(4,state, x0,x )
        x0(2) = -0.01
        res = adaptive_anderson_step(4,state, x0,x )
        call adaptive_anderson_end(state)

        write (*,*) "Sanity checks"
        call assert(.not. associated(adaptive_anderson_init(4,x0, &
           threshold=0d0)), "Sanity check of the threshold arg. failed")
        call assert(.not. associated(adaptive_anderson_init(4,x0, &
            alpha=-1d0)), "Sanity check of the alpha arg. failed")
        call assert(.not. associated(adaptive_anderson_init(4,x0, &
             history=-1)), "Sanity check of the history arg. failed")
        x0(2) = -1
        call assert(.not. associated(adaptive_anderson_init(4,x0, &
             weights=x0)), "Sanity check of the weights arg. failed")
        write (*,*) "All tests are ok"
      end
