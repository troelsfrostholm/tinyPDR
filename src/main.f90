program tinyPDR
  implicit none

  call init
  call run
  call cleanup

contains

  subroutine init
    use krome_main, only : krome_init
    use parameters, only : read_parameters
    use grid, only : init_grid
    use initcond, only : initial_condition
    use output, only : dump_header
    use extinction, only : init_extinction
    use rt, only : init_rt
    implicit none

    call krome_init
    call read_parameters
    call init_grid
    call initial_condition
    call init_extinction
    call init_rt
    call dump_header

  end subroutine

  subroutine run
    use parameters, only : tend, ntime
    use output
    use solver, only : step
    implicit none
    integer :: itime
    real*8 :: dt, t

    dt = tend/ntime
    t = 0d0
    do itime=1,ntime
      call dump_snapshot(t)
      call step(t, dt)
      t = t + dt
    end do
    call dump_snapshot(t)
#ifdef USE_HEATING
    call dump_heating(t)
#endif
#ifdef USE_COOLING
    call dump_cooling(t)
#endif

  end subroutine

  subroutine cleanup
    use grid, only : cleanup_grid
    implicit none

    call cleanup_grid

  end subroutine
end program tinyPDR