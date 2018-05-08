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
    use initcond, only : MC_in_ISM
    use output, only : dump_header
    use rt, only : init_rt
    implicit none

    call krome_init
    call read_parameters
    call init_grid
    call MC_in_ISM
    call init_rt
    call dump_header

  end subroutine

  subroutine run
    use parameters, only : tend, ntime
    use output, only : dump_snapshot
    use solver, only : step
    implicit none
    integer :: itime
    real*8 :: dt, t

    dt = tend/ntime
    t = 0d0
    do itime=1,ntime
      call dump_snapshot(t)
      call step(dt)
      t = t + dt
    end do
    call dump_snapshot(t)

  end subroutine

  subroutine cleanup
    use grid, only : cleanup_grid
    implicit none

    call cleanup_grid

  end subroutine
end program tinyPDR