module exponential_integral
  implicit none
  character(len=25), parameter :: table_file = "exponential-integral.dat"
  integer, parameter :: N_order = 3, N_value = 10000
  real*8, dimension(N_order+1,N_value) :: table
  real*8 :: dx, xmax, xmin

contains

  subroutine init
    implicit none
    integer, parameter :: unit=10
    integer i

    open(unit, file=trim(table_file), action='READ', form='FORMATTED')
    ! Skip header.. 
    read(unit,*)
    read(unit,*)
    ! Read table
    read(unit,*) table
    close(unit)

    xmax = maxval(table(1,:))
    xmin = minval(table(1,:))
    dx = xmax/N_value

    print*, xmin, xmax, dx
  end subroutine

  function expn(n, x)
    implicit none
    integer :: n
    real*8 :: x, expn
    real*8 :: x0, x1, y0, y1
    integer :: i

    if(x < xmin) then
      expn = xmin
      return
    end if

    if(x > xmax) then
      expn = xmax
      return
    end if

    ! Linear interpolation in table
    i = floor(x/dx) + 1
    x0 = table(1, i)
    x1 = table(1, i+1)
    y0 = table(n+1, i)
    y1 = table(n+1, i+1)

    expn = (y1 - y0)/(x1 - x0) * (x - x0) + y0

  end function
end module exponential_integral