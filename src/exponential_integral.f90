module exponential_integral
  implicit none
  character(len=255), parameter :: table_file = "../../dat/exponential-integral.dat"
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
    dx = (xmax-xmin)/(N_value-1)
  end subroutine

  pure function expn(n, x)
    implicit none
    integer, intent(in) :: n
    real*8, intent(in) :: x
    real*8 :: expn
    real*8 :: x0, x1, y0, y1
    integer :: i

    if(x < xmin) then
      expn = table(n+1,1)
      return
    end if

    if(x >= xmax) then
      expn = table(n+1,N_value)
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