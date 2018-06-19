module util
  implicit none

contains

  ! Returns an array with the cumulative sum of a
  function cumsum(a,n)
    implicit none
    real*8, intent(in) :: a(n)
    integer, intent(in) :: n
    real*8 :: cumsum(n)
    real*8 :: acc
    integer :: i

    acc = 0d0
    do i=1,n
      acc = acc + a(i)
      cumsum(i) = acc
    end do
  end function

  ! Prints msg and stops if expr is false
  subroutine assert(expr, msg, val)
    implicit none
    logical, intent(in) :: expr
    character(len=*), intent(in) :: msg
    real*8, optional, intent(in) :: val

    if(.not.expr) then
      if(present(val)) then
        print*, "Assertion failed: "//msg, val
      else
        print*, "Assertion failed: "//msg
      endif
      stop
    endif
  end subroutine
end module