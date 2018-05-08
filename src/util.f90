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
end module