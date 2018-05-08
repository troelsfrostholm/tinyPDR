program test_extinction
  use exponential_integral, only : init
  use extinction, only : extinction_spherical_cloud_uniform_incidence
  use exponential_integral, only : expn
  implicit none

  integer, parameter :: N=10000
  real*8, parameter :: r=1d2, dr=r/N
  real*8 :: x
  integer :: i

  call init

  do i=1,N
    x = i*dr
    write(1,*) x, extinction_spherical_cloud_uniform_incidence(x, r)
  end do

end program