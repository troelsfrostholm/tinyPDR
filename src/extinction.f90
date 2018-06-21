module extinction
  implicit none

contains
  
  subroutine init_extinction
    use exponential_integral, only : init
    implicit none

    call init

  end subroutine

  pure function extinction_single_direction_left(tau, tau_max)
    implicit none
    real*8, intent(in) :: tau, tau_max
    real*8 :: extinction_single_direction_left

    extinction_single_direction_left = 0.5d0*exp(-tau)
  end function

  pure function extinction_single_direction_right(tau, tau_max)
    implicit none
    real*8, intent(in) :: tau, tau_max
    real*8 :: extinction_single_direction_right

    extinction_single_direction_right = 0.5d0*exp(tau - tau_max)
  end function

  pure function extinction_spherical_cloud_uniform_incidence(tau, tau_max)
    use exponential_integral, only : expn
    implicit none
    real*8, intent(in) :: tau, tau_max
    real*8 :: extinction_spherical_cloud_uniform_incidence
    real*8 :: r
    real*8 :: x

    r = tau_max
    x = tau

    extinction_spherical_cloud_uniform_incidence = &
      & 1/(2d0*x)*( expn(3, r - x) - expn(3, r + x) + r*expn(2, r - x) - r*expn(2, r + x) )
  end function

end module