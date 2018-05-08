module solver
  implicit none

contains

  subroutine step(dt)
    use krome_user, only : krome_nmols, krome_nPhotoBins, krome_set_photoBinJ, krome_get_opacity_size_d2g
    use krome_main, only : krome
    use parameters, only : d2g, ngrid
    use grid, only : n, Tgas, tau, dr
    use rt, only : j0
    implicit none
    real*8, intent(in) :: dt
    real*8, dimension(krome_nPhotoBins) :: dtau, dtau_prev, tau_max
    real*8, dimension(krome_nPhotoBins,ngrid) :: mu
    real*8, dimension(krome_nPhotoBins,ngrid) :: j
    real*8, dimension(krome_nmols) :: nn
    integer :: i

    ! ------------------------------------------------------------
    ! Radiative transfer
    ! ------------------------------------------------------------

    ! Optical depth
    dtau_prev(:) = krome_get_opacity_size_d2g(n(:,1),Tgas(1),dr,d2g)
    tau(:,1) = 0.5_8*dtau_prev(:)
    do i=2,ngrid
      dtau(:) = krome_get_opacity_size_d2g(n(:,i),Tgas(i),dr,d2g)
      tau(:,i) = tau(:,i-1) + 0.5_8*(dtau_prev + dtau)
      dtau_prev(:) = dtau(:)
    end do
    do i=1,krome_nPhotoBins
      tau_max(i) = maxval(tau(i,:)) + 0.5_8*dtau(i)  ! tau(rmax)
    end do
    
    ! Total extinction (cm^-1)
    do i=1,ngrid
      mu(:,i) = exp(tau(:,i) - tau_max(:))
    end do

    ! Mean intensity (directional average)
    do i=1,ngrid
      j(:,i) = j0(:)*mu(:,i)
    end do

    ! ------------------------------------------------------------
    ! Chemistry
    ! ------------------------------------------------------------

    do i=1,ngrid
      ! Set flux in Krome
      call krome_set_photoBinJ(j(:,i))

      ! Set dissociation rates
      ! call krome_set_user_gamma_H2(gamma_thick(ibin_H2,i))
      ! call krome_set_user_gamma_CO(gamma_thick(ibin_CO,i))

      !call KROME to do chemistry
      nn = n(:,i)
      call krome(nn,Tgas(i),dt)
      n(:,i) = nn
    end do

  end subroutine
end module