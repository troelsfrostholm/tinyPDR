module solver
  implicit none

contains

  subroutine step(dt)
    use krome_user, only : krome_nmols, krome_nPhotoBins, krome_set_photoBinJ, krome_get_opacity_size_d2g, krome_idx_H2, krome_idx_CO, krome_set_user_gamma_H2, krome_set_user_gamma_CO
    use krome_main, only : krome
    use richtings_dissociation_rates, only : S_H2, S_H2_d, S_CO, S_CO_d, gamma_H2_thin, gamma_CO_thin
    use parameters, only : d2g, ngrid
    use grid, only : n, nHtot, Tgas, tau, dr
    use rt, only : j0
    use util, only : cumsum
    implicit none
    real*8, intent(in) :: dt
    real*8, dimension(krome_nPhotoBins) :: dtau, dtau_prev, tau_max
    real*8, dimension(krome_nPhotoBins,ngrid) :: mu
    real*8, dimension(krome_nPhotoBins,ngrid) :: j
    real*8, dimension(krome_nPhotoBins) :: j_H2, j_CO
    real*8, dimension(ngrid) :: N_H2, N_CO, N_Htot, gamma_H2, gamma_CO, mu_H2, mu_CO
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

    ! Self-shielding species
    !------------------------

    ! Column densities - to cell centers
    N_H2 = (cumsum(n(krome_idx_H2,:),ngrid) - 0.5_8*n(krome_idx_H2,:))*dr
    N_CO = (cumsum(n(krome_idx_CO,:),ngrid) - 0.5_8*n(krome_idx_CO,:))*dr
    N_Htot = (cumsum(nHtot(:),ngrid)) - 0.5_8*nHtot(:)*dr

    ! Reverse direction of column densities, so they measure the density from the end of array (cloud outer surface)
    N_H2 = maxval(N_H2(:)) + 0.5_8*n(krome_idx_H2,ngrid)*dr - N_H2(:)
    N_CO = maxval(N_CO(:)) + 0.5_8*n(krome_idx_CO,ngrid)*dr - N_CO(:)
    N_Htot = maxval(N_Htot(:)) + 0.5_8*nHtot(:)*dr - N_Htot(:)
    
    ! Look up shielding factors
    do i=1,ngrid
      mu_H2(i) = S_H2(N_H2(i), Tgas(i))*S_H2_d(N_Htot(i))
      mu_CO(i) = S_CO(N_CO(i), N_Htot(i),Tgas(i))*S_CO_d(N_Htot(i))
    end do

    ! Compute attenuated H2 and CO dissociation rates
    gamma_H2(:) = mu_H2(:)*gamma_H2_thin
    gamma_CO(:) = mu_CO(:)*gamma_CO_thin

    ! ------------------------------------------------------------
    ! Chemistry
    ! ------------------------------------------------------------

    do i=1,ngrid
      ! Set flux in Krome
      call krome_set_photoBinJ(j(:,i))

      ! Set dissociation rates
      call krome_set_user_gamma_H2(gamma_H2(i))
      call krome_set_user_gamma_CO(gamma_CO(i))

      !call KROME to do chemistry
      nn = n(:,i)
      call krome(nn,Tgas(i),dt)
      n(:,i) = nn
    end do

  end subroutine
end module