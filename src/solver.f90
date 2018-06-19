module solver
  implicit none

contains

  subroutine step(dt)
    use krome_user, only : krome_nmols, krome_nPhotoBins, krome_set_photoBinJ, krome_get_opacity_size_d2g, krome_idx_H2, krome_idx_CO, krome_set_user_gamma_H2, krome_set_user_gamma_CO
    use krome_main, only : krome
    use richtings_dissociation_rates, only : S_H2, S_H2_d, S_CO, S_CO_d, gamma_H2_thin, gamma_CO_thin
    use parameters, only : d2g, ngrid, extinction_type
    use grid, only : n, nHtot, Tgas, tau, dr
    use rt, only : j0
    use extinction, only : extinction_single_direction, extinction_spherical_cloud_uniform_incidence
    use util, only : cumsum
    implicit none
    real(kind=8), intent(in) :: dt
    real(kind=8), dimension(krome_nPhotoBins) :: dtau, dtau_prev, tau_max
    real(kind=8), dimension(krome_nPhotoBins,ngrid) :: mu, mu_spherical
    real(kind=8), dimension(krome_nPhotoBins,ngrid) :: j
    real(kind=8), dimension(krome_nPhotoBins) :: j_H2, j_CO
    real(kind=8), dimension(ngrid) :: N_H2, N_CO, N_Htot, T_col, gamma_H2, gamma_CO, mu_H2, mu_CO, tau_H2, tau_CO
    real(kind=8), dimension(krome_nmols) :: nn
    real(kind=8) :: N_H2_max, N_CO_max, N_Htot_max, T_col_max, tau_H2_0, tau_CO_0, Tgas_ss, mu_H2_0, mu_CO_0
    integer :: i, ibin

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
    select case (trim(extinction_type))
      case("single direction")
        do i=1,ngrid
          do ibin=1,krome_nPhotoBins
            mu(ibin,i) = 2d0*extinction_single_direction(tau(ibin,i), tau_max(ibin))
          end do
        end do
      case("spherical uniform")
        do i=1,ngrid
          do ibin=1,krome_nPhotoBins
            mu(ibin,i) = 2d0*extinction_spherical_cloud_uniform_incidence(tau(ibin,i), tau_max(ibin))
          end do
        end do
      case default
        print*, "Extinction type not supported: ", extinction_type
        stop
    end select

    ! Mean intensity (directional average)
    do i=1,ngrid
      j(:,i) = j0(:)*mu(:,i)
    end do

    ! Self-shielding species
    !------------------------

    ! Column densities - to cell centers
    N_H2 = (cumsum(n(krome_idx_H2,:),ngrid) - 0.5_8*n(krome_idx_H2,:))*dr
    N_CO = (cumsum(n(krome_idx_CO,:),ngrid) - 0.5_8*n(krome_idx_CO,:))*dr
    N_Htot = (cumsum(nHtot(:),ngrid) - 0.5_8*nHtot(:))*dr
    T_col = (cumsum(n(krome_idx_H2,:)*Tgas(:),ngrid) - 0.5_8*n(krome_idx_H2,:)*Tgas(:))*dr

    ! Column densities to rmax
    N_H2_max = maxval(N_H2(:)) + 0.5_8*n(krome_idx_H2,ngrid)*dr
    N_CO_max = maxval(N_CO(:)) + 0.5_8*n(krome_idx_CO,ngrid)*dr
    N_Htot_max = maxval(N_Htot(:)) + 0.5_8*nHtot(ngrid)*dr
    T_col_max = maxval(T_col(:)) + 0.5_8*n(krome_idx_H2,ngrid)*Tgas(ngrid)*dr
    
    ! Reverse direction of column densities, so they measure the density from the end of array (cloud outer surface)
    N_H2 = N_H2_max - N_H2(:)
    N_CO = N_CO_max - N_CO(:)
    N_Htot = N_Htot_max - N_Htot(:)
    T_col = T_col_max - T_col(:)
    
    ! Look up shielding factors
    do i=1,ngrid
      if(N_H2(i) < 1e-20) then
        Tgas_ss = Tgas(i)
      else
        Tgas_ss = T_col(i) / N_H2(i)
      endif
      mu_H2(i) = S_H2(N_H2(i), Tgas_ss)*S_H2_d(N_Htot(i))
      mu_CO(i) = S_CO(N_CO(i), N_Htot(i),Tgas_ss)*S_CO_d(N_Htot(i))
    end do

    mu_H2_0 = S_H2(N_H2_max, Tgas_ss)*S_H2_d(N_Htot_max)
    mu_CO_0 = S_CO(N_CO_max, N_Htot_max,Tgas_ss)*S_CO_d(N_Htot_max)

    if(trim(extinction_type) == "spherical uniform") then
      ! Go from extinction to optical depth
      tau_H2 = -log(mu_H2(:))
      tau_CO = -log(mu_CO(:))
      tau_H2_0 = -log(mu_H2_0)
      tau_CO_0 = -log(mu_CO_0)

      ! Reverse tau direction to measure from center of cloud
      tau_H2 = tau_H2_0 - tau_H2
      tau_CO = tau_CO_0 - tau_CO

      ! and back to extinction - now assuming a spherical cloud
      do i=1,ngrid
        mu_H2(i) = 2.0_8*extinction_spherical_cloud_uniform_incidence(tau_H2(i), tau_H2_0)
        mu_CO(i) = 2.0_8*extinction_spherical_cloud_uniform_incidence(tau_CO(i), tau_CO_0)
      end do
    endif

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