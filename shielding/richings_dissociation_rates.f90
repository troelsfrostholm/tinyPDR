  !
  ! Module for computing dissociation rates of H2 and CO
  ! following Richtings et. al. 2014 paper II (R14II in the following)
  !

module richtings_dissociation_rates
  real*8, parameter :: k_b = 1.38064852d-16         ! Boltzmann constant (erg K^-1)
  real*8, parameter :: m_H = 1.6737236d-24          ! Mass of hydrogen atom (g)
  real*8, parameter :: m_H2 = 2d0*m_H               ! Approx. mass of hydrogen molecule (g)
  real*8, parameter :: gamma_H2_thin = 7.5d-11      ! R14I, at photon number density 2.256 x 10^-4 cm^-3 btw. 12.24 and 13.51 eV
  real*8, parameter :: gamma_CO_thin = 2.6d-10      ! Visser 2009
  real*8, parameter :: fav = 4d-22                  ! N(H_tot) -> Av conversion factor from R14II eq. 3.4
  real*8, parameter :: gamma_H2_d = 3.74            ! Draine and Bertoldi 1996
  real*8, parameter :: gamma_CO_d = 3.53            ! R14II
  character(len=*),parameter :: data_dir="../../dat/"

  contains

  ! H2 self shielding factor
  function S_H2(N_H2, T)
    implicit none
    real*8, parameter :: one_over_one_point_three = 1d0/1.3d0
    
    real*8 :: N_H2         ! H2 column density
    real*8 :: T            ! Gas temperature
    real*8 :: S_H2         ! H2 self shielding extinction

    real*8 :: alpha, omega_H2, N_crit, x
    real*8 :: b_therm, b_turb, b_5

    ! Eq. 3.13 R14II
    omega_H2 = 0.013*(1d0 + (T/2700d0)**1.3)**one_over_one_point_three*exp(-(T/3900d0)**14.6)
    
    ! Eq. 3.14 - 3.15 R14II
    if (T < 3000d0) then
       alpha = 1.4d0
       N_crit = 1.3d0*(1d0 + (T/600d0)**0.8)
    else if (T < 4000d0) then
       alpha = (T/4500d0)**(-0.8)
       N_crit = (T/4760d0)**(-3.8)
    else
       alpha = 1.1d0
       N_crit = 2d0
    end if
    N_crit = N_crit * 1d14 ! cm^-2

    ! Eq. 3.16-3.17 R14II
    b_therm = sqrt(2d0*k_b*T/m_H2) ! cm s^-1
    b_turb = 7.1d5 ! cm s^-1
    b_5 = sqrt(b_therm**2 + b_turb**2)*1d-5

    ! print*, "omega_H2", omega_H2
    ! print*, "alpha", alpha
    ! print*, "N_crit", N_crit
    ! print*, "b_5", b_5

   ! omega_H2   1.3571033371816164E-002
   ! alpha   1.3999999999999999
   ! N_crit   204665392457850.22
   ! b_5   7.2721859922841201

    ! Overwrite everyting with values from Draine and Bertoldi
    ! omega_H2 = 0.035
    ! alpha = 2.0
    ! N_crit = 5d14 ! Reference number density of H2 (cm^-3)
    ! b_5 = 3d0

    x = N_H2 / N_crit        ! x' in the text
    
    ! Eq. 3.12 R14II
    S_H2 = (1d0 - omega_H2)/(1 + x/b_5)**alpha*exp(-5d-7*(1d0 + x)) &
         + omega_H2/(1d0 + x)**0.5*exp(-8.5d-4*(1d0 + x)**0.5)
  end function S_H2

  function S_H2_d(N_H_tot)
    implicit none
    real*8, intent(in) :: N_H_tot
    real*8 :: S_H2_d
    real*8 :: Av

    Av = fav*N_H_tot
    S_H2_d = exp(-gamma_H2_d*Av)
  end function S_H2_d

  function S_draine_bertoldi(N_H2)
    implicit none
    real*8             :: S_draine_bertoldi
    real*8, intent(in) :: N_H2
    real*8, parameter  :: omega_H = 0.035
    real*8, parameter  :: alpha = 2.0
    real*8, parameter  :: N_ref = 5e14 ! Reference number density of H2 (cm^-3)
    real*8, parameter  :: b = 3e5 ! Doppler broadening parameter (cm s^-1)
    real*8, parameter  :: b5 = b/1e5
    real*8             :: x

    x = N_H2/N_ref
    S_draine_bertoldi = (1 - omega_H)/(1+x/b5)**alpha + omega_H/(1+x)**0.5*exp(-8.5e-4*(1+x)**0.5)
  end function S_draine_bertoldi

  function gamma_H2(N_H2,N_H_tot,T)
    implicit none
    real*8, intent(in) :: N_H2, N_H_tot, T
    real*8 :: gamma_H2

    gamma_H2 = gamma_H2_thin*S_H2(N_H2,T)*S_H2_d(N_H_tot)
    !gamma_H2 = gamma_H2_thin*S_draine_bertoldi(N_H2)*S_H2_d(N_H_tot)
  end function gamma_H2

  ! CO self shielding factor, including shielding by H2
  ! For now ignores temperature and uses just one 
  ! CO extinction table
  function S_CO(N_CO,N_H2,T)
    use leiden_shielding
    implicit none
    real*8, intent(in) :: N_CO,  N_H2, T
    real*8 :: S_CO
    integer, parameter :: ispecies = 1  ! Use shielding values for 12C16O only
    character(len=*), parameter :: filename = data_dir//"shield.03.50.69-557-36.dat" ! at Tex(CO,H2) (K) = 50.00  353.55
    type(leiden_shielding_table), save :: table
    logical, save :: first_call = .true.

    if(first_call) then
      table = read_leiden_shielding(filename)
      first_call = .false.
    end if

    S_CO = interpolate(table, N_CO, N_H2, ispecies)
  end function

  function S_CO_d(N_H_tot)
    implicit none
    real*8, intent(in) :: N_H_tot
    real*8 :: S_CO_d
    real*8 :: Av

    Av = fav*N_H_tot
    S_CO_d = exp(-gamma_CO_d*Av)
  end function S_CO_d

  function gamma_CO(N_CO,N_H2,N_H_tot,T)
    implicit none
    real*8, intent(in) :: N_CO, N_H2, N_H_tot, T
    real*8 :: gamma_CO

    gamma_CO = gamma_CO_thin*S_CO(N_CO,N_H2,T)*S_CO_d(N_H_tot)
  end function gamma_CO
    
end module richtings_dissociation_rates
