module initcond

contains

  ! Set up initial condition
  subroutine initial_condition
    use parameters, only : initcond
    implicit none
    select case(trim(initcond))
      case ("MC_in_ISM")
        call MC_in_ISM
      case ("uniform")
        call uniform
      case default
        print*, "Unknown initial condition"
        stop
    end select
  end subroutine

  ! Power law molecular cloud profile with 
  ! constant low density outside and constant 
  ! high density inside
  ! This initial condition ignores rmin and overwrites rmax
  subroutine MC_in_ISM
    use krome_user
    use parameters, only : ngrid, inputfile, fav, rmax
    use grid, only : centered_uniform, n, nHtot, Tgas, Av, r, dr
    implicit none
    real*8, parameter :: T_floor = 10d0 ! temperature floor (K)
    real*8, parameter :: nT = 1d3       ! product of density and temperature
    real*8, parameter :: sigma = 1e-6
    real*8  :: Av_max
    real*8  :: r0
    real*8  :: r1
    real*8  :: n0
    real*8  :: n1
    real*8  :: alpha
    integer :: i
    namelist/initcond/Av_max,n0,n1

    ! Read initcond namelist
    Av_max = 1.0
    n0 = 1.0
    n1 = 2.0
    print*, "Reading namelist 'initcond' from : ", trim(inputfile)
    open(1, file=trim(inputfile))
    read(1,nml=initcond)
    close(1)
    write(*,initcond)

    ! Compute power law exponent and cloud size
    alpha = n0/(fav*Av_max*log(1d1))*(n1/n0 - 1d0)
    r1 = log10(n1/n0)/alpha
    r0 = 0.1*r1
    rmax = 1.25*(r0+r1)

    call centered_uniform

    ! Density and temperature
    do i=1,ngrid
      if(r(i) <= r0) then
        nHtot(i) = n1
      else if (r(i) >= r0+r1) then
        nHtot(i) = n0
      else
        nHtot(i) = n0*1d1**(alpha*(r0+r1-r(i)))
     end if
     Tgas(i) = max(nT/nHtot(i),T_floor)
    end do

    ! Compute visual extinction
    Av(1) = 0.5_8*nHtot(1)*dr / fav
    do i=2,ngrid
      Av(i) = Av(i-1) + 0.5_8*(nHtot(i-1)+nHtot(i))*dr / fav
    end do

    ! Composition
    n(:,:) = 0.0_8
    n(krome_idx_H,:) = nHtot(:)/3.
    n(krome_idx_H2,:) = nHtot(:)/3
    n(krome_idx_C,:) = 2.46d-4*nHtot(:) 
    n(krome_idx_O,:) = 4.9d-4*nHtot(:) 
    n(krome_idx_He,:) = 1d-1*nHtot(:)
    do i=1,ngrid
      if(r(i) > r0+r1) then
        n(krome_idx_H2,i) = 0.5*sigma*nHtot(i)
        n(krome_idx_H,i) = nHtot(i) - 2d0*n(krome_idx_H2,i)
      endif
      if(r(i) > r0+r1) then
        n(krome_idx_H,i) = sigma*nHtot(i)
        n(krome_idx_Hj,i) = nHtot(i) - (2d0*n(krome_idx_H2,i) + n(krome_idx_H,i))
        n(krome_idx_E,i) = n(krome_idx_Hj,i)
      endif
    end do
  end subroutine

  ! Uniform density, composition and temperature
  subroutine uniform
    use krome_user
    use util, only : assert
    use parameters, only : ngrid, inputfile, fav, rmax
    use grid, only : centered_uniform, n, nHtot, Tgas, Av, r, dr
    implicit none
    real*8 :: ntot, T, x(krome_nmols)
    integer :: i
    namelist/initcond/ntot,T,x

    ! Default values
    ntot = 1d0
    T = 10d0
    x(:) = 0d0

    ! Read namelist file
    print*, "Reading namelist 'initcond' from : ", trim(inputfile)
    open(1, file=trim(inputfile))
    read(1,nml=initcond)
    close(1)
    write(*,initcond)

    ! Validate input
    call assert(ntot > 0d0, "ntot must be > 0. It was ", ntot)
    call assert(T > 0d0, "T must be > 0. It was ", T)
    call assert(all(x(:) >= 0d0), "All x must be positive ")
    call assert(.not. all(x(:) == 0d0), "At least one x must be above 0")

    ! Normalize fractions
    x(:) = x(:) / sum(x(:))

    ! Compute number density of species
    do i=1,ngrid
      n(:,i) = x(:)*ntot
    end do

    ! Set gas temperature
    Tgas(:) = T

  end subroutine

end module