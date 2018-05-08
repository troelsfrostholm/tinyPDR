module initcond


contains

  subroutine MC_in_ISM
    use krome_user
    use parameters, only : ngrid, inputfile, fav
    use grid, only : centered_uniform, n, nHtot, Tgas, Av, rmax, r, dr
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
    r0 = 0.25*r1
    rmax = 1.1*(r0+r1)

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

end module