module grid
  implicit none

  real*8, allocatable :: r(:)           ! Grid point positions
  real*8, allocatable :: n(:,:)         ! Number density of chemical species = n(ispecies, igrid)
  real*8, allocatable :: nHtot(:)       ! Total number of H nuclei
  real*8, allocatable :: Tgas(:)        ! Gas temperature
  real*8, allocatable :: Av(:)          ! Visual extinction
  real*8, allocatable :: tau(:,:)        ! Optical depth
  real*8              :: dr             ! Cell size
  real*8              :: rmax           ! Simulation physical size

contains

  subroutine init_grid
    use parameters, only : ngrid
    use krome_user, only : krome_nmols, krome_nPhotoBins
    implicit none

    allocate(r(ngrid))
    allocate(n(krome_nmols, ngrid))
    allocate(nHtot(ngrid))
    allocate(Tgas(ngrid))
    allocate(Av(ngrid))
    allocate(tau(krome_nPhotoBins,ngrid))

    ! Initialize arrays to zero
    r(:) = 0d0
    n(:,:) = 0d0
    nHtot(:) = 0d0
    Tgas(:) = 0d0
    Av(:) = 0d0
    tau(:,:) = 0d0

  end subroutine

  subroutine cleanup_grid
    implicit none

    deallocate(r, n, nHtot, Tgas, Av, tau)
  end subroutine

  subroutine centered_uniform
    use parameters, only : ngrid
    implicit none
    integer :: igrid

    ! Build a centered uniform grid
    dr = rmax/ngrid
    do igrid=1,ngrid
      r(igrid) = (igrid - 0.5_8)*dr
    end do
  end subroutine

end module