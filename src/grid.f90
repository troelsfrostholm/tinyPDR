module grid
  implicit none

  real*8, allocatable :: r(:)           ! Grid point positions
  real*8, allocatable :: n(:,:)         ! Number density of chemical species = n(ispecies, igrid)
  real*8, allocatable :: nHtot(:)       ! Total number of H nuclei
  real*8, allocatable :: Tgas(:)        ! Gas temperature
  real*8, allocatable :: Tdust(:)       ! Dust temperature
  real*8, allocatable :: Av(:)          ! Visual extinction
  real*8, allocatable :: tau(:,:)       ! Optical depth
  real*8              :: dr             ! Cell size

contains

  subroutine init_grid
    use parameters, only : ngrid, grid_type
    use krome_user, only : krome_nmols, krome_nPhotoBins
    implicit none

    allocate(r(ngrid))
    allocate(n(krome_nmols, ngrid))
    allocate(nHtot(ngrid))
    allocate(Tgas(ngrid))
    allocate(Tdust(ngrid))
    allocate(Av(ngrid))
    allocate(tau(krome_nPhotoBins,ngrid))

    ! Initialize arrays to zero
    r(:) = 0d0
    n(:,:) = 0d0
    nHtot(:) = 0d0
    Tgas(:) = 0d0
    Tdust(:) = 0d0
    Av(:) = 0d0
    tau(:,:) = 0d0

    select case(trim(grid_type))
      case("centered")
        call centered_uniform
      case("log")
        call logarithmic
      case default
        print*, "Unknown grid type: ", grid_type
        stop
    end select

  end subroutine

  subroutine cleanup_grid
    implicit none

    deallocate(r, n, nHtot, Tgas, Av, tau)
  end subroutine

  ! Builds a centered uniform grid
  subroutine centered_uniform
    use parameters, only : ngrid, rmax
    implicit none
    integer :: igrid

    dr = rmax/ngrid
    do igrid=1,ngrid
      r(igrid) = (igrid - 0.5_8)*dr
    end do
  end subroutine

  ! Builds a logarithmic grid with first grid point at 
  ! rmin, and last at rmax
  subroutine logarithmic
    use parameters, only : ngrid, rmin, rmax
    use util, only : assert
    implicit none
    integer :: igrid

    call assert(rmin > 0d0, "rmin must be > 0 when using logarithmic grid. It was", rmin)

    do igrid=1,ngrid
      r(igrid) = 1d1**((igrid-1)*(log10(rmax)-log10(rmin))/(ngrid-1) &
          + log10(rmin))
    end do

  end subroutine

end module