module parameters
  use krome_user, only : krome_seconds_per_year, krome_nPhotoBins
  implicit none

  real*8, parameter :: fav = 1d0/4d-22
  real*8, parameter :: pc = 3.085d18 !1 pc in cm
  real*8, parameter :: spy = krome_seconds_per_year !s

  integer :: ngrid    ! Number of grid points
  integer :: ntime    ! Number of time steps
  real*8  :: tend
  real*8  :: rmax     ! Simulation physical size
  real*8  :: rmin     ! Location of least grid point for logarithmic grids
  real*8  :: d2g      ! Dust to gas mass ratio
  character(len=255) :: inputfile  ! Input namelist file
  character(len=255) :: outputdir
  character(len=255) :: datadir
  real*8, dimension(krome_nPhotoBins+1) :: photobin_limits
  character(len=255) :: sedfile  ! File name for loading spectral energy distribution of incident radiation
  character(len=255) :: opacityfile  ! File name for loading opacities
  character(len=255) :: unitenergy  ! Energy unit in opacity file (e.g. "micron" - see Krome's documentation)
  real*8 :: crate                   ! Cosmic ray rate
  character(len=255) :: extinction_type
  character(len=255) :: initcond
  character(len=255) :: grid_type

contains

  subroutine read_parameters
    implicit none
    namelist/params/ngrid,tend,ntime,rmax,rmin,outputdir,datadir,d2g,photobin_limits,sedfile,opacityfile,unitenergy,crate,extinction_type,initcond,grid_type

    call init_to_defaults

    if(iargc() == 0) then
      print*, "Usage: tinyPDR <namelist-file>"
      stop
    else if(iargc() == 1) then
      call getarg(1,inputfile)
    end if
  
    print*, "Reading namelist 'params' from : ", trim(inputfile)
    open(1, file=trim(inputfile))
    read(1,nml=params)
    close(1)

    tend = tend*spy

    write(*,params)

  end subroutine

  subroutine init_to_defaults
    implicit none

    ngrid = 128
    ntime = 100
    rmin = 0d0
    rmax = 1d0*pc
    outputdir="."
    datadir="../../dat"
    tend = 1d0
    d2g = 1d-2
    crate = 2.5d-17

    photobin_limits = 0d0
    sedfile="black87_eV.interp"
    opacityfile="opacityDraineR35.dat"
    unitenergy="micron"
    extinction_type="single direction"
    initcond="uniform"
    grid_type="centered"

  end subroutine

end module