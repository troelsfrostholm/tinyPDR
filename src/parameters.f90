module parameters
  use krome_user, only : krome_seconds_per_year
  implicit none

  real*8, parameter :: fav = 1d0/4d-22
  real*8, parameter :: pc = 3.085d18 !1 pc in cm
  real*8, parameter :: spy = krome_seconds_per_year !s

  integer :: ngrid    ! Number of grid points
  integer :: ntime    ! Number of time steps
  real*8  :: tend
  character(len=255) :: inputfile  ! Input namelist file
  character(len=255) :: datadir

contains

  subroutine read_parameters
    implicit none
    namelist/params/ngrid,tend,ntime,datadir

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

    write(*,params)

  end subroutine

  subroutine init_to_defaults
    implicit none

    ngrid = 128
    ntime = 100
    datadir="."
    tend = 1d0

  end subroutine

end module