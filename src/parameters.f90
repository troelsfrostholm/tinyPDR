module parameters
  use krome_user, only : krome_seconds_per_year, krome_nPhotoBins
  implicit none

  real*8, parameter :: fav = 1d0/4d-22
  real*8, parameter :: pc = 3.085d18 !1 pc in cm
  real*8, parameter :: spy = krome_seconds_per_year !s
  real*8, parameter :: h = 4.13566766d-15 ! plancks constant in eV s
  real*8, parameter :: pi = 3.14159265359 ! pi
  real*8, parameter :: sigma = 5.6704d-5 ! Stefan-boltzmann constant (erg K-1 cm-2)
  integer, parameter :: nmaxprintfluxesfor = 30


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
  character(len=16)  :: print_fluxes_for(nmaxprintfluxesfor)
  integer :: print_fluxes_for_ids(nmaxprintfluxesfor), nprintfluxesfor

contains

  subroutine read_parameters
    implicit none
    namelist/params/ngrid,tend,ntime,rmax,rmin,outputdir,datadir,d2g,photobin_limits,sedfile,opacityfile,unitenergy,crate,extinction_type,initcond,grid_type,print_fluxes_for

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

    print_fluxes_for_ids = mol_names_to_ids(print_fluxes_for)
    print*, print_fluxes_for_ids
    print*, nprintfluxesfor

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
    print_fluxes_for(:) = ""

    photobin_limits = 0d0
    sedfile="black87_eV.interp"
    opacityfile="opacityDraineR35.dat"
    unitenergy="micron"
    extinction_type="single direction"
    initcond="uniform"
    grid_type="uniform"

  end subroutine

  ! Maps a list of molecule names to krome ids
  function mol_names_to_ids(names)
    use util, only : indexof
    use krome_user, only : krome_get_names, krome_nmols
    implicit none
    character*16 :: names(nmaxprintfluxesfor), all_names(krome_nmols), mol
    integer :: mol_names_to_ids(nmaxprintfluxesfor), i

    mol_names_to_ids(:) = 0
    all_names(:) = krome_get_names()
    do i=1,nmaxprintfluxesfor
      mol = trim(names(i))
      if(mol .ne. "") then
        mol_names_to_ids(i) = indexof(mol, all_names, krome_nmols)
        if(mol_names_to_ids(i) == -1) then  ! Exit with error if not found
          print*, "Error in read_parameters: Molecule in namelist parameter print_fluxes_for not found in network. "
          print*, "Molecule: ", mol
          print*, "Molecules in network: ", all_names
          stop
        end if
      else
        nprintfluxesfor = i-1
        exit
      end if
    end do
  end function

end module