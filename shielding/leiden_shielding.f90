! This module provides functions to load shielding functions as provided on 
! http://home.strw.leidenuniv.nl/~ewine/photo/index.php?file=N2_photodissociation.php

module leiden_shielding

  type leiden_shielding_table
    integer :: nx, ny, nspecies
    real*8 :: logxmul, logymul                    ! Inverse logarithmic spacing
    real*8, allocatable :: x(:), y(:), table(:,:,:)
    real*8, allocatable :: logx(:), logy(:), logtable(:,:,:)
    character(len=20), allocatable :: species_names(:)
  end type

  contains

  ! Reads dimension metadata from row string. E.g. number of table columns. 
  function read_dimension(row_string)
    implicit none
    character(len=*), intent(in) :: row_string
    integer :: read_dimension

    integer :: ichar
    character(len=10) :: tail

    ! Read table dimensions
    ! first search for the '=' sign
    do ichar=1,len(row_string)
      if(row_string(ichar:ichar) == "=") exit
    end do
    if(ichar == len(row_string)) then
      print *,"ERROR: in read_leiden_shielding failed to read table dimensions"
      stop
    end if
    tail = trim(row_string(ichar+1:len(row_string)))
    read(tail,'(i4)') read_dimension
  end function read_dimension

  function format_string_exp(n)
    ! Construct formatting string
    ! To match n elements like : 1.000E+00
    implicit none
    integer, intent(in) :: n
    character(len=13) :: format_string_exp

    write(format_string_exp, '(a,i0.4,a)') '(', n, 'E10.3E3)'
  end function format_string_exp


  ! Reads leiden shielding functions from file
  function read_leiden_shielding(filename)
    implicit none
    character(len=*), intent(in) :: filename
    type(leiden_shielding_table) :: read_leiden_shielding

    integer, parameter :: row_string_len = 102
    type(leiden_shielding_table) :: t
    integer :: unit, ios, i, ix, iy, nmeta, nspec, nelem
    character(len=row_string_len)::row_string
    character(len=13) :: fmt

    ! Open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in read_leiden_shielding file ",trim(filename)," not found!"
      stop
    end if

    ! Skip metadata until specification of table dimensions
    nmeta = 0
    do
      read(unit,'(a)') row_string
      nmeta = nmeta + 1
      if(row_string(1:4)=="n[N(") exit
    end do

    ! Read table dimensions
    t%nx = read_dimension(row_string)
    read(unit,'(a)') row_string
    nmeta = nmeta + 1
    t%ny = read_dimension(row_string)

    ! Scan remaining lines and infer number of species
    nspec = -2
    do
      read(unit,'(a)',iostat=ios) row_string
      if(ios > 0) exit
      
      if(.not. row_string(3:3) == ".") then
        nspec = nspec + 1
      end if
    end do

    t%nspecies = nspec

    ! Allocate data structures
    allocate(t%x(t%nx))
    allocate(t%y(t%ny))
    allocate(t%species_names(nspec))
    allocate(t%table(t%nx, &
             t%ny, &
             nspec))

    ! Rewind and skip metadata again
    rewind(unit)
    do i=1,nmeta
      read(unit, *)
    end do

    ! Read x-axis
    fmt = format_string_exp(1)
    read(unit,*) ! Skip axis header
    do i=1,t%nx
      read(unit,fmt) t%x(i)
    end do

    ! Read y-axis
    fmt = format_string_exp(1)
    read(unit,*) ! Skip axis header
    do i=1,t%ny
      read(unit,fmt) t%y(i)
    end do

    ! Read tables
    do i=1,nspec
      read(unit,'(a)') row_string
      t%species_names(i) = trim(row_string)
      do iy=1,t%ny
        do ix=1,t%nx,10
          nelem = min(10,t%nx - ix + 1)
          fmt = format_string_exp(nelem)
          read(unit,fmt) t%table(ix:ix+nelem-1,iy,i)
        end do
      end do
    end do

    ! Store inverse logarithmic spacing
    t%logxmul = 1d0/(log(t%x(3)) - log(t%x(2)))
    t%logymul = 1d0/(log(t%y(3)) - log(t%y(2)))

    ! Store logarithmic table for fitting
    allocate(t%logx(t%nx))
    allocate(t%logy(t%ny))
    allocate(t%logtable(t%nx,t%ny,t%nspecies))
    t%logx(:) = log(t%x(:))
    t%logy(:) = log(t%y(:))
    t%logtable(:,:,:) = log(t%table(:,:,:))

    read_leiden_shielding = t

  end function read_leiden_shielding

  subroutine dump_table(table)
    implicit none
    type(leiden_shielding_table), intent(in) :: table
    integer :: ispec, irow
    character(len=13) :: fmt

    print*, "nx", table%nx
    print*, "ny", table%ny
    print*, "nspecies", table%nspecies
    print*
    print*, "x", table%x(:)
    print*, "y", table%y(:)
    do ispec=1,table%nspecies
      print*
      print*, table%species_names(ispec)
      print*, "==================="
      do irow=1,table%ny
        print*, table%table(:,irow,ispec)
      end do
    end do
  end subroutine dump_table

  ! Performs linear interpolation in loglog space
  function interpolate(table, xx, yy, ispecies)
    use krome_fit
    implicit none

    type(leiden_shielding_table), intent(in) :: table
    real*8, intent(in) :: xx, yy
    integer, intent(in) :: ispecies
    real*8 :: interpolate
    real*8 :: xmul, ymul
    real*8 :: x, y
    integer nx, ny

    x = xx
    y = yy

    nx = table%nx
    ny = table%ny

    ! Spacing differs at the table lower boundary and inside the table
    ! so there are the following edge and corner cases
    if(x >= table%x(2) .and. y >= table%y(2)) then
      ! The table interior has uniform log spacing
      interpolate = fit_anytab2D(table%logx(2:nx), &
                                 table%logy(2:ny), &
                                 table%logtable(2:nx,2:ny,ispecies), &
                                 table%logxmul, &
                                 table%logymul, &
                                 log(x), log(y))
    else if(x < table%x(2) .and. y < table%y(2)) then
      ! Optically thin corner has different spacing on both axes
      x = max(x, table%x(1))
      y = max(y, table%y(1))
      xmul = 1d0/(log(table%x(2)) - log(table%x(1)))
      ymul = 1d0/(log(table%y(2)) - log(table%y(1)))
      interpolate = fit_anytab2D(table%logx(:), &
                                 table%logy(:), &
                                 table%logtable(1:nx,1:ny,ispecies), &
                                 xmul, &
                                 ymul, &
                                 log(x), log(y))
    else if(x < table%x(2)) then
      ! Edge w. x optically thin has different spacing in x
      x = max(x, table%x(1))
      xmul = 1d0/(log(table%x(2)) - log(table%x(1)))
      interpolate = fit_anytab2D(table%logx(:), &
                                 table%logy(2:ny), &
                                 table%logtable(:,2:ny,ispecies), &
                                 xmul, &
                                 table%logymul, &
                                 log(x), log(y))
    else if(y < table%y(2)) then
      ! Edge w. y optically thin has different spacing in y
      y = max(y, table%y(1))
      ymul = 1d0/(log(table%y(2)) - log(table%y(1)))
      interpolate = fit_anytab2D(table%logx(2:nx), &
                                 table%logy(:), &
                                 table%logtable(2:nx,:,ispecies), &
                                 table%logxmul, &
                                 ymul, &
                                 log(x), log(y))
    end if

    interpolate = exp(interpolate)
  end function interpolate

  subroutine test_interpolate(filename)
    implicit none
    character(len=*) :: filename
    integer :: unit1, unit2
    type(leiden_shielding_table) :: table
    real*8 :: x, y
    integer :: i, j, ispecies
    character(len=10), parameter :: fmt = "(E10.3E2)"
    
    open(newunit=unit1,file=filename//".org",status="replace")
    open(newunit=unit2,file=filename//".fit",status="replace")

    table = read_leiden_shielding(filename)
    ispecies = 1
    do j=1,table%ny
      do i=1,table%nx
        x = table%x(i)
        y = table%y(j)
        write(unit1,*) x, y, table%table(i,j,ispecies)
        write(unit2,*) x, y, interpolate(table, x, y, ispecies)
      end do
    end do

    close(unit1)
    close(unit2)
  end subroutine test_interpolate

end module leiden_shielding