program test_leiden_shielding
  use leiden_shielding
  implicit none

  type(leiden_shielding_table) :: table
  real*8 :: x, y, z
  real*8, allocatable :: zz(:,:)
  integer i, j, ispecies, unit

  ! table = read_leiden_shielding("shield.03.50.69-557-36.dat")

  ! ! test interpolation by interpolating to half-points
  ! ispecies = 1
  ! allocate(zz(table%nx, table%ny))
  ! do j=1,table%ny
  !   do i=1,table%nx
  !     x = table%x(i)
  !     y = table%y(j)
  !     if(i<table%nx) x = 0.5*(x+table%x(i+1))
  !     if(j<table%ny) y = 0.5*(y+table%y(j+1))
  !     zz(i,j) = interpolate(table, x, y, ispecies)
  !   end do
  ! end do

  ! ! Save reconstructed table to file
  ! open(newunit=unit,file="shield.03.50.69-557-36.interp",status="old")
  ! write(unit,*) zz(:,:)
  ! close(unit)

  call test_interpolate("../data/shield.03.50.69-557-36.dat")

end program test_leiden_shielding
