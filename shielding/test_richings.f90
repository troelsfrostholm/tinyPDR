program test_richtings
  use richtings_dissociation_rates

  integer, parameter :: n=24
  real*8, parameter :: T=3d2, f_Av = 4d-22
  integer i, unit
  real*8 :: N_H2, S, Av

  open(newunit=unit,file="richtings_S_H2.txt")
  do i=1,n
    N_H2 = 1d1**i
    S = S_H2(N_H2,T)
    Av = f_Av*N_H2
    write(unit,*) N_H2, Av, S
  end do
  close(unit)

end program