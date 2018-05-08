module output
  implicit none

  character(len=255) :: outfile

contains

  subroutine dump_header
    use parameters, only : ngrid, ntime, outputdir
    use krome_user, only : krome_nmols, krome_nPhotoBins
    implicit none
    integer :: unit

    ! Write header for full history
    outfile = trim(outputdir)//"/output.dat"
    open(newunit=unit,file=trim(outfile), status="replace", action="write")
    write(unit,*) ngrid, ntime+1, krome_nmols, krome_nPhotoBins
    close(unit)
  end subroutine

  subroutine dump_snapshot(t)
    use parameters, only : ngrid, pc, spy
    use grid, only : r, n, nHtot, Tgas, Av, dr, tau
    implicit none
    real*8, intent(in) :: t
    integer :: i, unit

    open(newunit=unit,file=trim(outfile), access='append', status='old')
    do i=1,ngrid
      write(unit,'(200E17.8e3)') t/spy,r(i)/pc,Av(i),Tgas(i),n(:,i),tau(:,i)
    end do
    write(unit,*)
    close(unit)

  end subroutine

end module