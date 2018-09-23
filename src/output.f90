module output
  implicit none

  character(len=255) :: outfile
  character(len=255) :: outfile_fluxes

contains

  subroutine dump_header
    use parameters, only : ngrid, ntime, outputdir
    use krome_user, only : krome_nmols, krome_nPhotoBins, krome_nrea, krome_get_cooling_array, krome_get_cooling_names_header
    implicit none
    integer :: unit

    ! Write header for full history
    outfile = trim(outputdir)//"/output.dat"
    open(newunit=unit,file=trim(outfile), status="replace", action="write")
    write(unit,*) ngrid, ntime+1, krome_nmols, krome_nPhotoBins
    close(unit)

    ! Write header for fluxes
    outfile_fluxes = trim(outputdir)//"/fluxes.dat"
    open(newunit=unit,file=trim(outfile_fluxes), status="replace", action="write")
    write(unit,*) ngrid, ntime, krome_nrea
    close(unit)
  end subroutine

  subroutine dump_snapshot(t)
    use parameters, only : ngrid, pc, spy
    use grid, only : r, n, nHtot, Tgas, Tdust, Av, dr, tau, fluxes
    implicit none
    real*8, intent(in) :: t
    integer :: i, unit, unit_heatcool

    open(newunit=unit,file=trim(outfile), access='append', status='old')
    do i=1,ngrid
      write(unit,'(200E17.8e3)') t/spy,r(i)/pc,Av(i),Tgas(i),Tdust(i),n(:,i),tau(:,i)
    end do
    write(unit,*)
    close(unit)

    open(newunit=unit,file=trim(outfile_fluxes), access='append', action="write")
    do i=1,ngrid
      write(unit,*) t/spy,r(i)/pc,Av(i),fluxes(:,i)
    end do
    write(unit,*)
    close(unit)
  end subroutine

  subroutine dump_cooling(t)
    use parameters, only : ngrid, pc, outputdir
    use krome_user, only : krome_nmols, krome_nPhotoBins, krome_nrea, krome_get_cooling_names_header
    use grid, only : r, n, nHtot, Tgas, Tdust, Av, dr, tau, cooling
    implicit none
    integer :: unit, i
    character(len=128) :: outfile
    real*8::t
    outfile = trim(outputdir)//"/cooling.dat"
    open(newunit=unit,file=trim(outfile), status='replace', action="write")
    write(unit,*) krome_get_cooling_names_header()
    do i=1,ngrid
      write(unit,'(200E17.8e3)') r(i)/pc,Av(i),nHtot(i),cooling(:,i)
    end do
    close(unit)
  end subroutine

  subroutine dump_heating(t)
    use parameters, only : ngrid, pc, outputdir
    use krome_user, only : krome_nmols, krome_nPhotoBins, krome_nrea, krome_get_heating_names_header
    use grid, only : r, n, nHtot, Tgas, Tdust, Av, dr, tau, heating
    implicit none
    integer :: unit, i
    character(len=128) :: outfile
    real*8::t
    outfile = trim(outputdir)//"/heating.dat"
    open(newunit=unit,file=trim(outfile), status='replace', action="write")
    write(unit,*) krome_get_heating_names_header()
    do i=1,ngrid
      write(unit,'(200E17.8e3)') r(i)/pc,Av(i),nHtot(i),heating(:,i)
    end do
    close(unit)
  end subroutine

end module