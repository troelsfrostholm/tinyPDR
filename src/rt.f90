module rt
  use krome_user, only : krome_nPhotoBins
  implicit none

  real*8, dimension(krome_nPhotoBins) :: j0      ! Incident intensity (directional mean)

contains
  
  subroutine init_rt
    use krome_user, only : krome_set_photobinE_limits, krome_load_photoBin_file_2col, krome_load_opacity_table,krome_get_photoBinJ, krome_set_user_crate
    use parameters, only : photobin_limits, datadir, sedfile, opacityfile, unitenergy, crate
    use grid, only : Tgas
    implicit none

    call krome_set_photobinE_limits(photobin_limits, Tgas(1))
    call krome_load_photoBin_file_2col(trim(datadir)//"/"//trim(sedfile))
    call krome_load_opacity_table(trim(datadir)//"/"//trim(opacityfile),unitEnergy=trim(unitenergy))
    j0 = krome_get_photoBinJ()

    !set cosmic rays
    call krome_set_user_crate(crate)

  end subroutine
end module