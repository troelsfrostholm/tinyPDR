
!############### MODULE ##############
module krome_commons
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-10-24 12:49:09
  !  Changeset b21f657
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************
  integer,parameter::idx_E=1
  integer,parameter::idx_Hk=2
  integer,parameter::idx_Ck=3
  integer,parameter::idx_Ok=4
  integer,parameter::idx_H=5
  integer,parameter::idx_HE=6
  integer,parameter::idx_H2=7
  integer,parameter::idx_C=8
  integer,parameter::idx_O=9
  integer,parameter::idx_OH=10
  integer,parameter::idx_CO=11
  integer,parameter::idx_CH=12
  integer,parameter::idx_CH2=13
  integer,parameter::idx_C2=14
  integer,parameter::idx_HCO=15
  integer,parameter::idx_H2O=16
  integer,parameter::idx_O2=17
  integer,parameter::idx_CO_total=18
  integer,parameter::idx_H2O_total=19
  integer,parameter::idx_Hj=20
  integer,parameter::idx_HEj=21
  integer,parameter::idx_H2j=22
  integer,parameter::idx_Cj=23
  integer,parameter::idx_Oj=24
  integer,parameter::idx_HOCj=25
  integer,parameter::idx_HCOj=26
  integer,parameter::idx_H3j=27
  integer,parameter::idx_CHj=28
  integer,parameter::idx_CH2j=29
  integer,parameter::idx_COj=30
  integer,parameter::idx_CH3j=31
  integer,parameter::idx_OHj=32
  integer,parameter::idx_H2Oj=33
  integer,parameter::idx_H3Oj=34
  integer,parameter::idx_O2j=35
  integer,parameter::idx_HEjj=36
  integer,parameter::idx_CR=37
  integer,parameter::idx_g=38
  integer,parameter::idx_Tgas=39
  integer,parameter::idx_dummy=40
  integer,parameter::nrea=281
  integer,parameter::nmols=36
  integer,parameter::nspec=40
  integer,parameter::natoms=5
  integer,parameter::ndust=0
  integer,parameter::ndustTypes=0
  integer,parameter::nPhotoBins=7
  integer,parameter::nPhotoRea=17

  !cooling index
  integer,parameter::idx_cool_h2 = 1
  integer,parameter::idx_cool_h2gp = 2
  integer,parameter::idx_cool_atomic = 3
  integer,parameter::idx_cool_cen = 3
  integer,parameter::idx_cool_hd = 4
  integer,parameter::idx_cool_metal = 5
  integer,parameter::idx_cool_z = 5
  integer,parameter::idx_cool_dh = 6
  integer,parameter::idx_cool_enthalpic = 6
  integer,parameter::idx_cool_dust = 7
  integer,parameter::idx_cool_compton = 8
  integer,parameter::idx_cool_cie = 9
  integer,parameter::idx_cool_cont = 10
  integer,parameter::idx_cool_continuum = 10
  integer,parameter::idx_cool_expansion = 11
  integer,parameter::idx_cool_exp = 11
  integer,parameter::idx_cool_ff = 12
  integer,parameter::idx_cool_bss = 12
  integer,parameter::idx_cool_custom = 13
  integer,parameter::idx_cool_co = 14
  integer,parameter::idx_cool_zcie = 15
  integer,parameter::idx_cool_zcienouv = 16
  integer,parameter::idx_cool_zextend = 17
  integer,parameter::idx_cool_gh = 18
  integer,parameter::ncools = 18

  !heating index
  integer,parameter::idx_heat_chem = 1
  integer,parameter::idx_heat_compress = 2
  integer,parameter::idx_heat_compr = 2
  integer,parameter::idx_heat_photo = 3
  integer,parameter::idx_heat_dh = 4
  integer,parameter::idx_heat_enthalpic = 4
  integer,parameter::idx_heat_av = 5
  integer,parameter::idx_heat_photoav = 5
  integer,parameter::idx_heat_cr = 6
  integer,parameter::idx_heat_dust = 7
  integer,parameter::idx_heat_xray = 8
  integer,parameter::idx_heat_viscous = 9
  integer,parameter::idx_heat_visc = 9
  integer,parameter::idx_heat_custom = 10
  integer,parameter::idx_heat_zcie = 11
  integer,parameter::nheats = 11

  real*8::arr_k(nrea)

  !commons for rate tables
  !modify ktab_n according to the required precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx

  !thermo toggle (when >0 do cooling/heating)
  integer::krome_thermo_toggle
  !$omp threadprivate(krome_thermo_toggle)

  !debug bit flag, print and array with fallback values for extreme environments
  integer:: red_flag
  real*8::n_global(nspec)
  integer, save :: nprint_negative=10
  !$omp threadprivate(n_global,nprint_negative,red_flag)

  !commons for implicit RHS
  integer::arr_r1(nrea)
  integer::arr_r2(nrea)
  integer::arr_r3(nrea)
  integer::arr_p1(nrea)
  integer::arr_p2(nrea)
  integer::arr_p3(nrea)

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

  !commons for frequency bins
  real*8::photoBinJ(nPhotoBins) !intensity per bin, eV/sr/cm2
  real*8::photoBinJ_org(nPhotoBins) !intensity per bin stored, eV/sr/cm2
  real*8::photoBinEleft(nPhotoBins) !left limit of the freq bin, eV
  real*8::photoBinEright(nPhotoBins) !right limit of the freq bin, eV
  real*8::photoBinEmid(nPhotoBins) !middle point of the freq bin, eV
  real*8::photoBinEdelta(nPhotoBins) !size of the freq bin, eV
  real*8::photoBinEidelta(nPhotoBins) !inverse of the size of the freq bin, 1/eV
  real*8::photoBinJTab(nPhotoRea,nPhotoBins) !xsecs table, cm2
  real*8::photoBinRates(nPhotoRea) !photo rates, 1/s
  real*8::photoBinHeats(nPhotoRea) !photo heating, erg/s
  real*8::photoBinEth(nPhotoRea) !energy treshold, eV
  real*8::photoPartners(nPhotoRea) !index of the photoreactants
  real*8::opacityDust(nPhotoBins) !interpolated opacity from tables
  !$omp threadprivate(photoBinJ,photoBinJ_org,photoBinEleft,photoBinEright,photoBinEmid, &
      !$omp    photoBinEdelta,photoBinEidelta,photoBinJTab,photoBinRates,photoBinHeats,photoBinEth, &
      !$omp    photoPartners)

  ! Draine dust absorption data loaded from file, via load_kabs
  ! in krome_photo module
  real*8::find_Av_draine_kabs(nPhotoBins)

  !commons for H2 photodissociation (Solomon)
  ! note: paramters here are set depending on the data
  ! but if you have a different file you should modify them
  integer,parameter::H2pdData_nvibX=15
  integer,parameter::H2pdData_nvibB=37
  real*8::H2pdData_dE(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_pre(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_EX(H2pdData_nvibX)
  integer::H2pdData_binMap(H2pdData_nvibX,H2pdData_nvibB)

  !commons for dust optical properties

  !square of turbulence velocity for broadening
  real*8::broadeningVturb2

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0, krome_omp_thread
  !$omp threadprivate(krome_omp_thread)

  !user-defined commons variables from the reaction file
  real*8::user_crate,user_Av,user_G0,user_gamma_CO,user_gamma_H2
  !$omp threadprivate(user_crate,user_Av,user_G0,user_gamma_CO,user_gamma_H2)

  !commons for anytab

  !physical commons
  real*8::phys_Tcmb
  real*8::phys_zredshift
  real*8::phys_orthoParaRatio
  real*8::phys_metallicity
  real*8::phys_Tfloor
  !$omp threadprivate(phys_Tcmb)
  !$omp threadprivate(phys_zredshift)
  !$omp threadprivate(phys_orthoParaRatio)
  !$omp threadprivate(phys_metallicity)
  !$omp threadprivate(phys_Tfloor)

  !machine precision
  real*8::krome_epsilon

  !xrayJ21 for tabulated heating and rate
  real*8::J21xray

  !total metallicity relative to solar Z/Z_solar
  real*8::total_Z
  real*8::dust2gas_ratio

  !commons for dust tabs (cool,H2,Tdust)
  integer,parameter::dust_tab_imax=50, dust_tab_jmax=50
  real*8::dust_tab_ngas(dust_tab_imax)
  real*8::dust_tab_Tgas(dust_tab_jmax)
  real*8::dust_mult_Tgas,dust_mult_ngas
  real*8::dust_table_AvVariable_log

  real*8::dust_tab_cool(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_heat(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_Tdust(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_H2(dust_tab_imax, dust_tab_jmax)
  !$omp threadprivate(dust_tab_cool,dust_tab_Tdust,dust_tab_H2)

  !commons for exp(-a) table
  integer,parameter::exp_table_na=int(1d5)
  real*8,parameter::exp_table_aMax=1d4,exp_table_aMin=0d0
  real*8,parameter::exp_table_multa=(exp_table_na-1) &
      / (exp_table_aMax-exp_table_aMin)
  real*8,parameter::exp_table_da=1d0/exp_table_multa
  real*8::exp_table(exp_table_na)

  !stores the last evaluation of the rates in the fex
  real*8::last_coe(nrea)
  !$omp threadprivate(last_coe)

  !data for CO cooling
  integer,parameter::coolCOn1=40
  integer,parameter::coolCOn2=40
  integer,parameter::coolCOn3=40
  real*8::coolCOx1(coolCOn1),coolCOx2(coolCOn2),coolCOx3(coolCOn3)
  real*8::coolCOixd1(coolCOn1-1),coolCOixd2(coolCOn2-1),coolCOixd3(coolCOn3-1)
  real*8::coolCOy(coolCOn1,coolCOn2,coolCOn3)
  real*8::coolCOx1min,coolCOx1max
  real*8::coolCOx2min,coolCOx2max
  real*8::coolCOx3min,coolCOx3max
  real*8::coolCOdvn1,coolCOdvn2,coolCOdvn3

  !xsecs from file variables
  !xsec for C -> C+ + E
  real*8,allocatable::xsec217_val(:)
  real*8::xsec217_Emin
  real*8::xsec217_idE
  integer::xsec217_n

  !xsec for H2 -> H2+ + E
  real*8,allocatable::xsec218_val(:)
  real*8::xsec218_Emin
  real*8::xsec218_idE
  integer::xsec218_n

  !xsec for H- -> H + E
  real*8,allocatable::xsec219_val(:)
  real*8::xsec219_Emin
  real*8::xsec219_idE
  integer::xsec219_n

  !xsec for CH -> C + H
  real*8,allocatable::xsec220_val(:)
  real*8::xsec220_Emin
  real*8::xsec220_idE
  integer::xsec220_n

  !xsec for CH -> CH+ + E
  real*8,allocatable::xsec221_val(:)
  real*8::xsec221_Emin
  real*8::xsec221_idE
  integer::xsec221_n

  !xsec for C2 -> C + C
  real*8,allocatable::xsec222_val(:)
  real*8::xsec222_Emin
  real*8::xsec222_idE
  integer::xsec222_n

  !xsec for OH -> O + H
  real*8,allocatable::xsec223_val(:)
  real*8::xsec223_Emin
  real*8::xsec223_idE
  integer::xsec223_n

  !xsec for OH -> OH+ + E
  real*8,allocatable::xsec224_val(:)
  real*8::xsec224_Emin
  real*8::xsec224_idE
  integer::xsec224_n

  !xsec for H2O -> OH + H
  real*8,allocatable::xsec225_val(:)
  real*8::xsec225_Emin
  real*8::xsec225_idE
  integer::xsec225_n

  !xsec for H2O -> H2O+ + E
  real*8,allocatable::xsec226_val(:)
  real*8::xsec226_Emin
  real*8::xsec226_idE
  integer::xsec226_n

  !xsec for O2 -> O2+ + E
  real*8,allocatable::xsec227_val(:)
  real*8::xsec227_Emin
  real*8::xsec227_idE
  integer::xsec227_n

  !xsec for O2 -> O + O
  real*8,allocatable::xsec228_val(:)
  real*8::xsec228_Emin
  real*8::xsec228_idE
  integer::xsec228_n

  !xsec for H2 -> H+ + H + E
  real*8,allocatable::xsec229_val(:)
  real*8::xsec229_Emin
  real*8::xsec229_idE
  integer::xsec229_n

  ! Gibbs free energy data from file variables

  !partition function from file
  integer,parameter::zpart_nCO=641
  integer,parameter::zpart_nH2even=2000
  integer,parameter::zpart_nH2odd=2000
  real*8::zpart_CO(zpart_nCO),minpart_CO,partdT_CO
  real*8::zpart_H2even(zpart_nH2even),minpart_H2even,partdT_H2even
  real*8::zpart_H2odd(zpart_nH2odd),minpart_H2odd,partdT_H2odd

  !Habing flux for the photoelectric heating by dust
  ! and clumping factor for H2 formation
  ! on dust by Jura/Gnedin
  real*8::GHabing,Ghabing_thin,clump_factor
  !$omp threadprivate(GHabing,GHabing_thin)

  ! Photo reaction rates relevant for Gnedin-Hollon cooling/heating function
  real*8::PLW,PHI,PHeI,PCVI
  !$omp threadprivate(PLW,PHI,PHeI,PCVI)

  !partition functions common vars

  !verbatim reactions
  character*50::reactionNames(nrea)

  !store evaluate once rate
  real*8::rateEvaluateOnce(nrea)
  !$omp threadprivate(rateEvaluateOnce)

end module krome_commons

!############### MODULE ##############
module krome_constants
  implicit none

  !constants
  real*8,parameter::boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::iboltzmann_eV = 1d0/boltzmann_eV !K / eV
  real*8,parameter::iboltzmann_erg = 1d0/boltzmann_erg !K / erg
  real*8,parameter::planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::planck_J = 6.62606957d-34 !J s
  real*8,parameter::planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::iplanck_eV = 1d0/planck_eV !1 / eV / s
  real*8,parameter::iplanck_J = 1d0/planck_J !1 / J / s
  real*8,parameter::iplanck_erg = 1d0/planck_erg !1 / erg / s
  real*8,parameter::gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::e_mass = 9.10938188d-28 !g
  real*8,parameter::p_mass = 1.67262158d-24 !g
  real*8,parameter::n_mass = 1.674920d-24 !g
  real*8,parameter::ip_mass = 1d0/p_mass !1/g
  real*8,parameter::clight = 2.99792458e10 !cm/s
  real*8,parameter::pi = 3.14159265359d0 !#
  real*8,parameter::eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::km_to_cm = 1d5 !km -> cm
  real*8,parameter::cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::kvgas_erg = 8.d0*boltzmann_erg/pi/p_mass !
  real*8,parameter::pre_kvgas_sqrt = sqrt(8.d0*boltzmann_erg/pi) !
  real*8,parameter::pre_planck = 2.d0*planck_erg/clight**2 !erg/cm2*s3
  real*8,parameter::exp_planck = planck_erg / boltzmann_erg !s*K
  real*8,parameter::stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::N_avogadro = 6.0221d23 !#
  real*8,parameter::Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::hubble = 0.704d0 !dimensionless
  real*8,parameter::Omega0 = 1.0d0 !dimensionless
  real*8,parameter::Omegab = 0.0456d0 !dimensionless
  real*8,parameter::Hubble0 = 1.d2*hubble*km_to_cm*cm_to_Mpc !1/s

end module krome_constants

!############### MODULE ##############
module krome_fit
contains

  !*****************************
  subroutine init_anytab3D(filename,x,y,z,f,xmul,ymul,zmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8::rout(4)
    integer::i,j,k,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(f,1)) then
      print *,"ERROR: in init_anytab3D x size differs from f(x,y,z)"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(f,2)) then
      print *,"ERROR: in init_anytab3D y size differs from f(x,y,z)"
      stop
    end if

    !check the size of the Z input array
    if(size(z).ne.size(f,3)) then
      print *,"ERROR: in init_anytab3D z size differs from f(x,y,z)"
      stop
    end if

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab3D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of grid points"
      print *," per dimension in the format"
      print *,"  XX, YY, ZZ"
      print *,row_string
      stop
    end if

    !loop to read file (3rd dimension of f() is
    ! first in the tables. i.e. tables are z,x,y,
    ! while f() is x,y,z
    do i=1,size(z)
      do j=1,size(x)
        do k=1,size(y)
          read(unit,*,iostat=ios) rout(:)
          y(k) = rout(3)
          f(j,k,i) = rout(4)
        end do
        x(j) = rout(2)
        read(unit,*,iostat=ios) !skip blanks
      end do
      z(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))
    zmul = 1d0/(z(2)-z(1))

  end subroutine init_anytab3D

  !********************************************
  !load 2d tables from filename
  subroutine init_anytab2D(filename,x,y,z,xmul,ymul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:,:),xmul,ymul
    real*8::rout(3)
    integer::i,j,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(z,1)) then
      print *,"ERROR: in init_anytab2D x size differs from z"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(z,2)) then
      print *,"ERROR: in init_anytab2D y size differs from z"
      stop
    end if

    if (krome_mpi_rank<=1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab2D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of rows and "
      print *," columns in the format"
      print *,"  RR, CC"
      print *,row_string
      stop
    end if

    !loop to read file
    do i=1,size(x)
      do j=1,size(y)
        read(unit,*,iostat=ios) rout(:)
        y(j) = rout(2)
        z(i,j) = rout(3)
      end do
      x(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))

  end subroutine init_anytab2D

  !********************************************
  !load 1d tables from filename
  subroutine init_anytab1D(filename,x,y,xmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),xmul
    real*8::rout(2)
    integer::i,ios,unit

    !check the size of the X input array
    if(size(x) /= size(y)) then
      print *,"ERROR: in init_anytab1D x size differs from y"
      stop
    end if

    if (krome_mpi_rank <= 1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios /= 0) then
      print *,"ERROR: in init_anytab1D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    ! !check if first line is OK
    ! if(scan(row_string,",")==0) then
    !    print *,"ERROR: file "//filename//" should"
    !    print *," contain the number of rows and "
    !    print *," columns in the format"
    !    print *,"  RR, CC"
    !    print *,row_string
    !    stop
    ! end if

    !loop to read file
    do i=1,size(x)
      read(unit,*,iostat=ios) rout(:)
      y(i) = rout(2)
      x(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios /= 0) exit

    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))

  end subroutine init_anytab1D

  !******************************
  !test 2d fit and save to file
  subroutine test_anytab2D(fname,x,y,z,xmul,ymul)
    implicit none
    integer::i,j,unit1,unit2
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul
    real*8::xx,yy,zz
    character(len=*),intent(in)::fname

    open(newunit=unit1,file=fname//".fit",status="replace")
    open(newunit=unit2,file=fname//".org",status="replace")
    do i=1,size(x)
      do j=1,size(y)
        xx = x(i)
        yy = y(j)
        zz = fit_anytab2D(x(:),y(:),z(:,:),xmul,ymul,xx,yy)
        write(unit1,*) xx,yy,zz
        write(unit2,*) x(i),y(j),z(i,j)
      end do
      write(unit1,*)
      write(unit2,*)
    end do
    close(unit1)
    close(unit2)
    print *,"original file wrote in ",fname//".org"
    print *,"fit test file wrote in ",fname//".fit"

  end subroutine test_anytab2D

  !*****************************
  function fit_anytab3D(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D

  !******************************
  !return 2d fit at xx,yy
  function fit_anytab2D(x,y,z,xmul,ymul,xx,yy)
    implicit none
    real*8::fit_anytab2D
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D(x(:),zright(:),xmul,xx)

    fit_anytab2D = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D

  !*********************
  !return 1d fit at xx
  function fit_anytab1D(x,z,xmul,xx)
    real*8,intent(in)::x(:),z(:),xmul,xx
    real*8::fit_anytab1D,p
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    fit_anytab1D = p * (z(i2) - z(i1)) + z(i1)

  end function fit_anytab1D

  !*****************************
  function fit_anytab3D_linlinlog(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D_linlinlog,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D_linlog(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D_linlog(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D_linlinlog = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D_linlinlog

  !***************************
  function fit_anytab2D_linlog(x,y,z,xmul,ymul,xx,yy)
    real*8::fit_anytab2D_linlog,x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D_linlog(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D_linlog(x(:),zright(:),xmul,xx)

    fit_anytab2D_linlog = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D_linlog

  !*********************
  function fit_anytab1D_linlog(x,z,xmul,xx)
    real*8::fit_anytab1D_linlog,x(:),z(:),xmul,xx,p,z2,z1
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    z2 = z(i2)
    z1 = z(i1)
    if(z1<0d0 .and. z2<0d0) then
      z1 = log10(-z1)
      z2 = log10(-z2)
      fit_anytab1D_linlog = -1d1**(p * (z2 - z1) + z1)
      return
    end if

    if(z1>0d0 .and. z2>0d0) then
      z1 = log10(z1)
      z2 = log10(z2)
      fit_anytab1D_linlog = 1d1**(p * (z2 - z1) + z1)
      return
    end if

    fit_anytab1D_linlog = (p * (z2 - z1) + z1)

  end function fit_anytab1D_linlog

  !*****************************
  !spline interpolation at t using array  x,y (size n) as data
  function fspline(x,y,t)
    implicit none
    real*8::fspline,x(:),y(:),b(size(x)),c(size(x)),d(size(x)),t
    integer::n

    n = size(x)
    call spline(x(:),y(:),b(:),c(:),d(:),n)
    fspline = ispline(t,x(:),y(:),b(:),c(:),d(:),n)

  end function fspline

  !*******************************+
  subroutine spline(x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer::n
    real*8::x(n), y(n), b(n), c(n), d(n)
    integer::i, j, gap
    real*8::h

    gap = n-1

    !check input
    if(n<2) return
    if(n<3) then
      b(1) = (y(2)-y(1))/(x(2)-x(1)) !linear interpolation
      c(1) = 0d0
      d(1) = 0d0
      b(2) = b(1)
      c(2) = 0d0
      d(2) = 0d0
      return
    end if

    !step 1: preparation
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2d0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do

    ! step 2: end conditions
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0d0
    c(n) = 0d0
    if(n.ne.3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

    ! step 3: forward elimination
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do

    ! step 4: back substitution
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

    ! step 5: compute spline coefficients
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2d0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3d0*c(i)
    end do
    c(n) = 3d0*c(n)
    d(n) = d(n-1)
  end subroutine spline

  !*******************************
  function ispline(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    real*8::ispline
    integer::n
    real*8::u, x(n), y(n), b(n), c(n), d(n)
    integer::i, j, k
    real*8::dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u<=x(1)) then
      ispline = y(1)
      return
    end if

    if(u>=x(n)) then
      ispline = y(n)
      return
    end if

    ! binary search for for i, such that x(i) <= u <= x(i+1)
    i = 1
    j = n+1
    do while (j>i+1)
      k = (i+j)/2
      if(u<x(k)) then
        j=k
      else
        i=k
      end if
    end do

    ! evaluate spline interpolation
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

  end function ispline

end module krome_fit
!This module contains useful routines to get physical
! quantities, like mean molecular weight, mass density,
! mass, jeans length, etc. etc.

!############### MODULE ##############
module krome_getphys
contains

  !*****************************
  !get the mean molecular weight
  function get_mu(n)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:),get_mu,m(nspec)
    m(:) = get_mass()

    !ip_mass is 1/proton_mass_in_g
    get_mu = max(sum(n(1:nmols)*m(1:nmols)),1d-40) &
        / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu

  !***************************
  !get mean molecular weight
  function get_mu_rho(n,rhogas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_mu_rho,rhogas,n(:)

    !ip_mass is 1/proton_mass_in_g
    get_mu_rho = rhogas / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu_rho

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

    get_mass(1) = 9.10938188d-28	!E
    get_mass(2) = 1.67444345638d-24	!H-
    get_mass(3) = 2.00771060473d-23	!C-
    get_mass(4) = 2.67691710837d-23	!O-
    get_mass(5) = 1.67353251819d-24	!H
    get_mass(6) = 6.69206503638d-24	!HE
    get_mass(7) = 3.34706503638d-24	!H2
    get_mass(8) = 2.00761951091d-23	!C
    get_mass(9) = 2.67682601455d-23	!O
    get_mass(10) = 2.84417926637d-23	!OH
    get_mass(11) = 4.68444552546d-23	!CO
    get_mass(12) = 2.17497276273d-23	!CH
    get_mass(13) = 2.34232601455d-23	!CH2
    get_mass(14) = 4.01523902183d-23	!C2
    get_mass(15) = 4.85179877728d-23	!HCO
    get_mass(16) = 3.01153251819d-23	!H2O
    get_mass(17) = 5.3536520291d-23	!O2
    get_mass(18) = 0.d0	!CO_total
    get_mass(19) = 0.d0	!H2O_total
    get_mass(20) = 1.67262158d-24	!H+
    get_mass(21) = 6.69115409819d-24	!HE+
    get_mass(22) = 3.34615409819d-24	!H2+
    get_mass(23) = 2.00752841709d-23	!C+
    get_mass(24) = 2.67673492073d-23	!O+
    get_mass(25) = 4.85170768346d-23	!HOC+
    get_mass(26) = 4.85170768346d-23	!HCO+
    get_mass(27) = 5.01968661638d-24	!H3+
    get_mass(28) = 2.17488166891d-23	!CH+
    get_mass(29) = 2.34223492073d-23	!CH2+
    get_mass(30) = 4.68435443164d-23	!CO+
    get_mass(31) = 2.50958817255d-23	!CH3+
    get_mass(32) = 2.84408817255d-23	!OH+
    get_mass(33) = 3.01144142437d-23	!H2O+
    get_mass(34) = 3.17879467619d-23	!H3O+
    get_mass(35) = 5.35356093528d-23	!O2+
    get_mass(36) = 6.69024316d-24	!HE++
    get_mass(37) = 0.d0	!CR
    get_mass(38) = 0.d0	!g
    get_mass(39) = 0.d0	!Tgas
    get_mass(40) = 0.d0	!dummy

  end function get_mass

  !************************
  !get sqrt of the inverse of the masses (1/sqrt(g))
  function get_imass_sqrt()
    use krome_commons
    implicit none
    real*8::get_imass_sqrt(nspec)

    get_imass_sqrt(1) = 3.31326021505d+13	!E
    get_imass_sqrt(2) = 7.72795806394d+11	!H-
    get_imass_sqrt(3) = 2.23177004181d+11	!C-
    get_imass_sqrt(4) = 1.93278051341d+11	!O-
    get_imass_sqrt(5) = 7.73006102111d+11	!H
    get_imass_sqrt(6) = 3.86562679981d+11	!HE
    get_imass_sqrt(7) = 5.46597856701d+11	!H2
    get_imass_sqrt(8) = 2.23182067346d+11	!C
    get_imass_sqrt(9) = 1.93281339991d+11	!O
    get_imass_sqrt(10) = 1.87508740611d+11	!OH
    get_imass_sqrt(11) = 1.46106959624d+11	!CO
    get_imass_sqrt(12) = 2.14423849574d+11	!CH
    get_imass_sqrt(13) = 2.06621889668d+11	!CH2
    get_imass_sqrt(14) = 1.57813553259d+11	!C2
    get_imass_sqrt(15) = 1.43565011358d+11	!HCO
    get_imass_sqrt(16) = 1.82224271009d+11	!H2O
    get_imass_sqrt(17) = 1.36670546184d+11	!O2
    get_imass_sqrt(18) = 0.d0	!CO_total
    get_imass_sqrt(19) = 0.d0	!H2O_total
    get_imass_sqrt(20) = 7.732165696d+11	!H+
    get_imass_sqrt(21) = 3.86588992536d+11	!HE+
    get_imass_sqrt(22) = 5.46672253003d+11	!H2+
    get_imass_sqrt(23) = 2.23187130855d+11	!C+
    get_imass_sqrt(24) = 1.93284628808d+11	!O+
    get_imass_sqrt(25) = 1.43566359113d+11	!HOC+
    get_imass_sqrt(26) = 1.43566359113d+11	!HCO+
    get_imass_sqrt(27) = 4.463357746d+11	!H3+
    get_imass_sqrt(28) = 2.14428340044d+11	!CH+
    get_imass_sqrt(29) = 2.06625907582d+11	!CH2+
    get_imass_sqrt(30) = 1.46108380244d+11	!CO+
    get_imass_sqrt(31) = 1.99617572781d+11	!CH3+
    get_imass_sqrt(32) = 1.87511743463d+11	!OH+
    get_imass_sqrt(33) = 1.82227027061d+11	!H2O+
    get_imass_sqrt(34) = 1.77365342346d+11	!H3O+
    get_imass_sqrt(35) = 1.36671708942d+11	!O2+
    get_imass_sqrt(36) = 3.86615310465d+11	!HE++
    get_imass_sqrt(37) = 0.d0	!CR
    get_imass_sqrt(38) = 0.d0	!g
    get_imass_sqrt(39) = 0.d0	!Tgas
    get_imass_sqrt(40) = 0.d0	!dummy

  end function get_imass_sqrt

  !************************
  !get inverse of the species masses (1/g)
  function get_imass()
    use krome_commons
    implicit none
    real*8::get_imass(nspec)

    get_imass(1) = 1.09776932527d+27	!E
    get_imass(2) = 5.9721335838d+23	!H-
    get_imass(3) = 4.98079751954d+22	!C-
    get_imass(4) = 3.73564051301d+22	!O-
    get_imass(5) = 5.97538433901d+23	!H
    get_imass(6) = 1.49430705554d+23	!HE
    get_imass(7) = 2.9876921695d+23	!H2
    get_imass(8) = 4.98102351847d+22	!C
    get_imass(9) = 3.73576763885d+22	!O
    get_imass(10) = 3.51595278056d+22	!OH
    get_imass(11) = 2.13472436506d+22	!CO
    get_imass(12) = 4.59775872662d+22	!CH
    get_imass(13) = 4.26926052901d+22	!CH2
    get_imass(14) = 2.49051175924d+22	!C2
    get_imass(15) = 2.06109124864d+22	!HCO
    get_imass(16) = 3.32056849448d+22	!H2O
    get_imass(17) = 1.86788381943d+22	!O2
    get_imass(18) = 0.d0	!CO_total
    get_imass(19) = 0.d0	!H2O_total
    get_imass(20) = 5.97863863505d+23	!H+
    get_imass(21) = 1.4945104915d+23	!HE+
    get_imass(22) = 2.98850552203d+23	!H2+
    get_imass(23) = 4.98124953791d+22	!C+
    get_imass(24) = 3.73589477335d+22	!O+
    get_imass(25) = 2.0611299469d+22	!HOC+
    get_imass(26) = 2.0611299469d+22	!HCO+
    get_imass(27) = 1.99215623688d+23	!H3+
    get_imass(28) = 4.59795130141d+22	!CH+
    get_imass(29) = 4.2694265684d+22	!CH2+
    get_imass(30) = 2.13476587776d+22	!CO+
    get_imass(31) = 3.98471753628d+22	!CH3+
    get_imass(32) = 3.51606539365d+22	!OH+
    get_imass(33) = 3.32066893916d+22	!H2O+
    get_imass(34) = 3.14584646656d+22	!H3O+
    get_imass(35) = 1.86791560251d+22	!O2+
    get_imass(36) = 1.49471398286d+23	!HE++
    get_imass(37) = 0.d0	!CR
    get_imass(38) = 0.d0	!g
    get_imass(39) = 0.d0	!Tgas
    get_imass(40) = 0.d0	!dummy

  end function get_imass

  !************************
  !species binding energies (surface=BARE), K
  function get_EbindBare()
    use krome_commons
    implicit none
    real*8::get_EbindBare(nspec)

    get_EbindBare(:) = 1d99

    get_EbindBare(idx_H) = 500.0d0
    get_EbindBare(idx_H2) = 300.0d0
    get_EbindBare(idx_O) = 1700.0d0
    get_EbindBare(idx_OH) = 1360.0d0
    get_EbindBare(idx_CO) = 1100.0d0
    get_EbindBare(idx_HCO) = 1100.0d0
    get_EbindBare(idx_H2O) = 4800.0d0
    get_EbindBare(idx_O2) = 1250.0d0

  end function get_EbindBare

  !************************
  !species binding energies (surface=ICE), K
  function get_EbindIce()
    use krome_commons
    implicit none
    real*8::get_EbindIce(nspec)

    get_EbindIce(:) = 1d99

    get_EbindIce(idx_H) = 650.0d0
    get_EbindIce(idx_H2) = 300.0d0
    get_EbindIce(idx_O) = 1700.0d0
    get_EbindIce(idx_OH) = 3500.0d0
    get_EbindIce(idx_CO) = 1300.0d0
    get_EbindIce(idx_HCO) = 3100.0d0
    get_EbindIce(idx_H2O) = 4800.0d0
    get_EbindIce(idx_O2) = 900.0d0

  end function get_EbindIce

  !************************
  function get_kevap70()
    use krome_commons
    implicit none
    real*8::get_kevap70(nspec)

    get_kevap70(idx_E) = 0d0
    get_kevap70(idx_Hk) = 0d0
    get_kevap70(idx_Ck) = 0d0
    get_kevap70(idx_Ok) = 0d0
    get_kevap70(idx_H) = 790490323.12
    get_kevap70(idx_HE) = 0d0
    get_kevap70(idx_H2) = 13763786733.1
    get_kevap70(idx_C) = 0d0
    get_kevap70(idx_O) = 28.3692788833
    get_kevap70(idx_OH) = 3649.88043081
    get_kevap70(idx_CO) = 149751.929641
    get_kevap70(idx_CH) = 0d0
    get_kevap70(idx_CH2) = 0d0
    get_kevap70(idx_C2) = 0d0
    get_kevap70(idx_HCO) = 149751.929641
    get_kevap70(idx_H2O) = 1.65884938156e-18
    get_kevap70(idx_O2) = 17568.7715065
    get_kevap70(idx_CO_total) = 0d0
    get_kevap70(idx_H2O_total) = 0d0
    get_kevap70(idx_Hj) = 0d0
    get_kevap70(idx_HEj) = 0d0
    get_kevap70(idx_H2j) = 0d0
    get_kevap70(idx_Cj) = 0d0
    get_kevap70(idx_Oj) = 0d0
    get_kevap70(idx_HOCj) = 0d0
    get_kevap70(idx_HCOj) = 0d0
    get_kevap70(idx_H3j) = 0d0
    get_kevap70(idx_CHj) = 0d0
    get_kevap70(idx_CH2j) = 0d0
    get_kevap70(idx_COj) = 0d0
    get_kevap70(idx_CH3j) = 0d0
    get_kevap70(idx_OHj) = 0d0
    get_kevap70(idx_H2Oj) = 0d0
    get_kevap70(idx_H3Oj) = 0d0
    get_kevap70(idx_O2j) = 0d0
    get_kevap70(idx_HEjj) = 0d0
    get_kevap70(idx_CR) = 0d0
    get_kevap70(idx_g) = 0d0
    get_kevap70(idx_Tgas) = 0d0
    get_kevap70(idx_dummy) = 0d0

  end function get_kevap70

  !************************
  !get verbatim reaction names
  function get_rnames()
    use krome_commons
    implicit none
    character*50::get_rnames(nrea)

    !reaction names are loaded from file
    get_rnames(:) = reactionNames(:)

  end function get_rnames

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

    get_names(1) = "E"
    get_names(2) = "H-"
    get_names(3) = "C-"
    get_names(4) = "O-"
    get_names(5) = "H"
    get_names(6) = "HE"
    get_names(7) = "H2"
    get_names(8) = "C"
    get_names(9) = "O"
    get_names(10) = "OH"
    get_names(11) = "CO"
    get_names(12) = "CH"
    get_names(13) = "CH2"
    get_names(14) = "C2"
    get_names(15) = "HCO"
    get_names(16) = "H2O"
    get_names(17) = "O2"
    get_names(18) = "CO_total"
    get_names(19) = "H2O_total"
    get_names(20) = "H+"
    get_names(21) = "HE+"
    get_names(22) = "H2+"
    get_names(23) = "C+"
    get_names(24) = "O+"
    get_names(25) = "HOC+"
    get_names(26) = "HCO+"
    get_names(27) = "H3+"
    get_names(28) = "CH+"
    get_names(29) = "CH2+"
    get_names(30) = "CO+"
    get_names(31) = "CH3+"
    get_names(32) = "OH+"
    get_names(33) = "H2O+"
    get_names(34) = "H3O+"
    get_names(35) = "O2+"
    get_names(36) = "HE++"
    get_names(37) = "CR"
    get_names(38) = "g"
    get_names(39) = "Tgas"
    get_names(40) = "dummy"

  end function get_names

  !************************
  !get cooling names list (empty element if cooling not present)
  function get_cooling_names()
    use krome_commons
    implicit none
    character*16::get_cooling_names(ncools)

    get_cooling_names(:) = ""

    get_cooling_names(idx_cool_h2) = "H2"
    get_cooling_names(idx_cool_h2gp) = "H2GP"
    get_cooling_names(idx_cool_atomic) = "ATOMIC"
    get_cooling_names(idx_cool_cen) = "CEN"
    get_cooling_names(idx_cool_hd) = "HD"
    get_cooling_names(idx_cool_metal) = "METAL"
    get_cooling_names(idx_cool_z) = "Z"
    get_cooling_names(idx_cool_dh) = "DH"
    get_cooling_names(idx_cool_enthalpic) = "ENTHALPIC"
    get_cooling_names(idx_cool_dust) = "DUST"
    get_cooling_names(idx_cool_compton) = "COMPTON"
    get_cooling_names(idx_cool_cie) = "CIE"
    get_cooling_names(idx_cool_cont) = "CONT"
    get_cooling_names(idx_cool_continuum) = "CONTINUUM"
    get_cooling_names(idx_cool_expansion) = "EXPANSION"
    get_cooling_names(idx_cool_exp) = "EXP"
    get_cooling_names(idx_cool_ff) = "FF"
    get_cooling_names(idx_cool_bss) = "BSS"
    get_cooling_names(idx_cool_custom) = "CUSTOM"
    get_cooling_names(idx_cool_co) = "CO"
    get_cooling_names(idx_cool_zcie) = "ZCIE"
    get_cooling_names(idx_cool_zcienouv) = "ZCIENOUV"
    get_cooling_names(idx_cool_zextend) = "ZEXTEND"
    get_cooling_names(idx_cool_gh) = "GH"

  end function get_cooling_names

  !************************
  !get heating names list (empty element if heating not present)
  function get_heating_names()
    use krome_commons
    implicit none
    character*16::get_heating_names(nheats)

    get_heating_names(:) = ""

    get_heating_names(idx_heat_chem) = "CHEM"
    get_heating_names(idx_heat_compress) = "COMPRESS"
    get_heating_names(idx_heat_compr) = "COMPR"
    get_heating_names(idx_heat_photo) = "PHOTO"
    get_heating_names(idx_heat_dh) = "DH"
    get_heating_names(idx_heat_enthalpic) = "ENTHALPIC"
    get_heating_names(idx_heat_av) = "AV"
    get_heating_names(idx_heat_photoav) = "PHOTOAV"
    get_heating_names(idx_heat_cr) = "CR"
    get_heating_names(idx_heat_dust) = "DUST"
    get_heating_names(idx_heat_xray) = "XRAY"
    get_heating_names(idx_heat_viscous) = "VISCOUS"
    get_heating_names(idx_heat_visc) = "VISC"
    get_heating_names(idx_heat_custom) = "CUSTOM"
    get_heating_names(idx_heat_zcie) = "ZCIE"

  end function get_heating_names

  !******************************
  !get the total number of H nuclei
  function get_Hnuclei(n)
    use krome_commons
    real*8::n(:),get_Hnuclei,nH

    nH = n(idx_Hk) + &
        n(idx_H) + &
        n(idx_H2)*2d0 + &
        n(idx_OH) + &
        n(idx_CH) + &
        n(idx_CH2)*2d0 + &
        n(idx_HCO) + &
        n(idx_H2O_total)*2d0 + &
        n(idx_Hj) + &
        n(idx_H2j)*2d0 + &
        n(idx_HOCj) + &
        n(idx_HCOj) + &
        n(idx_H3j)*3d0 + &
        n(idx_CHj) + &
        n(idx_CH2j)*2d0 + &
        n(idx_CH3j)*3d0 + &
        n(idx_OHj) + &
        n(idx_H2Oj)*2d0 + &
        n(idx_H3Oj)*3d0
    get_Hnuclei = nH

  end function get_Hnuclei

  !***************************
  function get_zatoms()
    use krome_commons
    implicit none
    integer::get_zatoms(nspec)

    get_zatoms(1) = 0	!E
    get_zatoms(2) = 1	!H-
    get_zatoms(3) = 6	!C-
    get_zatoms(4) = 8	!O-
    get_zatoms(5) = 1	!H
    get_zatoms(6) = 2	!HE
    get_zatoms(7) = 2	!H2
    get_zatoms(8) = 6	!C
    get_zatoms(9) = 8	!O
    get_zatoms(10) = 9	!OH
    get_zatoms(11) = 14	!CO
    get_zatoms(12) = 7	!CH
    get_zatoms(13) = 8	!CH2
    get_zatoms(14) = 12	!C2
    get_zatoms(15) = 15	!HCO
    get_zatoms(16) = 10	!H2O
    get_zatoms(17) = 16	!O2
    get_zatoms(18) = 14	!CO_total
    get_zatoms(19) = 10	!H2O_total
    get_zatoms(20) = 1	!H+
    get_zatoms(21) = 2	!HE+
    get_zatoms(22) = 2	!H2+
    get_zatoms(23) = 6	!C+
    get_zatoms(24) = 8	!O+
    get_zatoms(25) = 15	!HOC+
    get_zatoms(26) = 15	!HCO+
    get_zatoms(27) = 3	!H3+
    get_zatoms(28) = 7	!CH+
    get_zatoms(29) = 8	!CH2+
    get_zatoms(30) = 14	!CO+
    get_zatoms(31) = 9	!CH3+
    get_zatoms(32) = 9	!OH+
    get_zatoms(33) = 10	!H2O+
    get_zatoms(34) = 11	!H3O+
    get_zatoms(35) = 16	!O2+
    get_zatoms(36) = 2	!HE++
    get_zatoms(37) = 0	!CR
    get_zatoms(38) = 0	!g
    get_zatoms(39) = 0	!Tgas
    get_zatoms(40) = 0	!dummy

  end function get_zatoms

  !******************************
  function get_qeff()
    use krome_commons
    implicit none
    real*8::get_qeff(nrea)

    get_qeff(:) = 0e0

  end function get_qeff

  !**************************
  function get_free_fall_time(n)
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),m(nspec)
    real*8::rhogas,get_free_fall_time

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols))
    get_free_fall_time = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time

  !**************************
  function get_free_fall_time_rho(rhogas)
    use krome_constants
    implicit none
    real*8::rhogas,get_free_fall_time_rho

    get_free_fall_time_rho = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time_rho

  !********************************
  function get_jeans_length(n,Tgas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::m(nspec),get_jeans_length
    m(:) = get_mass()
    rhogas = max(sum(n(1:nmols)*m(1:nmols)),1d-40)
    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length

  !********************************
  function get_jeans_length_rho(n,Tgas,rhogas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::get_jeans_length_rho

    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length_rho = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length_rho

  !***************************
  !number density to column density conversion
  function num2col(ncalc,n)
    use krome_commons
    implicit none
    real*8::num2col,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    num2col = 1.87d21*(max(ncalc,1d-40)*1d-3)**(2./3.)

  end function num2col

  !***********************
  !column density to number density conversion
  function col2num(ncalc,n)
    use krome_commons
    implicit none
    real*8::col2num,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    col2num = 1d3 * (max(ncalc,1d-40)/1.87d21)**1.5

  end function col2num

  !************************
  !get electrons by balancing charges
  function get_electrons(n)
    use krome_commons
    implicit none
    real*8::get_electrons,n(nspec)

    get_electrons =  - n(idx_Hk) &
        - n(idx_Ck) &
        - n(idx_Ok) &
        + n(idx_Hj) &
        + n(idx_HEj) &
        + n(idx_H2j) &
        + n(idx_Cj) &
        + n(idx_Oj) &
        + n(idx_HOCj) &
        + n(idx_HCOj) &
        + n(idx_H3j) &
        + n(idx_CHj) &
        + n(idx_CH2j) &
        + n(idx_COj) &
        + n(idx_CH3j) &
        + n(idx_OHj) &
        + n(idx_H2Oj) &
        + n(idx_H3Oj) &
        + n(idx_O2j) &
        + 2d0*n(idx_HEjj)
    get_electrons = max(get_electrons,0d0)

  end function get_electrons

  !************************
  !get species charges
  function get_charges()
    use krome_commons
    implicit none
    integer::get_charges(nspec)

    get_charges(1) = -1.d0 	!E
    get_charges(2) = -1.d0 	!H-
    get_charges(3) = -1.d0 	!C-
    get_charges(4) = -1.d0 	!O-
    get_charges(5) = 0.d0 	!H
    get_charges(6) = 0.d0 	!HE
    get_charges(7) = 0.d0 	!H2
    get_charges(8) = 0.d0 	!C
    get_charges(9) = 0.d0 	!O
    get_charges(10) = 0.d0 	!OH
    get_charges(11) = 0.d0 	!CO
    get_charges(12) = 0.d0 	!CH
    get_charges(13) = 0.d0 	!CH2
    get_charges(14) = 0.d0 	!C2
    get_charges(15) = 0.d0 	!HCO
    get_charges(16) = 0.d0 	!H2O
    get_charges(17) = 0.d0 	!O2
    get_charges(18) = 0.d0 	!CO_total
    get_charges(19) = 0.d0 	!H2O_total
    get_charges(20) = 1.d0 	!H+
    get_charges(21) = 1.d0 	!HE+
    get_charges(22) = 1.d0 	!H2+
    get_charges(23) = 1.d0 	!C+
    get_charges(24) = 1.d0 	!O+
    get_charges(25) = 1.d0 	!HOC+
    get_charges(26) = 1.d0 	!HCO+
    get_charges(27) = 1.d0 	!H3+
    get_charges(28) = 1.d0 	!CH+
    get_charges(29) = 1.d0 	!CH2+
    get_charges(30) = 1.d0 	!CO+
    get_charges(31) = 1.d0 	!CH3+
    get_charges(32) = 1.d0 	!OH+
    get_charges(33) = 1.d0 	!H2O+
    get_charges(34) = 1.d0 	!H3O+
    get_charges(35) = 1.d0 	!O2+
    get_charges(36) = 2.d0 	!HE++
    get_charges(37) = 0.d0 	!CR
    get_charges(38) = 0.d0 	!g
    get_charges(39) = 0.d0 	!Tgas
    get_charges(40) = 0.d0 	!dummy

  end function get_charges

  !*****************************
  ! get metallicity using C as reference
  function get_metallicityC(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityC,zC,nH

    nH = get_Hnuclei(n(:))

    zC = n(idx_Ck) &
        + n(idx_C) &
        + n(idx_CH) &
        + n(idx_CH2) &
        + 2d0*n(idx_C2) &
        + n(idx_HCO) &
        + n(idx_CO_total) &
        + n(idx_Cj) &
        + n(idx_HOCj) &
        + n(idx_HCOj) &
        + n(idx_CHj) &
        + n(idx_CH2j) &
        + n(idx_COj) &
        + n(idx_CH3j)

    zC = max(zC, 0d0)

    get_metallicityC = log10(zC/nH+1d-40) - (-3.57)

    phys_metallicity = get_metallicityC

  end function get_metallicityC

  !*****************************
  ! get metallicity using O as reference
  function get_metallicityO(n)
    use krome_commons
    implicit none
    real*8::n(:),get_metallicityO,zO,nH

    nH = get_Hnuclei(n(:))

    zO = n(idx_Ok) &
        + n(idx_O) &
        + n(idx_OH) &
        + n(idx_HCO) &
        + 2d0*n(idx_O2) &
        + n(idx_CO_total) &
        + n(idx_H2O_total) &
        + n(idx_Oj) &
        + n(idx_HOCj) &
        + n(idx_HCOj) &
        + n(idx_COj) &
        + n(idx_OHj) &
        + n(idx_H2Oj) &
        + n(idx_H3Oj) &
        + 2d0*n(idx_O2j)

    zO = max(zO, 0d0)

    get_metallicityO = log10(zO/nH+1d-40) - (-3.31)

    phys_metallicity = get_metallicityO

  end function get_metallicityO

end module krome_getphys
!This module contains the functions and subroutines
! needed to evaluate the adiabatic index.

!############### MODULE ##############
module krome_gadiab
contains

  !#KROME_header

  !**************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part
  function gamma_pop(array_part,dT_part,min_part,Tgasin)
    implicit none
    real*8::array_part(:),dT_part
    real*8::min_part,Tgas,gamma_pop,Tgas2,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !temperature above minimum data point
    inTgas = max(Tgasin,min_part)

    !data index
    idx = (inTgas-min_part)/dT_part+1
    !corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !ln of partition functions (3 points forward)
    logz = log(array_part(idx))
    logz1 = log(array_part(idx+1))
    logz2 = log(array_part(idx+2))

    !derivative for mean energy (2 points forward)
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv1 = (emed2-emed1)/dT_part

    !next point temperature
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of partition functions
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part(idx+3))

    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv2 = (emed2-emed1)/dT_part

    !interpolation for 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns result
    gamma_pop = Cv

  end function gamma_pop

  !*****************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part, for H2 with a ortho/para
  ! ratio of opratio. Needs even and odd partition functions.
  function gamma_pop_H2(array_part_even,array_part_odd,dT_part,&
        min_part,Tgasin,opratio)
    implicit none
    real*8::array_part_even(:),array_part_odd(:),dT_part,zcut(4)
    real*8::min_part,Tgas,opratio,gamma_pop_H2,Tgas2,a,b,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !Tgas above the data limit
    inTgas = max(Tgasin,min_part)

    !exponents for ortho/para ratio
    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    !index in the data for the given Tgas
    idx = (inTgas-min_part)/dT_part+1
    !get the corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !needed for ortho partition function (see Boley+2007)
    zcut(1) = exp(2d0*85.4/Tgas)
    zcut(2) = exp(2d0*85.4/(Tgas+dT_part))
    zcut(3) = exp(2d0*85.4/(Tgas+2d0*dT_part))
    zcut(4) = exp(2d0*85.4/(Tgas+3d0*dT_part))

    !ln of the composite partition function
    logz = log(array_part_even(idx)**b*(3d0*array_part_odd(idx)*zcut(1))**a)
    logz1 = log(array_part_even(idx+1)**b*(3d0*array_part_odd(idx+1)*zcut(2))**a)
    logz2 = log(array_part_even(idx+2)**b*(3d0*array_part_odd(idx+2)*zcut(3))**a)
    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the left point
    Cv1 = (emed2-emed1)/dT_part

    !Tgas of the right point
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of the composite function
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part_even(idx+3)**b*(3d0*array_part_odd(idx+3)*zcut(4))**a)
    !derivative for the mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the right point
    Cv2 = (emed2-emed1)/dT_part

    !interpolation of 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns the result
    gamma_pop_H2 = Cv
  end function gamma_pop_H2

  !**************************
  !function to get the partition function
  ! of H2 at Tgas with a orto-para ratio
  ! equal to opratio
  function zfop(Tgas,opratio)
    implicit none
    real*8::Tgas,zfop,brot,ibTgas
    real*8::a,b,zo,zp,opratio
    integer::j,jmax,j1
    brot = 85.4d0 !H2 rotational constant in K
    zo = 0d0 !sum for ortho partition function
    zp = 0d0 !sum for para partition function
    jmax = 10 !number of terms in sum

    ibTgas = brot/Tgas !pre-calc

    !loop over levels
    do j=0,jmax,2 !step 2
      j1 = j + 1
      zp = zp + (2d0*j+1d0) * exp(-j*(j+1d0)*ibTgas)
      zo = zo + 3d0 * (2d0*j1+1d0) * exp(-j1*(j1+1d0)*ibTgas)
    end do

    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    zfop = (zp**b * zo**a*exp(-2d0*ibTgas)) !final partition f

  end function zfop

  !*********************
  !get the partition function at Tgas
  ! of a diatom with rotational constant
  ! brot in K
  function zf(Tgas,brot)
    real*8::Tgas,zf,brot,z,ibTgas
    integer::j,jmax
    jmax = 10 !number of levels

    ibTgas = brot/Tgas !store
    z = 0d0
    !loop on levels
    do j=0,jmax
      z = z + (2d0*j+1d0)*exp(-j*(j+1d0)*ibTgas)
    end do

    zf = z

  end function zf

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of H2 with
  ! an ortho-para ratio of opratio
  function gamma_rotop(Tgas_in,opratio)
    implicit none
    real*8::gamma_rotop,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,opratio

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zfop(Tgas+dT,opratio)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zfop(Tgas,opratio)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zfop(Tgas+dT+dT,opratio))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rotop = (prot2-prot1)*idT

  end function gamma_rotop

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of a diatom
  ! with rotational constant brot in K
  function gamma_rot(Tgas_in,brot)
    implicit none
    real*8::gamma_rot,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,brot

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zf(Tgas+dT,brot)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zf(Tgas,brot)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zf(Tgas+dT+dT,brot))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rot = (prot2-prot1)*idT

  end function gamma_rot

  !*********************
  !get gamma
  function gamma_index(n)
    use krome_commons
    implicit none
    real*8::n(:),gamma_index,krome_gamma

    krome_gamma = 1.66666666667d0

    gamma_index = krome_gamma
  end function gamma_index

end module krome_gadiab
!This module contains functions and subroutines
! for the surface chemistry, including adsorption, desorption, chemisorption
! and icy grains.

!############### MODULE ##############
module krome_grfuncs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-10-24 12:49:09
  !  Changeset b21f657
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !**********************
  !get Tdust from tables, K
  function get_table_Tdust(n) result(Tdust)
    use krome_commons
    use krome_fit
    implicit none
    real*8,intent(in)::n(nspec)
    real*8::ntot,Tdust,Tgas

    Tgas = n(idx_Tgas)

    !default, K
    Tdust = 1d0

    !total densitym, cm-3
    ntot = sum(n(1:nmols))

    !zero density returns default
    if(ntot==0d0) return

    !get dust temperature from table, K
    Tdust = 1d1**fit_anytab2D(dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas, &
        log10(ntot), log10(Tgas))

  end function get_table_Tdust

  !**********************
  !adsorpion rate Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014
  function dust_adsorption_rate(nndust,ims,stick,adust2,sqrTgas)
    use krome_constants
    implicit none
    real*8::dust_adsorption_rate,nndust,ims,stick,adust2,sqrTgas

    dust_adsorption_rate = nndust * pi * adust2 &
        * pre_kvgas_sqrt * ims * sqrTgas &
        * stick

  end function dust_adsorption_rate

  !*****************************
  !desorption rate Cazaux+2010, Hocuk+2014
  function dust_desorption_rate(fice,expEice,expEbare)
    implicit none
    real*8::dust_desorption_rate
    real*8::fice,expEice,expEbare,nu0,fbare

    nu0 = 1d12 !1/s
    fbare = 1d0 - fice
    dust_desorption_rate = nu0 * (fbare * expEbare &
        + fice * expEice)

  end function dust_desorption_rate

  !**************************
  function dust_2body_rate(p,invphi,fice,expEice1,expEice2,&
        expEbare1,expEbare2,pesc_ice,pesc_bare)
    use krome_constants
    implicit none
    real*8::fice,expEice1,expEice2,expEbare1,expEbare2,invphi
    real*8::nu0,p,dust_2body_rate,fbare,pesc_ice,pesc_bare

    !no need to calculate this if the dust is not present
    dust_2body_rate = 0d0

    fbare = 1d0-fice
    nu0 = 1d12 ! 1/s
    dust_2body_rate = fbare * (expEbare1 + expEbare2) * pesc_bare &
        + fice * (expEice1 + expEice2) * pesc_ice
    dust_2body_rate = dust_2body_rate * p * nu0 * invphi

  end function dust_2body_rate

  !******************
  function krate_2bodySi(n,idx1,idx2,Ea,Tdust) result(krate)
    use krome_commons
    implicit none
    real*8,intent(in)::n(nspec),Ea,Tdust
    integer,intent(in)::idx1,idx2
    real*8::krate,amin,amax,pexp,d2g,rho0

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    krate = krate_2body(n(:),idx1,idx2,amin,amax,pexp,d2g,rho0,Ea,Tdust)

  end function krate_2bodySi

  !********************
  function krate_2body(n,idx1,idx2,amin,amax,pexp,d2g,rho0, &
        Ea,Tdust) result(krate)
    use krome_commons
    use krome_constants
    use krome_getphys
    implicit none
    integer,intent(in)::idx1,idx2
    real*8,intent(in)::n(nspec),amin,amax,pexp,d2g,rho0,Ea,Tdust
    real*8::rhog,p3,p4,ndns,krate,mred,fice,fbare,Preac
    real*8::iTd23,Ebare(nspec),Eice(nspec),mass(nspec)
    real*8,parameter::app2=(3d-8)**2 !cm^2 (Hocuk+2015)
    real*8,parameter::nu0=1d12 !1/s
    real*8,parameter::hbar=planck_erg/2d0/pi !erg*s
    real*8,parameter::ar=1d-8 !cm

    mass(:) = get_mass()

    !gas density, g/cm3
    rhog = sum(mass(1:nmols)*n(1:nmols))

    !exponentes
    p3 = pexp + 3d0
    p4 = pexp + 4d0

    !number of sites cm-3/mly
    ndns = rhog/(4d0/3d0*rho0*app2)*(amax**p3-amin**p3) &
        / (amax**p4-amin**p4) * p4 / p3

    !ice/bare fraction
    fbare = 1d0
    fice = (n(idx_H2O_total)-n(idx_H2O))/ndns
    fbare = 1d0 - fice

    !reduced mass
    mred = mass(idx1)*mass(idx2)/(mass(idx1)+mass(idx2))

    !tunneling probability
    Preac = exp(-2d0*ar/hbar*sqrt(2d0*mred*Ea*boltzmann_erg))

    !exponent
    iTd23 = 2d0/3d0/Tdust

    !get Ebind, K
    Ebare(:) = get_Ebind_bare()
    Eice(:) = get_Ebind_ice()

    !compute rate
    krate = fbare*(exp(-Ebare(idx1)*iTd23)+exp(-Ebare(idx2)*iTd23))
    krate = krate + fice*(exp(-Eice(idx1)*iTd23)+exp(-Eice(idx2)*iTd23))

    !rate in cm3/s
    krate = nu0*Preac/ndns*krate

  end function krate_2body

  !*************************
  function dust_get_inv_phi(asize2,nndust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::iapp2,dust_get_inv_phi(ndust),asize2(ndust)
    real*8::nndust(ndust),dephi
    integer::i

    iapp2 = (3d-8)**2 !1/cm2
    do i=1,ndust
      dust_get_inv_phi(i) = 0d0
      dephi = (4d0 * nndust(i) * pi * asize2(i))
      if(dephi.le.0d0) cycle
      dust_get_inv_phi(i) = iapp2 / dephi
    end do

  end function dust_get_inv_phi

  !****************************
  !returns an array with the sticking coefficient for each bin
  ! following Hollenbach+McKee 1979
  function dust_stick_array(Tgas,Tdust)
    use krome_commons
    implicit none
    real*8::dust_stick_array(ndust),Tgas,Tdust(ndust)
    real*8::Tg100,Td100
    integer::i

    Tg100 = Tgas * 1d-2
    do i=1,ndust
      Td100 = Tdust(i) * 1d-2
      dust_stick_array(i) = 1d0/(1d0+.4d0*sqrt(Tg100+Td100) &
          + .2d0*Tg100 + 0.08d0*Tg100**2)
    end do

  end function dust_stick_array

  !*************************
  function dust_stick(Tgas,Tdust)
    implicit none
    real*8,intent(in)::Tgas,Tdust
    real*8::dust_stick
    real*8::Tg100,Td100

    Tg100 = Tgas * 1d-2
    Td100 = Tdust * 1d-2
    dust_stick = 1d0/(1d0 + 0.4d0*sqrt(Tg100+Td100) &
        + 0.2d0*Tg100 + 0.08d0*Tg100**2)

  end function dust_stick

  !****************************
  !sticking rate (1/s), assuming power-law dust distribution
  ! example rate is
  !  @format:idx,R,P,rate
  !  1,CO,CO_ice,krate_stick(n(:),idx_CO,1d-7,1d-5,-3.5,3d0,1d-2)
  ! n(:): internal status array (number densities, temeperature, etc...)
  ! idx : index of the sticking species, e.g. idx_CO
  ! Tdust: dust temperature (assume same for all bins), K
  ! amin: min grain size, cm
  ! amax: max grain size, cm
  ! pexp: power-law exponent, usually -3.5
  ! rho0: bulk material density, g/cm3, e.g. 3 g/cm3 for silicates
  ! d2g: dust to gass mass ratio, usually 0.01
  function krate_stick(n,idx,Tdust,amin,amax,pexp,rho0,d2g) result(k)
    use krome_constants
    use krome_commons
    use krome_getphys
    implicit none
    real*8,intent(in)::n(nspec),Tdust,amin,amax,pexp,rho0,d2g
    real*8::k,imass(nspec),p4,p3,mass(nspec),rhod
    integer,intent(in)::idx

    !get inverse mass squared
    imass(:) = get_imass_sqrt()
    !get masses
    mass(:) = get_mass()
    !derived exponents
    p3 = pexp + 3.
    p4 = pexp + 4.

    !total dust density, g/cm3
    rhod = sum(n(1:nmols)*mass(1:nmols))*d2g

    !compute rate (1/s) coefficient assuming normalization
    k = pre_kvgas_sqrt*sqrt(n(idx_Tgas)) * imass(idx) &
        * rhod / (4./3.*rho0) * p4 / p3 &
        * (amax**p3-amin**p3) / (amax**p4-amin**p4) &
        * dust_stick(n(idx_Tgas),Tdust)

  end function krate_stick

  !********************************
  !compact version of krate_stick
  function krate_stickSi(n,idx,Tdust) result(k)
    use krome_commons
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,amin,amax,d2g,rho0,pexp

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    k = krate_stick(n(:),idx,Tdust,amin,amax,pexp,rho0,d2g)

  end function krate_stickSi

  !***************************
  !evaporation rate, 1/s
  function krate_evaporation(n,idx,Tdust) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,Ebind(nspec),nu0

    nu0 = 1d12 !1/s
    Ebind(:) = get_EbindBare()

    k = nu0 * exp(-Ebind(idx)/Tdust)

  end function krate_evaporation

  !***************************
  !non-thermal evaporation rate (1/s) following Hollenbach 2009,
  ! http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0809.1642
  !Gnot is the habing flux (1.78 is Draine)
  !Av is the visual extinction
  !crflux the ionization flux of cosmic rays, 1/s
  !yield is the efficiency of the photons to desorb the given molecule
  function krate_nonthermal_evaporation(idx, Gnot, Av, crflux, yield) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,parameter::crnot=1.3d-17
    real*8,parameter::Fnot=1d8 !desorbing photons flux, 1/s
    real*8,parameter::ap2=(3d-8)**2 !sites separation squared, cm2
    real*8,intent(in)::Gnot, Av, crflux, yield
    real*8::k,f70,kevap70(nspec)

    f70 = 3.16d-19*crflux/crnot
    kevap70(:) = get_kevap70()

    k = Gnot*Fnot*ap2*yield*exp(-1.8*Av)
    k = k + f70*kevap70(idx)

  end function krate_nonthermal_evaporation

  !***************************
  function dust_ice_fraction_array(invphi,nH2O)
    use krome_constants
    use krome_commons
    implicit none
    integer::i
    real*8::dust_ice_fraction_array(ndust)
    real*8::invphi(ndust),nH2O(ndust)

    do i=1,ndust
      dust_ice_fraction_array(i) = min(nH2O(i) * invphi(i), 1d0)
    end do

  end function dust_ice_fraction_array

  !*****************************
  function get_Ebareice_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice_exp_array(:) = 0d0

  end function get_Ebareice_exp_array

  !*****************************
  function get_Ebareice23_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice23_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice23_exp_array(:) = 0d0

  end function get_Ebareice23_exp_array

  !************************
  !returns the binding energy for ice coated grain (K)
  function get_Ebind_ice()
    use krome_commons
    implicit none
    real*8::get_Ebind_ice(nspec)

    get_Ebind_ice(:) = 0d0

  end function get_Ebind_ice

  !************************
  !returns the binding energy for bare grain (K)
  function get_Ebind_bare()
    use krome_commons
    implicit none
    real*8::get_Ebind_bare(nspec)

    get_Ebind_bare(:) = 0d0

  end function get_Ebind_bare

  !************************
  !returns the index of the parent dust bin (0 if none)
  function get_parent_dust_bin()
    use krome_commons
    implicit none
    integer::get_parent_dust_bin(nspec)

    get_parent_dust_bin(:) = 0

  end function get_parent_dust_bin

  !*****************************
  function get_exp_table(ain,invT)
    use krome_commons
    implicit none
    integer::ia
    real*8::get_exp_table,a,invT,ain
    real*8::x1a,f1,f2

    a = ain*invT
    a = min(a, exp_table_aMax - exp_table_da)

    ia = (a-exp_table_aMin) * exp_table_multa + 1
    ia = max(ia,1)

    x1a = (ia-1)*exp_table_da

    f1 = exp_table(ia)
    f2 = exp_table(ia+1)

    get_exp_table = (a-x1a) * exp_table_multa * (f2-f1) + f1

  end function get_exp_table

end module krome_grfuncs
!This module mainly contains shielding routine and
! function to initialize radiation background (e.g. Planck).

!############### MODULE ##############
module krome_phfuncs
contains

  !****************************
  !dust shielding factor
  function shield_dust(n,Tgas,gam)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::shield_dust,n(:),Tgas,gam,eff_d2g
    real*8::sigma_d,NHtot

    eff_d2g = dust2gas_ratio
    sigma_d = 2d-21*eff_d2g*gam !Richings et al. 2014
    !sigma_d = 2d-21 !Glover+2007
    !sigma_d = 4d-22 !Richings+ 2014
    !sigma_d = 4d-21 !Gnedin 2009

    NHtot = 0d0
    NHtot  = NHtot + num2col(n(idx_H),n(:))
    NHtot  = NHtot + num2col(n(idx_Hj),n(:))
    NHtot  = NHtot + 2d0 * num2col(n(idx_H2),n(:))

    shield_dust = exp(-sigma_d*NHtot)

  end function shield_dust

  !*******************
  !apply a shielding to Habing flux
  subroutine calcHabingThick(n,Tgas)
    use krome_commons
    implicit none
    real*8::getHabingThick,n(:),Tgas

    GHabing = GHabing_thin * shield_dust(n(:),Tgas,0.665d0)

  end subroutine calcHabingThick

  !*********************
  !return the ratio between the current flux an Draine's
  function get_ratioFluxDraine()
    implicit none
    real*8::get_ratioFluxDraine

    !7.95d-8 eV/cm2/sr is the integrated Draine flux
    get_ratioFluxDraine = get_integratedFlux()/7.95d-8

  end function get_ratioFluxDraine

  !**********************
  !return the curred integrated flux (eV/cm2/sr)
  ! as I(E)/E*dE
  function get_integratedFlux()
    use krome_commons
    implicit none
    integer::j
    real*8::get_integratedFlux,dE

    get_integratedFlux = 0d0
    do j=1,nPhotoBins
      dE = photoBinEdelta(j)
      get_integratedFlux = get_integratedFlux &
          + photoBinJ(j)*dE/photoBinEmid(j)
    end do

  end function get_integratedFlux

  !**********************
  !planck function in eV/s/cm2/Hz/sr
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB(x,Tbb)
    use krome_constants
    implicit none
    real*8::Tbb,x,xexp,planckBB

    !exponent
    xexp = x/boltzmann_eV/Tbb

    !default value
    planckBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
      planckBB = 2d0*x**3/planck_eV**2/clight**2 &
          / (exp(xexp)-1d0)
    end if

  end function planckBB

  !********************
  !planck function dTdust differential
  ! in eV/s/cm2/Hz/sr/K, where
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB_dT(x,Tbb)
    use krome_constants
    real*8::a,b,x,Tbb,xexp,planckBB_dT

    b = 1d0/boltzmann_eV
    xexp = b*x/Tbb

    planckBB_dT = 0d0

    if(xexp<3d2) then
      a = 2d0/planck_eV**2/clight**2
      planckBB_dT = a*b*x**4/Tbb/Tbb * exp(xexp)/(exp(xexp)-1d0)**2
    end if

  end function planckBB_dT

  !***********************
  !shielding function selected with -shield option
  function krome_fshield(n,Tgas)
    implicit none
    real*8::krome_fshield,n(:),Tgas

    krome_fshield = 1d0 !default shielding value

  end function krome_fshield

  !**************************
  !shielding function for H2O+ and H3O+
  ! following Glover+2010 MNRAS sect 2.2 eqn.4
  function fHnOj(Av)
    implicit none
    real*8::fHnOj,Av
    if(Av.le.15d0) then
      fHnOj = exp(-2.55*Av+0.0165*Av**2)
    else
      fHnOj = exp(-2.8*Av)
    end if
  end function fHnOj

  !******************************
  !self-shielding for H2
  ! following Glover+2010 MNRAS sect2.2 eqn.6
  ! N: column density (cm-2)
  ! b: doppler broadening (cm/s)
  function fselfH2(N, b)
    implicit none
    real*8::fselfH2,N,b,x,b5

    x = N*2d-15 !normalized column density (#)
    b5 = b*1d-5 !normalized doppler broadening (#)

    fselfH2 = 0.965d0/(1+x/b5)**2 + &
        0.035d0/sqrt(1d0+x) * &
        exp(max(-8.5d-4*sqrt(1+x),-250.))

  end function fselfH2

end module krome_phfuncs

!############### MODULE ##############
module krome_subs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2018-10-24 12:49:10
  !  Changeset b21f657
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !************************
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_constants
    use krome_user_commons
    use krome_getphys
    use krome_grfuncs
    use krome_phfuncs
    use krome_fit
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),kmax
    real*8::invTgas
    real*8::small,nmax
    integer::i
    real*8::ncr11  !preproc from coevar
    real*8::kHdissH2 !preproc from coevar
    real*8::invTe !preproc from coevar
    real*8:: Te  !preproc from coevar
    real*8::kHdiss !preproc from coevar
    real*8::Hnuclei  !preproc from coevar
    real*8:: lnTe  !preproc from coevar
    real*8:: T32  !preproc from coevar
    real*8::ncrdiss !preproc from coevar
    real*8::sumsav !preproc from coevar
    real*8::tabTdust  !preproc from coevar
    real*8::logT !preproc from coevar
    real*8::u3  !preproc from coevar
    real*8::kh11  !preproc from coevar
    real*8::u1  !preproc from coevar
    real*8:: invT  !preproc from coevar
    real*8::ncrdissH2 !preproc from coevar
    real*8:: invT32  !preproc from coevar
    real*8::asav2 !preproc from coevar
    real*8::asav3 !preproc from coevar
    real*8::asav0 !preproc from coevar
    real*8::asav1 !preproc from coevar
    real*8::asav6 !preproc from coevar
    real*8::asav7 !preproc from coevar
    real*8::asav4 !preproc from coevar
    real*8:: T  !preproc from coevar
    real*8:: HnOj  !preproc from coevar
    real*8::a11 !preproc from coevar
    real*8::a10 !preproc from coevar
    real*8::a12 !preproc from coevar
    real*8::kLdiss !preproc from coevar
    real*8::ntot !preproc from coevar
    real*8::a1 !preproc from coevar
    real*8::a3 !preproc from coevar
    real*8::a2 !preproc from coevar
    real*8::a5 !preproc from coevar
    real*8::ax11 !preproc from coevar
    real*8::a7 !preproc from coevar
    real*8::a6 !preproc from coevar
    real*8::a9 !preproc from coevar
    real*8::a8 !preproc from coevar
    real*8::kl11  !preproc from coevar
    real*8::sumsbv !preproc from coevar
    real*8::u2  !preproc from coevar
    real*8::invsqrT !preproc from coevar
    real*8::savthr !preproc from coevar
    real*8::bsav1 !preproc from coevar
    real*8::bsav0 !preproc from coevar
    real*8::bsav3 !preproc from coevar
    real*8::bsav2 !preproc from coevar
    real*8::bsav5 !preproc from coevar
    real*8::bsav4 !preproc from coevar
    real*8::bsav7 !preproc from coevar
    real*8::bsav6 !preproc from coevar
    real*8::kLdissH2 !preproc from coevar
    real*8::asav5 !preproc from coevar
    real*8::a4 !preproc from coevar
    !Tgas is in K
    Tgas = max(n(idx_Tgas), phys_Tcmb)
    Tgas = min(Tgas,1d9)

    !maxn initialization can be removed and small can be
    ! replaced with a proper value according to the environment
    nmax = max(maxval(n(1:nmols)),1d0)
    small = 1d-40/(nmax*nmax*nmax)

    invTgas = 1.d0/Tgas !inverse of T (1/K)

    Hnuclei = get_Hnuclei(n(:))
    ntot = sum(n(1:nmols))
    T32 = Tgas/3d2
    Te = Tgas*8.617343d-5
    invT = 1d0/Tgas
    lnTe = log(Te)
    T = Tgas
    invT32 = 1d0/T32
    invTe = 1d0/Te
    logT = log10(Tgas)
    invsqrT = 1d0/sqrt(Tgas)
    kl11 = 1d1**(-27.029d0+3.801d0*logT-29487d0*invT)
    kh11 = 1d1**(-2.729d0-1.75d0*logT-23474d0*invT)
    ncr11 = 1d1**(5.0792d0*(1d0-1.23d-5*(Tgas-2d3)))
    ax11 = 1.d0/(1.d0+(Hnuclei/(ncr11+1d-40)))
    a1 = 1.3500e-09
    a2 = 9.8493e-02
    a3 = 3.2852e-01
    a4 = 5.5610e-01
    a5 = 2.7710e-07
    a6 = 2.1826e+00
    a7 = 6.1910e-03
    a8 = 1.0461e+00
    a9 = 8.9712e-11
    a10 = 3.0424e+00
    a11 = 3.2576e-14
    a12 = 3.7741e+00
    asav0 = -1.9153214d2
    asav1 = 4.0129114d2
    asav2 = -3.7446991d2
    asav3 = 1.9078410d2
    asav4 = -5.7263467d1
    asav5 = 1.0133210d1
    asav6 = -9.8012853d-1
    asav7 = 4.0023414d-2
    bsav0 = -8.8755774d3
    bsav1 = 1.0081246d4
    bsav2 = -4.8606622d3
    bsav3 = 1.2889659d3
    bsav4 = -2.0319575d2
    bsav5 = 1.9057493d1
    bsav6 = -9.8530668d-1
    bsav7 = 2.1675387d-2
    sumsav = asav0+asav1*log10(Tgas)+asav2*(log10(Tgas))**2+asav3*(log10(Tgas))**3+asav4*(log10(Tgas))**4+asav5*(log10(Tgas))**5+asav6*(log10(Tgas))**6+asav7*(log10(Tgas))**7
    sumsbv = bsav0+bsav1*log10(Tgas)+bsav2*(log10(Tgas))**2+bsav3*(log10(Tgas))**3+bsav4*(log10(Tgas))**4+bsav5*(log10(Tgas))**5+bsav6*(log10(Tgas))**6+bsav7*(log10(Tgas))**7
    savthr = 2.12362754d4
    kHdiss = 3.52d-9*exp(-4.39d4*invT) + 1d-40
    kLdiss = 6.67d-12*sqrt(Tgas)*exp(-(1d0+63590.*invT)) + 1d-40
    ncrdiss = 1d1**(3. - 0.416*log10(Tgas/1d4) - 0.327*log10(Tgas/1d4)**2)
    kHdissH2 = 1.3d-9*exp(-5.33d4*invT) + 1d-40
    kLdissH2 = 5.996d-30*Tgas**4.1881/(1.+6.761d-6*Tgas)**5.6881 * exp(-5.46574d4*invT) + 1d-40
    ncrdissH2 = 1d1**(4.845 - 1.3*log10(Tgas/1d4) + 1.62*log10(Tgas/1d4)**2)
    u1 = 11.26d0*invTe
    u2 = 8.2d0*invTe
    u3 = 13.6*invTe
    tabTdust = get_table_Tdust(n(:))
    HnOj = fHnOj(user_Av)

    k(:) = small !inizialize coefficients

    !H + E -> H+ + E + E
    k(1) = small + (exp(-32.71396786d0+13.5365560d0&
        *lnTe-5.73932875d0*(lnTe**2)+1.56315498d0&
        *(lnTe**3)-0.28770560d0*(lnTe**4)+3.48255977d-2&
        *(lnTe**5)-2.63197617d-3*(lnTe**6)+1.11954395d-4&
        *(lnTe**7)-2.03914985d-6*(lnTe**8)))

    !H+ + E -> H
    if(Tgas.GE.2.73d0 .and. Tgas.LE.5.5d3) then
      k(2) = small + (3.92d-13&
          *invTe**0.6353d0)
    end if

    !H+ + E -> H
    if(Tgas.GT.5.5d3 .and. Tgas.LT.1d8) then
      k(3) = small + (exp(-28.61303380689232d0-0.7241125657826851d0&
          *lnTe-0.02026044731984691d0*lnTe**2-0.002380861877349834d0&
          *lnTe**3-0.0003212605213188796d0&
          *lnTe**4-0.00001421502914054107d0&
          *lnTe**5+4.989108920299513d-6*lnTe**6+5.755614137575758d-7&
          *lnTe**7-1.856767039775261d-8*lnTe**8-3.071135243196595d-9&
          *lnTe**9))
    end if

    !HE + E -> HE+ + E + E
    k(4) = small + (exp(-44.09864886d0+23.91596563d0&
        *lnTe-10.7532302d0*(lnTe**2)+3.05803875d0&
        *(lnTe**3)-0.56851189d0*(lnTe**4)+6.79539123d-2&
        *(lnTe**5)-5.00905610d-3*(lnTe**6)+2.06723616d-4&
        *(lnTe**7)-3.64916141d-6*(lnTe**8)))

    !HE+ + E -> HE
    if(Tgas.GE.2.73d0 .and. Tgas.LE.9.28d3) then
      k(5) = small + (3.92d-13&
          *invTe**0.6353d0)
    end if

    !HE+ + E -> HE
    if(Tgas.GT.9.28d3 .and. Tgas.LT.1d8) then
      k(6) = small + (1.54d-9&
          *(1.d0+0.3d0&
          /exp(8.099328789667d0&
          *invTe))&
          /(exp(40.49664394833662d0*invTe)&
          *Te**1.5d0)+3.92d-13&
          /Te**0.6353d0)
    end if

    !HE+ + E -> HE++ + E + E
    k(7) = small + (exp(-68.71040990212001d0+43.93347632635d0&
        *lnTe-18.48066993568d0*lnTe**2+4.701626486759002d0&
        *lnTe**3-0.7692466334492d0*lnTe**4+0.08113042097303d0&
        *lnTe**5-0.005324020628287001d0*lnTe**6+0.0001975705312221d0&
        *lnTe**7-3.165581065665d-6*lnTe**8))

    !HE+ + H -> HE + H+
    k(8) = small + (1.2d-15*(T32)*.25)

    !HE + H+ -> HE+ + H
    if(Tgas.LT.1d4) then
      k(9) = small + (1.26d-9&
          *Tgas**(-.75)*exp(-1.275d5*invT))
    end if

    !HE + H+ -> HE+ + H
    if(Tgas.GE.1d4) then
      k(10) = small + (4d-37&
          *Tgas**4.74)
    end if

    !H2 + HE -> H + H + HE
    k(11) = small + (kh11**(1.-ax11)&
        *kl11**ax11)

    !H2 + HE+ -> HE + H2+
    k(12) = small + (7.2d-15)

    !H2 + HE+ -> HE + H + H+
    k(13) = small + (3.7d-14*exp(-35d0&
        *invT))

    !H2 + HE+ -> HE+ + H + H
    k(14) = small + (3d-11*sqrt(T32)&
        *exp(-5.2d4*invT))

    !HE++ + E -> HE+
    k(15) = small + (1.891d-10/(sqrt(Tgas&
        /9.37)&
        *(1.+sqrt(Tgas/9.37))**0.2476&
        *(1.+sqrt(Tgas&
        /2.774d6))**1.7524))

    !H + E -> H-
    k(16) = small + (1.4d-18*Tgas**0.928&
        *exp(-Tgas&
        /16200.))

    !H- + H -> H2 + E
    k(17) = small + (a1*(Tgas**a2+a3&
        *Tgas**a4+a5*Tgas**a6)&
        /(1.+a7*Tgas**a8+a9*Tgas**a10+a11&
        *Tgas**a12))

    !H + H+ -> H2+
    if(Tgas.LT.3d1) then
      k(18) = small + (2.10e-20&
          *(Tgas&
          /30.)**(-0.15))
    end if

    !H + H+ -> H2+
    if(Tgas.GE.3d1) then
      k(19) = small + (10**(-18.20-3.194&
          *log10(Tgas)+1.786*(log10(Tgas))**2-0.2072&
          *(log10(Tgas))**3))
    end if

    !H2+ + H -> H2 + H+
    k(20) = small + (6d-10)

    !H2 + H+ -> H2+ + H
    if(Tgas.LT.1d5) then
      k(21) = small + (exp(-savthr&
          *invT)*1d1**sumsav)
    end if

    !H2 + H+ -> H2+ + H
    if(Tgas.GE.1d5) then
      k(22) = small + (exp(-savthr&
          *invT)*1d1**sumsbv)
    end if

    !H2 + E -> H + H + E
    k(23) = small + (4.38d-10*exp(-1.02d5&
        *invT)*Tgas**(0.35))

    !H2 + H -> H + H + H
    k(24) = small + (1d1**(log10(kHdiss)-log10(kHdiss&
        /kLdiss)/(1d0+ntot/ncrdiss)))

    !H- + E -> H + E + E
    k(25) = small + (exp(-18.01849334273d0+2.360852208681d0&
        *lnTe-0.2827443061704d0*lnTe**2+0.01623316639567d0&
        *lnTe**3-0.03365012031362999d0*lnTe**4+0.01178329782711d0&
        *lnTe**5-0.001656194699504d0*lnTe**6+0.0001068275202678d0&
        *lnTe**7-2.631285809207d-6*lnTe**8))

    !H- + H -> H + H + E
    if(Tgas.LE.1.16d3) then
      k(26) = small + (2.56d-9&
          *Te**1.78186)
    end if

    !H- + H -> H + H + E
    if(Tgas.GT.1.16d3) then
      k(27) = small + (exp(-20.37260896533324d0+1.139449335841631d0&
          *lnTe-0.1421013521554148d0*lnTe**2+0.00846445538663d0&
          *lnTe**3-0.0014327641212992d0*lnTe**4+0.0002012250284791d0&
          *lnTe**5+0.0000866396324309d0*lnTe**6-0.00002585009680264d0&
          *lnTe**7+2.4555011970392d-6*lnTe**8-8.06838246118d-8&
          *lnTe**9))
    end if

    !H- + H+ -> H + H
    k(28) = small + ((2.96d-6&
        /sqrt(Tgas) - 1.73d-9 + 2.50d-10&
        *sqrt(Tgas)))

    !H- + H+ -> H2+ + E
    k(29) = small + (1.d-8*Tgas**(-0.4))

    !H2+ + E -> H + H
    if(Tgas.LE.6.17d2) then
      k(30) = small + (1.d-8)
    end if

    !H2+ + E -> H + H
    if(Tgas.GT.6.17d2) then
      k(31) = small + (1.32d-6&
          *Tgas**(-0.76))
    end if

    !H2+ + H- -> H + H2
    k(32) = small + (5d-7*sqrt(1d2*invT))

    !H2 + H2 -> H2 + H + H
    k(33) = small + (1d1**(log10(kHdissH2)-log10(kHdissH2&
        /kLdissH2)/(1d0+ntot/ncrdissH2)))

    !H + H + HE -> H2 + HE
    k(34) = small + (6.9d-32&
        *Tgas**(-.4))

    !H + H + H -> H2 + H
    k(35) = small + (6d-32&
        *Tgas**(-.25) + 2d-31*Tgas**(-.5))

    !H2 + H + H -> H2 + H2
    k(36) = small + ((6d-32&
        *Tgas**(-0.25) + 2d-31*Tgas**(-.5)) &
        / 8d0)

    !C+ + E -> C
    if(Tgas.LE.7950d0) then
      k(37) = small + (4.67d-12&
          *(T32)**(-0.6))
    end if

    !C+ + E -> C
    if(Tgas.GT.7950d0 .and. Tgas.LE.21140d0) then
      k(38) = small + (1.23d-17&
          *(T32)**2.49*exp(21845.6d0*invT))
    end if

    !C+ + E -> C
    if(Tgas.GT.21140d0) then
      k(39) = small + (9.62d-8&
          *(T32)**(-1.37)*exp(-115786.2d0*invT))
    end if

    !O+ + E -> O
    if(Tgas.LE.4d2) then
      k(40) = small + (1.30d-10&
          *(Tgas)**(-0.64))
    end if

    !O+ + E -> O
    if(Tgas.GT.4d2) then
      k(41) = small + (1.41d-10&
          *(Tgas)**(-0.66) + 7.4d-4*(Tgas)**(-1.5)*exp(-1.75d5*invT)&
          *(1d0 + 0.062d0*exp(-1.45d5*invT)))
    end if

    !C + E -> C+ + E + E
    k(42) = small + (6.85d-8*u1**0.25&
        *exp(-u1)&
        /(0.193d0+u1))

    !O + E -> O+ + E + E
    k(43) = small + (3.59d-8*u3**0.34&
        *exp(-u3)&
        /(0.073d0+u3))

    !O+ + H -> O + H+
    k(44) = small + (4.99d-11&
        *Tgas**0.405 + 7.54d-10*invT**(0.458))

    !O + H+ -> O+ + H
    k(45) = small + ((1.08d-11&
        *Tgas**0.517 + 4d-10*Tgas**(0.00669))*exp(-2.27d2*invT))

    !O + HE+ -> O+ + HE
    k(46) = small + (4.991d-15*(Tgas&
        *1d-4)**0.3794*exp(-Tgas*8.9206d-7) + 2.78d-15*(Tgas&
        *1d-4)**(-0.2163)*exp(-Tgas*1.2258d-6))

    !C + H+ -> C+ + H
    k(47) = small + (3.9d-16*Tgas**(0.213))

    !C+ + H -> C + H+
    k(48) = small + (6.08d-14*(Tgas&
        *1d-4)**(1.96)*exp(-1.7d5*invT))

    !C + HE+ -> C+ + HE
    if(Tgas.LE.2d2) then
      k(49) = small + (8.58d-17&
          *Tgas**(0.757))
    end if

    !C + HE+ -> C+ + HE
    if(Tgas.GT.2d2 .and. Tgas.LE.2d3) then
      k(50) = small + (3.25d-17&
          *Tgas**(0.968))
    end if

    !C + HE+ -> C+ + HE
    if(Tgas.GT.2d3 .and. Tgas.LT.1d8) then
      k(51) = small + (2.77d-19&
          *Tgas**(1.597))
    end if

    !OH + H -> O + H + H
    k(52) = small + (6d-9*exp(-5.09d4&
        *invT))

    !HOC+ + H2 -> HCO+ + H2
    k(53) = small + (3.8d-10)

    !HOC+ + CO -> HCO+ + CO
    if(Tgas.LT.1d1) then
      k(54) = small + (1.604d-9)
    end if

    !HOC+ + CO -> HCO+ + CO
    if(Tgas.GE.1d1) then
      k(55) = small + (8.68d-10&
          *(1.+2.42717d-2*sqrt(3e2*invT)+7.1537*invT))
    end if

    !C + H2 -> CH + H
    k(56) = small + (6.64d-10*exp(-11700d0&
        *invT))

    !CH + H -> C + H2
    k(57) = small + (1.31d-10*exp(-8d1*invT))

    !CH + H2 -> CH2 + H
    k(58) = small + (5.46d-10*exp(-1943d0&
        *invT))

    !CH + C -> C2 + H
    k(59) = small + (2.40d-10)

    !CH + O -> CO + H
    k(60) = small + (1.02d-10*exp(-914d0&
        *invT))

    !CH + O -> HCO+ + E
    k(61) = small + (1.9d-11*(T32)**(-2.2)&
        *exp(-165.1d0*invT))

    !CH + O -> OH + C
    k(62) = small + (2.52d-11*exp(-2381d0&
        *invT))

    !CH2 + H -> CH + H2
    k(63) = small + (2.2d-10)

    !CH2 + O -> CO + H + H
    k(64) = small + (2.04d-10*exp(-270d0&
        *invT))

    !CH2 + O -> CO + H2
    k(65) = small + (1.36d-10*exp(-270d0&
        *invT))

    !CH2 + O -> HCO + H
    k(66) = small + (5.01d-11)

    !CH2 + O -> CH + OH
    k(67) = small + (4.98d-10*exp(-6d3&
        *invT))

    !C2 + O -> CO + C
    if(Tgas.LT.3d2) then
      k(68) = small + (2d-12&
          *(T32)**(-.12))
    end if

    !C2 + O -> CO + C
    if(Tgas.GE.3d2) then
      k(69) = small + (2d-12&
          *(T32)**(.757))
    end if

    !O + H2 -> OH + H
    k(70) = small + (1.46e-12*exp(-9650.&
        *invT))

    !OH + H -> O + H2
    if(Tgas.LT.280.d0) then
      k(71) = small + (6.99d-14&
          *T32**2.8*exp(-1950d0*invT))
    end if

    !OH + H -> O + H2
    if(Tgas.GE.280.d0) then
      k(72) = small + (5.45d-17)
    end if

    !H2 + OH -> H2O + H
    k(73) = small + (3.6d-16*T**(1.52)&
        *exp(-1.74d3*invT))

    !C + OH -> H + CO
    if(Tgas.LT.1d1) then
      k(74) = small + (7.051e-11)
    end if

    !C + OH -> H + CO
    if(Tgas.GE.1d1) then
      k(75) = small + (2.25d-11&
          *(T32)**(-.339)*exp(-.108d0*invT))
    end if

    !O + OH -> H + O2
    if(Tgas.GE.150.d0) then
      k(76) = small + (2.4d-11&
          *exp(110d0*invT))
    end if

    !O + OH -> H + O2
    if(Tgas.LT.150.d0) then
      k(77) = small + (4.997d-11)
    end if

    !OH + OH -> H2O + O
    k(78) = small + (1.65d-12*(T32)**1.14&
        *exp(-5d1*invT))

    !H2O + H -> H2 + OH
    k(79) = small + (1.59d-11*(T32)*1.2&
        *exp(-9610.*invT))

    !O2 + H -> OH + O
    k(80) = small + (2.61d-10*1.2*exp(-8156.&
        *invT))

    !O2 + H2 -> OH + OH
    k(81) = small + (3.16d-10*exp(-21890.d0&
        *invT))

    !O2 + C -> CO + O
    if(Tgas.LT.1052d0) then
      k(82) = small + (4.7d-11&
          *T32**(-.34))
    end if

    !O2 + C -> CO + O
    if(Tgas.GE.1052d0) then
      k(83) = small + (2.48d-12&
          *T32**1.54*exp(613d0*invT))
    end if

    !CO + H -> C + OH
    k(84) = small + (1.1d-10*T32**0.5&
        *exp(-77700d0*invT))

    !H2+ + H2 -> H3+ + H
    k(85) = small + (2.24d-9*T32**.042&
        *exp(-Tgas&
        /46600.))

    !H3+ + H -> H2+ + H2
    k(86) = small + (7.7d-9*exp(-17560d0&
        *invT))

    !C + H2+ -> CH+ + H
    k(87) = small + (2.4d-9)

    !C + H3+ -> CH+ + H2
    k(88) = small + ((1.0218d-9 + 7.2733d-11&
        *sqrt(Tgas) + 5.9203d-14*Tgas)&
        /(Tgas**0.1667 + 4.4914d-2&
        *sqrt(Tgas) - 5.9203d-14*Tgas + 2.6397d-6*Tgas**1.5))

    !C + H3+ -> CH2+ + H
    k(89) = small + ((8.5145d-10)&
        /(Tgas**(.1667) + 9.5666d-4&
        *sqrt(Tgas) - 4.404d-5*Tgas + 2.3496d-6 * Tgas**1.5))

    !C+ + H2 -> CH+ + H
    k(90) = small + (1d-10*exp(-4640d0&
        *invT))

    !CH+ + H -> C+ + H2
    k(91) = small + (7.5d-10)

    !CH+ + H2 -> CH2+ + H
    k(92) = small + (1.2d-9)

    !CH+ + O -> CO+ + H
    k(93) = small + (3.5d-10)

    !CH2+ + H -> CH+ + H2
    k(94) = small + (1d-9*exp(-7080d0&
        *invT))

    !CH2+ + H2 -> CH3+ + H
    k(95) = small + (1.6d-9)

    !CH2+ + O -> HCO+ + H
    k(96) = small + (7.5d-10)

    !CH3+ + H -> CH2+ + H2
    k(97) = small + (7.d-10*exp(-10560d0&
        *invT))

    !CH3+ + O -> HOC+ + H2
    k(98) = small + (2.5d-10)

    !CH3+ + O -> HCO+ + H2
    k(99) = small + (2.5d-10)

    !C2 + O+ -> CO+ + C
    k(100) = small + (4.8d-10)

    !O+ + H2 -> H + OH+
    k(101) = small + (1.69d-9)

    !O + H2+ -> H + OH+
    k(102) = small + (1.5d-9)

    !O + H3+ -> H2 + OH+
    k(103) = small + (7.98d-10&
        *T32**(-.156)*exp(-1.41d0*invT))

    !O + H3+ -> H + H2O+
    k(104) = small + (3.42d-10&
        *T32**(-.156)*exp(-1.41d0*invT))

    !OH + H3+ -> H2 + H2O+
    if(Tgas.LT.1d1) then
      k(105) = small + (2.277d-8)
    end if

    !OH + H3+ -> H2 + H2O+
    if(Tgas.GE.1d1) then
      k(106) = small + (1.52d-9&
          *(0.62d0 + 2.62185d0*(3d2*invT)**.5))
    end if

    !OH + C+ -> H + CO+
    if(Tgas.LT.1d1) then
      k(107) = small + (1.371d-08)
    end if

    !OH + C+ -> H + CO+
    if(Tgas.GE.1d1) then
      k(108) = small + (9.15d-10&
          *(0.62d0 + 2.62185d0*(3d2*invT)**.5))
    end if

    !OH+ + H2 -> H2O+ + H
    k(109) = small + (1.01d-9)

    !H2O+ + H2 -> H3O+ + H
    k(110) = small + (6.4d-10)

    !H2O + H3+ -> H2 + H3O+
    if(Tgas.LT.1d1) then
      k(111) = small + (2.55d-8)
    end if

    !H2O + H3+ -> H2 + H3O+
    if(Tgas.GE.1d1) then
      k(112) = small + (1.73d-9&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + C+ -> HOC+ + H
    k(113) = small + (1.8d-9)

    !H2O + C+ -> HCO+ + H
    if(Tgas.LT.1d1) then
      k(114) = small + (5.027d-9)
    end if

    !H2O + C+ -> HCO+ + H
    if(Tgas.GE.1d1) then
      k(115) = small + (3.4093e-10&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + C+ -> H2O+ + C
    k(116) = small + (2.4d-10)

    !H3O+ + C -> HCO+ + H2
    k(117) = small + (1d-11)

    !O2 + C+ -> CO+ + O
    k(118) = small + (3.42d-10)

    !O2 + C+ -> CO + O+
    k(119) = small + (4.53d-10)

    !O2 + CH2+ -> HCO+ + OH
    k(120) = small + (9.1d-10)

    !C + O2+ -> O + CO+
    k(121) = small + (5.2d-11)

    !C + O2+ -> O2 + C+
    k(122) = small + (5.2d-11)

    !CO + H3+ -> H2 + HCO+
    if(Tgas.LT.1d1) then
      k(123) = small + (2.468d-9)
    end if

    !CO + H3+ -> H2 + HCO+
    if(Tgas.GE.1d1) then
      k(124) = small + (1.88055d-9&
          *(1d0 + 0.02427d0*(3d2*invT)**.5 + 1.79558d0*invT))
    end if

    !CO + H3+ -> H2 + HOC+
    if(Tgas.LT.1d1) then
      k(125) = small + (1.421d-10)
    end if

    !CO + H3+ -> H2 + HOC+
    if(Tgas.GE.1d1) then
      k(126) = small + (1.08256d-10&
          *(1d0 + 0.02427d0*(3d2*invT)**.5 + 1.79558d0*invT))
    end if

    !HCO+ + C -> CO + CH+
    k(127) = small + (1.1d-9)

    !HCO+ + H2O -> CO + H3O+
    if(Tgas.LT.1d1) then
      k(128) = small + (7.279e-08)
    end if

    !HCO+ + H2O -> CO + H3O+
    if(Tgas.GE.1d1) then
      k(129) = small + (8.34d-10&
          *(1d0 + 0.5232d0*(3d2*invT)**.5 + 834.165880*invT))
    end if

    !CH + H+ -> CH+ + H
    if(Tgas.LT.1d1) then
      k(130) = small + (3.297d-8)
    end if

    !CH + H+ -> CH+ + H
    if(Tgas.GE.1d1) then
      k(131) = small + (3.54e-09&
          *(0.62d0 + 1.587411d0*(3d2*invT)**.5))
    end if

    !CH2 + H+ -> H2 + CH+
    if(Tgas.LT.1.5d2) then
      k(132) = small + (1.765d-9&
          *(0.62d0 + 0.672147d0*(3d2*invT)**.5))
    end if

    !CH2 + H+ -> H2 + CH+
    if(Tgas.GE.1.5d2) then
      k(133) = small + (1.765d-9&
          *(1d0 + 0.136347d0*(3d2*invT)**.5 + 56.66255d0*invT))
    end if

    !CH2 + H+ -> H + CH2+
    if(Tgas.LT.1.5d2) then
      k(134) = small + (1.765d-9&
          *(0.62d0 + 0.672147d0*(3d2*invT)**.5))
    end if

    !CH2 + H+ -> H + CH2+
    if(Tgas.GE.1.5d2) then
      k(135) = small + (1.765d-9&
          *(1d0 + 0.136347d0*(3d2*invT)**.5 + 56.66255d0*invT))
    end if

    !CH2 + HE+ -> HE + H2 + C+
    if(Tgas.LT.1.5d2) then
      k(136) = small + (9.65d-10&
          *(0.62d0 + 0.672147d0*(3d2*invT)**.5))
    end if

    !CH2 + HE+ -> HE + H2 + C+
    if(Tgas.GE.1.5d2) then
      k(137) = small + (9.65d-10&
          *(1d0 + 0.136347d0*(3d2*invT)**.5 + 56.6625498765d0&
          *invT))
    end if

    !CH2 + HE+ -> HE + H + CH+
    if(Tgas.LT.1.5d2) then
      k(138) = small + (9.65d-10&
          *(0.62d0 + 0.672147d0*(3d2*invT)**.5))
    end if

    !CH2 + HE+ -> HE + H + CH+
    if(Tgas.GE.1.5d2) then
      k(139) = small + (9.65d-10&
          *(1d0 + 0.136347d0*(3d2*invT)**.5 + 56.6625498765d0&
          *invT))
    end if

    !C2 + HE+ -> C+ + C + HE
    k(140) = small + (1.6d-9)

    !OH + H+ -> OH+ + H
    if(Tgas.LT.1d1) then
      k(141) = small + (3.745d-8)
    end if

    !OH + H+ -> OH+ + H
    if(Tgas.GE.1d1) then
      k(142) = small + (2.5d-9&
          *(0.62d0 + 2.62185d0*(3d2*invT)**.5))
    end if

    !OH + HE+ -> O+ + HE + H
    if(Tgas.LT.1d1) then
      k(143) = small + (2.022d-8)
    end if

    !OH + HE+ -> O+ + HE + H
    if(Tgas.GE.1d1) then
      k(144) = small + (1.35d-9&
          *(0.62d0 + 2.62185d0*(3d2*invT)**.5))
    end if

    !H2O + H+ -> H + H2O+
    if(Tgas.LT.1d1) then
      k(145) = small + (4.202d-8)
    end if

    !H2O + H+ -> H + H2O+
    if(Tgas.GE.1d1) then
      k(146) = small + (2.85d-9&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + HE+ -> HE + OH + H+
    if(Tgas.LT.1d1) then
      k(147) = small + (7.562d-9)
    end if

    !H2O + HE+ -> HE + OH + H+
    if(Tgas.GE.1d1) then
      k(148) = small + (5.1282d-10&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + HE+ -> HE + OH+ + H
    if(Tgas.LT.1d1) then
      k(149) = small + (7.562d-9)
    end if

    !H2O + HE+ -> HE + OH+ + H
    if(Tgas.GE.1d1) then
      k(150) = small + (5.1282d-10&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !H2O + HE+ -> HE + H2O+
    if(Tgas.LT.1d1) then
      k(151) = small + (7.56d-9)
    end if

    !H2O + HE+ -> HE + H2O+
    if(Tgas.GE.1d1) then
      k(152) = small + (5.1282d-10&
          *(0.62d0 + 2.578947d0*(3d2*invT)**.5))
    end if

    !O2 + H+ -> O2+ + H
    k(153) = small + (2d-9)

    !O2 + HE+ -> O2+ + HE
    k(154) = small + (3.3d-11)

    !O2 + HE+ -> O+ + HE + O
    k(155) = small + (1.1d-9)

    !CO + HE+ -> C+ + HE + O
    k(156) = small + (1.4d-9&
        *(T32)**(-.5))

    !CO + HE+ -> C + HE + O+
    k(157) = small + (1.4d-16&
        *(T32)**(-.5))

    !CO+ + H -> CO + H+
    k(158) = small + (7.5d-10)

    !C- + H+ -> C + H
    k(159) = small + (2.3d-7*(T32)**(-.5))

    !O- + H+ -> O + H
    k(160) = small + (2.3d-7*(T32)**(-.5))

    !HE+ + H- -> H + HE
    k(161) = small + (2.3d-7*T32**(-.5))

    !H3+ + E -> H2 + H
    k(162) = small + (2.34d-8*T32**(-.52))

    !H3+ + E -> H + H + H
    k(163) = small + (4.36d-8&
        *T32**(-.52))

    !CH+ + E -> C + H
    k(164) = small + (7d-8*T32**(-.5))

    !CH2+ + E -> CH + H
    k(165) = small + (1.6d-7*T32**(-.6))

    !CH2+ + E -> C + H2
    k(166) = small + (7.68d-8*T32**(-.6))

    !CH2+ + E -> C + H + H
    k(167) = small + (4.03d-7&
        *T32**(-.6))

    !CH3+ + E -> CH2 + H
    k(168) = small + (7.75d-8*T32**(-.5))

    !CH3+ + E -> CH + H2
    k(169) = small + (1.95d-7*T32**(-.5))

    !CH3+ + E -> CH + H + H
    k(170) = small + (2d-7*T32**(-.5))

    !OH+ + E -> O + H
    k(171) = small + (6.3d-9*T32**(-.48))

    !H2O+ + E -> O + H2
    k(172) = small + (3.9d-8*T32**(-.5))

    !H2O+ + E -> OH + H
    k(173) = small + (8.6d-8*T32**(-.5))

    !H2O+ + E -> O + H + H
    k(174) = small + (3.05d-7&
        *T32**(-.5))

    !H3O+ + E -> OH + H + H
    k(175) = small + (2.58d-7&
        *T32**(-.5))

    !H3O+ + E -> O + H + H2
    k(176) = small + (5.6d-9&
        *T32**(-.5))

    !H3O+ + E -> H + H2O
    k(177) = small + (1.08d-7*T32**(-.5))

    !H3O+ + E -> OH + H2
    k(178) = small + (6.02d-8*T32**(-.5))

    !O2+ + E -> O + O
    k(179) = small + (1.95d-7*T32**(-.7))

    !CO+ + E -> C + O
    k(180) = small + (2.75d-7*T32**(-.55))

    !HCO+ + E -> CO + H
    k(181) = small + (2.76d-7*T32**(-.64))

    !HCO+ + E -> OH + C
    k(182) = small + (2.4d-8*T32**(-.64))

    !HOC+ + E -> CO + H
    k(183) = small + (1.1d-7*invT32)

    !H- + C -> CH + E
    k(184) = small + (1d-9)

    !H- + O -> OH + E
    k(185) = small + (1d-10)

    !H- + OH -> H2O + E
    k(186) = small + (5d-10)

    !C- + H -> CH + E
    k(187) = small + (1d-13)

    !C- + H2 -> CH2 + E
    k(188) = small + (5d-10)

    !C- + O -> CO + E
    k(189) = small + (5d-10)

    !O- + H -> OH + E
    k(190) = small + (7d-10)

    !O- + H2 -> H2O + E
    k(191) = small + (7d-10)

    !O- + C -> CO + E
    k(192) = small + (5d-10)

    !H2 + H+ -> H + H + H+
    k(193) = small + (3d-11*T32**(.5)&
        *exp(-52000d0*invT))

    !H2 + H+ -> H3+
    k(194) = small + (1d-16)

    !C + E -> C-
    k(195) = small + (2.25d-15)

    !C + H -> CH
    k(196) = small + (1d-17)

    !C + H2 -> CH2
    k(197) = small + (1d-17)

    !C + C -> C2
    k(198) = small + (4.36d-18*T32**.35&
        *exp(-161.3d0*invT))

    !C + O -> CO
    k(199) = small + (3.09d-17*T32**.33&
        *exp(-1629d0*invT))

    !C+ + H -> CH+
    k(200) = small + (4.46d-16*Tgas**(-.5)&
        *exp(-4.93*Tgas**(-.6667)))

    !C+ + H2 -> CH2+
    k(201) = small + (2d-16*T32**(-1.3)&
        *exp(-23d0*invTgas))

    !C+ + O -> CO+
    if(Tgas.LT.3d2) then
      k(202) = small + (2.5d-18)
    end if

    !C+ + O -> CO+
    if(Tgas.GE.3d2) then
      k(203) = small + (3.14d-18&
          *T32**(-.15)*exp(-68d0*invT))
    end if

    !O + E -> O-
    k(204) = small + (1.5d-15)

    !O + H -> OH
    k(205) = small + (9.9d-19*T32**(-.38))

    !O + O -> O2
    k(206) = small + (4.9d-20*T32**(1.58))

    !OH + H -> H2O
    k(207) = small + (5.26d-18*T32**(-5.22)&
        *exp(-9d1*invT))

    !CO -> CO_ice
    k(208) = small + (krate_stickSi(n(:),idx_CO,tabTdust))

    !CO_ice -> CO
    k(209) = small + (krate_evaporation(n(:),idx_CO,tabTdust))

    !H2O -> H2O_ice
    k(210) = small + (krate_stickSi(n(:),idx_H2O,tabTdust))

    !H2O_ice -> H2O
    k(211) = small + (krate_evaporation(n(:),idx_H2O,tabTdust))

    !CO_ice -> CO
    k(212) = small + (krate_nonthermal_evaporation(idx_CO,1.78d0&
        *user_G0, user_Av, user_crate, 1d-3))

    !H -> H+ + E
    k(213) = small + (photoBinRates(1))

    !HE -> HE+ + E
    k(214) = small + (photoBinRates(2))

    !HE+ -> HE++ + E
    k(215) = small + (photoBinRates(3))

    !O -> O+ + E
    k(216) = small + (photoBinRates(4))

    !C -> C+ + E
    k(217) = small + (photoBinRates(5))

    !H2 -> H2+ + E
    k(218) = small + (photoBinRates(6))

    !H- -> H + E
    k(219) = small + (photoBinRates(7))

    !CH -> C + H
    k(220) = small + (photoBinRates(8))

    !CH -> CH+ + E
    k(221) = small + (photoBinRates(9))

    !C2 -> C + C
    k(222) = small + (photoBinRates(10))

    !OH -> O + H
    k(223) = small + (photoBinRates(11))

    !OH -> OH+ + E
    k(224) = small + (photoBinRates(12))

    !H2O -> OH + H
    k(225) = small + (photoBinRates(13))

    !H2O -> H2O+ + E
    k(226) = small + (photoBinRates(14))

    !O2 -> O2+ + E
    k(227) = small + (photoBinRates(15))

    !O2 -> O + O
    k(228) = small + (photoBinRates(16))

    !H2 -> H+ + H + E
    k(229) = small + (photoBinRates(17))

    !CO -> C + O
    k(230) = small + (user_gamma_CO)

    !H2 -> H + H
    k(231) = small + (user_gamma_H2)

    !H2+ -> H + H+
    k(232) = small + (user_G0*1.1d-9*exp(-1.9&
        *user_Av))

    !H3+ -> H2 + H+
    k(233) = small + (user_G0*4.9d-13*exp(-1.8&
        *user_Av))

    !H3+ -> H2+ + H
    k(234) = small + (user_G0*4.9d-13*exp(-2.3&
        *user_Av))

    !C- -> C + E
    k(235) = small + (user_G0*2.4d-7*exp(-.9&
        *user_Av))

    !CH+ -> C + H+
    k(236) = small + (user_G0*2.6d-10*exp(-2.5&
        *user_Av))

    !CH2 -> CH + H
    k(237) = small + (user_G0*7.1d-10*exp(-1.7&
        *user_Av))

    !CH2 -> CH2+ + E
    k(238) = small + (user_G0*5.9d-10*exp(-2.3&
        *user_Av))

    !CH2+ -> CH+ + H
    k(239) = small + (user_G0*4.6d-10*exp(-1.7&
        *user_Av))

    !CH3+ -> CH2+ + H
    k(240) = small + (user_G0*1d-9*exp(-1.7&
        *user_Av))

    !CH3+ -> CH+ + H2
    k(241) = small + (user_G0*1d-9*exp(-1.7&
        *user_Av))

    !O- -> O + E
    k(242) = small + (user_G0*2.4d-7*exp(-.5&
        *user_Av))

    !OH+ -> O + H+
    k(243) = small + (user_G0*1d-12*exp(-1.8&
        *user_Av))

    !H2O+ -> H2+ + O
    k(244) = small + (user_G0*5.d-11*HnOj)

    !H2O+ -> H+ + OH
    k(245) = small + (user_G0*5.d-11*HnOj)

    !H2O+ -> O+ + H2
    k(246) = small + (user_G0*5.d-11*HnOj)

    !H2O+ -> OH+ + H
    k(247) = small + (user_G0*1.5d-10*HnOj)

    !H3O+ -> H+ + H2O
    k(248) = small + (user_G0*2.5d-11*HnOj)

    !H3O+ -> H2+ + OH
    k(249) = small + (user_G0*2.5d-11*HnOj)

    !H3O+ -> H2O+ + H
    k(250) = small + (user_G0*7.5d-12*HnOj)

    !H3O+ -> OH+ + H2
    k(251) = small + (user_G0*2.5d-11*HnOj)

    !H -> H+ + E
    k(252) = rateEvaluateOnce(252)

    !HE -> HE+ + E
    k(253) = rateEvaluateOnce(253)

    !O -> O+ + E
    k(254) = rateEvaluateOnce(254)

    !CO -> C + O
    k(255) = rateEvaluateOnce(255)

    !CO -> CO+ + E
    k(256) = rateEvaluateOnce(256)

    !C2 -> C + C
    k(257) = rateEvaluateOnce(257)

    !H2 -> H + H
    k(258) = rateEvaluateOnce(258)

    !H2 -> H+ + H-
    k(259) = rateEvaluateOnce(259)

    !H2 -> H2+ + E
    k(260) = rateEvaluateOnce(260)

    !C -> C+ + E
    k(261) = rateEvaluateOnce(261)

    !CH -> C + H
    k(262) = rateEvaluateOnce(262)

    !O2 -> O + O
    k(263) = rateEvaluateOnce(263)

    !O2 -> O2+ + E
    k(264) = rateEvaluateOnce(264)

    !OH -> O + H
    k(265) = rateEvaluateOnce(265)

    !CH2 -> CH2+ + E
    k(266) = rateEvaluateOnce(266)

    !H2O -> OH + H
    k(267) = rateEvaluateOnce(267)

    !HCO -> CO + H
    k(268) = rateEvaluateOnce(268)

    !HCO -> HCO+ + E
    k(269) = rateEvaluateOnce(269)

    !H2 -> H + H+ + E
    k(270) = rateEvaluateOnce(270)

    !C + C -> C2
    if(Tgas.LT.5d3) then
      k(271) = small + (5.99d-33&
          *(Tgas&
          /5d3)**(-1.6)*ntot)
    end if

    !C + C -> C2
    if(Tgas.GE.5d3) then
      k(272) = small + (5.99d-33&
          *(Tgas&
          /5d3)**(-0.64)*exp(5255./Tgas)*ntot)
    end if

    !C + O -> CO
    if(Tgas.LT.2d3) then
      k(273) = small + (6.16d-29&
          *(Tgas&
          /3d2)**(-3.08)*ntot)
    end if

    !C + O -> CO
    if(Tgas.GE.2d3) then
      k(274) = small + (2.14d-29&
          *(Tgas&
          /3d2)**(-3.08)*exp(2114./Tgas)*ntot)
    end if

    !C+ + O -> CO+
    if(Tgas.LT.2d3) then
      k(275) = small + (6.16d-27&
          *(Tgas&
          /3d2)**(-3.08)*ntot)
    end if

    !C+ + O -> CO+
    if(Tgas.GE.2d3) then
      k(276) = small + (2.14d-27&
          *(Tgas&
          /3d2)**(-3.08)*exp(2114./Tgas)*ntot)
    end if

    !C + O+ -> CO+
    if(Tgas.LT.2d3) then
      k(277) = small + (6.16d-27&
          *(Tgas&
          /3d2)**(-3.08)*ntot)
    end if

    !C + O+ -> CO+
    if(Tgas.GE.2d3) then
      k(278) = small + (2.14d-27&
          *(Tgas&
          /3d2)**(-3.08)*exp(2114./Tgas)*ntot)
    end if

    !H + O -> OH
    k(279) = small + (4.33d-32*(T32)**(-1)*ntot)

    !OH + H -> H2O
    k(280) = small + (2.56d-31*(T32)**(-2)*ntot)

    !O + O -> O2
    k(281) = small + (9.2d-34*(T32)**(-1)*ntot)

    coe(:) = k(:) !set coefficients to return variable

    !!uncomment below to check coefficient values
    !kmax = 1d0
    !if(maxval(k)>kmax.or.minval(k)<0d0) then
    !   print *,"***************"
    !   do i=1,size(k)
    !      if(k(i)<0d0.or.k(i)>kmax) print *,i,k(i)
    !   end do
    !end if
  end function coe

  !*************************
  subroutine loadReactionsVerbatim()
    use krome_commons
    implicit none
    character*50::fname,line
    integer::ios,i,nunit

    fname = "reactions_verbatim.dat"

    !verbatim reactions are loaded from file
    ! to increase compilation speed
    open(newunit=nunit,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: "//trim(fname)//" file not present!"
      stop
    end if

    !load reactions from file
    do i=1,nrea
      read(nunit,'(a)',iostat=ios) line
      if(ios/=0) then
        print *,"ERROR: problem reading "//trim(fname)
        stop
      end if
      reactionNames(i) = trim(line)
    end do
    close(nunit)

  end subroutine loadReactionsVerbatim

  !*******************
  !The following functions compute the recombination rate
  ! on dust for H+, He+, C+, Si+, and O+. See Weingartner&Draine 2001
  ! dust2gas_ratio, D/D_sol, default is assumed equal to Z/Z_sol
  function H_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::H_recombination_on_dust

    H_recombination_on_dust = 0d0

    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    H_recombination_on_dust =  1.225d-13*dust2gas_ratio &
        /(1.d0+8.074d-6*psi**(1.378)*(1.d0+5.087d2 &
        *Tgas**(0.01586)*psi**(-0.4723-1.102d-5*log(Tgas))))

  end function H_recombination_on_dust

  !******************
  function He_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::He_recombination_on_dust

    He_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    He_recombination_on_dust = 5.572d-14*dust2gas_ratio&
        /(1.d0+3.185d-7*psi**(1.512)*(1.d0+5.115d3&
        *Tgas**(3.903d-7)*psi**(-0.4956-5.494d-7*log(Tgas))))

  end function He_recombination_on_dust

  !*******************
  function C_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::C_recombination_on_dust

    C_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    C_recombination_on_dust = 4.558d-13*dust2gas_ratio&
        /(1.d0+6.089d-3*psi**(1.128)*(1.d0+4.331d2&
        *Tgas**(0.04845)*psi**(-0.8120-1.333d-4*log(Tgas))))

  end function C_recombination_on_dust

  !******************
  function Si_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::Si_recombination_on_dust

    Si_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    Si_recombination_on_dust = 2.166d-14*dust2gas_ratio&
        /(1.d0+5.678d-8*psi**(1.874)*(1.d0+4.375d4&
        *Tgas**(1.635d-6)*psi**(-0.8964-7.538d-5*log(Tgas))))

  end function Si_recombination_on_dust

  !********************
  function O_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k_H
    real*8::O_recombination_on_dust

    k_H = H_recombination_on_dust(n(:),Tgas)
    O_recombination_on_dust = 0.25d0*k_H

  end function O_recombination_on_dust

  !*********************
  !This function returns the
  ! photorate of H2 occurring in the
  ! Lyman-Werner bands following the approximation
  ! provided by Glover&Jappsen 2007. Rate in 1/s.
  !Approximation valid at low-density, it assumes H2(nu = 0).
  !It also stores the rate as a common, needed for the photoheating
  function H2_solomonLW(myflux)
    use krome_commons
    use krome_constants
    implicit none
    real*8::H2_solomonLW,myflux

    !myflux is the radiation background at E = 12.87 eV
    !should be converted to erg
    H2_solomonLW = 1.38d9*myflux*eV_to_erg

  end function H2_solomonLW

  !****************************
  !tanh smoothing function that
  ! increses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_increase(xarg,xpos,slope)
    implicit none
    real*8::smooth_increase,xarg,xpos,slope

    smooth_increase = .5d0 * (tanh(slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_increase

  !****************************
  !tanh smoothing function that
  ! decreses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_decrease(xarg,xpos,slope)
    implicit none
    real*8::smooth_decrease,xarg,xpos,slope

    smooth_decrease = .5d0 * (tanh(-slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_decrease

  !*********************
  !sign: return 1d0 if x>=0d0,
  ! else return -1d0
  function get_sgn(x)
    implicit none
    real*8::x,get_sgn

    get_sgn = 1d0
    if(x==0d0) return
    get_sgn = x/abs(x)

  end function get_sgn

  !*********************
  function conserve(n,ni)
    use krome_commons
    implicit none
    real*8::conserve(nspec),n(nspec),ni(nspec),no(nspec)
    real*8::ntot,nitot,factor

    no(:) = n(:)

    conserve(:) = 0d0
    conserve(:) = no(:)

  end function conserve

  !*************************
  !this subroutine changes the x(:) mass fractions of the species
  ! to force conservation according to the reference ref(:)
  subroutine conserveLin_x(x,ref)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::x(nmols),ref(natoms)
    real*8::A(natoms,natoms),B(natoms),m(nspec)

    m(:) = get_mass()
    A(:,:) = 0d0

    B(:) = ref(:)

    !charge conservation
    x(idx_E) = m(idx_E)*(- 1d0*x(idx_Hk) / m(idx_Hk) &
        - 1d0*x(idx_Ck) / m(idx_Ck) &
        - 1d0*x(idx_Ok) / m(idx_Ok) &
        + 1d0*x(idx_Hj) / m(idx_Hj) &
        + 1d0*x(idx_HEj) / m(idx_HEj) &
        + 1d0*x(idx_H2j) / m(idx_H2j) &
        + 1d0*x(idx_Cj) / m(idx_Cj) &
        + 1d0*x(idx_Oj) / m(idx_Oj) &
        + 1d0*x(idx_HOCj) / m(idx_HOCj) &
        + 1d0*x(idx_HCOj) / m(idx_HCOj) &
        + 1d0*x(idx_H3j) / m(idx_H3j) &
        + 1d0*x(idx_CHj) / m(idx_CHj) &
        + 1d0*x(idx_CH2j) / m(idx_CH2j) &
        + 1d0*x(idx_COj) / m(idx_COj) &
        + 1d0*x(idx_CH3j) / m(idx_CH3j) &
        + 1d0*x(idx_OHj) / m(idx_OHj) &
        + 1d0*x(idx_H2Oj) / m(idx_H2Oj) &
        + 1d0*x(idx_H3Oj) / m(idx_H3Oj) &
        + 1d0*x(idx_O2j) / m(idx_O2j) &
        + 2d0*x(idx_HEjj) / m(idx_HEjj))
    !check if charge conservation goes wrong
    if(x(idx_E)<0d0) then
      print *,"ERROR in conserveLin, electrons < 0"
      stop
    end if

  end subroutine conserveLin_x

  !***************************
  !compute the total reference mass atom type by atom type
  function conserveLinGetRef_x(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::conserveLinGetRef_x(natoms),x(nmols)
    real*8::m(nspec)

    m(:) = get_mass()
    conserveLinGetRef_x(:) = 0d0

  end function conserveLinGetRef_x

  !***************************
  !Ref: Sasaki & Takahara (1993)
  !This function evaluate the recombination rate
  ! for H+ + e --> H + gamma and the same
  ! for D+ + e --> D + gamma
  function elec_recomb_ST93(nabund,nelec,ntot,nucleiH,Trad)
    use krome_commons
    use krome_constants
    implicit none
    real*8::nabund,nelec,Trad
    real*8::nucleiH,elec_recomb_ST93
    real*8::al,ak,rc2,r2c
    real*8::a0,b0,c0,d0,e0
    real*8::a1,b1,c1,d1,e1,f1,g1,h1
    real*8::ntot,ratio

    al = 8.227d0
    ak = 22.06d0 / (hubble  *(1d0 + phys_zredshift) &
        * sqrt(1d0 + Omega0 * phys_zredshift))
    !Rc2 evaluation
    rc2 = 8.76d-11 * (1d0 + phys_zredshift)**(-0.58)
    !R2c evaluation
    r2c = (1.80d10 * Trad)**(1.5) &
        * exp(-3.9472d4 / Trad) * rc2

    !coefficients
    a0 = nucleiH * rc2
    b0 = ak * al * nucleiH
    c0 = ak * rc2 * nucleiH * nucleiH
    d0 = r2c * exp(-1.18416d5/Trad)
    e0 = ak * r2c * nucleiH

    !polynomial terms
    a1 = -d0 * (1d0 + b0)
    b1 = d0 * (1d0 + 2d0 * b0)
    c1 = a0 + b0 * (a0 - d0)
    d1 = -a0 * b0
    e1 = a0 * c0
    f1 = 1d0 + b0 + e0
    g1 = -(b0 + e0)
    h1 = c0

    ratio = nabund / ntot

    elec_recomb_ST93 = ntot*(a1 + b1*ratio + c1*ratio**2 + d1*ratio**3 &
        + e1*ratio**4) / (f1 + g1*ratio + h1*ratio**2)

    elec_recomb_ST93 = elec_recomb_ST93 / (nabund * nelec)

  end function elec_recomb_ST93

  !********************
  subroutine load_parts()
    use krome_commons
    implicit none

  end subroutine load_parts

  !*************************
  subroutine load_part(fname,array_part,min_part,dT_part)
    character(len=*)::fname
    integer::ios,icount,i,cv
    real*8,allocatable::array_part(:),emed(:)
    real*8::min_part,dT_part,Told,array_tmp(int(1e5)),rout(2)

    open(33,file=trim(fname),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: partition function not found"
      print *," in file "//fname
      stop
    end if

    print *,"loading partition function from "//fname
    icount = 0
    min_part = 1d99
    Told = 0d0
    do
      read(33,*,iostat=ios) rout(:)
      if(ios<0) exit
      if(ios.ne.0) cycle
      icount = icount + 1
      min_part = min(min_part,rout(1))
      array_tmp(icount) = rout(2)
      dT_part = rout(1) - Told
      Told = rout(1)
    end do
    close(33)

    allocate(array_part(icount),emed(icount))
    array_part(:) = array_tmp(1:icount)

  end subroutine load_part

  !**********************
  function troe_falloff(k0,kinf,Fc,m)
    implicit none
    real*8::troe_falloff,k0,kinf,Fc,m,rm,xexp
    rm = k0*m/kinf
    xexp = 1d0/(1d0+log10(rm)**2)
    troe_falloff = k0*m/(1d0+rm)*Fc**xexp
  end function troe_falloff

  !*************************
  function k3body(k0,kinf,Fc,nM)
    implicit none
    real*8::k3body,k0,kinf,Fc,nM
    real*8::c,n,d,Pr,xexp,F

    c = -0.4d0-0.67d0*log10(Fc)
    n = 0.75d0-1.27d0*log10(Fc)
    d = 0.14d0
    Pr = k0*nM/kinf
    xexp = (log10(Pr)+c)/(n-d*(log10(Pr)+c))
    F = 1d1**(log10(Fc)/(1d0+xexp**2))
    k3body = kinf*(Pr/(1d0+Pr)) * F

  end function k3body

  !***********************
  !see http://kida.obs.u-bordeaux1.fr/help
  function KIDA3body(ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,&
        kcFc,kdFc,npart,Tgas,pmin,pmax)
    implicit none
    real*8::ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,kcFc,kdFc
    real*8::KIDA3body,kinf,p,f,npart,Tgas,fc,fexp,invT
    real*8::k0,cc,dd,nn,pmin,pmax

    KIDA3body = 0d0

    invT = 1d0/Tgas
    k0 = ka0*(Tgas/3d2)**kb0*exp(-kc0*invT)
    kinf = kainf*(Tgas/3d2)**kbinf*exp(-kcinf*invT)

    p = k0*npart/kinf
    if(p<pmin) return
    if(p>pmax) return

    fc = (1d0-kaFc)*exp(-Tgas/kbFc) + kaFc*exp(-Tgas/kbFc) &
        + exp(-kdFc*invT)

    cc = -0.4d0 - 0.67d0 *log10(fc)
    dd = 0.14d0
    nn = 0.75d0 - 1.27d0*log10(fc)
    fexp = 1d0 + ((log10(p)+cc)/(nn-dd*(log10(p)+cc)))**2

    f = fc**(1d0/fexp)

    KIDA3body = kinf*(p/(1d0+p))*f

  end function KIDA3body

  !******************************
  !collisional ionization rate from Verner+96
  ! unit: cm3/s
  function colion_v96(Tgas,dE,P,A,X,K)
    implicit none
    real*8::colion_v96,Tgas,dE,A,X,K,U,Te,P

    Te = Tgas * 8.621738d-5 !K to eV
    U = dE / Te
    colion_v96 = A * (1d0 + P*sqrt(U)) * U**K * exp(-U) / (X+U)

  end function colion_v96

  !****************************
  !radiative recombination rates from
  ! Verner routine, standard fit, cm3/s
  function recV96(Tgas,a,b)
    implicit none
    real*8::recV96,Tgas,a,b

    recV96 = a*(1d4/Tgas)**b

  end function recV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, new fit, cm3/s
  function recNewV96(Tgas,r1,r2,r3,r4)
    implicit none
    real*8::recNewV96,Tgas,r1,r2,r3,r4,tt

    tt = sqrt(Tgas/r3)
    recNewV96 = r1/(tt*(tt + 1d0)**(1.-r2) &
        * (1d0 + sqrt(Tgas/r4))**(1.+r2))

  end function recNewV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, iron only, cm3/s
  function recFeV96(Tgas,r1,r2,r3)
    implicit none
    real*8::recFeV96,Tgas,r1,r2,r3,tt

    tt = sqrt(Tgas*1d-4)
    recFeV96 = r1/tt**(r2 + r3 + log10(tt))

  end function recFeV96

  !******************************
  !radiative recombination rates from Verner+96
  ! unit: cm3/s
  function radrec_v96(Tgas,a,b,T0,T1)
    implicit none
    real*8::Tgas,a,b,T0,T1,radrec_v96,iT0

    iT0 = 1d0/T0
    radrec_v96 = a/(sqrt(Tgas*iT0) + (1d0*sqrt(Tgas*iT0))**(1.-b) &
        * (1d0+sqrt(Tgas/T1))**(1+b))

  end function radrec_v96

  !*******************************
  !radiative recombination rates low-temp fit, Verner+96
  ! unit: cm3/s
  function radrec_low_v96(Tgas,a,b,c,d,f)
    implicit none
    real*8::Tgas,a,b,c,d,f,radrec_low_v96,t,invt

    t = Tgas*1d-4
    invt = 1d0/t

    radrec_low_v96 = 1d-12 * (a*invt + b + c*t + d*t**2) &
        * t**(-1.5) * exp(-f*invt)

    radrec_low_v96 = max(0d0,radrec_low_v96)

  end function radrec_low_v96

  !***************************
  !Collisional dissociation rate (cm-3/s) by Martin et al. 1996
  ! H2+H->H+H+H
  !NOTE: the use of this rate is suggested
  ! for high-density regime and in the presence of UV backgrounds.
  ! if necessary it must be included in the reaction file as
  ! H2,H,,H,H,H,,NONE,NONE,dissH2_Martin96(n,Tgas)
  function dissH2_Martin96(n,Tgas)
    use krome_commons
    use krome_getphys
    integer::i
    real*8::n(nspec),Tgas,dissH2_Martin96
    real*8::CDrates,logTv(4),k_CIDm(21,2),k_CID,invT,logT,n_c1,n_c2,n_H
    real*8::logk_h1,logk_h2,logk_l1,logk_l2,logn_c1,logn_c2,p,logk_CID
    real*8::logT2,logT3

    !k_CID = collision-induced dissociation + dissociative tunneling

    !Collisional dissociation of H2
    k_CIDm(:,1) = (/-178.4239d0, -68.42243d0, 43.20243d0, -4.633167d0, &
        69.70086d0, 40870.38d0, -23705.70d0, 128.8953d0, -53.91334d0, &
        5.315517d0, -19.73427d0, 16780.95d0, -25786.11d0, 14.82123d0, &
        -4.890915d0, 0.4749030d0, -133.8283d0, -1.164408d0, 0.8227443d0,&
        0.5864073d0, -2.056313d0/)

    !Dissociative tunneling of H2
    k_CIDm(:,2) = (/-142.7664d0, 42.70741d0, -2.027365d0, -0.2582097d0, &
        21.36094d0, 27535.31d0, -21467.79d0, 60.34928d0, -27.43096d0, &
        2.676150d0, -11.28215d0, 14254.55d0, -23125.20d0, 9.305564d0, &
        -2.464009d0, 0.1985955d0, 743.0600d0, -1.174242d0, 0.7502286d0, &
        0.2358848d0, 2.937507d0/)

    n_H  = get_Hnuclei(n(:))
    logT = log10(Tgas)
    invT = 1.0d0/Tgas
    logT2 = logT*logT
    logT3 = logT2*logT
    logTv = (/1.d0, logT, logT2, logT3/)
    k_CID = 0.d0
    do i=1,2
      logk_h1 = k_CIDm(1,i)*logTv(1) + k_CIDm(2,i)*logTv(2) + &
          k_CIDm(3,i)*logTv(3) + k_CIDm(4,i)*logTv(4) + &
          k_CIDm(5,i)*log10(1.d0+k_CIDm(6,i)*invT)
      logk_h2 = k_CIDm(7,i)*invT
      logk_l1 = k_CIDm(8,i)*logTv(1) + k_CIDm(9,i)*logTv(2) + &
          k_CIDm(10,i)*logTv(3) + k_CIDm(11,i)*log10(1.d0+k_CIDm(12,i)*invT)
      logk_l2 = k_CIDm(13,i)*invT
      logn_c1 = k_CIDm(14,i)*logTv(1) + k_CIDm(15,i)*logTv(2) &
          + k_CIDm(16,i)*logTv(3) + k_CIDm(17,i)*invT
      logn_c2 = k_CIDm(18,i) + logn_c1
      p = k_CIDm(19,i) + k_CIDm(20,i)*exp(-Tgas/1.850d3) &
          + k_CIDm(21,i)*exp(-Tgas/4.40d2)
      n_c1 = 1d1**(logn_c1)
      n_c2 = 1d1**(logn_c2)
      logk_CID = logk_h1 - (logk_h1 - logk_l1) / (1.d0 + (n_H/n_c1)**p) &
          + logk_h2 - (logk_h2 - logk_l2) / (1.d0 + (n_H/n_c2)**p)
      k_CID = k_CID + 1.d1**logk_CID
    enddo

    dissH2_Martin96 = k_CID

  end function dissH2_Martin96

  !***********************************
  subroutine init_exp_table()
    use krome_commons
    implicit none
    integer::i
    real*8::a

    do i=1,exp_table_na
      a = (i-1)*(exp_table_aMax-exp_table_aMin)/(exp_table_na-1) + exp_table_aMin
      exp_table(i) = exp(-a)
    end do

  end subroutine init_exp_table

  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    use krome_getphys
    integer::get_index,i
    character*16::names(nspec)
    character*(*)::name
    names(:) = get_names()
    get_index = -1 !default index
    !loop on species to found the specie named name
    do i=1,nspec
      !when found store and break loop
      if(trim(names(i))== trim(name)) then
        get_index = i !store index
        exit
      end if
    end do

    !error if species not found
    if(get_index<0) then
      print *,"ERROR: can't find the index of ",name
      stop
    end if

  end function get_index

  !*****************************
  !computes revers kinetics from reaction and
  ! product indexes
  ! k_rev = k_for * revKc
  ! Note that reaction constant revKc is calculated with
  ! reactants and products from reverse reaction
  function revKc(Tgas,ridx,pidx)
    use krome_constants
    use krome_commons
    implicit none
    real*8::revKc,Tgas,dgibss,stoichiometricChange
    integer::ridx(:),pidx(:),i

    ! when considering forward reaction:
    ! Kc = (P°)**(p+p-r-r) * exp(-dGibss_forward°)
    ! where ° means at standard conditions of
    ! P° = 1 bar = (kb*T/1e6) dyn/cm^2 (cgs)
    ! when considering reverse:
    ! 1/Kc = revKc = (kb*T/1e6)**(p+p-r-r) * exp(-dGibss_reverse°)
    ! kb*T/1e6 is to go from 1 atm pressure to number density cm^-3
    ! When not at standard pressure this does not change:
    ! revKc = P**(p+p-r-r) *exp(-dGibss_reverse° - (p+p-r-r)*ln(P/P°))
    !       = (P°)**(p+p-r-r) * exp(-dGibss_reverse°)

    dgibss = 0.d0 ! Gibbs free energy/(R*T)
    stoichiometricChange = 0d0

    do i=1,size(pidx)
      dgibss = dgibss + revHS(Tgas,pidx(i))
      stoichiometricChange = stoichiometricChange + 1
    end do

    do i=1,size(ridx)
      dgibss = dgibss - revHS(Tgas,ridx(i))
      stoichiometricChange = stoichiometricChange - 1
    end do

    revKc = (boltzmann_erg * Tgas * 1e-6)**(-stoichiometricChange)&
        * exp(-dgibss)

  end function revKc

  !*****************************
  !compute H-S for species with index idx
  ! when temperature is Tgas
  function revHS(Tgas,idx)
    use krome_commons
    use krome_constants
    use krome_fit
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
    real*8::Tnist,Tnist2,Tnist3,Tnist4,invTnist,invTnist2,lnTnist
    real*8::p1_nasa(40,7), p2_nasa(40,7), Tlim_nasa(40,3), p(7)
    real*8::p1_nist(40,7), p2_nist(40,7), Tlim_nist(40,3)
    integer::idx

    p(:) = 0.d0
    p1_nasa(:,:) = 0.d0
    p2_nasa(:,:) = 0.d0
    Tlim_nasa(:,:) = 0.d0
    p1_nist(:,:) = 0.d0
    p2_nist(:,:) = 0.d0
    Tlim_nist(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)
    ! NIST polynomials are quite differernt
    ! it doesn't like easy stuff...
    Tnist = Tgas * 1.d-3
    Tnist2 = Tnist * Tnist
    Tnist3 = Tnist2 * Tnist
    Tnist4 = Tnist3 * Tnist2
    invTnist = 1d0/Tnist
    invTnist2 = invTnist * invTnist
    lnTnist = log(Tnist)

    p1_nasa(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p1_nasa(idx_Ck,:)  = (/2.50025151d0,&
        -1.19774349d-06,&
        2.28919443d-09,&
        -1.98276803d-12,&
        6.44398056d-16,&
        70064.893d0,&
        4.87847086d0/)
    p1_nasa(idx_Ok,:)  = (/2.90805921d0,&
        -0.00169804907d0,&
        2.98069955d-06,&
        -2.43835127d-09,&
        7.61229311d-13,&
        11435.7717d0,&
        2.80339097d0/)
    p1_nasa(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p1_nasa(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p1_nasa(idx_H2,:)  = (/2.34433112d0,&
        0.00798052075d0,&
        -1.9478151d-05,&
        2.01572094d-08,&
        -7.37611761d-12,&
        -917.935173d0,&
        0.683010238d0/)
    p1_nasa(idx_C,:)  = (/2.5542395d0,&
        -0.00032153772d0,&
        7.3379223d-07,&
        -7.3223487d-10,&
        2.6652144d-13,&
        85442.681d0,&
        4.5313085d0/)
    p1_nasa(idx_O,:)  = (/3.1682671d0,&
        -0.00327931884d0,&
        6.64306396d-06,&
        -6.12806624d-09,&
        2.11265971d-12,&
        29122.2592d0,&
        2.05193346d0/)
    p1_nasa(idx_OH,:)  = (/3.99198424d0,&
        -0.00240106655d0,&
        4.61664033d-06,&
        -3.87916306d-09,&
        1.36319502d-12,&
        3368.89836d0,&
        -0.103998477d0/)
    p1_nasa(idx_CO,:)  = (/3.5795335d0,&
        -0.00061035369d0,&
        1.0168143d-06,&
        9.0700586d-10,&
        -9.0442449d-13,&
        -14344.086d0,&
        3.5084093d0/)
    p1_nasa(idx_CH,:)  = (/3.4897583d0,&
        0.0003243216d0,&
        -1.6899751d-06,&
        3.162842d-09,&
        -1.4061803d-12,&
        70660.755d0,&
        2.0842841d0/)
    p1_nasa(idx_CH2,:)  = (/3.84261832d0,&
        -7.36676871d-06,&
        6.16970693d-06,&
        -6.96689962d-09,&
        2.64620979d-12,&
        45863.1528d0,&
        1.2758447d0/)
    p1_nasa(idx_HCO,:)  = (/4.36380907d0,&
        -0.00535204137d0,&
        2.31954508d-05,&
        -2.6610904d-08,&
        1.02711962d-11,&
        25010.8717d0,&
        2.98106307d0/)
    p1_nasa(idx_H2O,:)  = (/4.1986352d0,&
        -0.0020364017d0,&
        6.5203416d-06,&
        -5.4879269d-09,&
        1.771968d-12,&
        -30293.726d0,&
        -0.84900901d0/)
    p1_nasa(idx_O2,:)  = (/3.78245636d0,&
        -0.00299673416d0,&
        9.84730201d-06,&
        -9.68129509d-09,&
        3.24372837d-12,&
        -1063.94356d0,&
        3.65767573d0/)
    p1_nasa(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p1_nasa(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p1_nasa(idx_H2j,:)  = (/3.77256072d0,&
        -0.0019574659d0,&
        4.54812047d-06,&
        -2.82152141d-09,&
        5.33969209d-13,&
        178694.654d0,&
        -3.96609192d0/)
    p1_nasa(idx_Cj,:)  = (/2.61332254d0,&
        -0.000540148065d0,&
        1.03037233d-06,&
        -8.90092552d-10,&
        2.88500586d-13,&
        216862.274d0,&
        3.8345479d0/)
    p1_nasa(idx_Oj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        187935.284d0,&
        4.39337676d0/)
    p1_nasa(idx_H3j,:)  = (/4.1795698d0,&
        -0.000868875627d0,&
        -1.09017371d-07,&
        4.13349766d-09,&
        -2.37877027d-12,&
        132635.537d0,&
        -5.838001d0/)
    p1_nasa(idx_COj,:)  = (/3.77061642d0,&
        -0.00201773246d0,&
        4.61081738d-06,&
        -2.99175463d-09,&
        6.06065045d-13,&
        149006.795d0,&
        3.38129783d0/)
    p1_nasa(idx_OHj,:)  = (/3.50502572d0,&
        0.000241313747d0,&
        -1.42200948d-06,&
        2.64780232d-09,&
        -1.17038711d-12,&
        155210.676d0,&
        1.97907627d0/)
    p1_nasa(idx_H2Oj,:)  = (/4.02465912d0,&
        -0.00108851414d0,&
        5.13576558d-06,&
        -4.40027838d-09,&
        1.40726746d-12,&
        116895.616d0,&
        0.699968812d0/)
    p1_nasa(idx_H3Oj,:)  = (/3.79295251d0,&
        -0.000910852723d0,&
        1.16363521d-05,&
        -1.21364865d-08,&
        4.26159624d-12,&
        71402.7518d0,&
        1.47156927d0/)
    p1_nasa(idx_O2j,:)  = (/4.61017167d0,&
        -0.00635951952d0,&
        1.42425624d-05,&
        -1.20997923d-08,&
        3.70956878d-12,&
        139742.229d0,&
        -0.201326941d0/)
    p2_nasa(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p2_nasa(idx_Ck,:)  = (/2.50001597d0,&
        -1.71721376d-08,&
        6.9283294d-12,&
        -1.20607892d-15,&
        7.60308635d-20,&
        70064.9324d0,&
        4.87955907d0/)
    p2_nasa(idx_Ok,:)  = (/2.54474869d0,&
        -4.66695513d-05,&
        1.84912357d-08,&
        -3.18159223d-12,&
        1.98962956d-16,&
        11504.2089d0,&
        4.52131015d0/)
    p2_nasa(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p2_nasa(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p2_nasa(idx_H2,:)  = (/2.93286575d0,&
        0.000826608026d0,&
        -1.46402364d-07,&
        1.54100414d-11,&
        -6.888048d-16,&
        -813.065581d0,&
        -1.02432865d0/)
    p2_nasa(idx_C,:)  = (/2.605583d0,&
        -0.00019593434d0,&
        1.0673722d-07,&
        -1.642394d-11,&
        8.187058d-16,&
        85411.742d0,&
        4.1923868d0/)
    p2_nasa(idx_O,:)  = (/2.54363697d0,&
        -2.73162486d-05,&
        -4.1902952d-09,&
        4.95481845d-12,&
        -4.79553694d-16,&
        29226.012d0,&
        4.92229457d0/)
    p2_nasa(idx_OH,:)  = (/2.83853033d0,&
        0.00110741289d0,&
        -2.94000209d-07,&
        4.20698729d-11,&
        -2.4228989d-15,&
        3697.80808d0,&
        5.84494652d0/)
    p2_nasa(idx_CO,:)  = (/3.0484859d0,&
        0.0013517281d0,&
        -4.8579405d-07,&
        7.8853644d-11,&
        -4.6980746d-15,&
        -14266.117d0,&
        6.0170977d0/)
    p2_nasa(idx_CH,:)  = (/2.5209369d0,&
        0.0017653639d0,&
        -4.614766d-07,&
        5.9289675d-11,&
        -3.3474501d-15,&
        70994.878d0,&
        7.4051829d0/)
    p2_nasa(idx_CH2,:)  = (/3.11049513d0,&
        0.00373779517d0,&
        -1.37371977d-06,&
        2.23054839d-10,&
        -1.33567178d-14,&
        45971.5953d0,&
        4.62796405d0/)
    p2_nasa(idx_HCO,:)  = (/4.23892214d0,&
        0.0019657617d0,&
        -3.82075171d-07,&
        4.80137647d-11,&
        -3.11176347d-15,&
        24726.1645d0,&
        1.99698242d0/)
    p2_nasa(idx_H2O,:)  = (/2.6770389d0,&
        0.0029731816d0,&
        -7.7376889d-07,&
        9.4433514d-11,&
        -4.2689991d-15,&
        -29885.894d0,&
        6.88255d0/)
    p2_nasa(idx_O2,:)  = (/3.66096065d0,&
        0.000656365811d0,&
        -1.41149627d-07,&
        2.05797935d-11,&
        -1.29913436d-15,&
        -1215.97718d0,&
        3.41536279d0/)
    p2_nasa(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p2_nasa(idx_HEj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        285323.374d0,&
        1.62166556d0/)
    p2_nasa(idx_H2j,:)  = (/3.44204765d0,&
        0.000599083239d0,&
        6.69133685d-08,&
        -3.43574373d-11,&
        1.97626599d-15,&
        178650.236d0,&
        -2.79499055d0/)
    p2_nasa(idx_Cj,:)  = (/2.50827618d0,&
        -1.04354146d-05,&
        5.16160809d-09,&
        -1.14187475d-12,&
        9.43539946d-17,&
        216879.645d0,&
        4.3188599d0/)
    p2_nasa(idx_Oj,:)  = (/2.48542028d0,&
        2.56978695d-05,&
        -1.28833378d-08,&
        1.65525487d-12,&
        1.09933344d-16,&
        187940.874d0,&
        4.47425446d0/)
    p2_nasa(idx_H3j,:)  = (/2.01435718d0,&
        0.00415925769d0,&
        -1.42664877d-06,&
        2.22372739d-10,&
        -1.29346518d-14,&
        133230.507d0,&
        5.46168967d0/)
    p2_nasa(idx_COj,:)  = (/2.93062935d0,&
        0.00156033262d0,&
        -6.16246355d-07,&
        1.09957336d-10,&
        -6.66119284d-15,&
        149147.222d0,&
        7.3384673d0/)
    p2_nasa(idx_OHj,:)  = (/2.68358996d0,&
        0.00157006435d0,&
        -5.39972815d-07,&
        9.37643877d-11,&
        -5.70068067d-15,&
        155479.296d0,&
        6.44375894d0/)
    p2_nasa(idx_H2Oj,:)  = (/3.31570445d0,&
        0.00210648746d0,&
        -3.76341515d-07,&
        3.47525972d-11,&
        -1.70335643d-15,&
        117017.475d0,&
        4.03220514d0/)
    p2_nasa(idx_H3Oj,:)  = (/2.49647765d0,&
        0.0057284484d0,&
        -1.83953239d-06,&
        2.73577348d-10,&
        -1.54093917d-14,&
        71624.4227d0,&
        7.45850493d0/)
    p2_nasa(idx_O2j,:)  = (/3.31675922d0,&
        0.00111522244d0,&
        -3.83492556d-07,&
        5.72784687d-11,&
        -2.77648381d-15,&
        139876.823d0,&
        5.44726469d0/)
    Tlim_nasa(idx_Hk,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Ck,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Ok,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HE,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_C,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_OH,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CH,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_CH2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HCO,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2O,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Hj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HEj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Cj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H3j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_COj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_OHj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H3Oj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_O2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)

    ! pick NASA data if present for species
    if (Tlim_nasa(idx,2) /= 0.d0) then
      !select set of NASA polynomials using temperature
      if(Tlim_nasa(idx,1).le.Tgas .and. Tgas.le.Tlim_nasa(idx,2)) then
        p(:) = p1_nasa(idx,:)

      else if(Tlim_nasa(idx,2)<Tgas .and. Tgas.le.Tlim_nasa(idx,3)) then
        p(:) = p2_nasa(idx,:)

        ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NASA polynomials for enthalpy and enthropy (unitless)
      H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
          p(5)*Tgas4*0.2d0 + p(6)*invT
      S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
          p(5)*Tgas4*0.25d0 + p(7)

      revHS = H - S

      ! else pick NIST data (if present)
    else if (Tlim_nist(idx,2) /= 0.d0) then
      if (Tlim_nist(idx,1) < Tgas .and. Tgas < Tlim_nist(idx,2)) then
        p(:) = p1_nist(idx,:)

      else if (Tlim_nist(idx,2) < Tgas .and. Tgas < Tlim_nist(idx,3)) then
        p(:) = p2_nist(idx,:)

        ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NIST polynomials for enthalpy and enthropy
      ! H in (kJ/mol)
      H = p(1)*Tnist + p(2)*0.5d0*Tnist2 + p(3)*Tnist3/3.d0 + p(4)*Tnist4*0.25d0&
          - p(5)*invTnist + p(6)
      !  Unitsless
      H = H / (Rgas_kJ * Tgas)

      ! S in (J/mol*K)
      S = p(1)*lnTnist + p(2)*Tnist + p(3)*Tnist2*0.5d0 + p(4)*Tnist3/3.d0&
          - p(5)*invTnist2*0.5d0 + p(7)
      !  Unitless. Note: do not use Tnist
      S = S / Rgas_J

      revHS = H - S

      ! return zero is no data exists
    else
      print *, "No thermochemical data of species index", idx
      revHS = 0.d0

    end if

  end function revHS

  !******************************
  subroutine print_best_flux(n,Tgas,nbestin)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea)
    integer::nbest,idx(nrea),i,nbestin
    character*50::name(nrea)

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux

  !******************************
  subroutine print_best_flux_frac(n,Tgas,frac)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),frac
    integer::idx(nrea),i
    character*50::name(nrea)

    if(frac>1d0) then
      print *,"ERROR: fraction in krome_print_best_flux should be <=1!"
      stop
    end if

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nrea
      if(flux(idx(i))<flux(idx(1))*frac) exit
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux_frac

  !******************************
  subroutine print_best_flux_spec(n,Tgas,nbestin,idx_found)
    !print the first nbestin fluxes for the reactions
    ! that contains the species with index idx_found
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),maxflux
    integer::nbest,idx(nrea),i,nbestin,idx_found
    character*50::name(nrea)
    logical::found

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions
    maxflux = 0d0
    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names
    do i=1,nrea
      found = .false.
      if(arr_r1(i) == idx_found) found = .true.
      if(arr_r2(i) == idx_found) found = .true.
      if(arr_r3(i) == idx_found) found = .true.
      if(arr_p1(i) == idx_found) found = .true.
      if(arr_p2(i) == idx_found) found = .true.
      if(arr_p3(i) == idx_found) found = .true.
      maxflux = max(maxflux,flux(i))
      if(.not.found) flux(i) = 0d0
    end do

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,2E17.8)',idx(i)," ",name(idx(i)),flux(idx(i)),&
          flux(idx(i))/maxflux
    end do

  end subroutine print_best_flux_spec

  !*****************************
  function idx_sort(fin)
    !sorting algorithm: requires an array of real values fin
    ! and returns the sorted index list. descending.
    ! bubblesort: not very efficient, replace with what you prefer
    implicit none
    real*8::fin(:),f(size(fin)),ftmp
    integer::idx_sort(size(fin)),n,itmp,i
    logical::found

    f(:) = fin(:) !copy to local

    n = size(f)
    !init indexes
    do i=1,n
      idx_sort(i) = i
    end do

    !loop to sort
    do
      found = .false. !swapped something flag
      do i=2,n
        !> for descending, < for ascending
        if(f(i)>f(i-1)) then
          found = .true.
          !swap real value
          ftmp = f(i)
          f(i) = f(i-1)
          f(i-1) = ftmp
          !swap index
          itmp = idx_sort(i)
          idx_sort(i) = idx_sort(i-1)
          idx_sort(i-1) = itmp
        end if
      end do
      !if nothing swapped exit
      if(.not.found) exit
    end do

  end function idx_sort

  !******************************
  function get_flux(n,Tgas)
    !get the flux k*n*n*... of the rates
    use krome_commons
    implicit none
    integer::i
    integer::r1,r2,r3
    real*8::get_flux(nrea),n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      r3 = arr_r3(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)*n(r3)
    end do
    get_flux(:) = arr_flux(:)

  end function get_flux

  !*****************************
  subroutine load_arrays()
    !load the array containing reactants
    ! and product index
    use krome_commons

    arr_r1(1:281) = (/5,20,20,6,21,21,21,21,6,6,7,7,7,7,36,5,2,5&
        ,5,22,7,7,7,7,2,2,2,2,2,22,22,22,7,5,5,7,23,23,23,24,24,8,9&
        ,24,9,9,8,23,8,8,8,10,25,25,25,8,12,12,12,12,12,12,13,13,13&
        ,13,13,14,14,9,10,10,7,8,8,9,9,10,16,17,17,17,17,11,22,27,8,8&
        ,8,23,28,28,28,29,29,29,31,31,31,14,24,9,9,9,10,10,10,10,32&
        ,33,16,16,16,16,16,16,34,17,17,17,8,8,11,11,11,11,26,26,26,12&
        ,12,13,13,13,13,13,13,13,13,14,10,10,10,10,16,16,16,16,16,16&
        ,16,16,17,17,17,11,11,30,3,4,21,27,27,28,29,29,29,31,31,31,32&
        ,33,33,33,34,34,34,34,35,30,26,26,25,2,2,2,3,3,3,4,4,4,7,7,8&
        ,8,8,8,8,23,23,23,23,9,9,9,10,40,40,40,40,40,5,6,21,9,8,7,2&
        ,12,12,14,10,10,16,16,17,17,7,11,7,22,27,27,3,28,13,13,29,31&
        ,31,4,32,33,33,33,33,34,34,34,34,5,6,9,11,11,14,7,7,7,8,12,17&
        ,17,10,13,16,15,15,7,8,8,8,8,23,23,8,8,5,10&
        ,9/)
    arr_r2(1:281) = (/1,1,1,1,1,1,1,5,20,20,6,21,21,21,1&
        ,1,5,20,20,5,20,20,1,5,1,5,5,20,20,1,1,2,7,5,5,5,1,1,1,1,1,1&
        ,1,5,20,21,20,5,21,21,21,5,7,11,11,7,5,7,8,9,9,9,5,9,9,9,9,9&
        ,9,7,5,5,10,10,10,10,10,10,5,5,7,8,8,5,7,5,22,27,27,7,5,7,9,5&
        ,7,9,5,9,9,24,7,22,27,27,27,27,23,23,7,7,27,27,23,23,23,23,8&
        ,23,23,29,35,35,27,27,27,27,8,16,16,20,20,20,20,20,20,21,21&
        ,21,21,21,20,20,21,21,20,20,21,21,21,21,21,21,20,21,21,21,21&
        ,5,20,20,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,8,9,10&
        ,5,7,9,5,7,8,20,20,1,5,7,8,9,5,7,9,9,1,5,9,5,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,8,8,9&
        ,9,9,9,24,24,9,5,9/)
    arr_r3(1:281) = (/40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,6,5,5,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40/)
    arr_p1(1:281) = (/20,5,5,21,6,6,36,6,21,21,5,6,6,21&
        ,21,2,7,22,22,7,22,22,5,5,5,5,5,5,22,5,5,5,7,7,7,7,8,8,8,9,9&
        ,23,24,9,24,24,23,8,23,23,23,9,26,26,26,12,8,13,14,11,26,10&
        ,12,11,11,15,12,11,11,10,9,9,16,5,5,5,5,16,7,10,10,11,11,8,27&
        ,22,28,28,29,28,23,29,30,28,31,26,29,25,26,30,5,5,7,5,7,7,5,5&
        ,33,34,7,7,25,26,26,33,26,30,11,26,9,17,7,7,7,7,11,11,11,28&
        ,28,7,7,5,5,6,6,6,6,23,32,32,24,24,5,5,6,6,6,6,6,6,35,35,24&
        ,23,8,11,8,9,5,7,5,8,12,8,8,13,12,12,9,9,10,9,10,9,5,10,9,8&
        ,11,10,11,12,10,16,12,13,11,10,16,11,5,27,3,12,13,14,11,28,29&
        ,30,30,4,10,17,16,40,40,40,40,40,20,21,36,24,23,22,5,8,28,8,9&
        ,32,10,33,35,9,20,8,5,5,7,22,8,8,12,29,28,29,28,9,9,22,20,24&
        ,32,20,22,33,32,20,21,24,8,30,8,5,20,22,23,8,9,35,9,29,10,11&
        ,26,5,14,14,11,11,30,30,30,30,10,16&
        ,17/)
    arr_p2(1:281) = (/1,40,40,1,40,40,1,20,5,5,5,22,5,5&
        ,40,40,1,40,40,20,5,5,5,5,1,5,5,5,1,5,5,7,5,6,5,7,40,40,40,40&
        ,40,1,1,20,5,6,5,20,6,6,6,5,7,11,11,5,7,5,5,5,1,8,7,5,7,5,10&
        ,8,8,5,7,7,5,11,11,17,17,9,10,9,10,9,9,10,5,7,5,7,5,5,7,5,5,7&
        ,5,5,7,7,7,8,32,32,32,33,33,33,30,30,5,5,34,34,5,5,5,8,7,9,24&
        ,10,30,23,26,26,25,25,28,34,34,5,5,28,28,29,29,7,7,5,5,8,5,5&
        ,6,6,33,33,10,10,32,32,33,33,5,6,6,6,6,20,5,5,6,5,5,5,5,7,5,5&
        ,7,5,5,7,5,5,5,5,16,7,9,9,5,8,5,1,1,1,1,1,1,1,1,1,5,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,1,1,1,1,1,1&
        ,1,5,1,8,5,1,5,1,1,9,5,9,5,20,20,5,1,20,5,1,5,5,7,1,20,9,10,7&
        ,5,16,10,5,7,1,1,1,9,1,8,5,2,1,1,5,9,1,5,1,5,5,1,20,40,40,40&
        ,40,40,40,40,40,40,40,40/)
    arr_p3(1:281) = (/1,40,40,1,40&
        ,40,1,40,40,40,6,40,20,5,40,40,40,40,40,40,40,40,1,5,1,1,1,40&
        ,40,40,40,40,5,40,40,40,40,40,40,40,40,1,1,40,40,40,40,40,40&
        ,40,40,5,40,40,40,40,40,40,40,40,40,40,40,5,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,23,23,28,28,6,40,40,5,5,40,40,20,20,5,5,40,40&
        ,40,40,9,9,24,40,40,40,40,40,5,40,40,40,5,40,40,5,40,40,40,5&
        ,5,7,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,20,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,1,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40&
        ,40,40,40,40,40,40,40,40,40,40,40,40,40,40,1,40,40,40,40,40&
        ,40,40,40,40,40,40/)

  end subroutine load_arrays

  ! ************************************
  ! solves linear least squares
  subroutine llsq(n, x, y, a, b)

    !****************************************************
    !
    !! LLSQ solves a linear least squares problem matching a line to data.
    !
    !  Discussion:
    !
    !    A formula for a line of the form Y = A * X + B is sought, which
    !    will minimize the root-mean-square error to N data points
    !    ( X(I), Y(I) );
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 March 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    In: N, the number of data values.
    !
    !    In: X(N), Y(N), the coordinates of the data points.
    !
    !    Out: A, B, the slope and Y-intercept of the
    !    least-squares approximant to the data.
    !
    implicit none
    integer,intent(in)::n
    real*8,intent(out)::a, b
    real*8,intent(in)::x(n), y(n)
    real*8::bot, top, xbar, ybar

    ! special case
    if(n == 1) then
      a = 0d0
      b = y(1)
      return
    end if

    ! average X and Y
    xbar = sum(x) / n
    ybar = sum(y) / n

    ! compute beta
    top = dot_product(x(:) - xbar, y(:) - ybar)
    bot = dot_product(x(:) - xbar, x(:) - xbar)

    ! if top is zero a is zero
    if(top==0d0) then
      a = 0d0
    else
      a = top / bot
    end if

    b = ybar - a * xbar

  end subroutine llsq

end module krome_subs

!############### MODULE ##############
module krome_stars

end module krome_stars

!############### MODULE ##############
module krome_dust
contains

  subroutine load_int_JQabs_tab(integrand)
    use krome_commons
    use krome_constants
    !
    real*8,dimension(nPhotoBins) :: integrand ! tGamma in each photobin
    ! table related local variables
    character(len=255) :: fname
    integer :: nrec, lunit, i, j
    real(kind=8) :: energyL, energyR
    real(kind=8), allocatable, dimension(:) :: emid, el, eu, tGamma
    !
    ! open file, get nr of records, and read data
    fname = 'dust_table_absorption.dat'
    nrec = 0
    open(newunit=lunit, file=fname)
    do
      read(lunit,*, end=100)
      nrec = nrec + 1
    end do
    100 rewind(lunit)
    allocate(emid(nrec),el(nrec), eu(nrec), tGamma(nrec))
    do j=1,nrec
      read(lunit,*) emid(j), el(j), eu(j), tGamma(j)
    end do
    close(lunit)
    !
    ! use table to find tGamma integrated across each photo bin
    do i=1,nPhotoBins
      energyL = photoBinEleft(i) * eV_to_erg
      energyR = photoBinEright(i) * eV_to_erg !energy of the bin in erg
      integrand(i) = 0.0_8
      do j = 1, nrec
        if (el(j) <= energyR .and. eu(j) >= energyL) then
          if (el(j) >= energyL .and. eu(j) <= energyR) &
              integrand(i) = integrand(i) + tGamma(j)
          if (el(j) >= energyL .and. eu(j) > energyR) &
              integrand(i) = integrand(i) + tGamma(j) * (energyR - el(j)) / (eu(j) - el(j))
          if (el(j) < energyL .and. eu(j) <= energyR) &
              integrand(i) = integrand(i) + tGamma(j) * (eu(j) - energyL) / (eu(j) - el(j))
          if (el(j) < energyL .and. eu(j) > energyR) &
              integrand(i) = integrand(i) + tGamma(j) * (energyR - energyL) / (eu(j) - el(j))
        end if
      end do
    end do
    !
    deallocate(emid,el,eu,tGamma)
  end subroutine load_int_JQabs_tab
  !***********************
  !compute the absorbed radiation by the dust by integrating over the photobins
  !asuming a tabularised dust distribution
  function get_int_JQabs_tab()
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_int_JQabs_tab
    real*8,dimension(nPhotoBins), save :: integrand ! tGamma in each photobin
    logical, save :: first_call = .true.
    !
    ! Load tabularised data on first call.
    ! Make a double-nested threadprivate check of if this is first call
    ! to avoid race-condition without parallel overhead once table has been loaded
    !
    if (first_call) then
      !$omp critical
      if (first_call) call load_int_JQabs_tab(integrand)
      first_call = .false.
      !$omp end critical
    end if

    !loop over photo bins
    get_int_JQabs_tab = sum(photoBinJ * integrand) * eV_to_erg ! make sure result is in erg s^-1 cm^-3 as for table file
    !
  end function get_int_JQabs_tab
  !
  subroutine setup_2D_dust_tables
    use krome_commons
    implicit none
    integer       :: iAv
    real(kind=8)  :: lambda_abs, wl, wu, log10_Av
    logical,      save :: first_call=.true.
    integer,      save :: nAv, nrec
    real(kind=8), save :: log10_Av_lb, dlog10_Av
    real(kind=8), allocatable, dimension(:), save :: Av_tab, lambda_tab
    real(kind=8), allocatable, dimension(:,:,:), save :: dust_tab_H2_3D, dust_tab_cool_3D, dust_tab_Tdust_3D

    ! Compute absorption from radiation field
    lambda_abs = get_int_JQabs_tab()

    ! Make sure 3D tables are loaded
    !
    ! This will also define
    !    dust_tab_ngas, dust_tab_Tgas,
    !    dust_mult_ngas, dust_mult_Tgas
    ! for 2D table operation, and log10_Av_lb, dlog10_Av
    if (first_call) then
      !$omp critical
      if (first_call) call load_tables()
      first_call = .false.
      !$omp end critical
    end if

    ! Translate from absorption to an Av
    log10_Av = log10(convert_to_Av(lambda_abs))

    ! Find index (iAv) and weights(wl,wu; wl+wu=1)
    iAv = floor((log10_Av - log10_Av_lb)/dlog10_Av)+1
    iAv = min(max(iAv,1),nAv)
    wu  = min(max(log10_Av - ((iAv-1)*dlog10_Av + log10_Av_lb),0.0_8),1.0_8)
    wl = 1.0_8 - wu

    ! Interpolate in 3D table to generate 2D tables. Check if we are on edge point for ub (then wl==1.)
    if (iAv < nAv) then
      dust_tab_H2    = dust_tab_H2_3D(:,:,iAv)*wl + dust_tab_H2_3D(:,:,iAv+1)*wu
      dust_tab_cool  = dust_tab_cool_3D(:,:,iAv)*wl + dust_tab_cool_3D(:,:,iAv+1)*wu
      dust_tab_Tdust = dust_tab_Tdust_3D(:,:,iAv)*wl + dust_tab_Tdust_3D(:,:,iAv+1)*wu
    else
      dust_tab_H2    = dust_tab_H2_3D(:,:,nAv)
      dust_tab_cool  = dust_tab_cool_3D(:,:,nAv)
      dust_tab_Tdust = dust_tab_Tdust_3D(:,:,nAv)
    endif
  contains
    subroutine load_tables
      use krome_fit
      implicit none
      integer :: lunit, i, j, nTgas, nngas
      real(kind=8), dimension(:), allocatable :: dust_tab_Av
      character(len=255) :: fname
      ! Load Av-Lambda table
      fname = 'dust_table_Av_Lambda.dat'
      open(newunit=lunit, file=fname)
      nrec = 0
      do
        read(lunit,*, end=200)
        nrec = nrec + 1
      end do
      200 rewind(lunit)
      allocate(Av_tab(nrec), lambda_tab(nrec))
      do j=1,nrec
        read(lunit,*) Av_tab(j), lambda_tab(j)
      end do
      close(lunit)
      ! Load the three 3D tables
      nTgas = 50 ! FIXME, hardcoded!
      nngas = 50
      nAv = 20
      allocate(dust_tab_Tdust_3D(nTgas,nngas,nAv), &
          dust_tab_cool_3D(nTgas,nngas,nAv), &
          dust_tab_H2_3D(nTgas,nngas,nAv), &
          dust_tab_Av(nAv))
      !
      call init_anytab3D("dust_table_cool_3D.dat",dust_tab_ngas(:), &
          dust_tab_Tgas(:), dust_tab_Av(:), &
          dust_tab_cool_3D(:,:,:), dust_mult_ngas, &
          dust_mult_Tgas, dlog10_Av)
      call init_anytab3D("dust_table_Tdust_3D.dat",dust_tab_ngas(:), &
          dust_tab_Tgas(:), dust_tab_Av(:), &
          dust_tab_Tdust_3D(:,:,:), dust_mult_ngas, &
          dust_mult_Tgas, dlog10_Av)
      call init_anytab3D("dust_table_H2_3D.dat",dust_tab_ngas(:), &
          dust_tab_Tgas(:), dust_tab_Av(:), &
          dust_tab_H2_3D(:,:,:), dust_mult_ngas, &
          dust_mult_Tgas, dlog10_Av)

      log10_Av_lb = dust_tab_Av(1)
      dlog10_Av   = dust_tab_Av(2) - dust_tab_Av(1)

      deallocate(dust_tab_Av)
    end subroutine load_tables
    !
    function convert_to_Av(lambda_abs) result(Av)
      implicit none
      real(kind=8), intent(in) :: lambda_abs
      real(kind=8)             :: Av
      !
      real(kind=8) :: w
      integer :: i
      !
      if (lambda_abs >= lambda_tab(1)) then
        Av = Av_tab(1)
        return
      endif
      !
      if (lambda_abs <= lambda_tab(nrec)) then
        Av = Av_tab(nrec)
        return
      endif
      ! linear interpolation for Av in log of lambda_abs
      do i=2,nrec
        if (lambda_abs > lambda_tab(i)) then
          w = (log(lambda_abs) - log(lambda_tab(i))) / (log(lambda_tab(i-1)) - log(lambda_tab(i)))
          Av = Av_tab(i) * (1.0_8 - w) + Av_tab(i-1) * w
          return
        end if
      end do
      allocate(Av_tab(nrec), lambda_tab(nrec))
    end function convert_to_Av
  end subroutine setup_2D_dust_tables

  !***********************
  subroutine init_dust_tabs()
    use krome_commons
    use krome_fit
    implicit none

  end subroutine init_dust_tabs

end module krome_dust

!############### MODULE ##############
module krome_photo
contains

  !*******************************
  !load a frequency-dependent opacity table stored in fname file,
  ! column 1 is energy or wavelenght in un units of unitEnergy
  ! (default eV), column 2 is opacity in cm2/g.
  ! opacity is interpolated over the current photo-binning.
  subroutine load_opacity_table(fname, unitEnergy)
    use krome_commons
    use krome_constants
    implicit none
    integer,parameter::ntmp=int(1e5)
    character(len=*)::fname
    character(len=*),optional::unitEnergy
    character*10::eunit
    integer::ios,icount,iR,iL,i,j,fileUnit
    real*8::wl,opac,fL,fR,kk,dE
    real*8::wls(ntmp),opacs(ntmp)
    real*8,allocatable::energy(:),kappa(:)

    !read energy unit optional argument
    eunit = "eV" !default is eV
    if(present(unitEnergy)) then
      eunit = trim(unitEnergy)
    end if

    !read form file
    open(newunit=fileUnit,file=trim(fname),status="old",iostat=ios)
    !error if problems reading file
    if(ios/=0) then
      print *,"ERROR: problem while loading "//trim(fname)
      stop
    end if
    icount = 0
    !loop on file lines
    do
      !read wavelength and opacity
      read(fileUnit,*,iostat=ios) wl,opac
      if(ios/=0) exit
      icount = icount + 1
      wls(icount) = wl
      opacs(icount) = opac
    end do
    close(fileUnit)

    !allocate arrays
    allocate(energy(icount), kappa(icount))
    !copy temp arrays into allocated arrays, converting units
    if(trim(eunit)=="eV") then
      !eV->eV (default)
      kappa(:) = opacs(1:icount)
      energy(:) = wls(1:icount)
    elseif(trim(eunit)=="micron") then
      !micron->eV
      kappa(:) = opacs(1:icount)
      energy(:) = planck_eV*clight/(wls(1:icount)*1d-4)
    else
      print *,"ERROR: in load opacity table energy unit unknow",trim(eunit)
      stop
    end if

    !reverse array if necessary
    if(energy(2)<energy(1)) then
      energy(:) = energy(size(energy):1:-1)
      kappa(:) = kappa(size(kappa):1:-1)
    end if

    !check if photobins are intialized
    if(maxval(photoBinEleft)==0d0) then
      print *,"ERROR: empty photobins when interpolating dust Qabs"
      print *," from file "//trim(fname)
      print *,"You probably need to define a photobins metric before"
      print *," the call to krome_load_opacity_table"
      stop
    end if

    !check lower limit
    if(photoBinEleft(1)<energy(1)) then
      print *,"ERROR: dust table "//trim(fname)//" energy lower bound (eV)"
      print *,photoBinEleft(1), "<", energy(1)
      stop
    end if

    !check upper limit
    if(photoBinEright(nPhotoBins)>energy(size(energy))) then
      print *,"ERROR: dust table "//trim(fname)//" energy upper bound (eV)"
      print *,photoBinEright(nPhotoBins), ">", energy(size(energy))
      stop
    end if

    !interpolate on current energy distribution
    do j=1,nPhotoBins
      do i=2,size(energy)
        !find left bound position
        if(photoBinEleft(j)>energy(i-1) &
            .and. photoBinEleft(j)<energy(i)) then
        dE = energy(i)-energy(i-1)
        fL = (photoBinEleft(j)-energy(i-1))/dE &
            * (kappa(i)-kappa(i-1)) + kappa(i-1)
        iL = i
      end if

      !find right bound position
      if(photoBinEright(j)>energy(i-1) &
          .and. photoBinEright(j)<energy(i)) then
      dE = energy(i)-energy(i-1)
      fR = (photoBinEright(j)-energy(i-1))/dE &
          * (kappa(i)-kappa(i-1)) + kappa(i-1)
      iR = i
    end if
  end do

  !sum opacity for the given photo bin
  kk = 0d0
  !if there are other opacity points in between left and right limits
  if(iR-iL>0) then
    kk = kk + (energy(iL)-photoBinEleft(j))*(fL+kappa(iL))/2d0
    kk = kk + (photoBinEright(j)-energy(iR-1))*(fR+kappa(iR-1))/2d0
    !sum points in between
    do i=iL,iR-2
      kk = kk + (energy(i+1)-energy(i))*(kappa(i+1)+kappa(i))/2d0
    end do
  elseif(iR==iL) then
    !no opacity points in between
    kk = kk + (fL+fR)*(photoBinEright(j)-photoBinEleft(j))/2d0
  else
    print *,"ERROR: dust opacity interpolation error, iR-iL<0!"
    print *,"iR,iL:",iR,iL
    stop
  end if

  !copy to common and scale to bin size
  dE = photoBinEright(j)-photoBinEleft(j)
  opacityDust(j) = kk/dE

end do

!dump interpolated opacity
open(newunit=fileUnit,file="opacityDust.interp",status="replace")
do j=1,nPhotoBins
  write(fileUnit,*) photoBinEmid(j),opacityDust(j)
end do
close(fileUnit)

!dump original opacity file (as loaded by krome)
open(newunit=fileUnit,file="opacityDust.org",status="replace")
do i=1,size(energy)
  write(fileUnit,*) energy(i),kappa(i)
end do
close(fileUnit)

end subroutine load_opacity_table

!*************************
!get the intensity of the photon flux at
! a given energy in eV.
! returned value is in eV/cm2/s/Hz
function get_photoIntensity(energy)
use krome_commons
implicit none
real*8::get_photoIntensity,energy
integer::i

!check if requested energy is lower than the lowest limit
if(energy<photoBinEleft(1)) then
  get_photoIntensity = 0d0 !photoBinJ(1)
  return
end if

!check if requested energy is greater that the the largest limit
if(energy>photoBinEright(nPhotoBins)) then
  get_photoIntensity = 0d0 !photoBinJ(nPhotoBins)
  return
end if

!look for the interval
do i=1,nPhotoBins
  if(photoBinEleft(i).le.energy .and. photoBinEright(i).ge.energy) then
    get_photoIntensity = photoBinJ(i)
    return
  end if
end do

!error if nothing found
print *,"ERROR: no interval found in get_photoIntensity"
print *,"energy:",energy,"eV"
stop !halt program

end function get_photoIntensity

!*********************
!initialize/tabulate the bin-based xsecs
subroutine init_photoBins(Tgas)
use krome_constants
use krome_commons
use krome_dust
use krome_getphys
implicit none
integer::i,j
real*8::Tgas,imass(nspec),kt2
real*8::energy_eV,kk,energyL,energyR,dshift(nmols)

!rise error if photobins are not defined
if(photoBinEmid(nPhotoBins)==0d0) then
  print *,"ERROR: when using photo bins you must define"
  print *," the energy interval in bins!"
  stop
end if

!get inverse of mass
imass(:) = get_imass()

!precompute adimensional line broadening
dshift(:) = 0d0

call load_xsec("swri_C__C+_E.dat", xsec217_val, xsec217_Emin, xsec217_n, xsec217_idE)
call load_xsec("swri_H2__H2+_E.dat", xsec218_val, xsec218_Emin, xsec218_n, xsec218_idE)
call load_xsec("swri_H-__H_E.dat", xsec219_val, xsec219_Emin, xsec219_n, xsec219_idE)
call load_xsec("swri_CH__C_H.dat", xsec220_val, xsec220_Emin, xsec220_n, xsec220_idE)
call load_xsec("swri_CH__CH+_E.dat", xsec221_val, xsec221_Emin, xsec221_n, xsec221_idE)
call load_xsec("swri_C2__C_C.dat", xsec222_val, xsec222_Emin, xsec222_n, xsec222_idE)
call load_xsec("swri_OH__O_H.dat", xsec223_val, xsec223_Emin, xsec223_n, xsec223_idE)
call load_xsec("swri_OH__OH+_E.dat", xsec224_val, xsec224_Emin, xsec224_n, xsec224_idE)
call load_xsec("swri_H2O__OH_H.dat", xsec225_val, xsec225_Emin, xsec225_n, xsec225_idE)
call load_xsec("swri_H2O__H2O+_E.dat", xsec226_val, xsec226_Emin, xsec226_n, xsec226_idE)
call load_xsec("swri_O2__O2+_E.dat", xsec227_val, xsec227_Emin, xsec227_n, xsec227_idE)
call load_xsec("swri_O2__O_O.dat", xsec228_val, xsec228_Emin, xsec228_n, xsec228_idE)
call load_xsec("swri_H2__H+_H_E.dat", xsec229_val, xsec229_Emin, xsec229_n, xsec229_idE)

!tabulate the xsecs into a bin-based array
do j=1,nPhotoBins
  energyL = photoBinEleft(j)
  energyR = photoBinEright(j)
  energy_eV = photoBinEmid(j) !energy of the bin in eV

  !H -> H+ + E
  kk = 0d0
  if(energy_eV>1.360d+01.and.energy_eV<5.000d+04) kk =  sigma_v96(energy_ev, 4.298d-01, 5.475d+04, 3.288d+01, 2.963d+00, 0.000d+00, 0.000d+00, 0.000d+00)
  !$omp parallel
  photoBinJTab(1,j) = kk
  !$omp end parallel

  !HE -> HE+ + E
  kk = 0d0
  if(energy_eV>2.459d+01.and.energy_eV<5.000d+04) kk =  sigma_v96(energy_ev, 1.361d+01, 9.492d+02, 1.469d+00, 3.188d+00, 2.039d+00, 4.434d-01, 2.136d+00)
  !$omp parallel
  photoBinJTab(2,j) = kk
  !$omp end parallel

  !HE+ -> HE++ + E
  kk = 0d0
  if(energy_eV>5.442d+01.and.energy_eV<5.000d+04) kk =  sigma_v96(energy_ev, 1.720d+00, 1.369d+04, 3.288d+01, 2.963d+00, 0.000d+00, 0.000d+00, 0.000d+00)
  !$omp parallel
  photoBinJTab(3,j) = kk
  !$omp end parallel

  !O -> O+ + E
  kk = 0d0
  if(energy_eV>1.362d+01.and.energy_eV<5.380d+02) kk =  sigma_v96(energy_ev, 1.240d+00, 1.745d+03, 3.784d+00, 1.764d+01, 7.589d-02, 8.698d+00, 1.271d-01)
  !$omp parallel
  photoBinJTab(4,j) = kk
  !$omp end parallel

  !C -> C+ + E
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec217_val(:), xsec217_Emin,xsec217_idE, dshift(idx_C))
  !$omp parallel
  photoBinJTab(5,j) = kk
  !$omp end parallel

  !H2 -> H2+ + E
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec218_val(:), xsec218_Emin,xsec218_idE, dshift(idx_H2))
  !$omp parallel
  photoBinJTab(6,j) = kk
  !$omp end parallel

  !H- -> H + E
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec219_val(:), xsec219_Emin,xsec219_idE, dshift(idx_Hk))
  !$omp parallel
  photoBinJTab(7,j) = kk
  !$omp end parallel

  !CH -> C + H
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec220_val(:), xsec220_Emin,xsec220_idE, dshift(idx_CH))
  !$omp parallel
  photoBinJTab(8,j) = kk
  !$omp end parallel

  !CH -> CH+ + E
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec221_val(:), xsec221_Emin,xsec221_idE, dshift(idx_CH))
  !$omp parallel
  photoBinJTab(9,j) = kk
  !$omp end parallel

  !C2 -> C + C
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec222_val(:), xsec222_Emin,xsec222_idE, dshift(idx_C2))
  !$omp parallel
  photoBinJTab(10,j) = kk
  !$omp end parallel

  !OH -> O + H
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec223_val(:), xsec223_Emin,xsec223_idE, dshift(idx_OH))
  !$omp parallel
  photoBinJTab(11,j) = kk
  !$omp end parallel

  !OH -> OH+ + E
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec224_val(:), xsec224_Emin,xsec224_idE, dshift(idx_OH))
  !$omp parallel
  photoBinJTab(12,j) = kk
  !$omp end parallel

  !H2O -> OH + H
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec225_val(:), xsec225_Emin,xsec225_idE, dshift(idx_H2O))
  !$omp parallel
  photoBinJTab(13,j) = kk
  !$omp end parallel

  !H2O -> H2O+ + E
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec226_val(:), xsec226_Emin,xsec226_idE, dshift(idx_H2O))
  !$omp parallel
  photoBinJTab(14,j) = kk
  !$omp end parallel

  !O2 -> O2+ + E
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec227_val(:), xsec227_Emin,xsec227_idE, dshift(idx_O2))
  !$omp parallel
  photoBinJTab(15,j) = kk
  !$omp end parallel

  !O2 -> O + O
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec228_val(:), xsec228_Emin,xsec228_idE, dshift(idx_O2))
  !$omp parallel
  photoBinJTab(16,j) = kk
  !$omp end parallel

  !H2 -> H+ + H + E
  kk = 0d0
  if(energy_eV>0.0.and.energy_eV<1.d8) kk = xsec_interp(energyL, energyR, xsec229_val(:), xsec229_Emin,xsec229_idE, dshift(idx_H2))
  !$omp parallel
  photoBinJTab(17,j) = kk
  !$omp end parallel

end do

!save interpolated xsecs to file
call save_xsec("swri_C__C+_E.interp",5)
call save_xsec("swri_H2__H2+_E.interp",6)
call save_xsec("swri_H-__H_E.interp",7)
call save_xsec("swri_CH__C_H.interp",8)
call save_xsec("swri_CH__CH+_E.interp",9)
call save_xsec("swri_C2__C_C.interp",10)
call save_xsec("swri_OH__O_H.interp",11)
call save_xsec("swri_OH__OH+_E.interp",12)
call save_xsec("swri_H2O__OH_H.interp",13)
call save_xsec("swri_H2O__H2O+_E.interp",14)
call save_xsec("swri_O2__O2+_E.interp",15)
call save_xsec("swri_O2__O_O.interp",16)
call save_xsec("swri_H2__H+_H_E.interp",17)

!energy tresholds (eV)
!$omp parallel
photoBinEth(1) = 1.360d+01 !H -> H+ + E
!$omp end parallel
!$omp parallel
photoBinEth(2) = 2.459d+01 !HE -> HE+ + E
!$omp end parallel
!$omp parallel
photoBinEth(3) = 5.442d+01 !HE+ -> HE++ + E
!$omp end parallel
!$omp parallel
photoBinEth(4) = 1.362d+01 !O -> O+ + E
!$omp end parallel
!$omp parallel
photoBinEth(5) = 0.0 !C -> C+ + E
!$omp end parallel
!$omp parallel
photoBinEth(6) = 0.0 !H2 -> H2+ + E
!$omp end parallel
!$omp parallel
photoBinEth(7) = 0.0 !H- -> H + E
!$omp end parallel
!$omp parallel
photoBinEth(8) = 0.0 !CH -> C + H
!$omp end parallel
!$omp parallel
photoBinEth(9) = 0.0 !CH -> CH+ + E
!$omp end parallel
!$omp parallel
photoBinEth(10) = 0.0 !C2 -> C + C
!$omp end parallel
!$omp parallel
photoBinEth(11) = 0.0 !OH -> O + H
!$omp end parallel
!$omp parallel
photoBinEth(12) = 0.0 !OH -> OH+ + E
!$omp end parallel
!$omp parallel
photoBinEth(13) = 0.0 !H2O -> OH + H
!$omp end parallel
!$omp parallel
photoBinEth(14) = 0.0 !H2O -> H2O+ + E
!$omp end parallel
!$omp parallel
photoBinEth(15) = 0.0 !O2 -> O2+ + E
!$omp end parallel
!$omp parallel
photoBinEth(16) = 0.0 !O2 -> O + O
!$omp end parallel
!$omp parallel
photoBinEth(17) = 0.0 !H2 -> H+ + H + E
!$omp end parallel

!interpolate dust qabs

!map with X->B/C transition to bin corrspondence

end subroutine init_photoBins

!**********************
!save xsecs with index idx to file
subroutine save_xsec(fname,idx)
use krome_commons
implicit none
character(len=*)::fname
integer::idx,j
real*8::energyLeft,energyRight

open(22,file=trim(fname),status="replace")
do j=1,nPhotoBins
  energyLeft = photoBinELeft(j) !left bin energy, eV
  energyRight = photoBinERight(j) !right bin energy, eV
  write(22,*) energyLeft, energyRight, photoBinJTab(idx,j)
end do
close(22)

end subroutine save_xsec

!**********************
!compute integrals to derive phtorates (thin)
subroutine calc_photoBins()
use krome_commons
implicit none
real*8::n(nspec)

n(:) = 0d0
call calc_photoBins_thick(n)

end subroutine calc_photoBins

!**********************
!compute integrals to derive phtorates (thick)
subroutine calc_photoBins_thick(n)
use krome_commons
use krome_constants
use krome_subs
use krome_getphys
implicit none
integer::i,j
real*8::dE,kk,Jval,E,Eth,n(:),ncol(nmols),tau

!init rates and heating
photoBinRates(:) = 0d0 !1/s/Hz
photoBinHeats(:) = 0d0 !eV/s/Hz
GHabing_thin = 0d0 !habing flux
!loop on energy bins
do j=1,nPhotoBins
  dE = photoBinEdelta(j) !energy interval, eV
  E = photoBinEmid(j) !energy of the bin in eV
  Jval = photoBinJ(j) !radiation intensity eV/s/cm2/sr/Hz
  if(E>=6d0.and.E<=13.6)then
    GHabing_thin = GHabing_thin + Jval * dE
  endif
  tau = 0d0
  !loop on reactions
  do i=1,nPhotoRea
    Eth = photoBinEth(i) !reaction energy treshold, eV
    if(E>Eth) then
      !approx bin integral
      kk = photoBinJTab(i,j)*Jval/E*dE
      photoBinRates(i) = photoBinRates(i) + kk
      photoBinHeats(i) = photoBinHeats(i) + kk*(E-Eth)
    end if
  end do
end do

!Final Habing flux
GHabing_thin = GHabing_thin * 4d0 * pi / (1.6d-3) * iplanck_eV * eV_to_erg

!converts to 1/s
photoBinRates(:) = 4d0*pi*photoBinRates(:) * iplanck_eV

!converts to erg/s
photoBinHeats(:) = 4d0*pi*photoBinHeats(:) * iplanck_eV * eV_to_erg

end subroutine calc_photoBins_thick

!********************
!Verner+96 cross section fit (cm2)
function sigma_v96(energy_eV,E0,sigma_0,ya,P,yw,y0,y1)
implicit none
real*8::sigma_v96,energy_eV,sigma_0,Fy,yw,x,y,E0
real*8::y0,y1,ya,P
x = energy_eV/E0 - y0
y = sqrt(x**2 + y1**2)
Fy = ((x - 1.d0)**2 + yw**2) *  y**(0.5*P-5.5) &
    * (1.d0+sqrt(y/ya))**(-P)
sigma_v96 = 1d-18 * sigma_0 * Fy !cm2
end function sigma_v96

!********************
!Verner+96 cross section fit (cm2)
!Average by numerical integration
function sigma_v96_int(E_low,E_high,E0,sigma_0,ya,P,yw,y0,y1)
real*8::sigma_v96_int,E_low,E_high,sigma_0,integral,yw,x,y,E0
real*8::y0,y1,ya,P
real*8::binWidth,dE,E
integer::i
integer,parameter::N=100
integral = 0d0
binWidth = E_high-E_low
dE = binWidth/real(N,kind=8)
do i=1,N
  E = E_low + (i-0.5)*dE
  integral = integral + sigma_v96(E,E0,sigma_0,ya,P,yw,y0,y1)*dE
end do
sigma_v96_int = integral / binWidth !cm2
end function sigma_v96_int

!********************
function heat_v96(energy_eV,Eth,E0,sigma_0,ya,P,yw,y0,y1)
!Heating with Verner+96 cross section fit (cm2*eV)
use krome_constants
real*8::heat_v96,energy_eV,sigma_0,Fy,yw,x,y,E0,Eth
real*8::y0,y1,ya,P
x = energy_eV/E0 - y0
y = sqrt(x**2 + y1**2)
Fy = ((x - 1.d0)**2 + yw**2) *  y**(0.5*P-5.5) &
    * (1.d0+sqrt(y/ya))**(-P)
heat_v96 = 1d-18 * sigma_0 * Fy * (energy_eV - Eth) !cm2*eV
end function heat_v96

!************************
!load the xsecs from file and get limits
subroutine load_xsec(fname,xsec_val,xsec_Emin,xsec_n,xsec_idE)
implicit none
real*8,allocatable::xsec_val(:)
real*8::xsec_Emin,xsec_dE,xsec_val_tmp(int(1e6)),rout(2)
real*8::xsec_E_tmp(size(xsec_val_tmp)),xsec_idE,diff
integer::xsec_n,ios
character(*)::fname

!if file already loaded skip subroutine
if(allocated(xsec_val)) return

xsec_n = 0 !number of lines found
!open file
open(33,file=fname,status="old",iostat=ios)
!check if file exists
if(ios.ne.0) then
  print *,"ERROR: problems loading "//fname
  stop
end if

!read file line-by-line
do
  read(33,*,iostat=ios) rout(:) !read line
  if(ios<0) exit !eof
  if(ios/=0) cycle !skip blanks
  xsec_n = xsec_n + 1 !increase line number
  xsec_val_tmp(xsec_n) = rout(2) !read xsec value cm2
  xsec_E_tmp(xsec_n) = rout(1) !read energy value eV
  !compute the dE for the first interval
  if(xsec_n==2) xsec_dE = xsec_E_tmp(2)-xsec_E_tmp(1)
  !check if all the intervals have the same spacing
  if(xsec_n>2) then
    diff = xsec_E_tmp(xsec_n)-xsec_E_tmp(xsec_n-1)
    if(abs(diff/xsec_dE-1d0)>1d-6) then
      print *,"ERROR: spacing problem in file "//fname
      print *," energy points should be equally spaced!"
      print *,"Point number: ",xsec_n
      print *,"Found ",diff
      print *,"Should be",xsec_dE
      stop
    end if
  end if
end do
close(33)

!store the minimum energy
xsec_Emin = xsec_E_tmp(1)
!allocate the array with the values
allocate(xsec_val(xsec_n))
!copy the values from the temp array to the allocated one
xsec_val(:) = xsec_val_tmp(1:xsec_n)
!store the inverse of the delta energy
xsec_idE = 1d0 / xsec_dE

end subroutine load_xsec

!**********************
!return averaged xsec in the energy range [xL,xR]
! units: eV, cm2; broadening shift is adimensional
function xsec_interp(xL,xR,xsec_val,xsec_Emin,xsec_idE,dshift) result(xsecA)
use krome_user_commons
implicit none
real*8::xsecA,dE,dshift,dE_shift,eL,eR,dxi
real*8::energy,xsec_val(:),xsec_Emin,xsec_idE,xL,xR
integer::idx

!xsec energy step (regular grid)
dE = 1d0/xsec_idE
!store inverse of bin size
dxi = 1d0/(xR-xL)
xsecA = 0d0 !init integrated xsec
!loop on xsec vals
do idx=1,size(xsec_val)
  eL = (idx-1)*dE+xsec_Emin !left interval
  eR = eL + dE !right interval
  energy = (eL+eR)/2d0 !mid point

  !compute line broadening
  eL = eL - 0.5d0*dshift*energy
  eR = eR + 0.5d0*dshift*energy

  !if xsec energy in the interval compute area
  if(xR<eL.and.xL<eL) then
    xsecA = xsecA + 0d0
  elseif(xR>eL.and.xL>eL) then
    xsecA = xsecA + 0d0
  else
    !renormalize xsec area considering partial overlap
    xsecA = xsecA +xsec_val(idx) * (min(eR,xR)-max(eL,xL)) &
        * dxi
  end if
end do

end function xsec_interp

!**********************
!linear interpolation for the photo xsec
function xsec_interp_mid(energy,xsec_val,xsec_Emin,xsec_n,xsec_idE)
implicit none
real*8::xsec_interp_mid,E0
real*8::energy,xsec_val(:),xsec_Emin,xsec_idE
integer::xsec_n,idx

xsec_interp_mid = 0d0
!retrive index
idx = (energy-xsec_Emin) * xsec_idE + 1

!lower bound
E0 = xsec_Emin + (idx-1)/xsec_idE

!out of the limits is zero
if(idx<1.or.idx>xsec_n-1) return

!linear interpolation
xsec_interp_mid = (energy-E0) * xsec_idE &
    * (xsec_val(idx+1)-xsec_val(idx)) + xsec_val(idx)

!avoid negative xsec values when outside the limits
xsec_interp_mid = max(xsec_interp_mid,0d0)

end function xsec_interp_mid

!************************
!load photodissociation data from default file
subroutine kpd_H2_loadData()
use krome_commons
implicit none
integer::unit,ios,ii,jj
real*8::xE,dE,pre
character(len=20)::fname

!open file to read
fname = "H2pdB.dat"
open(newunit=unit,file=trim(fname),status="old",iostat=ios)
!check for errors
if(ios/=0) then
  print *,"ERROR: problem loading file "//trim(fname)
  stop
end if

!init data default
H2pdData_EX(:) = 0d0
H2pdData_dE(:,:) = 0d0
H2pdData_pre(:,:) = 0d0

!loop on file to read
do
  read(unit,*,iostat=ios) ii,jj,xE,dE,pre
  !skip comments
  if(ios==59.or.ios==5010) cycle
  !exit when eof
  if(ios/=0) exit
  !store data
  H2pdData_EX(ii+1) = xE !ground level energy, eV
  H2pdData_dE(ii+1,jj+1) = dE !Ej-Ei energy, eV
  H2pdData_pre(ii+1,jj+1) = pre !precomp (see file header)
end do

!check if enough data have been loaded (file size is expected)
if((ii+1/=H2pdData_nvibX).or.(jj+1/=H2pdData_nvibB)) then
  !print error message
  print *,"ERROR: missing data when loading "//fname
  print *,"found:",ii+1,jj+1
  print *,"expected:",H2pdData_nvibX,H2pdData_nvibB
  stop
end if

close(unit)

end subroutine kpd_H2_loadData

!************************
subroutine kpd_bin_map()
use krome_commons
implicit none
integer::i,j,k
logical::found

!loop on excited states (B)
do i=1,H2pdData_nvibB
  !loop on ground states (X)
  do j=1,H2pdData_nvibX
    !if prefactor is zero no need to check map
    ! default is set to 1 (be aware of it!)
    if(H2pdData_pre(j,i)==0d0) then
      H2pdData_binMap(j,i) = 1
      cycle
    end if

    found = .false.
    !loop on bins
    do k=1,nPhotoBins
      !find energy bin corresponding on the given dE
      if((photoBinEleft(k).le.H2pdData_dE(j,i)) &
          .and. (photoBinEright(k).ge.H2pdData_dE(j,i))) then
      H2pdData_binMap(j,i) = k
      found = .true.
    end if
  end do
  !error if outside bounds
  if(.not.found) then
    print *,"ERROR: problem when creating H2"
    print *," photodissociation map!"
    print *," min/max (eV):", minval(photoBinEleft), &
        maxval(photoBinEright)
    print *," transition:",j,i
    print *," corresponding energy (eV):",H2pdData_dE(j,i)
    print *," transitions min/max (eV):", &
        minval(H2pdData_dE, mask=((H2pdData_dE>0d0) .and. &
        (H2pdData_pre>0d0))), &
        maxval(H2pdData_dE, mask=(H2pdData_pre>0d0))
    stop
  end if
end do
end do

end subroutine kpd_bin_map

!************************
!compute vibrational partition function at given Tgas
! for all the loaded energies (for H2 Solomon)
function partitionH2_vib(Tgas) result(z)
use krome_constants
use krome_commons
implicit none
real*8::Tgas,z(H2pdData_nvibX),b
integer::j

!prepare partition function from ground (X) levels energies
b = iboltzmann_eV/Tgas
z(:) = exp(-H2pdData_EX(:)*b)

!normalize
z(:) = z(:)/sum(z)

end function partitionH2_vib

!************************
!compute H2 photodissociation rate (Solomon)
! state to state, using preloded data, 1/s
function kpd_H2(Tgas) result(kpd)
use krome_commons
implicit none
integer::i,j
real*8::Tgas,kpd,dE,z(H2pdData_nvibX)

!get partition for ground state X
z(:) = partitionH2_vib(Tgas)

!compute the rate, using preloaded data
kpd = 0d0
!loop on excited states (B)
do i=1,H2pdData_nvibB
!compute rate for ith state
kpd = kpd + sum(H2pdData_pre(:,i) &
    * photoBinJ(H2pdData_binMap(:,i)) * z(:))
end do

end function kpd_H2

!************************
!photodissociation H2 xsec from atomic data (for opacity)
function kpd_H2_xsec(Tgas) result(xsec)
use krome_constants
use krome_commons
implicit none
real*8::xsec(nPhotoBins),z(H2pdData_nvibX)
real*8::Tgas
integer::i

!get partition for ground state X
z(:) = partitionH2_vib(Tgas)

xsec(:) = 0d0
!loop on excited states (B)
do i=1,H2pdData_nvibB
xsec(H2pdData_binMap(:,i)) = &
    xsec(H2pdData_binMap(:,i)) &
    + H2pdData_pre(:,i)*z(:)
end do

!cm2
xsec(:) = xsec(:)*planck_eV

end function kpd_H2_xsec

!************************
!H2 direct photodissociation in the Lyman-Werner bands
! cross-section in cm^2 fit by Abel et al. 1997 of
! data by Allison&Dalgarno 1969
function H2_sigmaLW(energy_eV)
use krome_commons
implicit none
real*8::H2_sigmaLW,energy_eV
real*8::sL0,sW0,sL1,sW1,fact

!initialization
sL0 = 0d0
sL1 = 0d0
sW0 = 0d0
sW1 = 0d0

if(energy_eV>14.675.and.energy_eV<16.820)then
sL0 = 1d-18*1d1**(15.1289-1.05139*energy_eV)
elseif(energy_eV>16.820.and.energy_eV<17.6d0)then
sL0 = 1d-18*1d1**(-31.41d0+1.8042d-2*energy_eV**3-4.2339d-5*energy_eV**5)
endif

if(energy_eV>14.675d0.and.energy_eV<17.7d0)then
sW0 = 1d-18*1d1**(13.5311d0-0.9182618*energy_eV)
endif

if(energy_eV>14.159d0.and.energy_eV<15.302d0)then
sL1 = 1d-18*1d1**(12.0218406d0-0.819429*energy_eV)
elseif(energy_eV>15.302d0.and.energy_eV<17.2d0)then
sL1 = 1d-18*1d1**(16.04644d0-1.082438*energy_eV)
endif

if(energy_eV>14.159d0.and.energy_eV<17.2d0)then
sW1 = 1d-18*1d1**(12.87367-0.85088597*energy_eV)
endif

fact = 1d0/(phys_orthoParaRatio+1d0)

H2_sigmaLW = fact*(sL0+sW0)+(1d0-fact)*(sL1+sW1)

end function H2_sigmaLW

! *****************************
! load kabs from file
subroutine find_Av_load_kabs2(file_name)
use krome_commons
use krome_constants
implicit none
integer,parameter::imax=10000
character(len=*),intent(in),optional::file_name
character(len=200)::fname
integer::ios, unit, icount, i, j
real*8::tmp_energy(imax), tmp_data(imax), f1, f2, kavg, ksum
real*8,allocatable::Jdraine(:)

! check if energy bins are set
if(maxval(photoBinEleft)==0d0) then
print *, "ERROR: to load kabs for Av G0 finder you"
print *, " have to initialize some energy bins!"
stop
end if

! check if optional argument is present
fname = "kabs_draine_Rv31.dat"
if(present(file_name)) then
fname = trim(file_name)
end if

! open file to read
open(newunit=unit, file=fname, status="old", iostat=ios)
! check if file is there
if(ios/=0) then
print *, "ERROR: Kabs file not found!"
print *, trim(fname)
stop
end if

! loop on file lines
icount = 1
do
read(unit, *, iostat=ios) tmp_energy(icount), &
    tmp_data(icount)
if(ios/=0) exit
icount = icount + 1
end do
close(unit)

! convert microns to eV
tmp_energy(1:icount-1) = planck_eV * clight &
    / (tmp_energy(1:icount-1) * 1d-4)

! get corresponding draine flux
allocate(Jdraine(icount-1))
Jdraine(:) = get_draine(tmp_energy(1:icount-1))

! loop on photobins to get average kabs
do j=1,nPhotoBins
kavg = 0d0
ksum = 0d0
do i=1,icount-2
  ! integrate only in the bin range
  if(tmp_energy(i)>=photoBinEleft(j) &
      .and. tmp_energy(i+1)<=photoBinEright(j)) then
  ! numerator integral Jdraine(E)kabs(E)/E
  f1 = tmp_data(i)*Jdraine(i)/tmp_energy(i)
  f2 = tmp_data(i+1)*Jdraine(i+1)/tmp_energy(i+1)
  kavg = kavg + (f1+f2) / 2d0 &
      * (tmp_energy(i+1)-tmp_energy(i))

  ! denominator integral Jdraine(E)/E
  f1 = Jdraine(i)/tmp_energy(i)
  f2 = Jdraine(i+1)/tmp_energy(i+1)
  ksum = ksum + (f1+f2) / 2d0 &
      * (tmp_energy(i+1)-tmp_energy(i))
end if
!!$           if(tmp_energy(i)<photoBinEmid(j) &
    !!$                .and. tmp_energy(i+1)>photoBinEmid(j)) then
!!$              kavg = (photoBinEmid(j) - tmp_energy(i)) &
    !!$                   / (tmp_energy(i+1) - tmp_energy(i)) &
    !!$                   * (tmp_data(i+1) - tmp_data(i)) + tmp_data(i)
!!$              print *,photoBinEmid(j), kavg
!!$           end if
end do
! ratio of the integral is average absorption in the bin
find_Av_draine_kabs(j) = kavg / (ksum+1d-40)
end do

end subroutine find_Av_load_kabs2

! *****************************
! load kabs from file
subroutine find_Av_load_kabs(file_name)
use krome_commons
use krome_constants
implicit none
character(len=*),intent(in),optional::file_name
character(len=200)::fname
real*8::opacityDust_org(nPhotoBins)

! check if energy bins are set
if(maxval(photoBinEleft)==0d0) then
print *, "ERROR: to load kabs for Av G0 finder you"
print *, " have to initialize some energy bins!"
stop
end if

! check if optional argument is present
fname = "kabs_draine_Rv31.dat"
if(present(file_name)) then
fname = trim(file_name)
end if

opacityDust_org = opacityDust

call load_opacity_table(fname, "micron")

! ratio of the integral is average absorption in the bin
find_Av_draine_kabs = opacityDust

opacityDust = opacityDust_org

end subroutine find_Av_load_kabs

! *********************************
! given the current photo bin intensity distribution
! estimates G0 and Av using the bins in the Draine range
subroutine estimate_G0_Av(G0, Av, n, d2g)
use krome_constants
use krome_commons
use krome_getphys
use krome_subs
implicit none
real*8,intent(out)::G0, Av
real*8,intent(in)::d2g, n(nspec)
real*8::lnG0, mu, ntot
real*8::ydata(nPhotoBins)
integer::i
logical, save ::first_call=.true.
real*8,  save ::XH,Jdraine(nPhotoBins),xdata(nPhotoBins)
integer, save ::lb,ub,ndraine
!$omp threadprivate(first_call,XH,Jdraine,xdata,lb,ub,ndraine)

if (first_call) then
! get non-attenuated draine flux
Jdraine(:) = get_draine(photoBinEmid(:))

! only consider bins that have non-attenuated draine radiation
lb=0
do i=1,nPhotoBins
if (Jdraine(i)>0 .and. lb==0) lb=i
if (Jdraine(i)>0) ub=i
end do
ndraine = ub - lb + 1

! mean molecular weight and gas density
mu = get_mu(n(:))
ntot = sum(n(1:nspec))

! find mass fraction of H-nuclei as xH = mp * n_H / rho
xH = (p_mass * get_Hnuclei(n(:))) / (p_mass * mu * ntot)

! now we can calculate the xdata, which are constant:

! loop on photo bins
do i=lb,ub
! compute x in y = Av*x + ln(G0)
xdata(i-lb+1) = -find_Av_draine_kabs(i) * 1.8d21 * p_mass * d2g / xH
end do

! make sure we only do this once
! NOTICE: we assume hydrogen mass fraction is constant
! NOTICE: but this is implicitly the case anyway
! NOTICE: because below we translate between
! NOTICE: column density and Av
first_call = .false.
end if

! loop on photo bins
do i=lb,ub
! compute y in y = Av*x + ln(G0)
ydata(i - lb + 1) = log(photoBinJ(i) + 1d-200) - log(Jdraine(i))
end do

! needs at least one bin
if(ndraine<=1) then
print *,"ERROR: you want to estimate G0 and Av with less than 2 bins in the"
print *," Draine range, 5-13.6 eV! Nbins(Draine)=",ndraine
stop
end if

! call least squares to compute Av and ln(G0)
call llsq(ndraine, xdata(1:ndraine), ydata(1:ndraine), &
    Av, lnG0)

! Apply prior
if(lnG0 < -7d0 .or. lnG0 > 7d0) then
if(lnG0 < -7d0) lnG0 = -7d0
if(lnG0 > 7d0) lnG0 = 7d0
Av = sum(( ydata(1:ndraine)-lnG0)*xdata(1:ndraine))/sum(xdata(1:ndraine)**2)
end if

if(Av < 0d0) then
Av = 0d0
lnG0 = sum( ydata(1:ndraine))/ndraine
endif

! return G0
G0 = exp(lnG0)

end subroutine estimate_G0_Av

! ************************
function get_draine(energy_list) result(Jdraine)
use krome_commons
use krome_constants
implicit none
integer::i
real*8,intent(in)::energy_list(:)
real*8::x, Jdraine(size(energy_list))

do i=1,size(energy_list)
x = energy_list(i) !eV
!eV/cm2/sr
if(x<13.6d0.and.x>5d0) then
Jdraine(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
    * x *planck_eV
else
Jdraine(i) = 0d0
end if
end do
end function get_draine

end module krome_photo

!############### MODULE ##############
module krome_tabs
contains

!***********************+
function coe_tab(n)
!interface to tabs
use krome_subs
use krome_getphys
use krome_phfuncs
use krome_grfuncs
use krome_constants
use krome_commons
use krome_user_commons
implicit none
integer::idx,j
real*8::Tgas, coe_tab(nrea),n(nspec),small

Tgas = max(n(idx_Tgas),phys_Tcmb)
small = 0d0

coe_tab(:) = coe(n(:))

end function coe_tab

!**************************
!compute rates that remain constant during solver call
subroutine makeStoreOnceRates(n)
use krome_commons
implicit none
real*8,intent(in)::n(nspec)
real*8::small

small = 0d0
rateEvaluateOnce(:) = 0d0

!H -> H+ + E
rateEvaluateOnce(252) = small + (4.6d-1&
    *user_crate)

!HE -> HE+ + E
rateEvaluateOnce(253) = small + (5.d-1&
    *user_crate)

!O -> O+ + E
rateEvaluateOnce(254) = small + (2.8d0&
    *user_crate)

!CO -> C + O
rateEvaluateOnce(255) = small + (5d0&
    *user_crate)

!CO -> CO+ + E
rateEvaluateOnce(256) = small + (3d0&
    *user_crate)

!C2 -> C + C
rateEvaluateOnce(257) = small + (2.37d2&
    *user_crate)

!H2 -> H + H
rateEvaluateOnce(258) = small + (1d-1&
    *user_crate)

!H2 -> H+ + H-
rateEvaluateOnce(259) = small + (3d-4&
    *user_crate)

!H2 -> H2+ + E
rateEvaluateOnce(260) = small + (9.3d-1&
    *user_crate)

!C -> C+ + E
rateEvaluateOnce(261) = small + (1.02d3&
    *user_crate)

!CH -> C + H
rateEvaluateOnce(262) = small + (7.3d2&
    *user_crate)

!O2 -> O + O
rateEvaluateOnce(263) = small + (7.5d2&
    *user_crate)

!O2 -> O2+ + E
rateEvaluateOnce(264) = small + (1.17d2&
    *user_crate)

!OH -> O + H
rateEvaluateOnce(265) = small + (5.1d2&
    *user_crate)

!CH2 -> CH2+ + E
rateEvaluateOnce(266) = small + (5d2&
    *user_crate)

!H2O -> OH + H
rateEvaluateOnce(267) = small + (9.7d2&
    *user_crate)

!HCO -> CO + H
rateEvaluateOnce(268) = small + (4.21d2&
    *user_crate)

!HCO -> HCO+ + E
rateEvaluateOnce(269) = small + (1.17d3&
    *user_crate)

!H2 -> H + H+ + E
rateEvaluateOnce(270) = small + (9.3d-1&
    *user_crate)

end subroutine makeStoreOnceRates

end module krome_tabs

!############### MODULE ##############
module KROME_coolingGH
! ------------------------------------------------------------
!
! This module: Cooling and Heating Functions, table
! reader for Gnedin and Hollon cooling tables
! Language:  Fortran 77
!
!  UPDATED by Troels Haugboelle to Fortran 90, explicit kinds,
!  and encapsulated in a module
!
!  Copyright (c) 2012 Nick Gnedin
!  All rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions
!  are met:
!
!  Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!
!  Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!
!  Neither the name of Nick Gnedin nor the names of any contributors
!  may be used to endorse or promote products derived from this software
!  without specific prior written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR
!  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
!  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ------------------------------------------------------------
private
!  Table dimensions
integer,parameter::NT=81, NX=13, NP1=24, NP2=21, NP3=16, ND=3789
!  Number of components per T-D bin
integer,parameter:: NC=6, NICH=12, NRCH=13
!  Mode of table lookup
integer::mode
!  Boundaries
integer::np(3)
!  Data and index blocks
integer:: indx(NP1,NP2,NP3)
real*4::data(NC,NT,NX,ND)
! Indices
real*8::altval(NT), altmin, altstp, &
    xval(NX), xmin, xmax, xstp, &
    qmin(3), qmax(3), qstp(3)
public :: frtInitCF, frtCFCache, frtCFGetLn, frtGetCF
public :: NICH, NRCH
contains

subroutine frtInitCF(m,fname)
implicit none
integer                      :: m
character(len=*), intent(in) :: fname
!  Internally used unit number
integer, parameter :: IOCF=97
integer :: i, j, k, id, ix, it, ic, lt, ld, lp1, lp2, lp3, lp4, lx
real*4 :: q1, q2
real*4 :: altval4(NT), xmin4, xmax4, qmin4(3), qmax4(3)

mode = m

open(unit=IOCF, file=fname, status='old', form='unformatted', err=100)
read(IOCF,err=100) lt, ld, lp1, lp2, lp3, lp4, &
    (qmin4(j),j=1,3), q1, (qmax4(j),j=1,3), q2, lx, xmin4, xmax4

if(lt.ne.NT .or. ld.ne.ND .or. lx.ne.NX .or. lp1.ne.NP1 .or. &
    lp2.ne.NP2 .or. lp3.ne.NP3 .or. lp4.ne.1 .or. ld.eq.0) then
write(0,*) 'RT::InitCF: fatal error, corrupted table:'
write(0,*) '> NT= in file: ', lt, ' in code: ', NT
write(0,*) '> NX= in file: ', lx, ' in code: ', NX
write(0,*) '> ND= in file: ', ld, ' in code: ', ND
write(0,*) '> NP1= in file: ', lp1, ' in code: ', NP1
write(0,*) '> NP2= in file: ', lp2, ' in code: ', NP2
write(0,*) '> NP3= in file: ', lp3, ' in code: ', NP3
write(0,*) '> NP4= in file: ', lp4, ' in code: ', 1
close(IOCF)
m = -1
stop
end if

qmin=qmin4; qmax=qmax4; xmin=xmin4; xmax=xmax4

np(1) = lp1
np(2) = lp2
np(3) = lp3

do i=1,3
if(np(i) .gt. 1) then
qstp(i) = (qmax(i)-qmin(i))/(np(i)-1)
else
qstp(i) = 1.0
end if
end do

xstp = (xmax-xmin)/(NX-1)
do i=1,NX
xval(i) = xmin + xstp*(i-1)
end do

read(IOCF,err=100) (altval4(i),i=1,NT)
altval = altval4
!  Internally use natural log
do i=1,NT
altval(i) = altval(i)*log(10.0)
end do
altmin = altval(1)
altstp = altval(2) - altval(1)

read(IOCF,err=100) (((indx(i,j,k),i=1,lp1),j=1,lp2),k=1,lp3)

do id=1,ld
read(IOCF,err=100) (((data(ic,it,ix,id),ic=1,NC),it=1,NT),ix=1,NX)
end do

do id=1,ND
do ix=1,NX
do it=1,NT
do ic=1,NC
  data(ic,it,ix,id) = log(1d-37+abs(data(ic,it,ix,id)))
end do
end do
end do
end do

close(IOCF)

if(.false.) then
write(6,*) 'RT::InitCF: Table size = ', &
    NC*(NT*NX*ND/256/1024) + (NP1*NP2*NP3/256/1024), ' MB'
endif

m = 0

return

100   m = -1

end subroutine frtInitCF
!
!  Decode the interpolated function
!
!#define IXL      ich(3)
!#define IXU      ich(4)
!#define IPP(j)   ich(4+j)
!#define WXL      rch(6)
!#define WXU      rch(7)
!#define WPL(j)   rch(7+j)
!#define WPS(j)   rch(10+j)
!
subroutine frtCFPick(it,ich,rch,cfun,hfun)
implicit none
integer :: it
integer :: ich(:)
real*8 :: rch(:)
real*8 :: cfun, hfun
!
integer :: ic, j
real*8 :: v(NC), q(8), a0, a1, a2, Z
!
do ic=1,NC
do j=1,8
q(j) = rch(6)*data(ic,it,ich(3),ich(4+j)) + &
    rch(7)*data(ic,it,ich(4),ich(4+j))
end do
v(ic) = exp( &
    rch(10)*(rch(9)*(rch(8)*q( 1)+rch(11)*q( 2))+   &
    rch(12)*(rch(8)*q( 3)+rch(11)*q( 4))) + &
    rch(13)*(rch(9)*(rch(8)*q( 5)+rch(11)*q( 6))+   &
    rch(12)*(rch(8)*q( 7)+rch(11)*q( 8))))

end do

a0 = v(1)
a1 = v(2)
a2 = v(3)
v(2) = 2*a1 - 0.5*a2 - 1.5*a0
v(3) = 0.5*(a0+a2) - a1

a0 = v(4)
a1 = v(5)
a2 = v(6)
v(5) = 2*a1 - 0.5*a2 - 1.5*a0
v(6) = 0.5*(a0+a2) - a1

Z = rch(5)

if(mode .eq. 1) then
cfun = (Z*v(3)+v(2))*Z
hfun = (Z*v(6)+v(5))*Z
else
cfun = (Z*v(3)+v(2))*Z + v(1)
hfun = (Z*v(6)+v(5))*Z + v(4)
end if

end subroutine frtCFPick

!***********************
!  Cache some table information into arrays iCache and rCache
subroutine frtCFCache(den,Z,Plw,Ph1,Pg1,Pc6,ich,rch,ierr)
implicit none
real*8, intent(in)   :: den, Z, Plw, Ph1, Pg1, Pc6
real*8, dimension(:) :: rch
integer, dimension(:) :: ich
integer :: ierr
!
real*8  :: q(3), qh1, qg1, qc6, dl, w
integer :: j, il(3), is(3)

!  Convert from nb to nH from Cloudy models
dl = max(1.0e-10,den*(1.0-0.02*Z)/1.4)
ierr = 0

if(Plw .gt. 0.0) then
qh1 = log10(1.0e-37+Ph1/Plw)
qg1 = log10(1.0e-37+Pg1/Plw)
qc6 = log10(1.0e-37+Pc6/Plw)

!q(1) = log10(1.0e-37+Plw/dl) ! Should this be n_H or n_baryon ???
q(1) = log10(1.0e-37+Plw/den)
q(2) = 0.263*qc6 + 0.353*qh1 + 0.923*qg1
q(3) = 0.976*qc6 - 0.103*qh1 - 0.375*qg1

!  qmin, qstp, etc are boundaries of cells, not their centers
do j=1,3
w = 0.5 + (q(j)-qmin(j))/qstp(j)
il(j) = int(w) + 1
if(w .gt. il(j)-0.5) then
is(j) = il(j) + 1
else
is(j) = il(j) - 1
endif
rch(10+j) = abs(il(j)-0.5-w)
rch(7+j) = 1 - rch(10+j)

if(np(j) .gt. 1) then
if(max(il(j),is(j)) .gt. np(j)) ierr =  j
if(min(il(j),is(j)) .lt.     1) ierr = -j
endif

if(il(j) .lt. 1) il(j) = 1
if(is(j) .lt. 1) is(j) = 1
if(il(j) .gt. np(j)) il(j) = np(j)
if(is(j) .gt. np(j)) is(j) = np(j)
enddo

else

ierr = -1
do j=1,3
il(j) = 1
is(j) = 1
rch(7+j) = 1
rch(10+j) = 0
enddo

endif

!  Density interpolation is still CIC
w = (log10(dl)-xmin)/xstp
ich(3) = int(w) + 1
if(ich(3) .lt.  1) ich(3) = 1
if(ich(3) .ge. NX) ich(3) = NX-1
ich(4) = ich(3) + 1
rch(6) = max(0.0,min(1.0,ich(3)-w))
rch(7) = 1.0 - rch(6)

!  Do not forget C-to-F77 index conversion
ich(5) = 1 + indx(il(1),il(2),il(3))
ich(6) = 1 + indx(is(1),il(2),il(3))
ich(7) = 1 + indx(il(1),is(2),il(3))
ich(8) = 1 + indx(is(1),is(2),il(3))
ich(9) = 1 + indx(il(1),il(2),is(3))
ich(10) = 1 + indx(is(1),il(2),is(3))
ich(11) = 1 + indx(il(1),is(2),is(3))
ich(12) = 1 + indx(is(1),is(2),is(3))

rch(5) = Z

!  Clear temperature cache
ich(1) = 0
ich(2) = 0
end subroutine frtCFCache

! Get the cooling and heating functions for given
!  ln(T) from the cached dat
subroutine frtCFGetLn(alt,ich,rch,cfun,hfun)
real*8 :: alt
real*8,  dimension(:) :: rch
integer, dimension(:) :: ich
real*8 :: cfun, hfun
real*8, dimension(NC) :: v
integer      :: il, iu
real*8 :: ql, qu
!
il = int((alt-altmin)/altstp*0.99999) + 1
if(il .lt.  1) il = 1
if(il .ge. NT) il = NT-1
iu = il + 1
ql = max(0.0,min(1.0,(altval(iu)-alt)/altstp))
qu = 1.0 - ql
!
!  Shift cache lines as needed
!
if(ich(1) .eq. iu) then
ich(1) = 0
ich(2) = iu
rch(3) = rch(1)
rch(4) = rch(2)
endif
if(ich(2) .eq. il) then
ich(1) = il
ich(2) = 0
rch(1) = rch(3)
rch(2) = rch(4)
endif

!  Update the cache
if(ich(1) .ne. il) then
ich(1) = il
call frtCFPick(il,ich,rch,cfun,hfun)
rch(1) = cfun
rch(2) = hfun
endif
if(ich(2) .ne. iu) then
ich(2) = iu
call frtCFPick(iu,ich,rch,cfun,hfun)
rch(3) = cfun
rch(4) = hfun
endif

cfun = ql*rch(1) + qu*rch(3)
hfun = ql*rch(2) + qu*rch(4)
end subroutine frtCFGetLn

!  Get the cooling and heating functions for T
subroutine frtGetCF(tem,den,Z,Plw,Ph1,Pg1,Pc6,cfun,hfun,ierr)
implicit none
real*8 :: tem, den, Z, Plw, Ph1, Pg1, Pc6, cfun, hfun
integer :: ierr
! Cache arrays
real*8, dimension(NRCH) :: rch
integer,      dimension(NICH) :: ich

call frtCFCache(den,Z,Plw,Ph1,Pg1,Pc6,ich,rch,ierr)
call frtCFGetLn(log(max(1.0,tem)),ich,rch,cfun,hfun)

if (ierr > 0) then
print *,'frtGetCF: Problems with caching Gnedin and Hollon cooling table.'
print *,'frtGetCF:  Stopping. Error code :', ierr
print *,'frtGetCF: Input variables ntot, Tgas, PLW, PHI, PHeI, PCVI: ', den, tem, Plw, Ph1, Pg1, Pc6
stop
end if

end subroutine frtGetCF

end module KROME_coolingGH


!############### MODULE ##############
module KROME_cooling
! *************************************************************
!  This file has been generated with:
!  KROME 14.08.dev on 2018-10-24 12:49:10
!  Changeset b21f657
!  see http://kromepackage.org
!
!  Written and developed by Tommaso Grassi and Stefano Bovino
!
!  Contributors:
!  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
!  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
!  E.Tognelli
!  KROME is provided "as it is", without any warranty.
! *************************************************************
integer,parameter::coolTab_n=int(1e2)
integer,parameter::nZrate=50
real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
real*8::pop_level_CI(3)
real*8::pop_level_CII(2)
real*8::pop_level_OI(3)
!$omp threadprivate(pop_level_CI)
!$omp threadprivate(pop_level_CII)
!$omp threadprivate(pop_level_OI)
contains

!*******************
function cooling(n,inTgas)
use krome_commons
implicit none
real*8::n(:),inTgas,cooling,Tgas

Tgas = inTgas
cooling = sum(get_cooling_array(n(:),Tgas))

end function cooling

!*******************************
function get_cooling_array(n, Tgas)
use krome_commons
implicit none
real*8::n(:), Tgas
real*8::get_cooling_array(ncools),cools(ncools)
real*8::f1,f2,smooth

f1 = 1d0
f2 = 1d0

!returns cooling in erg/cm3/s
cools(:) = 0d0

cools(idx_cool_H2) = cooling_H2(n(:), Tgas)

cools(idx_cool_dust) = cooling_dust(n(:), Tgas)

cools(idx_cool_CO) = cooling_CO(n(:), Tgas)

cools(idx_cool_atomic) = cooling_Atomic(n(:), Tgas)

cools(idx_cool_ff) = cooling_ff(n(:), Tgas)

cools(idx_cool_cont) = cooling_continuum(n(:), Tgas)

!this parameter controls the smoothness of the
! merge between the two cooling functions
smooth = 1./500.

!smoothing functions | f1+f2=1
f1 = (tanh(smooth*(Tgas-2.5d3))+1.d0)*0.5d0
f2 = (tanh(smooth*(-Tgas+2.5d3))+1.d0)*0.5d0

! Only compute cooling terms if it really matters
if (f2 > 1d-20) then

cools(idx_cool_Z) = f2 * ( cooling_Z(n(:), Tgas)  )

endif

if (f1 > 1d-20) cools(idx_cool_GH) = f1 * cooling_GH(n(:), Tgas)

cools(idx_cool_custom) = cooling_custom(n(:),Tgas)

get_cooling_array(:) = cools(:)

end function get_cooling_array

!***************************
!CO cooling: courtesy of K.Omukai (Nov2014)
! method: Neufeld+Kaufman 1993 (bit.ly/1vnjcXV, see eqn.5).
! see also Omukai+2010 (bit.ly/1HIaGcn)
! H and H2 collisions
function cooling_CO(n,inTgas)
use krome_commons
use krome_subs
use krome_getphys
implicit none
integer,parameter::imax=coolCOn1
integer,parameter::jmax=coolCOn2
integer,parameter::kmax=coolCOn3
integer::i,j,k
real*8,parameter::eps=1d-5
real*8::cooling_CO,n(:),inTgas
real*8::v1,v2,v3,prev1,prev2,cH
real*8::vv1,vv2,vv3,vv4,vv12,vv34,xLd
real*8::x1(imax),x2(jmax),x3(kmax)
real*8::ixd1(imax-1),ixd2(jmax-1),ixd3(kmax-1)
real*8::v1min,v1max,v2min,v2max,v3min,v3max

!local copy of limits
v1min = coolCOx1min
v1max = coolCOx1max
v2min = coolCOx2min
v2max = coolCOx2max
v3min = coolCOx3min
v3max = coolCOx3max

!local copy of variables arrays
x1(:) = coolCOx1(:)
x2(:) = coolCOx2(:)
x3(:) = coolCOx3(:)

ixd1(:) = coolCOixd1(:)
ixd2(:) = coolCOixd2(:)
ixd3(:) = coolCOixd3(:)

!local variables
v3 = num2col(n(idx_CO),n(:)) !CO column density
cH = n(idx_H) + n(idx_H2)
if (n(idx_H) < 0. .or. n(idx_H2) < 0.) then
cH = n_global(idx_H) + n_global(idx_H2)
! Set red flag if not due to small excursion in n_H2
! Only relevant for low temperatures, bc otherwise H is ionised
! And CO cooling only relevant at high densities
if (abs(n(idx_H2)) / (abs(n(idx_H)) + 1d-40) > 1d-6 .and. &
    inTgas < 5d4 .and. &
    abs(n(idx_H)) > abs(n(idx_Hj)) .and. sum(n(1:nmols)) > 100.) &
    red_flag = ibset(red_flag,5)
endif

v2 = cH
v1 = inTgas !Tgas

!logs of variables
v1 = log10(v1)
v2 = log10(v2)
v3 = log10(v3)

!default value erg/s/cm3
cooling_CO = 0d0

!check limits
if(v1>=v1max) v1 = v1max*(1d0-eps)
if(v2>=v2max) v2 = v2max*(1d0-eps)
if(v3>=v3max) v3 = v3max*(1d0-eps)

if(v1<v1min) return
if(v2<v2min) return
if(v3<v3min) return

!gets position of variable in the array
i = (v1-v1min)*coolCOdvn1+1
j = (v2-v2min)*coolCOdvn2+1
k = (v3-v3min)*coolCOdvn3+1

!precompute shared variables
prev1 = (v1-x1(i))*ixd1(i)
prev2 = (v2-x2(j))*ixd2(j)

!linear interpolation on x1 for x2,x3
vv1 = prev1 * (coolCOy(k,j,i+1) - &
    coolCOy(k,j,i)) + coolCOy(k,j,i)
!linear interpolation on x1 for x2+dx2,x3
vv2 = prev1 * (coolCOy(k,j+1,i+1) - &
    coolCOy(k,j+1,i)) + coolCOy(k,j+1,i)
!linear interpolation on x2 for x3
vv12 = prev2 * (vv2 - vv1) + vv1

!linear interpolation on x1 for x2,x3+dx3
vv3 = prev1 * (coolCOy(k+1,j,i+1) - &
    coolCOy(k+1,j,i)) + coolCOy(k+1,j,i)
!linear interpolation on x1 for x2+dx2,x3+dx3
vv4 = prev1 * (coolCOy(k+1,j+1,i+1) - &
    coolCOy(k+1,j+1,i)) + coolCOy(k+1,j+1,i)
!linear interpolation on x2 for x3+dx3
vv34 = prev2 * (vv4 - vv3) + vv3

!linear interpolation on x3
xLd = (v3-x3(k))*ixd3(k)*(vv34 - &
    vv12) + vv12

!CO cooling in erg/s/cm3
cooling_CO = 1d1**xLd * cH * n(idx_CO)

end function cooling_CO

!************************
subroutine init_coolingCO()
use krome_commons
implicit none
integer::ios,iout(3),i
real*8::rout(4)

if (krome_mpi_rank<=1) print *,"load CO cooling..."
open(33,file="coolCO.dat",status="old",iostat=ios)
!check if file exists
if(ios.ne.0) then
print *,"ERROR: problems loading coolCO.dat!"
stop
end if

do
read(33,*,iostat=ios) iout(:),rout(:) !read line
if(ios<0) exit !eof
if(ios/=0) cycle !skip blanks
coolCOx1(iout(1)) = rout(1)
coolCOx2(iout(2)) = rout(2)
coolCOx3(iout(3)) = rout(3)
coolCOy(iout(3),iout(2),iout(1)) = rout(4)
end do

!store inverse of the differences
! to speed up interpolation
do i=1,coolCOn1-1
coolCOixd1(i) = 1d0/(coolCOx1(i+1)-coolCOx1(i))
end do
do i=1,coolCOn2-1
coolCOixd2(i) = 1d0/(coolCOx2(i+1)-coolCOx2(i))
end do
do i=1,coolCOn3-1
coolCOixd3(i) = 1d0/(coolCOx3(i+1)-coolCOx3(i))
end do

coolCOx1min = minval(coolCOx1)
coolCOx1max = maxval(coolCOx1)
coolCOx2min = minval(coolCOx2)
coolCOx2max = maxval(coolCOx2)
coolCOx3min = minval(coolCOx3)
coolCOx3max = maxval(coolCOx3)

coolCOdvn1 = (coolCOn1-1)/(coolCOx1max-coolCOx1min)
coolCOdvn2 = (coolCOn2-1)/(coolCOx2max-coolCOx2min)
coolCOdvn3 = (coolCOn3-1)/(coolCOx3max-coolCOx3min)

end subroutine init_coolingCO

!*****************************
function cooling_custom(n,Tgas)
use krome_commons
use krome_subs
use krome_constants
implicit none
real*8::n(:),Tgas,cooling_custom

cooling_custom = 0d0

end function cooling_custom

!**********************************
function kpla(n,Tgas)
!Planck opacity mean fit (Lenzuni+1996)
!only denisity dependent (note that the
! fit provided by Lenzuni is wrong)
! valid for T<3e3 K
!use krome_subs
use krome_commons
use krome_getphys
implicit none
real*8::kpla,rhogas,Tgas,n(:),y
real*8::a0,a1,m(nspec)

m(:) = get_mass()
rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3

kpla = 0.d0
!opacity is zero under 1e-12 g/cm3
if(rhogas<1d-12) return

!fit coefficients
a0 = 1.000042d0
a1 = 2.14989d0

!log density cannot exceed 0.5 g/cm3
y = log10(min(rhogas,0.5d0))

kpla = 1d1**(a0*y + a1) !fit density only

end function kpla

!*****************************
function coolingChem(n,Tgas)
implicit none
real*8::coolingChem,n(:),Tgas

!note that this function is a dummy.
! For chemical cooling you should see
! heatingChem function in krome_heating.f90

coolingChem = 0.d0

end function coolingChem

!**********************************
function cooling_Continuum(n,Tgas)
!cooling from continuum for a thin gas (no opacity)
!see Omukai+2000 for details
use krome_commons
use krome_constants
use krome_subs
use krome_getphys
implicit none
real*8::n(:),Tgas,cooling_Continuum,kgas,rhogas
real*8::lj,tau,beta,m(nspec)

m(:) = get_mass()
rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3
kgas = kpla(n(:),Tgas) !planck opacity cm2/g (Omukai+2000)
lj = get_jeans_length(n(:), Tgas) !cm
tau = lj * kgas * rhogas + 1d-40 !opacity
beta = min(1.d0,tau**(-2)) !beta escape (always <1.)
cooling_Continuum = 4.d0 * stefboltz_erg * Tgas**4 &
    * kgas * rhogas * beta !erg/s/cm3

end function cooling_Continuum

!****************************
function cooling_dust(n,Tgas)
use krome_commons
use krome_subs
use krome_fit
use krome_getphys
implicit none
real*8::n(:),Tgas,ntot,cooling_dust,coolFit
real*8::logn,logt

ntot = sum(n(1:nmols))
Tgas = n(idx_Tgas)

logn = log10(ntot)
logt = log10(Tgas)
!cooling fit from tables
coolFit = fit_anytab2D_linlog(dust_tab_ngas(:), dust_tab_Tgas(:), &
    dust_tab_cool(:,:), dust_mult_ngas, dust_mult_Tgas, &
    logn, logt)

cooling_dust = get_mu(n) * coolFit * ntot * ntot

end function cooling_dust

!*****************
!sigmoid function with x0 shift and s steepness
function sigmoid(x,x0,s)
implicit none
real*8::sigmoid,x,x0,s

sigmoid = 1d1/(1d1+exp(-s*(x-x0)))

end function sigmoid

!*******************
!window function for H2 cooling to smooth limits
function wCool(logTgas,logTmin,logTmax)
implicit none
real*8::wCool,logTgas,logTmin,logTmax,x

x = (logTgas-logTmin)/(logTmax-logTmin)
wCool = 1d1**(2d2*(sigmoid(x,-2d-1,5d1)*sigmoid(-x,-1.2d0,5d1)-1d0))
if(wCool<1d-199) wCool = 0d0
if(wCool>1d0) then
print *,"ERROR: wCool>1"
stop
end if

end function wCool

!ALL THE COOLING FUNCTIONS ARE FROM GLOVER & ABEL, MNRAS 388, 1627, 2008
!FOR LOW DENSITY REGIME: CONSIDER AN ORTHO-PARA RATIO OF 3:1
!UPDATED TO THE DATA REPORTED BY GLOVER 2015, MNRAS
!EACH SINGLE FUNCTION IS IN erg/s
!FINAL UNITS = erg/cm3/s
!*******************************
function cooling_H2(n, Tgas)
use krome_commons
use krome_subs
use krome_getphys
real*8::n(:),Tgas
real*8::temp,logt3,logt,cool,cooling_H2,T3
real*8::LDL,HDLR,HDLV,HDL
real*8::logt32,logt33,logt34,logt35,logt36,logt37,logt38
real*8::dump14,fH2H,fH2e,fH2H2,fH2Hp,fH2He,w14,w24
integer::i
character*16::names(nspec)

temp = Tgas
cooling_H2 = 0d0
!if(temp<2d0) return

T3 = temp * 1.d-3
logt3 = log10(T3)
logt = log10(temp)
cool = 0d0

logt32 = logt3 * logt3
logt33 = logt32 * logt3
logt34 = logt33 * logt3
logt35 = logt34 * logt3
logt36 = logt35 * logt3
logt37 = logt36 * logt3
logt38 = logt37 * logt3

w14 = wCool(logt, 1d0, 4d0)
w24 = wCool(logt, 2d0, 4d0)

!//H2-H
if(temp<=1d2) then
fH2H = 1.d1**(-16.818342D0 +3.7383713D1*logt3 &
    +5.8145166D1*logt32 +4.8656103D1*logt33 &
    +2.0159831D1*logt34 +3.8479610D0*logt35 )*n(idx_H)
elseif(temp>1d2 .and. temp<=1d3) then
fH2H = 1.d1**(-2.4311209D1 +3.5692468D0*logt3 &
    -1.1332860D1*logt32 -2.7850082D1*logt33 &
    -2.1328264D1*logt34 -4.2519023D0*logt35)*n(idx_H)
elseif(temp>1.d3.and.temp<=6d3) then
fH2H = 1d1**(-2.4311209D1 +4.6450521D0*logt3 &
    -3.7209846D0*logt32 +5.9369081D0*logt33 &
    -5.5108049D0*logt34 +1.5538288D0*logt35)*n(idx_H)
else
fH2H = 1.862314467912518E-022*wCool(logt,1d0,log10(6d3))*n(idx_H)
end if
cool = cool + fH2H

!//H2-Hp
if(temp>1d1.and.temp<=1d4) then
fH2Hp = 1d1**(-2.2089523d1 +1.5714711d0*logt3 &
    +0.015391166d0*logt32 -0.23619985d0*logt33 &
    -0.51002221d0*logt34 +0.32168730d0*logt35)*n(idx_Hj)
else
fH2Hp = 1.182509139382060E-021*n(idx_Hj)*w14
endif
cool = cool + fH2Hp

!//H2-H2
fH2H2 = w24*1d1**(-2.3962112D1 +2.09433740D0*logt3 &
    -.77151436D0*logt32 +.43693353D0*logt33 &
    -.14913216D0*logt34 -.033638326D0*logt35)*n(idx_H2) !&
    cool = cool + fH2H2

!//H2-e
fH2e = 0d0
if(temp<=5d2) then
fH2e = 1d1**(min(-2.1928796d1 + 1.6815730d1*logt3 &
    +9.6743155d1*logt32 +3.4319180d2*logt33 &
    +7.3471651d2*logt34 +9.8367576d2*logt35 &
    +8.0181247d2*logt36 +3.6414446d2*logt37 &
    +7.0609154d1*logt38,3d1))*n(idx_e)
elseif(temp>5d2)  then
fH2e = 1d1**(-2.2921189D1 +1.6802758D0*logt3 &
    +.93310622D0*logt32 +4.0406627d0*logt33 &
    -4.7274036d0*logt34 -8.8077017d0*logt35 &
    +8.9167183*logt36 + 6.4380698*logt37 &
    -6.3701156*logt38)*n(idx_e)
end if
cool = cool + fH2e*w24

!//H2-He
if(temp>1d1.and.temp<=1d4)then
fH2He = 1d1**(-2.3689237d1 +2.1892372d0*logt3&
    -.81520438d0*logt32 +.29036281d0*logt33 -.16596184d0*logt34 &
    +.19191375d0*logt35)*n(idx_He)
else
fH2He = 1.002560385050777E-022*n(idx_He)*w14
endif
cool = cool + fH2He

!check error
if(cool>1.d30) then
print *,"ERROR: cooling >1.d30 erg/s/cm3"
print *,"cool (erg/s/cm3): ",cool
names(:) = get_names()
do i=1,size(n)
print '(I3,a18,E11.3)',i,names(i),n(i)
end do
stop
end if

!this to avoid negative, overflow and useless calculations below
if(cool<=0d0) then
cooling_H2 = 0d0
return
end if

!high density limit from HM79, GP98 below Tgas = 2d3
!UPDATED USING GLOVER 2015 for high temperature corrections, MNRAS
!IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL
!the following mix of functions ensures the right behaviour
! at low (T<10 K) and high temperatures (T>2000 K) by
! using both the original Hollenbach and the new Glover data
! merged in a smooth way.
if(temp.lt.2d3)then
HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
    3.e-24*exp(-0.51/t3)) !erg/s
HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
HDL  = HDLR + HDLV !erg/s
elseif(temp>=2d3 .and. temp<=1d4)then
HDL = 1d1**(-2.0584225d1 + 5.0194035*logt3 &
    -1.5738805*logt32 -4.7155769*logt33 &
    +2.4714161*logt34 +5.4710750*logt35 &
    -3.9467356*logt36 -2.2148338*logt37 &
    +1.8161874*logt38)
else
dump14 = 1d0 / (1d0 + exp(min((temp-3d4)*2d-4,3d2)))
HDL = 5.531333679406485E-019*dump14
endif

LDL = cool !erg/s
if (HDL==0.) then
cooling_H2 = 0.d0
else
cooling_H2 = n(idx_H2)/(1.d0/HDL+1.d0/LDL)  !erg/cm3/s
endif

end function cooling_H2

!Atomic COOLING  Cen ApJS, 78, 341, 1992
!UNITS = erg/s/cm3
!*******************************
function cooling_Atomic(n, Tgas)
use krome_commons
use krome_subs
real*8::Tgas,cooling_atomic,n(:)
real*8::temp,T5,cool

temp = max(Tgas,10d0) !K
T5 = temp/1d5 !K
cool = 0d0 !erg/cm3/s

!COLLISIONAL IONIZATION: H, He, He+, He(2S)
cool = cool+ 1.27d-21*sqrt(temp)/(1.d0+sqrt(T5))&
    *exp(-1.578091d5/temp)*n(idx_e)*n(idx_H)

cool = cool+ 9.38d-22*sqrt(temp)/(1.d0+sqrt(T5))&
    *exp(-2.853354d5/temp)*n(idx_e)*n(idx_He)

cool = cool+ 4.95d-22*sqrt(temp)/(1.d0+sqrt(T5))&
    *exp(-6.31515d5/temp)*n(idx_e)*n(idx_Hej)
cool = cool+ 5.01d-27*temp**(-0.1687)/(1.d0+sqrt(T5))&
    *exp(-5.5338d4/temp)*n(idx_e)**2*n(idx_Hej)

!RECOMBINATION: H+, He+,He2+
cool = cool+ 8.7d-27*sqrt(temp)*(temp/1.d3)**(-0.2)&
    /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hj)

cool = cool+ 1.55d-26*temp**(0.3647)*n(idx_e)*n(idx_Hej)

cool = cool+ 3.48d-26*sqrt(temp)*(temp/1.d3)**(-0.2)&
    /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hejj)

!DIELECTRONIC RECOMBINATION: He
cool = cool+ 1.24d-13*temp**(-1.5)*exp(-4.7d5/temp)&
    *(1.d0+0.3d0*exp(-9.4d4/temp))*n(idx_e)*n(idx_Hej)

!COLLISIONAL EXCITATION:
!H(all n), He(n=2,3,4 triplets), He+(n=2)
cool = cool+ 7.5d-19/(1.d0+sqrt(T5))*exp(-1.18348d5/temp)*n(idx_e)*n(idx_H)

cool = cool+ 9.1d-27*temp**(-.1687)/(1.d0+sqrt(T5))&
    *exp(-1.3179d4/temp)*n(idx_e)**2*n(idx_Hej)
cool = cool+ 5.54d-17*temp**(-.397)/(1.d0+sqrt(T5))&
    *exp(-4.73638d5/temp)*n(idx_e)*n(idx_Hej)

cooling_atomic = max(cool, 0d0)  !erg/cm3/s

end function cooling_Atomic

!**************************
!free-free cooling (bremsstrahlung for all ions)
! using mean Gaunt factor value (Cen+1992)
function cooling_ff(n,Tgas)
use krome_commons
implicit none
real*8::n(:),Tgas,cool,cooling_ff,gaunt_factor,bms_ions

gaunt_factor = 1.5d0 !mean value

!BREMSSTRAHLUNG: all ions
bms_ions = +n(idx_Hj) +n(idx_HEj) +n(idx_Cj) +n(idx_Oj) +4.d0*n(idx_HEjj)
cool = 1.42d-27*gaunt_factor*sqrt(Tgas)&
    *bms_ions*n(idx_e)

cooling_ff = max(cool, 0.d0)  !erg/cm3/s

end function cooling_ff

!******************************
function cooling_GH(n, Tgas)
use krome_commons
use krome_constants, only : ip_mass
use krome_coolingGH
use krome_getphys, only : get_mass, get_metallicityC, get_metallicityO
implicit none
real*8::n(:),Tgas
real*8:: cooling_GH
real*8:: PLW_last=0d0,PHI_last=0d0,PHeI_last=0d0, &
    PCVI_last=0d0,nb_last=0d0
real*8,  save :: rch(NRCH)
integer, save :: ich(NICH)
!$omp threadprivate(PLW_last,PHI_last,PHeI_last,PCVI_last,nb_last,rch,ich)
real*8::log10Tgas,nb,Zmetal,ZAsplund,cfun,hfun
logical, save :: first_call=.true.
integer       :: ierr
! Parameters related to dust destruction at high temperatures
real*8, parameter :: Tdestruct=1e5, smooth=1./(0.1*Tdestruct)
! Check if we have to load table -- only done by one thread
! Do the double if-block so that we only enter the critical region
! if first_call has not happened or is in progress
!----------------------------------
if (first_call) then
!$omp critical
if (first_call) then
! NOTICE ierr is double duty.
! On input: sets mode of tables. =0: full heating nd cooling function. =1: only metals
! On output: error code
ierr = 1
call frtInitCF(ierr,'cf_table.I2.dat')
if(ierr .ne. 0) then
print *,'Error in reading Gnedin and Hollon'
print *,' cooling table data file: ', ierr
print *,'Maybe data file cf_table.I2.dat'
print *,' does not exist in current directory ?'
stop
endif
end if
first_call = .false.
!$omp end critical
end if

nb = sum(n(1:nspec)*get_mass()) * ip_mass

if (Tgas > 1d3) then
!check input values to see if we have to regenerate the cache
if (abs(nb-nb_last) > 1d-6*nb .or. &
    PLW  .ne. PLW_last  .or. &
    PHI  .ne. PHI_last  .or. &
    PHeI .ne. PHeI_last .or. &
    PCVI .ne. PCVI_last) then
! Gas-phase metallicity based on an average of C and O in case they differ wrt dust depletion
Zasplund = 10.**(0.5*(get_metallicityC(n)+get_metallicityO(n))) ! Normalised to Asplund 2009 value of Z=0.0134
Zmetal = Zasplund * 0.0134 / 0.02 ! Metallicity of gas in units of 0.02
! Above a certain temperature, Tdestruct, correct depletion by
! assuming dust-to-gas = 0.01 and all dust sublimated and therefore
! Zmetal = (Zmetal * 0.02 + 0.01) / 0.02
! At Tdestruct = 10^5 K 1% of protons have v > 50 km/s and can directly destruct dust grains
! 50 km/s is enough to destruct dust grains because everything is charged and gyrate around mag fields, which gives extra velocity boost
Zmetal = (Zmetal * 0.02 + &
    0.01 * ((tanh(smooth*(Tgas-Tdestruct))+1.d0)*0.5d0) &
    ) / 0.02
!
call frtCFCache(nb, Zmetal, PLW, PHI, PHeI, PCVI, ich, rch, ierr)
nb_last=nb
PLW_last=PLW
PHI_last=PHI
PHeI_last=PHeI
PCVI_last=PCVI
if (ierr==1 .or. ierr==3) then
print *,'Problems with caching Gnedin and Hollon cooling table.'
print *,' Stopping. Error code :', ierr
print *,'Input variables nb, Tgas, PLW, PHI, PHeI, PCVI: ', nb, Tgas, PLW, PHI, PHeI, PCVI
if(ierr==1) stop
end if
end if

call frtCFGetLn(log(max(1d0,Tgas)),ich,rch,cfun,hfun)
!call frtGetCF(Tgas, nb, 1d0, Plw, Phi, Phei, Pcvi, cfun, hfun, ierr)

else
cfun = 0d0
hfun = 0d0
endif

cooling_GH = nb**2 * (cfun-hfun)

end function cooling_GH

!*********************************************

!function for linear interpolation of f(x), using xval(:)
! and the corresponding yval(:) as reference values
! note: slow function, use only for initializations
function flin(xval,yval,x)
implicit none
real*8::xval(:),yval(:),x,flin
integer::i,n
logical::found
found = .false.
n = size(xval)
x = max(x,xval(1)) !set lower bound
x = min(x,xval(n)) !set upper bound
!loop to find interval (slow)
do i=2,n
if(x.le.xval(i)) then
!linear fit
flin = (yval(i) - yval(i-1)) / (xval(i) - xval(i-1)) * &
    (x - xval(i-1)) + yval(i-1)
found = .true. !found flag
exit
end if
end do
if(.not.found) flin = yval(n)

end function flin

!************************
!dump the level populations in a file
subroutine dump_cooling_pop(Tgas,nfile)
implicit none
integer::nfile,i
real*8::Tgas

!pop_level_CI(3)
do i=1,size(pop_level_CI)
write(nfile,'(a8,I5,3E17.8e3)') "CI", i, Tgas, pop_level_CI(i), sum(pop_level_CI(:))
end do

!pop_level_CII(2)
do i=1,size(pop_level_CII)
write(nfile,'(a8,I5,3E17.8e3)') "CII", i, Tgas, pop_level_CII(i), sum(pop_level_CII(:))
end do

!pop_level_OI(3)
do i=1,size(pop_level_OI)
write(nfile,'(a8,I5,3E17.8e3)') "OI", i, Tgas, pop_level_OI(i), sum(pop_level_OI(:))
end do

write(nfile,*)

end subroutine dump_cooling_pop

!***********************
!metal cooling as in Maio et al. 2007
! loaded from data file
function cooling_Z(n,inTgas)
use krome_commons
use krome_constants
implicit none
real*8::n(:), inTgas, cool, cooling_Z, k(nZrate), Tgas

Tgas = inTgas
k(:) = coolingZ_rate_tabs(Tgas)

cool = 0d0
cool = cool + coolingCI(n(:),inTgas,k(:))
cool = cool + coolingCII(n(:),inTgas,k(:))
cool = cool + coolingOI(n(:),inTgas,k(:))

cooling_Z = cool * boltzmann_erg

end function cooling_Z

!********************************
function coolingZ_rates(inTgas)
use krome_commons
use krome_subs
use krome_fit
implicit none
real*8::inTgas,coolingZ_rates(nZrate),k(nZrate)
real*8::Tgas,invT,logTgas
integer::i
real*8::T2,invTgas,lnT

Tgas = inTgas
invT = 1d0/Tgas
logTgas = log10(Tgas)

T2 = Tgas*1d-2
invTgas = 1d0/Tgas
lnT = log(Tgas)

if(Tgas.ge.1d4) return

!1->0, CI - H
k(1) = 1.6d-10*(T2)**(.14)

!2->0, CI - H
k(2) = 9.2d-11*(T2)**(.26)

!2->1, CI - H
k(3) = 2.9d-10*(T2)**(.26)

!1->0, CI - H+
k(4) = (9.6D-11 -1.8D-14*Tgas +1.9D-18*Tgas**2) *Tgas**(.45)

if(Tgas > 5d3) k(4) = 8.9D-10*Tgas**(.117)

!2->0, CI - H+
k(5) = (3.1D-12 -6.D-16*Tgas +3.9d-20*Tgas**2) *Tgas

if(Tgas > 5d3) k(5) = 2.3D-9*Tgas**(.0965)

!2->1, CI - H+
k(6) = (1.D-10 -2.2D-14*Tgas +1.7D-18*Tgas**2) *Tgas**(.70)

if(Tgas > 5d3) k(6) = 9.2D-9*Tgas**(.0535)

!1->0, CI - e
k(7) = 2.88D-6*Tgas**(-.5)*EXP(-9.25141 -7.73782D-1*lnT +3.61184D-1*lnT**2 -1.50892D-2*lnT**3 -6.56325D-4*lnT**4)

if(Tgas > 1D3) k(7) = 2.88D-6*Tgas**(-.5) *EXP(-4.446D2 -2.27913D2*lnT +4.2595D1*lnT**2 -3.4762*lnT**3 +1.0508D-1*lnT**4)

!2->0, CI - e
k(8) = 1.73D-6*Tgas**(-.5)*EXP(-7.69735 -1.30743*lnT +.697638*lnT**2 -.111338*lnT**3 +.705277D-2*lnT**4)

if(Tgas > 1D3) k(8) = 1.73D-6*Tgas**(-.5)*EXP(3.50609D2 -1.87474D2*lnT +3.61803D1*lnT**2 -3.03283*lnT**3 +9.38138D-2*lnT**4)

!2->1, CI - e
k(9) = 1.73D-6*Tgas**(-.5)*EXP(-7.4387 -.57443*lnT +.358264*lnT**2 -4.18166D-2*lnT**3 +2.35272D-3*lnT**4)

if(Tgas > 1D3) k(9) = 1.73D-6*Tgas**(-.5)*EXP(3.86186D2 -2.02192D2*lnT +3.85049D1*lnT**2 -3.19268*lnT**3 +9.78573D-2*lnT**4)

!1->0, CI - H2or
k(10) = 8.7d-11 -6.6d-11*exp(-Tgas/218.3) + 6.6d-11*exp(-2.*Tgas/218.3)

!2->0, CI - H2or
k(11) = 1.2d-10 -6.1d-11*exp(-Tgas/387.3)

!2->1, CI - H2or
k(12) = 2.9d-10 -1.9d-10*exp(-Tgas/348.9)

!1->0, CI - H2pa
k(13) = 7.9D-11 -8.7D-11*EXP(-Tgas/126.4) + 1.3D-10*EXP(-2.*Tgas/126.4)

!2->0, CI - H2pa
k(14) = 1.1D-10 -8.6D-11*EXP(-Tgas/223.) + 8.7D-11*EXP(-2.*Tgas/223.)

!2->1, CI - H2pa
k(15) = 2.7D-10 -2.6D-10*EXP(-Tgas/250.7) + 1.8D-10*EXP(-2.*Tgas/250.7)

!1->0, OI - H
k(16) = 9.2D-11*(T2)**(.67)

!2->0, OI - H
k(17) = 4.3D-11*(T2)**(.80)

!2->1, OI - H
k(18) = 1.1D-10*(T2)**(.44)

!1->0, OI - H+
k(19) = 6.38D-11*Tgas**(.4)

if(Tgas > 194.) k(19) = 7.75D-12*Tgas**(.8)

if(Tgas > 3686.) k(19) = 2.65D-10*Tgas**(.37)

!2->0, OI - H
k(20) = 6.1D-13*Tgas**(1.1)

if(Tgas > 511.) k(20) = 2.12D-12*Tgas**(.9)

if(Tgas > 7510.) k(20) = 4.49D-10*Tgas**(.3)

!2->1, OI - H
k(21) = 2.03D-11*Tgas**(.56)

if(Tgas > 2090.) k(21) = 3.43D-10*Tgas**(.19)

!1->0, OI - e
k(22) = 5.12D-10*Tgas**(-.075)

!2->0, OI - e
k(23) = 4.86D-10*Tgas**(-.026)

!2->1, OI - e
k(24) = 1.08D-14*Tgas**(.926)

!1->0, CII - e
k(25) = 2.8D-7*(Tgas/1d2)**(-.5)

!1->0, CII - H
k(26) = 8D-10*(Tgas/1d2)**(.07)

!2<-0, CI - H2or
k(27) = k(11) * 5.0d0 * exp(-6.300000d+01 * invT)

!1<-0, CI - e
k(28) = k(7) * 3.0d0 * exp(-2.400000d+01 * invT)

!2<-0, CI - e
k(29) = k(8) * 5.0d0 * exp(-6.300000d+01 * invT)

!1<-0, CI - H2or
k(30) = k(10) * 3.0d0 * exp(-2.400000d+01 * invT)

!1<-0, CI - H2pa
k(31) = k(13) * 3.0d0 * exp(-2.400000d+01 * invT)

!2<-0, CI - H2pa
k(32) = k(14) * 5.0d0 * exp(-6.300000d+01 * invT)

!1<-0, CI - H+
k(33) = k(4) * 3.0d0 * exp(-2.400000d+01 * invT)

!2<-1, CI - e
k(34) = k(9) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

!2<-1, CI - H
k(35) = k(3) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

!2<-0, CI - H
k(36) = k(2) * 5.0d0 * exp(-6.300000d+01 * invT)

!2<-1, CI - H+
k(37) = k(6) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

!2<-1, CI - H2pa
k(38) = k(15) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

!2<-1, CI - H2or
k(39) = k(12) * 1.66666666667d0 * exp(-3.900000d+01 * invT)

!1<-0, CI - H
k(40) = k(1) * 3.0d0 * exp(-2.400000d+01 * invT)

!2<-0, CI - H+
k(41) = k(5) * 5.0d0 * exp(-6.300000d+01 * invT)

!1<-0, CII - e
k(42) = k(25) * 2.0d0 * exp(-9.120000d+01 * invT)

!1<-0, CII - H
k(43) = k(26) * 2.0d0 * exp(-9.120000d+01 * invT)

!1<-0, OI - H
k(44) = k(16) * 0.6d0 * exp(-2.300000d+02 * invT)

!1<-0, OI - e
k(45) = k(22) * 0.6d0 * exp(-2.300000d+02 * invT)

!2<-1, OI - e
k(46) = k(24) * 0.333333333333d0 * exp(-1.000000d+02 * invT)

!2<-1, OI - H
k(47) = k(21) * 0.333333333333d0 * exp(-1.000000d+02 * invT)

!2<-0, OI - H
k(48) = k(20) * 0.2d0 * exp(-3.300000d+02 * invT)

!1<-0, OI - H+
k(49) = k(19) * 0.6d0 * exp(-2.300000d+02 * invT)

!2<-0, OI - e
k(50) = k(23) * 0.2d0 * exp(-3.300000d+02 * invT)

coolingZ_rates(:) = k(:)

!check rates > 1
if(maxval(k)>1d0) then
print *,"ERROR: found rate >1d0 in coolingZ_rates!"
print *," Tgas =",Tgas
do i=1,nZrate
if(k(i)>1d0) print *,i,k(i)
end do
stop
end if

!check rates <0
if(minval(k)<0d0) then
print *,"ERROR: found rate <0d0 in coolingZ_rates!"
print *," Tgas =",Tgas
do i=1,nZrate
if(k(i)<0d0) print *,i,k(i)
end do
stop
end if

end function coolingZ_rates

!**********************
function coolingZ_rate_tabs(inTgas)
use krome_commons
implicit none
real*8::inTgas,Tgas,coolingZ_rate_tabs(nZrate),k(nZrate)
integer::idx,j
Tgas = inTgas

idx = (log10(Tgas)-coolTab_logTlow) * inv_coolTab_idx + 1

idx = max(idx,1)
idx = min(idx,coolTab_n-1)

do j=1,nZrate
k(j) = (Tgas-coolTab_T(idx)) * inv_coolTab_T(idx) * &
    (coolTab(j,idx+1)-coolTab(j,idx)) + coolTab(j,idx)
k(j) = max(k(j), 0d0)
end do

coolingZ_rate_tabs(:) = k(:)

end function coolingZ_rate_tabs

!**********************
subroutine coolingZ_init_tabs()
use krome_commons
implicit none
integer::j,jmax,idx
real*8::Tgas,Tgasold

jmax = coolTab_n !size of the cooling tables (number of saples)

!note: change upper and lower limit for rate tables here
coolTab_logTlow = log10(2d0)
coolTab_logTup = log10(1d8)

!pre compute this value since used jmax times
inv_coolTab_idx = (jmax-1) / (coolTab_logTup-coolTab_logTlow)

!loop over the jmax interpolation points
do j=1,jmax
!compute Tgas for the given point
Tgas = 1d1**((j-1)*(coolTab_logTup-coolTab_logTlow) &
    /(jmax-1) + coolTab_logTlow)
!produce cooling rates for the given Tgas
coolTab(:,j) = coolingZ_rates(Tgas)
!store Tgas into the array
coolTab_T(j) = Tgas
!save 1/dT since it is known
if(j>1) inv_coolTab_T(j-1) = 1d0 / (Tgas-Tgasold)
Tgasold = Tgas
end do

end subroutine coolingZ_init_tabs

!*******************************
!this subroutine solves a non linear system
! with the equations stored in fcn function
! and a dummy jacobian jcn
subroutine nleq_wrap(x)
use krome_user_commons
integer,parameter::nmax=100 !problem size
integer,parameter::liwk=nmax+50 !size integer workspace
integer,parameter::lrwk=(nmax+13)*nmax+60 !real workspace
integer,parameter::luprt=6 !logical unit verbose output
integer::neq,iopt(50),ierr,niw,nrw,iwk(liwk),ptype,i
real*8::x(:),xscal(nmax),rtol,rwk(lrwk),idamp,mdamp,xi(size(x)),minx
real*8::store_invdvdz
neq = size(x)
niw = neq+50
nrw = (neq+13)*neq+60

ptype = 2 !initial problem type, 2=mildly non-linear
rtol = 1d-5 !realtive tolerance
xi(:) = x(:) !store initial guess
idamp = 1d-4 !initial damp (when ptype>=4, else default)
mdamp = 1d-8 !minimum damp (when ptype>=4, else default)
ierr = 0

!iterate until ierr==0 and non-negative solutions
do
if(ptype>50) then
print *,"ERROR in nleq1: can't find a solution after attempt",ptype
stop
end if

x(:) = xi(:) !restore initial guess

!if damping error or negative solutions
! prepares initial guess with the thin case
if(ptype>7.and.(ierr==3.or.ierr==0)) then
rtol = 1d-5
iwk(:) = 0
iopt(:) = 0
rwk(:) = 0d0
xscal(:) = 0d0
store_invdvdz = krome_invdvdz !store global variable
krome_invdvdz = 0d0 !this sets beta to 1
if(ierr.ne.0) then
print *,"ERROR in nleq for thin approx",ierr
stop
end if
krome_invdvdz = store_invdvdz !restore global variable
end if
xscal(:) = 0d0 !scaling factor
rtol = 1d-5 !relative tolerance
iwk(:) = 0 !default iwk
iwk(31) = int(1e8) !max iterations
iopt(:) = 0 !default iopt
iopt(31) = min(ptype,4) !problem type
rwk(:) = 0d0 !default rwk
!reduce damps if damping error
if(ptype>4.and.ierr==3) then
idamp = idamp * 1d-1 !reduce idamp
mdamp = mdamp * 1d-1 !reduce mdamp
end if
!if problem is extremely nonlinear use custom damps
if(ptype>4) then
rwk(21) = idamp !copy idamp to solver
rwk(22) = mdamp !copy mdamp to solver
end if

!check for errors
if(ierr.ne.0) then
!print *,"error",ierr
!problem with damping factor and/or problem type
if(ierr==3) then
ptype = ptype + 1 !change the problem type (non-linearity)
elseif(ierr==5) then
xi(:) = x(:)
else
!other type of error hence stop
print *,"ERROR in nleq1, ierr:",ierr
print *,"solutions found so far:"
do i=1,size(x)
print *,i,x(i)
end do
stop
end if
else
!if succesful search for negative results
minx = minval(x) !minimum value
!if minimum value is positive OK
if(minx.ge.0d0) then
exit
else
!if negative values are small set to zero
if(abs(minx)/maxval(x)<rtol) then
do i=1,neq
  x(i) = max(x(i),0d0)
end do
exit
else
!if large negative values increase non-linearity
ptype = ptype + 1
end if
end if
end if
end do
end subroutine nleq_wrap

!***************************
subroutine fcn(n,x,f,ierr)
implicit none
integer::n,ierr
real*8::x(n),f(n)

end subroutine fcn

!**********************************
!dummy jacobian for non linear equation solver
subroutine jcn()

end subroutine jcn

!************************************
function coolingCI(n,inTgas,k)
use krome_commons
use krome_photo
use krome_subs
implicit none
integer::i, hasnegative, nmax
real*8::coolingCI,n(:),inTgas,k(:)
real*8::A(3,3),Ain(3,3)
real*8::B(3),tmp(3)
real*8::coll_H2or
real*8::coll_Hj
real*8::coll_e
real*8::coll_H2pa
real*8::coll_H

!colliders should be >0
coll_H2or = max(n(idx_H2) * phys_orthoParaRatio / (phys_orthoParaRatio+1d0), 0d0)
coll_Hj = max(n(idx_Hj), 0d0)
coll_e = max(n(idx_e), 0d0)
coll_H2pa = max(n(idx_H2) / (phys_orthoParaRatio+1d0), 0d0)
coll_H = max(n(idx_H), 0d0)

!deafault cooling value
coolingCI = 0d0

if(n(idx_C)<1d-15) return

A(:,:) = 0d0

A(1,1) = 1d0
A(2,1) = + k(30) * coll_H2or &
    + k(33) * coll_Hj &
    + k(28) * coll_e &
    + k(31) * coll_H2pa &
    + k(40) * coll_H
A(3,1) = + k(41) * coll_Hj &
    + k(32) * coll_H2pa &
    + k(27) * coll_H2or &
    + k(29) * coll_e &
    + k(36) * coll_H
A(2,2) = - k(13) * coll_H2pa &
    - k(1) * coll_H &
    - k(7) * coll_e &
    - k(10) * coll_H2or &
    - k(4) * coll_Hj &
    - 7.900000d-08 &
    - k(38) * coll_H2pa &
    - k(39) * coll_H2or &
    - k(34) * coll_e &
    - k(37) * coll_Hj &
    - k(35) * coll_H
A(3,3) = - k(11) * coll_H2or &
    - k(14) * coll_H2pa &
    - k(2) * coll_H &
    - k(5) * coll_Hj &
    - k(8) * coll_e &
    - 2.100000d-14 &
    - k(9) * coll_e &
    - k(6) * coll_Hj &
    - k(15) * coll_H2pa &
    - k(3) * coll_H &
    - k(12) * coll_H2or &
    - 2.700000d-07
A(1,2) = 1d0
A(3,2) = + k(38) * coll_H2pa &
    + k(39) * coll_H2or &
    + k(34) * coll_e &
    + k(37) * coll_Hj &
    + k(35) * coll_H
A(1,3) = 1d0
A(2,3) = + k(9) * coll_e &
    + k(6) * coll_Hj &
    + k(15) * coll_H2pa &
    + k(3) * coll_H &
    + k(12) * coll_H2or &
    + 2.700000d-07

!build matrix B
B(:) = 0d0
B(1) = n(idx_C)

Ain(:,:) = A(:,:)

call mylin3(A(:,:), B(:))

!store population
pop_level_CI(:) = B(:)
!sanitize negative values
hasnegative = 0
do i=1,3
if(B(i)<0d0) then
if(abs(B(i)/n(idx_C))>1d-10) then
hasnegative = 1
else
B(i) = 1d-40
end if
end if
end do

!check if B has negative values
if(hasnegative>0)then
print *,"ERROR: minval(B)<0d0 in coolingCI"
print *,"ntot_CI =", n(idx_C)
print *,"Tgas =", inTgas
print *,"B(:) unrolled:"
do i=1,size(B)
print *, i, B(i)
end do
print *,"A(:,:) min/max:"
do i=1,size(B)
print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
end do

print *,"A(:,:)"
do i=1,size(B)
tmp(:) = Ain(i,:)
print '(I5,99E17.8)', i, tmp(:)
end do
stop
end if

coolingCI = B(2) * (7.900000d-08) * 2.400000d+01 &
    + B(3) * (2.700000d-07) * 3.900000d+01 &
    + B(3) * (2.100000d-14) * 6.300000d+01

end function coolingCI

!************************************
function coolingCII(n,inTgas,k)
use krome_commons
use krome_photo
use krome_subs
implicit none
integer::i, hasnegative, nmax
real*8::coolingCII,n(:),inTgas,k(:)
real*8::A(2,2),Ain(2,2)
real*8::B(2),tmp(2)
real*8::coll_e
real*8::coll_H

!colliders should be >0
coll_e = max(n(idx_e), 0d0)
coll_H = max(n(idx_H), 0d0)

!deafault cooling value
coolingCII = 0d0

if(n(idx_Cj)<1d-15) return

A(:,:) = 0d0

A(1,1) = 1d0
A(2,1) = + k(42) * coll_e &
    + k(43) * coll_H
A(2,2) = - k(25) * coll_e &
    - k(26) * coll_H &
    - 2.400000d-06
A(1,2) = 1d0

!build matrix B
B(:) = 0d0
B(1) = n(idx_Cj)

Ain(:,:) = A(:,:)

call mylin2(A(:,:), B(:))

!store population
pop_level_CII(:) = B(:)
!sanitize negative values
hasnegative = 0
do i=1,2
if(B(i)<0d0) then
if(abs(B(i)/n(idx_Cj))>1d-10) then
hasnegative = 1
else
B(i) = 1d-40
end if
end if
end do

!check if B has negative values
if(hasnegative>0)then
print *,"ERROR: minval(B)<0d0 in coolingCII"
print *,"ntot_CII =", n(idx_Cj)
print *,"Tgas =", inTgas
print *,"B(:) unrolled:"
do i=1,size(B)
print *, i, B(i)
end do
print *,"A(:,:) min/max:"
do i=1,size(B)
print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
end do

print *,"A(:,:)"
do i=1,size(B)
tmp(:) = Ain(i,:)
print '(I5,99E17.8)', i, tmp(:)
end do
stop
end if

coolingCII = B(2) * (2.400000d-06) * 9.120000d+01

end function coolingCII

!************************************
function coolingOI(n,inTgas,k)
use krome_commons
use krome_photo
use krome_subs
implicit none
integer::i, hasnegative, nmax
real*8::coolingOI,n(:),inTgas,k(:)
real*8::A(3,3),Ain(3,3)
real*8::B(3),tmp(3)
real*8::coll_e
real*8::coll_H
real*8::coll_Hj

!colliders should be >0
coll_e = max(n(idx_e), 0d0)
coll_H = max(n(idx_H), 0d0)
coll_Hj = max(n(idx_Hj), 0d0)

!deafault cooling value
coolingOI = 0d0

if(n(idx_O)<1d-15) return

A(:,:) = 0d0

A(1,1) = 1d0
A(2,1) = + k(45) * coll_e &
    + k(44) * coll_H &
    + k(49) * coll_Hj
A(3,1) = + k(50) * coll_e &
    + k(48) * coll_H
A(2,2) = - k(22) * coll_e &
    - k(19) * coll_Hj &
    - k(16) * coll_H &
    - 8.900000d-05 &
    - k(47) * coll_H &
    - k(46) * coll_e
A(3,3) = - k(23) * coll_e &
    - k(20) * coll_H &
    - 1.800000d-05 &
    - k(24) * coll_e &
    - k(21) * coll_H &
    - 1.300000d-10
A(1,2) = 1d0
A(3,2) = + k(47) * coll_H &
    + k(46) * coll_e
A(1,3) = 1d0
A(2,3) = + k(24) * coll_e &
    + k(21) * coll_H &
    + 1.300000d-10

!build matrix B
B(:) = 0d0
B(1) = n(idx_O)

Ain(:,:) = A(:,:)

call mylin3(A(:,:), B(:))

!store population
pop_level_OI(:) = B(:)
!sanitize negative values
hasnegative = 0
do i=1,3
if(B(i)<0d0) then
if(abs(B(i)/n(idx_O))>1d-10) then
hasnegative = 1
else
B(i) = 1d-40
end if
end if
end do

!check if B has negative values
if(hasnegative>0)then
print *,"ERROR: minval(B)<0d0 in coolingOI"
print *,"ntot_OI =", n(idx_O)
print *,"Tgas =", inTgas
print *,"B(:) unrolled:"
do i=1,size(B)
print *, i, B(i)
end do
print *,"A(:,:) min/max:"
do i=1,size(B)
print *, i, minval(Ain(i,:)), maxval(Ain(i,:))
end do

print *,"A(:,:)"
do i=1,size(B)
tmp(:) = Ain(i,:)
print '(I5,99E17.8)', i, tmp(:)
end do
stop
end if

coolingOI = B(2) * (8.900000d-05) * 2.300000d+02 &
    + B(3) * (1.800000d-05) * 3.300000d+02 &
    + B(3) * (1.300000d-10) * 1.000000d+02

end function coolingOI

!***********************
subroutine mylin2(a,b)
!solve Ax=B analytically for a 2-levels system
implicit none
integer,parameter::n=2
real*8::a(n,n),b(n),c(n),iab

!uncomment this: safer but slower function
!if(a(2,2)==a(2,1)) then
!   print *,"ERROR: a22=a21 in mylin2"
!   stop
!end if
iab = b(1)/(a(2,2)-a(2,1))
c(1) = a(2,2) * iab
c(2) = -a(2,1) * iab
b(:) = c(:)

end subroutine mylin2

!************************
subroutine mylin3(a,b)
!solve Ax=B analytically for a 3-levels system
implicit none
integer,parameter::n=3
real*8::iab,a(n,n),b(n),c(n)

!uncomment this: safer but slower function
!if(a(2,2)==a(2,3)) then
!   print *,"ERROR: a22=a23 in mylin3"
!   stop
!end if

!uncomment this: safer but slower
!if(a(2,1)*a(3,2)+a(2,2)*a(3,3)+a(2,3)*a(3,1) == &
    !     a(2,1)*a(3,3)+a(2,2)*a(3,1)+a(2,3)*a(3,2)) then
!   print *,"ERROR: division by zero in mylin3"
!   stop
!end if

iab = b(1) / (a(2,1)*(a(3,3)-a(3,2)) + a(2,2)*(a(3,1)-a(3,3)) &
    + a(2,3)*(a(3,2)-a(3,1)))
c(1) = (a(2,3)*a(3,2)-a(2,2)*a(3,3)) * iab
c(2) = -(a(2,3)*a(3,1)-a(2,1)*a(3,3)) * iab
c(3) = (a(3,1)*a(2,2)-a(2,1)*a(3,2)) * iab
b(:) = c(:)

end subroutine mylin3

!************************************
subroutine plot_cool(n)
!routine to plot cooling at runtime
real*8::n(:),Tgas,Tmin,Tmax
real*8::cool_atomic,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
real*8::cool_dH,cool_Z
integer::i,imax
imax = 1000
Tmin = log10(1d1)
Tmax = log10(1d8)
print *,"plotting cooling..."
open(33,file="KROME_cooling_plot.dat",status="replace")
do i=1,imax
Tgas = 1d1**(i*(Tmax-Tmin)/imax+Tmin)
cool_H2 = 0.d0
cool_H2GP = 0.d0
cool_HD = 0.d0
cool_atomic = 0.d0
cool_Z = 0.d0
cool_dH = 0.d0
cool_H2 = cooling_H2(n(:),Tgas)
cool_atomic = cooling_atomic(n(:),Tgas)
cool_Z = cooling_Z(n(:),Tgas)
cool_tot = cool_H2 + cool_atomic + cool_HD + cool_Z + cool_dH
cool_totGP = cool_H2GP + cool_atomic + cool_HD + cool_Z + cool_dH
write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, &
    cool_atomic, cool_HD, cool_H2GP, cool_Z, cool_dH
end do
close(33)
print *,"done!"

end subroutine plot_cool

!***********************************
!routine to dump cooling in unit nfile
subroutine dump_cool(n,Tgas,nfile)
use krome_commons
implicit none
real*8::Tgas,n(:),cools(ncools)
integer::nfile

cools(:) = get_cooling_array(n(:),Tgas)
write(nfile,'(99E14.5e3)') Tgas, sum(cools), cools(:)

end subroutine dump_cool

end module KROME_cooling

!############### MODULE ##############
module KROME_heating
contains

! *************************************************************
!  This file has been generated with:
!  KROME 14.08.dev on 2018-10-24 12:49:10
!  Changeset b21f657
!  see http://kromepackage.org
!
!  Written and developed by Tommaso Grassi and Stefano Bovino
!
!  Contributors:
!  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
!  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
!  E.Tognelli
!  KROME is provided "as it is", without any warranty.
! *************************************************************

!************************
function heating(n,inTgas,k,nH2dust)
implicit none
real*8::n(:), Tgas, inTgas, k(:), nH2dust
real*8::heating

Tgas = inTgas
heating = sum(get_heating_array(n(:),Tgas,k(:), nH2dust))

end function heating

!*******************************
function get_heating_array(n, Tgas, k, nH2dust)
use krome_commons
implicit none
real*8::n(:), Tgas, k(:), nH2dust
real*8::get_heating_array(nheats),heats(nheats)
real*8::smooth,f1,f2
!returns heating in erg/cm3/s

heats(:) = 0.d0

f2 = 1.
!this parameter controls the smoothness of the
! merge between the two cooling / heating functions
smooth = 1.d-3

!smoothing functions | f1+f2=1
!f1 = (tanh(smooth*(Tgas-1d4))+1.d0)*0.5d0
f2 = (tanh(smooth*(-Tgas+1d4))+1.d0)*0.5d0

!heating is already included in the GH cooling function (that thus can be negative), and f1 is not needed

heats(idx_heat_chem) = heatingChem(n(:), Tgas, k(:), nH2dust)

heats(idx_heat_photo) = photo_heating(n(:),f2)

heats(idx_heat_photoAv) = heat_photoAv(n(:),Tgas,k(:))

heats(idx_heat_CR) = heat_CR(n(:),Tgas,k(:))

heats(idx_heat_dust) = heat_photoDust(n(:),Tgas)

heats(idx_heat_custom) = heat_custom(n(:),Tgas)

get_heating_array(:) = heats(:)

end function get_heating_array

!*************************
function heat_custom(n,Tgas)
use krome_commons
use krome_subs
use krome_constants
implicit none
real*8::n(:),Tgas,heat_custom

heat_custom = 0d0

end function heat_custom

!***************************
function heat_photoDust(n,Tgas)
!photoelectric effect from dust in erg/s/cm3
!see Bakes&Tielens 1994 with a slight modification of Wolfire 2003
!on the amount of absorbed ultraviolet energy.
!This is for the local interstellar Habing flux and
!without considering the recombination (which at this
!radiation flux is indeed negligible)
use krome_commons
use krome_subs
use krome_getphys
implicit none
real*8::heat_photoDust,n(:),Tgas,ntot,eps
real*8::Ghab,z,psi

ntot = get_Hnuclei(n(:))
Ghab = 1.69d0 !habing flux, 1.69 is Draine78
Ghab  = user_G0
Ghab  = Ghab * exp(-2.5*user_av)
if(n(idx_e)>0d0) then
psi = Ghab * sqrt(Tgas) / n(idx_e)
else
psi = 1d99
end if
eps = 4.9d-2 / (1d0 + 4d-3 * psi**.73) + &
    3.7d-2 * (Tgas * 1d-4)**.7 / (1d0 + 2d-4 * psi)
z = 1d1**get_metallicityC(n(:)) !metallicty wrt solar
heat_photoDust = 1.3d-24*eps*Ghab*ntot*z

end function heat_photoDust

!******************************
function heat_photoAv(n,Tgas,k)
!heating from photoreactions using rate approximation (erg/s/cm3)
use krome_commons
use krome_user_commons
use krome_subs
use krome_getphys
implicit none
real*8::heat_photoAv,n(:),Tgas,k(:)
real*8::ncrn,ncrd1,ncrd2,yH,yH2,ncr,h2heatfac,dd,Rdiss

dd = get_Hnuclei(n(:))
ncrn  = 1.0d6*(Tgas**(-0.5d0))
ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**2)
ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

yH = n(idx_H)/dd   !dimensionless
yH2= n(idx_H2)/dd  !dimensionless

ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

Rdiss = k(231)

!photodissociation H2 heating
heat_photoAv = 6.4d-13*Rdiss*n(idx_H2)

!UV photo-pumping H2
heat_photoAv = heat_photoAv + 2.7d-11*Rdiss*h2heatfac*n(idx_H2)

end function heat_photoAv

!***************************
function heat_CR(n,Tgas,k)
!heating from cosmic rays erg/s/cm3
use krome_commons
implicit none
real*8::heat_CR,n(:),Tgas,Hfact,k(:)
real*8::logH2,QH2,QH,QHe,ev2erg

ev2erg = 1.60217662d-12
Hfact = 2d1*ev2erg !erg

!precompute log10(H2)
logH2 = log10(max(n(idx_H2),1d-40))

!init heating
heat_CR = 0d0

!heating per H ionization (eV)
QH = 4.3d0 * ev2erg

!heating per He ionization, same as H following Glassgold+2012
QHe = QH

!automatically generated fit function
! for H2 ionization, units: eV
! data from arXiv:1502.03380
if(n(idx_H2)>1d10) then
QH2 = 18.0217195233
else if((n(idx_H2).ge.1.03660174949e-05).and.(n(idx_H2)<1d5)) then
QH2 = 1.45108491351*logH2 + 7.23277032061
else if(n(idx_H2)<1.03660174949e-05) then
QH2 = 0d0
else
QH2 = 9.14621986426 &
    - 0.443120145068*logH2 &
    + 0.60127161756*logH2**2 &
    - 0.0710284101055*logH2**3 &
    + 0.00242079494592*logH2**4
end if

!convert eV to erg
QH2 = QH2 * ev2erg

!H -> H+ + E
heat_CR = heat_CR + k(252) * n(idx_H) * QH

!HE -> HE+ + E
heat_CR = heat_CR + k(253) * n(idx_HE) * QHe

!O -> O+ + E
heat_CR = heat_CR + k(254) * n(idx_O) * Hfact

!CO -> C + O
heat_CR = heat_CR + k(255) * n(idx_CO) * Hfact

!CO -> CO+ + E
heat_CR = heat_CR + k(256) * n(idx_CO) * Hfact

!C2 -> C + C
heat_CR = heat_CR + k(257) * n(idx_C2) * Hfact

!H2 -> H + H
heat_CR = heat_CR + k(258) * n(idx_H2) * Hfact

!H2 -> H+ + H-
heat_CR = heat_CR + k(259) * n(idx_H2) * Hfact

!H2 -> H2+ + E
heat_CR = heat_CR + k(260) * n(idx_H2) * QH2

!C -> C+ + E
heat_CR = heat_CR + k(261) * n(idx_C) * Hfact

!CH -> C + H
heat_CR = heat_CR + k(262) * n(idx_CH) * Hfact

!O2 -> O + O
heat_CR = heat_CR + k(263) * n(idx_O2) * Hfact

!O2 -> O2+ + E
heat_CR = heat_CR + k(264) * n(idx_O2) * Hfact

!OH -> O + H
heat_CR = heat_CR + k(265) * n(idx_OH) * Hfact

!CH2 -> CH2+ + E
heat_CR = heat_CR + k(266) * n(idx_CH2) * Hfact

!H2O -> OH + H
heat_CR = heat_CR + k(267) * n(idx_H2O) * Hfact

!HCO -> CO + H
heat_CR = heat_CR + k(268) * n(idx_HCO) * Hfact

!HCO -> HCO+ + E
heat_CR = heat_CR + k(269) * n(idx_HCO) * Hfact

!H2 -> H + H+ + E
heat_CR = heat_CR + k(270) * n(idx_H2) * Hfact

end function heat_CR

!**************************
function photo_heating(n,f2)
!photo heating in erg/cm3/s using bin-based
! approach. Terms are computed in the
! krome_photo module
use krome_commons
use krome_constants
implicit none
real*8::photo_heating,n(:),f2

photo_heating = 0.d0
photo_heating = photoBinHeats(1) * n(idx_H) &
    + photoBinHeats(2) * n(idx_HE) &
    + photoBinHeats(3) * n(idx_HEj) &
    + f2 *photoBinHeats(4) * n(idx_O) &
    + f2 *photoBinHeats(5) * n(idx_C) &
    + photoBinHeats(6) * n(idx_H2) &
    + photoBinHeats(7) * n(idx_Hk) &
    + photoBinHeats(8) * n(idx_CH) &
    + photoBinHeats(9) * n(idx_CH) &
    + photoBinHeats(10) * n(idx_C2) &
    + photoBinHeats(11) * n(idx_OH) &
    + photoBinHeats(12) * n(idx_OH) &
    + photoBinHeats(13) * n(idx_H2O) &
    + photoBinHeats(14) * n(idx_H2O) &
    + photoBinHeats(15) * n(idx_O2) &
    + photoBinHeats(16) * n(idx_O2) &
    + photoBinHeats(17) * n(idx_H2)

end function photo_heating

!H2 FORMATION HEATING and other exo/endothermic
! processes (including H2 on dust) in erg/cm3/s
!krome builds the heating/cooling term according
! to the chemical network employed
!*******************************
function heatingChem(n, Tgas, k, nH2dust)
use krome_constants
use krome_commons
use krome_dust
use krome_subs
use krome_getphys
implicit none
real*8::heatingChem, n(:), Tgas,k(:),nH2dust
real*8::h2heatfac,HChem,yH,yH2
real*8::ncr,ncrn,ncrd1,ncrd2,dd,n2H,small,nmax
dd = get_Hnuclei(n(:))

!replace small according to the desired enviroment
! and remove nmax if needed
nmax = maxval(n(1:nmols))
small = 1d-40/(nmax*nmax*nmax)

heatingChem = 0.d0

ncrn  = 1.0d6*(Tgas**(-0.5d0))
ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**2)
ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

yH = n(idx_H)/dd   !dimensionless
yH2= n(idx_H2)/dd  !dimensionless

ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

HChem = 0.d0 !inits chemical heating
n2H = n(idx_H) * n(idx_H)

!H- + H -> H2 + E (heating)
HChem = HChem + k(17) * (3.53d0*h2heatfac*n(idx_Hk)*n(idx_H))
!H2+ + H -> H2 + H+ (heating)
HChem = HChem + k(20) * (1.83d0*h2heatfac*n(idx_H2j)*n(idx_H))
!H2 + E -> H + H + E (cooling)
HChem = HChem + k(23) * (-4.48d0*n(idx_H2)*n(idx_E))
!H2 + H -> H + H + H (cooling)
HChem = HChem + k(24) * (-4.48d0*n(idx_H2)*n(idx_H))
!H2 + H2 -> H2 + H + H (cooling)
HChem = HChem + k(33) * (-4.48d0*n(idx_H2)*n(idx_H2))
!H + H + H -> H2 + H (heating)
HChem = HChem + k(35) * (4.48d0*h2heatfac*n(idx_H)*n(idx_H)*n(idx_H))
!H2 + H + H -> H2 + H2 (heating)
HChem = HChem + k(36) * (4.48d0*h2heatfac*n(idx_H2)*n(idx_H)*n(idx_H))

HChem = HChem + nH2dust * (4.2d0*h2heatfac + 0.2d0)

heatingChem = HChem * eV_to_erg  !erg/cm3/s

end function heatingChem

end module KROME_heating

!############### MODULE ##############
module krome_ode
contains

! *************************************************************
!  This file has been generated with:
!  KROME 14.08.dev on 2018-10-24 12:49:10
!  Changeset b21f657
!  see http://kromepackage.org
!
!  Written and developed by Tommaso Grassi and Stefano Bovino
!
!  Contributors:
!  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
!  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
!  E.Tognelli
!  KROME is provided "as it is", without any warranty.
! *************************************************************

subroutine fex(neq,tt,nin,dn)
use krome_commons
use krome_constants
use krome_subs
use krome_cooling
use krome_heating
use krome_tabs
use krome_photo
use krome_gadiab
use krome_getphys
use krome_phfuncs
use krome_fit
implicit none
integer::neq,idust
real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq)
real*8::dnChem_H2O,dnChem_CO
real*8::rr
integer::i,r1,r2,r3,p1,p2,p3

n(:) = nin(:)

nH2dust = 0.d0
n(idx_CR) = 1.d0
n(idx_g)  = 1.d0
n(idx_dummy) = 1.d0

dn(:) = 0.d0 !initialize differentials
n(idx_Tgas) = max(n(idx_tgas),2.73d0)
n(idx_Tgas) = min(n(idx_tgas),1d9)
Tgas = n(idx_Tgas) !get temperature

k(:) = coe_tab(n(:)) !compute coefficients

dnChem_H2O = 0d0  &
    +k(73)*n(idx_H2)*n(idx_OH)  &
    +k(78)*n(idx_OH)*n(idx_OH)  &
    -k(79)*n(idx_H2O)*n(idx_H)  &
    -k(111)*n(idx_H2O)*n(idx_H3j)  &
    -k(112)*n(idx_H2O)*n(idx_H3j)  &
    -k(113)*n(idx_H2O)*n(idx_Cj)  &
    -k(114)*n(idx_H2O)*n(idx_Cj)  &
    -k(115)*n(idx_H2O)*n(idx_Cj)  &
    -k(116)*n(idx_H2O)*n(idx_Cj)  &
    -k(128)*n(idx_HCOj)*n(idx_H2O)  &
    -k(129)*n(idx_HCOj)*n(idx_H2O)  &
    -k(145)*n(idx_H2O)*n(idx_Hj)  &
    -k(146)*n(idx_H2O)*n(idx_Hj)  &
    -k(147)*n(idx_H2O)*n(idx_HEj)  &
    -k(148)*n(idx_H2O)*n(idx_HEj)  &
    -k(149)*n(idx_H2O)*n(idx_HEj)  &
    -k(150)*n(idx_H2O)*n(idx_HEj)  &
    -k(151)*n(idx_H2O)*n(idx_HEj)  &
    -k(152)*n(idx_H2O)*n(idx_HEj)  &
    +k(177)*n(idx_H3Oj)*n(idx_E)  &
    +k(186)*n(idx_Hk)*n(idx_OH)  &
    +k(191)*n(idx_Ok)*n(idx_H2)  &
    +k(207)*n(idx_OH)*n(idx_H)  &
    -k(225)*n(idx_H2O)  &
    -k(226)*n(idx_H2O)  &
    +k(248)*n(idx_H3Oj)  &
    -k(267)*n(idx_H2O)  &
    +k(280)*n(idx_OH)*n(idx_H)
dnChem_CO = 0d0  &
    -k(54)*n(idx_HOCj)*n(idx_CO)  &
    +k(54)*n(idx_HOCj)*n(idx_CO)  &
    -k(55)*n(idx_HOCj)*n(idx_CO)  &
    +k(55)*n(idx_HOCj)*n(idx_CO)  &
    +k(60)*n(idx_CH)*n(idx_O)  &
    +k(64)*n(idx_CH2)*n(idx_O)  &
    +k(65)*n(idx_CH2)*n(idx_O)  &
    +k(68)*n(idx_C2)*n(idx_O)  &
    +k(69)*n(idx_C2)*n(idx_O)  &
    +k(74)*n(idx_C)*n(idx_OH)  &
    +k(75)*n(idx_C)*n(idx_OH)  &
    +k(82)*n(idx_O2)*n(idx_C)  &
    +k(83)*n(idx_O2)*n(idx_C)  &
    -k(84)*n(idx_CO)*n(idx_H)  &
    +k(119)*n(idx_O2)*n(idx_Cj)  &
    -k(123)*n(idx_CO)*n(idx_H3j)  &
    -k(124)*n(idx_CO)*n(idx_H3j)  &
    -k(125)*n(idx_CO)*n(idx_H3j)  &
    -k(126)*n(idx_CO)*n(idx_H3j)  &
    +k(127)*n(idx_HCOj)*n(idx_C)  &
    +k(128)*n(idx_HCOj)*n(idx_H2O)  &
    +k(129)*n(idx_HCOj)*n(idx_H2O)  &
    -k(156)*n(idx_CO)*n(idx_HEj)  &
    -k(157)*n(idx_CO)*n(idx_HEj)  &
    +k(158)*n(idx_COj)*n(idx_H)  &
    +k(181)*n(idx_HCOj)*n(idx_E)  &
    +k(183)*n(idx_HOCj)*n(idx_E)  &
    +k(189)*n(idx_Ck)*n(idx_O)  &
    +k(192)*n(idx_Ok)*n(idx_C)  &
    +k(199)*n(idx_C)*n(idx_O)  &
    -k(230)*n(idx_CO)  &
    -k(255)*n(idx_CO)  &
    -k(256)*n(idx_CO)  &
    +k(268)*n(idx_HCO)  &
    +k(273)*n(idx_C)*n(idx_O)  &
    +k(274)*n(idx_C)*n(idx_O)

ntot = sum(n(1:nmols))
nH2dust = get_mu(n(:)) * n(idx_H) * 1d1**fit_anytab2D(dust_tab_ngas(:),					dust_tab_Tgas(:), &
    dust_tab_H2(:,:), dust_mult_ngas, dust_mult_Tgas, &
    log10(ntot), log10(Tgas)) * ntot

!E
dn(idx_E) = &
    -k(1)*n(idx_H)*n(idx_E) &
    +2.d0*k(1)*n(idx_H)*n(idx_E) &
    -k(2)*n(idx_Hj)*n(idx_E) &
    -k(3)*n(idx_Hj)*n(idx_E) &
    -k(4)*n(idx_HE)*n(idx_E) &
    +2.d0*k(4)*n(idx_HE)*n(idx_E) &
    -k(5)*n(idx_HEj)*n(idx_E) &
    -k(6)*n(idx_HEj)*n(idx_E) &
    -k(7)*n(idx_HEj)*n(idx_E) &
    +2.d0*k(7)*n(idx_HEj)*n(idx_E) &
    -k(15)*n(idx_HEjj)*n(idx_E) &
    -k(16)*n(idx_H)*n(idx_E) &
    +k(17)*n(idx_Hk)*n(idx_H) &
    -k(23)*n(idx_H2)*n(idx_E) &
    +k(23)*n(idx_H2)*n(idx_E) &
    -k(25)*n(idx_Hk)*n(idx_E) &
    +2.d0*k(25)*n(idx_Hk)*n(idx_E) &
    +k(26)*n(idx_Hk)*n(idx_H) &
    +k(27)*n(idx_Hk)*n(idx_H) &
    +k(29)*n(idx_Hk)*n(idx_Hj) &
    -k(30)*n(idx_H2j)*n(idx_E) &
    -k(31)*n(idx_H2j)*n(idx_E) &
    -k(37)*n(idx_Cj)*n(idx_E) &
    -k(38)*n(idx_Cj)*n(idx_E) &
    -k(39)*n(idx_Cj)*n(idx_E) &
    -k(40)*n(idx_Oj)*n(idx_E) &
    -k(41)*n(idx_Oj)*n(idx_E) &
    -k(42)*n(idx_C)*n(idx_E) &
    +2.d0*k(42)*n(idx_C)*n(idx_E) &
    -k(43)*n(idx_O)*n(idx_E) &
    +2.d0*k(43)*n(idx_O)*n(idx_E) &
    +k(61)*n(idx_CH)*n(idx_O) &
    -k(162)*n(idx_H3j)*n(idx_E) &
    -k(163)*n(idx_H3j)*n(idx_E) &
    -k(164)*n(idx_CHj)*n(idx_E) &
    -k(165)*n(idx_CH2j)*n(idx_E) &
    -k(166)*n(idx_CH2j)*n(idx_E) &
    -k(167)*n(idx_CH2j)*n(idx_E) &
    -k(168)*n(idx_CH3j)*n(idx_E) &
    -k(169)*n(idx_CH3j)*n(idx_E) &
    -k(170)*n(idx_CH3j)*n(idx_E) &
    -k(171)*n(idx_OHj)*n(idx_E) &
    -k(172)*n(idx_H2Oj)*n(idx_E) &
    -k(173)*n(idx_H2Oj)*n(idx_E) &
    -k(174)*n(idx_H2Oj)*n(idx_E) &
    -k(175)*n(idx_H3Oj)*n(idx_E) &
    -k(176)*n(idx_H3Oj)*n(idx_E) &
    -k(177)*n(idx_H3Oj)*n(idx_E) &
    -k(178)*n(idx_H3Oj)*n(idx_E) &
    -k(179)*n(idx_O2j)*n(idx_E) &
    -k(180)*n(idx_COj)*n(idx_E) &
    -k(181)*n(idx_HCOj)*n(idx_E) &
    -k(182)*n(idx_HCOj)*n(idx_E) &
    -k(183)*n(idx_HOCj)*n(idx_E) &
    +k(184)*n(idx_Hk)*n(idx_C) &
    +k(185)*n(idx_Hk)*n(idx_O) &
    +k(186)*n(idx_Hk)*n(idx_OH) &
    +k(187)*n(idx_Ck)*n(idx_H) &
    +k(188)*n(idx_Ck)*n(idx_H2) &
    +k(189)*n(idx_Ck)*n(idx_O) &
    +k(190)*n(idx_Ok)*n(idx_H) &
    +k(191)*n(idx_Ok)*n(idx_H2) &
    +k(192)*n(idx_Ok)*n(idx_C) &
    -k(195)*n(idx_C)*n(idx_E) &
    -k(204)*n(idx_O)*n(idx_E) &
    +k(213)*n(idx_H) &
    +k(214)*n(idx_HE) &
    +k(215)*n(idx_HEj) &
    +k(216)*n(idx_O) &
    +k(217)*n(idx_C) &
    +k(218)*n(idx_H2) &
    +k(219)*n(idx_Hk) &
    +k(221)*n(idx_CH) &
    +k(224)*n(idx_OH) &
    +k(226)*n(idx_H2O) &
    +k(227)*n(idx_O2) &
    +k(229)*n(idx_H2) &
    +k(235)*n(idx_Ck) &
    +k(238)*n(idx_CH2) &
    +k(242)*n(idx_Ok) &
    +k(252)*n(idx_H) &
    +k(253)*n(idx_HE) &
    +k(254)*n(idx_O) &
    +k(256)*n(idx_CO) &
    +k(260)*n(idx_H2) &
    +k(261)*n(idx_C) &
    +k(264)*n(idx_O2) &
    +k(266)*n(idx_CH2) &
    +k(269)*n(idx_HCO) &
    +k(270)*n(idx_H2)
!H-
dn(idx_Hk) = &
    +k(16)*n(idx_H)*n(idx_E) &
    -k(17)*n(idx_Hk)*n(idx_H) &
    -k(25)*n(idx_Hk)*n(idx_E) &
    -k(26)*n(idx_Hk)*n(idx_H) &
    -k(27)*n(idx_Hk)*n(idx_H) &
    -k(28)*n(idx_Hk)*n(idx_Hj) &
    -k(29)*n(idx_Hk)*n(idx_Hj) &
    -k(32)*n(idx_H2j)*n(idx_Hk) &
    -k(161)*n(idx_HEj)*n(idx_Hk) &
    -k(184)*n(idx_Hk)*n(idx_C) &
    -k(185)*n(idx_Hk)*n(idx_O) &
    -k(186)*n(idx_Hk)*n(idx_OH) &
    -k(219)*n(idx_Hk) &
    +k(259)*n(idx_H2)
!C-
dn(idx_Ck) = &
    -k(159)*n(idx_Ck)*n(idx_Hj) &
    -k(187)*n(idx_Ck)*n(idx_H) &
    -k(188)*n(idx_Ck)*n(idx_H2) &
    -k(189)*n(idx_Ck)*n(idx_O) &
    +k(195)*n(idx_C)*n(idx_E) &
    -k(235)*n(idx_Ck)
!O-
dn(idx_Ok) = &
    -k(160)*n(idx_Ok)*n(idx_Hj) &
    -k(190)*n(idx_Ok)*n(idx_H) &
    -k(191)*n(idx_Ok)*n(idx_H2) &
    -k(192)*n(idx_Ok)*n(idx_C) &
    +k(204)*n(idx_O)*n(idx_E) &
    -k(242)*n(idx_Ok)
!H
dn(idx_H) = &
    -k(1)*n(idx_H)*n(idx_E) &
    +k(2)*n(idx_Hj)*n(idx_E) &
    +k(3)*n(idx_Hj)*n(idx_E) &
    -k(8)*n(idx_HEj)*n(idx_H) &
    +k(9)*n(idx_HE)*n(idx_Hj) &
    +k(10)*n(idx_HE)*n(idx_Hj) &
    +2.d0*k(11)*n(idx_H2)*n(idx_HE) &
    +k(13)*n(idx_H2)*n(idx_HEj) &
    +2.d0*k(14)*n(idx_H2)*n(idx_HEj) &
    -k(16)*n(idx_H)*n(idx_E) &
    -k(17)*n(idx_Hk)*n(idx_H) &
    -k(18)*n(idx_H)*n(idx_Hj) &
    -k(19)*n(idx_H)*n(idx_Hj) &
    -k(20)*n(idx_H2j)*n(idx_H) &
    +k(21)*n(idx_H2)*n(idx_Hj) &
    +k(22)*n(idx_H2)*n(idx_Hj) &
    +2.d0*k(23)*n(idx_H2)*n(idx_E) &
    -k(24)*n(idx_H2)*n(idx_H) &
    +3.d0*k(24)*n(idx_H2)*n(idx_H) &
    +k(25)*n(idx_Hk)*n(idx_E) &
    -k(26)*n(idx_Hk)*n(idx_H) &
    +2.d0*k(26)*n(idx_Hk)*n(idx_H) &
    -k(27)*n(idx_Hk)*n(idx_H) &
    +2.d0*k(27)*n(idx_Hk)*n(idx_H) &
    +2.d0*k(28)*n(idx_Hk)*n(idx_Hj) &
    +2.d0*k(30)*n(idx_H2j)*n(idx_E) &
    +2.d0*k(31)*n(idx_H2j)*n(idx_E) &
    +k(32)*n(idx_H2j)*n(idx_Hk) &
    +2.d0*k(33)*n(idx_H2)*n(idx_H2) &
    -2.d0*k(34)*n(idx_H)*n(idx_H)*n(idx_HE) &
    -3.d0*k(35)*n(idx_H)*n(idx_H)*n(idx_H) &
    +k(35)*n(idx_H)*n(idx_H)*n(idx_H) &
    -2.d0*k(36)*n(idx_H2)*n(idx_H)*n(idx_H) &
    -k(44)*n(idx_Oj)*n(idx_H) &
    +k(45)*n(idx_O)*n(idx_Hj) &
    +k(47)*n(idx_C)*n(idx_Hj) &
    -k(48)*n(idx_Cj)*n(idx_H) &
    -k(52)*n(idx_OH)*n(idx_H) &
    +2.d0*k(52)*n(idx_OH)*n(idx_H) &
    +k(56)*n(idx_C)*n(idx_H2) &
    -k(57)*n(idx_CH)*n(idx_H) &
    +k(58)*n(idx_CH)*n(idx_H2) &
    +k(59)*n(idx_CH)*n(idx_C) &
    +k(60)*n(idx_CH)*n(idx_O) &
    -k(63)*n(idx_CH2)*n(idx_H) &
    +2.d0*k(64)*n(idx_CH2)*n(idx_O) &
    +k(66)*n(idx_CH2)*n(idx_O) &
    +k(70)*n(idx_O)*n(idx_H2) &
    -k(71)*n(idx_OH)*n(idx_H) &
    -k(72)*n(idx_OH)*n(idx_H) &
    +k(73)*n(idx_H2)*n(idx_OH) &
    +k(74)*n(idx_C)*n(idx_OH) &
    +k(75)*n(idx_C)*n(idx_OH) &
    +k(76)*n(idx_O)*n(idx_OH) &
    +k(77)*n(idx_O)*n(idx_OH) &
    -k(79)*n(idx_H2O)*n(idx_H) &
    -k(80)*n(idx_O2)*n(idx_H) &
    -k(84)*n(idx_CO)*n(idx_H) &
    +k(85)*n(idx_H2j)*n(idx_H2) &
    -k(86)*n(idx_H3j)*n(idx_H) &
    +k(87)*n(idx_C)*n(idx_H2j) &
    +k(89)*n(idx_C)*n(idx_H3j) &
    +k(90)*n(idx_Cj)*n(idx_H2) &
    -k(91)*n(idx_CHj)*n(idx_H) &
    +k(92)*n(idx_CHj)*n(idx_H2) &
    +k(93)*n(idx_CHj)*n(idx_O) &
    -k(94)*n(idx_CH2j)*n(idx_H) &
    +k(95)*n(idx_CH2j)*n(idx_H2) &
    +k(96)*n(idx_CH2j)*n(idx_O) &
    -k(97)*n(idx_CH3j)*n(idx_H) &
    +k(101)*n(idx_Oj)*n(idx_H2) &
    +k(102)*n(idx_O)*n(idx_H2j) &
    +k(104)*n(idx_O)*n(idx_H3j) &
    +k(107)*n(idx_OH)*n(idx_Cj) &
    +k(108)*n(idx_OH)*n(idx_Cj) &
    +k(109)*n(idx_OHj)*n(idx_H2) &
    +k(110)*n(idx_H2Oj)*n(idx_H2) &
    +k(113)*n(idx_H2O)*n(idx_Cj) &
    +k(114)*n(idx_H2O)*n(idx_Cj) &
    +k(115)*n(idx_H2O)*n(idx_Cj) &
    +k(130)*n(idx_CH)*n(idx_Hj) &
    +k(131)*n(idx_CH)*n(idx_Hj) &
    +k(134)*n(idx_CH2)*n(idx_Hj) &
    +k(135)*n(idx_CH2)*n(idx_Hj) &
    +k(138)*n(idx_CH2)*n(idx_HEj) &
    +k(139)*n(idx_CH2)*n(idx_HEj) &
    +k(141)*n(idx_OH)*n(idx_Hj) &
    +k(142)*n(idx_OH)*n(idx_Hj) &
    +k(143)*n(idx_OH)*n(idx_HEj) &
    +k(144)*n(idx_OH)*n(idx_HEj) &
    +k(145)*n(idx_H2O)*n(idx_Hj) &
    +k(146)*n(idx_H2O)*n(idx_Hj) &
    +k(149)*n(idx_H2O)*n(idx_HEj) &
    +k(150)*n(idx_H2O)*n(idx_HEj) &
    +k(153)*n(idx_O2)*n(idx_Hj) &
    -k(158)*n(idx_COj)*n(idx_H) &
    +k(159)*n(idx_Ck)*n(idx_Hj) &
    +k(160)*n(idx_Ok)*n(idx_Hj) &
    +k(161)*n(idx_HEj)*n(idx_Hk) &
    +k(162)*n(idx_H3j)*n(idx_E) &
    +3.d0*k(163)*n(idx_H3j)*n(idx_E) &
    +k(164)*n(idx_CHj)*n(idx_E) &
    +k(165)*n(idx_CH2j)*n(idx_E) &
    +2.d0*k(167)*n(idx_CH2j)*n(idx_E) &
    +k(168)*n(idx_CH3j)*n(idx_E) &
    +2.d0*k(170)*n(idx_CH3j)*n(idx_E) &
    +k(171)*n(idx_OHj)*n(idx_E) &
    +k(173)*n(idx_H2Oj)*n(idx_E) &
    +2.d0*k(174)*n(idx_H2Oj)*n(idx_E) &
    +2.d0*k(175)*n(idx_H3Oj)*n(idx_E) &
    +k(176)*n(idx_H3Oj)*n(idx_E) &
    +k(177)*n(idx_H3Oj)*n(idx_E) &
    +k(181)*n(idx_HCOj)*n(idx_E) &
    +k(183)*n(idx_HOCj)*n(idx_E) &
    -k(187)*n(idx_Ck)*n(idx_H) &
    -k(190)*n(idx_Ok)*n(idx_H) &
    +2.d0*k(193)*n(idx_H2)*n(idx_Hj) &
    -k(196)*n(idx_C)*n(idx_H) &
    -k(200)*n(idx_Cj)*n(idx_H) &
    -k(205)*n(idx_O)*n(idx_H) &
    -k(207)*n(idx_OH)*n(idx_H) &
    -k(213)*n(idx_H) &
    +k(219)*n(idx_Hk) &
    +k(220)*n(idx_CH) &
    +k(223)*n(idx_OH) &
    +k(225)*n(idx_H2O) &
    +k(229)*n(idx_H2) &
    +2.d0*k(231)*n(idx_H2) &
    +k(232)*n(idx_H2j) &
    +k(234)*n(idx_H3j) &
    +k(237)*n(idx_CH2) &
    +k(239)*n(idx_CH2j) &
    +k(240)*n(idx_CH3j) &
    +k(247)*n(idx_H2Oj) &
    +k(250)*n(idx_H3Oj) &
    -k(252)*n(idx_H) &
    +2.d0*k(258)*n(idx_H2) &
    +k(262)*n(idx_CH) &
    +k(265)*n(idx_OH) &
    +k(267)*n(idx_H2O) &
    +k(268)*n(idx_HCO) &
    +k(270)*n(idx_H2) &
    -k(279)*n(idx_H)*n(idx_O) &
    -k(280)*n(idx_OH)*n(idx_H) - 2d0*nH2dust
!HE
dn(idx_HE) = &
    -k(4)*n(idx_HE)*n(idx_E) &
    +k(5)*n(idx_HEj)*n(idx_E) &
    +k(6)*n(idx_HEj)*n(idx_E) &
    +k(8)*n(idx_HEj)*n(idx_H) &
    -k(9)*n(idx_HE)*n(idx_Hj) &
    -k(10)*n(idx_HE)*n(idx_Hj) &
    -k(11)*n(idx_H2)*n(idx_HE) &
    +k(11)*n(idx_H2)*n(idx_HE) &
    +k(12)*n(idx_H2)*n(idx_HEj) &
    +k(13)*n(idx_H2)*n(idx_HEj) &
    -k(34)*n(idx_H)*n(idx_H)*n(idx_HE) &
    +k(34)*n(idx_H)*n(idx_H)*n(idx_HE) &
    +k(46)*n(idx_O)*n(idx_HEj) &
    +k(49)*n(idx_C)*n(idx_HEj) &
    +k(50)*n(idx_C)*n(idx_HEj) &
    +k(51)*n(idx_C)*n(idx_HEj) &
    +k(136)*n(idx_CH2)*n(idx_HEj) &
    +k(137)*n(idx_CH2)*n(idx_HEj) &
    +k(138)*n(idx_CH2)*n(idx_HEj) &
    +k(139)*n(idx_CH2)*n(idx_HEj) &
    +k(140)*n(idx_C2)*n(idx_HEj) &
    +k(143)*n(idx_OH)*n(idx_HEj) &
    +k(144)*n(idx_OH)*n(idx_HEj) &
    +k(147)*n(idx_H2O)*n(idx_HEj) &
    +k(148)*n(idx_H2O)*n(idx_HEj) &
    +k(149)*n(idx_H2O)*n(idx_HEj) &
    +k(150)*n(idx_H2O)*n(idx_HEj) &
    +k(151)*n(idx_H2O)*n(idx_HEj) &
    +k(152)*n(idx_H2O)*n(idx_HEj) &
    +k(154)*n(idx_O2)*n(idx_HEj) &
    +k(155)*n(idx_O2)*n(idx_HEj) &
    +k(156)*n(idx_CO)*n(idx_HEj) &
    +k(157)*n(idx_CO)*n(idx_HEj) &
    +k(161)*n(idx_HEj)*n(idx_Hk) &
    -k(214)*n(idx_HE) &
    -k(253)*n(idx_HE)
!H2
dn(idx_H2) = &
    -k(11)*n(idx_H2)*n(idx_HE) &
    -k(12)*n(idx_H2)*n(idx_HEj) &
    -k(13)*n(idx_H2)*n(idx_HEj) &
    -k(14)*n(idx_H2)*n(idx_HEj) &
    +k(17)*n(idx_Hk)*n(idx_H) &
    +k(20)*n(idx_H2j)*n(idx_H) &
    -k(21)*n(idx_H2)*n(idx_Hj) &
    -k(22)*n(idx_H2)*n(idx_Hj) &
    -k(23)*n(idx_H2)*n(idx_E) &
    -k(24)*n(idx_H2)*n(idx_H) &
    +k(32)*n(idx_H2j)*n(idx_Hk) &
    -2.d0*k(33)*n(idx_H2)*n(idx_H2) &
    +k(33)*n(idx_H2)*n(idx_H2) &
    +k(34)*n(idx_H)*n(idx_H)*n(idx_HE) &
    +k(35)*n(idx_H)*n(idx_H)*n(idx_H) &
    -k(36)*n(idx_H2)*n(idx_H)*n(idx_H) &
    +2.d0*k(36)*n(idx_H2)*n(idx_H)*n(idx_H) &
    -k(53)*n(idx_HOCj)*n(idx_H2) &
    +k(53)*n(idx_HOCj)*n(idx_H2) &
    -k(56)*n(idx_C)*n(idx_H2) &
    +k(57)*n(idx_CH)*n(idx_H) &
    -k(58)*n(idx_CH)*n(idx_H2) &
    +k(63)*n(idx_CH2)*n(idx_H) &
    +k(65)*n(idx_CH2)*n(idx_O) &
    -k(70)*n(idx_O)*n(idx_H2) &
    +k(71)*n(idx_OH)*n(idx_H) &
    +k(72)*n(idx_OH)*n(idx_H) &
    -k(73)*n(idx_H2)*n(idx_OH) &
    +k(79)*n(idx_H2O)*n(idx_H) &
    -k(81)*n(idx_O2)*n(idx_H2) &
    -k(85)*n(idx_H2j)*n(idx_H2) &
    +k(86)*n(idx_H3j)*n(idx_H) &
    +k(88)*n(idx_C)*n(idx_H3j) &
    -k(90)*n(idx_Cj)*n(idx_H2) &
    +k(91)*n(idx_CHj)*n(idx_H) &
    -k(92)*n(idx_CHj)*n(idx_H2) &
    +k(94)*n(idx_CH2j)*n(idx_H) &
    -k(95)*n(idx_CH2j)*n(idx_H2) &
    +k(97)*n(idx_CH3j)*n(idx_H) &
    +k(98)*n(idx_CH3j)*n(idx_O) &
    +k(99)*n(idx_CH3j)*n(idx_O) &
    -k(101)*n(idx_Oj)*n(idx_H2) &
    +k(103)*n(idx_O)*n(idx_H3j) &
    +k(105)*n(idx_OH)*n(idx_H3j) &
    +k(106)*n(idx_OH)*n(idx_H3j) &
    -k(109)*n(idx_OHj)*n(idx_H2) &
    -k(110)*n(idx_H2Oj)*n(idx_H2) &
    +k(111)*n(idx_H2O)*n(idx_H3j) &
    +k(112)*n(idx_H2O)*n(idx_H3j) &
    +k(117)*n(idx_H3Oj)*n(idx_C) &
    +k(123)*n(idx_CO)*n(idx_H3j) &
    +k(124)*n(idx_CO)*n(idx_H3j) &
    +k(125)*n(idx_CO)*n(idx_H3j) &
    +k(126)*n(idx_CO)*n(idx_H3j) &
    +k(132)*n(idx_CH2)*n(idx_Hj) &
    +k(133)*n(idx_CH2)*n(idx_Hj) &
    +k(136)*n(idx_CH2)*n(idx_HEj) &
    +k(137)*n(idx_CH2)*n(idx_HEj) &
    +k(162)*n(idx_H3j)*n(idx_E) &
    +k(166)*n(idx_CH2j)*n(idx_E) &
    +k(169)*n(idx_CH3j)*n(idx_E) &
    +k(172)*n(idx_H2Oj)*n(idx_E) &
    +k(176)*n(idx_H3Oj)*n(idx_E) &
    +k(178)*n(idx_H3Oj)*n(idx_E) &
    -k(188)*n(idx_Ck)*n(idx_H2) &
    -k(191)*n(idx_Ok)*n(idx_H2) &
    -k(193)*n(idx_H2)*n(idx_Hj) &
    -k(194)*n(idx_H2)*n(idx_Hj) &
    -k(197)*n(idx_C)*n(idx_H2) &
    -k(201)*n(idx_Cj)*n(idx_H2) &
    -k(218)*n(idx_H2) &
    -k(229)*n(idx_H2) &
    -k(231)*n(idx_H2) &
    +k(233)*n(idx_H3j) &
    +k(241)*n(idx_CH3j) &
    +k(246)*n(idx_H2Oj) &
    +k(251)*n(idx_H3Oj) &
    -k(258)*n(idx_H2) &
    -k(259)*n(idx_H2) &
    -k(260)*n(idx_H2) &
    -k(270)*n(idx_H2) + nH2dust
!C
dn(idx_C) = &
    +k(37)*n(idx_Cj)*n(idx_E) &
    +k(38)*n(idx_Cj)*n(idx_E) &
    +k(39)*n(idx_Cj)*n(idx_E) &
    -k(42)*n(idx_C)*n(idx_E) &
    -k(47)*n(idx_C)*n(idx_Hj) &
    +k(48)*n(idx_Cj)*n(idx_H) &
    -k(49)*n(idx_C)*n(idx_HEj) &
    -k(50)*n(idx_C)*n(idx_HEj) &
    -k(51)*n(idx_C)*n(idx_HEj) &
    -k(56)*n(idx_C)*n(idx_H2) &
    +k(57)*n(idx_CH)*n(idx_H) &
    -k(59)*n(idx_CH)*n(idx_C) &
    +k(62)*n(idx_CH)*n(idx_O) &
    +k(68)*n(idx_C2)*n(idx_O) &
    +k(69)*n(idx_C2)*n(idx_O) &
    -k(74)*n(idx_C)*n(idx_OH) &
    -k(75)*n(idx_C)*n(idx_OH) &
    -k(82)*n(idx_O2)*n(idx_C) &
    -k(83)*n(idx_O2)*n(idx_C) &
    +k(84)*n(idx_CO)*n(idx_H) &
    -k(87)*n(idx_C)*n(idx_H2j) &
    -k(88)*n(idx_C)*n(idx_H3j) &
    -k(89)*n(idx_C)*n(idx_H3j) &
    +k(100)*n(idx_C2)*n(idx_Oj) &
    +k(116)*n(idx_H2O)*n(idx_Cj) &
    -k(117)*n(idx_H3Oj)*n(idx_C) &
    -k(121)*n(idx_C)*n(idx_O2j) &
    -k(122)*n(idx_C)*n(idx_O2j) &
    -k(127)*n(idx_HCOj)*n(idx_C) &
    +k(140)*n(idx_C2)*n(idx_HEj) &
    +k(157)*n(idx_CO)*n(idx_HEj) &
    +k(159)*n(idx_Ck)*n(idx_Hj) &
    +k(164)*n(idx_CHj)*n(idx_E) &
    +k(166)*n(idx_CH2j)*n(idx_E) &
    +k(167)*n(idx_CH2j)*n(idx_E) &
    +k(180)*n(idx_COj)*n(idx_E) &
    +k(182)*n(idx_HCOj)*n(idx_E) &
    -k(184)*n(idx_Hk)*n(idx_C) &
    -k(192)*n(idx_Ok)*n(idx_C) &
    -k(195)*n(idx_C)*n(idx_E) &
    -k(196)*n(idx_C)*n(idx_H) &
    -k(197)*n(idx_C)*n(idx_H2) &
    -2.d0*k(198)*n(idx_C)*n(idx_C) &
    -k(199)*n(idx_C)*n(idx_O) &
    -k(217)*n(idx_C) &
    +k(220)*n(idx_CH) &
    +2.d0*k(222)*n(idx_C2) &
    +k(230)*n(idx_CO) &
    +k(235)*n(idx_Ck) &
    +k(236)*n(idx_CHj) &
    +k(255)*n(idx_CO) &
    +2.d0*k(257)*n(idx_C2) &
    -k(261)*n(idx_C) &
    +k(262)*n(idx_CH) &
    -2.d0*k(271)*n(idx_C)*n(idx_C) &
    -2.d0*k(272)*n(idx_C)*n(idx_C) &
    -k(273)*n(idx_C)*n(idx_O) &
    -k(274)*n(idx_C)*n(idx_O) &
    -k(277)*n(idx_C)*n(idx_Oj) &
    -k(278)*n(idx_C)*n(idx_Oj)
!O
dn(idx_O) = &
    +k(40)*n(idx_Oj)*n(idx_E) &
    +k(41)*n(idx_Oj)*n(idx_E) &
    -k(43)*n(idx_O)*n(idx_E) &
    +k(44)*n(idx_Oj)*n(idx_H) &
    -k(45)*n(idx_O)*n(idx_Hj) &
    -k(46)*n(idx_O)*n(idx_HEj) &
    +k(52)*n(idx_OH)*n(idx_H) &
    -k(60)*n(idx_CH)*n(idx_O) &
    -k(61)*n(idx_CH)*n(idx_O) &
    -k(62)*n(idx_CH)*n(idx_O) &
    -k(64)*n(idx_CH2)*n(idx_O) &
    -k(65)*n(idx_CH2)*n(idx_O) &
    -k(66)*n(idx_CH2)*n(idx_O) &
    -k(67)*n(idx_CH2)*n(idx_O) &
    -k(68)*n(idx_C2)*n(idx_O) &
    -k(69)*n(idx_C2)*n(idx_O) &
    -k(70)*n(idx_O)*n(idx_H2) &
    +k(71)*n(idx_OH)*n(idx_H) &
    +k(72)*n(idx_OH)*n(idx_H) &
    -k(76)*n(idx_O)*n(idx_OH) &
    -k(77)*n(idx_O)*n(idx_OH) &
    +k(78)*n(idx_OH)*n(idx_OH) &
    +k(80)*n(idx_O2)*n(idx_H) &
    +k(82)*n(idx_O2)*n(idx_C) &
    +k(83)*n(idx_O2)*n(idx_C) &
    -k(93)*n(idx_CHj)*n(idx_O) &
    -k(96)*n(idx_CH2j)*n(idx_O) &
    -k(98)*n(idx_CH3j)*n(idx_O) &
    -k(99)*n(idx_CH3j)*n(idx_O) &
    -k(102)*n(idx_O)*n(idx_H2j) &
    -k(103)*n(idx_O)*n(idx_H3j) &
    -k(104)*n(idx_O)*n(idx_H3j) &
    +k(118)*n(idx_O2)*n(idx_Cj) &
    +k(121)*n(idx_C)*n(idx_O2j) &
    +k(155)*n(idx_O2)*n(idx_HEj) &
    +k(156)*n(idx_CO)*n(idx_HEj) &
    +k(160)*n(idx_Ok)*n(idx_Hj) &
    +k(171)*n(idx_OHj)*n(idx_E) &
    +k(172)*n(idx_H2Oj)*n(idx_E) &
    +k(174)*n(idx_H2Oj)*n(idx_E) &
    +k(176)*n(idx_H3Oj)*n(idx_E) &
    +2.d0*k(179)*n(idx_O2j)*n(idx_E) &
    +k(180)*n(idx_COj)*n(idx_E) &
    -k(185)*n(idx_Hk)*n(idx_O) &
    -k(189)*n(idx_Ck)*n(idx_O) &
    -k(199)*n(idx_C)*n(idx_O) &
    -k(202)*n(idx_Cj)*n(idx_O) &
    -k(203)*n(idx_Cj)*n(idx_O) &
    -k(204)*n(idx_O)*n(idx_E) &
    -k(205)*n(idx_O)*n(idx_H) &
    -2.d0*k(206)*n(idx_O)*n(idx_O) &
    -k(216)*n(idx_O) &
    +k(223)*n(idx_OH) &
    +2.d0*k(228)*n(idx_O2) &
    +k(230)*n(idx_CO) &
    +k(242)*n(idx_Ok) &
    +k(243)*n(idx_OHj) &
    +k(244)*n(idx_H2Oj) &
    -k(254)*n(idx_O) &
    +k(255)*n(idx_CO) &
    +2.d0*k(263)*n(idx_O2) &
    +k(265)*n(idx_OH) &
    -k(273)*n(idx_C)*n(idx_O) &
    -k(274)*n(idx_C)*n(idx_O) &
    -k(275)*n(idx_Cj)*n(idx_O) &
    -k(276)*n(idx_Cj)*n(idx_O) &
    -k(279)*n(idx_H)*n(idx_O) &
    -2.d0*k(281)*n(idx_O)*n(idx_O)
!OH
dn(idx_OH) = &
    -k(52)*n(idx_OH)*n(idx_H) &
    +k(62)*n(idx_CH)*n(idx_O) &
    +k(67)*n(idx_CH2)*n(idx_O) &
    +k(70)*n(idx_O)*n(idx_H2) &
    -k(71)*n(idx_OH)*n(idx_H) &
    -k(72)*n(idx_OH)*n(idx_H) &
    -k(73)*n(idx_H2)*n(idx_OH) &
    -k(74)*n(idx_C)*n(idx_OH) &
    -k(75)*n(idx_C)*n(idx_OH) &
    -k(76)*n(idx_O)*n(idx_OH) &
    -k(77)*n(idx_O)*n(idx_OH) &
    -2.d0*k(78)*n(idx_OH)*n(idx_OH) &
    +k(79)*n(idx_H2O)*n(idx_H) &
    +k(80)*n(idx_O2)*n(idx_H) &
    +2.d0*k(81)*n(idx_O2)*n(idx_H2) &
    +k(84)*n(idx_CO)*n(idx_H) &
    -k(105)*n(idx_OH)*n(idx_H3j) &
    -k(106)*n(idx_OH)*n(idx_H3j) &
    -k(107)*n(idx_OH)*n(idx_Cj) &
    -k(108)*n(idx_OH)*n(idx_Cj) &
    +k(120)*n(idx_O2)*n(idx_CH2j) &
    -k(141)*n(idx_OH)*n(idx_Hj) &
    -k(142)*n(idx_OH)*n(idx_Hj) &
    -k(143)*n(idx_OH)*n(idx_HEj) &
    -k(144)*n(idx_OH)*n(idx_HEj) &
    +k(147)*n(idx_H2O)*n(idx_HEj) &
    +k(148)*n(idx_H2O)*n(idx_HEj) &
    +k(173)*n(idx_H2Oj)*n(idx_E) &
    +k(175)*n(idx_H3Oj)*n(idx_E) &
    +k(178)*n(idx_H3Oj)*n(idx_E) &
    +k(182)*n(idx_HCOj)*n(idx_E) &
    +k(185)*n(idx_Hk)*n(idx_O) &
    -k(186)*n(idx_Hk)*n(idx_OH) &
    +k(190)*n(idx_Ok)*n(idx_H) &
    +k(205)*n(idx_O)*n(idx_H) &
    -k(207)*n(idx_OH)*n(idx_H) &
    -k(223)*n(idx_OH) &
    -k(224)*n(idx_OH) &
    +k(225)*n(idx_H2O) &
    +k(245)*n(idx_H2Oj) &
    +k(249)*n(idx_H3Oj) &
    -k(265)*n(idx_OH) &
    +k(267)*n(idx_H2O) &
    +k(279)*n(idx_H)*n(idx_O) &
    -k(280)*n(idx_OH)*n(idx_H)
!CO_GAS
dn(idx_CO) = &
    dnChem_CO &
    -n(idx_CO)*(k(208)+k(212)) &
    +k(212)*n(idx_CO_total)
!CH
dn(idx_CH) = &
    +k(56)*n(idx_C)*n(idx_H2) &
    -k(57)*n(idx_CH)*n(idx_H) &
    -k(58)*n(idx_CH)*n(idx_H2) &
    -k(59)*n(idx_CH)*n(idx_C) &
    -k(60)*n(idx_CH)*n(idx_O) &
    -k(61)*n(idx_CH)*n(idx_O) &
    -k(62)*n(idx_CH)*n(idx_O) &
    +k(63)*n(idx_CH2)*n(idx_H) &
    +k(67)*n(idx_CH2)*n(idx_O) &
    -k(130)*n(idx_CH)*n(idx_Hj) &
    -k(131)*n(idx_CH)*n(idx_Hj) &
    +k(165)*n(idx_CH2j)*n(idx_E) &
    +k(169)*n(idx_CH3j)*n(idx_E) &
    +k(170)*n(idx_CH3j)*n(idx_E) &
    +k(184)*n(idx_Hk)*n(idx_C) &
    +k(187)*n(idx_Ck)*n(idx_H) &
    +k(196)*n(idx_C)*n(idx_H) &
    -k(220)*n(idx_CH) &
    -k(221)*n(idx_CH) &
    +k(237)*n(idx_CH2) &
    -k(262)*n(idx_CH)
!CH2
dn(idx_CH2) = &
    +k(58)*n(idx_CH)*n(idx_H2) &
    -k(63)*n(idx_CH2)*n(idx_H) &
    -k(64)*n(idx_CH2)*n(idx_O) &
    -k(65)*n(idx_CH2)*n(idx_O) &
    -k(66)*n(idx_CH2)*n(idx_O) &
    -k(67)*n(idx_CH2)*n(idx_O) &
    -k(132)*n(idx_CH2)*n(idx_Hj) &
    -k(133)*n(idx_CH2)*n(idx_Hj) &
    -k(134)*n(idx_CH2)*n(idx_Hj) &
    -k(135)*n(idx_CH2)*n(idx_Hj) &
    -k(136)*n(idx_CH2)*n(idx_HEj) &
    -k(137)*n(idx_CH2)*n(idx_HEj) &
    -k(138)*n(idx_CH2)*n(idx_HEj) &
    -k(139)*n(idx_CH2)*n(idx_HEj) &
    +k(168)*n(idx_CH3j)*n(idx_E) &
    +k(188)*n(idx_Ck)*n(idx_H2) &
    +k(197)*n(idx_C)*n(idx_H2) &
    -k(237)*n(idx_CH2) &
    -k(238)*n(idx_CH2) &
    -k(266)*n(idx_CH2)
!C2
dn(idx_C2) = &
    +k(59)*n(idx_CH)*n(idx_C) &
    -k(68)*n(idx_C2)*n(idx_O) &
    -k(69)*n(idx_C2)*n(idx_O) &
    -k(100)*n(idx_C2)*n(idx_Oj) &
    -k(140)*n(idx_C2)*n(idx_HEj) &
    +k(198)*n(idx_C)*n(idx_C) &
    -k(222)*n(idx_C2) &
    -k(257)*n(idx_C2) &
    +k(271)*n(idx_C)*n(idx_C) &
    +k(272)*n(idx_C)*n(idx_C)
!HCO
dn(idx_HCO) = &
    +k(66)*n(idx_CH2)*n(idx_O) &
    -k(268)*n(idx_HCO) &
    -k(269)*n(idx_HCO)
!H2O_GAS
dn(idx_H2O) = &
    dnChem_H2O &
    -n(idx_H2O)*(k(210)+k(211)) &
    +k(211)*n(idx_H2O_total)
!O2
dn(idx_O2) = &
    +k(76)*n(idx_O)*n(idx_OH) &
    +k(77)*n(idx_O)*n(idx_OH) &
    -k(80)*n(idx_O2)*n(idx_H) &
    -k(81)*n(idx_O2)*n(idx_H2) &
    -k(82)*n(idx_O2)*n(idx_C) &
    -k(83)*n(idx_O2)*n(idx_C) &
    -k(118)*n(idx_O2)*n(idx_Cj) &
    -k(119)*n(idx_O2)*n(idx_Cj) &
    -k(120)*n(idx_O2)*n(idx_CH2j) &
    +k(122)*n(idx_C)*n(idx_O2j) &
    -k(153)*n(idx_O2)*n(idx_Hj) &
    -k(154)*n(idx_O2)*n(idx_HEj) &
    -k(155)*n(idx_O2)*n(idx_HEj) &
    +k(206)*n(idx_O)*n(idx_O) &
    -k(227)*n(idx_O2) &
    -k(228)*n(idx_O2) &
    -k(263)*n(idx_O2) &
    -k(264)*n(idx_O2) &
    +k(281)*n(idx_O)*n(idx_O)
!CO_TOTAL
dn(idx_CO_total) = &
    dnChem_CO
!H2O_TOTAL
dn(idx_H2O_total) = &
    dnChem_H2O
!H+
dn(idx_Hj) = &
    +k(1)*n(idx_H)*n(idx_E) &
    -k(2)*n(idx_Hj)*n(idx_E) &
    -k(3)*n(idx_Hj)*n(idx_E) &
    +k(8)*n(idx_HEj)*n(idx_H) &
    -k(9)*n(idx_HE)*n(idx_Hj) &
    -k(10)*n(idx_HE)*n(idx_Hj) &
    +k(13)*n(idx_H2)*n(idx_HEj) &
    -k(18)*n(idx_H)*n(idx_Hj) &
    -k(19)*n(idx_H)*n(idx_Hj) &
    +k(20)*n(idx_H2j)*n(idx_H) &
    -k(21)*n(idx_H2)*n(idx_Hj) &
    -k(22)*n(idx_H2)*n(idx_Hj) &
    -k(28)*n(idx_Hk)*n(idx_Hj) &
    -k(29)*n(idx_Hk)*n(idx_Hj) &
    +k(44)*n(idx_Oj)*n(idx_H) &
    -k(45)*n(idx_O)*n(idx_Hj) &
    -k(47)*n(idx_C)*n(idx_Hj) &
    +k(48)*n(idx_Cj)*n(idx_H) &
    -k(130)*n(idx_CH)*n(idx_Hj) &
    -k(131)*n(idx_CH)*n(idx_Hj) &
    -k(132)*n(idx_CH2)*n(idx_Hj) &
    -k(133)*n(idx_CH2)*n(idx_Hj) &
    -k(134)*n(idx_CH2)*n(idx_Hj) &
    -k(135)*n(idx_CH2)*n(idx_Hj) &
    -k(141)*n(idx_OH)*n(idx_Hj) &
    -k(142)*n(idx_OH)*n(idx_Hj) &
    -k(145)*n(idx_H2O)*n(idx_Hj) &
    -k(146)*n(idx_H2O)*n(idx_Hj) &
    +k(147)*n(idx_H2O)*n(idx_HEj) &
    +k(148)*n(idx_H2O)*n(idx_HEj) &
    -k(153)*n(idx_O2)*n(idx_Hj) &
    +k(158)*n(idx_COj)*n(idx_H) &
    -k(159)*n(idx_Ck)*n(idx_Hj) &
    -k(160)*n(idx_Ok)*n(idx_Hj) &
    -k(193)*n(idx_H2)*n(idx_Hj) &
    +k(193)*n(idx_H2)*n(idx_Hj) &
    -k(194)*n(idx_H2)*n(idx_Hj) &
    +k(213)*n(idx_H) &
    +k(229)*n(idx_H2) &
    +k(232)*n(idx_H2j) &
    +k(233)*n(idx_H3j) &
    +k(236)*n(idx_CHj) &
    +k(243)*n(idx_OHj) &
    +k(245)*n(idx_H2Oj) &
    +k(248)*n(idx_H3Oj) &
    +k(252)*n(idx_H) &
    +k(259)*n(idx_H2) &
    +k(270)*n(idx_H2)
!HE+
dn(idx_HEj) = &
    +k(4)*n(idx_HE)*n(idx_E) &
    -k(5)*n(idx_HEj)*n(idx_E) &
    -k(6)*n(idx_HEj)*n(idx_E) &
    -k(7)*n(idx_HEj)*n(idx_E) &
    -k(8)*n(idx_HEj)*n(idx_H) &
    +k(9)*n(idx_HE)*n(idx_Hj) &
    +k(10)*n(idx_HE)*n(idx_Hj) &
    -k(12)*n(idx_H2)*n(idx_HEj) &
    -k(13)*n(idx_H2)*n(idx_HEj) &
    -k(14)*n(idx_H2)*n(idx_HEj) &
    +k(14)*n(idx_H2)*n(idx_HEj) &
    +k(15)*n(idx_HEjj)*n(idx_E) &
    -k(46)*n(idx_O)*n(idx_HEj) &
    -k(49)*n(idx_C)*n(idx_HEj) &
    -k(50)*n(idx_C)*n(idx_HEj) &
    -k(51)*n(idx_C)*n(idx_HEj) &
    -k(136)*n(idx_CH2)*n(idx_HEj) &
    -k(137)*n(idx_CH2)*n(idx_HEj) &
    -k(138)*n(idx_CH2)*n(idx_HEj) &
    -k(139)*n(idx_CH2)*n(idx_HEj) &
    -k(140)*n(idx_C2)*n(idx_HEj) &
    -k(143)*n(idx_OH)*n(idx_HEj) &
    -k(144)*n(idx_OH)*n(idx_HEj) &
    -k(147)*n(idx_H2O)*n(idx_HEj) &
    -k(148)*n(idx_H2O)*n(idx_HEj) &
    -k(149)*n(idx_H2O)*n(idx_HEj) &
    -k(150)*n(idx_H2O)*n(idx_HEj) &
    -k(151)*n(idx_H2O)*n(idx_HEj) &
    -k(152)*n(idx_H2O)*n(idx_HEj) &
    -k(154)*n(idx_O2)*n(idx_HEj) &
    -k(155)*n(idx_O2)*n(idx_HEj) &
    -k(156)*n(idx_CO)*n(idx_HEj) &
    -k(157)*n(idx_CO)*n(idx_HEj) &
    -k(161)*n(idx_HEj)*n(idx_Hk) &
    +k(214)*n(idx_HE) &
    -k(215)*n(idx_HEj) &
    +k(253)*n(idx_HE)
!H2+
dn(idx_H2j) = &
    +k(12)*n(idx_H2)*n(idx_HEj) &
    +k(18)*n(idx_H)*n(idx_Hj) &
    +k(19)*n(idx_H)*n(idx_Hj) &
    -k(20)*n(idx_H2j)*n(idx_H) &
    +k(21)*n(idx_H2)*n(idx_Hj) &
    +k(22)*n(idx_H2)*n(idx_Hj) &
    +k(29)*n(idx_Hk)*n(idx_Hj) &
    -k(30)*n(idx_H2j)*n(idx_E) &
    -k(31)*n(idx_H2j)*n(idx_E) &
    -k(32)*n(idx_H2j)*n(idx_Hk) &
    -k(85)*n(idx_H2j)*n(idx_H2) &
    +k(86)*n(idx_H3j)*n(idx_H) &
    -k(87)*n(idx_C)*n(idx_H2j) &
    -k(102)*n(idx_O)*n(idx_H2j) &
    +k(218)*n(idx_H2) &
    -k(232)*n(idx_H2j) &
    +k(234)*n(idx_H3j) &
    +k(244)*n(idx_H2Oj) &
    +k(249)*n(idx_H3Oj) &
    +k(260)*n(idx_H2)
!C+
dn(idx_Cj) = &
    -k(37)*n(idx_Cj)*n(idx_E) &
    -k(38)*n(idx_Cj)*n(idx_E) &
    -k(39)*n(idx_Cj)*n(idx_E) &
    +k(42)*n(idx_C)*n(idx_E) &
    +k(47)*n(idx_C)*n(idx_Hj) &
    -k(48)*n(idx_Cj)*n(idx_H) &
    +k(49)*n(idx_C)*n(idx_HEj) &
    +k(50)*n(idx_C)*n(idx_HEj) &
    +k(51)*n(idx_C)*n(idx_HEj) &
    -k(90)*n(idx_Cj)*n(idx_H2) &
    +k(91)*n(idx_CHj)*n(idx_H) &
    -k(107)*n(idx_OH)*n(idx_Cj) &
    -k(108)*n(idx_OH)*n(idx_Cj) &
    -k(113)*n(idx_H2O)*n(idx_Cj) &
    -k(114)*n(idx_H2O)*n(idx_Cj) &
    -k(115)*n(idx_H2O)*n(idx_Cj) &
    -k(116)*n(idx_H2O)*n(idx_Cj) &
    -k(118)*n(idx_O2)*n(idx_Cj) &
    -k(119)*n(idx_O2)*n(idx_Cj) &
    +k(122)*n(idx_C)*n(idx_O2j) &
    +k(136)*n(idx_CH2)*n(idx_HEj) &
    +k(137)*n(idx_CH2)*n(idx_HEj) &
    +k(140)*n(idx_C2)*n(idx_HEj) &
    +k(156)*n(idx_CO)*n(idx_HEj) &
    -k(200)*n(idx_Cj)*n(idx_H) &
    -k(201)*n(idx_Cj)*n(idx_H2) &
    -k(202)*n(idx_Cj)*n(idx_O) &
    -k(203)*n(idx_Cj)*n(idx_O) &
    +k(217)*n(idx_C) &
    +k(261)*n(idx_C) &
    -k(275)*n(idx_Cj)*n(idx_O) &
    -k(276)*n(idx_Cj)*n(idx_O)
!O+
dn(idx_Oj) = &
    -k(40)*n(idx_Oj)*n(idx_E) &
    -k(41)*n(idx_Oj)*n(idx_E) &
    +k(43)*n(idx_O)*n(idx_E) &
    -k(44)*n(idx_Oj)*n(idx_H) &
    +k(45)*n(idx_O)*n(idx_Hj) &
    +k(46)*n(idx_O)*n(idx_HEj) &
    -k(100)*n(idx_C2)*n(idx_Oj) &
    -k(101)*n(idx_Oj)*n(idx_H2) &
    +k(119)*n(idx_O2)*n(idx_Cj) &
    +k(143)*n(idx_OH)*n(idx_HEj) &
    +k(144)*n(idx_OH)*n(idx_HEj) &
    +k(155)*n(idx_O2)*n(idx_HEj) &
    +k(157)*n(idx_CO)*n(idx_HEj) &
    +k(216)*n(idx_O) &
    +k(246)*n(idx_H2Oj) &
    +k(254)*n(idx_O) &
    -k(277)*n(idx_C)*n(idx_Oj) &
    -k(278)*n(idx_C)*n(idx_Oj)
!HOC+
dn(idx_HOCj) = &
    -k(53)*n(idx_HOCj)*n(idx_H2) &
    -k(54)*n(idx_HOCj)*n(idx_CO) &
    -k(55)*n(idx_HOCj)*n(idx_CO) &
    +k(98)*n(idx_CH3j)*n(idx_O) &
    +k(113)*n(idx_H2O)*n(idx_Cj) &
    +k(125)*n(idx_CO)*n(idx_H3j) &
    +k(126)*n(idx_CO)*n(idx_H3j) &
    -k(183)*n(idx_HOCj)*n(idx_E)
!HCO+
dn(idx_HCOj) = &
    +k(53)*n(idx_HOCj)*n(idx_H2) &
    +k(54)*n(idx_HOCj)*n(idx_CO) &
    +k(55)*n(idx_HOCj)*n(idx_CO) &
    +k(61)*n(idx_CH)*n(idx_O) &
    +k(96)*n(idx_CH2j)*n(idx_O) &
    +k(99)*n(idx_CH3j)*n(idx_O) &
    +k(114)*n(idx_H2O)*n(idx_Cj) &
    +k(115)*n(idx_H2O)*n(idx_Cj) &
    +k(117)*n(idx_H3Oj)*n(idx_C) &
    +k(120)*n(idx_O2)*n(idx_CH2j) &
    +k(123)*n(idx_CO)*n(idx_H3j) &
    +k(124)*n(idx_CO)*n(idx_H3j) &
    -k(127)*n(idx_HCOj)*n(idx_C) &
    -k(128)*n(idx_HCOj)*n(idx_H2O) &
    -k(129)*n(idx_HCOj)*n(idx_H2O) &
    -k(181)*n(idx_HCOj)*n(idx_E) &
    -k(182)*n(idx_HCOj)*n(idx_E) &
    +k(269)*n(idx_HCO)
!H3+
dn(idx_H3j) = &
    +k(85)*n(idx_H2j)*n(idx_H2) &
    -k(86)*n(idx_H3j)*n(idx_H) &
    -k(88)*n(idx_C)*n(idx_H3j) &
    -k(89)*n(idx_C)*n(idx_H3j) &
    -k(103)*n(idx_O)*n(idx_H3j) &
    -k(104)*n(idx_O)*n(idx_H3j) &
    -k(105)*n(idx_OH)*n(idx_H3j) &
    -k(106)*n(idx_OH)*n(idx_H3j) &
    -k(111)*n(idx_H2O)*n(idx_H3j) &
    -k(112)*n(idx_H2O)*n(idx_H3j) &
    -k(123)*n(idx_CO)*n(idx_H3j) &
    -k(124)*n(idx_CO)*n(idx_H3j) &
    -k(125)*n(idx_CO)*n(idx_H3j) &
    -k(126)*n(idx_CO)*n(idx_H3j) &
    -k(162)*n(idx_H3j)*n(idx_E) &
    -k(163)*n(idx_H3j)*n(idx_E) &
    +k(194)*n(idx_H2)*n(idx_Hj) &
    -k(233)*n(idx_H3j) &
    -k(234)*n(idx_H3j)
!CH+
dn(idx_CHj) = &
    +k(87)*n(idx_C)*n(idx_H2j) &
    +k(88)*n(idx_C)*n(idx_H3j) &
    +k(90)*n(idx_Cj)*n(idx_H2) &
    -k(91)*n(idx_CHj)*n(idx_H) &
    -k(92)*n(idx_CHj)*n(idx_H2) &
    -k(93)*n(idx_CHj)*n(idx_O) &
    +k(94)*n(idx_CH2j)*n(idx_H) &
    +k(127)*n(idx_HCOj)*n(idx_C) &
    +k(130)*n(idx_CH)*n(idx_Hj) &
    +k(131)*n(idx_CH)*n(idx_Hj) &
    +k(132)*n(idx_CH2)*n(idx_Hj) &
    +k(133)*n(idx_CH2)*n(idx_Hj) &
    +k(138)*n(idx_CH2)*n(idx_HEj) &
    +k(139)*n(idx_CH2)*n(idx_HEj) &
    -k(164)*n(idx_CHj)*n(idx_E) &
    +k(200)*n(idx_Cj)*n(idx_H) &
    +k(221)*n(idx_CH) &
    -k(236)*n(idx_CHj) &
    +k(239)*n(idx_CH2j) &
    +k(241)*n(idx_CH3j)
!CH2+
dn(idx_CH2j) = &
    +k(89)*n(idx_C)*n(idx_H3j) &
    +k(92)*n(idx_CHj)*n(idx_H2) &
    -k(94)*n(idx_CH2j)*n(idx_H) &
    -k(95)*n(idx_CH2j)*n(idx_H2) &
    -k(96)*n(idx_CH2j)*n(idx_O) &
    +k(97)*n(idx_CH3j)*n(idx_H) &
    -k(120)*n(idx_O2)*n(idx_CH2j) &
    +k(134)*n(idx_CH2)*n(idx_Hj) &
    +k(135)*n(idx_CH2)*n(idx_Hj) &
    -k(165)*n(idx_CH2j)*n(idx_E) &
    -k(166)*n(idx_CH2j)*n(idx_E) &
    -k(167)*n(idx_CH2j)*n(idx_E) &
    +k(201)*n(idx_Cj)*n(idx_H2) &
    +k(238)*n(idx_CH2) &
    -k(239)*n(idx_CH2j) &
    +k(240)*n(idx_CH3j) &
    +k(266)*n(idx_CH2)
!CO+
dn(idx_COj) = &
    +k(93)*n(idx_CHj)*n(idx_O) &
    +k(100)*n(idx_C2)*n(idx_Oj) &
    +k(107)*n(idx_OH)*n(idx_Cj) &
    +k(108)*n(idx_OH)*n(idx_Cj) &
    +k(118)*n(idx_O2)*n(idx_Cj) &
    +k(121)*n(idx_C)*n(idx_O2j) &
    -k(158)*n(idx_COj)*n(idx_H) &
    -k(180)*n(idx_COj)*n(idx_E) &
    +k(202)*n(idx_Cj)*n(idx_O) &
    +k(203)*n(idx_Cj)*n(idx_O) &
    +k(256)*n(idx_CO) &
    +k(275)*n(idx_Cj)*n(idx_O) &
    +k(276)*n(idx_Cj)*n(idx_O) &
    +k(277)*n(idx_C)*n(idx_Oj) &
    +k(278)*n(idx_C)*n(idx_Oj)
!CH3+
dn(idx_CH3j) = &
    +k(95)*n(idx_CH2j)*n(idx_H2) &
    -k(97)*n(idx_CH3j)*n(idx_H) &
    -k(98)*n(idx_CH3j)*n(idx_O) &
    -k(99)*n(idx_CH3j)*n(idx_O) &
    -k(168)*n(idx_CH3j)*n(idx_E) &
    -k(169)*n(idx_CH3j)*n(idx_E) &
    -k(170)*n(idx_CH3j)*n(idx_E) &
    -k(240)*n(idx_CH3j) &
    -k(241)*n(idx_CH3j)
!OH+
dn(idx_OHj) = &
    +k(101)*n(idx_Oj)*n(idx_H2) &
    +k(102)*n(idx_O)*n(idx_H2j) &
    +k(103)*n(idx_O)*n(idx_H3j) &
    -k(109)*n(idx_OHj)*n(idx_H2) &
    +k(141)*n(idx_OH)*n(idx_Hj) &
    +k(142)*n(idx_OH)*n(idx_Hj) &
    +k(149)*n(idx_H2O)*n(idx_HEj) &
    +k(150)*n(idx_H2O)*n(idx_HEj) &
    -k(171)*n(idx_OHj)*n(idx_E) &
    +k(224)*n(idx_OH) &
    -k(243)*n(idx_OHj) &
    +k(247)*n(idx_H2Oj) &
    +k(251)*n(idx_H3Oj)
!H2O+
dn(idx_H2Oj) = &
    +k(104)*n(idx_O)*n(idx_H3j) &
    +k(105)*n(idx_OH)*n(idx_H3j) &
    +k(106)*n(idx_OH)*n(idx_H3j) &
    +k(109)*n(idx_OHj)*n(idx_H2) &
    -k(110)*n(idx_H2Oj)*n(idx_H2) &
    +k(116)*n(idx_H2O)*n(idx_Cj) &
    +k(145)*n(idx_H2O)*n(idx_Hj) &
    +k(146)*n(idx_H2O)*n(idx_Hj) &
    +k(151)*n(idx_H2O)*n(idx_HEj) &
    +k(152)*n(idx_H2O)*n(idx_HEj) &
    -k(172)*n(idx_H2Oj)*n(idx_E) &
    -k(173)*n(idx_H2Oj)*n(idx_E) &
    -k(174)*n(idx_H2Oj)*n(idx_E) &
    +k(226)*n(idx_H2O) &
    -k(244)*n(idx_H2Oj) &
    -k(245)*n(idx_H2Oj) &
    -k(246)*n(idx_H2Oj) &
    -k(247)*n(idx_H2Oj) &
    +k(250)*n(idx_H3Oj)
!H3O+
dn(idx_H3Oj) = &
    +k(110)*n(idx_H2Oj)*n(idx_H2) &
    +k(111)*n(idx_H2O)*n(idx_H3j) &
    +k(112)*n(idx_H2O)*n(idx_H3j) &
    -k(117)*n(idx_H3Oj)*n(idx_C) &
    +k(128)*n(idx_HCOj)*n(idx_H2O) &
    +k(129)*n(idx_HCOj)*n(idx_H2O) &
    -k(175)*n(idx_H3Oj)*n(idx_E) &
    -k(176)*n(idx_H3Oj)*n(idx_E) &
    -k(177)*n(idx_H3Oj)*n(idx_E) &
    -k(178)*n(idx_H3Oj)*n(idx_E) &
    -k(248)*n(idx_H3Oj) &
    -k(249)*n(idx_H3Oj) &
    -k(250)*n(idx_H3Oj) &
    -k(251)*n(idx_H3Oj)
!O2+
dn(idx_O2j) = &
    -k(121)*n(idx_C)*n(idx_O2j) &
    -k(122)*n(idx_C)*n(idx_O2j) &
    +k(153)*n(idx_O2)*n(idx_Hj) &
    +k(154)*n(idx_O2)*n(idx_HEj) &
    -k(179)*n(idx_O2j)*n(idx_E) &
    +k(227)*n(idx_O2) &
    +k(264)*n(idx_O2)
!HE++
dn(idx_HEjj) = &
    +k(7)*n(idx_HEj)*n(idx_E) &
    -k(15)*n(idx_HEjj)*n(idx_E) &
    +k(215)*n(idx_HEj)
!CR
dn(idx_CR) = &
    0.d0

!g
dn(idx_g) = 0.d0

!Tgas
dn(idx_Tgas) = 0.d0

!dummy
dn(idx_dummy) = 0.d0

last_coe(:) = k(:)

end subroutine fex

!***************************
subroutine jes(neq, tt, n, j, ian, jan, pdj)
use krome_commons
use krome_subs
use krome_tabs
use krome_cooling
use krome_heating
use krome_constants
use krome_gadiab
use krome_getphys
implicit none
integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

nH2dust = 0.d0
Tgas = n(idx_Tgas)

k(:) = last_coe(:) !get rate coefficients

if(j==1) then
elseif(j==1) then
pdj(1) =  &
    -k(6)*n(idx_HEj)  &
    -k(178)*n(idx_H3Oj)  &
    -k(168)*n(idx_CH3j)  &
    -k(174)*n(idx_H2Oj)  &
    -k(195)*n(idx_C)  &
    -k(169)*n(idx_CH3j)  &
    +2.d0*k(4)*n(idx_HE)  &
    -k(25)*n(idx_Hk)  &
    -k(170)*n(idx_CH3j)  &
    -k(1)*n(idx_H)  &
    -k(166)*n(idx_CH2j)  &
    -k(16)*n(idx_H)  &
    -k(39)*n(idx_Cj)  &
    -k(183)*n(idx_HOCj)  &
    -k(4)*n(idx_HE)  &
    -k(180)*n(idx_COj)  &
    -k(5)*n(idx_HEj)  &
    -k(181)*n(idx_HCOj)  &
    -k(162)*n(idx_H3j)  &
    -k(38)*n(idx_Cj)  &
    -k(43)*n(idx_O)  &
    +2.d0*k(42)*n(idx_C)  &
    -k(40)*n(idx_Oj)  &
    -k(7)*n(idx_HEj)  &
    -k(176)*n(idx_H3Oj)  &
    -k(163)*n(idx_H3j)  &
    -k(171)*n(idx_OHj)  &
    -k(175)*n(idx_H3Oj)  &
    +2.d0*k(7)*n(idx_HEj)  &
    -k(177)*n(idx_H3Oj)  &
    -k(2)*n(idx_Hj)  &
    -k(167)*n(idx_CH2j)  &
    -k(41)*n(idx_Oj)  &
    -k(182)*n(idx_HCOj)  &
    -k(173)*n(idx_H2Oj)  &
    -k(42)*n(idx_C)  &
    -k(37)*n(idx_Cj)  &
    -k(3)*n(idx_Hj)  &
    -k(164)*n(idx_CHj)  &
    -k(31)*n(idx_H2j)  &
    +2.d0*k(25)*n(idx_Hk)  &
    +k(23)*n(idx_H2)  &
    -k(23)*n(idx_H2)  &
    -k(172)*n(idx_H2Oj)  &
    +2.d0*k(1)*n(idx_H)  &
    -k(179)*n(idx_O2j)  &
    -k(204)*n(idx_O)  &
    -k(30)*n(idx_H2j)  &
    +2.d0*k(43)*n(idx_O)  &
    -k(165)*n(idx_CH2j)  &
    -k(15)*n(idx_HEjj)
pdj(2) =  &
    +k(16)*n(idx_H)  &
    -k(25)*n(idx_Hk)
pdj(3) =  &
    +k(195)*n(idx_C)
pdj(4) =  &
    +k(204)*n(idx_O)
pdj(5) =  &
    +k(183)*n(idx_HOCj)  &
    +k(3)*n(idx_Hj)  &
    +k(173)*n(idx_H2Oj)  &
    +k(171)*n(idx_OHj)  &
    +2.d0*k(174)*n(idx_H2Oj)  &
    +k(177)*n(idx_H3Oj)  &
    +k(165)*n(idx_CH2j)  &
    -k(1)*n(idx_H)  &
    -k(16)*n(idx_H)  &
    +k(168)*n(idx_CH3j)  &
    +3.d0*k(163)*n(idx_H3j)  &
    +2.d0*k(167)*n(idx_CH2j)  &
    +k(162)*n(idx_H3j)  &
    +k(2)*n(idx_Hj)  &
    +k(176)*n(idx_H3Oj)  &
    +2.d0*k(31)*n(idx_H2j)  &
    +2.d0*k(30)*n(idx_H2j)  &
    +k(164)*n(idx_CHj)  &
    +2.d0*k(170)*n(idx_CH3j)  &
    +2.d0*k(175)*n(idx_H3Oj)  &
    +k(181)*n(idx_HCOj)  &
    +k(25)*n(idx_Hk)  &
    +2.d0*k(23)*n(idx_H2)
pdj(6) =  &
    -k(4)*n(idx_HE)  &
    +k(6)*n(idx_HEj)  &
    +k(5)*n(idx_HEj)
pdj(7) =  &
    +k(166)*n(idx_CH2j)  &
    -k(23)*n(idx_H2)  &
    +k(172)*n(idx_H2Oj)  &
    +k(162)*n(idx_H3j)  &
    +k(169)*n(idx_CH3j)  &
    +k(178)*n(idx_H3Oj)  &
    +k(176)*n(idx_H3Oj)
pdj(8) =  &
    +k(38)*n(idx_Cj)  &
    -k(42)*n(idx_C)  &
    +k(167)*n(idx_CH2j)  &
    +k(180)*n(idx_COj)  &
    +k(37)*n(idx_Cj)  &
    +k(164)*n(idx_CHj)  &
    +k(166)*n(idx_CH2j)  &
    -k(195)*n(idx_C)  &
    +k(182)*n(idx_HCOj)  &
    +k(39)*n(idx_Cj)
pdj(9) =  &
    -k(204)*n(idx_O)  &
    +k(180)*n(idx_COj)  &
    +k(41)*n(idx_Oj)  &
    -k(43)*n(idx_O)  &
    +k(171)*n(idx_OHj)  &
    +2.d0*k(179)*n(idx_O2j)  &
    +k(40)*n(idx_Oj)  &
    +k(172)*n(idx_H2Oj)  &
    +k(174)*n(idx_H2Oj)  &
    +k(176)*n(idx_H3Oj)
pdj(10) =  &
    +k(182)*n(idx_HCOj)  &
    +k(178)*n(idx_H3Oj)  &
    +k(175)*n(idx_H3Oj)  &
    +k(173)*n(idx_H2Oj)
pdj(11) =  &
    +k(183)*n(idx_HOCj)  &
    +k(181)*n(idx_HCOj)
pdj(12) =  &
    +k(165)*n(idx_CH2j)  &
    +k(169)*n(idx_CH3j)  &
    +k(170)*n(idx_CH3j)
pdj(13) =  &
    +k(168)*n(idx_CH3j)
pdj(16) =  &
    +k(177)*n(idx_H3Oj)
pdj(20) =  &
    -k(2)*n(idx_Hj)  &
    -k(3)*n(idx_Hj)  &
    +k(1)*n(idx_H)
pdj(21) =  &
    +k(4)*n(idx_HE)  &
    -k(5)*n(idx_HEj)  &
    +k(15)*n(idx_HEjj)  &
    -k(7)*n(idx_HEj)  &
    -k(6)*n(idx_HEj)
pdj(22) =  &
    -k(30)*n(idx_H2j)  &
    -k(31)*n(idx_H2j)
pdj(23) =  &
    -k(38)*n(idx_Cj)  &
    -k(37)*n(idx_Cj)  &
    +k(42)*n(idx_C)  &
    -k(39)*n(idx_Cj)
pdj(24) =  &
    -k(41)*n(idx_Oj)  &
    -k(40)*n(idx_Oj)  &
    +k(43)*n(idx_O)
pdj(25) =  &
    -k(183)*n(idx_HOCj)
pdj(26) =  &
    -k(182)*n(idx_HCOj)  &
    -k(181)*n(idx_HCOj)
pdj(27) =  &
    -k(162)*n(idx_H3j)  &
    -k(163)*n(idx_H3j)
pdj(28) =  &
    -k(164)*n(idx_CHj)
pdj(29) =  &
    -k(166)*n(idx_CH2j)  &
    -k(167)*n(idx_CH2j)  &
    -k(165)*n(idx_CH2j)
pdj(30) =  &
    -k(180)*n(idx_COj)
pdj(31) =  &
    -k(168)*n(idx_CH3j)  &
    -k(170)*n(idx_CH3j)  &
    -k(169)*n(idx_CH3j)
pdj(32) =  &
    -k(171)*n(idx_OHj)
pdj(33) =  &
    -k(174)*n(idx_H2Oj)  &
    -k(173)*n(idx_H2Oj)  &
    -k(172)*n(idx_H2Oj)
pdj(34) =  &
    -k(176)*n(idx_H3Oj)  &
    -k(178)*n(idx_H3Oj)  &
    -k(175)*n(idx_H3Oj)  &
    -k(177)*n(idx_H3Oj)
pdj(35) =  &
    -k(179)*n(idx_O2j)
pdj(36) =  &
    +k(7)*n(idx_HEj)  &
    -k(15)*n(idx_HEjj)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(1)*1d-3
if(dnn>0.d0) then
nn(1) = n(1) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==2) then
pdj(1) =  &
    +k(219)  &
    +k(29)*n(idx_Hj)  &
    +k(184)*n(idx_C)  &
    +k(26)*n(idx_H)  &
    +k(27)*n(idx_H)  &
    -k(25)*n(idx_E)  &
    +k(186)*n(idx_OH)  &
    +k(17)*n(idx_H)  &
    +k(185)*n(idx_O)  &
    +2.d0*k(25)*n(idx_E)
pdj(2) =  &
    -k(27)*n(idx_H)  &
    -k(184)*n(idx_C)  &
    -k(161)*n(idx_HEj)  &
    -k(28)*n(idx_Hj)  &
    -k(25)*n(idx_E)  &
    -k(219)  &
    -k(29)*n(idx_Hj)  &
    -k(17)*n(idx_H)  &
    -k(32)*n(idx_H2j)  &
    -k(186)*n(idx_OH)  &
    -k(26)*n(idx_H)  &
    -k(185)*n(idx_O)
pdj(5) =  &
    -k(27)*n(idx_H)  &
    +k(219)  &
    +2.d0*k(26)*n(idx_H)  &
    +2.d0*k(27)*n(idx_H)  &
    +2.d0*k(28)*n(idx_Hj)  &
    -k(17)*n(idx_H)  &
    +k(161)*n(idx_HEj)  &
    +k(32)*n(idx_H2j)  &
    +k(25)*n(idx_E)  &
    -k(26)*n(idx_H)
pdj(6) =  &
    +k(161)*n(idx_HEj)
pdj(7) =  &
    +k(32)*n(idx_H2j)  &
    +k(17)*n(idx_H)
pdj(8) =  &
    -k(184)*n(idx_C)
pdj(9) =  &
    -k(185)*n(idx_O)
pdj(10) =  &
    +k(185)*n(idx_O)  &
    -k(186)*n(idx_OH)
pdj(12) =  &
    +k(184)*n(idx_C)
pdj(16) =  &
    +k(186)*n(idx_OH)
pdj(20) =  &
    -k(28)*n(idx_Hj)  &
    -k(29)*n(idx_Hj)
pdj(21) =  &
    -k(161)*n(idx_HEj)
pdj(22) =  &
    -k(32)*n(idx_H2j)  &
    +k(29)*n(idx_Hj)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(2)*1d-3
if(dnn>0.d0) then
nn(2) = n(2) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==3) then
pdj(1) =  &
    +k(189)*n(idx_O)  &
    +k(187)*n(idx_H)  &
    +k(235)  &
    +k(188)*n(idx_H2)
pdj(3) =  &
    -k(188)*n(idx_H2)  &
    -k(159)*n(idx_Hj)  &
    -k(189)*n(idx_O)  &
    -k(187)*n(idx_H)  &
    -k(235)
pdj(5) =  &
    -k(187)*n(idx_H)  &
    +k(159)*n(idx_Hj)
pdj(7) =  &
    -k(188)*n(idx_H2)
pdj(8) =  &
    +k(235)  &
    +k(159)*n(idx_Hj)
pdj(9) =  &
    -k(189)*n(idx_O)
pdj(11) =  &
    +k(189)*n(idx_O)
pdj(12) =  &
    +k(187)*n(idx_H)
pdj(13) =  &
    +k(188)*n(idx_H2)
pdj(20) =  &
    -k(159)*n(idx_Hj)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(3)*1d-3
if(dnn>0.d0) then
nn(3) = n(3) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==4) then
pdj(1) =  &
    +k(191)*n(idx_H2)  &
    +k(242)  &
    +k(192)*n(idx_C)  &
    +k(190)*n(idx_H)
pdj(4) =  &
    -k(242)  &
    -k(192)*n(idx_C)  &
    -k(190)*n(idx_H)  &
    -k(191)*n(idx_H2)  &
    -k(160)*n(idx_Hj)
pdj(5) =  &
    -k(190)*n(idx_H)  &
    +k(160)*n(idx_Hj)
pdj(7) =  &
    -k(191)*n(idx_H2)
pdj(8) =  &
    -k(192)*n(idx_C)
pdj(9) =  &
    +k(242)  &
    +k(160)*n(idx_Hj)
pdj(10) =  &
    +k(190)*n(idx_H)
pdj(11) =  &
    +k(192)*n(idx_C)
pdj(16) =  &
    +k(191)*n(idx_H2)
pdj(20) =  &
    -k(160)*n(idx_Hj)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(4)*1d-3
if(dnn>0.d0) then
nn(4) = n(4) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==5) then
pdj(1) =  &
    +k(26)*n(idx_Hk)  &
    -k(1)*n(idx_E)  &
    +k(27)*n(idx_Hk)  &
    +k(190)*n(idx_Ok)  &
    +2.d0*k(1)*n(idx_E)  &
    -k(16)*n(idx_E)  &
    +k(252)  &
    +k(187)*n(idx_Ck)  &
    +k(213)  &
    +k(17)*n(idx_Hk)
pdj(2) =  &
    -k(26)*n(idx_Hk)  &
    +k(16)*n(idx_E)  &
    -k(27)*n(idx_Hk)  &
    -k(17)*n(idx_Hk)
pdj(3) =  &
    -k(187)*n(idx_Ck)
pdj(4) =  &
    -k(190)*n(idx_Ok)
pdj(5) =  &
    -k(57)*n(idx_CH)  &
    -k(200)*n(idx_Cj)  &
    -k(94)*n(idx_CH2j)  &
    -k(26)*n(idx_Hk)  &
    +2.d0*k(52)*n(idx_OH)  &
    -k(205)*n(idx_O)  &
    -k(17)*n(idx_Hk)  &
    -k(187)*n(idx_Ck)  &
    -9.d0*k(35)*n(idx_H)*n(idx_H)  &
    +3.d0*k(35)*n(idx_H)*n(idx_H)  &
    -k(27)*n(idx_Hk)  &
    -k(213)  &
    -k(280)*n(idx_OH)  &
    -4.d0*k(36)*n(idx_H2)*n(idx_H)  &
    -k(207)*n(idx_OH)  &
    -k(44)*n(idx_Oj)  &
    -k(18)*n(idx_Hj)  &
    +2.d0*k(27)*n(idx_Hk)  &
    -k(158)*n(idx_COj)  &
    -k(24)*n(idx_H2)  &
    -k(8)*n(idx_HEj)  &
    -k(84)*n(idx_CO)  &
    -k(19)*n(idx_Hj)  &
    -k(48)*n(idx_Cj)  &
    +2.d0*k(26)*n(idx_Hk)  &
    -k(63)*n(idx_CH2)  &
    -k(1)*n(idx_E)  &
    -k(252)  &
    -k(52)*n(idx_OH)  &
    +3.d0*k(24)*n(idx_H2)  &
    -4.d0*k(34)*n(idx_H)*n(idx_HE)  &
    -k(279)*n(idx_O)  &
    -k(79)*n(idx_H2O)  &
    -k(97)*n(idx_CH3j)  &
    -k(72)*n(idx_OH)  &
    -k(71)*n(idx_OH)  &
    -k(91)*n(idx_CHj)  &
    -k(196)*n(idx_C)  &
    -k(86)*n(idx_H3j)  &
    -k(80)*n(idx_O2)  &
    -k(20)*n(idx_H2j)  &
    -k(16)*n(idx_E)  &
    -k(190)*n(idx_Ok)
pdj(6) =  &
    +2.d0*k(34)*n(idx_H)*n(idx_HE)  &
    -2.d0*k(34)*n(idx_H)*n(idx_HE)  &
    +k(8)*n(idx_HEj)
pdj(7) =  &
    -2.d0*k(36)*n(idx_H2)*n(idx_H)  &
    +k(63)*n(idx_CH2)  &
    +k(94)*n(idx_CH2j)  &
    +k(57)*n(idx_CH)  &
    +k(20)*n(idx_H2j)  &
    -k(24)*n(idx_H2)  &
    +k(97)*n(idx_CH3j)  &
    +k(72)*n(idx_OH)  &
    +k(71)*n(idx_OH)  &
    +k(86)*n(idx_H3j)  &
    +3.d0*k(35)*n(idx_H)*n(idx_H)  &
    +2.d0*k(34)*n(idx_H)*n(idx_HE)  &
    +4.d0*k(36)*n(idx_H2)*n(idx_H)  &
    +k(91)*n(idx_CHj)  &
    +k(79)*n(idx_H2O)  &
    +k(17)*n(idx_Hk)
pdj(8) =  &
    +k(84)*n(idx_CO)  &
    +k(48)*n(idx_Cj)  &
    +k(57)*n(idx_CH)  &
    -k(196)*n(idx_C)
pdj(9) =  &
    +k(44)*n(idx_Oj)  &
    +k(52)*n(idx_OH)  &
    -k(205)*n(idx_O)  &
    +k(80)*n(idx_O2)  &
    +k(72)*n(idx_OH)  &
    +k(71)*n(idx_OH)  &
    -k(279)*n(idx_O)
pdj(10) =  &
    +k(84)*n(idx_CO)  &
    -k(280)*n(idx_OH)  &
    -k(72)*n(idx_OH)  &
    -k(71)*n(idx_OH)  &
    +k(190)*n(idx_Ok)  &
    -k(52)*n(idx_OH)  &
    -k(207)*n(idx_OH)  &
    +k(80)*n(idx_O2)  &
    +k(205)*n(idx_O)  &
    +k(79)*n(idx_H2O)  &
    +k(279)*n(idx_O)
pdj(11) =  &
    +k(158)*n(idx_COj)  &
    -k(84)*n(idx_CO)
pdj(12) =  &
    +k(196)*n(idx_C)  &
    -k(57)*n(idx_CH)  &
    +k(187)*n(idx_Ck)  &
    +k(63)*n(idx_CH2)
pdj(13) =  &
    -k(63)*n(idx_CH2)
pdj(16) =  &
    -k(79)*n(idx_H2O)  &
    +k(280)*n(idx_OH)  &
    +k(207)*n(idx_OH)
pdj(17) =  &
    -k(80)*n(idx_O2)
pdj(20) =  &
    +k(158)*n(idx_COj)  &
    +k(44)*n(idx_Oj)  &
    +k(48)*n(idx_Cj)  &
    -k(18)*n(idx_Hj)  &
    +k(20)*n(idx_H2j)  &
    +k(1)*n(idx_E)  &
    -k(19)*n(idx_Hj)  &
    +k(252)  &
    +k(8)*n(idx_HEj)  &
    +k(213)
pdj(21) =  &
    -k(8)*n(idx_HEj)
pdj(22) =  &
    +k(18)*n(idx_Hj)  &
    +k(86)*n(idx_H3j)  &
    -k(20)*n(idx_H2j)  &
    +k(19)*n(idx_Hj)
pdj(23) =  &
    -k(48)*n(idx_Cj)  &
    +k(91)*n(idx_CHj)  &
    -k(200)*n(idx_Cj)
pdj(24) =  &
    -k(44)*n(idx_Oj)
pdj(27) =  &
    -k(86)*n(idx_H3j)
pdj(28) =  &
    +k(200)*n(idx_Cj)  &
    +k(94)*n(idx_CH2j)  &
    -k(91)*n(idx_CHj)
pdj(29) =  &
    -k(94)*n(idx_CH2j)  &
    +k(97)*n(idx_CH3j)
pdj(30) =  &
    -k(158)*n(idx_COj)
pdj(31) =  &
    -k(97)*n(idx_CH3j)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(5)*1d-3
if(dnn>0.d0) then
nn(5) = n(5) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==6) then
pdj(1) =  &
    +k(253)  &
    +k(214)  &
    +2.d0*k(4)*n(idx_E)  &
    -k(4)*n(idx_E)
pdj(5) =  &
    +k(9)*n(idx_Hj)  &
    +k(10)*n(idx_Hj)  &
    +2.d0*k(11)*n(idx_H2)  &
    -2.d0*k(34)*n(idx_H)*n(idx_H)
pdj(6) =  &
    -k(10)*n(idx_Hj)  &
    -k(11)*n(idx_H2)  &
    -k(214)  &
    -k(9)*n(idx_Hj)  &
    -k(253)  &
    -k(4)*n(idx_E)  &
    +k(11)*n(idx_H2)  &
    -k(34)*n(idx_H)*n(idx_H)  &
    +k(34)*n(idx_H)*n(idx_H)
pdj(7) =  &
    -k(11)*n(idx_H2)  &
    +k(34)*n(idx_H)*n(idx_H)
pdj(20) =  &
    -k(9)*n(idx_Hj)  &
    -k(10)*n(idx_Hj)
pdj(21) =  &
    +k(9)*n(idx_Hj)  &
    +k(10)*n(idx_Hj)  &
    +k(253)  &
    +k(214)  &
    +k(4)*n(idx_E)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(6)*1d-3
if(dnn>0.d0) then
nn(6) = n(6) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==7) then
pdj(1) =  &
    +k(23)*n(idx_E)  &
    +k(270)  &
    +k(191)*n(idx_Ok)  &
    -k(23)*n(idx_E)  &
    +k(260)  &
    +k(229)  &
    +k(188)*n(idx_Ck)  &
    +k(218)
pdj(2) =  &
    +k(259)
pdj(3) =  &
    -k(188)*n(idx_Ck)
pdj(4) =  &
    -k(191)*n(idx_Ok)
pdj(5) =  &
    +2.d0*k(23)*n(idx_E)  &
    +2.d0*k(258)  &
    +k(110)*n(idx_H2Oj)  &
    +k(90)*n(idx_Cj)  &
    +k(95)*n(idx_CH2j)  &
    +4.d0*k(33)*n(idx_H2)  &
    +2.d0*k(11)*n(idx_HE)  &
    +k(101)*n(idx_Oj)  &
    +k(22)*n(idx_Hj)  &
    +k(21)*n(idx_Hj)  &
    +k(92)*n(idx_CHj)  &
    +3.d0*k(24)*n(idx_H)  &
    +k(70)*n(idx_O)  &
    +k(109)*n(idx_OHj)  &
    +k(13)*n(idx_HEj)  &
    +k(85)*n(idx_H2j)  &
    +k(73)*n(idx_OH)  &
    +k(58)*n(idx_CH)  &
    +2.d0*k(14)*n(idx_HEj)  &
    +k(56)*n(idx_C)  &
    -2.d0*k(36)*n(idx_H)*n(idx_H)  &
    -k(24)*n(idx_H)  &
    +k(270)  &
    +2.d0*k(231)  &
    +k(229)  &
    +2.d0*k(193)*n(idx_Hj)
pdj(6) =  &
    +k(12)*n(idx_HEj)  &
    -k(11)*n(idx_HE)  &
    +k(13)*n(idx_HEj)  &
    +k(11)*n(idx_HE)
pdj(7) =  &
    -k(13)*n(idx_HEj)  &
    -k(270)  &
    -k(22)*n(idx_Hj)  &
    -k(260)  &
    +2.d0*k(33)*n(idx_H2)  &
    -k(218)  &
    -k(197)*n(idx_C)  &
    +2.d0*k(36)*n(idx_H)*n(idx_H)  &
    -k(70)*n(idx_O)  &
    -k(194)*n(idx_Hj)  &
    -k(259)  &
    -k(110)*n(idx_H2Oj)  &
    -k(188)*n(idx_Ck)  &
    -k(101)*n(idx_Oj)  &
    -k(191)*n(idx_Ok)  &
    -k(95)*n(idx_CH2j)  &
    -k(193)*n(idx_Hj)  &
    -k(201)*n(idx_Cj)  &
    -k(92)*n(idx_CHj)  &
    -k(11)*n(idx_HE)  &
    -k(14)*n(idx_HEj)  &
    -k(258)  &
    -k(12)*n(idx_HEj)  &
    -k(58)*n(idx_CH)  &
    -k(85)*n(idx_H2j)  &
    -k(53)*n(idx_HOCj)  &
    -k(36)*n(idx_H)*n(idx_H)  &
    -k(24)*n(idx_H)  &
    -k(73)*n(idx_OH)  &
    -k(109)*n(idx_OHj)  &
    -k(231)  &
    -k(90)*n(idx_Cj)  &
    -k(56)*n(idx_C)  &
    -k(23)*n(idx_E)  &
    -4.d0*k(33)*n(idx_H2)  &
    -k(229)  &
    -k(21)*n(idx_Hj)  &
    +k(53)*n(idx_HOCj)  &
    -k(81)*n(idx_O2)
pdj(8) =  &
    -k(197)*n(idx_C)  &
    -k(56)*n(idx_C)
pdj(9) =  &
    -k(70)*n(idx_O)
pdj(10) =  &
    +2.d0*k(81)*n(idx_O2)  &
    -k(73)*n(idx_OH)  &
    +k(70)*n(idx_O)
pdj(12) =  &
    -k(58)*n(idx_CH)  &
    +k(56)*n(idx_C)
pdj(13) =  &
    +k(58)*n(idx_CH)  &
    +k(197)*n(idx_C)  &
    +k(188)*n(idx_Ck)
pdj(16) =  &
    +k(191)*n(idx_Ok)  &
    +k(73)*n(idx_OH)
pdj(17) =  &
    -k(81)*n(idx_O2)
pdj(20) =  &
    +k(270)  &
    -k(22)*n(idx_Hj)  &
    -k(193)*n(idx_Hj)  &
    +k(259)  &
    +k(229)  &
    +k(193)*n(idx_Hj)  &
    -k(21)*n(idx_Hj)  &
    -k(194)*n(idx_Hj)  &
    +k(13)*n(idx_HEj)
pdj(21) =  &
    -k(13)*n(idx_HEj)  &
    -k(12)*n(idx_HEj)  &
    +k(14)*n(idx_HEj)  &
    -k(14)*n(idx_HEj)
pdj(22) =  &
    +k(22)*n(idx_Hj)  &
    +k(260)  &
    +k(21)*n(idx_Hj)  &
    -k(85)*n(idx_H2j)  &
    +k(12)*n(idx_HEj)  &
    +k(218)
pdj(23) =  &
    -k(90)*n(idx_Cj)  &
    -k(201)*n(idx_Cj)
pdj(24) =  &
    -k(101)*n(idx_Oj)
pdj(25) =  &
    -k(53)*n(idx_HOCj)
pdj(26) =  &
    +k(53)*n(idx_HOCj)
pdj(27) =  &
    +k(194)*n(idx_Hj)  &
    +k(85)*n(idx_H2j)
pdj(28) =  &
    -k(92)*n(idx_CHj)  &
    +k(90)*n(idx_Cj)
pdj(29) =  &
    -k(95)*n(idx_CH2j)  &
    +k(92)*n(idx_CHj)  &
    +k(201)*n(idx_Cj)
pdj(31) =  &
    +k(95)*n(idx_CH2j)
pdj(32) =  &
    +k(101)*n(idx_Oj)  &
    -k(109)*n(idx_OHj)
pdj(33) =  &
    +k(109)*n(idx_OHj)  &
    -k(110)*n(idx_H2Oj)
pdj(34) =  &
    +k(110)*n(idx_H2Oj)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(7)*1d-3
if(dnn>0.d0) then
nn(7) = n(7) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==8) then
pdj(1) =  &
    -k(195)*n(idx_E)  &
    +k(261)  &
    +k(192)*n(idx_Ok)  &
    +2.d0*k(42)*n(idx_E)  &
    +k(217)  &
    -k(42)*n(idx_E)  &
    +k(184)*n(idx_Hk)
pdj(2) =  &
    -k(184)*n(idx_Hk)
pdj(3) =  &
    +k(195)*n(idx_E)
pdj(4) =  &
    -k(192)*n(idx_Ok)
pdj(5) =  &
    +k(74)*n(idx_OH)  &
    +k(87)*n(idx_H2j)  &
    +k(89)*n(idx_H3j)  &
    -k(196)*n(idx_H)  &
    +k(75)*n(idx_OH)  &
    +k(59)*n(idx_CH)  &
    +k(47)*n(idx_Hj)  &
    +k(56)*n(idx_H2)
pdj(6) =  &
    +k(51)*n(idx_HEj)  &
    +k(50)*n(idx_HEj)  &
    +k(49)*n(idx_HEj)
pdj(7) =  &
    -k(197)*n(idx_H2)  &
    +k(88)*n(idx_H3j)  &
    -k(56)*n(idx_H2)  &
    +k(117)*n(idx_H3Oj)
pdj(8) =  &
    -k(59)*n(idx_CH)  &
    -k(87)*n(idx_H2j)  &
    -k(278)*n(idx_Oj)  &
    -k(196)*n(idx_H)  &
    -k(197)*n(idx_H2)  &
    -k(51)*n(idx_HEj)  &
    -k(42)*n(idx_E)  &
    -k(117)*n(idx_H3Oj)  &
    -k(89)*n(idx_H3j)  &
    -k(195)*n(idx_E)  &
    -k(274)*n(idx_O)  &
    -k(184)*n(idx_Hk)  &
    -4.d0*k(198)*n(idx_C)  &
    -k(277)*n(idx_Oj)  &
    -k(273)*n(idx_O)  &
    -k(217)  &
    -k(49)*n(idx_HEj)  &
    -k(88)*n(idx_H3j)  &
    -k(75)*n(idx_OH)  &
    -k(83)*n(idx_O2)  &
    -k(261)  &
    -k(127)*n(idx_HCOj)  &
    -k(50)*n(idx_HEj)  &
    -4.d0*k(272)*n(idx_C)  &
    -k(192)*n(idx_Ok)  &
    -4.d0*k(271)*n(idx_C)  &
    -k(199)*n(idx_O)  &
    -k(47)*n(idx_Hj)  &
    -k(122)*n(idx_O2j)  &
    -k(121)*n(idx_O2j)  &
    -k(74)*n(idx_OH)  &
    -k(56)*n(idx_H2)  &
    -k(82)*n(idx_O2)
pdj(9) =  &
    -k(274)*n(idx_O)  &
    -k(199)*n(idx_O)  &
    +k(83)*n(idx_O2)  &
    +k(121)*n(idx_O2j)  &
    +k(82)*n(idx_O2)  &
    -k(273)*n(idx_O)
pdj(10) =  &
    -k(74)*n(idx_OH)  &
    -k(75)*n(idx_OH)
pdj(11) =  &
    +k(74)*n(idx_OH)  &
    +k(192)*n(idx_Ok)  &
    +k(82)*n(idx_O2)  &
    +k(274)*n(idx_O)  &
    +k(83)*n(idx_O2)  &
    +k(127)*n(idx_HCOj)  &
    +k(75)*n(idx_OH)  &
    +k(199)*n(idx_O)  &
    +k(273)*n(idx_O)
pdj(12) =  &
    +k(196)*n(idx_H)  &
    -k(59)*n(idx_CH)  &
    +k(184)*n(idx_Hk)  &
    +k(56)*n(idx_H2)
pdj(13) =  &
    +k(197)*n(idx_H2)
pdj(14) =  &
    +2.d0*k(272)*n(idx_C)  &
    +k(59)*n(idx_CH)  &
    +2.d0*k(198)*n(idx_C)  &
    +2.d0*k(271)*n(idx_C)
pdj(17) =  &
    -k(83)*n(idx_O2)  &
    -k(82)*n(idx_O2)  &
    +k(122)*n(idx_O2j)
pdj(20) =  &
    -k(47)*n(idx_Hj)
pdj(21) =  &
    -k(49)*n(idx_HEj)  &
    -k(51)*n(idx_HEj)  &
    -k(50)*n(idx_HEj)
pdj(22) =  &
    -k(87)*n(idx_H2j)
pdj(23) =  &
    +k(261)  &
    +k(51)*n(idx_HEj)  &
    +k(42)*n(idx_E)  &
    +k(49)*n(idx_HEj)  &
    +k(50)*n(idx_HEj)  &
    +k(217)  &
    +k(122)*n(idx_O2j)  &
    +k(47)*n(idx_Hj)
pdj(24) =  &
    -k(278)*n(idx_Oj)  &
    -k(277)*n(idx_Oj)
pdj(26) =  &
    -k(127)*n(idx_HCOj)  &
    +k(117)*n(idx_H3Oj)
pdj(27) =  &
    -k(88)*n(idx_H3j)  &
    -k(89)*n(idx_H3j)
pdj(28) =  &
    +k(127)*n(idx_HCOj)  &
    +k(88)*n(idx_H3j)  &
    +k(87)*n(idx_H2j)
pdj(29) =  &
    +k(89)*n(idx_H3j)
pdj(30) =  &
    +k(121)*n(idx_O2j)  &
    +k(277)*n(idx_Oj)  &
    +k(278)*n(idx_Oj)
pdj(34) =  &
    -k(117)*n(idx_H3Oj)
pdj(35) =  &
    -k(121)*n(idx_O2j)  &
    -k(122)*n(idx_O2j)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(8)*1d-3
if(dnn>0.d0) then
nn(8) = n(8) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==9) then
pdj(1) =  &
    -k(43)*n(idx_E)  &
    -k(204)*n(idx_E)  &
    +k(189)*n(idx_Ck)  &
    +2.d0*k(43)*n(idx_E)  &
    +k(216)  &
    +k(185)*n(idx_Hk)  &
    +k(61)*n(idx_CH)  &
    +k(254)
pdj(2) =  &
    -k(185)*n(idx_Hk)
pdj(3) =  &
    -k(189)*n(idx_Ck)
pdj(4) =  &
    +k(204)*n(idx_E)
pdj(5) =  &
    +k(70)*n(idx_H2)  &
    +k(104)*n(idx_H3j)  &
    +k(76)*n(idx_OH)  &
    +k(60)*n(idx_CH)  &
    +k(96)*n(idx_CH2j)  &
    +k(93)*n(idx_CHj)  &
    +k(66)*n(idx_CH2)  &
    +k(102)*n(idx_H2j)  &
    -k(279)*n(idx_H)  &
    -k(205)*n(idx_H)  &
    +k(45)*n(idx_Hj)  &
    +2.d0*k(64)*n(idx_CH2)  &
    +k(77)*n(idx_OH)
pdj(6) =  &
    +k(46)*n(idx_HEj)
pdj(7) =  &
    -k(70)*n(idx_H2)  &
    +k(98)*n(idx_CH3j)  &
    +k(99)*n(idx_CH3j)  &
    +k(65)*n(idx_CH2)  &
    +k(103)*n(idx_H3j)
pdj(8) =  &
    -k(199)*n(idx_C)  &
    -k(274)*n(idx_C)  &
    +k(62)*n(idx_CH)  &
    +k(69)*n(idx_C2)  &
    +k(68)*n(idx_C2)  &
    -k(273)*n(idx_C)
pdj(9) =  &
    -k(199)*n(idx_C)  &
    -k(43)*n(idx_E)  &
    -k(216)  &
    -k(46)*n(idx_HEj)  &
    -k(202)*n(idx_Cj)  &
    -4.d0*k(281)*n(idx_O)  &
    -k(99)*n(idx_CH3j)  &
    -k(102)*n(idx_H2j)  &
    -k(77)*n(idx_OH)  &
    -k(205)*n(idx_H)  &
    -k(279)*n(idx_H)  &
    -k(69)*n(idx_C2)  &
    -k(65)*n(idx_CH2)  &
    -k(103)*n(idx_H3j)  &
    -k(70)*n(idx_H2)  &
    -k(45)*n(idx_Hj)  &
    -k(276)*n(idx_Cj)  &
    -k(185)*n(idx_Hk)  &
    -k(96)*n(idx_CH2j)  &
    -k(273)*n(idx_C)  &
    -k(203)*n(idx_Cj)  &
    -k(274)*n(idx_C)  &
    -k(60)*n(idx_CH)  &
    -k(62)*n(idx_CH)  &
    -k(66)*n(idx_CH2)  &
    -k(61)*n(idx_CH)  &
    -k(189)*n(idx_Ck)  &
    -k(98)*n(idx_CH3j)  &
    -k(67)*n(idx_CH2)  &
    -4.d0*k(206)*n(idx_O)  &
    -k(104)*n(idx_H3j)  &
    -k(275)*n(idx_Cj)  &
    -k(68)*n(idx_C2)  &
    -k(93)*n(idx_CHj)  &
    -k(204)*n(idx_E)  &
    -k(76)*n(idx_OH)  &
    -k(64)*n(idx_CH2)  &
    -k(254)
pdj(10) =  &
    +k(62)*n(idx_CH)  &
    +k(279)*n(idx_H)  &
    +k(70)*n(idx_H2)  &
    -k(76)*n(idx_OH)  &
    -k(77)*n(idx_OH)  &
    +k(185)*n(idx_Hk)  &
    +k(67)*n(idx_CH2)  &
    +k(205)*n(idx_H)
pdj(11) =  &
    +k(274)*n(idx_C)  &
    +k(60)*n(idx_CH)  &
    +k(189)*n(idx_Ck)  &
    +k(69)*n(idx_C2)  &
    +k(65)*n(idx_CH2)  &
    +k(273)*n(idx_C)  &
    +k(64)*n(idx_CH2)  &
    +k(199)*n(idx_C)  &
    +k(68)*n(idx_C2)
pdj(12) =  &
    -k(60)*n(idx_CH)  &
    -k(62)*n(idx_CH)  &
    -k(61)*n(idx_CH)  &
    +k(67)*n(idx_CH2)
pdj(13) =  &
    -k(64)*n(idx_CH2)  &
    -k(65)*n(idx_CH2)  &
    -k(66)*n(idx_CH2)  &
    -k(67)*n(idx_CH2)
pdj(14) =  &
    -k(69)*n(idx_C2)  &
    -k(68)*n(idx_C2)
pdj(15) =  &
    +k(66)*n(idx_CH2)
pdj(17) =  &
    +k(76)*n(idx_OH)  &
    +2.d0*k(281)*n(idx_O)  &
    +k(77)*n(idx_OH)  &
    +2.d0*k(206)*n(idx_O)
pdj(20) =  &
    -k(45)*n(idx_Hj)
pdj(21) =  &
    -k(46)*n(idx_HEj)
pdj(22) =  &
    -k(102)*n(idx_H2j)
pdj(23) =  &
    -k(275)*n(idx_Cj)  &
    -k(203)*n(idx_Cj)  &
    -k(202)*n(idx_Cj)  &
    -k(276)*n(idx_Cj)
pdj(24) =  &
    +k(46)*n(idx_HEj)  &
    +k(216)  &
    +k(45)*n(idx_Hj)  &
    +k(254)  &
    +k(43)*n(idx_E)
pdj(25) =  &
    +k(98)*n(idx_CH3j)
pdj(26) =  &
    +k(61)*n(idx_CH)  &
    +k(99)*n(idx_CH3j)  &
    +k(96)*n(idx_CH2j)
pdj(27) =  &
    -k(104)*n(idx_H3j)  &
    -k(103)*n(idx_H3j)
pdj(28) =  &
    -k(93)*n(idx_CHj)
pdj(29) =  &
    -k(96)*n(idx_CH2j)
pdj(30) =  &
    +k(275)*n(idx_Cj)  &
    +k(203)*n(idx_Cj)  &
    +k(93)*n(idx_CHj)  &
    +k(202)*n(idx_Cj)  &
    +k(276)*n(idx_Cj)
pdj(31) =  &
    -k(98)*n(idx_CH3j)  &
    -k(99)*n(idx_CH3j)
pdj(32) =  &
    +k(103)*n(idx_H3j)  &
    +k(102)*n(idx_H2j)
pdj(33) =  &
    +k(104)*n(idx_H3j)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(9)*1d-3
if(dnn>0.d0) then
nn(9) = n(9) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==10) then
pdj(1) =  &
    +k(224)  &
    +k(186)*n(idx_Hk)
pdj(2) =  &
    -k(186)*n(idx_Hk)
pdj(5) =  &
    -k(280)*n(idx_H)  &
    -k(207)*n(idx_H)  &
    +k(144)*n(idx_HEj)  &
    +k(223)  &
    +k(77)*n(idx_O)  &
    +k(141)*n(idx_Hj)  &
    +k(143)*n(idx_HEj)  &
    +k(142)*n(idx_Hj)  &
    -k(72)*n(idx_H)  &
    +k(76)*n(idx_O)  &
    -k(71)*n(idx_H)  &
    +k(107)*n(idx_Cj)  &
    +k(75)*n(idx_C)  &
    +k(108)*n(idx_Cj)  &
    +2.d0*k(52)*n(idx_H)  &
    +k(265)  &
    -k(52)*n(idx_H)  &
    +k(73)*n(idx_H2)  &
    +k(74)*n(idx_C)
pdj(6) =  &
    +k(143)*n(idx_HEj)  &
    +k(144)*n(idx_HEj)
pdj(7) =  &
    +k(72)*n(idx_H)  &
    +k(106)*n(idx_H3j)  &
    +k(105)*n(idx_H3j)  &
    -k(73)*n(idx_H2)  &
    +k(71)*n(idx_H)
pdj(8) =  &
    -k(75)*n(idx_C)  &
    -k(74)*n(idx_C)
pdj(9) =  &
    -k(77)*n(idx_O)  &
    +k(71)*n(idx_H)  &
    +k(223)  &
    -k(76)*n(idx_O)  &
    +k(72)*n(idx_H)  &
    +k(52)*n(idx_H)  &
    +k(265)  &
    +2.d0*k(78)*n(idx_OH)
pdj(10) =  &
    -k(280)*n(idx_H)  &
    -k(143)*n(idx_HEj)  &
    -k(77)*n(idx_O)  &
    -k(72)*n(idx_H)  &
    -k(142)*n(idx_Hj)  &
    -k(186)*n(idx_Hk)  &
    -k(224)  &
    -k(141)*n(idx_Hj)  &
    -k(76)*n(idx_O)  &
    -k(75)*n(idx_C)  &
    -4.d0*k(78)*n(idx_OH)  &
    -k(144)*n(idx_HEj)  &
    -k(73)*n(idx_H2)  &
    -k(106)*n(idx_H3j)  &
    -k(52)*n(idx_H)  &
    -k(207)*n(idx_H)  &
    -k(74)*n(idx_C)  &
    -k(223)  &
    -k(107)*n(idx_Cj)  &
    -k(108)*n(idx_Cj)  &
    -k(71)*n(idx_H)  &
    -k(265)  &
    -k(105)*n(idx_H3j)
pdj(11) =  &
    +k(75)*n(idx_C)  &
    +k(74)*n(idx_C)
pdj(16) =  &
    +k(73)*n(idx_H2)  &
    +k(207)*n(idx_H)  &
    +2.d0*k(78)*n(idx_OH)  &
    +k(280)*n(idx_H)  &
    +k(186)*n(idx_Hk)
pdj(17) =  &
    +k(77)*n(idx_O)  &
    +k(76)*n(idx_O)
pdj(20) =  &
    -k(142)*n(idx_Hj)  &
    -k(141)*n(idx_Hj)
pdj(21) =  &
    -k(144)*n(idx_HEj)  &
    -k(143)*n(idx_HEj)
pdj(23) =  &
    -k(107)*n(idx_Cj)  &
    -k(108)*n(idx_Cj)
pdj(24) =  &
    +k(143)*n(idx_HEj)  &
    +k(144)*n(idx_HEj)
pdj(27) =  &
    -k(105)*n(idx_H3j)  &
    -k(106)*n(idx_H3j)
pdj(30) =  &
    +k(108)*n(idx_Cj)  &
    +k(107)*n(idx_Cj)
pdj(32) =  &
    +k(224)  &
    +k(141)*n(idx_Hj)  &
    +k(142)*n(idx_Hj)
pdj(33) =  &
    +k(106)*n(idx_H3j)  &
    +k(105)*n(idx_H3j)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(10)*1d-3
if(dnn>0.d0) then
nn(10) = n(10) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==11) then
pdj(1) =  &
    +k(256)
pdj(5) =  &
    -k(84)*n(idx_H)
pdj(6) =  &
    +k(157)*n(idx_HEj)  &
    +k(156)*n(idx_HEj)
pdj(7) =  &
    +k(123)*n(idx_H3j)  &
    +k(124)*n(idx_H3j)  &
    +k(125)*n(idx_H3j)  &
    +k(126)*n(idx_H3j)
pdj(8) =  &
    +k(255)  &
    +k(84)*n(idx_H)  &
    +k(157)*n(idx_HEj)  &
    +k(230)
pdj(9) =  &
    +k(255)  &
    +k(156)*n(idx_HEj)  &
    +k(230)
pdj(10) =  &
    +k(84)*n(idx_H)
pdj(11) =  &
    +k(55)*n(idx_HOCj)  &
    -k(125)*n(idx_H3j)  &
    -k(157)*n(idx_HEj)  &
    -k(54)*n(idx_HOCj)  &
    -k(55)*n(idx_HOCj)  &
    -k(84)*n(idx_H)  &
    -k(230)  &
    -k(123)*n(idx_H3j)  &
    -k(156)*n(idx_HEj)  &
    +k(54)*n(idx_HOCj)  &
    -k(124)*n(idx_H3j)  &
    -k(256)  &
    -k(126)*n(idx_H3j)  &
    -k(255)
pdj(21) =  &
    -k(157)*n(idx_HEj)  &
    -k(156)*n(idx_HEj)
pdj(23) =  &
    +k(156)*n(idx_HEj)
pdj(24) =  &
    +k(157)*n(idx_HEj)
pdj(25) =  &
    -k(55)*n(idx_HOCj)  &
    +k(125)*n(idx_H3j)  &
    +k(126)*n(idx_H3j)  &
    -k(54)*n(idx_HOCj)
pdj(26) =  &
    +k(54)*n(idx_HOCj)  &
    +k(55)*n(idx_HOCj)  &
    +k(123)*n(idx_H3j)  &
    +k(124)*n(idx_H3j)
pdj(27) =  &
    -k(126)*n(idx_H3j)  &
    -k(124)*n(idx_H3j)  &
    -k(125)*n(idx_H3j)  &
    -k(123)*n(idx_H3j)
pdj(30) =  &
    +k(256)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(11)*1d-3
if(dnn>0.d0) then
nn(11) = n(11) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==12) then
pdj(1) =  &
    +k(61)*n(idx_O)  &
    +k(221)
pdj(5) =  &
    +k(131)*n(idx_Hj)  &
    +k(60)*n(idx_O)  &
    +k(220)  &
    -k(57)*n(idx_H)  &
    +k(262)  &
    +k(58)*n(idx_H2)  &
    +k(130)*n(idx_Hj)  &
    +k(59)*n(idx_C)
pdj(7) =  &
    +k(57)*n(idx_H)  &
    -k(58)*n(idx_H2)
pdj(8) =  &
    -k(59)*n(idx_C)  &
    +k(57)*n(idx_H)  &
    +k(220)  &
    +k(62)*n(idx_O)  &
    +k(262)
pdj(9) =  &
    -k(60)*n(idx_O)  &
    -k(62)*n(idx_O)  &
    -k(61)*n(idx_O)
pdj(10) =  &
    +k(62)*n(idx_O)
pdj(11) =  &
    +k(60)*n(idx_O)
pdj(12) =  &
    -k(60)*n(idx_O)  &
    -k(262)  &
    -k(131)*n(idx_Hj)  &
    -k(130)*n(idx_Hj)  &
    -k(57)*n(idx_H)  &
    -k(220)  &
    -k(61)*n(idx_O)  &
    -k(58)*n(idx_H2)  &
    -k(59)*n(idx_C)  &
    -k(62)*n(idx_O)  &
    -k(221)
pdj(13) =  &
    +k(58)*n(idx_H2)
pdj(14) =  &
    +k(59)*n(idx_C)
pdj(20) =  &
    -k(130)*n(idx_Hj)  &
    -k(131)*n(idx_Hj)
pdj(26) =  &
    +k(61)*n(idx_O)
pdj(28) =  &
    +k(130)*n(idx_Hj)  &
    +k(131)*n(idx_Hj)  &
    +k(221)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(12)*1d-3
if(dnn>0.d0) then
nn(12) = n(12) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==13) then
pdj(1) =  &
    +k(266)  &
    +k(238)
pdj(5) =  &
    +2.d0*k(64)*n(idx_O)  &
    +k(138)*n(idx_HEj)  &
    +k(66)*n(idx_O)  &
    -k(63)*n(idx_H)  &
    +k(139)*n(idx_HEj)  &
    +k(237)  &
    +k(134)*n(idx_Hj)  &
    +k(135)*n(idx_Hj)
pdj(6) =  &
    +k(139)*n(idx_HEj)  &
    +k(138)*n(idx_HEj)  &
    +k(137)*n(idx_HEj)  &
    +k(136)*n(idx_HEj)
pdj(7) =  &
    +k(133)*n(idx_Hj)  &
    +k(137)*n(idx_HEj)  &
    +k(136)*n(idx_HEj)  &
    +k(63)*n(idx_H)  &
    +k(132)*n(idx_Hj)  &
    +k(65)*n(idx_O)
pdj(9) =  &
    -k(66)*n(idx_O)  &
    -k(65)*n(idx_O)  &
    -k(67)*n(idx_O)  &
    -k(64)*n(idx_O)
pdj(10) =  &
    +k(67)*n(idx_O)
pdj(11) =  &
    +k(64)*n(idx_O)  &
    +k(65)*n(idx_O)
pdj(12) =  &
    +k(63)*n(idx_H)  &
    +k(67)*n(idx_O)  &
    +k(237)
pdj(13) =  &
    -k(66)*n(idx_O)  &
    -k(266)  &
    -k(237)  &
    -k(137)*n(idx_HEj)  &
    -k(136)*n(idx_HEj)  &
    -k(64)*n(idx_O)  &
    -k(134)*n(idx_Hj)  &
    -k(132)*n(idx_Hj)  &
    -k(67)*n(idx_O)  &
    -k(238)  &
    -k(63)*n(idx_H)  &
    -k(65)*n(idx_O)  &
    -k(135)*n(idx_Hj)  &
    -k(133)*n(idx_Hj)  &
    -k(139)*n(idx_HEj)  &
    -k(138)*n(idx_HEj)
pdj(15) =  &
    +k(66)*n(idx_O)
pdj(20) =  &
    -k(135)*n(idx_Hj)  &
    -k(132)*n(idx_Hj)  &
    -k(133)*n(idx_Hj)  &
    -k(134)*n(idx_Hj)
pdj(21) =  &
    -k(139)*n(idx_HEj)  &
    -k(138)*n(idx_HEj)  &
    -k(137)*n(idx_HEj)  &
    -k(136)*n(idx_HEj)
pdj(23) =  &
    +k(137)*n(idx_HEj)  &
    +k(136)*n(idx_HEj)
pdj(28) =  &
    +k(139)*n(idx_HEj)  &
    +k(138)*n(idx_HEj)  &
    +k(132)*n(idx_Hj)  &
    +k(133)*n(idx_Hj)
pdj(29) =  &
    +k(266)  &
    +k(134)*n(idx_Hj)  &
    +k(135)*n(idx_Hj)  &
    +k(238)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(13)*1d-3
if(dnn>0.d0) then
nn(13) = n(13) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==14) then
pdj(6) =  &
    +k(140)*n(idx_HEj)
pdj(8) =  &
    +2.d0*k(257)  &
    +2.d0*k(222)  &
    +k(140)*n(idx_HEj)  &
    +k(69)*n(idx_O)  &
    +k(68)*n(idx_O)  &
    +k(100)*n(idx_Oj)
pdj(9) =  &
    -k(68)*n(idx_O)  &
    -k(69)*n(idx_O)
pdj(11) =  &
    +k(69)*n(idx_O)  &
    +k(68)*n(idx_O)
pdj(14) =  &
    -k(222)  &
    -k(68)*n(idx_O)  &
    -k(140)*n(idx_HEj)  &
    -k(100)*n(idx_Oj)  &
    -k(69)*n(idx_O)  &
    -k(257)
pdj(21) =  &
    -k(140)*n(idx_HEj)
pdj(23) =  &
    +k(140)*n(idx_HEj)
pdj(24) =  &
    -k(100)*n(idx_Oj)
pdj(30) =  &
    +k(100)*n(idx_Oj)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(14)*1d-3
if(dnn>0.d0) then
nn(14) = n(14) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==15) then
pdj(1) =  &
    +k(269)
pdj(5) =  &
    +k(268)
pdj(11) =  &
    +k(268)
pdj(15) =  &
    -k(269)  &
    -k(268)
pdj(26) =  &
    +k(269)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(15)*1d-3
if(dnn>0.d0) then
nn(15) = n(15) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==16) then
pdj(1) =  &
    +k(226)
pdj(5) =  &
    +k(113)*n(idx_Cj)  &
    +k(115)*n(idx_Cj)  &
    +k(149)*n(idx_HEj)  &
    -k(79)*n(idx_H)  &
    +k(150)*n(idx_HEj)  &
    +k(114)*n(idx_Cj)  &
    +k(145)*n(idx_Hj)  &
    +k(225)  &
    +k(146)*n(idx_Hj)  &
    +k(267)
pdj(6) =  &
    +k(148)*n(idx_HEj)  &
    +k(149)*n(idx_HEj)  &
    +k(150)*n(idx_HEj)  &
    +k(151)*n(idx_HEj)  &
    +k(147)*n(idx_HEj)  &
    +k(152)*n(idx_HEj)
pdj(7) =  &
    +k(111)*n(idx_H3j)  &
    +k(112)*n(idx_H3j)  &
    +k(79)*n(idx_H)
pdj(8) =  &
    +k(116)*n(idx_Cj)
pdj(10) =  &
    +k(148)*n(idx_HEj)  &
    +k(267)  &
    +k(147)*n(idx_HEj)  &
    +k(225)  &
    +k(79)*n(idx_H)
pdj(11) =  &
    +k(128)*n(idx_HCOj)  &
    +k(129)*n(idx_HCOj)
pdj(16) =  &
    -k(129)*n(idx_HCOj)  &
    -k(148)*n(idx_HEj)  &
    -k(116)*n(idx_Cj)  &
    -k(115)*n(idx_Cj)  &
    -k(225)  &
    -k(79)*n(idx_H)  &
    -k(145)*n(idx_Hj)  &
    -k(267)  &
    -k(151)*n(idx_HEj)  &
    -k(150)*n(idx_HEj)  &
    -k(152)*n(idx_HEj)  &
    -k(128)*n(idx_HCOj)  &
    -k(147)*n(idx_HEj)  &
    -k(112)*n(idx_H3j)  &
    -k(111)*n(idx_H3j)  &
    -k(149)*n(idx_HEj)  &
    -k(226)  &
    -k(146)*n(idx_Hj)  &
    -k(114)*n(idx_Cj)  &
    -k(113)*n(idx_Cj)
pdj(20) =  &
    +k(148)*n(idx_HEj)  &
    -k(146)*n(idx_Hj)  &
    -k(145)*n(idx_Hj)  &
    +k(147)*n(idx_HEj)
pdj(21) =  &
    -k(152)*n(idx_HEj)  &
    -k(151)*n(idx_HEj)  &
    -k(150)*n(idx_HEj)  &
    -k(148)*n(idx_HEj)  &
    -k(147)*n(idx_HEj)  &
    -k(149)*n(idx_HEj)
pdj(23) =  &
    -k(116)*n(idx_Cj)  &
    -k(115)*n(idx_Cj)  &
    -k(114)*n(idx_Cj)  &
    -k(113)*n(idx_Cj)
pdj(25) =  &
    +k(113)*n(idx_Cj)
pdj(26) =  &
    -k(129)*n(idx_HCOj)  &
    +k(115)*n(idx_Cj)  &
    +k(114)*n(idx_Cj)  &
    -k(128)*n(idx_HCOj)
pdj(27) =  &
    -k(112)*n(idx_H3j)  &
    -k(111)*n(idx_H3j)
pdj(32) =  &
    +k(149)*n(idx_HEj)  &
    +k(150)*n(idx_HEj)
pdj(33) =  &
    +k(116)*n(idx_Cj)  &
    +k(151)*n(idx_HEj)  &
    +k(145)*n(idx_Hj)  &
    +k(146)*n(idx_Hj)  &
    +k(152)*n(idx_HEj)  &
    +k(226)
pdj(34) =  &
    +k(111)*n(idx_H3j)  &
    +k(128)*n(idx_HCOj)  &
    +k(112)*n(idx_H3j)  &
    +k(129)*n(idx_HCOj)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(16)*1d-3
if(dnn>0.d0) then
nn(16) = n(16) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==17) then
pdj(1) =  &
    +k(264)  &
    +k(227)
pdj(5) =  &
    +k(153)*n(idx_Hj)  &
    -k(80)*n(idx_H)
pdj(6) =  &
    +k(155)*n(idx_HEj)  &
    +k(154)*n(idx_HEj)
pdj(7) =  &
    -k(81)*n(idx_H2)
pdj(8) =  &
    -k(83)*n(idx_C)  &
    -k(82)*n(idx_C)
pdj(9) =  &
    +2.d0*k(263)  &
    +k(155)*n(idx_HEj)  &
    +k(118)*n(idx_Cj)  &
    +2.d0*k(228)  &
    +k(83)*n(idx_C)  &
    +k(82)*n(idx_C)  &
    +k(80)*n(idx_H)
pdj(10) =  &
    +2.d0*k(81)*n(idx_H2)  &
    +k(120)*n(idx_CH2j)  &
    +k(80)*n(idx_H)
pdj(11) =  &
    +k(82)*n(idx_C)  &
    +k(83)*n(idx_C)  &
    +k(119)*n(idx_Cj)
pdj(17) =  &
    -k(118)*n(idx_Cj)  &
    -k(263)  &
    -k(153)*n(idx_Hj)  &
    -k(264)  &
    -k(82)*n(idx_C)  &
    -k(119)*n(idx_Cj)  &
    -k(83)*n(idx_C)  &
    -k(80)*n(idx_H)  &
    -k(154)*n(idx_HEj)  &
    -k(81)*n(idx_H2)  &
    -k(228)  &
    -k(227)  &
    -k(155)*n(idx_HEj)  &
    -k(120)*n(idx_CH2j)
pdj(20) =  &
    -k(153)*n(idx_Hj)
pdj(21) =  &
    -k(155)*n(idx_HEj)  &
    -k(154)*n(idx_HEj)
pdj(23) =  &
    -k(119)*n(idx_Cj)  &
    -k(118)*n(idx_Cj)
pdj(24) =  &
    +k(155)*n(idx_HEj)  &
    +k(119)*n(idx_Cj)
pdj(26) =  &
    +k(120)*n(idx_CH2j)
pdj(29) =  &
    -k(120)*n(idx_CH2j)
pdj(30) =  &
    +k(118)*n(idx_Cj)
pdj(35) =  &
    +k(153)*n(idx_Hj)  &
    +k(264)  &
    +k(227)  &
    +k(154)*n(idx_HEj)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(17)*1d-3
if(dnn>0.d0) then
nn(17) = n(17) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==18) then
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(18)*1d-3
if(dnn>0.d0) then
nn(18) = n(18) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==19) then
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(19)*1d-3
if(dnn>0.d0) then
nn(19) = n(19) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==20) then
pdj(1) =  &
    +k(29)*n(idx_Hk)  &
    -k(3)*n(idx_E)  &
    -k(2)*n(idx_E)
pdj(2) =  &
    -k(29)*n(idx_Hk)  &
    -k(28)*n(idx_Hk)
pdj(3) =  &
    -k(159)*n(idx_Ck)
pdj(4) =  &
    -k(160)*n(idx_Ok)
pdj(5) =  &
    +k(145)*n(idx_H2O)  &
    +k(22)*n(idx_H2)  &
    +k(3)*n(idx_E)  &
    +k(159)*n(idx_Ck)  &
    +k(9)*n(idx_HE)  &
    +k(135)*n(idx_CH2)  &
    +k(153)*n(idx_O2)  &
    +k(21)*n(idx_H2)  &
    -k(19)*n(idx_H)  &
    +k(141)*n(idx_OH)  &
    +k(130)*n(idx_CH)  &
    -k(18)*n(idx_H)  &
    +2.d0*k(28)*n(idx_Hk)  &
    +k(146)*n(idx_H2O)  &
    +k(2)*n(idx_E)  &
    +k(10)*n(idx_HE)  &
    +k(47)*n(idx_C)  &
    +k(160)*n(idx_Ok)  &
    +k(142)*n(idx_OH)  &
    +k(134)*n(idx_CH2)  &
    +k(45)*n(idx_O)  &
    +2.d0*k(193)*n(idx_H2)  &
    +k(131)*n(idx_CH)
pdj(6) =  &
    -k(10)*n(idx_HE)  &
    -k(9)*n(idx_HE)
pdj(7) =  &
    -k(193)*n(idx_H2)  &
    -k(194)*n(idx_H2)  &
    +k(133)*n(idx_CH2)  &
    +k(132)*n(idx_CH2)  &
    -k(22)*n(idx_H2)  &
    -k(21)*n(idx_H2)
pdj(8) =  &
    +k(159)*n(idx_Ck)  &
    -k(47)*n(idx_C)
pdj(9) =  &
    -k(45)*n(idx_O)  &
    +k(160)*n(idx_Ok)
pdj(10) =  &
    -k(141)*n(idx_OH)  &
    -k(142)*n(idx_OH)
pdj(12) =  &
    -k(131)*n(idx_CH)  &
    -k(130)*n(idx_CH)
pdj(13) =  &
    -k(133)*n(idx_CH2)  &
    -k(132)*n(idx_CH2)  &
    -k(135)*n(idx_CH2)  &
    -k(134)*n(idx_CH2)
pdj(16) =  &
    -k(146)*n(idx_H2O)  &
    -k(145)*n(idx_H2O)
pdj(17) =  &
    -k(153)*n(idx_O2)
pdj(20) =  &
    -k(193)*n(idx_H2)  &
    -k(131)*n(idx_CH)  &
    -k(194)*n(idx_H2)  &
    -k(153)*n(idx_O2)  &
    -k(130)*n(idx_CH)  &
    -k(146)*n(idx_H2O)  &
    -k(134)*n(idx_CH2)  &
    -k(2)*n(idx_E)  &
    -k(19)*n(idx_H)  &
    -k(47)*n(idx_C)  &
    -k(145)*n(idx_H2O)  &
    -k(132)*n(idx_CH2)  &
    -k(18)*n(idx_H)  &
    -k(142)*n(idx_OH)  &
    -k(10)*n(idx_HE)  &
    -k(45)*n(idx_O)  &
    -k(160)*n(idx_Ok)  &
    -k(133)*n(idx_CH2)  &
    -k(22)*n(idx_H2)  &
    -k(29)*n(idx_Hk)  &
    -k(141)*n(idx_OH)  &
    -k(9)*n(idx_HE)  &
    -k(159)*n(idx_Ck)  &
    -k(28)*n(idx_Hk)  &
    -k(3)*n(idx_E)  &
    -k(135)*n(idx_CH2)  &
    +k(193)*n(idx_H2)  &
    -k(21)*n(idx_H2)
pdj(21) =  &
    +k(9)*n(idx_HE)  &
    +k(10)*n(idx_HE)
pdj(22) =  &
    +k(19)*n(idx_H)  &
    +k(29)*n(idx_Hk)  &
    +k(22)*n(idx_H2)  &
    +k(18)*n(idx_H)  &
    +k(21)*n(idx_H2)
pdj(23) =  &
    +k(47)*n(idx_C)
pdj(24) =  &
    +k(45)*n(idx_O)
pdj(27) =  &
    +k(194)*n(idx_H2)
pdj(28) =  &
    +k(130)*n(idx_CH)  &
    +k(133)*n(idx_CH2)  &
    +k(132)*n(idx_CH2)  &
    +k(131)*n(idx_CH)
pdj(29) =  &
    +k(135)*n(idx_CH2)  &
    +k(134)*n(idx_CH2)
pdj(32) =  &
    +k(141)*n(idx_OH)  &
    +k(142)*n(idx_OH)
pdj(33) =  &
    +k(145)*n(idx_H2O)  &
    +k(146)*n(idx_H2O)
pdj(35) =  &
    +k(153)*n(idx_O2)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(20)*1d-3
if(dnn>0.d0) then
nn(20) = n(20) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==21) then
pdj(1) =  &
    -k(7)*n(idx_E)  &
    -k(6)*n(idx_E)  &
    +2.d0*k(7)*n(idx_E)  &
    +k(215)  &
    -k(5)*n(idx_E)
pdj(2) =  &
    -k(161)*n(idx_Hk)
pdj(5) =  &
    +k(161)*n(idx_Hk)  &
    -k(8)*n(idx_H)  &
    +k(138)*n(idx_CH2)  &
    +k(150)*n(idx_H2O)  &
    +k(139)*n(idx_CH2)  &
    +k(149)*n(idx_H2O)  &
    +k(13)*n(idx_H2)  &
    +2.d0*k(14)*n(idx_H2)  &
    +k(143)*n(idx_OH)  &
    +k(144)*n(idx_OH)
pdj(6) =  &
    +k(157)*n(idx_CO)  &
    +k(138)*n(idx_CH2)  &
    +k(46)*n(idx_O)  &
    +k(50)*n(idx_C)  &
    +k(143)*n(idx_OH)  &
    +k(147)*n(idx_H2O)  &
    +k(156)*n(idx_CO)  &
    +k(8)*n(idx_H)  &
    +k(6)*n(idx_E)  &
    +k(49)*n(idx_C)  &
    +k(5)*n(idx_E)  &
    +k(140)*n(idx_C2)  &
    +k(152)*n(idx_H2O)  &
    +k(144)*n(idx_OH)  &
    +k(51)*n(idx_C)  &
    +k(151)*n(idx_H2O)  &
    +k(155)*n(idx_O2)  &
    +k(13)*n(idx_H2)  &
    +k(137)*n(idx_CH2)  &
    +k(148)*n(idx_H2O)  &
    +k(136)*n(idx_CH2)  &
    +k(161)*n(idx_Hk)  &
    +k(150)*n(idx_H2O)  &
    +k(139)*n(idx_CH2)  &
    +k(149)*n(idx_H2O)  &
    +k(12)*n(idx_H2)  &
    +k(154)*n(idx_O2)
pdj(7) =  &
    -k(13)*n(idx_H2)  &
    +k(136)*n(idx_CH2)  &
    -k(12)*n(idx_H2)  &
    +k(137)*n(idx_CH2)  &
    -k(14)*n(idx_H2)
pdj(8) =  &
    -k(49)*n(idx_C)  &
    -k(51)*n(idx_C)  &
    +k(140)*n(idx_C2)  &
    +k(157)*n(idx_CO)  &
    -k(50)*n(idx_C)
pdj(9) =  &
    +k(156)*n(idx_CO)  &
    -k(46)*n(idx_O)  &
    +k(155)*n(idx_O2)
pdj(10) =  &
    -k(143)*n(idx_OH)  &
    +k(147)*n(idx_H2O)  &
    +k(148)*n(idx_H2O)  &
    -k(144)*n(idx_OH)
pdj(11) =  &
    -k(157)*n(idx_CO)  &
    -k(156)*n(idx_CO)
pdj(13) =  &
    -k(137)*n(idx_CH2)  &
    -k(136)*n(idx_CH2)  &
    -k(139)*n(idx_CH2)  &
    -k(138)*n(idx_CH2)
pdj(14) =  &
    -k(140)*n(idx_C2)
pdj(16) =  &
    -k(150)*n(idx_H2O)  &
    -k(147)*n(idx_H2O)  &
    -k(149)*n(idx_H2O)  &
    -k(151)*n(idx_H2O)  &
    -k(148)*n(idx_H2O)  &
    -k(152)*n(idx_H2O)
pdj(17) =  &
    -k(155)*n(idx_O2)  &
    -k(154)*n(idx_O2)
pdj(20) =  &
    +k(148)*n(idx_H2O)  &
    +k(8)*n(idx_H)  &
    +k(147)*n(idx_H2O)  &
    +k(13)*n(idx_H2)
pdj(21) =  &
    -k(8)*n(idx_H)  &
    -k(136)*n(idx_CH2)  &
    -k(5)*n(idx_E)  &
    -k(143)*n(idx_OH)  &
    -k(139)*n(idx_CH2)  &
    -k(13)*n(idx_H2)  &
    -k(161)*n(idx_Hk)  &
    -k(46)*n(idx_O)  &
    -k(6)*n(idx_E)  &
    -k(140)*n(idx_C2)  &
    -k(49)*n(idx_C)  &
    -k(155)*n(idx_O2)  &
    -k(148)*n(idx_H2O)  &
    +k(14)*n(idx_H2)  &
    -k(150)*n(idx_H2O)  &
    -k(50)*n(idx_C)  &
    -k(147)*n(idx_H2O)  &
    -k(215)  &
    -k(138)*n(idx_CH2)  &
    -k(51)*n(idx_C)  &
    -k(157)*n(idx_CO)  &
    -k(154)*n(idx_O2)  &
    -k(152)*n(idx_H2O)  &
    -k(7)*n(idx_E)  &
    -k(149)*n(idx_H2O)  &
    -k(151)*n(idx_H2O)  &
    -k(144)*n(idx_OH)  &
    -k(12)*n(idx_H2)  &
    -k(137)*n(idx_CH2)  &
    -k(14)*n(idx_H2)  &
    -k(156)*n(idx_CO)
pdj(22) =  &
    +k(12)*n(idx_H2)
pdj(23) =  &
    +k(51)*n(idx_C)  &
    +k(156)*n(idx_CO)  &
    +k(49)*n(idx_C)  &
    +k(140)*n(idx_C2)  &
    +k(137)*n(idx_CH2)  &
    +k(50)*n(idx_C)  &
    +k(136)*n(idx_CH2)
pdj(24) =  &
    +k(157)*n(idx_CO)  &
    +k(155)*n(idx_O2)  &
    +k(46)*n(idx_O)  &
    +k(143)*n(idx_OH)  &
    +k(144)*n(idx_OH)
pdj(28) =  &
    +k(139)*n(idx_CH2)  &
    +k(138)*n(idx_CH2)
pdj(32) =  &
    +k(149)*n(idx_H2O)  &
    +k(150)*n(idx_H2O)
pdj(33) =  &
    +k(152)*n(idx_H2O)  &
    +k(151)*n(idx_H2O)
pdj(35) =  &
    +k(154)*n(idx_O2)
pdj(36) =  &
    +k(7)*n(idx_E)  &
    +k(215)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(21)*1d-3
if(dnn>0.d0) then
nn(21) = n(21) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==22) then
pdj(1) =  &
    -k(31)*n(idx_E)  &
    -k(30)*n(idx_E)
pdj(2) =  &
    -k(32)*n(idx_Hk)
pdj(5) =  &
    +2.d0*k(30)*n(idx_E)  &
    -k(20)*n(idx_H)  &
    +k(87)*n(idx_C)  &
    +k(32)*n(idx_Hk)  &
    +k(85)*n(idx_H2)  &
    +2.d0*k(31)*n(idx_E)  &
    +k(102)*n(idx_O)  &
    +k(232)
pdj(7) =  &
    +k(20)*n(idx_H)  &
    +k(32)*n(idx_Hk)  &
    -k(85)*n(idx_H2)
pdj(8) =  &
    -k(87)*n(idx_C)
pdj(9) =  &
    -k(102)*n(idx_O)
pdj(20) =  &
    +k(20)*n(idx_H)  &
    +k(232)
pdj(22) =  &
    -k(85)*n(idx_H2)  &
    -k(20)*n(idx_H)  &
    -k(87)*n(idx_C)  &
    -k(31)*n(idx_E)  &
    -k(30)*n(idx_E)  &
    -k(232)  &
    -k(102)*n(idx_O)  &
    -k(32)*n(idx_Hk)
pdj(27) =  &
    +k(85)*n(idx_H2)
pdj(28) =  &
    +k(87)*n(idx_C)
pdj(32) =  &
    +k(102)*n(idx_O)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(22)*1d-3
if(dnn>0.d0) then
nn(22) = n(22) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==23) then
pdj(1) =  &
    -k(37)*n(idx_E)  &
    -k(39)*n(idx_E)  &
    -k(38)*n(idx_E)
pdj(5) =  &
    +k(90)*n(idx_H2)  &
    +k(107)*n(idx_OH)  &
    +k(108)*n(idx_OH)  &
    -k(48)*n(idx_H)  &
    -k(200)*n(idx_H)  &
    +k(113)*n(idx_H2O)  &
    +k(115)*n(idx_H2O)  &
    +k(114)*n(idx_H2O)
pdj(7) =  &
    -k(201)*n(idx_H2)  &
    -k(90)*n(idx_H2)
pdj(8) =  &
    +k(39)*n(idx_E)  &
    +k(38)*n(idx_E)  &
    +k(37)*n(idx_E)  &
    +k(48)*n(idx_H)  &
    +k(116)*n(idx_H2O)
pdj(9) =  &
    -k(275)*n(idx_O)  &
    -k(276)*n(idx_O)  &
    -k(202)*n(idx_O)  &
    +k(118)*n(idx_O2)  &
    -k(203)*n(idx_O)
pdj(10) =  &
    -k(107)*n(idx_OH)  &
    -k(108)*n(idx_OH)
pdj(11) =  &
    +k(119)*n(idx_O2)
pdj(16) =  &
    -k(114)*n(idx_H2O)  &
    -k(115)*n(idx_H2O)  &
    -k(113)*n(idx_H2O)  &
    -k(116)*n(idx_H2O)
pdj(17) =  &
    -k(119)*n(idx_O2)  &
    -k(118)*n(idx_O2)
pdj(20) =  &
    +k(48)*n(idx_H)
pdj(23) =  &
    -k(39)*n(idx_E)  &
    -k(38)*n(idx_E)  &
    -k(275)*n(idx_O)  &
    -k(203)*n(idx_O)  &
    -k(114)*n(idx_H2O)  &
    -k(37)*n(idx_E)  &
    -k(107)*n(idx_OH)  &
    -k(276)*n(idx_O)  &
    -k(108)*n(idx_OH)  &
    -k(116)*n(idx_H2O)  &
    -k(202)*n(idx_O)  &
    -k(48)*n(idx_H)  &
    -k(90)*n(idx_H2)  &
    -k(201)*n(idx_H2)  &
    -k(113)*n(idx_H2O)  &
    -k(200)*n(idx_H)  &
    -k(115)*n(idx_H2O)  &
    -k(119)*n(idx_O2)  &
    -k(118)*n(idx_O2)
pdj(24) =  &
    +k(119)*n(idx_O2)
pdj(25) =  &
    +k(113)*n(idx_H2O)
pdj(26) =  &
    +k(115)*n(idx_H2O)  &
    +k(114)*n(idx_H2O)
pdj(28) =  &
    +k(90)*n(idx_H2)  &
    +k(200)*n(idx_H)
pdj(29) =  &
    +k(201)*n(idx_H2)
pdj(30) =  &
    +k(118)*n(idx_O2)  &
    +k(107)*n(idx_OH)  &
    +k(202)*n(idx_O)  &
    +k(276)*n(idx_O)  &
    +k(108)*n(idx_OH)  &
    +k(203)*n(idx_O)  &
    +k(275)*n(idx_O)
pdj(33) =  &
    +k(116)*n(idx_H2O)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(23)*1d-3
if(dnn>0.d0) then
nn(23) = n(23) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==24) then
pdj(1) =  &
    -k(40)*n(idx_E)  &
    -k(41)*n(idx_E)
pdj(5) =  &
    -k(44)*n(idx_H)  &
    +k(101)*n(idx_H2)
pdj(7) =  &
    -k(101)*n(idx_H2)
pdj(8) =  &
    +k(100)*n(idx_C2)  &
    -k(277)*n(idx_C)  &
    -k(278)*n(idx_C)
pdj(9) =  &
    +k(40)*n(idx_E)  &
    +k(41)*n(idx_E)  &
    +k(44)*n(idx_H)
pdj(14) =  &
    -k(100)*n(idx_C2)
pdj(20) =  &
    +k(44)*n(idx_H)
pdj(24) =  &
    -k(100)*n(idx_C2)  &
    -k(101)*n(idx_H2)  &
    -k(44)*n(idx_H)  &
    -k(278)*n(idx_C)  &
    -k(277)*n(idx_C)  &
    -k(40)*n(idx_E)  &
    -k(41)*n(idx_E)
pdj(30) =  &
    +k(100)*n(idx_C2)  &
    +k(278)*n(idx_C)  &
    +k(277)*n(idx_C)
pdj(32) =  &
    +k(101)*n(idx_H2)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(24)*1d-3
if(dnn>0.d0) then
nn(24) = n(24) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==25) then
pdj(1) =  &
    -k(183)*n(idx_E)
pdj(5) =  &
    +k(183)*n(idx_E)
pdj(7) =  &
    +k(53)*n(idx_H2)  &
    -k(53)*n(idx_H2)
pdj(11) =  &
    +k(54)*n(idx_CO)  &
    +k(183)*n(idx_E)  &
    -k(55)*n(idx_CO)  &
    +k(55)*n(idx_CO)  &
    -k(54)*n(idx_CO)
pdj(25) =  &
    -k(183)*n(idx_E)  &
    -k(55)*n(idx_CO)  &
    -k(53)*n(idx_H2)  &
    -k(54)*n(idx_CO)
pdj(26) =  &
    +k(53)*n(idx_H2)  &
    +k(54)*n(idx_CO)  &
    +k(55)*n(idx_CO)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(25)*1d-3
if(dnn>0.d0) then
nn(25) = n(25) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==26) then
pdj(1) =  &
    -k(182)*n(idx_E)  &
    -k(181)*n(idx_E)
pdj(5) =  &
    +k(181)*n(idx_E)
pdj(8) =  &
    +k(182)*n(idx_E)  &
    -k(127)*n(idx_C)
pdj(10) =  &
    +k(182)*n(idx_E)
pdj(11) =  &
    +k(127)*n(idx_C)  &
    +k(129)*n(idx_H2O)  &
    +k(128)*n(idx_H2O)  &
    +k(181)*n(idx_E)
pdj(16) =  &
    -k(129)*n(idx_H2O)  &
    -k(128)*n(idx_H2O)
pdj(26) =  &
    -k(129)*n(idx_H2O)  &
    -k(182)*n(idx_E)  &
    -k(127)*n(idx_C)  &
    -k(128)*n(idx_H2O)  &
    -k(181)*n(idx_E)
pdj(28) =  &
    +k(127)*n(idx_C)
pdj(34) =  &
    +k(129)*n(idx_H2O)  &
    +k(128)*n(idx_H2O)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(26)*1d-3
if(dnn>0.d0) then
nn(26) = n(26) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==27) then
pdj(1) =  &
    -k(162)*n(idx_E)  &
    -k(163)*n(idx_E)
pdj(5) =  &
    +3.d0*k(163)*n(idx_E)  &
    +k(104)*n(idx_O)  &
    -k(86)*n(idx_H)  &
    +k(234)  &
    +k(162)*n(idx_E)  &
    +k(89)*n(idx_C)
pdj(7) =  &
    +k(233)  &
    +k(86)*n(idx_H)  &
    +k(125)*n(idx_CO)  &
    +k(124)*n(idx_CO)  &
    +k(111)*n(idx_H2O)  &
    +k(123)*n(idx_CO)  &
    +k(126)*n(idx_CO)  &
    +k(105)*n(idx_OH)  &
    +k(106)*n(idx_OH)  &
    +k(88)*n(idx_C)  &
    +k(103)*n(idx_O)  &
    +k(162)*n(idx_E)  &
    +k(112)*n(idx_H2O)
pdj(8) =  &
    -k(88)*n(idx_C)  &
    -k(89)*n(idx_C)
pdj(9) =  &
    -k(104)*n(idx_O)  &
    -k(103)*n(idx_O)
pdj(10) =  &
    -k(105)*n(idx_OH)  &
    -k(106)*n(idx_OH)
pdj(11) =  &
    -k(126)*n(idx_CO)  &
    -k(125)*n(idx_CO)  &
    -k(124)*n(idx_CO)  &
    -k(123)*n(idx_CO)
pdj(16) =  &
    -k(111)*n(idx_H2O)  &
    -k(112)*n(idx_H2O)
pdj(20) =  &
    +k(233)
pdj(22) =  &
    +k(86)*n(idx_H)  &
    +k(234)
pdj(25) =  &
    +k(125)*n(idx_CO)  &
    +k(126)*n(idx_CO)
pdj(26) =  &
    +k(123)*n(idx_CO)  &
    +k(124)*n(idx_CO)
pdj(27) =  &
    -k(124)*n(idx_CO)  &
    -k(123)*n(idx_CO)  &
    -k(89)*n(idx_C)  &
    -k(88)*n(idx_C)  &
    -k(125)*n(idx_CO)  &
    -k(112)*n(idx_H2O)  &
    -k(233)  &
    -k(111)*n(idx_H2O)  &
    -k(86)*n(idx_H)  &
    -k(162)*n(idx_E)  &
    -k(163)*n(idx_E)  &
    -k(106)*n(idx_OH)  &
    -k(126)*n(idx_CO)  &
    -k(104)*n(idx_O)  &
    -k(234)  &
    -k(105)*n(idx_OH)  &
    -k(103)*n(idx_O)
pdj(28) =  &
    +k(88)*n(idx_C)
pdj(29) =  &
    +k(89)*n(idx_C)
pdj(32) =  &
    +k(103)*n(idx_O)
pdj(33) =  &
    +k(104)*n(idx_O)  &
    +k(105)*n(idx_OH)  &
    +k(106)*n(idx_OH)
pdj(34) =  &
    +k(112)*n(idx_H2O)  &
    +k(111)*n(idx_H2O)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(27)*1d-3
if(dnn>0.d0) then
nn(27) = n(27) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==28) then
pdj(1) =  &
    -k(164)*n(idx_E)
pdj(5) =  &
    +k(164)*n(idx_E)  &
    +k(92)*n(idx_H2)  &
    +k(93)*n(idx_O)  &
    -k(91)*n(idx_H)
pdj(7) =  &
    +k(91)*n(idx_H)  &
    -k(92)*n(idx_H2)
pdj(8) =  &
    +k(164)*n(idx_E)  &
    +k(236)
pdj(9) =  &
    -k(93)*n(idx_O)
pdj(20) =  &
    +k(236)
pdj(23) =  &
    +k(91)*n(idx_H)
pdj(28) =  &
    -k(236)  &
    -k(164)*n(idx_E)  &
    -k(93)*n(idx_O)  &
    -k(92)*n(idx_H2)  &
    -k(91)*n(idx_H)
pdj(29) =  &
    +k(92)*n(idx_H2)
pdj(30) =  &
    +k(93)*n(idx_O)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(28)*1d-3
if(dnn>0.d0) then
nn(28) = n(28) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==29) then
pdj(1) =  &
    -k(165)*n(idx_E)  &
    -k(166)*n(idx_E)  &
    -k(167)*n(idx_E)
pdj(5) =  &
    -k(94)*n(idx_H)  &
    +k(96)*n(idx_O)  &
    +k(95)*n(idx_H2)  &
    +k(165)*n(idx_E)  &
    +k(239)  &
    +2.d0*k(167)*n(idx_E)
pdj(7) =  &
    +k(94)*n(idx_H)  &
    -k(95)*n(idx_H2)  &
    +k(166)*n(idx_E)
pdj(8) =  &
    +k(166)*n(idx_E)  &
    +k(167)*n(idx_E)
pdj(9) =  &
    -k(96)*n(idx_O)
pdj(10) =  &
    +k(120)*n(idx_O2)
pdj(12) =  &
    +k(165)*n(idx_E)
pdj(17) =  &
    -k(120)*n(idx_O2)
pdj(26) =  &
    +k(120)*n(idx_O2)  &
    +k(96)*n(idx_O)
pdj(28) =  &
    +k(94)*n(idx_H)  &
    +k(239)
pdj(29) =  &
    -k(94)*n(idx_H)  &
    -k(95)*n(idx_H2)  &
    -k(166)*n(idx_E)  &
    -k(96)*n(idx_O)  &
    -k(167)*n(idx_E)  &
    -k(239)  &
    -k(165)*n(idx_E)  &
    -k(120)*n(idx_O2)
pdj(31) =  &
    +k(95)*n(idx_H2)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(29)*1d-3
if(dnn>0.d0) then
nn(29) = n(29) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==30) then
pdj(1) =  &
    -k(180)*n(idx_E)
pdj(5) =  &
    -k(158)*n(idx_H)
pdj(8) =  &
    +k(180)*n(idx_E)
pdj(9) =  &
    +k(180)*n(idx_E)
pdj(11) =  &
    +k(158)*n(idx_H)
pdj(20) =  &
    +k(158)*n(idx_H)
pdj(30) =  &
    -k(158)*n(idx_H)  &
    -k(180)*n(idx_E)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(30)*1d-3
if(dnn>0.d0) then
nn(30) = n(30) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==31) then
pdj(1) =  &
    -k(170)*n(idx_E)  &
    -k(168)*n(idx_E)  &
    -k(169)*n(idx_E)
pdj(5) =  &
    -k(97)*n(idx_H)  &
    +k(240)  &
    +k(168)*n(idx_E)  &
    +2.d0*k(170)*n(idx_E)
pdj(7) =  &
    +k(99)*n(idx_O)  &
    +k(241)  &
    +k(169)*n(idx_E)  &
    +k(98)*n(idx_O)  &
    +k(97)*n(idx_H)
pdj(9) =  &
    -k(99)*n(idx_O)  &
    -k(98)*n(idx_O)
pdj(12) =  &
    +k(169)*n(idx_E)  &
    +k(170)*n(idx_E)
pdj(13) =  &
    +k(168)*n(idx_E)
pdj(25) =  &
    +k(98)*n(idx_O)
pdj(26) =  &
    +k(99)*n(idx_O)
pdj(28) =  &
    +k(241)
pdj(29) =  &
    +k(240)  &
    +k(97)*n(idx_H)
pdj(31) =  &
    -k(169)*n(idx_E)  &
    -k(97)*n(idx_H)  &
    -k(99)*n(idx_O)  &
    -k(98)*n(idx_O)  &
    -k(170)*n(idx_E)  &
    -k(240)  &
    -k(241)  &
    -k(168)*n(idx_E)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(31)*1d-3
if(dnn>0.d0) then
nn(31) = n(31) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==32) then
pdj(1) =  &
    -k(171)*n(idx_E)
pdj(5) =  &
    +k(109)*n(idx_H2)  &
    +k(171)*n(idx_E)
pdj(7) =  &
    -k(109)*n(idx_H2)
pdj(9) =  &
    +k(243)  &
    +k(171)*n(idx_E)
pdj(20) =  &
    +k(243)
pdj(32) =  &
    -k(171)*n(idx_E)  &
    -k(109)*n(idx_H2)  &
    -k(243)
pdj(33) =  &
    +k(109)*n(idx_H2)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(32)*1d-3
if(dnn>0.d0) then
nn(32) = n(32) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==33) then
pdj(1) =  &
    -k(174)*n(idx_E)  &
    -k(173)*n(idx_E)  &
    -k(172)*n(idx_E)
pdj(5) =  &
    +k(173)*n(idx_E)  &
    +k(247)  &
    +2.d0*k(174)*n(idx_E)  &
    +k(110)*n(idx_H2)
pdj(7) =  &
    +k(246)  &
    +k(172)*n(idx_E)  &
    -k(110)*n(idx_H2)
pdj(9) =  &
    +k(172)*n(idx_E)  &
    +k(174)*n(idx_E)  &
    +k(244)
pdj(10) =  &
    +k(173)*n(idx_E)  &
    +k(245)
pdj(20) =  &
    +k(245)
pdj(22) =  &
    +k(244)
pdj(24) =  &
    +k(246)
pdj(32) =  &
    +k(247)
pdj(33) =  &
    -k(110)*n(idx_H2)  &
    -k(173)*n(idx_E)  &
    -k(245)  &
    -k(244)  &
    -k(174)*n(idx_E)  &
    -k(172)*n(idx_E)  &
    -k(247)  &
    -k(246)
pdj(34) =  &
    +k(110)*n(idx_H2)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(33)*1d-3
if(dnn>0.d0) then
nn(33) = n(33) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==34) then
pdj(1) =  &
    -k(178)*n(idx_E)  &
    -k(177)*n(idx_E)  &
    -k(176)*n(idx_E)  &
    -k(175)*n(idx_E)
pdj(5) =  &
    +2.d0*k(175)*n(idx_E)  &
    +k(250)  &
    +k(177)*n(idx_E)  &
    +k(176)*n(idx_E)
pdj(7) =  &
    +k(251)  &
    +k(117)*n(idx_C)  &
    +k(176)*n(idx_E)  &
    +k(178)*n(idx_E)
pdj(8) =  &
    -k(117)*n(idx_C)
pdj(9) =  &
    +k(176)*n(idx_E)
pdj(10) =  &
    +k(249)  &
    +k(175)*n(idx_E)  &
    +k(178)*n(idx_E)
pdj(16) =  &
    +k(248)  &
    +k(177)*n(idx_E)
pdj(20) =  &
    +k(248)
pdj(22) =  &
    +k(249)
pdj(26) =  &
    +k(117)*n(idx_C)
pdj(32) =  &
    +k(251)
pdj(33) =  &
    +k(250)
pdj(34) =  &
    -k(248)  &
    -k(178)*n(idx_E)  &
    -k(117)*n(idx_C)  &
    -k(175)*n(idx_E)  &
    -k(250)  &
    -k(251)  &
    -k(176)*n(idx_E)  &
    -k(177)*n(idx_E)  &
    -k(249)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(34)*1d-3
if(dnn>0.d0) then
nn(34) = n(34) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==35) then
pdj(1) =  &
    -k(179)*n(idx_E)
pdj(8) =  &
    -k(121)*n(idx_C)  &
    -k(122)*n(idx_C)
pdj(9) =  &
    +k(121)*n(idx_C)  &
    +2.d0*k(179)*n(idx_E)
pdj(17) =  &
    +k(122)*n(idx_C)
pdj(23) =  &
    +k(122)*n(idx_C)
pdj(30) =  &
    +k(121)*n(idx_C)
pdj(35) =  &
    -k(179)*n(idx_E)  &
    -k(121)*n(idx_C)  &
    -k(122)*n(idx_C)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(35)*1d-3
if(dnn>0.d0) then
nn(35) = n(35) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==36) then
pdj(1) =  &
    -k(15)*n(idx_E)
pdj(21) =  &
    +k(15)*n(idx_E)
pdj(36) =  &
    -k(15)*n(idx_E)
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(36)*1d-3
if(dnn>0.d0) then
nn(36) = n(36) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pdj(idx_Tgas) = (dn1-dn0)/dnn
end if

elseif(j==37) then
pdj(39) = 0.d0
elseif(j==38) then
pdj(39) = 0.d0
elseif(j==39) then
!use fex to compute temperature-dependent Jacobian
dnn = n(idx_Tgas)*1d-3
nn(:) = n(:)
nn(idx_Tgas) = n(idx_Tgas) + dnn
call fex(neq,tt,nn(:),dn(:))
do i=1,neq-1
pdj(i) = dn(i) / dnn
end do
elseif(j==40) then
pdj(39) = 0.d0
end if

return
end subroutine jes

!*************************
subroutine jex(neq,t,n,ml,mu,pd,npd)
use krome_commons
use krome_tabs
use krome_cooling
use krome_heating
use krome_constants
use krome_subs
use krome_gadiab
implicit none
real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
real*8::krome_gamma,nn(neq),nH2dust
integer::neq,ml,mu,npd

Tgas = n(idx_Tgas)
npd = neq
k(:) = coe_tab(n(:))
pd(:,:) = 0d0
krome_gamma = gamma_index(n(:))

!d[E_dot]/d[E]
pd(1,1) =  &
    -k(6)*n(idx_HEj)  &
    -k(178)*n(idx_H3Oj)  &
    -k(168)*n(idx_CH3j)  &
    -k(174)*n(idx_H2Oj)  &
    -k(195)*n(idx_C)  &
    -k(169)*n(idx_CH3j)  &
    +2.d0*k(4)*n(idx_HE)  &
    -k(25)*n(idx_Hk)  &
    -k(170)*n(idx_CH3j)  &
    -k(1)*n(idx_H)  &
    -k(166)*n(idx_CH2j)  &
    -k(16)*n(idx_H)  &
    -k(39)*n(idx_Cj)  &
    -k(183)*n(idx_HOCj)  &
    -k(4)*n(idx_HE)  &
    -k(180)*n(idx_COj)  &
    -k(5)*n(idx_HEj)  &
    -k(181)*n(idx_HCOj)  &
    -k(162)*n(idx_H3j)  &
    -k(38)*n(idx_Cj)  &
    -k(43)*n(idx_O)  &
    +2.d0*k(42)*n(idx_C)  &
    -k(40)*n(idx_Oj)  &
    -k(7)*n(idx_HEj)  &
    -k(176)*n(idx_H3Oj)  &
    -k(163)*n(idx_H3j)  &
    -k(171)*n(idx_OHj)  &
    -k(175)*n(idx_H3Oj)  &
    +2.d0*k(7)*n(idx_HEj)  &
    -k(177)*n(idx_H3Oj)  &
    -k(2)*n(idx_Hj)  &
    -k(167)*n(idx_CH2j)  &
    -k(41)*n(idx_Oj)  &
    -k(182)*n(idx_HCOj)  &
    -k(173)*n(idx_H2Oj)  &
    -k(42)*n(idx_C)  &
    -k(37)*n(idx_Cj)  &
    -k(3)*n(idx_Hj)  &
    -k(164)*n(idx_CHj)  &
    -k(31)*n(idx_H2j)  &
    +2.d0*k(25)*n(idx_Hk)  &
    +k(23)*n(idx_H2)  &
    -k(23)*n(idx_H2)  &
    -k(172)*n(idx_H2Oj)  &
    +2.d0*k(1)*n(idx_H)  &
    -k(179)*n(idx_O2j)  &
    -k(204)*n(idx_O)  &
    -k(30)*n(idx_H2j)  &
    +2.d0*k(43)*n(idx_O)  &
    -k(165)*n(idx_CH2j)  &
    -k(15)*n(idx_HEjj)

!d[H-_dot]/d[E]
pd(2,1) =  &
    +k(16)*n(idx_H)  &
    -k(25)*n(idx_Hk)

!d[C-_dot]/d[E]
pd(3,1) =  &
    +k(195)*n(idx_C)

!d[O-_dot]/d[E]
pd(4,1) =  &
    +k(204)*n(idx_O)

!d[H_dot]/d[E]
pd(5,1) =  &
    +k(183)*n(idx_HOCj)  &
    +k(3)*n(idx_Hj)  &
    +k(173)*n(idx_H2Oj)  &
    +k(171)*n(idx_OHj)  &
    +2.d0*k(174)*n(idx_H2Oj)  &
    +k(177)*n(idx_H3Oj)  &
    +k(165)*n(idx_CH2j)  &
    -k(1)*n(idx_H)  &
    -k(16)*n(idx_H)  &
    +k(168)*n(idx_CH3j)  &
    +3.d0*k(163)*n(idx_H3j)  &
    +2.d0*k(167)*n(idx_CH2j)  &
    +k(162)*n(idx_H3j)  &
    +k(2)*n(idx_Hj)  &
    +k(176)*n(idx_H3Oj)  &
    +2.d0*k(31)*n(idx_H2j)  &
    +2.d0*k(30)*n(idx_H2j)  &
    +k(164)*n(idx_CHj)  &
    +2.d0*k(170)*n(idx_CH3j)  &
    +2.d0*k(175)*n(idx_H3Oj)  &
    +k(181)*n(idx_HCOj)  &
    +k(25)*n(idx_Hk)  &
    +2.d0*k(23)*n(idx_H2)

!d[HE_dot]/d[E]
pd(6,1) =  &
    -k(4)*n(idx_HE)  &
    +k(6)*n(idx_HEj)  &
    +k(5)*n(idx_HEj)

!d[H2_dot]/d[E]
pd(7,1) =  &
    +k(166)*n(idx_CH2j)  &
    -k(23)*n(idx_H2)  &
    +k(172)*n(idx_H2Oj)  &
    +k(162)*n(idx_H3j)  &
    +k(169)*n(idx_CH3j)  &
    +k(178)*n(idx_H3Oj)  &
    +k(176)*n(idx_H3Oj)

!d[C_dot]/d[E]
pd(8,1) =  &
    +k(38)*n(idx_Cj)  &
    -k(42)*n(idx_C)  &
    +k(167)*n(idx_CH2j)  &
    +k(180)*n(idx_COj)  &
    +k(37)*n(idx_Cj)  &
    +k(164)*n(idx_CHj)  &
    +k(166)*n(idx_CH2j)  &
    -k(195)*n(idx_C)  &
    +k(182)*n(idx_HCOj)  &
    +k(39)*n(idx_Cj)

!d[O_dot]/d[E]
pd(9,1) =  &
    -k(204)*n(idx_O)  &
    +k(180)*n(idx_COj)  &
    +k(41)*n(idx_Oj)  &
    -k(43)*n(idx_O)  &
    +k(171)*n(idx_OHj)  &
    +2.d0*k(179)*n(idx_O2j)  &
    +k(40)*n(idx_Oj)  &
    +k(172)*n(idx_H2Oj)  &
    +k(174)*n(idx_H2Oj)  &
    +k(176)*n(idx_H3Oj)

!d[OH_dot]/d[E]
pd(10,1) =  &
    +k(182)*n(idx_HCOj)  &
    +k(178)*n(idx_H3Oj)  &
    +k(175)*n(idx_H3Oj)  &
    +k(173)*n(idx_H2Oj)

!d[CO_dot]/d[E]
pd(11,1) =  &
    +k(183)*n(idx_HOCj)  &
    +k(181)*n(idx_HCOj)

!d[CH_dot]/d[E]
pd(12,1) =  &
    +k(165)*n(idx_CH2j)  &
    +k(169)*n(idx_CH3j)  &
    +k(170)*n(idx_CH3j)

!d[CH2_dot]/d[E]
pd(13,1) =  &
    +k(168)*n(idx_CH3j)

!d[H2O_dot]/d[E]
pd(16,1) =  &
    +k(177)*n(idx_H3Oj)

!d[H+_dot]/d[E]
pd(20,1) =  &
    -k(2)*n(idx_Hj)  &
    -k(3)*n(idx_Hj)  &
    +k(1)*n(idx_H)

!d[HE+_dot]/d[E]
pd(21,1) =  &
    +k(4)*n(idx_HE)  &
    -k(5)*n(idx_HEj)  &
    +k(15)*n(idx_HEjj)  &
    -k(7)*n(idx_HEj)  &
    -k(6)*n(idx_HEj)

!d[H2+_dot]/d[E]
pd(22,1) =  &
    -k(30)*n(idx_H2j)  &
    -k(31)*n(idx_H2j)

!d[C+_dot]/d[E]
pd(23,1) =  &
    -k(38)*n(idx_Cj)  &
    -k(37)*n(idx_Cj)  &
    +k(42)*n(idx_C)  &
    -k(39)*n(idx_Cj)

!d[O+_dot]/d[E]
pd(24,1) =  &
    -k(41)*n(idx_Oj)  &
    -k(40)*n(idx_Oj)  &
    +k(43)*n(idx_O)

!d[HOC+_dot]/d[E]
pd(25,1) =  &
    -k(183)*n(idx_HOCj)

!d[HCO+_dot]/d[E]
pd(26,1) =  &
    -k(182)*n(idx_HCOj)  &
    -k(181)*n(idx_HCOj)

!d[H3+_dot]/d[E]
pd(27,1) =  &
    -k(162)*n(idx_H3j)  &
    -k(163)*n(idx_H3j)

!d[CH+_dot]/d[E]
pd(28,1) =  &
    -k(164)*n(idx_CHj)

!d[CH2+_dot]/d[E]
pd(29,1) =  &
    -k(166)*n(idx_CH2j)  &
    -k(167)*n(idx_CH2j)  &
    -k(165)*n(idx_CH2j)

!d[CO+_dot]/d[E]
pd(30,1) =  &
    -k(180)*n(idx_COj)

!d[CH3+_dot]/d[E]
pd(31,1) =  &
    -k(168)*n(idx_CH3j)  &
    -k(170)*n(idx_CH3j)  &
    -k(169)*n(idx_CH3j)

!d[OH+_dot]/d[E]
pd(32,1) =  &
    -k(171)*n(idx_OHj)

!d[H2O+_dot]/d[E]
pd(33,1) =  &
    -k(174)*n(idx_H2Oj)  &
    -k(173)*n(idx_H2Oj)  &
    -k(172)*n(idx_H2Oj)

!d[H3O+_dot]/d[E]
pd(34,1) =  &
    -k(176)*n(idx_H3Oj)  &
    -k(178)*n(idx_H3Oj)  &
    -k(175)*n(idx_H3Oj)  &
    -k(177)*n(idx_H3Oj)

!d[O2+_dot]/d[E]
pd(35,1) =  &
    -k(179)*n(idx_O2j)

!d[HE++_dot]/d[E]
pd(36,1) =  &
    +k(7)*n(idx_HEj)  &
    -k(15)*n(idx_HEjj)

!d[Tgas_dot]/d[E]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(1)*1d-3
if(dnn>0.d0) then
nn(1) = n(1) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,1) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H-]
pd(1,2) =  &
    +k(219)  &
    +k(29)*n(idx_Hj)  &
    +k(184)*n(idx_C)  &
    +k(26)*n(idx_H)  &
    +k(27)*n(idx_H)  &
    -k(25)*n(idx_E)  &
    +k(186)*n(idx_OH)  &
    +k(17)*n(idx_H)  &
    +k(185)*n(idx_O)  &
    +2.d0*k(25)*n(idx_E)

!d[H-_dot]/d[H-]
pd(2,2) =  &
    -k(27)*n(idx_H)  &
    -k(184)*n(idx_C)  &
    -k(161)*n(idx_HEj)  &
    -k(28)*n(idx_Hj)  &
    -k(25)*n(idx_E)  &
    -k(219)  &
    -k(29)*n(idx_Hj)  &
    -k(17)*n(idx_H)  &
    -k(32)*n(idx_H2j)  &
    -k(186)*n(idx_OH)  &
    -k(26)*n(idx_H)  &
    -k(185)*n(idx_O)

!d[H_dot]/d[H-]
pd(5,2) =  &
    -k(27)*n(idx_H)  &
    +k(219)  &
    +2.d0*k(26)*n(idx_H)  &
    +2.d0*k(27)*n(idx_H)  &
    +2.d0*k(28)*n(idx_Hj)  &
    -k(17)*n(idx_H)  &
    +k(161)*n(idx_HEj)  &
    +k(32)*n(idx_H2j)  &
    +k(25)*n(idx_E)  &
    -k(26)*n(idx_H)

!d[HE_dot]/d[H-]
pd(6,2) =  &
    +k(161)*n(idx_HEj)

!d[H2_dot]/d[H-]
pd(7,2) =  &
    +k(32)*n(idx_H2j)  &
    +k(17)*n(idx_H)

!d[C_dot]/d[H-]
pd(8,2) =  &
    -k(184)*n(idx_C)

!d[O_dot]/d[H-]
pd(9,2) =  &
    -k(185)*n(idx_O)

!d[OH_dot]/d[H-]
pd(10,2) =  &
    +k(185)*n(idx_O)  &
    -k(186)*n(idx_OH)

!d[CH_dot]/d[H-]
pd(12,2) =  &
    +k(184)*n(idx_C)

!d[H2O_dot]/d[H-]
pd(16,2) =  &
    +k(186)*n(idx_OH)

!d[H+_dot]/d[H-]
pd(20,2) =  &
    -k(28)*n(idx_Hj)  &
    -k(29)*n(idx_Hj)

!d[HE+_dot]/d[H-]
pd(21,2) =  &
    -k(161)*n(idx_HEj)

!d[H2+_dot]/d[H-]
pd(22,2) =  &
    -k(32)*n(idx_H2j)  &
    +k(29)*n(idx_Hj)

!d[Tgas_dot]/d[H-]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(2)*1d-3
if(dnn>0.d0) then
nn(2) = n(2) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,2) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[C-]
pd(1,3) =  &
    +k(189)*n(idx_O)  &
    +k(187)*n(idx_H)  &
    +k(235)  &
    +k(188)*n(idx_H2)

!d[C-_dot]/d[C-]
pd(3,3) =  &
    -k(188)*n(idx_H2)  &
    -k(159)*n(idx_Hj)  &
    -k(189)*n(idx_O)  &
    -k(187)*n(idx_H)  &
    -k(235)

!d[H_dot]/d[C-]
pd(5,3) =  &
    -k(187)*n(idx_H)  &
    +k(159)*n(idx_Hj)

!d[H2_dot]/d[C-]
pd(7,3) =  &
    -k(188)*n(idx_H2)

!d[C_dot]/d[C-]
pd(8,3) =  &
    +k(235)  &
    +k(159)*n(idx_Hj)

!d[O_dot]/d[C-]
pd(9,3) =  &
    -k(189)*n(idx_O)

!d[CO_dot]/d[C-]
pd(11,3) =  &
    +k(189)*n(idx_O)

!d[CH_dot]/d[C-]
pd(12,3) =  &
    +k(187)*n(idx_H)

!d[CH2_dot]/d[C-]
pd(13,3) =  &
    +k(188)*n(idx_H2)

!d[H+_dot]/d[C-]
pd(20,3) =  &
    -k(159)*n(idx_Hj)

!d[Tgas_dot]/d[C-]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(3)*1d-3
if(dnn>0.d0) then
nn(3) = n(3) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,3) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[O-]
pd(1,4) =  &
    +k(191)*n(idx_H2)  &
    +k(242)  &
    +k(192)*n(idx_C)  &
    +k(190)*n(idx_H)

!d[O-_dot]/d[O-]
pd(4,4) =  &
    -k(242)  &
    -k(192)*n(idx_C)  &
    -k(190)*n(idx_H)  &
    -k(191)*n(idx_H2)  &
    -k(160)*n(idx_Hj)

!d[H_dot]/d[O-]
pd(5,4) =  &
    -k(190)*n(idx_H)  &
    +k(160)*n(idx_Hj)

!d[H2_dot]/d[O-]
pd(7,4) =  &
    -k(191)*n(idx_H2)

!d[C_dot]/d[O-]
pd(8,4) =  &
    -k(192)*n(idx_C)

!d[O_dot]/d[O-]
pd(9,4) =  &
    +k(242)  &
    +k(160)*n(idx_Hj)

!d[OH_dot]/d[O-]
pd(10,4) =  &
    +k(190)*n(idx_H)

!d[CO_dot]/d[O-]
pd(11,4) =  &
    +k(192)*n(idx_C)

!d[H2O_dot]/d[O-]
pd(16,4) =  &
    +k(191)*n(idx_H2)

!d[H+_dot]/d[O-]
pd(20,4) =  &
    -k(160)*n(idx_Hj)

!d[Tgas_dot]/d[O-]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(4)*1d-3
if(dnn>0.d0) then
nn(4) = n(4) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,4) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H]
pd(1,5) =  &
    +k(26)*n(idx_Hk)  &
    -k(1)*n(idx_E)  &
    +k(27)*n(idx_Hk)  &
    +k(190)*n(idx_Ok)  &
    +2.d0*k(1)*n(idx_E)  &
    -k(16)*n(idx_E)  &
    +k(252)  &
    +k(187)*n(idx_Ck)  &
    +k(213)  &
    +k(17)*n(idx_Hk)

!d[H-_dot]/d[H]
pd(2,5) =  &
    -k(26)*n(idx_Hk)  &
    +k(16)*n(idx_E)  &
    -k(27)*n(idx_Hk)  &
    -k(17)*n(idx_Hk)

!d[C-_dot]/d[H]
pd(3,5) =  &
    -k(187)*n(idx_Ck)

!d[O-_dot]/d[H]
pd(4,5) =  &
    -k(190)*n(idx_Ok)

!d[H_dot]/d[H]
pd(5,5) =  &
    -k(57)*n(idx_CH)  &
    -k(200)*n(idx_Cj)  &
    -k(94)*n(idx_CH2j)  &
    -k(26)*n(idx_Hk)  &
    +2.d0*k(52)*n(idx_OH)  &
    -k(205)*n(idx_O)  &
    -k(17)*n(idx_Hk)  &
    -k(187)*n(idx_Ck)  &
    -9.d0*k(35)*n(idx_H)*n(idx_H)  &
    +3.d0*k(35)*n(idx_H)*n(idx_H)  &
    -k(27)*n(idx_Hk)  &
    -k(213)  &
    -k(280)*n(idx_OH)  &
    -4.d0*k(36)*n(idx_H2)*n(idx_H)  &
    -k(207)*n(idx_OH)  &
    -k(44)*n(idx_Oj)  &
    -k(18)*n(idx_Hj)  &
    +2.d0*k(27)*n(idx_Hk)  &
    -k(158)*n(idx_COj)  &
    -k(24)*n(idx_H2)  &
    -k(8)*n(idx_HEj)  &
    -k(84)*n(idx_CO)  &
    -k(19)*n(idx_Hj)  &
    -k(48)*n(idx_Cj)  &
    +2.d0*k(26)*n(idx_Hk)  &
    -k(63)*n(idx_CH2)  &
    -k(1)*n(idx_E)  &
    -k(252)  &
    -k(52)*n(idx_OH)  &
    +3.d0*k(24)*n(idx_H2)  &
    -4.d0*k(34)*n(idx_H)*n(idx_HE)  &
    -k(279)*n(idx_O)  &
    -k(79)*n(idx_H2O)  &
    -k(97)*n(idx_CH3j)  &
    -k(72)*n(idx_OH)  &
    -k(71)*n(idx_OH)  &
    -k(91)*n(idx_CHj)  &
    -k(196)*n(idx_C)  &
    -k(86)*n(idx_H3j)  &
    -k(80)*n(idx_O2)  &
    -k(20)*n(idx_H2j)  &
    -k(16)*n(idx_E)  &
    -k(190)*n(idx_Ok)

!d[HE_dot]/d[H]
pd(6,5) =  &
    +2.d0*k(34)*n(idx_H)*n(idx_HE)  &
    -2.d0*k(34)*n(idx_H)*n(idx_HE)  &
    +k(8)*n(idx_HEj)

!d[H2_dot]/d[H]
pd(7,5) =  &
    -2.d0*k(36)*n(idx_H2)*n(idx_H)  &
    +k(63)*n(idx_CH2)  &
    +k(94)*n(idx_CH2j)  &
    +k(57)*n(idx_CH)  &
    +k(20)*n(idx_H2j)  &
    -k(24)*n(idx_H2)  &
    +k(97)*n(idx_CH3j)  &
    +k(72)*n(idx_OH)  &
    +k(71)*n(idx_OH)  &
    +k(86)*n(idx_H3j)  &
    +3.d0*k(35)*n(idx_H)*n(idx_H)  &
    +2.d0*k(34)*n(idx_H)*n(idx_HE)  &
    +4.d0*k(36)*n(idx_H2)*n(idx_H)  &
    +k(91)*n(idx_CHj)  &
    +k(79)*n(idx_H2O)  &
    +k(17)*n(idx_Hk)

!d[C_dot]/d[H]
pd(8,5) =  &
    +k(84)*n(idx_CO)  &
    +k(48)*n(idx_Cj)  &
    +k(57)*n(idx_CH)  &
    -k(196)*n(idx_C)

!d[O_dot]/d[H]
pd(9,5) =  &
    +k(44)*n(idx_Oj)  &
    +k(52)*n(idx_OH)  &
    -k(205)*n(idx_O)  &
    +k(80)*n(idx_O2)  &
    +k(72)*n(idx_OH)  &
    +k(71)*n(idx_OH)  &
    -k(279)*n(idx_O)

!d[OH_dot]/d[H]
pd(10,5) =  &
    +k(84)*n(idx_CO)  &
    -k(280)*n(idx_OH)  &
    -k(72)*n(idx_OH)  &
    -k(71)*n(idx_OH)  &
    +k(190)*n(idx_Ok)  &
    -k(52)*n(idx_OH)  &
    -k(207)*n(idx_OH)  &
    +k(80)*n(idx_O2)  &
    +k(205)*n(idx_O)  &
    +k(79)*n(idx_H2O)  &
    +k(279)*n(idx_O)

!d[CO_dot]/d[H]
pd(11,5) =  &
    +k(158)*n(idx_COj)  &
    -k(84)*n(idx_CO)

!d[CH_dot]/d[H]
pd(12,5) =  &
    +k(196)*n(idx_C)  &
    -k(57)*n(idx_CH)  &
    +k(187)*n(idx_Ck)  &
    +k(63)*n(idx_CH2)

!d[CH2_dot]/d[H]
pd(13,5) =  &
    -k(63)*n(idx_CH2)

!d[H2O_dot]/d[H]
pd(16,5) =  &
    -k(79)*n(idx_H2O)  &
    +k(280)*n(idx_OH)  &
    +k(207)*n(idx_OH)

!d[O2_dot]/d[H]
pd(17,5) =  &
    -k(80)*n(idx_O2)

!d[H+_dot]/d[H]
pd(20,5) =  &
    +k(158)*n(idx_COj)  &
    +k(44)*n(idx_Oj)  &
    +k(48)*n(idx_Cj)  &
    -k(18)*n(idx_Hj)  &
    +k(20)*n(idx_H2j)  &
    +k(1)*n(idx_E)  &
    -k(19)*n(idx_Hj)  &
    +k(252)  &
    +k(8)*n(idx_HEj)  &
    +k(213)

!d[HE+_dot]/d[H]
pd(21,5) =  &
    -k(8)*n(idx_HEj)

!d[H2+_dot]/d[H]
pd(22,5) =  &
    +k(18)*n(idx_Hj)  &
    +k(86)*n(idx_H3j)  &
    -k(20)*n(idx_H2j)  &
    +k(19)*n(idx_Hj)

!d[C+_dot]/d[H]
pd(23,5) =  &
    -k(48)*n(idx_Cj)  &
    +k(91)*n(idx_CHj)  &
    -k(200)*n(idx_Cj)

!d[O+_dot]/d[H]
pd(24,5) =  &
    -k(44)*n(idx_Oj)

!d[H3+_dot]/d[H]
pd(27,5) =  &
    -k(86)*n(idx_H3j)

!d[CH+_dot]/d[H]
pd(28,5) =  &
    +k(200)*n(idx_Cj)  &
    +k(94)*n(idx_CH2j)  &
    -k(91)*n(idx_CHj)

!d[CH2+_dot]/d[H]
pd(29,5) =  &
    -k(94)*n(idx_CH2j)  &
    +k(97)*n(idx_CH3j)

!d[CO+_dot]/d[H]
pd(30,5) =  &
    -k(158)*n(idx_COj)

!d[CH3+_dot]/d[H]
pd(31,5) =  &
    -k(97)*n(idx_CH3j)

!d[Tgas_dot]/d[H]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(5)*1d-3
if(dnn>0.d0) then
nn(5) = n(5) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,5) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[HE]
pd(1,6) =  &
    +k(253)  &
    +k(214)  &
    +2.d0*k(4)*n(idx_E)  &
    -k(4)*n(idx_E)

!d[H_dot]/d[HE]
pd(5,6) =  &
    +k(9)*n(idx_Hj)  &
    +k(10)*n(idx_Hj)  &
    +2.d0*k(11)*n(idx_H2)  &
    -2.d0*k(34)*n(idx_H)*n(idx_H)

!d[HE_dot]/d[HE]
pd(6,6) =  &
    -k(10)*n(idx_Hj)  &
    -k(11)*n(idx_H2)  &
    -k(214)  &
    -k(9)*n(idx_Hj)  &
    -k(253)  &
    -k(4)*n(idx_E)  &
    +k(11)*n(idx_H2)  &
    -k(34)*n(idx_H)*n(idx_H)  &
    +k(34)*n(idx_H)*n(idx_H)

!d[H2_dot]/d[HE]
pd(7,6) =  &
    -k(11)*n(idx_H2)  &
    +k(34)*n(idx_H)*n(idx_H)

!d[H+_dot]/d[HE]
pd(20,6) =  &
    -k(9)*n(idx_Hj)  &
    -k(10)*n(idx_Hj)

!d[HE+_dot]/d[HE]
pd(21,6) =  &
    +k(9)*n(idx_Hj)  &
    +k(10)*n(idx_Hj)  &
    +k(253)  &
    +k(214)  &
    +k(4)*n(idx_E)

!d[Tgas_dot]/d[HE]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(6)*1d-3
if(dnn>0.d0) then
nn(6) = n(6) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,6) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H2]
pd(1,7) =  &
    +k(23)*n(idx_E)  &
    +k(270)  &
    +k(191)*n(idx_Ok)  &
    -k(23)*n(idx_E)  &
    +k(260)  &
    +k(229)  &
    +k(188)*n(idx_Ck)  &
    +k(218)

!d[H-_dot]/d[H2]
pd(2,7) =  &
    +k(259)

!d[C-_dot]/d[H2]
pd(3,7) =  &
    -k(188)*n(idx_Ck)

!d[O-_dot]/d[H2]
pd(4,7) =  &
    -k(191)*n(idx_Ok)

!d[H_dot]/d[H2]
pd(5,7) =  &
    +2.d0*k(23)*n(idx_E)  &
    +2.d0*k(258)  &
    +k(110)*n(idx_H2Oj)  &
    +k(90)*n(idx_Cj)  &
    +k(95)*n(idx_CH2j)  &
    +4.d0*k(33)*n(idx_H2)  &
    +2.d0*k(11)*n(idx_HE)  &
    +k(101)*n(idx_Oj)  &
    +k(22)*n(idx_Hj)  &
    +k(21)*n(idx_Hj)  &
    +k(92)*n(idx_CHj)  &
    +3.d0*k(24)*n(idx_H)  &
    +k(70)*n(idx_O)  &
    +k(109)*n(idx_OHj)  &
    +k(13)*n(idx_HEj)  &
    +k(85)*n(idx_H2j)  &
    +k(73)*n(idx_OH)  &
    +k(58)*n(idx_CH)  &
    +2.d0*k(14)*n(idx_HEj)  &
    +k(56)*n(idx_C)  &
    -2.d0*k(36)*n(idx_H)*n(idx_H)  &
    -k(24)*n(idx_H)  &
    +k(270)  &
    +2.d0*k(231)  &
    +k(229)  &
    +2.d0*k(193)*n(idx_Hj)

!d[HE_dot]/d[H2]
pd(6,7) =  &
    +k(12)*n(idx_HEj)  &
    -k(11)*n(idx_HE)  &
    +k(13)*n(idx_HEj)  &
    +k(11)*n(idx_HE)

!d[H2_dot]/d[H2]
pd(7,7) =  &
    -k(13)*n(idx_HEj)  &
    -k(270)  &
    -k(22)*n(idx_Hj)  &
    -k(260)  &
    +2.d0*k(33)*n(idx_H2)  &
    -k(218)  &
    -k(197)*n(idx_C)  &
    +2.d0*k(36)*n(idx_H)*n(idx_H)  &
    -k(70)*n(idx_O)  &
    -k(194)*n(idx_Hj)  &
    -k(259)  &
    -k(110)*n(idx_H2Oj)  &
    -k(188)*n(idx_Ck)  &
    -k(101)*n(idx_Oj)  &
    -k(191)*n(idx_Ok)  &
    -k(95)*n(idx_CH2j)  &
    -k(193)*n(idx_Hj)  &
    -k(201)*n(idx_Cj)  &
    -k(92)*n(idx_CHj)  &
    -k(11)*n(idx_HE)  &
    -k(14)*n(idx_HEj)  &
    -k(258)  &
    -k(12)*n(idx_HEj)  &
    -k(58)*n(idx_CH)  &
    -k(85)*n(idx_H2j)  &
    -k(53)*n(idx_HOCj)  &
    -k(36)*n(idx_H)*n(idx_H)  &
    -k(24)*n(idx_H)  &
    -k(73)*n(idx_OH)  &
    -k(109)*n(idx_OHj)  &
    -k(231)  &
    -k(90)*n(idx_Cj)  &
    -k(56)*n(idx_C)  &
    -k(23)*n(idx_E)  &
    -4.d0*k(33)*n(idx_H2)  &
    -k(229)  &
    -k(21)*n(idx_Hj)  &
    +k(53)*n(idx_HOCj)  &
    -k(81)*n(idx_O2)

!d[C_dot]/d[H2]
pd(8,7) =  &
    -k(197)*n(idx_C)  &
    -k(56)*n(idx_C)

!d[O_dot]/d[H2]
pd(9,7) =  &
    -k(70)*n(idx_O)

!d[OH_dot]/d[H2]
pd(10,7) =  &
    +2.d0*k(81)*n(idx_O2)  &
    -k(73)*n(idx_OH)  &
    +k(70)*n(idx_O)

!d[CH_dot]/d[H2]
pd(12,7) =  &
    -k(58)*n(idx_CH)  &
    +k(56)*n(idx_C)

!d[CH2_dot]/d[H2]
pd(13,7) =  &
    +k(58)*n(idx_CH)  &
    +k(197)*n(idx_C)  &
    +k(188)*n(idx_Ck)

!d[H2O_dot]/d[H2]
pd(16,7) =  &
    +k(191)*n(idx_Ok)  &
    +k(73)*n(idx_OH)

!d[O2_dot]/d[H2]
pd(17,7) =  &
    -k(81)*n(idx_O2)

!d[H+_dot]/d[H2]
pd(20,7) =  &
    +k(270)  &
    -k(22)*n(idx_Hj)  &
    -k(193)*n(idx_Hj)  &
    +k(259)  &
    +k(229)  &
    +k(193)*n(idx_Hj)  &
    -k(21)*n(idx_Hj)  &
    -k(194)*n(idx_Hj)  &
    +k(13)*n(idx_HEj)

!d[HE+_dot]/d[H2]
pd(21,7) =  &
    -k(13)*n(idx_HEj)  &
    -k(12)*n(idx_HEj)  &
    +k(14)*n(idx_HEj)  &
    -k(14)*n(idx_HEj)

!d[H2+_dot]/d[H2]
pd(22,7) =  &
    +k(22)*n(idx_Hj)  &
    +k(260)  &
    +k(21)*n(idx_Hj)  &
    -k(85)*n(idx_H2j)  &
    +k(12)*n(idx_HEj)  &
    +k(218)

!d[C+_dot]/d[H2]
pd(23,7) =  &
    -k(90)*n(idx_Cj)  &
    -k(201)*n(idx_Cj)

!d[O+_dot]/d[H2]
pd(24,7) =  &
    -k(101)*n(idx_Oj)

!d[HOC+_dot]/d[H2]
pd(25,7) =  &
    -k(53)*n(idx_HOCj)

!d[HCO+_dot]/d[H2]
pd(26,7) =  &
    +k(53)*n(idx_HOCj)

!d[H3+_dot]/d[H2]
pd(27,7) =  &
    +k(194)*n(idx_Hj)  &
    +k(85)*n(idx_H2j)

!d[CH+_dot]/d[H2]
pd(28,7) =  &
    -k(92)*n(idx_CHj)  &
    +k(90)*n(idx_Cj)

!d[CH2+_dot]/d[H2]
pd(29,7) =  &
    -k(95)*n(idx_CH2j)  &
    +k(92)*n(idx_CHj)  &
    +k(201)*n(idx_Cj)

!d[CH3+_dot]/d[H2]
pd(31,7) =  &
    +k(95)*n(idx_CH2j)

!d[OH+_dot]/d[H2]
pd(32,7) =  &
    +k(101)*n(idx_Oj)  &
    -k(109)*n(idx_OHj)

!d[H2O+_dot]/d[H2]
pd(33,7) =  &
    +k(109)*n(idx_OHj)  &
    -k(110)*n(idx_H2Oj)

!d[H3O+_dot]/d[H2]
pd(34,7) =  &
    +k(110)*n(idx_H2Oj)

!d[Tgas_dot]/d[H2]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(7)*1d-3
if(dnn>0.d0) then
nn(7) = n(7) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,7) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[C]
pd(1,8) =  &
    -k(195)*n(idx_E)  &
    +k(261)  &
    +k(192)*n(idx_Ok)  &
    +2.d0*k(42)*n(idx_E)  &
    +k(217)  &
    -k(42)*n(idx_E)  &
    +k(184)*n(idx_Hk)

!d[H-_dot]/d[C]
pd(2,8) =  &
    -k(184)*n(idx_Hk)

!d[C-_dot]/d[C]
pd(3,8) =  &
    +k(195)*n(idx_E)

!d[O-_dot]/d[C]
pd(4,8) =  &
    -k(192)*n(idx_Ok)

!d[H_dot]/d[C]
pd(5,8) =  &
    +k(74)*n(idx_OH)  &
    +k(87)*n(idx_H2j)  &
    +k(89)*n(idx_H3j)  &
    -k(196)*n(idx_H)  &
    +k(75)*n(idx_OH)  &
    +k(59)*n(idx_CH)  &
    +k(47)*n(idx_Hj)  &
    +k(56)*n(idx_H2)

!d[HE_dot]/d[C]
pd(6,8) =  &
    +k(51)*n(idx_HEj)  &
    +k(50)*n(idx_HEj)  &
    +k(49)*n(idx_HEj)

!d[H2_dot]/d[C]
pd(7,8) =  &
    -k(197)*n(idx_H2)  &
    +k(88)*n(idx_H3j)  &
    -k(56)*n(idx_H2)  &
    +k(117)*n(idx_H3Oj)

!d[C_dot]/d[C]
pd(8,8) =  &
    -k(59)*n(idx_CH)  &
    -k(87)*n(idx_H2j)  &
    -k(278)*n(idx_Oj)  &
    -k(196)*n(idx_H)  &
    -k(197)*n(idx_H2)  &
    -k(51)*n(idx_HEj)  &
    -k(42)*n(idx_E)  &
    -k(117)*n(idx_H3Oj)  &
    -k(89)*n(idx_H3j)  &
    -k(195)*n(idx_E)  &
    -k(274)*n(idx_O)  &
    -k(184)*n(idx_Hk)  &
    -4.d0*k(198)*n(idx_C)  &
    -k(277)*n(idx_Oj)  &
    -k(273)*n(idx_O)  &
    -k(217)  &
    -k(49)*n(idx_HEj)  &
    -k(88)*n(idx_H3j)  &
    -k(75)*n(idx_OH)  &
    -k(83)*n(idx_O2)  &
    -k(261)  &
    -k(127)*n(idx_HCOj)  &
    -k(50)*n(idx_HEj)  &
    -4.d0*k(272)*n(idx_C)  &
    -k(192)*n(idx_Ok)  &
    -4.d0*k(271)*n(idx_C)  &
    -k(199)*n(idx_O)  &
    -k(47)*n(idx_Hj)  &
    -k(122)*n(idx_O2j)  &
    -k(121)*n(idx_O2j)  &
    -k(74)*n(idx_OH)  &
    -k(56)*n(idx_H2)  &
    -k(82)*n(idx_O2)

!d[O_dot]/d[C]
pd(9,8) =  &
    -k(274)*n(idx_O)  &
    -k(199)*n(idx_O)  &
    +k(83)*n(idx_O2)  &
    +k(121)*n(idx_O2j)  &
    +k(82)*n(idx_O2)  &
    -k(273)*n(idx_O)

!d[OH_dot]/d[C]
pd(10,8) =  &
    -k(74)*n(idx_OH)  &
    -k(75)*n(idx_OH)

!d[CO_dot]/d[C]
pd(11,8) =  &
    +k(74)*n(idx_OH)  &
    +k(192)*n(idx_Ok)  &
    +k(82)*n(idx_O2)  &
    +k(274)*n(idx_O)  &
    +k(83)*n(idx_O2)  &
    +k(127)*n(idx_HCOj)  &
    +k(75)*n(idx_OH)  &
    +k(199)*n(idx_O)  &
    +k(273)*n(idx_O)

!d[CH_dot]/d[C]
pd(12,8) =  &
    +k(196)*n(idx_H)  &
    -k(59)*n(idx_CH)  &
    +k(184)*n(idx_Hk)  &
    +k(56)*n(idx_H2)

!d[CH2_dot]/d[C]
pd(13,8) =  &
    +k(197)*n(idx_H2)

!d[C2_dot]/d[C]
pd(14,8) =  &
    +2.d0*k(272)*n(idx_C)  &
    +k(59)*n(idx_CH)  &
    +2.d0*k(198)*n(idx_C)  &
    +2.d0*k(271)*n(idx_C)

!d[O2_dot]/d[C]
pd(17,8) =  &
    -k(83)*n(idx_O2)  &
    -k(82)*n(idx_O2)  &
    +k(122)*n(idx_O2j)

!d[H+_dot]/d[C]
pd(20,8) =  &
    -k(47)*n(idx_Hj)

!d[HE+_dot]/d[C]
pd(21,8) =  &
    -k(49)*n(idx_HEj)  &
    -k(51)*n(idx_HEj)  &
    -k(50)*n(idx_HEj)

!d[H2+_dot]/d[C]
pd(22,8) =  &
    -k(87)*n(idx_H2j)

!d[C+_dot]/d[C]
pd(23,8) =  &
    +k(261)  &
    +k(51)*n(idx_HEj)  &
    +k(42)*n(idx_E)  &
    +k(49)*n(idx_HEj)  &
    +k(50)*n(idx_HEj)  &
    +k(217)  &
    +k(122)*n(idx_O2j)  &
    +k(47)*n(idx_Hj)

!d[O+_dot]/d[C]
pd(24,8) =  &
    -k(278)*n(idx_Oj)  &
    -k(277)*n(idx_Oj)

!d[HCO+_dot]/d[C]
pd(26,8) =  &
    -k(127)*n(idx_HCOj)  &
    +k(117)*n(idx_H3Oj)

!d[H3+_dot]/d[C]
pd(27,8) =  &
    -k(88)*n(idx_H3j)  &
    -k(89)*n(idx_H3j)

!d[CH+_dot]/d[C]
pd(28,8) =  &
    +k(127)*n(idx_HCOj)  &
    +k(88)*n(idx_H3j)  &
    +k(87)*n(idx_H2j)

!d[CH2+_dot]/d[C]
pd(29,8) =  &
    +k(89)*n(idx_H3j)

!d[CO+_dot]/d[C]
pd(30,8) =  &
    +k(121)*n(idx_O2j)  &
    +k(277)*n(idx_Oj)  &
    +k(278)*n(idx_Oj)

!d[H3O+_dot]/d[C]
pd(34,8) =  &
    -k(117)*n(idx_H3Oj)

!d[O2+_dot]/d[C]
pd(35,8) =  &
    -k(121)*n(idx_O2j)  &
    -k(122)*n(idx_O2j)

!d[Tgas_dot]/d[C]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(8)*1d-3
if(dnn>0.d0) then
nn(8) = n(8) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,8) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[O]
pd(1,9) =  &
    -k(43)*n(idx_E)  &
    -k(204)*n(idx_E)  &
    +k(189)*n(idx_Ck)  &
    +2.d0*k(43)*n(idx_E)  &
    +k(216)  &
    +k(185)*n(idx_Hk)  &
    +k(61)*n(idx_CH)  &
    +k(254)

!d[H-_dot]/d[O]
pd(2,9) =  &
    -k(185)*n(idx_Hk)

!d[C-_dot]/d[O]
pd(3,9) =  &
    -k(189)*n(idx_Ck)

!d[O-_dot]/d[O]
pd(4,9) =  &
    +k(204)*n(idx_E)

!d[H_dot]/d[O]
pd(5,9) =  &
    +k(70)*n(idx_H2)  &
    +k(104)*n(idx_H3j)  &
    +k(76)*n(idx_OH)  &
    +k(60)*n(idx_CH)  &
    +k(96)*n(idx_CH2j)  &
    +k(93)*n(idx_CHj)  &
    +k(66)*n(idx_CH2)  &
    +k(102)*n(idx_H2j)  &
    -k(279)*n(idx_H)  &
    -k(205)*n(idx_H)  &
    +k(45)*n(idx_Hj)  &
    +2.d0*k(64)*n(idx_CH2)  &
    +k(77)*n(idx_OH)

!d[HE_dot]/d[O]
pd(6,9) =  &
    +k(46)*n(idx_HEj)

!d[H2_dot]/d[O]
pd(7,9) =  &
    -k(70)*n(idx_H2)  &
    +k(98)*n(idx_CH3j)  &
    +k(99)*n(idx_CH3j)  &
    +k(65)*n(idx_CH2)  &
    +k(103)*n(idx_H3j)

!d[C_dot]/d[O]
pd(8,9) =  &
    -k(199)*n(idx_C)  &
    -k(274)*n(idx_C)  &
    +k(62)*n(idx_CH)  &
    +k(69)*n(idx_C2)  &
    +k(68)*n(idx_C2)  &
    -k(273)*n(idx_C)

!d[O_dot]/d[O]
pd(9,9) =  &
    -k(199)*n(idx_C)  &
    -k(43)*n(idx_E)  &
    -k(216)  &
    -k(46)*n(idx_HEj)  &
    -k(202)*n(idx_Cj)  &
    -4.d0*k(281)*n(idx_O)  &
    -k(99)*n(idx_CH3j)  &
    -k(102)*n(idx_H2j)  &
    -k(77)*n(idx_OH)  &
    -k(205)*n(idx_H)  &
    -k(279)*n(idx_H)  &
    -k(69)*n(idx_C2)  &
    -k(65)*n(idx_CH2)  &
    -k(103)*n(idx_H3j)  &
    -k(70)*n(idx_H2)  &
    -k(45)*n(idx_Hj)  &
    -k(276)*n(idx_Cj)  &
    -k(185)*n(idx_Hk)  &
    -k(96)*n(idx_CH2j)  &
    -k(273)*n(idx_C)  &
    -k(203)*n(idx_Cj)  &
    -k(274)*n(idx_C)  &
    -k(60)*n(idx_CH)  &
    -k(62)*n(idx_CH)  &
    -k(66)*n(idx_CH2)  &
    -k(61)*n(idx_CH)  &
    -k(189)*n(idx_Ck)  &
    -k(98)*n(idx_CH3j)  &
    -k(67)*n(idx_CH2)  &
    -4.d0*k(206)*n(idx_O)  &
    -k(104)*n(idx_H3j)  &
    -k(275)*n(idx_Cj)  &
    -k(68)*n(idx_C2)  &
    -k(93)*n(idx_CHj)  &
    -k(204)*n(idx_E)  &
    -k(76)*n(idx_OH)  &
    -k(64)*n(idx_CH2)  &
    -k(254)

!d[OH_dot]/d[O]
pd(10,9) =  &
    +k(62)*n(idx_CH)  &
    +k(279)*n(idx_H)  &
    +k(70)*n(idx_H2)  &
    -k(76)*n(idx_OH)  &
    -k(77)*n(idx_OH)  &
    +k(185)*n(idx_Hk)  &
    +k(67)*n(idx_CH2)  &
    +k(205)*n(idx_H)

!d[CO_dot]/d[O]
pd(11,9) =  &
    +k(274)*n(idx_C)  &
    +k(60)*n(idx_CH)  &
    +k(189)*n(idx_Ck)  &
    +k(69)*n(idx_C2)  &
    +k(65)*n(idx_CH2)  &
    +k(273)*n(idx_C)  &
    +k(64)*n(idx_CH2)  &
    +k(199)*n(idx_C)  &
    +k(68)*n(idx_C2)

!d[CH_dot]/d[O]
pd(12,9) =  &
    -k(60)*n(idx_CH)  &
    -k(62)*n(idx_CH)  &
    -k(61)*n(idx_CH)  &
    +k(67)*n(idx_CH2)

!d[CH2_dot]/d[O]
pd(13,9) =  &
    -k(64)*n(idx_CH2)  &
    -k(65)*n(idx_CH2)  &
    -k(66)*n(idx_CH2)  &
    -k(67)*n(idx_CH2)

!d[C2_dot]/d[O]
pd(14,9) =  &
    -k(69)*n(idx_C2)  &
    -k(68)*n(idx_C2)

!d[HCO_dot]/d[O]
pd(15,9) =  &
    +k(66)*n(idx_CH2)

!d[O2_dot]/d[O]
pd(17,9) =  &
    +k(76)*n(idx_OH)  &
    +2.d0*k(281)*n(idx_O)  &
    +k(77)*n(idx_OH)  &
    +2.d0*k(206)*n(idx_O)

!d[H+_dot]/d[O]
pd(20,9) =  &
    -k(45)*n(idx_Hj)

!d[HE+_dot]/d[O]
pd(21,9) =  &
    -k(46)*n(idx_HEj)

!d[H2+_dot]/d[O]
pd(22,9) =  &
    -k(102)*n(idx_H2j)

!d[C+_dot]/d[O]
pd(23,9) =  &
    -k(275)*n(idx_Cj)  &
    -k(203)*n(idx_Cj)  &
    -k(202)*n(idx_Cj)  &
    -k(276)*n(idx_Cj)

!d[O+_dot]/d[O]
pd(24,9) =  &
    +k(46)*n(idx_HEj)  &
    +k(216)  &
    +k(45)*n(idx_Hj)  &
    +k(254)  &
    +k(43)*n(idx_E)

!d[HOC+_dot]/d[O]
pd(25,9) =  &
    +k(98)*n(idx_CH3j)

!d[HCO+_dot]/d[O]
pd(26,9) =  &
    +k(61)*n(idx_CH)  &
    +k(99)*n(idx_CH3j)  &
    +k(96)*n(idx_CH2j)

!d[H3+_dot]/d[O]
pd(27,9) =  &
    -k(104)*n(idx_H3j)  &
    -k(103)*n(idx_H3j)

!d[CH+_dot]/d[O]
pd(28,9) =  &
    -k(93)*n(idx_CHj)

!d[CH2+_dot]/d[O]
pd(29,9) =  &
    -k(96)*n(idx_CH2j)

!d[CO+_dot]/d[O]
pd(30,9) =  &
    +k(275)*n(idx_Cj)  &
    +k(203)*n(idx_Cj)  &
    +k(93)*n(idx_CHj)  &
    +k(202)*n(idx_Cj)  &
    +k(276)*n(idx_Cj)

!d[CH3+_dot]/d[O]
pd(31,9) =  &
    -k(98)*n(idx_CH3j)  &
    -k(99)*n(idx_CH3j)

!d[OH+_dot]/d[O]
pd(32,9) =  &
    +k(103)*n(idx_H3j)  &
    +k(102)*n(idx_H2j)

!d[H2O+_dot]/d[O]
pd(33,9) =  &
    +k(104)*n(idx_H3j)

!d[Tgas_dot]/d[O]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(9)*1d-3
if(dnn>0.d0) then
nn(9) = n(9) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,9) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[OH]
pd(1,10) =  &
    +k(224)  &
    +k(186)*n(idx_Hk)

!d[H-_dot]/d[OH]
pd(2,10) =  &
    -k(186)*n(idx_Hk)

!d[H_dot]/d[OH]
pd(5,10) =  &
    -k(280)*n(idx_H)  &
    -k(207)*n(idx_H)  &
    +k(144)*n(idx_HEj)  &
    +k(223)  &
    +k(77)*n(idx_O)  &
    +k(141)*n(idx_Hj)  &
    +k(143)*n(idx_HEj)  &
    +k(142)*n(idx_Hj)  &
    -k(72)*n(idx_H)  &
    +k(76)*n(idx_O)  &
    -k(71)*n(idx_H)  &
    +k(107)*n(idx_Cj)  &
    +k(75)*n(idx_C)  &
    +k(108)*n(idx_Cj)  &
    +2.d0*k(52)*n(idx_H)  &
    +k(265)  &
    -k(52)*n(idx_H)  &
    +k(73)*n(idx_H2)  &
    +k(74)*n(idx_C)

!d[HE_dot]/d[OH]
pd(6,10) =  &
    +k(143)*n(idx_HEj)  &
    +k(144)*n(idx_HEj)

!d[H2_dot]/d[OH]
pd(7,10) =  &
    +k(72)*n(idx_H)  &
    +k(106)*n(idx_H3j)  &
    +k(105)*n(idx_H3j)  &
    -k(73)*n(idx_H2)  &
    +k(71)*n(idx_H)

!d[C_dot]/d[OH]
pd(8,10) =  &
    -k(75)*n(idx_C)  &
    -k(74)*n(idx_C)

!d[O_dot]/d[OH]
pd(9,10) =  &
    -k(77)*n(idx_O)  &
    +k(71)*n(idx_H)  &
    +k(223)  &
    -k(76)*n(idx_O)  &
    +k(72)*n(idx_H)  &
    +k(52)*n(idx_H)  &
    +k(265)  &
    +2.d0*k(78)*n(idx_OH)

!d[OH_dot]/d[OH]
pd(10,10) =  &
    -k(280)*n(idx_H)  &
    -k(143)*n(idx_HEj)  &
    -k(77)*n(idx_O)  &
    -k(72)*n(idx_H)  &
    -k(142)*n(idx_Hj)  &
    -k(186)*n(idx_Hk)  &
    -k(224)  &
    -k(141)*n(idx_Hj)  &
    -k(76)*n(idx_O)  &
    -k(75)*n(idx_C)  &
    -4.d0*k(78)*n(idx_OH)  &
    -k(144)*n(idx_HEj)  &
    -k(73)*n(idx_H2)  &
    -k(106)*n(idx_H3j)  &
    -k(52)*n(idx_H)  &
    -k(207)*n(idx_H)  &
    -k(74)*n(idx_C)  &
    -k(223)  &
    -k(107)*n(idx_Cj)  &
    -k(108)*n(idx_Cj)  &
    -k(71)*n(idx_H)  &
    -k(265)  &
    -k(105)*n(idx_H3j)

!d[CO_dot]/d[OH]
pd(11,10) =  &
    +k(75)*n(idx_C)  &
    +k(74)*n(idx_C)

!d[H2O_dot]/d[OH]
pd(16,10) =  &
    +k(73)*n(idx_H2)  &
    +k(207)*n(idx_H)  &
    +2.d0*k(78)*n(idx_OH)  &
    +k(280)*n(idx_H)  &
    +k(186)*n(idx_Hk)

!d[O2_dot]/d[OH]
pd(17,10) =  &
    +k(77)*n(idx_O)  &
    +k(76)*n(idx_O)

!d[H+_dot]/d[OH]
pd(20,10) =  &
    -k(142)*n(idx_Hj)  &
    -k(141)*n(idx_Hj)

!d[HE+_dot]/d[OH]
pd(21,10) =  &
    -k(144)*n(idx_HEj)  &
    -k(143)*n(idx_HEj)

!d[C+_dot]/d[OH]
pd(23,10) =  &
    -k(107)*n(idx_Cj)  &
    -k(108)*n(idx_Cj)

!d[O+_dot]/d[OH]
pd(24,10) =  &
    +k(143)*n(idx_HEj)  &
    +k(144)*n(idx_HEj)

!d[H3+_dot]/d[OH]
pd(27,10) =  &
    -k(105)*n(idx_H3j)  &
    -k(106)*n(idx_H3j)

!d[CO+_dot]/d[OH]
pd(30,10) =  &
    +k(108)*n(idx_Cj)  &
    +k(107)*n(idx_Cj)

!d[OH+_dot]/d[OH]
pd(32,10) =  &
    +k(224)  &
    +k(141)*n(idx_Hj)  &
    +k(142)*n(idx_Hj)

!d[H2O+_dot]/d[OH]
pd(33,10) =  &
    +k(106)*n(idx_H3j)  &
    +k(105)*n(idx_H3j)

!d[Tgas_dot]/d[OH]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(10)*1d-3
if(dnn>0.d0) then
nn(10) = n(10) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,10) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[CO]
pd(1,11) =  &
    +k(256)

!d[H_dot]/d[CO]
pd(5,11) =  &
    -k(84)*n(idx_H)

!d[HE_dot]/d[CO]
pd(6,11) =  &
    +k(157)*n(idx_HEj)  &
    +k(156)*n(idx_HEj)

!d[H2_dot]/d[CO]
pd(7,11) =  &
    +k(123)*n(idx_H3j)  &
    +k(124)*n(idx_H3j)  &
    +k(125)*n(idx_H3j)  &
    +k(126)*n(idx_H3j)

!d[C_dot]/d[CO]
pd(8,11) =  &
    +k(255)  &
    +k(84)*n(idx_H)  &
    +k(157)*n(idx_HEj)  &
    +k(230)

!d[O_dot]/d[CO]
pd(9,11) =  &
    +k(255)  &
    +k(156)*n(idx_HEj)  &
    +k(230)

!d[OH_dot]/d[CO]
pd(10,11) =  &
    +k(84)*n(idx_H)

!d[CO_dot]/d[CO]
pd(11,11) =  &
    +k(55)*n(idx_HOCj)  &
    -k(125)*n(idx_H3j)  &
    -k(157)*n(idx_HEj)  &
    -k(54)*n(idx_HOCj)  &
    -k(55)*n(idx_HOCj)  &
    -k(84)*n(idx_H)  &
    -k(230)  &
    -k(123)*n(idx_H3j)  &
    -k(156)*n(idx_HEj)  &
    +k(54)*n(idx_HOCj)  &
    -k(124)*n(idx_H3j)  &
    -k(256)  &
    -k(126)*n(idx_H3j)  &
    -k(255)

!d[HE+_dot]/d[CO]
pd(21,11) =  &
    -k(157)*n(idx_HEj)  &
    -k(156)*n(idx_HEj)

!d[C+_dot]/d[CO]
pd(23,11) =  &
    +k(156)*n(idx_HEj)

!d[O+_dot]/d[CO]
pd(24,11) =  &
    +k(157)*n(idx_HEj)

!d[HOC+_dot]/d[CO]
pd(25,11) =  &
    -k(55)*n(idx_HOCj)  &
    +k(125)*n(idx_H3j)  &
    +k(126)*n(idx_H3j)  &
    -k(54)*n(idx_HOCj)

!d[HCO+_dot]/d[CO]
pd(26,11) =  &
    +k(54)*n(idx_HOCj)  &
    +k(55)*n(idx_HOCj)  &
    +k(123)*n(idx_H3j)  &
    +k(124)*n(idx_H3j)

!d[H3+_dot]/d[CO]
pd(27,11) =  &
    -k(126)*n(idx_H3j)  &
    -k(124)*n(idx_H3j)  &
    -k(125)*n(idx_H3j)  &
    -k(123)*n(idx_H3j)

!d[CO+_dot]/d[CO]
pd(30,11) =  &
    +k(256)

!d[Tgas_dot]/d[CO]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(11)*1d-3
if(dnn>0.d0) then
nn(11) = n(11) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,11) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[CH]
pd(1,12) =  &
    +k(61)*n(idx_O)  &
    +k(221)

!d[H_dot]/d[CH]
pd(5,12) =  &
    +k(131)*n(idx_Hj)  &
    +k(60)*n(idx_O)  &
    +k(220)  &
    -k(57)*n(idx_H)  &
    +k(262)  &
    +k(58)*n(idx_H2)  &
    +k(130)*n(idx_Hj)  &
    +k(59)*n(idx_C)

!d[H2_dot]/d[CH]
pd(7,12) =  &
    +k(57)*n(idx_H)  &
    -k(58)*n(idx_H2)

!d[C_dot]/d[CH]
pd(8,12) =  &
    -k(59)*n(idx_C)  &
    +k(57)*n(idx_H)  &
    +k(220)  &
    +k(62)*n(idx_O)  &
    +k(262)

!d[O_dot]/d[CH]
pd(9,12) =  &
    -k(60)*n(idx_O)  &
    -k(62)*n(idx_O)  &
    -k(61)*n(idx_O)

!d[OH_dot]/d[CH]
pd(10,12) =  &
    +k(62)*n(idx_O)

!d[CO_dot]/d[CH]
pd(11,12) =  &
    +k(60)*n(idx_O)

!d[CH_dot]/d[CH]
pd(12,12) =  &
    -k(60)*n(idx_O)  &
    -k(262)  &
    -k(131)*n(idx_Hj)  &
    -k(130)*n(idx_Hj)  &
    -k(57)*n(idx_H)  &
    -k(220)  &
    -k(61)*n(idx_O)  &
    -k(58)*n(idx_H2)  &
    -k(59)*n(idx_C)  &
    -k(62)*n(idx_O)  &
    -k(221)

!d[CH2_dot]/d[CH]
pd(13,12) =  &
    +k(58)*n(idx_H2)

!d[C2_dot]/d[CH]
pd(14,12) =  &
    +k(59)*n(idx_C)

!d[H+_dot]/d[CH]
pd(20,12) =  &
    -k(130)*n(idx_Hj)  &
    -k(131)*n(idx_Hj)

!d[HCO+_dot]/d[CH]
pd(26,12) =  &
    +k(61)*n(idx_O)

!d[CH+_dot]/d[CH]
pd(28,12) =  &
    +k(130)*n(idx_Hj)  &
    +k(131)*n(idx_Hj)  &
    +k(221)

!d[Tgas_dot]/d[CH]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(12)*1d-3
if(dnn>0.d0) then
nn(12) = n(12) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,12) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[CH2]
pd(1,13) =  &
    +k(266)  &
    +k(238)

!d[H_dot]/d[CH2]
pd(5,13) =  &
    +2.d0*k(64)*n(idx_O)  &
    +k(138)*n(idx_HEj)  &
    +k(66)*n(idx_O)  &
    -k(63)*n(idx_H)  &
    +k(139)*n(idx_HEj)  &
    +k(237)  &
    +k(134)*n(idx_Hj)  &
    +k(135)*n(idx_Hj)

!d[HE_dot]/d[CH2]
pd(6,13) =  &
    +k(139)*n(idx_HEj)  &
    +k(138)*n(idx_HEj)  &
    +k(137)*n(idx_HEj)  &
    +k(136)*n(idx_HEj)

!d[H2_dot]/d[CH2]
pd(7,13) =  &
    +k(133)*n(idx_Hj)  &
    +k(137)*n(idx_HEj)  &
    +k(136)*n(idx_HEj)  &
    +k(63)*n(idx_H)  &
    +k(132)*n(idx_Hj)  &
    +k(65)*n(idx_O)

!d[O_dot]/d[CH2]
pd(9,13) =  &
    -k(66)*n(idx_O)  &
    -k(65)*n(idx_O)  &
    -k(67)*n(idx_O)  &
    -k(64)*n(idx_O)

!d[OH_dot]/d[CH2]
pd(10,13) =  &
    +k(67)*n(idx_O)

!d[CO_dot]/d[CH2]
pd(11,13) =  &
    +k(64)*n(idx_O)  &
    +k(65)*n(idx_O)

!d[CH_dot]/d[CH2]
pd(12,13) =  &
    +k(63)*n(idx_H)  &
    +k(67)*n(idx_O)  &
    +k(237)

!d[CH2_dot]/d[CH2]
pd(13,13) =  &
    -k(66)*n(idx_O)  &
    -k(266)  &
    -k(237)  &
    -k(137)*n(idx_HEj)  &
    -k(136)*n(idx_HEj)  &
    -k(64)*n(idx_O)  &
    -k(134)*n(idx_Hj)  &
    -k(132)*n(idx_Hj)  &
    -k(67)*n(idx_O)  &
    -k(238)  &
    -k(63)*n(idx_H)  &
    -k(65)*n(idx_O)  &
    -k(135)*n(idx_Hj)  &
    -k(133)*n(idx_Hj)  &
    -k(139)*n(idx_HEj)  &
    -k(138)*n(idx_HEj)

!d[HCO_dot]/d[CH2]
pd(15,13) =  &
    +k(66)*n(idx_O)

!d[H+_dot]/d[CH2]
pd(20,13) =  &
    -k(135)*n(idx_Hj)  &
    -k(132)*n(idx_Hj)  &
    -k(133)*n(idx_Hj)  &
    -k(134)*n(idx_Hj)

!d[HE+_dot]/d[CH2]
pd(21,13) =  &
    -k(139)*n(idx_HEj)  &
    -k(138)*n(idx_HEj)  &
    -k(137)*n(idx_HEj)  &
    -k(136)*n(idx_HEj)

!d[C+_dot]/d[CH2]
pd(23,13) =  &
    +k(137)*n(idx_HEj)  &
    +k(136)*n(idx_HEj)

!d[CH+_dot]/d[CH2]
pd(28,13) =  &
    +k(139)*n(idx_HEj)  &
    +k(138)*n(idx_HEj)  &
    +k(132)*n(idx_Hj)  &
    +k(133)*n(idx_Hj)

!d[CH2+_dot]/d[CH2]
pd(29,13) =  &
    +k(266)  &
    +k(134)*n(idx_Hj)  &
    +k(135)*n(idx_Hj)  &
    +k(238)

!d[Tgas_dot]/d[CH2]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(13)*1d-3
if(dnn>0.d0) then
nn(13) = n(13) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,13) = (dn1-dn0)/dnn
end if

!d[HE_dot]/d[C2]
pd(6,14) =  &
    +k(140)*n(idx_HEj)

!d[C_dot]/d[C2]
pd(8,14) =  &
    +2.d0*k(257)  &
    +2.d0*k(222)  &
    +k(140)*n(idx_HEj)  &
    +k(69)*n(idx_O)  &
    +k(68)*n(idx_O)  &
    +k(100)*n(idx_Oj)

!d[O_dot]/d[C2]
pd(9,14) =  &
    -k(68)*n(idx_O)  &
    -k(69)*n(idx_O)

!d[CO_dot]/d[C2]
pd(11,14) =  &
    +k(69)*n(idx_O)  &
    +k(68)*n(idx_O)

!d[C2_dot]/d[C2]
pd(14,14) =  &
    -k(222)  &
    -k(68)*n(idx_O)  &
    -k(140)*n(idx_HEj)  &
    -k(100)*n(idx_Oj)  &
    -k(69)*n(idx_O)  &
    -k(257)

!d[HE+_dot]/d[C2]
pd(21,14) =  &
    -k(140)*n(idx_HEj)

!d[C+_dot]/d[C2]
pd(23,14) =  &
    +k(140)*n(idx_HEj)

!d[O+_dot]/d[C2]
pd(24,14) =  &
    -k(100)*n(idx_Oj)

!d[CO+_dot]/d[C2]
pd(30,14) =  &
    +k(100)*n(idx_Oj)

!d[Tgas_dot]/d[C2]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(14)*1d-3
if(dnn>0.d0) then
nn(14) = n(14) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,14) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[HCO]
pd(1,15) =  &
    +k(269)

!d[H_dot]/d[HCO]
pd(5,15) =  &
    +k(268)

!d[CO_dot]/d[HCO]
pd(11,15) =  &
    +k(268)

!d[HCO_dot]/d[HCO]
pd(15,15) =  &
    -k(269)  &
    -k(268)

!d[HCO+_dot]/d[HCO]
pd(26,15) =  &
    +k(269)

!d[Tgas_dot]/d[HCO]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(15)*1d-3
if(dnn>0.d0) then
nn(15) = n(15) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,15) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H2O]
pd(1,16) =  &
    +k(226)

!d[H_dot]/d[H2O]
pd(5,16) =  &
    +k(113)*n(idx_Cj)  &
    +k(115)*n(idx_Cj)  &
    +k(149)*n(idx_HEj)  &
    -k(79)*n(idx_H)  &
    +k(150)*n(idx_HEj)  &
    +k(114)*n(idx_Cj)  &
    +k(145)*n(idx_Hj)  &
    +k(225)  &
    +k(146)*n(idx_Hj)  &
    +k(267)

!d[HE_dot]/d[H2O]
pd(6,16) =  &
    +k(148)*n(idx_HEj)  &
    +k(149)*n(idx_HEj)  &
    +k(150)*n(idx_HEj)  &
    +k(151)*n(idx_HEj)  &
    +k(147)*n(idx_HEj)  &
    +k(152)*n(idx_HEj)

!d[H2_dot]/d[H2O]
pd(7,16) =  &
    +k(111)*n(idx_H3j)  &
    +k(112)*n(idx_H3j)  &
    +k(79)*n(idx_H)

!d[C_dot]/d[H2O]
pd(8,16) =  &
    +k(116)*n(idx_Cj)

!d[OH_dot]/d[H2O]
pd(10,16) =  &
    +k(148)*n(idx_HEj)  &
    +k(267)  &
    +k(147)*n(idx_HEj)  &
    +k(225)  &
    +k(79)*n(idx_H)

!d[CO_dot]/d[H2O]
pd(11,16) =  &
    +k(128)*n(idx_HCOj)  &
    +k(129)*n(idx_HCOj)

!d[H2O_dot]/d[H2O]
pd(16,16) =  &
    -k(129)*n(idx_HCOj)  &
    -k(148)*n(idx_HEj)  &
    -k(116)*n(idx_Cj)  &
    -k(115)*n(idx_Cj)  &
    -k(225)  &
    -k(79)*n(idx_H)  &
    -k(145)*n(idx_Hj)  &
    -k(267)  &
    -k(151)*n(idx_HEj)  &
    -k(150)*n(idx_HEj)  &
    -k(152)*n(idx_HEj)  &
    -k(128)*n(idx_HCOj)  &
    -k(147)*n(idx_HEj)  &
    -k(112)*n(idx_H3j)  &
    -k(111)*n(idx_H3j)  &
    -k(149)*n(idx_HEj)  &
    -k(226)  &
    -k(146)*n(idx_Hj)  &
    -k(114)*n(idx_Cj)  &
    -k(113)*n(idx_Cj)

!d[H+_dot]/d[H2O]
pd(20,16) =  &
    +k(148)*n(idx_HEj)  &
    -k(146)*n(idx_Hj)  &
    -k(145)*n(idx_Hj)  &
    +k(147)*n(idx_HEj)

!d[HE+_dot]/d[H2O]
pd(21,16) =  &
    -k(152)*n(idx_HEj)  &
    -k(151)*n(idx_HEj)  &
    -k(150)*n(idx_HEj)  &
    -k(148)*n(idx_HEj)  &
    -k(147)*n(idx_HEj)  &
    -k(149)*n(idx_HEj)

!d[C+_dot]/d[H2O]
pd(23,16) =  &
    -k(116)*n(idx_Cj)  &
    -k(115)*n(idx_Cj)  &
    -k(114)*n(idx_Cj)  &
    -k(113)*n(idx_Cj)

!d[HOC+_dot]/d[H2O]
pd(25,16) =  &
    +k(113)*n(idx_Cj)

!d[HCO+_dot]/d[H2O]
pd(26,16) =  &
    -k(129)*n(idx_HCOj)  &
    +k(115)*n(idx_Cj)  &
    +k(114)*n(idx_Cj)  &
    -k(128)*n(idx_HCOj)

!d[H3+_dot]/d[H2O]
pd(27,16) =  &
    -k(112)*n(idx_H3j)  &
    -k(111)*n(idx_H3j)

!d[OH+_dot]/d[H2O]
pd(32,16) =  &
    +k(149)*n(idx_HEj)  &
    +k(150)*n(idx_HEj)

!d[H2O+_dot]/d[H2O]
pd(33,16) =  &
    +k(116)*n(idx_Cj)  &
    +k(151)*n(idx_HEj)  &
    +k(145)*n(idx_Hj)  &
    +k(146)*n(idx_Hj)  &
    +k(152)*n(idx_HEj)  &
    +k(226)

!d[H3O+_dot]/d[H2O]
pd(34,16) =  &
    +k(111)*n(idx_H3j)  &
    +k(128)*n(idx_HCOj)  &
    +k(112)*n(idx_H3j)  &
    +k(129)*n(idx_HCOj)

!d[Tgas_dot]/d[H2O]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(16)*1d-3
if(dnn>0.d0) then
nn(16) = n(16) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,16) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[O2]
pd(1,17) =  &
    +k(264)  &
    +k(227)

!d[H_dot]/d[O2]
pd(5,17) =  &
    +k(153)*n(idx_Hj)  &
    -k(80)*n(idx_H)

!d[HE_dot]/d[O2]
pd(6,17) =  &
    +k(155)*n(idx_HEj)  &
    +k(154)*n(idx_HEj)

!d[H2_dot]/d[O2]
pd(7,17) =  &
    -k(81)*n(idx_H2)

!d[C_dot]/d[O2]
pd(8,17) =  &
    -k(83)*n(idx_C)  &
    -k(82)*n(idx_C)

!d[O_dot]/d[O2]
pd(9,17) =  &
    +2.d0*k(263)  &
    +k(155)*n(idx_HEj)  &
    +k(118)*n(idx_Cj)  &
    +2.d0*k(228)  &
    +k(83)*n(idx_C)  &
    +k(82)*n(idx_C)  &
    +k(80)*n(idx_H)

!d[OH_dot]/d[O2]
pd(10,17) =  &
    +2.d0*k(81)*n(idx_H2)  &
    +k(120)*n(idx_CH2j)  &
    +k(80)*n(idx_H)

!d[CO_dot]/d[O2]
pd(11,17) =  &
    +k(82)*n(idx_C)  &
    +k(83)*n(idx_C)  &
    +k(119)*n(idx_Cj)

!d[O2_dot]/d[O2]
pd(17,17) =  &
    -k(118)*n(idx_Cj)  &
    -k(263)  &
    -k(153)*n(idx_Hj)  &
    -k(264)  &
    -k(82)*n(idx_C)  &
    -k(119)*n(idx_Cj)  &
    -k(83)*n(idx_C)  &
    -k(80)*n(idx_H)  &
    -k(154)*n(idx_HEj)  &
    -k(81)*n(idx_H2)  &
    -k(228)  &
    -k(227)  &
    -k(155)*n(idx_HEj)  &
    -k(120)*n(idx_CH2j)

!d[H+_dot]/d[O2]
pd(20,17) =  &
    -k(153)*n(idx_Hj)

!d[HE+_dot]/d[O2]
pd(21,17) =  &
    -k(155)*n(idx_HEj)  &
    -k(154)*n(idx_HEj)

!d[C+_dot]/d[O2]
pd(23,17) =  &
    -k(119)*n(idx_Cj)  &
    -k(118)*n(idx_Cj)

!d[O+_dot]/d[O2]
pd(24,17) =  &
    +k(155)*n(idx_HEj)  &
    +k(119)*n(idx_Cj)

!d[HCO+_dot]/d[O2]
pd(26,17) =  &
    +k(120)*n(idx_CH2j)

!d[CH2+_dot]/d[O2]
pd(29,17) =  &
    -k(120)*n(idx_CH2j)

!d[CO+_dot]/d[O2]
pd(30,17) =  &
    +k(118)*n(idx_Cj)

!d[O2+_dot]/d[O2]
pd(35,17) =  &
    +k(153)*n(idx_Hj)  &
    +k(264)  &
    +k(227)  &
    +k(154)*n(idx_HEj)

!d[Tgas_dot]/d[O2]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(17)*1d-3
if(dnn>0.d0) then
nn(17) = n(17) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,17) = (dn1-dn0)/dnn
end if

!d[Tgas_dot]/d[CO_total]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(18)*1d-3
if(dnn>0.d0) then
nn(18) = n(18) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,18) = (dn1-dn0)/dnn
end if

!d[Tgas_dot]/d[H2O_total]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(19)*1d-3
if(dnn>0.d0) then
nn(19) = n(19) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,19) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H+]
pd(1,20) =  &
    +k(29)*n(idx_Hk)  &
    -k(3)*n(idx_E)  &
    -k(2)*n(idx_E)

!d[H-_dot]/d[H+]
pd(2,20) =  &
    -k(29)*n(idx_Hk)  &
    -k(28)*n(idx_Hk)

!d[C-_dot]/d[H+]
pd(3,20) =  &
    -k(159)*n(idx_Ck)

!d[O-_dot]/d[H+]
pd(4,20) =  &
    -k(160)*n(idx_Ok)

!d[H_dot]/d[H+]
pd(5,20) =  &
    +k(145)*n(idx_H2O)  &
    +k(22)*n(idx_H2)  &
    +k(3)*n(idx_E)  &
    +k(159)*n(idx_Ck)  &
    +k(9)*n(idx_HE)  &
    +k(135)*n(idx_CH2)  &
    +k(153)*n(idx_O2)  &
    +k(21)*n(idx_H2)  &
    -k(19)*n(idx_H)  &
    +k(141)*n(idx_OH)  &
    +k(130)*n(idx_CH)  &
    -k(18)*n(idx_H)  &
    +2.d0*k(28)*n(idx_Hk)  &
    +k(146)*n(idx_H2O)  &
    +k(2)*n(idx_E)  &
    +k(10)*n(idx_HE)  &
    +k(47)*n(idx_C)  &
    +k(160)*n(idx_Ok)  &
    +k(142)*n(idx_OH)  &
    +k(134)*n(idx_CH2)  &
    +k(45)*n(idx_O)  &
    +2.d0*k(193)*n(idx_H2)  &
    +k(131)*n(idx_CH)

!d[HE_dot]/d[H+]
pd(6,20) =  &
    -k(10)*n(idx_HE)  &
    -k(9)*n(idx_HE)

!d[H2_dot]/d[H+]
pd(7,20) =  &
    -k(193)*n(idx_H2)  &
    -k(194)*n(idx_H2)  &
    +k(133)*n(idx_CH2)  &
    +k(132)*n(idx_CH2)  &
    -k(22)*n(idx_H2)  &
    -k(21)*n(idx_H2)

!d[C_dot]/d[H+]
pd(8,20) =  &
    +k(159)*n(idx_Ck)  &
    -k(47)*n(idx_C)

!d[O_dot]/d[H+]
pd(9,20) =  &
    -k(45)*n(idx_O)  &
    +k(160)*n(idx_Ok)

!d[OH_dot]/d[H+]
pd(10,20) =  &
    -k(141)*n(idx_OH)  &
    -k(142)*n(idx_OH)

!d[CH_dot]/d[H+]
pd(12,20) =  &
    -k(131)*n(idx_CH)  &
    -k(130)*n(idx_CH)

!d[CH2_dot]/d[H+]
pd(13,20) =  &
    -k(133)*n(idx_CH2)  &
    -k(132)*n(idx_CH2)  &
    -k(135)*n(idx_CH2)  &
    -k(134)*n(idx_CH2)

!d[H2O_dot]/d[H+]
pd(16,20) =  &
    -k(146)*n(idx_H2O)  &
    -k(145)*n(idx_H2O)

!d[O2_dot]/d[H+]
pd(17,20) =  &
    -k(153)*n(idx_O2)

!d[H+_dot]/d[H+]
pd(20,20) =  &
    -k(193)*n(idx_H2)  &
    -k(131)*n(idx_CH)  &
    -k(194)*n(idx_H2)  &
    -k(153)*n(idx_O2)  &
    -k(130)*n(idx_CH)  &
    -k(146)*n(idx_H2O)  &
    -k(134)*n(idx_CH2)  &
    -k(2)*n(idx_E)  &
    -k(19)*n(idx_H)  &
    -k(47)*n(idx_C)  &
    -k(145)*n(idx_H2O)  &
    -k(132)*n(idx_CH2)  &
    -k(18)*n(idx_H)  &
    -k(142)*n(idx_OH)  &
    -k(10)*n(idx_HE)  &
    -k(45)*n(idx_O)  &
    -k(160)*n(idx_Ok)  &
    -k(133)*n(idx_CH2)  &
    -k(22)*n(idx_H2)  &
    -k(29)*n(idx_Hk)  &
    -k(141)*n(idx_OH)  &
    -k(9)*n(idx_HE)  &
    -k(159)*n(idx_Ck)  &
    -k(28)*n(idx_Hk)  &
    -k(3)*n(idx_E)  &
    -k(135)*n(idx_CH2)  &
    +k(193)*n(idx_H2)  &
    -k(21)*n(idx_H2)

!d[HE+_dot]/d[H+]
pd(21,20) =  &
    +k(9)*n(idx_HE)  &
    +k(10)*n(idx_HE)

!d[H2+_dot]/d[H+]
pd(22,20) =  &
    +k(19)*n(idx_H)  &
    +k(29)*n(idx_Hk)  &
    +k(22)*n(idx_H2)  &
    +k(18)*n(idx_H)  &
    +k(21)*n(idx_H2)

!d[C+_dot]/d[H+]
pd(23,20) =  &
    +k(47)*n(idx_C)

!d[O+_dot]/d[H+]
pd(24,20) =  &
    +k(45)*n(idx_O)

!d[H3+_dot]/d[H+]
pd(27,20) =  &
    +k(194)*n(idx_H2)

!d[CH+_dot]/d[H+]
pd(28,20) =  &
    +k(130)*n(idx_CH)  &
    +k(133)*n(idx_CH2)  &
    +k(132)*n(idx_CH2)  &
    +k(131)*n(idx_CH)

!d[CH2+_dot]/d[H+]
pd(29,20) =  &
    +k(135)*n(idx_CH2)  &
    +k(134)*n(idx_CH2)

!d[OH+_dot]/d[H+]
pd(32,20) =  &
    +k(141)*n(idx_OH)  &
    +k(142)*n(idx_OH)

!d[H2O+_dot]/d[H+]
pd(33,20) =  &
    +k(145)*n(idx_H2O)  &
    +k(146)*n(idx_H2O)

!d[O2+_dot]/d[H+]
pd(35,20) =  &
    +k(153)*n(idx_O2)

!d[Tgas_dot]/d[H+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(20)*1d-3
if(dnn>0.d0) then
nn(20) = n(20) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,20) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[HE+]
pd(1,21) =  &
    -k(7)*n(idx_E)  &
    -k(6)*n(idx_E)  &
    +2.d0*k(7)*n(idx_E)  &
    +k(215)  &
    -k(5)*n(idx_E)

!d[H-_dot]/d[HE+]
pd(2,21) =  &
    -k(161)*n(idx_Hk)

!d[H_dot]/d[HE+]
pd(5,21) =  &
    +k(161)*n(idx_Hk)  &
    -k(8)*n(idx_H)  &
    +k(138)*n(idx_CH2)  &
    +k(150)*n(idx_H2O)  &
    +k(139)*n(idx_CH2)  &
    +k(149)*n(idx_H2O)  &
    +k(13)*n(idx_H2)  &
    +2.d0*k(14)*n(idx_H2)  &
    +k(143)*n(idx_OH)  &
    +k(144)*n(idx_OH)

!d[HE_dot]/d[HE+]
pd(6,21) =  &
    +k(157)*n(idx_CO)  &
    +k(138)*n(idx_CH2)  &
    +k(46)*n(idx_O)  &
    +k(50)*n(idx_C)  &
    +k(143)*n(idx_OH)  &
    +k(147)*n(idx_H2O)  &
    +k(156)*n(idx_CO)  &
    +k(8)*n(idx_H)  &
    +k(6)*n(idx_E)  &
    +k(49)*n(idx_C)  &
    +k(5)*n(idx_E)  &
    +k(140)*n(idx_C2)  &
    +k(152)*n(idx_H2O)  &
    +k(144)*n(idx_OH)  &
    +k(51)*n(idx_C)  &
    +k(151)*n(idx_H2O)  &
    +k(155)*n(idx_O2)  &
    +k(13)*n(idx_H2)  &
    +k(137)*n(idx_CH2)  &
    +k(148)*n(idx_H2O)  &
    +k(136)*n(idx_CH2)  &
    +k(161)*n(idx_Hk)  &
    +k(150)*n(idx_H2O)  &
    +k(139)*n(idx_CH2)  &
    +k(149)*n(idx_H2O)  &
    +k(12)*n(idx_H2)  &
    +k(154)*n(idx_O2)

!d[H2_dot]/d[HE+]
pd(7,21) =  &
    -k(13)*n(idx_H2)  &
    +k(136)*n(idx_CH2)  &
    -k(12)*n(idx_H2)  &
    +k(137)*n(idx_CH2)  &
    -k(14)*n(idx_H2)

!d[C_dot]/d[HE+]
pd(8,21) =  &
    -k(49)*n(idx_C)  &
    -k(51)*n(idx_C)  &
    +k(140)*n(idx_C2)  &
    +k(157)*n(idx_CO)  &
    -k(50)*n(idx_C)

!d[O_dot]/d[HE+]
pd(9,21) =  &
    +k(156)*n(idx_CO)  &
    -k(46)*n(idx_O)  &
    +k(155)*n(idx_O2)

!d[OH_dot]/d[HE+]
pd(10,21) =  &
    -k(143)*n(idx_OH)  &
    +k(147)*n(idx_H2O)  &
    +k(148)*n(idx_H2O)  &
    -k(144)*n(idx_OH)

!d[CO_dot]/d[HE+]
pd(11,21) =  &
    -k(157)*n(idx_CO)  &
    -k(156)*n(idx_CO)

!d[CH2_dot]/d[HE+]
pd(13,21) =  &
    -k(137)*n(idx_CH2)  &
    -k(136)*n(idx_CH2)  &
    -k(139)*n(idx_CH2)  &
    -k(138)*n(idx_CH2)

!d[C2_dot]/d[HE+]
pd(14,21) =  &
    -k(140)*n(idx_C2)

!d[H2O_dot]/d[HE+]
pd(16,21) =  &
    -k(150)*n(idx_H2O)  &
    -k(147)*n(idx_H2O)  &
    -k(149)*n(idx_H2O)  &
    -k(151)*n(idx_H2O)  &
    -k(148)*n(idx_H2O)  &
    -k(152)*n(idx_H2O)

!d[O2_dot]/d[HE+]
pd(17,21) =  &
    -k(155)*n(idx_O2)  &
    -k(154)*n(idx_O2)

!d[H+_dot]/d[HE+]
pd(20,21) =  &
    +k(148)*n(idx_H2O)  &
    +k(8)*n(idx_H)  &
    +k(147)*n(idx_H2O)  &
    +k(13)*n(idx_H2)

!d[HE+_dot]/d[HE+]
pd(21,21) =  &
    -k(8)*n(idx_H)  &
    -k(136)*n(idx_CH2)  &
    -k(5)*n(idx_E)  &
    -k(143)*n(idx_OH)  &
    -k(139)*n(idx_CH2)  &
    -k(13)*n(idx_H2)  &
    -k(161)*n(idx_Hk)  &
    -k(46)*n(idx_O)  &
    -k(6)*n(idx_E)  &
    -k(140)*n(idx_C2)  &
    -k(49)*n(idx_C)  &
    -k(155)*n(idx_O2)  &
    -k(148)*n(idx_H2O)  &
    +k(14)*n(idx_H2)  &
    -k(150)*n(idx_H2O)  &
    -k(50)*n(idx_C)  &
    -k(147)*n(idx_H2O)  &
    -k(215)  &
    -k(138)*n(idx_CH2)  &
    -k(51)*n(idx_C)  &
    -k(157)*n(idx_CO)  &
    -k(154)*n(idx_O2)  &
    -k(152)*n(idx_H2O)  &
    -k(7)*n(idx_E)  &
    -k(149)*n(idx_H2O)  &
    -k(151)*n(idx_H2O)  &
    -k(144)*n(idx_OH)  &
    -k(12)*n(idx_H2)  &
    -k(137)*n(idx_CH2)  &
    -k(14)*n(idx_H2)  &
    -k(156)*n(idx_CO)

!d[H2+_dot]/d[HE+]
pd(22,21) =  &
    +k(12)*n(idx_H2)

!d[C+_dot]/d[HE+]
pd(23,21) =  &
    +k(51)*n(idx_C)  &
    +k(156)*n(idx_CO)  &
    +k(49)*n(idx_C)  &
    +k(140)*n(idx_C2)  &
    +k(137)*n(idx_CH2)  &
    +k(50)*n(idx_C)  &
    +k(136)*n(idx_CH2)

!d[O+_dot]/d[HE+]
pd(24,21) =  &
    +k(157)*n(idx_CO)  &
    +k(155)*n(idx_O2)  &
    +k(46)*n(idx_O)  &
    +k(143)*n(idx_OH)  &
    +k(144)*n(idx_OH)

!d[CH+_dot]/d[HE+]
pd(28,21) =  &
    +k(139)*n(idx_CH2)  &
    +k(138)*n(idx_CH2)

!d[OH+_dot]/d[HE+]
pd(32,21) =  &
    +k(149)*n(idx_H2O)  &
    +k(150)*n(idx_H2O)

!d[H2O+_dot]/d[HE+]
pd(33,21) =  &
    +k(152)*n(idx_H2O)  &
    +k(151)*n(idx_H2O)

!d[O2+_dot]/d[HE+]
pd(35,21) =  &
    +k(154)*n(idx_O2)

!d[HE++_dot]/d[HE+]
pd(36,21) =  &
    +k(7)*n(idx_E)  &
    +k(215)

!d[Tgas_dot]/d[HE+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(21)*1d-3
if(dnn>0.d0) then
nn(21) = n(21) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,21) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H2+]
pd(1,22) =  &
    -k(31)*n(idx_E)  &
    -k(30)*n(idx_E)

!d[H-_dot]/d[H2+]
pd(2,22) =  &
    -k(32)*n(idx_Hk)

!d[H_dot]/d[H2+]
pd(5,22) =  &
    +2.d0*k(30)*n(idx_E)  &
    -k(20)*n(idx_H)  &
    +k(87)*n(idx_C)  &
    +k(32)*n(idx_Hk)  &
    +k(85)*n(idx_H2)  &
    +2.d0*k(31)*n(idx_E)  &
    +k(102)*n(idx_O)  &
    +k(232)

!d[H2_dot]/d[H2+]
pd(7,22) =  &
    +k(20)*n(idx_H)  &
    +k(32)*n(idx_Hk)  &
    -k(85)*n(idx_H2)

!d[C_dot]/d[H2+]
pd(8,22) =  &
    -k(87)*n(idx_C)

!d[O_dot]/d[H2+]
pd(9,22) =  &
    -k(102)*n(idx_O)

!d[H+_dot]/d[H2+]
pd(20,22) =  &
    +k(20)*n(idx_H)  &
    +k(232)

!d[H2+_dot]/d[H2+]
pd(22,22) =  &
    -k(85)*n(idx_H2)  &
    -k(20)*n(idx_H)  &
    -k(87)*n(idx_C)  &
    -k(31)*n(idx_E)  &
    -k(30)*n(idx_E)  &
    -k(232)  &
    -k(102)*n(idx_O)  &
    -k(32)*n(idx_Hk)

!d[H3+_dot]/d[H2+]
pd(27,22) =  &
    +k(85)*n(idx_H2)

!d[CH+_dot]/d[H2+]
pd(28,22) =  &
    +k(87)*n(idx_C)

!d[OH+_dot]/d[H2+]
pd(32,22) =  &
    +k(102)*n(idx_O)

!d[Tgas_dot]/d[H2+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(22)*1d-3
if(dnn>0.d0) then
nn(22) = n(22) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,22) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[C+]
pd(1,23) =  &
    -k(37)*n(idx_E)  &
    -k(39)*n(idx_E)  &
    -k(38)*n(idx_E)

!d[H_dot]/d[C+]
pd(5,23) =  &
    +k(90)*n(idx_H2)  &
    +k(107)*n(idx_OH)  &
    +k(108)*n(idx_OH)  &
    -k(48)*n(idx_H)  &
    -k(200)*n(idx_H)  &
    +k(113)*n(idx_H2O)  &
    +k(115)*n(idx_H2O)  &
    +k(114)*n(idx_H2O)

!d[H2_dot]/d[C+]
pd(7,23) =  &
    -k(201)*n(idx_H2)  &
    -k(90)*n(idx_H2)

!d[C_dot]/d[C+]
pd(8,23) =  &
    +k(39)*n(idx_E)  &
    +k(38)*n(idx_E)  &
    +k(37)*n(idx_E)  &
    +k(48)*n(idx_H)  &
    +k(116)*n(idx_H2O)

!d[O_dot]/d[C+]
pd(9,23) =  &
    -k(275)*n(idx_O)  &
    -k(276)*n(idx_O)  &
    -k(202)*n(idx_O)  &
    +k(118)*n(idx_O2)  &
    -k(203)*n(idx_O)

!d[OH_dot]/d[C+]
pd(10,23) =  &
    -k(107)*n(idx_OH)  &
    -k(108)*n(idx_OH)

!d[CO_dot]/d[C+]
pd(11,23) =  &
    +k(119)*n(idx_O2)

!d[H2O_dot]/d[C+]
pd(16,23) =  &
    -k(114)*n(idx_H2O)  &
    -k(115)*n(idx_H2O)  &
    -k(113)*n(idx_H2O)  &
    -k(116)*n(idx_H2O)

!d[O2_dot]/d[C+]
pd(17,23) =  &
    -k(119)*n(idx_O2)  &
    -k(118)*n(idx_O2)

!d[H+_dot]/d[C+]
pd(20,23) =  &
    +k(48)*n(idx_H)

!d[C+_dot]/d[C+]
pd(23,23) =  &
    -k(39)*n(idx_E)  &
    -k(38)*n(idx_E)  &
    -k(275)*n(idx_O)  &
    -k(203)*n(idx_O)  &
    -k(114)*n(idx_H2O)  &
    -k(37)*n(idx_E)  &
    -k(107)*n(idx_OH)  &
    -k(276)*n(idx_O)  &
    -k(108)*n(idx_OH)  &
    -k(116)*n(idx_H2O)  &
    -k(202)*n(idx_O)  &
    -k(48)*n(idx_H)  &
    -k(90)*n(idx_H2)  &
    -k(201)*n(idx_H2)  &
    -k(113)*n(idx_H2O)  &
    -k(200)*n(idx_H)  &
    -k(115)*n(idx_H2O)  &
    -k(119)*n(idx_O2)  &
    -k(118)*n(idx_O2)

!d[O+_dot]/d[C+]
pd(24,23) =  &
    +k(119)*n(idx_O2)

!d[HOC+_dot]/d[C+]
pd(25,23) =  &
    +k(113)*n(idx_H2O)

!d[HCO+_dot]/d[C+]
pd(26,23) =  &
    +k(115)*n(idx_H2O)  &
    +k(114)*n(idx_H2O)

!d[CH+_dot]/d[C+]
pd(28,23) =  &
    +k(90)*n(idx_H2)  &
    +k(200)*n(idx_H)

!d[CH2+_dot]/d[C+]
pd(29,23) =  &
    +k(201)*n(idx_H2)

!d[CO+_dot]/d[C+]
pd(30,23) =  &
    +k(118)*n(idx_O2)  &
    +k(107)*n(idx_OH)  &
    +k(202)*n(idx_O)  &
    +k(276)*n(idx_O)  &
    +k(108)*n(idx_OH)  &
    +k(203)*n(idx_O)  &
    +k(275)*n(idx_O)

!d[H2O+_dot]/d[C+]
pd(33,23) =  &
    +k(116)*n(idx_H2O)

!d[Tgas_dot]/d[C+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(23)*1d-3
if(dnn>0.d0) then
nn(23) = n(23) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,23) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[O+]
pd(1,24) =  &
    -k(40)*n(idx_E)  &
    -k(41)*n(idx_E)

!d[H_dot]/d[O+]
pd(5,24) =  &
    -k(44)*n(idx_H)  &
    +k(101)*n(idx_H2)

!d[H2_dot]/d[O+]
pd(7,24) =  &
    -k(101)*n(idx_H2)

!d[C_dot]/d[O+]
pd(8,24) =  &
    +k(100)*n(idx_C2)  &
    -k(277)*n(idx_C)  &
    -k(278)*n(idx_C)

!d[O_dot]/d[O+]
pd(9,24) =  &
    +k(40)*n(idx_E)  &
    +k(41)*n(idx_E)  &
    +k(44)*n(idx_H)

!d[C2_dot]/d[O+]
pd(14,24) =  &
    -k(100)*n(idx_C2)

!d[H+_dot]/d[O+]
pd(20,24) =  &
    +k(44)*n(idx_H)

!d[O+_dot]/d[O+]
pd(24,24) =  &
    -k(100)*n(idx_C2)  &
    -k(101)*n(idx_H2)  &
    -k(44)*n(idx_H)  &
    -k(278)*n(idx_C)  &
    -k(277)*n(idx_C)  &
    -k(40)*n(idx_E)  &
    -k(41)*n(idx_E)

!d[CO+_dot]/d[O+]
pd(30,24) =  &
    +k(100)*n(idx_C2)  &
    +k(278)*n(idx_C)  &
    +k(277)*n(idx_C)

!d[OH+_dot]/d[O+]
pd(32,24) =  &
    +k(101)*n(idx_H2)

!d[Tgas_dot]/d[O+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(24)*1d-3
if(dnn>0.d0) then
nn(24) = n(24) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,24) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[HOC+]
pd(1,25) =  &
    -k(183)*n(idx_E)

!d[H_dot]/d[HOC+]
pd(5,25) =  &
    +k(183)*n(idx_E)

!d[H2_dot]/d[HOC+]
pd(7,25) =  &
    +k(53)*n(idx_H2)  &
    -k(53)*n(idx_H2)

!d[CO_dot]/d[HOC+]
pd(11,25) =  &
    +k(54)*n(idx_CO)  &
    +k(183)*n(idx_E)  &
    -k(55)*n(idx_CO)  &
    +k(55)*n(idx_CO)  &
    -k(54)*n(idx_CO)

!d[HOC+_dot]/d[HOC+]
pd(25,25) =  &
    -k(183)*n(idx_E)  &
    -k(55)*n(idx_CO)  &
    -k(53)*n(idx_H2)  &
    -k(54)*n(idx_CO)

!d[HCO+_dot]/d[HOC+]
pd(26,25) =  &
    +k(53)*n(idx_H2)  &
    +k(54)*n(idx_CO)  &
    +k(55)*n(idx_CO)

!d[Tgas_dot]/d[HOC+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(25)*1d-3
if(dnn>0.d0) then
nn(25) = n(25) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,25) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[HCO+]
pd(1,26) =  &
    -k(182)*n(idx_E)  &
    -k(181)*n(idx_E)

!d[H_dot]/d[HCO+]
pd(5,26) =  &
    +k(181)*n(idx_E)

!d[C_dot]/d[HCO+]
pd(8,26) =  &
    +k(182)*n(idx_E)  &
    -k(127)*n(idx_C)

!d[OH_dot]/d[HCO+]
pd(10,26) =  &
    +k(182)*n(idx_E)

!d[CO_dot]/d[HCO+]
pd(11,26) =  &
    +k(127)*n(idx_C)  &
    +k(129)*n(idx_H2O)  &
    +k(128)*n(idx_H2O)  &
    +k(181)*n(idx_E)

!d[H2O_dot]/d[HCO+]
pd(16,26) =  &
    -k(129)*n(idx_H2O)  &
    -k(128)*n(idx_H2O)

!d[HCO+_dot]/d[HCO+]
pd(26,26) =  &
    -k(129)*n(idx_H2O)  &
    -k(182)*n(idx_E)  &
    -k(127)*n(idx_C)  &
    -k(128)*n(idx_H2O)  &
    -k(181)*n(idx_E)

!d[CH+_dot]/d[HCO+]
pd(28,26) =  &
    +k(127)*n(idx_C)

!d[H3O+_dot]/d[HCO+]
pd(34,26) =  &
    +k(129)*n(idx_H2O)  &
    +k(128)*n(idx_H2O)

!d[Tgas_dot]/d[HCO+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(26)*1d-3
if(dnn>0.d0) then
nn(26) = n(26) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,26) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H3+]
pd(1,27) =  &
    -k(162)*n(idx_E)  &
    -k(163)*n(idx_E)

!d[H_dot]/d[H3+]
pd(5,27) =  &
    +3.d0*k(163)*n(idx_E)  &
    +k(104)*n(idx_O)  &
    -k(86)*n(idx_H)  &
    +k(234)  &
    +k(162)*n(idx_E)  &
    +k(89)*n(idx_C)

!d[H2_dot]/d[H3+]
pd(7,27) =  &
    +k(233)  &
    +k(86)*n(idx_H)  &
    +k(125)*n(idx_CO)  &
    +k(124)*n(idx_CO)  &
    +k(111)*n(idx_H2O)  &
    +k(123)*n(idx_CO)  &
    +k(126)*n(idx_CO)  &
    +k(105)*n(idx_OH)  &
    +k(106)*n(idx_OH)  &
    +k(88)*n(idx_C)  &
    +k(103)*n(idx_O)  &
    +k(162)*n(idx_E)  &
    +k(112)*n(idx_H2O)

!d[C_dot]/d[H3+]
pd(8,27) =  &
    -k(88)*n(idx_C)  &
    -k(89)*n(idx_C)

!d[O_dot]/d[H3+]
pd(9,27) =  &
    -k(104)*n(idx_O)  &
    -k(103)*n(idx_O)

!d[OH_dot]/d[H3+]
pd(10,27) =  &
    -k(105)*n(idx_OH)  &
    -k(106)*n(idx_OH)

!d[CO_dot]/d[H3+]
pd(11,27) =  &
    -k(126)*n(idx_CO)  &
    -k(125)*n(idx_CO)  &
    -k(124)*n(idx_CO)  &
    -k(123)*n(idx_CO)

!d[H2O_dot]/d[H3+]
pd(16,27) =  &
    -k(111)*n(idx_H2O)  &
    -k(112)*n(idx_H2O)

!d[H+_dot]/d[H3+]
pd(20,27) =  &
    +k(233)

!d[H2+_dot]/d[H3+]
pd(22,27) =  &
    +k(86)*n(idx_H)  &
    +k(234)

!d[HOC+_dot]/d[H3+]
pd(25,27) =  &
    +k(125)*n(idx_CO)  &
    +k(126)*n(idx_CO)

!d[HCO+_dot]/d[H3+]
pd(26,27) =  &
    +k(123)*n(idx_CO)  &
    +k(124)*n(idx_CO)

!d[H3+_dot]/d[H3+]
pd(27,27) =  &
    -k(124)*n(idx_CO)  &
    -k(123)*n(idx_CO)  &
    -k(89)*n(idx_C)  &
    -k(88)*n(idx_C)  &
    -k(125)*n(idx_CO)  &
    -k(112)*n(idx_H2O)  &
    -k(233)  &
    -k(111)*n(idx_H2O)  &
    -k(86)*n(idx_H)  &
    -k(162)*n(idx_E)  &
    -k(163)*n(idx_E)  &
    -k(106)*n(idx_OH)  &
    -k(126)*n(idx_CO)  &
    -k(104)*n(idx_O)  &
    -k(234)  &
    -k(105)*n(idx_OH)  &
    -k(103)*n(idx_O)

!d[CH+_dot]/d[H3+]
pd(28,27) =  &
    +k(88)*n(idx_C)

!d[CH2+_dot]/d[H3+]
pd(29,27) =  &
    +k(89)*n(idx_C)

!d[OH+_dot]/d[H3+]
pd(32,27) =  &
    +k(103)*n(idx_O)

!d[H2O+_dot]/d[H3+]
pd(33,27) =  &
    +k(104)*n(idx_O)  &
    +k(105)*n(idx_OH)  &
    +k(106)*n(idx_OH)

!d[H3O+_dot]/d[H3+]
pd(34,27) =  &
    +k(112)*n(idx_H2O)  &
    +k(111)*n(idx_H2O)

!d[Tgas_dot]/d[H3+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(27)*1d-3
if(dnn>0.d0) then
nn(27) = n(27) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,27) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[CH+]
pd(1,28) =  &
    -k(164)*n(idx_E)

!d[H_dot]/d[CH+]
pd(5,28) =  &
    +k(164)*n(idx_E)  &
    +k(92)*n(idx_H2)  &
    +k(93)*n(idx_O)  &
    -k(91)*n(idx_H)

!d[H2_dot]/d[CH+]
pd(7,28) =  &
    +k(91)*n(idx_H)  &
    -k(92)*n(idx_H2)

!d[C_dot]/d[CH+]
pd(8,28) =  &
    +k(164)*n(idx_E)  &
    +k(236)

!d[O_dot]/d[CH+]
pd(9,28) =  &
    -k(93)*n(idx_O)

!d[H+_dot]/d[CH+]
pd(20,28) =  &
    +k(236)

!d[C+_dot]/d[CH+]
pd(23,28) =  &
    +k(91)*n(idx_H)

!d[CH+_dot]/d[CH+]
pd(28,28) =  &
    -k(236)  &
    -k(164)*n(idx_E)  &
    -k(93)*n(idx_O)  &
    -k(92)*n(idx_H2)  &
    -k(91)*n(idx_H)

!d[CH2+_dot]/d[CH+]
pd(29,28) =  &
    +k(92)*n(idx_H2)

!d[CO+_dot]/d[CH+]
pd(30,28) =  &
    +k(93)*n(idx_O)

!d[Tgas_dot]/d[CH+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(28)*1d-3
if(dnn>0.d0) then
nn(28) = n(28) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,28) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[CH2+]
pd(1,29) =  &
    -k(165)*n(idx_E)  &
    -k(166)*n(idx_E)  &
    -k(167)*n(idx_E)

!d[H_dot]/d[CH2+]
pd(5,29) =  &
    -k(94)*n(idx_H)  &
    +k(96)*n(idx_O)  &
    +k(95)*n(idx_H2)  &
    +k(165)*n(idx_E)  &
    +k(239)  &
    +2.d0*k(167)*n(idx_E)

!d[H2_dot]/d[CH2+]
pd(7,29) =  &
    +k(94)*n(idx_H)  &
    -k(95)*n(idx_H2)  &
    +k(166)*n(idx_E)

!d[C_dot]/d[CH2+]
pd(8,29) =  &
    +k(166)*n(idx_E)  &
    +k(167)*n(idx_E)

!d[O_dot]/d[CH2+]
pd(9,29) =  &
    -k(96)*n(idx_O)

!d[OH_dot]/d[CH2+]
pd(10,29) =  &
    +k(120)*n(idx_O2)

!d[CH_dot]/d[CH2+]
pd(12,29) =  &
    +k(165)*n(idx_E)

!d[O2_dot]/d[CH2+]
pd(17,29) =  &
    -k(120)*n(idx_O2)

!d[HCO+_dot]/d[CH2+]
pd(26,29) =  &
    +k(120)*n(idx_O2)  &
    +k(96)*n(idx_O)

!d[CH+_dot]/d[CH2+]
pd(28,29) =  &
    +k(94)*n(idx_H)  &
    +k(239)

!d[CH2+_dot]/d[CH2+]
pd(29,29) =  &
    -k(94)*n(idx_H)  &
    -k(95)*n(idx_H2)  &
    -k(166)*n(idx_E)  &
    -k(96)*n(idx_O)  &
    -k(167)*n(idx_E)  &
    -k(239)  &
    -k(165)*n(idx_E)  &
    -k(120)*n(idx_O2)

!d[CH3+_dot]/d[CH2+]
pd(31,29) =  &
    +k(95)*n(idx_H2)

!d[Tgas_dot]/d[CH2+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(29)*1d-3
if(dnn>0.d0) then
nn(29) = n(29) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,29) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[CO+]
pd(1,30) =  &
    -k(180)*n(idx_E)

!d[H_dot]/d[CO+]
pd(5,30) =  &
    -k(158)*n(idx_H)

!d[C_dot]/d[CO+]
pd(8,30) =  &
    +k(180)*n(idx_E)

!d[O_dot]/d[CO+]
pd(9,30) =  &
    +k(180)*n(idx_E)

!d[CO_dot]/d[CO+]
pd(11,30) =  &
    +k(158)*n(idx_H)

!d[H+_dot]/d[CO+]
pd(20,30) =  &
    +k(158)*n(idx_H)

!d[CO+_dot]/d[CO+]
pd(30,30) =  &
    -k(158)*n(idx_H)  &
    -k(180)*n(idx_E)

!d[Tgas_dot]/d[CO+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(30)*1d-3
if(dnn>0.d0) then
nn(30) = n(30) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,30) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[CH3+]
pd(1,31) =  &
    -k(170)*n(idx_E)  &
    -k(168)*n(idx_E)  &
    -k(169)*n(idx_E)

!d[H_dot]/d[CH3+]
pd(5,31) =  &
    -k(97)*n(idx_H)  &
    +k(240)  &
    +k(168)*n(idx_E)  &
    +2.d0*k(170)*n(idx_E)

!d[H2_dot]/d[CH3+]
pd(7,31) =  &
    +k(99)*n(idx_O)  &
    +k(241)  &
    +k(169)*n(idx_E)  &
    +k(98)*n(idx_O)  &
    +k(97)*n(idx_H)

!d[O_dot]/d[CH3+]
pd(9,31) =  &
    -k(99)*n(idx_O)  &
    -k(98)*n(idx_O)

!d[CH_dot]/d[CH3+]
pd(12,31) =  &
    +k(169)*n(idx_E)  &
    +k(170)*n(idx_E)

!d[CH2_dot]/d[CH3+]
pd(13,31) =  &
    +k(168)*n(idx_E)

!d[HOC+_dot]/d[CH3+]
pd(25,31) =  &
    +k(98)*n(idx_O)

!d[HCO+_dot]/d[CH3+]
pd(26,31) =  &
    +k(99)*n(idx_O)

!d[CH+_dot]/d[CH3+]
pd(28,31) =  &
    +k(241)

!d[CH2+_dot]/d[CH3+]
pd(29,31) =  &
    +k(240)  &
    +k(97)*n(idx_H)

!d[CH3+_dot]/d[CH3+]
pd(31,31) =  &
    -k(169)*n(idx_E)  &
    -k(97)*n(idx_H)  &
    -k(99)*n(idx_O)  &
    -k(98)*n(idx_O)  &
    -k(170)*n(idx_E)  &
    -k(240)  &
    -k(241)  &
    -k(168)*n(idx_E)

!d[Tgas_dot]/d[CH3+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(31)*1d-3
if(dnn>0.d0) then
nn(31) = n(31) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,31) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[OH+]
pd(1,32) =  &
    -k(171)*n(idx_E)

!d[H_dot]/d[OH+]
pd(5,32) =  &
    +k(109)*n(idx_H2)  &
    +k(171)*n(idx_E)

!d[H2_dot]/d[OH+]
pd(7,32) =  &
    -k(109)*n(idx_H2)

!d[O_dot]/d[OH+]
pd(9,32) =  &
    +k(243)  &
    +k(171)*n(idx_E)

!d[H+_dot]/d[OH+]
pd(20,32) =  &
    +k(243)

!d[OH+_dot]/d[OH+]
pd(32,32) =  &
    -k(171)*n(idx_E)  &
    -k(109)*n(idx_H2)  &
    -k(243)

!d[H2O+_dot]/d[OH+]
pd(33,32) =  &
    +k(109)*n(idx_H2)

!d[Tgas_dot]/d[OH+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(32)*1d-3
if(dnn>0.d0) then
nn(32) = n(32) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,32) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H2O+]
pd(1,33) =  &
    -k(174)*n(idx_E)  &
    -k(173)*n(idx_E)  &
    -k(172)*n(idx_E)

!d[H_dot]/d[H2O+]
pd(5,33) =  &
    +k(173)*n(idx_E)  &
    +k(247)  &
    +2.d0*k(174)*n(idx_E)  &
    +k(110)*n(idx_H2)

!d[H2_dot]/d[H2O+]
pd(7,33) =  &
    +k(246)  &
    +k(172)*n(idx_E)  &
    -k(110)*n(idx_H2)

!d[O_dot]/d[H2O+]
pd(9,33) =  &
    +k(172)*n(idx_E)  &
    +k(174)*n(idx_E)  &
    +k(244)

!d[OH_dot]/d[H2O+]
pd(10,33) =  &
    +k(173)*n(idx_E)  &
    +k(245)

!d[H+_dot]/d[H2O+]
pd(20,33) =  &
    +k(245)

!d[H2+_dot]/d[H2O+]
pd(22,33) =  &
    +k(244)

!d[O+_dot]/d[H2O+]
pd(24,33) =  &
    +k(246)

!d[OH+_dot]/d[H2O+]
pd(32,33) =  &
    +k(247)

!d[H2O+_dot]/d[H2O+]
pd(33,33) =  &
    -k(110)*n(idx_H2)  &
    -k(173)*n(idx_E)  &
    -k(245)  &
    -k(244)  &
    -k(174)*n(idx_E)  &
    -k(172)*n(idx_E)  &
    -k(247)  &
    -k(246)

!d[H3O+_dot]/d[H2O+]
pd(34,33) =  &
    +k(110)*n(idx_H2)

!d[Tgas_dot]/d[H2O+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(33)*1d-3
if(dnn>0.d0) then
nn(33) = n(33) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,33) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[H3O+]
pd(1,34) =  &
    -k(178)*n(idx_E)  &
    -k(177)*n(idx_E)  &
    -k(176)*n(idx_E)  &
    -k(175)*n(idx_E)

!d[H_dot]/d[H3O+]
pd(5,34) =  &
    +2.d0*k(175)*n(idx_E)  &
    +k(250)  &
    +k(177)*n(idx_E)  &
    +k(176)*n(idx_E)

!d[H2_dot]/d[H3O+]
pd(7,34) =  &
    +k(251)  &
    +k(117)*n(idx_C)  &
    +k(176)*n(idx_E)  &
    +k(178)*n(idx_E)

!d[C_dot]/d[H3O+]
pd(8,34) =  &
    -k(117)*n(idx_C)

!d[O_dot]/d[H3O+]
pd(9,34) =  &
    +k(176)*n(idx_E)

!d[OH_dot]/d[H3O+]
pd(10,34) =  &
    +k(249)  &
    +k(175)*n(idx_E)  &
    +k(178)*n(idx_E)

!d[H2O_dot]/d[H3O+]
pd(16,34) =  &
    +k(248)  &
    +k(177)*n(idx_E)

!d[H+_dot]/d[H3O+]
pd(20,34) =  &
    +k(248)

!d[H2+_dot]/d[H3O+]
pd(22,34) =  &
    +k(249)

!d[HCO+_dot]/d[H3O+]
pd(26,34) =  &
    +k(117)*n(idx_C)

!d[OH+_dot]/d[H3O+]
pd(32,34) =  &
    +k(251)

!d[H2O+_dot]/d[H3O+]
pd(33,34) =  &
    +k(250)

!d[H3O+_dot]/d[H3O+]
pd(34,34) =  &
    -k(248)  &
    -k(178)*n(idx_E)  &
    -k(117)*n(idx_C)  &
    -k(175)*n(idx_E)  &
    -k(250)  &
    -k(251)  &
    -k(176)*n(idx_E)  &
    -k(177)*n(idx_E)  &
    -k(249)

!d[Tgas_dot]/d[H3O+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(34)*1d-3
if(dnn>0.d0) then
nn(34) = n(34) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,34) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[O2+]
pd(1,35) =  &
    -k(179)*n(idx_E)

!d[C_dot]/d[O2+]
pd(8,35) =  &
    -k(121)*n(idx_C)  &
    -k(122)*n(idx_C)

!d[O_dot]/d[O2+]
pd(9,35) =  &
    +k(121)*n(idx_C)  &
    +2.d0*k(179)*n(idx_E)

!d[O2_dot]/d[O2+]
pd(17,35) =  &
    +k(122)*n(idx_C)

!d[C+_dot]/d[O2+]
pd(23,35) =  &
    +k(122)*n(idx_C)

!d[CO+_dot]/d[O2+]
pd(30,35) =  &
    +k(121)*n(idx_C)

!d[O2+_dot]/d[O2+]
pd(35,35) =  &
    -k(179)*n(idx_E)  &
    -k(121)*n(idx_C)  &
    -k(122)*n(idx_C)

!d[Tgas_dot]/d[O2+]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(35)*1d-3
if(dnn>0.d0) then
nn(35) = n(35) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,35) = (dn1-dn0)/dnn
end if

!d[E_dot]/d[HE++]
pd(1,36) =  &
    -k(15)*n(idx_E)

!d[HE+_dot]/d[HE++]
pd(21,36) =  &
    +k(15)*n(idx_E)

!d[HE++_dot]/d[HE++]
pd(36,36) =  &
    -k(15)*n(idx_E)

!d[Tgas_dot]/d[HE++]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(36)*1d-3
if(dnn>0.d0) then
nn(36) = n(36) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,36) = (dn1-dn0)/dnn
end if

!d[Tgas_dot]/d[CR]
pd(39,37) = 0.d0

!d[Tgas_dot]/d[g]
pd(39,38) = 0.d0

!d[E_dot]/d[Tgas]
pd(1,39) = 0.d0

!d[H-_dot]/d[Tgas]
pd(2,39) = 0.d0

!d[C-_dot]/d[Tgas]
pd(3,39) = 0.d0

!d[O-_dot]/d[Tgas]
pd(4,39) = 0.d0

!d[H_dot]/d[Tgas]
pd(5,39) = 0.d0

!d[HE_dot]/d[Tgas]
pd(6,39) = 0.d0

!d[H2_dot]/d[Tgas]
pd(7,39) = 0.d0

!d[C_dot]/d[Tgas]
pd(8,39) = 0.d0

!d[O_dot]/d[Tgas]
pd(9,39) = 0.d0

!d[OH_dot]/d[Tgas]
pd(10,39) = 0.d0

!d[CO_dot]/d[Tgas]
pd(11,39) = 0.d0

!d[CH_dot]/d[Tgas]
pd(12,39) = 0.d0

!d[CH2_dot]/d[Tgas]
pd(13,39) = 0.d0

!d[C2_dot]/d[Tgas]
pd(14,39) = 0.d0

!d[HCO_dot]/d[Tgas]
pd(15,39) = 0.d0

!d[H2O_dot]/d[Tgas]
pd(16,39) = 0.d0

!d[O2_dot]/d[Tgas]
pd(17,39) = 0.d0

!d[CO_total_dot]/d[Tgas]
pd(18,39) = 0.d0

!d[H2O_total_dot]/d[Tgas]
pd(19,39) = 0.d0

!d[H+_dot]/d[Tgas]
pd(20,39) = 0.d0

!d[HE+_dot]/d[Tgas]
pd(21,39) = 0.d0

!d[H2+_dot]/d[Tgas]
pd(22,39) = 0.d0

!d[C+_dot]/d[Tgas]
pd(23,39) = 0.d0

!d[O+_dot]/d[Tgas]
pd(24,39) = 0.d0

!d[HOC+_dot]/d[Tgas]
pd(25,39) = 0.d0

!d[HCO+_dot]/d[Tgas]
pd(26,39) = 0.d0

!d[H3+_dot]/d[Tgas]
pd(27,39) = 0.d0

!d[CH+_dot]/d[Tgas]
pd(28,39) = 0.d0

!d[CH2+_dot]/d[Tgas]
pd(29,39) = 0.d0

!d[CO+_dot]/d[Tgas]
pd(30,39) = 0.d0

!d[CH3+_dot]/d[Tgas]
pd(31,39) = 0.d0

!d[OH+_dot]/d[Tgas]
pd(32,39) = 0.d0

!d[H2O+_dot]/d[Tgas]
pd(33,39) = 0.d0

!d[H3O+_dot]/d[Tgas]
pd(34,39) = 0.d0

!d[O2+_dot]/d[Tgas]
pd(35,39) = 0.d0

!d[HE++_dot]/d[Tgas]
pd(36,39) = 0.d0

!d[CR_dot]/d[Tgas]
pd(37,39) = 0.d0

!d[g_dot]/d[Tgas]
pd(38,39) = 0.d0

!d[Tgas_dot]/d[Tgas]
dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
nn(:) = n(:)
dnn = n(39)*1d-3
if(dnn>0.d0) then
nn(39) = n(39) + dnn
dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
pd(idx_Tgas,39) = (dn1-dn0)/dnn
end if

!d[dummy_dot]/d[Tgas]
pd(40,39) = 0.d0

!d[Tgas_dot]/d[dummy]
pd(39,40) = 0.d0

end subroutine jex

end module krome_ode

!############### MODULE ##############
module krome_user
implicit none

! *************************************************************
!  This file has been generated with:
!  KROME 14.08.dev on 2018-10-24 12:49:10
!  Changeset b21f657
!  see http://kromepackage.org
!
!  Written and developed by Tommaso Grassi and Stefano Bovino
!
!  Contributors:
!  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
!  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
!  E.Tognelli
!  KROME is provided "as it is", without any warranty.
! *************************************************************

integer,parameter::KROME_idx_E = 1	!E
integer,parameter::KROME_idx_Hk = 2	!H-
integer,parameter::KROME_idx_Ck = 3	!C-
integer,parameter::KROME_idx_Ok = 4	!O-
integer,parameter::KROME_idx_H = 5	!H
integer,parameter::KROME_idx_HE = 6	!HE
integer,parameter::KROME_idx_H2 = 7	!H2
integer,parameter::KROME_idx_C = 8	!C
integer,parameter::KROME_idx_O = 9	!O
integer,parameter::KROME_idx_OH = 10	!OH
integer,parameter::KROME_idx_CO = 11	!CO
integer,parameter::KROME_idx_CH = 12	!CH
integer,parameter::KROME_idx_CH2 = 13	!CH2
integer,parameter::KROME_idx_C2 = 14	!C2
integer,parameter::KROME_idx_HCO = 15	!HCO
integer,parameter::KROME_idx_H2O = 16	!H2O
integer,parameter::KROME_idx_O2 = 17	!O2
integer,parameter::KROME_idx_CO_total = 18	!CO_total
integer,parameter::KROME_idx_H2O_total = 19	!H2O_total
integer,parameter::KROME_idx_Hj = 20	!H+
integer,parameter::KROME_idx_HEj = 21	!HE+
integer,parameter::KROME_idx_H2j = 22	!H2+
integer,parameter::KROME_idx_Cj = 23	!C+
integer,parameter::KROME_idx_Oj = 24	!O+
integer,parameter::KROME_idx_HOCj = 25	!HOC+
integer,parameter::KROME_idx_HCOj = 26	!HCO+
integer,parameter::KROME_idx_H3j = 27	!H3+
integer,parameter::KROME_idx_CHj = 28	!CH+
integer,parameter::KROME_idx_CH2j = 29	!CH2+
integer,parameter::KROME_idx_COj = 30	!CO+
integer,parameter::KROME_idx_CH3j = 31	!CH3+
integer,parameter::KROME_idx_OHj = 32	!OH+
integer,parameter::KROME_idx_H2Oj = 33	!H2O+
integer,parameter::KROME_idx_H3Oj = 34	!H3O+
integer,parameter::KROME_idx_O2j = 35	!O2+
integer,parameter::KROME_idx_HEjj = 36	!HE++
integer,parameter::KROME_idx_CR = 37	!CR
integer,parameter::KROME_idx_g = 38	!g
integer,parameter::KROME_idx_Tgas = 39	!Tgas
integer,parameter::KROME_idx_dummy = 40	!dummy

integer,parameter::krome_idx_cool_h2 = 1
integer,parameter::krome_idx_cool_h2gp = 2
integer,parameter::krome_idx_cool_atomic = 3
integer,parameter::krome_idx_cool_cen = 3
integer,parameter::krome_idx_cool_hd = 4
integer,parameter::krome_idx_cool_metal = 5
integer,parameter::krome_idx_cool_z = 5
integer,parameter::krome_idx_cool_dh = 6
integer,parameter::krome_idx_cool_enthalpic = 6
integer,parameter::krome_idx_cool_dust = 7
integer,parameter::krome_idx_cool_compton = 8
integer,parameter::krome_idx_cool_cie = 9
integer,parameter::krome_idx_cool_cont = 10
integer,parameter::krome_idx_cool_continuum = 10
integer,parameter::krome_idx_cool_expansion = 11
integer,parameter::krome_idx_cool_exp = 11
integer,parameter::krome_idx_cool_ff = 12
integer,parameter::krome_idx_cool_bss = 12
integer,parameter::krome_idx_cool_custom = 13
integer,parameter::krome_idx_cool_co = 14
integer,parameter::krome_idx_cool_zcie = 15
integer,parameter::krome_idx_cool_zcienouv = 16
integer,parameter::krome_idx_cool_zextend = 17
integer,parameter::krome_idx_cool_gh = 18
integer,parameter::krome_ncools = 18

integer,parameter::krome_idx_heat_chem = 1
integer,parameter::krome_idx_heat_compress = 2
integer,parameter::krome_idx_heat_compr = 2
integer,parameter::krome_idx_heat_photo = 3
integer,parameter::krome_idx_heat_dh = 4
integer,parameter::krome_idx_heat_enthalpic = 4
integer,parameter::krome_idx_heat_av = 5
integer,parameter::krome_idx_heat_photoav = 5
integer,parameter::krome_idx_heat_cr = 6
integer,parameter::krome_idx_heat_dust = 7
integer,parameter::krome_idx_heat_xray = 8
integer,parameter::krome_idx_heat_viscous = 9
integer,parameter::krome_idx_heat_visc = 9
integer,parameter::krome_idx_heat_custom = 10
integer,parameter::krome_idx_heat_zcie = 11
integer,parameter::krome_nheats = 11

integer,parameter::krome_nrea=281
integer,parameter::krome_nmols=36
integer,parameter::krome_nspec=40
integer,parameter::krome_natoms=5
integer,parameter::krome_ndust=0
integer,parameter::krome_ndustTypes=0
integer,parameter::krome_nPhotoBins=7
integer,parameter::krome_nPhotoRates=17

real*8,parameter::krome_boltzmann_eV = 8.617332478d-5 !eV / K
real*8,parameter::krome_boltzmann_J = 1.380648d-23 !J / K
real*8,parameter::krome_boltzmann_erg = 1.380648d-16 !erg / K
real*8,parameter::krome_iboltzmann_eV = 1d0/krome_boltzmann_eV !K / eV
real*8,parameter::krome_iboltzmann_erg = 1d0/krome_boltzmann_erg !K / erg
real*8,parameter::krome_planck_eV = 4.135667516d-15 !eV s
real*8,parameter::krome_planck_J = 6.62606957d-34 !J s
real*8,parameter::krome_planck_erg = 6.62606957d-27 !erg s
real*8,parameter::krome_iplanck_eV = 1d0/krome_planck_eV !1 / eV / s
real*8,parameter::krome_iplanck_J = 1d0/krome_planck_J !1 / J / s
real*8,parameter::krome_iplanck_erg = 1d0/krome_planck_erg !1 / erg / s
real*8,parameter::krome_gravity = 6.674d-8 !cm3 / g / s2
real*8,parameter::krome_e_mass = 9.10938188d-28 !g
real*8,parameter::krome_p_mass = 1.67262158d-24 !g
real*8,parameter::krome_n_mass = 1.674920d-24 !g
real*8,parameter::krome_ip_mass = 1d0/krome_p_mass !1/g
real*8,parameter::krome_clight = 2.99792458e10 !cm/s
real*8,parameter::krome_pi = 3.14159265359d0 !#
real*8,parameter::krome_eV_to_erg = 1.60217646d-12 !eV -> erg
real*8,parameter::krome_ry_to_eV = 13.60569d0 !rydberg -> eV
real*8,parameter::krome_ry_to_erg = 2.179872d-11 !rydberg -> erg
real*8,parameter::krome_seconds_per_year = 365d0*24d0*3600d0 !yr -> s
real*8,parameter::krome_km_to_cm = 1d5 !km -> cm
real*8,parameter::krome_cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
real*8,parameter::krome_kvgas_erg = 8.d0*krome_boltzmann_erg/krome_pi/krome_p_mass !
real*8,parameter::krome_pre_kvgas_sqrt = sqrt(8.d0*krome_boltzmann_erg/krome_pi) !
real*8,parameter::krome_pre_planck = 2.d0*krome_planck_erg/krome_clight**2 !erg/cm2*s3
real*8,parameter::krome_exp_planck = krome_planck_erg / krome_boltzmann_erg !s*K
real*8,parameter::krome_stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
real*8,parameter::krome_N_avogadro = 6.0221d23 !#
real*8,parameter::krome_Rgas_J = 8.3144621d0 !J/K/mol
real*8,parameter::krome_Rgas_kJ = 8.3144621d-3 !kJ/K/mol
real*8,parameter::krome_hubble = 0.704d0 !dimensionless
real*8,parameter::krome_Omega0 = 1.0d0 !dimensionless
real*8,parameter::krome_Omegab = 0.0456d0 !dimensionless
real*8,parameter::krome_Hubble0 = 1.d2*krome_hubble*krome_km_to_cm*krome_cm_to_Mpc !1/s

contains

!*******************
subroutine krome_set_user_crate(argset)
use krome_commons
implicit none
real*8 :: argset
user_crate = argset
end subroutine krome_set_user_crate

!*******************
function krome_get_user_crate()
use krome_commons
implicit none
real*8 :: krome_get_user_crate
krome_get_user_crate = user_crate
end function krome_get_user_crate

!*******************
subroutine krome_set_user_Av(argset)
use krome_commons
implicit none
real*8 :: argset
user_Av = argset
end subroutine krome_set_user_Av

!*******************
function krome_get_user_Av()
use krome_commons
implicit none
real*8 :: krome_get_user_Av
krome_get_user_Av = user_Av
end function krome_get_user_Av

!*******************
subroutine krome_set_user_G0(argset)
use krome_commons
implicit none
real*8 :: argset
user_G0 = argset
end subroutine krome_set_user_G0

!*******************
function krome_get_user_G0()
use krome_commons
implicit none
real*8 :: krome_get_user_G0
krome_get_user_G0 = user_G0
end function krome_get_user_G0

!*******************
subroutine krome_set_user_gamma_CO(argset)
use krome_commons
implicit none
real*8 :: argset
user_gamma_CO = argset
end subroutine krome_set_user_gamma_CO

!*******************
function krome_get_user_gamma_CO()
use krome_commons
implicit none
real*8 :: krome_get_user_gamma_CO
krome_get_user_gamma_CO = user_gamma_CO
end function krome_get_user_gamma_CO

!*******************
subroutine krome_set_user_gamma_H2(argset)
use krome_commons
implicit none
real*8 :: argset
user_gamma_H2 = argset
end subroutine krome_set_user_gamma_H2

!*******************
function krome_get_user_gamma_H2()
use krome_commons
implicit none
real*8 :: krome_get_user_gamma_H2
krome_get_user_gamma_H2 = user_gamma_H2
end function krome_get_user_gamma_H2

!************************
!returns the Tdust averaged over the number density
! as computed in the tables
function krome_get_table_Tdust(x,Tgas)
use krome_commons
use krome_grfuncs
implicit none
real*8 :: Tgas
real*8 :: x(nmols), krome_get_table_Tdust
real*8::n(nspec)

n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas

krome_get_table_Tdust = get_table_Tdust(n(:))

end function krome_get_table_Tdust

!**********************
!convert from MOCASSIN abundances to KROME
! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
!  i=species, j=ionization level
! imap: matrix position index map, integer
! returns KROME abundances (cm-3, real*8)
function krome_convert_xmoc(xmoc,imap) result(x)
use krome_commons
use krome_subs
use krome_getphys
implicit none
real*4,intent(in):: xmoc(:,:)
real*8::x(nmols),n(nspec)
integer,intent(in)::imap(:)

x(:) = 0d0

x(idx_H) = xmoc(imap(1), 1)
x(idx_HE) = xmoc(imap(2), 1)
x(idx_C) = xmoc(imap(6), 1)
x(idx_O) = xmoc(imap(8), 1)
x(idx_Hj) = xmoc(imap(1), 2)
x(idx_HEj) = xmoc(imap(2), 2)
x(idx_Cj) = xmoc(imap(6), 2)
x(idx_Oj) = xmoc(imap(8), 2)
x(idx_HEjj) = xmoc(imap(2), 3)

n(1:nmols) = x(:)
n(nmols+1:nspec) = 0d0
x(idx_e) = get_electrons(n(:))

end function krome_convert_xmoc

!*************************
!convert from KROME abundances to MOCASSIN
! x: KROME abuances (cm-3, real*8)
! imap: matrix position index map, integer
! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
!  i=species, j=ionization level
subroutine krome_return_xmoc(x,imap,xmoc)
use krome_commons
implicit none
real*8,intent(in)::x(nmols)
real*4,intent(out)::xmoc(:,:)
integer,intent(in)::imap(:)

xmoc(:,:) = 0d0

xmoc(imap(1), 1) = x(idx_H)
xmoc(imap(2), 1) = x(idx_HE)
xmoc(imap(6), 1) = x(idx_C)
xmoc(imap(8), 1) = x(idx_O)
xmoc(imap(1), 2) = x(idx_Hj)
xmoc(imap(2), 2) = x(idx_HEj)
xmoc(imap(6), 2) = x(idx_Cj)
xmoc(imap(8), 2) = x(idx_Oj)
xmoc(imap(2), 3) = x(idx_HEjj)

end subroutine krome_return_xmoc

!**********************
!convert number density (cm-3) into column
! density (cm-2) using the specific density
! column method (see help for option
! -columnDensityMethod)
! num is the number density, x(:) is the species
! array, Tgas is the gas temperature
! If the method is not JEANS, x(:) and Tgas
! are dummy variables
function krome_num2col(num,x,Tgas)
use krome_subs
use krome_commons
use krome_getphys
implicit none
real*8 :: x(nmols),krome_num2col
real*8 :: Tgas,num
real*8::n(nspec)

n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas

krome_num2col = num2col(num,n(:))

end function krome_num2col

!***********************
!print on screen the current values of all phys variables
subroutine krome_print_phys_variables()
use krome_commons
implicit none

print *, "Tcmb:", phys_Tcmb
print *, "zredshift:", phys_zredshift
print *, "orthoParaRatio:", phys_orthoParaRatio
print *, "metallicity:", phys_metallicity
print *, "Tfloor:", phys_Tfloor

end subroutine krome_print_phys_variables

!*******************
subroutine krome_set_Tcmb(arg)
use krome_commons
implicit none
real*8 :: arg
phys_Tcmb = arg
end subroutine krome_set_Tcmb

!*******************
function krome_get_Tcmb()
use krome_commons
implicit none
real*8 :: krome_get_Tcmb
krome_get_Tcmb = phys_Tcmb
end function krome_get_Tcmb

!*******************
subroutine krome_set_zredshift(arg)
use krome_commons
implicit none
real*8 :: arg
phys_zredshift = arg
end subroutine krome_set_zredshift

!*******************
function krome_get_zredshift()
use krome_commons
implicit none
real*8 :: krome_get_zredshift
krome_get_zredshift = phys_zredshift
end function krome_get_zredshift

!*******************
subroutine krome_set_orthoParaRatio(arg)
use krome_commons
implicit none
real*8 :: arg
phys_orthoParaRatio = arg
end subroutine krome_set_orthoParaRatio

!*******************
function krome_get_orthoParaRatio()
use krome_commons
implicit none
real*8 :: krome_get_orthoParaRatio
krome_get_orthoParaRatio = phys_orthoParaRatio
end function krome_get_orthoParaRatio

!*******************
subroutine krome_set_metallicity(arg)
use krome_commons
implicit none
real*8 :: arg
phys_metallicity = arg
end subroutine krome_set_metallicity

!*******************
function krome_get_metallicity()
use krome_commons
implicit none
real*8 :: krome_get_metallicity
krome_get_metallicity = phys_metallicity
end function krome_get_metallicity

!*******************
subroutine krome_set_Tfloor(arg)
use krome_commons
implicit none
real*8 :: arg
phys_Tfloor = arg
end subroutine krome_set_Tfloor

!*******************
function krome_get_Tfloor()
use krome_commons
implicit none
real*8 :: krome_get_Tfloor
krome_get_Tfloor = phys_Tfloor
end function krome_get_Tfloor

!*******************
function krome_coolingCI(xin,inTgas)
use krome_commons
use krome_subs
use krome_cooling
use krome_constants
real*8 :: xin(nmols)
real*8 :: inTgas
real*8 :: krome_coolingCI
real*8::n(nspec),k(nZrate)
n(:) = 0d0
n(idx_Tgas) = inTgas
n(1:nmols) = xin(:)
k(:) = coolingZ_rate_tabs(inTgas)
krome_coolingCI = coolingCI(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
end function krome_coolingCI

!*******************
function krome_coolingCII(xin,inTgas)
use krome_commons
use krome_subs
use krome_cooling
use krome_constants
real*8 :: xin(nmols)
real*8 :: inTgas
real*8 :: krome_coolingCII
real*8::n(nspec),k(nZrate)
n(:) = 0d0
n(idx_Tgas) = inTgas
n(1:nmols) = xin(:)
k(:) = coolingZ_rate_tabs(inTgas)
krome_coolingCII = coolingCII(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
end function krome_coolingCII

!*******************
function krome_coolingOI(xin,inTgas)
use krome_commons
use krome_subs
use krome_cooling
use krome_constants
real*8 :: xin(nmols)
real*8 :: inTgas
real*8 :: krome_coolingOI
real*8::n(nspec),k(nZrate)
n(:) = 0d0
n(idx_Tgas) = inTgas
n(1:nmols) = xin(:)
k(:) = coolingZ_rate_tabs(inTgas)
krome_coolingOI = coolingOI(n(:),n(idx_Tgas),k(:)) *  boltzmann_erg
end function krome_coolingOI

!***************************
!dump the population of the Z cooling levels
! in the nfile file unit, using xvar as
! independent variable. alias of
! dump_cooling_pop subroutine
subroutine krome_popcool_dump(xvar,nfile)
use krome_cooling
implicit none
real*8 :: xvar
integer :: nfile

call dump_cooling_pop(xvar,nfile)

end subroutine krome_popcool_dump

!************************************
! when using 3D variables set the value of the 3rd dimension
! (usually Av). The other 2 dimensions are ntot and Tgas
! (take from internal KROME)
subroutine krome_set_dust_table_3D_variable(variable_value)
use krome_commons
implicit none
real*8,intent(in)::variable_value

dust_table_AvVariable_log = log10(variable_value)

end subroutine krome_set_dust_table_3D_variable

!*****************************
!get averaged Tdust from tables, with x(:) species array
! of size krome_nmols, and Tgas the gas temperature
function krome_get_Tdust(x,Tgas)
use krome_commons
use krome_fit
implicit none
real*8 :: x(nmols), krome_get_Tdust
real*8 :: Tgas
real*8 :: ntot

ntot = sum(x(1:nmols))
krome_get_Tdust = 1d1**fit_anytab2D(dust_tab_ngas(:), dust_tab_Tgas(:), &
    dust_tab_Tdust(:,:), dust_mult_ngas, dust_mult_Tgas, &
    log10(ntot), log10(Tgas))

end function krome_get_Tdust

!*****************************
!dump the data for restart (UNDER DEVELOPEMENT!)
!arguments: the species array and the gas temperature
subroutine krome_store(x,Tgas,dt)
use krome_commons
implicit none
integer::nfile,i
real*8 :: x(nmols)
real*8 :: Tgas,dt

nfile = 92

open(nfile,file="krome_dump.dat",status="replace")
!dump temperature
write(nfile,*) Tgas
write(nfile,*) dt
!dump species
do i=1,nmols
write(nfile,*) x(i)
end do
close(nfile)

end subroutine krome_store

!*****************************
!restore the data from a dump (UNDER DEVELOPEMENT!)
!arguments: the species array and the gas temperature
subroutine krome_restore(x,Tgas,dt)
use krome_commons
implicit none
integer::nfile,i
real*8 :: x(nmols)
real*8 :: Tgas,dt

nfile = 92

open(nfile,file="krome_dump.dat",status="old")
!restore temperature
read(nfile,*) Tgas
read(nfile,*) dt
!restore species
do i=1,nmols
read(nfile,*) x(i)
end do
close(nfile)

end subroutine krome_restore

!****************************
!switch on the thermal calculation
subroutine krome_thermo_on()
use krome_commons
krome_thermo_toggle = 1
end subroutine krome_thermo_on

!****************************
!switch off the thermal calculation
subroutine krome_thermo_off()
use krome_commons
krome_thermo_toggle = 0
end subroutine krome_thermo_off

!************************
! prepares tables for cross sections and
! photorates
subroutine krome_calc_photobins()
use krome_photo
call calc_photobins()
end subroutine krome_calc_photobins

!****************************
! set the energy per photo bin
! eV/cm2/sr
subroutine krome_set_photoBinJ(phbin)
use krome_commons
use krome_photo
implicit none
real*8 :: phbin(nPhotoBins)
photoBinJ(:) = phbin(:)
photoBinJ_org(:) = phbin(:) !for restore

!compute rates
call calc_photobins()

end subroutine krome_set_photoBinJ

!*************************
! set the energy (frequency) of the photobin
! as left-right limits in eV
subroutine krome_set_photobinE_lr(phbinleft,phbinright,Tgas)
use krome_commons
use krome_photo
implicit none
real*8 :: phbinleft(nPhotoBins),phbinright(nPhotoBins)
real*8,optional::Tgas
real*8::bTgas

!default Tgas for broadening
bTgas = 1d1
if(present(Tgas)) then
bTgas = Tgas
end if

!$omp parallel
photoBinEleft(:) = phbinleft(:)
photoBinEright(:) = phbinright(:)
photoBinEmid(:) = 0.5d0*(phbinleft(:)+phbinright(:))
photoBinEdelta(:) = phbinright(:)-phbinleft(:)
photoBinEidelta(:) = 1d0/photoBinEdelta(:)
!$omp end parallel

!initialize xsecs table
call init_photoBins(bTgas)

end subroutine krome_set_photobinE_lr

!*************************
! set the energy (frequency) of photobins
! when contiguous. Left and right limits are automatically
! extracted. Energy in eV
subroutine krome_set_photobinE_limits(phbinLimits,Tgas)
use krome_commons
use krome_photo
implicit none
real*8 :: phbinLimits(nPhotoBins+1)
real*8,optional::Tgas
real*8::phl(nPhotoBins),phr(nPhotoBins),bTgas

!default Tgas for broadening
bTgas = 1d1
if(present(Tgas)) then
bTgas = Tgas
end if
phl(:) = phbinLimits(1:nPhotoBins)
phr(:) = phbinLimits(2:nPhotoBins+1)

call krome_set_photobinE_lr(phl(:),phr(:),bTgas)

end subroutine krome_set_photobinE_limits

!*******************************
!set the energy (eV) of the photobin according
! to MOCASSIN way (position and width array)
subroutine krome_set_photobinE_moc(binPos,binWidth,Tgas)
use krome_commons
use krome_photo
implicit none
real*8 :: binPos(nPhotoBins),binWidth(nPhotoBins)
real*8,optional::Tgas
real*8::bTgas

bTgas = 1d1
if(present(Tgas)) then
bTgas = Tgas
end if

!$omp parallel
photoBinEleft(:) = binPos(:)-binWidth(:)/2d0
photoBinEright(:) = binPos(:)+binWidth(:)/2d0
photoBinEmid(:) = binPos(:)
photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
photoBinEidelta(:) = 1d0/photoBinEdelta(:)
!$omp end parallel

!initialize xsecs table
call init_photoBins(bTgas)

end subroutine krome_set_photobinE_moc

!********************************
! set the energy (eV) of the photobin
! linearly from lowest to highest energy value
! in eV
subroutine krome_set_photobinE_lin(lower,upper,Tgas)
use krome_commons
use krome_photo
implicit none
real*8 :: lower,upper
real*8,optional::Tgas
real*8::dE,bTgas
integer::i

bTgas = 1d1
if(present(Tgas)) then
bTgas = Tgas
end if

!$omp parallel
dE = abs(upper-lower)/nPhotoBins
!$omp end parallel
do i=1,nPhotoBins
!$omp parallel
photoBinEleft(i) = dE*(i-1) + lower
photoBinEright(i) = dE*i + lower
photoBinEmid(i) = 0.5d0*(photoBinEleft(i)+photoBinEright(i))
!$omp end parallel
end do
!$omp parallel
photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
photoBinEidelta(:) = 1d0/photoBinEdelta(:)
!$omp end parallel

!initialize xsecs table
call init_photoBins(bTgas)

end subroutine krome_set_photobinE_lin

!********************************
! set the energy (eV) of the photobin
! logarithmically from lowest to highest energy value
! in eV
subroutine krome_set_photobinE_log(lower,upper,Tgas)
use krome_commons
use krome_photo
implicit none
real*8 :: lower,upper
real*8,optional::Tgas
real*8::dE,logup,loglow,bTgas
integer::i

bTgas = 1d1
if(present(Tgas)) then
bTgas = Tgas
end if

if(lower.ge.upper) then
print *,"ERROR: in  krome_set_photobinE_log lower >= upper limit!"
stop
end if
loglow = log10(lower)
logup = log10(upper)
!$omp parallel
dE = 1d1**(abs(logup-loglow)/nPhotoBins)
!$omp end parallel
do i=1,nPhotoBins
!$omp parallel
photoBinEleft(i) = 1d1**((i-1)*(logup-loglow)/nPhotoBins + loglow)
photoBinEright(i) = 1d1**(i*(logup-loglow)/nPhotoBins + loglow)
photoBinEmid(i) = 0.5d0*(photoBinEleft(i)+photoBinEright(i))
!$omp end parallel
end do
!$omp parallel
photoBinEdelta(:) = photoBinEright(:)-photoBinEleft(:)
photoBinEidelta(:) = 1d0/photoBinEdelta(:)
!$omp end parallel

!initialize xsecs table
call init_photoBins(bTgas)

end subroutine krome_set_photobinE_log

!*********************************
!returns an array containing the flux for each photo bin
! in eV/cm2/sr
function krome_get_photoBinJ()
use krome_commons
real*8 :: krome_get_photoBinJ(nPhotoBins)
krome_get_photoBinJ(:) = photoBinJ(:)
end function krome_get_photoBinJ

!*********************************
!get an array containing all the left positions
! of the photobins, eV
function krome_get_photoBinE_left()
!returns an array of size krome_nPhotoBins with the
! left energy limits (eV)
use krome_commons
real*8 :: krome_get_photoBinE_left(nPhotoBins)
krome_get_photoBinE_left(:) = photoBinEleft(:)
end function krome_get_photoBinE_left

!*********************************
!returns an array of size krome_nPhotoBins with the
! right energy limits (eV)
function krome_get_photoBinE_right()
use krome_commons
real*8 :: krome_get_photoBinE_right(nPhotoBins)
krome_get_photoBinE_right(:) = photoBinEright(:)
end function krome_get_photoBinE_right

!*********************************
!returns an array of size krome_nPhotoBins with the
! middle energy values (eV)
function krome_get_photoBinE_mid()
use krome_commons
real*8 :: krome_get_photoBinE_mid(nPhotoBins)
krome_get_photoBinE_mid(:) = photoBinEmid(:)
end function krome_get_photoBinE_mid

!*********************************
!returns an array of size krome_nPhotoBins with the
! bin span (eV)
function krome_get_photoBinE_delta()
use krome_commons
real*8 :: krome_get_photoBinE_delta(nPhotoBins)
krome_get_photoBinE_delta(:) = photoBinEdelta(:)
end function krome_get_photoBinE_delta

!*********************************
!returns an array of size krome_nPhotoBins with the
! inverse of the bin span (1/eV)
function krome_get_photoBinE_idelta()
use krome_commons
real*8 :: krome_get_photoBinE_idelta(nPhotoBins)
krome_get_photoBinE_idelta(:) = photoBinEidelta(:)
end function krome_get_photoBinE_idelta

!*********************************
!returns an array of size krome_nPhotoBins with the
! integrated photo rates (1/s)
function krome_get_photoBin_rates()
use krome_commons
real*8 :: krome_get_photoBin_rates(nPhotoRea)
krome_get_photoBin_rates(:) = photoBinRates(:)
end function krome_get_photoBin_rates

!*********************************
!returns an array of size krome_nPhotoBins containing
! the cross section (cm2) of the idx-th photoreaction
function krome_get_xsec(idx)
use krome_commons
implicit none
real*8 :: krome_get_xsec(nPhotoBins)
integer :: idx

krome_get_xsec(:) = photoBinJTab(idx,:)

end function krome_get_xsec

!*********************************
!returns an array of size krome_nPhotoBins with the
! integrated photo heatings (erg/s)
function krome_get_photoBin_heats()
use krome_commons
implicit none
real*8 :: krome_get_photoBin_heats(nPhotoRea)
krome_get_photoBin_heats(:) = photoBinHeats(:)

end function krome_get_photoBin_heats

!****************************
!multiply all photobins by a factor real*8 xscale
subroutine krome_photoBin_scale(xscale)
use krome_commons
use krome_photo
implicit none
real*8 :: xscale

photoBinJ(:) = photoBinJ(:) * xscale

!compute rates
call calc_photobins()

end subroutine krome_photoBin_scale

!****************************
!multiply all photobins by a real*8 array xscale(:)
! of size krome_nPhotoBins
subroutine krome_photoBin_scale_array(xscale)
use krome_commons
use krome_photo
implicit none
real*8 :: xscale(nPhotoBins)

photoBinJ(:) = photoBinJ(:) * xscale(:)

!compute rates
call calc_photobins()

end subroutine krome_photoBin_scale_array

!********************************
!restore the original flux (i.e. undo any rescale).
! the flux is automatically stored by the functions
! that set the flux, or by the function
! krome_photoBin_store()
subroutine krome_photoBin_restore()
use krome_commons
implicit none

photoBinJ(:) = photoBinJ_org(:)

end subroutine krome_photoBin_restore

!**********************
!store flux to be restored with the subroutine
! krome_photoBin_restore later
subroutine krome_photoBin_store()
use krome_commons
implicit none

photoBinJ_org(:) = photoBinJ(:)

end subroutine krome_photoBin_store

!*********************
!load flux radiation from a two-columns file
! energy/eV, flux/(eV/cm2/sr)
! Flux is interpolated over the existing binning
! constant-area method
subroutine krome_load_photoBin_file_2col(fname, logarithmic)
use krome_commons
implicit none
integer,parameter::imax=int(1e4)
character(len=*) :: fname
logical, optional :: logarithmic
logical :: is_log
integer::unit,ios,icount,j,i
real*8::xtmp(imax),ftmp(imax),intA,eL,eR
real*8::xL,xR,pL,pR,fL,fR,Jflux(nPhotoBins)
real*8::a,b

if(present(logarithmic)) then
is_log = logarithmic
else
is_log = .false.
end if

!open file to read
open(newunit=unit,file=trim(fname),iostat=ios)
if(ios/=0) then
print *,"ERROR: problems reading "//trim(fname)
stop
end if

!read file line by line and store to temporary
ftmp(:) = 0d0
icount = 1
do
read(unit,*,iostat=ios) xtmp(icount), ftmp(icount)
if(ios/=0) exit
icount = icount + 1
end do
close(unit)
icount = icount - 1

if(is_log) ftmp = log(merge(ftmp,1d-40,ftmp>0d0))
!loop on photobins for interpolation
do j=1,nPhotoBins
intA = 0d0
!photobin limits
eL = photoBinEleft(j)
eR = photoBinEright(j)
!loop on flux bins
do i=1,icount-1
!flux bin limits
xL = xtmp(i)
xR = xtmp(i+1)
!if outside the bin skip
if((xR<eL).or.(xL>eR)) cycle
!get the interval limit (consider partial overlapping)
pL = max(xL,eL)
pR = min(xR,eR)
if(is_log) then
pL = log(max(pL,1d-40))
pR = log(max(pR,1d-40))
xL = log(max(xL,1d-40))
xR = log(max(xR,1d-40))
end if
!interpolate to get the flux at the interval limit
fL = (ftmp(i+1)-ftmp(i))*(pL-xL)/(xR-xL)+ftmp(i)
fR = (ftmp(i+1)-ftmp(i))*(pR-xL)/(xR-xL)+ftmp(i)
if(is_log) then
fL = exp(fL)
fR = exp(fR)
pL = exp(pL)
pR = exp(pR)
end if

!compute area of the overlapped area
intA = intA + (fL+fR)*(pR-pL)/2d0
end do
!distribute the flux in the photobin
Jflux(j) = intA / (eR-eL)
end do

!initialize intensity according to data
call krome_set_photoBinJ(Jflux(:))

end subroutine krome_load_photoBin_file_2col

!********************************
!load the radiation bins from the file fname
! data should be a 3-column file with
! energy Left (eV), energy Right (eV)
! intensity (eV/cm2/sr).
! This subroutine sets also the bin-size
subroutine krome_load_photoBin_file(fname) !! !! not yet callable from C
use krome_commons
implicit none
integer::ios,icount
character(len=*) :: fname
real*8::tmp_El(nPhotoBins),tmp_Er(nPhotoBins)
real*8::rout(3),tmp_J(nPhotoBins)

!open file and check for errors
open(33,file=fname,status="old",iostat=ios)
if(ios.ne.0) then
print *,"ERROR: problem opening "//fname//"!"
print *," (e.g. file not found)"
stop
end if

icount = 0 !count valid line
!loop on file
do
read(33,*,iostat=ios) rout(:)
if(ios==-1) exit !EOF
if(ios.ne.0) cycle !skip comments
icount = icount + 1
if(icount>nPhotoBins) exit !can't load more than nPhotoBins
tmp_El(icount) = rout(1) !energy L eV
tmp_Er(icount) = rout(2) !energy R eV
!check if left interval is before right
if(tmp_El(icount)>tmp_Er(icount)) then
print *,"ERROR: in file "//fname//" left"
print *, " interval larger than right one!"
print *,tmp_El(icount),tmp_Er(icount)
stop
end if
tmp_J(icount) = rout(3) !intensity eV/cm2/sr
end do
close(33)

!file data lines should be the same number of the photobins
if(icount/=nPhotoBins) then
print *,"ERROR: the number of data lines in the file"
print *," "//fname//" should be equal to the number of"
print *," photobins ",nPhotoBins
print *,"Found",icount
stop
end if

!initialize interval and intensity according to data
call krome_set_photobinE_lr(tmp_El(:),tmp_Er(:))
call krome_set_photoBinJ(tmp_J(:))

end subroutine krome_load_photoBin_file

!**********************************
!this subroutine sets an Hardt+Madau flux in the
! energy limits lower_in, upper_in (eV, log-spaced)
subroutine krome_set_photoBin_HMlog(lower_in,upper_in)
use krome_commons
use krome_photo
use krome_subs
use krome_fit
implicit none
real*8::z(59),energy(500),HM(59,500)
real*8::z_mul,energy_mul,x,lower,upper
real*8,parameter::limit_lower = 0.1237d0
real*8,parameter::limit_upper = 4.997d7
real*8,parameter::limit_redshift = 15.660d0
real*8,optional::lower_in,upper_in
integer::i

lower = limit_lower
upper = limit_upper
if(present(lower_in)) lower = lower_in
if(present(upper_in)) upper = upper_in

if(phys_zredshift>limit_redshift) then
print *,"ERROR: redshift out of range in HM"
print *,"redshift:",phys_zredshift
print *,"limit:",limit_redshift
stop
end if

if(lower<limit_lower .or. upper>limit_upper) then
print *,"ERROR: upper or lower limit out of range in HM."
print *,"lower limit (eV):",limit_lower
print *,"upper limit (eV):",limit_upper
stop
end if

call krome_set_photoBinE_log(lower,upper)

call init_anytab2D("krome_HMflux.dat", z(:), energy(:), &
    HM(:,:), z_mul, energy_mul)

do i=1,nPhotoBins
x = log10(photoBinEmid(i)) !log(eV)
photoBinJ(i) = 1d1**fit_anytab2D(z(:), energy(:), HM(:,:), &
    z_mul, energy_mul, phys_zredshift, x)
end do

photoBinJ_org(:) = photoBinJ(:)

call calc_photobins()

end subroutine krome_set_photoBin_HMlog

!**********************************
!this subroutine ADD an Hardt+Madau flux to the current radiation
! in the energy limits lower_in, upper_in (eV), It assumes
! the current binning
subroutine krome_set_photoBin_HMCustom(lower_in,upper_in,additive)
use krome_commons
use krome_photo
use krome_subs
use krome_fit
implicit none
real*8::z(59),energy(500),HM(59,500)
real*8::z_mul,energy_mul,x,lower,upper
real*8::photoTmpJ(nPhotoBins)
real*8,parameter::limit_lower = 0.1237d0
real*8,parameter::limit_upper = 4.997d7
real*8,parameter::limit_redshift = 15.660d0
logical,optional::additive
logical::add
real*8,optional :: lower_in,upper_in
integer::i

lower = limit_lower
upper = limit_upper
if(present(lower_in)) lower = lower_in
if(present(upper_in)) upper = upper_in

add = .false.
if(present(additive)) add = additive

if(phys_zredshift>limit_redshift) then
print *,"ERROR: redshift out of range in HM"
print *,"redshift:",phys_zredshift
print *,"limit:",limit_redshift
stop
end if

if(lower<limit_lower .or. upper>limit_upper) then
print *,"ERROR: upper or lower limit out of range in HM."
print *,"lower limit (eV):",limit_lower
print *,"upper limit (eV):",limit_upper
stop
end if

call init_anytab2D("krome_HMflux.dat", z(:), energy(:), &
    HM(:,:), z_mul, energy_mul)

do i=1,nPhotoBins
x = log10(photoBinEmid(i)) !log(eV)
photoTmpJ(i) = 1d1**fit_anytab2D(z(:), energy(:), HM(:,:), &
    z_mul, energy_mul, phys_zredshift, x)
end do

!add flux to already-present flux if optional argument
if(add) then
photoBinJ(:) = photoBinJ(:) + photoTmpJ(:)
else
photoBinJ(:) = photoTmpJ(:)
end if

photoBinJ_org(:) = photoBinJ(:)

call calc_photobins()

end subroutine krome_set_photoBin_HMCustom

!**********************************
!set the flux as a black body with temperature Tbb (K)
! in the range lower to upper (eV),  linear-spaced
subroutine krome_set_photoBin_BBlin(lower,upper,Tbb)
use krome_commons
use krome_constants
use krome_photo
use krome_subs
use krome_phfuncs
implicit none
real*8 :: lower,upper,Tbb
real*8::x
integer::i

call krome_set_photoBinE_lin(lower,upper)

!eV/cm2/sr
do i=1,nPhotoBins
x = photoBinEmid(i) !eV
photoBinJ(i) = planckBB(x,Tbb)
end do
photoBinJ_org(:) = photoBinJ(:)

call calc_photobins()

end subroutine krome_set_photoBin_BBlin

!**********************************
!set the flux as a black body with temperature Tbb (K)
! in the range lower to upper (eV), log-spaced
subroutine krome_set_photoBin_BBlog(lower,upper,Tbb)
use krome_commons
use krome_constants
use krome_photo
use krome_subs
use krome_phfuncs
implicit none
real*8 :: lower,upper,Tbb
real*8::x,xmax,xexp,Jlim
integer::i

!limit for the black body intensity to check limits
Jlim = 1d-3

call krome_set_photoBinE_log(lower,upper)

!eV/cm2/sr
do i=1,nPhotoBins
x = photoBinEmid(i) !eV
photoBinJ(i) = planckBB(x,Tbb)
end do
photoBinJ_org(:) = photoBinJ(:)

!uncomment this below for additional control
!!$    !find the maximum using Wien's displacement law
!!$    xmax = Tbb/2.8977721d-1 * clight * planck_eV !eV
!!$
!!$    if(xmax<lower) then
!!$       print *,"WARNING: maximum of the Planck function"
!!$       print *," is below the lowest energy bin!"
!!$       print *,"max (eV)",xmax
!!$       print *,"lowest (eV)",lower
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(xmax>upper) then
!!$       print *,"WARNING: maximum of the Planck function"
!!$       print *," is above the highest energy bin!"
!!$       print *,"max (eV)",xmax
!!$       print *,"highest (eV)",upper
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(photoBinJ(1)>Jlim) then
!!$       print *,"WARNING: lower bound of the Planck function"
!!$       print *," has a flux of (ev/cm2/s/Hz/sr)",photoBinJ(1)
!!$       print *," which is larger than the limit Jlim",Jlim
!!$       print *,"Tbb (K)",Tbb
!!$    end if
!!$
!!$    if(photoBinJ(nPhotoBins)>Jlim) then
!!$       print *,"WARNING: upper bound of the Planck function"
!!$       print *," has a flux of (ev/cm2/s/Hz/sr)",photoBinJ(nPhotoBins)
!!$       print *," which is larger than the limit Jlim",Jlim
!!$       print *,"Tbb (K)",Tbb
!!$    end if

!compute rates
call calc_photobins()

end subroutine krome_set_photoBin_BBlog

!*************************************
!set the BB spectrum and the limits using bisection
subroutine krome_set_photoBin_BBlog_auto(Tbb)
use krome_commons
use krome_subs
use krome_constants
use krome_phfuncs
implicit none
real*8 :: Tbb
real*8::xlow,xup,eps,xmax,J0,J1,x0,x1,xm,Jm
eps = 1d-6

!Rayleigh–Jeans approximation for the minimum energy
xlow = planck_eV*clight*sqrt(.5d0/Tbb/boltzmann_eV*eps)

!find energy of the Wien maximum (eV)
xmax = Tbb / 2.8977721d-1 * clight * planck_eV

!bisection to find the maximum
x0 = xmax
x1 = 2.9d2*Tbb*boltzmann_eV
J0 = planckBB(x0,Tbb) - eps
J1 = planckBB(x1,Tbb) - eps
if(J0<0d0.or.J1>0d0) then
print *,"ERROR: problems with auto planck bisection!"
stop
end if

do
xm = 0.5d0*(x0+x1)
Jm = planckBB(xm,Tbb) - eps
if(Jm>0d0) x0 = xm
if(Jm<0d0) x1 = xm
if(abs(Jm)<eps*1d-3) exit
end do
xup = xm

!initialize BB radiation using the values found
call krome_set_photoBin_BBlog(xlow,xup,Tbb)

end subroutine krome_set_photoBin_BBlog_auto

!*********************************
!return the ratio between the current flux an Draine's
function krome_get_ratioFluxDraine()
use krome_subs
use krome_phfuncs
implicit none
real*8::krome_get_ratioFluxDraine

krome_get_ratioFluxDraine = get_ratioFluxDraine()

end function krome_get_ratioFluxDraine

!**********************************
!set the flux as Draine's function
! in the range lower to upper (eV). the spacing is linear
subroutine krome_set_photoBin_draineLin(lower,upper)
use krome_commons
use krome_photo
use krome_constants
real*8 :: upper,lower
real*8::x
integer::i

call krome_set_photoBinE_lin(lower,upper)

do i=1,nPhotoBins
x = photoBinEmid(i) !eV
!eV/cm2/sr
if(x<13.6d0.and.x>5d0) then
photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
    * x *planck_eV
else
photoBinJ(i) = 0d0
end if
end do

photoBinJ_org(:) = photoBinJ(:)

!compute rates
call calc_photobins()

end subroutine krome_set_photoBin_draineLin

!**************************
!set the flux as Draine's function
! in the range lower to upper (eV). log-spaced
subroutine krome_set_photoBin_draineLog(lower,upper)
use krome_commons
use krome_photo
use krome_constants
real*8 :: upper,lower
real*8::x
integer::i

call krome_set_photoBinE_log(lower,upper)

do i=1,nPhotoBins
x = photoBinEmid(i) !eV
!eV/cm2/sr/s/Hz
if(x<13.6d0.and.x>5d0) then
photoBinJ(i) = (1.658d6*x - 2.152d5*x**2 + 6.919d3*x**3) &
    * x *planck_eV
else
photoBinJ(i) = 0d0
end if
end do

photoBinJ_org(:) = photoBinJ(:)

!compute rates
call calc_photobins()

end subroutine krome_set_photoBin_draineLog

!**************************
!set the flux as Draine's function with the current binning
! Note: you have to set the binning first
subroutine krome_set_photoBin_draineCustom()
use krome_commons
use krome_photo
use krome_constants
real*8::xL,xR,f1,f2
integer::i

!return error if binning is not set
if(maxval(photoBinEmid)==0d0) then
print *,"ERROR: not initialized binning in draineCustom!"
stop
end if

!loop on bins
do i=1,nPhotoBins
!eV/cm2/sr
if(xR<=13.6d0.and.xL>=5d0) then
xL = photoBinEleft(i) !eV
xR = photoBinEright(i) !eV
elseif(xL<5d0.and.xR>5d0) then
xL = 5d0 !eV
xR = photoBinEright(i) !eV
elseif(xL<13d0.and.xR>13.6d0) then
xL = photoBinEleft(i) !eV
xR = 13.6d0 !eV
else
xL = 0d0
xR = 0d0
end if
f1 = (1.658d6*xL - 2.152d5*xL**2 + 6.919d3*xL**3) &
    * planck_eV
f2 = (1.658d6*xR - 2.152d5*xR**2 + 6.919d3*xR**3) &
    * planck_eV
photoBinJ(i) = (f1+f2)*(xR-xL)/2d0

end do

photoBinJ_org(:) = photoBinJ(:)

!compute rates
call calc_photobins()

end subroutine krome_set_photoBin_draineCustom

!**************************
!set the flux as power-law (J21-style)
! in the range lower to upper (eV). linear-spaced
subroutine krome_set_photoBin_J21lin(lower,upper)
use krome_commons
use krome_photo
real*8 :: upper,lower

call krome_set_photoBinE_lin(lower,upper)

photoBinJ(:) = 6.2415d-10 * (13.6d0/photoBinEmid(:)) !eV/cm2/sr
photoBinJ_org(:) = photoBinJ(:)

!compute rates
call calc_photobins()

end subroutine krome_set_photoBin_J21lin

!**************************
!set the flux as power-law (J21-style)
! in the range lower to upper (eV). the spacing is logarithmic
subroutine krome_set_photoBin_J21log(lower,upper)
use krome_commons
use krome_photo
real*8 :: upper,lower

call krome_set_photoBinE_log(lower,upper)

photoBinJ(:) = 6.2415d-10 * (13.6d0/photoBinEmid(:)) !eV/cm2/sr
photoBinJ_org(:) = photoBinJ(:)

!compute rates
call calc_photobins()

end subroutine krome_set_photoBin_J21log

!*****************************
!get the opacity tau corresponding to x(:)
! chemical composition. The column density
! is computed using the expression in the
! num2col(x) function.
! An array of size krome_nPhotoBins is returned
function krome_get_opacity(x,Tgas)
use krome_commons
use krome_constants
use krome_photo
use krome_subs
use krome_getphys
implicit none
real*8 :: x(nmols),krome_get_opacity(nPhotoBins)
real*8,value :: Tgas
real*8::tau,n(nspec)
integer::i,j,idx

n(1:nmols) = x(:)
n(idx_Tgas) = Tgas

!loop on frequency bins
do j=1,nPhotoBins
tau = 0d0
!loop on species
do i=1,nPhotoRea
!calc opacity as column_density * cross_section
idx = photoPartners(i)
tau = tau + num2col(x(idx),n(:)) * photoBinJTab(i,j)
end do
krome_get_opacity(j) = tau !store
end do

end function krome_get_opacity

!*****************************
!get the opacity tau corresponding to the x(:)
! chemical composition. The column density
! is computed using the size of the cell (csize)
! An array of size krome_nPhotoBins is returned.
function krome_get_opacity_size(x,Tgas,csize)
use krome_commons
use krome_constants
use krome_photo
use krome_subs
use krome_dust
implicit none
real*8 :: x(nmols),krome_get_opacity_size(nPhotoBins)
real*8,value :: Tgas,csize
real*8::n(nspec),energy,tau
integer::i,j,idx

n(1:nmols) = x(:)
n(idx_Tgas) = Tgas

!loop on frequency bins
do j=1,nPhotoBins
tau = 0d0
!loop on species
do i=1,nPhotoRea
!calc opacity as column_density * cross_section
!where column_density is density*cell_size
idx = photoPartners(i)
tau = tau + x(idx) * photoBinJTab(i,j)
end do

krome_get_opacity_size(j) = tau * csize !store
end do

end function krome_get_opacity_size

!*****************************
!get the opacity tau corresponding to the x(:)
! chemical composition. The column density
! is computed using the size of the cell (csize).
! Dust is included using dust-to-gas mass ratio (d2g).
! You should load the dust tables with the subroutine
! krome_load_dust_opacity(fileName).
! An array of size krome_nPhotoBins is returned.
function krome_get_opacity_size_d2g(x,Tgas,csize,d2g)
use krome_commons
use krome_constants
use krome_photo
use krome_subs
use krome_dust
use krome_getphys
implicit none
real*8 :: x(nmols),krome_get_opacity_size_d2g(nPhotoBins)
real*8,value :: Tgas,csize,d2g
real*8::n(nspec),energy,tau,m(nspec),mgas
integer::i,j,idx

m(:) = get_mass()
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas
mgas = sum(n(1:nmols)*m(1:nmols))

!loop on frequency bins
do j=1,nPhotoBins
tau = 0d0
!loop on species
do i=1,nPhotoRea
!calc opacity as column_density * cross_section
!where column_density is density*cell_size
idx = photoPartners(i)
tau = tau + x(idx) * photoBinJTab(i,j)
end do

!sum dust opacity from interpolated table
tau = tau + d2g*mgas*opacityDust(j)

krome_get_opacity_size_d2g(j) = tau * csize !store
end do

end function krome_get_opacity_size_d2g

!*********************
!scale radiation intensity with opacity assuming a given
! cell size and gas composition
!  subroutine krome_opacity_scale_size(csize,n,Tgas)
!    use krome_commons
!    implicit none
!    real*8::csize,n(nmols),xscale(nPhotoBins),Tgas
!
!    xscale(:) = krome_get_opacity_size(n(:),Tgas,csize)
!    xscale(:) = exp(-xscale(:))
!    call krome_photoBin_scale_array(xscale(:))
!
!  end subroutine krome_opacity_scale_size

!*********************
!scale radiation intensity with opacity assuming a given
! cell size and gas composition
subroutine krome_opacity_scale_size(csize,n,Tgas)
use krome_commons
implicit none
real*8::csize,n(nmols),Tgas
real*8 :: xscale(nPhotoBins)

xscale(:) = krome_get_opacity_size(n(:),Tgas,csize)
xscale(:) = exp(-xscale(:))
call krome_photoBin_scale_array(xscale(:))
end subroutine krome_opacity_scale_size

!*******************************
!load a frequency-dependent opacity table stored in fname file,
! column 1 is energy or wavelenght in un units of unitEnergy
! (default eV), column 2 is opacity in cm2/g.
! opacity is interpolated over the current photo-binning.
subroutine krome_load_opacity_table(fname, unitEnergy)
use krome_commons
use krome_constants
use krome_photo
implicit none
integer,parameter::ntmp=int(1e5)
character(len=*)::fname
character(len=*),optional::unitEnergy
character*10::eunit
integer::ios,icount,iR,iL,i,j,fileUnit
real*8::wl,opac,fL,fR,kk,dE
real*8::wls(ntmp),opacs(ntmp)
real*8,allocatable::energy(:),kappa(:)

!read energy unit optional argument
eunit = "eV" !default is eV
if(present(unitEnergy)) then
eunit = trim(unitEnergy)
end if

call load_opacity_table(fname, eunit)

end subroutine krome_load_opacity_table

! ******************************
! load absorption data data from file, cm2/g
subroutine krome_load_average_kabs()
use krome_photo
implicit none

call find_Av_load_kabs()

end subroutine krome_load_average_kabs

! *******************************
! use linear least squares and the current Jflux distribution
! to return G0 and Av.
! x(:) are the abundances (use for mean molecular weight)
! and d2g is dust to gas mass ratio
subroutine krome_find_G0_Av(G0, Av, x, d2g)
use krome_commons
use krome_photo
implicit none
real*8,intent(out)::G0, Av
real*8,intent(in)::d2g, x(nmols)
real*8::n(nspec)

n(1:nmols) = x(:)
n(nmols+1:nspec) = 0d0

call estimate_G0_Av(G0, Av, n(:), d2g)

end subroutine krome_find_G0_Av

!*******************************
!dump the Jflux profile to the file
! with unit number nfile
subroutine krome_dump_Jflux(nfile)
use krome_commons
implicit none
integer::i
integer :: nfile

do i=1,nPhotoBins
write(nfile,*) photoBinEmid(i),photoBinJ(i)
end do

end subroutine krome_dump_Jflux

!**********************
!set the velocity for line broadening, cm/s
subroutine krome_set_lineBroadeningVturb(vturb)
use krome_constants
use krome_commons
implicit none
real*8::vturb

broadeningVturb2 = vturb**2

end subroutine krome_set_lineBroadeningVturb

!***************************
!alias for coe in krome_subs
! returns the coefficient array of size krome_nrea
! for a given Tgas
function krome_get_coef(Tgas,x)
use krome_commons
use krome_subs
use krome_tabs
real*8 :: krome_get_coef(nrea),x(nmols)
real*8,value:: Tgas
real*8::n(nspec)
n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas

call makeStoreOnceRates(n(:))

krome_get_coef(:) = coe(n(:))

end function krome_get_coef

!****************************
!get the mean molecular weight from
! mass fractions
function krome_get_mu_x(xin)
use krome_commons
implicit none
real*8 :: xin(nmols), krome_get_mu_x
real*8::n(nmols)
n(:) = krome_x2n(xin(:),1d0)
krome_get_mu_x = krome_get_mu(n(:))
end function krome_get_mu_x

!****************************
!return the adiabatic index from mass fractions
! and temperature in K
function krome_get_gamma_x(xin,inTgas)
use krome_commons
implicit none
real*8 :: inTgas
real*8 :: xin(nmols), krome_get_gamma_x
real*8::x(nmols),Tgas,rhogas

Tgas = inTgas
x(:) = krome_x2n(xin(:),1d0)
krome_get_gamma_x = krome_get_gamma(x(:),Tgas)

end function krome_get_gamma_x

!***************************
!normalize mass fractions and
! set charge to zero
subroutine krome_consistent_x(x)
use krome_commons
use krome_constants
implicit none
real*8 :: x(nmols)
real*8::isumx,sumx,xerr,imass(nmols),ee

!1. charge consistency
imass(:) = krome_get_imass()

x(idx_e) = 0.d0

ee = sum(krome_get_charges()*x(:)*imass(:))
ee = max(ee*e_mass,0d0)
x(idx_e) = ee

!2. mass fraction consistency
sumx = sum(x)

!NOTE: uncomment here if you want some additional control
!conservation error threshold: rise an error if above xerr
!xerr = 1d-2
!if(abs(sum-1d0)>xerr) then
!   print *,"ERROR: some problem with conservation!"
!   print *,"|sum(x)-1|=",abs(sum-1d0)
!   stop
!end if

isumx = 1d0/sumx
x(:) = x(:) * isumx

end subroutine krome_consistent_x

!*********************
!return an array sized krome_nmols containing
! the mass fractions (#), computed from the number
! densities (1/cm3) and the total density in g/cm3
function krome_n2x(n,rhogas)
use krome_commons
implicit none
real*8 :: n(nmols),krome_n2x(nmols)
real*8,value :: rhogas

krome_n2x(:) = n(:) * krome_get_mass() / rhogas

end function krome_n2x

!********************
!return an array sized krome_nmols containing
! the number densities (1/cm3), computed from the mass
! fractions and the total density in g/cm3
function krome_x2n(x,rhogas)
use krome_commons
implicit none
real*8 :: x(nmols),krome_x2n(nmols)
real*8,value :: rhogas

!compute densities from fractions
krome_x2n(:) = rhogas * x(:) * krome_get_imass()

end function krome_x2n

!******************
!returns free-fall time using the number density
! abundances of array x(:)
function krome_get_free_fall_time(x)
use krome_commons
use krome_getphys
implicit none
real*8::krome_get_free_fall_time
real*8::x(:),n(nspec)

n(1:nmols) = x(:)
n(nmols+1:nspec) = 0d0
krome_get_free_fall_time = get_free_fall_time(n(:))

end function krome_get_free_fall_time

!******************
!returns free-fall time using the total mass density
!  of gas, rhogas (g/cm3)
function krome_get_free_fall_time_rho(rhogas)
use krome_getphys
implicit none
real*8::krome_get_free_fall_time_rho
real*8::rhogas

krome_get_free_fall_time_rho = get_free_fall_time_rho(rhogas)

end function krome_get_free_fall_time_rho

!*******************
!do only cooling and heating
subroutine krome_thermo(x,Tgas,dt)
use krome_commons
use krome_cooling
use krome_heating
use krome_subs
use krome_tabs
use krome_constants
use krome_gadiab
implicit none
real*8 :: x(nmols)
real*8 :: Tgas,dt
real*8::n(nspec),nH2dust,dTgas,k(nrea),krome_gamma

nH2dust = 0d0
n(:) = 0d0
n(idx_Tgas) = Tgas
n(1:nmols) = x(:)
k(:) = coe_tab(n(:)) !compute coefficients
krome_gamma = gamma_index(n(:))

dTgas = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
    * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

Tgas = Tgas + dTgas*dt !update gas

end subroutine krome_thermo

!*************************
!get heating (erg/cm3/s) for a given species
! array x(:) and Tgas
function krome_get_heating(x,inTgas)
use krome_heating
use krome_subs
use krome_commons
implicit none
real*8 :: inTgas
real*8 :: x(nmols), krome_get_heating
real*8::Tgas,k(nrea),nH2dust,n(nspec)
n(1:nmols) = x(:)
Tgas = inTgas
n(idx_Tgas) = Tgas
k(:) = coe(n(:))
nH2dust = 0d0
krome_get_heating = heating(n(:),Tgas,k(:),nH2dust)
end function krome_get_heating

!*****************************
! get an array containing individual heatings (erg/cm3/s)
! the array has size krome_nheats. see heatcool.gps
! for index list
function krome_get_heating_array(x,inTgas)
use krome_heating
use krome_subs
use krome_commons
implicit none
real*8::n(nspec),Tgas,k(nrea),nH2dust
real*8 :: x(nmols),krome_get_heating_array(nheats)
real*8,value :: inTgas

n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = inTgas
!#KROME_Tdust_copy
k(:) = coe(n(:))
Tgas = inTgas
nH2dust = 0d0
krome_get_heating_array(:) = get_heating_array(n(:),Tgas,k(:),nH2dust)

end function krome_get_heating_array

!*************************
!get cooling (erg/cm3/s) for x(:) species array
! and Tgas
function krome_get_cooling(x,inTgas)
use krome_cooling
use krome_commons
implicit none
real*8 :: inTgas
real*8 :: x(nmols), krome_get_cooling
real*8::Tgas,n(nspec)
n(1:nmols) = x(:)
Tgas = inTgas
n(idx_Tgas) = Tgas
krome_get_cooling = cooling(n,Tgas)
end function krome_get_cooling

!*****************************
! get an array containing individual coolings (erg/cm3/s)
! the array has size krome_ncools. see heatcool.gps
! for index list
function krome_get_cooling_array(x,inTgas)
use krome_cooling
use krome_commons
implicit none
real*8::n(nspec),Tgas
real*8 :: x(nmols),krome_get_cooling_array(ncools)
real*8,value :: inTgas

n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = inTgas
!#KROME_Tdust_copy
Tgas = inTgas
krome_get_cooling_array(:) = get_cooling_array(n(:),Tgas)

end function krome_get_cooling_array

!******************
!alias of plot_cool
subroutine krome_plot_cooling(n)
use krome_cooling
implicit none
real*8 :: n(krome_nmols)

call plot_cool(n(:))

end subroutine krome_plot_cooling

!****************
!alias for dumping cooling in the unit nfile_in
subroutine krome_dump_cooling(n,Tgas,nfile_in)
use krome_cooling
use krome_commons
implicit none
real*8 :: n(nmols)
real*8 :: Tgas
real*8::x(nspec)
integer, optional :: nfile_in
integer::nfile
nfile = 31
x(:) = 0.d0
x(1:nmols) = n(:)
if(present(nfile_in)) nfile = nfile_in
call dump_cool(x(:),Tgas,nfile)

end subroutine krome_dump_cooling

!************************
!conserve the total amount of nucleii,
! alias for conserveLin_x in subs
subroutine krome_conserveLin_x(x,ref)
use krome_commons
use krome_subs
implicit none
real*8 :: x(nmols),ref(natoms)

call conserveLin_x(x(:),ref(:))

end subroutine krome_conserveLin_x

!************************
!conserve the total amount of nucleii,
! alias for conserveLin_x in subs
function krome_conserveLinGetRef_x(x)
use krome_commons
use krome_subs
implicit none
real*8 :: x(nmols),krome_conserveLinGetRef_x(natoms)

krome_conserveLinGetRef_x(:) = &
    conserveLinGetRef_x(x(:))

end function krome_conserveLinGetRef_x

!*************************
!force conservation to array x(:)
!using xi(:) as initial abundances.
!alias for conserve in krome_subs
function krome_conserve(x,xi)
use krome_subs
implicit none
real*8 :: x(krome_nmols),xi(krome_nmols),krome_conserve(krome_nmols)
real*8::n(krome_nspec),ni(krome_nspec)

n(:) = 0d0
ni(:) = 0d0
n(1:krome_nmols) = x(1:krome_nmols)
ni(1:krome_nmols) = xi(1:krome_nmols)
n(:) = conserve(n(:), ni(:))
krome_conserve(:) = n(1:krome_nmols)

end function krome_conserve

!***************************
!get the adiabatic index for x(:) species abundances
! and Tgas.
! alias for gamma_index in krome_subs
function krome_get_gamma(x,Tgas)
use krome_subs
use krome_commons
use krome_gadiab
real*8 :: Tgas
real*8 :: x(nmols), krome_get_gamma
real*8::n(nspec)
n(:) = 0.d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas
krome_get_gamma = gamma_index(n(:))
end function krome_get_gamma

!***************************
!get an integer array containing the atomic numbers Z
! of the spcecies.
! alias for get_zatoms
function krome_get_zatoms()
use krome_subs
use krome_commons
use krome_getphys
implicit none
integer :: krome_get_zatoms(nmols)
integer::zatoms(nspec)

zatoms(:) = get_zatoms()
krome_get_zatoms(:) = zatoms(1:nmols)

end function krome_get_zatoms

!****************************
!get the mean molecular weight from
! number density and mass density.
! alias for get_mu in krome_subs module
function krome_get_mu(x)
use krome_commons
use krome_subs
use krome_getphys
implicit none
real*8 :: x(nmols), krome_get_mu
real*8::n(1:nspec)
n(:) = 0d0
n(1:nmols) = x(:)
krome_get_mu = get_mu(n(:))
end function krome_get_mu

!***************************
!get the names of the reactions as a
! character*50 array of krome_nrea
! elements
!! !! cannot yet be called from C
function krome_get_rnames()
use krome_commons
use krome_subs
use krome_getphys
implicit none
character*50 :: krome_get_rnames(nrea)

krome_get_rnames(:) = get_rnames()

end function krome_get_rnames

!*****************
!get an array of double containing the masses in g
! of the species.
! alias for get_mass in krome_subs
function krome_get_mass()
use krome_subs
use krome_commons
use krome_getphys
implicit none
real*8::tmp(nspec)
real*8 :: krome_get_mass(nmols)
tmp(:) = get_mass()
krome_get_mass = tmp(1:nmols)
end function krome_get_mass

!*****************
!get an array of double containing the inverse
! of the mass (1/g) of the species
!alias for get_imass in krome_subs
function krome_get_imass()
use krome_subs
use krome_commons
use krome_getphys
implicit none
real*8::tmp(nspec)
real*8 :: krome_get_imass(nmols)
tmp(:) = get_imass()
krome_get_imass = tmp(1:nmols)
end function krome_get_imass

!***********************
!get the total number of H nuclei
function krome_get_Hnuclei(x)
use krome_commons
use krome_subs
use krome_getphys
real*8::n(nspec)
real*8 :: krome_get_Hnuclei, x(nmols)
n(:) = 0d0
n(1:nmols) = x(:)

krome_get_Hnuclei = get_Hnuclei(n(:))

end function krome_get_Hnuclei

!*****************
!get an array of size krome_nmols containing the
! charges of the species.
! alias for get_charges
function krome_get_charges()
use krome_subs
use krome_commons
use krome_getphys
implicit none
real*8::tmp(nspec)
real*8 :: krome_get_charges(nmols)
tmp(:) = get_charges()
krome_get_charges = tmp(1:nmols)
end function krome_get_charges

!*****************
!get an array of character*16 and size krome_nmols
! containing the names of all the species.
! alias for get_names
!!  !! cannot yet be called from C
function krome_get_names()
use krome_subs
use krome_commons
use krome_getphys
implicit none
character*16 :: krome_get_names(nmols)
character*16::tmp(nspec)
tmp(:) = get_names()
krome_get_names = tmp(1:nmols)
end function krome_get_names

!********************
!get space-separated header of chemical species
function krome_get_names_header()
use krome_commons
use krome_getphys
implicit none
character*141::krome_get_names_header
character*16::tmp(nspec)
integer::i

tmp(:) = get_names()

krome_get_names_header = ""
do i=1,nmols
krome_get_names_header = trim(krome_get_names_header)//" "//trim(tmp(i))
end do

end function krome_get_names_header

!********************
!get space-separated header of coolings
function krome_get_cooling_names_header()
use krome_commons
use krome_getphys
implicit none
character*130::krome_get_cooling_names_header
character*16::tmp(ncools)
integer::i

tmp(:) = get_cooling_names()

krome_get_cooling_names_header = ""
do i=1,ncools
if(trim(tmp(i))=="") cycle
krome_get_cooling_names_header = trim(krome_get_cooling_names_header)//" "//trim(tmp(i))
end do

end function krome_get_cooling_names_header

!********************
!get space-separated header of heatings
function krome_get_heating_names_header()
use krome_commons
use krome_getphys
implicit none
character*87::krome_get_heating_names_header
character*16::tmp(nheats)
integer::i

tmp(:) = get_heating_names()

krome_get_heating_names_header = ""
do i=1,nheats
if(trim(tmp(i))=="") cycle
krome_get_heating_names_header = trim(krome_get_heating_names_header)//" "//trim(tmp(i))
end do

end function krome_get_heating_names_header

!*****************
!get the index of the species with name name.
! alias for get_index
!! !! cannot yet be called from C
function krome_get_index(name)
use krome_subs
implicit none
integer :: krome_get_index
character*(*) :: name
krome_get_index = get_index(name)
end function krome_get_index

!*******************
!get the total density of the gas in g/cm3
! giving all the number densities n(:)
function krome_get_rho(n)
use krome_commons
real*8 :: krome_get_rho, n(nmols)
real*8::m(nmols)
m(:) = krome_get_mass()
krome_get_rho = sum(m(:)*n(:))
end function krome_get_rho

!*************************
!scale the abundances of the metals contained in n(:)
! to Z according to Asplund+2009.
! note that this applies only to neutral atoms.
subroutine krome_scale_Z(x,Z)
use krome_commons
use krome_getphys
real*8 :: x(nmols)
real*8 :: Z
real*8::Htot,n(nspec)

n(1:nmols) = x(:)
n(nmols+1:nspec) = 0d0

Htot = get_Hnuclei(n(:))
x(idx_C) = max(Htot * 1d1**(Z+(-3.57)), 1d-40)
x(idx_O) = max(Htot * 1d1**(Z+(-3.31)), 1d-40)

end subroutine krome_scale_Z

!*************************
!set the total metallicity
! in terms of Z/Z_solar
subroutine krome_set_Z(xarg)
use krome_commons
real*8 :: xarg

total_Z = xarg

end subroutine krome_set_Z

!*************************
!set D is in terms of D_solar (D/D_sol).
subroutine krome_set_dust_to_gas(xarg)
use krome_commons
real*8 :: xarg

dust2gas_ratio = xarg

end subroutine

!*************************
!set the clumping factor
subroutine krome_set_clump(xarg)
use krome_commons
real*8 :: xarg

clump_factor = xarg

end subroutine krome_set_clump

!***********************
!get the number of electrons assuming
! total neutral charge (cations-anions)
function krome_get_electrons(x)
use krome_commons
use krome_subs
use krome_getphys
real*8 :: x(nmols), krome_get_electrons
real*8::n(nspec)
n(1:nmols) = x(:)
n(nmols+1:nspec) = 0d0
krome_get_electrons = get_electrons(n(:))
end function krome_get_electrons

!**********************
!print on screen the first nbest highest reaction fluxes
subroutine krome_print_best_flux(xin,Tgas,nbest)
use krome_subs
use krome_commons
implicit none
real*8 :: xin(nmols)
real*8 :: Tgas
real*8::x(nmols),n(nspec)
integer :: nbest
n(1:nmols) = xin(:)
n(idx_Tgas) = Tgas
call print_best_flux(n,Tgas,nbest)

end subroutine krome_print_best_flux

!*********************
!print only the highest fluxes greater than a fraction frac
! of the maximum flux
subroutine krome_print_best_flux_frac(xin,Tgas,frac)
use krome_subs
use krome_commons
implicit none
real*8 :: xin(nmols)
real*8 :: Tgas,frac
real*8::n(nspec)
n(1:nmols) = xin(:)
n(idx_Tgas) = Tgas
call print_best_flux_frac(n,Tgas,frac)

end subroutine krome_print_best_flux_frac

!**********************
!print the highest nbest fluxes for reactions involving
!a given species using the index idx_find (e.g. krome_idx_H2)
subroutine krome_print_best_flux_spec(xin,Tgas,nbest,idx_find)
use krome_subs
use krome_commons
implicit none
real*8 :: xin(nmols)
real*8 :: Tgas
real*8::n(nspec)
integer :: nbest,idx_find
n(1:nmols) = xin(:)
n(idx_Tgas) = Tgas
call print_best_flux_spec(n,Tgas,nbest,idx_find)
end subroutine krome_print_best_flux_spec

!*******************************
!get an array of size krome_nrea with
! the fluxes of all the reactions in cm-3/s
function krome_get_flux(n,Tgas)
use krome_commons
use krome_subs
real*8 :: krome_get_flux(nrea),n(nmols)
real*8,value :: Tgas
real*8::x(nspec)
x(:) = 0.d0
x(1:nmols) = n(:)
x(idx_Tgas) = Tgas
krome_get_flux(:) = get_flux(x(:), Tgas)
end function krome_get_flux

!*****************************
!store the fluxes to the file unit ifile
! using the chemical composition x(:), and the
! gas temperature Tgas. xvar is th value of an
! user-defined independent variable that
! can be employed for plots.
! the file columns are as follow
! rate number, xvar, absolute flux,
!  flux/maxflux, flux fraction wrt total,
!  reaction name (*50 string)
subroutine krome_explore_flux(x,Tgas,ifile,xvar)
use krome_commons
use krome_subs
use krome_getphys
implicit none
real*8 :: x(nmols)
real*8 :: Tgas,xvar
real*8::flux(nrea),fluxmax,sumflux,n(nspec)
integer :: ifile
integer::i
character*50::rname(nrea)

!get reaction names
rname(:) = get_rnames()
n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas
!get fluxes
flux(:) = get_flux(n(:), Tgas)
fluxmax = maxval(flux) !maximum flux
sumflux = sum(flux) !sum of all the fluxes
!loop on reactions
do i=1,nrea
write(ifile,'(I8,5E17.8e3,a3,a50)') i,xvar,Tgas,flux(i),&
    flux(i)/fluxmax, flux(i)/sumflux," ",rname(i)
end do
write(ifile,*)

end subroutine krome_explore_flux

!*********************
!get nulcear qeff for the reactions
function krome_get_qeff()
use krome_commons
use krome_subs
use krome_getphys
implicit none
real*8 :: krome_get_qeff(nrea)

krome_get_qeff(:) = get_qeff()

end function krome_get_qeff

!************************
!dump the fluxes to the file unit nfile
subroutine krome_dump_flux(n,Tgas,nfile)
use krome_commons
real*8 :: n(nmols)
real*8 :: Tgas
real*8::flux(nrea)
integer :: nfile
integer::i

flux(:) = krome_get_flux(n(:),Tgas)
do i=1,nrea
write(nfile,'(I8,E17.8e3)') i,flux(i)
end do
write(nfile,*)

end subroutine krome_dump_flux

!************************
!dump all the evaluation of the coefficient rates in
! the file funit, in the range inTmin, inTmax, using
! imax points
subroutine krome_dump_rates(inTmin,inTmax,imax,funit)
use krome_commons
use krome_subs
implicit none
integer::i,j
integer :: funit,imax
real*8 :: inTmin,inTmax
real*8::Tmin,Tmax,Tgas,k(nrea),n(nspec)

Tmin = log10(inTmin)
Tmax = log10(inTmax)

n(:) = 1d-40
do i=1,imax
Tgas = 1d1**((i-1)*(Tmax-Tmin)/(imax-1)+Tmin)
n(idx_Tgas) = Tgas
k(:) = coe(n(:))
do j=1,nrea
write(funit,'(E17.8e3,I8,E17.8e3)') Tgas,j,k(j)
end do
write(funit,*)
end do

end subroutine krome_dump_rates

!************************
!print species informations on screen
subroutine krome_get_info(x, Tgas)
use krome_commons
use krome_subs
use krome_getphys
implicit none
integer::i,charges(nspec)
real*8 :: x(nmols)
real*8 :: Tgas
real*8::masses(nspec)
character*16::names(nspec)

names(:) = get_names()
charges(:) = get_charges()
masses(:) = get_mass()

print '(a4,a10,a11,a5,a11)',"#","Name","m (g)","Chrg","x"
do i=1,size(x)
print '(I4,a10,E11.3,I5,E11.3)',i," "//names(i),masses(i),charges(i),x(i)
end do
print '(a30,E11.3)'," sum",sum(x)

print '(a14,E11.3)',"Tgas",Tgas
end subroutine krome_get_info

!*****************************
subroutine krome_set_mpi_rank(xarg)
use krome_commons
implicit none
integer :: xarg
krome_mpi_rank=xarg
end subroutine krome_set_mpi_rank

end module krome_user

!############### MODULE ##############
module krome_reduction
contains

!**************************
function fex_check(n,Tgas)
use krome_commons
use krome_tabs
implicit none
integer::i
integer::r1,r2,r3
real*8::fex_check,n(nspec),k(nrea),rrmax,Tgas

k(:) = coe_tab(n(:))
rrmax = 0.d0
n(idx_dummy) = 1.d0
n(idx_g) = 1.d0
n(idx_CR) = 1.d0
do i=1,nrea
r1 = arr_r1(i)
r2 = arr_r2(i)
r3 = arr_r3(i)
arr_flux(i) = k(i)*n(r1)*n(r2)*n(r3)
rrmax = max(rrmax, arr_flux(i))
end do
fex_check = rrmax

end function fex_check

end module krome_reduction

!############### MODULE ##############
module krome_main

integer::krome_call_to_fex
!$omp threadprivate(krome_call_to_fex)

contains

! *************************************************************
!  This file has been generated with:
!  KROME 14.08.dev on 2018-10-24 12:49:10
!  Changeset b21f657
!  see http://kromepackage.org
!
!  Written and developed by Tommaso Grassi and Stefano Bovino
!
!  Contributors:
!  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
!  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
!  E.Tognelli
!  KROME is provided "as it is", without any warranty.
! *************************************************************

!********************************
!KROME main (interface to the solver library)

subroutine krome(x,Tgas,dt  )
use krome_commons
use krome_subs
use krome_ode
use krome_reduction
use krome_dust
use krome_getphys
use krome_tabs
implicit none
real*8 :: Tgas,dt
real*8 :: x(nmols)
real*8 :: rhogas

real*8::mass(nspec),n(nspec),tloc,xin
real*8::rrmax,totmass,n_old(nspec),ni(nspec),invTdust(ndust)
integer::icount,i,icount_max
integer:: ierr

!DLSODES variables
integer,parameter::meth=2 !1=adam, 2=BDF
integer::neq(1),itol,itask,istate,iopt,lrw,liw,mf
integer::iwork(641)
real*8::atol(nspec),rtol(nspec)
real*8::rwork(3957)
logical::got_error,equil

!****************************
!init DLSODES (see DLSODES manual)
call XSETF(0)!toggle solver verbosity
got_error = .false.
neq = nspec !number of eqns
liw = size(iwork)
lrw = size(rwork)
iwork(:) = 0
rwork(:) = 0d0
itol = 4 !both tolerances are scalar
rtol(:) = 1.000000d-04 !relative tolerance
atol(:) = 1.000000d-20 !absolute tolerance
icount_max = 100 !maximum number of iterations

itask = 1
iopt = 0

!MF=
!  = 222 internal-generated JAC and sparsity
!  = 121 user-provided JAC and internal generated sparsity
!  =  22 internal-generated JAC but sparsity user-provided
!  =  21 user-provided JAC and sparsity
MF = 222
!end init DLSODES
!****************************

ierr = 0 !error flag, zero==OK!
n(:) = 0d0 !initialize densities

n(1:nmols) = x(:)

n(idx_Tgas) = Tgas !put temperature in the input array

icount = 0 !count solver iterations
istate = 1 !init solver state
tloc = 0.d0 !set starting time

!store initial values
ni(:) = n(:)
n_global(:) = n(:)

call makeStoreOnceRates(n(:))

n_old(:) = -1d99
krome_call_to_fex = 0
do
icount = icount + 1
!solve ODE
CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, &
    ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
    LIW, JES, MF)

krome_call_to_fex = krome_call_to_fex + IWORK(12)
!check DLSODES exit status
if(istate==2) then
exit !sucsessful integration
elseif(istate==-1) then
istate = 1 !exceeded internal max iterations
elseif(istate==-5 .or. istate==-4) then
istate = 3 !wrong sparsity recompute
elseif(istate==-3) then
n(:) = ni(:)
istate = 1
else
got_error = .true.
end if

if(got_error.or.icount>icount_max) then
if (krome_mpi_rank>0) then
print *,krome_mpi_rank,"ERROR: wrong solver exit status!"
print *,krome_mpi_rank,"istate:",istate
print *,krome_mpi_rank,"iter count:",icount
print *,krome_mpi_rank,"max iter count:",icount_max
print *,krome_mpi_rank,"SEE KROME_ERROR_REPORT file"
else
print *,"ERROR: wrong solver exit status!"
print *,"istate:",istate
print *,"iter count:",icount
print *,"max iter count:",icount_max
print *,"SEE KROME_ERROR_REPORT file"
end if
call krome_dump(n(:), rwork(:), iwork(:), ni(:))
stop
end if

end do

!avoid negative species
do i=1,nspec
n(i) = max(n(i),0d0)
end do

!returns to user array
x(:) = n(1:nmols)

Tgas = n(idx_Tgas) !get new temperature

end subroutine krome

!*********************************
!integrates to equilibrium using constant temperature
subroutine krome_equilibrium(x,Tgas,verbosity)
use krome_ode
use krome_subs
use krome_commons
use krome_constants
use krome_getphys
use krome_tabs
implicit none
integer::mf,liw,lrw,itol,meth,iopt,itask,istate,neq(1)
integer::i,imax
integer,optional::verbosity
integer::verbose
real*8 :: Tgas
real*8 :: x(nmols)
real*8 :: rhogas
real*8::tloc,n(nspec),mass(nspec),ni(nspec)
real*8::dt,xin
integer::iwork(641)
real*8::atol(nspec),rtol(nspec)
real*8::rwork(3957)
real*8::ertol,eatol,max_time,t_tot,ntot_tol,err_species
logical::converged

integer, save :: ncall=0
integer, parameter :: ncall_print_frequency=20000
integer :: ncallp
integer::charges(nspec)
real*8::masses(nspec)
character*16::names(nspec)

!set verbosity from argument
verbose = 1 !default is verbose
if(present(verbosity)) verbose = verbosity

PLW = 2.11814e-13
PHI = 1.08928e-13
PHEI = 2.76947e-14
PCVI = 1.03070e-17

call XSETF(0)!toggle solver verbosity
meth = 2
neq = nspec !number of eqns
liw = size(iwork)
lrw = size(rwork)
iwork(:) = 0
rwork(:) = 0d0
itol = 4 !both tolerances are scalar
rtol(:) = 1d-6 !relative tolerance
atol(:) = 1d-20 !absolute tolerance

! Switches to decide when equilibrium has been reached
ertol = 1d-5  ! relative min change in a species
eatol = 1d-12 ! absolute min change in a species
max_time=seconds_per_year*5d8 ! max time we will be integrating for

!for DLSODES options see its manual
iopt = 0
itask = 1
istate = 1

mf = 222 !internally evaluated sparsity and jacobian
tloc = 0d0 !initial time

n(:) = 0d0 !initialize densities
!copy into array
n(nmols+1:) = 0d0
n(1:nmols) = x(:)

n(idx_Tgas) = Tgas

!store previous values
ni(:) = n(:)
n_global(:) = ni(:)

call makeStoreOnceRates(n(:))

imax = 1000

dt = seconds_per_year * 1d2
t_tot = dt
converged = .false.
do while (.not. converged)
do i=1,imax
!solve ODE
CALL DLSODES(fcn_tconst, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
    ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jcn_dummy, MF)
if(istate==2) then
exit
else
istate=1
end if
end do
!check errors
if(istate.ne.2) then
print *,"ERROR: no equilibrium found!"
stop
end if

!avoid negative species
do i=1,nspec
n(i) = max(n(i),0d0)
end do

! check if we have converged by comparing the error in any species with an relative abundance above eatol
converged = maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) .lt. ertol &
    .or. t_tot .gt. max_time

! Increase integration time by a reasonable factor
if(.not. converged) then
dt = dt * 3.
t_tot = t_tot + dt
ni = n
n_global = n
endif
enddo
!returns to user array
x(:) = n(1:nmols)

if(t_tot > max_time .and. &
    maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) > 0.2 .and. verbose>0) then
print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
print *, 'Tgas :', Tgas
names(:) = get_names()
charges(:) = get_charges()
masses(:) = get_mass()

print '(a4,a10,a11,a5,a16)',"#","Name","m (g)","Chrg","  Current / Last"
do i=1,nmols
print '(I4,a10,E11.3,I5,2E14.6,E11.3)',i," "//names(i),masses(i),charges(i),n(i),ni(i),abs(n(i) - ni(i)) / max(n(i),eatol*sum(n(1:nmols)))
end do
print '(a30,2E14.6)'," sum",sum(n(1:nmols)),sum(ni(1:nmols))
print *, 'Fractional error :', maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols))))
print *, 'Absolute and relative floors:', eatol, ertol
end if

! Print info ever so often
!$omp critical
ncall=ncall+1
ncallp = ncall
!$omp end critical

if(modulo(ncallp,ncall_print_frequency)==0 .and. verbose>0) then
print *, 'Found equilibrium for ', ncallp, ' cells.'
end if

end subroutine krome_equilibrium

!********************
!dummy jacobian
subroutine jcn_dummy()
implicit none
end subroutine jcn_dummy

!*******************
!dn/dt where dT/dt=0
subroutine fcn_tconst(n,tt,x,f)
use krome_commons
use krome_ode
implicit none
integer::n,ierr
real*8::x(n),f(n),tt
call fex(n,tt,x(:),f(:))
f(idx_Tgas) = 0d0
end subroutine fcn_tconst

!*******************************
subroutine krome_dump(n,rwork,iwork,ni)
use krome_commons
use krome_subs
use krome_tabs
use krome_reduction
use krome_ode
use krome_getphys
integer::fnum,i,iwork(:),idx(nrea),j
real*8::n(:),rwork(:),rrmax,k(nrea),kmax,rperc,kperc,dn(nspec),tt,ni(:)
character*16::names(nspec),FMTi,FMTr
character*50::rnames(nrea),fname,prex
integer,save::mx_dump=1000 ! max nr of reports before terminating
fnum = 99
if (krome_mpi_rank>0) then
write(fname,'(a,i5.5)') "KROME_ERROR_REPORT_",krome_mpi_rank
else
fname = "KROME_ERROR_REPORT"
endif
open(fnum,FILE=trim(fname),status="replace")
tt = 0d0
names(:) = get_names()
rnames(:) = get_rnames()
call fex(nspec,tt,n(:),dn(:))

write(fnum,*) "KROME ERROR REPORT"
write(fnum,*)
!SPECIES
write(fnum,*) "Species abundances"
write(fnum,*) "**********************"
write(fnum,'(a5,a20,3a12)') "#","name","qty","dn/dt","ninit"
write(fnum,*) "**********************"
do i=1,nspec
write(fnum,'(I5,a20,3E12.3e3)') i,names(i),n(i),dn(i),ni(i)
end do
write(fnum,*) "**********************"

!F90 FRIENDLY RESTART
write(fnum,*)
write(fnum,*) "**********************"
write(fnum,*) "F90-friendly species"
write(fnum,*) "**********************"
do i=1,nspec
write(prex,'(a,i3,a)') "x(",i,") = "
write(fnum,*) trim(prex),ni(i),"!"//names(i)
end do

write(fnum,*) "**********************"

!RATE COEFFIECIENTS
k(:) = coe_tab(n(:))
idx(:) = idx_sort(k(:))
kmax = maxval(k)
write(fnum,*)
write(fnum,*) "Rate coefficients (sorted) at Tgas",n(idx_Tgas)
write(fnum,*) "**********************"
write(fnum,'(a5,2a12,a10)') "#","k","k %","  name"
write(fnum,*) "**********************"
do j=1,nrea
i = idx(j)
kperc = 0.d0
if(kmax>0.d0) kperc = k(i)*1d2/kmax
write(fnum,'(I5,2E12.3e3,a2,a50)') i,k(i),kperc,"  ", rnames(i)
end do
write(fnum,*) "**********************"
write(fnum,*)

!FLUXES
call load_arrays
rrmax = fex_check(n(:), n(idx_Tgas))
idx(:) = idx_sort(arr_flux(:))
write(fnum,*)
write(fnum,*) "Reaction magnitude (sorted) [k*n1*n2*n3*...]"
write(fnum,*) "**********************"
write(fnum,'(a5,2a12,a10)') "#","flux","flux %","  name"
write(fnum,*) "**********************"
do j=1,nrea
i = idx(j)
rperc = 0.d0
if(rrmax>0.d0) rperc = arr_flux(i)*1d2/rrmax
write(fnum,'(I5,2E12.3e3,a2,a50)') i,arr_flux(i),rperc,"  ",rnames(i)
end do
write(fnum,*) "**********************"
write(fnum,*)

!SOLVER
FMTr = "(a30,E16.7e3)"
FMTi = "(a30,I10)"
write(fnum,*) "Solver-related information:"
write(fnum,FMTr) "step size last",rwork(11)
write(fnum,FMTr) "step size attempt",rwork(12)
write(fnum,FMTr) "time current",rwork(13)
write(fnum,FMTr) "tol scale factor",rwork(14)
write(fnum,FMTi) "numeber of steps",iwork(11)
write(fnum,FMTi) "call to fex",iwork(12)
write(fnum,FMTi) "call to jex",iwork(13)
write(fnum,FMTi) "last order used",iwork(14)
write(fnum,FMTi) "order attempt",iwork(15)
write(fnum,FMTi) "idx largest error",iwork(16)
write(fnum,FMTi) "RWORK size required",iwork(17)
write(fnum,FMTi) "IWORK size required",iwork(18)
write(fnum,FMTi) "NNZ in Jac",iwork(19)
write(fnum,FMTi) "extra fex to compute jac",iwork(20)
write(fnum,FMTi) "number of LU decomp",iwork(21)
write(fnum,FMTi) "base address in RWORK",iwork(22)
write(fnum,FMTi) "base address of IAN",iwork(23)
write(fnum,FMTi) "base address of JAN",iwork(24)
write(fnum,FMTi) "NNZ in lower LU",iwork(25)
write(fnum,FMTi) "NNZ in upper LU",iwork(21)
write(fnum,*) "See DLSODES manual for further details on Optional Outputs"
write(fnum,*)
write(fnum,*) "END KROME ERROR REPORT"
write(fnum,*)
close(fnum)

mx_dump = mx_dump - 1
if (mx_dump==0) stop

end subroutine krome_dump

!********************************
subroutine krome_init()
use krome_commons
use krome_tabs
use krome_subs
use krome_reduction
use krome_dust
use krome_cooling
use krome_photo
use krome_fit

!init phys common variables
!$omp parallel
phys_Tcmb = 2.73d0
phys_zredshift = 0d0
phys_orthoParaRatio = 3d0
phys_metallicity = 0d0
phys_Tfloor = 2.73d0
!$omp end parallel

!init metallicity default
!assuming solar
total_Z = 1d0

!default D/D_sol = Z/Z_sol
!assuming linear scaling
dust2gas_ratio = total_Z

!default broadening turubulence velocity
broadeningVturb2 = 0d0

!default clumping factor for
! H2 formation on dust by Jura/Gnedin
clump_factor = 1d0

!default for thermo toggle is ON
!$omp parallel
krome_thermo_toggle = 1
!$omp end parallel

!load arrays with ractants/products indexes
call load_arrays()

!initialize cooling tabel for metals
call coolingZ_init_tabs()

call init_dust_tabs()

!initialize CO cooling
call init_coolingCO()

!initialize the table for exp(-a/T) function
call init_exp_table()

call load_parts()

!init photo reactants indexes
photoPartners(1) = idx_H
photoPartners(2) = idx_HE
photoPartners(3) = idx_HEj
photoPartners(4) = idx_O
photoPartners(5) = idx_C
photoPartners(6) = idx_H2
photoPartners(7) = idx_Hk
photoPartners(8) = idx_CH
photoPartners(9) = idx_CH
photoPartners(10) = idx_C2
photoPartners(11) = idx_OH
photoPartners(12) = idx_OH
photoPartners(13) = idx_H2O
photoPartners(14) = idx_H2O
photoPartners(15) = idx_O2
photoPartners(16) = idx_O2
photoPartners(17) = idx_H2

!get machine precision
krome_epsilon = epsilon(0d0)

!load verbatim reactions
call loadReactionsVerbatim()

end subroutine krome_init

!****************************
function krome_get_coe(x,Tgas)
!krome_get_coe: public interface to obtain rate coefficients
use krome_commons
use krome_subs
use krome_tabs
implicit none
real*8 :: krome_get_coe(nrea), x(nmols), Tgas
real*8::n(nspec)

n(:) = 0d0
n(1:nmols) = x(:)
n(idx_Tgas) = Tgas
krome_get_coe(:) = coe_tab(n(:))

end function krome_get_coe

!****************************
function krome_get_coeT(Tgas)
!krome_get_coeT: public interface to obtain rate coefficients
! with argument Tgas only
use krome_commons
use krome_subs
use krome_tabs
implicit none
real*8 :: krome_get_coeT(nrea),Tgas
real*8::n(nspec)
n(idx_Tgas) = Tgas
krome_get_coeT(:) = coe_tab(n(:))
end function krome_get_coeT

end module krome_main
