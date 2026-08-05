!-----------------------------------------------------------------------------*
!= This is Fast-TUV FTUV4.2
!=  with 11 bins starting from 205nm
!=
!=  In general the error between TUV and FTUV is within 5%, and in
!=     some cases, it is above 5%
!=  There are some limitations to use this FTUV
!=   (1) J for stratospheric species (28-57) is not applied
!=   (2) CH3CHO have only one channel CH3CHO + hv --> CH3 + HCO
!=   (3) O2, N2O, and HO2 have large erros (due to absorption in short W)
!=
!=    Tropospheric Ultraviolet-Visible (TUV) radiation model                   =*
!=    Version 4.2                                                           =*
!=    May 2003                                                             =*
!-----------------------------------------------------------------------------*
!= Developed by Sasha Madronich with important contributions from:           =*
!= Chris Fischer, Siri Flocke, Julia Lee-Taylor, Bernhard Meyer,           =*
!= Irina Petropavlovskikh,  Xuexi Tie, and Jun Zen.                           =*
!= Special thanks to Knut Stamnes and co-workers for the development of the  =*
!= Discrete Ordinates code, and to Warren Wiscombe and co-workers for the    =*
!= development of the solar zenith angle subroutine. Citations for the many  =*
!= data bases (e.g. extraterrestrial irradiances, molecular spectra) may be  =*
!= found in the data files headers and/or in the subroutines that read them. =*
!=              To contact the author, write to:                             =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu  or tuv@acd.ucar.edu                     =*
!-----------------------------------------------------------------------------*
!= This program is free software; you can redistribute it and/or modify      =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                   =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                   =*
!= To obtain a copy of the GNU General Public License, write to:           =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96,97,98,99,2000,01,02,03  University Corporation   =*
!= for Atmospheric Research                                                   =*
!-----------------------------------------------------------------------------*
!

module ModTuv
    use mpas_log
    use, intrinsic :: iso_fortran_env, only: real64, real32, int32, int64

    implicit none

    private

!Public Variables
    public :: kw, kt, ks, kj, sw, ns
    public  :: initialized

!Public Subroutines
    public :: InitTuv
    public :: Tuv
    public :: slabel
    public :: jlabel
    public :: nj
    public :: tfiles
    public :: files
    public :: f
    public :: xRef
    public :: reverse
    public :: wl, wu, wc, nw
    public :: wbioStart, wbioEnd
    public :: nrad
    public :: get_nz

    logical, parameter :: isMpi = .true.

! BROADLY USED PARAMETERS:
!_________________________________________________
! i/o file unit numbers
    integer, parameter :: kout = 53
    integer, parameter :: kin = 12
!_________________________________________________
! altitude, wavelength, time (or solar zenith angle) grids
    integer, parameter :: kw = 20 ! wavelength  --- test Luiz Flavio
!INTEGER, PARAMETER :: kw=650 ! wavelength
    integer, parameter :: kt = 100 ! time/sza
!_________________________________________________
! number of weighting functions
    integer, parameter :: ks = 60 !  wavelength dependent
    integer, parameter :: kj = 80 !  wavelength and altitude dependent
    real, parameter :: deltax = 1.0e-4 ! delta for adding points at
! beginning or end of data grids
    real, parameter :: radius = 6.371e+3 ! radius of the earth:
    real, parameter :: largest = 1.0e+36 ! largest number of the machine:
    real, parameter :: pzero = +10./largest ! small numbers (positive and negative)
    real, parameter :: nzero = -10./largest
    real, parameter :: precis = 1.e-7  ! machine precision
! More physical constants:
!_________________________________________________________________
! Na = 6.022142E23  mol-1        = Avogadro constant
! kb = 1.38065E-23  J K-1        = Boltzmann constant
! R  = 8.31447      J mol-1 K-1 = molar gas constant
! h  = 6.626068E-34 J s        = Planck constant
! c  = 2.99792458E8 m s-1        = speed of light in vacuum
! G  = 6.673E-11    m3 kg-1 s-2 = Netwonian constant of gravitation
! sb = 5.67040E-8   W m-2 K-4   = Stefan-Boltzmann constant
!_________________________________________________________________
! (1) From NIST Reference on Constants, Units, and Uncertainty
! http://physics.nist.gov/cuu/index.html Oct. 2001.
! (2) These constants are not assigned to variable names;  in other
! words this is not Fortran code, but only a text table for quick
! reference.  To use, you must declare a variable name/type and
! assign the value to that variable. Or assign as parameter (see
! example for pi above).

    integer, parameter :: totReact = 115
    integer, parameter :: nmug = 10
    integer, parameter :: maxstr = 100
    integer, parameter :: maxtrm = 100
    integer, parameter :: maxsqt = 1000
    real, dimension(17), parameter :: xslod = (/ &
                                      6.2180730e-21, 5.8473627e-22, 5.6996334e-22, &
                                      4.5627094e-22, 1.7668250e-22, 1.1178808e-22, &
                                      1.2040544e-22, 4.0994668e-23, 1.8450616e-23, &
                                      1.5639540e-23, 8.7961075e-24, 7.6475608e-24, &
                                      7.6260556e-24, 7.5565696e-24, 7.6334338e-24, &
                                      7.4371992e-24, 7.3642966e-24/)
!mz = 5 zenith angle 0,20,40,60,80
!ms = 73 species
!mp = 5 pol coeff  0,1,2,3,4
    integer, parameter :: mz = 5
    integer, parameter :: ms = 73
    integer, parameter :: mp = 5
    real, parameter :: lbar = 206.214
    character(LEN=50), parameter :: jlabel(kj) = &! Photolysis coefficients labels
                                    (/ &
                                    'O2 -> O + O                                       ', & !1
                                    'O3 -> O2 + O(1D)                                  ', & !2 R2 CB07
                                    'O3 -> O2 + O(3P)                                  ', & !3 R3 CB07
                                    'NO2 -> NO + O(3P)                                 ', & !4 R1 CB07
                                    'NO3 -> NO + O2                                    ', & !5 R7 CB07
                                    'NO3 -> NO2 + O(3P)                                ', & !6 R8 CB07
                                    'N2O5 -> NO3 + NO + O(3P)                          ', & !7
                                    'N2O5 -> NO3 + NO2                                 ', & !8
                                    'N2O -> N2 + O(1D)                                 ', & !9
                                    'HO2 -> OH + O                                     ', & !10
                                    'H2O2 -> 2 OH                                      ', & !11 R9 CB07
                                    'HNO2 -> OH + NO                                   ', & !12 R4 CB07
                                    'HNO3 -> OH + NO2                                  ', & !13 R5 CB07
                                    'HNO4 -> HO2 + NO2                                 ', & !14
                                    'CH2O -> H + HCO                                   ', & !15
                                    'CH2O -> H2 + CO                                   ', & !16 R10 CB07
                                    'CH3CHO -> CH3 + HCO                               ', & !17
                                    'CH3CHO -> CH4 + CO                                ', & !18
                                    'CH3CHO -> CH3CO + H                               ', & !19
                                    'C2H5CHO -> C2H5 + HCO                             ', &
                                    'CHOCHO -> HCO + HCO                               ', &
                                    'CHOCHO -> CH2O + CO                               ', &
                                    'CH3COCHO -> CH3CO + HCO                           ', &
                                    'CH3COCH3 -> CH3CO + CH3                           ', &
                                    'CH3OOH -> CH3O + OH                               ', &
                                    'CH3ONO2 -> CH3O + NO2                             ', &
                                    'CH3CO(OONO2) -> Products                          ', &
                                    'ClOO -> Products                                  ', &
                                    'ClONO2 -> Cl + NO3                                ', &
                                    'ClONO2 -> ClO + NO2                               ', &
                                    'CH3Cl -> Products                                 ', &
                                    'CCl2O -> Products                                 ', &
                                    'CCl4 -> Products                                  ', &
                                    'CClFO -> Products                                 ', &
                                    'CF2O -> Products                                  ', &
                                    'CF2ClCFCl2 (CFC-113) -> Products                  ', &
                                    'CF2ClCF2Cl (CFC-114) -> Products                  ', &
                                    'CF3CF2Cl (CFC-115) -> Products                    ', &
                                    'CCl3F (CFC-11) -> Products                        ', &
                                    'CCl2F2 (CFC-12) -> Products                       ', &
                                    'CH3CCl3 -> Products                               ', &
                                    'CF3CHCl2 (HCFC-123) -> Products                   ', &
                                    'CF3CHFCl (HCFC-124) -> Products                   ', &
                                    'CH3CFCl2 (HCFC-141b) -> Products                  ', &
                                    'CH3CF2Cl (HCFC-142b) -> Products                  ', &
                                    'CF3CF2CHCl2 (HCFC-225ca) -> Products              ', &
                                    'CF2ClCF2CHFCl (HCFC-225cb) -> Products            ', &
                                    'CHClF2 (HCFC-22) -> Products                      ', &
                                    'BrONO2 -> BrO + NO2                               ', &
                                    'BrONO2 -> Br + NO3                                ', &
                                    'CH3Br -> Products                                 ', &
                                    'CHBr3 -> Products                                 ', &
                                    'CF3Br (Halon-1301) -> Products                    ', &
                                    'CF2BrCF2Br (Halon-2402) -> Products               ', &
                                    'CF2Br2 (Halon-1202) -> Products                   ', &
                                    'CF2BrCl (Halon-1211) -> Products                  ', &
                                    'Cl2 -> Cl + Cl                                    ', &
                                    'CH2(OH)CHO -> Products                            ', &
                                    'CH3COCOCH3 -> Products                            ', &
                                    'CH3COCHCH2 -> Products                            ', &
                                    'CH2C(CH3)CHO -> Products                          ', &
                                    'CH3COCO(OH) -> Products                           ', &
                                    'CH3CH2ONO2 -> CH3CH2O + NO2                       ', &
                                    'CH3CHONO2CH3 -> CH3CHOCH3 + NO2                   ', &
                                    'CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2          ', &
                                    'CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2              ', &
                                    'C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2                ', &
                                    'ClOOCl -> Cl + ClOO                               ', &
                                    'CH2(OH)COCH3 -> CH3CO + CH2(OH)                   ', &
                                    'CH2(OH)COCH3 -> CH2(OH)CO + CH3                   ', &
                                    'HOBr -> OH + Br                                   ', &
                                    'BrO -> Br + O                                     ', &
                                    'Br2 -> Br + Br                                    ', &
                                    '                                                  ', &
                                    '                                                  ', &
                                    '                                                  ', &
                                    '                                                  ', &
                                    '                                                  ', &
                                    '                                                  ', &
                                    '                                                  ' &
                                    /)

!_______________________________________________________________________
! select desired extra-terrestrial solar irradiance, using msun:
!  1 =   extsol.flx:  De Luisi, JGR 80, 345-354, 1975
!                     280-400 nm, 1 nm steps.
!  2 =   lowsun3.flx:  Lowtran (John Bahr, priv. comm.)
!                      173.974-500000 nm, ca. 0.1 nm steps in UV-B
!  3 =   modtran1.flx:  Modtran (Gail Anderson, priv. comm.)
!                       200.55-949.40, 0.05 nm steps
!  4 =   nicolarv.flx:  wvl<300 nm from Nicolet, Plan. Sp. Sci., 29,  951-974, 1981.
!                       wvl>300 nm supplied by Thekaekera, Arvesen Applied Optics 8,
!                       11, 2215-2232, 1969 (also see Thekaekera, Applied Optics, 13,
!                       3, 518, 1974) but with corrections recommended by:
!                       Nicolet, Plan. Sp. Sci., 37, 1249-1289, 1989.
!                       270.0-299.0 nm in 0.5 nm steps
!                       299.6-340.0 nm in ca. 0.4 nm steps
!                       340.0-380.0 nm in ca. 0.2 nm steps
!                       380.0-470.0 nm in ca. 0.1 nm steps
!  5 =  solstice.flx:  From:   MX%"ROTTMAN@virgo.hao.ucar.edu" 12-OCT-1994 13:03:01.62
!                      Original data gave Wavelength in vacuum
!                      (Converted to wavelength in air using Pendorf, 1967, J. Opt. Soc. Am.)
!                      279.5 to 420 nm, 0.24 nm spectral resolution, approx 0.07 nm steps
!  6 =  suntoms.flx: (from TOMS CD-ROM).  280-340 nm, 0.05 nm steps.
!  7 =  neckel.flx:  H.Neckel and D.Labs, "The Solar Radiation Between 3300 and 12500 A",
!                    Solar Physics v.90, pp.205-258 (1984).
!                    1 nm between 330.5 and 529.5 nm
!                    2 nm between 631.0 and 709.0 nm
!                    5 nm between 872.5 and 1247.4 nm
!                    Units: must convert to W m-2 nm-1 from photons cm-2 s-1 nm-1
!  8 =  atlas3.flx:  ATLAS3-SUSIM 13 Nov 94 high resolution (0.15 nm FWHM)
!                    available by ftp from susim.nrl.navy.mil
!                    atlas3_1994_317_a.dat, downloaded 30 Sept 98.
!                    150-407.95 nm, in 0.05 nm steps
!                    (old version from Dianne Prinz through Jim Slusser)
!                    orig wavelengths in vac, correct here to air.
!  9 =  solstice.flx:  solstice 1991-1996, average
!                    119.5-420.5 nm in 1 nm steps

! 10 =  susim_hi.flx:  SUSIM SL2 high resolution
!                      120.5-400.0 in 0.05 nm intervals (0.15 nm resolution)
! 11 =  wmo85.flx: from WMO 1995 Ozone Atmospheric Ozone (report no. 16)
!                  on variable-size bins.  Original values are per bin, not
!                  per nm.
! 12 = combine susim_hi.flx for .lt. 350 nm, neckel.flx for .gt. 350 nm.

! 13 = combine
!     for wl(iw) .lt. 150.01                                susim_hi.flx
!     for wl(iw) .ge. 150.01 and wl(iw) .le. 400            atlas3.flx
!     for wl(iw) .gt. 400                                   Neckel & Labs
    integer, parameter :: msun = 13
    integer, parameter :: nTotalFiles = 145
! quantum yield recommendation:
!    kjpl87:  JPL recommendation 1987                - JPL 87, 90, 92 do not "tail"
!    kjpl92:  JPL recommendations 1990/92 (identical) - still with no "tail"
!    kjpl97:  JPL recommendation 1997, includes tail, similar to Shetter et al.
!    kmich :  Michelsen et al., 1994
!    kshet :  Shetter et al., 1996
!    kjpl00:  JPL 2000
!    kmats:  Matsumi et al., 2002
    integer, parameter :: kmich = 1
    integer, parameter :: kjpl87 = 2
    integer, parameter :: kjpl92 = 3
    integer, parameter :: kshet = 4
    integer, parameter :: kjpl97 = 5
    integer, parameter :: kjpl00 = 6
    integer, parameter :: kmats = 7
    integer, dimension(20), parameter :: &
        mOption = (/1, 7, 1, 2, 6, 1, 3, 1, 1, 1, 4, 1, 2, 5, 1, 4, 2, 2, 2, 0/)
    integer, parameter, dimension(5) :: indSza = (/1, 2, 3, 4, 5/)
    real, parameter, dimension(5) :: angle = (/0.0, 20.0, 40.0, 60.0, 80.0/)
    real, parameter, dimension(5) :: fatSum = (/-0.2, -0.2, 0.0, 0.0, 0.0/)
    real, parameter, dimension(5) :: ca0 = &
                                     (/4.52372, 4.52372, 4.99378, 0.969867, 1.07801/)
    real, parameter, dimension(5) :: ca1 = &
                                     (/-5.94317, -5.94317, -7.92752, -0.841035, -2.39580/)
    real, parameter, dimension(5) :: ca2 = &
                                     (/2.63156, 2.63156, 3.94715, 0.878835, 2.32632/)
    real, parameter, dimension(5) :: cb0 = &
                                     (/2.43360, 2.43360, 3.98265, 3.49843, 3.06312/)
    real, parameter, dimension(5) :: cb1 = &
                                     (/-3.61363, -3.61363, -6.90516, -5.98839, -5.26281/)
    real, parameter, dimension(5) :: cb2 = &
                                     (/2.19018, 2.19018, 3.93602, 3.50262, 3.20980/)

!LFR - TUV problem with ozone column
! In original code the column of Ozone was obtained from O3 climatology
! over USA.
! Thefore for the new column (dynamic) from the CCATT model we have
! some negative values for Photolisys rate. The cause is the ammount
! of ozone at each level. The adjO3 is an attenuation factor for each
! reaction
    real, parameter :: adjO3(kj) = (/ &
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & !00
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, & !01
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & !02
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & !03
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & !04
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & !05
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, & !06
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 & !07
                       /)
!
    integer, allocatable, dimension(:) :: nrad
    real :: s226(kw), s263(kw), s298(kw)
    real :: s218(kw), s228(kw), s243(kw), s295(kw)
    real :: c0(kw), c1(kw), c2(kw)
    real :: yg(kw, totReact), yg1(kw, totReact), yg1n(kw, totReact), yg2(kw, totReact)
    real :: yg3(kw, totReact), yg4(kw, totReact), yg5(kw, totReact)
    real :: coeff(4, 3, totReact)
    real :: c(maxtrm)
    real :: hugeVar, powmax, powmin, tinyVar
    integer :: maxmsg, nummsg
    logical :: msglim
    integer :: maxmsg2, nummsg2

    real :: tbar(totReact)

!-------------------------------------------------------
! XS_COEFFS containing the Chebyshev
! polynomial coeffs necessary to calculate O2 effective
! cross-sections
!-------------------------------------------------------
    real(real64) :: ac(20, 17)
    real(real64) :: bc(20, 17) ! Chebyshev polynomial coeffs
    real(real32) :: wave_num(17)

    logical :: initialized
    logical :: isread
    logical :: firstcall
    logical :: call1
    logical :: pass1
    logical :: pass2
    logical :: pass3
    logical :: pass4
    logical :: pass5
    logical :: pass6
    logical :: doinit
    logical :: doinit2

    integer :: ila
    integer :: isrb

    real :: dither

    real :: sqt(maxsqt)

    real ::  pi
    real ::  twopi
    real ::  rpd

    double precision :: tol
    double precision :: dmach(4)
    real :: rmach(4)
    real :: gmu(nmug), gwt(nmug), ylmg(0:maxstr, nmug)

    real :: coef(mz, ms, mp)  ! pol coefficent for FTUV

! Wavelength grid:
    integer :: nw, ns, nj, iw
    real :: wl(kw), wc(kw), wu(kw)
    real :: f(kw) !extra terrestrial solar flux

! O2 absorption cross section
    real :: o2xs1(kw)

! SO2 absorption cross section
    real :: so2xs(kw)

! NO2 absorption cross section
    real :: no2xs(kw)

    real :: mm_o3xs(kw)

    integer :: wbioStart, wbioEnd
    real :: sw(ks, kw)
    character(LEN=50) :: slabel(ks) ! Spectral weighting functions labels

    integer           :: xRef(kj)
    logical           :: doReaction(kj)
    integer :: ii, jj !Temporary
! Parameters for shifting wavelengths air <--> vacuum
    integer           :: mrefr
    logical           :: lrefr
    real              :: airout

!Files
    type tFiles
        character(LEN=200) :: fileName
    end type tFiles
    type(tFiles), allocatable, dimension(:) :: files

    real, allocatable, dimension(:) :: x, y
    real, allocatable, dimension(:) :: x1, y1
    real, allocatable, dimension(:) :: x2, y2
    real, allocatable, dimension(:) :: x3, y3
    real, allocatable, dimension(:) :: x4, y4
    real, allocatable, dimension(:) :: x5, y5

    integer :: ierr

    integer :: noPr
    integer :: current_nz   ! número de níveis de altitude da chamada atual
!LFR   REAL :: rgasog,deltap

    integer,parameter :: nstr=0
    real,parameter    :: wstart=120.000
    real,parameter    :: wstop=735.000
    integer,parameter :: nwint=-11

    !real, allocatable :: jphoto(:,:,:)

    contains

    subroutine InitTuv(filesName, myNum, chemical_mechanism, nCells, nVertlevels, nr_photo)
        use mpi
        !include 'mpif.h'

        integer, intent(IN) :: myNum, nCells, nVertlevels, nr_photo
        character(LEN=*), intent(IN) :: filesName
        character(LEN=200) :: filesHome
        character(LEN=2) :: fname
        character(LEN=*), intent(IN) :: chemical_mechanism

        integer :: i, nFiles, nOfFile, iang, js, jp

        !if (.not. allocated(jphoto)) then
        !    allocate(jphoto(nVertlevels, nCells, nr_photo))
        !end if

        if (myNum == 0) then !Initialization just for myNum=0

            if (initialized) then
                print *, 'ERROR: TUV already initialized!'
                print *, '### Please, check your code ###'
                !call flush (6)
                stop
            end if
            filesHome = ''
            if (len(trim(filesName)) == 0) then
                !If is to use de deafault files name (test case)
                allocate (files(nTotalFiles))
                !Default values of files names
                files(1)%fileName = 'input/POL.out'
                files(2)%fileName = 'datae1/sun/extsol.flx'
                files(3)%fileName = 'datae1/sun/lowsun3.flx'
                files(4)%fileName = 'datae1/sun/modtran1.flx'
                files(5)%fileName = 'datae1/sun/nicolarv.flx'
                files(6)%fileName = 'datae2/sun/solstice.flx'
                files(7)%fileName = 'datae2/sun/suntoms.flx'
                files(8)%fileName = 'datae1/sun/neckel.flx'
                files(9)%fileName = 'datae1/sun/atlas3_1994_317_a.dat'
                files(10)%fileName = files(6)%fileName
                files(11)%fileName = 'datae1/sun/susim_hi.flx'
                files(12)%fileName = 'datae1/sun/wmo85.flx'
                files(13)%fileName = files(8)%fileName
                files(14)%fileName = files(9)%fileName
                files(15)%fileName = files(8)%fileName
                files(16)%fileName = 'datae1/o2/O2_brasseur.abs'
                files(17)%fileName = 'datae1/o2/O2_yoshino.abs'
                files(18)%fileName = 'datae1/so2/SO2xs.all'
                files(19)%fileName = 'datae1/no2/NO2_ncar_00.abs'
                files(20)%fileName = 'datas1/rbm.501'
                files(21)%fileName = 'datas1/dna.setlow.new'
                files(22)%fileName = 'datas1/SCUP-h'
                files(23)%fileName = 'datas1/ery.anders'
                files(24)%fileName = 'datas1/acgih.1992'
                files(25)%fileName = 'datas1/phaeo.bio'
                files(26)%fileName = 'datas1/proro.bio'
                files(27)%fileName = 'datas1/cataract_oriowo'
                files(28)%fileName = 'dataj1/yld/O3.param_jpl97.yld'
                files(29)%fileName = 'dataj1/yld/O3.param.yld'
                files(30)%fileName = 'dataj1/yld/O3_shetter.yld'
                files(31)%fileName = 'datae1/no2/NO2_jpl94.abs'
                files(32)%fileName = 'datae1/no2/NO2_Har.abs'
                files(33)%fileName = 'dataj1/yld/NO2_calvert.yld'
                files(34)%fileName = 'dataj1/abs/NO3_gj78.abs'
                files(35)%fileName = 'dataj1/abs/NO3_jpl94.abs'
                files(36)%fileName = 'dataj1/abs/N2O5_jpl97.abs'
                files(37)%fileName = 'dataj1/abs/HNO2_jpl92.abs'
                files(38)%fileName = 'dataj1/abs/HNO3_burk.abs'
                files(39)%fileName = 'dataj1/abs/HNO4_jpl92.abs'
                files(40)%fileName = 'dataj1/abs/H2O2_jpl94.abs'
                files(41)%fileName = 'dataj1/abs/CHBr3.abs'
                files(42)%fileName = 'dataj1/abs/CHBr3.jpl97'
                files(43)%fileName = 'dataj1/ch2o/CH2O_nbs.abs'
                files(44)%fileName = 'dataj1/CH2O_iupac1.abs'
                files(45)%fileName = 'dataj1/ch2o/CH2O_can_hr.abs'
                files(46)%fileName = 'dataj1/ch2o/CH2O_can_lr.abs'
                files(47)%fileName = 'dataj1/ch2o/CH2O_rog.abs'
                files(48)%fileName = 'dataj1/ch2o/CH2O_ncar.abs'
                files(49)%fileName = 'dataj1/ch2o/CH2O_i_mad.yld'
                files(50)%fileName = 'dataj1/ch2o/CH2O_ii_mad.yld'
                files(51)%fileName = 'dataj1/ch2o/CH2O_iupac.yld'
                files(52)%fileName = 'dataj1/ch2o/CH2O_jpl97.dat'
                files(53)%fileName = 'dataj1/ch3cho/CH3CHO_iup.abs'
                files(54)%fileName = 'dataj1/ch3cho/d021_cp.abs'
                files(55)%fileName = 'dataj1/ch3cho/CH3CHO_mar.abs'
                files(56)%fileName = 'dataj2/kfa/ch3cho.005'
                files(57)%fileName = 'dataj1/ch3cho/CH3CHO_iup.yld'
                files(58)%fileName = 'dataj1/ch3cho/d021_i.yld'
                files(59)%fileName = 'dataj1/ch3cho/d021_ii.yld'
                files(60)%fileName = 'dataj1/ch3cho/d021_iii.yld'
                files(61)%fileName = 'dataj1/ch3cho/CH3CHO_press.yld'
                files(62)%fileName = 'dataj1/c2h5cho/C2H5CHO_iup.abs'
                files(63)%fileName = 'dataj2/kfa/c2h5cho.001'
                files(64)%fileName = 'dataj1/c2h5cho/C2H5CHO_iup.yld'
                files(65)%fileName = 'dataj1/chocho/CHOCHO_iup.abs'
                files(66)%fileName = 'dataj2/kfa/chocho.001'
                files(67)%fileName = 'dataj1/chocho/glyoxal_orl.abs'
                files(68)%fileName = 'dataj1/chocho/glyoxal_horowitz.abs'
                files(69)%fileName = 'dataj1/ch3cocho/CH3COCHO_iup1.abs'
                files(70)%fileName = 'dataj1/ch3cocho/CH3COCHO_iup2.abs'
                files(71)%fileName = 'dataj1/ch3cocho/CH3COCHO_ncar.abs'
                files(72)%fileName = 'dataj2/kfa/ch3cocho.001'
                files(73)%fileName = 'dataj2/kfa/ch3cocho.002'
                files(74)%fileName = 'dataj2/kfa/ch3cocho.003'
                files(75)%fileName = 'dataj2/kfa/ch3cocho.004'
                files(76)%fileName = 'dataj1/ch3cococh3/biacetyl_plum.abs'
                files(77)%fileName = 'dataj1/chocho/glyoxal_orl.abs'
                files(78)%fileName = 'dataj1/ch3cocho/CH3COCHO_km.yld'
                files(79)%fileName = 'dataj1/ch3coch3/CH3COCH3_cp.abs'
                files(80)%fileName = 'dataj1/ch3coch3/CH3COCH3_iup.abs'
                files(81)%fileName = 'dataj1/ch3coch3/CH3COCH3_noaa.abs'
                files(82)%fileName = 'dataj1/ch3coch3/CH3COCH3_iup.yld'
                files(83)%fileName = 'dataj1/ch3ooh/CH3OOH_jpl94.abs'
                files(84)%fileName = 'dataj1/ch3ooh/CH3OOH_iup.abs'
                files(85)%fileName = 'dataj1/ch3ooh/CH3OOH_ct.abs'
                files(86)%fileName = 'dataj1/ch3ooh/CH3OOH_ma.abs'
                files(87)%fileName = 'dataj1/rono2/CH3ONO2_cp.abs'
                files(88)%fileName = 'dataj1/rono2/CH3ONO2_tal.abs'
                files(89)%fileName = 'dataj1/rono2/CH3ONO2_iup1.abs'
                files(90)%fileName = 'dataj1/rono2/CH3ONO2_iup2.abs'
                files(91)%fileName = 'dataj1/rono2/CH3ONO2_tay.abs'
                files(92)%fileName = 'dataj1/rono2/CH3ONO2_rat.abs'
                files(93)%fileName = 'dataj1/rono2/CH3ONO2_lib.abs'
                files(94)%fileName = 'dataj1/rono2/PAN_talukdar.abs'
                files(95)%fileName = 'dataj1/abs/CCl2O_jpl94.abs'
                files(96)%fileName = 'dataj1/abs/CCl4_jpl94.abs'
                files(97)%fileName = 'dataj1/abs/CClFO_jpl94.abs'
                files(98)%fileName = 'dataj1/abs/CF2O_jpl94.abs'
                files(99)%fileName = 'dataj1/abs/CFC-113_jpl94.abs'
                files(100)%fileName = 'dataj1/abs/CFC-114_jpl94.abs'
                files(101)%fileName = 'dataj1/abs/CFC-115_jpl94.abs'
                files(102)%fileName = 'dataj1/abs/CFC-11_jpl94.abs'
                files(103)%fileName = 'dataj1/abs/CFC-12_jpl94.abs'
                files(104)%fileName = 'dataj1/abs/CH3Br_jpl94.abs'
                files(105)%fileName = 'dataj1/abs/CH3CCl3_jpl94.abs'
                files(106)%fileName = 'dataj1/abs/CH3Cl_jpl94.abs'
                files(107)%fileName = 'dataj1/abs/ClOO_jpl94.abs'
                files(108)%fileName = 'dataj1/abs/HCFCs_orl.abs'
                files(109)%fileName = 'dataj1/abs/HCFCs_orl.abs'
                files(110)%fileName = 'dataj1/abs/HCFC-141b_jpl94.abs'
                files(111)%fileName = 'dataj1/abs/HCFCs_orl.abs'
                files(112)%fileName = 'dataj1/abs/HCFC-225ca_jpl94.abs'
                files(113)%fileName = 'dataj1/abs/HCFC-225cb_jpl94.abs'
                files(114)%fileName = 'dataj1/abs/HCFC-22_jpl94.abs'
                files(115)%fileName = 'dataj1/abs/HO2_jpl94.abs'
                files(116)%fileName = 'dataj1/abs/Halon-1202_jpl97.abs'
                files(117)%fileName = 'dataj1/abs/Halon-1211_jpl97.abs'
                files(118)%fileName = 'dataj1/abs/Halon-1301_jpl97.abs'
                files(119)%fileName = 'dataj1/abs/Halon-2402_jpl97.abs'
                files(120)%fileName = 'dataj1/abs/ClONO2_jpl97.abs'
                files(121)%fileName = 'dataj1/abs/BrONO2_jpl03.abs'
                files(122)%fileName = 'dataj1/abs/CL2_fpp.abs'
                files(123)%fileName = 'dataj1/ch2ohcho/glycolaldehyde.abs'
                files(124)%fileName = 'dataj1/ch3cococh3/biacetyl_plum.abs'
                files(125)%fileName = 'dataj1/ch3cococh3/biacetyl_horowitz.abs'
                files(126)%fileName = 'dataj1/abs/methylvinylketone.abs'
                files(127)%fileName = 'dataj1/abs/methacrolein.abs'
                files(128)%fileName = 'dataj1/ch3cocooh/pyruvic_horowitz.abs'
                files(129)%fileName = 'dataj1/rono2/RONO2_talukdar.abs'
                files(130)%fileName = 'dataj1/rono2/RONO2_talukdar.abs'
                files(131)%fileName = 'dataj1/abs/CLOOCL_jpl02.abs'
                files(132)%fileName = 'dataj1/abs/Hydroxyacetone.abs'
                files(133)%fileName = 'dataj1/abs/BrO.jpl03'
                files(134)%fileName = 'dataj1/abs/Br2.abs'
                files(135)%fileName = 'output/out_11'
                files(136)%fileName = 'datae1/grids/combined.grid'
                files(137)%fileName = 'datae1/grids/fast_tuv.grid'
                files(138)%fileName = 'datae1/o2/effxstex.txt'
                files(139)%fileName = 'datae1/wmo85'
                files(140)%fileName = 'datae1/o3/O3.molina.abs'
                files(141)%fileName = 'datae1/wmo85'
                files(142)%fileName = 'datae1/o3/o3absqs.dat'
                files(143)%fileName = 'datae1/wmo85'
                files(144)%fileName = 'datae1/o3/O3_bass.abs'
                !srf  files(145)%filename='input/xreference.dat'
                files(145)%filename = 'input/'//trim(chemical_mechanism)//'xreference.dat'
            else
                call mpas_log_write('TUV: Reading list of input files from file: '//trim(filesName//'listFiles.dat'))

                !begin RMF
                ! filesname already contains a relative directory link to
                ! ./tables/tuvData/listFiles.dat
                ! this was changed to ./tables/tuvData/ and filesHome var
                ! can be defined locally.
                !end RMF

                !Read all files name from filesName file.
                open (kin, FILE=trim(filesName//'listFiles.dat'))
                read (kin, *) !Header line
                read (kin, FMT='(I3.3)') nFiles
                call mpas_log_write(message='TUV: Number of files: ',intArgs=(/nFiles/))

                read (kin, *) !Header line
                read (kin, FMT='(A)') filesHome

                filesHome = trim(filesName)

                call mpas_log_write(message='TUV: Home directory of files: '//trim(filesHome))

                read (kin, *) !Header line
                allocate (files(nFiles))
                do i = 1, nFiles
                    read (kin, FMT='(I3,1X,A50)') nOfFile, files(i)%fileName
                    call mpas_log_write(message='TUV: File  = '//trim(files(i)%fileName),intArgs=(/i/))
                end do
                !-srf - special treatment for xreference.dat
                !-      since it depends on the chemical mechanism
                do i = 1, nFiles
                    if (files(i)%fileName == 'input/xreference.dat') then
                        files(i)%fileName = 'input/'//trim(chemical_mechanism)//'xreference.dat'
                        call mpas_log_write(message='TUV: xreference file='//trim(files(i)%fileName))
                    end if
                end do
                !-srf- end
            end if

            !Putting the prefix of files
            do i = 1, nFiles
                files(i)%fileName = trim(filesHome)//trim(files(i)%fileName)
            end do

            pi = 2.*asin(1.0)
            twopi = 2.*pi
            rpd = pi/180.0

            ! wavelengths (creates wavelength grid: lower, center, upper of each bin)
            ! NOTE:  Wavelengths are in vacuum.  To use wavelengths in air, see
            ! Section 3 below, where you must set lrefr= .TRUE.
            call gridw(wstart, wstop, nwint)
            call readpol(coef)
            call ReadAll(nw, wl)

            !**** Correction for air-vacuum wavelength shift:
            ! The TUV code assumes that all working wavelengths are strictly IN-VACUUM. This is assumed for ALL
            ! spectral data including extraterrestrial fluxes, ozone (and other) absorption cross sections,
            ! and various weighting functons (action spectra, photolysis cross sections, instrument spectral
            ! response functions).

            !  Occasionally, users may want their results to be given for wavelengths measured IN-AIR.
            ! The shift between IN-VACUUM and IN-AIR wavelengths depends on the index of refraction
            ! of air, which in turn depends on the local density of air, which in turn depends on
            ! altitude, temperature, etc.
            !  Here, we provide users with the option to use a wavelength grid IN-AIR, at the air density
            ! corresponding to the output altitude, airout = airden(izout), by setting the logical variable
            ! lrefr = .TRUE.  (default is lrefr = .FALSE.).  The wavelengths specified in gridw.f will be assumed
            ! to be IN-AIR, and will be shifted here to IN-VACUUM values to carry out the calculatons.
            ! The actual radiative transfer calculations will be done strictly with IN-VACUUM values.
            ! If this shift is applied (i.e., if lrefr = .TRUE.), the wavelength grid will be shifted back to air
            ! values just before the output is written.
            !  Note:  if this option is used (lref = .TRUE.), the wavelength values will be correct ONLY at the
            ! selected altitude, iz = iout.  The wavelength shift will be INCORRECT at all other altitudes.
            !  Note:  This option cannot be changed interactively in the input table.  It must be changed here.

            ! ___ SECTION 4: READ SPECTRAL DATA ____________________________
            ! read (and grid) extra terrestrial flux data:
            call rdetfl(nw, wl)

            ! read cross section data for
            !  O2 (will overwrite at Lyman-alpha and SRB wavelengths
            !               see subroutine la_srb.f)
            !  O3 (temperature-dependent)
            !  SO2
            !  NO2
            call rdo2xs(nw, wl, o2xs1)
            call rdso2xs(nw, wl, so2xs)
            call rdno2xs(nw, wl, no2xs)

            !***** Spectral weighting functions
            ! (Some of these depend on temperature T and pressure P, and therefore
            !  on altitude z.  Therefore they are computed only after the T and P profiles
            !  are set above with subroutines settmp and setair.)
            ! Photo-physical   set in swphys.f (transmission functions)
            ! Photo-biological set in swbiol.f (action spectra)
            ! Photo-chemical   set in swchem.f (cross sections x quantum yields)
            ! Physical and biological weigthing functions are assumed to depend
            ! only on wavelength.
            ! Chemical weighting functions (product of cross-section x quantum yield)
            ! for many photolysis reactions are known to depend on temperature
            ! and/or pressure, and therefore are functions of wavelength and altitude.
            ! Output:
            ! from pphys & pbiol:  s(ks,kw) - for each weighting function slabel(ks)
            ! from pchem:  sj(kj,kz,kw) - for each reaction jlabel(kj)
            ! For pchem, need to know temperature and pressure profiles.
            call swphys(nw, wl, wc, ns, sw, slabel)
            call swbiol(nw, wl, wc, ns, sw, slabel)

            initialized = .true.

            isread = .false.
            firstcall = .true.
            call1 = .true.
            pass1 = .true.
            pass2 = .true.
            pass3 = .true.
            pass4 = .true.
            pass5 = .true.
            pass5 = .true.
            doinit = .true.
            doinit2 = .true.

            tol = 10.*d1mach(4)
            maxmsg = 100
            nummsg = 0
            msglim = .false.
            maxmsg2 = 50
            nummsg2 = 0

            !------------------------------------------
            !      Loads Chebyshev polynomial Coeff.
            !------------------------------------------
            call init_xs

            dither = 10.*r1mach(4)
            !** Must dither more on Cray (14-digit prec)
            if (dither < 1.e-10) dither = 10.*dither

        end if !(just for processor myNum=0)

        if (isMpi) then
            !MPI broadcast
            !call MPI_Bcast(i, 1, MPI_INTEGERKIND, source, dminfo % comm, mpi_ierr)
            call MPI_BCAST(xRef, kj, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(doReaction, kj, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(mrefr, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(lrefr, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(airout, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(coef, mz*ms*mp, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wl, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wc, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wu, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(f, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            !call MPI_BCAST(o2xs, kz*kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(o2xs1, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(so2xs, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(no2xs, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(mm_o3xs, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(sqt, maxsqt, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(pi, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(twopi, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(rpd, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(tol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(dmach, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(rmach, 4, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(gmu, nmug, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(gwt, nmug, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(ylmg, (maxstr + 1)*nmug, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(s226, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(s263, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(s298, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            !call MPI_BCAST(o3xs, kz*kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(s218, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(s228, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(s243, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(s295, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(c0, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(c1, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(c2, kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(yg, kw*totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(yg1, kw*totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(yg1n, kw*totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(yg2, kw*totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(yg3, kw*totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(yg4, kw*totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(yg5, kw*totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(coeff, 4*3*totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(c, maxtrm, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(hugeVar, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(powmax, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(powmin, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(tinyVar, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(maxmsg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(nummsg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(msglim, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(maxmsg2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(nummsg2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(tbar, totReact, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(ac, 20*17, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(bc, 20*17, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(nw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(nj, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(iw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(sw, ks*kw, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(initialized, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        end if
        print *, '-------------------------------------------------'!; call flush (6)

    end subroutine InitTuv

! Radiative transfer scheme:
!   nstr = number of streams
!           If nstr < 2, will use 2-stream Delta Eddington
!           If nstr > 1, will use nstr-stream discrete ordinates
! Location (geographic):
!   lat = LATITUDE (degrees, North = positive)
!   lon = LONGITUDE (degrees, East = positive)
!   esfact = 1. (Earth-sun distance = 1.000 AU)
! Vertical grid:
!   zstart = surface elevation above sea level, km
!   zstop = top of the atmosphere (exospheric), km
!   nz = number of vertical levels, equally spaced
!         (nz will increase by +1 if zout does not match altitude grid)
! Wavlength grid:
!   wstart = starting wavelength, nm
!   wstop  = final wavelength, nm
!   nwint = number of wavelength intervals, equally spaced
!            if nwint < 0, the standard atmospheric wavelength grid, not
!            equally spaced, from 120 to 735 nm, will be used. In this
!            case, wstart and wstop values are ignored.
! Surface condition:
!   alsurf = surface albedo, wavelength independent
! Column amounts of absorbers (in Dobson Units, from surface to space):
!           Vertical profile for O3 from USSA76.  For SO2 and NO2, vertical
!           concentration profile is 2.69e10 molec cm-3 between 0 and
!           1 km above sea level, very small residual (10/largest) above 1 km.
!   so2col = sulfur dioxide (SO2)
!   no2col = nitrogen dioxide (NO2)
! Cloud, assumed horizontally uniform, total coverage, single scattering
!          albedo = 0.9999, asymmetry factor = 0.85, indep. of wavelength,
!          and also uniform vertically between zbase and ztop:
!   taucld = vertical optical depth, independent of wavelength
!   zbase = altitude of base, km above sea level
!   ztop = altitude of top, km above sea level
! Aerosols, assumed vertical provile typical of continental regions from
!          Elterman (1968):
!   tauaer = aerosol vertical optical depth at 550 nm, from surface to space.
!            If negative, will default to Elterman's values (ca. 0.235
!            at 550 nm).
!   ssaaer = single scattering albedo of aerosols, wavelength-independent.
!   alpha = Angstrom coefficient = exponent for wavelength dependence of
!            tauaer, so that  tauaer1/tauaer2  = (w2/w1)**alpha.
! Directional components of radiation, weighting factors:
!   dirsun = direct sun
!   difdn = down-welling diffuse
!   difup = up-welling diffuse
!         e.g. use:
!         dirsun = difdn = 1.0, difup = 0 for total down-welling irradiance
!         dirsun = difdn = difup = 1.0 for actinic flux from all directions
!         dirsun = difdn = 1.0, difup = -1 for net irradiance

    subroutine Tuv(mynum, &
                   nstr, nz, zLevel, sza, albedo, so2col, no2col, &
                   dtcld, & !
                   omcld, & !
                   gcld, & !
                   dtaer, & !
                   omaer, & !
                   gaer, & !
                   alpha, dirsun, difdn, &
                   difup, tlev, tlay, airden, cair, co3, tco3, esfact, &
                   valj)

        implicit none

        integer, intent(IN) :: myNum
        !! Number of current processor/core
        integer, intent(IN) :: nstr             
        !! number of radiation streams
        integer, intent(IN) :: nz               
        !! Altitude grid (number of levels)
        real, intent(IN) :: zLevel(nz)     
        !! vector of altitude levels (in km)
        real, intent(IN) :: sza            
        !! Solar zenith angle and azimuth
        real, intent(IN) :: albedo(kw)     
        !! surface albedo
        real, intent(IN) :: so2col         
        !! Total columns of SO2 (Dobson Units)
        real, intent(IN) :: no2col         
        !! Total columns of NO2 (Dobson Units)
        real, intent(IN) :: dtcld(nz, kw)  
        !!
        real, intent(IN) :: omcld(nz, kw)  
        !!
        real, intent(IN) :: gcld(nz, kw)   
        !!
        real, intent(IN) :: dtaer(nz, kw)  
        !!
        real, intent(IN) :: omaer(nz, kw)  
        !!
        real, intent(IN) :: gaer(nz, kw)   
        !!
        real, intent(IN) :: alpha               
        !! Angstrom alpha
        real, intent(IN) :: dirsun              
        !! direct sun
        real, intent(IN) :: difdn               
        !! down-welling diffuse
        real, intent(IN) :: difup               
        !! up-welling diffuse
        real, intent(IN) :: tlev(nz)       
        !! temperature (K) at each specified altitude level
        real, intent(IN) :: tlay(nz)       
        !! temperature (K) at each specified altitude layer
        real, intent(IN) :: airden(nz)     
        !! air density (molec/cc) at each specified altitude
        real, intent(IN) :: cair(nz)       
        !! number of air molecules per cm^2 at each altitude
        real, intent(IN) :: co3(nz)        
        !! Total of molec O3 cm-2
        real, intent(IN) :: tco3(nz)       
        !! Total column o3 - molec cm-2
        real, intent(IN) :: esfact              
        !! sun-earth distance UA
        !
        !Output variables:
        !real, intent(OUT) :: rate(ks, nz)  
        !! Weighted irradiances (dose rates) W m-2
        real, intent(OUT) :: valj(kj, nz)  
        !! Photolysis coefficients (j-values)
        !real, intent(OUT) :: sirrad(nz, kw)
        !! Spectral irradiance, [W m-2 nm-1]

        real :: saflux(nz, kw)
        !! Spectral actinic flux, quanta s-1 nm-1 cm-2


        logical, parameter :: is2print = .false.
        !! To print a debug on each column
! Altitude grid
        integer :: iz, izout !LFR>izout not defined
! Extra terrestrial solar flux
        real    :: etf(kw)
        integer :: is, iob
        real    :: drdw
        integer :: ij, iaux, jaux
        real    :: djdw
        integer :: iw
        integer :: iAng
! Other user-defined variables here:
        integer :: nn1, nn2             !XUEXI
        real    :: adjcoe1(kj, nz), adjcoe2(kj, nz) !adjcoe(kj,kz),

        real, dimension(nz)         :: edir
        real, dimension(nz)         :: edn
        real, dimension(nz)         :: eup
        real, dimension(nz)         :: fdir
        real, dimension(nz)         :: fdn
        real, dimension(nz)         :: fup
        real, dimension(nz)         :: scol
        real, dimension(nz)         :: vcol
        integer, dimension(0:nz)    :: nid
        real, dimension(nz, kw)     :: dtO2
        real, dimension(nz, kw)     :: dtO3
        real, dimension(nz, kw)     :: o3xs
        real, dimension(kj, nz)     :: adjcoe
        real, dimension(nz, kw)     :: o2xs
        real, dimension(nz, kw)     :: dtRl
        real, dimension(nz, kw)     :: dtSo2
        real, dimension(nz, kw)     :: dtNo2
        real, dimension(0:nz, nz)   :: dsdh
        real, dimension(kj, nz, kw) :: sj

        real :: valja

! print *,'LFR-DBG TUV albedo: ',minval(albedo(:)),maxval(albedo(:))
! print *,'LFR-DBG TUV sza: ',sza
! print *,'LFR-DBG TUV so2col: ',so2col
! print *,'LFR-DBG TUV no2col: ',no2col
! print *,'LFR-DBG TUV Tco3: ',minval(tco3(:)),maxval(tco3(:)) 
! print *,'LFR-DBG TUV airden: ',minval(airden(:)),maxval(airden(:))
! print *,'LFR-DBG TUV cAir: ',minval(cAir(:)),maxval(cAir(:))
! print *,'LFR-DBG TUV tLev: ',minval(tlev(:)),maxval(tlev(:))
! print *,'LFR-DBG TUV tLay: ',minval(tlay(:)),maxval(tlay(:))
! print *,'LFR-DBG TUV zLevel: ',minval(zLevel(:)),maxval(zLevel(:));call flush(6)
! print *,'LFR-DBG TUV dtcld: ',minval(dtCld(:,:)),maxval(dtCld(:,:))
! print *,'LFR-DBG TUV omcld: ',minval(omCld(:,:)),maxval(omCld(:,:))
! print *,'LFR-DBG TUV gcld: ',minval(gCld(:,:)),maxval(gCld(:,:))
! print *,'LFR-DBG TUV dtaer: ',minval(dtAer(:,:)),maxval(dtAer(:,:))
! print *,'LFR-DBG TUV omaer: ',minval(omAer(:,:)),maxval(omAer(:,:))
! print *,'LFR-DBG TUV gaer: ',minval(gAer(:,:)),maxval(gAer(:,:))
! print *,'LFR-DBG TUV alpha: ',alpha
! print *,'LFR-DBG TUV dirsun: ',dirsun
! print *,'LFR-DBG TUV difdn: ',difdn
! print *,'LFR-DBG TUV difup: ',difup
! print *,'LFR-DBG TUV airden(nbl,kz): ',minval(airden(:)),maxval(airden(:))
! print *,'LFR-DBG TUV cair(nbl,kz): ',minval(cAir(:)),maxval(cAir(:))
! print *,'LFR-DBG TUV co3(nbl,kz): ',minval(co3(:)),maxval(co3(:))
! print *,'LFR-DBG TUV tco3(nbl,kz): ',minval(tco3(:)),maxval(tco3(:))
! print *,'LFR-DBG TUV esfact: ',esfact
        current_nz = nz

        ! correction for earth-sun distance
        !default is 1.0 UA - must be set for each season
        do iw = 1, nw - 1
            etf(iw) = f(iw)*esfact
        end do

        !rate = 0.0
        valJ = 0.0
        !sirRad = 0.0
        saFlux = 0.0
        adjcoe = 1.0
!print *, 'LFR-DBG TUV: 001'
        ! Calculating sj values of the reaction
        call swchem(nw, wl, nz, nj, tLev, airDen, sj)
!print *, 'LFR-DBG TUV: 002'
        !Ozone molecular absorption cross section
!print *, 'LFR-DBG TUV: 002.5 nz, tlay=', nz, tLay(:)
!print *, 'LFR-DBG TUV: 002.5 nw=', nw, 'wl=', wl(:)
        call rdo3xs(nw, wl, nz, tLay, o3xs)
!print *, 'LFR-DBG TUV: 003'
        ! Rayleigh optical depth increments:
        call odrl(nz, nw, wl, cAir, dtRl)
!print *, 'LFR-DBG TUV: 004'
        ! O2 vertical profile and O2 absorption optical depths
        ! For now, O2 densitiy assumed as 20.95% of air density, can be change
        ! in subroutine.
        ! Optical depths in Lyman-alpha and SRB will be over-written
        ! in subroutine la_srb.f
        do iz = 1, nz
            do iw = 1, nw - 1
                if (iz > nz) cycle
                dtO2(iz, iw) = 0.2095*cAir(iz)*o2xs1(iw)
            end do
        end do
        ! Ozone optical depths
        do iw = 1, nw - 1
            do iz = 1, nz - 1
                dtO3(iz, iw) = co3(iz)*o3xs(iz, iw)
            end do
        end do
!print *, 'LFR-DBG TUV: 005'
        ! SO2 vertical profile and optical depths
        call setso2(nz, nw, so2xs, zLevel, so2col, cAir, dtSo2)
!print *, 'LFR-DBG TUV: 006'
        ! NO2 vertical profile and optical depths
        call setno2(nz, nw, no2xs, zLevel, no2col, cAir, dtNo2)
!print *, 'LFR-DBG TUV: 007'
        ! slant path lengths for spherical geometry
        call sphers(nz, zLevel, sza, nid, dsdh)
!print *, 'LFR-DBG TUV: 008'
        ! Calculate vertical and slant air columns
        call airmas(nz, cAir, scol, vcol, nid, dsdh)
!print *, 'LFR-DBG TUV: 009'
        !Obtain index for each angle
        iAng = int(sza/20.0 + 1)
        !Adjust for each angle
        if (iang == 1) then
            call setz(nz, nj, coef, adjcoe1, iAng, tLev(1), tcO3)
            do ij = 1, kj
                do iz = 1, nz
                    adjCoe(ij, iz) = adjcoe1(ij, iz)
                end do
            end do
        else if (iang < 5) then
            call setz(nz, nj, coef, adjcoe1, iAng, tLev(1), tcO3)
            call setz(nz, nj, coef, adjcoe2, iAng + 1, tLev(1), tcO3)
            do ij = 1, kj
                do iz = 1, nz
                    adjCoe(ij, iz) = adjcoe1(ij, iz) + &
                                      (adjcoe2(ij, iz) - adjcoe1(ij, iz)) &
                                      *(sza - angle(iAng))/(20.0)
                end do
            end do
        elseif (iang == 5) then
            call setz(nz, nj, coef, adjcoe1, iAng, tLev(1), tcO3)
            if (abs(sza) < 90) then
                do ij = 1, kj
                    do iz = 1, nz
                        adjCoe(ij, iz) = adjcoe1(ij, iz)
                    end do
                end do
            end if
        end if
        if (iAng > 5) then
            do ij = 1, kj
                do iz = 1, nz
                    adjCoe(ij, iz) = 1.0
                end do
            end do
        end if
!print *, 'LFR-DBG TUV: 010'
        ! Recalculate effective O2 optical depth and cross sections for Lyman-alpha
        ! and Schumann-Runge bands, must know zenith angle
        ! Then assign O2 cross section to sj(1,*,*)
        call la_srb(nz, nw, wl, o2xs1, tLev, dtO2, o2xs, scol, vcol)
        ! Update the weighting function (cross section x quantum yield) for O2
        ! photolysis.  The strong spectral variations in the O2 cross sections are
        ! parameterized into a few bands for Lyman-alpha (121.4-121.9 nm, one band)
        ! and Schumann-Runge (174.4-205.8, 17 bands) regions. The parameterizations
        ! depend on the overhead O2 column, and therefore on altitude and solar
        ! zenith angle, so they need to be updated at each time/zenith step.
!print *, 'LFR-DBG TUV: 011'
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(1, iz, iw) = o2xs(iz, iw)
            end do
        end do
!print *, 'LFR-DBG TUV: 012'
        !Main wavelength loop:
        do iw = 1, nw - 1
            !* monochromatic radiative transfer. Outputs are:
            !  normalized irradiances edir(iz), edn(iz), eup(iz)
            !  normalized actinic fluxes  fdir(iz), fdn(zi), fup(iz)
            !  where
            !  dir = direct beam, dn = down-welling diffuse,
            !   up = up-welling diffuse
!print *, 'LFR-DBG TUV: 012.0, iw=', iw,albedo(iw)
            call rtlink(nstr, nz, iw, sza, albedo, &
                        dtCld, omCld, gCld, dtAer, omAer, gAer, dtO2, dtO3, &
                        eDir, eDn, eUp, fDir, fDn, fUp, dtRl, dtSo2, dtNo2, &
                        nid, dsdh)
!print *, 'LFR-DBG TUV: 012.1'
            ! Spectral irradiance, W m-2 nm-1, down-welling:
            ! do iz = 1, nz
            !     sirRad(iz, iw) = etf(iw)* &
            !                      (dirsun*eDir(iz) + &
            !                       difdn*eDn(iz) + difup*eUp(iz))
            ! end do
!print *, 'LFR-DBG TUV: 012.2'
            ! Spectral actinic flux, quanta s-1 nm-1 cm-2, all directions:
            !     units conversion:  1.e-4 * (wc*1e-9) / (hc = 6.62E-34 * 2.998E8)
             do iz = 1, nz
                 saFlux(iz, iw) = &
                     etf(iw)*5.039e11*wc(iw)* &
                     (dirsun*fdir(iz) + &
                      difdn*fdn(iz) + difup*fup(iz))
            end do
!print *, 'LFR-DBG TUV: 012.3'
            !** Accumulate weighted integrals over wavelength, at all altitudes:
!            do iz = 1, nz
!                ! Weighted irradiances (dose rates) W m-2
!                do is = 1, ns
!                    drdw = sirRad(iz, iw)*sw(is, iw)
!                    rate(is, iz) = rate(is, iz) + drdw*(wu(iw) - wl(iw))
!                end do
!            end do

            do ij = 1, kj
!print *, 'LFR-DBG TUV: 012.4 , ij =', ij
                if (doReaction(ij)) then
                    do iz = 1, nz
!print *, 'LFR-DBG TUV: 012.5, ij=', ij, ' iz=', iz, valJ(ij, iz), saFlux(iz, iw), sj(ij, iz, iw), adjCoe(ij, iz),wu(iw), wl(iw)
                        ! Photolysis rate coefficients (J-values) s-1
                        djdw = saFlux(iz, iw) * sj(ij, iz, iw) * adjCoe(ij, iz)  ! XUEXI
                        valJ(ij, iz) = valJ(ij, iz) + djdw * (wu(iw) - wl(iw))
                        valJ(ij, iz) = valJ(ij, iz)/100.0 
!LFR-DBG only for test
                        write(51,fmt='(I3.3,1X,I3.3,1X,E15.5)') ij,iz,valJ(ij, iz)
!print *, 'LFR-DBG TUV: 012.6 valJ ',valJ(ij, iz)
                    end do
                end if
            end do
        end do
        !end wavelength loop
!print *, 'LFR-DBG TUV: 013'


        !* reset wavelength scale if needed:
        if (lrefr) then
            write (*, *) 'applying vacuum to air wavelength shift', airout
            mrefr = -mrefr
            call wshift(mrefr, nw, wl, airout)
            call wshift(mrefr, nwint, wc, airout)
            call wshift(mrefr, nwint, wu, airout)
        end if

901     format('zenith =   ', f10.1)
902     format(i10, a30)
903     format(i10, 4e13.3)

    end subroutine tuv

    subroutine calcoe(nz, ij, c, xzin, adjin, adjcoe)
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  ADJCOE - REAL, coross section adjust coefficients (in and out)           =*
!=  c(5,kj)- polynomal coef                                                   =*
!=  tt     - nomarlized temperature
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  11/2000 XUEXI                                                            =*
!-----------------------------------------------------------------------------*
!= This program is free software;  you can redistribute it and/or modify     =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                 =*
!= To obtain a copy of the GNU General Public License, write to:             =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= To contact the authors, please mail to:                                   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu                                            =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nz
        integer, intent(IN) :: ij
        real, intent(IN)    :: c(5, kj)
        real, intent(IN)    :: xzin(nz)
        real, intent(IN)    :: adjin
        real, intent(OUT)   :: adjcoe(kj, nz)
        real :: x2, x3, x4

        integer :: k
        real :: xz(nz)

        do k = 1, nz
            xz(k) = xzin(k)*adjO3(ij)
            x2 = xz(k)*xz(k)
            x3 = x2*xz(k)
            x4 = x3*xz(k)
            adjcoe(ij, k) = adjin*(1.0 + (c(1, ij) &
                                          + c(2, ij)*xz(k) + c(3, ij)*x2 &
                                          + c(4, ij)*x3 + c(5, ij)*x4)*0.01)
        end do

    end subroutine calcoe

    real function fery(w)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the action spectrum value for erythema at a given wavelength   =*
!=  according to: McKinlay, A.F and B.L.Diffey, A reference action spectrum  =*
!=  for ultraviolet induced erythema in human skin, CIE Journal, vol 6,      =*
!=  pp 17-22, 1987.                                                          =*
!=  Value at 300 nm = 0.6486                                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  W - REAL, wavelength (nm)                                             (I)=*
!-----------------------------------------------------------------------------*

        real, intent(IN)                     :: w

        if (w < 250.) then
            fery = 1.
! outside the ery spectrum range
        else if ((w >= 250.) .and. (w < 298)) then
            fery = 1.
        else if ((w >= 298.) .and. (w < 328.)) then
            fery = 10.**(0.094*(298.-w))
        else if ((w >= 328.) .and. (w < 400.)) then
            fery = 10.**(0.015*(139.-w))
        else
            fery = 1.e-36
! outside the ery spectrum range
        end if

    end function fery

!=============================================================================*

    real function fo3qy(w, t)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
! function to calculate the quantum yield O3 + hv -> O(1D) + O2,             =*
! according to JPL 2000 recommendation:                                      =*
!-----------------------------------------------------------------------------*

        real, intent(IN)  :: w
        real, intent(IN)  :: t
        real :: kt
        real, parameter :: a(3) = (/0.887, 2.35, 57.0/)
        real, parameter :: w0(3) = (/302.0, 311.1, 313.9/)
        real, parameter :: nu(3) = (/0.0, 820.0, 1190.0/)
        real, parameter :: om(3) = (/7.9, 2.2, 7.4/)

        fo3qy = 0.
        kt = 0.695*t

        if (w <= 300.) then
            fo3qy = 0.95
        else if (w > 300. .and. w <= 330.) then
            fo3qy = 0.06 + &
                    a(1)*exp(-((w - w0(1))/om(1))**4) + &
                    a(2)*(t/300.)**4*exp(-nu(2)/kt)*exp(-((w - w0(2))/om(2))**2) + &
                    a(3)*exp(-nu(3)/kt)*exp(-((w - w0(3))/om(3))**2)
        else if (w > 330. .and. w <= 345.) then
            fo3qy = 0.06
        else if (w > 345.) then
            fo3qy = 0.
        end if

    end function fo3qy

    real function fo3qy2(w, t)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
! function to calculate the quantum yield O3 + hv -> O(1D) + O2,             =*
! according to:
! Matsumi, Y., F. J. Comes, G. Hancock, A. Hofzumanhays, A. J. Hynes,
! M. Kawasaki, and A. R. Ravishankara, QUantum yields for production of O(1D)
! in the ultraviolet photolysis of ozone:  Recommendation based on evaluation
! of laboratory data, J. Geophys. Res., 107, 10.1029/2001JD000510, 2002.
!-----------------------------------------------------------------------------*

        real, intent(IN)  :: w
        real, intent(IN)  :: t
        real :: kt

        real, parameter :: a(3) = (/0.8036, 8.9061, 0.1192/)
        real, parameter :: x(3) = (/304.225, 314.957, 310.737/)
        real, parameter :: om(3) = (/5.576, 6.601, 2.187/)
        real :: q1, q2

        fo3qy2 = 0.0
        kt = 0.695*t
        q1 = 1.0
        q2 = exp(-825.518/kt)

        if (w <= 305.) then
            fo3qy2 = 0.90
        else if (w > 305. .and. w <= 328.) then

            fo3qy2 = 0.0765 + &
                     a(1)*(q1/(q1 + q2))*exp(-((x(1) - w)/om(1))**4) + &
                     a(2)*(t/300.)**2*(q2/(q1 + q2))*exp(-((x(2) - w)/om(2))**2) + &
                     a(3)*(t/300.)**1.5*exp(-((x(3) - w)/om(3))**2)

        else if (w > 328. .and. w <= 340.) then
            fo3qy2 = 0.08
        else if (w > 340.) then
            fo3qy2 = 0.
        end if

    end function fo3qy2

!=============================================================================*

    real function fsum(n, x)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute the sum of the first N elements of a floating point vector.      =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  N  - INTEGER, number of elements to sum                               (I)=*
!=  X  - REAL, vector whose components are to be summed                   (I)=*
!-----------------------------------------------------------------------------*
        integer, intent(IN) :: n
        real, intent(IN)    :: x(n)

! local:
        integer :: i

        fsum = 0.
        do i = 1, n
            fsum = fsum + x(i)
        end do

    end function fsum

!=============================================================================*

    real function futr(w)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the action spectrum value for skin cancer of albino hairless   =*
!=  mice at a given wavelength according to:  deGRuijl, F.R., H.J.C.M.Steren-=*
!=  borg, P.D.Forbes, R.E.Davies, C.Colse, G.Kelfkens, H.vanWeelden,         =*
!=  and J.C.van der Leun, Wavelength dependence of skin cancer induction by  =*
!=  ultraviolet irradiation of albino hairless mice, Cancer Research, vol 53,=*
!=  pp. 53-60, 1993                                                          =*
!=  (Action spectrum for carcinomas)                                         =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  W  - REAL, wavelength (nm)                                            (I)=*
!-----------------------------------------------------------------------------*

        real, intent(IN)                     :: w

! local:
        real :: a1, a2, a3, a4, a5, x1, x2, x3, x4, x5, &
                t1, t2, t3, t4, t5, b1, b2, b3, b4, b5, &
                p

        a1 = -10.91
        a2 = -0.86
        a3 = -8.60
        a4 = -9.36
        a5 = -13.15

        x1 = 270.
        x2 = 302.
        x3 = 334.
        x4 = 367.
        x5 = 400.

        t1 = (w - x2)*(w - x3)*(w - x4)*(w - x5)
        t2 = (w - x1)*(w - x3)*(w - x4)*(w - x5)
        t3 = (w - x1)*(w - x2)*(w - x4)*(w - x5)
        t4 = (w - x1)*(w - x2)*(w - x3)*(w - x5)
        t5 = (w - x1)*(w - x2)*(w - x3)*(w - x4)

        b1 = (x1 - x2)*(x1 - x3)*(x1 - x4)*(x1 - x5)
        b2 = (x2 - x1)*(x2 - x3)*(x2 - x4)*(x2 - x5)

        b3 = (x3 - x1)*(x3 - x2)*(x3 - x4)*(x3 - x5)
        b4 = (x4 - x1)*(x4 - x2)*(x4 - x3)*(x4 - x5)
        b5 = (x5 - x1)*(x5 - x2)*(x5 - x3)*(x5 - x4)

        p = a1*t1/b1 + a2*t2/b2 + a3*t3/b3 + a4*t4/b4 + a5*t5/b5

        futr = exp(p)

    end function futr

    subroutine gridw(wstart, wstop, nwint)

        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Create the wavelength grid for all interpolations and radiative transfer =*
!=  calculations.  Grid may be irregularly spaced.  Wavelengths are in nm.   =*
!=  No gaps are allowed within the wavelength grid.                          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW  - INTEGER, number of wavelength grid _points_                     (O)=*
!=  WL  - REAL, vector carrying the lower limit of each wavel. interval   (O)=*
!=  WC  - REAL, vector carrying the center wavel of each wavel. interval  (O)=*
!=              (wc(i) = 0.5*(wl(i)+wu(i), i = 1..NW-1)                      =*
!=  WU  - REAL, vector carrying the upper limit of each wavel. interval   (O)=*
!=
!=  MOPT- INTEGER OPTION for wave-length IF 3 good for JO2                (O)=*
!-----------------------------------------------------------------------------*

        real, intent(IN)     :: wstart
        real, intent(IN)     :: wstop
        integer, intent(IN)  :: nwint

        integer :: mopt
        real :: wincr
        integer :: iw

        character(LEN=200) :: fi
        character(LEN=20) :: wlabel

        real :: dum

        logical :: ok
!_______________________________________________________________________

!*** chose wavelength grid

! some pre-set options
!     mopt = 1    equal spacing
!     mopt = 2    grid defined in data table
!     mopt = 3    user-defined
!     mopt = 4    fast-TUV, troposheric wavelengths only
        mopt = 1
        if (nwint == -156) mopt = 2
        if (nwint <= -1 .and. nwint >= -20) mopt = 4     ! fast-J XUEXI

        select case (mopt)
        case (1)
            wlabel = 'equal spacing'
            nw = nwint + 1
            wincr = (wstop - wstart)/FLOAT(nwint)
            do iw = 1, nw - 1
                wl(iw) = wstart + wincr*FLOAT(iw - 1)
                wu(iw) = wl(iw) + wincr
                wc(iw) = (wl(iw) + wu(iw))/2.
            end do
            wl(nw) = wu(nw - 1)
        case (2)
! Input from table.  In this example:
! Wavelength grid will be read from a file.
! First line of table is:  nw = number of wavelengths (no. of intervals + 1)
! Then, nw wavelengths are read in, and assigned to wl(iw)
! Finally, wu(iw) and wc(iw) are computed from wl(iw)

!      wlabel = 'isaksen.grid'
!wlabel = 'combined.grid'

            fi = trim(files(136)%fileName)
            open (UNIT=kin, FILE=fi, STATUS='old')
            read (kin, *) nw
            do iw = 1, nw
                read (kin, *) wl(iw)
            end do
            close (kin)
            do iw = 1, nw - 1
                wu(iw) = wl(iw + 1)
                wc(iw) = 0.5*(wl(iw) + wu(iw))
            end do
        case (3)
! user-defined grid.  In this example, a single calculation is used to
! obtain results for two 1 nm wide intervals centered at 310 and 400 nm:
! interval 1 : 1 nm wide, centered at 310 nm
! interval 3 : 2 nm wide, centered at 400 nm
! (inteval 2 : 310.5 - 399.5 nm, required to connect intervals 1 & 3)

            nw = 4
            wl(1) = 309.5
            wl(2) = 310.5
            wl(3) = 399.5
            wl(4) = 400.5
            do iw = 1, nw - 1
                wu(iw) = wl(iw + 1)
                wc(iw) = 0.5*(wl(iw) + wu(iw))
            end do
        case (4)
            wlabel = 'fast-TUV tropospheric grid'

            fi = trim(files(137)%fileName)
            open (UNIT=kin, FILE=fi, STATUS='old')
            do iw = 1, 4
                read (kin, *)
            end do
!PRINT *, 'LFR->',fi
! skip wavelength shorter than 205 nm
            do iw = 1, 6
                read (kin, *)
            end do
            nw = abs(nwint) + 1
            do iw = 1, nw - 1
                read (kin, *) dum, wl(iw), dum, dum
            end do
            wl(nw) = dum
            do iw = 1, nw - 1
                wu(iw) = wl(iw + 1)
                wc(iw) = 0.5*(wl(iw) + wu(iw))
            end do
            close (kin)
        end select
        call gridck(kw, nw, wl, ok)

        if (.not. ok) then
            write (*, *) 'STOP in GRIDW:  The w-grid does not make sense'
            stop
        end if

!_______________________________________________________________________

    end subroutine gridw

    subroutine gridck(k, n, x, ok)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Check a grid X for various improperties.  The values in X have to comply =*
!=  with the following rules:                                                =*
!=  1) Number of actual points cannot exceed declared length of X            =*
!=  2) Number of actual points has to be greater than or equal to 2          =*
!=  3) X-values must be non-negative                                         =*
!=  4) X-values must be unique                                               =*
!=  5) X-values must be in ascending order                                   =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  K  - INTEGER, length of X as declared in the calling program          (I)=*
!=  N  - INTEGER, number of actual points in X                            (I)=*
!=  X  - REAL, vector (grid) to be checked                                (I)=*
!=  OK - LOGICAL, .TRUE. -> X agrees with rules 1)-5)                     (O)=*
!=                .FALSE.-> X violates at least one of 1)-5)                 =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN)    :: k
        integer, intent(IN)    :: n
        real, intent(IN)       :: x(k)
        logical, intent(OUT)   :: ok
! local:
        integer :: i
        integer, parameter :: kout = 6
!_______________________________________________________________________

        ok = .true.

! check if dimension meaningful and within bounds

        if (n > k) then
            ok = .false.
            write (kout, 100)
            return
        end if
100     format('Number of data exceeds dimension')

        if (n < 2) then
            ok = .false.
            write (kout, 101)
            return
        end if
101     format('Too few data, number of data points must be >= 2')

! disallow negative grid values

        if (x(1) < 0.) then
            ok = .false.
            write (kout, 105)
            return
        end if
105     format('Grid cannot start below zero')

! check sorting

        do i = 2, n
            if (x(i) <= x(i - 1)) then
                ok = .false.
                write (kout, 110)
                return
            end if
        end do
110     format('Grid is not sorted or contains multiple values')
!_______________________________________________________________________

    end subroutine gridck

    subroutine la_srb(nz, nw, wl, o2xs1, tLev, dtO2, o2xs, scol, vcol)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute equivalent optical depths for O2 absorption, and O2 effective    =*
!=  absorption cross sections, parameterized in the Lyman-alpha and SR bands =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  Z       - REAL, specified altitude working grid (km)                  (I)=*
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
!=            wavelength grid                                                =*
!=  WL      - REAL, vector of lxower limits of wavelength intervals in    (I)=*
!=            working wavelength grid                                        =*
!=  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)=*
!=            altitude layer                                                 =*
!=  ZEN     - REAL, solar zenith angle                                    (I)=*
!=                                                                           =*
!=  O2XS1   - REAL, O2 cross section from rdo2xs                          (I)=*
!=                                                                           =*
!=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer at each specified wavelength                    =*
!=  O2XS    - REAL, molecular absorption cross section in SR bands at     (O)=*
!=            each specified altitude and wavelength.  Includes Herzberg     =*
!=            continuum.                                                     =*
!-----------------------------------------------------------------------------*
        integer, intent(IN) :: nz
        integer, intent(IN) :: nw
        real, intent(IN)    :: wl(kw)
        real, intent(IN)    :: o2xs1(kw)
        real, intent(IN)    :: tLev(nz)
        real, intent(IN)    :: scol(nz)
        real, intent(IN)    :: vcol(nz)
        real, intent(INOUT) :: dtO2(nz, kw)
        real, intent(OUT)   :: o2xs(nz, kw)

        integer :: iz, iw
        real, allocatable, dimension(:) :: o2col
        real, allocatable, dimension(:) :: secchi

! Lyman-alpha variables
! O2 optical depth and equivalent cross section in the Lyman-alpha region

        integer, parameter :: nla = 1
        integer, parameter :: kla = 2
        real, dimension(kla), parameter :: wlla = (/121.4, 121.9/) ! Wavelengths for Lyman alpha and SRB parameterizations
        real, allocatable, dimension(:, :) :: dto2la
        real, allocatable, dimension(:, :) :: o2xsla

! grid on which Koppers' parameterization is defined
! O2 optical depth and equivalent cross section on Koppers' grid

        integer, parameter :: nsrb = 17
        integer, parameter :: ksrb = 18
        real, dimension(ksrb), parameter :: wlsrb = (/ &
                                            174.4, 177.0, 178.6, 180.2, 181.8, 183.5, 185.2, 186.9, &
                                            188.7, 190.5, 192.3, 194.2, 196.1, 198.0, 200.0, 202.0, &
                                            204.1, 205.8/)
        real, allocatable, dimension(:, :) :: dto2k
        real, allocatable, dimension(:, :) :: o2xsk

        integer :: i

!----------------------------------------------------------------------
! initalize O2 cross sections
!----------------------------------------------------------------------
        do iz = 1, nz
            do iw = 1, nw - 1
                if (iz > nz) cycle
                o2xs(iz, iw) = o2xs1(iw)
            end do
        end do

        if (wl(1) > wlsrb(nsrb)) return

!----------------------------------------------------------------------
! On first call, check that the user wavelength grid, WL(IW), is compatible
! with the wavelengths for the parameterizations of the Lyman-alpha and SRB.
! Also compute and save corresponding grid indices (ILA, ISRB)
!----------------------------------------------------------------------
        if (call1) then
!* locate Lyman-alpha wavelengths on grid
            ila = 0
            do iw = 1, nw
                if (abs(wl(iw) - wlla(1)) < 10.*precis) then
                    ila = iw
                    exit
                end if
            end do
! check
            if (ila == 0) stop ' Lyman alpha grid mis-match - 1'
            do i = 2, nla + 1
                if (abs(wl(ila + i - 1) - wlla(i)) > 10.*precis) then
                    write (*, *) 'Lyman alpha grid mis-match - 2'
                    stop
                end if
            end do
!* locate Schumann-Runge wavelengths on grid
            isrb = 0
            do iw = 1, nw
                if (abs(wl(iw) - wlsrb(1)) < 10.*precis) then
                    isrb = iw
                    exit
                end if
            end do
! check
            if (isrb == 0) stop ' SRB grid mis-match - 1'
            do i = 2, nsrb + 1
                if (abs(wl(isrb + i - 1) - wlsrb(i)) > 10.*precis) then
                    write (*, *) ' SRB grid mismatch - w'
                    stop
                end if
            end do
            call1 = .false.
        end if

!Local allocations
        allocate (o2col(nz))
        allocate (secchi(nz))
        allocate (dto2la(nz, kla - 1))
        allocate (o2xsla(nz, kla - 1))
        allocate (dto2k(nz, ksrb - 1))
        allocate (o2xsk(nz, ksrb - 1))
!----------------------------------------------------------------------
! Slant O2 column and x-sections.
!----------------------------------------------------------------------
        do iz = 1, nz
            o2col(iz) = 0.2095*scol(iz)
        end do

!----------------------------------------------------------------------
! Effective secant of solar zenith angle.
! Use 2.0 if no direct sun (value for isotropic radiation)
! For nz, use value at nz-1
!----------------------------------------------------------------------
        do iz = 1, nz - 1
            secchi(iz) = scol(iz)/vcol(iz)
            if (scol(iz) > largest*0.1) secchi(iz) = 2.
        end do
        secchi(nz) = secchi(nz - 1)
!---------------------------------------------------------------------
! Lyman-Alpha parameterization, output values of O2 optical depth
! and O2 effective (equivalent) cross section
!----------------------------------------------------------------------
        call lymana(nz, o2col, secchi, dto2la, o2xsla)
        do iw = ila, ila + nla - 1
            do iz = 1, nz
                dtO2(iz, iw) = dto2la(iz, iw - ila + 1)
                o2xs(iz, iw) = o2xsla(iz, iw - ila + 1)
            end do
        end do

!------------------------------------------------------------------------------
! Koppers' parameterization of the SR bands, output values of O2
! optical depth and O2 equivalent cross section
!------------------------------------------------------------------------------

        call schum(nz, o2col, secchi, dto2k, o2xsk, tLev)
        do iw = isrb, isrb + nsrb - 1
            do iz = 1, nz
                dto2(iz, iw) = dto2k(iz, iw - isrb + 1)
                o2xs(iz, iw) = o2xsk(iz, iw - isrb + 1)
            end do
        end do

    end subroutine la_srb

!=============================================================================*

    subroutine lymana(nz, o2col, secchi, dto2la, o2xsla)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha=*
!=  bands and an effective O2 optical depth at all altitudes.  Parameterized =*
!=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the  =*
!=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,  =*
!=  Vol.24, No.21, pp 2659-2662, 1997.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
!=            altitude                                                       =*
!=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer                                                 =*
!=  O2XSLA  - REAL, molecular absorption cross section in LA bands        (O)=*
!-----------------------------------------------------------------------------*
        integer, intent(IN) :: nz
        real, intent(IN)    :: o2col(nz)
        real, intent(IN)    :: secchi(nz)
        real, intent(OUT)   :: dto2la(nz, *)
        real, intent(OUT)   :: o2xsla(nz, *)

        double precision :: rm, ro2

        double precision, dimension(3), parameter :: bbb = (/6.8431d-01, 2.29841d-01, 8.65412d-02/)
        double precision, dimension(3), parameter :: ccc = (/8.22114d-21, 1.77556d-20, 8.22112d-21/)
        double precision, dimension(3), parameter :: ddd = (/6.0073d-21, 4.28569d-21, 1.28059d-20/)
        double precision, dimension(3), parameter :: eee = (/8.21666d-21, 1.63296d-20, 4.85121d-17/)

        integer :: iz, i
        real :: xsmin

!------------------------------------------------------------------------------*
!sm:  set minimum cross section
        xsmin = 1.e-20

        do iz = 1, nz
            rm = 0.0d+00
            ro2 = 0.0d+00
            do i = 1, 3
                rm = rm + bbb(i)*DEXP(-ccc(i)*dble(o2col(iz)))
                ro2 = ro2 + ddd(i)*DEXP(-eee(i)*dble(o2col(iz)))
            end do
            if (rm > 1.0d-100) then
                if (ro2 > 1.0d-100) then
                    o2xsla(iz, 1) = ro2/rm
                else
                    o2xsla(iz, 1) = xsmin
                end if
            else
                o2xsla(iz, 1) = xsmin
            end if
        end do

! calculate effective O2 optical depths and effective O2 cross sections
        do iz = 1, nz - 1
            if (rm > 1.0d-100) then
                if (rm > 0.) then
                    dto2la(iz, 1) = log(rm)/secchi(iz) - log(rm)/secchi(iz+1)
                else
                    dto2la(iz, 1) = 1000.
                end if
            else
                dto2la(iz, 1) = 1000.
            end if
        end do

! do top layer separately
        dto2la(nz, 1) = 0.

    end subroutine lymana

    subroutine schum(nz, o2col, secchi, dto2, o2xsk, tLev)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the equivalent absorption cross section of O2 in the SR bands. =*
!=  The algorithm is based on parameterization of G.A. Koppers, and          =*
!=  D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                         =*
!=  Final values do include effects from the Herzberg continuum.             =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
!=            altitude                                                       =*
!=  TLEV    - tmeperature at each level                                   (I)=*
!=  SECCHI  - ratio of slant to vertical o2 columns                       (I)=*
!=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
!=            vertical layer at each specified wavelength                    =*
!=  O2XSK  - REAL, molecular absorption cross section in SR bands at     (O)=*
!=            each specified wavelength.  Includes Herzberg continuum        =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nz
        real, intent(IN)    :: o2col(nz)
        real, intent(IN)    :: secchi(nz)
        real, intent(IN)    :: tLev(nz)

        real, intent(OUT)   :: dto2(nz, 17)
        real, intent(OUT)   :: o2xsk(nz, 17)

        real :: o2col1(nz)
        integer :: i, k
        integer :: ktop, ktop1, kbot
        real :: xs(17), x

!------------------------------------------
!sm  Initialize cross sections to values
!sm  at large optical depth
!------------------------------------------
        do i = 1, 17
            do k = 1, nz
                o2xsk(k, i) = xslod(i)
            end do
        end do

!------------------------------------------
!     Calculate cross sections
!sm:  Set smallest O2col = exp(38.) molec cm-2
!sm     to stay in range of parameterization
!sm     given by Koppers et al. at top of atm.
!------------------------------------------
        ktop = 121
        kbot = 0

        do k = 1, nz    !! loop for alt
            o2col1(k) = max(o2col(k), exp(38.))
            x = ALOG(o2col1(k))
            if (x < 38.0) then
                ktop1 = k - 1
                ktop = min(ktop1, ktop)
            else if (x > 56.0) then
                kbot = k
            else
                call effxs(x, tLev(k), xs)
                do i = 1, 17
                    o2xsk(k, i) = xs(i)
                end do
            end if
        end do                    !! finish loop for alt

!------------------------------------------
!  fill in cross section where X is out of range
!  by repeating edge table values
!------------------------------------------
!sm do not allow kbot = nz to avoid division by zero in
!   no light case.
        if (kbot == nz) kbot = nz - 1
        do k = 1, kbot
            do i = 1, 17
                o2xsk(k, i) = o2xsk(kbot + 1, i)
            end do
        end do
        do k = ktop + 1, nz
            do i = 1, 17
                o2xsk(k, i) = o2xsk(ktop, i)
            end do
        end do
!------------------------------------------
!  Calculate incremental optical depths
!------------------------------------------
        do i = 1, 17                   ! loop over wavelength
            do k = 1, nz - 1            ! loop for alt
!... calculate an optical depth weighted by density
!sm:  put in mean value estimate, if in shade
                if (abs(1.-o2col1(k + 1)/o2col1(k)) <= 2.*precis) then
                    dto2(k, i) = o2xsk(k + 1, i)*o2col1(k + 1)/(nz - 1)
                else
                    dto2(k, i) = abs((o2xsk(k + 1, i)*o2col1(k + 1) - &
                                       o2xsk(k, i)*o2col1(k))/ &
                                      (1.0 + ALOG(o2xsk(k + 1, i)/o2xsk(k, i))/ &
                                       ALOG(o2col1(k + 1)/o2col1(k))))
!... change to vertical optical depth
                    dto2(k, i) = 2.*dto2(k, i)/(secchi(k) + secchi(k + 1))
                end if
            end do
        end do

        do i = 1, 17                   ! loop over wavelength
            dto2(nz, i) = 0.0       ! set optical depth to zero at top
        end do

    end subroutine schum

!=============================================================================*
    subroutine effxs(x, t, xs)

!     Subroutine for evaluating the effective cross section
!     of O2 in the Schumann-Runge bands using parameterization
!     of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
!     68-79, 1996]

!     method:
!     ln(xs) = A(X)[T-220]+B(X)
!     X = log of slant column of O2
!     A,B calculated from Chebyshev polynomial coeffs
!     AC and BC using NR routine chebev.  Assume interval
!     is 38<ln(NO2)<56.

!     Revision History:

!     drm 2/97  initial coding

!-------------------------------------------------------------
        implicit none

        real(real32), intent(IN OUT) :: x
        real(real32), intent(IN)     :: t
        real(real32), intent(OUT)    :: xs(17)

        real(real32) a(17), b(17)
        integer :: i

        call calc_params(x, a, b)

        do i = 1, 17
            xs(i) = exp(a(i)*(t - 220.) + b(i))
        end do

    end subroutine effxs

!=============================================================================*
    subroutine calc_params(x, a, b)

!-------------------------------------------------------------

!       calculates coefficients (A,B), used in calculating the
! effective cross section, for 17 wavelength intervals
!       as a function of log O2 column density (X)
!       Wavelength intervals are defined in WMO1985

!-------------------------------------------------------------
        implicit none

        real(real32), intent(IN OUT) :: x
        real(real32), intent(OUT)    :: a(17)
        real(real32), intent(OUT)    :: b(17)

        integer :: i

!       call Chebyshev Evaluation routine to calc A and B from
! set of 20 coeficients for each wavelength
        do i = 1, 17
            a(i) = chebev(38.0, 56.0, ac(1, i), 20, x)
            b(i) = chebev(38.0, 56.0, bc(1, i), 20, x)
        end do

    end subroutine calc_params

!=============================================================================*

    subroutine init_xs()
        implicit none

!       locals
        integer(int32) in_lun ! file unit number
        integer(int32) i, j

        in_lun = 11

        open (UNIT=in_lun, FILE=trim(files(138)%fileName), FORM='FORMATTED')

        read (in_lun, 901)
        do i = 1, 20
            read (in_lun, 903) (ac(i, j), j=1, 17)
        end do
        read (in_lun, 901)
        do i = 1, 20
            read (in_lun, 903) (bc(i, j), j=1, 17)
        end do

901     format(/)
903     format(17(e23.14, 1x))

998     close (in_lun)

        do i = 1, 17
            wave_num(18 - i) = 48250.+(500.*i)
        end do

    end subroutine init_xs

!=============================================================================*
    real(real32) function chebev(a, b, c, m, x)
!     Chebyshev evaluation algorithm
!     See Numerical recipes p193
!-------------------------------------------------------------
        implicit none
        integer, intent(IN) :: m
        real(real32), intent(IN)  :: a
        real(real32), intent(IN)  :: b
        real(real64), intent(IN)  :: c(m)
        real(real32), intent(IN)  :: x

        integer :: j
        real :: d, dd, sv, y, y2

        if ((x - a)*(x - b) > 0.) then
            write (6, *) 'X NOT IN RANGE IN CHEBEV', x
            chebev = 0.0
            return
        end if

        d = 0.
        dd = 0.
        y = (2.*x - a - b)/(b - a)
        y2 = 2.*y
        do j = m, 2, -1
            sv = d
            d = y2*d - dd + c(j)
            dd = sv
        end do
        chebev = y*d - dd + 0.5*c(1)

    end function chebev

    subroutine inter1(ng, xg, yg, n, x, y)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on single, discrete points, onto a discrete target  =*
!=  grid.                                                                    =*
!=  The original input data are given on single, discrete points of an       =*
!=  arbitrary grid and are being linearly interpolated onto a specified      =*
!=  discrete target grid.  A typical example would be the re-gridding of a   =*
!=  given data set for the vertical temperature profile to match the speci-  =*
!=  fied altitude grid.                                                      =*
!=  Some caution should be used near the end points of the grids.  If the    =*
!=  input data set does not span the range of the target grid, the remaining =*
!=  points will be set to zero, as extrapolation is not permitted.           =*
!=  If the input data does not encompass the target grid, use ADDPNT to      =*
!=  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG  - INTEGER, number of points in the target grid                    (I)=*
!=  XG  - REAL, target grid (e.g. altitude grid)                          (I)=*
!=  YG  - REAL, y-data re-gridded onto XG                                 (O)=*
!=  N   - INTEGER, number of points in the input data set                 (I)=*
!=  X   - REAL, grid on which input data are defined                      (I)=*
!=  Y   - REAL, input y-data                                              (I)=*
!-----------------------------------------------------------------------------*

        integer, intent(IN)                      :: ng
        integer, intent(IN)                      :: n
        real, intent(IN)                         :: xg(ng)
        real, intent(IN)                         :: x(n)
        real, intent(IN)                         :: y(n)
        real, intent(OUT)                        :: yg(ng)

! local:
        real :: slope
        integer :: jsave, i, j
!_______________________________________________________________________

        jsave = 1
        do i = 1, ng
            yg(i) = 0.
            j = jsave
10          continue
            if ((x(j) > xg(i)) .or. (xg(i) >= x(j + 1))) then
                j = j + 1
                if (j <= n - 1) GO TO 10
!        ---- end of loop 10 ----
            else
                slope = (y(j + 1) - y(j))/(x(j + 1) - x(j))
                yg(i) = y(j) + slope*(xg(i) - x(j))
                jsave = j
            end if
        end do
!_______________________________________________________________________

    end subroutine inter1

!=============================================================================*

    subroutine inter2(ng, xg, yg, n, x, y, ierr)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on single, discrete points onto a set of target     =*
!=  bins.                                                                    =*
!=  The original input data are given on single, discrete points of an       =*
!=  arbitrary grid and are being linearly interpolated onto a specified set  =*
!=  of target bins.  In general, this is the case for most of the weighting  =*
!=  functions (action spectra, molecular cross section, and quantum yield    =*
!=  data), which have to be matched onto the specified wavelength intervals. =*
!=  The average value in each target bin is found by averaging the trapezoi- =*
!=  dal area underneath the input data curve (constructed by linearly connec-=*
!=  ting the discrete input values).                                         =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data set does not span the range of the target grid, an error      =*
!=  message is printed and the execution is stopped, as extrapolation of the =*
!=  data is not permitted.                                                   =*
!=  If the input data does not encompass the target grid, use ADDPNT to      =*
!=  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
!=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!=        bin i (i = 1..NG-1)                                                =*
!=  N   - INTEGER, number of points in input grid                         (I)=*
!=  X   - REAL, grid on which input data are defined                      (I)=*
!=  Y   - REAL, input y-data                                              (I)=*
!-----------------------------------------------------------------------------*

        integer, intent(IN)                      :: ng
        integer, intent(IN)                      :: n
        real, intent(IN)                         :: xg(ng)
        real, intent(OUT)                        :: yg(ng)
        real, intent(IN)                         :: x(n)
        real, intent(IN)                         :: y(n)
        integer, intent(OUT)                     :: ierr

! local:
        real :: area, xgl, xgu
        real :: darea, slope
        real :: a1, a2, b1, b2
        integer :: ngintv
        integer :: i, k, jstart

!_______________________________________________________________________

        ierr = 0

!  test for correct ordering of data, by increasing value of x
        do i = 2, n
            if (x(i) <= x(i - 1)) then
                ierr = 1
                write (*, *) 'data not sorted'
                return
            end if
        end do

        do i = 2, ng
            if (xg(i) <= xg(i - 1)) then
                ierr = 2
                write (0, *) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
                return
            end if
        end do

! check for xg-values outside the x-range
        if ((x(1) > xg(1)) .or. (x(n) < xg(ng))) then
            write (0, *) '>>> ERROR (inter2) <<<  Data do not span '//'grid.  '
            write (0, *) '                        Use ADDPNT to '// &
                'expand data and re-run.'
            stop
        end if

!  find the integral of each grid interval and use this to
!  calculate the average y value for the interval
!  xgl and xgu are the lower and upper limits of the grid interval
        jstart = 1
        ngintv = ng - 1
        do i = 1, ngintv

! initalize:
            area = 0.0
            xgl = xg(i)
            xgu = xg(i + 1)

!  discard data before the first grid interval and after the
!  last grid interval
!  for internal grid intervals, start calculating area by interpolating
!  between the last point which lies in the previous interval and the
!  first point inside the current interval

            k = jstart
            if (k <= n - 1) then

!  if both points are before the first grid, go to the next point
30              continue
                if (x(k + 1) <= xgl) then
                    jstart = k - 1
                    k = k + 1
                    if (k <= n - 1) GO TO 30
                end if

!  if the last point is beyond the end of the grid, complete and go to the next
!  grid
40              continue
                if ((k <= n - 1) .and. (x(k) < xgu)) then

                    jstart = k - 1

! compute x-coordinates of increment

                    a1 = max(x(k), xgl)
                    a2 = min(x(k + 1), xgu)

!  if points coincide, contribution is zero

                    if (x(k + 1) == x(k)) then
                        darea = 0.e0
                    else
                        slope = (y(k + 1) - y(k))/(x(k + 1) - x(k))
                        b1 = y(k) + slope*(a1 - x(k))
                        b2 = y(k) + slope*(a2 - x(k))
                        darea = (a2 - a1)*(b2 + b1)/2.
                    end if

!  find the area under the trapezoid from a1 to a2

                    area = area + darea

! go to next point

                    k = k + 1
                    GO TO 40

                end if

            end if

!  calculate the average y after summing the areas in the interval
            yg(i) = area/(xgu - xgl)

        end do
!_______________________________________________________________________

    end subroutine inter2

!=============================================================================*

    subroutine inter3(ng, xg, yg, n, x, y, foldin)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on a set of bins onto a different set of target     =*
!=  bins.                                                                    =*
!=  The input data are given on a set of bins (representing the integral     =*
!=  of the input quantity over the range of each bin) and are being matched  =*
!=  onto another set of bins (target grid).  A typical example would be an   =*
!=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
!=  vals, that has to be matched onto the working wavelength grid.           =*
!=  The resulting area in a given bin of the target grid is calculated by    =*
!=  simply adding all fractional areas of the input data that cover that     =*
!=  particular target bin.                                                   =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data do not span the full range of the target grid, the area in    =*
!=  the "missing" bins will be assumed to be zero.  If the input data extend =*
!=  beyond the upper limit of the target grid, the user has the option to    =*
!=  integrate the "overhang" data and fold the remaining area back into the  =*
!=  last target bin.  Using this option is recommended when re-gridding      =*
!=  vertical profiles that directly affect the total optical depth of the    =*
!=  model atmosphere.                                                        =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
!=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
!=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
!=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
!=           y-value for bin i (i = 1..NG-1)                                 =*
!=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
!=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
!=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
!=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
!=           y-value for bin i (i = 1..N-1)                                  =*
!=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
!=           FoldIn = 0 -> No folding of "overhang" data                     =*
!=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
!=                         last target bin                                   =*
!-----------------------------------------------------------------------------*
        integer, intent(IN)  :: ng
        integer, intent(IN)  :: n
        integer, intent(IN)  :: foldin
        real, intent(IN)     :: xg(ng)
        real, intent(IN)     :: x(n)
        real, intent(IN)     :: y(n)
        real, intent(OUT)    :: yg(ng)

! local:
        real :: a1, a2, sum
        real :: tail
        integer :: jstart, i, j, k
!_______________________________________________________________________

! check whether flag given is legal
        if ((foldin /= 0) .and. (foldin /= 1)) then
            write (0, *) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
            write (0, *) '                        Must be 0 or 1'
            stop
        end if

! do interpolation

        jstart = 1

        do i = 1, ng - 1

            yg(i) = 0.
            sum = 0.
            j = jstart

            if (j <= n - 1) then

20              continue

                if (x(j + 1) < xg(i)) then
                    jstart = j
                    j = j + 1
                    if (j <= n - 1) GO TO 20
                end if

25              continue

                if ((x(j) <= xg(i + 1)) .and. (j <= n - 1)) then

                    a1 = AMAX1(x(j), xg(i))
                    a2 = AMIN1(x(j + 1), xg(i + 1))

                    sum = sum + y(j)*(a2 - a1)/(x(j + 1) - x(j))
                    j = j + 1
                    GO TO 25

                end if

                yg(i) = sum

            end if

        end do

! if wanted, integrate data "overhang" and fold back into last bin

        if (foldin == 1) then

            j = j - 1
            a1 = xg(ng)     ! upper limit of last interpolated bin
            a2 = x(j + 1)     ! upper limit of last input bin considered

!        do folding only if grids don't match up and there is more input
            if ((a2 > a1) .or. (j + 1 < n)) then
                tail = y(j)*(a2 - a1)/(x(j + 1) - x(j))
                do k = j + 1, n - 1
                    tail = tail + y(k)*(x(k + 1) - x(k))
                end do
                yg(ng - 1) = yg(ng - 1) + tail
            end if

        end if
!_______________________________________________________________________

    end subroutine inter3

!=============================================================================*

    subroutine inter4(ng, xg, yg, n, x, y, foldin)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on a set of bins onto a different set of target     =*
!=  bins.                                                                    =*
!=  The input data are given on a set of bins (representing the integral     =*
!=  of the input quantity over the range of each bin) and are being matched  =*
!=  onto another set of bins (target grid).  A typical example would be an   =*
!=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
!=  vals, that has to be matched onto the working wavelength grid.           =*
!=  The resulting area in a given bin of the target grid is calculated by    =*
!=  simply adding all fractional areas of the input data that cover that     =*
!=  particular target bin.                                                   =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data do not span the full range of the target grid, the area in    =*
!=  the "missing" bins will be assumed to be zero.  If the input data extend =*
!=  beyond the upper limit of the target grid, the user has the option to    =*
!=  integrate the "overhang" data and fold the remaining area back into the  =*
!=  last target bin.  Using this option is recommended when re-gridding      =*
!=  vertical profiles that directly affect the total optical depth of the    =*
!=  model atmosphere.                                                        =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
!=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
!=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
!=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
!=           y-value for bin i (i = 1..NG-1)                                 =*
!=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
!=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
!=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
!=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
!=           y-value for bin i (i = 1..N-1)                                  =*
!=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
!=           FoldIn = 0 -> No folding of "overhang" data                     =*
!=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
!=                         last target bin                                   =*
!-----------------------------------------------------------------------------*
        integer, intent(IN) :: ng
        integer, intent(IN) :: n
        integer, intent(IN) :: foldin
        real, intent(IN)    :: xg(ng)
        real, intent(OUT)   :: yg(ng)
        real, intent(IN)    :: x(n)
        real, intent(IN)    :: y(n)

! local:
        real :: a1, a2, sum
        real :: tail
        integer :: jstart, i, j, k
!_______________________________________________________________________

! check whether flag given is legal
        if ((foldin /= 0) .and. (foldin /= 1)) then
            write (0, *) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
            write (0, *) '                        Must be 0 or 1'
            stop
        end if

! do interpolation

        jstart = 1

        do i = 1, ng - 1

            yg(i) = 0.
            sum = 0.
            j = jstart

            if (j <= n - 1) then

20              continue

                if (x(j + 1) < xg(i)) then
                    jstart = j
                    j = j + 1
                    if (j <= n - 1) GO TO 20
                end if

25              continue

                if ((x(j) <= xg(i + 1)) .and. (j <= n - 1)) then

                    a1 = AMAX1(x(j), xg(i))
                    a2 = AMIN1(x(j + 1), xg(i + 1))

                    sum = sum + y(j)*(a2 - a1)

                    j = j + 1
                    GO TO 25

                end if

                yg(i) = sum/(xg(i + 1) - xg(i))

            end if

        end do

! if wanted, integrate data "overhang" and fold back into last bin

        if (foldin == 1) then

            j = j - 1
            a1 = xg(ng)     ! upper limit of last interpolated bin
            a2 = x(j + 1)     ! upper limit of last input bin considered

!        do folding only if grids don't match up and there is more input
            if ((a2 > a1) .or. (j + 1 < n)) then
                tail = y(j)*(a2 - a1)/(x(j + 1) - x(j))
                do k = j + 1, n - 1
                    tail = tail + y(k)*(x(k + 1) - x(k))
                end do
                yg(ng - 1) = yg(ng - 1) + tail
            end if

        end if
!_______________________________________________________________________

    end subroutine inter4

!=============================================================================*

    subroutine addpnt(x, y, ld, n, xnew, ynew)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
!=  ascending order                                                          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
!=  Y    - REAL vector of length LD, y-values                            (IO)=*
!=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
!=         program                                                           =*
!=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
!=         N < LD.  On exit, N is incremented by 1.                          =*
!=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
!=  YNEW - REAL, y-value of point to be added                             (I)=*
!-----------------------------------------------------------------------------*
        integer, intent(IN)     :: ld
        integer, intent(IN OUT) :: n
        real, intent(IN OUT)    :: x(ld)
        real, intent(OUT)       :: y(ld)
        real, intent(IN)        :: xnew
        real, intent(IN)        :: ynew

! local variables

        integer :: insert
        integer :: i

!-----------------------------------------------------------------------

! check n<ld to make sure x will hold another point

        if (n >= ld) then
            write (0, *) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
            write (0, *) '                        All elements used.'
            stop
        end if

        insert = 1
        i = 2

! check, whether x is already sorted.
! also, use this loop to find the point at which xnew needs to be inserted
! into vector x, if x is sorted.

10      continue
        if (i < n) then
            if (x(i) < x(i - 1)) then
                write (0, *) '>>> ERROR (ADDPNT) <<<  x-data must be '// &
                    'in ascending order!'
                stop
            else
                if (xnew > x(i)) insert = i + 1
            end if
            i = i + 1
            GO TO 10
        end if

! if <xnew,ynew> needs to be appended at the end, just do so,
! otherwise, insert <xnew,ynew> at position INSERT

        if (xnew > x(n)) then

            x(n + 1) = xnew
            y(n + 1) = ynew

        else

! shift all existing points one index up

            do i = n, insert, -1
                x(i + 1) = x(i)
                y(i + 1) = y(i)
            end do

! insert new point

            x(insert) = xnew
            y(insert) = ynew

        end if

! increase total number of elements in x, y

        n = n + 1

    end subroutine addpnt

!=============================================================================*

    subroutine zero1(x, m)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Initialize all elements of a floating point vector with zero.            =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  X  - REAL, vector to be initialized                                   (O)=*
!=  M  - INTEGER, number of elements in X                                 (I)=*
!-----------------------------------------------------------------------------*
        integer, intent(IN)                      :: m
        real, intent(OUT)                        :: x(m)
        integer :: i

        do i = 1, m
            x(i) = 0.
        end do

    end subroutine zero1

!=============================================================================*

    subroutine zero2(x, m, n)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Initialize all elements of a 2D floating point array with zero.          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  X  - REAL, array to be initialized                                    (O)=*
!=  M  - INTEGER, number of elements along the first dimension of X,      (I)=*
!=       exactly as specified in the calling program                         =*
!=  N  - INTEGER, number of elements along the second dimension of X,     (I)=*
!=       exactly as specified in the calling program                         =*
!-----------------------------------------------------------------------------*

        integer, intent(IN)                      :: m
        integer, intent(IN)                      :: n
        real, intent(OUT)                        :: x(m, n)

! m,n : dimensions of x, exactly as specified in the calling program

        integer :: i, j

        do j = 1, n
            do i = 1, m
                x(i, j) = 0.0
            end do
        end do

    end subroutine zero2

    subroutine odrl(nz, nw, wl, cAir, dtRl)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute Rayleigh optical depths as a function of altitude and wavelength =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  Z       - REAL, specified altitude working grid (km)                  (I)=*
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
!=            wavelength grid                                                =*
!=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
!=            working wavelength grid                                        =*
!=  C       - REAL, number of air molecules per cm^2 at each specified    (O)=*
!=            altitude layer                                                 =*
!=  DTRL    - REAL, Rayleigh optical depth at each specified altitude     (O)=*
!=            and each specified wavelength                                  =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nz
        integer, intent(IN) :: nw
        real, intent(IN)    :: wl(kw)
        real, intent(IN)    :: cAir(nz)
        real, intent(OUT)   :: dtRl(nz, kw)
        real, parameter      :: oneover = 1.0/1.0e3

        real :: srayl, wc, wmicrn, xx
        integer :: iz, iw

!_______________________________________________________________________

! compute Rayleigh cross sections and depths:

        do iw = 1, nw - 1
            wc = (wl(iw) + wl(iw + 1))/2.

! Rayleigh scattering cross section from WMO 1985 (originally from
! Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
! An empirical formula for its calculation in the homoshpere, Planet.
! Space Sci., 32, 1467-1468, 1984.
            xx = 4.04
            wmicrn = wc*oneover
            if (wmicrn <= 0.55) xx = 3.6772 + 0.389*wmicrn + 0.09426/wmicrn
            srayl = 4.02e-28/(wmicrn)**xx
! alternate (older) expression from
! Frohlich and Shaw, Appl.Opt. v.11, p.1773 (1980).
!     xx = 3.916 + 0.074*wmicrn + 0.050/wmicrn
!     srayl(iw) = 3.90e-28/(wmicrn)**xx
            do iz = 1, nz
                if (iz > (nz - 1)) cycle
                dtRl(iz, iw) = cAir(iz)*srayl
            end do

        end do
!_______________________________________________________________________

    end subroutine odrl

! subroutines used for calculation of quantum yields for
! various photoreactions:
!     qyacet - q.y. for acetone, based on Blitz et al. (2004)

!*******************************************************************************

    subroutine qyacet(w, t, m, fco, fac)
! Compute acetone quantum yields according to the parameterization of:
! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield
!       (2004), Pressure and temperature-dependent quantum yields for the
!       photodissociation of acetone between 279 and 327.5 nm, Geophys.
!       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.
! input:
! w = wavelength, nm
! T = temperature, K
! m = air number density, molec. cm-3
! output
! fco = quantum yield for product CO
! fac = quantum yield for product CH3CO (acetyl radical)
        implicit none

        real, intent(IN)  :: w
        real, intent(IN)  :: t
        real, intent(IN)  :: m
        real, intent(OUT) :: fco
        real, intent(OUT) :: fac

        real :: a0, a1, a2, a3, a4
        real :: b0, b1, b2, b3, b4
        real :: c3
        real :: ca0, ca1, ca2, ca3, ca4

!** set out-of-range values:
! use low pressure limits for shorter wavelengths
! set to zero beyound 327.5

        if (w < 279.) then
            fco = 0.05
            fac = 0.95
            return
        end if

        if (w > 327.5) then
            fco = 0.
            fac = 0.
            return
        end if

!** CO (carbon monoxide) quantum yields:

        a0 = 0.350*(t/295.)**(-1.28)
        b0 = 0.068*(t/295.)**(-2.65)
        ca0 = exp(b0*(w - 248.))*a0/(1.-a0)

        fco = 1./(1 + ca0)

!** CH3CO (acetyl radical) quantum yields:

        if (w >= 279. .and. w < 302.) then

            a1 = 1.600e-19*(t/295.)**(-2.38)
            b1 = 0.55e-3*(t/295.)**(-3.19)
            ca1 = a1*exp(-b1*((1.e7/w) - 33113.))

            fac = (1.-fco)/(1 + ca1*m)

        end if

        if (w >= 302. .and. w < 327.5) then

            a2 = 1.62e-17*(t/295.)**(-10.03)
            b2 = 1.79e-3*(t/295.)**(-1.364)
            ca2 = a2*exp(-b2*((1.e7/w) - 30488.))

            a3 = 26.29*(t/295.)**(-6.59)
            b3 = 5.72e-7*(t/295.)**(-2.93)
            c3 = 30006*(t/295.)**(-0.064)
            ca3 = a3*exp(-b3*((1.e7/w) - c3)**2)

            a4 = 1.67e-15*(t/295.)**(-7.25)
            b4 = 2.08e-3*(t/295.)**(-1.16)
            ca4 = a4*exp(-b4*((1.e7/w) - 30488.))

            fac = (1.-fco)*(1.+ca3 + ca4*m)/ &
                  ((1.+ca3 + ca2*m)*(1.+ca4*m))

        end if

    end subroutine qyacet

    subroutine rdetfl(nw, wl)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and re-grid extra-terrestrial flux data.                            =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nw    ! input: (wavelength grid)
        real, intent(IN)    :: wl(kw)

        integer, parameter :: kdata = 20000
        integer :: iw

! work arrays for input data files:

        character(LEN=200) :: fil
        real :: x1(kdata)
        real :: y1(kdata)
        integer :: nhead, n, i, ierr
        real :: dum

! data gridded onto wl(kw) grid:

        real :: yg1(kw)
        real :: yg2(kw)
        real :: yg3(kw)

        real, parameter :: hc = 6.62e-34*2.998e8

! simple files are read and interpolated here in-line. Reading of
! more complex files may be done with longer code in a read#.f subroutine.

        select case (msun)
        case (1)
            fil = trim(files(2)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 3
            n = 121
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (2)
            fil = trim(files(3)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 3
            n = 4327
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (3)
            fil = trim(files(4)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 6
            n = 14980
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (4)
            fil = trim(files(5)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 8
            n = 1260
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (5)
! unofficial - do not use
            fil = trim(files(6)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 11
            n = 2047
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (6)
! unofficial - do not use
            fil = trim(files(7)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 3
            n = 1200
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-3
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (7)
            fil = trim(files(8)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 11
            n = 496
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) dum, y1(i)
                if (dum < 630.0) x1(i) = dum - 0.5
                if (dum > 630.0 .and. dum < 870.0) x1(i) = dum - 1.0
                if (dum > 870.0) x1(i) = dum - 2.5
                y1(i) = y1(i)*1.e4*hc/(dum*1.e-9)
            end do
            close (kin)
            x1(n + 1) = x1(n) + 2.5
            do i = 1, n
                y1(i) = y1(i)*(x1(i + 1) - x1(i))
            end do
            call inter3(nw, wl, yg2, n + 1, x1, y1, 0)
            do iw = 1, nw - 1
                yg1(iw) = yg1(iw)/(wl(iw + 1) - wl(iw))
            end do
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (8)
            nhead = 5
            fil = trim(files(9)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 13
            n = 5160
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-3
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (9)
            fil = trim(files(10)%filename)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 2
            n = 302
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (10)
!!WRITE(kout,*) 'datae1/sun/susim_hi.flx'
            call read1(nw, wl, yg1)
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (11)
!!WRITE(kout,*) 'datae1/sun/wmo85.flx'
            call read2(nw, wl, yg1)
            do iw = 1, nw - 1
                f(iw) = yg1(iw)
            end do
        case (12)
!!WRITE(kout,*) 'datae1/sun/susim_hi.flx'
            call read1(nw, wl, yg1)
            fil = trim(files(13)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 11
            n = 496
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) dum, y1(i)
                if (dum < 630.0) x1(i) = dum - 0.5
                if (dum > 630.0 .and. dum < 870.0) x1(i) = dum - 1.0
                if (dum > 870.0) x1(i) = dum - 2.5
                y1(i) = y1(i)*1.e4*hc/(dum*1.e-9)
            end do
            close (kin)
            x1(n + 1) = x1(n) + 2.5
            do i = 1, n
                y1(i) = y1(i)*(x1(i + 1) - x1(i))
            end do
            call inter3(nw, wl, yg2, n + 1, x1, y1, 0)
            do iw = 1, nw - 1
                yg2(iw) = yg2(iw)/(wl(iw + 1) - wl(iw))
            end do
            do iw = 1, nw - 1
                if (wl(iw) > 350.) then
                    f(iw) = yg2(iw)
                else
                    f(iw) = yg1(iw)
                end if
            end do
        case (13)
!!WRITE(kout,*) 'datae1/sun/susim_hi.flx'
            call read1(nw, wl, yg1)
            nhead = 5
            fil = trim(files(14)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            n = 5160
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-3
            end do
            close (kin)
            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg2, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, fil
                stop
            end if
            fil = trim(files(15)%fileName)
!WRITE(kout,*) fil
            open (UNIT=kin, FILE=trim(fil), STATUS='old')
            nhead = 11
            n = 496
            do i = 1, nhead
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) dum, y1(i)
                if (dum < 630.0) x1(i) = dum - 0.5
                if (dum > 630.0 .and. dum < 870.0) x1(i) = dum - 1.0
                if (dum > 870.0) x1(i) = dum - 2.5
                y1(i) = y1(i)*1.e4*hc/(dum*1.e-9)
            end do
            close (kin)

            x1(n + 1) = x1(n) + 2.5
            call inter4(nw, wl, yg3, n + 1, x1, y1, 0)

            do iw = 1, nw - 1
                if (wl(iw) < 150.01) then
                    f(iw) = yg1(iw)
                else if ((wl(iw) >= 150.01) .and. wl(iw) <= 400.) then
                    f(iw) = yg2(iw)
                else if (wl(iw) > 400.) then
                    f(iw) = yg3(iw)
                end if
            end do
        end select

    end subroutine rdetfl

    subroutine read1(nw, wl, f)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
!=  working wavelength grid.                                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nw    ! input: (wavelength grid)
        real, intent(IN)    :: wl(kw)
        real, intent(OUT)   :: f(kw) ! output: (extra terrestrial solar flux)

! local:

        real :: lambda_hi(10000), irrad_hi(10000)
        real :: lambda
        integer :: ierr
        integer :: i, j, n
        character(LEN=200) :: fil

!_______________________________________________________________________

!****** SUSIM irradiance
!_______________________________________________________________________
! VanHoosier, M. E., J.-D. F. Bartoe, G. E. Brueckner, and
! D. K. Prinz, Absolute solar spectral irradiance 120 nm -
! 400 nm (Results from the Solar Ultraviolet Spectral Irradiance
! Monitor - SUSIM- Experiment on board Spacelab 2),
! Astro. Lett. and Communications, 1988, vol. 27, pp. 163-168.
!     SUSIM SL2 high resolution (0.15nm) Solar Irridance data.
!     Irradiance values are given in milliwatts/m^2/nanomenters
!     and are listed at 0.05nm intervals.  The wavelength given is
!     the center wavelength of the 0.15nm triangular bandpass.
!     Normalized to 1 astronomical unit.
!  DATA for wavelengths > 350 nm are unreliable
! (Van Hoosier, personal communication, 1994).
!_______________________________________________________________________

!* high resolution

        fil = trim(files(11)%fileName)
        open (UNIT=kin, FILE=trim(fil), STATUS='old')
        do i = 1, 7
            read (kin, *)
        end do
        do i = 1, 559
            read (kin, *) lambda, (irrad_hi(10*(i - 1) + j), j=1, 10)
        end do
        close (kin)

! compute wavelengths, convert from mW to W

        n = 559*10
        do i = 1, n
            lambda_hi(i) = 120.5 + FLOAT(i - 1)*.05
            irrad_hi(i) = irrad_hi(i)/1000.
        end do
!_______________________________________________________________________

        call addpnt(lambda_hi, irrad_hi, 10000, n, lambda_hi(1)*(1.-deltax), 0.)
        call addpnt(lambda_hi, irrad_hi, 10000, n, 0., 0.)
        call addpnt(lambda_hi, irrad_hi, 10000, n, lambda_hi(n)*(1.+deltax), 0.)
        call addpnt(lambda_hi, irrad_hi, 10000, n, 1.e38, 0.)
        call inter2(nw, wl, f, n, lambda_hi, irrad_hi, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, fil
            stop
        end if

    end subroutine read1

!=============================================================================*

    subroutine read2(nw, wl, f)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
!=  working wavelength grid.                                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nw     ! input: (wavelength grid)
        real, intent(IN)    :: wl(kw)
        real, intent(OUT)   :: f(kw)  ! output: (extra terrestrial solar flux)

        real :: yg(kw)
        integer :: iw

! local:

        real :: x1(1000), y1(1000)
        real :: x2(1000)
        real :: x3(1000)
        integer :: i, n
        real :: dum
        integer :: idum

!_______________________________________________________________________

!********WMO 85 irradiance

        open (UNIT=kin, FILE=trim(files(12)%fileName), STATUS='old')
        do i = 1, 3
            read (kin, *)
        end do
        n = 158
        do i = 1, n
            read (kin, *) idum, x1(i), x2(i), y1(i), dum, dum, dum
            x3(i) = 0.5*(x1(i) + x2(i))

! average value needs to be calculated only if inter2 is
! used to interpolate onto wavelength grid (see below)
!        y1(i) =  y1(i) / (x2(i) - x1(i))

        end do
        close (kin)

        x1(n + 1) = x2(n)

! inter2: INPUT : average value in each bin
!         OUTPUT: average value in each bin
! inter3: INPUT : total area in each bin
!         OUTPUT: total area in each bin

        call inter3(nw, wl, yg, n + 1, x1, y1, 0)
!      CALL inter2(nw,wl,yg,n,x3,y1,ierr)

        do iw = 1, nw - 1
! from quanta s-1 cm-2 bin-1 to  watts m-2 nm-1
! 1.e4 * ([hc =] 6.62E-34 * 2.998E8)/(wc*1e-9)

! the scaling by bin width needs to be done only if
! inter3 is used for interpolation

            yg(iw) = yg(iw)/(wl(iw + 1) - wl(iw))
            f(iw) = yg(iw)*1.e4*(6.62e-34*2.998e8)/ &
                    (0.5*(wl(iw + 1) + wl(iw))*1.e-9)

        end do

    end subroutine read2

    subroutine rdno2xs(nw, wl, no2xs)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read NO2 molecular absorption cross section.  Re-grid data to match      =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  NO2XS  - REAL, molecular absoprtion cross section (cm^2) of NO2 at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nw
        real, intent(IN)    :: wl(kw)
        real, intent(OUT)   :: no2xs(kw)

        integer, parameter :: kdata = 1000

! input: (altitude working grid)

! local:
        real :: x1(kdata)
        real :: y1(kdata)
        real :: yg(kw)
        real :: dum
        integer :: ierr
        integer :: i, l, n, idum
        character(LEN=200) :: fil
!_______________________________________________________________________

!************ absorption cross sections:
!     measurements by:
! Davidson, J. A., C. A. Cantrell, A. H. McDaniel, R. E. Shetter,
! S. Madronich, and J. G. Calvert, Visible-ultraviolet absorption
! cross sections for NO2 as a function of temperature, J. Geophys.
! Res., 93, 7105-7112, 1988.
!  Values at 273K from 263.8 to 648.8 nm in approximately 0.5 nm intervals

        fil = trim(files(19)%fileName)

        open (UNIT=kin, FILE=trim(fil), STATUS='old')
        n = 750
        do i = 1, n
            read (kin, *) x1(i), y1(i), dum, dum, idum
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, fil
            stop
        end if

        do l = 1, nw - 1
            no2xs(l) = yg(l)
        end do

    end subroutine rdno2xs

    subroutine rdo2xs(nw, wl, o2xs1)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute equivalent O2 cross section, except                              =*
!=  the SR bands and the Lyman-alpha line.                                   =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
!=            wavelength grid                                                =*
!=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
!=            working wavelength grid
!=            vertical layer at each specified wavelength                    =*
!=  O2XS1   - REAL, O2 molecular absorption cross section                    =*
!=
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nw
        real, intent(IN)    :: wl(kw)
        real, intent(OUT)   :: o2xs1(kw) ! Output O2 xsect, temporary, will be
! over-written in Lyman-alpha and
! Schumann-Runge wavelength bands.

        integer :: i, n
        integer, parameter :: kdata = 200
        real :: x1(kdata), y1(kdata)
        real :: x, y
        integer :: ierr

!-----------------------------------------------------

! Read O2 absorption cross section data:
!  116.65 to 203.05 nm = from Brasseur and Solomon 1986
!  205 to 240 nm = Yoshino et al. 1988

! Note that subroutine la_srb.f will over-write values in the spectral regions
!   corresponding to:
! - Lyman-alpha (LA: 121.4-121.9 nm, Chabrillat and Kockarts parameterization)
! - Schumann-Runge bands (SRB: 174.4-205.8 nm, Koppers parameteriaztion)

        n = 0

        open (UNIT=kin, FILE=trim(files(16)%fileName))
        do i = 1, 7
            read (kin, *)
        end do
        do i = 1, 78
            read (kin, *) x, y
            if (x <= 204.) then
                n = n + 1
                x1(n) = x
                y1(n) = y
            end if
        end do
        close (kin)

        open (UNIT=kin, FILE=trim(files(17)%fileName), STATUS='old')
        do i = 1, 8
            read (kin, *)
        end do
        do i = 1, 36
            n = n + 1
            read (kin, *) x, y
            y1(n) = y*1.e-24
            x1(n) = x
        end do
        close (kin)

! Add termination points and interpolate onto the
!  user grid (set in subroutine gridw):

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, o2xs1, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'O2 -> O + O'
            stop
        end if

    end subroutine rdo2xs

    subroutine rdo3xs(nw, wl, nz, tLay, o3xs)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read ozone molecular absorption cross section.  Re-grid data to match    =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  O3XS   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
!=           each specified wavelength (WMO value at 273)                    =*
!=  S226   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
!=           each specified wavelength (value from Molina and Molina at 226K)=*
!=  S263   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
!=           each specified wavelength (value from Molina and Molina at 263K)=*
!=  S298   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
!=           each specified wavelength (value from Molina and Molina at 298K)=*
!=  opt    - if opt=0 read files
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nw
        integer, intent(IN) :: nz
        real, intent(IN)    :: wl(kw)
        real, intent(IN)    :: tLay(nz)
        real, intent(OUT)   :: o3xs(nz, kw)

        select case (mOption(1))
        case (1)
            call o3xs_mm(nw, wl, nz, o3xs, tLay)
        case (2)
            call o3xs_mal(nw, wl, nz, o3xs, tLay)
        case (3)
            call o3xs_bass(nw, wl, nz, o3xs, tLay)
        end select

    end subroutine rdo3xs

!=============================================================================*

    subroutine o3xs_mm(nw, wl, nz, xs, temp)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Interpolate the O3 cross section                                         =*
!=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
!=  175.439-847.5 nm) and:                                                   =*
!=  For Hartley and Huggins bands, use temperature-dependent values from     =*
!=  Molina, L. T., and M. J. Molina, Absolute absorption cross sections      =*
!=  of ozone in the 185- to 350-nm wavelength range, J. Geophys. Res.,       =*
!=  vol. 91, 14501-14508, 1986.                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
!=           at each defined wavelength and each defined altitude level      =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN)  :: nw
        integer, intent(IN)  :: nz
        real, intent(IN)     :: wl(kw)
        real, intent(IN)     :: temp(nz)
        real, intent(OUT)    :: xs(nz, kw)

        integer :: iw
        integer :: iz
        real, parameter :: div1 = 1.0/(263.-226.)
        real, parameter :: div2 = 1.0/(298.-263.)

        do iw = 1, nw - 1
            do iz = 1, nz
                xs(iz, iw) = mm_o3xs(iw)
                if (wl(iw) > 240.5 .and. wl(iw + 1) < 350.) then
                    if (temp(iz) < 263.) then
                        xs(iz, iw) = s226(iw) + (s263(iw) - s226(iw))* &
                                     (temp(iz) - 226.)*div1
                    else
                        xs(iz, iw) = s263(iw) + (s298(iw) - s263(iw))* &
                                     (temp(iz) - 263.)*div2
                    end if
                end if
            end do
        end do

    end subroutine o3xs_mm

!=============================================================================*

    subroutine o3xs_mal(nw, wl, nz, xs, temp)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and interpolate the O3 cross section                                =*
!=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
!=  175.439-847.5 nm) and:                                                   =*
!=  For Hartley and Huggins bands, use temperature-dependent values from     =*
!=  Malicet et al., J. Atmos. Chem.  v.21, pp.263-273, 1995.                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
!=           at each defined wavelength and each defined altitude level      =*
!-----------------------------------------------------------------------------*

        implicit none

        integer, intent(IN) :: nw
        integer, intent(IN) :: nz
        real, intent(IN)    :: wl(kw)
        real, intent(IN)     :: temp(nz)
        real, intent(OUT)   :: xs(nz, kw)

        integer :: iw
        integer :: iz
        integer, parameter :: kdata = 16000

! assign:
        do iw = 1, nw - 1
            do iz = 1, nz

                xs(iz, iw) = mm_o3xs(iw)
                if (wl(iw) > 195. .and. wl(iw + 1) < 345.) then
                    if (temp(iz) >= 243.) then
                        xs(iz, iw) = s243(iw) + (s295(iw) - s243(iw))* &
                                     (temp(iz) - 243.)/(295.-243.)
                    end if
                    if (temp(iz) < 254. .and. &
                        temp(iz) >= 228.) then
                        xs(iz, iw) = s228(iw) + (s243(iw) - s228(iw))* &
                                     (temp(iz) - 228.)/(243.-228.)
                    end if
                    if (temp(iz) < 228.) then
                        xs(iz, iw) = s218(iw) + (s228(iw) - s218(iw))* &
                                     (temp(iz) - 218.)/(228.-218.)
                    end if
                end if

            end do

        end do

    end subroutine o3xs_mal

!=============================================================================*

    subroutine o3xs_bass(nw, wl, nz, xs, temp)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and interpolate the O3 cross section                                =*
!=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
!=  175.439-847.5 nm) and:                                                   =*
!=  For Hartley and Huggins bands, use temperature-dependent values from     =*
!=  Bass
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
!=           at each defined wavelength and each defined altitude level      =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN)                      :: nw
        integer, intent(IN)                      :: nz
        real, intent(IN)                         :: wl(kw)
        real, intent(IN)     :: temp(nz)
        real, intent(OUT)                        :: xs(nz, kw)

        integer :: iw
        integer :: iz
        integer, parameter :: kdata = 2000
        real :: tc

        do iw = 1, nw - 1
            do iz = 1, nz

                tc = temp(iz) - 273.15

                xs(iz, iw) = mm_o3xs(iw)
                if (wl(iw) > 245. .and. wl(iw + 1) < 341.) then

                    xs(iz, iw) = 1.e-20*(c0(iw) + c1(iw)*tc + c2(iw)*tc*tc)

                end if

            end do

        end do

    end subroutine o3xs_bass

!=============================================================================*

    subroutine rdso2xs(nw, wl, so2xs)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read SO2 molecular absorption cross section.  Re-grid data to match      =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  SO2XS  - REAL, molecular absoprtion cross section (cm^2) of SO2 at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  02/97  Changed offset for grid-end interpolation to relative number      =*
!=         (x * (1 +- deltax)                                                =*
!-----------------------------------------------------------------------------*
!= This program is free software;  you can redistribute it and/or modify     =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                 =*
!= To obtain a copy of the GNU General Public License, write to:             =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= To contact the authors, please mail to:                                   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu                                            =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nw
        real, intent(IN)    :: wl(kw)
        real, intent(OUT)   :: so2xs(kw)

        integer, parameter :: kdata = 1000

! local:
        real :: x1(kdata)
        real :: y1(kdata)
        real :: yg(kw)
        integer :: ierr
        integer :: i, l, n!, idum
        character(LEN=200) :: fil
!_______________________________________________________________________

!************ absorption cross sections:
! SO2 absorption cross sections from J. Quant. Spectrosc. Radiat. Transfer
! 37, 165-182, 1987, T. J. McGee and J. Burris Jr.
! Angstrom vs. cm2/molecule, value at 221 K

        fil = 'data/McGee87'
        open (UNIT=kin, FILE=trim(files(18)%fileName), STATUS='old')
        do i = 1, 3
            read (kin, *)
        end do
!      n = 681
        n = 704
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            x1(i) = x1(i)/10.
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, fil
            stop
        end if

        do l = 1, nw - 1
            so2xs(l) = yg(l)
        end do

!_______________________________________________________________________

    end subroutine rdso2xs

    subroutine rtlink(nstr, nz, iw, sza, albedo, dtCld, omCld, gCld, &
                      dtAer, omAer, gAer, dtO2, dtO3, eDir, eDn, eUp, &
                      fDir, fDn, fUp, dtRl, dtSo2, dtNo2, nid, dsdh)
        implicit none

        integer, intent(IN)  :: nstr
        integer, intent(IN)  :: nz
        integer, intent(IN)  :: nid(0:nz)
        integer, intent(IN)  :: iw
        real   , intent(IN)  :: sza
        real   , intent(IN)  :: albedo(kw)
        real   , intent(IN)  :: dtCld(nz, kw)
        real   , intent(IN)  :: omCld(nz, kw)
        real   , intent(IN)  :: gCld(nz, kw)
        real   , intent(IN)  :: dtAer(nz, kw)
        real   , intent(IN)  :: omAer(nz, kw)
        real   , intent(IN)  :: gAer(nz, kw)
        real   , intent(IN)  :: dtO2(nz, kw)
        real   , intent(IN)  :: dtO3(nz, kw)
        real   , intent(IN)  :: dtRl(nz, kw)
        real   , intent(IN)  :: dtSo2(nz, kw)
        real   , intent(IN)  :: dtNo2(nz, kw)
        real   , intent(IN)  :: dsdh(0:nz, nz)
        !
        real   , intent(OUT) :: eDir(nz)
        real   , intent(OUT) :: eDn(nz)
        real   , intent(OUT) :: eUp(nz)
        real   , intent(OUT) :: fDir(nz)
        real   , intent(OUT) :: fDn(nz)
        real   , intent(OUT) :: fUp(nz)

        real :: dt(nz), om(nz), g(nz)
        real :: dtabs, dtsct, dscld, dacld, dsaer, daaer
        integer :: i, iii
! specific two ps2str
        real, dimension(nz) :: ediri, edni, eupi
        real, dimension(nz) :: fdiri, fdni, fupi
        logical, parameter :: delta = .true.
!  specific to psndo:
        real :: pmcld, pmray, pmaer
        real :: om1
        integer :: istr
        integer, parameter :: maxcly = 151
        integer, parameter :: maxulv = 151
        integer, parameter :: maxumu = 32
        integer, parameter :: maxcmu = 32
        integer, parameter :: maxphi = 3
        integer :: nlyr, numu
        real :: dtauc(maxcly), pmom(0:maxcmu, maxcly), &
                ssalb(maxcly), umu(maxumu), cwt(maxumu)
        real :: umu0
        real :: rfldir(maxulv), rfldn(maxulv), flup(maxulv), &
                u0u(maxumu, maxulv), &
                uavgso(maxulv), uavgup(maxulv), uavgdn(maxulv), &
                sindir(maxulv), sinup(maxulv), sindn(maxulv)
!print *,'LFR-DBG-RT 00'
        do i = 1, nz
            fdir(i) = 0.0
            fup(i) = 0.0
            fdn(i) = 0.0
            edir(i) = 0.0
            eup(i) = 0.0
            edn(i) = 0.0
            om(i) = 0.0
        end do

        umu0 = cos(sza*rpd)
!print *,'LFR-DBG-RT 01'
        do i = 1, nz - 1
!print *, 'LFR-DBG-RT 01.01',i
            dscld = dtCld(i, iw)*omCld(i, iw)
            dacld = dtCld(i, iw)*(1.-omCld(i, iw))
            dsaer = dtAer(i, iw)*omAer(i, iw)
            daaer = dtAer(i, iw)*(1.-omAer(i, iw))
            dtsct = dtRl(i, iw) + dscld + dsaer
            dtabs = dtso2(i, iw) + &
                    dtO2(i, iw) + dtO3(i, iw) + &
                    dtno2(i, iw) + dacld + daaer
            dtabs = AMAX1(dtabs, 1./largest)
            dtsct = AMAX1(dtsct, 1./largest)
            ! invert z-coordinate:
            iii = nz - i
            dt(iii) = dtsct + dtabs
            om(iii) = dtsct/(dtsct + dtabs)
            if (dtsct == 1./largest) om(iii) = 1./largest
!print *, 'LFR-DBG-RT 01.04', i, size(g),om(iii)
            g(iii) = (gCld(i, iw)*dscld + &
                      gAer(i, iw)*dsaer)/dtsct
!print *, 'LFR-DBG-RT 01.02', nstr
            if (nstr < 2) cycle

            ! DISORD parameters
            om1 = AMIN1(om(iii), 1.-precis)
            ssalb(iii) = AMAX1(om1, precis)
            dtauc(iii) = AMAX1(dt(iii), precis)
!print *, 'LFR-DBG-RT 01.03'
            !  phase function - assume Henyey-Greenstein for cloud and aerosol
            !  and Rayleigh for molecular scattering
            pmom(0, iii) = 1.0
            do istr = 1, nstr
                pmcld = gCld(i, iw)**(istr)
                pmaer = gAer(i, iw)**(istr)
                if (istr == 2) then
                    pmray = 0.1
                else
                    pmray = 0.0
                end if
                pmom(istr, iii) = (pmcld*dscld + pmaer*dsaer + &
                                    pmray*dtRl(i, iw))/dtsct
            end do
        end do
!print *,'LFR-DBG-RT 02'
        ! call rt routine:
        if (nstr < 2) then
            call ps2str(nz, sza, &
                        albedo(iw), dt, om, g, &
                        dsdh, nid, &
                        delta, fdiri, fupi, fdni, ediri, eupi, edni)
        else
            nlyr = nz - 1
            call psndo(nz, dsdh, &
                       nid, &
                       nlyr, dtauc, ssalb, pmom, &
                       albedo(iw), nstr, numu, umu, cwt, umu0, &
                       maxcly, maxulv, maxumu, maxcmu, maxphi, rfldir, &
                       rfldn, flup, u0u, &
                       uavgso, uavgup, uavgdn, sindir, sinup, sindn)
        end if ! output (invert z-coordinate)
!print *,'LFR-DBG-RT 03'
        if (nstr < 2) then
            do i = 1, nz
                iii = nz - i + 1
                fdir(i) = fdiri(iii)
                fup(i) = fupi(iii)
                fdn(i) = fdni(iii)
                edir(i) = ediri(iii)
                eup(i) = eupi(iii)
                edn(i) = edni(iii)
            end do
        else
            do i = 1, nz
                iii = nz - i + 1
                edir(i) = rfldir(iii)
                edn(i) = rfldn(iii)
                eup(i) = flup(iii)
                fdir(i) = 4.*pi*uavgso(iii)
                fdn(i) = 4.*pi*uavgdn(iii)
                fup(i) = 4.*pi*uavgup(iii)
            end do
        end if
!print *,'LFR-DBG-RT 04'
    end subroutine rtlink

!=============================================================================*

    subroutine ps2str(nlevel, zen, rsfc, tauu, omu, gu, dsdh, nid, delta, &
                      fdr, fup, fdn, edr, eup, edn)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Solve two-stream equations for multiple layers.  The subroutine is based =*
        !=  on equations from:  Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.=*
        !=  It contains 9 two-stream methods to choose from.  A pseudo-spherical     =*
        !=  correction has also been added.                                          =*
        !-----------------------------------------------------------------------------*
        !=  PARAMETERS:                                                              =*
        !=  NLEVEL  - INTEGER, number of specified altitude levels in the working (I)=*
        !=            grid                                                           =*
        !=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
        !=  RSFC    - REAL, surface albedo at current wavelength                  (I)=*
        !=  TAUU    - REAL, unscaled optical depth of each layer                  (I)=*
        !=  OMU     - REAL, unscaled single scattering albedo of each layer       (I)=*
        !=  GU      - REAL, unscaled asymmetry parameter of each layer            (I)=*
        !=  DSDH    - REAL, slant path of direct beam through each layer crossed  (I)=*
        !=            when travelling from the top of the atmosphere to layer i;     =*
        !=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
        !=  NID     - INTEGER, number of layers crossed by the direct beam when   (I)=*
        !=            travelling from the top of the atmosphere to layer i;          =*
        !=            NID(i), i = 0..NZ-1                                            =*
        !=  DELTA   - LOGICAL, switch to use delta-scaling                        (I)=*
        !=            .TRUE. -> apply delta-scaling                                  =*
        !=            .FALSE.-> do not apply delta-scaling                           =*
        !=  FDR     - REAL, contribution of the direct component to the total     (O)=*
        !=            actinic flux at each altitude level                            =*
        !=  FUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
        !=            the total actinic flux at each altitude level                  =*
        !=  FDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
        !=            the total actinic flux at each altitude level                  =*
        !=  EDR     - REAL, contribution of the direct component to the total     (O)=*
        !=            spectral irradiance at each altitude level                     =*
        !=  EUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
        !=            the total spectral irradiance at each altitude level           =*
        !=  EDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
        !=            the total spectral irradiance at each altitude level           =*
        !-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nlevel
        integer, intent(IN) :: nid(0:nlevel)
        real, intent(IN)    :: zen
        real, intent(IN)    :: rsfc
        real, intent(IN)    :: tauu(nlevel)
        real, intent(IN)    :: omu(nlevel)
        real, intent(IN)    :: gu(nlevel)
        real, intent(IN)    :: dsdh(0:nlevel, nlevel)
        logical, intent(IN) :: delta
        real, intent(OUT)   :: fdr(nlevel)
        real, intent(OUT)   :: fup(nlevel)
        real, intent(OUT)   :: fdn(nlevel)
        real, intent(OUT)   :: edr(nlevel)
        real, intent(OUT)   :: eup(nlevel)
        real, intent(OUT)   :: edn(nlevel)

        real :: tausla(0:nlevel), tauc(0:nlevel)
        real :: mu2(0:nlevel), mu, sum

        ! internal coefficients and matrix
        real :: lam(nlevel), taun(nlevel), bgam(nlevel)
        real :: e1(nlevel), e2(nlevel), e3(nlevel), e4(nlevel)
        real :: cup(nlevel), cdn(nlevel), cuptn(nlevel), cdntn(nlevel)
        real :: mu1(nlevel)
        integer :: row
        real, allocatable :: a(:), b(:), d(:), e(:), y(:)

        real :: pifs, fdn0
        real :: gi(nlevel), omi(nlevel), tempg
        real :: f, g, om
        real :: gam1, gam2, gam3, gam4


        ! For calculations of Associated Legendre Polynomials for GAMA1,2,3,4
        ! in delta-function, modified quadrature, hemispheric constant,
        ! Hybrid modified Eddington-delta function metods, p633,Table1.
        ! W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
        ! W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440,
        ! uncomment the following two lines and the appropriate statements further
        ! down.
        !     REAL YLM0, YLM2, YLM4, YLM6, YLM8, YLM10, YLM12, YLMS, BETA0,
        !    >     BETA1, BETAn, amu1, subd

        real :: expon, expon0, expon1, divisr, temp, up, dn
        real :: ssfc
        integer :: nlayer, mrows, lev, nrows

        integer :: i, j

        ! Some additional program constants:
        real, parameter :: eps = 1.e-3
!________
!print *,'LFR-DBG PS2 01'
        nrows = 2*nlevel
        allocate(a(nrows), b(nrows), d(nrows), e(nrows), y(nrows))

        ! MU = cosine of solar zenith angle
        ! RSFC = surface albedo
        ! TAUU =  unscaled optical depth of each layer
        ! OMU  =  unscaled single scattering albedo
        ! GU   =  unscaled asymmetry factor
        ! KLEV = max dimension of number of layers in atmosphere
        ! NLAYER = number of layers in the atmosphere
        ! NLEVEL = nlayer + 1 = number of levels

        ! initial conditions:  pi*solar flux = 1;  diffuse incidence = 0
        pifs = 1.
        fdn0 = 0.

        nlayer = nlevel - 1

        mu = cos(zen*rpd)
!print *,'LFR-DBG PS2 02', nlayer, nlevel, largest
        !************* compute coefficients for each layer:
        ! GAM1 - GAM4 = 2-stream coefficients, different for different approximations
        ! EXPON0 = calculation of e when TAU is zero
        ! EXPON1 = calculation of e when TAU is TAUN
        ! CUP and CDN = calculation when TAU is zero
        ! CUPTN and CDNTN = calc. when TAU is TAUN
        ! DIVISR = prevents division by zero
        do j = 0, nlevel
            tauc(j) = 0.
            tausla(j) = 0.
            mu2(j) = 1./sqrt(largest)
        end do
!print *,'LFR-DBG PS2 03', delta
        if (.not. delta) then
            do i = 1, nlayer
                gi(i) = gu(i)
                omi(i) = omu(i)
                taun(i) = tauu(i)
            end do
        else

            ! delta-scaling. Have to be done for delta-Eddington approximation,
            ! delta discrete ordinate, Practical Improved Flux Method, delta function,
            ! and Hybrid modified Eddington-delta function methods approximations
            do i = 1, nlayer
!print *,'LFR-DBG PS2 04',gu(i),omu(i)
                f = gu(i)*gu(i)
                gi(i) = (gu(i) - f)/(1 - f)
!print *,'LFR-DBG PS2 04.01',f,gi(i),(1 - f)*omu(i),(1 - omu(i)*f)
                omi(i) = (1 - f)*omu(i)/(1 - omu(i)*f)
                taun(i) = (1 - omu(i)*f)*tauu(i)
            end do
        end if
!print *,'LFR-DBG PS2 05'
        ! calculate slant optical depth at the top of the atmosphere when zen>90.
        ! in this case, higher altitude of the top layer is recommended which can
        ! be easily changed in gridz.f.
        if (zen > 90.0) then
            if (nid(0) < 0) then
                tausla(0) = largest
            else
                sum = 0.0
                do j = 1, nid(0)
                    sum = sum + 2.*taun(j)*dsdh(0, j)
                end do
                tausla(0) = sum
            end if
        end if
!print *,'LFR-DBG PS2 06', zen, nid(0), tausla(0)
        do i = 1, nlayer
            g = gi(i)
            om = omi(i)
            tauc(i) = tauc(i - 1) + taun(i)
            ! stay away from 1 by precision.  For g, also stay away from -1
            tempg = AMIN1(abs(g), 1.-precis)
            g = sign(tempg, g)
            om = AMIN1(om, 1.-precis)

            ! calculate slant optical depth
            if (nid(i) < 0) then
                tausla(i) = largest
            else
                sum = 0.0
                do j = 1, min(nid(i), i)
                    sum = sum + taun(j)*dsdh(i, j)
                end do
                do j = min(nid(i), i) + 1, nid(i)
                    sum = sum + 2.*taun(j)*dsdh(i, j)
                end do
                tausla(i) = sum
                if (tausla(i) == tausla(i - 1)) then
                    mu2(i) = sqrt(largest)
                else
                    mu2(i) = (tauc(i) - tauc(i - 1))/(tausla(i) - tausla(i - 1))
                    mu2(i) = sign(AMAX1(abs(mu2(i)), 1./sqrt(largest)), mu2(i))
                end if
            end if

            !** the following gamma equations are from pg 16,289, Table 1
            !** save mu1 for each approx. for use in converting irradiance to actinic flux
            ! Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):
            gam1 = (7.-om*(4.+3.*g))/4.
            gam2 = -(1.-om*(4.-3.*g))/4.
            gam3 = (2.-3.*g*mu)/4.
            gam4 = 1.-gam3
            mu1(i) = 0.5

! quadrature (Liou, 1973, JAS, 30, 1303-1326; 1974, JAS, 31, 1473-1475):
!               gam1 = 1.7320508*(2. - om*(1. + g))/2.
!               gam2 = 1.7320508*om*(1. - g)/2.
!               gam3 = (1. - 1.7320508*g*mu)/2.
!               gam4 = 1. - gam3
!               mu1(i) = 1./sqrt(3.)

! hemispheric mean (Toon et al., 1089, JGR, 94, 16287):

!               gam1 = 2. - om*(1. + g)
!               gam2 = om*(1. - g)
!               gam3 = (2. - g*mu)/4.
!               gam4 = 1. - gam3
!               mu1(i) = 0.5

! PIFM  (Zdunkovski et al.,1980, Conrib.Atmos.Phys., 53, 147-166):
!              GAM1 = 0.25*(8. - OM*(5. + 3.*G))
!              GAM2 = 0.75*OM*(1.-G)
!              GAM3 = 0.25*(2.-3.*G*MU)
!              GAM4 = 1. - GAM3
!              mu1(i) = 0.5

! delta discrete ordinates  (Schaller, 1979, Contrib.Atmos.Phys, 52, 17-26):
!              GAM1 = 0.5*1.7320508*(2. - OM*(1. + G))
!              GAM2 = 0.5*1.7320508*OM*(1.-G)
!              GAM3 = 0.5*(1.-1.7320508*G*MU)
!              GAM4 = 1. - GAM3
!              mu1(i) = 1./sqrt(3.)

! Calculations of Associated Legendre Polynomials for GAMA1,2,3,4
! in delta-function, modified quadrature, hemispheric constant,
! Hybrid modified Eddington-delta function metods, p633,Table1.
! W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
! W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440
!           YLM0 = 2.
!           YLM2 = -3.*G*MU
!           YLM4 = 0.875*G**3*MU*(5.*MU**2-3.)
!           YLM6=-0.171875*G**5*MU*(15.-70.*MU**2+63.*MU**4)
!          YLM8=+0.073242*G**7*MU*(-35.+315.*MU**2-693.*MU**4
!         *+429.*MU**6)
!          YLM10=-0.008118*G**9*MU*(315.-4620.*MU**2+18018.*MU**4
!         *-25740.*MU**6+12155.*MU**8)
!          YLM12=0.003685*G**11*MU*(-693.+15015.*MU**2-90090.*MU**4
!         *+218790.*MU**6-230945.*MU**8+88179.*MU**10)
!           YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
!           YLMS=0.25*YLMS
!           BETA0 = YLMS

!              amu1=1./1.7320508
!           YLM0 = 2.
!           YLM2 = -3.*G*amu1
!           YLM4 = 0.875*G**3*amu1*(5.*amu1**2-3.)
!           YLM6=-0.171875*G**5*amu1*(15.-70.*amu1**2+63.*amu1**4)
!          YLM8=+0.073242*G**7*amu1*(-35.+315.*amu1**2-693.*amu1**4
!         *+429.*amu1**6)
!          YLM10=-0.008118*G**9*amu1*(315.-4620.*amu1**2+18018.*amu1**4
!         *-25740.*amu1**6+12155.*amu1**8)
!          YLM12=0.003685*G**11*amu1*(-693.+15015.*amu1**2-90090.*amu1**4
!         *+218790.*amu1**6-230945.*amu1**8+88179.*amu1**10)
!           YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
!           YLMS=0.25*YLMS
!           BETA1 = YLMS

!              BETAn = 0.25*(2. - 1.5*G-0.21875*G**3-0.085938*G**5
!         *-0.045776*G**7)

! Hybrid modified Eddington-delta function(Meador and Weaver,1980,JAS,37,630):
!              subd=4.*(1.-G*G*(1.-MU))
!              GAM1 = (7.-3.*G*G-OM*(4.+3.*G)+OM*G*G*(4.*BETA0+3.*G))/subd
!              GAM2 =-(1.-G*G-OM*(4.-3.*G)-OM*G*G*(4.*BETA0+3.*G-4.))/subd
!              GAM3 = BETA0
!              GAM4 = 1. - GAM3
!              mu1(i) = (1. - g*g*(1.- mu) )/(2. - g*g)

!****
! delta function  (Meador, and Weaver, 1980, JAS, 37, 630):
!              GAM1 = (1. - OM*(1. - beta0))/MU
!              GAM2 = OM*BETA0/MU
!              GAM3 = BETA0
!              GAM4 = 1. - GAM3
!              mu1(i) = mu
!****
! modified quadrature (Meador, and Weaver, 1980, JAS, 37, 630):
!              GAM1 = 1.7320508*(1. - OM*(1. - beta1))
!              GAM2 = 1.7320508*OM*beta1
!              GAM3 = BETA0
!              GAM4 = 1. - GAM3
!              mu1(i) = 1./sqrt(3.)

! hemispheric constant (Toon et al., 1989, JGR, 94, 16287):
!              GAM1 = 2.*(1. - OM*(1. - betan))
!              GAM2 = 2.*OM*BETAn
!              GAM3 = BETA0
!              GAM4 = 1. - GAM3
!              mu1(i) = 0.5

!****

! lambda = pg 16,290 equation 21
! big gamma = pg 16,290 equation 22
! if gam2 = 0., then bgam = 0.

            lam(i) = sqrt(gam1*gam1 - gam2*gam2)

            if (gam2 /= 0.) then
                bgam(i) = (gam1 - lam(i))/gam2
            else
                bgam(i) = 0.
            end if

            expon = exp(-lam(i)*taun(i))

! e1 - e4 = pg 16,292 equation 44
            e1(i) = 1.+bgam(i)*expon
            e2(i) = 1.-bgam(i)*expon
            e3(i) = bgam(i) + expon
            e4(i) = bgam(i) - expon

! the following sets up for the C equations 23, and 24
! found on page 16,290
! prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
! which is approx equiv to shifting MU by 0.5*EPS* (MU)**3
            expon0 = exp(-tausla(i - 1))
            expon1 = exp(-tausla(i))

            divisr = lam(i)*lam(i) - 1./(mu2(i)*mu2(i))
            temp = AMAX1(eps, abs(divisr))
            divisr = sign(temp, divisr)

            up = om*pifs*((gam1 - 1./mu2(i))*gam3 + gam4*gam2)/divisr
            dn = om*pifs*((gam1 + 1./mu2(i))*gam4 + gam2*gam3)/divisr

! cup and cdn are when tau is equal to zero
! cuptn and cdntn are when tau is equal to taun
            cup(i) = up*expon0
            cdn(i) = dn*expon0
            cuptn(i) = up*expon1
            cdntn(i) = dn*expon1

        end do
!print *,'LFR-DBG PS2 07'
!**************** set up matrix ******
! ssfc = pg 16,292 equation 37  where pi Fs is one (unity).
!print *,'LFR-DBG PS2 07.01', rsfc, mu, tausla(nlayer), pifs

        ssfc = rsfc*mu*exp(-tausla(nlayer))*pifs

! MROWS = the number of rows in the matrix

        mrows = 2*nlayer

! the following are from pg 16,292  equations 39 - 43.
! set up first row of matrix:

        i = 1
        a(1) = 0.
        b(1) = e1(i)
        d(1) = -e2(i)
        e(1) = fdn0 - cdn(i)

        row = 1

! set up odd rows 3 thru (MROWS - 1):

        i = 0
        do row = 3, mrows - 1, 2
            i = i + 1
            a(row) = e2(i)*e3(i) - e4(i)*e1(i)
            b(row) = e1(i)*e1(i + 1) - e3(i)*e3(i + 1)
            d(row) = e3(i)*e4(i + 1) - e1(i)*e2(i + 1)
            e(row) = e3(i)*(cup(i + 1) - cuptn(i)) + e1(i)*(cdntn(i) - cdn(i + 1))
        end do
!print *,'LFR-DBG PS2 08'
! set up even rows 2 thru (MROWS - 2):

        i = 0
        do row = 2, mrows - 2, 2
            i = i + 1
            a(row) = e2(i + 1)*e1(i) - e3(i)*e4(i + 1)
            b(row) = e2(i)*e2(i + 1) - e4(i)*e4(i + 1)
            d(row) = e1(i + 1)*e4(i + 1) - e2(i + 1)*e3(i + 1)
            e(row) = (cup(i + 1) - cuptn(i))*e2(i + 1) - &
                     (cdn(i + 1) - cdntn(i))*e4(i + 1)
!print *,'LFR-DBG PS2 08.01', row, i, a(row), b(row), d(row), e(row)
        end do
! set up last row of matrix at MROWS:
!print *,'LFR-DBG PS2 09'
        row = mrows
        i = nlayer
!print *,'LFR-DBG PS2 09.01', row, i, e1(i), e2(i), e3(i), e4(i), cuptn(i), cdntn(i),ssfc, rsfc
        a(row) = e1(i) - rsfc*e3(i)
        b(row) = e2(i) - rsfc*e4(i)
        d(row) = 0.
        e(row) = ssfc - cuptn(i) + rsfc*cdntn(i)
!print *,'LFR-DBG PS2 10', row, mrows, a(row), b(row), d(row), e(row)
! solve tri-diagonal matrix:
        call tridag(a, b, d, e, y, mrows)
!print *,'LFR-DBG PS2 11 y=',y(:)
!*** unfold solution of matrix, compute output fluxes:
        row = 1
        lev = 1
        j = 1

! the following equations are from pg 16,291  equations 31 & 32
        fdr(lev) = exp(-tausla(0))
        edr(lev) = mu*fdr(lev)
        edn(lev) = fdn0
        eup(lev) = y(row)*e3(j) - y(row + 1)*e4(j) + cup(j)
        fdn(lev) = edn(lev)/mu1(lev)
        fup(lev) = eup(lev)/mu1(lev)
!print *,'LFR-DBG PS2 12'
        do lev = 2, nlayer + 1
!print *,'LFR-DBG PS2 12.01', lev, row, j, tausla(lev - 1), mu, e1(j),e2(j),e3(j),e4(j),cuptn(j),cdntn(j),y(row),y(row+1)
            fdr(lev) = exp(-tausla(lev - 1))
            edr(lev) = mu*fdr(lev)
            edn(lev) = y(row)*e3(j) + y(row + 1)*e4(j) + cdntn(j)
            eup(lev) = y(row)*e1(j) + y(row + 1)*e2(j) + cuptn(j)
!print *,'LFR-DBG PS2 12.02', lev, fdr(lev), edr(lev), edn(lev), eup(lev)
            if (fdr(lev) < 0) fdr(lev) = 0
            if (edr(lev) < 0) edr(lev) = 0
            if (edn(lev) < 0) edn(lev) = 0
            if (eup(lev) < 0) eup(lev) = 0
!print *,'LFR-DBG PS2 12.03', lev, fdr(lev), edr(lev), edn(lev), eup(lev)
            fdn(lev) = edn(lev)/mu1(j)
            fup(lev) = eup(lev)/mu1(j)
            row = row + 2
            j = j + 1
        end do
!print *,'LFR-DBG PS2 13'
!Debug-test LFR
        fup(nlayer + 1) = 0.0
        fdn(nlayer + 1) = 0.0
!End of debug
!_______________________________________________________________________

    end subroutine ps2str

!=============================================================================*

    subroutine tridag(a, b, c, r, u, n)
!_______________________________________________________________________
! solves tridiagonal system.  From Numerical Recipies, p. 40
!_______________________________________________________________________
        implicit none

        integer, intent(IN) :: n
        real, intent(IN)    :: a(n)
        real, intent(IN)    :: b(n)
        real, intent(IN)    :: c(n)
        real, intent(IN)    :: r(n)
        real, intent(OUT)   :: u(n)

        integer :: j

        real :: bet
        real, dimension(2*n) :: gam
!_______________________________________________________________________

        if (b(1) == 0.) stop 1001
        bet = b(1)
        u(1) = r(1)/bet
        do j = 2, n
            gam(j) = c(j - 1)/bet
            bet = b(j) - a(j)*gam(j)
            if (bet == 0.) stop 2002
            u(j) = (r(j) - a(j)*u(j - 1))/bet
        end do
        do j = n - 1, 1, -1
            u(j) = u(j) - gam(j + 1)*u(j + 1)
        end do
!_______________________________________________________________________

    end subroutine tridag

    subroutine psndo(nz, dsdh, nid, nlyr, dtauc, ssalb, pmom, &
                     albedo, nstr, numu, umu, cwt, umu0, &
                     maxcly, maxulv, maxumu, maxcmu, maxphi, rfldir, rfldn, flup, u0u, &
                     uavgso, uavgup, uavgdn, sindir, sinup, sindn)
        implicit none
! Improved handling of numerical instabilities. Bernhard Mayer on 5/3/99.
!  disort seems to produce unstable results for certain combinations
!  of single scattering albedo and phase function. A temporary fix has been
!  introduced to avoid this problem: The original instability check in
!  UPBEAM fails on certain compiler/machine combinations (e.g., gcc/LINUX,
!  or xlf/IBM RS6000). This check has therefore been replaced by a new one.
!  Whenever UPBEAM reports an instability, the single scattering albedo
!  of the respective layer is changed by a small amount, and the
!  calculation is repeated until numerically stable conditions are reached
!  (all the necessary changes are confined to the new subroutine SOLVEC
!  and the slighly changed subroutine UPBEAM). To check for potential
!  instabilities, the variable 'RCOND' returned by SGECO is compared to
!  a machine-dependent constant, 'MINRCOND'. The value of this constant
!  determines (a) if really all instabilities are caught; and (b) the
!  amount by which the single scattering albedo has to be changed. The
!  value of 'MINRCOND' is therefore a compromise between numerical
!  stability on the one hand and uncertainties introduced by changing
!  the atmospheric conditions and increased computational time on the
!  other hand (an increase of MINRCOND will lead to the detection of
!  more potential numerical instabilities, and thus to an increase in
!  computational time; by changing the atmospheric conditions, that is,
!  the single scattering albedo, the result might however be changed
!  unfavourably, if the change is too large). From a limited number
!  of experiments we found that 'MINRCOND = 5000. * R1MACH(4)' seems
!  to be a good choice if high accuracy is required (more tests are
!  definitely neccessary!). If an instability is encountered, a message
!  is printed telling about neccessary changes to the single scattering
!  albedo. This message may be switched off by setting 'DEBUG = .FALSE.'
!  in subroutine SOLVEC.

! modified to calculate sine-weighted intensities. Bernhard Mayer on 2/12/99.
! modified to handle some numerical instabilities. Chris Fischer on 1/22/99.
! modified by adding pseudo-spherical correction. Jun Zeng on 3/11/97.
! dsdh: slant path of direct beam through each layer crossed
!       when travelling from the top of the atmosphere to layer i;
!       dsdh(i,j), i = 0..nlyr, j = 1..nlyr;
! nid:  number of layers crossed by the direct beam when
!       travelling from the top of the atmosphere to layer i;
!       NID(i), i = 0..nlyr.
! uavgso, uvagup, and uvagdn are direct, downward diffuse, and upward
! diffuse actinic flux (mean intensity).
! u0u is the azimuthally averaged intensity, check DISORT.doc for details.
! *******************************************************************
!       Plane-parallel discrete ordinates radiative transfer program
!                      V E R S I O N    1.1
!             ( see DISORT.DOC for complete documentation )
! *******************************************************************

! +------------------------------------------------------------------+
!  Calling Tree (omitting calls to ERRMSG):
!  (routines in parentheses are not in this file)

!  DISORT-+-(R1MACH)
!         +-ZEROIT
!         +-CHEKIN-+-(WRTBAD)
!         |        +-(WRTDIM)
!         |        +-DREF
!         +-ZEROAL
!         +-SETDIS-+-QGAUSN (1)-+-(D1MACH)
!         +-PRTINP
!         +-LEPOLY see 2
!         +-SURFAC-+-QGAUSN see 1
!         |        +-LEPOLY see 2
!         |        +-ZEROIT
!         +-SOLEIG see 3
!         +-UPBEAM-+-(SGECO)
!         |        +-(SGESL)
!         +-TERPEV
!         +-TERPSO
!         +-SETMTX see 4
!         +-SOLVE0-+-ZEROIT
!         |        +-(SGBCO)
!         |        +-(SGBSL)
!         +-FLUXES--ZEROIT
!         +-PRAVIN
!         +-RATIO--(R1MACH)
!         +-PRTINT

! *** Intrinsic Functions used in DISORT package which take
!     non-negligible amount of time:

!    EXP :  Called by- ALBTRN, ALTRIN, CMPINT, FLUXES, SETDIS,
!                      SETMTX, SPALTR, USRINT, PLKAVG

!    SQRT : Called by- ASYMTX, LEPOLY, SOLEIG

! +-------------------------------------------------------------------+

!  Index conventions (for all DO-loops and all variable descriptions):

!     IU     :  for user polar angles

!  IQ,JQ,KQ  :  for computational polar angles ('quadrature angles')

!   IQ/2     :  for half the computational polar angles (just the ones
!               in either 0-90 degrees, or 90-180 degrees)

!     J      :  for user azimuthal angles

!     K,L    :  for Legendre expansion coefficients or, alternatively,
!               subscripts of associated Legendre polynomials

!     LU     :  for user levels

!     LC     :  for computational layers (each having a different
!               single-scatter albedo and/or phase function)

!    LEV     :  for computational levels

!    MAZIM   :  for azimuthal components in Fourier cosine expansion
!               of intensity and phase function

! +------------------------------------------------------------------+

!               I N T E R N A L    V A R I A B L E S

!   AMB(IQ/2,IQ/2)    First matrix factor in reduced eigenvalue problem
!                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG)

!   APB(IQ/2,IQ/2)    Second matrix factor in reduced eigenvalue problem
!                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG)

!   ARRAY(IQ,IQ)      Scratch matrix for SOLEIG, UPBEAM and UPISOT
!                     (see each subroutine for definition)

!   B()               Right-hand side vector of Eq. SC(5) going into
!                     SOLVE0,1;  returns as solution vector
!                     vector  L, the constants of integration

!   BDR(IQ/2,0:IQ/2)  Bottom-boundary bidirectional reflectivity for a
!                     given azimuthal component.  First index always
!                     refers to a computational angle.  Second index:
!                     if zero, refers to incident beam angle UMU0;
!                     if non-zero, refers to a computational angle.

!   BEM(IQ/2)         Bottom-boundary directional emissivity at compu-
!                     tational angles.

!   BPLANK            Intensity emitted from bottom boundary

!   CBAND()           Matrix of left-hand side of the linear system
!                     Eq. SC(5), scaled by Eq. SC(12);  in banded
!                     form required by LINPACK solution routines

!   CC(IQ,IQ)         C-sub-IJ in Eq. SS(5)

!   CMU(IQ)           Computational polar angles (Gaussian)

!   CWT(IQ)           Quadrature weights corresponding to CMU

!   DELM0             Kronecker delta, delta-sub-M0, where M = MAZIM
!                     is the number of the Fourier component in the
!                     azimuth cosine expansion

!   DITHER            Small quantity subtracted from single-scattering
!                     albedos of unity, in order to avoid using special
!                     case formulas;  prevents an eigenvalue of exactly
!                     zero from occurring, which would cause an
!                     immediate overflow

!   DTAUCP(LC)        Computational-layer optical depths (delta-M-scaled
!                     if DELTAM = TRUE, otherwise equal to DTAUC)

!   EMU(IU)           Bottom-boundary directional emissivity at user
!                     angles.

!   EVAL(IQ)          Temporary storage for eigenvalues of Eq. SS(12)

!   EVECC(IQ,IQ)      Complete eigenvectors of SS(7) on return from
!                     SOLEIG; stored permanently in  GC

!   EXPBEA(LC)        Transmission of direct beam in delta-M optical
!                     depth coordinates

!   FLYR(LC)          Truncated fraction in delta-M method

!   GL(K,LC)          Phase function Legendre polynomial expansion
!                     coefficients, calculated from PMOM by
!                     including single-scattering albedo, factor
!                     2K+1, and (if DELTAM=TRUE) the delta-M
!                     scaling

!   GC(IQ,IQ,LC)      Eigenvectors at polar quadrature angles,
!                     g  in Eq. SC(1)

!   GU(IU,IQ,LC)      Eigenvectors interpolated to user polar angles
!                     ( g  in Eqs. SC(3) and S1(8-9), i.e.
!                       G without the L factor )

!   HLPR()            Legendre coefficients of bottom bidirectional
!                     reflectivity (after inclusion of 2K+1 factor)

!   IPVT(LC*IQ)       Integer vector of pivot indices for LINPACK
!                     routines

!   KK(IQ,LC)         Eigenvalues of coeff. matrix in Eq. SS(7)

!   KCONV             Counter in azimuth convergence test

!   LAYRU(LU)         Computational layer in which user output level
!                     UTAU(LU) is located

!   LL(IQ,LC)         Constants of integration L in Eq. SC(1),
!                     obtained by solving scaled version of Eq. SC(5)

!   LYRCUT            TRUE, radiation is assumed zero below layer
!                     NCUT because of almost complete absorption

!   NAZ               Number of azimuthal components considered

!   NCUT              Computational layer number in which absorption
!                     optical depth first exceeds ABSCUT

!   OPRIM(LC)         Single scattering albedo after delta-M scaling

!   PASS1             TRUE on first entry, FALSE thereafter

!   PKAG(0:LC)        Integrated Planck function for internal emission

!   PSI(IQ)           Sum just after square bracket in  Eq. SD(9)

!   RMU(IU,0:IQ)      Bottom-boundary bidirectional reflectivity for a
!                     given azimuthal component.  First index always
!                     refers to a user angle.  Second index:
!                     if zero, refers to incident beam angle UMU0;
!                     if non-zero, refers to a computational angle.

!   TAUC(0:LC)        Cumulative optical depth (un-delta-M-scaled)

!   TAUCPR(0:LC)      Cumulative optical depth (delta-M-scaled if
!                     DELTAM = TRUE, otherwise equal to TAUC)

!   TPLANK            Intensity emitted from top boundary

!   UUM(IU,LU)        Expansion coefficients when the intensity
!                     (u-super-M) is expanded in Fourier cosine series
!                     in azimuth angle

!   U0C(IQ,LU)        Azimuthally-averaged intensity

!   UTAUPR(LU)        Optical depths of user output levels in delta-M
!                     coordinates;  equal to  UTAU(LU) if no delta-M

!   WK()              scratch array

!   XR0(LC)           X-sub-zero in expansion of thermal source func-
!                     tion preceding Eq. SS(14) (has no mu-dependence)

!   XR1(LC)           X-sub-one in expansion of thermal source func-
!                     tion;  see  Eqs. SS(14-16)

!   YLM0(L)           Normalized associated Legendre polynomial
!                     of subscript L at the beam angle (not saved
!                     as function of superscipt M)

!   YLMC(L,IQ)        Normalized associated Legendre polynomial
!                     of subscript L at the computational angles
!                     (not saved as function of superscipt M)

!   YLMU(L,IU)        Normalized associated Legendre polynomial
!                     of subscript L at the user angles
!                     (not saved as function of superscipt M)

!   Z()               scratch array used in  SOLVE0,1  to solve a
!                     linear system for the constants of integration

!   Z0(IQ)            Solution vectors Z-sub-zero of Eq. SS(16)

!   Z0U(IU,LC)        Z-sub-zero in Eq. SS(16) interpolated to user
!                     angles from an equation derived from SS(16)

!   Z1(IQ)            Solution vectors Z-sub-one  of Eq. SS(16)

!   Z1U(IU,LC)        Z-sub-one in Eq. SS(16) interpolated to user
!                     angles from an equation derived from SS(16)

!   ZBEAM(IU,LC)      Particular solution for beam source

!   ZJ(IQ)            Right-hand side vector  X-sub-zero in
!                     Eq. SS(19), also the solution vector
!                     Z-sub-zero after solving that system

!   ZZ(IQ,LC)         Permanent storage for the beam source vectors ZJ

!   ZPLK0(IQ,LC)      Permanent storage for the thermal source
!                     vectors  Z0  obtained by solving  Eq. SS(16)

!   ZPLK1(IQ,LC)      Permanent storage for the thermal source
!                     vectors  Z1  obtained by solving  Eq. SS(16)

! +-------------------------------------------------------------------+

!  LOCAL SYMBOLIC DIMENSIONS (have big effect on storage requirements):

!       MXCLY  = Max no. of computational layers
!       MXULV  = Max no. of output levels
!       MXCMU  = Max no. of computation polar angles
!       MXUMU  = Max no. of output polar angles
!       MXPHI  = Max no. of output azimuthal angles

! +-------------------------------------------------------------------+

        integer, intent(IN)   :: nz
        integer, intent(IN)   :: maxcly
        integer, intent(IN)   :: maxulv
        integer, intent(IN)   :: maxumu
        integer, intent(IN)   :: maxcmu
        integer, intent(IN)   :: maxphi
        integer, intent(IN)   :: nstr
        integer, intent(INOUT)   :: numu
        integer, intent(IN)   :: nlyr
        integer, intent(IN)   :: nid(0:nz)
        real, intent(IN)      :: dsdh(0:nz, nz)
        real, intent(IN)      :: dtauc(maxcly)
        real, intent(INOUT)      :: pmom(0:maxcmu, maxcly)
        real, intent(IN)      :: albedo
        real, intent(INOUT)      :: umu(maxumu)
        real, intent(INOUT)      :: cwt(maxcmu)
        real, intent(IN)      :: umu0
        real, intent(INOUT)   :: ssalb(maxcly)
        real, intent(OUT)     :: rfldir(maxulv)
        real, intent(OUT)     :: rfldn(maxulv)
        real, intent(OUT)     :: flup(maxulv)
        real, intent(OUT)     :: u0u(maxumu, maxulv)
        real, intent(OUT)     :: uavgso(maxulv)
        real, intent(OUT)     :: uavgup(maxulv)
        real, intent(OUT)     :: uavgdn(maxulv)
        real, intent(OUT)     :: sindir(maxulv)
        real, intent(OUT)     :: sinup(maxulv)
        real, intent(OUT)     :: sindn(maxulv)

        integer, parameter :: mxcly = 151
        integer, parameter :: mxulv = 151
        integer, parameter :: mxcmu = 32
        integer, parameter :: mxumu = 32
        integer, parameter :: mxphi = 3
        integer, parameter :: mi = mxcmu/2
        integer, parameter :: mi9m2 = 9*mi - 2
        integer, parameter :: nnlyri = mxcmu*mxcly
!     ..
!     .. Scalar Arguments ..

        integer :: ntau
        real :: btemp, temis, ttemp, wvnmhi, wvnmlo

!     sherical geometry

        real :: tausla(0:nz), tauslau(0:nz), mu2(0:nz)
!     ..
!     .. Array Arguments ..

        real :: albmed(maxumu), dfdt(maxulv), hl(0:maxcmu), phi(maxphi), &
                temper(0:maxcly), &
                trnmed(maxumu), uavg(maxulv), utau(maxulv), &
                uu(maxumu, maxulv, maxphi)

!     ..
!     .. Local Scalars ..

        logical :: lyrcut
        integer :: iq, iu, j, kconv, l, lc, lu, mazim, naz, ncol, ncos, ncut, nn
        real :: azerr, azterm, bplank, cosphi, delm0, &
                sgn, tplank
!     ..
!     .. Local Arrays ..
        real :: angcos(maxumu)
        integer :: ipvt(nnlyri), layru(mxulv)

        real :: amb(mi, mi), apb(mi, mi), array(mxcmu, mxcmu), &
                b(nnlyri), bdr(mi, 0:mi), bem(mi), &
                cband(mi9m2, nnlyri), cc(mxcmu, mxcmu), cmu(mxcmu), dtaucp(mxcly), &
                emu(mxumu), eval(mi), evecc(mxcmu, mxcmu), &
                expbea(0:mxcly), fldir(mxulv), fldn(mxulv), &
                flyr(mxcly), gc(mxcmu, mxcmu, mxcly), &
                gl(0:mxcmu, mxcly), gu(mxumu, mxcmu, mxcly), &
                hlpr(0:mxcmu), kk(mxcmu, mxcly), ll(mxcmu, mxcly), &
                oprim(mxcly), phirad(mxphi), pkag(0:mxcly), &
                psi(mxcmu), rmu(mxumu, 0:mi), tauc(0:mxcly), &
                taucpr(0:mxcly), u0c(mxcmu, mxulv), utaupr(mxulv), &
                uum(mxumu, mxulv), wk(mxcmu), xr0(mxcly), &
                xr1(mxcly), ylm0(0:mxcmu), ylmc(0:mxcmu, mxcmu), &
                ylmu(0:mxcmu, mxumu), z(nnlyri), z0(mxcmu), &
                z0u(mxumu, mxcly), z1(mxcmu), z1u(mxumu, mxcly), &
                zbeam(mxumu, mxcly), zj(mxcmu), &
                zplk0(mxcmu, mxcly), zplk1(mxcmu, mxcly), zz(mxcmu, mxcly)

!gy added glsave and dgl to allow adjustable dimensioning in SOLVEC
        real :: glsave(0:mxcmu), dgl(0:mxcmu)

        double precision :: aad(mi, mi), evald(mi), eveccd(mi, mi), wkd(mxcmu)
!     ..
!     .. External Functions ..

!  REAL :: plkavg!,  ratio
!LFR>   EXTERNAL  plkavg,  ratio
!     ..
!     .. External Subroutines ..

!LFR>   EXTERNAL  chekin, fluxes, lepoly, pravin, prtinp,  &
!LFR>       prtint, setdis, setmtx, soleig, solve0, surfac, upbeam, zeroal, zeroit
!     ..
!     .. Intrinsic Functions ..

        intrinsic ABS, ASIN, COS, LEN, MAX

! Discrete ordinate constants:
! For pseudo-spherical DISORT, PLANK, USRTAU and USRANG must be .FALSE.;
! ONLYFL must be .TRUE.; FBEAM = 1.; FISOT = 0.; IBCND = 0
        logical, parameter :: lamber = .true.
        logical, parameter :: usrtau = .false.
        logical, parameter :: plank = .false.
        logical, parameter :: usrang = .false.
        logical, parameter :: onlyfl = .true.
        logical, parameter :: deltam = .true. ! delat-M scaling option
        logical, dimension(7), parameter :: prnt = (/.false., .false., .false., .false., &
                                                     .false., .false., .false./)
        real, parameter :: accur = 0.0001
        character(LEN=127), parameter :: header = repeat(' ', 127)
        integer, parameter :: nphi = 0
        integer, parameter :: ibcnd = 0
        real, parameter :: fbeam = 1.0
        real, parameter :: fisot = 0.0
        real, parameter :: phi0 = 0.0

10      continue

!** Calculate cumulative optical depth
!   and dither single-scatter albedo
!   to improve numerical behavior of
!   eigenvalue/vector computation
        call zeroit(tauc, mxcly + 1)

        do lc = 1, nlyr

            if (ssalb(lc) == 1.0) ssalb(lc) = 1.0 - dither
            tauc(lc) = tauc(lc - 1) + dtauc(lc)

        end do
!                                ** Check input dimensions and variables

        call chekin(nlyr, dtauc, ssalb, pmom, temper, wvnmlo, wvnmhi, &
                    usrtau, ntau, utau, nstr, usrang, numu, umu, nphi, &
                    phi, ibcnd, fbeam, umu0, phi0, fisot, lamber, albedo, &
                    hl, btemp, ttemp, temis, plank, onlyfl, accur, tauc, &
                    maxcly, maxulv, maxumu, maxcmu, maxphi, mxcly, mxulv, mxumu, mxcmu, mxphi)

!                                 ** Zero internal and output arrays

        call zeroal(mxcly, expbea(1), flyr, oprim, taucpr(1), xr0, xr1, &
                    mxcmu, cmu, cwt, psi, wk, z0, z1, zj, mxcmu + 1, hlpr, ylm0, &
                    mxcmu**2, array, cc, evecc, (mxcmu + 1)*mxcly, gl, &
                    (mxcmu + 1)*mxcmu, ylmc, (mxcmu + 1)*mxumu, ylmu, &
                    mxcmu*mxcly, kk, ll, zz, zplk0, zplk1, mxcmu**2*mxcly, gc, &
                    mxulv, layru, utaupr, mxumu*mxcmu*mxcly, gu, &
                    mxumu*mxcly, z0u, z1u, zbeam, mi, eval, &
                    mi**2, amb, apb, nnlyri, ipvt, z, &
                    maxulv, rfldir, rfldn, flup, uavg, dfdt, maxumu, albmed, trnmed, &
                    maxumu*maxulv, u0u, maxumu*maxulv*maxphi, uu)

!                                 ** Perform various setup operations

        call setdis(nz, dsdh, nid, tausla, tauslau, mu2, &
                    cmu, cwt, deltam, dtauc, dtaucp, expbea, flyr, &
                    gl, hl, hlpr, ibcnd, lamber, layru, lyrcut, maxumu, &
                    maxcmu, mxcmu, ncut, nlyr, ntau, nn, nstr, plank, &
                    numu, onlyfl, oprim, pmom, ssalb, tauc, taucpr, utau, &
                    utaupr, umu, umu0, usrtau, usrang)

!                                 ** Print input information
        if (prnt(1)) call prtinp(nlyr, dtauc, dtaucp, ssalb, pmom, temper, &
                                 wvnmlo, wvnmhi, ntau, utau, nstr, numu, umu, &
                                 nphi, phi, ibcnd, fbeam, umu0, phi0, fisot, &
                                 lamber, albedo, hl, btemp, ttemp, temis, &
                                 deltam, plank, onlyfl, accur, flyr, lyrcut, &
                                 oprim, tauc, taucpr, maxcmu, prnt(7))

!                              ** Handle special case for getting albedo
!                                 and transmissivity of medium for many
!                                 beam angles at once
!                                   ** Calculate Planck functions

        bplank = 0.0
        tplank = 0.0
        call zeroit(pkag, mxcly + 1)

! ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
!           (EQ STWJ 5)

        kconv = 0
        naz = nstr - 1
!                                    ** Azimuth-independent case

        if (fbeam == 0.0 .or. (1.-umu0) < 1.e-5 .or. onlyfl .or. &
            (numu == 1 .and. (1.-umu(1)) < 1.e-5)) naz = 0

        do mazim = 0, naz

            if (mazim == 0) delm0 = 1.0
            if (mazim > 0) delm0 = 0.0

!                             ** Get normalized associated Legendre
!                                polynomials for
!                                (a) incident beam angle cosine
!                                (b) computational and user polar angle
!                                    cosines
            if (fbeam > 0.0) then

                ncos = 1
                angcos = -umu0
!Saulo: precisamos verificar isso aqui abaixo........
                call lepoly(ncos, mazim, mxcmu, nstr - 1, angcos, ylm0)

            end if

            if (.not. onlyfl .and. usrang) &
                call lepoly(numu, mazim, mxcmu, nstr - 1, umu, ylmu)

            call lepoly(nn, mazim, mxcmu, nstr - 1, cmu, ylmc)

!                       ** Get normalized associated Legendre polys.
!                          with negative arguments from those with
!                          positive arguments; Dave/Armstrong Eq. (15)
            sgn = -1.0

            do l = mazim, nstr - 1

                sgn = -sgn

                do iq = nn + 1, nstr
                    ylmc(l, iq) = sgn*ylmc(l, iq - nn)
                end do

            end do
!                                 ** Specify users bottom reflectivity
!                                    and emissivity properties
            if (.not. lyrcut) call surfac(albedo, delm0, fbeam, hlpr, lamber, &
                                          mi, mazim, mxcmu, mxumu, nn, numu, nstr, onlyfl, &
                                          umu, usrang, ylm0, ylmc, ylmu, bdr, emu, bem, rmu)

! ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

            do lc = 1, ncut

                call solvec(amb, apb, array, cmu, cwt, gl(0, lc), mi, &
                            mazim, mxcmu, nn, nstr, ylm0, ylmc, cc, &
                            evecc, eval, kk(1, lc), gc(1, 1, lc), aad, eveccd, &
                            evald, wk, wkd, delm0, fbeam, ipvt, pi, &
                            zj, zz(1, lc), oprim(lc), lc, dither, mu2(lc), glsave, dgl)
!gy added glsave and dgl to call to allow adjustable dimensioning

            end do

! ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============

!                      ** Set coefficient matrix of equations combining
!                         boundary and layer interface conditions

            call setmtx(bdr, cband, cmu, cwt, delm0, dtaucp, gc, kk, &
                        lamber, lyrcut, mi, mi9m2, mxcmu, ncol, ncut, nnlyri, nn, nstr, taucpr, wk)

!                      ** Solve for constants of integration in homo-
!                         geneous solution (general boundary conditions)

            call solve0(b, bdr, bem, bplank, cband, cmu, cwt, expbea, &
                        fbeam, fisot, ipvt, lamber, ll, lyrcut, mazim, mi, &
                        mi9m2, mxcmu, ncol, ncut, nn, nstr, nnlyri, pi, &
                        tplank, taucpr, umu0, z, zz, zplk0, zplk1)

!                                  ** Compute upward and downward fluxes

            if (mazim == 0) call fluxes(nz, tausla, tauslau, &
                                        cmu, cwt, fbeam, gc, kk, layru, ll, lyrcut, &
                                        maxulv, mxcmu, mxulv, ncut, nn, nstr, ntau, &
                                        pi, prnt, ssalb, taucpr, umu0, utau, utaupr, &
                                        xr0, xr1, zz, zplk0, zplk1, dfdt, flup, &
                                        fldn, fldir, rfldir, rfldn, uavg, u0c, uavgso, uavgup, uavgdn, &
                                        sindir, sinup, sindn)

            if (onlyfl) then

                if (maxumu >= nstr) then
!                                     ** Save azimuthal-avg intensities
!                                        at quadrature angles
                    do lu = 1, ntau

                        do iq = 1, nstr
                            u0u(iq, lu) = u0c(iq, lu)
                        end do

                    end do

                end if

                GO TO 170

            end if

            call zeroit(uum, mxumu*mxulv)

            if (mazim == 0) then
!                               ** Save azimuthally averaged intensities

                do lu = 1, ntau

                    do iu = 1, numu
                        u0u(iu, lu) = uum(iu, lu)

                        do j = 1, nphi
                            uu(iu, lu, j) = uum(iu, lu)
                        end do

                    end do

                end do
!                              ** Print azimuthally averaged intensities
!                                 at user angles

                if (prnt(4)) call pravin(umu, numu, maxumu, utau, ntau, u0u)
                if (naz > 0) then

                    call zeroit(phirad, mxphi)
                    do j = 1, nphi
                        phirad(j) = rpd*(phi(j) - phi0)
                    end do

                end if

            else
!                                ** Increment intensity by current
!                                   azimuthal component (Fourier
!                                   cosine series);  Eq SD(2)
                azerr = 0.0

                do j = 1, nphi

                    cosphi = cos(mazim*phirad(j))

                    do lu = 1, ntau

                        do iu = 1, numu
                            azterm = uum(iu, lu)*cosphi
                            uu(iu, lu, j) = uu(iu, lu, j) + azterm
                            azerr = max(azerr, ratio(abs(azterm), abs(uu(iu, lu, j))))
                        end do

                    end do

                end do

                if (azerr <= accur) kconv = kconv + 1

                if (kconv >= 2) GO TO 170

            end if

        end do

! ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============

!                                          ** Print intensities
170     continue
        if (prnt(5) .and. .not. onlyfl) call prtint(uu, utau, ntau, &
                                                    umu, numu, phi, nphi, maxulv, maxumu)

    end subroutine psndo

    subroutine asymtx(aa, evec, eval, m, ia, ievec, ier, wkd, aad, evecd, evald)

!    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

!       Solves eigenfunction problem for real asymmetric matrix
!       for which it is known a priori that the eigenvalues are real.

!       This is an adaptation of a subroutine EIGRF in the IMSL
!       library to use real instead of complex arithmetic, accounting
!       for the known fact that the eigenvalues and eigenvectors in
!       the discrete ordinate solution are real.  Other changes include
!       putting all the called subroutines in-line, deleting the
!       performance index calculation, updating many DO-loops
!       to Fortran77, and in calculating the machine precision
!       TOL instead of specifying it in a data statement.

!       EIGRF is based primarily on EISPACK routines.  The matrix is
!       first balanced using the Parlett-Reinsch algorithm.  Then
!       the Martin-Wilkinson algorithm is applied.

!       References:
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!             Sources and Development of Mathematical Software,
!             Prentice-Hall, Englewood Cliffs, NJ
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!             Clarendon Press, Oxford

!   I N P U T    V A R I A B L E S:

!       AA    :  input asymmetric matrix, destroyed after solved
!        M    :  order of  AA
!       IA    :  first dimension of  AA
!    IEVEC    :  first dimension of  EVEC

!   O U T P U T    V A R I A B L E S:

!       EVEC  :  (unnormalized) eigenvectors of  AA
!                   ( column J corresponds to EVAL(J) )

!       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )

!       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.

!   S C R A T C H   V A R I A B L E S:

!       WKD   :  work area ( dimension at least 2*M )
!       AAD   :  double precision stand-in for AA
!       EVECD :  double precision stand-in for EVEC
!       EVALD :  double precision stand-in for EVAL

!   Called by- SOLEIG
!   Calls- D1MACH, ERRMSG
! +-------------------------------------------------------------------+
        implicit none

        integer, intent(IN)   :: m
        integer, intent(IN)   :: ia
        integer, intent(IN)   :: ievec
        real, intent(IN)      :: aa(ia, m)
        real, intent(OUT)     :: evec(ievec, m)
        real, intent(OUT)     :: eval(m)
        integer, intent(OUT)  :: ier

        double precision, intent(OUT)  :: wkd(*)
        double precision, intent(OUT)  :: aad(ia, m)
        double precision, intent(OUT)  :: evecd(ia, m)
        double precision, intent(OUT)  :: evald(m)

        logical :: noconv, notlas
        integer :: i, iii, in, j, k, ka, kkk, l, lb, lll, n, n1, n2
        double precision :: col, discri, f, g, h, &
            p, q, r, repl, rnorm, row, s, scale, sgn, t, &
            uu, vv, w, x, y, z

        intrinsic ABS, DBLE, MIN, SIGN, SQRT

        double precision, parameter :: c1 = 0.4375d0
        double precision, parameter :: c2 = 0.5d0
        double precision, parameter :: c3 = 0.75d0
        double precision, parameter :: c4 = 0.95d0
        double precision, parameter :: c5 = 16.0d0
        double precision, parameter :: c6 = 256.0d0
        double precision, parameter :: zero = 0.0d0
        double precision, parameter :: one = 1.0d0

        ier = 0

        if (m < 1 .or. ia < m .or. ievec < m) &
            call errmsg('ASYMTX--bad input variable(s)', .true.)

!                           ** Handle 1x1 and 2x2 special cases

        if (m == 1) then

            eval(1) = aa(1, 1)
            evec(1, 1) = 1.0
            return

        else if (m == 2) then

            discri = (aa(1, 1) - aa(2, 2))**2 + 4.*aa(1, 2)*aa(2, 1)

            if (discri < 0.0) &
                call errmsg('ASYMTX--complex evals in 2x2 case', .true.)

            sgn = 1.0

            if (aa(1, 1) < aa(2, 2)) sgn = -1.0

            eval(1) = 0.5*(aa(1, 1) + aa(2, 2) + sgn*sqrt(discri))
            eval(2) = 0.5*(aa(1, 1) + aa(2, 2) - sgn*sqrt(discri))
            evec(1, 1) = 1.0
            evec(2, 2) = 1.0

            if (aa(1, 1) == aa(2, 2) .and. &
                (aa(2, 1) == 0.0 .or. aa(1, 2) == 0.0)) then

                rnorm = abs(aa(1, 1)) + abs(aa(1, 2)) + &
                        abs(aa(2, 1)) + abs(aa(2, 2))
                w = tol*rnorm
                evec(2, 1) = aa(2, 1)/w
                evec(1, 2) = -aa(1, 2)/w

            else

                evec(2, 1) = aa(2, 1)/(eval(1) - aa(2, 2))
                evec(1, 2) = aa(1, 2)/(eval(2) - aa(1, 1))

            end if

            return

        end if
!                               ** Put s.p. matrix into d.p. matrix
        do j = 1, m

            do k = 1, m
                aad(j, k) = dble(aa(j, k))
            end do

        end do

!                                ** Initialize output variables
        ier = 0

        do i = 1, m
            evald(i) = zero

            do j = 1, m
                evecd(i, j) = zero
            end do

            evecd(i, i) = one
        end do

!                  ** Balance the input matrix and reduce its norm by
!                     diagonal similarity transformation stored in WK;
!                     then search for rows isolating an eigenvalue
!                     and push them down
        rnorm = zero
        l = 1
        k = m

50      continue
        kkk = k

        do j = kkk, 1, -1

            row = zero

            do i = 1, k

                if (i /= j) row = row + abs(aad(j, i))

            end do

            if (row == zero) then

                wkd(k) = j

                if (j /= k) then

                    do i = 1, k
                        repl = aad(i, j)
                        aad(i, j) = aad(i, k)
                        aad(i, k) = repl
                    end do

                    do i = l, m
                        repl = aad(j, i)
                        aad(j, i) = aad(k, i)
                        aad(k, i) = repl
                    end do

                end if

                k = k - 1
                GO TO 50

            end if

        end do
!                                ** Search for columns isolating an
!                                   eigenvalue and push them left
100     continue
        lll = l

        do j = lll, k

            col = zero

            do i = l, k

                if (i /= j) col = col + abs(aad(i, j))

            end do

            if (col == zero) then

                wkd(l) = j

                if (j /= l) then

                    do i = 1, k
                        repl = aad(i, j)
                        aad(i, j) = aad(i, l)
                        aad(i, l) = repl
                    end do

                    do i = l, m
                        repl = aad(j, i)
                        aad(j, i) = aad(l, i)
                        aad(l, i) = repl
                    end do

                end if

                l = l + 1
                GO TO 100

            end if

        end do

!                           ** Balance the submatrix in rows L through K
        do i = l, k
            wkd(i) = one
        end do

160     continue
        noconv = .false.

        do i = l, k

            col = zero
            row = zero

            do j = l, k

                if (j /= i) then

                    col = col + abs(aad(j, i))
                    row = row + abs(aad(i, j))

                end if

            end do

            f = one
            g = row/c5
            h = col + row

180         continue
            if (col < g) then

                f = f*c5
                col = col*c6
                GO TO 180

            end if

            g = row*c5

190         continue
            if (col >= g) then

                f = f/c5
                col = col/c6
                GO TO 190

            end if
!                                                ** Now balance
            if ((col + row)/f < c4*h) then

                wkd(i) = wkd(i)*f
                noconv = .true.

                do j = l, m
                    aad(i, j) = aad(i, j)/f
                end do

                do j = 1, k
                    aad(j, i) = aad(j, i)*f
                end do

            end if

        end do

        if (noconv) GO TO 160
!                                   ** Is A already in Hessenberg form?
        if (k - 1 < l + 1) GO TO 370

!                                   ** Transfer A to a Hessenberg form
        do n = l + 1, k - 1

            h = zero
            wkd(n + m) = zero
            scale = zero
!                                                 ** Scale column
            do i = n, k
                scale = scale + abs(aad(i, n - 1))
            end do

            if (scale /= zero) then

                do i = k, n, -1
                    wkd(i + m) = aad(i, n - 1)/scale
                    h = h + wkd(i + m)**2
                end do

                g = -sign(sqrt(h), wkd(n + m))
                h = h - wkd(n + m)*g
                wkd(n + m) = wkd(n + m) - g
!                                            ** Form (I-(U*UT)/H)*A
                do j = n, m

                    f = zero

                    do i = k, n, -1
                        f = f + wkd(i + m)*aad(i, j)
                    end do

                    do i = n, k
                        aad(i, j) = aad(i, j) - wkd(i + m)*f/h
                    end do

                end do
!                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
                do i = 1, k

                    f = zero

                    do j = k, n, -1
                        f = f + wkd(j + m)*aad(i, j)
                    end do

                    do j = n, k
                        aad(i, j) = aad(i, j) - wkd(j + m)*f/h
                    end do

                end do

                wkd(n + m) = scale*wkd(n + m)
                aad(n, n - 1) = scale*g

            end if

        end do

        do n = k - 2, l, -1

            n1 = n + 1
            n2 = n + 2
            f = aad(n + 1, n)

            if (f /= zero) then

                f = f*wkd(n + 1 + m)

                do i = n + 2, k
                    wkd(i + m) = aad(i, n)
                end do

                if (n + 1 <= k) then

                    do j = 1, m

                        g = zero

                        do i = n + 1, k
                            g = g + wkd(i + m)*evecd(i, j)
                        end do

                        g = g/f

                        do i = n + 1, k
                            evecd(i, j) = evecd(i, j) + g*wkd(i + m)
                        end do

                    end do

                end if

            end if

        end do

370     continue

        n = 1

        do i = 1, m

            do j = n, m
                rnorm = rnorm + abs(aad(i, j))
            end do

            n = i

            if (i < l .or. i > k) evald(i) = aad(i, i)

        end do

        n = k
        t = zero

!                                      ** Search for next eigenvalues
400     continue
        if (n < l) GO TO 550

        in = 0
        n1 = n - 1
        n2 = n - 2
!                          ** Look for single small sub-diagonal element
410     continue

        do i = l, n
            lb = n + l - i

            if (lb == l) GO TO 430

            s = abs(aad(lb - 1, lb - 1)) + abs(aad(lb, lb))

            if (s == zero) s = rnorm

            if (abs(aad(lb, lb - 1)) <= tol*s) GO TO 430

        end do

430     continue
        x = aad(n, n)

        if (lb == n) then
!                                        ** One eigenvalue found
            aad(n, n) = x + t
            evald(n) = aad(n, n)
            n = n1
            GO TO 400

        end if

! next line has been included to avoid run time error caused by xlf

        if ((n1 <= 0) .or. (n <= 0)) then
            write (0, *) 'Subscript out of bounds in ASYMTX'
            stop 9999
        end if

        y = aad(n1, n1)
        w = aad(n, n1)*aad(n1, n)

        if (lb == n1) then
!                                        ** Two eigenvalues found
            p = (y - x)*c2
            q = p**2 + w
            z = sqrt(abs(q))
            aad(n, n) = x + t
            x = aad(n, n)
            aad(n1, n1) = y + t
!                                        ** Real pair
            z = p + sign(z, p)
            evald(n1) = x + z
            evald(n) = evald(n1)

            if (z /= zero) evald(n) = x - w/z

            x = aad(n, n1)
!                                  ** Employ scale factor in case
!                                     X and Z are very small
            r = sqrt(x*x + z*z)
            p = x/r
            q = z/r
!                                             ** Row modification
            do j = n1, m
                z = aad(n1, j)
                aad(n1, j) = q*z + p*aad(n, j)
                aad(n, j) = q*aad(n, j) - p*z
            end do
!                                             ** Column modification
            do i = 1, n
                z = aad(i, n1)
                aad(i, n1) = q*z + p*aad(i, n)
                aad(i, n) = q*aad(i, n) - p*z
            end do
!                                          ** Accumulate transformations
            do i = l, k
                z = evecd(i, n1)
                evecd(i, n1) = q*z + p*evecd(i, n)
                evecd(i, n) = q*evecd(i, n) - p*z
            end do

            n = n2
            GO TO 400

        end if

        if (in == 30) then

!                    ** No convergence after 30 iterations; set error
!                       indicator to the index of the current eigenvalue
            ier = n
            GO TO 700

        end if
!                                                  ** Form shift
        if (in == 10 .or. in == 20) then

            t = t + x

            do i = l, n
                aad(i, i) = aad(i, i) - x
            end do

            s = abs(aad(n, n1)) + abs(aad(n1, n2))
            x = c3*s
            y = x
            w = -c1*s**2

        end if

        in = in + 1

!                ** Look for two consecutive small sub-diagonal elements

! inhibit vectorization by CF77, as this will cause a run time error

!DIR$ NEXTSCALAR
        do j = lb, n2
            i = n2 + lb - j
            z = aad(i, i)
            r = x - z
            s = y - z
            p = (r*s - w)/aad(i + 1, i) + aad(i, i + 1)
            q = aad(i + 1, i + 1) - z - r - s
            r = aad(i + 2, i + 1)
            s = abs(p) + abs(q) + abs(r)
            p = p/s
            q = q/s
            r = r/s

            if (i == lb) GO TO 490

            uu = abs(aad(i, i - 1))*(abs(q) + abs(r))
            vv = abs(p)*(abs(aad(i - 1, i - 1)) + abs(z) + &
                         abs(aad(i + 1, i + 1)))

            if (uu <= tol*vv) GO TO 490

        end do

490     continue
        aad(i + 2, i) = zero

!                      ** fpp vectorization of this loop triggers
!                         array bounds errors, so inhibit
!FPP$ NOVECTOR L
        do j = i + 3, n
            aad(j, j - 2) = zero
            aad(j, j - 3) = zero
        end do

!             ** Double QR step involving rows K to N and columns M to N

        do ka = i, n1

            notlas = ka /= n1

            if (ka == i) then

                s = sign(sqrt(p*p + q*q + r*r), p)

                if (lb /= i) aad(ka, ka - 1) = -aad(ka, ka - 1)

            else

                p = aad(ka, ka - 1)
                q = aad(ka + 1, ka - 1)
                r = zero

                if (notlas) r = aad(ka + 2, ka - 1)

                x = abs(p) + abs(q) + abs(r)

                if (x == zero) cycle

                p = p/x
                q = q/x
                r = r/x
                s = sign(sqrt(p*p + q*q + r*r), p)
                aad(ka, ka - 1) = -s*x

            end if

            p = p + s
            x = p/s
            y = q/s
            z = r/s
            q = q/p
            r = r/p
!                                              ** Row modification
            do j = ka, m

                p = aad(ka, j) + q*aad(ka + 1, j)

                if (notlas) then

                    p = p + r*aad(ka + 2, j)
                    aad(ka + 2, j) = aad(ka + 2, j) - p*z

                end if

                aad(ka + 1, j) = aad(ka + 1, j) - p*y
                aad(ka, j) = aad(ka, j) - p*x
            end do
!                                                 ** Column modification
            do iii = 1, min(n, ka + 3)

                p = x*aad(iii, ka) + y*aad(iii, ka + 1)

                if (notlas) then

                    p = p + z*aad(iii, ka + 2)
                    aad(iii, ka + 2) = aad(iii, ka + 2) - p*r

                end if

                aad(iii, ka + 1) = aad(iii, ka + 1) - p*q
                aad(iii, ka) = aad(iii, ka) - p
            end do
!                                          ** Accumulate transformations
            do iii = l, k

                p = x*evecd(iii, ka) + y*evecd(iii, ka + 1)

                if (notlas) then

                    p = p + z*evecd(iii, ka + 2)
                    evecd(iii, ka + 2) = evecd(iii, ka + 2) - p*r

                end if

                evecd(iii, ka + 1) = evecd(iii, ka + 1) - p*q
                evecd(iii, ka) = evecd(iii, ka) - p
            end do

        end do

        GO TO 410
!                     ** All evals found, now backsubstitute real vector
550     continue

        if (rnorm /= zero) then

            do n = m, 1, -1
                n2 = n
                aad(n, n) = one

                do i = n - 1, 1, -1
                    w = aad(i, i) - evald(n)

                    if (w == zero) w = tol*rnorm

                    r = aad(i, n)

                    do j = n2, n - 1
                        r = r + aad(i, j)*aad(j, n)
                    end do

                    aad(i, n) = -r/w
                    n2 = i
                end do

            end do
!                      ** End backsubstitution vectors of isolated evals
            do i = 1, m

                if (i < l .or. i > k) then

                    do j = i, m
                        evecd(i, j) = aad(i, j)
                    end do

                end if

            end do
!                                   ** Multiply by transformation matrix
            if (k /= 0) then

                do j = m, l, -1

                    do i = l, k
                        z = zero

                        do n = l, min(j, k)
                            z = z + evecd(i, n)*aad(n, j)
                        end do

                        evecd(i, j) = z
                    end do

                end do

            end if

        end if

        do i = l, k

            do j = 1, m
                evecd(i, j) = evecd(i, j)*wkd(i)
            end do
        end do

!                           ** Interchange rows if permutations occurred
        do i = l - 1, 1, -1

            j = wkd(i)

            if (i /= j) then

                do n = 1, m
                    repl = evecd(i, n)
                    evecd(i, n) = evecd(j, n)
                    evecd(j, n) = repl
                end do

            end if

        end do

        do i = k + 1, m

            j = wkd(i)

            if (i /= j) then

                do n = 1, m
                    repl = evecd(i, n)
                    evecd(i, n) = evecd(j, n)
                    evecd(j, n) = repl
                end do

            end if

        end do

!                         ** Put results into output arrays
700     continue

        do j = 1, m

            eval(j) = evald(j)

            do k = 1, m
                evec(j, k) = evecd(j, k)
            end do

        end do

    end subroutine asymtx

    subroutine chekin(nlyr, dtauc, ssalb, pmom, temper, wvnmlo, &
                      wvnmhi, usrtau, ntau, utau, nstr, usrang, numu, &
                      umu, nphi, phi, ibcnd, fbeam, umu0, phi0, &
                      fisot, lamber, albedo, hl, btemp, ttemp, temis, &
                      plank, onlyfl, accur, tauc, maxcly, maxulv, &
                      maxumu, maxcmu, maxphi, mxcly, mxulv, mxumu, mxcmu, mxphi)

!           Checks the input dimensions and variables

!   Calls- WRTBAD, WRTDIM, DREF, ERRMSG
!   Called by- DISORT
! --------------------------------------------------------------------
        implicit none

        integer, intent(IN)  :: maxcly
        integer, intent(IN)  :: maxulv
        integer, intent(IN)  :: maxumu
        integer, intent(IN)  :: maxcmu
        integer, intent(IN)  :: maxphi
        integer, intent(IN)  :: mxcly
        integer, intent(IN)  :: mxulv
        integer, intent(IN)  :: mxumu
        integer, intent(IN)  :: mxcmu
        integer, intent(IN)  :: mxphi
        integer, intent(IN)  :: nlyr
        integer, intent(IN)  :: ntau
        integer, intent(IN)  :: nstr
        integer, intent(IN)  :: numu
        integer, intent(IN)  :: nphi
        integer, intent(IN)  :: ibcnd
        real, intent(IN)     :: dtauc(maxcly)
        real, intent(IN)     :: ssalb(maxcly)
        real, intent(IN)     :: pmom(0:maxcmu, maxcly)
        real, intent(IN)     :: temper(0:maxcly)
        real, intent(IN)     :: wvnmlo
        real, intent(IN)     :: wvnmhi
        real, intent(IN)     :: umu(maxumu)
        real, intent(IN)     :: phi(maxphi)
        real, intent(IN)     :: fbeam
        real, intent(IN)     :: umu0
        real, intent(IN)     :: phi0
        real, intent(IN)     :: fisot
        real, intent(IN)     :: albedo
        real, intent(IN)     :: hl(0:maxcmu)
        real, intent(IN)     :: btemp
        real, intent(IN)     :: ttemp
        real, intent(IN)     :: temis
        real, intent(IN)     :: accur
        real, intent(IN)     :: tauc(0:mxcly)
        real, intent(INOUT)  :: utau(maxulv)
        logical, intent(IN)  :: usrtau
        logical, intent(IN)  :: usrang
        logical, intent(IN)  :: lamber
        logical, intent(IN)  :: plank
        logical, intent(IN)  :: onlyfl

        logical :: inperr
        integer :: irmu, iu, j, k, lc, lu
        real :: flxalb, rmu

        intrinsic ABS, MOD

        inperr = .false.

        if (nlyr < 1) inperr = wrtbad('NLYR')

        if (nlyr > maxcly) inperr = wrtbad('MAXCLY')

        do lc = 1, nlyr

            if (dtauc(lc) < 0.0) inperr = wrtbad('DTAUC')

            if (ssalb(lc) < 0.0 .or. ssalb(lc) > 1.0) inperr = wrtbad('SSALB')

            if (plank .and. ibcnd /= 1) then

                if (lc == 1 .and. temper(0) < 0.0) inperr = wrtbad('TEMPER')

                if (temper(lc) < 0.0) inperr = wrtbad('TEMPER')

            end if

            do k = 0, nstr

                if (pmom(k, lc) < -1.0 .or. pmom(k, lc) > 1.0) &
                    inperr = wrtbad('PMOM')

            end do

        end do

        if (ibcnd == 1) then

            if (maxulv < 2) inperr = wrtbad('MAXULV')

        else if (usrtau) then

            if (ntau < 1) inperr = wrtbad('NTAU')

            if (maxulv < ntau) inperr = wrtbad('MAXULV')

            do lu = 1, ntau

                if (abs(utau(lu) - tauc(nlyr)) <= 1.e-4) utau(lu) = tauc(nlyr)

                if (utau(lu) < 0.0 .or. utau(lu) > tauc(nlyr)) &
                    inperr = wrtbad('UTAU')

            end do

        else

            if (maxulv < nlyr + 1) inperr = wrtbad('MAXULV')

        end if

        if (nstr < 2 .or. mod(nstr, 2) /= 0) inperr = wrtbad('NSTR')

!     IF( NSTR.EQ.2 )
!    &    CALL ERRMSG( 'CHEKIN--2 streams not recommended;'//
!    &                 ' use specialized 2-stream code instead',.False.)

        if (nstr > maxcmu) inperr = wrtbad('MAXCMU')

        if (usrang) then

            if (numu < 0) inperr = wrtbad('NUMU')

            if (.not. onlyfl .and. numu == 0) inperr = wrtbad('NUMU')

            if (numu > maxumu) inperr = wrtbad('MAXUMU')

            if (ibcnd == 1 .and. 2*numu > maxumu) inperr = wrtbad('MAXUMU')

            do iu = 1, numu

                if (umu(iu) < -1.0 .or. umu(iu) > 1.0 .or. &
                    umu(iu) == 0.0) inperr = wrtbad('UMU')

                if (ibcnd == 1 .and. umu(iu) < 0.0) inperr = wrtbad('UMU')

                if (iu > 1) then

                    if (umu(iu) < umu(iu - 1)) inperr = wrtbad('UMU')

                end if

            end do

        else

            if (maxumu < nstr) inperr = wrtbad('MAXUMU')

        end if

        if (.not. onlyfl .and. ibcnd /= 1) then

            if (nphi <= 0) inperr = wrtbad('NPHI')

            if (nphi > maxphi) inperr = wrtbad('MAXPHI')

            do j = 1, nphi

                if (phi(j) < 0.0 .or. phi(j) > 360.0) inperr = wrtbad('PHI')

            end do

        end if

        if (ibcnd < 0 .or. ibcnd > 1) inperr = wrtbad('ibcnd')

        if (ibcnd == 0) then

            if (fbeam < 0.0) inperr = wrtbad('FBEAM')

            if (fbeam > 0.0 .and. abs(umu0) > 1.0) inperr = wrtbad('UMU0')

            if (fbeam > 0.0 .and. (phi0 < 0.0 .or. phi0 > 360.0)) &
                inperr = wrtbad('PHI0')

            if (fisot < 0.0) inperr = wrtbad('FISOT')

            if (lamber) then

                if (albedo < 0.0 .or. albedo > 1.0) inperr = wrtbad('ALBEDO')

            else
!                    ** Make sure flux albedo at dense mesh of incident
!                       angles does not assume unphysical values

                do irmu = 0, 100
                    rmu = irmu*0.01
                    flxalb = dref(rmu, hl, nstr)

                    if (flxalb < 0.0 .or. flxalb > 1.0) inperr = wrtbad('HL')

                end do

            end if

        else if (ibcnd == 1) then

            if (albedo < 0.0 .or. albedo > 1.0) inperr = wrtbad('ALBEDO')

        end if

        if (plank .and. ibcnd /= 1) then

            if (wvnmlo < 0.0 .or. wvnmhi <= wvnmlo) inperr = wrtbad('WVNMLO,HI')

            if (temis < 0.0 .or. temis > 1.0) inperr = wrtbad('temis')

            if (btemp < 0.0) inperr = wrtbad('BTEMP')

            if (ttemp < 0.0) inperr = wrtbad('TTEMP')

        end if

        if (accur < 0.0 .or. accur > 1.e-2) inperr = wrtbad('accur')

        if (mxcly < nlyr) inperr = wrtdim('MXCLY', nlyr)

        if (ibcnd /= 1) then

            if (usrtau .and. mxulv < ntau) inperr = wrtdim('MXULV', ntau)

            if (.not. usrtau .and. mxulv < nlyr + 1) &
                inperr = wrtdim('MXULV', nlyr + 1)

        else

            if (mxulv < 2) inperr = wrtdim('MXULV', 2)

        end if

        if (mxcmu < nstr) inperr = wrtdim('MXCMU', nstr)

        if (usrang .and. mxumu < numu) inperr = wrtdim('MXUMU', numu)

        if (usrang .and. ibcnd == 1 .and. mxumu < 2*numu) &
            inperr = wrtdim('MXUMU', numu)

        if (.not. usrang .and. mxumu < nstr) inperr = wrtdim('MXUMU', nstr)

        if (.not. onlyfl .and. ibcnd /= 1 .and. mxphi < nphi) &
            inperr = wrtdim('MXPHI', nphi)

        if (inperr) call errmsg('DISORT--input and/or dimension errors', .true.)

        if (plank) then

            do lc = 1, nlyr

                if (abs(temper(lc) - temper(lc - 1)) > 20.0) &
                    call errmsg('CHEKIN--vertical temperature step may' &
                                //' be too large for good accuracy', .false.)
            end do

        end if

    end subroutine chekin

    subroutine fluxes(nz, tausla, tauslau, &
                      cmu, cwt, fbeam, gc, kk, layru, ll, lyrcut, &
                      maxulv, mxcmu, mxulv, ncut, nn, nstr, ntau, pi, &
                      prnt, ssalb, taucpr, umu0, utau, utaupr, xr0, &
                      xr1, zz, zplk0, zplk1, dfdt, flup, fldn, fldir, rfldir, rfldn, uavg, u0c, &
                      uavgso, uavgup, uavgdn, sindir, sinup, sindn)

!       Calculates the radiative fluxes, mean intensity, and flux
!       derivative with respect to optical depth from the m=0 intensity
!       components (the azimuthally-averaged intensity)

!    I N P U T     V A R I A B L E S:

!       CMU      :  Abscissae for Gauss quadrature over angle cosine
!       CWT      :  Weights for Gauss quadrature over angle cosine
!       GC       :  Eigenvectors at polar quadrature angles, SC(1)
!       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
!       LAYRU    :  Layer number of user level UTAU
!       LL       :  Constants of integration in Eq. SC(1), obtained
!                     by solving scaled version of Eq. SC(5);
!                     exponential term of Eq. SC(12) not included
!       LYRCUT   :  Logical flag for truncation of comput. layer
!       NN       :  Order of double-Gauss quadrature (NSTR/2)
!       NCUT     :  Number of computational layer where absorption
!                     optical depth exceeds ABSCUT
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
!       UTAUPR   :  Optical depths of user output levels in delta-M
!                     coordinates;  equal to UTAU if no delta-M
!       XR0      :  Expansion of thermal source function in Eq. SS(14)
!       XR1      :  Expansion of thermal source function Eqs. SS(16)
!       ZZ       :  Beam source vectors in Eq. SS(19)
!       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
!       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
!       (remainder are DISORT input variables)

!                   O U T P U T     V A R I A B L E S:

!       U0C      :  Azimuthally averaged intensities
!                   ( at polar quadrature angles )
!       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables)

!                   I N T E R N A L       V A R I A B L E S:

!       DIRINT   :  Direct intensity attenuated
!       FDNTOT   :  Total downward flux (direct + diffuse)
!       FLDIR    :  Direct-beam flux (delta-M scaled)
!       FLDN     :  Diffuse down-flux (delta-M scaled)
!       FNET     :  Net flux (total-down - diffuse-up)
!       FACT     :  EXP( - UTAUPR / UMU0 )
!       PLSORC   :  Planck source function (thermal)
!       ZINT     :  Intensity of m = 0 case, in Eq. SC(1)

!   Called by- DISORT
!   Calls- ZEROIT
! +-------------------------------------------------------------------+

        integer, intent(IN)  :: nz
        integer, intent(IN)  :: ncut
        integer, intent(IN)  :: nn
        integer, intent(IN)  :: nstr
        integer, intent(IN)  :: ntau
        integer, intent(IN)  :: maxulv
        integer, intent(IN)  :: mxcmu
        integer, intent(IN)  :: mxulv
        integer, intent(IN)  :: layru(mxulv)
        real, intent(IN)     :: tausla(0:nz)
        real, intent(IN)     :: tauslau(0:nz)
        real, intent(IN)     :: cmu(mxcmu)
        real, intent(IN)     :: cwt(mxcmu)
        real, intent(IN)     :: fbeam
        real, intent(IN)     :: gc(mxcmu, mxcmu, *)
        real, intent(IN)     :: kk(mxcmu, *)
        real, intent(IN)     :: ll(mxcmu, *)
        real, intent(IN)     :: pi
        real, intent(IN)     :: ssalb(*)
        real, intent(IN)     :: taucpr(0:*)
        real, intent(IN)     :: umu0
        real, intent(IN)     :: utau(maxulv)
        real, intent(IN)     :: utaupr(mxulv)
        real, intent(IN)     :: xr0(*)
        real, intent(IN)     :: xr1(*)
        real, intent(IN)     :: zz(mxcmu, *)
        real, intent(IN)     :: zplk0(mxcmu, *)
        real, intent(IN)     :: zplk1(mxcmu, *)
        logical, intent(IN)  :: lyrcut
        logical, intent(IN)  :: prnt(*)
        real, intent(INOUT)  :: uavg(maxulv)
        real, intent(INOUT)  :: flup(maxulv)
        real, intent(OUT)    :: dfdt(maxulv)
        real, intent(OUT)    :: fldn(mxulv)
        real, intent(OUT)    :: fldir(mxulv)
        real, intent(OUT)    :: rfldir(maxulv)
        real, intent(OUT)    :: rfldn(maxulv)
        real, intent(OUT)    :: u0c(mxcmu, mxulv)
        real, intent(OUT)    :: uavgso(*)
        real, intent(OUT)    :: uavgup(*)
        real, intent(OUT)    :: uavgdn(*)
        real, intent(OUT)    :: sindir(*)
        real, intent(OUT)    :: sinup(*)
        real, intent(OUT)    :: sindn(*)

        integer :: iq, jq, lu, lyu
        real :: ang1, ang2, dirint, fact, fdntot, fnet, plsorc, zint
!     ..

        intrinsic ACOS, EXP

        if (prnt(2)) write (*, 9000)
!                                          ** Zero DISORT output arrays
        call zeroit(u0c, mxulv*mxcmu)
        call zeroit(fldir, mxulv)
        call zeroit(fldn, mxulv)
        call zeroit(uavgso, maxulv)
        call zeroit(uavgup, maxulv)
        call zeroit(uavgdn, maxulv)
        call zeroit(sindir, maxulv)
        call zeroit(sinup, maxulv)
        call zeroit(sindn, maxulv)

!                                        ** Loop over user levels
        do lu = 1, ntau

            lyu = layru(lu)

            if (lyrcut .and. lyu > ncut) then
!                                                ** No radiation reaches
!                                                ** this level
                fdntot = 0.0
                fnet = 0.0
                plsorc = 0.0
                GO TO 70

            end if

            if (fbeam > 0.0) then

                fact = exp(-tausla(lu - 1))
                dirint = fbeam*fact
                fldir(lu) = umu0*(fbeam*fact)
                rfldir(lu) = umu0*fbeam*exp(-tauslau(lu - 1))
                sindir(lu) = sqrt(1.-umu0*umu0)*fbeam*exp(-tauslau(lu - 1))

            else

                dirint = 0.0
                fldir(lu) = 0.0
                rfldir(lu) = 0.0
                sindir(lu) = 0.0

            end if

            do iq = 1, nn

                zint = 0.0

                do jq = 1, nn
                    zint = zint + gc(iq, jq, lyu)*ll(jq, lyu)* &
                           exp(-kk(jq, lyu)*(utaupr(lu) - taucpr(lyu)))
                end do

                do jq = nn + 1, nstr
                    zint = zint + gc(iq, jq, lyu)*ll(jq, lyu)* &
                           exp(-kk(jq, lyu)*(utaupr(lu) - taucpr(lyu - 1)))
                end do

                u0c(iq, lu) = zint

                if (fbeam > 0.0) u0c(iq, lu) = zint + zz(iq, lyu)*fact

                u0c(iq, lu) = u0c(iq, lu) + zplk0(iq, lyu) + &
                              zplk1(iq, lyu)*utaupr(lu)
                uavg(lu) = uavg(lu) + cwt(nn + 1 - iq)*u0c(iq, lu)
                uavgdn(lu) = uavgdn(lu) + cwt(nn + 1 - iq)*u0c(iq, lu)
                sindn(lu) = sindn(lu) + cwt(nn + 1 - iq)* &
                            sqrt(1.-cmu(nn + 1 - iq)*cmu(nn + 1 - iq))*u0c(iq, lu)
                fldn(lu) = fldn(lu) + cwt(nn + 1 - iq)* &
                           cmu(nn + 1 - iq)*u0c(iq, lu)
            end do

            do iq = nn + 1, nstr

                zint = 0.0

                do jq = 1, nn
                    zint = zint + gc(iq, jq, lyu)*ll(jq, lyu)* &
                           exp(-kk(jq, lyu)*(utaupr(lu) - taucpr(lyu)))
                end do

                do jq = nn + 1, nstr
                    zint = zint + gc(iq, jq, lyu)*ll(jq, lyu)* &
                           exp(-kk(jq, lyu)*(utaupr(lu) - taucpr(lyu - 1)))
                end do

                u0c(iq, lu) = zint

                if (fbeam > 0.0) u0c(iq, lu) = zint + zz(iq, lyu)*fact

                u0c(iq, lu) = u0c(iq, lu) + zplk0(iq, lyu) + &
                              zplk1(iq, lyu)*utaupr(lu)
                uavg(lu) = uavg(lu) + cwt(iq - nn)*u0c(iq, lu)
                uavgup(lu) = uavgup(lu) + cwt(iq - nn)*u0c(iq, lu)
                sinup(lu) = sinup(lu) + cwt(iq - nn)*sqrt(1.-cmu(iq - nn)*cmu(iq - nn))* &
                            u0c(iq, lu)
                flup(lu) = flup(lu) + cwt(iq - nn)*cmu(iq - nn)*u0c(iq, lu)
            end do

            flup(lu) = 2.*pi*flup(lu)
            fldn(lu) = 2.*pi*fldn(lu)
            fdntot = fldn(lu) + fldir(lu)
            fnet = fdntot - flup(lu)
            rfldn(lu) = fdntot - rfldir(lu)
            uavg(lu) = (2.*pi*uavg(lu) + dirint)/(4.*pi)
            uavgso(lu) = dirint/(4.*pi)
            uavgup(lu) = (2.0*pi*uavgup(lu))/(4.*pi)
            uavgdn(lu) = (2.0*pi*uavgdn(lu))/(4.*pi)
            sindn(lu) = 2.*pi*sindn(lu)
            sinup(lu) = 2.*pi*sinup(lu)

            plsorc = xr0(lyu) + xr1(lyu)*utaupr(lu)
            dfdt(lu) = (1.-ssalb(lyu))*4.*pi*(uavg(lu) - plsorc)

70          continue
            if (prnt(2)) write (*, FMT=9010) utau(lu), lyu, &
                rfldir(lu), rfldn(lu), fdntot, flup(lu), fnet, &
                uavg(lu), plsorc, dfdt(lu)

        end do

        if (prnt(3)) then

            write (*, FMT=9020)

            do lu = 1, ntau

                write (*, FMT=9030) utau(lu)

                do iq = 1, nn
                    ang1 = 180./pi*acos(cmu(2*nn - iq + 1))
                    ang2 = 180./pi*acos(cmu(iq))
                    write (*, 9040) ang1, cmu(2*nn - iq + 1), u0c(iq, lu), &
                        ang2, cmu(iq), u0c(iq + nn, lu)
                end do

            end do

        end if

9000    format(//, 21x, &
                '<----------------------- FLUXES ----------------------->', /, &
                '   Optical  Compu    Downward    Downward    Downward     ', &
                ' Upward                    Mean      Planck   d(Net Flux)', /, &
                '     Depth  Layer      Direct     Diffuse       Total     ', &
                'Diffuse         Net   Intensity      Source   / d(Op Dep)',/)
9010    format(f10.4, i7, 1p, 7e12.3, e14.3)
9020    format(/, /, ' ******** AZIMUTHALLY AVERAGED INTENSITIES', &
                ' ( at polar quadrature angles ) *******')
9030    format(/, ' Optical depth =', f10.4, //, &
                '     Angle (deg)   cos(Angle)     Intensity', &
                '     Angle (deg)   cos(Angle)     Intensity')
9040    format(2(0p, f16.4, f13.5, 1p, e14.3))

    end subroutine fluxes

    subroutine lepoly(nmu, m, maxmu, twonm1, mu, ylm)

!       Computes the normalized associated Legendre polynomial,
!       defined in terms of the associated Legendre polynomial
!       Plm = P-sub-l-super-m as

!             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)

!       for fixed order m and all degrees from l = m to TWONM1.
!       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
!       from a prior call to the routine.

!       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
!                  High-Order Associated Legendre Polynomials,
!                  J. Quant. Spectrosc. Radiat. Transfer 10,
!                  557-562, 1970.  (hereafter D/A)

!       METHOD: Varying degree recurrence relationship.

!       NOTE 1: The D/A formulas are transformed by
!               setting  M = n-1; L = k-1.
!       NOTE 2: Assumes that routine is called first with  M = 0,
!               then with  M = 1, etc. up to  M = TWONM1.
!       NOTE 3: Loops are written in such a way as to vectorize.

!  I N P U T     V A R I A B L E S:

!       NMU    :  Number of arguments of YLM
!       M      :  Order of YLM
!       MAXMU  :  First dimension of YLM
!       TWONM1 :  Max degree of YLM
!       MU(i)  :  Arguments of YLM (i = 1 to NMU)

!       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist
!       from a prior call.

!  O U T P U T     V A R I A B L E:

!       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre
!                   polynomials evaluated at argument MU(i)

!   Called by- DISORT, ALBTRN, SURFAC
!   Calls- ERRMSG
! +-------------------------------------------------------------------+
        implicit none

        integer, intent(IN)    :: nmu
        integer, intent(IN)    :: m
        integer, intent(IN)    :: maxmu
        integer, intent(IN)    :: twonm1
        real, intent(IN)       :: mu(*)
        real, intent(OUT)      :: ylm(0:maxmu, *)

        integer :: i, l, ns
        real :: tmp1, tmp2

        if (pass1) then

            pass1 = .false.

            do ns = 1, maxsqt
                sqt(ns) = sqrt(FLOAT(ns))
            end do

        end if

        if (2*twonm1 > maxsqt) &
            call errmsg('LEPOLY--need to increase param MAXSQT', .true.)

        if (m == 0) then
!                             ** Upward recurrence for ordinary
!                                Legendre polynomials
            do i = 1, nmu
                ylm(0, i) = 1.0
                ylm(1, i) = mu(i)
            end do

            do l = 2, twonm1

                do i = 1, nmu
                    ylm(l, i) = ((2*l - 1)*mu(i)*ylm(l - 1, i) - &
                                 (l - 1)*ylm(l - 2, i))/l
                end do

            end do

        else

            do i = 1, nmu
!                               ** Y-sub-m-super-m; derived from
!                               ** D/A Eqs. (11,12)

                ylm(m, i) = -sqt(2*m - 1)/sqt(2*m)* &
                            sqrt(1.-mu(i)**2)*ylm(m - 1, i)

!                              ** Y-sub-(m+1)-super-m; derived from
!                              ** D/A Eqs.(13,14) using Eqs.(11,12)

                ylm(m + 1, i) = sqt(2*m + 1)*mu(i)*ylm(m, i)

            end do

!                                   ** Upward recurrence; D/A EQ.(10)
            do l = m + 2, twonm1

                tmp1 = sqt(l - m)*sqt(l + m)
                tmp2 = sqt(l - m - 1)*sqt(l + m - 1)

                do i = 1, nmu
                    ylm(l, i) = ((2*l - 1)*mu(i)*ylm(l - 1, i) - &
                                 tmp2*ylm(l - 2, i))/tmp1
                end do

            end do

        end if

    end subroutine lepoly

    subroutine pravin(umu, numu, maxumu, utau, ntau, u0u)

!        Print azimuthally averaged intensities at user angles

!   Called by- DISORT

!     LENFMT   Max number of polar angle cosines UMU that can be
!                printed on one line, as set in FORMAT statement
! --------------------------------------------------------------------
        implicit none

        integer, intent(IN)  :: numu
        integer, intent(IN)  :: maxumu
        integer, intent(IN)  :: ntau
        real, intent(IN)     :: umu(numu)
        real, intent(IN)     :: utau(ntau)
        real, intent(IN)     :: u0u(maxumu, ntau)

        integer :: iu, iumax, iumin, lenfmt, lu, np, npass

        intrinsic MIN

        if (numu < 1) return

        write (*, '(//,A)') ' *******  AZIMUTHALLY AVERAGED INTENSITIES '// &
            '(at user polar angles)  ********'

        lenfmt = 8
        npass = 1 + (numu - 1)/lenfmt

        write (*, '(/,A,/,A)') '   Optical   Polar Angle Cosines', '     Depth'

        do np = 1, npass

            iumin = 1 + lenfmt*(np - 1)
            iumax = min(lenfmt*np, numu)
            write (*, '(/,10X,8F14.5)') (umu(iu), iu=iumin, iumax)

            do lu = 1, ntau
                write (*, '(0P,F10.4,1P,8E14.4)') utau(lu), &
                    (u0u(iu, lu), iu=iumin, iumax)
            end do

        end do

    end subroutine pravin

    subroutine prtinp(nlyr, dtauc, dtaucp, ssalb, pmom, temper, &
                      wvnmlo, wvnmhi, ntau, utau, nstr, numu, umu, &
                      nphi, phi, ibcnd, fbeam, umu0, phi0, fisot, &
                      lamber, albedo, hl, btemp, ttemp, temis, &
                      deltam, plank, onlyfl, accur, flyr, lyrcut, &
                      oprim, tauc, taucpr, maxcmu, prtmom)

!        Print values of input variables

!   Called by- DISORT
! --------------------------------------------------------------------
        implicit none

        integer, intent(IN)  :: nlyr
        integer, intent(IN)  :: nstr
        integer, intent(IN)  :: numu
        integer, intent(IN)  :: ntau
        integer, intent(IN)  :: nphi
        integer, intent(IN)  :: ibcnd
        integer, intent(IN)  :: maxcmu
        real, intent(IN)     :: dtauc(*)
        real, intent(IN)     :: dtaucp(*)
        real, intent(IN)     :: ssalb(*)
        real, intent(IN)     :: pmom(0:maxcmu, *)
        real, intent(IN)     :: temper(0:*)
        real, intent(IN)     :: wvnmlo
        real, intent(IN)     :: wvnmhi
        real, intent(IN)     :: utau(*)
        real, intent(IN)     :: umu(*)
        real, intent(IN)     :: phi(*)
        real, intent(IN)     :: fbeam
        real, intent(IN)     :: umu0
        real, intent(IN)     :: phi0
        real, intent(IN)     :: fisot
        real, intent(IN)     :: albedo
        real, intent(IN)     :: hl(0:maxcmu)
        real, intent(IN)     :: btemp
        real, intent(IN)     :: ttemp
        real, intent(IN)     :: temis
        real, intent(IN)     :: accur
        real, intent(IN)     :: flyr(*)
        real, intent(IN)     :: oprim(*)
        real, intent(IN)     :: tauc(0:*)
        real, intent(IN)     :: taucpr(0:*)
        logical, intent(IN)  :: prtmom
        logical, intent(IN)  :: lamber
        logical, intent(IN)  :: deltam
        logical, intent(IN)  :: plank
        logical, intent(IN)  :: onlyfl
        logical, intent(IN)  :: lyrcut

        integer :: iu, j, k, lc, lu
        real :: yessct
!     ..

        write (*, '(/,A,I4,A,I4)') ' No. streams =', nstr, &
            '     No. computational layers =', nlyr

        if (ibcnd /= 1) write (*, '(I4,A,10F10.4,/,(26X,10F10.4))') &
            ntau, ' User optical depths :', (utau(lu), lu=1, ntau)

        if (.not. onlyfl) write (*, '(I4,A,10F9.5,/,(31X,10F9.5))') &
            numu, ' User polar angle cosines :', (umu(iu), iu=1, numu)

        if (.not. onlyfl .and. ibcnd /= 1) &
            write (*, '(I4,A,10F9.2,/,(28X,10F9.2))') &
            nphi, ' User azimuthal angles :', (phi(j), j=1, nphi)

        if (.not. plank .or. ibcnd == 1) write (*, '(A)') ' No thermal emission'

        write (*, '(A,I2)') ' Boundary condition flag: IBCND =', ibcnd

        if (ibcnd == 0) then

            write (*, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)') &
                '    Incident beam with intensity =', fbeam, &
                ' and polar angle cosine = ', umu0, '  and azimuth angle =', phi0, &
                '    plus isotropic incident intensity =', fisot

            if (lamber) write (*, '(A,0P,F8.4)') &
                '    Bottom albedo (Lambertian) =', albedo

            if (.not. lamber) write (*, '(A,/,(10X,10F9.5))') &
                '    Legendre coeffs of bottom bidirectional reflectivity :', &
                (hl(k), k=0, nstr)

            if (plank) write (*, '(A,2F14.4,/,A,F10.2,A,F10.2,A,F8.4)') &
                '    Thermal emission in wavenumber interval :', wvnmlo, wvnmhi, &
                '    Bottom temperature =', btemp, '    Top temperature =', ttemp, &
                '    Top emissivity =', temis

        else if (ibcnd == 1) then

            write (*, '(A)') '    Isotropic illumination from top and bottom'
            write (*, '(A,0P,F8.4)') '    Bottom albedo (Lambertian) =', albedo
        end if

        if (deltam) write (*, '(A)') ' Uses delta-M method'
        if (.not. deltam) write (*, '(A)') ' Does not use delta-M method'

        if (ibcnd == 1) then

            write (*, '(A)') ' Calculate albedo and transmissivity of'// &
                ' medium vs. incident beam angle'

        else if (onlyfl) then

            write (*, '(A)') ' Calculate fluxes and azim-averaged intensities only'

        else

            write (*, '(A)') ' Calculate fluxes and intensities'

        end if

        write (*, '(A,1P,E11.2)') &
            ' Relative convergence criterion for azimuth series =', accur

        if (lyrcut) write (*, '(A)') &
            ' Sets radiation = 0 below absorption optical depth 10'

!                                        ** Print layer variables
        if (plank) write (*, FMT=9180)
        if (.not. plank) write (*, FMT=9190)

        yessct = 0.0

        do lc = 1, nlyr

            yessct = yessct + ssalb(lc)

            if (plank) write (*, '(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)') &
                lc, dtauc(lc), tauc(lc), ssalb(lc), flyr(lc), &
                dtaucp(lc), taucpr(lc), oprim(lc), pmom(1, lc), temper(lc - 1)

            if (.not. plank) write (*, '(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)') &
                lc, dtauc(lc), tauc(lc), ssalb(lc), flyr(lc), &
                dtaucp(lc), taucpr(lc), oprim(lc), pmom(1, lc)
        end do

        if (plank) write (*, '(85X,F14.3)') temper(nlyr)

        if (prtmom .and. yessct > 0.0) then

            write (*, '(/,A)') ' Layer   Phase Function Moments'

            do lc = 1, nlyr

                if (ssalb(lc) > 0.0) write (*, '(I6,10F11.6,/,(6X,10F11.6))') &
                    lc, (pmom(k, lc), k=0, nstr)
            end do

        end if

!                ** (Read every other line in these formats)

9180    format(/, 37x, '<------------- Delta-M --------------->', /, &
                '                   Total    Single                           ', &
                'Total    Single', /, '       Optical   Optical   Scatter   Truncated   ', &
                'Optical   Optical   Scatter    Asymm', /, &
                '         Depth     Depth    Albedo    Fraction     ', &
                'Depth     Depth    Albedo   Factor   Temperature')
9190    format(/, 37x, '<------------- Delta-M --------------->', /, &
                '                   Total    Single                           ', &
                'Total    Single', /, '       Optical   Optical   Scatter   Truncated   ', &
                'Optical   Optical   Scatter    Asymm', /, &
                '         Depth     Depth    Albedo    Fraction     ', &
                'Depth     Depth    Albedo   Factor')

    end subroutine prtinp

    subroutine prtint(uu, utau, ntau, umu, numu, phi, nphi, maxulv, maxumu)

!         Prints the intensity at user polar and azimuthal angles

!     All arguments are DISORT input or output variables

!   Called by- DISORT

!     LENFMT   Max number of azimuth angles PHI that can be printed
!                on one line, as set in FORMAT statement
! +-------------------------------------------------------------------+
        implicit none

        integer, intent(IN) :: nphi
        integer, intent(IN) :: maxulv
        integer, intent(IN) :: maxumu
        integer, intent(IN) :: ntau
        integer, intent(IN) :: numu
        real, intent(IN)    :: uu(maxumu, maxulv, *)
        real, intent(IN)    :: utau(*)
        real, intent(IN)    :: umu(*)
        real, intent(IN)    :: phi(*)

        integer :: iu, j, jmax, jmin, lenfmt, lu, np, npass

        intrinsic MIN
!     ..

        if (nphi < 1) return

        write (*, '(//,A)') ' *********  I N T E N S I T I E S  *********'

        lenfmt = 10
        npass = 1 + (nphi - 1)/lenfmt

        write (*, '(/,A,/,A,/,A)') &
            '             Polar   Azimuth angles (degrees)', '   Optical   Angle', &
            '    Depth   Cosine'

        do lu = 1, ntau

            do np = 1, npass

                jmin = 1 + lenfmt*(np - 1)
                jmax = min(lenfmt*np, nphi)

                write (*, '(/,18X,10F11.2)') (phi(j), j=jmin, jmax)

                if (np == 1) write (*, '(F10.4,F8.4,1P,10E11.3)') &
                    utau(lu), umu(1), (uu(1, lu, j), j=jmin, jmax)
                if (np > 1) write (*, '(10X,F8.4,1P,10E11.3)') &
                    umu(1), (uu(1, lu, j), j=jmin, jmax)

                do iu = 2, numu
                    write (*, '(10X,F8.4,1P,10E11.3)') &
                        umu(iu), (uu(iu, lu, j), j=jmin, jmax)
                end do

            end do

        end do

    end subroutine prtint

    subroutine qgausn(m, gmu, gwt)

!       Compute weights and abscissae for ordinary Gaussian quadrature
!       on the interval (0,1);  that is, such that

!           sum(i=1 to M) ( GWT(i) f(GMU(i)) )

!       is a good approximation to

!           integral(0 to 1) ( f(x) dx )

!   INPUT :    M       order of quadrature rule

!   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
!             GWT(I)   array of weights (I = 1 TO M)

!   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
!                   Integration, Academic Press, New York, pp. 87, 1975

!   METHOD:  Compute the abscissae as roots of the Legendre
!            polynomial P-sub-M using a cubically convergent
!            refinement of Newton's method.  Compute the
!            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
!            that Newton's method can very easily diverge; only a
!            very good initial guess can guarantee convergence.
!            The initial guess used here has never led to divergence
!            even for M up to 1000.

!   ACCURACY:  relative error no better than TOL or computer
!              precision (machine epsilon), whichever is larger

!   INTERNAL VARIABLES:

!    ITER      : number of Newton Method iterations
!    MAXIT     : maximum allowed iterations of Newton Method
!    PM2,PM1,P : 3 successive Legendre polynomials
!    PPR       : derivative of Legendre polynomial
!    P2PRI     : 2nd derivative of Legendre polynomial
!    TOL       : convergence criterion for Legendre poly root iteration
!    X,XI      : successive iterates in cubically-convergent version
!                of Newtons Method (seeking roots of Legendre poly.)

!   Called by- SETDIS, SURFAC
!   Calls- D1MACH, ERRMSG
! +-------------------------------------------------------------------+
        implicit none

        integer, intent(IN) :: m
        real, intent(OUT)   :: gmu(m)
        real, intent(OUT)   :: gwt(m)

        integer :: iter, k, lim, nn, np1
        real :: cona, t
        double precision :: en, nnp1, p, p2pri, pm1, pm2, ppr, prod, &
            tmp, x, xi

        intrinsic ABS, ASIN, COS, FLOAT, MOD, TAN

        integer, parameter :: maxit = 1000
        double precision, parameter :: one = 1.0d0
        double precision, parameter :: two = 2.0d0

        if (m < 1) call errmsg('QGAUSN--Bad value of M', .true.)

        if (m == 1) then

            gmu(1) = 0.5
            gwt(1) = 1.0
            return

        end if

        en = m
        np1 = m + 1
        nnp1 = m*np1
        cona = FLOAT(m - 1)/(8*m**3)

        lim = m/2

        do k = 1, lim
!                                        ** Initial guess for k-th root
!                                           of Legendre polynomial, from
!                                           Davis/Rabinowitz (2.7.3.3a)
            t = (4*k - 1)*pi/(4*m + 2)
            x = cos(t + cona/tan(t))
            iter = 0
!                                        ** Upward recurrence for
!                                           Legendre polynomials
10          continue
            iter = iter + 1
            pm2 = one
            pm1 = x

            do nn = 2, m
                p = ((2*nn - 1)*x*pm1 - (nn - 1)*pm2)/nn
                pm2 = pm1
                pm1 = p
            end do
!                                              ** Newton Method
            tmp = one/(one - x**2)
            ppr = en*(pm2 - x*p)*tmp
            p2pri = (two*x*ppr - nnp1*p)*tmp
            xi = x - (p/ppr)*(one + (p/ppr)*p2pri/(two*ppr))

!                                              ** Check for convergence
            if (abs(xi - x) > tol) then

                if (iter > maxit) call errmsg('QGAUSN--max iteration count', .true.)

                x = xi
                GO TO 10

            end if
!                             ** Iteration finished--calculate weights,
!                                abscissae for (-1,1)
            gmu(k) = -x
            gwt(k) = two/(tmp*(en*pm2)**2)
            gmu(np1 - k) = -gmu(k)
            gwt(np1 - k) = gwt(k)
        end do
!                                    ** Set middle abscissa and weight
!                                       for rules of odd order
        if (mod(m, 2) /= 0) then

            gmu(lim + 1) = 0.0
            prod = one

            do k = 3, m, 2
                prod = prod*k/(k - 1)
            end do

            gwt(lim + 1) = two/prod**2
        end if

!                                        ** Convert from (-1,1) to (0,1)
        do k = 1, m
            gmu(k) = 0.5*gmu(k) + 0.5
            gwt(k) = 0.5*gwt(k)
        end do

    end subroutine qgausn

    subroutine setdis(nz, dsdh, nid, tausla, tauslau, mu2, &
                      cmu, cwt, deltam, dtauc, dtaucp, expbea, &
                      flyr, gl, hl, hlpr, ibcnd, lamber, layru, &
                      lyrcut, maxumu, maxcmu, mxcmu, ncut, nlyr, &
                      ntau, nn, nstr, plank, numu, onlyfl, oprim, &
                      pmom, ssalb, tauc, taucpr, utau, utaupr, umu, umu0, usrtau, usrang)

!          Perform miscellaneous setting-up operations

!       INPUT :  all are DISORT input variables (see DOC file)

!       OUTPUT:  NTAU,UTAU   if USRTAU = FALSE
!                NUMU,UMU    if USRANG = FALSE
!                CMU,CWT     computational polar angles and
!                               corresponding quadrature weights
!                EXPBEA      transmission of direct beam
!                FLYR        truncated fraction in delta-M method
!                GL          phase function Legendre coefficients multi-
!                              plied by (2L+1) and single-scatter albedo
!                HLPR        Legendre moments of surface bidirectional
!                              reflectivity, times 2K+1
!                LAYRU       Computational layer in which UTAU falls
!                LYRCUT      flag as to whether radiation will be zeroed
!                              below layer NCUT
!                NCUT        computational layer where absorption
!                              optical depth first exceeds  ABSCUT
!                NN          NSTR / 2
!                OPRIM       delta-M-scaled single-scatter albedo
!                TAUCPR      delta-M-scaled optical depth
!                UTAUPR      delta-M-scaled version of  UTAU

!   Called by- DISORT
!   Calls- QGAUSN, ERRMSG
! ----------------------------------------------------------------------
        implicit none

        integer, intent(IN)  :: nz
        integer, intent(IN)  :: maxumu
        integer, intent(IN)  :: maxcmu
        integer, intent(IN)  :: mxcmu
        integer, intent(IN)  :: nlyr
        integer, intent(IN)  :: nstr
        integer, intent(IN)  :: ibcnd
        integer, intent(IN)  :: nid(0:nz)
        real, intent(IN)     :: dsdh(0:nz, nz)
        real, intent(IN)     :: dtauc(*)
        real, intent(IN)     :: hl(0:maxcmu)
        real, intent(IN)     :: ssalb(*)
        real, intent(IN)     :: tauc(0:*)
        real, intent(IN)     :: umu0
        logical, intent(IN)  :: deltam
        logical, intent(IN)  :: lamber
        logical, intent(IN)  :: plank
        logical, intent(IN)  :: onlyfl
        logical, intent(IN)  :: usrtau
        logical, intent(IN)  :: usrang
        integer, intent(OUT) :: layru(*)
        integer, intent(OUT) :: numu
        integer, intent(OUT) :: ncut
        integer, intent(OUT) :: ntau
        integer, intent(OUT) :: nn
        real, intent(OUT)    :: tausla(0:nz)
        real, intent(OUT)    :: tauslau(0:nz)
        real, intent(OUT)    :: mu2(0:nz)
        real, intent(OUT)    :: cmu(mxcmu)
        real, intent(OUT)    :: cwt(mxcmu)
        real, intent(OUT)    :: dtaucp(*)
        real, intent(OUT)    :: expbea(0:*)
        real, intent(OUT)    :: flyr(*)
        real, intent(OUT)    :: gl(0:mxcmu, *)
        real, intent(OUT)    :: hlpr(0:mxcmu)
        real, intent(OUT)    :: oprim(*)
        real, intent(OUT)    :: pmom(0:maxcmu, *)
        real, intent(OUT)    :: taucpr(0:*)
        real, intent(OUT)    :: utau(*)
        real, intent(OUT)    :: utaupr(*)
        real, intent(OUT)    :: umu(maxumu)
        logical, intent(OUT) :: lyrcut

        real :: sum, sumu

        integer :: iq, iu, k, lc, lu, i
        real :: abstau, f
        real, parameter :: abscut = 10.0

        intrinsic ABS, EXP

        if (.not. usrtau) then
!                              ** Set output levels at computational
!                                 layer boundaries
            ntau = nlyr + 1

            do lc = 0, ntau - 1
                utau(lc + 1) = tauc(lc)
            end do

        end if
!                        ** Apply delta-M scaling and move description
!                           of computational layers to local variables
        expbea(0) = 1.0
        taucpr(0) = 0.0
        abstau = 0.0
        do i = 0, nz
            tausla(i) = 0.0
            tauslau(i) = 0.0
            mu2(i) = 1./largest
        end do

        do lc = 1, nlyr

            pmom(0, lc) = 1.0

            if (abstau < abscut) ncut = lc

            abstau = abstau + (1.-ssalb(lc))*dtauc(lc)

            if (.not. deltam) then

                oprim(lc) = ssalb(lc)
                dtaucp(lc) = dtauc(lc)
                taucpr(lc) = tauc(lc)

                do k = 0, nstr - 1
                    gl(k, lc) = (2*k + 1)*oprim(lc)*pmom(k, lc)
                end do

                f = 0.0

            else
!                                    ** Do delta-M transformation

                f = pmom(nstr, lc)
                oprim(lc) = ssalb(lc)*(1.-f)/(1.-f*ssalb(lc))
                dtaucp(lc) = (1.-f*ssalb(lc))*dtauc(lc)
                taucpr(lc) = taucpr(lc - 1) + dtaucp(lc)

                do k = 0, nstr - 1
                    gl(k, lc) = (2*k + 1)*oprim(lc)* &
                                (pmom(k, lc) - f)/(1.-f)
                end do

            end if

            flyr(lc) = f
            expbea(lc) = 0.0

        end do

! calculate slant optical depth

        if (umu0 < 0.0) then
            if (nid(0) < 0) then
                tausla(0) = largest
                tauslau(0) = largest
            else
                sum = 0.0
                sumu = 0.0
                do lc = 1, nid(0)
                    sum = sum + 2.*dtaucp(lc)*dsdh(0, lc)
                    sumu = sumu + 2.*dtauc(lc)*dsdh(0, lc)
                end do
                tausla(0) = sum
                tauslau(0) = sumu
            end if
        end if

        expbea(0) = exp(-tausla(0))

        do lc = 1, nlyr
            if (nid(lc) < 0) then
                tausla(lc) = largest
                tauslau(lc) = largest
            else
                sum = 0.0
                sumu = 0.0
                do lu = 1, min(nid(lc), lc)
                    sum = sum + dtaucp(lu)*dsdh(lc, lu)
                    sumu = sumu + dtauc(lu)*dsdh(lc, lu)
                end do
                do lu = min(nid(lc), lc) + 1, nid(lc)
                    sum = sum + 2.*dtaucp(lu)*dsdh(lc, lu)
                    sumu = sumu + 2.*dtauc(lu)*dsdh(lc, lu)
                end do
                tausla(lc) = sum
                tauslau(lc) = sumu
                if (tausla(lc) == tausla(lc - 1)) then
                    mu2(lc) = largest
                else
                    mu2(lc) = (taucpr(lc) - taucpr(lc - 1))/(tausla(lc) - tausla(lc - 1))
                    mu2(lc) = sign(AMAX1(abs(mu2(lc)), 1./largest), mu2(lc))
                end if
            end if
            expbea(lc) = exp(-tausla(lc))
        end do

!                      ** If no thermal emission, cut off medium below
!                         absorption optical depth = ABSCUT ( note that
!                         delta-M transformation leaves absorption
!                         optical depth invariant ).  Not worth the
!                         trouble for one-layer problems, though.
        lyrcut = .false.

        if (abstau >= abscut .and. .not. plank .and. ibcnd /= 1 .and. &
            nlyr > 1) lyrcut = .true.

        if (.not. lyrcut) ncut = nlyr

!                             ** Set arrays defining location of user
!                             ** output levels within delta-M-scaled
!                             ** computational mesh
        do lu = 1, ntau

            do lc = 1, nlyr

                if (utau(lu) >= tauc(lc - 1) .and. &
                    utau(lu) <= tauc(lc)) GO TO 60

            end do
            lc = nlyr

60          continue
            utaupr(lu) = utau(lu)
            if (deltam) utaupr(lu) = taucpr(lc - 1) + &
                                     (1.-ssalb(lc)*flyr(lc))*(utau(lu) - tauc(lc - 1))
            layru(lu) = lc

        end do
!                      ** Calculate computational polar angle cosines
!                         and associated quadrature weights for Gaussian
!                         quadrature on the interval (0,1) (upward)
        nn = nstr/2

        call qgausn(nn, cmu, cwt)
!                                  ** Downward (neg) angles and weights
        do iq = 1, nn
            cmu(iq + nn) = -cmu(iq)
            cwt(iq + nn) = cwt(iq)
        end do

!     IF( FBEAM.GT.0.0 ) THEN
!                               ** Compare beam angle to comput. angles
        do iq = 1, nn

!                      ** Dither mu2 if it is close to one of the
!                         quadrature angles.

            do lc = 1, nlyr
                if (abs(mu2(lc)) < 1.e5) then
                    if (abs(1.-abs(mu2(lc))/cmu(iq)) < 0.05) mu2(lc) = mu2(lc)*0.999
                end if
            end do

        end do

!     END IF

        if (.not. usrang .or. (onlyfl .and. maxumu >= nstr)) then

!                                   ** Set output polar angles to
!                                      computational polar angles
            numu = nstr

            do iu = 1, nn
                umu(iu) = -cmu(nn + 1 - iu)
            end do

            do iu = nn + 1, nstr
                umu(iu) = cmu(iu - nn)
            end do

        end if

        if (usrang .and. ibcnd == 1) then

!                               ** Shift positive user angle cosines to
!                                  upper locations and put negatives
!                                  in lower locations
            do iu = 1, numu
                umu(iu + numu) = umu(iu)
            end do

            do iu = 1, numu
                umu(iu) = -umu(2*numu + 1 - iu)
            end do

            numu = 2*numu

        end if

        if (.not. lyrcut .and. .not. lamber) then

            do k = 0, nstr
                hlpr(k) = (2*k + 1)*hl(k)
            end do

        end if

    end subroutine setdis

    subroutine setmtx(bdr, cband, cmu, cwt, delm0, dtaucp, gc, kk, &
                      lamber, lyrcut, mi, mi9m2, mxcmu, ncol, ncut, nnlyri, nn, nstr, taucpr, wk)

!        Calculate coefficient matrix for the set of equations
!        obtained from the boundary conditions and the continuity-
!        of-intensity-at-layer-interface equations;  store in the
!        special banded-matrix format required by LINPACK routines

!     I N P U T      V A R I A B L E S:

!       BDR      :  Surface bidirectional reflectivity
!       CMU      :  Abscissae for Gauss quadrature over angle cosine
!       CWT      :  Weights for Gauss quadrature over angle cosine
!       DELM0    :  Kronecker delta, delta-sub-m0
!       GC       :  Eigenvectors at polar quadrature angles, SC(1)
!       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
!       LYRCUT   :  Logical flag for truncation of comput. layer
!       NN       :  Number of streams in a hemisphere (NSTR/2)
!       NCUT     :  Total number of computational layers considered
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
!       (remainder are DISORT input variables)

!   O U T P U T     V A R I A B L E S:

!       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
!                      scaled by Eq. SC(12); in banded form required
!                      by LINPACK solution routines
!       NCOL     :  Counts of columns in CBAND

!   I N T E R N A L    V A R I A B L E S:

!       IROW     :  Points to row in CBAND
!       JCOL     :  Points to position in layer block
!       LDA      :  Row dimension of CBAND
!       NCD      :  Number of diagonals below or above main diagonal
!       NSHIFT   :  For positioning number of rows in band storage
!       WK       :  Temporary storage for EXP evaluations

!   Called by- DISORT, ALBTRN
!   Calls- ZEROIT
! +--------------------------------------------------------------------+
        implicit none

        integer, intent(IN)    :: mi
        integer, intent(IN)    :: mi9m2
        integer, intent(IN)    :: mxcmu
        integer, intent(IN)    :: ncut
        integer, intent(IN)    :: nnlyri
        integer, intent(IN)    :: nn
        integer, intent(IN)    :: nstr
        real, intent(IN)       :: bdr(mi, 0:mi)
        real, intent(IN)       :: cmu(mxcmu)
        real, intent(IN)       :: cwt(mxcmu)
        real, intent(IN)       :: delm0
        real, intent(IN)       :: dtaucp(*)
        real, intent(IN)       :: gc(mxcmu, mxcmu, *)
        real, intent(IN)       :: kk(mxcmu, *)
        real, intent(IN)       :: taucpr(0:*)
        logical, intent(IN)    :: lamber
        logical, intent(IN)    :: lyrcut
        integer, intent(OUT)   :: ncol
        real, intent(OUT)      :: wk(mxcmu)
        real, intent(OUT)      :: cband(mi9m2, nnlyri)

        integer :: iq, irow, jcol, jq, k, lc, lda, ncd, nncol, nshift
        real :: expa, sum

        intrinsic EXP

        call zeroit(cband, mi9m2*nnlyri)

        ncd = 3*nn - 1
        lda = 3*ncd + 1
        nshift = lda - 2*nstr + 1
        ncol = 0
!                         ** Use continuity conditions of Eq. STWJ(17)
!                            to form coefficient matrix in STWJ(20);
!                            employ scaling transformation STWJ(22)
        do lc = 1, ncut

            do iq = 1, nn
                wk(iq) = exp(kk(iq, lc)*dtaucp(lc))
            end do

            jcol = 0

            do iq = 1, nn

                ncol = ncol + 1
                irow = nshift - jcol

                do jq = 1, nstr
                    cband(irow + nstr, ncol) = gc(jq, iq, lc)
                    cband(irow, ncol) = -gc(jq, iq, lc)*wk(iq)
                    irow = irow + 1
                end do

                jcol = jcol + 1

            end do

            do iq = nn + 1, nstr

                ncol = ncol + 1
                irow = nshift - jcol

                do jq = 1, nstr
                    cband(irow + nstr, ncol) = gc(jq, iq, lc)*wk(nstr + 1 - iq)
                    cband(irow, ncol) = -gc(jq, iq, lc)
                    irow = irow + 1
                end do

                jcol = jcol + 1

            end do

        end do
!                  ** Use top boundary condition of STWJ(20a) for
!                     first layer

        jcol = 0

        do iq = 1, nn

            expa = exp(kk(iq, 1)*taucpr(1))
            irow = nshift - jcol + nn

            do jq = nn, 1, -1
                cband(irow, jcol + 1) = gc(jq, iq, 1)*expa
                irow = irow + 1
            end do

            jcol = jcol + 1

        end do

        do iq = nn + 1, nstr

            irow = nshift - jcol + nn

            do jq = nn, 1, -1
                cband(irow, jcol + 1) = gc(jq, iq, 1)
                irow = irow + 1
            end do

            jcol = jcol + 1

        end do
!                           ** Use bottom boundary condition of
!                              STWJ(20c) for last layer

        nncol = ncol - nstr
        jcol = 0

        do iq = 1, nn

            nncol = nncol + 1
            irow = nshift - jcol + nstr

            do jq = nn + 1, nstr

                if (lyrcut .or. (lamber .and. delm0 == 0)) then

!                          ** No azimuthal-dependent intensity if Lam-
!                             bert surface; no intensity component if
!                             truncated bottom layer

                    cband(irow, nncol) = gc(jq, iq, ncut)

                else

                    sum = 0.0

                    do k = 1, nn
                        sum = sum + cwt(k)*cmu(k)*bdr(jq - nn, k)* &
                              gc(nn + 1 - k, iq, ncut)
                    end do

                    cband(irow, nncol) = gc(jq, iq, ncut) - (1.+delm0)*sum
                end if

                irow = irow + 1

            end do

            jcol = jcol + 1

        end do

        do iq = nn + 1, nstr

            nncol = nncol + 1
            irow = nshift - jcol + nstr
            expa = wk(nstr + 1 - iq)

            do jq = nn + 1, nstr

                if (lyrcut .or. (lamber .and. delm0 == 0)) then

                    cband(irow, nncol) = gc(jq, iq, ncut)*expa

                else

                    sum = 0.0

                    do k = 1, nn
                        sum = sum + cwt(k)*cmu(k)*bdr(jq - nn, k)* &
                              gc(nn + 1 - k, iq, ncut)
                    end do

                    cband(irow, nncol) = (gc(jq, iq, ncut) - (1.+delm0)*sum)*expa
                end if

                irow = irow + 1

            end do

            jcol = jcol + 1

        end do

    end subroutine setmtx

    subroutine soleig(amb, apb, array, cmu, cwt, gl, mi, mazim, &
                      mxcmu, nn, nstr, ylmc, cc, evecc, eval, kk, gc, aad, eveccd, evald, wkd)

!         Solves eigenvalue/vector problem necessary to construct
!         homogeneous part of discrete ordinate solution; STWJ(8b)
!         ** NOTE ** Eigenvalue problem is degenerate when single
!                    scattering albedo = 1;  present way of doing it
!                    seems numerically more stable than alternative
!                    methods that we tried

!   I N P U T     V A R I A B L E S:

!       GL     :  Delta-M scaled Legendre coefficients of phase function
!                    (including factors 2l+1 and single-scatter albedo)
!       CMU    :  Computational polar angle cosines
!       CWT    :  Weights for quadrature over polar angle cosine
!       MAZIM  :  Order of azimuthal component
!       NN     :  Half the total number of streams
!       YLMC   :  Normalized associated Legendre polynomial
!                    at the quadrature angles CMU
!       (remainder are DISORT input variables)

!   O U T P U T    V A R I A B L E S:

!       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18)
!       EVAL   :  NN eigenvalues of Eq. SS(12) on return from ASYMTX
!                    but then square roots taken
!       EVECC  :  NN eigenvectors  (G+) - (G-)  on return
!                    from ASYMTX ( column j corresponds to EVAL(j) )
!                    but then  (G+) + (G-)  is calculated from SS(10),
!                    G+  and  G-  are separated, and  G+  is stacked on
!                    top of  G-  to form NSTR eigenvectors of SS(7)
!       GC     :  Permanent storage for all NSTR eigenvectors, but
!                    in an order corresponding to KK
!       KK     :  Permanent storage for all NSTR eigenvalues of SS(7),
!                    but re-ordered with negative values first ( square
!                    roots of EVAL taken and negatives added )

!   I N T E R N A L   V A R I A B L E S:

!       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced
!                    eigenvalue problem
!       ARRAY   :  Complete coefficient matrix of reduced eigenvalue
!                    problem: (alfa+beta)*(alfa-beta)
!       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11))
!       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11))
!       WKD     :  Scratch array required by ASYMTX

!   Called by- DISORT, ALBTRN
!   Calls- ASYMTX, ERRMSG
! +-------------------------------------------------------------------+
        implicit none

        integer, intent(IN) :: mi
        integer, intent(IN) :: mazim
        integer, intent(IN) :: mxcmu
        integer, intent(IN) :: nn
        integer, intent(IN) :: nstr
        real, intent(IN)    :: cmu(mxcmu)
        real, intent(IN)    :: cwt(mxcmu)
        real, intent(IN)    :: gl(0:mxcmu)
        real, intent(IN)    :: ylmc(0:mxcmu, mxcmu)
        real, intent(OUT)   :: amb(mi, mi)
        real, intent(OUT)   :: apb(mi, mi)
        real, intent(OUT)   :: array(mi, *)
        real, intent(OUT)   :: cc(mxcmu, mxcmu)
        real, intent(OUT)   :: evecc(mxcmu, mxcmu)
        real, intent(OUT)   :: eval(mi)
        real, intent(OUT)   :: kk(mxcmu)
        real, intent(OUT)   :: gc(mxcmu, mxcmu)
        double precision, intent(OUT) :: aad(mi, mi)
        double precision, intent(OUT) :: eveccd(mi, mi)
        double precision, intent(OUT) :: evald(mi)
        double precision, intent(OUT) :: wkd(mxcmu)

        integer :: ier, iq, jq, kq, l
        real :: alpha, beta, gpmigm, gpplgm, sum

        intrinsic ABS, SQRT

!                             ** Calculate quantities in Eqs. SS(5-6)
        do iq = 1, nn

            do jq = 1, nstr

                sum = 0.0
                do l = mazim, nstr - 1
                    sum = sum + gl(l)*ylmc(l, iq)*ylmc(l, jq)
                end do

                cc(iq, jq) = 0.5*sum*cwt(jq)

            end do

            do jq = 1, nn
!                             ** Fill remainder of array using symmetry
!                                relations  C(-mui,muj) = C(mui,-muj)
!                                and        C(-mui,-muj) = C(mui,muj)

                cc(iq + nn, jq) = cc(iq, jq + nn)
                cc(iq + nn, jq + nn) = cc(iq, jq)

!                                       ** Get factors of coeff. matrix
!                                          of reduced eigenvalue problem

                alpha = cc(iq, jq)/cmu(iq)
                beta = cc(iq, jq + nn)/cmu(iq)
                amb(iq, jq) = alpha - beta
                apb(iq, jq) = alpha + beta

            end do

            amb(iq, iq) = amb(iq, iq) - 1.0/cmu(iq)
            apb(iq, iq) = apb(iq, iq) - 1.0/cmu(iq)

        end do
!                      ** Finish calculation of coefficient matrix of
!                         reduced eigenvalue problem:  get matrix
!                         product (alfa+beta)*(alfa-beta); SS(12)
        do iq = 1, nn
            do jq = 1, nn
                sum = 0.
                do kq = 1, nn
                    sum = sum + apb(iq, kq)*amb(kq, jq)
                end do
                array(iq, jq) = sum
            end do
        end do
!                      ** Find (real) eigenvalues and eigenvectors

        call asymtx(array, evecc, eval, nn, mi, mxcmu, ier, wkd, aad, eveccd, evald)

        if (ier > 0) then

            write (*, FMT='(//,A,I4,A)') ' ASYMTX--eigenvalue no. ', &
                ier, '  didnt converge.  Lower-numbered eigenvalues wrong.'

            call errmsg('ASYMTX--convergence problems', .true.)

        end if

!DIR$ IVDEP
        do iq = 1, nn
            eval(iq) = sqrt(abs(eval(iq)))
            kk(iq + nn) = eval(iq)
!                                      ** Add negative eigenvalue
            kk(nn + 1 - iq) = -eval(iq)
        end do

!                          ** Find eigenvectors (G+) + (G-) from SS(10)
!                             and store temporarily in APB array
        do jq = 1, nn

            do iq = 1, nn

                sum = 0.
                do kq = 1, nn
                    sum = sum + amb(iq, kq)*evecc(kq, jq)
                end do

                apb(iq, jq) = sum/eval(jq)

            end do

        end do

        do jq = 1, nn
            do iq = 1, nn

                gpplgm = apb(iq, jq)
                gpmigm = evecc(iq, jq)
!                                ** Recover eigenvectors G+,G- from
!                                   their sum and difference; stack them
!                                   to get eigenvectors of full system
!                                   SS(7) (JQ = eigenvector number)

                evecc(iq, jq) = 0.5*(gpplgm + gpmigm)
                evecc(iq + nn, jq) = 0.5*(gpplgm - gpmigm)

!                                ** Eigenvectors corresponding to
!                                   negative eigenvalues (corresp. to
!                                   reversing sign of 'k' in SS(10) )
                gpplgm = -gpplgm
                evecc(iq, jq + nn) = 0.5*(gpplgm + gpmigm)
                evecc(iq + nn, jq + nn) = 0.5*(gpplgm - gpmigm)
                gc(iq + nn, jq + nn) = evecc(iq, jq)
                gc(nn + 1 - iq, jq + nn) = evecc(iq + nn, jq)
                gc(iq + nn, nn + 1 - jq) = evecc(iq, jq + nn)
                gc(nn + 1 - iq, nn + 1 - jq) = evecc(iq + nn, jq + nn)

            end do

        end do

    end subroutine soleig

    subroutine solve0(b, bdr, bem, bplank, cband, cmu, cwt, expbea, &
                      fbeam, fisot, ipvt, lamber, ll, lyrcut, mazim, &
                      mi, mi9m2, mxcmu, ncol, ncut, nn, nstr, nnlyri, &
                      pi, tplank, taucpr, umu0, z, zz, zplk0, zplk1)

!        Construct right-hand side vector B for general boundary
!        conditions STWJ(17) and solve system of equations obtained
!        from the boundary conditions and the continuity-of-
!        intensity-at-layer-interface equations.
!        Thermal emission contributes only in azimuthal independence.

!     I N P U T      V A R I A B L E S:

!       BDR      :  Surface bidirectional reflectivity
!       BEM      :  Surface bidirectional emissivity
!       BPLANK   :  Bottom boundary thermal emission
!       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
!                   scaled by Eq. SC(12); in banded form required
!                   by LINPACK solution routines
!       CMU      :  Abscissae for Gauss quadrature over angle cosine
!       CWT      :  Weights for Gauss quadrature over angle cosine
!       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
!       LYRCUT   :  Logical flag for truncation of comput. layer
!       MAZIM    :  Order of azimuthal component
!       ncol     :  Counts of columns in CBAND
!       NN       :  Order of double-Gauss quadrature (NSTR/2)
!       NCUT     :  Total number of computational layers considered
!       TPLANK   :  Top boundary thermal emission
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
!       ZZ       :  Beam source vectors in Eq. SS(19)
!       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
!       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
!       (remainder are DISORT input variables)

!   O U T P U T     V A R I A B L E S:

!       B        :  Right-hand side vector of Eq. SC(5) going into
!                   SGBSL; returns as solution vector of Eq. SC(12),
!                   constants of integration without exponential term

!      LL        :  Permanent storage for B, but re-ordered

!   I N T E R N A L    V A R I A B L E S:

!       IPVT     :  Integer vector of pivot indices
!       IT       :  Pointer for position in  B
!       NCD      :  Number of diagonals below or above main diagonal
!       RCOND    :  Indicator of singularity for CBAND
!       Z        :  Scratch array required by SGBCO

!   Called by- DISORT
!   Calls- ZEROIT, SGBCO, ERRMSG, SGBSL
! +-------------------------------------------------------------------+
        implicit none
        integer, intent(IN) :: mazim
        integer, intent(IN) :: mi
        integer, intent(IN) :: mi9m2
        integer, intent(IN) :: mxcmu
        integer, intent(IN) :: ncut
        integer, intent(IN) :: nn
        integer, intent(IN) :: nstr
        integer, intent(IN) :: nnlyri
        integer, intent(IN) :: ncol
        real, intent(IN)    :: bdr(mi, 0:mi)
        real, intent(IN)    :: bem(mi)
        real, intent(IN)    :: bplank
        real, intent(IN)    :: cmu(mxcmu)
        real, intent(IN)    :: cwt(mxcmu)
        real, intent(IN)    :: expbea(0:*)
        real, intent(IN)    :: fbeam
        real, intent(IN)    :: fisot
        real, intent(IN)    :: pi
        real, intent(IN)    :: tplank
        real, intent(IN)    :: taucpr(0:*)
        real, intent(IN)    :: umu0
        real, intent(IN)    :: zz(mxcmu, *)
        real, intent(IN)    :: zplk0(mxcmu, *)
        real, intent(IN)    :: zplk1(mxcmu, *)
        logical, intent(IN) :: lamber
        logical, intent(IN) :: lyrcut
        real, intent(INOUT) :: cband(mi9m2, nnlyri)
        integer, intent(OUT):: ipvt(*)
        real, intent(OUT)   :: ll(mxcmu, *)
        real, intent(OUT)   :: b(nnlyri)
        real, intent(OUT)   :: z(nnlyri)

        integer :: ipnt, iq, it, jq, lc, ncd
        real :: rcond, sum
!     ..

        call zeroit(b, nnlyri)
!                              ** Construct B,  STWJ(20a,c) for
!                                 parallel beam + bottom reflection +
!                                 thermal emission at top and/or bottom

        if (mazim > 0 .and. fbeam > 0.0) then

!                                         ** Azimuth-dependent case
!                                            (never called if FBEAM = 0)
            if (lyrcut .or. lamber) then

!               ** No azimuthal-dependent intensity for Lambert surface;
!                  no intensity component for truncated bottom layer

                do iq = 1, nn
!                                                  ** Top boundary
                    b(iq) = -zz(nn + 1 - iq, 1)*expbea(0)
!                                                  ** Bottom boundary

                    b(ncol - nn + iq) = -zz(iq + nn, ncut)*expbea(ncut)

                end do

            else

                do iq = 1, nn

                    b(iq) = -zz(nn + 1 - iq, 1)*expbea(0)

                    sum = 0.
                    do jq = 1, nn
                        sum = sum + cwt(jq)*cmu(jq)*bdr(iq, jq)* &
                              zz(nn + 1 - jq, ncut)*expbea(ncut)
                    end do

                    b(ncol - nn + iq) = sum
                    if (fbeam > 0.0) b(ncol - nn + iq) = sum + &
                                                         (bdr(iq, 0)*umu0*fbeam/pi - zz(iq + nn, ncut))*expbea(ncut)

                end do

            end if
!                             ** Continuity condition for layer
!                                interfaces of Eq. STWJ(20b)
            it = nn

            do lc = 1, ncut - 1

                do iq = 1, nstr
                    it = it + 1
                    b(it) = (zz(iq, lc + 1) - zz(iq, lc))*expbea(lc)
                end do

            end do

        else
!                                   ** Azimuth-independent case

            if (fbeam == 0.0) then

                do iq = 1, nn
!                                      ** Top boundary

                    b(iq) = -zplk0(nn + 1 - iq, 1) + fisot + tplank

                end do

                if (lyrcut) then
!                               ** No intensity component for truncated
!                                  bottom layer
                    do iq = 1, nn
!                                      ** Bottom boundary

                        b(ncol - nn + iq) = -zplk0(iq + nn, ncut) - &
                                            zplk1(iq + nn, ncut)*taucpr(ncut)
                    end do

                else

                    do iq = 1, nn

                        sum = 0.
                        do jq = 1, nn
                            sum = sum + cwt(jq)*cmu(jq)*bdr(iq, jq)* &
                                  (zplk0(nn + 1 - jq, ncut) + &
                                   zplk1(nn + 1 - jq, ncut)*taucpr(ncut))
                        end do

                        b(ncol - nn + iq) = 2.*sum + bem(iq)*bplank - &
                                            zplk0(iq + nn, ncut) - zplk1(iq + nn, ncut)* &
                                            taucpr(ncut)
                    end do

                end if
!                             ** Continuity condition for layer
!                                interfaces, STWJ(20b)
                it = nn
                do lc = 1, ncut - 1

                    do iq = 1, nstr
                        it = it + 1
                        b(it) = zplk0(iq, lc + 1) - zplk0(iq, lc) + &
                                (zplk1(iq, lc + 1) - zplk1(iq, lc))*taucpr(lc)
                    end do

                end do

            else

                do iq = 1, nn
                    b(iq) = -zz(nn + 1 - iq, 1)*expbea(0) - &
                            zplk0(nn + 1 - iq, 1) + fisot + tplank
                end do

                if (lyrcut) then

                    do iq = 1, nn
                        b(ncol - nn + iq) = -zz(iq + nn, ncut)*expbea(ncut) &
                                            - zplk0(iq + nn, ncut) - zplk1(iq + nn, ncut)*taucpr(ncut)
                    end do

                else

                    do iq = 1, nn

                        sum = 0.
                        do jq = 1, nn
                            sum = sum + cwt(jq)*cmu(jq)*bdr(iq, jq) &
                                  *(zz(nn + 1 - jq, ncut)*expbea(ncut) + zplk0(nn + 1 - jq, ncut) &
                                    + zplk1(nn + 1 - jq, ncut)*taucpr(ncut))
                        end do

                        b(ncol - nn + iq) = 2.*sum + (bdr(iq, 0)*umu0*fbeam/pi &
                                                      - zz(iq + nn, ncut))*expbea(ncut) + bem(iq)*bplank &
                                            - zplk0(iq + nn, ncut) - zplk1(iq + nn, ncut)*taucpr(ncut)
                    end do

                end if

                it = nn

                do lc = 1, ncut - 1

                    do iq = 1, nstr

                        it = it + 1
                        b(it) = (zz(iq, lc + 1) - zz(iq, lc))*expbea(lc) &
                                + zplk0(iq, lc + 1) - zplk0(iq, lc) + &
                                (zplk1(iq, lc + 1) - zplk1(iq, lc))*taucpr(lc)
                    end do

                end do

            end if

        end if
!                     ** Find L-U (lower/upper triangular) decomposition
!                        of band matrix CBAND and test if it is nearly
!                        singular (note: CBAND is destroyed)
!                        (CBAND is in LINPACK packed format)
        rcond = 0.0
        ncd = 3*nn - 1

        call sgbco(cband, mi9m2, ncol, ncd, ncd, ipvt, rcond, z)

        if (1.0 + rcond == 1.0) &
            call errmsg('SOLVE0--SGBCO says matrix near singular', .false.)

!                   ** Solve linear system with coeff matrix CBAND
!                      and R.H. side(s) B after CBAND has been L-U
!                      decomposed.  Solution is returned in B.

        call sgbsl(cband, mi9m2, ncol, ncd, ncd, ipvt, b, 0)

!                   ** Zero CBAND (it may contain 'foreign'
!                      elements upon returning from LINPACK);
!                      necessary to prevent errors

        call zeroit(cband, mi9m2*nnlyri)

        do lc = 1, ncut

            ipnt = lc*nstr - nn

            do iq = 1, nn
                ll(nn + 1 - iq, lc) = b(ipnt + 1 - iq)
                ll(iq + nn, lc) = b(iq + ipnt)
            end do

        end do

    end subroutine solve0

    subroutine surfac(albedo, delm0, fbeam, hlpr, lamber, mi, mazim, &
                      mxcmu, mxumu, nn, numu, nstr, onlyfl, umu, &
                      usrang, ylm0, ylmc, ylmu, bdr, emu, bem, rmu)

!       Specifies user's surface bidirectional properties, STWJ(21)

!   I N P U T     V A R I A B L E S:

!       DELM0  :  Kronecker delta, delta-sub-m0
!       HLPR   :  Legendre moments of surface bidirectional reflectivity
!                    (with 2K+1 factor included)
!       MAZIM  :  Order of azimuthal component
!       NN     :  Order of double-Gauss quadrature (NSTR/2)
!       YLM0   :  Normalized associated Legendre polynomial
!                 at the beam angle
!       YLMC   :  Normalized associated Legendre polynomials
!                 at the quadrature angles
!       YLMU   :  Normalized associated Legendre polynomials
!                 at the user angles
!       (remainder are DISORT input variables)

!    O U T P U T     V A R I A B L E S:

!       BDR :  Surface bidirectional reflectivity (computational angles)
!       RMU :  Surface bidirectional reflectivity (user angles)
!       BEM :  Surface directional emissivity (computational angles)
!       EMU :  Surface directional emissivity (user angles)

!    I N T E R N A L     V A R I A B L E S:

!       DREF      Directional reflectivity
!       NMUG   :  Number of angle cosine quadrature points on (0,1) for
!                   integrating bidirectional reflectivity to get
!                   directional emissivity (it is necessary to use a
!                   quadrature set distinct from the computational
!                   angles, because the computational angles may not be
!                   dense enough--NSTR may be too small--to give an
!                   accurate approximation for the integration).
!       GMU    :  The NMUG angle cosine quadrature points on (0,1)
!       GWT    :  The NMUG angle cosine quadrature weights on (0,1)
!       YLMG   :  Normalized associated Legendre polynomials
!                   at the NMUG quadrature angles

!   Called by- DISORT
!   Calls- QGAUSN, LEPOLY, ZEROIT, ERRMSG
! +-------------------------------------------------------------------+
        implicit none

        integer, intent(IN) :: mi
        integer, intent(IN) :: mazim
        integer, intent(IN) :: mxcmu
        integer, intent(IN) :: mxumu
        integer, intent(IN) :: nn
        integer, intent(IN) :: numu
        integer, intent(IN) :: nstr
        real, intent(IN)    :: albedo
        real, intent(IN)    :: delm0
        real, intent(IN)    :: fbeam
        real, intent(IN)    :: hlpr(0:mxcmu)
        real, intent(IN)    :: ylm0(0:mxcmu)
        real, intent(IN)    :: ylmc(0:mxcmu, mxcmu)
        real, intent(IN)    :: ylmu(0:mxcmu, mxumu)
        real, intent(IN)    :: umu(*)
        logical, intent(IN) :: lamber
        logical, intent(IN) :: onlyfl
        logical, intent(IN) :: usrang
        real, intent(OUT)   :: bdr(mi, 0:mi)
        real, intent(OUT)   :: emu(mxumu)
        real, intent(OUT)   :: bem(mi)
        real, intent(OUT)   :: rmu(mxumu, 0:mi)

        integer :: iq, iu, jg, jq, k
        real :: dref, sgn, sum

        if (pass4) then
            call qgausn(nmug, gmu, gwt)
            call lepoly(nmug, 0, maxstr, maxstr, gmu, ylmg)
! ** Convert Legendre polys. to negative GMU
            sgn = -1.0
            do k = 0, maxstr
                sgn = -sgn
                do jg = 1, nmug
                    ylmg(k, jg) = sgn*ylmg(k, jg)
                end do
            end do
            pass4 = .false.
        end if

        call zeroit(bdr, mi*(mi + 1))
        call zeroit(bem, mi)

        if (lamber .and. mazim == 0) then
            do iq = 1, nn
                bem(iq) = 1.-albedo
                do jq = 0, nn
                    bdr(iq, jq) = albedo
                end do
            end do
        else if (.not. lamber) then
            do iq = 1, nn
                do jq = 1, nn
                    sum = 0.0
                    do k = mazim, nstr - 1
                        sum = sum + hlpr(k)*ylmc(k, iq)*ylmc(k, jq + nn)
                    end do
                    bdr(iq, jq) = (2.-delm0)*sum
                end do
                if (fbeam > 0.0) then
                    sum = 0.0
                    do k = mazim, nstr - 1
                        sum = sum + hlpr(k)*ylmc(k, iq)*ylm0(k)
                    end do
                    bdr(iq, 0) = (2.-delm0)*sum
                end if
            end do

            if (mazim == 0) then

                if (nstr > maxstr) &
                    call errmsg('SURFAC--parameter MAXSTR too small', .true.)

!                              ** Integrate bidirectional reflectivity
!                                 at reflection polar angles CMU and
!                                 incident angles GMU to get
!                                 directional emissivity at
!                                 computational angles CMU.
                do iq = 1, nn

                    dref = 0.0

                    do jg = 1, nmug

                        sum = 0.0
                        do k = 0, nstr - 1
                            sum = sum + hlpr(k)*ylmc(k, iq)*ylmg(k, jg)
                        end do

                        dref = dref + 2.*gwt(jg)*gmu(jg)*sum

                    end do

                    bem(iq) = 1.-dref

                end do

            end if

        end if
!                                       ** Compute surface bidirectional
!                                          properties at user angles

        if (.not. onlyfl .and. usrang) then

            call zeroit(emu, mxumu)
            call zeroit(rmu, mxumu*(mi + 1))

            do iu = 1, numu

                if (umu(iu) > 0.0) then

                    if (lamber .and. mazim == 0) then

                        do iq = 0, nn
                            rmu(iu, iq) = albedo
                        end do

                        emu(iu) = 1.-albedo

                    else if (.not. lamber) then

                        do iq = 1, nn

                            sum = 0.0
                            do k = mazim, nstr - 1
                                sum = sum + hlpr(k)*ylmu(k, iu)*ylmc(k, iq + nn)
                            end do

                            rmu(iu, iq) = (2.-delm0)*sum

                        end do

                        if (fbeam > 0.0) then

                            sum = 0.0
                            do k = mazim, nstr - 1
                                sum = sum + hlpr(k)*ylmu(k, iu)*ylm0(k)
                            end do

                            rmu(iu, 0) = (2.-delm0)*sum

                        end if

                        if (mazim == 0) then

!                               ** Integrate bidirectional reflectivity
!                                  at reflection angles UMU and
!                                  incident angles GMU to get
!                                  directional emissivity at
!                                  user angles UMU.
                            dref = 0.0

                            do jg = 1, nmug

                                sum = 0.0
                                do k = 0, nstr - 1
                                    sum = sum + hlpr(k)*ylmu(k, iu)*ylmg(k, jg)
                                end do

                                dref = dref + 2.*gwt(jg)*gmu(jg)*sum

                            end do

                            emu(iu) = 1.-dref

                        end if

                    end if

                end if

            end do

        end if

    end subroutine surfac

!bm  SOLVEC calls SOLEIG and UPBEAM; if UPBEAM reports a potenially
!bm  unstable solution, the calculation is repeated with a slightly
!bm  changed single scattering albedo; this process is iterates
!bm  until a stable solution is found; as stable solutions may be
!bm  reached either by increasing or by decreasing the single
!bm  scattering albedo, both directions are explored ('upward' and
!bm  'downward' iteration); the solution which required the smaller
!bm  change in the single scattering albedo is finally returned
!bm  by SOLVEC.

    subroutine solvec(amb, apb, array, cmu, cwt, gl, mi, &
                      mazim, mxcmu, nn, nstr, ylm0, ylmc, cc, &
                      evecc, eval, kk, gc, aad, eveccd, evald, &
                      wk, wkd, delm0, fbeam, ipvt, pi, zj, zz, &
                      oprim, lc, dither, mu2, glsave, dgl)
        implicit none

!gy added glsave and dgl to call to allow adjustable dimensioning

        integer, intent(IN)  :: mi
        integer, intent(IN)  :: mazim
        integer, intent(IN)  :: mxcmu
        integer, intent(IN)  :: nn
        integer, intent(IN)  :: nstr
        integer, intent(IN)  :: lc
        real, intent(IN)     :: cmu(mxcmu)
        real, intent(IN)     :: cwt(mxcmu)
        real, intent(IN)     :: ylm0(0:mxcmu)
        real, intent(IN)     :: ylmc(0:mxcmu, mxcmu)
        real, intent(IN)     :: oprim
        real, intent(IN)     :: dither
        real, intent(IN)     :: mu2
        real, intent(IN)     :: delm0
        real, intent(IN)     :: fbeam
        real, intent(IN)     :: pi
        real, intent(INOUT)  :: gl(0:mxcmu)
        integer, intent(OUT) :: ipvt(*)
        real, intent(OUT)    :: evecc(mxcmu, mxcmu)
        real, intent(OUT)    :: eval(mi)
        real, intent(OUT)    :: kk(mxcmu)
        real, intent(OUT)    :: gc(mxcmu, mxcmu)
        real, intent(OUT)    :: amb(mi, mi)
        real, intent(OUT)    :: apb(mi, mi)
        real, intent(OUT)    :: array(mi, *)
        real, intent(OUT)    :: zj(mxcmu)
        real, intent(OUT)    :: zz(mxcmu)
        real, intent(OUT)    :: glsave(0:mxcmu)
        real, intent(OUT)    :: dgl(0:mxcmu)
        real, intent(OUT)    :: wk(mxcmu)
        real, intent(OUT)    :: cc(mxcmu, mxcmu)
        double precision, intent(OUT):: aad(mi, mi)
        double precision, intent(OUT):: eveccd(mi, mi)
        double precision, intent(OUT):: evald(mi)
        double precision, intent(OUT):: wkd(mxcmu)

!bm   Variables for instability fix

        integer :: uagain, dagain
        real :: minrcond, add, uadd, dadd, ssa, dssa, factor

        logical :: done, noup, nodn, debug, instab

        integer :: k

!bm   reset parameters

        done = .false.
        noup = .false.
        nodn = .false.

!bm   flag for printing debugging output
!      DEBUG  = .TRUE.
        debug = .false.

!bm   instability parameter; the solution is considered
!bm   unstable, if the RCOND reported by SGECO is smaller
!bm   than MINRCOND
        minrcond = 5000.*r1mach(4)

!bm   if an instability is detected, the single scattering albedo
!bm   is iterated downwards in steps of DADD and upwards in steps
!bm   of UADD; in practice, MINRCOND and -MINRCOND should
!bm   be reasonable choices for these parameters
        dadd = -minrcond
        uadd = minrcond

        uagain = 0
        dagain = 0
        add = dadd

!bm   save array GL( ) because it will be
!bm   changed if an iteration should be neccessary
        do k = mazim, nstr - 1
            glsave(k) = gl(k)
        end do

        ssa = oprim

!bm   in case of an instability reported by UPBEAM (INSTAB)
!bm   the single scattering albedo will be changed by a small
!bm   amount (ADD); this is indicated by DAGAIN or UAGAIN
!bm   being larger than 0; a change in the single scattering
!bm   albedo is equivalent to scaling the array GL( )

666     if (dagain > 0 .or. uagain > 0) then
            factor = (ssa + add)/ssa
            do k = mazim, nstr - 1
                gl(k) = gl(k)*factor
            end do

            ssa = ssa + add

!bm   if the single scattering albedo is now smaller than 0
!bm   the downward iteration is stopped and upward iteration
!bm   is forced instead

            if (ssa < dither) then
                nodn = .true.
                dagain = -1
                GO TO 778
            end if

!bm   if the single scattering albedo is now larger than its maximum
!bm   allowed value (1.0 - DITHER), the upward iteration is
!bm   stopped and downward iteration is forced instead

            if (ssa > 1.0 - dither) then
                noup = .true.
                uagain = -1
                GO TO 888
            end if
        end if

!     ** Solve eigenfunction problem in Eq. STWJ(8B);
!        return eigenvalues and eigenvectors

777     call soleig(amb, apb, array, cmu, cwt, gl, mi, &
                    mazim, mxcmu, nn, nstr, ylmc, cc, evecc, eval, kk, gc, aad, eveccd, evald, &
                    wkd)

!     ** Calculate particular solutions of
!        q.SS(18) for incident beam source

        if (fbeam > 0.0) then
            call upbeam(mu2, array, cc, cmu, delm0, fbeam, gl, &
                        ipvt, mazim, mxcmu, nn, nstr, pi, wk, &
                        ylm0, ylmc, zj, zz, minrcond, instab)
        end if

!     ** Calculate particular solutions of
!        Eq. SS(15) for thermal emission source
!        (not available in psndo.f)

!bm   finished if the result is stable on the first try
        if ((.not. instab) .and. (uagain == 0) .and. (dagain == 0)) then
            GO TO 999
        end if

!bm   downward iteration
        if (instab .and. uagain == 0) then
            dagain = dagain + 1
            GO TO 666
        end if

!bm   upward iteration
        if (instab .and. uagain > 0) then
            uagain = uagain + 1
            GO TO 666
        end if

!bm   ( DAGAIN .NE. 0 ) at this place means that the downward
!bm   iteration is finished

778     if (dagain /= 0 .and. uagain == 0) then

!bm   save downward iteration data for later use and
!bm   restore original input data
            do k = mazim, nstr - 1
                dgl(k) = gl(k)
                gl(k) = glsave(k)
            end do

            dssa = ssa
            ssa = oprim

!bm   start upward iteration
            add = uadd
            uagain = uagain + 1
            GO TO 666
        end if

!bm   both iterations finished
888     if (done) then
            GO TO 998
        end if

!bm  if neither upward nor downward iteration converged, the
!bm  original conditions are restored and SOLEIG/UPBEAM
!bm  is called for the last time

        if (noup .and. nodn) then

            do k = mazim, nstr - 1
                gl(k) = glsave(k)
            end do

            ssa = oprim

            if (debug) then
                write (*, *) '! *** Neither upward nor downward iteration'
                write (*, *) '! *** converged; using original result.'
            end if

            done = .true.
            GO TO 777
        end if

!bm  if upward iteration did not converge, the stable downward conditions
!bm  are restored and SOLEIG/UPBEAM is called for the last time
        if (noup) then
            do k = mazim, nstr - 1
                gl(k) = dgl(k)
            end do

            ssa = dssa

            if (debug) then
                write (*, *) '! *** The upward iteration did not converge.'
                write (*, *) '! *** Had to iterate ', dagain, &
                    ' times in layer LC =', lc, ';'
                write (*, *) '! *** changed SSA from ', oprim, ' to ', ssa, ','
                write (*, *) '! *** by a factor of ', ssa/oprim
            end if

            done = .true.
            GO TO 777
        end if

!bm  if downward iteration did not converge, we are done
!bm  (the result of the upward iteration will be used)
        if (nodn) then
            if (debug) then
                write (*, *) '! *** The downward iteration did not converge.'
                write (*, *) '! *** Had to iterate ', uagain, &
                    ' times in layer LC =', lc, ';'
                write (*, *) '! *** changed SSA from ', oprim, ' to ', ssa, ','
                write (*, *) '! *** by a factor of ', ssa/oprim
            end if

            done = .true.
            GO TO 998
        end if

!bm   if both iterations converged, and if the upward iteration
!bm   required more steps than the downward iteration, the stable
!bm   downward conditions are restored and SOLEIG/UPBEAM is
!bm   called for the last time

        if (uagain > dagain) then
            do k = mazim, nstr - 1
                gl(k) = dgl(k)
            end do

            ssa = dssa

            if (debug) then
                write (*, *) '! *** Both iterations converged;', ' using downward.'
                write (*, *) '! *** Had to iterate ', dagain, &
                    ' times in layer LC =', lc, ';'
                write (*, *) '! *** changed SSA from ', oprim, ' to ', ssa, ','
                write (*, *) '! *** by a factor of ', ssa/oprim
            end if

            done = .true.
            GO TO 777
        else

            if (debug) then
                write (*, *) '! *** Both iterations converged;', ' using upward.'
                write (*, *) '! *** Had to iterate ', uagain, &
                    ' times in layer LC =', lc, ';'
                write (*, *) '! *** changed SSA from ', oprim, ' to ', ssa, ','
                write (*, *) '! *** by a factor of ', ssa/oprim
            end if

            done = .true.
            GO TO 998
        end if

!bm   finally restore original input data
998     do k = mazim, nstr - 1
            gl(k) = glsave(k)
        end do

999     continue
    end subroutine solvec

    subroutine upbeam(mu2, array, cc, cmu, delm0, fbeam, gl, ipvt, mazim, &
                      mxcmu, nn, nstr, pi, wk, ylm0, ylmc, zj, zz, minrcond, instab)

!         Finds the incident-beam particular solution of SS(18)

!   I N P U T    V A R I A B L E S:

!       CC     :  C-sub-ij in Eq. SS(5)
!       CMU    :  Abscissae for Gauss quadrature over angle cosine
!       DELM0  :  Kronecker delta, delta-sub-m0
!       GL     :  Delta-M scaled Legendre coefficients of phase function
!                    (including factors 2L+1 and single-scatter albedo)
!       MAZIM  :  Order of azimuthal component
!       YLM0   :  Normalized associated Legendre polynomial
!                    at the beam angle
!       YLMC   :  Normalized associated Legendre polynomial
!                    at the quadrature angles
!       (remainder are DISORT input variables)

!   O U T P U T    V A R I A B L E S:

!       ZJ     :  Right-hand side vector X-sub-zero in SS(19); also the
!                 solution vector Z-sub-zero after solving that system

!       ZZ     :REAL, INTENT(IN)                         :: pi  Permanent storage for ZJ, but re-ordered

!   I N T E R N A L    V A R I A B L E S:

!       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19)
!       IPVT   :  Integer vector of pivot indices required by LINPACK
!       WK     :  Scratch array required by LINPACK

!   Called by- DISORT
!   Calls- SGECO, ERRMSG, SGESL
! +-------------------------------------------------------------------+
        implicit none
        integer, intent(IN)  :: mazim
        integer, intent(IN)  :: mxcmu
        integer, intent(IN)  :: nn
        integer, intent(IN)  :: nstr
        real, intent(IN)     :: mu2
        real, intent(IN)     :: cc(mxcmu, mxcmu)
        real, intent(IN)     :: cmu(mxcmu)
        real, intent(IN)     :: delm0
        real, intent(IN)     :: fbeam
        real, intent(IN)     :: gl(0:mxcmu)
        real, intent(IN)     :: pi
        real, intent(IN)     :: ylm0(0:mxcmu)
        real, intent(IN)     :: ylmc(0:mxcmu, *)
        real, intent(IN)     :: minrcond
        real, intent(INOUT)    :: array(mxcmu, mxcmu)
        integer, intent(OUT) :: ipvt(*)
        real, intent(OUT)    :: wk(mxcmu)
        real, intent(OUT)    :: zj(mxcmu)
        real, intent(OUT)    :: zz(mxcmu)
        logical, intent(OUT) :: instab

        integer :: iq, job, jq, k
        real :: rcond, sum

        do iq = 1, nstr

            do jq = 1, nstr
                array(iq, jq) = -cc(iq, jq)
            end do

            array(iq, iq) = 1.+cmu(iq)/mu2 + array(iq, iq)

            sum = 0.
            do k = mazim, nstr - 1
                sum = sum + gl(k)*ylmc(k, iq)*ylm0(k)
            end do

            zj(iq) = (2.-delm0)*fbeam*sum/(4.*pi)
        end do

!                  ** Find L-U (lower/upper triangular) decomposition
!                     of ARRAY and see if it is nearly singular
!                     (NOTE:  ARRAY is destroyed)
        rcond = 0.0

        call sgeco(array, mxcmu, nstr, ipvt, rcond, wk)

!bm      IF( 1.0 + RCOND.EQ.1.0 )
!bm     &    CALL ERRMSG('UPBEAM--SGECO says matrix near singular',.FALSE.)
!bm
!bm   replaced original check of RCOND by the following:

        instab = .false.
        if (abs(rcond) < minrcond) then
            instab = .true.
            return
        end if

!                ** Solve linear system with coeff matrix ARRAY
!                   (assumed already L-U decomposed) and R.H. side(s)
!                   ZJ;  return solution(s) in ZJ
        job = 0

        call sgesl(array, mxcmu, nstr, ipvt, zj, job)

!DIR$ IVDEP
        do iq = 1, nn
            zz(iq + nn) = zj(iq)
            zz(nn + 1 - iq) = zj(iq + nn)
        end do

    end subroutine upbeam

    subroutine zeroal(nd1, expbea, flyr, oprim, taucpr, xr0, xr1, &
                      nd2, cmu, cwt, psi, wk, z0, z1, zj, nd3, hlpr, ylm0, &
                      nd4, array, cc, evecc, nd5, gl, &
                      nd6, ylmc, nd7, ylmu, &
                      nd8, kk, ll, zz, zplk0, zplk1, nd9, gc, &
                      nd10, layru, utaupr, nd11, gu, &
                      nd12, z0u, z1u, zbeam, nd13, eval, &
                      nd14, amb, apb, nd15, ipvt, z, &
                      nd16, rfldir, rfldn, flup, uavg, dfdt, nd17, albmed, trnmed, &
                      nd18, u0u, nd19, uu)

!         ZERO ARRAYS; NDn is dimension of all arrays following
!         it in the argument list

!   Called by- DISORT
! --------------------------------------------------------------------
        implicit none

        integer, intent(IN) :: nd1
        integer, intent(IN) :: nd2
        integer, intent(IN) :: nd3
        integer, intent(IN) :: nd4
        integer, intent(IN) :: nd5
        integer, intent(IN) :: nd6
        integer, intent(IN) :: nd7
        integer, intent(IN) :: nd8
        integer, intent(IN) :: nd9
        integer, intent(IN) :: nd10
        integer, intent(IN) :: nd11
        integer, intent(IN) :: nd12
        integer, intent(IN) :: nd13
        integer, intent(IN) :: nd14
        integer, intent(IN) :: nd15
        integer, intent(IN) :: nd16
        integer, intent(IN) :: nd17
        integer, intent(IN) :: nd18
        integer, intent(IN) :: nd19
        integer, intent(OUT):: layru(*)
        integer, intent(OUT):: ipvt(*)
        real, intent(OUT)   :: expbea(*)
        real, intent(OUT)   :: flyr(*)
        real, intent(OUT)   :: oprim(*)
        real, intent(OUT)   :: taucpr(*)
        real, intent(OUT)   :: xr0(*)
        real, intent(OUT)   :: xr1(*)
        real, intent(OUT)   :: cmu(*)
        real, intent(OUT)   :: cwt(*)
        real, intent(OUT)   :: psi(*)
        real, intent(OUT)   :: wk(*)
        real, intent(OUT)   :: z0(*)
        real, intent(OUT)   :: z1(*)
        real, intent(OUT)   :: zj(*)
        real, intent(OUT)   :: hlpr(*)
        real, intent(OUT)   :: ylm0(*)
        real, intent(OUT)   :: array(*)
        real, intent(OUT)   :: cc(*)
        real, intent(OUT)   :: evecc(*)
        real, intent(OUT)   :: gl(*)
        real, intent(OUT)   :: ylmc(*)
        real, intent(OUT)   :: ylmu(*)
        real, intent(OUT)   :: kk(*)
        real, intent(OUT)   :: ll(*)
        real, intent(OUT)   :: zz(*)
        real, intent(OUT)   :: zplk0(*)
        real, intent(OUT)   :: zplk1(*)
        real, intent(OUT)   :: gc(*)
        real, intent(OUT)   :: utaupr(*)
        real, intent(OUT)   :: gu(*)
        real, intent(OUT)   :: z0u(*)
        real, intent(OUT)   :: z1u(*)
        real, intent(OUT)   :: zbeam(*)
        real, intent(OUT)   :: eval(*)
        real, intent(OUT)   :: amb(*)
        real, intent(OUT)   :: apb(*)
        real, intent(OUT)   :: z(*)
        real, intent(OUT)   :: rfldir(*)
        real, intent(OUT)   :: rfldn(*)
        real, intent(OUT)   :: flup(*)
        real, intent(OUT)   :: uavg(*)
        real, intent(OUT)   :: dfdt(*)
        real, intent(OUT)   :: albmed(*)
        real, intent(OUT)   :: trnmed(*)
        real, intent(OUT)   :: u0u(*)
        real, intent(OUT)   :: uu(*)

        integer :: n

        do n = 1, nd1
            expbea(n) = 0.0
            flyr(n) = 0.0
            oprim(n) = 0.0
            taucpr(n) = 0.0
            xr0(n) = 0.0
            xr1(n) = 0.0
        end do

        do n = 1, nd2
            cmu(n) = 0.0
            cwt(n) = 0.0
            psi(n) = 0.0
            wk(n) = 0.0
            z0(n) = 0.0
            z1(n) = 0.0
            zj(n) = 0.0
        end do

        do n = 1, nd3
            hlpr(n) = 0.0
            ylm0(n) = 0.0
        end do

        do n = 1, nd4
            array(n) = 0.0
            cc(n) = 0.0
            evecc(n) = 0.0
        end do

        do n = 1, nd5
            gl(n) = 0.0
        end do

        do n = 1, nd6
            ylmc(n) = 0.0
        end do

        do n = 1, nd7
            ylmu(n) = 0.0
        end do

        do n = 1, nd8
            kk(n) = 0.0
            ll(n) = 0.0
            zz(n) = 0.0
            zplk0(n) = 0.0
            zplk1(n) = 0.0
        end do

        do n = 1, nd9
            gc(n) = 0.0
        end do

        do n = 1, nd10
            layru(n) = 0
            utaupr(n) = 0.0
        end do

        do n = 1, nd11
            gu(n) = 0.0
        end do

        do n = 1, nd12
            z0u(n) = 0.0
            z1u(n) = 0.0
            zbeam(n) = 0.0
        end do

        do n = 1, nd13
            eval(n) = 0.0
        end do

        do n = 1, nd14
            amb(n) = 0.0
            apb(n) = 0.0
        end do

        do n = 1, nd15
            ipvt(n) = 0
            z(n) = 0.0
        end do

        do n = 1, nd16
            rfldir(n) = 0.
            rfldn(n) = 0.
            flup(n) = 0.
            uavg(n) = 0.
            dfdt(n) = 0.
        end do

        do n = 1, nd17
            albmed(n) = 0.
            trnmed(n) = 0.
        end do

        do n = 1, nd18
            u0u(n) = 0.
        end do

        do n = 1, nd19
            uu(n) = 0.
        end do

    end subroutine zeroal

    subroutine zeroit(a, length)

!         Zeros a real array A having LENGTH elements
! --------------------------------------------------------------------
        implicit none
        integer, intent(IN)  :: length
        real, intent(OUT)    :: a(length)

        integer :: l
!     ..

        do l = 1, length
            a(l) = 0.0
        end do

    end subroutine zeroit

    real function dref(mu, hl, nstr)

!        Exact flux albedo for given angle of incidence, given
!        a bidirectional reflectivity characterized by its
!        Legendre coefficients ( NOTE** these will only agree
!        with bottom-boundary albedos calculated by DISORT in
!        the limit as number of streams go to infinity, because
!        DISORT evaluates the integral 'CL' only approximately,
!        by quadrature, while this routine calculates it exactly.)

!  INPUT :   MU     Cosine of incidence angle
!            HL     Legendre coefficients of bidirectional reflectivity
!          NSTR     Number of elements of HL to consider

!  INTERNAL VARIABLES (P-sub-L is the L-th Legendre polynomial) :

!       CL      Integral from 0 to 1 of  MU * P-sub-L(MU)
!                   (vanishes for  L = 3, 5, 7, ... )
!       PL      P-sub-L
!       PLM1    P-sub-(L-1)
!       PLM2    P-sub-(L-2)

!   Called by- CHEKIN
!   Calls- ERRMSG
! +-------------------------------------------------------------------+
        implicit none
        integer, intent(IN)                      :: nstr
        real, intent(IN)                         :: mu
        real, intent(IN)                         :: hl(0:nstr)

        integer :: l
        real :: cl, pl, plm1, plm2

        if (pass5) then
            pass5 = .false.
            cl = 0.125
            c(2) = 10.*cl
            do l = 4, maxtrm, 2
                cl = -cl*(l - 3)/(l + 2)
                c(l) = 2.*(2*l + 1)*cl
            end do
        end if

        if (nstr < 2 .or. abs(mu) > 1.0) &
            call errmsg('DREF--input argument error(s)', .true.)

        if (nstr > maxtrm) call errmsg('DREF--parameter MAXTRM too small', .true.)

        dref = hl(0) - 2.*hl(1)*mu
        plm2 = 1.0
        plm1 = -mu

        do l = 2, nstr - 1
!                                ** Legendre polynomial recurrence

            pl = ((2*l - 1)*(-mu)*plm1 - (l - 1)*plm2)/l

            if (mod(l, 2) == 0) dref = dref + c(l)*hl(l)*pl

            plm2 = plm1
            plm1 = pl

        end do

        if (dref < 0.0 .or. dref > 1.0) &
            call errmsg('DREF--albedo value not in (0,1)', .false.)

    end function dref

    real function ratio(a, b)

!        Calculate ratio  A/B  with over- and under-flow protection
!        (thanks to Prof. Jeff Dozier for some suggestions here).
!        Since this routine takes two logs, it is no speed demon,
!        but it is invaluable for comparing results from two runs
!        of a program under development.

!        NOTE:  In Fortran90, built-in functions TINY and HUGE
!               can replace the R1MACH calls.
! ---------------------------------------------------------------
        implicit none
        real, intent(IN)                     :: a
        real, intent(IN)                     :: b

        real :: absa, absb, powa, powb

        intrinsic ABS, LOG10, SIGN
!     ..

        if (pass6) then

            tinyVar = r1mach(1)
            hugeVar = r1mach(2)
            powmax = log10(hugeVar)
            powmin = log10(tinyVar)
            pass6 = .false.

        end if

        if (a == 0.0) then

            if (b == 0.0) then

                ratio = 1.0

            else

                ratio = 0.0

            end if

        else if (b == 0.0) then

            ratio = sign(hugeVar, a)

        else

            absa = abs(a)
            absb = abs(b)
            powa = log10(absa)
            powb = log10(absb)

            if (absa < tinyVar .and. absb < tinyVar) then

                ratio = 1.0

            else if (powa - powb >= powmax) then

                ratio = hugeVar

            else if (powa - powb <= powmin) then

                ratio = tinyVar

            else

                ratio = absa/absb

            end if
!                      ** DONT use old trick of determining sign
!                      ** from A*B because A*B may (over/under)flow

            if ((a > 0.0 .and. b < 0.0) .or. &
                (a < 0.0 .and. b > 0.0)) ratio = -ratio

        end if

    end function ratio

    subroutine errmsg(messag, fatal)

!        Print out a warning or error message;  abort if error
!        after making symbolic dump (machine-specific)
        implicit none

        character(LEN=*), intent(IN)   :: messag
        logical, intent(IN)             :: fatal
!LOGICAL :: cray

        if (fatal) then
            write (*, '(//,2A,//)') ' ******* ERROR >>>>>>  ', messag
            stop
        end if

        nummsg = nummsg + 1
        if (msglim) return

        if (nummsg <= maxmsg) then
            write (*, '(/,2A,/)') ' ******* WARNING >>>>>>  ', messag
        else
            write (*, 99)
            msglim = .true.
        end if

        return

99      format(//, ' >>>>>>  TOO MANY WARNING MESSAGES --  ', &
                'They will no longer be printed  <<<<<<<', //)
    end subroutine errmsg

    logical function wrtbad(varnam)
!          Write names of erroneous variables and return 'TRUE'

!      INPUT :   VarNam = Name of erroneous variable to be written
!                         ( CHARACTER, any length )
        implicit none

        character(LEN=*), intent(IN)        :: varnam

        wrtbad = .true.
        nummsg2 = nummsg2 + 1
        write (*, '(3A)') ' ****  Input variable  ', varnam, '  in error  ****'
        if (nummsg2 == maxmsg2) &
            call errmsg('Too many input errors.  Aborting...', .true.)

    end function wrtbad

    logical function wrtdim(dimnam, minval)

!          Write name of too-small symbolic dimension and
!          the value it should be increased to;  return 'TRUE'

!      INPUT :  DimNam = Name of symbolic dimension which is too small
!                        ( CHARACTER, any length )
!               Minval = Value to which that dimension should be
!                        increased (at least)
        implicit none

        character(LEN=*), intent(IN)        :: dimnam
        integer, intent(IN)                  :: minval

        write (*, '(3A,I7)') ' ****  Symbolic dimension  ', dimnam, &
            '  should be increased to at least ', minval

        wrtdim = .true.
    end function wrtdim

    logical function tstbad(varnam, relerr)

!       Write name (VarNam) of variable failing self-test and its
!       percent error from the correct value;  return  'FALSE'.
        implicit none

        character(LEN=*), intent(IN)        :: varnam
        real, intent(IN)                     :: relerr

        tstbad = .false.
        write (*, '(/,3A,1P,E11.2,A)') &
            ' Output variable ', varnam, ' differed by ', 100.*relerr, &
            ' per cent from correct value.  Self-test failed.'

    end function tstbad

    subroutine sgbco(abd, lda, n, ml, mu, ipvt, rcond, z)

!         FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION
!         AND ESTIMATES THE CONDITION OF THE MATRIX.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.

!     INPUT:

!        ABD     REAL(LDA, N)
!                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
!                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
!                SEE THE COMMENTS BELOW FOR DETAILS.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!                LDA MUST BE .GE. 2*ML + MU + 1 .

!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.

!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!                0 .LE. ML .LT. N .

!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!                0 .LE. MU .LT. N .
!                MORE EFFICIENT IF  ML .LE. MU .

!     ON RETURN

!        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.

!        RCOND   REAL
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.

!        Z       REAL(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

!     BAND STORAGE

!           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
!           WILL SET UP THE INPUT.

!                   ML = (BAND WIDTH BELOW THE DIAGONAL)
!                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-MU)
!                      I2 = MIN0(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE

!           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
!           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
!           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
!           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
!           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
!           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.

!     EXAMPLE:  IF THE ORIGINAL MATRIX IS

!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66

!      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN

!            *  *  *  +  +  +  , * = NOT USED
!            *  * 13 24 35 46  , + = USED FOR PIVOTING
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *

!     ROUTINES CALLED:  FROM LINPACK: SGBFA
!                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
!                       FROM FORTRAN: ABS, AMAX1, MAX0, MIN0, SIGN

        implicit none

        integer, intent(IN) :: lda
        integer, intent(IN) :: n
        integer, intent(IN) :: ml
        integer, intent(IN) :: mu
        real, intent(INOUT)    :: abd(lda, *)
        integer, intent(OUT) :: ipvt(*)
        real, intent(OUT)   :: rcond
        real, intent(OUT)   :: z(*)

        real :: ek, t, wk, wkm
        real :: anorm, s, sm, ynorm
        integer :: is, info, j, ju, k, kb, kp1, l, la, lm, lz, m, mm

!                       ** COMPUTE 1-NORM OF A
        anorm = 0.0e0
        l = ml + 1
        is = l + mu
        do j = 1, n
            anorm = AMAX1(anorm, sasum(l, abd(is, j), 1))
            if (is > ml + 1) is = is - 1
            if (j <= mu) l = l + 1
            if (j >= n - ml) l = l - 1
        end do
!                                               ** FACTOR
        call sgbfa(abd, lda, n, ml, mu, ipvt, info)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!                     ** SOLVE TRANS(U)*W = E
        ek = 1.0e0
        do j = 1, n
            z(j) = 0.0e0
        end do

        m = ml + mu + 1
        ju = 0
        do k = 1, n
            if (z(k) /= 0.0e0) ek = sign(ek, -z(k))
            if (abs(ek - z(k)) > abs(abd(m, k))) then
                s = abs(abd(m, k))/abs(ek - z(k))
                call sscal(n, s, z, 1)
                ek = s*ek
            end if
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = abs(wk)
            sm = abs(wkm)
            if (abd(m, k) /= 0.0e0) then
                wk = wk/abd(m, k)
                wkm = wkm/abd(m, k)
            else
                wk = 1.0e0
                wkm = 1.0e0
            end if
            kp1 = k + 1
            ju = MIN0(MAX0(ju, mu + ipvt(k)), n)
            mm = m
            if (kp1 <= ju) then
                do j = kp1, ju
                    mm = mm - 1
                    sm = sm + abs(z(j) + wkm*abd(mm, j))
                    z(j) = z(j) + wk*abd(mm, j)
                    s = s + abs(z(j))
                end do
                if (s < sm) then
                    t = wkm - wk
                    wk = wkm
                    mm = m
                    do j = kp1, ju
                        mm = mm - 1
                        z(j) = z(j) + t*abd(mm, j)
                    end do
                end if
            end if
            z(k) = wk
        end do

        s = 1.0e0/sasum(n, z, 1)
        call sscal(n, s, z, 1)

!                         ** SOLVE TRANS(L)*Y = W
        do kb = 1, n
            k = n + 1 - kb
            lm = MIN0(ml, n - k)
            if (k < n) z(k) = z(k) + sdot(lm, abd(m + 1, k), 1, z(k + 1), 1)
            if (abs(z(k)) > 1.0e0) then
                s = 1.0e0/abs(z(k))
                call sscal(n, s, z, 1)
            end if
            l = ipvt(k)
            t = z(l)
            z(l) = z(k)
            z(k) = t
        end do

        s = 1.0e0/sasum(n, z, 1)
        call sscal(n, s, z, 1)

        ynorm = 1.0e0
!                         ** SOLVE L*V = Y
        do k = 1, n
            l = ipvt(k)
            t = z(l)
            z(l) = z(k)
            z(k) = t
            lm = MIN0(ml, n - k)
            if (k < n) call saxpy(lm, t, abd(m + 1, k), 1, z(k + 1), 1)
            if (abs(z(k)) > 1.0e0) then
                s = 1.0e0/abs(z(k))
                call sscal(n, s, z, 1)
                ynorm = s*ynorm
            end if
        end do

        s = 1.0e0/sasum(n, z, 1)
        call sscal(n, s, z, 1)
        ynorm = s*ynorm
!                           ** SOLVE  U*Z = W
        do kb = 1, n
            k = n + 1 - kb
            if (abs(z(k)) > abs(abd(m, k))) then
                s = abs(abd(m, k))/abs(z(k))
                call sscal(n, s, z, 1)
                ynorm = s*ynorm
            end if
            if (abd(m, k) /= 0.0e0) z(k) = z(k)/abd(m, k)
            if (abd(m, k) == 0.0e0) z(k) = 1.0e0
            lm = MIN0(k, m) - 1
            la = m - lm
            lz = k - lm
            t = -z(k)
            call saxpy(lm, t, abd(la, k), 1, z(lz), 1)
        end do
!                              ** MAKE ZNORM = 1.0
        s = 1.0e0/sasum(n, z, 1)
        call sscal(n, s, z, 1)
        ynorm = s*ynorm

        if (anorm /= 0.0e0) rcond = ynorm/anorm
        if (anorm == 0.0e0) rcond = 0.0e0

    end subroutine sgbco

    subroutine sgbfa(abd, lda, n, ml, mu, ipvt, info)

!         FACTORS A REAL BAND MATRIX BY ELIMINATION.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.

!     INPUT:  SAME AS 'SGBCO'

!     ON RETURN:

!        ABD,IPVT    SAME AS 'SGBCO'

!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
!                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.

!     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
!                       FROM FORTRAN: MAX0, MIN0
        implicit none

        integer, intent(IN)   :: n
        integer, intent(IN)   :: ml
        integer, intent(IN)   :: mu
        integer, intent(IN)   :: lda
        real, intent(INOUT)     :: abd(lda, *)
        integer, intent(OUT)  :: ipvt(*)
        integer, intent(OUT)  :: info

        real :: t
        integer :: i, i0, j, ju, jz, j0, j1, k, kp1, l, lm, m, mm, nm1

        m = ml + mu + 1
        info = 0
!                        ** ZERO INITIAL FILL-IN COLUMNS
        j0 = mu + 2
        j1 = MIN0(n, m) - 1
        do jz = j0, j1
            i0 = m + 1 - jz
            do i = i0, ml
                abd(i, jz) = 0.0e0
            end do
        end do
        jz = j1
        ju = 0

!                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
        nm1 = n - 1
        do k = 1, nm1
            kp1 = k + 1
!                                  ** ZERO NEXT FILL-IN COLUMN
            jz = jz + 1
            if (jz <= n) then
                do i = 1, ml
                    abd(i, jz) = 0.0e0
                end do
            end if
!                                  ** FIND L = PIVOT INDEX
            lm = MIN0(ml, n - k)
            l = isamax(lm + 1, abd(m, k), 1) + m - 1
            ipvt(k) = l + k - m

            if (abd(l, k) == 0.0e0) then
!                                      ** ZERO PIVOT IMPLIES THIS COLUMN
!                                      ** ALREADY TRIANGULARIZED
                info = k
            else
!                                ** INTERCHANGE IF NECESSARY
                if (l /= m) then
                    t = abd(l, k)
                    abd(l, k) = abd(m, k)
                    abd(m, k) = t
                end if
!                                   ** COMPUTE MULTIPLIERS
                t = -1.0e0/abd(m, k)
                call sscal(lm, t, abd(m + 1, k), 1)

!                               ** ROW ELIMINATION WITH COLUMN INDEXING

                ju = MIN0(MAX0(ju, mu + ipvt(k)), n)
                mm = m
                do j = kp1, ju
                    l = l - 1
                    mm = mm - 1
                    t = abd(l, j)
                    if (l /= mm) then
                        abd(l, j) = abd(mm, j)
                        abd(mm, j) = t
                    end if
                    call saxpy(lm, t, abd(m + 1, k), 1, abd(mm + 1, j), 1)
                end do

            end if

        end do

        ipvt(n) = n
        if (abd(m, n) == 0.0e0) info = n

    end subroutine sgbfa

    subroutine sgbsl(abd, lda, n, ml, mu, ipvt, b, job)

!         SOLVES THE REAL BAND SYSTEM
!            A * X = B  OR  TRANSPOSE(A) * X = B
!         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     INPUT:

!        ABD     REAL(LDA, N)
!                THE OUTPUT FROM SBGCO OR SGBFA.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .

!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.

!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.

!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.

!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SBGCO OR SGBFA.

!        B       REAL(N)
!                THE RIGHT HAND SIDE VECTOR.

!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
!                            TRANS(A)  IS THE TRANSPOSE.

!     ON RETURN

!        B       THE SOLUTION VECTOR  X .

!     ERROR CONDITION

!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
!        OR SGBFA HAS SET INFO .EQ. 0 .

!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
!                       FROM FORTRAN: MIN0
        implicit none

        integer, intent(IN)  :: lda
        integer, intent(IN)  :: n
        integer, intent(IN)  :: ml
        integer, intent(IN)  :: mu
        integer, intent(IN)  :: job
        integer, intent(IN)  :: ipvt(*)
        real, intent(IN)     :: abd(lda, *)
        real, intent(IN OUT) :: b(*)

        real :: t
        integer :: k, kb, l, la, lb, lm, m, nm1

        m = mu + ml + 1
        nm1 = n - 1
        if (job == 0) then
!                               ** JOB = 0 , SOLVE  A * X = B
!                               ** FIRST SOLVE L*Y = B
            if (ml /= 0) then
                do k = 1, nm1
                    lm = MIN0(ml, n - k)
                    l = ipvt(k)
                    t = b(l)
                    if (l /= k) then
                        b(l) = b(k)
                        b(k) = t
                    end if
                    call saxpy(lm, t, abd(m + 1, k), 1, b(k + 1), 1)
                end do
            end if
!                           ** NOW SOLVE  U*X = Y
            do kb = 1, n
                k = n + 1 - kb
                b(k) = b(k)/abd(m, k)
                lm = MIN0(k, m) - 1
                la = m - lm
                lb = k - lm
                t = -b(k)
                call saxpy(lm, t, abd(la, k), 1, b(lb), 1)
            end do

        else
!                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
!                                  ** FIRST SOLVE  TRANS(U)*Y = B
            do k = 1, n
                lm = MIN0(k, m) - 1
                la = m - lm
                lb = k - lm
                t = sdot(lm, abd(la, k), 1, b(lb), 1)
                b(k) = (b(k) - t)/abd(m, k)
            end do
!                                  ** NOW SOLVE TRANS(L)*X = Y
            if (ml /= 0) then
                do kb = 1, nm1
                    k = n - kb
                    lm = MIN0(ml, n - k)
                    b(k) = b(k) + sdot(lm, abd(m + 1, k), 1, b(k + 1), 1)
                    l = ipvt(k)
                    if (l /= k) then
                        t = b(l)
                        b(l) = b(k)
                        b(k) = t
                    end if
                end do
            end if

        end if

    end subroutine sgbsl

    subroutine sgeco(a, lda, n, ipvt, rcond, z)

!         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
!         AND ESTIMATES THE CONDITION OF THE MATRIX.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
!         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!     ON RETURN

!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.

!        RCOND   REAL
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.

!        Z       REAL(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

!     ROUTINES CALLED:  FROM LINPACK: SGEFA
!                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
!                       FROM FORTRAN: ABS, AMAX1, SIGN
        implicit none

        integer, intent(IN)   :: lda
        integer, intent(IN)   :: n
        real, intent(INOUT)   :: a(lda, *)
        integer, intent(OUT)  :: ipvt(*)
        real, intent(OUT)     :: rcond
        real, intent(OUT)     :: z(*)

        real :: ek, t, wk, wkm
        real :: anorm, s, sm, ynorm
        integer :: info, j, k, kb, kp1, l

!                        ** COMPUTE 1-NORM OF A
        anorm = 0.0e0
        do j = 1, n
            anorm = AMAX1(anorm, sasum(n, a(1, j), 1))
        end do
!                                      ** FACTOR
        call sgefa(a, lda, n, ipvt, info)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!                        ** SOLVE TRANS(U)*W = E
        ek = 1.0e0
        do j = 1, n
            z(j) = 0.0e0
        end do

        do k = 1, n
            if (z(k) /= 0.0e0) ek = sign(ek, -z(k))
            if (abs(ek - z(k)) > abs(a(k, k))) then
                s = abs(a(k, k))/abs(ek - z(k))
                call sscal(n, s, z, 1)
                ek = s*ek
            end if
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = abs(wk)
            sm = abs(wkm)
            if (a(k, k) /= 0.0e0) then
                wk = wk/a(k, k)
                wkm = wkm/a(k, k)
            else
                wk = 1.0e0
                wkm = 1.0e0
            end if
            kp1 = k + 1
            if (kp1 <= n) then
                do j = kp1, n
                    sm = sm + abs(z(j) + wkm*a(k, j))
                    z(j) = z(j) + wk*a(k, j)
                    s = s + abs(z(j))
                end do
                if (s < sm) then
                    t = wkm - wk
                    wk = wkm
                    do j = kp1, n
                        z(j) = z(j) + t*a(k, j)
                    end do
                end if
            end if
            z(k) = wk
        end do

        s = 1.0e0/sasum(n, z, 1)
        call sscal(n, s, z, 1)
!                                ** SOLVE TRANS(L)*Y = W
        do kb = 1, n
            k = n + 1 - kb
            if (k < n) z(k) = z(k) + sdot(n - k, a(k + 1, k), 1, z(k + 1), 1)
            if (abs(z(k)) > 1.0e0) then
                s = 1.0e0/abs(z(k))
                call sscal(n, s, z, 1)
            end if
            l = ipvt(k)
            t = z(l)
            z(l) = z(k)
            z(k) = t
        end do

        s = 1.0e0/sasum(n, z, 1)
        call sscal(n, s, z, 1)
!                                 ** SOLVE L*V = Y
        ynorm = 1.0e0
        do k = 1, n
            l = ipvt(k)
            t = z(l)
            z(l) = z(k)
            z(k) = t
            if (k < n) call saxpy(n - k, t, a(k + 1, k), 1, z(k + 1), 1)
            if (abs(z(k)) > 1.0e0) then
                s = 1.0e0/abs(z(k))
                call sscal(n, s, z, 1)
                ynorm = s*ynorm
            end if
        end do

        s = 1.0e0/sasum(n, z, 1)
        call sscal(n, s, z, 1)
!                                  ** SOLVE  U*Z = V
        ynorm = s*ynorm
        do kb = 1, n
            k = n + 1 - kb
            if (abs(z(k)) > abs(a(k, k))) then
                s = abs(a(k, k))/abs(z(k))
                call sscal(n, s, z, 1)
                ynorm = s*ynorm
            end if
            if (a(k, k) /= 0.0e0) z(k) = z(k)/a(k, k)
            if (a(k, k) == 0.0e0) z(k) = 1.0e0
            t = -z(k)
            call saxpy(k - 1, t, a(1, k), 1, z(1), 1)
        end do
!                                   ** MAKE ZNORM = 1.0
        s = 1.0e0/sasum(n, z, 1)
        call sscal(n, s, z, 1)
        ynorm = s*ynorm

        if (anorm /= 0.0e0) rcond = ynorm/anorm
        if (anorm == 0.0e0) rcond = 0.0e0

    end subroutine sgeco

    subroutine sgefa(a, lda, n, ipvt, info)

!         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .

!     INPUT:  SAME AS 'SGECO'

!     ON RETURN:

!        A,IPVT  SAME AS 'SGECO'

!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
        implicit none

        integer, intent(IN)  :: lda
        integer, intent(IN)  :: n
        real, intent(IN OUT) :: a(lda, *)
        integer, intent(OUT) :: ipvt(*)
        integer, intent(OUT) :: info

        real :: t
        integer :: j, k, kp1, l, nm1

!                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
        info = 0
        nm1 = n - 1
        do k = 1, nm1
            kp1 = k + 1
!                                            ** FIND L = PIVOT INDEX
            l = isamax(n - k + 1, a(k, k), 1) + k - 1
            ipvt(k) = l

            if (a(l, k) == 0.0e0) then
!                                     ** ZERO PIVOT IMPLIES THIS COLUMN
!                                     ** ALREADY TRIANGULARIZED
                info = k
            else
!                                     ** INTERCHANGE IF NECESSARY
                if (l /= k) then
                    t = a(l, k)
                    a(l, k) = a(k, k)
                    a(k, k) = t
                end if
!                                     ** COMPUTE MULTIPLIERS
                t = -1.0e0/a(k, k)
                call sscal(n - k, t, a(k + 1, k), 1)

!                              ** ROW ELIMINATION WITH COLUMN INDEXING
                do j = kp1, n
                    t = a(l, j)
                    if (l /= k) then
                        a(l, j) = a(k, j)
                        a(k, j) = t
                    end if
                    call saxpy(n - k, t, a(k + 1, k), 1, a(k + 1, j), 1)
                end do

            end if

        end do

        ipvt(n) = n
        if (a(n, n) == 0.0e0) info = n

    end subroutine sgefa

    subroutine sgesl(a, lda, n, ipvt, b, job)

!         SOLVES THE REAL SYSTEM
!            A * X = B  OR  TRANS(A) * X = B
!         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE OUTPUT FROM SGECO OR SGEFA.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SGECO OR SGEFA.

!        B       REAL(N)
!                THE RIGHT HAND SIDE VECTOR.

!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
!                            TRANS(A)  IS THE TRANSPOSE.

!     ON RETURN

!        B       THE SOLUTION VECTOR  X .

!     ERROR CONDITION

!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
!        OR SGEFA HAS SET INFO .EQ. 0 .

!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
        implicit none

        integer, intent(IN)  :: lda
        integer, intent(IN)  :: n
        integer, intent(IN)  :: job
        integer, intent(IN)  :: ipvt(*)
        real, intent(IN OUT) :: b(*)
        real, intent(IN)     :: a(lda, *)

        real :: t
        integer :: k, kb, l, nm1

        nm1 = n - 1
        if (job == 0) then
!                                 ** JOB = 0 , SOLVE  A * X = B
!                                     ** FIRST SOLVE  L*Y = B
            do k = 1, nm1
                l = ipvt(k)
                t = b(l)
                if (l /= k) then
                    b(l) = b(k)
                    b(k) = t
                end if
                call saxpy(n - k, t, a(k + 1, k), 1, b(k + 1), 1)
            end do
!                                    ** NOW SOLVE  U*X = Y
            do kb = 1, n
                k = n + 1 - kb
                b(k) = b(k)/a(k, k)
                t = -b(k)
                call saxpy(k - 1, t, a(1, k), 1, b(1), 1)
            end do

        else
!                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
!                                    ** FIRST SOLVE  TRANS(U)*Y = B
            do k = 1, n
                t = sdot(k - 1, a(1, k), 1, b(1), 1)
                b(k) = (b(k) - t)/a(k, k)
            end do
!                                    ** NOW SOLVE  TRANS(L)*X = Y
            do kb = 1, nm1
                k = n - kb
                b(k) = b(k) + sdot(n - k, a(k + 1, k), 1, b(k + 1), 1)
                l = ipvt(k)
                if (l /= k) then
                    t = b(l)
                    b(l) = b(k)
                    b(k) = t
                end if
            end do

        end if

    end subroutine sgesl

    real function sasum(n, sx, incx)

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
        implicit none

        integer, intent(IN) :: n
        integer, intent(IN) :: incx
        real, intent(IN)    :: sx(*)

        integer :: i, m

        sasum = 0.0
        if (n <= 0) return
        if (incx /= 1) then
!                                          ** NON-UNIT INCREMENTS
            do i = 1, 1 + (n - 1)*incx, incx
                sasum = sasum + abs(sx(i))
            end do
        else
!                                          ** UNIT INCREMENTS
            m = mod(n, 6)
            if (m /= 0) then
!                             ** CLEAN-UP LOOP SO REMAINING VECTOR
!                             ** LENGTH IS A MULTIPLE OF 6.
                do i = 1, m
                    sasum = sasum + abs(sx(i))
                end do
            end if
!                              ** UNROLL LOOP FOR SPEED
            do i = m + 1, n, 6
                sasum = sasum + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2)) &
                        + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
            end do
        end if

    end function sasum

    subroutine saxpy(n, sa, sx, incx, sy, incy)

!          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)

!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

! --OUTPUT--
!       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH
!                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
!            WHERE LX = 1          IF INCX .GE. 0,
!                     = (-INCX)*N  IF INCX .LT. 0
!            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
        implicit none

        integer, intent(IN)  :: n
        integer, intent(IN)  :: incx
        integer, intent(IN)  :: incy
        real, intent(IN)     :: sa
        real, intent(IN)     :: sx(*)
        real, intent(INOUT)    :: sy(*)

        integer :: i, m, ix, iy

        if (n <= 0 .or. sa == 0.0) return

        if (incx == incy .and. incx > 1) then

            do i = 1, 1 + (n - 1)*incx, incx
                sy(i) = sy(i) + sa*sx(i)
            end do

        else if (incx == incy .and. incx == 1) then

!                                        ** EQUAL, UNIT INCREMENTS
            m = mod(n, 4)
            if (m /= 0) then
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 4.
                do i = 1, m
                    sy(i) = sy(i) + sa*sx(i)
                end do
            end if
!                              ** UNROLL LOOP FOR SPEED
            do i = m + 1, n, 4
                sy(i) = sy(i) + sa*sx(i)
                sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
                sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
                sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
            end do

        else
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
            ix = 1
            iy = 1
            if (incx < 0) ix = 1 + (n - 1)*(-incx)
            if (incy < 0) iy = 1 + (n - 1)*(-incy)
            do i = 1, n
                sy(iy) = sy(iy) + sa*sx(ix)
                ix = ix + incx
                iy = iy + incy
            end do

        end if

    end subroutine saxpy

    real function sdot(n, sx, incx, sy, incy)

!          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'

!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

! --OUTPUT--
!     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
!            WHERE  LX = 1          IF INCX .GE. 0,
!                      = (-INCX)*N  IF INCX .LT. 0,
!            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
        implicit none

        integer, intent(IN) :: n
        integer, intent(IN) :: incx
        integer, intent(IN) :: incy
        real, intent(IN)    :: sx(*)
        real, intent(IN)    :: sy(*)

        integer :: i, m, ix, iy

        sdot = 0.0
        if (n <= 0) return

        if (incx == incy .and. incx > 1) then

            do i = 1, 1 + (n - 1)*incx, incx
                sdot = sdot + sx(i)*sy(i)
            end do

        else if (incx == incy .and. incx == 1) then

!                                        ** EQUAL, UNIT INCREMENTS
            m = mod(n, 5)
            if (m /= 0) then
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 4.
                do i = 1, m
                    sdot = sdot + sx(i)*sy(i)
                end do
            end if
!                              ** UNROLL LOOP FOR SPEED
            do i = m + 1, n, 5
                sdot = sdot + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) &
                       + sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
            end do

        else
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
            ix = 1
            iy = 1
            if (incx < 0) ix = 1 + (n - 1)*(-incx)
            if (incy < 0) iy = 1 + (n - 1)*(-incy)
            do i = 1, n
                sdot = sdot + sx(ix)*sy(iy)
                ix = ix + incx
                iy = iy + incy
            end do

        end if

    end function sdot

    subroutine sscal(n, sa, sx, incx)

!         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
!            SA  SINGLE PRECISION SCALE FACTOR
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX)
!                FOR I = 0 TO N-1
        implicit none

        integer, intent(IN)   :: n
        integer, intent(IN)   :: incx
        real, intent(IN)      :: sa
        real, intent(INOUT)   :: sx(*)

        integer :: i, m

        if (n <= 0) return

        if (incx /= 1) then

            do i = 1, 1 + (n - 1)*incx, incx
                sx(i) = sa*sx(i)
            end do

        else

            m = mod(n, 5)
            if (m /= 0) then
!                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                           ** IS A MULTIPLE OF 5.
                do i = 1, m
                    sx(i) = sa*sx(i)
                end do
            end if
!                             ** UNROLL LOOP FOR SPEED
            do i = m + 1, n, 5
                sx(i) = sa*sx(i)
                sx(i + 1) = sa*sx(i + 1)
                sx(i + 2) = sa*sx(i + 2)
                sx(i + 3) = sa*sx(i + 3)
                sx(i + 4) = sa*sx(i + 4)
            end do

        end if

    end subroutine sscal

    subroutine sswap(n, sx, incx, sy, incy)

!          INTERCHANGE S.P VECTORS  X  AND  Y

!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

! --OUTPUT--
!       SX  INPUT VECTOR SY (UNCHANGED IF N .LE. 0)
!       SY  INPUT VECTOR SX (UNCHANGED IF N .LE. 0)

!     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
!     WHERE LX = 1          IF INCX .GE. 0,
!              = (-INCX)*N  IF INCX .LT. 0
!     AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
        implicit none

        integer, intent(IN)  :: n
        integer, intent(IN)  :: incx
        integer, intent(IN)  :: incy
        real, intent(IN OUT) :: sx(*)
        real, intent(IN OUT) :: sy(*)

        real :: stemp1, stemp2, stemp3
        integer :: ix, iy, i, m

        if (n <= 0) return

        if (incx == incy .and. incx > 1) then

            do i = 1, 1 + (n - 1)*incx, incx
                stemp1 = sx(i)
                sx(i) = sy(i)
                sy(i) = stemp1
            end do

        else if (incx == incy .and. incx == 1) then

!                                        ** EQUAL, UNIT INCREMENTS
            m = mod(n, 3)
            if (m /= 0) then
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 3.
                do i = 1, m
                    stemp1 = sx(i)
                    sx(i) = sy(i)
                    sy(i) = stemp1
                end do
            end if
!                              ** UNROLL LOOP FOR SPEED
            do i = m + 1, n, 3
                stemp1 = sx(i)
                stemp2 = sx(i + 1)
                stemp3 = sx(i + 2)
                sx(i) = sy(i)
                sx(i + 1) = sy(i + 1)
                sx(i + 2) = sy(i + 2)
                sy(i) = stemp1
                sy(i + 1) = stemp2
                sy(i + 2) = stemp3
            end do

        else
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
            ix = 1
            iy = 1
            if (incx < 0) ix = 1 + (n - 1)*(-incx)
            if (incy < 0) iy = 1 + (n - 1)*(-incy)
            do i = 1, n
                stemp1 = sx(ix)
                sx(ix) = sy(iy)
                sy(iy) = stemp1
                ix = ix + incx
                iy = iy + incy
            end do

        end if

    end subroutine sswap

    integer function isamax(n, sx, incx)

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
!                         ABS(SX(1+(I-1)*INCX))
        implicit none

        integer, intent(IN)                      :: n
        integer, intent(IN)                      :: incx
        real, intent(IN OUT)                     :: sx(*)

        real :: smax, xmag
        integer :: iii, i

        if (n <= 0) then
            isamax = 0
        else if (n == 1) then
            isamax = 1
        else
            smax = 0.0
            iii = 1
            do i = 1, 1 + (n - 1)*incx, incx
                xmag = abs(sx(i))
                if (smax < xmag) then
                    smax = xmag
                    isamax = iii
                end if
                iii = iii + 1
            end do
        end if

    end function isamax

    double precision function d1mach(i)

!-----------------------------------------------------------------------------*
!= PURPOSE:                                                                  =*
!= D1MACH calculates various machine constants in single precision.          =*
!-----------------------------------------------------------------------------*
!= PARAMETERS:                                                               =*
!=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
!=   D1MACH  -  REAL, machine constant in single precision               (O) =*
!=      I=1     - the smallest non-vanishing normalized floating-point       =*
!=                power of the radix, i.e., D1MACH=FLOAT(IBETA)**MINEXP      =*
!=      I=2     - the largest finite floating-point number.  In              =*
!=                particular D1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
!=                Note - on some machines D1MACH will be only the            =*
!=                second, or perhaps third, largest number, being            =*
!=                too small by 1 or 2 units in the last digit of             =*
!=                the significand.                                           =*
!=      I=3     - A small positive floating-point number such that           =*
!=                1.0-D1MACH .NE. 1.0. In particular, if IBETA = 2           =*
!=                or  IRND = 0, D1MACH = FLOAT(IBETA)**NEGEPS.               =*
!=                Otherwise,  D1MACH = (IBETA**NEGEPS)/2.  Because           =*
!=                NEGEPS is bounded below by -(IT+3), D1MACH may not         =*
!=                be the smallest number that can alter 1.0 by               =*
!=                subtraction.                                               =*
!=      I=4     - the smallest positive floating-point number such           =*
!=                that  1.0+D1MACH .NE. 1.0. In particular, if either        =*
!=                IBETA = 2  or  IRND = 0, D1MACH=FLOAT(IBETA)**MACHEP.      =*
!=                Otherwise, D1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
!=  (see routine T665D for more information on different constants)          =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN)                      :: i

        select case (i)
        case (1)
            d1mach = tiny(1.0d0)
        case (2)
            d1mach = huge(1.0d0)
        case (3)
            d1mach = -epsilon(1.0d0)
        case (4)
            d1mach = epsilon(1.0d0)
        case default
            write (0, *) '>>> ERROR (D1MACH) <<<  invalid argument'
            stop
        end select

    end function d1mach

    real function r1mach(i)

!-----------------------------------------------------------------------------*
!= PURPOSE:                                                                  =*
!= R1MACH calculates various machine constants in single precision.          =*
!-----------------------------------------------------------------------------*
!= PARAMETERS:                                                               =*
!=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
!=   R1MACH  -  REAL, machine constant in single precision               (O) =*
!=      I=1     - the smallest non-vanishing normalized floating-point       =*
!=                power of the radix, i.e., R1MACH=FLOAT(IBETA)**MINEXP      =*
!=      I=2     - the largest finite floating-point number.  In              =*
!=                particular R1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
!=                Note - on some machines R1MACH will be only the            =*
!=                second, or perhaps third, largest number, being            =*
!=                too small by 1 or 2 units in the last digit of             =*
!=                the significand.                                           =*
!=      I=3     - A small positive floating-point number such that           =*
!=                1.0-R1MACH .NE. 1.0. In particular, if IBETA = 2           =*
!=                or  IRND = 0, R1MACH = FLOAT(IBETA)**NEGEPS.               =*
!=                Otherwise,  R1MACH = (IBETA**NEGEPS)/2.  Because           =*
!=                NEGEPS is bounded below by -(IT+3), R1MACH may not         =*
!=                be the smallest number that can alter 1.0 by               =*
!=                subtraction.                                               =*
!=      I=4     - the smallest positive floating-point number such           =*
!=                that  1.0+R1MACH .NE. 1.0. In particular, if either        =*
!=                IBETA = 2  or  IRND = 0, R1MACH=FLOAT(IBETA)**MACHEP.      =*
!=                Otherwise, R1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
!=  (see routine T665R for more information on different constants)          =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN)                      :: i

        select case (i)
        case (1)
            r1mach = tiny(1.0)
        case (2)
            r1mach = huge(1.0)
        case (3)
            r1mach = -epsilon(1.0)
        case (4)
            r1mach = epsilon(1.0)
        case default
            write (0, *) '>>> ERROR (D1MACH) <<<  invalid argument'
            stop
        end select

    end function r1mach

    subroutine swchem(nw, wl, nz, nj, tLev, airDen, sj)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Load various "weighting functions" (products of cross section and        =*
!=  quantum yield at each altitude and each wavelength).  The altitude       =*
!=  dependence is necessary to ensure the consideration of pressure and      =*
!=  temperature dependence of the cross sections or quantum yields.          =*
!=  The actual reading, evaluation and interpolation is done in separate     =*
!=  subroutines for ease of management and manipulation.  Please refer to    =*
!=  the inline documentation of the specific subroutines for detail          =*
!=  information.                                                             =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section * quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!=  OPT    - If opt=1 read, otherwise just calc
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN)   :: nw
        integer, intent(IN)   :: nz
        real, intent(IN)      :: wl(kw)
        real, intent(IN)      :: tLev(nz)
        real, intent(IN)      :: airDen(nz)
        real, intent(INOUT)   :: sj(kj, nz, kw)
        integer, intent(OUT)  :: nj

        real :: wc(kw)
        integer :: iw

        integer :: j

! complete wavelength grid
        do iw = 1, nw - 1
            wc(iw) = (wl(iw) + wl(iw + 1))/2.
        end do

!____________________________________________________________________________

! O2 + hv -> O + O
! reserve first position.  Cross section parameterization in Schumman-Runge and
! Lyman-alpha regions are zenith-angle dependent, will be written in
! subroutine seto2.f.

        j = 1
!jlabel(j) = 'O2 -> O + O'

! O3 + hv ->  (both channels)
        if (doReaction(2) .or. doReaction(3)) call r01(nw, wl, wc, nz, j, sj, tLev) !1,2

! NO2 + hv -> NO + O(3P)
        if (doReaction(4)) call r02(nw, nz, j, sj, tLev) !3

! NO3 + hv ->  (both channels)
        if (doReaction(5) .or. doReaction(6)) call r03(nw, wc, nz, j, sj, tLev) !4,5

! N2O5 + hv -> (both channels)
        if (doReaction(7) .or. doReaction(8)) call r04(nw, wc, nz, j, sj, tLev) !6,7

! N2O + hv -> N2 + O(1D)
        if (doReaction(9)) call r44(nw, wc, nz, j, sj, tLev) !8

! HO2 + hv -> OH + O
        if (doReaction(10)) call r39(nw, wc, nz, j, sj, tLev) !9

! H2O2 + hv -> 2 OH
        if (doReaction(11)) call r08(nw, wl, wc, nz, j, sj, tLev) !10

! HNO2 + hv -> OH + NO
        if (doReaction(12)) call r05(nw, nz, j, sj, tLev) !11

! HNO3 + hv -> OH + NO2
        if (doReaction(13)) call r06(nw, nz, j, sj, tLev) !12

! HNO4 + hv -> HO2 + NO2
        if (doReaction(14)) call r07(nw, nz, j, sj, tLev) !13

! CH2O + hv -> (both channels)
        if (doReaction(15) .or. doReaction(16)) &
            call r10(nw, wl, wc, nz, j, sj, tLev, airDen) !14,15

! CH3CHO + hv -> (all three channels)
        if (doReaction(17) .or. doReaction(18) .or. doReaction(19)) &
            call r11(nw, nz, j, sj, tLev, airDen) !16,17,18

! C2H5CHO + hv -> C2H5 + HCO
        if (doReaction(20)) call r12(nw, nz, j, sj, tLev, airDen) !19

! CHOCHO + hv -> Products
        if (doReaction(21) .or. doReaction(22)) call r13(nw, wc, nz, j, sj, tLev) !20,21

! CH3COCHO + hv -> Products
        if (doReaction(23)) call r14(nw, wc, nz, j, sj, tLev, airDen) !22

! CH3COCH3 + hv -> Products
        if (doReaction(24)) call r15(nw, wc, nz, j, sj, tLev, airDen) !23

! CH3OOH + hv -> CH3O + OH
        if (doReaction(25)) call r16(nw, nz, j, sj, tLev) !24

! CH3ONO2 + hv -> CH3O + NO2
        if (doReaction(26)) call r17(nw, nz, j, sj, tLev) !25

! PAN + hv -> Products
        if (doReaction(27)) call r18(nw, nz, j, sj, tLev) !26

! ClOO + hv -> Products
        if (doReaction(28)) call r31(nw, nz, j, sj, tLev) !27

! ClONO2 + hv -> Products
        if (doReaction(29) .or. doReaction(30)) call r45(nw, wc, nz, j, sj, tLev) !28,29

! CH3Cl + hv -> Products
        if (doReaction(31)) call r30(nw, nz, j, sj, tLev) !30

! CCl2O + hv -> Products
        if (doReaction(32)) call r19(nw, nz, j, sj, tLev) !31

! CCl4 + hv -> Products
        if (doReaction(33)) call r20(nw, nz, j, sj, tLev) !32

! CClFO + hv -> Products
        if (doReaction(34)) call r21(nw, nz, j, sj, tLev) !33

! CCF2O + hv -> Products
        if (doReaction(35)) call r22(nw, nz, j, sj, tLev) !34

! CF2ClCFCl2 (CFC-113) + hv -> Products
        if (doReaction(36)) call r23(nw, nz, j, sj, tLev) !35

! CF2ClCF2Cl (CFC-114) + hv -> Products
        if (doReaction(37)) call r24(nw, nz, j, sj, tLev) !36

! CF3CF2Cl (CFC-115) + hv -> Products
        if (doReaction(38)) call r25(nw, nz, j, sj, tLev) !37

! CCl3F (CFC-111) + hv -> Products
        if (doReaction(39)) call r26(nw, wc, nz, j, sj, tLev) !38

! CCl2F2 (CFC-112) + hv -> Products
        if (doReaction(40)) call r27(nw, nz, j, sj, tLev) !39

! CH3CCl3 + hv -> Products
        if (doReaction(41)) call r29(nw, nz, j, sj, tLev) !40

! CF3CHCl2 (HCFC-123) + hv -> Products
        if (doReaction(42)) call r32(nw, wc, nz, j, sj, tLev) !41

! CF3CHFCl (HCFC-124) + hv -> Products
        if (doReaction(43)) call r33(nw, wc, nz, j, sj, tLev) !42

! CH3CFCl2 (HCFC-141b) + hv -> Products
        if (doReaction(44)) call r34(nw, nz, j, sj, tLev) !43

! CH3CF2Cl (HCFC-142b) + hv -> Products
        if (doReaction(45)) call r35(nw, wc, nz, j, sj, tLev) !44

! CF3CF2CHCl2 (HCFC-225ca) + hv -> Products
        if (doReaction(46)) call r36(nw, nz, j, sj, tLev) !45

! CF2ClCF2CHFCl (HCFC-225cb) + hv -> Products
        if (doReaction(47)) call r37(nw, nz, j, sj, tLev) !46

! CHClF2 (HCFC-22) + hv -> Products
        if (doReaction(48)) call r38(nw, nz, j, sj, tLev) !47

! BrONO2 + hv -> Products
        if (doReaction(49) .or. doReaction(50)) call r46(nw, nz, j, sj, tLev) !48,49

! CH3Br + hv -> Products
        if (doReaction(51)) call r28(nw, nz, j, sj, tLev) !50

! CHBr3 + hv -> Products
        if (doReaction(52)) call r09(nw, wl, wc, nz, j, sj, tLev) !51

! CF3Br (Halon-1301) + hv -> Products
        if (doReaction(53)) call r42(nw, nz, j, sj, tLev) !52

! CF2BrCF2Br (Halon-2402) + hv -> Products
        if (doReaction(54)) call r43(nw, nz, j, sj, tLev) !53

! CF2Br2 (Halon-1202) + hv -> Products
        if (doReaction(55)) call r40(nw, nz, j, sj, tLev) !54

! CF2BrCl (Halon-1211) + hv -> Products
        if (doReaction(56)) call r41(nw, nz, j, sj, tLev) !55

! CL2 + hc -> CL + CL
        if (doReaction(57)) call r47(nw, nz, j, sj, tLev) !56

! CH2(OH)CH + hv -> Products
        if (doReaction(58)) call r101(nw, nz, j, sj, tLev) !57

! CH3COCOCH3 + hv -> Products
        if (doReaction(59)) call r102(nw, nz, j, sj, tLev) !58

! CH3COCHCH2 + hv -> Products
        if (doReaction(60)) call r103(nw, wc, nz, j, sj, tLev, airDen) !59

! CH2C(CH3)CHO + hv -> Products
        if (doReaction(61)) call r104(nw, nz, j, sj, tLev) !60

! CH3COCO(OH) + hv -> Products
        if (doReaction(62)) call r105(nw, nz, j, sj, tLev) !61

! CH3CH2ONO2 -> CH3CH2O + NO2
        if (doReaction(63)) call r106(nw, nz, j, sj, tLev) !62

! CH3CHONO2CH3 -> CH3CHOCH3 + NO2
        if (doReaction(64)) call r107(nw, nz, j, sj, tLev) !63

! CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2
        if (doReaction(65)) call r108(nw, wc, nz, j, sj, tLev) !64

! CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
        if (doReaction(66)) call r109(nw, wc, nz, j, sj, tLev) !65

! C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
        if (doReaction(67)) call r110(nw, wc, nz, j, sj, tLev) !66

! ClOOCl -> Cl + ClOO
        if (doReaction(68)) call r111(nw, nz, j, sj, tLev) !67

! CH2(OH)COCH3 -> CH3CO + CH2(OH)
! CH2(OH)COCH3 -> CH2(OH)CO + CH3
        if (doReaction(69) .or. doReaction(70)) call r112(nw, nz, j, sj, tLev) !68,69

! HOBr -> OH + Br'
        if (doReaction(71)) call r113(nw, wc, nz, j, sj, tLev) !70

! BrO -> Br + O'
        if (doReaction(72)) call r114(nw, nz, j, sj) !71

! Br2 -> Br + Br'
        if (doReaction(73)) call r115(nw, nz, j, sj, tLev) !72

!***************************************************************

        if (j > kj) stop '1002'

        nj = j

    end subroutine swchem

    subroutine r01(nw, wl, wc, nz, j, sj, tLev)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product of (cross section) x (quantum yield) for the two     =*
!=  O3 photolysis reactions:                                                 =*
!=             (a) O3 + hv -> O2 + O(1D)                                     =*
!=             (b) O3 + hv -> O2 + O(3P)                                     =*
!=  Cross section:  Combined data from WMO 85 Ozone Assessment (use 273K     =*
!=                  value from 175.439-847.5 nm) and data from Molina and    =*
!=                  Molina (use in Hartley and Huggins bans (240.5-350 nm)   =*
!=  Quantum yield:  Choice between                                           =*
!=                   (1) data from Michelsen et al, 1994                     =*
!=                   (2) JPL 87 recommendation                               =*
!=                   (3) JPL 90/92 recommendation (no "tail")                =*
!=                   (4) data from Shetter et al., 1996                      =*
!=                   (5) JPL 97 recommendation                               =*
!=                   (6) JPL 00 recommendation                               =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!=  OPT    - If opt=1 read, otherwise just calc
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 01 !1,2

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        real, intent(IN)                :: wl(kw)
        real, intent(IN)                :: wc(kw)
        real, intent(IN)                :: tLev(nz)
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 500
        real :: xs(nz, kw)
        real :: qy1d
        real :: qy3p
        real :: tau, tau2, tau3
        real :: a, b, c
        real :: a0, a1, a2, a3, a4, a5, a6 !, a7
        real :: xl, xl0
        integer :: i, iw

        j = j + 1
!jlabel(j) = 'O3 -> O2 + O(1D)'
        j = j + 1
!jlabel(j) = 'O3 -> O2 + O(3P)'

        select case (mOption(1))
        case (1)
            call o3xs_mm(nw, wl, nz, xs, tLev)
        case (2)
            call o3xs_mal(nw, wl, nz, xs, tLev)
        case (3)
            call o3xs_bass(nw, wl, nz, xs, tLev)
        end select

! compute cross sections and yields at different wavelengths, altitudes:

! quantum yields
        do iw = 1, nw - 1
            do i = 1, nz
! coefficients from jpl 87:
                if (mOption(2) == kjpl87) then
                    tau = tLev(i) - 230.
                    tau2 = tau*tau
                    tau3 = tau2*tau
                    xl = wc(iw)
                    xl0 = 308.2 + 4.4871e-2*tau + 6.938e-5*tau2 - 2.5452e-6*tau3
                    a = 0.9*(0.369 + 2.85e-4*tau + 1.28e-5*tau2 + 2.57e-8*tau3)
                    b = -0.575 + 5.59e-3*tau - 1.439e-5*tau2 - 3.27e-8*tau3
                    c = 0.9*(0.518 + 9.87e-4*tau - 3.94e-5*tau2 + 3.91e-7*tau3)
                    qy1d = a*atan(b*(xl - xl0)) + c
                    qy1d = AMAX1(0., qy1d)
                    qy1d = AMIN1(0.9, qy1d)
                end if

! from jpl90, jpl92:
! (caution: error in JPL92 for first term of a3)
                if (mOption(2) == kjpl92) then
                    tau = 298.-tLev(i)
                    tau2 = tau*tau
                    xl0 = wc(iw) - 305.
                    a0 = .94932 - 1.7039e-4*tau + 1.4072e-6*tau2
                    a1 = -2.4052e-2 + 1.0479e-3*tau - 1.0655e-5*tau2
                    a2 = 1.8771e-2 - 3.6401e-4*tau - 1.8587e-5*tau2
                    a3 = -1.4540e-2 - 4.7787e-5*tau + 8.1277e-6*tau2
                    a4 = 2.3287e-3 + 1.9891e-5*tau - 1.1801e-6*tau2
                    a5 = -1.4471e-4 - 1.7188e-6*tau + 7.2661e-8*tau2
                    a6 = 3.1830e-6 + 4.6209e-8*tau - 1.6266e-9*tau2
                    qy1d = a0 + a1*xl0 + a2*(xl0)**2 + a3*(xl0)**3 + &
                            a4*(xl0)**4 + a5*(xl0)**5 + a6*(xl0)**6
                    if (wc(iw) < 305.) qy1d = 0.95
                    if (wc(iw) > 320.) qy1d = 0.
                    if (qy1d < 0.02) qy1d = 0.
                end if

! from JPL'97
                if (mOption(2) == kjpl97) then
                    if (wc(iw) < 271.) then
                        qy1d = 0.87
                    else if (wc(iw) >= 271. .and. wc(iw) < 290.) then
                        qy1d = 0.87 + (wc(iw) - 271.)*(.95 - .87)/(290.-271.)
                    else if (wc(iw) >= 290. .and. wc(iw) < 305.) then
                        qy1d = 0.95
                    else if (wc(iw) >= 305. .and. wc(iw) <= 325.) then
                        if (i > nz) exit
                        qy1d = yg1(iw, nReact)*exp(-yg2(iw, nReact)/ &
                                                        tLev(i))
                    else
                        qy1d = 0.
                    end if
                end if

! from Michelsen, H. A., R.J. Salawitch, P. O. Wennber, and J. G. Anderson
! Geophys. Res. Lett., 21, 2227-2230, 1994.
                if (mOption(2) == kmich) then
                    if (wc(iw) < 271.) then
                        qy1d = 0.87
                    else if (wc(iw) >= 271. .and. wc(iw) < 305.) then
                        qy1d = 1.98 - 301./wc(iw)
                    else if (wc(iw) >= 305. .and. wc(iw) <= 325.) then
                        if (i > nz) exit
                        qy1d = yg1(iw, nReact)*exp(-yg2(iw, nReact)/ &
                                                        (0.6951*tLev(i)))
                    else
                        qy1d = 0.
                    end if
                end if

! Shetter et al.:
! phi = A * exp(-B/T), A and B are based on meas. at 298 and 230 K
! do linear interpolation between phi(298) and phi(230) for wavelengths > 321
! as phi(230)=0. for those wavelengths, so there are no A and B factors
                if (mOption(2) == kshet) then
                    if (wl(iw + 1) <= 321.) then
                        qy1d = yg1(iw, nReact)*exp(-1.*yg2(iw, nReact)/ &
                                                        tLev(i))
                    else
                        qy1d = (yg3(iw, nReact) - yg4(iw, nReact))/(298.-230.)* &
                                (tLev(i) - 230.) + yg4(iw, nReact)
                    end if
                end if

! JPL 2000:
                if (mOption(2) == kjpl00) then
                    qy1d = fo3qy(wc(iw), tLev(i))
                end if

! Matsumi et al.
                if (mOption(2) == kmats) then
                    qy1d = fo3qy2(wc(iw), tLev(i))
                end if

! compute product
                sj(2, i, iw) = qy1d*xs(i, iw) !1
                qy3p = 1.0 - qy1d
                sj(3, i, iw) = qy3p*xs(i, iw) !2
            end do
        end do

    end subroutine r01

!=============================================================================*
    subroutine r02(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for NO2            =*
!=  photolysis:                                                              =*
!=         NO2 + hv -> NO + O(3P)                                            =*
!=  Cross section from JPL94 (can also have Davidson et al.)                 =*
!=  Quantum yield from Gardiner, Sperry, and Calvert                         =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!=  OPT    - If opt=1 read, otherwise just calc
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 02 !3

        integer, intent(IN)             :: nw
        integer, intent(in)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)

        integer, parameter :: kdata = 200
        integer :: i, iw
        real :: r2_no2xs(nz, kw)

!*************** NO2 photodissociation
        j = j + 1
!jlabel(j) = 'NO2 -> NO + O(3P)'

        select case (mOption(3))
        case (1)
            do iw = 1, nw - 1
                do i = 1, nz
                    r2_no2xs(i, iw) = yg1(iw, nReact) + &
                                       yg2(iw, nReact)*(tLev(i) - 273.15)
                end do
            end do
        case (2)
            do iw = 1, nw - 1
                do i = 1, nz
                    r2_no2xs(i, iw) = yg1(iw, nReact)
                end do
            end do
        end select

! combine
        do iw = 1, nw - 1
            do i = 1, nz
                sj(4, i, iw) = r2_no2xs(i, iw)* &
                                yg1n(iw, nReact)         !3
            end do
        end do

    end subroutine r02

!=============================================================================*

    subroutine r03(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (absorptioon cross section) x (quantum yield) for    =*
!=  both channels of NO3 photolysis:                                         =*
!=          (a) NO3 + hv -> NO2 + O(3P)                                      =*
!=          (b) NO3 + hv -> NO + O2                                          =*
!=  Cross section combined from Graham and Johnston (<600 nm) and JPL 94     =*
!=  Quantum yield from Madronich (1988)                                      =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!=  OPT    - If opt=1 read, otherwise just calc
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 03 !4,5

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)

        integer, parameter :: kdata = 350
        real :: qyVal
        integer :: i, iw

!***************      jlabel(j) = 'NO3 -> NO2 + O(3P)'
!***************      jlabel(j) = 'NO3 -> NO + O2'
! quantum yield:
! from Madronich (1988) see CEC NO3 book.
        j = j + 1
!jlabel(j) = 'NO3 -> NO + O2'

! for   NO3 ->NO+O2
        do iw = 1, nw - 1
            if (wc(iw) < 584.) then
                qyVal = 0.
            else if (wc(iw) >= 640.) then
                qyVal = 0.
            else if (wc(iw) >= 595.) then
                qyVal = 0.35*(1.-(wc(iw) - 595.)/45.)
            else
                qyVal = 0.35*(wc(iw) - 584.)/11.
            end if
            do i = 1, nz
                sj(5, i, iw) = yg(iw, nReact)*qyVal !4
            end do
        end do

! for  NO3 ->NO2+O
        j = j + 1

!jlabel(j) = 'NO3 -> NO2 + O(3P)'
! for  NO3 ->NO2+O
        do iw = 1, nw - 1
            if (wc(iw) < 584.) then
                qyVal = 1.
            else if (wc(iw) > 640.) then
                qyVal = 0.
            else if (wc(iw) > 595.) then
                qyVal = 0.65*(1 - (wc(iw) - 595.)/45.)
            else
                qyVal = 1.-0.35*(wc(iw) - 584.)/11.
            end if
            do i = 1, nz
                sj(6, i, iw) = yg(iw, nReact)*qyVal !5
            end do
        end do

    end subroutine r03

!=============================================================================*

    subroutine r04(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yiels) for N2O5 photolysis =*
!=  reactions:                                                               =*
!=       (a) N2O5 + hv -> NO3 + NO + O(3P)                                   =*
!=       (b) N2O5 + hv -> NO3 + NO2                                          =*
!=  Cross section from JPL97: use tabulated values up to 280 nm, use expon.  =*
!=                            expression for >285nm, linearly interpolate    =*
!=                            between s(280) and s(285,T) in between         =*
!=  Quantum yield: Analysis of data in JPL94 (->dataj1/YLD/N2O5.qy)          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!=  OPT    - If opt=1 read, otherwise just calc
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 04 !6,7

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)

        integer, parameter :: kdata = 100
        real :: qy
        real :: xs, xs270, xs280, xst290, xst300
        real :: dum1, dum2
        real :: t
        integer :: i, iw

!*************** N2O5 photodissociation
        j = j + 1
!jlabel(j) = 'N2O5 -> NO3 + NO + O(3P)'
        j = j + 1
!jlabel(j) = 'N2O5 -> NO3 + NO2'
! quantum yield : see dataj1/YLD/N2O5.qy for explanation
! correct for T-dependence of cross section

        do iw = 1, nw - 1
            qy = min(1., 3.832441 - 0.012809638*wc(iw))
            qy = max(0., qy)
            do i = 1, nz
! temperature dependence only valid for 225 - 300 K.
                t = max(225., min(tLev(i), 300.))
! evaluation of exponential
                if (wc(iw) >= 285. .and. wc(iw) <= 380.) then
                    sj(7, i, iw) = qy*1.e-20* &
                                    exp(2.735 + (4728.5 - 17.127*wc(iw))/t) !6
                    sj(8, i, iw) = (1.-qy)*1.e-20* &
                                    exp(2.735 + (4728.5 - 17.127*wc(iw))/t) !7
! between 280 and 285 nm:  Extrapolate from both sides, then average.
                else if (wc(iw) >= 280. .and. wc(iw) < 285.) then
                    xst290 = 1.e-20*exp(2.735 + (4728.5 - 17.127*290.)/t)
                    xst300 = 1.e-20*exp(2.735 + (4728.5 - 17.127*300.)/t)
                    dum1 = xs270 + (wc(iw) - 270.)*(xs280 - xs270)/10.
                    dum2 = xst290 + (wc(iw) - 290.)*(xst300 - xst290)/10.
                    xs = 0.5*(dum1 + dum2)
                    sj(7, i, iw) = qy*xs !6
                    sj(8, i, iw) = (1.-qy)*xs !7
! for less than 280 nm, use tabulated values
                else if (wc(iw) < 280.) then
                    sj(7, i, iw) = qy*yg(iw, nReact) !6
                    sj(8, i, iw) = (1.-qy)*yg(iw, nReact) !7
! beyond 380 nm, set to zero
                else
                    sj(7, i, iw) = 0. !6
                    sj(8, i, iw) = 0. !7
                end if
            end do
        end do
    end subroutine r04

!=============================================================================*

    subroutine r05(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for HNO2 photolysis=*
!=     HNO2 + hv -> NO + OH                                                  =*
!=  Cross section:  from JPL97                                               =*
!=  Quantum yield:  assumed to be unity                                      =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 05 !11

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: i, iw

!*************** HNO2 photodissociation
! cross section from JPL92
! (from Bongartz et al., identical to JPL94, JPL97 recommendation)
        j = j + 1
!jlabel(j) = 'HNO2 -> OH + NO'
! quantum yield = 1
        qy = 1.
        do iw = 1, nw - 1
            do i = 1, nz
                sj(12, i, iw) = yg(iw, nReact)*qy !11
            end do
        end do

    end subroutine r05

!=============================================================================*

    subroutine r06(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yield) for HNO3 photolysis =*
!=        HNO3 + hv -> OH + NO2                                              =*
!=  Cross section: Burkholder et al., 1993                                   =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 06 !12

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)

        integer, parameter :: kdata = 100
        integer :: i, iw

!*************** HNO3 photodissociation
        j = j + 1
!jlabel(j) = 'HNO3 -> OH + NO2'
! quantum yield = 1
! correct for temperature dependence

        do iw = 1, nw - 1
            do i = 1, nz
                sj(13, i, iw) = yg1(iw, nReact)* &
                                 1.e-20*exp(yg2(iw, nReact)/1.e3* &
                                            (tLev(i) - 298.)) !12
            end do
        end do

    end subroutine r06

!=============================================================================*

    subroutine r07(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yield) for HNO4 photolysis =*
!=       HNO4 + hv -> HO2 + NO2                                              =*
!=  Cross section:  from JPL97                                               =*
!=  Quantum yield:  Assumed to be unity                                      =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 07 !13

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: i, iw

!*************** HNO4 photodissociation
! cross section from JPL85 (identical to JPL92 and JPL94 and JPL97)

        j = j + 1
!jlabel(j) = 'HNO4 -> HO2 + NO2'

! quantum yield = 1

        qy = 1.
        do iw = 1, nw - 1
            do i = 1, nz
                sj(14, i, iw) = yg(iw, nReact)*qy !13
            end do
        end do

    end subroutine r07

!=============================================================================*

    subroutine r08(nw, wl, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yield) for H2O2 photolysis =*
!=         H2O2 + hv -> 2 OH                                                 =*
!=  Cross section:  From JPL97, tabulated values @ 298K for <260nm, T-depend.=*
!=                  parameterization for 260-350nm                           =*
!=  Quantum yield:  Assumed to be unity                                      =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 08 !10

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wl(kw)
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)

        integer, parameter :: kdata = 100
        real :: qy
        real :: a0, a1, a2, a3, a4, a5, a6, a7
        real :: b0, b1, b2, b3, b4
        real :: xs
        real :: t
        integer :: i, iw
        real :: lambda
        real :: suma, sumb, chi

!*************** H2O2 photodissociation

! cross section from Lin et al. 1978

        j = j + 1
!jlabel(j) = 'H2O2 -> 2 OH'

        a0 = 6.4761e+04
        a1 = -9.2170972e+02
        a2 = 4.535649
        a3 = -4.4589016e-03
        a4 = -4.035101e-05
        a5 = 1.6878206e-07
        a6 = -2.652014e-10
        a7 = 1.5534675e-13

        b0 = 6.8123e+03
        b1 = -5.1351e+01
        b2 = 1.1522e-01
        b3 = -3.0493e-05
        b4 = -1.0924e-07

! quantum yield = 1

        qy = 1.

        do iw = 1, nw - 1

! Parameterization (JPL94)
! Range 260-350 nm; 200-400 K

            if ((wl(iw) >= 260.) .and. (wl(iw) < 350.)) then

                lambda = wc(iw)
                suma = ((((((a7*lambda + a6)*lambda + a5)*lambda + &
                           a4)*lambda + a3)*lambda + a2)*lambda + a1)*lambda + a0
                sumb = (((b4*lambda + b3)*lambda + b2)*lambda + b1)*lambda + b0

                do i = 1, nz
                    t = min(max(tLev(i), 200.), 400.)
                    chi = 1./(1.+exp(-1265./t))
                    xs = (chi*suma + (1.-chi)*sumb)*1e-21
                    sj(11, i, iw) = xs*qy !10
                end do
            else
                do i = 1, nz
                    sj(11, i, iw) = yg(iw, nReact)*qy !10
                end do
            end if

        end do

    end subroutine r08

!=============================================================================*

    subroutine r09(nw, wl, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yield) for CHBr3 photolysis=*
!=          CHBr3 + hv -> Products                                           =*
!=  Cross section: Choice of data from Atlas (?Talukdar???) or JPL97         =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 09 !51

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wl(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)
        real, intent(OUT)               :: wc(kw)

        integer, parameter :: kdata = 200
        real :: t
        real :: qy
        integer :: iw
        integer :: iz

        do iw = 1, nw - 1
            wc(iw) = (wl(iw) + wl(iw + 1))/2.
        end do

!*************** CHBr3 photodissociation
        j = j + 1
!jlabel(j) = 'CHBr3 -> Products'

! option:
        select case (mOption(4))
        case (1)
! quantum yield = 1
            qy = 1.
            do iw = 1, nw - 1
                do iz = 1, nz
                    t = tLev(iz)
                    if (t >= 296.) then
                        yg(iw, nReact) = yg1(iw, nReact)
                    else if (t >= 286.) then
                        yg(iw, nReact) = yg1(iw, nReact) + (t - 286.)* &
                                         (yg2(iw, nReact) - yg1(iw, nReact))/10.
                    else if (t >= 276.) then
                        yg(iw, nReact) = yg2(iw, nReact) + (t - 276.)* &
                                         (yg3(iw, nReact) - yg2(iw, nReact))/10.
                    else if (t >= 266.) then
                        yg(iw, nReact) = yg3(iw, nReact) + (t - 266.)* &
                                         (yg4(iw, nReact) - yg3(iw, nReact))/10.
                    else if (t >= 256.) then
                        yg(iw, nReact) = yg4(iw, nReact) + (t - 256.)* &
                                         (yg5(iw, nReact) - yg4(iw, nReact))/10.
                    else if (t < 256.) then
                        yg(iw, nReact) = yg5(iw, nReact)
                    end if
                    sj(52, iz, iw) = yg(iw, nReact)*qy !51
                end do
            end do
! jpl97, with temperature dependence formula,
!w = 290 nm to 340 nm,
!T = 210K to 300 K
!sigma, cm2 = exp((0.06183-0.000241*w)*(273.-T)-(2.376+0.14757*w))
        case (2)
! quantum yield = 1
            qy = 1.
            do iw = 1, nw - 1
                do iz = 1, nz
                    t = tLev(iz)
                    yg(iw, nReact) = yg1(iw, nReact)
                    if (wc(iw) > 290. .and. wc(iw) < 340. .and. &
                        t > 210 .and. t < 300) then
                        yg(iw, nReact) = exp((0.06183 - 0.000241*wc(iw))*(273.-t) - &
                                             (2.376 + 0.14757*wc(iw)))
                    end if
                    sj(52, iz, iw) = yg(iw, nReact)*qy !51
                end do
            end do
        end select

    end subroutine r09

!=============================================================================*

    subroutine r10(nw, wl, wc, nz, j, sj, tLev, AirDen)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
!=        (a) CH2O + hv -> H + HCO                                           =*
!=        (b) CH2O + hv -> H2 + CO                                           =*
!=  Cross section: Choice between                                            =*
!=                 1) Bass et al., 1980 (resolution: 0.025 nm)               =*
!=                 2) Moortgat and Schneider (resolution: 1 nm)              =*
!=                 3) Cantrell et al. (orig res.) for > 301 nm,              =*
!=                    IUPAC 92, 97 elsewhere                                 =*
!=                 4) Cantrell et al. (2.5 nm res.) for > 301 nm,            =*
!=                    IUPAC 92, 97 elsewhere                                 =*
!=                 5) Rogers et al., 1990                                    =*
!=                 6) new NCAR recommendation, based on averages of          =*
!=                    Cantrell et al., Moortgat and Schneider, and Rogers    =*
!=                    et al.                                                 =*
!=  Quantum yield: Choice between                                            =*
!=                 1) Evaluation by Madronich 1991 (unpublished)             =*
!=                 2) IUPAC 89, 92, 97                                       =*
!=                 3) Madronich, based on 1), updated 1998.                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 10

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wl(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(OUT)               :: wc(kw)
        real, intent(IN)                :: tLev(nz)
        real, intent(IN)                :: airDen(nz)

        integer, parameter :: kdata = 16000
        integer :: iw
        real :: x(kdata)
        real :: phi1, phi2, phi20, ak300, akt
        real :: qy1, qy2
        real :: sig, slope
        real :: t
        integer :: i

        do iw = 1, nw - 1
            wc(iw) = (wl(iw) + wl(iw + 1))/2.
        end do

!***************************************************************
!*************** CH2O photodissociatation

        j = j + 1
!jlabel(j) = 'CH2O -> H + HCO'

        j = j + 1
!jlabel(j) = 'CH2O -> H2 + CO'

! working grid arrays:
!     yg1(:,nReact) = cross section at a specific temperature
!     yg2(:,nReact), yg3(:,nReact) = cross sections at different temp or slope, for calculating
!                temperature depedence
!     yg4(:,nReact) = quantum yield data for radical channel
!     yg5(:,nReact) = quantum yield data for molecular channel

! combine
! y1 = xsect
! y2 = xsect(223), Cantrell et al.
! y3 = xsect(293), Cantrell et al.
! y4 = qy for radical channel
! y5 = qy for molecular channel
! pressure and temperature dependent for w > 330.
        do iw = 1, nw - 1
            if (mOption(5) == 6) then
                sig = yg2(iw, nReact)
            else
                sig = yg1(iw, nReact)
            end if
            do i = 1, nz
! correct cross section for temperature dependence for > 301. nm
                if (wl(iw) >= 301.) then
                    t = max(223.15, min(tLev(i), 293.15))
                    if (mOption(5) == 3 .or. mOption(5) == 6) then
                        sig = yg2(iw, nReact) + yg3(iw, nReact)*(t - 273.15)
                    else if (mOption(5) == 4) then
                        slope = (yg3(iw, nReact) - yg2(iw, nReact))/(293.-223.)
                        sig = yg2(iw, nReact) + slope*(t - 223.)
                    end if
                end if
                sig = max(sig, 0.)
! quantum yields:
! temperature and pressure dependence beyond 330 nm
                qy1 = yg4(iw, nReact)
                if ((wc(iw) >= 330.) .and. (yg5(iw, nReact) > 0.)) then
                    phi1 = yg4(iw, nReact)
                    phi2 = yg5(iw, nReact)
                    phi20 = 1.-phi1
                    ak300 = ((1./phi2) - (1./phi20))/2.54e+19
                    akt = ak300*(1.+61.69*(1.-tLev(i)/300.)* &
                                 (wc(iw)/329.-1.))
                    qy2 = 1./((1./phi20) + airDen(i)*akt)
                else
                    qy2 = yg5(iw, nReact)
                end if
                qy2 = max(0., qy2)
                qy2 = min(1., qy2)
                sj(15, i, iw) = sig*qy1 !14
                sj(16, i, iw) = sig*qy2 !15
            end do
        end do

    end subroutine r10

!=============================================================================*

    subroutine r11(nw, nz, j, sj, tLev, airDen)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: =*
!=      (a)  CH3CHO + hv -> CH3 + HCO                                        =*
!=      (b)  CH3CHO + hv -> CH4 + CO                                         =*
!=      (c)  CH3CHO + hv -> CH3CO + H                                        =*
!=  Cross section:  Choice between                                           =*
!=                   (1) IUPAC 97 data, from Martinez et al.                 =*
!=                   (2) Calvert and Pitts                                   =*
!=                   (3) Martinez et al., Table 1 scanned from paper         =*
!=                   (4) KFA tabulations                                     =*
!=  Quantum yields: Choice between                                           =*
!=                   (1) IUPAC 97, pressure correction using Horowith and    =*
!=                                 Calvert, 1982                             =*
!=                   (2) NCAR data file, from Moortgat, 1986                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 11 !16,17,18

        integer, intent(IN)    :: nw
        integer, intent(IN)    :: nz
        integer, intent(INOUT) :: j
        real, intent(INOUT)    :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)
        real, intent(IN) :: airDen(nz)

        integer, parameter :: kdata = 150
        integer :: i
        real :: qy1, qy2, qy3
        real :: sig
        integer :: iw

!CH3CHO photolysis
        j = j + 1
!jlabel(j) = 'CH3CHO -> CH3 + HCO'
        j = j + 1
!jlabel(j) = 'CH3CHO -> CH4 + CO'
        j = j + 1
!jlabel(j) = 'CH3CHO -> CH3CO + H'

! combine:
        do i = 1, nz
            do iw = 1, nw - 1
                sig = yg(iw, nReact)
! quantum yields:
                qy2 = yg2(iw, nReact)
                qy3 = yg3(iw, nReact)
                qy1 = yg1(iw, nReact)
! pressure correction for channel 1, CH3 + CHO
! based on Horowitz and Calvert 1982.
                qy1 = qy1*(1.+yg4(iw, nReact))/(1.+yg4(iw, nReact)* &
                                                airDen(i)/2.465e19)
                qy1 = min(1., qy1)
                qy1 = max(0., qy1)

                sj(17, i, iw) = sig*qy1 !16

                sj(18, i, iw) = sig*qy2 !17

                sj(19, i, iw) = sig*qy3   !18
            end do
        end do

    end subroutine r11

!=============================================================================*

    subroutine r12(nw, nz, j, sj, tLev, airDen)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for C2H5CHO             =*
!=  photolysis:                                                              =*
!=           C2H5CHO + hv -> C2H5 + HCO                                             =*
!=                                                                             =*
!=  Cross section:  Choice between                                             =*
!=                     (1) IUPAC 97 data, from Martinez et al.                     =*
!=                     (2) Calvert and Pitts, as tabulated by KFA              =*
!=  Quantum yield:  IUPAC 97 recommendation                                     =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working          (I)=*
!=             wavelength grid                                                     =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in          (I)=*
!=             working wavelength grid                                             =*
!=  WC     - REAL, vector of center points of wavelength intervals in          (I)=*
!=             working wavelength grid                                             =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level          (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J           - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=             photolysis reaction defined, at each defined wavelength and     =*
!=             at each defined altitude level                                     =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=             defined                                                             =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 12 !19

        integer, intent(IN)    :: nw
        integer, intent(IN)    :: nz
        integer, intent(INOUT) :: j
        real, intent(INOUT)    :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)
        real, intent(IN) :: airDen(nz)

        integer, parameter :: kdata = 150
        integer :: i, n
        integer :: n1
        real :: x1(kdata)
        real :: y1(kdata)
        real :: qy1
        real :: sig
        integer :: ierr
        integer :: iw

!************************ C2H5CHO photolysis
! 1:  C2H5 + HCO

        j = j + 1
!jlabel(j) = 'C2H5CHO -> C2H5 + HCO'

! combine:
        do iw = 1, nw - 1
            do i = 1, nz
                sig = yg(iw, nReact)
! quantum yields:
! use Ster-Volmer pressure dependence:
                if (yg1(iw, nReact) < pzero) then
                    qy1 = 0.
                else
                    qy1 = 1./(1.+(1./yg1(iw, nReact) - 1.)* &
                                   airDen(i)/2.45e19)
                end if
                qy1 = min(qy1, 1.)
                sj(20, i, iw) = sig*qy1 !19
            end do
        end do

    end subroutine r12

!=============================================================================*

    subroutine r13(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CHOCHO         =*
!=  photolysis:                                                              =*
!=              CHOCHO + hv -> Products                                      =*
!=                                                                           =*
!=  Cross section: Choice between                                            =*
!=                  (1) Plum et al., as tabulated by IUPAC 97                =*
!=                  (2) Plum et al., as tabulated by KFA.                    =*
!=                  (3) Orlando et al.                                       =*
!=                  (4) Horowitz et al., 2001                                =*
!=  Quantum yield: IUPAC 97 recommendation                                   =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 13 !20,21

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)

        integer, parameter :: kdata = 500
        integer :: i
        real :: qyii, qyiii
        real :: sig
        integer :: iw

!************************ CHOCHO photolysis
! see review by Madronich, Chapter VII in "The Mechansims of
!  Atmospheric Oxidation of the Alkanes, Calvert et al, Oxford U.
!  Press, 2000.
! Four possible channels:
!     I     H2 + 2 CO
!     II    2 HCO
!     III   HCHO + CO
!     IV    HCO + H + CO

!  Based on that review, the following quantum yield assignments are made:

!     qy_I = 0
!     qy_II = 0.63 for radiation between 280 and 380 nm
!     qy_III = 0.2  for radiation between 280 and 380 nm
!     qy_IV = 0
! The yields for channels II and III were determined by Bauerle et al. (personal
! communication from G. Moortgat, still unpublished as of Dec 2000).
! Bauerle et al. used broad-band irradiation 280-380 nm.
! According to Zhu et al., the energetic threshold (for II) is 417 nm.  Therefore,
! here the quantum yields were set to zero for wc > 417.  Furthermore, the
! qys of Bauerle et al. were reduced to give the same J values when using full solar
! spectrum.  The reduction factor was calculated by comparing the J-values (for
! high sun) using the 380 and 417 cut offs.  The reduction factor is 7.1

        j = j + 1
!jlabel(j) = 'CHOCHO -> HCO + HCO'
        j = j + 1
!jlabel(j) = 'CHOCHO -> CH2O + CO'

! combine:
        do iw = 1, nw - 1
            sig = yg(iw, nReact)
! quantum yields:
! Use values from Bauerle, but corrected to cutoff at 417 rather than 380.
! this correction is a reduction by 7.1.
! so that qyII = 0.63/7.1  and qyIII = 0.2/7.1
            if (wc(iw) < 417.) then
                qyii = 0.089
                qyiii = 0.028
            else
                qyii = 0.
                qyiii = 0.
            end if
            do i = 1, nz
                sj(21, i, iw) = sig*qyii !20
                sj(22, i, iw) = sig*qyiii !21
            end do
        end do

    end subroutine r13

!=============================================================================*

    subroutine r14(nw, wc, nz, j, sj, tLev, airDen)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3COCHO       =*
!=  photolysis:                                                              =*
!=           CH3COCHO + hv -> CH3CO + HCO                                    =*
!=                                                                           =*
!=  Cross section: Choice between                                            =*
!=                  (1) from Meller et al., 1991, as tabulated by IUPAC 97   =*
!=                         5 nm resolution (table 1) for < 402 nm            =*
!=                         2 nm resolution (table 2) for > 402 nm            =*
!=                  (2) average at 1 nm of Staffelbach et al., 1995, and     =*
!=                      Meller et al., 1991                                  =*
!=                  (3) Plum et al., 1983, as tabulated by KFA              =*
!=                  (4) Meller et al., 1991 (0.033 nm res.), as tab. by KFA  =*
!=                  (5) Meller et al., 1991 (1.0 nm res.), as tab. by KFA    =*
!=                  (6) Staffelbach et al., 1995, as tabulated by KFA        =*
!=  Quantum yield: Choice between                                            =*
!=                  (1) Plum et al., fixed at 0.107                          =*
!=                  (2) Plum et al., divided by 2, fixed at 0.0535           =*
!=                  (3) Staffelbach et al., 0.45 for < 300 nm, 0 for > 430 nm=*
!=                      linear interp. in between                            =*
!=                  (4) Koch and Moortgat, prv. comm., 1997                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 14

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: tLev(nz)
        real, intent(IN) :: airDen(nz)

        integer, parameter :: kdata = 500
        integer :: i
        real :: qy
        real :: sig
        integer :: iw
        real :: phi0, kq

!************************ CH3COCHO photolysis
! 1:  CH3COCHO

        j = j + 1
!jlabel(j) = 'CH3COCHO -> CH3CO + HCO'

! combine:
        do iw = 1, nw - 1
            sig = yg(iw, nReact)
            do i = 1, nz
! quantum yields:
                if (mOption(14) == 1) then
                    qy = 0.107
                else if (mOption(14) == 2) then
                    qy = 0.107/2.
                else if (mOption(14) == 3) then
                    if (wc(iw) <= 300.) then
                        qy = 0.45
                    else if (wc(iw) >= 430.) then
                        qy = 0.
                    else
                        qy = 0.45 + (0 - 0.45)*(wc(iw) - 300.)/(430.-300.)
                    end if
                else if (mOption(14) == 4) then
                    if (yg1(iw, nReact) > 0.) then
                        qy = yg2(iw, nReact)/(1.+ &
                                               (airDen(i)/2.465e19)* &
                                               ((yg2(iw, nReact)/yg1(iw, nReact)) - 1.))
                    else
                        qy = 0.
                    end if
                else if (mOption(14) == 5) then
! zero pressure yield:
! 1.0 for wc < 380 nm
! 0.0 for wc > 440 nm
! linear in between:
                    phi0 = 1.-(wc(iw) - 380.)/60.
                    phi0 = min(phi0, 1.)
                    phi0 = max(phi0, 0.)
! Pressure correction: quenching coefficient, torr-1
! in air, Koch and Moortgat:
                    kq = 1.36e8*exp(-8793/wc(iw))
! in N2, Chen et al:
!               kq = 1.93e4 * EXP(-5639/wc(iw))
                    if (phi0 > 0.) then
                        if (wc(iw) >= 380. .and. wc(iw) <= 440.) then
                            qy = phi0/(phi0 + kq*airDen(i)* &
                                            760./2.456e19)
                        else
                            qy = phi0
                        end if
                    else
                        qy = 0.
                    end if
                end if
                sj(23, i, iw) = sig*qy !22
            end do
        end do

    end subroutine r14

!=============================================================================*

    subroutine r15(nw, wc, nz, j, sj, tLev, airDen)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3COCH3 photolysis=*
!=          CH3COCH3 + hv -> Products                                        =*
!=                                                                           =*
!=  Cross section:  Choice between                                           =*
!=                   (1) Calvert and Pitts                                   =*
!=                   (2) Martinez et al., 1991, alson in IUPAC 97            =*
!=                   (3) NOAA, 1998, unpublished as of 01/98                 =*
!=  Quantum yield:  Choice between                                           =*
!=                   (1) Gardiner et al, 1984                                =*
!=                   (2) IUPAC 97                                            =*
!=                   (3) McKeen et al., 1997                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 15 !23

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: airDen(nz)

        integer, parameter :: kdata = 150
        integer :: i
        real :: qy
        real :: sig
        integer :: iw
        real :: a, b, t, m, w
        real :: fco, fac
        real, parameter :: deltat = 298.-235.

!*************** CH3COCH3 photodissociation
        j = j + 1
!jlabel(j) = 'CH3COCH3 -> CH3CO + CH3'

        do iw = 1, nw - 1
            do i = 1, nz
                sig = yg(iw, nReact)
                if (mOption(15) == 3) then
                    t = 298.-tLev(i)
                    t = min(t, deltaT)
                    t = max(t, 0.)
                    sig = yg(iw, nReact)*(1.+yg2(iw, nReact)*t + &
                                           yg3(iw, nReact)*t*t)
                end if
                if (mOption(16) == 1) then
                    qy = 0.0766 + 0.09415*exp(-airDen(i)/ &
                                                   3.222e18)
                else if (mOption(16) == 2) then
                    qy = yg1(iw, nReact)
                else if (mOption(16) == 3) then
                    if (wc(iw) <= 292.) then
                        qy = 1.
                    else if (wc(iw) >= 292. .and. wc(iw) < 308.) then
                        a = -15.696 + 0.05707*wc(iw)
                        b = exp(-88.81 + 0.15161*wc(iw))
                        qy = 1./(a + b*airDen(i))
                    else if (wc(iw) >= 308. .and. wc(iw) < 337.) then
                        a = -130.2 + 0.42884*wc(iw)
                        b = exp(-55.947 + 0.044913*wc(iw))
                        qy = 1./(a + b*airDen(i))
                    else if (wc(iw) >= 337.) then
                        qy = 0.
                    end if
                    qy = max(0., qy)
                    qy = min(1., qy)
                else if (mOption(16) == 4) then
                    w = wc(iw)
                    t = tLev(i)
                    m = airDen(i)
                    call qyacet(w, t, m, fco, fac)
                    qy = max(0., fac)
                    qy = min(1., fac)
                end if
                sj(24, i, iw) = sig*qy !23
            end do
        end do

    end subroutine r15

!=============================================================================*

    subroutine r16(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3OOH photolysis: =*
!=         CH3OOH + hv -> CH3O + OH                                          =*
!=                                                                           =*
!=  Cross section: Choice between                                            =*
!=                  (1) JPL 97 recommendation (based on Vaghjiana and        =*
!=                      Ravishankara, 1989), 10 nm resolution                =*
!=                  (2) IUPAC 97 (from Vaghjiana and Ravishankara, 1989),    =*
!=                      5 nm resolution                                      =*
!=                  (3) Cox and Tyndall, 1978; only for wavelengths < 280 nm =*
!=                  (4) Molina and Arguello, 1979;  might be 40% too high    =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 16 !24

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        integer :: i
        real :: qy
        integer :: iw

!*************** CH3OOH photodissociation
        j = j + 1
!jlabel(j) = 'CH3OOH -> CH3O + OH'

! quantum yield = 1
        qy = 1.
        do iw = 1, nw - 1

            do i = 1, nz
                sj(25, i, iw) = yg(iw, nReact)*qy !24
            end do
        end do

    end subroutine r16

!=============================================================================*

    subroutine r17(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3ONO2            =*
!=  photolysis:                                                              =*
!=          CH3ONO2 + hv -> CH3O + NO2                                       =*
!=                                                                           =*
!=  Cross section: Choice between                                            =*
!=                  (1) Calvert and Pitts, 1966                              =*
!=                  (2) Talukdar, Burkholder, Hunter, Gilles, Roberts,       =*
!=                      Ravishankara, 1997                                   =*
!=                  (3) IUPAC 97, table of values for 198K                   =*
!=                  (4) IUPAC 97, temperature-dependent equation             =*
!=                  (5) Taylor et al, 1980                                   =*
!=                  (6) fit from Roberts and Fajer, 1989                     =*
!=                  (7) Rattigan et al., 1992                                =*
!=                  (8) Libuda and Zabel, 1995                               =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 17 !25

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 2000
        integer :: i
        integer :: iw
        real :: qy
        real :: sig

!*************** CH3ONO2 photodissociation
        j = j + 1
!jlabel(j) = 'CH3ONO2 -> CH3O + NO2'

! quantum yield = 1
        qy = 1.

        do iw = 1, nw - 1
            sig = yg(iw, nReact)
            do i = 1, nz
                if (mOption(18) == 2) then
                    sig = yg(iw, nReact)*exp(yg1(iw, nReact)* &
                                              (tLev(i) - 298.))
                else if (mOption(18) == 4) then
                    sig = yg(iw, nReact)*10.**(yg1(iw, nReact)* &
                                                tLev(i))
                end if
                sj(26, i, iw) = qy*sig !25
            end do
        end do

    end subroutine r17

!=============================================================================*

    subroutine r18(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for PAN photolysis:    =*
!=       PAN + hv -> Products                                                =*
!=                                                                           =*
!=  Cross section: from Talukdar et al., 1995                                =*
!=  Quantum yield: Assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 18 !26

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100

        integer :: iw
        integer :: i
        real :: qy
        real :: sig

!*************** PAN photodissociation
        j = j + 1
!jlabel(j) = 'CH3CO(OONO2) -> Products'

! quantum yield
! yet unknown, but assumed to be 1.0 (Talukdar et al., 1995)
        qy = 1.0
        do iw = 1, nw - 1
            do i = 1, nz
                sig = yg(iw, nReact)*exp(yg2(iw, nReact)* &
                                         (tLev(i) - 298.))
                sj(27, i, iw) = qy*sig !26
            end do
        end do

    end subroutine r18

!=============================================================================*

    subroutine r19(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CCl2O photolysis:  =*
!=        CCl2O + hv -> Products                                             =*
!=                                                                           =*
!=  Cross section: JPL 94 recommendation                                     =*
!=  Quantum yield: Unity (Calvert and Pitts)                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 19 !31

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!************ CCl2O photodissociation
        j = j + 1
!jlabel(j) = 'CCl2O -> Products'

!** quantum yield unity (Calvert and Pitts)
        qy = 1.
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(32, iz, iw) = qy*yg(iw, nReact) !31
            end do
        end do

    end subroutine r19

!=============================================================================*

    subroutine r20(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CCl4 photolysis:   =*
!=      CCl4 + hv -> Products                                                =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 20

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!************ CCl4 photodissociation
        j = j + 1
!jlabel(j) = 'CCl4 -> Products'

!** quantum yield assumed to be unity
        qy = 1.
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(33, iz, iw) = qy*yg(iw, nReact) !32
            end do
        end do

    end subroutine r20

!=============================================================================*

    subroutine r21(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CClFO photolysis:  =*
!=         CClFO + hv -> Products                                            =*
!=  Cross section: from JPL 97                                               =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 21 !33

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!************ CClFO photodissociation
        j = j + 1
!jlabel(j) = 'CClFO -> Products'

!** quantum yield unity
        qy = 1.
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(34, iz, iw) = qy*yg(iw, nReact) !33
            end do
        end do

    end subroutine r21

    subroutine r22(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CF2O photolysis:   =*
!=        CF2O + hv -> Products                                              =*
!=  Cross section:  from JPL 97 recommendation                               =*
!=  Quantum yield:  unity (Nolle et al.)                                     =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 22 !34

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!************ CF2O photodissociation
        j = j + 1
!jlabel(j) = 'CF2O -> Products'

!** quantum yield unity (Nolle et al.)
        qy = 1.
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(35, iz, iw) = qy*yg(iw, nReact) !34
            end do
        end do

    end subroutine r22

!=============================================================================*

    subroutine r23(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-113 photolysis:=*
!=          CF2ClCFCl2 + hv -> Products                                      =*
!=  Cross section:  from JPL 97 recommendation, linear interp. between       =*
!=                  values at 210 and 295K                                   =*
!=  Quantum yield:  assumed to be unity                                      =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 23 !35

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        real :: t
        integer :: iw
        integer :: iz
        real :: slope

!************ CF2ClCFCl2 (CFC-113) photodissociation
        j = j + 1
!jlabel(j) = 'CF2ClCFCl2 (CFC-113) -> Products'

!** quantum yield assumed to be unity
        qy = 1.

        do iz = 1, nz
            t = max(210., min(tLev(iz), 295.))
            slope = (t - 210.)/(295.-210.)
            do iw = 1, nw - 1
                sj(36, iz, iw) = qy*(yg2(iw, nReact) + &
                                      slope*(yg1(iw, nReact) - yg2(iw, nReact))) !35
            end do
        end do

    end subroutine r23

!=============================================================================*

    subroutine r24(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-144 photolysis:=*
!=              CF2ClCF2Cl + hv -> Products                                  =*
!=  Cross section: from JPL 97 recommendation, linear interp. between values =*
!=                 at 210 and 295K                                           =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 24 !36

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        real :: t
        integer :: iw
        integer :: iz
        real :: slope

!************ CF2ClCF2Cl (CFC-114) photodissociation
        j = j + 1
!jlabel(j) = 'CF2ClCF2Cl (CFC-114) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!** quantum yield assumed to be unity
        qy = 1.

        do iz = 1, nz
            t = max(210., min(tLev(iz), 295.))
            slope = (t - 210.)/(295.-210.)
            do iw = 1, nw - 1
                sj(37, iz, iw) = qy*(yg2(iw, nReact) + &
                                      slope*(yg1(iw, nReact) - yg2(iw, nReact))) !36
            end do
        end do

    end subroutine r24

!=============================================================================*

    subroutine r25(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-115 photolysis =*
!=             CF3CF2Cl + hv -> Products                                     =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 25 !37

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ CF3CF2Cl (CFC-115) photodissociation

        j = j + 1
!jlabel(j) = 'CF3CF2Cl (CFC-115) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iw = 1, nw - 1
            do iz = 1, nz
                sj(38, iz, iw) = qy*yg(iw, nReact) !37
            end do
        end do

    end subroutine r25

!=============================================================================*

    subroutine r26(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-111 photolysis =*
!=          CCl3F + hv -> Products                                           =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 26 !38

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        real :: t
        integer :: iw
        integer :: iz

!************ CCl3F (CFC-11) photodissociation
        j = j + 1
!jlabel(j) = 'CCl3F (CFC-11) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iz = 1, nz
            t = 1e-04*(tLev(iz) - 298.)
            do iw = 1, nw - 1
                sj(39, iz, iw) = qy*yg(iw, nReact)* &
                                  exp((wc(iw) - 184.9)*t) !38
            end do
        end do

    end subroutine r26

!=============================================================================*

    subroutine r27(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CFC-112 photolysis:=*
!=         CCl2F2 + hv -> Products                                           =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 27 !39

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        real :: t
        integer :: iw
        integer :: iz

!************ CCl2F2 (CFC-12) photodissociation
        j = j + 1
!jlabel(j) = 'CCl2F2 (CFC-12) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)
!*** quantum yield assumed to be unity
        qy = 1.

        do iz = 1, nz
            t = 1e-04*(tLev(iz) - 298.)
            do iw = 1, nw - 1
                sj(40, iz, iw) = qy*yg(iw, nReact)* &
                                  exp((wc(iw) - 184.9)*t) !39
            end do
        end do

    end subroutine r27

!=============================================================================*

    subroutine r28(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3Br photolysis:  =*
!=         CH3Br + hv -> Products                                            =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 28 !50

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!************ CH3Br photodissociation
! data from JPL97 (identical to 94 recommendation)
        j = j + 1
!jlabel(j) = 'CH3Br -> Products'
!*** quantum yield assumed to be unity
        qy = 1.

        do iw = 1, nw - 1
            do iz = 1, nz
                sj(51, iz, iw) = qy*yg(iw, nReact) !50
            end do
        end do

    end subroutine r28

!=============================================================================*

    subroutine r29(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3CCl3 photolysis =*
!=           CH3CCl3 + hv -> Products                                        =*
!=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
!=                 of data at 210, 250, and 295K                             =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 29

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        real :: t
        integer :: iw
        integer :: iz
        real :: slope

!************ CH3CCl3 photodissociation
        j = j + 1
!jlabel(j) = 'CH3CCl3 -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iz = 1, nz
            t = min(295., max(tLev(iz), 210.))
            if (t <= 250.) then
                slope = (t - 210.)/(250.-210.)
                do iw = 1, nw - 1
                    sj(41, iz, iw) = qy*(yg3(iw, nReact) + &
                                          slope*(yg2(iw, nReact) - yg3(iw, nReact))) !40
                end do
            else
                slope = (t - 250.)/(295.-250.)
                do iw = 1, nw - 1
                    sj(41, iz, iw) = qy*(yg2(iw, nReact) + &
                                          slope*(yg1(iw, nReact) - yg2(iw, nReact))) !40
                end do
            end if
        end do

    end subroutine r29

!=============================================================================*

    subroutine r30(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for CH3Cl photolysis:  =*
!=            CH3Cl + hv -> Products                                         =*
!=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
!=                 from values at 255, 279, and 296K                         =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 30 !30

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        real :: t
        integer :: iw
        integer :: iz
        real :: slope

!************ CH3Cl photodissociation
        j = j + 1
!jlabel(j) = 'CH3Cl -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iz = 1, nz
            t = max(255., min(tLev(iz), 296.))
            if (t <= 279.) then
                slope = (t - 255.)/(279.-255.)
                do iw = 1, nw - 1
                    sj(31, iz, iw) = qy*(yg3(iw, nReact) + &
                                          slope*(yg2(iw, nReact) - yg3(iw, nReact))) !30
                end do
            else
                slope = (t - 279.)/(296.-279.)
                do iw = 1, nw - 1
                    sj(31, iz, iw) = qy*(yg2(iw, nReact) + slope* &
                                          (yg1(iw, nReact) - yg2(iw, nReact))) !30
                end do
            end if
        end do

    end subroutine r30

!=============================================================================*

    subroutine r31(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for ClOO photolysis:   =*
!=          ClOO + hv -> Products                                            =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 31 !27

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ ClOO photodissociation
        j = j + 1
!jlabel(j) = 'ClOO -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)
        qy = 1.
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(28, iz, iw) = qy*yg(iw, nReact) !27
            end do
        end do

    end subroutine r31

!=============================================================================*

    subroutine r32(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-123 photolysis=*
!=       CF3CHCl2 + hv -> Products                                           =*
!=  Cross section: from Orlando et al., 1991                                 =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 32

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real :: qy
        real :: t
        integer :: i, iw
        integer :: iz
        real :: lambda, sum_

!************ CF3CHCl2 (HCFC-123) photodissociation
        j = j + 1
!jlabel(j) = 'CF3CHCl2 (HCFC-123) -> Products'

!*** cross section from Orlando et al., 1991
!*** quantum yield assumed to be unity

        qy = 1.

        do iw = 1, nw - 1
            lambda = wc(iw)
! use parameterization only up to 220 nm, as the error bars associated with
! the measurements beyond 220 nm are very large (Orlando, priv.comm.)
            if (lambda >= 190. .and. lambda <= 220.) then
                do iz = 1, nz
                    t = min(295., max(tLev(iz), 203.)) - &
                         tbar(nReact)
                    sum_ = 0.
                    do i = 1, 4
                        sum_ = (coeff(i, 1, nReact) + t*(coeff(i, 2, nReact) + &
                                   t*coeff(i, 3, nReact)))* &
                                (lambda - lbar)**(i - 1) + sum_
                    end do
                    sj(42, iz, iw) = qy*exp(sum_) !41
                end do
            else
                do iz = 1, nz
                    sj(42, iz, iw) = 0. !41
                end do
            end if
        end do

    end subroutine r32

!=============================================================================*

    subroutine r33(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-124 photolysis=*
!=        CF3CHFCl + hv -> Products                                          =*
!=  Cross section: from Orlando et al., 1991                                 =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 33 !42

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real :: qy
        real :: t
        integer :: i, iw
        integer :: iz
        real :: lambda, sum_

!*************************************************************
!************ CF3CHFCl (HCFC-124) photodissociation

        j = j + 1
!jlabel(j) = 'CF3CHFCl (HCFC-124) -> Products'
!*** cross section from Orlando et al., 1991
!*** quantum yield assumed to be unity

        qy = 1.

        do iw = 1, nw - 1
            lambda = wc(iw)
            if (lambda >= 190. .and. lambda <= 230.) then
                do iz = 1, nz
                    t = min(295., max(tLev(iz), 203.)) - &
                         tbar(nReact)
                    sum_ = 0.0
                    do i = 1, 4
                        sum_ = (coeff(i, 1, nReact) + t*(coeff(i, 2, nReact) + t* &
                                   coeff(i, 3, nReact)))* &
                                (lambda - lbar)**(i - 1) + sum_
                    end do
                    sj(43, iz, iw) = qy*exp(sum_) !42
                end do
            else
                do iz = 1, nz
                    sj(43, iz, iw) = 0. !42
                end do
            end if
        end do

    end subroutine r33

!=============================================================================*

    subroutine r34(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-141b          =*
!=  photolysis:                                                              =*
!=         CH3CFCl2 + hv -> Products                                         =*
!=  Cross section: from JPL97 recommendation                                 =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 34 !43

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ CH3CFCl2 (HCFC-141b) photodissociation

        j = j + 1
!jlabel(j) = 'CH3CFCl2 (HCFC-141b) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(44, iz, iw) = qy*yg(iw, nReact) !43
            end do
        end do

    end subroutine r34

!=============================================================================*

    subroutine r35(nw, wc, nz, j, sj, tLev)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-142b          =*
!=  photolysis:                                                              =*
!=          CH3CF2Cl + hv -> Products                                        =*
!=  Cross section: from Orlando et al., 1991                                 =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 35 !44

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real :: qy
        real :: t
        integer :: i, iw
        integer :: iz
        real :: lambda, sum_

!*************************************************************
!************ CH3CF2Cl (HCFC-142b) photodissociation

        j = j + 1
!jlabel(j) = 'CH3CF2Cl (HCFC-142b) -> Products'

!*** cross section from Orlando et al., 1991

!*** quantum yield assumed to be unity

        qy = 1.
        do iw = 1, nw - 1
            lambda = wc(iw)
            if (lambda >= 190. .and. lambda <= 230.) then
                do iz = 1, nz
                    t = min(295., max(tLev(iz), 203.)) - &
                         tbar(nReact)
                    sum_ = 0.
                    do i = 1, 4
                        sum_ = (coeff(i, 1, nReact) + t*(coeff(i, 2, nReact) + t* &
                                   coeff(i, 3, nReact)))*(lambda - lbar)**(i - 1) + sum_
                    end do
! offeset exponent by 40 (exp(-40.) = 4.248e-18) to prevent exp. underflow errors
! on some machines.
!             sq(j,iz,iw) = qy * EXP(sum)
                    sj(45, iz, iw) = qy*4.248e-18*exp(sum_ + &
                                                           40.0) !44
                end do
            else
                do iz = 1, nz
                    sj(45, iz, iw) = 0.0 !44
                end do
            end if
        end do

    end subroutine r35

!=============================================================================*

    subroutine r36(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-225ca         =*
!=  photolysis:                                                              =*
!=           CF3CF2CHCl2 + hv -> Products                                    =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 36 !45

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ CF3CF2CHCl2 (HCFC-225ca) photodissociation

        j = j + 1
!jlabel(j) = 'CF3CF2CHCl2 (HCFC-225ca) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iw = 1, nw - 1
            do iz = 1, nz
                sj(46, iz, iw) = qy*yg(iw, nReact) !45
            end do
        end do

    end subroutine r36

!=============================================================================*

    subroutine r37(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-225cb         =*
!=  photolysis:                                                              =*
!=          CF2ClCF2CHFCl + hv -> Products                                   =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 37 !46

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ CF2ClCF2CHFCl (HCFC-225cb) photodissociation

        j = j + 1
!jlabel(j) = 'CF2ClCF2CHFCl (HCFC-225cb) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iw = 1, nw - 1
            do iz = 1, nz
                sj(47, iz, iw) = qy*yg(iw, nReact) !46
            end do
        end do

    end subroutine r37

!=============================================================================*

    subroutine r38(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HCFC-22 photolysis =*
!=          CHClF2 + hv -> Products                                          =*
!=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
!=                 from values at 210, 230, 250, 279, and 295 K              =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 38 !47

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        real :: t
        integer :: iw
        integer :: iz
        real :: slope

!*************************************************************
!************ CHClF2 (HCFC-22) photodissociation

        j = j + 1
!jlabel(j) = 'CHClF2 (HCFC-22) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iz = 1, nz
            t = min(295., max(tLev(iz), 210.))
            if (t <= 230.) then
                slope = (t - 210.)/(230.-210.)
                do iw = 1, nw - 1
                    sj(48, iz, iw) = qy*(yg5(iw, nReact) + slope* &
                                          (yg4(iw, nReact) - yg5(iw, nReact))) !47
                end do
            else if (t <= 250.) then
                slope = (t - 230.)/(250.-230.)
                do iw = 1, nw - 1
                    sj(48, iz, iw) = qy*(yg4(iw, nReact) + slope* &
                                          (yg3(iw, nReact) - yg4(iw, nReact))) !47
                end do
            else if (t <= 270.) then
                slope = (t - 250.)/(270.-250.)
                do iw = 1, nw - 1
                    sj(48, iz, iw) = qy*(yg3(iw, nReact) + slope* &
                                          (yg2(iw, nReact) - yg3(iw, nReact))) !47
                end do
            else
                slope = (t - 270.)/(295.-270.)
                do iw = 1, nw - 1
                    sj(48, iz, iw) = qy*(yg2(iw, nReact) + slope* &
                                          (yg1(iw, nReact) - yg2(iw, nReact))) !47
                end do
            end if
        end do

    end subroutine r38

!=============================================================================*

    subroutine r39(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for HO2 photolysis:    =*
!=          HO2 + hv -> OH + O                                               =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed shape based on work by Lee, 1982; normalized      =*
!=                 to unity at 248 nm                                        =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 39 !9

        integer, intent(IN)      :: nw
        real, intent(IN)         :: wc(kw)
        integer, intent(IN)      :: nz
        integer, intent(INOUT)   :: j
        real, intent(IN)         :: tLev(nz)
        real, intent(INOUT)      :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ HO2 photodissociation
        j = j + 1
!jlabel(j) = 'HO2 -> OH + O'
!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield:  absolute quantum yield has not been reported yet, but
!***                 Lee measured a quantum yield for O(1D) production at 248
!***                 nm that was 15 time larger than at 193 nm
!*** here:  a quantum yield of unity is assumed at 248 nm and beyond, for
!***        shorter wavelengths a linear decrease with lambda is assumed

        do iw = 1, nw - 1
            if (wc(iw) >= 248.) then
                qy = 1.
            else
                qy = 1./15.+(wc(iw) - 193.)*(14./15.)/(248.-193.)
                qy = max(qy, 0.)
            end if
            do iz = 1, nz
                sj(10, iz, iw) = qy*yg(iw, nReact) !9
            end do
        end do

    end subroutine r39

!=============================================================================*

    subroutine r40(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) Halon-1202 photolysis: =*
!=         CF2Br2 + hv -> Products                                           =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: unity (Molina and Molina)                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 40 !54

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ CF2Br2 (Halon-1202) photodissociation

        j = j + 1
!jlabel(j) = 'CF2Br2 (Halon-1202) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield unity (Molina and Molina)
        qy = 1.

        do iw = 1, nw - 1
            do iz = 1, nz
                sj(55, iz, iw) = qy*yg(iw, nReact) !54
            end do
        end do

    end subroutine r40

!=============================================================================*

    subroutine r41(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for Halon-1211         =*
!=  photolysis:                                                              =*
!=           CF2ClBr + hv -> Products                                        =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 41 !55

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ CF2BrCl (Halon-1211) photodissociation

        j = j + 1
!jlabel(j) = 'CF2BrCl (Halon-1211) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)
!*** quantum yield assumed to be unity
        qy = 1.

        do iw = 1, nw - 1
            do iz = 1, nz
                sj(56, iz, iw) = qy*yg(iw, nReact) !55
            end do
        end do

    end subroutine r41

!=============================================================================*

    subroutine r42(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for Halon-1301         =*
!=  photolysis:                                                              =*
!=         CF3Br + hv -> Products                                            =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 42 !52

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ CF3Br (Halon-1301) photodissociation

        j = j + 1
!jlabel(j) = 'CF3Br (Halon-1301) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iw = 1, nw - 1
            do iz = 1, nz
                sj(53, iz, iw) = qy*yg(iw, nReact) !52
            end do
        end do

    end subroutine r42

!=============================================================================*

    subroutine r43(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for Halon-2402         =*
!=  photolysis:                                                              =*
!=           CF2BrCF2Br + hv -> Products                                     =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 43  !53

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy
        integer :: iw
        integer :: iz

!*************************************************************
!************ CF2BrCF2Br (Halon-2402) photodissociation

        j = j + 1
!jlabel(j) = 'CF2BrCF2Br (Halon-2402) -> Products'

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)

!*** quantum yield assumed to be unity
        qy = 1.

        do iw = 1, nw - 1
            do iz = 1, nz
                sj(54, iz, iw) = qy*yg(iw, nReact) !53
            end do
        end do

    end subroutine r43

!=============================================================================*

    subroutine r44(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for N2O photolysis:    =*
!=              N2O + hv -> N2 + O(1D)                                       =*
!=  Cross section: from JPL 97 recommendation                                =*
!=  Quantum yield: assumed to be unity, based on Greenblatt and Ravishankara =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 44 !8

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real :: qy
        real :: a, b!, c
        real :: a0, a1, a2, a3, a4
        real :: b0, b1, b2, b3
        real :: t
        integer :: iw, iz
        real :: lambda

!*************************************************************
!************ N2O photodissociation
        j = j + 1
!jlabel(j) = 'N2O -> N2 + O(1D)'

!*** cross sections according to JPL97 recommendation (identical to 94 rec.)
!*** see file dataj1/abs/N2O_jpl94.abs for detail

        a0 = 68.21023
        a1 = -4.071805
        a2 = 4.301146e-02
        a3 = -1.777846e-04
        a4 = 2.520672e-07

        b0 = 123.4014
        b1 = -2.116255
        b2 = 1.111572e-02
        b3 = -1.881058e-05

!*** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
!*** Ravishankara), so quantum yield of O(1D) is assumed to be unity
        qy = 1.

        do iw = 1, nw - 1
            lambda = wc(iw)
            if (lambda >= 173. .and. lambda <= 240.) then
                do iz = 1, nz
                    t = max(194., min(tLev(iz), 320.))
                    a = (((a4*lambda + a3)*lambda + a2)*lambda + a1)*lambda + a0
                    b = (((b3*lambda + b2)*lambda + b1)*lambda + b0)
                    b = (t - 300.)*exp(b)
                    sj(9, iz, iw) = qy*exp(a + b) !8
                end do
            else
                do iz = 1, nz
                    sj(9, iz, iw) = 0. !8
                end do
            end if
        end do

    end subroutine r44

!=============================================================================*

    subroutine r45(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for ClONO2 photolysis: =*
!=        ClONO2 + hv -> Products                                            =*
!=                                                                           =*
!=  Cross section: JPL 97 recommendation                                     =*
!=  Quantum yield: JPL 97 recommendation                                     =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 45 !28,29

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 150

        real :: qy1, qy2
        real :: xs
        integer :: iw
        integer :: iz

!************ ClONO2 photodissociation

        j = j + 1
!jlabel(j) = 'ClONO2 -> Cl + NO3'

!** cross sections from JPL97 recommendation

        do iw = 1, nw - 1
!** quantum yields (from jpl97)
            if (wc(iw) < 308.) then
                qy1 = 0.6
            else if ((wc(iw) >= 308) .and. (wc(iw) <= 364.)) then
                qy1 = 7.143e-3*wc(iw) - 1.6
            else if (wc(iw) > 364.) then
                qy1 = 1.0
            end if
            qy2 = 1.-qy1
! compute T-dependent cross section
            do iz = 1, nz
                xs = yg1(iw, nReact)*(1.+yg2(iw, nReact)* &
                                      (tLev(iz) - 296) + &
                                      yg3(iw, nReact)*(tLev(iz) - 296)* &
                                      (tLev(iz) - 296))
                sj(29, iz, iw) = qy1*xs !28
                sj(30, iz, iw) = qy2*xs !29
            end do
        end do

        j = j + 1
!jlabel(j) = 'ClONO2 -> ClO + NO2'

    end subroutine r45

!=============================================================================*

    subroutine r46(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for BrONO2 photolysis: =*
!=        BrONO2 + hv -> Products                                            =*
!=                                                                           =*
!=  Cross section: JPL 03 recommendation                                     =*
!=  Quantum yield: JPL 03 recommendation                                     =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 46 !48,49

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 100
        real :: qy1, qy2
        integer :: iw
        integer :: iz

!************ BrONO2 photodissociation

        j = j + 1
!jlabel(j) = 'BrONO2 -> BrO + NO2'
        j = j + 1
!jlabel(j) = 'BrONO2 -> Br + NO3'

!** cross sections from JPL03 recommendation

!** quantum yields (from jpl97)

        qy1 = 0.71
        qy2 = 0.29
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(49, iz, iw) = qy1*yg1(iw, nReact) !48
                sj(50, iz, iw) = qy2*yg1(iw, nReact) !49
            end do
        end do

    end subroutine r46

!=============================================================================*

    subroutine r47(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide product (cross section) x (quantum yield) for Cl2 photolysis:    =*
!=        Cl2 + hv -> 2 Cl                                                   =*
!=                                                                           =*
!=  Cross section: JPL 97 recommendation                                     =*
!=  Quantum yield: 1     (Calvert and Pitts, 1966)                           =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 47 !56

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 150
        real :: qy
        integer :: iz, iw

!************ CL2 photodissociation

        j = j + 1
!jlabel(j) = 'Cl2 -> Cl + Cl'

!** cross sections from JPL97 recommendation (as tab by Finlayson-Pitts
! and Pitts, 1999.

!** quantum yield = 1 (Calvert and Pitts, 1966)

        qy = 1.
        do iw = 1, nw - 1
            do iz = 1, nz
                sj(57, iz, iw) = qy*yg(iw, nReact) !56
            end do
        end do

    end subroutine r47

!=============================================================================*

    subroutine r101(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH2(OH)CHO     =*
!=  (glycolaldehye, hydroxy acetaldehyde) photolysis:                        =*
!=           CH2(OH)CHO + hv -> Products                                     =*
!=                                                                           =*
!=  Cross section from                                                       =*
!= The Atmospheric Chemistry of Glycolaldehyde, C. Bacher, G. S. Tyndall     =*
!= and J. J. Orlando, J. Atmos. Chem., 39 (2001) 171-189.                    =*
!=                                                                           =*
!=  Quantum yield about 50%                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 101 !57

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 300
        integer :: i
        real :: qy
        integer :: iw

!************************ CH2(OH)CHO photolysis
! 1:  CH2(OH)CHO

        j = j + 1
!jlabel(j) = 'CH2(OH)CHO -> Products'

! combine:
        qy = 0.5

        do iw = 1, nw - 1
            do i = 1, nz
                sj(58, i, iw) = yg(iw, nReact)*qy !57
            end do
        end do

    end subroutine r101

!=============================================================================*

    subroutine r102(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3COCOCH3     =*
!=  (biacetyl) photolysis:                                                   =*
!=           CH3COCOCH3 + hv -> Products                                     =*
!=                                                                           =*
!=  Cross section from either                                                =*
!= 1.  Plum et al., Environ. Sci. Technol., Vol. 17, No. 8, 1983, p.480      =*
!= 2.  Horowitz et al., J. Photochem Photobio A, 146, 19-27, 2001.           =*
!=                                                                           =*
!=  Quantum yield =0.158                                                     =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 102 !58

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 300
        integer :: i
        real :: qy
        integer :: iw

!************************ CH3COCOCH3 photolysis
! 1:  CH3COCOCH3
        j = j + 1
!jlabel(j) = 'CH3COCOCH3 -> Products'

! quantum yield from Plum et al.

        qy = 0.158

        do iw = 1, nw - 1
            do i = 1, nz
                sj(59, i, iw) = yg(iw, nReact)*qy !58
            end do
        end do

    end subroutine r102

!=============================================================================*

    subroutine r103(nw, wc, nz, j, sj, tLev, airDen)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3COCHCH2     =*
!=  Methyl vinyl ketone photolysis:                                          =*
!=           CH3COCHCH2 + hv -> Products                                     =*
!=                                                                           =*
!=  Cross section from                                                       =*
!= W. Schneider and G. K. Moorgat, priv. comm, MPI Mainz 1989 as reported by =*
!= Roeth, E.-P., R. Ruhnke, G. Moortgat, R. Meller, and W. Schneider,        =*
!= UV/VIS-Absorption Cross Sections and QUantum Yields for Use in            =*
!= Photochemistry and Atmospheric Modeling, Part 2: Organic Substances,      =*
!= Forschungszentrum Julich, Report Jul-3341, 1997.                          =*
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 103 !59

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real, intent(IN) :: airDen(nz)

        integer, parameter :: kdata = 20000
        integer :: i
        real :: qy
        integer :: iw

!************************ CH3COCHCH2 photolysis

        j = j + 1
!jlabel(j) = 'CH3COCHCH2 -> Products'

! quantum yield from
! Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
! and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
! J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
! depends on pressure and wavelength, set upper limit to 1.0

        do iw = 1, nw - 1
            do i = 1, nz
                qy = exp(-0.055*(wc(iw) - 308.))/(5.5 + 9.2e-19* &
                                                  airDen(i))
                qy = min(qy, 1.)
                sj(60, i, iw) = yg(iw, nReact)*qy !59
            end do
        end do

    end subroutine r103

!=============================================================================*

    subroutine r104(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH2C(CH3)CHO   =*
!=  methacrolein        photolysis:                                          =*
!=           CH2C(CH3)CHO + hv -> Products                                   =*
!=                                                                           =*
!=  Cross section from                                                       =*
!= R. Meller, priv. comm, MPI Mainz 1990 as reported by =*
!= Roeth, E.-P., R. Ruhnke, G. Moortgat, R. Meller, and W. Schneider,        =*
!= UV/VIS-Absorption Cross Sections and QUantum Yields for Use in            =*
!= Photochemistry and Atmospheric Modeling, Part 2: Organic Substances,      =*
!= Forschungszentrum Julich, Report Jul-3341, 1997.                          =*
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 104 !60

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 20000
        integer :: i
        real :: qy
        integer :: iw

!************************ CH2C(CH3)CHO photolysis

        j = j + 1
!jlabel(j) = 'CH2C(CH3)CHO -> Products'

! quantum yield from
! Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
! and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
! J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
!   Upper limit, quantum yield < 0.01

        qy = 0.01

        do iw = 1, nw - 1
            do i = 1, nz
                sj(61, i, iw) = yg(iw, nReact)*qy !60
            end do
        end do
!PRINT *,'React 104. sj=',sj(:,61,:,:)
    end subroutine r104

!=============================================================================*

    subroutine r105(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3COCO(OH)    =*
!=  pyruvic acid        photolysis:                                          =*
!=           CH3COCO(OH) + hv -> Products                                    =*
!=                                                                           =*
!=  Cross section from                                                       =*
!= Horowitz, A., R. Meller, and G. K. Moortgat, The UV-VIS absorption cross  =*
!= section of the a-dicarbonyl compounds: pyruvic acid, biacetyl, and        =*
!= glyoxal. J. Photochem. Photobiol. A:Chemistry, v.146, pp.19-27, 2001.     =*
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 105 !61

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 20000
        integer :: i
        real :: qy
        integer :: iw

!************************ CH3COCO(OH) photolysis

        j = j + 1
!jlabel(j) = 'CH3COCO(OH) -> Products'

! quantum yield  = 1

        qy = 1.

        do iw = 1, nw - 1
            do i = 1, nz
                sj(62, i, iw) = yg(iw, nReact)*qy !61
            end do
        end do

    end subroutine r105

!=============================================================================*

    subroutine r106(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3CH2ONO2     =*
!=  ethyl nitrate       photolysis:                                          =*
!=           CH3CH2ONO2 + hv -> CH3CH2O + NO2                                =*
!=                                                                           =*
!= Absorption cross sections of several organic from                         =*
!= Talukdar, R. K., J. B. Burkholder, M. Hunter, M. K. Gilles,               =*
!= J. M Roberts, and A. R. Ravishankara, Atmospheric fate of several         =*
!= alkyl nitrates, J. Chem. Soc., Faraday Trans., 93(16) 2797-2805, 1997.    =*
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 106 !62

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 200
        integer :: i
        real :: qy, sig
        integer :: iw

!************************ CH3CH2ONO2 photolysis
        j = j + 1
!jlabel(j) = 'CH3CH2ONO2 -> CH3CH2O + NO2'

! quantum yield  = 1

        qy = 1.

        do iw = 1, nw - 1
            do i = 1, nz
                sig = yg1(iw, nReact)*exp(yg2(iw, nReact)* &
                                          (tLev(i) - 298.))
                sj(63, i, iw) = sig*qy !62
            end do
        end do

    end subroutine r106

!=============================================================================*

    subroutine r107(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for CH3CHONO2CH3   =*
!=  isopropyl nitrate   photolysis:                                          =*
!=           CH3CHONO2CH3 + hv -> CH3CHOCH3 + NO2                            =*
!=                                                                           =*
!= Absorption cross sections of several organic from                         =*
!= Talukdar, R. K., J. B. Burkholder, M. Hunter, M. K. Gilles,               =*
!= J. M Roberts, and A. R. Ravishankara, Atmospheric fate of several         =*
!= alkyl nitrates, J. Chem. Soc., Faraday Trans., 93(16) 2797-2805, 1997.    =*
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 107 !63

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 200
        integer :: i
        real :: qy, sig
        integer :: iw

!************************ CH3CHONO2CH3 photolysis
        j = j + 1
!jlabel(j) = 'CH3CHONO2CH3 -> CH3CHOCH3 + NO2'

! quantum yield  = 1

        qy = 1.

        do iw = 1, nw - 1
            do i = 1, nz
                sig = yg1(iw, nReact)*exp(yg2(iw, nReact)* &
                                          (tLev(i) - 298.))
                sj(64, i, iw) = sig*qy !63
            end do
        end do

    end subroutine r107

!=============================================================================*

    subroutine r108(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for                =*
!=   nitroxy ethanol CH2(OH)CH2(ONO2) + hv -> CH2(OH)CH2(O.) + NO2           =*
!=                                                                           =*
!=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
!=    sections of organic nitrates of potential atmospheric importance and   =*
!=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
!=    1989.
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 108 !64

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real :: qy, sig
        integer :: iw, i
        real :: a, b, c

!************************ CH2(OH)CH2(ONO2) photolysis

        j = j + 1
!jlabel(j) = 'CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2'
! coefficients from Roberts and Fajer 1989, over 270-306 nm
        a = -2.359e-3
        b = 1.2478
        c = -210.4

! quantum yield  = 1

        qy = 1.

        do iw = 1, nw - 1
            if (wc(iw) >= 270. .and. wc(iw) <= 306.) then
                sig = exp(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            else
                sig = 0.
            end if
            do i = 1, nz
                sj(65, i, iw) = sig*qy !64
            end do
        end do

    end subroutine r108

!=============================================================================*

    subroutine r109(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for                =*
!=   nitroxy acetone CH3COCH2(ONO2) + hv -> CH3COCH2(O.) + NO2               =*
!=                                                                           =*
!=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
!=    sections of organic nitrates of potential atmospheric importance and   =*
!=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
!=    1989.
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 109 !65

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real :: qy, sig
        integer :: iw, i
        real :: a, b, c

!************************ CH3COCH2(ONO2) photolysis

        j = j + 1
!jlabel(j) = 'CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2'

! coefficients from Roberts and Fajer 1989, over 284-335 nm
        a = -1.365e-3
        b = 0.7834
        c = -156.8

! quantum yield  = 1
        qy = 1.
        do iw = 1, nw - 1
            if (wc(iw) >= 284. .and. wc(iw) <= 335.) then
                sig = exp(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            else
                sig = 0.
            end if
            do i = 1, nz
                sj(66, i, iw) = sig*qy !65
            end do
        end do

    end subroutine r109

!=============================================================================*

    subroutine r110(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for                =*
!=  t-butyl nitrate C(CH3)3(ONO2) + hv -> C(CH3)(O.) + NO2                   =*
!=                                                                           =*
!=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
!=    sections of organic nitrates of potential atmospheric importance and   =*
!=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
!=    1989.
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 110 !66

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        real :: qy, sig
        integer :: iw, i
        real :: a, b, c

!************************ C(CH3)3(ONO2) photolysis
        j = j + 1
!jlabel(j) = 'C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2'

! coefficients from Roberts and Fajer 1989, over 270-330 nm
        a = -0.993e-3
        b = 0.5307
        c = -115.5

! quantum yield  = 1
        qy = 1.
        do iw = 1, nw - 1
            if (wc(iw) >= 270. .and. wc(iw) <= 330.) then
                sig = exp(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            else
                sig = 0.
            end if
            do i = 1, nz
                sj(67, i, iw) = sig*qy !66
            end do
        end do

    end subroutine r110

!=============================================================================*

    subroutine r111(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for ClOOCl         =*
!=  ClO dimer           photolysis:                                          =*
!=           ClOOCl + hv -> Cl + ClOO                                        =*
!=                                                                           =*
!=  Cross section from  JPL2002                                              =*
!=                                                                           =*
!=  Quantum yield assumed unity                                              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 111 !67

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 20000
        integer :: i
        real :: qy
        integer :: iw

!************************ ClOOCl photolysis
! from JPL-2002

        j = j + 1
!jlabel(j) = 'ClOOCl -> Cl + ClOO'
! quantum yield  = 1
        qy = 1.

        do iw = 1, nw - 1
            do i = 1, nz
                sj(68, i, iw) = yg(iw, nReact)*qy !67
            end do
        end do

    end subroutine r111

!=============================================================================*

    subroutine r112(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for hydroxyacetone =*
!=  CH2(OH)COCH3        photolysis:                                          =*
!=           CH2(OH)COCH3  -> CH3CO + CH2OH
!=                         -> CH2(OH)CO + CH3                                =*
!=                                                                           =*
!=  Cross section from Orlando et al. (1999)                                 =*
!=                                                                           =*
!=  Quantum yield assumed 0.325 for each channel (J. Orlando, priv.comm.2003)=*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 112 !68,69

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 20000
        integer :: i
        real :: qy
        integer :: iw

!************************ CH2(OH)COCH3 photolysis
! from Orlando et al. 1999

        j = j + 1
!jlabel(j) = 'CH2(OH)COCH3 -> CH3CO + CH2(OH)'
        j = j + 1
!jlabel(j) = 'CH2(OH)COCH3 -> CH2(OH)CO + CH3'

! Total quantum yield  = 0.65, equal for each of the two channels
        qy = 0.325
        do iw = 1, nw - 1
            do i = 1, nz
                sj(69, i, iw) = yg(iw, nReact)*qy !68
                sj(70, i, iw) = yg(iw, nReact)*qy   !69
            end do
        end do

    end subroutine r112

!=============================================================================*

    subroutine r113(nw, wc, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for HOBr           =*
!=  HOBr -> OH + Br                                                          =*
!=  Cross section from JPL 2003                                              =*
!=  Quantum yield assumed unity as in JPL2003                                =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 113 !70

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wc(kw)
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer :: i
        real :: qy, sig
        integer :: iw

!************************ HOBr photolysis
! from JPL2003

        j = j + 1
!jlabel(j) = 'HOBr -> OH + Br'

        qy = 1.
        do iw = 1, nw - 1
            sig = 24.77*exp(-109.80*(log(284.01/wc(iw)))**2) + &
                  12.22*exp(-93.63*(log(350.57/wc(iw)))**2) + &
                  2.283*exp(-242.40*(log(457.38/wc(iw)))**2)
            sig = sig*1.e-20
            if (wc(iw) < 250. .or. wc(iw) > 550.) sig = 0.
            do i = 1, nz
                sj(71, i, iw) = sig*qy !70
            end do
        end do

    end subroutine r113

!=============================================================================*

    subroutine r114(nw, nz, j, sj)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for BrO            =*
!=  BrO -> Br + O                                                            =*
!=  Cross section from JPL 2003                                              =*
!=  Quantum yield assumed unity as in JPL2003                                =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 114 !71

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer :: i
        integer :: iw
        real :: qy!, dum!, yg(kw)

!************************ HOBr photolysis
! from JPL2003

        j = j + 1
!jlabel(j) = 'BrO -> Br + O'
        qy = 1.
        do iw = 1, nw - 1
            do i = 1, nz
                sj(72, i, iw) = yg(iw, nReact)*qy !71
            end do
        end do

    end subroutine r114

!=============================================================================*

    subroutine r115(nw, nz, j, sj, tLev)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Provide the product (cross section) x (quantum yield) for BrO            =*
!=  Br2 -> Br + Br                                                           =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: nReact = 115 !72

        integer, intent(IN)             :: nw
        integer, intent(IN)             :: nz
        integer, intent(INOUT)          :: j
        real, intent(IN)                :: tLev(nz)
        real, intent(INOUT)             :: sj(kj, nz, kw)

        integer, parameter :: kdata = 50
        integer :: i
        integer :: iw
        real :: qy !, yg(kw)

!************************ Br2 photolysis

        j = j + 1
!jlabel(j) = 'Br2 -> Br + Br'

! Absorption cross section from:
! Seery, D.J. and D. Britton, The continuous absorption spectra of chlorine,
! bromine, bromine chloride, iodine chloride, and iodine bromide, J. Phys.
! Chem. 68, p. 2263 (1964).

        qy = 1.
        do iw = 1, nw - 1
            do i = 1, nz
                sj(73, i, iw) = yg(iw, nReact)*qy !72
            end do
        end do

    end subroutine r115

    subroutine atrim(a, aa, n)
        implicit none

! Trim blanks from character string
! Input: a
! Output:  aa, n
! Internal: i

        character(LEN=6), intent(IN)            :: a
        character(LEN=6), intent(OUT)           :: aa
        integer, intent(OUT)                     :: n

        integer :: i

        aa = ' '
        n = 0
        do i = 1, 6
            if (a(i:i) /= ' ') then
                n = n + 1
                aa(n:n) = a(i:i)
            end if
        end do

    end subroutine atrim

    subroutine setz(nz, nj, coef, adjcoe, iang, tLev, tcO3)
!-----------------------------------------------------------------------------*
!=  ZENITH = 0
!-----------------------------------------------------------------------------*
!=  INPUT:                                                                   =*
!=
!=  nz = height   = 121
!=  nj = species  = 73
!=  cz = overhead total ozone
!=  tlev=Temperature
!=  coef=POL coef(MZ,MS,MP)
!=
!=  OUTPUT:
!=
!=  ADJCOE - REAL, coross section adjust coefficients                        =*
!=                                                                           =*
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  08/2005 XUEXI                                                            =*
!-----------------------------------------------------------------------------*
!= This program is free software;  you can redistribute it and/or modify     =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                 =*
!= To obtain a copy of the GNU General Public License, write to:             =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= To contact the authors, please mail to:                                   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu                                            =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, parameter :: mz = 5
        integer, parameter :: ms = 73
        integer, parameter :: mp = 5

        integer, intent(IN) :: iang
        integer, intent(IN) :: nz
        integer, intent(IN) :: nj
        real, intent(IN)    :: coef(mz, ms, mp)
        real, intent(IN)    :: tLev
        real, intent(IN)    :: tco3(nz)

        real, intent(OUT)   :: adjcoe(kj, nz)

        integer :: k, js, jp

        real :: xz(nz), c(5, kj)

        real :: c0, c1, c2
        real :: adjin
        integer :: ij
        real :: tt

!===============================================
!  SET-UP
!===============================================

        do k = 1, nz
            do ij = 1, nj
                adjcoe(ij, k) = 1.00
            end do
        end do

        if (iang > 5) return

        do k = 1, nz
            xz(k) = tcO3(k)*1.e-18
        end do

        tt = tLev/281.0

        do js = 1, ms
            do jp = 1, mp
                c(jp, js) = coef(iang, js, jp)
            end do
        end do

! CAL ADJ COEF

!All species except tropospheric
        do ij = 1, 27
            adjin = 1.00
            call calcoe(nz, ij, c, xz, adjin, adjcoe)
        end do
        do ij = 58, 73
            adjin = 1.00
            call calcoe(nz, ij, c, xz, adjin, adjcoe)
        end do
!(2) O3 -> O2 + O(1D)
!----------------------------------------------------------------------
!      Temperature Modification
!      T0.9 (1.3) T0.95(1.25)  T1.0(1.2)  T1.15(1.18)  T1.1(1.16)
!----------------------------------------------------------------------
        c0 = ca0(iang)
        c1 = ca1(iang)
        c2 = ca2(iang)
        adjin = c0 + c1*tt + c2*tt*tt + fatSum(iAng)
        call calcoe(nz, 2, c, xz, adjin, adjcoe)
! 11 H2O2 -> 2 OH
!----------------------------------------------------------------------
!      Temperature Modification
!      T0.9(0.95)  T0.95(0.975)    T1.0(1.0)  T1.15(0.105)   T1.1(1.10)
!----------------------------------------------------------------------
        c0 = cb0(iAng)
        c1 = cb1(iAng)
        c2 = cb2(iAng)
        adjin = c0 + c1*tt + c2*tt*tt
        call calcoe(nz, 11, c, xz, adjin, adjcoe)

    end subroutine setz

    subroutine setno2(nz, nw, no2xs, zLevel, no2col, cAir, dtNo2)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Set up an altitude profile of NO2 molecules, and corresponding absorption=*
!=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
!=  that allows scaling of the entire profile to a given overhead NO2        =*
!=  column amount.                                                           =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NO2NEW - REAL, overhead NO2 column amount (molec/cm^2) to which       (I)=*
!=           profile should be scaled.  If NO2NEW < 0, no scaling is done    =*
!=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
!=           grid                                                            =*
!=  Z      - REAL, specified altitude working grid (km)                   (I)=*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  NO2XS  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
!=           each specified wavelength                                       =*
!=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
!=  DTNO2  - REAL, optical depth due to NO2 absorption at each            (O)=*
!=           specified altitude at each specified wavelength                 =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nz
        integer, intent(IN) :: nw
        real, intent(IN)    :: no2xs(kw)
        real, intent(IN)    :: zLevel(nz)
        real, intent(IN)    :: cAir(nz)
        real, intent(IN)    :: no2col
        real, intent(OUT)   :: dtNo2(nz, kw)

        integer, parameter :: kdata = 51
        real :: cz(nz)
! nitrogen dioxide profile data:
        real :: zd(kdata), no2(kdata)
        real :: cd(kdata)
        real :: hscale
        real :: colold
        real :: scale_
        real :: sno2

        integer :: i, l, nd

! Example:  set to 1 ppb in lowest 1 km, set to zero above that.
! - do by specifying concentration at 3 altitudes.

        nd = 3
        zd(1) = 0.
        no2(1) = 1.*2.69e10
        zd(2) = 1.
        no2(2) = 1.*2.69e10
        zd(3) = zd(2)*1.000001
        no2(3) = 10./largest
! compute column increments (alternatively, can specify these directly)
        do i = 1, nd - 1
            cd(i) = (no2(i + 1) + no2(i))*1.e5*(zd(i + 1) - zd(i))/2.
        end do
! Include exponential tail integral from top level to infinity.
! fold tail integral into top layer
! specify scale height near top of data (use ozone value)
        hscale = 4.50e5
        cd(nd - 1) = cd(nd - 1) + hscale*no2(nd)
!********** end data input.

! Compute column increments and total column on standard z-grid.
        call inter3(nz, zLevel, cz, nd, zd, cd, 1)
!*** Scaling of vertical profile by ratio of new to old column:
! If old column is near zero (less than 1 molec cm-2),
! use constant mixing ratio profile (nominal 1 ppt before scaling)
! to avoid numerical problems when scaling.

        if (fsum(nz - 1, cz) < 1.) then
            do i = 1, nz - 1
                cz(i) = 1.e-12*cAir(i)
            end do
        end if
        colold = fsum(nz - 1, cz)
        scale_ = 2.687e16*no2col/colold
        do i = 1, nz
            if (i > nz - 1) cycle
            cz(i) = cz(i)*scale_
        end do

!***********************************
! calculate optical depth for each layer.  Output: dtno2(nz,kw)
        do l = 1, nw - 1
            sno2 = no2xs(l)
            do i = 1, nz
                if (i > nz - 1) cycle
                dtNo2(i, l) = cz(i)*sno2
            end do
        end do
!_______________________________________________________________________

    end subroutine setno2
!=============================================================================*

    subroutine setso2(nz, nw, so2xs, zLevel, so2col, cAir, dtSo2)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Set up an altitude profile of SO2 molecules, and corresponding absorption=*
!=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
!=  that allows scaling of the entire profile to a given overhead SO2        =*
!=  column amount.                                                           =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  SO2NEW - REAL, overhead SO2 column amount (molec/cm^2) to which       (I)=*
!=           profile should be scaled.  If SO2NEW < 0, no scaling is done    =*
!=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
!=           grid                                                            =*
!=  Z      - REAL, specified altitude working grid (km)                   (I)=*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  SO2XS  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
!=           each specified wavelength                                       =*
!=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
!=  DTSO2  - REAL, optical depth due to SO2 absorption at each            (O)=*
!=           specified altitude at each specified wavelength                 =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nw
        integer, intent(IN) :: nz
        real, intent(IN)    :: so2xs(kw)
        real, intent(IN)    :: zLevel(nz)
        real, intent(IN)    :: cAir(nz)
        real, intent(IN)    :: so2col
        real, intent(OUT)   :: dtSo2(nz, kw)

        integer, parameter :: kdata = 51
        real :: cz(nz)
! sulfur dioxide profile data:
        real :: zd(kdata), so2(kdata)
        real :: cd(kdata)
        real :: hscale
        real :: colold
        real :: scale_
        real :: sso2

        integer :: i, l, nd

! Example:  set to 1 ppb in lowest 1 km, set to zero above that.
! - do by specifying concentration at 3 altitudes.

        nd = 3
        zd(1) = 0.
        so2(1) = 1.*2.69e10
        zd(2) = 1.
        so2(2) = 1.*2.69e10
        zd(3) = zd(2)*1.000001
        so2(3) = 10./largest
! compute column increments (alternatively, can specify these directly)
        do i = 1, nd - 1
            cd(i) = (so2(i + 1) + so2(i))*1.e5*(zd(i + 1) - zd(i))/2.
        end do
! Include exponential tail integral from top level to infinity.
! fold tail integral into top layer
! specify scale height near top of data (use ozone value)
        hscale = 4.50e5
        cd(nd - 1) = cd(nd - 1) + hscale*so2(nd)
!********** end data input.

! Compute column increments on standard z-grid.
        call inter3(nz, zLevel, cz, nd, zd, cd, 1)
!*** Scaling of vertical profile by ratio of new to old column:
! If old column is near zero (less than 1 molec cm-2),
! use constant mixing ratio profile (nominal 1 ppt before scaling)
! to avoid numerical problems when scaling.

        if (fsum(nz - 1, cz) < 1.) then
            do i = 1, nz - 1
                cz(i) = 1.e-12*cAir(i)
            end do
        end if
        colold = fsum(nz - 1, cz)
        scale_ = 2.687e16*so2col/colold
        do i = 1, nz
            if (i > nz - 1) cycle
            cz(i) = cz(i)*scale_
        end do

!***********************************
! calculate sulfur optical depth for each layer, with optional temperature
! correction.  Output, dtso2(nz,kw)

        do l = 1, nw - 1
            sso2 = so2xs(l)
            do i = 1, nz
! Leaving this part in in case i want to interpolate between
! the 221K and 298K data.
!            IF ( wl(l) .GT. 240.5  .AND. wl(l+1) .LT. 350. ) THEN
!               IF (tlay(i) .LT. 263.) THEN
!                  sso2 = s221(l) + (s263(l)-s226(l)) / (263.-226.) *
!     $                 (tlay(i)-226.)
!               ELSE
!                  sso2 = s263(l) + (s298(l)-s263(l)) / (298.-263.) *
!     $              (tlay(i)-263.)
!               ENDIF
!            ENDIF
                if (i > nz - 1) cycle
                dtSo2(i, l) = cz(i)*sso2
            end do
        end do
!_______________________________________________________________________

    end subroutine setso2

    subroutine sphers(nz, zLevel, sza, nid, dsdh)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate slant path over vertical depth ds/dh in spherical geometry.    =*
!=  Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model =*
!=  for computing the radiation field available for photolysis and heating   =*
!=  at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  Z       - REAL, specified altitude working grid (km)                  (I)=*
!=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
!=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
!=            when travelling from the top of the atmosphere to layer i;     =*
!=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
!=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
!=            travelling from the top of the atmosphere to layer i;          =*
!=            NID(i), i = 0..NZ-1                                            =*
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  double precision fix for shallow layers - Julia Lee-Taylor Dec 2000      =*
!-----------------------------------------------------------------------------*
!= This program is free software;  you can redistribute it and/or modify     =*
!= it under the terms of the GNU General Public License as published by the  =*
!= Free Software Foundation;  either version 2 of the license, or (at your   =*
!= option) any later version.                                                =*
!= The TUV package is distributed in the hope that it will be useful, but    =*
!= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
!= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
!= License for more details.                                                 =*
!= To obtain a copy of the GNU General Public License, write to:             =*
!= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
!-----------------------------------------------------------------------------*
!= To contact the authors, please mail to:                                   =*
!= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
!= send email to:  sasha@ucar.edu                                            =*
!-----------------------------------------------------------------------------*
!= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nz
        real, intent(IN) :: zLevel(nz)
        real, intent(IN) :: sza
        integer, intent(OUT) :: nid(0:nz)
        real, intent(OUT) :: dsdh(0:nz, nz)

        real :: re
        real :: ze(nz)
        double precision :: zenrad
        double precision :: rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
        integer :: i, j, k
        integer :: id

        integer :: nlayer
        real :: zd(0:nz - 1)

        zenrad = sza*rpd
! number of layers:
!nlayer = nz - 1

! include the elevation above sea level to the radius of the earth:
        re = radius + zLevel(1)
! correspondingly z changed to the elevation above earth surface:
        do k = 1, nz
            ze(k) = zLevel(k) - zLevel(1)
        end do

! inverse coordinate of z
        zd(0) = ze(nz)
        do k = 1, nz - 1
            zd(k) = ze(nz - k)
        end do

! initialize dsdh(i,j), nid(i)
        do i = 0, nz
            nid(i) = 0
            do j = 1, nz
                dsdh(i, j) = 0.
            end do
        end do

! calculate ds/dh of every layer
        do i = 0, nz - 1
            rpsinz = (re + zd(i))*sin(zenrad)
            nid(i) = -1
            if (.not. (sza > 90.0) .and. &
                (rpsinz < re)) then
!ELSE
! Find index of layer in which the screening height lies
                id = i
                if (sza > 90.0) then
                    do j = 1, nz - 1
                        if ((rpsinz < (zd(j - 1) + re)) .and. &
                            (rpsinz >= (zd(j) + re))) id = j
                    end do
                end if
                do j = 1, id
                    sm = 1.0
                    if (j == id .and. id == i .and. &
                        sza > 90.0) sm = -1.0
                    rj = re + zd(j - 1)
                    rjp1 = re + zd(j)
                    dhj = zd(j - 1) - zd(j)
                    if (dhj == 0.0) cycle !LFR for teste
                    ga = rj*rj - rpsinz*rpsinz
                    gb = rjp1*rjp1 - rpsinz*rpsinz
                    if (ga < 0.0) ga = 0.0
                    if (gb < 0.0) gb = 0.0
                    if (id > i .and. j == id) then
                        dsj = sqrt(ga)
                    else
                        dsj = sqrt(ga) - sm*sqrt(gb)
                    end if
                    dsdh(i, j) = dsj/dhj
                end do
                nid(i) = id
            end if
        end do

    end subroutine sphers

!=============================================================================*
    subroutine airmas(nz, cAir, scol, vcol, nid, dsdh)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate vertical and slant air columns, in spherical geometry, as a    =*
!=  function of altitude.                                                    =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
!=            grid                                                           =*
!=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
!=            when travelling from the top of the atmosphere to layer i;     =*
!=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
!=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
!=            travelling from the top of the atmosphere to layer i;          =*
!=            NID(i), i = 0..NZ-1                                            =*
!=  VCOL    - REAL, output, vertical air column, molec cm-2, above level iz  =*
!=  SCOL    - REAL, output, slant air column in direction of sun, above iz   =*
!=            also in molec cm-2                                             =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN) :: nz
        integer, intent(IN)  :: nid(0:nz)
        real, intent(IN)    :: cAir(nz)
        real, intent(IN)    :: dsdh(0:nz, nz)
        real, intent(OUT)   :: scol(nz)
        real, intent(OUT)   :: vcol(nz)

        integer :: id, j
        real :: sum, vsum

! calculate vertical and slant column from each level:
! work downward
        vsum = 0.
        do id = 0, nz - 1
            vsum = vsum + cAir(nz - id)
            vcol(nz - id) = vsum
            sum = 0.
            if (nid(id) < 0) then
                sum = largest
            else
! single pass layers:
                do j = 1, min(nid(id), id)
                    sum = sum + cAir(nz - j)*dsdh(id, j)
                end do
! double pass layers:
                do j = min(nid(id), id) + 1, nid(id)
                    sum = sum + 2.0*cAir(nz - j)*dsdh(id, j)
                end do
            end if
            scol(nz - id) = sum
        end do

    end subroutine airmas

    subroutine swbiol(nw, wl, wc, j, s, label)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Create or read various weighting functions, e.g. biological action       =*
!=  spectra, UV index, etc.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of central wavelength of wavelength intervals    I)=*
!=           in working wavelength grid                                      =*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  S      - REAL, value of each defined weighting function at each       (O)=*
!=           defined wavelength                                              =*
!=  LABEL  - CHARACTER*50, string identifier for each weighting function  (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN)               :: nw
        real, intent(IN)                  :: wl(kw)
        real, intent(IN)                  :: wc(kw)
        integer, intent(INOUT)            :: j
        real, intent(INOUT)               :: s(ks, kw)
        character(LEN=50), intent(INOUT) :: label(ks)

        integer, parameter :: kdata = 1000

! internal:
        real :: x1(kdata)
        real :: y1(kdata)
        real :: yg(kw)

!LFR>   REAL :: fery, futr
!EXTERNAL fery, futr
        integer :: i, iw, n

        integer :: ierr

        integer :: idum
        real :: dum1, dum2
        real :: em, a, b, c
!REAL :: sum

        real :: a0, a1, a2, a3

!_______________________________________________________________________

!******** Photosynthetic Active Radiation (400 < PAR < 700 nm)
! conversion to micro moles m-2 s-1:
!  s = s * (1e6/6.022142E23)(w/1e9)/(6.626068E-34*2.99792458E8)

        j = j + 1

        wbioStart = j
        label(j) = 'PAR, 400-700 nm, umol m-2 s-1'
        print *, j, label(j), wbioStart!; call flush (6)
        do iw = 1, nw - 1
            if (wc(iw) > 400. .and. wc(iw) < 700.) then
                s(j, iw) = 8.36e-3*wc(iw)
            else
                s(j, iw) = 0.
            end if
        end do

!********* unity raf constant slope:

        j = j + 1
        label(j) = 'Exponential decay, 14 nm/10'
        do iw = 1, nw - 1
            s(j, iw) = 10.**(-(wc(iw) - 300.)/14.)
        end do

!*********** DNA damage action spectrum
! from: Setlow, R. B., The wavelengths in sunlight effective in
!       producing skin cancer: a theoretical analysis, Proceedings
!       of the National Academy of Science, 71, 3363 -3366, 1974.
! normalize to unity at 300 nm
! Data read from original hand-drawn plot by Setlow
! received from R. Setlow in May 1995
! data is per quantum (confirmed with R. Setlow in May 1995).
! Therefore must put on energy basis if irradiance is is energy
! (rather than quanta) units.

        j = j + 1
        label(j) = 'DNA damage, in vitro (Setlow, 1974)'
        open (UNIT=kin, FILE=trim(files(21)%fileName), STATUS='old')
        do i = 1, 11
            read (kin, *)
        end do
        n = 55
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)/2.4e-02*x1(i)/300.
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, label(j)
            stop
        end if

        do iw = 1, nw - 1
            s(j, iw) = yg(iw)
        end do

!******** skin cancer in mice,  Utrecht/Phildelphia study
!from de Gruijl, F. R., H. J. C. M. Sterenborg, P. D. Forbes,
!     R. E. Davies, C. Cole, G. Kelfkens, H. van Weelden, H. Slaper,
!     and J. C. van der Leun, Wavelength dependence of skin cancer
!     induction by ultraviolet irradiation of albino hairless mice,
!     Cancer Res., 53, 53-60, 1993.
! normalize at 300 nm.

        j = j + 1
        label(j) = 'SCUP-mice (de Gruijl et al., 1993)'
        do iw = 1, nw - 1
            s(j, iw) = futr(wc(iw))/futr(300.)
        end do

!********** Utrecht/Philadelphia mice spectrum corrected for humans skin.
! From de Gruijl, F.R. and J. C. van der Leun, Estimate of the wavelength
! dependency of ultraviolet carcinogenesis and its relevance to the
! risk assessment of a stratospheric ozone depletion, Health Phys., 4,
! 317-323, 1994.

        j = j + 1
        label(j) = 'SCUP-human (de Gruijl and van der Leun, 1994)'
        open (UNIT=kin, FILE=trim(files(22)%fileName), STATUS='old')
        n = 28
        do i = 1, n
            read (kin, *) x1(i), y1(i)
        end do

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, label(j)
            stop
        end if

        do iw = 1, nw - 1
            s(j, iw) = yg(iw)
        end do
        close (kin)

!**************** CIE standard human erythema action spectrum
!from:
! McKinlay, A. F., and B. L. Diffey, A reference action spectrum for
! ultraviolet induced erythema in human skin, in Human Exposure to
! Ultraviolet Radiation: Risks and Regulations, W. R. Passchler
! and B. F. M. Bosnajokovic, (eds.), Elsevier, Amsterdam, 1987.

        j = j + 1
        label(j) = 'CIE human erythema (McKinlay and Diffey, 1987)'
        do iw = 1, nw - 1
            s(j, iw) = fery(wc(iw))
        end do

!**************** UV index (Canadian - WMO/WHO)
! from:
! Report of the WMO Meeting of experts on UV-B measurements, data quality
! and standardization of UV indices, World Meteorological Organization
! (WMO), report No. 95, Geneva, 1994.
! based on the CIE erythema weighting, multiplied by 40.

        j = j + 1
        label(j) = 'UV index (WMO, 1994)'
        do iw = 1, nw - 1
            s(j, iw) = 40.*fery(wc(iw))
        end do

!************ Human erythema - Anders et al.
! from:
! Anders, A., H.-J. Altheide, M. Knalmann, and H. Tronnier,
! Action spectrum for erythema in humands investigated with dye lasers,
! Photochem. and Photobiol., 61, 200-203, 1995.
! for skin types II and III, Units are J m-2.

        j = j + 1
        label(j) = 'Erythema, humans (Anders et al., 1995)'
        open (UNIT=kin, FILE=trim(files(23)%fileName), STATUS='old')
        do i = 1, 5
            read (kin, *)
        end do
        n = 28
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = 1./y1(i)
        end do

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, label(j)
            stop
        end if

        do iw = 1, nw - 1
            s(j, iw) = yg(iw)
        end do
        close (kin)

!******** 1991-92 ACGIH threshold limit values
! from
! ACGIH, 1991-1992 Threshold Limit Values, American Conference
!  of Governmental and Industrial Hygienists, 1992.

        j = j + 1
        label(j) = 'Occupational TLV (ACGIH, 1992)'
        open (UNIT=kin, FILE=trim(files(24)%fileName), STATUS='old')
        n = 56
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)
        end do

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, label(j)
            stop
        end if

        do iw = 1, nw - 1
            s(j, iw) = yg(iw)
        end do
        close (kin)

!******** phytoplankton, Boucher et al. (1994)
! from Boucher, N., Prezelin, B.B., Evens, T., Jovine, R., Kroon, B., Moline, M.A.,
! and Schofield, O., Icecolors '93: Biological weighting function for the ultraviolet
!  inhibition  of carbon fixation in a natural antarctic phytoplankton community,
! Antarctic Journal, Review 1994, pp. 272-275, 1994.
! In original paper, value of b and m (em below are given as positive.  Correct values
! are negative. Also, limit to positive values.

        j = j + 1
        label(j) = 'Phytoplankton (Boucher et al., 1994)'
        a = 112.5
        b = -6.223e-01
        c = 7.670e-04
        em = -3.17e-06
        do iw = 1, nw - 1
            if (wc(iw) > 290. .and. wc(iw) < 400.) then
                s(j, iw) = em + exp(a + b*wc(iw) + c*wc(iw)*wc(iw))
            else
                s(j, iw) = 0.
            end if
            s(j, iw) = max(s(j, iw), 0.)
        end do

!******** phytoplankton, Cullen et al.
! Cullen, J.J., Neale, P.J., and Lesser, M.P., Biological weighting function for the
!  inhibition of phytoplankton photosynthesis by ultraviolet radiation, Science, 25,
!  646-649, 1992.
! phaeo

        j = j + 1
        label(j) = 'Phytoplankton, phaeo (Cullen et al., 1992)'
        open (UNIT=kin, FILE=trim(files(25)%fileName), STATUS='old')
        n = 106
        do i = 1, n
            read (kin, *) idum, dum1, dum2, y1(i)
            x1(i) = (dum1 + dum2)/2.
        end do

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, label(j)
            stop
        end if

        do iw = 1, nw - 1
            s(j, iw) = yg(iw)
        end do
        close (kin)

! proro

        j = j + 1
        label(j) = 'Phytoplankton, proro (Cullen et al., 1992)'
        open (UNIT=kin, FILE=trim(files(26)%fileName), STATUS='old')
        n = 100
        do i = 1, n
            read (kin, *) idum, dum1, dum2, y1(i)
            x1(i) = (dum1 + dum2)/2.
        end do

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, label(j)
            stop
        end if

        do iw = 1, nw - 1
            s(j, iw) = yg(iw)
        end do
        close (kin)

!*** Damage to lens of pig eyes, from
! Oriowo, M. et al. (2001). Action spectrum for in vitro
! UV-induced cataract using whole lenses. Invest. Ophthalmol. & Vis. Sci. 42,
! 2596-2602.  For pig eyes. Last two columns computed by L.O.Bjorn.

        j = j + 1
        label(j) = 'Cataract, pig (Oriowo et al., 2001)'
        open (UNIT=kin, FILE=trim(files(27)%fileName), STATUS='old')
        do i = 1, 7
            read (kin, *)
        end do
        n = 18
        do i = 1, n
            read (kin, *) x1(i), dum1, dum1, y1(i)
        end do

! extrapolation to 400 nm (has very little effect on raf):
!      do i = 1, 30
!         n = n + 1
!         x1(n) = x1(n-1) + 1.
!         y1(n) = 10**(5.7666 - 0.0254*x1(n))
!      enddo

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, label(j)
            stop
        end if

        do iw = 1, nw - 1
            s(j, iw) = yg(iw)
        end do
        close (kin)

!***** Plant damage - Caldwell 1971
!  Caldwell, M. M., Solar ultraviolet radiation and the growth and
! development of higher plants, Photophysiology 6:131-177, 1971.

        j = j + 1
        label(j) = 'Plant damage (Caldwell, 1971)'

! Fit to Caldwell (1971) data by
! Green, A. E. S., T. Sawada, and E. P. Shettle, The middle
! ultraviolet reaching the ground, Photochem. Photobiol., 19,
! 251-259, 1974.

        do iw = 1, nw - 1
            s(j, iw) = 2.628*(1.-(wc(iw)/313.3)**2)*exp(-(wc(iw) - 300.)/31.08)
            if (s(j, iw) < 0. .or. wc(iw) > 313.) then
                s(j, iw) = 0.
            end if
        end do

! Alternative fit to Caldwell (1971) by
! Micheletti, M. I. and R. D. Piacentini, Photochem. Photobiol.,
! 76, pp.?, 2002.

        a0 = 570.25
        a1 = -4.70144
        a2 = 0.01274
        a3 = -1.13118e-5
        do iw = 1, nw - 1
            s(j, iw) = a0 + a1*wc(iw) + a2*wc(iw)**2 + a3*wc(iw)**3
            if (s(j, iw) < 0. .or. wc(iw) > 313.) then
                s(j, iw) = 0.
            end if
        end do

!***** Plant damage - Flint & Caldwell 2003
!  Flint, S. D. and M. M. Caldwell, A biological spectral weigthing
!  function for ozone depletion research with higher plants, Physiologia
!  Plantorum, in press, 2003.
!  Data available to 366 nm

        j = j + 1
        label(j) = 'Plant damage (Flint & Caldwell, 2003)'

        do iw = 1, nw - 1
            s(j, iw) = exp(4.688272*exp(-exp(0.1703411*(wc(iw) - 307.867)/1.15)) + &
                           ((390 - wc(iw))/121.7557 - 4.183832))

! put on per joule (rather than per quantum) basis:

            s(j, iw) = s(j, iw)*wc(iw)/300.

            if (s(j, iw) < 0. .or. wc(iw) > 366.) then
                s(j, iw) = 0.
            end if

        end do
        wbioEnd = j
!***************************************************************
!***************************************************************

!_______________________________________________________________________

        if (j > ks) stop '1001'
!_______________________________________________________________________

    end subroutine swbiol

    subroutine swphys(nw, wl, wc, j, s, label)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Create or read various spectral weighting functions, physically-based    =*
!=  e.g. UV-B, UV-A, visible ranges, instrument responses, etc.              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  WC     - REAL, vector of central wavelength of wavelength intervals    I)=*
!=           in working wavelength grid                                      =*
!=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
!=  S      - REAL, value of each defined weighting function at each       (O)=*
!=           defined wavelength                                              =*
!=  LABEL  - CHARACTER*40, string identifier for each weighting function  (O)=*
!=           defined                                                         =*
!-----------------------------------------------------------------------------*
        implicit none

        integer, intent(IN)             :: nw
        real, intent(IN)                :: wl(kw)
        real, intent(IN)                :: wc(kw)
        integer, intent(OUT)            :: j
        real, intent(OUT)               :: s(ks, kw)
        character(LEN=50), intent(OUT) :: label(ks)

        integer, parameter :: kdata = 1000

! internal:
        real :: x1(kdata)
        real :: y1(kdata)
        real :: yg(kw)

        integer :: i, iw, n

        integer :: ierr

!INTEGER :: idum
!REAL :: dum1, dum2
!REAL :: em, a, b, c
        real :: sum

!_______________________________________________________________________

        j = 0

!******** UV-B (280-315 nm)

        j = j + 1
        label(j) = 'UV-B, 280-315 nm'
        do iw = 1, nw - 1
            if (wc(iw) > 280. .and. wc(iw) < 315.) then
                s(j, iw) = 1.
            else
                s(j, iw) = 0.
            end if
        end do

!******** UV-B* (280-320 nm)

        j = j + 1
        label(j) = 'UV-B*, 280-320 nm'
        do iw = 1, nw - 1
            if (wc(iw) > 280. .and. wc(iw) < 320.) then
                s(j, iw) = 1.
            else
                s(j, iw) = 0.
            end if
        end do

!******** UV-A (315-400 nm)

        j = j + 1
        label(j) = 'UV-A, 315-400 nm'
        do iw = 1, nw - 1
            if (wc(iw) > 315. .and. wc(iw) < 400.) then
                s(j, iw) = 1.
            else
                s(j, iw) = 0.
            end if
        end do

!******** visible+ (> 400 nm)

        j = j + 1
        label(j) = 'vis+, > 400 nm'
        do iw = 1, nw - 1
            if (wc(iw) > 400.) then
                s(j, iw) = 1.
            else
                s(j, iw) = 0.
            end if
        end do

!*********  Gaussian transmission functions

        j = j + 1
        label(j) = 'Gaussian, 305 nm, 10 nm FWHM'
        sum = 0.
        do iw = 1, nw - 1
!srf -avoiding floating-point exception for single precision
            if ((log(2.)*((wc(iw) - 305.)/(5.))**2) < 80.) then
                s(j, iw) = exp(-(log(2.)*((wc(iw) - 305.)/(5.))**2))
            else
                s(j, iw) = 0.0
            end if

            sum = sum + s(j, iw)
        end do
        do iw = 1, nw - 1
            s(j, iw) = s(j, iw)/sum
        end do

        j = j + 1
        label(j) = 'Gaussian, 320 nm, 10 nm FWHM'
        sum = 0.
        do iw = 1, nw - 1
            if ((log(2.)*((wc(iw) - 320.)/(5.))**2) < 80.) then
                s(j, iw) = exp(-(log(2.)*((wc(iw) - 320.)/(5.))**2))
            else
                s(j, iw) = 0.0
            end if

            sum = sum + s(j, iw)
        end do
        do iw = 1, nw - 1
            s(j, iw) = s(j, iw)/sum
        end do

        j = j + 1
        label(j) = 'Gaussian, 340 nm, 10 nm FWHM'
        sum = 0.
        do iw = 1, nw - 1
            if ((log(2.)*((wc(iw) - 340.)/(5.))**2) < 80.) then
                s(j, iw) = exp(-(log(2.)*((wc(iw) - 340.)/(5.))**2))
            else
                s(j, iw) = 0.0
            end if
            sum = sum + s(j, iw)
        end do
        do iw = 1, nw - 1
            s(j, iw) = s(j, iw)/sum
        end do

        j = j + 1
        label(j) = 'Gaussian, 380 nm, 10 nm FWHM'
        sum = 0.
        do iw = 1, nw - 1
            if ((log(2.)*((wc(iw) - 380.)/(5.))**2) < 80.) then
                s(j, iw) = exp(-(log(2.)*((wc(iw) - 380.)/(5.))**2))
            else
                s(j, iw) = 0.0
            end if
            sum = sum + s(j, iw)
        end do
        do iw = 1, nw - 1
            s(j, iw) = s(j, iw)/sum
        end do

!********* RB Meter, model 501
!  private communication, M. Morys (Solar Light Co.), 1994.
! From: morys@omni.voicenet.com (Marian Morys)
! Received: from acd.ucar.edu by sasha.acd.ucar.edu (AIX 3.2/UCB 5.64/4.03)
!          id AA17274; Wed, 21 Sep 1994 11:35:44 -0600

        j = j + 1
        label(j) = 'RB Meter, model 501'
        open (UNIT=kin, FILE=trim(files(20)%fileName), STATUS='old')
        n = 57
        do i = 1, n
            read (kin, *) x1(i), y1(i)
        end do

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg, n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, label(j)
            stop
        end if

        do iw = 1, nw - 1
            s(j, iw) = yg(iw)
        end do
        close (kin)

!***************************************************************
!***************************************************************

!_______________________________________________________________________

        if (j > ks) stop '1001'
!_______________________________________________________________________

    end subroutine swphys

!_______________________________________________________________________
    real function refrac(w, airden)
! input vacuum wavelength, nm and air density, molec cm-3
! output refractive index for standard air
! (dry air at 15 deg. C, 101.325 kPa, 0.03% CO2)

        real, intent(IN) :: w
        real, intent(IN) :: airden

! internal

        real :: sig, dum

! from CRC Handbook, originally from Edlen, B., Metrologia, 2, 71, 1966.
! valid from 200 nm to 2000 nm
! beyond this range, use constant value

        sig = 1.e3/w

        if (w < 200.) sig = 1.e3/200.
        if (w > 2000.) sig = 1.e3/2000.

        dum = 8342.13 + 2406030./(130.-sig*sig) + 15997./(38.9 - sig*sig)

! adjust to local air density

        dum = dum*airden/(2.69e19*273.15/288.15)

! index of refraction:

        refrac = 1.+1.e-8*dum

    end function refrac

!_______________________________________________________________________
    subroutine wshift(mrefr, n, w, airden)

! Shift wavelength scale between air and vacuum.
! if mrefr = 1, shift input waveelengths in air to vacuum.
! if mrefr = -1, shift input wavelengths from vacuum to air
! if any other number, don't shift

        integer, intent(IN)   :: n
        integer, intent(IN)   :: mrefr
        real, intent(IN)      :: airden
        real, intent(INOUT)   :: w(n)

! internal
        integer :: i

        if (mrefr == 1) then
            do i = 1, n
                w(i) = w(i)/refrac(w(i), airden)
            end do
        else if (mrefr == -1) then
            do i = 1, n
                w(i) = w(i)*refrac(w(i), airden)
            end do
        end if

    end subroutine wshift

    subroutine readpol(coef)
!-----------------------------------------------------------------------------*
!=  PARAMETERS:  (XUEXI)                                                     =*
!=  coef(mz,ms,mp) - REAL, pol coefficent for FTUV                           =*
!=  mz = 5 zenith anagle 0,20,40,60,80                                       =*
!=  ms = 73 species
!=  mp = 5 pol coeff  0,1,2,3,4
        implicit none

        real, intent(OUT) :: coef(mz, ms, mp)

        integer :: iz, is, ip
        open (81, FILE=trim(files(1)%fileName))

        do iz = 1, mz
            read (81, *)
            do is = 1, ms
                read (81, *)
                do ip = 1, mp
                    read (81, *) coef(iz, is, ip)
                end do
            end do
        end do

    end subroutine readpol
    subroutine ChangeKData(kData)
        implicit none
        integer, intent(IN) :: kData

        if (allocated(x1)) deallocate (x1)
        if (allocated(y1)) deallocate (y1)
        if (allocated(x2)) deallocate (x2)
        if (allocated(y2)) deallocate (y2)
        if (allocated(x3)) deallocate (x3)
        if (allocated(y3)) deallocate (y3)
        if (allocated(x4)) deallocate (x4)
        if (allocated(y4)) deallocate (y4)
        if (allocated(x5)) deallocate (x5)
        if (allocated(y5)) deallocate (y5)
        if (allocated(x)) deallocate (x)
        if (allocated(y)) deallocate (y)

        allocate (x1(kData))
        allocate (y1(kData))
        allocate (x2(kData))
        allocate (y2(kData))
        allocate (x3(kData))
        allocate (y3(kData))
        allocate (x4(kData))
        allocate (y4(kData))
        allocate (x5(kData))
        allocate (y5(kData))
        allocate (x(kData))
        allocate (y(kData))

    end subroutine ChangeKData

    subroutine ReadAll(nw, wl)
        implicit none

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read the O3 cross section                                 =*
!=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
!=  175.439-847.5 nm) and:                                                   =*
!=  For Hartley and Huggins bands, use temperature-dependent values from         =*
!=  Molina, L. T., and M. J. Molina, Absolute absorption cross sections         =*
!=  of ozone in the 185- to 350-nm wavelength range, J. Geophys. Res.,         =*
!=  vol. 91, 14501-14508, 1986.                                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                                 =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
!=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
!=  XS     - REAL, cross section (cm^2) for O3                              (O)=*
!=           at each defined wavelength and each defined altitude level         =*
!-----------------------------------------------------------------------------*

        integer, intent(IN)  :: nw
        real, intent(IN)         :: wl(kw)

        integer :: iw

        integer :: kdata
        integer :: nReact

        integer :: i, n, k, idum
        real    :: a1, a2, dum

        integer :: ierr, irev
        integer :: n1, n2, n3, n4, n5
        real    :: yglocal(kw)
!REAL :: qyVal
        integer :: irow, icol
        real :: xs270
        real :: xs280
        character(LEN=120) :: inline

!Reading xreference file
        open (UNIT=kin, FILE=trim(files(145)%fileName), STATUS='old')
!PRINT *,'Reading crossreference file for reactions'
        do i = 1, 3
            read (kin, *)
        end do
        do i = 1, kj
            read (kin, FMT='(A57,I3.3,1X,L1)') inline, xRef(i), doReaction(i)
        end do
        close (UNIT=kin)

        kData = 250
        call ChangeKData(kData)
        select case (mOption(1))
        case (1)
!----------------------------------------------------------
! cross sections from WMO 1985 Ozone Assessment
! from 175.439 to 847.500 nm
! use value at 273 K
!PRINT *,'Reading cross sections from WMO 1985 Ozone Assessment'
!PRINT *,'LFR->'//trim(files(139)%fileName)//'!'
            open (UNIT=kin, FILE=trim(files(139)%fileName), STATUS='old')
            do i = 1, 3
                read (kin, *)
            end do
            n = 158
            do i = 1, n
                read (kin, *) idum, a1, a2, dum, dum, dum, dum, y1(i)
                x1(i) = (a1 + a2)/2.
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 cross section - WMO'
                stop
            end if

            do iw = 1, nw - 1
                mm_o3xs(iw) = yglocal(iw)
            end do

! For Hartley and Huggins bands, use temperature-dependent values from
! Molina, L. T., and M. J. Molina, Absolute absorption cross sections
! of ozone in the 185- to 350-nm wavelength range,
! J. Geophys. Res., vol. 91, 14501-14508, 1986.

            open (UNIT=kin, FILE=trim(files(140)%fileName), STATUS='old')
            do i = 1, 5
                read (kin, *)
            end do
            n1 = 220
            n2 = 220
            n3 = 220
            do i = 1, n1
                read (kin, *) x1(i), y1(i), y2(i), y3(i)
                x2(i) = x1(i)
                x3(i) = x1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - 226K Molina'
                stop
            end if
            do iw = 1, nw - 1
                s226(iw) = yglocal(iw)*1.e-20
            end do

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 0., 0.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - 263K Molina'
                stop
            end if
            do iw = 1, nw - 1
                s263(iw) = yglocal(iw)*1.e-20
            end do

            call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 0., 0.)
            call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n3, x3, y3, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - 298K Molina'
                stop
            end if
            do iw = 1, nw - 1
                s298(iw) = yglocal(iw)*1.e-20
            end do
        case (2)
!----------------------------------------------------------
! cross sections from WMO 1985 Ozone Assessment
! from 175.439 to 847.500 nm
! use value at 273 K
            print *, 'Reading cross sections from WMO 1985 Ozone Assessment'
            open (UNIT=kin, FILE=trim(files(141)%fileName), STATUS='old')
            do i = 1, 3
                read (kin, *)
            end do
            n = 158
            do i = 1, n
                read (kin, *) idum, a1, a2, dum, dum, dum, dum, y1(i)
                x1(i) = (a1 + a2)/2.
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 cross section - WMO'
                stop
            end if

            do iw = 1, nw - 1
                mm_o3xs(iw) = yglocal(iw)
            end do

!=  For Hartley and Huggins bands, use temperature-dependent values from     =*
!=  Malicet et al., J. Atmos. Chem.  v.21, pp.263-273, 1995.                   =*

            open (UNIT=kin, FILE=trim(files(142)%fileName), STATUS='old')
            do i = 1, 1
                read (kin, *)
            end do
            n1 = 15001
            n2 = 15001
            n3 = 15001
            n4 = 15001

            do i = 1, n1
                read (kin, *) x1(i), y1(i), y2(i), y3(i), y4(i)
                x2(i) = x1(i)
                x3(i) = x1(i)
                x4(i) = x1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - 295K Malicet'
                stop
            end if
            do iw = 1, nw - 1
                s295(iw) = yglocal(iw)
            end do

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 0., 0.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - 243K Malicet'
                stop
            end if
            do iw = 1, nw - 1
                s243(iw) = yglocal(iw)
            end do

            call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 0., 0.)
            call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n3, x3, y3, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - 228K Malicet'
                stop
            end if
            do iw = 1, nw - 1
                s228(iw) = yglocal(iw)
            end do

            call addpnt(x4, y4, kdata, n4, x4(1)*(1.-deltax), 0.)
            call addpnt(x4, y4, kdata, n4, 0., 0.)
            call addpnt(x4, y4, kdata, n4, x4(n4)*(1.+deltax), 0.)
            call addpnt(x4, y4, kdata, n4, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n4, x4, y4, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - 218K Malicet'
                stop
            end if
            do iw = 1, nw - 1
                s218(iw) = yglocal(iw)
            end do

        case (3)
!----------------------------------------------------------
! cross sections from WMO 1985 Ozone Assessment
! from 175.439 to 847.500 nm
! use value at 273 K

            open (UNIT=kin, FILE=trim(files(143)%fileName), STATUS='old')
            do i = 1, 3
                read (kin, *)
            end do
            n = 158
            do i = 1, n
                read (kin, *) idum, a1, a2, dum, dum, dum, dum, y1(i)
                x1(i) = (a1 + a2)/2.
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 cross section - WMO'
                stop
            end if

            do iw = 1, nw - 1
                mm_o3xs(iw) = yglocal(iw)
            end do

! For Hartley and Huggins bands, use temperature-dependent values from
!  Bass et al.

            open (UNIT=kin, FILE=trim(files(144)%fileName), STATUS='old')
            do i = 1, 8
                read (kin, *)
            end do
            n1 = 1915
            n2 = 1915
            n3 = 1915
            do i = 1, n1
                read (kin, *) x1(i), y1(i), y2(i), y3(i)
                x2(i) = x1(i)
                x3(i) = x1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - c0 Bass'
                stop
            end if
            do iw = 1, nw - 1
                c0(iw) = yglocal(iw)
            end do

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 0., 0.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - c1 Bass'
                stop
            end if
            do iw = 1, nw - 1
                c1(iw) = yglocal(iw)
            end do

            call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 0., 0.)
            call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 1.e+38, 0.)
            call inter2(nw, wl, yglocal, n3, x3, y3, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'O3 xsect - c2 Bass'
                stop
            end if
            do iw = 1, nw - 1
                c2(iw) = yglocal(iw)
            end do

        end select

!===============================================================
!R01
!===============================================================
! read parameters from JPL'97
        nReact = 01
        kData = 500
        call ChangeKData(kData)
        select case (mOption(2))
        case (kjpl97)
            open (UNIT=kin, FILE=trim(files(28)%fileName), STATUS='old')
            read (kin, *)
            read (kin, *)
            read (kin, *)
            n1 = 21
            n2 = n1
            do i = 1, n1
                read (kin, *) x1(i), y1(i), y2(i)
                x2(i) = x1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n1, 0., y1(1))
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), y1(n1))
            call addpnt(x1, y1, kdata, n1, 1.e+38, y1(n1))
            call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R01' ! jlabel(j)
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
            call addpnt(x2, y2, kdata, n2, 0., y2(1))
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), y2(n2))
            call addpnt(x2, y2, kdata, n2, 1.e+38, y2(n2))
            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R01'!, jlabel(j)
                stop
            end if
!       read parameters from Michelsen, H. A., R.J. Salawitch, P. O. Wennber,
!       and J. G. Anderson, Geophys. Res. Lett., 21, 2227-2230, 1994.
        case (kmich)
            open (UNIT=kin, FILE=trim(files(29)%fileName), STATUS='old')
            read (kin, *)
            read (kin, *)
            read (kin, *)
            n1 = 21
            n2 = n1
            do i = 1, n1
                read (kin, *) x1(i), y1(i), y2(i)
                x2(i) = x1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n1, 0., y1(1))
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), y1(n1))
            call addpnt(x1, y1, kdata, n1, 1.e+38, y1(n1))
            call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R01'!, jlabel(j)
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
            call addpnt(x2, y2, kdata, n2, 0., y2(1))
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), y2(n2))
            call addpnt(x2, y2, kdata, n2, 1.e+38, y2(n2))
            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R01'!, jlabel(j)
                stop
            end if
! quantum yield data from
! Shetter et al, J.Geophys.Res., v 101 (D9), pg. 14,631-14,641, June 20, 1996
        case (kshet)
            open (UNIT=kin, FILE=trim(files(30)%fileName), STATUS='OLD')
            read (kin, *) idum, n
            do i = 1, idum - 2
                read (kin, *)
            end do
            n = n - 2
            do i = 1, n
                read (kin, *) x1(i), y3(i), y4(i), y1(i), y2(i)
                x2(i) = x1(i)
                x3(i) = x1(i)
                x4(i) = x1(i)
            end do
            do i = n + 1, n + 2
                read (kin, *) x3(i), y3(i), y4(i)
                x4(i) = x3(i)
            end do
            close (kin)

            n1 = n
            n2 = n
            n3 = n + 2
            n4 = n + 2

! coefficients for exponential fit:

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n1, 0., y1(1))
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1e38, 0.)

            call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R01'!, jlabel(j)
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
            call addpnt(x2, y2, kdata, n2, 0., y2(1))
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1e38, 0.)

            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R01'!, jlabel(j)
                stop
            end if

! phi data at 298 and 230 K

            call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), y3(1))
            call addpnt(x3, y3, kdata, n3, 0., y3(1))
            call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 1e38, 0.)

            call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R01'!,jlabel(j)
                stop
            end if

            call addpnt(x4, y4, kdata, n4, x4(1)*(1.-deltax), y4(1))
            call addpnt(x4, y4, kdata, n4, 0., y4(1))
            call addpnt(x4, y4, kdata, n4, x4(n4)*(1.+deltax), 0.)
            call addpnt(x4, y4, kdata, n4, 1e38, 0.)

            call inter2(nw, wl, yg4(:, nReact), n4, x4, y4, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R01'!,jlabel(j)
                stop
            end if
        end select

!================================================================
!R02
!================================================================
        nReact = 02
        kData = 200
        call ChangeKData(kData)
! cross section
!------------NEED TO CHANGE kdata = 1000 FOR DAVIDSON ET AL. DATA---------
! measurements by:
! Davidson, J. A., C. A. Cantrell, A. H. McDaniel, R. E. Shetter,
! S. Madronich, and J. G. Calvert, Visible-ultraviolet absorption
! cross sections for NO2 as a function of temperature, J. Geophys.
! Res., 93, 7105-7112, 1988.
!     from 263.8 to 648.8 nm in approximately 0.5 nm intervals
!     OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_ncar_00.abs',STATUS='old')
!     n = 750
!     DO i = 1, n
!        READ(kin,*) x1(i), y1(i), dum, dum, idum
!     ENDDO
!     CLOSE(kin)

!     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,               0.,0.)
!     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!     CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
!     CALL inter2(nw,wl,yg,n,x1,y1,ierr)
!     IF (ierr .NE. 0) THEN
!        WRITE(*,*) ierr, jlabel(j)
!        STOP
!     ENDIF

! cross section data from JPL 94 recommendation
! JPL 97 recommendation is identical
        select case (mOption(3))
        case (1)
            open (UNIT=kin, FILE=trim(files(31)%fileName), STATUS='old')
            read (kin, *) idum, n
            do i = 1, idum - 2
                read (kin, *)
            end do
! read in wavelength bins, cross section at T0 and temperature correction
! coefficient a;  see input file for details.
! data need to be scaled to total area per bin so that they can be used with
! inter3

            do i = 1, n
                read (kin, *) x1(i), x3(i), y1(i), dum, y2(i)
                y1(i) = (x3(i) - x1(i))*y1(i)*1.e-20
                y2(i) = (x3(i) - x1(i))*y2(i)*1.e-22
                x2(i) = x1(i)
            end do
            close (kin)

            x1(n + 1) = x3(n)
            x2(n + 1) = x3(n)
            n = n + 1
            n1 = n

            call inter3(nw, wl, yg1(:, nReact), n, x1, y1, 0)
            call inter3(nw, wl, yg2(:, nReact), n1, x2, y2, 0)

! yg1(:,nReact), yg2(:,nReact) are per nm, so rescale by bin widths

            do iw = 1, nw - 1
                yg1(iw, nReact) = yg1(iw, nReact)/(wl(iw + 1) - wl(iw))
                yg2(iw, nReact) = yg2(iw, nReact)/(wl(iw + 1) - wl(iw))
            end do

        case (2)
            open (UNIT=kin, FILE=trim(files(32)%fileName), STATUS='old')
            do i = 1, 9
                read (kin, *)
            end do
            n = 135
            do i = 1, n
                read (kin, *) idum, y1(i)
                x1(i) = FLOAT(idum)
            end do

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n, 0., y1(1))
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n, x1, y1, ierr)
        end select

! quantum yield
! from Gardiner, Sperry, and Calvert
        open (UNIT=kin, FILE=trim(files(33)%fileName), STATUS='old')
        do i = 1, 8
            read (kin, *)
        end do
        n = 66
        do i = 1, n
            read (kin, *) x1(i), y1(i)
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
        call addpnt(x1, y1, kdata, n, 0., y1(1))
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg1n(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R02'!, jlabel(j)
            stop
        end if

!================================================================
!R03
!================================================================
        nReact = 03
        kData = 350
        call ChangeKData(kData)

! cross section
!     measurements of Graham and Johnston 1978
        open (UNIT=kin, FILE=trim(files(34)%fileName), STATUS='old')
        do i = 1, 9
            read (kin, *)
        end do
        n = 305
        do irow = 1, 30
            read (kin, *) (y1(10*(irow - 1) + icol), icol=1, 10)
        end do
        read (kin, *) (y1(300 + icol), icol=1, 5)
        close (kin)
        do i = 1, n
            y1(i) = y1(i)*1.e-19
            x1(i) = 400.+1.*FLOAT(i - 1)
        end do

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R03'!, jlabel(j)
            stop
        end if

!     cross section from JPL94:
        open (UNIT=kin, FILE=trim(files(35)%fileName), STATUS='old')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)
        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg1(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R03'!, jlabel(j)
            stop
        end if

! use JPL94 for wavelengths longer than 600 nm
        do iw = 1, nw - 1
            if (wl(iw) > 600.) yg(iw, nReact) = yg1(iw, nReact)
        end do

!================================================================
!R04
!================================================================
        nReact = 04
        kData = 100
        call ChangeKData(kData)
! cross section from jpl97, table up to 280 nm

        open (UNIT=kin, FILE=trim(files(36)%fileName), STATUS='old')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1.e-20
        end do
        xs270 = y1(n - 2)
        xs280 = y1(n)

        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), y1(n))
        call addpnt(x1, y1, kdata, n, 1.e36, y1(n))

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (0, *) ierr, 'Reading R04'!,jlabel(j)
            stop
        end if

!================================================================
!R05
!================================================================
        nReact = 05
        kData = 100
        call ChangeKData(kData)
        open (UNIT=kin, FILE=trim(files(37)%fileName), STATUS='old')
        do i = 1, 13
            read (kin, *)
        end do
        n = 91
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1.e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R05'!, jlabel(j)
            stop
        end if

!================================================================
!R06
!================================================================
        nReact = 06
        kData = 100
        call ChangeKData(kData)
!* cross section from JPL85

!      OPEN(UNIT=kin,FILE='dataj1/abs/HNO3.abs',STATUS='old')
!      DO i = 1, 9
!         READ(kin,*)
!      ENDDO
!      n = 29
!      DO i = 1, n
!         READ(kin,*) x1(i), y1(i)
!         y1(i) = y1(i) * 1.E-20
!      ENDDO
!      CLOSE (kin)

!      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,               0.,0.)
!      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
!      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
!      IF (ierr .NE. 0) THEN
!         WRITE(*,*) ierr, jlabel(j)
!         STOP
!      ENDIF

!* quantum yield = 1

!      qy = 1.
!      DO iw = 1, nw-1
!         DO i = 1, nz
!            sq(j,i,iw) = yg(iw)*qy
!         ENDDO
!      ENDDO

! HNO3 cross section parameters from Burkholder et al. 1993

        open (UNIT=kin, FILE=trim(files(38)%fileName), STATUS='old')
        do i = 1, 6
            read (kin, *)
        end do
        n1 = 83
        n2 = n1
        do i = 1, n1
            read (kin, *) y1(i), y2(i)
            x1(i) = 184.+i*2.
            x2(i) = x1(i)
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R06'!, jlabel(j)
            stop
        end if

        call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
        call addpnt(x2, y2, kdata, n2, 0., y2(1))
        call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), y2(n2))
        call addpnt(x2, y2, kdata, n2, 1.e+38, y2(n2))
        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R06'!, jlabel(j)
            stop
        end if

!================================================================
!R07
!================================================================
        nReact = 07
        kData = 100
        call ChangeKData(kData)
        open (UNIT=kin, FILE=trim(files(39)%fileName), STATUS='old')
        do i = 1, 4
            read (kin, *)
        end do
        n = 31
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1.e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R07'!, jlabel(j)
            stop
        end if

!================================================================
!R08
!================================================================
        nReact = 08
        kData = 100
        call ChangeKData(kData)
!     OPEN(UNIT=kin,FILE='dataj1/abs/H2O2_lin.abs',STATUS='old')
!     DO i = 1, 7
!        READ(kin,*)
!     ENDDO
!     n = 32
!     DO i = 1, n
!        READ(kin,*) x1(i), y1(i)
!        y1(i) = y1(i) * 1.E-20
!     ENDDO
!     CLOSE (kin)

!      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,               0.,0.)
!      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
!      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
!      IF (ierr .NE. 0) THEN
!         WRITE(*,*) ierr, jlabel(j)
!         STOP
!      ENDIF

! cross section from JPL94 (identical to JPL97)
! tabulated data up to 260 nm
        open (UNIT=kin, FILE=trim(files(40)%filename), STATUS='old')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1.e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R08'!, jlabel(j)
            stop
        end if

!================================================================
!R09
!================================================================
        nReact = 09
        kData = 100
        call ChangeKData(kData)
        select case (mOption(4))
        case (1)
            open (UNIT=kin, FILE=trim(files(41)%fileName), STATUS='old')
            do i = 1, 5
                read (kin, *)
            end do

            n5 = 25
            n4 = 27
            n3 = 29
            n2 = 31
            n1 = 39
            do i = 1, n5
                read (kin, *) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
            end do
            do i = n5 + 1, n4
                read (kin, *) x1(i), y1(i), y2(i), y3(i), y4(i)
            end do
            do i = n4 + 1, n3
                read (kin, *) x1(i), y1(i), y2(i), y3(i)
            end do
            do i = n3 + 1, n2
                read (kin, *) x1(i), y1(i), y2(i)
            end do
            do i = n2 + 1, n1
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            do i = 1, n1
                y1(i) = y1(i)*1.e-23
            end do
            do i = 1, n2
                x2(i) = x1(i)
                y2(i) = y2(i)*1.e-23
            end do
            do i = 1, n3
                x3(i) = x1(i)
                y3(i) = y3(i)*1.e-23
            end do
            do i = 1, n4
                x4(i) = x1(i)
                y4(i) = y4(i)*1.e-23
            end do
            do i = 1, n5
                x5(i) = x1(i)
                y5(i) = y5(i)*1.e-23
            end do

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n1, 0., y1(1))
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R09'!, jlabel(j)
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
            call addpnt(x2, y2, kdata, n2, 0., y2(1))
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R09'!, jlabel(j)
                stop
            end if

            call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), y3(1))
            call addpnt(x3, y3, kdata, n3, 0., y3(1))
            call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 1.e+38, 0.)
            call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R09'!, jlabel(j)

            end if

            call addpnt(x4, y4, kdata, n4, x4(1)*(1.-deltax), y4(1))
            call addpnt(x4, y4, kdata, n4, 0., y4(1))
            call addpnt(x4, y4, kdata, n4, x4(n4)*(1.+deltax), 0.)
            call addpnt(x4, y4, kdata, n4, 1.e+38, 0.)
            call inter2(nw, wl, yg4(:, nReact), n4, x4, y4, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R09'!, jlabel(j)
                stop
            end if

            call addpnt(x5, y5, kdata, n5, x5(1)*(1.-deltax), y5(1))
            call addpnt(x5, y5, kdata, n5, 0., y5(1))
            call addpnt(x5, y5, kdata, n5, x5(n5)*(1.+deltax), 0.)
            call addpnt(x5, y5, kdata, n5, 1.e+38, 0.)
            call inter2(nw, wl, yg5(:, nReact), n5, x5, y5, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R09'!, jlabel(j)
                stop
            end if
! jpl97, with temperature dependence formula,
!w = 290 nm to 340 nm,
!T = 210K to 300 K
!sigma, cm2 = exp((0.06183-0.000241*w)*(273.-T)-(2.376+0.14757*w))
        case (2)
            open (UNIT=kin, FILE=trim(files(42)%fileName), STATUS='old')
            do i = 1, 6
                read (kin, *)
            end do
            n1 = 87
            do i = 1, n1
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n1, 0., y1(1))
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R09'!, jlabel(j)
                stop
            end if
        end select

!================================================================
!R10
!================================================================
        nReact = 10
        kData = 16000
        call ChangeKData(kData)
        select case (mOption(5))
        case (1)
! read NBS/Bass data
            open (UNIT=kin, FILE=trim(files(43)%fileName), STATUS='old')
            n = 4032
            do i = 1, n
                read (kin, *) x(i), y(i)
            end do
            call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
            call addpnt(x, y, kdata, n, 0., 0.)
            call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
            call addpnt(x, y, kdata, n, 1.e+38, 0.)

            call inter2(nw, wl, yg1(:, nReact), n, x, y, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 10'
                stop
            end if
        case (2:4)
            open (UNIT=kin, FILE=trim(files(44)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 121
            do i = 1, n
                read (kin, *) x(i), y(i)
                y(i) = y(i)*1.e-20
            end do
            close (kin)
            call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
            call addpnt(x, y, kdata, n, 0., 0.)
            call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
            call addpnt(x, y, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n, x, y, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 10'
                stop
            end if
        end select
        select case (mOption(5))
        case (3)
! data are on wavenumber grid (cm-1), so convert to wavelength in nm:
! grid was on increasing wavenumbers, so need to reverse to get increasing
! wavelengths
! cross section assumed to be zero for wavelengths longer than 360 nm
! if y1 < 0, then make = 0 (some negative cross sections, actually 273 K intercepts
! are in the original data,  Here, make equal to zero)
            open (kin, FILE=trim(files(45)%fileName), STATUS='old')
            read (kin, *) idum, n
            do i = 1, idum - 2
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i), y2(i)
                x1(i) = 1./x1(i)*1e7
                if (x1(i) > 360.) then
                    y1(i) = 0.
                    y2(i) = 0.
                end if
            end do
            close (kin)

            do i = 1, n/2
                irev = n + 1 - i
                dum = x1(i)
                x1(i) = x1(irev)
                x1(irev) = dum
                dum = y1(i)
                y1(i) = y1(irev)
                y1(irev) = dum
                dum = y2(i)
                y2(i) = y2(irev)
                y2(irev) = dum
            end do
            do i = 1, n
                x2(i) = x1(i)
                y1(i) = max(y1(i), 0.)
            end do
            n1 = n
            n2 = n

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1e38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 0., 0.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1e38, 0.)
            call inter2(nw, wl, yg3(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if
        case (4)
            open (UNIT=kin, FILE=trim(files(46)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 23
            do i = 1, n
                read (kin, *) x2(i), y2(i), y3(i), dum, dum
                x3(i) = x2(i)
            end do
            close (kin)
            n2 = n
            n3 = n

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 0., 0.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1e38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if

            call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 0., 0.)
            call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 1e38, 0.)
            call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if
        case (5)
! read Rodgers data
            open (UNIT=kin, FILE=trim(files(47)%fileName), STATUS='old')
            do i = 1, 10
                read (kin, *)
            end do
            n = 261
            do i = 1, n
                read (kin, *) x(i), y(i), dum
                y(i) = y(i)*1.e-20
            end do
            call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
            call addpnt(x, y, kdata, n, 0., 0.)
            call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
            call addpnt(x, y, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n, x, y, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 10'
                stop
            end if
        case (6)
            open (UNIT=kin, FILE=trim(files(48)%fileName), STATUS='old')
            do i = 1, 3
                read (kin, *)
            end do
            n = 126
            do i = 1, n
                read (kin, *) x2(i), y2(i), y3(i)
                x3(i) = x2(i)
            end do
            close (kin)
            n2 = n
            n3 = n

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 0., 0.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1e38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if

            call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 0., 0.)
            call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 1e38, 0.)
            call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if
        end select
! quantum yield
        select case (mOption(6))
        case (1)
            open (UNIT=kin, FILE=trim(files(49)%fileName), STATUS='old')
            do i = 1, 11
                read (kin, *)
            end do
            n = 20
            do i = 1, n
                read (kin, *) x(i), y(i)
            end do
            close (kin)
            call addpnt(x, y, kdata, n, x(1)*(1.-deltax), y(1))
            call addpnt(x, y, kdata, n, 0., y(1))
            call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
            call addpnt(x, y, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg4(:, nReact), n, x, y, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 10'
                stop
            end if
            open (UNIT=kin, FILE=trim(files(50)%fileName), STATUS='old')
            do i = 1, 9
                read (kin, *)
            end do
            n = 33
            do i = 1, n
                read (kin, *) x(i), y(i)
            end do
            close (kin)
            call addpnt(x, y, kdata, n, x(1)*(1.-deltax), y(1))
            call addpnt(x, y, kdata, n, 0., y(1))
            call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
            call addpnt(x, y, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg5(:, nReact), n, x, y, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if
        case (2)
            open (UNIT=kin, FILE=trim(files(51)%fileName), STATUS='old')
            do i = 1, 7
                read (kin, *)
            end do
            n = 13
            do i = 1, n
                read (kin, *) x1(i), y1(i), y2(i)
                x2(i) = x1(i)
            end do
            close (kin)
            n1 = n
            n2 = n

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n1, 0., y1(1))
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg4(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
            call addpnt(x2, y2, kdata, n2, 0., y2(1))
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yg5(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if
        case (3)
            open (UNIT=kin, FILE=trim(files(52)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 23
            do i = 1, n
                read (kin, *) x1(i), dum, dum, dum, dum, y1(i), y2(i)
                x2(i) = x1(i)
            end do
            close (kin)
            n1 = n
            n2 = n

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n1, 0., y1(1))
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg4(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
            call addpnt(x2, y2, kdata, n2, 0., y2(1))
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yg5(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R10'
                stop
            end if
        end select

!================================================================
!R11
!================================================================
        nReact = 11
        kData = 150
        call ChangeKData(kData)
        select case (mOption(7))
        case (1)
            open (UNIT=kin, FILE=trim(files(53)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 106
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if
        case (2)
! cross section from Calvert and  Pitts
            open (UNIT=kin, FILE=trim(files(54)%fileName), STATUS='old')
            do i = 1, 14
                read (kin, *)
            end do
            n = 54
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                x1(i) = x1(i)/10.
                y1(i) = y1(i)*3.82e-21
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if
        case (3)
            open (UNIT=kin, FILE=trim(files(55)%fileName), STATUS='old')
            do i = 1, 3
                read (kin, *)
            end do
            n = 106
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if
        case (4)
!cross section from KFA tables
!ch3cho.001 - Calvert and Pitts 1966
!ch3cho.002 - Meyrahn thesis 1984
!ch3cho.003 - Schneider and Moortgat, priv comm. MPI Mainz 1989, 0.012 nm resol.
!ch3cho.004 - Schneider and Moortgat, priv comm. MPI Mainz 1989, 0.08  nm resol.
!ch3cho.005 - IUPAC'92
!ch3cho.006 - Libuda, thesis Wuppertal 1992

! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.001',STATUS='old')
! n = 217
! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.002',STATUS='old')
! n = 63
! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.003',STATUS='old')
! n = 13738
! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.004',STATUS='old')
! n = 2053
            open (UNIT=kin, FILE=trim(files(56)%fileName), STATUS='old')
            n = 18
! OPEN(UNIT=kin,FILE='dataj2/kfa/ch3cho.006',STATUS='old')
! n = 1705

            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if
        end select
! quantum yields
        select case (mOption(8))
        case (1)
            open (UNIT=kin, FILE=trim(files(57)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 12
            do i = 1, n
                read (kin, *) x1(i), y2(i), y1(i)
                x2(i) = x1(i)
            end do
            close (kin)
            n1 = n
            n2 = n

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 0., 0.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if

            do iw = 1, nw - 1
                yg3(:, nReact) = 0.
            end do
        case (2)
            open (UNIT=kin, FILE=trim(files(58)%fileName), STATUS='old')
            do i = 1, 18
                read (kin, *)
            end do
            n = 10
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n, 0., y1(1))
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if

            open (UNIT=kin, FILE=trim(files(59)%fileName), STATUS='old')
            do i = 1, 10
                read (kin, *)
            end do
            n = 9
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n, 0., y1(1))
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if

            open (UNIT=kin, FILE=trim(files(60)%fileName), STATUS='old')
            do i = 1, 10
                read (kin, *)
            end do
            n = 9
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), y1(1))
            call addpnt(x1, y1, kdata, n, 0., y1(1))
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg3(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R11'
                stop
            end if
        end select
!pressure-dependence parameters
        open (UNIT=kin, FILE=trim(files(61)%fileName), STATUS='old')
        do i = 1, 4
            read (kin, *)
        end do
        n = 5
        do i = 1, n
            read (kin, *) x1(i), dum, dum, y1(i)
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg4(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R11'
            stop
        end if

!================================================================
!R12
!================================================================
        nReact = 12
        kData = 150
        call ChangeKData(kData)

        select case (mOption(9))
        case (1)
            open (UNIT=kin, FILE=trim(files(62)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 106
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 12'
                stop
            end if
        case (2)
! cross section from KFA tables
! c2h5cho.001 - Calvert and Pitts 1966

            open (UNIT=kin, FILE=trim(files(63)%fileName), STATUS='old')
            n = 83

            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 12'
                stop
            end if
        end select
! quantum yields
        select case (mOption(10))
        case (1)
!PRINT *,'LFR->Abrindo "'//trim(files(64)%fileName)//'"'
            open (UNIT=kin, FILE=trim(files(64)%fileName), STATUS='old')
!PRINT *,'LFR->Abriu "'//trim(files(64)%fileName)//'"'
            do i = 1, 4
                read (kin, *)
            end do
            n = 5
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)
            n1 = n

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, 340., 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 12'
                stop
            end if
        case (2)
            stop
        end select

!================================================================
!R13
!================================================================
        nReact = 13
        kData = 500
        call ChangeKData(kData)

        select case (mOption(11))
        case (1)
            open (UNIT=kin, FILE=trim(files(65)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 110
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R13'
                stop
            end if
        case (2)
! cross section from KFA tables
! chocho.001 - Plum et al. 1983
            open (UNIT=kin, FILE=trim(files(66)%fileName), STATUS='old')
            n = 219

            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R13'
                stop
            end if
        case (3)
! cross section from Orlando et la.
! Orlando, J. J.; G. S. Tyndall, 2001:  The atmospheric chemistry of the
! HC(O)CO radical. Int. J. Chem. Kinet., 33, 149-156.
            open (UNIT=kin, FILE=trim(files(67)%fileName), STATUS='old')

            do i = 1, 6
                read (kin, *)
            end do
            n = 481
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R13'
                stop
            end if
        case (4)
            open (UNIT=kin, FILE=trim(files(68)%fileName), STATUS='old')

            do i = 1, 8
                read (kin, *)
            end do
            n = 270
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R13'
                stop
            end if
        end select

!================================================================
!R14
!================================================================
        nReact = 14
        kData = 500
        call ChangeKData(kData)

        select case (mOption(13))
        case (1)
            open (UNIT=kin, FILE=trim(files(69)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 38
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 14'
                stop
            end if

            open (UNIT=kin, FILE=trim(files(70)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 75
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 14'
                stop
            end if

            do iw = 1, nw - 1
                if (wc(iw) < 402.) then
                    yg(iw, nReact) = yg1(iw, nReact)
                else
                    yg(iw, nReact) = yg2(iw, nReact)
                end if
            end do
        case (2)
            open (UNIT=kin, FILE=trim(files(71)%fileName), STATUS='old')
            n = 271
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, 14), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 14'
                stop
            end if
        case (3:6)
!       cross section from KFA tables
!       ch3cocho.001 - Plum et al. 1983
!       ch3cocho.002 - Meller et al. 1991, 0.033 nm resolution
!       ch3cocho.003 - Meller et al. 1991, 1.0   nm resolution
!       ch3cocho.004 - Staffelbach et al. 1995
            select case (mOption(13))
            case (3)
                open (UNIT=kin, FILE=trim(files(72)%fileName), STATUS='old')
                n = 136
            case (4)
                open (UNIT=kin, FILE=trim(files(73)%fileName), STATUS='old')
                n = 8251
            case (5)
                open (UNIT=kin, FILE=trim(files(74)%fileName), STATUS='old')
                n = 275
            case (6)
                open (UNIT=kin, FILE=trim(files(75)%fileName), STATUS='old')
                n = 162
            end select

            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, 14), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 14'
                stop
            end if
        case (7)
            open (UNIT=kin, FILE=trim(files(76)%fileName), STATUS='old')
            do i = 1, 7
                read (kin, *)
            end do
            n = 55
            do i = 1, n
                read (kin, *) x(i), y(i)
                y(i) = y(i)*1.e-20
            end do
            close (kin)

            call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
            call addpnt(x, y, kdata, n, 0., 0.)
            call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
            call addpnt(x, y, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n, x, y, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 14'
                stop
            end if

            open (UNIT=kin, FILE=trim(files(77)%fileName), STATUS='old')
            do i = 1, 6
                read (kin, *)
            end do
            n = 481
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 14'
                stop
            end if

            do iw = 1, nw - 1
                yg(iw, nReact) = 0.5*(yg1(iw, nReact) + yg2(iw, nReact))
            end do
        end select
! quantum yields
        select case (mOption(14))
        case (4)
            open (UNIT=kin, FILE=trim(files(78)%fileName), STATUS='old')
            do i = 1, 5
                read (kin, *)
            end do
            n = 5
            do i = 1, n
                read (kin, *) x1(i), y1(i), y2(i)
                x2(i) = x1(i)
            end do
            close (kin)
            n1 = n
            n2 = n

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 1.)
            call addpnt(x1, y1, kdata, n1, 0., 1.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 14'
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 1.)
            call addpnt(x2, y2, kdata, n2, 0., 1.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 14'
                stop
            end if
        end select

!================================================================
!R15
!================================================================
        nReact = 15
        kData = 150
        call ChangeKData(kData)

        select case (mOption(15))
        case (1)
            open (UNIT=kin, FILE=trim(files(79)%fileName), STATUS='old')
            do i = 1, 6
                read (kin, *)
            end do
            n = 35
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*3.82e-21
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R15'
                stop
            end if
        case (2)
            open (UNIT=kin, FILE=trim(files(80)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 96
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R15'
                stop
            end if
        case (3)
            open (UNIT=kin, FILE=trim(files(81)%fileName), STATUS='old')
            do i = 1, 12
                read (kin, *)
            end do
            n = 135
            do i = 1, n
                read (kin, *) x1(i), y1(i), y2(i), y3(i)
                x2(i) = x1(i)
                x3(i) = x1(i)
            end do
            close (kin)
            n1 = n
            n2 = n
            n3 = n

            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R15'
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 0., 0.)
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
            call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
            call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R15'
                stop
            end if

            call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 0., 0.)
            call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
            call addpnt(x3, y3, kdata, n3, 1.e+38, 0.)
            call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R15'
                stop
            end if
        end select
        select case (mOption(16))
        case (2)
            open (UNIT=kin, FILE=trim(files(82)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 9
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg1(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R15'
                stop
            end if
        end select

!================================================================
!R16
!================================================================
        nReact = 16
        kData = 100
        call ChangeKData(kData)

        select case (mOption(17))
        case (1)
!            OPEN(UNIT=kin,FILE='dataj1/CH3OOH/CH3OOH_jpl85.abs',
!        $         STATUS='old')
!            OPEN(UNIT=kin,FILE='dataj1/CH3OOH/CH3OOH_jpl92.abs',
!        $         STATUS='old')
!            OPEN(UNIT=kin,FILE='dataj1/CH3OOH/CH3OOH_jpl94.abs',
!        $         STATUS='old')
            open (UNIT=kin, FILE=trim(files(83)%fileName), STATUS='old')
            read (kin, *) idum, n
            do i = 1, idum - 2
                read (kin, *)
            end do
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 16'
                stop
            end if
        case (2)
            open (UNIT=kin, FILE=trim(files(84)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 32
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 16'
                stop
            end if
        case (3)
            open (UNIT=kin, FILE=trim(files(85)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 12
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 16'
                stop
            end if
        case (4)
            open (UNIT=kin, FILE=trim(files(86)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 15
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n, 0., 0.)
            call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading 16'
                stop
            end if
        end select

!================================================================
!R17
!================================================================
        nReact = 17
        kData = 2000
        call ChangeKData(kData)

        select case (mOption(18))
        case (1)
            open (UNIT=kin, FILE=trim(files(87)%fileName), STATUS='old')
            do i = 1, 3
                read (kin, *)
            end do
            n = 15
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            n1 = n
            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if
        case (2)
!          sigma(T,lambda) = sigma(298,lambda) * exp(B * (T-298))
            open (UNIT=kin, FILE=trim(files(88)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 55
            do i = 1, n
                read (kin, *) x1(i), y1(i), y2(i)
                x2(i) = x1(i)
                y1(i) = y1(i)*1.e-20
            end do
            close (kin)

            n1 = n
            n2 = n
            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
            call addpnt(x2, y2, kdata, n2, 0., y2(1))
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), y2(n2))
            call addpnt(x2, y2, kdata, n2, 1.e+38, y2(n2))
            call inter2(nw, wl, yg1(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if
        case (3)
            open (UNIT=kin, FILE=trim(files(89)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 13
            do i = 1, n
                read (kin, *) x1(i), y1(i)
                y1(i) = y1(i)*1e-20
            end do
            close (kin)

            n1 = n
            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if
        case (4)
!          sigma(T,lambda) = sigma(298,lambda) * 10**(B * T)
            open (UNIT=kin, FILE=trim(files(90)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 7
            do i = 1, n
                read (kin, *) x1(i), y1(i), y2(i)
                x2(i) = x1(i)
                y1(i) = y1(i)*1.e-21
                y2(i) = y2(i)*1.e-3
            end do
            close (kin)

            n1 = n
            n2 = n
            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), -36.)
            call addpnt(x1, y1, kdata, n1, 0., -36.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), -36.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, -36.)
            call inter2(nw, wl, yg(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if

            call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), y2(1))
            call addpnt(x2, y2, kdata, n2, 0., y2(1))
            call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), y2(n2))
            call addpnt(x2, y2, kdata, n2, 1.e+38, y2(n2))
            call inter2(nw, wl, yg1(:, nReact), n2, x2, y2, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if
        case (5)
            open (UNIT=kin, FILE=trim(files(91)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 13
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            n1 = n
            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if
        case (6)
            do iw = 1, nw - 1
                if (wc(iw) > 284.) then
                    yg(iw, nReact) = exp(-1.044e-3*wc(iw)*wc(iw) + 0.5309*wc(iw) - 112.4)
                else
                    yg(iw, nReact) = 0.
                end if
            end do
        case (7)
            open (UNIT=kin, FILE=trim(files(92)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 24
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            n1 = n
            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if
        case (8)
            open (UNIT=kin, FILE=trim(files(93)%fileName), STATUS='old')
            do i = 1, 4
                read (kin, *)
            end do
            n = 1638
            do i = 1, n
                read (kin, *) x1(i), y1(i)
            end do
            close (kin)

            n1 = n
            call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 0., 0.)
            call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
            call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n1, x1, y1, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R17'
                stop
            end if
        end select

!================================================================
!R18
!================================================================
        nReact = 18
        kData = 100
        call ChangeKData(kData)

! cross section from Senum et al., 1984, J.Phys.Chem. 88/7, 1269-1270

!     OPEN(UNIT=kin,FILE='dataj1/RONO2/PAN_senum.abs',STATUS='OLD')
!     DO i = 1, 14
!        READ(kin,*)
!     ENDDO
!     n = 21
!     DO i = 1, n
!        READ(kin,*) x1(i), y1(i)
!        y1(i) = y1(i) * 1.E-20
!     ENDDO
!     CLOSE(kin)

!      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,               0.,0.)
!      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
!      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
!      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
!      IF (ierr .NE. 0) THEN
!         WRITE(*,*) ierr, 'Reading 18'
!         STOP
!      ENDIF

! cross section from
!      Talukdar et al., 1995, J.Geophys.Res. 100/D7, 14163-14174
        open (UNIT=kin, FILE=trim(files(94)%fileName), STATUS='OLD')
        do i = 1, 14
            read (kin, *)
        end do
        n = 78
        do i = 1, n
            read (kin, *) x1(i), y1(i), y2(i)
            y1(i) = y1(i)*1.e-20
            y2(i) = y2(i)*1e-3
            x2(i) = x1(i)
        end do
        n2 = n
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 18'
            stop
        end if

        call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 0., 0.)
        call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 1.e+38, 0.)
        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 18'
            stop
        end if

!================================================================
!R19
!================================================================
        nReact = 19
        kData = 100
        call ChangeKData(kData)

!** cross sections from JPL94 recommendation
        open (kin, FILE=trim(files(95)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R19'
            stop
        end if

!================================================================
!R20
!================================================================
        nReact = 20
        kData = 100
        call ChangeKData(kData)

!** cross sections from JPL97 recommendation (identical to 94 data)

        open (kin, FILE=trim(files(96)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R20'
            stop
        end if

!================================================================
!R21
!================================================================
        nReact = 21
        kData = 100
        call ChangeKData(kData)

!** cross sections from JPL97 recommendation (identical to 94 recommendation)

        open (kin, FILE=trim(files(97)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 21'
            stop
        end if

!================================================================
!R22
!================================================================
        nReact = 22
        kData = 100
        call ChangeKData(kData)

!*** cross sections from JPL97 recommendation (identical to 94 recommendation)
        open (kin, FILE=trim(files(98)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 22'
            stop
        end if

!================================================================
!R23
!================================================================
        nReact = 23
        kData = 100
        call ChangeKData(kData)

!** cross sections from JPL97 recommendation (identical to 94 recommendation)
        open (kin, FILE=trim(files(99)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i), y2(i)
            y1(i) = y1(i)*1e-20
            y2(i) = y2(i)*1e-20
            x2(i) = x1(i)
        end do
        close (kin)

        n1 = n
        n2 = n

!* sigma @ 295 K

        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1e38, 0.)

        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R23'
            stop
        end if

! sigma @ 210 K

        call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 0., 0.)
        call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 1e38, 0.)

        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R23'
            stop
        end if

!================================================================
!R24
!================================================================
        nReact = 24
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(100)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i), y2(i)
            y1(i) = y1(i)*1e-20
            y2(i) = y2(i)*1e-20
            x2(i) = x1(i)
        end do
        close (kin)

        n1 = n
        n2 = n

!* sigma @ 295 K

        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1e38, 0.)

        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R24'
            stop
        end if

! sigma @ 210 K

        call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 0., 0.)
        call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 1e38, 0.)

        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R24'
            stop
        end if

!================================================================
!R25
!================================================================
        nReact = 25
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(101)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 25'
            stop
        end if

!================================================================
!R26
!================================================================
        nReact = 26
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(102)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

!* sigma @ 298 K

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R26'
            stop
        end if

!================================================================
!R27
!================================================================
        nReact = 27
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(103)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

!* sigma @ 298 K

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R27'
            stop
        end if

!================================================================
!R28
!================================================================
        nReact = 28
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(104)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'REading R28'
            stop
        end if

!================================================================
!R29
!================================================================
        nReact = 29
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(105)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i), y2(i), y3(i)
            y1(i) = y1(i)*1e-20
            y2(i) = y2(i)*1e-20
            y3(i) = y3(i)*1e-20
            x2(i) = x1(i)
            x3(i) = x1(i)
        end do
        close (kin)

        n1 = n
        n2 = n
        n3 = n

!* sigma @ 295 K

        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1e38, 0.)

        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R29'
            stop
        end if

!* sigma @ 250 K

        call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 0., 0.)
        call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 1e38, 0.)

        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R29'
            stop
        end if

!* sigma @ 210 K

        call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
        call addpnt(x3, y3, kdata, n3, 0., 0.)
        call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
        call addpnt(x3, y3, kdata, n3, 1e38, 0.)

        call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R29'
            stop
        end if

!================================================================
!R30
!================================================================
        nReact = 30
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(106)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i), y2(i), y3(i)
            y1(i) = y1(i)*1e-20
            y2(i) = y2(i)*1e-20
            y3(i) = y3(i)*1e-20
            x2(i) = x1(i)
            x3(i) = x1(i)
        end do
        close (kin)

        n1 = n
        n2 = n
        n3 = n

!* sigma @ 296 K

        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1e38, 0.)

        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R30'
            stop
        end if

!* sigma @ 279 K

        call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 0., 0.)
        call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 1e38, 0.)

        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R30'
            stop
        end if

!* sigma @ 255 K

        call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
        call addpnt(x3, y3, kdata, n3, 0., 0.)
        call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
        call addpnt(x3, y3, kdata, n3, 1e38, 0.)

        call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R30'
            stop
        end if

!================================================================
!R31
!================================================================
        nReact = 31
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(107)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R31'
            stop
        end if

!================================================================
!R32
!================================================================
        nReact = 32
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(108)%fileName), STATUS='OLD')
        read (kin, *) idum
        do i = 1, idum - 2
            read (kin, *)
        end do
        read (kin, 100) inline
        read (inline(6:), *) tbar(nReact), i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        close (kin)

!================================================================
!R33
!================================================================
        nReact = 33
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(109)%fileName), STATUS='OLD')
        read (kin, *) idum
        idum = idum + 5
        do i = 1, idum - 2
            read (kin, *)
        end do
        read (kin, 100) inline
        read (inline(6:), *) tbar(nReact), i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        close (kin)

!================================================================
!R34
!================================================================
        nReact = 34
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(110)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'REading R34'
            stop
        end if

!================================================================
!R35
!================================================================
        nReact = 35
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(111)%fileName), STATUS='OLD')
        read (kin, *) idum
        idum = idum + 10
        do i = 1, idum - 2
            read (kin, *)
        end do
        read (kin, 101) inline
        read (inline(6:), *) tbar(nReact), i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        read (kin, *) i, (coeff(i, k, nReact), k=1, 3)
        close (kin)

!================================================================
!R36
!================================================================
        nReact = 36
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(112)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R36'
            stop
        end if

!================================================================
!R37
!================================================================
        nReact = 37
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(113)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R37'
            stop
        end if

!================================================================
!R38
!================================================================
        nReact = 38
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(114)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
            y1(i) = y1(i)*1e-20
            y2(i) = y2(i)*1e-20
            y3(i) = y3(i)*1e-20
            y4(i) = y4(i)*1e-20
            y5(i) = y5(i)*1e-20
            x2(i) = x1(i)
            x3(i) = x1(i)
            x4(i) = x1(i)
            x5(i) = x1(i)
        end do
        close (kin)

        n1 = n
        n2 = n
        n3 = n
        n4 = n
        n5 = n

!* sigma @ 295 K

        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1e38, 0.)

        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R38'
            stop
        end if

!* sigma @ 270 K

        call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 0., 0.)
        call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 1e38, 0.)

        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R38'
            stop
        end if

!* sigma @ 250 K

        call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
        call addpnt(x3, y3, kdata, n3, 0., 0.)
        call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
        call addpnt(x3, y3, kdata, n3, 1e38, 0.)

        call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R38'
            stop
        end if

!* sigma @ 230 K

        call addpnt(x4, y4, kdata, n4, x4(1)*(1.-deltax), 0.)
        call addpnt(x4, y4, kdata, n4, 0., 0.)
        call addpnt(x4, y4, kdata, n4, x4(n4)*(1.+deltax), 0.)
        call addpnt(x4, y4, kdata, n4, 1e38, 0.)

        call inter2(nw, wl, yg4(:, nReact), n4, x4, y4, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R38'
            stop
        end if

!* sigma @ 210 K

        call addpnt(x5, y5, kdata, n5, x5(1)*(1.-deltax), 0.)
        call addpnt(x5, y5, kdata, n5, 0., 0.)
        call addpnt(x5, y5, kdata, n5, x5(n5)*(1.+deltax), 0.)
        call addpnt(x5, y5, kdata, n5, 1e38, 0.)

        call inter2(nw, wl, yg5(:, nReact), n5, x5, y5, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R38'
            stop
        end if

!================================================================
!R39
!================================================================
        nReact = 39
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(115)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R39'
            stop
        end if

!================================================================
!R40
!================================================================
        nReact = 40
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(116)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R40'
            stop
        end if

!================================================================
!R41
!================================================================
        nReact = 41
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(117)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 41'
            stop
        end if

!================================================================
!R42
!================================================================
        nReact = 42
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(118)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R42'
            stop
        end if

!================================================================
!R43
!================================================================
        nReact = 43
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(119)%fileName), STATUS='OLD')
        read (kin, *) idum, n
        do i = 1, idum - 2
            read (kin, *)
        end do
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n) + deltax, 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)

        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)

        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 43'
            stop
        end if

!================================================================
!R44
!================================================================
! Fixed, not read.

!================================================================
!R45
!================================================================
        nReact = 45
        kData = 150
        call ChangeKData(kData)

        open (kin, FILE=trim(files(120)%fileName), STATUS='OLD')
        n = 119
        do i = 1, n
            read (kin, *) x1(i), y1(i), y2(i), y3(i)
            y1(i) = y1(i)*1e-20
            x2(i) = x1(i)
            x3(i) = x1(i)
        end do
        close (kin)

        n1 = n
        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1e38, 0.)
        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 45'
            stop
        end if

        n2 = n
        call addpnt(x2, y2, kdata, n2, x2(1)*(1.-deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 0., 0.)
        call addpnt(x2, y2, kdata, n2, x2(n2)*(1.+deltax), 0.)
        call addpnt(x2, y2, kdata, n2, 1e38, 0.)
        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 45'
            stop
        end if

        n3 = n
        call addpnt(x3, y3, kdata, n3, x3(1)*(1.-deltax), 0.)
        call addpnt(x3, y3, kdata, n3, 0., 0.)
        call addpnt(x3, y3, kdata, n3, x3(n3)*(1.+deltax), 0.)
        call addpnt(x3, y3, kdata, n3, 1e38, 0.)
        call inter2(nw, wl, yg3(:, nReact), n3, x3, y3, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 45'
            stop
        end if

!================================================================
!R46
!================================================================
        nReact = 46
        kData = 100
        call ChangeKData(kData)

        open (kin, FILE=trim(files(121)%fileName), STATUS='OLD')
        do i = 1, 13
            read (kin, *)
        end do
        n = 61
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        n1 = n
        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1e38, 0.)
        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 46'
            stop
        end if

!================================================================
!R47
!================================================================
        nReact = 47
        kData = 150
        call ChangeKData(kData)

        open (kin, FILE=trim(files(122)%fileName), STATUS='OLD')
        do i = 1, 5
            read (kin, *)
        end do
        n = 22
        do i = 1, n
            read (kin, *) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n, 0., 0.)
        call addpnt(x1, y1, kdata, n, x1(n)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n, 1e38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R47'
            stop
        end if

!================================================================
!R48 until R100
!================================================================
!Does not exist!

!================================================================
!R101
!================================================================
        nReact = 101
        kData = 300
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(123)%fileName), STATUS='old')
        do i = 1, 15
            read (kin, *)
        end do
        n = 131
        do i = 1, n
            read (kin, *) x(i), y(i)
        end do
        close (kin)

        call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
        call addpnt(x, y, kdata, n, 0., 0.)
        call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
        call addpnt(x, y, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R101'
            stop
        end if

!================================================================
!R102
!================================================================
        nReact = 102
        kData = 300
        call ChangeKData(kData)

        select case (mOption(19))
        case (1)
            open (UNIT=kin, FILE=trim(files(124)%fileName), STATUS='old')
            do i = 1, 7
                read (kin, *)
            end do
            n = 55
            do i = 1, n
                read (kin, *) x(i), y(i)
                y(i) = y(i)*1.e-20
            end do
            close (kin)

            call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
            call addpnt(x, y, kdata, n, 0., 0.)
            call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
            call addpnt(x, y, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R102'
                stop
            end if
        case (2)
            open (UNIT=kin, FILE=trim(files(125)%fileName), STATUS='old')
            do i = 1, 8
                read (kin, *)
            end do
            n = 287
            do i = 1, n
                read (kin, *) x(i), y(i)
                y(i) = y(i)*1.e-20
            end do
            close (kin)

            call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
            call addpnt(x, y, kdata, n, 0., 0.)
            call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
            call addpnt(x, y, kdata, n, 1.e+38, 0.)
            call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
            if (ierr /= 0) then
                write (*, *) ierr, 'Reading R102'
                stop
            end if
        end select

!================================================================
!R103
!================================================================
        nReact = 103
        kData = 20000
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(126)%fileName), STATUS='old')
        do i = 1, 9
            read (kin, *)
        end do
        n = 19682
        do i = 1, n
            read (kin, *) x(i), y(i)
        end do
        close (kin)

        call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
        call addpnt(x, y, kdata, n, 0., 0.)
        call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
        call addpnt(x, y, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R103'
            stop
        end if

!================================================================
!R104
!================================================================
        nReact = 104
        kData = 20000
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(127)%fileName), STATUS='old')
        do i = 1, 10
            read (kin, *)
        end do
        n = 15213
        do i = 1, n
            read (kin, *) x(i), y(i)
        end do
        close (kin)

        call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
        call addpnt(x, y, kdata, n, 0., 0.)
        call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
        call addpnt(x, y, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
!PRINT *,'React104 yg=',yg(:,nReact); CALL flush(6)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R104'
            stop
        end if

!================================================================
!R105
!================================================================
        nReact = 105
        kData = 20000
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(128)%fileName), STATUS='old')
        do i = 1, 8
            read (kin, *)
        end do
        n = 148
        do i = 1, n
            read (kin, *) x(i), y(i)
            y(i) = y(i)*1.e-20
        end do
        close (kin)

        call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
        call addpnt(x, y, kdata, n, 0., 0.)
        call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
        call addpnt(x, y, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R105'
            stop
        end if

!================================================================
!R106
!================================================================
        nReact = 106
        kData = 200
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(129)%fileName), STATUS='old')
        do i = 1, 10
            read (kin, *)
        end do
        n1 = 0
        n2 = 0
        do i = 1, 63
            read (kin, *) x1(i), dum, dum, y1(i), y2(i), dum, dum
            if (y1(i) > 0.) n1 = n1 + 1
            if (y2(i) > 0.) n2 = n2 + 1
            x2(i) = x1(i)
            y1(i) = y1(i)*1.e-20
            y2(i) = y2(i)*1.e-3
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R106'
            stop
        end if

        call addpnt(x2, y2, kdata, n2, 0., y2(1))
        call addpnt(x2, y2, kdata, n2, 1.e+38, y2(n2))
        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R106'
            stop
        end if

!================================================================
!R107
!================================================================
        nReact = 107
        kData = 200
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(130)%fileName), STATUS='old')
        do i = 1, 10
            read (kin, *)
        end do
        n1 = 0
        n2 = 0
        do i = 1, 63
            read (kin, *) x1(i), dum, dum, dum, dum, y1(i), y2(i)
            if (y1(i) > 0.) n1 = n1 + 1
            if (y2(i) > 0.) n2 = n2 + 1
            x2(i) = x1(i)
            y1(i) = y1(i)*1.e-20
            y2(i) = y2(i)*1.e-3
        end do
        close (kin)

        call addpnt(x1, y1, kdata, n1, x1(1)*(1.-deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 0., 0.)
        call addpnt(x1, y1, kdata, n1, x1(n1)*(1.+deltax), 0.)
        call addpnt(x1, y1, kdata, n1, 1.e+38, 0.)
        call inter2(nw, wl, yg1(:, nReact), n1, x1, y1, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R107'
            stop
        end if

        call addpnt(x2, y2, kdata, n2, 0., y2(1))
        call addpnt(x2, y2, kdata, n2, 1.e+38, y2(n2))
        call inter2(nw, wl, yg2(:, nReact), n2, x2, y2, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R107'
            stop
        end if

!================================================================
!R108 until R110
!================================================================
!Do not read - fixed

!================================================================
!R111
!================================================================
        nReact = 111
        kData = 20000
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(131)%fileName), STATUS='old')
        do i = 1, 25
            read (kin, *)
        end do
        n = 131
        do i = 1, n
            read (kin, *) x(i), y(i)
            y(i) = y(i)*1.e-20
        end do
        close (kin)

        call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
        call addpnt(x, y, kdata, n, 0., 0.)
        call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
        call addpnt(x, y, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R111'
            stop
        end if

!================================================================
!R112
!================================================================
        nReact = 112
        kData = 20000
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(132)%fileName), STATUS='old')
        do i = 1, 8
            read (kin, *)
        end do
        n = 101
        do i = 1, n
            read (kin, *) x(i), y(i)
        end do
        close (kin)

        call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
        call addpnt(x, y, kdata, n, 0., 0.)
        call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
        call addpnt(x, y, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading 112'
            stop
        end if

!================================================================
!R113
!================================================================
!Do not read - fixed

!================================================================
!R114
!================================================================
        nReact = 114
        kData = 20000
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(133)%fileName), STATUS='old')
        do i = 1, 14
            read (kin, *)
        end do
        n = 15
        do i = 1, n
            read (kin, *) x(i), dum, y(i)
            y(i) = y(i)*1.e-20
        end do
        n = n + 1
        x(n) = dum
        close (kin)

! use bin-to-bin interpolation

        call inter4(nw, wl, yg(:, nReact), n, x, y, 1)

!================================================================
!R115
!================================================================
        nReact = 115
        kData = 50
        call ChangeKData(kData)

        open (UNIT=kin, FILE=trim(files(134)%fileName), STATUS='old')

        do i = 1, 6
            read (kin, *)
        end do
        n = 29
        do i = 1, n
            read (kin, *) x(i), y(i)
        end do
        close (kin)

        call addpnt(x, y, kdata, n, x(1)*(1.-deltax), 0.)
        call addpnt(x, y, kdata, n, 0., 0.)
        call addpnt(x, y, kdata, n, x(n)*(1.+deltax), 0.)
        call addpnt(x, y, kdata, n, 1.e+38, 0.)
        call inter2(nw, wl, yg(:, nReact), n, x, y, ierr)
        if (ierr /= 0) then
            write (*, *) ierr, 'Reading R115'
            stop
        end if

100     format(a120)
101     format(a80)

    end subroutine ReadAll

    integer function reverse(pos, pmax)
        integer, intent(IN) :: pos, pmax
        reverse = pmax - pos + 1
    end function reverse

    subroutine alert(message)
        character(LEN=*), intent(IN) :: message

        write (noPr, *) message
        !call flush (noPr)

    end subroutine alert

    integer function get_nz()
        get_nz = current_nz
    end function

end module ModTuv