module ModTuv_driver

  use mpas_pool_routines
  use mpas_derived_types
  use mpas_constants, only : &
      rvord         ! huge(1.0_RKIND)
  use mpas_atmphys_constants,only: &
      degrad        ! conversion from degree to radiant
  use mpas_atmphys_manager, only: &
      gmt &         ! Greenwich mean time hour of model start (hr)
    , curr_julday & ! Current Julian day (= 0.0 at 0Z on January 1st)
    , julday &      ! Initial Julian day
    , year          ! Current year
  !
  use ModTuv, only : InitTuv, Tuv, nWavelengths => kw, wu, nw, ks, kj &
                    ,xRef
  use chem_list, only : nspecies, spc_name, nr_photo
  implicit none

!Include the parameters with arrays values
#include "tuvParam.inc"

  integer,dimension(:),allocatable :: tuv2carma
  integer ::  validCells
  real(kind=RKIND), dimension(33, 9, 6) :: mclat
  real(kind=RKIND), dimension(33, 6) :: mcol

  integer :: iposno2, iposso2, iposo3, maxNRad
  logical :: no2present, so2present, o3present

  private
  public :: Tuv_driver

  real(kind=RKIND), allocatable :: tlev(:,:)
  !! temperature at each specified altitude level (k, i) [K]
  real(kind=RKIND), allocatable :: tlay(:,:)
  !! temperature at each specified altitude layer (k, i) [K]
  real(kind=RKIND), allocatable :: zml(:,:)
  !! vector of altitude levels of middle (k, i) [m]
  real(kind=RKIND), allocatable :: ztl(:,:)
  !! vector of altitude levels of top (k, i) [m]
  real(kind=RKIND), allocatable :: prd(:,:)
  !! pressure at each altitude and cell (k, i) [Pa]
  real(kind=RKIND), allocatable :: temprd(:,:)
  !! temperature at each altitude and cell (k, i) [K]
  real(kind=RKIND), allocatable :: dair(:,:)
  !! Air density at each altitude and cell [molec/cm3]  
  real(kind=RKIND), allocatable :: dzl(:,:)
  !! thickness of each altitude layer (k, i) [km]
  real(kind=RKIND), allocatable :: rv(:,:)
  !! mixing ratio of water vapor at each altitude and cell [g/kg]
  real(kind=RKIND), allocatable :: o3l(:,:)
  !! Ozone concentration at each altitude and cell [molec/cm3]
  real(kind=RKIND),allocatable :: zLevel(:,:)
  !! vector of altitude levels (k, i) [km]
  real(kind=RKIND), allocatable :: albedo(:,:)
  !! Surface albedo at each cell and wavelength - Notice: using the same for all wavelengths for now, 
  !! but it can be changed to be wavelength dependent if needed
  real(kind=RKIND), allocatable :: dtcld(:,:,:) 
  !!  optical depth due to absorption by clouds at each altitude and wavelength
  real(kind=RKIND), allocatable :: omcld(:,:,:)   
  !! single scattering albedo due to clouds at each defined altitude and wavelength
  real(kind=RKIND), allocatable :: gcld (:,:,:)   
  !! cloud asymmetry factor at each defined altitude and wavelength
  real(kind=RKIND), allocatable :: dtaer(:,:,:)   
  !! optical depth due to absorption by aerosols at each altitude and wavelength
  real(kind=RKIND), allocatable :: omaer(:,:,:)  
  !!  single scattering albedo due to aerosols at each defined altitude and wavelength
  real(kind=RKIND), allocatable :: gaer (:,:,:)   
  !! aerosol asymmetry factor at each defined altitude and wavelength
  real(kind=RKIND), allocatable :: O3(:,:)
  !! ozone concentration at each altitude and cell [molec/cm3]
  real(kind=RKIND),allocatable :: cAir(:,:)
  !! number of air molecules per cm^2 at each altitude
  real(kind=RKIND), allocatable :: airden(:,:)
  !! AIr density at each altitude and cell [kg/m3]
  real(kind=RKIND), allocatable :: tCo3(:,:)
  !! Total of column O3 [cm-2] 
  real(kind=RKIND), allocatable :: so2col(:)    
  !! Total columns of SO2 (Dobson Units)
  real(kind=RKIND), allocatable :: no2col(:)   
  !! Total columns of NO2 (Dobson Units)
  real(kind=RKIND), allocatable :: zt(:)
  !! vector of altitude levels of top (k) [m]
  real(kind=RKIND), allocatable :: zm(:)
  !! vector of altitude levels of middle (k) [m]
  !real(kind=RKIND), allocatable :: rate(:,:,:)  
  !! Weighted irradiances (dose rates) W m-2
  real(kind=RKIND), allocatable :: valj(:,:,:)  
  !! Photolysis coefficients (j-values)
  !real(kind=RKIND), allocatable :: sirrad(:,:,:)
  !! Spectral irradiance, [W m-2 nm-1]
  !real(kind=RKIND), allocatable :: saflux(:,:,:)
  !! Spectral actinic flux, quanta s-1 nm-1 cm-2

contains        

    subroutine Tuv_driver(domain, iTimestep, press &
              , temp &
              , zt_ &
              , zm_ &
              , dzp &
              , rho &
              , pp &
              , coszr &
              , rlongup &
              , glat &
              , glon &
              , qv &
              , sfc_albedo &
              , nCells &
              , nVertLevels &
              ,myNum &
              ,jphoto &
              )
      !! ## Driver for TUV Radiation
      !!
      !! ## MONAN
      !!
      !! Author: rodrigues, L.F.
      !!
      !! E-mail: luflarois@gmail.com
      !!
      !! Date: 2026-03-04
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! Driver to adapt MONAN variables to the TUV radiation code.
      !!
      !! ** History**:
      !!
      !! - Itenizado_as_alterações_ao_longo_do_tempo (genérica)
      !!—
      !! ** Licence **:
      !!
      !! CC-GPL 3.0 License (https://creativecommons.org/licenses/GPL/3.0/)
      !!
      implicit none
      
      type(domain_type),intent(inout):: domain
      integer, intent(in) :: iTimestep
      integer,intent(in):: nCells, nVertLevels, mynum

      real,dimension(nVertLevels,nCells),intent(in) :: press
      !! Pressure field [Pa]
      real,dimension(nVertLevels,nCells),intent(in) :: temp
      !! Temperature field [K]
      real,dimension(nVertLevels,nCells),intent(in) :: zt_
      !! Geometric height of layer top [m]
      real,dimension(nVertLevels,nCells),intent(in) :: zm_
      !! Geometric height of layer midpoint [m]
      real,dimension(nVertLevels,nCells),intent(in) :: dzp
      !! Thickness of each layer [m]
      real,dimension(nVertLevels,nCells),intent(in) :: rho
      !! Air density at layer midpoint [kg/m3]
      real,dimension(nVertLevels,nCells),intent(in) :: pp
      !! Pressure at layer midpoint [Pa]  - SOundings
      real,dimension(nCells),intent(in)         :: coszr
      !! Cosine of solar zenith angle at each cell [unitless]
      real,dimension(nCells),intent(in)         :: rlongup
      !! !all-sky upwelling longwave flux at bottom-of-atmosphere
      real,dimension(nCells),intent(in)         :: glat
      !! latitude, south is negative [degrees]
      real,dimension(nCells),intent(in)         :: glon
      !! longitude, west is negative [degrees]
      real,dimension(nVertLevels,nCells),intent(in) :: qv
      !! !water vapor mixing ratio   [kg/kg
      real,dimension(nCells),intent(in) :: sfc_albedo

      real(kind=RKIND), intent(out) :: jphoto(nVertlevels, nCells, nr_photo)
 
      real(kind=RKIND) :: w1000, sig, sza
      integer :: n, i, j, k, k2, kvert, lv, lf, nz, j_ref
      real(kind=RKIND) :: tuv2carma(nWavelengths)  

      real(kind=RKIND), allocatable :: ssaer(:)
      !!
      real(kind=RKIND), allocatable :: co3FromTuv(:,:)
      !!
      integer, allocatable :: validCellsMap(:)
      !! Map of valid cells for radiation calculations - it depends on coszr

      integer :: nrad(nCells), narad(nCells)

      !!
      if (iTimestep == 1) then
          call InitTuv('./tuvData/', myNum, 'RELACS', nCells = nCells, nVertLevels = nVertLevels, nr_photo = nr_photo)
      end if
!      print *,'LFR-DBG - Inside TUV driver, timestep: ', iTimestep
      nz = nVertLevels - 1
      maxNRad = nVertLevels
      jphoto = 0.0

      if (iTimestep == 1) then
         ! Verify if NO2, SO2, and O3 are present in chem mechanism, and get their positions in the spc_name array
         call checkNo2Co2O3(nSpecies   = nspecies  , spc_name  = spc_name , no2present = no2present &
                       ,   so2present = so2present, o3Present = o3Present, iposno2    = iposNo2    &
                       ,      iposso2 = iposSo2   , iposo3    = iposO3)
         ! Map TUV wavelengths to Carma wavelengths (get the index of the Carma wavelength
           
         tuv2carma = 1
         do j_ref = 1, nw - 1
             do i = 1, nwave
                 w1000 = wave(i) * 1000.0
                 if (w1000 >= wu(j_ref)) exit 
                 tuv2carma(j_ref) = i
             end do
         end do

      end if

      do i = 1, nCells
         ! Get the number of rad levels for each column and max levels to be added (narad) above model top
         ! if model top pressure is greater than refPress.
         narad(i) = compute_rad_levs(prsnz = pp(nz-1,i), prsnzp = pp(nz,i), refPress = refPress)
         nrad(i) = nz + narad(i)
         maxNRad = max(maxNRad,nrad(i))
      end do
      call allocate_local_arrays(nCells, maxNRad, nWaveLengths, itimestep)
     
      zm(1) = 0.0
      do k= 2, nVertLevels
         zm(k) = zm_(k,1) - zm_(1,1)
      end do
      do k = 1, nVertLevels
         zt(k) = zm(k) - 0.5 * dzp(k,1)
      end do
      call compute_MC_Soundings(curr_julday)   
      do i = 1, nCells
         do k = 1, nVertLevels
            zml(k,i)  = zm(k)
            ztl(k,i)  = zt(k)
            if (k>1) then
               dzl(k,i) = zml(k,i) - zml(k-1,i)
            else
               dzl(k,i) = zt(k) - 0.0   
            end if
            prd(k,i)  = press(k,i)
            temprd(k,i) = temp(k,i)
            dair(k,i) = rho(k,i) 
            rv(k,i) = qv(k,i)
            o3l(k,i) = 0.0_RKIND !TBD: must be filled with chem O3
         end do
      end do  
      do i = 1, nCells
         if (nrad(i) == nVertLevels) cycle !If the top of model is above the max altitude for radiation calculations, 
                                      !skip the interpolation and column calculations, use from Model
         call compute_trop_ozone_column(nVertLevels = nVertLevels, nzrad = maxNRad, narad = narad(i), nrad = nrad(i) &
                                 ,   glat = glat(i), rlongup = rlongup(i), zm = zm, zt = zt  &
                                 ,   zml = zml(:,i), ztl = ztl(:,i), pl = prd(:,i), tl = temprd(:,i)  &
                                 ,    dl = dair(:,i), rl = rv(:,i), o3l = o3l(:,i), dzl = dzl(:,i))
      end do

      do i = 1, nCells
         so2col(i) = 0.0_RKIND
         no2col(i) = 0.0_RKIND
         do k = 1, nWaveLengths
            albedo(i,k) = sfc_albedo(i) !By now albedo doesn't depend on wavelength, but it need to be changed.
         end do   
         do k = 1, nrad(i)
            tlev(k,i) = temprd(k,i)
            zlevel(k,i) = zml(k,i) * 1.0e-3 !Convert to km
            dtcld(k,i,:) = 0.0_RKIND   !TBD: get cloud properties from model output when available
            omcld(k,i,:) = 0.0_RKIND
            gcld(k,i,:)  = 0.0_RKIND
            dtaer(k,i,:) = 0.0_RKIND
            omaer(k,i,:) = 0.0_RKIND
            gaer(k,i,:)  = 0.0_RKIND
            o3(k,i) = o3l(k,i)  * (fatmul / 48.00d0) * dzl(k,i) &
                                    * 1.d+3 * 1.d-4 !srf: ok
            cair(k,i)   = dair(k,i) * (fatmul / 28.96d0) * dzl(k,i) &
                                    * 1.d+3 * 1.d-4 !srf: ok
                 !- #molec[ar]/cm^3
            airden(k,i) = dair(k,i) * (fatmul / 28.96d0) &
                                    * 1.d+3 * 1.d-6  !srf: ok
         end do          
         do k = 1, nrad(i) - 1
            tlay(k,i) = 0.5 * (temprd(k,i) + temprd(k+1,i))
         end do      
         tlay(nrad(i),i) = tlay(nrad(i)-1,i) + (tlay(nrad(i)-1,i) - tlay(nrad(i)-2,i)) !Extrapolate temperature for the layer above nrad
         if (so2present) then
            do k = 1, nVertLevels
               so2col(i) =  so2col(i) + convert_kgkg_to_du(m = 0.0, rho = rho(k,i) & !TBD: mass must be gas value in kg/kg, currently 
                                                                                            !set to 0.0 for testing
                            , h = dzl(k,i), specie = 'SO2') ! Convert mixing ratio in kg/kg to Dobson Units
            end do
         end if 
         if (no2present) then
            do k = 1, nVertLevels
               no2col(i) = no2col(i) + convert_kgkg_to_du(m = 0.0, rho = rho(k,i) & !TBD: mass must be gas value in kg/kg, currently 
                                                                                            !set to 0.0 for testing
                            , h = dzl(k,i), specie = 'NO2') ! Convert mixing ratio in kg/kg to Dobson Units
            end do
         end if
         tco3(nrad(i),i) = o3(nrad(i),i)
         do k = nrad(i)-1,1,-1
            tco3(k,i)=tco3(k+1,i)+o3(k,i)
         end do          

         if (coszr(i) > 0.0_RKIND) then
            sza = acos(coszr(i))*f180PI
            call Tuv(mynum    = mynum                       &
                   , nstr     = zero                        & 
                   , alpha    = alpha                       &
                   , dirsun   = dirSun                      &
                   , difdn    = difdn                       &
                   , difup    = difup                       &
                   , esfact   = esFact                      &
                   , albedo   = albedo(i, :           )   & 
                   , nz       = nrad  (i              )   &  
                   , sza      = sza                      & 
                   , so2col   = so2col(i              )   & 
                   , no2col   = no2col(i              )   & 
                   , zLevel   = zlevel(1:nrad(i),i)   & 
                   , tlev     = tlev  (1:nrad(i),i)   & 
                   , tlay     = tlay  (1:nrad(i),i)   & 
                   , airden   = airden(1:nrad(i),i)   & 
                   , cair     = cAir  (1:nrad(i),i)   & 
                   , co3      = o3    (1:nrad(i),i)   & 
                   , tco3     = tco3  (1:nrad(i),i)   & 
                   , dtcld    = dtcld (1:nrad(i),i,:) & 
                   , omcld    = omcld (1:nrad(i),i,:) & 
                   , gcld     = gcld  (1:nrad(i),i,:) & 
                   , dtaer    = dtaer (1:nrad(i),i,:) & 
                   , omaer    = omaer (1:nrad(i),i,:) & 
                   , gaer     = gaer  (1:nrad(i),i,:) &
                   , valj     = valj  (1:nrad(i),i,:) &
                   )

         else
            valj  (1:nrad(i),i,:) = 0.0_RKIND
         end if
      end do !i loop

      do n=1, kj !loop over wavelength and altitude dependent
         if (xref(n) == 0) cycle
         do i = 1, nCells
            do k = 1, nVertLevels
               jphoto(k,i,xRef(n))=valj(k,i,n)
            end do
         end do
      end do

      !Adjust for relacs mechanism:
      jphoto(:,:,14) = 0.962055*jphoto(:,:,13)+0.0106247*jphoto(:,:,9)
      jphoto(:,:,16)=(12*jphoto(:,:,13))+(208*jphoto(:,:,15))

      call deallocate_local_arrays()

   end subroutine Tuv_driver


   subroutine allocate_local_arrays(nCells,maxNRad,nWaveLengths,itimestep)
      implicit none
      integer, intent(in) :: nCells, maxNRad, nWaveLengths, itimestep

      if (.not. allocated(tlev  )) allocate(tlev  (maxNRad,nCells))
      if (.not. allocated(tlay  )) allocate(tlay  (maxNRad,nCells))
      if (.not. allocated(zml   )) allocate(zml   (maxNRad,nCells))
      if (.not. allocated(ztl   )) allocate(ztl   (maxNRad,nCells))
      if (.not. allocated(prd   )) allocate(prd   (maxNRad,nCells))
      if (.not. allocated(temprd)) allocate(temprd(maxNRad,nCells))
      if (.not. allocated(dair  )) allocate(dair  (maxNRad,nCells))
      if (.not. allocated(rv    )) allocate(rv    (maxNRad,nCells))
      if (.not. allocated(dzl   )) allocate(dzl   (maxNRad,nCells))
      if (.not. allocated(o3l   )) allocate(o3l   (maxNRad,nCells))
      if (.not. allocated(zLevel)) allocate(zLevel(maxNRad,nCells))
      if (.not. allocated(albedo)) allocate(albedo(nCells        ,nWaveLengths))
      if (.not. allocated(dtcld )) allocate(dtcld (maxNRad,nCells,nWaveLengths))
      if (.not. allocated(omcld )) allocate(omcld (maxNRad,nCells,nWaveLengths))
      if (.not. allocated(gcld  )) allocate(gcld  (maxNRad,nCells,nWaveLengths))
      if (.not. allocated(dtaer )) allocate(dtaer (maxNRad,nCells,nWaveLengths))
      if (.not. allocated(omaer )) allocate(omaer (maxNRad,nCells,nWaveLengths))
      if (.not. allocated(gaer  )) allocate(gaer  (maxNRad,nCells,nWaveLengths))
      if (.not. allocated(cAir  )) allocate(cAir  (maxNRad,nCells))
      if (.not. allocated(o3    )) allocate(o3    (maxNRad,nCells))
      if (.not. allocated(airden)) allocate(airden(maxNRad,nCells))
      if (.not. allocated(tco3  )) allocate(tco3  (maxNRad,nCells))
      if (.not. allocated(so2col)) allocate(so2col(nCells))   
      if (.not. allocated(no2col)) allocate(no2col(nCells)) 
      if (.not. allocated(zt    )) allocate(zt    (maxNRad))
      if (.not. allocated(zm    )) allocate(zm    (maxNRad))
      if (itimestep == 1) then
         !if (.not.allocated(rate))   allocate(rate  (maxNRad,nCells,ks))
         if (.not.allocated(valj))   allocate(valj  (maxNRad,nCells,kj))
         !if (.not.allocated(sirrad)) allocate(sirrad(maxNRad,nCells,nWaveLengths))
         !if (.not.allocated(saflux)) allocate(saflux(maxNRad,nCells,nWaveLengths))
      end if

   end subroutine allocate_local_arrays

   subroutine deallocate_local_arrays()
      implicit none

      if (allocated(tlev  )) deallocate(tlev  )
      if (allocated(tlay  )) deallocate(tlay  )
      if (allocated(zml   )) deallocate(zml   )
      if (allocated(ztl   )) deallocate(ztl   )
      if (allocated(prd   )) deallocate(prd   )
      if (allocated(temprd)) deallocate(temprd)
      if (allocated(dair  )) deallocate(dair  )
      if (allocated(rv    )) deallocate(rv    )
      if (allocated(dzl   )) deallocate(dzl   )
      if (allocated(o3l   )) deallocate(o3l   )
      if (allocated(zLevel)) deallocate(zLevel)
      if (allocated(albedo)) deallocate(albedo)
      if (allocated(dtcld )) deallocate(dtcld )
      if (allocated(omcld )) deallocate(omcld )
      if (allocated(gcld  )) deallocate(gcld  )
      if (allocated(dtaer )) deallocate(dtaer )
      if (allocated(omaer )) deallocate(omaer )
      if (allocated(gaer  )) deallocate(gaer  )
      if (allocated(cAir  )) deallocate(cAir  )
      if (allocated(o3    )) deallocate(o3    )
      if (allocated(airden)) deallocate(airden)
      if (allocated(tco3  )) deallocate(tco3  )
      if (allocated(so2col)) deallocate(so2col)
      if (allocated(no2col)) deallocate(no2col)
      if (allocated(zt    )) deallocate(zt    )
      if (allocated(zm    )) deallocate(zm    )

   end subroutine deallocate_local_arrays

   function convert_kgkg_to_du(m, rho, H, specie) result(DU)
       implicit none
       real(kind = RKIND), intent(in) :: m       ! razão de mistura em massa (kg SO2 / kg ar)
       real(kind = RKIND), intent(in) :: rho     ! densidade do ar (kg/m³)
       real(kind = RKIND), intent(in) :: H       ! altura da coluna (m)
       real(kind = RKIND)            :: DU       ! resultado em Dobson Units
       character(len=*), intent(in) :: specie ! 'SO2' ou 'NO2'

       ! Constantes
       real(kind = RKIND), parameter :: NA = 6.02214076e23   ! Número de Avogadro (mol⁻¹)
       real(kind = RKIND), parameter :: M_SO2 = 64.0e-3      ! Massa molar do SO2 (kg/mol)
       real(kind = RKIND), parameter :: M_NO2 = 46.0e-3      ! Massa molar do NO2 (kg/mol)
       real(kind = RKIND), parameter :: DU_per_molec_cm2 = 2.69e16  ! moléculas/cm² por DU
       real(kind = RKIND), parameter :: cm2_per_m2 = 1.0e4   ! cm² por m²

       ! Variáveis intermediárias
       real(kind = RKIND) :: mass_gas_per_area   ! massa de SO2 por área (kg/m²)
       real(kind = RKIND) :: moles_gas_per_area  ! mols de SO2 por área (mol/m²)
       real(kind = RKIND) :: molecules_per_cm2   ! moléculas por cm²

       ! Cálculo da massa de SO₂ por unidade de área (kg/m²)
       mass_gas_per_area = m * rho * H

       ! Conversão para mols por área (mol/m²)
       if (specie == 'SO2') then
          moles_gas_per_area = mass_gas_per_area / M_SO2
       else if (specie == 'NO2') then
          moles_gas_per_area = mass_gas_per_area / M_NO2
       else
          print *, 'Error. Use only "SO2" or "NO2".'
          DU = -1.0  ! Valor de erro
          return
       end if

       ! Conversão para moléculas por cm²
       molecules_per_cm2 = moles_gas_per_area * NA / cm2_per_m2

       ! Conversão para Dobson Units
       DU = molecules_per_cm2 / DU_per_molec_cm2

   end function convert_kgkg_to_du

   subroutine compute_trop_ozone_column(nVertLevels, nzrad, narad, nrad, glat, rlongup, zm, zt &
                                       , zml, ztl, pl, tl, dl, rl, o3l, dzl)
      !! ## Compute the tropospheric ozone column and other variables needed for the TUV radiation code
      !!
      !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
      !! ## MONAN
      !!
      !! Author: Rodrigues, L.F.
      !!
      !! E-mail: luflarois@gmail.com
      !!
      !! Date: 2026-03-21
      !!
      !! #####Version: <>
      !!
      !! —
      !! **Full description**:
      !!
      !! At this point, subtropical, mid-latitude, sub-arctic,
      !! and arctic Mclatchy soundings have been interpolated between summer
      !! and winter values by time of year.  In this section of code,
      !! interpolate these 4 plus the all-year tropical sounding by latitude
      !!for the current i,j column in the grid.
      !!
      !! ** History**:
      !!
      !! - Itenizado_as_alterações_ao_longo_do_tempo (genérica)
      !!—
      !! ** Licence **:
      !!
      !! <img src='https://www.gnu.org/graphics/gplv3-127x51.png' width='63'>
      !!
      implicit none

      integer, intent(in) :: nVertLevels  
      integer, intent(in) :: nzrad, narad, nrad
      real(kind = RKIND), intent(in) :: glat
      real(kind = RKIND), intent(in) :: rlongup
      real(kind = RKIND), intent(in) :: zm(nzrad)
      real(kind = RKIND), intent(in) :: zt(nzrad)

      real(kind = RKIND), intent(inout) :: zml(nzrad)
      real(kind = RKIND), intent(inout) :: ztl(nzrad)
      real(kind = RKIND), intent(inout) :: pl(nzrad)
      real(kind = RKIND), intent(inout) :: tl(nzrad)
      real(kind = RKIND), intent(inout) :: dl(nzrad)
      real(kind = RKIND), intent(inout) :: rl(nzrad)
      real(kind = RKIND), intent(inout) :: o3l(nzrad)
      real(kind = RKIND), intent(inout) :: dzl(nzrad)

      integer, dimension(11), parameter :: latind = (/1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9/)
      real(kind=RKIND), dimension(12), parameter :: slat = (/ -90. , -70. , -60. , -45. &
                                                            , -25. , -15. ,  15. ,  25. &
                                                            ,  45. ,  60. ,  70. ,  90./)

      integer :: ind, index, lv, lf, k, lat
      real(kind = RKIND) :: wtnorth, wtsouth, wt, deltap

      do ind = 1, 11
          index = ind
          if (glat < slat(index + 1)) exit
      end do
      lat = latind(index)

      if (index == 1 .or. index == 6 .or. index == 11) then
         ! For pure arctic or tropical latitudes, assign sounding values without
          ! interpolation.
          do lv = 1, 33
             do lf = 1, 6
                mcol(lv, lf) = mclat(lv, lat, lf)
             end do
          end do
      else
         ! For other latitudes, linearly interpolate between soundings according
          ! to latitude bands defined in array `slat'.
          wtnorth = (glat - slat(index)) / (slat(index + 1) - slat(index))
          wtsouth = 1. - wtnorth
          do lv = 1, 33
             do lf = 1, 6
                mcol(lv, lf) = wtsouth * mclat(lv, lat, lf) &
                + wtnorth * mclat(lv, lat + 1, lf)
             end do
          end do
       end if

      ! Fill radiation column arrays (other than o3l) with variables from
       ! the model grid.
       ! Modif. by ALF-Sfreitas
       ! Original: deltap = (pl(nVertLevels-1) - 1500.) / float(narad)
      if (narad /= 0) then !IF (prsnz >= 3000.) THEN
         !-origdeltap = (pl(nVertLevels-1) - 1500.) / float(narad)
          deltap = (pl(nVertLevels - 1) - 200.) / float(narad)
      else
         deltap = 0.
      end if

      do k = nVertLevels, nRad
          pl(k) = pl(k - 1) - deltap
      end do
      ! Interpolate O3 from Mclatchy sounding to all levels in radiation column,
      ! and interpolate other variables (temperature, density,
      !vapor mixing ratio) to added levels.
      lv = 1
      do k = nVertLevels, nRad        
          do
             if (pl(k) > mcol(1, 2)) then
                o3l(k) = mcol(1, 5)
                if (k >= nVertLevels) then
                   tl(k) = mcol(1, 3)
                   rl(k) = mcol(1, 4)
                   dl(k) = mcol(1, 6)
                end if
                exit
             else if (pl(k) <= mcol(lv, 2) .and. pl(k) >= mcol(lv + 1, 2)) then
               wt = (pl(k) - mcol(lv, 2)) / (mcol(lv + 1, 2) - mcol(lv, 2))
                o3l(k) = mcol(lv, 5) + (mcol(lv + 1, 5) - mcol(lv, 5)) * wt
                if (k >= nVertLevels) then
                   tl(k) = mcol(lv, 3) + (mcol(lv + 1, 3) - mcol(lv, 3)) * wt
                   rl(k) = mcol(lv, 4) + (mcol(lv + 1, 4) - mcol(lv, 4)) * wt
                   dl(k) = mcol(lv, 6) + (mcol(lv + 1, 6) - mcol(lv, 6)) * wt
                end if
                exit
             else if (pl(k) < mcol(33, 2)) then
                o3l(k) = mcol(33, 5)
                if (k >= nVertLevels) then
                   tl(k) = mcol(33, 3)
                   rl(k) = mcol(33, 4)
                   dl(k) = mcol(33, 6)
                end if
                exit
             else if (pl(k) < mcol(lv + 1, 2)) then
                lv = lv + 1
                if (lv >= 33) then
                   print *, 'Lv greater or equal 33';
                   !call flush (6)
                   stop 'mclat1_fastjx'
                end if
                cycle
             else
                print *, 'pressure data improperly ordered';
                !call flush (6)
                stop 'mclat2_fastjx'
             end if
          end do
       end do

      ! Compute heights of added levels by hydrostatic integration.
       do k = nVertLevels, nRad
          ztl(k) = ztl(k - 1) + rgasog * (tl(k - 1) - tl(k) &
                  + (pl(k - 1) * tl(k) - pl(k) * tl(k - 1)) / (pl(k) - pl(k - 1)) &
                  * log(pl(k) / pl(k - 1)))
          zml(k - 1) = .5 * (ztl(k) + ztl(k - 1))
       end do
       zml(nRad) = 2. * ztl(nRad) - zml(nRad - 1)

      ! Compute dzl values.
       do k = 2, nRad
          dzl(k) = zml(k) - zml(k - 1)
       end do
       dzl(1) = zml(2)

      ! Fill surface values
       pl(1) = pl(2) + (zm(1) - zt(3)) / (zt(2) - zt(3)) * (pl(2) - pl(3))

       tl(1) = sqrt(sqrt(rlongup / stefan))

   end subroutine compute_trop_ozone_column


   subroutine compute_MC_Soundings(jday)
      !! ## Interpolate arctic, sub-arctic, mid-latitude, and subtropical, Mclatchy
      !!
      !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
      !! ## MONAN
      !!
      !! Author: Rodrigues, L.F.
      !!
      !! E-mail: luflarois@gmail.com
      !!
      !! Date: 2026-03-21
      !!
      !! #####Version: <>
      !!
      !! —
      !! **Full description**:
      !!
      !! Interpolate arctic, sub-arctic, mid-latitude, and subtropical, Mclatchy
      !! soundings between summer and winter values by time of year using cosine
      !! function.  Assume that extreme values occur on
      !! January 16 and 1/2 year later.
      !! (Done once per grid/node each radiation time.)
      !!
      !! ** History**:
      !!
      !! - Itenizado_as_alterações_ao_longo_do_tempo (genérica)
      !!—
      !! ** Licence **:
      !!
      !! <img src='https://www.gnu.org/graphics/gplv3-127x51.png' width='63'>
      !!
      implicit none

      real(kind=RKIND), intent(in) :: jday
      !! Julian day of the year (1-365)
      real(kind=RKIND) :: fjday, wtjan, wtjul
      integer :: lats, latn, isummer, iwinter, lv, lf

      fjday = jday
      wtjan = 0.5 * (1. + cos(6.283185 * (fjday-16.) / 365.))
      wtjul = 1. - wtjan
      do lats = 1,4
  	      latn = 10 - lats
  	      isummer = 2 * lats
  	      iwinter = isummer - 1
  	      do lv = 1,33
  	         do lf = 1,6
  		         mclat(lv,lats,lf) = wtjan * mcdat(lv,isummer,lf)  &
  		                             + wtjul * mcdat(lv,iwinter,lf)
  		         mclat(lv,latn,lf) = wtjan * mcdat(lv,iwinter,lf)  &
  		                           + wtjul * mcdat(lv,isummer,lf)
            end do
  	      end do
      end do

   end subroutine compute_MC_Soundings

   subroutine checkNo2Co2O3(nSpecies, spc_name, no2present, so2present, o3Present, iposno2, iposso2, iposo3)
      !! ## check if NO2 and CO2 are present in chem mechanism
      !!
      !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
      !! ## MONAN
      !!
      !! Author: USERRodrigues, L.F.
      !!
      !! E-mail: luflarois@gmail.com
      !!
      !! Date: 2026-03-20
      !!
      !! #####Version: <>
      !!
      !! —
      !! **Full description**:
      !!
      !! check if NO2 and CO2 are present in chem mechanism
      !!
      !! ** History**:
      !!
      !! - Itenizado_as_alterações_ao_longo_do_tempo (genérica)
      !!—
      !! ** Licence **:
      !!
      !! <img src='https://www.gnu.org/graphics/gplv3-127x51.png' width='63'>
      !!
      integer, intent(in) :: nSpecies
      character(len=*), dimension(nSpecies), intent(in) :: spc_name
      logical, intent(out) :: no2present, so2present, o3Present
      integer, intent(out) :: iposno2, iposso2, iposo3
      
      !!
      integer :: i

      no2present = .false.
      so2present = .false.
      o3Present  = .false.
      
      do i = 1, nSpecies
          if (trim(spc_name(i)) == 'NO2') then
              iposno2 = i
              no2present = .true.
          end if
          if (trim(spc_name(i)) == 'SO2') then
              iposso2 = i
              so2present = .true.
          end if
         if (trim(spc_name(i)) == 'O3') then
            iposo3 = i
            o3Present = .true.
         end if

      end do 

   end subroutine checkNo2Co2O3


   function compute_rad_levs(prsnz, prsnzp, refPress) result(narad)
      !! ## Verify if prsnz and prsnzp and compute naRad
      !!
      !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
      !! ## MONAN
      !!
      !! Author: rodrigues, L.F.
      !!
      !! E-mail: luflarois@gmail.com
      !!
      !! Date: 2026-03-19
      !!
      !! #####Version: 0.1.0
      !!
      !! —
      !! **Full description**:
      !!
      !! Verify if prsnz and prsnzp are less than refPress [Pa], and if so, set narad to zero. Otherwise, 
      !! compute the number of radiation levels to be added above. Also copy tropical sounding to mclat 
      !! array (done only once on each node).
      !!
      !! ** History**:
      !!
      !! - Itenizado_as_alterações_ao_longo_do_tempo (genérica)
      !!—
      !! ** Licence **:
      !!
      !! <img src='https://www.gnu.org/graphics/gplv3-127x51.png' width='63'>
      !!
      implicit none
      real(kind=RKIND), intent(in) :: prsnz, prsnzp
      real(kind=RKIND), intent(in) :: refPress
      ! Output
      integer :: narad

      integer :: lv, lf
      real(kind=RKIND) :: deltap
      

      ! Copy (since no time interpolation needed) tropical sounding to
      ! mclat array (done only once on each node).
      do lv = 1, 33
         do lf = 1, 6
            mclat(lv, 5, lf) = mcdat(lv, 9, lf)
         end do
      end do
!
      ! Compute number of levels to be added above model top from
      ! Mclatchy soundings.
      ! Base this number on prsnz and prsnzp, which was earlier computed from
      ! PI01DN(NNZP(1),1), so that it is the same for all compute nodes
      ! (done only once on each node).
      if (prsnz < refPress) then
!
         ! If prsnz, the pressure at the highest model prognostic level, is less
         ! than 3000 Pa, do not add any levels.
         narad = 1
!
      else
         ! If prsnzp is greater than 3000 Pa, add one or more radiation levels.
         ! Make the top level be 200 Pa.
!
         !srf-changed from 1500 to 200 Pa
         deltap = max(200. , &
         !mp
         (prsnz - prsnzp) / 2. , &
         !mp
         (prsnz - 200.1) / float(namax))
         narad = max(1, int((prsnz - 200.) / deltap))
!
      end if

   end function compute_rad_levs

end module ModTuv_driver