USE FULIOUMULTI
USE GENERATE_FULIOU_LEVELS ,only : gflq, generate_level_scheme
USE EXTRAS       ,only : getatmosphere, aer_scale_hgt
USE CALIPSO_OUTPUT, only : pack_sky,print_pack_sky,skyp,SKYP_TYPE

USE ICEDIRSFC,only: tau_uc !! Debug Diagnostic


! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    implicit none

    real  aot_by_band,aotf,wlf
    common /aotbyband/ aot_by_band(18) 
    common /tau_spline_aot/ aotf(18),wlf(18)
      
    TYPE (SKYP_TYPE) ut,tu

    real psfc
    !integer ii,jj,kk,nn,mm,i,is
    real psel(6)

    real                :: z1, z2

    integer             :: ii
    integer             :: jj
    integer             :: kk
    integer             :: mm
    integer             :: nn

    integer             :: i_num_cod
    real                :: min_cod
    real                :: max_cod
    real                :: delta_cod
    real                :: local_cod

    integer             :: i_num_sza
    real                :: min_sza
    real                :: max_sza
    real                :: delta_sza
    real                :: local_sza
  
    integer             :: i_num_cldfrac
    real                :: min_cldfrac
    real                :: max_cldfrac
    real                :: delta_cldfrac
    real                :: local_cldfrac
  
    integer             :: i_num_aod
    real                :: min_aod
    real                :: max_aod
    real                :: delta_aod
    real                :: local_aod
 
    integer             :: i_num_sfcalb

    integer             :: io8
    integer             :: istatus

    character(len = 80) :: filename
 
    real                :: alb(3)

    real(kind = 8), dimension(:), allocatable :: pval_ltop 
    real(kind = 8), dimension(:), allocatable :: pval_lbot
    real(kind = 8), dimension(:), allocatable :: zval_ltop 
    real(kind = 8), dimension(:), allocatable :: zval_lbot
    real(kind = 8), dimension(:), allocatable :: cldidx_top
    real(kind = 8), dimension(:), allocatable :: cldidx_bot
    real(kind = 8), dimension(:), allocatable :: aeridx_top
    real(kind = 8), dimension(:), allocatable :: aeridx_bot
    real(kind = 8), dimension(:), allocatable :: aprof

    real(kind = 8), dimension(:), allocatable :: alb_data 
    real(kind = 8), dimension(:), allocatable :: cod_data 
    real(kind = 8), dimension(:), allocatable :: sza_data 
    real(kind = 8), dimension(:), allocatable :: cdf_data 
    real(kind = 8), dimension(:), allocatable :: aod_data 
    real(kind = 8), dimension(:,:,:,:,:), allocatable :: swf_clr_data 
    real(kind = 8), dimension(:,:,:,:,:), allocatable :: swf_total_data 
    real(kind = 8), dimension(:,:,:,:,:), allocatable :: swf_pristine_data 
    real(kind = 8), dimension(:,:,:,:,:), allocatable :: swf_total_no_aer_data 

    !integer             :: itest

    !include 'test_head.h'

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

  !real                :: base_alb(3)
  !real                :: alb_interval(3)
  !real                :: cld(6)
  
  call set_default_options_fu ! Sets some of the more obsure inputs to reasonable values.
  
  fi%lscm   = .false. 
  fi%lscm(1:2)   = .true. 
  fi%lscm   = .true. 
  !InPut Profile Assignment
   !call getatmosphere('../../testatms/jmls.lay ', &
   !call getatmosphere('../../testatms/jmlw.lay ', &
   !call getatmosphere('../../testatms/jsas.lay ', &
  call getatmosphere('../../testatms/jsaw.lay ', &
  ! call getatmosphere('./testatms/jmls.lay ',&
    FI%VI%nlev,&
    FI%VI%pp,&
    FI%VI%pt,&
    FI%VI%ph,&
    FI%VI%po,&
    FI%pts)
  FI%VI%nlev = FI%VI%nlev+1  ! LAYER(getatm) to LEVEL

  FI%VI%hsfc = 0.00 !! SURFACE GEOPOTENTIAL OF FI%VI profile
! FI%VI%hsfc = 1600 !! SURFACE GEOPOTENTIAL OF FI%VI profile

 !gflq%hsfc = 1500. !Meters Surface elev. of ACTUAL FOV... to nearest 120m Multiple
  gflq%hsfc = 2.0! 
  

  gflq%mode = 'CALIP'
  gflq%mode = 'CERE3'

  gflq%nld =4
  gflq%internal_levels(1:4) = (/70.,200.,500.,850./)

  fi%HYBRID_SW_SOLVER =.true. !200130802 SYNI Ed4 $S Clear , 2S HOMO CLD , GWTSA INHOM Cld
  !fi%HYBRID_SW_SOLVER =.false. !checks isksolve,fourssl

  fi%isksolve= 1  ! Solver Method (0=fu 1=gwtsa) 
  
  fi%ss	   = 1365 ! Solar Constant wm-2
  !fi%u0      =  0.5 ! Cosine Solar Zenith Angle
  fi%ur      =  0.8 ! Cosine View Zenith Angle (for IR Radiance)
  
  
  
  !-------Cnd 2
  fi%wp_hgt_flag = 0  ! Constant lwc with height
  !fi%wp_hgt_flag = 1  ! Water Cloud  Top thicker than Base
  !fi%wp_hgt_flag = 2  ! Ice Cloud  Bottom thicker than top
  
  
  
  !Surface Properties --------------------------------------------------
  
  ! ocean albedo = 0.06
  ! ice albedo = 0.61
  ! land albedo = 0.25
  
  !Allow different albedos for Aerosol Vs. NO Aerosol cases , And for each Clear/Cloud Conditions
  !!#!fi%sfcalb(1:18,1,0)  = 0.06 ! Clear sky -Spectral Surface Albedo SW
  !!#!fi%sfcalb(1:18,2,0)  = 0.06 ! Pristine sky -Spectral Surface Albedo SW
  !!#!fi%sfcalb(1:18,1,1:)  = 0.06  ! CLOUDY w/AOT  sky -Spectral Surface Albedo SW
  !!#!fi%sfcalb(1:18,2,1:)  = 0.06  ! CLOUDY w/o AOT sky -Spectral Surface Albedo SW
  
  fi%ee(1:12)  = 0.99 ! Spectral Surface Emissivity LW
  
  !Aerosols ------------------------------------------------------------
  fi%nac      = 1   ! 2 aerosol types 
  !fi%itps(1)	     = 11	   ! Continental see types (1-18)
  fi%itps(1)	     = 11      ! Soot
  !fi%itps(2)	     = 1	   
  !fi%itps(2)	     = 11	   ! Soot	  see types (1-18)
  
  fi%n_atau	      = 1	   ! 1 wavelength input for aerosols
  fi%a_wli(1)	      = 0.550	   ! AOT wavelength(microns) of a_taus
  
  !fi%a_taus(1,2)	      =  0.81	   ! AOT for constituent 2
  
  !----------------------------------------------------------------------
 
  min_cod = 0.00
  max_cod = 50.0
  delta_cod = 2.5
  i_num_cod = (max_cod - min_cod + delta_cod) / delta_cod

  min_sza = 35.0
  max_sza = 85.0
  delta_sza = 2.5
  i_num_sza = (max_sza - min_sza + delta_sza) / delta_sza
  
  min_cldfrac = 0.0
  max_cldfrac = 1.0
  delta_cldfrac = 0.2
  i_num_cldfrac = (max_cldfrac - min_cldfrac + delta_cldfrac) / delta_cldfrac
  
  min_aod = 0.0
  max_aod = 2.0
  delta_aod = 0.1
  i_num_aod = (max_aod - min_aod + delta_aod) / delta_aod
 
  i_num_sfcalb = 3 
  alb(1) = 0.06
  alb(2) = 0.11
  alb(3) = 0.61

  ! Allocate dimension and data arrays
  ! ----------------------------------
  allocate(alb_data(i_num_sfcalb))
  allocate(cod_data(i_num_cod))
  allocate(sza_data(i_num_sza))
  allocate(cdf_data(i_num_cldfrac))
  allocate(aod_data(i_num_aod))

  ! alb, sza, cod, cldfrac, aod
  allocate(swf_clr_data(&
    i_num_aod, i_num_cldfrac, i_num_sza, i_num_cod, i_num_sfcalb))
  allocate(swf_total_data(&
    i_num_aod, i_num_cldfrac, i_num_sza, i_num_cod, i_num_sfcalb))
  allocate(swf_pristine_data(&
    i_num_aod, i_num_cldfrac, i_num_sza, i_num_cod, i_num_sfcalb))
  allocate(swf_total_no_aer_data(&
    i_num_aod, i_num_cldfrac, i_num_sza, i_num_cod, i_num_sfcalb))

  alb_data(:) = -9.
  cod_data(:) = -9.
  sza_data(:) = -9.
  cdf_data(:) = -9.
  aod_data(:) = -9.
  swf_clr_data(:,:,:,:,:) = -9.
  swf_total_data(:,:,:,:,:) = -9.
  swf_pristine_data(:,:,:,:,:) = -9.
  swf_total_no_aer_data(:,:,:,:,:) = -9.

  ! alb: mm
  !   sza: ii
  !     cod: jj
  !       cldfrac: kk
  !         aod: nn

  ! To access the cosine of the solar zenith angle:
  ! fi%u0 
  ! To access the cloud optical depth at cloud type number 3 with only 1 overlying layer:
  ! fi%fc(3)%tau_vis(1)
 
  !!#!base_alb(1) = 0.066
  !!#!base_alb(2) = 0.25
  !!#!base_alb(3) = 0.71
  !!#!alb_interval(1) = 0.003
  !!#!alb_interval(2) = 0.03
  !!#!alb_interval(3) = 0.02
  !!#!
  !!#!cld(1) = 0.0
  !!#!cld(2) = 0.2
  !!#!cld(3) = 0.4
  !!#!cld(4) = 0.6
  !!#!cld(5) = 0.8
  !!#!cld(6) = 1.0

  ! Structure of output SWUP values:
  !  
  ! ftoa(1)%swup - SWF_CLR
  ! ftoa(2)%swup - SWF_TOTAL
  ! ftoa(3)%swup - SWF_PRISTINE
  ! ftoa(4)%swup - SWF_TOTAL_NO_AER

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! BEGIN AEROSOL TEST CHUNK
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

  local_sza = 60.0
  local_cod = 10.0
  local_cldfrac = 0.4
  local_aod = 1.0

  mm = 1

  fi%fc(1)%novl = 1
  fi%fc(1)%cldfrac = local_cldfrac
  fi%fc(1)%tau_vis(1) = local_cod
  fi%fc(1)%sc(1)%mn_lin_tau = local_cod !fi%fc(1)%tau_vis(2) 
  fi%vd%cldpres(1:2, 1,1) = (/904,925/)

  fi%u0      = cos((local_sza * 3.14159265) / 180.  )   ! Cosine Solar Zenith Angle

  fi%a_taus(1,1)      =  local_aod   ! AOT for constituent 1

  fi%fc(1)%rphase(1)    =  1.0    ! Cloud Phase 1=Water 2=Ice
  fi%fc(1)%de(1) = 60.
  fi%fc(1)%re(1) = 15.

  fi%sfcalb(1:18,1,0)   = alb(mm) ! Clear sky -Spectral Surface Albedo SW
  fi%sfcalb(1:18,2,0)   = alb(mm) ! Pristine sky -Spectral Surface Albedo SW
  fi%sfcalb(1:18,1,1:)  = alb(mm)  ! CLOUDY w/AOT  sky -Spectral Surface Albedo SW
  fi%sfcalb(1:18,2,1:)  = alb(mm)  ! CLOUDY w/o AOT sky -Spectral Surface Albedo SW
  
  !!#!fi%fc(1)%cldfrac   = cld(ii)    ! Cloud Fraction (0-1) 
  
  call generate_level_scheme !! Define model Fixed layer structure pre-cloud by fixed DZ intervals...
       !     call print_vla_in 
       
  call prepare_model_profile_fu !! CALL After all FI%VD and FI%VI structures are defined.
  call vla_interface_fu     ! uses FI%VO !! Assign Model ATM Profile and CLD Levels
       !     call print_vla_out
       
       !    Aerosol Profile (after fi%pp is created )-----------------------------
       
  !!#!call aer_scale_hgt(fi%nv,fi%pp,9.0,fi%aprofs(1:fi%nv,1) )

  !call aer_scale_hgt(fi%nv,fi%pp,1.0,fi%aprofs(1:fi%nv,2) )
       !     RADIATVE TRANSFER --------------------------------------------------
       
            !call print_in_fu		   ! PRINTS INPUTS  AS ASCII 
            
  call rad_multi_fu  ! CALL THE CODE !!!

  do ii = 1, fi%nv
    write(*,*) ii, fi%pp(ii), fi%pp(ii+1), fi%aprofs(ii,1)
  enddo

  !write(*,*) 'HERE:',fi%aprofs(19:23,1)
 
  write(*,'(8(f8.3, 2x))') local_sza, local_cod, local_cldfrac, &
      local_aod, ftoa(1)%swup, ftoa(2)%swup, ftoa(3)%swup, &
      ftoa(4)%swup 
 
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! END AEROSOL TEST CHUNK
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

  ! Prep the profile variables
  ! --------------------------
  allocate(pval_ltop(fi%nv))
  allocate(pval_lbot(fi%nv))
  allocate(zval_ltop(fi%nv))
  allocate(zval_lbot(fi%nv))

  fi%aprofs(:,1) = 0.0
  fi%aprofs(21:25,1) = 20.0

  write(*,*) fi%nv

  do ii=1,fi%nv
  
    z1= 8.0* log( fi%pp(fi%nv+1) / fi%pp(ii))
    z2= 8.0* log( fi%pp(fi%nv+1) / fi%pp(ii+1))

    write(*,*) ii, z1, z2, fi%pp(ii), fi%pp(ii+1)

    pval_ltop(ii) = fi%pp(ii)
    pval_lbot(ii) = fi%pp(ii+1)

  enddo

  write(*,*) 
  write(*,*) fi%pp(:fi%nv)
  write(*,*) 
  write(*,*) fi%pp(2:fi%nv+1)
  write(*,*) 


  !!#!! Allocate dimension and data arrays
  !!#!! ----------------------------------
  !!#!
  !!#!do mm = 1, i_num_sfcalb

  !!#!  ! Insert into array
  !!#!  alb_data(mm) = alb(mm)

  !!#!  ! Open an output file
  !!#!  ! -------------------

  !!#!  !!#!if(fi%itps(1) > 9) then
  !!#!  !!#!  write(filename,'(a25,i2,a4,i1,a14)') 'fuliou_cloud_test_aertype',fi%itps(1),&
  !!#!  !!#!      '_alb', mm,'_loftplume.txt'
  !!#!  !!#!else
  !!#!  !!#!  write(filename,'(a25,i1,a4,i1,a14)') 'fuliou_cloud_test_aertype',fi%itps(1),&
  !!#!  !!#!      '_alb', mm,'_loftplume.txt'
  !!#!  !!#!endif

  !!#!  !!#!write(*,*) filename

  !!#!  !!#!open(io8, file = filename, iostat = istatus)
  !!#!  !!#!if(istatus /= 0) then
  !!#!  !!#!  write(*,*) "ERROR: Could not open output file"
  !!#!  !!#!  return
  !!#!  !!#!endif

  !!#!  !!#!write(io8,*) 'SolarZenith ','CloudOptDepth ', 'CloudFrac ','AOD ', &
  !!#!  !!#!    'SWF_CLR   ','SWF_TOTAL   ', 'SWF_PRISTINE   ','SWF_TOTAL_NO_AER'

  !!#!  ! Loop over the SZAs
  !!#!  do ii = 1, i_num_sza
  !!#!    write(*,*) ii
  !!#!    !fi%fc(1)%cldfrac   = 0.00000    ! Cloud Fraction (0-1) 
  !!#!    local_sza = min_sza + delta_sza * (ii - 1)

  !!#!    ! Insert into array
  !!#!    sza_data(ii) = local_sza

  !!#!    do jj = 1, i_num_cod
  !!#!      local_cod = min_cod + delta_cod * (jj - 1)

  !!#!      ! Insert into array
  !!#!      cod_data(jj) = local_cod

  !!#!      if(local_cod == 0) then
  !!#!        local_cod = 0.01
  !!#!      endif
  !!#!      do kk = 1, i_num_cldfrac
  !!#!        local_cldfrac = min_cldfrac + delta_cldfrac * (kk - 1)

  !!#!        ! Insert into array
  !!#!        cdf_data(kk) = local_cldfrac

  !!#!        do nn = 1, i_num_aod
  !!#!          local_aod = min_aod + delta_aod * (nn - 1)

  !!#!          if(local_aod == 0) then
  !!#!            local_aod = 0.01
  !!#!          endif

  !!#!          ! Insert into array
  !!#!          aod_data(nn) = local_aod

  !!#!          fi%fc(1)%novl = 1
  !!#!          fi%fc(1)%cldfrac = local_cldfrac
  !!#!          fi%fc(1)%tau_vis(1) = local_cod
  !!#!          fi%fc(1)%sc(1)%mn_lin_tau = local_cod !fi%fc(1)%tau_vis(2) 
  !!#!          fi%vd%cldpres(1:2, 1,1) = (/904,925/)

  !!#!          fi%u0      = cos((local_sza * 3.14159265) / 180.  )   ! Cosine Solar Zenith Angle

  !!#!          fi%a_taus(1,1)      =  local_aod   ! AOT for constituent 1

  !!#!          !!#!fi%fc(1)%dpi%ldpi = .false.
  !!#!          !!#!fi%fc(1)%cldfrac   = 0.00000    ! Cloud Fraction (0-1) 
  !!#!          !!#!fi%fc(1)%novl      =   2 
  !!#!          !!#!fi%fc(1)%novl      =   1 
  !!#!          !!#!
  !!#!          !!#!FI%VD%cldpres(1:2, 1,1) = (/200,400/)
  !!#!          !!#!!FI%VD%cldpres(1:2, 1,1) = (/400,800/)
  !!#!          !!#!
  !!#!          fi%fc(1)%rphase(1)    =  1.0    ! Cloud Phase 1=Water 2=Ice
  !!#!          fi%fc(1)%de(1) = 60.
  !!#!          fi%fc(1)%re(1) = 15.
  !!#!          !!#!
  !!#!          !!#!!fi%fc(1)%asp(1) = exp(iasp*0.1) !! Fu 20006 Ice AspectRatio !!!!! NEW FOR 20010130
  !!#!          !!#!
  !!#!          !!#!
  !!#!          !!#!fi%fc(1)%tau_vis(1)       = 10.00	    ! Cloud Visible Optical Depth ( Minnis)
  !!#!          !!#!fi%fc(1)%sc(1)%mn_lin_tau =  fi%fc(1)%tau_vis(1) *1.15
  !!#!          !!#!
  !!#!          !!#!
  !!#!          !!#!!-----
  !!#!          !!#!fi%fc(1)%rphase(2)    =  1.0    ! Cloud Phase 1=Water 2=Ice
  !!#!          !!#!fi%fc(1)%de(2) = 70.
  !!#!          !!#!fi%fc(1)%re(2) = 10.
  !!#!          !!#!
  !!#!          !!#!!fi%fc(1)%asp(1) = exp(iasp*0.1) !! Fu 20006 Ice AspectRatio !!!!! NEW FOR 20010130
  !!#!          !!#!
  !!#!          !!#!fi%fc(1)%tau_vis(2)       = 30	    ! Cloud Visible Optical Depth ( Minnis)
  !!#!          !!#!fi%fc(1)%sc(2)%mn_lin_tau =  fi%fc(1)%tau_vis(2) 
  !!#!          !!#!
  !!#!          !!#!fi%fc(1)%tau_vis(2)       = 1E-20	    ! Cloud Visible Optical Depth ( Minnis)
  !!#!          !!#!fi%fc(1)%sc(2)%mn_lin_tau = 1E-20 !fi%fc(1)%tau_vis(2) 

  !!#!          !Allow different albedos for Aerosol Vs. NO Aerosol cases , And for each Clear/Cloud Conditions
  !!#!          fi%sfcalb(1:18,1,0)   = alb(mm) ! Clear sky -Spectral Surface Albedo SW
  !!#!          fi%sfcalb(1:18,2,0)   = alb(mm) ! Pristine sky -Spectral Surface Albedo SW
  !!#!          fi%sfcalb(1:18,1,1:)  = alb(mm)  ! CLOUDY w/AOT  sky -Spectral Surface Albedo SW
  !!#!          fi%sfcalb(1:18,2,1:)  = alb(mm)  ! CLOUDY w/o AOT sky -Spectral Surface Albedo SW
  !!#!  
  !!#!          !!#!fi%fc(1)%cldfrac   = cld(ii)    ! Cloud Fraction (0-1) 
  !!#!    
  !!#!          call generate_level_scheme !! Define model Fixed layer structure pre-cloud by fixed DZ intervals...
  !!#!     !     call print_vla_in 
  !!#!     
  !!#!          call prepare_model_profile_fu !! CALL After all FI%VD and FI%VI structures are defined.
  !!#!          call vla_interface_fu     ! uses FI%VO !! Assign Model ATM Profile and CLD Levels
  !!#!     !     call print_vla_out
  !!#!     
  !!#!     !    Aerosol Profile (after fi%pp is created )-----------------------------
  !!#!     
  !!#!          !call aer_scale_hgt(fi%nv,fi%pp,10.0,fi%aprofs(1:fi%nv,1) )
  !!#!          !call aer_scale_hgt(fi%nv,fi%pp,10.0,fi%aprofs(1:fi%nv,2) )
  !!#!          fi%aprofs(:,1) = 0.0
  !!#!          fi%aprofs(21:25,1) = 20.0
  !!#!     !     RADIATVE TRANSFER --------------------------------------------------
  !!#!     
  !!#!          !call print_in_fu		   ! PRINTS INPUTS  AS ASCII 
  !!#!          
  !!#!          call rad_multi_fu  ! CALL THE CODE !!!

  !!#!          swf_clr_data(nn,kk,jj,ii,mm) = ftoa(1)%swup
  !!#!          swf_total_data(nn,kk,jj,ii,mm) = ftoa(2)%swup
  !!#!          swf_pristine_data(nn,kk,jj,ii,mm) = ftoa(3)%swup
  !!#!          swf_total_no_aer_data(nn,kk,jj,ii,mm) = ftoa(4)%swup

  !!#!          !!#!write(io8,'(8(f8.3, 2x))') local_sza, local_cod, local_cldfrac, &
  !!#!          !!#!    local_aod, ftoa(1)%swup, ftoa(2)%swup, ftoa(3)%swup, &
  !!#!          !!#!    ftoa(4)%swup 
  !!#!              
  !!#!          !write(*,*) ftoa(1)%swup, ftoa(2)%swup, ftoa(3)%swup, ftoa(4)%swup
 
  !!#!          !write(*,*) cld(ii), base_alb(jj), alb(jj) + alb_interval(jj) * (kk - 1), &
  !!#!          !    ftoa(1)%swup, ftoa(2)%swup, ftoa(3)%swup, ftoa(4)%swup
  !!#!          !write(*,*) cld(ii), alb(kk), ftoa(1)%swup, ftoa(2)%swup, ftoa(3)%swup, ftoa(4)%swup
  !!#!        enddo
  !!#!      enddo
  !!#!    enddo 
  !!#!  enddo

  !!#!  !!#!close(io8)

  !!#!enddo 

  !!#!write(*,*) alb_data
  !!#!write(*,*) sza_data
  !!#!write(*,*) cod_data
  !!#!write(*,*) cdf_data
  !!#!write(*,*) aod_data

  !!#!write(*,*) swf_total_data(1,3,3,3,3)

   !!!call generate_level_scheme !! Define model Fixed layer structure pre-cloud by fixed DZ intervals...
  !!!! call print_vla_in 
  
   !!!call prepare_model_profile_fu !! CALL After all FI%VD and FI%VI structures are defined.
   !!!call vla_interface_fu     ! uses FI%VO !! Assign Model ATM Profile and CLD Levels
  !!!! call print_vla_out
  
  !!!!Aerosol Profile (after fi%pp is created )-----------------------------
  
   !!!call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,1) )
   !!!call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,2) )
  !!!! RADIATVE TRANSFER --------------------------------------------------
  
   !!!!call print_in_fu		   ! PRINTS INPUTS  AS ASCII 
   !!!
   !!!call rad_multi_fu  ! CALL THE CODE !!!
  
   !!!!call print_out_fu		   ! PRINTS Lots of OUTPUTS  AS ASCII
  
   !!!write(*,*) ' MY OUTPUT'
   !!!write(*,*) fi%pts, ftoa(1)%swup, ftoa(2)%swup, ftoa(3)%swup, ftoa(4)%swup
  
   !call pack_sky
   !
   !call print_pack_sky

  deallocate(alb_data)
  deallocate(cod_data)
  deallocate(sza_data)
  deallocate(cdf_data)
  deallocate(aod_data)

  ! alb, sza, cod, cldfrac, aod
  deallocate(swf_clr_data)
  deallocate(swf_total_data)
  deallocate(swf_pristine_data)
  deallocate(swf_total_no_aer_data)
 
  deallocate(pval_ltop)
  deallocate(pval_lbot)
  deallocate(zval_ltop)
  deallocate(zval_lbot)
 
  stop ' Simple.f90 normal end'

end
