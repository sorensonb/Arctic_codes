USE FULIOUMULTI
USE GENERATE_FULIOU_LEVELS ,only : gflq, generate_level_scheme
USE EXTRAS       ,only : getatmosphere, aer_scale_hgt
USE CALIPSO_OUTPUT, only : pack_sky,print_pack_sky,skyp,SKYP_TYPE

USE ICEDIRSFC,only: tau_uc !! Debug Diagnostic
implicit none
 real  aot_by_band,aotf,wlf
 common /aotbyband/ aot_by_band(18) 
 common /tau_spline_aot/ aotf(18),wlf(18)
      

TYPE (SKYP_TYPE) ut,tu

real psfc
integer ii,jj,kk,i,is
real psel(6)

real                :: alb(3)
real                :: base_alb(3)
real                :: alb_interval(3)
real                :: cld(6)

 call set_default_options_fu ! Sets some of the more obsure inputs to reasonable values.
fi%lscm	   = .false. 
fi%lscm(1:2)	   = .true. 
fi%lscm	   = .true. 
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
fi%u0      =  0.5 ! Cosine Solar Zenith Angle
fi%ur      =  0.8 ! Cosine View Zenith Angle (for IR Radiance)



!-------Cnd 2
fi%wp_hgt_flag = 0  ! Constant lwc with height
!fi%wp_hgt_flag = 1  ! Water Cloud  Top thicker than Base
!fi%wp_hgt_flag = 2  ! Ice Cloud  Bottom thicker than top


fi%fc(1)%dpi%ldpi = .false.
fi%fc(1)%cldfrac   = 0.00000    ! Cloud Fraction (0-1) 
fi%fc(1)%novl      =   2 
fi%fc(1)%novl      =   1 

FI%VD%cldpres(1:2, 1,1) = (/200,400/)
FI%VD%cldpres(1:2, 1,2) = (/704,725/)
!FI%VD%cldpres(1:2, 1,1) = (/400,800/)

fi%fc(1)%rphase(1)    =  2.0    ! Cloud Phase 1=Water 2=Ice
fi%fc(1)%de(1) = 60.
fi%fc(1)%re(1) = 15.

!fi%fc(1)%asp(1) = exp(iasp*0.1) !! Fu 20006 Ice AspectRatio !!!!! NEW FOR 20010130


fi%fc(1)%tau_vis(1)       = 10.00	    ! Cloud Visible Optical Depth ( Minnis)
fi%fc(1)%sc(1)%mn_lin_tau =  fi%fc(1)%tau_vis(1) *1.15


!-----
fi%fc(1)%rphase(2)    =  1.0    ! Cloud Phase 1=Water 2=Ice
fi%fc(1)%de(2) = 70.
fi%fc(1)%re(2) = 10.

!fi%fc(1)%asp(1) = exp(iasp*0.1) !! Fu 20006 Ice AspectRatio !!!!! NEW FOR 20010130

fi%fc(1)%tau_vis(2)       = 30	    ! Cloud Visible Optical Depth ( Minnis)
fi%fc(1)%sc(2)%mn_lin_tau =  fi%fc(1)%tau_vis(2) 

fi%fc(1)%tau_vis(2)       = 1E-20	    ! Cloud Visible Optical Depth ( Minnis)
fi%fc(1)%sc(2)%mn_lin_tau = 1E-20 !fi%fc(1)%tau_vis(2) 


!Surface Properties --------------------------------------------------

! ocean albedo = 0.06
! ice albedo = 0.61
! land albedo = 0.25

!Allow different albedos for Aerosol Vs. NO Aerosol cases , And for each Clear/Cloud Conditions
fi%sfcalb(1:18,1,0)  = 0.06 ! Clear sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,2,0)  = 0.06 ! Pristine sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,1,1:)  = 0.06  ! CLOUDY w/AOT  sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,2,1:)  = 0.06  ! CLOUDY w/o AOT sky -Spectral Surface Albedo SW

fi%ee(1:12)  = 0.99 ! Spectral Surface Emissivity LW

!Aerosols ------------------------------------------------------------
fi%nac		     = 1	   ! 2 aerosol types 
fi%itps(1)	     = 11	   ! Continental see types (1-18)
!fi%itps(2)	     = 1	   
!fi%itps(2)	     = 11	   ! Soot	  see types (1-18)

fi%n_atau	      = 1	   ! 1 wavelength input for aerosols
fi%a_wli(1)	      = 0.641	   ! AOT wavelength(microns) of a_taus
fi%a_taus(1,1)	      =  0.01	   ! AOT for constituent 1

!fi%a_taus(1,2)	      =  0.81	   ! AOT for constituent 2

!----------------------------------------------------------------------

base_alb(1) = 0.066
base_alb(2) = 0.25
base_alb(3) = 0.71
alb(1) = 0.06
alb(2) = 0.11
alb(3) = 0.61
alb_interval(1) = 0.003
alb_interval(2) = 0.03
alb_interval(3) = 0.02

cld(1) = 0.0
cld(2) = 0.2
cld(3) = 0.4
cld(4) = 0.6
cld(5) = 0.8
cld(6) = 1.0


write(*,*) 'CloudFrac   ','Base_Albedo   ', 'Albedo   ','SWF_CLR   ','SWF_TOTAL   ', &
    'SWF_PRISTINE   ','SWF_TOTAL_NO_AER'
do ii = 1, 6
  fi%fc(1)%cldfrac   = 0.00000    ! Cloud Fraction (0-1) 
  do jj = 1, 3
    do kk = 1, 11
      !Allow different albedos for Aerosol Vs. NO Aerosol cases , And for each Clear/Cloud Conditions
      fi%sfcalb(1:18,1,0)   = alb(jj) + alb_interval(jj) * (kk - 1) ! Clear sky -Spectral Surface Albedo SW
      fi%sfcalb(1:18,2,0)   = alb(jj) + alb_interval(jj) * (kk - 1) ! Pristine sky -Spectral Surface Albedo SW
      fi%sfcalb(1:18,1,1:)  = alb(jj) + alb_interval(jj) * (kk - 1)  ! CLOUDY w/AOT  sky -Spectral Surface Albedo SW
      fi%sfcalb(1:18,2,1:)  = alb(jj) + alb_interval(jj) * (kk - 1)  ! CLOUDY w/o AOT sky -Spectral Surface Albedo SW

      fi%fc(1)%cldfrac   = cld(ii)    ! Cloud Fraction (0-1) 
  
      call generate_level_scheme !! Define model Fixed layer structure pre-cloud by fixed DZ intervals...
   !   call print_vla_in 
   
      call prepare_model_profile_fu !! CALL After all FI%VD and FI%VI structures are defined.
      call vla_interface_fu     ! uses FI%VO !! Assign Model ATM Profile and CLD Levels
   !   call print_vla_out
   
   !  Aerosol Profile (after fi%pp is created )-----------------------------
   
      call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,1) )
      call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,2) )
   !   RADIATVE TRANSFER --------------------------------------------------
   
      !call print_in_fu		   ! PRINTS INPUTS  AS ASCII 
      
      call rad_multi_fu  ! CALL THE CODE !!!
  
      write(*,*) cld(ii), base_alb(jj), alb(jj) + alb_interval(jj) * (kk - 1), &
          ftoa(1)%swup, ftoa(2)%swup, ftoa(3)%swup, ftoa(4)%swup
      !write(*,*) cld(ii), alb(kk), ftoa(1)%swup, ftoa(2)%swup, ftoa(3)%swup, ftoa(4)%swup
    enddo 
  enddo 
enddo

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

stop ' Simple.f90 normal end'

end
