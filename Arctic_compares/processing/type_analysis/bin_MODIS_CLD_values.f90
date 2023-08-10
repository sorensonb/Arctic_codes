subroutine bin_MODIS_CLD_values(ii,jj,kk)
!
! NAME:
!   bin_MODIS_CLD_values.f90
!
! PURPOSE:
!   Counts the number of individual swaths in the comp_grid_climo DTG file
!     AND the number of individual days in the file. 
! 
! CALLS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/08/08:
!     Written
!
!  ############################################################################

  use type_vars, only: &
    count_MODIS_CLD_0, &
    count_MODIS_CLD_1, &
    count_MODIS_CLD_2, &
    count_MODIS_CLD_3, &
    grid_MODIS_CLD

  implicit none

  integer               :: ii
  integer               :: jj
  integer               :: kk

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Examples of the cases for below:
  !   [xx,yy,zz,aa]
  !     - xx = counts of "cloudy"
  !     - yy = counts of "probably cloudy" 
  !     - zz = counts of "probably clear"
  !     - aa = counts of "clear"
  !
  ! "Cloudy" (out cloud value = 0.)
  !   - [2,0,0,0]
  !   - [2,0,0,1]
  !   - [2,0,1,1]
  !   - [2,1,0,0]
  !   - [2,1,1,0]
  !   - [2,1,0,1]
  ! "Probably Cloudy" (out cloud value = 1.)
  !   - [0,1,0,0]
  !   - [0,2,2,0]
  !   - [0,2,0,1]
  !   - [0,2,1,0]
  !   - [2,2,0,0]
  !
  ! "Probably Clear" (out cloud value = 2.)
  !   - [0,0,2,0]
  !   - [0,0,1,1]
  !   - [0,2,2,2]
  !   - [0,1,2,0]
  !   - [0,1,2,1]
  !   - [0,1,2,2]
  !
  ! "Clear" (out cloud value = 3.)
  !   - [0,0,0,2]
  !   - [0,0,1,2]
  !   - [0,1,0,2]
  !   - [0,1,1,2]
  !   - [1,0,0,1]
  !   - [1,0,1,1]
  !   - [2,2,2,2]
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

  ! "Clear" final pixels are only ones in which the 
  ! majority of the pixels in this OMI pixel are clear
  ! --------------------------------------------------
  if(  ((count_MODIS_CLD_3(ii,jj,kk) > count_MODIS_CLD_0(ii,jj,kk)) .and. &
        (count_MODIS_CLD_3(ii,jj,kk) > count_MODIS_CLD_1(ii,jj,kk)) .and. &
        (count_MODIS_CLD_3(ii,jj,kk) > count_MODIS_CLD_2(ii,jj,kk)))) then
    grid_MODIS_CLD(ii,jj,kk) = 3.
  ! "Probably clear" final pixels are only ones in which the 
  ! majority of the pixels in this OMI pixel are probably clear
  ! OR 
  ! the number of "probably clear" pixels is the same as the
  ! number of "clear" and "probably cloudy" pixels
  ! --------------------------------------------------
  else if(((count_MODIS_CLD_2(ii,jj,kk) > count_MODIS_CLD_3(ii,jj,kk)) .and. &
           (count_MODIS_CLD_2(ii,jj,kk) > count_MODIS_CLD_1(ii,jj,kk)) .and. &
           (count_MODIS_CLD_2(ii,jj,kk) > count_MODIS_CLD_0(ii,jj,kk))) .or. &
          ((count_MODIS_CLD_3(ii,jj,kk) > count_MODIS_CLD_0(ii,jj,kk)) .and. &
           (count_MODIS_CLD_2(ii,jj,kk) >= count_MODIS_CLD_3(ii,jj,kk)) .and. &
           (count_MODIS_CLD_2(ii,jj,kk) >= count_MODIS_CLD_1(ii,jj,kk))) .or. &
          ((count_MODIS_CLD_3(ii,jj,kk) == count_MODIS_CLD_0(ii,jj,kk)) .and. &
           (count_MODIS_CLD_2(ii,jj,kk) > count_MODIS_CLD_3(ii,jj,kk)) .and. &
           (count_MODIS_CLD_2(ii,jj,kk) == count_MODIS_CLD_1(ii,jj,kk)))) then
    grid_MODIS_CLD(ii,jj,kk) = 2.
  ! If there are the same number of "clear" and "cloudy"
  ! pixels, classify final pixel as "cloudy"
  ! ----------------------------------------------------
  else if(((count_MODIS_CLD_1(ii,jj,kk) > count_MODIS_CLD_1(ii,jj,kk)) .and. &
           (count_MODIS_CLD_1(ii,jj,kk) > count_MODIS_CLD_2(ii,jj,kk)) .and. &
           (count_MODIS_CLD_1(ii,jj,kk) > count_MODIS_CLD_0(ii,jj,kk))) .or. &
          ((count_MODIS_CLD_0(ii,jj,kk) > count_MODIS_CLD_3(ii,jj,kk)) .and. &
           (count_MODIS_CLD_1(ii,jj,kk) >= count_MODIS_CLD_0(ii,jj,kk)) .and. &
           (count_MODIS_CLD_2(ii,jj,kk) <= count_MODIS_CLD_1(ii,jj,kk))) .or. &
          ((count_MODIS_CLD_3(ii,jj,kk) == count_MODIS_CLD_1(ii,jj,kk)) .and. &
           (count_MODIS_CLD_3(ii,jj,kk) > count_MODIS_CLD_0(ii,jj,kk)) .and. &
           (count_MODIS_CLD_1(ii,jj,kk) > count_MODIS_CLD_2(ii,jj,kk)))) then
    grid_MODIS_CLD(ii,jj,kk) = 1.
  else if((count_MODIS_CLD_0(ii,jj,kk) >= count_MODIS_CLD_3(ii,jj,kk)) .and. &
          (count_MODIS_CLD_0(ii,jj,kk) >= count_MODIS_CLD_2(ii,jj,kk)) .and. &
          (count_MODIS_CLD_0(ii,jj,kk) >= count_MODIS_CLD_1(ii,jj,kk))) then
    grid_MODIS_CLD(ii,jj,kk) = 0.
  endif
  
end subroutine bin_MODIS_CLD_values
