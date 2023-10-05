subroutine bin_OMI_GPQF_values(ii,jj,kk)
!
! NAME:
!   bin_OMI_GPQF_values.f90
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

  use daily_vars, only: &
    count_sfland, count_seaice, count_permice, count_drysnow, count_ocean, &
    count_mixpixel, count_other, grid_GPQF

  implicit none

  integer               :: ii
  integer               :: jj
  integer               :: kk

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Value 1 - sfland (snow-free land)
  ! --------------------------------------------------
  if(  ((count_sfland(ii,jj,kk) > count_seaice(ii,jj,kk)) .and. &
        (count_sfland(ii,jj,kk) > count_permice(ii,jj,kk)) .and. &
        (count_sfland(ii,jj,kk) > count_drysnow(ii,jj,kk)) .and. &
        (count_sfland(ii,jj,kk) > count_ocean(ii,jj,kk)) .and. &
        (count_sfland(ii,jj,kk) > count_mixpixel(ii,jj,kk)) .and. &
        (count_sfland(ii,jj,kk) > count_other(ii,jj,kk)))) then
    grid_GPQF(ii,jj,kk) = 1
  ! Value 2 - seaice
  ! --------------------------------------------------
  else if(  ((count_seaice(ii,jj,kk) > count_sfland(ii,jj,kk)) .and. &
             (count_seaice(ii,jj,kk) > count_permice(ii,jj,kk)) .and. &
             (count_seaice(ii,jj,kk) > count_drysnow(ii,jj,kk)) .and. &
             (count_seaice(ii,jj,kk) > count_ocean(ii,jj,kk)) .and. &
             (count_seaice(ii,jj,kk) > count_mixpixel(ii,jj,kk)) .and. &
             (count_seaice(ii,jj,kk) > count_other(ii,jj,kk)))) then
    grid_GPQF(ii,jj,kk) = 2
  ! Value 3 - permanent ice (Greenland, Antarctica)
  ! ----------------------------------------------------
  else if(  ((count_permice(ii,jj,kk) > count_sfland(ii,jj,kk)) .and. &
             (count_permice(ii,jj,kk) > count_seaice(ii,jj,kk)) .and. &
             (count_permice(ii,jj,kk) > count_drysnow(ii,jj,kk)) .and. &
             (count_permice(ii,jj,kk) > count_ocean(ii,jj,kk)) .and. &
             (count_permice(ii,jj,kk) > count_mixpixel(ii,jj,kk)) .and. &
             (count_permice(ii,jj,kk) > count_other(ii,jj,kk)))) then
    grid_GPQF(ii,jj,kk) = 3
  ! Value 4 - dry snow (seasonal snow)
  ! ----------------------------------------------------
  else if(  ((count_drysnow(ii,jj,kk) > count_sfland(ii,jj,kk)) .and. &
             (count_drysnow(ii,jj,kk) > count_seaice(ii,jj,kk)) .and. &
             (count_drysnow(ii,jj,kk) > count_permice(ii,jj,kk)) .and. &
             (count_drysnow(ii,jj,kk) > count_ocean(ii,jj,kk)) .and. &
             (count_drysnow(ii,jj,kk) > count_mixpixel(ii,jj,kk)) .and. &
             (count_drysnow(ii,jj,kk) > count_other(ii,jj,kk)))) then
    grid_GPQF(ii,jj,kk) = 4
  ! Value 5 - ocean
  ! ----------------------------------------------------
  else if(  ((count_ocean(ii,jj,kk) > count_sfland(ii,jj,kk)) .and. &
             (count_ocean(ii,jj,kk) > count_seaice(ii,jj,kk)) .and. &
             (count_ocean(ii,jj,kk) > count_permice(ii,jj,kk)) .and. &
             (count_ocean(ii,jj,kk) > count_drysnow(ii,jj,kk)) .and. &
             (count_ocean(ii,jj,kk) > count_mixpixel(ii,jj,kk)) .and. &
             (count_ocean(ii,jj,kk) > count_other(ii,jj,kk)))) then
    grid_GPQF(ii,jj,kk) = 5
  ! Value 6 - mixed pixel
  ! ----------------------------------------------------
  else if(  ((count_mixpixel(ii,jj,kk) > count_sfland(ii,jj,kk)) .and. &
             (count_mixpixel(ii,jj,kk) > count_seaice(ii,jj,kk)) .and. &
             (count_mixpixel(ii,jj,kk) > count_permice(ii,jj,kk)) .and. &
             (count_mixpixel(ii,jj,kk) > count_drysnow(ii,jj,kk)) .and. &
             (count_mixpixel(ii,jj,kk) > count_ocean(ii,jj,kk)) .and. &
             (count_mixpixel(ii,jj,kk) > count_other(ii,jj,kk)))) then
    grid_GPQF(ii,jj,kk) = 6
  ! Value 7 - other
  ! ----------------------------------------------------
  else if(  ((count_other(ii,jj,kk) > count_sfland(ii,jj,kk)) .and. &
             (count_other(ii,jj,kk) > count_seaice(ii,jj,kk)) .and. &
             (count_other(ii,jj,kk) > count_permice(ii,jj,kk)) .and. &
             (count_other(ii,jj,kk) > count_drysnow(ii,jj,kk)) .and. &
             (count_other(ii,jj,kk) > count_ocean(ii,jj,kk)) .and. &
             (count_other(ii,jj,kk) > count_mixpixel(ii,jj,kk)))) then
    grid_GPQF(ii,jj,kk) = 7
  else
  ! Value 0 - none of the previous conditions met. Likely caused by
  ! grid boxes with the same counts across types
  ! ---------------------------------------------------------------
    grid_GPQF(ii,jj,kk) = 0
  endif
  
end subroutine bin_OMI_GPQF_values
