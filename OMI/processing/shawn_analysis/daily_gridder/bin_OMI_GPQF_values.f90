subroutine bin_OMI_GPQF_values(kk,jj,ii)
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

  integer               :: kk
  integer               :: jj
  integer               :: ii

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Value 1 - sfland (snow-free land)
  ! --------------------------------------------------
  if(  ((count_sfland(kk,jj,ii) > count_seaice(kk,jj,ii)) .and. &
        (count_sfland(kk,jj,ii) > count_permice(kk,jj,ii)) .and. &
        (count_sfland(kk,jj,ii) > count_drysnow(kk,jj,ii)) .and. &
        (count_sfland(kk,jj,ii) > count_ocean(kk,jj,ii)) .and. &
        (count_sfland(kk,jj,ii) > count_mixpixel(kk,jj,ii)) .and. &
        (count_sfland(kk,jj,ii) > count_other(kk,jj,ii)))) then
    grid_GPQF(kk,jj,ii) = 1
  ! Value 2 - seaice
  ! --------------------------------------------------
  else if(  ((count_seaice(kk,jj,ii) > count_sfland(kk,jj,ii)) .and. &
             (count_seaice(kk,jj,ii) > count_permice(kk,jj,ii)) .and. &
             (count_seaice(kk,jj,ii) > count_drysnow(kk,jj,ii)) .and. &
             (count_seaice(kk,jj,ii) > count_ocean(kk,jj,ii)) .and. &
             (count_seaice(kk,jj,ii) > count_mixpixel(kk,jj,ii)) .and. &
             (count_seaice(kk,jj,ii) > count_other(kk,jj,ii)))) then
    grid_GPQF(kk,jj,ii) = 2
  ! Value 3 - permanent ice (Greenland, Antarctica)
  ! ----------------------------------------------------
  else if(  ((count_permice(kk,jj,ii) > count_sfland(kk,jj,ii)) .and. &
             (count_permice(kk,jj,ii) > count_seaice(kk,jj,ii)) .and. &
             (count_permice(kk,jj,ii) > count_drysnow(kk,jj,ii)) .and. &
             (count_permice(kk,jj,ii) > count_ocean(kk,jj,ii)) .and. &
             (count_permice(kk,jj,ii) > count_mixpixel(kk,jj,ii)) .and. &
             (count_permice(kk,jj,ii) > count_other(kk,jj,ii)))) then
    grid_GPQF(kk,jj,ii) = 3
  ! Value 4 - dry snow (seasonal snow)
  ! ----------------------------------------------------
  else if(  ((count_drysnow(kk,jj,ii) > count_sfland(kk,jj,ii)) .and. &
             (count_drysnow(kk,jj,ii) > count_seaice(kk,jj,ii)) .and. &
             (count_drysnow(kk,jj,ii) > count_permice(kk,jj,ii)) .and. &
             (count_drysnow(kk,jj,ii) > count_ocean(kk,jj,ii)) .and. &
             (count_drysnow(kk,jj,ii) > count_mixpixel(kk,jj,ii)) .and. &
             (count_drysnow(kk,jj,ii) > count_other(kk,jj,ii)))) then
    grid_GPQF(kk,jj,ii) = 4
  ! Value 5 - ocean
  ! ----------------------------------------------------
  else if(  ((count_ocean(kk,jj,ii) > count_sfland(kk,jj,ii)) .and. &
             (count_ocean(kk,jj,ii) > count_seaice(kk,jj,ii)) .and. &
             (count_ocean(kk,jj,ii) > count_permice(kk,jj,ii)) .and. &
             (count_ocean(kk,jj,ii) > count_drysnow(kk,jj,ii)) .and. &
             (count_ocean(kk,jj,ii) > count_mixpixel(kk,jj,ii)) .and. &
             (count_ocean(kk,jj,ii) > count_other(kk,jj,ii)))) then
    grid_GPQF(kk,jj,ii) = 5
  ! Value 6 - mixed pixel
  ! ----------------------------------------------------
  else if(  ((count_mixpixel(kk,jj,ii) > count_sfland(kk,jj,ii)) .and. &
             (count_mixpixel(kk,jj,ii) > count_seaice(kk,jj,ii)) .and. &
             (count_mixpixel(kk,jj,ii) > count_permice(kk,jj,ii)) .and. &
             (count_mixpixel(kk,jj,ii) > count_drysnow(kk,jj,ii)) .and. &
             (count_mixpixel(kk,jj,ii) > count_ocean(kk,jj,ii)) .and. &
             (count_mixpixel(kk,jj,ii) > count_other(kk,jj,ii)))) then
    grid_GPQF(kk,jj,ii) = 6
  ! Value 7 - other
  ! ----------------------------------------------------
  else if(  ((count_other(kk,jj,ii) > count_sfland(kk,jj,ii)) .and. &
             (count_other(kk,jj,ii) > count_seaice(kk,jj,ii)) .and. &
             (count_other(kk,jj,ii) > count_permice(kk,jj,ii)) .and. &
             (count_other(kk,jj,ii) > count_drysnow(kk,jj,ii)) .and. &
             (count_other(kk,jj,ii) > count_ocean(kk,jj,ii)) .and. &
             (count_other(kk,jj,ii) > count_mixpixel(kk,jj,ii)))) then
    grid_GPQF(kk,jj,ii) = 7
  else
  ! Value 0 - none of the previous conditions met. Likely caused by
  ! grid boxes with the same counts across types
  ! ---------------------------------------------------------------
    grid_GPQF(kk,jj,ii) = 0
  endif
  
end subroutine bin_OMI_GPQF_values
