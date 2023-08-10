subroutine bin_NSIDC_ICE_values(ii,jj,kk)
!
! NAME:
!   bin_NSIDC_ICE_values.f90
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
    count_NSIDC_ICE_0100, &
    count_NSIDC_ICE_251, &
    count_NSIDC_ICE_252, &
    count_NSIDC_ICE_253, &
    count_NSIDC_ICE_254, &
    grid_NSIDC_ICE

  implicit none

  integer               :: ii
  integer               :: jj
  integer               :: kk

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Sea ice/ocean: only OMI grid boxes that contain only sea ice/ocean
  !                No pole hole, no coast, no land.
  if(  (count_NSIDC_ICE_0100(ii,jj,kk) /= 0) .and. &
       (count_NSIDC_ICE_251(ii,jj,kk) == 0) .and. &
       (count_NSIDC_ICE_252(ii,jj,kk) == 0) .and. &
       (count_NSIDC_ICE_253(ii,jj,kk) == 0) .and. &
       (count_NSIDC_ICE_254(ii,jj,kk) == 0)) then
    grid_NSIDC_ICE(ii,jj,kk) = &
      grid_NSIDC_ICE(ii,jj,kk) / count_NSIDC_ICE_0100(ii,jj,kk)
  ! "Pole hole" : if the OMI pixel contains any 251 values, classify the
  !               final pixel as pole hole. 
  ! NOTE: SHOULD NOT HAVE ANY PIXELS THAT CONTAIN POLE HOLE AND COAST/LAND.
  else if(count_NSIDC_ICE_251(ii,jj,kk) /= 0) then
    grid_NSIDC_ICE(ii,jj,kk) = 251.

  ! "Unused data". If any of these come up, just classify as 252
  else if(count_NSIDC_ICE_252(ii,jj,kk) /= 0) then
    grid_NSIDC_ICE(ii,jj,kk) = 252.

  ! "Coastline": As above: if the current pixel has more coastline than
  !              land, classify as land. Otherwise, classify as land
  else if(count_NSIDC_ICE_253(ii,jj,kk) >= count_NSIDC_ICE_254(ii,jj,kk)) then
    grid_NSIDC_ICE(ii,jj,kk) = 253.

  ! "Land": As above: if any coastline pixels here, class as 253
  else if(count_NSIDC_ICE_254(ii,jj,kk) > count_NSIDC_ICE_253(ii,jj,kk)) then
    grid_NSIDC_ICE(ii,jj,kk) = 254.

  ! ERROR CODE. Should not come here.
  else 
    grid_NSIDC_ICE(ii,jj,kk) = 888.
  endif
  
end subroutine bin_NSIDC_ICE_values
