subroutine check_row_xtrack(row_avgs,lat_thresh)
!
!  NAME:
!    check_row_xtrack
!
!  PURPOSE:
!    Loop over each value from the current OMI HDF5 file and
!    insert each AI value into the grid if it meets the criteria.
!
!  CALLS:
!   Modules:
!     - h5_vars (custom module, shared between count and climo analyses)
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/06/10: Written
!
!  ###########################################################################

  use h5_vars, only : AI_dims, LAT_dims, XTRACK_dims, &
                      AI_data, LAT_data, XTRACK_data

  implicit none

  integer                :: row_avgs(60)     ! one-degree AI grid
                                             ! values.
  real                   :: lat_thresh
  integer                :: ii                    ! loop counter
  integer                :: jj                    ! loop counter

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  ! Loop over the array contents
  ! -------------------------
  time_loop: do ii=1,AI_dims(2)
    row_loop: do jj=1,AI_dims(1) 

     if((LAT_data(jj,ii) >= lat_thresh) .and. &
        ((XTRACK_data(jj,ii) /= 0) .or. &
        (abs(AI_data(jj,ii)) <= -2e5))) then

        row_avgs(jj) = 1
     endif
    enddo row_loop
  enddo time_loop

end subroutine check_row_xtrack
