subroutine calc_grid_avgs(i_num_days, i_lat_size, i_lon_size)
!
! NAME:
!   write_output_file.f90
!
! PURPOSE:
! 
! CALLS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/08/08:
!     Written
!
!  ############################################################################

  use daily_vars, only: grid_AI, count_AI, day_values, lat_values, &
        lon_values, clear_daily_arrays

  implicit none

  integer                :: ii
  integer                :: jj
  integer                :: kk
  integer                :: i_num_days
  integer                :: i_lat_size
  integer                :: i_lon_size

  ! File read variables
  !!#!integer                :: error

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Calculate the average of each grid box in the grid 
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  day_loop: do ii = 1, i_num_days
    lat_loop: do jj = 1, i_lat_size
      lon_loop: do kk = 1, i_lon_size
        if(count_AI(kk,jj,ii) > 0) then
          grid_AI(kk,jj,ii) = grid_AI(kk,jj,ii) / count_AI(kk,jj,ii)
        else
          grid_AI(kk,jj,ii) = -9.
        endif
      enddo lon_loop
    enddo lat_loop
  enddo day_loop

  
end subroutine calc_grid_avgs
