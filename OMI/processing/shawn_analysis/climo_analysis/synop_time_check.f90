subroutine synop_time_check(synop_idx, int_hr, l_in_time)
!
!  NAME:
!    synop_time_check
!
!  PURPOSE:
!    Compare the hour from the current file timestamp to the current
!    analyzed synoptic time and see if the file hour falls within the
!    +/- 3 hr synoptic time window. The l_in_time logical variable is 
!    set to 'true' if the hour falls within the time window and 'false'
!    if the hour falls outside the time window.
!
!  CALLS:
!    None
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/06/10: Written
!
!  ###########################################################################

  implicit none

  integer                :: int_hr        ! integer variable for hour
  integer                :: synop_idx

  integer,dimension(4)   :: synop_times 
  logical                :: l_in_time

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  synop_times = [0,6,12,18] 

  if(synop_idx == 1) then
    if((int_hr < 3) .or. &
       (int_hr > 21)) then
      l_in_time = .true.
    else
      l_in_time = .false.
    endif
     
  else
    if((int_hr > (synop_times(synop_idx) - 3)) .and. &
       (int_hr < (synop_times(synop_idx) + 3))) then
      l_in_time = .true.
    else
      l_in_time = .false.
    endif
  endif

end subroutine synop_time_check
