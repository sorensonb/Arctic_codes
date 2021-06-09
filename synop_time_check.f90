subroutine synop_time_check(synop_idx, int_hr, l_in_time)
!
!  NAME:
!
!  PURPOSE:
!
!  CALLS:
!
!  MODIFICATIONS:
!
!  ###########################################################################

  implicit none

  integer                :: int_hr        ! integer variable for hour
  integer                :: synop_idx

  integer,dimension(4)   :: synop_times 
  logical                :: l_in_time
 
  synop_times = [0,6,12,18] 

  !write(*,*) "In synop_time_check"

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
