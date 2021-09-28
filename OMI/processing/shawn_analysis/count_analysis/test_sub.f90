subroutine test_sub()
!
!  NAME:
!    test_sub.f90
!
!  PURPOSE:
!    Calculate Q_ext, Q_sca, Q_abs, and g 
!
!  CALLS:
!    bessel.f90, newmann.f90, calc_deriv_A.f90, calc_deriv_B.f90
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>    - 2018/10/23:
!      Written
!
!  ###########################################################################

  implicit none

  integer                :: int_hr        ! integer variable for hour
  character(len = 12)    :: dtg
  integer,dimension(4)   :: synop_times 
  integer                :: synop_idx
  logical                :: l_in_time
 
  dtg = "200804220945" 
  read(dtg(9:10), *) int_hr
  synop_times = [0,6,12,18] 
  synop_idx = 2

  write(*,*) "In test_sub"

  write(*,*) int_hr,synop_times(synop_idx)
  
  call synop_time_check(synop_idx, int_hr, l_in_time)
  write(*,*) "After routine, l_in_time = ",l_in_time
  if(.not.l_in_time) then
    write(*,*) "Outside the time. Averaging"
  endif


end subroutine test_sub
