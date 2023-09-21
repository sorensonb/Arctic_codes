subroutine count_num_days_shawn(file_name_file, i_num_swath, i_num_days)
!
! NAME:
!   count_num_days_shawn.f90
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

  implicit none

  character(len = 255)   :: file_name_file  ! name of file containing the 
                                            ! list of shawn file names to 
                                            ! be analyzed
  integer                :: i_num_swath
  integer                :: i_num_days
  integer                :: i_int_day
  integer                :: i_current_day

  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer                :: istatus

  character(len = 255)   :: single_file_date        ! dtg from each line of file

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Loop over the file, counting both the total number of swaths AND
  ! the number of total days in the list of swath DTGs
  ! ----------------------------------------------------------------

  i_num_swath = 0
  i_num_days = 0
  i_current_day = 0

  ! Open the file name file
  ! -----------------------
  open(io8, file = trim(file_name_file), iostat = istatus)
  if(istatus > 0) then
    write(*,*) "ERROR: Problem reading "//trim(file_name_file)
    return 
  else
    ! Loop over the file name file
    ! ----------------------------
    file_loop: do
      ! Read the current total_file_name from the file
      ! ----------------------------------------------
      read(io8, '(A)', iostat=istatus) single_file_date
      if(istatus < 0) then 
        write(*,*) "End of "//trim(file_name_file)//" found"
        exit
      else if(istatus > 0) then
        write(*,*) "ERROR: problem reading total_file_name"
        cycle file_loop
      else
        i_num_swath = i_num_swath + 1
        read(single_file_date(7:8), '(i2)') i_int_day

        if( (i_int_day /= i_current_day)) then
          i_num_days = i_num_days + 1
        endif 

        i_current_day = i_int_day

      endif

    enddo file_loop  
  endif
 
  close(io8)
 
end subroutine count_num_days_shawn
