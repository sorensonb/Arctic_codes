subroutine check_bad_rows(errout,io10)
!
!  NAME:
!    check_bad_rows
!
!  PURPOSE:
!    Read a line from the OMI row anomaly bad row database and insert
!    any bad rows into the i_bad_list array. 
! 
!    NOTE: This code assumes that the first date in the bad row database
!          corresponds to the first date being analyzed. This could be coded
!          up better, but I am choosing not to right now.   
! 
!  CALLS:
!    None
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/06/10: Written
!
!  ###########################################################################

  use h5_vars, only : i_num_bad, i_bad_list

  implicit none

  integer                :: errout   ! file object for error file
  integer                :: io10     ! file object for bad row file

  character(len = 255)   :: dtg      ! dtg from each line in bad row file
  integer                :: ii       ! loop counter 
  integer                :: istatus  ! read status variable

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! if i_bad_list is allocated from a previous checking of the list, 
  ! then deallocate it before setting up a new list
  ! ----------------------------------------------------------------
  if(allocated(i_bad_list)) deallocate(i_bad_list)
 
  ! Read the current date and row count from the file
  ! -------------------------------------------------
  read(io10, '(a8,1x,i1)' , iostat=istatus,advance='no') &
      dtg,i_num_bad
  if(istatus < 0) then 
    write(*,*) "End of row file found"
    return
  else if(istatus > 0) then
    write(errout,*) "ERROR: problem reading line from bad row file"
    return
  else
    if(i_num_bad > 0) then
      allocate(i_bad_list(i_num_bad))

      ! Read each bad row number into the bad list array
      ! ------------------------------------------------
      do ii=1,i_num_bad
        read(io10,'(1x,i2)',iostat=istatus,advance='no') i_bad_list(ii)
      enddo
    endif

    ! Advance to the next line
    ! ------------------------
    read(io10,*)
  endif

end subroutine check_bad_rows
