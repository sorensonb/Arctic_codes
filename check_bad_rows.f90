subroutine check_bad_rows(c_in_name,errout,io10)
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

  character(len = 255)   :: c_in_name   ! HDF5 file name
  integer                :: errout      ! file object for error file
  integer                :: io10        ! file object for bad row file

  character(len = 255)   :: dtg         ! dtg from each line in bad row file
  integer                :: ii          ! loop counter 
  integer                :: istatus     ! read status variable

  integer                :: i_row_year  ! integer year from row file
  integer                :: i_row_month ! integer month from row file
  integer                :: i_row_day   ! integer day from row file

  integer                :: i_file_year   ! integer year from HDF5 file
  integer                :: i_file_month  ! integer month from HDF5 file
  integer                :: i_file_day    ! integer day from HDF5 file
 
  logical                :: l_continue   ! boolean used to handle the row loop

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  l_continue = .true.

  ! Extract day information from file name
  ! --------------------------------------
  read(c_in_name(44:47), *) i_file_year
  read(c_in_name(49:50), *) i_file_month
  read(c_in_name(51:52), *) i_file_day

  ! if i_bad_list is allocated from a previous checking of the list, 
  ! then deallocate it before setting up a new list
  ! ----------------------------------------------------------------
  if(allocated(i_bad_list)) deallocate(i_bad_list)

  check_loop: do while(l_continue)
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
      ! Extract the date information from the row date
      ! ----------------------------------------------
      read(dtg(1:4), *) i_row_year
      read(dtg(5:6), *) i_row_month
      read(dtg(7:8), *) i_row_day

      ! Compare the day from the file name file to the day from the bad
      ! row file. If not the same, must cycle the check loop
      ! ---------------------------------------------------------------
      if((i_row_year /= i_file_year) .or. &
         (i_row_month /= i_file_month) .or. &
         (i_row_day /= i_file_day)) then
        ! Cycle the row loop and read the next line
        ! -----------------------------------------
        read(io10,*)
        cycle check_loop
      endif
      
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

    ! Set the switch to exit the loop
    ! -------------------------------
    l_continue = .false.
  enddo check_loop 
end subroutine check_bad_rows
