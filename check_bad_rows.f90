subroutine check_bad_rows(errout,io10)
!
!  NAME:
!    check_bad_rows
!
!  PURPOSE:
!    Loop over each line from the current OMI AI perturbation text file and
!    insert each AI perturbation value into the grid if it meets the criteria.
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

  integer                :: errout
  integer                :: io10

  character(len = 255)   :: dtg
  integer                :: ii
  integer                :: istatus

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
    write(errout,*) "ERROR: problem reading line"
    return
  else
    if(i_num_bad > 0) then
      allocate(i_bad_list(i_num_bad))

      ! Read each bad row number
      ! ------------------------
      !write(*,*) trim(dtg),' has ',i_num_row, 'bad rows'
      do ii=1,i_num_bad
        read(io10,'(1x,i2)',iostat=istatus,advance='no') i_bad_list(ii)

        !write(*,*) '  ',temp_bad_row
      enddo
    endif

    ! Advance to the next line
    ! ------------------------
    read(io10,*)
  endif

end subroutine check_bad_rows
