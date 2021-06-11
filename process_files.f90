subroutine process_files(io8,errout,io6,data_path,grids,&
             i_counts,i_size,synop_idx,int_hr,lat_gridder,lat_thresh)
!
!  NAME:
!    process_files
!
!  PURPOSE:
!
!  CALLS:
!    None
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/06/10: Written
!
!  ###########################################################################

  implicit none

  integer                :: io8       ! File object for file name file
  integer                :: io6       ! Data output file
  integer                :: errout
  integer                :: synop_idx
  integer                :: i_size
  integer                :: int_hr        ! integer variable for hour
  character(len = 56)    :: c_total_file_name 
  character(len = 12)    :: dtg
  logical                :: l_in_time
  real                   :: grids(1440,i_size)
  integer                :: i_counts(1440,i_size)
  real                   :: lat_gridder
  real                   :: lat_thresh

  integer                :: istatus
  integer                :: index1        
  integer                :: index2        

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  ! Loop over the file
  file_loop: do
    ! Read the current dtg from the file
    read(io8, *, iostat=istatus) dtg
    if(istatus < 0) then 
      write(*,*) "End of omi_dates.txt found."
      exit
    else if(istatus > 0) then
      write(errout,*) "ERROR: problem reading dtg"
      cycle file_loop
    else

      ! Extract time information from dtg
      ! ---------------------------------
      read(dtg(9:10), *) int_hr

      ! See if the hour exceeds the current 6 hr assimilation window.
      ! If so, calculate averages and counts and reset variables.
      ! ------------------------------------------------------------
      call synop_time_check(synop_idx, int_hr, l_in_time)
      if(.not. l_in_time) then 
        call count_ai(io6,grids,i_counts,i_size,ai_thresh,synop_idx,&
                      ai_count,dtg)
        !! Loop over the grid and count up grids with high AI
        !do ii=1,i_size
        !  do jj=1,1440
        !    if(i_counts(jj,ii) > 0) then
        !      avg_ai = grids(jj,ii)/i_counts(jj,ii)
        !      if(avg_ai > ai_thresh) then
        !        ai_count = ai_count + 1
        !      endif 
        !    endif
        !  enddo  
        !enddo  

        !write(io6,*) dtg(1:8),synop_times(synop_idx), ai_count
        
        !! Reset grid arrays
        !synop_idx = synop_idx + 1
        !if(synop_idx == 5) synop_idx = 1
        !ai_count = 0     
        !grids(:,:) = 0.
        !i_counts(:,:) = 0

      endif  

      write(*,*) data_path//dtg

      ! Open the shawn file and look at contents
      open(io7, file = data_path//dtg, iostat = istatus)
      if(istatus /= 0) then
        write(errout, *) "ERROR: error opening file",data_path//dtg
        write(errout, *) "       cycling file_loop"
        cycle file_loop
      endif

      call process_files(io7,errout,data_path//dtg,grids,i_counts,&
                           i_size,lat_gridder,lat_thresh)
      ! Loop over the file
      ! -------------------------
      !data_loop: do
      !  read(io7, *, iostat = istatus)  &
      !          lat, lon, raw_ai, filter, clean_ai,v5,v6,v7,v8,v9,v10,&
      !            v11,v12,v13,v14
      !  if(istatus > 0) then
      !    write(errout, *) "ERROR: error reading data from ",data_path//dtg
      !    write(errout, *) "       cycling data_loop"
      !    cycle data_loop
      !  else if(istatus < 0) then
      !    write(errout, *) "End of data in file: ",data_path//dtg
      !    exit data_loop
      !  endif
      !  ! Read a line from the file

      !  if(lat > lat_thresh) then
      !    ! Average the data into the grid?
      !    ! -------------------------------
      !    index1 = floor(lat*4 - lat_gridder)
      !    index2 = floor(lon*4 + 720)

      !    if(index1 < 1) index1 = 1
      !    if(index1 > i_size) index1 = i_size
      !    if(index2 < 1) index2 = 1
      !    if(index2 > 1440) index2 = 1440

      !    grids(index2,index1) = ((grids(index2,index1) * &
      !        i_counts(index2,index1)) + clean_ai) / &
      !       (i_counts(index2,index1)+1)
      !    i_counts(index2,index1) = i_counts(index2,index1) + 1
      !  endif
      !enddo data_loop
        
      close(io7)
    endif
  enddo file_loop  

end subroutine process_files
