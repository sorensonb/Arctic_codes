subroutine read_shawn_file_climo(io7,errout,c_total_file_name,grids,i_counts,i_size,&
             lat_gridder,lat_thresh)
!
!  NAME:
!    read_shawn_file_climo
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

  implicit none

  integer                :: io7
  integer                :: errout
  integer                :: i_size
  character(len = 255)   :: c_total_file_name 
  real                   :: grids(360,i_size)
  integer                :: i_counts(360,i_size)
  real                   :: lat_gridder
  real                   :: lat_thresh

  integer                :: istatus
  integer                :: index1        
  integer                :: index2        

  ! Variables from each line in Shawn's files
  real                   :: lat
  real                   :: lon
  real                   :: raw_ai
  real                   :: filter
  real                   :: clean_ai
  real                   :: v5
  real                   :: v6
  real                   :: v7
  real                   :: v8
  real                   :: v9
  real                   :: v10
  real                   :: v11
  real                   :: v12
  real                   :: v13
  real                   :: v14

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Loop over the file
  ! -------------------------
  data_loop: do
    read(io7, *, iostat = istatus)  &
            lat, lon, raw_ai, filter, clean_ai,v5,v6,v7,v8,v9,v10,&
              v11,v12,v13,v14
    if(istatus > 0) then
      write(errout, *) "ERROR: error reading data from ", &
          trim(c_total_file_name)
      write(errout, *) "       cycling data_loop"
      cycle data_loop
    else if(istatus < 0) then
      !write(errout, *) "End of data in file: ",trim(c_total_file_name)
      exit data_loop
    endif
    ! Read a line from the file

    if(lat > lat_thresh) then
      ! Average the data into the grid?
      ! -------------------------------
      index1 = floor(lat - lat_gridder) + 1
      index2 = floor(lon + 180) + 1

      if(index1 < 1) index1 = 1
      if(index1 > i_size) index1 = i_size
      if(index2 < 1) index2 = 1
      if(index2 > 360) index2 = 360

      grids(index2,index1) = ((grids(index2,index1) * &
          i_counts(index2,index1)) + clean_ai) / &
         (i_counts(index2,index1)+1)
      i_counts(index2,index1) = i_counts(index2,index1) + 1
    endif
  enddo data_loop

end subroutine read_shawn_file_climo
