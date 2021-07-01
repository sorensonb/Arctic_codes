subroutine grid_raw_data(errout,c_total_file_name,grids,i_counts,i_size,&
             lat_gridder,lat_thresh)
!
!  NAME:
!    grid_raw_data
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

  use h5_vars, only : AI_dims, LAT_dims, LON_dims, XTRACK_dims, &
                      AI_data, LAT_data, LON_data, XTRACK_data

  implicit none

  integer                :: io7
  integer                :: errout
  integer                :: i_size
  character(len = 255)   :: c_total_file_name 
  real                   :: grids(1440,i_size)
  integer                :: i_counts(1440,i_size)
  real                   :: lat_gridder
  real                   :: lat_thresh

  integer                :: ii
  integer                :: jj
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
 
  ! Loop over the array contents
  ! -------------------------
  row_loop: do ii=1,AI_dims(2)
    do jj=1,AI_dims(1) 

      ! = = = = = = = = = = = = = = = = = = = =
      ! Account for bad rows here???????? 
      ! Cycle loop if this index in bad rows
      ! = = = = = = = = = = = = = = = = = = = =

      if((XTRACK_data(jj,ii) == 0) .and. &
          (LAT_data(jj,ii) > lat_thresh)) then
        ! Average the data into the grid
        ! -------------------------------
        index1 = floor(LAT_data(jj,ii)*4 - lat_gridder)
        index2 = floor(LON_data(jj,ii)*4 + 720)

        if(index1 < 1) index1 = 1
        if(index1 > i_size) index1 = i_size
        if(index2 < 1) index2 = 1
        if(index2 > 1440) index2 = 1440

        grids(index2,index1) = ((grids(index2,index1) * &
            i_counts(index2,index1)) + AI_data(jj,ii)) / &
           (i_counts(index2,index1)+1)
        i_counts(index2,index1) = i_counts(index2,index1) + 1
      endif
    enddo
  enddo row_loop

end subroutine grid_raw_data
