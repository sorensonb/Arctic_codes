subroutine grid_raw_data(errout,grids,i_counts,i_size,&
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

  use h5_vars, only : AI_dims, AZM_dims, GPQF_dims, LAT_dims, LON_dims, &
                      XTRACK_dims, &
                      AI_data, AZM_data, GPQF_data, LAT_data, LON_data, &
                      XTRACK_data, &
                      integer2binary

  implicit none

  integer                :: io7
  integer                :: errout
  integer                :: i_size
  real                   :: grids(1440,i_size)
  integer                :: i_counts(1440,i_size)
  real                   :: lat_gridder
  real                   :: lat_thresh

  integer                :: ii
  integer                :: jj
  integer                :: istatus
  integer                :: index1        
  integer                :: index2        
  integer,dimension(32)  :: bi

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  ! Loop over the array contents
  ! -------------------------
  row_loop: do ii=1,AI_dims(2)
    do jj=1,AI_dims(1) 

      ! = = = = = = = = = = = = = = = = = = = =
      ! Account for bad rows here???????? 
      ! Cycle loop if this index in bad rows
      ! = = = = = = = = = = = = = = = = = = = =

      ! Convert current GPQF value to binary
      if(GPQF_data(jj,ii) > -200) write(errout,*) GPQF_data(jj,ii)
      ! Then string
      ! Then extract indices
      ! Then convert it back to 

      if((XTRACK_data(jj,ii) == 0) .and. &
          (LAT_data(jj,ii) > lat_thresh) .and. &
          (AZM_data(jj,ii) > 100)) then
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
