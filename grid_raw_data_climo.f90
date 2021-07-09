subroutine grid_raw_data_climo(errout,grids,i_counts,i_size,&
             lat_gridder,lat_thresh)
!
!  NAME:
!    grid_raw_data_climo
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
                      get_ice_flags,i_num_bad,i_bad_list

  implicit none

  integer                :: io7
  integer                :: errout
  integer                :: i_size
  real                   :: grids(360,i_size)
  integer                :: i_counts(360,i_size)
  real                   :: lat_gridder
  real                   :: lat_thresh

  integer                :: ii
  integer                :: jj
  integer                :: kk
  integer                :: istatus
  integer                :: index1        
  integer                :: index2        

  integer                :: i_sfc_flag

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  ! Loop over the array contents
  ! -------------------------
  time_loop: do ii=1,AI_dims(2)
    row_loop: do jj=1,AI_dims(1) 

      ! Account for bad rows here
      ! Cycle loop if this index in bad rows
      ! ------------------------------------
      if(allocated(i_bad_list)) then
        do kk=1,i_num_bad
          if(jj == i_bad_list(kk)) cycle row_loop
        enddo
      endif

      ! Use the get_ice_flags function from h5_vars to extract the sfc type
      ! flag from the ground pixel quality flag
      ! -------------------------------------------------------------------
      i_sfc_flag = get_ice_flags(GPQF_data(jj,ii))

      if((XTRACK_data(jj,ii) == 0) .and. &
          (LAT_data(jj,ii) > lat_thresh) .and. &
          (AZM_data(jj,ii) > 100) .and. &
          (AI_data(jj,ii) > -2e5) .and. &
          ! VJZ2: no snow-free land either
          ( (i_sfc_flag >= 1 .and. i_sfc_flag <= 101) .or. &
            (i_sfc_flag == 104) )) then
        ! Average the data into the grid
        ! -------------------------------
        index1 = floor(LAT_data(jj,ii) - lat_gridder) + 1
        index2 = floor(LON_data(jj,ii) + 180) + 1

        if(index1 < 1) index1 = 1
        if(index1 > i_size) index1 = i_size
        if(index2 < 1) index2 = 1
        if(index2 > 360) index2 = 360

        grids(index2,index1) = ((grids(index2,index1) * &
            i_counts(index2,index1)) + AI_data(jj,ii)) / &
           (i_counts(index2,index1)+1)
        i_counts(index2,index1) = i_counts(index2,index1) + 1
      endif
    enddo row_loop
  enddo time_loop

end subroutine grid_raw_data_climo
