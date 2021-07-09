subroutine grid_raw_data_climo(grids,i_counts,i_size,lat_gridder,lat_thresh)
!
!  NAME:
!    grid_raw_data_climo
!
!  PURPOSE:
!    Loop over each value from the current OMI HDF5 file and
!    insert each AI value into the grid if it meets the criteria.
!
!  CALLS:
!   Modules:
!     - h5_vars (custom module, shared between count and climo analyses)
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

  real                   :: grids(360,i_size)     ! one-degree AI grid
                                                  ! values.
  integer                :: i_counts(360,i_size)  ! one-degree AI counts
  integer                :: i_size                ! lat array size
  real                   :: lat_gridder           ! lat grid value
  real                   :: lat_thresh            ! Latitude threshold value

  integer                :: ii                    ! loop counter
  integer                :: jj                    ! loop counter
  integer                :: kk                    ! loop counter    
  integer                :: index1                ! index in grid
  integer                :: index2                ! index in grid
  integer                :: i_sfc_flag            ! sfc QC flag from the GPQF
                                                  ! flag

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

        ! Insert the current AI value from the file into the running average
        ! ------------------------------------------------------------------ 
        grids(index2,index1) = ((grids(index2,index1) * &
            i_counts(index2,index1)) + AI_data(jj,ii)) / &
           (i_counts(index2,index1)+1)
        i_counts(index2,index1) = i_counts(index2,index1) + 1
      endif
    enddo row_loop
  enddo time_loop

end subroutine grid_raw_data_climo
