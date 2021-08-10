subroutine grid_raw_data_drift(avg_ai,i_count)
!
!  NAME:
!    grid_raw_data_drift
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
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/16: Written
!
!  ###########################################################################

  use h5_vars, only : AI_dims, AZM_dims, GPQF_dims, LAT_dims, LON_dims, &
                      XTRACK_dims, &
                      AI_data, AZM_data, GPQF_data, LAT_data, LON_data, &
                      XTRACK_data, &
                      get_ice_flags,i_num_bad,i_bad_list

  implicit none

  real                   :: avg_ai                ! avg AI value
  integer                :: i_count               ! AI counts

  integer                :: ii                    ! loop counter
  integer                :: jj                    ! loop counter
  integer                :: kk                    ! loop counter    
  integer                :: i_sfc_flag            ! sfc QC flag from the GPQF
                                                  ! flag

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  ! Loop over the array contents
  ! -------------------------
  time_loop: do ii=1,AI_dims(2)
    ! For BS2, use only rows 1 through 21, which should be unaffected through
    ! the entire dataset.
    ! -----------------------------------------------------------------------
    !row_loop: do jj=1,21 
    row_loop: do jj=1,AI_dims(1) 

      ! Account for bad rows here
      ! Cycle loop if this index in bad rows
      ! ------------------------------------
      if(allocated(i_bad_list)) then
        do kk=1,i_num_bad
          if(jj == i_bad_list(kk)) cycle row_loop
        enddo
      endif

      ! 
      if((jj > 20)) then
      !!#!if((jj == 49) .or. &
      !!#!   (jj == 50) .or. &
      !!#!   !(jj == 53) .or. &
      !!#!   (jj >= 56)) then 

        ! Use the get_ice_flags function from h5_vars to extract the sfc type
        ! flag from the ground pixel quality flag
        ! -------------------------------------------------------------------
        i_sfc_flag = get_ice_flags(GPQF_data(jj,ii))

        if((XTRACK_data(jj,ii) == 0) .and. &
            ((LAT_data(jj,ii) > -40.) .and. (LAT_data(jj,ii) <= 0.)) .and. &
            ((LON_data(jj,ii) > -180.) .and. (LON_data(jj,ii) <= -140.)) .and. &
            ! BS2: use rows 1-21, so must comment out AZM check
            ! -------------------------------------------------
            (AZM_data(jj,ii) > 100) .and. &
            (AI_data(jj,ii) > -2e5) .and. &
            ! VJZ2: no snow-free land either
            ( (i_sfc_flag >= 1 .and. i_sfc_flag <= 101) .or. &
              (i_sfc_flag == 104) )) then

          ! Insert the current AI value from the file into the running average
          ! ------------------------------------------------------------------ 
          avg_ai = ((avg_ai * i_count) + AI_data(jj,ii)) / (i_count+1)
          i_count = i_count + 1
        endif
      endif
    enddo row_loop
  enddo time_loop

end subroutine grid_raw_data_drift
