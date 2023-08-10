subroutine average_final_grid_data(i_num_days, i_lat_size, i_lon_size)
!
! NAME:
!   average_final_grid_data.f90
!
! PURPOSE:
!   Counts the number of individual swaths in the comp_grid_climo DTG file
!     AND the number of individual days in the file. 
! 
! CALLS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/08/08:
!     Written
!
!  ############################################################################

  use type_vars, only: &
    grid_OMI_AI, count_OMI_AI, &
    grid_TROP_AI, count_TROP_AI, &
    grid_MODIS_CH7, count_MODIS_CH7, &
    count_MODIS_CLD_0, count_MODIS_CLD_1, &
    count_MODIS_CLD_2, count_MODIS_CLD_3, &
    count_NSIDC_ICE_0100, count_NSIDC_ICE_251, &
    count_NSIDC_ICE_252, count_NSIDC_ICE_253, &
    count_NSIDC_ICE_254

  implicit none

  integer               :: i_num_days
  integer               :: i_lat_size
  integer               :: i_lon_size

  integer               :: ii
  integer               :: jj
  integer               :: kk

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  do ii = 1, i_num_days
    do jj = 1, i_lat_size
      do kk = 1, i_lon_size
   
        ! OMI AI raw data
        ! ---------------
        if(count_OMI_AI(ii,jj,kk) /= -9) then
          grid_OMI_AI(ii,jj,kk) = &
            grid_OMI_AI(ii,jj,kk) / count_OMI_AI(ii,jj,kk)
        endif

        ! TROP AI data
        ! ---------------
        if(count_TROP_AI(ii,jj,kk) /= -9) then
          grid_TROP_AI(ii,jj,kk) = &
            grid_TROP_AI(ii,jj,kk) / count_TROP_AI(ii,jj,kk)
        endif

        ! MODIS CH7 data
        ! ---------------
        if(count_MODIS_CH7(ii,jj,kk) /= -9) then
          !if(ii == 9) then
          !  write(*,*) 'Before', grid_MODIS_CH7(ii,jj,kk), count_MODIS_CH7(ii,jj,kk)
          !endif
          grid_MODIS_CH7(ii,jj,kk) = &
            grid_MODIS_CH7(ii,jj,kk) / count_MODIS_CH7(ii,jj,kk)
          !if(ii == 9) then
          !  write(*,*) 'After', grid_MODIS_CH7(ii,jj,kk), count_MODIS_CH7(ii,jj,kk)
          !endif
        endif

        ! Call subroutine to handle the MODIS cloud flag averaging
        ! --------------------------------------------------------
        if((count_MODIS_CLD_0(ii,jj,kk) > 0) .or. &
           (count_MODIS_CLD_1(ii,jj,kk) > 0) .or. &
           (count_MODIS_CLD_2(ii,jj,kk) > 0) .or. &
           (count_MODIS_CLD_3(ii,jj,kk) > 0)) then
          call bin_MODIS_CLD_values(ii,jj,kk)
        endif

        ! Call subroutine to handle the NSIDC surface type info
        ! -----------------------------------------------------
        if((count_NSIDC_ICE_0100(ii,jj,kk) > 0) .or. &
           (count_NSIDC_ICE_251(ii,jj,kk)  > 0) .or.  &
           (count_NSIDC_ICE_252(ii,jj,kk)  > 0) .or.  &
           (count_NSIDC_ICE_253(ii,jj,kk)  > 0) .or. &
           (count_NSIDC_ICE_254(ii,jj,kk)  > 0)) then
          call bin_NSIDC_ICE_values(ii,jj,kk)
        endif

      enddo
    enddo
  enddo

  
end subroutine average_final_grid_data
