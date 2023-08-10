subroutine grid_comp_data(single_file_date, lat_thresh, & 
                          i_day_idx, l_debug)
!
! NAME:
!   grid_comp_data.f90
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

  use hdf5
  use type_vars, only: &
    grid_OMI_AI, count_OMI_AI, &
    grid_TROP_AI, count_TROP_AI, &
    grid_MODIS_CH7, count_MODIS_CH7, &
    grid_NSIDC_ICE, &
    count_MODIS_CLD_0, count_MODIS_CLD_1, &
    count_MODIS_CLD_2, count_MODIS_CLD_3, &
    count_NSIDC_ICE_0100, count_NSIDC_ICE_251, &
    count_NSIDC_ICE_252,  count_NSIDC_ICE_253, &
    count_NSIDC_ICE_254, &
    lat_range, lon_range
  use comp_grid_vars, only : &
    clear_arrays, &
    MODIS_CH2_data, MODIS_CH2_dims, &
    MODIS_CH7_data, MODIS_CH7_dims, &
    MODIS_CLD_data, MODIS_CLD_dims, &
    MODIS_COD_data, MODIS_COD_dims, &
    !MODIS_LAT_data, MODIS_LAT_dims, &
    !MODIS_LON_data, MODIS_LON_dims, &
    NSIDC_data,     NSIDC_dims, &
    !NSIDC_LAT_data, NSIDC_LAT_dims, &
    !NSIDC_LON_data, NSIDC_LON_dims, &
    !CERES_CLD_data, CERES_CLD_dims, &
    CERES_LWF_data, CERES_LWF_dims, &
    CERES_SWF_data, CERES_SWF_dims, &
    !CERES_LAT_data, CERES_LAT_dims, &
    !CERES_LON_data, CERES_LON_dims, &
    OMI_AI_data,    OMI_AI_dims, &
    OMI_AI_raw_data,    OMI_AI_raw_dims, &
    OMI_LAT_data,   OMI_LAT_dims, &
    OMI_LON_data,   OMI_LON_dims, &
    !OMI_LATCRNR_data,   OMI_LATCRNR_dims, &
    !OMI_LONCRNR_data,   OMI_LONCRNR_dims, &
    OMI_SZA_data,   OMI_SZA_dims, &
    OMI_VZA_data,   OMI_VZA_dims, &
    OMI_AZM_data,   OMI_AZM_dims, &
    TROP_AI_data,   TROP_AI_dims, &
    TROP_SSA0_data, TROP_SSA0_dims, &
    TROP_SSA1_data, TROP_SSA1_dims, &
    TROP_SSA2_data, TROP_SSA2_dims
    !!#!MODIS_out_CH2_data, &
    !!#!MODIS_out_CH7_data, &
    !!#!MODIS_out_CLD_data, &
    !!#!!MODIS_out_LAT_data, &
    !!#!!MODIS_out_LON_data, &
    !!#!NSIDC_out_data,     &
    !!#!!NSIDC_out_LAT_data, &
    !!#!!NSIDC_out_LON_data, &
    !!#!CERES_out_LWF_data, &
    !!#!CERES_out_SWF_data   
    !!#!!CERES_out_LAT_data, &
    !!#!!CERES_out_LON_data

  implicit none

  character(len = 12)    :: single_file_date
  character(len = 255)   :: trop_name      ! filename
  character(len = 255)   :: total_file_name ! name of file containing the 
                                            ! list of shawn file names to 
                                            ! be analyzed

  real                   :: lat_thresh
  real                   :: local_trop_ai

  integer                :: ii
  integer                :: jj

  integer                :: i_day_idx
  integer                :: i_lat_idx
  integer                :: i_lon_idx

  logical                :: l_trop_found
  logical                :: l_debug
 
  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer                :: file_id
  integer                :: error

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


  ! Build the combined subset file name
  ! -----------------------------------
  total_file_name = '/home/bsorenson/OMI/arctic_comp/comp_data/'&
    //'colocated_subset_'//single_file_date//'.hdf5'

  ! Using the OMI name, determine if a matching TROPOMI colocation 
  ! file exists
  ! --------------------------------------------------------------
  l_trop_found = .false.
  trop_name = '/home/bsorenson/OMI/tropomi_colocate/coloc_data/'//&
    'colocated_tropomi_'//single_file_date//'.hdf5'
  
  inquire(FILE=trop_name, EXIST=l_trop_found)
  if(l_debug) then
    write(*,*) 'TROP file exist? :',l_trop_found,single_file_date
  endif


  ! Try to open the file
  ! --------------------
  ! Open the HDF5 file
  ! ------------------
  call h5fopen_f(total_file_name, H5F_ACC_RDWR_F, file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open file '//total_file_name
    return
  endif
       

  ! Open the HDF5 coloc file named total_file_name
  ! ----------------------------------------------
  call read_coloc_OMI_LAT(file_id)
  call read_coloc_OMI_LON(file_id)
  call read_coloc_OMI_AI_raw(file_id)
  call read_coloc_MODIS_CH7(file_id)
  call read_coloc_MODIS_CLD(file_id)
  call read_coloc_NSIDC_ICE(file_id)
  if(l_trop_found) then
    call read_coloc_TROP_AI(file_id)
  endif

  ! Loop over the coloc file, inserting the data into the grid
  ! ----------------------------------------------------------
  omi_loop1: do ii=1,OMI_AI_raw_dims(2)
    omi_loop2: do jj=1,OMI_AI_raw_dims(1) 

    if(l_trop_found) then
      local_trop_ai = TROP_AI_data(jj,ii)
    else
      local_trop_ai = -999.
    endif

    ! Check if the current pixel is missing
    ! -------------------------------------
    !if( (OMI_LAT_data(jj,ii) > lat_thresh) .and. &
    !    .not.isnan(OMI_AI_raw_data(jj,ii))) then
    if( (l_trop_found .and. ((.not.isnan(local_trop_ai) .and. &
         (local_trop_ai/= -999.)) .or. &
        .not.isnan(OMI_AI_raw_data(jj,ii)))) .or. &
        (.not.l_trop_found .and. .not.isnan(OMI_AI_raw_data(jj,ii)))) then

      ! Figure out which lat/lon bins match the current OMI pixel
      ! ---------------------------------------------------------
      i_lat_idx = minloc(abs(lat_range - OMI_LAT_data(jj,ii)), dim = 1)
      i_lon_idx = minloc(abs(lon_range - OMI_LON_data(jj,ii)), dim = 1)

      ! Grid the OMI AI data
      ! --------------------
      if(.not. isnan(OMI_AI_raw_data(jj,ii))) then
        if(count_OMI_AI(i_day_idx, i_lat_idx, i_lon_idx) == -9) then
          grid_OMI_AI(i_day_idx, i_lat_idx, i_lon_idx) = OMI_AI_raw_data(jj,ii)
          count_OMI_AI(i_day_idx, i_lat_idx, i_lon_idx) = 1
        else
          grid_OMI_AI(i_day_idx, i_lat_idx, i_lon_idx) = &
            grid_OMI_AI(i_day_idx, i_lat_idx, i_lon_idx) + OMI_AI_raw_data(jj,ii)
          count_OMI_AI(i_day_idx, i_lat_idx, i_lon_idx) = & 
            count_OMI_AI(i_day_idx, i_lat_idx, i_lon_idx) + 1
        endif
      endif

      if(l_trop_found) then
        if((count_TROP_AI(i_day_idx, i_lat_idx, i_lon_idx) == -9) .and. &
           (abs(local_trop_ai) < 20)) then
          grid_TROP_AI(i_day_idx, i_lat_idx, i_lon_idx) = local_trop_ai
          count_TROP_AI(i_day_idx, i_lat_idx, i_lon_idx) = 1
        else
          grid_TROP_AI(i_day_idx, i_lat_idx, i_lon_idx) = &
            grid_TROP_AI(i_day_idx, i_lat_idx, i_lon_idx) + local_trop_ai
          count_TROP_AI(i_day_idx, i_lat_idx, i_lon_idx) = & 
            count_TROP_AI(i_day_idx, i_lat_idx, i_lon_idx) + 1
        endif
      endif

      ! Add the MODIS data
      ! ------------------
      if(MODIS_CH7_data(jj,ii) > 0) then
        if(count_MODIS_CH7(i_day_idx, i_lat_idx, i_lon_idx) == -9) then
          grid_MODIS_CH7(i_day_idx, i_lat_idx, i_lon_idx) = MODIS_CH7_data(jj,ii)
          count_MODIS_CH7(i_day_idx, i_lat_idx, i_lon_idx) = 1
        else
          grid_MODIS_CH7(i_day_idx, i_lat_idx, i_lon_idx) = &
            grid_MODIS_CH7(i_day_idx, i_lat_idx, i_lon_idx) + &
            MODIS_CH7_data(jj,ii)
          count_MODIS_CH7(i_day_idx, i_lat_idx, i_lon_idx) = & 
            count_MODIS_CH7(i_day_idx, i_lat_idx, i_lon_idx) + 1
        endif
      endif

      if(MODIS_CLD_data(jj,ii) == 0.) then
        count_MODIS_CLD_0(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_MODIS_CLD_0(i_day_idx, i_lat_idx, i_lon_idx) + 1
      else if(MODIS_CLD_data(jj,ii) == 1.) then
        count_MODIS_CLD_1(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_MODIS_CLD_1(i_day_idx, i_lat_idx, i_lon_idx) + 1
      else if(MODIS_CLD_data(jj,ii) == 2.) then
        count_MODIS_CLD_2(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_MODIS_CLD_2(i_day_idx, i_lat_idx, i_lon_idx) + 1
      else if(MODIS_CLD_data(jj,ii) == 3.) then
        count_MODIS_CLD_3(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_MODIS_CLD_3(i_day_idx, i_lat_idx, i_lon_idx) + 1
      endif

      ! Add the NSIDC data
      ! ------------------
      if((NSIDC_data(jj,ii) >= 0.) .and. (NSIDC_data(jj,ii) <= 100.)) then
        if(count_NSIDC_ICE_0100(i_day_idx, i_lat_idx, i_lon_idx) == 0) then
          grid_NSIDC_ICE(i_day_idx, i_lat_idx, i_lon_idx) = &
            NSIDC_data(jj,ii)
        else
          grid_NSIDC_ICE(i_day_idx, i_lat_idx, i_lon_idx) = &
            grid_NSIDC_ICE(i_day_idx, i_lat_idx, i_lon_idx) + &
            NSIDC_data(jj,ii)
        endif
        count_NSIDC_ICE_0100(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_NSIDC_ICE_0100(i_day_idx, i_lat_idx, i_lon_idx) + 1
      else if(NSIDC_data(jj,ii) == 251.) then
        count_NSIDC_ICE_251(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_NSIDC_ICE_251(i_day_idx, i_lat_idx, i_lon_idx) + 1
      else if(NSIDC_data(jj,ii) == 252.) then
        count_NSIDC_ICE_252(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_NSIDC_ICE_252(i_day_idx, i_lat_idx, i_lon_idx) + 1
      else if(NSIDC_data(jj,ii) == 253.) then
        count_NSIDC_ICE_253(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_NSIDC_ICE_253(i_day_idx, i_lat_idx, i_lon_idx) + 1
      else if(NSIDC_data(jj,ii) == 254.) then
        count_NSIDC_ICE_254(i_day_idx, i_lat_idx, i_lon_idx) = &
          count_NSIDC_ICE_254(i_day_idx, i_lat_idx, i_lon_idx) + 1
      endif
    endif

    enddo omi_loop2
  enddo omi_loop1

  ! Deallocate the arrays holding this HDF5 coloc file's data
  ! ---------------------------------------------------------
  call clear_arrays

  ! Close file
  ! ----------
  call h5fclose_f(file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close file'
    return
  endif
  
end subroutine grid_comp_data
