program ceres_grid
!
! NAME:
!   ceres_grid.f90
!
! PURPOSE:
!   Calculate monthly averages of screened AI (following the JZ criteria) and
!   write the averages to an output file.
! 
! CALLS:
!   Modules:
!     - hdf5
!     - h5_vars (custom module, shared between count and climo analyses)
!   Subroutines:
!     - check_bad_rows
!     - clear_arrays (from h5_vars)
!     - grid_raw_data_climo
!     - print_climo
!     - read_h5_AI
!     - read_h5_LAT
!     - read_h5_LON
!     - read_h5_XTRACK
!     - read_h5_AZM
!     - read_h5_GPQF
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/09:
!     Written
!
!  ############################################################################

  use hdf5
  use ceres_vars, only : clear_arrays, i_bad_list, &
    find_distance_between_points, allocate_out_arrays, &
    pixel_in_box, &
    CERES_LAT_data, CERES_LAT_dims, &
    CERES_LON_data, CERES_LON_dims, &
    CERES_SWF_data, CERES_SWF_dims, &
    CERES_LWF_data, CERES_LWF_dims, &
    CERES_TIM_data, CERES_TIM_dims, &
    CERES_ALB_data, CERES_ALB_dims, &
    CERES_SZA_data, CERES_SZA_dims, &
    CERES_CLS_data, CERES_CLS_dims, &
    CERES_CLW_data, CERES_CLW_dims, &
    CERES_out_LWF_data, &
    CERES_out_SWF_data, &
    CERES_out_ALB_data, &
    CERES_out_swcnt_data, &
    CERES_out_lwcnt_data, &
    CERES_out_abcnt_data

  implicit none

  integer                :: ii            ! loop counter
  integer                :: jj            ! loop counter
  integer                :: nii            ! loop counter
  integer                :: njj            ! loop counter
  integer                :: index1
  integer                :: index2
  !!#!integer                :: int_month     ! integer variable for month
  !!#!integer                :: int_day       ! integer variable for day 
  integer                :: arg_count     ! Number of arguments passed to exec
  !!#!integer                :: work_month    ! currently analyzed month
  !!#!integer                :: work_day 
  integer                :: error         ! error flag
  integer                :: istatus
  !integer                :: num_nan
  integer                :: ceres_file_id       ! id for current HDF5 file
  integer                :: out_file_id       ! id for current HDF5 file

  integer                :: dspace_id_CLT  ! CERES LAT
  integer                :: dspace_id_CLN  ! CERES LON
  integer                :: dspace_id_CSW  ! CERES SWF
  integer                :: dspace_id_CLW  ! CERES LWF
  integer                :: dspace_id_CAL  ! CERES LWF
  integer                :: dspace_id_CSW_cnt  ! CERES SWF
  integer                :: dspace_id_CLW_cnt  ! CERES LWF
  integer                :: dspace_id_CAL_cnt  ! CERES LWF

  integer                :: dset_id_CLT  ! CERES LAT
  integer                :: dset_id_CLN  ! CERES LON
  integer                :: dset_id_CSW  ! CERES SWF
  integer                :: dset_id_CLW  ! CERES LWF
  integer                :: dset_id_CAL  ! CERES LWF
  integer                :: dset_id_CSW_cnt  ! CERES LWF
  integer                :: dset_id_CLW_cnt  ! CERES LWF
  integer                :: dset_id_CAL_cnt  ! CERES LWF

  integer                :: lats(30) = (/(ii, ii = 60, 89, 1)/)
  integer                :: lons(360) = (/(ii, ii = 0, 359, 1)/)

  integer                :: rank
  integer(hsize_t), dimension(2)    :: test_dims
  integer(hsize_t), dimension(1)    :: lat_dims
  integer(hsize_t), dimension(1)    :: lon_dims

  real                   :: arr_size        ! array size
  real                   :: lat_gridder
  !real                   :: distance
  !real                   :: closest_dist
  !real                   :: min_dist
  !real                   :: match_lat
  !real                   :: match_lon
  !real                   :: match_data1
  !real                   :: match_data2

  !real                   :: run_modis_total_ch2
  !real                   :: run_modis_total_ch7
  !integer                :: count_modis_total 

  character(len = 255)      :: ceres_name     ! filename
  character(len = 255)      :: out_file_name     ! filename

  !!#!! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 1) then
    write(*,*) 'SYNTAX: ./omi_exec ceres_aqua_name '! &
!        //'ceres_name omi_name'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1,ceres_name)
  !call get_command_argument(2,modis_name2)
  !call get_command_argument(3,nsidc_name)
  !call get_command_argument(4,ceres_name)
  !call get_command_argument(5,omi_name)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"

  ! Open the CERES file
  ! -------------------
  call h5fopen_f(trim(ceres_name), H5F_ACC_RDWR_F, ceres_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open CERES file'
    return
  endif

  call read_CERES_LAT(ceres_file_id)
  call read_CERES_LON(ceres_file_id)
  call read_CERES_SWF(ceres_file_id)
  call read_CERES_LWF(ceres_file_id)
  call read_CERES_SZA(ceres_file_id)
  call read_CERES_ALB(ceres_file_id)

  !do ii = 1, size(lons)
  !  write(*,*) lons(ii) + 0.5
  !enddo

  !do ii = 1, size(lats)
  !  write(*,*) lats(ii) + 0.5
  !enddo

  ! Average the data into the grid
  ! -------------------------------
  lat_gridder = 60.
  arr_size = 90. - lat_gridder
  ii = 333456

  write(*,*) CERES_LAT_data(ii), CERES_LON_data(ii)
  !write(*,*) index1, index2

  write(*,*) CERES_LAT_dims

  ! Allocate the output arrays
  ! --------------------------
  call allocate_out_arrays(size(lats), size(lons))


  ! Grid the data
  ! -------------
  ceres_loop: do ii = 1, CERES_LAT_dims(1)
  
    index1 = floor(CERES_LAT_data(ii) - lat_gridder) + 1
    index2 = floor(CERES_LON_data(ii) + 180) + 1

    if(index1 < 1) index1 = 1
    if(index1 > arr_size) index1 = arr_size
    if(index2 < 1) index2 = 1
    if(index2 > 360) index2 = 360
 
    ! Check the requirements for the pixel
    ! ------------------------------------
    if((CERES_SWF_data(ii) > 0.) .and. &
       (CERES_SWF_data(ii) < 5000.)) then

      CERES_out_SWF_data(index1, index2) = &
        CERES_out_SWF_data(index1, index2) + CERES_SWF_data(ii)
      CERES_out_swcnt_data(index1, index2) = &
        CERES_out_swcnt_data(index1, index2) + 1
    endif

    if((CERES_LWF_data(ii) > 0.) .and. &
       (CERES_LWF_data(ii) < 5000.)) then
   
      CERES_out_LWF_data(index1, index2) = &
        CERES_out_LWF_data(index1, index2) + CERES_LWF_data(ii) 
      CERES_out_lwcnt_data(index1, index2) = &
        CERES_out_lwcnt_data(index1, index2) + 1
    endif

    if((CERES_ALB_data(ii) > 0.) .and. &
       (CERES_ALB_data(ii) < 5000.)) then
   
      CERES_out_ALB_data(index1, index2) = &
        CERES_out_ALB_data(index1, index2) + CERES_ALB_data(ii) 
      CERES_out_abcnt_data(index1, index2) = &
        CERES_out_abcnt_data(index1, index2) + 1
    endif

  enddo ceres_loop

  CERES_out_SWF_data(:,:) = CERES_out_SWF_data(:,:) / CERES_out_swcnt_data(:,:)
  CERES_out_LWF_data(:,:) = CERES_out_LWF_data(:,:) / CERES_out_lwcnt_data(:,:)
  CERES_out_ALB_data(:,:) = CERES_out_ALB_data(:,:) / CERES_out_abcnt_data(:,:)

  !!#!! Look at the gridded output
  !!#!! --------------------------
  !!#!lat_loop: do ii = 1, size(lons)
  !!#!  lon_loop: do jj = 1, size(lats)
  !!#!    write(*,*) lats(jj), lons(ii), &
  !!#!      CERES_out_SWF_data(jj,ii), CERES_out_LWF_data(jj,ii), &
  !!#!      CERES_out_swcnt_data(jj,ii), CERES_out_lwcnt_data(jj,ii)
  !!#!      !CERES_out_SWF_data(jj,ii) / CERES_out_counts_data(jj,ii), &
  !!#!      !CERES_out_counts_data(jj,ii)
  !!#! 
  !!#!  enddo lon_loop
  !!#!enddo lat_loop


  ! Open the output file
  ! --------------------
  write(*,*) ceres_name(len(trim(ceres_name)) - 19:len(trim(ceres_name)) - 5)
  out_file_name = 'ceres_griddedL2_'//ceres_name(len(trim(ceres_name)) - &
    19:len(trim(ceres_name)) - 5)//'.hdf5'
  write(*,*) out_file_name

  call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open output file'
    return
  endif

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write OMI Latitude
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  rank = 2
  test_dims = (/30, 360/)
  lat_dims  = (/30/)
  lon_dims  = (/360/)

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, lat_dims, dspace_id_CLT, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'lat', H5T_NATIVE_INTEGER, dspace_id_CLT,  &
  !!#!                 dset_id_CLT, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'lat'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_CLT, H5T_NATIVE_INTEGER, lats, &
  !!#!                  lat_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_CLT, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_CLT, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI LAT'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write OMI Longitude
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, lon_dims, dspace_id_CLN, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'lon', H5T_NATIVE_INTEGER, dspace_id_CLN,  &
  !!#!                 dset_id_CLN, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'lon'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_CLN, H5T_NATIVE_INTEGER, lons, &
  !!#!                  lon_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_CLN, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_CLN, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close output file'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI LON'
  
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write CERES SWF data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_CSW, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'ceres_swf', H5T_NATIVE_DOUBLE, &
                   dspace_id_CSW,  dset_id_CSW, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_swf'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_CSW, H5T_NATIVE_DOUBLE, CERES_out_SWF_data, &
                  test_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_CSW, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_CSW, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote CERES SWF'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write CERES LWF data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_CLW, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'ceres_lwf', H5T_NATIVE_DOUBLE, &
                   dspace_id_CLW,  dset_id_CLW, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_lwf'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_CLW, H5T_NATIVE_DOUBLE, CERES_out_LWF_data, &
                  test_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_CLW, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_CLW, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote CERES LWF'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write CERES ALB data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_CAL, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'ceres_alb', H5T_NATIVE_DOUBLE, &
                   dspace_id_CAL,  dset_id_CAL, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_alb'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_CAL, H5T_NATIVE_DOUBLE, CERES_out_ALB_data, &
                  test_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_CAL, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_CAL, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote CERES ALB'





  ! Close output file
  ! -----------------
  call h5fclose_f(out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif

  write(*,*) 'Saved output file'//trim(out_file_name)
 
  ! Deallocate all the arrays for the next pass
  ! -------------------------------------------
  call clear_arrays

  ! Close file
  ! ----------
  call h5fclose_f(ceres_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close CERES file'
    return
  endif

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program ceres_grid
