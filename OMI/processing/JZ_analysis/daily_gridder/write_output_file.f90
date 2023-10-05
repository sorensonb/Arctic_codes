subroutine write_output_file(out_file_name, i_num_days, i_lat_size, i_lon_size)
!
! NAME:
!   write_output_file.f90
!
! PURPOSE:
! 
! CALLS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/08/08:
!     Written
!
!  ############################################################################

  use hdf5

  implicit none

  character(len = 255)   :: out_file_name

  integer                :: i_num_days
  integer                :: i_lat_size
  integer                :: i_lon_size

  integer                :: out_file_id

  integer(hsize_t), dimension(3)                   :: grid_data_dims
  integer(hsize_t), dimension(1)                   :: day_dims
  integer(hsize_t), dimension(1)                   :: lat_dims
  integer(hsize_t), dimension(1)                   :: lon_dims

  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer                :: error

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  write(*,*) trim(out_file_name)

  out_file_id = 2023

  ! Set up the dimensions for the output variables
  ! ----------------------------------------------

  !grid_data_dims = (/i_num_days,i_lat_size,i_lon_size/)
  grid_data_dims = (/i_lon_size,i_lat_size,i_num_days/)
  day_dims       = (/i_num_days/)
  lat_dims       = (/i_lat_size/)
  lon_dims       = (/i_lon_size/)

  ! Open the output file
  ! --------------------

  call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open output file'
    return
  endif

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Call subroutines to write variables to output file
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  call write_day_dates(out_file_id, day_dims)
  call write_day_lat(out_file_id, lat_dims)
  call write_day_lon(out_file_id, lon_dims)
  call write_day_grid_AI(out_file_id, grid_data_dims)
  call write_day_grid_GPQF(out_file_id, grid_data_dims)
 
  !!#!call write_day_values(out_file_id, day_dims)
  !!#!call write_lat_values(out_file_id, lat_dims)
  !!#!call write_lon_values(out_file_id, lon_dims)
  !!#!call write_grid_areas(out_file_id, grid_area_dims)
  !!#!call write_grid_OMI_AI_raw(out_file_id, grid_data_dims)
  !!#!call write_grid_TROP_AI(out_file_id, grid_data_dims)
  !!#!call write_grid_MODIS_CH7(out_file_id, grid_data_dims)
  !!#!call write_grid_MODIS_CLD(out_file_id, grid_data_dims)
  !!#!call write_grid_NSIDC_ICE(out_file_id, grid_data_dims)

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! End subroutine calls
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Close the output file
  ! ---------------------
  call h5fclose_f(out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif
  
end subroutine write_output_file
