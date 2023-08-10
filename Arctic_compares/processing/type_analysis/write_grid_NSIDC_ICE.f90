subroutine write_grid_NSIDC_ICE(out_file_id, grid_data_dims)
!
! NAME:
!   write_grid_NSIDC_ICE.f90
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
  use type_vars, only: &
    grid_NSIDC_ICE

  implicit none

  integer                :: out_file_id
  integer                :: dspace_id
  integer                :: dset_id

  integer(hsize_t), dimension(3)                   :: grid_data_dims

  ! File read variables
  integer                :: error

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! = = = = = = = = = = = 
  !
  ! Write the data
  !
  ! = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(3, grid_data_dims, dspace_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'nsidc_ice', H5T_NATIVE_REAL, &
                   dspace_id, dset_id, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_ice'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, grid_NSIDC_ICE, grid_data_dims, &
                      error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close dataset'
    return
  endif

  write(*,*) 'Wrote nsidc_ice'

  !!#!! = = = = = = = = = = = 
  !!#!!
  !!#!! Write the counts
  !!#!!
  !!#!! = = = = = = = = = = = 

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(3, grid_data_dims, dspace_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'modis_ch7_count', H5T_NATIVE_INTEGER, &
  !!#!                 dspace_id, dset_id, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_ch7_count'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, count_NSIDC_ICE, grid_data_dims, &
  !!#!                    error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote modis_ch7_count'

end subroutine write_grid_NSIDC_ICE
