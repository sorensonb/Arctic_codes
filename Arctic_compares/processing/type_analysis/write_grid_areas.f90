subroutine write_grid_areas(out_file_id, grid_area_dims)
!
! NAME:
!   write_grid_areas.f90
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
  use type_vars, only: grid_areas

  implicit none

  integer                         :: out_file_id
  integer                         :: dspace_id
  integer                         :: dset_id

  integer(hsize_t), dimension(2)  :: grid_area_dims

  ! File read variables
  integer                         :: error

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! = = = = = = = = = = = 
  !
  ! Write the data
  !
  ! = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(2, grid_area_dims, dspace_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'grid_areas', H5T_NATIVE_REAL, &
                   dspace_id, dset_id, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'grid_areas'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, grid_areas, grid_area_dims, &
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

  write(*,*) 'Wrote grid_areas'

end subroutine write_grid_areas
