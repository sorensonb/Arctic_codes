subroutine write_coloc_OMI_VZA(out_file_id)
!
! NAME:
!   write_grid_TROP_VZA.f90
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
  use colocate_vars, only: OMI_VZA_data,    OMI_AI_dims

  implicit none

  integer(hsize_t)                :: out_file_id
  integer(hsize_t)                :: dspace_id
  integer(hsize_t)                :: dset_id

  ! File read variables
  integer                :: rank
  integer                :: error

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! = = = = = = = = = = = 
  !
  ! Write the data
  !
  ! = = = = = = = = = = = 

  rank = 2

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, OMI_AI_dims, dspace_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_vza', H5T_NATIVE_DOUBLE, dspace_id,  &
                   dset_id, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_vza'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, OMI_VZA_data, OMI_AI_dims, &
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
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif

  write(*,*) 'Wrote OMI OVZ'
  

end subroutine write_coloc_OMI_VZA
