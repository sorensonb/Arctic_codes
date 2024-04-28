subroutine write_force_r2(out_file_id)
!
! NAME:
!   write_force_r2.f90
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
  use force_eff_vars, only: regress_r2, regress_dims

  implicit none

  integer(hsize_t)                :: out_file_id
  integer(hsize_t)                :: dspace_id
  integer(hsize_t)                :: dset_id

  ! File read variables
  integer                :: rank
  integer                :: error

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  !
  ! Write regression slope r2
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  rank = size(regress_dims)

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, regress_dims, dspace_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'regress_r2', H5T_NATIVE_REAL, &
                    dspace_id, dset_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, regress_r2, &
                   regress_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -------------------------------
  call h5dclose_f(dset_id, error)

  ! Close access to data space rank
  ! -------------------------------
  call h5sclose_f(dspace_id, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close dataset'
    return
  endif

  write(*,*) 'Wrote regress_r2'

end subroutine write_force_r2
