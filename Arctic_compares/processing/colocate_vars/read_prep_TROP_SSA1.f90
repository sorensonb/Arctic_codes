subroutine read_prep_TROP_SSA1(file_id)
!
! NAME:
!   read_prep_TROP_SSA1
!
! PURPOSE:
!   Read Latitude data from an HDF5 file pointed to by file_id 
!   and store the data in SSA1_SSA1 (and the dimensions in SSA1_dims) from
!   h5_vars.
!  
! CALLS:
!   Modules:
!     - hdf5
!     - h5_vars (custom module, shared between count and climo analyses)
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/09:
!     Written

!  ############################################################################

  use hdf5
  use colocate_vars, only : TROP_prep_SSA1_dims, TROP_prep_SSA1_data

  implicit none
  
  integer                         :: file_id     ! File id for the HDF5 file
  integer                         :: error       ! error flag

  integer                         :: ii          ! loop counter
  integer                         :: ds_id       ! dataset ID
  integer                         :: dspace      ! dataspace 
  integer                         :: ndims       ! number of dimensions in file
  integer(hsize_t), dimension(1)  :: dims        ! dimensions
  integer(hsize_t), dimension(4)  :: datadims    ! all data dimensions
  integer(hsize_t), dimension(4)  :: maxdatadims ! max data dims
 
  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Open data dataset
  ! ----------------
  !call h5dopen_f(file_id, 'HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation'&
  !      //' Fields/Latitude', ds_id, error)
  call h5dopen_f(file_id, 'ssa1', ds_id, error)
  !call h5dopen_f(file_id, 'latitude', ds_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset'
    return
  endif

  ! Determine the dataset dataspace
  ! -------------------------------
  call h5dget_space_f(ds_id, dspace, error)
  if (error /= 0) then
    write(*,*) " FATAL ERROR: Error determining dataspace"
    return
  endif

  ! Determine the number of dimensions in the dataset, allocate the data_dims
  ! array
  ! -------------------------------------------------------------------------
  call h5sget_simple_extent_ndims_f(dspace, ndims, error)
  if (error < 0) then
    write(*,*) " *** Error determining dataspace dimensionality"
    return
  endif

  allocate(TROP_prep_SSA1_dims(ndims),stat=error)
  if ( error /= 0 ) then
     write(*,*) " *** Error allocating dims"
     return
  endif

  ! Determine the dimensions in the dataset
  ! ---------------------------------------
  call h5sget_simple_extent_dims_f(dspace, datadims, maxdatadims,&
       error)
  if (error < 0) then
     write(*,*) " *** Error determining dataspace size"
     return
  endif

  ! Insert important dimensions into the data_dims array
  !TROP_prep_SSA1_dims(1) = datadims(1)
  do ii=1,ndims
    TROP_prep_SSA1_dims(ii) = datadims(ii)
  enddo

  ! Read the dataset and transfer the result to an allocated working array
  ! ----------------------------------------------------------------------
  allocate(TROP_prep_SSA1_data(TROP_prep_SSA1_dims(1)), stat=error)
  if ( error < 0 ) then
     write(*,*) " *** Error allocating H5dataset"
     return
  endif

  call h5dread_f(ds_id, H5T_NATIVE_DOUBLE, TROP_prep_SSA1_data, dims, &
                 error)
  if (error.lt.0) then
      write(*,*) " *** Error reading data"
      return
  endif

  ! Close dataset
  call h5dclose_f(ds_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close dataset'
    return
  endif


end subroutine read_prep_TROP_SSA1