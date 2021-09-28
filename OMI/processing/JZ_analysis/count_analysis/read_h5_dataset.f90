subroutine read_h5_AI(file_id,error)
!subroutine read_h5_dataset(file_id,data_path,AI_dims,&
!    H52DDoubledataset,error)
!
! NAME:
!
! PURPOSE:
! 
! callS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2018/10/24:
!     Written
!
!  ############################################################################

  use hdf5
  use h5_vars, only : AI_dims, AI_data

  implicit none
  
  integer :: file_id
  character(len=*) :: data_path
  !integer(hsize_t),dimension(:),allocatable :: AI_dims
  !!real(kind=8), dimension(:,:), allocatable, target, intent(out)     :: H52DDoubledataset
  integer :: error

  integer :: ii 
  integer :: ds_id
  integer :: dspace
  integer :: ndims
  integer(hsize_t), dimension(1)                 :: dims
  integer(hsize_t), dimension(4)                 :: datadims
  integer(hsize_t), dimension(4)                 :: maxdatadims
 
  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  write(*,*) "Inside read_h5_dataset"
 
  ! Open dataset
  call h5dopen_f(file_id, 'HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/'&
        //'UVAerosolIndex', ds_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset'
    return
  endif
  write(*,*) 'Dataset opened'

  call h5dget_space_f(ds_id, dspace, error)
  if (error /= 0) then
    write(*,*) " FATAL ERROR: Error determining dataspace"
    return
  endif

  ! Determine the number of dimensions in the dataset, allocate the AI_dims
  ! array
  ! -------------------------------------------------------------------------
  call h5sget_simple_extent_ndims_f(dspace, ndims, error)
  if (error < 0) then
    write(*,*) " *** Error determining dataspace dimensionality"
    return
  endif

  allocate(AI_dims(ndims),stat=error)
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

  ! Insert important dimensions into the AI_dims array
  do ii=1,ndims
    AI_dims(ii) = datadims(ii)
  enddo
  
  ! Read the dataset and transfer the result to an allocated working array
  ! ----------------------------------------------------------------------
  allocate(AI_data(AI_dims(1), AI_dims(2)), stat=error)
  if ( error < 0 ) then
     write(*,*) " *** Error allocating H5dataset"
     return
  endif

  call h5dread_f(ds_id, H5T_NATIVE_DOUBLE, AI_data, dims, &
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
  write(*,*) 'Dataset closed'


end subroutine read_h5_AI
