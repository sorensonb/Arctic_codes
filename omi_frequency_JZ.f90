program omi_frequency_JZ
!
! NAME:
!   omi_frequency_JZ.f90
!
! PURPOSE:
! 
! callS:
!   mie_calc.f90
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2018/10/24:
!     Written
!
!  ############################################################################

  use hdf5
  use h5_vars, only : AI_dims, LAT_dims, LON_dims, XTRACK_dims, &
                      AI_data, LAT_data, LON_data, XTRACK_data, &
                      clear_arrays

  implicit none

  integer                :: ii            ! loop counter
  integer                :: synop_idx
  integer                :: i_size        ! array size
  integer                :: ai_count      ! good AI counter
  integer,dimension(3)   :: ai_count2     ! good AI counter
  integer                :: int_hr        ! integer variable for hour

  real                   :: ai_thresh     ! threshold AI value
  real                   :: lat_thresh
  real                   :: lat_gridder

  logical                :: l_in_time

  real,dimension(:), allocatable :: lat_range
  real,dimension(:), allocatable :: lon_range   
  real,dimension(:,:), allocatable :: grids
  integer,dimension(:,:), allocatable :: i_counts

  integer,dimension(4)   :: synop_times 

  ! File read variables
  integer,parameter      :: io8    = 1042   ! File object for file name file
  integer,parameter      :: io6    = 1827 ! Data output file
  integer,parameter      :: errout = 1009 
  integer                :: istatus

  character(len = 255)   :: data_path
  character(len = 255)   :: out_file_name
  character(len = 255)   :: file_name_file
  character(len = 255)   :: total_file_name

  !! Variables from each line in Shawn's files
  !real                   :: lat
  !real                   :: lon
  !real                   :: raw_ai
  !real                   :: filter
  !real                   :: clean_ai
  !real                   :: v5
  !real                   :: v6
  !real                   :: v7
  !real                   :: v8
  !real                   :: v9
  !real                   :: v10
  !real                   :: v11
  !real                   :: v12
  !real                   :: v13
  !real                   :: v14

  integer                 :: arg_count

  integer :: error
  integer :: file_id
  integer :: gr_id

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./omi_exec out_file_name file_name_file'
    return
  endif

  ! Initialize the HDF5 interface
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"


  call get_command_argument(1,out_file_name)
  call get_command_argument(2,file_name_file)

  write(*,*) "out_file_name = ",trim(out_file_name)
  write(*,*) "file_name_file = ",trim(file_name_file)

  ai_count2(:) = 0

  synop_times = [0,6,12,18] 
  synop_idx = 1

  !data_path = "/home/bsorenson/OMI/shawn_analysis/test_dir/"
  data_path = "/Research/OMI/H5_files/"
  !out_file_name = "omi_counts_200501_200909.txt"
  !file_name_file = "omi_dates_200501_200909.txt"

  ! Set up lat/lon grids
  ! --------------------
  lat_thresh = 65.
  lat_gridder = lat_thresh * 4.
  i_size = (90. - lat_thresh) * 4.
  allocate(lat_range(i_size))
  do ii=1,i_size
    lat_range(ii) = lat_thresh + (ii-1)*0.25
  enddo

  allocate(lon_range(1440))
  do ii=1,1440
    lon_range(ii) = -180.0 + (ii-1)*0.25
  enddo


  ! Initialize grid arrays and set to -9 initially
  ! ----------------------------------------------

  allocate(grids(1440,i_size))
  allocate(i_counts(1440,i_size))

  grids(:,:) = 0.
  i_counts(:,:) = 0

  ! open debug file
  ! ---------------
  open(errout, file = "omi_jz_error.txt", iostat = istatus)
  if(istatus /= 0) then
    write(*,*) "error opening error file."
  endif

  write(errout,*) "TEST WRITE TO FILE"

  ! open output file
  ! ---------------
  write(errout,*) "Opening "//trim(out_file_name)
  open(io6, file = trim(out_file_name), iostat = istatus)
  if(istatus /= 0) then
    write(errout,*) "ERROR: error opening data count output file."
  endif
  write(io6,'(a10,5(a6))') 'Date','Cnt65','Cnt70','Cnt75','Cnt80','Cnt85'
 
  ! Set up count variables to count the number of grid boxes with
  ! high AI values
  ! -------------------------------------------------------------
  ai_thresh = 0.6
  ai_count  = 0
  
  ! Read the file names from the file name file
  open(io8, file = trim(file_name_file), iostat = istatus)
  if(istatus > 0) then
    write(errout,*) "ERROR: Problem opening file name containing files: "&
        //'titled '//trim(file_name_file)
    return 
  else
    !call process_files()
    ! Loop over the file
    file_loop: do
      ! Read the current total_file_name from the file
      read(io8, '(A)', iostat=istatus) total_file_name
      if(istatus < 0) then 
        write(*,*) "End of "//trim(file_name_file)//" found"
        exit
      else if(istatus > 0) then
        write(errout,*) "ERROR: problem reading total_file_name"
        cycle file_loop
      else

        ! Extract time information from total_file_name
        ! ---------------------------------
        read(total_file_name(54:55), *) int_hr

        ! See if the hour exceeds the current 6 hr assimilation window.
        ! If so, calculate averages and counts and reset variables.
        ! ------------------------------------------------------------
        call synop_time_check(synop_idx, int_hr, l_in_time)
        if(.not. l_in_time) then 
          call count_ai(io6,grids,i_counts,i_size,ai_thresh,synop_idx,&
                        total_file_name,lat_range)

        endif  

        write(*,*) trim(total_file_name)

        ! # # # # # # # # # # # # # # # #
        ! INSERT H5 READER CODE HERE
        ! # # # # # # # # # # # # # # # #

        ! Open the HDF5 file
        call h5fopen_f(total_file_name, H5F_ACC_RDWR_F, file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: could not open file'
          return
        endif
        !write(*,*) 'File opened'

        !!! Open group
        !!call h5gopen_f(file_id, 'HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/',&
        !!               gr_id, error)
        !!if(error /= 0) then
        !!  write(*,*) 'FATAL ERROR: could not open group'
        !!  return
        !!endif
        !!write(*,*) 'Group opened'

        ! Read in the necessary data
        call read_h5_AI(file_id,error)
        call read_h5_LAT(file_id,error)
        call read_h5_LON(file_id,error)
        call read_h5_XTRACK(file_id,error)   
        !call read_h5_AZM(file_id,error)         !NEED
        !call read_h5_GPQF(file_id,error)        !NEED

        ! Close file
        call h5fclose_f(file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: could not close file'
          return
        endif
        !write(*,*) 'File closed'


        !!#!! Open dataset
        !!#!call h5dopen_f(file_id, 'HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex', ds_id, error)
        !!#!if(error /= 0) then
        !!#!  write(*,*) 'FATAL ERROR: could not open dataset'
        !!#!  return
        !!#!endif
        !!#!write(*,*) 'Dataset opened'

        !!#!call h5dget_space_f(ds_id, dspace, error)
        !!#!if (error /= 0) then
        !!#!  write(*,*) " FATAL ERROR: Error determining dataspace"
        !!#!  return
        !!#!endif

        !!#!! Determine the number of dimensions in the dataset, allocate the main_dims
        !!#!! array
        !!#!! -------------------------------------------------------------------------
        !!#!call h5sget_simple_extent_ndims_f(dspace, ndims, error)
        !!#!if (error < 0) then
        !!#!  write(*,*) " *** Error determining dataspace dimensionality"
        !!#!  return
        !!#!endif

        !!#!allocate(main_dims(ndims),stat=error)
        !!#!if ( error /= 0 ) then
        !!#!   write(*,*) " *** Error allocating dims"
        !!#!   return
        !!#!endif

        !!#!! Determine the dimensions in the dataset
        !!#!! ---------------------------------------
        !!#!call h5sget_simple_extent_dims_f(dspace, datadims, maxdatadims,&
        !!#!     error)
        !!#!if (error < 0) then
        !!#!   write(*,*) " *** Error determining dataspace size"
        !!#!   return
        !!#!endif

        !!#!! Insert important dimensions into the main_dims array
        !!#!do ii=1,ndims
        !!#!  main_dims(ii) = datadims(ii)
        !!#!enddo
        
        !!#!! Read the dataset and transfer the result to an allocated working array
        !!#!! ----------------------------------------------------------------------
        !!#!allocate(H52DDoubledataset(main_dims(1), main_dims(2)), stat=error)
        !!#!if ( error < 0 ) then
        !!#!   write(*,*) " *** Error allocating H5dataset"
        !!#!   return
        !!#!endif

        !!#!call h5dread_f(ds_id, H5T_NATIVE_DOUBLE, H52DDoubledataset, dims, &
        !!#!               error)
        !!#!if (error.lt.0) then
        !!#!    write(*,*) " *** Error reading data"
        !!#!    return
        !!#!endif

        !!#!! Close dataset
        !!#!call h5dclose_f(ds_id, error)
        !!#!if(error /= 0) then
        !!#!  write(*,*) 'FATAL ERROR: could not close dataset'
        !!#!  return
        !!#!endif
        !!#!write(*,*) 'Dataset closed'


        ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
        ! RETURN MAIN_DIMS AND H52DDOUBLEDATASET TO MAIN FUNCTION
        ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

        !!! Close group
        !!call h5gclose_f(gr_id, error)
        !!if(error /= 0) then
        !!  write(*,*) 'FATAL ERROR: could not close group'
        !!  return
        !!endif
        !!write(*,*) 'Group closed'

        ! Insert this new data into the grid
        call grid_raw_data(errout,trim(data_path),grids,i_counts,i_size,&
                lat_gridder,lat_thresh)

        !!#!call read_shawn_file(io7,errout,trim(data_path)//total_file_name,grids,i_counts,&
        !!#!                     i_size,lat_gridder,lat_thresh)

        ! Deallocate all the arrays for the next pass
        ! -------------------------------------------
        call clear_arrays

      endif
    enddo file_loop  
  endif

  deallocate(grids)
  deallocate(i_counts)
  deallocate(lat_range)
  deallocate(lon_range)
  close(io8)
  close(io6)
  close(errout)  
  

  ! Close the HDF5 interface
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program omi_frequency_JZ
