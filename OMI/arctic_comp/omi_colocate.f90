program omi_colocate
!
! NAME:
!   omi_colocate.f90
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
  use comp_vars, only : clear_arrays, i_bad_list, &
    MODIS_CH2_data, MODIS_CH2_dims, &
    MODIS_LAT_data, MODIS_LAT_dims, &
    MODIS_LON_data, MODIS_LON_dims, &
    NSIDC_data,     NSIDC_dims, &
    NSIDC_LAT_data, NSIDC_LAT_dims, &
    NSIDC_LON_data, NSIDC_LON_dims, &
    CERES_LWF_data, CERES_LWF_dims, &
    CERES_SWF_data, CERES_SWF_dims, &
    CERES_LAT_data, CERES_LAT_dims, &
    CERES_LON_data, CERES_LON_dims, &
    OMI_AI_data,    OMI_AI_dims, &
    OMI_LAT_data,   OMI_LAT_dims, &
    OMI_LON_data,   OMI_LON_dims

  implicit none

  integer                :: ii            ! loop counter
  integer                :: jj            ! loop counter
  integer                :: nii            ! loop counter
  integer                :: njj            ! loop counter
  !!#!integer                :: i_size        ! array size
  !!#!integer                :: int_month     ! integer variable for month
  !!#!integer                :: int_day       ! integer variable for day 
  integer                :: arg_count     ! Number of arguments passed to exec
  !!#!integer                :: work_month    ! currently analyzed month
  !!#!integer                :: work_day 
  integer                :: error         ! error flag
  integer                :: istatus
  integer                :: modis_file_id       ! id for current HDF5 file
  integer                :: nsidc_file_id       ! id for current HDF5 file
  integer                :: ceres_file_id       ! id for current HDF5 file
  integer                :: omi_file_id         ! id for current HDF5 file

  !!#!real                   :: lat_thresh    ! latitude threshold. Only analyze
  !!#!                                        ! data north of this value.
  !!#!real                   :: lat_gridder   ! Used for setting up the latitude
  !!#!                                        ! grid

  !!#!real,dimension(:), allocatable      :: lat_range ! latitude grid
  !!#!real,dimension(:), allocatable      :: lon_range ! longitude grid
  !!#!real,dimension(:,:), allocatable    :: grids     ! quarter-degree AI grid
  !!#!                                                 ! values.
  !!#!integer,dimension(:,:), allocatable :: i_counts  ! quarter-degree AI counts

  !!#!! File read variables
  !!#!integer,parameter      :: io8    = 42   ! File object for file name file
  !!#!integer,parameter      :: io6    = 1827 ! Data output file
  !!#!integer,parameter      :: errout = 9    ! Error file
  !!#!integer,parameter      :: io10   = 1066 ! Row anomaly file

  !!#!character(len = 255)   :: out_file_name   ! output file name
  !!#!character(len = 255)   :: file_name_file  ! name of file containing the 
  !!#!                                          ! list of HDF5 file names to 
  !!#!                                          ! be analyzed
  !!#!character(len = 255)   :: total_file_name ! file name read from each line 
  !!#!                                          ! of file_name_file
  !!#!character(len = 4)     :: c_work_year     ! holds the previous year

  character(len = 255)      :: modis_name     ! filename
  character(len = 255)      :: nsidc_name     ! filename
  character(len = 255)      :: ceres_name     ! filename
  character(len = 255)      :: omi_name       ! filename

  !!#!! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  write(*,*) 'yayyy'

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 4) then
    write(*,*) 'SYNTAX: ./omi_exec modis_name nsidc_name ceres_name omi_name'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1,modis_name)
  call get_command_argument(2,nsidc_name)
  call get_command_argument(3,ceres_name)
  call get_command_argument(4,omi_name)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"

  write(*,*) trim(modis_name), trim(nsidc_name)

  ! Open the HDF5 file
  ! ------------------
  call h5fopen_f(trim(modis_name), H5F_ACC_RDWR_F, modis_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open file'
    return
  endif

  ! Open the HDF5 file
  ! ------------------
  call h5fopen_f(trim(nsidc_name), H5F_ACC_RDWR_F, nsidc_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open file'
    return
  endif

  ! Open the HDF5 file
  ! ------------------
  call h5fopen_f(trim(ceres_name), H5F_ACC_RDWR_F, ceres_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open file'
    return
  endif

  ! Open the HDF5 file
  ! ------------------
  call h5fopen_f(trim(omi_name), H5F_ACC_RDWR_F, omi_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open file'
    return
  endif

  call read_comp_MODIS_LAT(modis_file_id)
  call read_comp_MODIS_LON(modis_file_id)
  call read_comp_NSIDC_LAT(nsidc_file_id)
  call read_comp_NSIDC_LON(nsidc_file_id)
  call read_comp_CERES_LAT(ceres_file_id)
  call read_comp_CERES_LON(ceres_file_id)
  call read_comp_OMI_LAT(omi_file_id)
  call read_comp_OMI_LON(omi_file_id)

  write(*,*) MODIS_LAT_data(100,100), MODIS_LON_data(100,100)

  istatus = 0
  modis_loop1: do ii=1,MODIS_LAT_dims(2)
    modis_loop2: do jj=1,MODIS_LAT_dims(1) 
      ! Now, loop over the NSIDC data
      ! -----------------------------
      nsidc_loop1: do nii=1,NSIDC_LAT_dims(2)
        nsidc_loop2: do njj=1,NSIDC_LAT_dims(1) 

          !write(*,*) LAT_data(jj,ii), LON_data(jj,ii)
          istatus = 3
        enddo nsidc_loop2
      enddo nsidc_loop1

      ! Now, loop over the CERES data
      ! -----------------------------
      ceres_loop1: do nii=1,CERES_LAT_dims(2)
        ceres_loop2: do njj=1,CERES_LAT_dims(1) 

          !write(*,*) LAT_data(jj,ii), LON_data(jj,ii)
          istatus = istatus + 1
        enddo ceres_loop2
      enddo ceres_loop1

    enddo modis_loop2
  enddo modis_loop1
   
  write(*,*) istatus
 
  ! Deallocate all the arrays for the next pass
  ! -------------------------------------------
  call clear_arrays

  ! Close file
  ! ----------
  call h5fclose_f(modis_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close file'
    return
  endif

  call h5fclose_f(nsidc_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close file'
    return
  endif

  call h5fclose_f(ceres_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close file'
    return
  endif

  call h5fclose_f(omi_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close file'
    return
  endif

  !!#!! Initialize the working month to -1
  !!#!! ----------------------------------
  !!#!work_month = -1

  !!#!! Set up lat/lon grids
  !!#!! --------------------
  !!#!lat_thresh = 45.
  !!#!lat_gridder = lat_thresh
  !!#!i_size = (90. - lat_thresh)
  !!#!allocate(lat_range(i_size))
  !!#!do ii=1,i_size
  !!#!  lat_range(ii) = lat_thresh + (ii-1)
  !!#!enddo

  !!#!allocate(lon_range(360))
  !!#!do ii=1,360
  !!#!  lon_range(ii) = -180.0 + (ii-1)
  !!#!enddo

  !!#!! Initialize grid arrays and set to -9 initially
  !!#!! ----------------------------------------------

  !!#!allocate(grids(360,i_size))
  !!#!allocate(i_counts(360,i_size))

  !!#!grids(:,:) = 0.
  !!#!i_counts(:,:) = 0

  !!#!! open debug file
  !!#!! ---------------
  !!#!open(errout, file = "omi_jz_error_climo.txt", iostat = istatus)
  !!#!if(istatus /= 0) then
  !!#!  write(*,*) "error opening error file."
  !!#!endif

  !!#!! open row anomaly file
  !!#!! ---------------
  !!#!open(io10, file = "/home/bsorenson/OMI/"&
  !!#!  //"row_anomaly_dates_20050401_20201001.txt", iostat = istatus)
  !!#!if(istatus /= 0) then
  !!#!  write(*,*) "error opening row file."
  !!#!endif

  !!#!! open output file
  !!#!! ---------------
  !!#!open(io6, file = trim(out_file_name), iostat = istatus)
  !!#!if(istatus /= 0) then
  !!#!  write(errout,*) "ERROR: error opening climo output file."
  !!#!endif
  !!#!write(io6,'(a4,2x,3(a7))') 'Date','Lat Lon','Avg','#_obs'

  !!#!! Initialize the work_day value to -1
  !!#!! -----------------------------------
  !!#!work_day = -1

  !!#!! Open the file name file
  !!#!! -----------------------
  !!#!open(io8, file = trim(file_name_file), iostat = istatus)
  !!#!if(istatus > 0) then
  !!#!  write(errout,*) "ERROR: Problem reading "//trim(file_name_file)
  !!#!  return 
  !!#!else
  !!#!  ! Loop over the file name file
  !!#!  ! ----------------------------
  !!#!  file_loop: do
  !!#!    ! Read the current total_file_name from the file
  !!#!    ! ----------------------------------------------
  !!#!    read(io8, '(A)', iostat=istatus) total_file_name
  !!#!    if(istatus < 0) then 
  !!#!      write(*,*) "End of "//trim(file_name_file)//" found"
  !!#!      ! Print the final month of data to the output file using
  !!#!      ! print_climo
  !!#!      ! --------------------------------------------------------
  !!#!      call print_climo(io6,grids,i_counts,i_size,c_work_year,&
  !!#!                       work_month,lat_range,lon_range)
  !!#!      exit
  !!#!    else if(istatus > 0) then
  !!#!      write(errout,*) "ERROR: problem reading total_file_name"
  !!#!      cycle file_loop
  !!#!    else

  !!#!      ! Extract day information from file name
  !!#!      ! --------------------------------------
  !!#!      read(total_file_name(51:52), *) int_day

  !!#!      ! If the day of the new file is different than the current working
  !!#!      ! day, call check_bad_row and update the bad row list
  !!#!      ! ----------------------------------------------------------------
  !!#!      if(work_day /= int_day) then
  !!#!        call check_bad_rows(total_file_name,errout,io10)
  !!#!        work_day = int_day
  !!#!      endif

  !!#!      ! Extract month information from dtg
  !!#!      ! ----------------------------------
  !!#!      read(total_file_name(49:50), *) int_month

  !!#!      ! If the month of the new file is greater than the current working
  !!#!      ! month, call print_climo and print the grid values to the output
  !!#!      ! file.
  !!#!      ! ----------------------------------------------------------------
  !!#!      if(work_month == -1) then
  !!#!        work_month = int_month
  !!#!        c_work_year = total_file_name(44:47)  
  !!#!      else if(work_month /= int_month) then
  !!#!        call print_climo(io6,grids,i_counts,i_size,c_work_year,&
  !!#!                         work_month,lat_range,lon_range)
  !!#!        work_month = int_month
  !!#!        c_work_year = total_file_name(44:47)  
  !!#!      endif

  !!#!      ! Open the HDF5 file
  !!#!      ! ------------------
  !!#!      call h5fopen_f(total_file_name, H5F_ACC_RDWR_F, file_id, error)
  !!#!      if(error /= 0) then
  !!#!        write(*,*) 'FATAL ERROR: could not open file'
  !!#!        return
  !!#!      endif

  !!#!      ! Read in the necessary data using the read routines
  !!#!      ! --------------------------------------------------
  !!#!      call read_h5_AI(file_id)
  !!#!      call read_h5_LAT(file_id)
  !!#!      call read_h5_LON(file_id)
  !!#!      call read_h5_XTRACK(file_id)   
  !!#!      call read_h5_AZM(file_id) 
  !!#!      call read_h5_GPQF(file_id)

  !!#!      ! Close file
  !!#!      ! ----------
  !!#!      call h5fclose_f(file_id, error)
  !!#!      if(error /= 0) then
  !!#!        write(*,*) 'FATAL ERROR: could not close file'
  !!#!        return
  !!#!      endif

  !!#!      ! Insert this new data into the grid 
  !!#!      ! ----------------------------------
  !!#!      call grid_raw_data_climo(grids,i_counts,i_size,&
  !!#!              lat_gridder,lat_thresh)

  !!#!      ! Deallocate all the arrays for the next pass
  !!#!      ! -------------------------------------------
  !!#!      call clear_arrays
  !!#!    endif
  !!#!  enddo file_loop  
  !!#!endif

  !!#!! If, for some reason, the list of bad rows was not deallocated from
  !!#!! before, deallocate it.
  !!#!! ------------------------------------------------------------------
  !!#!if(allocated(i_bad_list)) deallocate(i_bad_list)

  !!#!! Deallocate the remaining allocated arrays and close all files
  !!#!! -------------------------------------------------------------
  !!#!close(io8)
  !!#!close(io6)
  !!#!close(io10)
  !!#!close(errout)  
  !!#!deallocate(grids)
  !!#!deallocate(i_counts)
  !!#!deallocate(lat_range)
  !!#!deallocate(lon_range)
  !!#!
  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program omi_colocate
