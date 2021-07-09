program omi_frequency_JZ
!
! NAME:
!   omi_frequency_JZ.f90
!
! PURPOSE:
! 
! CALLS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2018/10/24:
!     Written
!
!  ############################################################################

  use hdf5
  use h5_vars, only : clear_arrays, i_bad_list
 
  implicit none

  integer                :: ii            ! loop counter
  integer                :: synop_idx
  integer                :: i_size        ! array size
  integer                :: ai_count      ! good AI counter
  integer,dimension(3)   :: ai_count2     ! good AI counter
  integer                :: int_hr        ! integer variable for hour
  integer                :: int_day       ! integer variable for day 

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
  integer,parameter      :: io10   = 1066
  integer                :: istatus

  character(len = 255)   :: data_path
  character(len = 255)   :: out_file_name
  character(len = 255)   :: file_name_file
  character(len = 255)   :: total_file_name

  integer                :: arg_count
  integer                :: work_day 

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


  ! Initialize grid arrays and set to 0 initially
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

  ! open row anomaly file
  ! ---------------
  open(io10, file = "/home/bsorenson/OMI/"&
    //"row_anomaly_dates_20050401_20191001.txt", iostat = istatus)
  if(istatus /= 0) then
    write(*,*) "error opening row file."
  endif

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

        ! Extract day information from file name
        ! --------------------------------------
        read(total_file_name(51:52), *) int_day

        ! Extract time information from total_file_name
        ! ---------------------------------
        read(total_file_name(54:55), *) int_hr

        ! If the day of the new file is different than the current working
        ! day, call check_bad_row and update the bad row list
        ! ----------------------------------------------------------------
        if(work_day /= int_day) then
          call check_bad_rows(errout,io10)
          work_day = int_day
        endif

        ! See if the hour exceeds the current 6 hr assimilation window.
        ! If so, calculate averages and counts and reset variables.
        ! ------------------------------------------------------------
        call synop_time_check(synop_idx, int_hr, l_in_time)
        if(.not. l_in_time) then 
          call count_ai_JZ(io6,grids,i_counts,i_size,ai_thresh,synop_idx,&
                        total_file_name,lat_range)

        endif  

        !write(*,*) trim(total_file_name)

        ! Open the HDF5 file
        call h5fopen_f(total_file_name, H5F_ACC_RDWR_F, file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: could not open file'
          return
        endif
        !write(*,*) 'File opened'

        ! Read in the necessary data using the read routines
        ! --------------------------------------------------
        call read_h5_AI(file_id,error)
        call read_h5_LAT(file_id,error)
        call read_h5_LON(file_id,error)
        call read_h5_XTRACK(file_id,error)   
        call read_h5_AZM(file_id,error) 
        call read_h5_GPQF(file_id,error)

        ! Close file
        call h5fclose_f(file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: could not close file'
          return
        endif
        !write(*,*) 'File closed'


        ! Insert this new data into the grid 
        ! ----------------------------------
        call grid_raw_data(errout,grids,i_counts,i_size,&
                lat_gridder,lat_thresh)

        ! Deallocate all the arrays for the next pass
        ! -------------------------------------------
        call clear_arrays

      endif
    enddo file_loop  
  endif

  if(allocated(i_bad_list)) deallocate(i_bad_list)

  deallocate(grids)
  deallocate(i_counts)
  deallocate(lat_range)
  deallocate(lon_range)
  close(io10)
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
