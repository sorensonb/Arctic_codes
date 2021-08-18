program omi_frequency_JZ
!
! NAME:
!   omi_frequency_JZ.f90
!
! PURPOSE:
!   Calculate counts of high-AI (w.r.t. a user-defined threshold) quarter
!   degree grid-boxes for screened AI data averaged into 4 +/- 3hr windows 
!   (00, 06, 12, and 18). Since these codes read AI data directly from the 
!   HDF5 files, the hdf5 module is required.
! 
! CALLS:
!   Modules:
!     - hdf5
!     - h5_vars (custom module, shared between count and climo analyses)
!   Subroutines:
!     - check_bad_rows
!     - clear_arrays (from h5_vars)
!     - count_ai_JZ
!     - grid_raw_data
!     - read_h5_AI
!     - read_h5_LAT
!     - read_h5_LON
!     - read_h5_XTRACK
!     - read_h5_AZM
!     - read_h5_GPQF
!     - synop_time_check
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/09:
!     Written
!
!  ############################################################################

  use hdf5
  use h5_vars, only : clear_arrays, i_bad_list
 
  implicit none

  integer                :: ii            ! loop counter
  integer                :: synop_idx     ! index of current synoptic time
  integer                :: i_size        ! array size
  integer                :: int_hr        ! integer variable for hour
  integer                :: int_day       ! integer variable for day 
  integer                :: work_day      ! currently analyzed day
  integer                :: arg_count     ! Number of arguments passed to exec
  integer                :: error         ! error flag
  integer                :: istatus
  integer                :: file_id       ! id for current HDF5 file

  real                   :: ai_thresh     ! threshold AI value
  real                   :: lat_thresh    ! latitude threshold. Only analyze
                                          ! data north of this value.
  real                   :: lat_gridder   ! Used for setting up the latitude
                                          ! grid

  logical                :: l_in_time     ! used to determine if the current
                                          ! hour falls within the synoptic time
                                          ! range.

  real,dimension(:), allocatable      :: lat_range ! latitude grid
  real,dimension(:), allocatable      :: lon_range ! longitude grid
  real,dimension(:,:), allocatable    :: grids     ! quarter-degree AI grid
                                                   ! values.
  integer,dimension(:,:), allocatable :: i_counts  ! quarter-degree AI counts

  integer,dimension(4)   :: synop_times     ! list containing the 4 synoptic
                                            ! times

  character(len = 10)    :: c_work_dtg      ! dtg from each line of file
  character(len = 255)   :: out_file_name   ! output file name
  character(len = 255)   :: file_name_file  ! name of file containing the 
                                            ! list of HDF5 file names to 
                                            ! be analyzed
  character(len = 255)   :: total_file_name ! file name read from each line 
                                            ! of file_name_file

  ! File read variables
  integer,parameter      :: io8    = 1042 ! File object for file name file
  integer,parameter      :: io6    = 1827 ! Data output file
  integer,parameter      :: errout = 1009 ! Error file
  integer,parameter      :: io10   = 1066 ! Row anomaly file

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./omi_exec out_file_name file_name_file'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1,out_file_name)
  call get_command_argument(2,file_name_file)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"

  ! Set up the synoptic time list
  ! -----------------------------
  synop_times = [0,6,12,18] 
  synop_idx = 1

  ! Set up lat/lon grids
  ! --------------------
  lat_thresh = 60.
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
  write(io6,'(a10,6(a6))') 'Date','Cnt60','Cnt65','Cnt70','Cnt75','Cnt80',&
    'Cnt85'

  ! Initialize the work_day value to -1
  ! -----------------------------------
  work_day = -1
 
  ! Set up count variables to count the number of grid boxes with
  ! high AI values
  ! -------------------------------------------------------------
  ai_thresh = 0.6
  
  ! Open the file name file
  ! -----------------------
  open(io8, file = trim(file_name_file), iostat = istatus)
  if(istatus > 0) then
    write(errout,*) "ERROR: Problem opening file name containing files: "&
        //'titled '//trim(file_name_file)
    return 
  else
    ! Loop over the file name file
    ! ----------------------------
    file_loop: do
      ! Read the current total_file_name from the file
      ! ----------------------------------------------
      read(io8, '(A)', iostat=istatus) total_file_name
      if(istatus < 0) then 
        write(*,*) "End of "//trim(file_name_file)//" found"
        call count_ai_JZ(io6,grids,i_counts,i_size,ai_thresh,synop_idx,&
                      total_file_name,lat_range)
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
          call check_bad_rows(total_file_name,errout,io10)
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
        !!#!if(work_day == -1) then
        !!#!  call check_bad_rows(total_file_name,errout,io10)
        !!#!  work_day = int_day
        !!#!  c_work_dtg = total_file_name(44:47)//total_file_name(49:52)
        !!#!else if(work_day /= int_day) then
        !!#!  call count_ai_JZ(io6,grids,i_counts,i_size,ai_thresh,synop_idx,&
        !!#!                c_work_dtg,lat_range)
        !!#!  call check_bad_rows(total_file_name,errout,io10)
        !!#!  work_day = int_day
        !!#!  c_work_dtg = total_file_name(44:47)//total_file_name(49:52)
        !!#!endif

        ! Open the HDF5 file
        ! ------------------
        call h5fopen_f(total_file_name, H5F_ACC_RDWR_F, file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: could not open file'
          return
        endif

        ! Read in the necessary data using the read routines
        ! --------------------------------------------------
        call read_h5_AI(file_id)
        call read_h5_LAT(file_id)
        call read_h5_LON(file_id)
        call read_h5_XTRACK(file_id)   
        call read_h5_AZM(file_id) 
        call read_h5_GPQF(file_id)

        ! Close file
        ! ----------
        call h5fclose_f(file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: could not close file'
          return
        endif

        ! Insert this new data into the grid 
        ! ----------------------------------
        call grid_raw_data(grids,i_counts,i_size,lat_gridder,lat_thresh)

        ! Deallocate all the arrays for the next pass
        ! -------------------------------------------
        call clear_arrays

      endif
    enddo file_loop  
  endif

  ! If, for some reason, the list of bad rows was not deallocated from
  ! before, deallocate it.
  ! ------------------------------------------------------------------
  if(allocated(i_bad_list)) deallocate(i_bad_list)

  ! Deallocate the remaining allocated arrays and close all files
  ! -------------------------------------------------------------
  deallocate(grids)
  deallocate(i_counts)
  deallocate(lat_range)
  deallocate(lon_range)
  close(io10)
  close(io8)
  close(io6)
  close(errout)  
  

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program omi_frequency_JZ
