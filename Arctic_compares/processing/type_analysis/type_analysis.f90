program type_analysis
!
! NAME:
!   type_analysis.f90
!
! PURPOSE:
!
!   As of 2024/01/02: reads in the coloc files and generates daily
!   averages of AI, CH7, and NSIDC surface type.
! 
! CALLS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/09:
!     Written
!
!  ############################################################################

  use hdf5
  use type_vars, only: allocate_grid_arrays, clear_grid_arrays, &
    day_values, lat_range, lon_range, grid_areas

  implicit none

  integer                :: i_day_idx
  integer                :: i_lat_size        ! latitude array size
  integer                :: i_lon_size        ! longitude array size
  integer                :: int_whole_date    ! integer variable for day
  integer                :: int_day           ! integer variable for day
  integer                :: work_day 
  integer                :: arg_count     ! Number of arguments passed to exec
  integer                :: i_num_swath
  integer                :: i_num_days
  integer                :: error         ! error flag

  integer                :: test_day

  !real                   :: !uprev_lat 
  real                   :: lat_thresh    ! latitude threshold. Only analyze
                                          ! data north of this value.
  real                   :: lat_gridder
  real                   :: lon_gridder
 
  real                   :: resolution

  logical                :: l_debug

  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer,parameter      :: io7    = 22   ! File object for each data file 
  integer,parameter      :: io6    = 1827 ! Data output file
  integer,parameter      :: errout = 9 
  integer                :: istatus

  character(len = 120)   :: name_format
  character(len = 255)   :: out_file_name   
  character(len = 12)    :: single_file_date
  character(len = 12)    :: c_work_dtg        ! dtg from each line of file
  character(len = 255)   :: file_name_file  ! name of file containing the 
                                            ! list of shawn file names to 
                                            ! be analyzed

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  l_debug = .false.

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 1) then
    write(*,*) 'SYNTAX: ./type_exec file_name_file'
    return
  endif

  call get_command_argument(1,file_name_file)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! STEPS: 
  ! - Declare grid arrays to hold the daily averaged variables, dimensioned
  !       to:
  !     - # days
  !     - # lat grid boxes
  !     - # lon grid boxes
  ! - Loop over each of the given colocation files 
  ! - Grid all colocated variables (mainly OMI/TROP AI, COD, CLDFRAC, 
  !       NSIDC ICE, etc.) into a 0.25 x 0.25 / 0.5 x 0.5 / 1.0 x 1.0 deg.
  !       lat/lon grid for that day.
  ! - Once the new date string has a new day than before, increment 
  !       the date idx for the arrays? Average the data first?
  ! - Write the gridded arrays and areas to an HDF5 file
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


  ! Set up lat/lon grids
  ! --------------------
  resolution = 0.50
  lat_thresh = 70.0
  call fill_lat_lon_arrays(resolution, lat_thresh, i_lat_size, i_lon_size, &
                           lat_gridder, lon_gridder)

  write(*,*) lat_range

  ! Count the number of swaths and days in the file
  ! ----------------------------------------------- 
  call count_num_days(file_name_file, i_num_swath, i_num_days)

  write(*,*) 'In main:', i_num_swath, i_num_days, lat_gridder, lon_gridder

  allocate(day_values(i_num_days))
  day_values(:) = -99

  ! Initialize grid arrays and set to 0 initially
  ! ---------------------------------------------
  call allocate_grid_arrays(i_num_days, i_lat_size, i_lon_size)
 
  work_day = -1
  i_day_idx = 1 

  ! Read the file names from the file name file
  open(io8, file = trim(file_name_file), iostat = istatus)
  if(istatus > 0) then
    write(errout,*) "ERROR: Problem reading "//trim(file_name_file)
    return 
  else
    ! Loop over the dtg file
    ! ----------------------------
    file_loop: do
      ! Read the current dtg from the file
      ! ----------------------------------
      read(io8, *, iostat=istatus) single_file_date
      if(istatus < 0) then 
        write(*,*) "End of "//trim(file_name_file)//" found"
        exit
      else if(istatus > 0) then
        write(*,*) "ERROR: problem reading single_file_date"
        cycle file_loop
      else
        write(*,*) single_file_date

        ! Extract day information from file name
        ! --------------------------------------
        read(single_file_date(7:8), *) int_day

        ! If the day of the new file is different than the current working
        ! day, call count_ai and find the counts in the daily averages
        ! ----------------------------------------------------------------
        if(work_day == -1) then
          work_day = int_day
          c_work_dtg = single_file_date
        else if(work_day /= int_day) then
          work_day = int_day
          c_work_dtg = single_file_date
          i_day_idx = i_day_idx + 1
        endif

        read(single_file_date(1:8), '(i8)') day_values(i_day_idx)

        ! Read data from the current file and insert into the grids
        ! ---------------------------------------------------------
        call grid_comp_data(single_file_date, lat_thresh, &
                           i_day_idx, l_debug)

      endif
    enddo file_loop  
  endif

  ! Loop over the grid arrays and calculate averages
  ! ------------------------------------------------
  call average_final_grid_data(i_num_days, i_lat_size, i_lon_size)

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! DATA OUTPUT SECTION
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  if(resolution >= 1.0) then
    name_format = '(a19,i3,a5)'
    write(out_file_name,name_format) 'grid_coloc_test_res',int(resolution*100), '.hdf5'
  else 
    name_format = '(a20,i2,a5)'
    write(out_file_name,name_format) 'grid_coloc_test_res0',int(resolution*100), '.hdf5'
  endif

  call write_output_file(out_file_name, i_num_days, i_lat_size, i_lon_size)

  ! Write the data to the output HDF5 file
  ! --------------------------------------

  ! Deallocate the remaining allocated arrays and close all files
  ! -------------------------------------------------------------
  call clear_grid_arrays

  write(*,*) day_values
  
  deallocate(day_values)
  deallocate(lat_range)
  deallocate(lon_range)
  deallocate(grid_areas)

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"
  
end program type_analysis
