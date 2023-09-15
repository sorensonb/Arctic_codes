program omi_JZ_daily_climo
!
! NAME:
!   omi_JZ_daily_climo.f90
!
! PURPOSE:
!   Calculate daily averages of screened AI (following the JZ criteria) and
!   write the averages to an output file.
! 
! CALLS:
!   Modules:
!     - hdf5
!     - h5_vars (custom module, shared between count and climo analyses)
!   Subroutines:
!     - check_bad_rows
!     - clear_h5_arrays (from h5_vars)
!     - grid_raw_data_daily_climo
!     - print_climo
!     - read_h5_AI
!     - read_h5_LAT
!     - read_h5_LON
!     - read_h5_XTRACK
!     - read_h5_AZM
!     - read_h5_GPQF
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/09/11:
!     Written
!
!  ############################################################################

  use hdf5
  use h5_vars, only : clear_h5_arrays, i_bad_list
  use daily_vars, only: grid_AI, count_AI, day_values, lat_values, &
        lon_values, clear_daily_arrays

  implicit none

  integer                :: ii            ! loop counter
  integer                :: i_size        ! array size
  integer                :: i_day_idx
  integer                :: int_day       ! integer variable for day 
  integer                :: int_dtg       ! integer variable for day 
  integer                :: arg_count     ! Number of arguments passed to exec
  integer                :: work_month    ! currently analyzed month
  integer                :: work_day 
  integer                :: error         ! error flag
  integer                :: istatus
  integer                :: file_id       ! id for current HDF5 file
  integer                :: i_num_days
  integer                :: i_num_swath

  real                   :: lat_thresh    ! latitude threshold. Only analyze
                                          ! data north of this value.
  real                   :: lat_gridder   ! Used for setting up the latitude
                                          ! grid

  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer,parameter      :: io6    = 1827 ! Data output file
  integer,parameter      :: errout = 9    ! Error file
  integer,parameter      :: io10   = 1066 ! Row anomaly file

  character(len = 255)   :: out_file_start  ! output file name
  character(len = 255)   :: out_file_name   ! output file name
  character(len = 255)   :: file_name_file  ! name of file containing the 
                                            ! list of HDF5 file names to 
                                            ! be analyzed
  character(len = 255)   :: total_file_name ! file name read from each line 
                                            ! of file_name_file
  character(len = 8)     :: c_work_dtg      ! holds the previous year

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./omi_exec out_file_start file_name_file'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1,out_file_start)
  call get_command_argument(2,file_name_file)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"

  ! Initialize the working month to -1
  ! ----------------------------------
  work_month = -1

  ! Call subroutine to figure out how many days are in the file name file
  ! ---------------------------------------------------------------------
  call count_num_days(file_name_file, i_num_swath, i_num_days)
  write(*,*) 'Num days = ', i_num_days

  ! Allocate grid arrays
  ! --------------------
  lat_thresh = 65.
  lat_gridder = lat_thresh
  i_size = (90. - lat_thresh)
  allocate(lat_values(i_size))
  do ii=1,i_size
    lat_values(ii) = lat_thresh + (ii-1) + 0.5
  enddo

  allocate(lon_values(360))
  do ii=1,360
    lon_values(ii) = -180. + (ii-1) + 0.5
  enddo

  allocate(day_values(i_num_days))
  allocate(grid_AI(360, i_size, i_num_days))
  allocate(count_AI(360, i_size, i_num_days))

  ! Initialize grid arrays and set to 0 initially
  ! ----------------------------------------------
  day_values(:)   = 0
  grid_AI(:,:,:)  = 0.
  count_AI(:,:,:) = 0 

  ! Loop over each of the file names given in the file name file
  ! ------------------------------------------------------------

  ! open debug file
  ! ---------------
  open(errout, file = "omi_jz_error_climo.txt", iostat = istatus)
  if(istatus /= 0) then
    write(*,*) "error opening error file."
  endif

  ! open row anomaly file
  ! ---------------
  open(io10, file = "/home/bsorenson/OMI/"&
    //"row_anomaly_dates_20050401_20201001.txt", iostat = istatus)
  if(istatus /= 0) then
    write(*,*) "error opening row file."
  endif

  ! Initialize the work_day value to -1
  ! -----------------------------------
  work_day = -1

  i_day_idx = 0

  ! Open the file name file
  ! -----------------------
  open(io8, file = trim(file_name_file), iostat = istatus)
  if(istatus > 0) then
    write(errout,*) "ERROR: Problem reading "//trim(file_name_file)
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
        exit
      else if(istatus > 0) then
        write(errout,*) "ERROR: problem reading total_file_name"
        cycle file_loop
      else

        ! Extract day information from file name
        ! --------------------------------------
        read(total_file_name(51:52), *) int_day
        c_work_dtg = total_file_name(44:47)//total_file_name(49:52)
        read(c_work_dtg, *) int_dtg

        ! If the day of the new file is different than the current working
        ! day, call check_bad_row and update the bad row list
        ! ----------------------------------------------------------------
        if(work_day /= int_day) then
          call check_bad_rows(total_file_name,errout,io10)
          work_day = int_day
          i_day_idx = i_day_idx + 1
          day_values(i_day_idx) = int_dtg
          write(*,*) trim(total_file_name)
        endif

        !!#!! Extract month information from dtg
        !!#!! ----------------------------------
        !!#!read(total_file_name(49:50), *) int_month

        !!#!! If the month of the new file is greater than the current working
        !!#!! month, call print_climo and print the grid values to the output
        !!#!! file.
        !!#!! ----------------------------------------------------------------
        !!#!if(work_month == -1) then
        !!#!  work_month = int_month
        !!#!  c_work_year = total_file_name(44:47)  
        !!#!else if(work_month /= int_month) then
        !!#!  call print_climo(io6,grids,i_counts,i_size,c_work_year,&
        !!#!                   work_month,lat_range,lon_range)
        !!#!  work_month = int_month
        !!#!  c_work_year = total_file_name(44:47)  
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
        !!#!call read_h5_CLD(file_id) 
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
        call grid_raw_data_daily_climo(i_day_idx, i_size,lat_gridder, &
                lat_thresh)

        ! Deallocate all the arrays for the next pass
        ! -------------------------------------------
        call clear_h5_arrays
      endif
    enddo file_loop  
  endif


  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! DATA OUTPUT SECTION
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Average the grid values 
  ! -----------------------
  call calc_grid_avgs(i_num_days, i_size, 360)

  !!#!if(resolution >= 1.0) then
  !!#!  name_format = '(a19,i3,a5)'
  !!#!  write(out_file_name,name_format) 'grid_coloc_test_res',int(resolution*100), '.hdf5'
  !!#!else 
  !!#!  name_format = '(a20,i2,a5)'
  !!#!  write(out_file_name,name_format) 'grid_coloc_test_res0',int(resolution*100), '.hdf5'
  !!#!endif

  out_file_name = trim(out_file_start)//'.hdf5'
  call write_output_file(out_file_name, i_num_days, i_size, 360)

  ! If, for some reason, the list of bad rows was not deallocated from
  ! before, deallocate it.
  ! ------------------------------------------------------------------
  if(allocated(i_bad_list)) deallocate(i_bad_list)

  ! Deallocate the remaining allocated arrays and close all files
  ! -------------------------------------------------------------
  close(io8)
  !!#!close(io6)
  close(io10)
  close(errout)  
  !!#!deallocate(grids)
  !!#!deallocate(i_counts)

  call clear_daily_arrays
  
  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program omi_JZ_daily_climo
