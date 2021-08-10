program omi_JZ_drift
!
! NAME:
!   omi_JZ_drift.f90
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
!     - grid_raw_data_drift
!     - print_drift
!     - read_h5_AI
!     - read_h5_LAT
!     - read_h5_LON
!     - read_h5_XTRACK
!     - read_h5_AZM
!     - read_h5_GPQF
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/16:
!     Written
!
!  ############################################################################

  use hdf5
  use h5_vars, only : clear_arrays, i_bad_list

  implicit none

  integer                :: int_month     ! integer variable for month
  integer                :: int_day       ! integer variable for day 
  integer                :: arg_count     ! Number of arguments passed to exec
  integer                :: work_month    ! currently analyzed month
  integer                :: work_day 
  integer                :: error         ! error flag
  integer                :: istatus
  integer                :: file_id       ! id for current HDF5 file

  real                   :: avg_ai        ! quarter-degree AI grid
                                                   ! values.
  integer                :: i_count       ! quarter-degree AI counts

  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer,parameter      :: io6    = 1827 ! Data output file
  integer,parameter      :: errout = 9    ! Error file
  integer,parameter      :: io10   = 1066 ! Row anomaly file

  character(len = 255)   :: out_file_name   ! output file name
  character(len = 255)   :: file_name_file  ! name of file containing the 
                                            ! list of HDF5 file names to 
                                            ! be analyzed
  character(len = 255)   :: total_file_name ! file name read from each line 
                                            ! of file_name_file
  character(len = 4)     :: c_work_year     ! holds the previous year

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

  ! Initialize the working month and day to -1
  ! ------------------------------------------
  work_month = -1
  work_day = -1

  ! Initialize avg value / counts and set to 0 initially
  ! ----------------------------------------------------
  avg_ai     = 0.
  i_count    = 0

  ! open debug file
  ! ---------------
  open(errout, file = "omi_jz_error_climo.txt", iostat = istatus)
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
  open(io6, file = trim(out_file_name), iostat = istatus)
  if(istatus /= 0) then
    write(errout,*) "ERROR: error opening climo output file."
  endif
  write(io6,'(a4,2x,3(a8))') 'Date','Lat Lon','Avg','#_obs'

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
        ! Print the final month of data to the output file using
        ! print_climo
        ! --------------------------------------------------------
        call print_drift(io6,avg_ai,i_count,c_work_year,work_month)
        exit
      else if(istatus > 0) then
        write(errout,*) "ERROR: problem reading total_file_name"
        cycle file_loop
      else

        ! Extract day information from file name
        ! --------------------------------------
        read(total_file_name(51:52), *) int_day

        ! If the day of the new file is different than the current working
        ! day, call check_bad_row and update the bad row list
        ! ----------------------------------------------------------------
        if(work_day /= int_day) then
          call check_bad_rows(total_file_name,errout,io10)
          work_day = int_day
        endif

        ! Extract month information from dtg
        ! ----------------------------------
        read(total_file_name(49:50), *) int_month

        ! If the month of the new file is greater than the current working
        ! month, call print_climo and print the grid values to the output
        ! file.
        ! ----------------------------------------------------------------
        if(work_month == -1) then
          work_month = int_month
          c_work_year = total_file_name(44:47)  
        else if(work_month /= int_month) then
          call print_drift(io6,avg_ai,i_count,c_work_year,work_month)
          work_month = int_month
          c_work_year = total_file_name(44:47)  
        endif

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
        call grid_raw_data_drift(avg_ai,i_count)

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
  close(io8)
  close(io6)
  close(io10)
  close(errout)  
  
  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program omi_JZ_drift
