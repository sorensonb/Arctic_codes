program omi_shawn_climo
!
! NAME:
!   omi_shawn_climo.f90
!
! PURPOSE:
!   Calculate monthly averages of AI perturbations and
!   write the averages to an output file.
! 
! CALLS:
!   Subroutines:
!     - count_ai
!     - print_climo
!     - read_shawn_file
!     - synop_time_check
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2018/10/24:
!     Written
!
!  ############################################################################

  use hdf5
  use daily_vars, only: grid_AI, count_AI, day_values, lat_values, &
        lon_values, clear_daily_arrays, allocate_daily_arrays, &
        count_sfland, count_seaice, &
        count_permice, count_drysnow, count_ocean, count_mixpixel, &
        count_other   

  implicit none

  integer                :: ii            ! loop counter
  integer                :: i_size        ! array size
  integer                :: int_month     ! integer variable for month
  integer                :: arg_count     ! Number of arguments passed to exec
  integer                :: istatus       ! error flag
  integer                :: i_num_days
  integer                :: i_num_swath
  integer                :: i_day_idx
  integer                :: int_dtg
  integer                :: work_day
  integer                :: error

  real                   :: lat_thresh    ! latitude threshold. Only analyze
                                          ! data north of this value.
  real                   :: lat_gridder   ! Used for setting up the latitude
                                          ! grid
  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer,parameter      :: io7    = 22   ! File object for each data file 
  integer,parameter      :: errout = 9    ! File object for error file

  character(len = 255)   :: data_path       ! path to data files
  character(len = 12)    :: dtg             ! dtg from each line of file
  character(len = 255)   :: out_file_start  ! output file name
  character(len = 255)   :: date_file_name  ! name of file containing the 
                                            ! list of shawn file names to 
                                            ! be analyzed
  character(len = 255)   :: out_file_name   ! output file name
  character(len = 4)     :: c_work_year     ! holds the previous year

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./omi_exec out_file_start date_file_name'
    return
  endif

  call get_command_argument(1,out_file_start)
  call get_command_argument(2,date_file_name)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"

  data_path = "/Research/OMI/out_files-ltc4/"
  !data_path = "/Research/OMI/out_files-ltc3/new/"
  !data_path = "/Research/OMI/out_files-monthly.20210518/"

  ! Call subroutine to figure out how many days are in the file name file
  ! ---------------------------------------------------------------------
  call count_num_days_shawn(date_file_name, i_num_swath, i_num_days)
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

  ! Call routine to allocate the grid data arrays
  ! ---------------------------------------------
  call allocate_daily_arrays(360, i_size, i_num_days)

  ! open debug file
  ! ---------------
  open(errout, file = "omi_error_climo.txt", iostat = istatus)
  if(istatus /= 0) then
    write(*,*) "error opening error file."
  endif

  ! Initialize the work_day value to -1
  ! -----------------------------------
  work_day = -1

  i_day_idx = 0

  ! Read the file names from the file name file
  ! -------------------------------------------
  open(io8, file = trim(date_file_name), iostat = istatus)
  if(istatus > 0) then
    write(errout,*) "ERROR: Problem reading "//trim(date_file_name)
    return 
  else
    ! Loop over the file
    file_loop: do
      ! Read the current dtg from the file
      ! ----------------------------------
      read(io8, *, iostat=istatus) dtg
      if(istatus < 0) then 
        write(*,*) "End of "//trim(date_file_name)//" found"
        exit
      else if(istatus > 0) then
        write(errout,*) "ERROR: problem reading dtg: "//dtg
        cycle file_loop
      else

        ! Extract day information from file name
        ! --------------------------------------
        read(dtg(1:8), *) int_dtg

        ! If the day of the new file is different than the current working
        ! day, call check_bad_row and update the bad row list
        ! ----------------------------------------------------------------
        if(work_day /= int_dtg) then
          work_day = int_dtg
          i_day_idx = i_day_idx + 1
          day_values(i_day_idx) = int_dtg
          write(*,*) int_dtg
        endif

        ! Open the shawn file and look at contents
        open(io7, file = trim(data_path)//dtg, action = 'read', &
            iostat = istatus)
        if(istatus /= 0) then
          write(errout, *) "ERROR: error opening file",trim(data_path)//dtg
          write(errout, *) "       cycling file_loop"
          cycle file_loop
        endif

        ! Read data from the current file and insert into the grids
        ! ---------------------------------------------------------
        write(*,*) trim(data_path)//dtg
        call read_shawn_file_daily_climo(io7,errout,trim(data_path)//dtg,&
                i_day_idx, i_size,lat_gridder, lat_thresh)
        
        ! Close the current shawn file
        ! ----------------------------
        close(io7)
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

  out_file_name = trim(out_file_start)//'.hdf5'
  call write_output_file(out_file_name, i_num_days, i_size, 360)


  ! Deallocate the remaining allocated arrays and close all files
  ! -------------------------------------------------------------
  close(io8)
  close(errout)  
  
  call clear_daily_arrays

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program omi_shawn_climo
