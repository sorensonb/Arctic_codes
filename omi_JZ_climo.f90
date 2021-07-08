program omi_JZ_climo
!
! NAME:
!   omi_JZ_climo.f90
!
! PURPOSE:
! 
! CALLS:
!   mie_calc.f90
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2018/10/24:
!     Written
!
!  ############################################################################

  use hdf5
  use h5_vars, only : clear_arrays, temp_switch

  implicit none

  integer                :: ii            ! loop counter
  integer                :: jj            ! loop counter
  integer                :: index1        
  integer                :: index2        
  integer                :: i_size        ! array size
  integer                :: ai_count      ! good AI counter
  integer                :: int_month     ! integer variable for month

  real                   :: ai_thresh     ! threshold AI value
  real                   :: lat_thresh
  real                   :: lat_gridder
  real                   :: avg_ai

  logical                :: l_in_time

  real,dimension(:), allocatable :: lat_range
  real,dimension(:), allocatable :: lon_range   
  real,dimension(:,:), allocatable :: grids
  integer,dimension(:,:), allocatable :: i_counts

  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer,parameter      :: io7    = 22   ! File object for each data file 
  integer,parameter      :: io6    = 1827 ! Data output file
  integer,parameter      :: errout = 9 
  integer                :: istatus

  character(len = 255)   :: data_path
  character(len = 255)   :: out_file_name
  character(len = 255)   :: file_name_file
  character(len = 255)   :: total_file_name

  integer                :: arg_count
  integer                :: work_month 

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./omi_exec out_file_name date_file_name'
    return
  endif

  call get_command_argument(1,out_file_name)
  call get_command_argument(2,date_file_name)

  work_month = -1

  !data_path = "/home/bsorenson/OMI/JZ_analysis/test_dir/"
  data_path = "/Research/OMI/H5_files/"

  ! Set up lat/lon grids
  ! --------------------
  lat_thresh = 65.
  lat_gridder = lat_thresh
  i_size = (90. - lat_thresh)
  allocate(lat_range(i_size))
  do ii=1,i_size
    lat_range(ii) = lat_thresh + (ii-1)
  enddo

  allocate(lon_range(360))
  do ii=1,360
    lon_range(ii) = -180.0 + (ii-1)
  enddo

  ! Initialize grid arrays and set to -9 initially
  ! ----------------------------------------------

  allocate(grids(360,i_size))
  allocate(i_counts(360,i_size))

  grids(:,:) = 0.
  i_counts(:,:) = 0

  ! open debug file
  ! ---------------
  open(errout, file = "omi_jz_error_climo.txt", iostat = istatus)
  if(istatus /= 0) then
    write(*,*) "error opening error file."
  endif

  ! open output file
  ! ---------------
  open(io6, file = trim(out_file_name), iostat = istatus)
  if(istatus /= 0) then
    write(errout,*) "ERROR: error opening climo output file."
  endif
  write(io6,'(a4,2x,3(a7))') 'Date','Lat Lon','Avg','#_obs'

  ! Read the file names from the file name file
  open(io8, file = trim(date_file_name), iostat = istatus)
  if(istatus > 0) then
    write(errout,*) "ERROR: Problem reading "//trim(date_file_name)
    return 
  else
    ! Loop over the file
    file_loop: do
      ! Read the current dtg from the file
      read(io8, '(A)', iostat=istatus) total_file_name
      if(istatus < 0) then 
        write(*,*) "End of "//trim(date_file_name)//" found"
        call print_climo(io6,grids,i_counts,i_size,dtg(1:4),work_month,&
                         lat_range,lon_range)
        exit
      else if(istatus > 0) then
        write(errout,*) "ERROR: problem reading total_file_name"
        cycle file_loop
      else

        ! Extract month information from dtg
        ! ---------------------------------
        read(dtg(49:50), *) int_month

        ! If the month of the new file is greater than the current working
        ! month, call print_climo and print the grid values to the output
        ! file.
        if(work_month == -1) then
          work_month = int_month
        else if(work_month /= int_month) then
          call print_climo(io6,grids,i_counts,i_size,dtg(1:4),work_month,&
                           lat_range,lon_range)
          work_month = int_month
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
        call grid_raw_data_climo(errout,grids,i_counts,i_size,&
                lat_gridder,lat_thresh)

        ! Deallocate all the arrays for the next pass
        ! -------------------------------------------
        call clear_arrays
      endif
    enddo file_loop  
  endif

  close(io8)
  close(io6)
  close(errout)  
  deallocate(grids)
  deallocate(i_counts)
  deallocate(lat_range)
  deallocate(lon_range)
  

end program omi_JZ_climo
