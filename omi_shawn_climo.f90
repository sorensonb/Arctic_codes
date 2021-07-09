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

  implicit none

  integer                :: ii            ! loop counter
  integer                :: i_size        ! array size
  integer                :: int_month     ! integer variable for month
  integer                :: arg_count     ! Number of arguments passed to exec
  integer                :: work_month    ! currently analyzed month
  integer                :: istatus       ! error flag

  real                   :: lat_thresh    ! latitude threshold. Only analyze
                                          ! data north of this value.
  real                   :: lat_gridder   ! Used for setting up the latitude
                                          ! grid

  real,dimension(:), allocatable      :: lat_range ! latitude grid
  real,dimension(:), allocatable      :: lon_range ! longitude grid
  real,dimension(:,:), allocatable    :: grids     ! quarter-degree AI grid
                                                   ! values.
  integer,dimension(:,:), allocatable :: i_counts  ! quarter-degree AI counts

  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer,parameter      :: io7    = 22   ! File object for each data file 
  integer,parameter      :: io6    = 1827 ! Data output file
  integer,parameter      :: errout = 9    ! File object for error file

  character(len = 255)   :: data_path       ! path to data files
  character(len = 12)    :: dtg             ! dtg from each line of file
  character(len = 255)   :: out_file_name   ! output file name
  character(len = 255)   :: date_file_name  ! name of file containing the 
                                            ! list of shawn file names to 
                                            ! be analyzed

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./omi_exec out_file_name date_file_name'
    return
  endif

  call get_command_argument(1,out_file_name)
  call get_command_argument(2,date_file_name)

  ! Initialize the working month to -1
  ! ----------------------------------
  work_month = -1

  data_path = "/Research/OMI/out_files-monthly.20210518/"

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
  open(errout, file = "omi_error_climo.txt", iostat = istatus)
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
        ! Print the final month of data to the output file using
        ! print_climo
        ! --------------------------------------------------------
        call print_climo(io6,grids,i_counts,i_size,dtg(1:4),work_month,&
                         lat_range,lon_range)
        exit
      else if(istatus > 0) then
        write(errout,*) "ERROR: problem reading dtg: "//dtg
        cycle file_loop
      else

        ! Extract month information from dtg
        ! ---------------------------------
        read(dtg(5:6), *) int_month

        ! If the month of the new file is greater than the current working
        ! month, call print_climo and print the grid values to the output
        ! file.
        ! ----------------------------------------------------------------
        if(work_month == -1) then
          work_month = int_month
        else if(work_month /= int_month) then
          call print_climo(io6,grids,i_counts,i_size,dtg(1:4),work_month,&
                           lat_range,lon_range)
          work_month = int_month
        endif

        ! Open the shawn file and look at contents
        open(io7, file = trim(data_path)//dtg, iostat = istatus)
        if(istatus /= 0) then
          write(errout, *) "ERROR: error opening file",trim(data_path)//dtg
          write(errout, *) "       cycling file_loop"
          cycle file_loop
        endif

        ! Read data from the current file and insert into the grids
        ! ---------------------------------------------------------
        call read_shawn_file_climo(io7,errout,trim(data_path)//dtg,grids,&
                            i_counts,i_size,lat_gridder,lat_thresh)
        
        ! Close the current shawn file
        ! ----------------------------
        close(io7)
      endif
    enddo file_loop  
  endif

  ! Deallocate the remaining allocated arrays and close all files
  ! -------------------------------------------------------------
  close(io8)
  close(io6)
  close(errout)  
  deallocate(grids)
  deallocate(i_counts)
  deallocate(lat_range)
  deallocate(lon_range)
  

end program omi_shawn_climo
