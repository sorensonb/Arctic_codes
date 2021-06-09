program omi_frequency
!
! NAME:
!   omi_frequency.f90
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

  implicit none

  integer                :: ii            ! loop counter
  integer                :: jj            ! loop counter
  integer                :: index1        
  integer                :: index2        
  integer                :: synop_idx
  integer                :: i_size        ! array size
  integer                :: ai_count      ! good AI counter
  integer                :: int_hr        ! integer variable for hour

  real                   :: ai_thresh     ! threshold AI value
  real                   :: lat_thresh
  real                   :: lat_gridder
  real                   :: avg_ai

  logical                :: l_in_time

  real,dimension(:), allocatable :: lat_range
  real,dimension(:), allocatable :: lon_range   
  real,dimension(:,:), allocatable :: grids
  integer,dimension(:,:), allocatable :: i_counts

  integer,dimension(4)   :: synop_times 

  ! File read variables
  integer,parameter      :: io8    = 42
  integer,parameter      :: io7    = 22
  integer,parameter      :: errout = 9 
  integer                :: istatus

  character(len = 41)    :: data_path
  character(len = 12)    :: dtg

  ! Variables from each line in Shawn's files
  real                   :: lat
  real                   :: lon
  real                   :: raw_ai
  real                   :: filter
  real                   :: clean_ai
  real                   :: v5
  real                   :: v6
  real                   :: v7
  real                   :: v8
  real                   :: v9
  real                   :: v10
  real                   :: v11
  real                   :: v12
  real                   :: v13
  real                   :: v14

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  !call test_sub()
  !write(*,*) "hi"
  !return

  synop_times = [0,6,12,18] 
  synop_idx = 1

  data_path = "/Research/OMI/out_files-monthly.20210518/"

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
  open(errout, file = "omi_error.txt", iostat = istatus)
  if(istatus /= 0) then
    write(*,*) "error opening error file."
  endif

  ! Set up count variables to count the number of grid boxes with
  ! high AI values
  ! -------------------------------------------------------------
  ai_thresh = 0.3
  ai_count  = 0

  ! Read the file names from the file name file
  open(io8, file = "omi_dates.txt", iostat = istatus)
  if(istatus > 0) then
    write(errout,*) "ERROR: Problem reading 'omi_dates.txt'"
    return 
  else
    ! Loop over the file
    file_loop: do
      ! Read the current dtg from the file
      read(io8, *, iostat=istatus) dtg
      if(istatus < 0) then 
        write(*,*) "End of omi_dates.txt found."
        exit
      else if(istatus > 0) then
        write(errout,*) "ERROR: problem reading dtg"
        cycle file_loop
      else

        ! Extract time information from dtg
        ! ---------------------------------
        read(dtg(9:10), *) int_hr

        ! See if the hour exceeds the current 6 hr assimilation window
        ! ------------------------------------------------------------
        call synop_time_check(synop_idx, int_hr, l_in_time)
        if(.not. l_in_time) then 
          ! # # # # # # # # # # # 
          ! SUBROUTINE
          ! # # # # # # # # # # # 
          ! Loop over the grid and count up grids with high AI
          do ii=1,i_size
            do jj=1,1440
              if(i_counts(jj,ii) > 0) then
                avg_ai = grids(jj,ii)/i_counts(jj,ii)
                if(avg_ai > ai_thresh) then
                  ai_count = ai_count + 1
                endif 
              endif
            enddo  
          enddo  
        
          write(*,*) dtg(1:8),synop_times(synop_idx), ai_count
        
          ! Reset grid arrays
          synop_idx = synop_idx + 1
          if(synop_idx == 5) synop_idx = 1
          ai_count = 0     
          grids(:,:) = 0.
          i_counts(:,:) = 0
        endif  

        write(*,*) data_path//dtg

        ! Open the shawn file and look at contents
        open(io7, file = data_path//dtg, iostat = istatus)
        if(istatus /= 0) then
          write(errout, *) "ERROR: error opening file",data_path//dtg
          write(errout, *) "       cycling file_loop"
          cycle file_loop
        endif
        
        ! Loop over the file
        data_loop: do
          ! -------------------------
          read(io7, *, iostat = istatus)  &
            lat, lon, raw_ai, filter, clean_ai,v5,v6,v7,v8,v9,v10,&
              v11,v12,v13,v14
          if(istatus > 0) then
            write(errout, *) "ERROR: error reading data from ",data_path//dtg
            write(errout, *) "       cycling data_loop"
            cycle data_loop
          else if(istatus < 0) then
            write(errout, *) "End of data in file: ",data_path//dtg
            exit data_loop
          endif
          ! Read a line from the file

          if(lat > lat_thresh) then
            ! Average the data into the grid?
            ! -------------------------------
            index1 = floor(lat*4 - lat_gridder)
            index2 = floor(lon*4 + 720)

            if(index1 < 1) index1 = 1
            if(index1 > i_size) index1 = i_size
            if(index2 < 1) index2 = 1
            if(index2 > 1440) index2 = 1440

            grids(index2,index1) = ((grids(index2,index1) * &
                i_counts(index2,index1)) + clean_ai) / &
               (i_counts(index2,index1)+1)
            i_counts(index2,index1) = i_counts(index2,index1) + 1
          endif
          !! See if current value meets criteria
          !! -----------------------------------
          !if((clean_ai >= ai_thresh) .and. &
          !   (lat > lat_thresh)) then
          !  ai_count = ai_count + 1
          !  write(*,*) dtg, lat, lon
          !endif
           
        enddo data_loop
        close(io7)
      endif
    enddo file_loop  
  endif

  deallocate(grids)
  deallocate(i_counts)
  deallocate(lat_range)
  deallocate(lon_range)
  close(io8)
  close(errout)  
  

end program omi_frequency
