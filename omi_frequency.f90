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
  integer,parameter      :: io8    = 42
  integer,parameter      :: io7    = 22
  integer,parameter      :: errout = 9 
  integer                :: istatus

  character(len = 41)    :: data_path
  character(len = 12)    :: dtg

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

  data_path = "/Research/OMI/out_files-monthly.20210518/"
  dtg = "200804221345"

  ! Open debug file
  open(errout, file = "omi_error.txt", iostat = istatus)
  if(istatus /= 0) then
    write(*,*) "Error opening error file."
  endif

  ! Set up count variables to count the number of grid boxes with
  ! high AI values

  ! Loop over file names, once hour of file name goes outside the
  ! +/- 3 hrs from synoptic time, print counts and date

  ! Read the file names from the file name file
  open(io8, file = "omi_dates.txt", iostat = istatus)
  if(istatus /= 0) then
    write(errout,*) "ERROR: Problem reading 'omi_dates.txt'"
  else
    ! Loop over the file
    file_loop: do
      ! Read the current dtg from the file
      read(io8, *, iostat=istatus) dtg
      if(istatus /= 0) then 
        write(errout,*) "End of omi_dates.txt found."
        exit
      else
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
          read(io7, *, iostat = istatus)  &
            lat, lon, raw_ai, filter, clean_ai,v5,v6,v7,v8,v9,v10,&
              v11,v12,v13,v14
          if(istatus /= 0) then
            write(errout, *) "ERROR: error reading data from ",data_path//dtg
            write(errout, *) "       cycling data_loop"
            cycle data_loop
          endif

          if(lat > 70.0) then
            write(*,*) data_path//dtg, lat, lon
          endif
           
        enddo data_loop
      endif
    enddo file_loop  
  endif
   
  
  

end program omi_frequency
