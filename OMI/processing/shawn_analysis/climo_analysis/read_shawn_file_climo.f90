subroutine read_shawn_file_climo(io7,errout,c_total_file_name,grids,i_counts,i_size,&
             lat_gridder,lat_thresh)
!
!  NAME:
!    read_shawn_file_climo
!
!  PURPOSE:
!    Loop over each line from the current OMI AI perturbation text file and
!    insert each AI perturbation value into the grid if it meets the criteria.
!
!  CALLS:
!    None
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/06/10: Written
!
!  ###########################################################################

  implicit none

  integer                :: io7                   ! file pointer for shawn file
  integer                :: errout                ! error file pointer
  integer                :: i_size                ! array size
  character(len = 255)   :: c_total_file_name     ! total file name
  real                   :: grids(360,i_size)    ! quarter-degree AI grid
  integer                :: i_counts(360,i_size) ! AI counts
  real                   :: lat_gridder           !
  real                   :: lat_thresh            ! threshold lat value

  integer                :: istatus               ! read status
  integer                :: index1                ! omi file index
  integer                :: index2                ! omi file index

  ! Variables from each line in Shawn's files
  real                   :: lat                     ! Latitude
  real                   :: lon                     ! Longitude
  real                   :: raw_ai                  ! Original Value                            
  real                   :: filter                  ! Altered Value   ! WRONG
  real                   :: clean_ai                ! Climo Value
  real                   :: v5                      ! Time (of day)
  real                   :: v6                      ! Szea
  real                   :: v7                      ! Vzea
  real                   :: v8                      ! Raza
  real                   :: v9                      ! Albedo 1 (for bin)
  real                   :: v10                     ! Albedo 2 (not for bin)
  real                   :: v11                     ! Normalized Radiance 1
  real                   :: v12                     ! Normalized Radiance 2
  real                   :: cld_frc                 ! Cloud (flag?)
  real                   :: v14                     ! Ground flag
  real                   :: v15  ! added for SJ5    ! Row number (base 1)                                             
                                                    
  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Loop over the file
  ! -------------------------
  data_loop: do
    read(io7, *, iostat = istatus)  &
            lat, lon, raw_ai, filter, clean_ai,v5,v6,v7,v8,v9,v10,&
              v11,v12,cld_frc,v14,v15
    if(istatus > 0) then
      write(errout, *) "ERROR: error reading data from ", &
          trim(c_total_file_name)
      write(errout, *) "       cycling data_loop"
      cycle data_loop
    else if(istatus < 0) then
      exit data_loop
    endif
    ! Read a line from the file

    if(lat > lat_thresh) then
      ! SJ5: only use values from rows 55 - 60
      ! --------------------------------------
      !if(v15 >= 55.) then
      ! SJ42: add max cloud fraction of 0.2
      !if(cld_frc < 0.2) then
      ! SJ43: add max cloud fraction of 0.5
      if(cld_frc < 0.5) then
        ! Average the data into the grid?
        ! -------------------------------
        index1 = floor(lat - lat_gridder) + 1
        index2 = floor(lon + 180) + 1

        if(index1 < 1) index1 = 1
        if(index1 > i_size) index1 = i_size
        if(index2 < 1) index2 = 1
        if(index2 > 360) index2 = 360

        ! Average the current value into the grid
        ! ---------------------------------------
        grids(index2,index1) = ((grids(index2,index1) * &
            i_counts(index2,index1)) + clean_ai) / &
           (i_counts(index2,index1)+1)
        i_counts(index2,index1) = i_counts(index2,index1) + 1
      endif
      !endif
    endif
  enddo data_loop

end subroutine read_shawn_file_climo
