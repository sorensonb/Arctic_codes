subroutine count_ai_JZ(io6,grids,i_counts,i_size,ai_thresh,synop_idx,&
                        total_file_name,lat_range)
                        !c_work_dtg,lat_range)
!
!  NAME:
!    count_ai_JZ
!
!  PURPOSE:
!    Calculate the average AI for each quarter degree grid box
!    from the summed AI and counts. If the average AI perturbation exceeds
!    the desired threshold, increment a counter that is printed, along with
!    the dtg, at the end of the subroutine. Before leaving the subroutine,
!    the synoptic index is incremented, and the count variables are reset to 
!    0.
!
!  CALLS:
!    None
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/06/10: Written
!
!  ###########################################################################

  implicit none

  integer                :: io6                    ! output file object
  integer                :: synop_idx              ! synoptic time index
  integer                :: i_size                 ! lat array size
  real                   :: grids(1440,i_size)     ! gridded AI
  integer                :: i_counts(1440,i_size)  ! AI counts
  real                   :: ai_thresh              ! threshold AI value
  character(len = 255)   :: total_file_name        ! dtg from each line of file
  !character(len = 10)    :: c_work_dtg             ! dtg from each line of file
  real,dimension(i_size) :: lat_range              ! latitude grid

  integer        :: ai_count_60   ! good AI counter north of 60
  integer        :: ai_count_65   ! good AI counter north of 65
  integer        :: ai_count_70   ! good AI counter north of 70
  integer        :: ai_count_75   ! good AI counter north of 75
  integer        :: ai_count_80   ! good AI counter north of 80
  integer        :: ai_count_85   ! good AI counter north of 85
  integer        :: ii            ! loop counter 
  integer        :: jj            ! loop counter
  real           :: avg_ai

  integer,dimension(4)   :: synop_times 

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Re-define the synoptic times
  ! ---------------------------- 
  synop_times = [0,6,12,18] 

  ai_count_60 = 0 
  ai_count_65 = 0 
  ai_count_70 = 0 
  ai_count_75 = 0 
  ai_count_80 = 0 
  ai_count_85 = 0 

  ! Loop over the grid and count up grids with high average AI for this
  ! synoptic time.
  ! -------------------------------------------------------------------
  do ii=1,i_size
    do jj=1,1440
      if(i_counts(jj,ii) > 0) then
        avg_ai = grids(jj,ii)
        if(avg_ai > ai_thresh) then
          if(lat_range(ii) >= 60.) then
            ai_count_60 = ai_count_60 + 1
            if(lat_range(ii) >= 65.) then
              ai_count_65 = ai_count_65 + 1
              if(lat_range(ii) >= 70.) then
                ai_count_70 = ai_count_70 + 1
                if(lat_range(ii) >= 75.) then
                  ai_count_75 = ai_count_75 + 1
                  if(lat_range(ii) >= 80.) then
                    ai_count_80 = ai_count_80 + 1
                    if(lat_range(ii) >= 85.) then
                      ai_count_85 = ai_count_85 + 1
                    endif
                  endif
                endif
              endif
            endif
          endif
        endif 
      endif
    enddo  
  enddo  

  if(synop_times(synop_idx) < 12) then
    write(*,*) total_file_name(44:47)//total_file_name(49:52)//'0',&
        synop_times(synop_idx)
    write(io6,'(a9,i1,6(i6))') total_file_name(44:47)//total_file_name(49:52)//'0',synop_times(synop_idx),&
      ai_count_60,ai_count_65,ai_count_70,ai_count_75,ai_count_80,ai_count_85
  else
    write(*,*) total_file_name(44:47)//total_file_name(49:52),&
        synop_times(synop_idx)
    write(io6,'(a8,i2,6(i6))') total_file_name(44:47)//total_file_name(49:52),synop_times(synop_idx), &
      ai_count_60,ai_count_65,ai_count_70,ai_count_75,ai_count_80,ai_count_85
  endif    
  !!#!write(*,*) c_work_dtg
  !!#!write(io6,'(a10,2x,6(i6))') c_work_dtg,&
  !!#!  ai_count_60,ai_count_65,ai_count_70,ai_count_75,ai_count_80,ai_count_85

  ! Reset grid arrays
  ! -----------------
  synop_idx = synop_idx + 1
  if(synop_idx == 5) synop_idx = 1
  grids(:,:) = 0.
  i_counts(:,:) = 0

end subroutine count_ai_JZ
