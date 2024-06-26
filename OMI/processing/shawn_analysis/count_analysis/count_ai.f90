subroutine count_ai(io6,grids,i_counts,i_size,ai_thresh,synop_idx,&
                        dtg,lat_range,grid_areas)
!
!  NAME:
!    count_ai
!
!  PURPOSE:
!    Calculate the average AI perturbation for each quarter degree grid box
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

  !use omi_fort_lib, only : lat_lon_area

  implicit none

  integer                :: io6                    ! output file object
  integer                :: synop_idx              ! synoptic time index
  integer                :: i_size                 ! lat array size
  real                   :: grids(1440,i_size)     ! gridded AI
  integer                :: i_counts(1440,i_size)  ! AI counts
  real                   :: ai_thresh              ! threshold AI value
  character(len = 12)    :: dtg                    ! datetimegroup
  real                   :: lat_range(i_size)      ! latitude grid
  real                   :: grid_areas(i_size)     ! lat/lon box areas

  real                   :: ai_count_60   ! good AI counter north of 65
  real                   :: ai_count_65   ! good AI counter north of 65
  real                   :: ai_count_70   ! good AI counter north of 70
  real                   :: ai_count_75   ! good AI counter north of 75
  real                   :: ai_count_80   ! good AI counter north of 80
  real                   :: ai_count_85   ! good AI counter north of 85
  integer                :: ii            ! loop counter 
  integer                :: jj            ! loop counter
  real                   :: avg_ai

  integer,dimension(4)   :: synop_times 

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  ! Re-define the synoptic times
  ! ---------------------------- 
  synop_times = [0,6,12,18] 

  ai_count_60 = 0. 
  ai_count_65 = 0. 
  ai_count_70 = 0. 
  ai_count_75 = 0. 
  ai_count_80 = 0. 
  ai_count_85 = 0. 

  ! Loop over the grid and count up grids with high average AI for this
  ! synoptic time.
  ! -------------------------------------------------------------------
  do ii=1,i_size
    do jj=1,1440
      if(i_counts(jj,ii) > 0) then
        avg_ai = grids(jj,ii)
        if(avg_ai > ai_thresh) then
          if(lat_range(ii) >= 60.) then
            ai_count_60 = ai_count_60 + grid_areas(ii)
            !ai_count_60 = ai_count_60 + 1
            if(lat_range(ii) >= 65.) then
              ai_count_65 = ai_count_65 + grid_areas(ii)
              !ai_count_65 = ai_count_65 + 1
              if(lat_range(ii) >= 70.) then
                ai_count_70 = ai_count_70 + grid_areas(ii)
                !ai_count_70 = ai_count_70 + 1
                if(lat_range(ii) >= 75.) then
                  ai_count_75 = ai_count_75 + grid_areas(ii)
                  !ai_count_75 = ai_count_75 + 1
                  if(lat_range(ii) >= 80.) then
                    ai_count_80 = ai_count_80 + grid_areas(ii)
                    !ai_count_80 = ai_count_80 + 1
                    if(lat_range(ii) >= 85.) then
                      ai_count_85 = ai_count_85 + grid_areas(ii)
                      !ai_count_85 = ai_count_85 + 1
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

  ! For VSJ22, printing daily averages
  ! ----------------------------------
  !!#!if(synop_times(synop_idx) < 12) then
  !!#!  write(*,*) dtg(1:8)//'0',synop_times(synop_idx)
  !!#!  write(io6,'(a9,i1,6(i6))') dtg(1:8)//'0',synop_times(synop_idx),&
  !!#!    ai_count_60,ai_count_65,ai_count_70,ai_count_75,ai_count_80,ai_count_85
  !!#!else
  !!#!  write(*,*) dtg(1:8),synop_times(synop_idx)
  !!#!  write(io6,'(a8,i2,6(i6))') dtg(1:8),synop_times(synop_idx), &
  !!#!    ai_count_60,ai_count_65,ai_count_70,ai_count_75,ai_count_80,ai_count_85
  !!#!endif    
  write(*,*) dtg(1:8)
  write(io6,'(a8,6(1x,f10.0))') dtg(1:8),&
    ai_count_60,ai_count_65,ai_count_70,ai_count_75,ai_count_80,ai_count_85

  ! Reset grid arrays
  ! -----------------
  synop_idx = synop_idx + 1
  if(synop_idx == 5) synop_idx = 1
  grids(:,:) = 0.
  i_counts(:,:) = 0

end subroutine count_ai
