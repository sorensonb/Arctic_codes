subroutine count_ai(io6,grids,i_counts,i_size,ai_thresh,synop_idx,&
                        ai_count,dtg)
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

  implicit none

  integer        :: io6
  integer        :: synop_idx
  integer        :: i_size        ! array size
  integer        :: ai_count      ! good AI counter
  real           :: grids(1440,i_size)
  integer        :: i_counts(1440,i_size)
  real           :: ai_thresh     ! threshold AI value


  integer        :: ii
  integer        :: jj
  real           :: avg_ai

  character(len = 255)   :: out_string
  character(len = 12)    :: dtg
  integer,dimension(4)   :: synop_times 

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  synop_times = [0,6,12,18] 

  !write(*,*) "In count_ai"

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

  if(synop_times(synop_idx) < 12) then
    write(io6,'(a9,i1,i6)') dtg(1:8)//'0',synop_times(synop_idx), ai_count
  else
    write(io6,'(a8,i2,i6)') dtg(1:8),synop_times(synop_idx), ai_count
  endif    

  ! Reset grid arrays
  synop_idx = synop_idx + 1
  if(synop_idx == 5) synop_idx = 1
  ai_count = 0     
  grids(:,:) = 0.
  i_counts(:,:) = 0

end subroutine count_ai
