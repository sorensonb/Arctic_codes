subroutine count_ai(grids,i_counts,i_size,ai_thresh,synop_idx,&
                        ai_count,dtg)
!
!  NAME:
!
!  PURPOSE:
!
!  CALLS:
!
!  MODIFICATIONS:
!
!  ###########################################################################

  implicit none

  integer        :: synop_idx
  integer        :: i_size        ! array size
  integer        :: ai_count      ! good AI counter
  real           :: grids(1440,i_size)
  integer        :: i_counts(1440,i_size)
  real           :: ai_thresh     ! threshold AI value


  integer        :: ii
  integer        :: jj
  real           :: avg_ai

  character(len = 12)    :: dtg
  integer,dimension(4)   :: synop_times 
 
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
        
  write(*,*) dtg(1:8),synop_times(synop_idx), ai_count
        
  ! Reset grid arrays
  synop_idx = synop_idx + 1
  if(synop_idx == 5) synop_idx = 1
  ai_count = 0     
  grids(:,:) = 0.
  i_counts(:,:) = 0

end subroutine count_ai
