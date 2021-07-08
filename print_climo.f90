subroutine print_climo(io6,grids,i_counts,i_size,c_year,work_month,&
                       lat_range,lon_range)
!
!  NAME:
!
!  PURPOSE:
!
!  CALLS:
!    None
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/06/10: Written
!
!  ###########################################################################

  implicit none

  integer                :: io6
  real                   :: grids(360,i_size)
  integer                :: i_counts(360,i_size)
  integer                :: i_size        ! array size
  character(len = 4)     :: c_year
  integer                :: work_month
  real,dimension(i_size) :: lat_range
  real,dimension(360)    :: lon_range

  integer                :: ii
  integer                :: jj

  character(len = 6)     :: outstring
  character(len = 255)   :: out_fmt

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
  !write(*,*) "In print_climo"
  if(work_month < 10) then
    outstring = c_year//'0'
    out_fmt = '(a5,i1,i4,i5,1x,f9.5,i6)'
  else
    outstring = c_year 
    out_fmt = '(a4,i2,i4,i5,1x,f9.5,i6)'
  endif    

  ! Loop over the grid and count up grids with high AI
  write(*,*) outstring,work_month
  do ii=1,i_size
    do jj=1,360
      write(io6,trim(out_fmt)) outstring,work_month,int(lat_range(ii)),&
        int(lon_range(jj)),grids(jj,ii),i_counts(jj,ii)
    enddo  
  enddo  


  ! Reset grid arrays
  grids(:,:) = 0.
  i_counts(:,:) = 0

end subroutine print_climo
