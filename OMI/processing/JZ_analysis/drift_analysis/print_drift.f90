subroutine print_drift(io6,avg_ai,i_count,c_year,work_month)
!
!  NAME:
!    print_drift
!
!  PURPOSE:
!    Print the current monthly averages to the output file.
!
!  CALLS:
!    None
!
!  MODIFICATIONS:
!    Blake Sorenson <blake.sorenson@und.edu>     - 2021/06/10: Written
!
!  ###########################################################################

  implicit none

  integer                :: io6                  ! output file object
  real                   :: avg_ai               ! avg AI value
  integer                :: i_count              ! AI counts
  character(len = 4)     :: c_year               ! string year
  integer                :: work_month           ! current month

  character(len = 6)     :: outstring            ! working string for date
  character(len = 255)   :: out_fmt              ! format string for output

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Set up the output string differently if the month only has 1 digit
  ! ------------------------------------------------------------------ 
  if(work_month < 10) then
    outstring = c_year//'0'
    out_fmt = '(a5,i1,i4,i5,1x,f9.5,i9)'
  else
    outstring = c_year 
    out_fmt = '(a4,i2,i4,i5,1x,f9.5,i9)'
  endif    

  ! Loop over the grid and write the climatology
  ! --------------------------------------------
  write(*,*) outstring,work_month
  write(io6,trim(out_fmt)) outstring,work_month,-20,-160,avg_ai,i_count

  ! Reset avg values 
  ! ----------------
  avg_ai  = 0.
  i_count = 0

end subroutine print_drift
