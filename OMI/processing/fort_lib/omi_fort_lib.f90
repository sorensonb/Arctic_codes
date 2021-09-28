module omi_fort_lib
!
! NAME:
!
! PURPOSE:
! 
! CALLS:
!   None.
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/08/18:
!     Written
!
!  ############################################################################

  implicit none

  !real,dimension(:,:)

  contains

    ! -------------------------------------------------------------
    ! This function returns the surface area of a lat x lon box in
    ! square kilometers.
    ! -------------------------------------------------------------
    function lat_lon_area(lat1,lat0,lon1,lon0) result(area_out)
      real, intent(in)   :: lat1
      real, intent(in)   :: lat0
      real, intent(in)   :: lon1
      real, intent(in)   :: lon0

      real               :: area_out

      real               :: pi
   
      pi = 3.1415926535 

      area_out = (pi/180.) * (6371.**2.) * abs(sin(lat1 * pi/180.) - &
                 sin(lat0 * pi/180.)) * abs(lon1 - lon0)

    end function
  
    ! -------------------------------------------------------------
    ! This function returns 1 if the bit is True and 0 if the bit
    ! is False
    ! -------------------------------------------------------------
    function bit_check(i_gpqf,i_index) result(i_out)

      integer,intent(in)  :: i_gpqf
      integer,intent(in)  :: i_index

      integer :: i_out

      i_out = 0

      if(btest(i_gpqf,i_index)) i_out = 1

    end function

end module omi_fort_lib
