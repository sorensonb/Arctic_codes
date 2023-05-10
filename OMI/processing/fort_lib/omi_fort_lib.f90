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
    ! This function determines if lat/lon point is within a box
    ! made by 4 lat/lon pairs.
    ! -------------------------------------------------------------
    function pixel_in_box(lats, lons, plat, plon) result(l_inside)

      real(kind = 8), dimension(4), intent(in)    :: lats 
      real(kind = 8), dimension(4), intent(in)    :: lons   
      real(kind = 8), intent(in)                  :: plat
      real(kind = 8), intent(in)                  :: plon

      real(kind = 8), dimension(4)                :: local_lons
      real(kind = 8)                              :: local_lon
      integer                                     :: njj
      logical                                     :: l_inside


      
      ! Adjust the lon and lon corners to account for pixels
      ! that straddle the antimeridian
      do njj = 1, 4
        if(lons(njj) < 0) then
          local_lons(njj) = lons(njj) + 360
        else
          local_lons(njj) = lons(njj)     
        endif
      enddo 

      if(plon < 0) then
        local_lon = plon + 360
      else
        local_lon = plon
      endif 

      ! Handle case if lat/lon box straddles the prime meridian
      if((maxval(local_lons) - minval(local_lons)) > 180) then
        do njj = 1, 4
          if(local_lons(njj) <= 180) then
            local_lons(njj) = local_lons(njj) + 360
          !else
          !  local_lons(njj) = local_lons(njj) + 360
          endif 
        enddo

        if(local_lon <= 180) then
          local_lon = local_lon + 360
        endif
      endif

      if( (plat <= maxval(lats) .and. plat >= minval(lats)) .and. &
          (local_lon >= minval(local_lons) .and. &
           local_lon <= maxval(local_lons))) then
        !write(*,*) local_lon
        l_inside = .true.
      else
        ! Perform the check again with the adjusted lons
        if( (plat <= maxval(lats) .and. plat >= minval(lats)) .and. &
            (local_lon >= minval(local_lons) .and. &
             local_lon <= maxval(local_lons))) then
          l_inside = .true.
        else
          l_inside = .false.
        endif
        !else
        !  l_inside = .false.
        !endif
    
      endif

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
