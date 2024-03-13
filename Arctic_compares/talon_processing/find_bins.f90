subroutine find_bins(min_sza, max_sza, min_ice, max_ice, len_sza, len_ice)
!
! NAME:
!   find_bins.f90
!
! PURPOSE:
!   Based on the specified grid resolution and the minimum latitude,
!     allocates and fills the lat/lon grid arrays and fills the
!     grid box area array.
! 
! CALLS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/08/08:
!     Written
!
!  ############################################################################

  use omi_fort_lib, only : lat_lon_area
  use type_vars, only: lat_range, lon_range, grid_areas

  implicit none

  real                   :: resolution
  real                   :: lat_thresh
  real                   :: lat_gridder
  real                   :: lon_gridder

  integer                :: i_lat_size
  integer                :: i_lon_size

  integer                :: ii

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  lat_gridder = lat_thresh * (1 / resolution)
  lon_gridder = -180 * (1 / resolution)
  i_lat_size = (90. - lat_thresh) * (1 / resolution)
  i_lon_size = (360 * (1 / resolution))

  allocate(lat_range(i_lat_size))
  allocate(lon_range(i_lon_size))
  allocate(grid_areas(i_lat_size, i_lon_size))

  do ii=1,i_lon_size
    lon_range(ii) = -180.0 + (ii - 1) * resolution + (resolution / 2)
  enddo

  do ii=1,i_lat_size
    lat_range(ii) = lat_thresh + (ii - 1) * resolution + (resolution / 2)
    
    ! Since the distance between longitude lines is the same at any latitude
    ! band, only need a 1-d array to hold all the areas along each latitude
    ! line.
    ! ----------------------------------------------------------------------
    grid_areas(ii,:) = lat_lon_area(lat_range(ii) + resolution,&
      lat_range(ii),90.0 + resolution,90.0)
  enddo
 
end subroutine find_bins
