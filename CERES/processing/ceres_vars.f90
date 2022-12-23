module ceres_vars
!
! NAME:
!   ceres_vars
!
! PURPOSE:
!   Contain the data arrays and dimensions for all variable types, as well as
!   some functions.
! 
! CALLS:
!   None.
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/09:
!     Written
!
!  ############################################################################

  use hdf5

  implicit none

  integer                                           :: i_num_bad
  integer,dimension(:),allocatable                  :: i_bad_list

  ! Dimensions for input data
  ! ------------------------- 
  integer(hsize_t), dimension(:), allocatable       :: CERES_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_SWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_TIM_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_ALB_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_SZA_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_CLS_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_CLW_dims

  ! Dimensions for output data
  ! --------------------------
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_SWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_TIM_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_ALB_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_SZA_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_CLS_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_CLW_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_counts_dims
  

  ! Data from input data
  ! -------------------- 
  real(kind=8), dimension(:), allocatable, target   :: CERES_LAT_data
  real(kind=8), dimension(:), allocatable, target   :: CERES_LON_data
  real(kind=8), dimension(:), allocatable, target   :: CERES_SWF_data
  real(kind=8), dimension(:), allocatable, target   :: CERES_LWF_data
  real(kind=8), dimension(:), allocatable, target   :: CERES_TIM_data
  real(kind=8), dimension(:), allocatable, target   :: CERES_ALB_data
  real(kind=8), dimension(:), allocatable, target   :: CERES_SZA_data
  real(kind=8), dimension(:), allocatable, target   :: CERES_CLS_data
  real(kind=8), dimension(:), allocatable, target   :: CERES_CLW_data

  ! Data for output data
  ! --------------------
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_SWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_TIM_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_ALB_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_SZA_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_CLS_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_CLW_data
  integer(kind=8), dimension(:,:), allocatable, target :: CERES_out_swcnt_data
  integer(kind=8), dimension(:,:), allocatable, target :: CERES_out_lwcnt_data
  integer(kind=8), dimension(:,:), allocatable, target :: CERES_out_abcnt_data


  contains
 
    subroutine allocate_out_arrays(len1, len2)
      ! --------------------------------------------------------------
      ! This subroutine allocates the output data and dimension arrays
      ! --------------------------------------------------------------

      integer             :: len1, len2
      integer             :: error

      ! Dimensions for output data
      ! --------------------------
      allocate(CERES_out_SWF_data(len1, len2), stat = error)
      allocate(CERES_out_LWF_data(len1, len2), stat = error)
      allocate(CERES_out_swcnt_data(len1, len2), stat = error)
      allocate(CERES_out_lwcnt_data(len1, len2), stat = error)
      allocate(CERES_out_abcnt_data(len1, len2), stat = error)

      CERES_out_SWF_data(:,:) = 0.
      CERES_out_LWF_data(:,:) = 0.
      CERES_out_swcnt_data(:,:) = 0
      CERES_out_lwcnt_data(:,:) = 0

      if ( error < 0 ) then
         write(*,*) " *** Error allocating H5datasets in output data"
         return
      endif
        
      write(*,*) "Allocated out data to dimensions",len1, len2

    end subroutine allocate_out_arrays
 
    subroutine clear_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################

      !write(*,*) "Inside clear_arrays"

      ! Dimensions for input data
      ! ------------------------- 
      if(allocated(CERES_LAT_dims))   deallocate(CERES_LAT_dims)
      if(allocated(CERES_LON_dims))   deallocate(CERES_LON_dims)
      if(allocated(CERES_SWF_dims))   deallocate(CERES_SWF_dims)
      if(allocated(CERES_LWF_dims))   deallocate(CERES_LWF_dims)
      if(allocated(CERES_TIM_dims))   deallocate(CERES_TIM_dims)
      if(allocated(CERES_ALB_dims))   deallocate(CERES_ALB_dims)
      if(allocated(CERES_SZA_dims))   deallocate(CERES_SZA_dims)
      if(allocated(CERES_CLS_dims))   deallocate(CERES_CLS_dims)
      if(allocated(CERES_CLW_dims))   deallocate(CERES_CLW_dims)

      ! Dimensions for output data
      ! --------------------------
      if(allocated(CERES_out_LAT_dims))   deallocate(CERES_out_LAT_dims)
      if(allocated(CERES_out_LON_dims))   deallocate(CERES_out_LON_dims)
      if(allocated(CERES_out_SWF_dims))   deallocate(CERES_out_SWF_dims)
      if(allocated(CERES_out_LWF_dims))   deallocate(CERES_out_LWF_dims)
      if(allocated(CERES_out_TIM_dims))   deallocate(CERES_out_TIM_dims)
      if(allocated(CERES_out_ALB_dims))   deallocate(CERES_out_ALB_dims)
      if(allocated(CERES_out_SZA_dims))   deallocate(CERES_out_SZA_dims)
      if(allocated(CERES_out_CLS_dims))   deallocate(CERES_out_CLS_dims)
      if(allocated(CERES_out_CLW_dims))   deallocate(CERES_out_CLW_dims)

      ! Data from input data
      ! -------------------- 
      if(allocated(CERES_LAT_data))   deallocate(CERES_LAT_data)
      if(allocated(CERES_LON_data))   deallocate(CERES_LON_data)
      if(allocated(CERES_SWF_data))   deallocate(CERES_SWF_data)
      if(allocated(CERES_LWF_data))   deallocate(CERES_LWF_data)
      if(allocated(CERES_TIM_data))   deallocate(CERES_TIM_data)
      if(allocated(CERES_ALB_data))   deallocate(CERES_ALB_data)
      if(allocated(CERES_SZA_data))   deallocate(CERES_SZA_data)
      if(allocated(CERES_CLS_data))   deallocate(CERES_CLS_data)
      if(allocated(CERES_CLW_data))   deallocate(CERES_CLW_data)

      ! Data for output data
      ! --------------------
      if(allocated(CERES_out_LAT_data))   deallocate(CERES_out_LAT_data)
      if(allocated(CERES_out_LON_data))   deallocate(CERES_out_LON_data)
      if(allocated(CERES_out_SWF_data))   deallocate(CERES_out_SWF_data)
      if(allocated(CERES_out_LWF_data))   deallocate(CERES_out_LWF_data)
      if(allocated(CERES_out_TIM_data))   deallocate(CERES_out_TIM_data)
      if(allocated(CERES_out_ALB_data))   deallocate(CERES_out_ALB_data)
      if(allocated(CERES_out_SZA_data))   deallocate(CERES_out_SZA_data)
      if(allocated(CERES_out_CLS_data))   deallocate(CERES_out_CLS_data)
      if(allocated(CERES_out_CLW_data))   deallocate(CERES_out_CLW_data)
      if(allocated(CERES_out_swcnt_data))   deallocate(CERES_out_swcnt_data)
      if(allocated(CERES_out_lwcnt_data))   deallocate(CERES_out_lwcnt_data)
      if(allocated(CERES_out_abcnt_data))   deallocate(CERES_out_abcnt_data)

    end subroutine clear_arrays

    ! -------------------------------------------------------------
    ! This function converts degrees to radians
    ! -------------------------------------------------------------
    function degrees_to_radians(degree) result(radians)

      real(kind = 8), intent(in)    :: degree
      real                :: radians

      radians = degree * 3.14159265 / 180.

    end function degrees_to_radians

    ! -------------------------------------------------------------
    ! This function converts radians to degrees
    ! -------------------------------------------------------------
    function radians_to_degrees(radian) result(degrees)

      real(kind = 8), intent(in)    :: radian
      real                :: degrees

      degrees = radian * 180./ 3.14159265

    end function radians_to_degrees

    ! -------------------------------------------------------------
    ! This function finds the distance between 2 lat/lon pairs
    ! -------------------------------------------------------------
    function find_distance_between_points(lat1, lon1, lat2, lon2) &
            result(distance)

      real(kind = 8), intent(in)    :: lat1
      real(kind = 8), intent(in)    :: lon1
      real(kind = 8), intent(in)    :: lat2
      real(kind = 8), intent(in)    :: lon2

      real                :: r_lat1
      real                :: r_lon1
      real                :: r_lat2
      real                :: r_lon2

      real                :: r_e 
      real                :: distance

      real                :: dlat
      real                :: dlon
      real                :: const_a
      real                :: const_c

      r_e = 6371. ! km

      r_lat1 = degrees_to_radians(lat1)
      r_lon1 = degrees_to_radians(lon1)
      r_lat2 = degrees_to_radians(lat2)
      r_lon2 = degrees_to_radians(lon2)

      dlon = r_lon2 - r_lon1
      dlat = r_lat2 - r_lat1

      const_a = sin(dlat / 2)**2. + cos(r_lat1) * cos(r_lat2) * sin(dlon / 2)**2.
      const_c = 2. * atan2(sqrt(const_a), sqrt(1 - const_a))

      distance = r_e * const_c

    end function find_distance_between_points

    ! -------------------------------------------------------------
    ! This function determines if lat/lon point is within a box
    ! made by 4 lat/lon pairs.
    ! -------------------------------------------------------------
    function pixel_in_box(lats, lons, plat, plon) result(l_inside)

      real(kind = 8), dimension(4), intent(in)    :: lats 
      real(kind = 8), dimension(4), intent(in)    :: lons   
      real(kind = 8), intent(in)                  :: plat
      real(kind = 8), intent(in)                  :: plon

      logical                                     :: l_inside

      if( (plat <= maxval(lats) .and. plat >= minval(lats)) .and. &
          (plon >= minval(lons) .and. plon <= maxval(lons))) then
        l_inside = .true.
      else
        l_inside = .false.
      endif

    end function pixel_in_box
    
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

    end function bit_check

end module ceres_vars
