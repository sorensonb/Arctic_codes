module comp_vars
!
! NAME:
!   comp_vars
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
  integer(hsize_t), dimension(:), allocatable       :: OMI_AI_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_CH2_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_CH7_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_SWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LON_dims

  ! Dimensions for output data
  ! --------------------------
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_CH2_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_CH7_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_out_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_out_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_out_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_SWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LON_dims
  

  ! Data from input data
  ! -------------------- 
  real(kind=8), dimension(:,:), allocatable, target :: OMI_AI_data
  real(kind=8), dimension(:,:), allocatable, target :: OMI_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: OMI_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_CH2_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_CH7_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_LWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_SWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_LON_data

  ! Data for output data
  ! --------------------
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_CH2_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_CH7_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_out_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_out_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_out_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_SWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LON_data


  contains
 
    subroutine allocate_out_arrays
    ! --------------------------------------------------------------
    ! This subroutine allocates the output data and dimension arrays
    ! --------------------------------------------------------------

    integer             :: error


    ! Dimensions for output data
    ! --------------------------
    !!#!allocate(MODIS_out_CH2_dims)
    !!#!allocate(MODIS_out_CH7_dims)
    !!#!allocate(MODIS_out_LAT_dims)
    !!#!allocate(MODIS_out_LON_dims)
    !!#!allocate(NSIDC_out_dims)
    !!#!allocate(NSIDC_out_LAT_dims)
    !!#!allocate(NSIDC_out_LON_dims)
    !!#!allocate(CERES_out_LWF_dims)
    !!#!allocate(CERES_out_SWF_dims)
    !!#!allocate(CERES_out_LAT_dims)
    !!#!allocate(CERES_out_LON_dims)

    ! Dimensions for output data
    ! --------------------------
    allocate(MODIS_out_CH2_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(MODIS_out_CH7_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(MODIS_out_LAT_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(MODIS_out_LON_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(NSIDC_out_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(NSIDC_out_LAT_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(NSIDC_out_LON_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(CERES_out_LWF_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(CERES_out_SWF_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(CERES_out_LAT_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)
    allocate(CERES_out_LON_data(OMI_AI_dims(1), OMI_AI_dims(2)), stat = error)

    if ( error < 0 ) then
       write(*,*) " *** Error allocating H5datasets in output data"
       return
    endif
      
    write(*,*) "Allocated out data to dimensions",OMI_AI_dims

    end subroutine allocate_out_arrays
 
    subroutine clear_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################

      !write(*,*) "Inside clear_arrays"

      ! Dimensions for input data
      ! ------------------------- 
      if(allocated(OMI_AI_dims))      deallocate(OMI_AI_dims)
      if(allocated(OMI_LAT_dims))     deallocate(OMI_LAT_dims)
      if(allocated(OMI_LON_dims))     deallocate(OMI_LON_dims)
      if(allocated(MODIS_CH2_dims))   deallocate(MODIS_CH2_dims)
      if(allocated(MODIS_CH7_dims))   deallocate(MODIS_CH7_dims)
      if(allocated(MODIS_LAT_dims))   deallocate(MODIS_LAT_dims)
      if(allocated(MODIS_LON_dims))   deallocate(MODIS_LON_dims)
      if(allocated(NSIDC_dims))       deallocate(NSIDC_dims)
      if(allocated(NSIDC_LAT_dims))   deallocate(NSIDC_LAT_dims)
      if(allocated(NSIDC_LON_dims))   deallocate(NSIDC_LON_dims)
      if(allocated(CERES_LWF_dims))   deallocate(CERES_LWF_dims)
      if(allocated(CERES_SWF_dims))   deallocate(CERES_SWF_dims)
      if(allocated(CERES_LAT_dims))   deallocate(CERES_LAT_dims)
      if(allocated(CERES_LON_dims))   deallocate(CERES_LON_dims)

      ! Dimensions for output data
      ! --------------------------
      if(allocated(MODIS_out_CH2_dims))   deallocate(MODIS_out_CH2_dims)
      if(allocated(MODIS_out_CH7_dims))   deallocate(MODIS_out_CH7_dims)
      if(allocated(MODIS_out_LAT_dims))   deallocate(MODIS_out_LAT_dims)
      if(allocated(MODIS_out_LON_dims))   deallocate(MODIS_out_LON_dims)
      if(allocated(NSIDC_out_dims))       deallocate(NSIDC_out_dims)
      if(allocated(NSIDC_out_LAT_dims))   deallocate(NSIDC_out_LAT_dims)
      if(allocated(NSIDC_out_LON_dims))   deallocate(NSIDC_out_LON_dims)
      if(allocated(CERES_out_LWF_dims))   deallocate(CERES_out_LWF_dims)
      if(allocated(CERES_out_SWF_dims))   deallocate(CERES_out_SWF_dims)
      if(allocated(CERES_out_LAT_dims))   deallocate(CERES_out_LAT_dims)
      if(allocated(CERES_out_LON_dims))   deallocate(CERES_out_LON_dims)

      ! Data from input data
      ! -------------------- 
      if(allocated(OMI_AI_data))      deallocate(OMI_AI_data)
      if(allocated(OMI_LAT_data))     deallocate(OMI_LAT_data)
      if(allocated(OMI_LON_data))     deallocate(OMI_LON_data)
      if(allocated(MODIS_CH2_data))   deallocate(MODIS_CH2_data)
      if(allocated(MODIS_CH7_data))   deallocate(MODIS_CH7_data)
      if(allocated(MODIS_LAT_data))   deallocate(MODIS_LAT_data)
      if(allocated(MODIS_LON_data))   deallocate(MODIS_LON_data)
      if(allocated(NSIDC_data))       deallocate(NSIDC_data)
      if(allocated(NSIDC_LAT_data))   deallocate(NSIDC_LAT_data)
      if(allocated(NSIDC_LON_data))   deallocate(NSIDC_LON_data)
      if(allocated(CERES_LWF_data))   deallocate(CERES_LWF_data)
      if(allocated(CERES_SWF_data))   deallocate(CERES_SWF_data)
      if(allocated(CERES_LAT_data))   deallocate(CERES_LAT_data)
      if(allocated(CERES_LON_data))   deallocate(CERES_LON_data)

      ! Data for output data
      ! --------------------
      if(allocated(MODIS_out_CH2_data))   deallocate(MODIS_out_CH2_data)
      if(allocated(MODIS_out_CH7_data))   deallocate(MODIS_out_CH7_data)
      if(allocated(MODIS_out_LAT_data))   deallocate(MODIS_out_LAT_data)
      if(allocated(MODIS_out_LON_data))   deallocate(MODIS_out_LON_data)
      if(allocated(NSIDC_out_data))       deallocate(NSIDC_out_data)
      if(allocated(NSIDC_out_LAT_data))   deallocate(NSIDC_out_LAT_data)
      if(allocated(NSIDC_out_LON_data))   deallocate(NSIDC_out_LON_data)
      if(allocated(CERES_out_LWF_data))   deallocate(CERES_out_LWF_data)
      if(allocated(CERES_out_SWF_data))   deallocate(CERES_out_SWF_data)
      if(allocated(CERES_out_LAT_data))   deallocate(CERES_out_LAT_data)
      if(allocated(CERES_out_LON_data))   deallocate(CERES_out_LON_data)

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

    end function

    ! -------------------------------------------------------------
    ! This function extracts the ice flag from the whole GPQF
    ! -------------------------------------------------------------
    function get_ice_flags(i_gpqf) result(i_flag)

      integer,intent(in) :: i_gpqf
      integer            :: i_flag

      ! Check each bit and update i_flag accordingly
      i_flag = 64 * bit_check(i_gpqf,14) + &
               32 * bit_check(i_gpqf,13) + &
               16 * bit_check(i_gpqf,12) + &
               8  * bit_check(i_gpqf,11) + &
               4  * bit_check(i_gpqf,10) + &
               2  * bit_check(i_gpqf,9) + &
               1  * bit_check(i_gpqf,8) 

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

end module comp_vars
