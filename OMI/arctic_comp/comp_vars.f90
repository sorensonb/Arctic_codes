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
 
  integer(hsize_t), dimension(:), allocatable       :: OMI_AI_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_CH2_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_SWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LON_dims
  !!#!integer(hsize_t), dimension(:), allocatable       :: AZM_dims
  !!#!integer(hsize_t), dimension(:), allocatable       :: GPQF_dims
  !!#!integer(hsize_t), dimension(:), allocatable       :: XTRACK_dims
  real(kind=8), dimension(:,:), allocatable, target :: OMI_AI_data
  real(kind=8), dimension(:,:), allocatable, target :: OMI_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: OMI_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_CH2_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_LWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_SWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_LON_data
  !!#!real(kind=8), dimension(:,:), allocatable, target :: AZM_data
  !!#!integer, dimension(:,:), allocatable, target      :: GPQF_data
  !!#!integer, dimension(:,:), allocatable, target      :: XTRACK_data

  contains
  
    subroutine clear_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################

      !write(*,*) "Inside clear_arrays"

      if(allocated(OMI_AI_dims))      deallocate(OMI_AI_dims)
      if(allocated(OMI_LAT_dims))     deallocate(OMI_LAT_dims)
      if(allocated(OMI_LON_dims))     deallocate(OMI_LON_dims)
      if(allocated(MODIS_CH2_dims))   deallocate(MODIS_CH2_dims)
      if(allocated(MODIS_LAT_dims))   deallocate(MODIS_LAT_dims)
      if(allocated(MODIS_LON_dims))   deallocate(MODIS_LON_dims)
      if(allocated(NSIDC_dims))       deallocate(NSIDC_dims)
      if(allocated(NSIDC_LAT_dims))   deallocate(NSIDC_LAT_dims)
      if(allocated(NSIDC_LON_dims))   deallocate(NSIDC_LON_dims)
      if(allocated(CERES_LWF_dims))   deallocate(CERES_LWF_dims)
      if(allocated(CERES_SWF_dims))   deallocate(CERES_SWF_dims)
      if(allocated(CERES_LAT_dims))   deallocate(CERES_LAT_dims)
      if(allocated(CERES_LON_dims))   deallocate(CERES_LON_dims)
      !!#!if(allocated(AZM_dims))     deallocate(AZM_dims)
      !!#!if(allocated(GPQF_dims))    deallocate(GPQF_dims)
      !!#!if(allocated(XTRACK_dims))  deallocate(XTRACK_dims)

      if(allocated(OMI_AI_data))      deallocate(OMI_AI_data)
      if(allocated(OMI_LAT_data))     deallocate(OMI_LAT_data)
      if(allocated(OMI_LON_data))     deallocate(OMI_LON_data)
      if(allocated(MODIS_CH2_data))   deallocate(MODIS_CH2_data)
      if(allocated(MODIS_LAT_data))   deallocate(MODIS_LAT_data)
      if(allocated(MODIS_LON_data))   deallocate(MODIS_LON_data)
      if(allocated(NSIDC_data))       deallocate(NSIDC_data)
      if(allocated(NSIDC_LAT_data))   deallocate(NSIDC_LAT_data)
      if(allocated(NSIDC_LON_data))   deallocate(NSIDC_LON_data)
      if(allocated(CERES_LWF_data))   deallocate(CERES_LWF_data)
      if(allocated(CERES_SWF_data))   deallocate(CERES_SWF_data)
      if(allocated(CERES_LAT_data))   deallocate(CERES_LAT_data)
      if(allocated(CERES_LON_data))   deallocate(CERES_LON_data)
      !!#!if(allocated(AZM_data))     deallocate(AZM_data)
      !!#!if(allocated(GPQF_data))    deallocate(GPQF_data)
      !!#!if(allocated(XTRACK_data))  deallocate(XTRACK_data)

    end subroutine clear_arrays

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
