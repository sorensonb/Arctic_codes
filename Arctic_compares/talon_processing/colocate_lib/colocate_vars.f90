module colocate_vars
!
! NAME:
!   colocate_vars
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
  integer(hsize_t), dimension(:), allocatable       :: OMI_AI_raw_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_LATCRNR_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_LONCRNR_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_SZA_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_VZA_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_AZM_dims
  integer(hsize_t), dimension(:), allocatable       :: OMI_ALB_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_prep_AI_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_prep_SSA0_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_prep_SSA1_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_prep_SSA2_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_prep_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_prep_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_AI_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_SSA0_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_SSA1_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_SSA2_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_CH1_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_CH7_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_COD_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_CLD_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_CTP_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_SWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_CLD_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_ALB_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_LON_dims

  ! Dimensions for output data
  ! --------------------------
  integer(hsize_t), dimension(:), allocatable       :: TROP_out_AI_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_out_SSA0_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_out_SSA1_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_out_SSA2_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_out_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: TROP_out_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_CH1_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_CH7_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_COD_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_CLD_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_CTP_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: MODIS_out_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_out_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_out_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: NSIDC_out_LON_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_SWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LWF_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_CLD_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_ALB_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: CERES_out_LON_dims

  ! Data from input data
  ! -------------------- 
  real(kind=8), dimension(:,:), allocatable, target   :: OMI_AI_data
  real(kind=8), dimension(:,:), allocatable, target   :: OMI_AI_raw_data
  real(kind=8), dimension(:,:), allocatable, target   :: OMI_LAT_data
  real(kind=8), dimension(:,:), allocatable, target   :: OMI_LON_data
  real(kind=8), dimension(:,:,:), allocatable, target :: OMI_LATCRNR_data
  real(kind=8), dimension(:,:,:), allocatable, target :: OMI_LONCRNR_data
  real(kind=8), dimension(:,:), allocatable, target   :: OMI_SZA_data
  real(kind=8), dimension(:,:), allocatable, target   :: OMI_VZA_data
  real(kind=8), dimension(:,:), allocatable, target   :: OMI_AZM_data
  real(kind=8), dimension(:,:), allocatable, target   :: OMI_ALB_data
  real(kind=8), dimension(:),   allocatable, target   :: TROP_prep_AI_data
  real(kind=8), dimension(:),   allocatable, target   :: TROP_prep_SSA0_data
  real(kind=8), dimension(:),   allocatable, target   :: TROP_prep_SSA1_data
  real(kind=8), dimension(:),   allocatable, target   :: TROP_prep_SSA2_data
  real(kind=8), dimension(:),   allocatable, target   :: TROP_prep_LAT_data
  real(kind=8), dimension(:),   allocatable, target   :: TROP_prep_LON_data
  real(kind=8), dimension(:,:), allocatable, target   :: TROP_AI_data
  real(kind=8), dimension(:,:), allocatable, target   :: TROP_SSA0_data
  real(kind=8), dimension(:,:), allocatable, target   :: TROP_SSA1_data
  real(kind=8), dimension(:,:), allocatable, target   :: TROP_SSA2_data
  real(kind=8), dimension(:,:), allocatable, target   :: MODIS_CH1_data
  real(kind=8), dimension(:,:), allocatable, target   :: MODIS_CH7_data
  real(kind=8), dimension(:,:), allocatable, target   :: MODIS_COD_data
  real(kind=8), dimension(:,:), allocatable, target   :: MODIS_CLD_data
  real(kind=8), dimension(:,:), allocatable, target   :: MODIS_CTP_data
  real(kind=8), dimension(:,:), allocatable, target   :: MODIS_LAT_data
  real(kind=8), dimension(:,:), allocatable, target   :: MODIS_LON_data
  real(kind=8), dimension(:,:), allocatable, target   :: NSIDC_data
  real(kind=8), dimension(:,:), allocatable, target   :: NSIDC_LAT_data
  real(kind=8), dimension(:,:), allocatable, target   :: NSIDC_LON_data
  real(kind=8), dimension(:,:), allocatable, target   :: CERES_LWF_data
  real(kind=8), dimension(:,:), allocatable, target   :: CERES_SWF_data
  real(kind=8), dimension(:,:), allocatable, target   :: CERES_CLD_data
  real(kind=8), dimension(:,:), allocatable, target   :: CERES_ALB_data
  real(kind=8), dimension(:,:), allocatable, target   :: CERES_LAT_data
  real(kind=8), dimension(:,:), allocatable, target   :: CERES_LON_data

  ! Data for output data
  ! --------------------
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_CH1_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_CH7_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_COD_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_CLD_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_CTP_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: MODIS_out_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_out_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_out_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: NSIDC_out_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_SWF_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_CLD_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_ALB_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: CERES_out_LON_data
  real(kind=8), dimension(:,:), allocatable, target :: TROP_out_AI_data
  real(kind=8), dimension(:,:), allocatable, target :: TROP_out_SSA0_data
  real(kind=8), dimension(:,:), allocatable, target :: TROP_out_SSA1_data
  real(kind=8), dimension(:,:), allocatable, target :: TROP_out_SSA2_data
  real(kind=8), dimension(:,:), allocatable, target :: TROP_out_LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: TROP_out_LON_data


  contains
 
    subroutine allocate_out_arrays(dim1, dim2)
    ! --------------------------------------------------------------
    ! This subroutine allocates the output data and dimension arrays
    ! --------------------------------------------------------------

    integer(hsize_t)             :: dim1
    integer(hsize_t)             :: dim2
    integer             :: error


    ! Dimensions for output data
    ! --------------------------
    !!#!allocate(MODIS_out_CH1_dims)
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
    allocate(TROP_out_AI_data(dim1, dim2), stat = error)
    allocate(TROP_out_SSA0_data(dim1, dim2), stat = error)
    allocate(TROP_out_SSA1_data(dim1, dim2), stat = error)
    allocate(TROP_out_SSA2_data(dim1, dim2), stat = error)
    allocate(MODIS_out_CH1_data(dim1, dim2), stat = error)
    allocate(MODIS_out_CH7_data(dim1, dim2), stat = error)
    allocate(MODIS_out_COD_data(dim1, dim2), stat = error)
    allocate(MODIS_out_CLD_data(dim1, dim2), stat = error)
    allocate(MODIS_out_CTP_data(dim1, dim2), stat = error)
    allocate(NSIDC_out_data(dim1, dim2), stat = error)
    allocate(CERES_out_LWF_data(dim1, dim2), stat = error)
    allocate(CERES_out_SWF_data(dim1, dim2), stat = error)
    allocate(CERES_out_CLD_data(dim1, dim2), stat = error)
    allocate(CERES_out_ALB_data(dim1, dim2), stat = error)

    if ( error < 0 ) then
       write(*,*) " *** Error allocating H5datasets in output data"
       return
    endif
      
    write(*,*) "Allocated out data to dimensions",dim1,dim2

    end subroutine allocate_out_arrays
 
    subroutine clear_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################

      !write(*,*) "Inside clear_arrays"

      ! Dimensions for input data
      ! ------------------------- 
      if(allocated(OMI_AI_dims))           deallocate(OMI_AI_dims)
      if(allocated(OMI_AI_raw_dims))       deallocate(OMI_AI_raw_dims)
      if(allocated(OMI_LAT_dims))          deallocate(OMI_LAT_dims)
      if(allocated(OMI_LON_dims))          deallocate(OMI_LON_dims)
      if(allocated(OMI_LATCRNR_dims))      deallocate(OMI_LATCRNR_dims)
      if(allocated(OMI_LONCRNR_dims))      deallocate(OMI_LONCRNR_dims)
      if(allocated(OMI_SZA_dims))          deallocate(OMI_SZA_dims)
      if(allocated(OMI_VZA_dims))          deallocate(OMI_VZA_dims)
      if(allocated(OMI_AZM_dims))          deallocate(OMI_AZM_dims)
      if(allocated(OMI_ALB_dims))          deallocate(OMI_ALB_dims)
      if(allocated(TROP_prep_AI_dims))     deallocate(TROP_prep_AI_dims)
      if(allocated(TROP_prep_SSA0_dims))   deallocate(TROP_prep_SSA0_dims)
      if(allocated(TROP_prep_SSA1_dims))   deallocate(TROP_prep_SSA1_dims)
      if(allocated(TROP_prep_SSA2_dims))   deallocate(TROP_prep_SSA2_dims)
      if(allocated(TROP_prep_LAT_dims))    deallocate(TROP_prep_LAT_dims)
      if(allocated(TROP_prep_LON_dims))    deallocate(TROP_prep_LON_dims)
      if(allocated(TROP_AI_dims))          deallocate(TROP_AI_dims)
      if(allocated(TROP_SSA0_dims))        deallocate(TROP_SSA0_dims)
      if(allocated(TROP_SSA1_dims))        deallocate(TROP_SSA1_dims)
      if(allocated(TROP_SSA2_dims))        deallocate(TROP_SSA2_dims)
      if(allocated(MODIS_CH1_dims))        deallocate(MODIS_CH1_dims)
      if(allocated(MODIS_CH7_dims))        deallocate(MODIS_CH7_dims)
      if(allocated(MODIS_COD_dims))        deallocate(MODIS_COD_dims)
      if(allocated(MODIS_CLD_dims))        deallocate(MODIS_CLD_dims)
      if(allocated(MODIS_CTP_dims))        deallocate(MODIS_CTP_dims)
      if(allocated(MODIS_LAT_dims))        deallocate(MODIS_LAT_dims)
      if(allocated(MODIS_LON_dims))        deallocate(MODIS_LON_dims)
      if(allocated(NSIDC_dims))            deallocate(NSIDC_dims)
      if(allocated(NSIDC_LAT_dims))        deallocate(NSIDC_LAT_dims)
      if(allocated(NSIDC_LON_dims))        deallocate(NSIDC_LON_dims)
      if(allocated(CERES_LWF_dims))        deallocate(CERES_LWF_dims)
      if(allocated(CERES_SWF_dims))        deallocate(CERES_SWF_dims)
      if(allocated(CERES_CLD_dims))        deallocate(CERES_CLD_dims)
      if(allocated(CERES_ALB_dims))        deallocate(CERES_ALB_dims)
      if(allocated(CERES_LAT_dims))        deallocate(CERES_LAT_dims)
      if(allocated(CERES_LON_dims))        deallocate(CERES_LON_dims)

      ! Dimensions for output data
      ! --------------------------
      if(allocated(TROP_out_AI_dims))     deallocate(TROP_out_AI_dims)
      if(allocated(TROP_out_SSA0_dims))   deallocate(TROP_out_SSA0_dims)
      if(allocated(TROP_out_SSA1_dims))   deallocate(TROP_out_SSA1_dims)
      if(allocated(TROP_out_SSA2_dims))   deallocate(TROP_out_SSA2_dims)
      if(allocated(TROP_out_LAT_dims))    deallocate(TROP_out_LAT_dims)
      if(allocated(TROP_out_LON_dims))    deallocate(TROP_out_LON_dims)
      if(allocated(MODIS_out_CH1_dims))   deallocate(MODIS_out_CH1_dims)
      if(allocated(MODIS_out_CH7_dims))   deallocate(MODIS_out_CH7_dims)
      if(allocated(MODIS_out_COD_dims))   deallocate(MODIS_out_COD_dims)
      if(allocated(MODIS_out_CLD_dims))   deallocate(MODIS_out_CLD_dims)
      if(allocated(MODIS_out_CTP_dims))   deallocate(MODIS_out_CTP_dims)
      if(allocated(MODIS_out_LAT_dims))   deallocate(MODIS_out_LAT_dims)
      if(allocated(MODIS_out_LON_dims))   deallocate(MODIS_out_LON_dims)
      if(allocated(NSIDC_out_dims))       deallocate(NSIDC_out_dims)
      if(allocated(NSIDC_out_LAT_dims))   deallocate(NSIDC_out_LAT_dims)
      if(allocated(NSIDC_out_LON_dims))   deallocate(NSIDC_out_LON_dims)
      if(allocated(CERES_out_LWF_dims))   deallocate(CERES_out_LWF_dims)
      if(allocated(CERES_out_SWF_dims))   deallocate(CERES_out_SWF_dims)
      if(allocated(CERES_out_CLD_dims))   deallocate(CERES_out_CLD_dims)
      if(allocated(CERES_out_ALB_dims))   deallocate(CERES_out_ALB_dims)
      if(allocated(CERES_out_LAT_dims))   deallocate(CERES_out_LAT_dims)
      if(allocated(CERES_out_LON_dims))   deallocate(CERES_out_LON_dims)

      ! Data from input data
      ! -------------------- 
      if(allocated(OMI_AI_data))           deallocate(OMI_AI_data)
      if(allocated(OMI_AI_raw_data))       deallocate(OMI_AI_raw_data)
      if(allocated(OMI_LAT_data))          deallocate(OMI_LAT_data)
      if(allocated(OMI_LON_data))          deallocate(OMI_LON_data)
      if(allocated(OMI_LATCRNR_data))      deallocate(OMI_LATCRNR_data)
      if(allocated(OMI_LONCRNR_data))      deallocate(OMI_LONCRNR_data)
      if(allocated(OMI_SZA_data))          deallocate(OMI_SZA_data)
      if(allocated(OMI_VZA_data))          deallocate(OMI_VZA_data)
      if(allocated(OMI_AZM_data))          deallocate(OMI_AZM_data)
      if(allocated(OMI_ALB_data))          deallocate(OMI_ALB_data)
      if(allocated(TROP_prep_AI_data))     deallocate(TROP_prep_AI_data)
      if(allocated(TROP_prep_SSA0_data))   deallocate(TROP_prep_SSA0_data)
      if(allocated(TROP_prep_SSA1_data))   deallocate(TROP_prep_SSA1_data)
      if(allocated(TROP_prep_SSA2_data))   deallocate(TROP_prep_SSA2_data)
      if(allocated(TROP_prep_LAT_data))    deallocate(TROP_prep_LAT_data)
      if(allocated(TROP_prep_LON_data))    deallocate(TROP_prep_LON_data)
      if(allocated(TROP_AI_data))          deallocate(TROP_AI_data)
      if(allocated(TROP_SSA0_data))        deallocate(TROP_SSA0_data)
      if(allocated(TROP_SSA1_data))        deallocate(TROP_SSA1_data)
      if(allocated(TROP_SSA2_data))        deallocate(TROP_SSA2_data)
      if(allocated(MODIS_CH1_data))        deallocate(MODIS_CH1_data)
      if(allocated(MODIS_CH7_data))        deallocate(MODIS_CH7_data)
      if(allocated(MODIS_COD_data))        deallocate(MODIS_COD_data)
      if(allocated(MODIS_CLD_data))        deallocate(MODIS_CLD_data)
      if(allocated(MODIS_CTP_data))        deallocate(MODIS_CTP_data)
      if(allocated(MODIS_LAT_data))        deallocate(MODIS_LAT_data)
      if(allocated(MODIS_LON_data))        deallocate(MODIS_LON_data)
      if(allocated(NSIDC_data))            deallocate(NSIDC_data)
      if(allocated(NSIDC_LAT_data))        deallocate(NSIDC_LAT_data)
      if(allocated(NSIDC_LON_data))        deallocate(NSIDC_LON_data)
      if(allocated(CERES_LWF_data))        deallocate(CERES_LWF_data)
      if(allocated(CERES_SWF_data))        deallocate(CERES_SWF_data)
      if(allocated(CERES_CLD_data))        deallocate(CERES_CLD_data)
      if(allocated(CERES_ALB_data))        deallocate(CERES_ALB_data)
      if(allocated(CERES_LAT_data))        deallocate(CERES_LAT_data)
      if(allocated(CERES_LON_data))        deallocate(CERES_LON_data)

      ! Data for output data
      ! --------------------
      if(allocated(TROP_out_AI_data))     deallocate(TROP_out_AI_data)
      if(allocated(TROP_out_SSA0_data))   deallocate(TROP_out_SSA0_data)
      if(allocated(TROP_out_SSA1_data))   deallocate(TROP_out_SSA1_data)
      if(allocated(TROP_out_SSA2_data))   deallocate(TROP_out_SSA2_data)
      if(allocated(TROP_out_LAT_data))    deallocate(TROP_out_LAT_data)
      if(allocated(TROP_out_LON_data))    deallocate(TROP_out_LON_data)
      if(allocated(MODIS_out_CH1_data))   deallocate(MODIS_out_CH1_data)
      if(allocated(MODIS_out_CH7_data))   deallocate(MODIS_out_CH7_data)
      if(allocated(MODIS_out_COD_data))   deallocate(MODIS_out_COD_data)
      if(allocated(MODIS_out_CLD_data))   deallocate(MODIS_out_CLD_data)
      if(allocated(MODIS_out_CTP_data))   deallocate(MODIS_out_CTP_data)
      if(allocated(MODIS_out_LAT_data))   deallocate(MODIS_out_LAT_data)
      if(allocated(MODIS_out_LON_data))   deallocate(MODIS_out_LON_data)
      if(allocated(NSIDC_out_data))       deallocate(NSIDC_out_data)
      if(allocated(NSIDC_out_LAT_data))   deallocate(NSIDC_out_LAT_data)
      if(allocated(NSIDC_out_LON_data))   deallocate(NSIDC_out_LON_data)
      if(allocated(CERES_out_LWF_data))   deallocate(CERES_out_LWF_data)
      if(allocated(CERES_out_SWF_data))   deallocate(CERES_out_SWF_data)
      if(allocated(CERES_out_CLD_data))   deallocate(CERES_out_CLD_data)
      if(allocated(CERES_out_ALB_data))   deallocate(CERES_out_ALB_data)
      if(allocated(CERES_out_LAT_data))   deallocate(CERES_out_LAT_data)
      if(allocated(CERES_out_LON_data))   deallocate(CERES_out_LON_data)

    end subroutine clear_arrays

    ! -------------------------------------------------------------
    ! This function converts degrees to radians
    ! -------------------------------------------------------------
    function degrees_to_radians(degree) result(radians)

      real(kind = 8), intent(in)    :: degree
      real(kind = 8)                :: radians
      !real                :: radians

      radians = degree * 3.14159265 / 180.

    end function degrees_to_radians

    ! -------------------------------------------------------------
    ! This function converts radians to degrees
    ! -------------------------------------------------------------
    function radians_to_degrees(radian) result(degrees)

      real(kind = 8), intent(in)    :: radian
      real(kind = 8)                :: degrees
      !real                :: degrees

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

end module colocate_vars
