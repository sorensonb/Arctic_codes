module type_vars
!
! NAME:
!   type_vars
!
! PURPOSE:
!    Holds allocated arrays for generating a gridded, multidimensional
!    climatology of the arctic_comp colocated subset files. 
!
! CALLS:
!   None.
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/05/05:
!     Written
!
!  ############################################################################

  use hdf5

  implicit none

  integer, dimension(:), allocatable  :: day_values
  real, dimension(:), allocatable     :: lat_range
  real, dimension(:), allocatable     :: lon_range
  real, dimension(:,:), allocatable   :: grid_areas

  ! If HDF writing doesn't work, try adding (kind = 8) to the grid declaration
  real, dimension(:,:,:), allocatable      :: grid_OMI_AI
  integer, dimension(:,:,:), allocatable   :: count_OMI_AI
  real, dimension(:,:,:), allocatable      :: grid_TROP_AI
  integer, dimension(:,:,:), allocatable   :: count_TROP_AI
  real, dimension(:,:,:), allocatable      :: grid_MODIS_COD
  integer, dimension(:,:,:), allocatable   :: count_MODIS_COD
  real, dimension(:,:,:), allocatable      :: grid_MODIS_CH7
  integer, dimension(:,:,:), allocatable   :: count_MODIS_CH7
  real, dimension(:,:,:), allocatable      :: grid_MODIS_CLD
  integer, dimension(:,:,:), allocatable   :: count_MODIS_CLD_0
  integer, dimension(:,:,:), allocatable   :: count_MODIS_CLD_1
  integer, dimension(:,:,:), allocatable   :: count_MODIS_CLD_2
  integer, dimension(:,:,:), allocatable   :: count_MODIS_CLD_3
  real, dimension(:,:,:), allocatable      :: grid_NSIDC_ICE
  integer, dimension(:,:,:), allocatable   :: count_NSIDC_ICE_0100
  integer, dimension(:,:,:), allocatable   :: count_NSIDC_ICE_251
  integer, dimension(:,:,:), allocatable   :: count_NSIDC_ICE_252
  integer, dimension(:,:,:), allocatable   :: count_NSIDC_ICE_253
  integer, dimension(:,:,:), allocatable   :: count_NSIDC_ICE_254


  contains
 
    subroutine allocate_grid_arrays(dim1, dim2, dim3)

    ! --------------------------------------------------------------
    ! This subroutine allocates the output data and dimension arrays
    ! --------------------------------------------------------------

    !integer(hsize_t)             :: dim1
    !integer(hsize_t)             :: dim2
    integer                      :: dim1
    integer                      :: dim2
    integer                      :: dim3
    integer                      :: error
 
    error = 0

    ! Dimensions for output data
    ! --------------------------
    allocate(grid_OMI_AI(dim1, dim2, dim3), stat = error)
    allocate(grid_TROP_AI(dim1, dim2, dim3), stat = error)
    allocate(grid_NSIDC_ICE(dim1, dim2, dim3), stat = error)
    allocate(grid_MODIS_COD(dim1, dim2, dim3), stat = error)
    allocate(grid_MODIS_CLD(dim1, dim2, dim3), stat = error)
    allocate(grid_MODIS_CH7(dim1, dim2, dim3), stat = error)

    allocate(count_OMI_AI(dim1, dim2, dim3), stat = error)
    allocate(count_TROP_AI(dim1, dim2, dim3), stat = error)
    allocate(count_MODIS_COD(dim1, dim2, dim3), stat = error)
    allocate(count_MODIS_CH7(dim1, dim2, dim3), stat = error)
    allocate(count_MODIS_CLD_0(dim1, dim2, dim3), stat = error)
    allocate(count_MODIS_CLD_1(dim1, dim2, dim3), stat = error)
    allocate(count_MODIS_CLD_2(dim1, dim2, dim3), stat = error)
    allocate(count_MODIS_CLD_3(dim1, dim2, dim3), stat = error)
    allocate(count_NSIDC_ICE_0100(dim1, dim2, dim3), stat = error)
    allocate(count_NSIDC_ICE_251(dim1, dim2, dim3), stat = error)
    allocate(count_NSIDC_ICE_252(dim1, dim2, dim3), stat = error)
    allocate(count_NSIDC_ICE_253(dim1, dim2, dim3), stat = error)
    allocate(count_NSIDC_ICE_254(dim1, dim2, dim3), stat = error)

    grid_OMI_AI(:,:,:)     = -999.
    grid_TROP_AI(:,:,:)    = -999.
    grid_NSIDC_ICE(:,:,:)  = -999.
    grid_MODIS_COD(:,:,:)  = -999.
    grid_MODIS_CLD(:,:,:)  = -999.
    grid_MODIS_CH7(:,:,:)  = -999.

    count_OMI_AI(:,:,:)     = -9 
    count_TROP_AI(:,:,:)    = -9 
    count_MODIS_COD(:,:,:)  = -9
    count_MODIS_CH7(:,:,:)  = -9

    count_MODIS_CLD_0(:,:,:)  = 0
    count_MODIS_CLD_1(:,:,:)  = 0
    count_MODIS_CLD_2(:,:,:)  = 0
    count_MODIS_CLD_3(:,:,:)  = 0

    count_NSIDC_ICE_0100(:,:,:) = 0
    count_NSIDC_ICE_251(:,:,:)  = 0
    count_NSIDC_ICE_252(:,:,:)  = 0
    count_NSIDC_ICE_253(:,:,:)  = 0
    count_NSIDC_ICE_254(:,:,:)  = 0

    if ( error < 0 ) then
       write(*,*) " *** Error allocating H5datasets in output data"
       return
    endif
      
    write(*,*) "Allocated out data to dimensions",dim1,dim2,dim3

    end subroutine allocate_grid_arrays
 
    subroutine clear_grid_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################

      !write(*,*) "Inside clear_arrays"

      ! Dimensions for input data
      ! ------------------------- 
      if(allocated(grid_OMI_AI))              deallocate(grid_OMI_AI)
      if(allocated(grid_TROP_AI))             deallocate(grid_TROP_AI)
      if(allocated(grid_NSIDC_ICE))           deallocate(grid_NSIDC_ICE)
      if(allocated(grid_MODIS_COD))           deallocate(grid_MODIS_COD)
      if(allocated(grid_MODIS_CLD))           deallocate(grid_MODIS_CLD)
      if(allocated(grid_MODIS_CH7))           deallocate(grid_MODIS_CH7)

      if(allocated(count_OMI_AI))             deallocate(count_OMI_AI)
      if(allocated(count_TROP_AI))            deallocate(count_TROP_AI)
      if(allocated(count_MODIS_COD))          deallocate(count_MODIS_COD)
      if(allocated(count_MODIS_CH7))          deallocate(count_MODIS_CH7)
      if(allocated(count_MODIS_CLD_0))        deallocate(count_MODIS_CLD_0)
      if(allocated(count_MODIS_CLD_1))        deallocate(count_MODIS_CLD_1)
      if(allocated(count_MODIS_CLD_2))        deallocate(count_MODIS_CLD_2)
      if(allocated(count_MODIS_CLD_3))        deallocate(count_MODIS_CLD_3)

      if(allocated(count_NSIDC_ICE_0100))     deallocate(count_NSIDC_ICE_0100)
      if(allocated(count_NSIDC_ICE_251))      deallocate(count_NSIDC_ICE_251)
      if(allocated(count_NSIDC_ICE_252))      deallocate(count_NSIDC_ICE_252)
      if(allocated(count_NSIDC_ICE_253))      deallocate(count_NSIDC_ICE_253)
      if(allocated(count_NSIDC_ICE_254))      deallocate(count_NSIDC_ICE_254)

    end subroutine clear_grid_arrays

end module type_vars
