module daily_vars
!
! NAME:
!   daily_vars
!
! PURPOSE:
!   Contain the data arrays and dimensions for all variable types, as well as
!   some functions.
! 
! CALLS:
!   None.
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/09/11:
!     Written
!
!  ############################################################################

  use hdf5

  implicit none

  integer, dimension(:), allocatable, target            :: day_values
  real(kind = 8), dimension(:), allocatable, target     :: lat_values
  real(kind = 8), dimension(:), allocatable, target     :: lon_values
  real(kind = 8), dimension(:,:,:), allocatable, target :: grid_AI
  integer, dimension(:,:,:), allocatable, target        :: count_AI

  integer, dimension(:,:,:), allocatable                :: count_sfland
  integer, dimension(:,:,:), allocatable                :: count_seaice
  integer, dimension(:,:,:), allocatable                :: count_permice
  integer, dimension(:,:,:), allocatable                :: count_drysnow
  integer, dimension(:,:,:), allocatable                :: count_ocean  
  integer, dimension(:,:,:), allocatable                :: count_mixpixel
  integer, dimension(:,:,:), allocatable                :: count_other   

  integer, dimension(:,:,:), allocatable                :: grid_GPQF
    

  contains
 
    subroutine allocate_daily_arrays(dim1, dim2, dim3)

      ! ----------------------------------------------------------------------
      ! This subroutine allocates and initializes the grid arrays
      ! ----------------------------------------------------------------------

      integer               :: dim1  ! Longitude dimension
      integer               :: dim2  ! Latitude dimension
      integer               :: dim3  ! Day dimension
      integer               :: error

      ! ######################################################################

      error = 0

      allocate(day_values(dim3), stat = error)
      allocate(grid_AI(dim1, dim2, dim3), stat = error)
      allocate(count_AI(dim1, dim2, dim3), stat = error)
    
      allocate(count_sfland(dim1, dim2, dim3), stat = error)
      allocate(count_seaice(dim1, dim2, dim3), stat = error)
      allocate(count_permice(dim1, dim2, dim3), stat = error)
      allocate(count_drysnow(dim1, dim2, dim3), stat = error)
      allocate(count_ocean(dim1, dim2, dim3), stat = error)
      allocate(count_mixpixel(dim1, dim2, dim3), stat = error)
      allocate(count_other(dim1, dim2, dim3), stat = error)
 
      allocate(grid_GPQF(dim1, dim2, dim3), stat = error)
 
      ! Initialize grid arrays and set to 0 initially
      ! ----------------------------------------------
      day_values(:)   = 0
      grid_AI(:,:,:)  = 0.
      count_AI(:,:,:) = 0 

      count_sfland(:,:,:)   = 0
      count_seaice(:,:,:)   = 0
      count_permice(:,:,:)  = 0
      count_drysnow(:,:,:)  = 0
      count_ocean(:,:,:)    = 0
      count_mixpixel(:,:,:) = 0
      count_other(:,:,:)    = 0

      grid_GPQF(:,:,:)      = 0

      if ( error < 0) then
        write(*,*) '*** Error allocating H5datasets in output data'
        return
      endif

      write(*,*) 'Allocated out data to dimensions', dim1, dim2, dim3

    end subroutine allocate_daily_arrays
 
    subroutine clear_daily_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ######################################################################

      !write(*,*) "Inside clear_arrays"

      if(allocated(day_values))       deallocate(day_values)
      if(allocated(lat_values))       deallocate(lat_values)
      if(allocated(lon_values))       deallocate(lon_values)
      if(allocated(grid_AI))          deallocate(grid_AI)
      if(allocated(count_AI))         deallocate(count_AI)

      if(allocated(count_sfland))     deallocate(count_sfland)
      if(allocated(count_seaice))     deallocate(count_seaice)
      if(allocated(count_permice))    deallocate(count_permice)
      if(allocated(count_drysnow))    deallocate(count_drysnow)
      if(allocated(count_ocean))      deallocate(count_ocean)
      if(allocated(count_mixpixel))   deallocate(count_mixpixel)
      if(allocated(count_other))      deallocate(count_other)

      if(allocated(grid_GPQF))        deallocate(grid_GPQF)

    end subroutine clear_daily_arrays

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


end module daily_vars
