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

  !!#!integer                                             :: i_num_bad
  !!#!integer,dimension(:),allocatable                    :: i_bad_list
 
  !!#!integer(hsize_t), dimension(:), allocatable         :: AI_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: AZM_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: CLD_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: SZA_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: VZA_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: LAT_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: LON_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: LATCRNR_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: LONCRNR_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: GPQF_dims
  !!#!integer(hsize_t), dimension(:), allocatable         :: XTRACK_dims

  !!#!real(kind=8), dimension(:,:), allocatable, target   :: AI_data
  !!#!real(kind=8), dimension(:,:), allocatable, target   :: AZM_data
  !!#!real(kind=8), dimension(:,:), allocatable, target   :: CLD_data
  !!#!real(kind=8), dimension(:,:), allocatable, target   :: SZA_data
  !!#!real(kind=8), dimension(:,:), allocatable, target   :: VZA_data
  !!#!real(kind=8), dimension(:,:), allocatable, target   :: LAT_data
  !!#!real(kind=8), dimension(:,:), allocatable, target   :: LON_data
  !!#!real(kind=8), dimension(:,:,:), allocatable, target :: LATCRNR_data
  !!#!real(kind=8), dimension(:,:,:), allocatable, target :: LONCRNR_data
  !!#!integer, dimension(:,:), allocatable, target        :: GPQF_data
  !!#!integer, dimension(:,:), allocatable, target        :: XTRACK_data

  contains
  
    subroutine clear_daily_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################

      !write(*,*) "Inside clear_arrays"

      if(allocated(day_values))       deallocate(day_values)
      if(allocated(lat_values))       deallocate(lat_values)
      if(allocated(lon_values))       deallocate(lon_values)
      if(allocated(grid_AI))          deallocate(grid_AI)
      if(allocated(count_AI))         deallocate(count_AI)

    end subroutine clear_daily_arrays

end module daily_vars
