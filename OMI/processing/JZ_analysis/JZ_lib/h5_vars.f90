module h5_vars
!
! NAME:
!   h5_vars
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

  integer                                             :: i_num_bad
  integer,dimension(:),allocatable                    :: i_bad_list
 
  integer(hsize_t), dimension(:), allocatable         :: AI_dims
  integer(hsize_t), dimension(:), allocatable         :: AZM_dims
  integer(hsize_t), dimension(:), allocatable         :: CLD_dims
  integer(hsize_t), dimension(:), allocatable         :: SZA_dims
  integer(hsize_t), dimension(:), allocatable         :: VZA_dims
  integer(hsize_t), dimension(:), allocatable         :: LAT_dims
  integer(hsize_t), dimension(:), allocatable         :: LON_dims
  integer(hsize_t), dimension(:), allocatable         :: LATCRNR_dims
  integer(hsize_t), dimension(:), allocatable         :: LONCRNR_dims
  integer(hsize_t), dimension(:), allocatable         :: GPQF_dims
  integer(hsize_t), dimension(:), allocatable         :: XTRACK_dims

  real(kind=8), dimension(:,:), allocatable, target   :: AI_data
  real(kind=8), dimension(:,:), allocatable, target   :: AZM_data
  real(kind=8), dimension(:,:), allocatable, target   :: CLD_data
  real(kind=8), dimension(:,:), allocatable, target   :: SZA_data
  real(kind=8), dimension(:,:), allocatable, target   :: VZA_data
  real(kind=8), dimension(:,:), allocatable, target   :: LAT_data
  real(kind=8), dimension(:,:), allocatable, target   :: LON_data
  real(kind=8), dimension(:,:,:), allocatable, target :: LATCRNR_data
  real(kind=8), dimension(:,:,:), allocatable, target :: LONCRNR_data
  integer, dimension(:,:), allocatable, target        :: GPQF_data
  integer, dimension(:,:), allocatable, target        :: XTRACK_data

  contains
  
    subroutine clear_h5_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################

      !write(*,*) "Inside clear_arrays"

      if(allocated(AI_dims))          deallocate(AI_dims)
      if(allocated(AZM_dims))         deallocate(AZM_dims)
      if(allocated(CLD_dims))         deallocate(CLD_dims)
      if(allocated(SZA_dims))         deallocate(SZA_dims)
      if(allocated(VZA_dims))         deallocate(VZA_dims)
      if(allocated(LAT_dims))         deallocate(LAT_dims)
      if(allocated(LON_dims))         deallocate(LON_dims)
      if(allocated(LATCRNR_dims))     deallocate(LATCRNR_dims)
      if(allocated(LONCRNR_dims))     deallocate(LONCRNR_dims)
      if(allocated(GPQF_dims))        deallocate(GPQF_dims)
      if(allocated(XTRACK_dims))      deallocate(XTRACK_dims)

      if(allocated(AI_data))          deallocate(AI_data)
      if(allocated(AZM_data))         deallocate(AZM_data)
      if(allocated(CLD_data))         deallocate(CLD_data)
      if(allocated(SZA_data))         deallocate(SZA_data)
      if(allocated(VZA_data))         deallocate(VZA_data)
      if(allocated(LATCRNR_data))     deallocate(LATCRNR_data)
      if(allocated(LONCRNR_data))     deallocate(LONCRNR_data)
      if(allocated(GPQF_data))        deallocate(GPQF_data)
      if(allocated(XTRACK_data))      deallocate(XTRACK_data)

    end subroutine clear_h5_arrays

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

end module h5_vars
