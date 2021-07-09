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

  integer                                           :: i_num_bad
  integer,dimension(:),allocatable                  :: i_bad_list
 
  integer(hsize_t), dimension(:), allocatable       :: AI_dims
  integer(hsize_t), dimension(:), allocatable       :: AZM_dims
  integer(hsize_t), dimension(:), allocatable       :: GPQF_dims
  integer(hsize_t), dimension(:), allocatable       :: LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: LON_dims
  integer(hsize_t), dimension(:), allocatable       :: XTRACK_dims
  real(kind=8), dimension(:,:), allocatable, target :: AI_data
  real(kind=8), dimension(:,:), allocatable, target :: AZM_data
  integer, dimension(:,:), allocatable, target      :: GPQF_data
  real(kind=8), dimension(:,:), allocatable, target :: LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: LON_data
  integer, dimension(:,:), allocatable, target      :: XTRACK_data

  contains
  
    subroutine clear_arrays
    ! -------------------------------------------------------------
    ! This subroutine deallocates all the data and dimension arrays
    ! -------------------------------------------------------------
    
      ! ############################################################################

      !write(*,*) "Inside clear_arrays"

      if(allocated(AI_dims))      deallocate(AI_dims)
      if(allocated(AZM_dims))     deallocate(AZM_dims)
      if(allocated(GPQF_dims))    deallocate(GPQF_dims)
      if(allocated(LAT_dims))     deallocate(LAT_dims)
      if(allocated(LON_dims))     deallocate(LON_dims)
      if(allocated(XTRACK_dims))  deallocate(XTRACK_dims)

      if(allocated(AI_data))      deallocate(AI_data)
      if(allocated(AZM_data))     deallocate(AZM_data)
      if(allocated(GPQF_data))    deallocate(GPQF_data)
      if(allocated(LAT_data))     deallocate(LAT_data)
      if(allocated(LON_data))     deallocate(LON_data)
      if(allocated(XTRACK_data))  deallocate(XTRACK_data)

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

end module h5_vars
