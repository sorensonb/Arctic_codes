module h5_vars
!
! NAME:
!   omi_frequency_JZ.f90
!
! PURPOSE:
! 
! callS:
!   mie_calc.f90
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2018/10/24:
!     Written
!
!  ############################################################################

  use hdf5

  implicit none

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

    function integer2binary(i) result(b)
    ! -------------------------------------------------------------
    ! This function converts an integer to binary
    ! -------------------------------------------------------------
      integer,intent(in) :: i
      integer :: b(32)
      integer k,j
      b=0
      j=i
      do k=1,size(b)
        b(k)=mod(j,2)
        j=j/2
      enddo
    end function

end module h5_vars
