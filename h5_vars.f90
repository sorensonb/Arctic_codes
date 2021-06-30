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
  integer(hsize_t), dimension(:), allocatable       :: LAT_dims
  integer(hsize_t), dimension(:), allocatable       :: LON_dims
  integer(hsize_t), dimension(:), allocatable       :: XTRACK_dims
  real(kind=8), dimension(:,:), allocatable, target :: AI_data
  real(kind=8), dimension(:,:), allocatable, target :: LAT_data
  real(kind=8), dimension(:,:), allocatable, target :: LON_data
  integer, dimension(:,:), allocatable, target      :: XTRACK_data

end module h5_vars
