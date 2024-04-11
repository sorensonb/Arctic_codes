subroutine write_force_output_file(out_file_name)
!
! NAME:
!   write_output_file.f90
!
! PURPOSE:
! 
! CALLS:
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2023/08/08:
!     Written
!
!  ############################################################################

  use hdf5

  implicit none

  character(len = 255)   :: out_file_name

  integer(hsize_t)                :: out_file_id

  ! File read variables
  integer,parameter      :: io8    = 42   ! File object for file name file
  integer                :: error

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  write(*,*) trim(out_file_name)

  ! Open the output file
  ! --------------------
  call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open output file'
    write(*,*) '             Not writing to file'
    return
  endif

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Call subroutines to write variables to output file
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  call write_force_efficiency(out_file_id)
  call write_force_intercepts(out_file_id)
  call write_force_errors(out_file_id)
  !!#!call write_force_OMI_LAT(out_file_id) 
  !!#!call write_force_OMI_LON(out_file_id) 
  !!#!call write_force_OMI_SZA(out_file_id) 
  !!#!call write_force_OMI_VZA(out_file_id) 
  !!#!call write_force_OMI_AZM(out_file_id) 
  !!#!call write_force_OMI_AI(out_file_id) 
  !!#!call write_force_OMI_AI_RAW(out_file_id) 
  !!#!call write_force_TROP_AI(out_file_id) 
  !!#!call write_force_TROP_SSA0(out_file_id) 
  !!#!call write_force_TROP_SSA1(out_file_id) 
  !!#!call write_force_TROP_SSA2(out_file_id) 

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! End subroutine calls
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Close the output file
  ! ---------------------
  call h5fclose_f(out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif
  
end subroutine write_force_output_file
