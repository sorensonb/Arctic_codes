subroutine write_comp_output_file(out_file_name, l_trop_found)
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
  logical                :: l_trop_found

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

  call write_coloc_OMI_LAT(out_file_id) 
  call write_coloc_OMI_LON(out_file_id) 
  call write_coloc_OMI_SZA(out_file_id) 
  call write_coloc_OMI_VZA(out_file_id) 
  call write_coloc_OMI_AZM(out_file_id) 
  call write_coloc_OMI_ALB(out_file_id) 
  call write_coloc_OMI_AI(out_file_id) 
  call write_coloc_OMI_AI_raw(out_file_id) 
  call write_coloc_CERES_LWF(out_file_id) 
  call write_coloc_CERES_SWF(out_file_id) 
  call write_coloc_MODIS_CH1(out_file_id) 
  call write_coloc_MODIS_CH7(out_file_id) 
  call write_coloc_MODIS_COD(out_file_id) 
  call write_coloc_MODIS_CTP(out_file_id) 
  call write_coloc_MODIS_CLD(out_file_id) 
  call write_coloc_NSIDC_ICE(out_file_id) 
  if(l_trop_found) then
    call write_coloc_TROP_UVAI_comp(out_file_id) 
    !call write_coloc_TROP_SSA0(out_file_id) 
    !call write_coloc_TROP_SSA1(out_file_id) 
    !call write_coloc_TROP_SSA2(out_file_id) 
  endif

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
  
end subroutine write_comp_output_file
