program trop_colocate
!
! NAME:
!   trop_colocate.f90
!
! PURPOSE:
!   Calculate monthly averages of screened AI (following the JZ criteria) and
!   write the averages to an output file.
! 
! CALLS:
!   Modules:
!     - hdf5
!     - h5_vars (custom module, shared between count and climo analyses)
!   Subroutines:
!     - check_bad_rows
!     - clear_arrays (from h5_vars)
!     - grid_raw_data_climo
!     - print_climo
!     - read_h5_AI
!     - read_h5_LAT
!     - read_h5_LON
!     - read_h5_XTRACK
!     - read_h5_AZM
!     - read_h5_GPQF
!
! MODIFICATIONS:
!   Blake Sorenson <blake.sorenson@und.edu>     - 2021/07/09:
!     Written
!
!  ############################################################################

  use hdf5
  use colocate_vars, only : clear_arrays, i_bad_list, &
    find_distance_between_points, allocate_out_arrays, &
    pixel_in_box, &
    TROP_AI_data,  TROP_AI_dims, &
    TROP_LAT_data, TROP_LAT_dims, &
    TROP_LON_data, TROP_LON_dims, &
    OMI_AI_data,    OMI_AI_dims, &
    OMI_LAT_data,   OMI_LAT_dims, &
    OMI_LON_data,   OMI_LON_dims, &
    OMI_LATCRNR_data,   OMI_LATCRNR_dims, &
    OMI_LONCRNR_data,   OMI_LONCRNR_dims, &
    OMI_SZA_data,   OMI_SZA_dims, &
    OMI_VZA_data,   OMI_VZA_dims, &
    OMI_AZM_data,   OMI_AZM_dims, &
    TROP_out_AI_data

  implicit none

  integer                :: ii            ! loop counter
  integer                :: jj            ! loop counter
  integer                :: nii            ! loop counter
  integer                :: njj            ! loop counter
  !!#!integer                :: i_size        ! array size
  !!#!integer                :: int_month     ! integer variable for month
  !!#!integer                :: int_day       ! integer variable for day 
  integer                :: arg_count     ! Number of arguments passed to exec
  !!#!integer                :: work_month    ! currently analyzed month
  !!#!integer                :: work_day 
  integer                :: error         ! error flag
  integer                :: istatus
  integer                :: num_nan
  integer                :: omi_file_id         ! id for current HDF5 file
  integer                :: trop_file_id         ! id for current HDF5 file
  integer                :: out_file_id         ! id for current HDF5 file

  integer                :: dspace_id_OLT  ! OMI LAT
  integer                :: dspace_id_OLN  ! OMI LON
  integer                :: dspace_id_OSZ  ! OMI SZA
  integer                :: dspace_id_OVZ  ! OMI VZA
  integer                :: dspace_id_OAZ  ! OMI AZM
  integer                :: dspace_id_OAI  ! OMI AI
  integer                :: dspace_id_TAI  ! TROPOMI AI

  integer                :: dset_id_OLT  ! OMI LAT
  integer                :: dset_id_OLN  ! OMI LON
  integer                :: dset_id_OSZ  ! OMI SZA
  integer                :: dset_id_OVZ  ! OMI VZA
  integer                :: dset_id_OAZ  ! OMI AZM
  integer                :: dset_id_OAI  ! OMI AI
  integer                :: dset_id_TAI  ! TROPOMI AI

  integer                :: rank
  integer(hsize_t), dimension(2)    :: test_dims

  real                   :: distance
  real                   :: closest_dist
  real                   :: min_dist
  real                   :: match_lat
  real                   :: match_lon
  real                   :: min_lat
  real                   :: match_data2

  real                   :: run_trop_total_ai
  integer                :: count_trop_total 

  character(len = 255)      :: omi_name        ! filename
  character(len = 255)      :: trop_name       ! filename
  character(len = 255)      :: out_file_name   ! filename

  !!#!! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./trop_exec omi_name tropomi_name'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1, omi_name)
  call get_command_argument(2, trop_name)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"

  ! Open the TROPOMI AI HDF5 file
  ! ----------------------------
  call h5fopen_f(trim(trop_name), H5F_ACC_RDWR_F, trop_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open TROPOMI AI file'
    return
  endif

  ! Open the OMI file
  ! -----------------
  call h5fopen_f(trim(omi_name), H5F_ACC_RDWR_F, omi_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open OMI file'
    return
  endif

  call read_comp_OMI_AI(omi_file_id)
  call read_comp_OMI_LAT(omi_file_id)
  call read_comp_OMI_LON(omi_file_id)
  call read_comp_OMI_LATCRNR(omi_file_id)
  call read_comp_OMI_LONCRNR(omi_file_id)
  call read_comp_OMI_SZA(omi_file_id)
  call read_comp_OMI_VZA(omi_file_id)
  call read_comp_OMI_AZM(omi_file_id)

  call read_comp_TROP_AI(trop_file_id)
  call read_comp_TROP_LAT(trop_file_id)
  call read_comp_TROP_LON(trop_file_id)

  !test_dims = (/10, 20/)
  test_dims = (/OMI_AI_dims(1), OMI_AI_dims(2)/)
  !test_dims = OMI_AI_dims

  ! Allocate the output arrays
  ! --------------------------
  call allocate_out_arrays

  istatus = 0
  num_nan = 0 
  min_dist  = 50.
  closest_dist = 999999.
  match_lat  = -999.
  match_lon  = -999.
  min_lat = 65.
  match_data2 = -999.

  !!#!do ii = 1, 20
  !!#!  write(*,*) TROP_LAT_data(ii), TROP_LON_data(ii), TROP_AI_data(ii)
  !!#!enddo

  omi_loop1: do ii=1,OMI_LAT_dims(2)

    if(mod(ii, 50) == 0) then
      write(*,*) ii
    endif

    omi_loop2: do jj=1,OMI_LAT_dims(1) 

      ! Check if the current pixel is missing
      ! NEW: Check if the current pixel is above minlat
      ! -----------------------------------------------
      if(OMI_LAT_data(jj,ii) < min_lat) then
      !if(isnan(OMI_AI_data(jj,ii))) then
        TROP_out_AI_data(jj,ii) = -999.
        num_nan = num_nan + 1
      else
        
        count_trop_total = 0
        run_trop_total_ai = 0.0
        run_trop_total_ai = 0.0

        trop_loop1: do nii=1,TROP_AI_dims(1)
          !trop_loop2: do njj=1,MODIS_LAT_dims(1)

          ! Check if the current pixel is within the OMI pixel bounds
          ! ---------------------------------------------------------
          if(pixel_in_box(OMI_LATCRNR_data(:,jj,ii), OMI_LONCRNR_data(:,jj,ii), &
              TROP_LAT_data(nii), TROP_LON_data(nii))) then

              run_trop_total_ai = &
                  run_trop_total_ai + TROP_AI_data(nii)
              count_trop_total = count_trop_total + 1

          endif
          !enddo trop_loop2
        enddo trop_loop1

        if(count_trop_total == 0) then
          TROP_out_AI_data(jj,ii) = -999.
        else
          TROP_out_AI_data(jj,ii) = run_trop_total_ai / count_trop_total
        endif

        closest_dist = 999999.

      endif
    enddo omi_loop2
  enddo omi_loop1

  write(*,*) istatus, num_nan
  
  ! Open the output file
  ! --------------------
  write(*,*) trop_name(len(trim(trop_name)) - 16:len(trim(trop_name)) - 5)
  out_file_name = 'colocated_tropomi_'//trop_name(len(trim(trop_name)) - &
    16:len(trim(trop_name)) - 5)//'.hdf5'
  write(*,*) trim(out_file_name)

  call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open output file'
    return
  endif

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write OMI Latitude
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  rank = 2

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_lat', H5T_NATIVE_DOUBLE, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_lat'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_DOUBLE, OMI_LAT_data, OMI_AI_dims, &
                      error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_OLT, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_OLT, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close dataset'
    return
  endif

  write(*,*) 'Wrote OMI LAT'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write OMI Longitude
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_OLN, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_lon', H5T_NATIVE_DOUBLE, dspace_id_OLN,  &
                   dset_id_OLN, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_lon'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLN, H5T_NATIVE_DOUBLE, OMI_LON_data, OMI_AI_dims, &
                      error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_OLN, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_OLN, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif

  write(*,*) 'Wrote OMI LON'
  
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write OMI Solar Zenith Angle
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_OSZ, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'trop_sza', H5T_NATIVE_DOUBLE, dspace_id_OSZ,  &
                   dset_id_OSZ, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_sza'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OSZ, H5T_NATIVE_DOUBLE, OMI_SZA_data, OMI_AI_dims, &
                      error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_OSZ, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_OSZ, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif

  write(*,*) 'Wrote OMI OSZ'
  
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write OMI Viewing Zenith Angle
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_OVZ, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'trop_vza', H5T_NATIVE_DOUBLE, dspace_id_OVZ,  &
                   dset_id_OVZ, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_vza'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OVZ, H5T_NATIVE_DOUBLE, OMI_VZA_data, OMI_AI_dims, &
                      error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_OVZ, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_OVZ, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif

  write(*,*) 'Wrote OMI OVZ'
  
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write OMI Relative Azimuth Angle
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_OAZ, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_azm', H5T_NATIVE_DOUBLE, dspace_id_OAZ,  &
                   dset_id_OAZ, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_azm'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OAZ, H5T_NATIVE_DOUBLE, OMI_AZM_data, OMI_AI_dims, &
                      error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_OAZ, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_OAZ, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif

  write(*,*) 'Wrote OMI OAZ'
  
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write OMI AI data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_OAI, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_uvai_pert', H5T_NATIVE_DOUBLE, &
                   dspace_id_OAI,  dset_id_OAI, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_uvai_pert'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OAI, H5T_NATIVE_DOUBLE, OMI_AI_data, OMI_AI_dims, &
                      error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_OAI, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_OAI, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote OMI AI'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write TROPOMI AI data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_TAI, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'trop_ai', H5T_NATIVE_DOUBLE, &
                   dspace_id_TAI,  dset_id_TAI, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ai'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_TAI, H5T_NATIVE_DOUBLE, TROP_out_AI_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_TAI, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_TAI, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote MODIS CH2'

  ! Close output file
  ! -----------------
  call h5fclose_f(out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif

  write(*,*) 'Saved output file'//trim(out_file_name)
 
  ! Deallocate all the arrays for the next pass
  ! -------------------------------------------
  call clear_arrays

  call h5fclose_f(omi_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close OMI file'
    return
  endif

  call h5fclose_f(trop_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close TROPOMI file'
    return
  endif

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program trop_colocate
