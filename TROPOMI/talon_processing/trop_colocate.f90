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

  use omp_lib
  use hdf5
  use colocate_vars, only : clear_arrays, i_bad_list, &
    find_distance_between_points, allocate_out_arrays, &
    pixel_in_box, &
    TROP_prep_AI_data,    TROP_prep_AI_dims, &
    !!#!TROP_prep_SSA0_data,  TROP_prep_SSA0_dims, &
    !!#!TROP_prep_SSA1_data,  TROP_prep_SSA1_dims, &
    !!#!TROP_prep_SSA2_data,  TROP_prep_SSA2_dims, &
    TROP_prep_LAT_data, TROP_prep_LAT_dims, &
    TROP_prep_LON_data, TROP_prep_LON_dims, &
    TROP_out_AI_data!!#!, TROP_out_SSA0_data, TROP_out_SSA1_data, &
      !!#!TROP_out_SSA2_data
  use h5_vars, only: &
    clear_h5_arrays, &
    AI_data,         AI_dims, &
    LAT_data,        LAT_dims, &
    LON_data,        LON_dims, &
    LATCRNR_data,    LATCRNR_dims, &
    LONCRNR_data,    LONCRNR_dims, &
    SZA_data,        SZA_dims, &
    VZA_data,        VZA_dims, &
    AZM_data,        AZM_dims
    

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
  integer(hsize_t)                :: omi_file_id         ! id for current HDF5 file
  integer(hsize_t)                :: trop_file_id         ! id for current HDF5 file
  integer(hsize_t)                :: out_file_id         ! id for current HDF5 file

  integer(hsize_t)                :: dspace_id_OLT  ! OMI LAT
  integer(hsize_t)                :: dspace_id_OLN  ! OMI LON
  integer(hsize_t)                :: dspace_id_OSZ  ! OMI SZA
  integer(hsize_t)                :: dspace_id_OVZ  ! OMI VZA
  integer(hsize_t)                :: dspace_id_OAZ  ! OMI AZM
  integer(hsize_t)                :: dspace_id_OAI  ! OMI AI
  integer(hsize_t)                :: dspace_id_TAI  ! TROPOMI AI
  !!#!integer(hsize_t)                :: dspace_id_TS0  ! TROPOMI SSA0
  !!#!integer(hsize_t)                :: dspace_id_TS1  ! TROPOMI SSA1
  !!#!integer(hsize_t)                :: dspace_id_TS2  ! TROPOMI SSA2

  integer(hsize_t)                :: dset_id_OLT  ! OMI LAT
  integer(hsize_t)                :: dset_id_OLN  ! OMI LON
  integer(hsize_t)                :: dset_id_OSZ  ! OMI SZA
  integer(hsize_t)                :: dset_id_OVZ  ! OMI VZA
  integer(hsize_t)                :: dset_id_OAZ  ! OMI AZM
  integer(hsize_t)                :: dset_id_OAI  ! OMI AI
  integer(hsize_t)                :: dset_id_TAI  ! TROPOMI AI
  !!#!integer(hsize_t)                :: dset_id_TS0  ! TROPOMI SSA0
  !!#!integer(hsize_t)                :: dset_id_TS1  ! TROPOMI SSA1
  !!#!integer(hsize_t)                :: dset_id_TS2  ! TROPOMI SSA2

  integer                :: rank
  integer(hsize_t), dimension(2)    :: test_dims

  integer       :: thread_id
  integer       :: swaths_per_thread
  integer       :: begin_idx
  integer       :: end_idx
  integer       :: num_threads


  real                   :: min_lat

  real                   :: run_trop_total_ai
  integer                :: count_trop_ai
 
  real(kind = 8), dimension(4) :: local_lons
  real(kind = 8)               :: local_lon

  character(len = 12)       :: omi_date        ! filename
  character(len = 255)      :: omi_file_name   ! filename
  character(len = 255)      :: omi_just_name   ! filename
  character(len = 255)      :: omi_name_file   ! filename
  character(len = 255)      :: trop_file_name  ! filename
  character(len = 255)      :: out_file_name   ! filename

  integer                    :: io7  ! omi file name file

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./trop_exec trop_prep_file omi_name_file'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1, trop_file_name)
  call get_command_argument(2, omi_name_file)

  ! Open the OMI name file
  io7 = 42
  open(io7, file = trim(omi_name_file), iostat = istatus)
  if(istatus > 0) then
    write(*,*) "ERROR: Problem reading "//trim(omi_name_file)
    return 
  endif
 
  ! Read the OMI file name  
  read(io7, *, iostat = istatus) omi_just_name
 
  close(io7)

  ! NOTE: ON TALON, CHANGE THIS TO THE CORRECT PATH. WILL BE WORKING
  !       IN THE ONE DIRECTORY
  ! ----------------------------------------------------------------
  omi_file_name = '/Research/OMI/H5_files/'//trim(omi_just_name)

  write(*,*) trim(trop_file_name), '  ', trim(omi_file_name)


  ! Set up the out file name based on the OMI timestamp
  ! ---------------------------------------------------
  omi_date = omi_file_name(len(trim(omi_file_name)) - 46 : &
                           len(trim(omi_file_name)) - 43)//&
             omi_file_name(len(trim(omi_file_name)) - 41 : &
                           len(trim(omi_file_name)) - 38)//&
             omi_file_name(len(trim(omi_file_name)) - 36 : &
                           len(trim(omi_file_name)) - 33)
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
  call h5fopen_f(trim(trop_file_name), H5F_ACC_RDWR_F, trop_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open TROPOMI AI file at:'
    write(*,*) trim(trop_file_name)
    return
  endif

  ! Open the OMI file
  ! -----------------
  call h5fopen_f(trim(omi_file_name), H5F_ACC_RDWR_F, omi_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open OMI file at '
    write(*,*) trim(omi_file_name) 
    return
  endif



  call read_h5_AI(omi_file_id)
  call read_h5_LAT(omi_file_id)
  call read_h5_LON(omi_file_id)
  call read_h5_LATCRNR(omi_file_id)
  call read_h5_LONCRNR(omi_file_id)
  call read_h5_SZA(omi_file_id)
  call read_h5_VZA(omi_file_id)
  call read_h5_AZM(omi_file_id)

  call read_prep_TROP_AI(trop_file_id)
  !!#!call read_prep_TROP_SSA0(trop_file_id)
  !!#!call read_prep_TROP_SSA1(trop_file_id)
  !!#!call read_prep_TROP_SSA2(trop_file_id)
  call read_prep_TROP_LAT(trop_file_id)
  call read_prep_TROP_LON(trop_file_id)

  !test_dims = (/10, 20/)
  test_dims = (/AI_dims(1), AI_dims(2)/)
  !test_dims = AI_dims

  ! Allocate the output arrays
  ! --------------------------
  call allocate_out_arrays(AI_dims(1), AI_dims(2))
  write(*,*) "HERE2", LAT_data(20,20)

  istatus = 0
  num_nan = 0 
  min_lat = 65.

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! OPENMP CODE
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!$omp parallel private(thread_id, begin_idx, end_idx, run_trop_total_ai, count_trop_ai, local_lon)
  !!#!  thread_id   = omp_get_thread_num()
  !!#!  num_threads = omp_get_num_threads()
  !!#!  swaths_per_thread = (OMI_AI_dims(2) + num_threads - 1) / num_threads
  !!#!  
  !!#!  begin_idx = (swaths_per_thread * thread_id) + 1
  !!#!  end_idx   = swaths_per_thread * (thread_id + 1)
  !!#!  if(end_idx > OMI_AI_dims(2)) then
  !!#!    end_idx = end_idx - 1
  !!#!  endif
  !!#!  
  !!#!  omi_loop1: do ii = begin_idx, end_idx 
  !!#!  !!#!omi_loop1: do ii=1,LAT_dims(2)

  !!#!    !!#!if(mod(ii, 50) == 0) then
  !!#!    !!#!  write(*,*) ii
  !!#!    !!#!endif

  !!#!    omi_loop2: do jj=1,LAT_dims(1) 

  !!#!      ! Check if the current pixel is missing
  !!#!      ! NEW: Check if the current pixel is above minlat
  !!#!      ! -----------------------------------------------
  !!#!      if(LAT_data(jj,ii) < min_lat) then
  !!#!      !if(isnan(AI_data(jj,ii))) then
  !!#!        TROP_out_AI_data(jj,ii) = -999.
  !!#!        num_nan = num_nan + 1
  !!#!      else
  !!#!        
  !!#!        count_trop_ai     = 0
  !!#!        run_trop_total_ai = 0.0

  !!#!        ! Adjust the lon and lon corners to account for pixels
  !!#!        ! that straddle the antimeridian
  !!#!        do njj = 1, 4
  !!#!          if(LONCRNR_data(njj,jj,ii) < 0) then
  !!#!            local_lons(njj) = LONCRNR_data(njj,jj,ii) + 360
  !!#!          else
  !!#!            local_lons(njj) = LONCRNR_data(njj,jj,ii)     
  !!#!          endif
  !!#!        enddo
 
  !!#!        if((maxval(local_lons) - minval(local_lons)) > 180) then
  !!#!          do njj = 1, 4
  !!#!            if(local_lons(njj) <= 180) then
  !!#!              local_lons(njj) = local_lons(njj) + 360
  !!#!            !else
  !!#!            !  local_lons(njj) = local_lons(njj) + 360
  !!#!            endif 
  !!#!          enddo
  !!#!        endif

  !!#!        trop_loop1: do nii=1,TROP_prep_AI_dims(1)

  !!#!          if(TROP_prep_LON_data(nii) < 0) then
  !!#!            local_lon = TROP_prep_LON_data(nii) + 360
  !!#!          else
  !!#!            local_lon = TROP_prep_LON_data(nii)
  !!#!          endif 

  !!#!          ! Handle case if lat/lon box straddles the prime meridian
  !!#!          if((maxval(local_lons) - minval(local_lons)) > 180) then
  !!#!            if(local_lon <= 180) then
  !!#!              local_lon = local_lon + 360
  !!#!            endif
  !!#!          endif

  !!#!          !write(*,*) TROP_prep_AI_data(nii)
  !!#!          ! Check if the current pixel is within the OMI pixel bounds
  !!#!          ! ---------------------------------------------------------
  !!#!          if(pixel_in_box(LATCRNR_data(:,jj,ii), LONCRNR_data(:,jj,ii), &
  !!#!              TROP_prep_LAT_data(nii), TROP_prep_LON_data(nii))) then

  !!#!              run_trop_total_ai = &
  !!#!                  run_trop_total_ai + TROP_prep_AI_data(nii)
  !!#!              count_trop_ai = count_trop_ai + 1

  !!#!          endif
  !!#!          !enddo trop_loop2
  !!#!        enddo trop_loop1

  !!#!        if(count_trop_ai == 0) then
  !!#!          TROP_out_AI_data(jj,ii)   = -999.
  !!#!        else
  !!#!          TROP_out_AI_data(jj,ii) = run_trop_total_ai / count_trop_ai
  !!#!        endif
  !!#!      endif
  !!#!    enddo omi_loop2
  !!#!  enddo omi_loop1
  !!#!!$omp end parallel

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! NORMAL CODE
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  omi_loop1: do ii=1,LAT_dims(2)

    if(mod(ii, 50) == 0) then
      write(*,*) ii
    endif

    omi_loop2: do jj=1,LAT_dims(1) 

      ! Check if the current pixel is missing
      ! NEW: Check if the current pixel is above minlat
      ! -----------------------------------------------
      if(LAT_data(jj,ii) < min_lat) then
      !if(isnan(AI_data(jj,ii))) then
        TROP_out_AI_data(jj,ii) = -999.
        num_nan = num_nan + 1
      else
        
        count_trop_ai     = 0
        run_trop_total_ai = 0.0

        ! Adjust the lon and lon corners to account for pixels
        ! that straddle the antimeridian
        do njj = 1, 4
          if(LONCRNR_data(njj,jj,ii) < 0) then
            local_lons(njj) = LONCRNR_data(njj,jj,ii) + 360
          else
            local_lons(njj) = LONCRNR_data(njj,jj,ii)     
          endif
        enddo
 
        if((maxval(local_lons) - minval(local_lons)) > 180) then
          do njj = 1, 4
            if(local_lons(njj) <= 180) then
              local_lons(njj) = local_lons(njj) + 360
            !else
            !  local_lons(njj) = local_lons(njj) + 360
            endif 
          enddo
        endif

        trop_loop1: do nii=1,TROP_prep_AI_dims(1)

          if(TROP_prep_LON_data(nii) < 0) then
            local_lon = TROP_prep_LON_data(nii) + 360
          else
            local_lon = TROP_prep_LON_data(nii)
          endif 

          ! Handle case if lat/lon box straddles the prime meridian
          if((maxval(local_lons) - minval(local_lons)) > 180) then
            if(local_lon <= 180) then
              local_lon = local_lon + 360
            endif
          endif

          ! Check if the current pixel is within the OMI pixel bounds
          ! ---------------------------------------------------------
          if(pixel_in_box(LATCRNR_data(:,jj,ii), LONCRNR_data(:,jj,ii), &
              TROP_prep_LAT_data(nii), TROP_prep_LON_data(nii))) then

              run_trop_total_ai = &
                  run_trop_total_ai + TROP_prep_AI_data(nii)
              count_trop_ai = count_trop_ai + 1

          endif
          !enddo trop_loop2
        enddo trop_loop1

        if(count_trop_ai == 0) then
          TROP_out_AI_data(jj,ii)   = -999.
        else
          TROP_out_AI_data(jj,ii) = run_trop_total_ai / count_trop_ai

        endif


      endif
    enddo omi_loop2
  enddo omi_loop1

  write(*,*) istatus, num_nan
  
  ! Open the output file
  ! --------------------
  out_file_name = &
    !'/home/bsorenson/OMI/tropomi_colocate/coloc_data/colocated_tropomi_'//&
    '/home/blake.sorenson/OMI/tropomi_colocate/coloc_data/colocated_tropomi_'//&
    omi_date//'.hdf5'

  write(*,*) trim(out_file_name)

  call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open output file'
    return
  endif


  call write_output_file(out_file_name)


  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write OMI Latitude
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!rank = 2

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_OLT, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'omi_lat', H5T_NATIVE_DOUBLE, dspace_id_OLT,  &
  !!#!                 dset_id_OLT, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'omi_lat'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OLT, H5T_NATIVE_DOUBLE, LAT_data, AI_dims, &
  !!#!                    error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_OLT, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_OLT, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI LAT'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write OMI Longitude
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_OLN, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'omi_lon', H5T_NATIVE_DOUBLE, dspace_id_OLN,  &
  !!#!                 dset_id_OLN, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'omi_lon'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OLN, H5T_NATIVE_DOUBLE, LON_data, AI_dims, &
  !!#!                    error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_OLN, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_OLN, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close output file'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI LON'
  !!#!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write OMI Solar Zenith Angle
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_OSZ, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'omi_sza', H5T_NATIVE_DOUBLE, dspace_id_OSZ,  &
  !!#!                 dset_id_OSZ, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'omi_sza'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OSZ, H5T_NATIVE_DOUBLE, SZA_data, AI_dims, &
  !!#!                    error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_OSZ, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_OSZ, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close output file'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI OSZ'
  !!#!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write OMI Viewing Zenith Angle
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_OVZ, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'omi_vza', H5T_NATIVE_DOUBLE, dspace_id_OVZ,  &
  !!#!                 dset_id_OVZ, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'trop_vza'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OVZ, H5T_NATIVE_DOUBLE, VZA_data, AI_dims, &
  !!#!                    error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_OVZ, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_OVZ, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close output file'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI OVZ'
  !!#!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write OMI Relative Azimuth Angle
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_OAZ, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'omi_azm', H5T_NATIVE_DOUBLE, dspace_id_OAZ,  &
  !!#!                 dset_id_OAZ, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'omi_azm'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OAZ, H5T_NATIVE_DOUBLE, AZM_data, AI_dims, &
  !!#!                    error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_OAZ, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_OAZ, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close output file'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI OAZ'
  !!#!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write OMI AI data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_OAI, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'omi_uvai_pert', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_OAI,  dset_id_OAI, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'omi_uvai_pert'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OAI, H5T_NATIVE_DOUBLE, AI_data, AI_dims, &
  !!#!                    error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_OAI, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_OAI, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI AI'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write TROPOMI AI data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_TAI, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'trop_ai', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_TAI,  dset_id_TAI, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ai'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_TAI, H5T_NATIVE_DOUBLE, TROP_out_AI_data, &
  !!#!                AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_TAI, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_TAI, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote TROPOMI AI'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write TROPOMI SSA0 data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_TS0, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'trop_ssa0', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_TS0,  dset_id_TS0, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa0'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_TS0, H5T_NATIVE_DOUBLE, TROP_out_SSA0_data, &
  !!#!!!#!                AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_TS0, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_TS0, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote TROPOMI SSA0'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write TROPOMI SSA1 data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_TS1, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'trop_ssa1', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_TS1,  dset_id_TS1, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa1'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_TS1, H5T_NATIVE_DOUBLE, TROP_out_SSA1_data, &
  !!#!!!#!                AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_TS1, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_TS0, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote TROPOMI SSA1'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write TROPOMI SSA2 data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_TS2, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'trop_ssa2', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_TS2,  dset_id_TS2, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa2'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_TS2, H5T_NATIVE_DOUBLE, TROP_out_SSA2_data, &
  !!#!!!#!                AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_TS2, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_TS0, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote TROPOMI SSA2'


  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! END OF WRITING
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

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
  call clear_h5_arrays

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
