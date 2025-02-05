program trop_colocate_thin
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
    OMI_AI_data,    OMI_AI_dims, &
    OMI_AI_raw_data,    OMI_AI_raw_dims, &
    OMI_LAT_data,   OMI_LAT_dims, &
    OMI_LON_data,   OMI_LON_dims, &
    OMI_LATCRNR_data,   OMI_LATCRNR_dims, &
    OMI_LONCRNR_data,   OMI_LONCRNR_dims, &
    OMI_SZA_data,   OMI_SZA_dims, &
    OMI_VZA_data,   OMI_VZA_dims, &
    OMI_AZM_data,   OMI_AZM_dims, &
    TROP_prep_AI_data,    TROP_prep_AI_dims, &
    !!#!TROP_prep_SSA0_data,  TROP_prep_SSA0_dims, &
    !!#!TROP_prep_SSA1_data,  TROP_prep_SSA1_dims, &
    !!#!TROP_prep_SSA2_data,  TROP_prep_SSA2_dims, &
    TROP_prep_LAT_data, TROP_prep_LAT_dims, &
    TROP_prep_LON_data, TROP_prep_LON_dims, &
    TROP_out_AI_data
    !, TROP_out_SSA0_data
    !!#!, TROP_out_SSA1_data, &
    !!#!  TROP_out_SSA2_data

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
  integer(hsize_t)                :: dspace_id_OAI  ! OMI AI pert
  integer(hsize_t)                :: dspace_id_OAR  ! OMI AI raw
  integer(hsize_t)                :: dspace_id_TAI  ! TROPOMI AI
  integer(hsize_t)                :: dspace_id_TS0  ! TROPOMI SSA0
  integer(hsize_t)                :: dspace_id_TS1  ! TROPOMI SSA1
  integer(hsize_t)                :: dspace_id_TS2  ! TROPOMI SSA2

  integer(hsize_t)                :: dset_id_OLT  ! OMI LAT
  integer(hsize_t)                :: dset_id_OLN  ! OMI LON
  integer(hsize_t)                :: dset_id_OSZ  ! OMI SZA
  integer(hsize_t)                :: dset_id_OVZ  ! OMI VZA
  integer(hsize_t)                :: dset_id_OAZ  ! OMI AZM
  integer(hsize_t)                :: dset_id_OAI  ! OMI AI
  integer(hsize_t)                :: dset_id_OAR  ! OMI AI raw
  integer(hsize_t)                :: dset_id_TAI  ! TROPOMI AI
  integer(hsize_t)                :: dset_id_TS0  ! TROPOMI SSA0
  integer(hsize_t)                :: dset_id_TS1  ! TROPOMI SSA1
  integer(hsize_t)                :: dset_id_TS2  ! TROPOMI SSA2

  integer                :: rank
  integer(hsize_t), dimension(2)    :: test_dims

  integer       :: thread_id
  integer       :: swaths_per_thread
  integer       :: begin_idx
  integer       :: end_idx
  integer       :: num_threads

  real                   :: min_lat

  real                   :: run_trop_total_ai
  real                   :: run_trop_total_ssa0
  real                   :: run_trop_total_ssa1
  real                   :: run_trop_total_ssa2
  integer                :: count_trop_ai
  integer                :: count_trop_ssa0
  integer                :: count_trop_ssa1
  integer                :: count_trop_ssa2
 
  real(kind = 8), dimension(4) :: local_lons
  real(kind = 8)               :: local_lon

  character(len = 255)      :: omi_path        ! filename
  character(len = 255)      :: trop_path       ! filename
  character(len = 12)       :: omi_date        ! filename
  character(len = 12)       :: trop_date       ! filename
  character(len = 255)      :: omi_file_name   ! filename
  character(len = 255)      :: omi_just_name   ! filename
  character(len = 255)      :: trop_file_name  ! filename
  character(len = 255)      :: out_file_name   ! filename

  integer                    :: io7  ! omi file name file

  !logical                   :: test_logic
  !logical                   :: test_inbox
  !!#!real(kind = 8), dimension(4)  :: lats
  !!#!real(kind = 8), dimension(4)  :: lons
  !!#!real(kind = 8)                :: plat
  !!#!real(kind = 8)                :: plon

  !!#!! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
 
  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./trop_exec trop_prep_file omi_comp_file'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1, trop_file_name)
  call get_command_argument(2, omi_file_name)

  !! Set the paths to the data files
  !! -------------------------------
  !omi_path  = '/Research/OMI/H5_files/'
  !!omi_path  = '/home/bsorenson/OMI/arctic_comp/comp_data/'
  !trop_path = '/home/bsorenson/OMI/tropomi_colocate/prep_data/'

  !!#!! Open the OMI name file
  !!#!open(io7, file = trim(omi_name_file), iostat = istatus)
  !!#!if(istatus > 0) then
  !!#!  write(*,*) "ERROR: Problem reading "//trim(omi_name_file)
  !!#!  return 
  !!#!endif
 
  !!#!! Read the OMI file name  
  !!#!read(io7, *, iostat = istatus) omi_just_name
 
  !!#!close(io7)

  !!#!omi_file_name = '/Research/OMI/H5_files/'//trim(omi_just_name)

  !!#!write(*,*) trim(trop_file_name), '  ', trim(omi_file_name)


  !!#!! Set up the out file name based on the OMI timestamp
  !!#!! ---------------------------------------------------
  !!#!omi_date = omi_file_name(len(trim(omi_file_name)) - 46 : &
  !!#!                         len(trim(omi_file_name)) - 43)//&
  !!#!           omi_file_name(len(trim(omi_file_name)) - 41 : &
  !!#!                         len(trim(omi_file_name)) - 38)//&
  !!#!           omi_file_name(len(trim(omi_file_name)) - 36 : &
  !!#!                         len(trim(omi_file_name)) - 33)
  !!#!out_file_name = &
  !!#!  '/home/bsorenson/OMI/tropomi_colocate/coloc_data/colocated_tropomi_'//&
  !!#!  omi_date//'.hdf5'
  !!#!write(*,*) trim(out_file_name)

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


  ! Read data from OMI file
  ! -----------------------
  call read_comp_OMI_AI(omi_file_id)
  call read_comp_OMI_AI_raw(omi_file_id)
  call read_comp_OMI_LAT(omi_file_id)
  call read_comp_OMI_LON(omi_file_id)
  call read_comp_OMI_LATCRNR(omi_file_id)
  call read_comp_OMI_LONCRNR(omi_file_id)
  call read_comp_OMI_SZA(omi_file_id)
  call read_comp_OMI_VZA(omi_file_id)
  call read_comp_OMI_AZM(omi_file_id)

  ! Read data from TROPOMI file
  ! ----------------------------
  call read_prep_TROP_AI(trop_file_id)
  !!#!call read_prep_TROP_SSA0(trop_file_id)
  !!#!call read_prep_TROP_SSA1(trop_file_id)
  !!#!call read_prep_TROP_SSA2(trop_file_id)
  call read_prep_TROP_LAT(trop_file_id)
  call read_prep_TROP_LON(trop_file_id)

  ! Close the files after reading
  ! -----------------------------
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

  !test_dims = (/10, 20/)
  test_dims = (/OMI_AI_dims(1), OMI_AI_dims(2)/)
  !test_dims = AI_dims

  write(*,*) "DATA DIMS:",OMI_AI_dims

  ! Allocate the output arrays
  ! --------------------------
  call allocate_out_arrays(OMI_AI_dims(1), OMI_AI_dims(2))
  write(*,*) "HERE2", OMI_LAT_data(20,20)

  istatus = 0
  num_nan = 0 
  min_lat = 55.
  !min_lat = 65.

  !!#!do ii = 1, 20
  !!#!  write(*,*) TROP_prep_LAT_data(ii), TROP_prep_LON_data(ii), TROP_prep_AI_data(ii)
  !!#!enddo

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! OPENMP CODE
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !$omp parallel private(thread_id, begin_idx, end_idx, run_trop_total_ai, count_trop_ai, local_lon)
    thread_id   = omp_get_thread_num()
    num_threads = omp_get_num_threads()
    swaths_per_thread = (OMI_AI_dims(2) + num_threads - 1) / num_threads
    
    begin_idx = (swaths_per_thread * thread_id) + 1
    end_idx   = swaths_per_thread * (thread_id + 1)
    if(end_idx > OMI_AI_dims(2)) then
      end_idx = end_idx - (end_idx - OMI_AI_dims(2))
    endif
   
    write(*,*) thread_id, begin_idx, end_idx 

    ! Loop over swaths
    omi_loop1: do ii = begin_idx, end_idx 
    !!!omi_loop1: do ii=1,OMI_LAT_dims(2)

      !!#!if(mod(ii, 50) == 0) then
      !!#!  write(*,*) ii
      !!#!endif

      omi_loop2: do jj=1,OMI_LAT_dims(1) 

        ! Check if the current pixel is missing
        ! NEW: Check if the current pixel is above minlat
        ! -----------------------------------------------
        if(OMI_LAT_data(jj,ii) < min_lat) then
        !if(isnan(AI_data(jj,ii))) then
          TROP_out_AI_data(jj,ii) = -999.
          !!#!TROP_out_SSA0_data(jj,ii) = -999.
          !!#!TROP_out_SSA1_data(jj,ii) = -999.
          !!#!TROP_out_SSA2_data(jj,ii) = -999.
          !!#!num_nan = num_nan + 1
        else
          
          count_trop_ai     = 0
          !!#!count_trop_ssa0   = 0
          !!#!count_trop_ssa1   = 0
          !!#!count_trop_ssa2   = 0
          run_trop_total_ai = 0.0
          !!#!run_trop_total_ssa0 = 0.0
          !!#!run_trop_total_ssa1 = 0.0
          !!#!run_trop_total_ssa2 = 0.0

          !!#!if((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
          !!#!           (maxval(LONCRNR_data(:,jj,ii)) > 0)) then

          !!#!do njj = 1, 4
          !!#!  if(LONCRNR_data(njj,jj,ii) < 0) then
          !!#!    local_lon_crnr(njj) = LONCRNR_data(njj,jj,ii) + 360
          !!#!    !write(*,*) "HERE", LONCRNR_data(njj,jj,ii), &
          !!#!    !  local_lon_crnr(njj)
          !!#!  else
          !!#!    local_lon_crnr(njj) = LONCRNR_data(njj,jj,ii)     
          !!#!  endif
          !!#!enddo 

          !!#!if(LON_data(jj,ii) < 0) then
          !!#!  local_lon = LON_data(jj,ii) + 360
          !!#!else
          !!#!  local_lon = LON_data(jj,ii)
          !!#!endif 

          !!#!!local_lon_crnr = LONCRNR_data(:,jj,ii)
          !!#!!local_lon = LON_data(jj,ii)
          !!#!write(*,'(5(f8.3,2x),1(1l, 1x))') &
          !!#!           !LAT_data(jj,ii), LON_data(jj,ii), &
          !!#!           LAT_data(jj,ii), local_lon, LON_data(jj,ii), &
          !!#!           minval(local_lon_crnr) , &
          !!#!           !minval(LONCRNR_data(:,jj,ii)) , &
          !!#!           maxval(local_lon_crnr), &
          !!#!           !maxval(LONCRNR_data(:,jj,ii)) - 360, &
          !!#!    !!#!! Logical 1: LON corners straddle antimeridian
          !!#!    !!#!((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
          !!#!    !!#!    (maxval(LONCRNR_data(:,jj,ii)) > 0)), &
          !!#!    !!#!! Logical 2: pixel LAT within the range of the corner lats
          !!#!    !!#!(LAT_data(jj,ii) <= maxval(LATCRNR_data(:,jj,ii)) &
          !!#!    !!#!        .and. LAT_data(jj,ii) >= &
          !!#!    !!#!         minval(LATCRNR_data(:,jj,ii))), &
          !!#!    !!#!! Logical 3: pixel LON within the range of the corner lons
          !!#!    !!#!((local_lon >= minval(local_lon_crnr(:))) &
          !!#!    !!#!    .and. (local_lon <= maxval(local_lon_crnr(:)))), &
          !!#!    !!#!! Logical 4: combined
          !!#!    !!#!(((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
          !!#!    !!#!    (maxval(LONCRNR_data(:,jj,ii)) > 0)) .and. &
          !!#!    !!#!    ((LAT_data(jj,ii) <= maxval(LATCRNR_data(:,jj,ii)) &
          !!#!    !!#!    .and. LAT_data(jj,ii) >= &
          !!#!    !!#!     minval(LATCRNR_data(:,jj,ii))) .and. &
          !!#!    !!#!     (local_lon >= minval(local_lon_crnr) &
          !!#!    !!#!    .and. local_lon <= maxval(local_lon_crnr)))), &
          !!#!    !!#!! Logical 5: pixel in box
          !!#!    pixel_in_box(LATCRNR_data(:,jj,ii), &
          !!#!                 LONCRNR_data(:,jj,ii), &
          !!#!                 LAT_data(jj,ii), LON_data(jj,ii))
          !write(*,'(5x, 4(f8.3, 2x))')  LONCRNR_data(:,jj,ii)
          !write(*,'(5x, 4(f8.3, 2x))')  local_lon_crnr

          !endif

          !!#!if(((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
          !!#!    (maxval(LONCRNR_data(:,jj,ii)) > 0)) .and. &
          !!#!    ((plat <= maxval(lats) .and. plat >= minval(lats)) .and. &
          !!#!     (plon <= minval(lons) .or. plon >= maxval(lons)))) then
          ! Adjust the lon and lon corners to account for pixels
          ! that straddle the antimeridian
          do njj = 1, 4
            if(OMI_LONCRNR_data(njj,jj,ii) < 0) then
              local_lons(njj) = OMI_LONCRNR_data(njj,jj,ii) + 360
            else
              local_lons(njj) = OMI_LONCRNR_data(njj,jj,ii)     
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

            !!#!if((TROP_prep_LON_data(nii) > 170 .or. &
            !!#!    TROP_prep_LON_data(nii) < -170) .and. &
            !!#!   (LON_data(jj,ii) > 170 .or. &
            !!#!    LON_data(jj,ii) < -170)) then


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

            !!#!test_logic =  &
            !!#!    ((TROP_prep_LAT_data(nii) <= maxval(LATCRNR_data(:,jj,ii)) &
            !!#!    .and. TROP_prep_LAT_data(nii) >= &
            !!#!     minval(LATCRNR_data(:,jj,ii))) .and. &
            !!#!     (local_lon >= minval(local_lons) &
            !!#!    .and. local_lon <= maxval(local_lons)))
            !!#!test_inbox = pixel_in_box(LATCRNR_data(:,jj,ii), &
            !!#!                 LONCRNR_data(:,jj,ii), &
            !!#!                 TROP_prep_LAT_data(nii), &
            !!#!                 TROP_prep_LON_data(nii))

            !!#!if(test_logic .neqv. test_inbox) then
            !!#!  write(*,'(a8, 8(f8.3,2x),4(1l, 1x), 3x, 1l)') &
            !!#!             !LAT_data(jj,ii), LON_data(jj,ii), &
            !!#!             'ERROR', &
            !!#!             LAT_data(jj,ii), LON_data(jj,ii), &
            !!#!             TROP_prep_LAT_data(nii), local_lon, &
            !!#!             minval(LATCRNR_data(:,jj,ii)) , &
            !!#!             maxval(LATCRNR_data(:,jj,ii)), &
            !!#!             minval(local_lons(:)) , &
            !!#!             maxval(local_lons(:)), &
            !!#!             !maxval(LONCRNR_data(:,jj,ii)) - 360, &
            !!#!      ! Logical 1: LON corners straddle antimeridian
            !!#!      ((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
            !!#!          (maxval(LONCRNR_data(:,jj,ii)) > 0)), &
            !!#!      ! Logical 2: pixel LAT within the range of the corner lats
            !!#!      (TROP_prep_LAT_data(nii) <= maxval(LATCRNR_data(:,jj,ii)) &
            !!#!              .and. TROP_prep_LAT_data(nii) >= &
            !!#!               minval(LATCRNR_data(:,jj,ii))), &
            !!#!      ! Logical 3: pixel LON within the range of the corner lons
            !!#!      ((local_lon >= minval(local_lons(:))) &
            !!#!          .and. (local_lon <= maxval(local_lons(:)))), &
            !!#!      ! Logical 4: combined
            !!#!      test_logic, &
            !!#!      ! Logical 5: pixel in box
            !!#!      pixel_in_box(LATCRNR_data(:,jj,ii), &
            !!#!                   LONCRNR_data(:,jj,ii), &
            !!#!                   TROP_prep_LAT_data(nii), TROP_prep_LON_data(nii))
            !!#!endif
            !!#!endif
            !write(*,'(5x, 4(f8.3, 2x))')  LONCRNR_data(:,jj,ii)
            !write(*,'(5x, 4(f8.3, 2x))')  local_lon_crnr

            !!#!endif

            !write(*,*) TROP_prep_AI_data(nii)
            ! Check if the current pixel is within the OMI pixel bounds
            ! ---------------------------------------------------------
            if(pixel_in_box(OMI_LATCRNR_data(:,jj,ii), OMI_LONCRNR_data(:,jj,ii), &
                TROP_prep_LAT_data(nii), TROP_prep_LON_data(nii))) then

                run_trop_total_ai = &
                    run_trop_total_ai + TROP_prep_AI_data(nii)
                count_trop_ai = count_trop_ai + 1

                !!#!if(TROP_prep_SSA0_data(nii) > -2e5) then
                !!#!  run_trop_total_ssa0 = run_trop_total_ssa0 + TROP_prep_SSA0_data(nii)
                !!#!  count_trop_ssa0 = count_trop_ssa0 + 1
                !!#!endif
                !!#!if(TROP_prep_SSA1_data(nii) > -2e5) then
                !!#!  run_trop_total_ssa1 = run_trop_total_ssa1 + TROP_prep_SSA1_data(nii)
                !!#!  count_trop_ssa1 = count_trop_ssa1 + 1
                !!#!endif
                !!#!if(TROP_prep_SSA2_data(nii) > -2e5) then
                !!#!  run_trop_total_ssa2 = run_trop_total_ssa2 + TROP_prep_SSA2_data(nii)
                !!#!  count_trop_ssa2 = count_trop_ssa2 + 1
                !!#!endif

            endif
            !enddo trop_loop2
          enddo trop_loop1

          if(count_trop_ai == 0) then
            TROP_out_AI_data(jj,ii)   = -999.
            !!#!TROP_out_SSA0_data(jj,ii) = -999.
            !!#!TROP_out_SSA1_data(jj,ii) = -999.
            !!#!TROP_out_SSA2_data(jj,ii) = -999.
          else
            TROP_out_AI_data(jj,ii) = run_trop_total_ai / count_trop_ai

            !!#!if(count_trop_ssa0 == 0) then
            !!#!  TROP_out_SSA0_data(jj,ii) = -999.
            !!#!else
            !!#!  TROP_out_SSA0_data(jj,ii) = run_trop_total_ssa0 / count_trop_ssa0
            !!#!endif

            !!#!if(count_trop_ssa1 == 0) then
            !!#!  TROP_out_SSA1_data(jj,ii) = -999.
            !!#!else
            !!#!  TROP_out_SSA1_data(jj,ii) = run_trop_total_ssa1 / count_trop_ssa1
            !!#!endif

            !!#!if(count_trop_ssa2 == 0) then
            !!#!  TROP_out_SSA2_data(jj,ii) = -999.
            !!#!else
            !!#!  TROP_out_SSA2_data(jj,ii) = run_trop_total_ssa2 / count_trop_ssa2
            !!#!endif

          endif

        endif
      enddo omi_loop2
    enddo omi_loop1
  !$omp end parallel

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! NORMAL CODE
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Loop over swaths
  !!#!omi_loop1: do ii=1,OMI_LAT_dims(2)

  !!#!  !!#!if(mod(ii, 50) == 0) then
  !!#!  !!#!  write(*,*) ii
  !!#!  !!#!endif

  !!#!  omi_loop2: do jj=1,OMI_LAT_dims(1) 

  !!#!    ! Check if the current pixel is missing
  !!#!    ! NEW: Check if the current pixel is above minlat
  !!#!    ! -----------------------------------------------
  !!#!    if(OMI_LAT_data(jj,ii) < min_lat) then
  !!#!    !if(isnan(AI_data(jj,ii))) then
  !!#!      TROP_out_AI_data(jj,ii) = -999.
  !!#!      !!#!TROP_out_SSA0_data(jj,ii) = -999.
  !!#!      !!#!TROP_out_SSA1_data(jj,ii) = -999.
  !!#!      !!#!TROP_out_SSA2_data(jj,ii) = -999.
  !!#!      !!#!num_nan = num_nan + 1
  !!#!    else
  !!#!      
  !!#!      count_trop_ai     = 0
  !!#!      !!#!count_trop_ssa0   = 0
  !!#!      !!#!count_trop_ssa1   = 0
  !!#!      !!#!count_trop_ssa2   = 0
  !!#!      run_trop_total_ai = 0.0
  !!#!      !!#!run_trop_total_ssa0 = 0.0
  !!#!      !!#!run_trop_total_ssa1 = 0.0
  !!#!      !!#!run_trop_total_ssa2 = 0.0

  !!#!      !!#!if((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
  !!#!      !!#!           (maxval(LONCRNR_data(:,jj,ii)) > 0)) then

  !!#!      !!#!do njj = 1, 4
  !!#!      !!#!  if(LONCRNR_data(njj,jj,ii) < 0) then
  !!#!      !!#!    local_lon_crnr(njj) = LONCRNR_data(njj,jj,ii) + 360
  !!#!      !!#!    !write(*,*) "HERE", LONCRNR_data(njj,jj,ii), &
  !!#!      !!#!    !  local_lon_crnr(njj)
  !!#!      !!#!  else
  !!#!      !!#!    local_lon_crnr(njj) = LONCRNR_data(njj,jj,ii)     
  !!#!      !!#!  endif
  !!#!      !!#!enddo 

  !!#!      !!#!if(LON_data(jj,ii) < 0) then
  !!#!      !!#!  local_lon = LON_data(jj,ii) + 360
  !!#!      !!#!else
  !!#!      !!#!  local_lon = LON_data(jj,ii)
  !!#!      !!#!endif 

  !!#!      !!#!!local_lon_crnr = LONCRNR_data(:,jj,ii)
  !!#!      !!#!!local_lon = LON_data(jj,ii)
  !!#!      !!#!write(*,'(5(f8.3,2x),1(1l, 1x))') &
  !!#!      !!#!           !LAT_data(jj,ii), LON_data(jj,ii), &
  !!#!      !!#!           LAT_data(jj,ii), local_lon, LON_data(jj,ii), &
  !!#!      !!#!           minval(local_lon_crnr) , &
  !!#!      !!#!           !minval(LONCRNR_data(:,jj,ii)) , &
  !!#!      !!#!           maxval(local_lon_crnr), &
  !!#!      !!#!           !maxval(LONCRNR_data(:,jj,ii)) - 360, &
  !!#!      !!#!    !!#!! Logical 1: LON corners straddle antimeridian
  !!#!      !!#!    !!#!((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
  !!#!      !!#!    !!#!    (maxval(LONCRNR_data(:,jj,ii)) > 0)), &
  !!#!      !!#!    !!#!! Logical 2: pixel LAT within the range of the corner lats
  !!#!      !!#!    !!#!(LAT_data(jj,ii) <= maxval(LATCRNR_data(:,jj,ii)) &
  !!#!      !!#!    !!#!        .and. LAT_data(jj,ii) >= &
  !!#!      !!#!    !!#!         minval(LATCRNR_data(:,jj,ii))), &
  !!#!      !!#!    !!#!! Logical 3: pixel LON within the range of the corner lons
  !!#!      !!#!    !!#!((local_lon >= minval(local_lon_crnr(:))) &
  !!#!      !!#!    !!#!    .and. (local_lon <= maxval(local_lon_crnr(:)))), &
  !!#!      !!#!    !!#!! Logical 4: combined
  !!#!      !!#!    !!#!(((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
  !!#!      !!#!    !!#!    (maxval(LONCRNR_data(:,jj,ii)) > 0)) .and. &
  !!#!      !!#!    !!#!    ((LAT_data(jj,ii) <= maxval(LATCRNR_data(:,jj,ii)) &
  !!#!      !!#!    !!#!    .and. LAT_data(jj,ii) >= &
  !!#!      !!#!    !!#!     minval(LATCRNR_data(:,jj,ii))) .and. &
  !!#!      !!#!    !!#!     (local_lon >= minval(local_lon_crnr) &
  !!#!      !!#!    !!#!    .and. local_lon <= maxval(local_lon_crnr)))), &
  !!#!      !!#!    !!#!! Logical 5: pixel in box
  !!#!      !!#!    pixel_in_box(LATCRNR_data(:,jj,ii), &
  !!#!      !!#!                 LONCRNR_data(:,jj,ii), &
  !!#!      !!#!                 LAT_data(jj,ii), LON_data(jj,ii))
  !!#!      !write(*,'(5x, 4(f8.3, 2x))')  LONCRNR_data(:,jj,ii)
  !!#!      !write(*,'(5x, 4(f8.3, 2x))')  local_lon_crnr

  !!#!      !endif

  !!#!      !!#!if(((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
  !!#!      !!#!    (maxval(LONCRNR_data(:,jj,ii)) > 0)) .and. &
  !!#!      !!#!    ((plat <= maxval(lats) .and. plat >= minval(lats)) .and. &
  !!#!      !!#!     (plon <= minval(lons) .or. plon >= maxval(lons)))) then
  !!#!      ! Adjust the lon and lon corners to account for pixels
  !!#!      ! that straddle the antimeridian
  !!#!      do njj = 1, 4
  !!#!        if(OMI_LONCRNR_data(njj,jj,ii) < 0) then
  !!#!          local_lons(njj) = OMI_LONCRNR_data(njj,jj,ii) + 360
  !!#!        else
  !!#!          local_lons(njj) = OMI_LONCRNR_data(njj,jj,ii)     
  !!#!        endif
  !!#!      enddo
 
  !!#!      if((maxval(local_lons) - minval(local_lons)) > 180) then
  !!#!        do njj = 1, 4
  !!#!          if(local_lons(njj) <= 180) then
  !!#!            local_lons(njj) = local_lons(njj) + 360
  !!#!          !else
  !!#!          !  local_lons(njj) = local_lons(njj) + 360
  !!#!          endif 
  !!#!        enddo
  !!#!      endif

  !!#!      trop_loop1: do nii=1,TROP_prep_AI_dims(1)

  !!#!        !!#!if((TROP_prep_LON_data(nii) > 170 .or. &
  !!#!        !!#!    TROP_prep_LON_data(nii) < -170) .and. &
  !!#!        !!#!   (LON_data(jj,ii) > 170 .or. &
  !!#!        !!#!    LON_data(jj,ii) < -170)) then


  !!#!        if(TROP_prep_LON_data(nii) < 0) then
  !!#!          local_lon = TROP_prep_LON_data(nii) + 360
  !!#!        else
  !!#!          local_lon = TROP_prep_LON_data(nii)
  !!#!        endif 

  !!#!        ! Handle case if lat/lon box straddles the prime meridian
  !!#!        if((maxval(local_lons) - minval(local_lons)) > 180) then
  !!#!          if(local_lon <= 180) then
  !!#!            local_lon = local_lon + 360
  !!#!          endif
  !!#!        endif

  !!#!        !!#!test_logic =  &
  !!#!        !!#!    ((TROP_prep_LAT_data(nii) <= maxval(LATCRNR_data(:,jj,ii)) &
  !!#!        !!#!    .and. TROP_prep_LAT_data(nii) >= &
  !!#!        !!#!     minval(LATCRNR_data(:,jj,ii))) .and. &
  !!#!        !!#!     (local_lon >= minval(local_lons) &
  !!#!        !!#!    .and. local_lon <= maxval(local_lons)))
  !!#!        !!#!test_inbox = pixel_in_box(LATCRNR_data(:,jj,ii), &
  !!#!        !!#!                 LONCRNR_data(:,jj,ii), &
  !!#!        !!#!                 TROP_prep_LAT_data(nii), &
  !!#!        !!#!                 TROP_prep_LON_data(nii))

  !!#!        !!#!if(test_logic .neqv. test_inbox) then
  !!#!        !!#!  write(*,'(a8, 8(f8.3,2x),4(1l, 1x), 3x, 1l)') &
  !!#!        !!#!             !LAT_data(jj,ii), LON_data(jj,ii), &
  !!#!        !!#!             'ERROR', &
  !!#!        !!#!             LAT_data(jj,ii), LON_data(jj,ii), &
  !!#!        !!#!             TROP_prep_LAT_data(nii), local_lon, &
  !!#!        !!#!             minval(LATCRNR_data(:,jj,ii)) , &
  !!#!        !!#!             maxval(LATCRNR_data(:,jj,ii)), &
  !!#!        !!#!             minval(local_lons(:)) , &
  !!#!        !!#!             maxval(local_lons(:)), &
  !!#!        !!#!             !maxval(LONCRNR_data(:,jj,ii)) - 360, &
  !!#!        !!#!      ! Logical 1: LON corners straddle antimeridian
  !!#!        !!#!      ((minval(LONCRNR_data(:,jj,ii)) < 0) .and. &
  !!#!        !!#!          (maxval(LONCRNR_data(:,jj,ii)) > 0)), &
  !!#!        !!#!      ! Logical 2: pixel LAT within the range of the corner lats
  !!#!        !!#!      (TROP_prep_LAT_data(nii) <= maxval(LATCRNR_data(:,jj,ii)) &
  !!#!        !!#!              .and. TROP_prep_LAT_data(nii) >= &
  !!#!        !!#!               minval(LATCRNR_data(:,jj,ii))), &
  !!#!        !!#!      ! Logical 3: pixel LON within the range of the corner lons
  !!#!        !!#!      ((local_lon >= minval(local_lons(:))) &
  !!#!        !!#!          .and. (local_lon <= maxval(local_lons(:)))), &
  !!#!        !!#!      ! Logical 4: combined
  !!#!        !!#!      test_logic, &
  !!#!        !!#!      ! Logical 5: pixel in box
  !!#!        !!#!      pixel_in_box(LATCRNR_data(:,jj,ii), &
  !!#!        !!#!                   LONCRNR_data(:,jj,ii), &
  !!#!        !!#!                   TROP_prep_LAT_data(nii), TROP_prep_LON_data(nii))
  !!#!        !!#!endif
  !!#!        !!#!endif
  !!#!        !write(*,'(5x, 4(f8.3, 2x))')  LONCRNR_data(:,jj,ii)
  !!#!        !write(*,'(5x, 4(f8.3, 2x))')  local_lon_crnr

  !!#!        !!#!endif

  !!#!        !write(*,*) TROP_prep_AI_data(nii)
  !!#!        ! Check if the current pixel is within the OMI pixel bounds
  !!#!        ! ---------------------------------------------------------
  !!#!        if(pixel_in_box(OMI_LATCRNR_data(:,jj,ii), OMI_LONCRNR_data(:,jj,ii), &
  !!#!            TROP_prep_LAT_data(nii), TROP_prep_LON_data(nii))) then

  !!#!            run_trop_total_ai = &
  !!#!                run_trop_total_ai + TROP_prep_AI_data(nii)
  !!#!            count_trop_ai = count_trop_ai + 1

  !!#!            !!#!if(TROP_prep_SSA0_data(nii) > -2e5) then
  !!#!            !!#!  run_trop_total_ssa0 = run_trop_total_ssa0 + TROP_prep_SSA0_data(nii)
  !!#!            !!#!  count_trop_ssa0 = count_trop_ssa0 + 1
  !!#!            !!#!endif
  !!#!            !!#!if(TROP_prep_SSA1_data(nii) > -2e5) then
  !!#!            !!#!  run_trop_total_ssa1 = run_trop_total_ssa1 + TROP_prep_SSA1_data(nii)
  !!#!            !!#!  count_trop_ssa1 = count_trop_ssa1 + 1
  !!#!            !!#!endif
  !!#!            !!#!if(TROP_prep_SSA2_data(nii) > -2e5) then
  !!#!            !!#!  run_trop_total_ssa2 = run_trop_total_ssa2 + TROP_prep_SSA2_data(nii)
  !!#!            !!#!  count_trop_ssa2 = count_trop_ssa2 + 1
  !!#!            !!#!endif

  !!#!        endif
  !!#!        !enddo trop_loop2
  !!#!      enddo trop_loop1

  !!#!      if(count_trop_ai == 0) then
  !!#!        TROP_out_AI_data(jj,ii)   = -999.
  !!#!        !!#!TROP_out_SSA0_data(jj,ii) = -999.
  !!#!        !!#!TROP_out_SSA1_data(jj,ii) = -999.
  !!#!        !!#!TROP_out_SSA2_data(jj,ii) = -999.
  !!#!      else
  !!#!        TROP_out_AI_data(jj,ii) = run_trop_total_ai / count_trop_ai

  !!#!        !!#!if(count_trop_ssa0 == 0) then
  !!#!        !!#!  TROP_out_SSA0_data(jj,ii) = -999.
  !!#!        !!#!else
  !!#!        !!#!  TROP_out_SSA0_data(jj,ii) = run_trop_total_ssa0 / count_trop_ssa0
  !!#!        !!#!endif

  !!#!        !!#!if(count_trop_ssa1 == 0) then
  !!#!        !!#!  TROP_out_SSA1_data(jj,ii) = -999.
  !!#!        !!#!else
  !!#!        !!#!  TROP_out_SSA1_data(jj,ii) = run_trop_total_ssa1 / count_trop_ssa1
  !!#!        !!#!endif

  !!#!        !!#!if(count_trop_ssa2 == 0) then
  !!#!        !!#!  TROP_out_SSA2_data(jj,ii) = -999.
  !!#!        !!#!else
  !!#!        !!#!  TROP_out_SSA2_data(jj,ii) = run_trop_total_ssa2 / count_trop_ssa2
  !!#!        !!#!endif

  !!#!      endif
  !!#!    endif
  !!#!  enddo omi_loop2
  !!#!enddo omi_loop1

  !!#!write(*,*) istatus!, num_nan
  
  ! Open the output file
  ! --------------------
  out_file_name = &
    !'/home/bsorenson/OMI/tropomi_colocate/coloc_data/colocated_tropomi_'//&
    !'/home/blake.sorenson/OMI/tropomi_colocate/coloc_data/colocated_tropomi_'//&
    !'colocated_tropomi_'//&
    '/home/blake.sorenson/CSCI_544/final_project/coloc_data/colocated_tropomi_'//&
    omi_file_name(len(trim(omi_file_name)) - &
    16:len(trim(omi_file_name)) - 5)//'.hdf5'

  write(*,*) trim(out_file_name)

  call write_coloc_output_file(out_file_name)

  write(*,*) 'Saved output file'//trim(out_file_name)
 
  ! Deallocate all the arrays for the next pass
  ! -------------------------------------------
  call clear_arrays
  !call clear_h5_arrays

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program trop_colocate_thin
