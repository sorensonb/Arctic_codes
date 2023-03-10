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
    TROP_prep_AI_data,    TROP_prep_AI_dims, &
    TROP_prep_SSA0_data,  TROP_prep_SSA0_dims, &
    TROP_prep_SSA1_data,  TROP_prep_SSA1_dims, &
    TROP_prep_SSA2_data,  TROP_prep_SSA2_dims, &
    TROP_prep_LAT_data, TROP_prep_LAT_dims, &
    TROP_prep_LON_data, TROP_prep_LON_dims, &
    OMI_AI_data,    OMI_AI_dims, &
    OMI_LAT_data,   OMI_LAT_dims, &
    OMI_LON_data,   OMI_LON_dims, &
    OMI_LATCRNR_data,   OMI_LATCRNR_dims, &
    OMI_LONCRNR_data,   OMI_LONCRNR_dims, &
    OMI_SZA_data,   OMI_SZA_dims, &
    OMI_VZA_data,   OMI_VZA_dims, &
    OMI_AZM_data,   OMI_AZM_dims, &
    TROP_out_AI_data, TROP_out_SSA0_data, TROP_out_SSA1_data, &
      TROP_out_SSA2_data

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
  integer                :: dspace_id_TS0  ! TROPOMI SSA0
  integer                :: dspace_id_TS1  ! TROPOMI SSA1
  integer                :: dspace_id_TS2  ! TROPOMI SSA2

  integer                :: dset_id_OLT  ! OMI LAT
  integer                :: dset_id_OLN  ! OMI LON
  integer                :: dset_id_OSZ  ! OMI SZA
  integer                :: dset_id_OVZ  ! OMI VZA
  integer                :: dset_id_OAZ  ! OMI AZM
  integer                :: dset_id_OAI  ! OMI AI
  integer                :: dset_id_TAI  ! TROPOMI AI
  integer                :: dset_id_TS0  ! TROPOMI SSA0
  integer                :: dset_id_TS1  ! TROPOMI SSA1
  integer                :: dset_id_TS2  ! TROPOMI SSA2

  integer                :: rank
  integer(hsize_t), dimension(2)    :: test_dims

  real                   :: closest_dist
  real                   :: min_dist
  real                   :: match_lat
  real                   :: match_lon
  real                   :: min_lat
  real                   :: match_data2

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
  character(len = 255)      :: trop_file_name  ! filename
  character(len = 255)      :: out_file_name   ! filename

  !logical                   :: test_logic
  !logical                   :: test_inbox
  !!#!real(kind = 8), dimension(4)  :: lats
  !!#!real(kind = 8), dimension(4)  :: lons
  !!#!real(kind = 8)                :: plat
  !!#!real(kind = 8)                :: plon

  !!#!! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  !!#!lats = [73.5486, 73.4014, 73.2849, 73.4465]
  !!#!lons = [179.7175, 179.6764, 179.9431, -179.4345]

  !!#!! inside
  !!#!plat = 73.4117
  !!#!plon = -179.8443


  !!#!write(*,*) pixel_in_box(lats, lons, plat, plon) 

  !!#!! outside
  !!#!plat = 73.5925
  !!#!plon = -179.2505  
  !!#!write(*,*) pixel_in_box(lats, lons, plat, plon) 
 
  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 2) then
    write(*,*) 'SYNTAX: ./trop_exec omi_date tropomi_date'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1, omi_date)
  call get_command_argument(2, trop_date)

  ! Set the paths to the data files
  ! -------------------------------
  omi_path  = '/home/bsorenson/OMI/arctic_comp/comp_data/'
  trop_path = '/home/bsorenson/OMI/tropomi_colocate/prep_data/'

  omi_file_name  = trim(omi_path)//omi_date(1:8)//'/'//omi_date//'/omi_shawn_'//omi_date//'.hdf5'
  trop_file_name = trim(trop_path)//'tropomi_uvai_coloc_prep_'//trop_date//'.hdf5'

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


  call read_comp_OMI_AI(omi_file_id)
  call read_comp_OMI_LAT(omi_file_id)
  call read_comp_OMI_LON(omi_file_id)
  call read_comp_OMI_LATCRNR(omi_file_id)
  call read_comp_OMI_LONCRNR(omi_file_id)
  call read_comp_OMI_SZA(omi_file_id)
  call read_comp_OMI_VZA(omi_file_id)
  call read_comp_OMI_AZM(omi_file_id)

  call read_prep_TROP_AI(trop_file_id)
  call read_prep_TROP_SSA0(trop_file_id)
  call read_prep_TROP_SSA1(trop_file_id)
  call read_prep_TROP_SSA2(trop_file_id)
  call read_prep_TROP_LAT(trop_file_id)
  call read_prep_TROP_LON(trop_file_id)

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
  !!#!  write(*,*) TROP_prep_LAT_data(ii), TROP_prep_LON_data(ii), TROP_prep_AI_data(ii)
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
        TROP_out_SSA0_data(jj,ii) = -999.
        TROP_out_SSA1_data(jj,ii) = -999.
        TROP_out_SSA2_data(jj,ii) = -999.
        num_nan = num_nan + 1
      else
        
        count_trop_ai     = 0
        count_trop_ssa0   = 0
        count_trop_ssa1   = 0
        count_trop_ssa2   = 0
        run_trop_total_ai = 0.0
        run_trop_total_ssa0 = 0.0
        run_trop_total_ssa1 = 0.0
        run_trop_total_ssa2 = 0.0

        !!#!if((minval(OMI_LONCRNR_data(:,jj,ii)) < 0) .and. &
        !!#!           (maxval(OMI_LONCRNR_data(:,jj,ii)) > 0)) then

        !!#!do njj = 1, 4
        !!#!  if(OMI_LONCRNR_data(njj,jj,ii) < 0) then
        !!#!    local_lon_crnr(njj) = OMI_LONCRNR_data(njj,jj,ii) + 360
        !!#!    !write(*,*) "HERE", OMI_LONCRNR_data(njj,jj,ii), &
        !!#!    !  local_lon_crnr(njj)
        !!#!  else
        !!#!    local_lon_crnr(njj) = OMI_LONCRNR_data(njj,jj,ii)     
        !!#!  endif
        !!#!enddo 

        !!#!if(OMI_LON_data(jj,ii) < 0) then
        !!#!  local_lon = OMI_LON_data(jj,ii) + 360
        !!#!else
        !!#!  local_lon = OMI_LON_data(jj,ii)
        !!#!endif 

        !!#!!local_lon_crnr = OMI_LONCRNR_data(:,jj,ii)
        !!#!!local_lon = OMI_LON_data(jj,ii)
        !!#!write(*,'(5(f8.3,2x),1(1l, 1x))') &
        !!#!           !OMI_LAT_data(jj,ii), OMI_LON_data(jj,ii), &
        !!#!           OMI_LAT_data(jj,ii), local_lon, OMI_LON_data(jj,ii), &
        !!#!           minval(local_lon_crnr) , &
        !!#!           !minval(OMI_LONCRNR_data(:,jj,ii)) , &
        !!#!           maxval(local_lon_crnr), &
        !!#!           !maxval(OMI_LONCRNR_data(:,jj,ii)) - 360, &
        !!#!    !!#!! Logical 1: LON corners straddle antimeridian
        !!#!    !!#!((minval(OMI_LONCRNR_data(:,jj,ii)) < 0) .and. &
        !!#!    !!#!    (maxval(OMI_LONCRNR_data(:,jj,ii)) > 0)), &
        !!#!    !!#!! Logical 2: pixel LAT within the range of the corner lats
        !!#!    !!#!(OMI_LAT_data(jj,ii) <= maxval(OMI_LATCRNR_data(:,jj,ii)) &
        !!#!    !!#!        .and. OMI_LAT_data(jj,ii) >= &
        !!#!    !!#!         minval(OMI_LATCRNR_data(:,jj,ii))), &
        !!#!    !!#!! Logical 3: pixel LON within the range of the corner lons
        !!#!    !!#!((local_lon >= minval(local_lon_crnr(:))) &
        !!#!    !!#!    .and. (local_lon <= maxval(local_lon_crnr(:)))), &
        !!#!    !!#!! Logical 4: combined
        !!#!    !!#!(((minval(OMI_LONCRNR_data(:,jj,ii)) < 0) .and. &
        !!#!    !!#!    (maxval(OMI_LONCRNR_data(:,jj,ii)) > 0)) .and. &
        !!#!    !!#!    ((OMI_LAT_data(jj,ii) <= maxval(OMI_LATCRNR_data(:,jj,ii)) &
        !!#!    !!#!    .and. OMI_LAT_data(jj,ii) >= &
        !!#!    !!#!     minval(OMI_LATCRNR_data(:,jj,ii))) .and. &
        !!#!    !!#!     (local_lon >= minval(local_lon_crnr) &
        !!#!    !!#!    .and. local_lon <= maxval(local_lon_crnr)))), &
        !!#!    !!#!! Logical 5: pixel in box
        !!#!    pixel_in_box(OMI_LATCRNR_data(:,jj,ii), &
        !!#!                 OMI_LONCRNR_data(:,jj,ii), &
        !!#!                 OMI_LAT_data(jj,ii), OMI_LON_data(jj,ii))
        !write(*,'(5x, 4(f8.3, 2x))')  OMI_LONCRNR_data(:,jj,ii)
        !write(*,'(5x, 4(f8.3, 2x))')  local_lon_crnr

        !endif

        !!#!if(((minval(OMI_LONCRNR_data(:,jj,ii)) < 0) .and. &
        !!#!    (maxval(OMI_LONCRNR_data(:,jj,ii)) > 0)) .and. &
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
          !!#!   (OMI_LON_data(jj,ii) > 170 .or. &
          !!#!    OMI_LON_data(jj,ii) < -170)) then


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
          !!#!    ((TROP_prep_LAT_data(nii) <= maxval(OMI_LATCRNR_data(:,jj,ii)) &
          !!#!    .and. TROP_prep_LAT_data(nii) >= &
          !!#!     minval(OMI_LATCRNR_data(:,jj,ii))) .and. &
          !!#!     (local_lon >= minval(local_lons) &
          !!#!    .and. local_lon <= maxval(local_lons)))
          !!#!test_inbox = pixel_in_box(OMI_LATCRNR_data(:,jj,ii), &
          !!#!                 OMI_LONCRNR_data(:,jj,ii), &
          !!#!                 TROP_prep_LAT_data(nii), &
          !!#!                 TROP_prep_LON_data(nii))

          !!#!if(test_logic .neqv. test_inbox) then
          !!#!  write(*,'(a8, 8(f8.3,2x),4(1l, 1x), 3x, 1l)') &
          !!#!             !OMI_LAT_data(jj,ii), OMI_LON_data(jj,ii), &
          !!#!             'ERROR', &
          !!#!             OMI_LAT_data(jj,ii), OMI_LON_data(jj,ii), &
          !!#!             TROP_prep_LAT_data(nii), local_lon, &
          !!#!             minval(OMI_LATCRNR_data(:,jj,ii)) , &
          !!#!             maxval(OMI_LATCRNR_data(:,jj,ii)), &
          !!#!             minval(local_lons(:)) , &
          !!#!             maxval(local_lons(:)), &
          !!#!             !maxval(OMI_LONCRNR_data(:,jj,ii)) - 360, &
          !!#!      ! Logical 1: LON corners straddle antimeridian
          !!#!      ((minval(OMI_LONCRNR_data(:,jj,ii)) < 0) .and. &
          !!#!          (maxval(OMI_LONCRNR_data(:,jj,ii)) > 0)), &
          !!#!      ! Logical 2: pixel LAT within the range of the corner lats
          !!#!      (TROP_prep_LAT_data(nii) <= maxval(OMI_LATCRNR_data(:,jj,ii)) &
          !!#!              .and. TROP_prep_LAT_data(nii) >= &
          !!#!               minval(OMI_LATCRNR_data(:,jj,ii))), &
          !!#!      ! Logical 3: pixel LON within the range of the corner lons
          !!#!      ((local_lon >= minval(local_lons(:))) &
          !!#!          .and. (local_lon <= maxval(local_lons(:)))), &
          !!#!      ! Logical 4: combined
          !!#!      test_logic, &
          !!#!      ! Logical 5: pixel in box
          !!#!      pixel_in_box(OMI_LATCRNR_data(:,jj,ii), &
          !!#!                   OMI_LONCRNR_data(:,jj,ii), &
          !!#!                   TROP_prep_LAT_data(nii), TROP_prep_LON_data(nii))
          !!#!endif
          !!#!endif
          !write(*,'(5x, 4(f8.3, 2x))')  OMI_LONCRNR_data(:,jj,ii)
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

              if(TROP_prep_SSA0_data(nii) > -2e5) then
                run_trop_total_ssa0 = run_trop_total_ssa0 + TROP_prep_SSA0_data(nii)
                count_trop_ssa0 = count_trop_ssa0 + 1
              endif
              if(TROP_prep_SSA1_data(nii) > -2e5) then
                run_trop_total_ssa1 = run_trop_total_ssa1 + TROP_prep_SSA1_data(nii)
                count_trop_ssa1 = count_trop_ssa1 + 1
              endif
              if(TROP_prep_SSA2_data(nii) > -2e5) then
                run_trop_total_ssa2 = run_trop_total_ssa2 + TROP_prep_SSA2_data(nii)
                count_trop_ssa2 = count_trop_ssa2 + 1
              endif

          endif
          !enddo trop_loop2
        enddo trop_loop1

        if(count_trop_ai == 0) then
          TROP_out_AI_data(jj,ii)   = -999.
          TROP_out_SSA0_data(jj,ii) = -999.
          TROP_out_SSA1_data(jj,ii) = -999.
          TROP_out_SSA2_data(jj,ii) = -999.
        else
          TROP_out_AI_data(jj,ii) = run_trop_total_ai / count_trop_ai

          if(count_trop_ssa0 == 0) then
            TROP_out_SSA0_data(jj,ii) = -999.
          else
            TROP_out_SSA0_data(jj,ii) = run_trop_total_ssa0 / count_trop_ssa0
          endif

          if(count_trop_ssa1 == 0) then
            TROP_out_SSA1_data(jj,ii) = -999.
          else
            TROP_out_SSA1_data(jj,ii) = run_trop_total_ssa1 / count_trop_ssa1
          endif

          if(count_trop_ssa2 == 0) then
            TROP_out_SSA2_data(jj,ii) = -999.
          else
            TROP_out_SSA2_data(jj,ii) = run_trop_total_ssa2 / count_trop_ssa2
          endif

        endif


        closest_dist = 999999.

      endif
    enddo omi_loop2
  enddo omi_loop1

  write(*,*) istatus, num_nan
  
  ! Open the output file
  ! --------------------
  out_file_name = &
    '/home/bsorenson/OMI/tropomi_colocate/coloc_data/colocated_tropomi_'//&
    omi_date//'.hdf5'

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
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_lon'
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
  call h5dcreate_f(out_file_id, 'omi_sza', H5T_NATIVE_DOUBLE, dspace_id_OSZ,  &
                   dset_id_OSZ, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_sza'
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
  call h5dcreate_f(out_file_id, 'omi_vza', H5T_NATIVE_DOUBLE, dspace_id_OVZ,  &
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

  write(*,*) 'Wrote TROPOMI AI'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write TROPOMI SSA0 data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_TS0, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'trop_ssa0', H5T_NATIVE_DOUBLE, &
                   dspace_id_TS0,  dset_id_TS0, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa0'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_TS0, H5T_NATIVE_DOUBLE, TROP_out_SSA0_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_TS0, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_TS0, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote TROPOMI SSA0'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write TROPOMI SSA1 data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_TS1, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'trop_ssa1', H5T_NATIVE_DOUBLE, &
                   dspace_id_TS1,  dset_id_TS1, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa1'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_TS1, H5T_NATIVE_DOUBLE, TROP_out_SSA1_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_TS1, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_TS0, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote TROPOMI SSA1'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write TROPOMI SSA2 data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_TS2, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'trop_ssa2', H5T_NATIVE_DOUBLE, &
                   dspace_id_TS2,  dset_id_TS2, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa2'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_TS2, H5T_NATIVE_DOUBLE, TROP_out_SSA2_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_TS2, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_TS0, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote TROPOMI SSA2'


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
