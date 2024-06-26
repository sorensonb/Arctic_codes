program omi_colocate
!
! NAME:
!   omi_colocate.f90
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
  use comp_vars, only : clear_arrays, i_bad_list, &
    find_distance_between_points, allocate_out_arrays, &
    pixel_in_box, &
    MODIS_CH2_data, MODIS_CH2_dims, &
    MODIS_CH7_data, MODIS_CH7_dims, &
    MODIS_LAT_data, MODIS_LAT_dims, &
    MODIS_LON_data, MODIS_LON_dims, &
    NSIDC_data,     NSIDC_dims, &
    NSIDC_LAT_data, NSIDC_LAT_dims, &
    NSIDC_LON_data, NSIDC_LON_dims, &
    CERES_LWF_data, CERES_LWF_dims, &
    CERES_SWF_data, CERES_SWF_dims, &
    CERES_LAT_data, CERES_LAT_dims, &
    CERES_LON_data, CERES_LON_dims, &
    OMI_AI_data,    OMI_AI_dims, &
    OMI_LAT_data,   OMI_LAT_dims, &
    OMI_LON_data,   OMI_LON_dims, &
    OMI_LATCRNR_data,   OMI_LATCRNR_dims, &
    OMI_LONCRNR_data,   OMI_LONCRNR_dims, &
    OMI_SZA_data,   OMI_SZA_dims, &
    OMI_VZA_data,   OMI_VZA_dims, &
    OMI_AZM_data,   OMI_AZM_dims, &
    MODIS_out_CH2_data, &
    MODIS_out_CH7_data, &
    !MODIS_out_LAT_data, &
    !MODIS_out_LON_data, &
    NSIDC_out_data,     &
    !NSIDC_out_LAT_data, &
    !NSIDC_out_LON_data, &
    CERES_out_LWF_data, &
    CERES_out_SWF_data   
    !CERES_out_LAT_data, &
    !CERES_out_LON_data

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
  integer                :: modis1_file_id       ! id for current HDF5 file
  integer                :: modis2_file_id       ! id for current HDF5 file
  integer                :: nsidc_file_id       ! id for current HDF5 file
  integer                :: ceres_file_id       ! id for current HDF5 file
  integer                :: omi_file_id         ! id for current HDF5 file
  integer                :: out_file_id         ! id for current HDF5 file

  integer                :: dspace_id_OLT  ! OMI LAT
  integer                :: dspace_id_OLN  ! OMI LON
  integer                :: dspace_id_OSZ  ! OMI SZA
  integer                :: dspace_id_OVZ  ! OMI VZA
  integer                :: dspace_id_OAZ  ! OMI AZM
  integer                :: dspace_id_OAI  ! OMI AI
  !integer                :: dspace_id_CLT  ! CERES LAT
  !integer                :: dspace_id_CLN  ! CERES LON
  integer                :: dspace_id_CLW  ! CERES LWF
  integer                :: dspace_id_CSW  ! CERES SWF
  !integer                :: dspace_id_MLT  ! MODIS LAT
  !integer                :: dspace_id_MLN  ! MODIS LON
  integer                :: dspace_id_MC2  ! MODIS CH2
  integer                :: dspace_id_MC7  ! MODIS CH7
  !integer                :: dspace_id_NLT  ! NSIDC LAT
  !integer                :: dspace_id_NLN  ! NSIDC LON
  integer                :: dspace_id_NIC  ! NSIDC ICE

  integer                :: dset_id_OLT  ! OMI LAT
  integer                :: dset_id_OLN  ! OMI LON
  integer                :: dset_id_OSZ  ! OMI SZA
  integer                :: dset_id_OVZ  ! OMI VZA
  integer                :: dset_id_OAZ  ! OMI AZM
  integer                :: dset_id_OAI  ! OMI AI
  !integer                :: dset_id_CLT  ! CERES LAT
  !integer                :: dset_id_CLN  ! CERES LON
  integer                :: dset_id_CLW  ! CERES LWF
  integer                :: dset_id_CSW  ! CERES SWF
  !integer                :: dset_id_MLT  ! MODIS LAT
  !integer                :: dset_id_MLN  ! MODIS LON
  integer                :: dset_id_MC2  ! MODIS CH2
  integer                :: dset_id_MC7  ! MODIS CH7
  !integer                :: dset_id_NLT  ! NSIDC LAT
  !integer                :: dset_id_NLN  ! NSIDC LON
  integer                :: dset_id_NIC  ! NSIDC ICE

  integer                :: rank
  integer(hsize_t), dimension(2)    :: test_dims

  real                   :: distance
  real                   :: closest_dist
  real                   :: min_dist
  real                   :: match_lat
  real                   :: match_lon
  real                   :: match_data1
  real                   :: match_data2

  real                   :: run_modis_total_ch2
  real                   :: run_modis_total_ch7
  integer                :: count_modis_total 

  !!#!real                   :: lat_thresh    ! latitude threshold. Only analyze
  !!#!                                        ! data north of this value.
  !!#!real                   :: lat_gridder   ! Used for setting up the latitude
  !!#!                                        ! grid

  !!#!real,dimension(:), allocatable      :: lat_range ! latitude grid
  !!#!real,dimension(:), allocatable      :: lon_range ! longitude grid
  !!#!real,dimension(:,:), allocatable    :: grids     ! quarter-degree AI grid
  !!#!                                                 ! values.
  !!#!integer,dimension(:,:), allocatable :: i_counts  ! quarter-degree AI counts

  !!#!! File read variables
  !!#!integer,parameter      :: io8    = 42   ! File object for file name file
  !!#!integer,parameter      :: io6    = 1827 ! Data output file
  !!#!integer,parameter      :: errout = 9    ! Error file
  !!#!integer,parameter      :: io10   = 1066 ! Row anomaly file

  !!#!character(len = 255)   :: out_file_name   ! output file name
  !!#!character(len = 255)   :: file_name_file  ! name of file containing the 
  !!#!                                          ! list of HDF5 file names to 
  !!#!                                          ! be analyzed
  !!#!character(len = 255)   :: total_file_name ! file name read from each line 
  !!#!                                          ! of file_name_file
  !!#!character(len = 4)     :: c_work_year     ! holds the previous year

  character(len = 255)      :: modis_name1    ! filename
  character(len = 255)      :: modis_name2    ! filename
  character(len = 255)      :: nsidc_name     ! filename
  character(len = 255)      :: ceres_name     ! filename
  character(len = 255)      :: omi_name       ! filename
  character(len = 255)      :: out_file_name  ! filename

  !!#!! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 5) then
    write(*,*) 'SYNTAX: ./omi_exec modis_name1 modis_name2 nsidc_name ' &
        //'ceres_name omi_name'
    return
  endif

  ! Pull the output file name and input file list file from the command line
  ! ------------------------------------------------------------------------
  call get_command_argument(1,modis_name1)
  call get_command_argument(2,modis_name2)
  call get_command_argument(3,nsidc_name)
  call get_command_argument(4,ceres_name)
  call get_command_argument(5,omi_name)

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"

  ! Open the MODIS Ch2 HDF5 file
  ! ----------------------------
  call h5fopen_f(trim(modis_name1), H5F_ACC_RDWR_F, modis1_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open MODIS CH2 file'
    return
  endif

  ! Open the MODIS Ch7 HDF5 file
  ! ----------------------------
  call h5fopen_f(trim(modis_name2), H5F_ACC_RDWR_F, modis2_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open MODIS CH7 file'
    return
  endif

  ! Open the NSIDC file
  ! -------------------
  call h5fopen_f(trim(nsidc_name), H5F_ACC_RDWR_F, nsidc_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open NSIDC file'
    return
  endif

  ! Open the CERES file
  ! -------------------
  call h5fopen_f(trim(ceres_name), H5F_ACC_RDWR_F, ceres_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open CERES file'
    return
  endif

  ! Open the OMI file
  ! -----------------
  call h5fopen_f(trim(omi_name), H5F_ACC_RDWR_F, omi_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open OMI file'
    return
  endif

  call read_comp_MODIS_CH2(modis1_file_id)
  call read_comp_MODIS_CH7(modis2_file_id)
  call read_comp_MODIS_LAT(modis2_file_id)
  call read_comp_MODIS_LON(modis2_file_id)
  call read_comp_NSIDC_ICE(nsidc_file_id)
  call read_comp_NSIDC_LAT(nsidc_file_id)
  call read_comp_NSIDC_LON(nsidc_file_id)
  call read_comp_CERES_SWF(ceres_file_id)
  call read_comp_CERES_LWF(ceres_file_id)
  call read_comp_CERES_LAT(ceres_file_id)
  call read_comp_CERES_LON(ceres_file_id)
  call read_comp_OMI_AI(omi_file_id)
  call read_comp_OMI_LAT(omi_file_id)
  call read_comp_OMI_LON(omi_file_id)
  call read_comp_OMI_LATCRNR(omi_file_id)
  call read_comp_OMI_LONCRNR(omi_file_id)
  call read_comp_OMI_SZA(omi_file_id)
  call read_comp_OMI_VZA(omi_file_id)
  call read_comp_OMI_AZM(omi_file_id)

  !test_dims = (/10, 20/)
  test_dims = (/OMI_AI_dims(1), OMI_AI_dims(2)/)
  !test_dims = OMI_AI_dims

  ! Allocate the output arrays
  ! --------------------------
  call allocate_out_arrays

  !!#!write(*,*) MODIS_LAT_data(100,100), MODIS_LON_data(100,100)
  !!#!write(*,*) CERES_LAT_data(100,100), CERES_LON_data(100,100)

  !!#!distance = &
  !!#!  find_distance_between_points(MODIS_LAT_data(100,100), &
  !!#!                               MODIS_LON_data(100,100), &
  !!#!                               CERES_LAT_data(100,100), &
  !!#!                               CERES_LON_data(100,100))

  !!#!write(*,*) distance

  istatus = 0
  num_nan = 0 
  min_dist  = 50.
  closest_dist = 999999.
  match_lat  = -999.
  match_lon  = -999.
  match_data1 = -999.
  match_data2 = -999.

  omi_loop1: do ii=1,OMI_LAT_dims(2)

    if(mod(ii, 50) == 0) then
      write(*,*) ii
    endif

    omi_loop2: do jj=1,OMI_LAT_dims(1) 

      ! Check if the current pixel is missing
      ! -------------------------------------
      if(isnan(OMI_AI_data(jj,ii))) then
        !write(*,*) 'NAN VALUE'
        !NSIDC_out_LAT_data(jj,ii) = -999.
        !NSIDC_out_LON_data(jj,ii) = -999.
        NSIDC_out_data(jj,ii)     = -999.
        !CERES_out_LAT_data(jj,ii) = -999.
        !CERES_out_LON_data(jj,ii) = -999.
        CERES_out_LWF_data(jj,ii) = -999.
        CERES_out_SWF_data(jj,ii) = -999.
        !MODIS_out_LAT_data(jj,ii) = -999.
        !MODIS_out_LON_data(jj,ii) = -999.
        MODIS_out_CH2_data(jj,ii) = -999.
        MODIS_out_CH7_data(jj,ii) = -999.
        num_nan = num_nan + 1
      else

        ! Now, loop over the NSIDC data
        ! -----------------------------
        nsidc_loop1: do nii=1,NSIDC_LAT_dims(2)
          nsidc_loop2: do njj=1,NSIDC_LAT_dims(1) 

            ! Calculate the distance between the current pixels
            ! -------------------------------------------------
            distance = find_distance_between_points(OMI_LAT_data(jj,ii), &
                                                    OMI_LON_data(jj,ii), &
                                                    NSIDC_LAT_data(njj,nii), &
                                                    NSIDC_LON_data(njj,nii))

            if(distance < closest_dist) then
              closest_dist = distance
              !NSIDC_out_LAT_data(jj,ii) = NSIDC_LAT_data(njj, nii) 
              !NSIDC_out_LON_data(jj,ii) = NSIDC_LON_data(njj, nii) 
              NSIDC_out_data(jj,ii)     = NSIDC_data(njj, nii) 
            endif

          enddo nsidc_loop2
        enddo nsidc_loop1

        ! Check the distance requirement for the NSIDC pixel
        ! --------------------------------------------------
        if(closest_dist > min_dist) then
          !NSIDC_out_LAT_data(jj,ii) = -999.
          !NSIDC_out_LON_data(jj,ii) = -999.
          NSIDC_out_data(jj,ii)     = -999.
        endif

        closest_dist = 999999.

        ! Now, loop over the CERES data
        ! -----------------------------
        ceres_loop1: do nii=1,CERES_LAT_dims(2)
          ceres_loop2: do njj=1,CERES_LAT_dims(1) 

            ! Calculate the distance between the current pixels
            ! -------------------------------------------------
            distance = find_distance_between_points(OMI_LAT_data(jj,ii), &
                                                    OMI_LON_data(jj,ii), &
                                                    CERES_LAT_data(njj,nii), &
                                                    CERES_LON_data(njj,nii))

            if(distance < closest_dist) then
              closest_dist = distance
              !CERES_out_LAT_data(jj,ii) = CERES_LAT_data(njj, nii) 
              !CERES_out_LON_data(jj,ii) = CERES_LON_data(njj, nii) 
              CERES_out_LWF_data(jj,ii) = CERES_LWF_data(njj, nii) 
              CERES_out_SWF_data(jj,ii) = CERES_SWF_data(njj, nii) 
            endif

          enddo ceres_loop2
        enddo ceres_loop1

        ! Check the distance requirement for the CERES pixel
        ! --------------------------------------------------
        if(closest_dist > min_dist) then
          !CERES_out_LAT_data(jj,ii) = -999.
          !CERES_out_LON_data(jj,ii) = -999.
          CERES_out_LWF_data(jj,ii) = -999.
          CERES_out_SWF_data(jj,ii) = -999.
        endif

        closest_dist = 999999.

        count_modis_total = 0
        run_modis_total_ch2 = 0.0
        run_modis_total_ch7 = 0.0
        
        modis_loop1: do nii=1,MODIS_LAT_dims(2)
          modis_loop2: do njj=1,MODIS_LAT_dims(1) 

            ! Check if the current pixel is within the OMI pixel bounds
            ! ---------------------------------------------------------
            if(pixel_in_box(OMI_LATCRNR_data(:,jj,ii), OMI_LONCRNR_data(:,jj,ii), &
                MODIS_LAT_data(njj,nii), MODIS_LON_data(njj,nii))) then

                run_modis_total_ch2 = &
                    run_modis_total_ch2 + MODIS_CH2_data(njj,nii)
                run_modis_total_ch7 = &
                    run_modis_total_ch7 + MODIS_CH7_data(njj,nii)
                count_modis_total = count_modis_total + 1

            endif
            !!#!! Calculate the distance between the current pixels
            !!#!! -------------------------------------------------
            !!#!distance = find_distance_between_points(OMI_LAT_data(jj,ii), &
            !!#!                                        OMI_LON_data(jj,ii), &
            !!#!                                        MODIS_LAT_data(njj,nii), &
            !!#!                                        MODIS_LON_data(njj,nii))

            !!#!if(distance < closest_dist) then
            !!#!  closest_dist = distance
            !!#!  !MODIS_out_LAT_data(jj,ii) = MODIS_LAT_data(njj, nii) 
            !!#!  !MODIS_out_LON_data(jj,ii) = MODIS_LON_data(njj, nii) 
            !!#!  MODIS_out_CH2_data(jj,ii) = MODIS_CH2_data(njj, nii) 
            !!#!  MODIS_out_CH7_data(jj,ii) = MODIS_CH7_data(njj, nii) 
            !!#!endif

          enddo modis_loop2
        enddo modis_loop1

        !!#!! Check the distance requirement for the MODIS pixel
        !!#!! --------------------------------------------------
        !!#!if(closest_dist > min_dist) then
        !!#!  !MODIS_out_LAT_data(jj,ii) = -999.
        !!#!  !MODIS_out_LON_data(jj,ii) = -999.
        !!#!  MODIS_out_CH2_data(jj,ii) = -999.
        !!#!  MODIS_out_CH7_data(jj,ii) = -999.
        !!#!endif
        if(count_modis_total == 0) then
          !MODIS_out_LAT_data(jj,ii) = -999.
          !MODIS_out_LON_data(jj,ii) = -999.
          MODIS_out_CH2_data(jj,ii) = -999.
          MODIS_out_CH7_data(jj,ii) = -999.
        else
          MODIS_out_CH2_data(jj,ii) = run_modis_total_ch2 / count_modis_total
          MODIS_out_CH7_data(jj,ii) = run_modis_total_ch7 / count_modis_total
        endif

        closest_dist = 999999.

      endif

    enddo omi_loop2
  enddo omi_loop1

  write(*,*) istatus, num_nan
  
  ! Open the output file
  ! --------------------
  write(*,*) omi_name(len(trim(omi_name)) - 16:len(trim(omi_name)) - 5)
  out_file_name = 'colocated_subset_'//omi_name(len(trim(omi_name)) - &
    16:len(trim(omi_name)) - 5)//'.hdf5'
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
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_vza'
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

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write CERES LAT data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_CLT, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'ceres_lat', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_CLT,  dset_id_CLT, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_lat'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_CLT, H5T_NATIVE_DOUBLE, CERES_out_LAT_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_CLT, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_CLT, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote CERES LAT'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write CERES LON data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_CLN, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'ceres_lon', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_CLN,  dset_id_CLN, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_lon'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_CLN, H5T_NATIVE_DOUBLE, CERES_out_LON_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_CLN, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_CLN, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote CERES LON'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write CERES LWF data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_CLW, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'ceres_lwf', H5T_NATIVE_DOUBLE, &
                   dspace_id_CLW,  dset_id_CLW, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_lwf'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_CLW, H5T_NATIVE_DOUBLE, CERES_out_LWF_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_CLW, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_CLW, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote CERES LWF'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write CERES SWF data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_CSW, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'ceres_swf', H5T_NATIVE_DOUBLE, &
                   dspace_id_CSW,  dset_id_CSW, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_swf'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_CSW, H5T_NATIVE_DOUBLE, CERES_out_SWF_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_CSW, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_CSW, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote CERES SWF'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write MODIS LAT data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_MLT, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'modis_lat', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_MLT,  dset_id_MLT, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_lat'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_MLT, H5T_NATIVE_DOUBLE, MODIS_out_LAT_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_MLT, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_MLT, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote MODIS LAT'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write MODIS LON data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_MLN, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'modis_lon', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_MLN,  dset_id_MLN, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_lon'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_MLN, H5T_NATIVE_DOUBLE, MODIS_out_LON_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_MLN, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_MLN, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote MODIS LON'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write MODIS CH2 data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_MC2, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'modis_ch2', H5T_NATIVE_DOUBLE, &
                   dspace_id_MC2,  dset_id_MC2, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'modis_ch2'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_MC2, H5T_NATIVE_DOUBLE, MODIS_out_CH2_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_MC2, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_MC2, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote MODIS CH2'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write MODIS CH7 data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_MC7, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'modis_ch7', H5T_NATIVE_DOUBLE, &
                   dspace_id_MC7,  dset_id_MC7, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'modis_ch7'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_MC7, H5T_NATIVE_DOUBLE, MODIS_out_CH7_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_MC7, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_MC7, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote MODIS CH7'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write NSIDC LAT data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_NLT, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'nsidc_lat', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_NLT,  dset_id_NLT, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_lat'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_NLT, H5T_NATIVE_DOUBLE, NSIDC_out_LAT_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_NLT, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_NLT, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote NSIDC LAT'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write NSIDC LON data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_NLN, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'nsidc_lon', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_NLN,  dset_id_NLN, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_lon'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_NLN, H5T_NATIVE_DOUBLE, NSIDC_out_LON_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_NLN, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_NLN, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote NSIDC LON'

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Write NSIDC NIC data
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(rank, test_dims, dspace_id_NIC, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'nsidc_ice', H5T_NATIVE_DOUBLE, &
                   dspace_id_NIC,  dset_id_NIC, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_ice'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_NIC, H5T_NATIVE_DOUBLE, NSIDC_out_data, &
                  OMI_AI_dims, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  ! Close the dataset
  ! -----------------
  call h5dclose_f(dset_id_NIC, error)

  ! Close access to data space rank
  call h5sclose_f(dspace_id_NIC, error)

  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not write to dataset'
    return
  endif

  write(*,*) 'Wrote NSIDC ice'


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

  ! Close file
  ! ----------
  call h5fclose_f(modis1_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close MODIS CH2 file'
    return
  endif

  call h5fclose_f(modis2_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close MODIS CH7 file'
    return
  endif

  call h5fclose_f(nsidc_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close NSIDC file'
    return
  endif

  call h5fclose_f(ceres_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close CERES file'
    return
  endif

  call h5fclose_f(omi_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close OMI file'
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

end program omi_colocate
