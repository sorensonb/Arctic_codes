program gen_comp_grid_climo
!
! NAME:
!   gen_comp_grid_climo.f90
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
  !use comp_vars, only : clear_arrays, i_bad_list, &
  use comp_grid_vars, only : clear_arrays, i_bad_list, &
    find_distance_between_points, allocate_out_arrays, &
    pixel_in_box, &
    MODIS_CH2_data, MODIS_CH2_dims, &
    MODIS_CH7_data, MODIS_CH7_dims, &
    !MODIS_LAT_data, MODIS_LAT_dims, &
    !MODIS_LON_data, MODIS_LON_dims, &
    NSIDC_data,     NSIDC_dims, &
    !NSIDC_LAT_data, NSIDC_LAT_dims, &
    !NSIDC_LON_data, NSIDC_LON_dims, &
    CERES_LWF_data, CERES_LWF_dims, &
    CERES_SWF_data, CERES_SWF_dims, &
    !CERES_LAT_data, CERES_LAT_dims, &
    !CERES_LON_data, CERES_LON_dims, &
    OMI_AI_data,    OMI_AI_dims, &
    OMI_AI_raw_data,    OMI_AI_raw_dims, &
    OMI_LAT_data,   OMI_LAT_dims, &
    OMI_LON_data,   OMI_LON_dims, &
    !OMI_LATCRNR_data,   OMI_LATCRNR_dims, &
    !OMI_LONCRNR_data,   OMI_LONCRNR_dims, &
    OMI_SZA_data,   OMI_SZA_dims, &
    OMI_VZA_data,   OMI_VZA_dims, &
    OMI_AZM_data,   OMI_AZM_dims, &
    TROP_AI_data,   TROP_AI_dims, &
    TROP_SSA0_data, TROP_SSA0_dims, &
    TROP_SSA1_data, TROP_SSA1_dims, &
    TROP_SSA2_data, TROP_SSA2_dims, &
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
  integer                :: kk            ! loop counter
  integer                :: nn            ! loop counter
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
  integer                :: file_id
  integer                :: modis1_file_id       ! id for current HDF5 file
  integer                :: modis2_file_id       ! id for current HDF5 file
  integer                :: nsidc_file_id       ! id for current HDF5 file
  integer                :: ceres_file_id       ! id for current HDF5 file
  integer                :: omi_file_id         ! id for current HDF5 file
  integer                :: trop_file_id        ! id for current HDF5 file
  integer                :: out_file_id         ! id for current HDF5 file

  integer                :: dspace_id_OLT  ! OMI LAT
  integer                :: dspace_id_OLN  ! OMI LON
  integer                :: dspace_id_OSZ  ! OMI SZA
  integer                :: dspace_id_OVZ  ! OMI VZA
  integer                :: dspace_id_OAZ  ! OMI AZM
  integer                :: dspace_id_OAI  ! OMI AI pert
  integer                :: dspace_id_OAR  ! OMI AI raw
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
  integer                :: dspace_id_TAI  ! TROPOMI AI
  integer                :: dspace_id_TS0  ! TROPOMI SSA0
  integer                :: dspace_id_TS1  ! TROPOMI SSA1
  integer                :: dspace_id_TS2  ! TROPOMI SSA2

  integer                :: dset_id_OLT  ! OMI LAT
  integer                :: dset_id_OLN  ! OMI LON
  integer                :: dset_id_OSZ  ! OMI SZA
  integer                :: dset_id_OVZ  ! OMI VZA
  integer                :: dset_id_OAZ  ! OMI AZM
  integer                :: dset_id_OAI  ! OMI AI pert
  integer                :: dset_id_OAR  ! OMI AI raw
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
  integer                :: dset_id_TAI  ! TROPOMI AI
  integer                :: dset_id_TS0  ! TROPOMI SSA0
  integer                :: dset_id_TS1  ! TROPOMI SSA1
  integer                :: dset_id_TS2  ! TROPOMI SSA2

  real                   :: distance
  real                   :: closest_dist
  real                   :: min_dist
  real                   :: match_lat
  real                   :: match_lon
  real                   :: match_data1
  real                   :: match_data2
  real                   :: local_trop_ai

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

  !!#!character(len = 255)   :: file_name_file  ! name of file containing the 
  !!#!                                          ! list of HDF5 file names to 
  !!#!                                          ! be analyzed
  !!#!character(len = 255)   :: total_file_name ! file name read from each line 
  !!#!                                          ! of file_name_file
  !!#!character(len = 4)     :: c_work_year     ! holds the previous year

  character(len = 255)      :: file_name_file ! filename
  character(len = 255)      :: total_file_name
  character(len = 12)       :: single_file_date
  character(len = 255)      :: out_file_name   ! output file name

  logical                   :: l_trop_found

  real, dimension(:), allocatable       :: ai_bins
  real, dimension(:), allocatable       :: sza_bins
  real, dimension(:), allocatable       :: ice_bins
  real, dimension(:), allocatable       :: ch7_bins

  real, dimension(:), allocatable       :: ai_edges
  real, dimension(:), allocatable       :: sza_edges
  real, dimension(:), allocatable       :: ice_edges
  real, dimension(:), allocatable       :: ch7_edges

  real, dimension(:,:,:,:), allocatable          :: grid_swf_climo
  integer, dimension(:,:,:,:), allocatable       :: grid_swf_count

  real                      :: min_ai
  real                      :: max_ai
  real                      :: delta_ai
  real                      :: min_sza
  real                      :: max_sza
  real                      :: delta_sza
  real                      :: min_ice
  real                      :: max_ice
  real                      :: delta_ice
  real                      :: min_ch7
  real                      :: max_ch7
  real                      :: delta_ch7

  integer                   :: len_ai
  integer                   :: len_sza
  integer                   :: len_ice
  integer                   :: len_ch7
  integer                   :: num_data_grids

  integer                   :: ai_idx
  integer                   :: sza_idx
  integer                   :: ice_idx
  integer                   :: ch7_idx

  integer                           :: rank
  integer(hsize_t), dimension(4)    :: grid_dims
  integer(hsize_t), dimension(1)    :: ai_dims
  integer(hsize_t), dimension(1)    :: sza_dims
  integer(hsize_t), dimension(1)    :: ice_dims
  integer(hsize_t), dimension(1)    :: ch7_dims

  integer(hsize_t), dimension(1)    :: ai_edge_dims
  integer(hsize_t), dimension(1)    :: sza_edge_dims
  integer(hsize_t), dimension(1)    :: ice_edge_dims
  integer(hsize_t), dimension(1)    :: ch7_edge_dims

  integer                   :: io8
  integer                   :: io5

  ! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  io8 = 1999
  io5 = 2007
  
  ! Check command line arguments
  ! ----------------------------
  arg_count = command_argument_count()
  if(arg_count /= 1) then
    write(*,*) 'SYNTAX: ./comp_grid_exec file_name_file'
    write(*,*) '         file_name_file consists of lines of YYYYMMDDHHMM'
    return
  endif

  call get_command_argument(1,file_name_file)

  write(*,*) file_name_file

  ! Initialize the HDF5 interface
  ! -----------------------------
  call h5open_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not open HDF5 library'
    return
  endif
  write(*,*) "Interface opened"


  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Create gridded climatology bins
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Specify the bin ranges and sizes
  ! --------------------------------
  min_ai   = -2.0
  max_ai   = 12.0
  delta_ai = 0.25

  min_sza   = 35.0
  max_sza   = 85.0
  delta_sza = 2.5

  min_ice   = 0.0
  max_ice   = 105.0
  delta_ice = 5.0

  min_ch7   = 0.0
  max_ch7   = 0.350
  delta_ch7 = 0.025

  ! Determine the number of bins for each variable
  ! ----------------------------------------------
  len_ai = (max_ai - min_ai + delta_ai) / delta_ai
  len_sza = (max_sza - min_sza + delta_sza) / delta_sza
  len_ice = (max_ice - min_ice + delta_ice) / delta_ice
  len_ch7 = (max_ch7 - min_ch7 + delta_ch7) / delta_ch7

  !write(*,*) len_ai
  !write(*,*) len_sza
  !write(*,*) len_ice
  !write(*,*) len_ch7

  allocate(ai_bins(len_ai))
  allocate(sza_bins(len_sza))
  allocate(ice_bins(len_ice))
  allocate(ch7_bins(len_ch7))

  allocate(ai_edges(len_ai + 1))
  allocate(sza_edges(len_sza + 1))
  allocate(ice_edges(len_ice + 1))
  allocate(ch7_edges(len_ch7 + 1))
 
  ch7_bins(:) = -99.

  do ii = 1, len_ai
    ai_bins(ii) = min_ai + delta_ai*(ii-1)
  enddo

  do ii = 1, len_sza
    sza_bins(ii) = min_sza + delta_sza*(ii-1)
  enddo
 
  do ii = 1, len_ice
    ice_bins(ii) = min_ice + delta_ice*(ii-1)
  enddo

  do ii = 1, len_ch7
    ch7_bins(ii) = min_ch7 + delta_ch7*(ii-1)
  enddo

  do ii = 1, len_ai + 1 
    ai_edges(ii) = min_ai + delta_ai*((ii-1) - 0.5)
  enddo

  do ii = 1, len_sza + 1 
    sza_edges(ii) = min_sza + delta_sza*((ii-1) - 0.5)
  enddo
 
  do ii = 1, len_ice + 1 
    ice_edges(ii) = min_ice + delta_ice*((ii-1) - 0.5)
  enddo

  do ii = 1, len_ch7 + 1 
    ch7_edges(ii) = min_ch7 + delta_ch7*((ii-1) - 0.5)
  enddo

  write(*,*) 'total_size = ',len_ai * len_sza * len_ice * len_ch7

  ai_dims   = (/size(ai_bins)/)
  sza_dims  = (/size(sza_bins)/)
  ice_dims  = (/size(ice_bins)/)
  ch7_dims  = (/size(ch7_bins)/)
  grid_dims = (/size(ai_bins),size(sza_bins),size(ice_bins),size(ch7_bins)/)

  ai_edge_dims  = (/size(ai_edges)/)
  sza_edge_dims = (/size(sza_edges)/)
  ice_edge_dims = (/size(ice_edges)/)
  ch7_edge_dims = (/size(ch7_edges)/)

  write(*,*) ai_dims
  write(*,*) ai_bins

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Allocate output gridded climatology arrays
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  allocate(grid_swf_climo(len_ai, len_sza, len_ice, len_ch7))
  allocate(grid_swf_count(len_ai, len_sza, len_ice, len_ch7))

  grid_swf_climo(:,:,:,:) = -999.
  grid_swf_count(:,:,:,:) = -9


  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Insert data from combined_subset_YYYYMMDDHHMM.hdf5 files into the 
  ! climatology grid.
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! Open the file name file
  ! -----------------------
  open(io8, file = trim(file_name_file), iostat = istatus)
  if(istatus > 0) then
    write(*,*) "ERROR: Problem reading "//trim(file_name_file)
    return 
  else
    ! Loop over the file name file
    ! ----------------------------
    file_loop: do
      ! Read the current total_file_name from the file
      ! ----------------------------------------------
      read(io8, '(A)', iostat=istatus) single_file_date
      if(istatus < 0) then 
        write(*,*) "End of "//trim(file_name_file)//" found"
        exit
      else if(istatus > 0) then
        write(*,*) "ERROR: problem reading total_file_name"
        cycle file_loop
      else
        write(*,*) single_file_date

        ! Build the combined subset file name
        ! -----------------------------------
        total_file_name = '/home/bsorenson/OMI/arctic_comp/comp_data/'&
          //'colocated_subset_'//single_file_date//'.hdf5'

        ! Try to open the file
        ! --------------------
        ! Open the HDF5 file
        ! ------------------
        call h5fopen_f(total_file_name, H5F_ACC_RDWR_F, file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: could not open file '//total_file_name
          cycle file_loop
        endif
       
        !write(*,*) 'File ',trim(total_file_name),' has been opened'

        ! Read in the variables from this file
        ! As of 2023/05/09, these are the available variables from the 
        ! colocated_subset files:
        !  ceres_lwf
        !  x ceres_swf
        !  modis_ch2
        !  x modis_ch7
        !  x nsidc_ice
        !  omi_azm
        !  omi_lat
        !  omi_lon
        !  x omi_sza
        !  omi_uvai_pert
        !  x omi_uvai_raw
        !  omi_vza
        !  trop_ssa0
        !  trop_ssa1
        !  trop_ssa2
        !  trop_uvai
        ! ------------------------------------
        call read_coloc_MODIS_CH7(file_id)
        call read_coloc_CERES_SWF(file_id)
        call read_coloc_NSIDC_ICE(file_id)
        call read_coloc_OMI_AI_raw(file_id)
        call read_coloc_OMI_SZA(file_id)

        !!write(*,*) sza_bins, OMI_SZA_data(12,120)
        !ai_idx =  minloc(abs(sza_bins - OMI_SZA_data(12, 120)), dim = 1)    
        !write(*,*) 'SZA_IDX = ',ai_idx
        !write(*,*) minloc(abs(sza_bins - OMI_SZA_data(12, 120)), dim = 1), &
        !           sza_bins(minloc(abs(sza_bins - OMI_SZA_data(12, 120))))

        omi_loop1: do ii=1,OMI_AI_raw_dims(2)
          !if(mod(ii, 50) == 0) then
          !  write(*,*) ii
          !endif

          omi_loop2: do jj=1,OMI_AI_raw_dims(1) 

            ! Insert each pixel into the grid here
            ! ------------------------------------      
            ai_idx  = minloc(abs(ai_bins - OMI_AI_raw_data(jj,ii)), dim = 1)
            sza_idx = minloc(abs(sza_bins - OMI_SZA_data(jj,ii)), dim = 1)
            ice_idx = minloc(abs(ice_bins - NSIDC_data(jj,ii)), dim = 1)
            ch7_idx = minloc(abs(ch7_bins - MODIS_CH7_data(jj,ii)), dim = 1)

            if(CERES_SWF_data(jj,ii) > 0) then
              if(grid_swf_count(ai_idx, sza_idx,ice_idx,ch7_idx) == -9) then
                grid_swf_climo(ai_idx,sza_idx,ice_idx,ch7_idx) = &
                  CERES_SWF_data(jj,ii)
                grid_swf_count(ai_idx,sza_idx,ice_idx,ch7_idx) = 1

              else
                grid_swf_climo(ai_idx,sza_idx,ice_idx,ch7_idx) = &
                  grid_swf_climo(ai_idx,sza_idx,ice_idx,ch7_idx) + &
                  CERES_SWF_data(jj,ii)
                grid_swf_count(ai_idx,sza_idx,ice_idx,ch7_idx) = &
                  grid_swf_count(ai_idx,sza_idx,ice_idx,ch7_idx) + 1

              endif
            endif 
          enddo omi_loop2

        enddo omi_loop1
   
        ! Deallocate all the arrays for the next pass
        ! -------------------------------------------
        call clear_arrays

        ! Close file
        ! ----------
        call h5fclose_f(file_id, error)
        if(error /= 0) then
          write(*,*) 'FATAL ERROR: could not close file'
          return
        endif

      endif
    enddo file_loop  
  endif

  num_data_grids = 0 
  do ii = 1, len_ch7
    do jj = 1, len_ice
      do kk = 1, len_sza
        do nn = 1, len_ai
          if(grid_swf_count(nn,kk,jj,ii) >= 0) then
            !write(*,*) 'before', grid_swf_climo(nn,kk,jj,ii), grid_swf_count(nn,kk,jj,ii)
            grid_swf_climo(nn,kk,jj,ii) = &
              grid_swf_climo(nn,kk,jj,ii) / grid_swf_count(nn,kk,jj,ii)
            !write(*,*) 'after', grid_swf_climo(nn,kk,jj,ii)
            num_data_grids = num_data_grids + 1
          endif
        enddo
      enddo
    enddo
  enddo


  ! Open the output file
  ! --------------------
  out_file_name = 'comp_grid_climo_v1.hdf5'
  write(*,*) trim(out_file_name)

  call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open output file'
    return
  endif

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write OMI AI bins
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, ai_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_ai_bins', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_ai_bins'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ai_bins, ai_dims, &
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

  write(*,*) 'Wrote OMI ai bins'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write OMI SZA bins
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, sza_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_sza_bins', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_sza_bins'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, sza_bins, sza_dims, &
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

  write(*,*) 'Wrote OMI sza bins'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write NSIDC ICE bins
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, ice_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'nsidc_ice_bins', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_ice_bins'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ice_bins, ice_dims, &
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

  write(*,*) 'Wrote NSIDC ice bins'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write MODIS CH7 bins
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, ch7_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'modis_ch7_bins', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'modis_ch7_bins'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ch7_bins, ch7_dims, &
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

  write(*,*) 'Wrote modis ch7 bins'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write OMI AI edges
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, ai_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_ai_edges', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_ai_edges'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ai_edges, ai_dims, &
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

  write(*,*) 'Wrote OMI ai edges'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write OMI SZA edges
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, sza_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'omi_sza_edges', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'omi_sza_edges'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, sza_edges, sza_dims, &
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

  write(*,*) 'Wrote OMI sza edges'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write NSIDC ICE edges
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, ice_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'nsidc_ice_edges', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_ice_edges'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ice_edges, ice_dims, &
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

  write(*,*) 'Wrote NSIDC ice edges'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write MODIS CH7 edges
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, ch7_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'modis_ch7_edges', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'modis_ch7_edges'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ch7_edges, ch7_dims, &
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

  write(*,*) 'Wrote modis ch7 edges'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write gridded CERES values
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(4, grid_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'ceres_swf_climo', H5T_NATIVE_REAL, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_swf_climo'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, grid_swf_climo, grid_dims, &
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

  write(*,*) 'Wrote ceres_swf_climo'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write gridded CERES count
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(4, grid_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'ceres_swf_count', H5T_NATIVE_INTEGER, dspace_id_OLT,  &
                   dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_swf_count'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_INTEGER, grid_swf_count, grid_dims, &
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

  write(*,*) 'Wrote ceres_swf_count'







  ! Close output file
  ! -----------------
  call h5fclose_f(out_file_id, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not close output file'
    return
  endif


  !!#!write(*,*) 'pixels containing data:',num_data_grids

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write data to output file
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Write the following to the file:
  !!#!! - gridded SWF climatology
  !!#!! - counts in each bin of SWF climatology
  !!#!! - climo bins of AI, ICE, SZA, and CH7
  !!#!! ---------------------------------------

  !!#!! Open the output file
  !!#!! --------------------
  !!#!out_file_name = 'comp_grid_climo_v1.hdf5'
  !!#!write(*,*) trim(out_file_name)

  !!#!call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open output file'
  !!#!  return
  !!#!endif

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = 
  !!#!!
  !!#!! Write OMI AI bins
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = 

  !!#!rank = 1

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
  !!#!call h5dwrite_f(dset_id_OLT, H5T_NATIVE_DOUBLE, OMI_LAT_data, OMI_AI_dims, &
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


  !!#!! Close output file
  !!#!! -----------------
  !!#!call h5fclose_f(out_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close output file'
  !!#!  return
  !!#!endif

  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !
  ! Deallocate all used memory and clean up
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  close(io8)

  ! Clear out the allocated climo stuff
  ! -----------------------------------
  deallocate(grid_swf_climo)
  deallocate(grid_swf_count)

  ! Clear out the allocated bins
  ! ---------------------------- 
  deallocate(ai_bins)
  deallocate(sza_bins)
  deallocate(ice_bins)
  deallocate(ch7_bins)

  deallocate(ai_edges)
  deallocate(sza_edges)
  deallocate(ice_edges)
  deallocate(ch7_edges)

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

  !!#!! Check command line arguments
  !!#!! ----------------------------
  !!#!arg_count = command_argument_count()
  !!#!if(arg_count /= 5) then
  !!#!  write(*,*) 'SYNTAX: ./omi_exec modis_name1 modis_name2 nsidc_name ' &
  !!#!      //'ceres_name omi_name'
  !!#!  return
  !!#!endif

  !!#!! Pull the output file name and input file list file from the command line
  !!#!! ------------------------------------------------------------------------
  !!#!call get_command_argument(1,modis_name1)
  !!#!call get_command_argument(2,modis_name2)
  !!#!call get_command_argument(3,nsidc_name)
  !!#!call get_command_argument(4,ceres_name)
  !!#!call get_command_argument(5,omi_name)

  !!#!! Using the OMI name, determine if a matching TROPOMI colocation file exists
  !!#!! --------------------------------------------------------------------------
  !!#!l_trop_found = .false.
  !!#!trop_name = '/home/bsorenson/OMI/tropgen_comp_grid_climo/coloc_data/'//&
  !!#!  'colocated_tropomi_'//omi_name(len(trim(omi_name)) - &
  !!#!  16:len(trim(omi_name)) - 5)//'.hdf5'
  !!#!
  !!#!write(*,*) trim(trop_name)
  !!#!  
  !!#!inquire(FILE=trop_name, EXIST=l_trop_found)
  !!#!write(*,*) 'TROP file exist? :',l_trop_found 

  !!#!! Initialize the HDF5 interface
  !!#!! -----------------------------
  !!#!call h5open_f(error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: Could not open HDF5 library'
  !!#!  return
  !!#!endif
  !!#!write(*,*) "Interface opened"

  !!#!! Open the MODIS Ch2 HDF5 file
  !!#!! ----------------------------
  !!#!call h5fopen_f(trim(modis_name1), H5F_ACC_RDWR_F, modis1_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open MODIS CH2 file'
  !!#!  return
  !!#!endif

  !!#!! Open the MODIS Ch7 HDF5 file
  !!#!! ----------------------------
  !!#!call h5fopen_f(trim(modis_name2), H5F_ACC_RDWR_F, modis2_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open MODIS CH7 file'
  !!#!  return
  !!#!endif

  !!#!! Open the NSIDC file
  !!#!! -------------------
  !!#!call h5fopen_f(trim(nsidc_name), H5F_ACC_RDWR_F, nsidc_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open NSIDC file'
  !!#!  return
  !!#!endif

  !!#!! Open the CERES file
  !!#!! -------------------
  !!#!call h5fopen_f(trim(ceres_name), H5F_ACC_RDWR_F, ceres_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open CERES file'
  !!#!  return
  !!#!endif

  !!#!! Open the OMI file
  !!#!! -----------------
  !!#!call h5fopen_f(trim(omi_name), H5F_ACC_RDWR_F, omi_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open OMI file'
  !!#!  return
  !!#!endif

  !!#!! If TROPOMI file exists, open and read from the file
  !!#!! ---------------------------------------------------
  !!#!if(l_trop_found) then

  !!#!  ! Open the TROP file
  !!#!  ! -----------------
  !!#!  call h5fopen_f(trim(trop_name), H5F_ACC_RDWR_F, trop_file_id, error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open OMI file'
  !!#!    return
  !!#!  endif
 

  !!#!  write(*,*) "Opened TROPOMI file"

  !!#!  call read_comp_TROP_AI(trop_file_id)
  !!#!  call read_comp_TROP_SSA0(trop_file_id)
  !!#!  call read_comp_TROP_SSA1(trop_file_id)
  !!#!  call read_comp_TROP_SSA2(trop_file_id)

  !!#!endif

  !!#!call read_comp_MODIS_CH2(modis1_file_id)
  !!#!call read_comp_MODIS_CH7(modis2_file_id)
  !!#!call read_comp_MODIS_LAT(modis2_file_id)
  !!#!call read_comp_MODIS_LON(modis2_file_id)
  !!#!call read_comp_NSIDC_ICE(nsidc_file_id)
  !!#!call read_comp_NSIDC_LAT(nsidc_file_id)
  !!#!call read_comp_NSIDC_LON(nsidc_file_id)
  !!#!call read_comp_CERES_SWF(ceres_file_id)
  !!#!call read_comp_CERES_LWF(ceres_file_id)
  !!#!call read_comp_CERES_LAT(ceres_file_id)
  !!#!call read_comp_CERES_LON(ceres_file_id)
  !!#!call read_comp_OMI_AI(omi_file_id)
  !!#!call read_comp_OMI_AI_raw(omi_file_id)
  !!#!call read_comp_OMI_LAT(omi_file_id)
  !!#!call read_comp_OMI_LON(omi_file_id)
  !!#!call read_comp_OMI_LATCRNR(omi_file_id)
  !!#!call read_comp_OMI_LONCRNR(omi_file_id)
  !!#!call read_comp_OMI_SZA(omi_file_id)
  !!#!call read_comp_OMI_VZA(omi_file_id)
  !!#!call read_comp_OMI_AZM(omi_file_id)

  !!#!!test_dims = (/10, 20/)
  !!#!test_dims = (/OMI_AI_dims(1), OMI_AI_dims(2)/)
  !!#!!test_dims = OMI_AI_dims

  !!#!! Allocate the output arrays
  !!#!! --------------------------
  !!#!call allocate_out_arrays(OMI_AI_dims(1), OMI_AI_dims(2))

  !!#!!!#!write(*,*) distance

  !!#!istatus = 0
  !!#!num_nan = 0 
  !!#!min_dist  = 50.
  !!#!closest_dist = 999999.
  !!#!match_lat  = -999.
  !!#!match_lon  = -999.
  !!#!match_data1 = -999.
  !!#!match_data2 = -999.

  !!#!omi_loop1: do ii=1,OMI_LAT_dims(2)

  !!#!  if(mod(ii, 50) == 0) then
  !!#!    write(*,*) ii
  !!#!  endif

  !!#!  omi_loop2: do jj=1,OMI_LAT_dims(1) 

  !!#!    ! If TROPOMI data exists, put value in local variable. If not, set
  !!#!    ! local variable equal to -999.
  !!#!    if(l_trop_found) then
  !!#!      local_trop_ai = TROP_AI_data(jj,ii)
  !!#!    else
  !!#!      local_trop_ai = -999.
  !!#!    endif
  !!#!      
  !!#!    !!#!  (l_trop_found .and. ((.not.isnan(local_trop_ai) .and. &
  !!#!    !!#!  (local_trop_ai /=-999.)) .or. &
  !!#!    !!#!    .not.isnan(OMI_AI_data(jj,ii)))),  &
  !!#!    !!#!    (.not.l_trop_found .and. .not.isnan(OMI_AI_data(jj,ii)))
  !!#!    ! Check if the current pixel is missing
  !!#!    ! -------------------------------------
  !!#!    if( (l_trop_found .and. ((.not.isnan(local_trop_ai) .and. &
  !!#!         (local_trop_ai/= -999.)) .or. &
  !!#!        .not.isnan(OMI_AI_raw_data(jj,ii)))) .or. &
  !!#!        (.not.l_trop_found .and. .not.isnan(OMI_AI_raw_data(jj,ii)))) then
  !!#!      ! Now, loop over the NSIDC data
  !!#!      ! -----------------------------
  !!#!      nsidc_loop1: do nii=1,NSIDC_LAT_dims(2)
  !!#!        nsidc_loop2: do njj=1,NSIDC_LAT_dims(1) 

  !!#!          ! Calculate the distance between the current pixels
  !!#!          ! -------------------------------------------------
  !!#!          distance = find_distance_between_points(OMI_LAT_data(jj,ii), &
  !!#!                                                  OMI_LON_data(jj,ii), &
  !!#!                                                  NSIDC_LAT_data(njj,nii), &
  !!#!                                                  NSIDC_LON_data(njj,nii))

  !!#!          if(distance < closest_dist) then
  !!#!            closest_dist = distance
  !!#!            !NSIDC_out_LAT_data(jj,ii) = NSIDC_LAT_data(njj, nii) 
  !!#!            !NSIDC_out_LON_data(jj,ii) = NSIDC_LON_data(njj, nii) 
  !!#!            NSIDC_out_data(jj,ii)     = NSIDC_data(njj, nii) 
  !!#!          endif

  !!#!        enddo nsidc_loop2
  !!#!      enddo nsidc_loop1

  !!#!      ! Check the distance requirement for the NSIDC pixel
  !!#!      ! --------------------------------------------------
  !!#!      if(closest_dist > min_dist) then
  !!#!        !NSIDC_out_LAT_data(jj,ii) = -999.
  !!#!        !NSIDC_out_LON_data(jj,ii) = -999.
  !!#!        NSIDC_out_data(jj,ii)     = -999.
  !!#!      endif

  !!#!      closest_dist = 999999.

  !!#!      ! Now, loop over the CERES data
  !!#!      ! -----------------------------
  !!#!      ceres_loop1: do nii=1,CERES_LAT_dims(2)
  !!#!        ceres_loop2: do njj=1,CERES_LAT_dims(1) 

  !!#!          ! Calculate the distance between the current pixels
  !!#!          ! -------------------------------------------------
  !!#!          distance = find_distance_between_points(OMI_LAT_data(jj,ii), &
  !!#!                                                  OMI_LON_data(jj,ii), &
  !!#!                                                  CERES_LAT_data(njj,nii), &
  !!#!                                                  CERES_LON_data(njj,nii))

  !!#!          if(distance < closest_dist) then
  !!#!            closest_dist = distance
  !!#!            !CERES_out_LAT_data(jj,ii) = CERES_LAT_data(njj, nii) 
  !!#!            !CERES_out_LON_data(jj,ii) = CERES_LON_data(njj, nii) 
  !!#!            CERES_out_LWF_data(jj,ii) = CERES_LWF_data(njj, nii) 
  !!#!            CERES_out_SWF_data(jj,ii) = CERES_SWF_data(njj, nii) 
  !!#!          endif

  !!#!        enddo ceres_loop2
  !!#!      enddo ceres_loop1

  !!#!      ! Check the distance requirement for the CERES pixel
  !!#!      ! --------------------------------------------------
  !!#!      if(closest_dist > min_dist) then
  !!#!        !CERES_out_LAT_data(jj,ii) = -999.
  !!#!        !CERES_out_LON_data(jj,ii) = -999.
  !!#!        CERES_out_LWF_data(jj,ii) = -999.
  !!#!        CERES_out_SWF_data(jj,ii) = -999.
  !!#!      endif

  !!#!      closest_dist = 999999.

  !!#!      count_modis_total = 0
  !!#!      run_modis_total_ch2 = 0.0
  !!#!      run_modis_total_ch7 = 0.0
  !!#!      
  !!#!      modis_loop1: do nii=1,MODIS_LAT_dims(2)
  !!#!        modis_loop2: do njj=1,MODIS_LAT_dims(1) 

  !!#!          ! Check if the current pixel is within the OMI pixel bounds
  !!#!          ! ---------------------------------------------------------
  !!#!          if(pixel_in_box(OMI_LATCRNR_data(:,jj,ii), OMI_LONCRNR_data(:,jj,ii), &
  !!#!              MODIS_LAT_data(njj,nii), MODIS_LON_data(njj,nii))) then

  !!#!              run_modis_total_ch2 = &
  !!#!                  run_modis_total_ch2 + MODIS_CH2_data(njj,nii)
  !!#!              run_modis_total_ch7 = &
  !!#!                  run_modis_total_ch7 + MODIS_CH7_data(njj,nii)
  !!#!              count_modis_total = count_modis_total + 1

  !!#!          endif
  !!#!          !!#!! Calculate the distance between the current pixels
  !!#!          !!#!! -------------------------------------------------
  !!#!          !!#!distance = find_distance_between_points(OMI_LAT_data(jj,ii), &
  !!#!          !!#!                                        OMI_LON_data(jj,ii), &
  !!#!          !!#!                                        MODIS_LAT_data(njj,nii), &
  !!#!          !!#!                                        MODIS_LON_data(njj,nii))

  !!#!          !!#!if(distance < closest_dist) then
  !!#!          !!#!  closest_dist = distance
  !!#!          !!#!  !MODIS_out_LAT_data(jj,ii) = MODIS_LAT_data(njj, nii) 
  !!#!          !!#!  !MODIS_out_LON_data(jj,ii) = MODIS_LON_data(njj, nii) 
  !!#!          !!#!  MODIS_out_CH2_data(jj,ii) = MODIS_CH2_data(njj, nii) 
  !!#!          !!#!  MODIS_out_CH7_data(jj,ii) = MODIS_CH7_data(njj, nii) 
  !!#!          !!#!endif

  !!#!        enddo modis_loop2
  !!#!      enddo modis_loop1

  !!#!      !!#!! Check the distance requirement for the MODIS pixel
  !!#!      !!#!! --------------------------------------------------
  !!#!      !!#!if(closest_dist > min_dist) then
  !!#!      !!#!  !MODIS_out_LAT_data(jj,ii) = -999.
  !!#!      !!#!  !MODIS_out_LON_data(jj,ii) = -999.
  !!#!      !!#!  MODIS_out_CH2_data(jj,ii) = -999.
  !!#!      !!#!  MODIS_out_CH7_data(jj,ii) = -999.
  !!#!      !!#!endif
  !!#!      if(count_modis_total == 0) then
  !!#!        !MODIS_out_LAT_data(jj,ii) = -999.
  !!#!        !MODIS_out_LON_data(jj,ii) = -999.
  !!#!        MODIS_out_CH2_data(jj,ii) = -999.
  !!#!        MODIS_out_CH7_data(jj,ii) = -999.
  !!#!      else
  !!#!        MODIS_out_CH2_data(jj,ii) = run_modis_total_ch2 / count_modis_total
  !!#!        MODIS_out_CH7_data(jj,ii) = run_modis_total_ch7 / count_modis_total
  !!#!      endif

  !!#!      closest_dist = 999999.

  !!#!    else

  !!#!      !if(isnan(OMI_AI_data(jj,ii))) then
  !!#!      !write(*,*) 'NAN VALUE'
  !!#!      !NSIDC_out_LAT_data(jj,ii) = -999.
  !!#!      !NSIDC_out_LON_data(jj,ii) = -999.
  !!#!      NSIDC_out_data(jj,ii)     = -999.
  !!#!      !CERES_out_LAT_data(jj,ii) = -999.
  !!#!      !CERES_out_LON_data(jj,ii) = -999.
  !!#!      CERES_out_LWF_data(jj,ii) = -999.
  !!#!      CERES_out_SWF_data(jj,ii) = -999.
  !!#!      !MODIS_out_LAT_data(jj,ii) = -999.
  !!#!      !MODIS_out_LON_data(jj,ii) = -999.
  !!#!      MODIS_out_CH2_data(jj,ii) = -999.
  !!#!      MODIS_out_CH7_data(jj,ii) = -999.
  !!#!      num_nan = num_nan + 1

  !!#!    endif

  !!#!  enddo omi_loop2
  !!#!enddo omi_loop1

  !!#!write(*,*) istatus, num_nan
  !!#!
  !!#!! Open the output file
  !!#!! --------------------
  !!#!write(*,*) omi_name(len(trim(omi_name)) - 16:len(trim(omi_name)) - 5)
  !!#!out_file_name = 'colocated_subset_'//omi_name(len(trim(omi_name)) - &
  !!#!  16:len(trim(omi_name)) - 5)//'.hdf5'
  !!#!write(*,*) trim(out_file_name)

  !!#!call h5fcreate_f(trim(out_file_name), H5F_ACC_TRUNC_F, out_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open output file'
  !!#!  return
  !!#!endif

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
  !!#!call h5dwrite_f(dset_id_OLT, H5T_NATIVE_DOUBLE, OMI_LAT_data, OMI_AI_dims, &
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
  !!#!call h5dwrite_f(dset_id_OLN, H5T_NATIVE_DOUBLE, OMI_LON_data, OMI_AI_dims, &
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
  !!#!call h5dwrite_f(dset_id_OSZ, H5T_NATIVE_DOUBLE, OMI_SZA_data, OMI_AI_dims, &
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
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'omi_vza'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OVZ, H5T_NATIVE_DOUBLE, OMI_VZA_data, OMI_AI_dims, &
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
  !!#!call h5dwrite_f(dset_id_OAZ, H5T_NATIVE_DOUBLE, OMI_AZM_data, OMI_AI_dims, &
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
  !!#!! Write OMI AI pert data
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
  !!#!call h5dwrite_f(dset_id_OAI, H5T_NATIVE_DOUBLE, OMI_AI_data, OMI_AI_dims, &
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
  !!#!! Write OMI AI raw data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_OAR, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'omi_uvai_raw', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_OAR,  dset_id_OAR, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'omi_uvai_raw'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OAR, H5T_NATIVE_DOUBLE, OMI_AI_raw_data, OMI_AI_raw_dims, &
  !!#!                    error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_OAR, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_OAR, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote OMI AI raw'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write CERES LAT data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_CLT, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'ceres_lat', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_CLT,  dset_id_CLT, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_lat'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_CLT, H5T_NATIVE_DOUBLE, CERES_out_LAT_data, &
  !!#!!!#!                OMI_AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_CLT, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_CLT, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote CERES LAT'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write CERES LON data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_CLN, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'ceres_lon', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_CLN,  dset_id_CLN, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_lon'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_CLN, H5T_NATIVE_DOUBLE, CERES_out_LON_data, &
  !!#!!!#!                OMI_AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_CLN, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_CLN, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote CERES LON'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write CERES LWF data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_CLW, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'ceres_lwf', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_CLW,  dset_id_CLW, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_lwf'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_CLW, H5T_NATIVE_DOUBLE, CERES_out_LWF_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_CLW, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_CLW, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote CERES LWF'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write CERES SWF data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_CSW, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'ceres_swf', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_CSW,  dset_id_CSW, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_swf'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_CSW, H5T_NATIVE_DOUBLE, CERES_out_SWF_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_CSW, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_CSW, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote CERES SWF'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write MODIS LAT data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_MLT, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'modis_lat', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_MLT,  dset_id_MLT, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_lat'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_MLT, H5T_NATIVE_DOUBLE, MODIS_out_LAT_data, &
  !!#!!!#!                OMI_AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_MLT, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_MLT, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote MODIS LAT'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write MODIS LON data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_MLN, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'modis_lon', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_MLN,  dset_id_MLN, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_lon'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_MLN, H5T_NATIVE_DOUBLE, MODIS_out_LON_data, &
  !!#!!!#!                OMI_AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_MLN, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_MLN, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote MODIS LON'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write MODIS CH2 data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_MC2, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'modis_ch2', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_MC2,  dset_id_MC2, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_ch2'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_MC2, H5T_NATIVE_DOUBLE, MODIS_out_CH2_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_MC2, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_MC2, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote MODIS CH2'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write MODIS CH7 data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_MC7, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'modis_ch7', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_MC7,  dset_id_MC7, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_ch7'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_MC7, H5T_NATIVE_DOUBLE, MODIS_out_CH7_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_MC7, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_MC7, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote MODIS CH7'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write NSIDC LAT data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_NLT, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'nsidc_lat', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_NLT,  dset_id_NLT, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_lat'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_NLT, H5T_NATIVE_DOUBLE, NSIDC_out_LAT_data, &
  !!#!!!#!                OMI_AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_NLT, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_NLT, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote NSIDC LAT'

  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!!#!!
  !!#!!!#!! Write NSIDC LON data
  !!#!!!#!!
  !!#!!!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!!!#!! Create the dataspace
  !!#!!!#!! --------------------
  !!#!!!#!call h5screate_simple_f(rank, test_dims, dspace_id_NLN, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Create the dataset
  !!#!!!#!! ------------------
  !!#!!!#!call h5dcreate_f(out_file_id, 'nsidc_lon', H5T_NATIVE_DOUBLE, &
  !!#!!!#!                 dspace_id_NLN,  dset_id_NLN, error) 
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_lon'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Write to the dataset
  !!#!!!#!! --------------------
  !!#!!!#!call h5dwrite_f(dset_id_NLN, H5T_NATIVE_DOUBLE, NSIDC_out_LON_data, &
  !!#!!!#!                OMI_AI_dims, error)
  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!! Close the dataset
  !!#!!!#!! -----------------
  !!#!!!#!call h5dclose_f(dset_id_NLN, error)

  !!#!!!#!! Close access to data space rank
  !!#!!!#!call h5sclose_f(dspace_id_NLN, error)

  !!#!!!#!if(error /= 0) then
  !!#!!!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!!!#!  return
  !!#!!!#!endif

  !!#!!!#!write(*,*) 'Wrote NSIDC LON'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! Write NSIDC NIC data
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(rank, test_dims, dspace_id_NIC, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'nsidc_ice', H5T_NATIVE_DOUBLE, &
  !!#!                 dspace_id_NIC,  dset_id_NIC, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'nsidc_ice'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_NIC, H5T_NATIVE_DOUBLE, NSIDC_out_data, &
  !!#!                OMI_AI_dims, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!! Close the dataset
  !!#!! -----------------
  !!#!call h5dclose_f(dset_id_NIC, error)

  !!#!! Close access to data space rank
  !!#!call h5sclose_f(dspace_id_NIC, error)

  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!  return
  !!#!endif

  !!#!write(*,*) 'Wrote NSIDC ice'

  !!#!if(l_trop_found) then
  !!#!  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!  !
  !!#!  ! Write TROPOMI AI data
  !!#!  !
  !!#!  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!  ! Create the dataspace
  !!#!  ! --------------------
  !!#!  call h5screate_simple_f(rank, test_dims, dspace_id_TAI, error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!    return
  !!#!  endif

  !!#!  ! Create the dataset
  !!#!  ! ------------------
  !!#!  call h5dcreate_f(out_file_id, 'trop_uvai', H5T_NATIVE_DOUBLE, &
  !!#!                   dspace_id_TAI,  dset_id_TAI, error) 
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_uvai'
  !!#!    return
  !!#!  endif

  !!#!  ! Write to the dataset
  !!#!  ! --------------------
  !!#!  call h5dwrite_f(dset_id_TAI, H5T_NATIVE_DOUBLE, TROP_AI_data, TROP_AI_dims, &
  !!#!                      error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!    return
  !!#!  endif

  !!#!  ! Close the dataset
  !!#!  ! -----------------
  !!#!  call h5dclose_f(dset_id_TAI, error)

  !!#!  ! Close access to data space rank
  !!#!  call h5sclose_f(dspace_id_TAI, error)

  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!    return
  !!#!  endif

  !!#!  write(*,*) 'Wrote TROPOMI AI'

  !!#!  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!  !
  !!#!  ! Write TROPOMI SSA0 data
  !!#!  !
  !!#!  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!  ! Create the dataspace
  !!#!  ! --------------------
  !!#!  call h5screate_simple_f(rank, test_dims, dspace_id_TS0, error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!    return
  !!#!  endif

  !!#!  ! Create the dataset
  !!#!  ! ------------------
  !!#!  call h5dcreate_f(out_file_id, 'trop_ssa0', H5T_NATIVE_DOUBLE, &
  !!#!                   dspace_id_TS0,  dset_id_TS0, error) 
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa0'
  !!#!    return
  !!#!  endif

  !!#!  ! Write to the dataset
  !!#!  ! --------------------
  !!#!  call h5dwrite_f(dset_id_TS0, H5T_NATIVE_DOUBLE, TROP_SSA0_data, TROP_SSA0_dims, &
  !!#!                      error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!    return
  !!#!  endif

  !!#!  ! Close the dataset
  !!#!  ! -----------------
  !!#!  call h5dclose_f(dset_id_TS0, error)

  !!#!  ! Close access to data space rank
  !!#!  call h5sclose_f(dspace_id_TS0, error)

  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!    return
  !!#!  endif

  !!#!  write(*,*) 'Wrote TROPOMI SSA0'

  !!#!  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!  !
  !!#!  ! Write TROPOMI SSA1 data
  !!#!  !
  !!#!  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!  ! Create the dataspace
  !!#!  ! --------------------
  !!#!  call h5screate_simple_f(rank, test_dims, dspace_id_TS1, error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!    return
  !!#!  endif

  !!#!  ! Create the dataset
  !!#!  ! ------------------
  !!#!  call h5dcreate_f(out_file_id, 'trop_ssa1', H5T_NATIVE_DOUBLE, &
  !!#!                   dspace_id_TS1,  dset_id_TS1, error) 
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa1'
  !!#!    return
  !!#!  endif

  !!#!  ! Write to the dataset
  !!#!  ! --------------------
  !!#!  call h5dwrite_f(dset_id_TS1, H5T_NATIVE_DOUBLE, TROP_SSA1_data, TROP_SSA1_dims, &
  !!#!                      error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!    return
  !!#!  endif

  !!#!  ! Close the dataset
  !!#!  ! -----------------
  !!#!  call h5dclose_f(dset_id_TS1, error)

  !!#!  ! Close access to data space rank
  !!#!  call h5sclose_f(dspace_id_TS1, error)

  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!    return
  !!#!  endif

  !!#!  write(*,*) 'Wrote TROPOMI SSA1'

  !!#!  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!  !
  !!#!  ! Write TROPOMI SSA2 data
  !!#!  !
  !!#!  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!  ! Create the dataspace
  !!#!  ! --------------------
  !!#!  call h5screate_simple_f(rank, test_dims, dspace_id_TS2, error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!    return
  !!#!  endif

  !!#!  ! Create the dataset
  !!#!  ! ------------------
  !!#!  call h5dcreate_f(out_file_id, 'trop_ssa2', H5T_NATIVE_DOUBLE, &
  !!#!                   dspace_id_TS2,  dset_id_TS2, error) 
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not open dataset '//'trop_ssa2'
  !!#!    return
  !!#!  endif

  !!#!  ! Write to the dataset
  !!#!  ! --------------------
  !!#!  call h5dwrite_f(dset_id_TS2, H5T_NATIVE_DOUBLE, TROP_SSA2_data, TROP_SSA2_dims, &
  !!#!                      error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!    return
  !!#!  endif

  !!#!  ! Close the dataset
  !!#!  ! -----------------
  !!#!  call h5dclose_f(dset_id_TS2, error)

  !!#!  ! Close access to data space rank
  !!#!  call h5sclose_f(dspace_id_TS2, error)

  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not write to dataset'
  !!#!    return
  !!#!  endif

  !!#!  write(*,*) 'Wrote TROPOMI SSA2'

  !!#!endif

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  !!#!!
  !!#!! End Variable Writing
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  !!#!! Close output file
  !!#!! -----------------
  !!#!call h5fclose_f(out_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close output file'
  !!#!  return
  !!#!endif

  !!#!!!#!write(*,*) 'Saved output file'//trim(out_file_name)
 
  !!#!! Deallocate all the arrays for the next pass
  !!#!! -------------------------------------------
  !!#!call clear_arrays

  !!#!! Close file
  !!#!! ----------
  !!#!call h5fclose_f(modis1_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close MODIS CH2 file'
  !!#!  return
  !!#!endif

  !!#!call h5fclose_f(modis2_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close MODIS CH7 file'
  !!#!  return
  !!#!endif

  !!#!call h5fclose_f(nsidc_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close NSIDC file'
  !!#!  return
  !!#!endif

  !!#!call h5fclose_f(ceres_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close CERES file'
  !!#!  return
  !!#!endif

  !!#!call h5fclose_f(omi_file_id, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not close OMI file'
  !!#!  return
  !!#!endif

  !!#!if(l_trop_found) then
  !!#!  call h5fclose_f(trop_file_id, error)
  !!#!  if(error /= 0) then
  !!#!    write(*,*) 'FATAL ERROR: could not close TROPOMI file'
  !!#!    return
  !!#!  endif
  !!#!endif

  !!#!! Close the HDF5 interface
  !!#!! ------------------------
  !!#!call h5close_f(error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: Could not close HDF5 library'
  !!#!  return
  !!#!endif
  !!#!write(*,*) "Interface closed"

end program gen_comp_grid_climo
