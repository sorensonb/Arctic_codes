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
    MODIS_CLD_data, MODIS_CLD_dims, &
    !MODIS_LAT_data, MODIS_LAT_dims, &
    !MODIS_LON_data, MODIS_LON_dims, &
    NSIDC_data,     NSIDC_dims, &
    !NSIDC_LAT_data, NSIDC_LAT_dims, &
    !NSIDC_LON_data, NSIDC_LON_dims, &
    !CERES_CLD_data, CERES_CLD_dims, &
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
    MODIS_out_CLD_data, &
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

  integer                                       :: ii            ! loop counter
  integer                                       :: jj            ! loop counter
  integer                                       :: kk            ! loop counter
  integer                                       :: nn            ! loop counter
  integer                                       :: mm            ! loop counter
  integer                                       :: arg_count     ! Number of arguments
  integer                                       :: error         ! error flag
  integer                                       :: istatus
  integer                                       :: file_id
  integer                                       :: out_file_id   ! id for current H5 file
  integer                                       :: dspace_id_OLT ! OMI LAT
  integer                                       :: dset_id_OLT   ! OMI LAT

  integer                                       :: len_ai
  integer                                       :: len_sza
  integer                                       :: len_ice
  integer                                       :: len_ch7
  !integer                                       :: len_cld
  integer                                       :: num_data_grids

  integer                                       :: ai_idx
  integer                                       :: sza_idx
  integer                                       :: ice_idx
  integer                                       :: ch7_idx
  !integer                                       :: cld_idx

  integer                                       :: io8
  integer                                       :: io5

  real, dimension(:), allocatable                  :: ai_bins
  real, dimension(:), allocatable                  :: sza_bins
  real, dimension(:), allocatable                  :: ice_bins
  real, dimension(:), allocatable                  :: ch7_bins
  !real, dimension(:), allocatable                  :: cld_bins

  real, dimension(:), allocatable                  :: ai_edges
  real, dimension(:), allocatable                  :: sza_edges
  real, dimension(:), allocatable                  :: ice_edges
  real, dimension(:), allocatable                  :: ch7_edges
  !real, dimension(:), allocatable                  :: cld_edges

  integer, dimension(:,:,:,:), allocatable       :: grid_swf_count
  real, dimension(:,:,:,:), allocatable          :: grid_swf_climo

  real                                             :: min_ai
  real                                             :: max_ai
  real                                             :: delta_ai
  real                                             :: min_sza
  real                                             :: max_sza
  real                                             :: delta_sza
  real                                             :: min_ice
  real                                             :: max_ice
  real                                             :: delta_ice
  real                                             :: min_ch7
  real                                             :: max_ch7
  real                                             :: delta_ch7
  !real                                             :: min_cld
  !real                                             :: max_cld
  !real                                             :: delta_cld
  real                                             :: min_lat

  real                                             :: found_max_ai

  integer(hsize_t), dimension(4)                   :: grid_dims
  integer(hsize_t), dimension(1)                   :: ai_dims
  integer(hsize_t), dimension(1)                   :: sza_dims
  integer(hsize_t), dimension(1)                   :: ice_dims
  integer(hsize_t), dimension(1)                   :: ch7_dims
  !integer(hsize_t), dimension(1)                   :: cld_dims

  integer(hsize_t), dimension(1)                   :: ai_edge_dims
  integer(hsize_t), dimension(1)                   :: sza_edge_dims
  integer(hsize_t), dimension(1)                   :: ice_edge_dims
  integer(hsize_t), dimension(1)                   :: ch7_edge_dims
  !integer(hsize_t), dimension(1)                   :: cld_edge_dims

  character(len = 255)                             :: file_name_file 
  character(len = 255)                             :: total_file_name
  character(len = 12)                              :: single_file_date
  character(len = 255)                             :: out_file_name   

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

  !min_ch7   = 0.0
  !max_ch7   = 0.350
  !delta_ch7 = 0.025

  min_ch7   = 0.0
  max_ch7   = 3.0
  delta_ch7 = 1.

  !min_cld   = 0.0
  !max_cld   = 3.0
  !delta_cld = 1.

  min_lat = 70.

  ! Determine the number of bins for each variable
  ! ----------------------------------------------
  len_ai  = (max_ai  - min_ai  + delta_ai)  / delta_ai
  len_sza = (max_sza - min_sza + delta_sza) / delta_sza
  len_ice = (max_ice - min_ice + delta_ice) / delta_ice
  len_ch7 = (max_ch7 - min_ch7 + delta_ch7) / delta_ch7
  !len_cld = (max_cld - min_cld + delta_cld) / delta_cld

  allocate(ai_bins(len_ai))
  allocate(sza_bins(len_sza))
  allocate(ice_bins(len_ice))
  allocate(ch7_bins(len_ch7))
  !allocate(cld_bins(len_cld))

  allocate(ai_edges(len_ai + 1))
  allocate(sza_edges(len_sza + 1))
  allocate(ice_edges(len_ice + 1))
  allocate(ch7_edges(len_ch7 + 1))
  !allocate(cld_edges(len_cld + 1))
 
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

  !do ii = 1, len_cld
  !  cld_bins(ii) = min_cld + delta_cld*(ii-1)
  !enddo

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

  !do ii = 1, len_cld + 1 
  !  cld_edges(ii) = min_cld + delta_cld*((ii-1) - 0.5)
  !enddo

  write(*,*) 'total_size = ',len_ai * len_sza * len_ice * len_ch7

  ai_dims   = (/size(ai_bins)/)
  sza_dims  = (/size(sza_bins)/)
  ice_dims  = (/size(ice_bins)/)
  ch7_dims  = (/size(ch7_bins)/)
  !cld_dims  = (/size(cld_bins)/)
  grid_dims = (/size(ai_bins),size(sza_bins),size(ice_bins), &
    size(ch7_bins)/)

  ai_edge_dims  = (/size(ai_edges)/)
  sza_edge_dims = (/size(sza_edges)/)
  ice_edge_dims = (/size(ice_edges)/)
  ch7_edge_dims = (/size(ch7_edges)/)
  !cld_edge_dims = (/size(cld_edges)/)

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

  found_max_ai = -100.

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
        ! As of 2023/06/05, these are the available variables from the 
        ! colocated_subset files:
        !  ceres_cld
        !  ceres_lwf
        !  x ceres_swf
        !  modis_ch2
        !  x modis_ch7
        !  x modis_cld
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
        call read_coloc_MODIS_CLD(file_id)
        call read_coloc_CERES_SWF(file_id)
        !call read_coloc_CERES_CLD(file_id)
        call read_coloc_NSIDC_ICE(file_id)
        call read_coloc_OMI_AI_raw(file_id)
        call read_coloc_OMI_LAT(file_id)
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

        
            if(OMI_LAT_data(jj,ii) >= min_lat) then

              ! Insert each pixel into the grid here
              ! ------------------------------------      
              ai_idx  = minloc(abs(ai_bins - OMI_AI_raw_data(jj,ii)), dim = 1)
              sza_idx = minloc(abs(sza_bins - OMI_SZA_data(jj,ii)), dim = 1)
              ice_idx = minloc(abs(ice_bins - NSIDC_data(jj,ii)), dim = 1)
              ! NOTE: For "v7", am substituting the MODIS CH7 values for
              !       the MODIS CLD cloud mask values
              ch7_idx = minloc(abs(ch7_bins - MODIS_CLD_data(jj,ii)), dim = 1)
              !ch7_idx = minloc(abs(ch7_bins - MODIS_CH7_data(jj,ii)), dim = 1)
              !cld_idx = minloc(abs(cld_bins - MODIS_CLD_data(jj,ii)), dim = 1)

              if(((CERES_SWF_data(jj,ii) > 0) .and. &
                  (CERES_SWF_data(jj,ii) < 5000)) .and. &
                 (abs(OMI_AI_raw_data(jj,ii)) < max_ai)) then
                !if(grid_swf_count(ai_idx, sza_idx,ice_idx,ch7_idx,cld_idx) &
                if(grid_swf_count(ai_idx, sza_idx,ice_idx,ch7_idx) &
                    == -9) then
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

                if(OMI_AI_raw_data(jj,ii) > found_max_ai) then
                  found_max_ai = OMI_AI_raw_data(jj,ii)
                endif

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

  write(*,*) 'Found max AI', found_max_ai

  num_data_grids = 0 
  !do ii = 1, len_cld
  do jj = 1, len_ch7
    do kk = 1, len_ice
      do nn = 1, len_sza
        do mm = 1, len_ai
          !if(grid_swf_count(mm,nn,kk,jj,ii) > 0) then
          if(grid_swf_count(mm,nn,kk,jj) > 0) then
            !write(*,*) 'before', grid_swf_climo(mm,nn,kk,jj), &
            !  grid_swf_count(mm,nn,kk,jj)
            grid_swf_climo(mm,nn,kk,jj) = &
              grid_swf_climo(mm,nn,kk,jj) / grid_swf_count(mm,nn,kk,jj)
            !write(*,*) 'after', grid_swf_climo(mm,nn,kk,jj)
            num_data_grids = num_data_grids + 1
          endif
        enddo
      enddo
    enddo
  enddo
  !enddo

  ! Open the output file
  ! --------------------
  !out_file_name = 'comp_grid_climo_v6.hdf5'
  out_file_name = 'comp_grid_climo_v7.hdf5'
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

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = 
  !!#!!
  !!#!! Write MODIS CLD bins
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = 
      
  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(1, cld_dims, dspace_id_OLT, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif
      
  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'modis_cld_bins', H5T_NATIVE_REAL, dspace_id_OLT,  &
  !!#!                 dset_id_OLT, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_cld_bins'
  !!#!  return
  !!#!endif
      
  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, cld_bins, cld_dims, &
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
      
  !!#!write(*,*) 'Wrote modis cld bins'

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = 
  !!#!!
  !!#!! Write CERES CLD bins
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = 

  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(1, cld_dims, dspace_id_OLT, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif

  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'ceres_cld_bins', H5T_NATIVE_REAL, dspace_id_OLT,  &
  !!#!                 dset_id_OLT, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'ceres_cld_bins'
  !!#!  return
  !!#!endif

  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, cld_bins, cld_dims, &
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

  !!#!write(*,*) 'Wrote ceres cld bins'

  ! = = = = = = = = = = = = = = = = = = = = = = = 
  !
  ! Write OMI AI edges
  !
  ! = = = = = = = = = = = = = = = = = = = = = = = 

  ! Create the dataspace
  ! --------------------
  call h5screate_simple_f(1, ai_edge_dims, dspace_id_OLT, error)
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
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ai_edges, ai_edge_dims, &
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
  call h5screate_simple_f(1, sza_edge_dims, dspace_id_OLT, error)
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
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, sza_edges, sza_edge_dims, &
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
  call h5screate_simple_f(1, ice_edge_dims, dspace_id_OLT, error)
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
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ice_edges, ice_edge_dims, &
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
  call h5screate_simple_f(1, ch7_edge_dims, dspace_id_OLT, error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataspace'
    return
  endif

  ! Create the dataset
  ! ------------------
  call h5dcreate_f(out_file_id, 'modis_ch7_edges', H5T_NATIVE_REAL, &
                   dspace_id_OLT, dset_id_OLT, error) 
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: could not open dataset '//'modis_ch7_edges'
    return
  endif

  ! Write to the dataset
  ! --------------------
  call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, ch7_edges, ch7_edge_dims, &
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

  !!#!! = = = = = = = = = = = = = = = = = = = = = = = 
  !!#!!
  !!#!! Write MODIS CLD edges
  !!#!!
  !!#!! = = = = = = = = = = = = = = = = = = = = = = = 
      
  !!#!! Create the dataspace
  !!#!! --------------------
  !!#!call h5screate_simple_f(1, cld_edge_dims, dspace_id_OLT, error)
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataspace'
  !!#!  return
  !!#!endif
      
  !!#!! Create the dataset
  !!#!! ------------------
  !!#!call h5dcreate_f(out_file_id, 'modis_cld_edges', H5T_NATIVE_REAL, dspace_id_OLT,  &
  !!#!                 dset_id_OLT, error) 
  !!#!if(error /= 0) then
  !!#!  write(*,*) 'FATAL ERROR: could not open dataset '//'modis_cld_edges'
  !!#!  return
  !!#!endif
      
  !!#!! Write to the dataset
  !!#!! --------------------
  !!#!call h5dwrite_f(dset_id_OLT, H5T_NATIVE_REAL, cld_edges, cld_edge_dims, &
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
      
  !!#!write(*,*) 'Wrote modis cld edges'

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
  call h5dcreate_f(out_file_id, 'ceres_swf_climo', H5T_NATIVE_REAL, &
                   dspace_id_OLT, dset_id_OLT, error) 
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
  call h5dcreate_f(out_file_id, 'ceres_swf_count', H5T_NATIVE_INTEGER, &
                   dspace_id_OLT, dset_id_OLT, error) 
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
  !deallocate(cld_bins)

  deallocate(ai_edges)
  deallocate(sza_edges)
  deallocate(ice_edges)
  deallocate(ch7_edges)
  !deallocate(cld_edges)

  ! Close the HDF5 interface
  ! ------------------------
  call h5close_f(error)
  if(error /= 0) then
    write(*,*) 'FATAL ERROR: Could not close HDF5 library'
    return
  endif
  write(*,*) "Interface closed"

end program gen_comp_grid_climo
